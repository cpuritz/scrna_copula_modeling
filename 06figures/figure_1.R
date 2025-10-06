library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(ggsignif)
library(ggsci)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper")
source("Scripts/06figures/pairwise_wilcox_test.R")

if (!"res_all" %in% ls()) {
    res_all <- do.call(rbind, lapply(list.files("Results/comps", full.names = TRUE), readRDS))
}
res <- res_all
rownames(res) <- NULL
res <- res[grepl("norm|vine", res$family), ]

res <- res %>%
    dplyr::filter(norm == "two") %>%
    dplyr::select(-norm) %>%
    tidyr::pivot_wider(names_from = stat, values_from = err) %>%
    dplyr::group_by(ref, family, trial) %>%
    dplyr::summarize(
        ngene = unique(ngene),
        ncell = unique(ncell),
        pearson = mean(pearson),
        spearman = mean(spearman),
        kendall = mean(kendall),
        mi = mean(mi),
        bicor = mean(bicor),
        dCor = mean(dcor),
        .groups = "drop"
    ) %>%
    na.omit() %>%
    dplyr::group_by(ref, family) %>%
    dplyr::summarize(
        ngene = unique(ngene),
        ncell = unique(ncell),
        Pearson = mean(pearson),
        Kendall = mean(kendall),
        Spearman = mean(spearman),
        MI = mean(mi),
        Bicor = mean(bicor),
        dCor = mean(dCor),
        .groups = "drop"
    ) %>%
    dplyr::mutate(family = as.character(family)) %>%
    dplyr::mutate(family = case_match(
        family,
        "norm" ~ "Gaussian (count)",
        "vine" ~ "Vine (count)",
        "norm_jitter" ~ "Gaussian (jittered)",
        "vine_jitter" ~ "Vine (jittered)",
        .default = family
    ))

fv <- list(c("Gaussian (count)", "Gaussian (jittered)"),
           c("Vine (count)", "Vine (jittered)"))
plots <- list()
sbw <- c(0.10, 0.04, 0.04,
         0.05, 0.015, 0.08,
         0.02, 0.02, 0.02,
         0.02, 0.015, 0.02)
all_comp <- NULL
for (i in seq_along(fv)) {
    res_sub <- res[res$family %in% fv[[i]], ]
    vars <- setdiff(colnames(res), c("ref", "family", "ngene", "ncell"))
    bp <- boxplot_signif(res_sub,
                         cols = vars,
                         group_var = "family",
                         block_var = "ref",
                         var = "var",
                         adjust_method = "fdr",
                         alpha = NULL,
                         annot = "stars",
                         transform = log)
    all_comp <- rbind(all_comp, bp)
    df <- melt(res_sub, id.vars = c("family", "ref"), measure.vars = vars)
    for (j in seq_along(vars)) {
        data <- df[df$variable == vars[j], ]
        comp <- bp[bp$var == vars[j], ]
        comp$y <- comp$y * 0.98
        data$value <- log(data$value)
        
        if (i == 1) {
            colors <- ggsci::pal_npg()(10)[c(2, 5)]
        } else {
            colors <- ggsci::pal_npg()(10)[c(3, 4)]
        }

        ix <- (i - 1) * length(vars) + j
        plots[[ix]] <- 
            ggplot(data, aes(x = family, y = value, fill = family)) +
            geom_boxplot(outlier.shape = NA) +
            geom_point(aes(group = ref),
                       position = position_dodge(width = 0.22),
                       size = 0.6) +
            geom_line(aes(group = ref), position = position_dodge(width = 0.22),
                      linetype = "solid", color = "gray", linewidth = 0.15,
                      alpha = 0.7) +
            guides(fill = guide_legend(
                title = NULL,
                override.aes = list(shape = NA, linetype = 0)
            )) +
            xlab(NULL) +
            ylab(NULL) +
            ggtitle(vars[j]) +
            scale_fill_manual(values = colors) +
            theme_bw() +
            theme(panel.border = element_blank(),
                  axis.line = element_line(),
                  panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 13),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(size = 13),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size = 13, color = "black"),
                  axis.ticks.x = element_blank())
        if (j == 1) {
            plots[[ix]] <- plots[[ix]] +
                ylab("log(Frobenius error)") +
                theme(axis.title.y = element_text(size = 13))
        }
        
        vjust <- ifelse(comp$annot == "n.s.", -0.1, 0.4)
        plots[[ix]] <- plots[[ix]] +
            geom_signif(xmin = comp$xmin,
                        xmax = comp$xmax,
                        y_position = 0.925 * comp$y,
                        annotation = comp$annot,
                        tip_length = 0,
                        vjust = vjust,
                        textsize = 4) +
            scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                               limits = c(NA, 0.93 * max(comp$y)))
    }
}
all_comp <- all_comp %>% dplyr::select(V1, V2, var, padj, stat)
pplt_a <- wrap_plots(plots[1:6], nrow = 1) + plot_layout(guides = "collect")
pplt_b <- wrap_plots(plots[7:12], nrow = 1) + plot_layout(guides = "collect")

###############################################################################

get_times <- function(x) {
    files <- list.files("Results/04times", pattern = paste0("^", x, "_\\d+"),
                        full.names = TRUE)
    times <- do.call(rbind, lapply(files, readRDS))
    times <- times %>%
        dplyr::group_by(cells, genes) %>%
        dplyr::summarize(time = mean(time), .groups = "drop")
    return(times)
}
times <- list(
    list(count = get_times("norm"), jitter = get_times("norm_jitter")),
    list(count = get_times("vine"), jitter = get_times("vine_jitter"))
)
plots <- list()
pvals <- c()
for (i in seq_along(times)) {
    time <- dplyr::full_join(times[[i]]$count,
                             times[[i]]$jitter,
                             by = c("genes", "cells"),
                             suffix = c(".c", ".j")) %>%
        dplyr::mutate(ratio = time.c / time.j)
    
    # Paired Wilcox rank sum test comparing count and jitter times
    pvals[i] <- time %>%
        reshape2::melt(
            id.vars = c("cells", "genes"),
            measure.vars = c("time.c", "time.j"),
            variable.name = "group"
        ) %>%
        dplyr::mutate(block = factor(paste(cells, genes, sep = '_'))) %>%
        coin::wilcoxsign_test(
            value ~ group | block,
            data = .,
            distribution = "exact"
        ) %>%
        coin::pvalue()
    
    title <- ifelse(i == 1, "Gaussian", "Vine")
    binwidth <- ifelse(i == 1, 0.10, 50)
    
    plots[[i]] <- ggplot(time, aes(x = time.c, y = time.j)) +
        geom_point(size = 1.6) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("Time to fit count model (sec)") +
        ylab("Time to fit jittered model (sec)") +
        ggtitle(paste(title, "copula")) +
        theme_bw() +
        theme(axis.title = element_text(size = 13),
              plot.title = element_text(hjust = 0.5, size = 13),
              axis.text = element_text(size = 13, color = "black"),
              panel.grid = element_blank())
}
pplt_c <- plots[[1]]
pplt_d <- plots[[2]]

pplt_c <- pplt_c +
    annotate(geom = "text", x = 0.092, y = 0.40,
             label = paste("p =", sprintf("%.02f", pvals[1])),
             size = 4.3)
pplt_d <- pplt_d +
    annotate(geom = "text", x = 4.6, y = 180,
             label = paste("p =", sprintf("%.02e", pvals[2])),
             size = 4.3)

get_annot <- function(x) {
    plot_annotation(
        title = x,
        theme = theme(plot.title = element_text(face = 2, size = 20))
    )
}
pplt <- wrap_elements(pplt_a + get_annot("a")) +
    wrap_elements(pplt_b + get_annot("b")) +
    wrap_elements(pplt_c + get_annot("c")) +
    wrap_elements(pplt_d + get_annot("d")) +
    plot_layout(design = "1111\n222#\n34##",
                heights = c(0.32, 0.32, 0.38),
                widths = c(0.41, 0.41, 0.15, 0.03))
ggsave(plot = pplt, filename = "Figures/jitter.pdf", width = 11, height = 11)
