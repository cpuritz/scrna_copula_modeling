library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(ggsignif)
library(scales)
library(ggsci)
library(cowplot)
set.seed(0)
setwd("~/Documents/Graduate School/Copula Paper")
source("Scripts/Misc/pairwise_wilcox_test.R")

#res_all <- do.call(rbind, lapply(list.files("Results/comps", full.names = TRUE), readRDS))
res <- res_all
rownames(res) <- NULL
res <- res[!grepl("jitter", res$family), ]
res_raw <- res

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
    dplyr::mutate(
        family = dplyr::case_match(
            family,
            "norm" ~ "Gaussian",
            "nmle" ~ "ML-Gaussian",
            "vine" ~ "Vine",
            "t" ~ "t",
            .default = family)
    )

family_labels <- c("Gaussian", "ML-Gaussian", "t", "Vine")
names(family_labels) <- family_labels
res$family <- res$family %>%
    factor(levels = names(family_labels), ordered = TRUE) %>%
    droplevels()

vars <- setdiff(colnames(res), c("ref", "family", "ngene", "ncell"))
comps <- pairwise_wilcox_test(
    res,
    cols = vars,
    group_var = "family",
    block_var = "ref"
)
stats <- unique(comps$var)

pairs <- t(combn(unique(c(comps$V1, comps$V2)), 2))
effs <- c()
for (ix in seq_along(stats)) {
    stat <- stats[ix]
    vecs <- lapply(seq(dim(pairs)[1]), function(i) {
        res %>%
            filter(family %in% pairs[i, ]) %>%
            select(ref, family, stat) %>%
            pivot_wider(names_from = family, values_from = stat) %>%
            mutate(diff = .data[[pairs[i, 1]]] - .data[[pairs[i, 2]]]) %>%
            select(diff) %>%
            unlist()
    })
    effs <- c(effs, list(sapply(vecs, function(x) { abs(mean(x)) / sd(x) })))
}
names(effs) <- stats
palette <- colorRampPalette(c("#D9D9D9", "#008BB8"))(100)
colors <- palette[round(scales::rescale(unlist(effs), to = c(1, 100)))]

plots <- list()
binwidths <- c(20, 20, 20, 20, 20, 20,
               20, 20, 20, 20, 20, 20,
               20, 20, 20, 20, 20, 20,
               20, 20, 20, 20, 20, 20,
               20, 20, 20, 20, 20, 20,
               20, 20, 20, 20, 20, 20)
px <- c(NA, NA, NA, NA, NA, NA,
        NA, NA, NA, -0.08, -0.7, -0.7,
        NA, NA, NA, -0.08, -1, NA,
        0.25, 0.25, 0.3, 0.01, 0.15, 0.15,
        NA, NA, NA, -0.17, -0.25, NA,
        0.35, 0.4, NA, NA, NA, NA)
py <- c(NA, NA, NA, NA, NA, NA,
        NA, NA, NA, 25, 25, 25,
        NA, NA, NA, 26, 26, NA,
        16, 16, 16, 16, 16, 16,
        NA, NA, NA, 20, 20, NA,
        27, 27, NA, NA, NA, NA)
for (ix in seq_along(stats)) {
    stat <- stats[ix]
    pairs <- t(combn(unique(c(comps$V1, comps$V2)), 2))
    vecs <- lapply(seq(dim(pairs)[1]), function(i) {
        res %>%
            filter(family %in% pairs[i, ]) %>%
            select(ref, family, stat) %>%
            pivot_wider(names_from = family, values_from = stat) %>%
            mutate(diff = .data[[pairs[i, 1]]] - .data[[pairs[i, 2]]]) %>%
            select(diff) %>%
            unlist()
    })
    eff <- effs[[stat]]
    
    hists <- lapply(seq_along(vecs), function(i) {
        lix <- (ix - 1) * 6 + i
        padj <- comps[lix, "padj"]
        if (padj <= 0.05) {
            title <- paste("p =", sprintf("%.01e", padj))
        } else {
            title <- NULL
        }

        p <- ggplot(data.frame(x = vecs[[i]]), aes(x = x)) +
            geom_histogram(bins = binwidths[lix], color = "black",
                           fill = colors[lix]) +
            geom_vline(xintercept = 0, color = "lightgray", linewidth = 0.1,
                       alpha = 0.8) +
            xlab(NULL) +
            ylab(NULL) +
            ggtitle(title) +
            theme_minimal() +
            theme(panel.grid = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 8),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = "none",
                  plot.title = element_text(hjust = 0.5, size = 9))
 
        xmax <- max(abs(unlist(ggplot_build(p)$data[[1]][c("xmax", "xmin")])))
        p <- p + 
            scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
            coord_cartesian(xlim = c(-xmax, xmax)) +
            geom_blank(data = data.frame(x = c(-xmax, xmax)), aes(x = x))
        
        if (i == 1) {
            p <- p + theme(
                axis.text.y = element_text(size = 8),
                axis.ticks.y = element_line()
            )
        }
        if (ix == 1) {
            f1 <- pairs[i, 1]
            f2 <- pairs[i, 2]
            if (f1 == "ML-Gaussian") {
                f1 <- "MLG"
            } else if (f2 == "ML-Gaussian") {
                f2 <- "MLG"
            }
            ytitle <- paste0(f1, " - ", f2)
            label_panel <- ggplot() +
                annotate("text", x = 0, y = 0.5, label = ytitle, angle = 0,
                         size = 4, fontface = "bold") +
                theme_void() +
                theme(plot.background = element_rect(fill = "#D9D9D9",
                                                     color = NA))
            p <- patchwork::wrap_plots(list(label_panel, p), nrow = 2,
                                       heights = c(0.2, 0.8))
        }
        return(p)
    })
    
    ymax <- max(sapply(hists, function(x) {
        summarise(ggplot_build(x)$data[[1]], max_count = max(count))[[1]]
    }))
    for (i in seq_along(hists)) {
        hists[[i]] <- hists[[i]] +
            scale_y_continuous(limits = c(0, ymax),
                               breaks = scales::breaks_pretty(n = 3),
                               expand = expansion(mult = c(0, 0.05)))
    }
    
    slab <- ifelse(stat == "Pearson", "     Pearson", stat)
    label_panel <- ggplot() +
        annotate("text", x = 0, y = 0.5, label = slab, angle = 270,
                 size = 4, fontface = "bold") +
        theme_void() +
        theme(plot.background = element_rect(fill = "#D9D9D9", color = NA))
    plots[[ix]] <- patchwork::wrap_plots(c(hists, list(label_panel)),
                                  nrow = 1, widths = c(rep(1, 6), 0.15))
}
pplt <- wrap_plots(plots, ncol = 1, heights = c(1.31, rep(1, 5)))

eff_range <- seq(min(unlist(effs)), max(unlist(effs)),
                 length.out = length(palette))
legend_df <- data.frame(x = 1, y = eff_range, eff = eff_range)
legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = eff)) +
    geom_tile() +
    scale_fill_gradientn(colors = palette, name = "Effect\nsize") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.height = unit(0.8, "cm"))
pplt <- wrap_elements(pplt) +
    wrap_elements(full = cowplot::get_legend(legend_plot)) +
    plot_layout(ncol = 2, widths = c(1, 0.08))

ylab_panel <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = "Frequency", angle = 90,
             size = 4.3) +
    theme_void() +
    coord_cartesian(clip = "off")
pplt <- patchwork::wrap_plots(list(ylab_panel, plot_spacer(), pplt), nrow = 1,
                              widths = c(0.01, -0.01, 1))

xlab_panel <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = "Difference in error          ",
             angle = 0, size = 4.3) +
    theme_void() +
    coord_cartesian(clip = "off")
pplt <- patchwork::wrap_plots(list(pplt, plot_spacer(), xlab_panel), ncol = 1,
                              heights = c(1, -0.03, 0.01))

ggsave(plot = pplt, filename = "Figures/pairs.pdf", width = 10.2, height = 9)
