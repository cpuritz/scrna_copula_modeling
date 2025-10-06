library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(reshape2)
library(ggsignif)
library(ggtext)
library(grid)
library(ggsci)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper")
source("Scripts/06figures/pairwise_wilcox_test.R")

files <- list.files("Results/PCA", full.names = TRUE)
files <- files[grepl("pcaall", files)]
res_all <- do.call(rbind, lapply(files, readRDS))

family_labels <- c("Gaussian", "ML-Gaussian", "t", "Vine")
names(family_labels) <- family_labels

colors <- ggsci::pal_npg()(10)[c(2, 3, 6, 1)]
names(colors) <- c("Gaussian", "Vine", "ML-Gaussian", "t")

###############################################################################
## Supplemental figure ##

plots <- list()
for (i in seq(2, 10)) {
    res <- res_all %>%
        group_by(ref, family) %>%
        summarise(
            ff = mean(.data[[paste0("ff", i)]]),
            .groups = "drop"
        ) %>%
        mutate(family = factor(case_match(
            family,
            "norm" ~ "Gaussian",
            "vine" ~ "Vine",
            "nmle" ~ "ML-Gaussian",
            "t" ~ "t"
        ), levels = names(family_labels), ordered = TRUE))
    
    bp <- boxplot_signif(res,
                         cols = "ff",
                         group_var = "family",
                         block_var = "ref",
                         alpha = 0.05,
                         annot = "stars",
                         transform = log10)
    bp$y <- 0.88 * min(bp$y) + (order(bp$y) - 1) * diff(sort(bp$y))[1] * 0.15
    ix <- seq_along(bp$y)[-which.min(bp$y)]
    bp$y[ix] <- bp$y[ix] - diff(sort(bp$y))[1]
    df <- melt(res, id.vars = c("family", "ref"), measure.vars = "ff")
    df$value <- log10(df$value)
    
    ylabel <- NULL
    if (i %in% c(2, 5, 8)) {
        ylabel <- expression(log[10]~"(FF statistic)")
    }
    
    plots[[i - 1]] <- ggplot(df, aes(x = family, y = value, fill = family)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(aes(group = ref), position = position_dodge(width = 0.3),
                   size = 0.3) +
        geom_line(aes(group = ref), position = position_dodge(width = 0.3),
                  linetype = "solid", color = "gray", linewidth = 0.10,
                  alpha = 0.5) +
        guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                                   title = NULL)) +
        scale_fill_manual(values = colors, labels = family_labels) +
        xlab(NULL) +
        ylab(ylabel) +
        ggtitle(paste(i, "PCs")) +
        theme_bw() +
        theme(panel.border = element_blank(),
              axis.line = element_line(),
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.title = element_text(size = 13),
              legend.title = element_text(size = 13),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.ticks.x = element_blank(),
              legend.text = element_markdown(size = 13)) +
        ylim(c(NA, max(bp$y))) +
        geom_signif(xmin = bp$xmin,
                    xmax = bp$xmax,
                    y_position = bp$y,
                    annotation = bp$annot,
                    tip_length = 0,
                    vjust = 0.4,
                    textsize = 3.5)
}
pplt_sup <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect")
ggsave(filename = "Figures/pca_supp.pdf", width = 9, height = 9)

###############################################################################
## Main figure panel B ##

res <- res_all %>%
    group_by(ref, family) %>%
    summarise(ff = mean(ff2), .groups = "drop") %>%
    mutate(family = factor(case_match(
        family,
        "norm" ~ "Gaussian",
        "vine" ~ "Vine",
        "nmle" ~ "ML-Gaussian",
        "t" ~ "t"
    ), levels = names(family_labels), ordered = TRUE))
transform <- log10
bp <- boxplot_signif(res,
                     cols = "ff",
                     group_var = "family",
                     block_var = "ref",
                     alpha = 0.05,
                     annot = "stars",
                     transform = transform)
bp$y <- 0.88 * min(bp$y) + (order(bp$y) - 1) * diff(sort(bp$y))[1] * 0.15
ix <- seq_along(bp$y)[-which.min(bp$y)]
bp$y[ix] <- bp$y[ix] - diff(sort(bp$y))[1]
df_b <- melt(res, id.vars = c("family", "ref"), measure.vars = "ff")
df_b$value <- transform(df_b$value)

pplt_b <- ggplot(df_b, aes(x = family, y = value, fill = family)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group = ref), position = position_dodge(width = 0.30),
               size = 0.7) +
    geom_line(aes(group = ref), position = position_dodge(width = 0.30),
              linetype = "solid", color = "gray", linewidth = 0.15,
              alpha = 0.5) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                               title = NULL)) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab(expression(log[10]~"(FF statistic)")) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.ticks.x = element_blank(),
          legend.text = element_markdown(size = 13)) +
    ylim(c(NA, max(bp$y))) +
    geom_signif(xmin = bp$xmin,
                xmax = bp$xmax,
                y_position = bp$y,
                annotation = bp$annot,
                tip_length = 0,
                vjust = 0.4,
                textsize = 3.5)

###############################################################################
## Main figure panel C ##

res_ref <- res %>%
    group_by(ref) %>%
    summarise(ref_mean = mean(ff),
              ref_sd = sqrt(var(ff)),
              ref_max = max(ff),
              ref_min = min(ff))
ref_means <- setNames(res_ref$ref_mean, res_ref$ref)
ref_sdevs <- setNames(res_ref$ref_sd, res_ref$ref)
res$z <- (res$ff - ref_means[res$ref]) / ref_sdevs[res$ref]
df <- melt(res, id.vars = c("family", "ref"), measure.vars = "z")
pplt_c <- ggplot(df, aes(x = family, y = value, fill = family)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group = ref), position = position_dodge(width = 0.30),
               size = 0.7) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                               title = NULL)) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab("z-score") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.ticks.x = element_blank(),
          legend.text = element_markdown(size = 13))

###############################################################################
## Main figure panel A ##

pca_ex <- readRDS("Results/pca_example.rds")
pplt_a <- ggplot(pca_ex$df,
                 aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(colors, c("Reference" = "#B0B0B0"))) +
    xlab(NULL) +
    ylab("PC2") +
    facet_wrap(~plot, scales = "free", nrow = 1,
               labeller = labeller(plot = pca_ex$labels)) +
    ylim(c(min(pca_ex$df$PC2), max(pca_ex$df$PC2))) +
    annotation_custom(grob = textGrob("PC1", gp = gpar(fontsize = 13)),
                      xmin = -Inf, xmax = Inf, ymin = -80, ymax = -79) +
    coord_cartesian(clip = "off") +
    guides(color = guide_legend(override.aes = list(size = 2.2))) +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 13),
          strip.background = element_blank(),
          strip.text = element_text(size = 13),
          legend.title = element_blank())

###############################################################################
## Main figure ##

get_annot <- function(x) {
    plot_annotation(
        title = x,
        theme = theme(plot.title = element_text(face = 2, size = 20))
    )
}

pplt <- wrap_elements(pplt_a + get_annot("a")) +
        wrap_elements(pplt_b + get_annot("b")) +
        wrap_elements(pplt_c + get_annot("c")) +
        plot_layout(nrow = 3, heights = c(0.32, 0.34, 0.34))
ggsave(plot = pplt, filename = "Figures/pca.pdf", width = 8, height = 10)
