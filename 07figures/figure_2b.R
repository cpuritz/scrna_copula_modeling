library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(reshape2)
library(ggtext)
library(grid)
library(ggsci)
set.seed(0)

source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/04pca", full.names = TRUE)
files <- files[grepl("pcap", files)]
res_list <- lapply(files, readRDS)
for (i in seq_along(res_list)) {
    if (dim(res_list[[i]])[2] == 7) {
        res_list[[i]] <- res_list[[i]] %>%
            dplyr::select(ref, family, trial, sample, ngene, ncell, pval2 = pval)
    } else {
        res_list[[i]] <- res_list[[i]] %>%
            dplyr::select(ref, family, trial, sample, ngene, ncell, pval2)
    }
}
res <- do.call(rbind, res_list)

family_labels <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t")
names(family_labels) <- family_labels

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1)]
names(colors) <- family_labels

###############################################################################
## Main figure panel B ##

res <- res %>%
    group_by(ref, family) %>%
    summarise(pval = mean(pval2), .groups = "drop") %>%
    mutate(family = factor(case_match(
        family,
        "norm" ~ "Gaussian",
        "vine" ~ "Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "ind" ~ "Independence"
    ), levels = c("Independence", "Gaussian", "Jittered Gaussian",
                  "ML Gaussian", "t", "Vine", "Jittered Vine"))) %>%
    melt(id.vars = c("family", "ref"), measure.vars = "pval",
         value.name = "p") %>%
    mutate(p = -log10(p))

pplt_b <- ggplot(res, aes(x = family, y = p, fill = family)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(size = 1.1, position = position_dodge2(width = 0.24)) +
    geom_line(aes(group = ref), linewidth = 0.02, linetype = "dashed",
              color = "gray", position = position_dodge2(width = 0.24),
              alpha = 0.8) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                               title = NULL)) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab(expression(-log[10]~"(p)")) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.text.x = element_text(size = 13, angle = 30, hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          legend.position = "none",
          plot.margin = margin(5.5, 20, 5.5, 5.5))

###############################################################################
## Main figure panel A ##

pca_ex <- readRDS("Results/04pca/pca_ex.rds")
pplt_a <- ggplot(pca_ex$df,
                 aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c(colors, c("Reference" = "#424242"))) +
    xlab(NULL) +
    ylab("PC2") +
    facet_wrap(~plot, scales = "free", nrow = 1,
               labeller = labeller(plot = pca_ex$labels)) +
    ylim(c(min(pca_ex$df$PC2), max(pca_ex$df$PC2))) +
    annotation_custom(grob = textGrob("PC1", gp = gpar(fontsize = 13)),
                      xmin = -Inf, xmax = Inf, ymin = -70, ymax = -60) +
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
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.margin = margin(t = 3),
          plot.margin = margin(5.5, 5.5, 10.5, 5.5))

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
        plot_layout(nrow = 2, heights = c(0.47, 0.53))
ggsave(plot = pplt, filename = "Figures/figure_2.pdf", width = 8.25, height = 9)
