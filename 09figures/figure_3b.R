library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(grid)
library(ggsci)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper/")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/05pca", full.names = TRUE)
files <- files[grepl(".rds", files) & !grepl("ex", files)]
res_list <- lapply(files, readRDS)
res <- do.call(rbind, res_list)

family_labels <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t",
                   "ZINB-WaVE", "SPARSim")
names(family_labels) <- family_labels

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8, 9)]
names(colors) <- family_labels

###############################################################################
## Main figure panel B ##

res <- res %>%
    dplyr::filter(npc == 10) %>%
    dplyr::select(-npc) %>%
    group_by(ref, family) %>%
    summarise(pval = mean(p), .groups = "drop") %>%
    mutate(family = factor(recode_values(
        family,
        "norm" ~ "Gaussian",
        "vine" ~ "Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "ind" ~ "Independence",
        "zinbwave" ~ "ZINB-WaVE",
        "sparsim" ~ "SPARSim"
    ), levels = c("Independence", "Gaussian", "Jittered Gaussian",
                  "ML Gaussian", "t", "Vine", "Jittered Vine", "ZINB-WaVE",
                  "SPARSim")
    )) %>%
    melt(id.vars = c("family", "ref"), measure.vars = "pval",
         value.name = "p")

pplt_b <- ggplot(res, aes(x = family, y = p, fill = family)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(size = 1.1, position = position_dodge2(width = 0.24)) +
    geom_line(aes(group = ref), linewidth = 0.25, linetype = "dashed",
              color = "gray", position = position_dodge2(width = 0.24),
              alpha = 0.5) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                               title = NULL)) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab("p-value") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13, margin = margin(r = 5)),
          legend.title = element_text(size = 13),
          axis.text.x = element_text(size = 13, angle = 30, hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          legend.position = "none",
          plot.margin = margin(5.5, 20, 5.5, 5.5))

comp <- pairwise_wilcox_test(
    df = res[!res$family %in% c("Independence", "ZINB-WaVE", "SPARSim"), ],
    cols = "p",
    group_var = "family",
    block_var = "ref",
    adjust_method = "fdr"
)

# Swap families so that V1 always indicates the better performer
for (i in seq_len(dim(comp)[1])) {
    if (comp[i, "stat"] < 0) {
        comp[i, "stat"] <- -comp[i, "stat"]
        f1 <- comp[i, "V1"]
        comp[i, "V1"] <- comp[i, "V2"]
        comp[i, "V2"] <- f1
    }
}

effs <- apply(comp, 1, function(r) {
    r1 <- res[res$family == r["V1"], ]
    r2 <- res[res$family == r["V2"], ]
    assertthat::assert_that(all(r1$ref == r2$ref))
    return(effsize::cohen.d(r1$p, r2$p, paired = TRUE)$estimate)
})

comp <- comp %>%
    dplyr::rename(family1 = V1, family2 = V2) %>%
    dplyr::select(-var) %>%
    mutate(pval = sprintf("%.5e", pval),
           padj = sprintf("%.5e", padj),
           stat = sprintf("%.5f", stat)) %>%
    mutate(eff = sprintf("%.5f", abs(effs))) %>%
    arrange(desc(eff))

write.csv(
    x = comp,
    file = "Tables/pca_pvalues.csv",
    row.names = FALSE,
    quote = FALSE
)

###############################################################################
## Main figure panel A ##

pca_ex <- readRDS("Results/05pca/pca_ex.rds")
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
    theme(axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13, margin = margin(r = 5)),
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
          legend.box.margin = margin(t = 5),
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
ggsave(plot = pplt, filename = "Figures/figure_3.pdf", width = 8.25, height = 9)
