library(dplyr)
library(ggplot2)
library(ggsci)
library(tibble)
library(tidyr)
library(reshape2)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/04pca", full.names = TRUE)
files <- files[grepl("pcap", files)]
res_list <- lapply(files, readRDS)
res <- do.call(rbind, res_list)

family_labels <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t")
names(family_labels) <- family_labels
colors <- pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1)]
names(colors) <- family_labels

res <- res %>%
    group_by(ref, family) %>%
    summarise(
        pval2 = mean(pval2),
        pval3 = mean(pval3),
        pval4 = mean(pval4),
        pval5 = mean(pval5),
        pval6 = mean(pval6),
        pval7 = mean(pval7),
        pval8 = mean(pval8),
        pval9 = mean(pval9),
        pval10 = mean(pval10),
        .groups = "drop"
    ) %>%
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
    melt(id.vars = c("family", "ref"), variable.name = "npc",
         value.name = "p") %>%
    mutate(
        p = p,
        npc = factor(sapply(npc, function(x) {
            paste(as.integer(gsub("pval", "", x)), "PCs")
        }), levels = paste(seq(2, 10), "PCs"))
    )

pplt <- ggplot(res, aes(x = family, y = p, fill = family)) +
    geom_boxplot(outlier.size = 0.9) +
    facet_wrap(~ npc) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab("p-value") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          legend.position = "none",
          strip.text = element_text(size = 13, color = "black"),
          plot.margin = margin(5.5, 30, 5.5, 5.5),
          panel.spacing.x = unit(0.7, "cm"))
ggsave(plot = pplt, filename = "Figures/figure_s2.pdf", width = 9.25, height = 9.5)
