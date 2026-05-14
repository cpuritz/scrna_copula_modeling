library(dplyr)
library(ggplot2)
library(ggsci)
library(tibble)
library(tidyr)
library(reshape2)
set.seed(0)

source("scrna_copula_modeling/09figures/pairwise_wilcox_test.R")

files <- list.files("Results/05pca", full.names = TRUE)
files <- files[grepl(".rds", files) & !grepl("ex", files)]
res_list <- lapply(files, readRDS)
res <- do.call(rbind, res_list)

family_labels <- c("Independence", "Jittered Gaussian", "Sample Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t",
                   "ZINB-WaVE", "SPARSim")
names(family_labels) <- family_labels

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8, 9)]
names(colors) <- family_labels

res <- res %>%
    dplyr::group_by(ref, family, npc) %>%
    dplyr::summarise(p = mean(p), .groups = "drop") %>%
    dplyr::mutate(family = factor(dplyr::recode_values(
        family,
        "norm" ~ "Sample Gaussian",
        "vine" ~ "Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "ind" ~ "Independence",
        "sparsim" ~ "SPARSim",
        "zinbwave" ~ "ZINB-WaVE"
    ), levels = c("Independence", "Sample Gaussian", "Jittered Gaussian",
                  "ML Gaussian", "t", "Vine", "Jittered Vine",
                  "ZINB-WaVE","SPARSim"))) %>%
    dplyr::mutate(
        npc = factor(paste(npc, "PCs"), levels = paste(seq(2, 10), "PCs"))
    )

pplt <- ggplot(res, aes(x = family, y = p, fill = family)) +
    geom_boxplot(outlier.size = 1.2) +
    facet_wrap(~ npc) +
    scale_fill_manual(values = colors, labels = family_labels) +
    xlab(NULL) +
    ylab("FF p-value") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 13),
          legend.title = element_text(size = 13),
          axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                     color = "black"),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 13, margin = margin(r = 5)),
          legend.position = "none",
          strip.text = element_text(size = 13, color = "black"),
          plot.margin = margin(5.5, 30, 5.5, 5.5),
          panel.spacing.x = unit(0.7, "cm"))
ggsave(plot = pplt, filename = "Figures/figure_s2.pdf", width = 9.25, height = 9.5)
