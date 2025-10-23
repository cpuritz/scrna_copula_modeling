library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(tibble)
library(ggsignif)

set.seed(0)
setwd("~/Documents/Graduate School/Copula Paper")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

comp_all <- NULL
for (fig in c("c", "d")) {
    suffix <- ifelse(fig == "c", "", "_v2")
    res_norm  <- readRDS(paste0("Results/06wgcna/mod_norm", suffix, ".rds"))
    res_normj <- readRDS(paste0("Results/06wgcna/mod_norm_jitter", suffix, ".rds"))
    res_vinej <- readRDS(paste0("Results/06wgcna/mod_vine_jitter", suffix, ".rds"))
    res_boot  <- readRDS(paste0("Results/06wgcna/mod_boot", suffix, ".rds"))
    
    res <- rbind(
        melt(mutate(res_norm,  family = "norm"), id.vars = "family"),
        melt(mutate(res_normj, family = "norm_jitter"), id.vars = "family"),
        melt(mutate(res_vinej, family = "vine_jitter"), id.vars = "family"),
        melt(mutate(res_boot,  family = "boot"), id.vars = "family")
    )
    res <- res %>%
        dplyr::select(family, ARI = value) %>%
        mutate(family = case_match(
            family,
            "norm" ~ "Gaussian",
            "norm_jitter" ~ "Jittered Gaussian",
            "vine_jitter" ~ "Jittered Vine",
            "boot" ~ "Bootstrap"
        )) %>%
        mutate(family = factor(
            family,
            levels = c("Bootstrap", "Gaussian", "Jittered Gaussian",
                       "Jittered Vine")
        ))
    
    res <- res[!is.infinite(res$ARI), ]
    
    colors <- ggsci::pal_npg()(10)[c(8, 3, 5, 2)]
    names(colors) <- c("Bootstrap", "Jittered Vine", "Gaussian",
                       "Jittered Gaussian")
    
    comp <- pairwise_wilcox_test(
        df = res,
        cols = "ARI",
        group_var = "family",
        adjust_method = "fdr",
        alpha = NULL
    )
    comp$fig <- fig
    comp_all <- rbind(comp_all, comp)
    
    pplt <- ggplot(res, aes(x = family, y = ARI, fill = family)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 1.0, position = position_dodge2(width = 0.25)) +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(
            title = NULL,
            override.aes = list(shape = NA, linetype = 0)
        )) +
        xlab(NULL) +
        ylab("Adjusted Rand index") +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            axis.line = element_line(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 13),
            axis.text.x = element_text(size = 13, angle = 30, hjust = 1,
                                       color = "black"),
            axis.text.y = element_text(size = 13, color = "black"),
            legend.position = "none"
        )
    ggsave(plot = pplt, filename = paste0("Figures/figure_4", fig, ".pdf"),
           width = 7, height = 5)
}
comp_all$var <- NULL