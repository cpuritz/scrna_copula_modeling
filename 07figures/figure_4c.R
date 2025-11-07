library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("~/Documents/Graduate School/Copula Paper")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/06wgcna", full.names = TRUE)
files <- files[grepl(".rds", files)]
files <- files[!grepl("boot", files)]
res <- do.call(rbind, lapply(files, readRDS))
res <- res %>%
    group_by(family, ref, trial) %>%
    summarize(Z = mean(Z), .groups = "drop") %>%
    group_by(family, ref) %>%
    summarize(Z = mean(Z), .groups = "drop") %>%
    mutate(family = case_match(
        family,
        "norm" ~ "Gaussian",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "ind" ~ "Independence"
    )) %>%
    mutate(family = factor(
        family,
        levels = c("Independence", "Gaussian", "Jittered Gaussian",
                   "Jittered Vine")
    ))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3)]
names(colors) <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine")

comp <- boxplot_signif(
    df = res,
    cols = "Z",
    group_var = "family",
    block_var = "ref",
    annot = "stars",
    adjust_method = "fdr",
    offset = 0.03,
    scale = 0.04
)
comp$y <- min(comp$y) + c(0, 0.5, 3.5, 1.5, 2.5, 0) * diff(sort(comp$y))[1]

pplt <- ggplot(res, aes(x = family, y = Z)) +
    geom_boxplot(aes(fill = family), outliers = FALSE) +
    scale_fill_manual(values = colors) +
    geom_point(position = position_dodge2(width = 0.22), color = "black", size = 1.8) +
    geom_line(aes(group = ref), linewidth = 0.1, linetype = "solid",
              color = "gray", position = position_dodge2(width = 0.22),
              alpha = 0.8) +
    xlab(NULL) +
    ylab(expression(bar(Z)[summary])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    )
ggsave(pplt, filename = "Figures/figure_4c.pdf", width = 8, height = 5.25)