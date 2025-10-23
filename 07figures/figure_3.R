library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggsci)
library(scales)
library(viridis)
library(patchwork)
set.seed(0)

files <- list.files("Results/05times", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))

res <- res %>%
    group_by(family, genes, cells) %>%
    summarize(time = mean(time), .groups = "drop") %>%
    mutate(family = case_match(
        family,
        "norm" ~ "Gaussian",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine" ~ "Vine",
        "vine_jitter" ~ "Jittered Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t"
    )) %>%
    mutate(family = factor(family, levels = c(
        "Gaussian", "Jittered Gaussian", "ML Gaussian", "t",
        "Vine", "Jittered Vine"
    )))

pplt <- ggplot(res, aes(x = genes, y = time)) +
    geom_point(aes(color = cells), size = 2.2) +
    facet_wrap(~ family, scales = "free_y") +
    scale_color_viridis() +
    scale_y_continuous(transform = transform_log10(),
                       labels = scientific_format(digits = 2)) +
    xlab("Number of genes") +
    ylab("Time (seconds)") +
    labs(color = "Number\nof cells") +
    theme_bw() +
    theme(
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13, color = "black"),
        panel.spacing.x = unit(1.0, "cm"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, margin = margin(b = 10)),
        legend.key.height = unit(1.2, "cm")
    )
ggsave(plot = pplt, filename = "Figures/figure_3.pdf", width = 9.75, height = 6.5)
