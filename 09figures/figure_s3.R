library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggsci)
library(scales)
library(viridis)
library(patchwork)
set.seed(0)

files <- list.files("Results/07times", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))

res <- res %>%
    dplyr::group_by(family, genes, cells) %>%
    dplyr::summarize(time = mean(time), .groups = "drop") %>%
    dplyr::mutate(family = dplyr::recode_values(
        family,
        "norm" ~ "Sample Gaussian",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine" ~ "Vine",
        "vine_jitter" ~ "Jittered Vine",
        "nmle" ~ "ML-Gaussian",
        "t" ~ "t"
    )) %>%
    dplyr::mutate(family = factor(family, levels = c(
        "Sample Gaussian", "Vine", "ML-Gaussian", "Jittered Gaussian",
        "Jittered Vine", "t"
    )))

time_ratio <- function(x, y) {
    rx <- res[res$family == x, ]
    ry <- res[res$family == y, ]
    assertthat::assert_that(all(rx$cells == ry$cells), all(rx$genes == ry$genes))
    return(data.frame(x = rx$time / ry$time))
}

hists <- list()
hists[[1]] <- ggplot(time_ratio("Sample Gaussian", "Jittered Gaussian"), aes(x)) +
    geom_histogram(bins = 15, fill = "#4DBBD5FF", color = "black") +
    xlab("Ratio") +
    ylab("Frequency") +
    ggtitle("Sample Gaussian / Jittered Gaussian") +
    scale_y_continuous(expand = c(0, 0, 0, 0.4),
                       breaks = pretty_breaks(n = 4))
hists[[2]] <- plot_spacer()
hists[[3]] <- ggplot(time_ratio("Vine", "Jittered Vine"), aes(x)) +
    geom_histogram(bins = 15, fill = "#4DBBD5FF", color = "black") +
    xlab("Ratio") +
    ylab(NULL) +
    ggtitle("Vine / Jittered Vine") +
    scale_y_continuous(expand = c(0, 0, 0, 0.22))
hists[[4]] <- plot_spacer()
hists[[5]] <- ggplot(time_ratio("t", "ML-Gaussian"), aes(x)) +
    geom_histogram(bins = 11, fill = "#4DBBD5FF", color = "black") +
    xlab("Ratio") +
    ylab(NULL) +
    ggtitle("t / ML Gaussian") +
    scale_y_continuous(expand = c(0, 0, 0, 0.16))

pplt <- wrap_plots(hists, nrow = 1, widths = c(1, 0.01, 1, 0.01, 1)) &
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13, margin = margin(r = 7)),
        axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(size = 13, hjust = 0.5)
    )
ggsave(plot = pplt, filename = "Figures/figure_s3.pdf", width = 10.75, height = 4)
