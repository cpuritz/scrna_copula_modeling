library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
set.seed(0)
setwd("~/Documents/Graduate School/Copula Paper")

families <- c("norm" = "Gaussian",
              "vine" = "Vine",
              "nmle" = "ML-Gaussian",
              "t" = "t")
res <- do.call(rbind, lapply(names(families), function(f) {
    files <- list.files("Results/04times", pattern = paste0("^", f, "_\\d+"),
                        full.names = TRUE)
    do.call(rbind, lapply(files, readRDS))
}))
res$s_per_cg <- res$time / res$cells / res$genes
res$log_s_per_cg <- log10(res$s_per_cg)

plots <- list()
bins <- c(12, 12, 12, 12)
for (i in seq_along(families)) {
    df <- filter(res, family == !!names(families)[i])
    print(paste0(families[i], ": ", sprintf("%.0e", mean(df$s_per_cg))))
    plots[[i]] <- ggplot(df, aes(x = log_s_per_cg)) +
        geom_histogram(bins = bins[i], fill = "#4DBBD5FF",
                       color = "black") +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = expression(log[10]~paste("(second / (cell x gene))"))) +
        ylab("Frequency") +
        ggtitle(families[i]) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(),
              axis.title = element_text(size = 13),
              axis.text = element_text(size = 13, color = "black"),
              plot.title = element_text(size = 13, hjust = 0.5))
}

pplt <- wrap_plots(plots, ncol = 2)
ggsave(filename = "Figures/time.pdf", plot = pplt, width = 9, height = 7)
