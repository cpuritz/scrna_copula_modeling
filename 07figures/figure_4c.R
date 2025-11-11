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

comp <- pairwise_wilcox_test(
    df = res[res$family != "Independence", ],
    cols = "Z",
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
    return(effsize::cohen.d(r1$Z, r2$Z, paired = TRUE)$estimate)
})

comp <- comp %>%
    rename(family1 = V1, family2 = V2, measure = var) %>%
    mutate(pval = sprintf("%.5e", pval),
           padj = sprintf("%.5e", padj),
           stat = sprintf("%.5f", stat)) %>%
    mutate(eff = sprintf("%.5f", abs(effs))) %>%
    arrange(desc(eff))
