library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(ggsignif)
library(ggsci)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper/")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/04pairwise", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))
rownames(res) <- NULL

nna <- sum(is.na(res$err))
message(nna, "/", dim(res)[1], " samples had NA error and will be discarded.")
res <- na.omit(res)

res <- res %>%
    tidyr::pivot_wider(names_from = stats, values_from = err) %>%
    dplyr::group_by(ref, family, trial) %>%
    dplyr::summarize(
        pearson = mean(pearson, na.rm = TRUE),
        spearman = mean(spearman, na.rm = TRUE),
        kendall = mean(kendall, na.rm = TRUE),
        mi = mean(mi, na.rm = TRUE),
        bicor = mean(bicor, na.rm = TRUE),
        dcor = mean(dcor, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::group_by(ref, family) %>%
    dplyr::summarize(
        Pearson = mean(pearson),
        Spearman = mean(spearman),
        Kendall = mean(kendall),
        "Mutual Information" = mean(mi),
        "Biweight Midcorrelation" = mean(bicor),
        "Distance Correlation" = mean(dcor),
        .groups = "drop"
    ) %>%
    dplyr::mutate(family = as.character(family)) %>%
    dplyr::mutate(family = recode_values(
        family,
        "ind" ~ "Independence",
        "norm" ~ "Gaussian",
        "vine" ~ "Vine",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t",
        "zinbwave" ~ "ZINB-WaVE",
        "sparsim" ~ "SPARSim"
    )) %>%
    dplyr::mutate(family = factor(
        family,
        levels = c("Independence", "Gaussian", "Jittered Gaussian",
                   "ML Gaussian", "t", "Vine", "Jittered Vine",
                   "ZINB-WaVE", "SPARSim")
    ))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8, 9)]
names(colors) <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t",
                   "ZINB-WaVE", "SPARSim")

vars <- setdiff(colnames(res), c("ref", "family"))
df <- melt(res, id.vars = c("family", "ref"), measure.vars = vars)
df$value <- log(df$value)

pplt <- ggplot(data = df, aes(x = family, y = value, fill = family)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(size = 0.7, position = position_dodge2(width = 0.22)) +
    geom_line(aes(group = ref), linewidth = 0.1, linetype = "dashed",
              color = "gray", position = position_dodge2(width = 0.22),
              alpha = 0.6
              ) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ variable, nrow = 2, scales = "free_y") +
    xlab(NULL) +
    ylab(expression(log(Frobenius~Error))) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 13, hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 13, color = "black"),
        panel.spacing.x = unit(0.9, "cm")
    )

ggsave(plot = pplt, filename = "Figures/figure_2.pdf", width = 12, height = 7)

eff <- function(f1, f2, m) {
    r1 <- res[res$family == f1, ]
    r2 <- res[res$family == f2, ]
    assertthat::assert_that(all(r1$ref == r2$ref))
    return(effsize::cohen.d(r1[[m]], r2[[m]], paired = TRUE)$estimate)
}

cops <- c("Gaussian", "Jittered Gaussian", "ML Gaussian", "t", "Vine",
          "Jittered Vine")
comp <- pairwise_wilcox_test(
    df = res[res$family %in% cops, ],
    cols = vars,
    group_var = "family",
    block_var = "ref",
    var_name = "var",
    adjust_method = "fdr"
)

# Swap families so that V1 always indicates the worse performer
for (i in seq_len(dim(comp)[1])) {
    if (comp[i, "stat"] < 0) {
        comp[i, "stat"] <- -comp[i, "stat"]
        f1 <- comp[i, "V1"]
        comp[i, "V1"] <- comp[i, "V2"]
        comp[i, "V2"] <- f1
    }
}

effs <- apply(comp, 1, function(r) { eff(r["V1"], r["V2"], r["var"]) })

comp <- comp %>%
    dplyr::rename(family1 = V1, family2 = V2, measure = var) %>%
    dplyr::mutate(pval = sprintf("%.5e", pval),
                  padj = sprintf("%.5e", padj),
                  stat = sprintf("%.5f", stat)) %>%
    dplyr::mutate(
        measure = dplyr::recode_values(
            measure,
            "Pearson" ~ "pearson",
            "Kendall" ~ "kendall",
            "Spearman" ~ "spearman",
            "Distance Correlation" ~ "dcor",
            "Biweight Midcorrelation" ~ "bicor",
            "Mutual Information" ~ "MI"),
        eff = sprintf("%.5f", abs(effs))
    ) %>%
    dplyr::arrange(dplyr::desc(eff))

write.csv(
    x = comp,
    file = "Tables/pairwise_pvalues.csv",
    row.names = FALSE,
    quote = FALSE
)