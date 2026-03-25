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

files <- list.files("Results/03samples", full.names = TRUE)
files <- files[grepl(".rds", files)]
res <- do.call(rbind, lapply(files, readRDS))
rownames(res) <- NULL
res <- res[res$family != "boot", ]

res <- res %>%
    dplyr::filter(norm == "two" & family != "boot") %>%
    dplyr::mutate(rmse = err / ngene) %>%
    dplyr::select(-norm, -err) %>%
    tidyr::pivot_wider(names_from = stat, values_from = rmse) %>%
    dplyr::group_by(ref, family, trial) %>%
    dplyr::summarize(
        ngene = unique(ngene),
        ncell = unique(ncell),
        pearson = mean(pearson),
        spearman = mean(spearman),
        kendall = mean(kendall),
        mi = mean(mi),
        bicor = mean(bicor),
        dCor = mean(dcor),
        .groups = "drop"
    ) %>%
    na.omit() %>%
    dplyr::group_by(ref, family) %>%
    dplyr::summarize(
        ngene = unique(ngene),
        ncell = unique(ncell),
        Pearson = mean(pearson),
        Kendall = mean(kendall),
        Spearman = mean(spearman),
        'Mutual information' = mean(mi),
        'Biweight midcorrelation' = mean(bicor),
        'Distance correlation' = mean(dCor),
        .groups = "drop"
    ) %>%
    dplyr::mutate(family = as.character(family)) %>%
    dplyr::mutate(family = case_match(
        family,
        "ind" ~ "Independence",
        "norm" ~ "Gaussian",
        "vine" ~ "Vine",
        "norm_jitter" ~ "Jittered Gaussian",
        "vine_jitter" ~ "Jittered Vine",
        "nmle" ~ "ML Gaussian",
        "t" ~ "t",
        .default = family
    )) %>%
    dplyr::mutate(family = factor(
        family,
        levels = c("Independence", "Gaussian", "Jittered Gaussian",
                   "ML Gaussian", "t", "Vine", "Jittered Vine")
    ))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1)]
names(colors) <- c("Independence", "Jittered Gaussian", "Gaussian",
                   "Jittered Vine", "Vine", "ML Gaussian", "t")

vars <- c("Pearson", "Spearman", "Kendall", "Mutual information",
          "Biweight midcorrelation", "Distance correlation")

df <- melt(res, id.vars = c("family", "ref"), measure.vars = vars)
df$value <- log10(df$value)

pplt <- ggplot(data = df, aes(x = family, y = value, fill = family)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(size = 0.7, position = position_dodge2(width = 0.22)) +
    geom_line(aes(group = ref), linewidth = 0.018, linetype = "dashed",
              color = "gray", position = position_dodge2(width = 0.22),
              alpha = 0.7) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ variable, nrow = 2, scales = "free_y") +
    ylab(expression(log[10](RMSE))) +
    xlab(NULL) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 13, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13, color = "black"),
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 13, hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 13, color = "black"),
        panel.spacing.x = unit(0.9, "cm")
    )

ggsave(plot = pplt, filename = "Figures/figure_2.pdf", width = 10, height = 9)

eff <- function(f1, f2, m) {
    r1 <- res[res$family == f1, ]
    r2 <- res[res$family == f2, ]
    assertthat::assert_that(all(r1$ref == r2$ref))
    return(effsize::cohen.d(r1[[m]], r2[[m]], paired = TRUE)$estimate)
}

comp <- pairwise_wilcox_test(
    df = res[res$family != "Independence", ],
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
    rename(family1 = V1, family2 = V2, measure = var) %>%
    mutate(pval = sprintf("%.5e", pval),
           padj = sprintf("%.5e", padj),
           stat = sprintf("%.5f", stat)) %>%
    mutate(
        measure = case_match(
            measure,
            "Pearson" ~ "pearson",
            "Kendall" ~ "kendall",
            "Spearman" ~ "spearman",
            "Distance correlation" ~ "dcor",
            "Biweight midcorrelation" ~ "bicor",
            "Mutual information" ~ "MI"),
        eff = sprintf("%.5f", abs(effs))
    ) %>%
    arrange(desc(eff))

write.csv(
    x = comp,
    file = "Tables/pairwise_pvalues.csv",
    row.names = FALSE,
    quote = FALSE
)
