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

files_emp <- files[!grepl("par_", files) & grepl(".rds", files)]
res_emp <- do.call(rbind, lapply(files_emp, readRDS))
rownames(res_emp) <- NULL

files_par <- files[grepl("par_", files) & grepl(".rds", files)]
res_par <- do.call(rbind, lapply(files_par, readRDS))
rownames(res_par) <- NULL
res_par$family <- factor(paste(res_par$family, "par", sep = "_"))

res <- rbind(res_emp, res_par)
res$family <- factor(res$family)

res <- res %>%
	tidyr::pivot_wider(names_from = stats, values_from = err) %>%
	dplyr::group_by(ref, family, trial) %>%
	dplyr::summarize(
		ngene = unique(ngene),
		ncell = unique(ncell),
		pearson = mean(pearson),
		spearman = mean(spearman),
		# kendall = mean(kendall),
		mi = mean(mi),
		bicor = mean(bicor),
		# dcor = mean(dcor),
		.groups = "drop"
	) %>%
	na.omit() %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarize(
		ngene = unique(ngene),
		ncell = unique(ncell),
		Pearson = mean(pearson),
		Spearman = mean(spearman),
		# Kendall = mean(kendall),
		"Mutual Information" = mean(mi),
		Bicor = mean(bicor),
		# dCor = mean(dcor),
		.groups = "drop"
	) %>%
	dplyr::mutate(family = as.character(family)) %>%
	dplyr::mutate(family = recode_values(
		family,
		"ind" ~ "Independence",
		"ind_par" ~ "Independence (parametric)",
		"norm" ~ "Gaussian",
		"norm_par" ~ "Gaussian (parametric)",
		"norm_jitter" ~ "Jittered Gaussian",
		"norm_jitter_par" ~ "Jittered Gaussian (parametric)",
		"vine_jitter" ~ "Jittered Vine",
		"vine_jitter_par" ~ "Jittered Vine (parametric)",
	)) %>%
	dplyr::mutate(family = factor(
		family,
		levels = c("Independence", "Independence (parametric)",
				   "Gaussian", "Gaussian (parametric)",
				   "Jittered Gaussian", "Jittered Gaussian (parametric)",
				   "Jittered Vine", "Jittered Vine (parametric)")
	))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 1, 4, 6, 8)]
names(colors) <- c("Independence", "Jittered Gaussian", "Gaussian", "Jittered Vine",
				   "Independence (parametric)", "Jittered Gaussian (parametric)",
				   "Gaussian (parametric)", "Jittered Vine (parametric)")

vars <- c("Pearson", "Spearman", "Mutual Information", "Bicor")#, "Kendall", "dCor")

df <- melt(res, id.vars = c("family", "ref"), measure.vars = vars)
df$value <- log10(df$value)

pplt <- ggplot(data = df, aes(x = family, y = value, fill = family)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(size = 0.7, position = position_dodge2(width = 0.22)) +
	scale_fill_manual(values = colors) +
	facet_wrap(~ variable, nrow = 2, scales = "free_y") +
	ylab(expression(log[10](RMSE))) +
	xlab(NULL) +
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

#ggsave(plot = pplt, filename = "Figures/figure_2.pdf", width = 12, height = 7)
# 
# eff <- function(f1, f2, m) {
#     r1 <- res[res$family == f1, ]
#     r2 <- res[res$family == f2, ]
#     assertthat::assert_that(all(r1$ref == r2$ref))
#     return(effsize::cohen.d(r1[[m]], r2[[m]], paired = TRUE)$estimate)
# }
# 
# comp <- pairwise_wilcox_test(
#     df = res[res$family != "Independence", ],
#     cols = vars,
#     group_var = "family",
#     block_var = "ref",
#     var_name = "var",
#     adjust_method = "fdr"
# )
# 
# # Swap families so that V1 always indicates the worse performer
# for (i in seq_len(dim(comp)[1])) {
#     if (comp[i, "stat"] < 0) {
#         comp[i, "stat"] <- -comp[i, "stat"]
#         f1 <- comp[i, "V1"]
#         comp[i, "V1"] <- comp[i, "V2"]
#         comp[i, "V2"] <- f1
#     }
# }
# 
# effs <- apply(comp, 1, function(r) { eff(r["V1"], r["V2"], r["var"]) })
# 
# comp <- comp %>%
#     rename(family1 = V1, family2 = V2, measure = var) %>%
#     mutate(pval = sprintf("%.5e", pval),
#            padj = sprintf("%.5e", padj),
#            stat = sprintf("%.5f", stat)) %>%
#     mutate(
#         measure = recode_value(
#             measure,
#             "Pearson" ~ "pearson",
#             "Kendall" ~ "kendall",
#             "Spearman" ~ "spearman",
#             "Distance correlation" ~ "dcor",
#             "Biweight midcorrelation" ~ "bicor",
#             "Mutual information" ~ "MI"),
#         eff = sprintf("%.5f", abs(effs))
#     ) %>%
#     arrange(desc(eff))
# 
# write.csv(
#     x = comp,
#     file = "Tables/pairwise_pvalues.csv",
#     row.names = FALSE,
#     quote = FALSE
# )
