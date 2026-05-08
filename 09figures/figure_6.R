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

files <- list.files("Results/07dist", full.names = TRUE)
files <- files[grepl(".rds", files)]
res <- do.call(rbind, lapply(files, readRDS))
rownames(res) <- NULL

res <- res %>%
	tidyr::pivot_wider(names_from = stat, values_from = ks) %>%
	dplyr::group_by(ref, family, trial) %>%
	dplyr::summarize(
		ngene = unique(ngene),
		ncell = unique(ncell),
		pearson = mean(pearson),
		spearman = mean(spearman),
		bicor = mean(bicor),
		.groups = "drop"
	) %>%
	na.omit() %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarize(
		ngene = unique(ngene),
		ncell = unique(ncell),
		Pearson = mean(pearson),
		Spearman = mean(spearman),
		Bicor = mean(bicor),
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
		"zinbwave" ~ "zinbwave",
		.default = family
	)) %>%
	dplyr::mutate(family = factor(
		family,
		levels = c("Independence", "Gaussian", "Jittered Gaussian",
				   "ML Gaussian", "t", "Vine", "Jittered Vine", "zinbwave")
	))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8)]
names(colors) <- c("Independence", "Jittered Gaussian", "Gaussian",
				   "Jittered Vine", "Vine", "ML Gaussian", "t", "zinbwave")

vars <- c("Pearson", "Spearman", "Bicor")

df <- melt(res, id.vars = c("family", "ref"), measure.vars = vars)

pplt <- ggplot(data = df, aes(x = family, y = value, fill = family)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(size = 0.7, position = position_dodge2(width = 0.22)) +
	geom_line(aes(group = ref), linewidth = 0.018, linetype = "dashed",
			  color = "gray", position = position_dodge2(width = 0.22),
			  alpha = 0.7) +
	scale_fill_manual(values = colors) +
	facet_wrap(~ variable, nrow = 2, scales = "free_y") +
	ylab("KS") +
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