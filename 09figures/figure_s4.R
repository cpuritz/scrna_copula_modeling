library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(ggsignif)
library(ggsci)
set.seed(0)

source("09figures/pairwise_wilcox_test.R")

files <- list.files("Data/References/Subsets")
qc_df <- NULL
for (f in files) {
	sce <- readRDS(paste0("Data/References/Subsets/", f))
	qc_df <- rbind(qc_df, data.frame(
		ref = unlist(strsplit(f, ".rds"))[1],
		countsPerGene = median(rowSums(counts(sce))),
		countsPerCell = median(colSums(counts(sce))),
		totalCounts = sum(counts(sce)),
		nGene = dim(sce)[1],
		nCell = dim(sce)[2],
		pctZero = mean(counts(sce) == 0),
		pctNonzeroPerCell = median(colMeans(counts(sce) > 0))
	))
}

files <- list.files("Results/08margins", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))

##### Pairwise association ####

files <- list.files("Results/04pairwise", full.names = TRUE)
res2 <- do.call(rbind, lapply(files, readRDS)) %>%
	na.omit() %>%
	tidyr::pivot_wider(names_from = stats, values_from = err) %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarize(
        pearson = mean(pearson, na.rm = TRUE),
        spearman = mean(spearman, na.rm = TRUE),
        kendall = mean(kendall, na.rm = TRUE),
        mi = mean(mi, na.rm = TRUE),
        bicor = mean(bicor, na.rm = TRUE),
        dcor = mean(dcor, na.rm = TRUE),
        .groups = "drop"
    )

################
#### PCA ####

files <- list.files("Results/05pca", full.names = TRUE)
files <- files[grepl(".rds", files) & !grepl("ex", files)]
res_list <- lapply(files, readRDS)
res3 <- do.call(rbind, res_list) %>%
	dplyr::filter(npc == 10) %>%
	dplyr::select(-npc) %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarise(p = -mean(p), .groups = "drop")

################

files <- list.files("Results/06wgcna", full.names = TRUE)
res4 <- do.call(rbind, lapply(files, readRDS)) %>%
	na.omit() %>%
	group_by(family, ref, trial) %>%
	summarize(Z = mean(Z, na.rm = TRUE), .groups = "drop") %>%
	group_by(family, ref) %>%
	summarize(Z = -mean(Z), .groups = "drop")

################

files <- list.files("Results/08margins", full.names = TRUE)
res8 <- do.call(rbind, lapply(files, readRDS))
colnames(res8)[3] <- "pctMarginsRejected"

################

res <- dplyr::full_join(res2, res3, by = c("ref", "family")) %>%
	dplyr::full_join(res4, by = c("ref", "family")) %>%
	dplyr::full_join(res8, by = c("ref", "family")) %>%
	dplyr::full_join(qc_df, by = "ref")

stats <- c("pearson", "spearman", "kendall", "mi", "bicor", "dcor", "p", "Z")
qc_vars <- c("countsPerGene", "countsPerCell", "totalCounts", "nGene", "nCell",
			 "pctZero", "pctNonzeroPerCell", "pctMarginsRejected")
fams <- c("norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t")
combs <- expand.grid(qc_vars, stats, fams)
cors <- c()
for (i in seq_len(dim(combs)[1])) {
	m <- res[res$family == combs[i, 3], c(combs[i, 1], combs[i, 2])]
	m <- na.omit(m)
	if (dim(m)[1] > 0) {
		cors <- c(cors, cor(m)[1, 2])
	} else {
		cors <- c(cors, NA)
	}
}
qc_corr <- cbind(combs, cors)
colnames(qc_corr) <- c("qc", "Statistic", "Copula", "cor")

qc_corr <- qc_corr %>%
	dplyr::mutate(
		Copula = recode_values(
			Copula,
			"norm" ~ "Gaussian",
			"vine" ~ "Vine",
			"norm_jitter" ~ "Jittered Gaussian",
			"vine_jitter" ~ "Jittered Vine",
			"nmle" ~ "ML Gaussian",
			"t" ~ "t"
		),
		Statistic = recode_values(
			Statistic,
			"pearson" ~ "Pearson",
			"spearman" ~ "Spearman",
			"kendall" ~ "Kendall",
			"mi" ~ "MI",
			"bicor" ~ "Bicor",
			"dcor" ~ "dCor",
			"p" ~ "-p",
			"Z" ~ "-Z"
		)
	) %>%
	na.omit()

pplt <- ggplot(data = qc_corr, aes(x = qc, y = cor)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(aes(color = Statistic, shape = Copula), size = 3.25,
			   position = position_dodge2(width = 0.35)) +
	ylab("Pearson correlation") +
	xlab(NULL) +
	ylim(c(-1, 1)) +
	theme_bw() +
	geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
	theme(
		axis.title = element_text(size = 13),
		axis.text.x = element_text(size = 13, color = "black", angle = 30,
								   hjust = 1, vjust = 1),
		axis.text.y = element_text(size = 12, color = "black"),
		legend.text = element_text(size = 13),
		legend.title = element_text(size = 13)
	)
ggsave(plot = pplt, filename = "Figures/figure_s4.pdf", width = 13, height = 6.5)
