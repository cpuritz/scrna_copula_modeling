library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(grid)
library(ggsci)
set.seed(0)

source("scrna_copula_modeling/09figures/pairwise_wilcox_test.R")

files <- list.files("Results/08margins", full.names = TRUE)
files <- files[grepl("margins_", files)]
res_list <- lapply(files, readRDS)
res <- do.call(rbind, res_list)
res <- res %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarize(pval = mean(pval), .groups = "drop")
res$family2 <- res$family
res$family2[!res$family %in% c("zinbwave", "sparsim")] <- "NB/ZINB"
res$family2[res$family == "zinbwave"] <- "ZINB-WaVE"
res$family2[res$family == "sparsim"] <- "SPARSim"
res_a <- res %>%
	dplyr::group_by(ref, family2) %>%
	dplyr::summarize(pval = mean(pval), .groups = "drop")

colors <- ggsci::pal_npg()(10)[c(7, 8, 9)]
names(colors) <- c("NB/ZINB", "ZINB-WaVE", "SPARSim")

pplt_a <- ggplot(res_a, aes(x = family2, y = pval, fill = family2)) +
	geom_boxplot(outliers = TRUE) +
	scale_fill_manual(values = colors) +
	xlab(NULL) +
	ylab("KS p-value") +
	theme_bw() +
	theme(panel.border = element_blank(),
		  axis.line = element_line(),
		  panel.grid = element_blank(),
		  axis.title.x = element_text(size = 13),
		  axis.title.y = element_text(size = 13, margin = margin(r = 5)),
		  legend.title = element_text(size = 13),
		  axis.text.x = element_text(size = 13, angle = 30, hjust = 1,
		  						   color = "black"),
		  axis.text.y = element_text(size = 13, color = "black"),
		  legend.position = "none")

###########

res_ks <- res %>%
	dplyr::rename(ks_pval = pval) %>%
	dplyr::select(-family2)

files <- list.files("Results/05pca", full.names = TRUE)
files <- files[grepl(".rds", files) & !grepl("ex", files)]
res_ff <- do.call(rbind, lapply(files, readRDS)) %>%
	dplyr::filter(npc == 10) %>%
	dplyr::select(-npc) %>%
	dplyr::group_by(ref, family) %>%
	dplyr::summarise(pval = mean(p), .groups = "drop") %>%
	melt(id.vars = c("family", "ref"), measure.vars = "pval",
		 value.name = "ff_pval")

res <- res_ks %>%
	dplyr::full_join(res_ff, by = c("family", "ref")) %>%
	dplyr::mutate(
		family = factor(dplyr::recode_values(
			family,
			"norm" ~ "Sample Gaussian",
			"vine" ~ "Vine",
			"nmle" ~ "ML Gaussian",
			"t" ~ "t",
			"norm_jitter" ~ "Jittered Gaussian",
			"vine_jitter" ~ "Jittered Vine",
			"ind" ~ "Independence",
			"zinbwave" ~ "ZINB-WaVE",
			"sparsim" ~ "SPARSim"
		), levels = c("Independence", "Sample Gaussian", "Jittered Gaussian",
			  "ML Gaussian", "t", "Vine", "Jittered Vine", "ZINB-WaVE",
			  "SPARSim"))
	) %>%
	dplyr::rename(Model = family)

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8, 9)]
names(colors) <- c("Independence", "Jittered Gaussian", "Sample Gaussian",
				   "Jittered Vine", "Vine", "ML Gaussian", "t",
				   "ZINB-WaVE", "SPARSim")
pplt_b <- ggplot(res, aes(x = ks_pval, y = ff_pval, color = Model)) +
	geom_point(size = 1.6) +
	xlab("KS p-value") +
	ylab("FF p-value") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(panel.border = element_blank(),
		  axis.line = element_line(),
		  panel.grid = element_blank(),
		  axis.title.x = element_text(size = 13),
		  axis.title.y = element_text(size = 13, margin = margin(r = 5)),
		  legend.title = element_text(size = 13),
		  legend.text = element_text(size = 13),
		  axis.text.x = element_text(size = 13, angle = 30, hjust = 1,
		  						   color = "black"),
		  axis.text.y = element_text(size = 13, color = "black"))

get_annot <- function(x) {
	plot_annotation(
		title = x,
		theme = theme(plot.title = element_text(face = 2, size = 20))
	)
}

pplt <- wrap_elements(pplt_a + get_annot("a")) +
	wrap_elements(pplt_b + get_annot("b")) +
	plot_layout(nrow = 1, widths = c(0.42, 0.58)) +
	theme(plot.margin = margin(0, 0, 0, 0))
ggsave(plot = pplt, filename = "Figures/figure_s3.pdf", width = 11, height = 5)