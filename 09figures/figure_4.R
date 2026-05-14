library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(grid)
library(ggsci)
set.seed(0)

source("scrna_copula_modeling/09figures/pairwise_wilcox_test.R")

files <- list.files("Results/08margins", full.names = TRUE)
res_ls <- do.call(rbind, lapply(files[grepl("libsize_", files)], readRDS))
colnames(res_ls)[3] <- "ls"
res_zp <- do.call(rbind, lapply(files[grepl("zeroprop_", files)], readRDS))
colnames(res_zp)[3] <- "zp"
res <- res_ls %>%
	dplyr::full_join(res_zp, by = c("family", "ref")) %>%
	dplyr::mutate(family = recode_values(
		family,
		"ind" ~ "Independence",
		"norm" ~ "Sample Gaussian",
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
		levels = c("Independence", "Sample Gaussian", "Jittered Gaussian",
				   "ML Gaussian", "t", "Vine", "Jittered Vine",
				   "ZINB-WaVE", "SPARSim")
	))

colors <- ggsci::pal_npg()(10)[c(7, 2, 5, 3, 4, 6, 1, 8, 9)]
names(colors) <- c("Independence", "Jittered Gaussian", "Sample Gaussian",
				   "Jittered Vine", "Vine", "ML Gaussian", "t",
				   "ZINB-WaVE", "SPARSim")

res$label <- "Library size"
pplt_a <- ggplot(res, aes(x = family, y = ls, fill = family)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(size = 1.6, position = position_dodge2(width = 0.24)) +
	geom_line(aes(group = ref), linewidth = 0.10, linetype = "solid",
			  color = "gray", position = position_dodge2(width = 0.24),
			  alpha = 0.3) +
	scale_fill_manual(values = colors) +
	xlab(NULL) +
	ylab("KS p-value") +
	facet_wrap(~label) +
	theme_bw() +
	theme(panel.border = element_blank(),
		  axis.line = element_line(),
		  panel.grid = element_blank(),
		  axis.title.x = element_text(size = 13),
		  axis.title.y = element_text(size = 13, margin = margin(r = 5)),
		  legend.title = element_text(size = 13),
		  axis.text.x = element_text(size = 13, angle = 40, hjust = 1,
		  						   color = "black"),
		  axis.text.y = element_text(size = 13, color = "black"),
		  legend.position = "none",
		  strip.text = element_text(size = 13, color = "black"))
res$label <- "Zero proportion"
pplt_b <- ggplot(res, aes(x = family, y = zp, fill = family)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(size = 1.6, position = position_dodge2(width = 0.24)) +
	geom_line(aes(group = ref), linewidth = 0.10, linetype = "solid",
			  color = "gray", position = position_dodge2(width = 0.24),
			  alpha = 0.3) +
	scale_fill_manual(values = colors) +
	xlab(NULL) +
	ylab("KS p-value") +
	facet_wrap(~label) +
	theme_bw() +
	theme(panel.border = element_blank(),
		  axis.line = element_line(),
		  panel.grid = element_blank(),
		  axis.title.x = element_text(size = 13),
		  axis.title.y = element_text(size = 13, margin = margin(r = 5)),
		  legend.title = element_text(size = 13),
		  axis.text.x = element_text(size = 13, angle = 40, hjust = 1,
		  						   color = "black"),
		  axis.text.y = element_text(size = 13, color = "black"),
		  legend.position = "none",
		  strip.text = element_text(size = 13, color = "black"))
get_annot <- function(x) {
	plot_annotation(
		title = x,
		theme = theme(plot.title = element_text(face = 2, size = 20))
	)
}

pplt <- wrap_elements(pplt_a + get_annot("a")) +
	wrap_elements(pplt_b + get_annot("b")) +
	plot_layout(nrow = 1, widths = c(0.5, 0.5)) +
	theme(plot.margin = margin(0, 0, 0, 0))

ggsave(plot = pplt, filename = "Figures/figure_4.pdf", width = 10, height = 6)


comp <- pairwise_wilcox_test(
	df = res,
	cols = c("ls", "zp"),
	group_var = "family",
	block_var = "ref",
	var_name = "var",
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

eff <- function(f1, f2, m) {
	r1 <- res[res$family == f1, ]
	r2 <- res[res$family == f2, ]
	assertthat::assert_that(all(r1$ref == r2$ref))
	return(effsize::cohen.d(r1[[m]], r2[[m]], paired = TRUE)$estimate)
}
effs <- apply(comp, 1, function(r) { eff(r["V1"], r["V2"], r["var"]) })

comp <- comp %>%
	dplyr::rename(family1 = V1, family2 = V2, measure = var) %>%
	dplyr::mutate(eff = abs(effs)) %>%
	dplyr::arrange(dplyr::desc(eff)) %>%
	dplyr::mutate(
		eff = sprintf("%.5f", eff),
		measure = dplyr::recode_values(
			measure,
			"ls" ~ "Library size",
			"zp" ~ "Zero proportion"
		)
	)