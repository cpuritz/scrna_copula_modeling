library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(grid)
library(ggsci)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper/")
source("scrna_copula_modeling/07figures/pairwise_wilcox_test.R")

files <- list.files("Results/08margins", full.names = TRUE)
res_list <- lapply(files, readRDS)
res <- do.call(rbind, res_list)
res$family <- as.character(res$family)
res$family[!res$family %in% c("zinbwave", "sparsim")] <- "NB/ZINB"
res$family[res$family == "zinbwave"] <- "ZINB-WaVE"
res$family[res$family == "sparsim"] <- "SPARSim"
# Average pct_reject over all copula families for each dataset
res <- res %>%
	dplyr::group_by(family, ref) %>%
	dplyr::summarize(pct_reject = mean(pct_reject), .groups = "drop")

colors <- ggsci::pal_npg()(10)[c(7, 2, 5)]
names(colors) <- c("NB/ZINB", "ZINB-WaVE", "SPARSim")

pplt <- ggplot(res, aes(x = family, y = pct_reject, fill = family)) +
	geom_boxplot(outliers = FALSE) +
	geom_point(size = 1.6, position = position_dodge2(width = 0.24)) +
	scale_fill_manual(values = colors) +
	xlab(NULL) +
	ylab("Fraction rejected") +
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
		  legend.position = "none",
		  plot.margin = margin(5.5, 20, 5.5, 5.5))

ggsave(plot = pplt, filename = "Figures/margins.pdf", width = 8, height = 5)