library(dplyr)
library(parallel)
library(SingleCellExperiment)
set.seed(0, kind = "L'Ecuyer-CMRG")

setwd("/projects/b1167/connor_puritz/copgen")

args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1])
ref_ix <- as.integer(args[2])

refs <- list.dirs("Results/03samples", recursive = FALSE, full.names = FALSE)
ref <- refs[ref_ix]

message("Testing library size for ", ref)
ntrial <- 50L
nsample <- 20L

sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
N <- floor(dim(sce)[2] / 2)
ngene <- dim(sce)[1]

families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t", "zinbwave", "sparsim")
sims <- expand.grid(family = families, trial = seq(ntrial))

pct_reject <- c()
for (i in seq_len(dim(sims)[1])) {
	message(i, "/", dim(sims)[1])
	family <- sims[i, "family"]
	trial <- sims[i, "trial"]
	sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
	lib_test <- colSums(counts(sce_test))
	
	pvals <- mclapply(seq_len(nsample), function(j) {
		fname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, "_par.rds")
		if (!file.exists(fname)) {
			return(NULL)
		}
		sce_sim <- readRDS(fname)
		lib_sim <- colSums(counts(sce_sim))
		return(ks.test(lib_test, lib_sim, exact = TRUE)$p)
	}, mc.cores = cores)
	padj <- p.adjust(unlist(pvals), method = "fdr")
	pct_reject[i] <- mean(padj <= 0.05)
}
sims$pct_reject <- pct_reject
sims$ref <- ref
sims <- sims %>%
	dplyr::group_by(family, ref) %>%
	dplyr::summarize(pct_reject = mean(pct_reject), .groups = "drop")

saveRDS(sims, file = paste0("Results/08margins/libsize_", ref_ix, ".rds"))
message("Done!")
