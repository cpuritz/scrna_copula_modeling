library(dplyr)
library(parallel)
library(SingleCellExperiment)
set.seed(0, kind = "L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1])
ref_ix <- as.integer(args[2])

refs <- list.dirs("Results/03samples", recursive = FALSE, full.names = FALSE)
ref <- refs[ref_ix]

ntrial <- 50L
nsample <- 20L

sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
N <- floor(dim(sce)[2] / 2)
ngene <- dim(sce)[1]

families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t", "zinbwave", "sparsim")
sims <- expand.grid(family = families, trial = seq(ntrial))

pvals_all <- c()
for (i in seq_len(dim(sims)[1])) {
	message(i, "/", dim(sims)[1])
	family <- sims[i, "family"]
	trial <- sims[i, "trial"]
	sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
    zp_test <- apply(counts(sce_test) == 0, 2, mean)
	
	pvals <- mclapply(seq_len(nsample), function(j) {
		fname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, "_par.rds")
		if (!file.exists(fname)) {
			return(NULL)
		}
		sce_sim <- readRDS(fname)
        zp_sim <- apply(counts(sce_sim) == 0, 2, mean)
		return(ks.test(zp_test, zp_sim, exact = TRUE)$p)
	}, mc.cores = cores)
    pvals_all[i] <- mean(unlist(pvals))
}
sims$pval <- pvals_all
sims$ref <- ref
sims <- sims %>%
	dplyr::group_by(family, ref) %>%
	dplyr::summarize(pval = mean(pval), .groups = "drop")

saveRDS(sims, file = paste0("Results/08margins/zeroprop_", ref_ix, ".rds"))
