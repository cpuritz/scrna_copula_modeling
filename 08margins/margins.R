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

# The models for the different copula families have the same margins, so the only
# differences are due to sampling variability. So we'll only test one family here
# for simplicity.
families <- c("ind",  "zinbwave", "sparsim")
sims <- expand.grid(family = families, trial = seq(ntrial))

pvals_all <- c()
for (i in seq_len(dim(sims)[1])) {
    message(i, "/", dim(sims)[1])
	family <- sims[i, "family"]
	trial <- sims[i, "trial"]
	sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
    X_test <- counts(sce_test)
	
	pvals <- mclapply(seq_len(nsample), function(j) {
		fname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, "_par.rds")
		if (!file.exists(fname)) {
			return(NULL)
		}
		sce_sim <- readRDS(fname)
        X_sim <- counts(sce_sim)
        # Average p-value over genes
		mean(sapply(seq_len(ngene), function(k) {
			ks.test(X_sim[k, ], X_test[k, ], exact = TRUE)$p
		}))
	}, mc.cores = cores)
    # Average over samples
    pvals_all[i] <- mean(unlist(pvals))
}
sims$pval <- pvals_all
sims$ref <- ref
# Average over train-test splits
sims <- sims %>%
	dplyr::group_by(family, ref) %>%
	dplyr::summarize(pval = mean(pval), .groups = "drop")

saveRDS(sims, file = paste0("Results/08margins/margins_", ref_ix, ".rds"))
