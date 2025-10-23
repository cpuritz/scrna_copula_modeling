##############################################################################

#' Pairwise Wilcoxon Rank Sum Tests
#'
#' Perform pairwise Wilcoxon rank sum tests with FDR correction.
#'
#' @details The multiple-testing correction is performed over all response
#' variables simultaneously. Exact p-values are computed including in the
#' presence of ties.
#'
#' @param df A data.frame.
#' @param cols Names of response variable columns.
#' @param group_var Name of column indicating the grouping variable.
#' @param block_var Name of optional block factor. If not \code{NULL}, paired
#' tests are performed.
#' @param var_name Name of output variable.
#' @param adjust_method Method for adjusting p-values.
#' @param ... Additional arguments to pass to \code{coin::wilcox_test}.
#' 
#' @returns A data.frame.
pairwise_wilcox_test <- function(df,
                                 cols,
                                 group_var,
                                 block_var = NULL,
                                 var_name = "var",
                                 adjust_method = "fdr",
                                 ...) {
    all_cols <- c(cols, group_var, block_var)
    bad_cols <- all_cols[!all_cols %in% colnames(df)]
    if (length(bad_cols) > 0) {
        stop(bad_cols[1], " is not a column in 'df'.")
    }
    if (!adjust_method %in% stats::p.adjust.methods) {
        stop(adjust_method, " is not a valid method for adjusting p-values.")
    }

    groups <- sort(unique(as.character(df[[group_var]])))
    relevant_comps <- as.data.frame(t(utils::combn(groups, 2)))
    all_comps <- lapply(cols, function(x) {
        stats <- vector(mode = "numeric", length = nrow(relevant_comps))
        p_vals <- vector(mode = "numeric", length = nrow(relevant_comps))
        
        # Loop over all pairs of groups
        for (i in seq(nrow(relevant_comps))) {
            sub <- dplyr::filter(
                df, !!dplyr::sym(group_var) %in% relevant_comps[i, ]
            )
            if (is.null(block_var)) {
                # No block specified, don't run paired tests
                sub <- sub %>%
                    dplyr::select(dplyr::all_of(c(group_var, x))) %>%
                    dplyr::rename(value = !!dplyr::sym(x))
                # Refactor groups in alphabetical order
                group_factored <- droplevels(factor(sub[[group_var]],
                                                    levels = groups))
                sub_df <- data.frame(value = sub$value,
                                     group = group_factored)
                test <- coin::wilcox_test(value ~ group,
                                          data = sub_df,
                                          distribution = "exact",
                                          ...)
            } else {
                # Perform paired tests
                sub <- sub %>%
                    dplyr::select(
                        dplyr::all_of(c(group_var, x, block_var))
                    ) %>%
                    dplyr::rename(value = !!dplyr::sym(x))
                # Refactor groups in alphabetical order
                group_factored <- droplevels(factor(sub[[group_var]],
                                                    levels = groups))
                sub_df <- data.frame(value = sub$value,
                                     group = group_factored,
                                     block = factor(sub[[block_var]]))
                test <- coin::wilcoxsign_test(value ~ group | block,
                                              data = sub_df,
                                              distribution = "exact",
                                              ...)
            }
            stats[i] <- coin::statistic(test)
            p_vals[i] <- coin::pvalue(test)
        }
        return(dplyr::mutate(relevant_comps, !!var_name := x, pval = p_vals,
                             stat = stats))
    })

    # Merge tests and adjust p-values
    all_comps <- all_comps %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(padj = stats::p.adjust(pval, method = adjust_method))

    return(all_comps)
}

##############################################################################

#' Significance bars for boxplots
#'
#' Perform pairwise Wilcoxon rank sum tests with FDR correction and construct
#' significance bars for boxplots.
#'
#' @details The multiple-testing correction is performed over all response
#' variables simultaneously. Exact p-values are computed including in the
#' presence of ties.
#'
#' @param df A data.frame.
#' @param cols Names of response variable columns.
#' @param group_var Name of column indicating the grouping variable.
#' @param block_var Name of optional block factor. If not \code{NULL}, paired
#' tests are performed.
#' @param var_name Name of output variable.
#' @param adjust_method Method for adjusting p-values.
#' @param alpha Filter comparisons with adjusted p-values greater than this
#' threshold. If \code{NULL}, no filtering is done. Default is \code{NULL}.
#' @param annot Format for annotations. If set to \code{"value"}, adjusted
#' p-values are formatted in scientific notation with three digits. If set to
#' \code{"stars"}, standard significance codes are used (*** < 0.001,
#' ** < 0.01, * < 0.05).
#' @param offset How much to offset significance bars from top of boxplots.
#' @param scale How much space to place between stacked significance bars.
#' @param transform Optional data transformation (after testing, before plotting).
#' 
#' @returns A data.frame.
boxplot_signif <- function(df,
                           cols,
                           group_var,
                           block_var = NULL,
                           var_name = "var",
                           adjust_method = "fdr",
                           alpha = NULL,
                           annot = c("value", "stars"),
                           offset = 0.175,
                           scale = 0.375,
                           transform = NULL) {
    # Perform pairwise Wilcoxon rank sum tests
    all_comps <- pairwise_wilcox_test(
        df = df,
        cols = cols,
        group_var = group_var,
        block_var = block_var,
        var_name = var_name,
        adjust_method = adjust_method
    )
    
    # Filter by adjusted p-value
    if (is.null(alpha)) {
        alpha <- 1
    }
    all_comps <- dplyr::filter(all_comps, padj <= alpha)

    annot <- match.arg(annot)
    if (annot == "value") {
        annotations <- format(all_comps$padj, digits = 3, scientific = T)
    } else if (annot == "stars") {
        annotations <- sapply(all_comps$padj, function(x) {
            if (x <= 0.001) {
                return("***")
            } else if (x <= 0.01) {
                return("**")
            } else if (x <= 0.05) {
                return("*")
            } else {
                return("n.s.")
            }
        })
    }
    all_comps <- dplyr::mutate(all_comps, annot = annotations)

    # Add in x-coordinates for significance bars
    groups <- unique(df[[group_var]])
    xmin <- data.frame(V1 = sort(groups), xmin = seq_along(groups))
    xmax <- data.frame(V2 = sort(groups), xmax = seq_along(groups))
    all_comps <- all_comps %>%
        dplyr::left_join(xmin, by = "V1") %>%
        dplyr::left_join(xmax, by = "V2") %>%
        tibble::rownames_to_column(var = "ix")
    
    # Optional transform
    if (!is.null(transform)) {
        df[, cols] <- transform(df[, cols])
    }

    # Add in y-coordinates for significance bars
    if (dim(all_comps)[1] > 0) {
        all_comps$y <- apply(all_comps, 1, function(x) {
            yvals <- df %>%
                dplyr::filter(
                    !!dplyr::sym(group_var) %in% !!x[c("V1", "V2")]
                ) %>%
                dplyr::select(dplyr::all_of(x[[var_name]]))
            return((1 + offset) * max(yvals))
        })
    }

    # Adjust y-coordinates if there are multiple significance bars
    for (i in seq_along(cols)) {
        comps <- dplyr::filter(all_comps, !!dplyr::sym(var_name) == cols[i])
        if (nrow(comps) > 0) {
            diffs <- abs(comps$xmax - comps$xmin)
            ranks <- rank(diffs, ties.method = "first")
            for (j in 1:nrow(comps)) {
                yscale <- 1 + scale * (j - 1)
                all_comps[comps[ranks[j], ]$ix, ]$y <- yscale * max(comps$y)
            }
        }
    }
    return(all_comps)
}

##############################################################################
