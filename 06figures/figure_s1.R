library(SingleCellExperiment)
library(ggplot2)
library(ggsci)
set.seed(0)
setwd("~/Documents/Graduate School/Copula Paper")

files <- list.files("Data/References/Subsets")
dist <- NULL
for (i in seq_along(files)) {
    sce <- readRDS(paste0("Data/References/Subsets/", files[i]))
    df <- data.frame(
        ngene = dim(sce)[1],
        ncell = dim(sce)[2],
        type = sub(".*-(.*)\\.rds$", "\\1", files[i])
    )
    dist <- rbind(dist, df)
}
dist$type <- factor(
    dist$type,
    levels = c("HVG", "hsa05224", "hsa04915", "hsa04911", "hsa04260"),
    ordered = TRUE
)
colors <- setNames(
    ggsci::pal_npg()(10)[c(2, 6, 3, 5, 1)],
    c("HVG", "hsa04260", "hsa04915", "hsa05224", "hsa04911")
)

pplt <- ggplot(
    data = dist,
    aes(x = ncell, y = ngene, color = type, shape = type, fill = type,
        size = type)
    ) +
    geom_point() +
    xlab("Cells") +
    ylab("Genes") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(
        "HVG" = 16,
        "hsa04260" = 15,
        "hsa04915" = 17,
        "hsa05224" = 18,
        "hsa04911" = 25
    )) +
    scale_size_manual(values = c(
        "HVG" = 3.0,
        "hsa04260" = 2.5,
        "hsa04915" = 2.55,
        "hsa05224" = 3.5,
        "hsa04911" = 2.25
    )) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 13),
          axis.text = element_text(color = "black", size = 13),
          axis.title = element_text(size = 13))
