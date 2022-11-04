#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}
if (!exists("aoi_seurat")) {
    source("./plotting_utils.R")
    user_id <- args[3]
    setwd(file.path("temp", user_id))
    print(getwd())
    aoi_seurat <- readRDS("data/aoi_seurat.rds")
    sample_info <- suppressWarnings(read.csv(args[1],
            header = TRUE))
    sample_info <- sample_info[, -1]
    param <- colnames(sample_info)[2]
    goi <- read.table(args[2], sep = ",", header = FALSE)[, 1]
    #print(goi)
    output_dir <- "./GEC/"
}

# "./Gene_Expression_Comparison/"

gene_exp_comp_analysis <- function(aoi_seurat, goi,
                    param = "DISEASE", output_dir = "./GEC/") {
  # align format of each gene
    goi <- toupper(goi)
    goi <- gsub("\\-|\\_", "", goi)
    # check whether all goi in seurat features
    goi <- intersect(goi, rownames(aoi_seurat))
    feature_violin_folders <- list.files(output_dir,
            pattern = "featureDimPlot.png|ViolinComparePlot.png")
    # featureDim and ViolinCompare
    for (i in 1:length(goi)) {
        if (length(grep(goi[i], feature_violin_folders)) == 0) {
            p <- featureDim_plot(aoi_seurat = aoi_seurat, goi[i],
                            reduction = "umap", param = param)
            p_title <- paste0(output_dir, goi[i], "_",
                paste(param, collapse = "_"),
                "_featureDimPlot.png")
            ggsave(filename = p_title,
                    plot = p, width = 12, height = 12, units = "in")

            p <- ViolinComparePlot(aoi_seurat = aoi_seurat, goi[i],
                                reduction = "umap", param = param)
            p_title <- paste0(output_dir, goi[i],
                "_",
                paste(param, collapse = "_"),
                "_ViolinComparePlot.png")
            ggsave(filename = p_title,
                    plot = p, width = 12, height = 12, units = "in")
        }
    }
    # fraction table and fraction plot
    all_folders <- list.files(output_dir)
    frac_table_path <- paste0(output_dir, param,
            "_", paste(unlist(goi), collapse = "_"),
            "_fractionTable.csv")
    if (file.exists(frac_table_path)) {
        cat("Previous Fraction Table Calculation has already been stored \n")
        frac_table <- read.csv(frac_table_path)
    } else {
        frac_table <- fraction_param_table(aoi_seurat,
                                    threshold = 0.5,
                                    param = "macro_celltype", #will not change
                                    goi = goi)
        write.csv(frac_table, frac_table_path)
    }
    ps <- fraction_plot(frac_table = frac_table, goi = goi)
    p_title <- paste0(output_dir,
                paste(unlist(goi), collapse = "_"), "_fractionPlot.png")
    ggsave(filename = p_title,
            plot = wrap_plots(ps), width = 10, height = 8, units = "in")
    cat("Successfully run Gene Expression Comparison Analysis \n")
    # return(NULL)
}

suppressWarnings(gene_exp_comp_analysis(aoi_seurat,
                    param = param,
                    goi = goi,
                    output_dir = output_dir))