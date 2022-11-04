#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied \n", call. = FALSE)
}

if (!exists("aoi_seurat")) {
    source("./meta_utils.R")
    source("./plotting_utils.R")
    base_dir <- "db_data"
    all_folders <- list.files(base_dir, pattern = "harmony_seurat.rds",
                full.names = TRUE)
    all_folders <- file.path(getwd(), all_folders)
    user_id <- args[2]
    suppressWarnings(dir.create(file.path("temp", user_id), recursive = TRUE))
    # sample_info.txt
    sample_info <- suppressWarnings(read.table(file.path("temp",
            user_id, args[1]),
            header = TRUE))
    output_dir <- "./Step1/"
    setwd(file.path("temp", user_id))
}

acc_aoi_seurat <- function(all_folders, sample_info, threshold = 30000, user_id) {
    aoi_acc <- unique(sample_info[, "accession"])
    aoi_groups <- sample_info[, "groups"]
    ps <- list()
    threshold <- threshold / nrow(sample_info)
    for (i in 1:length(aoi_acc)) {
        ps[[i]] <- locate_seurat_per_acc(aoi_acc[i],
                    seurat_folders = all_folders,
                    group_info = aoi_groups[i],
                    threshold)
    }
    cat("Successfully Retrieve All Accession Seurat \n")
    ps <- ps[!sapply(ps, is.null)]
    seurat <- merge(ps[[1]], ps[-1])
    seurat <- harmony_merge(seurat)
    seurat <- seurat %>%
            FindNeighbors(reduction = "harmony", dims = 1:20) %>%
            RunUMAP(reduction = "harmony", dims = 1:20)
    # clean meta data
    meta <- seurat@meta.data
    meta <- meta_alignment(meta)
    seurat@meta.data <- meta
    suppressWarnings(dir.create("data", recursive = TRUE))
    readr::write_rds(seurat, "data/aoi_seurat.rds")
    aoi_table <- ct_disease_acc_table(seurat, 0)
    write.csv(aoi_table, "data/celltype_disease_acctable.csv")
    cat("Successfully generate AOI seurat \n")
    return(seurat)
}


plot_aoi_seurat <- function(output_dir = "Step1/", aoi_seurat, user_id) {
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
    disease_table <- table(aoi_seurat$DISEASE, aoi_seurat$orig.ident)
    write.csv(disease_table, "data/disease_acctable.csv")
    tissue_table <- table(aoi_seurat$TISSUE, aoi_seurat$orig.ident)
    write.csv(tissue_table, "data/tissue_acctable.csv")
    gender_table <- table(aoi_seurat$gender, aoi_seurat$orig.ident)
    write.csv(gender_table, "data/gender_acctable.csv")
    aging_table <- table(aoi_seurat$age, aoi_seurat$orig.ident)
    write.csv(gender_table, "data/age_acctable.csv")

    meta <- aoi_seurat@meta.data
    meta_param <- c("DISEASE", "TISSUE", "aging", "gender",
                "gse_alias", "macro_celltype")
    meta_param <- meta_param[meta_param %in% colnames(meta)]
    if_meta_dim <- metadata_dimplot(output_dir = output_dir,
        aoi_seurat = aoi_seurat, dimplot_param = meta_param)

    cat("Plotting Celltype Composition Bar plot ", "\n")
    compo_plot(aoi_seurat = aoi_seurat,
                    output_dir = output_dir)
    if (!file.exists(
            paste0(output_dir,
            paste(unique(aoi_seurat$groups), collapse = "_"),
            "_celltype_marker_table.csv"))) {
        cat("Calculating Top Marker and Plotting Violin and Heatmap ", "\n")
        deg_table <- top_markers_violin_heatmap(output_dir = output_dir,
                        aoi_seurat = aoi_seurat,
                        top_num = 10,
                        deg_table = FALSE)
        top_markers <- deg_table %>% group_by(cluster) %>%
                    slice_max(avg_log2FC, n = 5)
        p <- DotPlot(aoi_seurat,
                    group.by = "macro_celltype",
                    features = unique(top_markers$gene)) +
                    RotatedAxis()
        p_title <- paste0(output_dir, "celltype_TopMarkers_Dotplot.png")
        ggsave(filename = p_title,
            plot = p, width = 18, height = 12, units = "in")
    }
}


groupwise_deg_per_cluster <- function(aoi_seurat,
                    all_deg = FALSE,
                    celltype_options, output_dir = "Step1/") {
    deg_list <- list()
    Idents(aoi_seurat) <- aoi_seurat$groups
    for (i in 1:length(celltype_options)) {
        cat("Current celltyep:", celltype_options[i], "\n")
        tmp <- subset(aoi_seurat, macro_celltype == celltype_options[i])
        tmp <- tmp %>%
                FindNeighbors(reduction = "harmony", dims = 1:20)
        deg_table <- FindAllMarkers(tmp,
            logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
            min.diff.pct = 0,
            only.pos = FALSE,
            return.thresh = 0.05,
            assay = "RNA")
        if (nrow(deg_table) == 0) {
            next
        }
        deg_table <- subset(deg_table,
                p_val_adj < 0.05 & abs(avg_log2FC) >= 0.5)
        if (nrow(deg_table) == 0) {
            next
        }
        if (all_deg == FALSE) {
            sig_genes <- deg_table %>%
                group_by(cluster) %>% slice_max(abs(avg_log2FC), n = 20)
        } else {
            sig_genes <- deg_table
        }
        p <- DoHeatmap(tmp,
                    features = sig_genes$gene,
                    group.by = "groups",
                    assay = "RNA") +
                    scale_fill_gradientn(colors = c("blue", "white", "red"))
        p_title <- paste0(output_dir,
                "sig_", celltype_options[i], "_DEG_heatmap.png")
        ggsave(filename = p_title,
            plot = p, width = 28, height = 24, units = "in")
        deg_table$celltype <- celltype_options[i]
        deg_list[[i]] <- deg_table
    }
    write.csv(deg_table, "data/sig_groupwise_DEG_table.csv")
}

if ("data/aoi_seurat.rds" %!in% list.files(recursive = TRUE)) {
    aoi_seurat <- acc_aoi_seurat(all_folders,
            sample_info, threshold = 30000, user_id = user_id)
} else{
    aoi_seurat <- readRDS("data/aoi_seurat.rds")
}

print(getwd())
plot_aoi_seurat(output_dir = "Step1/", aoi_seurat, user_id = user_id)
celltype_options <- unique(aoi_seurat$macro_celltype)
groupwise_deg_per_cluster(aoi_seurat, all_deg = FALSE,
                    celltype_options = unique(aoi_seurat$macro_celltype))

cat("Successfully finish Step1 preparation \n")
