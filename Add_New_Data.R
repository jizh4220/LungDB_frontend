annotate_celltype <- function(seurat, ref_seurat) {
  seurat_anchors <- FindTransferAnchors(reference = ref_seurat,
                    query = seurat, dims = 1:30)
  predictions <- TransferData(anchorset = seurat_anchors,
                refdata = ref_seurat$celltype,
                dims = 1:30)
  seurat <- AddMetaData(object = seurat,
              metadata = predictions[1], col.name = "celltype")
  return(seurat)
}

cluster_macro_celltype <- function(meta) {
    epithelial_annotation <- c(
                        "Club",
                        "Neuroendocrine",
                        "Ionocyte",
                        "Serous",
                        "Mucous",
                        "Platelet/Megakaryocyte",
                        "Goblet")
    meta$macro_celltype <- meta$celltype
    meta$macro_celltype[grep("Fibroblast|Fibromyocyte",
        meta$macro_celltype, ignore.case = TRUE)] <- "Fibroblast"
    meta$macro_celltype[grep("Natural Killer|NK",
        meta$celltype, ignore.case = TRUE)] <- "NK"
    meta$macro_celltype[grep("Monocyte",
        meta$macro_celltype, ignore.case = TRUE)] <- "Monocytes"
    meta$macro_celltype[grep("Macrophage",
        meta$macro_celltype, ignore.case = TRUE)] <- "Macrophages"
    meta$macro_celltype[grep("Plasma",
        meta$macro_celltype, ignore.case = TRUE)] <- "Plasma"
    meta$macro_celltype[grep("Ciliated",
        meta$macro_celltype, ignore.case = TRUE)] <- "Ciliated"
    meta$macro_celltype[grep("Smooth Muscle|Pericyte",
        meta$macro_celltype, ignore.case = TRUE)] <- "Stromal"
    meta$macro_celltype[grep("Alveolar Epithelial Type 1|AT1",
        meta$macro_celltype, ignore.case = TRUE)] <- "AT1"
        meta$macro_celltype[grep("Alveolar Epithelial Type 2|AT2",
        meta$macro_celltype, ignore.case = TRUE)] <- "AT2"
    meta$macro_celltype[grep("Basal",
        meta$macro_celltype, ignore.case = TRUE)] <- "Basal"
    meta$macro_celltype[grep("Bronchial Vessel|Lymphatic|Capillary|Vein|Artery|Endothelial",
        meta$macro_celltype, ignore.case = TRUE)] <- "Endothelial"
    meta$macro_celltype[meta$macro_celltype %in%
        epithelial_annotation] <- "Epithelial"
    meta$macro_celltype[grep("Epithelial",
        meta$macro_celltype, ignore.case = TRUE)] <- "Epithelial"
    meta$macro_celltype[grep("Mesothelial",
        meta$macro_celltype, ignore.case = TRUE)] <- "Mesothelial"
    meta$macro_celltype[grep(" T|NK/T|T Cell",
        meta$macro_celltype, ignore.case = TRUE)] <- "T"
    meta$macro_celltype[grep(" T|NK/T|B Cell",
        meta$macro_celltype, ignore.case = TRUE)] <- "T"
    meta$macro_celltype[grep("Basophil|Mast",
        meta$macro_celltype,
        ignore.case = TRUE)] <- "Basophil_Mast"
    meta$macro_celltype[grep("DC|Dendritic",
        meta$macro_celltype,
        ignore.case = TRUE)] <- "DC"
    return(meta)
}

ref_seurat <- read_rds("/home/zhang_jiaxuan/Synapse/Ref_harmony_seurat.rds")
seurat <- annotate_celltype(seurat, ref_seurat)
seurat <- cluster_macro_celltype(seurat)

