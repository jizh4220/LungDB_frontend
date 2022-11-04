library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(tidyverse)


# Load Necessary Functions for Our Lung DB scRNA-seq Integration Demo
`%!in%` <- Negate(`%in%`)

library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 51539607552)

routine_QC <- function(seurat) {
  seurat$percent.mt <- PercentageFeatureSet(
                    object = seurat, pattern = "^MT-")
  seurat$log10GenesPerUMI <- log10(
    seurat$nFeature_RNA) / log10(seurat$nCount_RNA)

  seurat <- subset(seurat, percent.mt <= 10)

  seurat <- NormalizeData(seurat,
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

  return(seurat)
}


convert_ENSMBL <- function(counts) {
  # ensmbl_db <- read.csv("/home/zhang_jiaxuan/Lung_DB/utils/ENSMBL_gene_id_name.csv")
  ensmbl_db <- ensmbl_db[, -1]
  ensmbl_db$ENSGID <- gsub("\\.[0-9].*", "", ensmbl_db$ENSGID)
  ensmbl_db <- ensmbl_db[match(rownames(counts), ensmbl_db$ENSGID),]
  idx <- !is.na(ensmbl_db$ENSGID)
  f_ensmbl_db <- ensmbl_db[idx,]
  fcunts <- counts[idx,]
  rownames(fcunts) <- f_ensmbl_db$GENE
  return(fcunts)
}

sub_celltype <- function(seurat, req_celltype) {
  # Idents(seurat) <- seurat$celltype
  cluster_list <- unique(seurat$celltype)
  ctype <- grep(req_celltype, cluster_list, value = TRUE, ignore.case = TRUE)
  if (length(ctype) == 0) {
    print(paste0("No Requested ", req_celltype,
         " Cell Type Found In Current Seurat Object"))
    return(NULL)
  }
  ctype_seurat <- subset(seurat, celltype %in% ctype)
  print(paste0("Successfully Subset Cell Types: ", paste0(ctype, sep = ",")))
  return(ctype_seurat)
}

harmony_merge <- function(merged_seurat) {
  merged_seurat <- routine_QC(merged_seurat)
  library(harmony)
  merged_seurat <- RunHarmony(merged_seurat,
                    "orig.ident", assay.use = "RNA")
  return(merged_seurat)
}


# Fraction_plot necessary
cluster_frac_goi <- function(seurat = seurat,
                            cluster = cluster,
                            goi = goi,
                            param = param,
                            threshold = 0.5) {
    df <- FetchData(seurat,
                    cells = colnames(seurat)
                    [seurat[[param]][, 1] %in% cluster],
                    vars = goi)
    cluster_frac <- c()
    for (i in 1:ncol(df)) {
      tmp <- df[, i]
      exp_cells <- tmp[tmp >= threshold]
      cluster_frac <- c(cluster_frac, length(exp_cells) / length(tmp))
    }
    return(cluster_frac)
}

goi_frac_accession <- function(accession_frac_list,
              seurat, goi = goi,
              param = param, threshold = 0.5) {
    cluster_list <- unique(seurat[[param]])[, 1]
    ps <- lapply(cluster_list, cluster_frac_goi,
            seurat = seurat, goi = goi,
            param = param, threshold = threshold)
    frac_orig <- do.call(rbind, ps)
    frac_orig <- as.data.frame(frac_orig)
    rownames(frac_orig) <- cluster_list
    colnames(frac_orig) <- goi
    frac_orig$celltype <- rownames(frac_orig)
    #frac_orig$groups <- cluster_list
    frac_orig <- frac_orig %>%
                pivot_longer(1:ncol(frac_orig) - 1,
                names_to = "gene", values_to = "expressed_ratio")
    frac_orig$cell_counts <- ncol(seurat)
    frac_orig <- as.data.frame(frac_orig)
    return(frac_orig)
}

fraction_param_table <- function(seurat,
                                threshold, param = param,
                                goi = goi) {
  accession_frac_list <- list()
  orig <- unique(seurat$orig.ident)
  cat("Calculating Fraction Table per GOI and Parameter \n")
  for (i in 1:length(orig)) {
    tmp <- subset(seurat, orig.ident == orig[i])
    accession_frac_list[[i]] <- goi_frac_accession(accession_frac_list,
                                tmp,
                                threshold = threshold,
                                param = param,
                                goi = goi)
    accession_frac_list[[i]]$disease <- unique(tmp$DISEASE)
    accession_frac_list[[i]]$tissue <- unique(tmp$TISSUE)
    accession_frac_list[[i]]$gender <- unique(tmp$gender)
    accession_frac_list[[i]]$aging <- unique(tmp$aging)
  }
  print(accession_frac_list[[i]])
  frac_table <- do.call(rbind, accession_frac_list)
  return(frac_table)
}

fraction_plot <- function(frac_table, goi = goi,
                          reduction = "umap"
                          ) {
    ps <- list()
    for (i in 1:length(goi)) {
        p <- ggplot(frac_table[frac_table$gene == goi[i], ],
                    aes(x = celltype, y = expressed_ratio,
                    fill = disease)) +
                    geom_boxplot() +
                    # facet_wrap(~gene) +
                    labs(x = "Genes of Interest",
                        y = "Fraction of Cells Expressing Target Genes",
                        size = 100) +
                theme(axis.text = element_text(size = 12, angle = 45),
                        legend.key.size = unit(1.5, "cm"),
                        legend.text = element_text(size = 10),
                        plot.title = element_text(size = 30, hjust = 0.5))
        p <- p + labs(title = paste0("Fraction plot of GOI: ",
                                    goi[i]),
                      subtitle = paste0("Total cell counts: ",
                                  sum(unique(frac_table$cell_counts))),
                    )
        ps[[i]] <- p
    }
    return(ps)
}

fetch_by_cluster <- function(seurat_object = seurat, gene_list = gene_list,
                             gene1 = s_gene, clu = clu,
                             diag = "young", method = "pearson") {
  #TODO
  c_names <- gene_list
  #names of gene of interest
  c_names <- c(c_names, gene1)
  #specific cluster and specific Diagnosis: IPF/Control
  tmp <- FetchData(seurat_object,
                cells = colnames(seurat_object)
                [seurat_object$celltype == clu &
                seurat_object$DISEASE == diag],
                vars = c_names)
  #calculate correlation ignoring rows/columns with zeros
  c <- cor(tmp, method = method)
  x <- c[, gene1]
  return(x[names(x) != gene1])
}

single_gene_by_cluster <- function(seurat_object = seurat, gene1 = gene1,
                                   cluster = cluster, method = "pearson") {
  #specific cluster and specific Diagnosis: IPF/Control
  tmp <- FetchData(object = seurat_object,
                cells = colnames(seurat_object)
                [seurat_object$celltype == cluster],
                vars = gene1)
  #calculate correlation ignoring rows/columns with zeros
  return(sum(tmp))
}

#extract correlation matrix with s_gene per Ident
cor_gene_exp_per_Ident <- function(cluster_list, fetch_by_cluster,
                                  seurat_object = seurat, gene_list = gene_list,
                                  gene1 = s_gene, diag = diag,
                                  method = method) {
    #get correlations of each gene of interest in every cell
    m <- lapply(cluster_list, fetch_by_cluster,
                seurat_object = seurat_object,
                gene_list = gene_list,
                gene1 = gene1, diag = diag, method = method)
    whole_m <- do.call(rbind, m)
    rownames(whole_m) <- cluster_list
    return(whole_m)
}

bulk_violin <- function(seurat, goi = goi[1]) {
    p <- VlnPlot(seurat,
                features = goi, group.by = "pseudo_disease",
                pt.size = 0,
                #split.by = "DISEASE",
                same.y.lims = TRUE,
                )
    return(p)
}
bulk_feature <- function(seurat_object = ILD_seurat, goi,
                                 reduction = "umap", param = param, doi = doi) {
  ps <- list()
  #Featureplot with
  for (i in 1:length(doi)) {
    p1 <- FeaturePlot(seurat_object,
            cells = colnames(seurat_object)
            [seurat_object[[param]] == doi[i]],
            features = goi,
            reduction = reduction, pt.size = 1.0, repel = TRUE,
            label = TRUE, label.size = 6) +
            ggtitle(paste0(doi[i]))
    ps[[i]] <- p1
  }
  return(ps)
}

featureDim_plot <- function(aoi_seurat = seurat, goi,
                                 reduction = "umap", param = param) {
  #Featureplot by param
  Idents(aoi_seurat) <- aoi_seurat[[param]]
  p1 <- FeaturePlot(aoi_seurat,
            features = goi,
            reduction = reduction, pt.size = 0.2, repel = TRUE,
            label = TRUE, label.size = 10)
  p2 <- DimPlot(object = aoi_seurat,
              group.by = "macro_celltype",
              reduction = reduction, pt.size = 0.2, repel = TRUE,
              label = TRUE, label.size = 10)
  p <- p1 + p2
  return(p)
}

ViolinComparePlot <- function(aoi_seurat = seurat, goi,
                                 reduction = "umap", param = param) {
  #Violin Plot by param
  p1 <- VlnPlot(aoi_seurat,
                features = goi,
                group.by = param,
                pt.size = 0,
                same.y.lims = TRUE,
                )
  p2 <- VlnPlot(aoi_seurat,
                features = goi,
                group.by = "macro_celltype",
                pt.size = 0,
                same.y.lims = TRUE,
                )
  p <- p1 + p2
  return(p)
}

FractionComparePlot <- function(aoi_seurat = seurat,
                                goi,
                                reduction = "umap",
                                param = param) {
    p1 <- fraction_plot(seurat_object, goi = goi,
                        param = "DISEASE",
                        reduction = "umap"
                        )

}

metadata_dimplot <- function(output_dir, aoi_seurat, dimplot_param) {
  for (i in 1:length(dimplot_param)) {
      if (all(is.na(aoi_seurat[[dimplot_param[i]]][1]))) {
        next
      }
      cat("Plotting Dimplot by metadata:", dimplot_param[i], "\n")
      p <- DimPlot(aoi_seurat,
              group.by = dimplot_param[i],
              reduction = "umap", pt.size = 0.2, repel = TRUE,
              label = TRUE, label.size = 10)
      p_title <- paste0(output_dir, "meta_", dimplot_param[i], "_DimPlot.png")
      ggsave(filename = p_title,
          plot = p, width = 8, height = 8, units = "in")
  }
  return(TRUE)
}


# Celltype Composition
compo_plot <- function(
              aoi_seurat, 
              output_dir = output_dir
              ) {
    celltype_comp <- table(aoi_seurat$macro_celltype, aoi_seurat$groups)
    celltype_comp <- as.data.frame(celltype_comp)
    colnames(celltype_comp) <- c("celltype", "group", "counts")
    # Celltype Composition per group
    p <- ggplot(celltype_comp,
                aes(x = group,
                    y = counts,
                    fill = celltype)) +
        theme_classic(base_size = 15) +
        geom_col(position = "fill", width = 0.5) +
        labs(x = "Group of Interest",
            y = "Celltype Composition",
            title = "Celltype Composition Per Group") +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 20, angle = 45),
              legend.key.size = unit(3, "cm"),
              legend.text = element_text(size = 16),
              plot.title = element_text(size = 40, hjust = 0.5))
    p_title <- paste0(output_dir,
                paste(unique(aoi_seurat$groups),
                collapse = "_"), "_Compoplot.png")
    ggsave(filename = p_title,
          plot = p, width = 18, height = 24, units = "in")

    # Group Composition per Cell type
    p <- ggplot(celltype_comp,
                aes(x = celltype,
                y = counts, fill = group)) +
          theme_classic(base_size = 15) +
          geom_col(position = "fill", width = 0.5) +
          labs(x = "Group of Interest",
              y = "Celltype Composition",
              title = "Celltype Composition Per Group") +
          theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 20, angle = 45),
              legend.key.size = unit(3, "cm"),
              legend.text = element_text(size = 16),
              plot.title = element_text(size = 40, hjust = 0.5))
    p_title <- paste0(output_dir,
            paste(unique(aoi_seurat$groups),
            collapse = "_"), "_per_celltype_Compoplot.png")
    ggsave(filename = p_title,
          plot = p, width = 18, height = 24, units = "in")

    write.csv(celltype_comp, paste0(output_dir,
                paste(unique(aoi_seurat$groups),
                collapse = "_"), "_celltype_CompoTable.csv"))

    # Celltype Composition per Experiment
    celltype_comp <- table(aoi_seurat$macro_celltype, aoi_seurat$gse_alias)
    celltype_comp <- as.data.frame(celltype_comp)
    colnames(celltype_comp) <- c("celltype", "group", "counts")
    p <- ggplot(celltype_comp, aes(x = group, y = counts, fill = celltype)) +
            theme_classic(base_size = 15) +
            geom_col(position = "fill", width = 0.5) +
            labs(x = "Experiment",
                y = "Celltype Composition",
                title = "Celltype Composition Per Group") +
            theme(axis.title = element_text(size = 30),
                axis.text = element_text(size = 20, angle = 45),
                legend.key.size = unit(3, "cm"),
                legend.text = element_text(size = 16),
                plot.title = element_text(size = 40, hjust = 0.5))
    p_title <- paste0(output_dir,
          paste(unique(aoi_seurat$groups), collapse = "_"),
          "_per_experiment_Compoplot.png")
    ggsave(filename = p_title,
        plot = p, width = 36, height = 16, units = "in")
    write.csv(celltype_comp, paste0(output_dir,
                paste(unique(aoi_seurat$groups),
                collapse = "_"), "_per_experiment_CompoTable.csv"))
    # Celltype Composition per Accession
    celltype_comp <- table(aoi_seurat$macro_celltype, aoi_seurat$orig.ident)
    celltype_comp <- as.data.frame(celltype_comp)
    colnames(celltype_comp) <- c("celltype", "group", "counts")
    p <- ggplot(celltype_comp, aes(x = group,
            y = counts, fill = celltype)) +
            theme_classic(base_size = 15) +
            geom_col(position = "fill", width = 0.5) +
            labs(x = "Experiment",
                y = "Celltype Composition",
                title = "Celltype Composition Per Group") +
            theme(axis.title = element_text(size = 30),
                axis.text = element_text(size = 20, angle = 45),
                legend.key.size = unit(3, "cm"),
                legend.text = element_text(size = 16),
                plot.title = element_text(size = 40, hjust = 0.5))
    p_title <- paste0(output_dir,
          paste(unique(aoi_seurat$groups), collapse = "_"),
          "_per_accession_Compoplot.png")
    ggsave(filename = p_title,
        plot = p, width = 36, height = 16, units = "in")
    write.csv(celltype_comp, paste0(output_dir,
                  paste(unique(aoi_seurat$groups),
                  collapse = "_"), "_per_accession_CompoTable.csv"))
    cat("Finish Composition Table and Plot per Accession, Experiment, and Group \n")
    return(TRUE)
}

top_markers_violin_heatmap <- function(output_dir = output_dir,
            aoi_seurat, top_num = 5, deg_table = deg_table) {
    # Marker Genes Representation via Violinplot
    Idents(aoi_seurat) <- aoi_seurat$macro_celltype
    if (deg_table == FALSE) {
      deg_table <- FindAllMarkers(aoi_seurat,
        logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
        min.diff.pct = 0,
        only.pos = FALSE,
        return.thresh = 0.05,
        assay = "RNA")
      # if (deg_table)
      deg_table <- subset(deg_table,
          p_val_adj < 0.05 & abs(avg_log2FC) >= 0.5)
      write.csv(deg_table,
          paste0(output_dir,
          paste(unique(aoi_seurat$groups),
          collapse = "_"),
          "_celltype_marker_table.csv"))
    }

    #Marker Genes Representation via Heatmap
    top_markers <- deg_table %>% group_by(cluster) %>%
        slice_max(avg_log2FC, n = top_num)
    p <- DoHeatmap(aoi_seurat,
                features = top_markers$gene,
                group.by = "macro_celltype",
                assay = "RNA") +
                scale_fill_gradientn(colors = c("blue", "white", "red"))
    p_title <- paste0(output_dir,
          paste(unique(aoi_seurat$groups), collapse = "_"), "_DEG_heatmap.png")
    ggsave(filename = p_title,
        plot = p, width = 28, height = 24, units = "in")
    return(deg_table)
}
