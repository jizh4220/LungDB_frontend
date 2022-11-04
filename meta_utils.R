
extract_metadata <- function(meta = meta) {
    #filter metadata columns
    pat1 <- "dis|infection|age|tissue|sex|gender|ENA|sample|donor|origin|"
    pat2 <- "Feature|Count|umi|orig.ident|celltype|gene|gse|mt|percent|pct|mito|groups"
    pat <- paste0(pat1, pat2)
    idx <- grep(pat, colnames(meta), ignore.case = TRUE)
    meta <- meta[, idx]
    idx <- grep("ontology|term|score|CFauc|GOBPauc", colnames(meta), invert = T)
    meta <- meta[, idx]
    pseudo <- grep("pseudo_disease", colnames(meta))
    if(length(pseudo) > 0) {
        pseudo_col <- meta[, pseudo]
        meta[, pseudo] <- NULL
    }
    #nCount/nFeature
    if(length(grep("nCount_RNA|nFeature_RNA",
        colnames(meta), ignore.case = TRUE)) == 0) {
        col_idx <- grep("nUMIs|n_gene", colnames(meta), ignore.case = TRUE)
        colnames(meta)[col_idx] <- "nFeature_RNA"
        col_idx <- grep("n_count", colnames(meta), ignore.case = TRUE)
        colnames(meta)[col_idx] <- "nCount_RNA"
    }
    #adjust DISEASE
    col_idx <- grep("disease", colnames(meta), ignore.case = TRUE)
    #exist both
    if (length(col_idx) > 1) {
        s_idx <- grep("disease", colnames(meta))
        g_idx <- grep("DISEASE", colnames(meta))
        meta$DISEASE[!is.na(meta$disease)] <- meta$disease[!is.na(meta$disease)]
        meta[, s_idx] <- NULL
    } else if (length(col_idx) == 0) {
        meta$DISEASE <- "Normal"
    } else {
        colnames(meta)[col_idx] <- "DISEASE"
    }
    col_idx <- grep("disease", colnames(meta), ignore.case = TRUE)
    meta[, col_idx] <- as.character(meta[, col_idx])
    meta[, col_idx] <- disease_conversion(meta, col_idx)

    #adjust gender
    col_idx <- grep("gender|sex", colnames(meta), ignore.case = TRUE)
    if (length(col_idx) > 1) {
        s_idx <- grep("sex", colnames(meta), ignore.case = TRUE)
        g_idx <- grep("gender", colnames(meta), ignore.case = TRUE)
        meta[, g_idx] <- meta[, s_idx]
        meta[, s_idx] <- NULL
    } else if (length(col_idx) == 0) {
        meta$gender <- NA
    } else {
        colnames(meta)[col_idx] <- "gender"
    }
    col_idx <- grep("gender|sex", colnames(meta), ignore.case = TRUE)
    meta[, col_idx] <- as.character(meta[, col_idx])
    idx <- grep("female|f", meta[, col_idx], ignore.case = TRUE)
    meta[, col_idx][idx] <- "F"
    idx <- grep("mal|m", meta[, col_idx], ignore.case = TRUE)
    meta[, col_idx][idx] <- "M"

    #adjust tissue
    col_idx <- grep("TISSUE", colnames(meta), ignore.case = TRUE)
    if (length(col_idx) > 1) {
        s_idx <- grep("tissue", colnames(meta))
        g_idx <- grep("TISSUE", colnames(meta))
        meta[, g_idx] <- meta[, s_idx]
        meta[, s_idx] <- NULL
    } else if (length(col_idx) == 0) {
        meta$TISSUE <- NA
    } else {
        colnames(meta)[col_idx] <- "TISSUE"
    }

    #adjust age
    col_idx <- grep("age", colnames(meta), ignore.case = TRUE)
    if (length(col_idx) > 1) {
        s_idx <- grep("Age", colnames(meta))
        g_idx <- grep("age", colnames(meta))
        meta[, g_idx] <- meta[, s_idx]
        meta[, s_idx] <- NULL
    } else if (length(col_idx) == 0) {
        meta$age <- NA
    } else {
        colnames(meta)[col_idx] <- "age"
    }
    col_idx <- grep("age", colnames(meta), ignore.case = TRUE)
    meta[, col_idx] <- as.character(meta[, col_idx])
    meta[, col_idx] <- gsub("-year-old human stage", "", meta[, col_idx])
    meta[, col_idx] <- gsub("eighth.*", "80", meta[, col_idx])
    meta[, col_idx] <- gsub("seventh.*", "70", meta[, col_idx])
    meta[, col_idx] <- gsub("sixth.*", "60", meta[, col_idx])
    meta[, col_idx] <- gsub("fifth.*", "50", meta[, col_idx])
    meta[, col_idx] <- gsub("fourth.*", "40", meta[, col_idx])
    meta[, col_idx] <- gsub(" years", "", meta[, col_idx])
    meta[, col_idx] <- gsub("human middle.*", "medium", meta[, col_idx])
    meta[, col_idx] <- gsub("230", "organoid", meta[, col_idx])
    meta[, col_idx] <- gsub("<| y", "", meta[, col_idx])
    meta[, col_idx] <- gsub(" - .*", "", meta[, col_idx])
    meta <- cluster_pseudo_disease(meta)
    return(meta)
}



disease_conversion <- function(meta, col_idx) {
  library(stringr)
  idx <- grep("Normal|Control|Healthy|neg|ipsc|, lower lobe", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Normal"
  meta[, col_idx]["lung" == str_to_lower(meta[, col_idx] )] <- "Normal"
  idx <- grep("non-small cell lung carcinoma|NSC", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "NSCLC"
  idx <- grep("small cell lung carcinoma", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "SCLC"
  idx <- grep("lung adenocarcinoma|LUAD|lung cancer", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "LUAD"
  idx <- grep("COVID|SARS|COV", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "COVID-19"
  idx <- grep("tumor", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Tumor"
  idx <- grep("idiopathic pulmonary fibrosis|IPF", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "IPF"
  idx <- grep("lung fibrosis", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Lung fibrosis"
  idx <- grep(" carcinoma", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Carcinoma"
  idx <- grep("organoid", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "organoid"
  idx <- grep("chronic obstructive pulmonary disease|COPD",
              meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "COPD"
  idx <- grep("Non-specific interstitial pneumonia|Nonspecific interstitial pneumonia|NSIP",
              meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "NSIP"
  idx <- grep("interstitial lung disease|ILD", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "ILD"
  idx <- grep("LUSC", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "ILD"
  idx <- grep("hypersensitivity pneumonitis|HP", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "HP"
  idx <- grep("sacroidosis|sarcoidosis", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Sarcoidosis"
  idx <- grep("scleroderma", meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Scleroderma"
  return(meta[, col_idx])
}

cluster_disease <- function(meta) {
  library(stringr)
  col_idx <- grep("DISEASE", colnames(meta))
  idx <- grep("Normal|Control|Healthy|neg|ipsc|, lower lobe",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Normal"
  meta[, col_idx]["lung" == str_to_lower(meta[, col_idx])] <- "Normal"
  idx <- grep("non-small cell lung carcinoma|NSC",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "NSCLC"
  idx <- grep("small cell lung carcinoma",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "SCLC"
  idx <- grep("lung adenocarcinoma|LUAD",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "LUAD"
  idx <- grep("tumor",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Tumor"
  idx <- grep(" carcinoma",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Carcinoma"
  #COVID
  idx <- grep("COVID|SARS|COV",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "COVID-19"
  #Fibrosis
  idx <- grep("idiopathic pulmonary fibrosis|IPF",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "IPF"
  idx <- grep("chronic obstructive pulmonary disease|COPD",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "COPD"
  idx <- grep("organoid",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Unknown"
  idx <- grep("Non-specific interstitial pneumonia|Nonspecific interstitial pneumonia|NSIP",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "NSIP"
  idx <- grep("interstitial lung disease|ILD",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "ILD"
  idx <- grep("hypersensitivity pneumonitis|HP",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "HP"
  idx <- grep("sacroidosis|sarcoidosis",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Sarcoidosis"
  idx <- grep("Chronic Beryllium Disease",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "CBD"
  idx <- grep("pleuropulmonary blastoma",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "PPB"
  idx <- grep("scleroderma",
        meta[, col_idx], ignore.case = TRUE)
  meta[, col_idx][idx] <- "Scleroderma"
  meta[, col_idx][meta$gse_alias %in% c("GSE150263", "GSE152654")] <- "Normal"
  meta[, col_idx][meta$gse_alias %in% c("GSE183590")] <- "NSCLC"
  return(meta)
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

cluster_aging <- function(meta) {
    meta$aging <- meta$age
    meta$aging[grep("stage|<1|months", meta$aging)] <- "fetal"
    meta$aging[grep("- .* y|17", meta$aging)] <- "kid and teenager"
    meta$aging[meta$aging %in% 55:100] <- "old"
    meta$aging[meta$aging %in% 18:40] <- "adult"
    meta$aging[meta$aging %in% 41:54] <- "medium"
    return(meta)
}

cluster_tissue <- function(meta) {
    meta$TISSUE[grep("PBMC|Peripheral blood", meta$TISSUE)] <- "PBMC"
    meta$TISSUE[grep("peripheral", meta$TISSUE, ignore.case = TRUE)] <- "peripheral lung"
    meta$TISSUE[grep("Proximal lung", meta$TISSUE, ignore.case = TRUE)] <- "proximal lung"
    meta$TISSUE[grep("Whole|lung mesenchymal cells", meta$TISSUE, ignore.case = TRUE)] <- "whole lung"
    meta$TISSUE[grep("explant", meta$TISSUE, ignore.case = TRUE)] <- "lung explant"
    meta$TISSUE[grep("SCLC|cancer|Tumour|PPB|tumor", meta$TISSUE, ignore.case = TRUE)] <- "lung tumor"
    meta$TISSUE[grep("upper", meta$TISSUE)] <- "lung upper lobe"
    meta$TISSUE[grep("lower", meta$TISSUE)] <- "lung lower lobe"
    meta$TISSUE[grep("middle", meta$TISSUE)] <- "lung middle lobe"
    meta$TISSUE[grep("Brochoalvelolar|bronchoalveolar",
            meta$TISSUE, ignore.case = TRUE)] <- "BAL"
    meta$TISSUE[grep("Airway .* distal|distal airway",
            meta$TISSUE, ignore.case = TRUE)] <- "airway distal"
    meta$TISSUE[grep("Airway .* proximal|proximal airway",
            meta$TISSUE, ignore.case = TRUE)] <- "airway proximal"
    meta$TISSUE[meta$TISSUE == "Primary Pulmonary Endothelial"] <- "microvascular"
    meta$TISSUE[grep("Trachea", meta$TISSUE, ignore.case = TRUE)] <- "trachea"
    meta$TISSUE[grep("Pleura", meta$TISSUE, ignore.case = TRUE)] <- "pleura effusion"
    meta$TISSUE[meta$TISSUE %in% c("lung", "Lung")] <- "whole lung"
    meta$TISSUE[grep("MDA",
            meta$TISSUE, ignore.case = TRUE)] <- "CDX Models"
    meta$TISSUE[grep("Adjacent", meta$TISSUE, ignore.case = TRUE)] <- "Adjacent proximal"
    meta$TISSUE[meta$gse_alias %in% "GSE115982"] <- "lung explant"
    meta$TISSUE[meta$gse_alias %in% c("GSE132771", "GSE132914")] <- "whole lung"
    meta$TISSUE[meta$gse_alias == "GSE137026"] <- "lung fibroblast"
    meta$TISSUE[meta$gse_alias == "GSE160915"] <- "PaTu-8902/HUVECs"
    
    meta$TISSUE[meta$gse_alias %in% c("GSE138693", "GSE162045")] <- "PC9"
    meta$TISSUE[meta$gse_alias %in% c("GSE161934", "GSE178519")] <- "organoid"
    meta$TISSUE[meta$gse_alias %in% c("GSE154870", "GSE154869",
                    "GSE142286", "GSE126908",
                    "GSE126906", "GSE117450", "GSE183590")] <- "H2228, H1975, A549, H838, HCC827"
    meta$TISSUE[meta$orig.ident == "GSM5032894"] <- "lung tumor"
    meta$TISSUE[meta$orig.ident == "GSM4558305"] <- "liver"
    meta$TISSUE[meta$gse_alias == "GSE178360"] <- "airway distal"
    meta$TISSUE[meta$gse_alias == "GSE103918"] <- "hPSC"
    meta$TISSUE[grep("ipsc",
            meta$TISSUE, ignore.case = TRUE)] <- "iPSC"
    meta$TISSUE[grep("cell|P0|CF|COPD|COVID-19|Tumour|Normal|PF|LAM|measles virus|LUAD|mock",
            meta$TISSUE, ignore.case = TRUE)] <- "whole lung"
    meta$TISSUE[meta$gse_alias %in% c("GSE150263", "GSE152654")]<- "iPSC"
    return(meta)
}

cluster_pseudo_disease <- function(meta) {
    meta$pseudo_disease <- meta$DISEASE
    meta$pseudo_disease[grep("COPD|IPF|Sarcoidosis|Scleroderma|cHP|fibrosis|CF|Chronic Beryllium Disease|PF|ILD",
        meta$pseudo_disease)] <- "lung fibrosis"
    meta$pseudo_disease[grep("HP|LAM|Asthma|measles virus|NSIP",
        meta$pseudo_disease)] <- "lung inflammation"
    meta$pseudo_disease[grep("PAIVS|IPAH",
        meta$pseudo_disease)] <- "lung pulmonary artery"
    meta$pseudo_disease[grep("H9|Smoker|transplantation",
        meta$pseudo_disease)] <- "Unknown"
    meta$pseudo_disease[grep("cancer|carcinoma|tumor|LUAD|SCLC|PPB|TNBC|NSC",
        meta$pseudo_disease, ignore.case = TRUE)] <- "lung cancer/tumor"
    return(meta)
}

meta_alignment <- function(meta) {
    meta <- extract_metadata(meta)
    meta <- cluster_aging(meta)
    meta <- cluster_tissue(meta)
    meta <- cluster_disease(meta)
    ct_idx <- grep("celltype", colnames(meta))
    if (length(ct_idx != 0)) {
        meta <- cluster_macro_celltype(meta)
    }
    return(meta)
}

library(future)
plan()
plan("multicore", workers = 4)
options(future.globals.maxSize = 51539607552)

routine_QC <- function(seurat) {
    library(Seurat)
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

ct_disease_acc_table <- function(aoi_seurat, exp_acc = 1) {
    dis_ct_table <- table(aoi_seurat$macro_celltype, aoi_seurat$DISEASE)
    if (exp_acc == 1) {
        acc_ct_table <- table(aoi_seurat$macro_celltype, aoi_seurat$gse_alias)
    } else {
        acc_ct_table <- table(aoi_seurat$macro_celltype, aoi_seurat$orig.ident)
    }
    ct_table <- cbind(acc_ct_table, dis_ct_table)
    # write.table(ct_table, "celltype_disease_accession_table.csv")
    return(ct_table)
}

convert_ENSMBL <- function(counts) {
  ensmbl_db <- read.csv("/home/zhang_jiaxuan/Lung_DB/utils/ENSMBL_gene_id_name.csv")
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

locate_seurat <- function(collection, seurat_folders) {
  f <- grep(collection, seurat_folders, value = TRUE)
  seurat <- readRDS(f)
  return(seurat)
}

subset_aoi <- function(aoi, seurat) {
  aoi_seurat <- subset(seurat, orig.ident %in% aoi)
  return(aoi_seurat)
}

down_sampling_accession <- function(seurat, threshold) {
    cluster_list <- table(seurat$macro_celltype)
    if (length(cluster_list[cluster_list < ncol(seurat) * 0.005]) > 0) {
        non_filt <- subset(seurat,
            macro_celltype %in%
            names(cluster_list[cluster_list < ncol(seurat) * 0.005]))
        cluster_list <- cluster_list[cluster_list >= ncol(seurat) * 0.005]
    }
    print(cluster_list)
    ps <- list()
    for (i in 1:length(cluster_list)) {
        tmp <- subset(seurat, macro_celltype == names(cluster_list[i]))
        # make sure the sum of number of cell will not exceed 
        tmp_size <- threshold * ncol(tmp)/ncol(seurat)
        cat("Current Celltype",
            names(cluster_list[i]),
            "limit is:", tmp_size, "\n")
        ps[[i]] <- tmp[, sample(colnames(tmp),
                    size = tmp_size, replace = FALSE)]
    }
    down_seurat <- merge(ps[[1]], ps[-1])
    if (exists("non_filt")) {
        down_seurat <- merge(down_seurat, non_filt)
    }
    down_seurat <- harmony_merge(down_seurat)
    return(down_seurat)
}

locate_seurat_per_group <- function(collection, 
                param, ncollect,
                group, seurat_folders, threshold = 30000) {
    f <- grep(collection, seurat_folders, value = TRUE)
    seurat <- readRDS(f)
    seurat <- cluster_macro_celltype(seurat)
    if (param == "disease") {
        seurat <- subset(seurat, DISEASE %in% group)
        print(ncol(seurat))
    }
    if (param == "tissue") {
        seurat <- subset(seurat, TISSUE %in% group)
    }
    if (param == "gender") {
        seurat <- subset(seurat, gender %in% group)
    }
    if (param == "age") {
        seurat <- subset(seurat, age %in% group)
    }
    if (param == "celltype") {
        seurat <- subset(seurat, celltype %in% group)
    }
    if (param == FALSE) {
        # cat(("Down sampling accessions of interest"), collection, "\n")
        #seurat <- down_sampling_accession(seurat)
    }
    threshold <- threshold/ncollect
    cat("Each collection cell number threshold is:", threshold, "\n")
    # return()
    #check cell number of selected seurat
    if (ncol(seurat) > threshold) {
        cat("Down sampling accessions of interest", collection, "\n")
        seurat <- down_sampling_accession(seurat,
                threshold = threshold)
    }
    return(seurat)
}

locate_seurat_per_acc <- function(acc, seurat_folders,
                 group_info, threshold = 30000) {
    f <- grep(acc, seurat_folders, value = TRUE)
    if (length(f) == 0) {
        return(NULL)
    }
    seurat <- readRDS(f)
    seurat <- cluster_macro_celltype(seurat)
    seurat$groups <- group_info
    #check cell number of selected seurat
    if (ncol(seurat) > threshold) {
        cat("Down sampling accessions of interest", acc, "\n")
        seurat <- down_sampling_accession(seurat,
                threshold = threshold)
    }
    return(seurat)
}
