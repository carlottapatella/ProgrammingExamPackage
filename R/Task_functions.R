# ============================
# TASK 1a - Extract Gene Annotations from GTF
# ============================

#' Extract gene annotations from the GTF file you want to analyze
#'
#' Extracts gene_id and gene_type from a GTF file, keeping only entries where the feature type is "gene".
#'
#' @details "dplyr" package required
#' @param gtf_path Path to the GTF annotation file.
#' @return A data.frame with gene_id and gene_type of filtered elements annotated as genes
#' @examples
#' gene_info <- extract_gene_annotations_from_gtf("Homo_sapiens.GRCh38.111.gtf")
#' @export
extract_gene_annotations_from_gtf <- function(gtf_path) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install the 'dplyr' package.")
  gtf <- read.delim(gtf_path, header = FALSE, comment.char = "#")
  colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf_genes <- dplyr::filter(gtf, feature == "gene")

  # Helper function to extract specific attributes
  extract_attribute <- function(attr, key) {
    sub(paste0('.*', key, ' "([^"]+)".*'), "\1", attr)
  }

  gene_info <- data.frame(
    gene_id = extract_attribute(gtf_genes$attribute, "gene_id"),
    gene_type = extract_attribute(gtf_genes$attribute, "gene_type"),
    stringsAsFactors = FALSE
  )

  return(gene_info)
}

# ============================
# TASK 1b - Filter Seurat by Gene List
# ============================

#' Filter Seurat object using a list of gene IDs
#'
#' Filters the Seurat object to retain only genes specified in the provided list.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_ids Character vector of gene IDs to retain.
#' @return A filtered Seurat object containing only the selected genes.
#' @examples
#' seurat_filtered <- filter_seurat_by_genes(seurat_obj, gene_ids)
#' @export
filter_seurat_by_genes <- function(seurat_obj, gene_ids) {
  keep_genes <- rownames(seurat_obj)[rownames(seurat_obj) %in% gene_ids]
  seurat_obj_filtered <- subset(seurat_obj, features = keep_genes)
  return(seurat_obj_filtered)
}

# # ============================
# TASK 1c - Gene Annotation
# ============================

#' Filter Seurat object to keep only protein-coding genes
#'
#' This function reads a GTF file, extracts gene_id and gene_type using extract_attribute,
#' filters for protein-coding genes, and subsets the Seurat object accordingly.
#'
#' @details "dplyr" package required
#' @param seurat_obj A Seurat object.
#' @param gtf_path Path to the GTF annotation file.
#' @return A filtered Seurat object containing only protein-coding genes.
#' @examples
#' seurat_pc <- filter_protein_coding_genes(seurat_obj, "Homo_sapiens.GRCh38.111.gtf")
#' @export
filter_protein_coding_genes <- function(seurat_obj, gtf_path) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install the 'dplyr' package.")
  
  # Load and format GTF
  gtf <- read.delim(gtf_path, header = FALSE, comment.char = "#")
  colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  
  # Keep only gene-level annotations
  gtf_genes <- dplyr::filter(gtf, feature == "gene")
  
  # Helper function to extract specific attributes
  extract_attribute <- function(attr, key) {
    sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", attr)
  }
  
  # Extract gene_id and gene_type
  gene_ids <- extract_attribute(gtf_genes$attribute, "gene_id")
  gene_types <- extract_attribute(gtf_genes$attribute, "gene_type")
  
  # Keep only protein-coding genes
  protein_coding_ids <- gene_ids[gene_types == "protein_coding"]
  
  # Subset the Seurat object
  seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj) %in% protein_coding_ids)
  
  return(seurat_obj)
}


# ============================
# TASK 2 - Gene Expression Summary
# ============================

#' Plot number of genes with >=3 UMIs per cell
#'
#' @param seurat_obj A Seurat object.
#' @return A list with Seurat object with metadata column `genes_expressed_3UMI, and a violin plot.
#' @examples
#' result <- plot_gene_expression_summary(seurat_obj)
#' result$plot
#' @export
plot_gene_expression_summary <- function(seurat_obj) {
  umi_counts <- Matrix::colSums(seurat_obj@assays$RNA@counts >= 3)
  seurat_obj$genes_expressed_3UMI <- umi_counts
  plot <- VlnPlot(seurat_obj, features = "genes_expressed_3UMI") + NoLegend()
  return(list(seurat = seurat_obj, plot = plot))
}

# ============================
# TASK 3 - Remove Unwanted Gene Categories
# ============================
#'
#' Filters out genes based on name patterns: ribosomal, pseudogenes, and mitochondrial.
#'
#' @param seurat_obj A Seurat object.
#' @param ribo_pattern Regex for ribosomal proteins (default = "^RPS|^RPL").
#' @param pseudo_pattern Regex for pseudogenes (default = "pseudogene").
#' @param mito_prefix Prefix for mitochondrial genes (default = "MT-").
#' @return A list with filtered Seurat object and a summary table with counts of genes removed by category.
#' @examples
#' result <- filter_unwanted_gene_categories(seurat_obj)
#' result$summary_table
#' @export
filter_unwanted_gene_categories <- function(seurat_obj,
                                            ribo_pattern = "^RPS|^RPL",
                                            pseudo_pattern = "pseudogene",
                                            mito_prefix = "MT-") {
  gene_names <- rownames(seurat_obj)

  ribo_proteins <- grep(ribo_pattern, gene_names, value = TRUE)
  pseudo_genes <- grep(pseudo_pattern, gene_names, value = TRUE)
  mito_genes <- grep(mito_prefix, gene_names, value = TRUE)

  to_remove <- unique(c( ribo_proteins, pseudo_genes, mito_genes))
  filtered_seurat <- subset(seurat_obj, features = setdiff(gene_names, to_remove))

  summary <- data.frame(
    Category = c("Ribosomal", "Pseudogenes", "Mitochondrial"),
    Removed = c(length( ribo_proteins), length(pseudo_genes), length(mito_genes))
  )

  return(list(filtered_seurat = filtered_seurat, summary_table = summary))
}

# ============================
# TASK 4 - PCA
# ============================

#' Run PCA and plot variance explained
#'
#' Runs PCA on the Seurat object and plots the proportion of variance explained by each principal component.
#'
#' @param seurat_obj A Seurat object.
#' @param npcs Number of PCs to compute. Default 20.
#' @return A list with Seurat object and a variance barplot showing variance explained
#' @examples
#' result <- run_pca_and_plot(seurat_obj)
#' @export
run_pca_and_plot <- function(seurat_obj, npcs = 20) {
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs)
  stdev <- seurat_obj@reductions$pca@stdev[1:npcs]
  var_exp <- (stdev^2) / sum(stdev^2)
  plot <- barplot(var_exp, names.arg = 1:npcs, main = "Variance explained by PCs", xlab = "PC", ylab = "Proportion")
  return(list(seurat = seurat_obj, plot = plot))
}

# ============================
# TASK 5 - UMAP Visualization
# ============================

#' Run UMAP using selected principal componenet(PCs)
#'
#' Executes UMAP dimensionality reduction using user-defined principal components.
#'
#' @param seurat_obj A Seurat object.
#' @param dims PCs to use, e.g. 1:10.
#' @return A list with Seurat object and UMAP plot object with cluster labels.
#' @examples
#' result <- run_umap_custom(seurat_obj, dims = 1:10)
#' @export
run_umap_custom <- function(seurat_obj, dims) {
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  return(list(seurat = seurat_obj, plot = plot))
}

# ============================
# TASK 6 - Clustering
# ============================

#' Run clustering on Seurat object
#' 
#'Performs clustering using Seurat’s FindNeighbors and FindClusters functions, and returns a UMAP plot.
#'
#' @param seurat_obj A Seurat object.
#' @param dims Dimensions to use for clustering, e.g. 1:10.
#' @param resolution Clustering resolution. Default 0.5.
#' @return A list with Seurat object and clustering UMAP plot.
#' @examples
#' result <- run_clustering(seurat_obj, dims = 1:10)
#' @export
run_clustering <- function(seurat_obj, dims, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  return(list(seurat = seurat_obj, plot = plot))
}

# ============================
# TASK 7 - Cell Type Annotation
# ============================

#' Annotate cell types using SingleR
#'
#'Uses SingleR to assign cell type labels to each cell based on a reference.
#'
#' @param seurat_obj A Seurat object.
#' @param ref Reference matrix (genes × cell types).
#' @param labels Cell type labels corresponding to columns in "ref".
#' @return A list with Seurat object with cell type metadata and cell type UMAP plot labeled by cell type.
#' @examples
#' result <- annotate_cell_types(seurat_obj, ref_matrix, ref_labels)
#' @export
annotate_cell_types <- function(seurat_obj, ref, labels) {
  if (!requireNamespace("SingleR", quietly = TRUE)) stop("Install 'SingleR'.")
  pred <- SingleR::SingleR(test = GetAssayData(seurat_obj, slot = "data"), ref = ref, labels = labels)
  seurat_obj$cell_type <- pred$labels
  plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE)
  return(list(seurat = seurat_obj, plot = plot))
}

# ============================
# TASK 8 - Tissue Origin Inference
# ============================

#' Assign dominant cell type or tissue to each cluster
#'
#' Calculates the average expression of gene markers per cluster,
#' and assigns to each cluster the most likely tissue or cell type based on the highest marker signal.
#'
#' @param seurat_obj A Seurat object with clustering metadata.
#' @param marker_list A named list of marker genes (e.g. from derive_marker_list_from_reference()).
#' @return A data.frame with one row per cluster and the predicted cell type or tissue.
#' @examples
#' prediction <- assign_tissue_to_clusters(seurat_obj, marker_list)
#' @export
assign_tissue_to_clusters <- function(seurat_obj, marker_list) {
  avg_expr <- Seurat::AverageExpression(seurat_obj, return.seurat = FALSE)$RNA
  cluster_ids <- colnames(avg_expr)
  
  # Inizializza matrice marker score (tipo x cluster)
  score_matrix <- matrix(0, nrow = length(marker_list), ncol = length(cluster_ids))
  rownames(score_matrix) <- names(marker_list)
  colnames(score_matrix) <- cluster_ids
  
  for (type in names(marker_list)) {
    markers <- marker_list[[type]]
    markers <- intersect(markers, rownames(avg_expr))  # tieni solo marker presenti
    if (length(markers) > 0) {
      score_matrix[type, ] <- colMeans(avg_expr[markers, , drop = FALSE])
    }
  }
  
  # Per ogni cluster, assegna il tipo con punteggio massimo
  assignment <- apply(score_matrix, 2, function(x) names(which.max(x)))
  
  return(data.frame(
    cluster = cluster_ids,
    predicted_identity = assignment,
    stringsAsFactors = FALSE
  ))
}

