# Single-Cell RNA-Seq Analysis – Exam Project

This repository contains the full analysis pipeline for a single-cell RNA sequencing dataset (~10,000 cells, 10X Genomics 3’ technology) as part of the July exam project.

## Contents

- `exam_analysis.Rmd` – RMarkdown report with complete analysis and code  
- `functions.R` – Optional modular functions (if applicable)  
- `data/` – Folder containing input dataset and GTF file  
- `figures/` – All plots and visual outputs  

## Dataset

- **Technology**: 10X Genomics 3’ single-cell RNA-seq  
- **Cells**: ~10,000  
- **Input**:  
  - Raw gene expression matrix  
  - GTF file for gene annotation  

## Analysis Tasks

1. **Gene Annotation**  
   - Retained only protein-coding genes based on GTF.

2. **Gene Expression Summary**  
   - Counted genes per cell with ≥3 UMIs  
   - Visualized with violin plot

3. **Gene Filtering**  
   - Removed ribosomal, mitochondrial, and pseudogenes  
   - Provided summary table of filtered categories

4. **PCA**  
   - Computed first 20 PCs  
   - Visualized variance explained via histogram

5. **UMAP**  
   - Generated 2D embedding  
   - Selected components based on explained variance

6. **Clustering**  
   - Applied unsupervised clustering (method described in report)  
   - Discussed cluster resolution and biological interpretation

7. **Cell Type Annotation**  
   - Used automated annotation tool  
   - Compared results with clusters to evaluate biological consistency

8. **Tissue Origin Inference**  
   - Proposed tissue of origin  
   - Justified with marker genes and cluster identity

## 💻 How to Run Locally

1. Open `exam_analysis.Rmd` in RStudio  
2. Ensure required packages are installed (see report header)  
3. Click **"Knit"** to generate the HTML or PDF report  

## 🐳 Run with Docker

If you prefer a containerized environment, you can run the entire analysis via Docker.

docker run -v $(pwd):/app carlottapat/singlecell-project


