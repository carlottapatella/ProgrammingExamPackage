# STEP 1: Base image
FROM rocker/r-ver:4.3.1

# STEP 2: System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libzstd-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    build-essential \
    pandoc \
    pandoc-citeproc \
    libglpk-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# STEP 3: Set working directory
WORKDIR /app

# STEP 4: Copy all project files
COPY . /app/

# STEP 5: Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# STEP 6: Install CRAN packages
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'Matrix', 'patchwork', 'readr', 'uwot', 'Rtsne', 'tibble', 'stringr', 'scales', 'viridisLite', 'rmarkdown', 'knitr'), repos='https://cloud.r-project.org')"

# STEP 7: Install Bioconductor packages
RUN R -e "BiocManager::install(c('Seurat', 'SingleR', 'celldex', 'SummarizedExperiment', 'rtracklayer', 'org.Hs.eg.db', 'AnnotationDbi', 'scRNAseq'))"


# STEP 8: Knit the RMarkdown file to HTML automatically when the container runs
CMD ["Rscript", "-e", "rmarkdown::render('exam_analysis_FINAL.Rmd')"]
