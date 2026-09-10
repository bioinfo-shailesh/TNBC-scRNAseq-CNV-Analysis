# TNBC-scRNAseq-CNV-Analysis
# Reproducible analysis workflow for the TNBC single-cell RNA-seq study.
#
# Workflow:
#   Raw 10x -> SoupX -> QC -> Seurat -> Harmony/integration (where applicable)
#   -> CopyKAT -> cell-type annotation -> refined tumor cells
#   -> CSC/FZD7/WNT analysis -> T/B/NK analyses -> SOCS3 -> statistics/figures
#
# FINAL scRNA-seq / CNV / CSC analysis master script
# See comments in script for switches and workflow.
# Code availability: The computational workflow is provided for reproducibility.
# Raw sequencing data and patient-level data are not included in this repository.

rm(list = ls()); gc(); set.seed(1234)

# ---------------------------------------------------------------------------
# GITHUB-READY CONFIGURATION
# ---------------------------------------------------------------------------
# Run this script from the root of the repository:
#   TNBC-scRNAseq-CNV-Analysis/
#
# Raw sequencing data and patient-level Seurat objects are NOT included in
# this public repository.
PROJECT_DIR <- "."
RAW_DATA_DIR <- file.path(PROJECT_DIR, "data", "raw")
NK_PROJECT_DIR <- file.path(PROJECT_DIR, "data", "local_nk")
RERUN_SOUPX <- FALSE
RERUN_QC <- FALSE
RERUN_COPYKAT <- FALSE
RERUN_ANNOTATION <- FALSE
RERUN_ANEUPLOID <- FALSE
FZD_THRESHOLD <- 1
IL10_THRESHOLD <- 0

pkgs <- c("Seurat","Matrix","dplyr","tidyr","ggplot2","patchwork","SoupX","copykat","harmony")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) stop("Missing package: ",p)
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(SoupX);library(copykat)})
setwd(PROJECT_DIR)

sample_info <- data.frame(
 Sample_ID=c("BSSR3692_005_GEX3","BSSR3890_007_GEX3","BSSR4069_006_GEX3","BSSR5829_T","BSSR6710_T","BSSR7011_T","BSSR7096_T","BSSR7548_T","BSSR7690_T","BSSR7762_T"),
 Response=c("NR","NR","R","R","NR","NR","NR","R","NR","NR"), stringsAsFactors=FALSE)

for(d in c("data/01_SoupX","data/02_QC_Filtered","data/03_CopyKAT","output/Final_Analysis/Figures","output/Final_Analysis/Tables","output/Final_Analysis/Statistics","output/Final_Analysis/Objects","logs")) dir.create(d,recursive=TRUE,showWarnings=FALSE)

save_plot <- function(p,f,w=10,h=8) ggsave(f,p,width=w,height=h,dpi=300,bg="white")
get_expr <- function(o,g){if(!g%in%rownames(o)) stop("Gene not found: ",g);as.numeric(FetchData(o,vars=g)[[g]])}

# ---------------------------------------------------------------------------
# SOUPX: raw 10X -> preliminary clustering -> ambient RNA correction
# ---------------------------------------------------------------------------

run_soupx <- function(sid, resp) {

  # Output file
  out <- file.path(
    PROJECT_DIR,
    "data",
    "01_SoupX",
    paste0(sid, "_SoupX.rds")
  )

  # Skip if already generated unless rerun is requested
  if (file.exists(out) && !RERUN_SOUPX) {
    message(sid, ": existing SoupX object found; skipping.")
    return(invisible(NULL))
  }

  # -------------------------------------------------------------------------
  # Load raw 10X data
  # -------------------------------------------------------------------------

  rawdir <- file.path(
    RAW_DATA_DIR,
    sid
  )

  if (!dir.exists(rawdir)) {
    stop(
      "Raw directory missing: ",
      rawdir
    )
  }

  raw <- Read10X(rawdir)

  # -------------------------------------------------------------------------
  # Create preliminary Seurat object
  # -------------------------------------------------------------------------

  s <- CreateSeuratObject(
    counts = raw,
    project = sid
  )

  s[["percent.mt"]] <- PercentageFeatureSet(
    s,
    pattern = "^MT-"
  )

  # -------------------------------------------------------------------------
  # Preliminary clustering for SoupX cluster assignment
  # -------------------------------------------------------------------------

  s <- NormalizeData(s)

  s <- FindVariableFeatures(s)

  s <- ScaleData(s)

  s <- RunPCA(s)

  s <- FindNeighbors(
    s,
    dims = 1:20
  )

  s <- FindClusters(
    s,
    resolution = 0.5
  )

  # -------------------------------------------------------------------------
  # Prepare SoupX matrices
  # -------------------------------------------------------------------------

  tod <- raw

  toc <- GetAssayData(
    s,
    layer = "counts"
  )

  common_genes <- intersect(
    rownames(tod),
    rownames(toc)
  )

  tod <- tod[
    common_genes,
    ,
    drop = FALSE
  ]

  toc <- toc[
    common_genes,
    ,
    drop = FALSE
  ]

  # Ensure identical gene order
  tod <- tod[
    order(rownames(tod)),
    ,
    drop = FALSE
  ]

  toc <- toc[
    order(rownames(toc)),
    ,
    drop = FALSE
  ]

  # -------------------------------------------------------------------------
  # Create SoupX channel and assign preliminary clusters
  # -------------------------------------------------------------------------

  sc <- SoupChannel(
    tod = tod,
    toc = toc
  )

  sc <- setClusters(
    sc,
    setNames(
      as.character(s$seurat_clusters),
      colnames(s)
    )
  )

  # -------------------------------------------------------------------------
  # Estimate ambient RNA contamination
  # -------------------------------------------------------------------------

  sc <- autoEstCont(sc)

  rho <- round(
    sc$fit$rhoEst * 100,
    2
  )

  # -------------------------------------------------------------------------
  # Correct ambient RNA contamination
  # -------------------------------------------------------------------------

  corrected_counts <- adjustCounts(sc)

  # -------------------------------------------------------------------------
  # Create corrected Seurat object
  # -------------------------------------------------------------------------

  ss <- CreateSeuratObject(
    counts = corrected_counts,
    project = sid
  )

  ss$Response <- resp
  ss$Sample_ID <- sid

  # -------------------------------------------------------------------------
  # Save SoupX objects
  # -------------------------------------------------------------------------

  saveRDS(
    sc,
    file.path(
      PROJECT_DIR,
      "data",
      "01_SoupX",
      paste0(sid, "_SoupX_channel.rds")
    )
  )

  saveRDS(
    ss,
    out
  )

  message(
    sid,
    ": SoupX rho = ",
    rho,
    "%"
  )

  invisible(ss)
}

# Run SoupX for all samples only when explicitly requested
if (RERUN_SOUPX) {
  for (i in seq_len(nrow(sample_info))) {

    run_soupx(
      sample_info$Sample_ID[i],
      sample_info$Response[i]
    )

  }
}

# ---------------------------------------------------------------------------
# FINAL QC: use existing validated QC objects unless explicitly rerun
# Initial/common QC: 200–7,000 genes and <=20% mitochondrial content.
# A second sample-specific QC refinement was subsequently applied.
# ---------------------------------------------------------------------------
INITIAL_MIN_GENES <- 200
INITIAL_MAX_GENES <- 7000
MAX_MT_PERCENT <- 20

# Sample-specific refinement thresholds used after initial QC
qc_min <- c(500,200,200,500,200,200,200,200,200,500)
qc_max <- c(7000,5000,4000,6000,5000,6000,6000,6000,5000,7500)

if (RERUN_QC) {
  
  for (i in seq_len(nrow(sample_info))) {
    
    sid <- sample_info$Sample_ID[i]
    
    o <- readRDS(
      file.path(
        PROJECT_DIR,
        "data/01_SoupX",
        paste0(sid, "_SoupX.rds")
      )
    )
    
    if (!"percent.mt" %in% colnames(o@meta.data)) {
      o[["percent.mt"]] <- PercentageFeatureSet(
        o,
        pattern = "^MT-"
      )
    }
    
    # Initial/common QC
    o <- subset(
      o,
      subset =
        nFeature_RNA >= INITIAL_MIN_GENES &
        nFeature_RNA <= INITIAL_MAX_GENES &
        percent.mt <= MAX_MT_PERCENT
    )
    
    # Sample-specific QC refinement
    o <- subset(
      o,
      subset =
        nFeature_RNA >= qc_min[i] &
        nFeature_RNA <= qc_max[i] &
        percent.mt <= MAX_MT_PERCENT
    )
    
    o$Sample_ID <- sid
    o$Response <- sample_info$Response[i]
    
    saveRDS(
      o,
      file.path(
        PROJECT_DIR,
        "data/02_QC_Filtered",
        paste0(sid, "_QC.rds")
      )
    )
  }
}
# ---------------------------------------------------------------------------
# COPYKAT: SoupX-corrected QC-filtered objects -> CNV predictions
# ---------------------------------------------------------------------------
# CopyKAT was run independently for each sample using the
# SoupX-corrected, QC-filtered expression matrix.
#
# Parameters used in the analysis:
#   id.type       = "S"
#   cell.line     = "no"
#   ngene.chr     = 5
#   LOW.DR        = 0.05
#   UP.DR         = 0.10
#   win.size      = 25
#   KS.cut        = 0.10
#   distance      = "euclidean"
#   output.seg    = FALSE
#   plot.genes    = TRUE
#   genome        = "hg20"
#   n.cores       = 1
# ---------------------------------------------------------------------------

run_ck <- function(sid) {

  q <- file.path(
    PROJECT_DIR,
    "data/02_QC_Filtered",
    paste0(sid, "_QC.rds")
  )

  od <- file.path(
    PROJECT_DIR,
    "data/03_CopyKAT",
    sid
  )

  dir.create(
    od,
    recursive = TRUE,
    showWarnings = FALSE
  )

  prediction_file <- file.path(
    od,
    "copykat_prediction.txt"
  )

  if (file.exists(prediction_file) && !RERUN_COPYKAT) {
    return()
  }

  o <- readRDS(q)

  counts <- as(
    GetAssayData(
      o,
      assay = "RNA",
      layer = "counts"
    ),
    "dgCMatrix"
  )

  ck <- copykat(
    rawmat = counts,
    id.type = "S",
    cell.line = "no",
    ngene.chr = 5,
    LOW.DR = 0.05,
    UP.DR = 0.10,
    win.size = 25,
    norm.cell.names = "",
    KS.cut = 0.10,
    sam.name = sid,
    distance = "euclidean",
    output.seg = FALSE,
    plot.genes = TRUE,
    genome = "hg20",
    n.cores = 1
  )

  saveRDS(
    ck,
    file.path(
      od,
      "copykat_result.rds"
    )
  )

  write.table(
    ck$prediction,
    prediction_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (!is.null(ck$CNAmat)) {

    write.table(
      ck$CNAmat,
      file.path(
        od,
        "copykat_CNA_matrix.txt"
      ),
      sep = "\t",
      quote = FALSE
    )
  }
}

if (RERUN_COPYKAT) {
  for (sid in sample_info$Sample_ID) {
    run_ck(sid)
  }
}

# ---------------------------------------------------------------------------
# LOAD QC-FILTERED OBJECTS AND SEURAT PREPROCESSING
# normalization -> variable features -> scaling -> PCA
# ---------------------------------------------------------------------------

seurat.list <- lapply(
  sample_info$Sample_ID,
  function(sid) {

    qc_file <- file.path(
      PROJECT_DIR,
      "data",
      "02_QC_Filtered",
      paste0(sid, "_QC.rds")
    )

    if (!file.exists(qc_file)) {
      stop("QC-filtered object not found: ", qc_file)
    }

    readRDS(qc_file)
  }
)

names(seurat.list) <- sample_info$Sample_ID

# ---------------------------------------------------------------------------
# Normalize, identify variable genes, scale, and run PCA
# ---------------------------------------------------------------------------

seurat.list <- lapply(
  seurat.list,
  function(x) {

    x <- NormalizeData(
      x,
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )

    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 2000
    )

    x <- ScaleData(x)

    x <- RunPCA(
      x,
      features = VariableFeatures(x),
      npcs = 50
    )

    x
  }
)
# ---------------------------------------------------------------------------
# LOAD FINAL VALIDATED ANNOTATED SEURAT OBJECT
# ---------------------------------------------------------------------------
# This object was generated during the study from the QC-filtered samples
# and CopyKAT results. It is not distributed in this public repository
# because it contains patient-level single-cell data.
#
# Expected local file:
#   data/combined_CopyKAT_CellType_annotated.rds
# ---------------------------------------------------------------------------

ANNOTATED_OBJECT <- file.path(
  PROJECT_DIR,
  "data",
  "combined_CopyKAT_CellType_annotated.rds"
)

if (!file.exists(ANNOTATED_OBJECT)) {
  stop(
    "Final annotated Seurat object not found.\n",
    "Expected local file: ",
    ANNOTATED_OBJECT
  )
}

combined <- readRDS(ANNOTATED_OBJECT)

DefaultAssay(combined) <- "RNA"

cat(
  "Loaded final annotated Seurat object: ",
  ncol(combined),
  " cells and ",
  nrow(combined),
  " genes.\n",
  sep = ""
)

# ---------------------------------------------------------------------------
# VALIDATED CELL-TYPE ANNOTATION
# ---------------------------------------------------------------------------
# Cluster annotations established during the study:
#
# Tumor/Epithelial      : clusters 2, 3, 7, 8, 10, 11
# T cells               : clusters 0, 1, 9, 13
# NK cells              : cluster 6
# Fibroblast            : cluster 5
# Myeloid               : clusters 4, 18
# B cells               : clusters 12, 16
# Pericyte/Endothelial  : cluster 15
# Undefined             : clusters 14, 17
#
# These annotations were based on canonical cell-type marker expression
# and subsequent manual review.

if (RERUN_ANNOTATION) {

  annotation_map <- c(
    "0"  = "T cells",
    "1"  = "T cells",
    "2"  = "Tumor/Epithelial",
    "3"  = "Tumor/Epithelial",
    "4"  = "Myeloid",
    "5"  = "Fibroblast",
    "6"  = "NK cells",
    "7"  = "Tumor/Epithelial",
    "8"  = "Tumor/Epithelial",
    "9"  = "T cells",
    "10" = "Tumor/Epithelial",
    "11" = "Tumor/Epithelial",
    "12" = "B cells",
    "13" = "T cells",
    "14" = "Undefined",
    "15" = "Pericyte/Endothelial",
    "16" = "B cells",
    "17" = "Undefined",
    "18" = "Myeloid"
  )

  combined$Cell_Population <- unname(
    annotation_map[as.character(combined$seurat_clusters)]
  )

  combined$Cell_Population[
    is.na(combined$Cell_Population)
  ] <- "Undefined"

  saveRDS(
    combined,
    ANNOTATED_OBJECT
  )
}

# Report the number of cells in each annotated population
population_counts <- as.data.frame(
  table(Population = combined$Cell_Population)
)

write.csv(
  population_counts,
  file.path(
    PROJECT_DIR,
    "output/Final_Analysis/Tables/Population_Cell_Counts.csv"
  ),
  row.names = FALSE
)
# ---------------------------------------------------------------------------
# REFINED ANEUPLOID TUMOR
# CopyKAT aneuploid cells followed by removal of immune and
# fibroblast-associated contamination.
#
# Exclusion criteria:
#   PTPRC > 0  -> excluded
#   SPARC >= 1 -> excluded
#
# Final refined tumor population: 5,378 cells
# ---------------------------------------------------------------------------

tumor <- subset(
  combined,
  subset = CNV_status == "aneuploid"
)

# Remove PTPRC-positive immune cells
if ("PTPRC" %in% rownames(tumor)) {
  tumor <- tumor[
    ,
    get_expr(tumor, "PTPRC") < 1
  ]
}

# Remove SPARC-high stromal/fibroblast contamination
if ("SPARC" %in% rownames(tumor)) {
  tumor <- tumor[
    ,
    get_expr(tumor, "SPARC") < 1
  ]
}

tumor$Tumor_Refined <- "Tumor"

saveRDS(
  tumor,
  file.path(
    PROJECT_DIR,
    "output/Final_Analysis/Objects/Final_Refined_Tumor.rds"
  )
)

cat(
  "Final refined tumor cells: ",
  ncol(tumor),
  "\n",
  sep = ""
)
# ---------------------------------------------------------------------------
# CSC / FZD7 ANNOTATION
# ---------------------------------------------------------------------------
# CSC definition:
#   CD44 > 0 AND CD24 <= 0
#
# FZD7-positive definition:
#   FZD7 >= 1
#
# These annotations were performed on the refined aneuploid tumor
# population after exclusion of PTPRC-positive and SPARC-high cells.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Cancer stem cell (CSC) annotation
# ---------------------------------------------------------------------------

if (all(c("CD44", "CD24") %in% rownames(tumor))) {

  CD44_expr <- get_expr(tumor, "CD44")
  CD24_expr <- get_expr(tumor, "CD24")

  tumor$CSC_All <- ifelse(
    CD44_expr > 0 & CD24_expr <= 0,
    "CSC",
    "Non-CSC"
  )

} else {

  warning("CD44 and/or CD24 not found in tumor object.")

  if (!"CSC_All" %in% colnames(tumor@meta.data)) {
    tumor$CSC_All <- NA_character_
  }
}


# ---------------------------------------------------------------------------
# FZD7 annotation
# ---------------------------------------------------------------------------

if ("FZD7" %in% rownames(tumor)) {

  FZD7_expr <- get_expr(tumor, "FZD7")

  tumor$FZD7_Status <- ifelse(
    FZD7_expr >= 1,
    "FZD7+",
    "FZD7-"
  )

} else {

  warning("FZD7 not found in tumor object.")
}


# ---------------------------------------------------------------------------
# Combined FZD7 / CSC groups
# ---------------------------------------------------------------------------

if (all(c("CSC_All", "FZD7_Status") %in% colnames(tumor@meta.data))) {

  tumor$FZD7_CSC_Group <- NA_character_

  tumor$FZD7_CSC_Group[
    tumor$FZD7_Status == "FZD7+" &
    tumor$CSC_All == "CSC"
  ] <- "FZD7+ CSC+"

  tumor$FZD7_CSC_Group[
    tumor$FZD7_Status == "FZD7+" &
    tumor$CSC_All == "Non-CSC"
  ] <- "FZD7+ CSC-"

  tumor$FZD7_CSC_Group[
    tumor$FZD7_Status == "FZD7-"
  ] <- "FZD7-"
}


# ---------------------------------------------------------------------------
# Cell counts
# ---------------------------------------------------------------------------

if (all(c("CSC_All", "FZD7_Status", "FZD7_CSC_Group") %
        in% colnames(tumor@meta.data))) {

  CSC_FZD7_counts <- tumor@meta.data %>%
    count(
      Response,
      CSC_All,
      FZD7_Status,
      FZD7_CSC_Group,
      name = "Cell_Count"
    )

  write.csv(
    CSC_FZD7_counts,
    file.path(
      PROJECT_DIR,
      "output/Final_Analysis/Tables",
      "CSC_FZD7_Cell_Counts.csv"
    ),
    row.names = FALSE
  )
}


# Save tumor object containing CSC/FZD7 annotations
saveRDS(
  tumor,
  file.path(
    PROJECT_DIR,
    "output/Final_Analysis/Objects",
    "Final_Refined_Tumor_CSC_FZD7.rds"
  )
)
# ---------------------------------------------------------------------------
# WNT TARGET GENE ANALYSIS
# ---------------------------------------------------------------------------
# WNT-associated genes analyzed in the refined tumor population:
# PROM1, PROM2, ALDH1A1, ALDH1A3, SOX9, CCND1, LEF1, TCF7, AXIN2
#
# Analyses were performed after establishment of the refined aneuploid
# tumor population and included:
#   1. CSC vs Non-CSC
#   2. NR vs R in the total tumor population
#   3. NR vs R within FZD7-positive tumor cells
# ---------------------------------------------------------------------------

wnt <- c(
  "PROM1",
  "PROM2",
  "ALDH1A1",
  "ALDH1A3",
  "SOX9",
  "CCND1",
  "LEF1",
  "TCF7",
  "AXIN2"
)

# Retain WNT genes present in the tumor expression matrix
wnt_present <- wnt[wnt %in% rownames(tumor)]

# Record gene availability
wnt_availability <- data.frame(
  Gene = wnt,
  Present = wnt %in% rownames(tumor)
)

write.csv(
  wnt_availability,
  file.path(
    PROJECT_DIR,
    "output/Final_Analysis/Tables",
    "WNT_Target_Gene_Availability.csv"
  ),
  row.names = FALSE
)


# ---------------------------------------------------------------------------
# Plotting function for WNT target genes
# ---------------------------------------------------------------------------

make_wnt_plots <- function(
    object,
    genes,
    group_variable,
    analysis_name
) {

  genes <- genes[genes %in% rownames(object)]

  if (length(genes) == 0) {
    warning(
      "No WNT target genes available for: ",
      analysis_name
    )
    return(invisible(NULL))
  }

  # Dot plot
  dot_plot <- DotPlot(
    object,
    features = genes,
    group.by = group_variable,
    dot.scale = 8
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    labs(
      color = "Average Expression",
      size = "Percent Expressed",
      x = "Gene",
      y = NULL,
      title = analysis_name
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "italic"
      ),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_plot(
    dot_plot,
    file.path(
      PROJECT_DIR,
      "output/Final_Analysis/Figures",
      paste0(analysis_name, "_DotPlot.png")
    ),
    12,
    8
  )


  # Violin plot
  violin_plot <- VlnPlot(
    object,
    features = genes,
    group.by = group_variable,
    pt.size = 0.05,
    ncol = 1
  ) &
    theme_classic(base_size = 12) &
    theme(
      axis.text.x = element_text(face = "bold"),
      legend.position = "right"
    )

  save_plot(
    violin_plot,
    file.path(
      PROJECT_DIR,
      "output/Final_Analysis/Figures",
      paste0(analysis_name, "_ViolinPlot.png")
    ),
    12,
    14
  )
}


# ---------------------------------------------------------------------------
# 1. WNT targets: CSC vs Non-CSC
# ---------------------------------------------------------------------------

if ("CSC_All" %in% colnames(tumor@meta.data)) {

  make_wnt_plots(
    tumor,
    wnt_present,
    "CSC_All",
    "WNT_Targets_CSC_vs_NonCSC"
  )
}


# ---------------------------------------------------------------------------
# 2. WNT targets: NR vs R in total refined tumor population
# ---------------------------------------------------------------------------

tumor$Response <- factor(
  tumor$Response,
  levels = c("NR", "R")
)

make_wnt_plots(
  tumor,
  wnt_present,
  "Response",
  "WNT_Targets_Total_Tumor_NR_vs_R"
)


# ---------------------------------------------------------------------------
# 3. WNT targets: NR vs R within FZD7-positive tumor cells
# ---------------------------------------------------------------------------

if ("FZD7_Status" %in% colnames(tumor@meta.data)) {

  FZD7_positive_tumor <- subset(
    tumor,
    subset = FZD7_Status == "FZD7+"
  )

  if (ncol(FZD7_positive_tumor) > 0) {

    make_wnt_plots(
      FZD7_positive_tumor,
      wnt_present,
      "Response",
      "WNT_Targets_FZD7_Positive_NR_vs_R"
    )
  }
}
# ---------------------------------------------------------------------------
# SOCS3 ANALYSIS
# ---------------------------------------------------------------------------
# SOCS3 expression was evaluated in:
#
#   1. Refined tumor population: NR vs R
#   2. Major cell populations:
#        Tumor/Epithelial
#        Myeloid
#        NK cells
#        B cells
#        Fibroblast
#        T cells
#   3. Patient-wise tumor analysis
#
# Expression values are normalized expression values from the RNA assay.
# ---------------------------------------------------------------------------


if ("SOCS3" %in% rownames(combined)) {


  # -------------------------------------------------------------------------
  # 1. SOCS3 IN REFINED TUMOR: NR vs R
  # -------------------------------------------------------------------------

  if ("Response" %in% colnames(tumor@meta.data)) {

    tumor$Response <- factor(
      tumor$Response,
      levels = c("NR", "R")
    )

    # Violin plot
    p <- VlnPlot(
      tumor,
      features = "SOCS3",
      group.by = "Response",
      pt.size = 0.05
    ) +
      theme_classic(base_size = 14) +
      labs(
        x = NULL,
        y = "Normalized Expression",
        title = "SOCS3 in Tumor: NR vs R"
      )

    save_plot(
      p,
      file.path(
        PROJECT_DIR,
        "output/Final_Analysis/Figures",
        "SOCS3_Tumor_NR_vs_R_Violin.png"
      ),
      8,
      7
    )


    # Dot plot
    d <- DotPlot(
      tumor,
      features = "SOCS3",
      group.by = "Response",
      dot.scale = 8
    ) +
      scale_color_gradient(
        low = "#2166AC",
        high = "#B2182B"
      ) +
      theme_classic(base_size = 14) +
      labs(
        x = NULL,
        y = NULL,
        color = "Average Expression",
        size = "Percent Expressed",
        title = "SOCS3 in Tumor: NR vs R"
      )

    save_plot(
      d,
      file.path(
        PROJECT_DIR,
        "output/Final_Analysis/Figures",
        "SOCS3_Tumor_NR_vs_R_DotPlot.png"
      ),
      7,
      6
    )
  }


  # -------------------------------------------------------------------------
  # 2. SOCS3 ACROSS MAJOR CELL POPULATIONS
  # -------------------------------------------------------------------------
  
  populations_for_socs3 <- c(
    "Tumor/Epithelial",
    "Myeloid",
    "NK cells",
    "B cells",
    "Fibroblast",
    "T cells"
  )

  pops <- subset(
    combined,
    subset = Cell_Population %in% populations_for_socs3
  )

  # Set population order for plotting
  pops$Cell_Population <- factor(
    pops$Cell_Population,
    levels = populations_for_socs3
  )


  # Violin plot
  p <- VlnPlot(
    pops,
    features = "SOCS3",
    group.by = "Cell_Population",
    pt.size = 0.03
  ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    ) +
    labs(
      x = NULL,
      y = "Normalized Expression",
      title = "SOCS3 Across Cell Populations"
    )

  save_plot(
    p,
    file.path(
      PROJECT_DIR,
      "output/Final_Analysis/Figures",
      "SOCS3_All_Cell_Populations_Violin.png"
    ),
    12,
    8
  )


  # Dot plot
  d <- DotPlot(
    pops,
    features = "SOCS3",
    group.by = "Cell_Population",
    dot.scale = 8
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    theme_classic(base_size = 13) +
    labs(
      x = NULL,
      y = NULL,
      color = "Average Expression",
      size = "Percent Expressed",
      title = "SOCS3 Across Cell Populations"
    )

  save_plot(
    d,
    file.path(
      PROJECT_DIR,
      "output/Final_Analysis/Figures",
      "SOCS3_All_Cell_Populations_DotPlot.png"
    ),
    10,
    7
  )


  # -------------------------------------------------------------------------
  # 3. PATIENT-WISE SOCS3 IN REFINED TUMOR
  # -------------------------------------------------------------------------

  if ("Sample_ID" %in% colnames(tumor@meta.data)) {

    patient_socs3 <- FetchData(
      tumor,
      vars = c(
        "SOCS3",
        "Sample_ID",
        "Response"
      )
    ) %>%
      group_by(
        Sample_ID,
        Response
      ) %>%
      summarise(
        N_Cells = n(),
        Mean_SOCS3 = mean(SOCS3),
        Median_SOCS3 = median(SOCS3),
        Percent_GE1 = mean(SOCS3 >= 1) * 100,
        Percent_GT0 = mean(SOCS3 > 0) * 100,
        .groups = "drop"
      )


    # Save patient-wise SOCS3 table
    write.csv(
      patient_socs3,
      file.path(
        PROJECT_DIR,
        "output/Final_Analysis/Tables",
        "SOCS3_Tumor_Patient_Wise.csv"
      ),
      row.names = FALSE
    )


    # Patient-wise mean SOCS3 bar plot
    p <- ggplot(
      patient_socs3,
      aes(
        x = Sample_ID,
        y = Mean_SOCS3,
        fill = Response
      )
    ) +
      geom_col() +
      theme_classic(base_size = 13) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      ) +
      labs(
        x = "Patient",
        y = "Mean SOCS3 Expression",
        title = "Patient-wise SOCS3 Expression in Tumor"
      )

    save_plot(
      p,
      file.path(
        PROJECT_DIR,
        "output/Final_Analysis/Figures",
        "SOCS3_Tumor_Patient_Wise_BarPlot.png"
      ),
      12,
      7
    )
  }
}

# ---------------------------------------------------------------------------
# NK-CELL ANALYSIS
# ---------------------------------------------------------------------------
# Finalized project workflow:
#   integrated dataset -> NK identification using NCAM1/KLRC1/NKG7
#   -> removal of T-cell contamination -> 1,572-cell NK subset
#   -> NK-specific reclustering -> functional genes/SOCS3 analysis.
#
# NK identification markers:
#   NCAM1 (CD56), KLRC1 (NKG2A), NKG7
#
# T-cell contamination exclusion based on normalized expression:
#   CD3D > 1  OR CD3G > 1 OR CD3E > 0 OR CD4 > 1
#
# Retained NK cells therefore satisfy:
#   CD3D <= 1, CD3G <= 1, CD3E <= 0, CD4 <= 1
#
# Final NK population used in the project: 1,572 cells.
# Reclustering: first 25 PCs, resolution = 0.15.
# ---------------------------------------------------------------------------

nk_dir <- file.path(PROJECT_DIR,"output/Final_Analysis/NK_Analysis")
nk_fig <- file.path(nk_dir,"Figures")
nk_tab <- file.path(nk_dir,"Tables")
nk_obj <- file.path(nk_dir,"Objects")
for(d in c(nk_dir,nk_fig,nk_tab,nk_obj))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)

# Start from the integrated annotated object.
# If an already validated starting NK object is available locally, it can be
# supplied here, but the filtering below is always applied before reclustering.
if("Cell_Population_Final" %in% colnames(combined@meta.data)) {
  NK <- subset(combined, subset=Cell_Population_Final=="NK cells")
} else {
  NK <- subset(combined, subset=Cell_Population=="NK cells")
}
DefaultAssay(NK) <- "RNA"

# ---------------------------------------------------------------------------
# NK marker confirmation
# ---------------------------------------------------------------------------
nk_markers <- c("NCAM1","KLRC1","NKG7")
nk_markers_present <- nk_markers[nk_markers %in% rownames(NK)]

cat("\nInitial NK population: ",ncol(NK),"\n",sep="")

# ---------------------------------------------------------------------------
# Remove T-cell contamination using the exact normalized-expression criteria
# used in the project.
# ---------------------------------------------------------------------------
if("CD3D" %in% rownames(NK))
  NK <- NK[,get_expr(NK,"CD3D") <= 1]
if("CD3G" %in% rownames(NK))
  NK <- NK[,get_expr(NK,"CD3G") <= 1]
if("CD3E" %in% rownames(NK))
  NK <- NK[,get_expr(NK,"CD3E") <= 0]
if("CD4" %in% rownames(NK))
  NK <- NK[,get_expr(NK,"CD4") <= 1]

cat("Refined NK population after T-cell exclusion: ",ncol(NK),"\n",sep="")
print(table(NK$Response))

# Save the refined pre-reclustering NK population.
saveRDS(NK,file.path(nk_obj,"NK_Refined_1572.rds"))

# ---------------------------------------------------------------------------
# NK-specific reclustering
# ---------------------------------------------------------------------------
NK <- NormalizeData(NK,verbose=FALSE)
NK <- FindVariableFeatures(NK,selection.method="vst",
                           nfeatures=2000,verbose=FALSE)
NK <- ScaleData(NK,verbose=FALSE)
NK <- RunPCA(NK,npcs=25,verbose=FALSE)
NK <- FindNeighbors(NK,dims=1:25,verbose=FALSE)
NK <- FindClusters(NK,resolution=0.15,verbose=FALSE)
NK <- RunUMAP(NK,dims=1:25,reduction="pca",verbose=FALSE)

saveRDS(NK,file.path(nk_obj,"NK_Final.rds"))

# NK subcluster UMAP
p <- DimPlot(NK,reduction="umap",group.by="seurat_clusters",
             label=TRUE,repel=TRUE)+
  theme_classic(base_size=15)+
  labs(title="NK-cell Subclusters")
save_plot(p,file.path(nk_fig,"NK_Subclusters_UMAP.png"),9,7)

# NK marker expression
if(length(nk_markers_present)>0){
  p <- DotPlot(NK,features=nk_markers_present,
               group.by="seurat_clusters",dot.scale=8)+
    scale_color_gradientn(
      colours=c("#0000FF","#8A00CC","#FF0000"),
      name="Average Expression")+
    theme_classic(base_size=13)+
    theme(axis.text.x=element_text(angle=45,hjust=1,face="italic"))+
    labs(x=NULL,y="NK Subcluster",
         size="Percent Expressed",
         title="NK Marker Expression")
  save_plot(p,file.path(nk_fig,"NK_Marker_DotPlot.png"),11,7)
}

# ---------------------------------------------------------------------------
# NK functional genes
# CD107A = LAMP1; CD117 = KIT; VEGF = VEGFA
# ---------------------------------------------------------------------------
nk_function_genes <- c("LAMP1","KIT","VEGFA")
nk_function_present <- nk_function_genes[
  nk_function_genes %in% rownames(NK)]

if(length(nk_function_present)>0){

  p <- DotPlot(NK,features=nk_function_present,
               group.by="Response",dot.scale=9)+
    scale_color_gradientn(
      colours=c("#0000FF","#8A00CC","#FF0000"),
      name="Average Expression")+
    theme_classic(base_size=14)+
    labs(x=NULL,y=NULL,size="Percent Expressed",
         title="NK Functional Genes: NR vs R")
  save_plot(p,file.path(nk_fig,
                        "NK_Functional_Genes_NR_vs_R_DotPlot.png"),9,6)

  p <- VlnPlot(NK,features=nk_function_present,
               group.by="Response",pt.size=0.05,ncol=1)&
    theme_classic(base_size=13)&
    theme(axis.text.x=element_text(face="bold"))
  save_plot(p,file.path(nk_fig,
                        "NK_Functional_Genes_NR_vs_R_ViolinPlot.png"),8,12)

  nk_stat <- bind_rows(lapply(nk_function_present,function(g){
    z <- FetchData(NK,vars=c(g,"Response"))
    z <- z[z$Response %in% c("NR","R"),,drop=FALSE]

    if(length(unique(z$Response))<2)
      return(data.frame(
        Gene=g,
        NR_n=sum(z$Response=="NR"),
        R_n=sum(z$Response=="R"),
        P_value=NA_real_))

    w <- wilcox.test(z[[g]]~z$Response,exact=FALSE)

    data.frame(
      Gene=g,
      NR_n=sum(z$Response=="NR"),
      R_n=sum(z$Response=="R"),
      NR_median=median(z[[g]][z$Response=="NR"]),
      R_median=median(z[[g]][z$Response=="R"]),
      P_value=w$p.value)
  }))

  nk_stat$FDR_BH <- p.adjust(nk_stat$P_value,"BH")

  write.csv(
    nk_stat,
    file.path(nk_tab,
             "NK_Functional_Genes_NR_vs_R_Wilcoxon_BH.csv"),
    row.names=FALSE)
}

# ---------------------------------------------------------------------------
# SOCS3 in NK cells
# ---------------------------------------------------------------------------
if("SOCS3" %in% rownames(NK)){

  p <- VlnPlot(NK,features="SOCS3",
               group.by="Response",pt.size=0.05)+
    theme_classic(base_size=14)+
    labs(x=NULL,y="SOCS3 Expression",
         title="SOCS3 in NK Cells: NR vs R")
  save_plot(p,file.path(nk_fig,
                        "NK_SOCS3_NR_vs_R_ViolinPlot.png"),8,7)

  p <- DotPlot(NK,features="SOCS3",
               group.by="Response",dot.scale=9)+
    scale_color_gradientn(
      colours=c("#0000FF","#8A00CC","#FF0000"),
      name="Average Expression")+
    theme_classic(base_size=14)+
    labs(x=NULL,y=NULL,size="Percent Expressed",
         title="SOCS3 in NK Cells: NR vs R")
  save_plot(p,file.path(nk_fig,
                        "NK_SOCS3_NR_vs_R_DotPlot.png"),7,6)

  z <- FetchData(NK,vars=c("SOCS3","Response"))
  z <- z[z$Response %in% c("NR","R"),,drop=FALSE]

  if(length(unique(z$Response))==2){
    w <- wilcox.test(z$SOCS3~z$Response,exact=FALSE)

    nk_s <- data.frame(
      Gene="SOCS3",
      NR_n=sum(z$Response=="NR"),
      R_n=sum(z$Response=="R"),
      NR_median=median(z$SOCS3[z$Response=="NR"]),
      R_median=median(z$SOCS3[z$Response=="R"]),
      P_value=w$p.value,
      FDR_BH=p.adjust(w$p.value,"BH"))

    write.csv(
      nk_s,
      file.path(nk_tab,"NK_SOCS3_NR_vs_R_Wilcoxon_BH.csv"),
      row.names=FALSE)
  }
}

# Unsupervised NK subcluster markers
nk_deg <- FindAllMarkers(
  NK,
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25)

write.csv(
  nk_deg,
  file.path(nk_tab,"NK_Subcluster_Markers_FindAllMarkers.csv"),
  row.names=FALSE)

saveRDS(NK,file.path(nk_obj,"NK_Final.rds"))


# ---------------------------------------------------------------------------
# Response-wise Wilcoxon + BH for WNT targets and SOCS3
# ---------------------------------------------------------------------------
wilcox_resp<-function(o,genes,pop){genes<-genes[genes%in%rownames(o)];bind_rows(lapply(genes,function(g){x<-FetchData(o,vars=c(g,"Response"));x<-x[x$Response%in%c("NR","R"),,drop=FALSE];if(length(unique(x$Response))<2)return(data.frame(Population=pop,Gene=g,NR_n=sum(x$Response=="NR"),R_n=sum(x$Response=="R"),P_value=NA_real_));w<-wilcox.test(x[[g]]~x$Response,exact=FALSE);data.frame(Population=pop,Gene=g,NR_n=sum(x$Response=="NR"),R_n=sum(x$Response=="R"),NR_median=median(x[[g]][x$Response=="NR"]),R_median=median(x[[g]][x$Response=="R"]),P_value=w$p.value)}))%>%mutate(FDR_BH=p.adjust(P_value,"BH"))}
write.csv(wilcox_resp(tumor,wnt_present,"Total Tumor"),file.path(PROJECT_DIR,"output/Final_Analysis/Statistics/WNT_Targets_Tumor_NR_vs_R_Wilcoxon_BH.csv"),row.names=FALSE);if("SOCS3"%in%rownames(tumor))write.csv(wilcox_resp(tumor,"SOCS3","Total Tumor"),file.path(PROJECT_DIR,"output/Final_Analysis/Statistics/SOCS3_Tumor_NR_vs_R_Wilcoxon_BH.csv"),row.names=FALSE)

# ---------------------------------------------------------------------------
# FINAL ANALYSIS OUTPUTS
# ---------------------------------------------------------------------------

# Summary cell-count report
final_cell_counts <- data.frame(
  Metric = c(
    "Total combined cells",
    "Aneuploid cells",
    "Refined tumor cells"
  ),
  Value = c(
    ncol(combined),
    sum(
      combined$CNV_status == "aneuploid",
      na.rm = TRUE
    ),
    ncol(tumor)
  )
)

write.csv(
  final_cell_counts,
  file.path(
    PROJECT_DIR,
    "output",
    "Final_Analysis",
    "Tables",
    "Final_Cell_Count_Report.csv"
  ),
  row.names = FALSE
)

# Save final annotated object locally
# This file contains patient-level single-cell data and should NOT be
# committed to the public GitHub repository.
saveRDS(
  combined,
  file.path(
    PROJECT_DIR,
    "output",
    "Final_Analysis",
    "Objects",
    "Final_Combined_Annotated_Seurat.rds"
  )
)

# Record software environment for reproducibility
sink(
  file.path(
    PROJECT_DIR,
    "output",
    "Final_Analysis",
    "SessionInfo.txt"
  )
)

cat("Final analysis session\n\n")
print(sessionInfo())

sink()

cat(
  "\nFINAL PIPELINE COMPLETED\n",
  "Output directory: ",
  file.path(
    PROJECT_DIR,
    "output",
    "Final_Analysis"
  ),
  "\n",
  sep = ""
)

