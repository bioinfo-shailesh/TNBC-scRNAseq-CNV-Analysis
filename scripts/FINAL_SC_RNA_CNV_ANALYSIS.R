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

pkgs <- c("Seurat","Matrix","dplyr","tidyr","ggplot2","patchwork","SoupX","copykat")
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
run_soupx <- function(sid,resp){
 out <- file.path(PROJECT_DIR,"data/01_SoupX",paste0(sid,"_SoupX.rds"))
 if(file.exists(out)&&!RERUN_SOUPX)return()
 rawdir<-file.path(RAW_DATA_DIR,sid); if(!dir.exists(rawdir))stop("Raw directory missing: ",rawdir)
 raw<-Read10X(rawdir); s<-CreateSeuratObject(raw,project=sid)
 s[["percent.mt"]]<-PercentageFeatureSet(s,pattern="^MT-")
 s<-NormalizeData(s);s<-FindVariableFeatures(s);s<-ScaleData(s);s<-RunPCA(s);s<-FindNeighbors(s,dims=1:20);s<-FindClusters(s,resolution=.5)
 tod<-raw;toc<-GetAssayData(s,layer="counts");g<-intersect(rownames(tod),rownames(toc));tod<-tod[g,,drop=FALSE];toc<-toc[g,,drop=FALSE];tod<-tod[order(rownames(tod)),,drop=FALSE];toc<-toc[order(rownames(toc)),,drop=FALSE]
 sc<-SoupChannel(tod,toc);sc<-setClusters(sc,setNames(as.character(s$seurat_clusters),colnames(s)));sc<-autoEstCont(sc);rho<-round(sc$fit$rhoEst*100,2);corr<-adjustCounts(sc)
 ss<-CreateSeuratObject(corr,project=sid);ss$Response<-resp;ss$Sample_ID<-sid
 saveRDS(sc,file.path(PROJECT_DIR,"data/01_SoupX",paste0(sid,"_SoupX_channel.rds")));saveRDS(ss,out);message(sid,": SoupX rho = ",rho,"%")
}
if(RERUN_SOUPX) for(i in seq_len(nrow(sample_info))) run_soupx(sample_info$Sample_ID[i],sample_info$Response[i])

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
# LOAD FINAL VALIDATED ANNOTATED SEURAT OBJECT
# ---------------------------------------------------------------------------
# The final annotated Seurat object is a local analysis input and is not
# distributed in this public repository.
#
# Before running this section, place the locally generated annotated object at:
#   data/combined_CopyKAT_CellType_annotated.rds
#
# The object should contain the CopyKAT CNV classification and the integrated
# Seurat metadata required for downstream analyses.
# ---------------------------------------------------------------------------

ANNOTATED_OBJECT <- file.path(
  PROJECT_DIR,
  "data",
  "combined_CopyKAT_CellType_annotated.rds"
)

if (!file.exists(ANNOTATED_OBJECT)) {
  stop(
    "Final annotated Seurat object not found. ",
    "Please provide the local file: ",
    ANNOTATED_OBJECT
  )
}

combined <- readRDS(ANNOTATED_OBJECT)

DefaultAssay(combined) <- "RNA"

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
# REFINED ANEUPLOID TUMOR: CopyKAT aneuploid + PTPRC <=0 + SPARC <1
# ---------------------------------------------------------------------------
tumor<-subset(combined,subset=CNV_status=="aneuploid")
if("PTPRC"%in%rownames(tumor))tumor<-tumor[,get_expr(tumor,"PTPRC")<=0]
if("SPARC"%in%rownames(tumor))tumor<-tumor[,get_expr(tumor,"SPARC")<1]
tumor$Tumor_Refined<-"Tumor";saveRDS(tumor,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/Final_Refined_Tumor.rds"))

# ---------------------------------------------------------------------------
# CSC / FZD7 annotations
# CSC definition in project: CD44+ CD24-. FZD7+ threshold used in plots: >=1.
# ---------------------------------------------------------------------------
if(all(c("CD44","CD24")%in%rownames(tumor))){a<-get_expr(tumor,"CD44");b<-get_expr(tumor,"CD24");tumor$CSC_All<-ifelse(a>0&b<=0,"CSC","Non-CSC")}else if(!"CSC_All"%in%colnames(tumor@meta.data))tumor$CSC_All<-NA_character_
if("FZD7"%in%rownames(tumor)){f<-get_expr(tumor,"FZD7");tumor$FZD7_Status<-ifelse(f>=1,"FZD7+","FZD7-");tumor$FZD7_CSC_Group<-ifelse(f>=1&tumor$CSC_All=="CSC","FZD7+ CSC+",ifelse(f>=1&tumor$CSC_All!="CSC","FZD7+ CSC-","FZD7-"))}
if(all(c("CSC_All","FZD7_Status")%in%colnames(tumor@meta.data)))write.csv(tumor@meta.data%>%count(Response,CSC_All,FZD7_Status,FZD7_CSC_Group,name="Cell_Count"),file.path(PROJECT_DIR,"output/Final_Analysis/Tables/CSC_FZD7_Cell_Counts.csv"),row.names=FALSE)

# ---------------------------------------------------------------------------
# WNT targets used in the project
# ---------------------------------------------------------------------------
wnt<-c("PROM1","PROM2","ALDH1A1","ALDH1A3","SOX9","CCND1","LEF1","TCF7","AXIN2");wp<-wnt[wnt%in%rownames(tumor)];write.csv(data.frame(Gene=wnt,Present=wnt%in%rownames(tumor)),file.path(PROJECT_DIR,"output/Final_Analysis/Tables/WNT_Target_Gene_Availability.csv"),row.names=FALSE)

make_plots<-function(o,genes,grp,name){genes<-genes[genes%in%rownames(o)];if(!length(genes))return();dp<-DotPlot(o,features=genes,group.by=grp,dot.scale=8)+scale_color_gradient(low="#2166AC",high="#B2182B")+labs(color="Average Expression",size="Percent Expressed",x="Gene",y=NULL,title=name)+theme_classic(base_size=13)+theme(axis.text.x=element_text(angle=45,hjust=1,face="italic"),axis.text.y=element_text(face="bold"),plot.title=element_text(face="bold",hjust=.5));vp<-VlnPlot(o,features=genes,group.by=grp,pt.size=.05,ncol=1)&theme_classic(base_size=12)&theme(axis.text.x=element_text(face="bold"),legend.position="right");save_plot(dp,file.path(PROJECT_DIR,"output/Final_Analysis/Figures",paste0(name,"_DotPlot.png")),12,8);save_plot(vp,file.path(PROJECT_DIR,"output/Final_Analysis/Figures",paste0(name,"_ViolinPlot.png")),12,14)}
if("CSC_All"%in%colnames(tumor@meta.data))make_plots(tumor,wp,"CSC_All","WNT_Targets_CSC_vs_NonCSC")
tumor$Response<-factor(tumor$Response,levels=c("NR","R"));make_plots(tumor,wp,"Response","WNT_Targets_Total_Tumor_NR_vs_R")
if("FZD7_Status"%in%colnames(tumor@meta.data)){fp<-subset(tumor,subset=FZD7_Status=="FZD7+");if(ncol(fp)>0)make_plots(fp,wp,"Response","WNT_Targets_FZD7_Positive_NR_vs_R")}

# ---------------------------------------------------------------------------
# SOCS3: tumor, populations and patient-wise
# ---------------------------------------------------------------------------
if("SOCS3"%in%rownames(combined)){
 if("Response"%in%colnames(tumor@meta.data)){p<-VlnPlot(tumor,features="SOCS3",group.by="Response",pt.size=.05)+theme_classic(base_size=14)+labs(x=NULL,y="Normalized Expression",title="SOCS3 in Tumor: NR vs R");save_plot(p,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_Tumor_NR_vs_R_Violin.png"),8,7);d<-DotPlot(tumor,features="SOCS3",group.by="Response",dot.scale=8)+scale_color_gradient(low="#2166AC",high="#B2182B")+theme_classic(base_size=14)+labs(x=NULL,y=NULL,color="Average Expression",size="Percent Expressed",title="SOCS3 in Tumor: NR vs R");save_plot(d,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_Tumor_NR_vs_R_DotPlot.png"),7,6)}
 pops<-subset(combined,subset=Cell_Population%in%c("Tumor/Epithelial","Myeloid","NK cells","B cells","Fibroblast","T cells"));p<-VlnPlot(pops,features="SOCS3",group.by="Cell_Population",pt.size=.03)+theme_classic(base_size=13)+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y="Normalized Expression",title="SOCS3 Across Cell Populations");save_plot(p,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_All_Cell_Populations_Violin.png"),12,8);d<-DotPlot(pops,features="SOCS3",group.by="Cell_Population",dot.scale=8)+scale_color_gradient(low="#2166AC",high="#B2182B")+theme_classic(base_size=13)+labs(x=NULL,y=NULL,color="Average Expression",size="Percent Expressed",title="SOCS3 Across Cell Populations");save_plot(d,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_All_Cell_Populations_DotPlot.png"),10,7)
 if("Sample_ID"%in%colnames(tumor@meta.data)){ps<-FetchData(tumor,vars=c("SOCS3","Sample_ID","Response"))%>%group_by(Sample_ID,Response)%>%summarise(N_Cells=n(),Mean_SOCS3=mean(SOCS3),Median_SOCS3=median(SOCS3),Percent_GE1=mean(SOCS3>=1)*100,Percent_GT0=mean(SOCS3>0)*100,.groups="drop");write.csv(ps,file.path(PROJECT_DIR,"output/Final_Analysis/Tables/SOCS3_Tumor_Patient_Wise.csv"),row.names=FALSE);p<-ggplot(ps,aes(Sample_ID,Mean_SOCS3,fill=Response))+geom_col()+theme_classic(base_size=13)+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x="Patient",y="Mean SOCS3 Expression",title="Patient-wise SOCS3 Expression in Tumor");save_plot(p,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_Tumor_Patient_Wise_BarPlot.png"),12,7)}
}

# ---------------------------------------------------------------------------
# T and B cells
# ---------------------------------------------------------------------------
tc<-subset(combined,subset=Cell_Population=="T cells");bc<-subset(combined,subset=Cell_Population=="B cells");dir.create(file.path(PROJECT_DIR,"output/Final_Analysis/Objects/T_Cells"),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(PROJECT_DIR,"output/Final_Analysis/Objects/B_Cells"),recursive=TRUE,showWarnings=FALSE);saveRDS(tc,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/T_Cells/T_cells.rds"));saveRDS(bc,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/B_Cells/B_cells.rds"))
if(all(c("CD4","CD8A","FOXP3")%in%rownames(tc))){a<-get_expr(tc,"CD4");b<-get_expr(tc,"CD8A");c<-get_expr(tc,"FOXP3");tc$T_Cell_Subtype<-ifelse(c>0,"FOXP3+",ifelse(b>0,"CD8",ifelse(a>0,"CD4","Other")));if("SOCS3"%in%rownames(tc)){p<-VlnPlot(tc,features="SOCS3",group.by="T_Cell_Subtype",pt.size=.05)+theme_classic(base_size=13)+labs(x=NULL,y="Normalized Expression",title="SOCS3 in T-cell Subtypes");save_plot(p,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_T_Cell_Subtypes_Violin.png"),10,7);d<-DotPlot(tc,features="SOCS3",group.by="T_Cell_Subtype",dot.scale=8)+scale_color_gradient(low="#2166AC",high="#B2182B")+theme_classic(base_size=13)+labs(x=NULL,y=NULL,color="Average Expression",size="Percent Expressed",title="SOCS3 in T-cell Subtypes");save_plot(d,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_T_Cell_Subtypes_DotPlot.png"),8,6)};saveRDS(tc,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/T_Cells/T_cells_annotated.rds"))}
if("IL10"%in%rownames(bc)){bc$IL10_Status<-ifelse(get_expr(bc,"IL10")>IL10_THRESHOLD,"IL10+","IL10-");if("SOCS3"%in%rownames(bc)){p<-VlnPlot(bc,features="SOCS3",group.by="IL10_Status",pt.size=.08)+theme_classic(base_size=13)+labs(x=NULL,y="Normalized Expression",title="SOCS3 in IL10+ and IL10- B cells");save_plot(p,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_B_Cell_IL10_Status_Violin.png"),9,7);d<-DotPlot(bc,features="SOCS3",group.by="IL10_Status",dot.scale=8)+scale_color_gradient(low="#2166AC",high="#B2182B")+theme_classic(base_size=13)+labs(x=NULL,y=NULL,color="Average Expression",size="Percent Expressed",title="SOCS3 in IL10+ and IL10- B cells");save_plot(d,file.path(PROJECT_DIR,"output/Final_Analysis/Figures/SOCS3_B_Cell_IL10_Status_DotPlot.png"),8,6)};write.csv(as.data.frame(table(bc$IL10_Status,bc$Response)),file.path(PROJECT_DIR,"output/Final_Analysis/Tables/B_Cell_IL10_Response_Cell_Counts.csv"),row.names=FALSE);saveRDS(bc,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/B_Cells/B_cells_IL10_annotated.rds"))}

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
write.csv(wilcox_resp(tumor,wp,"Total Tumor"),file.path(PROJECT_DIR,"output/Final_Analysis/Statistics/WNT_Targets_Tumor_NR_vs_R_Wilcoxon_BH.csv"),row.names=FALSE);if("SOCS3"%in%rownames(tumor))write.csv(wilcox_resp(tumor,"SOCS3","Total Tumor"),file.path(PROJECT_DIR,"output/Final_Analysis/Statistics/SOCS3_Tumor_NR_vs_R_Wilcoxon_BH.csv"),row.names=FALSE)

saveRDS(combined,file.path(PROJECT_DIR,"output/Final_Analysis/Objects/Final_Combined_Annotated_Seurat.rds"))
sink(file.path(PROJECT_DIR,"output/Final_Analysis/SessionInfo.txt"));cat("Final analysis session\n\n");print(sessionInfo());sink()
write.csv(data.frame(Metric=c("Total combined cells","Aneuploid cells","Refined tumor cells"),Value=c(ncol(combined),sum(combined$CNV_status=="aneuploid",na.rm=TRUE),ncol(tumor))),file.path(PROJECT_DIR,"output/Final_Analysis/Tables/Final_Cell_Count_Report.csv"),row.names=FALSE)
cat("\nFINAL PIPELINE COMPLETED\nOutput: ",file.path(PROJECT_DIR,"output/Final_Analysis"),"\n")

