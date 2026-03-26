library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
#library(tidyverse)
library(stringr)
library(data.table)
#library(hdf5r)
#library(rhdf5)
library(SeuratData)
library(SeuratDisk)
library(BiocManager)
library(multtest)
#library(BPCells)
library(dittoSeq)
library(pacman)
library(Azimuth)



##### Library B23 - PBMCs were isolated with SiteCELL protocol

####Preprocessing
B23<-Read10X("~/B23-matched/filtered_feature_bc_matrix/")
B23<-CreateSeuratObject(counts = B23, min.cells = 3, min.features = 200)

#Add column of mitochondrial rna
B23$mitoPercent<- PercentageFeatureSet(B23, pattern = "^MT-")

#Add BARCODE column to the metadata
B23$BARCODE<-rownames(B23@meta.data)

#Export metadata to a file
metadata_B23<-B23@meta.data

View(metadata_B23) #17570 barcodes

#Import demux best-file
demuxB23<-read.table(file="~/B23/demuxlet-B23.best", header = T) 

#Check for matching barcodes
matching_rows<-merge(metadata_B23, demuxB23, by="BARCODE")
View(matching_rows)
View(matching_rows.2)

matching_rows.2 <- matching_rows %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

ggplot(matching_rows2, aes(x=RD.TOTL))+
  geom_histogram()+
  facet_grid(~BEST)+ 
  labs(title = "The total number of reads overlapping with variant sites for each droplet")+
  scale_x_continuous(limits = c(0, 75000), breaks = seq(0, 75000, 25000))

ggplot(matching_rows2, aes(x=N.SNP))+
  geom_histogram()+
  facet_grid(~ BEST) + 
  labs(title = "Total SNPS overlapping with any read in the droplet")+
  scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 25))


#To identify non matching barcodes 
non_matching1a<-setdiff(metadata_B23$BARCODE, matching_rows$BARCODE)
non_matching1a

[1] "AGAAGCGAGCTCCGAC-1" "AGGTCATTCGCTCTAC-1" "CAGTTCCAGTGCGTCC-1"
[4] "GCTGAATAGGACACTG-1" "TACGCTCAGGTCTGGA-1" "TTACTGTCAGCATGCC-1"

#Filter those non matching barcodes
filt.B23<-subset(B23, subset = BARCODE != "ACCAACAAGTTCCTGA-1") #1
filt.B23<-subset(filt.B23, subset = BARCODE != "ACCATTTGTCACGTGC-1")#2
filt.B23<-subset(filt.B23, subset = BARCODE != "ACGTAGTAGGTCGAGT-1")#3
filt.B23<-subset(filt.B23, subset = BARCODE != "AGGAAATAGGCACCAA-1")#4
filt.B23<-subset(filt.B23, subset = BARCODE != "ATTCATCGTTCTCGTC-1")#5
filt.B23<-subset(filt.B23, subset = BARCODE != "TCATGTTTCAACCGAT-1")#6
filt.B23<-subset(filt.B23, subset = BARCODE != "TCGGGACTCTCACCCA-1")#7
filt.B23<-subset(filt.B23, subset = BARCODE != "TGCCGAGTCGTAGTGT-1")#8
filt.B23<-subset(filt.B23, subset = BARCODE != "TGTTTGTGTAGGGTAC-1")#9

#Add BEST (demuxlet) column to the pool metadata
filt.B23@meta.data$Droplet <- matching_rows.2$Droplet
filt.B23@meta.data$Patient <- matching_rows.2$Patient

#Filter only SNG droplets
filt.B23<-subset(filt.B23, subset=Droplet=="SNG")


table(filt.B23@meta.data$Patient)
LCBR0093 LCBR0094 LCBR0095 LCBR0096 
2340     2120     2748     2339


#SUBSET BY PATIENT
LCBR0093<-subset(filt.B23, subset = Patient == 'LCBR0093')
saveRDS(LCBR0093, file = "~/LCBR0093-B23.rds")


#subset by patient
LCBR0094<-subset(filt.B23, subset = Patient == 'LCBR0094')
saveRDS(LCBR0094, file = "~/LCBR0094-B23.rds")

#subset by patient
LCBR0095<-subset(filt.B23, subset = Patient == 'LCBR0095')
saveRDS(LCBR0095, file = "~/LCBR0095-B23.rds")

#subset by patient
LCBR0096<-subset(filt.B23, subset = Patient == 'LCBR0096')
saveRDS(LCBR0096, file = "~/LCBR0096-B23.rds")

#MERGE PATIENTS
merged.b23<-merge(x=LCBR0093, y=c(LCBR0094,LCBR0095,LCBR0096), 
                 add.cell.ids=ls()[6:9])

#save merged
saveRDS(merged.b23, file = "~/matched-data/RDSs/merged-B23.rds")

merged.b23<-NormalizeData(merged.b23)
merged.b23<-FindVariableFeatures(merged.b23, selection.method = "vst", nfeatures = 3000)
merged.b23<-ScaleData(merged.b23)
merged.b23<-RunPCA(merged.b23)
ElbowPlot(merged.b23, ndims = 50)
merged.b23<-FindNeighbors(merged.b23, dims=1:30)
merged.b23<-FindClusters(merged.b23, resolution = 1) 
merged.b23<-RunUMAP(merged.b23, dims=1:30, return.model = T)

DimPlot(merged.b23, label=T, reduction = "umap", group.by = c("Patient","seurat_clusters"))

merged.b23

int.B23<-IntegrateLayers(merged.b23, method = CCAIntegration, 
                            orig.reduction = "pca", new.reduction="integrated.cca",
                            verbose=FALSE)

int.B23<-FindNeighbors(int.B23, reduction = "integrated.cca", dims=1:30)
int.B23<-FindClusters(int.B23, resolution = 1)

int.B23<-RunUMAP(int.B23, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(int.B23, reduction = "umap.cca",group.by = c("seurat_clusters","RNA_snn_res.0.8", "RNA_snn_res.1"), alpha = 0.5, pt.size = 1.2, label=TRUE)

#Join Layers
#Getting matrix
int.B23.j<-JoinLayers(int.B23)
int.B23.m<-as.matrix(GetAssayData(int.B23.j, layer = 'counts')[, WhichCells(int.B23.j)])
#Once created, RNA matrix needs to be converted into an RDS object.
saveRDS(int.B23.m, file = "~/matched-data/RDSs/B23-matched.rds")


#ANNOTATION AS POOL
#Origin list containing a list of barcodes of particulas interest that will be used as a reference
#to mark those same genes in another list
B23.o<-fread("~/matched-data/azimuth-annotations/azimuth_pred_B23-matched.tsv", stringsAsFactors = TRUE, header = TRUE)

write.table(int.B23.j@meta.data, file = "~/matched-data/B23_metadata.tsv", sep = '\t')
#Target list in which genes are going to be tagged if they have a matching barcode in a shorter list containing genes of interest
B23.t<-fread("~/matched-data/B23_metadata.tsv", stringsAsFactors = TRUE)

#Tagging barcode list
B23.o$Filtered<-"*"
B23.t$Tagged<-NA

B23.t$V1 <- factor(B23.t$V1, levels=c(levels(B23.t$V1), levels(B23.o$Filtered)))
B23.t$Tagged [na.omit(match(B23.o$cell, B23.t$V1))] <- B23.o$Filtered [which(B23.o$cell %in% B23.t$V1)]

#Add Tagged column into the Metadata of the seurat object 
int.B23.j@meta.data$Tagged<-B23.t$Tagged
#Subset seurat object using tagged column
int.B23.j<-subset(int.B23.j, subset=Tagged == "*") #569 cells filtered
#Add azimuth cell types prediction and score
int.B23.j@meta.data$celltype<-B23.o$predicted.celltype.l2
int.B23.j@meta.data$prediction.score<-B23.o$predicted.celltype.l2.score
#filter.pbmc1n2@meta.data$mapping.score<-Olist.joined.cca$mapping.score
#Last check for the metadata
View(int.B23.j@meta.data)





#CALCULATE proportion
table(int.B23.j@meta.data$celltype)

#Percentage of cells
int.B23.prop<-prop.table(table(int.B23.j@meta.data$celltype))

int.B23.j<-RenameIdents(object=int.B23.j,"0"="CD14 Mono","1"="CD4 TCM","2"="CD4 Naive","3"="CD8 TEM",
                         "4"="CD8 TCM","5"="CD8 Naive","6"="NK ","7"="CD4 Naive",
                         "8"="B Naive","9"="MAIT","10"="B Intemediate","11"="CD16 Mono",
                         "12"="Treg", "13"="CD4 Naive","14"="cDC", "15"="Plasmablast",
                         "16"="pDC", "17"="NK proliferative", "18"="NK_CD56bright")

#Add idents
int.B23.j$idents<-Idents(int.B23.j)


DimPlot(int.B23.j, reduction="umap.cca",group.by = "idents", label = TRUE) +NoLegend()

FeaturePlot(int.B23.j, features = c("HBD", "HBM","PPBP", "PF4"))

#Add idents

DefaultAssay(filt.B23)<-"RNA"

int.B23.j@meta.data$mitoPercent<- PercentageFeatureSet(int.B23.j, pattern = "^MT-")

agg.B23<-AggregateExpression(int.B23.j, return.seurat = T,
                                assays = "B23",
                                group.by = c("idents","Patient"))

#normalized expression WITH RNA ASSAY WITH SAMPLE AND CELL ANNOTATION
agg.B23.matrix<-as.matrix(GetAssayData(agg.B23, layer = 'counts')[, WhichCells(agg.B23)])

write.table(agg.B23.matrix, "agg-B23.tsv", sep = "\t")


#Proportions of celltypes

CO3.prop<-subset(int.B23.j, subset= Patient == "LC-CO-0003")
table(CO3.prop@meta.data$celltype)

CO25.prop<-subset(int.B23.j, subset= Patient == "LC-CO-0025")
table(CO25.prop@meta.data$celltype)

MX04.prop<-subset(int.B23.j, subset= Patient == "LCMX0004")
table(MX04.prop@meta.data$celltype)

MX147.prop<-subset(int.B23.j, subset= Patient == "LCMX0147")
table(MX147.prop@meta.data$celltype)

int.B23.ann<-RunAzimuth(int.B23.j, reference = "pbmcref")

View(int.B23.ann@meta.data)

DimPlot(int.B23.ann, reduction="umap.cca",group.by = "predicted.celltype.l2", label = TRUE) +NoLegend()
