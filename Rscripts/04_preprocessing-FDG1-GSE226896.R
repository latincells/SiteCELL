library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(tidyr)
library(dplyr)
#library(tidyverse)
library(stringr)
library(data.table)
#library(hdf5r)
#library(rhdf5)
#library(SeuratData)
#library(SeuratDisk)
library(BiocManager)
#library(BPCells)
library(dittoSeq)
library(Azimuth)

#CLINICAL IMMUNOLOGICAL = immuno

###Donor1
immunop1<-Read10X("../Data/cxcr4-recruitment-lymphoid-cells/uno/")
immunop1<-CreateSeuratObject(counts = immunop1, project = "Immunophenotyping-1", min.cell = 3, min.features = 200)
View(immunop1@meta.data)

###Donor2
immunop2<-Read10X("../Data/cxcr4-recruitment-lymphoid-cells/dos/")
immunop2<-CreateSeuratObject(counts = immunop2, project = "Immunophenotyping-2", min.cell = 3, min.features = 200)
View(immunop2@meta.data)

###Donor3
immunop3<-Read10X("../Data/cxcr4-recruitment-lymphoid-cells/tres/")
immunop3<-CreateSeuratObject(counts = immunop3, project = "Immunophenotyping-3", min.cell = 3, min.features = 200)
View(immunop3@meta.data)

###Donor4
immunop4<-Read10X("../Data/cxcr4-recruitment-lymphoid-cells/cuatro/")
immunop4<-CreateSeuratObject(counts = immunop4, project = "Immunophenotyping-4", min.cell = 3, min.features = 200)
View(immunop4@meta.data)

#Integration

#Starting with merge
merged.immuno<-merge(immunop1, y=c(immunop2, immunop3, immunop4),
                     add.cell.ids=ls()[59:62])

View(merged.immuno@meta.data)
# 60,982 cells

#mito calculation
merged.immuno$mitoPercent<- PercentageFeatureSet(merged.immuno, pattern = "^MT-")
#subset
merged.immuno<-subset(merged.immuno, subset=mitoPercent >=2 & mitoPercent <= 15)

View(merged.immuno@meta.data)
# 57,908

saveRDS(merged.immuno, "merged-immuno-raw.rds")

#Create new Patient column
merged.immuno$Sample<-rownames(merged.immuno@meta.data)

#Split first column into Patient 
split_ids<-strsplit(merged.immuno@meta.data$Sample, "_")
donor_ids<-sapply(split_ids, function(x) x[1])
merged.immuno@meta.data$Sample<-donor_ids
#CHange Sample for Patient
colnames(merged.immuno@meta.data)[4] <- "Patient"


#Object merged.immuno = 23404 cells
#Split by Patient column
#merged.immuno<-SplitObject(merged.immuno, split.by = 'Patient')
merged.immuno<-NormalizeData(merged.immuno)
merged.immuno<-FindVariableFeatures(merged.immuno, selection.method = "vst", nfeatures = 3000)
merged.immuno<-ScaleData(merged.immuno)
merged.immuno<-RunPCA(merged.immuno)
#ElbowPlot(merged.immuno, ndims = 50)
merged.immuno<-FindNeighbors(merged.immuno, dims=1:30)
merged.immuno<-FindClusters(merged.immuno) 
merged.immuno<-RunUMAP(merged.immuno, dims=1:30, return.model = T)
DimPlot(merged.immuno, label=T, reduction = "umap", group.by = "Patient")
View(merged.immuno@meta.data)

#integration CCA
int.immuno<-IntegrateLayers(merged.immuno, method = CCAIntegration, 
                            orig.reduction = "pca", new.reduction="integrated.cca.immuno",
                            verbose=FALSE)
int.immuno<-FindNeighbors(int.immuno, reduction = "integrated.cca.immuno", dims=1:25)
int.immuno<-FindClusters(int.immuno)

int.immuno<-RunUMAP(int.immuno, reduction = "integrated.cca.immuno", dims = 1:25, reduction.name = "umap.cca.immuno")

DimPlot(int.immuno, reduction = "umap.cca.immuno", group.by = c("Patient", "seurat_clusters"))

#Join Layers
int.immuno<-JoinLayers(int.immuno)

int.immuno.ann<-RunAzimuth(int.immuno, reference = "pbmcref")

#Until here, Seurat object will have only the merged-integrated Donor data but 
#does not include celltype annotation by azimuth.
View(int.immuno.ann@meta.data)
#56,810

#Plot UMAP in this case by Seurat clusters
DimPlot(int.immuno.ann, reduction ="umap.cca.immuno", label = T, pt.size = 2, label.size = 3, group.by = c("seurat_clusters","predicted.celltype.l2"))+NoLegend()

#Rename Idents
int.immuno.ann<-RenameIdents(object=int.immuno.ann,"0"="CD8 TM","1"="CD14 Mono","2"="CD4 TCM","3"="Naive",
                             "4"="NK","5"="MAIT","6"="CD16 Mono","7"="NK T","8"="B Naive",
                             "9"="CD4 TCM","10"="Treg","11"="CD14 Mono","12"="CD8 Effector",
                             "13"="B Memory", "14"="CD16 Mono","15"="cDC", "16"="NK_BrightCD56",
                             "17"="CD4 TM", "18"="Platelet", "19"="Proliferative","20"="pDC","21"="Eryth","22"="HSPC",
                             "23"="CD16 Mono 2","24"="Plasmablast")

#add idents
int.immuno.ann$idents<-Idents(int.immuno.ann)

#Dimplot
DimPlot(int.immuno.ann, reduction = "umap.cca.immuno", label = T, pt.size = 2, label.size = 5, group.by = "idents")

#Patient
int.immuno.ann$Sample<-rownames(int.immuno.ann@meta.data)
int.immuno.ann@meta.data<-separate(int.immuno.ann@meta.data, col='Sample', into = c('Patient', 'Barcode'), 
                                   sep='_')


#Creacion y conteo de tipos celulares
conteo.immuno<-table(int.immuno.ann$Patient, int.immuno.ann$predicted.celltype.l2)
immuno.df<-as.data.frame(conteo.immuno)

fwrite(immuno.df, "res/dataframes-t-count/immuno-df.txt")

saveRDS(int.immuno.ann, "res/int-annotated-datasets/int-immuno-ann.rds")