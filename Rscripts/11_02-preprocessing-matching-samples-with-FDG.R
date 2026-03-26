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
#library(multtest)
#library(BPCells)
library(dittoSeq)
library(pacman)
library(sctransform)



##### Library B24, PBMCS were isolated with FDG protocol. 

B24.seu<-Read10X("~/filtered_feature_bc_matrix/")
B24.seu.seu<-CreateSeuratObject(counts = B24.seu, min.cells = 3, min.features = 200)

#No crea las columnas de NCounts y nFeatures así que uso este código para forzar a que lo haga
B24.seu.seu[["nCount_RNA"]] <- Matrix::colSums(GetAssayData(B24.seu.seu, layer = "counts"))
B24.seu.seu[["nFeature_RNA"]] <- Matrix::colSums(GetAssayData(B24.seu.seu, layer = "counts") > 0)


#Add column of mitochondrial rna
B24.seu.seu$mitoPercent<- PercentageFeatureSet(B24.seu.seu, pattern = "^MT-")

#Add BARCODE column to the metadata
B24.seu$BARCODE<-rownames(B24.seu@meta.data)

#Export metadata to a file
metadata_B24.seu<-B24.seu@meta.data

View(metadata_B24.seu) #17570 barcodes

#Import demux best-file
demuxB24.seu<-read.table(file="~/demuxlet-B24.best", header = T) 

#Check for matching barcodes
matching_rows<-merge(metadata_B24.seu, demuxB24.seu, by="BARCODE")
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
non_matching1a<-setdiff(metadata_B24.seu$BARCODE, matching_rows$BARCODE)
non_matching1a

[1] "AGAAGCGAGCTCCGAC-1" "AGGTCATTCGCTCTAC-1" "CAGTTCCAGTGCGTCC-1"
[4] "GCTGAATAGGACACTG-1" "TACGCTCAGGTCTGGA-1" "TTACTGTCAGCATGCC-1"

#Filter those non matching barcodes
filt.B24<-subset(B24.seu, subset = BARCODE != "ACGTCCTGTCTTTCAT-1") #1
filt.B24<-subset(filt.B24, subset = BARCODE != "ATTACTCGTTAAGAAC-1")#2
filt.B24<-subset(filt.B24, subset = BARCODE != "CAATCGAGTCAGTTTG-1")#3
filt.B24<-subset(filt.B24, subset = BARCODE != "CATCCGTGTAGAGATT-1")#4
filt.B24<-subset(filt.B24, subset = BARCODE != "CATGCCTGTCAATGGG-1")#5
filt.B24<-subset(filt.B24, subset = BARCODE != "CTACGGGTCGCAGTTA-1")#6
filt.B24<-subset(filt.B24, subset = BARCODE != "CTATCCGTCGGCTGTG-1")#7
filt.B24<-subset(filt.B24, subset = BARCODE != "CTCCCAAGTCGACGCT-1")#8
filt.B24<-subset(filt.B24, subset = BARCODE != "GATCATGCAGGTGTGA-1")#9
filt.B24<-subset(filt.B24, subset = BARCODE != "GATCGTAAGTGGCCTC-1")#10
filt.B24<-subset(filt.B24, subset = BARCODE != "GATTGGTTCGAGCCTG-1")#11
filt.B24<-subset(filt.B24, subset = BARCODE != "GCATCGGTCTGGTCAA-1")#12
filt.B24<-subset(filt.B24, subset = BARCODE != "GCCAGGTTCTTCCAGC-1")#13
filt.B24<-subset(filt.B24, subset = BARCODE != "GTCCTCACATGACTCA-1")#14
filt.B24<-subset(filt.B24, subset = BARCODE != "TTGTTCACAGACCTGC-1")#15


#Add BEST (demuxlet) column to the pool metadata
filt.B24@meta.data$Droplet <- matching_rows.2$Droplet
filt.B24@meta.data$Patient <- matching_rows.2$Patient

#Filter only SNG droplets
filt.B24<-subset(filt.B24, subset=Droplet=="SNG")


table(filt.B24@meta.data$Patient)
LCBR0093 LCBR0094 LCBR0095 LCBR0096 
3663     3058     2758     2385 


##########                                            FREEMUXLET
freemp24<-read.table(file="~/B24.clust1.samples.gz", header = T)
View(freemp24) ###17075

matching_barcodes<-merge(matching_rows.2, freemp24, by="BARCODE", all=TRUE)
View(matching_barcodes)

matching_barcodes<-subset(matching_barcodes, subset=Droplet=="SNG")
View(matching_barcodes) # 12,390

LCBR0093<-subset(matching_barcodes, subset = Patient=="LCBR0093")
table(LCBR0093$BEST.GUESS) #========================= 1
#0,0  1,0  1,1  2,1  2,2  3,1 
#4    6 3639    4    5    5 

LCBR0094<-subset(matching_barcodes, subset = Patient=="LCBR0094")
table(LCBR0094$BEST.GUESS) #================================2
# 0,0  1,1  2,0  2,1  2,2  3,2  3,3 
#2    6    7    8 3029    4    2 

LCBR0095<-subset(matching_barcodes, subset = Patient=="LCBR0095")
table(LCBR0095$BEST.GUESS) #================================0
#0,0  1,0  1,1  2,0  2,2  3,0 
#2736    7    6    6    2    1

LCBR0096<-subset(matching_barcodes, subset = Patient=="LCBR0096")
table(LCBR0096$BEST.GUESS) #================================3
#0,0  1,1  2,2  3,0  3,1  3,3 
#2    8    2    3    5 2365 


#Merge metadata and freemux outputs by "BARCODE" column

matching1b.freemxmetadata<-merge(freemp24, metadata_B24.seu, by="BARCODE", all=TRUE)
matching1b.freemxmetadata<-subset(matching1b.freemxmetadata, subset=orig.ident =="SeuratProject")

#Merge previously merged data and demux4b by "BARCODE" column
matching.freemxmetaxdemux<-merge(matching1b.freemxmetadata, demuxB24.seu, by="BARCODE", all=TRUE)
matching.freemxmetaxdemux<-subset(matching.freemxmetaxdemux, subset=orig.ident=="SeuratProject")

#Merge only pool metadata and freemux file, keep only SNG and them assign genotypes
matching.freemxmetaxdemux<-matching.freemxmetaxdemux %>%
  mutate(Patient=NA) %>%
  mutate(Patient =ifelse(BEST.GUESS == "0,0","LCBR0095", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "1,1","LCBR0093", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "2,2","LCBR0094", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "3,3","LCBR0096", Patient))


table(matching.freemxmetaxdemux$Patient)
#LCBR0093 LCBR0094 LCBR0095 LCBR0096 
#3850     3187     2896     2466

#How many SNG do we have FROM FREEMUXLET?
table(matching.freemxmetaxdemux$DROPLET.TYPE)
#AMB   DBL   SNG 
#53   808 12346 

#Add columns to the metadata of the Seurat object
B24.seu@meta.data$Droplet<-matching.freemxmetaxdemux$DROPLET.TYPE
B24.seu@meta.data$Patient<-matching.freemxmetaxdemux$Patient
B24.seu@meta.data$demuxSNP<-matching.freemxmetaxdemux$N.SNP
B24.seu@meta.data$freemuxSNP<-matching.freemxmetaxdemux$NUM.SNPS

#Delete "NAs"
B24.seu@meta.data<-B24.seu@meta.data %>% drop_na() ###  12,384


#SUBSET BY PATIENT
LCBR0093<-subset(B24.seu, subset = Patient == 'LCBR0093')
saveRDS(LCBR0093, file = "~/LCBR0093-B24.seu.rds")


#subset by patient
LCBR0094<-subset(B24.seu, subset = Patient == 'LCBR0094')
saveRDS(LCBR0094, file = "~/LCBR0094-B24.seu.rds")

#subset by patient
LCBR0095<-subset(B24.seu, subset = Patient == 'LCBR0095')
saveRDS(LCBR0095, file = "~/LCBR0095-B24.seu.rds")

#subset by patient
LCBR0096<-subset(B24.seu, subset = Patient == 'LCBR0096')
saveRDS(LCBR0096, file = "~/LCBR0096-B24.seu.rds")

#MERGE PATIENTS
merged.B24.seu<-merge(x=LCBR0093, y=c(LCBR0094,LCBR0095,LCBR0096))


merged.B24.seu<-NormalizeData(merged.B24.seu)
merged.B24.seu<-FindVariableFeatures(merged.B24.seu, selection.method = "vst", nfeatures = 3000)
merged.B24.seu<-ScaleData(merged.B24.seu)
merged.B24.seu<-RunPCA(merged.B24.seu)
#ElbowPlot(merged.B24.seu, ndims = 50)
merged.B24.seu<-FindNeighbors(merged.B24.seu, dims=1:30)
merged.B24.seu<-FindClusters(merged.B24.seu, resolution = 1) 
merged.B24.seu<-RunUMAP(merged.B24.seu, dims=1:30, return.model = T)

DimPlot(merged.B24.seu, label=T, reduction = "umap", group.by = c("Patient","seurat_clusters"))

merged.B24.seu

int.B24.seu<-IntegrateLayers(merged.B24.seu, method = CCAIntegration, 
                         orig.reduction = "pca", new.reduction="integrated.cca",
                         verbose=FALSE)

int.B24.seu<-FindNeighbors(int.B24.seu, reduction = "integrated.cca", dims=1:30)
int.B24.seu<-FindClusters(int.B24.seu, resolution = 1.2)

int.B24.seu<-RunUMAP(int.B24.seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(int.B24.seu, reduction = "umap.cca",group.by =c("RNA_snn_res.0.8","RNA_snn_res.1","RNA_snn_res.1.2"), alpha = 0.5, pt.size = 1.2, label=TRUE)

FeaturePlot(int.B24.seu, features = c("HBD", "HBM","PPBP", "PF4"), reduction = "umap.cca")

#Join Layers
#Getting matrix
int.B24.seu.j<-JoinLayers(int.B24.seu)
int.B24.seu.ann<-RunAzimuth(int.B24.seu.j, reference = "pbmcref")

DimPlot(int.B24.seu.ann, reduction = "umap.cca",group.by =c("seurat_clusters","predicted.celltype.l2"), alpha = 0.5, pt.size = 1.2, label=TRUE)+NoLegend()

#CALCULATE proportion
table(int.B24.seu.j@meta.data$celltype)

#Percentage of cells
int.B24.seu.prop<-prop.table(table(int.B24.seu.j@meta.data$celltype))

int.B24.seu.j<-RenameIdents(object=int.B24.seu.j,"0"="CD4 TCM","1"="CD4 TCM","2"="CD4 Naive","3"="CD8 Naive",
                        "4"="CD14 Mono","5"="CD14 Mono","6"="NK","7"="CD8 TEM",
                        "8"="CD8 TEM","9"="B Naive","10"="B Intermediate","11"="MAIT",
                        "12"="CD8 TEM", "13"="CD16 Mono","14"="Treg", "15"="CD4 TCM",
                        "16"="Plasmablast", "17"="cDC", "18"="CD8 TEM","19"="CD4 TEM","20"="CD4 TCM", "21"="CD14 Mono", "22"="pDC",
                        "23"=" ","24"="HSCP")

#Add idents
int.B24.seu.j$idents<-Idents(int.B24.seu.j)


DimPlot(int.B24.seu.j, reduction="umap.cca", label = TRUE) +NoLegend()
