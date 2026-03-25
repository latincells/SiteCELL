#Import raw unprocessed Ficoll datasets

setwd("~/thesis")

merged.tr<-readRDS("~/time-resolving-raw/merge-tr.rds")

merged.immuno<-readRDS("merged-immuno-raw.rds")
merged.iphen<-readRDS("merged-iphen-raw.rds")
merged.kawasaki<-readRDS("merged-kawasaki-raw.rds")
merged.aging<-readRDS("merged-aging-raw.rds")

#Import raw unprocessed LatinCells datasets
merged.LC1<-readRDS("merged-pilot-raw.rds")
merged.LC2<-readRDS("merged-1a-raw.rds")
merged.LC3<-readRDS("czi1-merged-raw.rds")

#Adding Dataset ID column to metadata
merged.tr@meta.data<-merged.tr@meta.data%>%
  mutate(Dataset_ID="Ficoll3")

merged.immuno@meta.data<-merged.immuno@meta.data%>%
  mutate(Dataset_ID="Ficoll1")
merged.immuno$Patient<-rownames(merged.immuno@meta.data)
split_ids<-strsplit(merged.immuno@meta.data$Patient, "_")
donor_ids<-sapply(split_ids, function(x) x[1])
merged.immuno@meta.data$Patient<-donor_ids
merged.immuno$mitoPercent<- PercentageFeatureSet(merged.immuno, pattern = "^MT-")

merged.aging@meta.data<-merged.aging@meta.data%>%
  mutate(Dataset_ID="Ficoll4")
merged.aging$Patient<-rownames(merged.aging@meta.data)
split_ids<-strsplit(merged.aging@meta.data$Patient, "_")
donor_ids<-sapply(split_ids, function(x) x[1])
merged.aging@meta.data$Patient<-donor_ids
merged.aging$mitoPercent<- PercentageFeatureSet(merged.aging, pattern = "^MT-")

merged.kawasaki@meta.data<-merged.kawasaki@meta.data%>%
  mutate(Dataset_ID="Ficoll5")
merged.kawasaki$Patient<-rownames(merged.kawasaki@meta.data)
split_ids<-strsplit(merged.kawasaki@meta.data$Patient, "_")
donor_ids<-sapply(split_ids, function(x) x[1])
merged.kawasaki@meta.data$Patient<-donor_ids
merged.kawasaki$mitoPercent<- PercentageFeatureSet(merged.kawasaki, pattern = "^MT-")

merged.iphen@meta.data<-merged.iphen@meta.data%>%
  mutate(Dataset_ID="Ficoll2")

merged.LC1@meta.data<-merged.LC1@meta.data%>%
  mutate(Dataset_ID="LatinCells1")
merged.LC1@meta.data<-merged.LC1@meta.data%>%
  mutate(Patient = coalesce(Patient, Sample)) %>%  # merge into 'Patient'
  select(-Sample)

merged.LC2@meta.data<-merged.LC2@meta.data%>%
  mutate(Dataset_ID="LatinCells2")

merged.LC3@meta.data<-merged.LC3@meta.data%>%
  mutate(Dataset_ID="LatinCells3")
colnames(merged.LC3@meta.data)[7] <- "Patient"


#Merge all previously merged objects (from both Ficoll and LatinCells datasets)
merged.all<-merge(merged.tr, y=c(merged.immuno, merged.aging, merged.iphen, merged.LC1, merged.LC2, merged.kawasaki, merged.LC3))

#For integration

#Import raw unprocessed Ficoll datasets
int.tr<-readRDS("int-tr-raw.rds")
int.immuno<-readRDS("int-immuno-raw.rds")
int.iphen<-readRDS("int-iphen-raw.rds")
int.kawasaki<-readRDS("int-kawasaki-raw.rds")
int.aging<-readRDS("int-aging-raw.rds")

#Import raw unprocessed LatinCells datasets
int.LC1<-readRDS("int-LC1-raw.rds")
int.LC2<-readRDS("int-LC2-raw.rds")
int.LC3<-readRDS("int-LC3-raw.rds")




#Merge all previously integrated objects (from both Ficoll and LatinCells datasets)
merged.integration<-merge(int.LC3, y=c(int.aging, int.immuno.joined, int.iphen, int.kawasaki, int.LC1, int.LC2, int.tr)) 

#PROCESS MERGED
merged.integration<-NormalizeData(merged.integration)
merged.integration<-FindVariableFeatures(merged.integration)
merged.integration<-ScaleData(merged.integration)
merged.integration<-RunPCA(merged.integration)
#ElbowPlot(merged.integration, ndims = 50)
merged.integration<-FindNeighbors(merged.integration, reduction="pca", dims=1:35)
merged.integration<-FindClusters(merged.integration)
merged.integration<-RunUMAP(merged.integration, dims = 1:35, reduction = "pca")

#Integration 
int.all.objects<-IntegrateLayers(merged.integration, method = CCAIntegration, 
                             orig.reduction = "pca", new.reduction="integrated.cca", assay = "RNA")

int.all.objects<-FindNeighbors(int.all.objects, reduction = "integrated.cca", dims=1:35)
int.all.objects<-FindClusters(int.all.objects)

int.all.objects<-RunUMAP(int.all.objects, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")

DimPlot(int.all.objects, reduction = "umap.cca",group.by = c("ident","Dataset"))
View(int.all.objects@meta.data)

saveRDS(int.all.objects, "int.all.objects.rds")

