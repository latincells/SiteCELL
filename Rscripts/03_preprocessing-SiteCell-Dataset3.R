

#SiteCELL3
LC3<-Read10X(data.dir = "~/Data/SiteCELL/LC3")
LC3<-CreateSeuratObject(LC3, project = "LC3", min.cells = 3, min.features = 200)
View(LC3@meta.data)

LC3$mitoPercent<- PercentageFeatureSet(LC3, pattern = "^MT-")
LC3$BARCODE<-rownames(LC3@meta.data)
metadata_LC3<-LC3@meta.data         #12923

freemLC3<-read.table("~/Data/SiteCELL/LC3/best-files/freemuxlet.clust1.samples.gz", header = T)
table(freemLC3$DROPLET.TYPE)
#AMB   DBL   SNG 
#56   564 12362

#Matching barcodes 
matching_rows_LC3<-merge(metadata_LC3, freemLC3, by="BARCODE")
View(matching_rows_LC3)              #12923

#MERGE MATCHING LINES in BOTH FILES
#Add BEST (demuxlet) column to the pool metadata
LC3@meta.data$Droplet <- matching_rows_LC3$DROPLET.TYPE
LC3@meta.data$Sample <- matching_rows_LC3$BEST.GUESS

#Filter out ONLY SINGLETS 
filt.LC3<-subset(LC3, subset=Droplet =="SNG")

#12362

#Count number of cells per sample
table(filt.LC3@meta.data$Sample)

#0,0  1,1  2,2  3,3 
#3676 3298 2809 2556


#split by sample
LC3.1<-subset(filt.LC3, subset = Sample == '0,0')
LC3.2<-subset(filt.LC3, subset = Sample == '1,1')
LC3.3<-subset(filt.LC3, subset = Sample == '2,2')
LC3.4<-subset(filt.LC3, subset = Sample == '3,3')


#Merge all
LC3.merged<-merge(x=LC3.1, y=c(LC3.2, LC3.3, LC3.4), add.cell.ids=ls()[26:29])

#save
saveRDS(LC3.merged, "LC3-merged-raw.rds")
View(LC3.merged@meta.data)
#12,339 cells in total

#subset
LC3.merged<-subset(LC3.merged, subset= mitoPercent >= 2 & mitoPercent <=15)
#11,827 after filter total cells


#READ 
LC3.merged<-NormalizeData(LC3.merged)
LC3.merged<-FindVariableFeatures(LC3.merged)
LC3.merged<-ScaleData(LC3.merged)
LC3.merged<-RunPCA(LC3.merged)
#ElbowPlot(LC3.merged, ndims = 50)
LC3.merged<-FindNeighbors(LC3.merged, reduction="pca", dims=1:35)
LC3.merged<-FindClusters(LC3.merged)
LC3.merged<-RunUMAP(LC3.merged, dims = 1:35, reduction = "pca")

DimPlot(LC3.merged, reduction = "umap", pt.size = 1.8)

#Integration 
int.LC3<-IntegrateLayers(LC3.merged, method = CCAIntegration, 
                          orig.reduction = "pca", new.reduction="integrated.cca", assay = "RNA")

int.LC3<-FindNeighbors(int.LC3, reduction = "integrated.cca", dims=1:35)
int.LC3<-FindClusters(int.LC3)

int.LC3<-RunUMAP(int.LC3, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")


View(int.LC3@meta.data)
int.LC3@meta.data<-int.LC3@meta.data%>%
  rename(Patient = Sample) 

View(int.LC3@meta.data)
#11,756

DimPlot(int.LC3, reduction = "umap.cca",group.by = "Sample")

#JoinLayers
int.LC3.j<-JoinLayers(int.LC3)

#Run Azimuth
int.LC3.ann<-RunAzimuth(int.LC3.j, reference = "pbmcref")
View(int.LC3.ann@meta.data)


#Cell type counting
counting.czi<-table(int.LC3.ann$Patient, int.LC3.ann$predicted.celltype.l2)
czi.df<-as.data.frame(counting.czi)

