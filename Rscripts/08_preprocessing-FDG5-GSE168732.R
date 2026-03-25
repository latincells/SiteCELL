
#Importing data

#Donor1 
kawasaki1<-Read10X("../Data/Single-cell-sequencing-of-peripheral-blood-mononuclear-cells-in-acute-Kawasaki-disease/HC1/")
kawasaki1<-CreateSeuratObject(counts=kawasaki1, min.cells = 3, min.features = 200)
View(kawasaki1@meta.data)

#Donor2
kawasaki2<-Read10X("../Data/Single-cell-sequencing-of-peripheral-blood-mononuclear-cells-in-acute-Kawasaki-disease/HC2/")
kawasaki2<-CreateSeuratObject(counts=kawasaki2, min.cells = 3, min.features = 200)

#Donor3
kawasaki3<-Read10X("../Data/Single-cell-sequencing-of-peripheral-blood-mononuclear-cells-in-acute-Kawasaki-disease/HC3/")
kawasaki3<-CreateSeuratObject(counts=kawasaki3, min.cells = 3, min.features = 200)

###################
################### kawasaki1 kawasaki2 kawasaki3 
################### 5690      5077      5441 
###################

#Starting with merge
merged.kawasaki<-merge(kawasaki1, y=c(kawasaki2, kawasaki3),
                       add.cell.ids=ls()[81:83])

saveRDS(merged.kawasaki,"merged-kawasaki-raw.rds")

#calculate mitochondrial percent
merged.kawasaki$mitoPercent<-PercentageFeatureSet(merged.kawasaki, pattern = '^MT-')
merged.kawasaki$Sample<-rownames(merged.kawasaki@meta.data)
merged.kawasaki@meta.data<-separate(merged.kawasaki@meta.data, col='Sample', into = c('Patient', 'Barcode'), 
                                    sep='_')

#Object merged.kawasaki = 16,208 cells 

#subset
merged.kawasaki<-subset(merged.kawasaki, subset = mitoPercent >= 2 & mitoPercent <= 15)
#after subset
#12,864


merged.kawasaki<-NormalizeData(merged.kawasaki, )#Can be CPM, log-norm or centered-log ratio 
merged.kawasaki<-FindVariableFeatures(merged.kawasaki, selection.method = "vst", nfeatures = 3000)
merged.kawasaki<-ScaleData(merged.kawasaki)
merged.kawasaki<-RunPCA(merged.kawasaki)
#ElbowPlot(merged.kawasaki, ndims = 50)
merged.kawasaki<-FindNeighbors(merged.kawasaki, dims=1:30)
merged.kawasaki<-FindClusters(merged.kawasaki) 
merged.kawasaki<-RunUMAP(merged.kawasaki, dims=1:30, return.model = T)
DimPlot(merged.kawasaki, label=T, reduction = "umap", group.by = c("Patient","seurat_clusters"))

View(merged.kawasaki@meta.data)

#integration
int.kawasaki<-IntegrateLayers(merged.kawasaki, method = CCAIntegration, 
                              orig.reduction = "pca", new.reduction="integrated.cca.ka",
                              verbose=FALSE)
int.kawasaki<-FindNeighbors(int.kawasaki, reduction = "integrated.cca.ka", dims=1:30)
int.kawasaki<-FindClusters(int.kawasaki)
int.kawasaki<-RunUMAP(int.kawasaki, reduction = "integrated.cca.ka", dims = 1:30, reduction.name = "umap.cca.ka")
DimPlot(int.kawasaki, reduction = "umap.cca.ka", label = T ,group.by = c("seurat_clusters", "Patient" ))+NoLegend()


#Extract RNA matrix. Necessary to provide RNA assay to Azimuth
#Join Layers joins a merged object 
int.kawasaki.j<-JoinLayers(int.kawasaki)

#Azimuth annotation
int.kawasaki.ann<-RunAzimuth(int.kawasaki.j, reference = "pbmcref")

View(int.kawasaki.ann@meta.data)
#Plot UMAP in this case by Seurat clusters
DimPlot(int.kawasaki.ann, label = T, pt.size = 1.5, label.size = 5, group.by =c("seurat_clusters","predicted.celltype.l2"))

#Rename Idents
int.kawasaki.ann<-RenameIdents(object=int.kawasaki.ann,"0"="CD4 Naive","1"="CD8 Naive","2"="B Naive","3"="NK",
                               "4"="CD8 TEM","5"="CD4 TCM","6"="CD14 Mono","7"="B Memory",
                               "8"="gdT","9"="B Intermediate","10"="CD8 T effector","11"="CD16 Mono",
                               "12"="B Naive", "13"="Proliferative","14"="NK-CD56 Bright", "15"="dnT",
                               "16"="CD4 Naive", "17"="cDC","18"="B Naive", "19"="pDC")


#Add Idents to metadata
int.kawasaki.ann$idents<-Idents(int.kawasaki.ann)


#Cell type counting
conteo.kawa<-table(int.kawasaki.ann$Patient, int.kawasaki.ann$predicted.celltype.l2)
kawa.df<-as.data.frame(conteo.kawa)

DimPlot(int.kawasaki, reduction = "umap.cca.ka", label = T, pt.size = 2, label.size = 5, group.by = "idents")+NoLegend()

