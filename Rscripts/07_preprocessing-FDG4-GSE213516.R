
#Importing data

###Donor 1
HC1<-Read10X("../Data/Human-pbmc-scrna-seq-based-aging-clocks-reveal-ribosome-to-inflammation-balance-as-a-single-cell-aging/HC1/")
HC1<-CreateSeuratObject(counts = HC1,min.cells = 3, min.features = 200)

###Donor 2
HC2<-Read10X("../Data/Human-pbmc-scrna-seq-based-aging-clocks-reveal-ribosome-to-inflammation-balance-as-a-single-cell-aging/HC2/")
HC2<-CreateSeuratObject(counts = HC2, min.cells = 3, min.features = 200)

###Donor 3
HC3<-Read10X("../Data/Human-pbmc-scrna-seq-based-aging-clocks-reveal-ribosome-to-inflammation-balance-as-a-single-cell-aging/HC3/")
HC3<-CreateSeuratObject(counts = HC3, min.cells = 3, min.features = 200)

###Donor 4
HC4<-Read10X("../Data/Human-pbmc-scrna-seq-based-aging-clocks-reveal-ribosome-to-inflammation-balance-as-a-single-cell-aging/HC4/")
HC4<-CreateSeuratObject(counts = HC4, min.cells = 3, min.features = 200)

###Donor 5
HC5<-Read10X("../Data/Human-pbmc-scrna-seq-based-aging-clocks-reveal-ribosome-to-inflammation-balance-as-a-single-cell-aging/HC5/")
HC5<-CreateSeuratObject(counts = HC5, min.cells = 3, min.features = 200)

#Integration

#Starting with merge
merged.aging<-merge(HC1, y=c(HC2,HC3,HC4,HC5),
                    add.cell.ids=ls()[53:57])

View(merged.aging@meta.data)
saveRDS(merged.aging, "merged-aging-raw.rds")

#42,771 cells

#calculate mitopercent
merged.aging$mitoPercent<-PercentageFeatureSet(merged.aging, pattern = '^MT-')

#Patient
merged.aging$Sample<-rownames(merged.aging@meta.data)
merged.aging@meta.data<-separate(merged.aging@meta.data, col='Sample', into = c('Patient', 'Barcode'), 
                                 sep='_')

#subset
merged.aging<-subset(merged.aging, subset=mitoPercent >=2 & mitoPercent <= 15)
#39,291 cells after mito filter

merged.aging<-NormalizeData(merged.aging, )#Can be CPM, log-norm or centered-log ratio 
merged.aging<-FindVariableFeatures(merged.aging, selection.method = "vst", nfeatures = 3000)
merged.aging<-ScaleData(merged.aging)
merged.aging<-RunPCA(merged.aging)
#ElbowPlot(merged.aging, ndims = 50)
merged.aging<-FindNeighbors(merged.aging, dims=1:25)
merged.aging<-FindClusters(merged.aging) 
merged.aging<-RunUMAP(merged.aging, dims=1:25, return.model = T)
DimPlot(merged.aging, label=T, reduction = "umap", group.by = c("Patient","seurat_clusters"))

#integration CCA
int.aging<-IntegrateLayers(merged.aging, method = CCAIntegration, 
                           orig.reduction = "pca", new.reduction="integrated.cca",
                           verbose=FALSE)
int.aging<-FindNeighbors(int.aging, reduction = "integrated.cca", dims=1:25)
int.aging<-FindClusters(int.aging)
int.aging<-RunUMAP(int.aging, reduction = "integrated.cca", dims = 1:25, reduction.name = "umap.cca")
DimPlot(int.aging, reduction = "umap.cca", group.by = c("seurat_clusters", "celltype","Patient" ))+NoLegend()


#JoinLayers
int.aging.joined<-JoinLayers(int.aging)

#azimuth
int.aging.ann<-RunAzimuth(int.aging.joined, reference = "pbmcref")


#Plot UMAP in this case by Seurat clusters
DimPlot(int.aging.ann, reduction = "umap.cca", label = T, pt.size = 2, label.size = 5, group.by =c("predicted.celltype.l2","seurat_clusters"))


#Rename Idents#Correr de nuevo, modificaciones hechas
int.aging.ann<-RenameIdents(object=int.aging.ann,"0"="CD14 Mono","1"="NK","2"="TCM","3"="NK",
                            "4"="CD4 Naive","5"="CD8 Effector","6"="CD8 TEM","7"="B Naive",
                            "8"="CD8 Naive","9"="MAIT","10"="CD16 Mono","11"="CD14 Mono 2",
                            "12"="CD4 TCM", "13"="B Intermediate","14"="Platelet", "15"="CD14 Mono 3",
                            "16"="cDC", "17"="Platelet", "18"="NK_CD56bright","19"="CD14 Mono 3","20"="B Memory","21"="Proliferative",
                            "22"="pDC","23"="CD14 Mono 4")


#Add Idents to metadata
int.aging.ann$idents<-Idents(int.aging.ann)

#Cell type counting
conteo.aging<-table(int.aging.ann$Patient, int.aging.ann$predicted.celltype.l2)
aging.df<-as.data.frame(conteo.aging)

fwrite(aging.df, "res/dataframes-t-count/aging-df.txt")
saveRDS(int.aging.ann, "res/int-annotated-datasets/int-aging-ann.rds")
