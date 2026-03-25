####################################### Immunophenotyping 

iphen<-Read10X("../Data/Immunophenotyping-of-COVID19-and-influenza-highlights-the-role-of-type-I-interferons-in-development-of-severe-COVID-19/Data/")
iphen<-CreateSeuratObject(counts = iphen, project = "I-phen-covid-influenza",min.cells = 3, min.features = 200)

#Add column of mitochondrial rna
iphen$mitoPercent<- PercentageFeatureSet(iphen, pattern = "^MT-")

#Create Sample column
iphen$Sample<-rownames(iphen@meta.data)

#Split first column into Patient 
split_ids<-strsplit(iphen@meta.data$Sample, "-")
donor_ids<-sapply(split_ids, function(x) x[2])
iphen@meta.data$Sample<-donor_ids
#CHange Sample for Patient
colnames(iphen@meta.data)[5] <- "Patient"


#Subset each 
#Patient 5 5901
iphen5<-subset(iphen, subset= Patient ==5)

#Patient 13 5780
iphen13<-subset(iphen, subset= Patient ==13)

#Patient 14 6169
iphen14<-subset(iphen, subset= Patient ==14)

#Patient 19 5265
iphen19<-subset(iphen, subset= Patient ==19)

merged.iphen<-merge(iphen13, y=c(iphen14, iphen19, iphen5), #23 142 cells
                    add.cell.ids=ls()[80:83])

saveRDS(merged.iphen, "res/merged-iphen-raw.rds")

View(merged.iphen@meta.data)
#23,142 cells

VlnPlot(merged.iphen, features = c("nCount_RNA","nFeature_RNA","mitoPercent"))
merged.iphen<-subset(merged.iphen, subset=mitoPercent >=2 & mitoPercent <= 15 & nCount_RNA > 800)
View(merged.iphen@meta.data)
#17,022

#subset
#Processing
merged.iphen<-NormalizeData(merged.iphen)
merged.iphen<-FindVariableFeatures(merged.iphen, selection.method = "vst", nfeatures = 3000)
merged.iphen<-ScaleData(merged.iphen)
merged.iphen<-RunPCA(merged.iphen)
#ElbowPlot(merged.iphen, ndims = 50)
merged.iphen<-FindNeighbors(merged.iphen, dims=1:25)
merged.iphen<-FindClusters(merged.iphen) 
merged.iphen<-RunUMAP(merged.iphen, dims=1:25, return.model = T)
DimPlot(merged.iphen, label=T, reduction = "umap", group.by = "Patient")
View(merged.iphen@meta.data)

saveRDS(merged.iphen, "merged.iphen.rds")


#integration CCA
int.iphen<-IntegrateLayers(merged.iphen, method = CCAIntegration, 
                           orig.reduction = "pca", new.reduction="integrated.cca.iphen",
                           verbose=FALSE)
int.iphen<-FindNeighbors(int.iphen, reduction = "integrated.cca.iphen", dims=1:25)
int.iphen<-FindClusters(int.iphen)

int.iphen<-RunUMAP(int.iphen, reduction = "integrated.cca.iphen", dims = 1:30, reduction.name = "umap.cca.iphen")

DimPlot(int.iphen, reduction = "umap.cca.iphen", group.by = c("Patient", "seurat_clusters"))


#Plot UMAP in this case by Seurat clusters
DimPlot(int.iphen.ann, reduction ="umap.cca.iphen", label = T, pt.size = 2, label.size = 5, group.by = c("seurat_clusters","predicted.celltype.l2"))
View(int.iphen.ann@meta.data$predicted.celltype.l2)

#Azimuth

#Join Layers joins a merged object 
int.iphen.joined<-JoinLayers(int.iphen)

int.iphen.ann<-RunAzimuth(int.iphen.joined, reference = "pbmcref")
View(int.iphen.ann@meta.data)

#Rename Idents
int.iphen.ann<-RenameIdents(object=int.iphen.ann,"0"="NK","1"="CD14 Mono","2"="TCM","3"="CD8 TEM","4"="CD16 Mono",
                            "5"="CD14 Mono 2","6"="NK 2","7"="B Naive", "8"="MAIT",
                            "9"="B Memory","10"="T Naive","11"="gdT","12"="CD14 Mono",
                            "13"="cDC", "14"="NK T","15"="CD14 Mono 3", "16"="Proliferative",
                            "17"="Platelet", "18"="NK_CD56Bright", "19"="pDC")


int.iphen.ann$idents<-Idents(int.iphen.ann)


#Cell type 
counting.iphen<-table(int.iphen.ann$Patient, int.iphen.ann$predicted.celltype.l2)
iphen.df<-as.data.frame(counting.iphen)

fwrite(iphen.df, "res/dataframes-t-count/iphen-df.txt")
saveRDS(int.iphen.ann, "res/int-iphen-ann.rds")