
##Import data
lane1<-Read10X_h5(("GSM4929209_B3_10xlane01_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane2<-Read10X_h5(("GSM4929210_B3_10xlane02_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane3<-Read10X_h5(("GSM4929211_B3_10xlane03_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane4<-Read10X_h5(("GSM4929212_B3_10xlane04_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane5<-Read10X_h5(("GSM4929213_B3_10xlane05_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane6<-Read10X_h5(("GSM4929214_B3_10xlane06_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane7<-Read10X_h5(("GSM4929215_B3_10xlane07_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane8<-Read10X_h5(("GSM4929216_B3_10xlane08_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane9<-Read10X_h5(("GSM4929217_B3_10xlane09_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane10<-Read10X_h5(("GSM4929218_B3_10xlane10_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane11<-Read10X_h5(("GSM4929219_B3_10xlane11_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane12<-Read10X_h5(("GSM4929220_B3_10xlane12_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane13<-Read10X_h5(("GSM4929221_B3_10xlane13_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane14<-Read10X_h5(("GSM4929222_B3_10xlane14_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane15<-Read10X_h5(("GSM4929223_B3_10xlane15_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane16<-Read10X_h5(("GSM4929224_B3_10xlane16_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane17<-Read10X_h5(("GSM4929225_B3_10xlane17_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane18<-Read10X_h5(("GSM4929226_B3_10xlane18_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane19<-Read10X_h5(("GSM4929227_B3_10xlane19_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane20<-Read10X_h5(("GSM4929228_B3_10xlane20_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane21<-Read10X_h5(("GSM4929229_B3_10xlane21_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane22<-Read10X_h5(("GSM4929230_B3_10xlane22_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane23<-Read10X_h5(("GSM4929231_B3_10xlane23_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
lane24<-Read10X_h5(("GSM4929232_B3_10xlane24_RNA_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)

#Layer filtering
lane1[['Antibody Capture']]<-NULL
lane2[['Antibody Capture']]<-NULL
lane3[['Antibody Capture']]<-NULL
lane4[['Antibody Capture']]<-NULL
lane5[['Antibody Capture']]<-NULL
lane6[['Antibody Capture']]<-NULL
lane7[['Antibody Capture']]<-NULL
lane8[['Antibody Capture']]<-NULL
lane9[['Antibody Capture']]<-NULL
lane10[['Antibody Capture']]<-NULL
lane11[['Antibody Capture']]<-NULL
lane12[['Antibody Capture']]<-NULL
lane13[['Antibody Capture']]<-NULL
lane14[['Antibody Capture']]<-NULL
lane15[['Antibody Capture']]<-NULL
lane16[['Antibody Capture']]<-NULL
lane17[['Antibody Capture']]<-NULL
lane18[['Antibody Capture']]<-NULL
lane19[['Antibody Capture']]<-NULL
lane20[['Antibody Capture']]<-NULL
lane21[['Antibody Capture']]<-NULL
lane22[['Antibody Capture']]<-NULL
lane23[['Antibody Capture']]<-NULL
lane24[['Antibody Capture']]<-NULL


#Processing
lane1<-CreateSeuratObject(counts = lane1, project = "tr-raw1", min.cell = 3, min.features = 150, assay = 'RNA')
#lane1@assays$RNA$`counts.Antibody Capture`<-NULL

lane2<-CreateSeuratObject(counts = lane2, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane2@assays$RNA$`counts.Antibody Capture`<-NULL

lane3<-CreateSeuratObject(counts = lane3, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane3@assays$RNA$`counts.Antibody Capture`<-NULL

lane4<-CreateSeuratObject(counts = lane4, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane4@assays$RNA$`counts.Antibody Capture`<-NULL

lane5<-CreateSeuratObject(counts = lane5, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane5@assays$RNA$`counts.Antibody Capture`<-NULL

lane6<-CreateSeuratObject(counts = lane6, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane6@assays$RNA$`counts.Antibody Capture`<-NULL

lane7<-CreateSeuratObject(counts = lane7, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane7@assays$RNA$`counts.Antibody Capture`<-NULL

lane8<-CreateSeuratObject(counts = lane8, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane8@assays$RNA$`counts.Antibody Capture`<-NULL

lane9<-CreateSeuratObject(counts = lane9, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane9@assays$RNA$`counts.Antibody Capture`<-NULL

lane10<-CreateSeuratObject(counts = lane10, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane10@assays$RNA$`counts.Antibody Capture`<-NULL

lane11<-CreateSeuratObject(counts = lane11, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane11@assays$RNA$`counts.Antibody Capture`<-NULL

lane12<-CreateSeuratObject(counts = lane12, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane12@assays$RNA$`counts.Antibody Capture`<-NULL

lane13<-CreateSeuratObject(counts = lane13, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane13@assays$RNA$`counts.Antibody Capture`<-NULL

lane14<-CreateSeuratObject(counts = lane14, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane14@assays$RNA$`counts.Antibody Capture`<-NULL

lane15<-CreateSeuratObject(counts = lane15, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane15@assays$RNA$`counts.Antibody Capture`<-NULL

lane16<-CreateSeuratObject(counts = lane16, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane16@assays$RNA$`counts.Antibody Capture`<-NULL

lane17<-CreateSeuratObject(counts = lane17, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane17@assays$RNA$`counts.Antibody Capture`<-NULL

lane18<-CreateSeuratObject(counts = lane18, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane18@assays$RNA$`counts.Antibody Capture`<-NULL

lane19<-CreateSeuratObject(counts = lane19, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane19@assays$RNA$`counts.Antibody Capture`<-NULL

lane20<-CreateSeuratObject(counts = lane20, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane20@assays$RNA$`counts.Antibody Capture`<-NULL

lane21<-CreateSeuratObject(counts = lane21, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane21@assays$RNA$`counts.Antibody Capture`<-NULL

lane22<-CreateSeuratObject(counts = lane22, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane22@assays$RNA$`counts.Antibody Capture`<-NULL

lane23<-CreateSeuratObject(counts = lane23, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane23@assays$RNA$`counts.Antibody Capture`<-NULL

lane24<-CreateSeuratObject(counts = lane24, project = "tr-raw1", min.cell = 3, min.features = 150)
#lane24@assays$RNA$`counts.Antibody Capture`<-NULL
##########Demuxlet

#Demultiplexing
dem1<-fread("GSE161918_demuxletResults/B3_10xlane01_demuxed.best")
dem2<-fread("GSE161918_demuxletResults/B3_10xlane02_demuxed.best")
dem3<-fread("GSE161918_demuxletResults/B3_10xlane03_demuxed.best")
dem4<-fread("GSE161918_demuxletResults/B3_10xlane04_demuxed.best")
dem5<-fread("GSE161918_demuxletResults/B3_10xlane05_demuxed.best")
dem6<-fread("GSE161918_demuxletResults/B3_10xlane06_demuxed.best")
dem7<-fread("GSE161918_demuxletResults/B3_10xlane07_demuxed.best")
dem8<-fread("GSE161918_demuxletResults/B3_10xlane08_demuxed.best")
dem9<-fread("GSE161918_demuxletResults/B3_10xlane09_demuxed.best")
dem10<-fread("GSE161918_demuxletResults/B3_10xlane10_demuxed.best")
dem11<-fread("GSE161918_demuxletResults/B3_10xlane11_demuxed.best")
dem12<-fread("GSE161918_demuxletResults/B3_10xlane12_demuxed.best")
dem13<-fread("GSE161918_demuxletResults/B3_10xlane13_demuxed.best")
dem14<-fread("GSE161918_demuxletResults/B3_10xlane14_demuxed.best")
dem15<-fread("GSE161918_demuxletResults/B3_10xlane15_demuxed.best")
dem16<-fread("GSE161918_demuxletResults/B3_10xlane16_demuxed.best")
dem17<-fread("GSE161918_demuxletResults/B3_10xlane17_demuxed.best")
dem18<-fread("GSE161918_demuxletResults/B3_10xlane18_demuxed.best")
dem19<-fread("GSE161918_demuxletResults/B3_10xlane19_demuxed.best")
dem20<-fread("GSE161918_demuxletResults/B3_10xlane20_demuxed.best")
dem21<-fread("GSE161918_demuxletResults/B3_10xlane21_demuxed.best")
dem22<-fread("GSE161918_demuxletResults/B3_10xlane22_demuxed.best")
dem23<-fread("GSE161918_demuxletResults/B3_10xlane23_demuxed.best")
dem24<-fread("GSE161918_demuxletResults/B3_10xlane24_demuxed.best")



#List of lane object and demultiplexing files
lanes <- list(lane1, lane2, lane3, lane4, lane5, lane6, lane7, lane8, lane9, lane10,
              lane11, lane12, lane13, lane14, lane15, lane16, lane17, lane18, lane19,
              lane20, lane21, lane22, lane23, lane24)

demultiplexing_files <- list(dem1, dem2, dem3, dem4, dem5, dem6, dem7, dem8, dem9, dem10,
                             dem11, dem12, dem13, dem14, dem15, dem16, dem17, dem18, dem19,
                             dem20, dem21, dem22, dem23, dem24)

#List to storage results
matching_rows_all <- list()
metadata_lane<-list()

#Loop to process lane and demultiplexing files
for (i in 1:length(lanes)) {
  
  #Select lane and dem files
  lane <- lanes[[i]]
  lane$BARCODE<-rownames(lane@meta.data)
  lane$mitoPercent<-PercentageFeatureSet(lane, pattern = '^MT-')
  metadata_lane[[i]]<-lane@meta.data
  
  print(paste("Lane metadata", i))
  print(head(metadata_lane[[i]]))
  
  dem <- demultiplexing_files[[i]]
  
  #check content before merge
  print(paste("Lane", i, "metadata", nrow(metadata_lane[[i]]), "lines and", ncol(metadata_lane[[i]]), "columns"))
  print(paste("Lane", i, "demultiplexing", nrow(dem), "lines and", ncol(dem), "columns."))
  
  }
  
  
  #Check barcodes match from dem and metadata files
  matching_rows <- merge(metadata_lane[[i]], demultiplexing_files[[i]], by="BARCODE", all = TRUE)
  matching_rows_all[[i]]<-matching_rows
  

lane$BARCODE<-trimws(rownames(lane@meta.data))
dem$BARCODE<-trimws(dem$BARCODE)  

print(colnames(metadata_lane[[i]]))
print(colnames(dem))

#Check fot filtered lane 
lane1$mitoPercent<-matching_rows_all[[1]]$mitoPercent
lane1$Droplet<-matching_rows_all[[1]]$DROPLET.TYPE
lane1$Patient<-matching_rows_all[[1]]$BEST.GUESS

lane2$mitoPercent<-matching_rows_all[[2]]$mitoPercent
lane2$Droplet<-matching_rows_all[[2]]$DROPLET.TYPE
lane2$Patient<-matching_rows_all[[2]]$BEST.GUESS

lane3$mitoPercent<-matching_rows_all[[3]]$mitoPercent
lane3$Droplet<-matching_rows_all[[3]]$DROPLET.TYPE
lane3$Patient<-matching_rows_all[[3]]$BEST.GUESS

lane4$mitoPercent<-matching_rows_all[[4]]$mitoPercent
lane4$Droplet<-matching_rows_all[[4]]$DROPLET.TYPE
lane4$Patient<-matching_rows_all[[4]]$BEST.GUESS

lane5$mitoPercent<-matching_rows_all[[5]]$mitoPercent
lane5$Droplet<-matching_rows_all[[5]]$DROPLET.TYPE
lane5$Patient<-matching_rows_all[[5]]$BEST.GUESS

lane6$mitoPercent<-matching_rows_all[[6]]$mitoPercent
lane6$Droplet<-matching_rows_all[[6]]$DROPLET.TYPE
lane6$Patient<-matching_rows_all[[6]]$BEST.GUESS

lane7$mitoPercent<-matching_rows_all[[7]]$mitoPercent
lane7$Droplet<-matching_rows_all[[7]]$DROPLET.TYPE
lane7$Patient<-matching_rows_all[[7]]$BEST.GUESS

lane8$mitoPercent<-matching_rows_all[[8]]$mitoPercent
lane8$Droplet<-matching_rows_all[[8]]$DROPLET.TYPE
lane8$Patient<-matching_rows_all[[8]]$BEST.GUESS

lane9$mitoPercent<-matching_rows_all[[9]]$mitoPercent
lane9$Droplet<-matching_rows_all[[9]]$DROPLET.TYPE
lane9$Patient<-matching_rows_all[[9]]$BEST.GUESS

lane10$mitoPercent<-matching_rows_all[[10]]$mitoPercent
lane10$Droplet<-matching_rows_all[[10]]$DROPLET.TYPE
lane10$Patient<-matching_rows_all[[10]]$BEST.GUESS

lane11$mitoPercent<-matching_rows_all[[11]]$mitoPercent
lane11$Droplet<-matching_rows_all[[11]]$DROPLET.TYPE
lane11$Patient<-matching_rows_all[[11]]$BEST.GUESS

lane12$mitoPercent<-matching_rows_all[[12]]$mitoPercent
lane12$Droplet<-matching_rows_all[[12]]$DROPLET.TYPE
lane12$Patient<-matching_rows_all[[12]]$BEST.GUESS

lane13$mitoPercent<-matching_rows_all[[13]]$mitoPercent
lane13$Droplet<-matching_rows_all[[13]]$DROPLET.TYPE
lane13$Patient<-matching_rows_all[[13]]$BEST.GUESS

lane14$mitoPercent<-matching_rows_all[[14]]$mitoPercent
lane14$Droplet<-matching_rows_all[[14]]$DROPLET.TYPE
lane14$Patient<-matching_rows_all[[14]]$BEST.GUESS

lane15$mitoPercent<-matching_rows_all[[15]]$mitoPercent
lane15$Droplet<-matching_rows_all[[15]]$DROPLET.TYPE
lane15$Patient<-matching_rows_all[[15]]$BEST.GUESS

lane16$mitoPercent<-matching_rows_all[[16]]$mitoPercent
lane16$Droplet<-matching_rows_all[[16]]$DROPLET.TYPE
lane16$Patient<-matching_rows_all[[16]]$BEST.GUESS

lane17$mitoPercent<-matching_rows_all[[17]]$mitoPercent
lane17$Droplet<-matching_rows_all[[17]]$DROPLET.TYPE
lane17$Patient<-matching_rows_all[[17]]$BEST.GUESS

lane18$mitoPercent<-matching_rows_all[[18]]$mitoPercent
lane18$Droplet<-matching_rows_all[[18]]$DROPLET.TYPE
lane18$Patient<-matching_rows_all[[18]]$BEST.GUESS

lane19$mitoPercent<-matching_rows_all[[19]]$mitoPercent
lane19$Droplet<-matching_rows_all[[19]]$DROPLET.TYPE
lane19$Patient<-matching_rows_all[[19]]$BEST.GUESS

lane20$mitoPercent<-matching_rows_all[[20]]$mitoPercent
lane20$Droplet<-matching_rows_all[[20]]$DROPLET.TYPE
lane20$Patient<-matching_rows_all[[20]]$BEST.GUESS

lane21$mitoPercent<-matching_rows_all[[21]]$mitoPercent
lane21$Droplet<-matching_rows_all[[21]]$DROPLET.TYPE
lane21$Patient<-matching_rows_all[[21]]$BEST.GUESS

lane22$mitoPercent<-matching_rows_all[[22]]$mitoPercent
lane22$Droplet<-matching_rows_all[[22]]$DROPLET.TYPE
lane22$Patient<-matching_rows_all[[22]]$BEST.GUESS

lane23$mitoPercent<-matching_rows_all[[23]]$mitoPercent
lane23$Droplet<-matching_rows_all[[23]]$DROPLET.TYPE
lane23$Patient<-matching_rows_all[[23]]$BEST.GUESS

lane24$mitoPercent<-matching_rows_all[[24]]$mitoPercent
lane24$Droplet<-matching_rows_all[[24]]$DROPLET.TYPE
lane24$Patient<-matching_rows_all[[24]]$BEST.GUESS

#Merge lanes
all.lanes<-merge(lane1, c(lane10, lane11, lane12, lane13, lane14,lane15, lane16, lane17,
                          lane18, lane19, lane2, lane20, lane21, lane22, lane23, lane24,lane3, lane4,
                          lane5, lane6, lane7, lane8, lane9), add.cell.ids=ls()[43:66])

#Filter SINGLETS
all.lanes<-subset(all.lanes, subset=Droplet=="SNG")

# Filter matching datasets with BEST.GUESS (assigning cells to samples)
CHI014_tr <- subset(all.lanes, subset = Patient == "CHI014_B3,CHI014_B3,0.00")
#SHD1_NA,SHD1_NA,0.00 706
SHD1_tr <- subset(all.lanes, subset=Patient == "SHD1_NA,SHD1_NA,0.00")
#SHD3_NA,SHD3_NA,0.00  785
SHD3_tr <- subset(all.lanes, subset=Patient == "SHD3_NA,SHD3_NA,0.00")
#SHD5_NA,SHD5_NA,0.00  615            
SHD5_tr <- subset(all.lanes, subset=Patient == "SHD5_NA,SHD5_NA,0.00")
#SHD6_NA,SHD6_NA,0.00630
SHD6_tr <- subset(all.lanes, subset=Patient == "SHD6_NA,SHD6_NA,0.00")

#Merge samples (donors) 
merge.tr<-merge(CHI014_tr, c(SHD1_tr, SHD3_tr, SHD5_tr, SHD6_tr))
#77,908 cells in total 

#save
saveRDS(merge.tr, "../../data/merged-tr-raw.rds")

#clean mitoPercent
merge.tr<-subset(merge.tr, subset = mitoPercent >= 2 & mitoPercent <= 15)
#57,103 cells in total after mito filter


#processing
merge.tr<-NormalizeData(merge.tr)
merge.tr<-FindVariableFeatures(merge.tr, selection.method = "vst", nfeatures = 2000)
merge.tr<-ScaleData(merge.tr)
merge.tr<-RunPCA(merge.tr)
ElbowPlot(merge.tr, ndims = 50)
merge.tr<-FindNeighbors(merge.tr, dims=1:25)
merge.tr<-FindClusters(merge.tr) 
merge.tr<-RunUMAP(merge.tr, dims=1:25, return.model = T)
DimPlot(merge.tr, label=T, reduction = "umap", group.by = "Patient")
View(merge.tr@meta.data)
merge.tr.sct<-SCTransform(merge.tr, vars.to.regress = "mitoPercent")



#integration CCA
int.tr<-IntegrateLayers(merge.tr, method = CCAIntegration, 
                        orig.reduction = "pca", new.reduction="integrated.cca.immuno",
                        verbose=FALSE)
int.tr<-FindNeighbors(int.tr, reduction = "integrated.cca.immuno", dims=1:25)
int.tr<-FindClusters(int.tr)
int.tr<-RunUMAP(int.tr, reduction = "integrated.cca.immuno", dims = 1:25, reduction.name = "umap.cca.immuno")
int.tr.har<-IntegrateLayers(merge.tr, method=HarmonyIntegration, orig.reduction = "pca")
int.tr.har<-FindNeighbors(int.tr.har, dims=1:25)
int.tr.har<-FindClusters(int.tr.har)
int.tr.har<-RunUMAP(int.tr.har, dims = 1:25)


#Joinlayers
int.tr.j<-JoinLayers(int.tr.har)
int.tr.j<-JoinLayers(int.tr)

#Azimuth
int.tr.ann<-RunAzimuth(int.tr.j, reference="pbmcref")
int.tr<-RunAzimuth(int.tr.j, reference="pbmcref")

#Plot UMAP in this case by Seurat clusters
DimPlot(int.tr.ann, label = T , label.size = 3, group.by = c("seurat_clusters","predicted.celltype.l2"))+NoLegend()
DimPlot(int.tr, label = T, label.size = 3, group.by = c("seurat_clusters","predicted.celltype.l2"))+NoLegend()

#Rename Idents
int.tr.ann<-RenameIdents(object=int.tr.ann,"0"="CD4 TCM","1"="CD4 TM","2"="CD8 TEM","3"="CD4 TCM",
                         "4"="CD8 TCM","5"="CD14 Mono","6"="CD14 Mono","7"="B memory","8"="NK",
                         "9"="B naive","10"="CD14 Mono","11"="Platelet","12"="CD14 Mono",
                         "13"="Platelet", "14"="B Intermediate","15"="CD8 Naive", "16"="CD16 Mono",
                         "17"="CD4 TCM 2", "18"="MAIT", "19"="CD14 Mono","20"="Proliferative","21"="Plasmablast","22"="CD14 Mono 2",
                         "23"="cDC","24"="Eryth")

#add idents
int.tr.ann$idents<-Idents(int.tr.ann)

#Dimplot
DimPlot(int.immuno.joined, reduction = "umap.cca.immuno", label = T, pt.size = 2, label.size = 5, group.by = "idents")


#Cell type counting
conteo.tr<-table(int.tr.ann$Patient, int.tr.ann$predicted.celltype.l2)
tr.df<-as.data.frame(conteo.tr)

conteo.tr.2<-table(int.tr$Patient, int.tr$predicted.celltype.l2)
tr.df.2<-as.data.frame(conteo.tr.2)

fwrite(tr.df.2, "res/dataframes-t-count/timeres-df-2.txt")

saveRDS(int.tr.ann, "res/int-tr-ann.rds")