
####SiteCELL2
LC2<-Read10X("../Data/SiteCELL/LC2/")
LC2<-CreateSeuratObject(counts = LC2,min.cells = 3, min.features = 200)
View(LC2@meta.data) ##17,571 

#Add column of mitochondrial rna
LC2$mitoPercent<- PercentageFeatureSet(LC2, pattern = "^MT-")

#Add BARCODE column to the metadata
LC2$BARCODE<-rownames(LC2@meta.data)

#Export metadata to a file
metadata_LC2<-LC2@meta.data

#Import demux best-file
demux1a<-read.table(file="../Data/SiteCELL/best-files/LC2.best", header = T)
View(demux1a) ### 76,538

demux1a.2 <- demux1a %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(demux1a.2$Droplet)
#AMB   DBL   SNG 
#52729  1643 22166

demux1a.2<-subset(demux1a.2, subset = Droplet == "SNG")
table(demux1a.2$Patient)
#1112_LCMX0004 25_LC-CO-0025-ES  3_LC-CO-0003-ES      95_LCMX0147 
#5476             4675             3873             8142

#Check for matching barcodes
matching_rows<-merge(metadata_LC2, demux1a, by="BARCODE",all = TRUE)
View(matching_rows) #### 76,544

matching_rows<-subset(matching_rows, subset=orig.ident=="SeuratProject")
View(matching_rows) ##17,571


matching_rows.2 <- matching_rows %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(matching_rows.2$Droplet)
#AMB   DBL   SNG 
#455  1057 16053

matching_rows.2<-subset(matching_rows.2, subset=Droplet=="SNG")
table(matching_rows.2$Patient)
#112_LCMX0004 25_LC-CO-0025-ES  3_LC-CO-0003-ES      95_LCMX0147 
#3992             3332             2567             6162

freem1a<-read.table(file="../Data/SiteCELL/freemux-files/LC2.clust1.samples.gz", header = T)
View(freem1a) ###17,590

matching_barcodes<-merge(matching_rows.2, freem1a, by="BARCODE", all=TRUE)
View(matching_barcodes)

matching_barcodes<-subset(matching_barcodes, subset=Droplet=="SNG")
View(matching_barcodes) # 16,053 
table(matching_barcodes$Patient)

#112_LCMX0004 25_LC-CO-0025-ES  3_LC-CO-0003-ES      95_LCMX0147 
#3992             3332             2567             6162

MX04<-subset(matching_barcodes, subset = Patient=="112_LCMX0004")
table(MX04$BEST.GUESS) #========================= 1
#0,0  1,0  1,1  2,0  2,1  2,2  3,0  3,1  3,2  3,3 
#9    9 3955    2    2    2    2    5    1    5

CO25<-subset(matching_barcodes, subset = Patient=="25_LC-CO-0025-ES")
table(CO25$BEST.GUESS) #================================3
#0,0  1,0  1,1  2,0  2,1  2,2  3,0  3,1  3,2  3,3 
#6    5    1    3    1    6   11    8    7 3284

CO3<-subset(matching_barcodes, subset = Patient=="3_LC-CO-0003-ES")
table(CO3$BEST.GUESS)#=================================2
#0,0  1,0  1,1  2,0  2,1  2,2  3,0  3,2  3,3 
#5    4    4   15   10 2509    2   11    7

MX147<-subset(matching_barcodes, subset = Patient=="95_LCMX0147")
table(MX147$BEST.GUESS) #=================================0
#0,0  1,0  1,1  2,0  2,1  2,2  3,0  3,1  3,3 
#6092   18    4   19    1    6   12    1    9

#Merge metadata and freemux outputs by "BARCODE" column
matching1a.freemxmetadata<-merge(freem1a, metadata_LC2, by="BARCODE", all=TRUE)
matching1a.freemxmetadata<-subset(matching1a.freemxmetadata, subset=orig.ident =="SeuratProject")

#Merge previously merged data and demux4b by "BARCODE" column
matching.freemxmetaxdemux<-merge(matching1a.freemxmetadata, demux1a, by="BARCODE", all=TRUE)
matching.freemxmetaxdemux<-subset(matching.freemxmetaxdemux, subset=orig.ident=="SeuratProject") ##17,571

#Merge only pool metadata and freemux file, keep only SNG and them assign genotypes
matching.freemxmetaxdemux<-matching.freemxmetaxdemux %>%
  mutate(Patient=NA) %>%
  mutate(Patient =ifelse(BEST.GUESS == "0,0","95_LCMX0147", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "1,1","112_LCMX0004", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "2,2","3_LC-CO-0003-ES", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "3,3","25_LC-CO-0025-ES", Patient))

table(matching.freemxmetaxdemux$Patient)
#112_LCMX0004 25_LC-CO-0025-ES  3_LC-CO-0003-ES      95_LCMX0147 
#4132             3469             2722             6341

#How many SNG do we have FROM FREEMUXLET?
table(matching.freemxmetaxdemux$DROPLET.TYPE)
#AMB   DBL   SNG 
#16   907 16648

#Add columns to the metadata of the Seurat object
LC2@meta.data$Droplet<-matching.freemxmetaxdemux$DROPLET.TYPE
LC2@meta.data$Patient<-matching.freemxmetaxdemux$Patient
LC2@meta.data$demuxSNP<-matching.freemxmetaxdemux$N.SNP
LC2@meta.data$freemuxSNP<-matching.freemxmetaxdemux$NUM.SNPS

View(LC2@meta.data) ###17,571

#Delete "NAs"
LC2@meta.data<-LC2@meta.data %>% drop_na() ### 16,658

table(LC2@meta.data$Droplet)
#AMB SNG
#14  16644

#Subsetting
filt.LC2<-subset(LC2, subset = Droplet == 'SNG')

View(filt.LC2@meta.data)
#After subset 16,071

#Keep subsetting
LCCO0003<-subset(filt.LC2, subset = Patient == '3_LC-CO-0003-ES')
LCMX0147<-subset(filt.LC2, subset = Patient == '95_LCMX0147')
LCCO0025<-subset(filt.LC2, subset = Patient == '25_LC-CO-0025-ES')
LCMX0004<-subset(filt.LC2, subset = Patient == '112_LCMX0004')

#Merge
merged.1a<-merge(x=LCCO0003, y=c(LCCO0025,LCMX0004,LCMX0147), 
                 add.cell.ids=ls()[92:95])

#save
saveRDS(merged.1a, "res/merged-1a-raw.rds")


#FILTER MITORNA
merged.1a<-subset(merged.1a, subset = mitoPercent >=2 & mitoPercent <= 15)

merged.1a<-NormalizeData(merged.1a)
merged.1a<-FindVariableFeatures(merged.1a, selection.method = "vst", nfeatures = 3000)
merged.1a<-ScaleData(merged.1a)
merged.1a<-RunPCA(merged.1a)
#ElbowPlot(merged.1a, ndims = 50)
merged.1a<-FindNeighbors(merged.1a, dims=1:35)
merged.1a<-FindClusters(merged.1a) 
merged.1a<-RunUMAP(merged.1a, dims=1:35, return.model = T)

#integration
int.LC2<-IntegrateLayers(merged.1a, method = CCAIntegration, 
                            orig.reduction = "pca", new.reduction="integrated.cca",
                            verbose=FALSE)

int.LC2<-FindNeighbors(int.LC2, reduction = "integrated.cca", dims=1:35)
int.LC2<-FindClusters(int.LC2)

int.LC2<-RunUMAP(int.LC2, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")

DimPlot(int.LC2, reduction = "umap.cca", group.by = c("seurat_clusters","Patient"), label = TRUE)
DimPlot(int.LC2, reduction = "umap.cca", group.by = c("celltype", "seurat_clusters"), label = TRUE)

#Percentage of cells
table(int.LC2@meta.data$celltype)

#JoinLayers
int.LC2.j<-JoinLayers(int.LC2)

#Azimuth
int.LC2.ann<-RunAzimuth(int.LC2.j, reference="pbmcref")

int.LC2.ann<-RenameIdents(object=int.LC2.ann,"0"="NK","1"="B Naive","2"="CD4 TCM","3"="CD14 Mono",
                             "4"="CD8 TEM","5"="CD8 TE","6"="CD4 Naive","7"="CD4 TE",
                             "8"="CD4 TCM 2","9"="B Memory","10"="NK 2","11"="NK T effector",
                             "12"="B Intermediate", "13"="Treg","14"="CD16 Mono","15"="CD14 Mono 2",
                             "16"="CD8 Naive", "17"="cDC", "18"="NK Proliferative","19"="pDC")

#Add idents
int.LC2.ann$idents<-Idents(int.LC2.ann)

#cell type counting
counting.p1a<-table(int.LC2.ann$Patient, int.LC2.ann$predicted.celltype.l2)
p1a.df<-as.data.frame(counting.p1a)
