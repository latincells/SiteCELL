################################ LatinCells1 

#Importing data
LC1<-Read10X("../Data/LatinCells/PBMC_1n2_Otomi/")
LC1<-CreateSeuratObject(counts = LC1, min.cells = 3, min.features = 200)

#Add column of mitochondrial rna
LC1$mitoPercent<- PercentageFeatureSet(LC1, pattern = "^MT-")

#Add BARCODE column to the metadata
LC1$BARCODE<-rownames(LC1@meta.data)

#Export metadata to a file
metadata_LC1<-LC1@meta.data
View(metadata_LC1) #14144 barcodes

#Import demux best-file
demuxp1n2<-read.table(file="../Data/LatinCells/best-files/LC1.best", header = T)
View(demuxp1n2) ###  71,461

demuxp1n2.2 <- demuxp1n2 %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(demuxp1n2.2$Droplet)
#AMB   DBL   SNG 
#8783   272 10318 

demuxp1n2.2<-subset(demuxp1n2.2, subset = Droplet == "SNG")
table(demuxp1n2.2$Patient)
#90_LCMX0142 91_LCMX0143 92_LCMX0144 
#4387        5881          50

#Check for matching barcodes
matching_rows<-merge(metadata_LC1, demuxp1n2, by="BARCODE",all = TRUE)
View(matching_rows) ####  71,494 

matching_rows<-subset(matching_rows, subset=orig.ident=="SeuratProject")
View(matching_rows) ##17,067


matching_rows.2 <- matching_rows %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(matching_rows.2$Droplet)
#AMB  DBL  SNG 
#3993  269 9749

matching_rows.2<-subset(matching_rows.2, subset=Droplet=="SNG")
table(matching_rows.2$Patient)
#90_LCMX0142 91_LCMX0143 92_LCMX0144 
4012        5691          46

freemp1n2<-read.table(file="../Data/LatinCells/freemux-files/LC1.clust1.samples.gz", header = T)
View(freemp1n2) ###17075

matching_barcodes<-merge(matching_rows.2, freemp1n2, by="BARCODE", all=TRUE)
View(matching_barcodes)

matching_barcodes<-subset(matching_barcodes, subset=Droplet=="SNG")
View(matching_barcodes) # 12,390

MX142<-subset(matching_barcodes, subset = Patient=="90_LCMX0142")
table(MX142$BEST.GUESS) #========================= 0
#0,0  1,0  1,1 
#3864  118   30

MX143<-subset(matching_barcodes, subset = Patient=="91_LCMX0143")
table(MX143$BEST.GUESS) #================================1
#0,0  1,0  1,1 
#55   70 5566

MX144<-subset(matching_barcodes, subset = Patient=="92_LCMX0144")
table(MX144$BEST.GUESS) #================================NO
#0,0 1,0 1,1 
#26   4  16

#Merge metadata and freemux outputs by "BARCODE" column
matching1b.freemxmetadata<-merge(freemp1n2, metadata_LC1, by="BARCODE", all=TRUE)
matching1b.freemxmetadata<-subset(matching1b.freemxmetadata, subset=orig.ident =="SeuratProject")

#Merge previously merged data and demux4b by "BARCODE" column
matching.freemxmetaxdemux<-merge(matching1b.freemxmetadata, demuxp1n2, by="BARCODE", all=TRUE)
matching.freemxmetaxdemux<-subset(matching.freemxmetaxdemux, subset=orig.ident=="SeuratProject")

#Merge only pool metadata and freemux file, keep only SNG and then assign genotypes
matching.freemxmetaxdemux<-matching.freemxmetaxdemux %>%
  mutate(Patient=NA) %>%
  mutate(Patient =ifelse(BEST.GUESS == "0,0","90_LCMX0142", Patient))%>%
  mutate(Patient =ifelse(BEST.GUESS == "1,1","91_LCMX0143", Patient))

table(matching.freemxmetaxdemux$Patient)
#90_LCMX0142 91_LCMX0143 
#5246        8493

#How many SNG do we have FROM FREEMUXLET?
table(matching.freemxmetaxdemux$DROPLET.TYPE)
#AMB   DBL   SNG 
#5   405 13734

View(LC1@meta.data)  #17,067

#Add columns to the metadata of the Seurat object
LC1@meta.data$Droplet<-matching.freemxmetaxdemux$DROPLET.TYPE
LC1@meta.data$Patient<-matching.freemxmetaxdemux$Patient
LC1@meta.data$demuxSNP<-matching.freemxmetaxdemux$N.SNP
LC1@meta.data$freemuxSNP<-matching.freemxmetaxdemux$NUM.SNPS

View(LC1@meta.data) ###14144

#Delete "NAs"
LC1@meta.data<-LC1@meta.data %>% drop_na() ### 13,611

#SUBSET BY PATIENT  
p1n2MX0142<-subset(LC1, subset = Patient == '90_LCMX0142')
p1n2MX0143<-subset(LC1, subset = Patient == '91_LCMX0143')

#Import individual 3 LC1
p1n23MX0144<-Read10X("../Data/LatinCells/PBMC_3_Otomi/")
p1n23MX0144<-CreateSeuratObject(counts = p1n23MX0144, min.cells = 3, min.features = 200)

#Add column of mitochondrial rna
p1n23MX0144$mitoPercent<- PercentageFeatureSet(p1n23MX0144, pattern = "^MT-")

#Add BARCODE column to the metadata
p1n23MX0144$BARCODE<-rownames(p1n23MX0144@meta.data)

vect<-rep('SNG', 21378 )
p1n23MX0144@meta.data<-p1n23MX0144@meta.data %>% mutate(Droplet=vect)
vect2<-rep('92_LCMX0144',21378 )
p1n23MX0144@meta.data<-p1n23MX0144@meta.data %>% mutate(Sample = vect2)

View(p1n23MX0144@meta.data)
#21,378 cells


#Starting with merge
merged.LC1<-merge(p1n23MX0144, y=c(p1n2MX0142, p1n2MX0143),
                    add.cell.ids=ls()[138:140])


saveRDS(merged.LC1, "merged-LC1-raw.rds")
#34989 cells total

View(merged.LC1@meta.data)
saveRDS(merged.LC1, file = "merged-LC1.rds")


#Create new Patient column
merged.LC1$Sample<-rownames(merged.LC1@meta.data)
merged.LC1@meta.data<-separate(merged.LC1@meta.data, col='Sample', into = c('Patient', 'Barcode'), 
                                 sep='_')

merged.LC1<-subset(merged.LC1, subset = mitoPercent >= 2 & mitoPercent <= 15)
#After merged 34,023


#Object merged
merged.LC1<-NormalizeData(merged.LC1)
merged.LC1<-FindVariableFeatures(merged.LC1, selection.method = "vst", nfeatures = 3000)
merged.LC1<-ScaleData(merged.LC1)
merged.LC1<-RunPCA(merged.LC1)
#ElbowPlot(merged.LC1, ndims = 50)
merged.LC1<-FindNeighbors(merged.LC1, dims=1:30)
merged.LC1<-FindClusters(merged.LC1) 
merged.LC1<-RunUMAP(merged.LC1, dims=1:30, return.model = T)

int.LC1<-IntegrateLayers(merged.LC1, method = CCAIntegration, 
                           orig.reduction = "pca", new.reduction="integrated.cca",
                           verbose=FALSE)

int.LC1<-FindNeighbors(int.LC1, reduction = "integrated.cca", dims=1:35)
int.LC1<-FindClusters(int.LC1)
int.LC1<-RunUMAP(int.LC1, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")
DimPlot(int.LC1, reduction = "umap.cca", group.by = c("Patient","seurat_clusters"))

#join layers 
int.LC1.j<-JoinLayers(int.LC1)

#Azimuth
int.LC1.ann<-RunAzimuth(int.LC1.j, reference ="pbmcref")

int.LC1.ann<-RenameIdents(object=int.LC1.ann,"0"="CD14 Mono","1"="CD8 Effector 2","2"="CD8 TM","3"="CD8 TCM",
                            "4"="NK 2","5"="CD4 TCM","6"="CD8 TEM","7"="B Memory",
                            "8"="B Naive","9"="NK 1","10"="CD16 Mono","11"="CD4 Naive",
                            "12"="dgT", "13"="NK T","14"="MAIT","15"="CD8 Effector",
                            "16"="Treg", "17"="CD14 Mono 2", "18"="CD4","19"="cDC1","20"="CD8 Naive 2",
                            "21"="CD14 Mono 1","22"="B Naive","23"="Proliferative","24"="cDC2","25"="Plasmablast")

int.LC1.ann$idents<-Idents(int.LC1.ann)

#
saveRDS(int.LC1.ann, "res/int-LC1-ann.rds")

#Cell type counting
conteo.LC1<-table(int.LC1.ann$Patient, int.LC1.ann$predicted.celltype.l2)
LC1.df<-as.data.frame(conteo.LC1)

fwrite(LC1.df, "res/dataframes-t-count/LC1-df.txt")

