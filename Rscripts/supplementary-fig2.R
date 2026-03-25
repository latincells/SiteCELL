##########################Supplementary Figure 2 ###############################
#Brazil experiment using same samples with both ficoll and LatinCells protocols
#Latincells = SiteCELL

#Import table
brazil<-read.table("~/citometry_PBMCs.txt", header = TRUE)

#Delete "%" sign 
brazil$Percentage<-as.numeric(gsub("%", "", brazil$Percentage))

#Meand and SD calculation
brazil2<-brazil%>%
  group_by(Protocol, Celltype)%>%
  summarise(Mean_perc=mean(Percentage, na.rm=TRUE),
            SD_perc=sd(Percentage, na.rm=TRUE))



#Label sorting
celltype_order <- c("Leukocyte/CD45+",
                    "TcellCD4/CD45+/CD3+/CD4+", 
                    "TcellCD8/CD45+/CD3+/CD8a+",
                    "Tcell/CD45+/CD3+",
                    "Bcell/CD45+/CD19+",
                    "Neutrophil/CD16+CD66b+",
                    "MonocyteIntermediate/CD14++CD16++",
                    "MonocyteClassicalCD14++CD16-",
                    "MonocyteNonClassical/CD14+CD16++",
                    "NK/CD14-CD16+")

brazil2$Celltype <- factor(brazil2$Celltype, levels = celltype_order)

#filter
brazil2<-brazil2 %>%
  filter(Celltype!="Leukocyte/CD45+",
         Celltype!="Tcell/CD45+/CD3+")

brazil2$Celltype_clean <- recode(brazil2$Celltype, 
                                 "TcellCD4/CD45+/CD3+/CD4+" = "T cell CD4+",
                                 "TcellCD8/CD45+/CD3+/CD8a+" = "T cell CD8+", 
                                 "Bcell/CD45+/CD19+" = "B cell", 
                                 "Neutrophil/CD16+CD66b+" = "Neutrophil", 
                                 "MonocyteIntermediate/CD14++CD16++" = "Intermediate Monocyte", 
                                 "MonocyteClassicalCD14++CD16-" = "Classical Monocyte", 
                                 "MonocyteNonClassical/CD14+CD16++" = "Non Classical Monocyte", 
                                 "NK/CD14–CD16+" = "NK cell" ) 


# Plot
propsbrazil<-ggplot(brazil2, aes(x = fct_inorder(Celltype_clean), y = Mean_perc, fill = Protocol)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean_perc, ymax = Mean_perc + SD_perc),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "", y = "Frequency (%CD45+)", fill = "") +
  scale_fill_manual(values=colores.stress)+
  scale_y_continuous(limits = c(0, 50)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################### UMAPS ################

#Read integrated object
int.all.objects<-readRDS("~/data/int-meta-object.rds")

#Ficoll subsetting
Ficoll<-subset(int.all.objects, subset = Isolation_protocol == "Ficoll" & predicted.celltype.l2.score >= 0.5)

#LatinCells subsetting
LatCells<-subset(int.all.objects, subset = Isolation_protocol == "LatinCells" & predicted.celltype.l2.score >= 0.5)

#Reclustering
Ficoll<-FindNeighbors(Ficoll, reduction = "integrated.cca", dims=1:35)
Ficoll<-FindClusters(Ficoll, resolution = 1)
Ficoll<-RunUMAP(Ficoll, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")

#Reclustering
LatCells<-FindNeighbors(LatCells, reduction = "integrated.cca", dims=1:35)
LatCells<-FindClusters(LatCells, resolution = 1)
LatCells<-RunUMAP(LatCells, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")


#RENAME IDENTS
LatCells<-RenameIdents(object=LatCells,"0"="CD4 TEM","1"="CD8 TEM","2"="CD14 Mono","3"="CD4 TCM",
                       "4"="B Naive","5"="CD8 TCM","6"="CD4 Naive","7"="CD4 CTL","8"="B Memory",
                       "9"="NK","10"="CD14 Mono","11"="CD16 Mono","12"="gdT",
                       "13"="dnT", "14"="CD4 TCM","15"="CD8 Naive", "16"="MAIT",
                       "17"="Treg", "18"="CD14 Mono", "19"="NK CD56 Bright/ILC","20"="B Intermediate","21"="cDC","22"="CD14 Mono",
                       "23"="T/NK Proliferative","24"="B Naive","25"="pDC/ASDC","26"="Plasmablast","27"="HSPC","28"="Platelet")

#adding idents
LatCells$idents<-Idents(LatCells)


#RENAME IDENTS
Ficoll<-RenameIdents(object=Ficoll,"0"="CD4 TCM","1"="CD8 TCM","2"="CD4 CTL","3"="CD14 Mono",
                     "4"="CD4 Naive","5"="NK","6"="CD4 TEM","7"="CD8 TEM","8"="B Naive",
                     "9"="CD16 Mono","10"="CD14 Mono","11"="CD8 Naive","12"="gdT",
                     "13"="CD14 Mono", "14"="B Intermediate","15"="Treg", "16"="MAIT",
                     "17"="NK CD56 Bright/ILC", "18"="CD4 TCM", "19"="cDC","20"="Platelet","21"="B Memory","22"="CD14 Mono",
                     "23"="CD8 TCM","24"="T/NK Proliferative","25"="pDC/ASDC","26"="Plasmablast","27"="CD4 TCM","28"="Eryth",
                     "29"="CD14 Mono","30"="CD14 Mono","31"="CD14 Mono","32"="HSPC","33"="dnT","34"="B Naive")

#adding idents
Ficoll$idents<-Idents(Ficoll)

#Ficoll
umap2<-DimPlot(Ficoll, reduction = "umap.cca", group.by = "idents", raster = FALSE, alpha = 0.7)
umap3<-DimPlot(Ficoll, reduction = "umap.cca", group.by = "Dataset",raster = FALSE, alpha = 0.7)

ggsave(filename = "~/umapficoll-idents.jpg", plot = umap2, device = "jpg", width = 12, height = 7, dpi=300)
ggsave(filename = "~/umapficoll-datasets.jpg", plot = umap3, device = "jpg", width = 12, height = 7, dpi=300)

#LatinCells
umap4<-DimPlot(LatCells, reduction = "umap.cca", group.by = "idents", raster = FALSE, alpha = 0.7)
umap5<-DimPlot(LatCells, reduction = "umap.cca", group.by = "Dataset",raster = FALSE, alpha = 0.7)

umap4+umap5

##########################Supplementary Figure 2 ###############################
#Brazil experiment using same samples with both ficoll and LatinCells protocols

#Import table
brazil<-read.table("~/citometry_PBMCs.txt", header = TRUE)

#Delete "%" sign 
brazil$Percentage<-as.numeric(gsub("%", "", brazil$Percentage))

#Meand and SD calculation
brazil2<-brazil%>%
  group_by(Protocol, Celltype)%>%
  summarise(Mean_perc=mean(Percentage, na.rm=TRUE),
            SD_perc=sd(Percentage, na.rm=TRUE))



#Label sorting
celltype_order <- c("Leukocyte/CD45+",
                    "TcellCD4/CD45+/CD3+/CD4+", 
                    "TcellCD8/CD45+/CD3+/CD8a+",
                    "Tcell/CD45+/CD3+",
                    "Bcell/CD45+/CD19+",
                    "Neutrophil/CD16+CD66b+",
                    "MonocyteIntermediate/CD14++CD16++",
                    "MonocyteClassicalCD14++CD16-",
                    "MonocyteNonClassical/CD14+CD16++",
                    "NK/CD14-CD16+")

brazil2$Celltype <- factor(brazil2$Celltype, levels = celltype_order)

#filter
brazil2<-brazil2 %>%
  filter(Celltype!="Leukocyte/CD45+",
         Celltype!="Tcell/CD45+/CD3+")

brazil2$Celltype_clean <- recode(brazil2$Celltype, 
                                 "TcellCD4/CD45+/CD3+/CD4+" = "T cell CD4+",
                                 "TcellCD8/CD45+/CD3+/CD8a+" = "T cell CD8+", 
                                 "Bcell/CD45+/CD19+" = "B cell", 
                                 "Neutrophil/CD16+CD66b+" = "Neutrophil", 
                                 "MonocyteIntermediate/CD14++CD16++" = "Intermediate Monocyte", 
                                 "MonocyteClassicalCD14++CD16-" = "Classical Monocyte", 
                                 "MonocyteNonClassical/CD14+CD16++" = "Non Classical Monocyte", 
                                 "NK/CD14–CD16+" = "NK cell" ) 


# Plot
propsbrazil<-ggplot(brazil2, aes(x = fct_inorder(Celltype_clean), y = Mean_perc, fill = Protocol)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean_perc, ymax = Mean_perc + SD_perc),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "", y = "Frequency (%CD45+)", fill = "") +
  scale_fill_manual(values=colores.stress)+
  scale_y_continuous(limits = c(0, 50)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################### UMAPS ################

#Read integrated object
int.all.objects<-readRDS("~/data/int-meta-object.rds")

#Ficoll subsetting
Ficoll<-subset(int.all.objects, subset = Isolation_protocol == "Ficoll" & predicted.celltype.l2.score >= 0.5)

#LatinCells subsetting
LatCells<-subset(int.all.objects, subset = Isolation_protocol == "LatinCells" & predicted.celltype.l2.score >= 0.5)

#Reclustering
Ficoll<-FindNeighbors(Ficoll, reduction = "integrated.cca", dims=1:35)
Ficoll<-FindClusters(Ficoll, resolution = 1)
Ficoll<-RunUMAP(Ficoll, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")

#Reclustering
LatCells<-FindNeighbors(LatCells, reduction = "integrated.cca", dims=1:35)
LatCells<-FindClusters(LatCells, resolution = 1)
LatCells<-RunUMAP(LatCells, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")


#RENAME IDENTS
LatCells<-RenameIdents(object=LatCells,"0"="CD4 TEM","1"="CD8 TEM","2"="CD14 Mono","3"="CD4 TCM",
                       "4"="B Naive","5"="CD8 TCM","6"="CD4 Naive","7"="CD4 CTL","8"="B Memory",
                       "9"="NK","10"="CD14 Mono","11"="CD16 Mono","12"="gdT",
                       "13"="dnT", "14"="CD4 TCM","15"="CD8 Naive", "16"="MAIT",
                       "17"="Treg", "18"="CD14 Mono", "19"="NK CD56 Bright/ILC","20"="B Intermediate","21"="cDC","22"="CD14 Mono",
                       "23"="T/NK Proliferative","24"="B Naive","25"="pDC/ASDC","26"="Plasmablast","27"="HSPC","28"="Platelet")

#adding idents
LatCells$idents<-Idents(LatCells)


#RENAME IDENTS
Ficoll<-RenameIdents(object=Ficoll,"0"="CD4 TCM","1"="CD8 TCM","2"="CD4 CTL","3"="CD14 Mono",
                     "4"="CD4 Naive","5"="NK","6"="CD4 TEM","7"="CD8 TEM","8"="B Naive",
                     "9"="CD16 Mono","10"="CD14 Mono","11"="CD8 Naive","12"="gdT",
                     "13"="CD14 Mono", "14"="B Intermediate","15"="Treg", "16"="MAIT",
                     "17"="NK CD56 Bright/ILC", "18"="CD4 TCM", "19"="cDC","20"="Platelet","21"="B Memory","22"="CD14 Mono",
                     "23"="CD8 TCM","24"="T/NK Proliferative","25"="pDC/ASDC","26"="Plasmablast","27"="CD4 TCM","28"="Eryth",
                     "29"="CD14 Mono","30"="CD14 Mono","31"="CD14 Mono","32"="HSPC","33"="dnT","34"="B Naive")

#adding idents
Ficoll$idents<-Idents(Ficoll)

#Ficoll
umap2<-DimPlot(Ficoll, reduction = "umap.cca", group.by = "idents", raster = FALSE, alpha = 0.7)
umap3<-DimPlot(Ficoll, reduction = "umap.cca", group.by = "Dataset",raster = FALSE, alpha = 0.7)

ggsave(filename = "~/umapficoll-idents.jpg", plot = umap2, device = "jpg", width = 12, height = 7, dpi=300)
ggsave(filename = "~/umapficoll-datasets.jpg", plot = umap3, device = "jpg", width = 12, height = 7, dpi=300)

#LatinCells
umap4<-DimPlot(LatCells, reduction = "umap.cca", group.by = "idents", raster = FALSE, alpha = 0.7)
umap5<-DimPlot(LatCells, reduction = "umap.cca", group.by = "Dataset",raster = FALSE, alpha = 0.7)

umap4+umap5

ggsave(filename = "~/umaplatcell-idents.jpg", plot = umap4, device = "jpg", width = 12, height = 7, dpi=300)
ggsave(filename = "~/umaplatcell-dataset.jpg", plot = umap5, device = "jpg", width = 12, height = 7, dpi=300)



