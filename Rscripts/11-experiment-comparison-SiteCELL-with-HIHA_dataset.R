
###############################--------------MX0049 Matched sample experiment figure ------###########

fwrite(demux2c, "demux2c.txt")
fwrite(freemux2c, "freemux2c.txt")


fwrite(freemux4b, "pool4b.clust1.samples.gz")
fwrite(demux4b, "demux4b.best")

#LIBPOOL2B
LP2C<-Read10X(data.dir = "~/filtered_feature_bc_matrix/")
LP2C<-CreateSeuratObject(LP2C, project = "LP2C", min.cells = 3, min.features = 200)

LP4B<-Read10X(data.dir = "~/filtered_feature_bc_matrix/")
LP4B<-CreateSeuratObject(LP4B, project = "LP4B", min.cells = 3, min.features = 200)

#Add column of mitochondrial rna
LP2C$mitoPercent<- PercentageFeatureSet(LP2C, pattern = "^MT-")
LP4B$mitoPercent<- PercentageFeatureSet(LP4B, pattern = "^MT-")

#Add BARCODE column to the metadata
LP2C$BARCODE<-rownames(LP2C@meta.data)
LP4B$BARCODE<-rownames(LP4B@meta.data)

#Export metadata to a file
metadata_LP2C<-LP2C@meta.data
metadata_LP4B<-LP4B@meta.data

#Import demux best-file
demux2c<-read.table(file="demux2c.best", header = T, sep = ",")
demux4b<-read.table(file="demux4b.best", header = T, sep = ",")
View(demux2c) ### 38,000
View(demux4b) ### 23,000

demux2c.2 <- demux2c %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

demux4b.2 <- demux4b %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(demux2c.2$Droplet)
#  AMB   DBL   SNG 
# 25669 666 11820 

table(demux4b.2$Droplet)
#  AMB   DBL   SNG 
# 15989 551 7258 

demux2c.2<-subset(demux2c.2, subset = Droplet == "SNG")
table(demux2c.2$Patient)
# 11_LC_MX_0080    27_LC-MX-0049    29_LC_MX_0101 62_LC-CO-0061-ES 
# 3951             3936             2378             1555

demux4b.2<-subset(demux4b.2, subset = Droplet == "SNG")
table(demux4b.2$Patient)
# 27_LC-MX-0049 29_LC-CO-0029-ES  6_LC-CO-0006-ES 61_LC-CO-0060-ES      97_LCMX0149 
# 1635             1230             2129             1085             1179

##Change barcodes into BARCODE
demux2c.2 <- demux2c.2 %>%
  rename(BARCODE = barcodes)

demux4b.2 <- demux4b.2 %>%
  rename(BARCODE = barcodes)

matching_rows2c.2 <- matching_rows2c.2 %>%
  rename(BARCODE = barcodes)

matching_rows4b.2 <- matching_rows4b.2 %>%
  rename(BARCODE = barcodes)

#Check for matching barcodes
matching_rows2c<-merge(metadata_LP2C, demux2c.2, by="BARCODE",all = TRUE)
View(matching_rows2c) #### 18

matching_rows4b<-merge(metadata_LP4B, demux4b.2, by="BARCODE",all = TRUE)
View(matching_rows4b) #### 76,544

matching_rows2c<-subset(matching_rows2c, subset=orig.ident=="LP2C")
View(matching_rows2c) ##17,266
matching_rows4b<-subset(matching_rows4b, subset=orig.ident=="LP4B")
View(matching_rows4b) ##16,978


matching_rows2c <- matching_rows2c %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

matching_rows4b <- matching_rows4b %>%
  separate(BEST, into = c("Droplet", "Patient"), sep = "-", extra = "merge")

table(matching_rows2c$Droplet)
#SNG 
#10390

table(matching_rows4b$Droplet)
#SNG 
#6872

matching_rows2c.2<-subset(matching_rows2c, subset=Droplet=="SNG")
table(matching_rows2c.2$Patient)
# 11_LC_MX_0080    27_LC-MX-0049    29_LC_MX_0101 62_LC-CO-0061-ES 
# 3367             3356             2225             1442

matching_rows4b.2<-subset(matching_rows4b, subset=Droplet=="SNG")
table(matching_rows4b.2$Patient)
#  27_LC-MX-0049 29_LC-CO-0029-ES  6_LC-CO-0006-ES 61_LC-CO-0060-ES      97_LCMX0149 
#  1542             1190             2002             1030             1108 

######Fremuxlet ######
freem2c<-read.table(file="pool2c.clust1.samples.gz", header = T, sep = ",")
View(freem2c) ###17,293 

freem4b<-read.table(file="pool4b.clust1.samples.gz", header = T, sep = ",")
View(freem4b) ###17,007 

#change barcodes for BARCODE in freemuxlet files
freem2c <- freem2c %>%
  rename(BARCODE = barcodes)

freem4b <- freem4b %>%
  rename(BARCODE = barcodes)

matching_barcodes2c<-merge(matching_rows2c.2, freem2c, by="BARCODE", all=TRUE)
View(matching_barcodes2c) #17,293
matching_barcodes4b<-merge(matching_rows4b.2, freem4b, by="BARCODE", all=TRUE)
View(matching_barcodes4b) #17,007


matching_barcodes2c.2<-subset(matching_barcodes2c, subset=Droplet=="SNG")
View(matching_barcodes2c.2) # 10,390 
table(matching_barcodes2c.2$Patient)

#11_LC_MX_0080    27_LC-MX-0049    29_LC_MX_0101 62_LC-CO-0061-ES 
#3367             3356             2225             1442

matching_barcodes4b.2<-subset(matching_barcodes4b, subset=Droplet=="SNG")
View(matching_barcodes4b.2) # 6,872
table(matching_barcodes4b.2$Patient)

# 27_LC-MX-0049 29_LC-CO-0029-ES  6_LC-CO-0006-ES 61_LC-CO-0060-ES      97_LCMX0149 
# 1542             1190             2002             1030             1108

MX49.2c<-subset(matching_barcodes2c.2, subset = Patient=="27_LC-MX-0049")
table(MX49.2c$BEST.GUESS) #========================= 1
# 0,0  1,0  1,1  2,0  2,1  2,2  3,0  3,1  3,3 
# 3197   34   28   38    3   18   16    3   19

######Pool2C
#Merge metadata and freemux outputs by "BARCODE" column
matching2c.freemxmetadata<-merge(freem2c, metadata_LP2C, by="BARCODE", all=TRUE)
matching2c.freemxmetadata<-subset(matching2c.freemxmetadata, subset=orig.ident =="LP2C") #17225

#Merge previously merged data and demux4b by "BARCODE" column
matching2c.freemxmetaxdemux<-merge(matching2c.freemxmetadata, demux2c.2, by="BARCODE", all=TRUE) #total lines 18675
matching2c.freemxmetaxdemux<-subset(matching2c.freemxmetaxdemux, subset=orig.ident=="LP2C") # 17245

#Merge only pool metadata and freemux file, keep only SNG and them assign genotypes
matching2c.freemxmetaxdemux<-matching2c.freemxmetaxdemux %>%
  mutate(Patient=NA) %>%
  mutate(Patient =ifelse(BEST.GUESS == "0,0","27_LC-MX-0049", Patient))

table(matching2c.freemxmetaxdemux$Patient)
#27_LC-MX-0049 
#5468

#How many SNG do we have FROM FREEMUXLET?
table(matching2c.freemxmetaxdemux$DROPLET.TYPE)
# AMB   DBL   SNG 
# 10   752 16504 

#Add columns to the metadata of the Seurat object
LP2C@meta.data$Droplet<-matching2c.freemxmetaxdemux$DROPLET.TYPE
LP2C@meta.data$Patient<-matching2c.freemxmetaxdemux$Patient
LP2C@meta.data$demuxSNP<-matching2c.freemxmetaxdemux$N.SNP
LP2C@meta.data$freemuxSNP<-matching2c.freemxmetaxdemux$NUM.SNPS

View(LP2C@meta.data) ###17,266 

#Delete "NAs"
LP2C@meta.data<-LP2C@meta.data %>% drop_na() ### 16,658

table(LP2C@meta.data$Droplet)
#AMB SNG
#3327

#Subsetting
filt.LP2C<-subset(LP2C, subset = Droplet == 'SNG')

######################## ---------------- Pool4B---------------######################
MX49.4b<-subset(matching_barcodes4b.2, subset = Patient=="27_LC-MX-0049")
table(MX49.4b$BEST.GUESS) #================================3
# 0,0  1,0  1,1  2,0  2,1  2,2  3,1  3,2  4,0  4,1  4,2  4,3  4,4 
# 7    1    3   28   31 1394    1   19    3    6   39    1    9

#Merge metadata and freemux outputs by "BARCODE" column
matching4b.freemxmetadata<-merge(freem4b, metadata_LP4B, by="BARCODE", all=TRUE)
matching4b.freemxmetadata<-subset(matching4b.freemxmetadata, subset=orig.ident =="LP4B") #16937

#Merge previously merged data and demux4b by "BARCODE" column
matching4b.freemxmetaxdemux<-merge(matching4b.freemxmetadata, demux4b.2, by="BARCODE", all=TRUE) #total lines 17343
matching4b.freemxmetaxdemux<-subset(matching4b.freemxmetaxdemux, subset=orig.ident=="LP4B") # 16957

#Merge only pool metadata and freemux file, keep only SNG and them assign genotypes
matching4b.freemxmetaxdemux<-matching4b.freemxmetaxdemux %>%
  mutate(Patient=NA) %>%
  mutate(Patient =ifelse(BEST.GUESS == "2,2","27_LC-MX-0049", Patient))

table(matching4b.freemxmetaxdemux$Patient)
#27_LC-MX-0049 
#2990

#How many SNG do we have FROM FREEMUXLET?
table(matching4b.freemxmetaxdemux$DROPLET.TYPE)
# AMB   DBL   SNG 
# 20  1277 15681

#Add columns to the metadata of the Seurat object
LP4B@meta.data$Droplet<-matching4b.freemxmetaxdemux$DROPLET.TYPE
LP4B@meta.data$Patient<-matching4b.freemxmetaxdemux$Patient
LP4B@meta.data$demuxSNP<-matching4b.freemxmetaxdemux$N.SNP
LP4B@meta.data$freemuxSNP<-matching4b.freemxmetaxdemux$NUM.SNPS

#Delete "NAs"
LP4B@meta.data<-LP4B@meta.data %>% drop_na() ### 

table(LP4B@meta.data$Droplet)

#Subsetting
filt.LP4B<-subset(LP4B, subset = Droplet == 'SNG')

#########---------Preprocessing------------#######

lp2c<-readRDS("LP2C_LC-MX-0049.rds")
lp4b<-readRDS("LP4B_LC-MX-0049.rds")

lp2c@meta.data <- lp2c@meta.data %>%
  mutate(Pool = "LP2C")

lp4b@meta.data <- lp4b@meta.data %>%
  mutate(Pool = "LP4B")

#reset for merge
rna_counts_lp2c <- GetAssayData(lp2c, layer = "counts")
rna_counts_lp4b <- GetAssayData(lp4b, layer = "counts")

lp2c_clean <- CreateSeuratObject(
  counts = rna_counts_lp2c,
  meta.data = lp2c@meta.data)

lp4b_clean <- CreateSeuratObject(
  counts = rna_counts_lp4b,
  meta.data = lp4b@meta.data)

validObject(lp2c_clean)
validObject(lp4b_clean)

#merge
merged.MX49 <- merge(
  lp2c_clean,
  y = lp4b_clean,
  add.cell.ids = c("LP2C", "LP4B"))

##########--------------mito, nCount, nFeature -----------------############
metadata.mergedMX49<-merged.MX49@meta.data

#filtrado NA
metadata.mergedMX49 <- metadata.mergedMX49 %>%
  filter(!is.na(Pool), Pool != "")

############---------------------mitoPercent------------------------############

pool.colors <- c(
  "LP2C" = "#C6B7E2", 
  "LP4B" = "#F4B183"  
)

mergedmx49.mtplot<-ggplot(metadata.mergedMX49, aes(x = Pool, y = mitoPercent, 
                                         fill = Pool, color = Pool)) + 
  geom_violin(width = 0.7, color = "black", alpha = 0.9) +
  stat_summary(aes(group = Pool), fun = "mean", geom = "point", shape = 20,
               size = 1, fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = pool.colors) +
  scale_color_manual(values = pool.colors) +
  theme_classic(base_size = 14) +
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, by = 25)
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),  
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"),  
    axis.text.x = element_text(size = 6, face = "bold", family = "Arial"),  
    axis.text.y = element_text(size = 8, family = "Arial"),             
    legend.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Dataset",
    y = "mtRNA (%)"
  )

ggsave(filename = "mergedmx49-mtplot.jpg", plot = mergedmx49.mtplot, device = "jpg", width = 12, height = 9.5, dpi=600)

#########---------------------nCount plot-------------#######

mergedmx49.counts<-ggplot(metadata.mergedMX49, aes(x = Pool, y = nCount_RNA,  fill = Pool, color = Pool)) + 
  geom_violin(width = 0.7, color = "black", alpha = 0.9) +
  stat_summary(aes(group = Pool), fun = "mean", geom = "point", shape = 20,
               size = 1, fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = pool.colors) +
  scale_color_manual(values = pool.colors) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"), 
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"),  
    axis.text.x = element_text(size = 6, face = "bold", family = "Arial"),  
    axis.text.y = element_text(size = 8, family = "Arial"),               
    legend.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Dataset",
    y = "UMIs per cell"
  )

ggsave(filename = "mergedmx49-counts.jpg", plot = mergedmx49.counts, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "mergedmx49-counts.svg", plot = mergedmx49.counts, device = "svg", width = 12, height = 9.5, dpi=600)


#---------------------------------------------nFeatures

mergedmx49.features<-ggplot(metadata.mergedMX49, aes(x = Pool, y = nFeature_RNA, 
                                                   fill = Pool, color = Pool)) + 
  #geom_jitter(position = position_jitter(width = 0.25), size = 0.8, alpha = 0.8) +
  #geom_violin(width = 0.3, size = 0.7, col = "black", alpha = 0.9) +
  geom_violin(width = 0.7, color = "black", alpha = 0.9) +
  stat_summary(aes(group = Pool), fun = "mean", geom = "point", shape = 20,
               size = 1, fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = pool.colors) +
  scale_color_manual(values = pool.colors) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),  
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"), 
    axis.text.x = element_text(size = 6, face = "bold", family = "Arial"),   
    axis.text.y = element_text(size = 8, family = "Arial"),                 
    legend.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(x = "Dataset",  y = "Genes per cell")

ggsave(filename = "mergedmx49-features.jpg", plot = mergedmx49.features, device = "jpg", width = 12, height = 9.5, dpi=600)

##################-----------VIABILITY----------######
viab <- data.frame(
  Pool = factor(c("LP2C", "LP4B")),
  Viabilidad = c(96.6, 96.6))

viab.mx49.both<-ggplot(viab, aes(x = Pool, y = Viability, fill = Pool)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = pool.colors) +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 96.6, linetype = "dashed", color = "lightgray", linewidth = 0.7) +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Pool",
    y = "Viability (%)"
  ) +
  ylim(0, 100)


ggsave(filename = "viab-mx49-both.jpg", plot = viab.mx49.both, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "viab-mx49-both.svg", plot = viab.mx49.both, device = "svg", width = 12, height = 9.5, dpi=600)

#################-----------------processing merge---------------#####################

# Filter NA cells in metadata
merged.MX49 <- subset(merged.MX49, subset = !is.na(mitoPercent))
merged.MX49 <- subset(merged.MX49, subset = !is.na(nCount_RNA))
merged.MX49 <- subset(merged.MX49, subset = !is.na(nFeature_RNA))
merged.MX49 <- subset(merged.MX49, subset = !is.na(Pool))
merged.MX49 <- subset(merged.MX49, subset = !is.na(BARCODE))
merged.MX49 <- subset(merged.MX49, subset = !is.na(orig.ident))

#Processing
merged.MX49<-NormalizeData(merged.MX49)
merged.MX49<-FindVariableFeatures(merged.MX49, selection.method = "vst", nfeatures = 3000)
merged.MX49<-ScaleData(merged.MX49)
merged.MX49<-RunPCA(merged.MX49)
merged.MX49<-FindNeighbors(merged.MX49, dims=1:35)
merged.MX49<-FindClusters(merged.MX49, resolution = 1) 
merged.MX49<-RunUMAP(merged.MX49, dims=1:35, return.model = T)

DimPlot(merged.MX49, label=T, reduction = "umap", group.by = c("Pool","seurat_clusters"))


int.mx49<-IntegrateLayers(merged.MX49, method = CCAIntegration, 
                             orig.reduction = "pca", new.reduction="integrated.cca",
                             verbose=FALSE)

int.mx49<-FindNeighbors(int.mx49, reduction = "integrated.cca", dims=1:30)
int.mx49<-FindClusters(int.mx49, resolution = 1)

int.mx49<-RunUMAP(int.mx49, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

int.mx49.j<-JoinLayers(int.mx49)

int.mx49.ann<-RunAzimuth(int.mx49.j, reference="pbmcref")

int.mx49.ann@meta.data<-int.mx49.ann@meta.data%>%
  mutate(MajorGroup = case_when(
    predicted.celltype.l2 %in% c("CD8 TEM", "CD4 Naive",
                                 "CD8 TCM","CD4 TCM",
                                 "CD4 TEM","CD8 Naive","gdT","dnT","ILC","Treg","MAIT","CD4 CTL","CD4 Proliferative","CD8 Proliferating","NK Proliferating") ~  "T",
    predicted.celltype.l2 %in% c("NK","NK_CD56bright","NK Proliferating")~  "NK",
    predicted.celltype.l2 %in% c("B memory","B intermediate","B naive","Plasmablast")~  "B",
    predicted.celltype.l2 %in% c("CD14 Mono","CD16 Mono")~ "Monocyte",
    predicted.celltype.l2 %in% c("cDC1","cDC2","pDC","ASDC")~ "DC",
    predicted.celltype.l2 %in% c("Eryth")~ "Eryth",
    predicted.celltype.l2 %in% c("Platelet")~ "Platelet",
    predicted.celltype.l2 %in% c("HSPC")~ "HSPC",
    TRUE ~ "NO CLASS"
  ))

# 
int.mx49.ann.filt <- subset(int.mx49.ann, subset = MajorGroup != "NO CLASS")


umap.both.mx49<-DimPlot(int.mx49.ann.filt, reduction = "umap.cca", group.by = "MajorGroup", split.by="Pool",alpha = 0.5, pt.size = 1.2, label=FALSE)

ggsave(filename = "umap-mx49-both.jpg", plot = umap.both.mx49, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "umap-mx49-both.svg", plot = umap.both.mx49, device = "svg", width = 12, height = 9.5, dpi=600)


#######################-----------------COUNTS---------------####################
repmeta<-int.mx49.ann.filt@meta.data

prop.repmeta<-repmeta%>%
  add_count(Patient, Pool, predicted.celltype.l2, name = "subtype_count")

prop.repmeta<-prop.repmeta%>%
  add_count(Patient, Pool, MajorGroup, name = "MG_count")

proporciones_subtype <- prop.repmeta %>%
  group_by(Patient, Pool) %>%
  mutate(
    Total_in_Group = n(),
    n_subtype = subtype_count, 
    proportion_subtype = (n_subtype / Total_in_Group)*100,
    prop_mg_subtype=(MG_count/Total_in_Group)*100
  ) %>%
  ungroup() %>%
  select(Patient, Pool, predicted.celltype.l2, MajorGroup, Total_in_Group, n_subtype, proportion_subtype,prop_mg_subtype) %>%
  distinct()


df_agg <- proporciones_subtype %>%
  group_by(Pool, MajorGroup) %>%
  summarise(
    Percentage = mean(prop_mg_subtype)
  ) %>%
  ungroup()

#filter HSPC
df_agg1<-df_agg%>%
  filter(MajorGroup!="HSPC")

# Cell proportion
prop.mg.mx49.both<-ggplot(df_agg1, aes(x = MajorGroup, y = Percentage, fill = Pool)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", 
           size = 0.5) +
  scale_y_continuous(limits = c(0, max(df_agg$Percentage) * 1.1), 
                     expand = expansion(mult = 0)) + 
  scale_fill_manual(values = pool.colors) +
  scale_color_manual(values = pool.colors) +
  labs(
    x = "Major Group",
    y = "Cell percentage (%)")+
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )+
  theme_classic(base_size = 14)


ggsave(filename = "majorgroup-prop-mx49-both.jpg", plot = prop.mg.mx49.both, device = "jpg", width = 12, height = 9.5, dpi=600)

########################---------------Expression Matrix------------------##################
int.mx49.ann.filt@meta.data<-int.mx49.ann.filt@meta.data%>%
  mutate(Donor=case_when(
    Pool == "LP2C" ~ "MX491",
    Pool == "LP4B" ~ "MX492"))


###AGGREGATE EXPRESSION
mx49.exp<-AggregateExpression(int.mx49.ann.filt, return.seurat = T,
                               assays = "RNA",
                               group.by = c( "Pool","Donor","MajorGroup"))

#Getting matrix
pb.mx49<-as.matrix(GetAssayData(mx49.exp, layer = 'data')[, WhichCells(mx49.exp)])

#Write objects
fwrite(pb.mx49, "mx49-aggregatecounts.txt", row.names = TRUE, col.names = TRUE)


#Read object
aggregatecounts<-read.csv("mx49-aggregatecounts.txt", header = TRUE)

#extracting vector of azimuth 
azimuthgenes<-fread("exp_matrix_filt_brasil.txt")

#Unique genes found in azimuth as marker genes
interest.genes <- unique(azimuthgenes$azimuth.features)

##Filter
colnames_split <- strsplit(colnames(aggregatecounts)[-1], 
suffixes <- sapply(colnames_split, `[`, 2)
# Create a named list of columns by suffix
col_groups <- split(colnames(aggregatecounts)[-1], suffixes)

# Filter the data frame for the genes of interest
df_filt.markers <- aggregatecounts %>% filter(X %in% interest.genes)

# Initialize a data frame to store the results
df_stress.markers <- data.frame(interest.genes = df_filt.markers$X)


#expression scatter plots
# Pearson correlation
cor_res.exp <- cor.test(df_stress.markers$MX491,
                    df_stress.markers$MX492,
                    method = "pearson")

r_val.exp <- cor_res.exp$estimate
p_val.exp <- cor_res.exp$p.value

p <- ggplot(df_stress.markers, aes(x = MX491, y = MX492)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_fixed() +
  theme_classic() +
  annotate("text", x = Inf, y = -Inf,
    hjust = 1.1, vjust = -0.5,
    label = paste0("r = ", round(r_val.exp, 2), "\np = ", signif(p_val.exp, 2)))

ggsave(filename = "genexpcorr-mx49.jpg", plot = p, device = "jpg", width = 12, height = 9.5, dpi=600)

####### ---------- Cell type scatter plot--------- #########
prop_df <- prop.repmeta %>%
  distinct(Pool, MajorGroup, MG_count) %>%
  group_by(Pool) %>%
  mutate(
    total_cells = sum(MG_count),
    proportion = (MG_count / total_cells)*100
  ) %>%
  ungroup()

#Pivot wide
prop_wide <- prop_df %>%
  select(Pool, MajorGroup, proportion) %>%
  pivot_wider(names_from = Pool, values_from = proportion)

#conversion to 
pseudo <- 1e-4

prop_wide <- prop_wide %>%
  mutate(
    LP2C_log = LP2C + pseudo,
    LP4B_log = LP4B + pseudo)

###Correlation celltype
corr <- cor.test(prop_wide$LP2C_log, prop_wide$LP4B_log, method = "pearson")
r_val <- round(corr$estimate, 2)
p_val <- signif(corr$p.value, 3)

ctprop.mx49<-ggplot(prop_wide,
       aes(x = LP2C_log, y = LP4B_log, color = MajorGroup)) +
  geom_point(size = 9) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "lightgray") +
  coord_fixed() +
  theme_classic() +
  annotate("text",
           x = min(prop_wide$LP2C_log),
           y = max(prop_wide$LP4B_log),
           label = paste0("r = ", r_val, "\n", "p = ", p_val),
           hjust = 0, vjust = 1,
           size = 4)

ggsave(filename = "celltypecorr-mx49.jpg", plot = ctprop.mx49, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "celltypecorr-mx49.svg", plot = ctprop.mx49, device = "svg", width = 12, height = 9.5, dpi=600)

##################------HEATMAP-------###############

#From META OBJECT matrix
aggregatecounts<-read.csv("mx49-aggregatecounts.txt", header = TRUE, row.names = 1)

#Crossing feature.immune against matrix
matriz_filtrada.mx49 <- aggregatecounts[rownames(aggregatecounts) %in% features.stress, ]


log_mat.mx49<-apply(log1p(matriz_filtrada.mx49), 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
})



#matrix wrangling
matriz_filtrada.mx492 <- as.data.frame(matriz_filtrada.mx49)%>%
  rownames_to_column("Gene")%>%
  pivot_longer(cols=-Gene, names_to = "Pool", values_to = "Expression")%>%
  filter(Gene %in% features.stress)%>%
  mutate(MG = sub("_LC*","",Pool),
         Treatment=case_when(
           grepl("LP2C", Pool)~ "LP2C",
           grepl("LP4B",Pool)~"LP4B"
         ))

#filter names
df_clean <- matriz_filtrada.mx492 %>%
  mutate(
    MG = sub("^LP(2C|4B)_MX49[12]_", "", MG))


#filter HSPC
df_clean <- df_clean %>%
  filter(MG != "HSPC")

#wrangling
df_clean <- df_clean %>%
  mutate(Sample_CellType = paste(Treatment, MG, sep = "_"))

#Set color palette
palette_colors<-c("#762A83FF","white", "#FFA500") 

#Heatmap Plot
heatmap.mx49<-ggplot(df_clean, aes(x = Sample_CellType, y = Gene, fill = Expression)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  scale_fill_gradientn(
    colors=palette_colors)+
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12), 
    axis.text.y = element_text(face = "italic", size = 12), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)), 
    legend.position = "right", 
    legend.title = element_text(face = "bold", size = 12),  
    panel.grid = element_blank()
  ) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x))  

ggsave(filename = "heatmap-mx49.jpg", plot = heatmap.mx49, device = "jpg",width = 12, height = 8, dpi = 300)



##################---------------Comparison SiteCELL subcelltype proportions with annotated dataset---------------##################
##### Extra, not in the paper
Gong Q, Sharma M, Kuan EL, Glass MC, Chander A, Singh M, et al. Longitudinal Multi-omic Immune Profiling 
Reveals Age-Related Immune Cell Dynamics in Healthy Adults. Nature. 2025. doi:10.1038/s41586-025-09686-5

#Zellkonverter + SingleCellExperiment + Seurat
hiha<-readH5AD("HIHA/human_immune_health_atlas_full.h5ad")
hiha<-as.Seurat(hiha, counts = "counts", data="logcounts")
#convert assay a seurat
assayNames(hiha)[assayNames(hiha) == "X"] <- "counts"

seu <- CreateSeuratObject(counts = SummarizedExperiment::assay(hiha, "counts"))
seu <- AddMetaData(seu, metadata = as.data.frame(SummarizedExperiment::colData(hiha)))

seu.AIFI3<-subset(seu, subset=AIFI_L3)

#subsetting using healthy individuals (negative for ctyomegalovirus) and only T cells
seu.tcell<-subset(seu, subset = AIFI_L1=="T cell" & subject.cmv == "Negative")
#Number of lines = 651,818 T cells

########## LOAD HIHA database
seu <- LoadH5Seurat("HIHA/human_immune_health_atlas_full.h5seurat")

##############By Database
#Create metadata file
metadata<-seu.tcell@meta.data

table(metadata$AIFI_L3)

vec <- c("CD4 MAIT" = 1616, "CD8 MAIT" = 29007, "CD8aa" = 4079, "CM CD4 T cell" = 94097, "CM CD8 T cell" = 16761, "Core naive CD4 T cell" = 211784, "Core naive CD8 T cell" = 74157, "DN T cell" = 1477, "GZMB+ Vd2 gdT" = 11790, "GZMB- CD27+ EM CD4 T cell" = 38241, "GZMB- CD27- EM CD4 T cell" = 40133, "GZMK+ CD27+ EM CD8 T cell" = 29587, "GZMK+ Vd2 gdT" = 13247, "GZMK+ memory CD4 Treg" = 229, "GZMK- CD27+ EM CD8 T cell" = 3175, "ISG+ MAIT" = 515, "ISG+ memory CD4 T cell" = 2558, "ISG+ memory CD8 T cell" = 648,"ISG+ naive CD4 T cell" = 4118,"ISG+ naive CD8 T cell" = 459,"KLRB1+ memory CD4 Treg" = 1533,"KLRB1+ memory CD8 Treg" = 344,"KLRF1+ GZMB+ CD27- EM CD8 T cell" = 7303, "KLRF1+ effector Vd1 gdT" = 566, "KLRF1- GZMB+ CD27- EM CD8 T cell" = 13005, "KLRF1- GZMB+ CD27- memory CD4 T cell" = 748, "KLRF1- effector Vd1 gdT" = 155,  "Memory CD4 Treg" = 10235,  "Memory CD8 Treg" = 407, "Naive CD4 Treg" = 10252, "Naive Vd1 gdT" = 3050, "Proliferating T cell" = 1365, "SOX4+ Vd1 gdT" = 1219, "SOX4+ naive CD4 T cell" = 20206,  "SOX4+ naive CD8 T cell" = 3752)

table <- data.frame(celltype = names(vec), counts = as.numeric(vec))%>%
  mutate(proporcion = (counts / sum(counts))*100)

table <- table %>%
  mutate(
    celltype2 = case_when(
      grepl("CM CD8 T cell", celltype) ~ "CD8 TCM",
      celltype == "CD8aa" ~ "CD8aa",
      celltype == "DN T cell" ~ "DN T cell",
      celltype %in% c("CD8 MAIT", "CD4 MAIT", "ISG+ MAIT") ~ "MAIT",
      celltype == "CM CD4 T cell" ~ "Memory CD4 T cell",
      celltype == "GZMB- CD27+ EM CD4 T cell" ~ "Memory CD4 T cell",
      celltype == "GZMB- CD27- EM CD4 T cell" ~ "Memory CD4 T cell",
      celltype == "ISG+ memory CD4 T cell" ~ "Memory CD4 T cell",
      celltype == "KLRF1- GZMB+ CD27- memory CD4 T cell" ~ "Memory CD4 T cell",
      celltype %in% c("Core naive CD4 T cell",
                      "ISG+ naive CD4 T cell",
                      "SOX4+ naive CD4 T cell") ~ "Naive CD4 T cell",
      celltype %in% c("ISG+ naive CD8 T cell",
                      "Core naive CD8 T cell",
                      "SOX4+ naive CD8 T cell") ~ "Naive CD8 T cell",
      celltype == "Proliferating T cell" ~ "Proliferating T cell",
      celltype %in% c("GZMK+ memory CD4 Treg",
                      "KLRB1+ memory CD4 Treg",
                      "KLRB1+ memory CD8 Treg",
                      "Naive CD4 Treg",
                      "Memory CD4 Treg",
                      "Memory CD8 Treg") ~ "Treg",
      celltype %in% c("GZMK+ Vd2 gdT",
                      "Naive Vd1 gdT",
                      "SOX4+ Vd1 gdT",
                      "KLRF1- effector Vd1 gdT",
                      "KLRF1+ effector Vd1 gdT",
                      "GZMB+ Vd2 gdT") ~ "gdt",
      celltype %in% c("ISG+ memory CD8 T cell",
                      "KLRF1+ GZMB+ CD27- EM CD8 T cell",
                      "GZMK- CD27+ EM CD8 T cell",
                      "GZMK+ CD27+ EM CD8 T cell",
                      "KLRF1- GZMB+ CD27- EM CD8 T cell") ~ "Effector Memory CD8 T cell"))

table2<-table%>%
  group_by(celltype2)%>%
  summarise(counts2=sum(counts), .groups = "drop")%>%
  mutate(proporcion2=(counts2/sum(counts2)*100))

CD8TCMplot<-ggplot(table2, aes(x = reorder(celltype2, proporcion2), 
                               y = proporcion2)) +
  geom_bar(stat = "identity", fill = "#a6bddb") +
  geom_hline(yintercept = 5,  linetype = "dotted", color = "lightgray",  linewidth = 0.7)+
  theme_classic() +
  labs(x = "Cell type", 
       y = "Proportion (%)")

ggsave(filename = "~/revission-CD8TCM-plot.jpg", plot = CD8TCMplot, device = "jpg", width = 12, height = 7, dpi=300)

HumImmHealthA<-table2

############By batch #142157 cells
seu.t.sub<-subset(seu.tcell, subset=batch_id %in% c("B039", "B050","B063","B082","B094","B138","B145"))

batchmeta<-seu.t.sub@meta.data

prop.batchmeta<-batchmeta%>%
  add_count(batch_id, donor_id, AIFI_L2, name = "subtype_count")

cell_col<-"AIFI_L2"

meta_clean <- batchmeta %>%
  select(batch_id, donor_id, celltype = all_of(cell_col))

meta_clean_df <- as.data.frame(meta_clean, stringsAsFactors = FALSE)
meta_clean_df <- meta_clean_df[, c("batch_id", "donor_id", "celltype")]
counts <- meta_clean_df %>%
  count(batch_id, donor_id, celltype, name = "n")
counts <- counts %>%
  group_by(batch_id, donor_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


#by cell type
ggplot(counts, aes(x = reorder(celltype, proportion),    y = proportion)) +
  geom_bar(stat = "identity", fill = "#a6bddb") +
  theme_classic() +
  labs(x = "Cell type", 
       y = "Proportion (%)")


#per donor
ggplot(counts, aes(x = batch_id, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw()

#xtra highlighting Memory CD8
ggplot(counts, aes(x = batch_id, y = proportion, fill = celltype == "Memory CD8 T cell")) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("TRUE" = "#a6bddb", "FALSE" = "lightgray"),
                    labels = c("FALSE" = "Other celltypes", "TRUE"  = "Memory CD8 T cell"),
                    name = "Celltype") +  theme_bw()


#donor x batch x cell type
counts_plot <- counts %>%
  mutate(color_group = ifelse(celltype == "Memory CD8 T cell", "Memory CD8 T cell", "Other"),
         donor_batch = paste0(donor_id, "\n(", batch_id, ")"))


ggplot(counts_plot,
       aes(x = donor_batch, y = proportion,
           fill = color_group)) +
  geom_bar(stat = "identity",
           position = "fill",
           color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Memory CD8 T cell" = "#a6bddb",
    "Other" = "lightgray")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))+
  theme_bw() +
  labs(x = "Donor ID (Batch)", 
       y = "Proportion",
       fill = "Celltype",
       title = "Proportion of Memory CD8 T cells per Donor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
