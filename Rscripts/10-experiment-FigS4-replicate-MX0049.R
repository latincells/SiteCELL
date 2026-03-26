install.packages("ggdist")
install.packages("ggridges")
install.packages("ggbeeswarm")

library(ggdist)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggbeeswarm)


######## Matching samples processed with FDG and SiteCELL
######################## MITOPERCENT FIGURE

merged.matched@meta.data<-merged.matched@meta.data%>%
  mutate(Dataset_ID="merged-B23")

merged.B24.seu@meta.data<-merged.B24.seu@meta.data%>%
  mutate(Dataset_ID="merged-B24")

merged.matched<-merge(x = merged.matched, y=merged.B24.seu)

saveRDS(merged.matched, "~/RDSs/brmatched-merged.rds")
mtRNA.brmatched<-readRDS("~/brmatched-data/RDSs/brmatched-merged.rds")

#Create a meta file from merged object
meta.mtRNA.brmatched<-mtRNA.brmatched@meta.data

#Collapse mitoPercent + percent.mito 
meta.mtRNA.brmatched$mtRNA.collapsed <- ifelse(is.na(meta.mtRNA.brmatched$mitoPercent), meta.mtRNA.brmatched$percent.mito, meta.mtRNA.brmatched$mitoPercent)

#Save all merged metadata (for mtRNA analyses)
fwrite(meta.mtRNA.brmatched, "~/brmatched-data/mtRNA-brmatched.txt")

#Data wrangling

#Add Patient column based on Dataset ID values
meta.mtRNA.brmatched<-meta.mtRNA.brmatched%>%
  mutate(Protocol_isolation = case_when(
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B23" ~ "SiteCELL",
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B24" ~ "Ficoll",
    TRUE ~ "Sin clasificación"  # Por si hay algún valor fuera de los especificados
  ))

colores.stress <- c('SiteCELL' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


meta.mtRNA.brmatched$Dataset_ID <- factor(
  meta.mtRNA.brmatched$Dataset_ID,
  levels = c("merged-B23", "merged-B24"),
  labels = c("SiteCELL", "Ficoll")
)

#Order datasets
meta.mtRNA.brmatched <- meta.mtRNA.brmatched %>%
  mutate(Dataset_ID = factor(Dataset_ID, levels = Dataset_ID[order(Protocol_isolation)]))


br.mtplot<-ggplot(meta.mtRNA.brmatched, aes(x = Protocol, y = mitoPercent, 
                                      fill = Protocol_isolation, color = Protocol_isolation)) + 
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  stat_summary(aes(group = Dataset_ID), fun = "mean", geom = "point", shape = 20,
               size = 1, fill = "white", color = "black", stroke = 1.2) +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) +
  theme_classic(base_size = 14) +
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, by = 25)
  ) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),  # Eje X título
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"),  # Eje Y título
    axis.text.x = element_text(size = 6, face = "bold", family = "Arial"),   # Etiquetas x
    axis.text.y = element_text(size = 8, family = "Arial"),                  # Etiquetas y
    legend.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Dataset",
    y = "mtRNA (%)"
  )

ggsave(filename = "figs/mtRNA.jpg", plot = br.mtplot, device = "jpg", width = 12, height = 9.5, dpi=600)


##################UMIS and GENES per 

colores.stress <- c('SiteCELL' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


# Ncount
meta.mtRNA.brmatched$Protocol <- ifelse(grepl("Ficoll", meta.mtRNA.brmatched$Group), "Ficoll", "SiteCELL")


meta.mtRNA.brmatched <- meta.mtRNA.brmatched %>%
  arrange(Protocol, Group) %>%
  mutate(Group = factor(Group, levels = unique(Group)))


# 
ncount.br<-ggplot(meta.mtRNA.brmatched, aes(x = Group, y = nCount_RNA, fill = Protocol_isolation)) + 
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  scale_fill_manual(values = c('SiteCELL' = '#1c9099', 'Ficoll' = '#a6bddb')) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +
  scale_x_discrete(labels = function(x) gsub("_.*", "", x)) + 
  scale_y_log10(breaks=c(1000,10000,25000,100000), labels=scales::comma)+
  theme_classic(base_size = 14) +
  labs(x = "Patient", y = "nCount_RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(filename = "figs/ncount-br.jpg", plot = ncount.br, device = "jpg", width = 12, height = 9.5, dpi=600)


#NFeatures
nfeatu.br<-ggplot(meta.mtRNA.brmatched, aes(x = Group, y = nFeature_RNA, fill = Protocol_isolation)) + 
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  scale_fill_manual(values = c('SiteCELL' = '#1c9099', 'Ficoll' = '#a6bddb')) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +
  scale_x_discrete(labels = function(x) gsub("_.*", "", x)) +
  scale_y_log10(breaks=c(1000,5000,10000), labels=scales::comma) +
  theme_classic(base_size = 14) +
  labs(x = "Patient", y = "nFeature_RNA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")


ggsave(filename = "figs/nfeature-br.jpg", plot = nfeatu.br, device = "jpg", width = 12, height = 9.5, dpi=600)

############################### STRESS GENE FIGURE

#Add Patient column based on Dataset ID values
merged.matched@meta.data<-merged.matched@meta.data%>%
  mutate(Protocol_isolation = case_when(
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B23" ~ "SiteCELL",
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B24" ~ "Ficoll",
    TRUE ~ "NONE"  
  ))

merged.matched<-NormalizeData(merged.matched)
merged.matched<-FindVariableFeatures(merged.matched, selection.method = "vst", nfeatures = 3000)
merged.matched<-ScaleData(merged.matched)
merged.matched<-RunPCA(merged.matched)
#ElbowPlot(merged.matched, ndims = 50)
merged.matched<-FindNeighbors(merged.matched, dims=1:35)
merged.matched<-FindClusters(merged.matched, resolution = 1) 
merged.matched<-RunUMAP(merged.matched, dims=1:35, return.model = T)

DimPlot(merged.matched, label=T, reduction = "umap", group.by = c("Protocol_isolation","seurat_clusters"))


int.mbrmatched<-IntegrateLayers(merged.matched, method = CCAIntegration, 
                         orig.reduction = "pca", new.reduction="integrated.cca",
                         verbose=FALSE)

int.mbrmatched<-FindNeighbors(int.mbrmatched, reduction = "integrated.cca", dims=1:30)
int.mbrmatched<-FindClusters(int.mbrmatched, resolution = 1)

int.mbrmatched<-RunUMAP(int.mbrmatched, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(int.mbrmatched.ann, reduction = "umap.cca", group.by=c("seurat_clusters"),alpha = 0.5, pt.size = 1.2, label=TRUE)

int.mbrmatched.j<-JoinLayers(int.mbrmatched)

int.mbrmatched.ann<-RunAzimuth(int.mbrmatched.j, reference="pbmcref")
DimPlot(int.mbrmatched.ann, reduction = "umap.cca", group.by=c("seurat_clusters", "predicted.celltype.l2"),alpha = 0.5, pt.size = 1.2, label=TRUE)

umaps.by.protocol<-DimPlot(int.mbrmatched.ann, reduction = "umap.cca",group.by = "predicted.celltype.l2", split.by="Protocol_isolation",alpha = 0.5, pt.size = 1.2, label=FALSE)
ggsave(filename = "~/brmatched-data/figs/umaps-both-protocol-br.jpg", plot = umaps.by.protocol, device = "jpg",width = 12, height = 8, dpi = 600)


brmatched.SiteCELL<-subset(int.mbrmatched.ann, subset = Protocol_isolation =="SiteCELL")
brmatched.SiteCELL<-FindNeighbors(brmatched.SiteCELL,reduction = "integrated.cca", dims = 1:30)
brmatched.SiteCELL<-FindClusters(brmatched.SiteCELL, resolution = 1)
brmatched.SiteCELL<-RunUMAP(brmatched.SiteCELL, reduction = "integrated.cca", dims=1:30, reduction.name="umap.cca")

brmatched.ficoll<-subset(int.mbrmatched.ann, subset = Protocol_isolation =="Ficoll")
brmatched.ficoll<-FindNeighbors(brmatched.ficoll,reduction = "integrated.cca", dims = 1:30)
brmatched.ficoll<-FindClusters(brmatched.ficoll, resolution = 1)
brmatched.ficoll<-RunUMAP(brmatched.ficoll, reduction = "integrated.cca", dims=1:30, reduction.name="umap.cca")


umap.br.lc<-DimPlot(brmatched.SiteCELL, reduction = "umap.cca", group.by="predicted.celltype.l2",alpha = 0.5, pt.size = 1.2, label=TRUE)+NoLegend()
ggsave(filename = "~/brmatched-data/figs/umap-br-lc.jpg", plot = umap.br.lc, device = "jpg",width = 12, height = 8, dpi = 600)

umap.br.fic<-DimPlot(brmatched.ficoll, reduction = "umap.cca", group.by="predicted.celltype.l2",alpha = 0.5, pt.size = 1.2, label=TRUE)
ggsave(filename = "~/brmatched-data/figs/umap-br-fic.jpg", plot = umap.br.fic, device = "jpg",width = 12, height = 8, dpi = 600)

#Adding extra Majorgroup celltype

int.mbrmatched.ann@meta.data<-int.mbrmatched.ann@meta.data%>%
  mutate(MajorGroup = case_when(
    predicted.celltype.l2 %in% c("CD8 TEM", "CD4 Naive",
                  "CD8 TCM","CD4 TCM",
                  "CD4 TEM","CD8 Naive","gdT","dnT/ILC","Treg","MAIT","CD4 CTL","T/NK Proliferative") ~  "T",
    predicted.celltype.l2 %in% c("NK","NK_CD56bright","T/NK Proliferative")~  "NK",
    predicted.celltype.l2 %in% c("B Memory","B Intermediate","B Naive","Plasmablast")~  "B",
    predicted.celltype.l2 %in% c("CD14 Mono","CD16 Mono")~ "Monocyte",
    predicted.celltype.l2 %in% c("cDC","pDC/ASDC")~ "DC",
    predicted.celltype.l2 %in% c("Eryth")~ "Eryth",
    predicted.celltype.l2 %in% c("Platelet")~ "Platelet",
    predicted.celltype.l2 %in% c("HSPC")~ "HSPC",
    TRUE ~ "Sin clasificación"  # Por si hay algún valor fuera de los especificados
  ),
  # Agrega sufijo a los pacientes según el Dataset_ID
  Patient = case_when(
    Dataset_ID == "merged-B23" ~ paste0(Patient, "-LC"),
    Dataset_ID == "merged-B24" ~ paste0(Patient, "-Ficoll"),
    TRUE ~ as.character(Patient)  # Por si hubiera algún otro valor
  ))

###AGGREGATE EXPRESSION
pb.brmatched<-AggregateExpression(int.mbrmatched.ann, return.seurat = T,
                                    assays = "RNA",
                                    group.by = c("predicted.celltype.l2","Patient", "Protocol_isolation"))

#Getting matrix
matrix.brmatched<-as.matrix(GetAssayData(pb.brmatched, layer = 'data')[, WhichCells(pb.brmatched)])

#Write objects
fwrite(matrix.brmatched, "~/brmatched-data/brmatched-aggregatecounts.txt", row.names = TRUE, col.names = TRUE)

#Read gene vestors
#Genes from Massoni-Badossa et al.
features.stress<-c("CIRBP", "H3F3B", "EIF1", "NEAT1", "FTH1", "SRGN", "RBM3", "SAT1", "CXCR4")

#Genes from Savage et al.
features.time <- c("H3F3B",  "NEAT1", "SRGN", "SERF2","CFL1","MYL12A",
                   "CD52","PFN1","DUSP1","FOS","JUNB","JUND",
                   "KLF6","SAT1","GIMAP7")

#Read object
aggregatecounts<-read.csv("~/brmatched-data/brmatched-aggregatecounts.txt", header = TRUE)

# Exclude the 'column1' column
colnames_split <- strsplit(colnames(aggregatecounts)[-1], "_")  
suffixes <- sapply(colnames_split, `[`, 2)
# Create a named list of columns by suffix
col_groups <- split(colnames(aggregatecounts)[-1], suffixes)
# Filter the data frame for the genes of interest
aggregatecounts$Eryth_LCBR0094.Ficoll_Ficoll<-NULL
aggregatecounts$Platelet_LCBR0093.Ficoll_Ficoll<-NULL
aggregatecounts$Platelet_LCBR0094.Ficoll_Ficoll<-NULL
aggregatecounts$Platelet_LCBR0095.Ficoll_Ficoll<-NULL
aggregatecounts$Platelet_LCBR0096.Ficoll_Ficoll<-NULL


df_filt.stress <- aggregatecounts %>%
  filter(X %in% features.stress)

df_filt.time<- aggregatecounts %>%
  filter(X %in% features.time)

# Initialize a data frame to store the results
df_stress <- data.frame(features.stress = df_filt.stress$X)
df_time <- data.frame(features.time = df_filt.time$X)

# Loop through each group and sum the columns
for (suffix in names(col_groups)) {
  cols <- col_groups[[suffix]]
  valid_cols <- cols[cols %in% colnames(df_filt.stress)]
  if (length(valid_cols) == 0) {
    warning(paste("Not found:", suffix))
    next
  }
  df_stress[[suffix]] <- rowMeans(df_filt.stress[valid_cols])
}


for (suffix in names(col_groups)) {
  cols <- col_groups[[suffix]]
  valid_cols <- cols[cols %in% colnames(df_filt.time)]
  if (length(valid_cols) == 0) {
    warning(paste("Not found:", suffix))
    next
  }
  df_time[[suffix]] <- rowMeans(df_filt.time[valid_cols])
}

df_long.stress <- df_stress %>%
  pivot_longer(-features.stress, names_to = "Sample", values_to = "Expression")

#Using "Protocol_isolation" column 
df_long.stress<-df_long.stress%>%
  mutate(Isolation = case_when(
    Sample %in% c("LCBR0093.Ficoll", "LCBR0094.Ficoll",
                  "LCBR0095.Ficoll","LCBR0096.Ficoll") ~  "Ficoll",
    Sample %in% c("LCBR0093.LC","LCBR0094.LC","LCBR0095.LC","LCBR0096.LC")~  "SiteCELL",
    TRUE ~ "NA"  
  ))


df_long.time <- df_time %>%
  pivot_longer(-features.time, names_to = "Sample", values_to = "Expression")

#Creating new "Dataset" column 
df_long.time<-df_long.time%>%
  mutate(Isolation = case_when(
    Sample %in% c("LCBR0093.Ficoll", "LCBR0094.Ficoll",
                  "LCBR0095.Ficoll","LCBR0096.Ficoll") ~  "Ficoll",
    Sample %in% c("LCBR0093.LC","LCBR0094.LC","LCBR0095.LC","LCBR0096.LC")~  "SiteCELL",
    TRUE ~ "NA"  
  ))

df_long.stress2.1<-df_long.stress%>%
  group_by(Sample,Isolation)%>%
  summarise(sumexp_stress=mean(Expression))

df_long.time2.1<-df_long.time%>%
  group_by(Sample,Isolation)%>%
  summarise(sumexp_stress=mean(Expression))


#reorder factor levels
df_long.stress2.1<-df_long.stress2.1%>%
  mutate(Isolation=reorder(Isolation, sumexp_stress, FUN=mean))

#reorder factor levels
df_long.time2.1<-df_long.time2.1%>%
  mutate(Isolation=reorder(Isolation, sumexp_stress, FUN=mean))


colores.stress <- c('SiteCELL' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


stressg<-ggplot(df_long.stress2.1, aes(x = Isolation, y = sumexp_stress, fill = Isolation, color=Isolation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.8) +  # Boxplot
  geom_point(shape = 21, fill = "black", size = 1.2, alpha = 0.6, 
             position = position_jitter(width = 0.15)) +  
  #geom_smooth(aes(group = Protocol, color = Protocol), method = "lm", se = FALSE, size = 0.6, linetype = "solid", color="black") +
  labs(x = "Dataset", y = "Mean Stress Expression") +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) + 
  scale_y_continuous(limits = c(1.2,1.8), breaks = seq(1.2,1.8, by=0.2))+
  theme_classic() +
  theme(
  )

ggsave(filename = "~/brmatched-data/figs/stress-processing.svg", plot = stressg, device = "svg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brmatched-data/figs/stress-processing.jpg", plot = stressg, device = "jpg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brmatched-data/figs/stress-processing.pdf", plot = stressg, device = "pdf",width = 12, height = 8, dpi = 300)

stress.time<-ggplot(df_long.time2.1, aes(x = Isolation, y = sumexp_stress, fill = Isolation, color=Isolation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.8) +  
  geom_point(shape = 21, fill = "black", size = 1.2, alpha = 0.6, 
             position = position_jitter(width = 0.15)) +  
  #geom_smooth(aes(group = Protocol, color = Protocol), method = "lm", se = FALSE, linewidth = 0.6, linetype = "solid", color="black") +
  labs(x = "Dataset", y = "Mean Stress Expression") +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) + 
  theme_classic() +
  theme(
  )

ggsave(filename = "~/brmatched-data/figs/stress-time.jpg", plot = stress.time, device = "jpg",width = 12, height = 8, dpi = 300)

#Heat map
#HEATMAP

#From META OBJECT matrix
aggregatecounts<-read.csv("~/brmatched-data/brmatched-aggregatecounts.txt", header = TRUE)


# Extract expression matrix
matrix.brmatched<-as.matrix(pb.brmatched[["RNA"]]@layers$data)

#Normalize expression values between 0 and 1


#Crossing feature.immune against matrix
matriz_filtrada <- pb.brmatched[rownames(pb.brmatched) %in% features.stress, ]


log_mat.brmatched<-apply(log1p(matriz_filtrada), 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
})


#matrix wrangling
matriz_filtrada2 <- as.data.frame(matriz_normalizada)%>%
  rownames_to_column("Gene")%>%
  pivot_longer(cols=-Gene, names_to = "Sample", values_to = "Expression")%>%
  filter(Gene %in% features.stress)%>%
  mutate(CellType = sub("_LC*","",Sample),
         Treatment=case_when(
           grepl("Ficoll", Sample)~ "Ficoll",
           grepl("SiteCELL",Sample)~"SiteCELL"
         ))

matriz_filtrada2.1 <- matriz_filtrada2 %>%
  mutate(Sample = str_extract(Sample, "LCBR009[0-9]"))

df<-matriz_filtrada2.1%>%
  separate(CellType, into=c("CellType_clean", NA), sep="BR",remove=FALSE)%>%
  mutate(CellType =CellType_clean)%>%
  select(-CellType_clean)%>%
  mutate(MajorGroup = case_when(
    CellType %in% c("CD8 TEM", "CD4 Naive",
                    "CD8 TCM","CD4 TCM",
                    "CD4 TEM","CD8 Naive","gdT","dnT","ILC","Treg","MAIT","CD4 CTL","CD4 Proliferating","CD8 Proliferating","CD4 Naive") ~  "T",
    CellType %in% c("NK","NK-CD56bright","NK Proliferating")~  "NK",
    CellType %in% c("B memory","B intermediate","B naive","Plasmablast")~  "B",
    CellType %in% c("CD14 Mono","CD16 Mono")~ "Monocyte",
    CellType %in% c("cDC1","cDC2","pDC","ASDC")~ "DC",
    CellType %in% c("Eryth")~ "Eryth",
    CellType %in% c("Platelet")~ "Platelet",
    CellType %in% c("HSPC")~ "HSPC",
    TRUE ~ "NOCLASS"
  ))

df<-df%>%
  mutate(Sample_CellType=paste(MajorGroup, Sample, sep="_"))

df_ficoll<-df%>%
  filter(Treatment=="Ficoll")%>%
  filter(MajorGroup != "Platelet")%>%
  filter(MajorGroup != "Eryth")%>%
  filter(MajorGroup != "HSPC")


df_lc<-df%>%
  filter(Treatment=="SiteCELL")%>%
  filter(MajorGroup != "HSPC")
  

#Set color palette
palette_colors<-c("#762A83FF","white", "#FFA500") #darkviolet","#301E48FF","#F7AA14FF")

#Heatmap Plot
heatmap.br.fic<-ggplot(df_ficoll, aes(x = Sample_CellType, y = Gene, fill = Expression)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  scale_fill_gradientn(
    colors=palette_colors)+
  geom_vline(xintercept = split_index, color="black", linetype="dashed",linewidth=1)+
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

ggsave(filename = "~/brmatched-data/figs/heatmap-ficoll.jpg", plot = heatmap.br.fic, device = "jpg",width = 12, height = 8, dpi = 300)



heatmap.br.lc<-ggplot(df_lc, aes(x = Sample_CellType, y = Gene, fill = Expression)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  scale_fill_gradientn(
    colors=palette_colors,
    limits=c(0,4),
    oob=scales::squish) +  
  geom_vline(xintercept = split_index, color="black", linetype="dashed",linewidth=1)+
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

ggsave(filename = "~/brmatched-data/figs/heatmap-lc.jpg", plot = heatmap.br.lc, device = "jpg",width = 12, height = 8, dpi = 300)

######Proportion calculation
b23.93<-subset(int.B23.ann, subset = Patient == "LCBR0093")
b23.94<-subset(int.B23.ann, subset = Patient == "LCBR0094")
b23.95<-subset(int.B23.ann, subset = Patient == "LCBR0095")
b23.96<-subset(int.B23.ann, subset = Patient == "LCBR0096")
b24.93<-subset(int.B24.seu.ann, subset = Patient == "LCBR0093")
b24.94<-subset(int.B24.seu.ann, subset = Patient == "LCBR0094")
b24.95<-subset(int.B24.seu.ann, subset = Patient == "LCBR0095")
b24.96<-subset(int.B24.seu.ann, subset = Patient == "LCBR0096")

#list counts and props
sample<-c("b23.93","b23.94","b23.95","b23.96","b24.93","b24.94","b24.95","b24.96")

#Create lists for both counts and props
sam_list<-mget(sample)
names(sam_list)<-sample
list_counts<-lapply(sam_list, function(obj){
  table(obj$predicted.celltype.l2)
})

list_props<-lapply(list_counts, prop.table)

#create dataframe
df_summary<-bind_rows(lapply(names(list_counts), function(sample_id){
  data.frame(
    Sample=sample_id,
    CellType=names(list_counts[[sample_id]]),
    Count=as.vector(list_counts[[sample_id]]),
    Proportion=as.vector(list_props[[sample_id]])
  )
}))

df_summary<-df_summary%>%
  mutate(MajorGroup = case_when(
    CellType %in% c("CD8 TEM", "CD4 Naive",
                    "CD8 TCM","CD4 TCM",
                    "CD4 TEM","CD8 Naive","gdT","dnT","ILC","Treg","MAIT","CD4 CTL","CD4 Proliferating","CD8 Proliferating","CD4 Naive") ~  "T",
    CellType %in% c("NK","NK_CD56bright","NK Proliferating")~  "NK",
    CellType %in% c("B memory","B intermediate","B naive","Plasmablast")~  "B",
    CellType %in% c("CD14 Mono","CD16 Mono")~ "Monocyte",
    CellType %in% c("cDC1","cDC2","pDC","ASDC")~ "DC",
    CellType %in% c("Eryth")~ "Eryth",
    CellType %in% c("Platelet")~ "Platelet",
    CellType %in% c("HSPC")~ "HSPC",
    TRUE ~ "NOCLASS"
  ))

df_summary<-df_summary%>%
  mutate(Treatment = case_when(
    Sample %in% c("b23.93", "b23.94",
                  "b23.95","b23.96") ~  "SiteCELL",
    Sample %in% c("b24.93", "b24.94",
                  "b24.95","b24.96")~  "Ficoll",
    TRUE ~ "NOCLASS"
  ))


#Calculating subtype proportions
df_summary <- df_summary %>%
  mutate(Proportion = Proportion * 100)

df_summary_subtype <- df_summary %>%
  group_by(Treatment, CellType) %>%
  summarise(
    mean_prop = mean(Proportion),
    sd_prop   = sd(Proportion),
    n         = n(),
    se_prop   = sd_prop / sqrt(n),
    .groups = "drop")

df_summary_subtype <- df_summary_subtype %>%
  mutate(Treatment = recode(Treatment,
                            "Ficoll"     = "FDG",
                            "SiteCELL" = "SiteCELL"))

#plot subtypes in brmatched experiment
br.subtype.prop<-ggplot(df_summary_subtype, aes(x = reorder(CellType,mean_prop), y = mean_prop, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_classic() +
  #geom_errorbar(aes(ymin=mean_prop,ymax=mean_prop + se_prop),position = position_dodge(width=0.9),  width=0.2, linewidth=0.2)+
  geom_linerange(
    aes(
      ymin = mean_prop,
      ymax = (mean_prop + se_prop)
    ),
    position = position_dodge(width = 0.8),
    linewidth = 0.6
  ) +
  scale_fill_manual(values = c("FDG" = "#a6bddb", "SiteCELL" = "#1c9099")) +
  coord_flip() +
  labs(
    x = "Cell Proportion (%)",
    y = "Cell Subtype",
    fill = "Treatment"
  ) +
  theme(
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      family = "Arial",
      face = "plain",
      size = 11
    ),
    axis.title.x = element_text(
      family = "Arial",
      face = "bold",
      size = 12
    ),
    axis.title.y = element_text(
      family = "Arial",
      face = "bold",
      size = 12
    )
  )

ggsave(filename = "~/brmatched-data/figs/pbmc-subtype-br.jpg", plot = br.subtype.prop, device = "jpg",width = 8, height = 6, dpi = 300)



fwrite(df_summary, "counts-brmatched.txt")

#calculate counts and percentages for majorgroup
counts.bra<-df_summary %>%
  group_by(MajorGroup, Treatment, Sample)%>%
  summarise(count_mg=sum(Count),
            perc_mg=sum(Proportion))

#Calculate mean sd, se, var
counts.br.filt<-counts.bra%>%
  group_by(Treatment, MajorGroup)%>%
  summarise(Mean_perc_mg=mean(perc_mg, na.rm = TRUE),
            SD_perc_mg=sd(perc_mg, na.rm = TRUE),
            Var_perc_mg=var(perc_mg, na.rm = TRUE),
            SE_perc_mg=sd(perc_mg, na.rm = TRUE)/sqrt(sum(!is.na(perc_mg))),
            .groups="drop")


#Add eryth and platelet column for SiteCELL
counts.br.filt1<-counts.br.filt%>%
  bind_rows(
    tibble(
      MajorGroup="Eryth",
      Treatment="SiteCELL",
      Mean_perc_mg=0,
      SD_perc_mg=0,
      Var_perc_mg=0,
      SE_perc_mg=0))

counts.br.filt1<-counts.br.filt1%>%
  bind_rows(
    tibble(
      MajorGroup="Platelet",
      Treatment="SiteCELL",
      Mean_perc_mg=0,
      SD_perc_mg=0,
      Var_perc_mg=0,
      SE_perc_mg=0))

#filter out HSCP
counts.br.filt1<-counts.br.filt1%>%
  filter(MajorGroup != "HSPC")


###BARPLOTS
colores <- c('SiteCELL' = '#1c9099',  
             'Ficoll' = '#a6bddb')

counts.br.filt1.b<-counts.br.filt1 %>%
  filter(MajorGroup=="B")

b.overallmean<-mean(counts.br.filt1.b$Mean_perc_mg, na.mr=TRUE)

bcell<-ggplot(counts.br.filt1.b, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  #geom_point(aes(MajorGroup),size=2.2)+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = b.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5) +
  coord_flip() +
  theme_minimal()+
  scale_fill_manual(values=colores)+
  labs(x = "Dataset",
       y = "B Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6, face = "bold"),             
        axis.text = element_text(size = 6),               
        legend.text = element_text(size = 6)
  ) 

ggsave(filename = "~/brmatched-data/figs/Bcellproportion.svg", plot = bcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brmatched-data/figs/Bcellproportion-2.jpg", plot = bcell, device = "jpg",width = 6.15, height = 4.32, dpi = 300, units = "cm")



######NK CELL

counts.br.filt1.n<-counts.br.filt1 %>%
  filter(MajorGroup=="NK")

n.overallmean<-mean(counts.br.filt1.n$Mean_perc_mg, na.mr=TRUE)

nkcell<-ggplot(counts.br.filt1.n, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  #geom_point(aes(MajorGroup),size=2.2)+
  coord_flip()+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  
  geom_hline(yintercept = n.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5)+
  theme_minimal()+
  scale_fill_manual(values = colores)+
  labs(x = "Dataset",
       y = "NK Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),               
        axis.text = element_text(size = 6),                
        legend.text = element_text(size = 6)
  )


ggsave(filename = "~/brmatched-data/figs/NKcellproportion.jpg", plot = nkcell, device = "jpg",width = 8, height = 6, dpi = 300)


#####MONOCYTE
counts.br.filt1.m<-counts.br.filt1 %>%
  filter(MajorGroup=="Monocyte")

m.overallmean<-mean(counts.br.filt1.m$Mean_perc_mg, na.mr=TRUE)  


Monocell<-ggplot(counts.br.filt1.m, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(MajorGroup),size=2.2)+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = m.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5)+ 
  theme_minimal()+
  scale_fill_manual(values = colores)+
  labs(x = "Dataset",
       y = "Monocyte Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),               
        axis.text = element_text(size = 6),              
        legend.text = element_text(size = 6))

ggsave(filename = "~/brmatched-data/figs/Monocellproportion.jpg", plot = Monocell, device = "jpg",width = 8, height = 6,dpi=300)


############T CELL
counts.br.filt1.t<-counts.br.filt1 %>%
  filter(MajorGroup=="T")

t.overallmean<-mean(counts.br.filt1.t$Mean_perc_mg, na.mr=TRUE)

Tcell<-ggplot(counts.br.filt1.t, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(MajorGroup),size=2.2)+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = t.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5) +
  theme_minimal()+
  scale_fill_manual(values=colores)+
  labs(x = "Dataset",
       y = "T Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),              
        axis.text = element_text(size = 6),                
        legend.text = element_text(size = 6))

ggsave(filename = "~/brmatched-data/figs/Tcellproportion.jpg", plot = Tcell, device = "jpg",width = 8, height = 6, dpi = 300)

###########ERYTHROCYTE
counts.br.filt1.e<-counts.br.filt1 %>%
  filter(MajorGroup=="Eryth")
fwrite(counts.br.filt1.e, "~/data/table-eryth.txt")

#Some SiteCELL dataset does not have eryth cell type, we need to aggregate a representative column 
counts.br.filt1.e<-fread("~/data/table-eryth.txt")

e.overallmean<-mean(counts.br.filt1.e$Mean_perc_mg, na.mr=TRUE)

erythcell<-ggplot(counts.br.filt1.e, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = e.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5) +
  theme_minimal()+
  scale_fill_manual(values = colores)+
  labs(x = "Dataset",
       y = "Eryth Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),               
        axis.text = element_text(size = 6),              
        legend.text = element_text(size = 6))

ggsave(filename = "~/brmatched-data/figs/Erythcellproportion.jpg", plot = erythcell, device = "jpg",width = 8, height = 6, dpi = 300)

##############3PLATELETS
counts.br.filt1.p<-counts.br.filt1 %>%
  filter(MajorGroup=="Platelet")

#SiteCELL2 does not present any platelet count, we must to add a representative mock column
dt2 <- data.frame(Treatment = "SiteCELL2", MajorGroup = "Platelet", Mean_perc_mg=0, SD_perc_mg=0, Var_perc_mg=0, SE_perc_mg=0)
counts.br.filt1.p.2 <- rbind(counts.br.filt1.p, dt2)


p.overallmean<-mean(counts.br.filt1.p$Mean_perc_mg, na.mr=TRUE)

plateletcell<-ggplot(counts.br.filt1.p, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(MajorGroup),size=2.2)+
  geom_errorbar(aes(ymin=Mean_perc_mg, ymax=Mean_perc_mg+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = p.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5)+
  theme_minimal()+
  scale_fill_manual(values = colores)+
  labs(x = "Dataset",
       y = "Platelet Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),               
        axis.text = element_text(size = 6),                
        legend.text = element_text(size = 6))

ggsave(filename = "~/brmatched-data/figs/Plateletcellproportion.jpg", plot = plateletcell, device = "jpg",width = 8, height = 6, dpi = 300)

###ALL
#filter only platelets and eryth because of minor scale
counts.br.filt1.c.noep<-counts.br.filt1.c %>%
  filter(MajorGroup %in% c("Eryth", "Platelet"))

#filter out eryth and platelets
counts.br.filt1.c <- counts.br.filt1.c %>%
  filter(MajorGroup != "Platelet")%>%
  filter(MajorGroup != "Eryth")%>%
  group_by(MajorGroup) %>%
  mutate(MeanTotal = mean(Mean_perc_mg)) %>%
  ungroup() %>%
  mutate(MajorGroup = reorder(MajorGroup, -MeanTotal))  # el signo "-" ordena de mayor a menor

#plotting only PBMCs
br.pbmc.prop<-ggplot(counts.br.filt1.c, aes(x = MajorGroup, y = Mean_perc_mg, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(x = "Cell Type", y = "Mean Proportion (%)", fill = "Treatment") +
  theme_classic() +
  geom_errorbar(aes(ymin=Mean_perc_mg - SD_perc_mg, 
                    ymax=Mean_perc_mg + SD_perc_mg),
                position = position_dodge(width=0.9),
                width=0.2, linewidth=0.2)+
  scale_fill_manual(values = c("Ficoll" = "#a6bddb", "SiteCELL" = "#1c9099")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "~/brmatched-data/figs/pbmc-prop-br.jpg", plot = br.pbmc.prop, device = "jpg",width = 8, height = 6, dpi = 300)

#plotting only contamination
#assign NAs to zeros
counts.br.filt1.c.noep[is.na(counts.br.filt1.c.noep)]<-0

br.cont.prop<-ggplot(counts.br.filt1.c.noep, aes(x = MajorGroup, y = Mean_perc_mg, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(x = "Cell Type", y = "Mean Proportion (%)", fill = "Treatment") +
  scale_y_continuous(limits = c(0, 0.015),
                     breaks = seq(0, 0.015, by = 0.0050))+
  geom_errorbar(aes(ymin=Mean_perc_mg, 
                    ymax=Mean_perc_mg + SD_perc_mg),
                position = position_dodge(width=0.9),
                width=0.2, linewidth=0.2)+
  theme_classic() +
  scale_fill_manual(values = c("Ficoll" = "#a6bddb", "SiteCELL" = "#1c9099")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "~/brmatched-data/figs/cont-prop-br.jpg", plot = br.cont.prop, device = "jpg",width = 8, height = 6, dpi = 300)


meta.mtRNA.brmatched.1<-meta.mtRNA.brmatched %>%
  group_by(Protocol) %>%
  summarise(
    mean_mtRNA = mean(mitoPercent, na.rm = TRUE), 
    sd_mtRNA = sd(mitoPercent, na.rm = TRUE))

meta.mtRNA.brmatched.1<-meta.mtRNA.brmatched %>%
  group_by(Protocol) %>%
  summarise(n_mito_high=sum(mitoPercent>15, na.mr=TRUE))
    
