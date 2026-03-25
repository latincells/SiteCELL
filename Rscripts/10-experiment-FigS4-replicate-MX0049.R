install.packages("ggdist")
install.packages("ggridges")
install.packages("ggbeeswarm")

library(ggdist)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggbeeswarm)

######################## MITOPERCENT FIGURE

merged.brazil@meta.data<-merged.brazil@meta.data%>%
  mutate(Dataset_ID="merged-B23")

merged.B24.seu@meta.data<-merged.B24.seu@meta.data%>%
  mutate(Dataset_ID="merged-B24")

merged.brazil<-merge(x = merged.brazil, y=merged.B24.seu)

saveRDS(merged.brazil, "~/brazil-data/RDSs/brazil-merged.rds")
mtRNA.brazil<-readRDS("~/brazil-data/RDSs/brazil-merged.rds")

#Create a meta file from merged object
meta.mtRNA.brazil<-mtRNA.brazil@meta.data

#Collapse mitoPercent + percent.mito 
meta.mtRNA.brazil$mtRNA.collapsed <- ifelse(is.na(meta.mtRNA.brazil$mitoPercent), meta.mtRNA.brazil$percent.mito, meta.mtRNA.brazil$mitoPercent)

#Save all merged metadata (for mtRNA analyses)
fwrite(meta.mtRNA.brazil, "~/brazil-data/mtRNA-brazil.txt")

#Data wrangling

#Add Patient column based on Dataset ID values
meta.mtRNA.brazil<-meta.mtRNA.brazil%>%
  mutate(Protocol_isolation = case_when(
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B23" ~ "LatinCells",
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B24" ~ "Ficoll",
    TRUE ~ "Sin clasificación"  # Por si hay algún valor fuera de los especificados
  ))

colores.stress <- c('LatinCells' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


meta.mtRNA.brazil$Dataset_ID <- factor(
  meta.mtRNA.brazil$Dataset_ID,
  levels = c("merged-B23", "merged-B24"),
  labels = c("LatinCells", "Ficoll")
)

#Con líneas de tendencia
# 1. Ordenar Dataset_IDs según Protocol_isolation (Ficoll primero)
meta.mtRNA.brazil <- meta.mtRNA.brazil %>%
  mutate(Dataset_ID = factor(Dataset_ID, levels = Dataset_ID[order(Protocol_isolation)]))


br.mtplot<-ggplot(meta.mtRNA.brazil, aes(x = Protocol, y = mitoPercent, 
                                      fill = Protocol_isolation, color = Protocol_isolation)) + 
  #geom_jitter(position = position_jitter(width = 0.25), size = 0.8, alpha = 0.8) +
  #geom_violin(width = 0.3, size = 0.7, col = "black", alpha = 0.9) +
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  stat_summary(aes(group = Dataset_ID), fun = "mean", geom = "point", shape = 20,
               size = 1, fill = "white", color = "black", stroke = 1.2) +
  #geom_hline(yintercept = mean(meta.mtRNA.brazil$mitoPercent, na.rm = TRUE), linetype = "dashed", color = "black", alpha=0.7) +
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
ggsave(filename = "figs/mtRNA.svg", plot = br.mtplot, device = "svg", width = 12, height = 9.5, dpi=600)


##################UMIS and GENES per 

colores.stress <- c('LatinCells' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


# Ncount
meta.mtRNA.brazil$Protocol <- ifelse(grepl("Ficoll", meta.mtRNA.brazil$Group), "Ficoll", "LatinCells")


meta.mtRNA.brazil <- meta.mtRNA.brazil %>%
  arrange(Protocol, Group) %>%
  mutate(Group = factor(Group, levels = unique(Group)))


# Graficar
ncount.br<-ggplot(meta.mtRNA.brazil, aes(x = Group, y = nCount_RNA, fill = Protocol_isolation)) + 
  #geom_jitter(shape = 21,    size = 0.9,    alpha = 0.5,    stroke = 0.3,    position = position_jitter(width = 0.15))+
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  scale_fill_manual(values = c('LatinCells' = '#1c9099', 'Ficoll' = '#a6bddb')) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +
  #geom_hline(yintercept = mean(meta.mtRNA.brazil$nCount_RNA, na.rm = TRUE), linetype = "dashed", color = "black", alpha=0.7) +
  scale_x_discrete(labels = function(x) gsub("_.*", "", x)) + 
  scale_y_log10(breaks=c(1000,10000,25000,100000), labels=scales::comma)+
  theme_classic(base_size = 14) +
  labs(x = "Patient", y = "nCount_RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(filename = "figs/ncount-br.jpg", plot = ncount.br, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "figs/ncount-br.svg", plot = ncount.br, device = "svg", width = 12, height = 9.5, dpi=600)


#NFeatures
nfeatu.br<-ggplot(meta.mtRNA.brazil, aes(x = Group, y = nFeature_RNA, fill = Protocol_isolation)) + 
  #geom_jitter(shape = 21,size = 0.9, alpha = 0.5,stroke = 0.3, position = position_jitter(width = 0.15)) +
  geom_violin(width = 0.7, color = NA, alpha = 0.9) +
  scale_fill_manual(values = c('LatinCells' = '#1c9099', 'Ficoll' = '#a6bddb')) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +
  scale_x_discrete(labels = function(x) gsub("_.*", "", x)) +
  scale_y_log10(breaks=c(1000,5000,10000), labels=scales::comma) +
  theme_classic(base_size = 14) +
  labs(x = "Patient", y = "nFeature_RNA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")


ggsave(filename = "figs/nfeature-br.jpg", plot = nfeatu.br, device = "jpg", width = 12, height = 9.5, dpi=600)
ggsave(filename = "figs/nfeature-br.svg", plot = nfeatu.br, device = "svg", width = 12, height = 9.5, dpi=600)

############################### STRESS GENE FIGURE

#Add Patient column based on Dataset ID values
merged.brazil@meta.data<-merged.brazil@meta.data%>%
  mutate(Protocol_isolation = case_when(
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B23" ~ "LatinCells",
    Patient %in% c("LCBR0093", "LCBR0094", "LCBR0095", "LCBR0096") & Dataset_ID == "merged-B24" ~ "Ficoll",
    TRUE ~ "Sin clasificación"  # Por si hay algún valor fuera de los especificados
  ))

merged.brazil<-NormalizeData(merged.brazil)
merged.brazil<-FindVariableFeatures(merged.brazil, selection.method = "vst", nfeatures = 3000)
merged.brazil<-ScaleData(merged.brazil)
merged.brazil<-RunPCA(merged.brazil)
#ElbowPlot(merged.brazil, ndims = 50)
merged.brazil<-FindNeighbors(merged.brazil, dims=1:35)
merged.brazil<-FindClusters(merged.brazil, resolution = 1) 
merged.brazil<-RunUMAP(merged.brazil, dims=1:35, return.model = T)

DimPlot(merged.brazil, label=T, reduction = "umap", group.by = c("Protocol_isolation","seurat_clusters"))


int.mbrasil<-IntegrateLayers(merged.brazil, method = CCAIntegration, 
                         orig.reduction = "pca", new.reduction="integrated.cca",
                         verbose=FALSE)

int.mbrasil<-FindNeighbors(int.mbrasil, reduction = "integrated.cca", dims=1:30)
int.mbrasil<-FindClusters(int.mbrasil, resolution = 1)

int.mbrasil<-RunUMAP(int.mbrasil, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(int.mbrasil.ann, reduction = "umap.cca", group.by=c("seurat_clusters"),alpha = 0.5, pt.size = 1.2, label=TRUE)

int.mbrasil.j<-JoinLayers(int.mbrasil)

int.mbrasil.ann<-RunAzimuth(int.mbrasil.j, reference="pbmcref")
DimPlot(int.mbrasil.ann, reduction = "umap.cca", group.by=c("seurat_clusters", "predicted.celltype.l2"),alpha = 0.5, pt.size = 1.2, label=TRUE)

umaps.by.protocol<-DimPlot(int.mbrasil.ann, reduction = "umap.cca",group.by = "predicted.celltype.l2", split.by="Protocol_isolation",alpha = 0.5, pt.size = 1.2, label=FALSE)
ggsave(filename = "~/brazil-data/figs/umaps-both-protocol-br.jpg", plot = umaps.by.protocol, device = "jpg",width = 12, height = 8, dpi = 600)
ggsave(filename = "~/brazil-data/figs/umaps-both-protocol-br.svg", plot = umaps.by.protocol, device = "svg",width = 12, height = 8, dpi = 600)
ggsave(filename = "~/brazil-data/figs/umaps-both-protocol-br.pdf", plot = umaps.by.protocol, device = "pdf",width = 12, height = 8, dpi = 600)


brasil.latincells<-subset(int.mbrasil.ann, subset = Protocol_isolation =="LatinCells")
brasil.latincells<-FindNeighbors(brasil.latincells,reduction = "integrated.cca", dims = 1:30)
brasil.latincells<-FindClusters(brasil.latincells, resolution = 1)
brasil.latincells<-RunUMAP(brasil.latincells, reduction = "integrated.cca", dims=1:30, reduction.name="umap.cca")

brasil.ficoll<-subset(int.mbrasil.ann, subset = Protocol_isolation =="Ficoll")
brasil.ficoll<-FindNeighbors(brasil.ficoll,reduction = "integrated.cca", dims = 1:30)
brasil.ficoll<-FindClusters(brasil.ficoll, resolution = 1)
brasil.ficoll<-RunUMAP(brasil.ficoll, reduction = "integrated.cca", dims=1:30, reduction.name="umap.cca")


umap.br.lc<-DimPlot(brasil.latincells, reduction = "umap.cca", group.by="predicted.celltype.l2",alpha = 0.5, pt.size = 1.2, label=TRUE)+NoLegend()
ggsave(filename = "~/brazil-data/figs/umap-br-lc.jpg", plot = umap.br.lc, device = "jpg",width = 12, height = 8, dpi = 600)

umap.br.fic<-DimPlot(brasil.ficoll, reduction = "umap.cca", group.by="predicted.celltype.l2",alpha = 0.5, pt.size = 1.2, label=TRUE)
ggsave(filename = "~/brazil-data/figs/umap-br-fic.jpg", plot = umap.br.fic, device = "jpg",width = 12, height = 8, dpi = 600)

#Adding extra Majorgroup celltype

int.mbrasil.ann@meta.data<-int.mbrasil.ann@meta.data%>%
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
pb.brasil<-AggregateExpression(int.mbrasil.ann, return.seurat = T,
                                    assays = "RNA",
                                    group.by = c("predicted.celltype.l2","Patient", "Protocol_isolation"))

#Getting matrix
matrix.brasil<-as.matrix(GetAssayData(pb.brasil, layer = 'data')[, WhichCells(pb.brasil)])

#Write objects
fwrite(matrix.brasil, "~/brazil-data/brasil-aggregatecounts.txt", row.names = TRUE, col.names = TRUE)

#Read gene vestors
#Genes from Massoni-Badossa et al.
features.stress<-c("CIRBP", "H3F3B", "EIF1", "NEAT1", "FTH1", "SRGN", "RBM3", "SAT1", "CXCR4")

#Genes from Savage et al.
features.time <- c("H3F3B",  "NEAT1", "SRGN", "SERF2","CFL1","MYL12A",
                   "CD52","PFN1","DUSP1","FOS","JUNB","JUND",
                   "KLF6","SAT1","GIMAP7")

#Read object
aggregatecounts<-read.csv("~/brazil-data/brasil-aggregatecounts.txt", header = TRUE)

#IF getwd()
#"/data/users/aaron/tesis/Protocol/thesis"

#aggregatecounts<-read.table("pseudobulk_partial.tsv", sep = "\t", header=TRUE)
#View(aggregatecounts)

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
    Sample %in% c("LCBR0093.LC","LCBR0094.LC","LCBR0095.LC","LCBR0096.LC")~  "LatinCells",
    TRUE ~ "NA"  
  ))


df_long.time <- df_time %>%
  pivot_longer(-features.time, names_to = "Sample", values_to = "Expression")

#Creating new "Dataset" column 
df_long.time<-df_long.time%>%
  mutate(Isolation = case_when(
    Sample %in% c("LCBR0093.Ficoll", "LCBR0094.Ficoll",
                  "LCBR0095.Ficoll","LCBR0096.Ficoll") ~  "Ficoll",
    Sample %in% c("LCBR0093.LC","LCBR0094.LC","LCBR0095.LC","LCBR0096.LC")~  "LatinCells",
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


colores.stress <- c('LatinCells' = '#1c9099',  
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

ggsave(filename = "~/brazil-data/figs/stress-processing.svg", plot = stressg, device = "svg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brazil-data/figs/stress-processing.jpg", plot = stressg, device = "jpg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brazil-data/figs/stress-processing.pdf", plot = stressg, device = "pdf",width = 12, height = 8, dpi = 300)

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

ggsave(filename = "~/brazil-data/figs/stress-time.svg", plot = stress.time, device = "svg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brazil-data/figs/stress-time.jpg", plot = stress.time, device = "jpg",width = 12, height = 8, dpi = 300)

#Heat map
#HEATMAP

#From META OBJECT matrix
aggregatecounts<-read.csv("~/brazil-data/brasil-aggregatecounts.txt", header = TRUE)


# Extract expression matrix
matrix.brasil<-as.matrix(pb.brasil[["RNA"]]@layers$data)

#Normalize expression values between 0 and 1


#Crossing feature.immune against matrix
matriz_filtrada <- pb.brasil[rownames(pb.brasil) %in% features.stress, ]


log_mat.brasil<-apply(log1p(matriz_filtrada), 2, function(x) {
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
           grepl("LatinCells",Sample)~"LatinCells"
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
  filter(Treatment=="LatinCells")%>%
  filter(MajorGroup != "HSPC")
  

#Set color palette
palette_colors<-c("#762A83FF","white", "#FFA500") #darkviolet","#301E48FF","#F7AA14FF")

#Heatmap Plot
heatmap.bra.fic<-ggplot(df_ficoll, aes(x = Sample_CellType, y = Gene, fill = Expression)) +
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

ggsave(filename = "~/brazil-data/figs/heatmap-ficoll.svg", plot = heatmap.bra.fic, device = "svg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brazil-data/figs/heatmap-ficoll.jpg", plot = heatmap.bra.fic, device = "jpg",width = 12, height = 8, dpi = 300)



heatmap.bra.lc<-ggplot(df_lc, aes(x = Sample_CellType, y = Gene, fill = Expression)) +
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

ggsave(filename = "~/brazil-data/figs/heatmap-lc.svg", plot = heatmap.bra.lc, device = "svg",width = 12, height = 8, dpi = 300)
ggsave(filename = "~/brazil-data/figs/heatmap-lc.jpg", plot = heatmap.bra.lc, device = "jpg",width = 12, height = 8, dpi = 300)

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
                  "b23.95","b23.96") ~  "LatinCells",
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
                            "LatinCells" = "SiteCELL"))

#plot subtypes in brazil experiment
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

ggsave(filename = "~/brazil-data/figs/pbmc-subtype-br.svg", plot = br.subtype.prop, device = "svg",width = 8, height = 6, dpi = 300)
ggsave(filename = "~/brazil-data/figs/pbmc-subtype-br.jpg", plot = br.subtype.prop, device = "jpg",width = 8, height = 6, dpi = 300)
ggsave(filename = "~/brazil-data/figs/pbmc-subtype-br.pdf", plot = br.subtype.prop, device = cairo_pdf,width = 8, height = 6, dpi = 300)


fwrite(df_summary, "counts-brasil.txt")

#calculate counts and percentages for majorgroup
counts.bra<-df_summary %>%
  group_by(MajorGroup, Treatment, Sample)%>%
  summarise(count_mg=sum(Count),
            perc_mg=sum(Proportion))

#Calculate mean sd, se, var
counts.bra.filt<-counts.bra%>%
  group_by(Treatment, MajorGroup)%>%
  summarise(Mean_perc_mg=mean(perc_mg, na.rm = TRUE),
            SD_perc_mg=sd(perc_mg, na.rm = TRUE),
            Var_perc_mg=var(perc_mg, na.rm = TRUE),
            SE_perc_mg=sd(perc_mg, na.rm = TRUE)/sqrt(sum(!is.na(perc_mg))),
            .groups="drop")


#Add eryth and platelet column for latincells
counts.bra.filt1<-counts.bra.filt%>%
  bind_rows(
    tibble(
      MajorGroup="Eryth",
      Treatment="LatinCells",
      Mean_perc_mg=0,
      SD_perc_mg=0,
      Var_perc_mg=0,
      SE_perc_mg=0))

counts.bra.filt1<-counts.bra.filt1%>%
  bind_rows(
    tibble(
      MajorGroup="Platelet",
      Treatment="LatinCells",
      Mean_perc_mg=0,
      SD_perc_mg=0,
      Var_perc_mg=0,
      SE_perc_mg=0))

#filter out HSCP
counts.bra.filt1<-counts.bra.filt1%>%
  filter(MajorGroup != "HSPC")


###BARPLOTS
colores <- c('LatinCells' = '#1c9099',  
             'Ficoll' = '#a6bddb')

counts.bra.filt1.b<-counts.bra.filt1 %>%
  filter(MajorGroup=="B")

b.overallmean<-mean(counts.bra.filt1.b$Mean_perc_mg, na.mr=TRUE)

bcell<-ggplot(counts.bra.filt1.b, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
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

ggsave(filename = "~/brazil-data/figs/Bcellproportion.svg", plot = bcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/Bcellproportion-2.jpg", plot = bcell, device = "jpg",width = 6.15, height = 4.32, dpi = 300, units = "cm")



######NK CELL

counts.bra.filt1.n<-counts.bra.filt1 %>%
  filter(MajorGroup=="NK")

n.overallmean<-mean(counts.bra.filt1.n$Mean_perc_mg, na.mr=TRUE)

nkcell<-ggplot(counts.bra.filt1.n, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
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


ggsave(filename = "~/brazil-data/figs/NKcellproportion.svg", plot = nkcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/NKcellproportion.jpg", plot = nkcell, device = "jpg",width = 8, height = 6, dpi = 300)


#####MONOCYTE
counts.bra.filt1.m<-counts.bra.filt1 %>%
  filter(MajorGroup=="Monocyte")

m.overallmean<-mean(counts.bra.filt1.m$Mean_perc_mg, na.mr=TRUE)  


Monocell<-ggplot(counts.bra.filt1.m, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
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

ggsave(filename = "~/brazil-data/figs/Monocellproportion.svg", plot = Monocell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/Monocellproportion.jpg", plot = Monocell, device = "jpg",width = 8, height = 6,dpi=300)


############T CELL
counts.bra.filt1.t<-counts.bra.filt1 %>%
  filter(MajorGroup=="T")

t.overallmean<-mean(counts.bra.filt1.t$Mean_perc_mg, na.mr=TRUE)

Tcell<-ggplot(counts.bra.filt1.t, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
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

ggsave(filename = "~/brazil-data/figs/Tcellproportion.svg", plot = Tcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/Tcellproportion.jpg", plot = Tcell, device = "jpg",width = 8, height = 6, dpi = 300)

###########ERYTHROCYTE
counts.bra.filt1.e<-counts.bra.filt1 %>%
  filter(MajorGroup=="Eryth")
fwrite(counts.bra.filt1.e, "~/data/table-eryth.txt")

#Some LatinCells dataset does not have eryth cell type, we need to aggregate a representative column 
counts.bra.filt1.e<-fread("~/data/table-eryth.txt")

e.overallmean<-mean(counts.bra.filt1.e$Mean_perc_mg, na.mr=TRUE)

erythcell<-ggplot(counts.bra.filt1.e, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(MajorGroup),size=2.2)+
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

ggsave(filename = "~/brazil-data/figs/Erythcellproportion.svg", plot = erythcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/Erythcellproportion.jpg", plot = erythcell, device = "jpg",width = 8, height = 6, dpi = 300)

##############3PLATELETS
counts.bra.filt1.p<-counts.bra.filt1 %>%
  filter(MajorGroup=="Platelet")

#Latincells2 does not present any platelet count, we must to add a representative mock column
dt2 <- data.frame(Treatment = "LatinCells2", MajorGroup = "Platelet", Mean_perc_mg=0, SD_perc_mg=0, Var_perc_mg=0, SE_perc_mg=0)
counts.bra.filt1.p.2 <- rbind(counts.bra.filt1.p, dt2)


p.overallmean<-mean(counts.bra.filt1.p$Mean_perc_mg, na.mr=TRUE)

plateletcell<-ggplot(counts.bra.filt1.p, aes(x = reorder(Treatment, Mean_perc_mg), y = Mean_perc_mg, fill =Treatment, group = Treatment)) +
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

ggsave(filename = "~/brazil-data/figs/Plateletcellproportion.svg", plot = plateletcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/brazil-data/figs/Plateletcellproportion.jpg", plot = plateletcell, device = "jpg",width = 8, height = 6, dpi = 300)

###ALL
#filter only platelets and eryth because of minor scale
counts.bra.filt1.c.noep<-counts.bra.filt1.c %>%
  filter(MajorGroup %in% c("Eryth", "Platelet"))

#filter out eryth and platelets
counts.bra.filt1.c <- counts.bra.filt1.c %>%
  filter(MajorGroup != "Platelet")%>%
  filter(MajorGroup != "Eryth")%>%
  group_by(MajorGroup) %>%
  mutate(MeanTotal = mean(Mean_perc_mg)) %>%
  ungroup() %>%
  mutate(MajorGroup = reorder(MajorGroup, -MeanTotal))  # el signo "-" ordena de mayor a menor

#plotting only PBMCs
br.pbmc.prop<-ggplot(counts.bra.filt1.c, aes(x = MajorGroup, y = Mean_perc_mg, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(x = "Cell Type", y = "Mean Proportion (%)", fill = "Treatment") +
  theme_classic() +
  geom_errorbar(aes(ymin=Mean_perc_mg - SD_perc_mg, 
                    ymax=Mean_perc_mg + SD_perc_mg),
                position = position_dodge(width=0.9),
                width=0.2, linewidth=0.2)+
  scale_fill_manual(values = c("Ficoll" = "#a6bddb", "LatinCells" = "#1c9099")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "~/brazil-data/figs/pbmc-prop-br.svg", plot = br.pbmc.prop, device = "svg",width = 8, height = 6, dpi = 300)
ggsave(filename = "~/brazil-data/figs/pbmc-prop-br.jpg", plot = br.pbmc.prop, device = "jpg",width = 8, height = 6, dpi = 300)

#plotting only contamination
#assign NAs to zeros
counts.bra.filt1.c.noep[is.na(counts.bra.filt1.c.noep)]<-0

br.cont.prop<-ggplot(counts.bra.filt1.c.noep, aes(x = MajorGroup, y = Mean_perc_mg, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(x = "Cell Type", y = "Mean Proportion (%)", fill = "Treatment") +
  scale_y_continuous(limits = c(0, 0.015),
                     breaks = seq(0, 0.015, by = 0.0050))+
  geom_errorbar(aes(ymin=Mean_perc_mg, 
                    ymax=Mean_perc_mg + SD_perc_mg),
                position = position_dodge(width=0.9),
                width=0.2, linewidth=0.2)+
  theme_classic() +
  scale_fill_manual(values = c("Ficoll" = "#a6bddb", "LatinCells" = "#1c9099")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "~/brazil-data/figs/cont-prop-br.svg", plot = br.cont.prop, device = "svg",width = 8, height = 6, dpi = 300)
ggsave(filename = "~/brazil-data/figs/cont-prop-br.jpg", plot = br.cont.prop, device = "jpg",width = 8, height = 6, dpi = 300)


meta.mtRNA.brazil.1<-meta.mtRNA.brazil %>%
  group_by(Protocol) %>%
  summarise(
    mean_mtRNA = mean(mitoPercent, na.rm = TRUE), 
    sd_mtRNA = sd(mitoPercent, na.rm = TRUE))

meta.mtRNA.brazil.1<-meta.mtRNA.brazil %>%
  group_by(Protocol) %>%
  summarise(n_mito_high=sum(mitoPercent>15, na.mr=TRUE))
    