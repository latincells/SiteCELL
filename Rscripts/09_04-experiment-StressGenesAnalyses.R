library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(tidyr)
library(dplyr)
#library(tidyverse)
library(stringr)
library(data.table)
#library(hdf5r)
#library(rhdf5)
#library(SeuratData)
#library(SeuratDisk)
library(BiocManager)
#library(BPCells)
library(dittoSeq)
library(Azimuth)

#Read previously integrated object
int.all.objects<-readRDS("~/Data/int.all.objects.rds")

#Adding extra Majorgroup celltype
int.all.objects@meta.data<-int.all.objects@meta.data%>%
  mutate(MajorGroup = case_when(
    idents %in% c("CD8 TEM", "CD4 Naive",
                  "CD8 TCM","CD4 TCM",
                  "CD4 TEM","CD8 Naive","gdT","dnT/ILC","Treg","MAIT","CD4 CTL","T/NK Proliferative") ~  "T",
    idents %in% c("NK","NK_CD56bright","T/NK Proliferative")~  "NK",
    idents %in% c("B Memory","B Intermediate","B Naive","Plasmablast")~  "B",
    idents %in% c("CD14 Mono","CD16 Mono")~ "Monocyte",
    idents %in% c("cDC","pDC/ASDC")~ "DC",
    idents %in% c("Eryth")~ "Eryth",
    idents %in% c("Platelet")~ "Platelet",
    idents %in% c("HSPC")~ "HSPC",
    TRUE ~ "Sin clasificación"  
  ))

###AGGREGATE EXPRESSION
pb.meta.object<-AggregateExpression(int.all.objects, return.seurat = T,
                                     assays = "RNA",
                                     group.by = c("idents","Patient"))

#Getting matrix
matrix.meta.object<-as.matrix(GetAssayData(pb.meta.object, layer = 'data')[, WhichCells(pb.meta.object)])

#Write objects
fwrite(matrix.meta.object, "meta-object-aggregatecounts.txt", row.names = TRUE, col.names = TRUE)

#Read gene vectors
#Genes from Massoni-Badossa et al.
features.stress<-c("CIRBP", "H3F3B", "EIF1", "NEAT1", "FTH1", "SRGN", "RBM3", "SAT1", "CXCR4")

#Genes from Savage et al.
features.time <- c("H3F3B",  "NEAT1", "SRGN", "SERF2","CFL1","MYL12A",
                   "CD52","PFN1","DUSP1","FOS","JUNB","JUND",
                   "KLF6","SAT1","GIMAP7")

#Read object
aggregatecounts<-read.csv("~/meta-object-aggregatecounts.txt", header = TRUE)

# Exclude the 'column1' column
colnames_split <- strsplit(colnames(aggregatecounts)[-1], "_")  
suffixes <- sapply(colnames_split, `[`, 2)
# Create a named list of columns by suffix
col_groups <- split(colnames(aggregatecounts)[-1], suffixes)
# Filter the data frame for the genes of interest
df_filt.stress <- aggregatecounts %>% filter(X %in% features.stress)
df_filt.time<- aggregatecounts %>% filter(X %in% features.time)

# Initialize a data frame to store the results
df_stress <- data.frame(features.stress = df_filt.stress$X)
df_time <- data.frame(features.time = df_filt.time$X)

# Loop through each group and sum the columns
for (suffix in names(col_groups)) {
  cols <- col_groups[[suffix]]
  df_stress[[suffix]] <- rowMeans(df_filt.stress[cols])
}

for (suffix in names(col_groups)) {
  cols <- col_groups[[suffix]]
  df_time[[suffix]] <- rowMeans(df_filt.time[cols])
}

df_long.stress <- df_stress %>%
  pivot_longer(-features.stress, names_to = "Sample", values_to = "Expression")

#Creating new "Dataset" column 
df_long.stress<-df_long.stress%>%
  mutate(Dataset = case_when(
    Sample %in% c("SHD5.NA.SHD5.NA.0.00", "CHI014.B3.CHI014.B3.0.00",
                  "SHD6.NA.SHD6.NA.0.00","SHD3.NA.SHD3.NA.0.00",
                  "SHD1.NA.SHD1.NA.0.00") ~  "Ficoll3",
    Sample %in% c("HC1","HC2","HC3","HC4","HC5")~  "Ficoll4",
    Sample %in% c("p1n2MX0142","p1n2MX0143","p1n23MX0144")~  "LatinCells1",
    Sample %in% c("kawasaki1","kawasaki2","kawasaki3")~ "Ficoll5",
    Sample %in% c("immunop1","immunop2","immunop3","immunop4")~ "Ficoll1",
    Sample %in% c("Chile1","Chile2","Chile3","Chile4")~ "LatinCells3",
    Sample %in% c("5","13","14","15","19")~ "Ficoll2",
    Sample %in% c("25.LC.CO.0025.ES","3.LC.CO.0003.ES","112.LCMX0004","95.LCMX0147")~ "LatinCells2",
    TRUE ~ "NA"  
  ))


df_long.time <- df_time %>%
  pivot_longer(-features.time, names_to = "Sample", values_to = "Expression")

#Creating new "Dataset" column 
df_long.time<-df_long.time%>%
  mutate(Dataset = case_when(
    Sample %in% c("SHD5.NA.SHD5.NA.0.00", "CHI014.B3.CHI014.B3.0.00",
                  "SHD6.NA.SHD6.NA.0.00","SHD3.NA.SHD3.NA.0.00",
                  "SHD1.NA.SHD1.NA.0.00") ~  "Ficoll3",
    Sample %in% c("HC1","HC2","HC3","HC4","HC5")~  "Ficoll4",
    Sample %in% c("p1n2MX0142","p1n2MX0143","p1n23MX0144")~  "LatinCells1",
    Sample %in% c("kawasaki1","kawasaki2","kawasaki3")~ "Ficoll5",
    Sample %in% c("immunop1","immunop2","immunop3","immunop4")~ "Ficoll1",
    Sample %in% c("Chile1","Chile2","Chile3","Chile4")~ "LatinCells3",
    Sample %in% c("5","13","14","15","19")~ "Ficoll2",
    Sample %in% c("25.LC.CO.0025.ES","3.LC.CO.0003.ES","112.LCMX0004","95.LCMX0147")~ "LatinCells2",
  ))

df_long.stress2.1<-df_long.stress%>%
  group_by(Sample,Dataset)%>%
  summarise(sumexp_stress=mean(Expression))

df_long.time2.1<-df_long.time%>%
  group_by(Sample,Dataset)%>%
  summarise(sumexp_stress=mean(Expression))


#reorder factor levels
df_long.stress2.1<-df_long.stress2.1%>%
  mutate(Dataset=reorder(Dataset, sumexp_stress, FUN=mean))

#reorder factor levels
df_long.time2.1<-df_long.time2.1%>%
  mutate(Dataset=reorder(Dataset, sumexp_stress, FUN=mean))

df_long.stress2.1 <- df_long.stress2.1 %>%
  mutate(Dataset=as.character(Dataset))%>%
  mutate(Dataset = ifelse(Sample == "TimeRes", "Ficoll3", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Cxcr4", "Ficoll1", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "HumanAging", "Ficoll4", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Immunophenotyping", "Ficoll2", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Kawasaki", "Ficoll5", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells1", "LatinCells1", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells2", "LatinCells2", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells3", "LatinCells3", Dataset))

df_long.time2.1 <- df_long.time2.1 %>%
  mutate(Dataset=as.character(Dataset))%>%
  mutate(Dataset = ifelse(Sample == "TimeRes", "Ficoll3", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Cxcr4", "Ficoll1", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "HumanAging", "Ficoll4", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Immunophenotyping", "Ficoll2", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "Kawasaki", "Ficoll5", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells1", "LatinCells1", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells2", "LatinCells2", Dataset))%>%
  mutate(Dataset = ifelse(Sample == "LatinCells3", "LatinCells3", Dataset))

df_long.stress2.1$Protocol <- ifelse(grepl("LatinCells", df_long.stress2.1$Dataset), "LatinCells", "Ficoll")
df_long.time2.1$Protocol <- ifelse(grepl("LatinCells", df_long.time2.1$Dataset), "LatinCells", "Ficoll")

colores.stress <- c('LatinCells1' = '#1c9099', 
                    'LatinCells2' = '#1c9099',  
                    'LatinCells3' = '#1c9099',  
                    'Ficoll5' = '#a6bddb',
                    'Ficoll3' = '#a6bddb',
                    'Ficoll2' = '#a6bddb',
                    'Ficoll4' = '#a6bddb',
                    'Ficoll1' = '#a6bddb')

colores.stress <- c('LatinCells' = '#1c9099',  
                    'Ficoll' = '#a6bddb')


# Sort labels
orden_etiquetas_stress <- c("Ficoll5","Ficoll1","Ficoll4","Ficoll2","Ficoll3", "LatinCells3", "LatinCells1","LatinCells2") 
orden_etiquetas_time <- c("Ficoll5","Ficoll1","Ficoll4","Ficoll2","Ficoll3", "LatinCells3", "LatinCells1","LatinCells2") 

# Convert Dataset_ID into sorted factor 
df_long.stress2.1$Dataset <- factor(df_long.stress2.1$Dataset, levels = orden_etiquetas_stress)

df_long.time2.1$Dataset <- factor(df_long.time2.1$Dataset, levels = orden_etiquetas_stress)

##plot for STRESS dependent gene expression
stressg<-ggplot(df_long.stress2.1, aes(x = Dataset, y = sumexp_stress, fill = Protocol, color=Protocol)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.8) +  # Boxplot
  geom_point(shape = 21, fill = "black", size = 1.2, alpha = 0.6, 
             position = position_jitter(width = 0.15)) +  
    #geom_smooth(aes(group = Protocol, color = Protocol), method = "lm", se = FALSE, size = 0.6, linetype = "solid", color="black") +
  labs(x = "Dataset", y = "Mean Stress Expression") +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) + 
  theme_classic() +
  theme(
  )


ggsave(filename = "~/res/Figs/stress-genes.jpg", plot = stressg, device = "jpg",width = 12, height = 8, dpi=300)

##Plot for TIME dependent gene expression
ggplot(df_long.time2.1, aes(x = Dataset, y = sumexp_stress, fill = Protocol, color=Protocol)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.8) +  
  geom_point(shape = 21, fill = "black", size = 1.2, alpha = 0.6, 
             position = position_jitter(width = 0.15)) +  
  labs(x = "Dataset", y = "Mean Stress Expression") +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) + 
  theme_classic() +
  theme(
  )
