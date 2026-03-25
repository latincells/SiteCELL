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

#Read RDS
TSubset<-readRDS("~/Data/TSubset.rds")

#Recluster subset
TSubset<-FindNeighbors(TSubset, reduction = "integrated.cca", dims=1:35)
TSubset<-FindClusters(TSubset)
TSubset<-RunUMAP(TSubset, reduction = "integrated.cca", dims = 1:35, reduction.name = "umap.cca")

#UMAP
dimplotsubt<-DimPlot(TSubset, reduction = "umap.cca", group.by = "idents", raster = FALSE,label = FALSE, label.size = 2.5, alpha = 0.7)

#Save figures
ggsave(filename = "~/res/Figs/umapTsubtypes.jpg", plot = dimplotsubt, device = "jpg", width = 12, height = 7, dpi=300)


###########
#Bar plot looking for differences in subtype proportions among isolation protocols

#read object 
int.all.objects<-readRDS("~/int-meta-object.rds")
#Creating metadata from t subset
meta.int.t<-int.all.objects@meta.data

#Extract columns
meta.t.filter <- meta.int.t %>% 
  select(Patient, Protocol, Isolation_protocol, idents)

#data wrangling
summary_data <- meta.t.filter %>%
  filter(idents %in% c("CD8 TEM","CD8 TCM","CD4 TEM","CD4 TCM",
                       "CD4 Naive","CD8 Naive","gdT","dnT/ILC","Treg",
                       "MAIT","CD4 CTL","T/NK Proliferative"))%>%
  group_by(Protocol, idents) %>%
  summarise(count = n()) %>%
  group_by(Protocol) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

summary_data$idents <- as.character(summary_data$idents)


summary_data<-summary_data %>%
  rename(Dataset=Protocol)%>%
  mutate(Protocol=ifelse(grepl("Ficoll", Dataset),"Ficoll","LatinCells"))

#Create table for barplot
result<-tapply(summary_data$freq, list(summary_data$idents,summary_data$Protocol), mean)
result.sd<-tapply(summary_data$freq, list(summary_data$idents,summary_data$Protocol), sd)

fwrite(summary_data, "~/res/summary_data.txt")

result<-as.data.frame(result)%>%
  rownames_to_column("Celltype")

result.sd<-as.data.frame(result.sd)%>%
  rownames_to_column("Celltype")

result_long <- result %>%
  pivot_longer(cols = c(Ficoll, LatinCells), names_to = "protocol", values_to = "mean")
result_long.sd <- result.sd %>%
  pivot_longer(cols = c(Ficoll, LatinCells), names_to = "protocol", values_to = "SD")


#Add column 
result_long<-cbind(result_long, result_long.sd%>% select(SD))

# Set color palette
col.props <- c('LatinCells' = '#1c9099',  'Ficoll' = '#a6bddb')

#Add Standard error
Tsubcelltbarplot<-ggplot(result_long, aes(x = reorder(Celltype,mean), y = mean, fill = protocol)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=col.props)+
  geom_errorbar(aes(ymin = mean, ymax = mean + SD), 
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic()+labs(x="T Subcelltype", y="Proportion (%)")

ggsave(filename = "~/res/Figs/barplot2subtypet-diff.jpg", plot = Tsubcelltbarplot, device = "jpg", width = 12, height = 7, dpi=300)
