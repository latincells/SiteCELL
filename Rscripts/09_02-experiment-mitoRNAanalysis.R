merged.mtRNA<-readRDS("~/merged-global-mitoRNA.rds")

#Create a meta file from merged object
merged.all<-merged.mtRNA@meta.data

#Collapse mitoPercent + percent.mito 
merged.all$mtRNA.collapsed <- ifelse(is.na(merged.all$mitoPercent), merged.all$percent.mito, merged.all$mitoPercent)

#Save all merged metadata (for mtRNA analyses)
fwrite(merged.all, "mtRNA-partial-2.txt")

#Data wrangling

#Add Patient column based on Dataset ID values
merged.all<-merged.all%>%
  mutate(Dataset_ID = case_when(
    Patient %in% c("SHD5_NA,SHD5_NA,0.00", "CHI014_B3,CHI014_B3,0.00",
                   "SHD6_NA,SHD6_NA,0.00","SHD3_NA,SHD3_NA,0.00",
                   "SHD1_NA,SHD1_NA,0.00") ~  "Ficoll3",
    Patient %in% c("HC1","HC2","HC3","HC4","HC5")~  "Ficoll4",
    Patient %in% c("p1n2MX0142","p1n2MX0143","p1n23MX0144")~  "LatinCells1",
    Patient %in% c("kawasaki1","kawasaki2","kawasaki3")~ "Ficoll5",
    Patient %in% c("immunop1","immunop2","immunop3","immunop4")~ "Ficoll1",
    Patient %in% c("0,0","1,1","2,2","3,3")~ "LatinCells3",
    Patient %in% c("5","13","14","15","19")~ "Ficoll2",
    Patient %in% c("LCCO0025","LCCO0003","LCMX0004","LCMX0147")~ "LatinCells2",
    TRUE ~ "NO CLASS"  
  ))

#Modify Dataset ID column
merged.all<-merged.all%>%
  mutate(Isolation = case_when(
    Dataset_ID %in% c("Ficoll3") ~  "Ficoll",
    Dataset_ID %in% c("Ficoll4")~  "Ficoll",
    Dataset_ID %in% c("LatinCells1")~  "LatinCells",
    Dataset_ID %in% c("Ficoll5")~ "Ficoll",
    Dataset_ID %in% c("Ficoll1")~ "Ficoll",
    Dataset_ID %in% c("LatinCells3")~ "LatinCells",
    Dataset_ID %in% c("Ficoll2")~ "Ficoll",
    Dataset_ID %in% c("LatinCells2")~ "LatinCells",
    TRUE ~ "NO CLASS" 
  ))


#####mitoPercent
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


orden_etiquetas <- c("Ficoll5","Ficoll4","Ficoll3","Ficoll1","Ficoll2", "LatinCells3", "LatinCells1","LatinCells2") 

# Convertir Dataset_ID a factor con el orden deseado
merged.all$Dataset_ID <- factor(merged.all$Dataset_ID, levels = orden_etiquetas)

#Con líneas de tendencia
mtplot<-ggplot(merged.all, aes(x = Dataset_ID, y = mitoPercent, fill = Dataset_ID, color = Dataset_ID)) +
  geom_violin(width = 0.3, size = 0.7, alpha = 0.9) +
  geom_point(data = summary_df, aes(x = Dataset_ID, y = mean_mito),
             inherit.aes = FALSE, color = "black", size = 3) +
  scale_fill_manual(values = colores.stress) +
  scale_color_manual(values = colores.stress) +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(x = "Dataset", y = "mtRNA (%)")

ggsave(filename = "~/tesis/Protocol/thesis/meta-objects/meta-objects/Figures/UMAPs-reanalize/mtRNA-noline.jpg", 
       plot = mtplot, device = "jpg", width = 12, height = 7, dpi=600)
