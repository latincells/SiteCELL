

##############UMAPS datasets integration (global) and cell type (idents) Figure 2
int.all.objects<-readRDS("~/data/int.all.objects.rds")

umap1<-DimPlot(int.all.objects, reduction = "umap.cca", group.by = "Dataset",raster = FALSE,label = FALSE, alpha = 0.7)
umap2<-DimPlot(int.all.objects, reduction = "umap.cca", group.by = "idents",raster = FALSE,label = FALSE, label.size = 2.5, alpha = 0.7)
umap1+umap2

ggsave(filename = "~/data/Figures/umapglobaldatasets.jpg", plot = umap1, device = "jpg", width = 12, height = 7, dpi=600)
ggsave(filename = "~/data/Figures/umapglobalidents.jpg", plot = umap2, device = "jpg", width = 12, height = 7, dpi=300)
ggsave(filename = "~/data/Figures/umapglobaldatasets.svg", plot = umap1, device = "svg", width = 12, height = 7, dpi=600)
ggsave(filename = "~/data/Figures/umapglobalidents.svg", plot = umap2, device = "svg", width = 12, height = 7, dpi=300)


###############################
#Major group and all subcell type percentages. Figure 3

#Modified table with cell proportions obtained by calculating celltype.prediction.i2 in integrated object metadata
prop.table<-fread("~/data/proportiontable-integrated-object.txt")
prop.table<-fread("~/tesis/Protocol/Data/dataframes-cellcount/tabla-combinada-rawtimeres.txt")

p.table <- prop.table %>%
  mutate(majorgroup = if_else(celltype %in% c("cDC1", "cDC2", "pDC"), "DC", majorgroup))
  ))

##### Calculates cell percentge by major group
prop.table1<-p.table %>%
  group_by(majorgroup, Dataset_ID, patient)%>%
  summarise(Annotation_count_mg=sum(cellcount),
            Annotation_perc_mg=sum(percentage))

#Calculates the mean and variance of the major cell types proportion by dataset (merges samples)
##Mean and variance of majorgroups per dataset (merges individuals). Takes into account the PERCENTAGE value

#GENERATED FOR MAJOR GROUP PLOT
prop.table2<-prop.table1%>%
  group_by(Dataset_ID, majorgroup)%>%
  summarise(mean=mean(Annotation_perc_mg),
            SD_perc_mg=sd(Annotation_perc_mg),
            Var_perc_mg=var(Annotation_perc_mg),
            SE_perc_mg=sd(Annotation_perc_mg)/sqrt(n()))


fwrite(prop.table2, "prop-means.txt")

#FILTER HSPC  
prop.table3<-prop.table2 %>%
  filter(majorgroup!="HSPC")

prop.table3.1<-prop.table3 %>%
  group_by(majorgroup)%>%
  summarise(General_mean=mean(mean))

#####Celltypes props per dataset
prop.table3.2<-prop.table1 %>% 
  group_by(majorgroup, Dataset_ID) %>% 
  get_summary_stats(Annotation_perc_mg, type="mean_sd")


#Data wrangling prop.table_3
prop.table3.2<-prop.table3.2 %>%
  mutate(Dataset_ID=as.character(Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "TimeResolving", "Ficoll3", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Cxcr4", "Ficoll1", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Human-aging", "Ficoll4", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Immunophenotyping", "Ficoll2", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Kawasaki", "Ficoll5", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Batch1_LC", "LatinCells1", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Batch2_LC", "LatinCells2", Dataset_ID))%>%
  mutate(Dataset_ID = ifelse(Dataset_ID == "Batch3_LC", "LatinCells3", Dataset_ID))

#Set color palette
colores <- c('LatinCells1' = '#1c9099', 
             'LatinCells2' = '#1c9099',  
             'LatinCells3' = '#1c9099',  
             'Ficoll5' = '#a6bddb',
             'Ficoll3' = '#a6bddb',
             'Ficoll2' = '#a6bddb',
             'Ficoll4' = '#a6bddb',
             'Ficoll1' = '#a6bddb')
 

#############Subsetting celltypes
##########B CELL

prop.table_3.b<-prop.table3.2 %>%
  filter(majorgroup=="B")

b.overallmean<-mean(prop.table_3.b$mean, na.mr=TRUE)

bcell<-ggplot(prop.table_3.b, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
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

ggsave(filename = "~/data/Figures/Bcellproportion.svg", plot = bcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/Bcellproportion-2.jpg", plot = bcell, device = "jpg",width = 6.15, height = 4.32, dpi = 300, units = "cm")



######NK CELL

prop.table_3.n<-prop.table3.2 %>%
  filter(majorgroup=="NK")

n.overallmean<-mean(prop.table_3.n$mean, na.mr=TRUE)

nkcell<-ggplot(prop.table_3.n, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  #geom_point(aes(majorgroup),size=2.2)+
  coord_flip()+
  geom_errorbar(aes(ymin=mean, ymax=mean+SD_perc_mg), width=.2,
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


ggsave(filename = "~/data/Figures/NKcellproportion.svg", plot = nkcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/NKcellproportion.jpg", plot = nkcell, device = "jpg",width = 8, height = 6, dpi = 300)


#####MONOCYTE
prop.table_3.m<-prop.table3.1 %>%
  filter(majorgroup=="Mono")

m.overallmean<-mean(prop.table_3.m$mean, na.mr=TRUE)  


Monocell<-ggplot(prop.table_3.m, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+SD_perc_mg), width=.2,
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

ggsave(filename = "Monocellproportion.svg", plot = Monocell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "Monocellproportion.jpg", plot = Monocell, device = "jpg",width = 8, height = 6,dpi=300)


#####DC
prop.table_3.dc<-prop.table3.1 %>%
  filter(majorgroup=="DC")

dc.overallmean<-mean(prop.table_3.dc$mean, na.mr=TRUE)  


dccell<-ggplot(prop.table_3.dc, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+SD_perc_mg), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = dc.overallmean, linetype = "dashed", color = "#1c9099", linewidth = 0.5, alpha=0.5)+ 
  theme_minimal()+
  scale_fill_manual(values = colores)+
  labs(x = "Dataset",
       y = "Dendritic Cell Proportion")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 6,face = "bold"),               
        axis.text = element_text(size = 6),              
        legend.text = element_text(size = 6))

ggsave(filename = "~/data/Figures/dc-cellproportion.svg", plot = dccell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/dc-cellproportion.jpg", plot = dccell, device = "jpg",width = 8, height = 6,dpi=300)


############T CELL
prop.table_3.t<-prop.table3.1 %>%
  filter(majorgroup=="T")

t.overallmean<-mean(prop.table_3.t$mean, na.mr=TRUE)

Tcell<-ggplot(prop.table_3.t, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+SD_perc_mg), width=.2,
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

ggsave(filename = "~/data/Figures/Tcellproportion.svg", plot = Tcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/Tcellproportion.jpg", plot = Tcell, device = "jpg",width = 8, height = 6, dpi = 300)

###########ERYTHROCYTE
prop.table_3.e<-prop.table3.2 %>%
  filter(majorgroup=="Eryth")

fwrite(prop.table_3.e, "table-eryth.txt")

#Some LatinCells dataset does not have eryth cell type, we need to aggregate a representative column 
prop.table_3.e<-fread("~/data/table-eryth.txt")

e.overallmean<-mean(prop.table_3.e$mean, na.mr=TRUE)

erythcell<-ggplot(prop.table_3.e, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
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

ggsave(filename = "~/data/Figures/Erythcellproportion.svg", plot = erythcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/Erythcellproportion.jpg", plot = erythcell, device = "jpg",width = 8, height = 6, dpi = 300)

##############3PLATELETS
prop.table_3.p<-prop.table3.1 %>%
  filter(majorgroup=="Platelet")

#Latincells2 does not present any platelet count, we must to add a representative mock column
dt2 <- data.frame(Dataset_ID = "LatinCells2", majorgroup = "Platelet", mean=0, SD_perc_mg=0, Var_perc_mg=0, SE_perc_mg=0)
prop.table_3.p.2 <- rbind(prop.table_3.p, dt2)


p.overallmean<-mean(prop.table_3.p.2$mean, na.mr=TRUE)

plateletcell<-ggplot(prop.table_3.p.2, aes(x = reorder(Dataset_ID, mean), y = mean, fill =Dataset_ID, group = Dataset_ID)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip()+
  #geom_point(aes(majorgroup),size=2.2)+
  geom_errorbar(aes(ymin=mean, ymax=mean+SD_perc_mg), width=.2,
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

ggsave(filename = "~/data/Figures/Plateletcellproportion.svg", plot = plateletcell, device = "svg",width = 6.15, height = 4.32, dpi = 300, units = "cm")
ggsave(filename = "~/data/Figures/Plateletcellproportion.jpg", plot = plateletcell, device = "jpg",width = 8, height = 6, dpi = 300)




#######ALL SUBCELL TYPES BARPLOT FIGURE 3 
colores_personalizados<-c("#8dd3c7","#addfc1","#cdebbb","#edf8b6","#f6f6b8","#e4e2c3","#d2cfce","#bfbbd9","#cdabbf","#de9aa2",
                          "#f08a84","#ee857a","#cb9297","#a9a0b2","#86aece","#9cb1b8","#c0b299","#e3b379","#f7b762","#e2c364",
                          "#cdce66","#b8da66","#c1da82","#d6d5a5","#ebd0c8","#facde4","#f0d1e1","#e6d4dd","#dcd7da","#d3c9d3")

prop.table<-fread("~/data/proportiontable-integrated-object.txt")

prop.tablez <- prop.table %>%
  mutate(group = interaction(patient, Protocol, sep = " - "))


#BAR PLOT PROPORTIONS-all subtypes
barplot2subtypet<-ggplot(prop.tablez, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2, alpha=1) +  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +  
  scale_fill_manual(values =colores_personalizados) +  
  labs(
    x = "Isolation protocol",
    y = "Cell type proportion (%)",
    fill = "MajorGroup"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "~/data/Figures/barplot2subtypet.jpg", plot = barplot2subtypet, device = "jpg", width = 12, height = 7, dpi=300)
ggsave(filename = "~/data/Figures/barplot2subtypet.svg", plot = barplot2subtypet, device = "svg", width = 12, height = 7, dpi=300)





















