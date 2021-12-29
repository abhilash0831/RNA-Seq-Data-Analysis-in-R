#Load libraries

#Data translocation
library(writexl)
library(readxl)
library(xlsx)

#Plotting
library(ggplot2)
library(gplots)
library(ggpubr) #'ggplot2' Based Publication Ready Plots
library(grid)
library(ggvenn) #Draw Venn Diagram by 'ggplot2'
library(GOplot) #https://wencke.github.io/
library(RColorBrewer) #Provides color schemes for graphics
library(ggrepel)#provides geoms for ggplot2 to repel overlapping text labels
library(EnhancedVolcano) #For making volcano plots
library(pheatmap) #For making heatmaps
library(factoextra) #to extract and visualize the output of exploratory
#multivariate data analyses (Eg. PCA plot)

#Data wrangling
library(tidyverse) #Collection of packages data analysis packages
library(magrittr) #pipe like operator, %>%, to pipe a value forward into an 
#expression or function call
library(tidyr) #changing the shape (pivoting) and hierarchy 
#(nesting and 'unnesting') of a dataset
library(dplyr) #Set of functions for data manipulation
library(forcats) #for handling categorical variables
library(cluster) #Computes agglomerative hierarchical clustering of the dataset.
library(reshape2) #transform data between wide and long formats.

################################################################################
################################################################################
##**Bar plot of depicting the number of upregulated and downregulated genes**
################################################################################
################################################################################

##Figure 5A (Bar plot of biological processes in S1)
fig_5A <- read_excel("GO_count_BPMF_5_15.xlsx", sheet = "BP 5")
view(fig_5A)
ggplot(fig_5A, aes(fill=Regulation, y=`Biological Process`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,25,by = 5), limits = c(0,25))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("BP_5.png", width = 12, height = 5, dpi = 400, limitsize = FALSE) 


##Figure 5B (Bar plot of Molecular Function in S1)
fig_5B <- read_excel("GO_count_BPMF_5_15.xlsx", sheet = "MF 5")
view(fig_5B)
ggplot(fig_5B, aes(fill=Regulation, y=`Molecular Function`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,70,by = 10), limits = c(0,70))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("MF_5.png", width = 12, height = 5, dpi = 400, limitsize = FALSE)  


##Figure 5C (Bar plot of Cellular component in S1)
fig_5C <- read_excel("GO CC.xlsx", sheet = "GO CC 5")
View(fig_5C)
ggplot(fig_5C, aes(fill=Regulation, y=`Cellular component`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,110,by = 10), limits = c(0,110))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("CC_5.png", width = 12, height = 5, dpi = 400, limitsize = FALSE)

################################################################################
################################################################################
##Figure 5D (Bar plot of biological process in S2)
fig_5D <- read_excel("GO_count_BPMF_5_15.xlsx", sheet = "BP 15")
View(fig_5D)
ggplot(fig_5D, aes(fill=Regulation, y=`Biological Process`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,60,by = 10), limits = c(0,60))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("BP_15.png", width = 12, height = 5, dpi = 400, limitsize = FALSE) 


##Figure 5E (Bar plot of Molecular Function in S2)
fig_5E <- read_excel("GO_count_BPMF_5_15.xlsx", sheet = "MF 15")
View(fig_5E)
ggplot(fig_5E , aes(fill=Regulation, y=`Molecular Function`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,110,by = 10), limits = c(0,110))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("MF_15.png", width = 12, height = 5, dpi = 400, limitsize = FALSE)


##Figure 5F (Bar plot of Cellular component in S2)
fig_5F <- read_excel("GO CC.xlsx", sheet = "GO CC 15")
View(fig_5F)
ggplot(fig_5F, aes(fill=Regulation, y=`Cellular component`, x=`Gene Count`))+
  geom_bar(position="dodge", stat="identity", width = 0.5)+
  scale_x_continuous(breaks = seq(0,220,by = 20), limits = c(0,220))+
  theme(panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))
ggsave("CC_15.png", width = 12, height = 5, dpi = 400, limitsize = FALSE)

################################################################################
######****END****
################################################################################

























