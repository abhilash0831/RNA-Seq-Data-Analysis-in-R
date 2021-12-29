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
#**Grouped Bar plot to compare RNA-Seq and RT-PCR Data**
################################################################################


##Extracting Std Error of the mean for RNASeq data for S1

##Read RNA Seq data S1 (0 vs. 5uM Cu)
RNA_Seq_Data_S1 <- read_excel("DESeq2_results.xlsx", sheet = "DEG_0_vs_5uM")
typeof(RNA_Seq_Data_S1)
view(RNA_Seq_Data_S1)
dim(RNA_Seq_Data_S1) #No. of rows and columns for the data S1

#Remove NA values from data S1
RNA_Seq_Data_S1.1 <- RNA_Seq_Data_S1 %>% drop_na() 
view(RNA_Seq_Data_S1.1)

#Read RT-PCR data for S1
RNASeq_RT5 <- read_excel("RNASeq_RT_S1.xlsx", sheet = "0 vs 5")
view(RNASeq_RT5)

#Merge RNA-Seq and RT-PCR data
RNASeq_SE5 <- merge(RNASeq_RT5, RNA_Seq_Data_S1.1, by = "Gene_ID")
RNASeq_SE5 <- RNASeq_SE %>% select(Gene_ID, log2FC.y, StdErr, P_value.y)
View(RNASeq_SE5)

write.xlsx(RNASeq_SE5, file = "RNASeq_StdErr_5.xlsx", sheetName = "0vs5")

##Draw RTPCR-RNASeq grouped bar chart for S1
RNA_RT_dataS1 <- read_excel("RTPCR_Bar.xlsx", sheet = "Bar_Chart5")
view(RNA_RT_dataS1)

bar_plotS1 <- ggplot(RNA_RT_dataS1, aes(x=as.factor(Gene_ID), 
                                    y=log2FC, fill=Technique))+
  geom_bar(position = position_dodge(), stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=log2FC-SE, ymax=log2FC+SE), 
                width = 0.5, position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, 
                                   hjust = 1, size = 10, colour = "black"),
        axis.text.y = element_text(angle = 90, vjust = 0.2,
                                   hjust = 1, size = 10, colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black', size = 1),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"))+
  xlab("Gene ID") + ylab("log2FC")+
  scale_y_continuous(breaks = seq(-6,6, by = 2))+
  ggtitle("S1 (0 vs. 5uM)")
bar_plotS1


################################################################################
##Extracting Std Error of the mean for RNASeq data for S2

#Read RNA Seq data S2 
RNA_Seq_Data_S2 <- read_excel("DESeq2_results.xlsx", sheet = "DEG_0_vs_15uM")
typeof(RNA_Seq_Data_S2)
view(RNA_Seq_Data_S2)
dim(RNA_Seq_Data_S2) #No. of rows and columns for the data S1

#Remove NA values from data S1
RNA_Seq_Data_S2.1 <- RNA_Seq_Data_S2 %>% drop_na() 
view(RNA_Seq_Data_S2.1)

RNASeq_RT15 <- read_excel("RNASeq_RT_S1.xlsx", sheet = "0 vs 15")
view(RNASeq_RT15)

RNASeq_SE15 <- merge(RNASeq_RT15, RNA_Seq_Data_S2.1, by = "Gene_ID")
RNASeq_SE15 <- RNASeq_SE15 %>% select(Gene_ID, log2FC.y, StdErr, P_value.y)
View(RNASeq_SE15)

write.xlsx(RNASeq_SE15, file = "RNASeq_StdErr_15.xlsx", sheetName = "0vs15")

##Draw RTPCR-RNASeq bar chart for S2
RNA_RT_dataS2 <- read_excel("RTPCR_Bar.xlsx", sheet = "Bar_Chart15")
view(RNA_RT_dataS2)
bar_plotS2 <- ggplot(RNA_RT_dataS2, aes(x=as.factor(Gene_ID), 
                                      y=log2FC, fill=Technique))+
  geom_bar(position = position_dodge(), stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=log2FC-SE, ymax=log2FC+SE), 
                width = 0.5, position = position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, 
                                   hjust = 1, size = 10, colour = "black"),
        axis.text.y = element_text(angle = 90, vjust = 0.2,
                                   hjust = 1, size = 10, colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black', size = 1),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"))+
  xlab("Gene ID") + ylab("log2FC")+
  scale_y_continuous(breaks = seq(-6,6, by = 2))+
  ggtitle("S1 (0 vs. 15uM)")
bar_plotS2

################################################################################
######****END****
################################################################################