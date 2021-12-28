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
##**Observations are represented as rows** in our data frames or tibbles, while
##the **columns correspond to variables**
################################################################################
################################################################################

## Analyze gene statistics data for sample S1
#S1 refers to pairwise comparison 0 vs 5uM Cu(II) 

#Read DESeq2 output data for S1 (0_vs_5uM) from online galaxy server
RNA_Seq_Data_S1 <- read_excel("DESeq2_results.xlsx", sheet = "DEG_0_vs_5uM")
typeof(RNA_Seq_Data_S1)
view(RNA_Seq_Data_S1)
dim(RNA_Seq_Data_S1) #No. of rows and columns for the data S1

#Remove NA values from data S1
RNA_Seq_Data_S1.1 <- RNA_Seq_Data_S1 %>% drop_na() 
view(RNA_Seq_Data_S1.1)

#Read data as tibble 
RNA_Seq_Data_S1.2 <- tibble(RNA_Seq_Data_S1.1)
dim(RNA_Seq_Data_S1.2) 
view(RNA_Seq_Data_S1.2)

#Total genes removed after dropping NA values from data S1 (0vs5uM)
genes_removed_S1 = nrow(RNA_Seq_Data_S1) - nrow(RNA_Seq_Data_S1.2)
print(genes_removed_S1) #191 genes removed in data S1

#Total genes remaining after dropping NA values from data S1 (0vs5uM)
genes_remaining_S1 = nrow(RNA_Seq_Data_S1) - genes_removed_S1
print(genes_remaining_S1) #3067 genes remaining in data S1

#select_all_genes_with p-value < 0.05
RNA_Seq_Data_S1.3 <- RNA_Seq_Data_S1.2 %>% filter(P_value < 0.05) %>% 
  filter(log2FC >= 0 | log2FC <= 0 )
View(RNA_Seq_Data_S1.3)
dim(RNA_Seq_Data_S1.3)

write.csv(RNA_Seq_Data_S1.3, file = "All_genes_in_5(ONLY p-value).csv"
          , row.names = FALSE)

#Separate dataframe of only up and downregulated genes in data S1
RNA_Seq_Data_S1.3_only_upregulated <- RNA_Seq_Data_S1.3 %>% filter(log2FC > 0)
view(RNA_Seq_Data_S1.3_only_upregulated)
Total_no_upregulated_S1 = nrow(RNA_Seq_Data_S1.3_only_upregulated)
print(Total_no_upregulated_S1) #459 genes upregulated
write.csv(RNA_Seq_Data_S1.3_only_upregulated, 
          file = "up_genes_in_5(p-value).csv", row.names = FALSE)

RNA_Seq_Data_S1.3_only_downregulated <- RNA_Seq_Data_S1.3 %>% filter(log2FC < 0)
view(RNA_Seq_Data_S1.3_only_downregulated)
Total_no_downregulated_S1 = nrow(RNA_Seq_Data_S1.3_only_downregulated)
print(Total_no_downregulated_S1) #731 genes downregulated
write.csv(RNA_Seq_Data_S1.3_only_downregulated, 
          file = "down_genes_in_5(p-value).csv", row.names = FALSE)

#All down and top 50 downregulated genes in data S1 
All_down_0vs5 <- 
  RNA_Seq_Data_S1.3_only_downregulated[order
                                       (RNA_Seq_Data_S1.3_only_downregulated$
                                           log2FC),]
view(All_down_0vs5)
Top50_down_0vs5 <- All_down_0vs5[1:50,]
view(Top50_down_0vs5)
write.csv(Top50_down_0vs5, 
          file = "top_50_down_genes_in_0_vs_5.csv", row.names = FALSE)

#All up and top 50 upregulated genes in data S1 
All_up_0vs5 <- 
  RNA_Seq_Data_S1.3_only_upregulated[order
                                     (-RNA_Seq_Data_S1.3_only_upregulated$
                                         log2FC),]
view(All_up_0vs5)
Top50_up_0vs5 <- All_up_0vs5[1:50,]
view(Top50_up_0vs5)
write.csv(Top50_up_0vs5, 
          file = "top_50_up_genes_in_0_vs_5.csv", row.names = FALSE)



################################################################################
################################################################################
## Analyze gene statistics data for sample S2
#S2 refers to pairwise comparison 0 vs 15uM Cu(II) 

#Read DESeq2 output data for S2(0_vs_15uM) downloaded from online galaxy server
RNA_Seq_Data_S2 <- read_excel("DESeq2_results.xlsx", sheet = "DEG_0_vs_15uM")
view(RNA_Seq_Data_S2)
dim(RNA_Seq_Data_S2) #No. of rows and columns for the data S2

#Remove NA values from data S2
RNA_Seq_Data_S2.1 <- RNA_Seq_Data_S2 %>% drop_na() 
view(RNA_Seq_Data_S2.1)

#Read data as tibble 
RNA_Seq_Data_S2.2 <- tibble(RNA_Seq_Data_S2.1)
dim(RNA_Seq_Data_S2.2) 
view(RNA_Seq_Data_S2.2)

#Total genes removed after dropping NA values from data S2 (0vs15uM)
genes_removed_S2 = nrow(RNA_Seq_Data_S2) - nrow(RNA_Seq_Data_S2.2)
print(genes_removed_S2) #192 genes removed in data S2

#Total genes remaining after dropping NA values from data S2 (0vs15uM)
genes_remaining_S2 = nrow(RNA_Seq_Data_S2) - genes_removed_S2
print(genes_remaining_S2) #3066 genes remaining in data S2

#select_all_genes_with p-value < 0.05 from data S2 (0vs15uM)
RNA_Seq_Data_S2.3 <- RNA_Seq_Data_S2.2 %>% filter(P_value < 0.05) %>% 
  filter(log2FC >= 0 | log2FC <= 0 )
View(RNA_Seq_Data_S2.3)
dim(RNA_Seq_Data_S2.3)

write.csv(RNA_Seq_Data_S2.3, file = "All_genes_in_15(ONLY p-value).csv"
          , row.names = FALSE)

#Split dataframe into two separate dataframes consisting of up and 
#downregulated genes in S2
RNA_Seq_Data_S2.3_only_upregulated <- RNA_Seq_Data_S2.3 %>% filter(log2FC > 0)
view(RNA_Seq_Data_S2.3_only_upregulated)
Total_no_upregulated_S2 = nrow(RNA_Seq_Data_S2.3_only_upregulated)
print(Total_no_upregulated_S2) #987 genes upregulated
write.csv(RNA_Seq_Data_S2.3_only_upregulated, 
          file = "up_genes_in_15(p-value and log2FC).csv", row.names = FALSE)

RNA_Seq_Data_S2.3_only_downregulated <- RNA_Seq_Data_S2.3 %>% filter(log2FC < 0)
view(RNA_Seq_Data_S2.3_only_downregulated)
Total_no_downregulated_S2 = nrow(RNA_Seq_Data_S2.3_only_downregulated)
print(Total_no_downregulated_S2) #968 genes downregulated
write.csv(RNA_Seq_Data_S2.3_only_downregulated, 
          file = "down_genes_in_15(p-value and log2FC).csv", row.names = FALSE)

#Top 50 downregulated genes in data S2
All_down_0vs15 <- 
  RNA_Seq_Data_S2.3_only_downregulated[order
                                       (RNA_Seq_Data_S2.3_only_downregulated$
                                           log2FC),]
view(All_down_0vs15)
Top50_down_0vs15 <- All_down_0vs15[1:50,]
view(Top50_down_0vs15)
write.csv(Top50_down_0vs15, 
          file = "top_50_down_genes_in_0_vs_15.csv", row.names = FALSE)

#Top 50 upregulated genes in data S2 
All_up_0vs15 <- 
  RNA_Seq_Data_S2.3_only_upregulated[order
                                     (-RNA_Seq_Data_S2.3_only_upregulated$
                                         log2FC),]
view(All_up_0vs15)
Top50_up_0vs15 <- All_up_0vs15[1:50,]
view(Top50_up_0vs15)
write.csv(Top50_up_0vs15, 
          file = "top_50_up_genes_in_0_vs_15.csv", row.names = FALSE)

#################################################################################
#################################################################################

## Analyze gene statistics data for sample S3
#S3 refers to pairwise comparison 5 vs 15uM Cu(II) 

#Read DESeq2 output data for S3(3_vs_15uM) downloaded from online galaxy server
RNA_Seq_Data_S3 <- read_excel("DESeq2_results.xlsx", sheet = "DEG_5_vs_15uM")
view(RNA_Seq_Data_S3)
dim(RNA_Seq_Data_S3) #No. of rows and columns for the data S3

#Remove NA values from data S3
RNA_Seq_Data_S3.1 <- RNA_Seq_Data_S3 %>% drop_na() 
view(RNA_Seq_Data_S3.1)

#Read data as tibble 
RNA_Seq_Data_S3.2 <- tibble(RNA_Seq_Data_S3.1)
dim(RNA_Seq_Data_S3.2) 
view(RNA_Seq_Data_S3.2)

#Total genes removed after dropping NA values from data S3 (5vs15uM)
genes_removed_S3 = nrow(RNA_Seq_Data_S3) - nrow(RNA_Seq_Data_S3.2)
print(genes_removed_S3) #198 genes removed in data S3

#Total genes remaining after dropping NA values from data S3 (5vs15uM)
genes_remaining_S3 = nrow(RNA_Seq_Data_S3) - genes_removed_S3
print(genes_remaining_S3) #3060 genes remaining in data S3

#select_all_genes_with p-value < 0.05 from data S3 (5vs15uM)
RNA_Seq_Data_S3.3 <- RNA_Seq_Data_S3.2 %>% filter(P_value < 0.05) %>% 
  filter(log2FC >= 0 | log2FC <= 0 )
View(RNA_Seq_Data_S3.3)
dim(RNA_Seq_Data_S3.3)

write.csv(RNA_Seq_Data_S3.3, file = "All_genes_in_5vs15(ONLY p-value).csv"
          , row.names = FALSE)

#Split dataframe into two separate dataframes consisting of up and 
#downregulated genes in S3
RNA_Seq_Data_S3.3_only_upregulated <- RNA_Seq_Data_S3.3 %>% filter(log2FC > 0)
view(RNA_Seq_Data_S3.3_only_upregulated)
Total_no_upregulated_S3 = nrow(RNA_Seq_Data_S3.3_only_upregulated)
print(Total_no_upregulated_S3) #538 genes upregulated
write.csv(RNA_Seq_Data_S3.3_only_upregulated, 
          file = "up_genes_in_5vs15(p-value and log2FC).csv", row.names = FALSE)

RNA_Seq_Data_S3.3_only_downregulated <- RNA_Seq_Data_S3.3 %>% filter(log2FC < 0)
view(RNA_Seq_Data_S3.3_only_downregulated)
Total_no_downregulated_S3 = nrow(RNA_Seq_Data_S3.3_only_downregulated)
print(Total_no_downregulated_S3) #606 genes downregulated
write.csv(RNA_Seq_Data_S3.3_only_downregulated, 
          file = "down_genes_in_5vs15(p-value and log2FC).csv", row.names = FALSE)

#Top 50 downregulated genes in data S3
All_down_5vs15 <- 
  RNA_Seq_Data_S3.3_only_downregulated[order
                                       (RNA_Seq_Data_S3.3_only_downregulated$
                                           log2FC),]
view(All_down_5vs15)
Top50_down_5vs15 <- All_down_5vs15[1:50,]
view(Top50_down_5vs15)
write.csv(Top50_down_5vs15, 
          file = "top_50_down_genes_in_5_vs_15.csv", row.names = FALSE)

#Top 50 upregulated genes in data S3
All_up_5vs15 <- 
  RNA_Seq_Data_S3.3_only_upregulated[order
                                     (-RNA_Seq_Data_S3.3_only_upregulated$
                                         log2FC),]
view(All_up_5vs15)
Top50_up_5vs15 <- All_up_5vs15[1:50,]
view(Top50_up_5vs15)
write.csv(Top50_up_5vs15, 
          file = "top_50_up_genes_in_5_vs_15.csv", row.names = FALSE)

################################################################################
##**** Preliminary analysis of DEG in the data ****
##****Bar graph, Correlation analysis, Venn Diagram, Volcano plot****##
################################################################################

################################################################################
#Create a new dataframe to draw a bar chart for total count of genes that remain
#after filtering only on p-value (<0.05) in all pairwise samples (S1,S2,S3)
################################################################################
Total_genes_vs_condition <- tibble(Sample_0_vs_5uM = 1190, 
                                   Sample_0_vs_15uM = 1955, 
                                   Sample_5_vs_15uM = 1144)
Total_genes_vs_condition
view(Total_genes_vs_condition)

rownames(Total_genes_vs_condition) <- c("Number_of_genes")

Total_genes_vs_condition1 <- Total_genes_vs_condition %>% 
  pivot_longer(c
               ("Sample_0_vs_5uM",
                 "Sample_0_vs_15uM",
                 "Sample_5_vs_15uM"), 
               names_to = "Sample", values_to = "Total_genes") 

view(Total_genes_vs_condition1)
Total_genes_vs_condition1$Sample <- 
  factor(Total_genes_vs_condition1$Sample,
         levels=c("Sample_0_vs_5uM", "Sample_0_vs_15uM", "Sample_5_vs_15uM"))
ggplot(data = Total_genes_vs_condition1)+
  geom_col(mapping = aes(x=Sample, y=Total_genes, fill = Sample))+
  geom_text(aes(x=Sample, y=Total_genes, label = Total_genes), vjust = -0.5)+
  ylab("Total number of genes")+
  scale_fill_manual(values=c("#EFC000FF", "#868686FF", "#CD534CFF"), 
                    labels = c("SP1", "SP2", "SP3"))+
  theme(panel.background = element_blank())+
  theme_classic()+
  theme(axis.text.x=element_blank())+
  theme(axis.line = element_line(size = 1.2))+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
ggsave("Total_genes_vs_condition.png", width = 5, height = 5)

################################################################################
##****VennDiagram****
################################################################################

##Create a "list" for all samples (S1, S2 nd S3) only on Gene_ID column
options(max.print=2000) #Sets maxm print characters on console to 2000
RNA_Seq_Data_S1.3_vector <- as.vector(unlist(RNA_Seq_Data_S1.3$Gene_ID))
RNA_Seq_Data_S1.3_vector

RNA_Seq_Data_S2.3_vector <- as.vector(unlist(RNA_Seq_Data_S2.3$Gene_ID))
RNA_Seq_Data_S2.3_vector

RNA_Seq_Data_S3.3_vector <- as.vector(unlist(RNA_Seq_Data_S3.3$Gene_ID))
RNA_Seq_Data_S3.3_vector

RNA_Seq_Venn <- list(S1 = RNA_Seq_Data_S1.3_vector,
                     S2 = RNA_Seq_Data_S2.3_vector,
                     S3 = RNA_Seq_Data_S3.3_vector)

ggvenn(
  RNA_Seq_Venn, 
  fill_color = c("#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 8
)

################################################################################
###**Data verification for Venn Diagram**
################################################################################
##Common genes between samples S1 and S2 (should be 1000)
Common_genes_S1_and_S2 <- merge(RNA_Seq_Data_S1.3, RNA_Seq_Data_S2.3,
                                by.x = "Gene_ID", by.y = "Gene_ID")
nrow(Common_genes_S1_and_S2) #1000
View(Common_genes_S1_and_S2)
write.csv(Common_genes_S1_and_S2, file = "Common_genes_S1_and_S2.csv",
          row.names = FALSE)

##Common upregulated genes in S1 and S2
Common_genes_S1_and_S2_up <- Common_genes_S1_and_S2 %>% 
  filter(log2FC.x > 0 & log2FC.y > 0)
view(Common_genes_S1_and_S2_up) #373 commonly upregulated

##Common downregulated genes in S1 and S2
Common_genes_S1_and_S2_down <- Common_genes_S1_and_S2 %>% 
  filter(log2FC.x < 0 & log2FC.y < 0)
view(Common_genes_S1_and_S2_down) #613 commonly downregulated

##Downregulated in S1 and Upregulated in S2
Common_genes_S1down_and_S2up <- Common_genes_S1_and_S2 %>% 
  filter(log2FC.x < 0 & log2FC.y > 0)
view(Common_genes_S1down_and_S2up) #5 contra-regulated

##Upregulated in S1 and Downregulated in S2
Common_genes_S1up_and_S2down <- Common_genes_S1_and_S2 %>% 
  filter(log2FC.x > 0 & log2FC.y < 0)
view(Common_genes_S1up_and_S2down) #9 contra-regulated


##Common genes between samples S1 and S3 (should be 509)
Common_genes_S1_and_S3 <- merge(RNA_Seq_Data_S1.3, RNA_Seq_Data_S3.3,
                                by = "Gene_ID", all = FALSE)
nrow(Common_genes_S1_and_S3)
write.csv(Common_genes_S1_and_S3, file = "Common_genes_S1_and_S2.csv",
          row.names = FALSE)

##Common genes between samples S2 and S3 (should be 976)
Common_genes_S2_and_S3 <- merge(RNA_Seq_Data_S2.3, RNA_Seq_Data_S3.3,
                                by = "Gene_ID", all = FALSE)
nrow(Common_genes_S2_and_S3)
write.csv(Common_genes_S2_and_S3, file = "Common_genes_S1_and_S2.csv",
          row.names = FALSE)

################################################################################
##**Correlation between common DEG in S1 and S2
##**w.r.t control (0uM Cu)**
################################################################################

##Correlation between all common genes in S1 and S2
cor.test(Common_genes_S1_and_S2$log2FC.x, Common_genes_S1_and_S2$log2FC.y, 
         method = "pearson") #R=0.87 between all common genes in S1 and S2

ggplot(data = Common_genes_S1_and_S2, mapping = aes(x = log2FC.x, y=log2FC.y))+
  geom_point(colour = "red", size = 2, alpha = 0.3)+
  geom_smooth(method = "lm")+
  xlab("log2FC (S1)") + ylab("log2FC (S2)")+
  annotate(geom = "text", x = -2, y = 2.5, label = "R = 0.87", color = 'blue',
           size = 6)+
  coord_cartesian(xlim = c(-6,4), ylim = c(-7, 5))+
  theme(panel.background = element_blank())+
  theme_classic()+
  theme(axis.line = element_line(size = 1.2))+
  theme(axis.text = element_text(size = 26),
        axis.title = element_text(size = 26))+
  ggtitle("Correlation between all common \n genes in S1 and S2")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))
 
ggsave("correlation_5_vs_15(p-value and log2FC).png", width = 5, height = 5)

################################################################################
##Correlation between common downregulated genes in S1 and S2
cor.test(Common_genes_S1_and_S2_down$log2FC.x, 
         Common_genes_S1_and_S2_down$log2FC.y,
         method = "pearson") #R = 0.6017653 for common down genes in S1 and S2

ggplot(data = Common_genes_S1_and_S2_down, mapping = 
         aes(x = log2FC.x, y=log2FC.y))+
  geom_point(colour = "red", size = 2, alpha = 0.3)+
  geom_smooth(method = "lm")+
  xlab("log2FC (S1)") + ylab("log2FC (S2)")+
  annotate(geom = "text", x = -2, y = 2.5, label = "R = 0.60", color = 'blue',
           size = 6)+
  coord_cartesian(xlim = c(-6,4), ylim = c(-7, 5))+
  theme(panel.background = element_blank())+
  theme_classic()+
  theme(axis.line = element_line(size = 1.2))+
  theme(axis.text = element_text(size = 26),
        axis.title = element_text(size = 26))+
  ggtitle("Correlation between commonly downregulated \n genes in S1 and S2")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))

################################################################################
##Correlation between common upregulated genes in S1 and S2
cor.test(Common_genes_S1_and_S2_up$log2FC.x, 
         Common_genes_S1_and_S2_up$log2FC.y,
         method = "pearson") #R = 0.3378438 for common up genes in S1 and S2

ggplot(data = Common_genes_S1_and_S2_up, mapping = 
         aes(x = log2FC.x, y=log2FC.y))+
  geom_point(colour = "red", size = 2, alpha = 0.3)+
  geom_smooth(method = "lm")+
  xlab("log2FC (0 vs 5μM Cu)") + ylab("log2FC (0 vs 15μM Cu)")+
  annotate(geom = "text", x = 4, y = 1, label = "R = 0.34", 
           color = 'blue', size = 6)+
  coord_cartesian(xlim = c(0,5), ylim = c(0, 5))+
  ggtitle("Correlation between upregulated \n genes in S1 and S2")+
  theme(panel.background = element_blank())+
  theme_classic()+
  theme(axis.line = element_line(size = 1.2))+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))
ggsave("correlation_UP_S1_vs_S2.png", width = 5, height = 5)

################################################################################
##**Volcano plot for S1, S2 and S3**
################################################################################
#Volcano Plot of S1 (0uM vs 5uM Cu All genes)
volcano_S1 <- data_frame(RNA_Seq_Data_S1.2)
volcano_S1
view(volcano_S1)
typeof(volcano_S1)

volcano_S1.1 <- subset(volcano_S1, select = c(Gene_ID, log2FC, P_value))
view(volcano_S1.1)

Volcano_0vs5 <- EnhancedVolcano(toptable = data.frame(volcano_S1.1),
                                lab = "",
                                x = 'log2FC',
                                y = 'P_value',
                                title = 'S1',
                                xlim = c(-9, +9),
                                ylim = c(0, 150),
                                subtitle = paste0('0µM vs. 15µM Cu(II)'),
                                pCutoff = 0.05,
                                FCcutoff = 0,
                                pointSize = 3.0,
                                labSize = 6.0,
                                col = c('black', 'green', 'blue', 'red'),
                                colAlpha = 0.25,
                                legendLabels = c('Not sig.','Log2FC','p-value',
                                                 'p-value & Log2FC'),
                                legendLabSize = 12,
                                legendIconSize = 5.0,
                                legendPosition = "right",
                                caption = '',
                                cutoffLineType = 'twodash',
                                cutoffLineCol = 'cyan4',
                                cutoffLineWidth = 1,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE)

Volcano_0vs5
Volcano_0vs5 +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))

################################################################################
#Volcano Plot of S2 (0uM vs 15uM Cu All genes)
volcano_S2 <- data_frame(RNA_Seq_Data_S2.2)
view(volcano_S2)
typeof(volcano_S2)

volcano_S2.1 <- subset(volcano_S2, select = c(Gene_ID, log2FC, P_value))
view(volcano_S2.1)

Volcano_0vs15 <- EnhancedVolcano(toptable = data.frame(volcano_S2.1),
                                 lab = "",
                                 x = 'log2FC',
                                 y = 'P_value',
                                 title = 'S2',
                                 xlim = c(-9, +9),
                                 ylim = c(0, 150),
                                 subtitle = paste0('0µM vs. 15µM Cu(II)'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0,
                                 pointSize = 3.0,
                                 labSize = 6.0,
                                 col = c('black', 'green', 'blue', 'red'),
                                 colAlpha = 0.25,
                                 legendLabels = c('Not sig.','Log2FC','p-value',
                                                  'p-value & Log2FC'),
                                 legendLabSize = 12,
                                 legendIconSize = 5.0,
                                 legendPosition = "right",
                                 caption = '',
                                 cutoffLineType = 'twodash',
                                 cutoffLineCol = 'cyan4',
                                 cutoffLineWidth = 1,
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE)

Volcano_0vs15
Volcano_0vs15 +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))

################################################################################
#Volcano Plot of S3 (5uM vs 15uM Cu All genes)
volcano_S3 <- data_frame(RNA_Seq_Data_S3.2)
view(volcano_S3)
typeof(volcano_S3)

volcano_S3.1 <- subset(volcano_S3, select = c(Gene_ID, log2FC, P_value))
view(volcano_S3.1)

Volcano_5vs15 <- EnhancedVolcano(toptable = data.frame(volcano_S3.1),
                                 lab = "",
                                 x = 'log2FC',
                                 y = 'P_value',
                                 title = 'S2',
                                 xlim = c(-9, +9),
                                 ylim = c(0, 150),
                                 subtitle = paste0('0µM vs. 15µM Cu(II)'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0,
                                 pointSize = 3.0,
                                 labSize = 6.0,
                                 col = c('black', 'green', 'blue', 'red'),
                                 colAlpha = 0.25,
                                 legendLabels = c('Not sig.','Log2FC','p-value',
                                                  'p-value & Log2FC'),
                                 legendLabSize = 12,
                                 legendIconSize = 5.0,
                                 legendPosition = "right",
                                 caption = '',
                                 cutoffLineType = 'twodash',
                                 cutoffLineCol = 'cyan4',
                                 cutoffLineWidth = 1,
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE)

Volcano_5vs15
Volcano_5vs15 +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5))

################################################################################
################################################################################

################################################################################
##**Correlation heatmap between common genes in S1 and S2**
################################################################################
view(Common_genes_S1_and_S2)
heatMap_common_S1_and_S2 <- subset(Common_genes_S1_and_S2, 
                                   select = c(Gene_ID, log2FC.x, log2FC.y))
view(heatMap_common_S1_and_S2)

#Change column names of Gene IDs to row names
row.names(heatMap_common_S1_and_S2) <-heatMap_common_S1_and_S2$Gene_ID
view(heatMap_common_S1_and_S2)

#Rename Column names to S1 nad S2
x <- heatMap_common_S1_and_S2 %>% rename(S1 = log2FC.x, S2 = log2FC.y)
view(x)  

#Remove Gene_ID Column
y <- subset(x, select = -c(Gene_ID))
view(y)

#Compute the correlation matrix
cormat <- round(cor(y),2)
head(cormat)

#Create the correlation heatmap with ggplot2
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

################################################################################
######****END****
################################################################################