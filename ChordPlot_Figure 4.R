#Data translocation.
library(writexl)
library(readxl)

#Plotting
library(ggplot2)
library(gplots)

library(forcats)
library(magrittr)
library(tidyverse)

#Data wrangling
library(forcats)
library(magrittr)
library(tidyr)
library(dplyr)

#GOPlot
library(GOplot)
library(RColorBrewer)

################################################################################
##Chord plot for MF 15
################################################################################
##Read circ file
##In circ.xlsx logFC column is 'log2FC' and adj_pval is 'pval'

Circ_data <- read_xlsx("circ.xlsx", sheet = "MF_15")
Circ_data_1 <- data.frame(Circ_data)
typeof(Circ_data_1)
str(Circ_data_1)
view(Circ_data_1)

chord <- chord_dat(Circ_data_1)
GOChord(chord, title = 'Chord_plot_MF_S1', lfc.col=c('red','black','cyan'), 
        limit = c(3,5), space = 0.02, gene.order = 'logFC', 
        gene.space = 0.50, gene.size = 5)
GOChord(chord, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), 
        display = 'multiple', labels = 3)

GOBubble(Circ_data_1, title = 'Bubble plot', colour = 
           c('orange', 'darkred', 'gold'), 
         display = 'multiple', labels = 3)

Bar <- GOBar(subset(Circ_data_1, order.by.zscore = FALSE, Category == 'MF'))
Bar+
  geom_text(aes(label = adj_pval, vjust = 1))


#######################################################################
########################################################################
Circ_data123 <- read_xlsx("circ123.xlsx", sheet = "BP_15")
Circ_data_123 <- data.frame(Circ_data123)
typeof(Circ_data_123)
str(Circ_data_123)
view(Circ_data123)

All_15_pval <- read.csv("All_genes_in_15(ONLY p-value).csv", all = FALSE)
View(All_15_pval)

GO_pvalue <- merge(Circ_data123, All_15_pval, by = "Gene_ID", all = FALSE)
View(GO_pvalue)

write.csv(GO_pvalue, "Circ_data_pval15BP.csv", row.names = FALSE)

Circ123 <- read.csv(file = "Circ_data_pval15BP.csv")
view(Circ123)
Bar <- GOBar(subset(Circ123, order.by.zscore = FALSE, Category == 'BP'))
Bar+
  geom_text(aes(label = adj_pval, vjust = 1))


##ribbon.col = brewer.pal (6, "YlGnBu")
install.packages('GOplot')
DF1 <-read_xlsx("DF1.xlsx")
View(DF1)
DF2 <-read_xlsx("DF2.xlsx") 
view(DF2)
circ <- circle_dat(DF2, DF1)

########################################################################################
########################################################################################
##Sample Data from 
##https://wencke.github.io/#display-of-the-relationship-between-genes-and-terms-gochord
devtools::install_github("wencke/wencke.github.io")
library(GOplot)
library(RColorBrewer)

#######################################################################################
#######################################################################################
##GOBar plot for CC -- Creating the circ object fo 15CC

CC_15 <- read.csv("All_CC_15.csv", all = FALSE)
View(CC_15)

All_15_pval <- read.csv("All_genes_in_15(ONLY p-value).csv", all = FALSE)
View(All_15_pval)

CC_15_ALL <- merge(CC_15, All_15_pval, by = "Gene_ID", all = FALSE)
View(CC_15_ALL)

Transport_CC_15 <- GO_CC_15 %>% filter(grepl('Transport', GOCC))
View(Transport_CC_15)

write.csv(Transport_CC_15, "Transport_CC_15.csv", row.names = FALSE)

##Select specific GO terms from the last file created

GO_CC_15 <- read.csv("Circ_CC_15logFC.csv")
view(GO_CC_15)

ABCtrans_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0043190', GOCC))
View(ABCtrans_CC_15)

BFBB_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0009425', GOCC))
View(BFBB_CC_15)

ICM_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0016021', GOCC))
View(ICM_CC_15)
ICM_CC_15_rem_dup <-  ICM_CC_15 %>% distinct()
view(ICM_CC_15_rem_dup)
write.csv(ICM_CC_15_rem_dup, "ICM_CC_15_rem_dup.csv", row.names = FALSE)

PM_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0005886', GOCC))
View(PM_CC_15)
PM_CC_15_rem_dup <-  PM_CC_15 %>% distinct()
view(PM_CC_15_rem_dup)
write.csv(PM_CC_15_rem_dup, "PM_CC_15_rem_dup.csv", row.names = FALSE)

CY_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0005737', GOCC))
View(CY_CC_15)
CY_CC_15_rem_dup <-  CY_CC_15 %>% distinct()
view(CY_CC_15_rem_dup)
write.csv(CY_CC_15_rem_dup, "CY_CC_15_rem_dup.csv", row.names = FALSE)

EXCR_CC_15 <- GO_CC_15 %>% filter(grepl('GO:0005576', GOCC))
View(EXCR_CC_15)
EXCR_CC_15_rem_dup <-  EXCR_CC_15 %>% distinct()
view(EXCR_CC_15_rem_dup)
write.csv(EXCR_CC_15_rem_dup, "EXCR_CC_15_rem_dup.csv", row.names = FALSE)

Count_neg <- nrow(EXCR_CC_15_rem_dup[EXCR_CC_15_rem_dup$log2FC<0,])
Count_neg
Count_pos <- nrow(EXCR_CC_15_rem_dup[EXCR_CC_15_rem_dup$log2FC>0,])
Count_pos

#######################################################################################
#######################################################################################
##GOBar plot for CC -- Creating the circ object fo 5CC
CC_5 <- read.csv("All_CC_5.csv", all = FALSE)
View(CC_5)

All_5_pval <- read.csv("All_genes_in_5(ONLY p-value).csv", all = FALSE)
View(All_5_pval)

CC_5_ALL <- merge(CC_5, All_5_pval, by = "Gene_ID", all = FALSE, by.y = "Gene_ID", by.x = "Gene_ID")
View(CC_5_ALL)

write.csv(CC_5_ALL, "Circ_CC_5logFC.csv", row.names = FALSE)

GO_CC_5 <- read.csv("Circ_CC_5logFC.csv")
view(GO_CC_5)
Transporter_CC_5 <- GO_CC_5 %>% filter(grepl('Transport', Protein.names))
view(Transporter_CC_5)
write.csv(Transporter_CC_5, "TransporterONLY2_CC_5.csv", row.names = FALSE)


ICM_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0016021', GOCC))
View(ICM_CC_5)
ICM_CC_5_rem_dup <-  ICM_CC_5 %>% distinct()
view(ICM_CC_5_rem_dup)
write.csv(ICM_CC_5_rem_dup, "ICM_CC_5_rem_dup.csv", row.names = FALSE)

PM_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0005886', GOCC))
View(PM_CC_5)
PM_CC_5_rem_dup <-  PM_CC_5 %>% distinct()
view(PM_CC_5_rem_dup)
write.csv(PM_CC_5_rem_dup, "PM_CC_5_rem_dup.csv", row.names = FALSE)

CY_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0005737', GOCC))
View(CY_CC_5)
CY_CC_5_rem_dup <-  CY_CC_5 %>% distinct()
view(CY_CC_5_rem_dup)
write.csv(CY_CC_5_rem_dup, "CY_CC_5_rem_dup.csv", row.names = FALSE)

BFBB_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0009425', GOCC))
View(BFBB_CC_5)
BFBB_CC_5_rem_dup <-  BFBB_CC_5 %>% distinct()
view(BFBB_CC_5_rem_dup)
write.csv(BFBB_CC_5_rem_dup, "BFBB_CC_5_rem_dup.csv", row.names = FALSE)

EXCR_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0005576', GOCC))
View(EXCR_CC_5)
EXCR_CC_5_rem_dup <-  EXCR_CC_5 %>% distinct()
view(EXCR_CC_5_rem_dup)
write.csv(EXCR_CC_5_rem_dup, "EXCR_CC_5_rem_dup.csv", row.names = FALSE)

ABCtrans_CC_5 <- GO_CC_5 %>% filter(grepl('GO:0043190', GOCC))
View(ABCtrans_CC_5)
ABCtrans_CC_5_rem_dup <-  ABCtrans_CC_5 %>% distinct()
view(ABCtrans_CC_5_rem_dup)
write.csv(ABCtrans_CC_5_rem_dup, "ABCtrans_CC_5_rem_dup.csv", row.names = FALSE)

Count_neg <- nrow(ABCtrans_CC_5_rem_dup[ABCtrans_CC_5_rem_dup$log2FC<0,])
Count_neg
Count_pos <- nrow(ABCtrans_CC_5_rem_dup[ABCtrans_CC_5_rem_dup$log2FC>0,])
Count_pos

###################################################################################
###################################################################################
##BarPlot for CC (5 and 15) on pvalue and z-score
Circ_CC <- read_xlsx("Circ_CC.xlsx", sheet = "Circ_CC_5")
view(Circ_CC)
Bar <- GOBar(subset(Circ_CC, order.by.zscore = FALSE, category == 'CC'))
Bar+
  geom_text(aes(label = adj_pval, vjust = 1))+
  geom_hline(yintercept=1.30)







