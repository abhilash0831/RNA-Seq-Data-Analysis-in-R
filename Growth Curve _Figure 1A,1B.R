#Load Libraries
#Data translocation
library(writexl)
library(readxl)

#Plotting
library(ggplot2)
library(gplots)
library(ggpubr)
library(forcats)
library(magrittr)
library(tidyverse)

#Data wrangling
library(forcats)
library(magrittr)
library(tidyr)
library(dplyr)

################################################################################
###**Read Growth Data from working directory**
################################################################################
Growth_data <- read_xlsx("growth_curve.xlsx")
view(Growth_data)

#Reverse legend order to make graph in the order 0, 5, 15 and not 0, 15, 5
Growth_data$Sample <- factor(Growth_data$Sample, levels = 
                               c("0然","5然", "15然"))

Growth_plot <- ggplot(Growth_data, aes(x=Time, y=OD,
                                         color = Sample, fill = Sample))+
  geom_point(aes(shape = Sample, color = Sample), size = 3)+
  stat_smooth(aes(x=Time, y=OD, linetype = "Cubic Fit"), 
              method = "lm", formula = y ~ poly(x, 3), size = 0.5, se = FALSE, 
              color = "black")+
  #scale_shape_manual(values = c(1, 8, 2)) +
  #scale_color_manual(values = c("red", "blue", "green"))+
  geom_errorbar(data = Growth_data, mapping = aes(x=Time, ymax = SD_Max, ymin = SD_Min),
                width = 10)+
  xlab("\nTime (hours)")+
  ylab("\nO.D (600nm)")+
  theme(panel.background = element_blank())+
  theme(legend.position = "top")+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        axis.title.x = element_blank(),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  facet_wrap(~Sample, scales = "fixed", ncol = 3)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16))
growth_cur <- Growth_plot + guides(linetype = "none")
growth_cur

################################################################################
###**Read protein concentration data from working directory**
################################################################################

Protein_data <- read_xlsx("Protein_estimation.xlsx")
view(Protein_data)

#Reverse legend order to make graph in the order 0, 5, 15 and not 0, 15, 5
Protein_data$Sample <- factor(Protein_data$Sample, levels = 
                                c("0然","5然", "15然"))

Protein_plot <- ggplot(Protein_data, aes(x=Time, y=ProteinConcentration,
                                         color = Sample, fill = Sample))+
  geom_point(aes(shape = Sample, color = Sample), size = 3)+
  stat_smooth(aes(x=Time, y=ProteinConcentration, linetype = "Quadratic Fit"), 
              method = "lm", formula = y ~ poly(x, 2), size = 0.5, se = FALSE,
              color = "black")+
  #scale_shape_manual(values = c(1, 8, 2)) +
  #scale_color_manual(values = c("red", "blue", "green"))+
  geom_errorbar(aes(ymax = SD_Max, ymin = SD_Min), width = 5)+
  xlab("\nTime (hours)")+
  ylab("\nProtein Concentration (mg/L)")+
  theme(panel.background = element_blank())+
  theme(legend.position = "top")+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  facet_wrap(~Sample, scales = "fixed", ncol = 3)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16))
prot_conc <- Protein_plot + guides(linetype = FALSE)
prot_conc

################################################################################
##**Combine growth and protein plot in one frame**
##**Use the function ggarrange() from [ggpubr package]**
################################################################################

Combined_growth_and_protein <- ggarrange(growth_cur, 
                                         prot_conc, labels = c("A", "B"),
                                         ncol = 1, nrow = 2)
Combined_growth_and_protein
