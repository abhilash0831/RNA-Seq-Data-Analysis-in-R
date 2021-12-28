#Load Libraries.

#Data translocation.
library(writexl)
library(readxl)
library(xlsx)

#Plotting
library(ggplot2)
library(gplots)

#Data wrangling
library(tidyverse)
library(forcats)
library(magrittr)
library(tidyr)
library(dplyr)


################################################################################
##**Read data to draw bubble plot for S1 (0 vs. 5uM Cu(II))**
################################################################################
Bubble_plot_5 <- read_excel("GO_Bubble_plot.xlsx", sheet = "GO_5")
view(Bubble_plot_5)

#Grouping by GO_Category, 
#arranging GO_Term column in ascending order alphabetically
#and reordering that column as is.
Bubble_plot_5 %<>%
  group_by(`GO Category`) %>%
  arrange(GO_Term, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(GO_Term = forcats::fct_reorder(GO_Term, `GO Category`))

#Plotting.
ggplot(Bubble_plot_5, mapping = aes(x = `Gene Count`, 
                                    y = GO_Term, 
                                    size = `Gene Count`)) + 
  geom_tile(mapping = aes(width = Inf, y = GO_Term, fill = `GO Category`), alpha = 0.2) + 
  geom_point(mapping = aes(color = `GO Category`), alpha = 4.0) + 
  scale_fill_manual(values = c("darkcyan", "hotpink4", "darkgoldenrod3")) +
  scale_color_manual(values = c("darkcyan", "hotpink4", "darkgoldenrod3")) +
  theme(panel.background = element_blank()) +
  xlab("Gene Count") +
  ylab("GO Term") +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave("GO5Bubbleplot.png", width = 10, height = 8, dpi = 400, limitsize = FALSE)


################################################################################
##**Read data to draw bubble plot for S2 (0 vs. 15uM Cu(II))**
################################################################################
Bubble_plot_15 <- read_excel("GO_Bubble_plot.xlsx", sheet = "GO_15")
view(Bubble_plot_15)

#Grouping by GO_Category, 
#arranging GO_Term column in ascending order alphabetically
#and reordering that column as is.
Bubble_plot_15 %<>%
  group_by(`GO Category`) %>%
  arrange(GO_Term, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(GO_Term = forcats::fct_reorder(GO_Term, `GO Category`))

#Plotting.
ggplot(Bubble_plot_15, mapping = aes(x = `Gene Count`, 
                                     y = GO_Term, 
                                     size = `Gene Count`)) + 
  geom_tile(mapping = aes(width = Inf, y = GO_Term, fill = `GO Category`), alpha = 0.2) + 
  geom_point(mapping = aes(color = `GO Category`), alpha = 4.0) + 
  scale_fill_manual(values = c("darkcyan", "hotpink4", "darkgoldenrod3")) +
  scale_color_manual(values = c("darkcyan", "hotpink4", "darkgoldenrod3")) +
  theme(panel.background = element_blank()) +
  xlab("Gene Count") +
  ylab("GO Term") +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.y = element_text(angle = 0, size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

################################################################################
######****END****
################################################################################
