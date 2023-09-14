#new file
#Experimental evaluation of environmental DNA detection of a rare fish in turbid water
#Running title: eDNA experiments in turbid water
#Ann E. Holmes, Melinda R. Baerwald, Jeff Rodzen, Brian M. Schreier, Brian Mahardja, and Amanda J. Finger
#corresponding author aholmes@ucdavis.edu

#R code for box plots 

setwd("~/github/ann_holmes/DSM_eDNA")
library(ggplot2)
library(ggpattern) #to add a pattern to the box plots
library(tidyverse)
library(cowplot) #for saving plot
#Running R version 4.1.3 in RStudio 1.4.1717

#read in qPCR results
#Ct values for non-detects have been imputed using non-detects R package (McCall et al. 2014)
#DNA copy numbers were generated using the standard curve following standardized eDNA methodology (Klymus et al. 2020; Merkes et al. 2020)
#binning the number of replicates that amplified (out of 6)

filter_data <-read.csv("S5_eDNA_expt_qPCR_results.csv")

filter_data$Filter_type_f = factor(filter_data$Filter_type, levels=c('ST','PC5','PC10','GF'))
filter_data$Tank_water = as.factor(filter_data$Tank_water)
filter_data$Treatment = as.factor(filter_data$Treatment)
filter_data$Det_rate = as.factor(filter_data$Det_rate)
filter_data$Detections = factor(filter_data$Det_rate_bin, levels=c('6 reps','4-5 reps','≤3 reps'))

#create boxplots to show results in DNA copies detected
#x axis for each box is the amount of tank water added Tank_water
#DNA copies per liter of water filtered (calculated from: DNA copies per reaction/DNA extract per reaction (6.1 ul) * Volume of DNA extract (100 ul) / Volume of water filtered (0.1 to 1 L depending on the sample)) is the response variable
#plotted by Filter_type and Treatment: 5 NTU (non-turbid), 50 NTU (turbid without prefilter), and 50 NTU - PF (turbid with prefilter)

#color for non-turbid water: blue #bbe4ee (matches point in map figure)
#color for turbid water: brown #c1b090  (matches point in map figure)
#color for turbid water with prefilter: yellow #fefd98

#retain only rows that are detections (replicates that amplified with a Ct value less than the LOD)
det_data <- filter(filter_data, Imputed == "No")

#boxplots for copies detected per L of water filtered
#LOD (red line) and LOQ (dash black line) have been calculated for copies/L so their value changes depending on the volume of water filtered for the sample
copies_box <- ggplot(det_data, 
                     aes(x = Tank_water, 
                         y = Copies_per_L_water, 
                         group=Comb, 
                         pattern=Detections)) +  
  geom_rect(data=det_data, 
            aes(fill = Treatment), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  facet_wrap(~Treatment~Filter_type_f, scales = "free") +  
  geom_boxplot() +
  geom_boxplot_pattern(color = "black", 
                       pattern_fill = "white", 
                       pattern_angle = 45, 
                       pattern_density = 0.1, 
                       pattern_spacing = 0.025, 
                       pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c('none', 'stripe', 'crosshatch')) +
  scale_fill_manual(values = alpha(c("#bbe4ee", #light blue for nonturbid
                                     "#c1b090", #brown for turbid
                                     "#fefd98" ))) + #yellow for prefiltered turbid
  theme_bw() +
  labs(x = "Tank water (mL) added to estuary water (1L)", y = "Delta smelt DNA copies per L") +
  geom_hline(data=filter_data_bins, #LOQ line
             aes(yintercept=eff_LOQ_per_L_water), 
             linetype="dashed", 
             color = "black") + 
  geom_hline(data=filter_data_bins, #LOD line
             aes(yintercept=eff_6repLOD_per_L_water),
             color = "red") +
  theme(axis.text.y = element_text(size=13), 
        axis.text.x = element_text(size=13, angle = 45, vjust = 1, hjust = 1), 
        axis.title.x = element_text(size=15, face="bold"), 
        axis.title.y = element_text(size=15, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "bottom") +
  coord_flip()

save_plot("Figure_4_qPCR_results.png",
          copies_box,
          dpi=300,
          base_height = 10,
          base_width = 8,
          bg="white")

#same boxplots but response variable is Ct value
#note that LOD (red line) is now on the right and LOQ (dash black line) since larger Ct values indicate fewer copies detected
Ct_box <- ggplot(det_data, aes(x = Tank_water, y = Ct_raw_plot, group=Comb, pattern=Detections)) +  
  scale_y_continuous(limits=c(32,41),expand = c(0, 0)) +
  geom_rect(data=det_data, aes(fill = Treatment), xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
  facet_wrap(~Treatment~Filter_type_f) +  
  geom_boxplot() +
  geom_boxplot_pattern(color = "black", 
                       pattern_fill = "white", 
                       pattern_angle = 45, 
                       pattern_density = 0.1, 
                       pattern_spacing = 0.025, 
                       pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c('none', 'stripe', 'crosshatch')) +
  scale_fill_manual(values = alpha(c("#bbe4ee", #light blue for nonturbid
                                     "#c1b090", #brown for turbid
                                     "#fefd98" ))) + #yellow for prefiltered turbid
  theme_bw() +
  labs(x = "Tank water (mL) added to estuary water (1L)", y = "Ct value") +
  geom_hline(yintercept=35.445, linetype="dashed", color = "black") + #LOQ
  geom_hline(yintercept=40.378, color = "red") + #LOD
  theme(axis.text.x = element_text(size=13), 
        axis.text.y = element_text(size=13,), 
        axis.title.x = element_text(size=13, face="bold"), 
        axis.title.y = element_text(size=13, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "bottom") +
  coord_flip()

Ct_box