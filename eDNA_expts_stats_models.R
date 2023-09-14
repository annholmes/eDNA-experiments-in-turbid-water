#Experimental evaluation of environmental DNA detection of a rare fish in turbid water
#Running title: eDNA experiments in turbid water
#Ann E. Holmes, Melinda R. Baerwald, Jeff Rodzen, Brian M. Schreier, Brian Mahardja, and Amanda J. Finger
#corresponding author aholmes@ucdavis.edu

#R code for stats and models comparing detection by filtration methods and turbidity

setwd("~/github/ann_holmes/DSM_eDNA")

library(dplyr)
library(ggplot2)
library(FSA) #dunn test
library(DHARMa) #test for overdispersion
library(lme4) #glmm
library(MASS) #fit negative binomial glm
library(bbmle) #AICc comparison
library(openxlsx) #save model summary to excel

#need to reconcile file names

#entire data set with results for each technical replicate
filter_data <-read.csv("S5_eDNA_expt_qPCR_results.csv")

filter_data$Filter_type_f = factor(filter_data$Filter_type, levels=c('GF','PC5','PC10','ST'))
filter_data$Tank_water = as.factor(filter_data$Tank_water)
filter_data$Treatment = as.factor(filter_data$Treatment)
filter_data$Det_rate = as.factor(filter_data$Det_rate)
filter_data$Detections = factor(filter_data$Det_rate_bin, levels=c('6 reps','4-5 reps','≤3 reps'))

#check for normality Shapiro-Wilk normality test
#hypothesis is that data is normally distributed
shapiro.test(filter_data$Copies_per_L_water)
#significant p value indicates a rejected hypothesis hypothesis; data is not normally distributed

#considering non-turbid water as the control level
nonturb <- filter(filter_data, Turbidity == "5")
shapiro.test(nonturb$Copies_per_L_water)
#checking for normality just in case; also not normal 

#Kruskal-Wallis rank sum test to determine if copies captured in nonturbid water differ by filter type
#hypothesis is that samples are not different
kruskal.test(Copies_per_L_water ~ Filter_type_f, data = nonturb)
kruskal.test(Copies_per_L_water ~ Filter_type_f, data = filter_data)
#significant p value indicates a rejected hypothesis; samples collected by different filters are different

dunnTest(Copies_per_L_water ~ Filter_type_f,
         data=nonturb,
         method="bonferroni")

rem_out_nonturb <-nonturb %>% filter(Copies_per_L_water < 12000)

kruskal.test(Copies_per_L_water ~ Filter_type_f, data = rem_out_nonturb)

dunnTest(Copies_per_L_water ~ Filter_type_f,
         data=rem_out_nonturb,
         method="bonferroni")

#results by sample
data_by_sample <-read.csv("S5_qPCR_results_imputed_by_sample.csv")

data_by_sample$Filter_type_f = factor(data_by_sample$Filter_type, levels=c('ST','PC5','PC10','GF'))
data_by_sample$Tank_water = as.factor(data_by_sample$Tank_water)
data_by_sample$Treatment = as.factor(data_by_sample$Treatment)
data_by_sample$Det_rate = as.factor(data_by_sample$Det_rate)

#check for normality Shapiro-Wilk normality test
#hypothesis is that data is normally distributed
shapiro.test(data_by_sample$Mean_copies_per_L_by_sample)
#significant p value indicates a rejected hypothesis hypothesis; data is not normally distributed

#considering non-turbid water as the control level
nonturb_by_sample <- filter(data_by_sample, Turbidity == "5")
shapiro.test(nonturb_by_sample$Mean_copies_per_L_by_sample)
#checking for normality just in case; also not normal 

#Kruskal-Wallis rank sum test to determine if copies captured in nonturbid water differ by filter type
#null hypothesis is that filter types do not generate different values for copies/L
kruskal.test(Mean_copies_per_L_by_sample ~ Filter_type_f, data = nonturb_by_sample)
#significant p value indicates a rejected hypothesis; samples collected by different filters are different

dunnTest(Mean_copies_per_L_by_sample ~ Filter_type_f,
         data=nonturb_by_sample,
         method="bonferroni")

#null hypothesis is that treatments (5 NTU, 50 NTU, 50 NTU+) do not generate different values for copies/L
kruskal.test(Mean_copies_per_L_by_sample ~ Treatment, data = data_by_sample)

dunnTest(Mean_copies_per_L_by_sample ~ Treatment,
         data=data_by_sample,
         method="bonferroni")

GF <- filter(data_by_sample, Filter_type_f == "GF")
kruskal.test(Mean_copies_per_L_by_sample ~ Treatment, data = GF)
dunnTest(Mean_copies_per_L_by_sample ~ Treatment,
         data=GF,
         method="bonferroni")

ST <- filter(data_by_sample, Filter_type_f == "ST")
kruskal.test(Mean_copies_per_L_by_sample ~ Treatment, data = ST)
dunnTest(Mean_copies_per_L_by_sample ~ Treatment,
         data=ST,
         method="bonferroni")

PC5 <- filter(data_by_sample, Filter_type_f == "PC5")
kruskal.test(Mean_copies_per_L_by_sample ~ Treatment, data = PC5)
dunnTest(Mean_copies_per_L_by_sample ~ Treatment,
         data=PC5,
         method="bonferroni")

PC10 <- filter(data_by_sample, Filter_type_f == "PC10")
kruskal.test(Mean_copies_per_L_by_sample ~ Treatment, data = PC10)
dunnTest(Mean_copies_per_L_by_sample ~ Treatment,
         data=PC10,
         method="bonferroni")

#testing for a statistical difference in detection rate by sample
#a sample is considered positive if at least half of the replicates produce a Ct value below the Ct value corresponding to the LOD
#samples with zero, 1 or 2 replicates amplifying are considered negative

#for comparisons where 1 of 6 bottles is negative
five <- data.frame(
  "detection_yes" = c(6, 0),
  "detection_no" = c(5, 1),
  row.names = c("nonturbid", "turbid"),
  stringsAsFactors = FALSE)

fisher.test(five)
#p-value = 1

#for comparisons where 3 of 6 bottles are negative
three <- data.frame(
  "detection_yes" = c(6, 0),
  "detection_no" = c(3, 3),
  row.names = c("nonturbid", "turbid"),
  stringsAsFactors = FALSE)

fisher.test(three)
#p-value = 0.1818

#consider the same data, but set the threshold for positive to 4 of 6 replicates (instead of 3 of 6 replicates)
#now there are cases where only 1 bottle of 6 meets the threshold for positive
one <- data.frame(
  "detection_yes" = c(6, 0),
  "detection_no" = c(1, 5),
  row.names = c("nonturbid", "turbid"),
  stringsAsFactors = FALSE)

fisher.test(one)
#p-value = 0.01515


#glm to examine which variables influence eDNA copy number (average copies by biological replicate with imputed values for nondetects)

eDNA_det <- read.csv("S5_eDNA_expt_qPCR_results_det.csv")
eDNA_det$Filter_type_f = factor(eDNA_det$Filter_type, levels=c('GF','PC5','PC10','ST'))
eDNA_det$Tank_water = as.factor(eDNA_det$Tank_water)
eDNA_det$Treatment = as.factor(eDNA_det$Treatment)
eDNA_det$Det_rate = as.factor(eDNA_det$Det_rate)

mod_poisson <-glm(Copies_per_L_water ~ Filter_type*Treatment, 
                  data=eDNA_det, family = poisson)

sim_poisson <- simulateResiduals(mod_poisson, refit=T) 
testDispersion(sim_poisson)
#dispersion test significant

#make plots showing results of test
png("sim_poisson_plot.png",
    res = 300,
    width = 20,
    height = 10,
    units = "cm")

plot(sim_poisson)

dev.off()

#now using negative binomial distribution
mod_nb <-glm.nb(Copies_per_L_water ~ Filter_type*Treatment, 
                data=eDNA_det)

sim_nb <- simulateResiduals(mod_nb, refit=T) 
testDispersion(sim_nb)
#dispersion test significant

png("sim_nb_plot.png",
    res = 300,
    width = 20,
    height = 10,
    units = "cm")

sim_nb_plot <- plot(sim_nb) 

dev.off()

summary(mod_nb)

#glm for copy number/L by sample (biological replicate)
#negative binomial distribution

#full model
full_mod_av_copies=glm.nb(Copies_per_L_water ~ 
                Filter_type_f + Prefilter + Turbidity + Tank_water +
                Filter_type_f:Prefilter + Filter_type_f:Turbidity + Filter_type_f:Tank_water + Prefilter:Tank_water + Turbidity:Tank_water, 
              data=eDNA_det)

summary(full_mod_av_copies)

#write full model summary results to excel file
sum_full_mod_av_copies <- coef(summary(full_mod_av_copies))
full_mod_summary<- createWorkbook()
addWorksheet(full_mod_summary, "full_mod_av_copies")
writeData(full_mod_summary, "full_mod_av_copies", 
          matrix(sum_full_mod_av_copies, ncol = 4),
          colNames = FALSE)
saveWorkbook(full_mod_summary, file = "full_mod_av_copies.xlsx")

#best model
best_mod_av_copies=glm.nb(Copies_per_L_water ~ Filter_type_f + Turbidity +
                          Filter_type_f:Prefilter,
                          data=eDNA_det)

summary(best_mod_av_copies)

#AICc comparison in bbmle
AICctab(full_mod_av_copies, best_mod_av_copies,
        weights = TRUE, base = TRUE)

#write best model summary results to excel file
sum_best_mod_av_copies <- coef(summary(best_mod_av_copies))
best_mod_summary<- createWorkbook()
addWorksheet(best_mod_summary, "best_mod_av_copies")
writeData(best_mod_summary, "best_mod_av_copies", 
          matrix(sum_best_mod_av_copies, ncol = 4),
          colNames = FALSE)
saveWorkbook(best_mod_summary, file = "best_mod_av_copies.xlsx")

#glm for presence/absence at sample level using 3 reps as the threshold for presence
#binomial distribution

mod_det_3reps=glm(Det_3reps ~ 
             Filter_type_f + Prefilter + Turbidity + Filtr_vol_L + Tank_water, 
           data=eDNA_det, 
           family = binomial)
summary(mod_det_3reps)
#no significant effects

#glm for presence/absence at sample level using 4 reps as the threshold for presence
mod_det_4reps=glm(Det_4reps ~ 
                    Filter_type_f + Prefilter + Turbidity + Filtr_vol_L + Tank_water, 
                  data=eDNA_det, 
                  family = binomial)
summary(mod_det_4reps)
#no significant effects

#mixed effects model for presence/absence using qPCR replicate presence/absence 
#binomial distribution

mod_rep_det_full=glmer(Detection ~ Filter_type_f + Prefilter + Turbidity + Tank_water + 
             (1|Treatment/Filter_type_f/Bio_rep), 
           data=filter_data, 
           family = binomial)

summary(mod_rep_det_full)

sum_mod_rep_det_full <- coef(summary(mod_rep_det_full))
mod_rep_det_full_summary<- createWorkbook()
addWorksheet(mod_rep_det_full_summary, "mod_rep_det_full")
writeData(mod_rep_det_full_summary, "mod_rep_det_full", 
          matrix(sum_mod_rep_det_full, ncol = 4),
          colNames = FALSE)
saveWorkbook(mod_rep_det_full_summary, file = "mod_rep_det_full.xlsx")

#questions: 
#no fixed effects are significant when I run the model with filtration volume, but when I take it out turbidity is significant
#filtration vol is too correlated with filter type in the copies model, and therefore not included in this model, but it does run with filtration vol in the model
#interactions cause the model to not run (maybe because 5 NTU and many 50NTU/PF are 100% detection?)

mod_rep_det_best=glmer(Detection ~ Turbidity + Prefilter + 
                         (1|Treatment/Filter_type_f/Bio_rep), 
                       data=filter_data, 
                       family = binomial)

summary(mod_rep_det_best)
#turbidity is significant but prefilter is not, however there is a slight (0.4) improvement in AICc when I add prefilter
mod_rep_det_best1=glmer(Detection ~ Turbidity + 
                         (1|Treatment/Filter_type_f/Bio_rep), 
                       data=filter_data, 
                       family = binomial)
summary(mod_rep_det_best1)
#in this model the p-value of turbidity goes to 0.017 from 0.003

sum_mod_rep_det_best <- coef(summary(mod_rep_det_best))
mod_rep_det_best_summary<- createWorkbook()
addWorksheet(mod_rep_det_best_summary, "mod_rep_det_best")
writeData(mod_rep_det_best_summary, "mod_rep_det_best", 
          matrix(sum_mod_rep_det_best, ncol = 4),
          colNames = FALSE)
saveWorkbook(mod_rep_det_best_summary, file = "mod_rep_det_best.xlsx")

AICctab(mod_rep_det_full, mod_rep_det_best, mod_rep_det_best1,
        weights = TRUE, base = TRUE)
