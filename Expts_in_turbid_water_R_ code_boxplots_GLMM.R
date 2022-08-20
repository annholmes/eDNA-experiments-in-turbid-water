#code for box plots and GLMM

setwd("~/github/ann_holmes/DSM_eDNA")

filter_data <-read.csv("S3. Filter sample metadata and qPCR results .csv")

library(ggplot2)

#color for non-turbid water: blue #bbe4ee (matches Fig 2)
#color for turbid water: brown #c1b090  (matches Fig 2)
#color for turbid water with prefilter: yellow #fefd98

#x axis is the amount of tank water added Tank_water
#y axis is the DNA concentration given in Log(copies+1)
#plotted by Filter_type and Treatment: 5 NTU (non-turbid), 50 NTU (turbid without prefilter), and 50 NTU - PF (turbid with prefilter)

#order filter type
filter_data$Filter_type_f = factor(filter_data$Filter_type, levels=c('ST','PC5','PC10','GF'))

#plot with DNA copies (represented as Log(Copies+1)) as response variable
ggplot(filter_data, aes(x = Tank_water, y = Log_copies_plus_1, group=Tank_water)) +   
  geom_rect(data=filter_data, aes(fill = Treatment), 
            xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
  facet_wrap(~Treatment~Filter_type_f) +  
  geom_boxplot() +
  scale_fill_manual(values = alpha(c("#bbe4ee", "#c1b090", "#fefd98"))) +
  geom_hline(yintercept=1.8325, linetype="dashed", color = "black") + #Log of LOQ+1
  geom_hline(yintercept=0.3053, color = "red") + #Log of LOD+1
  labs(x = "Tank water added to estuary water (mL)", y = "Log(copies + 1)") + 
  theme_bw()
#saved as jpeg 1000x1000

#plot with detection rate (represented as % positive replicates out of 6 for each of 3 bottles) as response variable
ggplot(filter_data, aes(x = Tank_water, y = Det_rate, group=Tank_water)) +   
  geom_rect(data=filter_data, aes(fill = Treatment), 
            xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.3) +
  facet_wrap(~Treatment~Filter_type_f) +  
  geom_boxplot() +
  scale_fill_manual(values = alpha(c("#bbe4ee", "#c1b090", "#fefd98"))) +
  labs(x = "Tank water added to estuary water (mL)", y = "Detection rate by biological replicate") + 
  theme_bw()
#saved as jpeg 1000x1000

#GLMM
library(lme4)
library(bbmle)

head(filter_data)

filter_data$Biological_rep <-factor(filter_data$Biological_rep)
filter_data$Filter_type <-factor(filter_data$Filter_type)
filter_data$Prefilter <-factor(filter_data$Prefilter)
filter_data$Treatment <-factor(filter_data$Treatment)

#data exploration
mod1 <- lm(Log_copies_plus_1 ~ Filter_type, data = filter_data)

summary(mod1)
plot(resid(mod1) ~ filter_data$Turbidity)
plot(resid(mod1) ~ filter_data$Prefilter)

hist(filter_data$Log_copies_plus_1) #histogram of the response variable- skewed because of zeros
dotchart(filter_data$Log_copies_plus_1)

plot(Log_copies_plus_1~Filter_type, data = filter_data) #plot response variable by each of the predictors
plot(Log_copies_plus_1~Prefilter, data = filter_data)
plot(Log_copies_plus_1~Turbidity, data = filter_data)
plot(Log_copies_plus_1~Filtr_vol, data = filter_data)
plot(Log_copies_plus_1~Tank_water, data = filter_data)
plot(Log_copies_plus_1~Treatment, data = filter_data)

#Log(eDNAcopies+1) is response variable
#linear model no random effect
mod.full1 <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#linear model with Biological replicate (each 1 L bottle) as random effect
mod.random1 <- lmer(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter + (1|Biological_rep), data=filter_data)
#rank deficiency warning is a result of prefilter only applied to turbid samples
#model result still valid

summary(mod.full1)
summary(mod.random1)

AICctab(mod.full1, mod.random1, weights = TRUE)

#result of AICc
#mod.full1 0.0 with 14 df, model weight 1
#mod.random1 dAICc  68.3 with 15 df, model weight <0.001
#select model without random effect

mod1a <- lm(Log_copies_plus_1 ~ 1, data=filter_data) 

#removing interactions
mod1b <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water, data=filter_data)
mod1c <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity, data=filter_data)
mod1d <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Prefilter, data=filter_data)

#remove filtration volume
mod1e <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Tank_water, data=filter_data)
mod1f <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Turbidity, data=filter_data)
mod1g <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Prefilter, data=filter_data)
mod1h <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#remove tank water added
mod1i <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol, data=filter_data)
mod1j <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Turbidity, data=filter_data)
mod1k <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Prefilter, data=filter_data)
mod1l <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#remove turbidity
mod1m <- lm(Log_copies_plus_1 ~ Filter_type + Prefilter + Filtr_vol + Tank_water, data=filter_data)
mod1n <- lm(Log_copies_plus_1 ~ Filter_type + Prefilter + Filtr_vol + Tank_water+ Filter_type*Prefilter, data=filter_data)

#remove filter type
mod1o <- lm(Log_copies_plus_1 ~ Turbidity + Prefilter + Filtr_vol + Tank_water, data=filter_data)

#remove prefilter
mod1p <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Filtr_vol + Tank_water, data=filter_data)
mod1q <- lm(Log_copies_plus_1 ~ Turbidity + Filter_type + Filtr_vol + Tank_water+ Filter_type*Turbidity, data=filter_data)

AICctab(mod.full1, mod1a, mod1b, mod1c, mod1d, mod1e, mod1f, mod1g, mod1i, mod1h, mod1j, mod1k, mod1l, mod1m, mod1n, mod1o, mod1p, mod1q, weights = TRUE)

#with detection as response variable
#linear model no random effect
mod.full2 <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#linear model with Biological replicate (each 1 L bottle) as random effect
mod.random2 <- lmer(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter + (1|Biological_rep), data=filter_data)
#rank deficiency warning is a result of prefilter only applied to turbid samples
#model result still valid

summary(mod.full2)
summary(mod.random2)

AICctab(mod.full2, mod.random2, weights = TRUE)

#result of dAICc
#mod.full1 0.0 with 14 df, model weight 1
#mod.random1 dAICc 91 with 15 df, model weight <0.001
#select model without random effect

mod2a <- lm(Detection ~ 1, data=filter_data) 

#removing interactions
mod2b <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water, data=filter_data)
mod2c <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Turbidity, data=filter_data)
mod2d <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Tank_water + Filter_type*Prefilter, data=filter_data)

#remove filtration volume
mod2e <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Tank_water, data=filter_data)
mod2f <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Turbidity, data=filter_data)
mod2g <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Prefilter, data=filter_data)
mod2h <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Tank_water + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#remove tank water added
mod2i <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol, data=filter_data)
mod2j <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Turbidity, data=filter_data)
mod2k <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Prefilter, data=filter_data)
mod2l <- lm(Detection ~ Turbidity + Filter_type + Prefilter + Filtr_vol + Filter_type*Turbidity + Filter_type*Prefilter, data=filter_data)

#remove turbidity
mod2m <- lm(Detection ~ Filter_type + Prefilter + Filtr_vol + Tank_water, data=filter_data)
mod2n <- lm(Detection ~ Filter_type + Prefilter + Filtr_vol + Tank_water+ Filter_type*Prefilter, data=filter_data)

#remove filter type
mod2o <- lm(Detection ~ Turbidity + Prefilter + Filtr_vol + Tank_water, data=filter_data)

#remove prefilter
mod2p <- lm(Detection ~ Turbidity + Filter_type + Filtr_vol + Tank_water, data=filter_data)
mod2q <- lm(Detection ~ Turbidity + Filter_type + Filtr_vol + Tank_water+ Filter_type*Turbidity, data=filter_data)

AICctab(mod.full2, mod2a, mod2b, mod2c, mod2d, mod2e, mod2f, mod2g, mod2i, mod2h, mod2j, mod2k, mod2l, mod2m, mod2n, mod2o, mod2p, mod2q, weights = TRUE)
        
###########

##Assemble Full Starting (reference) Model
mod.full <- lme(Log_copies ~ Turbidity + Filter + Prefilter + Volume + Tank + Turbidity*Filter + Turbidity*Prefilter + Turbidity*Volume + Turbidity*Tank + Filter*Prefilter + Filter*Volume + Filter*Tank + Prefilter*Volume + Prefilter*Tank + Volume*Tank, 
                method = "REML",
                data = experiment)
summary(mod.full) #look at summary output



## 2) DETERMINE THE PROPER RANDOM STRUCTURE (USING method="REML")
mod.full <- lme(Log_copies ~ Turbidity + Filter + Prefilter + Volume + Tank + Turbidity*Filter + Turbidity*Prefilter + Turbidity*Volume + Turbidity*Tank + Filter*Prefilter + Filter*Volume + Filter*Tank + Prefilter*Volume + Prefilter*Tank + Volume*Tank, 
                method = "REML",
                data = experiment)
summary(mod.full) #look at summary output

#need to compare with gls() b/c no random effect in this model
mod.date <- gls(Log_copies ~ Turbidity + Filter + Prefilter + Volume + Tank + Turbidity*Filter + Turbidity*Prefilter + Turbidity*Volume + Turbidity*Tank + Filter*Prefilter + Filter*Volume + Filter*Tank + Prefilter*Volume + Prefilter*Tank + Volume*Tank, 
                method = "REML", #drop random effect
                data = experiment)
summary(mod.date)
anova(mod.full, mod.date) 


#4 Jan 22 analysis starts here

#### INFORMATION ABOUT EXPERIMENTAL DESIGN & PARAMETERS MEASURED
#response variable (s)- Copies, Detection
#fixed effects (categorical): 
#Filter (Filter type used for filtration: material, pore size and shape)
#Turbidity (Estuarine water at 5 NTU or 50 NTU)
#Tank (Amount of water from tank of delta smelt added to container)
#fixed effects (continuous): 
#Volume (Total volume of water filtered before filter clogged)
#random effect:
#Bottle (individual containers i.e. biological replicates) 

turbidity<- read.csv("Filtration_turbidity.csv")
turbidity$Bottle <-factor(turbidity$Bottle)
turbidity$Filter <-factor(turbidity$Filter)

mod0 <- lmer(Copies ~ 1  + (1|Bottle), data=turbidity)   
mod1 <- lmer(Copies ~ Turbidity + Filter + Tank + Volume + (1|Bottle), data=turbidity) 
mod2 <- lmer(Copies ~ Turbidity + Filter + Tank + (1|Bottle), data=turbidity) 
mod3 <- lmer(Copies ~ Turbidity + Filter + Volume + (1|Bottle), data=turbidity) 
mod4 <- lmer(Copies ~ Turbidity + Tank + Volume + (1|Bottle), data=turbidity) 
mod5 <- lmer(Copies ~ Filter + Tank + Volume + (1|Bottle), data=turbidity) 


mod6 <- lmer(Copies ~ Turbidity  + Filter + (1|Bottle), data=turbidity) #AICc 3118.2
mod7 <- lmer(Copies ~ Turbidity  + Volume + (1|Bottle), data=turbidity) 
mod8 <- lmer(Copies ~ Turbidity  + Tank + (1|Bottle), data=turbidity) 
mod9 <- lmer(Copies ~ Filter  + Tank + (1|Bottle), data=turbidity) 
mod10 <- lmer(Copies ~ Filter  + Volume + (1|Bottle), data=turbidity) 
mod11 <- lmer(Copies ~ Volume  + Tank + (1|Bottle), data=turbidity) 

mod12 <- lmer(Copies ~ Turbidity  + (1|Bottle), data=turbidity) 
mod13 <- lmer(Copies ~ Tank + (1|Bottle), data=turbidity) 
mod14 <- lmer(Copies ~ Filter + (1|Bottle), data=turbidity) 
mod15 <- lmer(Copies ~ Volume + (1|Bottle), data=turbidity) 

mod16 <- lmer(Copies ~ Turbidity  + Filter + Turbidity*Filter + (1|Bottle), data=turbidity)

mod17 <- lmer(Copies ~ Turbidity  + Filter + Turbidity*Filter + Volume + Tank + (1|Bottle), data=turbidity)
mod18 <- lmer(Copies ~ Turbidity  + Filter + Turbidity*Filter + Tank + (1|Bottle), data=turbidity)
mod19 <- lmer(Copies ~ Turbidity  + Filter + Turbidity*Filter + Volume + (1|Bottle), data=turbidity)

summary(mod0)
summary(mod1)
summary(mod2)
AICctab(mod0, mod1, mod2, mod3, mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16, mod17,mod18, mod19, weights=TRUE, base=TRUE)

#500ml added
turbidity_500<- read.csv("Filtration_turbidity_500.csv")
turbidity_500$Bottle <-factor(turbidity_500$Bottle)
turbidity_500$Filter <-factor(turbidity_500$Filter)

mod0a <- lmer(Copies ~ 1  + (1|Bottle), data=turbidity_500)   
mod1a <- lmer(Copies ~ Turbidity + Filter + Volume + (1|Bottle), data=turbidity_500) 

mod2a <- lmer(Copies ~ Turbidity  + Filter + (1|Bottle), data=turbidity_500) #AICc 3118.2
mod3a <- lmer(Copies ~ Turbidity  + Volume + (1|Bottle), data=turbidity_500) 
mod4a <- lmer(Copies ~ Filter  + Volume + (1|Bottle), data=turbidity_500) 

mod5a <- lmer(Copies ~ Turbidity  + (1|Bottle), data=turbidity_500) 
mod6a <- lmer(Copies ~ Filter + (1|Bottle), data=turbidity_500) 
mod7a <- lmer(Copies ~ Volume + (1|Bottle), data=turbidity_500) 

mod8a <- lmer(Copies ~ Turbidity  + Filter + Turbidity*Filter + (1|Bottle), data=turbidity_500)
mod9a <- lmer(Copies ~ Turbidity  + Filter + Volume + Turbidity*Filter + (1|Bottle), data=turbidity_500)
mod10a <- lmer(Copies ~ Turbidity  + Filter + Volume + Volume*Filter + (1|Bottle), data=turbidity_500)
mod11a <- lmer(Copies ~ Turbidity  + Filter + Volume + Volume*Filter + Turbidity*Filter + (1|Bottle), data=turbidity_500)
mod12a <- lmer(Copies ~ Turbidity  + Filter + Volume + Volume*Filter + Turbidity*Filter + Turbidity*Volume + (1|Bottle), data=turbidity_500)

AICctab(mod0a, mod1a, mod2a, mod3a, mod4a,mod5a,mod6a,mod7a, mod8a, mod9a, mod10a,mod11a, mod12a, weights=TRUE, base=TRUE)


#prefilter
prefilter <- read.csv("Filtration_prefilter.csv")
prefilter$Bottle <-factor(prefilter$Bottle)
prefilter$Filter <-factor(prefilter$Filter)

mod0b <- lmer(Copies ~ 1  + (1|Bottle), data=prefilter)   
mod1b <- lmer(Copies ~ Prefilter + Filter + Tank + Volume + (1|Bottle), data=prefilter) 
mod2b <- lmer(Copies ~ Prefilter + Filter + Tank + (1|Bottle), data=prefilter) 
mod3b <- lmer(Copies ~ Prefilter + Filter + Volume + (1|Bottle), data=prefilter) 
mod4b <- lmer(Copies ~ Prefilter + Tank + Volume + (1|Bottle), data=prefilter) 
mod5b <- lmer(Copies ~ Filter + Tank + Volume + (1|Bottle), data=prefilter) 


mod6b <- lmer(Copies ~ Prefilter + Filter + (1|Bottle), data=prefilter) #AICc 3118.2
mod7b <- lmer(Copies ~ Prefilter + Volume + (1|Bottle), data=prefilter) 
mod8b <- lmer(Copies ~ Prefilter + Tank + (1|Bottle), data=prefilter) 
mod9b <- lmer(Copies ~ Filter  + Tank + (1|Bottle), data=prefilter) 
mod10b <- lmer(Copies ~ Filter  + Volume + (1|Bottle), data=prefilter) 
mod11b <- lmer(Copies ~ Volume  + Tank + (1|Bottle), data=prefilter) 

mod12b <- lmer(Copies ~ Prefilter + (1|Bottle), data=prefilter) 
mod13b <- lmer(Copies ~ Tank + (1|Bottle), data=prefilter) 
mod14b <- lmer(Copies ~ Filter + (1|Bottle), data=prefilter) 
mod15b <- lmer(Copies ~ Volume + (1|Bottle), data=prefilter) 

mod16b <- lmer(Copies ~ Prefilter + Filter + Prefilter*Filter + (1|Bottle), data=prefilter)

mod17b <- lmer(Copies ~ Prefilter + Filter + Prefilter*Filter + Volume + Tank + (1|Bottle), data=prefilter)
mod18b <- lmer(Copies ~ Prefilter + Filter + Prefilter*Filter + Tank + (1|Bottle), data=prefilter)
mod19b <- lmer(Copies ~ Prefilter + Filter + Prefilter*Filter + Volume + (1|Bottle), data=prefilter)

AICctab(mod0b, mod1, mod2, mod3b, mod4b,mod5b,mod6b,mod7b,mod8b,mod9b,mod10b,mod11b,mod12b,mod13b,mod14b,mod15b,mod16b, mod17b, mod18b, mod19b, weights=TRUE, base=TRUE)



detection<- read.csv("Prefiltration_detection.csv")


mod0c <- glm(cbind(Detections, Total-Detections) ~ 1, data=detection, family = "binomial")   
mod1c <- glm(cbind(Detections, Total-Detections) ~ Turbidity + Filter + Tank + Volume, data=detection, family = "binomial") 
mod2c <- glm(cbind(Detections, Total-Detections) ~ Turbidity + Filter + Tank, data=detection, family = "binomial") 
mod3c <- glm(cbind(Detections, Total-Detections) ~ Turbidity + Filter + Volume, data=detection, family = "binomial") 
mod4c <- glm(cbind(Detections, Total-Detections) ~ Turbidity + Tank + Volume, data=detection, family = "binomial") 
mod5c <- glm(cbind(Detections, Total-Detections) ~ Filter + Tank + Volume, data=detection, family = "binomial") 


mod6c <- glm(cbind(Detections, Total-Detections) ~ Turbidity  + Filter, data=detection, family = "binomial") 
mod7c <- glm(cbind(Detections, Total-Detections) ~ Turbidity  + Volume, data=detection, family = "binomial") 
mod8c <- glm(cbind(Detections, Total-Detections) ~ Turbidity  + Tank, data=detection, family = "binomial") 
mod9c <- glm(cbind(Detections, Total-Detections) ~ Filter  + Tank, data=detection, family = "binomial") 
mod10c <- glm(cbind(Detections, Total-Detections) ~ Filter  + Volume, data=detection, family = "binomial") 
mod11c <- glm(cbind(Detections, Total-Detections) ~ Volume  + Tank, data=detection, family = "binomial") 

mod12c <- glm(cbind(Detections, Total-Detections) ~ Turbidity, data=detection, family = "binomial") 
mod13c <- glm(cbind(Detections, Total-Detections) ~ Tank, data=detection, family = "binomial") 
mod14c <- glm(cbind(Detections, Total-Detections) ~ Filter, data=detection, family = "binomial") 
mod15c <- glm(cbind(Detections, Total-Detections) ~ Volume, data=detection, family = "binomial") 

mod16c <- glm(cbind(Detections, Total-Detections) ~ Volume  + Filter + Volume*Filter, data=detection, family = "binomial")

AICctab(mod0c, mod1c, mod2c, mod3c, mod4c,mod5c,mod6c,mod7c,mod8c,mod9c,mod10c,mod11c,mod12c,mod13c,mod14c,mod15c,mod16c, weights=TRUE, base=TRUE)



#detection for turbidity

detection<- read.csv("Prefiltration_detection.csv")


mod0c <- glm(cbind(Detections, Total-Detections) ~ 1, data=detection, family = "binomial")   
mod1c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Filter + Tank + Volume, data=detection, family = "binomial") 
mod2c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Filter + Tank, data=detection, family = "binomial") 
mod3c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Filter + Volume, data=detection, family = "binomial") 
mod4c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Tank + Volume, data=detection, family = "binomial") 
mod5c <- glm(cbind(Detections, Total-Detections) ~ Filter + Tank + Volume, data=detection, family = "binomial") 


mod6c <- glm(cbind(Detections, Total-Detections) ~ Prefilter  + Filter, data=detection, family = "binomial") 
mod7c <- glm(cbind(Detections, Total-Detections) ~ Prefilter  + Volume, data=detection, family = "binomial") 
mod8c <- glm(cbind(Detections, Total-Detections) ~ Prefilter  + Tank, data=detection, family = "binomial") 
mod9c <- glm(cbind(Detections, Total-Detections) ~ Filter  + Tank, data=detection, family = "binomial") 
mod10c <- glm(cbind(Detections, Total-Detections) ~ Filter  + Volume, data=detection, family = "binomial") 
mod11c <- glm(cbind(Detections, Total-Detections) ~ Volume  + Tank, data=detection, family = "binomial") 

mod12c <- glm(cbind(Detections, Total-Detections) ~ Prefilter, data=detection, family = "binomial") 
mod13c <- glm(cbind(Detections, Total-Detections) ~ Tank, data=detection, family = "binomial") 
mod14c <- glm(cbind(Detections, Total-Detections) ~ Filter, data=detection, family = "binomial") 
mod15c <- glm(cbind(Detections, Total-Detections) ~ Volume, data=detection, family = "binomial") 

mod16c <- glm(cbind(Detections, Total-Detections) ~ Prefilter  + Filter + Prefilter*Filter, data=detection, family = "binomial")


AICctab(mod0c, mod1c, mod2c, mod3c, mod4c,mod5c,mod6c,mod7c,mod8c,mod9c,mod10c,mod11c,mod12c,mod13c,mod14c,mod15c,mod16c, weights=TRUE, base=TRUE)


###############

detection<- read.csv("Prefiltration_detection.csv")
str(detection)

mod1c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Tank + Filter + Volume + Prefilter*Filter, data=detection, family="binomial")  #AICc 140.5
mod2c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Tank + Filter + Volume, data=detection, family="binomial") #AICc 177.1
mod3c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Tank  + Filter + Prefilter*Filter, data=detection, family="binomial") #AICc 140.5
mod4c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Tank, data=detection, family="binomial") #AICc 187.4
mod5c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Filter, data=detection, family="binomial") #AICc 187.3
mod6c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Volume + Filter, data=detection, family="binomial") #AIC174.8
mod7c <- glm(cbind(Detections, Total-Detections) ~ Volume + Filter, data=detection, family="binomial")#AIC174.8
mod8c <- glm(cbind(Detections, Total-Detections) ~ Prefilter + Filter + Prefilter*Filter, data=detection, family="binomial") #138.1




  
  
  
    

