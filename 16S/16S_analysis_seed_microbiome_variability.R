#################################### Bean microbiomes variability_16S #####################################
# Date: September 26th 2019
# By : AF. Bintarti
# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("sjmisc") 
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages('mvtnorm', dep = TRUE)
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("userfriendlyscience")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
#install.packages("cowplot")
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("car")
install.packages("multcomp")
library(multcomp)
library(car)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
#library(cowplot)
library(ggplot2)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(userfriendlyscience)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)

# SET THE WORKING DIRECTORY
# setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/16S')
wd <- print(getwd())
otu <- read.table('OTU_table_tax_filt.txt', sep='\t', header=T, row.names = 1)
otu
tax <- otu[,'taxonomy']
str(tax)
#write.csv(tax, file = "taxonomy.csv")
dim(otu)
otu <- otu[,-51]
dim(otu) # otu= 273, otu table still has Mock, NC, and PC in the sample
sort(rowSums(otu, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
taxonomy = read.csv("taxonomy.edit.csv", header=T)
rownames(taxonomy) <- rownames(otu)
dim(taxonomy)
#read the metadata
map <- read.csv("bean.var.map.csv")
#select only biological sample from otu table
otu.bac <- otu[,1:47] #unselect Mock, NC, and PC from the otu table
dim(otu.bac) 
sort(rowSums(otu.bac, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in sample
otu.bac1=otu.bac[which(rowSums(otu.bac) > 0),]
dim(otu.bac1) # otu= 229, otu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(otu.bac1, na.rm = FALSE, dims = 1), decreasing = F)


#otu table of the negative control
otu.NC <- otu[,"NC",drop=FALSE]#only negative control
otu.NC

############## remove contaminant reads/otus from otu bac using microDecon package ##################

#install and load microDecon package
#install.packages("devtools")
#devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

#merge otu.NC to otu bac and put otu.NC as the first column
otu.NC <- rownames_to_column(otu.NC, var = "OTU.ID")
otu.bac1 <- rownames_to_column(otu.bac1, var = "OTU.ID")
otu.NC.bac <- merge(otu.NC, otu.bac1, by="OTU.ID")
dim(otu.NC.bac)

#decontamination
otu.decon <- decon(data = otu.NC.bac, numb.blanks = 1, numb.ind = 47, taxa = F)
otu.decon$OTUs.removed # remove 18 OTUs
otu.tab <- as.data.frame(otu.decon$decon.table)

#remove NC from the otu.dec.table 
otu.tab <- otu.tab[,c(1,3:49)]
dim(otu.tab) #there are 211 otus and 47 samples

#merge otu.dec.table with taxonomy to have match taxonomy table
taxonomy <- rownames_to_column(taxonomy, var = "OTU.ID")
otu.tab.tax <- merge(otu.tab, taxonomy, by = "OTU.ID")
otu.tab.tax

#separate the otu table and taxonomy table
#otu table
decon.otu <- data.frame(otu.tab.tax[,c(1:48)])
dim(decon.otu) #211 otus, 47samples
tax <- data.frame(otu.tab.tax[,c(1,49:55)], row.names = 1)
dim(tax) #211 otus, 7columns/levels
decon.otu <- column_to_rownames(decon.otu, var = "OTU.ID")
sort(colSums(decon.otu, na.rm = FALSE, dims = 1), decreasing = F)
sort(rowSums(decon.otu, na.rm = FALSE, dims = 1), decreasing = F)

#remove OTUs that do not present in sample: "rowSums(MRcounts(bac.data))=0"
#otu1=otu[which(rowSums(otu) > 0),]
#sort(rowSums(otu1, na.rm = FALSE, dims = 1), decreasing = F)
#dim(otu1) #211 otus, 47samples, clean otu table after decontamination before reads normalization

############## Reads normalization using metagnomeSeq ####################

#install and load metagenomeSeq package
#BiocManager::install("metagenomeSeq", dependencies=TRUE)
library(metagenomeSeq)

#Creating a MRexperiment object

#loading count data (otu)
decon.otu

#loading metadata
map <- read.csv("bean.var.map.csv")
map <- column_to_rownames(map, var="sample.id")

#create MRexperiment object 
phenotypeData <- AnnotatedDataFrame(map)
phenotypeData

OTUdata <- AnnotatedDataFrame(tax)
OTUdata

#create model
model <- newMRexperiment(decon.otu, phenoData = phenotypeData, featureData = OTUdata)
model

sort(rowSums(MRcounts(model), na.rm = FALSE, dims = 1), decreasing = T)

#normalising the data

#normalise the data to account for differences due to uneven sequencing depth
#metagenomeSeq uses Cumulative Sum Scaling (CSS) normalisation instead of rarefaction
#cumNormStatFast=Calculates the percentile for which to sum counts up to and scale by.

p <- cumNormStatFast(model, pFlag = T)

#cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.

bac.norm <- cumNorm(model, p = p)
bac.norm

#export count matrix
otu.norm <- MRcounts(bac.norm, norm = TRUE, log = F)

otu.norm <- as.data.frame(otu.norm)

head(sort(colSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otu.norm, na.rm = FALSE, dims = 1), decreasing = FALSE))

otu.norm <- rownames_to_column(otu.norm, var="OTU.ID")
#write.table(otu.norm, "otu_norm.txt", sep = "\t", quote = F, row.names = F)
dim(otu.norm) #there are 211 otus and 47 samples

### calculate alpha diversity ###


# calculate richness
otu.norm <- column_to_rownames(otu.norm, var="OTU.ID")
s <- specnumber(otu.norm, MARGIN = 2) # richness
rich <- as.data.frame(s)
bean.map <- map
bean.map$Richness <- s
bean.map

#add Faith's PD index
pd.index <- read.table('PD.txt', sep='\t', header=T, row.names = 1)
pd.index <- rownames_to_column(pd.index, "sample.id")
#join Faith's PD index to the map 
bean.map <- rownames_to_column(bean.map, var = "sample.id")
bean.map <- merge(bean.map, pd.index, by="sample.id", all = T)

###############################################################################################################################################
###############################################################################################################################################

### Compare Richness among plants and pods ###

# 1. Model fitting with linear mixed effect modelling with lme() function 

bean.map$plant <- as.factor(bean.map$plant)
bean.map$pod <- as.factor(bean.map$pod)

set.seed(13)
# Fit random intercept model
model = lme(Richness ~ plant, random = ~ 1|pod, data=bean.map, method="REML")
summary(model)
aov.mod = anova.lme(model, 
          type="marginal", 
          adjustSigma = FALSE)
aov.mod
anova(model, type='marginal')
# Fit random intercept and slope model
model2 <- lme(Richness ~ plant, random = ~ plant | pod, data = bean.map, method = "REML")
AIC(model, model2)
anova(model, model2)
# Random intercepts and slope model does not fit the data better Use simpler random intercepts model
model <- update(model, method = "REML")
summary(model)
aov.mod = anova.lme(model, 
          type="marginal", 
          adjustSigma = FALSE)
aov.mod

# Result: plant has significant effect to the bacterial/archaeal richness.

# Test the significance of the random effect in the mixed effects model
model.fixed = gls(Richness ~ plant, 
                  data=bean.map, 
                 method="REML")
anova(model, 
      model.fixed) # meanwhile pod has no significant effect.

# Checking assumptions of the model
hist(residuals(model), 
     col="darkgray") # the model is skewed

# Using the aov function for a nested anova
fit = aov(Richness ~ plant + Error(pod), data=bean.map)
summary(fit)

# Using Mean Sq and Df values to get p-value for H = plant and Error = pod
pf(q=656.9/95.0,
   df1=2,
   df2=9,
   lower.tail=FALSE)
# Using Mean Sq and Df values to get p-value for H = plant and Error = Within
summary(fit)
pf(q=95.0/58.38,
   df1=9,
   df2=35,
   lower.tail=F)

# model evaluation
# normality
plot(model)
qqnorm(resid(model))
qqline(resid(model)) #There is some deviation from from the expected normal line towards the tails, but overall the line looks straight and therefore pretty normal and suggests that the assumption is not violated
library(sjPlot)
plot_grid(plot_model(model, type = "diag"))
kurtosis(resid(model),method = 'sample')
# homoscedasticity
leveneTest(residuals(model) ~ bean.map$plant) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
boxplot(residuals(model) ~ bean.map$plant) 

# posthoc
set.seed(13)
posthoc = glht(model,
               linfct = mcp(plant="Tukey"))
posthoc

mcs = summary(posthoc,
              test=adjusted("BH"))
mcs



### RESULT: There are significant differences of bacterial and archaeal richness among plants but not among pods

# 2. Model fitting with linear mixed effect modelling with lmer function
data = bean.map %>% group_by(plant, pod) %>% summarize(Richness = mean(Richness))
data
boxplot(Richness ~ plant, data)

library(lme4)
library(lmerTest)

model.lmer = lmer(Richness ~ plant + (1|pod),
             data=bean.map,
             REML=TRUE)

model.lmer1 <- lmer(Richness ~ plant + (plant|pod), 
               data=bean.map, 
               REML = TRUE)

anova(model.lmer, model.lmer1)

plot(model.lmer)
plot(fitted(model.lmer), residuals(model.lmer, type = "pearson",
    scaled = TRUE))
qqnorm(resid(model.lmer))
qqline(resid(model.lmer))

summary(model.lmer)
anova(model.lmer)

mod1 = update(model.lmer, REML = FALSE)
mod2 = update(model.lmer, ~. - plant, REML = FALSE)
anova(mod1, mod2)

rand(model.lmer)
difflsmeans(model.lmer, 
            test.effs="plant")


# posthoc
posthoc.lmer = glht(model.lmer,
               linfct = mcp(plant="Tukey"))

mcs.lmer = summary(posthoc.lmer,
              test=adjusted("BH"))

mcs.lmer

# get the significant letter
label <- cld(mcs,
    level=0.05,
    decreasing=TRUE)

label
# make the dataframe of the group and the significant letter
label.df <- as.data.frame(label$mcletters$Letters)
names(label.df)[names(label.df) == "label$mcletters$Letters"] <- "Letter"
label.df <- rownames_to_column(label.df, "plant")
label.df
# calculate the max value of the richness and put them in the same dataframe of the group and the significant letter
sum_rich_plant <- bean.map %>%
  group_by(plant) %>% 
  summarize(max.rich=max(Richness))
sum_rich_plant_new=left_join(label.df,sum_rich_plant, by='plant')
sum_rich_plant_new

# 1. Power analysis for richness using lmer

#install.packages("simr")
library(simr)
set.seed(13)
fixef(model.lmer)["plantC"]
power <- powerSim(model.lmer, nsim=1000, test=fixed("plant"), seed=13) #92.20% (90.36, 93.79)
print(power, alpha = power$alpha, level = 0.95)
#The 3 level categorical predictor (plant) is split into 2 dummy variables when entering the model. 
#If you are interested in the effect of one specific dummy variable, you can run a z-test on it
powerB <- powerSim(model.lmer,fixed("plantB",'z'), nsim=200, seed=13)
print(powerB, alpha = power$alpha, level = 0.95)
powerC <- powerSim(model.lmer,fixed("plantC",'z'), nsim=200, seed=13) 
print(powerC, alpha = power$alpha, level = 0.95)
#The power to reject the null hypothesis of zero trend in Richness is about > 80%, 
#traditionally 80% power is considered adequate (although this arbitrary threshold is not always appropriate

# do test for random variable
doTest(model.lmer, random())
powerSim(model.lmer,random(), nsim=100, seed=13)
# To find out about the dummy variables, you can take a look at the model summary:
summary(model.lmer)$coef

# power curve
pc1 <- powerCurve(model.lmer, nsim=1000, seed=13)
summary(pc1)
print(pc1)
plot(pc1)

bean.map$obs <- 1:47
pc <- powerCurve(model.lmer, along="obs")
plot(pc)
summary(pc)

######################################################################################################################################################
######################################################################################################################################################

### Analysis of bacterial/archaeal richness using GLMM ###
### why using it? GLMM is designed for analyzing data that not normally distributed and have unbalanced design
### GLMM is also used if you have mixed variables (fixed and random variables)

# install the package and load it
# install.packages("glmmTMB")
library("glmmTMB")

# since we deal with discrete count data, then try to analyze the data using poisson and negative binomial distribution
poismod <- glmmTMB(Richness ~ plant + (1|pod), data = bean.map, family="poisson") #try poisson
nbmod <- glmmTMB(Richness ~ plant  + (1|pod), data = bean.map, family="nbinom2") #try negative binomial 
nbmod1 <- glmmTMB(Richness ~ plant + (1|pod), data = bean.map, family="nbinom1") #nbinom1 is similar to quasi-poisson

# compare the model to see which one is the most fit
library(bbmle)
AICtab(poismod, nbmod, nbmod1) # nbmod1 is the most fit model (it has the lowest AIC value)
anova(poismod, nbmod, nbmod1) 

# try to add another fixed factor " dryseed_weight_g" into the model
nbmod2 <- glmmTMB(Richness ~ plant + dryseed_weight_g + (1|pod), data = bean.map, family="nbinom1")
nbmod3 <- glmmTMB(Richness ~ plant * dryseed_weight_g + (1|pod), data = bean.map, family="nbinom1")

# compare the model to see which one is the most fit
AICtab(nbmod1, nbmod2, nbmod3) # nbmod1 is the most fit model (it has the lowest AIC value)
anova(nbmod1, nbmod2, nbmod3)

# model inference
summary(nbmod1) #there are differences of bacterial/archaeal richness among plants, plant B and C have lower richness than plant A
car::Anova(nbmod1) #get the p-value of the fixed factor

# check the significance of the random factor
m2 <- glmmTMB(Richness ~ plant, data = bean.map, family="nbinom1")
anova(nbmod1, m2) # there are no differences of bacterial/archaeal richness among pods within plant

G2 = -2 * logLik(nbmod1) + 2 * logLik(m2)
pchisq(as.numeric(G2), df=1, lower.tail=F)

# check the model 
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = nbmod1, plot = T)
#plot(simulationOutput, quantreg = T, asFactor = F)
residuals(simulationOutput)
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)
hist(simulationOutput)
testResiduals(simulationOutput)
testDispersion(nbmod1)

# posthoc
library(emmeans)
tukey.rich <- emmeans(nbmod1, "plant")
tukey.rich
pairs(tukey.rich)
pwpm(tukey.rich)
eff_size(tukey.rich, sigma = sigma(nbmod1), edf = Inf)
plot(tukey.rich, comparisons = TRUE)

# get the signif letters
rich.label <- multcomp::cld(tukey.rich,alpha = 0.05, Letters=letters)
rich.label
# calculate the max value of the richness and put them in the same dataframe of the group and the significant letter
sum_rich_plant <- bean.map %>%
  group_by(plant) %>% 
  summarize(max.rich=max(Richness))
sum_rich_plant_new=left_join(rich.label,sum_rich_plant, by='plant')
sum_rich_plant_new

################################################################################################################################
#plot richness among plants
library(viridis)
rich.plant <- ggplot(bean.map, aes(x=plant, y=Richness, fill=plant))+
                    geom_boxplot() +
                    #scale_fill_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FF","#3CBC75F","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
                    scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    geom_text(data=sum_rich_plant_new, aes(x=plant,y=2+max.rich,label=Letter), vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
rich.plant
# save the plot
ggsave("rich.plant.tiff",
       rich.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

#plot richness among pods
#install.packages("viridis")  # Install
library("viridis")           # Load
rich.p <- ggplot(bean.map, aes(x=pod, y=Richness, fill = pod))+
                    geom_boxplot() +
                    #scale_fill_manual(labels = c("A1","A2", "A3", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2", "C3"),values=c("#ffc6c4","#cc607d","#672044","#fef6b5", "#ffd795", "#ffb77f", "#fd9576", "#f37378", "#e15383", "#A9DC67", "#6EC574", "#39AB7E"))+
                    scale_fill_viridis(discrete = T)+
                    geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #labs(title = "A. Bacteria/archaea")+
                    #geom_text(data=new.rich_pod.summarized,aes(x=pod,y=0.5+max.rich,label=new.rich_pod.summarized$groups),vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
rich.p
rich.pod <- rich.p +
  facet_grid(. ~ plant, scales="free_x")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
rich.pod

ggsave("rich.pod.tiff",
       rich.pod, device = "tiff",
       width = 6, height =4.5, 
       units= "in", dpi = 600)


######################################################################################################################
######################################################################################################################

### Compare PD Whole Tree among plants ###

# 1. Model fitting with linear mixed effect modelling with lme() function 
pd.model = lme(PD_whole_tree ~ plant, random = ~ 1|pod, data=bean.map, method="REML")
summary(pd.model)
aov.pd.mod = anova.lme(pd.model, 
          type="marginal", 
          adjustSigma = FALSE)
aov.pd.mod
anova(pd.model, type='marginal')
# Fit random intercept and slope model
pd.model2 <- lme(PD_whole_tree ~ plant, random = ~ plant | pod, data = bean.map, method = "REML")
AIC(pd.model, pd.model2)
anova(pd.model, pd.model2)
# check the normality
shapiro.test(resid(pd.model)) # normal
# Result: plant has significant effect to the bacterial/archaeal richness.

# Test the significance of the random effect in the mixed effects model
pd.model.fixed = gls(PD_whole_tree ~ plant, 
                  data=bean.map, 
                 method="REML")
anova(pd.model, 
      pd.model.fixed) # meanwhile pod has no significant effect.

# Checking assumptions of the model
hist(residuals(pd.model), 
     col="darkgray") # the model is skewed

# model evaluation
# normality
plot(pd.model)
qqnorm(resid(pd.model))
qqline(resid(pd.model)) #There is some deviation from from the expected normal line towards the tails, but overall the line looks straight and therefore pretty normal and suggests that the assumption is not violated
library(sjPlot)
plot_grid(plot_model(pd.model, type = "diag"))

# homoscedasticity
leveneTest(residuals(pd.model) ~ bean.map$plant) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
boxplot(residuals(pd.model) ~ bean.map$plant) 

# posthoc
set.seed(13)
bean.map$plant <- as.factor(bean.map$plant)

pd.posthoc = glht(pd.model,
               linfct = mcp(plant="Tukey"))
pd.posthoc
pd.mcs = summary(pd.posthoc,
              test=adjusted("BH"))
pd.mcs


# 2. Model fitting with linear mixed effect modelling with lmer() function 

bean.map$plant <- as.factor(bean.map$plant)
bean.map$pod <- as.factor(bean.map$pod)

set.seed(13)
library(lme4)
library(lmerTest)
pd.model.lmer = lmer(PD_whole_tree ~ plant + (1|pod),
             data=bean.map,
             REML=TRUE)

pd.model.lmer1 <- lmer(PD_whole_tree ~ plant + (plant|pod), 
               data=bean.map, 
               REML = TRUE)

anova(pd.model.lmer, pd.model.lmer1)


# Checking assumptions of the model
plot(pd.model.lmer)
plot(fitted(pd.model.lmer), residuals(pd.model.lmer, type = "pearson",
    scaled = TRUE))
qqnorm(resid(pd.model.lmer))
qqline(resid(pd.model.lmer))
qqnorm(pd.model.lmer, ~ranef(., level=2))

# model inference
summary(pd.model.lmer)
anova(pd.model.lmer)

rand(pd.model.lmer)
difflsmeans(pd.model.lmer, 
            test.effs="plant")


# 2. Power analysis for PD using lmer

#install.packages("simr")
library(simr)
set.seed(13)

power.pd <- powerSim(pd.model.lmer, nsim=1000, test=fixed("plant"), seed=13) 
print(power.pd, alpha = power$alpha, level = 0.95) # 89.00% (86.89, 90.87)
#The 3 level categorical predictor (plant) is split into 2 dummy variables when entering the model. 
#If you are interested in the effect of one specific dummy variable, you can run a z-test on it
powerB.pd <- powerSim(pd.model.lmer,fixed("plantB",'z'), nsim=200, seed=13)
print(powerB.pd, alpha = power$alpha, level = 0.95)
powerC.pd <- powerSim(pd.model.lmer,fixed("plantC",'z'), nsim=200, seed=13) 
print(powerC.pd, alpha = power$alpha, level = 0.95)
#The power to reject the null hypothesis of zero trend in Richness is about > 80%, 
#traditionally 80% power is considered adequate (although this arbitrary threshold is not always appropriate

# do test for random variable
doTest(model.lmer, random())
powerSim(model.lmer,random(), nsim=100, seed=13)
# To find out about the dummy variables, you can take a look at the model summary:
summary(model.lmer)$coef



### RESULT: There are significant differences of bacterial and archaeal phylogenetic diversity among plants but not among pods

# posthoc

set.seed(13)
bean.map$plant <- as.factor(bean.map$plant)
pd.posthoc = glht(pd.model,
               linfct = mcp(plant="Tukey"))
pd.posthoc
pd.mcs = summary(pd.posthoc,
              test=adjusted("BH"))
pd.mcs
# get the significant letter
pd.label <- cld(pd.mcs,
    level=0.05,
    decreasing=TRUE)
pd.label
# make the dataframe of the group and the significant letter
pd.label.df <- as.data.frame(pd.label$mcletters$Letters)
names(pd.label.df)[names(pd.label.df) == "pd.label$mcletters$Letters"] <- "Letter"
pd.label.df <- rownames_to_column(pd.label.df, "plant")
pd.label.df
# calculate the max value of the phylogenetic diversity and put them in the same dataframe of the group and the significant letter
sum_pd_plant <- bean.map %>%
  group_by(plant) %>% 
  summarize(max.pd=max(PD_whole_tree))
sum_pd_plant_new=left_join(pd.label.df,sum_pd_plant, by='plant')
sum_pd_plant_new
#plot phylogenetic diversity among plant
pd.plant <- ggplot(bean.map, aes(x=plant, y=PD_whole_tree, fill=plant))+
                    geom_boxplot() +
                    #scale_fill_manual(labels = c("A", "B", "C"),values=c("#CC6677", "#DDCC77","#117733"))+
                    scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(y="Faith's Phylogenetic Diversity")+
                    geom_text(data=sum_pd_plant_new, aes(x=plant,y=0.4+max.pd,label=Letter), vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
pd.plant
# save the plot
ggsave("pd.plant.tiff",
       pd.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

#plot PD among pods

#install.packages("viridis")  # Install
library("viridis")           # Load
pd.p <- ggplot(bean.map, aes(x=pod, y=PD_whole_tree, fill = pod))+
                    geom_boxplot() +
                    #scale_fill_manual(labels = c("A1","A2", "A3", "B1", "B2", "B3", "B4", "B5", "B6", "C5", "C6", "C7"),values=c("#ffc6c4","#cc607d","#672044","#fef6b5", "#ffd795", "#ffb77f", "#fd9576", "#f37378", "#e15383", "#A9DC67", "#6EC574", "#39AB7E"))+
                    scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(y="Faith's Phylogenetic Diversity")+
                    #geom_text(data=new.rich_pod.summarized,aes(x=pod,y=0.5+max.rich,label=new.rich_pod.summarized$groups),vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
pd.pod <- pd.p +
  facet_grid(. ~ plant, scales="free_x")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
pd.pod

ggsave("pd.pod.tiff",
       pd.pod, device = "tiff",
       width = 6, height =4, 
       units= "in", dpi = 600)



### RESULT: There are no significant differences of bacterial and archaeal richness among pods after post hoc Dunn test

###########################################################################################################################
### Compare bacterial and archaeal richness among different pods of the same plants using one-way ANOVA######

plantA <- bean.map[c(1:12),]
plantA$plant <- factor(plantA$plant)
plantA$pod <- factor(plantA$pod)
plantA$sample.id <- factor(plantA$sample.id)
plantB <- bean.map[c(13:36),]
plantB$plant <- factor(plantB$plant)
plantB$pod <- factor(plantB$pod)
plantB$sample.id <- factor(plantB$sample.id)
plantC <- bean.map[c(37:47),]
plantC$plant <- factor(plantC$plant)
plantC$pod <- factor(plantC$pod)
plantC$sample.id <- factor(plantC$sample.id)

#. 1. Plant A
rich.podA <- aov(plantA$Richness ~ pod, data = plantA)
summary(rich.podA) #3.156 0.0915
# testing assumptions
# Generate residual and predicted values
rich.podA.resids <- residuals(rich.podA)
rich.podA.preds <- predict(rich.podA)
# Look at a plot of residual vs. predicted values
plot(rich.podA.resids ~ rich.podA.preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(rich.podA.resids) # p-value = 0.33, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ pod, data=plantA, na.action=na.exclude) # p-val=0.26 variances among group are homogenous
plot(density(rich.podA.resids)) 
qqnorm(rich.podA.resids)
qqline(rich.podA.resids)
hist(rich.podA.resids)
skew_xts <- skewness(rich.podA.resids)
kurtosis(rich.podA.resids,method = 'sample')
boxplot(Richness ~ pod,
        ylab="Richness", xlab="pod", data = plantA) 
#### RESULT: There are no differences of bacterial and archaeal richness among pods
(rich.podA.plot <- ggplot(plantA, aes(x=pod, y=Richness))+
                    geom_boxplot() +
                    geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "A. Plant A")+
                    # geom_text(data=new.richness.summarized,aes(x=Site,y=32+max.Richness,label=new.richness.summarized$groups),vjust=0)+
                    theme(axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()))
ggsave("rich.podA.norm.tiff",
       rich.podA.plot, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

#. 2. Plant B
rich.podB <- aov(plantB$Richness ~ pod, data = plantB)
summary(rich.podB) #3.202 0.0306 *
# testing assumptions
# Generate residual and predicted values
rich.podB.resids <- residuals(rich.podB)
rich.podB.preds <- predict(rich.podB)
# Look at a plot of residual vs. predicted values
plot(rich.podB.resids ~ rich.podB.preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(rich.podB.resids) # GOOD!! p-value = 0.72, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ pod, data=plantB, na.action=na.exclude) # p-val=0.96 variances among group are homogenous
plot(density(rich.podB.resids)) 
qqnorm(rich.podB.resids)
qqline(rich.podB.resids)
hist(rich.podB.resids)
skew_xts <- skewness(rich.podB.resids)
kurtosis(rich.podB.resids,method = 'sample')
boxplot(Richness ~ pod,
        ylab="Richness", xlab="pod", data = plantB) 
### RESULT: There are differences of bacterial and archaeal richness among pods
# Do Tukey's HSD Post Hoc Test
rich_podB.hsd <- HSD.test(rich.podB, "pod", alpha = 0.05,group = T ,main = NULL,console=TRUE)
rich_podB.hsd <- HSD.test(rich.podB, "pod",alpha = 0.05, group = FALSE, main = NULL,console=TRUE)
# Do Plot
# add significance letters from HSD.test into box plot   
rich_podB.hsd.summarized <- plantB %>% group_by(pod) %>% summarize(max.rich=max(Richness))
rich_podB.letter <- rich_podB.hsd$groups
rich_podB.letter$pod <- rownames(rich_podB.letter)
rich_podB.hsd.summ <- left_join(rich_podB.letter,rich_podB.hsd.summarized, by='pod')  
(rich.podB.plot <- ggplot(plantB, aes(x=pod, y=Richness))+
                    geom_boxplot() +
                    geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "B. Plant B")+
                    geom_text(data=rich_podB.hsd.summ,aes(x=pod,y=2+max.rich,label=rich_podB.hsd.summ$groups),vjust=0)+
                    theme(axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()))

ggsave("rich.podB.norm.tiff",
       rich.podB.plot, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

#. 3. Plant C
rich.podC <- aov(plantC$Richness ~ pod, data = plantC)
summary(rich.podC) #0.327   0.73
# testing assumptions
# Generate residual and predicted values
rich.podC.resids <- residuals(rich.podC)
rich.podC.preds <- predict(rich.podC)
# Look at a plot of residual vs. predicted values
plot(rich.podC.resids ~ rich.podC.preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(rich.podC.resids) # GOOD!! p-value = 0.13, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ pod, data=plantC) # p-val=0.66 variances among group are homogenous
plot(density(rich.podC.resids)) 
plot(rich.podC,1)
qqnorm(rich.podC.resids)
qqline(rich.podC.resids)
hist(rich.podC.resids)
skew_xts <- skewness(rich.podC.resids)
kurtosis(rich.podC.resids,method = 'sample')
boxplot(Richness ~ pod,
        ylab="Richness", xlab="pod", data = plantC) 
### RESULT: There are no differences of bacterial and archaeal richness among pods
(rich.podC.plot <- ggplot(plantC, aes(x=pod, y=Richness))+
                    geom_boxplot() +
                    geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
                    theme_bw()+
                     expand_limits(x = 0, y = 0)+
                    labs(title = "C. Plant C")+
                    # geom_text(data=new.richness.summarized,aes(x=Site,y=32+max.Richness,label=new.richness.summarized$groups),vjust=0)+
                    theme(axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()))
ggsave("rich.podC.norm.tiff",
       rich.podC.plot, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

#######################################################################################################################################################3

### Rarefaction curves ######
# using GlobalPatterns

library(phyloseq)

# 1. rarefaction curve for decontaminated non-normalized OTU table

decon.otu # decontaminated non-normalized OTU table
# make phyloseq otu table and taxonomy
otu1.phyl = otu_table(decon.otu, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object

phyl.obj1 <- merge_phyloseq(otu1.phyl,tax.phyl,map.phyl)
phyl.obj1

#setEPS()
#postscript("bacterial rarecurve.eps", height = 5, width = 5)
#rarecurve(t(otu_table(phyl.obj1)), 
                           #step=1000, cex=0.5,
                           #xlab = "Reads", 
                           #ylab = "Bacterial/archaeal OTUs", col = curve_colors)
#dev.off()

map_16S <- bean.map[bean.map$sample.id%in%colnames(decon.otu),]
curve_colors <- rep("darkgreen", ncol(decon.otu))
curve_colors[map_16S$plant=="A"] <- "#440154FF"
curve_colors[map_16S$plant=="B"] <- "#1F968BFF"
curve_colors[map_16S$plant=="C"] <- "#FDE725FF"

lty <- c("solid", "solid", "solid")

setEPS()
postscript("bacterial rarecurve decontaminated non-normalized.eps", height = 5, width = 5)
rarecurve(t(decon.otu), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", 
                           ylab = "Bacterial/archaeal OTUs")
# Add a legend
legend(2500, 20, legend=c("A", "B", "C"),
       col=c("#440154FF", "#1F968BFF","#FDE725FF"), lty = lty, lwd=2, title = "Plant")
dev.off()
graphics.off()

# 2. rarefaction curve for undecontaminated non-normalized OTU table

# read the otu table
otu.bac1 # undecontaminated non-normalized OTU table

setEPS()
postscript("bacterial rarecurve undecontaminated non-normalized.eps", height = 5, width = 5)
rarecurve(t(otu.bac1), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", 
                           ylab = "Bacterial/archaeal OTUs")
# Add a legend
legend(2500, 20, legend=c("A", "B", "C"),
       col=c("#440154FF", "#1F968BFF","#FDE725FF"), lty = lty, lwd=2, title = "Plant")
dev.off()

#######################################################################################################################################################

# BACTERIA COMPOSITION

#BiocManager::install("phyloseq")
library(phyloseq)

# load normalized otu table
otu.norm
head(otu.norm)
#otu.norm <- column_to_rownames(otu.norm, var = "OTU.ID")

# bacterial taxonomy
head(tax)
#tax <- column_to_rownames(tax, var = "OTU.ID")
rownames(tax) <- rownames(otu.norm)
# make phyloseq otu table and taxonomy
otu.phyl = otu_table(otu.norm, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj

# merge taxa by phylum

# 1. phylum - Bacteria
bac.phylum <- tax_glom(phyl.obj, taxrank = "Phylum", NArm = F)
bac.phylum.ra <- transform_sample_counts(bac.phylum, function(x) x/sum(x))
bac.phylum.ra

df.phylum <- psmelt(bac.phylum.ra) %>%
  group_by(plant, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

# OR


###create dataframe from phyloseq object
library(data.table)
df.bac.phyl <- data.table(psmelt(bac.phylum.ra))
dim(df.bac.phyl)
###convert phylum to character vector from a factor
df.bac.phyl$Phylum <- as.character(df.bac.phyl$Phylum)

# calculating mean relative abundance per plant
mean.phylrelabund.perplant= df.bac.phyl %>%
  group_by(plant, Phylum) %>%
  summarise(ra=mean(Abundance))



# barplot of bacterial/archaeal composition across pods at Phylum level
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
my_colors = carto_pal(12, "Safe")
my_colors

# New facet label names for plant variable
plant.labs <- c("Plant A", "Plant B", "Plant C")
names(plant.labs) <- c("A", "B", "C")

# Create the plot

pod.phylum <- ggplot(data=df.phylum, aes(x=pod, y=Mean, fill=Phylum))
plot.pod.phylum <- pod.phylum + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     scale_fill_manual(values= my_colors)+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Mean Relative Abundance", x="Pod")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           
plot.pod.phylum

plot.pod.phylum1 <- plot.pod.phylum +
  facet_grid(. ~ plant, labeller = labeller(plant=plant.labs), scales="free_x")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
plot.pod.phylum1

#ggsave("plot.pod.phylum.eps",
      #plot.pod.phylum1, device = "eps",
      #width = 9.5, height =5, 
      #units= "in", dpi = 600)

# merge taxa by class

# 1. class - Bacteria
bac.cl <- tax_glom(phyl.obj, taxrank = "Class", NArm = F)
bac.cl.ra <- transform_sample_counts(bac.cl, function(x) x/sum(x))
bac.cl.ra

df.cl <- psmelt(bac.cl.ra) %>%
  group_by(sample.id, plant, pod, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.cl$Class <- as.character(df.cl$Class)
df.cl$Class[df.cl$Mean < 0.1] <- "Other"

# barplot of bacterial/archaeal composition across pods at Phylum level
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
my_colors = carto_pal(12, "Safe")
my_colors

# New facet label names for plant variable
plant.labs <- c("Plant A", "Plant B", "Plant C")
names(plant.labs) <- c("A", "B", "C")

# Create the plot

pod.cl <- ggplot(data=df.cl, aes(x=pod, y=Mean, fill=Class))
plot.pod.cl <- pod.cl + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Mean Relative Abundance", x="Pod")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           #axis.ticks.x = element_blank(),
                           #axis.title.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           
plot.pod.cl

plot.pod.cl1<- plot.pod.cl +
  facet_grid(. ~ plant, labeller = labeller(plant=plant.labs), scales="free_x")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
plot.pod.cl1

ggsave("plot.pod.class.eps",
      plot.pod.cl1, device = "eps",
       width = 8.6, height =5, 
       units= "in", dpi = 600)
# subset genus Rhizobia from the phyloseq object
#subset_taxa(phyl.obj, Phylum=="Bacteroidetes")
bac.genus <- tax_glom(phyl.obj, taxrank = "Genus", NArm = F)
bac.genus.ra <- transform_sample_counts(bac.genus, function(x) x/sum(x))
bac.genus.ra
###create dataframe from phyloseq object
library(data.table)
df.bac.genus <- data.table(psmelt(bac.genus.ra))
dim(df.bac.genus)
# calculating mean relative abundance per genus
mean.relabund.genus= df.bac.genus %>%
  group_by(Genus) %>%
  summarise(ra=mean(Abundance)) 



##############################################################################################################################

# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA

# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
bacnorm_PA <- 1*(otu.norm>0)
bacnorm_PA
otu_dist <- vegdist(t(bacnorm_PA), binary = TRUE, method = "jaccard") #Sorensen

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
#env <- map[,c(11:22, 24:36)]

# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)

# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
map <- cbind(bean.map,ax1.scores,ax2.scores)

# simple plot
pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))

# PCoA Plot 
require("ggrepel")
library(ggrepel)
library(viridis)
set.seed(13)

pod.pcoa <- ggplot(data = map, aes(x=ax1.scores, y=ax2.scores))+
            theme_bw()+
            geom_point(data = map, aes(x = ax1.scores, y = ax2.scores, col=factor(plant)),size=5, alpha =0.7)+
            #scale_color_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FFF","#3CBB75FF","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
            coord_fixed() + 
            labs(colour = "Plant")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = rel(1), face="bold"),
            axis.text=element_text(size=10), 
            axis.title=element_text(size=12,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))

set.seed(13)
pod.pcoa2 <- pod.pcoa + geom_text_repel(aes(label = pod),size = 3) 
pod.pcoa2
ggsave("pcoa.jac2.tiff",
       pod.pcoa2, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

set.seed(13)
#adonis <- adonis(otu_dist ~ map$plant/map$pod, 
                 #permutation=999,
                 #method="jaccard", 
                 #strata = NULL)
#adonis
#install.packages("BiodiversityR")
library(BiodiversityR)
set.seed(13)
map$plant <- as.factor(map$plant)
map$pod <- as.factor(map$pod)
nested.npmanova(otu_dist ~ plant + pod, 
                data = map, 
                method = "jac", 
                permutations = 999)

###############################################################################################################################

## Betadisper 
groups.plant <- factor(c(rep("A",12),rep("B",24), rep("C",11)))
otu_dist <- vegdist(t(bacnorm_PA), binary = TRUE, method = "jaccard") #Sorensen
mod <- betadisper(otu_dist, groups.plant)
mod
mod$distances
dispersion <- as.data.frame(mod$distance)
names(dispersion)[names(dispersion) == "mod$distance"] <- "Dispersion"
#add dispersion index
dispersion <- rownames_to_column(dispersion, "sample.id")
#join dispersion index to the map 
bean.map <- merge(bean.map, dispersion, by="sample.id", all = T)


boxplot(mod)
# Null hypothesis of no difference in dispersion between groups
anova(mod) # there is significant differences in dispersion between groups
# the variances among groups are not homogenous,
hsd=TukeyHSD(mod) #which groups differ in relation to their variances
hsd
plot(hsd)

hsd.group=hsd$group
df.hsd.group=as.data.frame(hsd.group)
df.hsd.group=rownames_to_column(df.hsd.group, var = "Comparison")
names(df.hsd.group)[names(df.hsd.group) == "p adj"] <- "P.adj"
df.hsd.group
# get the significant letter
dis.summ.plant <- bean.map %>% group_by(plant) %>% summarize(max.dis=max(Dispersion))

hsd.letter=cldList(P.adj ~ Comparison,
        data      = df.hsd.group,
        threshold = 0.05)
names(hsd.letter)[names(hsd.letter) == "Group"] <- "plant"

new.dis.sum <- left_join(hsd.letter,dis.summ.plant,by='plant')  
new.dis.sum
#plot betadisper among plant
dis.plant <- ggplot(bean.map, aes(x=plant, y=Dispersion, fill=plant))+
                    geom_boxplot() +
                    #scale_fill_manual(labels = c("A", "B", "C"),values=c("#CC6677", "#DDCC77","#117733"))+
                    scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(y="Dispersion")+
                    geom_text(data=new.dis.sum, aes(x=plant,y=0.05+max.dis,label=Letter), vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
dis.plant
# save the plot
ggsave("dis.plant.tiff",
       dis.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

################################################################################################################################################################
## Rank-abundance Curves

#install.packages("goeveg")
library(goeveg)
racurve(t(otu.norm), main = "Rank-abundance diagram")
racurves(t(decon.otu), main = "Rank-abundance diagram", bw=F)

# OR
#install.packages("BiodiversityR")
library(BiodiversityR)
RankAbun.1 <- rankabundance(t(decon.otu))
RankAbun.1
rankabunplot(RankAbun.1,scale='abundance', addit=FALSE, specnames=c(1,2,3))
bean.map$plant <- as.factor(bean.map$plant)
bean.map$pod <- as.factor(bean.map$pod)
rankabuncomp(t(decon.otu), y=bean.map, factor='plant', type = "l",
    scale='proportion', legend=T, rainbow = T)

# OR
#install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2",force = T)
library(ampvis2)

# make ampvis data
# read the otu table and tax table
decon.otu.tax <- cbind(decon.otu, tax)
dim(decon.otu.tax)
# change "Domain" to "Kingdom"
names(decon.otu.tax)[names(decon.otu.tax) == "Domain"] <- "Kingdom"
# load the amp data
amp.data <- amp_load(decon.otu.tax,bean.map,
  fasta = NULL,
  tree = NULL,
  pruneSingletons = FALSE)
# generate rank abundance curve
amp_rankabundance(amp.data, group_by="plant", showSD = TRUE, log10_x = F)

# using phyloseq

# 1. make phyloseq object of decontaminated non-normalized otu table

library(phyloseq)
decon.otu # decontaminated non-normalized OTU table
# make phyloseq otu table and taxonomy
otu1.phyl = otu_table(decon.otu, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object

phyl.obj1 <- merge_phyloseq(otu1.phyl,tax.phyl,map.phyl)
phyl.obj1
phyl.obj1 <- merge_phyloseq(otu1.phyl,tax.phyl,map.phyl)
phyl.obj1

# this converts taxa counts in each sample to a percentage
phyloTemp <-  transform_sample_counts(phyl.obj1, function(x) x/sum(x))
clusterData <-  psmelt(phyloTemp)

# calculating mean relative abundance per plant
relabund.perplant= clusterData %>%
  group_by(plant, OTU) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(desc(ra),.by_group = TRUE) %>%
  mutate(rank = row_number(-ra))

relabund.perplant$plant <- as.factor(relabund.perplant$plant)

# plotting
rank.abund <- ggplot(relabund.perplant,aes(x=rank,y=ra, colour=plant)) +
  geom_line(aes(group = plant), size=1.2) +
  scale_colour_viridis(discrete = T)+
  labs(x="Rank", y="Relative Abundance")+
  theme_bw()+
  theme(axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title = element_text(size=15,face="bold"),
        legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# 2. make phyloseq object of decontaminated normalized otu table

# load normalized otu table
otu.norm
head(otu.norm)
#otu.norm <- column_to_rownames(otu.norm, var = "OTU.ID")

# bacterial taxonomy
head(tax)
#tax <- column_to_rownames(tax, var = "OTU.ID")
rownames(tax) <- rownames(otu.norm)
# make phyloseq otu table and taxonomy
otu.phyl = otu_table(otu.norm, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object  

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj

# this converts taxa counts in each sample to a percentage
phyloTemp2 <-  transform_sample_counts(phyl.obj, function(x) x/sum(x))
clusterData2 <-  psmelt(phyloTemp2)

# calculating mean relative abundance per plant
relabund.perplant2= clusterData2 %>%
  group_by(plant, OTU) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(desc(ra),.by_group = TRUE) %>%
  mutate(rank = row_number(-ra))

relabund.perplant2$plant <- as.factor(relabund.perplant2$plant)

# plotting
rank.abund2 <- ggplot(relabund.perplant2,aes(x=rank,y=ra, colour=plant)) +
  geom_line(aes(group = plant), size=1.2) +
  scale_colour_viridis(discrete = T)+
  labs(x="Rank", y="Relative Abundance")+
  theme_bw()+
  theme(axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title = element_text(size=15,face="bold"),
        legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
rank.abund2
ggsave("rank.abund.normalized.eps",
       rank.abund2, device = "eps",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

# 3. make phyloseq object of undecontaminated non-normalized OTU table

# read the otu table
otu.bac1 # undecontaminated non-normalized OTU table
head(otu.bac1)
otu.bac1 <- rownames_to_column(otu.bac1, var = "OTU.ID")
# bacterial taxonomy
taxonomy
head(taxonomy)
taxonomy <- rownames_to_column(taxonomy, var = "OTU.ID")
otu.bac1.tax <- left_join(otu.bac1,taxonomy, by='OTU.ID')
tax.bac1 <- data.frame(otu.bac1.tax[,c(1,49:55)], row.names = 1)
head(tax.bac1)

# make phyloseq otu table and taxonomy
otu.bac1 <- column_to_rownames(otu.bac1, var = "OTU.ID")
#tax.bac1 <- column_to_rownames(tax.bac1, var = "OTU.ID")
otu.bac1.phyl = otu_table(otu.bac1, taxa_are_rows = TRUE)
tax.bac1.phyl = tax_table(as.matrix(tax.bac1))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object  

phyl.bac1.obj <- merge_phyloseq(otu.bac1.phyl,tax.bac1.phyl,map.phyl)
phyl.bac1.obj

# this converts taxa counts in each sample to a percentage
phyloTemp3 <-  transform_sample_counts(phyl.bac1.obj, function(x) x/sum(x))
clusterData3 <-  psmelt(phyloTemp3)

# calculating mean relative abundance per plant
relabund.perplant3= clusterData3 %>%
  group_by(plant, OTU) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(desc(ra),.by_group = TRUE) %>%
  mutate(rank = row_number(-ra))
# plotting
rank.abund3 <- ggplot(relabund.perplant3,aes(x=rank,y=ra, colour=plant)) +
  geom_line(aes(group = plant), size=1.2) +
  scale_colour_viridis(discrete = T)+
  labs(x="Rank", y="Relative Abundance")+
  theme_bw()+
  theme(axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title = element_text(size=15,face="bold"),
        legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
rank.abund3
ggsave("rank.abund.nondecontaminatednonnormalized.eps",
       rank.abund3, device = "eps",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

#####################################################################################################################################
#####################################################################################################################################
# Distribution of Chloroplast and mitochondria per smaple


#read the unfiltered otu table
otu.unfil <- read.table('OTU_table_tax.txt', sep='\t', header=T, row.names = 1)
otu.unfil
tax.unfil <- otu.unfil[,'taxonomy']
tax.unfil
#write.csv(tax.unfil, file = "taxonomy.unfil.csv")
dim(otu.unfil)
otu.unfil <- otu.unfil[,-51]
dim(otu.unfil) # otu= 311, otu table still has Mock, NC, and PC in the sample
sort(rowSums(otu.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
taxonomy.unfil = read.csv("taxonomy.unfil.edit.csv", header=T)
rownames(taxonomy.unfil) <- rownames(otu.unfil)
dim(taxonomy.unfil)
#read the metadata
map <- read.csv("bean.var.map.csv")
#select only biological sample from otu table
otu.bac.unfil <- otu.unfil[,1:47] #unselect Mock, NC, and PC from the otu table
dim(otu.bac.unfil)
sort(rowSums(otu.bac.unfil, na.rm = FALSE, dims = 1), decreasing = F)
# remove OTUs that do not present in sample
otu.bac1.unfil=otu.bac.unfil[which(rowSums(otu.bac.unfil) > 0),]
dim(otu.bac1.unfil) # otu= 267, otu table before normalization using metagenomeSeq package and before decontamination
sort(rowSums(otu.bac1.unfil, na.rm = FALSE, dims = 1), decreasing = F)

library(phyloseq)

# load the otu table
head(otu.bac1.unfil)
otu.bac1.unfil <- rownames_to_column(otu.bac1.unfil, var = "OTU.ID")

# bacterial taxonomy
head(taxonomy.unfil)
taxonomy.unfil <- rownames_to_column(taxonomy.unfil, var = "OTU.ID")
otu.tax.unfil <- merge(otu.bac1.unfil, taxonomy.unfil, by="OTU.ID")
dim(otu.tax.unfil)
#separate the taxonomy table
tax.unf <- data.frame(otu.tax.unfil[,c(1,49:55)], row.names = 1)
dim(tax.unf) #267 otus, 7columns/levels
otu.bac1.unfil <- column_to_rownames(otu.bac1.unfil, var = "OTU.ID")
rownames(tax.unf) <- rownames(otu.bac1.unfil)

# make phyloseq otu table and taxonomy
otu.phyl.unf = otu_table(otu.bac1.unfil, taxa_are_rows = TRUE)
tax.phyl.unf = tax_table(as.matrix(tax.unf))

# make phyloseq map
rownames(bean.map) <- bean.map$sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object

phyl.obj.unf <- merge_phyloseq(otu.phyl.unf,tax.phyl.unf,map.phyl)
phyl.obj.unf

# merge taxa by order

# 1. order 
bac.o <- tax_glom(phyl.obj.unf, taxrank = "Order", NArm = F)
bac.o.ra <- transform_sample_counts(bac.o, function(x) x/sum(x))
bac.o.ra

df.o <- psmelt(bac.o.ra) %>%
  group_by(sample.id, plant, Order) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

#df.o$percent <- df.o$Mean*100
#df.o$Order <- as.character(df.o$Order)
#df.o$Order[df.o$Order != "Chloroplast" & df.o$Order !="NA" ] <- "Other"

#o.ra <- ggplot(data=df.o, aes(x=sample.id, y=percent, fill=Order))
#plot.o <- o.ra + 
                     geom_bar(aes(), stat="identity") + 
                     scale_fill_manual(values = my_colors3)+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Relative abundance (%)", x="Sample")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           axis.line.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
#plot.o

# OR

###create dataframe from phyloseq object
library(data.table)
df.bac.o <- data.table(psmelt(bac.o.ra))
dim(df.bac.o)
###convert phylum to character vector from a factor
df.bac.o$Order <- as.character(df.bac.o$Order)
# calculating mean relative abundance per sample
mean.bac.o.perplant= df.bac.o %>%
  group_by(Sample, Order) %>%
  summarise(ra=mean(Abundance))

# select only order Chloroplast

df.ch <- df.o %>%
  filter(Order == "Chloroplast")
df.ch$percent <- df.ch$Mean*100

# Create the plot
chlo.ra <- ggplot(data=df.ch, aes(x=sample.id, y=percent))
plot.chlo <- chlo.ra + 
                     geom_bar(aes(), stat="identity") + 
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Abundance of chloroplast (%)", x="Sample")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text.y =element_text(size=12, face = "bold"),
                           axis.text.x =element_text(size=12, face = "bold", angle = 90, hjust = 0.5, vjust=0.5),
                           axis.line.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           
plot.chlo
ggsave("plot.chlo.eps",
      plot.chlo, device = "eps",
      width = 8, height =5, 
      units= "in", dpi = 600)

# 2. family
bac.f <- tax_glom(phyl.obj.unf, taxrank = "Family", NArm = F)
bac.f.ra <- transform_sample_counts(bac.f, function(x) x/sum(x))
bac.f.ra

df.f <- psmelt(bac.f.ra) %>%
  group_by(sample.id, plant, Family) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

# select only order Mitochondria

df.mit <- df.f %>%
  filter(Family == "Mitochondria")
df.mit$percent <- df.mit$Mean*100

# Create the plot
mit.ra <- ggplot(data=df.mit, aes(x=sample.id, y=percent))
plot.mit <- mit.ra + 
                     geom_bar(aes(), stat="identity") + 
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Abundance of mitochondria (%)", x="Sample")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text.y =element_text(size=12, face = "bold"),
                           axis.text.x =element_text(size=12, face = "bold", angle = 90, hjust = 0.5, vjust=0.5),
                           axis.line.x = element_blank(),
                           axis.title =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           
plot.mit
ggsave("plot.mit.eps",
      plot.mit, device = "eps",
      width = 8, height =5, 
      units= "in", dpi = 600)












