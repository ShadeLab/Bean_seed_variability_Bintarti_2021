#################################### Bean microbiomes variability_16S #####################################
##
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
setwd('/Users/arifinabintarti/Documents/PAPER/Bean_seed_variability_Bintarti_2020/16S')
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
# merge otu.bac1 with taxonomy to have match taxonomy table
head(otu.bac1)
otu.bac1 <- rownames_to_column(otu.bac1,var = "OTU.ID")
taxonomy <- rownames_to_column(taxonomy,var = "OTU.ID")
otu.bac1.tax <- merge(otu.bac1, taxonomy, by="OTU.ID")
dim(otu.bac1.tax)
#  separate  the otu table
#otu table
otu.bac.plantfil <- data.frame(otu.bac1.tax[,c(1:48)])
head(otu.bac.plantfil)
otu.bac.plantfil <- column_to_rownames(otu.bac.plantfil,var="OTU.ID")
sum(otu.bac.plantfil)

#otu table of the negative control
otu.NC <- otu[,"NC",drop=FALSE]#only negative control
otu.NC

######################################################################################################################################
######################################################################################################################################

############## remove contaminant reads/otus from otu bac using microDecon package ##################

#install and load microDecon package
#install.packages("devtools")
#devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

#merge otu.NC to otu bac and put otu.NC as the first column
head(otu.NC)
otu.NC <- rownames_to_column(otu.NC, var = "OTU.ID")
head(otu.bac1)
#otu.bac1 <- rownames_to_column(otu.bac1, var = "OTU.ID")
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
head(taxonomy)
head(otu.tab)
#taxonomy <- rownames_to_column(taxonomy, var = "OTU.ID")
otu.tab.tax <- merge(otu.tab, taxonomy, by = "OTU.ID")
otu.tab.tax

#separate the otu table and taxonomy table
#otu table
decon.otu <- data.frame(otu.tab.tax[,c(1:48)])
dim(decon.otu) #211 otus, 47samples
tax <- data.frame(otu.tab.tax[,c(1,49:55)], row.names = 1)
dim(tax) #211 otus, 7columns/levels
head(decon.otu)
decon.otu <- column_to_rownames(decon.otu, var = "OTU.ID")
sort(colSums(decon.otu, na.rm = FALSE, dims = 1), decreasing = F)
sort(rowSums(decon.otu, na.rm = FALSE, dims = 1), decreasing = F)

#remove OTUs that do not present in sample: "rowSums(MRcounts(bac.data))=0"
#otu1=otu[which(rowSums(otu) > 0),]
#sort(rowSums(otu1, na.rm = FALSE, dims = 1), decreasing = F)
#dim(otu1) #211 otus, 47samples, clean otu table after decontamination before reads normalization

######################################################################################################################################
######################################################################################################################################

############## Reads normalization using metagnomeSeq ####################

#install and load metagenomeSeq package
#BiocManager::install("metagenomeSeq", dependencies=TRUE)
library(metagenomeSeq)

#Creating a MRexperiment object

#loading count data (otu)
decon.otu

#loading metadata
map <- read.csv("bean.var.map.csv")
head(map)
map <- column_to_rownames(map, var="Sample.id")
View(map)

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

head(otu.norm)
otu.norm <- rownames_to_column(otu.norm, var="OTU.ID")
#write.table(otu.norm, "otu_norm.txt", sep = "\t", quote = F, row.names = F)
dim(otu.norm) #there are 211 otus and 47 samples

######################################################################################################################################
######################################################################################################################################

######### FIGURE 1A & 1B ###############

### Rarefaction curves ######
# using GlobalPatterns

library(phyloseq)

# 1. rarefaction curve for decontaminated non-normalized OTU table

decon.otu # decontaminated non-normalized OTU table
# make phyloseq otu table and taxonomy
otu1.phyl = otu_table(decon.otu, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
bean.map <- map
head(bean.map)
bean.map <- rownames_to_column(bean.map, var = "Sample.id")
bean.map$Plant <- as.factor(bean.map$Plant)
bean.map$Sample.id <- as.factor(bean.map$Sample.id)
head(bean.map)
rownames(bean.map) <- bean.map$Sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object

phyl.obj1 <- merge_phyloseq(otu1.phyl,tax.phyl,map.phyl)
phyl.obj1

#bean.map <- rownames_to_column(map, var="Sample.id")
#map_16S <- bean.map[bean.map$Sample.id%in%colnames(decon.otu),]
#View(map_16S)
#curve_colors <- rep("darkgreen", ncol(decon.otu))
#curve_colors[map_16S$Plant=="A"] <- "#440154FF"
#curve_colors[map_16S$Plant=="B"] <- "#1F968BFF"
#curve_colors[map_16S$Plant=="C"] <- "#FDE725FF"

#lty <- c("solid", "solid", "solid")

#View(curve_colors)
#setEPS()
#postscript("bacterial rarecurve decontaminated non-normalized.eps", height = 5, width = 5)
#rarecurve(t(decon.otu), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", ylab = "Bacterial/archaeal OTUs")
# Add a legend
#legend(2500, 20, legend=c("A", "B", "C"),col=c("#440154FF", "#1F968BFF","#FDE725FF"), lty = lty, lwd=2, title = "Plant")
#dev.off()
#graphics.off()

#set seed
set.seed(42)

#rarefy the data
# make sure to run ggrare function in the "generating_rarecurfe.r" file
# data = phyloseq object of decontaminated non normalized otu table

# run the ggrare function attached in the file "generating_rarecurve.r"
p.rare <- ggrare(phyl.obj1, step = 1, color = "Plant", label = "Sample", se = FALSE)

#set up your own color palette
Palette <- c("#440154FF","#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(phyl.obj1)$Plant)
Palette

#plot the rarecurve
#p <- ggrare(psdata, step = 1000, color = "SampleType", label = "Sample", se = FALSE)
library(ggtext)
p.bac <- p.rare + 
 #facet_wrap(~Plant, labeller = label_both)+
 theme_bw()+
 scale_color_manual(values = Palette)+
 scale_size_manual(values = 60)+
 labs(title = "(a) Bacteria/archaea")+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.y = element_text(size=15,face="bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        #legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Number of OTUs")
 
plot(p.bac)

#setEPS()
#postscript("Fig.1ed.eps", height =4.5 , width = 8)
#par(oma = c(1.5,0,0,0))
#par(mfrow=c(1,2))
#set.seed(115)
#rarecurve(t(decon.otu), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", ylab = "Bacterial/archaeal OTUs", font.lab=2)
#title("(a)", adj=0, font.main= 2,  cex.main=2)
#rarecurve(t(otu.fil), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", ylab = "Fungal OTUs", font.lab=2)
#title("(b)", adj=0, font.main= 2,  cex.main=2)
#par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("bottom", c("A", "B", "C"), xpd = TRUE, horiz = TRUE, inset = c(0, 
    #0), bty = "n", col=c("#440154FF", "#1F968BFF","#FDE725FF"), cex = 0.8, 
    #lty = lty, lwd = 2, title = "Plant")
#dev.off()

#install.packages("iNEXT")
#library(iNEXT)

#col <- as.vector(colSums(decon.otu))

#t=decon.otu
#colSums(decon.otu)
#col <- as.vector(colSums(decon.otu))

#next.obj=iNEXT(t, q=0, datatype="abundance")
#next.obj$DataInfo
#next.obj$iNextEst

#df.next <- fortify(next.obj, type = 1)

#df.point <- df.next[which(df.next$method=="observed"),]
#df.line <- df.next[which(df.next$method!="observed"),]
#df.line$method <- factor(df.line$method, c("interpolated", "extrapolated"),c("interpolation", "extrapolation"))

#df.next.ed <- df.next %>%  mutate(Plant = case_when(startsWith(site, "A") ~ "A",startsWith(site, "B") ~ "B",startsWith(site, "C") ~ "C"))

#max(df.next.ed$x, na.rm = TRUE)
#max(df.next.ed$y, na.rm = TRUE)


#ggplot(df.next.ed, aes(x, y, group =site, colour=Plant)) + 
  #geom_line(aes(linetype=method), lwd=1.5, data=df.next.ed)+
  #scale_colour_viridis(discrete = T)+
  #labs(x="Reads", y="Bacterial/archaeal OTUs") +
  #theme(legend.position = "none") 
        #legend.title=element_blank(),
        #text=element_text(size=18),
        #legend.box = "vertical") 

  #geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  
#ggiNEXT(next.obj, type=1, facet.var = "none", color.var = "site")

#dat <- colnames(t) %>% data.frame()  
#colnames(dat)[1] <- "norep_filter_name" 
#sub_metadata <- bean.map[,1:8] %>%  
  #tibble::rownames_to_column(var = "norep_filter_name")

#dat_iNEXT <- dat %>%    
  #left_join(sub_metadata, by = "norep_filter_name") %>%
  #mutate(color = case_when(
   #startsWith(norep_filter_name, "A") ~ "#440154FF",
    #startsWith(norep_filter_name, "B") ~ "#1F968BFF",
    #startsWith(norep_filter_name, "C") ~ "#FDE725FF"))

#p_iNEXT <- ggiNEXT(next.obj, type=1, facet.var="none", se=F) + 
  #scale_color_manual(values = dat_iNEXT$color,  guide = FALSE) +
  #scale_fill_manual(values = dat_iNEXT$color, guide = FALSE)+
  #theme(legend.position = "none") +
  #theme(plot.title = element_text())

######################################################################################################################################
######################################################################################################################################

## Rank-abundance Curves

#install.packages("goeveg")
#library(goeveg)
#racurve(t(otu.norm), main = "Rank-abundance diagram")
#racurves(t(decon.otu), main = "Rank-abundance diagram", bw=F)

# OR
#install.packages("BiodiversityR")
#library(BiodiversityR)
#RankAbun.1 <- rankabundance(t(decon.otu))
#RankAbun.1
#rankabunplot(RankAbun.1,scale='abundance', addit=FALSE, specnames=c(1,2,3))
#bean.map$Plant <- as.factor(bean.map$Plant)
#bean.map$Pod <- as.factor(bean.map$Pod)
#rankabuncomp(t(decon.otu), y=bean.map, factor='Plant', type = "l",scale='proportion', legend=T, rainbow = T)

# OR
#install.packages("remotes")
#remotes::install_github("MadsAlbertsen/ampvis2",force = T)
#library(ampvis2)

# make ampvis data
# read the otu table and tax table
#decon.otu.tax <- cbind(decon.otu, tax)
#dim(decon.otu.tax)
# change "Domain" to "Kingdom"
#names(decon.otu.tax)[names(decon.otu.tax) == "Domain"] <- "Kingdom"
# load the amp data
#amp.data <- amp_load(decon.otu.tax,bean.map,fasta = NULL,tree = NULL,pruneSingletons = FALSE)
# generate rank abundance curve
#amp_rankabundance(amp.data, group_by="plant", showSD = TRUE, log10_x = F)


# OR
# using phyloseq

# make phyloseq object of decontaminated normalized otu table

# load normalized otu table
otu.norm
head(otu.norm)
otu.norm <- column_to_rownames(otu.norm, var = "OTU.ID")

# bacterial taxonomy
head(tax)
#tax <- column_to_rownames(tax, var = "OTU.ID")
rownames(tax) <- rownames(otu.norm)
# make phyloseq otu table and taxonomy
otu.phyl = otu_table(otu.norm, taxa_are_rows = TRUE)
tax.phyl = tax_table(as.matrix(tax))

# make phyloseq map
rownames(bean.map) <- bean.map$Sample.id
map.phyl <- sample_data(bean.map)

# make phyloseq object  

phyl.obj <- merge_phyloseq(otu.phyl,tax.phyl,map.phyl)
phyl.obj

# this converts taxa counts in each sample to a percentage
phyloTemp2 <-  transform_sample_counts(phyl.obj, function(x) x/sum(x))
clusterData2 <-  psmelt(phyloTemp2)

# calculating mean relative abundance per plant
relabund.perplant2= clusterData2 %>%
  group_by(Plant, Pod, OTU) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(desc(ra),.by_group = TRUE) %>%
  mutate(rank = row_number(-ra))

relabund.perplant2$Plant <- as.factor(relabund.perplant2$Plant)
relabund.perplant2$Pod <- as.factor(relabund.perplant2$Pod)
View(relabund.perplant2)

# plotting
library(viridis)

#set up your own color palette
col.pal <- c("#440154FF","#1F968BFF","#FDE725FF")
names(col.pal) <- levels(relabund.perplant2$Plant)
col.pal

rank.abund2 <- ggplot(relabund.perplant2,aes(x=rank,y=ra, colour=Plant)) +
  geom_line(aes(group = Pod), size=1) +
  scale_colour_manual(values = col.pal)+
  labs(title = "(c)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_markdown(size=15,face="bold"),
        axis.text.x=element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = 20, face='bold'),
        axis.title.x = element_blank(),
        legend.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=12, face = 'bold'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        ylab("Relative Abundance")

rank.abund2

data.rank.abund2 <- ggplot_build(rank.abund2)
gtable.rank.abund2 <- ggplot_gtable(data.rank.abund2)

#ggsave("rank.abund.normalized.eps",
      #rank.abund2, device = "eps",
      #width = 5, height =4.5, 
      #units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

##### Fig. 1 Make Rarefaction and Rank-Abundance Plots of Bacteria and Fungi in the Same Panel######

setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
library(patchwork)

#rarecurve
p.bac
p.its

Rare <- (p.bac | p.its)+xlab(label = "Sequencing Depth (Reads)")+
  theme(axis.title.x = element_text(size=15, face='bold',hjust=-2))
#Rare

#rank abundance curve
rank.abund2
rank.abund2.its

Rank <- (rank.abund2 | rank.abund2.its)+xlab(label = "Rank")+
  theme(axis.title.x = element_text(size=15, face='bold',hjust=-0.2))
#Rank

Rare.Rank <- (Rare / Rank) +
  plot_layout(guides = "collect")  & 
  theme(legend.position = 'right') &
        guides(colour=guide_legend(title.position = "top", keywidth = 2, ncol=1, bycol = T,title.hjust =0.5,
                                   title.theme = element_text(size=12, face='bold'),
                                   label.theme = element_text(size = 12, face = 'bold')))
Rare.Rank

Rare.Rank
ggsave("Fig.1.eps",
       Rare.Rank, device = "eps",
       width = 9, height = 8, 
       units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

### calculate alpha diversity ###


# calculate richness
head(otu.norm)
#otu.norm <- column_to_rownames(otu.norm, var="OTU.ID")
s <- specnumber(otu.norm, MARGIN = 2) # richness
rich <- as.data.frame(s)
bean.map <- map
dim(bean.map)
bean.map$Richness <- s
bean.map

#add Faith's PD index
pd.index <- read.table('PD.txt', sep='\t', header=T, row.names = 1)
pd.index <- rownames_to_column(pd.index, "Sample.id")
#join Faith's PD index to the map 
bean.map <- rownames_to_column(bean.map, var = "Sample.id")
bean.map <- merge(bean.map, pd.index, by="Sample.id", all = T)

###############################################################################################################################################
###############################################################################################################################################

### Compare Richness among plants and pods ###

# 1. Model fitting with linear mixed effect modelling with lme() function 

bean.map$Plant <- as.factor(bean.map$Plant)
bean.map$Pod <- as.factor(bean.map$Pod)

set.seed(13)
# Fit random intercept model
model = lme(Richness ~ Plant, random = ~ 1|Pod, data=bean.map, method="REML")
summary(model)
aov.mod = anova.lme(model, 
          type="marginal", 
          adjustSigma = FALSE)
aov.mod
anova(model, type='marginal') #  for the unbalanced data
# Fit random intercept and slope model
model2 <- lme(Richness ~ Plant, random = ~ Plant | Pod, data = bean.map, method = "REML")
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
model.fixed = gls(Richness ~ Plant, 
                  data=bean.map, 
                 method="REML")
anova(model, 
      model.fixed) # meanwhile pod has no significant effect.

# Checking assumptions of the model
hist(residuals(model), 
     col="darkgray") # the model is skewed

# Using the aov function for a nested anova
fit = aov(Richness ~ Plant + Error(Pod), data=bean.map)
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
leveneTest(residuals(model) ~ bean.map$Plant) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
boxplot(residuals(model) ~ bean.map$Plant) 

# posthoc
set.seed(13)
posthoc = glht(model,
               linfct = mcp(Plant="Tukey"))
posthoc

mcs = summary(posthoc,
              test=adjusted("BH"))
mcs

### RESULT: There are significant differences of bacterial and archaeal richness among plants but not among pods

# get the significant letter
label <- cld(mcs,
    level=0.05,
    decreasing=TRUE)

label

# make the dataframe of the group and the significant letter
label.df <- as.data.frame(label$mcletters$Letters)
names(label.df)[names(label.df) == "label$mcletters$Letters"] <- "Letter"
label.df <- rownames_to_column(label.df, "Plant")
label.df

# calculate the max value of the richness and put them in the same dataframe of the group and the significant letter
sum_rich_plant <- bean.map %>%
  group_by(Plant) %>% 
  summarize(max.rich=max(Richness))
sum_rich_plant_new=left_join(label.df,sum_rich_plant, by='Plant')
sum_rich_plant_new

##################################################################################################################################
##################################################################################################################################

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
  group_by(Plant) %>% 
  summarize(max.rich=max(Richness))
sum_rich_plant_new=left_join(rich.label,sum_rich_plant, by='Plant')
sum_rich_plant_new


##################################################################################################################################
##################################################################################################################################

#plot richness among plants
library(viridis)
#install.packages("ggtext")
library(ggtext)
rich.plant <- ggplot(bean.map, aes(x=Plant, y=Richness, fill=Plant))+
                    scale_fill_viridis(discrete = T)+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #scale_fill_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FF","#3CBC75F","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
                    #scale_colour_viridis(discrete = T)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5, aes(colour=factor(Plant)))+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=14.5+max.rich,label=Letter), size=5,vjust=0)+ 
                    geom_signif(comparisons = list(c("A", "B")), 
                                annotations = "**", textsize = 6,
                                y_position = 51, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    geom_signif(comparisons = list(c("A", "C")), 
                                annotations = "**", textsize = 6,
                                y_position =60, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    labs(title = "(a)")+
                    ylab("Bacterial/archaeal<br>Richness\n")+
                    theme(legend.position="none",
                          axis.text.x=element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20, face = 'bold'),
                          #axis.title.y=element_text(size=13,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=13, color="red", shape=95)
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          
                          
rich.plant
# save the plot
ggsave("rich.plant.tiff",
       rich.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

#plot richness among pods

#install.packages("viridis")  # Install
library("viridis")           # Load
rich.p <- ggplot(bean.map, aes(x=Pod, y=Richness, fill = Pod))+
                    #geom_boxplot() +
                    #scale_fill_manual(labels = c("A1","A2", "A3", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2", "C3"),values=c("#ffc6c4","#cc607d","#672044","#fef6b5", "#ffd795", "#ffb77f", "#fd9576", "#f37378", "#e15383", "#A9DC67", "#6EC574", "#39AB7E"))+
                    scale_fill_viridis(discrete = T)+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "(b)")+
                    #geom_text(data=new.rich_pod.summarized,aes(x=pod,y=0.5+max.rich,label=new.rich_pod.summarized$groups),vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y = element_blank(),
                          plot.title = element_text(size = 20, face = 'bold'),
                          axis.title=element_blank(),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=6, color="red", shape=95)
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          
rich.p
rich.pod <- rich.p + 
  facet_grid(. ~ Plant, scales="free_x",labeller = label_both)+
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
pd.model = lme(PD_whole_tree ~ Plant, random = ~ 1|Pod, data=bean.map, method="REML")
summary(pd.model)
aov.pd.mod = anova.lme(pd.model, 
          type="marginal", 
          adjustSigma = FALSE)
aov.pd.mod
anova(pd.model, type='marginal')

# Fit random intercept and slope model
pd.model2 <- lme(PD_whole_tree ~ Plant, random = ~ Plant | Pod, data = bean.map, method = "REML")
AIC(pd.model, pd.model2)
anova(pd.model, pd.model2)

# check the normality
shapiro.test(resid(pd.model)) # normal
# Result: plant has significant effect to the bacterial/archaeal richness.

# Test the significance of the random effect in the mixed effects model
pd.model.fixed = gls(PD_whole_tree ~ Plant, 
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
leveneTest(residuals(pd.model) ~ bean.map$Plant) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
boxplot(residuals(pd.model) ~ bean.map$Plant) 

# posthoc
set.seed(13)
bean.map$Plant <- as.factor(bean.map$Plant)

pd.posthoc = glht(pd.model,
               linfct = mcp(Plant="Tukey"))
pd.posthoc
pd.mcs = summary(pd.posthoc,
              test=adjusted("BH"))
pd.mcs


##################################################################################################################################
##################################################################################################################################

# 2. Power analysis for PD using lmer

#install.packages("simr")
#library(simr)
#set.seed(13)

#power.pd <- powerSim(pd.model.lmer, nsim=1000, test=fixed("plant"), seed=13) 
#print(power.pd, alpha = power$alpha, level = 0.95) # 89.00% (86.89, 90.87)
#The 3 level categorical predictor (plant) is split into 2 dummy variables when entering the model. 
#If you are interested in the effect of one specific dummy variable, you can run a z-test on it
#powerB.pd <- powerSim(pd.model.lmer,fixed("plantB",'z'), nsim=200, seed=13)
#print(powerB.pd, alpha = power$alpha, level = 0.95)
#powerC.pd <- powerSim(pd.model.lmer,fixed("plantC",'z'), nsim=200, seed=13) 
#print(powerC.pd, alpha = power$alpha, level = 0.95)
#The power to reject the null hypothesis of zero trend in Richness is about > 80%, 
#traditionally 80% power is considered adequate (although this arbitrary threshold is not always appropriate

# do test for random variable
#doTest(model.lmer, random())
#powerSim(model.lmer,random(), nsim=100, seed=13)
# To find out about the dummy variables, you can take a look at the model summary:
#summary(model.lmer)$coef

##################################################################################################################################
##################################################################################################################################


### RESULT: There are significant differences of bacterial and archaeal phylogenetic diversity among plants but not among pods

# posthoc

set.seed(13)
bean.map$Plant <- as.factor(bean.map$Plant)
pd.posthoc = glht(pd.model,
               linfct = mcp(Plant="Tukey"))
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
pd.label.df <- rownames_to_column(pd.label.df, "Plant")
pd.label.df

# calculate the max value of the phylogenetic diversity and put them in the same dataframe of the group and the significant letter
sum_pd_plant <- bean.map %>%
  group_by(Plant) %>% 
  summarize(max.pd=max(PD_whole_tree))
sum_pd_plant_new=left_join(pd.label.df,sum_pd_plant, by='Plant')
sum_pd_plant_new
#plot phylogenetic diversity among plant
pd.plant <- ggplot(bean.map, aes(x=Plant, y=PD_whole_tree, fill=Plant))+
                    #geom_boxplot() +
                    #scale_fill_manual(labels = c("A", "B", "C"),values=c("#CC6677", "#DDCC77","#117733"))+
                    scale_fill_viridis(discrete = T)+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "(c)")+
                    ylab("Bacterial/archaeal<br>Faith's PD\n")+
                    #geom_text(data=sum_pd_plant_new, aes(x=Plant,y=1.8+max.pd,label=Letter), vjust=0, size=5)+
                    geom_signif(comparisons = list(c("A", "B")), 
                                annotations = "**", textsize = 6,
                                y_position = 7.2, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    geom_signif(comparisons = list(c("A", "C")), 
                                annotations = "*", textsize = 6,
                                y_position =7.8, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    theme(legend.position="none",
                          axis.text.x=element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20,face = 'bold'),
                          axis.title.x =element_blank(),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=13, color="red", shape=95)
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          

pd.plant
# save the plot
ggsave("pd.plant.tiff",
       pd.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)

#plot PD among pods

#install.packages("viridis")  # Install
library("viridis")           # Load
pd.p <- ggplot(bean.map, aes(x=Pod, y=PD_whole_tree, fill = Pod))+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #geom_boxplot() +
                    #scale_fill_manual(labels = c("A1","A2", "A3", "B1", "B2", "B3", "B4", "B5", "B6", "C5", "C6", "C7"),values=c("#ffc6c4","#cc607d","#672044","#fef6b5", "#ffd795", "#ffb77f", "#fd9576", "#f37378", "#e15383", "#A9DC67", "#6EC574", "#39AB7E"))+
                    scale_fill_viridis(discrete = T)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title="(d)")+
                    #geom_text(data=new.rich_pod.summarized,aes(x=pod,y=0.5+max.rich,label=new.rich_pod.summarized$groups),vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_blank(),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20, face = 'bold'),
                          axis.title.x=element_blank(),
                          axis.title.y = element_blank(),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=6, color="red", shape=95)
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
pd.p
pd.pod <- pd.p +
  facet_grid(. ~ Plant, scales="free_x")+
  theme(strip.text = element_blank())
  #theme(strip.background =element_rect(fill="grey"))+
  #theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
pd.pod

ggsave("pd.pod.tiff",
       pd.pod, device = "tiff",
       width = 6, height =4, 
       units= "in", dpi = 600)



### RESULT: There are no significant differences of bacterial and archaeal richness among pods after post hoc Dunn test


#######################################################################################################################################################
#######################################################################################################################################################

##### Fig. 2 Make Bacterial and Fungal Richness and PD Plots Among Plants and Pods in the Same Panel######

rich.plant
rich.pod
pd.plant
pd.pod
fg.rich.plant
fg.rich.pod1
setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
#RichPD <- ggarrange(rich.plant,rich.pod,pd.plant,pd.pod, legend=NULL, nrow=2,ncol = 2, align = "hv")
#library(cowplot)
library(patchwork)
RichPD <- (rich.plant | rich.pod ) / (pd.plant | pd.pod) / (fg.rich.plant | fg.rich.pod1)
#RichPD <- plot_grid(rich.plant,rich.pod,pd.plant,pd.pod, align = "hv", axis = "btrl")
ggsave("Fig.2.tiff",
       RichPD, device = "tiff",
       width = 14, height = 12.5, 
       units= "in", dpi = 600)

#######################################################################################################################################################
#######################################################################################################################################################

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
            geom_point(data = map, aes(x = ax1.scores, y = ax2.scores, col=factor(Plant)),size=5, alpha =0.7)+
            #scale_color_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FFF","#3CBB75FF","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Plant",  title = "(a) Bacteria/archaea")+
            theme(legend.position="none",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            #plot.subtitle = element_text(size = 20, face = 'bold'),
            axis.text=element_text(size=14), 
            axis.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=12),
            legend.title = element_text(size = 12),
            legend.spacing.x = unit(0.05, 'cm'))

set.seed(13)
pod.pcoa2 <- pod.pcoa + geom_text_repel(aes(label = Pod),size = 3, max.overlaps = Inf) 
pod.pcoa2
ggsave("pcoa.jac2.tiff",
       pod.pcoa2, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

set.seed(13)

######## Calculated the statistical analysis of beta diversity using nested permanova #########
set.seed(13)
otu_dist
adonis <- adonis(otu_dist ~ map$Plant/map$Pod, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis
#install.packages("BiodiversityR")
library(BiodiversityR)
set.seed(13)
map$Plant <- as.factor(map$Plant)
map$Pod <- as.factor(map$Pod)
nested.npmanova(otu_dist ~ Plant + Pod, 
                data = map, 
                method = "jac", 
                permutations = 999)

###############################################################################################################################

## Betadisper 
set.seed(13)
groups.plant <- factor(c(rep("A",12),rep("B",24), rep("C",11)))
otu_dist <- vegdist(t(bacnorm_PA), binary = TRUE, method = "jaccard") #Sorensen
mod <- betadisper(otu_dist, groups.plant)
mod
plot(mod)
mod$distances
dispersion <- as.data.frame(mod$distance)
names(dispersion)[names(dispersion) == "mod$distance"] <- "Dispersion"
#add dispersion index
dispersion <- rownames_to_column(dispersion, "Sample.id")
#join dispersion index to the map 
bean.map <- merge(bean.map, dispersion, by="Sample.id", all = T)

set.seed(13)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod <- permutest(mod, permutations = 999, pairwise = T)
permod # there is significant differences in dispersion between groups

hsd=TukeyHSD(mod) #which groups differ in relation to their variances
hsd
plot(hsd)

hsd.group=hsd$group
df.hsd.group=as.data.frame(hsd.group)
df.hsd.group=rownames_to_column(df.hsd.group, var = "Comparison")
names(df.hsd.group)[names(df.hsd.group) == "p adj"] <- "P.adj"
df.hsd.group
# get the significant letter
detach(package:plyr)
library(dplyr)
dis.summ.plant <- bean.map %>% group_by(Plant) %>% summarize(max.dis=max(Dispersion))

hsd.letter=cldList(P.adj ~ Comparison,
        data      = df.hsd.group,
        threshold = 0.05)
names(hsd.letter)[names(hsd.letter) == "Group"] <- "Plant"

new.dis.sum <- left_join(hsd.letter,dis.summ.plant,by='Plant')  
new.dis.sum
#plot betadisper among plant
dis.plant <- ggplot(bean.map, aes(x=Plant, y=Dispersion, fill=Plant))+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #scale_fill_manual(labels = c("A", "B", "C"),values=c("#CC6677", "#DDCC77","#117733"))+
                    scale_fill_viridis(discrete = T)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(y="Dispersion", title = "(c)")+
                    #geom_text(data=new.dis.sum, aes(x=Plant,y=0.05+max.dis,label=Letter), vjust=0)+
                    geom_signif(comparisons = list(c("A", "B")), 
                                annotations = "***", textsize = 6,
                                y_position = 0.82, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    geom_signif(comparisons = list(c("A", "C")), 
                                annotations = "***", textsize = 6,
                                y_position =0.86, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20,  face = 'bold'),
                          axis.title.y =element_text(size=18,face="bold"),
                          axis.title.x = element_blank(),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=13, color="red", shape=95)
dis.plant
# save the plot
ggsave("dis.plant.tiff",
       dis.plant, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)


groups.pod <- factor(c(rep("A1",4),rep("A2",4), rep("A3",4), rep("B1",4), rep("B2",4), rep("B3",4), rep("B4",4),rep("B5",4), rep("B6",4), rep("C5",3), rep("C6",4), rep("C7",4)))
otu_dist <- vegdist(t(bacnorm_PA), binary = TRUE, method = "jaccard") #Sorensen
mod.pod <- betadisper(otu_dist, groups.pod)
mod.pod
mod.pod$distances
dispersion.pod <- as.data.frame(mod.pod$distance)
names(dispersion.pod)[names(dispersion.pod) == "mod.pod$distance"] <- "Dispersion"
#add dispersion index
dispersion.pod <- rownames_to_column(dispersion.pod, "Sample.id")
#join dispersion index to the map 
bean.map <- merge(bean.map, dispersion.pod, by="Sample.id", all = T)

boxplot(mod.pod)
# Null hypothesis of no difference in dispersion between groups
anova(mod.pod) # there is significant differences in dispersion between groups
# the variances among groups are not homogenous,
hsd.pod=TukeyHSD(mod.pod) #which groups differ in relation to their variances
hsd.pod
plot(hsd.pod)


#######################################################################################################################################################
#######################################################################################################################################################

##### Fig. 3 Make Bacterial and Fungal PCoA Plots and Dispersion Plots in the Same Panel######

pod.pcoa2
pod.pcoa.its2
dis.plant
disperplot.its

setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')

library(patchwork)


PCoA.Beta <- (pod.pcoa2 | pod.pcoa.its2 ) / (dis.plant | disperplot.its)
PCoA.Beta

PCoA.Beta2 <- patchwork::patchworkGrob(PCoA.Beta)
PCoA.Beta3 <- gridExtra::grid.arrange(PCoA.Beta2, bottom =textGrob("Plant", gp=gpar(fontsize=18,fontface='bold')))

ggsave("Fig.3.tiff",
       PCoA.Beta3, device = "tiff",
       width = 10, height = 9, 
       units= "in", dpi = 600)

#######################################################################################################################################################
#######################################################################################################################################################
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
rownames(bean.map) <- bean.map$Sample.id
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
  group_by(Plant, Pod, Phylum) %>%
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
  group_by(Plant, Pod, Phylum) %>%
  summarise(ra=mean(Abundance))%>%
  arrange(-ra)



# barplot of bacterial/archaeal composition across pods at Phylum level
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
my_colors = carto_pal(12, "Safe")
my_colors

# New facet label names for plant variable
plant.labs <- c("Plant A", "Plant B", "Plant C")
names(plant.labs) <- c("A", "B", "C")

# Create the plot

pod.phylum <- ggplot(data=df.phylum, aes(x=Pod, y=Mean, fill=Phylum))
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
  facet_grid(. ~ Plant, labeller = labeller(Plant=plant.labs), scales="free_x")+
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
  group_by(Sample.id, Plant, Pod, Class) %>%
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
plant.labs <- c("Plant: A", "Plant: B", "Plant: C")
names(plant.labs) <- c("A", "B", "C")

# Create the plot

pod.cl <- ggplot(data=df.cl, aes(x=Pod, y=Mean, fill=Class))
plot.pod.cl <- pod.cl + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888"))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     #labs(y= "Mean Relative Abundance", x="Plant")+
                     labs(y= "Mean Relative Abundance", x="Pod", title = "(a) Bacteria/archaea")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
                           axis.line.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y =  element_markdown(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           #strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           #facet_grid(~plant, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                           
plot.pod.cl

plot.pod.cl1<- plot.pod.cl +
  facet_grid(. ~ Plant, labeller = labeller(Plant=plant.labs), scales="free_x")+
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
  summarise(ra=mean(Abundance)) %>%
  arrange(-ra)

sum(mean.relabund.genus$ra)


df.bac.cl <- data.table(psmelt(bac.cl.ra))
mean.relabund.cl= df.bac.cl %>%
  group_by(Plant, Class) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(-ra)



##############################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

##### Fig. 4 Make Bacterial and Fungal Relative Abundance Plots in the Same Panel######

plot.pod.cl1
plot.fg.pod.cl1

setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')

library(patchwork)

Fig.4 <- plot.pod.cl1 / plot.fg.pod.cl1
Fig.4

ggsave("Fig.4.eps",
       Fig.4, device = "eps",
       width = 8.8, height = 7.5, 
       units= "in", dpi = 600)

#######################################################################################################################################################
#######################################################################################################################################################












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
dim(otu.bac1.unfil)
sum(otu.bac1.unfil)
rownames(tax.unf) <- rownames(otu.bac1.unfil)

# make phyloseq otu table and taxonomy
otu.phyl.unf = otu_table(otu.bac1.unfil, taxa_are_rows = TRUE)
tax.phyl.unf = tax_table(as.matrix(tax.unf))

# make phyloseq map
rownames(bean.map) <- bean.map$Sample.id
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
  group_by(Sample.id, Plant, Order) %>%
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
chlo.ra <- ggplot(data=df.ch, aes(x=Sample.id, y=percent))
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

######################################################################################################################################################
######################################################################################################################################################
# Read Proportion of Chloroplast and mitochondria 


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

# load the otu table
head(otu.bac1.unfil)
otu.bac1.unfil <- rownames_to_column(otu.bac1.unfil, var = "OTU.ID")

# merge the taxonomy with otu table
head(taxonomy.unfil)
taxonomy.unfil <- rownames_to_column(taxonomy.unfil, var = "OTU.ID")
otu.tax.unfil <- merge(otu.bac1.unfil, taxonomy.unfil, by="OTU.ID")
dim(otu.tax.unfil)

#select only the otu table and "Order"  & "Family"
otu.tax.unfil.ed <- otu.tax.unfil[,c(1:48,52,53)]
colnames(otu.tax.unfil.ed)

#edit the taxonomy
otu.tax.unfil.ed1 <- otu.tax.unfil.ed %>%
    mutate(Taxonomy = case_when(Order == "Chloroplast" ~ 'Chloroplast',
                                  Family == "Mitochondria" ~ 'Mitochondria',
                                  TRUE ~ 'Bacteria/archaea')) %>%
    mutate(Kingdom = case_when(Order == "Chloroplast" ~ 'Plant',
                                  Family == "Mitochondria" ~ 'Plant',
                                  TRUE ~ 'Bacteria/archaea'))

tail(otu.tax.unfil.ed1)
otu.tax.unfil.ed2 <- otu.tax.unfil.ed1[,c(1:48,51:52)]
colnames(otu.tax.unfil.ed2)
tail(otu.tax.unfil.ed2)

long.dat <- gather(otu.tax.unfil.ed2, Sample, Read, A11:C74, factor_key = T)
long.dat


df.unfil <- long.dat %>%
  group_by(Sample, Kingdom) %>%
  summarize(read.number = sum(Read))
df.unfil1 <- df.unfil %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

with(df.unfil1, sum(percent[Sample ==  "A11"]))

library(ggbeeswarm)
plot.unfil.king <- ggplot(df.unfil1, aes(x=Kingdom, y=percent, fill=Kingdom))+
                    geom_violin(trim = F, scale="width") +
                    #geom_beeswarm(dodge.width = 1, alpha = 0.3)+
                    scale_fill_manual(labels = c("Bacteria/archaea","Plant"),values=c("#88CCEE", "#117733"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=2+max.rich,label=Letter), vjust=0)+
                    labs(title = "A")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          #axis.text.x=element_blank(),
                          #axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          #axis.title.y=element_text(size=13,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
                          #width=1, position=position_dodge(),show.legend = FALSE)

plot.unfil.king
setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
ggsave("plant.tiff",
       plot.unfil.king, device = "tiff",
       width = 5, height =5, 
       units= "in", dpi = 600)


df.unfil.tax <- long.dat %>%
  group_by(Sample, Taxonomy) %>%
  summarize(read.number = sum(Read))
df.unfil.tax1 <- df.unfil.tax %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

install.packages("ggbeeswarm")
library(ggbeeswarm)
plot.unfil.tax <- ggplot(df.unfil.tax1, aes(x=Taxonomy, y=percent, fill=Taxonomy))+
                    geom_violin(trim = F, scale="width") +
                    #geom_beeswarm(dodge.width = 1, alpha = 0.3)+
                    #scale_fill_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FF","#3CBC75F","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
                    #scale_fill_viridis(discrete = T)+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.3)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=2+max.rich,label=Letter), vjust=0)+
                    labs(title = "A")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          #axis.text.x=element_blank(),
                          #axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          #axis.title.y=element_text(size=13,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
                          #width=1, position=position_dodge(),show.legend = FALSE)

plot.unfil.tax
setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
ggsave("chlomito.tiff",
       plot.unfil.tax, device = "tiff",
       width = 5.5, height =5, 
       units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

##### Fig.  Make plot of Plant Microbes sequence proportion in the Same Panel######

plot.unfil.king
plot.unfil.its
setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
library(patchwork)
plantmicro.plot <- plot.unfil.king | plot.unfil.its
plantmicro.plot
ggsave("Fig.tiff",
       plantmicro.plot, device = "tiff",
       width = 7, height = 4, 
       units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

#### 1. Occupancy across all seed samples #####
setwd('/Users/arifinabintarti/Documents/PAPER/Bean_seed_variability_Bintarti_2020/16S')
wd <- print(getwd())

# load normalized otu table
otu.norm
head(otu.norm)
dim(otu.norm)
otu.norm <- column_to_rownames(otu.norm, var = "OTU.ID")

# load bacterial taxonomy
head(tax)
dim(tax)
tax <- rownames_to_column(tax, var = "OTU.ID")


# load metadata
map <- read.csv("bean.var.map.csv")
head(map)

##build a long data frame joining rarefied otu table, map root, and tax root
longDF.seed <- data.frame(OTU.ID=as.factor(rownames(otu.norm)), otu.norm) %>%
  gather(Sample.id, abun, -OTU.ID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'Sample.id' column)
  left_join(tax)  %>% #adding the taxonomy info (grouped by the 'OTU.ID' column)
  group_by(OTU.ID, Sample.id) %>%
  #BioProject, plant_family, plant_genus, host, compartment
  summarise(n=sum(abun))
  #mutate(n=if_else(n>0, 1,0))

##build the new table: OTU.ID as rownames and Sample.id as colnames
wideDF.seed <- as.data.frame(spread(longDF.seed, OTU.ID, n, fill=0))
rownames(wideDF.seed) <-  wideDF.seed[,1]
wideDF.seed <- wideDF.seed[,-1]
wideDF.seed <- t(wideDF.seed)

##calculate the occupancy of each OTU.ID across Sample.id/Seed
wideDF.seed.PA <- 1*((wideDF.seed>0)==1)
Occ.seed <- rowSums(wideDF.seed.PA)/ncol(wideDF.seed.PA)
df.Occ.seed <- as.data.frame(Occ.seed)
df.Occ.seed <- rownames_to_column(df.Occ.seed, var = "OTU.ID")
df.Occ.seed.tax <- merge(df.Occ.seed, tax, by="OTU.ID")
sort.df.Occ.seed.tax <- df.Occ.seed.tax[order(df.Occ.seed.tax$Occ.seed, decreasing = TRUE),]
write.csv(sort.df.Occ.seed.tax, file = "sort.df.Occ.seed.tax.csv")

### 2. Occupancy across 3 plants

longDF.plant <- data.frame(OTU.ID=as.factor(rownames(otu.norm)), otu.norm) %>%
  gather(Sample.id, abun, -OTU.ID) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(map) %>%  #will add the info form mapping file (grouped by the 'Sample.id' column)
  left_join(tax)  %>% #adding the taxonomy info (grouped by the 'OTU.ID' column)
  group_by(OTU.ID, Plant) %>%
  summarise(n=sum(abun))
  #mutate(n=if_else(n>0, 1,0))

##build the new table: OTU.ID as rownames and Plant as colnames
wideDF.plant <- as.data.frame(spread(longDF.plant, OTU.ID, n, fill=0))
rownames(wideDF.plant) <-  wideDF.plant[,1]
wideDF.plant <- wideDF.plant[,-1]
wideDF.plant <- t(wideDF.plant)
dim(wideDF.plant) #211 taxa, 3 plant 

##calculate the occupancy of each otu across plants
wideDF.plant.PA <- 1*((wideDF.plant>0)==1)
Occ.plant <- rowSums(wideDF.plant.PA)/ncol(wideDF.plant.PA)
df.Occ.plant <- as.data.frame(Occ.plant)
df.Occ.plant <- rownames_to_column(df.Occ.plant, var = "OTU.ID")
df.Occ.plant.tax <- merge(df.Occ.plant, tax, by="OTU.ID")
sort.df.Occ.plant.tax <- df.Occ.plant.tax[order(df.Occ.plant.tax$Occ.plant, decreasing = TRUE),]
write.csv(sort.df.Occ.plant.tax, file = "sort.df.Occ.plant.tax.csv")

### 3. Occupancy across seeds within plant

#  1. Plant A
otu.normA <- data.frame(otu.norm[,c(1:12)])
##calculate the occupancy of each otu across seeds within plant
otu.normA.PA <- 1*((otu.normA>0)==1)
Occ.A <- rowSums(otu.normA.PA)/ncol(otu.normA.PA)
df.Occ.A <- as.data.frame(Occ.A)
df.Occ.A <- rownames_to_column(df.Occ.A, var = "OTU.ID")
df.Occ.A.tax <- merge(df.Occ.A, tax, by="OTU.ID")
sort.df.Occ.A.tax <- df.Occ.A.tax[order(df.Occ.A.tax$Occ.A, decreasing = TRUE),]
write.csv(sort.df.Occ.A.tax, file = "sort.df.Occ.A.tax.csv")

#  2. Plant B
otu.normB <- data.frame(otu.norm[,c(13:36)])
##calculate the occupancy of each otu across seeds within plant
otu.normB.PA <- 1*((otu.normB>0)==1)
Occ.B <- rowSums(otu.normB.PA)/ncol(otu.normB.PA)
df.Occ.B <- as.data.frame(Occ.B)
df.Occ.B <- rownames_to_column(df.Occ.B, var = "OTU.ID")
df.Occ.B.tax <- merge(df.Occ.B, tax, by="OTU.ID")
sort.df.Occ.B.tax <- df.Occ.B.tax[order(df.Occ.B.tax$Occ.B, decreasing = TRUE),]
write.csv(sort.df.Occ.B.tax, file = "sort.df.Occ.B.tax.csv")

#  3. Plant C
otu.normC <- data.frame(otu.norm[,c(37:47)])
##calculate the occupancy of each otu across seeds within plant
otu.normC.PA <- 1*((otu.normC>0)==1)
Occ.C <- rowSums(otu.normC.PA)/ncol(otu.normC.PA)
df.Occ.C <- as.data.frame(Occ.C)
df.Occ.C <- rownames_to_column(df.Occ.C, var = "OTU.ID")
df.Occ.C.tax <- merge(df.Occ.C, tax, by="OTU.ID")
sort.df.Occ.C.tax <- df.Occ.C.tax[order(df.Occ.C.tax$Occ.C, decreasing = TRUE),]
write.csv(sort.df.Occ.C.tax, file = "sort.df.Occ.C.tax.csv")




