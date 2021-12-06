#################################### Bean seed microbiomes variability (new fungal ITS analysis from Nejc) #####################################
##
# Date: FJanuary 16th 2021
# By : Ari Fina Bintarti
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

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/GitHub/Bean_seed_variability_Bintarti_2021/ITS/ITSolddata/')
wd <- print(getwd())
otu.its <- read.table('otu_table_ITS_UPARSE.txt', sep='\t', header=T, row.names = 1)
dim(otu.its) ## total otu= 87, otu table still has Mock, NC, and PC in the sample

#select only biological sample from otu table
otu.its.unfil <- otu.its[,1:45] #unselect Mock, NC, and PC from the otu table
dim(otu.its.unfil)
#otu.its.unfil <- column_to_rownames(otu.its.unfil, var = "OTUid")
sort(rowSums(otu.its.unfil, na.rm = FALSE, dims = 1), decreasing = F)

# remove OTUs that do not present in sample
otu.its1.unfil=otu.its.unfil[which(rowSums(otu.its.unfil) > 0),]
dim(otu.its1.unfil) # otu= 76
sort(rowSums(otu.its1.unfil, na.rm = FALSE, dims = 1), decreasing = F)

#read taxonomy
taxonomy.its = read.csv("ITSotu_taxonomy_RDP.csv", header=T)
dim(taxonomy.its) #83 otus, still contains plantae taxa and contaminant microbial taxa

#read the metadata
map.its <- read.csv("bean.var.map.its.csv")

#merge 'otu.its' with 'taxonomy.its' to have match otu and taxonomy table
otu.its1.unfil <- rownames_to_column(otu.its1.unfil, var = "OTUid")
otu.tax.its <- merge(otu.its1.unfil,taxonomy.its, by = "OTUid")
dim(otu.tax.its) #72 fungal otus before plant and microbial contaminants removal

#filter out the plantae taxa from the 'otu.tax.its' table
otu.filter <- filter(otu.tax.its, Kingdom != "Plantae")
dim(otu.filter) #59 otus, filter out 13 plantae taxa

#separate the otu table and taxonomy table
otu.tab <- data.frame(otu.filter[,c(1:46)])
dim(otu.tab)# 59otus,
tax.tab <- data.frame(otu.filter[,c(1,47:53)])
head(tax.tab)

#fungal otu
otu.fg1 <- column_to_rownames(otu.tab, var = "OTUid")
sort(rowSums(otu.fg1, na.rm = FALSE, dims = 1), decreasing = F)
head(otu.fg1)
dim(otu.fg1) #59 plant filtered otus; otu table before normalization using metagenomeSeq package and before decontamination

#otu table of the negative control
#otu.its <- column_to_rownames(otu.its, var = "OTUid")
otu.its.NC <- otu.its[,"NC",drop=FALSE]#only negative control
otu.its.NC

######################################################################################################################################
######################################################################################################################################


############## remove contaminant reads/otus from otu fg using microDecon package ##################

#install and load microDecon package
#devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

#merge otu.NC to otu bac and put otu.NC as the first column
otu.its.NC <- rownames_to_column(otu.its.NC, var = "OTUid")
otu.fg1 <- rownames_to_column(otu.fg1, var = "OTUid")
otu.fg1.NC <- merge(otu.its.NC, otu.fg1, by="OTUid")
otu.fg1.NC
dim(otu.fg1.NC)

#decontamination
otu.its.decon <- decon(data = otu.fg1.NC, numb.blanks = 1, numb.ind = 45, taxa = F)
otu.its.decon$OTUs.removed # remove 2 OTUs (OTU_2 and OTU_3)
otu.its.dec.table <- as.data.frame(otu.its.decon$decon.table)
dim(otu.its.dec.table) # otu = 57

#remove NC from the otu.dec.table 
otu.its.dec.table <- otu.its.dec.table[,c(1,3:47)]
otu.its.dec.table
dim(otu.its.dec.table)#there are 57 otus and 45 samples
head(otu.its.dec.table)
#use first column as a rownames
#otu.its.dec.table <- data.frame(otu.its.dec.table, row.names = 1)
#otu.its.dec.table
#head(sort(colSums(otu.its.dec.table, na.rm = FALSE, dims = 1), decreasing = T))
#head(sort(rowSums(otu.its.dec.table, na.rm = FALSE, dims = 1), decreasing = F))

#merge 'otu.its.dec.table' with 'tax.tab' to have match otu and taxonomy table
#tax.tab <- rownames_to_column(tax.tab, var = "OTUid")
otu.tax.fil <- merge(otu.its.dec.table,tax.tab, by = "OTUid")
dim(otu.tax.fil) #57 otus

#separate the taxonomy from the otu table
otu.tax.fil <- column_to_rownames(otu.tax.fil, var = "OTUid")
otu.fil <- data.frame(otu.tax.fil[,c(1:45)])
dim(otu.fil)
tax.fil <- data.frame(otu.tax.fil[,c(46:52)])

#write to  excel to edit the blank with the last taxonomic level possible 
#write.csv(tax.fil, file = "tax.fil.csv")
#make the taxonomic csv file "tax.fil.edited.csv"


######################################################################################################################################
######################################################################################################################################


############# Reads normalization using metagnomeSeq ####################
library(metagenomeSeq)

#Creating a MRexperiment object
head(otu.its.dec.table) #make sure you have "OTUid" as the first column
dim(otu.its.dec.table) # 57 OTUs
#write.table(otu.its.dec.table, file = "otu.its.dec.table.txt", col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")

#load count data (otu.dec.its)
otu.dec.its <- loadMeta("otu.its.dec.table.txt", sep = "\t")
otu.dec.its$counts
dim(otu.dec.its$counts) #[1] 57 45

#load taxonomy
tax.fil.edited = read.csv("tax.fil.edited.csv", sep=',', header=T)
#write.table(tax.fil.edited, file = "tax.fil.edited.txt", col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")
#loading taxonomy table
tax.dec.its <- read.delim("tax.fil.edited.txt", stringsAsFactors = FALSE, sep = "\t")
tax.dec.its <- tax.dec.its %>%
     remove_rownames() %>%
     column_to_rownames(var = 'OTUid')
dim(tax.dec.its) #[1] 57 otus  7
#loading metadata
its.map <- read.csv("bean.var.map.its.csv")
#write.table(its.map, file = "bean.var.map.its.txt",col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")
meta.its <- loadPhenoData("bean.var.map.its.txt", sep="\t",tran = TRUE)
meta.its
#create one MRexperiment object 
phenotypeData <- AnnotatedDataFrame(meta.its)
phenotypeData
OTUdata <- AnnotatedDataFrame(tax.dec.its)
OTUdata
#create model
fg.data <- newMRexperiment(otu.dec.its$counts, phenoData = phenotypeData, featureData = OTUdata)
fg.data

#normalising the data

#normalise the data to account for differences due to uneven sequencing depth
#metagenomeSeq uses Cumulative Sum Scaling (CSS) normalisation instead of rarefaction

#cumNormStatFast=Calculates the percentile for which to sum counts up to and scale by.
p <- cumNormStat(fg.data)
p
#cumNorm=Calculates each column’s quantile and calculates the sum up to and including that quantile.
fg.data.norm <- cumNorm(fg.data, p = p)
fg.data.norm
#export count matrix
fgnorm <- MRcounts(fg.data.norm, norm = TRUE, log = F)
fgnorm <- as.data.frame(fgnorm)
fgnorm
head(sort(colSums(fgnorm, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(fgnorm, na.rm = FALSE, dims = 1), decreasing = FALSE))
#write.table(fgnorm, "normalised_logged_fg_otu.txt", sep = "\t", quote = F)
dim(fgnorm) #there are 57 otus and 45 samples



######################################################################################################################################
######################################################################################################################################

### Rarefaction curves ######
# using GlobalPatterns

library(phyloseq)
library(graphics)


# 1. rarefaction curve for decontaminated non-normalized OTU table



# decontaminated non-normalized OTU table
otu.fil
dim(otu.fil)

# fungal taxonomy
head(tax.fil.edited)
dim(tax.fil.edited)
rownames(tax.fil.edited) <- rownames(otu.fil)
# make phyloseq otu table and taxonomy
otu.fil.phyl = otu_table(otu.fil, taxa_are_rows = TRUE)
tax.fil.phyl = tax_table(as.matrix(tax.fil.edited))

# make phyloseq map
rownames(its.map) <- its.map$Sample.id
its.map.phyl <- sample_data(its.map)

# make phyloseq object  
fil.obj <- merge_phyloseq(otu.fil.phyl,tax.fil.phyl,its.map.phyl)
fil.obj


# OR


#its.map.rare <- its.map[its.map$Sample.id%in%colnames(otu.fil),]
#curve_colors <- rep("darkgreen", ncol(otu.fil))
#curve_colors[its.map.rare$Plant=="A"] <- "#440154FF"
#curve_colors[its.map.rare$Plant=="B"] <- "#1F968BFF"
#curve_colors[its.map.rare$Plant=="C"] <- "#FDE725FF"

#lty <- c("solid", "solid", "solid")

#setEPS()
#postscript("Fig.1A Fungal rarecurve decontaminated non-normalized.eps", height = 5, width = 5)
#rarecurve(t(otu.fil), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", ylab = "Fungal OTUs", main="A", adj=0)
# Add a legend
#legend(4000, 5, legend=c("A", "B", "C"),col=c("#440154FF", "#1F968BFF","#FDE725FF"), lty = lty, lwd=2, title = "Plant")
#dev.off()
#graphics.off()


# 2. rarefaction curve for undecontaminated non-normalized OTU table

# read the otu table
#otu.fg1 # undecontaminated non-normalized OTU table
#otu.fg1 <- column_to_rownames(otu.fg1, var = "OTUid")

#setEPS()
#postscript("Fig.1B Fungal rarecurve undecontaminated non-normalized.eps", height = 5, width = 5)
#rarecurve(t(otu.fg1), step=1, label=F, col = curve_colors, lty=lty, lwd=2, xlab = "Reads", ylab = "Fungal OTUs", main="B", adj=0)
# Add a legend
#legend(7000, 6, legend=c("A", "B", "C"),col=c("#440154FF", "#1F968BFF","#FDE725FF"), lty = lty, lwd=2, title = "Plant")
#dev.off()

#set seed
set.seed(42)

#rarefy the data
# data = phyloseq object of decontaminated non normalized otu table
p.rare.its <- ggrare(fil.obj, step = 1, color = "Plant", label = "Sample", se = FALSE)

#set up your own color palette
Palette <- c("#440154FF","#1F968BFF","#FDE725FF")
names(Palette) <- levels(sample_data(fil.obj)$Plant)
Palette

#plot the rarecurve
#p <- ggrare(psdata, step = 1000, color = "SampleType", label = "Sample", se = FALSE)
p.its <- p.rare.its + 
 #facet_wrap(~Plant)+
 theme_bw()+
 scale_color_manual(values = Palette)+
 labs(title = "(b) Fungi")+
 expand_limits(x = 0, y = 0)+
 theme( strip.text.x = element_text(size=14, face='bold'),
        axis.text=element_text(size = 14),
        #axis.text.y = element_text(size = 13),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(size=15,face="bold"),
        #axis.title.y = element_markdown(),
        legend.position = "none",
        #legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        #xlab("Sequencing depth (Reads)") +  ylab("Fungal<br>OTU Number")
 

plot(p.its)
######################################################################################################################################
######################################################################################################################################


# Make Rank Abundance Curve 

library(phyloseq)

# load normalized otu table
fgnorm
dim(fgnorm)
head(fgnorm)
#otu.norm <- column_to_rownames(otu.norm, var = "OTUid")

# fungal taxonomy
head(tax.fil.edited)
tax.fil.edited <- column_to_rownames(tax.fil.edited, var = "OTUid")
dim(tax.fil.edited)
rownames(tax.fil.edited) <- rownames(fgnorm)
# make phyloseq otu table and taxonomy
otu.its.phyl = otu_table(fgnorm, taxa_are_rows = TRUE)
tax.its.phyl = tax_table(as.matrix(tax.fil.edited))

# make phyloseq map
rownames(its.map) <- its.map$Sample.id
its.map.phyl <- sample_data(its.map)

# make phyloseq object  

its.phyl.obj <- merge_phyloseq(otu.its.phyl,tax.its.phyl,its.map.phyl)
its.phyl.obj

# this converts taxa counts in each sample to a percentage
phyloTemp2.its <-  transform_sample_counts(its.phyl.obj, function(x) x/sum(x))
clusterData2.its <-  psmelt(phyloTemp2.its)

# calculating mean relative abundance per plant
relabund.perplant2.its= clusterData2.its %>%
  group_by(Plant,Pod, OTU) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(desc(ra),.by_group = TRUE) %>%
  mutate(rank = row_number(-ra))

relabund.perplant2.its$Plant <- as.factor(relabund.perplant2.its$Plant)
relabund.perplant2.its$Pod <- as.factor(relabund.perplant2.its$Pod)

#set up your own color palette
col.pal.its <- c("#440154FF","#1F968BFF","#FDE725FF")
names(col.pal.its) <- levels(relabund.perplant2.its$Plant)
col.pal.its

# plotting
rank.abund2.its <- ggplot(relabund.perplant2.its,aes(x=rank,y=ra, colour=Plant)) +
  geom_line(aes(group = Pod), size=1) +
  scale_colour_manual(values = col.pal.its)+
  labs(title = "(d)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = 14),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size =20 ,face='bold'),
        #axis.title = element_text(size=14,face="bold"),
        legend.position = "none",
        #legend.title = element_text(size=15),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        #xlab("\nRank") +  ylab("Fungal<br>Relative Abundance")
        #guides(colour=guide_legend(title.position = "top", keywidth = 2, nrow = 1, byrow = T,title.hjust =0.5,
                                   #title.theme = element_text(size=13),
                                   #label.theme = element_text(size = 12, face = 'bold')))
rank.abund2.its


data.rank.abund2.its <- ggplot_build(rank.abund2.its)
gtable.rank.abund2.its <- ggplot_gtable(data.rank.abund2.its)


ggsave("Fig.2.Fungal rankabund curve.eps",
       rank.abund2.its, device = "eps",
       width = 5, height =4.5, 
       units= "in", dpi = 600)


library(ggplot2)
library(scales)



######################################################################################################################################
######################################################################################################################################

###### calculate alpha diversity  ##########

# calculate richness
s <- specnumber(fgnorm, MARGIN = 2) # richness
rich <- as.data.frame(s)
# calculate shannon index
h <- diversity(t(fgnorm), index = 'shannon') # Shannon index
shannon <- as.data.frame(h)
# calculate evenness
j <- h/log(s) # Pielou's evenness
pie <- as.data.frame(j)
# gather the alpha diversity into the metadata
its.map$Richness <- s
its.map$Shannon <- h
its.map$Pielou <- j


################################################################################################################################################
###############################################################################################################################################

### Compare Fungal Richness among plants and pods ###

# 1. Model fitting with linear mixed effect modelling with lme() function 

its.map$Plant <- as.factor(its.map$Plant)
its.map$Pod <- as.factor(its.map$Pod)

set.seed(13)
library(nlme)
# Fit random intercept model
model.its = lme(Richness ~ Plant, random = ~ 1|Pod, data=its.map, method="REML")
summary(model.its)
aov.mod.its = anova.lme(model.its, 
          type="marginal", 
          adjustSigma = FALSE)
aov.mod.its
anova(model.its, type='marginal')
# Result: plant does not has significant effect to the fungal richness, 0.3713

# Test the significance of the random effect in the mixed effects model
model.fixed.its = gls(Richness ~ Plant, 
                  data=its.map, 
                 method="REML")
anova(model.its, 
      model.fixed.its) # pod also has no significant effect to the fungal richness

# Checking assumptions of the model
hist(residuals(model.its), 
     col="darkgray") # the model is almost normal
# Using the aov function for a nested anova
fit.its.x = aov(Richness ~ Plant + Error(Pod), data=its.map)
summary(fit.its.x)

# Using Mean Sq and Df values to get p-value for H = plant and Error = pod
pf(q=4.007/3.502,
   df1=2,
   df2=9,
   lower.tail=FALSE)
# Using Mean Sq and Df values to get p-value for H = plant and Error = Within
summary(fit)
pf(q=3.502/2.657,
   df1=9,
   df2=33,
   lower.tail=F)

# model evaluation

# normality
plot(model.its)
qqnorm(resid(model.its))
qqline(resid(model.its)) #There is some deviation from from the expected normal line towards the tails, but overall the line looks straight and therefore pretty normal and suggests that the assumption is not violated
library(sjPlot)
plot_grid(plot_model(model.its, type = "diag"))
kurtosis(resid(model.its),method = 'sample')
# homoscedasticity
leveneTest(residuals(model.its) ~ its.map$Plant) #Since the p value is greater than 0.05, we can say that the variance of the residuals is equal and therefore the assumption of homoscedasticity is met
boxplot(residuals(model.its) ~ its.map$Plant) 

#1. plot fungal richness among plants

library(viridis)
library(ggtext)
#Fig.1A
set.seed(1)
fg.rich.plant <- ggplot(its.map, aes(x=Plant, y=Richness, fill=Plant))+
                    #geom_boxplot() +
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    scale_fill_viridis(discrete = T)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "(e)")+
                    ylab("Fungal<br>Richness")+
                    xlab("\nPlant")+
                    #geom_text(data=sum_rich_plant_new, aes(x=plant,y=2+max.rich,label=Letter), vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20, face = 'bold'),
                          axis.title.x=element_text(size=15,face="bold"),
                          axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=13, color="red", shape=95)
fg.rich.plant

#2. plot fungal richness among pods

set.seed(1)
fg.rich.pod <- ggplot(its.map, aes(x=Pod, y=Richness, fill = Pod))+
                    #geom_boxplot() +
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    scale_fill_viridis(discrete = T)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title = "(f)")+
                    xlab("\nPod")+
                    #geom_text(data=new.rich_pod.summarized,aes(x=pod,y=0.5+max.rich,label=new.rich_pod.summarized$groups),vjust=0)+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_blank(),
                          strip.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          plot.title = element_text(size = 20,face =  'bold'),
                          axis.title.x =element_text(size=15,face="bold"),
                          axis.title.y = element_blank(),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=6, color="red", shape=95)
fg.rich.pod
fg.rich.pod1 <- fg.rich.pod +
  facet_grid(. ~ Plant, scales="free_x")+
  theme(strip.text = element_blank())
  #theme(strip.background =element_rect(fill="grey"))+
  #theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
fg.rich.pod1

##################################### FUNGAL COMPOSITION ####################################################################

#BiocManager::install("phyloseq")
library(phyloseq)

# load normalized otu table
fgnorm
dim(fgnorm)
head(fgnorm)
#otu.norm <- column_to_rownames(otu.norm, var = "OTUid")

# fungal taxonomy
head(tax.fil.edited)
tax.fil.edited <- column_to_rownames(tax.fil.edited, var = "OTUid")
dim(tax.fil.edited)
rownames(tax.fil.edited) <- rownames(fgnorm)
# make phyloseq otu table and taxonomy
otu.its.phyl = otu_table(fgnorm, taxa_are_rows = TRUE)
tax.its.phyl = tax_table(as.matrix(tax.fil.edited))

# make phyloseq map
rownames(its.map) <- its.map$Sample.id
its.map.phyl <- sample_data(its.map)

# make phyloseq object

its.phyl.obj <- merge_phyloseq(otu.its.phyl,tax.its.phyl,its.map.phyl)
its.phyl.obj

# merge taxa by phylum

# 1. phylum - Fungi
fg.phylum <- tax_glom(its.phyl.obj, taxrank = "Phylum", NArm = F)
fg.phylum.ra <- transform_sample_counts(fg.phylum, function(x) x/sum(x))
fg.phylum.ra

df.fg.phylum <- psmelt(fg.phylum.ra) %>%
  group_by(Sample.id,Plant,Pod,Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

# OR

###create dataframe from phyloseq object
library(data.table)
df.fg.phyl <- data.table(psmelt(fg.phylum.ra))
dim(df.fg.phyl)
###convert phylum to character vector from a factor
df.fg.phyl$Phylum <- as.character(df.fg.phyl$Phylum)

# calculating mean relative abundance per plant
fg.mean.phylrelabund.perplant= df.fg.phyl %>%
  group_by(Plant, Phylum) %>%
  summarise(ra=mean(Abundance))

# barplot of fungal composition across pods at Phylum level

library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
my_colors = carto_pal(12, "Safe")
my_colors

# New facet label names for plant variable
plant.labs <- c("Plant A", "Plant B", "Plant C")
names(plant.labs) <- c("A", "B", "C")

# Create the plot

fg.pod.phylum <- ggplot(data=df.fg.phylum, aes(x=Pod, y=Mean, fill=Phylum))
plot.fg.pod.phylum <- fg.pod.phylum + 
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
                           
plot.fg.pod.phylum

plot.fg.pod.phylum1 <- plot.fg.pod.phylum +
  facet_grid(. ~ Plant, labeller = labeller(Plant=plant.labs), scales="free_x")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black', size = 14, face = 'bold'))
plot.fg.pod.phylum1


# Merge taxa by Class

# 2. CLass - Fungi
fg.cl <- tax_glom(its.phyl.obj, taxrank = "Class", NArm = F)
fg.cl.ra <- transform_sample_counts(fg.cl, function(x) x/sum(x))
fg.cl.ra

df.fg.cl <- psmelt(fg.cl.ra) %>%
  group_by(Sample.id, Plant, Pod, Class) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.fg.cl$Class <- as.character(df.fg.cl$Class)
df.fg.cl$Class[df.fg.cl$Mean < 0.1] <- "Other"

# Create the plot

#Fig.2. Fungal Barplot

fg.pod.cl <- ggplot(data=df.fg.cl, aes(x=Pod, y=Mean, fill=Class))
plot.fg.pod.cl <- fg.pod.cl + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     scale_fill_manual(values=c("#44AA99", "#332288", "#117733","#CC6677","#DDCC77", "#88CCEE","#661100","#AA4499" ,"#888888",'tomato','#ffd8b1'))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(y= "Mean Relative Abundance", x="Pod", title = "(b) Fungi")+
                     theme(plot.title = element_text(size = 20, face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=14),
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
                           guides(fill=guide_legend(ncol=1,bycol=TRUE))
                           
plot.fg.pod.cl

plot.fg.pod.cl1<- plot.fg.pod.cl +
  facet_grid(. ~ Plant, scales="free_x")+
  theme(strip.text.x = element_blank())
plot.fg.pod.cl1


pd.pod <- pd.p +
  facet_grid(. ~ Plant, scales="free_x")+
  theme(strip.text = element_blank())


ggsave("Fig.4.Fungal relative abund.eps",
      plot.fg.pod.cl1, device = "eps",
       width = 8.6, height =5, 
       units= "in", dpi = 600)

fg.genus <- tax_glom(its.phyl.obj, taxrank = "Genus", NArm = F)
fg.genus.ra <- transform_sample_counts(fg.genus, function(x) x/sum(x))
fg.genus.ra
###create dataframe from phyloseq object
library(data.table)
df.fg.genus <- data.table(psmelt(fg.genus.ra))
dim(df.fg.genus)
# calculating mean relative abundance per genus
mean.its.relabund.genus= df.fg.genus %>%
  group_by(Genus) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(-ra)


df.fg.cl <- data.table(psmelt(fg.cl.ra))
mean.its.relabund.cl= df.fg.cl %>%
  group_by(Class) %>%
  summarise(ra=mean(Abundance)) %>%
  arrange(-ra)


# other: Fungal classes with mean relative abundances of less than 10 %

##############################################################################################################################

# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR FUNGI

# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
fgnorm_PA <- 1*(fgnorm>0)
fgnorm_PA
otu_dist.its <- vegdist(t(fgnorm_PA), binary = TRUE, method = "jaccard") #Sorensen

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa.its <- cmdscale(otu_dist.its, eig=T)

# scores of PC1 and PC2
ax1.scores.its=otu_pcoa.its$points[,1]
ax2.scores.its=otu_pcoa.its$points[,2] 

# calculate percent variance explained, then add to plot
ax1.its <- otu_pcoa.its$eig[1]/sum(otu_pcoa.its$eig)
ax2.its <- otu_pcoa.its$eig[2]/sum(otu_pcoa.its$eig)
its.map.pcoa <- cbind(its.map,ax1.scores.its,ax2.scores.its)

# simple plot
pcoa_plot.its <- plot(ax1.scores.its, ax2.scores.its, xlab=paste("PCoA1: ",round(ax1.its,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.its,3)*100,"% var. explained", sep=""))

# PCoA Plot 
require("ggrepel")
library(ggrepel)
library(viridis)
set.seed(13)

# Fig.3. Fungal PCoA Plot
set.seed(13)
pod.pcoa.its <- ggplot(data = its.map.pcoa, aes(x=ax1.scores.its, y=ax2.scores.its))+
            theme_bw()+
            geom_point(data = its.map.pcoa, aes(x = ax1.scores.its, y = ax2.scores.its, col=factor(Plant)),size=5, alpha =0.7)+
            scale_color_manual(name = "Plant and Pod", labels = c("A (Pod A1:A3)", "B (Pod B1:B6)", "C (Pod C5:C7)"), values=c("#440154FF", "#287D8EFF","#FDE725FF"))+
            #scale_color_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FFF","#3CBB75FF","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
            #scale_color_viridis(discrete = T, ) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.its,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste(round(ax2.its,3)*100,"% var. explained", sep=""))+
            #coord_fixed() + 
            labs(colour = "Plant", title = "(b) Fungi")+
            theme(legend.position="right",
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 20, face="bold"),
            axis.text=element_text(size=14), 
            axis.title=element_text(size=15,face="bold"),
            legend.text=element_text(size=14),
            legend.title = element_text(size = 14, face = 'bold'),
            legend.spacing.x = unit(0.05, 'cm'))

set.seed(13)
pod.pcoa.its2 <- pod.pcoa.its + geom_text_repel(aes(label = Pod),size = 3, max.overlaps = Inf) 
pod.pcoa.its2

#ggsave("Fig.5.Fungal PCoA plot tiff",
       #pod.pcoa.its2, device = "tiff",
       #width = 5, height =4, 
       #units= "in", dpi = 600)



######## Calculated the statistical analysis of beta diversity using nested permanova #########

set.seed(13)
adonis.its <- adonis(otu_dist.its ~ its.map$Plant/its.map$Pod, 
                 permutation=999,
                 method="jaccard", 
                 strata = NULL)
adonis.its
#install.packages("BiodiversityR")
library(BiodiversityR)

its.map$Plant <- as.factor(its.map$Plant)
its.map$Pod <- as.factor(its.map$Pod)

set.seed(13)
nested.npmanova(otu_dist.its ~ Plant + Pod, 
                data = its.map, 
                method = "jac", 
                permutations = 999)

###############################################################################################################################

## Betadisper
set.seed  (13)
groups.plant.its <- factor(c(rep("A",11),rep("B",23), rep("C",11)))
otu_dist.its <- vegdist(t(fgnorm_PA), binary = TRUE, method = "jaccard") #Sorensen
mod.its <- betadisper(otu_dist.its, groups.plant.its)
mod.its
mod.its$distances
dispersion.its <- as.data.frame(mod.its$distance)
names(dispersion.its)[names(dispersion.its) == "mod.its$distance"] <- "Dispersion"
#add dispersion index
dispersion.its <- rownames_to_column(dispersion.its, "Sample.id")
#join dispersion index to the map 
its.map <- merge(its.map, dispersion.its, by="Sample.id", all = T)

set.seed  (13)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.its <- permutest(mod.its, permutations = 999, pairwise = T)
permod.its # there is no significant differences in dispersion between groups
# the variances among groups are homogenous,

#plot betadisper among plant

# Fig. 4 Betadispersion Plot
library(viridis)
set.seed(1)
disperplot.its <- ggplot(its.map, aes(x=Plant, y=Dispersion, fill=Plant))+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    #scale_fill_manual(labels = c("A", "B", "C"),values=c("#CC6677", "#DDCC77","#117733"))+
                    scale_fill_viridis(discrete = T)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    #geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title= "(d)", y="Dispersion")+
                    theme(legend.position="none",
                          axis.text.x=element_text(size = 14),
                          axis.text.y = element_text(size = 14),
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20, face = 'bold'),
                          #axis.title.x =element_text(size=18,face="bold"),
                          axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          stat_summary(fun="median",geom="point", size=13, color="red", shape=95)
disperplot.its
# save the plot
ggsave("Fig.6.Fungal betadisper.tiff",
       disperplot.its, device = "tiff",
       width = 5, height =4.5, 
       units= "in", dpi = 600)
# A non-significant result in betadisper is not necessarily related to a significant/non-significant result in adonis.

#####################################################################################################################

######################################################################################################################################################
######################################################################################################################################################
# Read Proportion of Chloroplast and mitochondria 


#read the unfiltered otu table and tax table
otu.tax.its

#select only the otu table and "Order"  & "Family"
otu.its <- otu.tax.its[,c(1:46)]
otu.its <- column_to_rownames(otu.its, var = "OTUid")
sum(otu.its)
otu.its.unfil <- otu.tax.its[,c(1:47)]
colnames(otu.its.unfil)

#edit the taxonomy
otu.its.unfil.ed <- otu.its.unfil %>%
    mutate(Taxonomy = case_when(Kingdom == "Plantae" ~ 'Plant',
                                TRUE ~ 'Fungi'))

long.dat.its <- gather(otu.its.unfil.ed, Sample, Read, A12:C74, factor_key = T)
long.dat.its


df.unfil.its <- long.dat.its %>%
  group_by(Sample, Taxonomy) %>%
  summarize(read.number = sum(Read))
df.unfil.its1 <- df.unfil.its %>% 
  group_by(Sample) %>% 
  mutate(percent= prop.table(read.number) * 100)

with(df.unfil.its1, sum(percent[Sample ==  "A12"]))

plot.unfil.its <- ggplot(df.unfil.its1, aes(x=Taxonomy, y=percent, fill=Taxonomy))+
                    geom_violin(trim = F, scale="width") +
                    #scale_fill_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FF","#3CBC75F","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
                    scale_fill_manual(labels = c("Fungi","Plant"),values=c("#88CCEE", "#117733"))+
                    geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    #geom_text(data=sum_rich_plant_new, aes(x=Plant,y=2+max.rich,label=Letter), vjust=0)+
                    labs(title = "B")+
                    ylab("Read Proportion (%)")+
                    theme(legend.position="none",
                          axis.text.y=element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x= element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 14, face = 'bold'),
                          axis.title.y=element_blank(),
                          #axis.title.y = element_markdown(size=15,face="bold"),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
                          #plot.margin = unit(c(0, 0, 0, 0), "cm"))
                          stat_summary(fun="median",geom="point", size=7, color="red", shape=95)
plot.unfil.its

setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
ggsave("plant.its.tiff",
       plot.unfil.its, device = "tiff",
       width = 5, height =5, 
       units= "in", dpi = 600)

#####################################################################################################################################
######################################################################################################################################

#### 1.  Occupancy-mean relative abundance across metadata #####
setwd('/Users/arifinabintarti/Documents/GitHub/Bean_seed_variability_Bintarti_2021/ITS/')
wd <- print(getwd())

# load normalized otu table
fgnorm
dim(fgnorm)
head(fgnorm)
#otu.norm <- column_to_rownames(otu.norm, var = "OTUid")

# load fungal taxonomy
head(tax.fil.edited)
tax.fil.edited <- rownames_to_column(tax.fil.edited, var = "OTUid")
dim(tax.fil.edited)

# load metadata
its.map <- read.csv("bean.var.map.its.csv")
head(its.map)

##calculate the occupancy of each OTUid across Sample.id/Seed
fgnorm.seed.PA <- 1*((fgnorm>0)==1)
Occ.seed.its <- rowSums(fgnorm.seed.PA)/ncol(fgnorm.seed.PA)
df.Occ.seed.its <- as.data.frame(Occ.seed.its)
df.Occ.seed.its <- rownames_to_column(df.Occ.seed.its, var = "OTUid")
df.Occ.seed.its.tax <- merge(df.Occ.seed.its, tax.fil.edited, by="OTUid")
sort.df.Occ.seed.its.tax <- df.Occ.seed.its.tax[order(df.Occ.seed.its.tax$Occ.seed.its, decreasing = TRUE),]
#write.csv(sort.df.Occ.seed.its.tax, file = "sort.df.Occ.seed.its.tax.csv")

### 2. Occupancy across plants

longDF.plant.its <- data.frame(OTUid=as.factor(rownames(fgnorm)), fgnorm) %>%
  gather(Sample.id, abun, -OTUid) %>%  #keep same column nameing as in mapping file, calling counts as "abun" (abundance)
  left_join(its.map) %>%  #will add the info form mapping file (grouped by the 'Sample.id' column)
  left_join(tax.fil.edited)  %>% #adding the taxonomy info (grouped by the 'OTU.ID' column)
  group_by(OTUid, Plant) %>%
  #BioProject, plant_family, plant_genus, host, compartment
  summarise(n=sum(abun))
  #mutate(n=if_else(n>0, 1,0))

##build the new table: OTUid as rownames and Plant as colnames
wideDF.plant.its <- as.data.frame(spread(longDF.plant.its, OTUid, n, fill=0))
rownames(wideDF.plant.its) <-  wideDF.plant.its[,1]
wideDF.plant.its <- wideDF.plant.its[,-1]
wideDF.plant.its <- t(wideDF.plant.its)
dim(wideDF.plant.its) #57 taxa, 3 plant 

##calculate the occupancy of each otu across plants
wideDF.plant.PA.its <- 1*((wideDF.plant.its>0)==1)
Occ.plant.its <- rowSums(wideDF.plant.PA.its)/ncol(wideDF.plant.PA.its)
df.Occ.plant.its <- as.data.frame(Occ.plant.its)
df.Occ.plant.its <- rownames_to_column(df.Occ.plant.its, var = "OTUid")
df.Occ.plant.tax.its <- merge(df.Occ.plant.its, tax.fil.edited, by="OTUid")
sort.df.Occ.plant.tax.its <- df.Occ.plant.tax.its[order(df.Occ.plant.tax.its$Occ.plant.its, decreasing = TRUE),]
write.csv(sort.df.Occ.plant.tax.its, file = "sort.df.Occ.plant.tax.its.csv")

### 3. Occupancy across seeds within plant

# 1.Plant A
fgnormA <- data.frame(fgnorm[,c(1:11)])
##calculate the occupancy of each otu across seeds within plant
fgnormA.PA <- 1*((fgnormA>0)==1)
fgOcc.A <- rowSums(fgnormA.PA)/ncol(fgnormA.PA)
df.fgOcc.A <- as.data.frame(fgOcc.A)
df.fgOcc.A <- rownames_to_column(df.fgOcc.A, var = "OTUid")
df.fgOcc.A.tax <- merge(df.fgOcc.A, tax.fil.edited, by="OTUid")
sort.df.fgOcc.A.tax <- df.fgOcc.A.tax[order(df.fgOcc.A.tax$fgOcc.A, decreasing = TRUE),]
write.csv(sort.df.fgOcc.A.tax, file = "sort.df.fgOcc.A.tax.csv")

# 2.Plant B
fgnormB <- data.frame(fgnorm[,c(12:34)])
##calculate the occupancy of each otu across seeds within plant
fgnormB.PA <- 1*((fgnormB>0)==1)
fgOcc.B <- rowSums(fgnormB.PA)/ncol(fgnormB.PA)
df.fgOcc.B <- as.data.frame(fgOcc.B)
df.fgOcc.B <- rownames_to_column(df.fgOcc.B, var = "OTUid")
df.fgOcc.B.tax <- merge(df.fgOcc.B, tax.fil.edited, by="OTUid")
sort.df.fgOcc.B.tax <- df.fgOcc.B.tax[order(df.fgOcc.B.tax$fgOcc.B, decreasing = TRUE),]
write.csv(sort.df.fgOcc.B.tax, file = "sort.df.fgOcc.B.tax.csv")

# 3.Plant C
fgnormC <- data.frame(fgnorm[,c(35:45)])
##calculate the occupancy of each otu across seeds within plant
fgnormC.PA <- 1*((fgnormC>0)==1)
fgOcc.C <- rowSums(fgnormC.PA)/ncol(fgnormC.PA)
df.fgOcc.C <- as.data.frame(fgOcc.C)
df.fgOcc.C <- rownames_to_column(df.fgOcc.C, var = "OTUid")
df.fgOcc.C.tax <- merge(df.fgOcc.C, tax.fil.edited, by="OTUid")
sort.df.fgOcc.C.tax <- df.fgOcc.C.tax[order(df.fgOcc.C.tax$fgOcc.C, decreasing = TRUE),]
write.csv(sort.df.fgOcc.C.tax, file = "sort.df.fgOcc.C.tax.csv")

