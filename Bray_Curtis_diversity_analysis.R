######################################################################################
##### CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA USING BRAY CURTIS METHOD #####
######################################################################################

setwd('/Users/arifinabintarti/Documents/GitHub/Bean_seed_variability_Bintarti_2021/16S')
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis dissimilarity between samples)
otu_dist_bc <- vegdist(t(otu.norm), binary = F, method = "bray")
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa_bc <- cmdscale(otu_dist_bc, eig=T)
# scores of PC1 and PC2
ax1.scores.bc=otu_pcoa_bc$points[,1]
ax2.scores.bc=otu_pcoa_bc$points[,2] 
# calculate percent variance explained, then add to plot
ax1.bc <- otu_pcoa_bc$eig[1]/sum(otu_pcoa_bc$eig)
ax2.bc <- otu_pcoa_bc$eig[2]/sum(otu_pcoa_bc$eig)

#loading metadata
map <- read.csv("bean.var.map.csv")
head(map)
map <- column_to_rownames(map, var="Sample.id")
View(map)

map.pcoa.bc <- cbind(map,ax1.scores.bc,ax2.scores.bc)
map.pcoa.bc
# PCoA Plot 
set.seed(1)
pod.pcoa.bc <- ggplot(data = map.pcoa.bc, aes(x=ax1.scores.bc, y=ax2.scores.bc))+
            theme_bw()+
            geom_point(data = map.pcoa.bc, aes(x = ax1.scores.bc, y = ax2.scores.bc, col=factor(Plant)),size=5, alpha =0.7)+
            #scale_color_manual(labels = c("A1","A2", "A3","B1","B2","B3","B4","B5","B6","C5","C6","C7"),values=c("#440154FF", "#482677FF","#3F4788FF","#238A8DFF","#1F968BFF","#20A386FF","#29AF7FFF","#3CBB75FF","#56C667FF","#B8DE29FF","#DCE318FF","#FDE725FF"))+
            scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.bc,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste("PCoA2:\n",round(ax2.bc,3)*100,"% var. explained", sep=""))+
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

pod.pcoa2.bc <- pod.pcoa.bc + geom_text_repel(aes(label = Pod),size = 3, max.overlaps = Inf) 
pod.pcoa2.bc

######## Calculated the statistical analysis of beta diversity using nested permanova #########
#install.packages("BiodiversityR")
library(BiodiversityR)
map.pcoa.bc$Plant <- as.factor(map.pcoa.bc$Plant)
map.pcoa.bc$Pod <- as.factor(map.pcoa.bc$Pod)
set.seed(1)
nested.npmanova(unname(otu_dist_bc) ~ Plant + Pod, 
                data = map.pcoa.bc, 
                method = "bray", 
                permutations = 999)
###############################################################################################################################

## Betadisper grouped by Plant

set.seed(13)
groups.plant <- factor(c(rep("A",12),rep("B",24), rep("C",11)))
otu_dist_bc <- vegdist(t(otu.norm), binary = F, method = "bray")
mod.bc <- betadisper(otu_dist_bc, groups.plant)
mod.bc
plot(mod.bc)
mod.bc$distances
dispersion.bc <- as.data.frame(mod.bc$distance)
names(dispersion.bc)[names(dispersion.bc) == "mod.bc$distance"] <- "Dispersion"
#add dispersion index
dispersion.bc <- rownames_to_column(dispersion.bc, "Sample.id")
#join dispersion index to the map 
map.pcoa.bc.mod  = rownames_to_column(map,var = "Sample.id")
map.pcoa.bc.mod <- merge(map.pcoa.bc.mod, dispersion.bc, by="Sample.id", all = T)
map.pcoa.bc.mod


set.seed(1)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.bc <- permutest(mod.bc, permutations = 999, pairwise = T)
permod.bc # there is significant differences in dispersion between groups

hsd.bc=TukeyHSD(mod.bc) #which groups differ in relation to their variances
hsd.bc
plot(hsd.bc)

hsd.group.bc=hsd.bc$group
df.hsd.group.bc=as.data.frame(hsd.group.bc)
df.hsd.group.bc=rownames_to_column(df.hsd.group.bc, var = "Comparison")
names(df.hsd.group.bc)[names(df.hsd.group.bc) == "p adj"] <- "P.adj"
df.hsd.group.bc
# get the significant letter
detach(package:plyr)
library(dplyr)
dis.summ.plant.bc <- map.pcoa.bc.mod %>% group_by(Plant) %>% summarize(max.dis=max(Dispersion))

hsd.letter.bc = cldList(P.adj ~ Comparison,
         data = df.hsd.group.bc,
         threshold = 0.05)
names(hsd.letter.bc)[names(hsd.letter.bc) == "Group"] <- "Plant"

new.dis.sum.bc <- left_join(hsd.letter.bc,dis.summ.plant.bc,by='Plant')  
new.dis.sum.bc
#plot betadisper among plant
set.seed(1)
dis.plant.bc <- ggplot(map.pcoa.bc.mod, aes(x=Plant, y=Dispersion, fill=Plant))+
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
                                y_position = 0.85, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
                    geom_signif(comparisons = list(c("A", "C")), 
                                annotations = "***", textsize = 6,
                                y_position =0.895, tip_length = 0,
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
dis.plant.bc

## Betadisper grouped by Pod

groups.pod <- factor(c(rep("A1",4),rep("A2",4), rep("A3",4), rep("B1",4), rep("B2",4), rep("B3",4), rep("B4",4),rep("B5",4), rep("B6",4), rep("C5",3), rep("C6",4), rep("C7",4)))
mod.pod.bc <- betadisper(otu_dist_bc, groups.pod)
mod.pod.bc
mod.pod.bc$distances
dispersion.pod.bc <- as.data.frame(mod.pod.bc$distance)
names(dispersion.pod.bc)[names(dispersion.pod.bc) == "mod.pod.bc$distance"] <- "Dispersion"
#add dispersion index
dispersion.pod.bc <- rownames_to_column(dispersion.pod.bc, "Sample.id")
#join dispersion index to the map 
bean.map.bc <- merge(bean.map, dispersion.pod.bc, by="Sample.id", all = T)

boxplot(mod.pod.bc)
# Null hypothesis of no difference in dispersion between groups
anova(mod.pod.bc) # there is significant differences in dispersion between groups
# the variances among groups are not homogenous,
hsd.pod.bc=TukeyHSD(mod.pod.bc) #which groups differ in relation to their variances
hsd.pod.bc
plot(hsd.pod.bc)


######################################################################################
###### CALCULATE BETA DIVERSITY (PCoA PLOT) FOR FUNGI USING BRAY CURTIS METHOD #######
######################################################################################

# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
otu_dist.its_bc <- vegdist(t(fgnorm), binary = F, method = "bray")

# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa.its.bc <- cmdscale(otu_dist.its_bc, eig=T)

# scores of PC1 and PC2
ax1.scores.its.bc=otu_pcoa.its.bc$points[,1]
ax2.scores.its.bc=otu_pcoa.its.bc$points[,2] 

# calculate percent variance explained, then add to plot
ax1.its.bc <- otu_pcoa.its.bc$eig[1]/sum(otu_pcoa.its.bc$eig)
ax2.its.bc <- otu_pcoa.its.bc$eig[2]/sum(otu_pcoa.its.bc$eig)

setwd('/Users/arifinabintarti/Documents/GitHub/Bean_seed_variability_Bintarti_2021/ITS/')
its.map <- read.csv("bean.var.map.its.csv")
head(its.map)

its.map.pcoa.bc <- cbind(its.map,ax1.scores.its.bc,ax2.scores.its.bc)
its.map.pcoa.bc
# simple plot
pcoa_plot.its.bc <- plot(ax1.scores.its.bc, ax2.scores.its.bc, xlab=paste("PCoA1: ",round(ax1.its.bc,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.its.bc,3)*100,"% var. explained", sep=""))

# PCoA Plot 
require("ggrepel")
library(ggrepel)
library(viridis)
set.seed(13)

# Fig.3. Fungal PCoA Plot

set.seed(1)
pod.pcoa.its.bc <- ggplot(data = its.map.pcoa.bc, aes(x=ax1.scores.its.bc, y=ax2.scores.its.bc))+
            theme_bw()+
            geom_point(data = its.map.pcoa.bc, aes(x = ax1.scores.its.bc, y = ax2.scores.its.bc, col=factor(Plant)),size=5, alpha =0.7)+
            scale_color_manual(name = "Plant and Pod", labels = c("A (Pod A1:A3)", "B (Pod B1:B6)", "C (Pod C5:C7)"), values=c("#440154FF", "#287D8EFF","#FDE725FF"))+
            #scale_color_viridis(discrete = T) +
            scale_x_continuous(name=paste("PCoA1:\n",round(ax1.its.bc,3)*100,"% var. explained", sep=""))+
            scale_y_continuous(name=paste(round(ax2.its.bc,3)*100,"% var. explained", sep=""))+
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

pod.pcoa.its2.bc <- pod.pcoa.its.bc + geom_text_repel(aes(label = Pod),size = 3, max.overlaps = Inf) 
pod.pcoa.its2.bc

######## Calculated the statistical analysis of beta diversity using nested permanova #########

#install.packages("BiodiversityR")
library(BiodiversityR)

its.map.pcoa.bc$Plant <- as.factor(its.map.pcoa.bc$Plant)
its.map.pcoa.bc$Pod <- as.factor(its.map.pcoa.bc$Pod)

set.seed(1)
nested.npmanova(unname(otu_dist.its_bc) ~ Plant + Pod, 
                data = its.map.pcoa.bc, 
                method = "bray", 
                permutations = 999)

###############################################################################################################################

## Betadisper grouped by  plant

set.seed  (13)
groups.plant.its <- factor(c(rep("A",10),rep("B",19), rep("C",11)))
mod.its.bc <- betadisper(otu_dist.its_bc, groups.plant.its)
mod.its.bc
mod.its.bc$distances
dispersion.its.bc <- as.data.frame(mod.its.bc$distance)
names(dispersion.its.bc)[names(dispersion.its.bc) == "mod.its.bc$distance"] <- "Dispersion"
#add dispersion index
dispersion.its.bc <- rownames_to_column(dispersion.its.bc, "Sample.id")
dispersion.its.bc

#join dispersion index to the map 
its.map
its.map.bc.mod <- merge(its.map, dispersion.its.bc, by="Sample.id", all = T)
its.map.bc.mod


set.seed  (1)
#permutation-based test for multivariate homogeneity of group dispersion (variances)
permod.its.bc <- permutest(mod.its.bc, permutations = 999, pairwise = T)
permod.its.bc # there is marginal differences in dispersion between groups
# the variances among groups are not homogenous,

set.seed  (1)
hsd.its.bc=TukeyHSD(mod.its.bc) #which groups differ in relation to their variances
hsd.its.bc
plot(hsd.its.bc)

hsd.group.its.bc=hsd.its.bc$group
df.hsd.group.its.bc=as.data.frame(hsd.group.its.bc)
df.hsd.group.its.bc=rownames_to_column(df.hsd.group.its.bc, var = "Comparison")
names(df.hsd.group.its.bc)[names(df.hsd.group.its.bc) == "p adj"] <- "P.adj"
df.hsd.group.its.bc
# get the significant letter
detach(package:plyr)
library(dplyr)
dis.summ.plant.its.bc <- its.map.bc.mod %>% group_by(Plant) %>% summarize(max.dis=max(Dispersion))

hsd.letter.its.bc = cldList(P.adj ~ Comparison,
         data = df.hsd.group.its.bc,
         threshold = 0.05)
names(hsd.letter.its.bc)[names(hsd.letter.its.bc) == "Group"] <- "Plant"

new.dis.sum.its.bc <- left_join(hsd.letter.its.bc,dis.summ.plant.its.bc,by='Plant')  
new.dis.sum.its.bc

#plot betadisper among plant
library(viridis)
set.seed(1)
disperplot.its.bc <- ggplot(its.map.bc.mod, aes(x=Plant, y=Dispersion, fill=Plant))+
                    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.5, trim = F) +
                    scale_fill_viridis(discrete = T)+
                    geom_point(shape = 21,size=2, position = position_jitterdodge(),alpha=1)+
                    theme_bw()+
                    expand_limits(x = 0, y = 0)+
                    labs(title= "(d)", y="Dispersion")+
                    geom_signif(comparisons = list(c("B", "C")), 
                                annotations = "*", textsize = 6,
                                y_position = 0.85, tip_length = 0,
                                map_signif_level=TRUE, vjust = 0.5) +
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
disperplot.its.bc

# A non-significant result in betadisper is not necessarily related to a significant/non-significant result in adonis.

##### Supplementary Fig.2 Make Bacterial and Fungal PCoA Plots and Dispersion Plots in the Same Panel######

pod.pcoa2.bc
pod.pcoa.its2.bc
dis.plant.bc
disperplot.its.bc

setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures/NewFigures')

library(patchwork)
PCoA.Beta.bc <- (pod.pcoa2.bc | pod.pcoa.its2.bc ) / (dis.plant.bc | disperplot.its.bc)
PCoA.Beta.bc

PCoA.Beta2.bc <- patchwork::patchworkGrob(PCoA.Beta.bc)
PCoA.Beta3.bc <- gridExtra::grid.arrange(PCoA.Beta2.bc, bottom =textGrob("Plant", gp=gpar(fontsize=18,fontface='bold')))

ggsave("Fig.S2.eps",
       PCoA.Beta3.bc, device=cairo_ps,
       width = 10, height = 9, 
       units= "in", fallback_resolution = 600)



