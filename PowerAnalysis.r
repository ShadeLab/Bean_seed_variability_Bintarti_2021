################################################################################################################################
################################  Determining Power for Bean Drought Experiment   #######################################

install.packages("paramtest")
install.packages("pwr")
install.packages("knitr")
install.packages("nlme")
install.packages("lavaan")

library('paramtest')
library('pwr')
library('ggplot2')
library('knitr')
library('nlme')
library('lavaan')
library('dplyr')

# We want to determine how many samples needed to allow us to detecan effect of given size with a given degree of confidence
# there will be control and drought treatments, so we want to see the effect of drought to seed microbiome structure
# we will use t-test for the analysis

pwr.t.test(n = , d = , sig.level = , power = , type = c("two.sample", "one.sample", "paired"))

# d = effect size
# n = sample size

#load the  data
bean.map
data.rich <- bean.map[,c(2,6,7)]

#calculate mean and sd of richness and PD
data.bac <- bean.map %>%
 summarise(mean.rich = mean(Richness),
          mean.pd = mean(PD_whole_tree), 
          sd.rich = sd(Richness),
          sd.pd = sd(PD_whole_tree))

#group by plant
data.bac.plant <- bean.map %>%
 group_by(Plant) %>%
summarise(mean.rich = mean(Richness),
          sd.rich = sd(Richness),
          var.rich = var(Richness),
          mean.pd = mean(PD_whole_tree), 
          sd.pd = sd(PD_whole_tree),
          var.pd = var(PD_whole_tree))

data.bac.plant

###  do some simulation for t-test two means / t-test two samples for richness ###

#plant A
n1 <- 12
mean.rich.a <- 30.58
sd.rich.a <- 6.41

#plant B
n2 <- 24
mean.rich.b <- 18.20
sd.rich.b <- 7.35

#add size to the data frame
data.bac.plant$Size <- c(12,24,11)

install.packages("utilities")
library(utilities)

# calculating pooled standard deviation
#sample.decomp(n = data.bac.plant$Size, sample.mean = data.bac.plant$mean.rich, sample.sd = data.bac.plant$sd.rich, include.sd = TRUE)
#sd.rich.pool <- 9.568392

#calculating delta manually
#delta <- (mean.rich.a-mean.rich.b)/sd.rich.pool
#delta

# Calculating effect size (Cohen's d) using the max and min means out of the 3 groups using linear model
group <- 3
data.bac.plant

fit.rich <- lm(Richness~Plant, data = bean.map)
anova(fit.rich)

# Error or within-group variance or Means Squared Error
anova(fit.rich)["Residuals", "Mean Sq"] #[1] 65.85873 or MSE

# Calculating effect size (Cohen's d)
# d= mean.largest-mean.smallest / sqrt(MSE) 

d=((30.58333-18.20833))/(sqrt(65.85873))
d # 1.52489

(30.58333-18.20833)/8.114

# sqrt(MSE)  OR RMSE also can be calculated  below:
sqrt(sum(fit.rich$residuals^2) / fit.rich$df)

# corrected Cohen's  (Hedges's g)
g = d * (1-(3/((4*36)-9)))
g




# simulation
n <- seq(0, 100, 1)
n

power <- sapply(seq_along(n), function(i)
 power.t.test(n = n[i], delta =g, sig.level = 0.05, type = 'two.sample')$power)
power

power.df <- data.frame(n=n, power=power)
power.rich.plot <- ggplot(power.df, aes(x=n, y=power))+
 geom_line(size=1.5)+
 #geom_hline(yintercept = 0.8, linetype = 2, color = "gray30")+
 geom_hline(yintercept = 0.8, linetype = 2, color = "gray30")+
 geom_vline(xintercept = 10, linetype = 2, color = "gray30")+
 scale_x_continuous("Sample Size", breaks = seq(0,100, 20))+
 scale_y_continuous("Power", breaks = seq(0,1,.2))+
 labs(title="(a)")+
 expand_limits(x = 0, y = 0)+
 theme_bw(base_size = 14)+
 theme(axis.text = element_text(size = 14),
       axis.title.x = element_blank(),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20,face = 'bold'),
                          axis.title.y = element_text(size=15,face="bold"))

power.rich.plot

pwr.t.test(d=g,power=.95,sig.level=.05,type="two.sample")

#pwrt<-pwr.t.test(d=d,n=c(2,5,10,15,20,25,30,40,50,60,70,80,90,100),sig.level=.05,type="two.sample")
#pwrt
#plot(pwrt$n,pwrt$power,type="b",xlab="sample size",ylab="power")




###  do some simulation for t-test two means / t-test two samples for PD ###



#plant A
n1 <- 12
mean.pd.a <- 4.17
sd.pd.a <- 0.88

#plant B
n2 <- 24
mean.pd.b <- 2.92
sd.pd.b <- 0.82

data.bac.plant$Size <- c(12,24,11)

# Calculating effect size (Cohen's d) using the max and min means out of the 3 groups using linear model

fit.pd <- lm(PD_whole_tree~Plant, data = bean.map)
anova(fit.pd)

# Error or within-group variance or Means Squared Error
anova(fit.pd)["Residuals", "Mean Sq"] #[1] 0.9937546 or MSE

# Calculating effect size
# d= mean.largest-mean.smallest / sqrt(MSE) 

d.pd=((4.174887-2.924510))/(sqrt(0.9937546))
d.pd # 1.2543

# corrected Cohen's  (Hedges's g)
g.pd = d.pd * (1-(3/((4*36)-9)))
g.pd


n <- seq(0, 100, 1)
n

power.pd <- sapply(seq_along(n), function(i)
 power.t.test(n = n[i], delta =g.pd, sig.level = 0.05, type = 'two.sample')$power)
power.pd

power.pd.df <- data.frame(n=n, power=power.pd)
power.pd.plot <- ggplot(power.pd.df, aes(x=n, y=power))+
 geom_line(size=1.5)+
 geom_hline(yintercept = 0.8, linetype = 2, color = "gray30")+
 geom_vline(xintercept = 10, linetype = 2, color = "gray30")+
 scale_x_continuous("Sample Size", breaks = seq(0,100, 20))+
 scale_y_continuous("Power", breaks = seq(0,1,.2))+
 expand_limits(x = 0, y = 0)+
 #scale_x_discrete(labels= label)+
 labs(title="(b)")+
 theme_bw(base_size = 14)+
 theme(axis.text.y = element_blank(),
       axis.ticks.y = element_blank(),
       axis.title = element_blank(),
       axis.text.x = element_text(size = 14),
                          strip.text = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = 20,face = 'bold'))
                          #axis.title = element_text(size=15,face="bold"))
power.pd.plot 

pwr.t.test(d=g.pd,power=.8,sig.level=.05,type="two.sample")
pwr.t.test(d=g.pd,power=.95,sig.level=.05,type="two.sample")


#pwrt.pd<-pwr.t.test(d=d.pd,n=c(2,5,10,15,20,25,30,40,50,60,70,80,90,100),sig.level=.05,type="two.sample")
#pwrt.pd
#plot(pwrt.pd$n,pwrt.pd$power,type="b",xlab="sample size",ylab="power")

effect_sizes <- c(.25, .5, .8)
conditions <- expand.grid(n = n, effect_sizes = effect_sizes)

power_curve <- sapply(seq_len(nrow(conditions)), function(i) 
 power.t.test(n = conditions[i, 'n'],
              delta =  conditions[i, 'effect_sizes'], 
              type = 'two.sample')$power)
power_curve

power_curve_df <- bind_cols( conditions,power = power_curve)
ggplot(power_curve_df, aes(x = n, y = power)) +
 geom_line(aes(color = factor(effect_sizes)), size = 2) + 
 geom_hline(yintercept = 0.8, linetype = 2, color = 'gray30') + 
 scale_x_continuous("Sample Size", breaks = seq(0, 100, 20)) +
 scale_y_continuous("Power", breaks = seq(0, 1, .2)) + 
 scale_color_grey("Effect Size") +
 theme_bw(base_size = 14)


#######################################################################################################################################################
#######################################################################################################################################################

groupmeans <- data.bac.plant$mean.rich
p.rich <- power.anova.test(groups = length(groupmeans), 
between.var = var(groupmeans), within.var = 65.85873 , 
power=0.95,sig.level=0.05,n=NULL)
p.rich

groupmeans <- data.bac.plant$mean.rich
n.aov <- c(seq(2,10,by=1),seq(12,20,by=2),seq(25,50,by=5))
p.rich <- power.anova.test(groups = length(groupmeans), 
between.var = var(groupmeans), within.var = 65.85873, 
power=NULL, sig.level=0.05,n=n.aov)
p.rich

p.rich.df <- data.frame(n=n, power=power)
plot(n.aov,p.rich$power)




install.packages("DescTools")
library(DescTools)

set.seed(13)
nested.rich <- aov(bean.map$Richness~bean.map$Plant/factor(bean.map$Pod))

EtaSq(nested.rich, type=1, anova=TRUE)
EtaSq(fit.rich, type=1, anova=TRUE)



set.seed(13)
pwr.anova.test(k = 2, f =0.719, sig.level =0.05 , power =0.8 )

pwr.anova.test(k = 2, f =0.719, sig.level =0.05 , power =0.95 )


setwd('/Users/arifinabintarti/Documents/Bean_seed_variability_Bintarti_2020/Figures')
library(patchwork)

#rarecurve
power.rich.plot
power.pd.plot

power <- (power.rich.plot | power.pd.plot)+xlab(label = "Sample Size")+
  theme(axis.title.x = element_text(size=15, face='bold',hjust=-0.4))
power
ggsave("Fig.5.eps",
       power, device = "eps",
       width = 7, height = 4, 
       units= "in", dpi = 600)

