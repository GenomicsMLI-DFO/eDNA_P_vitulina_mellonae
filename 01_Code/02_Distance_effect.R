#Author: Marion Chevrinais
#date: August 2024

## Required packages 
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(car)
library(ggbreak)
library(lmerTest)
library(lsmeans)
library(readxl)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(DHARMa)
library('ggeffects')

#load data
data.path <- file.path(here::here(), "00_Data/02_Distance_effect/Phoca_dist_ADNe.xlsx")
phoca.dist <- readxl::read_excel(data.path, sheet = "Sheet1")

#Summary stats
phoca.dist<-phoca.dist %>% mutate(DNA_copies = as.numeric(as.character(DNA_copies)),
                                  Numero_unique_extrait = factor(Numero_unique_extrait))

stats.extraits<-phoca.dist%>% group_by(Numero_unique_extrait) %>% dplyr::summarise(count_extrait = n(), mean_extrait = mean(DNA_copies, na.rm = TRUE), sd_extrait=sd(DNA_copies, na.rm=TRUE),  se_extrait= sd_extrait / sqrt(count_extrait))

phoca.dist <- phoca.dist %>% left_join(stats.extraits, by ='Numero_unique_extrait')

stats.station<-phoca.dist%>% group_by(Station) %>% dplyr::summarise(count_station = n(), mean_station = mean(mean_extrait, na.rm = TRUE), sd_station=sd(mean_extrait, na.rm=TRUE),  se_station= sd_station / sqrt(count_station))

phoca.graph <- phoca.dist %>% left_join(stats.station, by =('Station'))

#figure

phoca.graph$Percentage_positive <- as.factor(phoca.graph$Percentage_positive) #depending on the value in column Type the percentage is for qPCR replicates or samples
phoca.graph <- phoca.graph %>% filter(phoca.graph$Type=='sample')

Figure <- ggplot(phoca.graph, aes(x=dist_phoca_m, group=Station)) +
  geom_point(aes(y=mean_station, size = Percentage_positive))+
  geom_errorbar(aes(ymin=mean_station-se_station, ymax=mean_station+se_station), width=0
  )+
  scale_y_continuous(
    # Features of the first axis
    name = "DNA copies per qPCR reaction") + 
  theme_bw()

Figure

ggsave(file.path(here::here(), "02_Results/Figure3.tiff"),
       Figure,
       height = 5, width = 6, scale = 1)

#Modelization
#data formatting
xtabs(~ dist_phoca_m + Numero_unique_extrait + qPCR_rep, phoca.dist)
hist(phoca.dist$DNA_copies) #exploratory graph showing a high frequency of 0 DNA copy

phoca<-as_tibble(phoca.dist)
phoca$group<-paste(phoca$dist_phoca_m, phoca$Numero_unique_extrait, sep="_")
#summary stats by group DNA extract_phoca distance
test<-phoca %>% group_by(group) %>% summarise(n = n(),mean = mean(DNA_copies), sd=sd(DNA_copies))
limits<-aes(x=group, ymin=mean-sd, ymax=mean+sd) #set limits
dodge <- position_dodge(width = 0.9)
# histogram with stanadr deviation for each group
g1 <-ggplot(test) +
  geom_bar(aes(x=group, y=mean),position = dodge, stat="identity", fill="lightblue", alpha=0.9)+
  geom_errorbar( limits,position = dodge, width=0.4, colour="black", alpha=0.9, size=.5)  +
  ggtitle("histogram showing standard deviation by group ")

g1

#zero-inflated model
m5 <- glmmTMB(DNA_copies ~ factor(dist_phoca_m) + (1|Numero_unique_extrait), zi = ~ factor(dist_phoca_m), family = "nbinom2", data = phoca.dist) #model construction, pseudo replication of qPCR in Numero_unique_extrait
plot(residuals(m5) ~ fitted(m5))
summary(m5) #significance of the model
testDispersion(m5)
simulationOutput <- simulateResiduals(fittedModel = m5, plot = F)
plot(simulationOutput)

plot(ggpredict(m5, terms = "dist_phoca_m"))
predictions <- ggpredict(m5, c("dist_phoca_m"))

g2<- ggplot(
  predictions,
  aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, colour = group, fill = group)
) +
  geom_line() +
  geom_ribbon(alpha = .1, colour = NA) +
  theme_minimal()

ggsave(file.path(here::here(), "02_Results/Model_predictions.tiff"),
       g2,
       height = 5, width = 6, scale = 1)