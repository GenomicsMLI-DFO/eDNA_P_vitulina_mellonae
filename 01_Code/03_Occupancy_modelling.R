#Author: Marion Chevrinais
#Date: July 2024

## Required packages 
library(msocc)
library(tidyverse)
library(gridExtra)
library(rlang)
library(ggplot2)
library(tidyr)
library(gtable)

#Load data
resp.path<-file.path(here::here(), "00_Data/03_Occupancy_modelling/resp_msocc2.csv")
resp <-read.csv(resp.path, header = T,sep = ';')
head(resp)

site.path<-file.path(here::here(), "00_Data/03_Occupancy_modelling/site_msocc2.csv")
site <-read.csv(site.path, header = T,sep = ';')
head(site)

sample.path<-file.path(here::here(), "00_Data/03_Occupancy_modelling/sample_msocc3.csv")
sample <-read.csv(sample.path, header = T,sep = ';')
head(sample)

rep.path<-file.path(here::here(), "00_Data/03_Occupancy_modelling/rep_msocc2.csv")
rep <-read.csv(rep.path, header = T,sep = ';')
head(rep)

#Fit a model to the data with only intercept modeled which means that we have probabilities of presence, occurrence and detection
mod <- msocc_mod(wide_data=resp,
                 site=list(model=~1,cov_tbl=site),
                 sample=list(model=~1,cov_tbl=sample),
                 rep=list(model=~1,cov_tbl=rep),
                 progress=F)

#numerical summaries of the model, values of psi, theta and p - unique combinations
posterior_summary(mod,print=T)

#For a more in-depth look at each of the site, sample, and replicate levels and a description of uncertainty, we can specify the level in posterior_summary.
posterior_summary(mod,level='site',print=T)
posterior_summary(mod,level='sample',print=T)
posterior_summary(mod,level='rep',print=T)

#constant psi, theta as a function of covariates, constant p
# make sure that covariates are numerical
site<-site%>% mutate(volume=as.numeric(volume))
sample<-sample%>% mutate(month=as.numeric(month), species=as.numeric(species))

#Volume
mod_cov1 <- msocc_mod(wide_data=resp,
                      site=list(model=~1,cov_tbl=site),
                      sample=list(model=~volume,cov_tbl=sample),
                      rep=list(model=~1,cov_tbl=rep),
                      progress=F)

#numerical summaries of the model, values of psi, theta and p - unique combinations
posterior_summary(mod_cov1,print=T)
psum<-posterior_summary(mod_cov1,level='sample',print=T)
psum
str(psum)

#add a site_sample column to psum and to sample
psum_u<-unite(psum, col='site_sample', c('site', 'sample'), sep='_')
sample_u<-unite(sample, col = 'site_sample', c('site','sample'),sep='_')

#change name of intervals columns
names(psum_u)[5]<-paste("low")
names(psum_u)[6]<-paste("high")
data.graph<- sample_u %>% left_join(psum_u %>% dplyr::select(site_sample, mean, median, low, high))
data.graph<-na.omit(data.graph)

ggplot(data.graph,aes(x=mean, y=site_sample, shape= factor(volume)
))+
  geom_errorbar(aes(xmin=low, xmax=high), width=0.3)+
  geom_point(aes(x=mean), cex=2)+
  ylab ('sample') +
  xlab('psi')  + 
  ggtitle ('95% credibility intervals for theta by site with volume as covariable')+
  theme_bw()


#Month
mod_cov2 <- msocc_mod(wide_data=resp,
                      site=list(model=~1,cov_tbl=site),
                      sample=list(model=~month,cov_tbl=sample),
                      rep=list(model=~1,cov_tbl=rep),
                      progress=F)

#numerical summaries of the model, values of psi, theta and p - unique combinations
posterior_summary(mod_cov2,print=T)
psum<-posterior_summary(mod_cov2,level='sample',print=T)
psum
str(psum)

#add a site_sample column to psum and to sample
psum_u<-unite(psum, col='site_sample', c('site', 'sample'), sep='_')
sample_u<-unite(sample, col = 'site_sample', c('site','sample'),sep='_')

#change name of intervals columns
names(psum_u)[5]<-paste("low")
names(psum_u)[6]<-paste("high")

data.graph<- sample_u %>% left_join(psum_u %>% dplyr::select(site_sample, mean, median, low, high))
data.graph<-na.omit(data.graph)

ggplot(data.graph,aes(x=mean, y=site_sample, shape=factor(month)))+
  geom_errorbar(aes(xmin=low, xmax=high))+
  geom_point(aes(x=mean), cex=2)+
  ylab ('site_sample') +
  xlab('theta')  + 
  ggtitle ('95% credibility intervals for theta by site_sample with month as covariable')+
  theme_bw()

#species presence
mod_cov3 <- msocc_mod(wide_data=resp,
                      site=list(model=~1,cov_tbl=site),
                      sample=list(model=~species,cov_tbl=sample),
                      rep=list(model=~1,cov_tbl=rep),
                      progress=F)

#numerical summaries of the model, values of psi, theta and p - unique combinations
posterior_summary(mod_cov3,print=T)
psum<-posterior_summary(mod_cov3,level='sample',print=T)
psum
str(psum)

#add a site_sample column to psum and to sample
psum_u<-unite(psum, col='site_sample', c('site', 'sample'), sep='_')
sample_u<-unite(sample, col = 'site_sample', c('site','sample'),sep='_')

#change name of intervals columns
names(psum_u)[5]<-paste("low")
names(psum_u)[6]<-paste("high")

data.graph<- sample_u %>% left_join(psum_u %>% dplyr::select(site_sample, mean, median, low, high))
data.graph<-na.omit(data.graph)

ggplot(data.graph,aes(x=median,y=site_sample, shape=factor(species)))+
  geom_errorbar(aes(xmin=low, xmax=high))+
  geom_point(aes(x=mean), cex=2)+
  ylab ('site_sample') +
  xlab('theta')  + 
  ggtitle ('95% credibility intervals for theta by site-sample with species presence as covariable')+
  theme_bw()

#3 levels volume and species presence and month
mod_cov4 <- msocc_mod(wide_data=resp,
                      site=list(model=~1,cov_tbl=site),
                      sample=list(model=~species+month+volume,cov_tbl=sample),
                      rep=list(model=~1,cov_tbl=rep),
                      progress=F)

#numerical summaries of the model, values of psi, theta and p - unique combinations
#psum_site<-posterior_summary(mod_cov4,level='site',print=T)
#psum_site
psum_sample<-posterior_summary(mod_cov4,level='sample',print=T)
sta_mod4_cov<-posterior_summary(mod_cov4,print=T)
psum_sample
str(psum_sample)

#add a site_sample column to psum and to sample
psum_sample<-unite(psum_sample, col='site_sample', c('site', 'sample'), sep='_')
sample_u<-unite(sample, col = 'site_sample', c('site','sample'),sep='_')

#change name of intervals columns
#names(psum_site)[4]<-paste("low")
#names(psum_site)[5]<-paste("high")

names(psum_sample)[5]<-paste("low")
names(psum_sample)[6]<-paste("high")

data.graph<- sample_u %>% left_join(psum_sample %>% dplyr::select(site_sample, mean, median, low, high))
data.graph<-na.omit(data.graph) #%>% transform(data.graph, volume = as.numeric(volume))

data.graph$volume <- as.factor(data.graph$volume)

g1 <- ggplot(data.graph,aes(x=median,y=site_sample, shape=factor(species), fill=factor(month)
))+
  geom_errorbar(aes(xmin=low, xmax=high#, color=factor(volume)
  ))+
  scale_shape_manual(values = c(21, 22) )+
  # scale_color_manual(values = c("black","white", 'grey50','red') )+
  scale_fill_manual(values=c('grey80','grey50','black'))+
  geom_point(aes(x=mean, size=data.graph$volume))+
  scale_size_manual(values = c(1,3,5,7))+
  ylab ('site_sample') +
  #geom_line(aes(linetype=factor(month), color=factor(volume)))+
  xlab('theta')  + 
  ggtitle ('95% credibility intervals for theta by site-sample with species presence, month and volume as covariable')+
  theme_bw()

ggsave(file.path(here::here(), "02_Results/Figure4.tiff"),
       g1,
       height = 6, width = 7, scale = 1)

#comparison of models
waic(mod_cov1,type=2)
waic(mod_cov2,type=2)
waic(mod_cov3,type=2)
waic(mod_cov4,type=2)
