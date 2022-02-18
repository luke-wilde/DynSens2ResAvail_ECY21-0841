# Author: Luke Wilde
# Date compiled: 08/14/2021
# Desc: Novel code to select top model from global Bayesian HM using Bayes Factors weighting for manuscript: "Dynamic sensitivity in early life to resource availability influences population responses to mismatches in a shorebird"

# Background: We identified which factors influenced survival of Hudsonian godwit chicks from hatch to pre-fledge (0-28 d). This was a known fate procedure, with code written with the help of Emily Weiser (USGS Anchorage).

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#

#load the packages 
x <- c("runjags", "rjags", "arm"); lapply(x, require, character.only = T)

#
#

# ensures that all parameters would be output. Otherwise, runjags has a ceiling on the number of parameters  summarized
runjags.options(force.summary=T)

#
#

#define functions

# function to calculate the last day alive:
get.last <- function(m) max(which(m %in% c(0,1)))
get.first <- function(m) min(which(m %in% c(0,1)))

# logit and inverse logit functions:
fn_logit <- function(X) log(X / (1 - X))
fn_invlogit <- function(X)	1/(1+exp(-(X)))

# 
# 

#set working directory
#this will need to be updated on each machine
#always reload wd after running the script

#
init <- Sys.time()	# record the start time for the example script

#
setwd("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/Individual Mismatch/Data")

#
outfolder <- "HUGO_DSR_output/"
dir.create(outfolder)


#now we define the parameters
#load the data
#for those that do not change across days (constant), we unlist to get them into vector form, then convert them to numerics using the coercive 'as.numeric()'
#for those that change over time (dynamic), we instead keep as a matrix

#
NH <- as.matrix(read.csv("HUGO_chick_hist.v3try.csv", fileEncoding="UTF-8-BOM", header = F));head(NH) # survival history, 1 row for each individual, 1 for all days survival, 0 for dead. 2 trailing NAs before 0 when not confirmed mortality
#
year <- as.matrix(unlist(read.csv("HUGO_year.csv", fileEncoding="UTF-8-BOM", header = F, ))) #study year
year <- as.numeric(as.factor(year))
n.year <- length((unique(year)))
#
insect <- as.matrix(read.csv("HUGO_insect_v1.csv", fileEncoding="UTF-8-BOM", header = F)) #invertebrate biomass (see manuscript for details)
insect#

#
size <- as.matrix(read.csv("HUGO_size.csv", fileEncoding="UTF-8-BOM", header = F)) # invertebrate body mass (see manuscript for details)
size#

age <- as.matrix(read.csv("HUGO_chick_age.full.csv", fileEncoding="UTF-8-BOM", header = F)) #chick age

#
zage <- as.matrix(read.csv("HUGO_chick_zage.full.csv", fileEncoding="UTF-8-BOM", header = F)) # standardized age
#

brood <- as.matrix(unlist(read.csv("HUGO_broods.csv", fileEncoding="UTF-8-BOM", header = F))) # chick brood ID
brood <- as.numeric(as.factor(brood))
n.brood <- length((unique(brood)))
#
hatch <- as.matrix(read.csv("HUGO_hatch_dates_rescaled.csv", fileEncoding="UTF-8-BOM", header = F)) # chick hatch date
#
plot <- as.matrix(unlist(read.csv("HUGO_plot.csv", fileEncoding="UTF-8-BOM", header = F))) # study plot: North or South plot (see manuscript for details)
plot <- as.numeric(as.factor(plot))
n.plot <- length((unique(plot)))


#this shouldnt really be f. f should come from a function to calculate the day of first encounter, but for me that is 1 for all
# first day the chick was monitored; any nest history before that is ignored 
# in the analysis):

#
f <- as.numeric(apply(NH, 1, get.first))	
l <- as.numeric(apply(NH, 1, get.last))	# last day that each nest was alive

# rescale all variables by 2 SD
#Note: the 'age' variable is created below and hatch is brought in from a standardized .csv file
#
zinsect <- rescale(insect)
zsize <- rescale(size)


#
## 2. Define parameters from metadata

#these are the parameters from my MARK model
m <- 0.96	# mean DSR (at midpoint of pre-fledge stage)
B <- 0.077 	# effect of chick age on DSR on the logit scale, from Mark
d <- 21	# duration of the nesting stage in days. Used an odd number so 
# that centered age had a mean of exactly 0. 
ad <- 1	# age at which nests were discovered. Age must never be 0 for 
# the purposes of the data simulation and analysis (indexing), so ad = 1 
# corresponds to age at discovery = 0 as described in the accompanying paper.
adcv <- 0.33	# coefficient of variation (CV) of the age at discovery # this is from our captures of chicks both on the nests and found off nest

# Begin nest survival model, first assume constant effect of insects, then assume insect*Age

#
#
# for each day of the pre-fledge, calculate the centered value for age. Tidier than the absolute value
#age is 0 and thus the mean DSR is m = 0.96 applies to the midpoint of the pre-fledge
#mean age is the midpoint
age_mean <- round(d/2+0.01,0)
#
#zage <- c(1:d) - age_mean
#
#
# calculate the age when chicks were active. Center the age on the day of the chick histories to the observed mean (will not be the midpoint of the chick stage before some fail early)
ages <- numeric(0)
for(i in 1:nrow(NH)) ages <- c(ages, which(NH[i,]==1))
age_c <- as.vector(1:ncol(NH) - mean(ages))
n.age <- length(unique(age_c))

#
#
#

#first, define input to JAGS:
adapt <- 600 #number of iterations to discard
bi <- 1000 #number of burn in, discarded before chains stabilize
ni <- 10000 #number of iterations to run
nt <- 3 #number of thinning iterations, refers to every 3rd character saved for space???
nc <- 3 #number of chains to run, my CORE i7 has 6 cores

#
setwd("C:/Users/14064/Dropbox/Chapter 2/Data/Data for Analyses/Individual Mismatch/Data/Outputs")
#

# run global model with BF selection
sink("HUGO_match.model_global.txt")
cat("model {
	for(i in 1:nnest){
		for(j in f[i]:(l[i]-1)){
			logit(phi[i,j]) <- intercept + w.c*insecteff*zinsect[i,j] + w.b*eff_size*zsize[i,j] + w.e*eff_int_s*zage[i,j]*zsize[i,j] + w.d*hatcheff*zhatch[i,1] + w.a*eff_age*zage[i,j]+ alpha[plot[j]] + alpha2[year[j]] + alpha1[brood[j]]
		}
	}
	# likelihood:
	for(i in 1:nnest){
		for(j in (f[i]+1):l[i]){
			mu[i,j]<-phi[i,j-1]*NH[i,j-1]
			NH[i,j]~dbern(mu[i,j])
		}
	}
	# set priors:
	#model priors
  tau.total ~ dgamma(3.29,7.8)
  K <- w.a + w.b + w.c + w.d + w.e
  tau.model <- K*tau.total

  w.a ~ dbern(0.5)
  w.b ~ dbern(0.5)
  w.c ~ dbern(0.5)
  w.d ~ dbern(0.5)
	w.e ~ dbern(0.5)
	
	#priors are either diffuse, or in the case of eff_age, and hatcheff are from a MARK known fate model run in the same study
	
	# parameter priors
	intercept~dnorm(0,0.01)
  eff_age~dnorm(0,25)
  eff_size~dnorm(0,25)
  eff_int_s~dnorm(0,25)
	insecteff~dnorm(0,25)
	hatcheff~dnorm(0,25)
	
	#random effect priors
  for(j in 1:n.plot){
    alpha[j]~dnorm(0,tau.plot) ###I(-2,2)
  }
  #int.plot ~ dnorm(0,0.001) #initial intercept of plot effect
  tau.plot<-pow(sigma.plot, -2)
  sigma.plot~dunif(0,25)
  
  
  #
   for(j in 1:n.brood){
    alpha1[j]~dnorm(0,tau.brood) ###I(-98,98)
  }
  #int.brood ~ dnorm(0,0.001) 
  tau.brood<-pow(sigma.brood, -2)
  sigma.brood~dunif(0,25)
  
  #
   for(j in 1:n.year){
    alpha2[j]~dnorm(0,tau.year) ###I(-7,7)
  }
  #int.year ~ dnorm(0,0.001) 
  tau.year<-pow(sigma.year, -2) 
sigma.year~dunif(0,25)
  
  
	# back-transform expected DSR and nest survival to fledging
	
	dsr_mean <- 1/(1+exp(-(intercept)))
	nsurv_mean <- dsr_mean^d
	
		# or, use age-specific estimates of DSR to fully account for the age effect:
	# 
#	for(j in 1:length(age)){
	#	dsr_age[j] <- 1/(1+exp(-(intercept + eff_age*age[j])))
#	}
#	nsurv_byage <- prod(dsr_age)
	
}",fill=TRUE)
sink()


#run the model with constant insect effect 
global_BF <- run.jags(model = "HUGO_match.model_global.txt",	data=list(NH=NH, zinsect=zinsect, zsize=zsize, plot=plot,year=year,brood=brood, l=l, f=f, d=d, zage=zage,  zhatch=hatch, nnest=nrow(NH),n.year=n.year,n.plot=n.plot,n.brood=n.brood), monitor=c("dsr_mean", "nsurv_mean", "eff_int_s", "eff_size", "insecteff" ,"eff_age", "hatcheff","intercept", "w.a", "w.b", "w.c", "w.d", "w.e"),inits=function(){list(intercept=rnorm(1,0,0.2),insecteff=runif(1), eff_int_s=runif(1),eff_age=runif(1),hatcheff=runif(1),alpha=runif(n.plot),eff_size=runif(1), alpha1=runif(n.brood), alpha2=runif(n.year))},thin=nt, n.chains=nc, burnin=bi, adapt=adapt, sample=ni,psrf.target=1.1, method='rjparallel')

#dic cannot be assessed by parrellel chains
summary(global_BF)

BF_biomass <- (global_BF$summary$statistics[9,1]/((1-global_BF$summary$statistics[9,1])))/(0.5/((1-global_BF$summary$statistics[9,1])))
BF_size <- (global_BF$summary$statistics[10,1]/((1-global_BF$summary$statistics[10,1])))/(0.5/((1-global_BF$summary$statistics[10,1])))
BF_AS <- (global_BF$summary$statistics[11,1]/((1-global_BF$summary$statistics[11,1])))/(0.5/((1-global_BF$summary$statistics[11,1])))
BF_H <- (global_BF$summary$statistics[12,1]/((1-global_BF$summary$statistics[12,1])))/(0.5/((1-global_BF$summary$statistics[12,1])))


BF_biomass; BF_size; BF_AS; BF_H



BF_A/(1+BF_A); BF_AS/(1+BF_AS); BF_S/(1+BF_S); BF_R/(1+BF_R); BF_H/(1+BF_H)






#top model: probability of support
sink("HUGO_match.model_topmodel.txt")
cat("model {
	for(i in 1:nnest){
		for(j in f[i]:(l[i]-1)){
			logit(phi[i,j]) <- intercept + w.a*insecteff*zinsect[i,j] + w.b*eff_int_s*zage[i,j]*zsize[i,j] + w.c*hatcheff*zhatch[i,1] + w.d*eff_size*zsize[i,j] + alpha[plot[j]] + alpha2[year[j]] + alpha1[brood[j]]
		}
	}
	# likelihood:
	for(i in 1:nnest){
		for(j in (f[i]+1):l[i]){
			mu[i,j]<-phi[i,j-1]*NH[i,j-1]
			NH[i,j]~dbern(mu[i,j])
		}
	}
	# set priors:
	#model priors
  tau.total ~ dgamma(3.29,7.8)
  K <- w.a + w.b + w.c + w.d
  tau.model <- K*tau.total

  w.a ~ dbern(0.5)
  w.b ~ dbern(0.5)
  w.c ~ dbern(0.5)
  w.d ~ dbern(0.5)
	
	#priors are either diffuse, or in the case of eff_age, and hatcheff are from a MARK known fate model run in the same study
	
	# parameter priors
	intercept~dnorm(0,0.01)
  eff_size~dnorm(0, 25)
  eff_int_s~dnorm(0, 25)
  insecteff~dnorm(0, 25)
  eff_age~dnorm(0, 25)
	hatcheff~dnorm(0, 25)
	
	#random effect priors
  for(j in 1:n.plot){
    alpha[j]~dnorm(0,tau.plot) ###I(-2,2)
  }
  #int.plot ~ dnorm(0,0.001) #initial intercept of plot effect
  tau.plot<-pow(sigma.plot, -2)
  sigma.plot~dunif(0,25)
  
  
  #
   for(j in 1:n.brood){
    alpha1[j]~dnorm(0,tau.brood) ###I(-98,98)
  }
  #int.brood ~ dnorm(0,0.001) 
  tau.brood<-pow(sigma.brood, -2)
  sigma.brood~dunif(0,25)
  
  #
   for(j in 1:n.year){
    alpha2[j]~dnorm(0,tau.year) ###I(-7,7)
  }
  #int.year ~ dnorm(0,0.001) 
  tau.year<-pow(sigma.year, -2) 
sigma.year~dunif(0,25)
  
  
	# back-transform expected DSR and nest survival to fledging
	
	dsr_mean <- 1/(1+exp(-(intercept)))
	nsurv_mean <- dsr_mean^d
	
		# or, use age-specific estimates of DSR to fully account for the age effect:
	# 
#	for(j in 1:length(age)){
	#	dsr_age[j] <- 1/(1+exp(-(intercept + eff_age*age[j])))
#	}
#	nsurv_byage <- prod(dsr_age)
	
}",fill=TRUE)
sink()


#run the model with constant insect effect 
topmodel <- run.jags(model = "HUGO_match.model_topmodel.txt",	data=list(NH=NH, zinsect=zinsect, zsize=zsize, plot=plot,year=year,brood=brood, l=l, f=f, d=d, zhatch=hatch,zage=zage,nnest=nrow(NH),n.year=n.year,n.plot=n.plot,n.brood=n.brood), monitor=c("dsr_mean", "nsurv_mean", "hatcheff", "eff_int_s","eff_size","insecteff" ,"intercept", "w.a", "w.b", "w.c", "w.d"),inits=function(){list(intercept=rnorm(1,0,0.2),insecteff=runif(1), eff_int_s=runif(1),hatcheff=runif(1),alpha=runif(n.plot),eff_size=runif(1), alpha1=runif(n.brood), alpha2=runif(n.year))},thin=nt, n.chains=nc, burnin=bi, adapt=adapt, sample=ni,psrf.target=1.1, method='rjparallel')


summary(topmodel)


BF_R <- (topmodel$summary$statistics[8,1]/((1-topmodel$summary$statistics[8,1])))/(0.5/((1-topmodel$summary$statistics[8,1])))
BF_AS <- (topmodel$summary$statistics[9,1]/((1-topmodel$summary$statistics[9,1])))/(0.5/((1-topmodel$summary$statistics[9,1])))
BF_H <- (topmodel$summary$statistics[10,1]/((1-topmodel$summary$statistics[10,1])))/(0.5/((1-topmodel$summary$statistics[10,1])))
BF_S <- (topmodel$summary$statistics[11,1]/((1-topmodel$summary$statistics[11,1])))/(0.5/((1-topmodel$summary$statistics[11,1])))

#bayes factors
BF_R; BF_AS; BF_H; BF_S
exp(BF_R); exp(BF_AS); exp(BF_H); exp(BF_S)

#probability of inclusion (aka support)
BF_R/(1+BF_R); BF_AS/(1+BF_AS); BF_H/(1+BF_H); BF_S/(1+BF_S)



options(max.print=1000000); options(scipen = 999)



#top model - calculate the beta coeffecients - best if run with only top covariates, but BF not in the model to save effective degrees of freedom
sink("HUGO_match.model_topmodel_noBF.txt")
cat("model {
	for(i in 1:nnest){
		for(j in f[i]:(l[i]-1)){
			logit(phi[i,j]) <- intercept + w.a*insecteff*zinsect[i,j] + w.b*hatcheff*zhatch[i,1] + w.c*eff_int_s*zage[i,j]*zsize[i,j] + alpha[plot[j]] + alpha2[year[j]] + alpha1[brood[j]]
		}
	}
	# likelihood:
	for(i in 1:nnest){
		for(j in (f[i]+1):l[i]){
			mu[i,j]<-phi[i,j-1]*NH[i,j-1]
			NH[i,j]~dbern(mu[i,j])
		}
	}
	# set priors:
	#priors are either diffuse, or in the case of eff_age, and hatcheff are from a MARK known fate model run in the same study
	tau.total ~ dgamma(3.29,7.8)
  K <- w.a + w.b + w.c
  tau.model <- K*tau.total

  w.a ~ dbern(0.5)
  w.b ~ dbern(0.5)
  w.c ~ dbern(0.5)

	
	# parameter priors
	intercept~dnorm(0,0.01)
  eff_int_s~dnorm(0, 25)
  insecteff~dnorm(0, 25)
	hatcheff~dnorm(0, 25)
	
	#random effect priors
  for(j in 1:n.plot){
    alpha[j]~dnorm(0,tau.plot) ###I(-2,2)
  }
  #int.plot ~ dnorm(0,0.001) #initial intercept of plot effect
  tau.plot<-pow(sigma.plot, -2)
  sigma.plot~dunif(0,25)
  
  
  #
   for(j in 1:n.brood){
    alpha1[j]~dnorm(0,tau.brood) ###I(-98,98)
  }
  #int.brood ~ dnorm(0,0.001) 
  tau.brood<-pow(sigma.brood, -2)
  sigma.brood~dunif(0,25)
  
  #
   for(j in 1:n.year){
    alpha2[j]~dnorm(0,tau.year) ###I(-7,7)
  }
  #int.year ~ dnorm(0,0.001) 
  tau.year<-pow(sigma.year, -2) 
sigma.year~dunif(0,25)
  
  
	# back-transform expected DSR and nest survival to fledging
	
	dsr_mean <- 1/(1+exp(-(intercept)))
	nsurv_mean <- dsr_mean^d
	
		# or, use age-specific estimates of DSR to fully account for the age effect:
	# 
#	for(j in 1:length(age)){
	#	dsr_age[j] <- 1/(1+exp(-(intercept + eff_age*age[j])))
#	}
#	nsurv_byage <- prod(dsr_age)
	
}",fill=TRUE)
sink()


#run the model with constant insect effect 
topmodel_noBF <- run.jags(model = "HUGO_match.model_topmodel_noBF.txt",	data=list(NH=NH, zinsect=zinsect, zsize=zsize, plot=plot,year=year,zhatch=hatch,brood=brood, l=l, f=f, d=d, zage=zage,nnest=nrow(NH),n.year=n.year,n.plot=n.plot,n.brood=n.brood), monitor=c("dsr_mean", "nsurv_mean", "eff_int_s", "insecteff","hatcheff" ,"intercept","w.a", "w.b", "w.c"),inits=function(){list(intercept=rnorm(1,0,0.2),insecteff=runif(1), eff_int_s=runif(1),alpha=runif(n.plot), hatcheff=runif(1), alpha1=runif(n.brood), alpha2=runif(n.year))},thin=nt, n.chains=nc, burnin=bi, adapt=adapt, sample=ni,psrf.target=1.1, method='rjparallel')


summary(topmodel_noBF)


BF_R <- (topmodel_noBF$summary$statistics[7,1]/((1-topmodel_noBF$summary$statistics[7,1])))/(0.5/((1-topmodel_noBF$summary$statistics[7,1])))
BF_H <- (topmodel_noBF$summary$statistics[8,1]/((1-topmodel_noBF$summary$statistics[8,1])))/(0.5/((1-topmodel_noBF$summary$statistics[8,1])))
BF_AS <- (topmodel_noBF$summary$statistics[9,1]/((1-topmodel_noBF$summary$statistics[9,1])))/(0.5/((1-topmodel_noBF$summary$statistics[9,1])))

#bayes factors
BF_R; BF_H; BF_AS
exp(BF_R);exp(BF_H); exp(BF_AS)

#probability of inclusion (aka support)
BF_R/(1+BF_R); BF_H/(1+BF_H); BF_AS/(1+BF_AS)
