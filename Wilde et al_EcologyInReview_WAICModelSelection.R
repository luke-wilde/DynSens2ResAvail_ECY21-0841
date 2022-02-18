# Author: Luke Wilde
# Date compiled: 08/14/2021
# Desc: Novel code to perform intital model selection of candidate Bayesian HMs with age interaction for manuscript: "Dynamic sensitivity in early life to resource availability influences population responses to mismatches in a shorebird"

# Background: We tested the hypothesis that either biomass (abundance) or median body mass (quality) changed with chick age, using Hudsonian godwit chicks from hatch to pre-fledge (0-28 d).



#### Set up ####
#create function to load and install (missing) packages
foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

#load or install packages
foo( c("R2jags", "loo"))


#see outputs in standard (not scientific notation)
options(scipen = 999)


set.seed(123)


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
setwd("C:/Users/14064/Desktop/Ecology Submission/Data/Data/Data for Analyses/Individual Mismatch/Data")

# create destination
outfolder <- "HUGO_DSR_output_3/"
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
adapt <- 10 #number of iterations to discard
bi <- 10 #number of burn in, discarded before chains stabilize
ni <- 50 #number of iterations to run
nt <- 3 #number of thinning iterations, refers to every 3rd character saved for space???
nc <- 1 #number of chains to run, my CORE i7 has 6 cores

#
setwd("C:/Users/14064/Desktop/Ecology Submission/Data/Data/Data for Analyses/Individual Mismatch/Data/HUGO_DSR_output_3/")
#

# Example Jags code -----------------------------https://github.com/andrewcparnell/jags_examples/blob/master/R%20Code/jags_WAIC.R----------------------------------

# create models and compare using Watanabe's Akaike Information Criterion (WAIC)

with_size <- "
model 
{
	for(i in 1:nnest){
		for(j in f[i]:(l[i]-1)){
			logit(phi[i,j]) <- intercept + insecteff*zinsect[i,j] + eff_int_s*zage[i,j]*zsize[i,j] + eff_size*zsize[i,j] + hatcheff*zhatch[i,1] + eff_age*zage[i,j]
		}
	}
	# likelihood:
	for(i in 1:nnest){
		for(j in (f[i]+1):l[i]){
			mu[i,j]<-phi[i,j-1]*NH[i,j-1]
			NH[i,j]~dbern(mu[i,j])
      log_lik[i,j] <- logdensity.norm(mu[i,j], NH[i,j], tau)
		}
	}
	# set priors:
	#priors are either diffuse, or in the case of eff_age, and hatcheff are from a MARK known fate model run in the same study
	
	# parameter priors
	intercept~dnorm(0,0.01)
  eff_age~dnorm(-0.0094,0.07)
  eff_int_s~dnorm(0,25)
  eff_size~dnorm(0,25)
	insecteff~dnorm(0,25)
	hatcheff~dnorm(-0.09,0.0689)
	tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
# 	#random effect priors
#   for(j in 1:n.plot){
#     alpha[j]~dnorm(0,tau.plot) ###I(-2,2)
#   }
#   #int.plot ~ dnorm(0,0.001) #initial intercept of plot effect
#   tau.plot<-pow(sigma.plot, -2)
#   sigma.plot~dunif(0,25)
#   
#   
#   #
#    for(j in 1:n.brood){
#     alpha1[j]~dnorm(0,tau.brood) ###I(-98,98)
#   }
#   #int.brood ~ dnorm(0,0.001) 
#   tau.brood<-pow(sigma.brood, -2)
#   sigma.brood~dunif(0,25)
#   
#   #
#    for(j in 1:n.year){
#     alpha2[j]~dnorm(0,tau.year) ###I(-7,7)
#   }
#   #int.year ~ dnorm(0,0.001) 
#   tau.year<-pow(sigma.year, -2) 
# sigma.year~dunif(0,25)
	
}
"



data <- list(NH=NH, zinsect=zinsect, zsize=zsize, plot=plot,year=year,brood=brood, l=l, f=f, zage=zage,  zhatch=hatch, nnest=nrow(NH),n.year=n.year,n.plot=n.plot,n.brood=n.brood)

params <- c("eff_age", "eff_int_s", "eff_size", "insecteff","log_lik")

#run the model with constant insect effect 
# with_size <- run.jags(model = with_size,	data=data, monitor = params,inits=function(){list(intercept=rnorm(1,0,0.2),insecteff=runif(1),eff_int_s=runif(1),eff_age=runif(1),hatcheff=runif(1),alpha=runif(n.plot),eff_size=runif(1), alpha1=runif(n.brood), alpha2=runif(n.year))},thin=nt, n.chains=nc, burnin=bi, adapt=adapt, sample=ni,psrf.target=1.1)

size_run <- jags(data=data, parameters.to.save = params,model.file = textConnection(with_size),
                  n.chains = 4, # Number of different starting positions
                  n.iter = 500, # Number of iterations
                  n.burnin = 200, # Number of iterations to remove at start
                  n.thin = 2)

# Get the log likelihood
log_lik_size <- size_run$BUGSoutput$sims.list$log_lik
waic(log_lik_size)



with_abund <- "
model 
{
	for(i in 1:nnest){
		for(j in f[i]:(l[i]-1)){
			logit(phi[i,j]) <- intercept + insecteff*zinsect[i,j] + eff_int*zage[i,j]*zinsect[i,j] + eff_size*zsize[i,j] + hatcheff*zhatch[i,1] + eff_age*zage[i,j]
		}
	}
	# likelihood:
	for(i in 1:nnest){
		for(j in (f[i]+1):l[i]){
			mu[i,j]<-phi[i,j-1]*NH[i,j-1]
			NH[i,j]~dbern(mu[i,j])
      log_lik[i,j] <- logdensity.norm(mu[i,j], NH[i,j], tau)
		}
	}
	# set priors:
	#priors are either diffuse, or in the case of eff_age, and hatcheff are from a MARK known fate model run in the same study
	
	# parameter priors
	intercept~dnorm(0,0.01)
  eff_age~dnorm(-0.0094,0.07)
  eff_int~dnorm(0,25)
  eff_size~dnorm(0,25)
	insecteff~dnorm(0,25)
	hatcheff~dnorm(-0.09,0.0689)
	tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
# 	#random effect priors
#   for(j in 1:n.plot){
#     alpha[j]~dnorm(0,tau.plot) ###I(-2,2)
#   }
#   #int.plot ~ dnorm(0,0.001) #initial intercept of plot effect
#   tau.plot<-pow(sigma.plot, -2)
#   sigma.plot~dunif(0,25)
#   
#   
#   #
#    for(j in 1:n.brood){
#     alpha1[j]~dnorm(0,tau.brood) ###I(-98,98)
#   }
#   #int.brood ~ dnorm(0,0.001) 
#   tau.brood<-pow(sigma.brood, -2)
#   sigma.brood~dunif(0,25)
#   
#   #
#    for(j in 1:n.year){
#     alpha2[j]~dnorm(0,tau.year) ###I(-7,7)
#   }
#   #int.year ~ dnorm(0,0.001) 
#   tau.year<-pow(sigma.year, -2) 
# sigma.year~dunif(0,25)
	
}
"



data <- list(NH=NH, zinsect=zinsect, zsize=zsize, plot=plot,year=year,brood=brood, l=l, f=f, zage=zage,  zhatch=hatch, nnest=nrow(NH),n.year=n.year,n.plot=n.plot,n.brood=n.brood)

params <- c("eff_age", "eff_int", "eff_size", "insecteff","log_lik")

#run the model with constant insect effect 
# with_size <- run.jags(model = with_size,	data=data, monitor = params,inits=function(){list(intercept=rnorm(1,0,0.2),insecteff=runif(1),eff_int_s=runif(1),eff_age=runif(1),hatcheff=runif(1),alpha=runif(n.plot),eff_size=runif(1), alpha1=runif(n.brood), alpha2=runif(n.year))},thin=nt, n.chains=nc, burnin=bi, adapt=adapt, sample=ni,psrf.target=1.1)

abund_run <- jags(data=data, parameters.to.save = params,model.file = textConnection(with_abund),
                 n.chains = 4, # Number of different starting positions
                 n.iter = 500, # Number of iterations
                 n.burnin = 200, # Number of iterations to remove at start
                 n.thin = 2)

# Get the log likelihood
log_lik_abund <- abund_run$BUGSoutput$sims.list$log_lik
waic(log_lik_abund)


#compare the two models
waic(log_lik_abund)
waic(log_lik_size)
