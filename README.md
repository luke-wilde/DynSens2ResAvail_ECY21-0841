# DynSens2ResAvail_ECY21-0841.R1


##

Overview
Novel code developed for manuscript in review: "Dynamic sensitivity in early life to resource availability influences population responses to mismatches in a shorebird". We studied the effect of daily invertebrate biomass and daily median invertebrate body mass on Hudsonian godwit (Limosa heamastica) growth and survival. Presented code performed three steps of the analysis: (1) code to calculate the Watanabe's Akaike Information Criterion for selection of the top model from candidate Bayesian heirarchical models, (2) the creation of the dataset and parametirization of a global model using Just Another Gibs Sampler (JAGS), (3) use of Link and Barker's (2006) approach to calculate Bayes Factors and weights for probability of including covariates in top model, (4) population level overlap models to compare among peak and whole demand models (see manuscript for details).

##

Methods
All scrips were written in Program R (v.4.1.0). Dependencies are loaded within scripts.

##

Files List:

Code file 1 name: Wilde et al._EcologyInReview_WAICModelSelection_.R
Code file 1 description:

Code file 2 name: Wilde et al._EcologyInReview_JAGSmodel_.R
Code file 2 description:

Code file 3 name: Wilde et al._EcologyInReview_BayesFactorsTopModelSelection_.R
Code file 3 description: 

Code file 4 name: Wilde et al._EcologyInReview_OverlapMismatchModels_.R
Code file 4 description: 

##

Code specfic information: Wilde et al._EcologyInReview_WAIC_.R

Code specific information: Wilde et al._EcologyInReview_JAGSmodel_.R

Code specific information: Wilde et al._EcologyInReview_BayesFactorsTopModelSelection_.R

Code specific information: Wilde et al._EcologyInReview_OverlapModel_.R





