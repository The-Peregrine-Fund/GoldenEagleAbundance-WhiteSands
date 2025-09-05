############################
# Supplemental materials for:
# B. W. Rolek, D. J. Harrison, D. W. Linden,  C. S. Loftin, 
# P. B. Wood. 2021. Associations among breeding 
# conifer-associated birds, forestry treatments, 
# years-since-harvest, and vegetation characteristics in 
# regenerating stands. 
#############################
## ---- basic Poisson --------
#############################################
# This file is not needed to replicate analyses.
# We provide the simplest version of the model used
# as a template for other studies
#############################################
# software used
# JAGS 4.3.0 
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
library (jagsUI) # v1.5.1
library (abind)
load ("./data/example-data.Rdata")
# data manipulation
# TOM YOU CAN PROBABLY IGNORE THIS PART
# BUT LOOK AT THE FINAL DATA LIST 
datalfoc$SPP <- length(spp.list.foc)
yr <- array(NA, dim=c(dim (ab)[1], 9) )
yr[,1:3] <- 1; yr[,4:6] <- 2; yr[,7:9] <- 3
datalfoc$yr <- yr
s.year <- array(NA, dim=c(114, 9))
s.year[,1:3] <- 1; s.year[,4:6] <- 2; s.year[,7:9] <- 3
datalfoc$s.year <- s.year

# Add visit and year
# for removal model
datalfoc$visit <- ifelse(datalfoc$yr_rot %in% c(1,4,7), 1, 
                         ifelse(datalfoc$yr_rot %in% c(2,5,8), 2,
                                3))
datalfoc$year <- ifelse(datalfoc$yr_rot %in% c(1,2,3), 1, 
                        ifelse(datalfoc$yr_rot %in% c(4,5,6), 2,
                               3))

nobs <- datalfoc$nobs
int <- datalfoc$int
site <- datalfoc$site
yr_rot <- datalfoc$yr_rot
year <- datalfoc$year
visit <- datalfoc$visit

datalfoc$nvisit <- 3
datalfoc$YR <- 3
# print sample sizes 
apply(ab2[,1:2,,,dimnames(ab2)[[5]] %in% spp.list.foc], c(5), sum, na.rm=T)

# Subset data to species=Yellow-bellied flycatcher spp.list.foc[2]
spp <- spp.list.foc[1]
spp.num<- which(dimnames(nobs)[[3]]==spp)
Nav <- apply(ab2[,1:2,,,spp], c(1,4), sum, na.rm=T)

no <- abind(Nav[, c(1,2,3)],
            Nav[, c(4,5,6)],
            Nav[, c(7,8,9)],
            along=3 )
dimnames(no) <- list(site= dimnames(Nav)[[1]],                       
                     visitnum= 1:3,
                     year= 2013:2015)

#***********************
#* Run the above code but
#* START HERE TOM
#* These are all the data needed
#* to run the model. You can add covariates. 
#* See comments about dimensions for each data
datalist <- list()
datalist$nobs <- no # an array of 3 dims: nsites x nvisits x nyears
# the following are all nL in length, ie the number of detections
datalist$int <- datalfoc$int[datalfoc$species==spp.num] # the time interval of first detection for each eagle observed during a count
datalist$site <- datalfoc$site[datalfoc$species==spp.num] # the site of detection for each individual detection
datalist$visit <- datalfoc$visit[datalfoc$species==spp.num] # the visit of detection
datalist$year <- datalfoc$year[datalfoc$species==spp.num] # the year detected
# dimensions
datalist$nL <- length(datalfoc$species[datalfoc$species==spp.num]) # the number of detections
datalist$nsites <- datalfoc$nsites
datalist$nyr <- 3
datalist$nvisits <- 3
datalist$nR <- datalfoc$R
# Examine full dataset
str(datalist)
# These two should be equal to the number of detections
# and each other. Not quite the case here.
# But do better than I did! 
length(datalist$int)
sum(datalist$nobs, na.rm=T)
####################################
# (1) basic Poisson model
####################################
cat("
    model {
    # p.pa.beta: prob. of availability
    # p.pem.beta: prob. of temporary emigration
    # lam.beta: abundance (log-link scale, use exp function to backtransform)
    ##### PRIORS ###############################################
    pa.beta <- logit(p.pa.beta)
    p.pa.beta ~ dbeta(1, 1)
    pem.beta <- logit(p.pem.beta) # backtransform to probability
    p.pem.beta ~ dbeta(1, 1) # probability of emigration
    lam.beta ~ dnorm(0, 0.01)
    
    ##### REMOVAL #####################################
    # Removal submodel
    for (l in 1:nL) {
    int[l] ~ dcat(pi.pa.c[site[l], visit[l], year[l], ]) # removal class frequencies
    } # L

    for (i in 1:nsites) {
    for (t in 1:nyr) {
    # Abundance (inserted here to reduce number of loops)
    N[i,t] ~ dpois(lambda[i,t])
    log(lambda[i,t]) <- lam.beta # add covariates for abundance here
    
    # More removal submodel
    for (j in 1:nvisits) {
    for (r in 1:nR){
    pi.pa[i,j,t,r] <- p.a[i,j,t]*pow(1-p.a[i,j,t], (r-1))
    pi.pa.c[i,j,t,r] <- pi.pa[i,j,t,r] / pcap[i,j,t]
    }  #nR
    pcap[i,j,t] <- sum(pi.pa[i,j,t,1:nR])
    
    # Detection  and emigration models 
    pmarg[i,j,t] <-  pcap[i,j,t] * pem[i,j,t]
    logit(p.a[i,j,t]) <- pa.beta # add covariates for availability (time-removal) here
    logit(pem[i,j,t]) <- pem.beta # add covariates for temporary emigration
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,j,t] ~ dbin(pmarg[i,j,t], N[i,t])  

    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,j,t] ~ dbin(pmarg[i,j,t], N[i,t]) # create new realization of model
    e.p[i,j,t] <- pmarg[i,j,t] * N[i,t] # original model prediction
    E.p[i,j,t] <- pow((nobs[i,j,t]- e.p[i,j,t]),2)/(e.p[i,j,t]+0.5)
    E.New.p[i,j,t]<- pow((nobs.fit[i,j,t]-e.p[i,j,t]),2)/(e.p[i,j,t]+0.5)
    }}} #nyr #nsites 
    fit.p <- sum(E.p[1:nsites, 1:nvisits, 1:nyr])
    fit.new.p <- sum(E.New.p[1:nsites, 1:nvisits, 1:nyr])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:nyr){
    Ntot[t] <- sum(N[1:nsites,t])
    # Can calculate density if you have a definable area of survey
    # Example for circular point counts
    # D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } # nyr
    } # End model
    ",file="./model-basic-Poisson.txt")
  
  inits <- function(){  list(
    N = apply(no, c(1,3), max, na.rm=T),
    lam.beta = mean(Nav, na.rm=T),
    p.pa.beta= runif(1, 0.01, 0.99),
    p.em.beta= runif(1, 0.01, 0.99)
  )  }
  
  params <- c("pa.beta", "p.pa.beta",
              "pem.beta", "p.pem.beta",
              "lam.beta",  
              "Ntot",  "bayesp"
  )
  
  # MCMC settings
  ni <- 10000  ;   nb <- 5000   ;   nt <- 5   ;   nc <- 3 ; na=100
  #ni <- 200  ;   nb <- 100   ;   nt <- 1   ;   nc <- 3 ; na=100
  # Run JAGS
  out <- jags(datalist, inits=inits, params, 
              "./model-basic-Poisson.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm")
  )
  
  fn<- paste( "./", spp, "_basic-pois.RData", sep="" )
  save(list= c("out"), file=fn)

#*********************************
#* Below is a zero-inflated POisson
#* model in case you have poor model fit. 
#* I don't think we'll need this because it
#* deals with excess zeroes and most of your 
#* sites were occupied. 
#* Requires modification to work 
#********************************* 
# ## ---- basic ZIP --------
# # software used
# # JAGS 4.3.0 
# # R version 4.0.2 (2020-06-22) -- "Taking Off Again"
# library (jagsUI) # v1.5.1
# load ("./DATA.Rdata")
# # set up data
# datalfoc$SPP <- length(spp.list.foc)
# yr <- array(NA, dim=c(dim (ab)[1], 9) )
# yr[,1:3] <- 1; yr[,4:6] <- 2; yr[,7:9] <- 3
# datalfoc$yr <- yr
# s.year <- array(NA, dim=c(114, 9))
# s.year[,1:3] <- 1; s.year[,4:6] <- 2; s.year[,7:9] <- 3
# datalfoc$s.year <- s.year
# 
# datalfoc$ba <- datalfoc$CovsLam[, "ba"]
# nobs <- datalfoc$nobs
# dclass <- datalfoc$dclass
# int <- datalfoc$int
# site <- datalfoc$site
# yr_rot <- datalfoc$yr_rot
# 
# # print sample sizes 
# apply(ab2[,1:2,,,dimnames(ab2)[[5]] %in% spp.list.foc], c(5), sum, na.rm=T)
# 
# ###################################
# # (2) basic zero-inflated Poisson model
# ###################################
# cat("
#     model {
#     ##### PRIORS ###############################################
#     pa.beta <- logit(p.pa.beta) # backtransform to probability
#     p.pa.beta ~ dbeta(1, 1) # probability of availability
#     pem.beta <- logit(pem.beta) # backtransform to probability
#     pem.beta ~ dbeta(1, 1) # probability of emigration
#     psi.beta <- logit(p.psi.beta) # backtransform to probability
#     p.psi.beta ~ dbeta(1, 1) probability of habitat suitability
#     lam.beta ~ dnorm(0, 0.01) # abundance log scale
#     
#     ##### REMOVAL #####################################
#     # Removal submodel
#     for (l in 1:L) {
#     int[l] ~ dcat(pi.pa.c[site[l], visit[l], yr[l], ]) # removal class frequencies
#     } # L
#     
#     for (i in 1:nsites) {
#     for (j in 1:nvisits) {
#     for (t in 1:YR) {
#     for (r in 1:R){
#     pi.pa[i,j,t,r] <- p.a[i,j,t]*pow(1-p.a[i,j,t], (r-1))
#     pi.pa.c[i,j,t,r] <- pi.pa[i,j,t,r] / pcap[i,j,t]
#     }  #R
#     pcap[i,j,t] <- sum(pi.pa[i,j,t,1:R])
#     
#     # Detection submodels 
#     pmarg[i,j,t] <-  pcap[i,j,t]  * pem[i,j,t]
#     logit(p.a[i,j,t]) <- pa.beta # add covariates for availability (time-removal) here
#     logit(pem[i,j,t]) <- pem.beta
#     
#     ##### POINT-LEVEL ABUNDANCE ###########################
#     # Abundance/suitability submodels
#     # N is the abundance which includes any GOEA using the site, ie emigrants
#     nobs[i,j,t] ~ dbin(pmarg[i,j,t], N[i,t]) 
#     }}} # i j t 
#     
#     N[i,t] ~ dpois(lambda.eff[i,t])
#     lambda.eff[i,t] <- lambda[i,t] * w.lam[i,t]
#     w.lam[i,t] ~  dbern(psi[i,t])
#     log(lambda[i,t]) <- lam.beta # add covariates for abundance
#     logit(psi[i,t]) <- psi.beta # add covariates for habitat suitability
#     # If not running set inits near psi=1 on logit scale. For example psi.beta=10
#     
#     ##### GOODNESS OF FIT #######################################
#     nobs.fit[i,j,t] ~ dbin(pmarg[i,j,t], N[i,t]) # create new realization of model
#     e.p[i,j,t] <- pmarg[i,j,t] * N[i,t] # original model prediction
#     E.p[i,j,t] <- pow((nobs[i,j,t]- e.p[i,j,t]),2)/(e.p[i,j,t]+0.5)
#     E.New.p[i,j,t]<- pow((nobs.fit[i,j,t]-e.p[i,j,t]),2)/(e.p[i,j,t]+0.5)
#     }} #YR #nsites 
#     fit.p <- sum(E.p[1:nsites, 1:nvisits, 1:YR])
#     fit.new.p <- sum(E.New.p[1:nsites, 1:nvisits, 1:YR])
#     bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
#     
#     ##### DERIVED QUANTITIES ####################################
#     for(t in 1:YR){
#     Ntot[t] <- sum(N[1:nsites,t])
#     # You can calculate density if you have a quantifiable area of survey
#     # D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
#     } #YR
#     } # End model
#     ",file="./model-basic-ZIP.txt")
# 
# ###################
# # CORRECTION
# ###################
# # add more bins across distance to better 
# # approximate the integral
# midpt.temp <- seq(0,50, length=21)
# midpt2 <- c(1.25)
# for (z in 2:length(midpt.temp)){ midpt2[z] <-  midpt2[z-1]+2.5 }
# datalfoc$midpt <- midpt2 <- midpt2[-21]
# datalfoc$delta <- rep(2.5, length(midpt2))
# datalfoc$nD <- length(midpt2)
# # END CORRECTION 
# 
# for (i in 1:19){ 
#   spp <- spp.list.foc[i]
#   spp.num<- which(dimnames(nobs)[[3]]==spp)
#   datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
#   Mst <- apply(Nav, c(1), max, na.rm=T) +1
#   
#   inits <- function(){  list(
#     N = Nav, p.psi.beta= 0.99, # setting psi near 1 helps model run
#     p.pa.beta0= runif(1, 0.1, 0.9),
#     pp.beta= runif(1, 5, 100)
#   )  }
#   
#   params <- c("pa.beta", "p.pa.beta",
#               "pem.beta", "p.pem.beta",
#               "lam.beta", "psi.beta", "p.psi.beta",
#               "Ntot", "D",  "bayesp"
#   )
#   
#   # MCMC settings
#   ni <- 2000  ;   nb <- 1000   ;   nt <- 1   ;   nc <- 3 ; na=100
#   # Run JAGS
#   out <- jags(datalfoc, inits=inits, params, 
#                "./model-basic-ZIP.txt", 
#               n.thin=nt, n.chains=nc, 
#               n.burnin=nb, n.iter=ni, n.adapt=na,
#               parallel = T, modules=c("glm")
#               )
#   
#   fn<- paste( "./", spp, "_basic-ZIP.RData", sep="" )
#   save(list= c("out"), file=fn)
# }
# 
