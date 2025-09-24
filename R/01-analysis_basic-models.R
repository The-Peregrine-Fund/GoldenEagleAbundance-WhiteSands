
# This file is not needed to replicate analyses.
# We provide the simplest version of the model used
# as a template for other studies
#############################################
# software used
# JAGS 4.3.0 
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
library (jagsUI) # v1.5.1
# load the wide data from GOEA White Sands project
#df.w <- read.csv("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\GoldenEagleDensity-WhiteSands\\data\\surveydatawide.csv")
df.ab <- read.csv("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\GoldenEagleDensity-WhiteSands\\data\\abundancenmixclean.csv")

counts <- with(df.ab, tapply(count, list(site, rotation, survey_year), sum, na.rm=T) )
tri <- with(df.ab, tapply(tri, list(site), mean ) )
obs <- with(df.ab, tapply(multi_observer, list(site, rotation, survey_year), sum) )
day <- with(df.ab, tapply(survey_day, list(site, rotation, survey_year), mean) )

datl <- list(counts = counts[,1:5,3:5],
             mult_surveyors = ifelse(obs[,1:5,3:5]==1, 1, -1 ),
             day.sc = ((day-mean(day, na.rm=T))/sd(day, na.rm=T))[,1:5,3:5],
             tri.sc= (tri-mean(tri, na.rm=T))/sd(tri, na.rm=T),
             nvisits = 5,
             nyears = 3, 
             nsites = nrow(counts))
datl$day.sc[is.na(datl$day.sc)] <- 0
datl$tri.sc[is.na(datl$tri.sc)] <- 0
datl$mult_surveyors[is.na(datl$mult_surveyors)] <- -1
####################################
# (1) N-mixture model for abundance
# repeated survey 
####################################
cat("
    model {
    # pbeta: prob. detection
    # lam.beta: abundance (log-link scale, use exp function to backtransform)
    ##### PRIORS ###############################################
    lpbeta <- logit(pbeta) # backtransform to probability pbeta
    pbeta ~ dbeta(1, 1) # probability of emigration
    for (t in 1:nyears){
        lam.beta[t] ~ dnorm(0, 0.01) # average abundance on log scale
    }
    beta1 ~ dnorm(0, 0.01) # coefficient from num_surveyors
    beta2 ~ dnorm(0, 0.01) # coefficient from day of year
    beta3 ~ dnorm(0, 0.01) # coefficient from topographical roughness index
    sd.site.p ~ dnorm(0, 1/(2*2))T(0,) # random effect for site
    sd.site.N ~ dnorm(0, 1/(2*2))T(0,) # random effect for site
    
    for (i in 1:nsites) {
    for (t in 1:nyears) {
    # Abundance (inserted here to reduce number of loops)
    N[i,t] ~ dpois(lambda[i,t])
    log(lambda[i,t]) <- lam.beta[t] + # average abundance varies by year
                          beta3*tri.sc[i] + # tri effect
                          eps[i] # random effect for site
    
    # Separate detection and abundance    
    for (j in 1:nvisits) {
    counts[i,j,t] ~ dbin(p[i,j,t], N[i,t])
    # Detection submodel
    logit(p[i,j,t]) <- lpbeta + # average detection
                     beta1*mult_surveyors[i,j,t] + # mult observer effect
                     beta2*day.sc[i,j,t] + # day effect
                     eta[i] # random effect for site
    }}}
    # random effects for site
    for (i in 1:nsites){ 
    eta[i] ~ dnorm(0, sd.site.p) 
    eps[i] ~ dnorm(0, sd.site.N)
    }
    
    ##### GOODNESS OF FIT #######################################
    for (i in 1:nsites) {
    for (t in 1:nyears) {
    for (j in 1:nvisits) {
    counts.fit[i,j,t] ~ dbin(p[i,j,t], N[i,t]) # simulate new data realization from the model
    e.p[i,j,t] <- p[i,j,t] * N[i,t] # original model prediction
    E.p[i,j,t] <- pow((counts[i,j,t]- e.p[i,j,t]),2)/(e.p[i,j,t]+0.5) # calculate discrepancy between data and e.p
    E.New.p[i,j,t]<- pow((counts.fit[i,j,t]-e.p[i,j,t]),2)/(e.p[i,j,t]+0.5) # calculate discrepancy between simulations and e.p
    }}} #nyr #nsites
    fit.p <- sum(E.p[1:nsites, 1:nvisits, 1:nyears])
    fit.new.p <- sum(E.New.p[1:nsites, 1:nvisits, 1:nyears])
    bayesp <- step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    } # End model
    ",file="./N-mix-repeatedsurveys.txt")

# We need good initial values for 
# N to get the model running, so
# we estimate it from the data.
N.inits <- apply(datl$counts, c(1, 3), max, na.rm=T) 
N.inits[N.inits=="-Inf"] <- 2

# Create afunction to provide initial values
# for the model. Each chain will randomly 
# generate initial values.
inits <- function(){  list(
  N = N.inits,
  lam.beta = rep(mean(datl$counts, na.rm=T) |> log(),3),
  pbeta= runif(1, 0.01, 0.99),
  beta1=runif(1,-1,1), 
  beta2=runif(1,-1,1),
  beta3=runif(1,-1,1),
  sd.site.p=runif(1),
  sd.site.N=runif(1)
)  }

params <- c("pbeta", "lpbeta",
            "lam.beta",
            "beta1", "beta2", "beta3",
            "bayesp", # goodness of fit measure, should be near 0.5
            "sd.site.p", "sd.site.N",
             "eta", "eps", "N"
)

# MCMC settings
ni <- 300000  ;   nb <- 200000   ;   nt <- 100   ;   nc <- 3 ; na=1000
#ni <- 200  ;   nb <- 100   ;   nt <- 1   ;   nc <- 3 ; na=100
# Run JAGS
out <- jags(datl, inits=inits, params, 
            "./N-mix-repeatedsurveys.txt", 
            n.thin=nt, n.chains=nc, 
            n.burnin=nb, n.iter=ni, n.adapt=na,
            parallel = T, modules=c("glm")
)

traceplot(out, c("pbeta","lambda", 
                 "beta1", "beta2", "beta3",
                 "sd.site.p", "sd.site.N") )

# View average abundance ests at sites
Nmn<- apply(out$sims.list$N, c(2,3), mean)
hist(Nmn, main="Histogram of model-predicted abundance")

fn<- paste( "./outputs/N-mix-repeatedsurveys.RData", sep="" )
save(list= c("out"), file=fn)


####################################
# (2) N-mixture model for abundance
# removal/ repeated survey 
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

#********************
#* Royle-Nichols Model
#* Abundance from detection/nondeteciton data
#**********************
  cat("
    model {
    # p.beta: prob. of detection
    # lam.beta: abundance (log-link scale, use exp function to backtransform)
    ##### PRIORS ###############################################

    lp.beta <- logit(p.beta)
    p.beta ~ dbeta(1,1)
    lam.beta ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    sd.year ~ dnorm(0, 1/(2*2))T(0,) # translates tau to sd=2
    sd.site ~ dnorm(0, 1/(2*2))T(0,)
    sd.both ~ dnorm(0, 1/(2*2))T(0,)
    
    # Royle Nichols model
    # observation model
      for (i in 1:ndata) {
        for (j in 1:nvisits) {
          y[i,j] ~ dbern(Pstar[i,j])
          Pstar[i,j] <- 1-(1-p[i,j])^N[i]
          logit(p[i,j]) <- lp.beta + 
                            beta1*num_surveyors[i,j] +
                            beta2*oday.sc[i,j]
    } # j
    # abundance model
        N[i] ~ dpois(lambda[i])
        log(lambda[i]) <- lam.beta + 
                          eps[ year[i] ] +
                          eta[ site[i] ] +
                          nu[ year[i], site[i]]
    } # i
    # random effect for year
    for (t in 1:nyears){
    eps[t] ~ dnorm(0, sd.year)
    }
    for (s in 1:nsites){
    eta[s] ~ dnorm(0, sd.site)
        for (t in 1:nyears){
          nu[ t, s] ~ dnorm(0, sd.both)
    }}
    } # End model
    ",file="./model-Royle-Nichols.txt")
  
  inits <- function(){  list(
    # inits for lam.beta and p.beta
    # taken from run of model with 
    # no covariates
    # to help run
    lam.beta = 1.168, 
    p.beta= 0.388,
    beta1= 0, 
    beta2= 0,
    sd.year=runif(1, 0.01, 0.1),
    sd.site=runif(1, 0.01, 0.1),
    sd.both=runif(1, 0.01, 0.1)
  )  }
  
  params <- c("p.beta", "lp.beta",
              "lam.beta", 
              "beta1", "beta2",
              "sd.year", "sd.site", "sd.both",
              "eps", "eta", "nu"
  )
  
  # MCMC settings
  ni <- 100000  ;   nb <- 50000   ;   nt <- 100   ;   nc <- 3 ; na=1000
  # Run JAGS
  out <- jags(datl, inits=inits, params, 
              "./model-Royle-Nichols.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm"))
  
  params2 <- c("p.beta", 
              "lam.beta", 
              "beta1", "beta2",
              "sd.year", "sd.site", "sd.both"
  )
  traceplot(out, params2)
  

# Postprocess to estimate abundance
# at each territory, for each year
N <- array(NA, dim=list(datl$nsite, 
                        datl$nyear, 
                        length(out$sims.list$lam.beta)))
for (t in 1:datl$nyear){
  for (s in 1:datl$nsite){
  N[s, t, ] <- exp(out$sims.list$lam.beta + 
             out$sims.list$eps[,t] +
             out$sims.list$eta[,s] +
             out$sims.list$nu[,t,s]) 
  }
}

N.mean <- apply(N, c(1,2), mean)
hist(N.mean)
N.hdis <- apply(N, c(1,2), HDInterval::hdi)

fn<- paste( "./outputs/Royle-Nichols.RData", sep="" )
save(list= c("out"), file=fn)
  