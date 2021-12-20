
############################################################################
# Simulations -------------------------------------------------------------

rm(list=ls())

# Simulations baseline (Poisson & rfixed) ---------------------------------

set.seed(221)

design = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)#design of simulated surveys

Ab = c(0.05, 0.10, 0.20, 1, 2, 4, 6)#range of abundances
1 - dpois(0, Ab)
r = c(0.15, 0.45, 0.75)#range of detectability
np = 2
spp_params = expand.grid(Ab, r, np)

params = NULL
for(i in 1:nrow(design))
{
  tmp = matrix(rep(design[i,], nrow(spp_params)), ncol = 2, byrow = T)
  out= cbind(tmp, spp_params)
  params = rbind(params, out)
}
params = data.frame(params)
colnames(params) = c("n.su", "nv", "Ab", "r", "np")

#outputs
Ys <- vector("list", nrow(params))#Observations
Xs <- vector("list", nrow(params))#Covariates
Betas<- vector("list", nrow(params))#Species-specific responses
Ns<- vector("list", nrow(params))#Number of individuals in each site
Lambdas <- vector("list", nrow(params))#Expected number of individuals
Rs  <- vector("list", nrow(params))#Individual detection
Ps<- vector("list", nrow(params))#Probability of species detection
Psis_L<- vector("list", nrow(params))#Probability of occurrence (based on lambda)
As<- vector("list", nrow(params))#Baseline abundances


nreps = 100
for(pp in 1:nrow(params))
{
  Ab = params$Ab[pp]
  n.su = params$n.su[pp]
  nv = params$nv[pp]
  r = params$r[pp]
  np = params$np[pp]
  #outputs per paramset
  Y = array(NA, c(n.su, nv, nreps))
  X = array(NA, c(n.su, np, nreps))
  B = array(NA, c(nreps, np))
  N = array(NA, c(n.su, nreps))
  L = array(NA, c(n.su, nreps))
  P = array(NA, c(n.su, nreps))
  R =array(NA, c(n.su, nreps))
  Psi_L = array(NA, c(n.su, nreps))
  A =rep(NA,  nreps)
  for (rr in 1:nreps)
  {
    y = array(NA, c(n.su, nv))
    tmp = c(rep(1, n.su), rnorm(n.su*(np-1), 0, 1))
    x = matrix(tmp, nrow = n.su, ncol = np)
    if(runif(1, 0,1) <= 0.5)#mixed-normal (mean -1/1, sd = 1)
    {
      betas = c(log(Ab), rnorm(1, 1, 1)) 
    }else{
      betas = c(log(Ab), rnorm(1, -1, 1))
    }
    lambda = exp(as.vector(betas %*% t(x)))
    n = rep(NA, n.su)
    p = rep(NA, n.su)
    psi_l = rep(NA, n.su)
    for(i in 1:n.su)
    {
      n[i]  = rpois(1, lambda[i]) #sample indiviuals
      psi_l[i] = 1 - dpois(0, lambda[i])#calculate probability of occurence
      p [i] = 1 - (1-r)^n[i]#probability of species detection
      y[i,1:nv] = rbinom(nv, 1, p[i])#observed data
    }
    #save data
    Y[,,rr] = y
    X[,,rr] = x
    B[rr,] = betas
    N[,rr] = n
    L[,rr] = lambda
    P[,rr] = p
    R[,rr] = r
    Psi_L[,rr] = psi_l
    A[rr] = Ab
  }
  Ys[[pp]] = Y
  Xs[[pp]] = X
  Betas[[pp]] = B
  Ns[[pp]] = N
  Lambdas[[pp]] = L
  Ps[[pp]] = P
  Rs[[pp]] = R
  Psis_L[[pp]] = Psi_L
  As[[pp]] = A
}

save.image("Simulations_baseline.RData")


# Simulations scenario 2 :Poisson & r beta ---------------------------------
rm(list=ls())

phi = function(CV, mu){
  A = 1 - mu
  B = CV^2*mu
  C = B-A
  phi = - (C/B)
  return(phi)
}#Function to find phi according to a CV and a mu (beta distribution for individual detection)

CV = 0.5

set.seed(221)

design = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)


Ab = c(0.05, 0.10, 0.20, 1, 2, 4, 6)
1 - dpois(0, Ab)
r = c(0.15, 0.45, 0.75)
np = 2
spp_params = expand.grid(Ab, r, np)

params = NULL
for(i in 1:nrow(design))
{
  tmp = matrix(rep(design[i,], nrow(spp_params)), ncol = 2, byrow = T)
  out= cbind(tmp, spp_params)
  params = rbind(params, out)
}
params = data.frame(params)
colnames(params) = c("n.su", "nv", "Ab", "r", "np")

#outputs

Ys <- vector("list", nrow(params))
Xs <- vector("list", nrow(params))
Betas<- vector("list", nrow(params))
Ns<- vector("list", nrow(params))
Lambdas <- vector("list", nrow(params))
Rs  <- vector("list", nrow(params))
Ps<- vector("list", nrow(params))
Psis_L<- vector("list", nrow(params))
As<- vector("list", nrow(params))

nreps = 100
for(pp in 1:nrow(params))
{
  Ab = params$Ab[pp]
  n.su = params$n.su[pp]
  nv = params$nv[pp]
  r = params$r[pp]
  np = params$np[pp]
  #outputs per paramset
  Y = array(NA, c(n.su, nv, nreps))
  X = array(NA, c(n.su, np, nreps))
  B = array(NA, c(nreps, np))
  N = array(NA, c(n.su, nreps))
  L = array(NA, c(n.su, nreps))
  P = array(NA, c(n.su, nreps))
  R = array(NA, c(n.su, nreps))
  Psi_L = array(NA, c(n.su, nreps))
  A =rep(NA,  nreps)
  for (rr in 1:nreps)
  {
    y = array(NA, c(n.su, nv))
    tmp = c(rep(1, n.su), rnorm(n.su*(np-1), 0, 1))
    x = matrix(tmp, nrow = n.su, ncol = np)
    if(runif(1, 0,1) <= 0.5)
    {
      betas = c(log(Ab), rnorm(1, 1, 1)) 
    }else{
      betas = c(log(Ab), rnorm(1, -1, 1))
    }
    lambda = exp(as.vector(betas %*% t(x)))
    
    phi_ = phi(CV = CV, mu = r)#Find phi according to the CV and r
    alpha = r*phi_
    beta = (1 - r)*phi_
    r_ = rbeta(n.su, alpha, beta)#sample site-specific detections
    n = rep(NA, n.su)
    p = rep(NA, n.su)
    psi_l = rep(NA, n.su)
    for(i in 1:n.su)
    {
      n[i]  = rpois(1, lambda[i]) 
      psi_l[i] = 1 - dpois(0, lambda[i])
      p [i] = 1 - (1-r_[i])^n[i]
      y[i,1:nv] = rbinom(nv, 1, p[i])
    }
    #save data
    Y[,,rr] = y
    X[,,rr] = x
    B[rr,] = betas
    N[,rr] = n
    L[,rr] = lambda
    P[,rr] = p
    R[,rr] = r_
    Psi_L[,rr] = psi_l
    A[rr] = Ab
  }
  Ys[[pp]] = Y
  Xs[[pp]] = X
  Betas[[pp]] = B
  Ns[[pp]] = N
  Lambdas[[pp]] = L
  Ps[[pp]] = P
  Rs[[pp]] = R
  Psis_L[[pp]] = Psi_L
  As[[pp]] = A
}

save.image("Simulations_Poisson_beta.RData")



# Scenario 3: Negbinom ---------------------------------------------


rm(list=ls())
phi_nb = function(N, mu){
  A = mu^2
  B = mu*(N-1)
  phi =  A/B
  
  return(phi)
}#find phi of negbinom according Var(X) = N*mu

Nmu = 2#Overdisperssion for local abundances


set.seed(221)

design = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)
#Best , Intermediate and "Common" case scenarios of sampling

Ab = c(0.05, 0.10, 0.20, 1, 2, 4, 6)
1 - dpois(0, Ab)
r = c(0.15, 0.45, 0.75)
np = 2
spp_params = expand.grid(Ab, r, np)

params = NULL
for(i in 1:nrow(design))
{
  tmp = matrix(rep(design[i,], nrow(spp_params)), ncol = 2, byrow = T)
  out= cbind(tmp, spp_params)
  params = rbind(params, out)
}
params = data.frame(params)
colnames(params) = c("n.su", "nv", "Ab", "r", "np")

#Simulations

Ys <- vector("list", nrow(params))
Xs <- vector("list", nrow(params))
Betas<- vector("list", nrow(params))
Ns<- vector("list", nrow(params))
Lambdas <- vector("list", nrow(params))
Rs  <- vector("list", nrow(params))
Ps<- vector("list", nrow(params))
Psis_L<- vector("list", nrow(params))
As<- vector("list", nrow(params))

nreps = 100
for(pp in 1:nrow(params))
{
  Ab = params$Ab[pp]
  n.su = params$n.su[pp]
  nv = params$nv[pp]
  r = params$r[pp]
  np = params$np[pp]
  #outputs per paramset
  Y = array(NA, c(n.su, nv, nreps))
  X = array(NA, c(n.su, np, nreps))
  B = array(NA, c(nreps, np))
  N = array(NA, c(n.su, nreps))
  L = array(NA, c(n.su, nreps))
  P = array(NA, c(n.su, nreps))
  R = array(NA, c(n.su, nreps))
  Psi_L = array(NA, c(n.su, nreps))
  A =rep(NA,  nreps)
  for (rr in 1:nreps)
  {
    y = array(NA, c(n.su, nv))
    tmp = c(rep(1, n.su), rnorm(n.su*(np-1), 0, 1))
    x = matrix(tmp, nrow = n.su, ncol = np)
    if(runif(1, 0,1) <= 0.5)#
    {
      betas = c(log(Ab), rnorm(1, 1, 1)) 
    }else{
      betas = c(log(Ab), rnorm(1, -1, 1))
    }
    lambda = exp(as.vector(betas %*% t(x)))
    phi_nb_ = phi_nb(N =Nmu, mu= Ab)#dispersion parameter according to overdisperssion and mean values
    #Observation process   
    n = rep(NA, n.su)
    p = rep(NA, n.su)
    psi_l = rep(NA, n.su)
    for(i in 1:n.su)
    {
      n[i]  = rnbinom(1, mu = lambda[i], size = phi_nb_) 
      psi_l[i] = 1 - dnbinom(0, size = phi_nb_, mu = lambda[i])#probability of occurence
      p [i] = 1 - (1-r)^n[i]
      y[i,1:nv] = rbinom(nv, 1, p[i])
    }
    #save data
    Y[,,rr] = y
    X[,,rr] = x
    B[rr,] = betas
    N[,rr] = n
    L[,rr] = lambda
    P[,rr] = p
    R[,rr] = r
    Psi_L[,rr] = psi_l
    A[rr] = Ab
  }
  Ys[[pp]] = Y
  Xs[[pp]] = X
  Betas[[pp]] = B
  Ns[[pp]] = N
  Lambdas[[pp]] = L
  Ps[[pp]] = P
  Rs[[pp]] = R
  Psis_L[[pp]] = Psi_L
  As[[pp]] = A
}

save.image("Simulations_Nbinom.RData")

##########################################################################
# Adjust data -------------------------------------------------------------

# Models ------------------------------------------------------------------


# Royle y Nichols 2003 ----------------------------------------------------

cat(file = "RN_model.bug","
model {
  # Likelihood
  for (i in 1:nSite) {
    # Process model
    
    lambda[i] <- exp(b0 + b1*x[i,2])
    N[i] ~ dpois(lambda[i])	# latent abundance of species i
    z[i] <- step(N[i] - 1) # z=1 if N>0, ie. site is occupied
    psi_L[i] <- 1 - exp(-lambda[i])#site-specific occupancy
    # Observation model
    p[i] <- 1 - pow(1-r, N[i])
    y[i] ~ dbin(p[i], nRep)
  }
  # Priors
  r ~ dunif(0, 1)
  b0 ~ dnorm(0, 0.01) 
  b1 ~ dnorm(0, 0.01) 
  # Derived quantities
  psi.sample <- mean(z[])   # Proportion of occupied sites in the sample
}
"
)



# Bernoulli ---------------------------------------------------------------

cat(file = "bernoulli_model.bug", "
 #Likelihood:
    model { 
     for(i in 1: M){ 
       y[i] ~ dbin(pex[i], J[i]) 
       pex[i] <- z[i] * p   
       z[i] ~ dbin(psi[i], 1)      
       logit(psi[i]) <- b0 + b1 * X[i] 
    } 
  #Previas: 
    p  ~ dunif(0,1)                
    b0 ~ dnorm(0, 0.01) 
    b1 ~ dnorm(0, 0.01) 
    psi.sample = sum(z)/M
    }
    ")



# Adjust data ---------------------------------------------------------

rm(list=ls())
require("jagsUI")
load("Simulations_baseline.RData")#To change if Pois-beta or Negative-binom


nmod = 2


for(pp in 1:nrow(params))
{
  print(paste("param", pp/nrow(params)))
  sum_RN = vector("list", nreps)##Save summary RN
  sum_bern = vector("list", nreps)##Save summary MC
  for (rr in 1:nreps)
  {
    tmp = (rr/nreps)*100
    print(paste (tmp, "%"))
    ##load data
    Y = Ys[[pp]][,,rr]
    X = Xs[[pp]][,,rr]
    N = Ns[[pp]][,rr]
    tmp = length(which(N > 100))#Avoid simulations with sites N>100
    if(sum(N) > 0 | length(tmp) == 0)
    {
      ##Chain properties (change according to needs)
      ni <- 50000  # number of iterations
      nt <- 5 # thining
      nb <- 5000  #'burn in'
      nc <- 3  # number of chains
      
      
      # Adjust RN  --------------------------------------------------------------
      data <- list(y = rowSums(Y), nRep = ncol(Y), nSite = nrow(Y), x = X)
      
      Nst <- as.numeric(rowSums(Y) > 0)  # N must be >0 if y > 0, so start at 1
      inits <- function() list(N = Nst,  b0 = rnorm(1), b1 = rnorm(1), 
                               r=runif(1, 0, 1))     
      parameters <- c("b0", "b1", "lambda", "N", "psi_L", "p", "r", "psi.sample")
      
      
      mod_RN<- jags(data, inits, parameters, "RN_model.bug", n.chains = nc, n.thin = nt, 
                    n.iter = ni, n.burnin = nb, parallel = T)
      summary_RN = as.data.frame(mod_RN$summary)
      tmp = which(row.names(summary_RN) %in% c("b0", "b1", "r"))#select model parameters (not derived quantities)
      summary_RN2 = summary_RN [tmp,]
      if(length(which(summary_RN2$Rhat > 1.1)) == 0 & #If converged and Neff>=100
         length(which(summary_RN2$n.eff < 100)) == 0)
      {
        sum_RN[[rr]] = summary_RN
      }
      
      ###########################################################################
      # Adjust Bernoulli model-------------------------------------------
      data <- list (y = rowSums(Y), J = rep(ncol(Y), nrow(Y)), M = nrow(Y), X = X[,2])
      inits <- function() {
        list(z = rbinom(nrow(Y), 1, 1), b0 = rnorm(1), p = runif(1), b1 = rnorm(1))
      }
      
      parameters <- c("p", "b0", "b1", "psi.sample", "psi")
      mod_bern <- jags(data, inits, parameters, "bernoulli_model.bug", n.chains = nc, n.thin = nt, 
                       n.iter = ni, n.burnin = nb, parallel = T)
      summary_mod_bern = as.data.frame(mod_bern$summary)
      tmp = which(row.names(summary_mod_bern) %in% c("p", "b0", "b1"))
      summary_mod_bern2 = summary_mod_bern [tmp,]
      if(length(which(summary_mod_bern2$Rhat > 1.1)) == 0 &
         length(which(summary_mod_bern2$n.eff < 100)) == 0)
      {
        sum_bern[[rr]] = summary_mod_bern
      }
      
      ##Save.RData
      if(rr %% nreps == 0)
      {
        file = paste("fit_data_baseline", pp, "_", rr, ".RData", sep = "")##To change when Pois-beta and Negbinom
        save(sum_RN, sum_bern,
             file = file)#Save summaries of fitted models
      }
    }#If N<100
  }#end repetitions loop
}#end of parameters loop



##########################################################################
# Calculate errors --------------------------------------------------------

rm(list = ls())

load("Simulations_baseline.RData")#Change if Pois-beta or Nbinom
library("forecast")

OUT = array(NA, c(nrow(params), 13))
for(pp in 1:nrow(params))
{
  file = paste("fit_data_baseline", pp, "_100.RData",  sep = "")#Change if Pois-beta or Nbinom
  load(file)
  # Output vector and matrices  ---------------------------------------------
  #1- RN model, 2- MC model
  Power1 = rep(NA, nreps)
  Power2 = rep(NA, nreps)
  RMSE1 = rep(NA, nreps)
  RMSE2 = rep(NA, nreps)
  Reg1 = array(NA, c(nreps,2))#intercept and slope
  Reg2 = array(NA, c(nreps,2))
  #which simulations have N < 100 in some plots (realistic values of abundance)
  id= apply(Ns[[pp]], 2, function(x) length(which(x > 100)))
  id = which(id == 0)
  for(rr in id)
  {
    # Model fit ---------------------------------------------------------------
    # RN  -------------------------------------------------------------
    summary_mod = sum_RN[[rr]]
    if(!is.null(summary_mod))#if it converged
    {
      # Power --------------------------------------------------------------------
      true_beta = round(Betas[[pp]][rr,2], 2)
      est_beta = summary_mod[grepl("b1", rownames(summary_mod)),]
      if(true_beta != 0 & est_beta$overlap0 == 1)
      {Power1[rr] = 1#Type II error
      }else{
        Power1[rr] = 0#
      }
      # Site-specific occupancy (RMSE and regression) --------------------------------------------
      true_psis = Psis_L[[pp]][,rr]
      est_psis = summary_mod[grepl("psi_L", rownames(summary_mod)),]
      tmp = accuracy(est_psis$mean, true_psis)
      RMSE1[rr] = tmp[[2]] 
      mod = lm(est_psis$mean~true_psis)
      Reg1[rr,] = mod$coefficients
    }#end of RN (loop for those that converged and had enough Neff)
    
    # MC ----------------------------------------------------------------------
    summary_mod = sum_bern[[rr]]
    if(!is.null(summary_mod))#if it converged
    {
      # Power --------------------------------------------------------------------
      true_beta = round(Betas[[pp]][rr,2], 2)
      est_beta = summary_mod[grepl("b1", rownames(summary_mod)),]
      if(true_beta != 0 & est_beta$overlap0 == 1)
      {Power2[rr] = 1#Type II error
      }else{
        Power2[rr] = 0#
      }
      # Site-specific occupancy (RMSE and regression) --------------------------------------------
      true_psis = Psis_L[[pp]][,rr]
      est_psis = summary_mod[grepl("psi", rownames(summary_mod)),]
      tmp = which(rownames(est_psis) == "psi.sample" | rownames(est_psis) == "psi0")
      est_psis = est_psis[-tmp,]
      tmp = accuracy(est_psis$mean, true_psis)
      RMSE2[rr] = tmp[[2]] 
      mod = lm(est_psis$mean~true_psis)
      Reg2[rr,] = mod$coefficients
    }#end of RN (loop for those that converged and had enough Neff)
  }#rep loop
  
  # Save outputs ------------------------------------------------------------
  type2_1 = sum(na.omit(Power1))/length(which(!is.na(Power1)))
  power1 = 1-type2_1
  type2_2 = sum(na.omit(Power2))/length(which(!is.na(Power2)))
  power2 = 1-type2_2
  rmse1 = median(na.omit(RMSE1))
  rmse2 = median(na.omit(RMSE2))
  reg1 = apply(Reg1, 2, function(x) quantile(na.omit(x), probs = 0.5))
  reg2 = apply(Reg2, 2, function(x) quantile(na.omit(x), probs = 0.5))
  out = c(params$n.su[pp], params$nv[pp], params$Ab[pp], params$r[pp], "baseline", #Change if Pois-beta or Nbinom
          power1, power2, rmse1, rmse2, reg1, reg2)
  OUT[pp,] = out
}#params loop
OUT = data.frame(OUT)
colnames(OUT) = c("sites", "visits", "Abundance", "r", "Model_type", 
                  "power_RN", "power_MC",
                  "Occ_rmse_RN", "Occ_rmse_MC", 
                  "int_RN", "slope_RN", "int_MC", "slope_MC")
write.csv(OUT, "error_baseline.csv", row.names= F) #Change if Pois-beta or Nbinom


###########################################################################
###MC extension############################################################
###########################################################################

#MC model in which site-specific detection depends on the same environmental
#covariates affecting abundances

# MC extension---------------------------------------------------------------

cat(file = "bernoulli_model_plus.bug", "
 #Likelihood:
    model { 
     for(i in 1: M){ 
       y[i] ~ dbin(pex[i], J[i]) 
       pex[i] <- z[i] * p [i]  
       z[i] ~ dbin(psi[i], 1)      
       logit(psi[i]) <- b0 + b1 * X[i] 
       logit(p[i]) = a0+  a1 * X[i] 
    } 
  #Prios: 
    a0 ~ dnorm(0, 0.01) 
    a1 ~ dnorm(0, 0.01) 
    b0 ~ dnorm(0, 0.01) 
    b1 ~ dnorm(0, 0.01) 
  #Derived quantities
    logit(psi0) <- b0#probabilidad basal de presencia
    psi.sample = sum(z)/M#proporción de sitios ocupados
    }
    ")

# Adjust data ---------------------------------------------------------

rm(list=ls())
require("jagsUI")
load("Simulations_baseline.RData")##To change if Pois-beta or Nbinom


for(pp in 1:nrow(params))
{
  sum_bern = vector("list", nreps)#save sumaries
  for (rr in 1:nreps)
  {
    tmp = (rr/nreps)*100
    print(paste (tmp, "%"))
    ##load data
    Y = Ys[[pp]][,,rr]#detection/not detection data [sites, visits]
    X = Xs[[pp]][,,rr]#covariates (intercept+covariate)
    tmp = length(which(N > 100))#More than 100 individuals per site seem unrealistic
    if(sum(N) > 0 | length(tmp) == 0)
    {
      ##Chain properties
      ni <- 50000  # n?mero de iteraciones
      nt <- 5 # descartamos 'thin' algunos valores
      nb <- 5000  # cuantas iteraciones usamos de 'burn in'
      nc <- 3  # y cuantas cadenas corremos
      
      
      ###########################################################################
      # Adjust  model-------------------------------------------
      data <- list (y = rowSums(Y), J = rep(ncol(Y), nrow(Y)), M = nrow(Y), X = X[,2])
      inits <- function() {
        list(z = rbinom(nrow(Y), 1, 1), b0 = rnorm(1), b1 = rnorm(1), a0 = rnorm(1), a1 = rnorm(1))
      }
      
      parameters <- c("p", "b0", "b1", "a0", "a1", "psi.sample", "psi")
      mod_bern <- jags(data, inits, parameters, "bernoulli_model_plus.bug", n.chains = nc, n.thin = nt, 
                       n.iter = ni, n.burnin = nb, parallel = T)
      summary_mod_bern = as.data.frame(mod_bern$summary)
      tmp = which(row.names(summary_mod_bern) %in% c("p", "b0", "b1", "a0", "a1"))
      summary_mod_bern2 = summary_mod_bern [tmp,]
      if((length(which(summary_mod_bern2$Rhat > 1.1)) == 0 &#converged and adequate Neff 
          length(which(summary_mod_bern2$n.eff < 100)) == 0))
      {
        sum_bern[[rr]] = summary_mod_bern
      }
      
      ##Save .RData
      if(rr %% nreps == 0)
      {
        file = paste("fit_data_baseline_McPlus", pp, "_", rr, ".RData", sep = "")#Cambiar si es PB o Nbinom
        save(sum_bern,
             file = file)
      }
    }#site-specific N > 0 and N<100
  }#nreps loop
}#params loop


# Calculate errors --------------------------------------------------------


rm(list = ls())

load("Simulations_baseline.RData")#Change if Poisson-beta or Nbinom
library("forecast")

OUT = array(NA, c(nrow(params), 13))
for(pp in 1:nrow(params))
{
  file = paste("fit_data_baseline_McPlus", pp, "_100.RData",  sep = "")#Change if Poisson-beta or Nbinom
  load(file)
  # Output vector and matrices  ---------------------------------------------
  Power = rep(NA, nreps)#Power to detect covariate effects
  RMSE = rep(NA, nreps)#RMSE site-specific occupancy
  Reg = array(NA, c(nreps,2))#intercept and slope estimated ~ true site-specific occupancy
  RMSE_p = rep(NA, nreps)#RMSE species detection
  Reg_p = array(NA, c(nreps,2))#intercept and slope est~true species detection (site-specific)
  #which simulations have N < 100 in some plots (realistic values of abundance)
  id= apply(Ns[[pp]], 2, function(x) length(which(x > 100)))
  id = which(id == 0)
  count = 0#count sample size (to increase niter of chains if necessary)
  for(rr in id)
  {
    # Model fit ---------------------------------------------------------------
    # MC plus  -------------------------------------------------------------
    summary_mod = sum_bern[[rr]]
    if(!is.null(summary_mod))#if it converged
    {
      count = count+1
      # Power --------------------------------------------------------------------
       true_beta = round(Betas[[pp]][rr,2], 2)
      est_beta = summary_mod[grepl("b1", rownames(summary_mod)),]
      if(true_beta != 0 & est_beta$overlap0 == 1)
      {Power[rr] = 1#Type II error
      }else{
        Power[rr] = 0#
      }
      
      # Site-specific occupancy (RMSE and regression) --------------------------------------------
      true_psis = Psis_L[[pp]][,rr]
      est_psis = summary_mod[grepl("psi", rownames(summary_mod)),]#select outputs with psis
      tmp = which(rownames(est_psis) == "psi.sample" | rownames(est_psis) == "psi0")#remove outputs of basal and landscape occupancy.
      est_psis = est_psis[-tmp,]
      tmp = accuracy(est_psis$mean, true_psis)#estimate accuracy
      RMSE[rr] = tmp[[2]] #select RMSE
      mod = lm(est_psis$mean~true_psis)#regression of site specific occupancies
      Reg[rr,] = mod$coefficients#save coefficients
      
      # Probability of detection ------------------------------------------------
      true_ps = Ps[[pp]][,rr]
      est_ps = summary_mod[grepl("p", rownames(summary_mod)),]
      tmp = which(grepl("psi", rownames(est_ps)) == TRUE)#remove those with psi
      est_ps = est_ps[-tmp,]
      tmp = accuracy(est_ps$mean, true_ps)
      RMSE_p[rr] = tmp[[2]]#save RMSE
      mod = lm(est_ps$mean ~true_ps)#regress values
      Reg_p[rr,] = mod$coefficients
    }#if converged and Neff good enough
  }#end loop repetition
  type2 = sum(na.omit(Power))/length(which(!is.na(Power)))#calculate typeII error rate
  power = 1-type2
  #calculate median values of the rest of error rates
  rmse = median(na.omit(RMSE))#
  reg = apply(Reg, 2, function(x) quantile(na.omit(x), probs = 0.5))
  rmse_p = median(na.omit(RMSE_p))
  reg_p = apply(Reg_p, 2, function(x) quantile(na.omit(x), probs = 0.5))
  
  out = c(params$n.su[pp], params$nv[pp], params$Ab[pp], params$r[pp], "baseline", #change if PB o Nbinom
          power, rmse, reg, rmse_p, reg_p, count)
  OUT[pp,] = out
  
  # Save outputs ------------------------------------------------------------
}#end parameter loop

OUT = data.frame(OUT)
colnames(OUT) = c("sites", "visits", "Abundance", "r", "Model_type", 
                  "power","Occ_rmse", "int", "slope", "P_rmse", "int_p", "slope_p", "sample_size")
write.csv(OUT, "error_baseline_MCplus.csv", row.names= F)#Change if Poisson-beta or Negative binomial


