
# Simulate data -----------------------------------------------------------
rm(list=ls())
#setwd("D:/tere/11_06_20")
library(ape)
library(mvtnorm)
library(truncnorm)
set.seed(1234)
library(emulator)
library(matrixcalc)
library(Matrix)


params = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)
colnames(params)= c("nsites", "n.v")


n.s <- 20 # number of bird spp
n.p <- 2   # number of parameters that model abundance (b0, b1)
n.t <- 1# number of traits
n.t.p <- 1 #number of traits prob. of detection
n.p.p <- 1 #number of parameters that model prob of detection

nreps = 20

for(pp in 1:nrow(params))
{
  n.su = params[pp,1]
  n.v = params[pp,2]
  n = n.s*n.su
  N_m = array(NA, c(n, nreps))
  L_m = array(NA, c(n, nreps))
  Occ_m = array(NA, c(n.s, nreps))
  P_det = array(NA, c(n.s, nreps))
  Ys = vector("list", nreps)
  Xs = vector("list", nreps)
  Ns = vector ("list", nreps)
  #Parameters that model multi-spececies occupancy
  Rho = rep(NA, nreps)
  Ls = vector("list", nreps)
  Thetas = vector ("list", nreps)
  Rs  =  vector ("list", nreps)
  Zs = vector("list", nreps)
  Traits = vector("list", nreps)
  Cs = vector("list", nreps)
  Sigmas = vector("list", nreps)
  Omegas = vector("list", nreps)
  Taus = vector("list", nreps)
  ##Parameters that model probability of detection 
  P_det = array(NA, c(n.s, nreps))
  Traits2 = vector("list", nreps)
  Z2 = vector("list", nreps)
  Tau2 = rep(NA, nreps)
  Rho2 = rep(NA, nreps)
  count = 1
  while(count <= nreps)
  { print(count)
 # Simulate species environmental responses-------------------

    #traits that modify response covariates
    t1 <- rnorm(n.s, 0, 1)  # 
    # create trait matrix including ones for the intercept
    TT <- as.matrix(cbind(rep(1, n.s), scale(t1)))
    #traits that modify probability of detection
    t2 <- rnorm(n.s, 0, 1)  # 
    # define trait matrix includding ones for the intercept
    TT2 <- as.matrix(cbind(rep(1, n.s), scale(t2)))
    
    # simulate bird phylogeny 
    tree <- rtree(n = n.s)    
    CC <- vcv(tree, corr = TRUE) # species correlation based on phylogeny
    # sort species and re-arrange phylogenetic correlation
    tmp <- dimnames(CC)
    ids <- as.numeric(as.factor(tmp[[1]]))
    C <- matrix(NA, ncol(CC), ncol(CC))
    for(i in 1:ncol(CC)){
      for(j in 1:ncol(CC)){
        C[ids[i],ids[j]] <- CC[i,j]
      }
    }
    
    #Trait effects on abundance & thetas
    if((runif(1, 0, 1) <= 0.5))#Trait effects can be either positive or negative
    {
      tmp = rnorm(1, -1, 0.5)
      tmp2 =  rnorm(1, -1, 0.5)
      Z= matrix(c(log(0.3), 0, tmp2, tmp), nrow= 2, ncol = 2, byrow = FALSE)#
    }else{
      tmp = rnorm(1, 1, 0.5)
      tmp2 =  rnorm(1, 1, 0.5)
      Z= matrix(c(log(0.3), 0, tmp2, tmp), nrow= 2, ncol = 2, byrow = FALSE)#
    }
    
    M <- TT %*% Z#Obtain expected responses
    cor = rnorm(1, 0, 0.25)#Correlation between parameters
    omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= n.p, nrow= n.p)
    taus = rep(runif(1, 0.1, 0.5) , 2)
    Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= n.p, nrow= n.p)
    Sigma = round(Sigma, 4)
    if(is.symmetric.matrix(Sigma) == FALSE)#Evaluate if it symmetric and positive-definite (is a variance-covariance matrix)
    {
      while(is.symmetric.matrix(Sigma) == FALSE )
      { cor = rnorm(1, 0, 0.5)
      omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= n.p, nrow= n.p)
      taus = runif(2, 0.1, 0.5)
      Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= n.p, nrow= n.p)   
      }#symmetric condition
      if(is.positive.definite(Sigma) == FALSE)
      {
        while(is.positive.definite(Sigma) == FALSE )
        { cor = rnorm(1, 0, 0.5)
        omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= n.p, nrow= n.p)
        taus = runif(2, 0.1, 0.5)
        Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= n.p, nrow= n.p)   
        }
      }#positive.definite condition
    }else{
      if(is.positive.definite(Sigma) == FALSE)
      {
        while(is.positive.definite(Sigma) == FALSE )
        { cor = rnorm(1, 0, 0.5)
        omega = matrix(c(1, cor, cor, 1), byrow = T, ncol= n.p, nrow= n.p)
        taus = runif(2, 0.1, 0.5)
        Sigma = matrix(as.numeric(quad.form(omega, diag(taus))),  byrow = T, ncol= n.p, nrow= n.p)   
        }
      }#positive.definite condition
    }#end if-else
    

# Simulate species detectaiblity ------------------------------------------
    rho = runif(1, 0.2, 0.8)  # sample phylogeny effects
    thetas <- rmvnorm(1, mean = as.vector(M), kronecker(Sigma, rho*C + (1-rho) * diag(n.s)))#sample species-specific responses
    Theta <- matrix(thetas[1,], n.s, n.p)
    
    #Trait effects on probability of detection 
    if((runif(1, 0, 1) <= 0.5))#They can be positive or negative
    {
      tmp2 =  rnorm(1, -1, 0.5)
    }else{
      tmp2 =  rnorm(1, 1, 0.5)
    }
    Zp= matrix(c(qlogis(0.3) , tmp2), nrow= 2, ncol = 1, byrow = FALSE)#Matrix of trait effects
    M2 =  TT2 %*% Zp
    taup = runif(1, 0.1, 0.5)#variance of probability of detection
    rhop = runif(1, 0.2, 0.8) #phylogeny effects
    p.ds = rmvnorm(1, mean = as.vector(M2), taup*(rhop*C+(1-rhop)*diag(n.s)))#sample species-specific detection
    p.ds = as.vector(p.ds)
    r = plogis(p.ds)
    
    
    # Simulate data ------------------------------------------------------
    
    x = array(NA, c(n.su, n.p))
    x[,n.p] = rnorm(n.su, 0, 1)#Covariate effects
    x[,1] = rep(1, nrow(x))#X matrix
    
    # species indicator & extended x matrix
    j <- rep(1:n.s, each = n.su) 
    x.s = NULL
    for(i in 1:n.s)
    {
      x.s = rbind(x.s, x)
    }
    n <- length(j) # sample size
    lambda <- exp(Theta[j,1] + Theta[j,2] * x.s[,2])##
    r2 = r[j]
    
    #Simulate observation process
    N = rep(NA, n)
    p = rep(NA, n)
    y <- matrix(NA, n, n.v)
    for(i in 1:n)
    {
      N[i]  = rpois(1, lambda[i]) 
      p [i] = 1 - (1-r2[i])^N[i]
      y[i,1:n.v] = rbinom(n.v, 1, p[i])
    }
    Occ = rep(NA, n.s)#summarize species occupancies (to select rare and common)
    for(i in 1:n.s)
    {
      init = ((i-1)*n.su)+1
      fin = i*n.su
      tmp = N[init:fin]
      Occ[i] = length(which(tmp>0))/n.su
    }
    tmp = which(round(Occ,1) <= 0.10)
    tmp2 = which(round(Occ,1) >= 0.5)
    if(length(tmp)>0 & length(tmp2)>0 & max(N)<100)#If there are common and rare species and abundances are realistic
    {
      
      N_m[,count] = N
      L_m[,count] = lambda
      Occ_m[,count] = Occ
      P_det[,count] = p.ds
      Ys[[count]] = y
      Xs[[count]] = x.s
      Thetas[[count]] = Theta
      Rho[count] = rho
      Zs[[count]] = Z
      Traits[[count]] = TT
      Cs[[count]] = C
      Sigmas[[count]] = Sigma
      Omegas[[count]]= omega
      Taus[[count]] = taus
      P_det[,count] = r
      Traits2[[count]] = TT2
      Z2[[count]] = Zp
      Tau2[count] = taup
      Rho2[count] = rhop
      count = count+1
    }
    
  }#end loop nreps
  file = paste("sim_", pp, "F.RData", sep = "")
  save(N_m, L_m, Occ_m, P_det, Ys, Xs, Thetas, Rho, Zs, Traits, Cs,
       Sigmas, Omegas, Taus, params, 
       n.s, n.su, n.v, P_det, Traits2, Z2, Tau2, Rho2, file = file)
}#end loop params


# Adjust JSM --------------------------------------------------------------

rm(list = ls())
library("rstan")

# Stan model --------------------------------------------------------------

cat(file = "RNJSM.stan", "
functions { 
    /* compute the kronecker product
    * Args: 
    *   A,B: matrices 
    * Returns: 
    *   kronecker product of A and B
    */ 
    matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:cols(A)) { 
    for (j in 1:rows(A)) { 
    kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
    } 
    } 
    return kron; 
    } 
} 
    data { 
  int<lower=1> N;                 // total number of observations 
  int<lower=1> T;                 // Number of (potential) replicates at each site
  int<lower=0,upper=1> Y[N,T];    // response variable 
  int<lower=1> K;                 // number of sample-level predictors 
  int<lower=1> N_J;               // num of groups (spp) 
  int<lower=1> L_J;               // num group level predictors (traits)
  int<lower=1,upper=N_J> J[N];    // group id 
  matrix[N, K] X;                // obs-level design matrix 
  matrix[N_J, L_J] TT;       // group-level traits for presence
  matrix[N_J, N_J] C;             // phylogenetic correlation matrix
  vector[N_J] ones;               // vector on 1s
  int<lower=0,upper=1> yz[N];     // at least one detection
  int<lower=0> n_max;             // Upper bound of population size
  int<lower=1> Kp;                // 1 (detection probability)
  int<lower=1> L_Jp;              // num group level predictors (traits) for detection probability
  matrix[N_J, L_Jp] TTp;          // group-level traits for detection probability
}


parameters {
  corr_matrix[K] Omega;     // correlation matrix for var-covar of betas 
  vector<lower=0>[K] tau;   // scales for the variance covariance of betas 
  vector[L_J * K] z;        // coeffs for trait effects on species environmental responses
  vector[N_J * K] betas;    // species environmental responses
  real<lower=0,upper=1> rho;  // phylogenetic effects on environmental responses
  real<lower=0> taup;         // scales for the variance covariance of p
  vector[N_J*Kp] r;          //vector with species-specific probability of detection
  vector[L_Jp * Kp] zp;       // coeffs for traits effects on detectability 
  real<lower=0,upper=1> rhop; // phylogenetic effects on detectability
}


transformed parameters { 
  matrix[K, K] Sigma = quad_form_diag(Omega, tau);//Variance-covariance matrix of parameters
  matrix[N_J*K, N_J*K] S = kronecker(Sigma, rho * C + (1-rho) * diag_matrix(ones));//Variance-covariance multivar-ate normal
  matrix[L_J, K] Z = to_matrix(z, L_J, K); //
  vector[N_J * K] m = to_vector(TT * Z);          // Expected species environmental reponses 
  matrix[N_J, K] b_m = to_matrix(betas, N_J, K);  // 
 matrix[N_J, N_J] Sp = taup *(rhop * C + (1-rhop) * diag_matrix(ones));//variance-covariance for detectability
  matrix[L_Jp, Kp] Zp = to_matrix(zp, L_Jp, Kp);  // 
  vector[N_J * Kp] mp = to_vector(TTp * Zp);      // Expected species detectability
  matrix[N_J, Kp] p_m = to_matrix(r, N_J, Kp);   //
  real lambda[N];
  for (n in 1:N){
    lambda[n]  = exp(b_m[J[n], 1] +b_m[J[n],2]* X[n,2]); // Espected local abundances per species
  }
  }
  

model {
  // priors & sampling species-specific responses
   Omega ~ lkj_corr(2);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal() (other possibilities)
  betas ~ multi_normal(m, S); // sampling species-specific environmental responses
  z ~ normal(0,1);
  taup ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()(other possibilities)
  r~multi_normal(mp, Sp); // sample species-specific detectability
  zp ~ normal(0,1);
  //prior for rho & rhop
  target += log_sum_exp(log(0.5) +  beta_lpdf(rho|1, 100), log(0.5) +  beta_lpdf(rho|2, 2));
  target += log_sum_exp(log(0.5) +  beta_lpdf(rhop|1, 100), log(0.5) +  beta_lpdf(rhop|2, 2));
//Likelihood
   for (n in 1:N){
    vector[n_max - yz[n] + 1] lp;
    if(yz[n] == 0){
      lp[1] = poisson_lpmf(0 | lambda[n]) + 1;
    }
    else lp[1] = poisson_lpmf(1 | lambda[n]) + binomial_lpmf(Y[n] | 1, inv_logit(p_m[J[n],1]));
    for (j in 2:(n_max - yz[n] + 1)){
      lp[j] = poisson_lpmf(yz[n] + j - 1 | lambda[n]) 
      + binomial_lpmf(Y[n] | 1, 1 - (1 - inv_logit(p_m[J[n],1]))^(yz[n] + j - 1) );
    }
    target += log_sum_exp(lp);
  }
  }
  
")


# Adjust simulated data ---------------------------------------------------

rm(list=ls()) 
library("rstan")
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)

params = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)
colnames(params)= c("nsites", "n.v")



nreps = 20
for(pp in 1:nrow(params))
{
  file = paste("sim_", pp, "F.RData", sep = "")
  load(file)
  for(rr in 1:nreps)
  {
    y = Ys[[rr]]
    x = Xs[[rr]]
    TT = Traits[[rr]]
    C = Cs[[rr]]
    yz = rowSums(y)
    yz[yz>0] = 1
    j <- rep(1:n.s, each = n.su) 
    TT_det = Traits2[[rr]]
    stan_dat = list(
      N  = nrow(y),
      T = ncol(y),
      Y = y,
      K = dim(x)[2],
      N_J = max(j), 
      L_J = dim(TT)[2],
      J = j , 
      X =x, 
      TT = TT, 
      C = C, 
      ones  = numeric(max(j)) + 1,
      yz = yz,
      n_max = 10,
      Kp = 1,
      L_Jp = dim(TT_det)[2], 
      TTp = TT_det
    )
    pars <- c("Omega", "tau", "betas", "rho", "z",
              "r", "zp", "taup", "rhop", "Sigma")
    fit <- stan(file = 'RNJSM.stan', 
                data = stan_dat,
                pars = pars, 
                iter = 10000, 
                chains = 3)
    file = paste("fit_JSM_", pp, "_", rr, "F.RData", sep = "")
    save(fit, file = file)
  }
  }
  
  
  # RN basic ----------------------------------------------------------------
  
  rm(list=ls()) 
  require("jagsUI")
  
  # Royle y Nichols 2003 ----------------------------------------------------
  
  cat(file = "RN_model_data.bug","
model {
  # Likelihood
  for (i in 1:M) {
    # Process model
    lambda[i] <- exp(b0 * X[i,1] + b1* X[i,2])
    N[i] ~ dpois(lambda[i])	# latent abundance of species i
    z[i] <- step(N[i] - 1) # z=1 if N>0, ie. site is occupied
    psi_L[i] <- 1 - exp(-lambda[i])
    # Observation model
    p[i] <- 1 - pow(1-r, N[i])
    y[i] ~ dbin(p[i], J[i])
  }
  # Priors
  b0 ~ dnorm(0, 0.01) 
  b1 ~ dnorm(0, 0.01) 
  r ~ dunif(0, 1)
  # Derived quantities
  psi.sample <- mean(z[])   # Proportion of occupied sites in the sample
}
"
  )
  
  ni <- 10000  # iterations
  nt <- 10  # thining
  nb <- 1000  #burin
  nc <- 3  # chains
  nreps = 20
  
  params = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)
  colnames(params)= c("nsites", "n.v")
  
  for(pp in 1:nrow(params))
  {
    file = paste("sim_", pp, "F.RData", sep = "")
    load(file)
    SUM_RN_common = vector("list", nreps)
    SUM_RN_rare = vector("list", nreps)
    rare_spp = rep(NA, nreps)
    com_spp = rep(NA, nreps)
    for(rr in 1:nreps)
    {
      Y = Ys[[rr]]
      X = Xs[[rr]]
      Occ = Occ_m[,rr] 
      rs = P_det[,rr] 
      tmp = which(round(Occ,1) <= 0.1)#Identify and select a rare species
      if(length(tmp)>0)
      {
        if(length(tmp) == 1)
        {
          rare_spp[rr] = tmp
        }else{
          rare_spp[rr] = tmp[which(Occ[tmp] == min(Occ[tmp]))][1]
        }
      }
      tmp = which(round(Occ,1) >= 0.50)#Identify and select a common species
      if(length(tmp)>0)
      {
        if(length(tmp) == 1)
        {
          com_spp[rr] = tmp
        }else{
          com_spp[rr] = tmp[which(Occ[tmp] == max(Occ[tmp]))][1]
        } 
        
      }
      #Adjust rare species
      sp = rare_spp[rr]
      init = ((sp-1)*n.su)+1
      fin = init + (n.su-1)
      Y_sp  = Y[init:fin,]
      y_sp = rowSums(Y_sp)
      visits = rep(ncol(Y_sp), nrow(Y_sp))
      x  = X[init:fin,]
      # Royle-Nichols -----------------------------------------------------------
      data <- list (y = y_sp, J = visits, M = n.su, X = x)
      Nst <- as.numeric(y_sp > 0)  # N must be >0 if y > 0, so start at 1
      inits <- function() list(N = Nst, 
                               b0 = rnorm(1), b1 = rnorm(1),
                               r=runif(1, 0, 1))  
      parameters <- c("b0", "b1",
                      "lambda", "N", "psi_L", "p", "r", "psi.sample")
      mod_RN_rare<- jags(data, inits, parameters, "RN_model_data.bug", n.chains = nc, n.thin = nt, 
                         n.iter = ni, n.burnin = nb)#, parallel = T)
      SUM_RN_rare[[rr]] = as.data.frame(mod_RN_rare$summary)
      
      #Adjust common species
      sp = com_spp[rr]
      init = ((sp-1)*n.su)+1
      fin = init + (n.su-1)
      Y_sp  = Y[init:fin,]
      y_sp = rowSums(Y_sp)
      visits = rep(ncol(Y_sp), nrow(Y_sp))
      x  = X[init:fin,]
      # Royle-Nichols -----------------------------------------------------------
      data <- list (y = y_sp, J = visits, M = n.su, X = x)
      Nst <- as.numeric(y_sp > 0)  # N must be >0 if y > 0, so start at 1
      inits <- function() list(N = Nst, 
                               b0 = rnorm(1), b1 = rnorm(1),
                               r=runif(1, 0, 1))  
      parameters <- c("b0", "b1",
                      "lambda", "N", "psi_L", "p", "r", "psi.sample")
      mod_RN_common<- jags(data, inits, parameters, "RN_model_data.bug", n.chains = nc, n.thin = nt, 
                           n.iter = ni, n.burnin = nb)#, parallel = T)
      SUM_RN_common[[rr]] = as.data.frame(mod_RN_common$summary)
    }
    
    file = paste("fit_RN_rare_", pp, "F.RData", sep = "")
    save(SUM_RN_rare, rare_spp, file = file)
    file = paste("fit_RN_common_", pp, "F.RData", sep = "")
    save(SUM_RN_common, com_spp, file = file)
    
  }#End parameters loop

  
############################################################################  
  # Analyze parameterizations -----------------------------------------------
  rm(list=ls()) 
  library(ape)
  library(mvtnorm)
  library(truncnorm)
  set.seed(1234)
  library("rstan")
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  rstan_options(auto_write = TRUE)
 # load("id_targets.RData")
  
  params = matrix(c(150,3, 75,6,150,6,75,3), ncol = 2, nrow = 4, byrow=T)
  colnames(params)= c("nsites", "n.v")
  
  
  n.s <- 20 # number of bird spp
  n.p <- 2   # number of parameters that model abundance (b0, b1)
  n.t <- 1# number of traits
  #n.t.p <- 1 #number of traits prob. of detection
  #n.p.p <- 1 #number of parameteres that model prob of detection
  nreps = 20
  
  # Detect which simulations did not converge (avoid) -----------------------
  
  CONV = vector("list", nrow(params))
  
  for(p in 1:nrow(params))
  {
    file = paste("fit_RN_common_", p, "F.RData", sep = "")
    load(file)
    file = paste("fit_RN_rare_", p, "F.RData", sep = "")
    load(file)
    conv = array(0, c(nreps, 2))#com?n y rara
    for(r in 1:nreps)
    {
      id_com = com_spp[r]
      RN_com = SUM_RN_common[[r]]
      tmp = c("b0", "b1", "r")
      id1 = which(rownames(RN_com) %in% tmp)
      tmp = c("mean", "2.5%", "97.5%", "Rhat")
      id2 = which(colnames(RN_com) %in% tmp)
      RN_com= RN_com[id1,id2]
      if(length(which(RN_com$Rhat > 1.2)) == 0)
      {conv[r, 1] = 1}
      id_rare = rare_spp[r]
      RN_rare = SUM_RN_rare[[r]]
      tmp = c("b0", "b1", "r")
      id1 = which(rownames(RN_rare) %in% tmp)
      tmp = c("mean", "2.5%", "97.5%", "Rhat")
      id2 = which(colnames(RN_rare) %in% tmp)
      RN_rare= RN_rare[id1,id2]
      if(length(which(RN_rare$Rhat > 1.2)) == 0)
      {conv[r, 2] = 1}
    }
    CONV[[p]] = conv
  }
  
  
  # RN basic -----------------------------------------------------------

    
  OUT = NULL
  for(p in 1:nrow(params))
  {
    file = paste("sim_", p, "F.RData", sep = "")
    load(file)
    file = paste("fit_RN_common_", p, "F.RData", sep = "")
    load(file)
    file = paste("fit_RN_rare_", p, "F.RData", sep = "")
    load(file)
    out_p =NULL
    conv = CONV[[p]]
    #  tmp = params[p,]
    #  id_par = c(tmp$n.v, tmp$r)
    for(r in 1:nreps)
    {
        id_com = com_spp[r]#id common spp
        RN_com = SUM_RN_common[[r]]#single-spp RN parameterization
        tmp = c("b0", "b1", "r")
        id1 = which(rownames(RN_com) %in% tmp)
        tmp = c("mean", "2.5%", "97.5%")
        id2 = which(colnames(RN_com) %in% tmp)
        RN_com= RN_com[id1,id2]
        t_com = Thetas[[r]][id_com,]#true species environmental responses
        p.d = P_det[id_com,r]#spp detectability
        out1 = cbind(rep(p, 3), rep(params[p,2],3), rep(r, 3), rep("com", 3), rep(id_com,3), c("b0", "b1", "r"),
                     c(t_com, p.d), RN_com, rep(conv[r,1],3))
        colnames(out1) = c("pset", "nv", "rep", "type_spp", "id_spp", "param", "true", "m_RN", "L_RN", "U_RN", "Conv")
      #rare
        id_rare = rare_spp[r]#identify rare spp
        RN_rare = SUM_RN_rare[[r]]#single spp RN for rare
        tmp = c("b0", "b1", "r")
        id1 = which(rownames(RN_rare) %in% tmp)
        tmp = c("mean", "2.5%", "97.5%")
        id2 = which(colnames(RN_rare) %in% tmp)
        RN_rare= RN_rare[id1,id2]
        t_rare = Thetas[[r]][id_rare,]#environmental response of rare
        p.d = P_det[id_rare, r]
        out2 = cbind(rep(p, 3), rep(params[p,2],3), rep(r, 3), rep("rare", 3), rep(id_rare,3), c("b0", "b1", "r"),
                     c(t_rare, p.d), RN_rare, rep(conv[r,2],3))
        colnames(out2) = c("pset", "nv", "rep", "type_spp", "id_spp", "param", "true", "m_RN", "L_RN", "U_RN", "Conv")
      
      out_p = rbind(out_p, out1, out2)
    }
    OUT = rbind(OUT, out_p)
  }
  
  write.csv(OUT, "fitted_params.csv", row.names = F)#true and fitted parameters for single-spp RN
  
  
  # Add JSM results -------------------------------------------------
 
   XXX = NULL
  for(p in 1:nrow(params)) #nrow(params))
  {
    OUT_p = NULL
    for(r in 1:nreps)
    {
      file = paste("fit_RN_common_", p, "F.RData", sep = "")
      load(file)
      file = paste("fit_RN_rare_", p, "F.RData", sep = "")
      load(file)
      id_com = com_spp[r]
      id_rare = rare_spp[r]
      RN_com = SUM_RN_common[[r]]
      file = paste("fit_JSM_", p, "_", r, "F.RData", sep = "")
      load(file)
      fit_summary <- summary(fit,  probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
      bs <- fit_summary[grepl("betas", rownames(fit_summary)),]
      pp = 1
      init = ((pp-1)*n.s)+1#identify rows corresponding to the intercept
      fin = init+(n.s-1)
      b0 = data.frame(bs[init:fin,])
      tmp = c("mean", "2.5%", "97.5%")
      id2 = which(colnames(fit_summary) %in% tmp)
      b0_com = b0[id_com,id2]#select fitted values for common spp
      b0_rare = b0[id_rare, id2]#select fitted values for rare species
      pp = 2
      init = ((pp-1)*n.s)+1#
      fin = init+(n.s-1)
      b1 = data.frame(bs[init:fin,])#select b1k parameters
      b1_com = b1[id_com,id2]#select common
      b1_rare = b1[id_rare, id2]#select rare
      pdet <- fit_summary[which(rownames(fit_summary) == "r"),]
      pdet = pdet[id2]#
      out = rbind(b0_com, b1_com, pdet, b0_rare, b1_rare, pdet)
      tmp = nrow(out)
      id_type = c(rep("com", 3), rep("rare",3))
      id_spp = c(rep(id_com, 3), rep(id_rare,3))
      id = cbind(rep(p, tmp), rep(params[p,2], tmp), rep(r, tmp),
                 id_type, id_spp, rep(c("b0", "b1", "r"),2))
      out2 = cbind(id, out)
      OUT_p = rbind(OUT_p, out2)
    }
    XXX = rbind(XXX, OUT_p)
  }
  View(XXX)
  jsm = XXX
  colnames(jsm) =  c("pset", "nv", "rep", "type_spp", "id_spp", "param", "m_JSM", "L_JSM", "U_JSM")
  

# Join RN and JSM datasets ------------------------------------------------

  #Ensure that they have the same order
  fitted = read.csv("fitted_params.csv")
  XXX = NULL
  for(i in 1:nrow(fitted))
  {
    tmp2 = fitted[i,1:6]
    tmp3 = jsm[i, 1:6]
    out = rbind(tmp2, tmp3)
    out = cbind(rep(i, 2), out)
    XXX = rbind(XXX, out)
  }
  View(XXX)
  #Ok
  tmp = c( "m_JSM"  ,  "L_JSM"  ,  "U_JSM")
  tmp2 = which(colnames(jsm) %in% tmp)
  jsm2 = jsm[,tmp2]
  fitted_tot = cbind(fitted, jsm2)#join single spp RN (fitted) with RN-JSM results
  write.csv(fitted_tot, "fitted_params.csv", row.names = F)
  
  
 