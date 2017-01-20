## --------------------------------------------------------------##
## --This program was written by Wenbao Yu.----------------------##
## --note:-----x--maker; y--disease status; R-- missing indicator;
## --------cova1 -- the covariates in disease model
## --------cova2 -- the covariates in verification model
##---------------------------------------------------------------##


library(nleqslv)      # solving nonlinear equations 
library(compiler)
library(pROC)    # calculate ROC indices
library(Rcpp)
library(MASS)
library(maxLik)

sourceCpp('estAUC_nonIg_Rcpp.cpp')  ## call some C codes 


data_simu <- function(n = 200, scenario){
  
  x = runif(n, -1, 1)
  
  if(scenario == "scenarioI"){
    work.v = v = rbinom(n, 1, 0.5)
    tmu = c(2, -2.5, -1)   ## true mu
    tphi = c(1.2, -1, 0)   ## true phi
    beta <- -1.5
  }
  if(scenario == "scenarioII"){
    ## ignorable case
    work.v = v = rbinom(n, 1, 0.5)
    tmu = c(2, -2.5, -1)
    tphi = c(2, -1, -1)
    beta <- 0
  }
  if(scenario == 'scenarioIII'){
    work.v = v = rnorm(n)
    tmu = c(2, -2.5, -1)
    tphi = c(1, -1, 0)
    beta <- -1.5
    work.v = sign(v) * abs(v)^(1/3)
  }
  if(scenario == 'scenarioIV'){
    work.v = v = rnorm(n)
    #work.v = v = rbinom(n, 1, 0.5)
    tmu = c(0.5, -2.5, -1)   ## true mu
    tphi = c(2, -1, 0.8)   ## true phi
    beta <- -2
  }
  if(scenario == 'scenarioV'){
    v = rnorm(n)
    tmu = c(0.5, -2.5, -1)
    tphi = c(2, -1, 0.8)
    beta = -2
    work.v = sign(v) * abs(v)^(1/3)
  } 
  if(grepl(scenario, pattern = 'scenarioVI')){
    v1 = rbinom(n, 1, 0.5)
    v2 = rnorm(n)
    work.v = v = cbind(v1, v2)
    tmu = c(1, -1.5, 0.5, -0.5)
    tphi = c(2, -1, 0, 0)
    beta = -1.5
  }
  
  if(scenario == 'scenarioIII-DP'){
    work.v = v = rnorm(n)
    tmu = c(1, -2.5, -1)
    tphi = c(1.5, -1, 0)
    beta <- -2
    work.v = sign(v) * abs(v)^(1/3)
  }
  
  if(scenario == 'scenarioIII-DV'){
    work.v = v = rnorm(n)
    tmu = c(2, -2.5, -1)
    tphi = c(1.8, -1, 0)
    beta <- -2
    work.v = sign(v) * abs(v)^(1/3)
  }
  
  if(scenario == 'Multiple'){
    ## try using more covariates
    v1 = rbinom(n, 1, 0.5)
    v2 = rnorm(n)
    v3 = runif(n, 0, 1)
    work.v = v = cbind(v1, v2, v3)
    tmu = c(0.6, -1.5, 0.5, -0.5, 0.5)
    tphi = c(1, -1, 0.5, -0.5, 0.5)
    beta = -1.5
  }
  
  cova = cbind(1, x, v)
  py = 1/(1 + exp(cova %*% tmu))  ## true disease model
  y = rbinom(n, 1, py)
  
  ps = 1/(1 + exp(cova %*% tphi + beta * y))
  R = rbinom(n, 1, ps)
  
  v = as.matrix(work.v)
  #true_phi <<- tphi
  return(list('y' = y, 'x' = x, 'v' = v, 'R' = R))
}
data_simu = cmpfun(data_simu)


## estimate parameters mu in disease model
## using partial likelihood
mu_est <- function(yr, R, cova){
  # the disease model p1 = pr(y=1|R=1, 1, v)
  # p1 = 1/(1 + exp(mu[1] + mu[2] * x + mu[3] * v))
  
  # use glm, logistic regression directly
  ind = which(R == 1)
  cova.obs = cova[ind, ]
  yobs = yr[ind]
  
  result = glm(yobs ~ cova.obs , family="binomial") 
  mu <- -coef(result)  
  
  mu[is.na(mu)] = 0  ## avoid constant covariate
  return(mu)
}
mu_est = cmpfun(mu_est)


## estimate parameters phi in verification model 

# using nleqslv package, faster
phi_EM_nleqslv <- function(init_phi, mu, yr, R, cova1, cova2, em_run = 100){
  n = length(R)
  
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  
  old_phi = init_phi
  
  # objective score function
  sbarfun <- function(z){
    sfun <- function(ty){
      s = matrix(0, n, k2 + 2)
      temps = temp0 + z[k2 + 2] * ty
      
      tps = 1/(1 + exp(temps))
      tt = tps - R
      s[, 1] = tt
      s[, k2+2] = ty * tt
      
      for(i in 1:k2){
        s[, i + 1] =  cova2[, i] * tt
      }
      
      return(s)
    }
    temp0 = z[1]
    for(i in 1:k2) temp0 =  temp0 + z[1 + i] * cova2[, i]
    
    s1 = sfun(1)
    s0 = sfun(0)
    E0s = p0 * s1 + (1 - p0) * s0
    sbar = colSums(R * sfun(yr)) + colSums((1 - R) * E0s)
    return(sbar/n)
  }
  
  # jacobian for the objective function
  jac_fun <- function(z){
    DS2 = matrix(0, k2+2, k2+2)
    b0 = c(rep(0, k2 + 1), 1)
    pi_hat = rep(0, n)
    for(i in 1:n){
      ay = c(1, cova2[i, ], yr[i])
      a1 = c(1, cova2[i, ], 1)
      a0 = c(1, cova2[i, ], 0)
      pi1 = 1/(1 + exp(sum(z * a1)))      ## p(R=1|y=1)
      pi0 = 1/(1 + exp(sum(z * a0)))      ## p(R=1|y=0)
      
      pi_hat[i] = 1/(1 + exp(sum(z * ay)))
      ds1 = -pi1 * (1 - pi1) * outer(a1, a1)
      ds0 = -pi0 * (1 - pi0) * outer(a0, a0)
      DS2 = DS2 - R[i] * pi_hat[i] * (1 - pi_hat[i]) * outer(ay, ay) + 
        (1 - R[i]) * (p0[i] * (1 - p0[i]) * outer((pi1 - pi0) * a0 + (pi1 - R[i]) * b0, b0 ) +
                        p0[i] * (ds1 - ds0) + ds0 ) 
    }
    return(DS2)
  }
  
  temp = mu[1]
  for(i in 1:k1){
    temp = temp + mu[1+i] * cova1[, i]
  }
  p1 = 1/(1 + exp(temp))   ## p1=pr(y=1|R=1, x, v)
  
  # using EM
  run = 0
  repeat{
    run = run + 1
    
    p0 = p1 * exp(old_phi[k2+2])/(1 - p1 + p1 * exp(old_phi[k2+2]))
    
    result = nleqslv(old_phi, fn = sbarfun) 
    phi <- result$x
    if(max(abs(old_phi - phi)) < 0.0001 || run == em_run) break
    
    old_phi = phi
    
  }
  cat(run,'\n')
  
  return(list('phi' = phi, 'fval' = sum(abs(result$fvec)), 'convN' = run))
}
phi_EM_nleqslv = cmpfun(phi_EM_nleqslv)


## solve the equivalent optimization problem by EM
phi_EM_optim <- function(init_phi, mu, yr, R, cova1, cova2, em_run = 100){
  n = length(R)
  
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  
  old_phi = init_phi
  
  # objective score function
  opt_fun <- function(z){
    temp0 = z[1]
    for(i in 1:k2) temp0 =  temp0 + z[1 + i] * cova2[, i]
    pi = 1/(1 + exp(temp0 + z[k2 + 2] * yr))
    pi_1 = 1/(1 + exp(temp0 + z[k2 + 2]))
    pi_0 = 1/(1 + exp(temp0))
    
    E0s = p0 * log(1 - pi_1) + (1 - p0) * log(1 - pi_0)
    sbar = sum(R * log(pi)) + sum((1 - R) * E0s)
    return(sbar/n)
  }
  
  temp = mu[1]
  for(i in 1:k1){
    temp = temp + mu[1+i] * cova1[, i]
  }
  p1 = 1/(1 + exp(temp))   ## p1=pr(y=1|R=1, x, v)
  
  # using EM
  run = 0
  repeat{
    run = run + 1
    
    p0 = p1 * exp(old_phi[k2+2])/(1 - p1 + p1 * exp(old_phi[k2+2]))
    
    result = maxLik(start = old_phi, logLik = opt_fun) 
    phi <- result$estimate
    if(max(abs(old_phi - phi)) < 0.001 || run == em_run) break
    
    old_phi = phi
    
  }
  cat(run,'\n')
  
  return(list('phi' = phi, 'convN' = run))
}
phi_EM_optim = cmpfun(phi_EM_optim)


## parameter estimation for liu and zhou's (2010) method, note asAUC_fp in our paper
## output variance of parameters too, using optim function
par_est_lzh <- function(init.par, yr, R, cova1, cova2){
  ## cova1 -- the covariates in disease model
  ## cova2 -- the covariates in verification model
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  obj <- function(z){
    temp1 = z[1] 
    temp2 = z[k1+2] 
    for(i in 1:k1) temp1 = temp1 + z[i+1] * cova1[, i]
    for(i in 1:k2) temp2 = temp2 + z[k1+2 + i] * cova2[, i]
    
    rho = 1 - 1/(1 + exp(temp1))
    pi1 = 1 - 1/(1 + exp(temp2 + z[k1 + k2 + 3]))
    pi0 = 1 - 1/(1 + exp(temp2))
    lvalues = sum(yr * R * log(pi1 * rho)) + sum((1 - yr) * R * log(pi0 * (1 - rho))) +
      sum((1 - R) * log(1 - pi1 * rho - pi0 * (1 - rho)))
    return(-lvalues)
  }
  
  opt_res = optim(par = init.par, fn = obj, hessian = TRUE)
  pars = opt_res$par

  hessian = -opt_res$hessian

  vars = diag(-ginv(hessian))
  
  temp = pars[k1+2] + pars[k1 + k2 + 3] * yr
  for(i in 1:k2) temp = temp + pars[k1+2+i] * cova2[, i]
  
  ps_lzh =  1 - 1/(1 + exp(temp))
  mu = pars[1:(k1+1)]
  phi = pars[-c(1:(k1+1))]
  
  n = length(R)
  vars_mu = vars[1: (k1+1)]
  vars_phi = vars[-(1:(k1+1))]
  
  pv_mu = 2 * pnorm(-abs(mu)/sqrt(vars_mu))
  pv_phi = 2 * pnorm(-abs(phi)/sqrt(vars_phi))
  
  return(list('ps' = ps_lzh, 'phi' = -phi, 'mu' = -mu,
              'vars_mu' = vars_mu, 'vars_phi' = vars_phi,
              'pv_mu' = pv_mu, 'pv_phi' = pv_phi, 'hessian' = hessian))
}
par_est_lzh = cmpfun(par_est_lzh)


## alternatively, using maxLik package
par_est_lzh_mle <- function(init.par, yr, R, cova1, cova2){
  ## cova1 -- the covariates in disease model
  ## cova2 -- the covariates in verification model
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  obj <- function(z){
    temp1 = z[1] 
    temp2 = z[k1+2] 
    for(i in 1:k1) temp1 = temp1 + z[i+1] * cova1[, i]
    for(i in 1:k2) temp2 = temp2 + z[k1+2 + i] * cova2[, i]
    
    rho = 1 - 1/(1 + exp(temp1))
    pi1 = 1 - 1/(1 + exp(temp2 + z[k1 + k2 + 3]))
    pi0 = 1 - 1/(1 + exp(temp2))
    lvalues = sum(yr * R * log(pi1 * rho)) + sum((1 - yr) * R * log(pi0 * (1 - rho))) +
      sum((1 - R) * log(1 - pi1 * rho - pi0 * (1 - rho)))
    return(lvalues)
  }
  
  opt_res = maxLik(start = init.par,  logLik = obj)
  pars = opt_res$estimate

  hessian = opt_res$hessian
  #grad = opt_res$gradient

  vars = diag(-ginv(hessian))
  
  temp = pars[k1+2] + pars[k1 + k2 + 3] * yr
  for(i in 1:k2) temp = temp + pars[k1+2+i] * cova2[, i]
  
  ps_lzh =  1 - 1/(1 + exp(temp))
  mu = pars[1:(k1+1)]
  phi = pars[-c(1:(k1+1))]
  
  n = length(R)
  vars_mu = vars[1: (k1+1)]
  vars_phi = vars[-(1:(k1+1))]
  
  pv_mu = 2 * pnorm(-abs(mu)/sqrt(vars_mu))
  pv_phi = 2 * pnorm(-abs(phi)/sqrt(vars_phi))
  

  # make phi be consistent as in our notation
  return(list('ps' = ps_lzh, 'phi' = -phi, 'mu' = -mu,
              'vars_mu' = vars_mu, 'vars_phi' = vars_phi,
              'pv_mu' = pv_mu, 'pv_phi' = pv_phi, 'hessian' = hessian))
}
par_est_lzh_mle = cmpfun(par_est_lzh_mle)


## estimate variance of AUC_fp
var_auc_lzh <- function(pi_hat, mu, phi, yr, R, cova1, cova2, xid = 1, auc1, hessian){ 
  n = length(R)
  
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  

  
  ## calculate the first term of variance, using c
  x = cova1[, xid]
  F1 = cal_F1(pi_hat, auc1, yr, x, R)
    
  # calculate Gamma1 (sum of gradient of Dij), using c
  gamma1 = cal_gamma1_v2(pi_hat, auc1, yr, R, cova2)  ## length(phi) vector
  gamma1 = c(rep(0, length(mu)), gamma1)

  ## gradient of the likelihood function
  grad = matrix(0, n, k1 + k2 + 3) 
  for(i in 1:n){
    ax = c(1, cova1[i, ])
    rho = 1/(1 + exp(sum(mu * ax)))
    a1 = c(1, cova2[i, ], 1)
    a0 = c(1, cova2[i, ], 0)
    pi1 = 1/(1 + exp(sum(phi * a1)))      ## p(R=1|y=1)
    pi0 = 1/(1 + exp(sum(phi * a0)))      ## p(R=1|y=0)
    denom = (1 - pi0 + rho * (pi0 - pi1))
    grad[i, 1:(k1 + 1)] = (R[i] * yr[i] * (rho - 1) + R[i] * (1 - yr[i]) * rho + 
                           (1 - R[i]) * rho * (rho -1) * (pi0 - pi1) /denom) * ax
    grad[i, -(1:(k1 + 1))] = R[i] * (yr[i] * (pi1 - 1) * a1 + (1 -yr[i]) * (pi0 - 1) * a0) + 
                             (1- R[i]) * ((1 - rho)*pi0*(1-pi0)*a0 + rho*pi1*(1-pi1) * a1)/denom
  }

  vq = var(F1 - grad %*% ginv(hessian) %*% gamma1/n)
 
  dp1 = sum(R * yr / pi_hat)/n    ## estimation of pr(y=1)
  dp0 = sum(R * (1 - yr)/ pi_hat)/n  ## another way of estimating pr(y=0)
  #dp0 = 1- dp1

  vA = vq/(dp1 * dp0) ^2
  
    
  return(vA/n)
}
var_auc_lzh = cmpfun(var_auc_lzh)


## estimate variance of AUC_IV, the proposed one
## make use of c function
var_auc_iv <- function(pi_hat, mu, phi, yr, R, cova1, cova2, xid = 1, auc1){ 
  n = length(R)
  
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  
  temp = mu[1]
  
  for(i in 1:k1) temp = temp + mu[i+1] * cova1[, i]
  
  
  p1 = 1/(1 + exp(temp))    ## pr(y=1|R=1)
  
  ebeta = exp(phi[k2 + 2])
  p0 = p1 * exp(phi[k2 + 2])/(1 - p1 + p1 * ebeta)  ## pr(y=1|R=0) 
  
  ## calculate the first term of variance, using c
  x = cova1[, xid]
  F1 = cal_F1(pi_hat, auc1, yr, x, R)
    
  # calculate Gamma1 and Gamma2, using c
  
  gamma1 = cal_gamma1_v2(pi_hat, auc1, yr, R, cova2)
   
  # calculate s2 and derivative of s2
  s2 = matrix(0, n, k2 + 2)
  s1 = matrix(0, n, k1 + 1)
  FIM2 = matrix(0, k2+2, k2+2)
  FIM1 = matrix(0, k1+1, k1+1)
  b0 = c(rep(0, k2 + 1), 1)
  ds2_mu = matrix(0, k1+1, k2+2)
  for(i in 1:n){
    ay = c(1, cova2[i, ], yr[i])
    a1 = c(1, cova2[i, ], 1)
    a0 = c(1, cova2[i, ], 0)
    pi1 = 1/(1 + exp(sum(phi * a1)))      ## p(R=1|y=1)
    pi0 = 1/(1 + exp(sum(phi * a0)))      ## p(R=1|y=0)
    s2[i, ] = R[i] * (pi_hat[i] - R[i]) * ay + (1 - R[i]) * 
      (p0[i] * (pi1 - R[i]) * a1 + (1 - p0[i]) * (pi0 - R[i]) * a0)
    
    
    ds2_1 = -pi1 * (1 - pi1) * outer(a1, a1)
    ds2_0 = -pi0 * (1 - pi0) * outer(a0, a0)
    FIM2 = FIM2 - R[i] * pi_hat[i] * (1 - pi_hat[i]) * outer(ay, ay) + 
      (1 - R[i]) * (p0[i] * (1 - p0[i]) * outer((pi1 - pi0) * a0 + (pi1 - R[i]) * b0, b0 ) +
                      p0[i] * (ds2_1 - ds2_0) + ds2_0 ) 

    ## add variance term comes from estimation mu (added since 10/23/2016)
    ax = c(1, cova1[i, ])
    s1[i, ] = R[i] * (p1[i] - yr[i]) * ax 
    FIM1 = FIM1 - R[i] * p1[i] * (1 - p1[i]) * outer(ax, ax)
    ds2_mu = ds2_mu - (1 - R[i]) * p1[i] * (1 - p1[i]) * ebeta/(1-p1[i]+p1[i]*ebeta)^2 *
             (pi1 * outer(ax, a1) - pi0 * outer(ax, a0)) 
  }
  ds2_mu = ds2_mu/n
  FIM1 = FIM1/n
  FIM2 = FIM2/n
  gamma1 = gamma1/n/n
  vq = var(F1 - (s2 - s1 %*% ginv(FIM1) %*% ds2_mu) %*% ginv(FIM2) %*% gamma1)
 
  dp1 = sum(R * yr / pi_hat)/n    ## estimation of pr(y=1)
  dp0 = sum(R * (1 - yr)/ pi_hat)/n  ## another way of estimating pr(y=0)
  #dp0 = 1- dp1

  vA = vq/(dp1 * dp0) ^2
  
    
  return(vA/n)
}
var_auc_iv = cmpfun(var_auc_iv)


## estimate variance for AUC_ig
var_auc_ig <- function(pi_hat, phi, yr, R, cova1, xid  = 1, auc1){ 
  n = length(R)
  
  k1 = ncol(cova1)
  
  ## calculate the first term of variance, using c
  x = cova1[, xid]
  F1 = cal_F1(pi_hat, auc1, yr, x, R)
    
  # calculate Gamma1  using c
  gamma1 = cal_gamma1_ig(pi_hat, auc1, yr, R, cova1)
  s1 = matrix(0, n, k1 + 1) 
  FIM1 = matrix(0, k1+1, k1+1)
  for(i in 1:n){
    
   ax = c(1, cova1[i, ])
  
   pi1 = 1/(1 + exp(sum(phi * ax)))      ## p(R=1|y=1) = p(R=1|y=0)
    
   s1[i, ] = (pi1 - R[i]) * ax 
   FIM1 = FIM1 - pi1 * (1 - pi1) * outer(ax, ax)
   
  }
  
  FIM1 = FIM1/n
  
  gamma1 = gamma1/n/n
  vq = var(F1 - s1 %*% ginv(FIM1) %*% gamma1)
 
  dp1 = sum(R * yr / pi_hat)/n    ## estimation of pr(y=1)
  dp0 = sum(R * (1 - yr)/ pi_hat)/n  ## another way of estimating pr(y=0)
  #dp0 = 1- dp1

  vA = vq/(dp1 * dp0) ^2
  
    
  return(vA/n)
}
var_auc_ig = cmpfun(var_auc_ig)


# estimate auc_iv and its variance
est_auc_iv <- function(init_phi, yr, R, cova1, cova2, xid = 1, em_run){  
    k1 = ncol(cova1)
    k2 = ncol(cova2)
    
    mu = mu_est(yr, R, cova1)
    temp_res = phi_EM_nleqslv(init_phi, mu, yr, R, cova1, cova2, em_run)
    #temp_res = phi_EM_optim(init_phi, mu, yr, R, cova1, cova2, em_run)
    
    phi = temp_res$phi
  
    temp = phi[1] + phi[k2 + 1] * yr
    for(j in 1:k2) temp = temp + phi[1 + j] * cova2[, j]
    pi_hat = 1/(1 + exp(temp))
    
    x = cova1[, xid]
  
    ## calculate Auc_iv
    auc_iv = cal_auc1(pi_hat, yr, x, R)
    var_auc = var_auc_iv(pi_hat, mu, phi, yr, R, cova1, cova2, xid, auc_iv)
    
    return(c(auc_iv, var_auc))
}


est_auc_lzh <- function(inits, yr, R, cova1, cova2, xid = 1){
 # inits = rep(0,  ncol(cova1) + ncol(cova1) + 3)
 # inits[1] = -1
  x = cova1[, xid]
  res = par_est_lzh(inits, yr, R, cova1, cova1)
  hessian = res$hessian
  ps_lzh = res$ps
  mu = res$mu
  phi = res$phi

  auc_fp= cal_auc1(ps_lzh, yr, x, R)
 
  var_auc = var_auc_lzh(ps_lzh, mu, phi, yr, R, cova1, cova2, xid = 1, auc_fp, hessian)
  return(c(auc_fp, var_auc))
}


est_auc_ig <- function(yr, R, cova1, xid = 1){
  result = glm(R ~ cova1 , family = "binomial")
  phi_ig = - coef(result)
  temp = phi_ig[1] 
  for(i in 1:ncol(cova1)) temp = temp + phi_ig[1 + i] * cova1[, i]
  ps_ig = 1/(1 + exp(temp))
  x = cova1[, xid]
  auc_ig = cal_auc1(ps_ig, yr, x, R)
  var_auc = var_auc_ig(ps_ig, phi_ig, yr, R, cova1, xid, auc_ig)
  return(c(auc_ig, var_auc))
}


est_auc_v <- function(yr, R, x){
  ## auc using verified samples only
  obs.id = which(R == 1)
  res = roc(response = yr[obs.id], predictor = x[obs.id])
  auc_v = res$auc
  var_auc = var(res)
  return(c(auc_v, var_auc))
}


est_auc_f <- function(y, x){
  ## auc all samples suppose no missing
  res = roc(response = y, predictor = x)
  auc_f = res$auc
  var_auc = var(res)
  return(c(auc_f, var_auc))
}


## estimate Omega -- for getting variance of phi, not used
est_Omega <- function(pi_hat, R, yr, mu, phi, cova1, cova2){
  k = length(phi)
  n = length(R)
  kk = ncol(cova1) + 1
  u = s = s2 = matrix(0, n, k)
  s1 = matrix(0, n, kk)
  T11 = matrix(0, kk, kk)
  T21 = kabba = matrix(0, k, kk)
  T22 = matrix(0, k, k)
  
  #tempO = phi[1]
  #for(i in 1:ncol(cova2)) tempO = tempO + phi[i+1] * cova2[, i]
  #O = as.vector((1 + exp(tempO + phi[k]))/
  #                (1 + exp(tempO)))
  O = 1
  ## p1=pr(y=1|R=1, cova1), p0=pr(y=1|R=0, cova1)
  temp = mu[1] 
  for(i in 1:ncol(cova1)) temp = temp + mu[i+1] * cova1[, i]
  p1 = 1/(1 + O * exp(temp))
  
  p0 = p1 * exp(phi[k])/(1 - p1 + p1 * exp(phi[k]))
  
  s1[, 1] =  R * (yr - p1)
  for(i in 1:ncol(cova1)) s1[, 1+i] = s1[, 1] * cova1[, i]
  
  s2 = matrix(0, n, k)
  
  b0 = c(rep(0, k-1), 1)
  for(i in 1:n){
    ay = c(1,  cova2[i, ], yr[i])
    av = c(1,  cova1[i, ])
    a1 = c(1,  cova2[i, ], 1)
    a0 = c(1,  cova2[i, ], 0)
    pi1 = 1/(1 + exp(sum(phi * a1)))      ## pi(R=1|y=1)
    pi0 = 1/(1 + exp(sum(phi * a0)))      ## pi(R=1|y=0)
    s2[i, ] = R[i] * (pi_hat[i] - R[i]) * ay + (1 - R[i]) * 
      (p0[i] * (pi1 - R[i]) * a1 + (1 - p0[i]) * (pi0 - R[i]) * a0)
    
    T11 = T11 - R[i] * (1 - p1[i]) * p1[i] * outer(av, av)
    temp = (pi1 - R[i]) * a1 - (pi0 - R[i]) * a0
    T21 = T21 + (1 - R[i]) * p0[i]^2 * (1 - p1[i]) * outer(temp, av)
    
    ds1 = -pi1 * (1 - pi1) * outer(a1, a1)
    ds0 = -pi0 * (1 - pi0) * outer(a0, a0)
    T22 = T22 + R[i] * pi_hat[i] * (1 - pi_hat[i]) * outer(ay, ay) - 
      (1 - R[i]) * (p0[i] * (1 - p0[i]) * outer((pi1 - pi0) * a0 + (pi1 - R[i]) * b0, b0 ) +
                      p0[i] * (ds1 - ds0) + ds0) 
  }
  T11 = T11/n
  T21 = T21/n
  T22 = T22/n
  ivT11 = solve(T11)
  kabba = T21 %*% ivT11
  u = s2 - s1 %*% t(kabba) 
  ivT22 = solve(T22)
  Omega = ivT22 %*% cov(u) %*% t(ivT22)
  return(Omega)
}
est_Omega = cmpfun(est_Omega)


