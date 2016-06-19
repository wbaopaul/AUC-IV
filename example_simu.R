## --------------------------------------------------------------##
## --This program was written by Wenbao Yu.----------------------##
## --It is a revised version based on Biometric resivision ------##
## --not include A2 ---------------------------------------------##
## --note:-----x--maker; y--disease status; R-- missing indicator;
## --------cova1 -- the covariates in disease model
## --------cova2 -- the covariates in verification model
##---------------------------------------------------------------##


source("auc_IV_source.R")


## run simulation B times ####


rept_run_simu <- function(n = 200,  B=500, scenario, 
                          cova.id1 = 1, cova.id2 = 1, em_run){
  ## cova.id1 -- the covaraites ids used in disease model
  ## cova.id2 -- the covariates ids used in verification model
  
  out_auc = matrix(0, B, 5)
  asym_var = rep(0, B)
  set.seed(1)
  i = 0
  repeat{
    i = i + 1
    dat = data_simu(n, scenario)
    ##dat = dat_simu_lzh(n)
    y = dat$y  ## only used for calculate the full estimator
    x = dat$x
    v = as.matrix(dat$v)
    R = dat$R
    yr <- ifelse(R == 1, y, -1)  
    ## to avoid NA, recode NA in the known y by another number
    ## note that when R==0, y is not accurally used in any piece of the following programming, 
    ## so how to code those y does not matter
    
    cova = cbind(x, v)
    cova1 = as.matrix(cova[, cova.id1])
    cova2 = as.matrix(cova[, cova.id2])
    k1 = ncol(cova1)
    k2 = ncol(cova2)
    
    # try different initial parameters for phi
    
    
    init_phi = rep(0, k2 + 2)
    init_phi[1] = 1
    mu = mu_est(yr, R, cova1)
    temp_res = phi_EM_nleqslv(init_phi, mu, yr, R, cova1, cova2, em_run)

    
    phi = temp_res$phi
    
    temp = phi[1] + phi[k2 + 1] * yr
    for(j in 1:k2) temp = temp + phi[1 + j] * cova2[, j]
    pi_hat = 1/(1 + exp(temp))
    
    out_auc[i, 1:4] = auc_est(pi_hat, yr, R, xid=1, cova1)
    
    ## the full estimator (as if all y are observed)
    out_auc[i, 5] = auc(y, x)
    
    asym_var[i] = Cvar_est_alter(pi_hat, mu, phi, yr, R, cova1, cova2,
                                 out_auc[i, 1])
    if(i == B)  break
  }
  
  colnames(out_auc) = c('A-iv', 'A-ig', 'A-v', 'A-fp', 'A-f')
  
  #save results in file
  #filename1 = paste('./', 'result1_for_n=', n,  '_scneario=', scenario,  '.Rdata', sep='')
  #save(out_auc, asym_var, file = filename1)
  
  bias = apply(out_auc, 2, mean) - mean(out_auc[, 5])
  vars = apply(out_auc, 2, var)
  return(rbind(bias, vars))
}
rept_run_simu = cmpfun(rept_run_simu)



## run scenario  ####

# specify the scenario you want to run, for example
scenario = 'scenarioI'

if( !grepl(scenario, pattern = 'scenarioVI') & 
      !grepl(scenario, pattern = 'scenarioM')){
  cova.id1 = 1:2     
  cova.id2 = cova.id1[-2]
}

if( scenario == 'scenarioVI' ){
  cova.id1 = 1:3     
  cova.id2 = cova.id1[-2]
}


if(scenario == 'Multiple'){
  cova.id1 = 1:4     
  cova.id2 = cova.id1[-4]
}


B = 500       ## number of data sets for each scenario
em_run = 500  ## maximum number of repeats for EM
rept_run_simu(200,  B, scenario, cova.id1, cova.id2, em_run)
rept_run_simu(2000,  B, scenario, cova.id1, cova.id2, em_run)
