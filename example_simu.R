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
  asym_var = out_auc
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
    
    tmp_iv = est_auc_iv(init_phi, yr, R, cova1, cova2, xid = 1, em_run)

    inits = rep(0, 2*k1 + 3)
    inits[1] = -1
    tmp_lzh = est_auc_lzh(inits, yr, R, cova1, cova1, xid = 1)
    tmp_ig = est_auc_ig(yr, R, cova1, xid = 1)
    tmp_v = est_auc_v(yr, R, cova1[, 1])
    tmp_f = est_auc_f(y, cova1[, 1])

    out_auc[i, 1] = tmp_iv[1]
    asym_var[i, 1] = tmp_iv[2]

    out_auc[i, 2] = tmp_ig[1]
    asym_var[i, 2] = tmp_ig[2]

    out_auc[i, 3] = tmp_v[1]
    asym_var[i, 3] = tmp_v[2]
    
    out_auc[i, 4] = tmp_lzh[1]
    asym_var[i, 4] = tmp_lzh[2]

    out_auc[i, 5] = tmp_f[1]
    asym_var[i, 5] = tmp_f[2]

    if(i == B)  break
  }
  
  colnames(out_auc) = c('A-iv', 'A-ig', 'A-v', 'A-fp', 'A-f')
  
  
  bias = apply(out_auc, 2, mean) - mean(out_auc[, 5])
  vars = apply(out_auc, 2, var)
  
  ## save results in file
  filename1 = paste('./', 'result1_for_n=', n,  '_scneario=', scenario,  '.Rdata', sep='')
  save(out_auc, asym_var, file = filename1)
  
  ## calculate coverage probability
  
  
  return(rbind(bias, vars))
}
rept_run_simu = cmpfun(rept_run_simu)



## run scenario  ####

# specify the scenario you want to run, for example
args = commandArgs(T)

scenario = args[1]

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

