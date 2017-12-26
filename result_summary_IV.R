## make results into latex table
library(tables)

summary_dat <- function(scenario, n){
  filename1 = paste('./simu_result_for_n=', n,  '_scneario=', 
                    scenario,  '.Rdata', sep='')
  load(file = filename1)
  ests = colnames(out_auc)
  
  ## make all the result in a data frame
  dat = data.frame( 'methods' = ests,
                    stringsAsFactors = FALSE)
  
  
  ## calculate bias and sample variance for each scenario
  if(scenario %in% c('scenarioI', 'scenarioII')) aucs = 0.8137
  if(scenario %in% c('scenarioIII', 'scenarioIII-DP', "scenarioIII-DV", 'scenarioIV', 'scenarioV')) aucs = 0.792
  if(scenario %in% c('scenarioVI')) aucs = 0.713
 
  
  aucs = mean(out_auc[, 5])
  bias = apply(out_auc, 2, mean) - aucs
  
  vars = apply(out_auc, 2, var)
  
  
  ## calculate coverage probability
  B = 500
  cfun <- function(ind){
    #se = (sqrt(asym_var))
    se = median((sqrt(asym_var)))
    if(ind == 5 ) se = sd(out_auc[, ind])
    CIc = (ifelse(aucs >= (out_auc[, ind] - 1.96 * se) & aucs <= (out_auc[,ind] + 1.96 * se), 1, 0))
    #mean(CIc[-which(asym_var)>1])
    mean(CIc)
  }
  cp = sapply(1:5, cfun)
  
  dat$bias = bias
  dat$cp = cp
  dat$n = rep(n, 5)
  
  dat$svar = vars/vars[5]
  dat$smse = (bias^2 + vars)/(bias[5]^2 + vars[5])
  
  return(dat)
}

summary_vars <- function(scenario, n){
  filename1 = paste('./simu_result_for_n=', n,  '_scneario=',
                    scenario,  '.Rdata', sep='')
  load(file = filename1)
  
  vars = apply(out_auc, 2, var)
  
  avar = median(asym_var)
  
  var12 = c(n, vars, avar)
  return(var12)
}



scenario = commandArgs(TRUE)   ## specify a scenario

dats = summary_dat(scenario, 200)
datl = summary_dat(scenario, 2000)

outdat = rbind(dats, datl)

## write latex table


dat1 = data.frame('n'= as.factor(outdat$n), 
                  'methods' = outdat$methods, 'bias' = outdat$bias, 'cp' = outdat$cp,
                  'var' = outdat$svar, 'mse' = outdat$smse, stringsAsFactors=T)

dat1


tt = tabular(data = dat1,   methods ~ 
               n *(round(bias, 3) + round(var, 3) + round(mse, 3) + Format(digits = 3) * cp) * mean)
latex(tt)
tt

tt = tabular(data = dat1,   methods ~ 
               n *(round(bias, 3) + round(var, 3) + round(mse, 3)) * mean)
latex(tt)
tt

## extract variance ####
vars = round(rbind(summary_vars(scenario, 200), summary_vars(scenario, 2000)), 4)
colnames(vars)[1] = 'n'
colnames(vars)[7] = 'asym-var'
vars


out_file = 'summary_results_IV.txt'
write(scenario, file = out_file, append=T)

vars = format(vars, digits = 4, justify='left')

dat1 = format(dat1, digits = 4, justify='left')

#write.table(vars, file = out_file, append=T, row.names=F, quote = F)
#write.table(dat1, file = out_file, append=T, row.names=F, quote = F)

dat1
