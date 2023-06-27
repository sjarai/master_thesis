########################################################################################################
### This R-code performs the simulation study belonging to the thesis: Correcting for Linkage Errors ###
### in Contingency Tables: New Methods to Improve the Correction Approach by Sjarai Eikenhout.       ###
### The code is paralleled and takes approximately 3 days to run.                                    ### 
### Comments are added to the code to explain what is happening.                                     ###
### The corresponding formulas in the thesis are also given as comments for each function.           ###
########################################################################################################


###########################
### Start with clear memory
rm(list = ls())
gc()


#########################
### Set working directory
setwd("S:/Code")


#################
### Load packages
library(gtools)
library(dummies)
library(foreach)
library(doParallel)


###################################
### Preparation for parallel coding
cores = detectCores()
cl = makeCluster(cores[1]-1)
registerDoParallel(cl)


########################
### Generating variables

set.seed(12345)

n = 300 # number of records of the data sets that will be linked
n_tab = 10 # number of tables we will be generating

## auxiliary variable x
x = rnorm(n*n_tab, 20, 10)

## target variable z
z = rnorm(n*n_tab, 20, 15)

## target variable y
y = 20 + 2*x + 1*z + rnorm(n *n_tab, 0, 4)

## variable w
w = rnorm(n*n_tab, 40, 30)


#################################
### Set the number of simulations

n_sim = 10000


########################
### Discretize variables

n_cat1 = 5 # number of categories for variable y in data set 1
n_cat2 = 5 # number of categories for variable z & w in data set 2

Ylist = vector(mode = "list", length = n_tab)
Zlist = vector(mode = "list", length = n_tab)
Wlist = vector(mode = "list", length = n_tab)

freq1list = vector(mode = "list", length = n_tab)
freq2list = vector(mode = "list", length = n_tab)

ass_deplist = vector(mode = "list", length = n_tab)
ass_indeplist = vector(mode = "list", length = n_tab)

for (i in 1:n_tab){
  yg = cut(y[(i-1)*n + (1:n)], breaks = n_cat1, labels = 1:n_cat1) # categorize y into n_cat1 categories
  Ylist[[i]] = dummy(yg) # discretize y
  
  zg = cut(z[(i-1)*n + (1:n)], breaks = n_cat2, labels = 1:n_cat2) # categorize z into n_cat2 categories
  Zlist[[i]] = dummy(zg) # discretize z
  
  wg = cut(w[(i-1)*n + (1:n)], breaks = n_cat2, labels = 1:n_cat2) # categorize w into n_cat2 categories
  Wlist[[i]] = dummy(wg) # discretize w
  
  
  #########################################################################
  ## Compute true contingency table & linkage error probabilities matrix Q
  
  freq1list[[i]] = t(Ylist[[i]]) %*% Zlist[[i]] # list with 10 tables t(Y_g)*Z_g with dependent attributes
  freq2list[[i]] = t(Ylist[[i]]) %*% Wlist[[i]] # list with 10 tables t(Y_g)*W_g with independent attributes
  
  
  ##########################
  ## Determine associations
  
  ## For dependent attributes
  R1 = matrix(rep(rowSums(freq1list[[i]]), each = n_cat2), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
  S1 = matrix(rep(colSums(freq1list[[i]]), times = n_cat1), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
  
  ass_dep = ifelse(n*freq1list[[i]] - R1*S1 < 0, 'neg', ifelse(n*freq1list[[i]] - R1*S1 > 0, 'pos', 'equal'))
  ass_deplist[[i]] = as.vector(t(ass_dep)) # column vector containing the associations (over the rows)
  
  ## For independent attributes
  R2 = matrix(rep(rowSums(freq2list[[i]]), each = n_cat2), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
  S2 = matrix(rep(colSums(freq2list[[i]]), times = n_cat1), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
  
  ass_indep = ifelse(n*freq2list[[i]] - R2*S2 < 0, 'neg', ifelse(n*freq2list[[i]] - R2*S2 > 0, 'pos', 'equal'))
  ass_indeplist[[i]] = as.vector(t(ass_indep)) # column vector containing the associations (over the rows)
}


################################################
### set all values for q that will be considered

qvec = seq(0.9, 0.1, by = -0.1)


###################################### 
### Start of parallel simulation study

res_parallel = foreach(iq = 1:length(qvec),
                       .packages = c("gtools")) %dopar% {
                         
                         ################################
                         ## Set seed for reproducibility
                         
                         set.seed(iq)  
                         
                         
                         ###############################
                         ### Set q en compute Q and Qinv
                         
                         q = qvec[iq]
                         # compute error matrix Q (see formula (4))
                         Q = matrix(rep((1-q)/(n-1), n^2), nrow = n, ncol = n)
                         diag(Q) = q
                         # compute the inverse of error matrix Q (see formula (7))
                         Qinv = (1/(n*q-1)) * (((n-1)*diag(n)) - ((1-q)*c(rep(1, n))%*%t(c(rep(1, n)))))
                         
                         
                         ###################
                         ### Load functions
                         setwd("S:/Code") # set to the right directory that contains the file Functions.R
                         source("Functions.R")
                         
                         
                         ##################
                         ### Compute lambda
                         
                         lambda = lambdaExact(n = n, q = q)
                         
                         
                         #####################################
                         ### Generating permutation matrices C
                         
                         Clist = list(n_sim)
                         for (i in 1:n_sim) {
                           Clist[[i]] = drawC(lambda = lambda, n = n)
                         }
                         
                         
                         ######################################################
                         ### Generating results for each table & each estimator
                         
                         dev_list1 = vector(mode = "list", length = n_tab) # to store results for the dependent tables
                         dev_list2 = vector(mode = "list", length = n_tab) # to store results for the independent tables
                         
                         for (i in 1:n_tab){
                           freq1 = freq1list[[i]]
                           R1 = matrix(rep(rowSums(freq1), each = n_cat2), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
                           S1 = matrix(rep(colSums(freq1), times = n_cat1), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
                           
                           freq2 = freq2list[[i]]
                           R2 = matrix(rep(rowSums(freq2), each = n_cat2), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
                           S2 = matrix(rep(colSums(freq2), times = n_cat1), nrow = n_cat1, ncol = n_cat2, byrow = TRUE)
                           
                           Yg = Ylist[[i]]
                           Zg = Zlist[[i]]
                           Wg = Wlist[[i]]
                           
                           
                           #########################################
                           ### Computations for dependent attributes
                           
                           v.prag = v.prag_fun(n = n, q = q) # (see formula (17))
                           v.opt.true = v.opt_fun(n = n, q = q, freq = freq1, R = R1, S = S1) # (see formula (16))
                           
                           u.prag = u.prag_fun(n = n, q = q) # (see formula (A.34)) !not included in thesis!
                           u.opt.true = u.opt_fun(n = n, q = q, freq = freq1, R = R1, S = S1) # (see formula (A.33)) !not included in thesis!
                           
                           dev = sapply(Clist, function(C) {
                             
                             ## Simulate linkage errors
                             Zstar = C %*% Zg 
                             
                             ## Existing estimators considered by Scholtus, Shlomo and De Waal (2022)
                             freq.obs = t(Yg) %*% Zstar # (see formula (8))
                             freq.Q = t(Q %*% Yg) %*% Zstar # == reg.est(0) (see formula (9))
                             freq.Qinv = t(Yg) %*% (Qinv %*% Zstar) # == reg.est(1) (see formula (10))
                             
                             # Set negative values to 0 for Qinv
                             freq.Qinv = ifelse(freq.Qinv < 0, 0, freq.Qinv)
                             
                             ## Expected value by means of probabilities using alternative prior where all possible values have the same probability
                             probabilities1 = rev_probabilities_fun(n, lambda, freq.obs, R1, S1)$alternative1 # (see formula (27) in combination with (18))
                             freq.exp1 = expected_values_fun(n_cat1, n_cat2, probabilities1) # (see formula (28))
                             
                             ## Expected value by means of probabilities using proposed prior by Scholtus & de Waal (2020-2022)
                             probabilities2 = rev_probabilities_fun(n, lambda, freq.obs, R1, S1)$proposed # (see formula (27) in combination with (19))
                             freq.exp2 = expected_values_fun(n_cat1, n_cat2, probabilities2) # (see formula (28))
                             
                             ## Expected value by means of probabilities using alternative prior by de Waal (2023) where a truncated binomial model is used
                             probabilities3 = rev_probabilities_fun(n, lambda, freq.obs, R1, S1)$alternative2 # (see formula (27) in combination with (21))
                             freq.exp3 = expected_values_fun(n_cat1, n_cat2, probabilities3) # (see formula (28))
                             
                             
                             ## compute the observed optimal nu for the regularised estimators (see formula (16) in combination with (8))
                             v.opt.obs = v.opt_fun(n = n, q = q, freq = freq.obs, R = R1, S = S1)
                             
                             ## compute the expected optimal nu for the regularised estimators (using prior 1) (see formula (16) in combination with (28) and (18))
                             v.opt.exp1 = v.opt_fun(n = n, q = q, freq = freq.exp1, R = R1, S = S1)
                             
                             ## compute the expected optimal nu for the regularised estimators (using prior 2) (see formula (16) in combination with (28) and (19))
                             v.opt.exp2 = v.opt_fun(n = n, q = q, freq = freq.exp2, R = R1, S = S1)
                             
                             # compute the expected optimal nu for the regularised estimators (using prior 3) (see formula (16) in combination with (28) and (21))
                             v.opt.exp3 = v.opt_fun(n = n, q = q, freq = freq.exp3, R = R1, S = S1)
                             
                             # compute the observed optimal mu for the regularised estimators (see formula (A.33) in combination with (8)) !not included in thesis!
                             u.opt.obs = u.opt_fun(n = n, q = q, freq = freq.obs, R = R1, S = S1)
                             
                             # compute the observed optimal mu for the regularized estimators (using prior 1) (see formula (A.33) in combination with (28) and (18)) !not included in thesis! 
                             u.opt.exp1 = u.opt_fun(n = n, q = q, freq = freq.exp1, R = R1, S = S1)
                             
                             # compute the expected optimal mu for the regularized estimators (using prior 2) (see formula (A.33) in combination with (28) and (19)) !not included in thesis!
                             u.opt.exp2 = u.opt_fun(n = n, q = q, freq = freq.exp2, R = R1, S = S1)
                             
                             # compute the observed optimal mu for the regularized estimators  (using prior 3) (see formula (A.33) in combination with (28) and (21)) !not included in thesis!
                             u.opt.exp3 = u.opt_fun(n = n, q = q, freq = freq.exp3, R = R1, S = S1)
                             
                             # regularised estimators proposed by Frank Pijpers (2021) (see formula (15))
                             freq.reg.prag = reg.est1(n = n, q = q, Yg = Yg, v = v.prag, Zstar = Zstar)
                             freq.reg.opt.true = reg.est_opt1(n = n, q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.true, Zstar = Zstar)
                             freq.reg.opt.obs = reg.est_opt1(n = n, q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.obs, Zstar = Zstar)
                             freq.reg.opt.exp1 = reg.est_opt1(n = n, q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp1, Zstar = Zstar)
                             freq.reg.opt.exp2 = reg.est_opt1(n = n, q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp2, Zstar = Zstar)
                             freq.reg.opt.exp3 = reg.est_opt1(n = n, q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp3, Zstar = Zstar)
                             
                             # !the following six use the second variant of the regularised estimator and are not included in the thesis! (see formula (A.30))
                             freq.reg2.prag = reg.est2(n = n,  q = q, Yg = Yg, mu = u.prag, Zstar = Zstar)
                             freq.reg2.opt.true = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg, mu = u.opt.true, Zstar = Zstar)
                             freq.reg2.opt.obs = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg, mu = u.opt.obs, Zstar = Zstar)
                             freq.reg2.opt.exp1 = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg, mu = u.opt.exp1, Zstar = Zstar)
                             freq.reg2.opt.exp2 = reg.est_opt2(n = n, q = q,n_cat1, n_cat2, Yg = Yg, mu = u.opt.exp2, Zstar = Zstar)
                             freq.reg2.opt.exp3 = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg, mu = u.opt.exp3, Zstar = Zstar)
                             
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (18), and (14))
                             e_var_exp1 = var_fun(n = n, q = q, freq = freq.exp1, R = R1, S = S1)
                             emse_obs1 = ((-(1-q)/(n-1))*(n*freq.exp1-(R1*S1)))^2 + e_var_exp1
                             emse_Q1 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp1-(R1*S1)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp1
                             emse_Qinv1 = ((n-1)^2/(n*q-1)) * e_var_exp1
                             freq.weighed.emse1 = (((1/emse_Q1)*freq.Q) + ((1/emse_obs1)*freq.obs) + ((1/emse_Qinv1)*freq.Qinv)) / ((1/emse_Q1) + (1/emse_obs1) + (1/emse_Qinv1))
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (19), and (14))
                             e_var_exp2 = var_fun(n = n, q = q, freq = freq.exp2, R = R1, S = S1)
                             emse_obs2 = ((-(1-q)/(n-1))*(n*freq.exp2-(R1*S1)))^2 + e_var_exp2
                             emse_Q2 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp2-(R1*S1)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp2
                             emse_Qinv2 = ((n-1)^2/(n*q-1)) * e_var_exp2
                             freq.weighed.emse2 = (((1/emse_Q2)*freq.Q) + ((1/emse_obs2)*freq.obs) + ((1/emse_Qinv2)*freq.Qinv)) / ((1/emse_Q2) + (1/emse_obs2) + (1/emse_Qinv2))
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (21), and (14))
                             e_var_exp3 = var_fun(n = n, q = q, freq = freq.exp3, R = R1, S = S1)
                             emse_obs3 = ((-(1-q)/(n-1))*(n*freq.exp3-(R1*S1)))^2 + e_var_exp3
                             emse_Q3 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp3-(R1*S1)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp3
                             emse_Qinv3 = ((n-1)^2/(n*q-1)) * e_var_exp3
                             freq.weighed.emse3 = (((1/emse_Q3)*freq.Q) + ((1/emse_obs3)*freq.obs) + ((1/emse_Qinv3)*freq.Qinv)) / ((1/emse_Q3) + (1/emse_obs3) + (1/emse_Qinv3))
                             
                             
                             # compute differences of contingency table per estimator with the true contingency table freq1
                             dev.res = array(NA, dim = c(dim(freq1), 21))
                             dev.res[,,1] = freq.obs - freq1
                             dev.res[,,2] = freq.Q - freq1
                             dev.res[,,3] = freq.Qinv - freq1
                             dev.res[,,4] = freq.exp1 - freq1
                             dev.res[,,5] = freq.exp2 - freq1
                             dev.res[,,6] = freq.exp3 - freq1
                             dev.res[,,7] = freq.reg.prag - freq1
                             dev.res[,,8] = freq.reg.opt.true - freq1
                             dev.res[,,9] = freq.reg.opt.obs - freq1
                             dev.res[,,10] = freq.reg.opt.exp1 - freq1
                             dev.res[,,11] = freq.reg.opt.exp2 - freq1
                             dev.res[,,12] = freq.reg.opt.exp3 - freq1
                             dev.res[,,13] = freq.reg2.prag - freq1
                             dev.res[,,14] = freq.reg2.opt.true - freq1
                             dev.res[,,15] = freq.reg2.opt.obs - freq1
                             dev.res[,,16] = freq.reg2.opt.exp1 - freq1
                             dev.res[,,17] = freq.reg2.opt.exp2 - freq1
                             dev.res[,,18] = freq.reg2.opt.exp3 - freq1
                             dev.res[,,19] = freq.weighed.emse1 - freq1
                             dev.res[,,20] = freq.weighed.emse2 - freq1
                             dev.res[,,21] = freq.weighed.emse3 - freq1
                             
                             return(dev.res)
                           }, simplify = 'array')
                           
                           
                           ## add names to first 3 dimensions with the categories of y, the categories of z and the names of the estimators, respectively.
                           dimnames(dev) = list(row.names(freq1),
                                                colnames(freq1),
                                                c('obs','Q','Qinv',
                                                  'exp1','exp2','exp3',
                                                  'reg1.prag','reg1.opt.true','reg1.opt.obs',
                                                  'reg1.opt.exp1','reg1.opt.exp2','reg1.opt.exp3',
                                                  'reg2.prag','reg2.opt.true','reg2.opt.obs',
                                                  'reg2.opt.exp1','reg2.opt.exp2','reg2.opt.exp2',
                                                  'weighted.emse1','weighted.emse2','weighted.emse3'),
                                                NULL)
                           
                           
                           #########################################
                           ### Computations for independent attributes
                           
                           v.prag = v.prag_fun(n = n, q = q) # (see formula (17))
                           v.opt.true = v.opt_fun(n = n, q = q, freq = freq2, R = R2, S = S2) # (see formula (16))
                           
                           u.prag = u.prag_fun(n = n, q = q) # (see formula (A.34)) !not included in thesis!
                           u.opt.true = u.opt_fun(n = n, q = q, freq = freq2, R = R2, S = S2) # (see formula (A.33)) !not included in thesis!
                           
                           dev2 = sapply(Clist, function(C) {
                             
                             # Simulate linkage errors
                             Wstar = C %*% Wg  
                             
                             # estimators considered by Scholtus, Shlomo and De Waal (2022)
                             freq.obs = t(Yg) %*% Wstar # (see formula (8))
                             freq.Q = t(Q %*% Yg) %*% Wstar # == reg.est(0) (see formula (9))
                             freq.Qinv = t(Yg) %*% (Qinv %*% Wstar) # == reg.est(1) (see formula (10))
                          
                             # Set negative values to 0 for Qinv
                             freq.Qinv = ifelse(freq.Qinv < 0, 0, freq.Qinv)
                             
                             # Expected value by means of probabilities using alternative prior where all possible values have the same probability
                             probabilities1 = rev_probabilities_fun(n, lambda, freq.obs, R2, S2)$alternative1 # (see formula (27) in combination with (18))
                             freq.exp1 = expected_values_fun(n_cat1, n_cat2, probabilities1) # (see formula (28))
                             
                             # Expected value by means of probabilities using proposed prior by Scholtus & de Waal (2020-2022)
                             probabilities2 = rev_probabilities_fun(n, lambda, freq.obs, R2, S2)$proposed  # (see formula (27) in combination with (19))
                             freq.exp2 = expected_values_fun(n_cat1, n_cat2, probabilities2) # (see formula (28))
                             
                             # Expected value by means of probabilities using alternative prior by de Waal (2023) where a truncated binomial model is used
                             probabilities3 = rev_probabilities_fun(n, lambda, freq.obs, R2, S2)$alternative2  # (see formula (27) in combination with (21))
                             freq.exp3 = expected_values_fun(n_cat1, n_cat2, probabilities3) # (see formula (28))
                             
                             
                             ## compute the observed optimal nu for the regularised estimators (see formula (16) in combination with (8))
                             v.opt.obs = v.opt_fun(n = n, q = q, freq = freq.obs, R = R2, S = S2)
                             
                             ## compute the expected optimal nu for the regularised estimators (using prior 1) (see formula (16) in combination with (28) and (18))
                             v.opt.exp1 = v.opt_fun(n = n, q = q, freq = freq.exp1, R = R2, S = S2)
                             
                             ## compute the expected optimal nu for the regularised estimators (using prior 2) (see formula (16) in combination with (28) and (19))
                             v.opt.exp2 = v.opt_fun(n = n, q = q, freq = freq.exp2, R = R2, S = S2)
                             
                             # compute the expected optimal nu for the regularised estimators (using prior 3) (see formula (16) in combination with (28) and (21))
                             v.opt.exp3 = v.opt_fun(n = n, q = q, freq = freq.exp3, R = R2, S = S2)
                             
                             # compute the observed optimal mu for the regularised estimators (see formula (A.33) in combination with (8)) !not included in thesis!
                             u.opt.obs = u.opt_fun(n = n, q = q, freq = freq.obs, R = R2, S = S2)
                             
                             # compute the observed optimal mu for the regularized estimators (using prior 1) (see formula (A.33) in combination with (28) and (18)) !not included in thesis! 
                             u.opt.exp1 = u.opt_fun(n = n, q = q, freq = freq.exp1, R = R2, S = S2)
                             
                             # compute the expected optimal mu for the regularized estimators (using prior 2) (see formula (A.33) in combination with (28) and (19)) !not included in thesis!
                             u.opt.exp2 = u.opt_fun(n = n, q = q, freq = freq.exp2, R = R2, S = S2)
                             
                             # compute the observed optimal mu for the regularized estimators  (using prior 3) (see formula (A.33) in combination with (28) and (21)) !not included in thesis!
                             u.opt.exp3 = u.opt_fun(n = n, q = q, freq = freq.exp3, R = R2, S = S2)
                             
                             # regularised estimators proposed by Frank Pijpers (2021) (see formula (15))
                             freq.reg.prag = reg.est1(n = n, q = q, Yg = Yg, v = v.prag, Zstar = Wstar)
                             freq.reg.opt.true = reg.est_opt1(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.true, Zstar = Wstar)
                             freq.reg.opt.obs = reg.est_opt1(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.obs, Zstar = Wstar)
                             freq.reg.opt.exp1 = reg.est_opt1(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp1, Zstar = Wstar)
                             freq.reg.opt.exp2 = reg.est_opt1(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp2, Zstar = Wstar)
                             freq.reg.opt.exp3 = reg.est_opt1(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, v = v.opt.exp3, Zstar = Wstar)
                             
                             # !the following six use the second variant of the regularised estimator and are not included in the thesis! (see formula (A.30))
                             freq.reg2.prag = reg.est2(n = n,  q = q, Yg = Yg, mu = u.prag, Zstar = Wstar)
                             freq.reg2.opt.true = reg.est_opt2(n = n,  q = q, n_cat1, n_cat2, Yg = Yg, mu = u.opt.true, Zstar = Wstar)
                             freq.reg2.opt.obs = reg.est_opt2(n = n,  q = q, n_cat1, n_cat2, Yg = Yg,mu = u.opt.obs, Zstar = Wstar)
                             freq.reg2.opt.exp1 = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg,mu = u.opt.exp1, Zstar = Wstar)
                             freq.reg2.opt.exp2 = reg.est_opt2(n = n,  q = q,n_cat1, n_cat2, Yg = Yg,mu = u.opt.exp2, Zstar = Wstar)
                             freq.reg2.opt.exp3 = reg.est_opt2(n = n, q = q, n_cat1, n_cat2, Yg = Yg,mu = u.opt.exp3, Zstar = Wstar)
                             
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (18), and (14))
                             e_var_exp1 = var_fun(n = n, q = q, freq = freq.exp1, R = R2, S = S2)
                             emse_obs1 = ((-(1-q)/(n-1))*(n*freq.exp1-(R2*S2)))^2 + e_var_exp1
                             emse_Q1 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp1-(R2*S2)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp1
                             emse_Qinv1 = ((n-1)^2/(n*q-1)) * e_var_exp1
                             freq.weighed.emse1 = (((1/emse_Q1)*freq.Q) + ((1/emse_obs1)*freq.obs) + ((1/emse_Qinv1)*freq.Qinv)) / ((1/emse_Q1) + (1/emse_obs1) + (1/emse_Qinv1))
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (19), and (14))
                             e_var_exp2 = var_fun(n = n, q = q, freq = freq.exp2, R = R2, S = S2)
                             emse_obs2 = ((-(1-q)/(n-1))*(n*freq.exp2-(R2*S2)))^2 + e_var_exp2
                             emse_Q2 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp2-(R2*S2)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp2
                             emse_Qinv2 = ((n-1)^2/(n*q-1)) * e_var_exp2
                             freq.weighed.emse2 = (((1/emse_Q2)*freq.Q) + ((1/emse_obs2)*freq.obs) + ((1/emse_Qinv2)*freq.Qinv)) / ((1/emse_Q2) + (1/emse_obs2) + (1/emse_Qinv2))
                             
                             # Weight the three correction methods Q, Qinv and naive with their MSE (see formula (29) in combination with (28), (21), and (14))
                             e_var_exp3 = var_fun(n = n, q = q, freq = freq.exp3, R = R2, S = S2)
                             emse_obs3 = ((-(1-q)/(n-1))*(n*freq.exp3-(R2*S2)))^2 + e_var_exp3
                             emse_Q3 = (-(1+((n*q-1)/(n-1))) * ((1-q)/(n-1)) * (n*freq.exp3-(R2*S2)))^2 *  + ((n*q-1)^2/(n-1)^2) * e_var_exp3
                             emse_Qinv3 = ((n-1)^2/(n*q-1)) * e_var_exp3
                             freq.weighed.emse3 = (((1/emse_Q3)*freq.Q) + ((1/emse_obs3)*freq.obs) + ((1/emse_Qinv3)*freq.Qinv)) / ((1/emse_Q3) + (1/emse_obs3) + (1/emse_Qinv3))
                             
                             
                             # compute differences of contingency table per estimator with the true contingency table freq2
                             dev.res = array(NA, dim = c(dim(freq2), 21))
                             dev.res[,,1] = freq.obs - freq2
                             dev.res[,,2] = freq.Q - freq2
                             dev.res[,,3] = freq.Qinv - freq2
                             dev.res[,,4] = freq.exp1 - freq2
                             dev.res[,,5] = freq.exp2 - freq2
                             dev.res[,,6] = freq.exp3 - freq2
                             dev.res[,,7] = freq.reg.prag - freq2
                             dev.res[,,8] = freq.reg.opt.true - freq2
                             dev.res[,,9] = freq.reg.opt.obs - freq2
                             dev.res[,,10] = freq.reg.opt.exp1 - freq2
                             dev.res[,,11] = freq.reg.opt.exp2 - freq2
                             dev.res[,,12] = freq.reg.opt.exp3 - freq2
                             dev.res[,,13] = freq.reg2.prag - freq2
                             dev.res[,,14] = freq.reg2.opt.true - freq2
                             dev.res[,,15] = freq.reg2.opt.obs - freq2
                             dev.res[,,16] = freq.reg2.opt.exp1 - freq2
                             dev.res[,,17] = freq.reg2.opt.exp2 - freq2
                             dev.res[,,18] = freq.reg2.opt.exp3 - freq2
                             dev.res[,,19] = freq.weighed.emse1 - freq2
                             dev.res[,,20] = freq.weighed.emse2 - freq2
                             dev.res[,,21] = freq.weighed.emse3 - freq2
                             
                             return(dev.res)
                           }, simplify = 'array')
                           
                           ## add names to first 3 dimensions with the categories of y, the categories of z and the names of the estimators, respectively.
                           dimnames(dev2) = list(row.names(freq2),
                                                 colnames(freq2),
                                                 c('obs','Q','Qinv',
                                                   'exp1','exp2','exp3',
                                                   'reg1.prag','reg1.opt.true','reg1.opt.obs',
                                                   'reg1.opt.exp1','reg1.opt.exp2','reg1.opt.exp3',
                                                   'reg2.prag','reg2.opt.true','reg2.opt.obs',
                                                   'reg2.opt.exp1','reg2.opt.exp2','reg2.opt.exp3',
                                                   'weighted.emse1','weighted.emse2','weighted.emse3'),
                                                 NULL)
                           
                           dev_list1[[i]] = dev
                           dev_list2[[i]] = dev2
                         }
                         return(list(dev_list1,dev_list2))
                       }
stopCluster(cl) 
####################################
### End of parallel simulation study


###############
## Save results

setwd("S:/Results") # set directory to where you want to save the results
save('res_parallel',
     file = sprintf("Simulation_results.Rdata"))

