############################################################################################################
### This R-code makes data frames and graphs from the simulation study belonging to the thesis:          ###
### Correcting for Linkage Errors in Contingency Tables: New Methods to Improve the Correction Approach  ###
### by Sjarai Eikenhout. These visualizations and results can be found in section 4.2 and Appendix E.    ###           
### Comments are added to the code to explain what is happening.                                         ###
############################################################################################################


#################
### Load packages

library(dummies)
library(ggplot2) 
library(latex2exp)


#############################################################################################
### Load the results and set the following variables the same as in the performed simulation.
### Note: this is not necessary if the simulation is performed in the same R-session and all
### the results are still in the global environment. In that case, lines 20-92 can be skipped.

load("S:/Results/Simulation_results.Rdata") # Copy the path to the data here

n_tab = 10 # number of tables generated
n = 300 # number of records in data sets
n_sim = 10000 # number of simulations
n_cat1 = 5 # number of categories in data set 1
n_cat2 = 5 # number of categories in data set 2
qvec = seq(0.9, 0.1, by = -0.1) # values of q

set.seed(12345) # use same seed as in Simulation_study.R

## auxiliary variable x
x = rnorm(n*n_tab, 20, 10)

## target variable z
z = rnorm(n*n_tab, 20, 15)

## target variable y
y = 20 + 2*x + 1*z + rnorm(n *n_tab, 0, 4)

## variable w
w = rnorm(n*n_tab, 40, 30)


## Discretize variables
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
  
  
  ##################################
  ### Compute true contingency table
  
  freq1list[[i]] = t(Ylist[[i]]) %*% Zlist[[i]] # dependent attributes
  freq2list[[i]] = t(Ylist[[i]]) %*% Wlist[[i]] # independent attributes
  
  
  ##########################
  ### Determine associations
  
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


###################################################################################
### Analyse results as done in Scholtus et al. (2021)  by looking at the following:
### 1) The percentages that the naive approach outperforms the alternative approach
### 2) The total relative difference from the naive empirical RMSE 

## Cell positions
cell_pos = c()
for (i in 1:n_cat1){
  for (j in 1:n_cat2){
    cell_pos = rbind(cell_pos, paste0("R",i," x C",j))
  }
}

## Make multi-dimensional arrays to store the results in n_tab tables
bias_dependent = mse_dependent = var_dependent = rmse_dependent = array(NA, dim = (c(dim(freq1list[[1]]), 21, length(qvec), n_tab)))
bias_independent = mse_independent = var_independent = rmse_independent = array(NA, dim = (c(dim(freq1list[[1]]), 21, length(qvec), n_tab)))
RMSE_diff_dep = RMSE_diff_indep = array(NA, dim = (c(dim(freq1list[[1]]), 20, length(qvec), n_tab)))
averages_dependent = differences_dependent = averages_independent = differences_independent = 
  array(NA, dim = c(length(qvec), 20, n_tab), dimnames = list(paste("q =",qvec),
                                                              c('Q','Qinv',
                                                                'exp1','exp2','exp3',
                                                                'reg1.prag','reg1.opt.true','reg1.opt.obs',
                                                                'reg1.opt.exp1','reg1.opt.exp2','reg1.opt.exp3',
                                                                'reg2.prag','reg2.opt.true','reg2.opt.obs',
                                                                'reg2.opt.exp1','reg2.opt.exp2','reg2.opt.exp3',
                                                                'weighted.emse1','weighted.emse2','weighted.emse3'),
                                                              paste("Table", 1:n_tab)))

for (i in 1:n_tab){
  for (j in 1:length(qvec)){
    
    ## Compute bias, mse, var, and rmse for q = j, table i with dependent attributes
    bias_dependent[,,,j,i] = apply(res_parallel[[j]][[1]][[i]], 1:3, mean)
    mse_dependent[,,,j,i] = apply(res_parallel[[j]][[1]][[i]]^2, 1:3, mean)
    var_dependent[,,,j,i] = mse_dependent[,,,j,i] - bias_dependent[,,,j,i]
    rmse_dependent[,,,j,i] = sqrt(mse_dependent[,,,j,i])
    
    ## Compute bias, mse, var, and rmse for q = j, table i with independent attributes
    bias_independent[,,,j,i] = apply(res_parallel[[j]][[2]][[i]], 1:3, mean)
    mse_independent[,,,j,i] = apply(res_parallel[[j]][[2]][[i]]^2, 1:3, mean)
    var_independent[,,,j,i] = mse_independent[,,,j,i] - bias_independent[,,,j,i]
    rmse_independent[,,,j,i] = sqrt(mse_independent[,,,j,i])
    
    ## Percentages that the naive approach outperforms the alternative approach for tables with dependent attributes and q = j
    performance_Q = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,2,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2))/ n_sim * 100 
    performance_Qinv = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,3,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,4,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,5,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,6,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.prag = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,7,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.true = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,8,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.obs = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,9,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,10,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,11,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,12,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.prag = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,13,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.true = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,14,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.obs = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,15,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,16,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,17,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,18,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse1 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,19,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse2 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,20,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse3 = rowSums(matrix(c(abs(res_parallel[[j]][[1]][[i]][,,21,]) > abs(res_parallel[[j]][[1]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100
    performances_dep = data.frame(cell_pos, ass_deplist[[i]], 
                                  performance_Q, performance_Qinv,
                                  performance_exp1, performance_exp2, performance_exp3,
                                  performance_reg1.prag, performance_reg1.opt.true, performance_reg1.opt.obs,
                                  performance_reg1.opt.exp1, performance_reg1.opt.exp2, performance_reg1.opt.exp3,
                                  performance_reg2.prag, performance_reg2.opt.true, performance_reg2.opt.obs,
                                  performance_reg2.opt.exp1, performance_reg2.opt.exp2, performance_reg2.opt.exp3,
                                  performance_weighed.emse1, performance_weighed.emse2, performance_weighed.emse3)
    averages_dependent[j,,i] = colSums(performances_dep[,3:22])/(n_cat1*n_cat2) 
    
    ## Total relative difference from the naive empirical RMSE for dependent tables and q = j
    diff = ((rmse_dependent[,,2:21,j,i] - rep(rmse_dependent[,,1,j,i],20)))/ rep(rmse_dependent[,,1,j,i],20)
    RMSE_diff_dep[,,,j,i] = diff
    differences_dependent[j,,i] = apply(diff, 3, sum)
    
    ## Percentages that the naive approach outperforms the alternative approach for tables with independent attributes and q = j
    performance_Q = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,2,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2))/ n_sim * 100 
    performance_Qinv = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,3,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,4,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,5,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,6,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.prag = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,7,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.true = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,8,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.obs = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,9,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,10,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,11,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg1.opt.exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,12,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.prag = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,13,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.true = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,14,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.obs = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,15,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp1 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,16,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp2 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,17,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_reg2.opt.exp3 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,18,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse1 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,19,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse2 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,20,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100 
    performance_weighed.emse3 = rowSums(matrix(c(abs(res_parallel[[j]][[2]][[i]][,,21,]) > abs(res_parallel[[j]][[2]][[i]][,,1,])), nrow = n_cat1 * n_cat2)) / n_sim * 100
    performances_indep = data.frame(cell_pos, ass_indeplist[[i]], 
                                    performance_Q, performance_Qinv,
                                    performance_exp1, performance_exp2, performance_exp3,
                                    performance_reg1.prag, performance_reg1.opt.true, performance_reg1.opt.obs,
                                    performance_reg1.opt.exp1, performance_reg1.opt.exp2, performance_reg1.opt.exp3,
                                    performance_reg2.prag, performance_reg2.opt.true, performance_reg2.opt.obs,
                                    performance_reg2.opt.exp1, performance_reg2.opt.exp2, performance_reg2.opt.exp3,
                                    performance_weighed.emse1, performance_weighed.emse2, performance_weighed.emse3)
    averages_independent[j,,i] = colSums(performances_indep[,3:22])/(n_cat1*n_cat2) 
    
    ## Total relative difference from the naive empirical RMSE for independent tables  and q = j
    diff = ((rmse_independent[,,2:21,j,i] - rep(rmse_independent[,,1,j,i],20)))/ rep(rmse_independent[,,1,j,i],20)
    RMSE_diff_indep[,,,j,i] = diff
    differences_independent[j,,i] = apply(diff, 3, sum)
  }
}


#####################################################################################################
### Make data frame with the average percentages where the naive approach outperforms the alternative
### approach over the 10 tables with dependent attributes for each value of q

means_dep = apply(averages_dependent, 1:2, mean) # compute average percentage over the 10 tables with dependent attributes (freq1)
means_sd = apply(averages_dependent, 1:2, sd) # compute standard deviation over the 10 tables with dependent attributes (freq1)
rownames(means_dep) = rownames(means_sd) = qvec # add row names
df_means_dep = as.data.frame.table(means_dep)
df_sd_dep =  as.data.frame.table(means_sd)
df_means_dep = cbind(df_means_dep, df_sd_dep[,3])
colnames(df_means_dep) = c("q", "Estimator", "average", "sd") # add column names
drop = c("reg2.prag", "reg2.opt.true", "reg2.opt.obs", 
         "reg2.opt.exp1", "reg2.opt.exp2", "reg2.opt.exp3") # remove results second variant of the regularised estimators as they are not included in thesis
df_means_dep = df_means_dep[!(df_means_dep$Estimator %in% drop),]
df_means_dep$q = as.numeric(qvec)


#####################################################################################################
### Make data frame with the average percentages where the naive approach outperforms the alternative
### approach over the 10 tables with independent attributes for each value of q

means_indep = apply(averages_independent, 1:2, mean) # compute average percentage over the 10 tables with independent attributes (freq2)
means_sd = apply(averages_independent, 1:2, sd) # compute standard deviation over the 10 tables with independent attributes (freq2)
rownames(means_dep) = rownames(means_sd) = qvec # add row names
df_means_indep = as.data.frame.table(means_indep)
df_sd_indep =  as.data.frame.table(means_sd)
df_means_indep = cbind(df_means_indep, df_sd_indep[,3])
colnames(df_means_indep) = c("q", "Estimator", "average", "sd") # add column names
drop = c("reg2.prag", "reg2.opt.true", "reg2.opt.obs", 
         "reg2.opt.exp1", "reg2.opt.exp2", "reg2.opt.exp3")  # remove second variant of the regularised estimators
df_means_indep = df_means_indep[!(df_means_indep$Estimator %in% drop),]
df_means_indep$q = as.numeric(qvec)


############################################################################################
### Make data frame with the average total relative difference from the naive empirical RMSE
### over the 10 tables with dependent attributes for each value of q

diff_dep = apply(differences_dependent, 1:2, mean) # compute average total relative difference over the 10 tables with dependent attributes (freq1)
diff_sd = apply(differences_dependent, 1:2, sd) # compute standard deviation over the 10 tables with dependent attributes (freq1)
rownames(diff_dep) = rownames(diff_sd) = qvec # add row names
df_diff_dep = as.data.frame.table(diff_dep)
df_sd_dep =  as.data.frame.table(diff_sd)
df_diff_dep = cbind(df_diff_dep, df_sd_dep[,3])
colnames(df_diff_dep) = c("q", "Estimator", "diff", "sd") # add column names
drop = c("reg2.prag", "reg2.opt.true", "reg2.opt.obs", 
         "reg2.opt.exp1", "reg2.opt.exp2", "reg2.opt.exp3") # remove second variant of the regularised estimators
df_diff_dep = df_diff_dep[!(df_diff_dep$Estimator %in% drop),]
df_diff_dep$q = as.numeric(qvec)


############################################################################################
### Make data frame with the average total relative difference from the naive empirical RMSE
### over the 10 tables with independent attributes for each value of q

diff_indep = apply(differences_independent, 1:2, mean) # compute average total relative difference over the 10 tables with independent attributes (freq2)
diff_sd = apply(differences_independent, 1:2, sd) # compute standard deviation over the 10 tables with independent attributes (freq2)
rownames(diff_indep) = rownames(diff_sd) = qvec # add row names
df_diff_indep = as.data.frame.table(diff_indep)
df_sd_indep =  as.data.frame.table(diff_sd)
df_diff_indep = cbind(df_diff_indep, df_sd_indep[,3])
colnames(df_diff_indep) = c("q", "Estimator", "diff", "sd") # add column names
drop = c("reg2.prag", "reg2.opt.true", "reg2.opt.obs", 
         "reg2.opt.exp1", "reg2.opt.exp2", "reg2.opt.exp3") # remove second variant of the regularised estimators
df_diff_indep = df_diff_indep[!(df_diff_indep$Estimator %in% drop),]
df_diff_indep$q = as.numeric(qvec)


###################
### Plots in thesis


## set colors
colors = c("#FF0000", "#0008FF", 
           "#00FFD8", "#FFF700", "#FF6BF8", 
           "#1FFF00", "#32CAFF", "#FF9B00", "#FF0064", "#D49DFF", "#F7FF00",
           "#E400FF", "#00960E", "#CC0000")

## set labels
labels = unname(TeX(c("$\\textbf{Q}$", "$\\textbf{Q}^{-1}$",
                      "$E_1$", "$E_2$", "$E_3$",
                      "$Reg_{prag}$", "$Reg_{opt(true)}$", "$Reg_{opt(obs)}$", "$Reg_{opt(E_1)}$", "$Reg_{opt(E_2)}$", "$Reg_{opt(E_3)}$",
                      "$W_{\\widehat{MSE}_1}$", "$W_{\\widehat{MSE}_2}$", "$W_{\\widehat{MSE}_3}$")))

## Figure 2
ggplot(df_means_dep, aes(x = q, y = average)) + 
  ylab("Average percentage (%)") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = average, ymin = average - sd, ymax = average + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  xlim(0,1) + 
  ylim(0,100) +
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))

## Figure 3
ggplot(df_means_indep, aes(x = q, y = average)) + 
  ylab("Average percentage (%)") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = average, ymin = average - sd, ymax = average + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  xlim(0,1) + 
  ylim(0,100) +
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))

## Figure 4
ggplot(df_diff_dep, aes(x = q, y = diff)) + 
  ylab("Average total relative difference") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = diff, ymin = diff - sd, ymax = diff + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  coord_cartesian(ylim = c(-10, 30)) +
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))

## Figure 5
ggplot(df_diff_indep, aes(x = q, y = diff)) + 
  ylab("Average total relative difference") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = diff, ymin = diff - sd, ymax = diff + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  coord_cartesian(ylim = c(-10, 30)) +
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))

## Figure E.6
ggplot(df_diff_dep, aes(x = q, y = diff)) + 
  ylab("Average total relative difference") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = diff, ymin = diff - sd, ymax = diff + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))

## Figure E.7
ggplot(df_diff_indep, aes(x = q, y = diff)) + 
  ylab("Average total relative difference") +
  geom_line(aes(color = Estimator)) +
  geom_point(aes(color = Estimator, shape = Estimator, size = Estimator)) +
  theme_bw() +
  scale_color_manual(values = colors, labels = labels, name = "Kleur", aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(15, 15, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18), labels = labels, name = "Kleur") +
  scale_size_manual(values = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), labels = labels, name = "Kleur") +
  geom_ribbon(aes(y = diff, ymin = diff - sd, ymax = diff + sd, fill = Estimator),          
              alpha = 0.05, color = NA) +
  guides(fill = "none") + 
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20))
