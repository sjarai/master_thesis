##################################################################################################
### This R-code contains all the functions that are needed for the simulation study performed  ###
### in the thesis: Correcting for Linkage Errors in Contingency Tables: New Methods to Improve ###
### the Correction Approach by Sjarai Eikenhout.                                               ###
### To be able to reproduce the simulation study in the thesis, this code should be ran first. ###
### This is automatically done in the parallel simulation study in Simulation_study.R          ###
### Comments are added to the code to explain what is happening.                               ###
### The corresponding formulas in the thesis are also given as comments for each function.     ###
##################################################################################################


#######################################################################################
### Functions to generate permutations matrices with proposed method by Scholtus (2020)
### Described in Appendix B

## Compute lambda (see formula (B.35))
lambdaExact = function(n, q) { 
  ((n-1)/2) * (log(n-1) - log(n*q - 1))
}

## Draw permutation matrix C
drawPerm = function(K, n) {
  C0 = diag(n)
  transSet = permutations(n = n, r = 2)
  L = nrow(transSet)
  
  if (K <= 0) return(C0)
  
  C = C0
  for (k in 1:K) {
    s = sample(x = L, size = 1)
    trans = transSet[s,]
    Ctrans = C0
    Ctrans[trans[1], trans] = rev(Ctrans[trans[1], trans])
    Ctrans[trans[2], trans] = rev(Ctrans[trans[2], trans])
    C = C %*% Ctrans
  }
  
  return(C)
}

## Draw K first and thereafter C
drawC = function(lambda, n) {
  K = rpois(n = 1, lambda = lambda) 
  drawPerm(K = K, n = n)
}


###############################################################
### Variance of contingency table - function (see formula (14))

var_fun = function(n, q, freq, R, S){
  q * (1-q) * freq + (((1-q) / (n-1)) * (1 - ((1-q) / (n-1))) * (R*S - freq)) +
    (-q) * ((1-q)/(n-1)) * 2 * freq * (R + S - 2) +
    ((((2*q) / (n-1)) + ((1 - 2/(n-1))*((n*q-1)/(n-1))^(2-2/n)) - q^2 - ((1-q)/(n-1))^2) * freq * (freq - 1)) +
    ((1/(n-2)) * ((((n-3)/(n-1)) * (q^2 - ((n*q-1)/(n-1))^(2-2/n))) + ((1-q)/(n-1))^2)) * 2 * freq * ((R-1)*(S-1)-(freq-1)) +
    (-(((1-q)/(n-1))^2 * (R*S - 2*freq)*(R+S-2))) +
    (((1/(n-2)) * (((1/(n-1))*(((n*q-1)/(n-1))^(2-2/n)-q^2)) + ((1-q)/(n-1))^2)) * ((R*S*(R-1)*(S-1)) - (2*freq*(2*(R-1)*(S-1)-(freq-1)))))
}


#################################################################################
### First variant of the regularised estimator by Pijpers (2021) from section 3.2
### Note: nu is coded as "v"

## pragmatic nu (see formula (17))
v.prag_fun = function(n, q){
  (n*q - 1)^2/((n-1)^2 + (n*q - 1)^2) 
}

## optimal nu (see formula (16))
v.opt_fun = function(n, q, freq, R, S){
  var.obs = var_fun(n = n, q = q, freq = freq, R = R, S = S)
  v = 1 - (((n^2 * (n-1)^2)/((n-1)^2 - (n*q-1)^2)) * (var.obs / (R*S - n*freq)^2))
  v[v < 0] = 0 # since 0 <= nu <= 1
  v[v > 1] = 1 # since 0 <= nu <= 1
  return(v)
}

## family of regularized estimators variant 1 where nu is a scalar
## (see formulas (15) and (17) in thesis in combination with formulas (15) and (17) in Pijpers (2021))
reg.est1 = function(n, q, Yg, v, Zstar){
  lambda = ((n-1)^2-(n*q-1)^2)/(n*(n-1)*(1-q))
  u = c(rep(1, n))
  t(Yg) %*% ((((n*q)-1)/(n-1-(n*v*lambda*(1-q))) * diag(n)) + (((1-q)/(n-1-(n*v*lambda*(1-q))))*(1-(v*lambda))*u%*%t(u))) %*% Zstar
}

## family of regularized estimators variant 1 where nu is a matrix with unique values for each element 
## (see formula (15) and (16) in combination with formulas (15) and (17) in Pijpers (2021))
reg.est_opt1 = function(n, q, n_cat1, n_cat2, Yg, v, Zstar){
  out = matrix(nrow = n_cat1, ncol = n_cat2)
  lambda = ((n-1)^2-(n*q-1)^2)/(n*(n-1)*(1-q))
  u = c(rep(1, n))
  for (i in 1:n_cat1){
    for (j in 1:n_cat2){
      out[i,j] = Yg[ ,i] %*% ((((n*q)-1)/(n-1-(n*v[i,j]*lambda*(1-q))) * diag(n)) +(((1-q)/(n-1-(n*v[i,j]*lambda*(1-q))))*(1-(v[i,j]*lambda))*u%*%t(u))) %*% Zstar[ ,j]
    }
  }
  return(out)
}


###################################################################################################
### Functions for the second variant of the regularised estimator by Pijpers (2021) from appendix A
### Note: mu is coded as "u" + not included in study and thesis

## pragmatic mu - function (see formula (A.34))
u.prag_fun = function(n, q){
  (n*q - 1)/(n*q + n - 2) 
}

## optimal mu - function (see formula (A.33))
u.opt_fun = function(n, q, freq, R, S){
  var.obs = var_fun(n = n, q = q, freq = freq, R = R, S = S)
  u = 1 - (((n-1))^2/((n*q-1)*(1-q)^2) * (var.obs / (R*S - n*freq)^2))
  u[u < 0] = 0 # since 0 <= mu <= 1
  u[u > 1] = 1 # since 0 <= mu <= 1
  return(u)
}

## family of regularized estimators variant 2 where mu is a scalar
## (see formula (A.30) in combination with formula (A.34))
reg.est2 = function(n, q, Yg, mu, Zstar) {
  u = c(rep(1, n))
  t(Yg) %*% ((((n-1)/((mu*n*(q-1)) + (n-1))) * diag(n)) - (((mu*(1-q))/((mu*n*(q-1))+(n-1)))*u%*%t(u))) %*% Zstar
}

## family of regularized estimators variant 2 where mu is a matrix with unique values for each element 
## (see formula (A.30) in combination with formula (A.33))
reg.est_opt2 = function(n, q, n_cat1, n_cat2, Yg, mu, Zstar){
  out = matrix(nrow = n_cat1, ncol = n_cat2)
  u = c(rep(1, n))
  for (i in 1:n_cat1){
    for (j in 1:n_cat2){
      out[i,j] = Yg[ ,i] %*% ((((n-1)/(mu[i,j]*n*(q-1)+(n-1))) * diag(n)) - ((mu[i,j]*(1-q))/(mu[i,j]*n*(q-1)+(n-1))*u%*%t(u))) %*% Zstar[ ,j]
    }
  }
  return(out)
}


####################################################################################################
### Functions needed for the computations of the probabilities Pr(t*_jk | t = t_jk) from section 3.3

## + function (see formula (23))
plus_fun = function(n, r, s, t_star){
  2*((r-t_star)*(s-t_star))/(n*(n-1))
}

## - function (see formula (24))
minus_fun = function(n, r, s, t_star){
  2*(t_star*(n-r-s+t_star))/(n*(n-1))
}

## = function (see formula (25))
equal_fun = function(n, r, s, t_star){
  1 - 2 *((r*s + (n*t_star) - (2*(r+s)*t_star) + (2*t_star^2))/(n*(n-1)))
}

#####################################################################################
### Compute probabilities using Bayes' rule for each individual cell (see appendix D)
### In one step -> P is computed once and 3 different priors are used
### Note that the order is different here than in the thesis: 
### The proposed prior in the code is the second prior in the thesis.
### The first alternative prior in the code is the first prior in the thesis.
### The second alternative prior in the code is the third prior in the thesis.

rev_probabilities_fun = function(n, lambda, freq, R, S, range = 19){
  probabilities_prop = array(NA, dim = c(unique(sort(c(R,S)))[length(unique(sort(c(R,S))))-1]+1, 2, dim(freq)))
  probabilities_alt1 = array(NA, dim = c(unique(sort(c(R,S)))[length(unique(sort(c(R,S))))-1]+1, 2, dim(freq)))
  probabilities_alt2 = array(NA, dim = c(unique(sort(c(R,S)))[length(unique(sort(c(R,S))))-1]+1, 2, dim(freq)))
  
  for (ii in 1:nrow(freq)){
    for (jj in 1:ncol(freq)){
      t_obs = freq[ii,jj]
      r = R[ii,jj]
      s = S[ii,jj]
      
      L = max(0, r + s - n)
      U = min(r, s)
      t_star = seq(L, U)
      
      if (length(t_star) > 2*range+1){
        t_star_min = max(0, t_obs - range)
        t_star_max = min(U, t_obs + range)
        
        t_star = seq(t_star_min, t_star_max)
      }
      
      trans = data.frame(t_star = t_star)
      trans$plus = plus_fun(n, r, s, trans$t_star)
      trans$equal = equal_fun(n, r, s, trans$t_star)
      trans$min = minus_fun(n, r, s, trans$t_star)
      
      # Make matrix A (see formula (D.36))
      A = matrix(0, nrow = length(t_star), ncol = length(t_star))
      for (x in 2:(length(t_star))) A[x, x-1] = trans$plus[x-1]
      for (x in 1:(length(t_star))) A[x, x] = trans$equal[x]
      for (x in 1:(length(t_star)-1)) A[x, x+1] = trans$min[x+1]
      A[,1] = A[,1]/sum(A[,1])
      A[,length(t_star)] = A[,length(t_star)]/sum(A[,length(t_star)])
      
      # Compute eigenvalues (d) and eigenvectors (V)
      eig = eigen(A)
      d = eig$values # eigenvalues
      V = eig$vectors # eigenvectors
      
      # compute P in one step (see formula (D.39))
      P = tryCatch(
        # try inverting V first by using solve
        { Vinv = solve(V)
        if (is.complex(Vinv)) stop("Complex eigenvalue decomposition")
        V %*% diag(exp(lambda * (d-1))) %*% Vinv},
        
        #if an error occurs:
        error = function(e){
          P = diag(nrow(A)) * exp(-lambda)
          A_power = diag(nrow(A))
          count = 1
          finished = FALSE
          while (!finished & count < 10000){
            A_power = A_power %*% A
            Pold = P
            P = P + A_power * exp(-lambda + count * log(lambda) - lfactorial(count))
            
            finished = (max(abs(P-Pold)) < 1e-16)
            count = count + 1
          }
          return(P)
        })
      colnames(P) = paste("t =",t_star) 
      rownames(P) = paste("t* =",t_star)
      
      # using the proposed prior by Scholtus & de Waal (2020-2022)
      rev_probabilities1 = multiply_P_prob(n, t_obs, r, s, P, t_star, prior = "proposed")
      colnames(rev_probabilities1) = paste("t* =",t_star) 
      rownames(rev_probabilities1) = paste("t =",t_star)
      # Save the probabilities corresponding to the specific observed cell
      probabilities_prop[1:length(t_star),1,ii,jj] = t_star
      probabilities_prop[1:length(t_star),2,ii,jj] = rev_probabilities1[,which(t_star == t_obs,arr.ind=TRUE)]
      
      # using alternative prior 1 (all equal probabilities)
      rev_probabilities2 =  multiply_P_prob(n, t_obs, r, s, P, t_star, prior = "alternative1")
      # Save the probabilities corresponding to the specific observed cell
      probabilities_alt1[1:length(t_star),1,ii,jj] = t_star
      probabilities_alt1[1:length(t_star),2,ii,jj] = rev_probabilities2[,which(t_star == t_obs,arr.ind=TRUE)]
      
      # using alternative prior 2 by de Waal (2023)
      rev_probabilities3 =  multiply_P_prob(n, t_obs, r, s, P, t_star, prior = "alternative2")
      # Save the probabilities corresponding to the specific observed cell
      probabilities_alt2[1:length(t_star),1,ii,jj] = t_star
      probabilities_alt2[1:length(t_star),2,ii,jj] = rev_probabilities3[,which(t_star == t_obs,arr.ind=TRUE)]
    } 
  }
  return(list("proposed" = probabilities_prop,
              "alternative1" = probabilities_alt1,
              "alternative2" = probabilities_alt2))
}

## Reverse probabilities using specific prior distribution
multiply_P_prob = function(n, t_obs, r, s, P, t_star, prior = "proposed"){
  if (prior == "proposed"){
    priors = prior_fun(n = n, r = r, s = s, t = t_star)
  } else if (prior == "alternative1"){
    priors = rep(1/(length(t_star)), length(t_star))
  } else if (prior == "alternative2"){
    priors = alt_prior_fun(U = max(t_star), t_obs = t_obs, t = t_star)
  }
  
  prior_matrix = matrix(rep(priors, length(t_star)), byrow = TRUE, nrow = length(t_star), ncol = length(t_star))
  denominators = rowSums(prior_matrix * P)
  rev_probabilities = t(P * prior_matrix / denominators)
  colnames(rev_probabilities) = paste("t* =", t_star) 
  rownames(rev_probabilities) = paste("t =", t_star)
  return(rev_probabilities)
}

## Proposed prior probabilities Pr_0(t_jk) by Scholtus & de Waal (2020-2022) (see formula (19))
## Uses Stirlings approximation (lfactorial) because of big values of n
prior_fun = function(n, r, s, t){
  prior = choose(r, t) * choose(s, t) * exp((lfactorial(t)+lfactorial(n-r)+lfactorial(n-s)-lfactorial(n-r-s+t))-lfactorial(n))
  prior/sum(prior) # to make sure the prior probabilities sum up to 1
}

## Alternative prior probabilities Pr_0(t_jk) by de Waal(2023) (see formula (21))
alt_prior_fun = function(U, t_obs, t){
  p_hat = t_obs/U
  prior = choose(U, t)* p_hat^t * (1-p_hat)^(U-t)
  prior/sum(prior) # to make sure the prior probabilities sum up to 1
}


################################################################
### Expected values for the contingency table from section 3.3.4

expected_values_fun = function(n_cat1, n_cat2, probabilities){
  out = matrix(NA, nrow = n_cat1, ncol = n_cat2)
  for (i in 1:nrow(out)){
    for (j in 1:ncol(out)){
      out[i,j] = sum(na.omit(probabilities[,1,i,j])*na.omit(probabilities[,2,i,j]))
    }
  }
  return(out)
}
