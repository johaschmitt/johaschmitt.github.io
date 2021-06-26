# Simulations for the Alternative score

library(data.table)
library(ggplot2)
library(emulator)
library(gridExtra)
library(optimr)
library(EstimationTools)
library(kernlab)
library(igraph)
library(rje)
library(purrr)
library(beepr)


#setwd("/home/johnny/MT_ROOT/WD")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#------------------- GROUND TRUTH SECTION -------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------- Ground Truth Setup ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# number of variables
n_var <- 3
all_var <- rep(1:n_var)
all_var_names <- c("x_1", "x_2", "x_3")

# choose ground truth causal graph
g_true <- graph(edges = c(1,3, 2,3), n = n_var, directed = T)

# plot graph
plot(g_true)

# adjaceny matrix, rows: outgoing edeges, columns: incoming edges ... of that node
adjm_true <- as.matrix(g_true[])

# choose ground truth SCM
scm_true <- list(X1 = function(n) rnorm(n, 0, 1),
                 X2 = function(n) rnorm(n, 1, 1),
                 X3 = function(x1, x2, n) 2*tanh(x1) + 4*tanh(x2) + rnorm(n, 0, 0.1))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#---------- Available functions ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' generate n samples (x_1. ...) from the joint underlying distribution
#' 
#' @param input SCM (either true observational or some interventional)
#' @param n sample size
#' 
#' @return  the samples (x_1. ...) as data table  
joint_sample <- function(scm, n){
  
  # evaluation of variables given the input SCM 
  x1 = scm$X1(n)
  x2 = scm$X2(n)
  
  x <- list(x_1 = x1,
            x_2 = x2,
            x_3 = scm$X3(x1, x2, n))
  
  return(as.data.table(x))
}


#' generate n samples (default 1) from the interventional distribution P(.|do(X_i = x))
#' 
#' @param scm_intv the SCM that's scheduled for intervention
#' @param i_do index of variable that is being intervened on
#' @param x_do intervention value on X_i
#' @param n number of returned samples from the interventional SCM
#' 
#' @return the samples (x_1. ...) as data table  
Do_sample <- function(scm_intv, i_do, x_do, n = 1){
  
  # set the intervention variable to the intervention value
  body(scm_intv[[i_do]]) <- substitute(x_do)
  
  # get new samples
  smp <- joint_sample(scm = scm_intv, n = n)
  
  # tag the samples
  smp <- smp[,doX := i_do,]
  
  return(smp)
}




#' generate n(=1) sample(s) from the interventional distribution P(.|do(X_i = xi, X_j = xj))
#' 
#' @param scm_intv the SCM that's scheduled for intervention
#' @param i_do index of 1st variable that is being intervened on
#' @param xi_do intervention value on X_i
#' @param j_do index of 2nd variable that is being intervened on
#' @param xj_do intervention value on X_j
#' @param n number of returned samples from the interventional SCM
#' 
#' @return the samples (x_1. ...) as data table  
Do2_sample <- function(scm_intv, i_do, xi_do, j_do, xj_do, n = 1){
  
  # set the intervention variables to the intervention values
  body(scm_intv[[i_do]]) <- substitute(xi_do)
  body(scm_intv[[j_do]]) <- substitute(xj_do)
  
  # get new samples
  smp <- joint_sample(scm = scm_intv, n = n)
  
  # tag the samples
  smp <- smp[,doX := NA,][,doXi := i_do,][,doXj := j_do,]
  
  return(smp)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------------------- DISCOVERY --------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#---------- Useful general functions ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' elementary symmetric polynomials
#' 
#' @param z vector of inputs
#' @param o_p order of polynomial
#' 
#' @return value of elementary symmetric polynomial
es_poly <- function(z, o_p){
  
  # input dimension
  d = length(z)
  
  # compute the powerset of z up to order o_p
  pS <- powerSet(z, m = o_p)
  
  # filter out the subsets not of that order
  S <- Filter(Negate(is.null), lapply(pS, function(i) if(length(i) == o_p) i))
  
  # return the sum of products of the respective interactions
  return(sum(unlist(lapply(S, function(i) prod(i)))))
}


#' customized base kernel: squared exponential (+ linear, turned off) 
#' 
#' @param v outer coefficient of exponential part
#' @param w outer coefficient of exponential part
#' @param a coefficient of linear part
#' @param x first input
#' @param y second input
#' 
#' @return kernel value
k_base <- function(v, w, a, x, y) return(v*exp(-w*((abs(x-y))^2)/2))

# define class from kernlab package 
class(k_base) <- "kernel"


#' generic kernel/covariance function with customized kernel with noise
#' 
#' @param lhp list of log parameters v, w, a, k_add_o, vv
#' @param X_l data in design matrix form (columns = features resp. dimensions)
#' @param X_r secondary data for K(X_l,X_r) in design matrix form
#' @param vv_incl Boolean for inclusion of variance in resulting kernel matrix, only sensible for X_l == X_r
#' 
#' @return the kernel/covariance matrix
cov_cust <- function(lhp, X_l, X_r, vv_incl = F){
  
  # number of input dimensions
  D = ncol(X_l)
  
  # additive kernel
  ker_cust <- function(x,y){
    
    # base kernel values (see Duvenaud)
    z <- k_base(v = exp(lhp$ex.v), w = exp(lhp$ex.w), a = exp(lhp$lin.a), x = x, y = y)
    
    # full additive kernel w/ resp. weights up to full order D 
    
    # (not recursive)
    k_add <- sum(exp(lhp$k_add_w) * sapply(1:D, function(o) es_poly(z, o)))
    
    # switch to  alternative recursive, if D >= 5
    # list of (o-1)th order kernels, 1st elemnt dummy for the sum in the for loop
    #k_add_l <- list(1)
    
    # recursive formula "Newton-Girard" to efficiently compute the elementary symmetric polynomials
    #for(o in 1:D){
    #      k_add_l[[o + 1]] <- 1/o * sum(sapply(1:o, function(k) (-1)^(k-1) * k_add_l[[(o-k) + 1]] * sum(z^k))) 
    #}
    
    #k_add <- sum(exp(lhp$k_add_w) * unlist(k_add_l)[2:(D+1)])
    
    return(k_add)
  } 
  
  if(vv_incl){
    
    # case of X_l == X_r, number of observations: 
    n <- nrow(X_l)
    return(kernelMatrix(kernel = ker_cust, x = as.matrix(X_l), y = as.matrix(X_r)) + diag(exp(lhp$vv), n))
  } 
  
  else return(kernelMatrix(kernel = ker_cust, x = as.matrix(X_l), y = as.matrix(X_r)))
}


#' derivative of generic kernel/covariance function with customized kernel with noise
#' 
#' @param lhp list of log parameters v, w, a, vv
#' @param X data 
#' 
#' @return "gradient" of the kernel/covariance matrix (order of hyperparameters in lhp maintained)
#'         i.e. list of derivative matrices
D_cov_cust <- function(lhp, X){
  
  # number of observations and dimensions
  n = nrow(X) 
  D = ncol(X)
  
  
  # derivatives ...
  
  # ... wrt ex.v
  ker_cust_dv <- function(x, y){
    
    # base kernel values (see Duvenaud)
    z <- k_base(v = exp(lhp$ex.v), w = exp(lhp$ex.w), a = exp(lhp$lin.a), x = x, y = y)
    
    # apply chain rule: df/dz * dz/dv, where z = c(z_1, z_2, ..., z_D) and f = k_add_o
    # dk_add_o/dz: D 1xD matrices (see Duvenaud)
    dk_add_o_diff_dz <- t(sapply(1:D, function(o) sapply(1:D, function(d) es_poly(z[-d],o-1)))) 
    
    # dz/dv
    dz_diff_dv <- as.matrix(k_base(v = 1, w = exp(lhp$ex.w), a = 0, x = x, y = y))
    
    # column vector partial order derivatives w/o additive weights
    dk_add_o_diff_dv <- dk_add_o_diff_dz %*% dz_diff_dv
    
    # entire partial derivative
    return(sum(exp(lhp$k_add_w) * dk_add_o_diff_dv))
  }
  
  cov_dv <- kernelMatrix(kernel = ker_cust_dv, x = as.matrix(X), y = as.matrix(X))
  
  
  # ... wrt ex.w
  ker_cust_dw <- function(x, y){
    
    # base kernel values (see Duvenaud)
    z <- k_base(v = exp(lhp$ex.v) * (-((abs(x-y))^2)/2), w = exp(lhp$ex.w), a = 0, x = x, y = y)
    
    # apply chain rule: df/dz * dz/dv, where z = c(z_1, z_2, ..., z_D) and f = k_add_o
    # dk_add_o/dz: D 1xD matrices (see Duvenaud)
    dk_add_o_diff_dz <- t(sapply(1:D, function(o) sapply(1:D, function(d) es_poly(z[-d],o-1)))) 
    
    # dz/dv
    dz_diff_dv <- as.matrix(k_base(v = 1, w = exp(lhp$ex.w), a = 0, x = x, y = y))
    
    # column vector partial order derivatives w/o additive weights
    dk_add_o_diff_dv <- dk_add_o_diff_dz %*% dz_diff_dv
    
    # entire partial derivative
    return(sum(exp(lhp$k_add_w) * dk_add_o_diff_dv))
  }
  
  cov_dw <- kernelMatrix(kernel = ker_cust_dw, x = as.matrix(X), y = as.matrix(X))
  
  
  # ... wrt lin.a
  ker_cust_da <- function(x, y){
    
    # base kernel values (see Duvenaud)
    z <- k_base(v = 0, w = 0, a = 1, x = x, y = y)
    
    # apply chain rule: df/dz * dz/dv, where z = c(z_1, z_2, ..., z_D) and f = k_add_o
    # dk_add_o/dz: D 1xD matrices (see Duvenaud)
    dk_add_o_diff_dz <- t(sapply(1:D, function(o) sapply(1:D, function(d) es_poly(z[-d],o-1)))) 
    
    # dz/dv
    dz_diff_dv <- as.matrix(k_base(v = 1, w = exp(lhp$ex.w), a = 0, x = x, y = y))
    
    # column vector partial order derivatives w/o additive weights
    dk_add_o_diff_dv <- dk_add_o_diff_dz %*% dz_diff_dv
    
    # entire partial derivative
    return(sum(exp(lhp$k_add_w) * dk_add_o_diff_dv))
  }
  
  cov_da <- kernelMatrix(kernel = ker_cust_da, x = as.matrix(X), y = as.matrix(X))
  
  
  # ... wrt k_add_w[]: list of ker_cust functions and then list of covariance matrices
  ker_cust_dadd_w <- lapply(1:D, function(d) return(function(x, y){
    
    # base kernel values (see Duvenaud)
    z <- k_base(v = exp(lhp$ex.v), w = exp(lhp$ex.w), a = exp(lhp$lin.a), x = x, y = y)
    
    # derivative of full additive kernel w/ respect to weights, not recursive
    dk_add <- es_poly(z, d)
    
    return(dk_add)
  }))
  
  cov_dadd_w <- lapply(1:D, function(d) kernelMatrix(kernel = ker_cust_dadd_w[[d]], x = as.matrix(X), y = as.matrix(X)))
  
  
  # ... wrt vv
  cov_dvv <- diag(n)
  
  # list w/ all cov matrices
  cov_list <- list(cov_dv, cov_dw, cov_da)
  for(d in 1:D) cov_list[[d+3]] <- cov_dadd_w[[d]]
  cov_list[[D+4]] <- cov_dvv
  
  return(cov_list)
}


#' source node estimation w/ ML and simple univariate Gaussians
#' 
#' @param xtrain univariate samples (expected still in data.table format)
#' 
#' @return list of estimated mean and variance
est_source <- function(xtrain){
  
  # simple ML scheme for Gaussian
  f_Xloglike <- maxlogL(x = as.double(unlist(xtrain)), dist = 'dnorm', start=c(0, 1), lower=c(-4, 0.01), upper=c(4, 4))
  
  # estimates for mean and variance
  return(list(X_mea = f_Xloglike$fit$par[1], X_var = f_Xloglike$fit$par[2]^2))
}


#' dependant node estimation: GPR w/ ML-II scheme for hyperparameter optimization
#' 
#' @param ytrain univariate data of the dependant node
#' @param xtrain possibly multivariate data of parents
#' @param dep_par list of parameters for the regression situation (level of dep_param[[1]][[1]])
#' 
#' @return list of estimates in standard form: 
est_dependant <- function(ytrain, xtrain, dep_par){
  
  # ensure canonical form: column vectors for the various dimensions  
  ytrain <- t(t(ytrain))
  xtrain <- t(t(xtrain))
  
  # number of observations and dimensions
  n = nrow(xtrain)
  D = ncol(xtrain)
  
  # get initial hyperparameter
  log_hp <- dep_par$para
  
  
  #------ functions 
  
  #' (marginal) GP loglikelohood 
  #' 
  #' @param lpara log hyperparameter 
  #' 
  #' @return likelihood value
  minus_maloglike <- function(lpara){
    
    lhp <- list("ex.v" = lpara[1], "ex.w" = lpara[2], "lin.a" = lpara[3], "k_add_w" = lpara[4:(D+3)], "vv" = lpara[(D+4)])
    Cv <- cov_cust(lhp, X_l = xtrain, X_r = xtrain, vv_incl = T)
    
    return(0.5*quad.form(solve(Cv), unlist(ytrain)) + 0.5*log(det(Cv)) + 0.5*n*log(2*pi))
  }
  
  #' gradients of those  minus GP loglikelihoods 
  #' 
  #' @param lpara log hyperparameter 
  #' 
  #' @return likelihood value
  D_minus_maloglike <- function(lpara){
    
    lhp <- list("ex.v" = lpara[1], "ex.w" = lpara[2], "lin.a" = lpara[3], "k_add_w" = lpara[4:(D+3)], "vv" = lpara[(D+4)])
    Cv <- cov_cust(lhp, X_l = xtrain, X_r = xtrain, vv_incl = T)
    Cv_inv <- solve(Cv)
    D_Cv <- D_cov_cust(lhp, X = xtrain)
    
    alph <- Cv_inv %*% as.matrix(ytrain)
    
    # return derivative of minus loglikelihood 
    return(unlist(lapply( D_Cv, function(i) return(- 0.5 * tr(((alph %*% t(alph)) - Cv_inv) %*% i)))))
  }
  
  
  #------ main
  
  lhp_init <- list("ex.v" = log(5), "ex.w" = log(10), "lin.a" = log(10), 
                   "k_add_w" = log(rep(1/D, D)), "vv" = log(1))
  
  o_res <- optimr(fn = minus_maloglike, gr = D_minus_maloglike, par = as.double(unlist(lhp_init)))
  
  
  dep_par$para = list("ex.v" = o_res$par[1], "ex.w" = o_res$par[2], "lin.a" = o_res$par[3], 
                      "k_add_w" = o_res$par[4:(D+3)], "vv" = o_res$par[(D+4)])
  
  dep_par$conv = o_res$convergence
  
  # calculate covariance matrices with optimized parameters
  dep_par$Cov <- cov_cust(dep_par$para, X_l = xtrain, X_r = xtrain, vv_incl = T)
  
  return(dep_par)
}


#' pseudo source node estimation: GPR w/ ML-II scheme for hyperparameter optimization
#' 
#' @param data data table with observational data and fixed value intervention on all other d-1 variables
#' @param psou_par list of parameters for the regression situation (level of psou_param[[1]][[1]])
#' 
#' @return list of estimates in standard form: 
est_pseudo_source <- function(psou_data, psou_par){
  
  # select correct (regular + additional) observational sample data, restricted to max. 3 variable case
  if(psou_par$sou == 1){
    psou_data <- psou_data[,-c(2,3)]
    sou_data_obs_add <- data_obs_add[,1]
  } 
  if(psou_par$sou == 2){
    psou_data <-psou_data[,-c(1,3)]
    sou_data_obs_add <- data_obs_add[,2]
  } 
  if(psou_par$sou == 3){
    psou_data <- psou_data[,-c(1,2)]
    sou_data_obs_add <- data_obs_add[,3]
  } 

  # calculate observational sample variance (corrected), also from additional observational data
  smp_var_obs = as.double(var(rbind(psou_data[doX == 0][,1],sou_data_obs_add)))
  
  # calculate intventional sample variance (corrected)
  smp_var_intv = as.double(var(psou_data[is.na(doX)][,1]))
  
  # set base noise variance equal to the noise variance around functional values in GPR settings (base level now fixed sd = 0.1, var = 0.01)
  eps_base <- 0.01
  
  # set additional noise according to linear function: 
  #     if intervetional sample variance is high (probably source with fixed interventions) -> small additional noise, low penalty
  #     else -> additional variance as high as spread of entire variable 'smp_var_obs', high penalty
  eps_add <- max(0,-smp_var_intv + smp_var_obs)
  
  # calculate total noise as sum of base and additional noise
  eps_total <- eps_add + eps_base
  
  
  
  # generate pseudo GP fitting data, first the observational data with just base noise
  data_pseudofit_obs <- list(x_1 = psou_data[doX == 0][,1],
                             x_2 = rnorm(n_obs, 0, sqrt(eps_base)))
  
  # then the for the "interventional data" randomly distributed points with full noise eps_total, change with different "support interval"
  data_pseudofit_intv <- list(x_1 = runif(((n_var-1)*n_intv), min = -2, max = 2),
                              x_2 = rnorm(((n_var-1)*n_intv), 0, sqrt(eps_total)))
  
  # entire data to do GP pseudo fit on 
  data_pseudofit <- rbind(as.data.table(data_pseudofit_obs), as.data.table(data_pseudofit_intv))
  
  #plot(data_pseudofit)
  
  # fit pseudo GP
  res_psou_par <- est_dependant(ytrain = data_pseudofit$x_2, xtrain = data_pseudofit$x_1, dep_par = psou_par)
  
  res_psou_par[["ydat"]] <- data_pseudofit$x_2
  
  return(res_psou_par)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------- Learning Setup ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# considered models
models <- list(
  # v-structures
  graph(edges = c(1,3, 2,3), n = 3, directed = T),
  graph(edges = c(1,2, 3,2), n = 3, directed = T),
  graph(edges = c(2,1, 3,1), n = 3, directed = T),
  # reversed v-structures
  graph(edges = c(1,2, 1,3), n = 3, directed = T),
  graph(edges = c(2,1, 2,3), n = 3, directed = T),
  graph(edges = c(3,1, 3,2), n = 3, directed = T),
  # hierarchical lines
  graph(edges = c(1,2, 2,3), n = 3, directed = T),
  graph(edges = c(1,3, 3,2), n = 3, directed = T),
  graph(edges = c(2,1, 1,3), n = 3, directed = T),
  graph(edges = c(2,3, 3,1), n = 3, directed = T),
  graph(edges = c(3,1, 1,2), n = 3, directed = T),
  graph(edges = c(3,2, 2,1), n = 3, directed = T),
  # v-structures plus parent dependency
  graph(edges = c(1,2, 1,3, 2,3), n = 3, directed = T),
  graph(edges = c(1,3, 1,2, 3,2), n = 3, directed = T),
  graph(edges = c(2,1, 2,3, 1,3), n = 3, directed = T),
  graph(edges = c(2,3, 2,1, 3,1), n = 3, directed = T),
  graph(edges = c(3,1, 3,2, 1,2), n = 3, directed = T),
  graph(edges = c(3,2, 3,1, 2,1), n = 3, directed = T),
  # one independent variable
  graph(edges = c(1,2), n = 3, directed = T),
  graph(edges = c(2,1), n = 3, directed = T),
  graph(edges = c(1,3), n = 3, directed = T),
  graph(edges = c(3,1), n = 3, directed = T),
  graph(edges = c(2,3), n = 3, directed = T),
  graph(edges = c(3,2), n = 3, directed = T),
  # all independent
  graph(edges = c(), n = 3, directed = T)
)

# number of models
n_models <- length(models)

# list of respective adjacency matrices
adjm <- lapply(models, function(g) as.matrix(g[]))

# distribution of considered graphs
prob_g <- rep(1/n_models, n_models)

# current sample size
n_obs <- 5 

# generate 5 samples
data <- joint_sample(scm = scm_true, n_obs)

# intervention indicator: 0 = no (obs), 1 = doX_1, 2 = doX_2, ...
data <- data[,doX := 0,]



# allow to draw 15 additional observational samples just to estimate the sample variances for pseudo nodes later on equal grounds,
# not used otherwise
data_obs_add <- joint_sample(scm = scm_true, 15)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#---------- Causal Discovery ----------  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


######## for observational and interventional data


#****************************************************** 
# number of interventions
n_intv <- 10

# total number of values that GPs are fitted on
n_total <- ((n_var-1)*n_intv) + n_obs

# Experiments: get interventional data of K_1
for(i in 1:n_intv){
  
  # choose interventions at random uniformly in [-1, 1]
  x_intv <- runif(n_var, min = -1, max = 1)
  
  # get interventional data
  for(i_v in 1:n_var) data = rbind(data, Do_sample(scm_intv = scm_true, i_do = i_v, x_do = x_intv[i_v], n = 1))
  
}



# ... and get interventional data of K_2
# generate intervention data for fixed-value-interventions
data <- data[,doXi := NA,][,doXj := NA,] 

# intervene on every variable w/ fixed intv value 0
for(not_v in 1:n_var){
  
  i_v <- setdiff(1:n_var, not_v)[1]
  j_v <- setdiff(1:n_var, not_v)[2]
  
  # sample size (n_var-1)*n_intv to match sample size for GP =  n_obs + (n_var-1)*n_intv   
  smp <- Do2_sample(scm_intv = scm_true, i_do = i_v, xi_do = 0, j_do = j_v, xj_do = 0, n = ((n_var-1)*n_intv))
  
  data = rbind(data, smp)
} 




#******************************************************   
# for source nodes: ML estimation
# get source nodes (no incoming edges means zero sum columns)
sources <- lapply(adjm, function(m) which(colSums(m) == 0))


#******************************************************          
# for non-source nodes or dependants 
# identify the dependants for each model
dependants <- lapply(1:n_models, function(i_m) setdiff(all_var, sources[[i_m]]))

# identify parents for each dependant
# result: list[models] of lists[dependant nodes] of lists[dependant node (1) or his parents (2)]
dep_parent <- lapply(1:n_models, function(i_m)  lapply(dependants[[i_m]], function(i_d)  list("dep" = i_d, "parent" = which(adjm[[i_m]][, i_d] == 1))))

# list to carry the GPR parameters to be estimated
# result: list[models] of lists[dependant nodes] of lists[dependant node (1)or the list of parameter of the additive kernel GPR (2)]



print("estimating dependant parameters...")
# INITIALIZE all dependant parameters for GPR
#' @param ex.v log linear weight for (exp)^2 part
#' @param ex.w log exp weight for (exp)^2 part
#' @param lin.a log linear weight lin part
#' @param k_add_w log vector of weights for the o-th order additive kernel respectively
#' @param vv log variance of GPR w/ normal noise 
# i_d now list in style of "dep_parent[[1]][[1]]"
dep_param <- lapply(1:n_models, function(i_m)  lapply(dep_parent[[i_m]], 
                                                      function(i_d) {
                                                        i_d[["para"]] <- list("ex.v" = log(5), "ex.w" = log(10), "lin.a" = log(10), 
                                                                              "k_add_w" = log(rep(1/length(i_d$parent), length(i_d$parent))), "vv" = log(1))
                                                        i_d[["conv"]] <- 0
                                                        i_d[["Cov"]] <- matrix()
                                                        return(i_d)
                                                      }))

# ESTIMATE all dependant parameters via the Gaussian process likelihoods
dep_param <- lapply(1:(n_models-1), function(i_m)  lapply(dep_param[[i_m]], 
  function(i_d) est_dependant(as.matrix(data[doX != i_d$dep & !is.na(doX)])[,i_d$dep], as.matrix(data[doX != i_d$dep & !is.na(doX)])[,i_d$parent], i_d)))

dep_param[[1]][[1]]



print("estimating pseudo source parameters...")
# INITIALIZE all pseudo source parameters for GPR
#' @param ex.v log linear weight for (exp)^2 part
#' @param ex.w log exp weight for (exp)^2 part
#' @param lin.a log linear weight lin part
#' @param k_add_w log vector of weights for the o-th order additive kernel respectively
#' @param vv log variance of GPR w/ normal noise 
# i_d now list in style of "dep_parent[[1]][[1]]"
psou_param <- lapply(1:n_models, function(i_m)  lapply(sources[[i_m]], 
                                                      function(i_s) {
                                                        l <- list()
                                                        l[["sou"]] <- i_s
                                                        l[["para"]] <- list("ex.v" = log(5), "ex.w" = log(10), "lin.a" = log(10), 
                                                                            "k_add_w" = log(1), "vv" = log(1))
                                                        l[["conv"]] <- 0
                                                        l[["Cov"]] <- matrix()
                                                        l[["ydat"]] <- c()
                                                        return(l)
                                                      }))

# ESTIMATE all pseudo source parameter via the Gaussian process likelihoods
psou_param <- lapply(1:n_models, function(i_m)  lapply(psou_param[[i_m]], 
  function(i_s) est_pseudo_source(data[(doXi != i_s$sou & doXj != i_s$sou & is.na(doX)) | doX == 0], i_s)))




# log likelihoods using Markov factorization
# ... for sources
maloglike_psources <- lapply(1:n_models, function(i_m) sum(sapply(psou_param[[i_m]], function(i_s){
  return(-0.5*quad.form(solve(i_s$Cov), t(t(i_s$ydat))) - 0.5*log(det(i_s$Cov)) - 0.5*n_total*log(2*pi))
})))

# .. for dependants, fix n...
maloglike_dependants <- lapply(1:(n_models-1), function(i_m) sum(sapply(dep_param[[i_m]], function(i_d){
  return(-0.5*quad.form(solve(i_d$Cov), as.matrix(as.matrix(data[doX != i_d$dep & !is.na(doX)])[,i_d$dep])) - 0.5*log(det(i_d$Cov)) - 0.5*n_total*log(2*pi))
})))

# nothing for last model of independent variables
maloglike_dependants[[25]] <- 0



# ... total 
totallike <- exp(unlist(maloglike_psources) + unlist(maloglike_dependants))

# nornalization constant
norm_const <- sum(prob_g*totallike)

# probablities for the individual graphs 
prob_g = prob_g*totallike/norm_const  


beep(sound = 5)


