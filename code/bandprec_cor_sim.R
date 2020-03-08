options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(toString(args[1]))
P <- as.numeric(toString(args[2]))
count_missing <- as.numeric(toString(args[3]))

N=500
P=50
count_missing=25
prop_missing = count_missing/100

library(Matrix)
library(CVXR)
library(Robocov)
library(CorShrink)
library(emdbook)
library(MASS)
library(methods)

banded_prec_sim = function(N, P){
  diags <- list()
  diags[[1]] <- rep(1, P)
  diags[[2]] <- rep(-0.5, P)
  Kinv <- bandSparse(P, k = -(0:1), diag = diags, symm = TRUE)
  K <- solve(Kinv)
  corSigma <- cov2cor(K)
  data <- MASS::mvrnorm(N,rep(0,P),corSigma)
  ll = list("dat" = data, "cor" = corSigma)
  return(ll)
}


corSigma = as.matrix(banded_prec_sim(N, P)$cor)
pcorSigma = -cov2cor(solve(corSigma))
diag(pcorSigma) = 1

nloglik = function(data, cormat){
  llik = 0
  for(m in 1:nrow(data)){
    idx = which(!is.na(data[m,]))
    if(length(idx) > 2){
      llik = llik + emdbook::dmvnorm(data[m, idx], rep(0, length(idx)), cormat[idx, idx], log = T)
    }
  }
  return(-llik)
}

angle_norm = function(S, Sigma){
  dist = 1 - (tr(as.matrix(cov2cor(S)%*%cov2cor(Sigma))))/(norm(cov2cor(S), type = "F")* norm(Sigma, type = "F"))
  return(dist)
}


NUM_SIM=1
ll <- vector(mode="list", length=NUM_SIM)

for(nsim in 1:NUM_SIM){
  data = banded_prec_sim(N, P)$dat

  data_missing = t(apply(data, 1, function(x){
    if(prop_missing > 0){
      rand = sample(1:length(x), floor(prop_missing*length(x)), replace = F)
      y = x
      y[rand] = NA
      return(y)
    }else{
      return(x)
    }}))

  robocov_precision = Robocov_box(data_with_missing = data_missing, loss = "lasso")
  image(robocov_precision)

  gg = glasso::glasso(cor(data_missing, use = "pairwise.complete.obs"), rho=1)
  pglasso = -cov2cor(gg$wi)
  diag(pglasso) = 1
  image(pglasso)


  df1 = c(norm(robocov_precision - pcorSigma, type = "2"),
          norm(pglasso - pcorSigma, type = "2"))

  df2 = c(norm(robocov_precision[which(pcorSigma==0)] - pcorSigma[which(pcorSigma==0)], type = "2"),
          norm(pglasso[which(pcorSigma==0)] - pcorSigma[which(pcorSigma==0)], type = "2"))

  df4 = c(length(which(robocov_precision[which(pcorSigma==0)] > 0.1))/length(robocov_precision[which(pcorSigma==0)]),
          length(which(pglasso[which(pcorSigma==0)] > 0.1))/length(pglasso[which(pcorSigma==0)]))

  df6 = c(length(which(robocov_precision[which(pcorSigma!=0)] > 0.1))/length(robocov_precision[which(pcorSigma!=0)]),
          length(which(pglasso[which(pcorSigma!=0)] > 0.1))/length(pglasso[which(pcorSigma!=0)]))

  df = rbind(df1, df2, 1-df4, df6)

  colnames(df) = c("Robocov_pcor",
                   "GLASSO")
  rownames(df) = c("2-norm", "2-norm0", "specificity-0.1", "sensitivity-0.1")

  ll[[nsim]] = df
  cat("We are at simulation trial:", nsim, "\n")
}




save(ll, file = paste0("/n/groups/price/kushal/Robocov/output_Nov26_prec/robocov_prec_sim_bandprec_n_",
                       N, "_p_", P, "_prop_", count_missing, ".rda"))

