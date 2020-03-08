options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(toString(args[1]))
P <- as.numeric(toString(args[2]))
count_missing <- as.numeric(toString(args[3]))

N=50
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

DM_toeplitz = function(n,P){
  library("MASS")
  index1=sort(sample(seq(1:n),(n/2)))
  index2=seq(1:n)[-index1]
  Sigmatp=function(P){
    a=array(0,dim=c(P,P))
    for(i in 1:P){
      for(j in 1:P){
        a[i,j]=max(1-0.1*(abs(i-j)),0)
      }
    }
    return(a)
  }
  Sigma = Sigmatp(P)
  data = mvrnorm(n,rep(0,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = solve(Sigma)
  return(list(Xtrain = Xtrain, Xtest = Xtest, Sigma = Sigma))
}


toeplitz_sim = function(N, P){
  ll <- DM_toeplitz(n=N, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)
  ll = list("dat" = data, "cor" = corSigma)
  return(ll)
}

corSigma = as.matrix(toeplitz_sim(N, P)$cor)
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
  data = toeplitz_sim(N, P)$dat

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

  cormat = cor(data_missing, use = "pairwise.complete.obs")

  gg = glasso::glasso(cormat, rho=1)
  pglasso = -cov2cor(gg$wi)
  diag(pglasso) = 1
  image(pglasso)


  df1 = c(norm(robocov_precision - pcorSigma, type = "2"),
          norm(pglasso - pcorSigma, type = "2"))

  df2 = c(norm(robocov_precision[which(pcorSigma < 0.01)] - pcorSigma[which(pcorSigma < 0.01)], type = "2"),
          norm(pglasso[which(pcorSigma< 0.01)] - pcorSigma[which(pcorSigma < 0.01)], type = "2"))

  df4 = c(length(which(robocov_precision[which(pcorSigma< 0.01)] > 0.1))/length(robocov_precision[which(pcorSigma< 0.01)]),
          length(which(pglasso[which(pcorSigma< 0.01)] > 0.1))/length(pglasso[which(pcorSigma< 0.01)]))

  df6 = c(length(which(robocov_precision[which(pcorSigma!=0)] > 0.1))/length(robocov_precision[which(pcorSigma!=0)]),
          length(which(pglasso[which(pcorSigma!=0)] > 0.1))/length(pglasso[which(pcorSigma!=0)]))

  df = rbind(df1, df2, 1-df4, df6)

  colnames(df) = c("Robocov_pcor",
                   "GLASSO")
  rownames(df) = c("2-norm", "2-norm0", "specificity-0.1", "sensitivity-0.1")

  ll[[nsim]] = df
  cat("We are at simulation trial:", nsim, "\n")
}


save(ll, file = paste0("/n/groups/price/kushal/Robocov/output_Nov26_prec/robocov_prec_sim_toeplitz_n_",
                       N, "_p_", P, "_prop_", count_missing, ".rda"))

