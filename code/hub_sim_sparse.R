options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(toString(args[1]))
P <- as.numeric(toString(args[2]))
prop_missing <- as.numeric(toString(args[3]))

N=500
P=100
prop_missing = 0

library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(CVXR)
library(Robocov)
library(CorShrink)
library(emdbook)

hub_sim = function(n, p, block){
  mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
  Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
  corSigma <- cov2cor(Sigma)
  data <- MASS::mvrnorm(n,rep(0,p),corSigma)
  ll = list("dat" = data, "cor" = corSigma)
  return(ll)
}
corSigma = as.matrix(hub_sim(N, P, 10)$cor)


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
  data = hub_sim(N, P, 10)$dat

  #######################   Turn some of the entries to NA   ###################################

  data_missing = t(apply(data, 1, function(x){
    if(prop_missing > 0){
      rand = sample(1:length(x), floor(prop_missing*length(x)), replace = F)
      y = x
      y[rand] = NA
      return(y)
    }else{
      return(x)
    }}))



  standard_cor = cor(data_missing, use = "pairwise.complete.obs")

  strimmer_sample <- corpcor::cor.shrink(data)
  strimmer_cor = as.matrix(cov2cor(strimmer_sample[1:P,1:P]))

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
  pdsoft_cor = as.matrix(cov2cor(pdsoft_sample$sigma))

  corrplot::corrplot(as.matrix(standard_cor), diag = FALSE,
                     col = colorRampPalette(c("blue", "white", "red"))(200),
                     tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                     rect.col = "white",na.label.col = "white",
                     method = "color", type = "upper")
  cov_sample_ML <-  CorShrinkData(data_missing, sd_boot = FALSE,
                                  ash.control = list())
  corshrink_cor = cov2cor(cov_sample_ML$cor)
  corrplot::corrplot(as.matrix(corshrink_cor), diag = FALSE,
                     col = colorRampPalette(c("blue", "white", "red"))(200),
                     tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                     rect.col = "white",na.label.col = "white",
                     method = "color", type = "upper")

  robocov_box_cor = Robocov_box(data_with_missing = data_missing)
  corrplot::corrplot(as.matrix(robocov_box_cor), diag = FALSE,
                     col = colorRampPalette(c("blue", "white", "red"))(200),
                     tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                     rect.col = "white",na.label.col = "white",
                     method = "color", type = "upper")
  alpha_vec = c(1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  nlik_list = rep(0, length(alpha_vec))
  for(m in 1:length(alpha_vec)){
    temp_cor = Robocov_box_slack(data_missing, alpha=alpha_vec[m])
    nlik_list[m] = nloglik(data_missing, as.matrix(nearPD(temp_cor)$mat))
    cat("Passed alpha value: ", alpha_vec[m], "\n")
  }
  final_alpha = alpha_vec[which.min(nlik_list)]
  robocov_box_slack_cor =  Robocov_box_slack(data_missing, alpha=final_alpha)
  corrplot::corrplot(as.matrix(robocov_box_slack_cor), diag = FALSE,
                     col = colorRampPalette(c("blue", "white", "red"))(200),
                     tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                     rect.col = "white",na.label.col = "white",
                     method = "color", type = "upper")

  alpha_vec = c(1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  nlik_list = rep(0, length(alpha_vec))
  for(m in 1:length(alpha_vec)){
    temp_cor = Robocov_local(data_missing, alpha=alpha_vec[m])
    nlik_list[m] = nloglik(data_missing, as.matrix(nearPD(temp_cor)$mat))
    cat("Passed alpha value: ", alpha_vec[m], "\n")
  }
  final_alpha = alpha_vec[which.min(nlik_list)]
  robocov_local_cor =  Robocov_local(data_missing, alpha=5)
  corrplot::corrplot(as.matrix(robocov_local_cor), diag = FALSE,
                     col = colorRampPalette(c("blue", "white", "red"))(200),
                     tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                     rect.col = "white",na.label.col = "white",
                     method = "color", type = "upper")


  df2 = c(norm(robocov_box_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(robocov_local_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(robocov_box_slack_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(corshrink_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(standard_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(strimmer_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"),
          norm(pdsoft_cor[which(corSigma==0)] - corSigma[which(corSigma==0)], type = "2"))

  df3 = c(max(abs(robocov_box_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(robocov_local_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(robocov_box_slack_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(corshrink_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(standard_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(strimmer_cor[which(corSigma==0)] - corSigma[which(corSigma==0)])),
          max(abs(pdsoft_cor[which(corSigma==0)] - corSigma[which(corSigma==0)]))
          )

  df4 = c(length(which(robocov_box_cor[which(corSigma==0)] > 0.1))/length(robocov_box_cor[which(corSigma==0)]),
          length(which(robocov_local_cor[which(corSigma==0)] > 0.1))/length(robocov_local_cor[which(corSigma==0)]),
          length(which(robocov_box_slack_cor[which(corSigma==0)] > 0.1))/length(robocov_box_slack_cor[which(corSigma==0)]),
          length(which(corshrink_cor[which(corSigma==0)] > 0.1))/length(corshrink_cor[which(corSigma==0)]),
          length(which(standard_cor[which(corSigma==0)] > 0.1))/length(standard_cor[which(corSigma==0)]),
          length(which(strimmer_cor[which(corSigma==0)] > 0.1))/length(standard_cor[which(corSigma==0)]),
          length(which(pdsoft_cor[which(corSigma==0)] > 0.1))/length(standard_cor[which(corSigma==0)])
          )

  df6 = c(length(which(robocov_box_cor[which(corSigma!=0)] > 0.1))/length(robocov_box_cor[which(corSigma!=0)]),
          length(which(robocov_local_cor[which(corSigma!=0)] > 0.1))/length(robocov_local_cor[which(corSigma!=0)]),
          length(which(robocov_box_slack_cor[which(corSigma!=0)] > 0.1))/length(robocov_box_slack_cor[which(corSigma!=0)]),
          length(which(corshrink_cor[which(corSigma!=0)] > 0.1))/length(corshrink_cor[which(corSigma!=0)]),
          length(which(standard_cor[which(corSigma!=0)] > 0.1))/length(standard_cor[which(corSigma!=0)]),
          length(which(strimmer_cor[which(corSigma!=0)] > 0.1))/length(standard_cor[which(corSigma!=0)]),
          length(which(pdsoft_cor[which(corSigma!=0)] > 0.1))/length(standard_cor[which(corSigma!=0)])
  )

  df = rbind(df2, df3, 1-df4, 1-df5, df6)

  colnames(df) = c("Robocov_box",
                "Robocov_local",
                "Robocov_box_slack_cor",
                "CorShrink",
                "standard",
                "strimmer",
                "pdsoft")
  rownames(df) = c("2-norm-0", "Inf-norm-0", "specificity-0.1", "specificity-0.05", "sensitivity-0.1")

  ll[[nsim]] = df
  cat("We are at simulation trial:", nsim, "\n")
}

xtable2(t((df[c(1, 3, 5),c(1,3,2,4,6,7, 5)])))

save(df, file = paste0("/n/groups/price/kushal/Robocov/output_sparse/robocov_sim_hub_n_",
                       N, "_p_", P, "_prop_", prop_missing, ".rda"))



