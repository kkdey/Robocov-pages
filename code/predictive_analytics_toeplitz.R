options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(toString(args[1]))
P <- as.numeric(toString(args[2]))
count_missing <- as.numeric(toString(args[3]))
prop_missing = count_missing/100

library(Matrix)
library(CVXR)
library(Robocov)
library(CorShrink)

N=500
P=50
count_missing=0
prop_missing = count_missing/100

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

data = toeplitz_sim(N, P)$dat

#######################   Turn some of the entries to NA (missing)   ###################################

data_missing = t(apply(data, 1, function(x){
  if(prop_missing > 0){
    rand = sample(1:length(x), floor(prop_missing*length(x)), replace = F)
    y = x
    y[rand] = NA
    return(y)
  }else{
    return(x)
  }}))

gene_data = data_missing
measure = c()
measure2 = c()

for(num_iter in 1:5){
  train_sample_id = sample(1:nrow(gene_data), floor(nrow(gene_data)/2), replace = FALSE)
  predict_sample_id = setdiff(1:nrow(gene_data), train_sample_id)
  train_datamat = gene_data[train_sample_id,]

  empirical_cor = cor(train_datamat, method = "pearson", use = "pairwise.complete.obs")
  empirical_cor[is.na(empirical_cor)] = 0

  corshrink_out = CorShrinkData(train_datamat, sd_boot = FALSE, image = "null",
                                image.control = list(tl.cex = 0.2))
  corshrink_cor = corshrink_out$cor

  robocov_box_cor = Robocov_box(data_with_missing = train_datamat, loss = "lasso")

  predict_datamat = gene_data[predict_sample_id,]
  cormat2 = cor(predict_datamat, method = "pearson", use = "pairwise.complete.obs")
  cormat2[is.na(cormat2)] = 0
  measure = rbind(measure, c(mean(abs(cormat2 - empirical_cor)),
                             mean(abs(cormat2 - corshrink_cor)),
                             mean(abs(cormat2 - robocov_box_cor)))
  )
  measure2 = rbind(measure2, c(sqrt(mean((cormat2 - empirical_cor)^2)),
                               sqrt(mean((cormat2 - corshrink_cor)^2)),
                               sqrt(mean((cormat2 - robocov_box_cor)^2))))
  cat("We finished iteration:", num_iter)
}

df = cbind(apply(measure, 2, mean),
           apply(measure, 2, sd),
           apply(measure2, 2, mean),
           apply(measure2, 2, sd))
colnames(df) = c("L1", "se(L1)", "L2", "se(L2")
rownames(df) = c("Sample-Est", "CorShrink", "Robocov")

ll = list("L1" = measure, "L2" = measure2)
save(ll, file = paste0("/n/groups/price/kushal/Robocov/output/Predictive/toeplitz_N_",N,
                       "_P_", P, "_pi_",count_missing, ".rda"))
