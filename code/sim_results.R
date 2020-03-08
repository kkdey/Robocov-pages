
N=500
P=50
count_missing=0
tmp = get(load(paste0("/n/groups/price/kushal/Robocov/output_Nov26_prec/robocov_prec_sim_hub_n_",
                      N, "_p_", P, "_prop_", count_missing, ".rda")))
summ=matrix(0, nrow(tmp[[1]]), ncol(tmp[[1]]))
for(mm in 1:10){
  summ = summ+tmp[[mm]]
}

xtable2(t(summ/10))
