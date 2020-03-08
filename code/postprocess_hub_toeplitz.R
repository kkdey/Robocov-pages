
#################################  postprocess runs from cluster  #####################################

Nvec=c(50, 100, 200, 500)
Pvec=100
propvec = c(0, 0.25, 0.5, 0.75)

for (n in 1:length(Nvec)){
  for(p in 1:length(propvec)){
    if(file.exists(paste0("/Users/kushaldey/Documents/Robocov-pages/output/output_sparse/robocov_sim_toeplitz_n_",
                          Nvec[n], "_p_", Pvec, "_prop_", propvec[p], ".rda"))){
      dat = get(load(paste0("/Users/kushaldey/Documents/Robocov-pages/output/output_sparse/robocov_sim_toeplitz_n_",
                            Nvec[n], "_p_", Pvec, "_prop_", propvec[p], ".rda")))
    }
  }
}

dat = get(load(paste0("/Users/kushaldey/Documents/Robocov-pages/output/output_sparse/robocov_sim_toeplitz_n_",
                      500, "_p_", 100, "_prop_", 0.5, ".rda")))
xtable2(t(dat[[1]][c(1, 3, 5), c(1,3,2,4,5)]))
