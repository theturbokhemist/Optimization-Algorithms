#Create topology
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
# topology_AB <- data.frame(Source = c("Z", "Z", "C"), Target = c("B", "D", "D"), Type = c(1,1, 1))
# test6 <- sracipeSimulate(topology_AB, integrate = FALSE, genParams = TRUE)
# names(test6)
# sracipeConfig(reference_rset)
# colnames(sracipeParams(test6))
# sracipePlotCircuit(test6, plotToFile = F)
# threshold_genes1 <- grep(pattern = "TH", colnames(sracipeParams(reference_rset)))
# threshold_medians1 <- sracipeConfig(reference_rset)$thresholds
# threshold_names1 <- names(reference_rset)
# 
# sracipeParams(reference_rset)
# 
# threshold_values1 <- vector(length = length(threshold_genes1))

# for (k in 1:length(threshold_names1)) {
#   
#   names_index <- grep(paste0("TH_", threshold_names1[k]), colnames(sracipeParams(reference_rset))[threshold_genes1])
#   threshold_values1[c(names_index)] <- threshold_medians1[k]
#   
# }
# 
# threshold_df1 <- data.frame(interactions = c(colnames(sracipeParams(reference_rset))[threshold_genes1]), indices = threshold_genes1, medians = threshold_values1)


score_distribution <- function(topology, numModels = 2000, ref_bounds = c(1,10), sim_bounds = c(1,10), method = "PCs", number_of_iterations = 500, step_size = 0.02) {
  
  scores_individual <- vector(mode = "list", length = number_of_iterations)
  
  scores_sum <- vector(length = number_of_iterations)
  
  scores_object <- list(scores_sep = scores_individual, scores_sum = scores_sum)
  
  #generate unintegrated rset object
  rset_reference <- sracipeSimulate(topology, integrate = FALSE, numModels = numModels, genParams = TRUE)
  
  #replace params
  sracipeParams(rset_reference)[1] <- runif(min = ref_bounds[1], max = ref_bounds[2], n = numModels)
  # sracipeParams(rset_TS_10000_reference)[2] <- runif(min = 70, max = 95, n = 10000)
  # sracipeParams(rset_TS_10000_reference)[3] <- runif(min = 0.2, max = 0.6, n = 10000)
  # sracipeParams(rset_TS_10000_reference)[4] <- runif(min = 0.3, max = 0.7, n = 10000)
  # sracipeParams(rset_TS_10000_reference)[5] <- runif(min = 0.5, max = 1.5, n = 10000)*sracipeConfig(rset_TS_10000_reference)$thresholds[2]
  # sracipeParams(rset_TS_10000_reference)[6] <- runif(min = 0.02, max = 1.98, n = 10000)*sracipeConfig(rset_TS_10000_reference)$thresholds[1]
  # sracipeParams(rset_TS_10000_reference)[7] <- sample(1:3, size = 10000, replace = TRUE)
  # sracipeParams(rset_TS_10000_reference)[8] <- sample(4:5, size = 10000, replace = TRUE)
  # sracipeParams(rset_TS_10000_reference)[9] <- runif(min = 15, max = 45, n = 10000)
  # sracipeParams(rset_TS_10000_reference)[10] <- runif(min = 30, max = 60, n = 10000)
  
  #simulate reference
  rset_reference <- sracipeSimulate(rset_reference, integrate = TRUE, genParams = FALSE, integrateStepSize = step_size)
  
  #log2 transform expresion data
  expr_reference_log <- log2(t(assay(rset_reference)))
  
  means_reference <- colMeans(expr_reference_log)
  
  sds_reference <- apply(expr_reference_log, 2, sd)
  
  #z normalize
  expr_reference_log_z <- sweep(expr_reference_log, 
                                2, means_reference, FUN = "-")
  
  expr_reference_log_z <- sweep(expr_reference_log_z, 
                                2, sds_reference, FUN = "/")
  
  #find eigenvector,loading scores
  pc_reference <- prcomp(expr_reference_log_z, center = FALSE, scale. = FALSE)
  
  
  
  list_fvs <- vector(mode = "list", length = number_of_iterations)
  
  for (i in 1:number_of_iterations) {
    
    
    rset_sim <- sracipeSimulate(topology, integrate = FALSE, numModels = numModels, genParams = TRUE)
    
    sracipeParams(rset_sim)[1] <- runif(min = sim_bounds[1], max = sim_bounds[2], n = numModels)
    # sracipeParams(rset_TS_10000_simulated)[2] <- runif(min = 70, max = 95, n = 10000)
    # sracipeParams(rset_TS_10000_simulated)[3] <- runif(min = 0.2, max = 0.6, n = 10000)
    # sracipeParams(rset_TS_10000_simulated)[4] <- runif(min = 0.3, max = 0.7, n = 10000)
    # sracipeParams(rset_TS_10000_simulated)[5] <- runif(min = 0.5, max = 1.5, n = 10000)*sracipeConfig(rset_TS_10000_reference)$thresholds[2]
    # sracipeParams(rset_TS_10000_simulated)[6] <- runif(min = 0.02, max = 1.98, n = 10000)*sracipeConfig(rset_TS_10000_reference)$thresholds[1]
    # sracipeParams(rset_TS_10000_simulated)[7] <- sample(1:3, size = 10000, replace = TRUE)
    # sracipeParams(rset_TS_10000_simulated)[8] <- sample(4:5, size = 10000, replace = TRUE)
    # sracipeParams(rset_TS_10000_simulated)[9] <- runif(min = 15, max = 45, n = 10000)
    # sracipeParams(rset_TS_10000_simulated)[10] <- runif(min = 30, max = 60, n = 10000)
    
    rset_sim <- sracipeSimulate(rset_sim, integrate = TRUE, genParams = FALSE, integrateStepSize = step_size)
    
    expr_sim_log <- log2(t(assay(rset_sim)))
    
    expr_sim_log_z <- sweep(expr_sim_log, 2, means_reference, FUN = "-")
    
    expr_sim_log_z <- sweep(expr_sim_log_z, 2, sds_reference, FUN = "/")
    
    
    if(method == "PCs") {
      
      simulated_expression_PCs <- expr_sim_log_z %*% pc_reference$rotation
      
      #set up vector to store stat for each gene of the network
      stats_vec <- vector(length = ncol(simulated_expression_PCs))
      
      
      for (l in 1:length(stats_vec)) {
        
        ks <- dgof::ks.test(x = pc_reference$x[,l], y = simulated_expression_PCs[,l])
        stats_vec[l] <- as.numeric(ks$statistic)
        
      }
      
    } else {
      
      #set up vector to store stat for each gene of the network
      stats_vec <- vector(length = ncol(expr_sim_log_z))
      
      
      for (l in 1:length(stats_vec)) {
        
        ks <- dgof::ks.test(x = expr_reference_log_z[,l], y = expr_sim_log_z[,l])
        stats_vec[l] <- as.numeric(ks$statistic)
        
      }
      
    }
    
    scores_object$scores_sep[[i]] <- stats_vec
    
    scores_object$scores_sum[i] <- sum(stats_vec)
    
    
  }
  
  scores_object
}

#score_distribution(topology = topology_TS, numModels = 2000, ref_bounds = c(1,10), sim_bounds = c(1,10), method = "PCs", number_of_iterations = 2, step_size = 0.02)

#############################################################################################################################################################################################
topology = topology_TS
numModels = 1000
ref_bounds = c(1,10)
sim_bounds = c(1,10)
method = "Default"
number_of_iterations = 2
step_size = 0.02

arguments <- list(topology = topology_TS,
                  numModels = numModels,
                  ref_bounds =ref_bounds,
                  sim_bounds = sim_bounds,
                  method = method,
                  number_of_iterations = number_of_iterations,
                  step_size = step_size)

list_of_scores_objects_TS_10000_default_0.02 <- vector(mode = "list", length = 11)
list_of_scores_objects_TS_10000_default_0.02[[11]] <- arguments
names(list_of_scores_objects_TS_10000_default_0.02)[[11]] <- "arguments"

for (j in 1:10) {
  
  names(list_of_scores_objects_TS_10000_default_0.02)[j] <- paste0("sim_bounds_", sim_bounds[1], "_", sim_bounds[2])
  
  list_of_scores_objects_TS_10000_default_0.02[[j]] <- score_distribution(topology = topology, numModels = numModels, ref_bounds = ref_bounds, sim_bounds = sim_bounds, method = method, number_of_iterations = number_of_iterations, step_size = step_size)
  
  sim_bounds <- sim_bounds + 10
  
  save.image("~/Desktop/Modeling/Metropolis/score_term_refinement_env3.RData")
  
  print(j)
  
  
}

########################################################################################################################################################################
fplot(density(unlist(lapply(fv_vec_test_stat_list_twicenorm_onceintegrate, sum))))
range(unlist(lapply(fv_vec_test_stat_list_twicenorm_twiceintegrate, sum)))
range(unlist(lapply(fv_vec_test_stat_list_twicenorm_onceintegrate, sum)))
range(unlist(lapply(fv_vec_test_stat_list_oncenorm_onceintegrate, sum)))
range(unlist(lapply(fv_vec_pcomp_stat_list_oncenorm_onceintegrate, sum)))
plot(density(unlist(lapply(fv_vec_pcomp_stat_list_oncenorm_onceintegrate, sum))))
range(unlist(lapply(b, sum)))
range(unlist(lapply(b_10000, sum)))

plot(density(unlist(lapply(fv_vec_pcomp_stat_list_oncenorm_onceintegrate, sum))))

list_of_scores_objects_TS_10000_default_0.02


