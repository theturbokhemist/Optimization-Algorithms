

individual <- function(topology, sample_rset, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  number_of_genes <- ncol(sracipeParams(sample_rset))
  parameter_names <-colnames(sracipeParams(sample_rset))
  
  ind <- vector(mode = "list", length = number_of_genes)
  names(ind) <- parameter_names
  
  production_genes <- grep(pattern = "G", parameter_names)
  degradation_genes <- grep(pattern = "K", parameter_names)
  threshold_genes <- grep(pattern = "TH", parameter_names)
  hill_genes <- grep(pattern = "N", parameter_names)
  FC_genes <- grep(pattern = "FC", parameter_names)
  
  bounds <- vector(mode = "list", length = 2)
  
  #Production
  for (i in 1:length(production_genes)) {
    lower_limit_rnd <- runif(G_lower, G_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (G_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    ind[[production_genes[i]]] <- bounds
  }
  
  #Degradation
  for (j in 1:length(degradation_genes)) {
    lower_limit_rnd <- runif(K_lower, K_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (K_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    ind[[degradation_genes[j]]] <- bounds
  }
  
  #Threshold
  for (l in 1:length(threshold_genes)) {
    lower_limit_rnd <- runif(TH_lower, TH_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (TH_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    ind[[threshold_genes[l]]] <- bounds
  }
  
  
  #Hill
  for (m in 1:length(hill_genes)) {
    lower_limit_rnd <- sample(N_lower:N_upper, size = 1, replace = TRUE)
    upper_limit_rnd <- lower_limit_rnd + round((N_upper - lower_limit_rnd)*runif(0, 1, n = 1))
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    ind[[hill_genes[m]]] <- bounds
  }
  
  #FC
  for (n in 1:length(FC_genes)) {
    lower_limit_rnd <- runif(FC_lower, FC_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (FC_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    ind[[FC_genes[n]]] <- bounds
  }
  
  ind
  
}

population <- function(topology, sample_rset, population_size = 8, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  population_list <- vector(mode = "list", length = population_size)
  for (i in 1:population_size) {
    
    population_list[[i]] <- individual(topology = topology , sample_rset = sample_rset,  G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                                       N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
  }
  population_list
  
}


fv <- function(individual, reference_expression, sample_rset, threshold_df, means_ref, sds_ref, pc_ref, numModels = 2000, scoring_method = "PCs") {
  
  simulated_rset <- sample_rset
  
  hill_genes <- grep(pattern = "N", colnames(sracipeParams(simulated_rset)))
  
  
  for (i in 1:length(individual)) {
    
    sracipeParams(simulated_rset)[,i] <- runif(numModels, min= as.numeric(individual[[i]][1]), max= as.numeric(individual[[i]][2]))  
    
  }
  
  for (j in 1:nrow(threshold_df)) {

    sracipeParams(simulated_rset)[,threshold_df[j, "indices"]] <- sracipeParams(simulated_rset)[,threshold_df[j, "indices"]]*threshold_df[j, "medians"]  
    
  }
  
  for (k in 1:length(hill_genes)) {
    
    sracipeParams(simulated_rset)[,hill_genes[k]] <- round(sracipeParams(simulated_rset)[,hill_genes[k]])
    
  }


  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, numModels = numModels, genParams = FALSE)

  simulated_expression <- log2(t(assay(simulated_rset)))

  simulated_expression_z <- sweep(simulated_expression,
                                  2, means_ref, FUN = "-")

  simulated_expression_z <- sweep(simulated_expression_z,
                                  2, sds_ref, FUN = "/")

if (scoring_method == "PCs") {
  
  simulated_expression_PCs <- simulated_expression_z %*% pc_ref$rotation
  
  fv_vec <-  vector(length = ncol(pc_ref$x))
  
  for (l in 1:ncol(pc_ref$x)) {
    
    ks_test <- dgof::ks.test(x = pc_ref$x[,l], y = simulated_expression_PCs[,l])
    fv_vec[l] <- as.numeric(ks_test$statistic)
    
  }
  
  
} else {
  fv_vec <-  vector(length = ncol(reference_expression))

  for (l in 1:ncol(reference_expression)) {

    ks_test <- dgof::ks.test(x = reference_expression[,l], y = simulated_expression_z[,l])
    fv_vec[l] <- as.numeric(ks_test$statistic)

  }

}
  fv_vec

}

selection <- function(generation, reference_expression, sample_rset, threshold_df, means_ref, sds_ref, pc_ref, numModels = 2000, scoring_method = "PCs") {
  fitness_values_list <- vector(mode = "list", length = length(generation))
  fitness_values <- vector(length = length(generation))
  
  for (i in 1:length(generation)) {

    fitness_values_list[[i]] <- fv(generation[[i]], reference_expression = reference_expression, threshold_df = threshold_df, sample_rset = sample_rset, numModels = numModels,
                                   means_ref = means_ref, sds_ref = sds_ref, pc_ref = pc_ref, scoring_method = scoring_method)
    fitness_values[i] <- sum(fitness_values_list[[i]])
  }
  
  fv_index_ordered <- order(fitness_values, decreasing = F)
  fitness_values_list <<- fitness_values_list[c(fv_index_ordered)]
  generation <- generation[c(fv_index_ordered)]
  selected <- generation[1:(length(generation)/2)]
  selected
}


mating <- function(mating_pair) {
  genes <- length(mating_pair[[1]])
  child1 <- mating_pair[[1]]
  child2 <- mating_pair[[2]]
  
  offspring <- vector(mode = "list", length = 2)
  pivot_point <- sample(2:genes, 1)
  
  if (genes == 2) {
    
    pivot_point <- 2
    
  }
  
  if (genes == 1) {
    
    pivot_point <- 1
    
  }
  
  child1[pivot_point:genes] <- child2[pivot_point:genes]
  child2[pivot_point:genes] <- mating_pair[[1]][pivot_point:genes]
  
  offspring[[1]] <- child1
  offspring[[2]] <- child2
  offspring
}

mutation <- function(individual, sample_rset, mutation_rate = 1, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  number_of_genes <- ncol(sracipeParams(sample_rset))
  parameter_names <-colnames(sracipeParams(sample_rset))
  
  ind <- vector(mode = "list", length = number_of_genes)
  names(ind) <- parameter_names
  
  production_genes <- grep(pattern = "G", parameter_names)
  degradation_genes <- grep(pattern = "K", parameter_names)
  threshold_genes <- grep(pattern = "TH", parameter_names)
  hill_genes <- grep(pattern = "N", parameter_names)
  FC_genes <- grep(pattern = "FC", parameter_names)
  
  mutation_indexes <- sample(1:length(individual), mutation_rate)
  
  for (l in 1:mutation_rate) {
    
    bound <- round(runif(1, 0, 1))

    
    if (mutation_indexes[l] %in% production_genes) {
      
      if (bound == 0) {
        
        individual[[mutation_indexes[l]]][1] <- runif(min = G_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)

        
      } else {
        
        individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = G_upper, n = 1)
        
      }
    }
    
    
    if (mutation_indexes[l] %in% degradation_genes) {
      
      if (bound == 0) {
        
        individual[[mutation_indexes[l]]][1] <- runif(min = K_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
        
        
      } else {
        
        individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = K_upper, n = 1)
        
      }
      
    }
    
    
    if (mutation_indexes[l] %in% threshold_genes) {
      
      if (bound == 0) {
        
        individual[[mutation_indexes[l]]][1] <- runif(min = TH_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
        
        
      } else {
        
        individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = TH_upper, n = 1)
        
      }
      
    }
    
    if (mutation_indexes[l] %in% hill_genes) {
      
      if (bound == 0) {
        
        individual[[mutation_indexes[l]]][1] <- round(runif(min = N_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1))
        
        
      } else {
        
        individual[[mutation_indexes[l]]][2] <- round(runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = N_upper, n = 1))
        
      }
      
    }
    
    if (mutation_indexes[l] %in% FC_genes) {
      
      if (bound == 0) {
        
        individual[[mutation_indexes[l]]][1] <- runif(min = FC_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
        
        
      } else {
        
        individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = FC_upper, n = 1)
        
      }
      
    }
    

  }
  
  individual
  
}

reproduction <- function(selected, sample_rset, mutation_rate = 1, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  pairs <- seq(1, length(selected), 2)
  unmutated <- list()
  mutated <- vector(mode = "list", length = length(selected))
  
  for (i in 1:length(pairs)) {
    
    unmutated <- append(unmutated, mating(selected[pairs[[i]]:(pairs[[i]] + 1)]))
    # unmutated <- do.call(c, list(unmutated, mating(selected[pairs[[i]]:(pairs[[i]] + 1)])))
    
  }
  for (j in 1:length(unmutated)) {
    
    mutated[[j]] <- mutation(unmutated[[j]], sample_rset = sample_rset,  G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                             N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
  }
  
  new_generation <- append(selected, mutated)
  new_generation
  
}
 

GA_R_RACIPE <- function(topology, reference_rset, population_size = 16, number_of_generations = 10, numModels = 2000, mutation_rate = 1, scoring_method = "PCs", 
                        G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100 ) {
  
  reference_rset <- sracipeSimulate(reference_rset, integrate = TRUE, numModels = numModels, genParams = FALSE)
  
  #transform/normalize/pc reference data
  reference_expression <- log2(t(assay(reference_rset)))
  
  means_reference <- colMeans(reference_expression)
  
  sds_reference <- apply(reference_expression, 2, sd)
  
  reference_expression_log_z <- sweep(reference_expression, 
                                                   2, means_reference, FUN = "-")
  
  reference_expression_log_z <- sweep(reference_expression_log_z, 
                                                   2, sds_reference, FUN = "/")
  
  reference_expression <- reference_expression_log_z
  pc_reference <- prcomp(reference_expression, center = FALSE, scale. = FALSE)
  
  #thresholds
  threshold_genes <- grep(pattern = "TH", colnames(sracipeParams(reference_rset)))
  threshold_medians <- sracipeConfig(reference_rset)$thresholds
  threshold_names <- names(reference_rset)
  
  threshold_values <- vector(length = length(threshold_genes))
  
  for (k in 1:length(threshold_names)) {
    
    names_index <- grep(paste0("TH_", threshold_names[k]), colnames(sracipeParams(reference_rset))[threshold_genes])
    threshold_values[c(names_index)] <- threshold_medians[k]
    
  }
  
  threshold_df <- data.frame(interactions = c(colnames(sracipeParams(reference_rset))[threshold_genes]), indices = threshold_genes, medians = threshold_values)
  
  #generate intial population
  generations_list <<- vector(mode = "list", length = number_of_generations)
  fvs_list <<- vector(mode = "list", length = number_of_generations)

  generation <- population(topology = topology, sample_rset = reference_rset, population_size = population_size, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                           K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
 

  for (i in 1:number_of_generations) {
    print(paste("Generation:", i))
    generations_list[[i]] <<- generation
    selected <- selection(generation = generation, reference_expression = reference_expression, sample_rset = reference_rset, threshold_df = threshold_df, numModels = numModels,
                          means_ref = means_reference, sds_ref = sds_reference, pc_ref = pc_reference, scoring_method = scoring_method)
  
    fvs_list[[i]] <<- fitness_values_list


    generation <- reproduction(selected = selected, mutation_rate = mutation_rate, sample_rset = reference_rset, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                               K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    generation

  }
  
  
  
  
}

#######
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE)
sracipeParams(reference_rset)[1] <- runif(min = 5, max = 30, n = 10000)
sracipeParams(reference_rset)[2] <- runif(min = 70, max = 95, n = 10000)
#reference_rset <- sracipeSimulate(reference_rset, integrate = TRUE, numModels = 4000, genParams = FALSE)

###
GA_R_RACIPE(topology_TS, reference_rset = reference_rset, population_size = 40, number_of_generations = 40, mutation_rate = 2, numModels = 10000, N_upper = 5)
generations_list[[40]][[1]]
fvs_list[[40]]

assay(reference_rset)

generations_list[[2]][[1]]
generations_list[[10]][[1]]
fvs_list[[178]]
generations_list[[178]][[1]]

# generations_list1_100_16_TS <- generations_list
# fv_list1_100_16_TS <- fvs_list
# 
# generations_list3_100_16_TS <- generations_list
# fv_list3_100_16_TS <- fvs_list

test_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 2000, genParams = TRUE)

poptest <- population(topology_TS, sample_rset = test_rset, population_size = 8)
fv(individual = poptest[[1]], reference_expression = rSet_TS_2000_custom_expression_log_z, sample_rset = rset_TS_2000_custom, threshold_df = test_t_df, numModels = 2000)
seltest <- selection(generation = poptest, reference_expression = rSet_TS_2000_custom_expression_log_z, sample_rset = rset_TS_2000_custom, threshold_df = test_t_df, numModels = 2000)
gentest <- reproduction(seltest, sample_rset = test_rset, mutation_rate = 1)

gentest[c(1,5)]

umut <- poptest[[1]]
mut <- mutation(individual = poptest[[1]], sample_rset = test_rset, mutation_rate = 3)

for (i in 1:length(umut)) {
  
  if ((as.numeric(umut[[i]][1]) != as.numeric(mut[[i]][1])) | (as.numeric(umut[[i]][2]) != as.numeric(mut[[i]][2]))) {
    
    print(umut[[i]])
    print(mut[[i]])
  }
  
}

umut[[1]] == mut[[1]]

####
rset_TS_2000_custom <- rSet_TS_2000
colnames(sracipeParams(rset_TS_2000_custom))

t_genes <- grep(pattern = "TH", colnames(sracipeParams(rset_TS_2000_custom)))
t_medians <- sracipeConfig(rset_TS_2000_custom)$thresholds
test_t_df <- data.frame(interactions = c(colnames(sracipeParams(rset_TS_2000_custom))[t_genes]), indices = t_genes, medians = t_medians)

#######
fv_vec_test_list <- vector(mode = "list", length = 20)
fv_vec_test_sum <- vector(length = 20)

for (i in 1:20) {

rset_TS_2000_custom <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 2000, genParams = TRUE)
sracipeParams(rset_TS_2000_custom)[1] <- runif(min = 5, max = 30, n = 2000)
sracipeParams(rset_TS_2000_custom)[2] <- runif(min = 70, max = 95, n = 2000)
sracipeParams(rset_TS_2000_custom)[3] <- runif(min = 0.2, max = 0.6, n = 2000)
sracipeParams(rset_TS_2000_custom)[4] <- runif(min = 0.3, max = 0.7, n = 2000)
sracipeParams(rset_TS_2000_custom)[5] <- runif(min = 0.5, max = 1.5, n = 2000)*sracipeConfig(rset_TS_2000_custom)$thresholds[2]
sracipeParams(rset_TS_2000_custom)[6] <- runif(min = 0.02, max = 1.98, n = 2000)*sracipeConfig(rset_TS_2000_custom)$thresholds[1]
sracipeParams(rset_TS_2000_custom)[7] <- sample(1:3, size = 2000, replace = TRUE)
sracipeParams(rset_TS_2000_custom)[8] <- sample(4:5, size = 2000, replace = TRUE)
sracipeParams(rset_TS_2000_custom)[9] <- runif(min = 15, max = 45, n = 2000)
sracipeParams(rset_TS_2000_custom)[10] <- runif(min = 30, max = 60, n = 2000)

rset_TS_2000_custom <- sracipeSimulate(rset_TS_2000_custom, integrate = TRUE, numModels = 2000, genParams = FALSE)

rSet_TS_2000_custom_expression_log <- log2(t(assay(rset_TS_2000_custom)))

rSet_TS_2000_custom_expression_log_z <- sweep(rSet_TS_2000_custom_expression_log, 
                                        2, colMeans(rSet_TS_2000_custom_expression_log), FUN = "-")

rSet_TS_2000_custom_expression_log_z <- sweep(rSet_TS_2000_custom_expression_log_z, 
                                        2, apply(rSet_TS_2000_custom_expression_log, 2, sd), FUN = "/")



rset_TS_2000_test <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 2000, genParams = TRUE)

sracipeParams(rset_TS_2000_test)[1] <- runif(min = 5, max = 30, n = 2000)
sracipeParams(rset_TS_2000_test)[2] <- runif(min = 70, max = 95, n = 2000)
sracipeParams(rset_TS_2000_test)[3] <- runif(min = 0.2, max = 0.6, n = 2000)
sracipeParams(rset_TS_2000_test)[4] <- runif(min = 0.3, max = 0.7, n = 2000)
sracipeParams(rset_TS_2000_test)[5] <- runif(min = 0.5, max = 1.5, n = 2000)*sracipeConfig(rset_TS_2000_test)$thresholds[1]
sracipeParams(rset_TS_2000_test)[6] <- runif(min = 0.02, max = 1.98, n = 2000)*sracipeConfig(rset_TS_2000_test)$thresholds[2]
sracipeParams(rset_TS_2000_test)[7] <- sample(1:3, size = 2000, replace = TRUE)
sracipeParams(rset_TS_2000_test)[8] <- sample(4:5, size = 2000, replace = TRUE)
sracipeParams(rset_TS_2000_test)[9] <- runif(min = 15, max = 45, n = 2000)
sracipeParams(rset_TS_2000_test)[10] <- runif(min = 30, max = 60, n = 2000)

rset_TS_2000_custom <- sracipeSimulate(rset_TS_2000_test, integrate = TRUE, numModels = 2000, genParams = FALSE)

rset_TS_2000_test_expression_log <- log2(t(assay(rset_TS_2000_custom)))

rset_TS_2000_test_expression_log_z <- sweep(rset_TS_2000_test_expression_log, 
                                              2, colMeans(rset_TS_2000_test_expression_log), FUN = "-")

rset_TS_2000_test_expression_log_z <- sweep(rset_TS_2000_test_expression_log_z, 
                                              2, apply(rset_TS_2000_test_expression_log, 2, sd), FUN = "/")


fv_vec_test <-  vector(length = ncol(rSet_TS_2000_custom_expression_log_z))

for (l in 1:ncol(rSet_TS_2000_custom_expression_log_z)) {
  
  ks_test1 <- dgof::ks.test(x = rSet_TS_2000_custom_expression_log_z[,l], y = rset_TS_2000_test_expression_log_z[,l])
  fv_vec_test[l] <- as.numeric(10000- 10000*ks_test1$statistic)
  
  
}
fv_vec_test_list[[i]] <- fv_vec_test
fv_vec_test_sum[i] <- 20000 - sum(fv_vec_test)

}

range(unlist(lapply(fv_vec_test_list, function(x) {20000 - sum(x)})))

#####
prcomp_rset_TS_2000_test_expression_log_z <- prcomp(x = rset_TS_2000_test_expression_log_z, center = F, scale. = F)
prcomp_rset_TS_2000_test_expression_log_z$rotation

screeplot_2000_DP_t_filtered_log_z <- prcomp_2000_DP_t_filtered_log_z$sdev^2/sum(prcomp_2000_DP_t_filtered_log_z$sdev^2)*100
screeplot_2000_DP_t_filtered_log_z <- data.frame(PCs = factor(colnames(prcomp_2000_DP_t_filtered_log_z$x), levels = colnames(prcomp_2000_DP_t_filtered_log_z$x) ),
                                                 Percent_Variance = screeplot_2000_DP_t_filtered_log_z)

screeplot_2000_DP_t_filtered_log_z <- ggplot2::ggplot(data = screeplot_2000_DP_t_filtered_log_z, aes(x=PCs, y = Percent_Variance)) + geom_bar(stat = "identity")

PCA_2000_DP_t_filtered_log_z <- ggplot(data = as.data.frame(prcomp_2000_DP_t_filtered_log_z$x), aes(x = PC1, y = PC2)) + geom_point()
prcomp_2000_DP_t_filtered_log_z$x

den3d <- kde2d(x = prcomp_2000_DP_t_filtered_log_z$x[, 1], y = prcomp_2000_DP_t_filtered_log_z$x[, 2])
plot_ly(x = den3d$x, y = den3d$y, z = den3d$z) %>% add_surface()
#####

plot_density <- ggplot2::ggplot(data = data.frame(expression = rSet_TS_2000_custom_expression_log_z[,1]), aes(x = expression)) + geom_density() +
  geom_density(data = data.frame(expression = rset_TS_2000_test_expression_log_z[,1]), aes(x = expression), color = "red")
print(plot_density)

list1 <- list()
list2 <- list(list(1,2), list(list(3,4), list(5,6)))
list3 <- list(1,2,3)
list4 <- append(list1, list3)
append(list4, list4)

a <- sample(x = 6:8, size = 500, replace = TRUE)
b <- sample(x = 8:10, size = 200, replace = TRUE)
c <- sample(x = 1:5, size = 300, replace = TRUE)

a <- runif(n = 500, min = 6, max = 8)
b <- runif(n = 200, min = 8, max = 10)
c <- runif(n = 300, min = 1, max = 6)

plot(density(c(a,b,c)), main="PDF for kevins level of prickness at any given moment", 
     xlab="level of prickness (out of 10)", ylab="probability density")
plot(c(a,b,c))

plot(ecdf(c(a,b,c)), main="CDF for kevins level of prickness at any given moment", 
     xlab="level of prickness (out of 10)", ylab="probability that kevins prickness will take any value less than the x-value" )
