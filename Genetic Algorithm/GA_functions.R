library(sRACIPE)
library(dgof)
library(rlist)
library(plyr)


###################Reference_Info###########################################################################################################################################################################################
# topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
# reference_rset_10000_0.02 <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE, integrateStepSize = 0.02)
# 
# ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
# names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")
###################generate_ref_individual_object###########################################################################################################################################################################################


generate_ref_individual_object <- function(topology, reference_individual, numModels = 10000, integrateStepSize = 0.02) {
  
  reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  #hill gene indexes
  hill_genes <- grep(pattern = "N", colnames(sracipeParams(reference_rset)))
  #assign vaues from uniform samples distributions to each model of reference rset.
  for (i in 1:length(reference_individual)) {
    sracipeParams(reference_rset)[,i] <- runif(numModels, min= as.numeric(reference_individual[[i]][1]), max= as.numeric(reference_individual[[i]][2]))
  }
  for (k in 1:length(hill_genes)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(reference_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(reference_individual[[hill_genes[k]]][1]) - 0.5, max= as.numeric(reference_individual[[hill_genes[k]]][2]) + 0.5))
  }
  
  ref_expression_rset <- sracipeSimulate(reference_rset, integrate = TRUE, genParams = FALSE)
  expr_reference_log <- log2(t(assay(ref_expression_rset)))
  means_reference <- colMeans(expr_reference_log)
  sds_reference <- apply(expr_reference_log, 2, sd)
  #z normalize
  expr_reference_log_z <- sweep(expr_reference_log,
                                2, means_reference, FUN = "-")
  expr_reference_log_z <- sweep(expr_reference_log_z,
                                2, sds_reference, FUN = "/")
  #find eigenvector,loading scores
  prcomp_reference <- prcomp(expr_reference_log_z, center = FALSE, scale. = FALSE)
  #list of all the ref information
  ref_individual_object <- list(reference_rset = ref_expression_rset, reference_individual = reference_individual,
                                expression_ref = expr_reference_log_z, log2_expression_ref = expr_reference_log,
                                prcomp_ref = prcomp_reference,
                                means_ref = means_reference, sds_ref = sds_reference, topology = topology)
  ref_individual_object
}

reference_individual_object_10000_0.02 <- generate_ref_individual_object(topology = topology_TS, reference_individual = ref_state1, numModels = 10000, integrateStepSize = 0.02)
# reference_individual_object_10000_0.02

##########################generate_fitness_value_object####################################################################################################################################################################################

generate_fitness_value_object <- function(individual, reference_individual_object, normalize_with_reference = FALSE, scoring_method = "PCs",
                                          numModels = 10000, integrateStepSize = 0.02) {
  
  if(scoring_method == "Toy") {
    fv_vec <- vector(length = length(individual))
    for (i in 1:length(individual)) {
      fv_vec[[i]] <- as.numeric(individual[[i]][1]) - as.numeric(individual[[i]][2])
    }
    
    fitness_value_object <- list(fitness_values_sum = sum(fv_vec), simulated_rset = reference_individual_object$reference_rset, individual = individual)
    
  } else {
    
    simulated_rset <- sracipeSimulate(reference_individual_object$topology, integrate = FALSE,
                                      numModels = numModels, genParams = TRUE,
                                      integrateStepSize = integrateStepSize)
    #hill gene indexes
    hill_genes <- grep(pattern = "N", colnames(sracipeParams(simulated_rset)))
    #assign vaues from uniform samples distributions to each model of simulated rset.
    for (i in 1:length(individual)) {
      sracipeParams(simulated_rset)[,i] <- runif(numModels, min= as.numeric(individual[[i]][1]), max= as.numeric(individual[[i]][2]))
    }
    #sample and round all the hill coefficients
    for (k in 1:length(hill_genes)) {
      #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
      sracipeParams(simulated_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(individual[[hill_genes[k]]][1]) - 0.5, max= as.numeric(individual[[hill_genes[k]]][2]) + 0.5))
    }
    
    sim_expression_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE,
                                           integrateStepSize = integrateStepSize )
    expr_simulated_log <- log2(t(assay(sim_expression_rset)))
    
    
    if (normalize_with_reference == TRUE) {
      
      #z normalize
      expr_simulated_log_z <- sweep(expr_simulated_log,
                                    2, reference_individual_object$means_ref, FUN = "-")
      expr_simulated_log_z <- sweep(expr_simulated_log_z,
                                    2, reference_individual_object$sds_ref, FUN = "/")
      #rotate
      simulated_expression_PCs <- expr_simulated_log_z %*% reference_individual_object$prcomp_ref$rotation
      
      #fitness_value vector
      fv_vec <- vector(length = ncol(expr_simulated_log_z))
      
      #scoring
      if (scoring_method == "PCs") {
        
        #ks score for distribution of each gene
        for (l in 1:length(fv_vec)) {
          ks_test <- dgof::ks.test(x = reference_individual_object$prcomp_ref$x[,l], y = simulated_expression_PCs[,l])
          fv_vec[l] <- as.numeric(ks_test$statistic)
        }
        
      } else {
        
        for (l in 1:length(fv_vec)) {
          ks_test <- dgof::ks.test(x = reference_individual_object$expression_ref[,l], y = expr_simulated_log_z[,l])
          fv_vec[l] <- as.numeric(ks_test$statistic)
        }
      }
      
    } else {
      
      #calculate means/sds for sim data
      means_sim <- colMeans(expr_simulated_log)
      sds_sim <- apply(expr_simulated_log, 2, sd)
      
      #z normalize
      expr_simulated_log_z <- sweep(expr_simulated_log,
                                    2, means_sim, FUN = "-")
      expr_simulated_log_z <- sweep(expr_simulated_log_z,
                                    2, sds_sim, FUN = "/")
      #calculate PC's
      prcomp_simulated <- prcomp(expr_simulated_log_z, center = FALSE, scale. = FALSE)
      simulated_expression_PCs <- prcomp_simulated$x
      
      #fitness_value vector
      fv_vec <- vector(length = ncol(expr_simulated_log_z))
      
      #scoring
      if (scoring_method == "PCs") {
        
        #ks score for distribution of each gene
        for (l in 1:length(fv_vec)) {
          ks_test <- dgof::ks.test(x = reference_individual_object$prcomp_ref$x[,l], y = simulated_expression_PCs[,l])
          fv_vec[l] <- as.numeric(ks_test$statistic)
        }
        
      } else {
        
        for (l in 1:length(fv_vec)) {
          ks_test <- dgof::ks.test(x = reference_individual_object$expression_ref[,l], y = expr_simulated_log_z[,l])
          fv_vec[l] <- as.numeric(ks_test$statistic)
        }
      }
    }
    
    
    fitness_value_object <- list(fitness_values = fv_vec, fitness_values_sum = sum(fv_vec), simulated_rset = sim_expression_rset, individual = individual, expression_simulated_znormalized = expr_simulated_log_z, PCs_simulated = simulated_expression_PCs)
    
  }
}

#########################individual#####################################################################################################################################################################################
individual <- function(sample_rset, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  number_of_genes <- ncol(sracipeParams(sample_rset))
  parameter_names <-colnames(sracipeParams(sample_rset))
  
  #list that is length of number of genes
  ind <- vector(mode = "list", length = number_of_genes)
  names(ind) <- parameter_names
  
  #index of each parameter type in vector of parameter names 
  production_genes <- grep(pattern = "G", parameter_names)
  degradation_genes <- grep(pattern = "K", parameter_names)
  threshold_genes <- grep(pattern = "TH", parameter_names)
  hill_genes <- grep(pattern = "N", parameter_names)
  FC_genes <- grep(pattern = "FC", parameter_names)
  
  #list that will contain lower and upper bound for each parameter
  bounds <- vector(mode = "list", length = 2)
  
  #Production
  for (i in 1:length(production_genes)) {
    
    #sample uniformly from absolute bounds
    lower_limit_rnd <- runif(G_lower, G_upper, n = 1)
    upper_limit_rnd <- runif(G_lower, G_upper, n = 1)
    
    #swap bounds if upper is greater than lower
    if (upper_limit_rnd < lower_limit_rnd) {
      
      bounds[[2]] <- lower_limit_rnd
      bounds[[1]] <- upper_limit_rnd
      
    } else {
      
      bounds[[1]] <- lower_limit_rnd
      bounds[[2]] <- upper_limit_rnd
    }
    
    ind[[production_genes[i]]] <- bounds
  }
  
  #Degradation
  for (i in 1:length(degradation_genes)) {
    
    lower_limit_rnd <- runif(K_lower, K_upper, n = 1)
    upper_limit_rnd <- runif(K_lower, K_upper, n = 1)
    
    if (upper_limit_rnd < lower_limit_rnd) {
      
      bounds[[2]] <- lower_limit_rnd
      bounds[[1]] <- upper_limit_rnd
      
    } else {
      
      bounds[[1]] <- lower_limit_rnd
      bounds[[2]] <- upper_limit_rnd
    }
    
    ind[[degradation_genes[i]]] <- bounds
  }
  
  
  #Threshold
  for (l in 1:length(threshold_genes)) {
    lower_limit_rnd <- runif(TH_lower, TH_upper, n = 1)
    upper_limit_rnd <- runif(TH_lower, TH_upper, n = 1)
    
    if (upper_limit_rnd < lower_limit_rnd) {
      
      bounds[[2]] <- lower_limit_rnd
      bounds[[1]] <- upper_limit_rnd
      
    } else {
      
      bounds[[1]] <- lower_limit_rnd
      bounds[[2]] <- upper_limit_rnd
    }
    ind[[threshold_genes[l]]] <- bounds
  }
  
  
  #Hill
  for (m in 1:length(hill_genes)) {
    lower_limit_rnd <- sample(N_lower:N_upper, size = 1, replace = TRUE)
    upper_limit_rnd <- sample(N_lower:N_upper, size = 1, replace = TRUE)
    
    if (upper_limit_rnd < lower_limit_rnd) {
      
      bounds[[2]] <- lower_limit_rnd
      bounds[[1]] <- upper_limit_rnd
      
    } else {
      
      bounds[[1]] <- lower_limit_rnd
      bounds[[2]] <- upper_limit_rnd
    }
    ind[[hill_genes[m]]] <- bounds
  }
  
  #FC
  for (n in 1:length(FC_genes)) {
    lower_limit_rnd <- runif(FC_lower, FC_upper, n = 1)
    upper_limit_rnd <- runif(FC_lower, FC_upper, n = 1)
    
    if (upper_limit_rnd < lower_limit_rnd) {
      
      bounds[[2]] <- lower_limit_rnd
      bounds[[1]] <- upper_limit_rnd
      
    } else {
      
      bounds[[1]] <- lower_limit_rnd
      bounds[[2]] <- upper_limit_rnd
    }
    
    ind[[FC_genes[n]]] <- bounds
  }
  
  ind
  
}

# a <- individual(sample_rset = reference_rset)
# unlist(a)
############################population##################################################################################################################################################################################
population <- function(sample_rset, population_size = 8, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  population_list <- vector(mode = "list", length = population_size)
  for (i in 1:population_size) {
    
    population_list[[i]] <- individual(sample_rset = sample_rset,  G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                                       N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
  }
  population_list
  
}
# a_pop <- population(sample_rset = reference_rset)
# a_pop
# 
# b_pop <- population(sample_rset = reference_rset_10000_0.02, population_size = 7)
#################selection#############################################################################################################################################################################################
selection <- function(generation, reference_individual_object, normalize_with_reference = TRUE, numModels = 10000, integrateStepSize = 0.02, scoring_method = "PCs", counter, generations_object) {
  
  #intialize list of fitness value objects
  fitness_value_objects_list <- vector(mode = "list", length = length(generation))
  
  #intialize vector of fitness values
  fitness_values_vec <- vector(length = length(fitness_value_objects_list))
  
  #to prevent redundant RACIPE analyses
  skip <- matrix(nrow = length(generation), ncol = length(generation))
  
  if (!missing(counter) & !missing(generations_object)) {
    
    if (counter > 2) {
      
      for (j in 1:length(generation)) {
        
        for (k in 1:length(generation)) {
          
          #are the genes of the kth individual in the previous generation identical to the genes of the jth individual of the current generation?
          skip[j,k] <- as.logical(identical(generation[[j]], generations_object[[counter-2]][[k]]$individual))
        }
      }
    }
    
  }
  
  
  if (!missing(counter) & !missing(generations_object)) {
    
    if (counter > 2) {
      
      for (i in 1:length(generation)) {
        
        if(sum(skip[i,]) > 0) {
          
          fitness_value_objects_list[[i]] <- generations_object[[counter-2]][[(which(skip[i,] == TRUE)[1])]]
          
          fitness_values_vec[[i]] <- generations_object[[counter-2]][[(which(skip[i,] == TRUE)[1])]]$fitness_values_sum
          
        } else {
          
          #compute fv object
          fitness_value_objects_list[[i]] <- generate_fitness_value_object(individual = generation[[i]], reference_individual_object = reference_individual_object,
                                                                           normalize_with_reference = normalize_with_reference, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize)
          
          #sums of fv scores of each individual
          fitness_values_vec[i] <- fitness_value_objects_list[[i]]$fitness_values_sum
          
        }
      }
    }
  } else {
    
    for (i in 1:length(generation)) {
      
      #compute fv object
      fitness_value_objects_list[[i]] <- generate_fitness_value_object(individual = generation[[i]], reference_individual_object = reference_individual_object,
                                                                       normalize_with_reference = normalize_with_reference, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize)
      
      #sums of fv scores of each individual
      fitness_values_vec[i] <- fitness_value_objects_list[[i]]$fitness_values_sum
      
    }
  }
  
  #order individuals in generation from lowest to highest scores
  fvs_ordered <- order(fitness_values_vec, decreasing = F)
  
  fitness_value_objects_list <- fitness_value_objects_list[c(fvs_ordered)]
  generation_ordered <- generation[c(fvs_ordered)]
  
  names(fitness_value_objects_list) <- paste("ind_", 1:length(generation))
  
  #depending on the chosen population size for all the generations, the individuals selected will be different. If the populations are a multiple of 4, or 4 + 1, we take half rounded down.
  #If 4 + 2 or 4 + 3, we take half and then round up after first adding a small number to the half.
  
  if (length(generation)%%4 == 0 | length(generation)%%4 == 1) {
    
    selected <- generation_ordered[1:floor(length(generation)/2)]
    
  } else {
    
    selected <- generation_ordered[1:ceiling(length(generation)/2 + 0.1)]
    
  }
  
  #create selections object to use later  
  selections_object <- list(selected = selected, generation_ordered = generation_ordered, fitness_value_objects_list = fitness_value_objects_list, skip_matrix_sums = rowSums(skip))
  
  selections_object
  
}

# fv_objects_list_test <- selection(a_pop, reference_individual_object = reference_individual_object_10000_0.02, integrateStepSize = 0.1, scoring_method = "Default")
# fv_objects_list_test$fitness_value_objects_list

# fv_objects_list_tes1t <- selection(b_pop, reference_individual_object = reference_individual_object_10000_0.02, integrateStepSize = 0.1, scoring_method = "Default")
# length(unlist(fv_objects_list_tes1t$selected))/2

# fv_objects_list_toy <- selection(a_pop, reference_individual_object = reference_individual_object_10000_0.02, integrateStepSize = 0.1, scoring_method = "Toy")
# fv_objects_list_toy

######################mating########################################################################################################################################################################################
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

# mating_test <- mating(a_pop[1:2])
# unlist(a_pop[[1]])
# unlist(mating_test[[1]])
# unlist(mating_test[[2]])
########################mutation######################################################################################################################################################################################
mutation <- function(individual, sample_rset, mutation_rate = 1, mutation_intensity = 0.1, mutation_method = "Random", G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  number_of_genes <- ncol(sracipeParams(sample_rset))
  parameter_names <-colnames(sracipeParams(sample_rset))
  
  mutated_individual <- individual
  
  production_genes <- grep(pattern = "G", parameter_names)
  degradation_genes <- grep(pattern = "K", parameter_names)
  threshold_genes <- grep(pattern = "TH", parameter_names)
  hill_genes <- grep(pattern = "N", parameter_names)
  FC_genes <- grep(pattern = "FC", parameter_names)
  
  mutation_indexes <- sample(1:length(individual), size = mutation_rate, replace = FALSE)
  
  if (mutation_method == "Random") {
    
    for (l in 1:mutation_rate) {
      
      #which bound will be changed, lower(0) or upper(1)
      bound <- round(runif(1, 0, 1))
      
      #which direction will be bound be changed, down (0) or up(1)
      direction <- round(runif(1, 0, 1))
      
      #random number between 0-1 to multiply mutation-size with
      random <- runif(1, 0, 1)
      
      #if the random index is one of the indexes for production parameters
      if (mutation_indexes[l] %in% production_genes) {
        
        #calculate the size of the mutation
        mutation <- mutation_intensity*(G_upper - G_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - mutation
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) - mutation
            
            #if temp lbound is less than global lbound
            if (new_bound < G_lower) {
              
              #new lbound becomes global lbound
              mutated_individual[[mutation_indexes[l]]][1] <- G_lower
              
            } else {
              
              #new lbound becomes temp lbound
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + mutation 
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) + mutation
            
            #if temp lbound > global lbound
            if (new_bound > G_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- G_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(individual[[mutation_indexes[l]]][2])) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][2])
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) - mutation #temp ubound is current ubound - mutation 
            
            if (new_bound < G_lower) {
              
              new_bound <- G_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(individual[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][1])
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) + mutation #if uound will be increased: temporary ubound is current ubound + mutation 
            
            if (new_bound > G_upper) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- G_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      
      #if the random index is one of the indexes for degradation parameters
      if (mutation_indexes[l] %in% degradation_genes) {
        
        #calculate the size of the mutation
        mutation <- mutation_intensity*(K_upper - K_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - mutation
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) - mutation
            
            #if temp lbound is less than global lbound
            if (new_bound < K_lower) {
              
              #new lbound becomes global lbound
              mutated_individual[[mutation_indexes[l]]][1] <- K_lower
              
            } else {
              
              #new lbound becomes temp lbound
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + mutation 
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) + mutation
            
            #if temp lbound > global lbound
            if (new_bound > K_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- K_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(individual[[mutation_indexes[l]]][2])) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][2])
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) - mutation #temp ubound is current ubound - mutation 
            
            if (new_bound < K_lower) {
              
              new_bound <- K_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(individual[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][1])
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) + mutation #if uound will be increased: temporary ubound is current ubound + mutation 
            
            if (new_bound > K_upper) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- K_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      
      #if the random index is one of the indexes for threshold parameters
      if (mutation_indexes[l] %in% threshold_genes) {
        
        #calculate the size of the mutation
        mutation <- mutation_intensity*(TH_upper - TH_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - mutation
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) - mutation
            
            #if temp lbound is less than global lbound
            if (new_bound < TH_lower) {
              
              #new lbound becomes global lbound
              mutated_individual[[mutation_indexes[l]]][1] <- TH_lower
              
            } else {
              
              #new lbound becomes temp lbound
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + mutation 
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) + mutation
            
            #if temp lbound > global lbound
            if (new_bound > TH_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- TH_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(individual[[mutation_indexes[l]]][2])) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][2])
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) - mutation #temp ubound is current ubound - mutation 
            
            if (new_bound < TH_lower) {
              
              new_bound <- TH_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(individual[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][1])
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) + mutation #if uound will be increased: temporary ubound is current ubound + mutation 
            
            if (new_bound > TH_upper) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- TH_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      #if the random index is one of the indexes for hill parameters
      if (mutation_indexes[l] %in% hill_genes) {
        
        if (bound == 0) {
          
          if (direction == 0) {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) - 1
            
            if (new_bound < N_lower) {
              
              new_bound <- N_lower
            }
            
            mutated_individual[[mutation_indexes[l]]][1] <- new_bound
            
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) + 1
            
            if (new_bound > N_upper) {
              
              new_bound <- N_upper
            }
            
            if (new_bound > as.numeric(individual[[mutation_indexes[l]]][2])) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- individual[[mutation_indexes[l]]][2]
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
            }
          }
        } else {
          
          if (direction == 0) {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) - 1
            
            if (new_bound < N_lower) {
              
              new_bound <- N_lower
            }
            
            if (new_bound < as.numeric(individual[[mutation_indexes[l]]][1])) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- individual[[mutation_indexes[l]]][1]
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
            }
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) + 1
            
            if (new_bound > N_upper) {
              
              new_bound <- N_upper
            }
            
            mutated_individual[[mutation_indexes[l]]][2] <- new_bound
            
          }
          
        }
      }
      
      #if the random index is one of the indexes for foldchange parameters
      if (mutation_indexes[l] %in% FC_genes) {
        
        #calculate the size of the mutation
        mutation <- mutation_intensity*(FC_upper - FC_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - mutation
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) - mutation
            
            #if temp lbound is less than global lbound
            if (new_bound < FC_lower) {
              
              #new lbound becomes global lbound
              mutated_individual[[mutation_indexes[l]]][1] <- FC_lower
              
            } else {
              
              #new lbound becomes temp lbound
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + mutation 
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][1]) + mutation
            
            #if temp lbound > global lbound
            if (new_bound > FC_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- FC_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(individual[[mutation_indexes[l]]][2])) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][2])
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) - mutation #temp ubound is current ubound - mutation 
            
            if (new_bound < FC_lower) {
              
              new_bound <- FC_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(individual[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][1])
              
              mutated_individual[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(individual[[mutation_indexes[l]]][2]) + mutation #if uound will be increased: temporary ubound is current ubound + mutation 
            
            if (new_bound > FC_upper) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- FC_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
    }    
    
    
  } else if (mutation_method == "Bounded") {
    
    for (l in 1:mutation_rate) {
      
      bound <- round(runif(1, 0, 1))
      
      if (mutation_indexes[l] %in% production_genes) {
        
        if (bound == 0) {
          
          mutated_individual[[mutation_indexes[l]]][1] <- runif(min = G_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
          
          
        } else {
          
          mutated_individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = G_upper, n = 1)
          
        }
      }
      
      
      if (mutation_indexes[l] %in% degradation_genes) {
        
        if (bound == 0) {
          
          mutated_individual[[mutation_indexes[l]]][1] <- runif(min = K_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
          
          
        } else {
          
          mutated_individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = K_upper, n = 1)
          
        }
        
      }
      
      
      if (mutation_indexes[l] %in% threshold_genes) {
        
        if (bound == 0) {
          
          mutated_individual[[mutation_indexes[l]]][1] <- runif(min = TH_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
          
          
        } else {
          
          mutated_individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = TH_upper, n = 1)
          
        }
        
      }
      
      if (mutation_indexes[l] %in% hill_genes) {
        
        direction <- round(runif(1, 0, 1))
        
        #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
        if ((bound == 0 | as.numeric(individual[[mutation_indexes[l]]][2]) == N_upper) & as.numeric(individual[[mutation_indexes[l]]][1]) != N_lower ) {
          
          #if the two bounds are equal to each other...lower bound must mutate down, OR if direction = 0
          if (as.numeric(individual[[mutation_indexes[l]]][1]) == as.numeric(individual[[mutation_indexes[l]]][2]) | direction == 0) {
            
            #new lower bound = previous lower - 1
            
            if ((as.numeric(individual[[mutation_indexes[l]]][1]) - 1) < N_lower) {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][1])
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][1]) - 1
            }
            
          } else {
            
            #new lower bound = previous lower + 1
            mutated_individual[[mutation_indexes[l]]][1] <- as.numeric(individual[[mutation_indexes[l]]][1]) + 1
          }
          
        } else {
          
          #if the two bounds are equal to each upper bound must mutate up, OR if direction = 1
          if (as.numeric(individual[[mutation_indexes[l]]][2]) == as.numeric(individual[[mutation_indexes[l]]][2]) | direction == 1) {
            
            #new upper bound = previous upper + 1
            
            if ((as.numeric(individual[[mutation_indexes[l]]][2]) + 1) > N_upper) {
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][2])
              
            } else {
              
              mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][2]) + 1
            }
            
          } else {
            
            #new upper bound = previous upper - (previous_upper - 1
            mutated_individual[[mutation_indexes[l]]][2] <- as.numeric(individual[[mutation_indexes[l]]][2]) - 1
          }
        }
      }
      
      if (mutation_indexes[l] %in% FC_genes) {
        
        if (bound == 0) {
          
          mutated_individual[[mutation_indexes[l]]][1] <- runif(min = FC_lower, max = as.numeric(individual[[mutation_indexes[l]]][2]), n = 1)
          
          
        } else {
          
          mutated_individual[[mutation_indexes[l]]][2] <- runif(min = as.numeric(individual[[mutation_indexes[l]]][1]), max = FC_upper, n = 1)
          
        }
        
      }
      
    }
  }
  
  mutated_individual
  
}

# unlist(a)
# mutation_test_random <- mutation(individual = a, sample_rset = reference_rset_10000_0.02, mutation_rate = 10, mutation_intensity = 0.5, mutation_method = "Random")
# unlist(mutation_test_random)
# mutation_test_bounded <- mutation(individual = a, sample_rset = reference_rset_10000_0.02, mutation_rate = 10, mutation_intensity = 0.5, mutation_method = "Bounded")
# unlist(mutation_test_bounded)
#############################reproduction#################################################################################################################################################################################
reproduction <- function(selections_object, mutation_rate = 1, mutation_intensity = 0.1, mating_method = "Random", mutation_method = "Random", G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  #make temp objects for downstream use/readability of code
  population_size <- length(selections_object$fitness_value_objects_list)
  
  selected <- selections_object$selected
  
  #order selected individuals based on the mating method
  if (mating_method == "Random") {
    
    random <- sample(1:length(selected), size = length(selected), replace = FALSE)
    selected_for_mating <- selected[random]
    
  } else if (mating_method == "Best") {
    
    selected_for_mating <- selected
    
  }
  
  sample_rset <- selections_object$fitness_value_objects_list[[1]]$simulated_rset
  
  #create sequence of every other index from 1 to the length of selected (the number of selected individuals for mating)
  pairs <- seq(1, length(selected), 2)
  
  #initilialize list of children 
  children <- list()
  
  #mate
  for (i in 1:length(pairs)) {
    
    children <- append(children, mating(selected_for_mating[pairs[i]:(pairs[i] + 1)]))
    # children <- do.call(c, list(children, mating(selected[pairs[[i]]:(pairs[[i]] + 1)])))
    
  }
  
  #mutate children
  mutated <- vector(mode = "list", length = length(children)) #initialize list
  
  for (j in 1:length(mutated)) { 
    
    mutated[[j]] <- mutation(children[[j]], sample_rset = sample_rset,
                             mutation_rate = mutation_rate, mutation_intensity = mutation_intensity, mutation_method = mutation_method,
                             G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                             N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
  }
  
  
  
  #new generation depends on population size. In the cases of odd population sizes, a totally new individual is generated
  if (population_size %% 4 == 0) {
    
    new_generation <- append(selected, mutated)
    
  } else if (population_size %% 4 == 1) {
    
    new_generation <- append(selected, mutated)
    
    new_individual <- individual(sample_rset = sample_rset, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                                 N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
    new_generation <- append(new_generation, list(new_individual))
    
  } else if (population_size %% 4 == 2) {
    
    new_generation <- append(selected[1:(length(selected)-2)], mutated)
    
  } else if (population_size %% 4 == 3) {
    
    new_generation <- append(selected[1:(length(selected)-2)], mutated)
    
    new_individual <- individual(sample_rset = sample_rset, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                                 N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
    new_generation <- append(new_generation, list(new_individual))
    
  }
  
  new_generation
  
}

##########

####################GA_R_RACIPE##########################################################################################################################################################################################
GA_R_RACIPE <- function(reference_individual_object, population_size = 16, max_generations = 10, numModels = 10000, integrateStepSize = 0.02,  mutation_rate = 1, mutation_intensity = 0.25, mutation_method = "Random",
                        scoring_method = "Default", normalize_with_reference = TRUE, min_score = "NA", saving_frequency = 2, mating_method = "Random",
                        G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100 ) {
  
  
  #make temp objects for downstream use/readability of code
  sample_rset <- reference_individual_object$reference_rset
  
  #initialize results lists/vectors
  generations_object <- vector(mode = "list", length = max_generations)
  names(generations_object) <- paste0("generation_", 1:max_generations)
  
  generations_object_saved <- vector(mode = "list", length = max_generations)
  names(generations_object_saved) <- paste0("generation_", 1:max_generations)
  
  generations_list <- vector(mode = "list", length = max_generations)
  names(generations_list) <- paste0("generation_", 1:max_generations)
  
  top_fv_vec <- vector(length = max_generations)
  names(top_fv_vec) <- paste0("generation_", 1:max_generations)
  
  skips_object <- vector(mode = "list", length = max_generations)
  names(skips_object) <- paste0("generation_", 1:max_generations)
  
  #generate intial population
  generation <- population(sample_rset = sample_rset, population_size = population_size, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                           K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
  
  #Analyze 1st two generations out of loop
  #1
  
  selections_object <- selection(generation = generation, reference_individual_object = reference_individual_object, 
                                 normalize_with_reference = normalize_with_reference, numModels = numModels, integrateStepSize = integrateStepSize, scoring_method = scoring_method)
  
  
  generations_object[[1]] <- selections_object$fitness_value_objects_list
  generations_object_saved[[1]] <- selections_object$fitness_value_objects_list
  generations_list[[1]] <- selections_object$generation_ordered
  
  top_fv_vec[1] <- selections_object$fitness_value_objects_list[[1]]$fitness_values_sum
  skips_object[[1]] <- selections_object$skip_matrix_sums
  
  generation <- reproduction(selections_object = selections_object, mutation_rate = mutation_rate, mutation_intensity = mutation_intensity, mutation_method = mutation_method, mating_method = mating_method, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                             K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
  
  #2
  selections_object <- selection(generation = generation, reference_individual_object = reference_individual_object, 
                                 normalize_with_reference = normalize_with_reference, numModels = numModels, integrateStepSize = integrateStepSize, scoring_method = scoring_method)
  
  
  generations_object[[2]] <- selections_object$fitness_value_objects_list
  generations_object_saved[[2]] <- selections_object$fitness_value_objects_list
  generations_list[[2]] <- selections_object$generation_ordered
  
  top_fv_vec[2] <- selections_object$fitness_value_objects_list[[1]]$fitness_values_sum
  skips_object[[2]] <- selections_object$skip_matrix_sums
  
  generation <- reproduction(selections_object = selections_object, mutation_rate = mutation_rate, mutation_intensity = mutation_intensity, mutation_method = mutation_method, mating_method = mating_method, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                             K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
  
  
  i <- 3
  
  while (i <= max_generations) {
    
    counter <- i
    
    selections_object <- selection(generation = generation, reference_individual_object = reference_individual_object, 
                                   normalize_with_reference = normalize_with_reference, numModels = numModels, integrateStepSize = integrateStepSize, scoring_method = scoring_method, counter = counter, generations_object = generations_object)
    
    generations_object[[i]] <- selections_object$fitness_value_objects_list
    
    #saving generations object
    if (counter%%saving_frequency == 0 | counter == max_generations) {
      
      generations_object_saved[[i]] <- selections_object$fitness_value_objects_list
      
    }
    
    generations_list[[i]] <- selections_object$generation_ordered
    
    top_fv_vec[i] <- selections_object$fitness_value_objects_list[[1]]$fitness_values_sum
    
    
    skips_object[[i]] <- selections_object$skip_matrix_sums
    
    generation <- reproduction(selections_object = selections_object, mutation_rate = mutation_rate, mutation_intensity = mutation_intensity, mutation_method = mutation_method, mating_method = mating_method, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower,
                               K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
    
    #break loop if min_score exists and if it is greater than the current score
    if(!missing(min_score) & is.numeric(min_score)) {
      
      if (top_fv_vec[i] <= min_score) {
        
        generations_object_saved[[i]] <- selections_object$fitness_value_objects_list
        
        break 
      }
      
    }
    
    #increment generation number
    i <- i + 1
    
  }
  top_fv_vec <- top_fv_vec[1:(i-1)]
  generations_object_saved <- plyr::compact(generations_object_saved)
  generations_list <- plyr::compact(generations_list)
  skips_object <- plyr::compact(skips_object)
  
  GA_results <- list(top_fv_vec = top_fv_vec, generations_list = generations_list, generations_object = generations_object_saved, skips_object = skips_object)
  
}

#####################Testing#########################################################################################################################################################################################

# GA_R_RACIPE_toy <- GA_R_RACIPE(reference_individual_object = reference_individual_object_10000_0.02, population_size = 4, max_generations = 1000, mutation_rate = 4, mutation_intensity = 0.25, mutation_method = "Random",
#                                min_score = -505, scoring_method = "Toy", normalize_with_reference = FALSE, numModels = 2000, integrateStepSize = 0.02)
# 
# 
# GA_R_RACIPE_test <- GA_R_RACIPE(reference_individual_object = reference_individual_object_10000_0.02, population_size = 8, max_generations = 10, mutation_rate = 2, mutation_intensity = 0.6, mutation_method = "Random", mating_method = "Random", saving_frequency = 2,
#                                  min_score = "NA", scoring_method = "Default", normalize_with_reference = TRUE, numModels = 10000, integrateStepSize = 0.02)


