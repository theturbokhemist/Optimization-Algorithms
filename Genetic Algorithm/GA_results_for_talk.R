##############################generate_ref_object########################################################################################################################################################
generate_ref_object <- function(topology, reference_state, numModels = 10000, integrateStepSize = 0.02, global_parameter) {
  
  reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  state <- reference_state
  
  #hill gene indexes
  hill_genes <- grep(pattern = "N", colnames(sracipeParams(reference_rset)))
  
  if (!missing(global_parameter)) {
    
    parametric_variation_reference_state <- reference_state
    
    for (w in 1:length(parametric_variation_reference_state)) {
      
      new_bounds <- parametric_variation_bounds(parameter_bounds = reference_state[[w]], global_parameter = global_parameter)
      
      parametric_variation_reference_state[[w]][[1]] <- new_bounds$bound_min
      parametric_variation_reference_state[[w]][[2]] <- new_bounds$bound_max
      
    }
    
    for (u in 1:length(hill_genes)) {
      
      parametric_variation_reference_state[[hill_genes[u]]][[1]] <- round(parametric_variation_reference_state[[hill_genes[u]]][[1]])
      parametric_variation_reference_state[[hill_genes[u]]][[2]] <- round(parametric_variation_reference_state[[hill_genes[u]]][[2]]) 
      
    }
    
    state <- parametric_variation_reference_state
    
  }
  
  
  #assign vaues from uniform samples distributions to each model of reference rset.  
  for (i in 1:length(state)) {
    
    sracipeParams(reference_rset)[,i] <- runif(numModels, min= as.numeric(state[[i]][1]), max= as.numeric(state[[i]][2]))  
    
  }
  
  for (k in 1:length(hill_genes)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(reference_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(state[[hill_genes[k]]][1]) - 0.5, max= as.numeric(state[[hill_genes[k]]][2]) + 0.5))  
    
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
  reference_rset_object <- list(reference_rset = ref_expression_rset, reference_state = reference_state, expression_ref = expr_reference_log_z, prcomp_ref = prcomp_reference, means_ref = means_reference, sds_ref = sds_reference, topology = topology)
  
  if (!missing(global_parameter)) {
    
    reference_rset_object <- list(reference_rset = ref_expression_rset, reference_state = reference_state, expression_ref = expr_reference_log_z, prcomp_ref = prcomp_reference, means_ref = means_reference, sds_ref = sds_reference, topology = topology,
                                  parametric_variation_reference_state = parametric_variation_reference_state, global_parameter = global_parameter)
    
  }
  
  reference_rset_object
  
}

########################################################




cost_function_GA <- function(state, reference_rset_object, global_bounds, scoring_method = "Default", numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE, skipped_params) {
  
  #######Pre-Processing######
  ##########################################transform state_vector into list
  state_current <- state
  
  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)
  
  ##########################################append reference values for skipped parameters
  if (!missing(skipped_params)) {
    
    if (sum(skipped_params %in% c("G")) > 0) {
      
      state_current <- append(state_current, state_reference[production_params], after = production_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("K")) > 0) {
      
      state_current <- append(state_current, state_reference[degradation_params], after = degradation_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("TH")) > 0) {
      
      state_current <- append(state_current, state_reference[threshold_params], after = threshold_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("N")) > 0) {
      
      state_current <- append(state_current, state_reference[hill_params], after = hill_params[1] - 1)   
      
    }
    
    if (sum(skipped_params %in% c("FC")) > 0) {
      
      state_current <- append(state_current, state_reference[FC_params], after = FC_params[1] - 1)   
      
    }
  }
  
  names(state_current) <- parameter_names
  
  ##########################################global bounds for each parameter
  global_bounds_per_parameter <- list()
  
  for (i in 1:length(global_bounds)) {
    
    global_bounds_per_parameter <- append(global_bounds_per_parameter, (rep(global_bounds[i], length(params_list[[i]]))))
    
  }
  
  names(global_bounds_per_parameter) <- parameter_names
  
  ##########################################check if current state bounds are within global bounds for each paramter
  
  for (i in 1:length(state_current)) {
    
    if (state_current[[i]][[1]] < global_bounds_per_parameter[[i]][[1]]) {
      
      state_current[[i]][[1]] <- global_bounds_per_parameter[[i]][[1]]
      
    }
    
    if (state_current[[i]][[2]] > global_bounds_per_parameter[[i]][[2]]) {
      
      state_current[[i]][[2]] <- global_bounds_per_parameter[[i]][[2]]
      
    }
    
  }
  
  #######Pre-Simulation and Simulation######
  ##########################################create rset object
  state <- state_current
  
  simulated_rset <- sracipeSimulate(reference_rset_object$topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  ##########################################sample parameter ranges from current state *numModels* times, add those results into the rset objects parameter values matrix
  for (t in 1:length(state)) {
    
    sracipeParams(simulated_rset)[,t] <- runif(numModels, min = state[[t]][[1]], max = state[[t]][[2]])  
    
  }
  
  #sample and round all the hill coefficients
  for (u in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(simulated_rset)[,hill_params[u]] <- round(runif(numModels, min= round(state[[hill_params[u]]][[1]]) - 0.5, max= round(state[[hill_params[u]]][[2]]) + 0.5))
    
  }
  
  ##########################################simulate
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
  
  #######Calculations######
  
  ##########################################log transform
  simulated_rset_expr_log <- log2(t(assay(simulated_rset)))
  
  ##########################################z-score normalize
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log,
                                     2, reference_rset_object$means_ref, FUN = "-")
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log_z,
                                     2, reference_rset_object$sds_ref, FUN = "/")
  
  ##########################################rotate (transformation)
  simulated_rset_expr_PCs <- simulated_rset_expr_log_z %*% reference_rset_object$prcomp_ref$rotation
  
  #######Scoring######
  
  #cost_function vector
  cf_vec <- vector(length = ncol(simulated_rset_expr_log_z))
  
  ##########################################kolmogorov-smirnov
  if (scoring_method == "PCs") {
    
    simulated_expression <- simulated_rset_expr_PCs
    
    #ks score for distribution of each gene
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,w], y = simulated_rset_expr_PCs[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
    
  } else {
    
    simulated_expression <- simulated_rset_expr_log_z
    
    for (w in 1:length(cf_vec)) {
      
      ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,w], y = simulated_rset_expr_log_z[,w])
      cf_vec[w] <- as.numeric(ks_test$statistic)
      
    }
  }
  
  costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec))
  
  if (save_rset == TRUE) {
    
    costs_object <- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset, simulated_expression = simulated_expression,
                         arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state, params_list = params_list, reference_state = reference_state, reference_rset_object = reference_rset_object))
    
  }
  
  gv_current_state_scores_object <<- list(cf_vec = cf_vec, cf_vec_sum = sum(cf_vec), simulated_rset = simulated_rset, simulated_expression = simulated_expression,
                                          arguments = list(scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, state = state, params_list = params_list, reference_state = reference_state, reference_rset_object = reference_rset_object))
  
  costs_object
  
  
}


cost_function_wrapper(vector = state_vector)




####################################
cost_function_wrapper_GA <- function(vector) {
  
  print(vector)
  
  differences_vector <- vector(length = length(vector)/2)
  
  counter <- 0
  
  for (i in 1:length(differences_vector)) {
    
    if (vector[i + counter] > vector[i + counter + 1]) {
      
      differences_vector[i] <- vector[i + counter + 1] - vector[i + counter] 
      
    } else {
      
      differences_vector[i] <- 0
      
    }
    
    counter <- counter + 1
    
  }
  
  print(differences_vector)
  
  if (sum(differences_vector) < 0) {
    
    output <- sum(differences_vector)
    
  } else {
    
    
    iteration_cf <<- iteration_cf + 1
    print(paste("Cost Function Iteration: ", iteration_cf))
    
    state_list <- split(vector, ceiling(seq_along(vector)/2))
    
    cost_function_input <<- append(cost_function_input, list(state_list))
    
    output_vec <- vector(length = number_of_iterations_gv)
    
    for (k in 1:number_of_iterations_gv) {
      
      costs_object <- cost_function_GA(state = state_list, reference_rset_object = reference_object, global_bounds = g_bounds, scoring_method = scoring_method, numModels = numModels,
                                       integrateStepSize = integrateStepSize, save_rset = FALSE, skipped_params = skip)
      
      output_vec[k] <- costs_object$cf_vec_sum
      
    }
    
    cf_output <<- append(cf_output, list(output_vec))
    
    output <- -1*mean(output_vec)
    
  }
  
  print(output)
  
  output
  
}

##############
reference_object <- generate_ref_object(topology = topology_TS, reference_state = list(c(1,10), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)),
                                 numModels = 10000, integrateStepSize = 0.02)


skip <- c("K","TH", "N", "FC")
global_bounds <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)
numModels <- 5000
integrateStepSize <- 0.02
g_bounds <- list(G = c(1, 100), K = c(0.1, 1), TH = c(0, 55), N = c(1, 5), FC = c(1, 100))
scoring_method <- "Default"

iteration_cf <- 0

cost_function_input <- list()

number_of_iterations_gv <- 4

cf_output <- list()


a <- cost_function_GA(state = list(c(1,100), c(1,100)), reference_rset_object = reference_object, global_bounds = g_bounds, scoring_method = scoring_method,
                 numModels = 5000, integrateStepSize = integrateStepSize, save_rset = TRUE, skipped_params = skip)

a$cf_vec_sum

b <- cost_function_wrapper_GA(c(1,100,1,100))


GA_final <- ga(type = "real-valued", 
         fitness =  cost_function_wrapper_GA,
         lower = c(1,1, 1, 1), upper = c(100,100, 100, 100), 
         popSize = 50, maxiter = 100)

GA_final@solution
plot(GA_final)
GA_final@population
##########

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))

standard_state_TS <- list(c(1,100), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #standard

standard_object_TS <- generate_standard_object(topology = topology_TS, standard_state = standard_state_TS, numModels = 10000)

standard_expr_A <- standard_object_TS$expression_ref[, "A"]

standard_expression_TS.kmeans <- kmeans(standard_object_TS$expression_ref, centers = 2, nstart = 25)


experimental_state_TS.G_A.l <- list(c(1,10), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #experimental

experimental_object_TS.G_A.l <- generate_experimental_object(experimental_state = experimental_state_TS.G_A.l, standard_rset_object = standard_object_TS, numModels = 1000,
                                                             integrateStepSize = 0.02, save_rset = TRUE)
experimental_expr_A <- experimental_object_TS.G_A.l$simulated_rset_expr_log_z[, "A"]
#################
