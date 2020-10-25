#####################Globals##################################################################################################################################################################
current_state_og <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(2,4), list(2,4), list(20,50), list(20,50))
current_state_og <- split(unlist(current_state_og), ceiling(seq_along(unlist(current_state_og))/2))

current_state_og2 <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(1.95,4), list(2,4.05), list(20,50), list(20,50))

#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50), list(0.2,0.5), list(0.2, 0.5), list(10,40), list(10,40), list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("N")
#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50), list(10,40), list(10,40), list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("K","N")
#######################################################################################################################################################################################
current_state <- list(list(20,50), list(20,50))
state_vector <- unlist(current_state)
skip <- c("G","K","TH", "N")


#######################################################################################################################################################################################
reference_state <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
numModels <- 10000
integrateStepSize <- 0.02
reference <- generate_ref_object(topology = topology_TS, reference_state = reference_state, numModels = numModels, integrateStepSize = integrateStepSize)

global_bounds <- list(G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100)
g_bounds <- list(G = c(1, 100), K = c(0.1, 1), TH = c(0, 55), N = c(1, 5), FC = c(1, 100))
scoring_method <- "Default"
step_size <- 0.02
gv_numModels_FC <- 2

iteration_pdf <- 0
iteration_cf <- 0

cost_function_input <- list()
partial_derivs_input <- list()

number_of_iterations_gv <- 4

gv_current_state_scores_object <- c()

cf_output <- list()


cost_function_wrapper_GA(state_vector)
test <- calc_score(state = list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(2,100), list(2,100)), reference_rset_object = reference)
test$cf_vec_sum
#####################cost function wrapper############################################################################################################################################################################################################################
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
    
    costs_object <- cost_function_GA(state = state_list, reference_rset_object = reference, global_bounds = g_bounds, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = FALSE, skipped_params = skip)
    
    output_vec[k] <- costs_object$cf_vec_sum
    
  }

  cf_output <<- append(cf_output, list(output_vec))
  
  output <- -1*mean(output_vec)
  
  }
  
  print(output)
  
  output
  
}

#####################cost_function############################################################################################################################################################################################################################
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



#####################ga testing############################################################################################################################################################################################################################


GA <- ga(type = "real-valued", 
         fitness =  cost_function_wrapper_GA,
         lower = c(1,1, 1, 1), upper = c(100,100, 100, 100), 
         popSize = 50, maxiter = 10)

GA@solution

#####################hessian############################################################################################################################################################################################################################
hessian <- function(state, reference_rset_object, scoring_method = "Default", numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE) {
  
  
  #######Pre-Processing######
  ##########################################transform state_vector into list
  state_current <- state
  
  state_reference <- split(unlist(reference_rset_object$reference_state), ceiling(seq_along(unlist(reference_rset_object$reference_state))/2))
  
  ##########################################extract parameter names and identify index of each parameter type
  parameter_names <- colnames(sracipeParams(reference_rset_object$reference_rset))
  
  parameter_bounds_names <- paste(rep(parameter_names, each = 2),
                                  rep(c("lower", "upper"), times = length(parameter_names)), sep = "_")
  
  production_params <- grep(pattern = "G", parameter_names)
  
  degradation_params <- grep(pattern = "K", parameter_names)
  
  threshold_params <- grep(pattern = "TH", parameter_names)
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  FC_params <- grep(pattern = "FC", parameter_names)
  
  params_list <- list(production_params = production_params, degradation_params = degradation_params, threshold_params = threshold_params, hill_params = hill_params, FC_params = FC_params)
  
  ##########################################index of parameters which should be skipped from full parameter list
  if (!missing(skipped_params)) {
    
    params_skipped <- c()
    bounds_skipped <- c()
    
    for (i in 1:length(skipped_params)) {
      
      params_skipped <- c(params_skipped, grep(pattern = skipped_params[i], parameter_names) )
      
    }
    
    for (i in 1:length(skipped_params)) {
      
      bounds_skipped <- c(bounds_skipped, grep(pattern = skipped_params[i], parameter_bounds_names) )
      
    }
    
  }
  
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
  
  # if (!missing(skipped_params)) {
  #   
  #   global_bounds_per_parameter <- global_bounds_per_parameter[-params_skipped]
  #   
  # }
  
  
  ##########################################check if current state bounds are within global bounds for each paramter
  
  for (i in 1:length(state_current)) {
    
    if (state_current[[i]][[1]] < global_bounds_per_parameter[[i]][[1]]) {
      
      state_current[[i]][[1]] <- global_bounds_per_parameter[[i]][[1]]
      
    }
    
    if (state_current[[i]][[2]] > global_bounds_per_parameter[[i]][[2]]) {
      
      state_current[[i]][[2]] <- global_bounds_per_parameter[[i]][[2]]
      
    }
    
  }
  
  #######Generate partial states list######
  
  ##########################################number of models to simulate to add
  numModels_partial <- numModels*step_size*numModels_FC
  
  ##########################################generate vector of the step magnitude for each parameter
  step_magnitudes <- vector(length = length(parameter_names))
  
  counter <- 0
  
  for (i in 1:length(step_magnitudes)) {
    
    step_magnitudes[i] <- (state_current[[i]][[2]] - state_current[[i]][[1]])*step_size
    
  }
  
  step_magnitudes_rep <- rep(step_magnitudes, each = 2)
  
  ##########################################generate list of partial states
  partial_states_list <- vector(mode = "list", length = length(parameter_names)*2) #*2 because subtracting from lower bound and adding to upper bound
  
  counter <- 0
  
  for (j in 1:length(state_current)) {
    
    #first lower/upper refers to which bound, second minus/plus refers t
    
    partial_state_lower_bound <- state_current
    partial_state_upper_bound <- state_current
    
    partial_state_lower_bound[[j]][[1]] <- state_current[[j]][[1]] - step_magnitudes[j]
    partial_state_lower_bound[[j]][[2]] <- state_current[[j]][[1]]
    
    partial_state_upper_bound[[j]][[1]] <- state_current[[j]][[2]]
    partial_state_upper_bound[[j]][[2]] <- state_current[[j]][[2]] + step_magnitudes[j]
    
    
    partial_states_list[[j + counter]] <- partial_state_lower_bound
    partial_states_list[[j + counter + 1]] <- partial_state_upper_bound
    
    counter <- counter + 1
    
  }
  
  #######Simulate partial states and Calculate Partial Derivatives######
  ##########################################Simulate current state
  
  # current_state_scores_object <- calc_score(state = state_current, reference_rset_object = reference_rset_object, scoring_method = scoring_method, numModels = numModels, integrateStepSize = integrateStepSize, save_rset = TRUE)
  
  ##########################################Duplicate current state expression matrix *FC* times
  current_state_expr_rep <- do.call(rbind, replicate(numModels_FC, current_state_scores_object$simulated_expression, simplify = FALSE))
  
  
  
  
  
  
}
