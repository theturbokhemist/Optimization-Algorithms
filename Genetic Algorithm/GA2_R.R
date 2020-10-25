topology_TS4 <- rbind(topology_TS3, data.frame(Source = "D", Target = "B", Type = 1))

simulated_rset4 <- sracipeSimulate(topology_TS4, integrate = FALSE, numModels = 100, genParams = TRUE)
sracipeConfig(simulated_rset4)
colnames(sracipeParams(simulated_rset4))
names(simulated_rset4)
grep(pattern = "G", names)
names <- colnames(sracipeParams(simulated_rset4))
names

sracipeConfig(simulated_rset)
names(simulated_rset)
individual(topology_TS)

individual <- function(topology, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  simulated_rset <- sracipeSimulate(topology, integrate = FALSE, numModels = 100, genParams = TRUE)
  number_of_genes <- ncol(sracipeParams(simulated_rset))
  parameter_names <-colnames(sracipeParams(simulated_rset))

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
  
  thresholds <- sracipeConfig(simulated_rset)$thresholds
  threshold_names <- names(simulated_rset)

  
  threshold_values <- vector(length = length(threshold_genes))
  
  for (k in 1:length(names(simulated_rset))) {
    
    names_index <- grep(paste0("TH_", threshold_names[k]), parameter_names[threshold_genes])
    threshold_values[c(names_index)] <- thresholds[k]
    
  }

  threshold_df <- data.frame(interactions = c(parameter_names[threshold_genes]), medians = threshold_values)
  
  
  for (l in 1:length(threshold_genes)) {
    lower_limit_rnd <- runif(TH_lower, TH_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (TH_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd*threshold_df[l, 2]
    bounds[[2]] <- upper_limit_rnd*threshold_df[l, 2]
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

testpop <- population(topology_TS4)
testpop[[2]]

population <- function(topology, population_size = 8, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 6, FC_lower = 1, FC_upper = 100) {
  
  population_list <- vector(mode = "list", length = population_size)
  for (i in 1:population_size) {
    
    population_list[[i]] <- individual(topology, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper,
                                       N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
  }
  population_list
  
}

fv <- function(individual, reference_expression, topology, numModels = 2000) {
  
  simulated_rset <- sracipeSimulate(topology, integrate = FALSE, numModels = numModels, genParams = TRUE)
  
  for (i in 1:length(individual)) {
    
    
    
  }
  
  sracipeParams(simulated_rset)[,1] <- runif(numModels, min= as.numeric(individual[[1]][1]), max= as.numeric(individual[[1]][2]))
  
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, numModels = numModels, genParams = FALSE)
  
  simulated_expression <- log2(t(assay(simulated_rset)))
  
  simulated_expression_z <- sweep(simulated_expression, 
                                  2, colMeans(simulated_expression), FUN = "-")
  
  simulated_expression_z <- sweep(simulated_expression_z, 
                                  2, apply(simulated_expression, 2, sd), FUN = "/")
  
  # cdf <- ggplot2::ggplot(data = data.frame(expression = c(rSet_TS_2000_GA_25max_expression_log_z[,1], rSet_TS_2000_GA_25max_expression_log_z[,2])), aes(x = expression)) + stat_ecdf(geom = "step", color = "red") +
  #   stat_ecdf(data = data.frame(expression = c(simulated_expression_z[,1], simulated_expression_z[,2])),color = "green", geom = "step")
  # print(cdf)
  
  fv_vec <-  vector(length = ncol(reference_expression))
  
  for (j in 1:ncol(reference_expression)) {
    
    ks_test <- dgof::ks.test(x = reference_expression[,j], y = simulated_expression_z[,j])
    fv_vec[j] <- as.numeric(1000- 1000*ks_test$statistic)
    
  }
  
  
  fv <- sum(fv_vec)
  fv 
  
}

fv(b, reference_expression = rSet_TS_2000_GA_25max_expression_log_z, topology = topology_TS)

a = population(number_of_genes = 1, population_size = 8)
length(a)
length(a[[1]])
b = individual(number_of_genes = 1)
b[[1]][1]


selection <- function(generation, reference_expression, topology) {
  fitness_values <- vector(length = length(generation))
  
  for (i in 1:length(generation)) {
    
    fitness_values[i] <- fv(generation[[i]], reference_expression = reference_expression, topology = topology)
    
  }
  fv_index_ordered <- order(fitness_values, decreasing = T)
  fitness_values <<- fitness_values[c(fv_index_ordered)]
  generation <- generation[c(fv_index_ordered)]
  selected <- generation[1:(length(generation)/2)]
  selected
}
t <- selection(a, reference_expression = rSet_TS_2000_GA_25max_expression_log_z, topology = topology_TS)

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

mutation <- function(individual, lower_limit, upper_limit, mutation_rate = 1) {
  
  mutation_indexes <- sample(1:length(individual), mutation_rate)
  
  for (l in 1:mutation_rate) {
    
    bound <- round(runif(1, 0, 1))
    if (bound == 0) {
      
      individual[[mutation_indexes[l]]][1] <- sample(lower_limit:as.numeric(individual[[mutation_indexes[l]]][2]), 1)
      
    } else {
      
      individual[[mutation_indexes[l]]][2] <- sample(as.numeric(individual[[mutation_indexes[l]]][1]):upper_limit, 1)
      
    }
    
  }
  
  individual
  
}

reproduction <- function(selected, lower_limit, upper_limit, mutation_rate = 1) {
  
  pairs <- seq(1, length(selected), 2)
  unmutated <- list()
  mutated <- vector(mode = "list", length = length(selected))
  
  for (i in 1:length(pairs)) {
    
    unmutated <- append(unmutated, mating(selected[pairs[[i]]:(pairs[[i]] + 1)]))
    # unmutated <- do.call(c, list(unmutated, mating(selected[pairs[[i]]:(pairs[[i]] + 1)])))
    
  }
  for (j in 1:length(unmutated)) {
    
    mutated[[j]] <- mutation(unmutated[[j]], lower_limit = lower_limit, upper_limit = upper_limit, mutation_rate = mutation_rate)
    
  }
  
  new_generation <- append(selected, mutated)
  new_generation
  
}

GA_R_RACIPE <- function(topology, reference_expression, population_size = 4, number_of_genes = 1, lower_limit = 1, upper_limit = 100, number_of_generations = 10, numModels = 2000, mutation_rate = 1 ) {
  
  generations_list <<- vector(mode = "list", length = number_of_generations)
  fvs_list <<- vector(mode = "list", length = number_of_generations)
  
  generation <- population(population_size = population_size, number_of_genes = number_of_genes, lower_limit = lower_limit, upper_limit = upper_limit)
  
  for (i in 1:number_of_generations) {
    print(paste("Generation:", i))
    generations_list[[i]] <<- generation
    selected <- selection(generation = generation, reference_expression = reference_expression, topology = topology)
    fvs_list[[i]] <<- fitness_values
    
    
    generation <- reproduction(selected = selected, lower_limit = lower_limit, upper_limit = upper_limit, mutation_rate = mutation_rate)
    generation
    
  }
  
  
  
  
}



h <- list()
g <- list("a", "b")

append(h, g)
l
h[1:2]
for (k in 1:length(h)) {
  print(h[[k]])
  k <- k + 1
  
}
seq(1, 50, 2)
