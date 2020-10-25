##############################starting_point########################################################################################################################################################
starting_point <- function(sample_rset, G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0, TH_upper = 55, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  #number of parameters
  number_of_params <- ncol(sracipeParams(sample_rset))
  
  #vector of paramater names
  parameter_names <-colnames(sracipeParams(sample_rset))
  
  #list that is length of number of params
  state <- vector(mode = "list", length = number_of_params)
  names(state) <- parameter_names
  
  #index of each parameter type in vector of parameter names
  production_params <- grep(pattern = "G", parameter_names)
  degradation_params <- grep(pattern = "K", parameter_names)
  threshold_params <- grep(pattern = "TH", parameter_names)
  hill_params <- grep(pattern = "N", parameter_names)
  FC_params <- grep(pattern = "FC", parameter_names)
  
  #list that will contain lower and upper bound for each parameter
  bounds <- vector(mode = "list", length = 2)
  
  #Randomize lower and upper bounds for each Production parameter
  for (i in 1:length(production_params)) {
    
    #sample uniformly from absolute bounds to get lower bound
    lower_limit_rnd <- runif(G_lower, G_upper, n = 1)
    
    #sample uniformly from lower bound and absolute upper bound to get upper bound
    upper_limit_rnd <- lower_limit_rnd + (G_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    
    #assign bounds of this specific parameter
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    
    state[[production_params[i]]] <- bounds
  }
  
  #Degradation
  for (j in 1:length(degradation_params)) {
    lower_limit_rnd <- runif(K_lower, K_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (K_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    state[[degradation_params[j]]] <- bounds
  }
  
  #Threshold
  for (l in 1:length(threshold_params)) {
    lower_limit_rnd <- runif(TH_lower, TH_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (TH_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    state[[threshold_params[l]]] <- bounds
  }
  
  
  #Hill
  for (m in 1:length(hill_params)) {
    lower_limit_rnd <- sample(N_lower:N_upper, size = 1, replace = TRUE)
    upper_limit_rnd <- lower_limit_rnd + round((N_upper - lower_limit_rnd)*runif(0, 1, n = 1))
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    state[[hill_params[m]]] <- bounds
  }
  
  #FC
  for (n in 1:length(FC_params)) {
    lower_limit_rnd <- runif(FC_lower, FC_upper, n = 1)
    upper_limit_rnd <- lower_limit_rnd + (FC_upper - lower_limit_rnd)*runif(0, 1, n = 1)
    bounds[[1]] <- lower_limit_rnd
    bounds[[2]] <- upper_limit_rnd
    state[[FC_params[n]]] <- bounds
  }
  
  state
  
}
##############################update_state########################################################################################################################################################


update_state <- function(current_state, params_changed = 1, update_method = "Random", step_size = 0.1,
                         G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  new_state <- current_state
  
  #number of parameters
  number_of_params <- length(current_state)
  
  #vector of paramater names
  parameter_names <-names(current_state)
  
  #index of each parameter type in vector of parameter names
  production_params <- grep(pattern = "G", parameter_names)
  degradation_params <- grep(pattern = "K", parameter_names)
  threshold_params <- grep(pattern = "TH", parameter_names)
  hill_params <- grep(pattern = "N", parameter_names)
  FC_params <- grep(pattern = "FC", parameter_names)
  
  #index of which params to change
  mutation_indexes <- sample(1:length(new_state), params_changed)
  
  if (update_method == "Random") {
    
    for (l in 1:params_changed) {
      
      #which bound will be changed, lower(0) or upper(1)
      bound <- round(runif(1, 0, 1))
      
      #which direction will be bound be changed, down (0) or up(1)
      direction <- round(runif(1, 0, 1))
      
      #random number between 0-1 to multiply step-size with
      random <- runif(1, 0, 1)

      #if the random index is one of the indexes for production parameters
      if (mutation_indexes[l] %in% production_params) {
        
        #calculate the size of the step
        step <- step_size*(G_upper - G_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - step
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) - step
            
            #if temp lbound is less than global lbound
            if (new_bound < G_lower) {
              
              #new lbound becomes global lbound
              new_state[[mutation_indexes[l]]][1] <- G_lower
              
            } else {
              
              #new lbound becomes temp lbound
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + step 
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) + step
            
            #if temp lbound > global lbound
            if (new_bound > G_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- G_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(current_state[[mutation_indexes[l]]][2])) {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][2])
              
              new_state[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) - step #temp ubound is current ubound - step 
            
            if (new_bound < G_lower) {
              
              new_bound <- G_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(current_state[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][1])
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
        
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) + step #if uound will be increased: temporary ubound is current ubound + step 
            
            if (new_bound > G_upper) {
              
              new_state[[mutation_indexes[l]]][2] <- G_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      
      #if the random index is one of the indexes for degradation parameters
      if (mutation_indexes[l] %in% degradation_params) {
        
        #calculate the size of the step
        step <- step_size*(K_upper - K_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - step
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) - step
            
            #if temp lbound is less than global lbound
            if (new_bound < K_lower) {
              
              #new lbound becomes global lbound
              new_state[[mutation_indexes[l]]][1] <- K_lower
              
            } else {
              
              #new lbound becomes temp lbound
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + step 
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) + step
            
            #if temp lbound > global lbound
            if (new_bound > K_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- K_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(current_state[[mutation_indexes[l]]][2])) {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][2])
              
              new_state[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) - step #temp ubound is current ubound - step 
            
            if (new_bound < K_lower) {
              
              new_bound <- K_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(current_state[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][1])
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) + step #if uound will be increased: temporary ubound is current ubound + step 
            
            if (new_bound > K_upper) {
              
              new_state[[mutation_indexes[l]]][2] <- K_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      
      #if the random index is one of the indexes for threshold parameters
      if (mutation_indexes[l] %in% threshold_params) {
        
        #calculate the size of the step
        step <- step_size*(TH_upper - TH_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - step
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) - step
            
            #if temp lbound is less than global lbound
            if (new_bound < TH_lower) {
              
              #new lbound becomes global lbound
              new_state[[mutation_indexes[l]]][1] <- TH_lower
              
            } else {
              
              #new lbound becomes temp lbound
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + step 
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) + step
            
            #if temp lbound > global lbound
            if (new_bound > TH_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- TH_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(current_state[[mutation_indexes[l]]][2])) {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][2])
              
              new_state[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) - step #temp ubound is current ubound - step 
            
            if (new_bound < TH_lower) {
              
              new_bound <- TH_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(current_state[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][1])
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) + step #if uound will be increased: temporary ubound is current ubound + step 
            
            if (new_bound > TH_upper) {
              
              new_state[[mutation_indexes[l]]][2] <- TH_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
      #if the random index is one of the indexes for hill parameters
      if (mutation_indexes[l] %in% hill_params) {
        
        #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
        if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == N_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != N_lower ) {
          
          #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
          if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
            
            #new lower bound = previous lower - 1
            
            if ((as.numeric(current_state[[mutation_indexes[l]]][1]) - 1) < N_lower) {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1])
              
            } else {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - 1
            }
            
          } else {
            
            #new lower bound = previous lower + 1
            new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + 1
            
          }
          
          
        } else {
          
          #if the two bounds are equal to each upper bound must step up, OR if direction = 1
          if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
            
            #new upper bound = previous upper + 1
            
            if ((as.numeric(current_state[[mutation_indexes[l]]][2]) + 1) > N_upper) {
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2])
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + 1
              
            }
            
          } else {
            
            #new upper bound = previous upper - (previous_upper - 1
            new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - 1
            
          }
          
        }
        
      }
      
      #if the random index is one of the indexes for foldchange parameters
      if (mutation_indexes[l] %in% FC_params) {
        
        #calculate the size of the step
        step <- step_size*(FC_upper - FC_lower)*random
        
        #if lower bound
        if (bound == 0) {
          
          #if lbound will be decreased
          if (direction == 0) {
            
            #temporary lbound is current lbound - step
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) - step
            
            #if temp lbound is less than global lbound
            if (new_bound < FC_lower) {
              
              #new lbound becomes global lbound
              new_state[[mutation_indexes[l]]][1] <- FC_lower
              
            } else {
              
              #new lbound becomes temp lbound
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
            
          } else {
            #if lbound will be increased: temporary lbound is current lbound + step 
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][1]) + step
            
            #if temp lbound > global lbound
            if (new_bound > FC_upper) {
              
              #temp lbound becomes global lbound
              new_bound <- FC_upper
              
            }
            #if temp lbound > current lbound: current ubound becomes new lbound, and temp lbound becomes new ubound 
            if (new_bound > as.numeric(current_state[[mutation_indexes[l]]][2])) {
              
              new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][2])
              
              new_state[[mutation_indexes[l]]][2] <- new_bound
              
              #otherwise, lbound becomes temp lbound
            } else {
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            }
          }
          
          #if upper bound
        } else {
          #if ubound will be decreased
          if (direction == 0) { 
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) - step #temp ubound is current ubound - step 
            
            if (new_bound < FC_lower) {
              
              new_bound <- FC_lower #temp ubound becomes global ubound IF temp ubound < global ubound
              
            }
            
            if (new_bound < as.numeric(current_state[[mutation_indexes[l]]][1])) { #if temp ubound < current lbound: current lbound becomes new ubound, and temp ubound becomes new lbound
              
              new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][1])
              
              new_state[[mutation_indexes[l]]][1] <- new_bound
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp ubound
              
            }
            
          } else {
            
            new_bound <- as.numeric(current_state[[mutation_indexes[l]]][2]) + step #if uound will be increased: temporary ubound is current ubound + step 
            
            if (new_bound > FC_upper) {
              
              new_state[[mutation_indexes[l]]][2] <- FC_upper #new ubound becomes global upper if temp bound > global_upper
              
            } else {
              
              new_state[[mutation_indexes[l]]][2] <- new_bound #otherwise, new ubound becomes temp bound
              
            }
            
          }
          
        }
      }
      
    }    
      
    
  } else if(update_method == "Difference_Between_Bounds") {
  
  #for each time a param needs to be changed
  for (l in 1:params_changed) {
    
    #which bound will be changed, lower(0) or upper(1)
    bound <- round(runif(1, 0, 1))
    
    #which direction will be bound be changed, down (0) or up(1)
    direction <- round(runif(1, 0, 1))
    
    #if the random index is one of the indexes for production parameters
    if (mutation_indexes[l] %in% production_params) {
      
      #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
      if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == G_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != G_lower ) {
        
        #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
        if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
          
          #new lower bound = previous lower - (previous_lower - global_lower)*stepsize, AKA goes down 10% of the current difference
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - (as.numeric(current_state[[mutation_indexes[l]]][1]) - G_lower)*step_size
          
        } else {
          
          #new lower bound = previous lower + (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
        
      } else {
        
        #if the two bounds are equal to each upper bound must step up, OR if direction = 1
        if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
          
          #new upper bound = previous upper + (global upper - previous upper)*stepsize, AKA goes up 10% of the current difference between upper and global
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + (G_upper - as.numeric(current_state[[mutation_indexes[l]]][2]))*step_size
          
        } else {
          
          #new upper bound = previous upper - (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
      }
    }
    
    
    if (mutation_indexes[l] %in% degradation_params) {
      
      #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
      if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == K_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != K_lower ) {
        
        #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
        if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
          
          #new lower bound = previous lower - (previous_lower - global_lower)*stepsize, AKA goes down 10% of the current difference
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - (as.numeric(current_state[[mutation_indexes[l]]][1]) - K_lower)*step_size
          
        } else {
          
          #new lower bound = previous lower + (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
        
      } else {
        
        #if the two bounds are equal to each upper bound must step up, OR if direction = 1
        if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
          
          #new upper bound = previous upper + (global upper - previous upper)*stepsize, AKA goes up 10% of the current difference between upper and global
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + (K_upper - as.numeric(current_state[[mutation_indexes[l]]][2]))*step_size
          
        } else {
          
          #new upper bound = previous upper - (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
      }
    }
    
      
    
    
    
    if (mutation_indexes[l] %in% threshold_params) {
      
      #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
      if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == TH_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != TH_lower ) {
        
        #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
        if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
          
          #new lower bound = previous lower - (previous_lower - global_lower)*stepsize, AKA goes down 10% of the current difference
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - (as.numeric(current_state[[mutation_indexes[l]]][1]) - TH_lower)*step_size
          
        } else {
          
          #new lower bound = previous lower + (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
        
      } else {
        
        #if the two bounds are equal to each upper bound must step up, OR if direction = 1
        if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
          
          #new upper bound = previous upper + (global upper - previous upper)*stepsize, AKA goes up 10% of the current difference between upper and global
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + (TH_upper - as.numeric(current_state[[mutation_indexes[l]]][2]))*step_size
          
        } else {
          
          #new upper bound = previous upper - (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
      }
      
    }
    
    if (mutation_indexes[l] %in% hill_params) {
      
      #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
      if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == N_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != N_lower ) {
        
        #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
        if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
          
          #new lower bound = previous lower - 1
          
          if ((as.numeric(current_state[[mutation_indexes[l]]][1]) - 1) < N_lower) {
            
            new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1])
            
          } else {
          
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - 1
          }
          
        } else {
          
          #new lower bound = previous lower + 1
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + 1
          
        }
        
        
      } else {
        
        #if the two bounds are equal to each upper bound must step up, OR if direction = 1
        if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
          
          #new upper bound = previous upper + 1
          
          if ((as.numeric(current_state[[mutation_indexes[l]]][2]) + 1) > N_upper) {
            
            new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2])
            
          } else {
          
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + 1
          
          }
          
        } else {
          
          #new upper bound = previous upper - (previous_upper - 1
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - 1
          
        }
        
      }
      
    }
    
    if (mutation_indexes[l] %in% FC_params) {
      
      #if (bound = 0 OR upper bound = global ubound) AND (lower bound is NOT equal to global lbound)
      if ((bound == 0 | as.numeric(current_state[[mutation_indexes[l]]][2]) == FC_upper) & as.numeric(current_state[[mutation_indexes[l]]][1]) != FC_lower ) {
        
        #if the two bounds are equal to each other...lower bound must step down, OR if direction = 0
        if (as.numeric(current_state[[mutation_indexes[l]]][1]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 0) {
          
          #new lower bound = previous lower - (previous_lower - global_lower)*stepsize, AKA goes down 10% of the current difference
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) - (as.numeric(current_state[[mutation_indexes[l]]][1]) - FC_lower)*step_size
          
        } else {
          
          #new lower bound = previous lower + (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][1] <- as.numeric(current_state[[mutation_indexes[l]]][1]) + (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
        
      } else {
        
        #if the two bounds are equal to each upper bound must step up, OR if direction = 1
        if (as.numeric(current_state[[mutation_indexes[l]]][2]) == as.numeric(current_state[[mutation_indexes[l]]][2]) | direction == 1) {
          
          #new upper bound = previous upper + (global upper - previous upper)*stepsize, AKA goes up 10% of the current difference between upper and global
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) + (FC_upper - as.numeric(current_state[[mutation_indexes[l]]][2]))*step_size
          
        } else {
          
          #new upper bound = previous upper - (previous_upper - previous_lower)*stepsize, AKA goes up 10% of the difference between two bounds
          new_state[[mutation_indexes[l]]][2] <- as.numeric(current_state[[mutation_indexes[l]]][2]) - (as.numeric(current_state[[mutation_indexes[l]]][2]) - as.numeric(current_state[[mutation_indexes[l]]][1]))*step_size
          
        }
        
      }
    
    
  }
  }
  }
  new_state
  
}

##############################parametric_variation_functions########################################################################################################################################################
parametric_variation_bounds <- function(parameter_bounds, global_parameter) {
  
  bound_min <- (parameter_bounds[[2]] +  parameter_bounds[[1]])/2 - ((parameter_bounds[[2]] -  parameter_bounds[[1]])/2)*global_parameter
  bound_max <- (parameter_bounds[[2]] +  parameter_bounds[[1]])/2 + ((parameter_bounds[[2]] -  parameter_bounds[[1]])/2)*global_parameter
  
  bounds <- list(bound_min = bound_min, bound_max = bound_max)
  bounds
}
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
ref_state0 <- list(list(1,100), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
test_pv <- generate_ref_object(topology = topology_TS, reference_state = ref_state0, numModels = 2000, integrateStepSize = 0.02, global_parameter = 0.8)
test_pv$parametric_variation_reference_state

############################cost_function##########################################################################################################################################################
cost_function <- function(state, reference_rset_object, scoring_method = "PCs") {
  
  if(scoring_method == "Toy") {
    
    cf_vec <- vector(length = length(state))
    
    for (i in 1:length(state)) {
      
      cf_vec[[i]] <- as.numeric(state[[i]][1]) - as.numeric(state[[i]][2])
      
    }
    
    sum(cf_vec)
    
  } else {
    
    simulated_rset <- reference_rset_object$reference_rset
    numModels <- nrow(reference_rset_object$expression_ref)
    
    #hill gene indexes
    hill_genes <- grep(pattern = "N", colnames(sracipeParams(simulated_rset)))
    
    #assign vaues from uniform samples distributions to each model of simulated rset.  
    for (i in 1:length(state)) {
      
      sracipeParams(simulated_rset)[,i] <- runif(numModels, min= as.numeric(state[[i]][1]), max= as.numeric(state[[i]][2]))  
      
    }
    
    #sample and round all the hill coefficients
    for (k in 1:length(hill_genes)) {
      #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
      sracipeParams(simulated_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(state[[hill_genes[k]]][1]) - 0.5, max= as.numeric(state[[hill_genes[k]]][2]) + 0.5))  


    }
    
    # print(state)
    # 
    # plot(density(sracipeParams(simulated_rset)[,1]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,1]))
    # 
    # plot(density(sracipeParams(simulated_rset)[,3]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,3]))
    # 
    # plot(density(sracipeParams(simulated_rset)[,7]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,7]))
    
    #simulate/transform
    sim_expression_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
    
    expr_simulated_log <- log2(t(assay(sim_expression_rset)))
    
    #z normalize
    expr_simulated_log_z <- sweep(expr_simulated_log, 
                                  2, reference_rset_object$means_ref, FUN = "-")
    
    expr_simulated_log_z <- sweep(expr_simulated_log_z, 
                                  2, reference_rset_object$sds_ref, FUN = "/")
    
    #rotate
    simulated_expression_PCs <- expr_simulated_log_z %*% reference_rset_object$prcomp_ref$rotation
    
    #cost_function vector
    cf_vec <- vector(length = ncol(reference_rset_object$expression_ref))
    
    #scoring
    if (scoring_method == "PCs") {
      
      #ks score for distribution of each gene
      for (l in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,l], y = simulated_expression_PCs[,l])
        cf_vec[l] <- as.numeric(ks_test$statistic)
        
      }
    
    } else {
      
      for (l in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,l], y = expr_simulated_log_z[,l])
        cf_vec[l] <- as.numeric(ks_test$statistic)
        
      }
  }

    cf_vec
}
}

#####################Acceptance functions##################################################################################################################################################################
#Acceptance functions
normal_PDF <- function(x, mean, sd) {
  
  y <- 1/(sd*sqrt(2*pi))*exp(-((x-mean)^2)/(2*sd^2))
  y
  
}


Maxwell.Boltzmann_PDF <- function(x, T = 1) {
  
  y <- sqrt(2/pi)*(x^2)*(exp((-x^2)/(2*(T^2)))/(T^3))
  y
  
}

exponential <- function(x, T) {
  
  y <- exp(-x/T)
  y
  
}

T = 0
curve(exponential(x, T = T), from = 0, to = 1, col = "orange"); print(paste(exponential(0.5, T), "(x=0.5)_", exponential(2, T), "(x=2)_", exponential(10, T), "(x=10)"))


T = 1
curve(exponential(x, T = T), from = 0, to = 15, col = "red"); print(paste(exponential(0.5, T), "(x=0.5)_", exponential(2, T), "(x=2)_", exponential(10, T), "(x=10)"))
T = 0.1
curve(exponential(x, T = T), from = 0, to = 15, col = "orange", add = TRUE); print(paste(exponential(0.5, T), "(x=0.5)_", exponential(2, T), "(x=2)_", exponential(10, T), "(x=10)"))
T = 10
curve(exponential(x, T = T), from = 0, to = 15, col = "blue", add = TRUE); print(paste(exponential(0.5, T), "(x=0.5)_", exponential(2, T), "(x=2)_", exponential(10, T), "(x=10)"))
T = 100
curve(exponential(x, T = T), from = 0, to = 15, col = "green", add = TRUE); print(paste(exponential(0.5, T), "(x=0.5)_", exponential(2, T), "(x=2)_", exponential(10, T), "(x=10)"))

T = 1
curve(Maxwell.Boltzmann_PDF(x, T = T), from = 0, to = 5, col = "red")
T = 5
curve(Maxwell.Boltzmann_PDF(x, T = T), from = 0, to = 5, col = "blue", add = TRUE)
T = 10
curve(Maxwell.Boltzmann_PDF(x, T = T), from = 0, to = 5, col = "green", add = TRUE)


##########################MCMC_Metropolis############################################################################################################################################################

MCMC_Metropolis <- function(reference_rset, reference_state, max_iterations = 1000, params_changed = 1, update_method = "Random", step_size = 0.1, scoring_method = "Toy", acceptance_fun = "None", temperature = "NA", annealing_rate = "NA", 
                            min_score = "NA", G_lower = 1, G_upper = 100, K_lower = 0.1, K_upper = 1, TH_lower = 0.02, TH_upper = 1.98, N_lower = 1, N_upper = 5, FC_lower = 1, FC_upper = 100) {
  
  if (scoring_method == "Toy") {
  
  MCMC_object <- list(states_list = vector(mode = "list", length = max_iterations), costs_vec = vector(length = max_iterations) , acceptance = 0)


  MCMC_object$states_list[[1]] <- starting_point(sample_rset = reference_rset, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper,
                                                 TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper) 
  
  MCMC_object$costs_vec[1] <- cost_function(state = MCMC_object$states_list[[1]], scoring_method = scoring_method)

  w <- 2
  
while (w <= max_iterations +1) {

  #cost of current state
  # cost_current <- cost_function(state = MCMC_object$states_list[[w-1]], scoring_method = scoring_method)
  cost_current <- MCMC_object$costs_vec[w-1]
  
  #generate new state
  state_new <- update_state(current_state = MCMC_object$states_list[[w-1]], params_changed = params_changed, update_method = update_method, step_size = step_size,
                            G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
  
  #cost of new state
  cost_new <- cost_function(state = state_new, scoring_method = scoring_method)
  
  #difference in costs
  cost_difference <- cost_new - cost_current
  
  #if cost_new < cost_current, accept new change, else reject
  if (cost_new < cost_current) {
    
    MCMC_object$states_list[[w]] <- state_new
    
    MCMC_object$costs_vec[w] <- cost_new
    
    MCMC_object$acceptance <- MCMC_object$acceptance + 1
    
  } else {
    
    #reject the rejection if random number between 0 and 1 is less than the "acceptance" PDF evaluated at the cost difference between states
      if (acceptance_fun == "Exponential") {
    
        L <- exponential(x = cost_difference, T = temperature)
    
        chance <- runif(1, min = 0, max = 1)
    
        if (chance < L) {
      
          MCMC_object$states_list[[w]] <- state_new
          
          MCMC_object$costs_vec[w] <- cost_new
          
          MCMC_object$acceptance <- MCMC_object$acceptance + 1
  
        } else {
          
          MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
          
          MCMC_object$costs_vec[w] <- cost_current
      
        }
    
  } else if(acceptance_fun == "Boltzmann") {
    
      L <- Maxwell.Boltzmann_PDF(x = cost_difference, T = temperature)
    
      chance <- runif(1, min = 0, max = 1)
    
      if (chance < L) {
      
        MCMC_object$states_list[[w]] <- state_new
      
        MCMC_object$costs_vec[w] <- cost_new
      
        MCMC_object$acceptance <- MCMC_object$acceptance + 1
      
    } else {
      
        MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
      
        MCMC_object$costs_vec[w] <- cost_current
      
    }
    
    
  } else {
    
    MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
    
    MCMC_object$costs_vec[w] <- cost_current
    
    
  }
    
}
#increment iteration umber
  w <- w + 1
  
  
  #break loop if min_score exists and if it is greater than the current score
  if(!missing(min_score) & is.numeric(min_score)) {
    
    if (MCMC_object$costs_vec[w-1] <= min_score) {
     break 
    }
  
  }
  
  #annealing step: if annealing is present and if temperature is greater than 1
  if(!missing(annealing_rate) & is.numeric(temperature) & is.numeric(annealing_rate)) {
    
    if (temperature >= 1) {
    
    temperature <- temperature - annealing_rate
    }
  }

  
}
  MCMC_object$costs_vec <- MCMC_object$costs_vec[1:w-1]
  
  MCMC_object$states_list <- plyr::compact(MCMC_object$states_list)
  
  MCMC_object
  
  } 
  else {
    
    reference_rset_object <- generate_ref_object(reference_rset = reference_rset, reference_state = reference_state)
    
    arguments <- list(reference_rset_object = reference_rset_object, params_changed = params_changed, update_method = update_method, update_step_size = step_size, scoring_method = scoring_method, acceptance_function = acceptance_fun,
                      temperature = temperature, annealing_rate = annealing_rate, min_score = min_score, max_iterations = max_iterations, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper,
                      TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
    MCMC_object <- list(states_list = vector(mode = "list", length = max_iterations), costs_vec = vector(length = max_iterations), costs_list = vector(mode = "list", length = max_iterations), acceptance = 0, arguments = arguments)
    
    MCMC_object$states_list[[1]] <- starting_point(sample_rset = reference_rset, G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper,
                                                   TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
    
    MCMC_object$costs_list[[1]] <- cost_function(state = MCMC_object$states_list[[1]], scoring_method = scoring_method, reference_rset_object = reference_rset_object)
    
    MCMC_object$costs_vec[1] <- sum(MCMC_object$costs_list[[1]])
    
    w <- 2
    
    while (w <= max_iterations +1) {
      
      #cost of current state
      # cost_current <- cost_function(state = MCMC_object$states_list[[w-1]], scoring_method = scoring_method)
      cost_current <- MCMC_object$costs_vec[w-1]
      cost_scores_current <- MCMC_object$costs_list[[w-1]]
      
      #generate new state
      state_new <- update_state(current_state = MCMC_object$states_list[[w-1]], params_changed = params_changed, update_method = update_method, step_size = step_size,
                                G_lower = G_lower, G_upper = G_upper, K_lower = K_lower, K_upper = K_upper, TH_lower = TH_lower, TH_upper = TH_upper, N_lower = N_lower, N_upper = N_upper, FC_lower = FC_lower, FC_upper = FC_upper)
      
      #cost of new state
      cost_scores_new <- cost_function(state = state_new, scoring_method = scoring_method, reference_rset_object = reference_rset_object)
      cost_new <-sum(cost_scores_new)
      
      #difference in costs
      cost_difference <- cost_new - cost_current
      
      #if cost_new < cost_current, accept new change, else reject
      if (cost_new < cost_current) {
        
        MCMC_object$states_list[[w]] <- state_new
        
        MCMC_object$costs_list[[w]] <- cost_scores_new
        
        MCMC_object$costs_vec[w] <- cost_new
        
        MCMC_object$acceptance <- MCMC_object$acceptance + 1
        
      } else {
        
        #reject the rejection if random number between 0 and 1 is less than the "acceptance" PDF evaluated at the cost difference between states
        if (acceptance_fun == "Exponential") {
          
          L <- exponential(x = cost_difference, T = temperature)
          print(L)
          chance <- runif(1, min = 0, max = 1)
          
          if (chance <= L) {
            
            MCMC_object$states_list[[w]] <- state_new
            
            MCMC_object$costs_list[[w]] <- cost_scores_new
            
            MCMC_object$costs_vec[w] <- cost_new
            
            MCMC_object$acceptance <- MCMC_object$acceptance + 1
            
          } else {
            
            MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
            
            MCMC_object$costs_list[[w]] <- cost_scores_current
            
            MCMC_object$costs_vec[w] <- cost_current
            
          }
          
        } else if(acceptance_fun == "Boltzmann") {
          
          L <- Maxwell.Boltzmann_PDF(x = cost_difference, T = temperature)
          
          chance <- runif(1, min = 0, max = 1)
          
          if (chance < L) {
            
            MCMC_object$states_list[[w]] <- state_new
            
            MCMC_object$costs_list[[w]] <- cost_scores_new
            
            MCMC_object$costs_vec[w] <- cost_new
            
            MCMC_object$acceptance <- MCMC_object$acceptance + 1
            
          } else {
            
            MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
            
            MCMC_object$costs_list[[w]] <- cost_scores_current
            
            MCMC_object$costs_vec[w] <- cost_current
            
          }
          
          
        } else {
          
          MCMC_object$states_list[[w]] <- MCMC_object$states_list[[w-1]]
          
          MCMC_object$costs_list[[w]] <- cost_scores_current
          
          MCMC_object$costs_vec[w] <- cost_current
          
          
        }
        
      }
      #increment iteration umber
      w <- w + 1
      
      
      #break loop if min_score exists and if it is greater than the current score
      if(!missing(min_score) & is.numeric(min_score)) {
        
        if (MCMC_object$costs_vec[w-1] <= min_score) {
          break 
        }
        
      }
      
      #annealing step: if annealing is present and if temperature is greater than 1
      if(!missing(annealing_rate) & is.numeric(temperature) & is.numeric(annealing_rate)) {
        
        if (temperature > 0) {
          
          temperature <- temperature - annealing_rate
        }
        
        if (temperature <= 0) {
          
          temperature <- 0
        }
        
      }
      
      
    }
    MCMC_object$costs_vec <- MCMC_object$costs_vec[1:w-1]
    
    MCMC_object$costs_list <- plyr::compact(MCMC_object$costs_list)
    
    MCMC_object$states_list <- plyr::compact(MCMC_object$states_list)
    
    MCMC_object

  }
  
}

##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset_10000_0.02 <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE, integrateStepSize = 0.02)

ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")

ref_object_10000_0.02 <- generate_ref_object(reference_rset = reference_rset, reference_state = ref_state1)


MCMC_toy_test <- MCMC_Metropolis(reference_rset = ref_object_10000_0.02$reference_rset, reference_state = ref_object_10000_0.02$reference_state, scoring_method = "Toy", max_iterations = 2000, min_score = -604, params_changed = 8,
                              update_method = "Random", step_size = 0.25, acceptance_fun = "None", temperature = "NA", annealing_rate = "NA", TH_lower = 0, TH_upper = 100)

MCMC_toy_test$costs_vec
##########################################################################################################################################################################################################################################
MCMC_test1 <- MCMC_Metropolis(reference_rset = reference_rset, scoring_method = "PCs", max_iterations = 2000, min_score = 0.04, params_changed = 1,
                              update_method = "Random", step_size = 0.25, acceptance_fun = "None", temperature = "NA", annealing_rate = "NA", TH_lower = 0, TH_upper = 100)

plot(MCMC_test1$costs_vec)
tail(MCMC_test1$states_list, n = 1)
tail(MCMC_test1$costs_vec)
##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE)
sracipeParams(reference_rset)[,1] <- runif(min = 1, max = 10, n = 10000)

MCMC_test2 <- MCMC_Metropolis(reference_rset = reference_rset, scoring_method = "PCs", max_iterations = 2000, min_score = 0.04, params_changed = 3,
                              update_method = "Random", step_size = 0.5, acceptance_fun = "None", temperature = "NA", annealing_rate = "NA", TH_lower = 0, TH_upper = 100)

plot(MCMC_test2$costs_vec)
tail(MCMC_test2$states_list, n= 1)
tail(MCMC_test2$costs_vec)

##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE)

ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")

ref_object_10000_0.02 <- generate_ref_object(reference_rset = reference_rset, reference_state = ref_state1)


MCMC_test3 <- MCMC_Metropolis(reference_rset = ref_object_10000_0.02$reference_rset, reference_state = ref_object_10000_0.02$reference_state, scoring_method = "PCs", max_iterations = 2000, min_score = 0.025, params_changed = 8,
                              update_method = "Random", step_size = 0.75, acceptance_fun = "None", temperature = "NA", annealing_rate = "NA", TH_lower = 0, TH_upper = 100)

plot(MCMC_test3$costs_vec)
tail(MCMC_test3$states_list, n = 1)
MCMC_test3$costs_vec[[2000]]
ref_object <- generate_ref_object(reference_rset = reference_rset)

names(ref_state)
ref_state <- sim_state_final
ref_state[[1]] <- list(1,10)
ref_state[[2]] <- list(1,100)
ref_state[[3]] <- list(0.1,1)
ref_state[[4]] <- list(0.1,1)
ref_state[[5]] <- list(0,55)
ref_state[[6]] <- list(0,55)
ref_state[[7]] <- list(1,5)
ref_state[[8]] <- list(1,5)
ref_state[[9]] <- list(1,100)
ref_state[[10]] <- list(1,100)
ref_state_test <- cost_function(state = ref_state, scoring_method = "PCs", reference_rset_object = ref_object)
sum(ref_state_test)

sim_state_final <- MCMC_test2$states_list[[2000]]
sim_state1 <- sim_state_final
sim_state1[[9]] <- list(1, 100)
sim_state_test1 <- cost_function(state = sim_state1, scoring_method = "PCs", reference_rset_object = ref_object)

sim_state2 <- sim_state_final
sim_state2[[8]] <- list(1, 5)
sim_state_test2 <- cost_function(state = sim_state2, scoring_method = "PCs", reference_rset_object = ref_object)
sum(sim_state_test2)

##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE, integrateStepSize = 0.02)

ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")

MCMC_test4 <- MCMC_Metropolis(reference_rset = reference_rset, reference_state = ref_state1, scoring_method = "PCs", max_iterations = 2000, min_score = 0.025, params_changed = 4,
                              update_method = "Random", step_size = 0.6, acceptance_fun = "Exponential", temperature = 1, annealing_rate = 0.00049, TH_lower = 0, TH_upper = 100)

plot(MCMC_test4$costs_vec)
which(MCMC_test4$costs_vec == min(MCMC_test4$costs_vec))
MCMC_test2$states_list[[1961]]
MCMC_test4$acceptance/2000

##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE, integrateStepSize = 0.02)

ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")

MCMC_test5 <- MCMC_Metropolis(reference_rset = reference_rset, reference_state = ref_state1, scoring_method = "PCs", max_iterations = 2000, min_score = 0.025, params_changed = 5,
                              update_method = "Random", step_size = 0.95, acceptance_fun = "Exponential", temperature = 1, annealing_rate = 0.001, TH_lower = 0, TH_upper = 100)

plot(MCMC_test5$costs_vec)
which(MCMC_test5$costs_vec == min(MCMC_test5$costs_vec))
MCMC_test5$states_list[[2000]]

##########################################################################################################################################################################################################################################
cost_function_modified <- function(state, reference_rset_object, scoring_method = "PCs") {
  
  if(scoring_method == "Toy") {
    
    cf_vec <- vector(length = length(state))
    
    for (i in 1:length(state)) {
      
      cf_vec[[i]] <- as.numeric(state[[i]][1]) - as.numeric(state[[i]][2])
      
    }
    
    sum(cf_vec)
    
  } else {
    
    simulated_rset <- reference_rset_object$reference_rset
    numModels <- nrow(reference_rset_object$expression_ref)
    
    #hill gene indexes
    hill_genes <- grep(pattern = "N", colnames(sracipeParams(simulated_rset)))
    
    #assign vaues from uniform samples distributions to each model of simulated rset.  
    for (i in 1:length(state)) {
      
      sracipeParams(simulated_rset)[,i] <- runif(numModels, min= as.numeric(state[[i]][1]), max= as.numeric(state[[i]][2]))  
      
    }
    
    #sample and round all the hill coefficients
    for (k in 1:length(hill_genes)) {
      #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
      sracipeParams(simulated_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(state[[hill_genes[k]]][1]) - 0.5, max= as.numeric(state[[hill_genes[k]]][2]) + 0.5))  
      
      
    }
    
    # print(state)
    # 
    # plot(density(sracipeParams(simulated_rset)[,1]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,1]))
    # 
    # plot(density(sracipeParams(simulated_rset)[,3]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,3]))
    # 
    # plot(density(sracipeParams(simulated_rset)[,7]))
    # plot(density(sracipeParams(reference_rset_object$reference_rset)[,7]))
    
    #simulate/transform
    sim_expression_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
    
    expr_simulated_log <- log2(t(assay(sim_expression_rset)))
    
    #z normalize
    expr_simulated_log_z <- sweep(expr_simulated_log, 
                                  2, reference_rset_object$means_ref, FUN = "-")
    
    expr_simulated_log_z <- sweep(expr_simulated_log_z, 
                                  2, reference_rset_object$sds_ref, FUN = "/")
    
    #rotate
    simulated_expression_PCs <- expr_simulated_log_z %*% reference_rset_object$prcomp_ref$rotation
    
    #cost_function vector
    cf_vec <- vector(length = ncol(reference_rset_object$expression_ref))
    
    #cost_function_results_list
    cf_results <- list(cf_vec = cf_vec, simulated_rset = sim_expression_rset, simulated_expression = 1)
    
    #scoring
    if (scoring_method == "PCs") {
      
      cf_results$simulated_expression <- simulated_expression_PCs
      
      #ks score for distribution of each gene
      for (l in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$prcomp_ref$x[,l], y = simulated_expression_PCs[,l])
        cf_vec[l] <- as.numeric(ks_test$statistic)
        
      }
      
    } else {
      
      cf_results$simulated_expression <- expr_simulated_log_z
      
      for (l in 1:length(cf_vec)) {
        
        ks_test <- dgof::ks.test(x = reference_rset_object$expression_ref[,l], y = expr_simulated_log_z[,l])
        cf_vec[l] <- as.numeric(ks_test$statistic)
        
      }
    }
    cf_results$cf_vec <- cf_vec
    cf_results
  }
}



cf_results_test <- cost_function_modified(state = MCMC_test2$states_list[[1961]], reference_rset_object = MCMC_test5$arguments$reference_rset_object, scoring_method = "PCs")
cf_results_test$simulated_expression

MCMC_test2$states_list[[1961]]

plot(density(MCMC_test5$arguments$reference_rset_object$expression_ref[,1]))
lines(density(cf_results_test$simulated_expression[,1]), col = "blue")

##########################################################################################################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE, integrateStepSize = 0.02)

ref_state1 <- list(list(1,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
names(ref_state1) <- c("G_A", "G_B", "K_A", "K_B", "TH_B_A", "TH_A_B", "N_B_A", "N_A_B", "FC_B_A", "FC_A_B")

MCMC_test6 <- MCMC_Metropolis(reference_rset = reference_rset, reference_state = ref_state1, scoring_method = "PCs", max_iterations = 2000, min_score = 0.025, params_changed = 1,
                              update_method = "Random", step_size = 0.4, acceptance_fun = "Exponential", temperature = 0.03, annealing_rate = 0.00003, TH_lower = 0, TH_upper = 100)

ref_state2 <- list(list(5,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state3 <- list(list(5,15), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state4 <- list(list(10,15), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state5 <- list(list(3,12), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state6 <- list(list(2,11), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state7 <- list(list(1,11), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))
ref_state8 <- list(list(2,10), list(1,100), list(0.1,1), list(0.1, 1), list(0,55), list(0,55), list(1,5), list(1,5), list(1,100), list(1,100))

test_ref_state8 <- cost_function(state = ref_state8, reference_rset_object = ref_rset_object, scoring_method = "Default")
sum(test_ref_state8)
sum(test_ref_state7)
sum(test_ref_state6)
sum(test_ref_state5)
sum(test_ref_state4)
sum(test_ref_state3)
sum(test_ref_state2)
sum(test_ref_state1)

ref_rset_object <- generate_ref_object(reference_rset = reference_rset, reference_state = ref_state1)

plot(MCMC_test6$costs_vec)
which(MCMC_test6$costs_vec == min(MCMC_test6$costs_vec))
MCMC_test6$states_list[[2000]]

