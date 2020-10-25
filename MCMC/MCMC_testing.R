###################
topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))
reference_rset <- sracipeSimulate(topology_TS, integrate = FALSE, numModels = 10000, genParams = TRUE)

hist(sracipeParams(reference_rset)[,10])
range(sracipeParams(reference_rset)[,10])

sracipeParams(reference_rset)[,1] <- runif(min = 5, max = 30, n = 10000)
sracipeParams(reference_rset)[,2] <- runif(min = 70, max = 95, n = 10000)

start <- starting_point(reference_rset)
cost_function(start, scoring_method = "Test")
update1 <- update_state(current_state = start, params_changed = 10, update_method = "Random", step_size = 0.2)

for (q in 1:length(start)) {
  
  print(start[[q]])
  print(update1[[q]])
  print("separate")
  
}
######################################################################################################################################################################################
##TESTING with Toy Model##
test_length <- 100
max_iterations <- 5000
params_changed <- 1
update_method <- "Random"
step_size <- 0.1
min_score <- -405
acceptance_fun <- "Exponential"
temperature <- 1
annealing_rate <- "NA"

arguments_MCMC <- list(test_length = test_length,
                       max_iterations = max_iterations,
                       params_changed = params_changed,
                       update_method = update_method,
                       step_size = step_size,
                       min_score = min_score,
                       acceptance_fun = acceptance_fun,
                       temperature = temperature,
                       annealing_rate = annealing_rate)


list_of_MCMC_objects <- vector(mode = "list", length = test_length)

name <- paste0("MasterMCMCObject__", "UpdateMethod.", update_method, "__", "StepSize.", step_size, "__", "ParamsChanged.", params_changed, "__",
               "AcceptanceFun.", acceptance_fun, "__", "Temperature.", temperature, "__", "AnnealingRate.", annealing_rate)

assign(name, list(MCMC_objects_list = list_of_MCMC_objects, arguments = arguments_MCMC))

for (k in 1:test_length) {
  
  .GlobalEnv[[name]]$MCMC_objects_list[[k]] <- MCMC_Metropolis(topology = topology_TS, reference = "None", numModels = 1000, scoring_method = "Test", max_iterations = max_iterations, min_score = min_score, params_changed = params_changed,
                                                               update_method = update_method, step_size = step_size, acceptance_fun = acceptance_fun, temperature = temperature, annealing_rate = annealing_rate)
  
}



######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA


i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
update_method.random_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  update_method.random_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  update_method.random_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}
######################################################################################################################################################################################
#MasterMCMCObject__UpdateMethod.Difference_Between_Bounds__StepSize.0.1__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Difference_Between_Bounds__StepSize.0.1__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA


i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
update_method.difference_between_bounds_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  update_method.difference_between_bounds_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  update_method.difference_between_bounds_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.2__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.2__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA


i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.2_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.2_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.2_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.3__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.3__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.3_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.3_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.3_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.4__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.4__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.4_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.4_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.4_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.8__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.8__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.8_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.8_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.8_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.10__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.10__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.10_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.10_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.10_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.7__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.7__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.7_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.7_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.7_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.9__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.1__ParamsChanged.9__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
ParamsChanged.9_list <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  ParamsChanged.9_list$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  ParamsChanged.9_list$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.25__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.25__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
StepSize.0.25 <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  StepSize.0.25$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  StepSize.0.25$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.5__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.5__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA

i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
StepSize.0.5 <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  StepSize.0.5$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  StepSize.0.5$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.25__ParamsChanged.8__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.25__ParamsChanged.8__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
StepSize.0.25__ParamsChanged.8 <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  StepSize.0.25__ParamsChanged.8$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  StepSize.0.25__ParamsChanged.8$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################

#MasterMCMCObject__UpdateMethod.Random__StepSize.0.75__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
temp <- MasterMCMCObject__UpdateMethod.Random__StepSize.0.75__ParamsChanged.1__AcceptanceFun.Random_Walk__Temperature.NA__AnnealingRate.NA
i <- 1
vec_lengths <- vector(length = length(temp$MCMC_objects_list))
vec_accept <- vector(length = length(temp$MCMC_objects_list))


#CHANGE THIS
StepSize.0.75 <- list(lengths = vec_lengths, accept = vec_accept)

for (i in 1:length(temp$MCMC_objects_list)) {
  
  StepSize.0.75$lengths[i] <- length(temp$MCMC_objects_list[[i]]$costs_vec)
  StepSize.0.75$accept[i] <- temp$MCMC_objects_list[[i]]$acceptance/(length(temp$MCMC_objects_list[[i]]$costs_vec))
}

######################################################################################################################################################################################


lengths_mean <- c(mean(update_method.random_list$lengths), mean(update_method.difference_between_bounds_list$lengths), mean(ParamsChanged.2_list$lengths), mean(ParamsChanged.3_list$lengths), mean(ParamsChanged.4_list$lengths),
                  mean(ParamsChanged.7_list$lengths), mean(ParamsChanged.8_list$lengths), mean(ParamsChanged.9_list$lengths), mean(ParamsChanged.10_list$lengths), mean(StepSize.0.25$lengths), mean(StepSize.0.5$lengths), mean(StepSize.0.75$lengths), mean(StepSize.0.25__ParamsChanged.8$lengths))

lengths_sd <- c(sd(update_method.random_list$lengths), sd(update_method.difference_between_bounds_list$lengths), sd(ParamsChanged.2_list$lengths), sd(ParamsChanged.3_list$lengths), sd(ParamsChanged.4_list$lengths),
                sd(ParamsChanged.7_list$lengths), sd(ParamsChanged.8_list$lengths), sd(ParamsChanged.9_list$lengths), sd(ParamsChanged.10_list$lengths), sd(StepSize.0.25$lengths), sd(StepSize.0.5$lengths), sd(StepSize.0.75$lengths), sd(StepSize.0.25__ParamsChanged.8$lengths))

accept_mean <- c(mean(update_method.random_list$accept), mean(update_method.difference_between_bounds_list$accept), mean(ParamsChanged.2_list$accept), mean(ParamsChanged.3_list$accept), mean(ParamsChanged.4_list$accept),
                 mean(ParamsChanged.7_list$accept), mean(ParamsChanged.8_list$accept), mean(ParamsChanged.9_list$accept), mean(ParamsChanged.10_list$accept), mean(StepSize.0.25$accept), mean(StepSize.0.5$accept), mean(StepSize.0.75$accept), mean(StepSize.0.25__ParamsChanged.8$accept))

accept_sd <- c(sd(update_method.random_list$accept), sd(update_method.difference_between_bounds_list$accept), sd(ParamsChanged.2_list$accept), sd(ParamsChanged.3_list$accept), sd(ParamsChanged.4_list$accept),
               sd(ParamsChanged.7_list$accept), sd(ParamsChanged.8_list$accept), sd(ParamsChanged.9_list$accept), sd(ParamsChanged.10_list$accept), sd(StepSize.0.25$accept), sd(StepSize.0.5$accept), sd(StepSize.0.75$accept), sd(StepSize.0.25__ParamsChanged.8$accept))

names <- as.factor(c("update_method.random_list", "update_method.difference_between_bounds_list", "ParamsChanged.2_list", "ParamsChanged.3_list", "ParamsChanged.4_list", "ParamsChanged.7_list",  "ParamsChanged.8_list",  "ParamsChanged.9_list", "ParamsChanged.10_list", "StepSize.0.25",
                     "StepSize.0.5", "StepSize.0.75", "StepSize.0.25__ParamsChanged.8" ))


MCMC_tests_df <- data.frame(lengths_mean = lengths_mean, lengths_sd = lengths_sd,
                            accept_mean = accept_mean, accept_sd = accept_sd, names = names)


MCMC_tests_lengths_plot <-ggplot(data=MCMC_tests_df, aes(x=reorder(names, lengths_mean), y=lengths_mean, fill = names)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=lengths_mean-lengths_sd, ymax=lengths_mean+lengths_sd), width=.2,
                position=position_dodge(.9))
MCMC_tests_lengths_plot

MCMC_tests_accept_plot <-ggplot(data=MCMC_tests_df, aes(x=reorder(names, accept_mean), y=accept_mean, fill = names)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=accept_mean-accept_sd, ymax=accept_mean+accept_sd), width=.2,
                position=position_dodge(.9))
MCMC_tests_accept_plot

#change T to 0.05
#change annealing rate to 0.001
#with exponential as the acceptance function

