#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:00:46 2020

@author: gordid
"""

import numpy as np
from numpy.random import randint
from random import random as rnd
from random import gauss, randrange
import copy
import random

class GA:
     "This class is for recording all the data produced from the genetic algorithm"
     def __init__(self):
        self.generations = []
        self.fitness_values = []
        self.selected = []
        self.number_of_genes = None
        self.population_size = None
        
class Error(Exception):
   """Base class for other exceptions"""
   pass

class MutationRateTooHigh(Error):
   """Mutation Rate cannot be greater than the number of genes"""
   pass

def individual(number_of_genes, lower_limit, upper_limit):
    individual = [0]*number_of_genes
    for i in range(number_of_genes):
        lower_limit_rnd = round(rnd()*upper_limit)
        upper_limit_rnd = lower_limit_rnd + round(rnd()*(upper_limit-lower_limit_rnd))
        individual[i] = [lower_limit_rnd, upper_limit_rnd]
    return(individual)

def population(population_size,
               number_of_genes, lower_limit, upper_limit):
    population_list = [0]*population_size
    for i in range(population_size):
        population_list[i] = individual(number_of_genes, lower_limit, upper_limit)
    return(population_list)

def fitness_value(individual):
    ranges_list = [0]*len(individual)
    for i in range(len(individual)):
        ranges_list[i] = individual[i][1] - individual[i][0]
    fv = sum(ranges_list)
    return(fv)
        
def selection(generation, method = 'Fittest Half'):
    generations_list = [0]*len(generation)
    fitness_values = [0]*len(generation) 
    
    fv_list = [0]*len(generation) 
    for i in range(len(generation)):
        fv_list[i] = fitness_value(generation[i])
        
    # normalized_fv = [0]*len(fv_list)
    # for j in range(len(fv_list)):
    #     normalized_fv[j] = fv_list[j]/sum(fv_list)
        
    if method == 'Fittest Half':
        if len(generation) % 4 == 0:
            fv_index_sorted = sorted(range(len(fv_list)), key = lambda k: fv_list[k], reverse = True)
            for l in range(len(fv_index_sorted)):
                fitness_values[l] = fv_list[fv_index_sorted[l]]
                generations_list[l] = generation[fv_index_sorted[l]]
                
            
            selected = generations_list[:int(len(generation)/2)]
            return(selected)
                   
    elif method == 'Roulette Wheel':
        print(fv_list)
    elif method == 'Random':
        print(fv_list)              
        
def pairing(selected, method = 'Fittest'):
    if method == 'Fittest':
        if len(selected) % 4 == 0:
            paired = selected
        return(paired)                   
    elif method == 'Random':
        return(paired)
    elif method == 'Weighted Random':
        return(paired)

def mating(mating_pair, method = 'Single Point', hotspot = False):
    if method == 'Single Point':
        if hotspot == False:
            pivot_point = randint(1, len(mating_pair[0]))
            offspring = [mating_pair[0][0:pivot_point]+mating_pair[1][pivot_point:]]
            offspring.append(mating_pair[1][0:pivot_point]+mating_pair[0][pivot_point:])
            return(offspring)
        else:
            return(offspring)
    elif method == 'Double Point':
       return(offspring)
    elif method == 'Random':
       return(offspring)
   
 
def mutation(individual, upper_limit, lower_limit, mutation_rate = 1, method = 'Reset', sd = 0.001):
    if mutation_rate > len(individual):
        print("Mutation Rate cannot be greater than the number of genes")
        raise MutationRateTooHigh
    else:
        mutated_individual = copy.deepcopy(individual)
        mutation_indexes = random.sample(range(len(individual)), mutation_rate)
        
        if method == "Reset":
            for i in range(mutation_rate):
                lower = round(rnd()*upper_limit)
                upper = lower + round(rnd()*(upper_limit-lower))                
                mutated_individual[mutation_indexes[i]][0] = lower
                mutated_individual[mutation_indexes[i]][1] = upper     
            return(mutated_individual)
        
        elif method == "Reset Bound":
            for i in range(mutation_rate):
                bound = round(rnd())
                if bound == 0:
                    mutated_individual[mutation_indexes[i]][0] = round(rnd()*individual[mutation_indexes[i]][1])
                else:
                    mutated_individual[mutation_indexes[i]][1] = upper_limit - round(rnd()*(upper_limit-individual[mutation_indexes[i]][0]))      
            return(mutated_individual)
        elif method == "Gauss":
            return(mutated_individual)


def reproduction(paired, mating_method = 'Single Point', hotspot = False,
                 mutation_rate = 1, mutation_method = 'Reset', sd = 0.001, upper_limit = 100, lower_limit = 1):
    unmutated = []
    next_generation = copy.deepcopy(paired)
    paired_off_list = list(zip(paired[0::2], paired[1::2]))
    for i in range(len(paired_off_list)):
        unmutated.extend(mating(paired_off_list[i], method = mating_method, hotspot = hotspot))
    for j in range(len(unmutated)):
        next_generation.append(mutation(unmutated[j], upper_limit = upper_limit, lower_limit = lower_limit, mutation_rate = mutation_rate,
                                        sd = sd, method = mutation_method))
    return(next_generation)

  

def GA_RACIPE(GA_object, population_size = 8, number_of_genes = 4, number_of_generations = 10, mutation_rate = 1,
              lower_limit = 1, upper_limit = 100, selection_method = 'Fittest Half', pairing_method = 'Fittest',
              mating_method = 'Single Point', hotspot = False, mutation_method = 'Reset', sd = 0.001, termination = True, termination_percent = 0.9):
    GA_object.number_of_genes = number_of_genes
    GA_object.population_size = population_size
    
    generation = population(population_size, number_of_genes, lower_limit, upper_limit)
    fvs = [0]*len(generation) 
    for j in range(len(generation)):
        fvs[j] = fitness_value(generation[j])
    
    if termination == False:
    
        for i in range(number_of_generations):
            GA_object.generations.append(generation)
            
            fvs = [0]*len(generation) 
            for j in range(len(generation)):
                fvs[j] = fitness_value(generation[j])
            GA_object.fitness_values.append(fvs)
        
            selected = selection(generation, method = selection_method)    
            GA_object.selected.append(selected)
        
            paired = pairing(selected, method = pairing_method)
        
            generation = reproduction(paired, mating_method = mating_method, hotspot = hotspot,
                                       mutation_rate = mutation_rate, mutation_method = mutation_method,
                                       sd = sd, upper_limit = upper_limit, lower_limit = lower_limit)
            print(np.mean(fvs))
    else:
        while fvs[0] < termination_percent*upper_limit*number_of_genes:
            print(fvs[0])
            GA_object.generations.append(generation)
            
            fvs = [0]*len(generation) 
            for j in range(len(generation)):
                fvs[j] = fitness_value(generation[j])
            GA_object.fitness_values.append(fvs)
        
            selected = selection(generation, method = selection_method)    
            GA_object.selected.append(selected)
        
            paired = pairing(selected, method = pairing_method)
        
            generation = reproduction(paired, mating_method = mating_method, hotspot = hotspot,
                                       mutation_rate = mutation_rate, mutation_method = mutation_method,
                                       sd = sd, upper_limit = upper_limit, lower_limit = lower_limit)
        


random.sample(range(1,5),4)
    
GA_test = GA()
GA_RACIPE(GA_test, number_of_genes = 10, number_of_generations = 20000, mutation_rate = 1, termination = True, lower_limit = 1, upper_limit = 100,
          termination_percent = .8, mutation_method = 'Reset Bound')
print(GA_test.fitness_values[-1])
len(GA_test.fitness_values)
GA_test.generations[-1]

GA_test.__dict__

data = [2,4]
np.mean(data)

b = (population(8,2,1,100))
print(b)
c = selection(b)
print(c)
d = reproduction(c)
print(d)   
1+1
