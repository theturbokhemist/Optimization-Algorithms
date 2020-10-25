#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:22:33 2020

@author: gordid
"""
import numpy as np
from numpy.random import randint
from random import random as rnd
from random import gauss, randrange
import copy


class GA:
     "This class is for recording all the data produced from the genetic algorithm"
     def __init__(self):
        self.generations = []
        self.fitness_values = []
        self.selected = []
        self.number_of_genes = None
        self.number_of_individuals = None

GA_ranges = GA()
vars(GA)    

def individual(number_of_genes, lower_limit, upper_limit):
    parameters_matrix = [0]*number_of_genes
    for i in range(number_of_genes):
        lower_limit_rnd = round(rnd()*upper_limit)
        upper_limit_rnd = lower_limit_rnd + round(rnd()*(upper_limit-lower_limit_rnd))
        parameters_matrix[i] = [lower_limit_rnd, upper_limit_rnd]
    GA_ranges.number_of_genes = number_of_genes
    return(parameters_matrix)

 
a = (individual(5,1,100))
GA_ranges.__dict__
print(a)


def population(number_of_individuals,
               number_of_genes, lower_limit, upper_limit):
    population_list = [0]*number_of_individuals
    for i in range(number_of_individuals):
        population_list[i] = individual(number_of_genes, lower_limit, upper_limit)
    GA_ranges.number_of_genes = number_of_genes
    GA_ranges.number_of_individuals = number_of_individuals
    return(population_list)

b = (population(8,2,1,100))
print(b)
d = b
GA_ranges.__dict__

def fitness_value(individual):
    ranges_list = [0]*len(individual)
    for i in range(len(individual)):
        ranges_list[i] = individual[i][1] - individual[i][0]
    fv = sum(ranges_list)
    return(fv)

print(fitness_value(a))

        
def selection(generation, method = 'Fittest Half'):
    generations_list = [0]*len(generation)
    fitness_values = [0]*len(generation) 
    
    fv_list = [0]*len(generation) 
    for i in range(len(generation)):
        fv_list[i] = fitness_value(generation[i])
        
    normalized_fv = [0]*len(fv_list)
    for j in range(len(fv_list)):
        normalized_fv[j] = fv_list[j]/sum(fv_list)
        
    if method == 'Fittest Half':
        if len(generation) % 4 == 0:
            fv_index_sorted = sorted(range(len(fv_list)), key = lambda k: fv_list[k], reverse = True)
            for l in range(len(fv_index_sorted)):
                fitness_values[l] = fv_list[fv_index_sorted[l]]
                generations_list[l] = generation[fv_index_sorted[l]]
                
            GA_ranges.fitness_values.append(fitness_values)
            GA_ranges.generations.append(generations_list)
            selected = generations_list[:int(len(generation)/2)]
            GA_ranges.selected.append(selected)
            return(selected)                   
    elif method == 'Roulette Wheel':
        print(fv_list)
    elif method == 'Random':
        print(fv_list)              
        
c = selection(b)
print(c)    
selection(b)
GA_ranges.__dict__

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
        if method == "Reset":
            for i in range(mutation_rate):
                bound = round(rnd())
                if bound == 0:
                    mutated_individual[i][0] = round(rnd()*individual[i][1])
                else:
                    mutated_individual[i][1] = upper_limit - round(rnd()*(upper_limit-individual[i][0]))      
            return(mutated_individual)
        elif method == "Gauss":
            return(mutated_individual)

            
f = mutation(d[0], 100, 1, 1)
f
    

def reproduction(paired, mating_method = 'Single Point', upper_limit = 100, lower_limit = 1):
    unmutated = []
    next_generation = copy.deepcopy(paired)
    paired_off_list = list(zip(paired[0::2], paired[1::2]))
    for i in range(len(paired_off_list)):
        unmutated.extend(mating(paired_off_list[i]))
    for j in range(len(unmutated)):
        next_generation.append(mutation(unmutated[j], upper_limit = upper_limit, lower_limit = lower_limit))
    return(next_generation)
 
c       
d = reproduction(c)        
d        
d[0]

def termination()

class Error(Exception):
   """Base class for other exceptions"""
   pass

class MutationRateTooHigh(Error):
   """Mutation Rate cannot be greater than the number of genes"""
   pass


rnd()
round(rnd())
int(rnd())
int(0.9)
z = zip(c[0::2], c[1:2])
list(z)
      

testlist = [2,3,1,4,10,9,7, 8]
testlist
testlist[0:2]
testlist.append(1)
testlist = testlist[0,2,5,4,1]
print(len(testlist))
print(range(len(testlist)))
print(sorted(range(len(testlist))))
sorted(range(len(testlist)), key=lambda k: testlist[k])
int(4.0)

print(testlist[:int(len(testlist)/2)])
print(GA_ranges.fitness_values[0])
print(GA_ranges.fitness_values[0].sort())
int(5/2)
test3 = [1,2] + [3,4]
print(test3)

data = c
data = [[[1,2],[3,4]],[[5,6], [7,8]],[[10,20],[30,40]],[[50,60], [70,80]]]
data[0]
q = zip(data[0::2], data[1::2])
r = list(q)
len(r)
r[0]
for i,k in zip(data[0::2], data[1::2]):
    print(i + k)
    
a = zip([1,2,1], ["a", "b", "a"], [True, False, True])
next(a)
print(a)
list(a)
set(a)
d = set([1,2,3,3])
print(d)
list(d)
