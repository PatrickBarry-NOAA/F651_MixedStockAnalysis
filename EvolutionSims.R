#### Allele Freq Change Simulations ####

library(tidyverse)

source("./F651_SimFunctions.R")

set.seed(123) # Always set a seed if you want to be able to replicate your results

#### Genetic Drift ####
# Start with an allele frequency of 0.5
# Population size of 100 (so 200 diploid alleles)
# And 100 Generations

simulate_drift(p = 0.5, N = 100, G = 100, replicates = 20)

#What happens as you decrease the starting allele frequency?

#What happens as you change the population size?


#### Mutation with Drift ####
#u controls the mutation from p to q

#First lets ensure that allele freq are stable when the pop size is 
#super high and there is no mutation
simulate_drift_two_way_mutation(p = 0.5, N = 100000, G = 100, u = 0, v = 0, replicates = 10)

#Now lets introduce some mutation
simulate_drift_two_way_mutation(p = 0.5, N = 1000, G = 200, u = 0.05, v = 0.01, replicates = 10)

#What happens when we introduce a reasonable mutation rate 10^-8?

#### Migration #### 
# We'll simulate a simple island model with a population that receives migrants
# So now we'll have pop sizes for the source (N_source), sink (N_sink),
# allele freq for the source (p_source) and sink (p_sink),
# migration rate from the source to sink


simulate_source_sink(p_source = 0.9, p_sink = 0.1, 
                     N_source = 10000, N_sink = 500, 
                     m = 0.1, G = 50, replicates = 10)

# How many migrants are going from the source and into the sink here? 

# What if only ~1% of the sink pop in any given year was migrants

