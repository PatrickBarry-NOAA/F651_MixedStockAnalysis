#### Sim Functions ####

simulate_drift <- function(p = 0.5, N = 100, G = 100, replicates = 10) {
  
  # Total number of alleles in a diploid population
  total_alleles <- 2 * N 
  
  # Create an empty data frame to store results
  results <- data.frame(
    Generation = integer(),
    Frequency = numeric(),
    Replicate = factor()
  )
  
  # Run the simulation for the specified number of replicates
  for (rep in 1:replicates) {
    
    # Initialize a vector to hold the allele frequencies for this replicate
    freqs <- numeric(G)
    freqs[1] <- p
    
    # Loop through each generation
    for (gen in 2:G) {
      # Draw the number of alleles for the next generation
      # from a binomial distribution
      alleles_drawn <- rbinom(n = 1, size = total_alleles, prob = freqs[gen - 1])
      
      # Calculate the new allele frequency
      freqs[gen] <- alleles_drawn / total_alleles
    }
    
    # Append the replicate's data to the main results dataframe
    rep_data <- data.frame(
      Generation = 1:G,
      Frequency = freqs,
      Replicate = as.factor(rep)
    )
    
    results <- rbind(results, rep_data)
  }
  
  # Plot the results
  p_plot <- ggplot(results, aes(x = Generation, y = Frequency, color = Replicate)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = paste("Genetic Drift Simulation (N =", N, ")"),
      subtitle = paste("Starting allele frequency (p) =", p, "| Generations (G) =", G),
      x = "Generation",
      y = "Allele Frequency"
    )
  
  print(p_plot)
  return(invisible(results))
}


# Genetic Drift and mutation
simulate_drift_mutation <- function(p = 0.5, N = 100, G = 100, u = 0.001, replicates = 10) {
  
  # Total number of alleles in a diploid population
  total_alleles <- 2 * N 
  
  # Create an empty data frame to store results
  results <- data.frame(
    Generation = integer(),
    Frequency = numeric(),
    Replicate = factor()
  )
  
  # Run the simulation for the specified number of replicates
  for (rep in 1:replicates) {
    
    freqs <- numeric(G)
    freqs[1] <- p
    
    for (gen in 2:G) {
      current_p <- freqs[gen - 1]
      
      # 1. MUTATION: p mutates to q at rate u
      # The probability of remaining p is (1 - u)
      mutated_p <- current_p * (1 - u)
      
      # 2. DRIFT: Draw alleles based on the mutated probability
      alleles_drawn <- rbinom(n = 1, size = total_alleles, prob = mutated_p)
      
      # Calculate the new allele frequency
      freqs[gen] <- alleles_drawn / total_alleles
    }
    
    rep_data <- data.frame(
      Generation = 1:G,
      Frequency = freqs,
      Replicate = as.factor(rep)
    )
    
    results <- rbind(results, rep_data)
  }
  
  # Plot the results
  p_plot <- ggplot(results, aes(x = Generation, y = Frequency, color = Replicate)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = paste("Drift & Mutation Simulation (N =", N, ")"),
      subtitle = paste("Start p =", p, "| Gens =", G, "| Mutation (p to q) u =", u),
      x = "Generation",
      y = "Frequency of Allele 'p'"
    )
  
  print(p_plot)
  return(invisible(results))
}

simulate_drift_two_way_mutation <- function(p = 0.5, N = 100, G = 100, u = 0.001, v = 0.001, replicates = 10) {
  
  total_alleles <- 2 * N 
  
  results <- data.frame(
    Generation = integer(),
    Frequency = numeric(),
    Replicate = factor()
  )
  
  for (rep in 1:replicates) {
    
    freqs <- numeric(G)
    freqs[1] <- p
    
    for (gen in 2:G) {
      current_p <- freqs[gen - 1]
      
      # MUTATION: 
      # p alleles that don't mutate: current_p * (1 - u)
      # q alleles that mutate to p: (1 - current_p) * v
      expected_p <- (current_p * (1 - u)) + ((1 - current_p) * v)
      
      # DRIFT: Draw alleles based on the expected probability after mutation
      alleles_drawn <- rbinom(n = 1, size = total_alleles, prob = expected_p)
      
      # Calculate the new allele frequency
      freqs[gen] <- alleles_drawn / total_alleles
    }
    
    rep_data <- data.frame(
      Generation = 1:G,
      Frequency = freqs,
      Replicate = as.factor(rep)
    )
    
    results <- rbind(results, rep_data)
  }
  
  # Plot the results
  p_plot <- ggplot(results, aes(x = Generation, y = Frequency, color = Replicate)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = paste("Drift & Two-Way Mutation (N =", N, ")"),
      subtitle = sprintf("Start p = %g | Gens = %d | p->q (u) = %g | q->p (v) = %g", p, G, u, v),
      x = "Generation",
      y = "Frequency of Allele 'p'"
    )
  
  print(p_plot)
  return(invisible(results))
}

library(ggplot2)
library(dplyr)
library(tidyr)

# Function to simulate genetic drift with source-to-sink migration
simulate_source_sink <- function(p_source = 0.8, p_sink = 0.1, 
                                 N_source = 1000, N_sink = 50, 
                                 m = 0.05, G = 100, replicates = 5) {
  
  total_alleles_source <- 2 * N_source
  total_alleles_sink <- 2 * N_sink
  
  results <- data.frame(
    Generation = integer(),
    Source_Freq = numeric(),
    Sink_Freq = numeric(),
    Replicate = factor()
  )
  
  for (rep in 1:replicates) {
    
    freq_source <- numeric(G)
    freq_sink <- numeric(G)
    
    freq_source[1] <- p_source
    freq_sink[1] <- p_sink
    
    for (gen in 2:G) {
      curr_source <- freq_source[gen - 1]
      curr_sink <- freq_sink[gen - 1]
      
      # MIGRATION: 
      # Source frequency is unaffected by the sink (unidirectional migration)
      exp_source <- curr_source
      
      # Sink frequency is a mix of its native alleles (1-m) and source migrants (m)
      exp_sink <- (curr_sink * (1 - m)) + (curr_source * m)
      
      # DRIFT: Draw alleles based on expected probabilities
      alleles_drawn_source <- rbinom(n = 1, size = total_alleles_source, prob = exp_source)
      alleles_drawn_sink <- rbinom(n = 1, size = total_alleles_sink, prob = exp_sink)
      
      freq_source[gen] <- alleles_drawn_source / total_alleles_source
      freq_sink[gen] <- alleles_drawn_sink / total_alleles_sink
    }
    
    rep_data <- data.frame(
      Generation = 1:G,
      Source_Freq = freq_source,
      Sink_Freq = freq_sink,
      Replicate = as.factor(rep)
    )
    
    results <- rbind(results, rep_data)
  }
  
  # Reshape data for plotting both populations easily
  plot_data <- results %>%
    pivot_longer(cols = c(Source_Freq, Sink_Freq), 
                 names_to = "Population", 
                 values_to = "Frequency")
  
  # Plot the results
  p_plot <- ggplot(plot_data, aes(x = Generation, y = Frequency, 
                                  color = Population, group = interaction(Replicate, Population))) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    scale_color_manual(values = c("Source_Freq" = "#1f77b4", "Sink_Freq" = "#d62728"),
                       labels = c("Sink", "Source")) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    labs(
      title = "Source-Sink Migration Model",
      subtitle = sprintf("Source: N=%d, p0=%g | Sink: N=%d, p0=%g | Migration (m)=%g", 
                         N_source, p_source, N_sink, p_sink, m),
      x = "Generation",
      y = "Allele Frequency"
    )
  
  print(p_plot)
  return(invisible(results))
}