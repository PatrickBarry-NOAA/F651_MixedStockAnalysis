#### Chinook Salmon Baseline #####

# Go ahead and download the chinook baseline from
# https://datadryad.org/dataset/doi:10.5061/dryad.8063
# You should have ChinoookSNPGeno.csv and a readme
# Or better yet, just install MsGSI cause Bobby includeds the Templin baselin already formatted 
# Note if you want the full baseline use the dryad link, Bobby only supplies 28 loci

#install Ms.GSI and rubias - you can use devtools:: or remotes:: to install from github
devtools::install_github("boppingshoe/Ms.GSI", build_vignettes = F)
devtools::install_github("eriqande/rubias")

library(Ms.GSI)

#Print out the reporting group summary
templin_pops211 %>%
  group_by(repunit,grpvec) %>%
  summarize(nPops = n()) %>%
  arrange(grpvec)

#Now lets look at the baseline itself
base_templin

# Let's see how well we can assign fish from their reference population if 
# we remove them from the baseline

sa_chinook <- self_assign(reference = base_templin, gen_start_col = 5)

sa_to_repu <- sa_chinook %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

sa_to_repu %>%
  group_by(indiv,repunit) %>%
  slice_max(repu_scaled_like, with_ties = FALSE) %>% #pick the highest prob of assignment to repuint
  ungroup() %>%
  group_by(repunit) %>%
  summarize(nSamps = n(),
            nCorrect = sum(repunit == inferred_repunit),
            propCorrect = mean(repunit == inferred_repunit) 
  )

# What does this mean about assigning individuals from these groups?
# How did the lower yukon do?

#Lets simulate a mixture and then try to estimate the proportion
# we'll use the assess_reference_loo function in rubias

chin_sims <- assess_reference_loo(reference = base_templin, 
                                  gen_start_col = 5, 
                                  alpha_repunit = .05,
                                  reps = 50, 
                                  mixsize = 200)

chin_sims_sum <- chin_sims %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

ggplot(chin_sims_sum, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.1, slope = 1,lty=2,col='red') +
  geom_abline(intercept = -0.1, slope = 1,lty=2,col='red') +
  facet_wrap(~ repunit)+
  xlab("True Proportion")+
  ylab("Estimated Mean Proportion")

# We like to see 95% of the simulations within 10% of the true value (the red dashed lines)

# What if we put the lower yukon river in the Coastal Western Alaska Group?

base_templin_mod <- base_templin %>%
  mutate(repunit = ifelse(repunit=='Lower Yukon','Coastal West Alaska',repunit))

#Now lets redo those simulations and see what the effect is on the Coastal Western Alaska Group

# <Insert your code here>

