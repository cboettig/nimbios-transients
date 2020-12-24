Simulation study for Bayesian hierarchical model.

R package longTransients must be installed:

  install.packages("path/to/longTransients_1.1.tar.gz", repos = NULL)

simulate_trajectories.R creates the simulated data stored in /data

The simulation_study.R script fits both models to several subsets of the 10 simulated trajectories for a single value of the measurement error variance. The script also creates Figure 7.

single_combinations.R fits both models to one subset of the simulated trajectories.

kl_combinations.R fits the parametric model to subsets of the simulated trajectories used to create Figure 5.

kl_divergence.R summarizes the results of the fits in kl_combinations.R through numerical approximations of marginal KL divergences to create Figure 5.

neg_grad_potential_plot.R creates Figure 6.