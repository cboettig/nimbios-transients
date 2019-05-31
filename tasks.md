


- Long-running replicate vs ensemble dynamics (Jody)
- Perturbing the system multiple times

- Likelihood / posterior distributions from a single trajectory vs from the ensemble (Carl) 

- How divergent can the posterior distribution of a parameter from a given model be when estimated from two different realizations of the same model? (Chris)

- Conditions for overconfident learning? (posterior is much narrower than it should be) (Henry)

- Quantify learning from the data (KL divergence prior -> posterior) (Kai)

- First passage time problem for ghost attractor? (Carl -> Sebastian)
- Hierarchical estimation of a over replicates (replicates draw their a's from the distribution rather than fixed)  (Lizzie)

- Black-box / ML (Jorge, Jiang)
- Mixture model, change-point analysis, linear regression

----

Narrative for Reimer et al (2020)

Motivation: 
- change points
- typical analysis (multiple or self citations): 
  - change-point analysis, 
  - EWS prediction
  - mixture model
  - ...

- model misspecification can be bad
- even having the right model might not help if the data come from a transient
- estimation on the right model can even be misleading or overconfident (?)

- illustrate consequences of poor inference

So what can we do? 
- replicate data (or perturbation experiments) may be essential
- mcmc is hard for some problems...
- manage for multiple outcomes
  
  







