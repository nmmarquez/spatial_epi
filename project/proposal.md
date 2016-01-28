# BYM2, Blending Error, and Covariate Bias

### Background  
In a recent study by Riebler et al[1], an analysis was conducted on a new
parameterization of a recently proposed model, the BYM2 model. The authors 
first state the importance of using scaling in a model in order for the 
interpretation and setting of hyper priors to be intuitive and then propose 
a new framework for model selection, the penalized complexity(PC) prior. 

The BYM2 model has two hyper priors which a modeler needs to select for, 
$\phi$ and $\tau_b$. When $\phi$ is $0$ the model is one that accounts for 
pure overdispersion of the relative risk attributed to random error. When 
$\phi$ is $1$ the model random effects are correlated via an intrinsic 
Gaussian Markov random field (referred to as the Besag model). The PC 
prior is one in which simplicity is the assumption, models need not 
introduce complexity if there is no proof of it in the data. This prior 
causes the BYM2 model to reduce itself in complexity from the Besag model to 
the overdispersion model to constant relative risk. This model is then tested 
against simulated data where the expected values are fluctuated along side 
different scenarios of overdispersion, i.e. constant relative risk, spatially 
uncorrelated and correlated error.  

### Proposal  
There are at least two remaining scenarios where the BYM2 model with the PC 
prior should be tested.

1) The effect of blended sources of variation ie spatially and non spatially correlated. 
2) The effects that the PC prior has on potential bias in $\beta$ estimation.

These two points still need to be tested in the BYM2 model against more 
traditional models and the same methods of evaluation used by Riebler et al[1]  
could be used to test these effects. It will be important to see the 
effects of a PC prior when simulated data is a mixture of spatially correlated 
and uncorrelated overdispersion. 


### References  
[1] Andrea Riebler, Sigrunn H. Sorbye, Daniel Simpson, Harvard Rue. An 
intuitive Bayesian spatial model for disease mapping that accounts for scaling 
(2016).