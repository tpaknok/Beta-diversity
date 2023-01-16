# Environmental availablity and uniqueness 

This repository contains all codes required to reproduce results from the paper Tsang et al. (XXXXX)

This function is likely to evolve based on users' feedback, and the latest version can be found in a separate file.

site_dist.R is the function required to compute Uniche. A compositional or distance matrix can be provided to the function. Formula can be specificed as in GAM / GAMM function in mgcv, but always start with pair_dist as the LHS of the formula.

For random effect and correlation structure specifications, please refer to the documentation of gamm. 

Further arguments (e.g. optimizer, maximum number of iterations) will also be passed on to gam/gamm. 
