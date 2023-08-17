# Environmental availablity and uniqueness
# (Latest version - 2023-08-17)

This repository contains all codes required to reproduce results in Tsang et al. (2024)

Cont.R - Results based on a simulated continuous variable - as described in Supplementary Text S1
empirical_tree.R - Empirical analyses using BCI data - as described in the main text
simulation.R - Results based on a simulated categorical variable - as described in the main text 
site_dist.R - Function to correct ecological uniqueness under uneven environmental availability

The site_dist.R function is likely to evolve based on users' feedback. Both the latest and old version will be kept here

# Installing site_dist.R

The easiest way is to use devtools and the function source_url. Alternatively, you can simply download the script and load it in R using source.

```
library(devtools)
source_url("https://raw.githubusercontent.com/tpaknok/Beta-diversity/main/site_dist.R")
```

# Description 

site_dist.R calculates compositional uniqueness within a niche space using Generalized Additive Model or Generalized Additive Mixed Model

empircal_tree.R and simulation.R are the scripts to reproduce results described in Tsang et al. (2023, Methods in Ecology & Evolution)
```
site_dist(formula,random=NULL,correlation=NULL,spatial_var = NULL,
          mat,env_data,dissim="jaccard",C=0.5,
          dissim_mat = NULL,
          env_space = NULL,
          family="gaussian",
          weights=NULL,
          length.cont=25,...)
```

# Arguments

|Argument|Description|
|---|---|
|formula| A GAM/GAMM formula. Note that the LHS must start with pair_dist (see Example). Syntax as in mgcv.|
|random| Optional; A random effect structure, as in mgcv. Default is NULL.|
|correlation| Optional; A correlation structure, as in mgcv. Default is NULL.|
|spatial_var| Optional; Name of the spatial variable, specificed as a vector (e.g. c("Long","Lat"). Default is NULL.|
|mat| A community data (species as column, site as row) as a data frame. Ignored if dissimilarity matrix is provided in dissim_mat.|
|env_data| A data frame containing environmental data. Must be in data frame format.|
|dissim| Default is jaccard dissimilarity, but other indices can be used. Values will be passed on to vegdist in vegan. See ?vegdist for a list of possible inputs.|
|C| For developmental purposes. Used for rarefaction pairwise beta diversity to obtain pairwise beta at specific completeness.|
|dissim_mat| Optional; A dissimilarity matrix as a dist object or data frame.|
|env_space| Optional; A data frame with environmental conditions specified in the model. If not provided the function will create a niche space automatically based on minima and maxima of each variable, and the number of virtual site for each axis equates to length.cont.|
|family| Default is gaussian. Other distribution can be used. See gamm in mgcv.|
|weights| Weight of each observation. Default is NULL - all observations have the same weight.|
|length.cont| Number of evenly-distributed sites on each environmental gradient. Default is 25.|
|...| Other arguments that will be passed on to gam/gamm in mgcv|

# Notes

If factors are ordered, the function will convert them into numeric variables before the regression. The order will follow how the factors were ordered originally.
For both ordinal and nominal variables, in the formula simply set pair_dist~Factor1+Factor2+Factor3....etc. Note that no splines can be added for all factors. 

You can also add continuous variable (e.g. pair_dist~Factor1 + s(Cont.var, k=3)).

The users are recommended to read the documentation of gam and gamm in mgcv.

This function will generate N evenly-distributed virtual sites along each environmental gradient, and exhaust all possible combinations of them. Therefore, if 25 evenly-distributed sites were generated for each gradient, the total number of sites would be 25^N

For total beta diversity, users can choose to provide the compositional data or its distance matrix. For other beta diversity, such as turnover component of beta diversity, functional trait and phylogenetic beta diversity, users should supply a distance matrix

# Author
Toby P.N. Tsang (paknok.tsang@utoronto.ca)

# Example

```
library(BiodiversityR)
library(vegan)

data(BCI) #load BCI data
data(BCI.env) #load topographic variable

#may take some time to run!
Start <- Sys.time()
predict_pb <- site_dist(pair_dist ~ s(elevation,k=3) + s(convex,k=3) + s(slope,k=3),
                        spatial_var = c("UTM.EW","UTM.NS"),
                        correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),
                        mat=BCI,env_data=BCI.env,
                        length.cont=25,
                        control=lmeControl(opt = 'optim',maxIter = 1e8, msMaxIter = 1e8)) #increased the number of iterations, and specified optimizer

End <- Sys.time()
End-Start

m_niche <- gls(avg_dis~elevation+convex+slope,
               correlation=corExp(1,form=~UTM.EW+UTM.NS,nugget=T),
               data=predict_pb)
summary(m_niche)
```
