# AdaPTGMM


## Overview
This package implements AdaPTGMM: [link to be added](https://github.com/patrickrchao/AdaPTGMM/). AdaPTGMM is a flexible multiple testing method that uses arbitrary covariates to model the local false discovery rate.


The main method `adapt_gmm` allows flexibility in covariates, testing type, input, and classification model. We include default classification model implementations of using neural networks, multinomial logistic regression, glmnet, Generalized Additive Models (GAM).

Users may use their own custom beta models, read the vignette for more details.


## Installation         

```
# install.packages("devtools")
devtools::install_github("patrickrchao/AdaPTGMM")
```


### Simulation Example
We create simulations similar to Section 4.4 of the paper.

```
# Load package
library("AdaPTGMM")
library("splines")

# Generate data
n <- 3000

x <- rnorm(n)
expit <- function(x) 1/(1+exp(-x))
pi1 <-   expit(2*(x-1))
H <- 1 * (runif(n) < pi1)
theta <- H * rlogis(n, location=2, scale=1)
z <- rnorm(n,mean=theta)
# Two sided testing
pvals <- 2 * pnorm(abs(z), lower.tail = FALSE)



# Run adapt_gmm
x <- data.frame(x = x)
formulas = paste("splines::ns(x,df=", c(2,3,4), ")")
alphas <- c(0.01,0.05,0.1,0.2)

# Can use p-values
res <- adapt_gmm(x=x,pvals=pvals, alphas=alphas,
        beta_formulas = formulas, model_type = "neural", nclasses= c(3,4))
        
# Can use z-values
zres <- adapt_gmm(x=x,z=z, alphas=alphas,
        beta_formulas = formulas, model_type = "neural", nclasses= c(3,4))
        
      
# Fitted Gaussian means
print(paste("Means:", zres$params$mu, collapse=" "))

# Fitted Gaussian variances
print(paste("Variances:", zres$params$var, collapse=" "))

# Selected beta formula from model selection
print(paste("Selected formula:", format(zres$args$beta_formula), collapse=" "))

# Selected number of classes from model selection
print(paste("Selected number of classes:", zres$args$nclasses, collapse=" "))
```
