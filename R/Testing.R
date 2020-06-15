# n = 200
# z = c(rnorm(n*0.8,mean=0),rnorm(n*0.1,mean=1.5,sd=1.3),rnorm(n*0.1,mean=3,sd = 1.5))
# p_values=pnorm(z,lower.tail = FALSE)
# output <- adapt_gmm(x=runif(n),z=z,testing = "interval", rendpoint=2,selection = "BIC")
