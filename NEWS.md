### 0.1.1 

* Added check for infinity and NAN in the log-likelihood function. Now, "Inf" and "NAN" are replaced by 2.2e+10 and "-Inf" is replaced by -2.2e+10 avoiding issues with the optim function.

* Added an internal function to calculate the numerical hessian. This function is only used if the inversion of the output of numDeriv::hessian fails.


### 0.1.0

* First commit.
* Release date: 	2022-01-13