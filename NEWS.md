### 0.1.2 

* Fixed an error: Inside the function predict.ptsr, the argument xreg was not being passed to the internal function ".ptsr.predict.". 

* Fixed an error: in the internal function ".ptsr.predict.", the MA lags were incorrectly specified. This was causing an error in the prediction for some models.


### 0.1.1 

* Added check for infinity and NAN in the log-likelihood function. Now, "Inf" and "NAN" are replaced by 2.2e+10 and "-Inf" is replaced by -2.2e+10 avoiding issues with the optim function.

* Added an internal function to calculate the numerical hessian. This function is only used if the inversion of the output of numDeriv::hessian fails.

* Release date: 	2022-01-13

### 0.1.0

* First commit.
* Release date: 	2022-01-13