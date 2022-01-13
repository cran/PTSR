#  internal function
#  extracts mu(t) and eps(t)
.extract <- function(par, a, xreg, xregar, b, yt, p, q, arlag, malag,
                     ddist, g1, g1.inv, g2) {

  n <- length(yt)
  mut <- numeric(n)
  epst  <- numeric(n + q)

  # adding initial values for yt.
  # Using the unconditional mean
  if (p > 0) yt <- c(rep(mean(yt[1:p]), p), yt)
  g2y <- g2(yt)

  alpha <- 0
  beta <- 0
  l <- 0       # counter

  # alpha
  if (a == 1) {alpha <- par[1]; l <- l + 1}
  # beta
  if (b > 0) {beta <- par[(l + 1):(l + b)]; l <- l + b}
  # ar
  if (p > 0) {
    if (is.null(arlag)) {
      phi <- par[(l + 1):(l + p)]
      l <- l + p
    }
    else {
      phi <- numeric(p)
      phi[arlag] <- par[(l + 1):(l + length(arlag))]
      l <- l + length(arlag)
    }
  }
  # ma
  if (q > 0) {
    if (is.null(arlag)) {
      theta <- par[(l + 1):(l + q)]
      l <- l + q
    }
    else {
      theta <- numeric(q)
      theta[malag] <- par[(l + 1):(l + length(malag))]
      l <- l + length(malag)
    }
  }

  # X*beta
  XB <- numeric(n)
  if (!is.null(xreg)) XB <- xreg %*% beta
  if(p > 0){
    # setting the initial p values as the unconditional mean
    if(xregar)	XB = c(rep(mean(XB[1:p]),p), XB)
  }

  # starting the loop to calculate mu(t) and eps(t)
  for (t in 1:n) {
    ls <- alpha + XB[t+p]

    # shift on t to compensate for "p" and "q"
    if (p > 0){
      xr = 0
      if(xregar) xr = XB[(p + t) - c(1:p)]
      ls <- ls + sum(phi * (g2y[(p + t) - c(1:p)] - xr))
    }
    if (q > 0) ls <- ls + sum(theta * epst[(q + t) - (1:q)])

    mut[t] <- g1.inv(ls)
    epst[q + t] <- yt[p + t] - mut[t]
  }

  return(list(mut = mut, epst = epst[(q+1):(q+n)]))

}

# internal function
# calculates the sum of the -log-likelihood
# for a given distribution
.ll <- function(x, mu, varphi, ddist){
  ll <- sum(ddist(x, mu, varphi, log = TRUE))
  return(-ll)
}


#  internal function
#  sum of the -log-likelihood for the selected model
.loglike <- function(par, a, xreg, xregar, b, yt, p, q, arlag, malag,
                     ddist, g1, g1.inv, g2){
  temp <- .extract(par = par, a = a, xreg = xreg, xregar = xregar, b = b,
                   yt = yt, p = p, q = q, arlag = arlag, malag = malag,
                   ddist = ddist, g1 = g1, g1.inv = g1.inv, g2 = g2)
  # if ddist does not have a varphi parameter, this value will be
  # ignored by the function so there is no problem on setting varphi = par[length(par)]
  sll <- .ll(yt, mu = temp$mut, varphi = par[length(par)], ddist = ddist)
  return(sll)
}


#' Title Function to fit a PTSR model
#'
#' @param start a vector with the starting values for the non-fixed coefficients
#'  of the model.
#' @param yt the time series
#' @param xreg optionally, a vector or matrix of external regressors. Default is \code{NULL}
#' @param xregar logical, if \code{FALSE}, the regressors are not included in the
#' AR component of the model. Default is \code{TRUE}.
#' @param fit.alpha logical, if FALSE, alpha is set to zero. Default is \code{TRUE}
#' @param p order of the AR polinomial
#' @param q order of the MA polinomial
#' @param arlag the lags to be included in the AR polinomial. Default is \code{NULL}, meaning that all lags will be included.
#' @param malag the lags to be included in the MA polinomial. Default is \code{NULL}, meaning that all lags will be included.
#' @param ddist function, the density function to be used
#' @param link1 character indicating which link must be used for \eqn{\mu_t}.  See \code{\link{ptsr.link}} for available links. Default is \sQuote{log}.
#' @param link2 character indicating which link must be used for \eqn{y_t} in the AR recursion.  See \code{\link{ptsr.link}} for available links. Default is \sQuote{identity}.
#' @param g1 optionally, a link function to be used  for \eqn{\mu_t}. Default is \code{NULL}, so that it is calculated internally, using \code{link1}.
#' @param g1.inv optionally, a the inverse link function to be used  for \eqn{\eta_t}. It must the the ivnerse of \code{g1}. Default is \code{NULL}, so that it is calculated internally, using \code{link1}.
#' @param g2 optionally, a link function to be used  for \eqn{y_t}. Default is \code{NULL}, so that it is calculated internally, using \code{link2}.
#' @param method  The method to be used. See [optim][stats::optim] for details.
#' @param ... Further arguments to be passed to \code{optim}.
#'
#' @return
#' The same arguments return by \code{optim} plus a the  following arguments
#' \itemize{
#'    \item \code{coefficients}: a vector with the estimated coefficients;
#'    \item \code{sll}: the sum of the log-likelihood for the fitted model;
#'    \item \code{series}: the original time series;
#'    \item \code{xreg}: the regressors (if any);
#'    \item \code{fitted.values}:  the conditional mean, which corresponds to
#' the in-sample forecast, also denoted fitted values;
#'    \item \code{residuals}: the observed minus the fitted values;
#'    \item \code{model}: a list with the configurations used to fit the model.
#' }
#'
#' @examples
#'
#' #-------------------------------------------------------------------
#' # Gamma-ARMA(1,1) model with no regressors
#' #-------------------------------------------------------------------
#'
#' simu = ptsr.sim(n = 3000, burn = 50,
#'                 varphi = 20, alpha = 0,
#'                 phi = 0.35, theta = 0.2,
#'                 seed = 1234, rdist = r.gamma,
#'                 link1 = "log", link2 = "log")
#'
#' fit1 = ptsr.fit(start =  c(0,0,0,10), yt = simu$yt,
#'                fit.alpha = TRUE, p = 1, q = 1,
#'                ddist = d.gamma, link1 = "log",
#'                link2 = "log", method = "L-BFGS-B")
#' summary(fit1)
#'
#' # removing alpha from the model
#' fit2 = ptsr.fit(start =  c(0,0,10), yt = simu$yt,
#'                fit.alpha = FALSE, p = 1, q = 1,
#'                ddist = d.gamma, link1 = "log",
#'                link2 = "log", method = "L-BFGS-B")
#' summary(fit2)
#'
#' @export
ptsr.fit <- function(start, yt,  xreg = NULL, xregar = TRUE,
                     fit.alpha = TRUE, p = 0, q = 0,
                     arlag = NULL, malag = NULL,
                     ddist = d.gamma, link1 = "log", link2 = "identity",
                     g1 = NULL, g1.inv = NULL, g2 = NULL,
                     method = "L-BFGS-B",...){

  # par = c(alpha, beta, phi, theta, varphi) if ddist has 2 parameters
  # par = c(alpha, beta, phi, theta) if ddist only has mu

  # function to be applied to mu(t) and eta(t)
  if(is.null(g1) | is.null(g1.inv)){
    lk1 <- ptsr.link(link = link1)
    g1 <- lk1$linkfun
    g1.inv <- lk1$linkinv
  }

  # function to be applied to y(t)
  if(is.null(g2)) {
		lk2 <- ptsr.link(link = link2)
		g2 <- lk2$linkfun
		}

  b <- 0
  if(!is.null(xreg)){
    xreg <- as.matrix(xreg)
    b <- ncol(xreg)
  }

  temp <- stats::optim(par = start, fn = .loglike,
                       a = as.integer(fit.alpha),
                       xreg = xreg, xregar = xregar, b = b, yt = yt,
                       p = p, q = q, arlag = arlag, malag = malag,
                       ddist = ddist, g1 = g1, g1.inv = g1.inv, g2 = g2,
                       method = method,...)

  final <- .extract(par = temp$par, a = as.integer(fit.alpha),
                    xreg = xreg, xregar = xregar, b = b, yt = yt,
                    p = p, q = q, arlag = arlag, malag = malag,
                    ddist = ddist, g1 = g1, g1.inv = g1.inv, g2 = g2)

  out <- c()
  out$optim <- temp
  out$model <- list(a = as.integer(fit.alpha),
                    b = b, xregar = xregar, p = p, q = q,
                    arlag = arlag, malag = malag,
                    g1 = g1, g1.inv = g1.inv, g2 = g2,
                    ddist = ddist)

  out$coefficients <- temp$par
  out$sll <- -temp$value
  out$series <- yt
  out$xreg <- NULL
  if(b > 0) out$xreg <- xreg
  out$fitted.values <- final$mut
  out$residuals <- final$epst

  class(out) <- c("ptsr", class(out))
  invisible(out)

}



#' @title Summary Method of class PTSR
#'
#' @description \code{summary} method for class \code{"ptsr"}.
#'
#' @name summary
#'
#' @aliases summary.ptsr
#' @aliases print.summary.ptsr
#'
#' @param object object of class \code{"ptsr"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'
#' @return
#' The function \code{summary.ptsr} computes and returns a list
#' of summary statistics of the fitted model given in \code{object}.
#' Returns a list of class \code{summary.ptsr}, which contains the
#' following components:
#'
#'  \item{residuals}{the residuals of the model.}
#'
#'  \item{coefficients}{a \eqn{k \times 4}{k x 4} matrix with columns for
#'  the estimated coefficient, its standard error, z-statistic and corresponding
#'  (two-sided) p-value.}
#'
#'  \item{sigma.res}{the square root of the estimated variance of the random
#'  error \deqn{\hat\sigma^2 = \frac{1}{n-k}\sum_i{r_i^2},}{\sigma^2 = \frac{1}{n-k} \sum_i r[i]^2,}
#'  where \eqn{r_i}{r[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
#'
#'  \item{df}{degrees of freedom, a 3-vector \eqn{(k, n-k, k*)}, the first
#'  being the number of non-aliased coefficients, the last being the total
#'  number of coefficients.}
#'
#'  \item{vcov}{a \eqn{k \times k}{k \times k} matrix of (unscaled) covariances.
#'  The inverse ov the information matrix.}
#'
#'  \item{loglik}{the sum of the log-likelihood values}
#'
#'  \item{aic}{the AIC value. \eqn{AIC = -2*loglik+2*k}.}
#'
#'  \item{bic}{the BIC value. \eqn{BIC = -2*loglik + log(n)*k}.}
#'
#'  \item{hqc}{the HQC value. \eqn{HQC = -2*loglik + log(log(n))*k}.}
#'
#' @importFrom stats pnorm
#'
#' @export
#' @md
summary.ptsr <- function(object,...) {

  if(!"ptsr" %in% class(object))
    stop("the argument 'object' must be a 'ptsr' object")

  temp <- object$model

  # calculates the numerical hessian
  hs <- numDeriv::hessian(func = .loglike, x = object$coefficients,
                          yt = object$series, a = temp$a,
                          xreg = object$xreg, xregar = temp$xregar,
                          b = temp$b, p = temp$p, q = temp$q,
                          arlag = temp$arlag, malag = temp$malag,
                          ddist = temp$ddist, g1 = temp$g1,
                          g1.inv = temp$g1.inv, g2 = temp$g2)
  sv <- try(solve(hs))
  if ("try-error" %in% class(sv)) {
    print("Error - Returning the hessian")
    return(hs)
  }

  npar <- length(object$coefficients)
  ans <- c()
  ans$residuals <- object$residuals
  n <- length(object$residuals)
  rdf <- ans$df.residuals <- n-npar
  ans$sigma.res <- sqrt(sum(ans$residuals^2)/rdf)
  class(ans) <- "summary.ptsr"

  ans$df = c(npar, rdf)
  ans$vcov <- sv
  stderror <- sqrt(diag(abs(ans$vcov)))
  zstat <- abs(object$coefficients/stderror)
  ans$coefficients <- cbind(Estimate = object$coefficients,
                            `Std. Error` = stderror,
                            `z value` = zstat,
                            `Pr(>|t|)` = 2*(1 - pnorm(zstat)))

  nms = NULL
  if (temp$a > 0 ) nms = c("alpha")
  if (temp$b > 0 ) nms = c(nms, paste("beta(",1:temp$b,")", sep = ""))
  ptemp <- 0
  if (temp$p > 0) {
    lags = 1:temp$p
    ptemp <- temp$p
    if (!is.null(temp$arlag)){
      lags <- temp$arlag
      ptemp <- length(lags)
    }
    nms <- c(nms, paste("ar(",lags,")", sep = ""))
  }
  qtemp <- 0
  if (temp$q > 0) {
    lags <- 1:temp$q
    qtemp <- temp$q
    if (!is.null(temp$malag)){
      lags <- temp$malag
      qtemp <- length(lags)
    }
    nms <- c(nms, paste("ma(",lags,")", sep = ""))
  }
  np <- temp$a + temp$b + ptemp + qtemp
  nphi = 0
  if (length(object$coefficients) > np){
    nms <- c(nms, "varphi")
    nphi = 1
  }
  np <- temp$a + temp$b + temp$p + temp$q + nphi
  ans$df = c(npar, rdf, np)
  rownames(ans$coefficients) <- nms

  ans$loglik <- object$sll
  ans$aic <- -2*ans$loglik+2*npar
  ans$bic <- -2*ans$loglik + log(n)*npar
  ans$hq <- -2*ans$loglik + log(log(n))*npar
  ans$coefficients
  return(ans)

}


#' Users are not encouraged to call these internal functions directly.
#' Internal functions for package PTSR.
#'
#' @title Print Method of class PTSR
#'
#' @description
#' Print method for objects of class \code{ptsr}.
#'
#' @param x object of class \code{ptsr}.
#' @param digits  minimal number of significant digits, see
#' \code{\link{print.default}}.
#' @param ... further arguments to be passed to or from other methods.
#' They are ignored in this function
#'
#' @return Invisibly returns its argument, \code{x}.
#'
#' @importFrom stats coef
#'
#' @export
#'
print.ptsr <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{

  if(length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}


#-----------------------------------------------
# Internal function for printing the summary
#-----------------------------------------------
#' @rdname summary
#' @importFrom stats quantile printCoefmat
#'
#' @param x an object of class \code{"summary.ptsr"},
#' usually, a result of a call to \code{summary.ptsr}.
#' @param digits  minimal number of significant digits, see
#' \code{\link{print.default}}.
#' @param signif.stars logical. If \code{TRUE},
#' \sQuote{significance stars} are printed for each coefficient.
#'
#' @details
#' \code{print.summary.btsr} tries to be smart about formatting the
#' coefficients, standard errors, etc. and additionally provides
#' \sQuote{significance stars}.
#'
#' @export
print.summary.ptsr <- function (x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),	...)
{

  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat("\n")
  cat("-----------------------------------------------")
  cat("\nCall:\n",
      "Positive Time Serie Regression Model", "\n\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L)
      structure(apply(t(resid), 1L, quantile),
                dimnames = list(nam, dimnames(resid)[[2L]]))
    else  {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  }
  else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  } else { # rdf == 0 : perfect fit!
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  ##
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  cat("-----------------------------------------------\n")
  cat("\n")
  invisible(x)
}
