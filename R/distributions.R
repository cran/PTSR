############################################################
# Reparametrized distributions
# Beta-Prime:           d.betap         r.betap
# F:                    d.F             r.F
# Gamma:                d.gamma         r.gamma
# Inverse Gaussian:     d.invGauss      r.invGauss
# Log-Logistic:         d.logLogis      r.logLogis
# Log-Normal:           d.logNorm       r.logNorm
# Chi-Squared:          d.chi           r.chi
# Rayleigh:             d.ray           r.ray
############################################################
#----------------------------------------------------------
#' @title
#' Reparametrized Distributions
#'
#' @name ddist
#' @aliases rdist
#'
#' @order 1
#'
#' @description
#' Density function and random numbers generation for models with support on the positive real line.
#'
#' @details
#'
#' @return
#' For any avaliable \code{dist}, \code{ddist} gives the density and \code{rdist} generates random deviates.
#'
#' The length of the result is determined by \code{n} for \code{rdist}, and is the maximum of the lengths of the numerical arguments for \code{rdist}.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @md
NULL
#> NULL



############################################################# Reparametrized Beta-Prime Distribution
############################################################
#' @rdname ddist
#' @aliases d.betap
#' @order 2
#'
#' @param x vector of real values
#' @param mu non-negative parameter (the distribution's mean. See \sQuote{Details})
#' @param varphi non-negative parameter
#' @param log logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}.
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Beta-Prime distribution, the functions [dbetapr][extraDistr::dbetapr] and [rbetapr][extraDistr::rbetapr] are imported from the package \code{\link{extraDistr}}. The following holds
#'   \deqn{shape1 = mu*varphi}
#'   \deqn{shape2 = varphi + 1}
#'   \deqn{scale = 1}
#'   }
#'
#' @export
#' @md
d.betap <- function(x, mu, varphi, log = FALSE){
  shape1 = mu*varphi
  shape2 = varphi + 1
  return(extraDistr::dbetapr(x = x, shape1 = shape1, shape2 = shape2, scale = 1, log = log))
}

#' @rdname ddist
#' @aliases r.betap
#' @order 3
#' @param n sample size
#' @export
#' @md
r.betap <- function(n, mu, varphi){
  shape1 = mu*varphi
  shape2 = varphi + 1
  return(extraDistr::rbetapr(n = n, shape1 = shape1, shape2 = shape2, scale = 1))
}

############################################################
# Reparametrized F Distribution
############################################################
#' @rdname ddist
#' @aliases d.F
#' @order 4
#'
#' @details
#' \itemize{
#'   \item For the reparametrized F distribution, the functions [df][stats::df] and  [rf][stats::rf] are imported from \code{\link{stats}}. The following holds
#'   \deqn{df1 = varphi}
#'   \deqn{df2 = 2*mu/(mu - 1)}
#'   so that the parameter \eqn{\mu} must satisfy \eqn{\mu > 1}.
#'   }
#'
#' @export
#' @md
d.F <- function(x, mu, varphi, log = FALSE){
  df1 = varphi
  df2 = 2*mu/(mu - 1)
  return(stats::df(x = x, df1 = df1, df2 = df2, log = log))
}

#' @rdname ddist
#' @aliases r.F
#' @order 5
#' @export
#' @md
r.F <- function(n, mu, varphi){
  df1 = varphi
  df2 = 2*mu/(mu - 1)
  return(stats::rf(n = n, df1 = df1, df2 = df2))
}


############################################################
# Reparametrized Gamma Distribution
############################################################
#' @rdname ddist
#' @aliases d.gamma
#' @order 6
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Gamma distribution, the functions [dgamma][stats::dgamma] and  [rgamma][stats::rgamma] are imported from \code{\link{stats}}. The following holds
#'   \deqn{shape = varphi}
#'   \deqn{rate = varphi/mu}
#'   }
#'
#' @export
#' @md
d.gamma <- function(x, mu, varphi, log = FALSE){
  shape = varphi
  rate = varphi/mu
  return(stats::dgamma(x = x, shape = shape, rate = rate, log = log))
}

#' @rdname ddist
#' @aliases r.gamma
#' @order 7
#' @export
#' @md
r.gamma <- function(n, mu, varphi) {
  shape = varphi
  rate = varphi/mu
  return(stats::rgamma(n = n, shape = shape, rate = rate))
}


############################################################
# Reparametrized Inverse Gaussian Distribution
############################################################
#' @rdname ddist
#' @aliases d.invGauss
#' @order 8
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Inverse Gaussian distribution, the functions [dinvGauss][SuppDists::dinvGauss] and [rinvGauss][SuppDists::rinvGauss] are imported from \code{SuppDists}. The following holds
#'   \deqn{nu = mu}
#'   \deqn{lambda = 1/varphi}
#'   }
#'
#' @export
#' @md
d.invGauss <- function(x, mu, varphi, log = FALSE){
  nu = mu
  lambda = 1/varphi
  return(SuppDists::dinvGauss(x = x, nu = nu, lambda = lambda, log = log))
}

#' @rdname ddist
#' @aliases r.invGauss
#' @order 9
#' @export
#' @md
r.invGauss <- function(n, mu, varphi){
  nu = mu
  lambda = 1/varphi
  return(SuppDists::rinvGauss(n = n, nu = nu, lambda = lambda))
}


############################################################
# Reparametrized Log-logistic Distribution
############################################################
#' @rdname ddist
#' @aliases d.logLogis
#' @order 10
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Log-logistic distribution, the functions [dllogis][actuar::dllogis] and [rllogis][actuar::rllogis] a are imported from \code{actuar}. The following holds
#'   \deqn{shape = varphi}
#'   \deqn{rate = (pi/varphi)/(mu*sin(pi/varphi))}
#'   }
#'
#' @export
#' @md
d.logLogis <- function(x, mu, varphi, log = FALSE){
  shape = varphi
  rate = (pi/varphi)/(mu*sin(pi/varphi))
  return(actuar::dllogis(x = x, shape = shape, rate = rate, log = log))
}

#' @rdname ddist
#' @aliases r.logLogis
#' @order 11
#' @export
#' @md
r.logLogis <- function(n, mu, varphi){
  shape = varphi
  rate = (pi/varphi)/(mu*sin(pi/varphi))
  return(actuar::rllogis(n = n, shape = shape, rate = rate))
}


############################################################
# Reparametrized Log-Normal Distribution
############################################################
#' @rdname ddist
#' @aliases d.logNorm
#' @order 12
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Log-Normal distribution, the functions [dlnorm][stats::dlnorm] and [rlnorm][stats::rlnorm] are imported from \code{\link{stats}}. The following holds
#'   \deqn{meanlog = log(mu) - varphi^2/2}
#'   \deqn{sdlog = varphi}
#'   }
#'
#' @export
#' @md
d.logNorm <- function(x, mu, varphi, log = FALSE){
  meanlog = log(mu) - varphi^2/2
  sdlog = varphi
  return(stats::dlnorm(x = x, meanlog = meanlog, sdlog = sdlog, log = log))
}

#' @rdname ddist
#' @aliases r.logNorm
#' @order 13
#' @export
#' @md
r.logNorm <- function(n, mu, varphi){
  meanlog = log(mu) - varphi^2/2
  sdlog = varphi
  return(stats::rlnorm(n = n, meanlog = meanlog, sdlog = sdlog))
}


############################################################
# Reparametrized Chi-squared Distribution
############################################################
#' @rdname ddist
#' @aliases d.chi
#' @order 14
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Chi-squared F distribution, the functions [dchisq][stats::dchisq] and [rchisq][stats::rchisq] are imported from \code{\link{stats}}. The following holds
#'   \deqn{df = mu}
#'   }
#'
#' @export
#' @md
d.chi <- function(x, mu, log = FALSE,...){
  df = mu
  return(stats::dchisq(x = x, df = df, log = log))
}

#' @rdname ddist
#' @aliases d.chi
#' @order 15
#' @param ... for compatibility with other functions
#' @export
#' @md
r.chi <- function(n, mu, ...){
  df = mu
  return(stats::rchisq(n = n, df = df))
}


############################################################
# Reparametrized Rayleigh Distribution
############################################################
#' @rdname ddist
#' @aliases d.ray
#' @order 16
#'
#' @details
#' \itemize{
#'   \item For the reparametrized Rayleigh distribution, the functions [drayleigh][extraDistr::drayleigh] and [rrayleigh][extraDistr::rrayleigh] are imported from \code{\link{extraDistr}}. The following holds
#'   \deqn{sigma = mu/sqrt(pi/2)}
#'   }
#'
#' @export
#' @md
d.ray <- function(x, mu, log = FALSE,...){
  sigma = mu/sqrt(pi/2)
  return(extraDistr::drayleigh(x = x, sigma = sigma, log = log))
}

#' @rdname ddist
#' @aliases r.ray
#' @order 17
#' @export
#' @md
r.ray <- function(n, mu, ...){
  sigma = mu/sqrt(pi/2)
  return(extraDistr::rrayleigh(n = n, sigma = sigma))
}

