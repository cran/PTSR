#' @title
#' Function to simulate a PTSR model
#'
#' @details
#' The function \code{ptsr.sim} generates a random sample from a positive time
#' series regression model, with a given distribution.
#'
#' @param n a strictly positive integer. The sample size of yt (after burn-in).
#'  Default is 1.
#' @param burn a non-negative integer. length of "burn-in" period. Default is 0.
#' @param xreg optionally, a vector or matrix of external regressors.
#'  For simulation purposes, the length of xreg must be \code{n+burn}.
#'  Default is \code{NULL}.
#' @param xregar logical, if \code{FALSE}, the regressors are not included in the
#' AR component of the model. Default is \code{TRUE}.
#' @param varphi non-negative parameter. Default is 1.
#' @param alpha a numeric value corresponding to the intercept.  Default is 0.
#' @param beta optionally, a vector of coefficients corresponding to the
#'   regressors in \code{xreg}. Default is \code{NULL}.
#' @param phi optionally, for the simulation function this must be a vector
#'   of size \eqn{p}, corresponding to the autoregressive coefficients
#'   (including the ones that are zero), where \eqn{p} is the AR order.
#'   Default is \code{NULL}.
#' @param theta optionally, for the simulation function this must be a vector
#'   of size \eqn{q}, corresponding to the moving average coefficients
#'   (including the ones that are zero), where \eqn{q} is the MA order.
#'   Default is \code{NULL}.
#'  that \eqn{g_2(y_t) = 0}, for \eqn{t < 1}.
#' @param seed optionally, an integer which gives the value of the fixed
##'  seed to be used by the random number generator. If missing, a random integer
##'  is chosen uniformly from 1,000 to 10,000.
#' @param rdist function, the random number generator to be used
#' @param link1 character indicating which link must be used for \eqn{\mu_t}.  See \code{\link{ptsr.link}} for available links. Default is \sQuote{log}.
#' @param link2 character indicating which link must be used for \eqn{y_t} in the AR recursion.  See \code{\link{ptsr.link}} for available links. Default is \sQuote{identity}.
#' @param g1 optionally, a link function to be used  for \eqn{\mu_t}. Default is \code{NULL}, so that it is calculated internally, using \code{link1}.
#' @param g1.inv optionally, a the inverse link function to be used  for \eqn{\eta_t}. It must the the ivnerse of \code{g1}. Default is \code{NULL}, so that it is calculated internally, using \code{link1}.
#' @param g2 optionally, a link function to be used  for \eqn{y_t}. Default is \code{NULL}, so that it is calculated internally, using \code{link2}.
#'
#' @return
#' Returns a list with the following components
#' \itemize{
#'
#' \item \code{yt}: the simulated time series
#'
#' \item \code{mut}: the conditional mean
#'
#' \item \code{etat}: the linear predictor \eqn{g(\mu_t)}
#'
#' \item \code{error}: the error term.
#'
#' }
#'
#' @examples
#'
#' #-------------------------------------------------------------------
#' # Generating a sample of a Gamma-ARMA(1,1) model with no regressors
#' #-------------------------------------------------------------------
#'
#' simu = ptsr.sim(n = 300, burn = 50,
#'                 varphi = 20, alpha = 0,
#'                 phi = 0.35, theta = 0.2,
#'                 seed = 1234, rdist = r.gamma,
#'                 link1 = "log", link2 = "log")
#'
#' names(simu)
#' plot.ts(simu$yt)
#' lines(simu$mut, col= "red")
#'
#' @export
ptsr.sim <- function(n = 1, burn = 0, xreg = NULL, xregar = TRUE,
                     varphi = 1, alpha = 0, beta = NULL,
                     phi = NULL, theta = NULL,
                     seed = stats::runif(1, 1000, 10000), rdist = r.gamma,
                     link1 = "log", link2 = "identity",
                     g1 = NULL, g1.inv = NULL, g2 = NULL) {

  # function to be applied to mu and eta
  if (is.null(g1) | is.null(g1.inv)) {
    lk1 <- ptsr.link(link = link1)
    g1 <- lk1$linkfun
    g1.inv <- lk1$linkinv
  }

  # function to be applied to r2
  if (is.null(g2)) {
    lk2 <- ptsr.link(link = link2)
    g2 <- lk2$linkfun
  }

  set.seed(seed)

  p <- length(phi)
  q <- length(theta)

  n <- n + burn
  # vectors to save the time series
  yt <- numeric(n+p)
  gyt <- numeric(n+p)
  mut <- numeric(n)
  epst <- numeric(n+q)
  etat <- numeric(n)

  # checking beta and xreg
  if (is.null(beta) | is.null(xreg)) {
    XB <- matrix(0, nrow = n, ncol = 1)
  }
  else {
    if (is.null(xreg)) stop("missing xreg")
    xreg <- as.matrix(xreg)
    if (is.null(beta)) stop("missing beta")
    if (length(beta) != ncol(xreg)) stop("beta is not compatible with xreg")
    XB <- xreg %*% beta
  }
  if (nrow(XB) < n) stop("xreg must have n + burn rows")
  if(p > 0){
    # setting the initial p values as the unconditional mean
    if(xregar)	XB = c(rep(mean(XB[1:p]),p), XB)
  }

  # starting the loop
  for (t in 1:n) {

    etat[t] = alpha + XB[t]
    if(p > 0){
      iAR = (t+p) - (1:p)
      xr = 0
      if(xregar) xr = XB[iAR]
      etat[t] = etat[t] + sum(phi * (gyt[iAR] - xr))
    }
    if(q > 0) etat[t] = etat[t] + sum(theta * epst[(t+q) - (1:q)])

    mut[t] = g1.inv(etat[t])
    yt[t+p] = rdist(n = 1, mu = mut[t], varphi = varphi)
    gyt[t+p] = g2(yt[t+p])

    epst[t+q] = yt[t+p] - mut[t]
  }

	if(sum(yt[(p+burn + 1):(n+p)] <= 0) > 0) warning("Error generating the process. Zeros produced")

  invisible(list(yt = yt[(p+burn + 1):(n+p)],
                 mut = mut[(burn + 1):n],
                 etat = etat[(burn + 1):n],
                 epst = epst[(q+burn + 1):(n+q)]))
}
