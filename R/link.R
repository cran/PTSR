#------------------------------------------------------------
#  Similar to make.link()
#------------------------------------------------------------
##'  Given the name of a link, this function returns a link function,
##'  an inverse link function, the derivative   \eqn{d\eta / d\mu}{deta/dmu}
##'  and the derivative \eqn{d\mu / d\eta}{dmu/deta}.
##'
##' @title Create a Link for PTSR models
##'
##' @param link  character; one of  \code{"log"}, \code{"log1"}. See \sQuote{Details}.
##'
##' @return An object of class \code{"link-ptsr"}, a list with components
##'
##'  \item{linkfun}{Link function \code{function(mu)}}
##'  \item{linkinv}{Inverse link function \code{function(eta)}}
##'  \item{linkdif}{Derivative \code{function(mu)} \eqn{d\eta / d\mu}{deta/dmu}}
##'  \item{mu.eta}{Derivative \code{function(eta)} \eqn{d\mu / d\eta}{dmu/deta}}
##'  \item{name}{a name to be used for the link}
##'
##'@details The available links are:
##'
##'  log:    \eqn{f(x) = log(x)}
##'
##'  log1: \eqn{f(x) = log(x-1)}
##'
##' @export
##' @md

ptsr.link <- function(link = "log"){
  ##--------------------------------------------------
  ## linkfun:  Link function function(mu)
  ## linkinv:  Inverse link function function(eta)
  ## mu.eta:   Derivative function(eta) dmu/deta
  ## diflink:  Derivative function(mu) deta/dmu
  ##--------------------------------------------------

	if(!(link %in% c("log", "log1", "identity")))
		stop(cat("Unknown link:", link, "\n"))

  # defines g(mu)
  linkfun <- function(mu){
    switch(link,
           log = log(mu),
           log1 = log(mu - 1),
           identity = mu
    )
  }

  # defines g^{-1}(eta)
  linkinv <- function(eta){
    switch(link,
           log = exp(eta),
           log1 = exp(eta) + 1,
           identity = eta
    )
  }

  # defines dg/dmu
  diflink <- function(mu){
    switch(link,
           log = 1/mu,
           log1 = 1/(mu - 1),
           identity = 1
    )
  }

  # defines dmu/deta = 1/g'(ginv(eta))
  mu.eta <- function(eta){
    1/diflink(mu = linkinv(eta = eta))
  }

  structure(list(linkfun = linkfun,
                 linkinv = linkinv,
                 diflink = diflink,
                 mu.eta = mu.eta,
                 name = link), class = "link-ptsr")
}
