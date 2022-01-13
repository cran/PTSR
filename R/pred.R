#' @title Predict method for PTSR
#'
#' @description Predicted values based on ptsr object.
#'
#' @param object Object of class inheriting from \code{"ptsr"}
#' @param newdata A matrix with new values for the regressors.   If omitted
#' and \code{"xreg"} is present in the model, the fitted values are returned.
#' If the model does not include regressors, the functions will use
#' the value of \code{nnew}.
#' @param nnew number of out-of-sample forecasts required. If \code{newdata} is
#' provided, \code{nnew} is ignored.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'   \code{predict.ptsr} produces predicted values, obtained by evaluating
#'   the regression function in the frame \code{newdata}.
#'
#'   If \code{newdata} is omitted the predictions are based on the data
#'   used for the fit.
#'
#'   For now, prediction intervals are not provided.
#'
#' @return A list with the following arguments
#'
#'  \item{series}{The original time series yt.}
#'
#'  \item{xreg}{The original regressors (if any).}
#'
#'  \item{fitted.values}{The in-sample forecast given by \eqn{\mu_t}.}
#'
#'  \item{etat}{In-sample values of \eqn{g(\mu[t])}.}
#'
#'  \item{error}{The error term}
#'
#'  \item{residuals}{The (in-sample) residuals, that is, the observed minus the predicted values.}
#'
#'  \item{forecast}{The predicted values for yt.}
#'
##' @export
predict.ptsr <-
  function(object, newdata, nnew = 0,...){

    out <- list()
    nms.out <- c("series", "xreg", "fitted.values", "residuals")

    if(missing(newdata)) newdata = NULL

    if(is.null(newdata) & nnew <= 0){
      #------------------------------------------------------
      # New data was not provided.
      # Extracting existing components and returning
      #------------------------------------------------------
      out[nms.out] <- object[nms.out]
    }else{

      if(is.null(newdata) & object$model$b > 0)
        stop("Please, provide the new values for the regressors")

      #------------------------------------------------------
      # New data was provided.
      # Making the necessary calculations
      #------------------------------------------------------
      xnew = NULL
      if(!is.null(newdata)){
        xnew = as.matrix(newdata)
        nnew = nrow(xnew)
      }

      temp <- .ptsr.predict(par = object$coefficients,
                            h = nnew, xnew = xnew, yt = object$series,
                            mut = object$fitted.values,
                            epst = object$residuals,
                            model = object$model)


      out[nms.out] <- object[nms.out]
      out$forecast <- temp
      out$xnew <- NULL
      if(object$model$b > 0) out$xnew <- xnew
    }
    out
  }


# Internal function.
# Calculates out-of-sample predicted values
.ptsr.predict <- function(par, h, xreg, xnew, yt, mut, epst, model) {

  n = length(yt)
  mut = c(mut, numeric(h))
  epst  = c(epst, numeric(h))

  p = model$p
  q = model$q
  a = model$a
  b = model$b
  arlag = model$arlag
  malag = model$arlag
  g1 = model$g1
  g1.inv = model$g1.inv
  g2 = model$g2

  g2y = c(g2(yt), numeric(h))
  yt = c(yt, numeric(h))

  alpha = 0
  beta = 0
  l = 0
  if (a == 1) {alpha = par[1]; l = l + 1}
  if (b > 0) {beta = par[(l + 1):(l + b)]; l = l + b}
  if (p > 0) {
    if (is.null(arlag)) {
      phi = par[(l + 1):(l + p)]
      l = l + p
    }
    else {
      phi = numeric(p)
      phi[arlag] = par[(l + 1):(l + length(arlag))]
      l = l + length(arlag)
    }
  }
  if (q > 0) {
    if (is.null(arlag)) {
      theta = par[(l + 1):(l + q)]
      l = l + q
    }
    else{
      theta = numeric(q)
      theta[malag] = par[(l + 1):(l + length(malag))]
      l = l + length(malag)
    }
  }

  XB = numeric(n+h)
  if (!is.null(xnew)){
	  xtemp = rbind(as.matrix(xreg), as.matrix(xnew))
		XB = xtemp %*% beta
	}

  for (t in 1:h) {
    ls = alpha + XB[n+t]
    if (p > 0) {
      xr = 0
      if (model$xregar) xr = XB[(n + t) - c(1:p)]
      temp = sum(phi * (g2y[(n + t) - c(1:p)] - xr))
      ls = ls + temp
    }
    if (q > 0) ls = ls + sum(theta * epst[(n + t) - (1:q)])
    mut[n + t] = g1.inv(ls)
    yt[n + t] = mut[n + t]
    g2y[n+t] = g2(mut[n + t])
  }

  ynew = yt[(n + 1):(n + h)]

  invisible(ynew)
}

