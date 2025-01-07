#' Temporal Disaggregation Methods
#'
#' This function contains the traditional standard-dimensional temporal disaggregation methods proposed by \insertCite{denton1971adjustment;textual}{DisaggregateTS}, \insertCite{dagum2006benchmarking;textual}{DisaggregateTS},
#' \insertCite{chow1971best;textual}{DisaggregateTS}, \insertCite{fernandez1981methodological;textual}{DisaggregateTS} and \insertCite{litterman1983random;textual}{DisaggregateTS},
#' and the high-dimensional methods of \insertCite{10.1111/rssa.12952;textual}{DisaggregateTS}.
#'
#' Takes in a \eqn{n_l \times 1} low-frequency series to be disaggregated \eqn{Y} and a \eqn{n \times p} high-frequency matrix of p indicator series \eqn{X}. If \eqn{n > n_l \times aggRatio} where \eqn{aggRatio}
#' is the aggregation ratio (e.g. \eqn{aggRatio = 4} if annual-to-quarterly disagg, or \eqn{aggRatio = 3} if quarterly-to-monthly disagg) then extrapolation is done
#' to extrapolate up to \eqn{n}.
#'
#' @param Y  		   The low-frequency response series (\eqn{n_l \times 1} matrix).
#' @param X  		   The high-frequency indicator series (\eqn{n \times p} matrix).
#' @param aggMat 	 Aggregation matrix according to 'first', 'sum', 'average', 'last' (default is 'sum').
#' @param aggRatio Aggregation ratio e.g. 4 for annual-to-quarterly, 3 for quarterly-to-monthly (default is 4).
#' @param method   Disaggregation method using 'Denton', 'Denton-Cholette', 'Chow-Lin', 'Fernandez', 'Litterman', 'spTD' or 'adaptive-spTD' (default is 'Chow-Lin').
#' @param Denton	 Type of differencing for Denton method: 'simple-diff', 'additive-first-diff', 'additive-second-diff', 'proportional-first-diff' and 'proportional-second-diff' (default is 'additive-first-diff'). For instance, 'simple-diff' differencing refers to the differences between the original and revised values, whereas 'additive-first-diff' differencing refers to the differences between the first differenced original and revised values.
#' @return \code{y_Est}:	  Estimated high-frequency response series (output is an \eqn{n \times 1} matrix).
#' @return \code{beta_Est}:	Estimated coefficient vector (output is a \eqn{p \times 1} matrix).
#' @return \code{rho_Est}:	Estimated residual AR(1) autocorrelation parameter.
#' @return \code{ul_Est}:	  Estimated aggregate residual series (output is an \eqn{n_l \times 1} matrix).
#' @keywords Denton Denton-Cholette Chow-Lin Fernandez Litterman spTD adaptive-spTD lasso temporal-disaggregation
#' @import Matrix lars
#' @export
#' @examples
#' data <- TempDisaggDGP(n_l=25,n=100,p=10,rho=0.5)
#' X <- data$X_Gen
#' Y <- data$Y_Gen
#' fit_chowlin <- disaggregate(Y=Y,X=X,method='Chow-Lin')
#' y_hat = fit_chowlin$y_Est
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm optimize




disaggregate <- function(Y, X = matrix(data = rep(1, times = (nrow(Y)*aggRatio)), nrow = (nrow(Y)*aggRatio)), aggMat = 'sum', aggRatio = 4, method = 'Chow-Lin', Denton = 'additive-first-diff'){


  if(is.matrix(X) == FALSE || is.matrix(Y) == FALSE){

    stop("X and Y must be a matrices! \n")

  }

  if(!(method == 'Denton' || method == 'Denton-Cholette' || method == 'Chow-Lin' || method == 'Fernandez' || method == 'Litterman' || method == 'spTD' || method == 'adaptive-spTD')) {

    stop("Wrong method inputted \n")

  }

  n_l = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  nfull = aggRatio*n_l
  extr = n - nfull # number of extrapolations

  if(nfull > n) {

    stop("X does not have enough observations. \n")

  }


  # Generate the aggregation matrix C

  if(aggMat == 'sum'){

    C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))

  }else if(aggMat == 'average'){

    C <- kronecker(diag(n_l), matrix(data = 1/aggRatio, nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))

  }else if(aggMat == 'first'){

    C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = aggRatio-1)), nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))

  }else if(aggMat == 'last'){

    C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = aggRatio-1), 1), nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))

  }

  X_l = C %*% X


  if(method == 'Denton-Cholette' || method == 'Denton'){

    if(dim(as.matrix(X))[2] > 1){
      stop("X has more than 1 indicator. The Denton/Denton-Cholette methods requires only one indicator. \n")
    }

    if(!(Denton == 'simple-diff' || Denton == 'additive-first-diff' || Denton == 'additive-second-diff' || Denton == 'proportional-first-diff' || Denton == 'proportional-second-diff')){
      stop("Wrong Denton method inputted.")
    }

    # Difference between the low-frequency and the transformed high-frequency series.

    u_l <- Y - X_l

    # First difference matrix

    Delta_0 <- diag(n)

    diags <- list(rep(1, times = n), rep(-1, times = n-1))
    Delta_t <- as.matrix(bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE))  # Convert to dense matrix
    Delta <- t(Delta_t)

    X_inv <- diag(1 / (as.numeric(X) / mean(X)))

    if(Denton == 'simple-diff'){
      Delta_0 = Delta_0
    }else if(Denton == 'additive-first-diff'){
      Delta_0 = Delta
    }else if(Denton == 'additive-second-diff'){
      Delta_0 = t(Delta) %*% Delta
    }else if(Denton == 'proportional-first-diff'){
      Delta_0 = Delta %*% X_inv
    }else{
      Delta_0 = t(Delta) %*% Delta %*% X_inv
    }

    if(method == 'Denton'){

      vcov = solve(t(Delta_0) %*% Delta_0)

      # The distribution matrix

      D <- vcov %*% t(C) %*% solve(C %*% vcov %*% t(C))

      # Generate the high-frequency series

      y <- X + (D %*% u_l)

      # Unnecessary parameter outputs

      rho_opt <- NaN
      betaHat <- NaN

    }else if(method == 'Denton-Cholette'){

      if(Denton == 'simple'){
        Delta_1 = Delta_0
      }else if(Denton == 'additive-first-diff' || Denton == 'proportional-first-diff'){
        Delta_1 = Delta_0[-1,]
      }else{
        Delta_1 = Delta_0[-(1:2),]
      }

      D1 = t(Delta_1) %*% Delta_1

      # Eq. (2.2) from Denton (1971); Eq (6.8) from Cholette and Dagum (2006)
      y <- solve(
        rbind(cbind(D1, t(C)), cbind(C, matrix(0, nrow = n_l, ncol = n_l)))
      ) %*% rbind(
        cbind(D1, matrix(0, nrow = n, ncol = n_l)),
        cbind(C, diag(1, nrow = n_l, ncol = n_l))
      ) %*% matrix(c(X, u_l))

      y <- y[1:n]

      # Unnecessary parameter outputs

      rho_opt <- NaN
      betaHat <- NaN

    }

  }else if(method == 'Fernandez') {


    diags <- list(rep(1, times = n), rep(-1, times = n-1))
    Delta_t <- as.matrix(bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE))  # Convert to dense matrix
    Delta <- t(Delta_t)


    # H(rho) matrix is first an nxn identity matrix

    vcov <- solve(Delta_t %*% Delta)
    vcov_agg <- C %*% vcov %*% t(C)

    # Simplification and Cholesky factorization of the Sigma_opt

    Uchol <- chol(vcov_agg)
    Lchol <- t(Uchol)

    # Preconditioning the variables

    X_F <- solve(Lchol) %*% X_l
    Y_F <- solve(Lchol) %*% Y

    # First estimate betaHat_opt using OLS assuming Sigma = (Delta'Delta)^{-1}

    betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F

    # Obtain the residuals using betaHat_1

    u_l <- Y - X_l %*% betaHat

    # The distribution matrix

    D <- vcov %*% t(C) %*% solve(vcov_agg)

    # Generate the high-frequency series

    y <- X %*% betaHat + (D %*% u_l)

    rho_opt <- NaN

  }else {

    if(method == 'Chow-Lin') {

      Objective <- function(rho) {
        -chowlin_likelihood(
          Y = Y, X = X_l, vcov = C %*% ARcov(rho,n) %*% t(C)
        )
      }
    }else if(method == 'Litterman') {

      Objective <- function(rho) {
        -as.numeric(chowlin_likelihood(
          Y = Y, X = X_l, vcov = forceSymmetric(C %*% ARcov_lit(rho,n) %*% t(C))
        ))
      }

    }else if(method == 'spTD' || method == 'adaptive-spTD') {

      Objective <- function(rho) {
        sptd_BIC(
          Y = Y, X = X_l, vcov = C %*% ARcov(rho,n) %*% t(C)
        )
      }

    }

    # optimise for best rho

    optimise_rho <- optimize(Objective,
                             lower = 0, upper = 0.999, tol = 1e-16,
                             maximum = FALSE
    )

    rho_opt = optimise_rho$minimum

    # Generate the optimal Toeplitz covariance matrix

    if(method == 'Chow-Lin'){

      fit = chowlin(Y = Y, X = X, rho = rho_opt, aggMat = aggMat, aggRatio = aggRatio, litterman = FALSE)
      betaHat = fit$betaHat
      y = fit$y
      u_l = fit$u_l

    }else if(method == 'Litterman'){

      fit = chowlin(Y = Y, X = X, rho = rho_opt, aggMat = aggMat, aggRatio = aggRatio, litterman = TRUE)
      betaHat = fit$betaHat
      y = fit$y
      u_l = fit$u_l

    }else if(method == 'spTD'){

      fit = sptd(Y = Y, X = X, rho = rho_opt, aggMat = aggMat, aggRatio = aggRatio, adaptive = FALSE)
      betaHat = fit$betaHat
      y = fit$y
      u_l = fit$u_l

    }else if(method == 'adaptive-spTD'){

      fit = sptd(Y = Y, X = X, rho = rho_opt, aggMat = aggMat, aggRatio = aggRatio, adaptive = TRUE)
      betaHat = fit$betaHat
      y = fit$y
      u_l = fit$u_l

    }

  }

  data_list <- list(y, betaHat, rho_opt, u_l)
  names(data_list) <- c("y_Est", "beta_Est", "rho_Est","ul_Est")

  return(data_list)
}
