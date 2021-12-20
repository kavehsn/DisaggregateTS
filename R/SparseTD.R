#' High-dimensional temporal disaggregation  
#'
#' This function provides the Sparse Temporal Disaggregation (spTD) methods proposed by \insertCite{mosley2021sparse;textual}{DisaggregateTS}
#' to perform temporal disaggregation of time series data in both standard and high-dimensional settings. Variable selection is also
#' performed by a LASSO penalty \insertCite{tibshirani1996regression;textual}{DisaggregateTS} or an Adaptive LASSO penalty 
#' \insertCite{zou2006adaptive;textual}{DisaggregateTS}. 
#' 
#' @param Y       The low-frequency response vector. 
#' @param X       The high-frequency indicator matrix. 
#' @param penalty Nominates the choice of regularisation ('lasso' or 'adalasso').
#' @param aggMat  Aggregation matrix according to 'first', 'sum', 'average', 'last'. 
#' @keywords Sparse Temporal Disaggregation Lasso Time Series Disaggregation
#' @import lars 
#' @export
#' @examples 
#' data = TempDisaggDGP(n_l = 50, m = 4, p = 10, beta = 3, sparsity = 0.5, method = 'Chow-Lin', rho = 0.5)
#' X = data$X_Gen
#' Y = data$Y_Gen
#' fit_spTD = SparseTD(Y = Y, X = X, penalty = 'lasso', aggMat = 'sum')
#' y_hat = fit_spTD$y
#' @references 
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

### SparseTD main function --------------------

SparseTD <- function(Y, X = matrix(data = rep(1, times = nrow(Y)), nrow = nrow(Y)), penalty = 'lasso', aggMat = 'sum') {
  
  if(is.matrix(X) == FALSE || is.matrix(Y) == FALSE){
    
    stop("X and Y must be a matrices! \n")
    
  }else{
    
    
    if(nrow(X) %% nrow(Y) != 0){
      
      
      stop("The high-frequency series must be a integer multiple of the low-frequency series! \n")
      
      
    }else{
      
      
      # Find the integer multiple
      
      m <- nrow(X) / nrow(Y)
      
      n <- nrow(X)
      n_l <- nrow(Y)
      
      # Generate the aggregation matrix C
      
      if(aggMat == 'sum'){
        
        C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = m))
        
      }else if(aggMat == 'avg'){
        
        C <- kronecker(diag(n_l), matrix(data = n_l/n, nrow = 1, ncol = m))
        
      }else if(aggMat == 'first'){
        
        C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = m-1)), nrow = 1, ncol = m))
        
      }else if(aggMat == 'last'){
        
        C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = m-1), 1), nrow = 1, ncol = m))
        
      }      
      
      # Aggregate X 
      
      X_l <- C %*% X
      
      # Grid for AR parameters rho
      
      grid <- seq(0.01, 0.99, by = 0.01)
      
      ## Loop over grid and store BIC for each rho
      
      # Storage       
      min_bic <- lambdahat <- c()
      betahat <- list()
      
      for(rho in 1:length(grid)) {
        
        
        # Generate AR auto-covariance matrix 
        sqnc <- grid[rho]^seq(0, n, by = 1)
        Omega <- toeplitz(sqnc[1: n])
        V <- (1/(1-grid[rho]^2)) * Omega
        
        # Aggregate V 
        V_l <- C %*% tcrossprod(V, C)
        
        # Simplification and Cholesky factorization of the auto-covariance matrix 
        Uchol <- chol(V_l)
        Lchol <- t(Uchol)
        
        # Preconditioning the variables
        X_F <- solve(Lchol) %*% X_l
        Y_F <- solve(Lchol) %*% Y
        
        # Fit LARS algorithm to the data 
        lars.fit <- lars(X_F, Y_F, intercept = F, normalize = F)
        betamat <- lars.fit$beta 
        
        # Don't allow support to be bigger than n_l/2
        npath <- k.index(betamat, n_l)
        
        # Find BIC for each re-fitted betahat 
        beta_refit <- list()
        BIC <- c()
        BIC[1] <- hdBIC(X_F, Y_F, V_l, betamat[1,])
        beta_refit[[1]] <- betamat[1,]
        
        for(lam in 2:npath) {
          
          beta_refit[[lam]] <- refit(X_F, Y_F, betamat[lam, ])
          BIC[lam] <- hdBIC(X_F, Y_F, V_l, beta_refit[[lam]])
          
        }
        
        # Store best BIC, betahat and lambdahat 
        min_bic[rho] <- min(BIC)
        min_bic_idx <- which.min(BIC)
        betahat[[rho]] <- beta_refit[[min_bic_idx]]
        lambdahat[rho] <- lars.fit$lambda[min_bic_idx]
        
        
      }
      
      # Store the final parameter estimates for lasso regularisation
      rho_idx <- which.min(min_bic)
      rhohat_final <- grid[rho_idx]
      betahat_final <- betahat[[rho_idx]]
      lambdahat_final <- lambdahat[rho_idx]
      
      # Generate AR auto-covariance matrix 
      sqnc <- rhohat_final^seq(0, n, by = 1)
      Omega <- toeplitz(sqnc[1: n])
      V <- (1/(1-rhohat_final^2)) * Omega
      
      # Aggregate V
      V_l <- C %*% V %*% t(C)
      
      yhat_final <- X %*% betahat_final + V %*% t(C) %*% solve(V_l) %*% (Y - X_l %*% betahat_final)
      
      # If adaptive lasso, do the following 
      if(penalty == 'adalasso') {
        
        # get adaptive lasso weight
        ada_weight <- abs(betahat_final)
        
        
        # Simplification and Cholesky factorization of the auto-covariance matrix 
        Uchol <- chol(V_l)
        Lchol <- t(Uchol)
        
        # Preconditioning the variables
        X_F <- solve(Lchol) %*% X_l
        Y_F <- solve(Lchol) %*% Y
        
        # Scale the data A*|betahat_lasso|
        X_F_scaled <- scale(X_F, center = F, scale = 1/ada_weight)
        
        # Apply lars to scaled data 
        lars.fit <- lars(X_F_scaled, Y_F, intercept = F, normalize = F)
        betamat <- lars.fit$beta 
        
        # Scale beta estimates back to original 
        betamat_scaled <- scale(betamat, center = F, scale = 1/ada_weight)
        
        # Tune with BIC 
        npath <- k.index(betamat_scaled, n_l)
        beta_refit <- list()
        BIC <- c()
        BIC[1] <- hdBIC(X_F, Y_F, V_l, betamat_scaled[1,])
        beta_refit[[1]] <-  betamat_scaled[1,]
        
        for(lam in 2:npath) {
          
          beta_refit[[lam]] <- refit(X_F, Y_F, betamat_scaled[lam, ])
          BIC[lam] <- hdBIC(X_F, Y_F, V_l, beta_refit[[lam]])
          
        }
        
        # Store best BIC, betahat and lambdahat 
        min_bic_idx <- which.min(BIC)
        betahat_final <- beta_refit[[min_bic_idx]]
        lambdahat_final <- lars.fit$lambda[min_bic_idx]
        
        yhat_final <- X %*% betahat_final + V %*% t(C) %*% solve(V_l) %*% (Y - X_l %*% betahat_final)
        
        
      }
      
      
    }
    
  }   
  
  output <- list('y' = yhat_final, 'beta' = betahat_final, 'rho' = rhohat_final, 'lambda' = lambdahat_final)
  
  return(output)
  
}
