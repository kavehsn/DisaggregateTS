#' Low-dimensional temporal dissaggregation toolbox
#' 
#' This function contains the movement preservation and regression-based low-dimensional temporal disaggregation methods proposed by \insertCite{denton1971adjustment;textual}{DisaggregateTS}, \insertCite{dagum2006benchmarking;textual}{DisaggregateTS}
#' \insertCite{chow1971best;textual}{DisaggregateTS}, \insertCite{fernandez1981methodological;textual}{DisaggregateTS} and \insertCite{litterman1983random;textual}{DisaggregateTS}.
#' 
#' @param Y  		The low-frequency response vector.
#' @param X  		The high-frequency indicator matrix.
#' @param method 	Disaggregation using 'Denton', 'Denton-Cholette', 'Chow-Lin', 'Fernandez', 'Litterman'.
#' @param aggMat 	Aggregation matrix according to 'first', 'sum', 'average', 'last'.
#' @param Denton	The 'absolute', 'first', 'second' and 'proportional' difference Sigma for the Denton method. 
#' @keywords Denton Denton-Cholette Chow-Lin Fernandez Litterman temporal-disaggregation
#' @import Matrix
#' @export
#' @examples
#' data = TempDisaggDGP(n_l=50,m=4,p=4,method='Chow-Lin',rho=0.5)
#' X = data$X_Gen
#' Y = data$Y_Gen
#' fit_chowlin = TempDisaggToolbox(Y=Y,X=X,method='Chow-Lin')
#' y_hat = fit_chowlin$y
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt	
#' @importFrom stats lm rbinom rnorm

TempDisaggToolbox <- function(Y, X = matrix(data = rep(1, times = nrow(Y)), nrow = nrow(Y)), method = 'Denton-Cholette', aggMat = 'sum', Denton = 'first'){


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
			
			# Generate the disaggregation matrix C

			if(aggMat == 'sum'){

				C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = m))

			}else if(aggMat == 'avg'){

				C <- kronecker(diag(n_l), matrix(data = n_l/n, nrow = 1, ncol = m))

			}else if(aggMat == 'first'){

				C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = m-1)), nrow = 1, ncol = m))

			}else if(aggMat == 'last'){

				C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = m-1), 1), nrow = 1, ncol = m))

			}

			if(method == 'Denton-Cholette' || method == 'Denton'){

				# Difference between the low-frequency and the transformed high-frequency series.

				u_l <- Y - C %*% X

				# First difference matrix

				diags <- list(rep(1, times = n), rep(-1, times = n-1))
				Delta_t <-  bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
				Delta <- t(Delta_t)

				if(method == 'Denton'){

					# Absolute difference matrix

					if(Denton == 'abs'){

						Sigma <- diag(n)

					# First difference Sigma

					}else if(Denton == 'first'){

						Sigma <- solve(Delta_t %*% Delta)

					# Second difference Sigma

					}else if(Denton == 'second'){

						Sigma <- solve(Delta_t %*% Delta_t %*% Delta %*% Delta) 

					# Proportional difference Sigma

					}else if(Denton == 'prop'){

						Sigma <- solve(solve(Diagonal(n, X)) %*% Delta_t %*% Delta %*% solve(Diagonal(n, X)))

					}

					# The distribution matrix

					D <- Sigma %*% t(C) %*% solve(C %*% Sigma %*% t(C))

					# Generate the high-frequency series

					y <- X + (D %*% u_l) 

					# Unnecessary parameter outputs

					rho_opt <- NaN
					betaHat_opt <- NaN

				}else if(method == 'Denton-Cholette'){

					# Removed the first roq of the Delta matrix

					Delta_DC <- Delta[2: nrow(Delta), ]

					# The distribution matrix

					D <- Delta_DC

					# Generate the high-frequency series

					y <- X + (D %*% u_l) 

					# Unnecessary parameter outputs

					rho_opt <- NaN
					betaHat_opt <- NaN

				}

			}else if(method == 'Chow-Lin'){


				# Find the AR(1) parameter corresponding that maximizes the likelihood function 

				rho <- 0
				counter <- 1

				LF <- rep(0, times = 99)

				while(rho < 1){

					# Generate the Toeplitz covariance matrix

					sqnc <- rho^seq(0, n, by = 1)
					Omega <- toeplitz(sqnc[1: n])
					Sigma <- (1/(1-rho^2)) * Omega

					# Simplification and Cholesky factorization of the Sigma 

					Uchol <- chol(C %*% Sigma %*% t(C))
					Lchol <- t(Uchol)

					# Preconditioning the variables

					X_F <- solve(Lchol) %*% C %*% X
					Y_F <- solve(Lchol) %*% Y  

					# Estimate betaHat_0 using GLS assuming Sigma with rho

					betaHat_0 <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F 

					# Obtain the residuals using betaHat_0

					u_l_sim <- Y - C %*% X %*% betaHat_0

					# Preconditioning for the LF function

					u_l_sim_F <- solve(Lchol) %*% u_l_sim

					# Calculate the likelihood function

					LF[counter] <- (-0.5 * t(u_l_sim_F) %*% u_l_sim_F)-(n_l/2)*log(2*pi)-1/2*log(det(Lchol %*% Uchol))
					
					counter <- counter + 1

					rho <- counter * 0.01

				}

				# Generate the optimal Toeplitz covariance matrix

				rho_opt <- (which.max(LF)-1) * 0.01

				sqnc_opt <- rho_opt ^ seq(0, n, by = 1)
				Omega_opt <- toeplitz(sqnc_opt[1: n])
				Sigma_opt <- (1/(1-rho_opt ^2)) * Omega_opt

				# Simplification and Cholesky factorization of the Sigma_opt 

				Uchol_opt <- chol(C %*% Sigma_opt %*% t(C))
				Lchol_opt <- t(Uchol_opt)

				# Preconditioning the variables

				X_F_opt <- solve(Lchol_opt) %*% C %*% X
				Y_F_opt <- solve(Lchol_opt) %*% Y  

				# Estimate betaHat_opt using GLS assuming Sigma is Toeplitz

				betaHat_opt <- solve(t(X_F_opt) %*% X_F_opt) %*% t(X_F_opt) %*% Y_F_opt 

				# The distribution matrix

				D <- Sigma_opt %*% t(C) %*% solve(C %*% Sigma_opt %*% t(C))

				# Obtain the residuals using betaHat_1

				u_l <- Y - C %*% X %*% betaHat_opt  

				# Generate the high-frequency series

				y <- X %*% betaHat_opt + (D %*% u_l)


			}else if(method == 'Litterman'|| method == 'Fernandez'){

				# First difference matrix

				diags <- list(rep(1, times = n), rep(-1, times = n-1))
				Delta_t <- bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
				Delta <- t(Delta_t) 

				if(method == 'Fernandez'){

					# H(rho) matrix is first an nxn identity matrix

					Sigma_opt <- solve(Delta_t %*% Delta)

					# Simplification and Cholesky factorization of the Sigma_opt 

					Uchol_opt <- chol(C %*% Sigma_opt %*% t(C))
					Lchol_opt <- t(Uchol_opt)

					# Preconditioning the variables

					X_F_opt <- solve(Lchol_opt) %*% C %*% X
					Y_F_opt <- solve(Lchol_opt) %*% Y  

					# First estimate betaHat_opt using OLS assuming Sigma = (Delta'Delta)^{-1}

					betaHat_opt <- solve(t(X_F_opt) %*% X_F_opt) %*% t(X_F_opt) %*% Y_F_opt 

					# Obtain the residuals using betaHat_1

					u_l <- Y - C %*% X %*% betaHat_opt

					# The distribution matrix

					D <- Sigma %*% t(C) %*% solve(C %*% Sigma %*% t(C))

					# Generate the high-frequency series

					y <- X %*% betaHat_opt + (D %*% u_l)

					rho_opt <- NaN

				}else if(method == 'Litterman'){

					# Find the AR(1) parameter corresponding to the errors of the random walk

					rho <- 0
					counter <- 1

					LF <- rep(0, times = 99)

					while(rho < 1){

						diags <- list(rep(1, times = n), rep(-rho, times = n-1))
						H_r_t <-  bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
						H_r <- t(H_r_t) 

						Sigma <- solve(Delta_t %*% H_r_t %*% H_r %*% Delta)

						# Simplification and Cholesky factorization of the Sigma 

						Uchol <- chol(C %*% Sigma %*% t(C))
						Lchol <- t(Uchol)

						# Preconditioning the variables

						X_F <- solve(Lchol) %*% C %*% X
						Y_F <- solve(Lchol) %*% Y   

						# Estimate betaHat_0 using OLS assuming Sigma with rho

						betaHat_0 <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F 

						# Obtain the residuals using betaHat_0

						u_l_sim <- Y - C %*% X %*% betaHat_0

						# Preconditioning for the LF function

						u_l_sim_F <- solve(Lchol) %*% u_l_sim

						# Calculate the likelihood function

						LF[counter] <- (-0.5 * t(u_l_sim_F) %*% u_l_sim_F)-(n_l/2)*log(2*pi)-1/2*log(det(Lchol %*% Uchol))

						counter <- counter + 1

						rho <- counter * 0.01

					}

					# Generate the optimal covariance matrix

					rho_opt <- (which.max(LF)-1) * 0.01

					diags_opt <- list(rep(1, times = n), rep(-rho_opt, times = n-1))
					H_r_t_opt <-  bandSparse(n, k = 0:1, diagonals = diags_opt, symmetric = FALSE)
					H_r_opt <- t(H_r_t_opt) 

					Sigma_opt <- solve(Delta_t %*% H_r_t_opt %*% H_r_opt %*% Delta)

					# Simplification and Cholesky factorization of the Sigma_opt 

					Uchol_opt <- chol(C %*% Sigma_opt %*% t(C))
					Lchol_opt <- t(Uchol_opt)

					# Preconditioning the variables

					X_F_opt <- solve(Lchol_opt) %*% C %*% X
					Y_F_opt <- solve(Lchol_opt) %*% Y  

					# Estimate betaHat_1 using GLS assuming optimal Sigma

					betaHat_opt <- solve(t(X_F_opt) %*% X_F_opt) %*% t(X_F_opt) %*% Y_F_opt 

					# The distribution matrix

					D <- Sigma %*% t(C) %*% solve(C %*% Sigma %*% t(C))

					# Obtain the residuals using betaHat_1

					u_l <- Y - C %*% X %*% betaHat_opt  

					# Generate the high-frequency series

					y <- X %*% betaHat_opt + (D %*% u_l)

				}


			}


		}


	}

	data_list <- list(y, betaHat_opt, rho_opt, u_l)
	names(data_list) <- c("y_Gen", "betaHat_Opt", "rho_Opt","u_Gen")
	
	return(data_list)
}
