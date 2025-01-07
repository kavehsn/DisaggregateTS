## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(DisaggregateTS)

## -----------------------------------------------------------------------------
# Load the combined data from the package
data(Data)

# Extract Data_Y and Data_X from the combined data
Data_Y <- Data$Data_Y
Data_X <- Data$Data_X

# Select IBM GHG data and dates for Q3 2005 - Q3 2021
Dates <- Data_Y$Dates[c(7:23)]
Y <- Data_Y$IBM[c(7:23)]
Y <- as.matrix(as.numeric(Y))

# HF data available from 12-2004 (observation 21) up to 09-2021 (observation 88)
Dates_Q <- Data_X$Dates[c(21:88)]
X <- Data_X[c(21:88),]
X <- sapply(X, as.numeric)

# Remove columns containing NAs
X <- X[ , colSums(is.na(X))==0] 


# Remove highly correlated variables (pairwise correlation >= 0.99)
tmp <- cor(X)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

X2 <- X[, !apply(tmp, 2, function(x) any(abs(x) >= 0.99, na.rm = TRUE))]

## -----------------------------------------------------------------------------
C_sparse <- disaggregate(
  as.matrix(Y),
  as.matrix(X2),
  aggMat   = "sum",
  aggRatio = 4,
  method   = "adaptive-spTD")

# Temporally disaggregated time series
Y_HF <- C_sparse$y_Est

## ----plot-results, fig.width=8, fig.height=5, echo=TRUE-----------------------
par(mar = c(5, 6, 4, 5) + 0.1)  # Adjust margins for better spacing

# Plot the temporal disaggregated data
plot(Dates_Q, Y_HF, type = "b", pch = 19, ylab = "GHG emissions", xlab = "Time",
     lwd = 2, cex.lab = 1.4, cex.axis = 1.2, main = "Temporal Disaggregation of GHG Emissions")

# Add a legend with adjusted font size and position
legend("bottomleft", inset = 0.05, 
       legend = "Temporal disaggregated observations",
       col = "black", lty = 1, lwd = 2, pch = 19, 
       cex = 1.2, pt.cex = 1.2)

