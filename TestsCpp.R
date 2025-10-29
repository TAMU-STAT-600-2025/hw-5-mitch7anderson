
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
testthat::expect_equal(soft_c(2, 1), soft(2, 1))
testthat::expect_equal(soft_c(5, 7), soft(5, 7))

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
# Set up test data
set.seed(123)
X <- matrix(rnorm(20), 5, 4)
Y <- rnorm(5)
beta <- rnorm(4)
lambda <- 0.1

# Compute R and cpp functions
out_R <- lasso(X, Y, beta, lambda)
out_cpp <- lasso(X, Y, beta, lambda)

# Test equality
testthat::expect_equal(out_R, out_cpp)

# Different lambda
lambda <- 2.2

# Compute R and cpp functions
out_R <- lasso(X, Y, beta, lambda)
out_cpp <- lasso(X, Y, beta, lambda)

# Test equality
testthat::expect_equal(out_R, out_cpp)

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
