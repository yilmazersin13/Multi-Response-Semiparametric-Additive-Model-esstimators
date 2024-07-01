# Load required libraries
library(svMisc)
library(pspline)
library(splines)
library(psych)
library(splines2)
library(Rlab)
library(pracma)
library(matrixcalc)
library(mvtnorm)
library(rockchalk)
library(optR)
library(condSURV)

# Data Generation Function
datagen <- function(n, CL, cov) {
  # n: sample size
  # CL: Censoring level
  # no.z: number of responses (3 or 5)
  #-----------------------------------------------------------------------------
  x <- matrix(c(rnorm(n / 3), rnorm(n / 3)), n / 3, 2)
  beta1 <- c(1, 2)
  beta2 <- c(0.5, -2)
  beta3 <- c(-1, 0.2)
  m <- 3
  #------------------------------------------------------------------------------
  # Creating Block Design Matrix
  zeros <- matrix(0, n / 3, ncol(x))
  X_B1 <- matrix(c(x, zeros, zeros), n / 3, 6)
  X_B2 <- matrix(c(zeros, x, zeros), n / 3, 6)
  X_B3 <- matrix(c(zeros, zeros, x), n / 3, 6)
  X_mat <- rbind(X_B1, X_B2, X_B3)
  #------------------------------------------------------------------------------ 
  t1 <- 0
  for (j in 1:(n)) {
    t1[j] <- 2.5 * (j - 0.5) / (n)
  }
  #----------------------------------------------------------------------------
  g <- t1^2 * (sin(0.95 * -t1^2))
  #------------------------------------------------------------------------------ 
  mu <- c(0, 0, 0)  # Mean
  Sigma <- matrix(c(1, cov, cov, cov, 1, cov, cov, cov, 1), 3)  # Covariance matrix
  is.positive.definite(Sigma, tol = 1e-100)
  error <- mvrnorm(n / 3, mu = mu, Sigma = Sigma)
  ERR <- matrix(c(error[, 1], error[, 2], error[, 3]), n, 1)
  p <- 2
  BETA <- matrix(c(beta1, beta2, beta3), p * 3, 1)
  Z <- X_mat %*% BETA + g + ERR
  #------------------------------------------------------------------------------  
  # Generation of the right-censored data
  censoring <- function(n, CL, z) {
    delta <- 1 - rbern(n, CL)
    c <- matrix(0, n, 1)
    y <- 0
    for (i in 1:n) {
      if (delta[i] == 0) {
        while (Z[i] <= c[i]) {
          c[i] <- rnorm(1, mean(z), sd = sd(z))
        }
      } else {
        c[i] <- z[i]
      }
    }
    for (j in 1:n) {
      if (z[j] <= c[j]) {
        y[j] <- z[j]
      } else {
        y[j] <- c[j]
      }
    }
    return(y)
  }
  
  syndata <- function(y, delta) {
    library(pracma)
    # where y: right-censored observations, delta: censorship indicator
    n <- length(y)
    M <- 0
    yg <- 0
    M <- 0
    delta2 <- 0
    # Synthetic data transformation
    y1 <- cbind(y, delta)
    y2 <- sortrows(y1, k = 2)
    delta2 <- y2[, 2]
    delta2[n] <- 1
    sy <- sort(y)
    for (i1 in 1:n) {
      Ma <- 1
      for (i2 in 1:n) {
        mGt <- ((n - i2) / (n - i2 + 1))^((delta2[i2] == 0) * (sy[i2] <= sy[i1]))
        Ma <- Ma * mGt
      }
      M[i1] = Ma
      yg[i1] <- (delta[i1] * y[i1]) / M[i1]
    }
    return(yg)
  }
  delta <- 1 - rbern(n, CL)
  
  Y <- censoring(n, CL, Z)
  syn.Y <- syndata(Z, delta)
  
  #------------------------------------------------------------------------------
  dat <- new.env()
  dat$X <- X_mat
  dat$otherx <- x
  dat$Z <- Z
  dat$Y <- Y
  dat$syn.Y <- syn.Y
  dat$nc <- t1
  dat$BETA <- BETA
  dat$g <- g
  return(dat)
}

# GCV Function
gcvfunc <- function(y, yhat, p) {
  y <- matrix(y)
  yhat <- matrix(yhat)
  n <- length(y)
  score <- (1 / n) * (norm(y - yhat)^2) / ((1 / n) * p)
  return(score)
}

# Local Polynomial Estimation Function
local_poly_estimation <- function(data, sn, bw_seq, LO) {
  x <- data$X
  y <- data$syn.Y
  t <- data$nc
  n <- length(t)
  
  ones <- matrix(1, n, 1)
  index1 <- seq(min(t) - 0.1, max(t) + 0.1, length.out = n)
  T1 <- matrix(c(ones, (t - index1), (t - index1)^2), n, 3)
  GCV_g <- numeric(LO)
  
  for (i in 1:LO) {
    W1 <- matrix(0, n, n)
    for (j in 1:n) {
      W1[, j] <- NWW(t, index1[j], bw = bw_seq[i])
    }
    S_g1 <- T1 %*% solve(t(T1) %*% W1 %*% T1) %*% t(T1) %*% W1
    xtil <- (diag(n) - S_g1) %*% x
    ytil <- (diag(n) - S_g1) %*% y
    var.m <- var(ytil)
    cov_mat <- var.m[1] * t(x) %*% W1 %*% x
    BETAHAT <- solve(t(xtil) %*% xtil %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(xtil) %*% ytil
    g1hat <- W1 %*% (y - x %*% BETAHAT)
    pg1 <- tr(S_g1)
    GCV_g[i] <- gcvfunc(data$g, g1hat[, 1], pg1)
  }
  
  bwg <- bw_seq[which.min(GCV_g)]
  
  # Use optimal bandwidth bwg to get final estimates
  W1 <- matrix(0, n, n)
  for (j in 1:n) {
    W1[, j] <- NWW(t, index1[j], bw = bwg)
  }
  S_g1 <- T1 %*% solve(t(T1) %*% W1 %*% T1) %*% t(T1) %*% W1
  xtil <- (diag(n) - S_g1) %*% x
  ytil <- (diag(n) - S_g1) %*% y
  var.m <- var(ytil)
  cov_mat <- var.m[1] * t(x) %*% W1 %*% x
  BETAHAT <- solve(t(xtil) %*% xtil %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(xtil) %*% ytil
  g1hat <- W1 %*% (y - x %*% BETAHAT)
  
  results <- list(BETAHAT = BETAHAT, g1hat = g1hat, GCV_g = GCV_g, bw_seq = bw_seq, bwg = bwg)
  return(results)
}

# Local Constant Estimation Function
local_constant_estimation <- function(data, sn, bw_seq, LO) {
  x <- data$X
  y <- data$syn.Y
  t <- data$nc
  n <- length(t)
  
  GCV_g <- numeric(LO)
  
  for (i in 1:LO) {
    K <- function(u) { return(dnorm(u)) }
    S_hk0 <- matrix(0, n, n)
    for (j in 1:n) {
      weights <- K((t - t[j]) / bw_seq[i])
      S_hk0[j, ] <- weights / sum(weights)
    }
    S_h0 <- diag(rep(1, n)) - S_hk0
    Y_tilde_G <- S_h0 %*% y
    X_tilde <- S_h0 %*% x
    var.m <- var(Y_tilde_G)
    cov_mat <- var.m[1] * t(x) %*% S_hk0 %*% x
    BETAHAT_LC <- solve(t(X_tilde) %*% X_tilde %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(X_tilde) %*% Y_tilde_G
    g_hat_LC <- S_hk0 %*% (y - x %*% BETAHAT_LC)
    pg1 <- tr(S_hk0)
    GCV_g[i] <- gcvfunc(data$g, g_hat_LC[, 1], pg1)
  }
  
  bwg <- bw_seq[which.min(GCV_g)]
  
  # Use optimal bandwidth bwg to get final estimates
  K <- function(u) { return(dnorm(u)) }
  S_hk0 <- matrix(0, n, n)
  for (j in 1:n) {
    weights <- K((t - t[j]) / bwg)
    S_hk0[j, ] <- weights / sum(weights)
  }
  S_h0 <- diag(rep(1, n)) - S_hk0
  Y_tilde_G <- S_h0 %*% y
  X_tilde <- S_h0 %*% x
  var.m <- var(Y_tilde_G)
  cov_mat <- var.m[1] * t(x) %*% S_hk0 %*% x
  BETAHAT_LC <- solve(t(X_tilde) %*% X_tilde %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(X_tilde) %*% Y_tilde_G
  g_hat_LC <- S_hk0 %*% (y - x %*% BETAHAT_LC)
  
  results <- list(BETAHAT = BETAHAT_LC, g1hat = g_hat_LC, GCV_g = GCV_g, bw_seq = bw_seq, bwg = bwg)
  return(results)
}

# Local Linear Estimation Function
local_linear_estimation <- function(data, sn, bw_seq, LO) {
  x <- data$X
  y <- data$syn.Y
  t <- data$nc
  n <- length(t)
  
  ones <- matrix(1, n, 1)
  index1 <- seq(min(t) - 0.1, max(t) + 0.1, length.out = n)
  T1 <- matrix(c(ones, (t - index1)), n, 2)
  GCV_g <- numeric(LO)
  
  for (i in 1:LO) {
    W1 <- matrix(0, n, n)
    for (j in 1:n) {
      W1[, j] <- NWW(t, index1[j], bw = bw_seq[i])
    }
    S_g1 <- T1 %*% solve(t(T1) %*% W1 %*% T1) %*% t(T1) %*% W1
    xtil <- (diag(n) - S_g1) %*% x
    ytil <- (diag(n) - S_g1) %*% y
    var.m <- var(ytil)
    cov_mat <- var.m[1] * t(x) %*% W1 %*% x
    BETAHAT <- solve(t(xtil) %*% xtil %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(xtil) %*% ytil
    g1hat <- W1 %*% (y - x %*% BETAHAT)
    pg1 <- tr(S_g1)
    GCV_g[i] <- gcvfunc(data$g, g1hat[, 1], pg1)
  }
  
  bwg <- bw_seq[which.min(GCV_g)]
  
  # Use optimal bandwidth bwg to get final estimates
  W1 <- matrix(0, n, n)
  for (j in 1:n) {
    W1[, j] <- NWW(t, index1[j], bw = bwg)
  }
  S_g1 <- T1 %*% solve(t(T1) %*% W1 %*% T1) %*% t(T1) %*% W1
  xtil <- (diag(n) - S_g1) %*% x
  ytil <- (diag(n) - S_g1) %*% y
  var.m <- var(ytil)
  cov_mat <- var.m[1] * t(x) %*% W1 %*% x
  BETAHAT <- solve(t(xtil) %*% xtil %*% solve(cov_mat)) %*% solve(cov_mat) %*% t(xtil) %*% ytil
  g1hat <- W1 %*% (y - x %*% BETAHAT)
  
  results <- list(BETAHAT = BETAHAT, g1hat = g1hat, GCV_g = GCV_g, bw_seq = bw_seq, bwg = bwg)
  return(results)
}

# Function to perform one simulation
perform_simulation <- function(sn, CL, cov, bw_seq, LO) {
  data <- datagen(sn * 3, CL, cov)
  
  # Local Polynomial Estimation
  poly_results <- local_poly_estimation(data, sn, bw_seq, LO)
  beta_hat_poly <- poly_results$BETAHAT
  g_hat_poly <- poly_results$g1hat
  
  # Local Constant Estimation
  const_results <- local_constant_estimation(data, sn, bw_seq, LO)
  beta_hat_const <- const_results$BETAHAT
  g_hat_const <- const_results$g1hat
  
  # Local Linear Estimation
  linear_results <- local_linear_estimation(data, sn, bw_seq, LO)
  beta_hat_linear <- linear_results$BETAHAT
  g_hat_linear <- linear_results$g1hat
  
  # Compute Biases
  bias_poly <- data$BETA - beta_hat_poly
  bias_const <- data$BETA - beta_hat_const
  bias_linear <- data$BETA - beta_hat_linear
  
  list(
    beta_poly = beta_hat_poly,
    beta_const = beta_hat_const,
    beta_linear = beta_hat_linear,
    bias_poly = bias_poly,
    bias_const = bias_const,
    bias_linear = bias_linear,
    g_hat_poly = g_hat_poly,
    g_hat_const = g_hat_const,
    g_hat_linear = g_hat_linear
  )
}

# Function to perform multiple simulations for different configurations
perform_all_simulations <- function() {
  sample_sizes <- c(50, 100, 200)
  censoring_levels <- c(0.05, 0.15)
  covariance_levels <- c(0, 0.9)
  num_simulations <- 100
  LO <- 5
  bw_seq <- seq(0.035, 0.3, length.out = LO)
  
  overall_summary_table <- data.frame()
  beta_summary_table <- data.frame()
  mse_summary_table <- data.frame()
  ss <- 1
  for (sn in sample_sizes) {
    for (CL in censoring_levels) {
      for (cov in covariance_levels) {
        set.seed(123 + ss)  # For reproducibility
        ss <- ss + 1
        results <- replicate(num_simulations, perform_simulation(sn, CL, cov, bw_seq, LO), simplify = FALSE)
        data <- datagen(sn * 3, CL, cov)
        
        # Extract and analyze results
        beta_poly_all <- do.call(rbind, lapply(results, function(res) res$beta_poly))
        beta_const_all <- do.call(rbind, lapply(results, function(res) res$beta_const))
        beta_linear_all <- do.call(rbind, lapply(results, function(res) res$beta_linear))
        
        bias_poly_all <- do.call(rbind, lapply(results, function(res) res$bias_poly))
        bias_const_all <- do.call(rbind, lapply(results, function(res) res$bias_const))
        bias_linear_all <- do.call(rbind, lapply(results, function(res) res$bias_linear))
        
        g_hat_poly_all <- do.call(cbind, lapply(results, function(res) res$g_hat_poly))
        g_hat_const_all <- do.call(cbind, lapply(results, function(res) res$g_hat_const))
        g_hat_linear_all <- do.call(cbind, lapply(results, function(res) res$g_hat_linear))
        
        # Calculate variances of betas
        var_beta_poly <- apply(beta_poly_all, 2, var)
        var_beta_const <- apply(beta_const_all, 2, var)
        var_beta_linear <- apply(beta_linear_all, 2, var)
        
        # Calculate Mean Squared Errors (MSE) for nonparametric parts
        mse_g_poly <- mean((rowMeans(g_hat_poly_all) - data$g)^2)
        mse_g_const <- mean((rowMeans(g_hat_const_all) - data$g)^2)
        mse_g_linear <- mean((rowMeans(g_hat_linear_all) - data$g)^2)
        
        # Create summary table for current configuration
        current_overall_summary <- data.frame(
          Sample_Size = sn,
          Censoring_Level = CL,
          Covariance_Level = cov,
          Estimator = c("LP", "LC", "LL"),
          MSE_Nonparametric = c(mse_g_poly, mse_g_const, mse_g_linear),
          Bias_Avg = c(mean(abs(bias_poly_all)), mean(abs(bias_const_all)), mean(abs(bias_linear_all))),
          Variance_Betas = c(mean(var_beta_poly), mean(var_beta_const), mean(var_beta_linear))
        )
        
        # Append to the overall summary table
        overall_summary_table <- rbind(overall_summary_table, current_overall_summary)
        
        # Bias and variance table for betas
        current_beta_summary <- data.frame(
          Sample_Size = sn,
          Censoring_Level = CL,
          Covariance_Level = cov,
          Estimator = rep(c("LP", "LC", "LL"), each = 6),
          Beta = rep(1:6, 3),
          Bias = c(colMeans(matrix(bias_poly_all, nrow = num_simulations, byrow = TRUE)),
                   colMeans(matrix(bias_const_all, nrow = num_simulations, byrow = TRUE)),
                   colMeans(matrix(bias_linear_all, nrow = num_simulations, byrow = TRUE))),
          Variance = c(apply(matrix(beta_poly_all, nrow = num_simulations, byrow = TRUE), 2, var),
                       apply(matrix(beta_const_all, nrow = num_simulations, byrow = TRUE), 2, var),
                       apply(matrix(beta_linear_all, nrow = num_simulations, byrow = TRUE), 2, var))
        )
        
        # Append to the beta summary table
        beta_summary_table <- rbind(beta_summary_table, current_beta_summary)
        
        # MSE table for different parts of g functions
        mse_g1_poly <- mean((rowMeans(g_hat_poly_all[1:sn, ]) - data$g[1:sn])^2)
        mse_g2_poly <- mean((rowMeans(g_hat_poly_all[(sn + 1):(2 * sn), ]) - data$g[(sn + 1):(2 * sn)])^2)
        mse_g3_poly <- mean((rowMeans(g_hat_poly_all[(2 * sn + 1):(3 * sn), ]) - data$g[(2 * sn + 1):(3 * sn)])^2)
        
        mse_g1_const <- mean((rowMeans(g_hat_const_all[1:sn, ]) - data$g[1:sn])^2)
        mse_g2_const <- mean((rowMeans(g_hat_const_all[(sn + 1):(2 * sn), ]) - data$g[(sn + 1):(2 * sn)])^2)
        mse_g3_const <- mean((rowMeans(g_hat_const_all[(2 * sn + 1):(3 * sn), ]) - data$g[(2 * sn + 1):(3 * sn)])^2)
        
        mse_g1_linear <- mean((rowMeans(g_hat_linear_all[1:sn, ]) - data$g[1:sn])^2)
        mse_g2_linear <- mean((rowMeans(g_hat_linear_all[(sn + 1):(2 * sn), ]) - data$g[(sn + 1):(2 * sn)])^2)
        mse_g3_linear <- mean((rowMeans(g_hat_linear_all[(2 * sn + 1):(3 * sn), ]) - data$g[(2 * sn + 1):(3 * sn)])^2)
        
        current_mse_summary <- data.frame(
          Sample_Size = sn,
          Censoring_Level = CL,
          Covariance_Level = cov,
          Estimator = rep(c("LP", "LC", "LL"), each = 3),
          g_Part = rep(c("g1", "g2", "g3"), 3),
          MSE = c(mse_g1_poly, mse_g2_poly, mse_g3_poly,
                  mse_g1_const, mse_g2_const, mse_g3_const,
                  mse_g1_linear, mse_g2_linear, mse_g3_linear)
        )
        
        # Append to the MSE summary table
        mse_summary_table <- rbind(mse_summary_table, current_mse_summary)
        
        # Save plots
        plot_dir <- paste0("plots/sn_", sn, "_CL_", CL, "_cov_", cov)
        dir.create(plot_dir, showWarnings = TRUE, recursive = TRUE)
        
        # Plot boxplots of biases for each beta
        jpeg(file = paste0(plot_dir, "/boxplots_biases.jpg"), width = 600, height = 350, quality = 100)
        par(mfrow = c(2, 3))
        for (i in 1:6) {
          start_idx <- (i - 1) * num_simulations + 1
          end_idx <- i * num_simulations
          boxplot(
            bias_poly_all[start_idx:end_idx],
            bias_const_all[start_idx:end_idx],
            bias_linear_all[start_idx:end_idx],
            names = c("LP", "LC", "LL"),
            main = paste("Bias of Beta", i),
            ylab = "Bias", col = c("red", "green", "blue"),cex=1.3,
            cex.axis = 1.5,
            cex.lab = 1.3,   # Increase size of axis labels
            cex.main = 1.3   # Increase size of main title
          )
          grid()
        }
        dev.off()
        
        # Plot nonparametric function estimations
        plot_nonparametric <- function(g_hat_all, true_g, title) {
          matplot(data$nc, g_hat_all, type = "l", col = rgb(0.7, 0.7, 0.7, 0.5), lty = 1,
                  main = title, xlab = "t", ylab = "g(t)", ylim = c(-4, 8))
          lines(data$nc, rowMeans(g_hat_all), col = "blue", lwd = 2, ylim = c(-4, 8))
          lines(data$nc, true_g, col = "red", lwd = 2, ylim = c(-4, 8))
          abline(v = data$nc[sn], lty = 2, col = "gray", lwd = 2)
          abline(v = data$nc[2 * sn], lty = 2, col = "gray", lwd = 2)
          text(0.4, 2.5, "g(y1)", cex = 1.6)
          text(1.3, 2.5, "g(y2)", cex = 1.6)
          text(2.3, 0.5, "g(y3)", cex = 1.6)
          legend("topright", legend = c("Fitted curve", "True g(t)"),
                 col = c("blue", "red"), lwd = 2,cex=1.5)
        }
        
        jpeg(file = paste0(plot_dir, "/nonparametric_estimations.jpg"), width = 1000, height = 600, quality = 100)
        par(mfrow = c(1, 3))
        plot_nonparametric(g_hat_poly_all, data$g, "LP Fit")
        plot_nonparametric(g_hat_const_all, data$g, "LC Fit")
        plot_nonparametric(g_hat_linear_all, data$g, "LL Fit")
        dev.off()
        message(ss, "ends")
      }
    }
  }
  
  # Save the summary tables to CSV files
  write.csv(overall_summary_table, file = "overall_summary_table.csv", row.names = FALSE)
  write.csv(beta_summary_table, file = "beta_summary_table.csv", row.names = FALSE)
  write.csv(mse_summary_table, file = "mse_summary_table.csv", row.names = FALSE)
}

# Run all simulations
perform_all_simulations()
