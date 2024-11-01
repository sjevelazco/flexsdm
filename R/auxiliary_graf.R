## %######################################################%##
#                                                          #
####                     GRaF codes                     ####
#                                                          #
## %######################################################%##
# This codes belogn to GRaF package available at https://github.com/goldingn/GRaF

## citation
# Golding, N., & Purse, B. V. (2016). Fast and flexible Bayesian species distribution modelling
# using Gaussian processes. Methods in Ecology and Evolution, 7, 598–608.
# https://doi.org/10.1111/2041-210X.12523


# density and gradient of the prior over the log hyperparameters
theta.prior <- function(theta, pars) {
  if (any(is.na(pars))) {
    return(list(
      density = 0,
      gradient = rep(0, length(theta))
    ))
  }

  density <- sum(stats::dnorm(theta, pars[1], pars[2], log = TRUE))

  gradient <- (pars[1] - theta) / pars[2]^2

  return(list(
    density = density,
    gradient = gradient
  ))
}

# psiline
psiline <-
  function(s, adiff, a, K, y, d0, mn = 0, wt) {
    a <- a + s * as.vector(adiff)
    f <- K %*% a + mn
    psi(a, f, mn, y, d0, wt)
  }

# psi
psi <-
  function(a, f, mn, y, d0, wt) {
    0.5 * t(a) %*% (f - mn) - sum(wt * d0(f, y))
  }


# moments of the standard normal distribution for EP on the probit GP
ep.moments <- function(y, sigma2, mu) {
  # 1 + sigma^{2}
  sigma2p1 <- 1 + sigma2

  # get the likelihood
  z <- y * (mu / sqrt(sigma2p1))

  # get cdf and pdf of the likelihood
  cdf_z <- stats::pnorm(z)
  pdf_z <- stats::dnorm(z)

  # new \hat{\mu} (mean; second moment)
  muhat <- mu + (y * sigma2 * pdf_z) / cdf_z

  # new \hat{\sigma^{2}} (variance; third moment)
  sigma2hat <- sigma2 -
    (sigma2^2 * pdf_z) /
      (sigma2p1 * cdf_z) *
      (z + pdf_z / cdf_z)

  # log of the 0th moment (cdf)
  logM0 <- log(cdf_z)

  # return the three elements
  return(list(
    logM0 = logM0,
    muhat = muhat,
    sigma2hat = sigma2hat
  ))
}

# update.ep
update.ep <-
  function(i, y, mn, lis) {
    # update at index \code{i} for the EP approximation given \code{y}
    # (on +1, -1 scale), \code{mn} (on the Gaussian scale) and the current
    # list of parameters \code{lis}.

    Sigma <- lis$Sigma
    ttau <- lis$ttau
    tnu <- lis$tnu
    mu <- lis$mu

    # calculate approximate cavity parameters \nu_{-i} and \tau_{-i}
    # \sigma^{2}_i = \Sigma_{i, i}
    # \tau_{-i} = \sigma^{-2}_i - \tilde{\tau}_i
    # \nu_{-i} = \mu_i / \sigma^{2}_i + m_i * \tau_{-i} - \tilde{\nu}_i
    sig2_i <- Sigma[i, i]
    tau_ni <- 1 / sig2_i - ttau[i]
    nu_ni <- mu[i] / sig2_i + mn[i] * tau_ni - tnu[i]

    # compute marginal moments \hat{\mu}_i and \hat{\sigma}_i^2
    # calculate required derivatives of individual log partition function
    z <- nu_ni / tau_ni / (sqrt(1 + 1 / tau_ni))
    yz <- y[i] * z
    lZ <- stats::pnorm(yz, log.p = TRUE)
    n_p <- stats::dnorm(yz) / exp(lZ)
    dlZ <- y[i] * n_p / sqrt(1 + 1 / tau_ni)
    d2lZ <- -n_p * (yz + n_p) / (1 + 1 / tau_ni)

    # save old \tilde{\tau} before finding new one
    ttau_old <- ttau[i]

    # update \tilde{\tau}_i, forcing it non-negative
    # then update \tilde{\nu}_i
    ttau[i] <- max(-d2lZ / (1 + d2lZ / tau_ni), 0)
    tnu[i] <- (dlZ + (mn[i] - nu_ni / tau_ni) * d2lZ) / (1 + d2lZ / tau_ni)

    # rank-1 update \Sigma and \mathbf{\mu}
    # \delta\sigma^2
    ds2 <- ttau[i] - ttau_old
    # get column i
    #   si <- Sigma[, i]
    # recompute\Sigma \& \mu
    lis <- update.sigma.mu(Sigma, ds2, i, tnu)

    # return update list of parameters
    return(list(
      Sigma = lis$Sigma,
      ttau = ttau,
      tnu = tnu,
      mu = lis$mu
    ))
  }

# update.sigma.mu
update.sigma.mu <-
  function(Sigma, ds2, i, tnu) {
    # this takes 70% of total time in gpml Matlab code
    # tried re-coding in C++, but all the cost is the matrx algebra,
    # so little improvement. The approach is inherently slow.

    # get column i
    si <- Sigma[, i]
    # recalculate Sigma
    Sigma <- Sigma - ds2 / (1 + ds2 * si[i]) * si %*% t(si)
    # and mu
    mu <- Sigma %*% tnu

    return(list(
      Sigma = Sigma,
      mu = mu
    ))
  }

# cov.SE
cov.SE <- function(x1, x2 = NULL, e1 = NULL, e2 = NULL, l) {
  n1 <- nrow(x1)
  n2 <- ifelse(is.null(x2), n1, nrow(x2))
  n3 <- ncol(x1)

  # distance matrices
  if (is.null(x2)) {
    e2 <- e1
    # if no second matrix do with distance matrices for speed up
    dists <- lapply(1:n3, function(i, x) dist(x[, i])^2, x1)
  } else {
    dists <- list()
    for (i in 1:n3) {
      dists[[i]] <- x1[, i]^2 %*% t(rep(1, n2)) +
        rep(1, n1) %*% t(x2[, i]^2) - 2 * x1[, i] %*% t(x2[, i])
    }
  }

  # with error matrices
  if (!is.null(e1)) {
    E1 <- list()
    ones <- t(rep(1, n2))
    for (i in 1:n3) {
      E1[[i]] <- e1[, i] %*% ones
    }

    if (!is.null(e2)) {
      E2 <- list()
      ones <- t(rep(1, n1))
      for (i in 1:n3) {
        E2[[i]] <- t(e2[, i] %*% ones)
      }
    } else {
      E2 <- as.list(rep(0, n3))
    }

    # run through each covariate

    sumdiffs <- 0
    denom <- 1
    lower <- lower.tri(E1[[1]])
    for (i in 1:n3) {
      err <- E1[[i]] + E2[[i]]
      if (is.null(x2)) {
        err <- err[lower] # save only lower portion for speed up
      }
      sumdiffs <- sumdiffs + dists[[i]] / (err + l[i])
      denom <- denom * (1 + err / l[i])
    }
    # inverse kronecker delta
    ikds <- as.numeric(sumdiffs > 0)
    diag(ikds <- 1)
    denom <- sqrt(denom) * ikds
    K <- exp(-0.5 * sumdiffs) / denom
  } else {
    # without error matrices
    sumdiffs <- 0
    for (i in 1:n3) {
      sumdiffs <- sumdiffs + dists[[i]] / l[i]
    }
    K <- exp(-0.5 * sumdiffs) # to matrix?
  }

  if (class(sumdiffs)[1] == "dist") {
    K <- as.matrix(K)
    diag(K) <- 1
  }
  K
}

# cov.SE.d1
cov.SE.d1 <- function(x, e = NULL, l) {
  # get gradients (matrices) of the kernel wrt. the parameters
  # CURRENTLY IGNORES e!!

  # number of parameters
  n <- length(l)

  # assign vector for gradients
  grads <- list()

  # get full covariance matrix
  K <- cov.SE(x1 = x, e1 = e, l = l)

  # loop through them
  for (i in 1:n) {
    # squared distances
    d2_i <- as.matrix(dist(x[, i])^2)

    # gradient for each parameter
    grads[[i]] <- K * (1 / l[i]^2) * d2_i / 2
  }
  # return as a list
  return(grads)
}

# d0
d0 <-
  function(z, y) {
    if (length(y) != length(z)) y <- rep(y, length(z))
    pr <- y > 0 & y < 1
    npr <- !pr
    ans <- vector("numeric", length(y))
    y[pr] <- stats::qnorm(y[pr])
    y[npr] <- 2 * y[npr] - 1
    ans[pr] <- stats::dnorm(y[pr], z[pr], log = TRUE)
    ans[npr] <- stats::pnorm(y[npr] * z[npr], log.p = TRUE)
    ans
  }

# d1
d1 <-
  function(z, y) {
    pr <- y > 0 & y < 1
    npr <- !pr
    ans <- vector("numeric", length(y))
    y[pr] <- stats::qnorm(y[pr])
    y[npr] <- 2 * y[npr] - 1
    ans[pr] <- y[pr] - z[pr]
    ans[npr] <- y[npr] * stats::dnorm(z[npr]) / stats::pnorm(y[npr] * z[npr])
    ans
  }

# d2
d2 <-
  function(z, y) {
    pr <- y > 0 & y < 1
    npr <- !pr
    ans <- vector("numeric", length(y))
    y[npr] <- 2 * y[npr] - 1
    ans[pr] <- -1
    a <- stats::dnorm(z[npr])^2 / stats::pnorm(y[npr] * z[npr])^2
    b <- y[npr] * z[npr] * stats::dnorm(z[npr]) / stats::pnorm(y[npr] * z[npr])
    ans[npr] <- -a - b
    ans
  }

# d3
d3 <- function(z, y) {
  pr <- y > 0 & y < 1
  npr <- !pr
  ans <- vector("numeric", length(y))
  y[npr] <- 2 * y[npr] - 1
  ans[pr] <- 0
  n <- stats::dnorm(z[npr])
  p <- stats::pnorm(y[npr] * z[npr])
  f <- z[npr]
  a <- n / p^3
  b <- f * (y[npr]^2 + 2) * n * p
  c <- y[npr] * (f^2 - 1) * p^2 + 2 * y[npr] * n^2
  ans[npr] <- a * (b + c)
  ans
}


# define the objective and gradient functions to optimise the hyperparameters
# of a GRaF model
objective <- function(theta, prior.pars, isfac, args, fun) {
  # unpack theta
  l <- ifelse(isfac, 0.01, NA)
  l[!isfac] <- exp(theta)

  # set the lengthscales
  args$l <- l

  # run the model
  m <- do.call(fun, args)

  # log likelihood and prior
  llik <- -m$mnll
  lpri <- theta.prior(theta, prior.pars)$density

  # log posterior
  lpost <- llik + lpri

  # and objective
  objective <- -lpost

  return(objective)
}


gradient <- function(theta, prior.pars, isfac, args, fun) {
  # unpack theta
  l <- ifelse(isfac, 0.01, NA)
  l[!isfac] <- exp(theta)

  # set the lengthscales
  args$l <- l

  # run the model
  m <- do.call(fun, args)

  # gradient of llik w.r.t. l
  dLdl <- m$l_grads[!isfac]

  # gradient of l w.r.t. theta
  dldtheta <- exp(theta)

  # gradient of lpri w.r.t. theta
  dpdtheta <- theta.prior(theta, prior.pars)$gradient

  # gradient of llik w.r.t. theta
  dLdtheta <- dLdl * dldtheta

  # gradient of lpost w.r.t. theta
  dPdtheta <- dLdtheta + dpdtheta

  # gradient of objective w.r.t. lpost
  dOdP <- -1

  # gradient of objective w.r.t. theta
  dOdtheta <- dOdP * dPdtheta

  return(dOdtheta)
}

optimise.graf <- function(args) {
  # pass all the arguments of a call to graf, memoize and
  # optimise the model, and return afitted version

  # set opt.l to FALSE
  args[["opt.l"]] <- FALSE

  # memoise graf
  #   mgraf <- memoise(graf)
  mgraf <- graf

  # set up initial lengthscales
  k <- ncol(args$x)
  l <- rep(1, k)

  # find factors and drop them from theta
  notfacs <- 1:k

  facs <- which(unlist(lapply(args$x, is.factor)))

  if (length(facs) > 0) {
    notfacs <- notfacs[-facs]
    l[facs] <- 0.01
  }

  # log them
  theta <- log(l[notfacs])

  # logical vector of factors
  isfac <- 1:k %in% facs

  # optimisation arguments
  if (args$method == "Laplace") {
    meth <- "L-BFGS-B"
    grad <- gradient
  } else {
    meth <- "BFGS"
    grad <- NULL
  }

  low <- -Inf
  up <- Inf

  opt <- stats::optim(theta,
    fn = objective,
    gr = grad,
    prior.pars = args$theta.prior.pars,
    isfac = isfac,
    args = args,
    fun = mgraf,
    hessian = args$hessian,
    lower = low,
    upper = up,
    method = meth,
    control = args$opt.control
  )

  # get the resultant lengthscales
  l[notfacs] <- exp(opt$par)

  args$l <- l

  # fit the final model and return
  m <- do.call(mgraf, args)

  # replace hessian with the hessian matrix or NULL
  if (args$hessian) {
    m$hessian <- opt$hessian
  } else {
    m$hessian <- NULL
  }

  # un-memoize graf
  # forget(mgraf)

  return(m)
}


# graf.fit.laplace
graf.fit.laplace <-
  function(y, x, mn, l, wt, e = NULL, tol = 10^-6, itmax = 50,
           verbose = FALSE) {
    if (is.vector(x)) x <- as.matrix(x)
    mn <- stats::qnorm(mn)
    n <- length(y)

    # create the covariance matrix
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)

    # an identity matrix for the calculations
    eye <- diag(n)

    # initialise
    a <- rep(0, n)
    f <- mn
    obj.old <- Inf
    obj <- -sum(wt * d0(f, y))
    it <- 0

    # start newton iterations
    while (obj.old - obj > tol & it < itmax) {
      it <- it + 1
      obj.old <- obj
      W <- -(wt * d2(f, y))
      rW <- sqrt(W)
      cf <- f - mn
      mat1 <- rW %*% t(rW) * K + eye
      L <- tryCatch(chol(mat1),
        error = function(x) {
          return(NULL)
        }
      )
      b <- W * cf + wt * d1(f, y)
      mat2 <- rW * (K %*% b)
      adiff <- b - rW * backsolve(L, forwardsolve(t(L), mat2)) - a
      dim(adiff) <- NULL

      # find optimum step size using Brent's method
      res <- stats::optimise(psiline, c(0, 2), adiff, a, as.matrix(K), y, d0, mn, wt)
      a <- a + res$minimum * adiff
      f <- K %*% a + mn
      obj <- psi(a, f, mn, y, d0, wt)
    }

    # recompute key components
    cf <- f - mn
    W <- -(wt * d2(f, y))
    rW <- sqrt(W)
    mat1 <- rW %*% t(rW) * K + eye
    L <- tryCatch(chol(mat1),
      error = function(x) {
        return(NULL)
      }
    )

    # return marginal negative log-likelihood
    mnll <- (a %*% cf)[1, 1] / 2 + sum(log(diag(L)) - (wt * d0(f, y)))

    # get partial gradients of the objective wrt l

    # gradient components
    W12 <- matrix(rep(rW, n), n)
    R <- W12 * backsolve(L, forwardsolve(t(L), diag(rW)))
    C <- forwardsolve(t(L), (W12 * K))

    # partial gradients of the kernel
    dK <- cov.SE.d1(x, e, l)

    # rate of change of likelihood w.r.t. the mode
    s2 <- (diag(K) - colSums(C^2)) / 2 * d3(f, y)

    # vector to store gradients
    l_grads <- rep(NA, length(l))

    for (i in 1:length(l)) {
      grad <- sum(R * dK[[i]]) / 2
      grad <- grad - (t(a) %*% dK[[i]] %*% a) / 2
      b <- dK[[i]] %*% d1(f, y)
      grad <- grad - t(s2) %*% (b - K %*% (R %*% b))
      l_grads[i] <- -grad
    }

    if (verbose) cat(paste("  ", it, "Laplace iterations\n"))
    if (it == itmax) print("timed out, don't trust the inference!")
    return(list(
      y = y, x = x, MAP = f, ls = l, a = a, W = W, L = L, K = K,
      e = e, obsx = x, obsy = y, mnll = mnll, wt = wt,
      l_grads = l_grads
    ))
  }

# graf.fit.ep
graf.fit.ep <-
  function(y, x, mn, l, wt, e = NULL, tol = 1e-6, itmax = 50, itmin = 2,
           verbose = FALSE) {
    # fit a GRaF model using expectation-propagation
    # as implemented in the gpml matlab library
    # If parallel  = TRUE, the EP uses the parallel update
    # as implemented in GPStuff

    # whether to use parallel EP (rather than sequential)
    parallel <- FALSE

    if (is.vector(x)) {
      x <- as.matrix(x)
    }

    # mn to probability scale
    mn <- stats::qnorm(mn)
    n <- length(y)
    # covariance matrix
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    # identity matrix
    eye <- diag(n)

    # convert observations to +1, -1, save 0, 1 version
    oldy <- y
    y <- y * 2 - 1

    # initialise

    # \tilde{\mathbf{\nu}} = \mathbf{0}
    tnu <- rep(0, n)
    # \tilde{\mathbf{\tau}} = \mathbf{0}
    ttau <- rep(0, n)
    # \mathbf{\mu} = \mathbf{0}
    mu <- rep(0, n)
    # \Sigma = \mathbf{K}  (only used in sequential EP)
    Sigma <- K

    # calculate marginal negative log likelihood at ttau = tnu = mu = 0s
    z <- mu / sqrt(1 + diag(K))
    mnll <- -sum(stats::pnorm(y * z), log.p = TRUE)

    # ~~~~~~~~~~~~~~~~~~~~~
    # set up for EP algorithm
    # set up damping factor (for parallel EP) as in GPStuff
    df <- 0.8

    # set old mnll to Inf & start iteration counter
    mnll_old <- Inf
    logM0_old <- logM0 <- rep(0, n)
    it <- 1
    converged <- FALSE

    while (!converged) {
      if (parallel) {
        # calculate key parameters
        dSigma <- diag(Sigma)
        tau <- 1 / dSigma - ttau
        nu <- 1 / dSigma * mn - tnu
        mu <- nu / tau
        sigma2 <- 1 / tau

        # get marginal moments of posterior
        lis <- ep.moments(y, sigma2, mu)

        # recalculate parameters from these moments
        # \delta\tilde{\tau} (change in \tilde{\tau})
        dttau <- 1 / lis$sigma2hat - tau - ttau

        # \tilde{\tau}
        ttau <- ttau + df * dttau

        # \delta\tilde{\nu} (change in \tilde{\nu})
        dtnu <- 1 / lis$sigma2hat * lis$muhat - nu - tnu

        # \tilde{\nu}
        tnu <- tnu + df * dtnu
      } else {
        # otherwise use sequential EP

        lis <- list(
          Sigma = Sigma,
          ttau = ttau,
          tnu = tnu,
          mu = mu
        )

        # cycle through in random order
        for (i in sample(1:n)) {
          lis <- update.ep(i, y, mn, lis)
        } # end permuted for loop

        Sigma <- lis$Sigma
        ttau <- lis$ttau
        tnu <- lis$tnu
        mu <- lis$mu

        sigma2 <- diag(Sigma)
      } # end parallel / sequential

      # recompute the approximate posterior parameters \Sigma and \mathbf{\mu}
      # using eq. 3.53 and eq. 3.68

      # sW = \tilde{S}^{\frac{1}{2}} = \sqrt{\mathbf{\tilde{\tau}}}
      sW <- sqrt(ttau)

      # L = cholesky(I_n  + \tilde{S}^{\frac{1}{2}} * K * \tilde{S}^{\frac{1}{2}})
      L <- chol(sW %*% t(sW) * K + eye)

      # V = L^T \\ \tilde{S}^{\frac{1}{2}} * K
      sWmat <- matrix(rep(sW, n), n) # byrow = TRUE?
      V <- backsolve(L, sWmat * K, transpose = TRUE)

      # \Sigma = \mathbf{K} - \mathbf{V}^T \mathbf{V}
      Sigma <- K - t(V) %*% V

      # \mathbf{\mu} = \Sigma \tilde{\mathbf{\nu}}
      mu <- Sigma %*% tnu # + mn

      # calculate new mnll and assess convergence
      # compute logZ_{EP} using eq. 3.65, 3.73 and 3.74 and the existing L
      # \mathbf{\sigma}^2 = diag(\Sigma)
      sigma2 <- diag(Sigma)
      tau_n <- 1 / sigma2 - ttau
      nu_n <- mu / sigma2 - tnu + mn * tau_n

      z <- nu_n / tau_n / (sqrt(1 + 1 / tau_n))
      lZ <- stats::pnorm(y * z, log.p = TRUE)

      # split the final equation up into 5 bits...
      mnll.a <- sum(log(diag(L))) - sum(lZ) - t(tnu) %*% Sigma %*% tnu / 2
      mnll.b <- t(nu_n - mn * tau_n)
      mnll.c <- ((ttau / tau_n * (nu_n - mn * tau_n) - 2 * tnu) / (ttau + tau_n)) / 2
      mnll.d <- sum(tnu^2 / (tau_n + ttau)) / 2
      mnll.e <- sum(log(1 + ttau / tau_n)) / 2

      mnll <- as.numeric(mnll.a - mnll.b %*% mnll.c + mnll.d - mnll.e)

      # improvement in negative log marginal likelihood
      dmnll <- abs(mnll - mnll_old)

      # improvement in log of the 0th moment
      dlogM0 <- max(abs(logM0 - logM0_old))

      # both under tolerance?
      sub_tol <- dmnll < tol & dlogM0 < tol

      if ((sub_tol & it >= itmin) | it >= itmax) {
        # stop if there was little improvement and there have been at least
        # itmin iterations or if there have been itmax or more
        # iterations
        converged <- TRUE
      } else {
        it <- it + 1
        mnll_old <- mnll
      }
    } # end while loop

    # throw an error if the iterations maxed out before convergence
    if (it >= itmax) {
      stop(paste0("maximum number of iterations (", itmax, ") reached
                  without convergence"))
    }

    # calculate posterior parameters
    sW <- sqrt(ttau)
    alpha <- tnu - sW * backsolve(L, backsolve(L, sW * (K %*% tnu), transpose = TRUE))
    f <- crossprod(K, alpha) + mn

    # return relevant parameters
    return(list(
      y = oldy,
      x = x,
      MAP = f,
      ls = l,
      a = alpha,
      W = ttau,
      L = L,
      K = K,
      e = e,
      obsx = x,
      obsy = oldy,
      mnll = mnll,
      wt = wt,
      l_grads = NULL
    ))
  }

graf <-
  function(y, x, error = NULL, weights = NULL, prior = NULL, l = NULL, opt.l = FALSE,
           theta.prior.pars = c(log(10), 1), hessian = FALSE, opt.control = list(),
           verbose = FALSE, method = c("Laplace", "EP")) {
    method <- match.arg(method)

    # optionally optimise graf (by recursively calling this function)
    if (opt.l) {
      # get all visible object as a list
      args <- capture.all()

      # get the expected objects
      expected_args <- names(formals(graf))

      # remove any unexpected arguments
      args <- args[names(args) %in% expected_args]

      # pass this to optimiser
      fit <- optimise.graf(args)

      # skip out of this function and return
      return(fit)
    }


    if (!is.data.frame(x)) stop("x must be a dataframe")

    # convert any ints to numerics
    for (i in 1:ncol(x)) if (is.integer(x[, i])) x[, i] <- as.numeric(x[, i])

    obsx <- x
    k <- ncol(x)
    n <- length(y)

    if (is.null(weights)) {
      # if weights aren't provided
      weights <- rep(1, n)
    } else {
      # if they are, run some checks
      # throw an error if weights are specified with EP
      if (method == "EP") {
        stop("weights are not implemented for the EP algorithm (yet)")
      }
      # or if any are negative
      if (any(weights < 0)) {
        stop("weights must be positive or zero")
      }
    }

    # find factors and convert them to numerics
    notfacs <- 1:k
    facs <- which(unlist(lapply(x, is.factor)))
    if (length(facs) > 0) notfacs <- notfacs[-facs]
    for (fac in facs) {
      x[, fac] <- as.numeric(x[, fac])
    }
    x <- as.matrix(x)

    # scale the matrix, retaining scaling
    scaling <- apply(as.matrix(x[, notfacs]), 2, function(x) c(mean(x), stats::sd(x)))
    for (i in 1:length(notfacs)) {
      x[, notfacs[i]] <- (x[, notfacs[i]] - scaling[1, i]) / scaling[2, i]
    }

    # set up the default prior, if not specified
    exp.prev <- sum(weights[y == 1]) / sum(weights)
    if (is.null(prior)) {
      mnfun <- function(x) rep(exp.prev, nrow(x))
    } else {
      mnfun <- prior
    }

    # give an approximation to l, if not specified (or optimised)
    if (is.null(l)) {
      l <- rep(0.01, k)
      l[notfacs] <- apply(x[y == 1, notfacs, drop = FALSE], 2, sd) * 8
      l[l == 0] <- 1
    }

    # calculate mean (on unscaled data and probability scale)
    mn <- mnfun(obsx)

    # fit model
    if (method == "Laplace") {
      # by Laplace approximation
      fit <- graf.fit.laplace(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = verbose)
    } else {
      # or using the expectation-propagation algorithm
      fit <- graf.fit.ep(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = FALSE)
    }

    fit$mnfun <- mnfun
    fit$obsx <- obsx
    fit$facs <- facs
    fit$hessian <- hessian
    fit$scaling <- scaling
    fit$peak <- obsx[which(fit$MAP == max(fit$MAP))[1], ]
    class(fit) <- "graf"
    fit
  }

capture.all <- function() {
  # capture all visible objects in the parent environment and pass to a list
  env <- parent.frame()
  object_names <- objects(env)
  objects <- lapply(object_names,
    get,
    envir = env
  )
  names(objects) <- object_names
  return(objects)
}

# pred
pred <-
  function(predx, fit, mn, std = TRUE, maxn = 250, same = FALSE) {
    predx <- as.matrix(predx)
    n <- nrow(predx)
    if (n > maxn & !same) {
      inds <- split(1:n, ceiling((1:n) / maxn))
      fun <- function(ind, X, fit, std, maxn) {
        pred(X[ind, , drop = FALSE], fit, mn[ind], std, maxn, same)
      }
      prediction <- lapply(inds, fun, predx, fit, std, maxn)
      prediction <- do.call("rbind", prediction)
    } else {
      if (same) {
        # if predicting back to input data re-use covariance matrix
        Kx <- fit$K
        prediction <- fit$MAP
      } else {
        Kx <- cov.SE(x1 = fit$x, x2 = predx, e1 = fit$e, e2 = NULL, l = fit$ls)
        mpred <- stats::qnorm(mn)
        prediction <- crossprod(Kx, fit$a) + mpred
      }

      if (std) {
        v <- backsolve(fit$L, sqrt(as.vector(fit$W)) * Kx, transpose = T)
        # using correlation matrix, so diag(kxx) is all 1s, no need to compute kxx
        predvar <- 1 - crossprod(v)
        prediction <- cbind(prediction, sqrt(diag(predvar)))
        colnames(prediction) <- c("MAP", "std")
      }
    }
    prediction
  }

# predict.graf
predict.graf <-
  function(object, newdata = NULL, type = c("response", "latent"),
           CI = 0.95, maxn = NULL, ...) {
    type <- match.arg(type)
    if (is.null(maxn)) maxn <- round(nrow(object$x) / 10)
    # set up data
    if (is.null(newdata)) {
      # use already set up inference data if not specified
      newdata <- object$x
      # get mean on raw data
      mn <- object$mnfun(object$obsx)
    } else {
      # convert any ints to numerics
      for (i in 1:ncol(newdata)) if (is.integer(newdata[, i])) newdata[, i] <- as.numeric(newdata[, i])

      if (is.data.frame(newdata) & all(sapply(object$obsx, class) == sapply(newdata, class))) {
        # get mean on raw data
        mn <- object$mnfun(newdata)

        k <- ncol(newdata)
        # numericize factors
        for (fac in object$facs) {
          newdata[, fac] <- as.numeric(newdata[, fac])
        }
        # convert to a matrix
        newdata <- as.matrix(newdata)
        # scale, if needed
        if (!is.null(object$scaling)) {
          notfacs <- (1:k)
          if (length(object$facs) > 0) notfacs <- notfacs[-object$facs]
          for (i in 1:length(notfacs)) {
            newdata[, notfacs[i]] <- (newdata[, notfacs[i]] - object$scaling[1, i]) / object$   scaling[2, i]
          }
        }
      } else {
        stop("newdata must be either a dataframe with the same elements as used for inference, or NULL")
      }
    }

    # check CI
    if (!is.null(CI)) {
      if (!(CI == "std" & type == "latent")) {
        if (CI >= 1 | CI <= 0) {
          stop("CI must be a number between 0 and 1, or NULL")
        }
        err <- stats::qnorm(1 - (1 - CI) / 2)
      }
    }
    # latent case
    if (type == "latent") {
      if (is.null(CI)) {
        # if CIs aren't wanted
        ans <- pred(predx = newdata, fit = object, mn = mn, std = FALSE, maxn = maxn)
        colnames(ans) <- "posterior mean"
      } else if (CI == "std") { # if standard deviations are wanted instead
        ans <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        colnames(ans) <- c("posterior mean", "posterior std")
      } else {
        # if they are
        pred <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        upper <- pred[, 1] + err * pred[, 2]
        lower <- pred[, 1] - err * pred[, 2]
        ans <- cbind(pred[, 1], lower, upper)
        colnames(ans) <- c(
          "posterior mean", paste("lower ", round(100 * CI), "% CI", sep = ""),
          paste("upper ", round(100 * CI), "% CI", sep = "")
        )
      }
    } else {
      # response case
      if (is.null(CI)) {
        # if CIs aren't required
        ans <- stats::pnorm(pred(newdata, object, mn, std = FALSE, maxn = maxn))
        colnames(ans) <- "posterior mode"
      } else {
        # if CIs are required
        pred <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        upper <- pred[, 1] + err * pred[, 2]
        lower <- pred[, 1] - err * pred[, 2]
        ans <- stats::pnorm(cbind(pred[, 1], lower, upper))
        colnames(ans) <- c(
          "posterior mode", paste("lower ", round(100 * CI), "% CI", sep = ""),
          paste("upper ", round(100 * CI), "% CI", sep = "")
        )
      }
    }
    ans
  }
