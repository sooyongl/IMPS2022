#' S3 generic for Outcome data generation
#'
genOutcome <- function(Data, ...) { # Data = sim_info
  UseMethod("genOutcome", Data)
}

#' generate outcome (Y)
#'
genOutcome.default <- function(Data) {

  N      <- Data$N
  nsec   <- Data$nsec
  theta  <- Data$theta
  xdata  <- Data$x
  nfac   <- Data$nfac

  R2Y    <- Data$R2Y
  linear <- Data$linear
  ydist  <- Data$ydist

  omega  <- Data$omega # round(runif(nfac, 0.1, 0.3),3)
  tau0   <- Data$tau0  # round(runif(1, 0.2, 0.4),3)
  tau1   <- Data$tau1  # round(runif(nfac, -0.2, -0.1),3)

  section  <- Data$section
  studentM <- Data$studentM
  grad     <- Data$grad
  lv.resp  <- Data$lv.resp

  a_idx <- gen_a_idx(nsec, nfac)
  fi_idx <- detect_firstitem(a_idx)

  Y.R2  <- R2Y
  eta   <- theta
  x1    <- xdata[,"x1"]
  x2    <- xdata[,"x2"]
  Z     <- rep(c(1,0), each=N/2)

  n.eta <- ifelse(!is.null(dim(eta)),  ncol(eta), 1)

  Y <-
    tau0*Z +
    matrix(eta, ncol=n.eta)%*%matrix(omega) +
    (Z*matrix(eta, ncol=n.eta))%*%tau1 +
    1*x1 +
    0.5*x2

  exp_var  <- var(Y)
  unex_var <- (1 - Y.R2)*exp_var / Y.R2
  Y <- Y + MASS::mvrnorm(N, 0, sqrt(unex_var), empirical = T)

  Y <- Y - mean(Y)


  Data$stan_dt <- list(
    # data info
    nsecWorked = length(section),
    nstud      = N,
    nsec       = nsec,
    nfac       = nfac,
    min_k      = min(grad),
    max_k      = max(grad),
    ncov       = ncol(xdata),
    # index
    studentM     = studentM,
    section      = section,
    factoridx    = a_idx,
    firstitem    = fi_idx,
    # data
    grad   = grad,
    X      = xdata,
    Z      = Z,
    Y      = c(Y)
  )

  return(Data)
}
