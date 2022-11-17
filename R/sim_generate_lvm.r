#' Set latent variable model information
#'
genLVinfo <- function(sim_info) {

  nitem   <- sim_info$nsec
  nfac    <- sim_info$nfac
  lvmodel <- sim_info$lvmodel

  if(sim_info$lvmodel %in% c("rasch","2pl","3pl","gpcm","pcm","grm","ggrm")){

    ipar <- genIRTpar(nitem, nfac, lvmodel, ncat = 4)

  } else {
    stop("LV models are incorrect")
  }

  sim_info$lvinfo$ipar <- ipar

  return(sim_info)
}

#' Generate IRT parameters
#'
genIRTpar <- function(nitem=25, nfac=1, lvmodel, ncat = 4) {

  lvmodel <- tolower(lvmodel)

  if(ncat <= 1) {
    stop("the number of cateories should be at least 2")
  } else if(ncat == 2 & lvmodel %in% c("grm","gpcm")) {
    stop("For GRM and GPCM, cateories should be at least 3")
  }

  a_list <- gen_a(nitem, nfac)
  a <- a_list$a
  a_idx <- a_list$a_idx

  if(lvmodel %in% c("grm","gpcm","ggrm","pcm")) {
    if(lvmodel == "gpcm") {
      diffs <- t(apply(matrix(runif(nitem * (ncat-1), 0.5, 1), nitem), 1, cumsum))
      d <- -(diffs - rowMeans(diffs))
      d <- -1*d
    }
    if(lvmodel == "grm") {
      diffs <- t(apply(matrix(runif(nitem * (ncat-1), 0.5, 1), nitem), 1, cumsum))
      d <- -(diffs - rowMeans(diffs))
    }

    colnames(d) <- paste0("d",1:ncol(d))
    ipar <- data.frame(a, d)

  } else if(lvmodel %in% c("rasch","1pl","2pl","3pl")) {

    g <- 0
    d <- rnorm(nitem)

    if(lvmodel %in% c("1pl","rasch")) {
      a[which(a_idx == 1)] <- 1
    } else if(lvmodel == "3pl") {
      g = runif(nitem, 0, 0.2)
    }

    ipar <- data.frame(a, d, g)
  }

  return(ipar)
}

#' Generate LV model data
#'
genLVM <- function(info) { # info = sim_info
  N      <- info$N
  nsec   <- info$nsec
  nfac   <- info$nfac
  lambda <- info$lambda

  lv.gen.dt <- generateLV(info)

  lv.par <- lv.gen.dt$lv.par
  lv.resp <- lv.gen.dt$resp

  total_N <- N/2

  nworked <- rep(floor(nsec * lambda), total_N)

  studentM <- do.call("c", lapply(seq(total_N),
                                  function(n) rep(n,each=nworked[n])))
  section <- do.call("c", lapply(seq(total_N),
                                 function(n) {
                                   sort(sample(1:nsec, nworked[n],
                                               replace = FALSE))}))
  ss <- cbind(studentM, section)
  grad <- sapply(1:dim(ss)[1], function(n) lv.resp[ss[n,1], ss[n,2]] )

  res <- list(
    lv.par   = lv.par,
    lv.resp  = lv.resp,
    grad     = grad,
    studentM = studentM,
    section  = section
  )

  info <- structure(append(info, res), class = attr(info, "class"))

  return(info)
}

#' S3 generic for latent model data generation
#'
generateLV <- function(info, ...) {
  UseMethod("generateLV", info)
}

#' methods for rasch model
#'
generateLV.rasch <- function(.x, ...) {.Class <- "irt"; NextMethod()}

#' method for 1PL model
#'
generateLV.1pl   <- function(.x, ...) {.Class <- "irt"; NextMethod()}

#' method for 2PL model
#'
generateLV.2pl   <- function(.x, ...) {.Class <- "irt"; NextMethod()}

#' method for 3PL model
#'
generateLV.3pl   <- function(.x, ...) {.Class <- "irt"; NextMethod()}

#' method for GPCM model
#'
generateLV.gpcm   <- function(.x, ...) {.Class <- "irt"; NextMethod()}
generateLV.pcm   <- function(.x, ...) {.Class <- "irt"; NextMethod()}
#' method for GRM
#'
generateLV.ggrm   <- function(.x, ...) {.Class <- "irt"; NextMethod()}
generateLV.grm   <- function(.x, ...) {.Class <- "irt"; NextMethod()}

#' method for all IRT model
#'
#' @examples
#'
#' lvmodel <- "gpcm"
#' ipar <- genIRTpar(20, ncat = 3, 2, lvmodel)
#' eta <- MASS::mvrnorm(100, rep(0, 2), matrix(c(1,0,0,1),ncol=2))
#' generateLV.irt(lvmodel, eta, ipar)
#'
generateLV.irt <- function(info) {

  theta <- as.matrix(info$theta, ncol = info$nfac);
  nitem <- info$nsec;
  lv_info <- info$lvinfo
  lvmodel <- tolower(info$lvmodel)
  ipar <- lv_info$ipar

  a <- ipar[grep("a",names(ipar))]
  d <- ipar[grep("d|b",names(ipar))]
  guess <- 0

  stopifnot(is.data.frame(ipar))
  stopifnot(ncol(a) == ncol(theta))

  N    <- nrow(theta)
  nfac <- ncol(theta)

  lvmodel <- switch(lvmodel,
                    "rasch" = "dich",
                    "1pl" = "dich",
                    "2pl" = "dich",
                    "3pl" = "dich",
                    "gpcm" = "gpcm",
                    "pcm" = "gpcm",
                    "grm" = "graded",
                    "ggrm" = "graded")

  if(lvmodel == "dich") { guess <- ipar[,grep("g",names(ipar))] }

  resp <- simIRTdata(a = as.matrix(a),
                  d = as.matrix(d),
                  guess = as.vector(guess),
                  N = N,
                  theta = theta,
                  itemtype = lvmodel)

  if(lvmodel != "dich")
    resp <- resp + 1


  return(list(resp = data.frame(resp), lv.par = ipar))
}

#' Generate IRT data
#'
simIRTdata <- function(a, d, guess, N, theta, itemtype) {

  if(itemtype != "gpcm"){
    resp <- mirt::simdata(
      a = a,
      d = d,
      guess = guess,
      N = N,
      Theta = theta,
      itemtype = itemtype)

  } else {
    resp <- simData.pcm(
      a = a,
      d = d,
      theta = theta
    )
  }

  # if(itemtype == "dich"){
  #   resp <- simData.dich(a, d, guess, N, theta)
  #
  # } else if(itemtype == "gpcm") {
  #   resp <- simData.gpcm(a, d, N, theta)
  #
  # } else if(itemtype == "graded") {
  #
  # }

  resp
}

#' Generate GPCM data
#'
simData.pcm <- function(a,d,theta) {

  N     <- nrow(theta)
  nfac  <- ncol(a)
  nitem <- nrow(a)
  ncat <- ncol(d) + 1

  difficulty = rowMeans(d)
  steps <- apply(d, 2, function(x) x - difficulty)
  steps <- cbind(difficulty, steps)

  # Create an empty response matrix
  res <- matrix(NA, nrow = N, ncol = nitem)
  # for(h in 1:nfac) { # h = 1
  for(k in 1:N) { # k = 1
    for(i in 1:nitem) { # i = 1
      measure=0
      p <- vector()
      p[1] <- 1

      for(j in 2:ncat) {
        measure <- measure + a[i, 1]*theta[k, 1] - d[i, j-1]
        p[j] <- p[(j-1)] + exp(measure)
      }

      U <- runif(1, 0, 1)
      U = U * p[ncat]

      for(j in 1:ncat) {
        if(U <= p[j]) {
          res[k,i] <- (j-1)
          break
        }
      }
    }
  }
  # }
  return(res)
}
