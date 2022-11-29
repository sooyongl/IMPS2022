#' Make stan codes
#' @export
mk_stanmodel <- function(lvmodel = c("2pl","grm","gpcm")) {

  lvmodel　<- tolower(lvmodel)
  lvmodel  <- match.arg(lvmodel)

  if(lvmodel == "2pl") {
    stanmodel <- "data{
//Sample sizes
  int<lower=1> nsecWorked;
  int<lower=1> ncov;
  int<lower=1> nstud;
  int<lower=1> nsec;
  int<lower=1> nfac;
  int<lower=0> min_k;
  int<lower=1> max_k;

// Data indices
  int<lower=1,upper=nstud> studentM[nsecWorked];
  int<lower=1,upper=nsec> section[nsecWorked];

// Index for factor loadings
  matrix[nsec, nfac] factoridx;
  int<lower=0> firstitem[nsec];

// Input data
  int<lower=min_k,upper=max_k> grad[nsecWorked];
  matrix[nstud,ncov] X;
  int<lower=0,upper=1> Z[nstud];
  real Y[nstud];
}

parameters{
  // IRT model
  vector[nstud] eta;  // Latent traits
  real<lower=0> sigU; // Latent variable variance

  matrix<lower=0, upper=10>[nsec, nfac] a1_free; // Item Slopes
  real d[nsec];                 // Item intercepts

  vector[ncov] betaU; // Covariate effects on latent variable
  vector[ncov] betaY; // Covariate effects on outcome

  real omega;
  real yint;
  real tau0;
  real tau1;

  real<lower=0> sigY[2]; // Outcome variance
}

transformed parameters {
  matrix<lower=0, upper=10>[nsec, nfac] a1;

  // Factor loading constraints
  for(jjj in 1:nfac) {
    for(jj in 1:nsec) {
      if(factoridx[jj, jjj] != 0) {
        if(firstitem[jj] == 1) { // first loading per factor constrained to 1.
        a1[jj, jjj] = 1;
        } else {
          a1[jj, jjj] = a1_free[jj, jjj];
        }
      } else {
        a1[jj, jjj] = 0;
      }
    }
  };
}

model{
  vector[nstud] muEta;
  vector[nstud] muY;
  real sigYI[nstud];

// Fully Latent Principal Stratification model
// Structural part -----------------
   for(i in 1:nstud){
     muEta[i] = X[i, ]*betaU;
	 muY[i]  = yint+ omega*eta[i] + Z[i] * (tau0 + tau1*eta[i]) + X[i,]*betaY;
	 sigYI[i]=sigY[Z[i]+1];

	 eta[i] ~ normal(muEta[i], sigU);
	 Y[i] ~ normal(muY[i], sigYI[i]);
   };

// Measurement part -----------------
   for(j in 1:nsecWorked) {
     grad[j] ~ bernoulli_logit(d[section[j]] + a1[section[j],1] * eta[studentM[j]]);
   };

// Priors ------------------
// IRT priors
   d ~ normal(0, 1);
   for(i in 1:nsec) {
     for(j in 1:nfac) {
       a1_free[i, j] ~ lognormal(0, 1);
     };
   };

// Priors for structural model
  betaY ~ normal(0, 1);
  betaU ~ normal(0, 1);
  omega ~ normal(0, 1);
  yint ~ normal(0, 1);
  tau0  ~ normal(0, 1);
  tau1 ~ normal(0, 1);
}
"

  } else if(lvmodel == "grm") {
    stanmodel <- "data{
//Sample sizes
  int<lower=1> nsecWorked;
  int<lower=1> nstud;
  int<lower=1> nsec;
  int<lower=0> min_k;
  int<lower=2> max_k;
  int<lower=0> ncov;
  int<lower=1> nfac;

// Data indices
  int<lower=1,upper=nstud> studentM[nsecWorked];
  int<lower=1,upper=nsec> section[nsecWorked];

// Index for factor loadings
  matrix[nsec, nfac] factoridx;
  int<lower=0> firstitem[nsec];

// Input data
  int<lower=min_k,upper=max_k> grad[nsecWorked];
  matrix[nstud,ncov] X;
  int<lower=0,upper=1> Z[nstud];
  real Y[nstud];
}

parameters{
  vector[nstud] eta;
  real<lower=0> sigU;

  matrix<lower=0, upper=10>[nsec, nfac] a1_free; // Item Slopes
  ordered[max_k-1] d[nsec];                        // Item category intercept

  vector[ncov] betaU;
  vector[ncov] betaY;

  real omega;
  real yint;
  real tau0;
  real tau1;

  real<lower=0> sigY[2];
}

transformed parameters{
 matrix<lower=0, upper=10>[nsec, nfac] a1;

// Factor loading constraints
  for(jjj in 1:nfac) {
    for(jj in 1:nsec) {
	  if(factoridx[jj, jjj] != 0) {
        if(firstitem[jj] == 1) {  // first loading per factor constrained to 1.
          a1[jj, jjj] = 1;
        } else {
          a1[jj, jjj] = a1_free[jj, jjj];
        }
      } else {
        a1[jj, jjj] = 0;
      }
    }
  };
}

model{
  vector[nstud] muEta;
  vector[nstud] muY;
  real sigYI[nstud];

// Fully Latent Principal Stratification model
// Structural model ----------------
  for(i in 1:nstud){
 	muEta[i] = X[i, ]*betaU;
	muY[i]  = yint+ omega*eta[i] + Z[i] * (tau0 + tau1*eta[i]) + X[i,]*betaY;
	sigYI[i]=sigY[Z[i]+1];

	eta[i] ~ normal(muEta[i], sigU);
	Y[i] ~ normal(muY[i], sigYI[i]);
  };

// Measurement model ----------------
  for (i in 1:nsecWorked){
    grad[i]~ordered_logistic(a1[section[i],1]*eta[studentM[i]],d[section[i]]);
  };

// Priors ------------------
// IRT priors
  for(i in 1:nsec) {
    for(ii in 1:(max_k-1)) {
        if(ii == 1) {
           d[i, ii] ~ normal(-1, 1);
        }
        if(ii == 2) {
           d[i, ii] ~ normal(0,1);
        }
        if(ii == 3) {
           d[i, ii] ~ normal(1, 1);
        }
    };
    for(j in 1:nfac) {
        a1_free[i, j] ~ lognormal(0, 1);
    };
  };

// Priors for structural model
  betaY ~ normal(0, 1);
  betaU ~ normal(0, 1);
  omega ~ normal(0, 1);
  yint ~ normal(0, 1);
  tau0  ~ normal(0, 1);
  tau1 ~ normal(0, 1);
}
// last line"

  } else if(lvmodel == "gpcm") {
    stanmodel <- "data{
//Sample sizes
  int<lower=1> nsecWorked;
  int<lower=1> nstud;
  int<lower=1> nsec;
  int<lower=0> min_k;
  int<lower=2> max_k;
  int<lower=0> ncov;
  int<lower=1> nfac;

// Data indices
  int<lower=1,upper=nstud> studentM[nsecWorked];
  int<lower=1,upper=nsec> section[nsecWorked];

// Index for factor loadings
  matrix[nsec, nfac] factoridx;
  int<lower=0> firstitem[nsec];

// Input data
  int<lower=min_k,upper=max_k> grad[nsecWorked];
  matrix[nstud,ncov] X;
  int<lower=0,upper=1> Z[nstud];
  real Y[nstud];
}

parameters{
  vector[nstud] eta;
  real<lower=0> sigU;

  matrix<lower=0, upper=10>[nsec, nfac] a1_free; // Item Slopes
  vector[max_k-1] d[nsec];                         // // Item category intercept

  vector[ncov] betaU;
  vector[ncov] betaY;

  real omega;
  real yint;
  real tau0;
  real tau1;

  real<lower=0> sigY[2];
}

transformed parameters{
  matrix<lower=0, upper=10>[nsec, nfac] a1;

// Factor loading constraints
  for(jjj in 1:nfac) {
    for(jj in 1:nsec) {
	  if(factoridx[jj, jjj] != 0) {
        if(firstitem[jj] == 1) {   // first loading per factor constrained to 1.
          a1[jj, jjj] = 1;
        } else {
          a1[jj, jjj] = a1_free[jj, jjj];
        }
      } else {
        a1[jj, jjj] = 0;
      }
    }
  };
}

model{
  matrix[max_k, nsecWorked] p;  // probs of reponse
  matrix[max_k, nsecWorked] s;  // logits of reponse

  vector[nstud] muEta;
  vector[nstud] muY;
  real sigYI[nstud];
  real measure;

// Fully Latent Principal Stratification model
// Structural model -------------------
  for(i in 1:nstud){
 	  muEta[i] = X[i, ]*betaU;
	  muY[i]  = yint + omega*eta[i] + Z[i] * (tau0 + tau1*eta[i]) + X[i,]*betaY;
	  sigYI[i] = sigY[Z[i]+1];

	  eta[i] ~ normal(muEta[i], sigU);
	  Y[i] ~ normal(muY[i], sigYI[i]);
  };

// Measurement model -------------------
  for(i in 1:nsecWorked) {
    s[1,i] = 1; //reference
	measure = 0;
    for(k in 2:max_k) {
		measure = measure + a1[section[i], 1]*eta[studentM[i]] + d[section[i], k-1];
		s[k, i] = s[k-1, i] + exp(measure);
	}
    for(kk in 1:max_k) {
	   p[kk,i] = s[kk,i] / sum(s[,i]);
	}
    grad[i] ~ categorical(p[,i]);
  }

// Priors ------------------
// IRT priors
  for(i in 1:nsec) {
    for(ii in 1:(max_k-1)) {
        if(ii == 1) {
           d[i, ii] ~ normal(-1, 1);
        }
        if(ii == 2) {
           d[i, ii] ~ normal(0,1);
        }
        if(ii == 3) {
           d[i, ii] ~ normal(1, 1);
        }
    };
    for(j in 1:nfac) {
        a1_free[i, j] ~ lognormal(0, 1);
    };
  };

// Priors for structural model
  betaY ~ normal(0,1);
  betaU ~ normal(0,1);
  omega ~ normal(0,1);
  yint ~ normal(0,1);
  tau0  ~ normal(0,1);
  tau1 ~ normal(0,1);
}
// last line"

  } else {
    stop("not supported")
  }
  stanmodel
}
