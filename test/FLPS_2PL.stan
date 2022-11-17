data{
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
// last line blank
