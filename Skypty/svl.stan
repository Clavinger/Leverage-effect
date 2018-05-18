data {
  int<lower=0> T;   // # time points (equally spaced)
    vector[T] y;      // mean corrected return at time t
}
parameters {
  real<lower=-1,upper=1> rho;                     // asymmetry 
  real mu;                     // mean log volatility
  real<lower=0,upper=1> phiStar;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h_std;             // std log volatility time t
}

transformed parameters {
  vector[T] h;   // log volatility at time t
  real<lower=-1,upper=1> phi; 
  phi=2*phiStar-1;
  h = h_std * sqrt(sigma)*sqrt(1 - rho * rho);
  h[1] = h[1]+mu;
  for (t in 2:T){
    h[t] = h[t]+mu + phi * (h[t-1] - mu) +  sqrt(sigma)*rho *y[t-1] *exp(-h[t-1]/2) ;
  }
}


model {
  sigma ~ inv_gamma(2.5, 0.025); 
  phiStar ~ beta(20,1.5);   
  mu ~ normal(0,100);  
  h_std ~ normal(0,1);
  rho ~ uniform(-1,1);
  y ~ normal(0, exp(h / 2));
}
