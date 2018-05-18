data {
  int<lower=0> T;   // # time points (equally spaced)
    vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h_std;             // std log volatility time t
}
transformed parameters {
  vector[T] h;   // log volatility at time t
  h = h_std * sigma;
  h[1] = h[1] / sqrt(1 - phi * phi);
  h = h + mu;
  for (t in 2:T){
    h[t] = h[t] + phi * (h[t-1] - mu);
  }
}

model {
  sigma*sigma ~ inv_gamma(2.5, 0.025); //do zmiany
  (phi+1)/2 ~ beta(20,1.5);   
  mu ~ normal(0,100);  
  h_std ~ normal(0,1);
  y ~ normal(0, exp(h / 2));
}
