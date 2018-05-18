data {
  int<lower=0> T;   // # time points (equally spaced)
    vector[T] y;      // mean corrected return at time t
}
parameters {
  real G_0;                     // initial
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  real<lower=0> nu;  
  vector[T] h_std;             // std log volatility time t
  vector[T] G_std;   
}

transformed parameters {
  vector[T] h;   // log volatility at time t
  vector[T] G;
  G[1] = G_0 +sqrt(nu)*G_std[1];
  h[1] = mu +sqrt(sigma)*h_std[1] / sqrt(1 - phi * phi);
  for (t in 2:T){
    G[t] = G[t-1] + sqrt(nu)*G_std[t];
    h[t] = mu + phi * (h[t-1] - mu) +  tanh(G[t]) * y[t-1] *exp(-h[t-1]/2)*sqrt(sigma)                      + h_std[t] * sqrt(sigma)*sqrt(1 - tanh(G[t]) * tanh(G[t])); 
  }
}


model {
  sigma ~ inv_gamma(2.5, 0.025); 
  nu ~ inv_gamma(2.5, 0.025); 
  phi~uniform(-1,1);                            //(phi+1)/2 ~ beta(20,1.5);   
  mu ~ normal(0,100);  
  G_0 ~ normal(0,100);  
  h_std ~ normal(0,1);
  G_std ~ normal(0,1);
  y ~ normal(0, exp(h / 2));
}