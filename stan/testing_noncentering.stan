functions {
    real a_theor(real ao, real k) {
        /*return ao + k*(100);*/
        return ao * exp(k*100); // look @ end of growth for now
    }
}

parameters {
    real<lower=0> mu_ao;
    real<lower=0> sig_ao;
    real<lower=0> sig_k;
    real<lower=0> mu_sigma;
    real<lower=0> sig_sigma;
    
    real ao_tilde;     // bounding tildes? -> probs not
    real k_tilde;
    real mu_k_tilde;
    real sigma_tilde;
    real a_tilde;
}

transformed parameters {
    real ao;
    real mu_k;
    real k;
    real sigma;
    real mu;
    real a;
    
    ao    = mu_ao    + ao_tilde * sig_ao;
    mu_k  = 0.0116   + 0.0032  * mu_k_tilde;
    k     = mu_k     + k_tilde * sig_k;
    sigma = mu_sigma + sigma_tilde * sig_sigma;
    mu    = a_theor(ao, k);
    a     = mu       + a_tilde * sigma;
}

model {
    ao_tilde     ~ normal(0,1);
    mu_ao        ~ lognormal(0.15,0.15);
    sig_ao       ~ exponential(20);
    
    k_tilde      ~ normal(0,1);
    mu_k_tilde   ~ normal(0,1);
    sig_k        ~ exponential(1000);
    
    sigma_tilde  ~ normal(0,1);
    mu_sigma     ~ exponential(100);
    sig_sigma    ~ exponential(500);
    
    a_tilde      ~ normal(0, 1);
}






