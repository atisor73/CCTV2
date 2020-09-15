functions {
    real a_theor_exp(real ao, real k, real t) {
        return ao * exp(k*t);}
    
    real a_theor_lin(real ao, real k, real t) {
        return ao + k*t;}
}

data {
    int len_df1;
    int len_df2;
    int nI_cycles;
    int nII_cycles;
    
    real timesI[len_df1];
    real areasI[len_df1];
    int cyclesI[len_df1];
    
    real timesII[len_df2];
    real areasII[len_df2];
    int cyclesII[len_df2];
}

parameters {
    real sigma_tilde;
    real<lower=0, upper=0.06> mu_sigma;
    real<lower=0, upper=0.1> sig_sigma;
    
    real<lower=0, upper=0.3> sig_ao;
    real<lower=0, upper=0.002> sig_k;
    
    real<lower=0> mu_aoI;
    real<lower=0> mu_aoII;
    real mu_kI_tilde;
    real mu_kII_tilde;
    
    real aoI_tilde[nI_cycles];
    real aoII_tilde[nII_cycles];
    real kI_tilde[nI_cycles];
    real kII_tilde[nII_cycles];
}

transformed parameters {
    real sigma;
    real mu_kI;
    real mu_kII;
    real aoI[nI_cycles];
    real aoII[nII_cycles];
    real kI[nI_cycles];
    real kII[nII_cycles];
    real muI[len_df1];
    real muII[len_df2];
    
    
    sigma = mu_sigma + sigma_tilde * sig_sigma;
    mu_kI = 0.0116 + 0.0032 * mu_kI_tilde;
    mu_kII = 0.0116 + 0.0032 * mu_kII_tilde;
    
    for (i in 1:nI_cycles) {
        aoI[i] = mu_aoI + aoI_tilde[i] * sig_ao;
        kI[i] = mu_kI + kI_tilde[i] * sig_k;}
    for (i in 1:nII_cycles) {
        aoII[i] = mu_aoII + aoII_tilde[i] * sig_ao;
        kII[i] = mu_kII + kII_tilde[i] * sig_k;}
    
    for (i in 1:len_df1) {
        muI[i] = a_theor_lin(aoI[cyclesI[i]+1], kI[cyclesI[i]+1], timesI[i]);}
    for (i in 1:len_df2) {
        muII[i] = a_theor_lin(aoII[cyclesII[i]+1], kII[cyclesII[i]+1], timesII[i]);}
}

model {
    sigma_tilde ~ normal(0,1);
    mu_sigma ~ exponential(100);
    sig_sigma ~ exponential(500);
    
    sig_ao ~ exponential(20);
    sig_k ~ exponential(1200);
    
    mu_aoI ~ lognormal(0.15,0.15);
    mu_aoII ~ lognormal(0.15,0.15);
    mu_kI_tilde ~ normal(0,1);
    mu_kII_tilde ~ normal(0,1);
    
    for (i in 1:nI_cycles) {
        aoI_tilde[i] ~ normal(0,1);
        kI_tilde[i] ~ normal(0,1);}
    
    for (i in 1:nII_cycles) {
        aoII_tilde[i] ~ normal(0,1);
        kII_tilde[i] ~ normal(0,1);}
    
    areasI ~ normal(muI, sigma);
    areasII ~ normal(muII, sigma);
    
}


generated quantities {
    real areasI_post[len_df1];
    real areasII_post[len_df2];
    
    areasI_post = normal_rng(muI, sigma);
    areasII_post = normal_rng(muII, sigma);
}



































