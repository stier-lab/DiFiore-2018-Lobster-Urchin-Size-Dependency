
    data {
    int<lower=1> Nind;
    int<lower=0> Ndensities;
    int killed[Nind, Ndensities];
    int initial[Nind, Ndensities];
    real<lower=0> mc[Nind];
    }
    
    parameters {
    real<lower=0> alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real<lower=0> alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    
    }
    
    transformed parameters{
    // likelihood
    
    vector[Nind] loga;
    vector[Nind] logh;
    vector[Nind] a;
    vector[Nind] h;
    real<lower=0, upper=1> prob[Nind, Ndensities];
    
    for(i in 1:Nind){
    loga[i] = alphaa + beta1a*log(mc[i]);
    a[i] = exp(loga[i]);
    logh[i] = alphah + beta1h*log(mc[i]);
    h[i] = exp(logh[i]);
    
    for(j in 1:Ndensities){
    prob[i,j] = 1/(1/a[i] + h[i]*initial[i,j]);
    }
    
    }
    
    }
    
    
    model {
    
    //priors
    
    // Population level prior on attack rate intercept
    alphaa ~ normal(0, varalphaa);
    //log(alphaa) ~ normal(0, varalphaa);
    //target += -log(alphaa);
    varalphaa ~ exponential(1);
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(0, varalphah);
    //log(alphah) ~ normal(0, varalphah);
    //target += -log(alphah);
    varalphah ~ exponential(1);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    // likelihood
    
    for(i in 1:Nind){
    for(j in 1:Ndensities){
    killed[i,j] ~ binomial(initial[i,j], prob[i,j]);
    }
    }
    
    }
