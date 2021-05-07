
    data {
    int<lower=1> N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int<lower=1> Nind;
    int<lower=1> id[N]; // this will need to be character???
    real<lower=0> mc[N];
    real<lower=0> mr[N];
    }
    
    parameters {
    real<lower=0> alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real<lower=0> alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    }
    
    transformed parameters{
    // likelihood

    real<lower=0, upper=1> prob[N];
    
    for(i in 1:N){
    prob[i] = 1/((1/(alphaa*(mc[i]^beta1a)*(mr[i]^beta2a))) + (alphah*(mc[i]^beta1h)*(mr[i]^beta2h)*initial[i]));
    }
    
    }
    
    
    model {
    
    //priors
    
    // Population level prior on attack rate intercept
    log(alphaa) ~ normal(0, varalphaa);
    target += -log(alphaa);
    varalphaa ~ uniform(0,10);
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    log(alphah) ~ normal(0, varalphah);
    target += -log(alphah);
    varalphah ~ uniform(0,10);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(2,1);
    
    // likelihood
    
    for(i in 1:N){
    killed[i] ~ binomial(initial[i], prob[i]);
    }
    
    }
