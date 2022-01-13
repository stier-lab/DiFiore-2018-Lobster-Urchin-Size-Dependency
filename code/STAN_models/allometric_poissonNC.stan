
    data {
    int N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int Nind;
    int id[N]; 
    real<lower=0> mc[N];
    real<lower=0> mr[N];
    }
    
    parameters {

    // These are all the hyperparameters estimated across all individuals
    real alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    // These are parameters associated with variation in attack rates and handling times between individuals. In other words they are the random effect of individual. 
    real<lower=0> sigmaa;
    vector[Nind] zza; 
    real<lower=0> sigmah; 
    vector[Nind] zzh;
    }
    
    transformed parameters{
    
    vector[N] a;
    vector[N] h; 
    vector[N] lambda; // This is the estimated number of prey eaten in any trial i. Lambda is the parameters for a poisson distribution.
    
    for(i in 1:N){
    a[i] = exp(alphaa + beta1a*log(mc[i]) + beta2a*log(mr[i]) + zza[id[i]]*sigmaa);
    h[i] = exp(alphah + beta1h*log(mc[i]) + beta2h*log(mr[i]) + zzh[id[i]]*sigmah);
    lambda[i] = a[i]*initial[i] / (1 + a[i]*h[i]*initial[i]);
    }
    

    }
    
    
    model {
    
    //priors
    
    // individual variation
    sigmaa ~ gamma(2,1);
    for(i in 1:Nind){
    zza[i] ~ normal(0, sigmaa);
    }
    
    sigmah ~ gamma(2,1);
    for(i in 1:Nind){
    zzh[i] ~ normal(0, sigmah);
    }
    
    // Population level prior on attack rate intercept
    alphaa ~ normal(0, varalphaa);
    varalphaa ~ gamma(1,1);    
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(1,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(1,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(0, varalphah);
    varalphah ~ gamma(1,1);    
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(1,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(1,1);
    
    // likelihood
    
    killed ~ poisson(lambda);
    
    }
