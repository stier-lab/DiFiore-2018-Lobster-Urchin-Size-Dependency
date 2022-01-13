sink("code/STAN_models/allometric.stan")
cat("
    data {
    int N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int Nind;
    int<lower=1> id[N]; // this will need to be character???
    real<lower=0> mc[Nind];
    real<lower=0> mr[Nind];
    }
    
    parameters {
    real mualphaa;
    real<lower=0> varmualphaa;
    real alphaa[Nind];
    real<lower=0> varaind;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real mualphah;
    real<lower=0> varmualphah;
    real alphah[Nind];
    real<lower=0> varhind;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    }
    
    transformed parameters{
    // likelihood
    
    real<lower=0> a[Nind];
    real<lower=0> h[Nind]; 
    real loga[Nind];
    real logh[Nind];
    real<lower=0, upper=1> prob[N];
    
    for(i in 1:Nind){
    loga[i] = alphaa[i] + beta1a*log(mc[i]) + beta2a*log(mr[i]);
    a[i] = exp(loga[i]);
    logh[i] = alphah[i] + beta1h*log(mc[i]) + beta2h*log(mr[i]);
    h[i] = exp(logh[i]);
    }
    
    for(i in 1:N){
    prob[i] = 1/(1/a[id[i]] + h[id[i]]*initial[i]);
    //prob[i] = a[id[i]] / (1 + a[id[i]]*h[id[i]]*initial[i]);
    }
    
    }
    
    
    
    model {
    
    //priors
    
    // Population level prior on attack rate intercept
    mualphaa ~ normal(-4.5, varmualphaa);
    //varmualphaa ~ uniform(0,10);
    varmualphaa ~ gamma(1,1);    
    
    // Individual level variation in alphaa
    for(i in 1:Nind){
    alphaa[i] ~ normal(mualphaa, varaind);
    }
    //varaind ~ uniform(0,10);
    varaind ~ gamma(2,1);
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    mualphah ~ normal(11, varmualphah);
    //varmualphah ~ uniform(0,10);
    varmualphah ~ gamma(1,1);    
    
    //Individual level variation in handling time interpect 
    for(i in 1:Nind){
    alphah[i] ~ normal(mualphah, varhind);
    }
    varhind ~ uniform(0,10);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(2,1);
    
    // likelihood
    
    for(i in 1:N){
    killed[i] ~ binomial(initial[i], prob[i]);
    }
    
    
    }",fill = TRUE)
sink()