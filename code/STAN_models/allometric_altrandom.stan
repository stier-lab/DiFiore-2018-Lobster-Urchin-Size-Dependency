
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
    real alphaa;
    //real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real alphah;
    //real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    real<lower=0> sigma;
    vector[Nind] zz;

    
    }
    
    transformed parameters{

    vector[N] a;
    vector[N] h; 
    vector[N] prob;

    for(i in 1:N){
    a[i] = exp(alphaa + beta1a*log(mc[i]) + beta2a*log(mr[i]) + zz[id[i]]);
    h[i] = exp(alphah + beta1h*log(mc[i]) + beta2h*log(mr[i]) + zz[id[i]]);
    prob[i] = 1/(1/a[i] + h[i]*initial[i]);
    }

    }
    

    model {
    
    //priors

    // individual variation
    sigma ~ gamma(2,1);
    for(i in 1:Nind){
    zz[i] ~ normal(0, sigma);
    }
    
    // Population level prior on attack rate intercept
    alphaa ~ normal(-4.5, 1);
    //varalphaa ~ gamma(1,1);    
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(11, 1);
    //varalphah ~ gamma(1,1);    
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(2,1);
    
    // likelihood
 
    killed ~ binomial(initial, prob);
    
    }
