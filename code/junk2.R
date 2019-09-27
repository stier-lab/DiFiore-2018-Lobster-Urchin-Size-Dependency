jagsscript = cat("

                 model{
                 
                 # This model assumed 4 temperature treatments an 11, 16, 21, 26. This model has three levels of heirarchy: a population level estiamte of the FR, a treatment level estimate of the FR, and an indivual level estiamte of the FR.
                 
                 PRIORS
                 
                 # # Global population estiamte
                 #     #a
                 #     logit_mu.a ~ dnorm(0, 1)
                 #     mu.a <- exp(logit_mu.a)/(1+exp(logit_mu.a))
                 # 
                 #     #h
                 #     log_mu.h ~ dnorm(0, 1)
                 #     mu.h <- log(log_mu.h)
                 
                 # Treatment level effects
                 for(i in 1:Ntreats){
                 #a
                 t.logit.a[i] ~ dnorm(mu.a, t.tau.a)
                 t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))
                 
                 #h
                 t.log.h[i] ~ dnorm(mu.h, t.tau.h)
                 t.h[i] <- log(t.log.h[i])
                 
                 }
                 
                 
                 # Individual level effects
                 for(i in 1:num.ind){
                 
                 # a
                 loga[i] ~ dnorm(t.a[tind[i]], tau_int.a)
                 a[i] <- exp(loga[i])/(1+exp(loga[i]))
                 
                 # h
                 logh[i] ~ dnorm(t.h[tind[i]], tau_int.h)
                 h[i] <- log(logh[i])
                 }
                 
                 
                 # for(i in 1:n){
                 # #a
                 # t.logit.a[t[i]] ~ dnorm(0, t.tau.a)
                 # t.a[t[i]] <- exp(t.logit.a[t[i]])/(1+exp(t.logit.a[t[i]]))
                 # logit.a[id[i]] ~ dnorm(t.a[t[i]], tau_int.a)
                 # a[id[i]] <- exp(logit.a[id[i]])/(1+exp(logit.a[id[i]]))
                 # 
                 # 
                 # #h
                 # t.logit.h[t[i]] ~ dnorm(0, t.tau.h)
                 # t.h[t[i]] <- exp(t.logit.h[t[i]])/(1+exp(t.logit.h[t[i]]))
                 # logit.h[id[i]] ~ dnorm(t.h[t[i]], tau_int.h)
                 # h[id[i]] <- exp(logit.h[id[i]])/(1+exp(logit.h[id[i]]))
                 #                  
                 # }
                 
                 # for(i in 1:n){
                 # t.a[t[i]] ~ dnorm(0, t.tau.a)
                 # a[id[i]] ~ dnorm(t.a[t[i]], tau_int.a)
                 # 
                 # t.h[t[i]] ~ dnorm(0, t.tau.h)
                 # h[id[i]] ~ dnorm(t.h[t[i]], tau_int.h)
                 # }
                 
                 # functional response likelihood
                 
                 # # For accurate estimates of treatment level effects
                 # for(i in 1:n){
                 #   
                 #    prob[i] <- max(0.0001,min(0.9999,1/(1/t.a[t[i]] + t.h[t[i]]*initial[i])))
                 #    
                 #    killed[i] ~ dbin(prob[i],initial[i])
                 #    
                 # }
                 
                 # For accurate estimates of individual level effects
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 
                 # Variances for all levels
                 sigma_int.a ~dunif(0,10)
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a)
                 sigma_int.h ~dunif(0,10)
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)
                 
                 sigma_t.a ~ dunif(0,10)
                 t.tau.a <- 1/(sigma_t.a*sigma_t.a)
                 sigma_t.h ~ dunif(0,10)
                 t.tau.h <- 1/(sigma_t.h*sigma_t.h)
                 
                 }
                 
                 ", file = here("code", "heirarchical_jagsSAM.txt"))
