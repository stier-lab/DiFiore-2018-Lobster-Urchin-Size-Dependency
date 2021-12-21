

// Decision about data:
  // * Filter by number of cases? 
  // * Not all counties start reporting the same day
// * Throw out counties that have <= 5 days pre-lockdown
// * Are we going to be able to calculate active cases pre-lockdown if we don't have data for 12 days before lockdown?
//     - Tentative assumption: In very early parts of data set (pre-lockdown), assume that unrecorded cases, means 0 cases
//     - Reasoning: there were already false 0s in the data that will have a strong effect
//     - And it makes things a lot easier because we don't need counties to have 17 days of data pre-lockdown to be included in the dataset

data {
  // Information we need for coding purposes
  // Stan can't read in the data until it knows the size of the matrices
  int<lower=1> counties; // The number of counties in the data set
  int<lower=1> timepointsPre; // The number of time points (days with new cases data) pre-lockdown
  int<lower=1> timepointsPost; // Number of time points post-lockdown
  
  // Statistical data
  
  // This reads in a matrix with rows 
  int<lower=0> activeCasesPre[counties, timepointsPre];  // Active cases = Cases for the last 12 days - deaths
  int<lower=0> newCasesPre[counties, timepointsPre];     // Case counts on that day
  
  int<lower=0> activeCasesPost[counties, timepointsPost];  // Active cases = Cases for the last 12 days - deaths
  int<lower=0> newCasesPost[counties, timepointsPost];     // Case counts on that day
  
}

parameters {
// parameters with beta in their name are mean slopes (on log scale) of the effect of active cases on new cases
//vector[counties] betaPre; // mean effect of active cases for each county pre-lockdown
//vector[counties] betaPost; // mean effect of active cases for each county post-lockdown
matrix[counties, 2] beta;
vector[2] betaBar;
//real betaPreBar; // mean effect of active cases across all counties, pre-lockdown
//real betaPostBar; // mean effect of active cases across all counties, post-lockdown

// parameters with sigma and R in their name are related to standard deviations/covariances
real<lower=0> sigmaPre; // standard deviation of the distribution of betaPre (mean effect of active cases pre-lockdown) between counties
real<lower=0> sigmaPost; // standard deviation the distribution of betaPost between counties
matrix[2, 2] R; // used to get the covariance terms between sigmaPre and sigmaPost for each county (cov=R[1,2] * sigmaPre * sigmaPost)

// parameters with eta in their name refer to noise in the relationship between active cases and new cases
matrix[counties, timepointsPre] etaPre; // noise in the effect of active cases on new cases, varies with time point
matrix[counties, timepointsPost] etaPost; // noise, same distribution as etaPre. Separated into etaPre (pre-lockown) and etaPost for convenience
real<lower=0> sigmaEta; // standard deviation of noise parameter
}

// There are a couple things I'm not sure I did right in this block
// * equation for lamdba
// * decision to define beta, betaBar, and sigma here
transformed parameters {
  
  // lambda refers to Poisson parameters for the distribution of active cases.
  // There is a different value of lambda for each county and timepoint, because lambda depends both on the county
  // and on the number of active cases at that time point.
  // lambda also depends whether lockdown has occurred
  matrix<lower=0>[counties, timepointsPre] lambdaPre; // poisson parameter for each county and each pre-lockdown timepoint
  matrix<lower=0>[counties, timepointsPost] lambdaPost; // poisson parameter for each county and each post-lockdown timepoint
  
  // Here I've decided (perhaps wrongly) to make some variables that are matrices or vectors containing other parameters
  // I did this because later on we need to write a multivariate normal distribution using these parameters, and I couldn't
  // figure out how to do that without writing the things involved in the distribution as vectors/matrices.
  // In particular, we need:
    // * a matrix with the values of beta for each county pre- and post-lockdown
  // * a vector of the statewide mean beta value pre- and post-lockdown
  // * a matrix whose diagonals are the pre- and post-lockdown standard deviation in the distribution from which the county-level beta values are drawn
  //matrix[counties, 2] beta; // matrix of betaPre (1st column) and betaPost (2nd column) for all counties
  // because I can't figure out how to write the model block without it
  //vector[2] betaBar; //matrix of betaPreBar and betaPostBar
  matrix[2,2] sigma; // matrix of [sigmaPre 0,  0 sigmaPost] for use in generating Sigma for multivariate normal in model
  matrix[2, 2] S; // VCV matrix for multivariate normal distribution from which betaPre and betaPost are drawn for each county
  
  // lambda (poisson parameter) for each county and timepoint
  // exp(county-level effect of active cases on new cases + noise for each observation of cases)*activecases that day
  for (i in 1:counties) {
  for (j in 1:timepointsPre) {
  //lambdaPre[i, j] = exp(betaPre[i] + etaPre[i, j]) * activeCasesPre[i, j];
  lambdaPre[i, j] = exp(beta[i,1] + etaPre[i, j]) * activeCasesPre[i, j];
  }
  
  for (k in 1:timepointsPost) {
  //lambdaPost[i, k] = exp(betaPost[i] + etaPost[i, k]) * activeCasesPost[i, k];
  lambdaPost[i, k] = exp(beta[i,2] + etaPost[i, k]) * activeCasesPost[i, k];
  }
  }
  
  // making vectors and matrices that compile several parameters into one variable
  //for (i in 1:counties) {
  //    beta[i, 1] = betaPre[i];
  //    beta[i, 2] = betaPost[i];
  //}
  
  //betaBar[1] = betaPreBar;
  //betaBar[2] = betaPostBar;
  
  sigma[1,1] = sigmaPre;
  sigma[1, 2] = 0;
  sigma[2, 1] = 0;
  sigma[2, 2] = sigmaPost;
  
  S = sigma * R * sigma;
  
}

model {

//Priors go here
// standard deviations of the statewide distribution of beta values
sigmaPre ~ exponential(1);
sigmaPost ~ exponential(1);

// statewide mean beta
betaBar[1] ~ normal(0, 2);
betaBar[2] ~ normal(0, 2);

// off-diagonals give correlation between pre-lockdown and post-lockdown beta
R ~ lkj_corr(2);

// standard deviation of noise
sigmaEta ~ exponential(1);

// distribution of beta
for (i in 1:counties) {
beta[i]' ~ multi_normal(betaBar, S);
}

//I guess technically these are priors too!
  for (i in 1:counties) {
    for (j in 1:timepointsPre) {
      etaPre[i, j] ~ normal(0, sigmaEta);
    }
    
    for (k in 1:timepointsPost) {
      etaPost[i, k] ~ normal(0, sigmaEta);
    }
  }

//And now we need the likelihood which is Poisson distributed with the calculated lambda
for (i in 1:counties) {
  for (j in 1:timepointsPre) {
    activeCasesPre[i, j] ~ poisson(lambdaPre[i, j]);
  }
  
  for (k in 1:timepointsPost) {
    activeCasesPost[i, k] ~ poisson(lambdaPost[i, k]);
  }
}

}

   
