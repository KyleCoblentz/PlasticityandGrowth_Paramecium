//
// Stan program to perform hierarchical TPC fitting 
// of the Lactin 2 model
//
//


// functions{
//   
//   real lactin2(real T, real rho, real Tmax, real delta, real lambda){
//     
//     real r;
//     r = exp(rho*T) - exp(rho*Tmax - (Tmax - T)/delta) + lambda;
//     return r;
//     
//   }
//   
// }

// data block

data{
  int N;
  int n_genotypes;
  array[N] real growth;
  array[N] real temp;
  array[N] int genotype;
}

// parameters

parameters{
  array[n_genotypes] real<lower = 0> rho_i;
  real<lower = 0> mu_rho;
  real<lower = 0> sigma_rho;
  array[n_genotypes] real<lower = 0> Tmax_i;
  real<lower = 0> mu_Tmax;
  real<lower = 0> sigma_Tmax;
  array[n_genotypes] real delta_i;
  real mu_delta;
  real<lower = 0> sigma_delta;
  array[n_genotypes] real lambda_i;
  real mu_lambda;
  real<lower = 0> sigma_lambda;
  real<lower = 0> sigma;
}

// model likelihood and priors

model{
  
  // priors
  mu_rho ~ normal(0, 25);
  sigma_rho ~ exponential(1);
  mu_Tmax ~ normal(40,25);
  sigma_Tmax ~ exponential(1);
  mu_delta ~ normal(0, 25);
  sigma_delta ~ exponential(1);
  mu_lambda ~ normal(0, 25);
  sigma_lambda ~ exponential(1);
  
    sigma ~ exponential(1);
  
  // likelihoods
  
  for(i in 1:N){
    
    growth[i] ~ normal(exp(rho_i[genotype[i]]*temp[i]) - exp(rho_i[genotype[i]]*Tmax_i[genotype[i]] - ((Tmax_i[genotype[i]] - temp[i])/delta_i[genotype[i]])) + lambda_i[genotype[i]], sigma);
    
  }
  
  for(i in 1:n_genotypes) {
  
    rho_i[i] ~ normal(mu_rho, sigma_rho);
    
    Tmax_i[i] ~ normal(mu_Tmax, sigma_Tmax);
    
    delta_i[i] ~ normal(mu_delta, sigma_delta);
    
    lambda_i[i] ~ normal(mu_lambda, sigma_lambda);
  
  }
  
}


