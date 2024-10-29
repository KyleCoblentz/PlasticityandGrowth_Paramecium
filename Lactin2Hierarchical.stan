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
  array[n_genotypes] real<lower = 0,upper = 1> rho_i;
  real<lower = 0> mu_rho;
  real<lower = 0> sigma_rho;
  array[n_genotypes] real Tmax_i;
  real mu_Tmax;
  real<lower = 0> sigma_Tmax;
  array[n_genotypes] real<lower = 0> delta_i;
  real<lower = 0> mu_delta;
  real<lower = 0> sigma_delta;
  array[n_genotypes] real lambda_i;
  real mu_lambda;
  real<lower = 0> sigma_lambda;
  array[n_genotypes] real<lower = 0> sigma_i;
}

// model likelihood and priors

model{
  
  // priors
  mu_rho ~ uniform(0, 1);
  sigma_rho ~ normal(0, 0.5);
  mu_Tmax ~ normal(40, 5);
  sigma_Tmax ~ normal(0, 10);
  mu_delta ~ normal(0, 10);
  sigma_delta ~ normal(0, 5);
  mu_lambda ~ normal(0, 2.5);
  sigma_lambda ~ normal(0, 1);
  
  for(i in 1:n_genotypes){
  sigma_i[i] ~ normal(0, 10);
  } 
  
  // likelihoods
  
    for(i in 1:n_genotypes) {
  
    rho_i[i] ~ normal(mu_rho, sigma_rho);
    
    Tmax_i[i] ~ normal(mu_Tmax, sigma_Tmax);
    
    delta_i[i] ~ normal(mu_delta, sigma_delta);
    
    lambda_i[i] ~ normal(mu_lambda, sigma_lambda);
  
  }
  
  for(i in 1:N){
    
    growth[i] ~ normal(exp(rho_i[genotype[i]]*temp[i]) - exp(rho_i[genotype[i]]*Tmax_i[genotype[i]] - ((Tmax_i[genotype[i]] - temp[i])/delta_i[genotype[i]])) + lambda_i[genotype[i]], sigma_i[genotype[i]]);
    
  }
  
}


