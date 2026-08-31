data {
  int<lower=0> I; //number of individuals 
  int<lower=0> T; //number of traits
  int<lower=0> F; //number of fixed effects
  matrix[I,F] X; //fixed effect matrix  
  matrix[I,I] A; //relatedness matrix
  array[I] vector[T] Z; //mv normal traits
  vector[T] Z_med; //medians for brms-style prior
  vector[T] Z_mad; //MADs for prior
  
}

transformed data{
  //lower Cholesky (improves efficiency)
  matrix[I,I] LA = cholesky_decompose(A); 
  
  //used to specify default brms priors
  vector[T] prior_scale = fmax(rep_vector(2.5, T), Z_mad);
}

parameters {
  //fixed effects
  matrix[F,T] B; //fixed effects
  
  //random effects
  vector<lower=0>[T] Gsd; //genetic sd
  vector<lower=0>[T] Esd; //genetic sd
  
  cholesky_factor_corr[T] LGr; //Cholesky genetic r
  cholesky_factor_corr[T] LEr; //Cholesky environmental r
  matrix[I,T] av; //unscaled breeding values
}

transformed parameters {
  //scaled breeding values (matrix normal method)
  matrix[I,T] a = LA * av * diag_pre_multiply(Gsd, LGr)' ; 
}

model {
  //linear predictor
  matrix[I,T] mu = X * B + a;
  
  //likelihood
  for(i in 1:I){
  Z[i] ~ multi_normal_cholesky(mu[i], //breeding values + fixed effects
          diag_pre_multiply(Esd, LEr)); //Cholesky enviromental cov
  }

  //diffuse priors for demonstration (following brms defaults)
  //weakly regularizing priors recommended in practice
  for (t in 1:T){
    B[1,t] ~ student_t(3, Z_med[t], prior_scale[t]);
    Gsd[t] ~ student_t(3, 0, prior_scale[t]);
    Esd[t] ~ student_t(3, 0, prior_scale[t]);
  }
  to_vector(av) ~ std_normal(); //z-scores (don't change)
  LGr ~ lkj_corr_cholesky(1);
  LEr ~ lkj_corr_cholesky(1); 
}

generated quantities {
  //variances
  vector[T] Gv = Gsd .* Gsd; 
  vector[T] Ev = Esd .* Esd; 
  
  //correlations
  matrix[T,T] Gr = LGr * LGr'; 
  matrix[T,T] Er = LEr * LEr'; 
  
  //covariances S * R * S
  matrix[T,T] G = diag_matrix(Gsd)* Gr * diag_matrix(Gsd); 
  matrix[T,T] E = diag_matrix(Esd)* Er * diag_matrix(Esd); 
  
  //heritability
  vector[T] h2 = Gv ./ (Gv + Ev);
}
