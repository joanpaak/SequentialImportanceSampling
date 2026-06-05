/*
N = Number of observations, 
K = Number of categories
*/
data{
  int N;
  int K;
  
  vector[K] prior;
  array[N] int y;
}

parameters{
  simplex[K] phi;
}

model{
  phi ~ dirichlet(prior);
  y ~ categorical(phi);
}
