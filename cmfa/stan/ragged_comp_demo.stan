/* A basic model implementing compositional regression on ragged data.

  The ragged data comes in two containers:

    - the vector stacked_y contains all the values in one long vector.
    - the array y_sizes contains the size of each element of the ragged array.

  Transformation from and to simplexes is handled by the centered log ratio
  transformation, implemented by functions clr and clr_inv in the functions
  block. See the python package compositional_stats for reference implentations
  in Python.

*/

functions {
  real geometric_mean(vector x){
   return exp(inv(size(x) * sum(log(x))));
  }
  vector clr(vector x){
   return log(x / geometric_mean(x));
  }
  vector clr_inv(vector x){
   vector[size(x)] expx = exp(x);
   return expx / sum(expx);
  }
  tuple(int, int) get_ragged_start_and_end(int pos, array[] int sizes){
    tuple(int, int) out;
    out.1 = pos == 1 ? 1 : 1 + sum(sizes[1:pos-1]);
    out.2 = out.1 + sizes[pos] - 1;
   return out;
  }
  vector extract_ragged(int pos, vector long, array[] int sizes){
   int start, end;
   (start, end) = get_ragged_start_and_end(pos, sizes);
   return long[start: end];
  }
}
data {
 int<lower=1> N;
 int<lower=1> N_measurement;
 array[N_measurement] int<lower=1> y_sizes;
 vector[N] stacked_y;
}
parameters {
 vector[N] stacked_yhat_clr;
 real<lower=0> sigma;
}
model {
 sigma ~ normal(0, 1);
 for (n in 1:N_measurement){
   int k = y_sizes[n];
   vector[k] y = extract_ragged(n, stacked_y, y_sizes);
   vector[k] yhat_clr = extract_ragged(n, stacked_yhat_clr, y_sizes);
   clr(y) ~ normal(yhat_clr, sigma);
 }
}
generated quantities {
 vector[N] stacked_yhat;
 vector[N] stacked_yrep;
 for (n in 1:N_measurement){
  int start, end;
  (start, end) = get_ragged_start_and_end(n, y_sizes);
  int k = y_sizes[n];
  vector[k] yhat_clr = stacked_yhat_clr[start: end];
  vector[k] yrep_clr = to_vector(normal_rng(yhat_clr, sigma));
  stacked_yhat[start: end] = clr_inv(yhat_clr);
  stacked_yrep[start: end] = clr_inv(yhat_clr);
 }
}
