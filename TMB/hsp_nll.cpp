// simple half-sibling CKMR with exponential pop growth and known ages
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(n_UP);
  DATA_VECTOR(n_HSP);
  DATA_VECTOR(born_year);
  DATA_VECTOR(age_diff);
  PARAMETER(N_init);
  PARAMETER(lambda);
  PARAMETER(phiA);
  ADREPORT(N_init);
  ADREPORT(lambda);
  ADREPORT(phiA);

  Type N;    // for storing the pop size at time of younger kids birth
  Type hsp;  // for storing the half-sibling pair prob

  Type nll = 0;

  for(int i=0;i<n_UP.size(); i++) {
    N = N_init * exp(lambda * born_year[i]);
    hsp = (4 / N) * pow(phiA, age_diff[i]);
    nll -= log(hsp) * n_HSP[i] + log(1 - hsp) * n_UP[i];
  }

  return nll;
}
