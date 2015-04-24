// Define common functions for later use
#include <Rcpp.h>
using namespace Rcpp;

double exp_int(double k);

// RTMs
NumericVector prospect4_def(NumericVector param);
NumericVector prospect5_def(NumericVector param);
NumericVector prospect5b_def(NumericVector param);

// Priors
double prospect4_priors(int param, double value);
double prospect5_priors(int param, double value);
double prospect5b_priors(int param, double value);


// Truncated normal distribution functions
double rtnorm(double mu, double sd, double MIN);
double dtnorm(double X, double mu, double sd, double MIN);
double rtnorm_c(double mu, double sd, double MIN);
double dtnorm_c(double X, double mu, double sd, double MIN);

// MCMC support functions
NumericMatrix SpecError(NumericVector Model, NumericMatrix Observed);
NumericMatrix SpecError(NumericMatrix Model, NumericMatrix Observed);
double Likelihood(NumericMatrix Error, double rsd);
double Likelihood(NumericVector Error, double rsd);

// RTM selection functions
typedef NumericVector (*select_model)(NumericVector);
typedef double (*select_prior)(int, double);
select_model MODEL(std::string RTM);
select_prior PRIOR(std::string RTM);
NumericVector PMIN(std::string RTM);

// MCMC samplers
void sampler_MH(
        NumericVector &inits,
        double rsd,
        NumericVector &Jump,
        NumericMatrix &Observed,
        NumericMatrix &PrevError,
        NumericVector &ar,
        NumericVector (*Model)(NumericVector),
        double (*Prior)(int, double),
        NumericVector &pmin);
