#include <Rcpp.h>
using namespace Rcpp;

// Fortran definitions
void prospect_4_(
        double &N, 
        double &Cab, 
        double &Cw,
        double &Cm,
        NumericMatrix &RT);

void prospect_5_(
        double &N, 
        double &Cab, 
        double &Car, 
        double &Cw,
        double &Cm,
        NumericMatrix &RT);

void prospect_5b_(
        double &N, 
        double &Cab, 
        double &Car, 
        double &Cbrown,
        double &Cw,
        double &Cm,
        NumericMatrix &RT);

