#include <Rcpp.h>
using namespace Rcpp;

// Fortran definitions
void prospect_4(
        double N, 
        double Cab, 
        double Cw,
        double Cm,
        NumericMatrix RT);

void prospect_5(
        double N, 
        double Cab, 
        double Car, 
        double Cw,
        double Cm,
        NumericMatrix RT);

void prospect_5B(
        double N, 
        double Cab, 
        double Car, 
        double Cbrown,
        double Cw,
        double Cm,
        NumericMatrix RT);

