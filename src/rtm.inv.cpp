#include "fortran.h"

NumericVector prospect4_def(NumericVector param){
    double N = param[0], Cab = param[1],
           Cw = param[2], Cm = param[3];
    NumericMatrix RT(2101,2);
    NumericVector R(2101);

    prospect_4_(N, Cab, Cw, Cm, RT);
    R = RT(_,0);
    return R;
}

NumericVector prospect5_def(NumericVector param){
    double N = param[0], Cab = param[1], Car = param[2],
           Cw = param[3], Cm = param[4];
    NumericMatrix RT(2101,2);
    NumericVector R(2101);

    prospect_5_(N, Cab, Car, Cw, Cm, RT);
    R = RT(_,0);
    return R;
}

NumericVector prospect5b_def(NumericVector param){
    double N = param[0], Cab = param[1], Car = param[2],
           Cbrown = param[3], Cw = param[3], Cm = param[4];
    NumericMatrix RT(2101,2);
    NumericVector R(2101);

    prospect_5b_(N, Cab, Car, Cbrown, Cw, Cm, RT);
    R = RT(_,0);
    return R;
}
