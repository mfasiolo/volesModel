// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// volesFullCpp
Rcpp::List volesFullCpp(SEXP nMon_, SEXP nSimul_, SEXP nBurn_, SEXP params_, SEXP nSteps_, SEXP T0_, SEXP randInit_, SEXP startVole_, SEXP startWeas_, SEXP addObsNoise_);
RcppExport SEXP volesModel_volesFullCpp(SEXP nMon_SEXP, SEXP nSimul_SEXP, SEXP nBurn_SEXP, SEXP params_SEXP, SEXP nSteps_SEXP, SEXP T0_SEXP, SEXP randInit_SEXP, SEXP startVole_SEXP, SEXP startWeas_SEXP, SEXP addObsNoise_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type nMon_(nMon_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nSimul_(nSimul_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nBurn_(nBurn_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type params_(params_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nSteps_(nSteps_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type T0_(T0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type randInit_(randInit_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type startVole_(startVole_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type startWeas_(startWeas_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type addObsNoise_(addObsNoise_SEXP);
    __result = Rcpp::wrap(volesFullCpp(nMon_, nSimul_, nBurn_, params_, nSteps_, T0_, randInit_, startVole_, startWeas_, addObsNoise_));
    return __result;
END_RCPP
}
// volesStdCpp
Rcpp::List volesStdCpp(SEXP nMon_, SEXP nSimul_, SEXP nBurn_, SEXP params_, SEXP nSteps_, SEXP T0_, SEXP randInit_, SEXP startVole_, SEXP startWeas_, SEXP addObsNoise_);
RcppExport SEXP volesModel_volesStdCpp(SEXP nMon_SEXP, SEXP nSimul_SEXP, SEXP nBurn_SEXP, SEXP params_SEXP, SEXP nSteps_SEXP, SEXP T0_SEXP, SEXP randInit_SEXP, SEXP startVole_SEXP, SEXP startWeas_SEXP, SEXP addObsNoise_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type nMon_(nMon_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nSimul_(nSimul_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nBurn_(nBurn_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type params_(params_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nSteps_(nSteps_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type T0_(T0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type randInit_(randInit_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type startVole_(startVole_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type startWeas_(startWeas_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type addObsNoise_(addObsNoise_SEXP);
    __result = Rcpp::wrap(volesStdCpp(nMon_, nSimul_, nBurn_, params_, nSteps_, T0_, randInit_, startVole_, startWeas_, addObsNoise_));
    return __result;
END_RCPP
}
