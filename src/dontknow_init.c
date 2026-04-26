#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP logLik_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP score_rho_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP hess_rho_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP score_mu1_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP score_mu2_dontknow(SEXP, SEXP, SEXP,SEXP, SEXP);
SEXP hess_mu1_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP hess_mu2_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP score_alpha_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP hess_alpha_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP z_weights_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"logLik_dontknow", (DL_FUNC) &logLik_dontknow, 7},
  {"score_rho_dontknow", (DL_FUNC) &score_rho_dontknow, 5},
  {"hess_rho_dontknow",  (DL_FUNC) &hess_rho_dontknow,  5},
  {"score_mu1_dontknow", (DL_FUNC) &score_mu1_dontknow, 5},
  {"score_mu2_dontknow", (DL_FUNC) &score_mu2_dontknow, 5},
  {"hess_mu1_dontknow",  (DL_FUNC) &hess_mu1_dontknow,  5},
  {"hess_mu2_dontknow",  (DL_FUNC) &hess_mu2_dontknow,  5},
  {"score_alpha_dontknow",(DL_FUNC) &score_alpha_dontknow,6},
  {"hess_alpha_dontknow", (DL_FUNC) &hess_alpha_dontknow, 6},
  {"z_weights_dontknow",  (DL_FUNC) &z_weights_dontknow,  7},
  {NULL, NULL, 0}
};

void R_init_dontknow(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
