#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP logLik_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP score_rho_dontknow(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP hess_rho_dontknow(SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"logLik_dontknow", (DL_FUNC) &logLik_dontknow, 7},
  {"score_rho_dontknow", (DL_FUNC) &score_rho_dontknow, 5},
  {"hess_rho_dontknow",  (DL_FUNC) &hess_rho_dontknow,  4},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

