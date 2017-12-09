// nonlin_c_test.h
#ifndef NONLIN_C_TEST_INCLUDED_
#define NONLIN_C_TEST_INCLUDED_

#include "nonlin.h"
#include "c_test_core.h"
#include <math.h>

// Macros
#define SQR(x) ((x) * (x))

// Routines containing functions to solve
void fcn1(int neqn, int nvar, const double *x, double *f);
void jac1(int neqn, int nvar, const double *x, double *jac);

// Solver Test Routines
bool test_quasinewton();
bool test_newton();
bool test_least_squares();
bool test_nelder_mead();
bool test_bfgs();

// Polynomial Test Routines
bool test_poly_roots();
bool test_poly_multiply();

#endif
