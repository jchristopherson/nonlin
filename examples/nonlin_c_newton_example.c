// nonlin_c_newton_example.c

#include <stdio.h>
#include <stdlib.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))
#define INDEX(i, j, m) ((m) * (j) + (i))

void eqns(int neqn, int nvar, const double *x, double *f);
void jacobian(int neqn, int nvar, const double *x, double *jac);

int main() {
    // Local Variables
    const int n = 2;
    int flag;
    double x[2], f[2];
    iteration_controls cntrls;
    line_search_controls linesearch;
    iteration_process stats;

    // Defin the iteration and line search controls - use defaults
    c_set_default_solver_settings(&cntrls);
    c_set_default_line_search_settings(&linesearch);

    // Define initial conditions
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution - analytical Jacobian
    flag = c_solver_newton(eqns, jacobian, n, x, f, &cntrls, 
        &linesearch, &stats);
    printf("Analytical Jacobian Results:\n");
    printf("Solution: (%f, %f)\n", x[0], x[1]);
    printf("Residual: (%e, %e)\n", f[0], f[1]);
    printf("Iterations: %i\n", stats.iteration_count);
    printf("Jacobian Evaluations: %i\n", stats.jacobian_eval_count);

    // Reset the initial conditions
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution - numerical Jacobian
    flag = c_solver_newton(eqns, NULL, n, x, f, &cntrls, 
        &linesearch, &stats);
    printf("\nNumerical Jacobian Results:\n");
    printf("Solution: (%f, %f)\n", x[0], x[1]);
    printf("Residual: (%e, %e)\n", f[0], f[1]);
    printf("Iterations: %i\n", stats.iteration_count);
    printf("Jacobian Evaluations: %i\n", stats.jacobian_eval_count);

    // End
    return 0;
}

// Equations:
// x^2 + y^2 = 34
// x^2 - 2 * y^2 = 7
void eqns(int neqn, int nvar, const double *x, double *f) {
    f[0] = SQR(x[0]) + SQR(x[1]) - 34.0;
    f[1] = SQR(x[0]) - 2.0 * SQR(x[1]) - 7.0;
}

// Jacobian:
// | df1/dx  df1/dy |   | 2*x   2*y |
// |                | = |           |
// | df2/dx  df2/dy |   | 2*x  -4*y |
void jacobian(int neqn, int nvar, const double *x, double *jac) {
    jac[INDEX(0,0,neqn)] = 2.0 * x[0];      jac[INDEX(0,1,neqn)] = 2.0 * x[1];
    jac[INDEX(1,0,neqn)] = 2.0 * x[0];      jac[INDEX(1,1,neqn)] = -4.0 * x[1];
}
