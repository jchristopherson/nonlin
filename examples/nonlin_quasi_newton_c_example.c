// nonlin_quasi_newton_c_example.c

#include <stdio.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))
void fcn(int neqn, int nvar, const double *x, double *f);

int main() {
    // Local Variables
    solver_control tol;
    line_search_control ls;
    iteration_behavior ib;
    double f[2], x[2] = {1.0, 1.0};

    // Define tolerances - use default values
    set_nonlin_defaults(&tol);
    set_nonlin_ls_defaults(&ls);

    // Compute the solution
    solve_quasi_newton(fcn, NULL, 2, x, f, &tol, &ls, &ib, NULL);

    // Display the results
    printf("Solution: (%f, %f)\n", x[0], x[1]);
    printf("Residual: (%e, %e)\n", f[0], f[1]);
    printf("Iterations: %i\nFunction Evaluations: %i\nJacobian Evaluations: %i\n",
        ib.iter_count, ib.fcn_count, ib.jacobian_count);

    // End
    return 0;
}

void fcn(int neqn, int nvar, const double *x, double *f) {
    // x^2 + y^2 = 34
    // x^2 - 2 * y^2 = 7
    f[0] = SQR(x[0]) + SQR(x[1]) - 34.0;
    f[1] = SQR(x[0]) - 2.0 * SQR(x[1]) - 7.0;
}
