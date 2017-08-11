// nonlin_bfgs_c_example.c

#include <stdio.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))
double fcn(int nvar, const double *x);

int main() {
    // Local Variables
    iteration_behavior ib;
    solver_control tol;
    line_search_control ls;
    double f, x[2] = {0.0, 0.0};

    // Define the default tolerances
    set_nonlin_defaults(&tol);
    set_nonlin_ls_defaults(&ls);

    // Compute the solution.  Let the solver compute the gradient vector
    // via numerical means.
    bfgs(fcn, NULL, 2, x, &f, &tol, &ls, &ib, NULL);

    // Display the results
    printf("Solution: (%6.4f, %6.4f)\nFunction Value: %6.4f\nIterations: %i\n",
        x[0], x[1], f, ib.iter_count);
}

// Rosenbrock's Function:
double fcn(int nvar, const double *x) {
    return 1.0e2 * SQR(x[1] - SQR(x[0])) + SQR(x[0] - 1.0);
}
