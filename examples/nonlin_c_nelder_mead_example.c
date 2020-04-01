// nonlin_c_nelder_mead_example.c

#include <stdio.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))

double fun(int n, const double *x);

int main() {
    // Local Variables
    const int nvar = 2;
    int flag;
    double x[2], f;
    iteration_controls cntrls;
    iteration_process stats;

    // Define the iteration controls - use defaults
    c_set_default_solver_settings(&cntrls);

    // Define an initial guess
    x[0] = x[1] = 0.0;

    // Compute the solution
    flag = c_solver_nelder_mead(fun, nvar, x, &f, &cntrls, &stats);

    // Print the solution
    printf("Solution: (%f, %f)\n", x[0], x[1]);
    printf("Function Value: %e\n", f);
    printf("Iterations: %i\n", stats.iteration_count);
    printf("Function Evaluations: %i\n", stats.function_eval_count);

    // End
    return 0;
}

// Rosenbrock's Function
// f(x) = 1e2 * (x2 - x1^2)^2 + (x1 - 1)^2
double fun(int n, const double *x) {
    return 1.0e2 * SQR(x[1] - SQR(x[0])) + SQR(x[0] - 1.0);
}
