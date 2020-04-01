// nonlin_c_bfgs_example.c

#include <stdlib.h>
#include <stdio.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))

double fun(int n, const double *x);
void grad(int n, const double *x, double *g);

int main() {
    // Local Variables
    const int nvar = 2;
    int flag;
    double x[2], f;
    iteration_controls cntrls;
    line_search_controls linesearch;
    iteration_process stats;

    // Define the iteration controls - use defaults
    c_set_default_solver_settings(&cntrls);
    c_set_default_line_search_settings(&linesearch);

    // Define an initial guess
    x[0] = x[1] = 0.0;

    // Compute the solution - analytical gradient
    flag = c_solver_bfgs(fun, grad, nvar, x, &f, &cntrls, &linesearch, &stats);

    // Print the solution
    printf("\nAnalytical Gradient:\n");
    printf("Solution: (%f, %f)\n", x[0], x[1]);
    printf("Function Value: %e\n", f);
    printf("Iterations: %i\n", stats.iteration_count);
    printf("Function Evaluations: %i\n", stats.function_eval_count);

    // Reset the initial conditions
    x[0] = x[1] = 0.0;

    // Compute the solution - numerical gradient
    flag = c_solver_bfgs(fun, NULL, nvar, x, &f, &cntrls, &linesearch, &stats);

    // Print the solution
    printf("\nNumerical Gradient:\n");
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

// Gradient Vector
void grad(int n, const double *x, double *g) {
    g[0] = 2.0 * (x[0] - 1.0) - 4.0e2 * x[0] * (x[1] - SQR(x[0]));
    g[1] = 2.0e2 * (x[1] - SQR(x[0]));
}
