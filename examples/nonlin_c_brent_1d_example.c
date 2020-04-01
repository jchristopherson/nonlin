// nonlin_c_brent_1d_example.c

#include <stdio.h>
#include "nonlin.h"

double eqn(double x);

int main() {
    // Local Variables
    int flag;
    double xo, fo;
    value_pair limits;
    iteration_controls cntrls;
    iteration_process stats;

    // Define the allowable solution domain
    limits.x1 = -2.0;
    limits.x2 = 2.0;

    // Define an initial guess
    xo = 0.0;

    // Define the iteration controls - use defaults
    c_set_default_solver_settings(&cntrls);

    // Compute the solution
    flag = c_solver_brent_1var(eqn, &xo, &fo, limits, &cntrls, &stats);

    // Print the solution
    printf("Solution: %f\nFunction Value: %f\n", xo, fo);

    // Print the iteration information
    printf("Iteration Count: %i\n", stats.iteration_count);
    printf("Function Evaluation Count: %i\n", stats.function_eval_count);
    printf("Converge on Function Value: %s\n", 
        stats.converge_on_function ? "true" : "false");
    printf("Converge on Change in Varible: %s\n",
        stats.converge_on_solution_change ? "true" : "false");
    printf("Converge on Derivative: %s\n",
        stats.converge_on_gradient ? "true" : "false");

    // End
    return 0;
}

// The equation to solve:
// f(x) = x^3 - 2 * x - 1
//
// Solution:
// x = -1, -0.618034, 1.618034
double eqn(double x) {
    return (x * x * x) - 2.0 * x - 1.0;
}
