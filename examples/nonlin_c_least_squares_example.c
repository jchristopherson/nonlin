// nonlin_c_least_squares_example.c

#include <stdlib.h>
#include <stdio.h>
#include "nonlin.h"

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))

void fun(int neqn, int nvar, const double *x, double *f);

int main() {
    // Local Variables
    const int neqn = 21;
    const int nvar = 4;
    int flag;
    double x[4], f[21];
    iteration_controls cntrls;
    iteration_process stats;

    // Defin the iteration controls - use defaults
    c_set_default_solver_settings(&cntrls);

    // Define an initial guess
    x[0] = x[1] = x[2] = x[3] = 1.0;

    // Compute the solution
    flag = c_solver_least_squares(fun, NULL, neqn, nvar, x, f, &cntrls, &stats);

    // Print out the solution
    printf("Solution: (%f, %f, %f, %f)\n", x[0], x[1], x[2], x[3]);
    printf("Residual: (%f, %f, %f, %f)\n", f[0], f[1], f[2], f[3]);
    printf("Iterations: %i\n", stats.iteration_count);
    printf("Function Evaluations: %i\n", stats.function_eval_count);
    printf("Jacobian Evaluations: %i\n", stats.jacobian_eval_count);

    // End
    return 0;
}

void fun(int neqn, int nvar, const double *x, double *f) {
    // Solution Points
    double xp[] = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
        0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
        1.8, 1.9, 2.0
    };
    double yp[] = {
        1.216737514, 1.250032542, 1.305579195, 1.040182335, 
        1.751867738, 1.109716707, 2.018141531, 1.992418729, 
        1.807916923, 2.078806005, 2.698801324, 2.644662712, 
        3.412756702, 4.406137221, 4.567156645, 4.999550779, 
        5.652854194, 6.784320119, 8.307936836, 8.395126494, 
        10.30252404
    };
    int i;

    // Apply the cubic polynomial to each equation
    for (i = 0; i < neqn; ++i) {
        f[i] = x[0] * CUBE(xp[i]) + 
            x[1] * SQR(xp[i]) + 
            x[2] * xp[i] + 
            x[3] - yp[i];
    }
}
