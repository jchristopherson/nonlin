// nonlin_c_test_solve.c

#include "nonlin_c_test.h"

// System of equations #1
//
// x^2 + y^2 = 34
// x^2 - 2 * y^2 = 7;
void fcn1(int neqn, int nvar, const double *x, double *f) {
    f[0] = SQR(x[0]) + SQR(x[1]) - 34.0;
    f[1] = SQR(x[0]) - 2.0 * SQR(x[1]) - 7.0;
}

// Jacobian for system #1
//
//     | 2x  2y |
// J = |        |
//     | 2x  -4y|
void jac1(int neqn, int nvar, const double *x, double *jac) {
    jac[0] = 2.0 * x[0];
    jac[1] = 2.0 * x[0];
    jac[2] = 2.0 * x[1];
    jac[3] = -4.0 * x[1];
}

// Tests against the answer for system #1
bool is_ans_1(const double *x, double tol) {
    bool rst = true;
    const double x1 = 5.0;
    const double x2 = 3.0;
    double ax1, ax2;
    ax1 = abs(x[0]) - x1;
    ax2 = abs(x[1]) - x2;
    printf("AX1 = %f\nAX2 = %f\nX[0] = %f\nX[1] = %f\nX1 = %f\nX2 = %f\n", ax1, ax2, x[0], x[1], x1, x2);
    if (abs(ax1) > tol || abs(ax2) > tol) rst = false;
    return rst;
}



bool test_quasinewton() {
    // Local Variables
    const double test = 1.0e-6;
    const double ans1 = 5.0;
    const double ans2 = 3.0;

    bool rst = true;
    iteration_behavior ib;
    solver_control tol;
    double x[2], f[2];

    // Set up tolerances
    tol.max_evals = 500;
    tol.fcn_tolerance = 1.0e-8;
    tol.var_tolerances = 1.0e-12;
    tol.grad_tolerances = 1.0e-12;
    tol.print_status = false;

    // Define an initial guess
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution
    solve_quasi_newton(fcn1, NULL, 2, x, f, &tol, NULL, &ib, NULL);

    // Test
    if (!is_ans_1(x, test)) {
        rst = false;
        printf("Test Failed: Quasi-Newton, Sytem #1\nExpected: +/-(%f, %f)\nReceived: (%f, %f)\n",
            ans1, ans2, x[0], x[1]);
    }

    // End
    return rst;
}