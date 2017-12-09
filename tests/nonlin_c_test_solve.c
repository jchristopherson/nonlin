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
    double ax1, ax2;
    ax1 = fabs(x[0]) - 5.0;
    ax2 = fabs(x[1]) - 3.0;
    if (abs(ax1) > tol || abs(ax2) > tol) rst = false;
    return rst;
}

double rosenbrock(int n, const double *x) {
    return 1.0e2 * SQR(x[1] - SQR(x[0])) + SQR(x[0] - 1.0);
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
    tol.var_tolerance = 1.0e-12;
    tol.grad_tolerance = 1.0e-12;
    tol.print_status = false;

    // Define an initial guess
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution
    solve_quasi_newton(fcn1, jac1, 2, x, f, &tol, NULL, &ib, NULL);

    // Test
    if (!is_ans_1(x, test)) {
        rst = false;
        printf("Test Failed: Quasi-Newton, Sytem #1\nExpected: +/-(%f, %f)\nReceived: (%f, %f)\n",
            ans1, ans2, x[0], x[1]);
    }

    // End
    return rst;
}



bool test_newton() {
    // Local Variables
    const double test = 1.0e-6;
    const double ans1 = 5.0;
    const double ans2 = 3.0;

    bool rst = true;
    iteration_behavior ib;
    solver_control tol;
    line_search_control ls;
    double x[2], f[2];

    // Set up tolerances
    set_nonlin_defaults(&tol);
    set_nonlin_ls_defaults(&ls);

    // Define an initial guess
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution
    solve_newton(fcn1, jac1, 2, x, f, &tol, &ls, &ib, NULL);

    // Test
    if (!is_ans_1(x, test)) {
        rst = false;
        printf("Test Failed: Newton, Sytem #1\nExpected: +/-(%f, %f)\nReceived: (%f, %f)\n",
            ans1, ans2, x[0], x[1]);
    }

    // End
    return rst;
}



bool test_least_squares() {
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
    tol.var_tolerance = 1.0e-12;
    tol.grad_tolerance = 1.0e-12;
    tol.print_status = false;

    // Define an initial guess
    x[0] = 1.0;
    x[1] = 1.0;

    // Compute the solution
    solve_nl_least_squares(fcn1, jac1, 2, 2, x, f, &tol, &ib, NULL);

    // Test
    if (!is_ans_1(x, test)) {
        rst = false;
        printf("Test Failed: Levenberg-Marquardt, Sytem #1\nExpected: +/-(%f, %f)\nReceived: (%f, %f)\n",
            ans1, ans2, x[0], x[1]);
    }

    // End
    return rst;
}


bool test_nelder_mead() {
    // Local Variables
    const double check = 1.0e-5;
    bool rst = true;
    iteration_behavior ib;
    solver_control tol;
    double f, s[6], x[2] = {0.0, 0.0};

    // Initialization
    set_nonlin_defaults(&tol);
    tol.max_evals = 500; // This routine can take a lot of evaluations to converge
    tol.fcn_tolerance = 1.0e-12;
    
    // Define a simplex (this is an optional step, but can help convergence)
    s[0] = 0.0;     s[2] = 0.0;     s[4] = 2.0;
    s[1] = 0.0;     s[3] = 2.0;     s[5] = 0.0;

    // Compute the solution
    nelder_mead(rosenbrock, 2, x, &f, s, &tol, &ib, NULL);

    // Test (solution is [1, 1])
    if (fabs(x[0] - 1.0) > check || fabs(x[1] - 1.0) > check) {
        rst = false;
        printf("Test Failed: Nelder-Mead, Rosebrock Function\nExpected: (1, 1)\nReceived: (%f, %f)\n",
            x[0], x[1]);
    }

    // End
    return rst;
}


bool test_bfgs() {
    // Local Variables
    const double check = 1.0e-5;
    bool rst = true;
    iteration_behavior ib;
    solver_control tol;
    line_search_control ls;
    double f, x[2] = {0.0, 0.0};

    // Initialization
    set_nonlin_defaults(&tol);
    set_nonlin_ls_defaults(&ls);

    // Compute the solution
    bfgs(rosenbrock, NULL, 2, x, &f, &tol, &ls, &ib, NULL);

    // Test (solution is [1, 1])
    if (fabs(x[0] - 1.0) > check || fabs(x[1] - 1.0) > check) {
        rst = false;
        printf("Test Failed: BFGS, Rosebrock Function\nExpected: (1, 1)\nReceived: (%f, %f)\n",
            x[0], x[1]);
    }

    // End
    return rst;
}
