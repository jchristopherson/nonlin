// nonlin_c_test.c

#include <stdio.h>
#include "nonlin_c_test.h"

int main() {
    // Local Variables
    bool rst, overall = true;

    // Tests
    rst = test_quasinewton();
    if (!rst) overall = false;

    rst = test_newton();
    if (!rst) overall = false;

    rst = test_least_squares();
    if (!rst) overall = false;

    rst = test_nelder_mead();
    if (!rst) overall = false;

    rst = test_bfgs();
    if (!rst) overall = false;

    rst = test_poly_roots();
    if (!rst) overall = false;

    rst = test_poly_multiply();
    if (!rst) overall = false;

    // End
    if (overall) printf("NONLIN C TEST STATUS: PASS\n");
    else printf("NONLIN C TEST STATUS: FAILED\n");
    return 0;
}