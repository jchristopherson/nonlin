// nonlin_c_test_poly.c

#include "nonlin_c_test.h"

bool test_poly_roots() {
    // Local Variables
    const double tol = 1.0e-8;
    const int order = 10;

    int i;
    polynomial p;
    double coeff[order + 1];
    double complex rts[order], sol[order];
    bool rst = true;

    // Initialization
    make_rand_mtx(order + 1, 1, coeff);
    alloc_polynomial(&p, order);
    for (i = 0; i <= order; ++i) 
        set_polynomial_coefficient(&p, i + 1, coeff[i], NULL);
    
    // Compute the roots
    polynomial_roots(&p, order, rts, NULL);

    // Compute the value of the polynomial at each root, and ensure it is 
    // sufficiently close to zero
    evaluate_polynomial_cmplx(&p, order, rts, sol);
    for (i = 0; i < order; ++i) {
        if (cabs(sol[i]) > tol) {
            rst = false;
            printf("Test Failed: Polynomial Roots Test\n");
            break;
        }
    }

    // End
    free_polynomial(&p);
    return rst;
}

bool test_poly_multiply() {
    // Local Variables
    const double tol = 1.0e-8;
    bool rst = true;

    int i;
    polynomial x, y, z;
    double sol[6], ans[6] = {5.0, 10.0, 30.0, 26.0, 52.0, 24.0};

    // Initialization
    alloc_polynomial(&x, 3);
    alloc_polynomial(&y, 2);
    alloc_polynomial(&z, 5);

    set_polynomial_coefficient(&x, 1, 5.0, NULL);
    set_polynomial_coefficient(&x, 2, 0.0, NULL);
    set_polynomial_coefficient(&x, 3, 10.0, NULL);
    set_polynomial_coefficient(&x, 4, 6.0, NULL);

    set_polynomial_coefficient(&y, 1, 1.0, NULL);
    set_polynomial_coefficient(&y, 2, 2.0, NULL);
    set_polynomial_coefficient(&y, 3, 4.0, NULL);

    // Compute the product of the two polynomials
    polynomial_multiply(&x, &y, &z);

    // Test
    for (i = 0; i < 6; ++i)
        sol[i] = get_polynomial_coefficient(&z, i + 1, NULL);
    if (!is_dbl_mtx_equal(6, 1, sol, ans, tol)) {
        rst = false;
        printf("Test Failed: Polynomial Multiplication:\n");
    }


    // End
    free_polynomial(&x);
    free_polynomial(&y);
    free_polynomial(&z);
    return rst;
}
