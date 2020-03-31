// nonlin_c_polynomial_example.c

#include <stdio.h>
#include "nonlin.h"

int main() {
    // Local Variables
    const int neqn = 21;
    const int norder = 3;
    c_polynomial poly;
    int i, flag;
    double coeffs[4], y[21];
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

    // Create a new third order polynomial
    flag = c_init_polynomial(norder, &poly);
    if (flag != NL_NO_ERROR) printf("Error Code: %i\n", flag);

    // Check the polynomial order
    flag = c_get_polynomial_order(&poly);
    printf("Polynomial Order: %i\n", flag);

    // Fit data to the polynomial
    flag = c_fit_polynomial(&poly, neqn, xp, yp, false);

    // Print out the coefficients
    flag = c_get_polynomial_coefficients(&poly, norder + 1, coeffs);
    for (i = 0; i <= norder; ++i) printf("c[%i] = %f\n", i, coeffs[i]);

    // Evaluate the fitted polynomial and compare to the actual values
    c_evaluate_polynomial_real(&poly, neqn, xp, y);
    printf("\n");
    for (i = 0; i < neqn; ++i) 
    printf("Actual: %f\tFitted: %f\n", yp[i], y[i]);

    // End
    c_free_polynomial(&poly);
    return 0;
}