// nonlin.h
#ifndef NONLIN_H_DEFINED
#define NONLIN_H_DEFINED

#include <stdbool.h>

#ifndef __cplusplus
extern "C" {
#endif

/** @brief Defines a set of solver control information. */
typedef struct {
    /** @brief The maximum number of function evaluations allowed. */
    int max_evals;
    /** @brief The convergence criteria on function values. */
    double fcn_tolerance;
    /** @brief The convergence criteria on change in variable values. */
    double var_tolerances;
    /** @brief The convergence criteria for the slope of the gradient vector. */
    double grad_tolerances;
    /** @brief Controls whether iteration status is printed. */
    bool print_status;
} solver_control;

/** @brief Defines a set of line search controls. */
typedef struct {
    /** @brief The maximum number of function evaluations allowed per search. */
    int max_evals;
    /** @brief Defines the scaling of the product of the gradient and direction 
     * vectors such that F(X + LAMBDA * P) <= 
     * F(X) + LAMBDA * ALPHA * P**T * G, where P is the search direction 
     * vector, G is the gradient vector, and LAMBDA is the scaling factor.  
     * The parameter must exist on the set (0, 1).  A value of 1e-4 is 
     * typically a good starting point.
     */
    double alpha;
    /** @brief Defines a minimum factor X used to determine a minimum value
     * LAMBDA such that MIN(LAMBDA) = X * LAMBDA, where LAMBDA defines the
     * distance along the line search direction assuming a value of one means
     * the full length of the direction vector is traversed.  As such, the value
     * must exist on the set (0, 1); however, for practical considerations,
     * the minimum value should be limited to 0.1 such that the value must
     * exist on the set [0.1, 1).
     */
    double factor;
} line_search_control;

#ifndef __cplusplus
}
#endif
#endif // NONLIN_H_DEFINED