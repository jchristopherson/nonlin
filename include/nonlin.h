// nonlin.h
#ifndef NONLIN_H_DEFINED
#define NONLIN_H_DEFINED

#include <stdbool.h>
#include "external/linalg/ferror/include/ferror.h"

/** @brief Defines the signature of a function of one variable.
 *
 * @param x The independent variable.
 *
 * @return The value of the function at @p x.
 */
typedef double (*fcn1var)(double x);

/** @brief Describes an M-element vector-valued function of N-variables.
 *
 * @param[in] x An N-element array containing the independent variables.
 * @param[out] f An M-element array that, on output, contains the values
 *  of the M functions.
 */
typedef void (*vecfcn)(const double *x, double *f);

/** @brief Describes a routine capable of computing the Jacobian matrix
 * of M functions of N unknowns.
 *
 * @param[in] x An N-element array containing the independent variables.
 * @param[out] jac An M-by-N matrix where the Jacobian will be written.
 */
typedef void (*jacobianfcn)(const double *x, double *jac);


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

/** @brief Defines a pair of numeric values. */
typedef struct {
    /** @brief Value 1. */
    double x1;
    /** @brief Value 2. */
    double x2;
} value_pair;

/** @brief Defines a set of parameters that describe the behavior of the
 * iteration process.
 */
typedef struct {
    /** @brief Specifies the number of iterations performed. */
    int iter_count;
    /** @brief Specifies the number of function evaluations performed. */
    int fcn_count;
    /** @brief Specifies the number of Jacobian evaluations performed. */
    int jacobian_count;
    /** @brief True if the solution converged as a result of a zero-valued
     * function; else, false.
     */
    bool converge_on_fcn;
    /** @brief True if the solution converged as a result of no appreciable
     * change in solution points between iterations; else, false.
     */
    bool converg_on_chng;
    /** @brief True if the solution appears to have settled on a stationary
     * point such that the gradient of the function is zero-valued; else,
     * false.
     */
    bool converge_on_zero_diff;
} iteration_behavior;


#ifdef __cplusplus
extern "C" {
#endif

/** @brief Solves an equation of one variable using Brent's method.
 *
 * @param[in] fcn A pointer to the routine containing the function to solve.
 * @param[in] lim A value_pair object defining the search limits.
 * @param[out] x On output, the solution.
 * @param[out] f On output, the residual as computed at @p x.
 * @param[in] tol A solver_control object defining the solver control
 *  parameters.
 * @param[out] ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param[in] err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 *  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
 *      than the number of variables.
 *  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
 *      the allowed number of iterations.
 */
void solve_brent(fcn1var fcn, value_pair lim, double *x, double *f,
                 const solver_control *tol, iteration_behavior *ib,
                 errorhandler err);

/** @brief Applies the quasi-Newton's method developed by Broyden in
 * conjunction with a backtracking type line search to solve N equations
 * of N unknowns.
 *
 * @param[in] fcn A pointer to the routine containing the system of
 *  equations to solve.
 * @param[in] jac A pointer to a routine used to compute the Jacobian of
 *  the system of equations.  To let the program compute the Jacobian
 *  numerically, simply pass NULL.
 * @param[in] n The number of equations, and the number of unknowns.
 * @param[in,out] x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 *  N is the number of variables.
 * @param[out] fvec An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param[in] tol A solver_control object defining the solver control
 *  parameters.
 * @param[in] lsearch A pointer to a line_search_control object defining
 *  the line search control parameters.  If no line search is desired,
 *  simply pass NULL.
 * @param[out] ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param[in] err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 *  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
 *      than the number of variables.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      correctly.
 *  - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
 *      pointing in an apparent uphill direction.
 *  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
 *      the allowed number of iterations.
 *  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
 *      gradient vector becomes sufficiently close to zero.
 */
void solve_quasi_newton(vecfcn fcn, jacobianfcn jac, int n, double *x,
                        double *fvec, const solver_control *tol,
                        line_search_control *lsearch,
                        iteration_behavior *ib, errorhandler err);

#ifdef __cplusplus
}
#endif
#endif // NONLIN_H_DEFINED
