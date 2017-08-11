// nonlin.h
#ifndef NONLIN_H_DEFINED
#define NONLIN_H_DEFINED

#include <stdbool.h>
#include <complex.h>
#include "ferror.h"

/** An error flag denoting an invalid input. */
#define NL_INVALID_INPUT_ERROR = 201
/** An error flag denoting an improperly sized array. */
#define NL_ARRAY_SIZE_ERROR = 202
/** An error denoting that there is insufficient memory available. */
#define NL_OUT_OF_MEMORY_ERROR = 203
/** An error resulting from an invalid operation. */
#define NL_INVALID_OPERATION_ERROR = 204
/** An error resulting from a lack of convergence. */
#define NL_CONVERGENCE_ERROR = 205
/** An error resulting from a divergent condition. */
#define NL_DIVERGENT_BEHAVIOR_ERROR = 206
/** An error indicating a possible spurious convergence condition. */
#define NL_SPURIOUS_CONVERGENCE_ERROR = 207
/** An error indicating the user-requested tolerance is too small to be
 * practical for the problem at hand. */
#define NL_TOLERANCE_TOO_SMALL_ERROR = 208

/** @brief Defines the signature of a function of one variable.
 *
 * @param x The independent variable.
 *
 * @return The value of the function at @p x.
 */
typedef double (*fcn1var)(double x);

/** @brief Describes an M-element vector-valued function of N-variables.
 *
 * @param neqn The number of equations.
 * @param nvar The number of variables.
 * @param x An NVAR-element array containing the independent variables.
 * @param f An NEQN-element array that, on output, contains the values
 *  of the M functions.
 */
typedef void (*vecfcn)(int neqn, int nvar, const double *x, double *f);

/** @brief Describes a routine capable of computing the Jacobian matrix
 * of M functions of N unknowns.
 *
 * @param neqn The number of equations.
 * @param nvar The number of variables.
 * @param x An NVAR-element array containing the independent variables.
 * @param jac An NEQN-by-NVAR matrix where the Jacobian will be written.
 */
typedef void (*jacobianfcn)(int neqn, int nvar, const double *x, double *jac);

/** @brief Describes a function of N variables.
!!
!! @param[in] nvar The number of variables.
!! @param[in] x An NVAR-element array containing the independent variables.
!! @return The value of the function at @p x.
 */
typedef double (*fcnnvar)(int nvar, const double *x);

/** @brief Describes a routine capable of computing the gradient vector
!! of an equation of N variables.
!!
!! @param[in] nvar The number of variables.
!! @param[in] x An NVAR-element array containing the independent variables.
!! @param[out] g An NVAR-element array where the gradient vector will be
!!  written as output.
 */
typedef void (*gradientfcn)(int nvar, const double *x, double *g);


/** @brief Defines a set of solver control information. */
typedef struct {
    /** @brief The maximum number of function evaluations allowed. */
    int max_evals;
    /** @brief The convergence criteria on function values. */
    double fcn_tolerance;
    /** @brief The convergence criteria on change in variable values. */
    double var_tolerance;
    /** @brief The convergence criteria for the slope of the gradient vector. */
    double grad_tolerance;
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
    /** Specifies the number of gradient vector evaluations performed. */
    int gradient_count;
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

/** @brief A C compatible type encapsulating a polynomial object. */
typedef struct {
    /** @brief A pointer to the polynomial object. */
    void *ptr;
    /** @brief The size of the polynomial object, in bytes. */
    int n;
} polynomial;


#ifdef __cplusplus
extern "C" {
#endif

/** @brief Solves an equation of one variable using Brent's method.
 *
 * @param fcn A pointer to the routine containing the function to solve.
 * @param lim A value_pair object defining the search limits.
 * @param x On output, the solution.
 * @param f On output, the residual as computed at @p x.
 * @param tol A solver_control object defining the solver control
 *  parameters.
 * @param ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param err A pointer to the C error handler object.  If no error
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
 * @param fcn A pointer to the routine containing the system of
 *  equations to solve.
 * @param jac A pointer to a routine used to compute the Jacobian of
 *  the system of equations.  To let the program compute the Jacobian
 *  numerically, simply pass NULL.
 * @param n The number of equations, and the number of unknowns.
 * @param x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 *  N is the number of variables.
 * @param fvec An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param tol A solver_control object defining the solver control
 *  parameters.
 * @param lsearch A pointer to a line_search_control object defining
 *  the line search control parameters.  If no line search is desired,
 *  simply pass NULL.
 * @param ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param err A pointer to the C error handler object.  If no error
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
                        const line_search_control *lsearch,
                        iteration_behavior *ib, errorhandler *err);

/** @brief Applies Newton's method in conjunction with a backtracking type
 * line search to solve N equations of N unknowns.
 *
 * @param fcn A pointer to the routine containing the system of
 *  equations to solve.
 * @param jac A pointer to a routine used to compute the Jacobian of
 *  the system of equations.  To let the program compute the Jacobian
 *  numerically, simply pass NULL.
 * @param n The number of equations, and the number of unknowns.
 * @param x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 *  N is the number of variables.
 * @param fvec An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param tol A solver_control object defining the solver control
 *  parameters.
 * @param lsearch A pointer to a line_search_control object defining
 *  the line search control parameters.  If no line search is desired,
 *  simply pass NULL.
 * @param ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param err A pointer to the C error handler object.  If no error
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
void solve_newton(vecfcn fcn, jacobianfcn jac, int n, double *x, double *fvec,
                  const solver_control *tol, const line_search_control *lsearch,
                  iteration_behavior *ib, errorhandler *err);

/** @brief Applies the Levenberg-Marquardt method to solve the nonlinear
 * least-squares problem.
 *
 * @param fcn A pointer to the routine containing the system of
 *  equations to solve.
 * @param jac A pointer to a routine used to compute the Jacobian of
 *  the system of equations.  To let the program compute the Jacobian
 *  numerically, simply pass NULL.
 * @param neqn The number of equations.
 * @param nvar The number of unknowns.  This must be less than or equal
 *  to @p neqn.
 * @param x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 *  N is the number of variables.
 * @param fvec An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param tol A solver_control object defining the solver control
 *  parameters.
 * @param ib On output, an iteration_behavior object containing the
 *  iteration performance statistics.
 * @param err A pointer to the C error handler object.  If no error
 *  handling is desired, simply pass NULL, and errors will be dealt with
 *  by the default internal error handler.  Possible errors that may be
 *  encountered are as follows.
 *  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 *  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
 *      than the number of variables.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      correctly.
 *  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
 *      the allowed number of iterations.
 *  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 *  - NL_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
 *      to small to be practical for the problem at hand.
 */
void solve_nl_least_squares(vecfcn fcn, jacobianfcn jac, int neqn, int nvar,
                            double *x, double *fvec, const solver_control *tol,
                            iteration_behavior *ib, errorhandler *err);

/** @brief Sets defaults for the solver_control type.
 *
 * @param tol The solver_control object.
 */
void set_nonlin_defaults(solver_control *tol);

/** @brief Sets defaults for the line_search_control type.
 *
 * @param ls The line_search_control object.
 */
void set_nonlin_ls_defaults(line_search_control *ls);

/** @brief Utilizes the Nelder-Mead simplex method for finding a minimum
!! value of the specified function.
!!
!! @param[in] fcn A pointer to the routine containing the function on which
!!  to operate.
!! @param[in] nvar The dimension of the problem (number of variables).
!! @param[in,out] x On input, the initial guess at the optimal point.
!!  On output, the updated optimal point estimate.
!! @param[out] f An optional output, that if provided, returns the
!!  value of the function at @p x.
!! @param[in] smplx An optional NVAR-by-(NVAR + 1) matrix, that if supplied
!!  provides an initial simplex geometry (each column is a vertex location).
!!  If not provided (NULL), the solver generates its own estimate of a
!!  starting simplex geometry.
!! @param[in] tol A solver_control object defining the solver control
!!  parameters.
!! @param[out] ib On output, an iteration_behavior object containing the
!!  iteration performance statistics.
!! @param[in] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
!!  - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
!!      the problem as defined in @p fcn.
!!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
!!      available.
!!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
!!      the allowed number of iterations.
 */
void nelder_mead(fcnnvar fcn, int nvar, double *x, double *f, 
                 const double *smplx, const solver_control *tol, 
                 iteration_behavior *ib, errorhandler *err);

/** @brief Utilizes the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
!! for finding a minimum value of the specified function.
!!
!! @param[in] fcn A pointer to the routine containing the function on which
!!  to operate.
!! @param[in] grad An optional pointer to a routine capable of computing
!!  the gradient of the function contained within @p fcn.  If no routine
!!  is supplied (NULL), the solver will numerically estimate the gradient.
!! @param[in] nvar The dimension of the problem (number of variables).
!! @param[in,out] x On input, the initial guess at the optimal point.
!!  On output, the updated optimal point estimate.
!! @param[out] f An optional output, that if provided, returns the
!!  value of the function at @p x.
!! @param[in] tol A solver_control object defining the solver control
!!  parameters.
!! @param[out] ib On output, an iteration_behavior object containing the
!!  iteration performance statistics.
!! @param[in] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
!!  - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
!!      the problem as defined in @p fcn.
!!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
!!      available.
!!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
!!      the allowed number of iterations.
 */
void bfgs(fcnnvar fcn, gradientfcn grad, int nvar, double *x, double *f,
          const solver_control *tol, const line_search_control *lsearch,
          iteration_behavior *ib, errorhandler *err);


/** @brief Initializes a new polynomial object.
 *
 * @param poly The polynomial object to initialize.
 * @param order The order of the polynomial.  This value must be > 0.
 */
void alloc_polynomial(polynomial *poly, int order);

/** @brief Frees resources held by a polynomial object.
 *
 * @param obj The polynomial object.
 */
void free_polynomial(polynomial *poly);

/** @brief Gets the order of the polynomial.
 *
 * @param poly A pointer to the polynomial object.
 * @return The order of the polynomial object.
 */
int get_polynomial_order(const polynomial *poly);

/** @brief Fits a polynomial of the specified order to a data set.
 *
 * @param poly The polynomial object to initialize.
 * @param n The size of the arrays.
 * @param x An N-element array containing the independent variable data
 *  points.  Notice, must be N > @p order.
 * @param y On input, an N-element array containing the dependent
 *  variable data points.  On output, the contents are overwritten.
 * @param order The order of the polynomial (must be >= 1).
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
 *      was specified, or if order is too large for the data set.
 *  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
 */
void fit_polynomial(polynomial *poly, int n, const double *x, double *y,
                    int order, errorhandler *err);

/** @brief Fits a polynomial of the specified order that passes through zero
 * to a data set.
 *
 * @param poly The c_polynomial object to initialize.
 * @param n The size of the arrays.
 * @param x An N-element array containing the independent variable data
 *  points.  Notice, must be N > @p order.
 * @param y On input, an N-element array containing the dependent
 *  variable data points.  On output, the contents are overwritten.
 * @param order The order of the polynomial (must be >= 1).
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
 *      was specified, or if order is too large for the data set.
 *  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
 */
void fit_polynomial_thru_zero(polynomial *poly, int n, const double *x,
                              double *y, int order, errorhandler *err);

/** @brief Evaluates a polynomial at the specified points.
 *
 * @param poly A pointer to the polynomial object.
 * @param n The number of points to evaluate.
 * @param x An N-element array containing the points at which to
 *  evaluate the polynomial.
 * @param y An N-element array where the resulting polynomial outputs
 *  will be written.
 */
void evaluate_polynomial(const polynomial *poly, int n, const double *x,
                         double *y);

/** @brief Evaluates a polynomial at the specified points.
 *
 * @param poly A pointer to the polynomial object.
 * @param n The number of points to evaluate.
 * @param x An N-element array containing the points at which to
 *  evaluate the polynomial.
 * @param y An N-element array where the resulting polynomial outputs
 *  will be written.
 */
void evaluate_polynomial_cmplx(const polynomial *poly, int n,
                               const double complex *x, double complex *y);

/** @brief Computes all the roots of a polynomial by computing the
 * eigenvalues of the polynomial companion matrix.
 *
 * @param poly A pointer to the polynomial object.
 * @param n The size of @p rts.  This value should be the same as the
 *  order of the polynomial.
 * @param rts An N-element array where the roots of the polynomial
 *  will be written.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
 *      there is insufficient memory available.
 *  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
 */
void polynomial_roots(const polynomial *poly, int n, double complex *rts,
                      errorhandler *err);

/** @brief Gets the requested polynomial coefficient by index.  The
 * coefficient index is established as follows: c(1) + c(2) * x +
 * c(3) * x**2 + ... c(n) * x**n-1.
 *
 * @param poly A pointer to the polynomial object.
 * @param ind The polynomial coefficient index (0 < ind <= order + 1).
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 * - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
 *      equal to zero, or if the requested index exceeds the number of
 *      polynomial coefficients.
 */
double get_polynomial_coefficient(const polynomial *poly, int ind,
                                  errorhandler *err);

/** @brief Sets the requested polynomial coefficient by index.  The
 * coefficient index is established as follows: c(1) + c(2) * x +
 * c(3) * x**2 + ... c(n) * x**n-1.
 *
 * @param poly A pointer to the polynomial object.
 * @param ind The polynomial coefficient index (0 < ind <= order + 1).
 * @param x The polynomial coefficient.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 * - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
 *      equal to zero, or if the requested index exceeds the number of
 *      polynomial coefficients.
 */
void set_polynomial_coefficient(polynomial *poly, int ind, double x,
                                errorhandler *err);

/** @brief Adds two polynomials.
 *
 * @param p1 The left-hand-side argument.
 * @param p2 The right-hand-side argument.
 * @param rst The resulting polynomial.
 */
void polynomial_add(const polynomial *p1, const polynomial *p2, 
                    polynomial *rst);

/** @brief Subtracts two polynomials.
 *
 * @param p1 The left-hand-side argument.
 * @param p2 The right-hand-side argument.
 * @param rst The resulting polynomial.
 */
void polynomial_subtract(const polynomial *p1, const polynomial *p2, 
                         polynomial *rst);

/** @brief Multiplies two polynomials.
 *
 * @param p1 The left-hand-side argument.
 * @param p2 The right-hand-side argument.
 * @param rst The resulting polynomial.
 */
void polynomial_multiply(const polynomial *p1, const polynomial *p2,
                         polynomial *rst);

/** @brief Copies the contents of one polynomial object to another.
 *
 * @param src The source polynomial object.
 * @param dst The destination polynomial.
 */
void polynomial_copy(const polynomial *src, polynomial *dst);

#ifdef __cplusplus
}
#endif
#endif // NONLIN_H_DEFINED
