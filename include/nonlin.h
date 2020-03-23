#ifndef NONLIN_H_
#define NONLIN_H_

#include <stdbool.h>

/** A type for providing a set of iteration control parameters. */
typedef struct {
    /** Defines the maximum number of allowable function evaluations. */
    int max_function_evals;
    /** Defines convergence criteria based upon function value. */
    double function_tolerance;
    /** Defines convergence criteria based upon the solution value. */
    double solution_tolerance;
    /**
     * Defines convergence criteria based upon the slope of the gradient 
     * vector. 
     * */
    double gradient_tolerance;
    /**
     * Gets a logical value determining if iteration status should be printed.
     */
    bool print_status;
} iteration_controls;

/** A type providing information on the iteration process. */
typedef struct {
    /** The number of iteration performed. */
    int iteration_count;
    /**
     * The number of function evaluations performed.  This typically
     * does not include derivative evaulations.
     */
    int function_eval_count;
    /** The number of Jacobian evaluations performed. */
    int jacobian_eval_count;
    /** The number of gradient vector evaluations performed. */
    int gradient_eval_count;
    /**
     * True if the solution converged as a result of a zero-valued function;
     * else, false.
     */
    bool converge_on_function;
    /**
     * True if the solution converged as a result of no appreciable
     * change in the solution between iterations; else, false.
     */
    bool converge_on_solution_change;
    /**
     * True if the solution appears to have settled on a stationary
     * point such that the gradient of the function is zero-valued; else,
     * false.
     */
    bool converge_on_gradient;
} iteration_process;

/** A type for providing controls on line search parameters. */
typedef struct {
    /** Defines whether line-searching should be employed. */
    bool enable;
    /**
     * Defines the maximum number of allowable function evaluations
     * during a single line search.
     */
    int max_function_evals;
    /**
     * Defines the scaling of the product of the gradient and 
     * direction vectors such that F(X + LAMBDA * P) <=
     * F(X) + LAMBDA * ALPHA * P**T * G, where P is the search direction
     * vector, G is the gradient vector, and LAMBDA is the scaling factor.
     * The parameter must exist on the set (0, 1).  A value of 1e-4 is
     * typically a good starting point.
     */
    double alpha;
    /**
     * Defines a minimum factor X used to determine a minimum value 
     * LAMBDA such that MIN(LAMBDA) = X * LAMBDA, where LAMBDA defines the 
     * distance along the line search direction assuming a value of one 
     * means the full length of the direction vector is traversed.  As such,
     * the value must exist on the set (0, 1); however, for practical 
     * considerations, the minimum value should be limited to 0.1 such that 
     * the value must exist on the set [0.1, 1).
     */
    double factor;
} line_search_controls;

/** Defines a pair of values. */
typedef struct {
    /** The first value. */
    double x1;
    /** The second value. */
    double x2;
} value_pair;

/**
 * Describes a function of one variable.
 * 
 * @param x The independent variable.
 * @return The value of the function at @p x.
 */
typedef double (*fcn1var)(double x);

/*
 * Describes a function of N variables.
 * 
 * @param n The number of independent variables.
 * @param x An N-element array containing the independent variables.
 * @return The value of the function at @p x.
 */
typedef double (*fcnnvar)(int n, const double *x);

/**
 * Describes a vector-valued function.
 *
 * @param neqn The number of equations.
 * @param nvar The number of independent variables.
 * @param x An NVAR-element array containing the independent
 *  variables.
 * @param f An NEQN-element array containing the values of the
 *  function at @p x.
 */
typedef void (*vecfcn)(int neqn, int nvar, const double *x, double *f);

/*
 * Describes a routin capable of computing the gradient vector
 * of an equation of N variables.
 *
 * @param n The number of independent variables.
 * @param x An N-element array containing the independent variables.
 * @param g An N-element array where the gradient vector will be
 *  written as output.
 */
typedef void (*gradientfcn)(int n, const double *x, double *g);

/*
 * Describes a routine capable of computing the Jacobian matrix
 * of a system of equations.
 *
 * @param neqn The number of equations.
 * @param nvar The number of independent variables.
 * @param x An NVAR-element array containing the independent 
 *  variables.
 * @param jac An NEQN-by-NVAR matrix where the Jacobian will be
 *  written.
 */
typedef void (*jacobianfcn)(int neqn, int nvar, const double *x, double *jac);

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Defines standard solver settings.
 *
 * @param x The iteration_controls item to populate.
 */
void c_set_default_solver_settings(iteration_controls *x);

/**
 * Defines standard line search settings.
 *
 * @param x The line_search_controls item to populate.
 */
void c_set_default_line_search_settings(line_search_controls *x);

/**
 * Utilizes Newton's method to sovle an equation of one variable.  
 * The derivative calculation will be numerically estimated.
 *
 * @param fcn The function to solve.
 * @param limits A set of limits on the solution bounding the solution
 *  range thereby preventing the solver from wandering too far off course.
 * @param x On output, the solution.
 * @param f The value of the function at the solution.
 
 * @param cntrls The iteration controls.
 * @param stats The iteration status.
 *
 * @return An error flag with the following possible values.
 * - NL_NO_ERROR: No error has occurred - successful execution.
 * - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 * - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within 
 *      the allowed number of iterations.
 */
int c_solver_newton_1var(fcn1var fcn, value_pair limits, double *x, 
    double* f, const iteration_controls *cntrls, iteration_process *stats);

/**
 * Utilizes Brent's method to solve an equation of one variable.
 *
 * @param fcn The function to solve.
 * @param x On input, the initial guess.  On output, the solution.
 * @param f The value of the function at the solution.
 * @param limits A set of limits on the solution bounding the solution
 *  range thereby preventing the solver from wandering too far off course.
 * @param cntrls The iteration controls.
 * @param stats The iteration status.
 *
 * @return An error flag with the following possible values.
 * - NL_NO_ERROR: No error has occurred - successful execution.
 * - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 * - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within 
 *      the allowed number of iterations.
 */
int c_solver_brent_1var(fcn1var fcn, double *x, double *f, value_pair limits,
    const iteration_controls *cntrls, iteration_process *stats);

/**
 * Utilizes Broyden's Quasi-Newton method to solve a system of N
 * equations of N unknowns.  A backtracking type line search is also 
 * employed.
 *
 * @param fcn The function to solve.
 * @param jac A function for evaluating the Jacobian.  If null, the
 *  Jacobian is estimated numerically.
 * @param n The number of equations.
 * @param x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 * @param f An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param cntrls The iteration controls.
 * @param ls The line search controls.
 * @param stats The iteration status.
 *
 * @return An error flag with the following possible values.
 * - NL_NO_ERROR: No error has occurred - successful execution.
 * - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      correctly.
 * - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
 *      pointing in an apparent uphill direction.
 * - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
 *      the allowed number of iterations.
 * - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 * - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
 *      gradient vector becomes sufficiently close to zero.
 */
int c_solver_quasi_newton(vecfcn fcn, jacobianfcn jac, int n, double *x,
    double *f, const iteration_controls *cntrls, const line_search_controls *ls,
    iteration_process *stats);

/**
 * Utilizes Newton's method to solve a system of N equations of N 
 * unknowns in conjuction with a backtracking type line search.
 *
 * @param fcn The function to solve.
 * @param jac A function for evaluating the Jacobian.  If null, the
 *  Jacobian is estimated numerically.
 * @param n The number of equations.
 * @param x On input, an N-element array containing an initial
 *  estimate to the solution.  On output, the updated solution estimate.
 * @param f An N-element array that, on output, will contain
 *  the values of each equation as evaluated at the variable values
 *  given in @p x.
 * @param cntrls The iteration controls.
 * @param ls The line search controls.
 * @param stats The iteration status.
 *
 * @return An error flag with the following possible values.
 * - NL_NO_ERROR: No error has occurred - successful execution.
 * - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 *  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
 *      correctly.
 * - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
 *      pointing in an apparent uphill direction.
 * - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
 *      the allowed number of iterations.
 * - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 * - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
 *      gradient vector becomes sufficiently close to zero.
 */
int c_solver_newton(vecfcn fcn, jacobianfcn jac, int n, double *x, double *f, 
    const iteration_controls *cntrls, const line_search_controls *ls,
    iteration_process *stats);

/**
 * Utilizes the Levenberg-Marquardt method to solve a least-squares
 * problem of M equations of N unknowns.  There must be at least as many
 * equations as unknowns for this solver.
 *
 * @param fcn The function to solve.
 * @param jac A function for evaluating the Jacobian.  If null, the
 *  Jacobian is estimated numerically.
 * @param neqn The number of equations.
 * @param nvar The number of variables.
 * @param x On input, an NVAR element array containing the initial 
 *  estimate to the solution.  On output, the solution.
 * @param f An NEQN-element array that, on output, will contain the
 *  values of each equation as evaluated at the output solution given in
 *  @p x.
 * @param cntrls The iteration controls.
 * @param stats The iteration status.
 *
 * @return An error flag with the following possible values.
 * - NL_NO_ERROR: No error has occurred - successful execution.
 * - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
 * - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
 *      than the number of variables.
 * - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
 *      the allowed number of iterations.
 * - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
 *      available.
 * - NL_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
 *      to small to be practical for the problem at hand.
 */
int c_solver_least_squares(vecfcn fcn, jacobianfcn jac, int neqn, int nvar,
    double *x, double *f, const iteration_controls *cntrls, 
    iteration_process *stats);

#ifdef __cplusplus
}
#endif
#endif
