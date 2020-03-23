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

/** Defines a pair of values. */
typedef struct {
    /** The first value. */
    double x1;
    /** The second value. */
    double x2;
} value_pair;

/**
 * Defines the signature of a function used for defining an equation of one
 * variable.
 * 
 * @param x The value at which to evaluate the function.
 * @return The value of the function.
 */
typedef double (*fcn1var)(double x);

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
 * Utilizes Newton's method to sovle an equation of one variable.  
 * The derivative calculation will be numerically estimated.
 *
 * @param fcn The function to solve.
 * @param xi On input, the initial guess at the solution.  On
 *  output, the solution.
 * @param fo The value of the function at the solution.
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
int c_solver_newton_1var(fcn1var fcn, double *xi, double *fo, 
    const value_pair *limits, const iteration_controls *cntrls,
    iteration_process *stats);


#ifdef __cplusplus
}
#endif
#endif
