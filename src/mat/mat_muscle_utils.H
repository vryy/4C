/*-----------------------------------------------------------*/
/*! \file

\brief Utility methods for skeletal muscle


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef MAT_MUSCLE_UTILS_H
#define MAT_MUSCLE_UTILS_H

#include <vector>
#include "linalg_fixedsizematrix.H"

namespace MAT
{
  namespace UTILS
  {
    namespace MUSCLE
    {
      /*!
       * @brief Structure holding the values of a function f(x) and its derivative dfdx(x) at x
       */
      struct ValuesFunctAndFunctDeriv
      {
        double val_funct;
        double val_deriv_funct;
      };

      /*!
       * @brief Structure holding the values of a function f(x) and its first derivative dfdx(x) and
       * second derivative ddfddx(x) at x
       */
      struct ValuesFunctAndFunctDerivs
      {
        double val_funct;
        double val_deriv_funct;
        double val_deriv_deriv_funct;
      };

      /*!
       * @brief Find a zero of a function using the Newton-Raphson method.
       *
       * Find a zero of the scalar-valued function func given a nearby scalar starting point x_init.
       * The input function needs to be defined as a function that accepts the argument x and that
       * returns a ValuesFunctAndFunctDeriv holding f(x) and dfdx(x)
       *
       * @param[in]     func Scalar function f(x) whose zero should be found
       * @param[in]     x_init Starting point of x
       * @param[in]     tol Tolerance for zero value f(x)
       * @param[in]     maxiter Maximal number of iterations
       * @return        x Point x where f(x) is zero
       */
      double NewtonScalar(const std::function<ValuesFunctAndFunctDeriv(double)> &func,
          const double x_init, const double tol, const int maxiter);

      /*!
       * @brief Find a zero of a continuous function using the bisection method
       *
       * Find a zero of the scalar-valued continuous function func defined on an interval [a, b]
       * given the starting endpoints of the interval a_init and b_init. Note that f(a) and f(b)
       * need to have opposite signs.
       *
       * The input function needs to be defined as a function that accepts the argument x and that
       * returns a ValuesFunctAndFunctDeriv holding f(x) and dfdx(x)
       *
       * @param[in]     func Scalar function f(x) whose zero should be found
       * @param[in]     a_init Starting point of a for the interval [a,b]
       * @param[in]     b_init Starting point of b for the interval [a,b]
       * @param[in]     tol Tolerance for the zero value f(c)
       * @param[in]     maxiter Maximal number of iterations
       * @return        c Midpoint c where f(c) is zero
       */
      double Bisection(const std::function<ValuesFunctAndFunctDeriv(double)> &func,
          const double a_init, const double b_init, const double tol, const int maxiter);

      /*!
       * @brief Evaluate a function f(x) and the first and second derivatives w.r.t. x using central
       * differences
       *
       * @param[in]     func Scalar function f(x) to be evaluated
       * @param[in]     x Point x where f(x) and derivatives should be evaluated
       * @param[in]     delta_x Step size x_{i+1}-x_i
       */
      ValuesFunctAndFunctDerivs EvaluateFunctionAndDerivativesCentralDifferences(
          const std::function<double(double)> &func, const double x, const double delta_x);

      /*!
       * @brief Evaluate the first derivative of f(x) w.r.t. x using central differences
       *
       * @param[in]     f_i_minus_1 Solution f(x_{i-1})
       * @param[in]     f_i_plus_1 Solution f(x_{i+1})
       * @param[in]     delta_x Step size x_{i+1}-x_i
       * @param[out]    dfdx First derivative of f at x_i
       */
      double FirstDerivativeCentralDifferences(
          const double f_i_minus_1, const double f_i_plus_1, const double delta_x);

      /*!
       * @brief Evaluate the second derivative of f(x) w.r.t. x using central differences
       *
       * @param[in]     f_i_minus_1 Solution f(x_{i-1})
       * @param[in]     f_i Solution f(x_{i})
       * @param[in]     f_i_plus_1 Solution f(x_{i+1})
       * @param[in]     delta_x Step size x_{i+1}-x_i
       * @param[out]    ddfddx Second derivative of f at x_i
       */
      double SecondDerivativeCentralDifferences(const double f_i_minus_1, const double f_i,
          const double f_i_plus_1, const double delta_x);

      /*!
       * @brief Evaluate Lambert W function with Halley's method
       *
       * Solution of Lambert W function is functional inverse of xi = W_0*exp(W_0)
       * Computation here restricted to principal branch W_0
       * Use of Halley's method according to:
       * https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-funsolutiction/
       *
       * @param[in]     xi Argument of Lambert W function W(xi)
       * @param[in,out] WO Solution of principal branch of Lambert W function
       * @param[in]     tol Tolerance for Halley's approximation
       * @param[in]     maxiter Maximal number of iterations
       */
      void EvaluateLambert(const double xi, double &W0, const double tol, const int maxiter);

      /*!
       * @brief Evaluate the force-stretch dependency according to Ehret et al.
       *
       * Reference: A. E. Ehret, M. Boel, M. Itskova, 'A continuum constitutive model for the active
       * behaviour of skeletal muscle', Journal of the Mechanics and Physics of Solids, vol. 59, pp.
       * 625-636, 2011, doi: 10.1016/j.jmps.2010.12.008.
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaMin Minimal fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[out]    fxi Force-stretch function result between zero and one
       */
      double EvaluateForceStretchDependencyEhret(
          const double lambdaM, const double lambdaMin, const double lambdaOpt);

      /*!
       * @brief Evaluate the derivative of the force-stretch dependency according to Ehret et al.
       *        w.r.t. the fiber stretch
       *
       * Reference: A. E. Ehret, M. Boel, M. Itskova, 'A continuum constitutive model for the active
       * behaviour of skeletal muscle', Journal of the Mechanics and Physics of Solids, vol. 59, pp.
       * 625-636, 2011, doi: 10.1016/j.jmps.2010.12.008.
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaMin Minimal fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[out]    dFxidLamdaM Derivative of the force-stretch function w.r.t. the fiber
       *                            stretch
       */
      double EvaluateDerivativeForceStretchDependencyEhret(
          const double lambdaM, const double lambdaMin, const double lambdaOpt);

      /*!
       * @brief Evaluate the integral of the force-stretch dependency according to Giantesio et al.
       *        w.r.t. the fiber stretch
       *
       * Reference: G. Giantesio, A. Musesti, 'Strain-dependent internal parameters in hyperelastic
       * biological materials', International Journal of Non-Linear Mechanics, vol. 95, pp. 162-167,
       * 2017, doi:10.1016/j.ijnonlinmec.2017.06.012.
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaMin Minimal fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[out]    intFxi Integral of the force-stretch function w.r.t. the fiber stretch
       */
      double EvaluateIntegralForceStretchDependencyEhret(
          const double lambdaM, const double lambdaMin, const double lambdaOpt);

      /*!
       * @brief Evaluate the force-velocity dependency in the style of Boel et al.
       *
       * Force-velocity dependecy according to Boel et al. with some slight modifications. Boel et.
       * al. scale and shift the eccentric branch of the velocity dependency using the parameter d.
       * To be able to turn the velocity dependency on and off, this concept is adopted for the
       * concentric branch as well. Choosing de=dc=0 leads to a neglection of the velocity
       * dependence. Choosing dc=1 and de=d-1 reproduces the function presented by Boel et. al.
       *
       * Reference: M. Boel and S. Reese, 'Micromechanical modelling of skeletal muscles based on
       * the finite element method', Computer Methods in Biomechanics and Biomedical Engineering,
       * vol. 11, no. 5, pp. 489-504, 2008, doi: 10.1080/10255840701771750.
       *
       * @param[in]     dotLambdaM Contraction velocity / stretch rate
       * @param[in]     dotLambdaMMin Minimal stretch rate
       * @param[in]     de Amplitude of the eccentric velocity dependency function
       * @param[in]     dc Amplitude of the concentric velocity dependency function
       * @param[in]     ke Curvature of the eccentric velocity dependency function
       * @param[in]     kc Curvature of the concentric velocity dependency function
       * @param[out]    fv Force-velocity function result between zero and one
       */
      double EvaluateForceVelocityDependencyBoel(const double dotLambdaM,
          const double dotLambdaMMin, const double de, const double dc, const double ke,
          const double kc);

      /*!
       * @brief Evaluate the derivative of the force-velocity dependency in the style of Boel et al.
       *        w.r.t. the fiber stretch
       *
       * Force-velocity dependecy according to Boel et al. with some slight modifications. Boel et.
       * al. scale and shift the eccentric branch of the velocity dependency using the parameter d.
       * To be able to turn the velocity dependency on and off, this concept is adopted for the
       * concentric branch as well. Choosing de=dc=0 leads to a neglection of the velocity
       * dependence. Choosing dc=1 and de=d-1 reproduces the function presented by Boel et. al.
       *
       * Reference: M. Boel and S. Reese, 'Micromechanical modelling of skeletal muscles based on
       * the finite element method', Computer Methods in Biomechanics and Biomedical Engineering,
       * vol. 11, no. 5, pp. 489-504, 2008, doi: 10.1080/10255840701771750.
       *
       * @param[in]     dotLambdaM Contraction velocity / stretch rate
       * @param[in]     dDotLambdaMdLambdaM Derivative of contraction velocity w.r.t. the fiber
       *                                    stretch
       * @param[in]     dotLambdaMMin Minimal stretch rate
       * @param[in]     de Amplitude of the eccentric velocity dependency function
       * @param[in]     dc Amplitude of the concentric velocity dependency function
       * @param[in]     ke Curvature of the eccentric velocity dependency function
       * @param[in]     kc Curvature of the concentric velocity dependency function
       * @param[out]    dFvdLamdaM Derivative of the force-velocity function w.r.t. the fiber
       *                           stretch
       */
      double EvaluateDerivativeForceVelocityDependencyBoel(const double dotLambdaM,
          const double dDotLambdaMdLambdaM, const double dotLambdaMMin, const double de,
          const double dc, const double ke, const double kc);

      /*!
       * @brief Evaluate the time-dependent optimal active stress by summation single twitches
       * according to Ehret et al.
       *
       * Reference: A. E. Ehret, M. Boel, M. Itskova, 'A continuum constitutive model for the active
       * behaviour of skeletal muscle', Journal of the Mechanics and Physics of Solids, vol. 59, pp.
       * 625-636, 2011, doi: 10.1016/j.jmps.2010.12.008.
       *
       * @param[in]     Na Number of active motor units (MU) per undeformed muscle cross-sectional
       *                   area
       * @param[in]     muTypesNum Number of motor unit (MU) types
       * @param[in]     rho Fraction of MU types
       * @param[in]     I Interstimulus intervals of MU types
       * @param[in]     F Twitch forces of MU types
       * @param[in]     T Twitch contraction times of MU types
       * @param[in]     actIntervalsNum Number of time intervals where activation is prescribed
       * @param[in]     actTimes Time boundaries between intervals
       * @param[in]     actValues Scaling factor in intervals (1=full activation, 0=no activation)
       * @param[in]     currentTime Current time
       * @param[out]    Poptft Time-dependent optimal active stress at currentTime
       */
      double EvaluateTimeDependentActiveStressEhret(const double Na, const int muTypesNum,
          const std::vector<double> &rho, const std::vector<double> &I,
          const std::vector<double> &F, const std::vector<double> &T, const int actIntervalsNum,
          const std::vector<double> &actTimes, const std::vector<double> &actValues,
          const double currentTime);

      /*!
       * @brief Evaluate the active force-stretch dependency according to Blemker et al.
       *
       * Reference: S. S. Blemker, P. M. Pinsky und S. L. Delp, 'A 3D model of muscle reveals the
       * causes of nonuniform strains in the biceps brachii', Journal of biomechanics, vol. 38, no.
       * 4, pp. 657-665, 2005. doi: 10.1016/j.jbiomech.2004.04.009
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[out]    fxi Force-stretch function result between zero and one
       */
      double EvaluateActiveForceStretchDependencyBlemker(
          const double lambdaM, const double lambdaOpt);

      /*!
       * @brief Evaluate the derivative of the active force-stretch dependency according to
       *        Blemker et al. w.r.t. the fiber stretch
       *
       * Reference: S. S. Blemker, P. M. Pinsky und S. L. Delp, 'A 3D model of muscle reveals the
       * causes of nonuniform strains in the biceps brachii', Journal of biomechanics, vol. 38, no.
       * 4, pp. 657-665, 2005. doi: 10.1016/j.jbiomech.2004.04.009
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[out]    dFxidLamdaM Derivative of the force-stretch function w.r.t. the fiber
       *                            stretch
       */
      double EvaluateDerivativeActiveForceStretchDependencyBlemker(
          const double lambdaM, const double lambdaOpt);

      /*!
       * @brief Evaluate the passive force-stretch dependency according to Blemker et al.
       *
       * Reference: S. S. Blemker, P. M. Pinsky und S. L. Delp, 'A 3D model of muscle reveals the
       * causes of nonuniform strains in the biceps brachii', Journal of biomechanics, vol. 38, no.
       * 4, pp. 657-665, 2005. doi: 10.1016/j.jbiomech.2004.04.009
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[in]     lambdaStar Fiber strech where normalized passive fiber force becomes linear
       * @param[in]     P1 Linear material parameter for along-fiber response
       * @param[in]     P2 Exponential material parameter for along-fiber response
       * @param[out]    fxi Force-stretch function result between zero and one
       */
      double EvaluatePassiveForceStretchDependencyBlemker(const double lambdaM,
          const double lambdaOpt, const double lambdaStar, const double P1, const double P2);

      /*!
       * @brief Evaluate the derivative of the passive force-stretch dependency according to
       *        Blemker et al. w.r.t. the fiber stretch
       *
       * Reference: S. S. Blemker, P. M. Pinsky und S. L. Delp, 'A 3D model of muscle reveals the
       * causes of nonuniform strains in the biceps brachii', Journal of biomechanics, vol. 38, no.
       * 4, pp. 657-665, 2005. doi: 10.1016/j.jbiomech.2004.04.009
       *
       * @param[in]     lambdaM Fiber stretch
       * @param[in]     lambdaOpt Optimal fiber stretch
       * @param[in]     lambdaStar Fiber strech where normalized passive fiber force becomes linear
       * @param[in]     P1 Linear material parameter for along-fiber response
       * @param[in]     P2 Exponential material parameter for along-fiber response
       * @param[out]    dFxidLamdaM Derivative of the force-stretch function w.r.t. the fiber
       *                            stretch
       */
      double EvaluateDerivativePassiveForceStretchDependencyBlemker(const double lambdaM,
          const double lambdaOpt, const double lambdaStar, const double P1, const double P2);

      /*!
       * @brief Evaluate the time-dependent optimal (i.e. maximal) active stress using a
       * tanh-function
       *
       * The time-dependent activation function is computed for any t_tot>t_act_start as
       * ft = alpha*tanh(beta*(t_tot-t_act_start)).
       *
       * The time-dependent optimal active stress is obtained by sigma_opt = sigma_max * ft
       *
       * @param[in]     sigma_max optimal (i.e. maximal) active stress
       * @param[in]     alpha tetanised activation level
       * @param[in]     beta scaling factor
       * @param[in]     t_act_start time of start of activation
       * @param[in]     t_current current time
       * @param[out]    sigma_max_ft Time-dependent optimal active stress at t_current
       */
      double EvaluateTimeDependentActiveStressTanh(const double sigma_max, const double alpha,
          const double beta, const double t_act_start, const double t_current);

      /*!
       *  @brief Returns the fiber stretch in the current configuration
       *
       *  @param[in] C Cauchy-Green strain tensor
       *  @param[in] M Structural tensor of fiber directions
       *  @param[out] lambdaM Fiber stretch
       */
      double FiberStretch(const LINALG::Matrix<3, 3> &C, const LINALG::Matrix<3, 3> &M);

      /*!
       *  @brief Returns the derivative of the fiber stretch w.r.t. the Cauchy-Green strain
       *
       *  @param[in] lambdaM Fiber stretch
       *  @param[in] C Cauchy-Green strain tensor
       *  @param[in] M Structural tensor of fiber directions
       *  @param[out] dlambdaMdC Derivative of the fiber stretch w.r.t. the Cauchy-Green strains
       */
      LINALG::Matrix<3, 3> DFiberStretch_DC(
          const double lambdaM, const LINALG::Matrix<3, 3> &C, const LINALG::Matrix<3, 3> &M);

      /*!
       *  @brief Returns the contraction velocity computed by a Backward Euler approximation
       *
       *  @param[in] lambdaM Fiber stretch at time t_n
       *  @param[in] lambdaMOld Fiber stretch at time t_{n-1}
       *  @param[in] timeStepSize Timestep size t_n - t_{n-1}
       *  @return dotLambdaM Contraction velocity, i.e. time derivative of the fiber stretch
       */
      double ContractionVelocityBWEuler(
          const double lambdaM, const double lambdaMOld, const double timeStepSize);
    }  // namespace MUSCLE
  }    // namespace UTILS
}  // namespace MAT
#endif
