/*----------------------------------------------------------------------------*/
/*!

\brief Line search based on polynomials

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_POLYNOMIAL_H
#define LINESEARCH_POLYNOMIAL_H

/*----------------------------------------------------------------------------*/
/* headers */

// baci
#include "linesearch_base.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class LineSearchPolynomial
   *  \brief Line search based on polynomials
   *
   *  Create a series of second order polynomials to model the L2-norm of the
   *  residual along the search direction. Try minimum of this polynomial model
   *  as next line search parameter. Build polynomial models recursively until
   *  the sufficient decrease condition is fulfilled.
   *
   *  <h3> References </h3>
   *  <ul>
   *  <li>[1] Kelley CT: Iterative Methods for Linear and Nonlinear Equations,
   *          SIAM, Philadelphia (1995)
   *  <li>[2] Kelley CT: Iterative Methods for Optimization, SIAM,
   *          Philadelphia (1999)
   *  <li>[3] Kelley CT: Solving Nonlinear Equations with Newton's Method,
   *          SIAM, Philadelphia (2003)
   *  </ul>
   *
   *  \sa NLNSOL::LineSearchBase::IsSufficientDecrease()
   *
   *  \author mayr.mt \date 05/2014
   */
  class LineSearchPolynomial : public NLNSOL::LineSearchBase
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchPolynomial();

    //! Destructor
    virtual ~LineSearchPolynomial(){};

    //@}

    //! @name Setup
    //@{

    //! Setup the line search object
    virtual void Setup();

    //@}

    //! @name Computation
    //@{

    /*! \brief Compute the step length / line search parameter
     *
     *  As suggested in C.T. Kelley's books, a polynomial line search consists
     *  of three steps:
     *
     *  <h3> Outline of the Algorithm </h3>
     *
     *  1. Try a full step: Evaluate residual f(x_i + p) with the previous
     *  solution \f$x_i\f$ and the search direction \f$p\f$. Test for sufficient
     *  decrease (cf. IsSufficientDecrease()).
     *
     *  2. If Step 1 doesn't satisfy sufficient decrease, try a half step, i.e.
     *  \f$f(x_i+0.5p)\f$. Test for sufficient decrease
     *  (cf. IsSufficientDecrease()).
     *
     *  3. If Step 2 doesn't satisfy sufficient decrease, build a second order
     *  polynomial \f$p(\lambda) = a\lambda^2 + b\lambda +c\f$ to model the
     *  L2-norm of the residual along the search direction \f$p\f$. Use the
     *  known pairings \f$(\lambda,p(\lambda))\f$ from the previous iterate
     *  (i.e. \f$(0.0,f(x_i))\f$), from Step 1 (i.e. \f$(1.0,f(x_i+p_i))\f$ and
     *  from Step 2 (i.e. \f$(0.5,f(x_i+0.5p_i))\f$. Then, recursively
     *  further polynomials are built, where we keep the pairing
     *  \f$(0.0,f(x_i))\f$), but always drop the 'oldest' point and add the most
     *  recent one as new data point.
     *
     *  <h3> How to determine the coefficients \f$a,b\f$ and \f$c\f$ of the
     *  quadratic polynomial? </h3>
     *
     *  Using the three data pairings \f$(l_0,y_0),(l_1,y_1),(l_2,y_2)\f$, the
     *  coefficients are computed as
     *  \f[
     *  a = \frac{y_2 - \frac{y_1-y_0}{l_1-l_0}l_2 - y_0 + \frac{y_1-y_0}{l_1-l_0}l_0}
     *      {l_2^2 - \frac{l_1^2-l_0^2}{l_1-l_0}l_2 - l_0^2 + \frac{l_1^2-l_0^2}{l_1-l_0}l_0}
     *  ,\quad
     *  b = \frac{y_1-y_0}{l_1-l_0} - \frac{l_1^2-l_0^2}{l_1-l_0}a;
     *  ,\quad
     *  c = y_0 - l_0^2a - l_0b;
     *  \f]
     *
     *  \note Coefficient \f$c\f$ is not necessary to determine the minimizer
     *  of the polynomial. Therefore, \f$a\f$ and \f$b\f$ are sufficient.
     *
     *  <h3> How to determine the minimizer \f$l_{min}\f$ of the quadratic
     *  polynomial? </h3>
     *
     *  The minimizer is computed as the root of the first derivative, i.e.
     *  \f[
     *  l_{min} = -\frac{b}{2a}
     *  \f]
     *
     *  \note See [Kelley (1995), p. 143] for treatment of special cases.
     *
     *  Line search is finished as soon as the sufficient decrease condition is
     *  satisfied (cf. IsSufficientDecrease()). We stop after #itermax_
     *  evaluations of the polynomial model.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeLSParam(double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr  ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! @name Description of this class
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::LineSearchPolynomial"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "LS_Polynomial"; }

    //@}

   protected:
   private:
    //! maximum number of recursive polynomial evaluations
    int itermax_;

  };  // class LineSearchPolynomial
}  // namespace NLNSOL

#endif /* LINESEARCH_POLYNOMIAL_H */
