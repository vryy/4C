/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.H

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ACOU_H
#define INPAR_ACOU_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  // Input parameters for acoustical simulations
  namespace ACOU
  {
    // Type of time integrator
    enum DynamicType
    {
      acou_impleuler,   // implicit Euler
      acou_trapezoidal, // trapezoidal rule
      acou_dirk23,      // diagonally implicit Runge-Kutta with 2 stages of order 3
      acou_dirk33,      // diagonally implicit Runge-Kutta with 2 stages of order 3
      acou_dirk34,      // diagonally implicit Runge-Kutta with 2 stages of order 3
      acou_dirk54,      // diagonally implicit Runge-Kutta with 2 stages of order 3
      acou_bdf2,        // backward differentiation formula of order 2
      acou_bdf3,        // backward differentiation formula of order 3
      acou_bdf4         // backward differentiation formula of order 4
    };

    // physical type
    enum PhysicalType
    {
      acou_lossless,
      acou_viscous
    };

    /// Calculation of the error
    enum CalcError
    {
      calcerror_no,
      calcerror_1d
    };

    /// Type of Inverse Analysis
    enum InvAnalysisType
    {
      inv_none,          // no inverse analysis
      inv_pat            // photoacoustic inverse analysis
    };

    /// Paramterization of material parameters
    enum InvAnalysisMatParametrization
    {
      inv_mat_none,        // nothing defined
      inv_mat_elementwise, // parameters per element
      inv_mat_uniform      // uniform i.e. spatially constant parameters
    };

  } // namespace ACOU
} // namespace INPAR


#endif // #ifndef INPAR_ACOU_H
