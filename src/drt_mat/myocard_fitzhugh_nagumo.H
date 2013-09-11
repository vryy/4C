/*----------------------------------------------------------------------*/
/*!
\file myocard_fitzhugh_nagumo.H
\brief

<pre>
Maintainer: Cristobal Bertoglio
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 09/13 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_FITZHUGH_NAGUMO_H
#define MYOCARD_FITZHUGH_NAGUMO_H

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 09/13 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "myocard_tools.H"
#include "myocard_general.H"

/*----------------------------------------------------------------------*/
/// Myocard material according to [1]
///
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac muscle tissue
///
/// d(r)/dt =  b*d*(\phi/d-r)
/// d(\phi)/dt + c1*\phi*(\phi-a)*(\phi-1.0) + c2*\phi*r_ = 0
///
/// with \phi the potential variable and r the internal state
/// </ul>
///
/// \author cbert


/// \date 09/13

class Myocard_Fitzhugh_Nagumo : public Myocard_General

{

public:

  /// construct empty material object
  Myocard_Fitzhugh_Nagumo();

  /// construct empty material object
  explicit Myocard_Fitzhugh_Nagumo(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~Myocard_Fitzhugh_Nagumo() { ; }

  /// compute reaction coefficient
  double ComputeReactionCoeff(const double phi, const double dt);

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val);

  /// time update for this material
  void Update(const double phi, const double dt);

private:

  Myocard_Tools tools_;

  /// perturbation for numerical approximation of the derivative
  double eps_deriv_;

   /// last gating variables MV
  double r0_; /// fast inward current

  /// current gating variables MV
  double r_; /// fast inward current

  /// model parameters
  double a_;
  double b_;
  double c1_;
  double c2_;
  double d_;

}; // Myocard_Fitzhugh_Nagumo


#endif
