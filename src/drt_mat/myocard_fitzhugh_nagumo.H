/*----------------------------------------------------------------------*/
/*!
\file myocard_fitzhugh_nagumo.H

\brief Fitzhugh Nagumo model for myocard material

\level 2

<pre>
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
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
  explicit Myocard_Fitzhugh_Nagumo(const double eps_deriv_myocard, const std::string tissue, int num_gp);

  /// destructor
  virtual ~Myocard_Fitzhugh_Nagumo() { ; }

  /// compute reaction coefficient
  double ReaCoeff(const double phi, const double dt);

  /// compute reaction coefficient for multiple points per element
  double ReaCoeff(const double phi, const double dt, int gp);

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const;

  ///  returns current internal state of the material for multiple points per element
  double GetInternalState(const int k, int gp) const;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val);

  ///  set internal state of the material for multiple points per element
  void SetInternalState(const int k, const double val, int gp);

  ///  return number of ionic currents
  int GetNumberOfIonicCurrents() const;

  ///  return ionic currents
  double GetIonicCurrents(const int k) const;

  ///  return ionic currents for multiple points per element
  double GetIonicCurrents(const int k, int gp) const;

  /// time update for this material
  void Update(const double phi, const double dt);

  /// get number of Gauss points
  int GetNumberOfGP() const
  {
    return r0_.size();
  };

  /// resize internal state variables if number of Gauss point changes
  void ResizeInternalStateVariables(int gp)
  {
    r0_.resize(gp);
    r_.resize(gp);
    J1_.resize(gp);
    J2_.resize(gp);
    mechanical_activation_.resize(gp);
  }

private:

  Myocard_Tools tools_;

  /// perturbation for numerical approximation of the derivative
  double eps_deriv_;

   /// last gating variables MV
  std::vector<double> r0_; /// fast inward current

  /// current gating variables MV
  std::vector<double> r_; /// fast inward current

  /// ionic currents
  std::vector<double> J1_;
  std::vector<double> J2_;

  /// model parameters
  double a_;
  double b_;
  double c1_;
  double c2_;
  double d_;

  // Variables for electromechanical coupling
  std::vector<double> mechanical_activation_; // to store the variable for activation (phi in this case=)
  double act_thres_ ; // activation threshold (so that activation = 1.0 if mechanical_activation_ >= act_thres_)


}; // Myocard_Fitzhugh_Nagumo


#endif
