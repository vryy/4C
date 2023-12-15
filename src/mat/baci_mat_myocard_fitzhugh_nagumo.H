/*----------------------------------------------------------------------*/
/*! \file
\brief Fitzhugh Nagumo model for myocard material

\level 2

*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 09/13 |
 *----------------------------------------------------------------------*/
#ifndef BACI_MAT_MYOCARD_FITZHUGH_NAGUMO_H
#define BACI_MAT_MYOCARD_FITZHUGH_NAGUMO_H

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 09/13 |
 *----------------------------------------------------------------------*/
#include "baci_comm_parobjectfactory.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_mat_material.H"
#include "baci_mat_myocard_general.H"
#include "baci_mat_myocard_tools.H"
#include "baci_mat_par_parameter.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/// Myocard material according to [1]
///
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac
/// muscle tissue
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
  explicit Myocard_Fitzhugh_Nagumo(
      const double eps_deriv_myocard, const std::string tissue, int num_gp);

  /// compute reaction coefficient
  double ReaCoeff(const double phi, const double dt) override;

  /// compute reaction coefficient for multiple points per element
  double ReaCoeff(const double phi, const double dt, int gp) override;

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const override;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const override;

  ///  returns current internal state of the material for multiple points per element
  double GetInternalState(const int k, int gp) const override;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val) override;

  ///  set internal state of the material for multiple points per element
  void SetInternalState(const int k, const double val, int gp) override;

  ///  return number of ionic currents
  int GetNumberOfIonicCurrents() const override;

  ///  return ionic currents
  double GetIonicCurrents(const int k) const override;

  ///  return ionic currents for multiple points per element
  double GetIonicCurrents(const int k, int gp) const override;

  /// time update for this material
  void Update(const double phi, const double dt) override;

  /// get number of Gauss points
  int GetNumberOfGP() const override { return r0_.size(); };

  /// resize internal state variables if number of Gauss point changes
  void ResizeInternalStateVariables(int gp) override
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
  std::vector<double> r0_;  /// fast inward current

  /// current gating variables MV
  std::vector<double> r_;  /// fast inward current

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
  std::vector<double>
      mechanical_activation_;  // to store the variable for activation (phi in this case=)
  double act_thres_;  // activation threshold (so that activation = 1.0 if mechanical_activation_ >=
                      // act_thres_)


};  // Myocard_Fitzhugh_Nagumo


BACI_NAMESPACE_CLOSE

#endif
