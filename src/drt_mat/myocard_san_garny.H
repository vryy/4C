/*----------------------------------------------------------------------*/
/*! \file
\brief san garny myocard material model

\level 3

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef Myocard_SAN_Garny_H
#define Myocard_SAN_Garny_H

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 08/13 |
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
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac
/// muscle tissue
///
/// <h3>References</h3>
/// <ul>
/// <li> [1] A Garny et. al., "One-dimensional Rabbit Sinoatrial Node Models: Benefits and
/// Limitations", The Journal of Cardiovascular Electrophysiology 14 (2003) 121-132
/// </ul>
///
/// \author ljag


/// \date 09/13

class Myocard_SAN_Garny : public Myocard_General

{
 public:
  /// construct empty material object
  Myocard_SAN_Garny();

  /// construct empty material object
  explicit Myocard_SAN_Garny(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~Myocard_SAN_Garny() { ; }

  /// compute reaction coefficient
  double ReaCoeff(const double phi, const double dt);

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val);

  ///  return number of ionic currents
  int GetNumberOfIonicCurrents() const;

  ///  return ionic currents
  double GetIonicCurrents(const int k) const;

  /// time update for this material
  void Update(const double phi, const double dt);

 private:
  Myocard_Tools tools_;

  /// perturbation for numerical approximation of the derivative
  double eps_deriv_;

  /// gating variables SAN_Garny
  std::vector<double> s0_;
  std::vector<double> s_;
  std::vector<double> r_;
  std::vector<double> a_;
  std::vector<double> c_;

};  // Myocard_SAN_Garny


#endif
