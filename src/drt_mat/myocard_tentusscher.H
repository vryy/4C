/*----------------------------------------------------------------------*/
/*!
\file myocard_tentusscher.H
\brief ten tusscher myocard material model

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef Myocard_TenTusscher_H
#define Myocard_TenTusscher_H

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
/// <li> [1] KH ten Tusscher et. al., "Alternans and spiral breakup in a human ventricular tissue
/// model", Am J Physiol Heart Circ Physiol 291 (2006) 1088-1100
/// </ul>
///
/// \author ljag


/// \date 09/13

class Myocard_TenTusscher : public Myocard_General

{
 public:
  /// construct empty material object
  Myocard_TenTusscher();

  /// construct empty material object
  explicit Myocard_TenTusscher(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~Myocard_TenTusscher() { ; }

  /// compute reaction coefficient
  double ReaCoeff(const double phi, const double dt);

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const;

  ///  return current internal state of the material
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

  /// gating variables Inada
  std::vector<double> s0_;
  std::vector<double> s_;
  std::vector<double> r_;
  std::vector<double> a_;
  std::vector<double> c_;

  double VOI_;  // current time (for debugging with CellML)

};  // Myocard_TenTusscher


#endif
