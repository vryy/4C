/*----------------------------------------------------------------------*/
/*!
\file myocard.H
\brief

<pre>
Maintainer: Cristobal Bertoglio
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_MINIMAL_H
#define MYOCARD_MINIMAL_H

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
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac muscle tissue
///
/// <h3>References</h3>
/// <ul>
/// <li> [1] A Bueno-Orovio et. al., "Minimal model for human ventricular action potentials in tissue", Journal of Theoretical Biology 253 (2008) 544-560
/// </ul>
///
/// \author cbert


/// \date 08/13

class Myocard_Minimal : public Myocard_General

{

public:

  /// construct empty material object
  Myocard_Minimal();

  /// construct empty material object
  explicit Myocard_Minimal(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~Myocard_Minimal() { ; }

  /// compute reaction coefficient
  double ComputeReactionCoeff(const double phi, const double dt);

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

   /// last gating variables MV
  double v0_; /// fast inward current
  double w0_; /// slow outward current
  double s0_; /// slow inward current

  /// current gating variables MV
  double v_; /// fast inward current
  double w_; /// slow outward current
  double s_; /// slow inward current

  /// ionic currents MV
  double Jfi_; /// fast inward current
  double Jso_; /// slow outward current
  double Jsi_; /// slow inward current

  /// model parameters
  double u_o_;// = 0.0;
  double u_u_;// = 1.55;//1.58;
  double Theta_v_;// = 0.3;
  double Theta_w_;// = 0.13;//0.015;
  double Theta_vm_;// = 0.006;//0.015;
  double Theta_o_;// = 0.006;
  double Tau_v1m_;// = 60.0;
  double Tau_v2m_;// = 1150.0;
  double Tau_vp_;// = 1.4506;
  double Tau_w1m_;// = 60.0;//70.0;
  double Tau_w2m_;// = 15.0;//20.0;
  double k_wm_;// = 65.0;
  double u_wm_;// = 0.03;
  double Tau_wp_;// = 200.0;//280.0;
  double Tau_fi_;// = 0.11;
  double Tau_o1_;// = 400.0;//6.0;
  double Tau_o2_;// = 6.0;
  double Tau_so1_;// = 30.0181;//43.0;
  double Tau_so2_;// = 0.9957;//0.2;
  double k_so_;// = 2.0458;//2.0;
  double u_so_;// = 0.65;
  double Tau_s1_;// = 2.7342;
  double Tau_s2_;// = 16.0;//3.0;
  double k_s_;// = 2.0994;
  double u_s_;// = 0.9087;
  double Tau_si_;// = 1.8875;//2.8723;
  double Tau_winf_;// = 0.07;
  double w_infs_;// = 0.94;

  // Variables for electromechanical coupling
  double mechanical_activation_; // to store the variable for activation (phi in this case=)
  double act_thres_ ; // activation threshold (so that activation = 1.0 if mechanical_activation_ >= act_thres_)

}; // Myocard_Minimal


#endif
