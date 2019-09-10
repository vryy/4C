/*----------------------------------------------------------------------*/
/*! \file
\brief Fitzhugh Nagumo model for myocard material

\level 2

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*
 |  headers                                                  ljag 07/12 |
 *----------------------------------------------------------------------*/

#include <vector>
#include "myocard_fitzhugh_nagumo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Fitzhugh_Nagumo::Myocard_Fitzhugh_Nagumo() {}


/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Fitzhugh_Nagumo::Myocard_Fitzhugh_Nagumo(
    const double eps_deriv_myocard, const std::string tissue, int num_gp)
    : tools_(), r0_(num_gp), r_(num_gp), J1_(num_gp), J2_(num_gp), mechanical_activation_(num_gp)
{
  // Initial condition
  for (int i = 0; i < num_gp; ++i)
  {
    r0_[i] = 0.0;
    r_[i] = r0_[i];
    mechanical_activation_[i] = 0.0;  // to store the variable for activation
  }


  eps_deriv_ = eps_deriv_myocard;

  // initialization of the material parameters
  a_ = 0.13;
  b_ = 0.013;
  c1_ = 0.26;
  c2_ = 0.1;
  d_ = 1.0;


  // Variables for electromechanical coupling
  act_thres_ = 0.2;  // activation threshold (so that activation = 1.0 if mechanical_activation_ >=
                     // act_thres_)
}

double Myocard_Fitzhugh_Nagumo::ReaCoeff(const double phi, const double dt)
{
  return Myocard_Fitzhugh_Nagumo::ReaCoeff(phi, dt, 0);
}

double Myocard_Fitzhugh_Nagumo::ReaCoeff(const double phi, const double dt, int gp)
{
  double reacoeff;
  r_[gp] = tools_.GatingVarCalc(dt, r0_[gp], phi / d_, 1.0 / (b_ * d_));
  J1_[gp] = c1_ * phi * (phi - a_) * (phi - 1.0);
  J2_[gp] = c2_ * phi * r_[gp];
  reacoeff = J1_[gp] + J2_[gp];

  // For electromechanics
  mechanical_activation_[gp] = phi;

  return reacoeff;
}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Fitzhugh_Nagumo::GetNumberOfInternalStateVariables() const { return 1; }

/*----------------------------------------------------------------------*
 |  returns current internal state of the material          cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetInternalState(const int k) const
{
  return GetInternalState(k, 0);
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material       hoermann 09/19 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetInternalState(const int k, int gp) const
{
  double val = 0.0;
  switch (k)
  {
    case -1:
    {
      val = mechanical_activation_[gp];
      break;
    }
    case 0:
    {
      val = r_[gp];
      break;
    }
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                     cbert 08/13 |
 *----------------------------------------------------------------------*/
void Myocard_Fitzhugh_Nagumo::SetInternalState(const int k, const double val)
{
  SetInternalState(k, val, 0);
  return;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                  hoermann 09/16 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
void Myocard_Fitzhugh_Nagumo::SetInternalState(const int k, const double val, int gp)
{
  switch (k)
  {
    case -1:
    {
      mechanical_activation_[gp] = val;
      break;
    }
    case 0:
    {
      r0_[gp] = val;
      r_[gp] = val;
      break;
    }
    default:
    {
      dserror("There are only 1 internal variables in this material!");
      break;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Fitzhugh_Nagumo::GetNumberOfIonicCurrents() const { return 2; }

/*----------------------------------------------------------------------*
 |  returns current internal currents                       cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetIonicCurrents(const int k) const
{
  return GetIonicCurrents(k, 0);
}

/*----------------------------------------------------------------------*
 |  returns current internal currents                    hoermann 09/16 |
 |  for multiple points per element                                     |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetIonicCurrents(const int k, int gp) const
{
  double val = 0.0;
  switch (k)
  {
    case 0:
    {
      val = J1_[gp];
      break;
    }
    case 1:
    {
      val = J2_[gp];
      break;
    }
  }
  return val;
}


/*----------------------------------------------------------------------*
 |  update of material at the end of a time step             ljag 07/12 |
 *----------------------------------------------------------------------*/
void Myocard_Fitzhugh_Nagumo::Update(const double phi, const double dt)
{
  // update initial values for next time step
  r0_ = r_;

  return;
}
