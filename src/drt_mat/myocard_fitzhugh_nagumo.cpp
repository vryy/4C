/*!----------------------------------------------------------------------
\file myocard_fitzhugh_nagumo.cpp

<pre>
Maintainer: Cristobal Bertogloi
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
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
Myocard_Fitzhugh_Nagumo::Myocard_Fitzhugh_Nagumo()
{

}


/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Fitzhugh_Nagumo::Myocard_Fitzhugh_Nagumo(const double eps_deriv_myocard,const std::string tissue)
:
  tools_()
{
  // Initial condition
  r0_ = 0.0; r_ = r0_ ;
  mechanical_activation_ = 0.0;

  eps_deriv_ = eps_deriv_myocard;

  // initialization of the material parameters
  a_=0.13;
  b_=0.013;
  c1_=0.26;
  c2_=0.1;
  d_=1.0;


  // Variables for electromechanical coupling
   mechanical_activation_ = 0.0;; // to store the variable for activation (phi in this case=)
   act_thres_ = 0.2; // activation threshold (so that activation = 1.0 if mechanical_activation_ >= act_thres_)

}

double Myocard_Fitzhugh_Nagumo::ReaCoeff(const double phi, const double dt)
{

     double reacoeff;
     r_ = tools_.GatingVarCalc(dt, r0_, phi/d_, 1.0/(b_*d_));
     J1_ = c1_*phi*(phi-a_)*(phi-1.0);
     J2_ = c2_*phi*r_;
     reacoeff = J1_+J2_;

     // For electromechanics
      mechanical_activation_ = phi;

     return reacoeff;

}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Fitzhugh_Nagumo::GetNumberOfInternalStateVariables() const
{
  return 1;
}

/*----------------------------------------------------------------------*
 |  returns current internal state of the material          cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetInternalState(const int k) const
{
  double val = 0.0;
  switch(k){
    case -1: {val =mechanical_activation_; break;}
    case 0: {val = r_; break;}
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                     cbert 08/13 |
 *----------------------------------------------------------------------*/
void Myocard_Fitzhugh_Nagumo::SetInternalState(const int k, const double val)
{
  switch(k){
    case -1: {mechanical_activation_ = val; break;}
    case 0: {r0_ = val; r_ = val; break;}
    default: {dserror("There are only 1 internal variables in this material!"); break;}
  }
  return;
}

/*----------------------------------------------------------------------*
 |  returns number of internal state variables of the material  cbert 08/13 |
 *----------------------------------------------------------------------*/
int Myocard_Fitzhugh_Nagumo::GetNumberOfIonicCurrents() const
{
  return 2;
}

/*----------------------------------------------------------------------*
 |  returns current internal currents          cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Fitzhugh_Nagumo::GetIonicCurrents(const int k) const
{
  double val=0.0;
  switch(k){
    case 0: {val=J1_; break;}
    case 1: {val=J2_; break;}
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
