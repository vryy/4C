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

  eps_deriv_ = eps_deriv_myocard;

  // initialization of the material parameters
  a_=0.13;
  b_=0.013;
  c1_=0.26;
  c2_=0.1;
  d_=1.0;

}

double Myocard_Fitzhugh_Nagumo::ComputeReactionCoeff(const double phi, const double dt)
{

     double reacoeff;
     r_ = tools_.GatingVarCalc(dt, r0_, phi/d_, 1.0/(b_*d_));
     reacoeff = c1_*phi*(phi-a_)*(phi-1.0)+c2_*phi*r_;
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
  double val=0.0;
  switch(k){
    case 0: {val=r_; break;}
  }
  return val;
}

/*----------------------------------------------------------------------*
 |  set  internal state of the material                     cbert 08/13 |
 *----------------------------------------------------------------------*/
void Myocard_Fitzhugh_Nagumo::SetInternalState(const int k, const double val)
{
  switch(k){
    case 0: {r0_ = val; r_ = val; break;}
    default: {dserror("There are only 1 internal variables in this material!"); break;}
  }
  return;
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
