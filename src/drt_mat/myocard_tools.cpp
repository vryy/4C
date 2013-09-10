/*!----------------------------------------------------------------------
\file myocard.cpp

<pre>
Maintainer: Cristobal Bertogloi
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 08/13 |
 *----------------------------------------------------------------------*/
#include "myocard_tools.H"
#include <math.h>       /* tanh, log */

/*----------------------------------------------------------------------*
 |  Constructor                                    (public)  cbert 08/13 |
 *----------------------------------------------------------------------*/
Myocard_Tools::Myocard_Tools()
{

}

/*----------------------------------------------------------------------*
 |                                                           cbert 08/13 |
 *----------------------------------------------------------------------*/
double Myocard_Tools::GatingFunction(const double Gate1, const double Gate2, const double p, const double var, const double thresh) const
{
  double Erg = Gate1+(Gate2-Gate1)*(1.0+tanh(p*(var-thresh)))/2;
  return Erg;
}


/*----------------------------------------------------------------------*
 |                                                           cbert 08/13 |
 *----------------------------------------------------------------------*/
double  Myocard_Tools::GatingVarCalc(const double dt, double y_0, const double y_inf, const double y_tau) const
{
  double Erg =  1.0/(1.0/dt + 1.0/y_tau)*(y_0/dt + y_inf/y_tau);
  return Erg;
}
