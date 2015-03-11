/*----------------------------------------------------------------------*/
/*!
\file acou_utils.cpp

\brief utility functions for acoustic problems

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/

#include <cmath>

#include "acou_utils.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/


void ACOU::FillDIRKValues(INPAR::ACOU::DynamicType scheme,
                          double (&a)[6][6],
                          double (&b)[6],
                          double (&c)[6],
                          int &q)
{
  // TODO: no fill, but object with values
  switch(scheme)
  {
  case INPAR::ACOU::acou_dirk23:
  {
    a[0][0] = 0.5 + 0.5 / sqrt(3.0);  a[0][1] = 0.0;
    a[1][0] = -1.0 / sqrt(3.0);       a[1][1] = 0.5 + 0.5 / sqrt(3.0);
    b[0] = 0.5;
    b[1] = 0.5;
    c[0] = 0.5 + 0.5 / sqrt(3.0);
    c[1] = 0.5 - 0.5 / sqrt(3.0);
    q = 2;
    break;
  }
  case INPAR::ACOU::acou_dirk33:
  {
    double alpha = 0.435866521508459;
    a[0][0] = alpha;                                 a[0][1] = 0.0;                                   a[0][2] = 0.0;
    a[1][0] = 0.5 - alpha / 2.0;                     a[1][1] = alpha;                                 a[1][2] = 0.0;
    a[2][0] = -(6.0*alpha*alpha-16.0*alpha+1.0)/4.0; a[2][1] = (6.0*alpha*alpha-20.0*alpha+5.0)/4.0;  a[2][2] = alpha;
    b[0] = -(6.0*alpha*alpha-16.0*alpha+1.0)/4.0;
    b[1] = (6.0*alpha*alpha-20.0*alpha+5.0)/4.0;
    b[2] = alpha;
    c[0] = alpha;
    c[1] = (1.0 + alpha) / 2.0;
    c[2] = 1.0;
    q = 3;
    break;
  }
  case INPAR::ACOU::acou_dirk34:
  {
    double alpha = 1.1371580426;
    a[0][0] = (1.0 + alpha) / 2.0; a[0][1] = 0.0;                 a[0][2] = 0.0;
    a[1][0] = -alpha / 2.0;        a[1][1] = (1.0 + alpha) / 2.0; a[1][2] = 0.0;
    a[2][0] = 1.0 + alpha;         a[2][1] = -1.0 - 2.0 * alpha;  a[2][2] = (1.0 + alpha) / 2.0;
    b[0] = 1.0 / 6.0 / alpha / alpha;
    b[1] = 1.0 - 1.0 / 3.0 / alpha / alpha;
    b[2] = 1.0 / 6.0 / alpha / alpha;
    c[0] = (1.0 + alpha) / 2.0;
    c[1] = 0.5;
    c[2] = (1.0 - alpha) / 2.0;
    q = 3;
    break;
  }
  case INPAR::ACOU::acou_dirk54: // from "Diagonally Implicit Runge-Kutta Formulae with Error Estimates", J.R. Cash, J. Inst. MAths Applics, 1979, 24, pp. 293-301, Equation (2.3)
  {
    a[0][0] = 0.4358665215;   a[0][1] = 0.0;             a[0][2] = 0.0;              a[0][3] = 0.0;             a[0][4] = 0.0;
    a[1][0] = -1.13586652150; a[1][1] = 0.4358665215;    a[1][2] = 0.0;              a[1][3] = 0.0;             a[1][4] = 0.0;
    a[2][0] = 1.08543330679;  a[2][1] = -0.721299828287; a[2][2] = 0.4358665215;     a[2][3] = 0.0;             a[2][4] = 0.0;
    a[3][0] = 0.416349501547; a[3][1] = 0.190984004184;  a[3][2] = -0.118643265417;  a[3][3] = 0.4358665215;    a[3][4] = 0.0;
    a[4][0] = 0.896869652944; a[4][1] = 0.0182725272734; a[4][2] = -0.0845900310706; a[4][3] = -0.266418670647; a[4][4] = 0.4358665215;
    b[0]    = 0.896869652944; b[1]    = 0.0182725272734; b[2]    = -0.0845900310706; b[3]    = -0.266418670647; b[4]    = 0.4358665215;
    c[0]    = 0.4358665215;   c[1]    = -0.7;            c[2]    = 0.8;              c[3]    = 0.924556761814;  c[4]    = 1.0;
    q = 5;
    break;
  }
  default:
    dserror("unknown DIRK scheme");
    break;
  }

  return;
};

std::string ACOU::DIRKTypeToString(INPAR::ACOU::DynamicType scheme)
{
  std::string s = "";
  switch(scheme)
  {
  case INPAR::ACOU::acou_dirk23:      s = "DIRK23";      break;
  case INPAR::ACOU::acou_dirk33:      s = "DIRK33";      break;
  case INPAR::ACOU::acou_dirk34:      s = "DIRK34";      break;
  case INPAR::ACOU::acou_dirk54:      s = "DIRK54";      break;
  default: dserror("no string for DirkType defined"); break;
  }
  return s;
}





