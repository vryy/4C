/*!----------------------------------------------------------------------
\file brownianrod.H

\brief utility functions for application of finite element method to rod-like continua

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Christian J. Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BROWNIANROD_H_
#define BROWNIANROD_H_

#include "../linalg/linalg_fixedsizematrix.H"


/************************************************************************/
/* Large rotations functions                                            */
/************************************************************************/
/*!
\brief BROWNIANROD: namespace for utility functions which allow for upgrading
any finite rod-element in such a way that in captures Brownian motion according
to the fluctuation dissipation theorem. The implementation was performed with 
respect to
C.J.Cyron, W.A.Wall, 2009, Finite element approach to the Brownian dynamics of polymers, Physical Review E
C.J.Cyron, W.A.Wall, 2010

*/
namespace BROWNIANROD
{
 
  //!computing from quaternion q the Rodrigues parameters omega, Crisfield, Vol. 2, equation (16.79)
  void BrownianRodClassic3D(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& omega);


} //BROWNIANROD

#endif /*BROWNIANROD_H_*/
#endif  // #ifdef CCADISCRET
