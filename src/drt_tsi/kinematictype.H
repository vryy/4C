/*----------------------------------------------------------------------*/
/*!
\file kinematictype.H

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 11/12 |
 *----------------------------------------------------------------------*/
#ifndef KINEMATIC_H
#define KINEMATIC_H


//! kinematic types of the underlying elements
//! default: total Lagrange, i.e. nonlinear
enum GenKinematicType
{
  geo_linear,
  geo_nonlinear,
  geo_none
};


/*----------------------------------------------------------------------*/
#endif  // KINEMATICTYPE_H

