/*----------------------------------------------------------------------*/
/*!
\file kinematictype.H

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 11/12 |
 *----------------------------------------------------------------------*/
#ifndef KINEMATIC_H
#define KINEMATIC_H


//! kinematic types of the underlying elements
//! default: total Lagrange, i.e. nonlinear
// TODO 2012-11-12 typedef enum KinematicType
enum GenKinematicType
{
  geo_linear,
  geo_nonlinear,
  geo_none
};


/*----------------------------------------------------------------------*/
#endif  // KINEMATICTYPE_H
