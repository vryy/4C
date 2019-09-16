/*----------------------------------------------------------------------*/
/*! \file

\brief Entry routines for FSI problems and some other problem types as well

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_DYN_H
#define FSI_DYN_H

void fluid_ale_drt();
void fluid_xfem_drt();
void fluid_fluid_fsi_drt();
void fluid_freesurf_drt();

/*! \brief Entry routine to all ALE-based FSI algorithms
 *
 *  All ALE-based FSI algorithms rely on a ceratin DOF ordering, namely
 *
 *               structure dof < fluid dof < ale dof
 *
 *  We establish this ordering by calling FillComplete() on the three
 *  discretizations in the order (1) structure (2) fluid (3) ALE.
 */
void fsi_ale_drt();
void xfsi_drt();
void xfpsi_drt();


#endif
