/*!----------------------------------------------------------------------
\file beam2_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
/*----------------------------------------------------------------------*
|                                                        cyron 01/08     |
| vector of material laws                                                |
| defined in global_control.c
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

extern "C" 
{
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         mgit 03/07
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}


#include "beam2.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2::ReadElement()
{
  // read element's nodes; in case of a beam element always line2 shape
  int ierr=0;
  int nnode=0;
  int nodes[2]; 
  
  frchk("LINE2",&ierr);
  frint_n("LINE2",nodes,nnode,&ierr);
  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;
  
  SetNodeIds(nnode,nodes);
  
  
  // read number of material model and assing young's modulus to ym 
  // note: ym is the only material parameter necessary for Bernoulli beams
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
 

  // read beam cross section
  cross_section_ = 1.0;
  frdouble("CROSS",&cross_section_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  
   // read beam cross section
  double shear_correction = 0.0;
  frdouble("SHEARCORR",&shear_correction,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  cross_section_corr_ = cross_section_ * shear_correction;
  
  // read beam moment of inertia of area
  moment_inertia_ = 1.0;
  frdouble("INERMOM",&moment_inertia_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  

  return true;
} // Beam2::ReadElement()




#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
