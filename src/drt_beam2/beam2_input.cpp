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

#include "beam2.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*
|                                                        cyron 01/08     |
| vector of material laws                                                |
| defined in global_control.c
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2::ReadElement()
{
  // read element's nodes; in case of a beam element always line2 shape
  int ierr=0;

  //note: BACI intern type is LINE2, but gid input files work with LIN2
  frchk("LIN2",&ierr);
  // two figures have to be read by frint
  int nnode=2;
  // provide an array of length two in order to store the two figures read
  int nodes[2];
  frint_n("LIN2",nodes,nnode,&ierr);
  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");

  // read beam cross section
  crosssec_ = 1.0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");

   // read beam cross section
  double shear_correction = 0.0;
  frdouble("SHEARCORR",&shear_correction,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  crosssecshear_ = crosssec_ * shear_correction;

  // read beam moment of inertia of area
  mominer_ = 1.0;
  frdouble("INERMOM",&mominer_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
  
  // element can use consistent or lumped mass matrix 
  lumpedflag_ = 0;
  frint("LUMPED",&lumpedflag_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2 element failed");
   
  // reading thermal energy kT responsible for statistical thermodynamical forces
    thermalenergy_ = 0;
    frdouble("THERMEN",&thermalenergy_,&ierr);
    if (ierr!=1) dserror("Reading of Beam2 element failed"); 
   
  return true;
} // Beam2::ReadElement()




#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
