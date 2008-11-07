/*!----------------------------------------------------------------------
\file truss3_input.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TRUSS3
#ifdef CCADISCRET

#include "truss3.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
/*----------------------------------------------------------------------*
| vector of material laws defined in global_control.c        cyron 08/08|
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3::ReadElement()
{
  // read element's nodes; in case of a truss element always line2 shape
  int ierr=0;

  //note: BACI intern type is LINE2, but gid input files work with LIN2
  frchk("LIN2",&ierr);
  // two figures have to be read by frint
  const int nnode=2;
  // provide an array of length two in order to store the two node IDs read by frint_n
  int nodes[nnode];
  frint_n("LIN2",nodes,nnode,&ierr);
  
  // if that does not work try LINE2, in case .dat file was created with pre_exodus
  if (ierr != 1)
  {
    ierr=0;
    frchk("LINE2",&ierr);
    frint_n("LINE2",nodes,nnode,&ierr);
  }
  
  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one for meeting BACI intern standard
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of Truss3 element failed");

  // read truss cross section
  crosssec_ = 0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Truss3 element failed");
  
  // we expect kintype to be total lagrangian
  kintype_ = tr3_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically non-linear with Total Lagrangean approach
   if (strncmp(buffer,"totlag",6)==0)    kintype_ = tr3_totlag;
   // geometrically non-linear approach with engineering strains
   else if (strncmp(buffer,"engstr",6)==0)   kintype_ = tr3_engstrain;
   
   else dserror("Reading of Truss3 element failed because of unknown kinematic type!");
  }  
  return true;
} // Truss3::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
