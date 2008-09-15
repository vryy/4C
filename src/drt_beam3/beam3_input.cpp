/*!----------------------------------------------------------------------
\file beam3_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#ifdef CCADISCRET

#include "beam3.H"
#include "../drt_lib/standardtypes_cpp.H"
/*----------------------------------------------------------------------*
|                                                        cyron 03/08     |
| vector of material laws                                                |
| defined in global_control.c
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 03/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3::ReadElement()
{
  // read element's nodes; in case of a beam element always line2 shape
  int ierr=0;

  //note: BACI intern type is LINE2, but gid input files work with LIN2
  frchk("LIN2",&ierr);
  // two figures have to be read by frint
  const int nnode=2;
  // provide an array of length two in order to store the two node IDs read by frint_n
  int nodes[nnode];
  frint_n("LIN2",nodes,nnode,&ierr);
  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one for meeting BACI intern standard
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of Beam3 element failed");

  // read beam cross section
  crosssec_ = 0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Beam3 element failed");

  // read beam shear correction and compute corrected cross section
  double shear_correction = 0.0;
  frdouble("SHEARCORR",&shear_correction,&ierr);
  if (ierr!=1) dserror("Reading of Beam3 element failed");
  crosssecshear_ = crosssec_ * shear_correction;

  /*read beam moments of inertia of area; currently the beam3 element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases; the only
   * point which has to be made sure is that the nodal triad T_ is initialized in the registration process 
   * (->beam3.cpp) in such a way that t1 is the unit vector along the beam axis and t2 and t3 are the principal
   * axes with moment of inertia of area Iyy_ and Izz_, respectively; so a modification to more general kinds of
   * cross sections can be done easily by allowing for more complex input right here and by calculating an approxipate
   * initial nodal triad in the frame of the registration; */
  
  frdouble("MOMIN",&Iyy_,&ierr);
  frdouble("MOMIN",&Izz_,&ierr);
  frdouble("MOMINPOL",&Irr_,&ierr);
  if (ierr!=1) dserror("Reading of Beam3 element failed"); 
   
   
  return true;
} // Beam3::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
