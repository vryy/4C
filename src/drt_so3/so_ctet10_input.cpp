/*!----------------------------------------------------------------------**###
\file so_ctet10_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

extern "C"
{
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}
#include "so_ctet10.H" //**
#include "../drt_lib/dstrc.H"

/*----------------------------------------------------------------------**########
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_ctet10::ReadElement()
{
  DSTraceHelper dst("So_ctet10::ReadElement");


  int ierr=0;
  const int nnode=10;
  int nodes[10];
  frchk("SOLIDCTET10",&ierr);
  if (ierr==1)
  {
    frint_n("CTET10",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDCTET10 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_CTET10 element material failed");
  SetMaterial(material);

  // read gaussian points
  frint_n("GP",ngp_,3,&ierr);
  /*if (ierr!=1) dserror("Reading of SO_CTET10 element gp failed");
  for (int i=0; i<3; ++i) if (ngp_[i]!=4) dserror("Only 2 GP for TET10");
*/
  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = so_ctet10_geolin;
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = so_ctet10_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = so_ctet10_updlag;
       dserror("Updated Lagrange for SO_CTET10 is not implemented!");
   }
   else dserror("Reading of SO_CTET10 element failed");
  }

  // read stress evaluation/output type
  frchar("STRESS",buffer,&ierr);
  if (ierr!=1) dserror("reading of SO_CTET10 stress failed");
  if (strncmp(buffer,"none",4)==0)  stresstype_= so_ctet10_stress_none;
  if (strncmp(buffer,"Gpxyz",5)==0) stresstype_= so_ctet10_stress_gpxyz;
  if (strncmp(buffer,"Gprst",5)==0) stresstype_= so_ctet10_stress_gprst;
  if (strncmp(buffer,"Gp123",5)==0) stresstype_= so_ctet10_stress_gp123;
  if (strncmp(buffer,"Ndxyz",5)==0) stresstype_= so_ctet10_stress_ndxyz;
  if (strncmp(buffer,"Ndrst",5)==0) stresstype_= so_ctet10_stress_ndrst;
  if (strncmp(buffer,"Nd123",5)==0) stresstype_= so_ctet10_stress_nd123;
  // set default: no stresses
  else stresstype_= so_ctet10_stress_none;

  return true;
} // So_ctet10::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET
