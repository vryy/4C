/*!----------------------------------------------------------------------
\file so_hex8_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

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
#include "so_hex8.H"
#include "../drt_lib/dstrc.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::So_hex8::ReadElement()
{
  DSTraceHelper dst("So_hex8::ReadElement");

  // read element's nodes
  int ierr=0;
  const int nnode=8;
  int nodes[8];
  frchk("SOLIDH8",&ierr);
  if (ierr==1)
  {
    frint_n("HEX8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDH8 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of SO_HEX8 element material failed");

  // read gaussian points
  frint_n("GP",ngp_,3,&ierr);
  if (ierr!=1) dserror("Reading of SOLID3 element gp failed");
  for (int i=0; i<3; ++i) if (ngp_[i]!=2) dserror("Only 2 GP for HEX8 yet");

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = soh8_geolin;
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = soh8_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = soh8_updlag;
       dserror("Updated Lagrange for SO_HEX8 is not implemented!");
   }
   else dserror("Reading of SO_HEX8 element failed");
  }
  
//  // read EAS technology flag
//  eastype_ = soh8_easnone;     // default: no EAS
//  frchar("EAS",buffer,&ierr);
//  if (ierr){
//    // full EAS technology
//    if      (strncmp(buffer,"full",4)==0){
//      eastype_ = soh8_easfull;
//      neas_ = 9;               // number of eas parameters for full EAS
//      soh8_easinit();
//    }
//    // mild EAS technology
//    else if (strncmp(buffer,"mild",4)==0){
//      eastype_ = soh8_easmild;
//      neas_ = 9;               // number of eas parameters for mild EAS
//      soh8_easinit();
//    }
//    // no EAS technology
//    else if (strncmp(buffer,"none",4)==0) eastype_ = soh8_easnone;
//    else dserror("Reading of SO_HEX8 EAS technology failed");
//  } 

  // read stress evaluation/output type
  frchar("STRESS",buffer,&ierr);
  if (ierr!=1) dserror("reading of SO_HEX8 stress failed");
  if (strncmp(buffer,"none",4)==0)  stresstype_= soh8_stress_none;
  if (strncmp(buffer,"Gpxyz",5)==0) stresstype_= soh8_stress_gpxyz;
  if (strncmp(buffer,"Gprst",5)==0) stresstype_= soh8_stress_gprst;
  if (strncmp(buffer,"Gp123",5)==0) stresstype_= soh8_stress_gp123;
  if (strncmp(buffer,"Ndxyz",5)==0) stresstype_= soh8_stress_ndxyz;
  if (strncmp(buffer,"Ndrst",5)==0) stresstype_= soh8_stress_ndrst;
  if (strncmp(buffer,"Nd123",5)==0) stresstype_= soh8_stress_nd123;
  // set default: no stresses
  else stresstype_= soh8_stress_none;

  return true;
} // So_hex8::ReadElement()


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
