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

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

extern "C"
{
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

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8::ReadElement()
{
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
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_HEX8 element material failed");
  SetMaterial(material);

  // read possible gaussian points, obsolete for computation
  int ngp[3];
  frint_n("GP",ngp,3,&ierr);
  if (ierr==1) for (int i=0; i<3; ++i) if (ngp[i]!=2) dserror("Only 2 GP for So_SH8");

  // we expect kintype to be total lagrangian
  kintype_ = soh8_totlag;
   
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

  // read EAS technology flag
  eastype_ = soh8_easnone;     // default: no EAS
  frchar("EAS",buffer,&ierr);
  if (ierr){
    // full EAS technology
    if      (strncmp(buffer,"full",4)==0){
      eastype_ = soh8_easfull;
      neas_ = 21;               // number of eas parameters for full EAS
      soh8_easinit();
    }
    // mild EAS technology
    else if (strncmp(buffer,"mild",4)==0){
      eastype_ = soh8_easmild;
      neas_ = 9;               // number of eas parameters for mild EAS
      soh8_easinit();
    }
    // no EAS technology
    else if (strncmp(buffer,"none",4)==0) eastype_ = soh8_easnone;
    else dserror("Reading of SO_HEX8 EAS technology failed");
  }

  //Initialize fiber vector
  fiberdirection_.resize(3);

  // Initialize winding flags
  donerewinding_ = false;

  return true;
} // So_hex8::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
