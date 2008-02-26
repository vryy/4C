/*!----------------------------------------------------------------------
\file so_sh8_input.cpp
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
#include "so_sh8.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh8::ReadElement()
{
  // read element's nodes
  int ierr=0;
  const int nnode=8;
  int nodes[8];
  frchk("SOLIDSH8",&ierr);
  if (ierr==1)
  {
    frint_n("HEX8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDSH8 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i){
    nodes[i]--;
    inp_nodeIds_[i] = nodes[i];
  }

  SetNodeIds(nnode,nodes);


  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_SH8 element material failed");
  SetMaterial(material);

  // read possible gaussian points, obsolete for computation
  int ngp[3];
  frint_n("GP",ngp,3,&ierr);
  if (ierr==1) for (int i=0; i<3; ++i) if (ngp[i]!=2) dserror("Only 2 GP for So_SH8");

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
       dserror("Updated Lagrange for SO_SH8 is not implemented!");
   }
   else dserror("Reading of SO_SH8 element failed");
  }

  // read EAS technology flag
  eastype_ = soh8_eassosh8;     // default: EAS for Solid-Shell8
  frchar("EAS",buffer,&ierr);
  if (ierr){
    // full EAS technology
    if      (strncmp(buffer,"sosh8",5)==0){
      eastype_ = soh8_eassosh8;
      neas_ = 7;               // number of eas parameters for EAS_SOSH8
      soh8_easinit();
    }
    // no EAS technology
    else if (strncmp(buffer,"none",4)==0){
      cout << "Warning: Solid-Shell8 without EAS" << endl;
      eastype_ = soh8_easnone;
    }
    else dserror("Reading of SO_SH8 EAS technology failed");
  }

  // read global coordinate of shell-thickness direction
  thickdir_ = autoj;           // default: auto by Jacobian
  nodes_rearranged_ = false;
  frchar("THICKDIR",buffer,&ierr);
  if (ierr)
  {
   // global X
   if      (strncmp(buffer,"XDIR",4)==0)    thickdir_ = globx;
   // global Y
   else if (strncmp(buffer,"YDIR",4)==0)    thickdir_ = globy;
   // global Z
   else if (strncmp(buffer,"ZDIR",4)==0)    thickdir_ = globz;
   // find automatically through Jacobian of Xrefe
   else if (strncmp(buffer,"AUTO",4)==0)    thickdir_ = autoj;
   // no noderearrangement
   else if (strncmp(buffer,"NONE",4)==0){
     thickdir_ = none;
     nodes_rearranged_ = true;
   }
   else dserror("Reading of SO_SH8 thickness direction failed");
  }

  return true;
} // So_sh8::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
