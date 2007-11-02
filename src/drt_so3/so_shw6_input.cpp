/*!----------------------------------------------------------------------
\file so_shw6_input.cpp
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
#include "so_shw6.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::So_shw6::ReadElement()
{
  // read element's nodes
  int ierr=0;
  const int nnode=6;
  int nodes[6];
  frchk("SOLIDSHW6",&ierr);
  if (ierr==1)
  {
    frint_n("WEDGE6",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDSHW6 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i){
    nodes[i]--;
  }

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_WEG6 element material failed");
  SetMaterial(material);

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = sow6_geolin;
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = sow6_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = sow6_updlag;
       dserror("Updated Lagrange for SOLIDSHW6 is not implemented!");
   }
   else dserror("Reading of SOLIDSHW6 element failed");
  }

  
  // read stress evaluation/output type
  frchar("STRESS",buffer,&ierr);
  if (ierr!=1) dserror("reading of SOLIDSHW6 stress failed");
  if (strncmp(buffer,"none",4)==0)  stresstype_= sow6_stress_none;
  if (strncmp(buffer,"Gpxyz",5)==0) stresstype_= sow6_stress_gpxyz;
  if (strncmp(buffer,"Gprst",5)==0) stresstype_= sow6_stress_gprst;
  if (strncmp(buffer,"Gp123",5)==0) stresstype_= sow6_stress_gp123;
  if (strncmp(buffer,"Ndxyz",5)==0) stresstype_= sow6_stress_ndxyz;
  if (strncmp(buffer,"Ndrst",5)==0) stresstype_= sow6_stress_ndrst;
  if (strncmp(buffer,"Nd123",5)==0) stresstype_= sow6_stress_nd123;
  // set default: no stresses
  else stresstype_= sow6_stress_none;

  return true;
} // So_weg6::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
