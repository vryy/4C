/*!----------------------------------------------------------------------
\file wall1_input.cpp
\brief

<pre>
Maintainer: Markus Gitterle 
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
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

  <pre>                                                         mgit 03/07
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}
#include "wall1.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              mwgit 03/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::Wall1::ReadElement()
{
  // read element's nodes
  int ierr=0;
  int nnode=0;
  int nodes[9]; 
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    nnode = 4;
    frint_n("QUAD4",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    nnode = 8;
    frint_n("QUAD8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    nnode = 9;
    frint_n("QUAD9",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    nnode = 3;
    frint_n("TRI3",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI6",&ierr);
  if (ierr==1)
  {
    nnode = 6;
    frint_n("TRI6",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;
  
  SetNodeIds(nnode,nodes);
  
  // read number of material model
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");
  
  // read wall thickness
  thickness_ = 1.0;
  frdouble("THICK",&thickness_,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");

  // read gaussian points
  frint_n("GP",ngp_,2,&ierr);
  if (ierr!=1) dserror("Reading of WALL1 element failed");

//  // read gaussian points for triangle element
//  frint("GP_TRI",&ngptri_,&ierr);
//  if (ierr!=1) dserror("Reading of WALL1 element failed");

  //read 2D problem type
  frchk("PLANE_STRESS",&ierr);
  if (ierr==1) wtype_ = plane_stress; 
  frchk("PLANE_STRAIN",&ierr);
  if (ierr==1) wtype_ = plane_strain;
  
  //read model (EAS or not)
  frchk("EAS_Model",&ierr);
  if (ierr==1) iseas_=true;

  //read lokal or global stresses
  char buffer [50];
  frchar("STRESSES",buffer,&ierr);
  if (ierr)
  {
     if      (strncmp(buffer,"XY",2)==0)       stresstype_ = w1_xy;
     else if (strncmp(buffer,"RS",2)==0)       stresstype_ = w1_rs;
     else dserror ("Reading of WALL1 element failed");
  }

  //read TSI
#if 0
  tsi_couptyp = tsi_coup_none;  /* default */
  char buffer[50];
  frchar("TSI_COUPTYP",buffer,&ierr);
  if (ierr)
  {
    if (strncmp(buffer,"None",4)==0) tsi_couptyp_ = tsi_coup_none;
    if (strncmp(buffer,"Thermconf",9)==0) tsi_couptyp_ = tsi_coup_thermconf;
    if (strncmp(buffer,"Thermcreate",11)==0) tsi_couptyp_ = tsi_coup_thermcreate;
  }
#endif
        
  return true;
} // Wall1::ReadElement()




#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
