#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include <mpi.h>
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

#include "ale3.H"


bool DRT::Elements::Ale3::ReadElement()
{
  // read element's nodes
  int   ierr = 0;
  int   nnode = 0;
  int   nodes[27];

  frchk("HEX8",&ierr);
  if (ierr==1)
  {
    nnode=8;
    frint_n("HEX8",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("HEX20",&ierr);
  if (ierr==1)
  {
    nnode=20;
    frint_n("HEX20",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("HEX27",&ierr);
  if (ierr==1)
  {
    nnode=27;
    frint_n("HEX27",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("TET4",&ierr);
  if (ierr==1)
  {
    nnode=4;
    frint_n("TET4",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("TET10",&ierr); /* rearrangement??????? */
  if (ierr==1)
  {
    dserror("TET10 element not yet tested!!!\n");
    nnode=10;
    frint_n("TET10",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of ALE3 element failed\n");
  if (material_==0) dserror("No material defined for ALE3 element\n");

  return true;
}


#endif
#endif
#endif
