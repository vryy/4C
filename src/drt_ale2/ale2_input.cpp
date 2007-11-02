#ifdef D_ALE
#ifdef CCADISCRET

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

#include "ale2.H"


bool DRT::Elements::Ale2::ReadElement()
{
  // read element's nodes
  int   ierr = 0;
  int   nnode = 0;
  int   nodes[9];

  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    nnode=4;
    frint_n("QUAD4",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    nnode=8;
    frint_n("QUAD8",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    nnode=9;
    frint_n("QUAD9",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    nnode=3;
    frint_n("TRI3",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("TRI6",&ierr); /* rearrangement??????? */
  if (ierr==1)
  {
    nnode=6;
    frint_n("TRI6",nodes,nnode,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of ALE2 element failed\n");
  if (material==0) dserror("No material defined for ALE2 element\n");
  SetMaterial(material);

  return true;
}


#endif
#endif
