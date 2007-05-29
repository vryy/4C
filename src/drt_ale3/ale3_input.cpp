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
  char  buffer[50];

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

  // read gaussian points

  if (nnode==8 || nnode==20 || nnode==27)
  {
    frint_n("GP",ngp_,3,&ierr);
    if (ierr!=1) dserror("Reading of ALE3 element failed: GP\n");
  }

  // read number of gaussian points for tetrahedral elements */
  if (nnode==4 || nnode==10)
  {
    frint("GP_TET",&ngp_[0],&ierr);
    if (ierr!=1) dserror("Reading of ALE3 element failed: GP_TET\n");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of ALE3 element failed: GP_ALT\n");
    /*
     * integration for TET-elements is distinguished into different cases. This is
     * necessary to get the right integration parameters from ALE_DATA.
     * The flag for the integration case is saved in nGP[1]. For detailed informations
     * see /ale3/f3_intg.c.
     */

    switch (ngp_[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=0;
        else
          dserror("Reading of ALE3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ngp_[1]=2;
        else
          dserror("Reading of ALE3 element failed: GP_ALT\n");
        break;
      case 5:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=3;
        else
          dserror("Reading of ALE3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      default:
        dserror("Reading of ALE3 element failed: integration points\n");
    }
  }

  return true;
}


#endif
#endif
#endif
