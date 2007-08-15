/*!----------------------------------------------------------------------
\file condif2_input.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
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
#include "condif2.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              gammi 04/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::Condif2::ReadElement()
{
  // read element's nodes
  int   ierr = 0;
  int   nnode = 0;
  int   nodes[27];
  char  buffer[50];
  int   ngp_[2];

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
  if (ierr!=1) dserror("Reading Material for CONDIF2 element failed");
  if (material==0) dserror("No material defined for CONDIF2 element");
  SetMaterial(material);

  // read gaussian points

  if (nnode==4 || nnode==8 || nnode==9)
  {
    frint_n("GP",ngp_,2,&ierr);
    if (ierr!=1) dserror("Reading of CONDIF2 element failed: GP\n");
  }

  // read number of gaussian points for triangle elements */
  if (nnode==3 || nnode==6)
  {
    frint("GP_TRI",&ngp_[0],&ierr);
    if (ierr!=1) dserror("Reading of CONDIF2 element failed: GP_TRI\n");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of CONDIF2 element failed: GP_ALT\n");
    /*
     * integration for TRI-elements is distinguished into different cases.
     * This is necessary to get the right integration parameters from
     * FLUID_DATA.
     * The flag for the integration case is saved in nGP[1].
     * For detailed informations see /fluid3/f3_intg.c.
     */

    switch(ngp_[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=0;
        else
          dserror("Reading of CONDIF2 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      case 3:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ngp_[1]=2;
        else
          dserror("Reading of CONDIF2 element failed: GP_ALT\n");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=3;
        else
          dserror("Reading of CONDIF2 element failed: gauss-radau not possible!\n");
        break;
      case 6:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=4;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ngp_[1]=5;
        else
          dserror("Reading of CONDIF2 element failed: GP_ALT\n");
        break;
      case 7:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=6;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ngp_[1]=7;
        else
        dserror("Reading of CONDIF2 element failed: GP_ALT\n");
      case 9:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=8;
        else
          dserror("Reading of CONDIF2 element failed: gauss-radau not possible!\n");
        break;
      case 12:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=9;
        else
          dserror("Reading of CONDIF2 element failed: gauss-radau not possible!\n");
        break;
      case 13:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=10;
        else
          dserror("Reading of CONDIF2 element failed: gauss-radau not possible!\n");
        break;
        
      default:
        dserror("Reading of CONDIF2 element failed: integration points\n");
    }
  } // end reading gaussian points for tetrahedral elements


  return true;

} // Condif2::ReadElement()


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
