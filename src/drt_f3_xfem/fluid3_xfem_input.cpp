/*!----------------------------------------------------------------------
\file fluid3_input.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
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
#include "fluid3_xfem.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              g.bau 03/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::XFluid3::ReadElement()
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
  if (ierr!=1) dserror("Reading of FLUID3 element failed\n");
  if (material_==0) dserror("No material defined for FLUID3 element\n");

  // read gaussian points

   if (nnode==8 || nnode==20 || nnode==27)
  {
    frint_n("GP",ngp_,3,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP\n");
  }

  // read number of gaussian points for tetrahedral elements */
  if (nnode==4 || nnode==10)
  {
    frint("GP_TET",&ngp_[0],&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_TET\n");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_ALT\n");
    /*
     * integration for TET-elements is distinguished into different cases. This is
     * necessary to get the right integration parameters from FLUID_DATA.
     * The flag for the integration case is saved in nGP[1]. For detailed informations
     * see /fluid3/f3_intg.c.
     */

    switch(ngp_[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=0;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ngp_[1]=2;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT\n");
        break;
      case 5:
        if (strncmp(buffer,"standard",8)==0)
          ngp_[1]=3;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      default:
        dserror("Reading of FLUID3 element failed: integration points\n");
    }
  } // end reading gaussian points for tetrahedral elements


  // read net algo
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
        strncmp(buffer,"ALE",3)==0 ||
        strncmp(buffer,"Ale",3)==0 )
    {
      is_ale_=true;
    }

    else if (strncmp(buffer,"euler",5)==0 ||
             strncmp(buffer,"EULER",5)==0 ||
             strncmp(buffer,"Euler",5)==0 )
      is_ale_=false;
    else
      dserror("Reading of FLUID3 element failed: Euler/Ale\n");
  }
  else
    dserror("Reading of FLUID3 element net algorithm failed: NA\n");




  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!
  // see src/fluid3/f3_inpele.c for missing details


  return true;

} // Fluid3::ReadElement()


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
