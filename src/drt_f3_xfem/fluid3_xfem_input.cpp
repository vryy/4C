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
  dsassert(ierr!=1, "Reading of material for FLUID3 element failed\n");
  dsassert(material_==0, "No material defined for FLUID3 element\n");

  // read/set gaussian rule
  gaussrule_ = get_optimal_gaussrule(Shape());

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
