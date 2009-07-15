/*!----------------------------------------------------------------------**###
\file so_ptet_input.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_ptet.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                gee 05/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Ptet::ReadElement()
{
  int ierr=0;
  const int nnode=4;
  int nodes[4];
  frchk("PTET4",&ierr);
  if (ierr==1)
  {
    frint_n(" TET4 ",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  }
  else
  {
    dserror ("Reading of PTET4 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of PTET4 element material failed");
  SetMaterial(material);

  return true;
} // Ptet::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
