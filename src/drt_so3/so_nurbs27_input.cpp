/*!----------------------------------------------------------------------
\file so_nurbs27_input.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nurbs27.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::NURBS::So_nurbs27::ReadElement()
{
  // read element's nodes
  int ierr=0;
  const int nnode=27;
  int nodes[27];
  frchk("SONURBS27",&ierr);
  if (ierr==1)
  {
    frint_n("NURBS27",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SONURBS27 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_NURBS27 element material failed");
  SetMaterial(material);

  // read possible gaussian points, obsolete for computation
  int ngp[3];
  frint_n("GP",ngp,3,&ierr);
  if (ierr==1) for (int i=0; i<3; ++i) if (ngp[i]!=3) dserror("Only version with 3 GP for So_H27 implemented");

  // we expect kintype to be total lagrangian
  kintype_ = sonurbs27_totlag;

  return true;
} // So_nurbs27::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
