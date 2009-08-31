/*!----------------------------------------------------------------------
\file truss3_input.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TORSION3
#ifdef CCADISCRET

#include "torsion3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Torsion3::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  if (buffer=="totlag")
    kintype_ = tr3_totlag;

  // geometrically non-linear approach with engineering strains
  else if (buffer=="engstr")
    kintype_ = tr3_engstrain;

  else
    dserror("Reading of Torsion2 element failed because of unknown kinematic type!");

  return true;
}


#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Torsion3::ReadElement()
{
  // read element's nodes; in case of a truss element always line2 shape
  int ierr=0;

  //note: BACI intern type is LINE2, but gid input files work with LIN2
  frchk("LIN2",&ierr);
  // two figures have to be read by frint
  const int nnode=2;
  // provide an array of length two in order to store the two node IDs read by frint_n
  int nodes[nnode];
  frint_n("LIN2",nodes,nnode,&ierr);

  // if that does not work try LINE2, in case .dat file was created with pre_exodus
  if (ierr != 1)
  {
    ierr=0;
    frchk("LINE2",&ierr);
    frint_n("LINE2",nodes,nnode,&ierr);
  }

  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one for meeting BACI intern standard
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of Torsion3 element failed");
  SetMaterial(material_);

  // read truss cross section
  crosssec_ = 0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Torsion3 element failed");

  // we expect kintype to be total lagrangian
  kintype_ = tr3_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically non-linear with Total Lagrangean approach
   if (strncmp(buffer,"totlag",6)==0)    kintype_ = tr3_totlag;
   // geometrically non-linear approach with engineering strains
   else if (strncmp(buffer,"engstr",6)==0)   kintype_ = tr3_engstrain;

   else dserror("Reading of Torsion3 element failed because of unknown kinematic type!");
  }
  return true;
} // Torsion3::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION3
