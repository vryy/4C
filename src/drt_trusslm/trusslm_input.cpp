/*!-----------------------------------------------------------------------------------------------------------
 \file trusslm_input.cpp
 \brief three dimensional total Lagrange truss element

<pre>
Maintainer: Kei Mueller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "trusslm.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::TrussLm::ReadElement(const std::string& eletype,
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

  // geometrically non-linear with Total Lagrangean approach
  if (buffer=="totlag")
    kintype_ = trlm_totlag;

  /*/ geometrically non-linear approach with engineering strains
  else if (buffer=="engstr")
    kintype_ = trlm_engstrain;*/

  else
    dserror("Reading of Torsion3 element failed because of unknown kinematic type!");

  return true;
}


#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::TrussLm::ReadElement()
{
  // read element's nodes; in case of a trusslm element line4 shape
  int ierr=0;

  //note: BACI intern type is LINE4, but gid input files work with LIN4
  frchk("LIN4",&ierr); //originally LIN2
  // two figures have to be read by frint
  const int nnode=4; //originally (truss3): const int nnode=2;
  // provide an array of length two in order to store the two node IDs read by frint_n
  int nodes[nnode];
  frint_n("LIN4",nodes,nnode,&ierr);//originally LIN2

  // if that does not work try LINE2, in case .dat file was created with pre_exodus
  if (ierr != 1)
  {
    ierr=0;
    frchk("LINE4",&ierr);
    frint_n("LINE4",nodes,nnode,&ierr);
  }

  if (ierr != 1) dserror("Reading of ELEMENT Topology failed");

  // reduce node numbers by one for meeting BACI intern standard
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of TrussLm element failed");
  SetMaterial(material_);

  // read truss cross section
  crosssec_ = 0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of TrussLm element failed");

  // we expect kintype to be total lagrangian
  kintype_ = trlm_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically non-linear with Total Lagrangean approach
   if (strncmp(buffer,"totlag",6)==0)    kintype_ = trlm_totlag;
   // geometrically non-linear approach with engineering strains
   else if (strncmp(buffer,"engstr",6)==0)   kintype_ = trlm_engstrain;

   else dserror("Reading of TrussLm element failed because of unknown kinematic type!");
  }
  return true;
} // TrussLm::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
