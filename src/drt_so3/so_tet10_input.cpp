/*!----------------------------------------------------------------------**###
\file so_tet10_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet10.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if      (buffer=="Geolin")    kintype_ = so_tet10_geolin;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="Totlag")    kintype_ = so_tet10_totlag;
  // geometrically non-linear with Updated Lagrangean approach
  else if (buffer=="Updlag")
  {
    kintype_ = so_tet10_updlag;
    dserror("Updated Lagrange for SO_TET10 is not implemented!");
  }
  else dserror("Reading of SO_TET10 element failed");

  return true;
}


#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::ReadElement()
{
  int ierr=0;
  const int nnode=10;
  int nodes[10];
  frchk("SOLIDT10",&ierr);
  if (ierr==1)
  {
    frint_n("TET10",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDTET10 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_TET10 element material failed");
  SetMaterial(material);

  // we expect kintype to be total lagrangian
  kintype_ = so_tet10_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = so_tet10_geolin;
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = so_tet10_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = so_tet10_updlag;
       dserror("Updated Lagrange for SO_TET10 is not implemented!");
   }
   else dserror("Reading of SO_TET10 element failed");
  }

  return true;
} // So_tet10::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
