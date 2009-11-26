/*!----------------------------------------------------------------------
\file so_q1p0hex8_input.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "so_hex8p1j1.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_Hex8P1J1::ReadElement(const std::string& eletype,
                                             const std::string& distype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_plneohooke)
  {
    MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(Material().get());
    plastic->Setup(NUMGPT_SOH8);
  }
    
  return true;
}

#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                                 lw 12/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_Hex8P1J1::ReadElement()
{
  // read element's nodes
  int ierr = 0;
  const int nnode = NUMNOD_SOH8;
  int nodes[nnode];
  frchk("SOLIDH8P1J1",&ierr);
  if (ierr==1)
  {
    frint_n("HEX8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of solid Q1P0 hex8 element failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i)
  {
    nodes[i]--;
  }

  SetNodeIds(nnode,nodes);


  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of solid Q1P0 hex8 element material failed");
  SetMaterial(material);

  return true;
} // So_Hex8P1J1::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

