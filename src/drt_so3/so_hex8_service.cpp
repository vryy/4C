/*!----------------------------------------------------------------------
\file so_hex8_service.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_utils.H"
using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  return Center Coords in Reference System                   maf 11/07|
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_hex8::soh8_ElementCenterRefeCoords()
{
  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  Epetra_SerialDenseVector funct(NUMNOD_SOH8);
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  Epetra_SerialDenseMatrix midpoint(1,NUMDIM_SOH8);
  midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  vector<double> centercoords(3);
  centercoords[0] = midpoint(0,0);
  centercoords[1] = midpoint(0,1);
  centercoords[2] = midpoint(0,2);
  return centercoords;
}

#endif
#endif
