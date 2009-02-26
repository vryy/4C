/*!----------------------------------------------------------------------
\file so_hex27_service.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef D_SOLID3
#include "so_hex27.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  return Center Coords in Reference System                   maf 11/07|
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_hex27::soh27_ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOH27,NUMDIM_SOH27> xrefe;  // material coord. of element
  for (int i=0; i<NUMNOD_SOH27; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOH27,1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  LINALG::Matrix<1,NUMDIM_SOH27> midpoint;
  //midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  vector<double> centercoords(3);
  centercoords[0] = midpoint(0,0);
  centercoords[1] = midpoint(0,1);
  centercoords[2] = midpoint(0,2);
  return centercoords;
}

#endif
#endif
