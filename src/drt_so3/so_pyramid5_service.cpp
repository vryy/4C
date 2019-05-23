/*!----------------------------------------------------------------------
\brief pyramid shaped solid element
\level 2
\maintainer Christoph Meier

*----------------------------------------------------------------------*/


#include "so_pyramid5.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_node.H"


/*----------------------------------------------------------------------*
 |  return Center Coords in Reference System                            |
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::So_pyramid5::sop5_ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOP5, NUMDIM_SOP5> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOP5; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOP5, 1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.25, distype);
  LINALG::Matrix<1, NUMDIM_SOP5> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}
