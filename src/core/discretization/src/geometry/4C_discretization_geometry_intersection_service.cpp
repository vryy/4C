/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations

\level 3

*----------------------------------------------------------------------*/


#include "4C_discretization_geometry_intersection_service.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_geometry_element_coordtrafo.hpp"
#include "4C_discretization_geometry_intersection_service_templates.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ICS:    checks if an element is CARTESIAN, LINEAR and    u.may 07/08|
 |          HIGHERORDER                                                 |
 *----------------------------------------------------------------------*/
void CORE::GEO::checkGeoType(const CORE::Elements::Element* element,
    const CORE::LINALG::SerialDenseMatrix& xyze_element, EleGeoType& eleGeoType)
{
  bool cartesian = true;
  int CartesianCount = 0;
  const int dimCoord = 3;
  const CORE::FE::CellType distype = element->Shape();
  const int eleDim = CORE::FE::getDimension(distype);

  if (CORE::FE::getOrder(distype) == 1)
    eleGeoType = LINEAR;
  else if (CORE::FE::getOrder(distype) == 2)
    eleGeoType = HIGHERORDER;
  else
    FOUR_C_THROW("order of element shapefuntion is not correct");

  // check if cartesian
  if (eleDim == 3)
  {
    const std::vector<std::vector<int>> eleNodeNumbering =
        CORE::FE::getEleNodeNumberingSurfaces(distype);
    std::vector<Teuchos::RCP<CORE::Elements::Element>> surfaces =
        (const_cast<CORE::Elements::Element*>(element))->Surfaces();
    for (int i = 0; i < element->NumSurface(); i++)
    {
      CartesianCount = 0;
      const CORE::Elements::Element* surfaceP = surfaces[i].get();

      for (int k = 0; k < dimCoord; k++)
      {
        int nodeId = eleNodeNumbering[i][0];
        const double nodalcoord = xyze_element(k, nodeId);
        for (int j = 1; j < surfaceP->num_node(); j++)
        {
          nodeId = eleNodeNumbering[i][j];
          if (fabs(nodalcoord - xyze_element(k, nodeId)) > TOL7)
          {
            CartesianCount++;
            break;
          }
        }
      }
      if (CartesianCount > 2)
      {
        cartesian = false;
        break;
      }
    }  // for xfem surfaces
  }    // if eleDim == 3
  else if (eleDim == 2 || eleDim == 1)
  {
    CartesianCount = 0;
    for (int k = 0; k < dimCoord; k++)
    {
      const double nodalcoord = xyze_element(k, 0);
      for (int j = 1; j < element->num_node(); j++)
      {
        if (fabs(nodalcoord - xyze_element(k, j)) > TOL7)
        {
          CartesianCount++;
          break;
        }
      }
    }
    if (CartesianCount > 2) cartesian = false;
  }
  else
    FOUR_C_THROW("dimension of element is not correct");



  if (cartesian) eleGeoType = CARTESIAN;
}


/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          u.may 12/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
std::map<int, CORE::LINALG::Matrix<3, 2>> CORE::GEO::getCurrentXAABBs(
    const DRT::Discretization& dis,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions)
{
  std::map<int, CORE::LINALG::Matrix<3, 2>> currentXAABBs;
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j = 0; j < dis.NumMyColElements(); ++j)
  {
    const CORE::Elements::Element* element = dis.lColElement(j);
    const CORE::LINALG::SerialDenseMatrix xyze_element(
        CORE::GEO::getCurrentNodalPositions(element, currentpositions));
    CORE::GEO::EleGeoType eleGeoType(CORE::GEO::HIGHERORDER);
    CORE::GEO::checkGeoType(element, xyze_element, eleGeoType);
    const CORE::LINALG::Matrix<3, 2> xaabbEle =
        CORE::GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
    currentXAABBs[element->Id()] = xaabbEle;
  }
  return currentXAABBs;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if two 18DOPs intersect                   u.may 12/08| |
 *----------------------------------------------------------------------*/
bool CORE::GEO::intersectionOfKDOPs(
    const CORE::LINALG::Matrix<9, 2>& cutterDOP, const CORE::LINALG::Matrix<9, 2>& xfemDOP)
{
  // check intersection of 18 kdops
  for (int i = 0; i < 9; i++)
    if (!(((cutterDOP(i, 0) > (xfemDOP(i, 0) - CORE::GEO::TOL7)) &&
              (cutterDOP(i, 0) < (xfemDOP(i, 1) + CORE::GEO::TOL7))) ||
            ((cutterDOP(i, 1) > (xfemDOP(i, 0) - CORE::GEO::TOL7)) &&
                (cutterDOP(i, 1) < (xfemDOP(i, 1) + CORE::GEO::TOL7))) ||
            ((xfemDOP(i, 0) > (cutterDOP(i, 0) - CORE::GEO::TOL7)) &&
                (xfemDOP(i, 0) < (cutterDOP(i, 1) + CORE::GEO::TOL7))) ||
            ((xfemDOP(i, 1) > (cutterDOP(i, 0) - CORE::GEO::TOL7)) &&
                (xfemDOP(i, 1) < (cutterDOP(i, 1) + CORE::GEO::TOL7)))))
      return false;

  return true;
}


/*----------------------------------------------------------------------*
 |  checks the intersection between two bounding volumes (AABB)         |
 |                                                          wirtz 08/14 |
 *----------------------------------------------------------------------*/
bool CORE::GEO::intersectionOfBVs(
    const CORE::LINALG::Matrix<3, 2>& currentBV, const CORE::LINALG::Matrix<3, 2>& queryBV)
{
  return (overlap(currentBV(0, 0), currentBV(0, 1), queryBV(0, 0), queryBV(0, 1)) and
          overlap(currentBV(1, 0), currentBV(1, 1), queryBV(1, 0), queryBV(1, 1)) and
          overlap(currentBV(2, 0), currentBV(2, 1), queryBV(2, 0), queryBV(2, 1)));

  return 0;
}


/*----------------------------------------------------------------------*
 |  checks the overlap of two intervals in one coordinate               |
 |                                                          wirtz 08/14 |
 *----------------------------------------------------------------------*/
bool CORE::GEO::overlap(double smin, double smax, double omin, double omax)
{
  return ((omax > smin - CORE::GEO::TOL7 and omin < smax + CORE::GEO::TOL7) or
          (smax > omin - CORE::GEO::TOL7 and smin < omax + CORE::GEO::TOL7));
}

FOUR_C_NAMESPACE_CLOSE
