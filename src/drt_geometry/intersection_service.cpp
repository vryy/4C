/*!----------------------------------------------------------------------
\file intersection_service.cpp

\brief collection of service methods for intersection computations


<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>
*----------------------------------------------------------------------*/


#include "element_coordtrafo.H"
#include "element_normals.H"
#include "intersection_interfacepoint.H"
#include "intersection_service.H"
#include "intersection_service_templates.H"
#include "position_array.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"



/*----------------------------------------------------------------------*
 |  ML:     computes the cross product                       u.may 08/07|
 |          of 2 vectors c = a x b                                      |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::computeCrossProduct(
    const LINALG::Matrix<3, 1>& a, const LINALG::Matrix<3, 1>& b)
{
  LINALG::Matrix<3, 1> c;

  c(0) = a(1) * b(2) - a(2) * b(1);
  c(1) = a(2) * b(0) - a(0) * b(2);
  c(2) = a(0) * b(1) - a(1) * b(0);

  return c;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if an element is CARTESIAN, LINEAR and    u.may 07/08|
 |          HIGHERORDER                                                 |
 *----------------------------------------------------------------------*/
void GEO::checkGeoType(const DRT::Element* element, const LINALG::SerialDenseMatrix& xyze_element,
    EleGeoType& eleGeoType)
{
  bool cartesian = true;
  int CartesianCount = 0;
  const int dimCoord = 3;
  const DRT::Element::DiscretizationType distype = element->Shape();
  const int eleDim = DRT::UTILS::getDimension(distype);

  if (DRT::UTILS::getOrder(distype) == 1)
    eleGeoType = LINEAR;
  else if (DRT::UTILS::getOrder(distype) == 2)
    eleGeoType = HIGHERORDER;
  else
    dserror("order of element shapefuntion is not correct");

  // check if cartesian
  if (eleDim == 3)
  {
    const std::vector<std::vector<int>> eleNodeNumbering =
        DRT::UTILS::getEleNodeNumberingSurfaces(distype);
    std::vector<Teuchos::RCP<DRT::Element>> surfaces =
        (const_cast<DRT::Element*>(element))->Surfaces();
    for (int i = 0; i < element->NumSurface(); i++)
    {
      CartesianCount = 0;
      const DRT::Element* surfaceP = surfaces[i].get();

      for (int k = 0; k < dimCoord; k++)
      {
        int nodeId = eleNodeNumbering[i][0];
        const double nodalcoord = xyze_element(k, nodeId);
        for (int j = 1; j < surfaceP->NumNode(); j++)
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
      for (int j = 1; j < element->NumNode(); j++)
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
    dserror("dimension of element is not correct");



  if (cartesian) eleGeoType = CARTESIAN;
}



/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          u.may 12/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
const std::map<int, LINALG::Matrix<3, 2>> GEO::getCurrentXAABBs(
    const DRT::Discretization& dis, const std::map<int, LINALG::Matrix<3, 1>>& currentpositions)
{
  std::map<int, LINALG::Matrix<3, 2>> currentXAABBs;
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j = 0; j < dis.NumMyColElements(); ++j)
  {
    const DRT::Element* element = dis.lColElement(j);
    const LINALG::SerialDenseMatrix xyze_element(
        GEO::getCurrentNodalPositions(element, currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkGeoType(element, xyze_element, eleGeoType);
    const LINALG::Matrix<3, 2> xaabbEle =
        GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
    currentXAABBs[element->Id()] = xaabbEle;
  }
  return currentXAABBs;
}



/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          u.may 02/09|
 | triangle list                                                        |
 *----------------------------------------------------------------------*/
const std::map<int, LINALG::Matrix<3, 2>> GEO::getTriangleXAABBs(
    const std::vector<std::vector<int>>& triangleList,
    const std::vector<GEO::InterfacePoint>& pointList)
{
  std::map<int, LINALG::Matrix<3, 2>> triangleXAABBs;
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int i = 0; i < (int)triangleList.size(); ++i)
  {
    LINALG::SerialDenseMatrix xyze_triElement(3, 3);
    for (int j = 0; j < 3; j++)
    {
      LINALG::Matrix<3, 1> node = pointList[triangleList[i][j]].getCoord();
      for (int k = 0; k < 3; k++) xyze_triElement(k, j) = node(k);
    }
    LINALG::Matrix<3, 2> xaabbEle =
        GEO::computeFastXAABB(DRT::Element::tri3, xyze_triElement, GEO::EleGeoType(LINEAR));
    triangleXAABBs[i] = xaabbEle;
  }
  return triangleXAABBs;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes 18Dops                                  u.may 12/08|
 |          (only the slabs which are not present in an XAABB)          |
 *----------------------------------------------------------------------*/
LINALG::Matrix<6, 2> GEO::computeContact18Dop(
    const DRT::Element* element, const LINALG::SerialDenseMatrix& xyze)
{
  // consider only remaining slabs
  LINALG::Matrix<6, 2> slabs;

  for (int j = 0; j < 6; j++)
    slabs(j, 0) = slabs(j, 1) =
        (GEO::Dop18Normals[j][0] * xyze(0, 0) + GEO::Dop18Normals[j][1] * xyze(1, 0) +
            GEO::Dop18Normals[j][2] * xyze(2, 0)) /
        sqrt(2.0);

  // remaining element nodes
  for (int k = 1; k < element->NumNode(); k++)
    for (int j = 0; j < 6; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      const double dcurrent =
          (GEO::Dop18Normals[j][0] * xyze(0, k) + GEO::Dop18Normals[j][1] * xyze(1, k) +
              GEO::Dop18Normals[j][2] * xyze(2, k)) /
          sqrt(2.0);
      if (dcurrent > slabs(j, 1)) slabs(j, 1) = dcurrent;
      if (dcurrent < slabs(j, 0)) slabs(j, 0) = dcurrent;
    }

  // CONTACT::CElement* selement = static_cast<CONTACT::CElement*>(element);
  // double area = selement->Area();
  double area = 0.01;
  // inflation of slabs by tolerance
  for (int j = 0; j < 6; j++)
  {
    slabs(j, 0) -= 0.1 * area;
    slabs(j, 1) += 0.1 * area;
  }

  // if higher order include proper extension
  return slabs;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if two 18DOPs intersect                   u.may 12/08| |
 *----------------------------------------------------------------------*/
bool GEO::intersectionOfKDOPs(
    const LINALG::Matrix<9, 2>& cutterDOP, const LINALG::Matrix<9, 2>& xfemDOP)
{
  // check intersection of 18 kdops
  for (int i = 0; i < 9; i++)
    if (!(((cutterDOP(i, 0) > (xfemDOP(i, 0) - GEO::TOL7)) &&
              (cutterDOP(i, 0) < (xfemDOP(i, 1) + GEO::TOL7))) ||
            ((cutterDOP(i, 1) > (xfemDOP(i, 0) - GEO::TOL7)) &&
                (cutterDOP(i, 1) < (xfemDOP(i, 1) + GEO::TOL7))) ||
            ((xfemDOP(i, 0) > (cutterDOP(i, 0) - GEO::TOL7)) &&
                (xfemDOP(i, 0) < (cutterDOP(i, 1) + GEO::TOL7))) ||
            ((xfemDOP(i, 1) > (cutterDOP(i, 0) - GEO::TOL7)) &&
                (xfemDOP(i, 1) < (cutterDOP(i, 1) + GEO::TOL7)))))
      return false;

  return true;
}


/*----------------------------------------------------------------------*
 |  checks the intersection between two bounding volumes (AABB)         |
 |                                                          wirtz 08/14 |
 *----------------------------------------------------------------------*/
bool GEO::intersectionOfBVs(
    const LINALG::Matrix<3, 2>& currentBV, const LINALG::Matrix<3, 2>& queryBV)
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
bool GEO::overlap(double smin, double smax, double omin, double omax)
{
  return ((omax > smin - GEO::TOL7 and omin < smax + GEO::TOL7) or
          (smax > omin - GEO::TOL7 and smin < omax + GEO::TOL7));
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within a given element   u.may 06/07|
 *----------------------------------------------------------------------*/
bool GEO::checkPositionWithinElement(const DRT::Element* element,
    const LINALG::SerialDenseMatrix& xyze, const LINALG::Matrix<3, 1>& x)
{
  dsassert(
      DRT::UTILS::getDimension(element->Shape()) == 3, "only valid for 3 dimensional elements");
  LINALG::Matrix<3, 1> xsi(true);
  bool nodeWithinElement = currentToVolumeElementCoordinates(element->Shape(), xyze, x, xsi);
  // printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol =
  // %20.16f\n", xsi(0),xsi(1),xsi(2), residual, TOL14);

  nodeWithinElement = checkPositionWithinElementParameterSpace(xsi, element->Shape());

  return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  RQI:    searches the nearest point on a surface          u.may 02/08|
 |          element for a given point in physical coordinates           |
 *----------------------------------------------------------------------*/
bool GEO::searchForNearestPointOnSurface(const DRT::Element* surfaceElement,
    const LINALG::SerialDenseMatrix& xyze_surfaceElement, const LINALG::Matrix<3, 1>& physCoord,
    LINALG::Matrix<2, 1>& eleCoord, LINALG::Matrix<3, 1>& normal, double& distance)
{
  CurrentToSurfaceElementCoordinates(
      surfaceElement->Shape(), xyze_surfaceElement, physCoord, eleCoord);

  const bool pointWithinElement =
      checkPositionWithinElementParameterSpace(eleCoord, surfaceElement->Shape());

  // normal vector at position xsi
  static LINALG::Matrix<3, 1> eleNormalAtXsi;
  computeNormalToSurfaceElement(
      surfaceElement->Shape(), xyze_surfaceElement, eleCoord, eleNormalAtXsi);

  LINALG::Matrix<3, 1> x_surface_phys;
  elementToCurrentCoordinates(
      surfaceElement->Shape(), xyze_surfaceElement, eleCoord, x_surface_phys);
  // normal pointing away from the surface towards physCoord
  normal.Update(1.0, physCoord, -1.0, x_surface_phys);
  // absolute distance between point and surface
  distance = normal.Norm2();


  if (fabs(distance) > GEO::TOL7)
  {
    // compute distance with sign
    const double scalarproduct = eleNormalAtXsi(0) * normal(0) + eleNormalAtXsi(1) * normal(1) +
                                 eleNormalAtXsi(2) * normal(2);
    const double teiler = eleNormalAtXsi.Norm2() * normal.Norm2();
    const double cosphi = scalarproduct / teiler;
    const double vorzeichen = cosphi / abs(cosphi);
    distance *= vorzeichen;
  }
  else
    distance = 0.0;

  return pointWithinElement;
}
