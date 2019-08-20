/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations

\level 3

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef INTERSECTION_SERVICE_H_
#define INTERSECTION_SERVICE_H_


#include "geo_utils.H"
#include "intersection_math.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace GEO
{
  class InterfacePoint;


  const double Dop18Normals[6][3] = {{1.0, 1.0, 0.0},  //  x  y
      {1.0, 0.0, 1.0},                                 //  x  z
      {0.0, 1.0, 1.0},                                 //  y  z
      {1.0, 0.0, -1.0},                                //  x -z
      {-1.0, 1.0, 0.0},                                // -x  y
      {0.0, -1.0, 1.0}};                               // -y  z


  /*!
  \brief compares two points

  \param point1       (in)    : first point
  \param point2       (in)    : second point
  \return true if both points equal each other, false otherwise
  */
  template <int length>
  bool comparePoints(
      const LINALG::Matrix<length, 1>& point1, const LINALG::Matrix<length, 1>& point2)
  {
    for (int i = 0; i < length; i++)
      if (fabs(point1(i) - point2(i)) > TOL7)
      {
        return false;
      }
    return true;
  }


  const std::map<int, LINALG::Matrix<3, 2>> getCurrentXAABBs(
      const DRT::Discretization& dis, const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);


  const std::map<int, LINALG::Matrix<3, 2>> getTriangleXAABBs(
      const std::vector<std::vector<int>>& triangleList,
      const std::vector<GEO::InterfacePoint>& pointList);


  /*!
  \brief Computes a rough overestimating extended
         axis-aligned bounding box for an element (XAABB)
  \param element        (in)  element
  \param xyze           (in)  nodal position array (3,numnodes)
  \return extended axis-aligned bounding box  (XAABB) for an element
   */
  LINALG::Matrix<6, 2> computeContact18Dop(
      const DRT::Element* element, const LINALG::SerialDenseMatrix& xyze);


  /*!
  \brief checks if two 18 Dops are intersecting (note : for efficiency it only checks slabs
         which are not present for XAABBs)

  \param cutterDOP   (in)         : DOP of the cutting element
  \param xfemDOP     (in)         : DOP of the xfem element
  \return true if the DOP's intersect or false otherwise
   */
  bool intersectionOfKDOPs(
      const LINALG::Matrix<9, 2>& cutterDOP, const LINALG::Matrix<9, 2>& xfemDOP);

  /*!
  \brief checks the intersection between two bounding volumes (AABB)
  \param currentBV   (in)         : AABB of the current element
  \param queryBV     (in)         : AABB of the query element
  \return true if the AABB's intersect or false otherwise
   */
  bool intersectionOfBVs(
      const LINALG::Matrix<3, 2>& currentBV, const LINALG::Matrix<3, 2>& queryBV);

  /*!
  \brief checks the overlap of two intervals in one coordinate
  \param smin     (in)         : minimum value of the current interval
  \param smax     (in)         : maximum value of the current interval
  \param omin     (in)         : minimum value of the query interval
  \param omax     (in)         : maximum value of the query interval
  \return true if the intervals's overlap or false otherwise
   */
  bool overlap(double smin, double smax, double omin, double omax);


  /*!
  \brief checks if an element is Cartesian, linear or higherorder

  \param element        (in)         : element
  \param xyze_element   (in)         : coordinates of the element
  \param eleGeoType     (out)        : element geometric type CARTESIAN LINEAR or HIGHERORDER
   */
  void checkGeoType(const DRT::Element* element, const LINALG::SerialDenseMatrix& xyze_element,
      EleGeoType& eleGeoType);

  /*!
  \brief checks if an element is Cartesian, linear or higherorder

  \tparam DISTYPE       (in)         : shape of element
  \param element        (in)         : element
  \param xyze_element   (in)         : coordinates of the element
  \return eleGeoType    (out)        : element geometric type CARTESIAN LINEAR or HIGHERORDER
   */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  EleGeoType checkGeoTypeT(const DRT::Element* element, const M& xyze_element)
  {
    if (element->Shape() != DISTYPE) dserror("shape mismatch!");

    EleGeoType eleGeoType = HIGHERORDER;
    bool cartesian = true;
    int CartesianCount = 0;
    const int dimCoord = 3;
    const int eleDim = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    if (DRT::UTILS::DisTypeToEdgeOrder<DISTYPE>::order == 1)
      eleGeoType = LINEAR;
    else if (DRT::UTILS::DisTypeToEdgeOrder<DISTYPE>::order == 2)
      eleGeoType = HIGHERORDER;
    else
      dserror("order of element shapefunction is not correct");

    // check if cartesian
    if (eleDim == 3)
    {
      const std::vector<std::vector<int>> eleNodeNumbering =
          DRT::UTILS::getEleNodeNumberingSurfaces(DISTYPE);
      const std::vector<Teuchos::RCP<DRT::Element>> surfaces =
          (const_cast<DRT::Element*>(element))->Surfaces();
      for (int i = 0; i < element->NumSurface(); i++)
      {
        CartesianCount = 0;
        const Teuchos::RCP<DRT::Element>& surfaceP = surfaces[i];

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

    return eleGeoType;
  }

  /*!
  \brief checks if a position in current coordinates lies within a certain element

  \param element              (in)        : element
  \param xyze                 (in)        : element nodal positions (3,numnode)
  \param x                    (in)        : node in current coordinates
  */
  bool checkPositionWithinElement(const DRT::Element* element,
      const LINALG::SerialDenseMatrix& xyze, const LINALG::Matrix<3, 1>& x);

  /*!
  \brief  finds the nearest point for a given point on a surface element
          the element coordinates of this point are stored in eleCoord, the normal
          point from the given point to the nearest point on the surface,
          distance returns the absolut length of the normal.
          If the computed point lies outside the element , the distance which is retuned equals -1
          and the normal is set to zero.
  \param surfaceElement        (in) :   surface element
  \param xyze_surfaceElement   (in) :   nodal position array (3,numnode)
  \param physCoord             (in) :   node in physical coordinates (x, y, z)
  \param eleCoord              (out):   node in element coordinates (r, s)
  \param normal                (out):   normal point from a given point to the surface point
  \param distance              (out):   length of the normal
  */
  bool searchForNearestPointOnSurface(const DRT::Element* surfaceElement,
      const LINALG::SerialDenseMatrix& xyze_surfaceElement, const LINALG::Matrix<3, 1>& physCoord,
      LINALG::Matrix<2, 1>& eleCoord, LINALG::Matrix<3, 1>& normal, double& distance);

}  // namespace GEO


#endif  // INTERSECTION_MATH_H_
