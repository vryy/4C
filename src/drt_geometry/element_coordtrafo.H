/*----------------------------------------------------------------------*/
/*! \file

\brief computes coordinate transformation between current and element
       domains (both ways)

\level 2

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#ifndef ELEMENT_COORDTRAFO_H_
#define ELEMENT_COORDTRAFO_H_


#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_gauss.H"


namespace GEO
{
  /*----------------------------------------------------------------------*
   | GM:      transforms a node in element coordinates         u.may 07/07|
   |          into current coordinates                                    |
   *----------------------------------------------------------------------*/
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static inline void elementToCurrentCoordinatesInPlaceT(
      const M& xyze, LINALG::Matrix<3, 1>& eleCoord)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<numNodes, 1> funct;

    DRT::UTILS::shape_function<DISTYPE>(eleCoord, funct);

    eleCoord.Clear();
    for (int i = 0; i < numNodes; i++)
      for (int j = 0; j < 3; j++) eleCoord(j) += xyze(j, i) * funct(i);

    return;
  }



  /*----------------------------------------------------------------------*
   | GM:      transforms a node in element coordinates       u.may 07/07|
   |          into current coordinates                                    |
   *----------------------------------------------------------------------*/
  template <class M>
  static inline void elementToCurrentCoordinatesInPlace(
      const DRT::Element::DiscretizationType distype, const M& xyze, LINALG::Matrix<3, 1>& eleCoord)
  {
    switch (distype)
    {
      case DRT::Element::line2:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::line2>(xyze, eleCoord);
        break;
      case DRT::Element::line3:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::line3>(xyze, eleCoord);
        break;
      case DRT::Element::tri3:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::tri3>(xyze, eleCoord);
        break;
      case DRT::Element::tri6:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::tri6>(xyze, eleCoord);
        break;
      case DRT::Element::quad4:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::quad4>(xyze, eleCoord);
        break;
      case DRT::Element::quad8:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::quad8>(xyze, eleCoord);
        break;
      case DRT::Element::quad9:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::quad9>(xyze, eleCoord);
        break;
      case DRT::Element::hex8:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::hex8>(xyze, eleCoord);
        break;
      case DRT::Element::hex20:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::hex20>(xyze, eleCoord);
        break;
      case DRT::Element::hex27:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::hex27>(xyze, eleCoord);
        break;
      case DRT::Element::tet4:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::tet4>(xyze, eleCoord);
        break;
      case DRT::Element::tet10:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::tet10>(xyze, eleCoord);
        break;
      case DRT::Element::wedge6:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::wedge6>(xyze, eleCoord);
        break;
      case DRT::Element::wedge15:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::wedge15>(xyze, eleCoord);
        break;
      case DRT::Element::pyramid5:
        elementToCurrentCoordinatesInPlaceT<DRT::Element::pyramid5>(xyze, eleCoord);
        break;
      default:
        std::cout << DistypeToString(distype) << std::endl;
        dserror("add your 3D distype to this switch!");
        break;
    }
    return;
  }



  /*!
  \brief transforms a point in element coordinates to a point
         in current coordinates with respect to a given element
  \param element              (in)        : element
  \param xyze                 (in)        : nodal positions of element
  \param xsi                  (in)        : element coordinates
  \param x                    (out)       : position in physical coordinates (x, y, z)
  */
  template <DRT::Element::DiscretizationType DISTYPE, class V, class M, unsigned probDim>
  static inline void elementToCurrentCoordinatesT(
      const M& xyze, const V& xsi, LINALG::Matrix<probDim, 1>& x)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<numNodes, 1> funct;

    DRT::UTILS::shape_function<DISTYPE>(xsi, funct);

    // set to zero
    x.Clear();
    for (int i = 0; i < numNodes; i++)
      for (unsigned j = 0; j < probDim; j++) x(j) += xyze(j, i) * funct(i);

    return;
  }



  /*!
  \brief transforms a point in element coordinates to a point
         in current coordinates with respect to a given element
  \param element              (in)        : element
  \param xyze                 (in)        : nodal positions of element
  \param xsi                  (in)        : element coordinates
  \param x                    (out)       : position in physical coordinates (x, y, z)
  */
  template <class V, class M, unsigned probDim>
  static inline void elementToCurrentCoordinates(const DRT::Element::DiscretizationType distype,
      const M& xyze, const V& xsi, LINALG::Matrix<probDim, 1>& x)
  {
    switch (distype)
    {
      case DRT::Element::line2:
        elementToCurrentCoordinatesT<DRT::Element::line2>(xyze, xsi, x);
        break;
      case DRT::Element::line3:
        elementToCurrentCoordinatesT<DRT::Element::line3>(xyze, xsi, x);
        break;
      case DRT::Element::tri3:
        elementToCurrentCoordinatesT<DRT::Element::tri3>(xyze, xsi, x);
        break;
      case DRT::Element::tri6:
        elementToCurrentCoordinatesT<DRT::Element::tri6>(xyze, xsi, x);
        break;
      case DRT::Element::quad4:
        elementToCurrentCoordinatesT<DRT::Element::quad4>(xyze, xsi, x);
        break;
      case DRT::Element::quad8:
        elementToCurrentCoordinatesT<DRT::Element::quad8>(xyze, xsi, x);
        break;
      case DRT::Element::quad9:
        elementToCurrentCoordinatesT<DRT::Element::quad9>(xyze, xsi, x);
        break;
      case DRT::Element::hex8:
        elementToCurrentCoordinatesT<DRT::Element::hex8>(xyze, xsi, x);
        break;
      case DRT::Element::hex20:
        elementToCurrentCoordinatesT<DRT::Element::hex20>(xyze, xsi, x);
        break;
      case DRT::Element::hex27:
        elementToCurrentCoordinatesT<DRT::Element::hex27>(xyze, xsi, x);
        break;
      case DRT::Element::tet4:
        elementToCurrentCoordinatesT<DRT::Element::tet4>(xyze, xsi, x);
        break;
      case DRT::Element::tet10:
        elementToCurrentCoordinatesT<DRT::Element::tet10>(xyze, xsi, x);
        break;
      case DRT::Element::wedge6:
        elementToCurrentCoordinatesT<DRT::Element::wedge6>(xyze, xsi, x);
        break;
      case DRT::Element::wedge15:
        elementToCurrentCoordinatesT<DRT::Element::wedge15>(xyze, xsi, x);
        break;
      case DRT::Element::pyramid5:
        elementToCurrentCoordinatesT<DRT::Element::pyramid5>(xyze, xsi, x);
        break;
      default:
        std::cout << DistypeToString(distype) << std::endl;
        dserror("add your 3D distype to this switch!");
        break;
    }
    return;
  }



  /*!
  \brief computes starting value for current to volume element computation
  \param xsi                 (in)   starting vector
  \param distype               (in)   shape of the element
  \return  true if point lies in the element, false otherwise
  */
  template <DRT::Element::DiscretizationType distype, class V>
  static inline void startingValueCurrentToElementCoords(V& xsi)
  {
    switch (distype)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        xsi.Clear();
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        xsi.PutScalar(0.3);
        break;
      }
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        xsi.Clear();
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        xsi.PutScalar(0.3);
        break;
      }
      case DRT::Element::line2:
      case DRT::Element::line3:
      {
        xsi.Clear();
        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      case DRT::Element::pyramid5:
      {
        xsi.Clear();
        break;
      }
      default:
        dserror("distype not yet implemented");
        break;
    }
    return;
  }


  /*!
  \brief computes starting value for current to volume element computation with a fast guess
  \param xyze                 (in)        : nodal positions of element
  \param x                    (in)        : position in physical coordinates (x, y, z)
  \param xsi                  (out)       : initial guess for xsi
  \param distype              (in)        : shape of the element
  \return  true if point lies in the element, false otherwise (ATTENTION: this is just an
  approximation and can be incorrect)
  */
  template <DRT::Element::DiscretizationType distype, class M1, class V3>
  static inline bool FastInitialGuess(const M1& xyze,  ///< nodal position array
      const V3& x,                                     ///< (x,y,z)
      LINALG::Matrix<3, 1>& xsi)                       ///< initial guess for xsi
  {
    switch (distype)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        // split hex element into 5 tetrahedra in order to obtain a good initial guess
        static LINALG::Matrix<4, 4> A_tet;
        static LINALG::Matrix<4, 1> b_tet;
        static LINALG::Matrix<4, 1> x_tet;


        // first tet 1 3 4 6 + fourth row filled with ones
        A_tet(0, 0) = xyze(0, 1);
        A_tet(1, 0) = xyze(1, 1);
        A_tet(2, 0) = xyze(2, 1);
        A_tet(3, 0) = 1.0;

        A_tet(0, 1) = xyze(0, 3);
        A_tet(1, 1) = xyze(1, 3);
        A_tet(2, 1) = xyze(2, 3);
        A_tet(3, 1) = 1.0;

        A_tet(0, 2) = xyze(0, 4);
        A_tet(1, 2) = xyze(1, 4);
        A_tet(2, 2) = xyze(2, 4);
        A_tet(3, 2) = 1.0;

        A_tet(0, 3) = xyze(0, 6);
        A_tet(1, 3) = xyze(1, 6);
        A_tet(2, 3) = xyze(2, 6);
        A_tet(3, 3) = 1.0;

        // same rhs for all tets
        b_tet(0) = x(0);
        b_tet(1) = x(1);
        b_tet(2) = x(2);
        b_tet(3) = 1.0;

        double det = LINALG::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) dserror("determinant is near zero %d", det);

        bool inside = true;
        for (int i = 0; i < 4; ++i)
        {
          if (std::min(x_tet(i), 1.0 - x_tet(i)) < -TOL7)
          {
            inside = false;
            break;
          }
        }

        // if inside: map xsi of tet into xsi of hex and leave
        if (inside)
        {
          xsi(0) = 1.0 - 2 * (x_tet(1) + x_tet(2));
          xsi(1) = -1.0 + 2 * (x_tet(1) + x_tet(3));
          xsi(2) = -1.0 + 2 * (x_tet(2) + x_tet(3));
          return inside;
        }


        // second tet 0 1 3 4 + fourth row filled with ones
        A_tet(0, 0) = xyze(0, 0);
        A_tet(1, 0) = xyze(1, 0);
        A_tet(2, 0) = xyze(2, 0);
        A_tet(3, 0) = 1.0;

        A_tet(0, 1) = xyze(0, 1);
        A_tet(1, 1) = xyze(1, 1);
        A_tet(2, 1) = xyze(2, 1);
        A_tet(3, 1) = 1.0;

        A_tet(0, 2) = xyze(0, 3);
        A_tet(1, 2) = xyze(1, 3);
        A_tet(2, 2) = xyze(2, 3);
        A_tet(3, 2) = 1.0;

        A_tet(0, 3) = xyze(0, 4);
        A_tet(1, 3) = xyze(1, 4);
        A_tet(2, 3) = xyze(2, 4);
        A_tet(3, 3) = 1.0;

        // same rhs for all tets
        b_tet(0) = x(0);
        b_tet(1) = x(1);
        b_tet(2) = x(2);
        b_tet(3) = 1.0;

        det = LINALG::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) dserror("determinant is near zero %d", det);

        inside = true;
        for (int i = 0; i < 4; ++i)
        {
          if (std::min(x_tet(i), 1.0 - x_tet(i)) < -TOL7)
          {
            inside = false;
            break;
          }
        }

        // if inside: map xsi of tet into xsi of hex and leave
        if (inside)
        {
          xsi(0) = -1.0 + 2 * x_tet(1);
          xsi(1) = -1.0 + 2 * x_tet(2);
          xsi(2) = -1.0 + 2 * x_tet(3);
          return inside;
        }


        // third tet 2 3 1 6 + fourth row filled with ones
        A_tet(0, 0) = xyze(0, 2);
        A_tet(1, 0) = xyze(1, 2);
        A_tet(2, 0) = xyze(2, 2);
        A_tet(3, 0) = 1.0;

        A_tet(0, 1) = xyze(0, 3);
        A_tet(1, 1) = xyze(1, 3);
        A_tet(2, 1) = xyze(2, 3);
        A_tet(3, 1) = 1.0;

        A_tet(0, 2) = xyze(0, 1);
        A_tet(1, 2) = xyze(1, 1);
        A_tet(2, 2) = xyze(2, 1);
        A_tet(3, 2) = 1.0;

        A_tet(0, 3) = xyze(0, 6);
        A_tet(1, 3) = xyze(1, 6);
        A_tet(2, 3) = xyze(2, 6);
        A_tet(3, 3) = 1.0;

        // same rhs for all tets
        b_tet(0) = x(0);
        b_tet(1) = x(1);
        b_tet(2) = x(2);
        b_tet(3) = 1.0;

        det = LINALG::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) dserror("determinant is near zero %d", det);

        inside = true;
        for (int i = 0; i < 4; ++i)
        {
          if (std::min(x_tet(i), 1.0 - x_tet(i)) < -TOL7)
          {
            inside = false;
            break;
          }
        }

        // if inside: map xsi of tet into xsi of hex and leave
        if (inside)
        {
          xsi(0) = 1.0 - 2 * x_tet(1);
          xsi(1) = 1.0 - 2 * x_tet(2);
          xsi(2) = -1.0 + 2 * x_tet(3);
          return inside;
        }


        // fourth tet 5 4 6 1 + fourth row filled with ones
        A_tet(0, 0) = xyze(0, 5);
        A_tet(1, 0) = xyze(1, 5);
        A_tet(2, 0) = xyze(2, 5);
        A_tet(3, 0) = 1.0;

        A_tet(0, 1) = xyze(0, 4);
        A_tet(1, 1) = xyze(1, 4);
        A_tet(2, 1) = xyze(2, 4);
        A_tet(3, 1) = 1.0;

        A_tet(0, 2) = xyze(0, 6);
        A_tet(1, 2) = xyze(1, 6);
        A_tet(2, 2) = xyze(2, 6);
        A_tet(3, 2) = 1.0;

        A_tet(0, 3) = xyze(0, 1);
        A_tet(1, 3) = xyze(1, 1);
        A_tet(2, 3) = xyze(2, 1);
        A_tet(3, 3) = 1.0;

        // same rhs for all tets
        b_tet(0) = x(0);
        b_tet(1) = x(1);
        b_tet(2) = x(2);
        b_tet(3) = 1.0;

        det = LINALG::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) dserror("determinant is near zero %d", det);

        inside = true;
        for (int i = 0; i < 4; ++i)
        {
          if (std::min(x_tet(i), 1.0 - x_tet(i)) < -TOL7)
          {
            inside = false;
            break;
          }
        }

        // if inside: map xsi of tet into xsi of hex and leave
        if (inside)
        {
          xsi(0) = 1.0 - 2 * x_tet(1);
          xsi(1) = -1.0 + 2 * x_tet(2);
          xsi(2) = 1.0 - 2 * x_tet(3);
          return inside;
        }


        // fifth tet 7 6 4 3 + fourth row filled with ones
        A_tet(0, 0) = xyze(0, 7);
        A_tet(1, 0) = xyze(1, 7);
        A_tet(2, 0) = xyze(2, 7);
        A_tet(3, 0) = 1.0;

        A_tet(0, 1) = xyze(0, 6);
        A_tet(1, 1) = xyze(1, 6);
        A_tet(2, 1) = xyze(2, 6);
        A_tet(3, 1) = 1.0;

        A_tet(0, 2) = xyze(0, 4);
        A_tet(1, 2) = xyze(1, 4);
        A_tet(2, 2) = xyze(2, 4);
        A_tet(3, 2) = 1.0;

        A_tet(0, 3) = xyze(0, 3);
        A_tet(1, 3) = xyze(1, 3);
        A_tet(2, 3) = xyze(2, 3);
        A_tet(3, 3) = 1.0;

        // same rhs for all tets
        b_tet(0) = x(0);
        b_tet(1) = x(1);
        b_tet(2) = x(2);
        b_tet(3) = 1.0;

        det = LINALG::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) dserror("determinant is near zero %d", det);

        inside = true;
        for (int i = 0; i < 4; ++i)
        {
          if (std::min(x_tet(i), 1.0 - x_tet(i)) < -TOL7)
          {
            inside = false;
            break;
          }
        }

        // if inside: map xsi of tet into xsi of hex and leave
        if (inside)
        {
          xsi(0) = -1.0 + 2 * x_tet(1);
          xsi(1) = 1.0 - 2 * x_tet(2);
          xsi(2) = 1.0 - 2 * x_tet(3);
          return inside;
        }

        // NOTE: showing up here is equivalent to x being outside the 5 tetrahedrons
        // which is not necessarily identical to being outside the hex element
        return false;
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        dserror("no fast initial guess implemented");
        break;
      }
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        dserror("no fast initial guess implemented");
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        dserror("no fast initial guess implemented");
        break;
      }
      case DRT::Element::line2:
      case DRT::Element::line3:
      {
        dserror("no fast initial guess implemented");
        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      case DRT::Element::pyramid5:
      {
        dserror("no fast initial guess implemented");
        break;
      }
      default:
        dserror("distype not yet implemented");
        break;
    }
    return false;
  }


  /*!
  \brief checks if a position in element coordinates lies within a certain Element parameter space
  \param eleCoord              (in)   node in element coordinates (r,s,t)
  \param distype               (in)   shape of the element
  \return  true if point lies in the element, false otherwise
  */
  template <class V>
  bool checkPositionWithinElementParameterSpace(
      const V& eleCoord, const DRT::Element::DiscretizationType distype)
  {
    switch (distype)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      case DRT::Element::line2:
      case DRT::Element::line3:
      {
        for (int i = 0; i < DRT::UTILS::getDimension(distype); i++)
          if (eleCoord(i) > (1.0 + GEO::TOL7) || eleCoord(i) < (-1) * (1.0 + GEO::TOL7))
            return false;

        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        if (eleCoord(0) > (1.0 + GEO::TOL7) || eleCoord(0) < (-1) * GEO::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(0) + 2 * GEO::TOL7) ||
            eleCoord(1) < (-1) * GEO::TOL7)  // s = 0 ... 1 -r
          return false;
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        if (eleCoord(0) > (1.0 + GEO::TOL7) || eleCoord(0) < (-1) * GEO::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 + GEO::TOL7) || eleCoord(1) < (-1) * GEO::TOL7)  // s = 0 ... 1
          return false;
        if (eleCoord(2) > (1.0 - eleCoord(0) - eleCoord(1) + 3 * GEO::TOL7) ||
            eleCoord(2) < (-1) * GEO::TOL7)  // t = 0 ... 1 - r -s
          return false;                      // TODO check tolerance
        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
        if (eleCoord(0) > (1.0 + GEO::TOL7) || eleCoord(0) < (-1) * GEO::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(0) + 2 * GEO::TOL7) ||
            eleCoord(1) < (-1) * GEO::TOL7)  // s = 0 ... 1 -r
          return false;
        if (eleCoord(2) > (1.0 + GEO::TOL7) || eleCoord(2) < (-1) * (1.0 + GEO::TOL7)) return false;
        break;
      case DRT::Element::pyramid5:
        if (eleCoord(2) > (1.0 + GEO::TOL7) || eleCoord(0) < (-1) * GEO::TOL7)  // t = 0 ... 1
          return false;
        if (eleCoord(0) > (1.0 - eleCoord(2) + GEO::TOL7) ||
            eleCoord(0) < (-1.0 + eleCoord(2) - GEO::TOL7))  // r
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(2) + GEO::TOL7) ||
            eleCoord(1) < (-1.0 + eleCoord(2) - GEO::TOL7))  // s
          return false;
        break;
      default:
        dserror("distype not yet implemented");
        break;
    }
    return true;
  }



  /*!
  \brief updates the system matrix at the corresponding element coordinates for the
         computation if a node in current coordinates lies within an element

  \tparam dim   dimension of the problem

  */
  template <DRT::Element::DiscretizationType DISTYPE, int dim, class M1, class V, class M2>
  static inline void updateAForNWE(M1& A,  ///< system matrix
      const V& xsi,                        ///< vector of element coordinates
      const M2& xyze                       ///< nodal position array (3,numNodes)
  )
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<dim, numNodes> deriv1;
    DRT::UTILS::shape_function_deriv1<DISTYPE>(xsi, deriv1);

    A.Clear();
    for (int isd = 0; isd < 3; ++isd)
      for (int jsd = 0; jsd < dim; ++jsd)
        for (int inode = 0; inode < numNodes; ++inode)
          A(isd, jsd) += xyze(isd, inode) * deriv1(jsd, inode);
  }



  /*!
  \brief updates the rhs at the corresponding element coordinates for the
         computation whether a node in current coordinates lies within an element

  */
  template <DRT::Element::DiscretizationType DISTYPE,  ///< shape of the element
      int dim, class V1, class V2, class V3, class M1>
  static inline void updateRHSForNWE(V1& b,  ///< right-hand-sid
      const V2& xsi,                         ///< vector of element coordinates
      const V3& x,                           ///< node in physical coordinates (x,y,z)
      const M1& xyze                         ///< nodal coordinates of an element with shape DISTYPE
  )
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<numNodes, 1> funct;
    funct.Clear();
    DRT::UTILS::shape_function<DISTYPE>(xsi, funct);

    for (int isd = 0; isd < 3; ++isd)
    {
      b(isd) = x(isd);
      for (int inode = 0; inode < numNodes; ++inode) b(isd) -= xyze(isd, inode) * funct(inode);
    }
  }


  /*!
  \brief transforms a point in current coordinates to a point
         in element coordinates with respect to a given element
         The nonlinear system of equation is solved with help of the Newton-method.
         Fast templated version

  \param xyze                 (in)        : element nodal positions (3,numnode)
  \param x                    (in)        : node in current coordinates \f$(x, y, z)\f$
  \param xsi                  (inout)     : node in element coordinates
  \param fastinitguess        (in)        : if true, fast (but approximate) initial guess is used
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M1, class V3>
  static inline bool currentToVolumeElementCoordinatesT(const M1& xyze,  ///< nodal position array
      const V3& x,                                                       ///< (x,y,z)
      LINALG::Matrix<3, 1>& xsi,                                         ///< (r,s,t)
      bool fastinitguess)  ///< enhanced intitial guess
  {
    // REMARK: This function seemed to deliver wrong results!
    // It seems that some matrices (especially dx) were not cleared at the beginning!
    // we fixed this at 11.09.2012
    // if you want to use this function you can remove the dserror, but please check carefully if it
    // works properly REMARK: tested this method for quite a while and seems to be fine now
    // (28.11.2012 ghamm)
    //  dserror("This function seems to deliver wrong results! Check carefully before use or use a
    //  different function, e.g. GEO::CUT::Position");
    // REMARK: Added a scaling of the function -> biggest entry is 1 and thus this function should
    // work for elements far away
    //         from the origin as well. (27.07.2015 winter & seitz)

    // number space dimensions for element
    const size_t dim = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // number of nodes of element
    const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    bool nodeWithinElement = true;
    const int maxiter = 20;  // 40;
    double residual = 1.0;

    static LINALG::Matrix<3, 3> A;
    static LINALG::Matrix<3, 1> b;
    static LINALG::Matrix<3, 1> dx;

    A.Clear();
    b.Clear();
    dx.Clear();

    // initial guess
    if (fastinitguess == false)
    {
      startingValueCurrentToElementCoords<DISTYPE>(xsi);
    }
    else
    {
      // either get good initial guess and continue or
      // leave here when outside (which is just an approximation and can be incorrect)
      nodeWithinElement = FastInitialGuess<DISTYPE>(xyze, x, xsi);
      if (nodeWithinElement == false) return nodeWithinElement;
    }

    // Scale the problem before attempting Newton iteration.
    //  This is necessary if elements are "far" away from the origin.
    // Eq that is solved:
    //  (x_GP  - N(xi) * xyze )*Row_Scaling = 0

    // Row_Scaling is a diagonal scaling matrix to transform the system
    //   to scale each row to [-1,1]

    //               |s_1  0   0 |
    // Row_Scaling = | 0  s_2  0 |
    //               | 0   0  s_3|

    LINALG::Matrix<3, 1> scale_vector;
    //  for(int i=0; i<(int)xyze.N(); ++i)
    //    for(int j=0; j<(int)xyze.M(); ++j)
    for (int i = 0; i < (int)nen; ++i)
      for (int j = 0; j < (int)dim; ++j)
        scale_vector(j) = std::max(scale_vector(j), fabs(xyze(j, i)));

    V3 x_scaled(x);
    M1 xyze_scaled(xyze);

    for (int j = 0; j < (int)dim; ++j)
    {
      scale_vector(j) = 1 / scale_vector(j);
      for (int i = 0; i < (int)nen; ++i)
      {
        xyze_scaled(j, i) *= scale_vector(j);
      }
      x_scaled(j) *= scale_vector(j);
    }

    updateRHSForNWE<DISTYPE, dim>(b, xsi, x_scaled, xyze_scaled);

    int iter = 0;
    while (residual > GEO::TOL14)
    {
      updateAForNWE<DISTYPE, dim>(A, xsi, xyze_scaled);

      double det = LINALG::gaussElimination<true, 3>(A, b, dx);

      if (fabs(det) < GEO::TOL14)
      {
        dserror("determinant is near zero %d", det);
      }

      xsi += dx;
      updateRHSForNWE<DISTYPE, dim>(b, xsi, x_scaled, xyze_scaled);

      residual = b.Norm2();
      iter++;

      if (iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8)
      {
        nodeWithinElement = false;
        break;
      }
    }

    // final check whether node is within element
    if (nodeWithinElement == true)
      nodeWithinElement = checkPositionWithinElementParameterSpace(xsi, DISTYPE);

    return nodeWithinElement;
  }



  /*! \brief computes element coordinates of point x within element
   *
   *  \param distype              (in)        : shape of the element
   *  \param xyze                 (in)        : nodal positions of element
   *  \param x                    (in)        : position in physical coordinates (x, y, z)
   *  \param xsi                  (out)       : element coordinates xsi
   *  \param fastinitguess        (in)        : if true, approximate xsi first (default: false)
   *  \return  true if point lies in the element, false otherwise
   *
   *  \note If \c fastinitguess == TRUE this is just an approximation and can be incorrect */
  template <class M1, class V3>
  bool currentToVolumeElementCoordinates(const DRT::Element::DiscretizationType distype,
      const M1& xyze, const V3& x, LINALG::Matrix<3, 1>& xsi, bool fastinitguess = false)
  {
    bool nodeWithinElement = false;

    switch (distype)
    {
      case DRT::Element::hex8:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::hex8>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::hex20:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::hex20>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::hex27:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::hex27>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::tet4:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::tet4>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::tet10:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::tet10>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::wedge6:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::wedge6>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::wedge15:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::wedge15>(xyze, x, xsi, fastinitguess);
        break;
      case DRT::Element::pyramid5:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<DRT::Element::pyramid5>(xyze, x, xsi, fastinitguess);
        break;
      default:
        std::cout << DistypeToString(distype) << std::endl;
        dserror("add your 3D distype to this switch!");
        nodeWithinElement = false;
        break;
    }

    return nodeWithinElement;
  }



  /*----------------------------------------------------------------------*
   | GM:  transforms a node in current coordinates            u.may 12/07 |
   |      into element coordinates                                        |
   *----------------------------------------------------------------------*/
  template <class M1, class V3>
  LINALG::Matrix<3, 1> currentToVolumeElementCoordinatesExact(
      const DRT::Element::DiscretizationType distype, const M1& xyze, const V3& x, const double tol)
  {
    LINALG::Matrix<3, 1> xsi;
    currentToVolumeElementCoordinates(distype, xyze, x, xsi);

    switch (distype)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        // rounding 1 and -1 to be exact for the CDT
        for (int j = 0; j < 3; j++)
        {
          if (fabs((fabs(xsi(j)) - 1.0)) < tol && xsi(j) < 0) xsi(j) = -1.0;
          if (fabs((fabs(xsi(j)) - 1.0)) < tol && xsi(j) > 0) xsi(j) = 1.0;
        }
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        // round r on 0 or 1
        // round s on 0 or 1
        // round t on 0 or 1-r-s
        if (fabs(xsi(0)) < tol) xsi(0) = 0.0;
        if (fabs(xsi(0) - 1.0) < tol) xsi(0) = 1.0;

        if (fabs(xsi(1)) < tol) xsi(1) = 0.0;
        if (fabs(xsi(1) - 1.0) < tol) xsi(1) = 1.0;

        if (fabs(xsi(2) - (1.0 - xsi(0) - xsi(1))) < 3 * tol) xsi(2) = 1.0 - xsi(0) - xsi(1);
        if (fabs(xsi(2)) < 2 * tol) xsi(2) = 0.0;

        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      {
        if (fabs(xsi(0)) < tol) xsi(0) = 0.0;
        if (fabs(xsi(0) - 1.0) < tol) xsi(0) = 1.0;

        if (fabs(xsi(1)) < tol) xsi(1) = 0.0;
        if (fabs(xsi(1) - 1.0) < tol) xsi(1) = 1.0;

        if (fabs((fabs(xsi(2)) - 1.0)) < tol && xsi(2) < 0) xsi(2) = -1.0;
        if (fabs((fabs(xsi(2)) - 1.0)) < tol && xsi(2) > 0) xsi(2) = 1.0;

        break;
      }
      case DRT::Element::pyramid5:
      {
        if (fabs(xsi(2)) < tol) xsi(2) = 0.0;
        if (fabs(xsi(2) - 1.0) < tol) xsi(2) = 1.0;

        if (fabs(xsi(0)) < -1.0 + xsi(2) - tol) xsi(0) = -1.0 + xsi(2);
        if (fabs(xsi(0)) > 1.0 - xsi(2) + tol) xsi(0) = 1.0 - xsi(2);

        if (fabs(xsi(1)) < -1.0 + xsi(2) - tol) xsi(1) = -1.0 + xsi(2);
        if (fabs(xsi(1)) > 1.0 - xsi(2) + tol) xsi(1) = 1.0 - xsi(2);

        break;
      }
      default:
        dserror("element type not implemented");
        break;
    };

    return xsi;
  }



  /*!
  \brief  updates the 3x2 Jacobian matrix for the computation of the
          surface element coordinates for a point in physical coordinates
  \param Jacobi               (in) :  Jacobi
  \param xsi                  (in) :  node in element coordinates (r, s)
  \param xyze_surfaceElement  (in) :  element nodal coordinates
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static inline void updateJacobianForMap3To2(
      LINALG::Matrix<3, 2>& Jacobi, const LINALG::Matrix<2, 1>& xsi, const M& xyze_surfaceElement)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    static LINALG::Matrix<2, numNodes> deriv1;
    DRT::UTILS::shape_function_2D_deriv1(deriv1, xsi(0), xsi(1), DISTYPE);

    Jacobi.Clear();
    for (int isd = 0; isd < 3; isd++)
      for (int jsd = 0; jsd < 2; jsd++)
        for (int inode = 0; inode < numNodes; inode++)
          Jacobi(isd, jsd) += xyze_surfaceElement(isd, inode) * deriv1(jsd, inode);

    return;
  }



  /*!
  \brief  updates the nonlinear equations for the computation of the
          surface element coordinates for a surface point in physical coordinates
  \param F                    (out):  right hand sied of Xp -X(r,s)
  \param xsi                  (in) :  node in element coordinates (r, s)
  \param x                    (in) :  node in physical coordinates
  \param xyze_surfaceElement  (in) :  element nodal coordinates
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static inline void updateFForMap3To2(LINALG::Matrix<3, 1>& F, const LINALG::Matrix<2, 1>& xsi,
      const LINALG::Matrix<3, 1>& x, const M& xyze_surfaceElement)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<numNodes, 1> funct;
    DRT::UTILS::shape_function_2D(funct, xsi(0), xsi(1), DISTYPE);

    F.Clear();
    for (int isd = 0; isd < 3; ++isd)
      for (int inode = 0; inode < numNodes; inode++)
        F(isd) += xyze_surfaceElement(isd, inode) * funct(inode);

    F -= x;
    return;
  }



  /*!
  \brief  updates the nonlinear equations for the computation of the
          surface element coordinates for a surface point in physical coordinates
  \param A                    (out):  system matrix
  \param Jacobi               (in) :  Jacobian matrix
  \param F                    (in) :  right hand sied of Xp -X(r,s)
  \param xsi                  (in) :  node in element coordinates (r, s)
  \param xyze_surfaceElement  (in) :  element nodal coordinates
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static void updateAForMap3To2(LINALG::Matrix<2, 2>& A, const LINALG::Matrix<3, 2>& Jacobi,
      const LINALG::Matrix<3, 1>& F, const LINALG::Matrix<2, 1>& xsi, const M& xyze_surfaceElement)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<3, numNodes> deriv2;
    DRT::UTILS::shape_function_2D_deriv2(deriv2, xsi(0), xsi(1), DISTYPE);

    // third order tensor 3 x 2 x 2 stored as 3x2 and 3x2 so 3x4
    static LINALG::Matrix<3, 4> tensor3order;
    tensor3order.Clear();

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          for (int inode = 0; inode < numNodes; inode++)
            tensor3order(i, k * 2 + j) += xyze_surfaceElement(i, inode) * deriv2(j, inode);


    A.Clear();
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 3; ++k)
          A(i, j) += Jacobi(k, i) * Jacobi(k, j) + F(k) * tensor3order(k, j * 2 + i);
  }



  /*!
  \brief compute element coordinates from a given point
         in the 3-dim physical space lies on a given surface element
         Template version
  \param xyze_surfaceElement      (in) :  element nodal coordinates
  \param physCoord                (in) :  physical coordinates (x,y,z)
  \param eleCoord                 (in) :  element coordinates (r,s)

   RQI
   */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static void currentToSurfaceElementCoordinatesT(const M& xyze_surfaceElement,
      const LINALG::Matrix<3, 1>& physCoord, LINALG::Matrix<2, 1>& eleCoord)
  {
    startingValueCurrentToElementCoords<DISTYPE>(eleCoord);

    const int maxiter = 20;
    int iter = 0;
    while (iter < maxiter)
    {
      iter++;

      // compute Jacobian, f and b
      static LINALG::Matrix<3, 2> Jacobi;
      static LINALG::Matrix<3, 1> F;
      updateJacobianForMap3To2<DISTYPE>(Jacobi, eleCoord, xyze_surfaceElement);
      updateFForMap3To2<DISTYPE>(F, eleCoord, physCoord, xyze_surfaceElement);
      static LINALG::Matrix<2, 1> b;
      b.Clear();

      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j) b(i) -= Jacobi(j, i) * F(j);

      const double residual = b.Norm2();
      if (residual < GEO::TOL14)
      {
        break;
      }

      // compute system matrix A
      static LINALG::Matrix<2, 2> A;
      updateAForMap3To2<DISTYPE>(A, Jacobi, F, eleCoord, xyze_surfaceElement);

      static LINALG::Matrix<2, 1> dx;
      dx = 0.0;

      double det = LINALG::gaussElimination<true, 2>(A, b, dx);

      if (fabs(det) < GEO::TOL14)
      {
        break;
      }

      eleCoord += dx;
    }

    //  std::cout << "nodeWithinElement: " << nodeWithinElement << std::endl;
    //  printf("iter = %d\n", iter);
    //  printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n",
    //  eleCoord(0),eleCoord(1), residual, GEO::TOL14); exit(1);
    return;
  }



  /*----------------------------------------------------------------------*
   |  RQI:    compute element coordinates from a given point              |
   |          in the 3-dim physical space lies on a given surface element |
   *----------------------------------------------------------------------*/
  template <class M>
  static void CurrentToSurfaceElementCoordinates(const DRT::Element::DiscretizationType distype,
      const M& xyze_surfaceElement, const LINALG::Matrix<3, 1>& physCoord,
      LINALG::Matrix<2, 1>& eleCoord)
  {
    switch (distype)
    {
      case DRT::Element::quad4:
        currentToSurfaceElementCoordinatesT<DRT::Element::quad4>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case DRT::Element::quad8:
        currentToSurfaceElementCoordinatesT<DRT::Element::quad8>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case DRT::Element::quad9:
        currentToSurfaceElementCoordinatesT<DRT::Element::quad9>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case DRT::Element::tri3:
        currentToSurfaceElementCoordinatesT<DRT::Element::tri3>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case DRT::Element::tri6:
        currentToSurfaceElementCoordinatesT<DRT::Element::tri6>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      default:
        dserror("please add your surface element type here");
        break;
    }
  }



  /*!
  \brief compute element coordinates from a given point
         in the 3-dim physical space lies on a given line element
         Template version
  \param xyze_lineElement         (in) :  element nodal coordinates
  \param physCoord                (in) :  physical coordinates (x,y,z)
  \param eleCoord                 (in) :  element coordinates (r)

   RQI
   */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  static void currentToLineElementCoordinatesT(const M& xyze_lineElement,
      const LINALG::Matrix<3, 1>& physCoord, LINALG::Matrix<1, 1>& eleCoord)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    const int maxiter = 20;
    int iter = 0;
    double residual = 1.0;

    // starting value
    startingValueCurrentToElementCoords<DISTYPE>(eleCoord);

    while (residual > GEO::TOL13)
    {
      iter++;

      // determine shapefunction, 1. and 2. derivative at current solutiom
      static LINALG::Matrix<numNodes, 1> funct;
      DRT::UTILS::shape_function_1D(funct, eleCoord(0), DISTYPE);

      static LINALG::Matrix<1, numNodes> deriv1;
      DRT::UTILS::shape_function_1D_deriv1(deriv1, eleCoord(0), DISTYPE);

      static LINALG::Matrix<1, numNodes> deriv2;
      DRT::UTILS::shape_function_1D_deriv2(deriv2, eleCoord(0), DISTYPE);

      // compute nonlinear system
      static LINALG::Matrix<3, 1> F;
      // compute first derivative of r
      static LINALG::Matrix<3, 1> F_deriv1;
      // compute first derivative of r
      static LINALG::Matrix<3, 1> F_deriv2;

      F.Clear();
      F_deriv1.Clear();
      F_deriv2.Clear();

      for (int i = 0; i < 3; i++)
        for (int inode = 0; inode < numNodes; inode++)
        {
          F(i) += xyze_lineElement(i, inode) * funct(inode);
          F_deriv1(i) += xyze_lineElement(i, inode) * deriv1(0, inode);
          F_deriv2(i) += xyze_lineElement(i, inode) * deriv2(0, inode);
        }

      F -= physCoord;

      // std::cout << "F " << F << std::endl;
      // std::cout << "F_deriv1 " << F_deriv1 << std::endl;
      // std::cout << "F_deriv2 " << F_deriv2 << std::endl;

      // determine system matrix A and rhs b
      double A = 0.0;
      double b = 0.0;
      // update system matrix A and rhs b
      for (int i = 0; i < 3; i++)
      {
        A += F(i) * F_deriv2(i) + F_deriv1(i) * F_deriv1(i);
        b += F_deriv1(i) * F(i);
      }

      if (fabs(A) < GEO::TOL14)
      {
        // printf("A is equal to zero 0");
        break;
      }
      // solve scalar linear equation  Delta_x = -b/A
      eleCoord(0) += (-1) * b / A;

      // std::cout << " A = " << A << std::endl;
      // std::cout << " b = " << b << std::endl;
      // std::cout << " x = " << eleCoord(0) << std::endl;

      if (iter >= maxiter || fabs(eleCoord(0)) > GEO::TOLPLUS8) break;

      residual = fabs(b);
    }
  }



  /*----------------------------------------------------------------------*
   |  RQI:    compute element coordinates from a given point              |
   |          in the 3-dim physical space lies on a given line element    |
   *----------------------------------------------------------------------*/
  template <class M>
  static void CurrentToLineElementCoordinates(const DRT::Element::DiscretizationType distype,
      const M& xyze_lineElement, const LINALG::Matrix<3, 1>& physCoord,
      LINALG::Matrix<1, 1>& eleCoord)
  {
    switch (distype)
    {
      case DRT::Element::line2:
        currentToLineElementCoordinatesT<DRT::Element::line2>(
            xyze_lineElement, physCoord, eleCoord);
        break;
      case DRT::Element::line3:
        currentToLineElementCoordinatesT<DRT::Element::line3>(
            xyze_lineElement, physCoord, eleCoord);
        break;
      default:
        dserror("please add your line element type here");
        break;
    }
  }



  /*----------------------------------------------------------------------*
   | GM:      transforms a point in global coordinates        schott 07/12|
   |          into local element coordinates                              |
   *----------------------------------------------------------------------*/
  // template <DRT::Element::DiscretizationType DISTYPE>
  template <DRT::Element::DiscretizationType DISTYPE>
  static bool ComputeLocalCoordinates(
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim,
          DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement>& xyze,
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim, 1>& x,
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<DISTYPE>::dim, 1>& xsi)
  {
    bool inelement = true;

    const int nsd_ = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    const int maxiter = 20;
    double residual = 1.0;

    LINALG::Matrix<nsd_, nsd_> A(true);
    LINALG::Matrix<nsd_, nsd_> A_inv(true);

    LINALG::Matrix<nsd_, 1> b(true);
    LINALG::Matrix<nsd_, 1> dx(true);

    // initialize != 0
    dx(0) = 1.0;

    LINALG::Matrix<numNodes, 1> funct(true);
    LINALG::Matrix<nsd_, numNodes> deriv1(true);


    // initial guess
    startingValueCurrentToElementCoords<DISTYPE>(xsi);

    // update rhs b= -(x(xi)-x_point)
    DRT::UTILS::shape_function<DISTYPE>(xsi, funct);

    for (int isd = 0; isd < nsd_; ++isd)
    {
      b(isd) = x(isd);
      for (int inode = 0; inode < numNodes; ++inode) b(isd) -= xyze(isd, inode) * funct(inode);
    }


    // newton loop
    int iter = 0;
    while (residual > GEO::TOL12 or dx.Norm2() > GEO::TOL12)
    {
      DRT::UTILS::shape_function_deriv1<DISTYPE>(xsi, deriv1);

      // get the jacobian
      A.Clear();
      for (int isd = 0; isd < nsd_; ++isd)
        for (int jsd = 0; jsd < nsd_; ++jsd)
          for (int inode = 0; inode < numNodes; ++inode)
            A(isd, jsd) += xyze(isd, inode) * deriv1(jsd, inode);

      const double det = A_inv.Invert(A);

      if (fabs(det) < GEO::TOL14)
      {
        dserror("determinant is near zero %d", det);
      }

      dx.Clear();
      dx.Multiply(A_inv, b);

      // update of local coordinates
      for (int i = 0; i < nsd_; i++) xsi(i, 0) += dx(i);

      // update rhs b= -(x(xi)-x_point)
      DRT::UTILS::shape_function<DISTYPE>(xsi, funct);
      for (int isd = 0; isd < nsd_; ++isd)
      {
        b(isd) = x(isd);
        for (int inode = 0; inode < numNodes; ++inode) b(isd) -= xyze(isd, inode) * funct(inode);
      }

      residual = b.Norm2();
      iter++;

      if (iter >= maxiter || xsi.Norm2() > GEO::TOLPLUS8)
      {
        inelement = false;
        break;
      }
    }

    inelement = checkPositionWithinElementParameterSpace(xsi, DISTYPE);

    return inelement;
  }



}  // namespace GEO


#endif /*ELEMENT_COORDTRAFO_H_*/
