/*----------------------------------------------------------------------*/
/*! \file

\brief computes coordinate transformation between current and element
       domains (both ways)

\level 2

*----------------------------------------------------------------------*/


#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_ELEMENT_COORDTRAFO_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_ELEMENT_COORDTRAFO_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_geometry_intersection_math.hpp"
#include "4C_linalg_gauss.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::Geo
{
  /*!
  \brief transforms a point in element coordinates to a point
         in current coordinates with respect to a given element
  \param element              (in)        : element
  \param xyze                 (in)        : nodal positions of element
  \param xsi                  (in)        : element coordinates
  \param x                    (out)       : position in physical coordinates (x, y, z)
  */
  template <FE::CellType DISTYPE, class V, class M, unsigned probDim>
  static inline void elementToCurrentCoordinatesT(
      const M& xyze, const V& xsi, Core::LinAlg::Matrix<probDim, 1>& x)
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;
    static Core::LinAlg::Matrix<numNodes, 1> funct;

    Core::FE::shape_function<DISTYPE>(xsi, funct);

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
  static inline void elementToCurrentCoordinates(
      const FE::CellType distype, const M& xyze, const V& xsi, Core::LinAlg::Matrix<probDim, 1>& x)
  {
    switch (distype)
    {
      case FE::CellType::line2:
        elementToCurrentCoordinatesT<FE::CellType::line2>(xyze, xsi, x);
        break;
      case FE::CellType::line3:
        elementToCurrentCoordinatesT<FE::CellType::line3>(xyze, xsi, x);
        break;
      case FE::CellType::tri3:
        elementToCurrentCoordinatesT<FE::CellType::tri3>(xyze, xsi, x);
        break;
      case FE::CellType::tri6:
        elementToCurrentCoordinatesT<FE::CellType::tri6>(xyze, xsi, x);
        break;
      case FE::CellType::quad4:
        elementToCurrentCoordinatesT<FE::CellType::quad4>(xyze, xsi, x);
        break;
      case FE::CellType::quad8:
        elementToCurrentCoordinatesT<FE::CellType::quad8>(xyze, xsi, x);
        break;
      case FE::CellType::quad9:
        elementToCurrentCoordinatesT<FE::CellType::quad9>(xyze, xsi, x);
        break;
      case FE::CellType::hex8:
        elementToCurrentCoordinatesT<FE::CellType::hex8>(xyze, xsi, x);
        break;
      case FE::CellType::hex20:
        elementToCurrentCoordinatesT<FE::CellType::hex20>(xyze, xsi, x);
        break;
      case FE::CellType::hex27:
        elementToCurrentCoordinatesT<FE::CellType::hex27>(xyze, xsi, x);
        break;
      case FE::CellType::tet4:
        elementToCurrentCoordinatesT<FE::CellType::tet4>(xyze, xsi, x);
        break;
      case FE::CellType::tet10:
        elementToCurrentCoordinatesT<FE::CellType::tet10>(xyze, xsi, x);
        break;
      case FE::CellType::wedge6:
        elementToCurrentCoordinatesT<FE::CellType::wedge6>(xyze, xsi, x);
        break;
      case FE::CellType::wedge15:
        elementToCurrentCoordinatesT<FE::CellType::wedge15>(xyze, xsi, x);
        break;
      case FE::CellType::pyramid5:
        elementToCurrentCoordinatesT<FE::CellType::pyramid5>(xyze, xsi, x);
        break;
      default:
        std::cout << Core::FE::CellTypeToString(distype) << std::endl;
        FOUR_C_THROW("add your 3D distype to this switch!");
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
  template <FE::CellType distype, class V>
  static inline void startingValueCurrentToElementCoords(V& xsi)
  {
    switch (distype)
    {
      case FE::CellType::hex8:
      case FE::CellType::hex20:
      case FE::CellType::hex27:
      {
        xsi.Clear();
        break;
      }
      case FE::CellType::tet4:
      case FE::CellType::tet10:
      {
        xsi.PutScalar(0.3);
        break;
      }
      case FE::CellType::quad4:
      case FE::CellType::quad8:
      case FE::CellType::quad9:
      {
        xsi.Clear();
        break;
      }
      case FE::CellType::tri3:
      case FE::CellType::tri6:
      {
        xsi.PutScalar(0.3);
        break;
      }
      case FE::CellType::line2:
      case FE::CellType::line3:
      {
        xsi.Clear();
        break;
      }
      case FE::CellType::wedge6:
      case FE::CellType::wedge15:
      case FE::CellType::pyramid5:
      {
        xsi.Clear();
        break;
      }
      default:
        FOUR_C_THROW("distype not yet implemented");
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
  template <FE::CellType distype, class M1, class V3>
  static inline bool FastInitialGuess(const M1& xyze,  ///< nodal position array
      const V3& x,                                     ///< (x,y,z)
      Core::LinAlg::Matrix<3, 1>& xsi)                 ///< initial guess for xsi
  {
    switch (distype)
    {
      case FE::CellType::hex8:
      case FE::CellType::hex20:
      case FE::CellType::hex27:
      {
        // split hex element into 5 tetrahedra in order to obtain a good initial guess
        static Core::LinAlg::Matrix<4, 4> A_tet;
        static Core::LinAlg::Matrix<4, 1> b_tet;
        static Core::LinAlg::Matrix<4, 1> x_tet;


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

        double det = Core::LinAlg::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) FOUR_C_THROW("determinant is near zero %d", det);

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

        det = Core::LinAlg::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) FOUR_C_THROW("determinant is near zero %d", det);

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

        det = Core::LinAlg::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) FOUR_C_THROW("determinant is near zero %d", det);

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

        det = Core::LinAlg::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) FOUR_C_THROW("determinant is near zero %d", det);

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

        det = Core::LinAlg::gaussElimination<true, 4>(A_tet, b_tet, x_tet);
        if (fabs(det) < 1E-14) FOUR_C_THROW("determinant is near zero %d", det);

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
      case FE::CellType::tet4:
      case FE::CellType::tet10:
      {
        FOUR_C_THROW("no fast initial guess implemented");
        break;
      }
      case FE::CellType::quad4:
      case FE::CellType::quad8:
      case FE::CellType::quad9:
      {
        FOUR_C_THROW("no fast initial guess implemented");
        break;
      }
      case FE::CellType::tri3:
      case FE::CellType::tri6:
      {
        FOUR_C_THROW("no fast initial guess implemented");
        break;
      }
      case FE::CellType::line2:
      case FE::CellType::line3:
      {
        FOUR_C_THROW("no fast initial guess implemented");
        break;
      }
      case FE::CellType::wedge6:
      case FE::CellType::wedge15:
      case FE::CellType::pyramid5:
      {
        FOUR_C_THROW("no fast initial guess implemented");
        break;
      }
      default:
        FOUR_C_THROW("distype not yet implemented");
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
  bool checkPositionWithinElementParameterSpace(const V& eleCoord, const FE::CellType distype)
  {
    switch (distype)
    {
      case FE::CellType::hex8:
      case FE::CellType::hex20:
      case FE::CellType::hex27:
      case FE::CellType::quad4:
      case FE::CellType::quad8:
      case FE::CellType::quad9:
      case FE::CellType::line2:
      case FE::CellType::line3:
      {
        for (int i = 0; i < Core::FE::getDimension(distype); i++)
          if (eleCoord(i) > (1.0 + Core::Geo::TOL7) || eleCoord(i) < (-1) * (1.0 + Core::Geo::TOL7))
            return false;

        break;
      }
      case FE::CellType::tri3:
      case FE::CellType::tri6:
      {
        if (eleCoord(0) > (1.0 + Core::Geo::TOL7) ||
            eleCoord(0) < (-1) * Core::Geo::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(0) + 2 * Core::Geo::TOL7) ||
            eleCoord(1) < (-1) * Core::Geo::TOL7)  // s = 0 ... 1 -r
          return false;
        break;
      }
      case FE::CellType::tet4:
      case FE::CellType::tet10:
      {
        if (eleCoord(0) > (1.0 + Core::Geo::TOL7) ||
            eleCoord(0) < (-1) * Core::Geo::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 + Core::Geo::TOL7) ||
            eleCoord(1) < (-1) * Core::Geo::TOL7)  // s = 0 ... 1
          return false;
        if (eleCoord(2) > (1.0 - eleCoord(0) - eleCoord(1) + 3 * Core::Geo::TOL7) ||
            eleCoord(2) < (-1) * Core::Geo::TOL7)  // t = 0 ... 1 - r -s
          return false;                            // TODO check tolerance
        break;
      }
      case FE::CellType::wedge6:
      case FE::CellType::wedge15:
        if (eleCoord(0) > (1.0 + Core::Geo::TOL7) ||
            eleCoord(0) < (-1) * Core::Geo::TOL7)  // r = 0 ... 1
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(0) + 2 * Core::Geo::TOL7) ||
            eleCoord(1) < (-1) * Core::Geo::TOL7)  // s = 0 ... 1 -r
          return false;
        if (eleCoord(2) > (1.0 + Core::Geo::TOL7) || eleCoord(2) < (-1) * (1.0 + Core::Geo::TOL7))
          return false;
        break;
      case FE::CellType::pyramid5:
        if (eleCoord(2) > (1.0 + Core::Geo::TOL7) ||
            eleCoord(0) < (-1) * Core::Geo::TOL7)  // t = 0 ... 1
          return false;
        if (eleCoord(0) > (1.0 - eleCoord(2) + Core::Geo::TOL7) ||
            eleCoord(0) < (-1.0 + eleCoord(2) - Core::Geo::TOL7))  // r
          return false;
        if (eleCoord(1) > (1.0 - eleCoord(2) + Core::Geo::TOL7) ||
            eleCoord(1) < (-1.0 + eleCoord(2) - Core::Geo::TOL7))  // s
          return false;
        break;
      default:
        FOUR_C_THROW("distype not yet implemented");
        break;
    }
    return true;
  }

  /*!
  \brief updates the system matrix at the corresponding element coordinates for the
         computation if a node in current coordinates lies within an element

  \tparam dim   dimension of the problem

  */
  template <FE::CellType DISTYPE, int dim, class M1, class V, class M2>
  static inline void updateAForNWE(M1& A,  ///< system matrix
      const V& xsi,                        ///< vector of element coordinates
      const M2& xyze                       ///< nodal position array (3,numNodes)
  )
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;
    static Core::LinAlg::Matrix<dim, numNodes> deriv1;
    Core::FE::shape_function_deriv1<DISTYPE>(xsi, deriv1);

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
  template <FE::CellType DISTYPE,  ///< shape of the element
      int dim, class V1, class V2, class V3, class M1>
  static inline void updateRHSForNWE(V1& b,  ///< right-hand-sid
      const V2& xsi,                         ///< vector of element coordinates
      const V3& x,                           ///< node in physical coordinates (x,y,z)
      const M1& xyze                         ///< nodal coordinates of an element with shape DISTYPE
  )
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;
    static Core::LinAlg::Matrix<numNodes, 1> funct;
    funct.Clear();
    Core::FE::shape_function<DISTYPE>(xsi, funct);

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
  template <FE::CellType DISTYPE, class M1, class V3>
  static inline bool currentToVolumeElementCoordinatesT(const M1& xyze,  ///< nodal position array
      const V3& x,                                                       ///< (x,y,z)
      Core::LinAlg::Matrix<3, 1>& xsi,                                   ///< (r,s,t)
      bool fastinitguess)  ///< enhanced intitial guess
  {
    // REMARK: This function seemed to deliver wrong results!
    // It seems that some matrices (especially dx) were not cleared at the beginning!
    // we fixed this at 11.09.2012
    // if you want to use this function you can remove the FOUR_C_THROW, but please check carefully
    // if it works properly REMARK: tested this method for quite a while and seems to be fine now
    // (28.11.2012 ghamm)
    //  FOUR_C_THROW("This function seems to deliver wrong results! Check carefully before use or
    //  use a different function, e.g. Core::Geo::Cut::Position");
    // REMARK: Added a scaling of the function -> biggest entry is 1 and thus this function should
    // work for elements far away from the origin as well. (27.07.2015 winter & seitz)

    // number space dimensions for element
    const size_t dim = Core::FE::dim<DISTYPE>;

    // number of nodes of element
    const size_t nen = Core::FE::num_nodes<DISTYPE>;
    bool nodeWithinElement = true;
    const int maxiter = 20;  // 40;
    double residual = 1.0;

    static Core::LinAlg::Matrix<3, 3> A;
    static Core::LinAlg::Matrix<3, 1> b;
    static Core::LinAlg::Matrix<3, 1> dx;

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

    Core::LinAlg::Matrix<3, 1> scale_vector;
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
    while (residual > Core::Geo::TOL14)
    {
      updateAForNWE<DISTYPE, dim>(A, xsi, xyze_scaled);

      double det = Core::LinAlg::gaussElimination<true, 3>(A, b, dx);

      if (fabs(det) < Core::Geo::TOL14)
      {
        FOUR_C_THROW("determinant is near zero %d", det);
      }

      xsi += dx;
      updateRHSForNWE<DISTYPE, dim>(b, xsi, x_scaled, xyze_scaled);

      residual = b.Norm2();
      iter++;

      if (iter >= maxiter || xsi.Norm1() > Core::Geo::TOLPLUS8)
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
  bool currentToVolumeElementCoordinates(const FE::CellType distype, const M1& xyze, const V3& x,
      Core::LinAlg::Matrix<3, 1>& xsi, bool fastinitguess = false)
  {
    bool nodeWithinElement = false;

    switch (distype)
    {
      case FE::CellType::hex8:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::hex8>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::hex20:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::hex20>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::hex27:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::hex27>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::tet4:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::tet4>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::tet10:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::tet10>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::wedge6:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::wedge6>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::wedge15:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::wedge15>(xyze, x, xsi, fastinitguess);
        break;
      case FE::CellType::pyramid5:
        nodeWithinElement =
            currentToVolumeElementCoordinatesT<FE::CellType::pyramid5>(xyze, x, xsi, fastinitguess);
        break;
      default:
        std::cout << Core::FE::CellTypeToString(distype) << std::endl;
        FOUR_C_THROW("add your 3D distype to this switch!");
        nodeWithinElement = false;
        break;
    }

    return nodeWithinElement;
  }

  /*!
  \brief  updates the 3x2 Jacobian matrix for the computation of the
          surface element coordinates for a point in physical coordinates
  \param Jacobi               (in) :  Jacobi
  \param xsi                  (in) :  node in element coordinates (r, s)
  \param xyze_surfaceElement  (in) :  element nodal coordinates
  */
  template <FE::CellType DISTYPE, class M>
  static inline void updateJacobianForMap3To2(Core::LinAlg::Matrix<3, 2>& Jacobi,
      const Core::LinAlg::Matrix<2, 1>& xsi, const M& xyze_surfaceElement)
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;

    static Core::LinAlg::Matrix<2, numNodes> deriv1;
    Core::FE::shape_function_2D_deriv1(deriv1, xsi(0), xsi(1), DISTYPE);

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
  template <FE::CellType DISTYPE, class M>
  static inline void updateFForMap3To2(Core::LinAlg::Matrix<3, 1>& F,
      const Core::LinAlg::Matrix<2, 1>& xsi, const Core::LinAlg::Matrix<3, 1>& x,
      const M& xyze_surfaceElement)
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;
    static Core::LinAlg::Matrix<numNodes, 1> funct;
    Core::FE::shape_function_2D(funct, xsi(0), xsi(1), DISTYPE);

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
  template <FE::CellType DISTYPE, class M>
  static void updateAForMap3To2(Core::LinAlg::Matrix<2, 2>& A,
      const Core::LinAlg::Matrix<3, 2>& Jacobi, const Core::LinAlg::Matrix<3, 1>& F,
      const Core::LinAlg::Matrix<2, 1>& xsi, const M& xyze_surfaceElement)
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;
    static Core::LinAlg::Matrix<3, numNodes> deriv2;
    Core::FE::shape_function_2D_deriv2(deriv2, xsi(0), xsi(1), DISTYPE);

    // third order tensor 3 x 2 x 2 stored as 3x2 and 3x2 so 3x4
    static Core::LinAlg::Matrix<3, 4> tensor3order;
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
  \brief compute element coordinates from a given point in the 3-dim physical space
         lies on a given surface element
  \param xyze_surfaceElement      (in) :  element nodal coordinates
  \param physCoord                (in) :  physical coordinates (x,y,z)
  \param eleCoord                 (in) :  element coordinates (r,s)
   */
  template <FE::CellType DISTYPE, class M>
  static void currentToSurfaceElementCoordinatesT(const M& xyze_surfaceElement,
      const Core::LinAlg::Matrix<3, 1>& physCoord, Core::LinAlg::Matrix<2, 1>& eleCoord)
  {
    startingValueCurrentToElementCoords<DISTYPE>(eleCoord);

    const int maxiter = 20;
    int iter = 0;
    while (iter < maxiter)
    {
      iter++;

      // compute Jacobian, f and b
      static Core::LinAlg::Matrix<3, 2> Jacobi;
      static Core::LinAlg::Matrix<3, 1> F;
      updateJacobianForMap3To2<DISTYPE>(Jacobi, eleCoord, xyze_surfaceElement);
      updateFForMap3To2<DISTYPE>(F, eleCoord, physCoord, xyze_surfaceElement);
      static Core::LinAlg::Matrix<2, 1> b;
      b.Clear();

      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j) b(i) -= Jacobi(j, i) * F(j);

      const double residual = b.Norm2();
      if (residual < Core::Geo::TOL14)
      {
        break;
      }

      // compute system matrix A
      static Core::LinAlg::Matrix<2, 2> A;
      updateAForMap3To2<DISTYPE>(A, Jacobi, F, eleCoord, xyze_surfaceElement);

      static Core::LinAlg::Matrix<2, 1> dx;
      dx = 0.0;

      double det = Core::LinAlg::gaussElimination<true, 2>(A, b, dx);

      if (fabs(det) < Core::Geo::TOL14)
      {
        break;
      }

      eleCoord += dx;
    }

    return;
  }

  /*!
  \brief compute element coordinates from a given point in the 3-dim physical space
         lies on a given line element
  \param xyze_lineElement         (in) :  element nodal coordinates
  \param physCoord                (in) :  physical coordinates (x,y,z)
  \param eleCoord                 (in) :  element coordinates (r)
   */
  template <class M>
  static void CurrentToSurfaceElementCoordinates(const FE::CellType distype,
      const M& xyze_surfaceElement, const Core::LinAlg::Matrix<3, 1>& physCoord,
      Core::LinAlg::Matrix<2, 1>& eleCoord)
  {
    switch (distype)
    {
      case FE::CellType::quad4:
        currentToSurfaceElementCoordinatesT<FE::CellType::quad4>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case FE::CellType::quad8:
        currentToSurfaceElementCoordinatesT<FE::CellType::quad8>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case FE::CellType::quad9:
        currentToSurfaceElementCoordinatesT<FE::CellType::quad9>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case FE::CellType::tri3:
        currentToSurfaceElementCoordinatesT<FE::CellType::tri3>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      case FE::CellType::tri6:
        currentToSurfaceElementCoordinatesT<FE::CellType::tri6>(
            xyze_surfaceElement, physCoord, eleCoord);
        break;
      default:
        FOUR_C_THROW("please add your surface element type here");
        break;
    }
  }

  /*!
  \brief compute element coordinates from a given point
         in the 3-dim physical space lies on a given line element
  \param xyze_lineElement         (in) :  element nodal coordinates
  \param physCoord                (in) :  physical coordinates (x,y,z)
  \param eleCoord                 (in) :  element coordinates (r)
   */
  template <FE::CellType DISTYPE, class M>
  static void currentToLineElementCoordinatesT(const M& xyze_lineElement,
      const Core::LinAlg::Matrix<3, 1>& physCoord, Core::LinAlg::Matrix<1, 1>& eleCoord)
  {
    const int numNodes = Core::FE::num_nodes<DISTYPE>;

    const int maxiter = 20;
    int iter = 0;
    double residual = 1.0;

    // starting value
    startingValueCurrentToElementCoords<DISTYPE>(eleCoord);

    while (residual > Core::Geo::TOL13)
    {
      iter++;

      // determine shapefunction, 1. and 2. derivative at current solutiom
      static Core::LinAlg::Matrix<numNodes, 1> funct;
      Core::FE::shape_function_1D(funct, eleCoord(0), DISTYPE);

      static Core::LinAlg::Matrix<1, numNodes> deriv1;
      Core::FE::shape_function_1D_deriv1(deriv1, eleCoord(0), DISTYPE);

      static Core::LinAlg::Matrix<1, numNodes> deriv2;
      Core::FE::shape_function_1D_deriv2(deriv2, eleCoord(0), DISTYPE);

      // compute nonlinear system
      static Core::LinAlg::Matrix<3, 1> F;
      // compute first derivative of r
      static Core::LinAlg::Matrix<3, 1> F_deriv1;
      // compute first derivative of r
      static Core::LinAlg::Matrix<3, 1> F_deriv2;

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

      // determine system matrix A and rhs b
      double A = 0.0;
      double b = 0.0;
      // update system matrix A and rhs b
      for (int i = 0; i < 3; i++)
      {
        A += F(i) * F_deriv2(i) + F_deriv1(i) * F_deriv1(i);
        b += F_deriv1(i) * F(i);
      }

      if (fabs(A) < Core::Geo::TOL14)
      {
        // printf("A is equal to zero 0");
        break;
      }
      // solve scalar linear equation  Delta_x = -b/A
      eleCoord(0) += (-1) * b / A;

      if (iter >= maxiter || fabs(eleCoord(0)) > Core::Geo::TOLPLUS8) break;

      residual = fabs(b);
    }
  }

  /*!
  \brief compute element coordinates from a given point in the 3-dim physical space lies on a given
  line element
  */
  template <class M>
  static void CurrentToLineElementCoordinates(const FE::CellType distype, const M& xyze_lineElement,
      const Core::LinAlg::Matrix<3, 1>& physCoord, Core::LinAlg::Matrix<1, 1>& eleCoord)
  {
    switch (distype)
    {
      case FE::CellType::line2:
        currentToLineElementCoordinatesT<FE::CellType::line2>(
            xyze_lineElement, physCoord, eleCoord);
        break;
      case FE::CellType::line3:
        currentToLineElementCoordinatesT<FE::CellType::line3>(
            xyze_lineElement, physCoord, eleCoord);
        break;
      default:
        FOUR_C_THROW("please add your line element type here");
        break;
    }
  }

  /*!
  \brief transforms a point in global coordinates into local element coordinates

  \param xyze (in) : global coordinates of the element nodes
  \param x (out) : global coordinates of the point
  \param xsi (in) : local element coordinates of the point (parameter space)
  */
  template <Core::FE::CellType DISTYPE>
  static bool ComputeLocalCoordinates(
      Core::LinAlg::Matrix<Core::FE::dim<DISTYPE>, Core::FE::num_nodes<DISTYPE>>& xyze,
      Core::LinAlg::Matrix<Core::FE::dim<DISTYPE>, 1>& x,
      Core::LinAlg::Matrix<Core::FE::dim<DISTYPE>, 1>& xsi)
  {
    bool inelement = true;

    const int numDim = Core::FE::dim<DISTYPE>;
    const int numNodes = Core::FE::num_nodes<DISTYPE>;

    const int maxiter = 20;
    double residual = 1.0;

    Core::LinAlg::Matrix<numDim, numDim> A(true);
    Core::LinAlg::Matrix<numDim, numDim> A_inv(true);

    Core::LinAlg::Matrix<numDim, 1> b(true);
    Core::LinAlg::Matrix<numDim, 1> dx(true);

    // initialize != 0
    dx(0) = 1.0;

    Core::LinAlg::Matrix<numNodes, 1> funct(true);
    Core::LinAlg::Matrix<numDim, numNodes> deriv1(true);

    // initial guess
    startingValueCurrentToElementCoords<DISTYPE>(xsi);

    // update rhs b= -(x(xi)-x_point)
    Core::FE::shape_function<DISTYPE>(xsi, funct);

    // newton loop
    int iter = 0;
    while (residual > Core::Geo::TOL12 or dx.Norm2() > Core::Geo::TOL12)
    {
      Core::FE::shape_function_deriv1<DISTYPE>(xsi, deriv1);

      // get the jacobian
      A.Clear();
      for (int isd = 0; isd < numDim; ++isd)
        for (int jsd = 0; jsd < numDim; ++jsd)
          for (int inode = 0; inode < numNodes; ++inode)
            A(isd, jsd) += xyze(isd, inode) * deriv1(jsd, inode);

      const double det = A_inv.Invert(A);

      if (fabs(det) < Core::Geo::TOL14)
      {
        FOUR_C_THROW("determinant is near zero %d", det);
      }

      // update rhs b= -(x(xi)-x_point)
      Core::FE::shape_function<DISTYPE>(xsi, funct);
      for (int isd = 0; isd < numDim; ++isd)
      {
        b(isd) = x(isd);
        for (int inode = 0; inode < numNodes; ++inode) b(isd) -= xyze(isd, inode) * funct(inode);
      }

      // update of local coordinates
      dx.Clear();
      dx.Multiply(A_inv, b);
      for (int i = 0; i < numDim; i++) xsi(i, 0) += dx(i);

      residual = b.Norm2();
      iter++;

      if (iter >= maxiter || xsi.Norm2() > Core::Geo::TOLPLUS8)
      {
        inelement = false;
        break;
      }
    }

    inelement = checkPositionWithinElementParameterSpace(xsi, DISTYPE);

    return inelement;
  }

}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
