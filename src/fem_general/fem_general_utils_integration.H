/*----------------------------------------------------------------------*/
/*! \file

\brief structs to store integration points and weights together

\level 0

*----------------------------------------------------------------------*/
#ifndef FEM_GENERAL_UTILS_INTEGRATION_H
#define FEM_GENERAL_UTILS_INTEGRATION_H

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

#include "lib_dserror.H"
#include "lib_element.H"

namespace DRT
{
  namespace UTILS
  {
    // forward declarations
    struct IntegrationPoints1D;
    struct IntegrationPoints2D;
    struct IntegrationPoints3D;


    //! supported 3d integration rules
    //!
    //! integration rules can be distinguished into open and closed rules
    //! hereby, open rules are rules where all integration points are inside the domain
    //! and none on the boundary of the integrated domain
    enum class GaussRule3D
    {
      undefined,    ///< use this to initialize a gaussrule, if you don't know the value,
                    ///< yet!
      hex_1point,   ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
      hex_8point,   ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
      hex_18point,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 in xy-plane, 3 in z-direction
                    ///< (open)
      hex_27point,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 (open)
      hex_64point,
      hex_125point,
      hex_216point,
      hex_343point,
      hex_512point,
      hex_729point,
      hex_1000point,
      tet_1point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
      tet_4point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (open)
      tet_4point_gauss_radau,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (closed)
      tet_5point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
      tet_10point,             ///< KEAST3           , DEG.OF PRECISION 3 (closed)
      tet_11point,             ///< KEAST5           , DEG.OF PRECISION 4 (open)
      tet_15point,             ///< KEAST6           , DEG.OF PRECISION 5 (closed)
      tet_24point,             ///< KEAST7           , DEG.OF PRECISION 6 (open)
      tet_64point_peano,       ///< GAUSS-LOBATTO    , DEG.OF PRECISION 6 (open)   (Peano, 1982)
      tet_45point,             ///< KEAST9           , DEG.OF PRECISION 8 (open)
      tet_125point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 9 (open)   (Peano, 1982)
      tet_343point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 13 (open)  (Peano, 1982)
      tet_729point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 16 (open)  (Peano, 1982)
      wedge_1point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      wedge_6point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      wedge_9point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      pyramid_1point,          ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      pyramid_8point           ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    };

    inline std::string GaussRuleToString(GaussRule3D gauss_rule)
    {
      switch (gauss_rule)
      {
        case GaussRule3D::undefined:
          return "undefined";
        case GaussRule3D::hex_1point:
          return "hex_1point";
        case GaussRule3D::hex_8point:
          return "hex_8point";
        case GaussRule3D::hex_18point:
          return "hex_18point";
        case GaussRule3D::hex_27point:
          return "hex_27point";
        case GaussRule3D::hex_64point:
          return "hex_64point";
        case GaussRule3D::hex_125point:
          return "hex_125point";
        case GaussRule3D::hex_216point:
          return "hex_216point";
        case GaussRule3D::hex_343point:
          return "hex_343point";
        case GaussRule3D::hex_512point:
          return "hex_512point";
        case GaussRule3D::hex_729point:
          return "hex_729point";
        case GaussRule3D::hex_1000point:
          return "hex_1000point";
        case GaussRule3D::tet_1point:
          return "tet_1point";
        case GaussRule3D::tet_4point:
          return "tet_4point";
        case GaussRule3D::tet_4point_gauss_radau:
          return "tet_4point_gauss_radau";
        case GaussRule3D::tet_5point:
          return "tet_5point";
        case GaussRule3D::tet_10point:
          return "tet_10point";
        case GaussRule3D::tet_11point:
          return "tet_11point";
        case GaussRule3D::tet_15point:
          return "tet_15point";
        case GaussRule3D::tet_24point:
          return "tet_24point";
        case GaussRule3D::tet_64point_peano:
          return "tet_64point_peano";
        case GaussRule3D::tet_45point:
          return "tet_45point";
        case GaussRule3D::tet_125point_peano:
          return "tet_125point_peano";
        case GaussRule3D::tet_343point_peano:
          return "tet_343point_peano";
        case GaussRule3D::tet_729point_peano:
          return "tet_729point_peano";
        case GaussRule3D::wedge_1point:
          return "wedge_1point";
        case GaussRule3D::wedge_6point:
          return "wedge_6point";
        case GaussRule3D::wedge_9point:
          return "wedge_9point";
        case GaussRule3D::pyramid_1point:
          return "pyramid_1point";
        case GaussRule3D::pyramid_8point:
          return "pyramid_8point";
      }
      dserror("Your GaussRule3D (%d) has not beed added into the list above",
          static_cast<int>(gauss_rule));
      return "unknown";
    }


    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule3D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case GaussRule3D::hex_1point:
          open = true;
          break;
        case GaussRule3D::hex_8point:
          open = true;
          break;
        case GaussRule3D::hex_27point:
          open = true;
          break;
        case GaussRule3D::hex_64point:
          open = true;
          break;
        case GaussRule3D::hex_125point:
          open = true;
          break;
        case GaussRule3D::hex_216point:
          open = true;
          break;
        case GaussRule3D::hex_343point:
          open = true;
          break;
        case GaussRule3D::hex_512point:
          open = true;
          break;
        case GaussRule3D::hex_729point:
          open = true;
          break;
        case GaussRule3D::hex_1000point:
          open = true;
          break;
        case GaussRule3D::tet_1point:
          open = true;
          break;
        case GaussRule3D::tet_4point:
          open = true;
          break;
        case GaussRule3D::tet_4point_gauss_radau:
          open = false;
          break;
        case GaussRule3D::tet_5point:
          open = true;
          break;
        case GaussRule3D::tet_10point:
          open = false;
          break;
        case GaussRule3D::tet_11point:
          open = true;
          break;
        case GaussRule3D::tet_15point:
          open = false;
          break;
        case GaussRule3D::tet_24point:
          open = true;
          break;
        case GaussRule3D::tet_64point_peano:
          open = true;
          break;
        case GaussRule3D::tet_45point:
          open = true;
          break;
        case GaussRule3D::tet_125point_peano:
          open = true;
          break;
        case GaussRule3D::tet_343point_peano:
          open = true;
          break;
        case GaussRule3D::tet_729point_peano:
          open = true;
          break;
        case GaussRule3D::wedge_1point:
          open = true;
          break;
        case GaussRule3D::wedge_6point:
          open = true;
          break;
        case GaussRule3D::wedge_9point:
          open = true;
          break;
        case GaussRule3D::pyramid_1point:
          open = true;
          break;
        case GaussRule3D::pyramid_8point:
          open = true;
          break;
        default:
          std::cout << "Integration rule not categorized" << std::endl;
          std::exit(1);
      }

      return open;
    }

    /*!
     * @brief Evaluates the Gauss rule given the number of Gauss points
     *
     * @note Only Gauss-Legendre integration rules are returned.
     *
     * @tparam distype (in) : Discretization type
     * @param numgp (in): number of Gauss points
     * @return GaussRule1D / GaussRule2D / GaussRule3D depending on the dimension of the
     * discretization
     */
    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_hex_v<distype>, int> = 0>
    GaussRule3D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule3D::hex_1point;
        case 8:
          return GaussRule3D::hex_8point;
        case 18:
          return GaussRule3D::hex_18point;
        case 27:
          return GaussRule3D::hex_27point;
        case 64:
          return GaussRule3D::hex_64point;
        case 125:
          return GaussRule3D::hex_125point;
        case 216:
          return GaussRule3D::hex_216point;
        case 343:
          return GaussRule3D::hex_343point;
        case 512:
          return GaussRule3D::hex_512point;
        case 729:
          return GaussRule3D::hex_729point;
        case 1000:
          return GaussRule3D::hex_1000point;
        default:
          dserror("I don't know the GaussRule3D for hex elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule3D::undefined;
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_tet_v<distype>, int> = 0>
    GaussRule3D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule3D::tet_1point;
        case 4:
          return GaussRule3D::tet_4point;
        case 11:
          return GaussRule3D::tet_11point;
        case 24:
          return GaussRule3D::tet_24point;
        case 45:
          return GaussRule3D::tet_45point;
        default:
          dserror("I don't know an open GaussRule3D for tet elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule3D::undefined;
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_wedge_v<distype>, int> = 0>
    GaussRule3D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule3D::wedge_1point;
        case 6:
          return GaussRule3D::wedge_6point;
        case 9:
          return GaussRule3D::wedge_9point;
        default:
          dserror("I don't know the GaussRule3D for wedge elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule3D::undefined;
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_pyramid_v<distype>, int> = 0>
    GaussRule3D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule3D::pyramid_1point;
        case 8:
          return GaussRule3D::pyramid_8point;
        default:
          dserror("I don't know the GaussRule3D for pyramid elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule3D::undefined;
    }

    /// get total number of GPs of the 3D Gauss rule at compile time
    template <GaussRule3D rule>
    struct GaussRule3DToNumGaussPoints
    {
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::undefined>
    {
      static const unsigned value = 0;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_8point>
    {
      static const unsigned value = 8;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_18point>
    {
      static const unsigned value = 18;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_27point>
    {
      static const unsigned value = 27;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_64point>
    {
      static const unsigned value = 64;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_125point>
    {
      static const unsigned value = 125;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_216point>
    {
      static const unsigned value = 216;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_343point>
    {
      static const unsigned value = 343;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_512point>
    {
      static const unsigned value = 512;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_729point>
    {
      static const unsigned value = 729;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_1000point>
    {
      static const unsigned value = 1000;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_4point>
    {
      static const unsigned value = 4;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_4point_gauss_radau>
    {
      static const unsigned value = 4;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_5point>
    {
      static const unsigned value = 5;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_10point>
    {
      static const unsigned value = 10;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_11point>
    {
      static const unsigned value = 11;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_15point>
    {
      static const unsigned value = 15;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_24point>
    {
      static const unsigned value = 24;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_64point_peano>
    {
      static const unsigned value = 64;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_45point>
    {
      static const unsigned value = 45;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_125point_peano>
    {
      static const unsigned value = 125;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_343point_peano>
    {
      static const unsigned value = 343;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_729point_peano>
    {
      static const unsigned value = 729;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_6point>
    {
      static const unsigned value = 6;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_9point>
    {
      static const unsigned value = 9;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::pyramid_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<GaussRule3D::pyramid_8point>
    {
      static const unsigned value = 8;
    };

    //! supported 2d integration rules
    //!
    //! integration rules can be distinguished into open and closed rules
    //! hereby, open rules are rules where all integration points are inside the domain
    //! and none on the boundary of the integrated domain
    enum class GaussRule2D
    {
      undefined,            ///< use this to initialize a gaussrule, if you don't know the value,
                            ///< yet!
      quad_1point,          ///< degree of precision: 1  (open)
      quad_4point,          ///< degree of precision: 3  (open)
      quad_6point,          ///< degree of precision: 3/5  (open) (2 by 3 gauss points for quad6
                            ///< element
      quad_9point,          ///< degree of precision: 5  (open)
      quad_16point,         ///< degree of precision: 7  (open)
      quad_25point,         ///< degree of precision: 9  (open)
      quad_lobatto25point,  ///< degree of precision: 7  (closed)
      quad_36point,         ///< degree of precision: 11 (open)
      quad_49point,         ///< degree of precision: 13 (open)
      quad_64point,         ///< degree of precision: 15 (open)
      quad_81point,         ///< degree of precision: 17 (open)
      quad_100point,        ///< degree of precision: 19 (open)
      quad_256point,        ///< degree of precision: 31 (open)
      quad_400point,        ///< degree of precision: 39 (open)
      quad_1024point,       ///< degree of precision: 63 (open)
      tri_1point,           ///< degree of precision: 1  (open)
      tri_3point,  ///< degree of precision: 2  (open)   (Hughes, The Finite Element Method)
      tri_3point_gauss_radau,  ///< degree of precision: 2  (closed) (Hughes, The Finite
                               ///< Element Method)
      tri_4point,   ///< degree of precision: 3  (open)   (Hughes, The Finite Element Method)
      tri_6point,   ///< degree of precision: 4  (open)   (Hughes, The Finite Element Method)
      tri_7point,   ///< degree of precision: 5  (open)   (Gilbert Strang, George Fix, An
                    ///< Analysis of the Finite Element Method, Cambridge, 1973,ISBN:
                    ///< 096140888X,LC: TA335.S77.)
      tri_12point,  ///< degree of precision: 6  (open)   (Hughes, The Finite Element
                    ///< Method)
      tri_16point,  ///< degree of precision: 8  (open)   (Linbo Zhang, Tao Cui and Hui Liu:
                    ///< A SET OF SYMMETRIC QUADRATURE RULES ON TRIANGLES
                    ///<                                   AND TETRAHEDRA, Journal of
                    ///<                                   Computational Mathematics,
                    ///<                                   Vol.27, No.1, 2009, 89-96.)
      tri_37point,  ///< degree of precision: 13 (open)   (rule from ACM TOMS algorithm
                    ///< #706)
      tri_64point   ///< degree of precision: 15 (open)   (essentially a product of two 8
                    ///< point 1D Gauss-Legendre rules)
    };



    inline std::string GaussRuleToString(GaussRule2D gauss_rule)
    {
      switch (gauss_rule)
      {
        case GaussRule2D::undefined:
          return "undefined";
        case GaussRule2D::quad_1point:
          return "quad_1point";
        case GaussRule2D::quad_4point:
          return "quad_4point";
        case GaussRule2D::quad_6point:
          return "quad_6point";
        case GaussRule2D::quad_9point:
          return "quad_9point";
        case GaussRule2D::quad_16point:
          return "quad_16point";
        case GaussRule2D::quad_25point:
          return "quad_25point";
        case GaussRule2D::quad_lobatto25point:
          return "quad_lobatto25point";
        case GaussRule2D::quad_36point:
          return "quad_36point";
        case GaussRule2D::quad_49point:
          return "quad_49point";
        case GaussRule2D::quad_64point:
          return "quad_64point";
        case GaussRule2D::quad_81point:
          return "quad_81point";
        case GaussRule2D::quad_100point:
          return "quad_100point";
        case GaussRule2D::quad_256point:
          return "quad_256point";
        case GaussRule2D::quad_400point:
          return "quad_400point";
        case GaussRule2D::quad_1024point:
          return "quad_1024point";
        case GaussRule2D::tri_1point:
          return "tri_1point";
        case GaussRule2D::tri_3point:
          return "tri_3point";
        case GaussRule2D::tri_3point_gauss_radau:
          return "tri_3point_gauss_radau";
        case GaussRule2D::tri_4point:
          return "tri_4point";
        case GaussRule2D::tri_6point:
          return "tri_6point";
        case GaussRule2D::tri_7point:
          return "tri_7point";
        case GaussRule2D::tri_12point:
          return "tri_12point";
        case GaussRule2D::tri_16point:
          return "tri_16point";
        case GaussRule2D::tri_37point:
          return "tri_37point";
        case GaussRule2D::tri_64point:
          return "tri_64point";
      }
      dserror("Your GaussRule2D (%d) has not beed added into the list above",
          static_cast<int>(gauss_rule));
      return "unknown";
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_quad_v<distype>, int> = 0>
    GaussRule2D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule2D::quad_1point;
        case 4:
          return GaussRule2D::quad_4point;
        case 6:
          return GaussRule2D::quad_6point;
        case 9:
          return GaussRule2D::quad_9point;
        case 16:
          return GaussRule2D::quad_16point;
        case 25:
          return GaussRule2D::quad_25point;
        case 49:
          return GaussRule2D::quad_49point;
        case 64:
          return GaussRule2D::quad_64point;
        case 81:
          return GaussRule2D::quad_81point;
        case 100:
          return GaussRule2D::quad_100point;
        case 256:
          return GaussRule2D::quad_256point;
        case 400:
          return GaussRule2D::quad_400point;
        case 1024:
          return GaussRule2D::quad_1024point;
        default:
          dserror("I don't know the GaussRule2D for quad elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule2D::undefined;
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_tri_v<distype>, int> = 0>
    GaussRule2D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule2D::tri_1point;
        case 3:
          return GaussRule2D::tri_3point;
        case 4:
          return GaussRule2D::tri_4point;
        case 6:
          return GaussRule2D::tri_6point;
        case 7:
          return GaussRule2D::tri_7point;
        case 12:
          return GaussRule2D::tri_12point;
        case 16:
          return GaussRule2D::tri_16point;
        case 37:
          return GaussRule2D::tri_37point;
        case 64:
          return GaussRule2D::tri_64point;
        default:
          dserror("I don't know the GaussRule2D for tri elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule2D::undefined;
    }

    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule2D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case GaussRule2D::quad_1point:
          open = true;
          break;
        case GaussRule2D::quad_4point:
          open = true;
          break;
        case GaussRule2D::quad_6point:
          open = true;
          break;
        case GaussRule2D::quad_9point:
          open = true;
          break;
        case GaussRule2D::quad_16point:
          open = true;
          break;
        case GaussRule2D::quad_25point:
          open = true;
          break;
        case GaussRule2D::quad_lobatto25point:
          open = false;
          break;
        case GaussRule2D::quad_36point:
          open = true;
          break;
        case GaussRule2D::quad_49point:
          open = true;
          break;
        case GaussRule2D::quad_64point:
          open = true;
          break;
        case GaussRule2D::quad_81point:
          open = true;
          break;
        case GaussRule2D::quad_100point:
          open = true;
          break;
        case GaussRule2D::quad_256point:
          open = true;
          break;
        case GaussRule2D::quad_400point:
          open = true;
          break;
        case GaussRule2D::quad_1024point:
          open = true;
          break;
        case GaussRule2D::tri_1point:
          open = true;
          break;
        case GaussRule2D::tri_3point:
          open = true;
          break;
        case GaussRule2D::tri_3point_gauss_radau:
          open = false;
          break;
        case GaussRule2D::tri_6point:
          open = true;
          break;
        case GaussRule2D::tri_7point:
          open = true;
          break;
        case GaussRule2D::tri_12point:
          open = true;
          break;
        case GaussRule2D::tri_16point:
          open = true;
          break;
        case GaussRule2D::tri_37point:
          open = true;
          break;
        case GaussRule2D::tri_64point:
          open = true;
          break;
        default:
          std::cout << "Integration rule not categorized" << std::endl;
          std::exit(1);
      }

      return open;
    }

    //! supported 1d integration rules
    enum class GaussRule1D
    {
      undefined,            ///< use this to initialize a gaussrule, if you don't know the value,
                            ///< yet!
      line_1point,          ///< degree of precision: 1  (open)
      line_2point,          ///< degree of precision: 3  (open)
      line_3point,          ///< degree of precision: 5  (open)
      line_4point,          ///< degree of precision: 7  (open)
      line_5point,          ///< degree of precision: 9  (open)
      line_6point,          ///< degree of precision: 11 (open)
      line_7point,          ///< degree of precision: 13 (open)
      line_8point,          ///< degree of precision: 15 (open)
      line_9point,          ///< degree of precision: 17 (open)
      line_10point,         ///< degree of precision: 19 (open)
      line_16point,         ///< degree of precision: 31 (open)
      line_20point,         ///< degree of precision: 39 (open)
      line_32point,         ///< degree of precision: 63 (open)
      line_50point,         ///< degree of precision: 99 (open)
      line_lobatto2point,   ///< degree of precision: 1  (closed)
      line_lobatto3point,   ///< degree of precision: 3  (closed)
      line_lobatto4point,   ///< degree of precision: 5  (closed)
      line_lobatto5point,   ///< degree of precision: 7  (closed)
      line_lobatto6point,   ///< degree of precision: 9  (closed)
      line_lobatto7point,   ///< degree of precision: 11 (closed)
      line_lobatto8point,   ///< degree of precision: 13 (closed)
      line_lobatto9point,   ///< degree of precision: 15 (closed)
      line_lobatto10point,  ///< degree of precision: 17 (closed)
      line_lobatto11point,  ///< degree of precision: 19 (closed)
      line_lobatto12point,  ///< degree of precision: 21 (closed)
      line_lobatto13point,  ///< degree of precision: 23 (closed)
      line_lobatto14point   ///< degree of precision: 25 (closed)
    };

    inline std::string GaussRuleToString(GaussRule1D gauss_rule)
    {
      switch (gauss_rule)
      {
        case GaussRule1D::undefined:
          return "undefined";
        case GaussRule1D::line_1point:
          return "line_1point";
        case GaussRule1D::line_2point:
          return "line_2point";
        case GaussRule1D::line_3point:
          return "line_3point";
        case GaussRule1D::line_4point:
          return "line_4point";
        case GaussRule1D::line_5point:
          return "line_5point";
        case GaussRule1D::line_6point:
          return "line_6point";
        case GaussRule1D::line_7point:
          return "line_7point";
        case GaussRule1D::line_8point:
          return "line_8point";
        case GaussRule1D::line_9point:
          return "line_9point";
        case GaussRule1D::line_10point:
          return "line_10point";
        case GaussRule1D::line_16point:
          return "line_16point";
        case GaussRule1D::line_20point:
          return "line_20point";
        case GaussRule1D::line_32point:
          return "line_32point";
        case GaussRule1D::line_50point:
          return "line_50point";
        case GaussRule1D::line_lobatto2point:
          return "line_lobatto2point";
        case GaussRule1D::line_lobatto3point:
          return "line_lobatto3point";
        case GaussRule1D::line_lobatto4point:
          return "line_lobatto4point";
        case GaussRule1D::line_lobatto5point:
          return "line_lobatto5point";
        case GaussRule1D::line_lobatto6point:
          return "line_lobatto6point";
        case GaussRule1D::line_lobatto7point:
          return "line_lobatto7point";
        case GaussRule1D::line_lobatto8point:
          return "line_lobatto8point";
        case GaussRule1D::line_lobatto9point:
          return "line_lobatto9point";
        case GaussRule1D::line_lobatto10point:
          return "line_lobatto10point";
        case GaussRule1D::line_lobatto11point:
          return "line_lobatto11point";
        case GaussRule1D::line_lobatto12point:
          return "line_lobatto12point";
        case GaussRule1D::line_lobatto13point:
          return "line_lobatto13point";
        case GaussRule1D::line_lobatto14point:
          return "line_lobatto14point";
      }
      dserror("Your GaussRule1D (%d) has not beed added into the list above",
          static_cast<int>(gauss_rule));
      return "unknown";
    }

    template <DRT::Element::DiscretizationType distype,
        std::enable_if_t<DRT::is_line_v<distype>, int> = 0>
    GaussRule1D NumGaussPointsToGaussRule(unsigned numgp)
    {
      switch (numgp)
      {
        case 1:
          return GaussRule1D::line_1point;
        case 2:
          return GaussRule1D::line_2point;
        case 3:
          return GaussRule1D::line_3point;
        case 4:
          return GaussRule1D::line_4point;
        case 5:
          return GaussRule1D::line_5point;
        case 6:
          return GaussRule1D::line_6point;
        case 7:
          return GaussRule1D::line_7point;
        case 8:
          return GaussRule1D::line_8point;
        case 9:
          return GaussRule1D::line_9point;
        case 10:
          return GaussRule1D::line_10point;
        case 16:
          return GaussRule1D::line_16point;
        case 20:
          return GaussRule1D::line_20point;
        case 32:
          return GaussRule1D::line_32point;
        case 50:
          return GaussRule1D::line_50point;
        default:
          dserror("I don't know the GaussRule1D for line elements with %d Gauss points", numgp);
      }

      // to make the compiler happy
      return GaussRule1D::undefined;
    }

    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule1D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case GaussRule1D::line_1point:
          open = true;
          break;
        case GaussRule1D::line_2point:
          open = true;
          break;
        case GaussRule1D::line_3point:
          open = true;
          break;
        case GaussRule1D::line_4point:
          open = true;
          break;
        case GaussRule1D::line_5point:
          open = true;
          break;
        case GaussRule1D::line_6point:
          open = true;
          break;
        case GaussRule1D::line_7point:
          open = true;
          break;
        case GaussRule1D::line_8point:
          open = true;
          break;
        case GaussRule1D::line_9point:
          open = true;
          break;
        case GaussRule1D::line_10point:
          open = true;
          break;
        case GaussRule1D::line_16point:
          open = true;
          break;
        case GaussRule1D::line_20point:
          open = true;
          break;
        case GaussRule1D::line_32point:
          open = true;
          break;
        case GaussRule1D::line_50point:
          open = true;
          break;
        case GaussRule1D::line_lobatto2point:
          open = false;
          break;
        case GaussRule1D::line_lobatto3point:
          open = false;
          break;
        case GaussRule1D::line_lobatto4point:
          open = false;
          break;
        case GaussRule1D::line_lobatto5point:
          open = false;
          break;
        case GaussRule1D::line_lobatto6point:
          open = false;
          break;
        case GaussRule1D::line_lobatto7point:
          open = false;
          break;
        case GaussRule1D::line_lobatto8point:
          open = false;
          break;
        case GaussRule1D::line_lobatto9point:
          open = false;
          break;
        case GaussRule1D::line_lobatto10point:
          open = false;
          break;
        case GaussRule1D::line_lobatto11point:
          open = false;
          break;
        case GaussRule1D::line_lobatto12point:
          open = false;
          break;
        case GaussRule1D::line_lobatto13point:
          open = false;
          break;
        case GaussRule1D::line_lobatto14point:
          open = false;
          break;
        default:
          std::cout << "Integration rule not categorized" << std::endl;
          std::exit(1);
      }

      return open;
    }

    /*!
     \brief a point with a position and an associated weight

     */
    template <int dim, class Vec>
    class IntegrationPoint
    {
     public:
      //! Standard constructor
      explicit IntegrationPoint(const Vec& pos, const double weight)
          : position_(pos), weight_(weight)
      {
        return;
      }

      //! copy constructor
      explicit IntegrationPoint(const IntegrationPoint& other)
          : position_(other.position_), weight_(other.weight_)
      {
        return;
      }

      ~IntegrationPoint();

      //! return string representation of this class
      std::string toString() const
      {
        std::stringstream s;
        s << "IntPoint, Position: " << position_ << ", weight: " << weight_;
        return s.str();
      }

      //! get integration point position
      Vec Position() const { return position_; }

      //! get integration point weight
      double Weight() const { return weight_; }

     private:
      //! 1d, 2d or 3d position of integration point in element coordinates
      const Vec position_;

      //! integration weight
      const double weight_;

      //! Empty constructor
      explicit IntegrationPoint() : position_(), weight_(0.0) { return; }
    };  // class IntegrationPoint


    /*!
    \brief integration parameters
    \author Axel Gerstenberger

    In this structure the coordinates and weights used by gauss integration are stored.

    \note tetrahedral integration rule weights are stored such that after integration over the
    parameter space, we get the real volume of 0.1666..., NOT the unit volume of 1.0. Integration
    over the hex rules gives an volume of 8.0.
    */
    struct IntegrationPoints3D
    {
      /// constructor
      explicit IntegrationPoints3D(const GaussRule3D intrule  ///< 3d rule
      );

     private:
      static const int max_nquad = 1000;  ///< size of c array
      DRT::UTILS::GaussRule3D intrule_;   ///< associated integration rule name
     public:
      int nquad;                 ///< number of gausspoints
      double qxg[max_nquad][3];  ///< coordinates
      double qwgt[max_nquad];    ///< weights

      /// gauss point coordinates
      const double* Point(int point) const
      {
        dsassert(point < max_nquad, "Index out of range");
        return qxg[point];
      };

      [[nodiscard]] int NumPoints() const { return nquad; }

      /// get used integration rule
      DRT::UTILS::GaussRule3D GetIntRule() const { return intrule_; };
    };

    /*!
    \brief integration parameters
    \author Axel Gerstenberger

    In this structure the coordinates and weights used by gauss integration are stored.

    \note triangular integration rule weights are stored such that after integration over the
    parameter space, we get the real area of 0.5, NOT the unit area of 1.0. Integration over the
    quad rules gives an area of 4.0.
    */
    struct IntegrationPoints2D
    {
      /// constructor
      explicit IntegrationPoints2D(const GaussRule2D intrule  ///< 2d rule
      );

     private:
      static const int max_nquad = 1024;  ///< size of c array
      DRT::UTILS::GaussRule2D intrule_;   ///< associated integration rule name
     public:
      int nquad;                 ///< number of gausspoints
      double qxg[max_nquad][2];  ///< coordinates
      double qwgt[max_nquad];    ///< weights

      /// gauss point coordinates
      const double* Point(int point) const
      {
        dsassert(point < max_nquad, "Index out of range");
        return qxg[point];
      };

      [[nodiscard]] int NumPoints() const { return nquad; }

      /// get used integration rule
      DRT::UTILS::GaussRule2D GetIntRule() const { return intrule_; };
    };

    /*!
    \brief integration parameters

    \author gerstenberger \date 06/07

    In this structure the coordinates and weights used by numerical integration
    are stored.

        */
    struct IntegrationPoints1D
    {
      /// constructor
      explicit IntegrationPoints1D(const GaussRule1D intrule  ///< 1d rule
      );

     private:
      static const int max_nquad = 50;   ///< size of c array
      DRT::UTILS::GaussRule1D intrule_;  ///< associated integration rule name
     public:
      int nquad;                 ///< number of integration points
      double qxg[max_nquad][1];  ///< coordinates
      double qwgt[max_nquad];    ///< weights

      /// gauss point coordinates
      const double* Point(int point) const
      {
        dsassert(point < max_nquad, "Index out of range");
        return qxg[point];
      };

      [[nodiscard]] int NumPoints() const { return nquad; }

      /// get used integration rule
      DRT::UTILS::GaussRule1D GetIntRule() const { return intrule_; };
    };


    namespace DETAIL
    {
      template <unsigned dim>
      struct DimToIntegrationPoints
      {
      };

      template <>
      struct DimToIntegrationPoints<1>
      {
        using type = DRT::UTILS::IntegrationPoints1D;
      };

      template <>
      struct DimToIntegrationPoints<2>
      {
        using type = DRT::UTILS::IntegrationPoints2D;
      };

      template <>
      struct DimToIntegrationPoints<3>
      {
        using type = DRT::UTILS::IntegrationPoints3D;
      };
    }  // namespace DETAIL

    template <unsigned dim>
    using IntegrationPoints = typename DETAIL::DimToIntegrationPoints<dim>::type;



    /*!
    \brief integration points and weights
    \author Georg Bauer
    \date 12/2008

    This class unifies the access of integration points and corresponding weights
    w.r.t. to the number of space dimensions.
    For each number of space dimensions (i.e., 1,2 and 3) a specialization of this class
    is defined. This avoids performing a switch over the number of space dimensions
    during runtime!
    The already existing structs IntegrationPoints3D,... are reused in the specializations
    in form of member variables of the class.

    */
    template <const int NSD>
    class IntPointsAndWeights
    {
     public:
      /// Constructor
      template <class G>
      explicit IntPointsAndWeights(const G gaussrule);
    };


    //! specialization for 3D
    template <>
    class IntPointsAndWeights<3>
    {
     public:
      /// Constructor for 3D specialization
      explicit IntPointsAndWeights(const GaussRule3D gaussrule)
          : intpoints_(IntegrationPoints3D(gaussrule))
      {
        return;
      };

      //! get reference to the integration points and weights
      const IntegrationPoints3D& IP() const { return intpoints_; };

      inline unsigned NumPoints() const { return IP().nquad; };
      inline unsigned NumDimension() const { return 3; }
      inline const double* Point(int point) const { return IP().qxg[point]; };
      inline double Weight(int point) const { return IP().qwgt[point]; };

     private:
      //! integration points and weights for 3D
      IntegrationPoints3D intpoints_;
    };


    //! specialization for 2D
    template <>
    class IntPointsAndWeights<2>
    {
     public:
      /// Constructor for 2D specialization
      explicit IntPointsAndWeights(const GaussRule2D gaussrule)
          : intpoints_(IntegrationPoints2D(gaussrule))
      {
        return;
      };

      //! get reference to the integration points and weights
      const IntegrationPoints2D& IP() const { return intpoints_; };

      inline unsigned NumPoints() const { return IP().nquad; };
      inline unsigned NumDimension() const { return 2; }
      inline const double* Point(int point) const { return IP().qxg[point]; };
      inline double Weight(int point) const { return IP().qwgt[point]; };

     private:
      //! integration points and weights for 3D
      IntegrationPoints2D intpoints_;
    };


    //! specialization for 1D
    template <>
    class IntPointsAndWeights<1>
    {
     public:
      /// Constructor for 1D specialization
      explicit IntPointsAndWeights(const GaussRule1D gaussrule)
          : intpoints_(IntegrationPoints1D(gaussrule))
      {
        return;
      };

      //! get reference to the integration points and weights
      const IntegrationPoints1D& IP() const { return intpoints_; };

      inline unsigned NumPoints() const { return IP().nquad; };
      inline unsigned NumDimension() const { return 1; }
      inline const double* Point(int point) const { return IP().qxg[point]; };
      inline double Weight(int point) const { return IP().qwgt[point]; };

     private:
      //! integration points and weights for 3D
      IntegrationPoints1D intpoints_;
    };


  }  // namespace UTILS
}  // namespace DRT

#endif  // FEM_GENERAL_UTILS_INTEGRATION_H
