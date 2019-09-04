/*----------------------------------------------------------------------*/
/*! \file

\brief structs to store integration points and weights together

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/
#ifndef UTILS_INTEGRATION_H
#define UTILS_INTEGRATION_H

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

#include "../drt_lib/drt_dserror.H"

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
    enum GaussRule3D
    {
      intrule3D_undefined,  ///< use this to initialize a gaussrule, if you don't know the value,
                            ///< yet!
      intrule_hex_1point,   ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
      intrule_hex_8point,   ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
      intrule_hex_18point,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 in xy-plane, 3 in z-direction
                            ///< (open)
      intrule_hex_27point,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 (open)
      intrule_hex_64point,
      intrule_hex_125point,
      intrule_hex_216point,
      intrule_hex_343point,
      intrule_hex_512point,
      intrule_hex_729point,
      intrule_hex_1000point,
      intrule_tet_1point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
      intrule_tet_4point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (open)
      intrule_tet_4point_gauss_radau,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (closed)
      intrule_tet_5point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
      intrule_tet_10point,             ///< KEAST3           , DEG.OF PRECISION 3 (closed)
      intrule_tet_11point,             ///< KEAST5           , DEG.OF PRECISION 4 (open)
      intrule_tet_15point,             ///< KEAST6           , DEG.OF PRECISION 5 (closed)
      intrule_tet_24point,             ///< KEAST7           , DEG.OF PRECISION 6 (open)
      intrule_tet_64point_peano,   ///< GAUSS-LOBATTO    , DEG.OF PRECISION 6 (open)   (Peano, 1982)
      intrule_tet_45point,         ///< KEAST9           , DEG.OF PRECISION 8 (open)
      intrule_tet_125point_peano,  ///< GAUSS-LOBATTO    , DEG.OF PRECISION 9 (open)   (Peano, 1982)
      intrule_tet_343point_peano,  ///< GAUSS-LOBATTO    , DEG.OF PRECISION 13 (open)  (Peano, 1982)
      intrule_tet_729point_peano,  ///< GAUSS-LOBATTO    , DEG.OF PRECISION 16 (open)  (Peano, 1982)
      intrule_wedge_1point,        ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      intrule_wedge_6point,        ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      intrule_wedge_9point,        ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      intrule_pyramid_1point,      ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
      intrule_pyramid_8point       ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    };

    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule3D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case intrule_hex_1point:
          open = true;
          break;
        case intrule_hex_8point:
          open = true;
          break;
        case intrule_hex_27point:
          open = true;
          break;
        case intrule_hex_64point:
          open = true;
          break;
        case intrule_hex_125point:
          open = true;
          break;
        case intrule_hex_216point:
          open = true;
          break;
        case intrule_hex_343point:
          open = true;
          break;
        case intrule_hex_512point:
          open = true;
          break;
        case intrule_hex_729point:
          open = true;
          break;
        case intrule_hex_1000point:
          open = true;
          break;
        case intrule_tet_1point:
          open = true;
          break;
        case intrule_tet_4point:
          open = true;
          break;
        case intrule_tet_4point_gauss_radau:
          open = false;
          break;
        case intrule_tet_5point:
          open = true;
          break;
        case intrule_tet_10point:
          open = false;
          break;
        case intrule_tet_11point:
          open = true;
          break;
        case intrule_tet_15point:
          open = false;
          break;
        case intrule_tet_24point:
          open = true;
          break;
        case intrule_tet_64point_peano:
          open = true;
          break;
        case intrule_tet_45point:
          open = true;
          break;
        case intrule_tet_125point_peano:
          open = true;
          break;
        case intrule_tet_343point_peano:
          open = true;
          break;
        case intrule_tet_729point_peano:
          open = true;
          break;
        case intrule_wedge_1point:
          open = true;
          break;
        case intrule_wedge_6point:
          open = true;
          break;
        case intrule_wedge_9point:
          open = true;
          break;
        case intrule_pyramid_1point:
          open = true;
          break;
        case intrule_pyramid_8point:
          open = true;
          break;
        default:
          std::cout << "Integration rule not categorized" << std::endl;
          std::exit(1);
      }

      return open;
    }

    /// get total number of GPs of the 3D Gauss rule at compile time
    template <GaussRule3D rule>
    struct GaussRule3DToNumGaussPoints
    {
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule3D_undefined>
    {
      static const unsigned value = 0;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_8point>
    {
      static const unsigned value = 8;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_18point>
    {
      static const unsigned value = 18;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_27point>
    {
      static const unsigned value = 27;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_64point>
    {
      static const unsigned value = 64;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_125point>
    {
      static const unsigned value = 125;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_216point>
    {
      static const unsigned value = 216;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_343point>
    {
      static const unsigned value = 343;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_512point>
    {
      static const unsigned value = 512;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_729point>
    {
      static const unsigned value = 729;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_hex_1000point>
    {
      static const unsigned value = 1000;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_4point>
    {
      static const unsigned value = 4;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_4point_gauss_radau>
    {
      static const unsigned value = 4;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_5point>
    {
      static const unsigned value = 5;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_10point>
    {
      static const unsigned value = 10;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_11point>
    {
      static const unsigned value = 11;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_15point>
    {
      static const unsigned value = 15;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_24point>
    {
      static const unsigned value = 24;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_64point_peano>
    {
      static const unsigned value = 64;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_45point>
    {
      static const unsigned value = 45;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_125point_peano>
    {
      static const unsigned value = 125;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_343point_peano>
    {
      static const unsigned value = 343;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_tet_729point_peano>
    {
      static const unsigned value = 729;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_wedge_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_wedge_6point>
    {
      static const unsigned value = 6;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_wedge_9point>
    {
      static const unsigned value = 9;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_pyramid_1point>
    {
      static const unsigned value = 1;
    };
    template <>
    struct GaussRule3DToNumGaussPoints<intrule_pyramid_8point>
    {
      static const unsigned value = 8;
    };

    //! supported 2d integration rules
    //!
    //! integration rules can be distinguished into open and closed rules
    //! hereby, open rules are rules where all integration points are inside the domain
    //! and none on the boundary of the integrated domain
    enum GaussRule2D
    {
      intrule2D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value,
                             ///< yet!
      intrule_quad_1point,   ///< degree of precision: 1  (open)
      intrule_quad_4point,   ///< degree of precision: 3  (open)
      intrule_quad_6point,   ///< degree of precision: 3/5  (open) (2 by 3 gauss points for quad6
                             ///< element
      intrule_quad_9point,   ///< degree of precision: 5  (open)
      intrule_quad_16point,  ///< degree of precision: 7  (open)
      intrule_quad_25point,  ///< degree of precision: 9  (open)
      intrule_quad_lobatto25point,  ///< degree of precision: 7  (closed)
      intrule_quad_36point,         ///< degree of precision: 11 (open)
      intrule_quad_49point,         ///< degree of precision: 13 (open)
      intrule_quad_64point,         ///< degree of precision: 15 (open)
      intrule_quad_81point,         ///< degree of precision: 17 (open)
      intrule_quad_100point,        ///< degree of precision: 19 (open)
      intrule_quad_256point,        ///< degree of precision: 31 (open)
      intrule_quad_400point,        ///< degree of precision: 39 (open)
      intrule_quad_1024point,       ///< degree of precision: 63 (open)
      intrule_tri_1point,           ///< degree of precision: 1  (open)
      intrule_tri_3point,  ///< degree of precision: 2  (open)   (Hughes, The Finite Element Method)
      intrule_tri_3point_gauss_radau,  ///< degree of precision: 2  (closed) (Hughes, The Finite
                                       ///< Element Method)
      intrule_tri_4point,  ///< degree of precision: 3  (open)   (Hughes, The Finite Element Method)
      intrule_tri_6point,  ///< degree of precision: 4  (open)   (Hughes, The Finite Element Method)
      intrule_tri_7point,  ///< degree of precision: 5  (open)   (Gilbert Strang, George Fix, An
                           ///< Analysis of the Finite Element Method, Cambridge, 1973,ISBN:
                           ///< 096140888X,LC: TA335.S77.)
      intrule_tri_12point,  ///< degree of precision: 6  (open)   (Hughes, The Finite Element
                            ///< Method)
      intrule_tri_16point,  ///< degree of precision: 8  (open)   (Linbo Zhang, Tao Cui and Hui Liu:
                            ///< A SET OF SYMMETRIC QUADRATURE RULES ON TRIANGLES
                            ///<                                   AND TETRAHEDRA, Journal of
                            ///<                                   Computational Mathematics,
                            ///<                                   Vol.27, No.1, 2009, 89–96.)
      intrule_tri_37point,  ///< degree of precision: 13 (open)   (rule from ACM TOMS algorithm
                            ///< #706)
      intrule_tri_64point   ///< degree of precision: 15 (open)   (essentially a product of two 8
                            ///< point 1D Gauss-Legendre rules)
    };

    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule2D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case intrule_quad_1point:
          open = true;
          break;
        case intrule_quad_4point:
          open = true;
          break;
        case intrule_quad_6point:
          open = true;
          break;
        case intrule_quad_9point:
          open = true;
          break;
        case intrule_quad_16point:
          open = true;
          break;
        case intrule_quad_25point:
          open = true;
          break;
        case intrule_quad_lobatto25point:
          open = false;
          break;
        case intrule_quad_36point:
          open = true;
          break;
        case intrule_quad_49point:
          open = true;
          break;
        case intrule_quad_64point:
          open = true;
          break;
        case intrule_quad_81point:
          open = true;
          break;
        case intrule_quad_100point:
          open = true;
          break;
        case intrule_quad_256point:
          open = true;
          break;
        case intrule_quad_400point:
          open = true;
          break;
        case intrule_quad_1024point:
          open = true;
          break;
        case intrule_tri_1point:
          open = true;
          break;
        case intrule_tri_3point:
          open = true;
          break;
        case intrule_tri_3point_gauss_radau:
          open = false;
          break;
        case intrule_tri_6point:
          open = true;
          break;
        case intrule_tri_7point:
          open = true;
          break;
        case intrule_tri_12point:
          open = true;
          break;
        case intrule_tri_16point:
          open = true;
          break;
        case intrule_tri_37point:
          open = true;
          break;
        case intrule_tri_64point:
          open = true;
          break;
        default:
          std::cout << "Integration rule not categorized" << std::endl;
          std::exit(1);
      }

      return open;
    }

    //! supported 1d integration rules

    /*note: the following enum should be renamed the following way: GaussRule1D -> IntRule1D,
     *  intrule1D_line*point -> gauss1D_*point, intrule_line_lobatto*point -> lobatto1D_*point*/
    enum GaussRule1D
    {
      intrule1D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value,
                             ///< yet!
      intrule_line_1point,   ///< degree of precision: 1  (open)
      intrule_line_2point,   ///< degree of precision: 3  (open)
      intrule_line_3point,   ///< degree of precision: 5  (open)
      intrule_line_4point,   ///< degree of precision: 7  (open)
      intrule_line_5point,   ///< degree of precision: 9  (open)
      intrule_line_6point,   ///< degree of precision: 11 (open)
      intrule_line_7point,   ///< degree of precision: 13 (open)
      intrule_line_8point,   ///< degree of precision: 15 (open)
      intrule_line_9point,   ///< degree of precision: 17 (open)
      intrule_line_10point,  ///< degree of precision: 19 (open)
      intrule_line_16point,  ///< degree of precision: 31 (open)
      intrule_line_20point,  ///< degree of precision: 39 (open)
      intrule_line_32point,  ///< degree of precision: 63 (open)
      intrule_line_50point,  ///< degree of precision: 99 (open)
      intrule_line_lobatto2point,   ///< degree of precision: 1  (closed)
      intrule_line_lobatto3point,   ///< degree of precision: 3  (closed)
      intrule_line_lobatto4point,   ///< degree of precision: 5  (closed)
      intrule_line_lobatto5point,   ///< degree of precision: 7  (closed)
      intrule_line_lobatto6point,   ///< degree of precision: 9  (closed)
      intrule_line_lobatto7point,   ///< degree of precision: 11 (closed)
      intrule_line_lobatto8point,   ///< degree of precision: 13 (closed)
      intrule_line_lobatto9point,   ///< degree of precision: 15 (closed)
      intrule_line_lobatto10point,  ///< degree of precision: 17 (closed)
      intrule_line_lobatto11point,  ///< degree of precision: 19 (closed)
      intrule_line_lobatto12point,  ///< degree of precision: 21 (closed)
      intrule_line_lobatto13point,  ///< degree of precision: 23 (closed)
      intrule_line_lobatto14point   ///< degree of precision: 25 (closed)

    };

    /// allow to verify the openness of an integration rule
    ///
    /// \note: Open means that no integration point is on the boundary of the integrated domain
    inline bool IntRuleIsOpen(DRT::UTILS::GaussRule1D intrule)
    {
      bool open = false;

      switch (intrule)
      {
        case intrule_line_1point:
          open = true;
          break;
        case intrule_line_2point:
          open = true;
          break;
        case intrule_line_3point:
          open = true;
          break;
        case intrule_line_4point:
          open = true;
          break;
        case intrule_line_5point:
          open = true;
          break;
        case intrule_line_6point:
          open = true;
          break;
        case intrule_line_7point:
          open = true;
          break;
        case intrule_line_8point:
          open = true;
          break;
        case intrule_line_9point:
          open = true;
          break;
        case intrule_line_10point:
          open = true;
          break;
        case intrule_line_16point:
          open = true;
          break;
        case intrule_line_20point:
          open = true;
          break;
        case intrule_line_32point:
          open = true;
          break;
        case intrule_line_50point:
          open = true;
          break;
        case intrule_line_lobatto2point:
          open = false;
          break;
        case intrule_line_lobatto3point:
          open = false;
          break;
        case intrule_line_lobatto4point:
          open = false;
          break;
        case intrule_line_lobatto5point:
          open = false;
          break;
        case intrule_line_lobatto6point:
          open = false;
          break;
        case intrule_line_lobatto7point:
          open = false;
          break;
        case intrule_line_lobatto8point:
          open = false;
          break;
        case intrule_line_lobatto9point:
          open = false;
          break;
        case intrule_line_lobatto10point:
          open = false;
          break;
        case intrule_line_lobatto11point:
          open = false;
          break;
        case intrule_line_lobatto12point:
          open = false;
          break;
        case intrule_line_lobatto13point:
          open = false;
          break;
        case intrule_line_lobatto14point:
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

      /// get used integration rule
      DRT::UTILS::GaussRule1D GetIntRule() const { return intrule_; };
    };


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

#endif  // #ifndef UTILS_INTEGRATION_H
