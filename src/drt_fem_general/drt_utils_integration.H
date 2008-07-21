/*!
\file drt_utils_integration.H

\brief structs to store integration points and weights together

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef UTILS_INTEGRATION_H
#define UTILS_INTEGRATION_H

#include "../drt_lib/drt_element.H"

namespace DRT
{
namespace UTILS
{

//! supported 3d integration rules
enum GaussRule3D
{
    intrule3D_undefined,    ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_hex_1point,
    intrule_hex_8point,
    intrule_hex_27point,
    intrule_tet_1point,     ///< GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 0
    intrule_tet_4point,     ///< GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 1
    intrule_tet_4point_alternative,  ///< ALT.GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION ? (closed)
    intrule_tet_5point,     ///< GAUSS INTEGRATION        5 SAMPLING POINTS, DEG.OF PRECISION 2 (open)
    intrule_tet_10point,    ///< GAUSS INTEGRATION        10 SAMPLING POINTS, DEG.OF PRECISION 3 (closed)
    intrule_tet_11point,    ///< KEAST5,        11 SAMPLING POINTS, DEG.OF PRECISION 4 (open)
    intrule_tet_24point,    ///< KEAST7,        24 SAMPLING POINTS, DEG.OF PRECISION 6 (open)
    intrule_wedge_1point,
    intrule_wedge_6point,   ///< GAUSS INTEGRATION         6 SAMPLING POINTS, degree of precission ???
    intrule_wedge_9point,   ///< GAUSS INTEGRATION         9 SAMPLING POINTS, degree of precission ???
    intrule_pyramid_1point,
    intrule_pyramid_8point
};


//! supported 2d integration rules
enum GaussRule2D
{
    intrule2D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_quad_1point,
    intrule_quad_4point,
    intrule_quad_9point,
    intrule_tri_1point,
    intrule_tri_3point,    ///< gp on triangular lines/edges
    intrule_tri_3point_on_corners,
    intrule_tri_6point,
    intrule_tri_12point,   ///< degree of precision: 6 (Hughes, The Finite Element Method)
    intrule_tri_37point    ///< degree of precision: 13 (rule from ACM TOMS algorithm #706)
};

//! supported 1d integration rules
enum GaussRule1D
{
    intrule1D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_line_1point,
    intrule_line_2point,
    intrule_line_3point,
    intrule_line_4point,
    intrule_line_5point
};

/*!
 \brief a point with a position an an associated weight

 */
template <int dim, class Vec>
class IntegrationPoint
{
public:
    //! Default constructor
    explicit IntegrationPoint(
            const Vec&       pos,
            const double     weight) : 
                position_(pos),
                weight_(weight)
    {
        return;
    }
    
    //! copy constructor
    IntegrationPoint(
            const IntegrationPoint&     other) :
                position_(other.position_),
                weight_(other.weight_)
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
    inline Vec Position() const { return position_; }
    //! get integration point weight
    inline double Weight() const { return weight_; }
    
private:
    //! 1d, 2d or 3d position of integration point in element coordinates
    const Vec position_;
    //! integration weight
    const double weight_;
    
    //! Empty constructor
    explicit IntegrationPoint() :
                position_(),
                weight_(0.0)
    {
        return;
    }
}; // class IntegrationPoint


/*!
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.

\note tetrahedral integration rule weights are stored such that after integration over the parameter space,
we get the real volume of 0.1666..., NOT the unit volume of 1.0. Integration over the hex rules gives an volume of 8.0.
*/
struct IntegrationPoints3D
{
  /// constructor
  IntegrationPoints3D(
      const GaussRule3D gaussrule   ///< 3d rule
      );
  int           nquad;         ///< number of gausspoints
  double        qxg[27][3];    ///< coordinates
  double        qwgt[27];      ///< weights
};

/*!
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.

\note triangular integration rule weights are stored such that after integration over the parameter space,
we get the real area of 0.5, NOT the unit area of 1.0. Integration over the quad rules gives an area of 4.0.
*/
struct IntegrationPoints2D
{
  /// constructor
  IntegrationPoints2D(
      const GaussRule2D gaussrule   ///< 2d rule
      );
  int           nquad;        ///< number of gausspoints
  double        qxg[37][2];    ///< coordinates
  double        qwgt[37];      ///< weights
};

/*!
\brief integration parameters

<pre>                                                gerstenberger 06/07

In this structure the coordinates and weights used by numerical integration
are stored.

</pre>
*/
struct IntegrationPoints1D
{
  /// constructor
  IntegrationPoints1D(
      const GaussRule1D gaussrule   ///< 1d rule
      );
  int           nquad;        ///< number of integration points
  double        qxg[5];       ///< coordinates */
  double        qwgt[5];      ///< weights */
};

/*!
\brief return object with 3D integration points and weights
*/
inline IntegrationPoints3D getIntegrationPoints3D(
  const  GaussRule3D gaussrule   ///< used to switch between the implemented integration rules
  )
{
  return IntegrationPoints3D(gaussrule);
}

/*!
\brief return object with 2D integration points and weights
*/
inline IntegrationPoints2D getIntegrationPoints2D(
  const  GaussRule2D gaussrule   ///< used to switch between the implemented integration rules
  )
{
  return IntegrationPoints2D(gaussrule);
}

/*!
\brief return object with 1D integration points and weights
*/
inline IntegrationPoints1D getIntegrationPoints1D(
  const  GaussRule1D gaussrule   ///< used to switch between the implemented integration rules
  )
{
  return IntegrationPoints1D(gaussrule);
}

/*!
 *  returns the element size in local coordinates
 */    
double getSizeInLocalCoordinates( 
    const DRT::Element::DiscretizationType     distype);


} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_INTEGRATION_H
#endif  // #ifdef CCADISCRET
