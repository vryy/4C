/*!
\file drt_utils_integration.H

\brief structs to store integration points and weights together

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

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

#include <string>
#include <sstream>

namespace DRT
{
namespace UTILS
{

//! supported 3d integration rules
//!
//! integration rules can be distinguished into open and closed rules
//! hereby, open rules are rules where all integration points are inside the domain
//! and none on the boundary of the integrated domain
enum GaussRule3D
{
    intrule3D_undefined,            ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_hex_1point,             ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
    intrule_hex_8point,             ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
    intrule_hex_27point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 (open)
    intrule_tet_1point,             ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
    intrule_tet_4point,             ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (open)
    intrule_tet_4point_gauss_radau, ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (closed)
    intrule_tet_5point,             ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
    intrule_tet_10point,            ///< KEAST3           , DEG.OF PRECISION 3 (closed)
    intrule_tet_11point,            ///< KEAST5           , DEG.OF PRECISION 4 (open)
    intrule_tet_15point,            ///< KEAST6           , DEG.OF PRECISION 5 (closed)
    intrule_tet_24point,            ///< KEAST7           , DEG.OF PRECISION 6 (open)
    intrule_tet_64point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 6 (open)   (Peano, 1982)
    intrule_tet_45point,            ///< KEAST9           , DEG.OF PRECISION 8 (open)
    intrule_tet_125point_peano,     ///< GAUSS-LOBATTO    , DEG.OF PRECISION 9 (open)   (Peano, 1982)
    intrule_tet_343point_peano,     ///< GAUSS-LOBATTO    , DEG.OF PRECISION 13 (open)  (Peano, 1982)
    intrule_tet_729point_peano,     ///< GAUSS-LOBATTO    , DEG.OF PRECISION 16 (open)  (Peano, 1982)
    intrule_wedge_1point,           ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    intrule_wedge_6point,           ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    intrule_wedge_9point,           ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    intrule_pyramid_1point,         ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    intrule_pyramid_8point          ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
};


//! supported 2d integration rules
//!
//! integration rules can be distinguished into open and closed rules
//! hereby, open rules are rules where all integration points are inside the domain
//! and none on the boundary of the integrated domain
enum GaussRule2D
{
    intrule2D_undefined,             ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_quad_1point,             ///< degree of precision: 1  (open)
    intrule_quad_4point,             ///< degree of precision: 3  (open)
    intrule_quad_9point,             ///< degree of precision: 5  (open)
    intrule_tri_1point,              ///< degree of precision: 1  (open)
    intrule_tri_3point,              ///< degree of precision: 2  (open)   (Hughes, The Finite Element Method)
    intrule_tri_3point_gauss_radau,  ///< degree of precision: 2  (closed) (Hughes, The Finite Element Method)
    intrule_tri_6point,              ///< degree of precision: 4  (open)   (Hughes, The Finite Element Method)
    intrule_tri_7point,              ///< degree of precision: 5  (open)   (Gilbert Strang, George Fix, An Analysis of the Finite Element Method, Cambridge, 1973,ISBN: 096140888X,LC: TA335.S77.)
    intrule_tri_12point,             ///< degree of precision: 6  (open)   (Hughes, The Finite Element Method)
    intrule_tri_37point,             ///< degree of precision: 13 (open)   (rule from ACM TOMS algorithm #706)
    intrule_tri_64point              ///< degree of precision: 15 (open)   (essentially a product of two 8 point 1D Gauss-Legendre rules)
};

//! supported 1d integration rules
enum GaussRule1D
{
    intrule1D_undefined,           ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_line_1point,           ///< degree of precision: 1  (open)
    intrule_line_2point,           ///< degree of precision: 3  (open)
    intrule_line_3point,           ///< degree of precision: 5  (open)
    intrule_line_4point,           ///< degree of precision: 7  (open)
    intrule_line_5point,           ///< degree of precision: 9  (open)
    intrule_line_6point,           ///< degree of precision: 11 (open)
    intrule_line_7point,           ///< degree of precision: 13 (open)
    intrule_line_8point,           ///< degree of precision: 15 (open)
    intrule_line_9point,           ///< degree of precision: 17 (open)
    intrule_line_10point           ///< degree of precision: 19 (open)

};

/*!
 \brief a point with a position and an associated weight

 */
template <int dim, class Vec>
class IntegrationPoint
{
public:
    //! Standard constructor
    explicit IntegrationPoint(
            const Vec&       pos,
            const double     weight) :
                position_(pos),
                weight_(weight)
    {
        return;
    }

    //! copy constructor
    explicit IntegrationPoint(
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
    Vec Position() const { return position_; }

    //! get integration point weight
    double Weight() const { return weight_; }

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
  explicit IntegrationPoints3D(
      const GaussRule3D gaussrule   ///< 3d rule
      );
private:
  static const int max_nquad = 729;  ///< size of c array
public:
  int           nquad;               ///< number of gausspoints
  double        qxg[max_nquad][3];   ///< coordinates
  double        qwgt[max_nquad];     ///< weights
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
  explicit IntegrationPoints2D(
      const GaussRule2D gaussrule   ///< 2d rule
      );
private:
  static const int max_nquad = 64;    ///< size of c array
public:
  int           nquad;                ///< number of gausspoints
  double        qxg[max_nquad][2];    ///< coordinates
  double        qwgt[max_nquad];      ///< weights
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
  explicit IntegrationPoints1D(
      const GaussRule1D gaussrule   ///< 1d rule
      );
private:
  static const int max_nquad = 5;     ///< size of c array
public:
  int           nquad;                ///< number of integration points
  double        qxg[max_nquad][1];    ///< coordinates
  double        qwgt[max_nquad];      ///< weights
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
template<const int NSD>
class IntPointsAndWeights
{
public:
  /// Constructor
  template<class G>
  explicit IntPointsAndWeights(const G gaussrule);
};


//! specialization for 3D
template<>
class IntPointsAndWeights<3>
{
public:
  /// Constructor for 3D specialization
  explicit IntPointsAndWeights(const GaussRule3D gaussrule) :
  intpoints_(IntegrationPoints3D(gaussrule))
  {return; };

  //! get reference to the integration points and weights
  const IntegrationPoints3D& IP() const {return intpoints_;};

private:
  //! integration points and weights for 3D
  IntegrationPoints3D intpoints_;
};


//! specialization for 2D
template<>
class IntPointsAndWeights<2>
{
public:
  /// Constructor for 2D specialization
  explicit IntPointsAndWeights(const GaussRule2D gaussrule) :
  intpoints_(IntegrationPoints2D(gaussrule))
  {return; };

  //! get reference to the integration points and weights
  const IntegrationPoints2D& IP() const {return intpoints_;};

private:
  //! integration points and weights for 3D
  IntegrationPoints2D intpoints_;
};


//! specialization for 1D
template<>
class IntPointsAndWeights<1>
{
public:
  /// Constructor for 1D specialization
  explicit IntPointsAndWeights(const GaussRule1D gaussrule) :
  intpoints_(IntegrationPoints1D(gaussrule))
  {return; };

  //! get reference to the integration points and weights
  const IntegrationPoints1D& IP() const {return intpoints_;};

private:
  //! integration points and weights for 3D
  IntegrationPoints1D intpoints_;
};


} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_INTEGRATION_H
#endif  // #ifdef CCADISCRET
