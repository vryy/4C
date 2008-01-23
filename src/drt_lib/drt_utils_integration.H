/*!
\file drt_utils_integration.H

\brief structs to store gausspoints and weights together

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

#include "drt_element.H"

#include "Epetra_Vector.h"
using namespace std;

namespace DRT
{
namespace UTILS
{

//! supported 3d gauss rules
enum GaussRule3D
{
    intrule3D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_hex_1point,
    intrule_hex_8point,
    intrule_hex_27point,
    intrule_tet_1point,
    intrule_tet_4point,
    intrule_tet_4point_alternative,
    intrule_tet_5point,
    intrule_tet_10point,
    intrule_tet_11point,
    intrule_tet_24point,
    intrule_wedge_1point,
    intrule_wedge_6point,
    intrule_wedge_9point,
    intrule_pyramid_1point,
    intrule_pyramid_8point
};


//! supported 2d gauss rules
enum GaussRule2D
{
    intrule2D_undefined,   ///< use this to initialize a gaussrule, if you don't know the value, yet!
    intrule_quad_1point,
    intrule_quad_4point,
    intrule_quad_9point,
    intrule_tri_1point,
    intrule_tri_3point,
    intrule_tri_3point_on_corners,
    intrule_tri_6point
};

//! supported 1d gauss rules
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
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.
*/
struct IntegrationPoints3D
{
  IntegrationPoints3D(const GaussRule3D gaussrule);
  int           nquad;         ///< number of gausspoints
  double        qxg[84][3];    ///< coordinates
  double        qwgt[84];      ///< weights
};

/*!
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.
*/
struct IntegrationPoints2D
{
  IntegrationPoints2D(const GaussRule2D gaussrule);
  int           nquad;        ///< number of gausspoints
  double        qxg[9][2];    ///< coordinates
  double        qwgt[9];      ///< weights
};

/*!
\brief integration parameters

<pre>                                                gerstenberger 06/07

In this structure the coordinates and weights used by gauss integration
are stored.

</pre>
*/
struct IntegrationPoints1D
{
  IntegrationPoints1D(const GaussRule1D gaussrule);
  int           nquad;        ///< number of gausspoints
  double        qxg[5];       ///< coordinates */
  double        qwgt[5];      ///< weights */
};

/*!
\brief return object with 3D integration points and weights
*/
inline IntegrationPoints3D getIntegrationPoints3D(
  const  GaussRule3D gaussrule   ///< used to switch between the implemented gaussrules
  )
{
  return IntegrationPoints3D(gaussrule);
}

/*!
\brief return object with 2D integration points and weights
*/
inline IntegrationPoints2D getIntegrationPoints2D(
  const  GaussRule2D gaussrule   ///< used to switch between the implemented gaussrules
  )
{
  return IntegrationPoints2D(gaussrule);
}

/*!
\brief return object with 1D integration points and weights
*/
inline IntegrationPoints1D getIntegrationPoints1D(
  const  GaussRule1D gaussrule   ///< used to switch between the implemented gaussrules
  )
{
  return IntegrationPoints1D(gaussrule);
}

/*----------------------------------------------------------------------*
 |  returns the elementsize in local coordinates             a.ger 11/07|
 *----------------------------------------------------------------------*/    
double getSizeInLocalCoordinates( 
    const DRT::Element::DiscretizationType     distype);


} // namespace UTILS
} // namespace DRT

#endif // #ifndef UTILS_INTEGRATION_H
#endif  // #ifdef CCADISCRET
