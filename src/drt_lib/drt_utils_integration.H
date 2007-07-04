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
#ifdef TRILINOS_PACKAGE
#ifndef UTILS_INTEGRATION_H
#define UTILS_INTEGRATION_H


#include "Epetra_Vector.h"
using namespace std;

namespace DRT
{
namespace Utils
{

//! supported 3d gauss rules
enum GaussRule3D
{
    intrule_hex_1point,
    intrule_hex_8point,
    intrule_hex_27point,
    intrule_tet_1point,
    intrule_tet_4point,
    intrule_tet_4point_alternative,
    intrule_tet_10point,
    intrule_weg_9point
};


//! supported 2d gauss rules
enum GaussRule2D
{
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
    intrule_line_1point,
    intrule_line_2point,
    intrule_line_3point
};

/*!
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.
*/
typedef struct
{
  int           nquad;               ///< number of gausspoints
  double        qxg[MAXGAUSS][3];    ///< coordinates
  double        qwgt[MAXGAUSS];      ///< weights
} IntegrationPoints3D;

/*!
\brief integration parameters
\author Axel Gerstenberger

In this structure the coordinates and weights used by gauss integration are stored.
*/
typedef struct
{
  int           nquad;               ///< number of gausspoints
  double        qxg[MAXGAUSS][2];    ///< coordinates
  double        qwgt[MAXGAUSS];      ///< weights
} IntegrationPoints2D;

/*!
\brief integration parameters

<pre>                                                gerstenberger 06/07

In this structure the coordinates and weights used by gauss integration
are stored.

</pre>
*/
typedef struct
{
  int           nquad;               ///< number of gausspoints
  double        qxg[MAXGAUSS];       ///< coordinates */
  double        qwgt[MAXGAUSS];      ///< weights */
} IntegrationPoints1D;

/*!
\brief return object with 3D integration points and weights
*/
IntegrationPoints3D getIntegrationPoints3D(
  const  GaussRule3D gaussrule   ///< used to switch between the implemented gaussrules
  );

/*!
\brief return object with 2D integration points and weights
*/
IntegrationPoints2D getIntegrationPoints2D(
  const  GaussRule2D gaussrule   ///< used to switch between the implemented gaussrules
  );

/*!
\brief return object with 1D integration points and weights
*/
IntegrationPoints1D getIntegrationPoints1D(
  const  GaussRule1D gaussrule   ///< used to switch between the implemented gaussrules
  );

} // namespace Utils
} // namespace DRT

#endif // #ifndef UTILS_INTEGRATION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
