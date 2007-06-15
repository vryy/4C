/*!----------------------------------------------------------------------
\file fluid3_xfem.H

\class DRT::Elements::XFluid3

\brief A XFEM enabled copy of the fluid3 element

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef FLUID3_XFEM_INTEGRATION_H
#define FLUID3_XFEM_INTEGRATION_H


#include "Epetra_Vector.h"
using namespace std;


//! GaussRules
enum GaussRule3D  
  {
    hex_1point,
    hex_8point,
    hex_27point,
    tet_1point,
    tet_4point,
    tet_4point_alternative,
    tet_10point
};


//! GaussRules
enum GaussRule2D  
  {
    quad_1point,
    quad_4point,
    quad_9point,
    tri_1point,
    tri_3point,
    tri_6point
};

/*!----------------------------------------------------------------------
\brief integration parameters

<pre>                                                gerstenberger 06/07

In this structure the coordinates and weights used by gauss integration
are stored.

</pre>

------------------------------------------------------------------------*/
typedef struct _INTEGRATION_POINTS_2D
{
  int           nquad;
  double        qxg[MAXGAUSS][2];    /*!< coordinates */
  double        qwgt[MAXGAUSS];      /*!< weights */
} INTEGRATION_POINTS_2D;

/*!----------------------------------------------------------------------
\brief integration parameters

<pre>                                                gerstenberger 06/07

In this structure the coordinates and weights used by gauss integration
are stored.

</pre>

------------------------------------------------------------------------*/
typedef struct _INTEGRATION_POINTS_3D
{
  int           nquad;
  double        qxg[MAXGAUSS][3];    /*!< coordinates */
  double        qwgt[MAXGAUSS];      /*!< weights */
} INTEGRATION_POINTS_3D;


// evaluate element integration points (and storing them in intpoints)
void integration_points_3d(struct _INTEGRATION_POINTS_3D& intpoints,
                           const  GaussRule3D gaussrule);

// evaluate element integration points (and storing them in intpoints)
void integration_points_2d(struct _INTEGRATION_POINTS_2D& intpoints,
                           const  GaussRule2D gaussrule);


#endif // #ifndef FLUID3_SHAPE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
