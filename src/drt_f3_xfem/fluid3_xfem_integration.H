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


//! GaussRules supported by fluid3
enum GaussRule  
  {
    hex_1point,
    hex_8point,
    hex_27point,
    tet_1point,
    tet_4point,
    tet_4point_alternative,
    tet_10point
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


// evaluate element integration points (and storing them in data)
void integration_points_3d(struct _INTEGRATION_POINTS_3D& data,
                           const  GaussRule gaussrule);




#endif // #ifndef FLUID3_SHAPE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
