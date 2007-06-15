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
#ifndef FLUID3_XFEM_SHAPE_H
#define FLUID3_XFEM_SHAPE_H


#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;



void shape_function_3D_deriv0( Epetra_SerialDenseVector&                  funct,
                               const double&                              r,
                               const double&                              s,
                               const double&                              t,
                               const DRT::Element::DiscretizationType&    distype);

void shape_function_3D_deriv1( Epetra_SerialDenseMatrix&                  deriv,
                               const double&                              r,
                               const double&                              s,
                               const double&                              t,
                               const DRT::Element::DiscretizationType&    distype);
                            
void shape_function_3D_deriv2( Epetra_SerialDenseMatrix&                  deriv2,
                               const double&                              r,
                               const double&                              s,
                               const double&                              t,
                               const DRT::Element::DiscretizationType&    distype);


#endif // #ifndef FLUID3_SHAPE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
