/*!
 \file drt_utils_dgfem_basisfunctions.H

 \brief Provide function values of basis-functions used in STDG-FEM

 Provided are 1D, 2D and 3D shape functions with special functions to be used 
 for the FEM discretization in time.
 
 Use is made of Legendre polynomial basis-functions using the Rodrigues representation 
 for which the $p$-th order polynomial is given by $\phi_p=\frac{1}{2^p p!}\frac{d^p}{dx^p}(x^2-1)^p$
 and the first realizations are given by $1$,$x$,$(3x^2-1)/2$,$\ldots$.
 
 \author Fedderik van der Bos
 bos@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15253
 */

#ifdef CCADISCRET
#ifndef DRT_UTILS_DGFEM_BASISFUNCTIONS
#define DRT_UTILS_DGFEM_BASISFUNCTIONS

#include <blitz/array.h>

#include "drt_element.H"
#include "Epetra_Vector.h"

using namespace std;
using namespace Teuchos;

namespace DRT
{
namespace Utils
{
  /*!
  \brief Fill a blitz-vector with with 3D shape function
  */
  void DGFEM_BasisFunction_3D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with with 3D first shape function derivatives
  */
  void DGFEM_BasisFunction_3D_deriv1(
    blitz::Array<double, 2>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with with 3D second shape function derivatives
  */
  void DGFEM_BasisFunction_3D_deriv2(
    blitz::Array<double, 2>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-vector with 2D shape function
  */
  void DGFEM_BasisFunction_2D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with first 2D shape function derivative
  */
  void DGFEM_BasisFunction_2D_deriv1(
    blitz::Array<double, 2>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with second 2D shape function derivative
  */
  void DGFEM_BasisFunction_2D_deriv2(
    blitz::Array<double, 2>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-vector with 1D shape function
  */
  void DGFEM_BasisFunction_1D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with first 1D shape function derivative
  */
  void DGFEM_BasisFunction_1D_deriv1(
    blitz::Array<double, 1>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a blitz-matrix with second 1D shape function derivative
  */
  void DGFEM_BasisFunction_1D_deriv2(
    blitz::Array<double, 1>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  shape,   ///< distuinguish between shape of the element
    const int&                               npol     ///< Number of polynomials
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with with 3D shape function
  */
  void DGFEM_BasisFunction_3D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D first shape function derivatives
  */
  void DGFEM_BasisFunction_3D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D second shape function derivatives
  */
  void DGFEM_BasisFunction_3D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with 2D shape function
  */
  void DGFEM_BasisFunction_2D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 2D shape function derivative
  */
  void DGFEM_BasisFunction_2D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 2D shape function derivative
  */
  void DGFEM_BasisFunction_2D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with 1D shape function
  */
  void DGFEM_BasisFunction_1D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 1D shape function derivative
  */
  void DGFEM_BasisFunction_1D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 1D shape function derivative
  */
  void DGFEM_BasisFunction_1D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

} // namespace Utils
} // namespace DRT

#endif  // #ifdef DRT_UTILS_DGFEM_BASISFUNCTIONS
#endif  // #ifdef CCADISCRET


