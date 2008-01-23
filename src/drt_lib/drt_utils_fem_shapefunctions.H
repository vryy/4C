/*!
 \file drt_utils_fem_shapefunctions.H

 \brief Provide a node numbering scheme together with a set of shape functions

 Provided are 1D, 2D and 3D shape functions

 The surface mapping gives the node numbers such that the 2D shape functions can be used
 Nodal mappings describe the relation between volume, surface and line node numbering.
 They should be used as the only reference for such relationships.
 The corresponding graphics and a detailed description can be found in the Baci guide in the Convention chapter.
 The numbering of lower order elements is included in the higher order element, such that
 e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping

 !!!!
 The corresponding graphics and a detailed description can be found
 in the Baci guide in the Convention chapter.
 !!!!

 \author Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15236

 */

#ifdef CCADISCRET
#ifndef DRT_UTILS_FEM_SHAPEFUNCTION
#define DRT_UTILS_FEM_SHAPEFUNCTION

#include <blitz/array.h>

#include "drt_element.H"
#include "Epetra_Vector.h"

using namespace std;
using namespace Teuchos;

namespace DRT
{
namespace UTILS
{

  /*!
  \brief Fill blitz-vector with with 3D shape function
  */
  void shape_function_3D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 3D first shape function derivatives
  */
  void shape_function_3D_deriv1(
    blitz::Array<double, 2>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 3D second shape function derivatives
  */
  void shape_function_3D_deriv2(
    blitz::Array<double, 2>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with 2D shape function
  */
  void shape_function_2D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with first 2D shape function derivative
  */
  void shape_function_2D_deriv1(
    blitz::Array<double, 2>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with second 2D shape function derivative
  */
  void shape_function_2D_deriv2(
    blitz::Array<double, 2>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with 1D shape function
  */
  void shape_function_1D(
    blitz::Array<double, 1>&                 funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with first 1D shape function derivative
  */
  void shape_function_1D_deriv1(
    blitz::Array<double, 2>&                 deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with second 1D shape function derivative
  */
  void shape_function_1D_deriv2(
    blitz::Array<double, 2>&                 deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with with 3D shape function
  */
  void shape_function_3D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D first shape function derivatives
  */
  void shape_function_3D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D second shape function derivatives
  */
  void shape_function_3D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const double&                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with 2D shape function
  */
  void shape_function_2D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 2D shape function derivative
  */
  void shape_function_2D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 2D shape function derivative
  */
  void shape_function_2D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const double&                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with 1D shape function
  */
  void shape_function_1D(
    Epetra_SerialDenseVector&                funct,   ///< to be filled with shape function values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 1D shape function derivative
  */
  void shape_function_1D_deriv1(
    Epetra_SerialDenseMatrix&                deriv1,  ///< to be filled with shape function derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 1D shape function derivative
  */
  void shape_function_1D_deriv2(
    Epetra_SerialDenseMatrix&                deriv2,  ///< to be filled with shape function 2-nd derivative values
    const double&                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType&  distype  ///< distuinguish between DiscretizationType
    );

} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_FEM_SHAPEFUNCTION
#endif  // #ifdef CCADISCRET


