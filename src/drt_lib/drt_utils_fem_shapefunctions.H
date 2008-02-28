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
#include "drt_utils_local_connectivity_matrices.H"

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

  /*!
  \brief Fill a blitz-vector with with 3D shape function
  
  \return shape function values
  */
  blitz::Array<double, 1> shape_function_3D(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const double                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 3D first shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_3D_deriv1(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const double                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 3D second shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_3D_deriv2(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const double                            t,       ///< t coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );
  
  /*!
  \brief Fill a blitz-vector with with 2D shape function
  
  \return shape function values
  */
  blitz::Array<double, 1> shape_function_2D(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 2D first shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_2D_deriv1(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-matrix with with 2D second shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_2D_deriv2(
    const double                            r,       ///< r coordinate
    const double                            s,       ///< s coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );
  
  /*!
  \brief Fill a blitz-vector with with 1D shape function
  
  \return shape function values
  */
  blitz::Array<double, 1> shape_function_1D(
    const double                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-vector with with 1D first shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_1D_deriv1(
    const double                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );

  /*!
  \brief Fill a blitz-vector with with 1D second shape function derivatives
  
  \return shape function derivative values
  */
  blitz::Array<double, 2> shape_function_1D_deriv2(
    const double                            r,       ///< r coordinate
    const DRT::Element::DiscretizationType  distype  ///< distuinguish between DiscretizationType
    );
  
  //////////////////////////////////////////////////////////////////
  // versions that return the shape function and derivatives as return argument
  // allows to declare them const when they are used

  //
  // shape functions
  //
  template<DRT::Element::DiscretizationType distype, class T>
  blitz::Array<double, 1> shape_function(
          const T& xsi)
  {
    static blitz::Array<double, 1> f(DRT::UTILS::getNumberOfElementNodes<distype>());
    switch(DRT::UTILS::getDimension<distype>())
    {
        case 1:
        {
            shape_function_1D(f, xsi(0), distype);
            break;
        }
        case 2:
        {
            shape_function_2D(f, xsi(0), xsi(1), distype);
            break;
        }
        case 3:
        {
            shape_function_3D(f, xsi(0), xsi(1), xsi(2), distype);
            break;
        }
        default:
            dserror("dimension of the element is not correct");
    }
    return f;
  }

  //
  // first natural derivative of shape functions
  //
  template<DRT::Element::DiscretizationType distype, class T>
  blitz::Array<double, 2> shape_function_deriv1(
          const T& xsi)
  {
    switch(DRT::UTILS::getDimension<distype>())
    {
        case 1:
        {
            static blitz::Array<double, 2> d(1, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_1D_deriv1(d, xsi(0), distype);
            return d;
            break;
        }
        case 2:
        {
            static blitz::Array<double, 2> d(2, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_2D_deriv1(d, xsi(0), xsi(1), distype);
            return d;
            break;
        }
        case 3:
        {
            static blitz::Array<double, 2> d(3, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_3D_deriv1(d, xsi(0), xsi(1), xsi(2), distype);
            return d;
            break;
        }
        default:
            dserror("dimension of the element is not correct");
            exit(1);
    }
  }

  //
  // Second natural derivative of shape functions
  //
  template<DRT::Element::DiscretizationType distype, class T>
  blitz::Array<double, 2> shape_function_deriv2(
          const T& xsi)
  {
    switch(DRT::UTILS::getDimension<distype>())
    {
        case 1:
        {
            static blitz::Array<double, 2> d2(1, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_1D_deriv2(d2, xsi(0), distype);
            return d2;
            break;
        }
        case 2:
        {
            static blitz::Array<double, 2> d2(3, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_2D_deriv2(d2, xsi(0), xsi(1), distype);
            return d2;
            break;
        }
        case 3:
        {
            static blitz::Array<double, 2> d2(6, DRT::UTILS::getNumberOfElementNodes<distype>(), blitz::ColumnMajorArray<2>());
            shape_function_3D_deriv2(d2, xsi(0), xsi(1), xsi(2), distype);
            return d2;
            break;
        }
        default:
            dserror("dimension of the element is not correct");
            exit(1);
    }
  }

} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_FEM_SHAPEFUNCTION
#endif  // #ifdef CCADISCRET


