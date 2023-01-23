/*----------------------------------------------------------------------*/
/*! \file

\brief Postprocessing utilities for Gauss point quantities

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_H
#define FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_H

#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseMatrix.h>

namespace DRT
{
  class Element;
}

namespace DRT::ELEMENTS
{

  /*!
   * @brief Extrapolation of Gauss point quantities given in @data to the nodes of the Element @ele
   * using the shape functions of the element and assembly to the global nodal data @nodal_data.
   *
   * @note On shared nodes, the values of all participating elements will be averaged
   *
   * @param ele (in) : Reference to finite element
   * @param data (in) : Gauss point data in a Matrix (numgp x numdim of vector)
   * @param nodal_data (out) : Assembled data
   */
  void ExtrapolateGaussPointQuantityToNodes(
      DRT::Element& ele, const Epetra_SerialDenseMatrix& data, Epetra_MultiVector& nodal_data);

  /*!
   * @brief Averaging of all Gauss point quantities in @data within the element @ele and assembly to
   * the element vector @element_data
   *
   * @param ele (in) : Reference to finite element
   * @param data (in) : Gauss point data in a Matrix (numgp x numdim of vector)
   * @param element_data (out) : Assembled data
   */
  void EvaluateGaussPointQuantityAtElementCenter(
      DRT::Element& ele, const Epetra_SerialDenseMatrix& data, Epetra_MultiVector& element_data);

}  // namespace DRT::ELEMENTS

#endif