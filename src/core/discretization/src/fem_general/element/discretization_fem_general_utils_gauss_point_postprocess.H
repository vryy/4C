/*----------------------------------------------------------------------*/
/*! \file

\brief Postprocessing utilities for Gauss point quantities

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef DISCRETIZATION_FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_H
#define DISCRETIZATION_FEM_GENERAL_UTILS_GAUSS_POINT_POSTPROCESS_H

#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include "lib_element.H"

namespace DRT
{
  class Element;
}

namespace CORE::DRT::ELEMENTS
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
      ::DRT::Element& ele, const Epetra_SerialDenseMatrix& data, Epetra_MultiVector& nodal_data);

  /*!
   * @brief Averaging of all Gauss point quantities in @data within the element @ele and assembly to
   * the element vector @element_data
   *
   * @param ele (in) : Reference to finite element
   * @param data (in) : Gauss point data in a Matrix (numgp x numdim of vector)
   * @param element_data (out) : Assembled data
   */
  void EvaluateGaussPointQuantityAtElementCenter(
      ::DRT::Element& ele, const Epetra_SerialDenseMatrix& data, Epetra_MultiVector& element_data);


  /*!
   * \brief Assemble averaged data. The data at the Gauss points are averaged within the element.
   *
   * \tparam T Type of the data, either SerialDenseMatrix or CORE::LINALG::Matrix
   * \param global_data Global cell data
   * \param gp_data (numgp x size) matrix of the Gauss point data
   * \param ele element
   */
  template <class T>
  void AssembleAveragedElementValues(
      Epetra_MultiVector& global_data, const T& gp_data, const ::DRT::Element& ele);


  // --- template and inline functions --- //

  template <class T>
  void DRT::ELEMENTS::AssembleAveragedElementValues(
      Epetra_MultiVector& global_data, const T& gp_data, const ::DRT::Element& ele)
  {
    const Epetra_BlockMap& elemap = global_data.Map();
    int lid = elemap.LID(ele.Id());
    if (lid != -1)
    {
      for (decltype(gp_data.N()) i = 0; i < gp_data.N(); ++i)
      {
        double& s = (*(global_data(i)))[lid];  // resolve pointer for faster access
        s = 0.;
        for (decltype(gp_data.M()) j = 0; j < gp_data.M(); ++j)
        {
          s += gp_data(j, i);
        }
        s *= 1.0 / gp_data.M();
      }
    }
  }

}  // namespace CORE::DRT::ELEMENTS

#endif