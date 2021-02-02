/*----------------------------------------------------------------------*/
/*! \file
\brief Collection of free functions to reduce code duplication between elements

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef SO_ELEMENT_SERVICE_H
#define SO_ELEMENT_SERVICE_H

#include <Epetra_MultiVector.h>
#include <Epetra_IntVector.h>
#include "so_base.H"

namespace DRT
{
  namespace ELEMENTS
  {
    /*!
     * \brief Assemble nodal data (that comes from an extrapolation from Gauss points) to the global
     * data
     *
     * If the last parameter is true, the nodal data of each node is divided by the number of
     * elements that share the specific node to get an averaged value at each node. If the parameter
     * is false, the division by the number of elements has to be done later.
     *
     * \tparam T Type of the data, either SerialDenseMatrix or LINALG::Matrix
     * \param global_data Global nodal data vector
     * \param nodal_data (numnod x size) matrix with data at the nodes
     * \param ele element
     * \param nodal_average whether or not to compute the nodal average based on the element that
     * share the nodes
     */
    template <class T>
    void AssembleExtrapolatedNodalValues(Epetra_MultiVector& global_data, const T& nodal_data,
        const DRT::ELEMENTS::So_base* ele, bool nodal_average = true);

    /*!
     * \brief Assemble averaged date. The data at the Gauss points are averaged within the element.
     *
     * \tparam T Type of the data, either SerialDenseMatrix or LINALG::Matrix
     * \param global_data Global cell data
     * \param gp_data (numgp x size) matrix of the Gauss point data
     * \param ele element
     */
    template <class T>
    void AssembleAveragedElementValues(
        Epetra_MultiVector& global_data, const T& gp_data, const DRT::ELEMENTS::So_base* ele);

    /*!
     * \brief Assemble nodal element count
     *
     * \param global_data Add a 1 to all nodes belonging to this element
     * \param ele element
     */
    void AssembleNodalElementCount(
        Epetra_IntVector& global_data, const DRT::ELEMENTS::So_base* ele);

    /*!
     * \brief Assemble Gauss point data into an array of global cell data
     *
     * \param global_data array of global cell data (length at least number of gauss points)
     * \param gp_data (numgp x size) matrix of the Gauss point data
     * \param ele element
     */
    void AssembleGaussPointValues(std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
        const Epetra_SerialDenseMatrix& gp_data, const DRT::ELEMENTS::So_base* ele);

  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SO_ELEMENT_SERVICE_H