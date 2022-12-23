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
  namespace UTILS
  {
    class GaussIntegration;
  }
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
        const DRT::Element& ele, bool nodal_average = true);

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
        Epetra_MultiVector& global_data, const T& gp_data, const DRT::Element& ele);

    /*!
     * \brief Assemble nodal element count
     *
     * \param global_data Add a 1 to all nodes belonging to this element
     * \param ele element
     */
    void AssembleNodalElementCount(Epetra_IntVector& global_data, const DRT::Element& ele);

    /*!
     * \brief Assemble Gauss point data into an array of global cell data
     *
     * \param global_data array of global cell data (length at least number of gauss points)
     * \param gp_data (numgp x size) matrix of the Gauss point data
     * \param ele element
     */
    void AssembleGaussPointValues(std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
        const Epetra_SerialDenseMatrix& gp_data, const DRT::Element& ele);

    template <DRT::Element::DiscretizationType distype>
    void ExtrapolateGPQuantityToNodesAndAssemble(const DRT::Element& ele,
        const LINALG::SerialDenseMatrix& gp_data, Epetra_MultiVector& global_data,
        bool nodal_average, const DRT::UTILS::GaussIntegration&);

    template <DRT::Element::DiscretizationType distype>
    LINALG::SerialDenseMatrix EvaluateGaussPointsToNodesExtrapolationMatrix(
        const DRT::UTILS::GaussIntegration&);

    /*!
     * @brief calculate and return the value of the quantity at position xi based on the
     * quantity node vector
     *
     * @tparam distype        discretization type of element
     * @param xi              position to project to in local coordinates
     * @param nodal_quantity  nodal vector of the quantity to be projected
     * @return quantities projected to position xi
     */
    template <DRT::Element::DiscretizationType distype>
    std::vector<double> ProjectNodalQuantityToXi(
        const LINALG::Matrix<3, 1>& xi, const std::vector<double>& nodal_quantity);
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SO_ELEMENT_SERVICE_H