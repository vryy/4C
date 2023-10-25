/*----------------------------------------------------------------------*/
/*! \file
\brief Collection of free functions to reduce code duplication between elements

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SO3_ELEMENT_SERVICE_H
#define BACI_SO3_ELEMENT_SERVICE_H

#include "baci_so3_base.H"

#include <Epetra_IntVector.h>
#include <Epetra_MultiVector.h>

namespace DRT
{
  namespace UTILS
  {
    class GaussIntegration;
  }
  namespace ELEMENTS
  {
    /*!
     * \brief Assemble nodal element count
     *
     * \param global_count Add a 1 to all nodes belonging to this element
     * \param ele element
     */
    void AssembleNodalElementCount(Epetra_IntVector& global_count, const DRT::Element& ele);

    /*!
     * \brief Assemble Gauss point data into an array of global cell data
     *
     * \param global_data array of global cell data (length at least number of gauss points)
     * \param gp_data (numgp x size) matrix of the Gauss point data
     * \param ele element
     */
    void AssembleGaussPointValues(std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
        const CORE::LINALG::SerialDenseMatrix& gp_data, const DRT::Element& ele);

    /*!
     * @brief calculate and return the value of the quantity at position xi based on the
     * quantity node vector
     *
     * @tparam distype        discretization type of element
     * @param xi              position to project to in local coordinates
     * @param nodal_quantity  nodal vector of the quantity to be projected
     * @return quantities projected to position xi
     */
    template <CORE::FE::CellType distype>
    std::vector<double> ProjectNodalQuantityToXi(
        const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& nodal_quantity);
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SO3_ELEMENT_SERVICE_H