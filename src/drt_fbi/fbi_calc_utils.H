/*----------------------------------------------------------------------*/
/*! \file
 *
 *\brief Utility functions for fluid beam interaction related calculations
 *
 *\level 3
 *
 */
/*----------------------------------------------------------------------*/

#ifndef SRC_DRT_FBI_FBI_CALC_UTILS_H
#define SRC_DRT_FBI_FBI_CALC_UTILS_H

#include <vector>
#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG

namespace FBI
{
  namespace UTILS
  {
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Pointer to the fluid or structure discretization.
     * @param ele (in) Pointer to the fluid or beam element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     *
     */
    void GetFBIElementCenterlineDOFIndices(DRT::Discretization const& discret,
        const DRT::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
        unsigned int& num_dof);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assembles centerline DOF coupling contributions into element coupling matrices and
     * force vectors
     *
     * The function expects a vector object of length 2 containing discretizations
     *
     * \param[in] discretizations vector to the structure and then fluid discretization
     * \param[in] elegid global id for the beam and fluid element
     * \param[in] eleforce_centerlineDOFs force vectors containing contributions from the centerline
     * dofs in the case of beams, otherwise just the classical element DOFs
     * \param[in]
     * elestiffcenterlineDOFs stiffness matrices containing contributions from the centerline dofs
     * in the case of beams, otherwise just the classical element DOFs
     * \param[inout] eleforce
     * element force vector
     * \param[inout] elestiff element stiffness matrix
     *
     */

    void AssembleCenterlineDofForceStiffIntoFBIElementForceStiff(
        const DRT::Discretization& discretization1, const DRT::Discretization& discretization2,
        std::vector<int> const& elegid,
        std::vector<LINALG::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<LINALG::SerialDenseVector>* eleforce,
        std::vector<std::vector<LINALG::SerialDenseMatrix>>* elestiff);

    /**
     * \brief Assembles element coupling contributions into global coupling matrices and force
     * vectors needed for partitioned algorithms
     *
     * \param[in] discretization1 discretization to the first field
     * \param[in] discretization2 discretization to the second field
     * \param[in] elegids vector of length 2 containing the global IDs of the interacting elements
     * \param[in] elevec vector of length 2 containing the discrete element residual vectors of the
     * interacting elements
     * \param[in, out] c22 coupling matrix relating DOFs in the second discretization to each other
     * \param[in, out] c11 coupling matrix relating DOFs in the first discretization to each other
     * \param[in, out] c12 coupling matrix relating DOFs in the second discretization to DOFs in the
     * first discretization \param[in, out] c21 coupling matrix relating DOFs in the first
     * discretization to DOFs in the second discretization
     *
     */
    void FEAssembleEleForceStiffIntoSystemVectorMatrices(const DRT::Discretization& discretization1,
        const DRT::Discretization& discretization2, std::vector<int> const& elegid,
        std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector>& f1, Teuchos::RCP<Epetra_FEVector>& f2,
        Teuchos::RCP<LINALG::SparseMatrix>& c11, Teuchos::RCP<LINALG::SparseOperator> c22,
        Teuchos::RCP<LINALG::SparseMatrix>& c12, Teuchos::RCP<LINALG::SparseMatrix>& c21);

  }  // namespace UTILS
}  // namespace FBI



#endif /* SRC_DRT_FBI_FBI_CALC_UTILS_H */
