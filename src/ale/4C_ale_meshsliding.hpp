/*--------------------------------------------------------------------------*/
/*! \file

\brief Mesh sliding for ale problems


\level 3
*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_ALE_MESHSLIDING_HPP
#define FOUR_C_ALE_MESHSLIDING_HPP

#include "4C_config.hpp"

#include "4C_ale_meshtying.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class CouplingNonLinMortar;
}

namespace ALE
{
  class Meshsliding : public Meshtying
  {
   public:
    //! Constructor
    Meshsliding(Teuchos::RCP<DRT::Discretization> dis,          ///> actual discretisation
        CORE::LINALG::Solver& solver,                           ///> solver
        int msht,                                               ///> meshting parameter list
        int nsd,                                                ///> number space dimensions
        const UTILS::MapExtractor* surfacesplitter = nullptr);  ///> surface splitter

    //! Set up mesh sliding framework
    Teuchos::RCP<CORE::LINALG::SparseOperator> Setup(
        std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp) override;

   private:
    //! Call the constructor and the setup of the mortar coupling adapter
    void AdapterMortar(std::vector<int> coupleddof) override;

    //! Compare the size of the slave and master dof row map
    void CompareNumDof() override;

    //! Get function for the slave and master dof row map
    void DofRowMaps() override;

    //! Get function for the P matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMortarMatrixP() override;

    //! Condensation operation for a block matrix
    void condensation_operation_block_matrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator>&
            sysmat,                             ///> sysmat established by the element routine
        Teuchos::RCP<Epetra_Vector>& residual,  ///> residual established by the element routine
        Teuchos::RCP<Epetra_Vector>& dispnp) override;  ///> current displacement vector

    //! Get functions for the mortar matrices
    void GetMortarMatrices(Teuchos::RCP<CORE::LINALG::SparseMatrix>& Aco_mm,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& Aco_ms,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& Aco_sm,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& Aco_ss,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& N_m,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& N_s);

    //! Split the mortar matrix into its slave and its master part
    void SplitMortarMatrix(Teuchos::RCP<CORE::LINALG::SparseMatrix>& MortarMatrix,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& MasterMatrix,
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& SlaveMatrix,
        Teuchos::RCP<const Epetra_Map>& dofrowmap);

    //! Compute and update the increments of the slave node (do nothing in the mesh sliding case)
    void UpdateSlaveDOF(
        Teuchos::RCP<Epetra_Vector>& inc, Teuchos::RCP<Epetra_Vector>& dispnp) override{};

    //! Recover method for Lagrange multipliers
    void Recover(Teuchos::RCP<Epetra_Vector>& inc) override;

    //! Solve ALE mesh sliding problem
    int SolveMeshtying(CORE::LINALG::Solver& solver,
        Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector>& disi,
        Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<Epetra_Vector>& dispnp) override;

    //! adapter to nonlinear mortar coupling framework
    Teuchos::RCP<ADAPTER::CouplingNonLinMortar> adaptermeshsliding_;

    Teuchos::RCP<Epetra_Vector> lm_;  // current vector of Lagrange multipliers at t_n+1

    Teuchos::RCP<CORE::LINALG::SparseMatrix> a_ss_;   // stiffness block A_ss (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> a_sm_;   // stiffness block A_sm (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> a_sn_;   // stiffness block A_sn (needed for LM)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> d_inv_;  // inverse of Mortar matrix D (needed for LM)
    Teuchos::RCP<Epetra_Vector> rs_;                  // slave side effective forces (needed for LM)
  };

}  // namespace ALE

FOUR_C_NAMESPACE_CLOSE

#endif
