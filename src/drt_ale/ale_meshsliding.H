/*--------------------------------------------------------------------------*/
/*!
\file ale_meshsliding.H

\brief Mesh sliding for ale problems

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_ALE_ALE_MESHSLIDING_H_
#define SRC_DRT_ALE_ALE_MESHSLIDING_H_

#include "ale_meshtying.H"

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
  Meshsliding(Teuchos::RCP<DRT::Discretization>      dis,             ///> actual discretisation
            LINALG::Solver&               solver,                   ///> solver
            int                           msht,                   ///> meshting parameter list
            int                           nsd,                      ///> number space dimensions
            const UTILS::MapExtractor*    surfacesplitter = NULL);  ///> surface splitter

  virtual ~Meshsliding(){return;};

  //! Set up mesh sliding framework
  Teuchos::RCP<LINALG::SparseOperator> Setup(std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp);

private:

  //! Call the constructor and the setup of the mortar coupling adapter
  void AdapterMortar(std::vector<int> coupleddof);

  //! Compare the size of the slave and master dof row map
  void CompareNumDof();

  //! Get function for the slave and master dof row map
  void DofRowMaps();

  //! Get function for the P matrix
  Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo();

  //! Condensation operation for a block matrix
  void CondensationOperationBlockMatrix(
      Teuchos::RCP<LINALG::SparseOperator>&  sysmat,         ///> sysmat established by the element routine
      Teuchos::RCP<Epetra_Vector>&           residual,       ///> residual established by the element routine
      Teuchos::RCP<Epetra_Vector>&           dispnp);        ///> current displacement vector

  //! Get functions for the mortar matrices
  void GetMortarMatrices(
      Teuchos::RCP<LINALG::SparseMatrix>& Aco_mm,
      Teuchos::RCP<LINALG::SparseMatrix>& Aco_ms,
      Teuchos::RCP<LINALG::SparseMatrix>& Aco_sm,
      Teuchos::RCP<LINALG::SparseMatrix>& Aco_ss,
      Teuchos::RCP<LINALG::SparseMatrix>& N_m,
      Teuchos::RCP<LINALG::SparseMatrix>& N_s);

  //! Split the mortar matrix into its slave and its master part
  void SplitMortarMatrix(
      Teuchos::RCP<LINALG::SparseMatrix>& MortarMatrix,
      Teuchos::RCP<LINALG::SparseMatrix>& MasterMatrix,
      Teuchos::RCP<LINALG::SparseMatrix>& SlaveMatrix,
      Teuchos::RCP<const Epetra_Map>&           dofrowmap);

  //! Compute and update the increments of the slave node (do nothing in the mesh sliding case)
  void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>&   inc,
                          Teuchos::RCP<Epetra_Vector>&   dispnp){};

  //! Recover method for Lagrange multipliers
  void Recover(Teuchos::RCP<Epetra_Vector>& inc);

  //! Solve ALE mesh sliding problem
  int SolveMeshtying(
    LINALG::Solver&                       solver,
    Teuchos::RCP<LINALG::SparseOperator>  sysmat,
    Teuchos::RCP<Epetra_Vector>&          disi,
    Teuchos::RCP<Epetra_Vector>           residual,
    Teuchos::RCP<Epetra_Vector>&          dispnp);

  //! adapter to nonlinear mortar coupling framework
  Teuchos::RCP<ADAPTER::CouplingNonLinMortar> adaptermeshsliding_;

  Teuchos::RCP<Epetra_Vector>  lm_;          // current vector of Lagrange multipliers at t_n+1

  Teuchos::RCP<LINALG::SparseMatrix> A_ss_;  // stiffness block A_ss (needed for LM)
  Teuchos::RCP<LINALG::SparseMatrix> A_sm_;  // stiffness block A_sm (needed for LM)
  Teuchos::RCP<LINALG::SparseMatrix> A_sn_;  // stiffness block A_sn (needed for LM)
  Teuchos::RCP<LINALG::SparseMatrix> D_inv_; // inverse of Mortar matrix D (needed for LM)
  Teuchos::RCP<Epetra_Vector>        rs_;    // slave side effective forces (needed for LM)

};

}

#endif /* SRC_DRT_ALE_ALE_MESHSLIDING_H_ */
