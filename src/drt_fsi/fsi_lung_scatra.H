#ifdef CCADISCRET

#ifndef FSI_LUNG_SCATRA_H
#define FSI_LUNG_SCATRA_H

#include "../drt_scatra/passive_scatra_algorithm.H"

namespace FSI
{
  class LungScatra
  {

  public:

    LungScatra(Teuchos::RCP<FSI::Monolithic> fsi);

    ~LungScatra(){};

    Teuchos::RCP<Epetra_Vector> Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv);

    Teuchos::RCP<Epetra_Vector> Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv);

    /// read restart
    void ReadRestart();

    void SetupFSISystem();

    void Timeloop();

    void DoFsiStep();

    void DoScatraStep();

    void PrepareTimeStep();

    void EvaluateScatraFields();

    bool AbortScatraNonlinIter(const int itnum);

    void SetupCoupledScatraSystem();

    void SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector> globalvec,
                                  const Teuchos::RCP<const Epetra_Vector> vec1,
                                  const Teuchos::RCP<const Epetra_Vector> vec2);

    void ExtractScatraFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
                                   Teuchos::RCP<const Epetra_Vector>& vec1,
                                   Teuchos::RCP<const Epetra_Vector>& vec2);

    void SetupCoupledScatraMatrix();

    void LinearSolveScatra();

    void ScatraMergeSolve();

    void SetVelocityFields();

    /// extract fluid convective and structure grid velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel);

    /// update concentrations in each iteration
    void FieldUpdateIter();

    void UpdateAndOutput();

    /// make sure potential Dirichlet conditions at the scatra coupling
    /// interface are defined on both discretizations
    void CheckInterfaceDirichletBC();

  private:

    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra algorithms
    std::vector<Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> > scatravec_;

    /// vector of scatra field map extractors (coupled vs. uncoupled dofs)
    std::vector<LINALG::MultiMapExtractor> scatrafieldexvec_;

    /// coupling of dofs at the scatra interface
    ADAPTER::Coupling scatracoup_;

    /// dof row map of scatra problems splitted in (field) blocks
    LINALG::MultiMapExtractor scatraglobalex_;

    /// scatra rhs vector
    Teuchos::RCP<Epetra_Vector> scatrarhs_;

    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> scatraincrement_;

    /// scatra block matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrasystemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks
    UTILS::MatrixRowColTransform sbbtransform_;
    UTILS::MatrixRowTransform sbitransform_;
    UTILS::MatrixColTransform sibtransform_;
    ///@}
  };

}

#endif
#endif
