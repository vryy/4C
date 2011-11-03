#ifdef CCADISCRET

#ifndef FS3I_1WC_H
#define FS3I_1WC_H

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "fs3i.H"
#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_fsi/fsi_matrixtransform.H"

namespace FS3I
{
  class FS3I_1WC : public FS3I_Base
  {
  public:

    FS3I_1WC(Epetra_Comm& comm);

    virtual ~FS3I_1WC(){};

    Teuchos::RCP<Epetra_Vector> Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv);

    Teuchos::RCP<Epetra_Vector> Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv);

    virtual void ReadRestart();

    virtual void TestResults(const Epetra_Comm& comm);

    virtual void SetupSystem();

    virtual void Timeloop();

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

    void SetupCoupledScatraRHS();

    void SetupCoupledScatraMatrix();

    void LinearSolveScatra();

    void CoupledScatraSolve();

    void SetVelocityFields();

    void SetMeshDisp();

    /// extract fluid convective and structure grid velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel);

    /// update concentrations in each iteration
    void FieldUpdateIter();

    virtual void Output();

    void UpdateScatraFields();

    /// make sure potential Dirichlet conditions at the scatra coupling
    /// interface are defined on both discretizations
    void CheckInterfaceDirichletBC();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }


  protected:

    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra algorithms
    std::vector<Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> > scatravec_;


  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

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

    /// coupling forces (in case of surface permeability)
    std::vector<Teuchos::RCP<Epetra_Vector> > scatracoupforce_;

    /// coupling matrices (in case of surface permeability)
    std::vector<Teuchos::RCP<LINALG::SparseMatrix> > scatracoupmat_;

    /// zero vector (needed for application of Dirichlet BC on
    /// coupling vector)
    std::vector<Teuchos::RCP<Epetra_Vector> > scatrazeros_;

    /// scatra solver
    Teuchos::RCP<LINALG::Solver> scatrasolver_;

    /// consideration of surface permeability
    bool permeablesurf_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks
    ::FSI::UTILS::MatrixRowColTransform sbbtransform_;
    ::FSI::UTILS::MatrixRowTransform sbitransform_;
    ::FSI::UTILS::MatrixColTransform sibtransform_;

    ::FSI::UTILS::MatrixRowTransform fbitransform_;
    ///@}
  };

}

#endif
#endif
