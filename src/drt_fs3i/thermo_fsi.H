#ifdef CCADISCRET

#ifndef THERMOFSI_H
#define THERMOFSI_H

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "fs3i.H"

namespace FS3I
{
  class ThermoFSI : public FS3I_Base
  {
  public:

    ThermoFSI(const Epetra_Comm& comm);

    virtual ~ThermoFSI(){};

    virtual void SetupSystem();

    virtual void Timeloop();

    virtual void TestResults(const Epetra_Comm& comm) { dserror("not implemented"); return; };

    void DoScatraStep();

    void SetupCoupledScatraMatrix();

    void SetupCoupledScatraRHS();

    void ExtractScatraFieldVectors(Teuchos::RCP<const Epetra_Vector>  globalvec,
                                   Teuchos::RCP<const Epetra_Vector>& vec1,
                                   Teuchos::RCP<const Epetra_Vector>& vec2);

    void IterUpdate();

    void UpdateScatraFields();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }


  protected:

  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// vector of scatra field map extractors (coupled vs. uncoupled dofs)
    std::vector<LINALG::MultiMapExtractor> scatrafieldexvec_;

    /// dof row map of scatra problems splitted in (field) blocks
    LINALG::MultiMapExtractor scatraglobalex_;

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

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks
    FSI::UTILS::MatrixRowColTransform sbbtransform_;
    FSI::UTILS::MatrixRowTransform sbitransform_;
    FSI::UTILS::MatrixColTransform sibtransform_;
    ///@}
  };

}

#endif
#endif
