#ifdef CCADISCRET

#ifndef FSI_MONOLITHIC_H
#define FSI_MONOLITHIC_H

#include "fsi_algorithm.H"

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_mapextractor.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>

namespace FSI
{

  /// base class of all monolithic FSI algorithms
  /*!

    Monolithic FSI is a Netwon solver on a block matrix with field blocks.

    \author u.kue
    \date 02/08
   */
  class Monolithic : public Algorithm,
                     public NOX::Epetra::Interface::Required,
                     public NOX::Epetra::Interface::Jacobian,
                     public NOX::Epetra::Interface::Preconditioner
  {
  public:
    Monolithic(Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Teuchos::RCP<Epetra_Vector> f) const = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat) const = 0;

    //@}

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return mat_; }

    /// Evaluate all fields at position x
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> x) const;

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) = 0;

  protected:

    /// setup solver for global block system
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils) = 0;

    /// setup of NOX convergance tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp) = 0;

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) const = 0;

    //! @name Access methods for subclasses

    Teuchos::RCP<NOX::Utils> Utils() const { return utils_; }

    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //@}

  private:

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    /// block matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat_;

    /// output utilities
    Teuchos::RCP<NOX::Utils> utils_;
  };
}


#endif
#endif
