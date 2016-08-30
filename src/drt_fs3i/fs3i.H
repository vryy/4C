/*!----------------------------------------------------------------------
\file fs3i.H
\brief H-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I) and fluid-porous-structure-scalar-scalar
       interaction (FPS3I).

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


*----------------------------------------------------------------------*/

#ifndef FS3I_H
#define FS3I_H


#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


// forward declarations
namespace ADAPTER
{
  class Coupling;
  class ScaTraBaseAlgorithm;
}

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace LINALG
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}

namespace FS3I
{
  class FS3I_Base
  {

  public:

    /// constructor of base class
    FS3I_Base();

    /// destructor of base class
    virtual ~FS3I_Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) = 0;

    /// read restart
    virtual void ReadRestart() = 0;

    /// needed for redistribution of FPS3I interface, if running on parallel
    virtual void RedistributeInterface() = 0;

    //! make sure potential Dirichlet conditions at the scatra coupling
    //! interface are defined on both discretizations
    void CheckInterfaceDirichletBC();

    //! Check FS3I specific inputs
    void CheckFS3IInputs();

    //! output of scalars and mean scalars
    void ScatraOutput();

    //! increment step and time
    void IncrementTimeAndStep();

    //! update ScaTra solution vectors (new time step)
    void UpdateScatraFields();

    //! evaluate, solve and iteratively update coupled ScaTra problem
    void ScatraEvaluateSolveIterUpdate();

    //! @name monolithic ScaTra problem
    //@{

    //! evaluate ScaTra fields
    virtual void EvaluateScatraFields();

    //! set Membrane concentration in scatra fields
    void SetMembraneConcentration() const;

    //! set-up of global matrix and rhs of the monolithic ScaTra problem
    void SetupCoupledScatraSystem();

    //! set-up of global rhs of the monolithic ScaTra problem
    void SetupCoupledScatraVector
    (
      Teuchos::RCP<Epetra_Vector>  globalvec, //!< resulting global vector
      Teuchos::RCP<const Epetra_Vector>& vec1,//!< vector in fluid ScaTra map
      Teuchos::RCP<const Epetra_Vector>& vec2 //!< vector in solid ScaTra map
    );

    //! set-up of global rhs of the monolithic ScaTra problem
    void SetupCoupledScatraRHS();

    //! set-up of global matrix of the monolithic ScaTra problem
    void SetupCoupledScatraMatrix();

    Teuchos::RCP<Epetra_Vector> Scatra2ToScatra1 (Teuchos::RCP<const Epetra_Vector> iv) const;

    Teuchos::RCP<Epetra_Vector> Scatra1ToScatra2 (Teuchos::RCP<const Epetra_Vector> iv) const;

    //! linear solution of monolithic ScaTra problem
    void LinearSolveScatra();

    //! iterative update of ScaTra solution vectors
    void ScatraIterUpdate();

    //! extraction of field-specific vectors from global ScaTra vector
    void ExtractScatraFieldVectors
    (
      Teuchos::RCP<const Epetra_Vector>  globalvec, //!< global vector
      Teuchos::RCP<const Epetra_Vector>& vec1,      //!< resulting vector in fluid ScaTra map
      Teuchos::RCP<const Epetra_Vector>& vec2       //!< resulting vector in solid ScaTra map
    );

  private:

    ///extracts membrane concentration in membrane (interface)
    void ExtractMembraneConcentration(std::vector<Teuchos::RCP<Epetra_Vector> >& MeanConc) const;

    ///Calculation of membane concentration in the membrane between fluid-scatra and structure-scatra
    Teuchos::RCP<Epetra_Vector> CalcMembraneConcentration() const;

  protected:
    /// vector of scatra algorithms
    std::vector<Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> > scatravec_;

    /// scatra rhs vector
    Teuchos::RCP<Epetra_Vector> scatrarhs_;

    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> scatraincrement_;

    /// dof row map of scatra problems splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> scatraglobalex_;

    /// vector of scatra field map extractors (coupled vs. uncoupled dofs)
    std::vector<Teuchos::RCP<LINALG::MultiMapExtractor> > scatrafieldexvec_;

    /// coupling of dofs at the scatra interface
    Teuchos::RCP<ADAPTER::Coupling> scatracoup_;

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrasystemmatrix_;

    /// coupling forces (in case of surface permeability)
    std::vector<Teuchos::RCP<Epetra_Vector> > scatracoupforce_;

    /// coupling matrices (in case of surface permeability)
    std::vector<Teuchos::RCP<LINALG::SparseMatrix> > scatracoupmat_;

    /// zero vector (needed for application of Dirichlet BC on coupling vector)
    std::vector<Teuchos::RCP<Epetra_Vector> > scatrazeros_;

    /// scatra solver
    Teuchos::RCP<LINALG::Solver> scatrasolver_;

    /// flag for infinite surface permeability
    const bool infperm_;

    /// @name  control parameters for time-integration scheme

    /// maximal simulation time
    const double timemax_;

    /// number of steps to simulate
    const int numstep_;

    /// timestep
    const double dt_;

    /// current time
    double time_;

    ///current step
    int step_;

    //@}
  private:

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> sbbtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> sbitransform_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> sibtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> fbitransform_;
    ///@}

  };
}

#endif
