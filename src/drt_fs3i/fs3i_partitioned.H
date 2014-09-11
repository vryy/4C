/*!----------------------------------------------------------------------
\file fs3i_partitioned.H
\brief H-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I), that is, algorithmic routines not specifically
       related to partitioned solution approaches to one -or
       two-way-coupled problem configurations, respectively

<pre>
Maintainers: Lena Yoshihara & Volker Gravemeier
             {yoshihara,vgravem}@lnm.mw.tum.de
             089/289-15303,-15245
</pre>

*----------------------------------------------------------------------*/


#ifndef PARTFS3I_H
#define PARTFS3I_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "fs3i.H"

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
  class PartFS3I : public FS3I_Base
  {

  public:

    //! constructor of base class for partitioned FS3I
    PartFS3I(const Epetra_Comm& comm);

    //! destructor of base class for partitioned FS3I
    virtual ~PartFS3I(){};

    //! @name overall FS3I system
    //@{

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    virtual void Timeloop() = 0;

    //! flag whether time loop should be finished
    bool NotFinished() { return step_ < numstep_ and time_ <= timemax_; };

    //! increment step and time
    void IncrementTimeAndStep();

    //! read and set fields needed for restart
    virtual void ReadRestart();

    //! set-up of FSI and ScaTra systems
    virtual void SetupSystem();

    //! test results for individual fields
    virtual void TestResults(const Epetra_Comm& comm);

    //! information transfer FSI -> ScaTra
    void SetFSISolution();

    /// set scatra solution on structure field
    void SetStructScatraSolution() ;

    //! evalute, solver and iteratively update coupled ScaTra problem
    void ScatraEvaluateSolveIterUpdate();

    //@}


    //! @name monolithic ScaTra problem
    //@{

    //! evaluate ScaTra fields
    void EvaluateScatraFields();

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

    //! linear solution of monolithic ScaTra problem
    void LinearSolveScatra();

    //! extraction of field-specific vectors from global ScaTra vector
    void ExtractScatraFieldVectors
    (
      Teuchos::RCP<const Epetra_Vector>  globalvec, //!< global vector
      Teuchos::RCP<const Epetra_Vector>& vec1,      //!< resulting vector in fluid ScaTra map
      Teuchos::RCP<const Epetra_Vector>& vec2       //!< resulting vector in solid ScaTra map
    );

    //! iterative update of ScaTra solution vectors
    void ScatraIterUpdate();

    //! check convergence of monolithic ScaTra problem (depends on which
    //coupling is considered)
    virtual bool ScatraConvergenceCheck(int itnum) = 0;

    //! update ScaTra solution vectors (new time step)
    void UpdateScatraFields();

    //! return communicator
    const Epetra_Comm& Comm() const { return comm_; }


    Teuchos::RCP<Epetra_Vector> Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv);

    Teuchos::RCP<Epetra_Vector> Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv);

    //! make sure potential Dirichlet conditions at the scatra coupling
    //! interface are defined on both discretizations
    void CheckInterfaceDirichletBC();

    /// extract fluid convective and structure convective velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel,
                    std::vector<Teuchos::RCP<const Epetra_Vector> >& convel);

    void SetVelocityFields();

    //! routine for preparing time step to be defined in inherited classes
    //! (structure depends on coupling, that is, either one- or two-way)
    virtual void PrepareTimeStep() = 0;

    void SetMeshDisp();

    //! output of scalars and mean scalars
    void ScatraOutput();
    //@}

  protected:

    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra algorithms
    std::vector<Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> > scatravec_;

    /// scatra rhs vector
    Teuchos::RCP<Epetra_Vector> scatrarhs_;

    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> scatraincrement_;

    /// dof row map of scatra problems splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> scatraglobalex_;

    //! @name  control parameters for time-integration scheme
    int step_;
    int numstep_;

    double dt_;
    double time_;
    double timemax_;
    //@}

  private:

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of dofs at the scatra interface
    Teuchos::RCP<ADAPTER::Coupling> scatracoup_;

    /// vector of scatra field map extractors (coupled vs. uncoupled dofs)
    std::vector<Teuchos::RCP<LINALG::MultiMapExtractor> > scatrafieldexvec_;

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

    /// flag for infinite surface permeability
    bool infperm_;

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
