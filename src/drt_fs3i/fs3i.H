#ifdef CCADISCRET

#ifndef FS3I_H
#define FS3I_H

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_fsi/fsi_monolithic_nox.H"
#include "../drt_fsi/fsi_matrixtransform.H"

namespace FS3I
{
  class FS3I_Base
  {

  public:

    /// constructor of base class
    FS3I_Base(){};

    /// destructor of base class
    virtual ~FS3I_Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) = 0;

    Teuchos::RCP<Epetra_Vector> Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv);

    Teuchos::RCP<Epetra_Vector> Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv);

    /// read restart
    void ReadRestart();

    /// make sure potential Dirichlet conditions at the scatra coupling
    /// interface are defined on both discretizations
    void CheckInterfaceDirichletBC();

    /// extract fluid convective and structure convective velocities
    void ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel,
                    std::vector<Teuchos::RCP<const Epetra_Vector> >& convel);

    void SetVelocityFields();

    void PrepareTimeStep();

    void SetMeshDisp();

    void DoFsiStep();

    bool ConvergenceCheck(int itnum);

    /// output of scalars and mean scalars
    void ScatraOutput();

    /// fsi algorithm
    Teuchos::RCP< ::FSI::MonolithicNOX> fsi_;

    /// vector of scatra algorithms
    std::vector<Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> > scatravec_;

    /// coupling of dofs at the scatra interface
    ADAPTER::Coupling scatracoup_;

    /// scatra rhs vector
    Teuchos::RCP<Epetra_Vector> scatrarhs_;

    /// scatra increment vector
    Teuchos::RCP<Epetra_Vector> scatraincrement_;

  };
}

#endif
#endif
