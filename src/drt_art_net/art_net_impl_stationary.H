/*!----------------------------------------------------------------------
\brief Control routine for arterial network stationary formulation.

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/

#ifndef SRC_DRT_ART_NET_ART_NET_IMPL_STATIC_H_
#define SRC_DRT_ART_NET_ART_NET_IMPL_STATIC_H_

#include "art_net_timint.H"

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace ART
{
  /*!
  \brief stationary formulation for arterial network problems

  \author kremheller
  */
  class ArtNetImplStationary : public TimInt
  {
   public:
    /*!
    \brief Standard Constructor

    */
    ArtNetImplStationary(Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        FILE* errfile, IO::DiscretizationWriter& output);

    // destructor
    virtual ~ArtNetImplStationary();

    // initialization
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname);

    // test results
    void TestResults();

    // create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// setup the variables to do a new time step
    void TimeUpdate();

    /// prepare time step
    void PrepareTimeStep();

    /// setup Dirichlet Boundary conditions
    void ApplyDirichletBC();

    //! Apply Neumann boundary conditions
    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    /// add neumann BC to residual
    void AddNeumannToResidual();

    /// initialization
    void InitSaveState() { dserror("InitSaveState() not available for stationary formulation"); }

    // restart
    void ReadRestart(int step, bool CoupledTo3D = false);

    /// save state
    void SaveState() { dserror("SaveState() not available for stationary formulation"); }

    void LoadState() { dserror("LoadState() not available for stationary formulation"); }

    // output
    void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams);

    //! output of element radius
    void OutputRadius();

    //! set the initial field on the artery discretization
    virtual void SetInitialField(const INPAR::ARTDYN::InitialField init,  //!< type of initial field
        const int startfuncno  //!< number of spatial function
    );

    // prepare the loop
    virtual void PrepareTimeLoop();

    // solve artery system of equation
    void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    // prepare linear solve (apply DBC)
    void PrepareLinearSolve();

    // Assembling of the RHS Vector and the LHS Matrix
    void AssembleMatAndRHS();

    // Solve the Linear System of equations
    void LinearSolve();

    // Solve Scatra equations
    void SolveScatra();

    // get solution vector = pressure
    Teuchos::RCP<const Epetra_Vector> Pressurenp() const { return pressurenp_; }

    //! iterative update of primary variable
    virtual void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
    {
      pressurenp_->Update(1.0, *inc, 1.0);
      return;
    }


   private:
    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;
    //! pressure at time n+1
    Teuchos::RCP<Epetra_Vector> pressurenp_;
    //! pressure increment at time n+1
    Teuchos::RCP<Epetra_Vector> pressureincnp_;
    //! the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;
    /// underlying scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

  };  // class ArtNetImplStationary
}  // namespace ART


#endif /* SRC_DRT_ART_NET_ART_NET_IMPL_STATIC_H_ */
