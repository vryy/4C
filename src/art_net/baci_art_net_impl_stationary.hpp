/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network stationary formulation.


\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_IMPL_STATIONARY_HPP
#define FOUR_C_ART_NET_IMPL_STATIONARY_HPP

#include "baci_config.hpp"

#include "baci_art_net_timint.hpp"

FOUR_C_NAMESPACE_OPEN

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
        IO::DiscretizationWriter& output);


    // initialization
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    // test results
    void TestResults() override;

    // create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override;

    /// setup the variables to do a new time step
    void TimeUpdate() override;

    /// prepare time step
    void PrepareTimeStep() override;

    /// setup Dirichlet Boundary conditions
    void ApplyDirichletBC();

    /// reset artery diameter of previous time step
    void ResetArteryDiamPreviousTimeStep();

    //! Apply Neumann boundary conditions
    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    /// add neumann BC to residual
    void AddNeumannToResidual();

    /// initialization
    void InitSaveState() override
    {
      dserror("InitSaveState() not available for stationary formulation");
    }

    // restart
    void ReadRestart(int step, bool CoupledTo3D = false) override;

    /// save state
    void SaveState() override { dserror("SaveState() not available for stationary formulation"); }

    void LoadState() override { dserror("LoadState() not available for stationary formulation"); }

    // output
    void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) override;

    //! output of element radius
    void OutputRadius();

    //! output of element volumetric flow
    void OutputFlow();

    //! set the initial field on the artery discretization
    void SetInitialField(const INPAR::ARTDYN::InitialField init,  //!< type of initial field
        const int startfuncno                                     //!< number of spatial function
        ) override;

    // prepare the loop
    void PrepareTimeLoop() override;

    // solve artery system of equation
    void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) override;

    // prepare linear solve (apply DBC)
    void PrepareLinearSolve() override;

    // Assembling of the RHS Vector and the LHS Matrix
    void AssembleMatAndRHS() override;

    // Solve the Linear System of equations
    void LinearSolve();

    // Solve Scatra equations
    void SolveScatra() override;

    // get solution vector = pressure
    Teuchos::RCP<const Epetra_Vector> Pressurenp() const override { return pressurenp_; }

    //! get element volume flow
    Teuchos::RCP<const Epetra_Vector> EleVolflow() const { return ele_volflow_; }

    //! get element radius
    Teuchos::RCP<const Epetra_Vector> EleRadius() const { return ele_radius_; }

    //! iterative update of primary variable
    void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc) override
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
    //! volumetric flow (for output)
    Teuchos::RCP<Epetra_Vector> ele_volflow_;
    //! element radius (for output)
    Teuchos::RCP<Epetra_Vector> ele_radius_;
    /// underlying scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

  };  // class ArtNetImplStationary
}  // namespace ART


FOUR_C_NAMESPACE_CLOSE

#endif
