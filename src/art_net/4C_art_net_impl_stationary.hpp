/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network stationary formulation.


\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_IMPL_STATIONARY_HPP
#define FOUR_C_ART_NET_IMPL_STATIONARY_HPP

#include "4C_config.hpp"

#include "4C_art_net_timint.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class ScaTraBaseAlgorithm;
}

namespace Arteries
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
    ArtNetImplStationary(Teuchos::RCP<Core::FE::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Core::IO::DiscretizationWriter& output);


    // initialization
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    // test results
    void TestResults() override;

    // create field test
    Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() override;

    /// setup the variables to do a new time step
    void TimeUpdate() override;

    /// prepare time step
    void prepare_time_step() override;

    /// setup Dirichlet Boundary conditions
    void apply_dirichlet_bc();

    /// reset artery diameter of previous time step
    void reset_artery_diam_previous_time_step();

    //! Apply Neumann boundary conditions
    void apply_neumann_bc(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    /// add neumann BC to residual
    void add_neumann_to_residual();

    /// initialization
    void InitSaveState() override
    {
      FOUR_C_THROW("InitSaveState() not available for stationary formulation");
    }

    // restart
    void read_restart(int step, bool CoupledTo3D = false) override;

    /// save state
    void SaveState() override
    {
      FOUR_C_THROW("SaveState() not available for stationary formulation");
    }

    void LoadState() override
    {
      FOUR_C_THROW("LoadState() not available for stationary formulation");
    }

    // output
    void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) override;

    //! output of element radius
    void OutputRadius();

    //! output of element volumetric flow
    void OutputFlow();

    //! set the initial field on the artery discretization
    void SetInitialField(const Inpar::ArtDyn::InitialField init,  //!< type of initial field
        const int startfuncno                                     //!< number of spatial function
        ) override;

    // prepare the loop
    void prepare_time_loop() override;

    // solve artery system of equation
    void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) override;

    // prepare linear solve (apply DBC)
    void PrepareLinearSolve() override;

    // Assembling of the RHS Vector and the LHS Matrix
    void assemble_mat_and_rhs() override;

    // Solve the Linear System of equations
    void linear_solve();

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
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra_;

  };  // class ArtNetImplStationary
}  // namespace Arteries


FOUR_C_NAMESPACE_CLOSE

#endif
