/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for implicit time integration schemes for
        multiphas porous fluid problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_TIMINT_IMPLICIT_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_TIMINT_IMPLICIT_HPP



#include "4C_config.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_adapter_porofluidmultiphase.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_inpar_porofluidmultiphase.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace Discret
{
  class ResultTest;
}  // namespace Discret

namespace Core::DOFSets
{
  class DofSet;
}  // namespace Core::DOFSets


namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
  class KrylovProjector;
}  // namespace Core::LinAlg

namespace Adapter
{
  class ArtNet;
}


namespace POROFLUIDMULTIPHASE
{
  // forward declaration
  class MeshtyingStrategyBase;

  /*!
   * \brief implicit time integration for porous multiphase flow problems
   */

  class TimIntImpl : public Adapter::PoroFluidMultiphase
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntImpl(Teuchos::RCP<Core::FE::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! initialize time integration
    void init(bool isale, int nds_disp, int nds_vel, int nds_solidpressure, int nds_scalar,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! add global state vectors specific for time-integration scheme
    virtual void add_time_integration_specific_vectors() = 0;

    //! prepare time loop
    void prepare_time_loop() override;

    //! setup the variables to do a new time step
    void prepare_time_step() override;

    //! initialization procedure prior to evaluation of first time step
    virtual void prepare_first_time_step();

    //! initialization procedure prior to evaluation of first time step
    virtual void calc_initial_time_derivative();

    //! read restart data
    void read_restart(int step) override;

    /// create result test for porous fluid field
    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override;

    //! finite difference check for system matrix
    void fd_check();

    /*--- calculate and update -----------------------------------------------*/

    //! do time integration (time loop)
    void time_loop() override;

    //! general solver call for coupled algorithms
    void solve() override;

    //! update the solution after convergence of the nonlinear iteration.
    void update() override;

    ///  compute time derivative
    virtual void compute_time_derivative() = 0;

    ///  compute intermediate values if necessary
    virtual void compute_intermediate_values() = 0;

    //! apply moving mesh data
    void apply_mesh_movement(Teuchos::RCP<const Epetra_Vector> dispnp  //!< displacement vector
        ) override;

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    void set_velocity_field(Teuchos::RCP<const Epetra_Vector> vel  //!< velocity vector
        ) override;

    //! set state on discretization
    void set_state(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state) override;

    //! calculate error compared to analytical solution
    void evaluate_error_compared_to_analytical_sol() override;

    /*--- query and output ---------------------------------------------------*/

    //! print information about current time step to screen
    virtual void print_time_step_info();

    //! iterative update of phinp
    void update_iter(const Teuchos::RCP<const Epetra_Vector> inc  //!< increment vector for phi
        ) override;

    //! build linear system tangent matrix, rhs/force residual
    void evaluate() override;

    //! apply Dirichlet Boundary Condition
    void prepare_system_for_newton_solve();

    //! direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_);
    };

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() const override
    {
      return dbcmaps_with_volfracpress_;
    }

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> rhs() const override { return residual_; }

    //! right-hand side alias the dynamic force residual for coupled system
    Teuchos::RCP<const Epetra_Vector> artery_porofluid_rhs() const override;

    //! return discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() const override { return discret_; }

    //! access dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) const override;

    //! access dof row map
    Teuchos::RCP<const Epetra_Map> artery_dof_row_map() const override;

    //! direct access to block system matrix of artery poro problem
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat() const override;

    //! output solution and restart data to file
    void output() override;

    //! output solution and restart data to file
    virtual void print_header();

    //! write additional data required for restart
    virtual void output_restart() = 0;

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- query and output ---------------------------------------------------*/

    //! return current time value
    double time() const { return time_; }

    //! return current step number
    int step() const { return step_; }

    //! return number of newton iterations in last timestep
    double iter_num() const { return iternum_; }

    //! return time step size
    double dt() const { return dt_; }

    /*========================================================================*/
    //! @name degrees of freedom and related
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! set the initial scalar field phi
    virtual void set_initial_field(
        const Inpar::POROFLUIDMULTIPHASE::InitialField init,  //!< type of initial field
        const int startfuncno                                 //!< number of spatial function
    );

    /*--- query and output ---------------------------------------------------*/

    //! return pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> phinp() const override { return phinp_; }

    //! return scalar field phi at time n
    Teuchos::RCP<const Epetra_Vector> phin() const override { return phin_; }

    //! return time derivative of scalar field phi at time n
    Teuchos::RCP<const Epetra_Vector> phidtn() const { return phidtn_; }

    //! return time derivative of scalar field phi at time n+1
    Teuchos::RCP<const Epetra_Vector> phidtnp() const { return phidtnp_; }

    //! return scalar field history
    Teuchos::RCP<const Epetra_Vector> hist() const { return hist_; }

    //! return solid pressure field
    Teuchos::RCP<const Epetra_Vector> solid_pressure() const override
    {
      if (!output_solidpress_)
        FOUR_C_THROW("solid pressure requested but flag OUTPUT_SOLIDPRESS set to no");
      return solidpressure_;
    }

    //! return pressure field
    Teuchos::RCP<const Epetra_Vector> pressure() const override
    {
      if (!output_satpress_)
        FOUR_C_THROW("pressure requested but flag OUTPUT_SATANDPRESS set to no");
      return pressure_;
    }

    //! return saturation field
    Teuchos::RCP<const Epetra_Vector> saturation() const override
    {
      if (!output_satpress_)
        FOUR_C_THROW("saturation requested but flag OUTPUT_SATANDPRESS set to no");
      return saturation_;
    }

    //! return phase flux field at time n+1
    Teuchos::RCP<const Epetra_MultiVector> flux() const override { return flux_; }

    //! return phase velocity at time n+1
    Teuchos::RCP<const Epetra_MultiVector> phase_velocity() const { return phase_velocities_; }

    //! return number of dof set associated with solid pressure
    int get_dof_set_number_of_solid_pressure() const override { return nds_solidpressure_; };

    //! return valid volume fraction species
    Teuchos::RCP<const Epetra_Vector> valid_vol_frac_spec_dofs() const override
    {
      return valid_volfracspec_dofs_;
    }

    //! return number of domain integral functions
    int num_domain_int_functions() const { return num_domainint_funct_; }

    //! return the values of the domain integrals
    Teuchos::RCP<const Core::LinAlg::SerialDenseVector> domain_int_values() const
    {
      return domain_integrals_;
    }

    //! return the meshtying strategy
    Teuchos::RCP<POROFLUIDMULTIPHASE::MeshtyingStrategyBase> mesh_tying_strategy() const
    {
      return strategy_;
    }

   protected:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! don't want copy constructor
    TimIntImpl(const TimIntImpl& old);

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! Set element time step parameters (varying every time step)
    virtual void set_element_time_step_parameter() const = 0;

    //! set time for evaluation of Neumann boundary conditions
    virtual void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) = 0;

    //! Set general element parameters
    void set_element_general_parameters() const;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! Set the part of the righthandside belonging to the last timestep.
    virtual void set_old_part_of_righthandside() = 0;

    /*--- calculate and update -----------------------------------------------*/

    //! Apply Dirichlet boundary conditions on provided state vector
    void apply_dirichlet_bc(const double time,  //!< evaluation time
        Teuchos::RCP<Epetra_Vector> prenp,      //!< pressure (may be = null)
        Teuchos::RCP<Epetra_Vector> predt       //!< first time derivative (may be = null)
    );

    //! potential residual scaling and potential addition of Neumann terms
    void scaling_and_neumann();

    //! add actual Neumann loads multipl. with time factor to the residual
    virtual void add_neumann_to_residual() = 0;

    //! Apply Neumann boundary conditions
    void apply_neumann_bc(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    //! call elements to calculate system matrix and rhs and assemble
    virtual void assemble_mat_and_rhs();

    //! call elements to find the valid volume frac pressures and species
    virtual void evaluate_valid_volume_frac_press_and_spec();

    //! apply the additional volume fraction pressures as DBC
    virtual void apply_additional_dbc_for_vol_frac_press();

    //! apply the starting Dirichlet boundary condition
    virtual void apply_starting_dbc();

    //! call elements to calculate fluid coupling matrix with structure and assemble
    void assemble_fluid_struct_coupling_mat(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs) override;

    //! call elements to calculate fluid coupling matrix with scatra and assemble
    void assemble_fluid_scatra_coupling_mat(
        Teuchos::RCP<Core::LinAlg::SparseOperator> k_pfs) override;

    //! return the right time-scaling-factor for the true residual
    virtual double residual_scaling() const = 0;

    //! contains the nonlinear iteration loop
    virtual void nonlinear_solve();

    //! check convergence (or divergence) of nonlinear iteration
    bool abort_nonlin_iter(const int itnum,  //!< current value of iteration step counter
        const int itemax,                    //!< maximum number of iteration steps
        const double abstolres,              //!< absolute tolerance for the residual norm
        double& actresidual                  //!< return value of the current residual
    );

    //! linear solve
    virtual void linear_solve(bool isadapttol, double actresidual, double adaptolbetter);

    //! reconstruct pressures and saturation from current solution
    void reconstruct_pressures_and_saturations() override;

    //! reconstruct solid pressures from current solution
    void reconstruct_solid_pressures();

    //! reconstruct fluxes from current solution
    void reconstruct_flux() override;

    //! calculate phase velocities from current solution
    void calculate_phase_velocities() override;

    //! reconstruct porosity from current solution
    void reconstruct_porosity();

    //! evaluate domain integrals
    void evaluate_domain_integrals();

    /*--- query and output ---------------------------------------------------*/

    //! is output needed for the current time step?
    bool do_output() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    //! write state vectors prenp to BINIO
    virtual void output_state();

    //! print header of convergence table to screen
    virtual void print_convergence_header();

    //! print first line of convergence table to screen
    virtual void print_convergence_values_first_iter(
        const int& itnum,                      //!< current Newton-Raphson iteration step
        const int& itemax,                     //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,                   //!< relative tolerance for Newton-Raphson scheme
        const std::vector<double>& preresnorm  //!< L2 norm of pressure residual
    );

    //! print current line of convergence table to screen
    virtual void print_convergence_values(
        const int& itnum,     //!< current Newton-Raphson iteration step
        const int& itemax,    //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,  //!< relative tolerance for Newton-Raphson scheme
        const std::vector<double>& preresnorm,  //!< norm of pressure residual
        const std::vector<double>& incprenorm,  //!< norm of pressure increment
        const std::vector<double>& prenorm      //!< norm of pressure state vector
    );

    //! print finish line of convergence table to screen
    virtual void print_convergence_finish_line();

    // return arterial network time integrator
    Teuchos::RCP<Adapter::ArtNet> art_net_tim_int() override;

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! increment time and step value
    void increment_time_and_step();

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! linear solver
    Teuchos::RCP<Core::LinAlg::Solver> solver_;

    //! solver number in input file
    const int linsolvernumber_;

    //! parameter list of global control problem
    const Teuchos::ParameterList& params_;

    //! parameter list of poro fluid multiphase problem
    const Teuchos::ParameterList& poroparams_;

    //! processor id
    int myrank_;

    //! number of space dimensions
    int nsd_;

    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for Eulerian or ALE formulation of equation(s)
    bool isale_;

    //! flag if initial time derivative should be skipped
    bool skipinitder_;

    //! flag if saturations and pressures should be output
    bool output_satpress_;

    //! flag if solid pressure should be output
    bool output_solidpress_;

    //! flag if porosity should be output
    bool output_porosity_;

    //! flag if phase velocities should be written to output
    bool output_phase_velocities_;

    //! flag if blood vessel volume fraction should be output (for 1D-3D coupling)
    bool output_bloodvesselvolfrac_;

    //! flag for biot stabilization
    bool stab_biot_;

    /*--- query and output ---------------------------------------------------*/

    //! parameters for domain integration
    std::vector<int> domainint_funct_;
    int num_domainint_funct_;

    //! values of domain integrals
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> domain_integrals_;

    //! flag for error calculation
    const Inpar::POROFLUIDMULTIPHASE::CalcError calcerr_;

    //! flag for flux reconstruction
    const Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod fluxrecon_;

    //! solver number for flux reconstruction
    const int fluxreconsolvernum_;

    //! what to do when nonlinear solution fails
    enum Inpar::POROFLUIDMULTIPHASE::DivContAct divcontype_;

    //! flag for finite difference check
    const Inpar::POROFLUIDMULTIPHASE::FdCheck fdcheck_;

    //! perturbation magnitude for finite difference check
    const double fdcheckeps_;

    //! relative tolerance for finite difference check
    const double fdchecktol_;

    //! scaling factor for biot stabilization
    double stab_biot_scaling_;

    /*========================================================================*/
    //! @name Time, time-step, and iteration variables
    /*========================================================================*/

    //! actual time
    double time_;

    //! maximum simulation time
    double maxtime_;

    //! actual step number
    int step_;

    //! maximum number of steps
    const int stepmax_;

    //! time step size
    double dt_;

    //! time measurement element
    double dtele_;

    //! time measurement solve
    double dtsolve_;

    //! number of newton iterations in actual timestep
    int iternum_;

    //! maximum number of newton iterations
    const int itemax_;

    //! write results every upres_ steps ? writesolutionevery_
    const int upres_;

    //! write restart data every uprestart_ steps ? writesolutioneveryrestart_
    const int uprestart_;

    // vector norm for residuals
    enum Inpar::POROFLUIDMULTIPHASE::VectorNorm vectornormfres_;
    // vector norm for increments
    enum Inpar::POROFLUIDMULTIPHASE::VectorNorm vectornorminc_;

    //! convergence tolerance for increments
    double ittolres_;
    //! convergence tolerance for residuals
    double ittolinc_;

    //! flag if artery coupling is active
    bool artery_coupling_active_;

    /*========================================================================*/
    //! @name degrees of freedom variables
    /*========================================================================*/

    //! phi at time n
    Teuchos::RCP<Epetra_Vector> phin_;
    //! phi at time n+1
    Teuchos::RCP<Epetra_Vector> phinp_;

    //! time derivative of phi at time n
    Teuchos::RCP<Epetra_Vector> phidtn_;
    //! time derivative of phi at time n+1
    Teuchos::RCP<Epetra_Vector> phidtnp_;

    //! histvector --- a linear combination of phinm, phin (BDF)
    //!                or phin, phidtn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector> hist_;

    /*========================================================================*/
    //! @name degrees of freedom and related
    /*========================================================================*/

    //! pressure at time n+1
    Teuchos::RCP<Epetra_Vector> pressure_;

    //! saturation at time n+1
    Teuchos::RCP<Epetra_Vector> saturation_;

    //! solid pressure at time n+1
    Teuchos::RCP<Epetra_Vector> solidpressure_;

    //! porosity at time n+1
    Teuchos::RCP<Epetra_Vector> porosity_;

    //! vector with valid volume fraction pressure dofs, this vector identifies volume fraction
    //! pressure DOFs,
    //  which actually have to be evaluated with a double >= 1.0, see also
    //  EvaluatorValidVolFracPressures: if at least one nodal volume fraction value of an element is
    //  bigger than a threshold (min volfrac), the volume fraction pressure is a valid (physically
    //  meaningful) quantity in this element and the respective Darcy equation has to be solved
    //  for volume fraction species we only evaluate if all nodal volume fraction values of the
    //  element are bigger than the threshold (min volfrac), this turned out to be the most stable
    //  approach
    Teuchos::RCP<Epetra_Vector> valid_volfracpress_dofs_;
    Teuchos::RCP<Epetra_Vector> valid_volfracspec_dofs_;

    //! flux of each phase at time n+1 (post-processed from pressure solution)
    Teuchos::RCP<Epetra_MultiVector> flux_;

    //! velocity of each phase at time n+1 (post-processed from pressure solution)
    Teuchos::RCP<Epetra_MultiVector> phase_velocities_;

    //! number of dofset associated with displacement dofs
    int nds_disp_;

    //! number of dofset associated with velocity related dofs
    int nds_vel_;

    //! number of dofset associated with solid pressure dofs
    int nds_solidpressure_;

    //! number of dofset associated with scatra dofs
    int nds_scatra_;

    /*========================================================================*/
    //! @name Galerkin discretization, boundary conditions, and related
    /*========================================================================*/

    //! the porous multiphase flow discretization
    Teuchos::RCP<Core::FE::Discretization> discret_;

    //! the discretization writer
    Teuchos::RCP<Core::IO::DiscretizationWriter> output_;

    //! system matrix (either sparse matrix or block sparse matrix)
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat_;

    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;

    //! maps for extracting Dirichlet and free DOF sets, here the additional dofs have been added
    //! which have to be zeroed out for the volume fraction pressure since it is not defined if the
    //! corresponding volume fraction is equal to zero (or smaller than minvolfrac)
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_with_volfracpress_;

    //! maps for extracting Dirichlet and free DOF sets with additional starting Dirichlet boundary
    //! condition
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_starting_condition_;

    //! the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;

    //! residual vector
    Teuchos::RCP<Epetra_Vector> residual_;

    //! true (rescaled) residual vector without zeros at Dirichlet conditions
    Teuchos::RCP<Epetra_Vector> trueresidual_;

    //! nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector> increment_;

    //! meshtying strategy (includes standard case without meshtying)
    Teuchos::RCP<POROFLUIDMULTIPHASE::MeshtyingStrategyBase> strategy_;

    //! end time point when to switch off the starting Dirichlet boundary condition
    double starting_dbc_time_end_;

    //! switch the starting Dirichlet boundary condition on or off for the different DOFs
    std::vector<bool> starting_dbc_onoff_;

    //! function ID prescribing the starting Dirichlet boundary condition
    std::vector<int> starting_dbc_funct_;

    /*========================================================================*/

  };  // class TimIntImpl
}  // namespace POROFLUIDMULTIPHASE


FOUR_C_NAMESPACE_CLOSE

#endif
