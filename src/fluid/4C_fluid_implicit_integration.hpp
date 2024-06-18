/*-----------------------------------------------------------*/
/*! \file

\brief Control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o a one-step-theta time-integration scheme,

     o a two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm),

     o two variants of a generalized-alpha time-integration scheme

     and a stationary solver.


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_IMPLICIT_INTEGRATION_HPP
#define FOUR_C_FLUID_IMPLICIT_INTEGRATION_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fluid_timint.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_turbulence.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
  class DiscretizationFaces;
}  // namespace Core::FE
namespace Core::IO
{
  class DiscretizationWriter;
}
namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace Core::Conditions
{
  class LocsysManager;
}


namespace Adapter
{
  class CouplingMortar;
}

namespace FLD
{
  // forward declarations
  class TurbulenceStatisticManager;
  class ForcingInterface;
  class DynSmagFilter;
  class Vreman;
  class Boxfilter;
  class Meshtying;
  class XWall;
  class TransferTurbulentInflowCondition;
  namespace UTILS
  {
    class FluidInfNormScaling;
    class FluidImpedanceWrapper;
    class StressManager;
  }  // namespace UTILS

  /*!
  \brief time integration for fluid problems

  */
  class FluidImplicitTimeInt : public TimInt
  {
    friend class TurbulenceStatisticManager;
    friend class HomIsoTurbInitialField;
    friend class HomIsoTurbForcing;
    friend class PeriodicHillForcing;
    friend class FluidResultTest;

   public:
    /*!
    \brief Standard Constructor

    */
    FluidImplicitTimeInt(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void init() override;

    /*!
    \brief initialization of nonlinear BCs

    */
    virtual void InitNonlinearBC();

    /*!
    \brief start time loop for starting algorithm, normal problems and restarts

    */
    void Integrate() override;

    /*!
    \brief Do time integration (time loop)

    */
    virtual void TimeLoop();

    /*!
    \brief Print information about current time step to screen

    */
    virtual void print_time_step_info()
    {
      FOUR_C_THROW("you are in the base class");
      return;
    }

    /*!
    \brief Set theta_ to its value, dependent on integration method for GenAlpha and BDF2

    */
    virtual void SetTheta() { return; }

    /*!
    \brief Set the part of the righthandside belonging to the last
           timestep for incompressible or low-Mach-number flow

       for low-Mach-number flow: distinguish momentum and continuity part
       (continuity part only meaningful for low-Mach-number flow)

       Stationary/af-generalized-alpha:

                     mom: hist_ = 0.0
                    (con: hist_ = 0.0)

       One-step-Theta:

                     mom: hist_ = veln_  + dt*(1-Theta)*accn_
                    (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

       BDF2: for constant time step:

                     mom: hist_ = 4/3 veln_  - 1/3 velnm_
                    (con: hist_ = 4/3 densn_ - 1/3 densnm_)


    */
    virtual void set_old_part_of_righthandside() = 0;

    /*!
    \brief Set gamma to a value

    */
    virtual void SetGamma(Teuchos::ParameterList& eleparams) = 0;

    /*!
    \brief Initialize function which is called after that the constructor of the time integrator has
    been called

    */
    virtual void CompleteGeneralInit();

    /*!
     * \brief Create internal faces extension
     */
    void create_faces_extension();


    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    virtual void SetStateTimInt() = 0;

    /*!
    \brief Set time factor in GenAlpha

    */
    virtual double SetTimeFac() { return 1.0; }

    /*!
    \brief Scale separation

    */
    virtual void Sep_Multiply() = 0;

    /*!
    \brief Update velaf_ for GenAlpha

    */
    virtual void UpdateVelafGenAlpha() {}

    /*!
    \brief Update velam_ for GenAlpha

    */
    virtual void UpdateVelamGenAlpha() {}

    /*!
    \brief Insert Womersley condition

    */
    virtual void insert_volumetric_surface_flow_cond_vector(
        Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> res)
    {
      return;
    }

    /*!
    \brief treat turbulence models in assemble_mat_and_rhs

    */
    virtual void treat_turbulence_models(Teuchos::ParameterList& eleparams);

    /*!
    \brief Evaluate for AVM3 Separation

    */
    virtual void av_m3_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams);

    /*!
    \brief Get scale separation matrix

    */
    virtual void av_m3_get_scale_separation_matrix();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void set_custom_ele_params_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams) {}

    /*!
    \brief Call discret_->ClearState() after assembly (HDG needs to read from state vectors...)

    */
    virtual void clear_state_assemble_mat_and_rhs() { discret_->ClearState(); }

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void set_custom_ele_params_apply_nonlinear_boundary_conditions(
        Teuchos::ParameterList& eleparams)
    {
    }

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void set_custom_ele_params_linear_relaxation_solve(Teuchos::ParameterList& eleparams) {}

    /*!
    \brief Prepare calculation of acceleration

    */
    virtual void tim_int_calculate_acceleration();

    /*!
    \brief Additional function for RedModels in linear_relaxation_solve

    */
    virtual void CustomSolve(Teuchos::RCP<Epetra_Vector> relax) {}

    /*!
    \brief Call statistics manager (special case in TimIntLoma)

    */
    virtual void call_statistics_manager();

    /*!
    \brief return thermpressaf_ in TimIntLoma

    */
    virtual double ReturnThermpressaf() { return 0.0; }

    /*!
    \brief Calculate time derivatives for
           stationary/one-step-theta/BDF2/af-generalized-alpha time integration
           for incompressible and low-Mach-number flow
    */
    virtual void calculate_acceleration(
        const Teuchos::RCP<const Epetra_Vector> velnp,  ///< velocity at     n+1
        const Teuchos::RCP<const Epetra_Vector> veln,   ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector> velnm,  ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector> accn,   ///< acceleration at n-1
        const Teuchos::RCP<Epetra_Vector> accnp         ///< acceleration at n+1
        ) = 0;

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void set_element_general_fluid_parameter();

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void set_element_turbulence_parameters();

    //! @name Set general parameter in class fluid_ele_parameter_intface
    /*!

    \brief parameter (fix over all time step) are set in this method.

    */
    virtual void set_face_general_fluid_parameter();

    /// initialize vectors and flags for turbulence approach
    virtual void set_general_turbulence_parameters();

    /*!
    \brief do explicit predictor step to start nonlinear iteration from
           a better initial value
                          +-                                      -+
                          | /     dta \          dta  veln_-velnm_ |
     velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
                          | \     dtp /          dtp     dtp       |
                          +-                                      -+
    */
    virtual void explicit_predictor();

    /// setup the variables to do a new time step
    void prepare_time_step() override;

    /*!
    \brief (multiple) corrector

    */
    void Solve() override;

    /*!
  \brief solve linearised fluid

  */
    Teuchos::RCP<Core::LinAlg::Solver> LinearSolver() override { return solver_; };

    /*!
    \brief preparatives for solver

    */
    void PrepareSolve() override;

    /*!
    \brief preparations for Krylov space projection

    */
    virtual void init_krylov_space_projection();
    void setup_krylov_space_projection(Core::Conditions::Condition* kspcond) override;
    void update_krylov_space_projection() override;
    void check_matrix_nullspace() override;

    /*!
    \brief update within iteration

    */
    void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) override;

    /*!
   \brief convergence check

    */
    bool convergence_check(int itnum, int itmax, const double velrestol, const double velinctol,
        const double presrestol, const double presinctol) override;

    /*!
      \brief build linear system matrix and rhs

      Monolithic FSI needs to access the linear fluid problem.
    */
    void evaluate(Teuchos::RCP<const Epetra_Vector>
            stepinc  ///< solution increment between time step n and n+1
        ) override;


    /*!
    \brief Update the solution after convergence of the nonlinear
           iteration. Current solution becomes old solution of next
           timestep.
    */
    virtual void TimeUpdate();

    /*!
    \ time update of stresses
    */
    virtual void TimeUpdateStresses();

    virtual void time_update_nonlinear_bc();

    /*!
    \ Update of external forces

    */
    virtual void time_update_external_forces();

    /// Implement Adapter::Fluid
    void update() override { TimeUpdate(); }

    //! @name Time step size adaptivity in monolithic FSI
    //@{

    //! access to time step size of previous time step
    virtual double DtPrevious() const { return dtp_; }

    //! set time step size
    void set_dt(const double dtnew) override;

    //! set time and step
    void SetTimeStep(const double time,  ///< time to set
        const int step) override;        ///< time step number to set

    /*!
    \brief Reset time step

    In case of time step size adaptivity, time steps might have to be repeated.
    Therefore, we need to reset the solution back to the initial solution of the
    time step.

    \author mayr.mt
    \date 08/2013
    */
    void reset_step() override
    {
      accnp_->Update(1.0, *accn_, 0.0);
      velnp_->Update(1.0, *veln_, 0.0);
      dispnp_->Update(1.0, *dispn_, 0.0);

      return;
    }

    /*! \brief Reset time and step in case that a time step has to be repeated
     *
     *  Fluid field increments time and step at the beginning of a time step. If a
     *  time step has to be repeated, we need to take this into account and
     *  decrease time and step beforehand. They will be incremented right at the
     *  beginning of the repetition and, thus, everything will be fine.
     *
     *  Currently, this is needed for time step size adaptivity in FSI.
     *
     *  \author mayr.mt \date 08/2013
     */
    void reset_time(const double dtold) override { SetTimeStep(Time() - dtold, Step() - 1); }

    //! Give order of accuracy
    virtual int method_order_of_accuracy() const
    {
      return std::min(method_order_of_accuracy_vel(), method_order_of_accuracy_pres());
    }

    //! Give local order of accuracy of velocity part
    virtual int method_order_of_accuracy_vel() const
    {
      FOUR_C_THROW("Not implemented in base class. May be overwritten by derived class.");
      return 0;
    }

    //! Give local order of accuracy of pressure part
    virtual int method_order_of_accuracy_pres() const
    {
      FOUR_C_THROW("Not implemented in base class. May be overwritten by derived class.");
      return 0;
    }
    //! Return linear error coefficient of velocity
    virtual double method_lin_err_coeff_vel() const
    {
      FOUR_C_THROW("Not implemented in base class. May be overwritten by derived class.");
      return 0;
    }

    //@}

    /*!
    \brief lift'n'drag forces, statistics time sample and
           output of solution and statistics

    */
    void StatisticsAndOutput() override;

    /*!
    \brief statistics time sample and
           output of statistics

    */
    void StatisticsOutput() override;

    /*!
    \brief update configuration and output to file/screen

    */
    void Output() override;

    /*
     * \brief Write fluid runtime output
     */
    virtual void WriteRuntimeOutput();

    virtual void OutputNonlinearBC();

    virtual void output_to_gmsh(const int step, const double time, const bool inflow) const;

    /*!
    \output of external forces for restart

    */
    virtual void output_external_forces();

    /*!
    \brief get access to map extractor for velocity and pressure

    */
    Teuchos::RCP<Core::LinAlg::MapExtractor> GetVelPressSplitter() override
    {
      return velpressplitter_;
    };

    /*!
    \brief set initial flow field for analytical test problems

    */
    void SetInitialFlowField(
        const Inpar::FLUID::InitialField initfield, const int startfuncno) override;

    /// Implement Adapter::Fluid
    Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres) override;

    /// Implement Adapter::Fluid
    Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres) override;

    /// Reset state vectors
    void Reset(bool completeReset = false, int numsteps = 1, int iter = -1) override;

    /*!
    \brief calculate error between a analytical solution and the
           numerical solution of a test problems

    */
    virtual Teuchos::RCP<std::vector<double>> evaluate_error_compared_to_analytical_sol();

    /*!
    \brief evaluate divergence of velocity field

    */
    virtual Teuchos::RCP<double> EvaluateDivU();

    /*!
    \brief calculate adaptive time step with the CFL number

    */
    virtual double evaluate_dt_via_cfl_if_applicable();

    /*!
    \brief read restart data

    */
    void read_restart(int step) override;

    /*!
    \brief get restart data in case of turbulent inflow computation

    */
    void SetRestart(const int step, const double time, Teuchos::RCP<const Epetra_Vector> readvelnp,
        Teuchos::RCP<const Epetra_Vector> readveln, Teuchos::RCP<const Epetra_Vector> readvelnm,
        Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn) override;

    //! @name access methods for composite algorithms
    /// monolithic FSI needs to access the linear fluid problem

    Teuchos::RCP<const Epetra_Vector> initial_guess() override { return incvel_; }

    /// return implemented residual (is not an actual force in Newton [N])
    virtual Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }

    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS() override { return Residual(); }

    /// Return true residual, ie the actual force in Newton [N]
    Teuchos::RCP<const Epetra_Vector> TrueResidual() override { return trueresidual_; }

    Teuchos::RCP<const Epetra_Vector> Velnp() override { return velnp_; }
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override { return velnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf() override { return velaf_; }
    virtual Teuchos::RCP<const Epetra_Vector> Velam() { return velam_; }
    Teuchos::RCP<const Epetra_Vector> Veln() override { return veln_; }
    Teuchos::RCP<const Epetra_Vector> Velnm() override { return velnm_; }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp() { return accnp_; }
    Teuchos::RCP<const Epetra_Vector> Accnp() override { return accnp_; }
    Teuchos::RCP<const Epetra_Vector> Accn() override { return accn_; }
    Teuchos::RCP<const Epetra_Vector> Accnm() override { return accnm_; }
    Teuchos::RCP<const Epetra_Vector> Accam() override { return accam_; }
    Teuchos::RCP<const Epetra_Vector> Scaaf() override { return scaaf_; }
    Teuchos::RCP<const Epetra_Vector> Scaam() override { return scaam_; }
    Teuchos::RCP<const Epetra_Vector> Hist() override { return hist_; }
    Teuchos::RCP<const Epetra_Vector> GridVel() override { return gridv_; }
    Teuchos::RCP<const Epetra_Vector> GridVeln() override { return gridvn_; }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessGridVel() { return gridv_; }
    Teuchos::RCP<const Epetra_Vector> FsVel() override
    {
      // get fine-scale part of velocity at time n+alpha_F or n+1
      if (Sep_ != Teuchos::null)
      {
        Sep_Multiply();
      }

      // set fine-scale velocity for parallel nigthly tests
      // separation matrix depends on the number of proc here
      if (turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales and
          (Core::UTILS::IntegralValue<int>(
              params_->sublist("MULTIFRACTAL SUBGRID SCALES"), "SET_FINE_SCALE_VEL")))
        fsvelaf_->PutScalar(0.01);

      return fsvelaf_;
    }

    /// access to Dirichlet maps
    Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() override
    {
      return dbcmaps_;
    }

    /// Expand the Dirichlet DOF set
    ///
    /// The method expands the DOF set (map) which contains the DOFs
    /// subjected to Dirichlet boundary conditions. For instance, the method is
    /// called by the staggered FSI in which the velocities on the FSI
    /// interface are prescribed by the other fields.
    void add_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoadd) override;

    /// Contract the Dirichlet DOF set
    ///
    /// !!! Be careful using this! You might delete dirichlet values set in the input file !!!
    /// !!! So make sure you are only touching the desired dofs.                           !!!
    ///
    /// The method contracts the DOF set (map) which contains the DOFs
    /// subjected to Dirichlet boundary conditions. This method is
    /// called solely by immersed FSI to remove the Dirichlet values from
    /// the previous solution step before a new set is prescribed.
    void remove_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoremove) override;

    /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    virtual Teuchos::RCP<const Epetra_Vector> Dirichlet();

    /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
    ///
    /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
    /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
    /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
    virtual Teuchos::RCP<const Epetra_Vector> InvDirichlet();

    //! Return locsys manager
    virtual Teuchos::RCP<Core::Conditions::LocsysManager> LocsysManager() { return locsysman_; }

    //! Return wss manager
    virtual Teuchos::RCP<FLD::UTILS::StressManager> StressManager() { return stressmanager_; }

    //! Return impedance BC
    virtual Teuchos::RCP<FLD::UTILS::FluidImpedanceWrapper> ImpedanceBC_() { return impedancebc_; }

    //! Evaluate Dirichlet and Neumann boundary conditions
    virtual void set_dirichlet_neumann_bc();

    //! Apply Dirichlet boundary conditions on provided state vectors
    virtual void apply_dirichlet_bc(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector,    //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> systemvectord,   //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> systemvectordd,  //!< (may be Teuchos::null)
        bool recreatemap                             //!< recreate mapextractor/toggle-vector
                                                     //!< which stores the DOF IDs subjected
                                                     //!< to Dirichlet BCs
                                                     //!< This needs to be true if the bounded DOFs
                                                     //!< have been changed.
    );

    Teuchos::RCP<const Epetra_Vector> Dispnp() override { return dispnp_; }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() { return dispnp_; }

    //! Create mesh displacement at time level t_{n+1}
    virtual Teuchos::RCP<Epetra_Vector> CreateDispnp()
    {
      const Epetra_Map* aledofrowmap = discret_->dof_row_map(ndsale_);
      dispnp_ = Core::LinAlg::CreateVector(*aledofrowmap, true);
      return dispnp_;
    }

    Teuchos::RCP<const Epetra_Vector> Dispn() override { return dispn_; }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispn() { return dispn_; }

    //! Create mesh displacement at time level t_{n}
    virtual Teuchos::RCP<Epetra_Vector> CreateDispn()
    {
      const Epetra_Map* aledofrowmap = discret_->dof_row_map(ndsale_);
      dispn_ = Core::LinAlg::CreateVector(*aledofrowmap, true);
      return dispn_;
    }
    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix() override
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_);
    }
    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemSparseMatrix() override
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_)->Merge();
    }
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_);
    }
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> ShapeDerivatives() override
    {
      return shapederivatives_;
    }

    virtual Teuchos::RCP<Core::LinAlg::MapExtractor> VelPresSplitter() { return velpressplitter_; };
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() override;
    Teuchos::RCP<const Epetra_Map> PressureRowMap() override;
    //  virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);
    //  double TimeScaling() const;

    /// Use residual_scaling() to convert the implemented fluid residual to an actual force with
    /// unit Newton [N]
    /*! In order to avoid division by time step size \f$\Delta t\f$
     *  the fluid balance of linear momentum is implemented in a way
     *  that the residual does not have the unit Newton [N].
     *  By multiplication with residual_scaling() the residual is
     *  converted to the true residual in unit Newton [N], ie a real force.
     *
     *  \sa trueresidual_
     *  \sa TrueResidual()
     */
    double residual_scaling() const override = 0;

    /*!
    \brief return scheme-specific time integration parameter

    */
    double TimIntParam() const override = 0;

    /*!
    \brief compute values at intermediate time steps for gen.-alpha
           for given vectors and store result in given vectors
           Helper method which can be called from outside fluid (e.g. for coupled problems)

    */
    virtual void gen_alpha_intermediate_values(
        Teuchos::RCP<Epetra_Vector>& vecnp, Teuchos::RCP<Epetra_Vector>& vecn)
    {
      return;
    }

    /// update velocity increment after Newton step
    void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel) override;

    //  int Itemax() const { return params_->get<int>("max nonlin iter steps"); }
    void SetItemax(int itemax) override { params_->set<int>("max nonlin iter steps", itemax); }


    /*!
    \brief set scalar fields within outer iteration loop

    */
    void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<Core::FE::Discretization> scatradis, int dofset) override;

    /*!
    \brief set scalar fields

    */
    void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp, const double thermpressnp,
        Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<Core::FE::Discretization> scatradis, const int whichscalar = -1) override;

    /*!
    \brief set velocity field obtained by separate computation

    */
    void set_velocity_field(Teuchos::RCP<const Epetra_Vector> setvelnp) override
    {
      velnp_->Update(1.0, *setvelnp, 0.0);
      return;
    }

    /// provide access to turbulence statistics manager
    Teuchos::RCP<FLD::TurbulenceStatisticManager> turbulence_statistic_manager() override;
    /// provide access to the box filter for dynamic Smagorinsky model
    Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() override;
    /// provide access to the box filter for Vreman model
    Teuchos::RCP<FLD::Vreman> Vreman() override;

    /// introduce surface split extractor object
    /*!
      This method must (and will) be called during setup with a properly
      initialized extractor object if we are on an ale mesh.
     */
    virtual void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter)
    {
      surfacesplitter_ = surfacesplitter;
    }

    /// determine grid velocity
    virtual void UpdateGridv();

    /// prepare AVM3-based scale separation
    virtual void av_m3_preparation();

    /// AVM3-based scale separation
    virtual void av_m3_separation();

    /// compute flow rate
    virtual void ComputeFlowRates() const;

    /// integrate shape functions at nodes marked by condition
    /*!
      Needed for Mortar coupling at the FSI interface
     */
    virtual Teuchos::RCP<Epetra_Vector> integrate_interface_shape(std::string condname);

    /// switch fluid field to block matrix
    virtual void use_block_matrix(
        Teuchos::RCP<std::set<int>> condelements,  ///< conditioned elements of fluid
        const Core::LinAlg::MultiMapExtractor&
            domainmaps,  ///< domain maps for split of fluid matrix
        const Core::LinAlg::MultiMapExtractor& rangemaps,  ///< range maps for split of fluid matrix
        bool splitmatrix = true                            ///< flag for split of matrices
    );

    /// switch fluid field to block matrix (choose maps for shape derivatives separately)
    virtual void use_block_matrix(
        Teuchos::RCP<std::set<int>> condelements,  ///< conditioned elements of fluid
        const Core::LinAlg::MultiMapExtractor&
            domainmaps,  ///< domain maps for split of fluid matrix
        const Core::LinAlg::MultiMapExtractor& rangemaps,  ///< range maps for split of fluid matrix
        Teuchos::RCP<std::set<int>> condelements_shape,    ///< conditioned elements
        const Core::LinAlg::MultiMapExtractor&
            domainmaps_shape,  ///< domain maps for split of shape deriv. matrix
        const Core::LinAlg::MultiMapExtractor&
            rangemaps_shape,     ///< domain maps for split of shape deriv. matrix
        bool splitmatrix = true  ///< flag for split of matrices
    );

    /// linear solve with prescribed dirichlet conditions and without history
    /*!
      This is the linear solve as needed for steepest descent FSI.
     */
    virtual void linear_relaxation_solve(Teuchos::RCP<Epetra_Vector> relax);

    //@}

    //! @name methods for turbulence models

    virtual void apply_scale_separation_for_les();

    virtual void OutputofFilteredVel(
        Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec) = 0;

    virtual void print_turbulence_model();
    //@}

    /// set the initial porosity field
    void set_initial_porosity_field(
        const Inpar::PoroElast::InitialField,  ///< type of initial field
                                               // const int, ///< type of
                                               // initial field
        const int startfuncno                  ///< number of spatial function
        ) override
    {
      FOUR_C_THROW("not implemented in base class");
    }

    virtual void update_iter_incrementally(
        Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
    );

    //! @name methods for fsi
    /// Extrapolation of vectors from mid-point to end-point t_{n+1}
    virtual Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
        Teuchos::RCP<Epetra_Vector> vecn,  ///< vector at time level t_n
        Teuchos::RCP<Epetra_Vector> vecm   ///< vector at time level of equilibrium
    );
    //@}

    /// apply external forces to the fluid
    void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext) override;

    /// create field test
    Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() override;

    Teuchos::RCP<const Epetra_Vector> ConvectiveVel() override;

    /*! \brief Calculate a integrated divergence operator in vector form
     *
     *   The vector valued operator \f$B\f$ is constructed such that
     *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
     */
    virtual Teuchos::RCP<Epetra_Vector> CalcDivOp();

    //! @name Biofilm methods
    //@{

    // set fluid displacement vector due to biofilm growth
    void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp) override;
    //@}

    /*!
    \brief evaluate and update problem-specific boundary conditions

    */
    virtual void do_problem_specific_boundary_conditions() { return; }

    ///< Print stabilization details to screen
    virtual void print_stabilization_details() const;

    ///< Add contribution to external load vector ( add possibly pre-existing external_loads_);
    void add_contribution_to_external_loads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector) override;

    ///< Update slave dofs for multifield simulations with fluid mesh tying
    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& f);

    /** \brief Here additional contributions to the system matrix can be set
     *
     * To enforce weak dirichlet conditions as they arise from meshtying for example, such
     * contributions can be set here and will be assembled into the system matrix
     *
     * \param[in] matrix (size fluid_dof x fluid_dof) linear matrix containing entries that need
     * to be assembled into the overall fluid system matrix
     */
    virtual void set_coupling_contributions(
        Teuchos::RCP<const Core::LinAlg::SparseOperator> matrix);

    void ResetExternalForces();

    Teuchos::RCP<const FLD::Meshtying> GetMeshtying() { return meshtying_; }

   protected:
    // don't want = operator and cctor
    // FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
    FluidImplicitTimeInt(const FluidImplicitTimeInt& old);

    /*!
    \brief timeloop break criterion
     */
    virtual bool not_finished() { return step_ < stepmax_ and time_ < maxtime_; }

    /*!
    \brief  increment time and step value

    */
    void increment_time_and_step() override
    {
      step_ += 1;
      time_ += dta_;
    }

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    virtual void assemble_mat_and_rhs();

    /*!
    \brief call elements to calculate system matrix/rhs and assemble, called from
    assemble_mat_and_rhs

    */
    virtual void evaluate_mat_and_rhs(Teuchos::ParameterList& eleparams);

    /*!
    \brief calculate intermediate solution

    */
    void calc_intermediate_solution() override;

    /*!
    \brief apply Dirichlet boundary conditions to system of equations

    */
    virtual void apply_dirichlet_to_system();

    /*!
    \brief apply weak or mixed hybrid Dirichlet boundary conditions to system of equations

    */
    virtual void apply_nonlinear_boundary_conditions();

    /*!
    \brief update acceleration for generalized-alpha time integration

    */
    virtual void gen_alpha_update_acceleration() { return; }

    /*!
    \brief compute values at intermediate time steps for gen.-alpha

    */
    virtual void gen_alpha_intermediate_values() { return; }

    //! Predict velocities which satisfy exactly the Dirichlet BCs
    //! and the linearised system at the previously converged state.
    //!
    //! This is an implicit predictor, i.e. it calls the solver once.
    virtual void predict_tang_vel_consist_acc();

    /*!
    \brief update surface tension (free surface flow only)

    */
    virtual void free_surface_flow_surface_tension_update();

    /*!
    \brief Update of an Ale field based on the fluid state

    */
    virtual void ale_update(std::string condName);

    /*!
    \brief For a given node, obtain local indices of dofs in a vector (like e.g. velnp)

    */
    void get_dofs_vector_local_indicesfor_node(int nodeGid, Teuchos::RCP<Epetra_Vector> vec,
        bool withPressure, std::vector<int>* dofsLocalInd);

    /*!
    \brief add mat and rhs of edge-based stabilization

    */
    virtual void assemble_edge_based_matand_rhs();

    /*!
    \brief Setup meshtying

    */
    virtual void setup_meshtying();

    /*!
    \brief velocity required for evaluation of related quantites required on element level

    */
    virtual Teuchos::RCP<const Epetra_Vector> evaluation_vel() = 0;

    /*!
      \brief add problem dependent vectors

     */
    virtual void add_problem_dependent_vectors() { return; };

    /*!
    \brief Initialize forcing

    */
    virtual void init_forcing();

    /*!
    \brief calculate lift&drag forces and angular momenta

    */
    void LiftDrag() const override;

    /**
     * \brief Here, the coupling contributions collected in the linear matrix couplingcontributions_
     * will be added to the system matrix
     */
    virtual void assemble_coupling_contributions();

    //! @name general algorithm parameters

    //! do we move the fluid mesh and calculate the fluid on this moving mesh?
    bool alefluid_;
    //! do we have a turbulence model?
    enum Inpar::FLUID::TurbModelAction turbmodel_;

    //@}

    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    int numstasteps_;  ///< number of steps for starting algorithm
    //@}


    /// gas constant (only for low-Mach-number flow)
    double gasconstant_;

    //! use (or not) linearisation of reactive terms on the element
    Inpar::FLUID::LinearisationAction newton_;

    //! kind of predictor used in nonlinear iteration
    std::string predictor_;

    //! @name restart variables
    int writestresses_;
    int write_wall_shear_stresses_;
    int write_eledata_everystep_;
    //@}

    int write_nodedata_first_step_;

    //! @name time step sizes
    double dtp_;  ///< time step size of previous time step
    //@}

    //! @name time-integration-scheme factors
    /// declaration required here in base class
    double theta_;

    //@}

    //! @name parameters for sampling/dumping period
    int samstart_;
    int samstop_;
    int dumperiod_;
    //@}

    std::string statistics_outfilename_;

    //! @name cfl number for adaptive time step
    Inpar::FLUID::AdaptiveTimeStepEstimator cfl_estimator_;  ///< type of adaptive estimator
    double cfl_;                                             ///< cfl number
    //@}

    //! @name norms for convergence check
    double incvelnorm_L2_;
    double incprenorm_L2_;
    double velnorm_L2_;
    double prenorm_L2_;
    double vresnorm_;
    double presnorm_;
    //@}

    //! flag to skip calculation of residual after solution has converged
    bool inconsistent_;

    //! flag to reconstruct second derivative for fluid residual
    bool reconstructder_;

    /// flag for special turbulent flow
    std::string special_flow_;

    /// flag for potential nonlinear boundary conditions
    bool nonlinearbc_;

    /// form of convective term
    std::string convform_;

    /// fine-scale subgrid-viscosity flag
    Inpar::FLUID::FineSubgridVisc fssgv_;

    /// cpu-time measures
    double dtele_;
    double dtfilter_;
    double dtsolve_;

    /// (standard) system matrix
    Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat_;

    /// linearization with respect to mesh motion
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> shapederivatives_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    /// the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;

    /// the vector containing external loads
    Teuchos::RCP<Epetra_Vector> external_loads_;

    /// the vector containing volume force externally computed
    Teuchos::RCP<Epetra_Vector> forcing_;

    /// the vector containing potential Neumann-type outflow terms
    //  Teuchos::RCP<Epetra_Vector>    outflow_;

    /// a vector containing the integrated traction in boundary normal direction for slip boundary
    /// conditions (Unit: Newton [N])
    Teuchos::RCP<Epetra_Vector> slip_bc_normal_tractions_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector> residual_;

    /// true (rescaled) residual vector without zeros at dirichlet positions (Unit: Newton [N])
    Teuchos::RCP<Epetra_Vector> trueresidual_;

    /// Nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector> incvel_;

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n) and
    //! n-1
    //@{
    Teuchos::RCP<Epetra_Vector> accnp_;  ///< acceleration at time \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> accn_;   ///< acceleration at time \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> accam_;  ///< acceleration at time \f$t^{n+\alpha_M}\f$
    Teuchos::RCP<Epetra_Vector> accnm_;  ///< acceleration at time \f$t^{n-1}\f$
    //@}

    //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F (and n+alpha_M for
    //! weakly_compressible)
    //@{
    Teuchos::RCP<Epetra_Vector> velnp_;  ///< velocity at time \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;   ///< velocity at time \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> velaf_;  ///< velocity at time \f$t^{n+\alpha_F}\f$
    Teuchos::RCP<Epetra_Vector> velam_;  ///< velocity at time \f$t^{n+\alpha_M}\f$
    Teuchos::RCP<Epetra_Vector> velnm_;  ///< velocity at time \f$t^{n-1}\f$
    //@}

    //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
    Teuchos::RCP<Epetra_Vector> scaaf_;
    Teuchos::RCP<Epetra_Vector> scaam_;
    //@}

    //! @name displacements at time n+1, n and n-1
    //@{
    Teuchos::RCP<Epetra_Vector> dispnp_;  ///< displacement at time \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> dispn_;   ///< displacement at time \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> dispnm_;  ///< displacement at time \f$t^{n-1}\f$
    //@}

    //! @name flow rate and volume at time n+1 (i+1), n+1 (i), n and n-1 for flow-dependent pressure
    //! boundary conditions
    std::vector<double> flowratenp_;
    std::vector<double> flowratenpi_;
    std::vector<double> flowraten_;
    std::vector<double> flowratenm_;

    std::vector<double> flowvolumenp_;
    std::vector<double> flowvolumenpi_;
    std::vector<double> flowvolumen_;
    std::vector<double> flowvolumenm_;

    //@}

    /// only necessary for AVM3: scale-separation matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Sep_;

    /// only necessary for AVM3: fine-scale solution vector
    Teuchos::RCP<Epetra_Vector> fsvelaf_;

    /// only necessary for LES models including filtered quantities: filter type
    enum Inpar::FLUID::ScaleSeparation scale_sep_;

    /// fine-scale scalar: only necessary for multifractal subgrid-scale modeling in loma
    Teuchos::RCP<Epetra_Vector> fsscaaf_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector> gridv_;

    /// grid velocity at time step n (set from the adapter!)
    Teuchos::RCP<Epetra_Vector> gridvn_;

    /// histvector --- a linear combination of velnm, veln (BDF)
    ///                or veln, accn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector> hist_;


    //! manager for turbulence statistics
    Teuchos::RCP<FLD::TurbulenceStatisticManager> statisticsmanager_;

    //! forcing for homogeneous isotropic turbulence
    Teuchos::RCP<FLD::ForcingInterface> forcing_interface_;

    //! @name Dynamic Smagorinsky model: methods and variables
    //        -------------------------

    //! one instance of the filter object
    Teuchos::RCP<FLD::DynSmagFilter> DynSmag_;
    //! one instance of the filter object
    Teuchos::RCP<FLD::Vreman> Vrem_;
    Teuchos::RCP<FLD::Boxfilter> Boxf_;

    //@}

    //! Extractor to split velnp_ into velocities and pressure DOFs
    //!
    //! velocities  = OtherVector
    //! pressure    = CondVector
    Teuchos::RCP<Core::LinAlg::MapExtractor> velpressplitter_;

    /// row dof map extractor
    const UTILS::MapExtractor* surfacesplitter_;

    /// a manager doing the transfer of boundary data for
    /// turbulent inflow profiles from a separate (periodic) domain
    Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

    /// @name special relaxation state

    bool inrelaxation_;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> dirichletlines_;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> meshmatrix_;

    /// coupling of fluid-fluid at an internal interface
    Teuchos::RCP<FLD::Meshtying> meshtying_;

    /// coupling of fluid-fluid at an internal interface
    Teuchos::RCP<FLD::XWall> xwall_;

    /// flag for mesh-tying
    enum Inpar::FLUID::MeshTying msht_;

    /// face discretization (only initialized for edge-based stabilization)
    Teuchos::RCP<Core::FE::DiscretizationFaces> facediscret_;

    //@}

    // possible inf-norm scaling of linear system / fluid matrix
    Teuchos::RCP<FLD::UTILS::FluidInfNormScaling> fluid_infnormscaling_;

    //! @name Biofilm specific stuff
    //@{
    Teuchos::RCP<Epetra_Vector> fldgrdisp_;
    //@}

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<Core::Conditions::LocsysManager> locsysman_;

    /// windkessel (outflow) boundaries
    Teuchos::RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<FLD::UTILS::StressManager> stressmanager_;

    /// flag for windkessel outflow condition
    bool isimpedancebc_;

    /// flag for windkessel outflow condition
    bool off_proc_assembly_;

    /// number of dofset for ALE quantities (0 by default and 2 in case of HDG discretizations)
    unsigned int ndsale_;


    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void set_element_time_parameter() = 0;

   private:
    //! @name Special method for turbulent variable-density flow at low Mach number with
    //! multifractal subgrid-scale modeling
    /*!

    \brief adaptation of CsgsD to CsgsB
           Since CsgsB depends on the resolution if the near-wall limit is included,
           CsgsD is adapted accordingly by using the mean value of the near-wall  correction.

    */
    virtual void recompute_mean_csgs_b();


    /*! \brief Prepares the locsys manager by calculating the node normals
     *
     */
    void setup_locsys_dirichlet_bc(const double time);

    /// prepares and evalutes egde-based internal face integrals
    void evaluate_fluid_edge_based(Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::ParameterList edgebasedparams);

    /*! \brief Compute kinetic energy and write it to file
     *
     *  Kinetic energy of the system is calculated as \f$E_{kin} = \frac{1}{2}u^TMu\f$
     *  with the velocity vector \f$u\f$ and the mass matrix \f$M\f$. Then, it is
     *  written to an output file.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void write_output_kinetic_energy();

    ///< Evaluate mass matrix
    virtual void evaluate_mass_matrix();

    /// mass matrix (not involved in standard evaluate() since it is invluded in #sysmat_)
    Teuchos::RCP<Core::LinAlg::SparseOperator> massmat_;

    /// output stream for energy-file
    Teuchos::RCP<std::ofstream> logenergy_;

    /** \brief This matrix can be used to hold contributions to the system matrix like such that
     * arise from meshtying methods or in general weak Dirichlet conditions
     *
     */
    Teuchos::RCP<const Core::LinAlg::SparseOperator> couplingcontributions_;
    double meshtyingnorm_;

  };  // class FluidImplicitTimeInt

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
