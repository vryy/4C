/*----------------------------------------------------------------------*/
/*! \file


\brief H-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches

\level 3


----------------------------------------------------------------------*/


#ifndef FOUR_C_FS3I_AC_FSI_HPP
#define FOUR_C_FS3I_AC_FSI_HPP

#include "4C_config.hpp"

#include "4C_fs3i_partitioned_1wc.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class MapExtractor;
}

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace FS3I
{
  class MeanManager;

  /*!
   * What does the problem type Atherosclerosis_Fluid_Structure_Interaction?
   * Short answer: cool stuff!
   * And here is the long answer:
   * Its doing a multiscale (in time) approach with an full FS3I simulation at the small time scale
   * (seconds) and a scalar transport simulation at the larger time scale (days). The solving
   * strategy is as follows: We start with the full small time scale FS3i simulation (including
   * fluid Windkessel and WSS permeability). After each FSI cycle we check if the FSI problem is
   * periodic, by looking if the Windkessel does produce peridodic results. Afterwards we continue
   * the small time scale but do not solve the fsi subproblem anymore. Instead we peridically repeat
   * it by calling suitable restarts. When the fluid scatra subproblem gets periodic at the FS3I
   * interface we stop the small time scale and switch to the large time scale Now we higher dt_ and
   * only solve the structure scatra problem. We thereby use the WSS and interface concentrations of
   * the small time scale. Each time when there as been 'created' enough growth inducing mass we do
   * a growth update. If we have finally grew to much, we go back to the small time scale. And so
   * on, and so on,...
   */
  class ACFSI : public PartFS3I
  {
   public:
    /// constructor
    ACFSI(const Epetra_Comm& comm);

    /// initialize this class
    void init() override;

    /// setup this class
    void setup() override;

    /// Read restart
    void read_restart() override;

    /// timeloop
    void timeloop() override;

    /// timeloop for small time scales
    void small_time_scale_loop();

    /// flag whether small time scale time loop should be finished
    bool small_time_scale_loop_not_finished();

    /// Prepare small time scale time step
    void small_time_scale_prepare_time_step();

    /// Prepare time step
    void prepare_time_step() override
    {
      FOUR_C_THROW(
          "This function is not implemented! Use small_time_scale_prepare_time_step() or "
          "large_time_scale_prepare_time_step() instead!");
    };

    /// outer_loop
    void small_time_scale_outer_loop();

    /// outer_loop for sequentially staggered FS3I scheme
    void small_time_scale_outer_loop_sequ_stagg();

    /// outer_loop for iterative staggered FS3I scheme
    void small_time_scale_outer_loop_iter_stagg();

    /// Do a single fsi step (including all subcycles)
    void do_fsi_step();

    void is_small_time_scale_periodic();

    /// Decide if fsi problem is already periodic
    void is_fsi_periodic();

    /// provide wall shear stresses from FS3I subproblem for scatra subproblem
    void set_wall_shear_stresses() const override;

    /// Decide if fluid scatra problem is periodic
    void is_scatra_periodic();

    /// Do a standard fsi step
    void do_fsi_step_standard();

    /// Do a fsi step with subcycling
    void do_fsi_step_subcycled(const int subcyclingsteps);

    /// Get fsi solution from one period before
    void do_fsi_step_periodic();

    /// Get step number of on cycle ago
    double get_step_of_one_period_ago_and_prepare_reading(const int actstep, const double time);

    /// Get step number of the beginning of this cycle
    double get_step_of_beginn_of_this_period_and_prepare_reading(
        const int actstep, const double acttime, const double dt);

    /// Get filename in which the equivalent step of the last period is written
    std::string get_file_name(const int step);

    /// Set time and step in FSI and all subfields
    void set_time_and_step_in_fsi(const double time, const int step);

    /// Do a single scatra step
    void small_time_scale_do_scatra_step();

    /// Update and output the small time scale
    void small_time_scale_update_and_output();

    /// Write FSI output
    void fsi_output();

    /// check convergence of scatra fields
    bool scatra_convergence_check(const int itnum) override;

    /// Convergence check for iterative staggered FS3I scheme
    bool part_fs3i_convergence_ckeck(const int itnum);

    //--------------------------------------------------------------------
    /// @name  control routines for the large time scale
    //--------------------------------------------------------------------

    /// timeloop for large time scales
    void large_time_scale_loop();

    /// Prepare the large time scale loop
    void prepare_large_time_scale_loop();

    /// Set mean wall shear stresses in scatra fields
    void set_mean_wall_shear_stresses() const;

    /// Set mean concentration of the fluid scatra field
    void set_mean_fluid_scatra_concentration();

    /// Set zero velocity field in scatra fields
    void set_zero_velocity_field();

    /// Evaluate surface permeability condition for structure scatra field
    void evaluateith_scatra_surface_permeability(const int i  ///< id of scalar to evaluate
    );

    /// Finish the large time scale loop
    void finish_large_time_scale_loop();

    /// flag whether large time scale time loop should be finished
    bool large_time_scale_loop_not_finished();

    /// Prepare large time scale time step
    void large_time_scale_prepare_time_step();

    /// outer_loop for sequentially staggered FS3I scheme
    void large_time_scale_outer_loop();

    /// Do a large time scale structe scatra step
    void do_struct_scatra_step();

    /// evaluate, solver and iteratively update structure scalar problem
    void struct_scatra_evaluate_solve_iter_update();

    /// check convergence of structure scatra field
    bool struct_scatra_convergence_check(const int itnum  ///< current iteration number
    );

    /// Do the structure scatra displacments need to update
    bool does_growth_needs_update();

    /// update the structure scatra displacments due to growth
    void large_time_scale_do_growth_update();

    /// outer_loop for large time scale iterative staggered FS3I scheme
    void large_time_scale_outer_loop_iter_stagg();

    /// set mean FSI values in scatra fields (only to be used in large time scale!!)
    void large_time_scale_set_fsi_solution();

    /// Update and output the large time scale
    void large_time_scale_update_and_output();
    //@}

    /// Build map extractor which extracts the j-th dof
    std::vector<Teuchos::RCP<Core::LinAlg::MapExtractor>> build_map_extractor();

    /// optional safety check for times and dt's of all fields
    void check_if_times_and_steps_and_dts_match();

    /// Compare if two doubles are relatively zero
    bool is_realtive_equal_to(const double A,  ///< first value
        const double B,                        ///< second value
        const double Ref = 1.0                 ///< reference value
    );

    /// Compare if A mod B is relatively equal to zero
    bool modulo_is_realtive_zero(const double value,  ///< value to mod
        const double modulo,                          ///< mod value
        const double Ref = 1.0                        ///< reference value
    );

   private:
    /// structure increment vector
    Teuchos::RCP<Epetra_Vector> structureincrement_;
    /// fluid increment vector
    Teuchos::RCP<Epetra_Vector> fluidincrement_;
    /// ale increment vector
    Teuchos::RCP<Epetra_Vector> aleincrement_;
    /// mean fluid phinp vector of the last period
    Teuchos::RCP<Epetra_Vector> fluidphinp_lp_;
    /// structurephinp vector at the beginning of the large time scale loop
    Teuchos::RCP<Epetra_Vector> structurephinp_blts_;
    /// growth update counter
    int growth_updates_counter_;

    /// mean WSS vector of the last period
    Teuchos::RCP<Epetra_Vector> wall_shear_stress_lp_;

    /// time of one fsi period, e.g. time of a heart cycle
    const double fsiperiod_;

    /// time step for the large time scale problem
    const double dt_large_;

    /// flag iff fsi subproblem is periodic
    bool fsiisperiodic_;

    /// flag iff fluid scatra subproblem is periodic
    bool scatraisperiodic_;

    /// flag iff fluid scatra subproblem is periodic
    bool fsineedsupdate_;

    /// Extract the j-th out of numscal_ dof
    std::vector<Teuchos::RCP<Core::LinAlg::MapExtractor>> extractjthstructscalar_;

    /// pointer to mean manager object
    Teuchos::RCP<FS3I::MeanManager> meanmanager_;
  };

  class MeanManager
  {
   public:
    /// constructor
    MeanManager(const Epetra_Map& wssmap, const Epetra_Map& phimap, const Epetra_Map& pressuremap);

    /// destructor
    virtual ~MeanManager() = default;

    /// add value into the mean manager
    void add_value(
        const std::string type, const Teuchos::RCP<const Epetra_Vector> value, const double dt);

    /// reset mean manager
    void reset();

    /// get some mean value
    Teuchos::RCP<const Epetra_Vector> get_mean_value(const std::string type) const;

    /// Write restart of mean manager
    void write_restart(Teuchos::RCP<Core::IO::DiscretizationWriter> fluidwriter) const;

    /// Read restart of mean manager
    void read_restart(Core::IO::DiscretizationReader& fluidreader);

   private:
    /// weighted sum of all prior wall shear stresses
    Teuchos::RCP<Epetra_Vector> sum_wss_;
    /// weighted sum of all prior concentrations
    Teuchos::RCP<Epetra_Vector> sum_phi_;
    /// weighted sum of all prior pressures
    Teuchos::RCP<Epetra_Vector> sum_pres_;

    double sum_dt_wss_;
    double sum_dt_phi_;
    double sum_dt_pres_;
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
