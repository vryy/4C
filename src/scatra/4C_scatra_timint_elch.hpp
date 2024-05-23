/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_ELCH_HPP
#define FOUR_C_SCATRA_TIMINT_ELCH_HPP

#include "4C_config.hpp"

#include "4C_inpar_elch.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace FLD
{
  class Meshtying;
}

namespace ADAPTER
{
  class Coupling;
  class ScaTraBaseAlgorithm;
}  // namespace ADAPTER

namespace SCATRA
{
  class CCCVCondition;

  class ScaTraTimIntElch : public virtual ScaTraTimIntImpl
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScaTraTimIntElch(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! initialize algorithm
    void Init() override;

    //! initialize algorithm
    void Setup() override;

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/

    //! Setup splitter for concentration and potential dofs
    void SetupSplitter() override;

    //! additional, to standard partitioning in scatra, the global system matrix in elch can be
    //! partitioned into concentration and potential dofs
    void BuildBlockMaps(
        const std::vector<Teuchos::RCP<CORE::Conditions::Condition>>& partitioningconditions,
        std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const override;

    void build_block_null_spaces(
        Teuchos::RCP<CORE::LINALG::Solver> solver, int init_block_number) const override;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! Set elch-specific parameters
    void set_element_specific_sca_tra_parameters(Teuchos::ParameterList& eleparams) const override;

    //! contains the nonlinear iteration loop
    void NonlinearSolve() override;

    //! calculate error compared to analytical solution
    void evaluate_error_compared_to_analytical_sol() override;

    Teuchos::RCP<CORE::UTILS::ResultTest> create_sca_tra_field_test() override;

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! setup natural convection
    void SetupNatConv() override;

    /*--- calculate and update -----------------------------------------------*/

    [[nodiscard]] bool NotFinished() const override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

    /*--- query and output ---------------------------------------------------*/

    void check_and_write_output_and_restart() override;

    //! problem-specific outputs
    void output_problem_specific() override;

    //! read problem-specific restart data
    void read_restart_problem_specific(int step, IO::DiscretizationReader& reader) override;

    //! output electrode domain status information to screen and file
    void output_electrode_info_domain();

    //! output electrode boundary status information to screen and file
    void output_electrode_info_boundary();

    //! \brief evaluate status information on single line or surface electrode
    //!
    //! \param condid      ID of condition to be evaluated
    //! \param condstring  name of condition to be evaluated
    //! \return            evaluated scalars
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> evaluate_single_electrode_info(
        int condid, const std::string& condstring);

    //! \brief evaluate status information on single point electrode
    //!
    //! \param condition   condition to be evaluated
    //! \return            evaluated scalars
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> evaluate_single_electrode_info_point(
        Teuchos::RCP<CORE::Conditions::Condition> condition);

    //! \brief post-process status information on single electrode
    //!
    //! \param scalars        scalar quantities associated with electrode status information (in)
    //! \param                id electrode ID (in)
    //! \param print          flag for output to screen and file (in)
    //! \param currentsum     net current involving all conditions (out)
    //! \param currtangent    tangent of current w.r.t. electrode potential (out)
    //! \param currresidual   negative residual of current equation (out)
    //! \param electrodeint   physical dimensions of the electrode region (out)
    //! \param electrodepot   electrode potential on electrode side (out)
    //! \param meanoverpot    mean overpotential (out)
    void post_process_single_electrode_info(CORE::LINALG::SerialDenseVector& scalars, int id,
        bool print, double& currentsum, double& currtangent, double& currresidual,
        double& electrodeint, double& electrodepot, double& meanoverpot);

    //! output electrode interior status information to screen and files
    void output_electrode_info_interior();

    //! output cell voltage to screen and file
    void OutputCellVoltage();

    void WriteRestart() const override;

    //! output type of closing equation for electric potential
    INPAR::ELCH::EquPot EquPot() const { return equpot_; }

    //! return constant F/RT
    double FRT() const { return fr_ / temperature_; }

    //! current temperature is determined and returned
    double get_current_temperature() const;

    //! return elch parameter list
    Teuchos::RCP<const Teuchos::ParameterList> ElchParameterList() const { return elchparams_; }

    //! return states of charge of resolved electrodes
    const std::map<int, double>& ElectrodeSOC() const { return electrodesoc_; };

    //! return C rates with respect to resolved electrodes
    const std::map<int, double>& ElectrodeCRates() const { return electrodecrates_; };

    //! return mean reactant concentrations at electrode boundaries
    const std::map<int, double>& ElectrodeConc() const { return electrodeconc_; };

    //! return mean electric overpotentials at electrode boundaries
    const std::map<int, double>& ElectrodeEta() const { return electrodeeta_; };

    //! return total electric currents at electrode boundaries
    const std::map<int, double>& ElectrodeCurr() const { return electrodecurr_; };

    //! return cell voltage
    const double& CellVoltage() const { return cellvoltage_; };

    //! return map extractor for macro scale in multi-scale simulations
    const Teuchos::RCP<const CORE::LINALG::MultiMapExtractor>& SplitterMacro() const
    {
      return splitter_macro_;
    };

    //! prepare time integrator specific things before calculation of initial potential field
    virtual void pre_calc_initial_potential_field() = 0;

    //! clean up settings from pre_calc_initial_potential_field() after initial potential field is
    //! calculated
    virtual void post_calc_initial_potential_field() = 0;

   protected:
    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! add parameters depending on the problem, i.e., loma, level-set, ...
    void add_problem_specific_parameters_and_vectors(Teuchos::ParameterList& params) override;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! assemble global system of equations
    void AssembleMatAndRHS() override;

    //! prepare time loop
    void PrepareTimeLoop() override;

    void PrepareTimeStep() override;

    void prepare_first_time_step() override;

    //! initialize meshtying strategy (including standard case without meshtying)
    void CreateScalarHandler() override;

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! compute time step size
    void ComputeTimeStepSize(double& dt) final;

    //! temperature is computed based on function prescribed in input file
    double compute_temperature_from_function() const;

    //! evaluate SOC and c-rate of electrode
    void evaluate_electrode_info_interior();

    //! evaluate cell voltage of electrode
    void EvaluateCellVoltage();

    //! Evaluate cccv phase based on c rate and cell voltage
    void EvaluateCCCVPhase();

    //! extrapolate current state and adapt time step
    double extrapolate_state_adapt_time_step(double dt);

    //! Parameter check for diffusion-conduction formulation
    void valid_parameter_diff_cond();

    //! Initialize Nernst-BC
    void InitNernstBC();

    //! initialize meshtying strategy (including standard case without meshtying)
    void create_meshtying_strategy() override;

    //! set up concentration-potential splitter
    void SetupConcPotSplit();

    //! set up concentration-potential-potential splitter for macro scale in multi-scale
    //! simulations
    void setup_conc_pot_pot_split();

    //! reduces the dimension of the null space by one (if system matrix is partitioned according to
    //! concentration and potential). The original full null space was computed for all degrees of
    //! freedom on the discretization, such that the reduced null spaces still have the full
    //! dimension. Thus, the dimension of each null space is decreased by one, and the corresponding
    //! zero null space vector is removed from the null space.
    void reduce_dimension_null_space_blocks(
        Teuchos::RCP<CORE::LINALG::Solver> solver, int init_block_number) const;

    /*--- calculate and update -----------------------------------------------*/

    //! calculate initial electric potential field
    virtual void calc_initial_potential_field();

    //! \brief computes different conductivity expressions for electrolyte solutions (ELCH)
    //!
    //! \param sigma       result vector
    //! \param effCond     flag for computation of effective conductivity
    //! \param specresist  flag for computation of specific electrolyte resistance
    //! \return            specific resistance
    double ComputeConductivity(
        CORE::LINALG::SerialDenseVector& sigma, bool effCond = false, bool specresist = false);

    //! apply galvanostatic control (update electrode potential) (ELCH)
    bool apply_galvanostatic_control();

    //! \brief evaluate domain or boundary conditions for electrode kinetics
    //!
    //! \param systemmatrix    global system matrix
    //! \param rhs             global right-hand side vector
    //! \param condstring      name of condition to be evaluated
    void evaluate_electrode_kinetics_conditions(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs,
        const std::string& condstring);

    //! \brief evaluate point boundary conditions for electrode kinetics
    //!
    //! \param systemmatrix  global system matrix
    //! \param rhs           global right-hand side vector
    void evaluate_electrode_boundary_kinetics_point_conditions(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs);

    //! Add Linearization for Nernst-BC
    void linearization_nernst_condition();

    //! update time-dependent electrode state variables at the end of an time step
    virtual void electrode_kinetics_time_update() = 0;

    void evaluate_solution_depending_conditions(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<Epetra_Vector> rhs) override;

    void ApplyDirichletBC(
        double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt) override;

    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads) override;

    void perform_aitken_relaxation(
        Epetra_Vector& phinp, const Epetra_Vector& phinp_inc_diff) override;

    /*--- query and output ---------------------------------------------------*/

    //! check for negative values of concentrations (ELCH)
    void check_concentration_values(
        Teuchos::RCP<Epetra_Vector> vec  //!< current phi vector to be checked
    );

    void OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux, const std::string& fluxtype) override;

    /*========================================================================*/
    //! @name ELCH variables
    /*========================================================================*/

    //! the parameter list for elch problems
    Teuchos::RCP<const Teuchos::ParameterList> elchparams_;

    //! type of closing equation for electric potential
    INPAR::ELCH::EquPot equpot_;

    //! ELCH-specific parameter F/R
    const double fr_;

    //! function number describing the temporal temperature curve
    const int temperature_funct_num_;

    //! homogeneous temperature within the scalar transport field (can be time dependent)
    double temperature_;

    //! number of iterations in galvanostatic mode (ELCH)
    int gstatnumite_;

    //! value of electric potential increment in galvanostatic mode (ELCH)
    double gstatincrement_;

    //! flag for (de)activation of double layer capacity
    bool dlcapexists_;

    //! electro-kinetics toggle
    //! Toggle which defines dof's with Nernst-BC or Dirichlet condition
    Teuchos::RCP<Epetra_Vector> ektoggle_;

    //! dirichlet toggle
    //! Toggle which defines dof's with a Dirichlet condition
    Teuchos::RCP<Epetra_Vector> dctoggle_;

    //! initial volumes of resolved electrodes
    std::map<int, double> electrodeinitvols_;

    //! states of charge of resolved electrodes
    std::map<int, double> electrodesoc_;

    //! C rates with respect to resolved electrodes
    std::map<int, double> electrodecrates_;

    //! mean reactant concentrations at electrode boundaries
    std::map<int, double> electrodeconc_;

    //! mean electric overpotentials at electrode boundaries
    std::map<int, double> electrodeeta_;

    //! total electric currents at electrode boundaries
    std::map<int, double> electrodecurr_;

    //! voltage at both conditions
    std::map<int, double> electrodevoltage_;

    //! cell voltage
    double cellvoltage_;

    //! cell voltage from previous time step
    double cellvoltage_old_;

    Teuchos::RCP<SCATRA::CCCVCondition> cccv_condition_;

    //! cell C rate
    double cellcrate_;

    //! cell C rate from previous time step
    double cellcrate_old_;

    //! modified time step size for CCCV cell cycling
    const double cycling_timestep_;

    //! flag indicating modified time step size for CCCV cell cycling
    bool adapted_timestep_active_;

    //! adapted time step
    double dt_adapted_;

    //! time step number of last modification of time step size
    int last_dt_change_;

    //! map extractor for macro scale in multi-scale simulations
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> splitter_macro_;

    //! CSV writers for SOC, c-rate and cell voltage
    std::map<int, std::optional<IO::RuntimeCsvWriter>> runtime_csvwriter_soc_;
    std::optional<IO::RuntimeCsvWriter> runtime_csvwriter_cell_voltage_;
  };  // class ScaTraTimIntElch


  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Helper class for managing different number of degrees of freedom per node
   */
  class ScalarHandlerElch : public ScalarHandler
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScalarHandlerElch();

    //! initialize time integration
    void Setup(const ScaTraTimIntImpl* scatratimint) override;

    /*========================================================================*/
    //! @name Access and Query methods
    /*========================================================================*/

    //! return maximum number of dofs per node
    int NumDofPerNode() const override
    {
      CheckIsSetup();
      return *(numdofpernode_.rbegin());
    };

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    int NumScal() const override
    {
      CheckIsSetup();
      return *(numscal_.rbegin());
    }

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    int NumScalInCondition(const CORE::Conditions::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const override;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
   protected:
    //! number of transported scalars (without potential and current density)
    std::set<int> numscal_;

  };  // class ScalarHandlerElch

}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif