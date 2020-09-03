/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_ELCH_H
#define SCATRA_TIMINT_ELCH_H

#include "../drt_inpar/inpar_elch.H"

#include "../drt_scatra/scatra_timint_implicit.H"

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace FLD
{
  class Meshtying;
}

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
    ScaTraTimIntElch(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
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
    void BuildBlockMaps(const std::vector<Teuchos::RCP<DRT::Condition>>& partitioningconditions,
        std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const override;

    void BuildBlockNullSpaces(
        Teuchos::RCP<LINALG::Solver> solver, int init_block_number) const override;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! Set elch-specific parameters
    void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const override;

    //! contains the nonlinear iteration loop
    void NonlinearSolve() override;

    //! calculate error compared to analytical solution
    void EvaluateErrorComparedToAnalyticalSol() override;

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! setup natural convection
    void SetupNatConv() override;

    /*--- calculate and update -----------------------------------------------*/

    //! determine whether there are still time steps to be evaluated
    bool NotFinished() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(int num) override;

    /*--- query and output ---------------------------------------------------*/

    //! output solution and restart data to file
    void Output(int num) final;

    //! problem-specific outputs
    void OutputProblemSpecific() override;

    //! read problem-specific restart data
    void ReadRestartProblemSpecific(int step, IO::DiscretizationReader& reader) override;

    //! output electrode domain status information to screen and file
    void OutputElectrodeInfoDomain();

    //! output electrode boundary status information to screen and file
    void OutputElectrodeInfoBoundary();

    //! evaluate status information on single line or surface electrode
    Teuchos::RCP<Epetra_SerialDenseVector> EvaluateSingleElectrodeInfo(
        int condid,                    //!< ID of condition to be evaluated
        const std::string& condstring  //!< name of condition to be evaluated
    );

    //! evaluate status information on single point electrode
    Teuchos::RCP<Epetra_SerialDenseVector> EvaluateSingleElectrodeInfoPoint(
        Teuchos::RCP<DRT::Condition> condition  //!< condition to be evaluated
    );

    //! post-process status information on single electrode
    void PostProcessSingleElectrodeInfo(
        Epetra_SerialDenseVector&
            scalars,  //!< scalar quantities associated with electrode status information (in)
        unsigned id,  //!< electrode ID (in)
        bool print,   //!< flag for output to screen and file (in)
        double& currentsum,    //!< net current involving all conditions (out)
        double& currtangent,   //!< tangent of current w.r.t. electrode potential (out)
        double& currresidual,  //!< negative residual of current equation (out)
        double& electrodeint,  //!< physical dimensions of the electrode region (out)
        double& electrodepot,  //!< electrode potential on electrode side (out)
        double& meanoverpot    //!< mean overpotential (out)
    );

    //! output electrode interior status information to screen and files
    void OutputElectrodeInfoInterior();

    //! output cell voltage to screen and file
    void OutputCellVoltage();

    //! output restart data
    void OutputRestart() const override;

    //! output type of closing equation for electric potential
    INPAR::ELCH::EquPot EquPot() const { return equpot_; }

    //! return constant F/RT
    double FRT() const { return frt_; }

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
    const Teuchos::RCP<const LINALG::MultiMapExtractor>& SplitterMacro() const
    {
      return splitter_macro_;
    };

   protected:
    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! add parameters depending on the problem, i.e., loma, level-set, ...
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params) override;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! assemble global system of equations
    void AssembleMatAndRHS() override;

    //! prepare time loop
    void PrepareTimeLoop() override;

    //! initialization procedure prior to evaluation of first time step
    void PrepareFirstTimeStep() override;

    //! initialize meshtying strategy (including standard case without meshtying)
    void CreateScalarHandler() override;

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! compute time step size
    void ComputeTimeStepSize(double& dt) final;

    //! extrapolate current state and adapt time step
    double ExtrapolateStateAdaptTimeStep(double dt);

    //! Parameter check for diffusion-conduction formulation
    void ValidParameterDiffCond();

    //! Initialize Nernst-BC
    void InitNernstBC();

    //! initialize meshtying strategy (including standard case without meshtying)
    void CreateMeshtyingStrategy() override;

    //! set up concentration-potential splitter
    void SetupConcPotSplit();

    //! set up concentration-potential-potential splitter for macro scale in multi-scale
    //! simulations
    void SetupConcPotPotSplit();

    //! reduces the dimension of the null space by one (if system matrix is partitioned according to
    //! concentration and potential). The original full null space was computed for all degrees of
    //! freedom on the discretization, such that the reduced null spaces still have the full
    //! dimension. Thus, the dimension of each null space is decreased by one, and the corresponding
    //! zero null space vector is removed from the null space.
    void ReduceDimensionNullSpaceBlocks(
        Teuchos::RCP<LINALG::Solver> solver, int init_block_number) const;

    /*--- calculate and update -----------------------------------------------*/

    //! calculate initial electric potential field
    //! this function is never called directly, but only via overloading
    virtual void CalcInitialPotentialField();

    //! computes different conductivity expressions for electrolyte solutions (ELCH)
    double ComputeConductivity(Epetra_SerialDenseVector& sigma,  //! result vector
        bool effCond = false,    //! flag for computation of effective conductivity
        bool specresist = false  //! flag for computation of specific electrolyte resistance
    );

    //! apply galvanostatic control (update electrode potential) (ELCH)
    bool ApplyGalvanostaticControl();

    //! evaluate domain or boundary conditions for electrode kinetics
    void EvaluateElectrodeKineticsConditions(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< global system matrix
        Teuchos::RCP<Epetra_Vector> rhs,                    //!< global right-hand side vector
        const std::string& condstring                       //!< name of condition to be evaluated
    );

    //! evaluate point boundary conditions for electrode kinetics
    void EvaluateElectrodeBoundaryKineticsPointConditions(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< global system matrix
        Teuchos::RCP<Epetra_Vector> rhs                     //!< global right-hand side vector
    );

    //! Add Linearization for Nernst-BC
    void LinearizationNernstCondition();

    //! update time-dependent electrode state variables at the end of an time step
    virtual void ElectrodeKineticsTimeUpdate() = 0;

    void EvaluateSolutionDependingConditions(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<Epetra_Vector> rhs) override;

    void ApplyDirichletBC(
        double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt) override;

    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads) override;

    void PerformAitkenRelaxation(
        Epetra_Vector& phinp, const Epetra_Vector& phinp_inc_diff) override;

    /*--- query and output ---------------------------------------------------*/

    //! check for negative values of concentrations (ELCH)
    void CheckConcentrationValues(
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

    //! ELCH-specific parameter F/RT (default value for the temperature is 298K)
    double frt_;

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

    //! map extractor for macro scale in multi-scale simulations
    Teuchos::RCP<const LINALG::MultiMapExtractor> splitter_macro_;
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
    int NumScalInCondition(const DRT::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const override;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
   protected:
    //! number of transported scalars (without potential and current density)
    std::set<int> numscal_;

  };  // class ScalarHandlerElch

}  // namespace SCATRA
#endif