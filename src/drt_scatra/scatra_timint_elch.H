/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_elch.H

\brief scatra time integration for elch

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
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

/*==========================================================================*/
//! Algorithm and time integration classes for scalar transport problems: elch
/*==========================================================================*/
namespace SCATRA
{
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

    //! virtual destructor
    virtual ~ScaTraTimIntElch() { return; };

    //! initialize algorithm
    virtual void Init();

    //! initialize algorithm
    virtual void Setup();

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/

    void SetupSplitter();

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! Set elch-specific parameters
    void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const;

    //! contains the nonlinear iteration loop
    void NonlinearSolve();

    //! calculate error compared to analytical solution
    void EvaluateErrorComparedToAnalyticalSol();

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! setup natural convection
    void SetupNatConv();

    /*--- calculate and update -----------------------------------------------*/

    //! determine whether there are still time steps to be evaluated
    bool NotFinished();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    virtual void Update(const int num = 0);

    /*--- query and output ---------------------------------------------------*/

    //! output solution and restart data to file
    void Output(const int num = 0) final;

    //! problem-specific outputs
    void OutputProblemSpecific();

    //! read problem-specific restart data
    void ReadRestartProblemSpecific(const int step, IO::DiscretizationReader& reader);

    //! output electrode domain status information to screen and file
    void OutputElectrodeInfoDomain();

    //! output electrode boundary status information to screen and file
    void OutputElectrodeInfoBoundary();

    //! evaluate status information on single line or surface electrode
    Teuchos::RCP<Epetra_SerialDenseVector> EvaluateSingleElectrodeInfo(
        const unsigned condid,        //!< ID of condition to be evaluated
        const std::string condstring  //!< name of condition to be evaluated
    );

    //! evaluate status information on single point electrode
    Teuchos::RCP<Epetra_SerialDenseVector> EvaluateSingleElectrodeInfoPoint(
        Teuchos::RCP<DRT::Condition> condition  //!< condition to be evaluated
    );

    //! post-process status information on single electrode
    void PostProcessSingleElectrodeInfo(
        Epetra_SerialDenseVector&
            scalars,        //!< scalar quantities associated with electrode status information (in)
        const unsigned id,  //!< electrode ID (in)
        const bool print,   //!< flag for output to screen and file (in)
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
    virtual void OutputRestart() const;

    //! output type of closing equation for electric potential
    INPAR::ELCH::EquPot EquPot() const { return equpot_; }

    //! return constant F/RT
    double FRT() const { return frt_; }

    //! return elch parameter list
    const Teuchos::RCP<const Teuchos::ParameterList> ElchParameterList() const
    {
      return elchparams_;
    }

    //! return states of charge of resolved electrodes
    const Teuchos::RCP<std::vector<double>>& ElectrodeSOC() const { return electrodesoc_; };

    //! return C rates with respect to resolved electrodes
    const Teuchos::RCP<std::vector<double>>& ElectrodeCRates() const { return electrodecrates_; };

    //! return mean reactant concentrations at electrode boundaries
    const Teuchos::RCP<std::vector<double>>& ElectrodeConc() const { return electrodeconc_; };

    //! return mean electric overpotentials at electrode boundaries
    const Teuchos::RCP<std::vector<double>>& ElectrodeEta() const { return electrodeeta_; };

    //! return total electric currents at electrode boundaries
    const Teuchos::RCP<std::vector<double>>& ElectrodeCurr() const { return electrodecurr_; };

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
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params  //!< parameter list
    );

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! prepare time loop
    void PrepareTimeLoop();

    //! initialization procedure prior to evaluation of first time step
    void PrepareFirstTimeStep();

    //! initialize meshtying strategy (including standard case without meshtying)
    virtual void CreateScalarHandler();

    /*========================================================================*/
    //! @name ELCH methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! compute time step size
    void ComputeTimeStepSize(double& dt) final;

    //! Parameter check for diffusion-conduction formulation
    void ValidParameterDiffCond();

    //! Initialize Nernst-BC
    void InitNernstBC();

    //! initialize meshtying strategy (including standard case without meshtying)
    void CreateMeshtyingStrategy();

    //! set up concentration-potential splitter
    void SetupConcPotSplit();

    //! set up concentration-potential-potential splitter for macro scale in multi-scale simulations
    void SetupConcPotPotSplit();

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
        const std::string condstring                        //!< name of condition to be evaluated
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

    //! compute contributions of solution-depending boundary and interface conditions to global
    //! system of equations
    void EvaluateSolutionDependingConditions(
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                     //!< rhs vector
    );

    //! apply Dirichlet boundary conditions to state vectors
    void ApplyDirichletBC(const double time,  //!< time
        Teuchos::RCP<Epetra_Vector> phinp,  //!< state vector of scalar transport degrees of freedom
        Teuchos::RCP<Epetra_Vector>
            phidt  //!< state vector of time derivatives of scalar transport degrees of freedom
    );

    //! evaluate Neumann boundary conditions
    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    //! perform Aitken relaxation
    void PerformAitkenRelaxation(Epetra_Vector& phinp,  //!< state vector to be relaxed
        const Epetra_Vector&
            phinp_inc_diff  //!< difference between current and previous state vector increments
    );

    /*--- query and output ---------------------------------------------------*/

    //! check for negative values of concentrations (ELCH)
    void CheckConcentrationValues(
        Teuchos::RCP<Epetra_Vector> vec  //!< current phi vector to be checked
    );

    //! output flux vectors
    void OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux,  //!< flux vector
        const std::string& fluxtype                         //!< flux type ("domain" or "boundary")
    );

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
    Teuchos::RCP<std::vector<double>> electrodeinitvols_;

    //! states of charge of resolved electrodes
    Teuchos::RCP<std::vector<double>> electrodesoc_;

    //! C rates with respect to resolved electrodes
    Teuchos::RCP<std::vector<double>> electrodecrates_;

    //! mean reactant concentrations at electrode boundaries
    Teuchos::RCP<std::vector<double>> electrodeconc_;

    //! mean electric overpotentials at electrode boundaries
    Teuchos::RCP<std::vector<double>> electrodeeta_;

    //! total electric currents at electrode boundaries
    Teuchos::RCP<std::vector<double>> electrodecurr_;

    //! cell voltage
    double cellvoltage_;

    //! cell voltage from previous time step
    double cellvoltage_old_;

    //! cutoff voltage for CCCV cell cycling
    double cutoff_voltage_;

    //! cell C rate
    double cellcrate_;

    //! flag indicating whether cell is currently being charged or discharged
    bool charging_;

    //! flag indicating whether cell is currently being operated in constant-current (CC),
    //! constant-voltage (CV), or relaxation (RX) mode
    INPAR::ELCH::CCCVCyclingMode mode_cccv_;

    //! ID of constant-current constant-voltage (CCCV) half-cycle condition currently in effect
    int condid_cccv_;

    //! total number of charge and discharge half-cycles
    int nhalfcycles_;

    //! number of current charge or discharge half-cycle
    int ihalfcycle_;

    //! time step at the start of current charge or discharge half-cycle
    int ihalfcycle_startstep_;

    //! end time of current relaxation phase
    double relax_endtime_;

    //! modified time step size for CCCV cell cycling
    const double cycling_timestep_;

    //! flag indicating modified time step size for CCCV cell cycling
    bool cycling_timestep_active_;

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

    //! Destructor
    virtual ~ScalarHandlerElch(){};

    //! initialize time integration
    void Setup(const ScaTraTimIntImpl* const scatratimint);

    /*========================================================================*/
    //! @name Access and Query methods
    /*========================================================================*/

    //! return maximum number of dofs per node
    int NumDofPerNode() const
    {
      CheckIsSetup();
      return *(numdofpernode_.rbegin());
    };

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    int NumScal() const
    {
      CheckIsSetup();
      return *(numscal_.rbegin());
    }

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    int NumScalInCondition(const DRT::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
   protected:
    //! number of transported scalars (without potential and current density)
    std::set<int> numscal_;

  };  // class ScalarHandlerElch

};  // namespace SCATRA
#endif
