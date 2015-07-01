/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_elch.H

\brief scatra time integration for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
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
  ScaTraTimIntElch(
      Teuchos::RCP<DRT::Discretization>        dis,
      Teuchos::RCP<LINALG::Solver>             solver,
      Teuchos::RCP<Teuchos::ParameterList>     params,
      Teuchos::RCP<Teuchos::ParameterList>     sctratimintparams,
      Teuchos::RCP<Teuchos::ParameterList>     extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>   output);

  //! virtual destructor
  virtual ~ScaTraTimIntElch(){return;};

  //! initialize algorithm
  virtual void Init();

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

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  virtual void Update(const int num=0);

  /*--- query and output ---------------------------------------------------*/

  //! problem-specific outputs
  void OutputProblemSpecific();

  //! output of all electrode boundary status information to screen and file
  void OutputElectrodeInfoBoundary();

  //! output of electrode information for a single boundary condition
  void OutputSingleElectrodeInfoBoundary(
    DRT::Condition* condition, //!< pointer to the condition
    const int  condid,         //!< condition id
    const bool print,          //!< flag for output to screen and file
    double& currentsum,        //!< net current over all boundaries (out)
    double& pottangent,        //!< tangent of current w.r.t. electrode potential (out)
    double& currresidual,      //!< neg. residual of current equation (out)
    double& electrodesurface,  //!< surface of the electrode
    double& electrodepot,      //!< electrode potential (metal side)
    double& meanoverpot        //!< mean surface overpotential
  );

  //! output electrode interior status information to screen and files
  void OutputElectrodeInfoInterior();

  //! output cell voltage to screen
  void OutputCellVoltage();

  //! output type of closing equation for electric potential
  const INPAR::ELCH::EquPot EquPot() const { return equpot_; }

  //! return constant F/RT
  const double FRT() const {return frt_;}

  //! return elch parameter list
  Teuchos::RCP<Teuchos::ParameterList> ElchParameterList() { return elchparams_;}

  //! return states of charge of resolved electrodes
  const Teuchos::RCP<std::vector<double> >& ElectrodeSOC() const {return electrodesoc_;};

  //! return mean reactant concentrations at electrode boundaries
  const Teuchos::RCP<std::vector<double> >& ElectrodeConc() const {return electrodeconc_;};

  //! return mean electric overpotentials at electrode boundaries
  const Teuchos::RCP<std::vector<double> >& ElectrodeEta() const {return electrodeeta_;};

  //! return total electric currents at electrode boundaries
  const Teuchos::RCP<std::vector<double> >& ElectrodeCurr() const {return electrodecurr_;};

  //! return cell voltage
  const double& CellVoltage() const {return cellvoltage_;};

  protected:

  /*========================================================================*/
  //! @name set element parameters
  /*========================================================================*/

  //! add parameters depending on the problem, i.e., loma, level-set, ...
  void AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params //!< parameter list
  );

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! Calculate problem specific norm
  void CalcProblemSpecificNorm(
      double& conresnorm,
      double& incconnorm_L2,
      double& connorm_L2,
      double& incpotnorm_L2,
      double& potnorm_L2,
      double& potresnorm,
      double& conresnorminf);

  //! call evaluate routine
  void Evaluate(Teuchos::ParameterList eleparams);

  /*========================================================================*/
  //! @name ELCH methods
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! Parameter check for diffusion-conduction formulation
  void ValidParameterDiffCond();

  //! Initialize Nernst-BC
  void InitNernstBC();

  //! initialize meshtying strategy (including standard case without meshtying)
  void CreateMeshtyingStrategy();

  //! adapt number of transported scalars
  void AdaptNumScal();

  /*--- calculate and update -----------------------------------------------*/

  //! calculate initial field for electric potential based on ion concentrations (ELCH)
  void CalcInitialPotentialField();

  //! computes different conductivity expressions for electrolyte solutions (ELCH)
  double  ComputeConductivity(bool effCond=false, bool specresist=false);

  //! apply galvanostatic control (update electrode potential) (ELCH)
  bool ApplyGalvanostaticControl();

  //! evaluate domain or boundary conditions for electrode kinetics
  void EvaluateElectrodeKineticsConditions(
      Teuchos::RCP<LINALG::SparseOperator>   systemmatrix,   //!< global system matrix
      Teuchos::RCP<Epetra_Vector>            rhs,            //!< global right-hand side vector
      const std::string                      condstring      //!< name of condition to be evaluated
  );

  //! Add Linearization for Nernst-BC
  void  LinearizationNernstCondition();

  //! update time-dependent electrode state variables at the end of an time step
  virtual void ElectrodeKineticsTimeUpdate() = 0;

  //! compute contributions of solution-depending boundary and interface conditions to global system of equations
  void EvaluateSolutionDependingConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,      //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs                //!< rhs vector
  );

  /*--- query and output ---------------------------------------------------*/

  //! check for negative values of concentrations (ELCH)
  void CheckConcentrationValues(
    Teuchos::RCP<Epetra_Vector> vec //!< current phi vector to be checked
  );

  //! print header of convergence table to screen
  void PrintConvergenceHeader();

  //! print first line of convergence table to screen
  void PrintConvergenceValuesFirstIter(
      const int&              itnum,          //!< current Newton-Raphson iteration step
      const int&              itemax,         //!< maximum number of Newton-Raphson iteration steps
      const double&           ittol,          //!< relative tolerance for Newton-Raphson scheme
      const double&           conresnorm,     //!< L2 norm of concentration residual
      const double&           potresnorm,     //!< L2 norm of potential residual (only relevant for electrochemistry)
      const double&           conresnorminf   //!< infinity norm of concentration residual
  );

  //! print current line of convergence table to screen
  void PrintConvergenceValues(
      const int&              itnum,           //!< current Newton-Raphson iteration step
      const int&              itemax,          //!< maximum number of Newton-Raphson iteration steps
      const double&           ittol,           //!< relative tolerance for Newton-Raphson scheme
      const double&           conresnorm,      //!< L2 norm of concentration residual
      const double&           potresnorm,      //!< L2 norm of potential residual (only relevant for electrochemistry)
      const double&           incconnorm_L2,   //!< L2 norm of concentration increment
      const double&           connorm_L2,      //!< L2 norm of concentration state vector
      const double&           incpotnorm_L2,   //!< L2 norm of potential increment
      const double&           potnorm_L2,      //!< L2 norm of potential state vector
      const double&           conresnorminf    //!< infinity norm of concentration residual
  );

  //! print finish line of convergence table to screen
  void PrintConvergenceFinishLine();

  /*========================================================================*/
  //! @name ELCH variables
  /*========================================================================*/

  //! the parameter list for elch problems
  Teuchos::RCP<Teuchos::ParameterList> elchparams_;

  //! type of closing equation for electric potential
  INPAR::ELCH::EquPot equpot_;

  //! ELCH-specific parameter F/RT (default value for the temperature is 298K)
  double frt_;

  //! number of iterations in galvanostatic mode (ELCH)
  int gstatnumite_;

  //! value of electric potential increment in galvanostatic mode (ELCH)
  double gstatincrement_;

  //! conductivity of the single species and the whole electrolyte
  Teuchos::RCP<Epetra_SerialDenseVector> sigma_;

  //! flag for (de)activation of double layer capacity
  bool dlcapexists_;

  //! electro-kinetics toggle
  //! Toggle which defines dof's with Nernst-BC or Dirichlet condition
  Teuchos::RCP<Epetra_Vector> ektoggle_;

  //! dirichlet toggle
  //! Toggle which defines dof's with a Dirichlet condition
  Teuchos::RCP<Epetra_Vector> dctoggle_;

  //! states of charge of resolved electrodes
  Teuchos::RCP<std::vector<double> > electrodesoc_;

  //! mean reactant concentrations at electrode boundaries
  Teuchos::RCP<std::vector<double> > electrodeconc_;

  //! mean electric overpotentials at electrode boundaries
  Teuchos::RCP<std::vector<double> > electrodeeta_;

  //! total electric currents at electrode boundaries
  Teuchos::RCP<std::vector<double> > electrodecurr_;

  //! cell voltage
  double cellvoltage_;
};
};
#endif
