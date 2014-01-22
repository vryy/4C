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

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_inpar/inpar_elch.H"

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

  //! Initialization of system matrix
  virtual void InitSystemMatrix();

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! set elch-specific parameter
  void SetElementGeneralScaTraParameter();

  //! contains the nonlinear iteration loop
  void NonlinearSolve();

  //! calculate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /*========================================================================*/
  //! @name ELCH methods
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  //! compute density field from concentration(s) for ELCH natural convection
  void ComputeDensity();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  virtual void Update(const int num=0);

  //! update density field for ELCH natural convection
  virtual void UpdateDensityElch() = 0;

  /*--- query and output ---------------------------------------------------*/

  //! return density field at time n+1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNp() { return elchdensnp_; }

  //! return density field at time n for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchN() { return elchdensn_; }

  //! return density field at time n-1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNm() { return elchdensnm_; }

  //! output of all electrode information to screen and/or file (default: 2x true)
  Teuchos::RCP< std::vector<double> > OutputElectrodeInfo(
    bool printtoscreen = true,  //!< flag for printing information to screen
    bool printtofile = true     //!< flag for printing information to file
  );

  //! output of electrode information for a single boundary condition
  Teuchos::RCP< std::vector<double> > OutputSingleElectrodeInfo(
    DRT::Condition* condition, //!< pointer to the condition
    const int  condid,         //!< condition id
    const bool printtoscreen,  //!< flag for output to screen
    const bool printtofile,    //!< flag for output to file
    double& currentsum,        //!< net current over all boundaries (out)
    double& pottangent,        //!< tangent of current w.r.t. electrode potential (out)
    double& currresidual,      //!< neg. residual of current equation (out)
    double& electrodesurface,  //!< surface of the electrode
    double& electrodepot,      //!< electrode potential (metal side)
    double& meanoverpot        //!< mean surface overpotential
  );

  protected:

  /*========================================================================*/
  //! @name set element parameters
  /*========================================================================*/

  //! add parameters depending on the problem, i.e., loma, level-set, ...
  void AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params //!< parameter list
  );

  //! add parameters depending on the problem, i.e., loma, level-set, ... with special modifications for
  //! calculation fo initial phidt
  void AddProblemSpecificParametersAndVectorsForCalcInitialPhiDt(
    Teuchos::ParameterList& params //!< parameter list
  ){return;};

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

  /*========================================================================*/
  //! @name ELCH methods
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! setup of natural convection
  void SetupElchNatConv();

  //! Parameter check for diffusion-conduction formulation
  void ValidParameterDiffCond();

  //! Initialize Nernst-BC
  void InitNernstBC();

  /*--- calculate and update -----------------------------------------------*/

  //! calculate initial field for electric potential based on ion concentrations (ELCH)
  void CalcInitialPotentialField();

  //! computes different conductivity expressions for electrolyte solutions (ELCH)
  Teuchos::RCP<Epetra_SerialDenseVector>  ComputeConductivity();

  //! apply galvanostatic control (update electrode potential) (ELCH)
  bool ApplyGalvanostaticControl();

  //! compute contribution of electrode kinetics to eq. system
  void EvaluateSolutionDependingBC(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  );

  //! Add Linearization for Nernst-BC
  void  LinearizationNernstCondition();

  //! update time-dependent electrode state variables at the end of an time step
  virtual void ElectrodeKineticsTimeUpdate() = 0;

  /*--- query and output ---------------------------------------------------*/

  //! check for negative values of concentrations (ELCH)
  void CheckConcentrationValues(
    Teuchos::RCP<Epetra_Vector> vec //!< current phi vector to be checked
  );

  /*========================================================================*/
  //! @name ELCH variables
  /*========================================================================*/

  //! the parameter list for elch problems
  Teuchos::RCP<Teuchos::ParameterList> elchparams_;

  //! elch type
  INPAR::ELCH::ElchType elchtype_;

  //! ELCH-specific parameter F/RT (default value for the temperature is 298K)
  double frt_;

  //! density (nodebased!) at times n-1 for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensnm_;
  //! density (nodebased!) at times n for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensn_;
  //! density (nodebased!) at times n+1 for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensnp_;

  //! number of iterations in galvanostatic mode (ELCH)
  int gstatnumite_;

  //! value of electric potential increment in galvanostatic mode (ELCH)
  double gstatincrement_;

  //! initial concentrations of the electrolyte solution (ELCH)
  std::vector<double> c0_;

  //! conductivity of the single species and the whole electrolyte
  Teuchos::RCP<Epetra_SerialDenseVector> sigma_;

  //! densification coefficients for ELCH natural convection
  std::vector<double> densific_;

  //! flag for (de)activation of double layer capacity
  bool dlcapexists_;

  //! electro-kinetics toggle
  // Toggle defining the dof's with Nernst-BC
  Teuchos::RCP<Epetra_Vector> ektoggle_;
};

};

#endif
