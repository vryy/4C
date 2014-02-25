/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.H
\brief Implicit time integration for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_IMPL_H
#define STRTIMINT_IMPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>
#include <NOX_Epetra_Group.H>

#include "strtimint.H"

// forward declarations
namespace DRT
{
  class Condition;
}

namespace LINALG
{
  class MultiMapExtractor;
  class KrylovProjector;
}

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  namespace AUX
  {
    class MapExtractor;
  }

  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics
   *        with \b implicit time integration
   *
   * <h3> About </h3>
   * The implicit time integrator object is a derivation of the base time integrators with an eye towards implicit time integration. #TimIntImpl provides the environment needed to execute implicit integrators. This is chiefly the non-linear solution technique, e.g., Newton-Raphson iteration. These iterative solution techniques require a set of control parameters which are stored within this object. It is up to derived object to implement the time-space discretised residuum and its tangent. This object provides some utility functions to obtain various force vectors necessary in the calculation of the force residual in the derived time integrators.
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntImpl : public TimInt,
                     public NOX::Epetra::Interface::Required,
                     public NOX::Epetra::Interface::Jacobian,
                     public NOX::Epetra::Interface::Preconditioner
  {

  public:

    //! @name Construction
    //@{

    //! Constructor
    TimIntImpl
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntImpl() { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //! return time integration factor
    virtual const double TimIntParam() const = 0;

    //! return Newmark parameter \f$beta\f$
    virtual const double GetNewmarkBeta() const = 0;

    //@}

    //! Do time integration of single step
    int IntegrateStep();

    //! @name Prediction
    //@{

    //! Predict target solution and identify residual
    void Predict();

    //! Identify residual
    //! This method does not predict the target solution but
    //! evaluates the residual and the stiffness matrix.
    //! In partitioned solution schemes, it is better to keep the current
    //! solution instead of evaluating the initial guess (as the predictor)
    //! does.
    void PreparePartitionStep();

    //! Predict constant displacements, velocities and accelerations,
    //! i.e. the initial guess is equal to the last converged step
    //! except Dirichlet BCs
    void PredictConstDisVelAcc();

    //! Predict constant displacements, however the velocities
    //! and accelerations are consistent to the time integration
    //! if the constant displacements are taken as correct displacement
    //! solution.
    //! This method has to be implemented by the individual time
    //! integrator.
    virtual void PredictConstDisConsistVelAcc() = 0;

    //! Predict displacements based on the assumption of constant
    //! velocities. Calculate consistent velocities and accelerations
    //! afterwards.
    //! This method has to be implemented by the individual time
    //! integrator.
    virtual void PredictConstVelConsistAcc() = 0;

    //! Predict displacements based on the assumption of constant
    //! accelerations. Calculate consistent velocities and accelerations
    //! afterwards.
    //! This method has to be implemented by the individual time
    //! integrator.
    virtual void PredictConstAcc() = 0;

    //! Predict displacements which satisfy exactly the Dirichlet BCs
    //! and the linearised system at the previously converged state.
    //!
    //! This is an implicit predictor, i.e. it calls the solver once.
    void PredictTangDisConsistVelAcc();

    //!
    void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
    //!
    void UpdateKrylovSpaceProjection();

    //@}

    //! @name Forces and stiffnesses
    //@{

    //! Do residual force due to global balance of momentum
    //! and its tangent with respect to the current
    //! displacements \f$D_{n+1}\f$
    //!
    //! This is <i>the</i> central method which is different for each
    //! derived implicit time integrator. The time integrator implementation
    //! is expected to set members #fres_ and #stiff_.
    //! The residual #fres_ is expected to follow the <i>same</i> sign
    //! convention like its tangent #stiff_, i.e. to use
    //! Newton--Raphson's method the residual will be scaled by -1.
    //!
    //! For structural contact simulations, it is of importance
    //! whether this method is called from a predictor step or
    //! from a regular Newton step. Thus we hand in a boolean flag
    //! which reflects the respective situation.
    virtual void EvaluateForceStiffResidual(Teuchos::ParameterList& params) = 0;

    //! Evaluate forces due to constraints
    void ApplyForceStiffConstraint
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< last evaluated displacements
      const Teuchos::RCP<Epetra_Vector> disn,  //!< evaluation displacements
      Teuchos::RCP<Epetra_Vector>& fint,  //!< forces are added onto
      Teuchos::RCP<LINALG::SparseOperator>& stiff,  //!< stiffness is added onto
      Teuchos::ParameterList pcon   //!< Parameterlist containing scale factors for matrix entries
    );

    //! Evaluate forces due to Windkessel bcs
    void ApplyForceStiffWindkessel
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< last evaluated displacements
      const Teuchos::RCP<Epetra_Vector> disn,  //!< evaluation displacements
      Teuchos::ParameterList pcon   //!< Parameterlist containing scale factors for matrix entries
    );

    //! Evaluate _certain_ surface stresses and stiffness
    //!
    //! The Evaluation happens internal-force-like.
    void ApplyForceStiffSurfstress
    (
      const double time,                         //!< evaluation time
      const double dt,                           //!< time step size
      const Teuchos::RCP<Epetra_Vector> disn,    //!< displacement state at n+1
      Teuchos::RCP<Epetra_Vector>& fint,         //!< internal force
      Teuchos::RCP<LINALG::SparseOperator>& stiff  //! stiffness
    );

    //! evaluate certain potential forces and its stiffness
    //!
    //! Evaluation happens internal-force-like.
    //! Which potentials are actually meant can be found
    //! at DRT::PotentialManager
    void ApplyForceStiffPotential
    (
      const double time, //!< total time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force vector
      Teuchos::RCP<LINALG::SparseOperator>& stiff  //!< stiffness
    );


    //! evaluate certain potential forces and its stiffness
    //! and computes an analytical solution
    //!
    //! Evaluation happens internal-force-like.
    //! Which potentials are actually meant can be found
    //! at DRT::PotentialManager
    void TestForceStiffPotential
    (
      const double                      time,
      const Teuchos::RCP<Epetra_Vector> dis,
      const int                         step
    );

    //! Evaluate forces and stiffness due to contact / meshtying
    void ApplyForceStiffContactMeshtying
    (
      Teuchos::RCP<LINALG::SparseOperator>& stiff,  //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector>& fres,  //!< residual forces are modified
      Teuchos::RCP<Epetra_Vector>& dis,   //!< current displacement state
      bool predict                        //!< flag indicating predictor step
    );

    //! Evaluate forces and stiffness due to beam contact
    void ApplyForceStiffBeamContact
    (
      Teuchos::RCP<LINALG::SparseOperator>& stiff,  //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector>& fres,  //!< residual forces are modified
      Teuchos::RCP<Epetra_Vector>& dis,   //!< current displacement state
      bool predict                        //!< flag indicating predictor step
    );

    //! Evaluate forces and stiffness due to spring dashpot boundary condition
    void ApplyForceStiffSpringDashpot
    (
      Teuchos::RCP<LINALG::SparseOperator> stiff,  //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal forces are modified
      Teuchos::RCP<Epetra_Vector> dis,   //!< current displacement state
      Teuchos::RCP<Epetra_Vector> vel,   //!< current velocity state
      bool predict,                      //!< flag indicating predictor step
      Teuchos::ParameterList psprdash    //!< Parameterlist containing scale factors for matrix entries
    );

    //@}

    //! @name Solution
    //@{

    //! determine characteristic norms for relative
    //! error checks of residual displacements
    //! \author lw  \date 12/07
    virtual double CalcRefNormDisplacement() = 0;

    //! determine characteristic norms for relative
    //! error checks of residual forces
    //! \author lw  \date 12/07
    virtual double CalcRefNormForce() = 0;

    //! Is convergence reached of iterative solution technique?
    //! Keep your fingers crossed...
    //! \author lw  \date 12/07
    bool Converged();

    //! Solve dynamic equilibrium
    //!
    //! This is a general wrapper around the specific techniques.
    int Solve();

    //! Do full Newton-Raphson iteration
    //!
    //! This routines expects a prepared negative residual force #fres_
    //! and associated effective stiffness matrix #stiff_
    //void NewtonFull();
    int NewtonFull();

    //! check for success of linear solve otherwise return error code
    int LinSolveErrorCheck(int linerror);

    //! check for success of nonlinear solve otherwise return error code
    int NewtonFullErrorCheck(int linerror);

    //! Do (so-called) modified Newton-Raphson iteration in which
    //! the initial tangent is kept and not adapted to the current
    //! state of the displacement solution
    void NewtonModified() { dserror("Not impl."); }

    //! Do classical augmented lagrange for volume constraint
    //!
    //! Potential is linearized wrt displacements
    //! keeping Lagrange multiplier fixed.
    //! Until convergence Lagrange multiplier
    //! is increased by Uzawa_param*(Vol_err)
    //!
    //! \author tk (originally)
    int UzawaNonLinearNewtonFull();

    //! do full Newton iteration respecting volume constraint
    //!
    //! Potential is linearized wrt displacements
    //! and Lagrange multipliers
    //! Linear problem is solved with Uzawa algorithm.
    //!
    //! \author tk (originally) \date 11/07
    int UzawaLinearNewtonFull();

    //! check for success of nonlinear solve otherwise return error code
    int UzawaLinearNewtonFullErrorCheck(int linerror);

    //! Do pseudo transient continuation non-linear iteration
    //!
    //! Pseudo transient continuation is a variant of full newton that has a
    //! larger convergence radius than newton and is therefore more stable
    //! and/or can do larger time steps
    //!
    //! \author mwgee (originally) \date 03/12
    int PTC();

    //! Do nonlinear iteration for contact / meshtying
    //!
    //! \author popp (originally) \date 03/10
    int CmtNonlinearSolve();

    //! Call linear solver for contact / meshtying
    //!
    //! The point is that of course we usually just want to perform a
    //! standard solver call, except for the standard Lagrange strategy
    //! where we explicitly have to set up the saddle point system blocks
    //! and then solve this extended system.
    //!
    //! \author popp (originally) \date 03/10
    void CmtLinearSolve();

    int CmtWindkConstrNonlinearSolve();
    void CmtWindkConstrLinearSolve();

    //! Do nonlinear iteration for beam contact
    //!
    //! \author popp (originally) \date 11/11
    int BeamContactNonlinearSolve();

    //@}

    //! @name STC business
    //@{
    //! Precondition system and return operator
    void STCPreconditioning();

    //! Build up STC Matrix
    void ComputeSTCMatrix();

    //! recover standard solution
    void RecoverSTCSolution();
    //@}

    //! @name NOX solution
    //@{

    //! Compute the residual of linear momentum
    //!
    //! Computes the set of nonlinear equations, \f$F(x)=0\f$, to be solved by NOX.
    //! This method must be supplied by the user.
    bool computeF
    (
      const Epetra_Vector& x,  //!< solution vector \f$x\f$ specified from NOX
                               //!< i.e. total displacements
      Epetra_Vector& RHS,  //!< to be filled with the \f$F\f$ values that correspond
                           //!< to the input solution vector \f$x\f$.
      const NOX::Epetra::Interface::Required::FillType flag  //!< enumerated
                                        //!< type (see NOX::Epetra::FillType)
                                        //!< that tells a users interface why
                                        //!< computeF() was called. NOX has
                                        //!< the ability to generate
                                        //!< Jacobians based on numerical
                                        //!< differencing using calls to
                                        //!< computeF(). In this case,
                                        //!< the user may want to compute
                                        //!< an inexact
                                        //!< (and hopefully cheaper) \f$F\f$
                                        //!< since it
                                        //!< is only used in the Jacobian (or
                                        //!< preconditioner).
    );

    //! Compute effective dynamic stiffness matrix
    bool computeJacobian
    (
      const Epetra_Vector& x,  //!< solution vector \f$x\f$ specified from NOX
                               //!< i.e. total displacements
      Epetra_Operator& Jac  //!< a reference to the Jacobian operator
                            //!< \f$\frac{\partial F}{\partial x}\f$
                            //!< that the user supplied in the
                            //!< NOX::Epetra::Group constructor.
    );

    //! Compute preconditioner
    //!
    //! This method allows a user to supply their own preconditioner.
    //! The method should compute a preconditioner based upon the
    //! solution vector x and store it in the Epetra_Operator M.
    //! Subsequent calls to the Epetra_Operator::Apply method will apply
    //! this user supplied preconditioner to epetra vectors.
    //! The Epetra_Operator M can also be an Epetra_RowMatrix that can
    //! be used by internally constructed preconditioner objects such as AztecOO, Ifpack, and ML
    bool computePreconditioner
    (
      const Epetra_Vector& x,  //!< solution vector \f$x\f$ specified from NOX
                               //!< i.e. total displacements
      Epetra_Operator& M,  //!< a reference to the operator that is to be filled.
                           //!< This operator should be updated to correspond to the current
                           //!< solution vector x.
      Teuchos::ParameterList* precParams  //!< a pointer to a parameter list
                                          //!< that can be used by the preconditioner.
    );

    //! Setup for solution with NOX
    void NoxSetup();

    //! Setup for solution with NOX
    void NoxSetup(
      const Teuchos::ParameterList& noxparams  //!< NOX parameters from read-in
    );

    //! Create status test for non-linear solution with NOX
    Teuchos::RCP<NOX::StatusTest::Combo> NoxCreateStatusTest
    (
      Teuchos::RCP<NOX::Abstract::Group> grp  //!< NOX group
    );

    //! Create solver parameters for  non-linear solution with NOX
    Teuchos::RCP<Teuchos::ParameterList> NoxCreateSolverParameters();

    //! Create printing parameters for non-linear solution with NOX
    Teuchos::RCP<Teuchos::ParameterList> NoxCreatePrintParameters
    (
      const bool verbose = false  //!< verbosity level
    ) const;

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem> NoxCreateLinearSystem
    (
      Teuchos::ParameterList& nlParams,
      NOX::Epetra::Vector& noxSoln,
      Teuchos::RCP<NOX::Utils> utils
    );

    //! Do non-linear solve with NOX
    int NoxSolve();

    //! moved from the ADAPTER to here, to be checked carefully
//    Teuchos::RCP<Epetra_Vector> RelaxationSolve(
//      Teuchos::RCP<Epetra_Vector> iforce
//    )
//    {
//      Teuchos::RCP<Epetra_Vector> relax = interface_->InsertFSICondVector(iforce);
//      SetForceInterface(relax);
//      Teuchos::RCP<Epetra_Vector> idisi = SolveRelaxationLinear();
//
//      // we are just interested in the incremental interface displacements
//      return interface_->ExtractFSICondVector(idisi);
//    }
    // check for success of nonlinear solve otherwise return error code
    int NoxErrorCheck(NOX::StatusTest::StatusType status, Teuchos::RCP<NOX::Solver::Generic> solver);
    //@}

    //! @name Updates
    //@{

    //! Update iteration
    //!
    //! This handles the iterative update of the current
    //! displacements \f$D_{n+1}\f$ with the residual displacements
    //! The velocities and accelerations follow on par.
    void UpdateIter
    (
      const int iter  //!< iteration counter
    );

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #accn_
    //! from scratch by using the newly updated #disn_ according to the time
    //! integration scheme.
    //! The method respects the Dirichlet DOFs which are not touched.
    //! We have to use this update routine if we are not sure whether velocities
    //! and accelerations have been computed consistently to the displacements
    //! based on time integration scheme specific formulas. Hence, this method
    //! is necessary for certain predictors (like #PredictConstDisVelAcc).
    virtual void UpdateIterIncrementally() = 0;

    //! Update iteration incrementally with prescribed residual
    //! displacements
    void UpdateIterIncrementally
    (
      const Teuchos::RCP<const Epetra_Vector> disi  //!< input residual displacements
    );

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    //! We can use this update routine if we are sure that the velocities
    //! and accelerations have been computed consistently to the displacements
    //! based on time integration scheme specific formulas.
    virtual void UpdateIterIteratively() = 0;

    //@}


    //! Update configuration after time step
    //!
    //! This means, the state set
    //! \f$D_{n} := D_{n+1}\f$, \f$V_{n} := V_{n+1}\f$ and \f$A_{n} := A_{n+1}\f$
    //! Thus the 'last' converged state is lost and a reset
    //! of the time step becomes impossible.
    //! We are ready and keen awaiting the next time step.
    virtual void UpdateStepState() = 0;

    //! Update Element
    virtual void UpdateStepElement() = 0;

    //! Update step for surface stress
    void UpdateStepSurfstress();

    //! Update step for constraints
    void UpdateStepConstraint();

    //! Update step for Windkessel
    void UpdateStepWindkessel();

    //! @name Output
    //@{

    //! Print to screen predictor informations about residual norm etc.
    //! \author lw (originally) \date 12/07
    void PrintPredictor();

    //! Print to screen information about residual forces and displacements
    //! \author lw (originally) \date 12/07
    void PrintNewtonIter();

    //! Contains text to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterText
    (
      FILE* ofile  //!< output file handle
    );

    //! Contains header to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterHeader
    (
      FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void PrintNewtonConv();

    //! print summary after step
    void PrintStep();

    //! The text for summary print, see #PrintStep
    void PrintStepText
    (
      FILE* ofile  //!< output file handle
    );

    // export contact integration time and active set into textfiles
    // xxx.time and xxx.active
    void ExportContactQuantities();
    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    virtual enum INPAR::STR::DynamicType MethodName() const = 0;

    //! These time integrators are all implicit (mark their name)
    bool MethodImplicit() { return true; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() = 0;

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() = 0;

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() = 0;

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() = 0;

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() = 0;

    //! Return bool indicating if constraints are defined
    bool HaveConstraint();
    //! Return bool indicating if Windkessel bcs are defined
    bool HaveWindkessel();

    //! Return Teuchos::rcp to ConstraintManager conman_
    const Teuchos::RCP< ::UTILS::ConstrManager> GetConstraintManager()
    {
      return conman_;
    }

    //! Return Teuchos::rcp to WindkesselManager windkman_
    const Teuchos::RCP< ::UTILS::WindkesselManager> GetWindkesselManager()
    {
      return windkman_;
    }

    //! Get type of thickness scaling for thin shell structures
    INPAR::STR::STC_Scale GetSTCAlgo()
    {
      return stcscale_;
    }

    //! Access to scaling matrix for STC
    Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    {
      if (!stccompl_)
        ComputeSTCMatrix();
      stccompl_ = true;
      return stcmat_;
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr
    (
      Teuchos::RCP<Epetra_Vector> lagrincr ///< Lagrange multiplier increment
    );

    //! Update iteration
    //! Add residual increment to pressures stored in Windkessel manager
    virtual void UpdateIterIncrWindkessel
    (
      Teuchos::RCP<Epetra_Vector> presincr ///< pressure increment
    );

    //@}

    //! @name Fluid-structure-interaction specific methods
    //@{

    //! switch structure field to block matrix in fsi simulations
    void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    virtual void EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params) = 0;

    //! Linear structure solve with just an interface load
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear();

    //! update all EpetraVectors as discretization changes after crack insertion
    void updateEpetraVectorsCrack( std::map<int,int>& oldnew );

	  /// update all the field vectors to account for new node introduced by crack
    virtual void updateMethodSpecificEpetraCrack( std::map<int,int>& oldnew ) = 0;

        //! Prepare system for solving with Newton's method
    //!
    //! - negative residual
    //! - blank residual on Dirichlet DOFs
    //! - apply Dirichlet boundary conditions on system
    void PrepareSystemForNewtonSolve();

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return external force \f$F_{ext,n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> FextNew() = 0;

    //! Return reaction forces
    //!
    //! This is a vector of length holding zeros at
    //! free DOFs and reaction force component at DOFs on DBCs.
    //! Mark, this is not true for DBCs with local co-ordinate
    //! systems in which the non-global reaction force
    //! component is stored in global Cartesian components.
    //! The reaction force resultant is not affected by
    //! this operation.
    Teuchos::RCP<Epetra_Vector> Freact() { return freact_; }

    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;

    //! Write internal and external forces for restart
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! FORMERLY: Return residual displacements \f$\Delta D_{n+1}^{<k>}\f$
    //! Called from the previous adapters as InitialGuess()
    Teuchos::RCP<const Epetra_Vector> InitialGuess() { return disi_; }

    //! Prepare time step for coupled problems with monolithic approach
    void PrepareTimeStep()
    {
      // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
      SetTimen( (*time_)[0]+(*dt_)[0] );

      // prepare contact for new time step
      PrepareStepContact();

      // predict
      Predict();
    }

    //! Evaluate routine for coupled problems with monolithic approach
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> disiterinc
    )
    {
      UpdateIterIncrementally(disiterinc);

      Teuchos::ParameterList params;

      // builds tangent, residual and applies DBC
      EvaluateForceStiffResidual(params);
      PrepareSystemForNewtonSolve();
      return;
    }

    //! Update routine for coupled problems with monolithic approach
    virtual void Update()
    {
      UpdateStepState();
      UpdateStepTime();
      UpdateStepElement();
      return;
    }

    //! Update routine for coupled problems with monolithic approach with time adaptivity
    virtual void Update(const double endtime)
    {
      UpdateStepState();

      timen_ = endtime;

      //Update
      time_->UpdateSteps(timen_);
      step_ = stepn_;
      stepn_ += 1;

      UpdateStepElement();
      return;
    }



    //! Output routine for coupled problems with monolithic approach
    void Output(bool forced_writerestart = false)
    {
      OutputStep(forced_writerestart);
    }

    //! Set residual displacements \f$\Delta D_{n+1}^{<k>}\f$
    void SetDisResidual
    (
      const Teuchos::RCP<const Epetra_Vector> disi  //!< input residual displacements
    )
    {
      if (disi != Teuchos::null) disi_->Update(1.0, *disi, 0.0);
    }

    //! Return the rhs-vector (negative sign for Newton is already included.)
    Teuchos::RCP<const Epetra_Vector> RHS() { return fres_; }

    virtual Teuchos::RCP<Epetra_Vector> Velaf()
    {
      dserror("mid velocities not available");
      return Teuchos::null;
    }

    //@}

  protected:

    //! copy constructor is NOT wanted
    TimIntImpl(const TimIntImpl& old);

    //! @name General purpose algorithm parameters
    //@{
    enum INPAR::STR::PredEnum pred_;  //!< predictor
    //@}

    //! @name Iterative solution technique
    //@{
    enum INPAR::STR::NonlinSolTech itertype_;  //!< kind of iteration technique
                                               //!< or non-linear solution technique
    enum INPAR::STR::ConvNorm normtypedisi_; //!< convergence check for residual displacements
    enum INPAR::STR::ConvNorm normtypefres_; //!< convergence check for residual forces
    enum INPAR::STR::ConvNorm normtypepres_; //!< convergence check for residual pressure
    enum INPAR::STR::ConvNorm normtypepfres_; //!< convergence check for residual pressure forces
    enum INPAR::STR::ConvNorm normtypecontconstr_; //!< convergence check for contact constraints (saddlepoint formulation only)
    enum INPAR::STR::ConvNorm normtypeplagrincr_;  //!< convergence check for Lagrange multiplier increment (saddlepoint formulation only)

    enum INPAR::STR::BinaryOp combdispre_; //!< binary operator to combine field norms
    enum INPAR::STR::BinaryOp combfrespfres_; //!< binary operator to combine field norms
    enum INPAR::STR::BinaryOp combdisifres_; //!< binary operator to combine displacement and forces
    enum INPAR::STR::BinaryOp combfrescontconstr_; //!< binary operator to combine field norms (forces and contact constraints, contact/meshtying in saddlepoint formulation only)
    enum INPAR::STR::BinaryOp combdisilagr_; //!< binary operator to combine field norms (displacement increments and LM increments, contact/meshtying in saddlepoint formulation only)
    enum INPAR::STR::BinaryOp combfresplconstr_; //!< binary operator to combine field norms (forces and plastic constraints, semi-smooth plasticity only)
    enum INPAR::STR::BinaryOp combdisiLp_; //!< binary operator to combine field norms (displacement increments and Lp increments, semi-smooth plasticity only)

    enum INPAR::STR::VectorNorm iternorm_;  //!< vector norm to check with
    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iterations

    enum INPAR::STR::DivContAct divcontype_ ; // what to do when nonlinear solution fails
    double toldisi_;  //!< tolerance residual displacements
    double tolfres_;  //!< tolerance force residual
    double tolpfres_; //!< norm of residual pressure forces
    double tolpres_; //!< norm of residual pressures
    double tolcontconstr_; //!< norm of rhs for contact constraints (saddlepoint formulation only)
    double tollagr_; //!< tolerance of LM multiplier increments (saddlepoint formulation only)

    double uzawaparam_;  //!< Parameter for Uzawa algorithm dealing
                         //!< with Lagrange multipliers
    int uzawaitermax_;  //!< maximally permitted Uzawa iterations
    double tolcon_;  //!< tolerance constraint
    double tolwindk_;  //!< tolerance Windkessel
    int iter_;  //!< iteration step
    double normcharforce_;  //!< characteristic norm for residual force
    double normchardis_;  //!< characteristic norm for residual displacements
    double normfres_;  //!< norm of residual forces
    double normdisi_;  //!< norm of residual displacements
    double normcon_;  //!< norm of constraint
    double normwindk_;  //!< norm of Windkessel
    double normpfres_; //!< norm of residual pressure forces
    double normpres_; //!< norm of residual pressures
    double normcontconstr_;  //!< norm of contact/meshtying constraint rhs (contact/meshtying in saddlepoint formulation only)
    double normlagr_;     //!< norm of lagrange multipliers
    double normw_;     //!< norm of wear
    double normwrhs_;
    double normwm_;     //!< norm of wear master
    double normwmrhs_;

    Teuchos::RCP<Epetra_Vector> disi_;  //!< residual displacements (and pressure)
                                        //!< \f$\Delta{D}^{<k>}_{n+1}\f$
    //! @name Various global forces
    //@{
    Teuchos::RCP<Epetra_Vector> fres_;  //!< force residual used for solution
    Teuchos::RCP<Epetra_Vector> freact_;  //!< reaction force
    //@}

    //! @name NOX variables
    //@{
    Teuchos::RCP<NOX::StatusTest::Combo> noxstatustest_;
    Teuchos::RCP<Teuchos::ParameterList> noxparams_;
    Teuchos::RCP<NOX::Utils> noxutils_;
    //@}

    //! @name Krylov projection variables
    bool updateprojection_;                             //!< bool triggering update of Krylov projection
    Teuchos::RCP<LINALG::KrylovProjector>   projector_; //!< Krylov projector himself
    //@}

    //! @name SDC Scaling for thin shell structures
    //@{
    enum INPAR::STR::STC_Scale stcscale_; //!< scale thickness of shells?
    double stcfact_;  //!< scaling factor for STC
    int stclayer_;  //! number of layers for multilayered case
    Teuchos::RCP<LINALG::SparseMatrix> stcmat_;  //!< scaling matrix for STC
    bool stccompl_; ///< keep track if stc matrix has been evaluated
    //@}

    //! @name Pseudo Transient Continuation Parameters
    //@{
    double ptcdt_; //!< pseudo time step size for PTC
    double dti_;   //!< scaling factor for PTC (initially 1/ptcdt_, then adapted)
    //@}

  };  // class TimIntImpl

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_IMPL_H
