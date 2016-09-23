/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.H
\brief Implicit time integration for structural dynamics

<pre>
\level 1

\maintainer Alexander Popp
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
#include "../drt_lib/drt_discret.H"

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
   * The implicit time integrator object is a derivation of the base time
   * integrators with an eye towards implicit time integration. #TimIntImpl
   * provides the environment needed to execute implicit integrators. This is
   * chiefly the non-linear solution technique, e.g., Newton-Raphson iteration.
   * These iterative solution techniques require a set of control parameters
   * which are stored within this object. It is up to derived object to
   * implement the time-space discretised residuum and its tangent. This object
   * provides some utility functions to obtain various force vectors necessary
   * in the calculation of the force residual in the derived time integrators.
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
      const Teuchos::ParameterList& timeparams,
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

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    virtual void Init(
        const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<LINALG::Solver> solver );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();


    //! Resize \p TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //! return time integration factor
    virtual double TimIntParam() const = 0;

    //@}

    //! Do time integration of single step
    int IntegrateStep();

    //! Create Edges of for discrete shell elements
    void InitializeEdgeElements();

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

    //! Check if line search is applied in combination with elements
    //! that perform a local condensation (e.g. EAS)
    void PrepareLineSearch();

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

    /*! @name Forces
     *
     *  Apply all sets of forces (external, internal, damping, inertia, ...)
     *  and corresponding stiffnesses based on the current solution state.
     *
     *  On the level of STR::TimIntImpl, we deal with forces and their stiffness
     *  contributions since an implicit time integration in BACI often requires
     *  full linearization.
     *
     *  For application of forces only (without stiffness), see STR::TimInt.
     *
     *  \sa STR::TimInt
     */
    //@{

    /*! \brief Evaluate forces, stiffness #stiff_ and residual #fres_
     *
     *  Do residual force due to global balance of momentum
     *  and its tangent with respect to the current
     *  displacements \f$D_{n+1}\f$
     *
     *  This is <i>the</i> central method which is different for each
     *  derived implicit time integrator. The time integrator implementation
     *  is expected to set members #fres_ and #stiff_.
     *  The residual #fres_ is expected to follow the <i>same</i> sign
     *  convention like its tangent #stiff_, i.e. to use
     *  Newton--Raphson's method the residual will be scaled by -1.
     *
     *  A parameter list is used to pass on additional information
     *  from from the non-linear solution procedure. For instance,
     *  in structural contact simulations and semi-smooth newton plasticity,
     *  it is of importance whether this method is called from a predictor step
     *  or from a regular Newton step.
     */
    virtual void EvaluateForceStiffResidual(Teuchos::ParameterList& params) = 0;

    /*! \brief Evaluate forces and residual #fres_
     *
     *  Do residual force due to global balance of momentum.
     *
     *  This is <i>the</i> central method which is different for each
     *  derived implicit time integrator. The time integrator implementation
     *  is expected to set members #fres_.
     *
     *  For now, this only provides basic functionalities, i.e. pure solid
     *  dynamics without contact, plasticity, constraints, ...
     *
     *  \f[ f_{res} = f_{int} - f_{ext} \f]
     *
     *  \author mayr.mt \date 09/2014
     */
    virtual void EvaluateForceResidual() = 0;

    //! Apply external force, its stiffness at state
    void ApplyForceStiffExternal
    (
      const double time, //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis, //!< old displacement state
      const Teuchos::RCP<Epetra_Vector> disn, //!< new displacement state
      const Teuchos::RCP<Epetra_Vector> vel, // velocity state
      Teuchos::RCP<Epetra_Vector>& fext, //!< external force
      Teuchos::RCP<LINALG::SparseOperator>& fextlin //!<linearization of external force
    );

    //! Evaluate ordinary internal force, its stiffness at state
    void ApplyForceStiffInternal
    (
      const double time, //!< evaluation time
      const double dt, //!< step size
      const Teuchos::RCP<Epetra_Vector> dis, //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi, //!< residual displacements
      const Teuchos::RCP<Epetra_Vector> vel, // velocity state
      Teuchos::RCP<Epetra_Vector> fint, //!< internal force
      Teuchos::RCP<LINALG::SparseOperator> stiff, //!< stiffness matrix
      Teuchos::ParameterList& params, //!< parameters from nonlinear solver
      Teuchos::RCP<LINALG::SparseOperator> damp = Teuchos::null //!< material damping matrix
    );

    //! Evaluate internal and inertia forces and their linearizations
    void ApplyForceStiffInternalAndInertial
    (
      const double time, //!< evaluation time
      const double dt, //!< step size
      const double timintfac_dis, //!< time integration factor for d(Res) / d dis
      const double timintfac_vel, //!< time integration factor for d(Res) / d vel
      const Teuchos::RCP<Epetra_Vector> dis, //!< displacement state
      const Teuchos::RCP<Epetra_Vector> disi, //!< residual displacements
      const Teuchos::RCP<Epetra_Vector> vel, //! velocity state
      const Teuchos::RCP<Epetra_Vector> acc, //! acceleration state
      Teuchos::RCP<Epetra_Vector> fint, //!< internal force
      Teuchos::RCP<Epetra_Vector> finert, //!< inertia force
      Teuchos::RCP<LINALG::SparseOperator> stiff, //!< stiffness matrix
      Teuchos::RCP<LINALG::SparseOperator> mass, //!< mass matrix
      Teuchos::ParameterList& params, //!< parameters from nonlinear solver
      const double beta=1000.0, //!< time integration parameters for element-wise time integration (necessary for time integration of rotations)
      const double gamma=1000.0, //!< time integration parameters for element-wise time integration (necessary for time integration of rotations)
      const double alphaf=1000.0, //!< time integration parameters for element-wise time integration (necessary for time integration of rotations)
      const double alpham=1000.0 //!< time integration parameters for element-wise time integration (necessary for time integration of rotations)
    );

    //! Evaluate forces due to constraints
    void ApplyForceStiffConstraint
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis, //!< last evaluated displacements
      const Teuchos::RCP<Epetra_Vector> disn, //!< evaluation displacements
      Teuchos::RCP<Epetra_Vector>& fint, //!< forces are added onto
      Teuchos::RCP<LINALG::SparseOperator>& stiff, //!< stiffness is added onto
      Teuchos::ParameterList pcon  //!< parameter list containing scale factors for matrix entries
    );

    //! Evaluate forces due to Cardiovascular0D BCs
    void ApplyForceStiffCardiovascular0D
    (
      const double time, //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> disn, //!< evaluation displacements
      Teuchos::RCP<Epetra_Vector>& fint, //!< forces are added onto
      Teuchos::RCP<LINALG::SparseOperator>& stiff, //!< stiffness is added onto
      Teuchos::ParameterList pcon //!< parameter list containing scale factors for matrix entries
    );

    //! Evaluate _certain_ surface stresses and stiffness
    //!
    //! The Evaluation happens internal-force-like.
    void ApplyForceStiffSurfstress
    (
      const double time, //!< evaluation time
      const double dt, //!< time step size
      const Teuchos::RCP<Epetra_Vector> disn, //!< displacement state at n+1
      Teuchos::RCP<Epetra_Vector>& fint, //!< internal force
      Teuchos::RCP<LINALG::SparseOperator>& stiff //! stiffness
    );

    //! evaluate certain potential forces and its stiffness
    //!
    //! Evaluation happens internal-force-like.
    //! Which potentials are actually meant can be found
    //! at DRT::PotentialManager
    void ApplyForceStiffPotential
    (
      const double time, //!< total time
      const Teuchos::RCP<Epetra_Vector> dis, //!< displacement state
      Teuchos::RCP<Epetra_Vector>& fint, //!< internal force vector
      Teuchos::RCP<LINALG::SparseOperator>& stiff //!< stiffness
    );

    //! evaluate certain potential forces and its stiffness
    //! and computes an analytical solution
    //!
    //! Evaluation happens internal-force-like.
    //! Which potentials are actually meant can be found
    //! at DRT::PotentialManager
    void TestForceStiffPotential
    (
      const double time,
      const Teuchos::RCP<Epetra_Vector> dis,
      const int step
    );

    //! Evaluate forces and stiffness due to contact / meshtying
    void ApplyForceStiffContactMeshtying
    (
      Teuchos::RCP<LINALG::SparseOperator>& stiff, //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector>& fres, //!< residual forces are modified
      Teuchos::RCP<Epetra_Vector>& dis, //!< current displacement state
      bool predict //!< flag indicating predictor step
    );

    //! Evaluate forces and stiffness due to beam contact
    void ApplyForceStiffBeamContact
    (
      Teuchos::RCP<LINALG::SparseOperator>& stiff, //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector>& fres, //!< residual forces are modified
      Teuchos::RCP<Epetra_Vector>& dis, //!< current displacement state
      bool predict //!< flag indicating predictor step
    );

    //! Check residual displacement and scale it if necessary
    void LimitStepsizeBeamContact
    (
      Teuchos::RCP<Epetra_Vector>& disi //!< residual displacement vector
    );

    //! Check residual displacement and scale it if necessary
    void LimitStepsizeBeam
    (
      Teuchos::RCP<Epetra_Vector>& disi //!< residual displacement vector
    );

    //! Evaluate forces and stiffness due to spring dash-pot boundary condition
    void ApplyForceStiffSpringDashpot
    (
      Teuchos::RCP<LINALG::SparseOperator> stiff, //!< stiffness is modified
      Teuchos::RCP<Epetra_Vector> fint, //!< internal forces are modified
      Teuchos::RCP<Epetra_Vector> dis, //!< current displacement state
      Teuchos::RCP<Epetra_Vector> vel, //!< current velocity state
      bool predict, //!< flag indicating predictor step
      Teuchos::ParameterList psprdash //!< parameter list containing scale factors for matrix entries
    );

    //@}

    //! @name Solution
    //@{

    //! determine characteristic norms for relative
    //! error checks of residual displacements
    //! \author lw  \date 12/07
    virtual double CalcRefNormDisplacement();

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
    //! Do the nonlinear solve, i.e. (multiple) corrector,
    //! for the time step. All boundary conditions have
    //! been set.
    INPAR::STR::ConvergenceStatus Solve();

    //! Do full Newton-Raphson iteration
    //!
    //! This routines expects a prepared negative residual force #fres_
    //! and associated effective stiffness matrix #stiff_
    //void NewtonFull();
    int NewtonFull();

    //! check for success of element evaluation in that no negative Jacobian
    //! determinant occured, otherwise return error code
    int ElementErrorCheck(bool evalerr);

    //! check for success of linear solve otherwise return error code
    int LinSolveErrorCheck(int linerror);

    //! check for success of nonlinear solve otherwise return error code
    int NewtonFullErrorCheck(int linerror, int eleerror);

    //! Do (so-called) modified Newton-Raphson iteration in which
    //! the initial tangent is kept and not adapted to the current
    //! state of the displacement solution
    void NewtonModified() { dserror("Not impl."); }

    //! Do Line search iteration
    //!
    //! This routines expects a prepared negative residual force #fres_
    //! and associated effective stiffness matrix #stiff_
    int NewtonLS();

    //! Solver call (line search)
    int LsSolveNewtonStep();

    //! Update structural RHS and stiffness matrix (line search)
    void LsUpdateStructuralRHSandStiff(bool& isexcept,double& merit_fct);

    //! Evaluate the specified merit function (line search)
    //! (for pure structural problems this routine is rather short.
    //! However, to incorporate contact problems more easily the
    //! evaluation of the merit function is performed here.)
    //! return 0 if successful
    int LsEvalMeritFct(double& merit_fct);

    //! Check the inner linesearch loop for convergence (line search)
    bool LsConverged(double* mf_value, double step_red);

    //! Print information concerning the last line search step (line search)
    void LsPrintLineSearchIter(double* mf_value, int iter_ls, double step_red);

    //! Contains text to LsPrintLineSearchIter
    void LsPrintLsIterText
    (
      FILE* ofile  //!< output file handle
    );

    //! Contains header to LsPrintLineSearchIter
    void LsPrintLsIterHeader
    (
      FILE* ofile  //!< output file handle
    );

    /*! \brief Use nonlinear solver package NLNSOL
     *
     *  \author mayr.mt \date 06/2014
     */
    int NlnSolver();

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
    int UzawaLinearNewtonFullErrorCheck(int linerror, int eleerror);

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
    int CmtWindkConstrLinearSolve();

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

    /*! \brief Compute the residual of linear momentum
     *
     *  Computes the set of nonlinear equations, \f$F(x)=0\f$, to be solved by
     *  NOX. This method must be supplied by the user.
     */
    bool computeF(
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

    /*! \brief Compute preconditioner
     *
     *  This method allows a user to supply their own preconditioner.
     *  The method should compute a preconditioner based upon the solution
     *  vector x and store it in the Epetra_Operator M. Subsequent calls to the
     *  Epetra_Operator::Apply method will apply this user supplied
     *  preconditioner to epetra vectors. The Epetra_Operator M can also be an
     *  Epetra_RowMatrix that can be used by internally constructed
     *  preconditioner objects such as AztecOO, Ifpack, and ML.
     */
    bool computePreconditioner(
      const Epetra_Vector& x, //!< solution vector \f$x\f$ specified from NOX
                              //!< i.e. total displacements
      Epetra_Operator& M, //!< a reference to the operator that is to be filled.
                          //!< This operator should be updated to correspond to the current
                          //!< solution vector x.
      Teuchos::ParameterList* precParams //!< a pointer to a parameter list
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
    Teuchos::RCP<Teuchos::ParameterList> NoxCreatePrintParameters(
      const bool verbose = false //!< verbosity level
      ) const;

    //! Create the linear system that is passed into NOX
    Teuchos::RCP<NOX::Epetra::LinearSystem> NoxCreateLinearSystem(
      Teuchos::ParameterList& nlParams, ///< parameter list
      NOX::Epetra::Vector& noxSoln, ///< solution vector to operate on
      Teuchos::RCP<NOX::Utils> utils ///< printing utilities
      );

    //! Do non-linear solve with NOX
    int NoxSolve();

    //! check for success of nonlinear solve otherwise return error code
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

    //! Update step for Cardiovascular0D
    void UpdateStepCardiovascular0D();

    //! Update step for SpringDashpot
    void UpdateStepSpringDashpot();

    //@}

    // export contact integration time and active set into text files
    // xxx.time and xxx.active
    void ExportContactQuantities();

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    virtual enum INPAR::STR::DynamicType MethodName() const = 0;

    //! These time integrators are all implicit (mark their name)
    bool MethodImplicit() { return true; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() const = 0;

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const = 0;

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const = 0;

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const = 0;

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const = 0;

    //! Return bool indicating if constraints are defined
    bool HaveConstraint();

    //! Return bool indicating if Cardiovascular0D bcs are defined
    bool HaveCardiovascular0D();

    //! Return bool indicating if spring dashpot BCs are defined
    bool HaveSpringDashpot();

    //! Return Teuchos::rcp to ConstraintManager conman_
    const Teuchos::RCP< ::UTILS::ConstrManager> GetConstraintManager()
    {
      return conman_;
    }

    //! Return Teuchos::rcp to Cardiovascular0DManager cardvasc0dman_
    const Teuchos::RCP< ::UTILS::Cardiovascular0DManager> GetCardiovascular0DManager()
    {
      return cardvasc0dman_;
    }

    //! Return Teuchos::rcp to SpringDashpotManager springman_
    const Teuchos::RCP< ::UTILS::SpringDashpotManager> GetSpringDashpotManager()
    {
      return springman_;
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
    //! Add residual increment to pressures stored in Cardiovascular0D manager
    virtual void UpdateIterIncrCardiovascular0D
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

    //! check, if according to divercont flag time step size can be increased
    void CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status);

    /// update all the field vectors to account for new node introduced by crack
    virtual void updateMethodSpecificEpetraCrack( std::map<int,int>& oldnew ) = 0;

    /*! \brief Prepare system for solving with Newton's method
     *
     *  Blank DOFs with Dirichlet BCs in the residual. By default
     *  (preparejacobian = true), apply Dirichlet BCs to #stiff_ as well. This
     *  can be switched off when only the residual has been evaluated
     *  (\sa EvaluateForceResidual()).
     */
    void PrepareSystemForNewtonSolve(
        const bool preparejacobian = true ///< prepare Jacobian matrix
        );

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return external force \f$F_{ext,n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> FextNew() = 0;

    /*! \brief Return reaction forces
     *
     *  This is a vector of length holding zeros at free DOFs and reaction force
     *  component at DOFs on DBCs.
     *
     *  \note This is not true for DBCs with local coordinate systems in which
     *  the non-global reaction force  component is stored in global Cartesian
     *  components. The reaction force resultant is not affected by this
     *  operation.
     */
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
      // safety checks
      CheckIsInit();
      CheckIsSetup();

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
      PreUpdate();
      UpdateStepState();
      UpdateStepTime();
      UpdateStepElement();
      PostUpdate();
      return;
    }

    //! Update routine for coupled problems with monolithic approach with time adaptivity
    virtual void Update(const double endtime)
    {
      PreUpdate();
      UpdateStepState();

      timen_ = endtime;

      //Update
      time_->UpdateSteps(timen_);
      step_ = stepn_;
      stepn_ += 1;

      UpdateStepElement();
      PostUpdate();
      return;
    }

    //! Output routine for coupled problems with monolithic approach
    void Output(bool forced_writerestart = false)
    {
      PreOutput();
      OutputStep(forced_writerestart);
      // write Gmsh output
      writeGmshStrucOutputStep();
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

    //@}

  protected:

    //! @name Output to file or screen
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

    //@}

    //! copy constructor is NOT wanted
    TimIntImpl(const TimIntImpl& old);

    //! @name General purpose algorithm parameters
    //@{
    enum INPAR::STR::PredEnum pred_;  //!< predictor
    //@}

    //! @name Iterative solution technique
    //@{
    enum INPAR::STR::NonlinSolTech itertype_;       //!< kind of iteration technique or non-linear solution technique

    enum INPAR::STR::ConvNorm normtypedisi_;        //!< convergence check for residual displacements
    enum INPAR::STR::ConvNorm normtypefres_;        //!< convergence check for residual forces
    enum INPAR::STR::ConvNorm normtypepres_;        //!< convergence check for residual pressure
    enum INPAR::STR::ConvNorm normtypepfres_;       //!< convergence check for residual pressure forces
    enum INPAR::STR::ConvNorm normtypecontconstr_;  //!< convergence check for contact constraints (saddlepoint formulation only)
    enum INPAR::STR::ConvNorm normtypeplagrincr_;   //!< convergence check for Lagrange multiplier increment (saddlepoint formulation only)
    enum INPAR::STR::BinaryOp combfresplconstr_; //!< binary operator to combine field norms (forces and plastic constraints, semi-smooth plasticity only)
    enum INPAR::STR::BinaryOp combdisiLp_; //!< binary operator to combine field norms (displacement increments and Lp increments, semi-smooth plasticity only)
    enum INPAR::STR::BinaryOp combfresEasres_; //!< binary operator to combine field norms (forces and EAS residuals, semi-smooth plasticity only)
    enum INPAR::STR::BinaryOp combdisiEasIncr_; //!< binary operator to combine field norms (displacement increments and EAS increments, semi-smooth plasticity only)

    enum INPAR::STR::BinaryOp combdispre_;          //!< binary operator to combine field norms
    enum INPAR::STR::BinaryOp combfrespfres_;       //!< binary operator to combine field norms
    enum INPAR::STR::BinaryOp combdisifres_;        //!< binary operator to combine displacement and forces
    enum INPAR::STR::BinaryOp combfrescontconstr_;  //!< binary operator to combine field norms (forces and contact constraints, contact/meshtying in saddlepoint formulation only)
    enum INPAR::STR::BinaryOp combdisilagr_;        //!< binary operator to combine field norms (displacement increments and LM increments, contact/meshtying in saddlepoint formulation only)

    enum INPAR::STR::VectorNorm iternorm_;  //!< vector norm to check with
    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iterations

    double toldisi_;        //!< tolerance residual displacements
    double tolfres_;        //!< tolerance force residual
    double tolpfres_;       //!< norm of residual pressure forces
    double tolpres_;        //!< norm of residual pressures
    double tolcontconstr_;  //!< norm of rhs for contact constraints (saddlepoint formulation only)
    double tollagr_;        //!< tolerance of LM multiplier increments (saddlepoint formulation only)

    double uzawaparam_;  //!< Parameter for Uzawa algorithm dealing
                         //!< with Lagrange multipliers
    int uzawaitermax_;  //!< maximally permitted Uzawa iterations
    double tolcon_;  //!< tolerance constraint
    double tolcardvasc0d_;  //!< tolerance for 0D cardiovascular residual
    double tolcardvasc0ddofincr_;  //!< tolerance for 0D cardiovascular dof increment
    int iter_;  //!< iteration step
    double normcharforce_;  //!< characteristic norm for residual force
    double normchardis_;  //!< characteristic norm for residual displacements
    double normfres_;  //!< norm of residual forces
    double normdisi_;  //!< norm of residual displacements
    double normcon_;  //!< norm of constraint
    double normcardvasc0d_;  //!< norm of 0D cardiovascular residual
    double normcardvasc0ddofincr_;  //!< norm of 0D cardiovascular dof increment
    double normpfres_; //!< norm of residual pressure forces
    double normpres_; //!< norm of residual pressures
    double normcontconstr_;  //!< norm of contact/meshtying constraint rhs (contact/meshtying in saddlepoint formulation only)
    double normlagr_;     //!< norm of lagrange multipliers
    double normw_;     //!< norm of wear
    double normwrhs_;
    double normwm_;     //!< norm of wear master
    double normwmrhs_;
    double alpha_ls_; //!< line search step reduction
    double sigma_ls_; //!< line search sufficient descent factor
    double ls_maxiter_; //!< maximum number of line search steps
    double cond_res_; //!< residual norm of condensed variables (e.g. EAS) needed for line search




    Teuchos::RCP<Epetra_Vector> disi_;  //!< residual displacements (and pressure) \f$\Delta{D}^{<k>}_{n+1}\f$

    //@}

    //! @name Various global forces
    //@{
    Teuchos::RCP<Epetra_Vector> fres_;   //!< force residual used for solution
    Teuchos::RCP<Epetra_Vector> freact_; //!< reaction force
    //@}

    //! @name NOX variables
    //@{
    Teuchos::RCP<NOX::StatusTest::Combo> noxstatustest_;  //!< NOX status test for convergence check
    Teuchos::RCP<Teuchos::ParameterList> noxparams_;      //!< NOX parameter list to configure the NOX solver
    Teuchos::RCP<NOX::Utils> noxutils_;                   //!< NOX utils for printing
    //@}

    //! @name Krylov projection variables
    bool updateprojection_;                           //!< bool triggering update of Krylov projection
    Teuchos::RCP<LINALG::KrylovProjector> projector_; //!< Krylov projector himself
    //@}

    //! @name STC Scaling for thin shell structures
    //@{
    enum INPAR::STR::STC_Scale stcscale_;       //!< scale thickness of shells?
    double stcfact_;                            //!< scaling factor for STC
    int stclayer_;                              //! number of layers for multilayered case
    Teuchos::RCP<LINALG::SparseMatrix> stcmat_; //!< scaling matrix for STC
    bool stccompl_;                             //!< keep track if stc matrix has been evaluated
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
