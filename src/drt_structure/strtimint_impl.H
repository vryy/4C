/*----------------------------------------------------------------------*/
/*!
\file strtimint_impl.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_IMPL_H
#define STRTIMINT_IMPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>
#include <NOX_Epetra_Group.H>

#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics
   *        with \b implicit time integration
   *
   * <h3> About </h3>
   * The implicit time integrator object is a derivation of the base time integrators with an eye towards implicit time integration. #TimIntImpl provides the environment needed to execute implicit integrators. This is chiefly the non-linear solution technique, e.g., Newton-Raphson iteration. These iterative solution techniques require a set of control parameters which are stored within this object. It is up to derived object to implement the time-space discretised residuum an d its tangent. This object provides some utility functions to obtain various force vectors necessary in the calculation of the force residual in the derived tim integrators.
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
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntImpl() { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //@}

    //! Do time integration of single step
    void IntegrateStep();

    //! @name Prediction
    //@{

    //! Predict target solution and identify residual
    void Predict();

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

    //@}

    //! @name Forces and stiffnesses
    //@{

    //! Do residual force due to global balance of momentum
    //! and its tangent with respect to the current 
    //! displacements \f$D_{n+1}\f$
    //!
    //! The residual #fres_ is determined readily negative
    virtual void EvaluateForceStiffResidual() = 0;

    //! Evaluate forces due to constraints
    void ApplyForceStiffConstraint
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< last evaluated displacements
      const Teuchos::RCP<Epetra_Vector> disn,  //!< evaluation displacements
      Teuchos::RCP<Epetra_Vector>& fint,  //!< forces are added onto
      Teuchos::RCP<LINALG::SparseMatrix>& stiff,  //!< stiffness is added onto
      Teuchos::ParameterList pcon   //!< Parameterlist containing scale factors for matrix entries
    );

    //! Evaluate _certain_ surface stresses and stiffness
    //!
    //! The Evaluation happens internal-force-like.
    void ApplyForceStiffSurfstress
    (
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force
      Teuchos::RCP<LINALG::SparseMatrix>& stiff  //! stiffness 
    );

    //! evaluate certain potential forces and its stiffness
    //!
    //! Evaluation happens internal-force-like.
    //! Which potentials are actually meant can be found 
    //! at DRT::PotentialManager
    void ApplyForceStiffPotential
    (
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force vector
      Teuchos::RCP<LINALG::SparseMatrix>& stiff  //!< stiffness
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
    void Solve();

    //! Do full Newton-Raphson iteration
    //! 
    //! This routines expects a prepared negative reisdual force #fres_
    //! and associated effective stiffness matrix #stiff_
    void NewtonFull();

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
    //!
    //! DEAR SIRS AND MADAMS, PLEASE PROVIDE REFERENCE
    void UzawaNonLinearNewtonFull();

    //! do full Newton iteration respecting volume constraint
    //!
    //! Potential is linearized wrt displacements 
    //! and Lagrange multipliers
    //! Linear problem is solved with Uzawa algorithm.
    //!
    //! \author tk (originally) \date 11/07
    //!
    //! DEAR SIRS AND MADAMS, PLEASE PROVIDE REFERENCE
    void UzawaLinearNewtonFull();

    //! Do non-linear CG iteration using NOX
    //!
    //! \author mwgee (originally) \date 03/07
    void NonlinearCG() { dserror("Not impl."); }

    //! Do pseudo transient continuation non-linear iteration
    //! 
    //! Pseudo transient continuation is a variant of full newton that has a
    //! larger convergence radius than newton and is therefore more stable 
    //! and/or can do larger time steps
    //!
    //! \author mwgee (originally) \date 03/07
    void PTC() { dserror("Not impl."); }

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
    //! be used by internally constructed preconditioner objects such as AztecOO, Ipfack, and ML
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
    void NoxSolve();



    //@}

    //! @name Updates
    //@{

    //! Update iteration
    //!
    //! This handles the iterative update of the current 
    //! displacements \f$D_{n+1}\f$ with the residual dispplacements
    //! The velocities and accelerations follow on par.
    void UpdateIter
    (
      const int iter  //!< iteration counter
    );

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #accn_
    //! from scratch by using the newly updated #disn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors 
    //! (like #PredictConstDisConsistVelAcc)
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
    virtual void UpdateIterIteratively() = 0;

    //@}


    //! Update configuration after time step
    //!
    //! This means, the state set
    //! \f$D_{n} := D_{n+1}\f$, \f$V_{n} := V_{n+1}\f$ and \f$A_{n} := A_{n+1}\f$
    //! Thus the 'last' converged state is lost and a reset 
    //! of the time step becomes impossible.
    //! We are ready and keen awating the next time step.
    virtual void UpdateStepState() = 0;

    //! Update step for surface stress
    void UpdateStepSurfstress()
    {
      if (surfstressman_ != Teuchos::null) surfstressman_->Update();
    }

    //! Update step for constraints
    void UpdateStepConstraint()
    {
      if (conman_ -> HaveConstraint()) conman_->Update();
    }
    
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

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    virtual enum NameEnum MethodName() = 0;

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
    virtual bool HaveConstraint()
    {
      return conman_->HaveConstraint();
    }
    
    //! Return rcp to ConstraintManager conman_
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager()
    {
      return conman_;
    }
    
    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr
    (
        Teuchos::RCP<Epetra_Vector> lagrincr ///< Lagrange multiplier increment
    )
    {
      conman_->UpdateLagrMult(lagrincr);
    }
    
    
    //@}

    //! @name Fluid-structure-interaction specific methods
    //@{

    //! introduce (robin) fsi surface extractor object
    //!
    //! Introduce the FSI surface to the fluid algorithm. This way special case
    //! algorithms (robin coupling) can work on the surface.
    void SetSurfaceFSI
    (
      const LINALG::MapExtractor* fsisurface  //!< the FSI surface
    );

    //! Set forces due to interface with fluid
    void SetForceInterface
    (
      const LINALG::MapExtractor& extractor,
      Teuchos::RCP<Epetra_Vector> iforce  ///< the force on interface
    );

    //! Set forces due to interface with fluid,
    //! the force is expected external-force-like
    void SetForceInterface
    (
      Teuchos::RCP<Epetra_Vector> iforce  //!< the force on interface
    );

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    virtual void EvaluateForceStiffResidualRelax() = 0;

    //! Linear structure solve with just an interface load 
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear();

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

    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;
   
    //! Return residual displacements \f$\Delta D_{n+1}^{<k>}\f$
    Teuchos::RCP<const Epetra_Vector> DisRes() const { return disi_; }

    //! Set residual displacements \f$\Delta D_{n+1}^{<k>}\f$
    void SetDisResidual
    (
      const Teuchos::RCP<const Epetra_Vector> disi  //!< input residual displacements
    )
    {  
      if (disi != Teuchos::null) disi_->Update(1.0, *disi, 0.0);
    }

    //! Return effective residual force \f$R_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> ForceRes() const { return fres_; }

    //@}

  protected:

    //! copy constructor is NOT wanted
    TimIntImpl(const TimIntImpl& old);

    //! @name FSI stuff
    //@{
    const LINALG::MapExtractor* fsisurface_;  //!< surface of structural
                                              //!< which is in contact
                                              //!< with fluid
    //@}

    //! @name General purpose algorithm parameters
    //@{
    enum INPAR::STR::PredEnum pred_;  //!< predictor
    //@}

    //! @name Iterative solution technique
    //@{
    enum INPAR::STR::NonlinSolTech itertype_;  //!< kind of iteration technique 
                                               //!< or non-linear solution technique
    enum INPAR::STR::ConvCheck itercnvchk_;  //!< convergence check
    enum INPAR::STR::VectorNorm iternorm_;  //!< vector norm to 
    //!< check with
    int itermax_;  //!< maximally permitted iterations
    double toldisi_;  //!< tolerance residual displacements
    double tolfres_;  //!< tolerance force residual
    double uzawaparam_;  //!< Parameter for Uzawa algorithm dealing 
    //!< with Lagrange multipliers
    int uzawaitermax_;  //!< maximally permitted Uzawa iterations
    double tolcon_;  //!< tolerance constraint
    int iter_;  //!< iteration step
    double normcharforce_;  //!< characteristic norm for residual force
    double normchardis_;  //!< characteristic norm for residual displacements
    double normfres_;  //!< norm of residual forces
    double normdisi_;  //!< norm of residual displacements
    double normcon_;  //!< norm of constraint
    Teuchos::RCP<Epetra_Vector> disi_;  //!< residual displacements
                                        //!< \f$\Delta{D}^{<k>}_{n+1}\f$
    Epetra_Time timer_;  //!< timer for solution technique
    //@}

    //! @name Various global forces
    //@{
    Teuchos::RCP<Epetra_Vector> fres_;  //!< force residual used for solution
    //@}

    //! @name Interface stuff
    //@{
    Teuchos::RCP<Epetra_Vector> fifc_;  //!< interface force in e.g. FSI
    //@}

    //! @name NOX variables
    //@{
    Teuchos::RCP<NOX::StatusTest::Combo> noxstatustest_;
    Teuchos::RCP<Teuchos::ParameterList> noxparams_;
    Teuchos::RCP<NOX::Utils> noxutils_;
    //@}

  };  // class TimIntImpl

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_IMPL_H
#endif  // #ifdef CCADISCRET
