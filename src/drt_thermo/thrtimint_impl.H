/*----------------------------------------------------------------------*/
/*! \file
\brief base class for all implicit time integrators in thermo-field
\level 1
*/

/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef THRTIMINT_IMPL_H
#define THRTIMINT_IMPL_H

/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "thrtimint.H"
#include "thr_aux.H"
#include <Epetra_Time.h>

// forward declaration
namespace ADAPTER
{
  class CouplingMortar;
}

/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  //!
  //! \brief Front-end for thermal dynamics
  //!        with \b implicit time integration
  //!
  //! <h3> About </h3>
  //! The implicit time integrator object is a derivation of the base time integrators with an eye
  //! towards implicit time integration. #TimIntImpl provides the environment needed to execute
  //! implicit integrators. This is chiefly the non-linear solution technique, e.g., Newton-Raphson
  //! iteration. These iterative solution techniques require a set of control parameters which are
  //! stored within this object. It is up to derived object to implement the time-space discretised
  //! residuum an d its tangent. This object provides some utility functions to obtain various force
  //! vectors necessary in the calculation of the force residual in the derived time integrators.
  //!
  //! \author bborn
  //! \date 06/08
  class TimIntImpl : public TimInt
  {
   public:
    //! @name Construction
    //@{

    //! Constructor
    TimIntImpl(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,       //!< input parameters
        const Teuchos::ParameterList& xparams,          //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,       //!< current discretization
        Teuchos::RCP<LINALG::Solver> solver,            //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output   //!< the output
    );

    //! Destructor
    virtual ~TimIntImpl() { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //@}

    //! Do time integration of single step
    virtual void IntegrateStep();

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> tempi);

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate();

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
    virtual void PreparePartitionStep();

    //! Predict constant temperature, temperature rate,
    //! i.e. the initial guess is equal to the last converged step
    //! except Dirichlet BCs
    void PredictConstTempRate();

    //! Predict constant temperature, however the rate
    //! is consistent to the time integration
    //! if the constant temperature is taken as correct temperature
    //! solution.
    //! This method has to be implemented by the individual time
    //! integrator.
    virtual void PredictConstTempConsistRate() = 0;

    //! Predict temperature which satisfy exactly the Dirichlet BCs
    //! and the linearised system at the previously converged state.
    //!
    //! This is an implicit predictor, i.e. it calls the solver once.
    void PredictTangTempConsistRate();

    //! prepare time step
    virtual void PrepareTimeStep();

    // finite difference check for the tangent K_TT
    void FDCheck();

    //@}

    //! @name Forces and tangents
    //@{

    //! Do residual force due to global balance of energy
    //! and its tangent with respect to the current
    //! temperatures \f$T_{n+1}\f$
    //!
    //! This is <i>the</i> central method which is different for each
    //! derived implicit time integrator. The time integrator implementation
    //! is expected to set members #fres_ and #tang_.
    //! The residual #fres_ is expected to follow the <i>same</i> sign
    //! convention like its tangent #tang_, i.e. to use
    //! Newton--Raphson's method the residual will be scaled by -1.
    virtual void EvaluateRhsTangResidual() = 0;

    //@}

    //! @name Solution
    //@{

    //! determine characteristic norms for relative
    //! error checks of residual temperatures
    //! \author lw  \date 12/07
    virtual double CalcRefNormTemperature() = 0;

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
    virtual INPAR::THR::ConvergenceStatus Solve();

    //! Do full Newton-Raphson iteration
    //!
    //! This routines expects a prepared negative reisdual force #fres_
    //! and associated effective tangent matrix #tang_
    virtual INPAR::THR::ConvergenceStatus NewtonFull();


    void ApplyLatentHeatIntegration();

    //! Blank Dirichlet dofs form residual and reactions
    //! calculate norms for convergence checks
    void BlankDirichletAndCalcNorms();

    // check for success of nonlinear solve
    INPAR::THR::ConvergenceStatus NewtonFullErrorCheck();

    //! Do (so-called) modified Newton-Raphson iteration in which
    //! the initial tangent is kept and not adapted to the current
    //! state of the temperature solution
    void NewtonModified() { dserror("Not impl."); }

    //! Prepare system for solving with Newton's method
    //!
    //! - negative residual
    //! - blank residual on Dirichlet DOFs
    //! - apply Dirichlet boundary conditions on system
    void PrepareSystemForNewtonSolve();

    //@}

    //! @name Updates
    //@{

    //! Update iteration
    //!
    //! This handles the iterative update of the current
    //! temperature \f$T_{n+1}\f$ with the residual temperature
    //! The temperature rate follow on par.
    void UpdateIter(const int iter  //!< iteration counter
    );

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #tempn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    virtual void UpdateIterIncrementally() = 0;

    //! Update iteration incrementally with prescribed residual
    //! temperatures
    void UpdateIterIncrementally(
        const Teuchos::RCP<const Epetra_Vector> tempi  //!< input residual temperatures
    );

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #tempn_ and #raten_ by
    //! incrementing these vector proportional to the residual
    //! temperatures #tempi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual temperatures #tempi_ are blanked at these DOFs.
    virtual void UpdateIterIteratively() = 0;

    //! Update configuration after time step
    //!
    //! This means, the state set
    //! \f$T_{n} := T_{n+1}\f$ and \f$R_{n} := R_{n+1}\f$
    //! Thus the 'last' converged state is lost and a reset
    //! of the time step becomes impossible.
    //! We are ready and keen awaiting the next time step.
    virtual void UpdateStepState() = 0;

    //! Update Element
    virtual void UpdateStepElement() = 0;

    //! update time step
    virtual void Update();

    //! update Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> tempi);

    //@}


    //! @name Output
    //@{

    //! Print to screen predictor informations about residual norm etc.
    //! \author lw (originally) \date 12/07
    void PrintPredictor();

    //! Print to screen information about residual forces and temperatures
    //! \author lw (originally) \date 12/07
    void PrintNewtonIter();

    //! Contains text to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterText(FILE* ofile  //!< output file handle
    );

    //! Contains header to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterHeader(FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void PrintNewtonConv();

    //! print summary after step
    void PrintStep();

    //! The text for summary print, see #PrintStep
    void PrintStepText(FILE* ofile  //!< output file handle
    );

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    virtual enum INPAR::THR::DynamicType MethodName() const = 0;

    //! These time integrators are all implicit (mark their name)
    bool MethodImplicit() { return true; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() = 0;

    //! Give local order of accuracy of temperature part
    virtual int MethodOrderOfAccuracy() = 0;

    //! Return linear error coefficient of temperatures
    virtual double MethodLinErrCoeff() = 0;

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
    //! Mark, this is not true for DBCs with local coordinate
    //! systems in which the non-global reaction force
    //! component is stored in global Cartesian components.
    //! The reaction force resultant is not affected by
    //! this operation.
    Teuchos::RCP<Epetra_Vector> Freact() { return freact_; }

    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;

    //! Write internal and external forces for restart
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! Return residual temperatures \f$\Delta T_{n+1}^{<k>}\f$
    Teuchos::RCP<const Epetra_Vector> TempRes() const { return tempi_; }

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() { return tempi_; }

    //! Set residual temperatures \f$\Delta T_{n+1}^{<k>}\f$
    void SetTempResidual(
        const Teuchos::RCP<const Epetra_Vector> tempi  //!< input residual temperatures
    )
    {
      if (tempi != Teuchos::null) tempi_->Update(1.0, *tempi, 0.0);
    }

    //! Return effective residual force \f$R_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> ForceRes() const { return fres_; }

    //! right-hand side alias the dynamic force residual
    virtual Teuchos::RCP<const Epetra_Vector> RHS() { return fres_; }

    //@}

   protected:
    //! copy constructor is NOT wanted
    TimIntImpl(const TimIntImpl& old);

    // called when unconverged AND dicvont_halve_step
    void HalveTimeStep();

    void CheckForTimeStepIncrease();

    //! @name General purpose algorithm parameters
    //@{
    enum INPAR::THR::PredEnum pred_;  //!< predictor
    //@}

    //! @name Iterative solution technique
    //@{
    enum INPAR::THR::NonlinSolTech itertype_;  //!< kind of iteration technique
                                               //!< or non-linear solution technique
    enum INPAR::THR::ConvNorm normtypetempi_;  //!< convergence check for residual temperatures
    enum INPAR::THR::ConvNorm normtypefres_;   //!< convergence check for residual forces

    enum INPAR::THR::BinaryOp
        combtempifres_;  //!< binary operator to combine temperatures and forces

    enum INPAR::THR::VectorNorm iternorm_;    //!< vector norm to check with
    int itermax_;                             //!< maximally permitted iterations
    int itermin_;                             //!< minimally requested iterations
    enum INPAR::THR::DivContAct divcontype_;  // what to do when nonlinear solution fails
    int divcontrefinelevel_;                  //!< refinement level of adaptive time stepping
    int divcontfinesteps_;  //!< number of time steps already performed at current refinement level
    double toltempi_;       //!< tolerance residual temperatures
    double tolfres_;        //!< tolerance force residual
    int iter_;              //!< iteration step
    int resetiter_;  //<! number of iterations already performed in resets of the current step
    double normcharforce_;                 //!< characteristic norm for residual force
    double normchartemp_;                  //!< characteristic norm for residual temperatures
    double normfres_;                      //!< norm of residual forces
    double normtempi_;                     //!< norm of residual temperatures
    Teuchos::RCP<Epetra_Vector> tempi_;    //!< residual temperatures
                                           //!< \f$\Delta{T}^{<k>}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> tempinc_;  //!< sum of temperature vectors already applied,
                                           //!< i.e. the incremental temperature
    Epetra_Time timer_;                    //!< timer for solution technique
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;  //!< mortar coupling adapter
    //@}

    //! @name Various global forces
    //@{
    Teuchos::RCP<Epetra_Vector> fres_;    //!< force residual used for solution
    Teuchos::RCP<Epetra_Vector> freact_;  //!< reaction force
    //@}

    //! @name melting stuff
    //@{
    Teuchos::RCP<Epetra_Vector> fmelt_;  //!< artificial source term accounting for latent heat
    const bool heatint_;                 //!< should heat integration be performed
    const double melttol_;  //!< tolerance until which heat source increments are applied

    //@}

  };  // class TimIntImpl

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THRTIMINT_IMPL_H
