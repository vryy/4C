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
#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics
   *        with \b implicit time integration
   * \author bborn
   * \date 06/08
   */
  class TimIntImpl : public TimInt
  {

  public:

    //! NOX interface is friend of this class
    //friend class NoxInterface;

    //! type of solution techniques
    enum SolTechEnum
    {
      soltech_vague,  //!< undefined
      soltech_newtonfull,  //!< full Newton-Raphson iteration
      soltech_newtonmod,  //!< modified Newton-Raphson iteration
      soltech_uzawalinnewton,  //!< linear Uzawa iteration for 
      //!< constraint system
      soltech_uzawanonlinnewton  //!< non-linear Uzawa iteration 
      //!< for constraint system
    };

    //! map solution technique identification string to enum term
    static enum SolTechEnum MapSolTechStringToEnum
    (
      const std::string name  //!< name identification string
    );

    //! map solution technique enum term identification string
    static std::string MapSolTechEnumToString
    (
      const enum SolTechEnum name  //!< identifying enum term
    );

    //! kind of iteration convergence checks
    //! \author lw \date 12/07
    enum ConvCheckEnum
    {
      convcheck_absres_or_absdis,  //!< absolute norms of residual forces
                                   //!< OR iterative displacement increments
      convcheck_absres_and_absdis,  //!< absolute norms of residual forces
                                    //!< AND iterative displacement increments
      convcheck_relres_or_absdis,  //!< relative norm of residual forces
                                   //!< OR absolute norm if iterative
                                   //!< displacement increments
      convcheck_relres_and_absdis,  //!< relative norm of residual forces
                                    //!< AND absolute norm if iterative
                                    //!< displacement increments
      convcheck_relres_or_reldis,  //!< relative norms of residual forces
                                   //!< OR iterative displacement increments
      convcheck_relres_and_reldis,  //!< relative norms of residual forces
                                    //!< AND iterative displacement increments
      convcheck_vague  //!< undetermined
    };

    //! map convergence check to enum term
    static enum ConvCheckEnum MapConvCheckStringToEnum
    (
      const std::string name  //!< name identification string
    );

    //! Type of predictor
    enum PredEnum
    {
      pred_vague,  //!< undetermined
      pred_constdis,  //!< constant displacements
      pred_constdisvelacc  //!< constant displacements, 
                           //!< velocities and accelerations
    };

    //! Map enum term to input string
    static enum PredEnum MapPredictorStringToEnum
    (
      const std::string name  //!< identifier
    );

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
      Teuchos::ParameterList pcon   //!<Parameterlist containing scale factors for matrix entries
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
      const Teuchos::RCP<Epetra_Vector> disi  //!< input residual displacements
    )
    {
      // select residual displacements
      if (disi != Teuchos::null)
        disi_->Update(1.0, *disi, 0.0);  // set the new solution we just got
      else
        disi_->PutScalar(0.0);
      // Update
      UpdateIterIncrementally();
    }

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
    //! This means, the state set
    //! \f$D_{n} := D_{n+1}\f$, \f$V_{n} := V_{n+1}\f$ and \f$A_{n} := A_{n+1}\f$
    //! Thus the 'last' converged state is lost and a reset 
    //! of the time step becomes impossible.
    //! We are ready and keen awating the next time step.
    virtual void UpdateStep() = 0;

    //! Update configuration and time after time step
    void UpdateStepAndTime()
    {
      // system state
      UpdateStep();
      // update time and step
      time_->UpdateSteps(timen_);
      step_ = stepn_;
      // 
      timen_ += (*dt_)[0];
      stepn_ += 1;
    }

    //! Update step for surface stress
    void UpdateStepSurfstress()
    {
      if (surfstressman_ != Teuchos::null) surfstressman_->Update();
    }

    //! Update step for potential forces
    void UpdateStepPotential()
    {
      if (potman_ != Teuchos::null) potman_->Update();
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
    Teuchos::RCP<Epetra_Vector> DisRes() { return disi_; }

    //! Set residual displacements \f$\Delta D_{n+1}^{<k>}\f$
    void SetDisResidual
    (
      const Teuchos::RCP<Epetra_Vector> disi  //!< input residual displacements
    )
    {  
      if (disi != Teuchos::null) disi_->Update(1.0, *disi, 0.0);
    }

    //! Return effective residual force \f$R_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ForceRes() { return fres_; }

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
    enum PredEnum pred_;  //!< predictor
    //@}

    //! @name Iterative solution technique
    //@{
    enum SolTechEnum itertype_;  //!< kind of iteration technique
    enum ConvCheckEnum itercnvchk_;  //!< convergence check
    enum TimIntVector::NormEnum iternorm_;  //!< vector norm to 
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

  };  // class TimIntImpl

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_IMPL_H
#endif  // #ifdef CCADISCRET
