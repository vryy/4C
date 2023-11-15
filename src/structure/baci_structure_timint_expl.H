/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for structural dynamics

\level 1

*/

/*----------------------------------------------------------------------*/
#ifndef BACI_STRUCTURE_TIMINT_EXPL_H
#define BACI_STRUCTURE_TIMINT_EXPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_structure_timint.H"

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
   *        with \b explicit time integrators
   *
   * <h3> About </h3>
   * This object bridges the gap between the base time integator STR::TimInt
   * and the specific implementation of explicit time integrators.
   *
   * \author bborn
   * \date 07/08
   */
  class TimIntExpl : public TimInt
  {
   public:
    //! @name Life
    //@{

    //! constructor
    TimIntExpl(const Teuchos::ParameterList& timeparams,   //!< time parameters
        const Teuchos::ParameterList& ioparams,            //!< ioflags
        const Teuchos::ParameterList& sdynparams,          //!< input parameters
        const Teuchos::ParameterList& xparams,             //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,          //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,         //!< the solver
        Teuchos::RCP<CORE::LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
        Teuchos::RCP<IO::DiscretizationWriter> output      //!< the output
    );

    //! Empty constructor
    TimIntExpl() : TimInt() { ; }

    //! Copy constructor
    TimIntExpl(const TimIntExpl& old) : TimInt(old) { ; }

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
    void Init(const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<CORE::LINALG::Solver> solver) override;

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
    void Setup() override;

    //@}

    //! @name Actions
    //@{

    //! Resize \p TimIntMStep<T> multi-step quantities
    void ResizeMStep() override = 0;

    //! Do time integration of single step
    int IntegrateStep() override = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next time step.
    void UpdateStepState() override = 0;

    //! Update Element
    void UpdateStepElement() override = 0;
    /*
        //! Update configuration and time after time step
        void UpdateStepAndTime()
        {
          // system state
          UpdateStepState();
          // update time and step
          time_->UpdateSteps(timen_);
          step_ = stepn_;
          //
          timen_ += (*dt_)[0];
          stepn_ += 1;
          // element update
          UpdateStepElement();
        }
    */
    //@}

    //! @name Output
    //@{

    //! print summary after step
    void PrintStep() override;

    //! The text for summary print, see #PrintStep
    void PrintStepText(FILE* ofile  //!< output file handle
    );

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum INPAR::STR::DynamicType MethodName() const override = 0;

    //! These time integrators are all explicit (mark their name)
    bool MethodImplicit() override { return false; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() const override = 0;

    //! Give local order of accuracy of displacement part
    int MethodOrderOfAccuracyDis() const override = 0;

    //! Give local order of accuracy of velocity part
    int MethodOrderOfAccuracyVel() const override = 0;

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const override = 0;

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const override = 0;

    //! return time integration factor
    double TimIntParam() const override { return 0.0; }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override = 0;

    //! Return reaction forces
    Teuchos::RCP<Epetra_Vector> Freact() override
    {
      dserror("Not impl.");
      return Teuchos::null;
    };

    //! Read and set external forces from file
    void ReadRestartForce() override = 0;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) override = 0;

    //! InitialGuess is not available for explicit time integrators
    Teuchos::RCP<const Epetra_Vector> InitialGuess() override
    {
      dserror("InitialGuess() is not available for explicit time integrators");
      return Teuchos::null;
    }

    //! RHS() is not available for explicit time integrators
    Teuchos::RCP<const Epetra_Vector> RHS() override
    {
      dserror("RHS() is not available for explicit time integrators");
      return Teuchos::null;
    }

    //! Prepare time step
    void PrepareTimeStep() override
    {
      // safety checks
      CheckIsInit();
      CheckIsSetup();

      // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
      SetTimen((*time_)[0] + (*dt_)[0]);

      // things that need to be done before Predict
      PrePredict();

      // prepare contact for new time step
      PrepareStepContact();

      return;
    }

    //!  Update displacement state in case of coupled problems
    void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override
    {
      dserror(
          "All monolithically coupled problems work with implicit time "
          "integration schemes. Thus, calling UpdateStateIncrementally() in an explicit scheme "
          "is not possible.");
    }

    //!  Evaluate routine for coupled problems with monolithic approach
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc  ///< iterative solution increment
        ) override
    {
      dserror(
          "All monolithically coupled problems work with implicit time "
          "integration schemes. Thus, calling Evaluate() in an explicit scheme "
          "is not possible.");
    }

    //! Apply external force
    void ApplyForceExternal(const double time,  //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
        const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
        Teuchos::RCP<Epetra_Vector>& fext       //!< external force
    );

    /// has to be renamed either here or PrintStep()
    void Output(bool forced_writerestart) override
    {
      OutputStep(forced_writerestart);
      // write Gmsh output
      WriteGmshStrucOutputStep();
      return;
    }

    /// has to be renamed either here or UpdateStepState() /UpdateStepStateElement()
    void Update() override
    {
      PreUpdate();
      UpdateStepState();
      UpdateStepTime();
      UpdateStepElement();
      PostUpdate();
      return;
    }

    //! Update routine for coupled problems with monolithic approach with time adaptivity
    void Update(const double endtime) override
    {
      dserror("Not implemented. No time adaptivity available for explicit time integration.");
    }


    /* Linear structure solve with just an interface load */
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() override
    {
      dserror("SolveRelaxationLinear() not implemented");
      return Teuchos::null;
    }

    /// are there any algebraic constraints?
    bool HaveConstraint() override
    {
      dserror("HaveConstraint() has not been tested for explicit time integrators");
      return false;
    };

    /// are there any Cardiovascular0D bcs?
    virtual bool HaveCardiovascular0D()
    {
      dserror("HaveCardiovascular0D() has not been tested for explicit time integrators");
      return false;
    };

    /// are there any spring dashpot BCs?
    bool HaveSpringDashpot() override
    {
      dserror("HaveSpringDashpot() has not been tested for explicit time integrators");
      return false;
    };

    //! Return Teuchos::rcp to ConstraintManager conman_
    Teuchos::RCP<::UTILS::ConstrManager> GetConstraintManager() override
    {
      dserror("GetConstraintManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Return Teuchos::rcp to Cardiovascular0DManager windkman_
    virtual Teuchos::RCP<::UTILS::Cardiovascular0DManager> GetCardiovascular0DManager()
    {
      dserror("GetCardiovascular0DManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Return Teuchos::rcp to SpringDashpotManager springman_
    Teuchos::RCP<::UTILS::SpringDashpotManager> GetSpringDashpotManager() override
    {
      dserror("GetSpringDashpotManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Get type of thickness scaling for thin shell structures
    INPAR::STR::STC_Scale GetSTCAlgo() override
    {
      dserror("GetSTCAlgo() has not been tested for explicit time integrators");
      return INPAR::STR::stc_none;
    };

    //! Access to scaling matrix for STC
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override
    {
      dserror("GetSTCMat() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    void UpdateIterIncrConstr(
        Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
        ) override
    {
      dserror("UpdateIterIncrConstr() has not been tested for explicit time integrators");
      return;
    }

    void UpdateIterIncrCardiovascular0D(
        Teuchos::RCP<Epetra_Vector> presincr  ///< pressure increment
        ) override
    {
      dserror("UpdateIterIncrCardiovascular0D() has not been tested for explicit time integrators");
      return;
    }

    /// Do the nonlinear solve, i.e. (multiple) corrector,
    /// for the time step. All boundary conditions have
    /// been set.
    INPAR::STR::ConvergenceStatus Solve() final
    {
      IntegrateStep();
      return INPAR::STR::conv_success;
    }

    //! prepare partiton step
    void PreparePartitionStep() override
    {
      // do nothing for explicit time integrators
      return;
    }

    void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override
    {
      dserror("UseBlockMatrix() not implemented");
    }
    //@}
  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // STRUCTURE_TIMINT_EXPL_H
