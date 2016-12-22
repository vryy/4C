/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.H

\brief Explicit time integration for structural dynamics

\level 1

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238

*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_EXPL_H
#define STRTIMINT_EXPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint.H"

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
    TimIntExpl
    (
      const Teuchos::ParameterList& timeparams,  //!< time parameters
      const Teuchos::ParameterList& ioparams,    //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,     //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,       //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    TimIntExpl() : TimInt() { ; }

    //! Copy constructor
    TimIntExpl(const TimIntExpl& old) : TimInt(old) { ; }

    //! Destructor
    virtual ~TimIntExpl() { ; };

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

    //@}

    //! @name Actions
    //@{

    //! Resize \p TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next time step.
    virtual void UpdateStepState() = 0;

    //! Update Element
    virtual void UpdateStepElement() = 0;
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
    virtual enum INPAR::STR::DynamicType MethodName() const = 0;

    //! These time integrators are all explicit (mark their name)
    bool MethodImplicit() { return false; }

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

    //! return time integration factor
    virtual double TimIntParam() const { dserror("TimeFac not yet implemented"); return -1.0; }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return reaction forces
    Teuchos::RCP<Epetra_Vector> Freact()
    {
      dserror("Not impl.");
      return Teuchos::null;
    };

    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;

    //! Write internal and external forces for restart
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! InitialGuess is not available for explicit time integrators
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    {
      dserror("InitialGuess() is not available for explicit time integrators");
      return Teuchos::null;
    }

    //! RHS() is not available for explicit time integrators
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    {
      dserror("RHS() is not available for explicit time integrators");
      return Teuchos::null;
    }

    //! Prepare time step
    virtual void PrepareTimeStep()
    {
      // safety checks
      CheckIsInit();
      CheckIsSetup();

      // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
      SetTimen( (*time_)[0]+(*dt_)[0] );

      // things that need to be done before Predict
      PrePredict();

      // prepare contact for new time step
      PrepareStepContact();

      // things that need to be done after Predict
      PostPredict();

      return;
    }

    //!  Evaluate routine for coupled problems with monolithic approach
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> disiterinc ///< iterative solution increment
    )
    {
      dserror("All monolithically coupled problems work with implicit time "
          "integration schemes. Thus, calling Evaluate() in an explicit scheme "
          "is not possible.");
    }

    //! Apply external force
    void ApplyForceExternal
    (
      const double time,  //!< evaluation time
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector>& fext  //!< external force
    );

    /// has to be renamed either here or PrintStep()
    virtual void Output(bool forced_writerestart)
    {
      PreOutput();
      OutputStep(forced_writerestart);
      // write Gmsh output
      writeGmshStrucOutputStep();
      return;
    }

    /// has to be renamed either here or UpdateStepState() /UpdateStepStateElement()
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
      dserror("Not implemented. No time adaptivity available for explicit time integration.");
    }


    /* Linear structure solve with just an interface load */
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear()
    {
      dserror("SolveRelaxationLinear() not implemented");
      return Teuchos::null;
    }

    /// are there any algebraic constraints?
    virtual bool HaveConstraint()
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
    virtual bool HaveSpringDashpot()
    {
      dserror("HaveSpringDashpot() has not been tested for explicit time integrators");
      return false;
    };

    //! Return Teuchos::rcp to ConstraintManager conman_
    virtual const Teuchos::RCP< ::UTILS::ConstrManager> GetConstraintManager()
    {
      dserror("GetConstraintManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Return Teuchos::rcp to Cardiovascular0DManager windkman_
    virtual const Teuchos::RCP< ::UTILS::Cardiovascular0DManager> GetCardiovascular0DManager()
    {
      dserror("GetCardiovascular0DManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Return Teuchos::rcp to SpringDashpotManager springman_
    virtual const Teuchos::RCP< ::UTILS::SpringDashpotManager> GetSpringDashpotManager()
    {
      dserror("GetSpringDashpotManager() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    //! Get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo()
    {
      dserror("GetSTCAlgo() has not been tested for explicit time integrators");
      return INPAR::STR::stc_none;
    };

    //! Access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    {
      dserror("GetSTCMat() has not been tested for explicit time integrators");
      return Teuchos::null;
    };

    virtual void UpdateIterIncrConstr
    (
      Teuchos::RCP<Epetra_Vector> lagrincr ///< Lagrange multiplier increment
    )
    {
      dserror("UpdateIterIncrConstr() has not been tested for explicit time integrators");
      return;
    }

    virtual void UpdateIterIncrCardiovascular0D
    (
      Teuchos::RCP<Epetra_Vector> presincr ///< pressure increment
    )
    {
      dserror("UpdateIterIncrCardiovascular0D() has not been tested for explicit time integrators");
      return;
    }

    /// Do the nonlinear solve, i.e. (multiple) corrector,
    /// for the time step. All boundary conditions have
    /// been set.
    virtual INPAR::STR::ConvergenceStatus Solve()
    {
      IntegrateStep();
      return INPAR::STR::conv_success;
    }

    //! prepare partiton step
    virtual void PreparePartitionStep()
    {
      // do nothing for explicit time integrators
      return;
    }

    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                                Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
    {
      dserror("UseBlockMatrix() not implemented");
    }
    //@}

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_EXPL_H
