/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for thermal dynamics
\level 3
\maintainer Sebastian Proell
*/

/*----------------------------------------------------------------------*/
#ifndef THRTIMINT_EXPL_H
#define THRTIMINT_EXPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "thrtimint.H"

/*----------------------------------------------------------------------*/
/* belongs to thermal dynamics namespace */
namespace THR
{
  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics
   *        with \b explicit time integrators
   *
   * <h3> About </h3>
   * This object bridges the gap between the base time integator THR::TimInt
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
    TimIntExpl(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,       //!< input parameters
        const Teuchos::ParameterList& xparams,          //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,       //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,            //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output   //!< the output
    );

    //! Empty constructor
    TimIntExpl() : TimInt() { ; }

    //! Copy constructor
    TimIntExpl(const TimIntExpl& old) : TimInt(old) { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //@}

    //! @name Actions
    //@{

    //! Do time integration of single step
    virtual void IntegrateStep() = 0;

    //! Solve dynamic equilibrium
    //! This is a general wrapper around the specific techniques.
    virtual INPAR::THR::ConvergenceStatus Solve()
    {
      IntegrateStep();
      return INPAR::THR::conv_success;
    }

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> tempi)
    {
      dserror("not implemented for explicit time integration");
      return;
    }

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate()
    {
      dserror("not implemented for explicit time integration");
      return;
    }

    //! prepare time step
    virtual void PrepareTimeStep()
    {
      // do nothing
      return;
    }

    //! for implicit partitioned schemes
    virtual void PreparePartitionStep()
    {
      // do nothing
      return;
    }

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next time step.
    virtual void UpdateStepState() = 0;

    //! Update Element
    virtual void UpdateStepElement() = 0;

    //! update at time step end
    virtual void Update();

    //! update Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> tempi)
    {
      dserror("not needed for explicit time integration");
      return;
    }
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
    void PrintStepText(FILE* ofile  //!< output file handle
    );

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    virtual enum INPAR::THR::DynamicType MethodName() const = 0;

    //! These time integrators are all explicit (mark their name)
    bool MethodImplicit() { return false; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual int MethodSteps() = 0;

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracy() = 0;

    //! Return linear error coefficient of temperatures
    virtual double MethodLinErrCoeff() = 0;

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

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    {
      dserror("not needed for explicit time integration");
      return Teuchos::null;
    }

    //! right-hand side alias the dynamic force residual
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    {
      dserror("not needed for explicit time integration");
      return Teuchos::null;
    }

    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;

    //! Write internal and external forces for restart
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //@}

   protected:
    // currently nothing
  };

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THRTIMINT_EXPL_H
