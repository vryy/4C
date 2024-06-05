/*----------------------------------------------------------------------*/
/*! \file
\brief Explicit time integration for thermal dynamics
\level 3
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_THERMO_TIMINT_EXPL_HPP
#define FOUR_C_THERMO_TIMINT_EXPL_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_thermo_timint.hpp"

FOUR_C_NAMESPACE_OPEN

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
    TimIntExpl(const Teuchos::ParameterList& ioparams,       //!< ioflags
        const Teuchos::ParameterList& tdynparams,            //!< input parameters
        const Teuchos::ParameterList& xparams,               //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,            //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,           //!< the solver
        Teuchos::RCP<CORE::IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    TimIntExpl() : TimInt() { ; }

    //! Copy constructor
    TimIntExpl(const TimIntExpl& old) : TimInt(old) { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    void ResizeMStep() override = 0;

    //@}

    //! @name Actions
    //@{

    //! Do time integration of single step
    void IntegrateStep() override = 0;

    //! Solve dynamic equilibrium
    //! This is a general wrapper around the specific techniques.
    INPAR::THR::ConvergenceStatus Solve() override
    {
      IntegrateStep();
      return INPAR::THR::conv_success;
    }

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    void Evaluate(Teuchos::RCP<const Epetra_Vector> tempi) override
    {
      FOUR_C_THROW("not implemented for explicit time integration");
      return;
    }

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    void Evaluate() override
    {
      FOUR_C_THROW("not implemented for explicit time integration");
      return;
    }

    //! prepare time step
    void prepare_time_step() override
    {
      // do nothing
      return;
    }

    //! for implicit partitioned schemes
    void prepare_partition_step() override
    {
      // do nothing
      return;
    }

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next time step.
    void UpdateStepState() override = 0;

    //! Update Element
    void UpdateStepElement() override = 0;

    //! update at time step end
    void Update() override;

    //! update Newton step
    void UpdateNewton(Teuchos::RCP<const Epetra_Vector> tempi) override
    {
      FOUR_C_THROW("not needed for explicit time integration");
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
    void PrintStep() override;

    //! The text for summary print, see #PrintStep
    void print_step_text(FILE* ofile  //!< output file handle
    );

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum INPAR::THR::DynamicType MethodName() const override = 0;

    //! These time integrators are all explicit (mark their name)
    bool MethodImplicit() override { return false; }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps() override = 0;

    //! Give local order of accuracy of displacement part
    int method_order_of_accuracy() override = 0;

    //! Return linear error coefficient of temperatures
    double MethodLinErrCoeff() override = 0;

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override = 0;

    //! Return reaction forces
    Teuchos::RCP<Epetra_Vector> Freact() override
    {
      FOUR_C_THROW("Not impl.");
      return Teuchos::null;
    };

    //! initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> initial_guess() override
    {
      FOUR_C_THROW("not needed for explicit time integration");
      return Teuchos::null;
    }

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() override
    {
      FOUR_C_THROW("not needed for explicit time integration");
      return Teuchos::null;
    }

    //! Read and set external forces from file
    void ReadRestartForce() override = 0;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<CORE::IO::DiscretizationWriter> output) override = 0;

    //@}

   protected:
    // currently nothing
  };

}  // namespace THR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
