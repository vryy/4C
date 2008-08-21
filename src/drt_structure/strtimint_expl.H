/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_EXPL_H
#define STRTIMINT_EXPL_H

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
   *        with \b explicit time integrators
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
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
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

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step 
    //! becomes impossible. We are ready and keen awating the next time step.
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
    virtual enum NameEnum MethodName() = 0;

    //! These time integrators are all explicit (mark their name)
    bool MethodImplicit() { return false; }

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

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;
    
    //! Read and set external forces from file
    virtual void ReadRestartForce() = 0;
  

    //@}

  protected:
  
    // currently nothing

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_EXPL_H
#endif  // #ifdef CCADISCRET
