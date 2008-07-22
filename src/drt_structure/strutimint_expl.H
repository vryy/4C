/*----------------------------------------------------------------------*/
/*!
\file strutimint_expl.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_EXPL_H
#define STRUTIMINT_EXPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for structural dynamics
   *        with *explicit* time integrators
   * \author bborn
   * \date 07/08
   */
  class StruTimIntExpl : public StruTimInt
  {

  public:

    //! @name Life
    //@{

    //! constructor
    StruTimIntExpl
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    StruTimIntExpl() : StruTimInt() { ; }

    //! Copy constructor
    StruTimIntExpl(const StruTimIntExpl& old) : StruTimInt(old) { ; }

    //! Resize #StruTimIntMStep<T> multi-step quantities
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

    //! Output to file
    //! This routine always prints the last converged state, i.e.
    //! \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    void OutputStep();

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    void OutputRestart
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Output displacements, velocities and accelerations
    //! \author mwgee (originally) \date 03/07
    void OutputState
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
    );

    //! Stress & strain output
    //! \author lw (originally)
    void OutputStressStrain
    (
      bool& datawritten  //!< (in/out) read and append if
                         //!< it was written at this time step
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

  protected:
  
    // currently nothing

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_EXPL_H
#endif  // #ifdef CCADISCRET
