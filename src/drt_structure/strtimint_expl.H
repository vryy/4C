/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.H
\brief Explicit time integration for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
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
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
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
      const AUX::MapExtractor* fsisurface  //!< the FSI surface
    );

    //! Set forces due to interface with fluid
    void SetForceInterface
    (
      const AUX::MapExtractor& extractor,
      Teuchos::RCP<Epetra_Vector> iforce  ///< the force on interface
    );

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

    virtual Teuchos::RCP<Epetra_Vector> Velaf()
    {
      dserror("mid velocities not available");
      return Teuchos::null;
    }

    //@}

  protected:

    //! @name FSI stuff
    //@{
    const AUX::MapExtractor* fsisurface_;  //!< surface of structural
                                             //!< which is in contact
                                             //!< with fluid
    //@}

    //! @name Interface stuff
    //@{
    Teuchos::RCP<Epetra_Vector> fifc_;  //!< interface force in e.g. FSI
    //@}

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_EXPL_H
#endif  // #ifdef CCADISCRET
