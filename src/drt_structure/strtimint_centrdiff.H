/*----------------------------------------------------------------------*/
/*!
\file strtimint_centrdiff.H
\brief Structural time integration with central difference scheme 2nd order (explicit)

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_CENTRDIFF_H
#define STRTIMINT_CENTRDIFF_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_expl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Central differences: 2nd order accurate,
   *                             explicit time integrator,
   *                             1-step method
   *
   *  The displacement update reads
   *  \f[
   *  d_{n+1} = d_n + \Delta t \dot{d}_n + \frac{1}{2}\Delta t^2 \ddot{d}_n
   *  \f]
   *
   *  The velocity update reads
   *  \f[
   *  \dot{d}_{n+1} = \dot{d}_n
   *                + \frac{1}{2}\Delta t \ddot{d}_n
   *                + \frac{1}{2}\Delta t \ddot{d}_{n+1}
   *  \f]
   *
   * \author popp
   * \date 01/12
   */
  class TimIntCentrDiff : public TimIntExpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntCentrDiff
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output,  //!< the output
      bool initmassdampconsistaccel = true  //!< bool whether DetermineMassDampConsistAccel() is called
    );

    //! Empty constructor
    TimIntCentrDiff() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntCentrDiff(const TimIntCentrDiff& old) : TimIntExpl(old) { ; }

    //! Destructor
    virtual ~TimIntCentrDiff() { ; }

    //@}

    //! @name Actions
    //@{

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awaiting the next time step.
    void UpdateStepState();

    //! Update Element
    void UpdateStepElement();

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum INPAR::STR::DynamicType MethodName() const
    {
      return INPAR::STR::dyna_centrdiff;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    virtual const int MethodSteps() const { return 1; }

    //! Give local order of accuracy of displacement part
    virtual const int MethodOrderOfAccuracyDis() const { return 2; }

    //! Give local order of accuracy of velocity part
    virtual const int MethodOrderOfAccuracyVel() const { return 2; }

    /*! Return linear error coefficient of displacements
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx \frac{1}{6}\Delta t_n^3 \dddot{d_n} + HOT(\Delta t_n^4)
     *  \f]
     */
    virtual const double MethodLinErrCoeffDis() const { return 1./6.; }

    /*! Return linear error coefficient of velocities
     *
     *  The local discretization error reads
     *  \f[
     *  e \approx -\frac{1}{12}\Delta t_n^3 d_n^{(4)} + HOT(\Delta t_n^4)
     *  \f]
     */
    virtual const double MethodLinErrCoeffVel() const { return -1./12.; }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      return fextn_;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      dserror("FextNew() not available in CentrDiff");
      return Teuchos::null;
    }

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);

    virtual Teuchos::RCP<Epetra_Vector> Velaf()
    {
      dserror("mid velocities not available");
      return Teuchos::null;
    }

    //@}


  protected:

    //! @name Global forces at \f$t_{n+1}\f$
    //@{
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force
                                         //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Epetra_Vector> fviscn_; //!< Rayleigh viscous forces
                                         //!< \f$C \cdot V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fcmtn_;  //!< contact or meshtying forces
                                         //!< \f$F_{cmt;n+1}\f$
    Teuchos::RCP<Epetra_Vector> frimpn_; //!< time derivative of
                                         //!< linear momentum
                                         //!< (temporal rate of impulse)
                                         //!< \f$\dot{P}_{n+1} = M \cdot \dot{V}_{n+1}\f$
    //@}

  };  // class TimIntCentrDiff

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_CENTRDIFF_H
