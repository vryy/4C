/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.H
\brief Statics analysis

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_STATICS_H
#define STRTIMINT_STATICS_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_impl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Static analysis
   *
   * This static analysis inside the structural dynamics methods appears
   * slightly displaced, however, it comes in handy in case of
   * fluid-structure-interaction, which is built upon structural
   * dynamics.
   *
   * Regarding this matter, please direct any complaints to Axel Gerstenberger.
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntStatics : public TimIntImpl
  {

  public:

    //! @name Construction
    //@{

    //! Constructor
    TimIntStatics
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here exept when doing optimization
    void ResizeMStep(){ ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const
    {
      return INPAR::STR::dyna_statics;
    }

    //! Provide number of steps, a single-step method returns 1
    virtual const int MethodSteps() const { return 1; }

    //! Give local order of accuracy of displacement part
    virtual const int MethodOrderOfAccuracyDis() const
    {
      dserror("Sensible to ask?");
      return 0;
    }

    //! Give local order of accuracy of velocity part
    virtual const int MethodOrderOfAccuracyVel() const
    {
      dserror("Sensible to ask?");
      return 0;
    }

    //! Return linear error coefficient of displacements
    virtual const double MethodLinErrCoeffDis() const
    {
      dserror("Sensible to ask?");
      return 0.0;
    }

    //! Return linear error coefficient of velocities
    virtual const double MethodLinErrCoeffVel() const
    {
      dserror("Sensible to ask?");
      return 0.0;
    }

    //! return time integration factor
    virtual const double TimIntParam() const { return 0.0; }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void PredictConstDisConsistVelAcc();

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! For quasi-static problems this is equivalent to
    //! a linear extrapolation of the displacement field.
    //! In the first step we do TangDis
    void PredictConstVelConsistAcc() ;

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! For quasi-static problems this is equivalent to
    //! a quadratic extrapolation of the displacement field.
    //! In the first step we do TangDis, in the second ConstVel
    void PredictConstAcc() ;

    //! Create force residual #fres_ and its stiffness #stiff_
    void EvaluateForceStiffResidual(Teuchos::ParameterList& params);

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    void EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params);

    //! Evaluate residual #fres_
    void EvaluateForceResidual();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_
    //! from scratch by using the newly updated #disn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstDisConsistVelAcc)
    void UpdateIterIncrementally();

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void UpdateIterIteratively();

    //! Update step
    void UpdateStepState();

    //! Update element
    void UpdateStepElement();

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);
    //@}

    virtual void ApplyDirichletBC
    (
      const double time,  //!< at time
      Teuchos::RCP<Epetra_Vector> dis,  //!< displacements
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> vel,  //!< velocities
                                        //!< (may be Teuchos::null)
      Teuchos::RCP<Epetra_Vector> acc,  //!< accelerations
                                       //!< (may be Teuchos::null)
      bool recreatemap  //!< recreate mapextractor/toggle-vector
                        //!< which stores the DOF IDs subjected
                        //!< to Dirichlet BCs
                        //!< This needs to be true if the bounded DOFs
                        //!< have been changed.
    );

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      dserror("Statics: no external forces at t_n available");
      return Teuchos::null;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      return fextn_;
    }

    //@}

    /// update all the field vectors to account for new node introduced by crack
    void updateMethodSpecificEpetraCrack( std::map<int,int>& oldnew );

  protected:

    //! equal operator is hidden
    TimIntStatics operator=(const TimIntStatics& old);

    //! copy constructor is hidden
    TimIntStatics(const TimIntStatics& old);

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$

    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;  //!< internal force at \f$t_n\f$

    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$
    //@}

  };  // class TimIntStatics

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_STATICS_H
