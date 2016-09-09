/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.H
\brief Statics analysis
\level 1
\maintainer Alexander Popp
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
      const Teuchos::ParameterList& timeparams,  //!< ioflags
      const Teuchos::ParameterList& ioparams,    //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,     //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,       //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

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
    virtual int MethodSteps() const { return 1; }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const
    {
      dserror("Sensible to ask?");
      return 0;
    }

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const
    {
      dserror("Sensible to ask?");
      return 0;
    }

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const
    {
      dserror("Sensible to ask?");
      return 0.0;
    }

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const
    {
      dserror("Sensible to ask?");
      return 0.0;
    }

    //! return time integration factor
    virtual double TimIntParam() const { return 0.0; }

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
