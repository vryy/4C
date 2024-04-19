/*----------------------------------------------------------------------*/
/*! \file
\brief Statics analysis, i.e. no capacity terms
\level 1
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 06/09 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_THERMO_TIMINT_STATICS_HPP
#define FOUR_C_THERMO_TIMINT_STATICS_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 06/09 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_thermo_timint_impl.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  /*!
   * \brief Static analysis
   *
   * This static analysis inside the thermal dynamics methods appears
   * slightly displaced, however, it comes in handy in case of
   * thermo-structure-interaction, which is built upon thermal
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
    TimIntStatics(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,          //!< input parameters
        const Teuchos::ParameterList& xparams,             //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,          //!< current discretisation
        Teuchos::RCP<CORE::LINALG::Solver> solver,         //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output      //!< the output
    );

    //! Destructor
    // ....

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() override { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::THR::DynamicType MethodName() const override { return INPAR::THR::dyna_statics; }

    //! Provide number of steps, a single-step method returns 1
    int MethodSteps() override { return 1; }

    //! Give local order of accuracy of temperature part
    int MethodOrderOfAccuracy() override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0;
    }

    //! Return linear error coefficient of temperature
    // virtual double MethodLinErrCoeffTemp()
    double MethodLinErrCoeff() override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0.0;
    }

    //! Consistent predictor with constant temperatures
    //! and consistent temperature rates and temperatures
    void PredictConstTempConsistRate() override;

    //! Evaluate ordinary internal force, its tangent at state
    void ApplyForceTangInternal(const double time,     //!< evaluation time
        const double dt,                               //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,        //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,       //!< residual temperatures
        Teuchos::RCP<Epetra_Vector> fint,              //!< internal force
        Teuchos::RCP<CORE::LINALG::SparseMatrix> tang  //!< tangent matrix
    );

    //! Evaluate ordinary internal force
    void ApplyForceInternal(const double time,    //!< evaluation time
        const double dt,                          //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,   //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
        Teuchos::RCP<Epetra_Vector> fint          //!< internal force
    );

    //! Evaluate a convective boundary condition
    // (nonlinear --> add term to tangent)
    void ApplyForceExternalConv(const double time,     //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> tempn,       //!< temperature state T_n
        const Teuchos::RCP<Epetra_Vector> temp,        //!< temperature state T_n+1
        Teuchos::RCP<Epetra_Vector> fext,              //!< internal force
        Teuchos::RCP<CORE::LINALG::SparseMatrix> tang  //!< tangent matrix
    );

    //! Create force residual #fres_ and its tangent #tang_
    void EvaluateRhsTangResidual() override;

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    // void EvaluateForceTangResidualRelax();

    //! Determine characteristic norm for temperatures
    //! \author lw (originally)
    double CalcRefNormTemperature() override;

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce() override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #tempn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    void UpdateIterIncrementally() override;

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #tempn_ and #raten_ by
    //! incrementing these vector proportional to the residual
    //! temperatures #tempi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual temperatures #tempi_ are blanked at these DOFs.
    void UpdateIterIteratively() override;

    //! Update step
    void UpdateStepState() override;

    //! Update element
    void UpdateStepElement() override;

    //! Read and set restart for forces
    void ReadRestartForce() override;

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) override;

    //@}

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext() override { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() override { return fextn_; }

    //@}

   protected:
    //! equal operator is hidden
    TimIntStatics operator=(const TimIntStatics& old);

    //! copy constructor is hidden
    TimIntStatics(const TimIntStatics& old);

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;   //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;   //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$
    //@}

  };  // class TimIntStatics

}  // namespace THR


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
