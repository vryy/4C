/*-----------------------------------------------------------*/
/*!
\file str_impl_statics.H

\brief Static (time) integrator.

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_STATICS_H_
#define STR_IMPL_STATICS_H_

#include "str_impl_generic.H"
// necessary due to enum usage
#include <NOX_Abstract_Vector.H>

namespace STR {
namespace IMPLICIT {
  class Statics : public Generic
  {
  public:
    //! constructor
    Statics();

    //! destructor
    virtual ~Statics() {};

    //! Setup the class variables
    virtual void Setup();

    //! Reset state variables (derived)
    virtual void SetState(const Epetra_Vector& x);

    //! Apply the rhs only (derived)
    virtual bool ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f);

    //! Apply the stiffness only (derived)
    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac);

    //! Apply force and stiff at once (derived)
    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac);

    //! (derived)
    virtual bool AssembleForce( Epetra_Vector& f,
        const std::vector<INPAR::STR::ModelType>* without_these_models = NULL ) const;

    //! (derived)
    virtual void WriteRestart(IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! (derived)
    virtual void ReadRestart(IO::DiscretizationReader& ioreader);

    //! (derived)
    virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type);

    //! return time integration factor (derived)
    virtual double GetIntParam() const;

    //! @name Monolithic update routines
    //! @{
    //! things that should be done before updating (derived)
    virtual void PreUpdate();

    //! Update configuration after time step (derived)
    virtual void UpdateStepState();

    //! Update everything on element level after time step and after output (derived)
    virtual void UpdateStepElement();
    //! @}

    //! @name Predictor routines (dependent on the implicit integration scheme)
    //! @{
    /*! Predict constant displacements, consistent velocities and accelerations (derived) */
    virtual void PredictConstDisConsistVelAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const;

    /*! Predict displacements based on constant velocities and consistent accelerations (derived) */
    virtual bool PredictConstVelConsistAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const;

    /*! Predict displacements based on constant accelerations and consistent velocities (derived) */
    virtual bool PredictConstAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const;
    //! @}

  protected:
    /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
     *
     * \remark Nothing needs to be done in the static case. */
    virtual void AddViscoMassContributions(Epetra_Vector& f) const {return;};

    /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
     *
     * \remark Nothing needs to be done in the static case. */
    virtual void AddViscoMassContributions(LINALG::SparseOperator& jac) const {return;};

    //! reset the time step dependent parameters for the element evaluation [derived]
    virtual void ResetEvalParams();
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_STATICS_H_ */
