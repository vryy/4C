/*-----------------------------------------------------------*/
/*!
\file str_impl_generic.H

\maintainer Michael Hiermeier

\date Aug 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_GENERIC_H_
#define STR_IMPL_GENERIC_H_

#include "str_integrator.H"
// required because of used enums
#include <NOX_Abstract_Vector.H>

// forward declaration ...
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace STR {
namespace IMPLICIT {
  class Generic : public Integrator
  {
  public:
    //! constructor
    Generic() {};

    //! destructor
    virtual ~Generic() {};

    //! Setup (has to be implemented by the derived classes)
    virtual void Setup() = 0;

    //! Set state variables (derived)
    virtual void SetState(const Epetra_Vector& x) = 0;

    //! derived
    virtual bool ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f) = 0;

    /*! \brief Apply the stiffness only
     *
     * Normally this one is unnecessary, since it makes more sense
     * to evaluate the stiffness and right hand side at once, because of
     * the lower computational overhead. */
    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac) = 0;

    /*! \brief Apply force and stiff at once
     *
     *  Only one loop over all elements. Especially in the contact case,
     *  the difference between this call and first call ApplyForce and
     *  then ApplyStiff is mentionable because of the projection operations. */
    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac) = 0;

    /*! \brief Calculate characteristic/reference norms for forces
     *
     *  The reference norms are used to scale the calculated iterative
     *  displacement norm and/or the residual force norm. For this
     *  purpose we only need the right order of magnitude, so we don't
     *  mind evaluating the corresponding norms at possibly different
     *  points within the time-step (end point, generalized midpoint). */
    virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type) = 0;

    //! return time integration factor
    virtual double GetIntParam() const = 0;

    //! @name Monolithic update routines
    //! @{
    //! things that should be done before updating (derived)
    virtual void PreUpdate() {/* do nothing for now */};

    //! Update configuration after time step (derived)
    virtual void UpdateStepState() = 0;

    //! Update everything on element level after time step and after output (derived)
    virtual void UpdateStepElement() = 0;

    //! Output routine
    virtual void OutputStepState() = 0;

    //! things that should be done after updating (derived)
    virtual void PostUpdate() {/* do nothing for now */};
    //! @}

    //! @name Predictor routines (dependent on the implicit integration scheme)
    //! @{
    /*! Predict constant displacements, however the velocities
     *  and accelerations are consistent to the time integration
     *  if the constant displacements are taken as correct displacement
     *  solution.
     *  This method has to be implemented by the individual time
     *  integrator. */
    virtual void PredictConstDisConsistVelAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const = 0;

    /*! Predict displacements based on the assumption of constant
     *  velocities. Calculate consistent velocities and accelerations
     *  afterwards.
     *  This method has to be implemented by the individual time
     *  integrator. */
    virtual bool PredictConstVelConsistAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const = 0;

    /*! Predict displacements based on the assumption of constant
     *  accelerations. Calculate consistent velocities and accelerations
     *  afterwards.
     *  This method has to be implemented by the individual time
     *  integrator. */
    virtual bool PredictConstAcc(
        Epetra_Vector& disnp,
        Epetra_Vector& velnp,
        Epetra_Vector& accnp) const = 0;
    //! @}

  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_GENERIC_H_ */
