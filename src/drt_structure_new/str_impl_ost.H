/*-----------------------------------------------------------*/
/*!
\file str_impl_ost.H

\maintainer Philipp Farah

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_OST_H_
#define STR_IMPL_OST_H_

#include "str_impl_generic.H"

namespace STR {
namespace IMPLICIT {
  class OneStepTheta : public Generic
  {
  public:
    //! constructor
    OneStepTheta();

    //! destructor
    virtual ~OneStepTheta() {};

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
    virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type);

    //! (derived)
    virtual double GetIntParam() const;

    //! @name Monolithic update routines
    //! @{
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
  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_OST_H_ */
