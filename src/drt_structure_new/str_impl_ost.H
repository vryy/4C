/*-----------------------------------------------------------*/
/*!
\file str_impl_ost.H

\brief One step theta time integrator.

\maintainer Philipp Farah

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_OST_H_
#define STR_IMPL_OST_H_

#include "str_impl_generic.H"

//forward declarations
class Epetra_MultiVector;

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
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! (derived)
    virtual void ReadRestart(IO::DiscretizationReader& ioreader);

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

    /*! \brief things that should be done after updating [derived]
     *
     *  We use in the OneStepTheta case to update constant contributions (during one time step)
     *  of the SetState routine.*/
    virtual void PostUpdate();
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

    /*! \brief Update constant contributions of the current state for the new time step \f$t_{n+1}\f$
     *
     *  The constant contributions are stored in the const_vel_acc_update_ptr_ multi-vector
     *  pointer. The 1st entry represents the velocity, and the 2nd the acceleration.
     *
     *  Calculate the consistent velocity update part:
     *  \f[
     *    V_{n+1}^{0} =   \frac{\beta-\gamma}{\beta} V_{n}
     *                  + \frac{(2 \beta - \gamma) \Delta t}{2 \beta} A_{n}
     *  \f]
     *
     *  Calculate the consistent acceleration update part:
     *  \f[
     *    A_{n+1}^{0} = - \frac{1}{\beta (\Delta t)^{2}} V_{n}
     *                  + \frac{2 \beta - 1}{2 \beta} A_{n}
     *  \f]
     *
     *  Set the new DBC values at time \f$t=t_{n+1}\f$:
     *  \f[
     *    D_{n+1}^{0} \leftarrow DBC(D_{n+1}^{0}, t_{n+1}) \\
     *    V_{n+1}^{0} \leftarrow DBC(V_{n+1}^{0}, t_{n+1}) \\
     *    A_{n+1}^{0} \leftarrow DBC(A_{n+1}^{0}, t_{n+1})
     *  \f]
     *
     *  Finally, add the constant displacement update part to the constant velocity and
     *  acceleration vectors:
     *  \f[
     *    \tilde{V}_{n+1} = V_{n+1}^{0} - \frac{1}{\theta \Delta t} D_{n+1}^{0} \\
     *    \tilde{A}_{n+1} = A_{n+1}^{0} - \frac{1}{theta^{2} (\Delta t)^{2}} D_{n+1}^{0}
     *  \f]
     *
     *  See the SetState() routine for the iterative update of the current state.
     *
     *  \author vuong
     *  \date 04/2016 */
    void UpdateConstantStateContributions();

  protected:
    //! run directly after the setup routine has been finished
    virtual void PostSetup();

  private:
    /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
     *
     * \remark The remaining contributions have been considered in the corresponding model
     *         evaluators. This is due to the fact, that some models use a different
     *         time integration scheme for their terms (e.g. GenAlpha for the structure
     *         and OST for the remaining things).
     *
     *  \f[
     *    Res = M . [\theta  * A_{n+1} + (1-\theta) * A_{n}]
     *        + C . [\theta  * V_{n+1} + \(1-\theta) * V_{n}]
     *        + \theta * Res_{\mathrm{statics},n+1} + (1-\theta) * Res_{\mathrm{statics},n}
     *  \f] */
    void AddViscoMassContributions(Epetra_Vector& f) const;

    /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
     *
     *  \remark The remaining blocks have been considered in the corresponding model
     *          evaluators. This is due to the fact, that some models use a different
     *          time integration scheme for their terms (e.g. GenAlpha for the structure
     *          and OST for the remaining things). Furthermore, constraint/Lagrange
     *          multiplier blocks need no scaling anyway.
     *
     *  \f[
     *    \boldsymbol{K}_{T,effdyn} = (1 - \frac{1}{\theta (\Delta t)^{2}} \boldsymbol{M}
     *                + (1 - \frac{1}{\Delta t} \boldsymbol{C}
     *                + theta  \boldsymbol{K}_{T}
     *  \f] */
    void AddViscoMassContributions(LINALG::SparseOperator& jac) const;

  private:
    //! theta factor: feasible interval (0,1]
    double theta_;

    /*! @name New vectors for internal use only
     *
     *  If an external use seems necessary, move these vectors to the
     *  global state data container and just store a pointer to the global
     *  state variable. */
    //! @{

    //! viscous mid-point force vector F_viscous F_{viscous;n+1}
    Teuchos::RCP<Epetra_Vector> fvisconp_ptr_;

    //! viscous mid-point force vector F_viscous F_{viscous;n}
    Teuchos::RCP<Epetra_Vector> fviscon_ptr_;

    /*! \brief Holds the during a time step constant contributions to
     *  the velocity and acceleration state update.
     *
     *  entry (0): constant velocity contribution \f$\tilde{V}_{n+1}\f$
     *  entry (1): constant acceleration contribution \f$\tilde{A}_{n+1}\f$ */
    Teuchos::RCP<Epetra_MultiVector> const_vel_acc_update_ptr_;
    //! @}

    //! @name pointers to the global state data container content
    //! @{

    //! pointer to inertial force vector F_{inertial,n} at last time
    Teuchos::RCP<Epetra_Vector> finertian_ptr_;

    //! pointer to inertial force vector F_{inertial,n+1} at new time
    Teuchos::RCP<Epetra_Vector> finertianp_ptr_;
    //! @}

  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_OST_H_ */
