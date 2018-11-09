/*-----------------------------------------------------------*/
/*!
\file str_impl_genalpha.H

\brief Generalized Alpha time integrator.

\maintainer Michael Hiermeier

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_GENALPHA_H_
#define STR_IMPL_GENALPHA_H_

#include "str_impl_generic.H"

class Epetra_MultiVector;

namespace STR
{
  namespace IMPLICIT
  {
    /*! Generalized-\f$\alpha\f$ time integration for 2nd-order ODEs
     *
     * <h3> References </h3>
     * - Chung J, Hulbert GM:
     *   A Time Integration Algorithm for Structural Dynamics With Improved Numerical Dissipation:
     *   The Generalized-\f$\alpha\f$ Method
     *   Journal of Applied Mechanics, 60(2):371--375 (1993)
     */
    class GenAlpha : public Generic
    {
      // forwards declaration
      struct Coefficients;

     public:
      //! constructor
      GenAlpha();

      //! destructor
      virtual ~GenAlpha(){};

      //! Setup the class variables
      virtual void Setup();

      //! Reset state variables [derived]
      virtual void SetState(const Epetra_Vector& x);

      //! Apply the rhs only [derived]
      virtual bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f);

      //! Apply the stiffness only [derived]
      virtual bool ApplyStiff(const Epetra_Vector& x, LINALG::SparseOperator& jac);

      //! Apply force and stiff at once [derived]
      virtual bool ApplyForceStiff(
          const Epetra_Vector& x, Epetra_Vector& f, LINALG::SparseOperator& jac);

      //! [derived]
      virtual bool AssembleForce(Epetra_Vector& f,
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL) const;

      virtual bool AssembleJac(LINALG::SparseOperator& jac,
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL) const;

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! [derived]
      virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type);

      //! [derived]
      virtual double GetIntParam() const;

      /// access the parameter for the accelerations at \f$t_{n}\f$, i.e. alpham
      double GetAccIntParam() const override;

      //! @name Monolithic update routines
      //!@{

      //! Update configuration after time step [derived]
      virtual void UpdateStepState();

      //! Update everything on element level after time step and after output [derived]
      virtual void UpdateStepElement();

      /*! \brief things that should be done after updating [derived]
       *
       *  We use in the GenAlpha case to update constant contributions (during one time step)
       *  of the SetState routine.*/
      virtual void PostUpdate();

      //!@}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //!@{

      //! Predict constant displacements, consistent velocities and accelerations [derived]
      virtual void PredictConstDisConsistVelAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      //! Predict displacements based on constant velocities and consistent accelerations [derived]
      virtual bool PredictConstVelConsistAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      //! Predict displacements based on constant accelerations and consistent velocities [derived]
      virtual bool PredictConstAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      //!@}

     protected:
      //! run directly after the setup routine has been finished
      virtual void PostSetup();

      //! reset the time step dependent parameters for the element evaluation [derived]
      virtual void ResetEvalParams();

     private:
      /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
       *
       * \remark The remaining contributions have been considered in the corresponding model
       *         evaluators. This is due to the fact, that some models use a different
       *         time integration scheme for their terms (e.g. GenAlpha for the structure
       *         and OST for the remaining things).
       *
       *  \f[
       *    Res = M . [(1-\alpha_m) * A_{n+1} + \alpha_m * A_{n}]
       *        + C . [(1-\alpha_f) * V_{n+1} + \alpha_f * V_{n}]
       *        + (1-\alpha_f) * Res_{\mathrm{statics},n+1} + \alpha_f * Res_{\mathrm{statics},n}
       *  \f]
       *
       *  \param[in/out] f Right-hand side vector
       *
       *  \author hiermeier
       *  \date 03/2016 */
      virtual void AddViscoMassContributions(Epetra_Vector& f) const;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       *
       *  \remark The remaining blocks have been considered in the corresponding model
       *          evaluators. This is due to the fact, that some models use a different
       *          time integration scheme for their terms (e.g. GenAlpha for the structure
       *          and OST for the remaining things). Furthermore, constraint/Lagrange
       *          multiplier blocks need no scaling anyway.
       *
       *  \f[
       *    \boldsymbol{K}_{T,effdyn} = (1 - \frac{\alpha_m}{\beta (\Delta t)^{2}} \boldsymbol{M}
       *                + (1 - \frac{\alpha_f \gamma}{\beta \Delta t} \boldsymbol{C}
       *                + (1 - \alpha_f)  \boldsymbol{K}_{T}
       *  \f]
       *
       * \param[in/out] jac Jacobian matrix
       *
       *  \author hiermeier
       *  \date 03/2016 */
      virtual void AddViscoMassContributions(LINALG::SparseOperator& jac) const;

      /*! \brief Update constant contributions of the current state for the new time step
       * \f$t_{n+1}\f$
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
       *    \tilde{V}_{n+1} = V_{n+1}^{0} - \frac{\gamma}{\beta \Delta t} D_{n+1}^{0} \\
       *    \tilde{A}_{n+1} = A_{n+1}^{0} - \frac{1}{\beta (\Delta t)^{2}} D_{n+1}^{0}
       *  \f]
       *
       *  See the SetState() routine for the iterative update of the current state.
       *
       *  \author hiermeier
       *  \date 03/2016 */
      virtual void UpdateConstantStateContributions();

      /// set the time integration coefficients
      void SetTimeIntegrationCoefficients(Coefficients& coeffs) const;

      /// Return a reliable model value which can be used for line search
      double GetModelValue(const Epetra_Vector& x) override;

     private:
      /*! \brief Time integration coefficients container
       *
       * \note For bounds on the individual parameters, please consult the
       * original publication by Chung and Hulbert (1993). In practive however,
       * it is advised to choose the spectral radius as \f$\rho_\infty \in [0.5,1]\f$.
       */
      struct Coefficients
      {
        /// copy operator (currently default is enough)
        Coefficients& operator=(const Coefficients& source) = default;

        //! Parameter \f$\beta \in (0,1/2]\f$
        double beta_ = -1.0;

        //! Parameter \f$\gamma \in (0,1]\f$
        double gamma_ = -1.0;

        //! Parameter \f$\alpha_f \in [0,1)\f$
        double alphaf_ = -1.0;

        //! Parameter \f$\alpha_m \in [-1,1)\f$
        double alpham_ = -1.0;

        //! Spectral radius \f$\rho_\infty \in [0,1]\f$
        double rhoinf_ = -1.0;
      };
      Coefficients coeffs_;

     protected:
      /** @name Generalized alpha parameters
       *
       *  \note redirection to the content of the coefficient struct */
      //!@{

      //! Parameter \f$\beta \in (0,1/2]\f$
      double& beta_;

      //! Parameter \f$\gamma \in (0,1]\f$
      double& gamma_;

      //! Parameter \f$\alpha_f \in [0,1)\f$
      double& alphaf_;

      //! Parameter \f$\alpha_m \in [-1,1)\f$
      double& alpham_;

      //! Spectral radius \f$\rho_\infty \in [0,1]\f$
      double& rhoinf_;

      //!@}

      /*! @name New vectors for internal use only
       *
       *  If an external use seems necessary, move these vectors to the
       *  global state data container and just store a pointer to the global
       *  state variable. */
      //!@{

      /*! \brief Holds the during a time step constant contributions to
       *  the velocity and acceleration state update.
       *
       *  entry (0): constant velocity contribution \f$\tilde{V}_{n+1}\f$
       *  entry (1): constant acceleration contribution \f$\tilde{A}_{n+1}\f$ */
      Teuchos::RCP<Epetra_MultiVector> const_vel_acc_update_ptr_;

      //!@}

      //! @name pointers to the global state data container content
      //!@{

      //! viscous force vector F_viscous F_{viscous;n+1}
      Teuchos::RCP<Epetra_Vector> fvisconp_ptr_;

      //! viscous force vector F_viscous F_{viscous;n}
      Teuchos::RCP<Epetra_Vector> fviscon_ptr_;

      //! pointer to inertial force vector F_{inertial,n+1} at new time
      Teuchos::RCP<Epetra_Vector> finertianp_ptr_;

      //! pointer to inertial force vector F_{inertial,n} at last time
      Teuchos::RCP<Epetra_Vector> finertian_ptr_;

      //!@}
    };  // namespace IMPLICIT
  }     // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_GENALPHA_H_ */
