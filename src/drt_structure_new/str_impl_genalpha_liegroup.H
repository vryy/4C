/*-----------------------------------------------------------*/
/*!

\brief Modified Generalized Alpha time integrator for Lie groups
       (can handle non-additive rotational pseudo-vector DoFs)

\maintainer Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_GENALPHA_LIEGROUP_H_
#define STR_IMPL_GENALPHA_LIEGROUP_H_

#include "str_impl_genalpha.H"

class Epetra_MultiVector;
namespace STR
{
  namespace IMPLICIT
  {
    /*! \brief Modified GenAlpha time integrator for Lie groups
     *
     *  We need this for all beam elements with rotation vector
     *  DoFs (currently beam3r, beam3k) because those are non-
     *  additive and the corresponding rotation tensors are
     *  members of the special orthogonal group SO3 which is no
     *  vector space!
     *
     *  The implementation is based on [Arnold, Brüls (2007)],
     *  [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)].
     *  See also Christoph Meier's dissertation for explanations,
     *  governing equations, differences to classical GenAlpha,
     *  comments on why we chose this time discretization, ...
     *
     *  In theory, this should be a generalization of classical
     *  GenAlpha, which yields the same results if applied to
     *  problems whose primary DoFs are elements of vector spaces.
     *
     *  The main difference to classical GenAlpha is that a
     *  modified acceleration vector is introduced (see accn_mod_)
     *  and that the weak form is evaluated at t_{n+1} only.
     *
     *  Furthermore, the update of non-additive rotation vector
     *  DoFs needs a special (multiplicative) treatment. However,
     *  for the implemented beam elements with rotation vector DoFs,
     *  the calculation of angular velocities and accelerations
     *  from temporal derivatives of the primary rotation vector DoFs
     *  is intricate (because of non-linear triad interpolation
     *  schemes). Therefore, these beam elements first discretize
     *  in time and apply the Newmark scheme to space-continuous
     *  angular velocities and accelerations. All this (including
     *  storage of old state variables) happens inside the beam
     *  elements. The beam elements are therefore independent from
     *  the velocity and acceleration vectors on time integrator level.
     *
     *  However, in mixed discretizations (beams and other elements),
     *  we need the global velocity and acceleration vector to evaluate
     *  the inertia terms of all other than beam elements in the usual manner.
     *  Note: because of the non-linear inertia forces and mass matrix,
     *        the beams require "MASSLIN = rotations" which evaluates
     *        inertia force and mass matrix in every iteration by calling
     *        Evaluate with action type "nlnstiffmass". Because this
     *        setting (ml_rotations) is true for entire discretization,
     *        all other than beam elements must adapt to this and also
     *        evaluate the inertia force and mass matrix inside the element.
     *
     *  Possible future work (to make this even nicer):
     *  We could use the global velocity and acceleration vectors
     *  for all translational (and hence additive) DoFs to compute
     *  the inertia terms inside the element but the rotation vector
     *  DoFs anyway need the element-internal treatment described above.
     *
     *  Note: - viscous Rayleigh damping is supported in principle, but
     *          never applied or tested yet!
     *
     *        - this is not yet rigorously tested with other than beam
     *          and rigidsphere elements and thus experimental code
     *
     *  \author grill
     *  \date 07/16*/

    class GenAlphaLieGroup : public GenAlpha
    {
     public:
      //! constructor
      GenAlphaLieGroup();

      //! destructor
      virtual ~GenAlphaLieGroup(){};

      //! Setup the class variables [derived]
      virtual void Setup();

      //! Reset state variables [derived]
      virtual void SetState(const Epetra_Vector& x);

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! [derived]
      virtual double GetIntParam() const;

      //! @name Monolithic update routines
      //! @{
      //! Update configuration after time step [derived]
      virtual void UpdateStepState();
      //! @}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //! @{
      /*! Predict constant displacements, consistent velocities and accelerations [derived] */
      virtual void PredictConstDisConsistVelAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      /*! Predict displacements based on constant velocities
       *  and consistent accelerations [derived] */
      virtual bool PredictConstVelConsistAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      /*! Predict displacements based on constant accelerations
       *  and consistent velocities [derived] */
      virtual bool PredictConstAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;
      //! @}
     protected:
      //! run directly after the setup routine has been finished [derived]
      virtual void PostSetup();

      //! reset the time step dependent parameters for the element evaluation [derived]
      virtual void ResetEvalParams();

     private:
      /*! \brief Calculate the right-hand-side vector at \f$ t_{n+1} \f$  [derived]
       *
       *  \f[
       *    Res = F_{inert;n+1}
       *        + F_{visco;n+1}
       *        + F_{int;n+1}
       *        - F_{ext;n+1}
       *  \f]
       *     * Remark: In the case of a Lie group Gen-Alpha time integration scheme,
       *               all forces are evaluated at the end point n+1. */
      virtual void AddViscoMassContributions(Epetra_Vector& f) const;

      /*! \brief Calculate the structural stiffness block at \f$ t_{n+1} \f$  [derived]
       *
       *  \f[
       *    \boldsymbol{K}_{T,effdyn} = \frac{1 - \alpha_m}{\beta (\Delta t)^{2} (1 - \alpha_f) }
       * \boldsymbol{M}
       *                                + \frac{\gamma}{\beta \Delta t} \boldsymbol{C}
       *                                + \boldsymbol{K}_{T}
       *  \f] */
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
       *                  + \frac{(2 \beta - \gamma) \Delta t}{2 \beta} A_{mod,n}
       *  \f]
       *
       *  Calculate the consistent acceleration update part:
       *  \f[
       *    A_{n+1}^{0} = - \frac{1 - \alpha_m}{\beta \Delta t (1-\alpha_f)} V_{n}
       *                  - \frac{\alpha_f}{1-\alpha_f} A_{n}
       *                  + ( \frac{\alpha_m}{1-\alpha_f} - \frac{(1-\alpha_m)(0.5-\beta)}{\beta
       * (1-\alpha_f)} ) A_{mod,n} \f]
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
       *    \tilde{A}_{n+1} = A_{n+1}^{0} - \frac{1-\alpha_m}{\beta (\Delta t)^{2} (1-\alpha_f)}
       * D_{n+1}^{0} \f]
       *
       *  See the SetState() routine for the iterative update of the current state.
       *
       *  A_{mod,n} is updated in method UpdateStepState() according to
       *  \f[
       *    A_{mod,n+1} = \frac{1}{1-\alpha_m} ( (1-\alpha_f) A_{n+1} + \alpha_f A_n - \alpha_m
       * A_{mod,n} ) \f]
       *
       *  [derived]
       *
       *  \author grill
       *  \date 07/2016 */
      virtual void UpdateConstantStateContributions();

     private:
      /*! @name New vectors for internal use only
       *
       *  If an external use seems necessary, move these vectors to the
       *  global state data container and just store a pointer to the global
       *  state variable. */
      //! @{
      //! modified acceleration vector acc_{mod;n}
      Teuchos::RCP<Epetra_Vector> accn_mod_;

      //! @}
    };
  }  // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_GENALPHA__LIEGROUP_H_ */
