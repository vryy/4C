/*-----------------------------------------------------------*/
/*! \file

\brief Static (time) integrator.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_IMPL_STATICS_HPP
#define FOUR_C_STRUCTURE_NEW_IMPL_STATICS_HPP

#include "4C_config.hpp"

#include "4C_structure_new_impl_generic.hpp"

#include <NOX_Abstract_Vector.H>

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace IMPLICIT
  {
    class Statics : public Generic
    {
     public:
      //! constructor
      Statics();

      //! Setup the class variables
      void setup() override;

      //! (derived)
      void post_setup() override;

      //! Reset state variables (derived)
      void set_state(const Epetra_Vector& x) override;

      //! Apply the rhs only (derived)
      bool apply_force(const Epetra_Vector& x, Epetra_Vector& f) override;

      //! Apply the stiffness only (derived)
      bool apply_stiff(const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac) override;

      //! Apply force and stiff at once (derived)
      bool apply_force_stiff(
          const Epetra_Vector& x, Epetra_Vector& f, Core::LinAlg::SparseOperator& jac) override;

      //! (derived)
      bool assemble_force(
          Epetra_Vector& f, const std::vector<Inpar::Solid::ModelType>* without_these_models =
                                nullptr) const override;

      //! (derived)
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! (derived)
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! (derived)
      double calc_ref_norm_force(const enum ::NOX::Abstract::Vector::NormType& type) const override;

      //! return time integration factor (derived)
      [[nodiscard]] double get_int_param() const override;

      //! derived
      double get_model_value(const Epetra_Vector& x) override;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      void pre_update() override;

      //! Update configuration after time step (derived)
      void update_step_state() override;

      //! Update everything on element level after time step and after output (derived)
      void update_step_element() override;
      //! @}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //! @{
      /*! predict constant displacements, consistent velocities and accelerations (derived) */
      void predict_const_dis_consist_vel_acc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;

      /*! predict displacements based on constant velocities and consistent accelerations (derived)
       */
      bool predict_const_vel_consist_acc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;

      /*! predict displacements based on constant accelerations and consistent velocities (derived)
       */
      bool predict_const_acc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;
      //! @}

      //! @name Attribute access functions
      //@{

      //! Return name
      enum Inpar::Solid::DynamicType method_name() const override
      {
        return Inpar::Solid::dyna_statics;
      }

      //! Provide number of steps, a single-step method returns 1
      int method_steps() const override { return 1; }

      //! Give local order of accuracy of displacement part
      int method_order_of_accuracy_dis() const override
      {
        FOUR_C_THROW("Sensible to ask?");
        return 0;
      }

      //! Give local order of accuracy of velocity part
      int method_order_of_accuracy_vel() const override
      {
        FOUR_C_THROW("Sensible to ask?");
        return 0;
      }

      //! Return linear error coefficient of displacements
      double method_lin_err_coeff_dis() const override
      {
        FOUR_C_THROW("Sensible to ask?");
        return 0.0;
      }

      //! Return linear error coefficient of velocities
      double method_lin_err_coeff_vel() const override
      {
        FOUR_C_THROW("Sensible to ask?");
        return 0.0;
      }

      //@}

     protected:
      /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
       *
       * \remark Nothing needs to be done in the static case. */
      void add_visco_mass_contributions(Epetra_Vector& f) const override{};

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       *
       * \remark Nothing needs to be done in the static case. */
      void add_visco_mass_contributions(Core::LinAlg::SparseOperator& jac) const override{};

      //! reset the time step dependent parameters for the element evaluation [derived]
      void reset_eval_params() override;
    };
  }  // namespace IMPLICIT
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
