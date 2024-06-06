/*-----------------------------------------------------------*/
/*! \file

\brief Some special solution strategy

\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_STEEPEST_ASCENT_SP_STRATEGY_HPP
#define FOUR_C_CONTACT_AUG_STEEPEST_ASCENT_SP_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    class LagrangeMultiplierFunction;
    class PenaltyUpdate;
    namespace SteepestAscent
    {
      /*--------------------------------------------------------------------------*/
      class DataContainer
      {
       public:
        /// constructor
        DataContainer();

        /// destructor
        virtual ~DataContainer() = default;

        void set_penalty_correction_parameter(const double correction_param)
        {
          if (correction_param < 0.0)
            FOUR_C_THROW(
                "The correction parameter must be equal to 0.0 (default) or "
                "larger! A value of %.4e is not allowed!",
                correction_param);

          penalty_correction_parameter_ = correction_param;
        }

        double get_penalty_correction_parameter() const { return penalty_correction_parameter_; }

        void set_penalty_decrease_correction_parameter(const double correction_param)
        {
          if (correction_param < get_penalty_correction_parameter())
            FOUR_C_THROW(
                "The DECREASE CORRECTION PARAMETER must be larger than the"
                "penalty correction parameter [currently set to %.4e]! A value of "
                "%.4e is not allowed!",
                get_penalty_correction_parameter(), correction_param);

          penalty_decrease_correction_parameter_ = correction_param;
        }

        double get_penalty_decrease_correction_parameter() const
        {
          return penalty_decrease_correction_parameter_;
        }

        Teuchos::RCP<const LagrangeMultiplierFunction> lagrange_multiplier_func_ptr() const
        {
          return lm_func_ptr_.getConst();
        }

        Teuchos::RCP<LagrangeMultiplierFunction>& lagrange_multiplier_func_ptr()
        {
          return lm_func_ptr_;
        }

        LagrangeMultiplierFunction& lagrange_multiplier_func()
        {
          if (lm_func_ptr_.is_null())
            FOUR_C_THROW("The lm_func_ptr_ is not initialized correctly!");

          return *lm_func_ptr_;
        }

        Teuchos::RCP<const Aug::PenaltyUpdate> PenaltyUpdatePtr() const
        {
          return penalty_update_ptr_.getConst();
        }

        Teuchos::RCP<Aug::PenaltyUpdate>& PenaltyUpdatePtr() { return penalty_update_ptr_; }

        Aug::PenaltyUpdate& PenaltyUpdate()
        {
          if (penalty_update_ptr_.is_null())
            FOUR_C_THROW("The lm_func_ptr_ is not initialized correctly!");

          return *penalty_update_ptr_;
        }

        void SetStepLength(double steplength) { steplength_ = steplength; }

        double StepLength() const { return steplength_; }

       protected:
        // don't want assignment operator and copy constructor
        DataContainer operator=(const DataContainer& old);
        DataContainer(const DataContainer& old);

       private:
        /// correction control parameter for a possible increase of the cn value
        double penalty_correction_parameter_ = 0.0;

        /// correction control parameter for a possible decrease of the cn value
        double penalty_decrease_correction_parameter_ = 1.0;

        /// current step length
        double steplength_ = 1.0;

        /// Lagrange multiplier function pointer
        Teuchos::RCP<LagrangeMultiplierFunction> lm_func_ptr_ = Teuchos::null;

        /// pointer to the cn-correction object
        Teuchos::RCP<Aug::PenaltyUpdate> penalty_update_ptr_ = Teuchos::null;

      };  // class DataContainer
    }     // namespace SteepestAscent

    namespace SteepestAscentSaddlePoint
    {
      /*--------------------------------------------------------------------------*/
      /** \brief Saddle-point variant of the modified Newton approach.
       *
       * \author hiermeier \date 08/18 */
      class Strategy : public Aug::Strategy
      {
        /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
         *  strategies and needs access to all methods. */
        friend class CONTACT::Aug::ComboStrategy;

       public:
        /// constructor
        Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
            const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
            const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
            const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

        /// derived
        Inpar::CONTACT::SolvingStrategy Type() const override
        {
          return Inpar::CONTACT::solution_steepest_ascent_sp;
        }

       protected:
        /// scale factor similar to the std. Lagrangian strategy
        double inactive_scale_factor() const override { return 1.0; }

        /** \brief Assemble the structural contact rhs [derived]
         *
         *  In contradiction to the base class only the the Lagrange multiplier
         *  forces are considered.
         *
         *  \author hiermeier \date 03/17 */
        void eval_str_contact_rhs() override;

        /// derived
        void post_setup(bool redistributed, bool init) override;

        /// derived
        void run_post_apply_jacobian_inverse(const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
            const NOX::Nln::Group& grp) override;

        /// initiate the actual update of the cn-value (increase or decrease)
        void run_post_iterate(const CONTACT::ParamsInterface& cparams) override;

        /// call the penalty update object for a possible cn increase
        void update_cn(const CONTACT::ParamsInterface& cparams);

        /// call the penalty update object for a possible cn decrease
        void decrease_cn(const CONTACT::ParamsInterface& cparams);

        /// set current (default-step) solution state in the penalty update method
        void set_penalty_update_state(const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& xold, const Epetra_Vector& dir);

        /// used to add the modified diagonal entries
        void add_contributions_to_matrix_block_lm_lm(
            Core::LinAlg::SparseMatrix& kzz) const override;

        /// Returns the diagonal modification vector for the active dofs
        Teuchos::RCP<Epetra_Vector> get_kzz_diag_modification() const;
      };
    }  // namespace SteepestAscentSaddlePoint
  }    // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
