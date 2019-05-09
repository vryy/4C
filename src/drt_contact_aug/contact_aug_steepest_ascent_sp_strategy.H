/*-----------------------------------------------------------*/
/*!
\maintainer Matthias Mayr

\brief Some special solution strategy

\level 3
*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_SP_STRATEGY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_SP_STRATEGY_H_

#include "contact_augmented_strategy.H"

namespace CONTACT
{
  namespace AUG
  {
    class LagrangeMultiplierFunction;
    class PenaltyUpdate;
    namespace STEEPESTASCENT
    {
      /*--------------------------------------------------------------------------*/
      class DataContainer
      {
       public:
        /// constructor
        DataContainer();

        /// destructor
        virtual ~DataContainer(){};

        void SetPenaltyCorrectionParameter(const double correction_param)
        {
          if (correction_param < 0.0)
            dserror(
                "The correction parameter must be equal to 0.0 (default) or "
                "larger! A value of %.4e is not allowed!",
                correction_param);

          penalty_correction_parameter_ = correction_param;
        }

        double GetPenaltyCorrectionParameter() const { return penalty_correction_parameter_; }

        void SetPenaltyDecreaseCorrectionParameter(const double correction_param)
        {
          if (correction_param < GetPenaltyCorrectionParameter())
            dserror(
                "The DECREASE CORRECTION PARAMETER must be larger than the"
                "penalty correction parameter [currently set to %.4e]! A value of "
                "%.4e is not allowed!",
                GetPenaltyCorrectionParameter(), correction_param);

          penalty_decrease_correction_parameter_ = correction_param;
        }

        double GetPenaltyDecreaseCorrectionParameter() const
        {
          return penalty_decrease_correction_parameter_;
        }

        Teuchos::RCP<const LagrangeMultiplierFunction> LagrangeMultiplierFuncPtr() const
        {
          return lm_func_ptr_.getConst();
        }

        Teuchos::RCP<LagrangeMultiplierFunction>& LagrangeMultiplierFuncPtr()
        {
          return lm_func_ptr_;
        }

        LagrangeMultiplierFunction& LagrangeMultiplierFunc()
        {
          if (lm_func_ptr_.is_null()) dserror("The lm_func_ptr_ is not initialized correctly!");

          return *lm_func_ptr_;
        }

        Teuchos::RCP<const AUG::PenaltyUpdate> PenaltyUpdatePtr() const
        {
          return penalty_update_ptr_.getConst();
        }

        Teuchos::RCP<AUG::PenaltyUpdate>& PenaltyUpdatePtr() { return penalty_update_ptr_; }

        AUG::PenaltyUpdate& PenaltyUpdate()
        {
          if (penalty_update_ptr_.is_null())
            dserror("The lm_func_ptr_ is not initialized correctly!");

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
        Teuchos::RCP<AUG::PenaltyUpdate> penalty_update_ptr_ = Teuchos::null;

      };  // class DataContainer
    }     // namespace STEEPESTASCENT

    namespace STEEPESTASCENT_SP
    {
      /*--------------------------------------------------------------------------*/
      /** \brief Saddle-point variant of the modified Newton approach.
       *
       * \author hiermeier \date 08/18 */
      class Strategy : public AUG::Strategy
      {
        /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
         *  strategies and needs access to all methods. */
        friend class CONTACT::AUG::ComboStrategy;

       public:
        /// constructor
        Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
            const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
            const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
            const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

        /// destructor
        virtual ~Strategy(){};

        /// derived
        INPAR::CONTACT::SolvingStrategy Type() const override
        {
          return INPAR::CONTACT::solution_steepest_ascent_sp;
        }

       protected:
        /// scale factor similar to the std. Lagrangian strategy
        double InactiveScaleFactor() const override { return 1.0; }

        /** \brief Assemble the structural contact rhs [derived]
         *
         *  In contradiction to the base class only the the Lagrange multiplier
         *  forces are considered.
         *
         *  \author hiermeier \date 03/17 */
        void EvalStrContactRHS() override;

        /// derived
        void PostSetup(bool redistributed, bool init) override;

        /// derived
        void RunPostApplyJacobianInverse(const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
            const NOX::NLN::Group& grp) override;

        /// initiate the actual update of the cn-value (increase or decrease)
        void RunPostIterate(const CONTACT::ParamsInterface& cparams) override;

        /// call the penalty update object for a possible cn increase
        void UpdateCn(const CONTACT::ParamsInterface& cparams);

        /// call the penalty update object for a possible cn decrease
        void DecreaseCn(const CONTACT::ParamsInterface& cparams);

        /// set current (default-step) solution state in the penalty update method
        void SetPenaltyUpdateState(const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& xold, const Epetra_Vector& dir);

        /// used to add the modified diagonal entries
        void AddContributionsToMatrixBlockLmLm(LINALG::SparseMatrix& kzz) const override;

        /// Returns the diagonal modification vector for the active dofs
        Teuchos::RCP<Epetra_Vector> GetKzzDiagModification() const;
      };
    }  // namespace STEEPESTASCENT_SP
  }    // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_SP_STRATEGY_H_ */
