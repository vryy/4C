/*---------------------------------------------------------------------*/
/*!
\file contact_aug_steepest_ascent_strategy.H

\brief Steepest ascent solution strategy based on the augmented contact
       formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 7, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_

#include "contact_augmented_strategy.H"

namespace CONTACT {
namespace AUG {
  class LagrangeMultiplierFunction;
  class PenaltyUpdate;
namespace STEEPESTASCENT {
  // forward declarations
  class Interface;

  /*--------------------------------------------------------------------------*/
  class DataContainer
  {
  public:
    /// constructor
    DataContainer();

    /// destructor
    virtual ~DataContainer() {};

    void SetOldInfeasibilityMeasure( const double old_infeasibility_measure )
    {
      old_infeasibility_measure_ = old_infeasibility_measure;
    }

    double GetOldInfeasibilityMeasure() const
    {
      return old_infeasibility_measure_;
    }

    void SetCnUpperBound( const double cn_upper_bound )
    {
      cn_upper_bound_ = cn_upper_bound;
    }

    double GetCnUpperBound() const
    {
      return cn_upper_bound_;
    }

    void SetPenaltyCorrectionParameter( const double correction_param )
    {
      if ( correction_param < 0.0 )
        dserror("The correction parameter must be equal to 0.0 (default) or "
            "larger! A value of %.4e is not allowed!", correction_param );

      penalty_correction_parameter_ = correction_param;
    }

    double GetPenaltyCorrectionParameter() const
    {
      return penalty_correction_parameter_;
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
      if ( lm_func_ptr_.is_null() )
        dserror( "The lm_func_ptr_ is not initialized correctly!" );

      return *lm_func_ptr_;
    }

    Teuchos::RCP<const AUG::PenaltyUpdate> PenaltyUpdatePtr() const
    {
      return penalty_update_ptr_.getConst();
    }

    Teuchos::RCP<AUG::PenaltyUpdate>& PenaltyUpdatePtr()
    {
      return penalty_update_ptr_;
    }

    AUG::PenaltyUpdate& PenaltyUpdate()
    {
      if ( penalty_update_ptr_.is_null() )
        dserror( "The lm_func_ptr_ is not initialized correctly!" );

      return *penalty_update_ptr_;
    }

    void SetStepLength( double steplength )
    {
      steplength_ = steplength;
    }

    double StepLength() const
    {
      return steplength_;
    }

  protected:
    // don't want assignment operator and copy constructor
    DataContainer operator =(const DataContainer& old);
    DataContainer(const DataContainer& old);

  private:
    // old constraint violation from the last converged Newton step
    double old_infeasibility_measure_ = 0.0;

    // user provided upper bound for the cn value
    double cn_upper_bound_ = 0.0;

    double penalty_correction_parameter_ = 0.0;

    double steplength_ = 1.0;

    Teuchos::RCP<LagrangeMultiplierFunction> lm_func_ptr_ = Teuchos::null;

    Teuchos::RCP<AUG::PenaltyUpdate> penalty_update_ptr_ = Teuchos::null;

  };  // class DataContainer

  /*--------------------------------------------------------------------------*/
  class Strategy : public ::CONTACT::AUG::Strategy
  {
    /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
     *  strategies and needs access to all methods. */
    friend class CONTACT::AUG::ComboStrategy;

  public:
    /// constructor
    Strategy(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const plain_interface_set& interfaces,
        int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        int maxdof );

    /// destructor
    virtual ~Strategy() {};

    INPAR::CONTACT::SolvingStrategy Type() const override
    {
      return INPAR::CONTACT::solution_steepest_ascent;
    }

  protected:

    /// derived
    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtrForNormCheck(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    /// derived
    void PostSetup(
        bool redistributed,
        bool init) override;

    /// derived
    void AddContributionsToConstrRHS(
        Epetra_Vector& augConstrRhs ) const override;

    /// derived
    Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = NULL) const override;

    /// derived
    void AddContributionsToMatrixBlockDisplDispl(
        LINALG::SparseMatrix& kdd,
        const CONTACT::ParamsInterface* cparams = NULL ) const override;

    /*! \brief Derived function, which is called at the very beginning of a call to
     *  NOX::NLN::Group::computeX()
     *
     *  This method is used to get access to the current direction vector and
     *  to augment/modify the direction vector before the actual state update is
     *  performed. One possible scenario is the steepest ascent method, where
     *  we calculate the Lagrange multiplier increment in a post-processing
     *  step.
     *
     *  \author hiermeier \date 03/17 */
    void RunPreComputeX(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable ) override;

    void RunPostApplyJacobianInverse(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& rhs,
        Epetra_Vector& result,
        const Epetra_Vector& xold,
        const NOX::NLN::Group& grp ) override;

    void RunPostIterate(
        const CONTACT::ParamsInterface& cparams ) override;

    /** \brief parameter correction during a more sophisticated non-linear
     *  solver approach (e.g. the filter method) [derived]
     *
     *  \author hiermeier \date 12/17 */
    void CorrectParameters(
        CONTACT::ParamsInterface& cparams,
        const NOX::NLN::CorrectionType type ) override;

  private:

    void UpdateCn( const CONTACT::ParamsInterface& cparams );

    void AugmentDirection(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable );

    Teuchos::RCP<Epetra_Vector> ComputeActiveLagrangeIncrInNormalDirection(
        const Epetra_Vector& displ_incr );

    Teuchos::RCP<Epetra_Vector> ComputeInactiveLagrangeIncrInNormalDirection(
        const Epetra_Vector& displ_incr, const Epetra_Vector& zold );

    void PostAugmentDirection(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir );
  };  //  class Strategy

} // namespace STEEPESTASCENT
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_ */
