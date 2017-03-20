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
namespace STEEPESTASCENT {
  // forward declarations
  class Interface;

  /*--------------------------------------------------------------------------*/
  class DataContainer : public ::CONTACT::AUG::DataContainer
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

  protected:
    // don't want assignment operator and copy constructor
    DataContainer operator =(const DataContainer& old);
    DataContainer(const DataContainer& old);

  private:
    // old constraint violation from the last converged Newton step
    double old_infeasibility_measure_;

    // user provided upper bound for the cn value
    double cn_upper_bound_;

  };  // class DataContainer

  /*--------------------------------------------------------------------------*/
  class Strategy : public ::CONTACT::AUG::Strategy
  {
  public:
    /// constructor
    Strategy(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        int maxdof );

    /// destructor
    virtual ~Strategy() {};

  protected:
    /// derived
    virtual void AddContributionsToConstrRHS(
        Epetra_Vector& augConstrRhs ) const;

    /// derived
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const;

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
    virtual void RunPreComputeX(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable );

    void RunPostIterate(
        const CONTACT::ParamsInterface& cparams );

    /// derived
    virtual bool AssembleContactRHS();

  private:
    /** \brief Get access to the internal data container of the strategy (mutable)
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    CONTACT::AUG::STEEPESTASCENT::DataContainer& Data()
    {
      if (saDataPtr_.is_null())
        dserror("The steepest-ascent strategy data container is not initialized!");
      return *saDataPtr_;
    }

    /** \brief Get access to the internal data container of the strategy (read-only)
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    const CONTACT::AUG::STEEPESTASCENT::DataContainer& Data() const
    {
      if (saDataPtr_.is_null())
        dserror("The steepest-ascent strategy data container is not initialized!");
      return *saDataPtr_;
    }

    void UpdateCn( const CONTACT::ParamsInterface& cparams );

    void AugmentDirection(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable );

    Teuchos::RCP<Epetra_Vector> ComputeActiveLagrangeIncrInNormalDirection(
        const Epetra_Vector& displ_incr );

    Teuchos::RCP<Epetra_Vector> ComputeInactiveLagrangeIncrInNormalDirection(
        const Epetra_Vector& displ_incr, const Epetra_Vector& zold_inactive  );

    void PostAugmentDirection();

  private:
    Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::DataContainer> saDataPtr_;
    std::vector<Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::Interface> > interface_;

    typedef std::vector<Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::Interface> > InterfaceVector;
  };  //  class Strategy

} // namespace STEEPESTASCENT
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_STRATEGY_H_ */
