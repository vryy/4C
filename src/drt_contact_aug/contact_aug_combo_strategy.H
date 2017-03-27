/*---------------------------------------------------------------------*/
/*!
\file contact_aug_combo_strategy.H

\brief This strategy allows the combination of an arbitrary number of
       augmented contact solving strategies.

\level 3

\maintainer Michael Hiermeier

\date Mar 20, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_COMBO_STRATEGY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_COMBO_STRATEGY_H_

#include "contact_aug_utils.H"
#include "../drt_contact/contact_abstract_strategy.H"

//#define DEBUG_COMBO_STRATEGY

namespace CONTACT {
namespace AUG {
  class Strategy;
  class Interface;
  class DataContainer;

  class ComboStrategy : public CoAbstractStrategy
  {
  public:
    /// create a new combo strategy object
    static Teuchos::RCP<CoAbstractStrategy> Create(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data,
        const Epetra_Map* dof_row_map,
        const Epetra_Map* node_row_map,
        const Teuchos::ParameterList& params,
        const plain_interface_set& ref_interfaces,
        const int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        const int maxdof );

  private:
    static void CreateStrategyInterfaces(
        const enum INPAR::CONTACT::SolvingStrategy strat_type,
        const plain_interface_set& ref_interfaces,
        plain_interface_set& strat_interfaces );

  public:
    /// constructor
    ComboStrategy(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& stratData,
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const plain_strategy_set& strategies,
        const int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        const int maxdof );

    /// destructor
    virtual ~ComboStrategy() {};

    /* Hide the public derived methods as well, since we do not want to access
     * any of these methods directly.                         hiermeier 03/17 */
  protected:
    /// @name Derived public methods from the base class ( a.k.a. CoAbstractStrategy )
    /// @{

    virtual void Reset(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& dispnp,
        const Epetra_Vector& xnew );

    /// function wrapper
    virtual bool IsSaddlePointSystem() const;

    /// function wrapper
    virtual void ResetActiveSet();

    /// function wrapper
    virtual double ConstraintNorm() const;

    /// function wrapper
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis);

    /// function wrapper
    virtual bool ActiveSetConverged();

    /// function wrapper
    virtual int ActiveSetSteps();

    /// function wrapper
    virtual void UpdateActiveSet();

    /// function wrapper
    virtual void EvaluateRelMovPredict();

    /// function wrapper
    virtual bool ActiveSetSemiSmoothConverged() const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const;

    /// function wrapper
    virtual const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const;

    /// function wrapper
    virtual const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
            const enum DRT::UTILS::VecBlockType& bt) const;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr() const;

    /// function wrapper
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
            const enum DRT::UTILS::MatBlockType& bt) const;

    /// function wrapper
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetCondensedMatrixBlockPtr() const;

    /// function wrapper
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs();

    /// function wrapper
    virtual void Initialize();

    /// function wrapper
    virtual void EvalConstrRHS();

    /// function wrapper
    virtual void UpdateActiveSetSemiSmooth();

    /// @}

    /// @{

    /// @name accessors to public augmented Lagrangian methods
    /// @{

    /// function wrapper
    bool WasInContactLastIter() const;

    /// @}

  protected:
    /// @name Derived protected methods from the base class ( a.k.a. CoAbstractStrategy )
    /// @{

    /// Get the set of currently active interfaces
    virtual std::vector<Teuchos::RCP<CONTACT::CoInterface> >& Interfaces();

    /// Get the set of currently active interfaces
    virtual const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& Interfaces() const;

    /// function wrapper
    virtual void OutputStresses();

    /** \brief function wrapper: Redistribute and setup augmented Lagrangian members
     *
     *  \author hiermeier \date 03/17 */
    virtual void PostSetup( bool redistributed, bool init );

    /** \brief function wrapper: Compute force terms
     *
     *  \author hiermeier \date 03/17 */
    virtual void EvalForce(CONTACT::ParamsInterface& cparams);

    /** \brief function wrapper: Compute force and stiffness terms
     *
     *  \author hiermeier \date 03/17 */
    virtual void EvalForceStiff(CONTACT::ParamsInterface& cparams);

    /// function wrapper
    virtual void RunPostEvaluate( CONTACT::ParamsInterface& cparams );

    /// function wrapper
    virtual void RunPostIterate(
          const CONTACT::ParamsInterface& cparams );

    /// function wrapper
    virtual void RunPreComputeX(
          const CONTACT::ParamsInterface& cparams,
          const Epetra_Vector& xold,
          Epetra_Vector& dir_mutable );

    /// function wrapper
    virtual void RecoverState(
            const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& xold,
            const Epetra_Vector& dir,
            const Epetra_Vector& xnew);

    /// function wrapper
    virtual void ResetLagrangeMultipliers(
            const CONTACT::ParamsInterface& cparams,
            const Epetra_Vector& xnew);

    virtual void PostStoreDirichletStatus(
        Teuchos::RCP<const LINALG::MapExtractor> dbcmaps );

    /// @}
  private:

    /// access the currently active contact strategy
    CONTACT::AUG::Strategy& Get();

    /// access the currently active contact strategy (read-only)
    const CONTACT::AUG::Strategy& Get() const;

    void RunPostEvalForce( CONTACT::ParamsInterface& cparams );

    void RunPostEvalForceStiff( CONTACT::ParamsInterface& cparams );

    //! @name Unsupported derived routines (dead-end)
    //! @{

    //! @name Deprecated methods
    //! @{
    virtual void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Deprecated function call!"); };
    virtual void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Deprecated function call!"); };
    virtual void BuildSaddlePointSystem(
        Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd,
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
        int numiter,Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs)
    { dserror("Deprecated function call!"); };
    virtual void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<Epetra_Vector> blocksol)
    { dserror("Deprecated function call!"); };
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi)
    { dserror("Deprecated function call! Replaced by RecoverState()."); };
    //! @}

    /*! @name Dead-end for penalty and Uzawa methods (wrong strategy)
     *
     * Please note, that the definition of these functions seems completely unnecessary here.
     * Actually it would be a much better idea to cast the object to the right strategy at the
     * place where it is needed.                                            hiermeier 05/16 */
    //! @{
    virtual double InitialPenalty()
    { dserror("Wrong strategy!"); exit(EXIT_FAILURE); };
    virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Wrong strategy!"); };
    virtual void ResetPenalty()
    { dserror("Wrong strategy!"); };
    virtual void ModifyPenalty()
    { dserror("Wrong strategy!"); };
    virtual void UpdateUzawaAugmentedLagrange()
    { dserror("Wrong strategy!"); };
    virtual void UpdateConstraintNorm(int uzawaiter = 0)
    { dserror("Wrong strategy!"); };
    //! @}
    //! @}

  private:
    class Switching;


  private:
    /// wrapped contact solution strategies
    plain_strategy_set strategies_;

    /// interface of the wrapped contact solution strategies
    plain_interface_sets interface_sets_;

    CONTACT::AUG::DataContainer& data_;

    enum MORTAR::ActionType eval_action_type_;

    bool is_assemblecontactrhs_;

    struct GlobalNoDbc
    {
      /// constructor
      GlobalNoDbc()
          : slMaMap_( Teuchos::null ),
            slMap_( Teuchos::null ),
            maMap_( Teuchos::null ),
            slMaForce_( Teuchos::null ),
            slForce_( Teuchos::null ),
            maForce_( Teuchos::null )
      { /*do nothing*/ }

      void Assemble( const Epetra_Map& dbcmap,
          const CONTACT::AUG::DataContainer& data );

      void Redistribute(  const CONTACT::AUG::DataContainer& data );

      void Reset( const Epetra_Map& slMaMap,
          const CONTACT::AUG::DataContainer& data );

      Teuchos::RCP<Epetra_Map> slMaMap_;
      Teuchos::RCP<Epetra_Map> slMap_;
      Teuchos::RCP<Epetra_Map> maMap_;

      Teuchos::RCP<Epetra_Vector> slMaForce_;
      Teuchos::RCP<Epetra_Vector> slForce_;
      Teuchos::RCP<Epetra_Vector> maForce_;
    };

    GlobalNoDbc no_dbc_;

    Teuchos::RCP<Switching> switch_;
  };

  /*--------------------------------------------------------------------------*/
  class ComboStrategy::Switching
  {
  public:
    static Teuchos::RCP<Switching> Create( ComboStrategy& combo );

    /// constructor
    Switching( ComboStrategy& combo,
        const Teuchos::ParameterList& p_combo );

    /// destructor
    virtual ~Switching() {};

    /// return the active ID
    unsigned Id() const;

    /// return the ID of the given solving strategy
    unsigned Id( enum INPAR::CONTACT::SolvingStrategy sol_type ) const;

    /// check the switching condition and update the member variables
    void Update( CONTACT::ParamsInterface& cparams );

  private:
    unsigned FindId( INPAR::CONTACT::SolvingStrategy sol_type ) const;

    /// detect strategy types of the wrapped strategies
    void GetStrategyTypes(
        const plain_strategy_set& strategies,
        plain_strattype_set& strat_types ) const;

    bool CheckPenetration() const;

    bool CheckResidual(
        CONTACT::ParamsInterface& cparams );

    double GetPenetrationBound() const;

    Teuchos::RCP<Epetra_Vector> GetStructuralForceWithoutDbcDofs(
        const CONTACT::ParamsInterface& cparams );

    bool CheckContactResidualNorm(
        const Epetra_Vector& str_slmaforce,
        const Epetra_Vector& constr_slmaforce ) const;

    bool CheckAngleBetweenStrForceAndContactForce(
        const Epetra_Vector& str_slmaforce,
        const Epetra_Vector& constr_slmaforce ) const;

    void GetActiveSlMaForces(
        const Epetra_Vector& str_force,
        Teuchos::RCP<Epetra_Vector>& str_slmaforce,
        Teuchos::RCP<Epetra_Vector>& constr_slmaforce ) const;

    void GetGlobalSlMaActiveForceMaps(
        const Epetra_Vector& slforce,
        const Epetra_Vector& maforce,
        Teuchos::RCP<Epetra_Map>& gSlActiveForceMap,
        Teuchos::RCP<Epetra_Map>& gMaActiveForceMap ) const;

  private:
    /// reference to the combo strategy
    ComboStrategy& combo_;

    /// id of the pre-asymptotic solution strategy
    unsigned preasymptotic_id_;

    /// id of the asymptotic solution strategy
    unsigned asymptotic_id_;

    /// which phase is currently active?
    bool is_asymptotic_;

    plain_strattype_set strat_types_;
  };  // class SwitchingStrategy
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_COMBO_STRATEGY_H_ */
