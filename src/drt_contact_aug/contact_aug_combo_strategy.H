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

    /// return the type of currently active wrapped contact strategy
    INPAR::CONTACT::SolvingStrategy Type() const override;

    void Reset(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& dispnp,
        const Epetra_Vector& xnew ) override;

    /// function wrapper
    bool IsSaddlePointSystem() const override;

    /// function wrapper
    void ResetActiveSet() override;

    /// function wrapper
    double ConstraintNorm() const override;

    /// function wrapper
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override;

    /// function wrapper
    bool ActiveSetConverged() override;

    /// function wrapper
    int ActiveSetSteps() override;

    /// function wrapper
    void UpdateActiveSet() override;

    /// function wrapper
    void EvaluateRelMovPredict() override;

    /// function wrapper
    bool ActiveSetSemiSmoothConverged() const override;

    /// function wrapper
    Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const override;

    /// function wrapper
    Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override;

    /// function wrapper
    Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const override;

    /// function wrapper
    const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const override;

    /// function wrapper
    Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const override;

    /// function wrapper
    const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const override;

    /// function wrapper
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
            const enum DRT::UTILS::VecBlockType& bt) const;

    /// function wrapper
    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtrForNormCheck(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    /// function wrapper
    Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr(
        Epetra_Vector& f, const double & timefac_np ) const override;

    /// function wrapper
    Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams=NULL ) const override;

    /// function wrapper
    Teuchos::RCP<LINALG::SparseMatrix> GetCondensedMatrixBlockPtr(
        Teuchos::RCP<LINALG::SparseMatrix>& kteff,
        const double & timefac_np ) const override;

    /// function wrapper
    Teuchos::RCP<Epetra_Vector> ConstrRhs() override;

    /// function wrapper
    void Initialize() override;

    /// function wrapper
    void EvalConstrRHS() override;

    /// function wrapper
    void UpdateActiveSetSemiSmooth() override;

    /// function wrapper
    void DoReadRestart(
        IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis,
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr) override;

    /// function wrapper
    void Update( Teuchos::RCP<const Epetra_Vector> dis ) override;

    /// function wrapper
    double GetPotentialValue(
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type ) const override;

    /// function wrapper
    double GetLinearizedPotentialValueTerms(
        const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype ) const override;

    /// function wrapper
    void WriteOutput( IO::DiscretizationWriter& writer ) const override;

    /// function wrapper
    void EvaluateReferenceState(
        Teuchos::RCP<const Epetra_Vector> vec) override;

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
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& Interfaces() override;

    /// Get the set of currently active interfaces
    const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& Interfaces() const override;

    /// function wrapper
    void OutputStresses() override;

    /** \brief function wrapper: Redistribute and setup augmented Lagrangian members
     *
     *  \author hiermeier \date 03/17 */
    void PostSetup( bool redistributed, bool init ) override;

    /** \brief function wrapper: Compute force terms
     *
     *  \author hiermeier \date 03/17 */
    void EvalForce(CONTACT::ParamsInterface& cparams) override;

    /** \brief function wrapper: Compute force and stiffness terms
     *
     *  \author hiermeier \date 03/17 */
    void EvalForceStiff(CONTACT::ParamsInterface& cparams) override;

    /// function wrapper
    void EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams) override;

    /// function wrapper
    void RunPreEvaluate( CONTACT::ParamsInterface& cparams ) override;

    /// function wrapper
    void RunPostEvaluate( CONTACT::ParamsInterface& cparams ) override;

    /// function wrapper
    void RunPostIterate(
        const CONTACT::ParamsInterface& cparams ) override;

    /// function wrapper
    void RunPreComputeX(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable ) override;

    /// function wrapper
    void RunPostComputeX(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew) override;

    /// function wrapper
    void RunPostApplyJacobianInverse(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& rhs,
        Epetra_Vector& result,
        const Epetra_Vector& xold,
        const NOX::NLN::Group& grp ) override;

    /// function wrapper
    void RemoveCondensedContributionsFromRhs(
        Epetra_Vector& str_rhs ) const override;

    /// function wrapper
    void CorrectParameters(
        CONTACT::ParamsInterface& cparams,
        const NOX::NLN::CorrectionType type ) override;

    /// function wrapper
    void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xnew) override;

    void PostStoreDirichletStatus(
        Teuchos::RCP<const LINALG::MapExtractor> dbcmaps ) override;

    /// @}
  private:
    /// @name class member functions
    /// @{

    /// access the currently active contact strategy
    CONTACT::AUG::Strategy& Get();

    /// access the currently active contact strategy (read-only)
    const CONTACT::AUG::Strategy& Get() const;

    void SwitchUpdate( CONTACT::ParamsInterface& cparams );

    void RunPostEvalForce( CONTACT::ParamsInterface& cparams );

    void RunPostEvalForceStiff( CONTACT::ParamsInterface& cparams );

    void RunPostEvalStaticConstraintRHS( CONTACT::ParamsInterface& cparams );

    /// @}

    //! @name Unsupported derived routines (dead-end)
    //! @{

    //! @name Deprecated methods
    //! @{
    void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    { dserror("Deprecated function call!"); };
    void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override
    { dserror("Deprecated function call!"); };
    void BuildSaddlePointSystem(
        Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd,
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
        int numiter,Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs) override
    { dserror("Deprecated function call!"); };
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<Epetra_Vector> blocksol) override
    { dserror("Deprecated function call!"); };
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override
    { dserror("Deprecated function call! Replaced by RunPostComputeX()."); };
    //! @}

    /*! @name Dead-end for penalty and Uzawa methods (wrong strategy)
     *
     * Please note, that the definition of these functions seems completely unnecessary here.
     * Actually it would be a much better idea to cast the object to the right strategy at the
     * place where it is needed.                                            hiermeier 05/16 */
    //! @{
    double InitialPenalty() override
    { dserror("Wrong strategy!"); exit(EXIT_FAILURE); };
    void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) override
    { dserror("Wrong strategy!"); };
    void ResetPenalty() override
    { dserror("Wrong strategy!"); };
    void ModifyPenalty() override
    { dserror("Wrong strategy!"); };
    void UpdateUzawaAugmentedLagrange() override
    { dserror("Wrong strategy!"); };
    void UpdateConstraintNorm(int uzawaiter = 0) override
    { dserror("Wrong strategy!"); };
    //! @}
    //! @}

  private:
    // nested classes
    class Switching;
    class PreAsymptoticSwitching;

  private:
    /// wrapped contact solution strategies
    plain_strategy_set strategies_;

    /// interface of the wrapped contact solution strategies
    plain_interface_sets interface_sets_;

    CONTACT::AUG::DataContainer& data_;

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

    class Output
    {
    public:
      Output() = default;
      ~Output() = default;

      void initScreenOutput( bool print2screen );

      std::ostream& oscreen() const;

    private:
      std::ostream* oscreen_ = NULL;

      const Teuchos::RCP<Teuchos::oblackholestream> blackhole_ =
          Teuchos::rcp( new Teuchos::oblackholestream );
    };

    GlobalNoDbc no_dbc_;
    Output output_;

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

    /// return the ID of the given solving strategy
    unsigned Id( enum INPAR::CONTACT::SolvingStrategy sol_type ) const;

    /// return the active ID
    virtual unsigned Id() const = 0;

    /// check the switching condition and update the member variables
    virtual void Update( CONTACT::ParamsInterface& cparams,
        std::ostream& os ) = 0;

  protected:
    unsigned FindId( INPAR::CONTACT::SolvingStrategy sol_type ) const;

    /// detect strategy types of the wrapped strategies
    void GetStrategyTypes(
        const plain_strategy_set& strategies,
        plain_strattype_set& strat_types ) const;

  protected:
    /// reference to the combo strategy
    ComboStrategy& combo_;

    plain_strattype_set strat_types_;
  };  // class SwitchingStrategy

  /*--------------------------------------------------------------------------*/
  class ComboStrategy::PreAsymptoticSwitching : public ComboStrategy::Switching
  {
  public:
    /// constructor
    PreAsymptoticSwitching( ComboStrategy& combo,
        const Teuchos::ParameterList& p_combo );

    /// destructor
    virtual ~PreAsymptoticSwitching() {};

  private:
    /// derived
    unsigned Id() const override;

    /// derived
    void Update( CONTACT::ParamsInterface& cparams,
        std::ostream& os ) override;

    bool CheckPenetration( std::ostream& os );

    bool CheckResidual(
        CONTACT::ParamsInterface& cparams,
        std::ostream& os );

    double GetPenetrationBound() const;

    Teuchos::RCP<Epetra_Vector> GetStructuralForceWithoutDbcDofs(
        const CONTACT::ParamsInterface& cparams );

    bool CheckContactResidualNorm(
        const Epetra_Vector& str_slmaforce,
        const Epetra_Vector& constr_slmaforce,
        std::ostream& os ) const;

    bool CheckAngleBetweenStrForceAndContactForce(
        const Epetra_Vector& str_slmaforce,
        const Epetra_Vector& constr_slmaforce,
        std::ostream& os ) const;

    void GetActiveSlMaForces(
        const Epetra_Vector& str_force,
        Teuchos::RCP<Epetra_Vector>& str_slmaforce,
        Teuchos::RCP<Epetra_Vector>& constr_slmaforce ) const;

    void GetGlobalSlMaActiveForceMaps(
        const Epetra_Vector& slforce,
        const Epetra_Vector& maforce,
        Teuchos::RCP<Epetra_Map>& gSlActiveForceMap,
        Teuchos::RCP<Epetra_Map>& gMaActiveForceMap ) const;

    void PrintUpdateHead( std::ostream& os ) const;

  private:
    struct MinAWGap
    {
      MinAWGap( PreAsymptoticSwitching& switching )
          : switch_( switching ),
            asymptotic_( 1.0e12 ),
            pre_asymptotic_( 1.012 )
      {};

      inline void Update( const double min_awgap )
      {
        if ( switch_.is_asymptotic_ )
          asymptotic_ = min_awgap;
        else
          pre_asymptotic_ = min_awgap;
      }

      const PreAsymptoticSwitching& switch_;
      double asymptotic_;
      double pre_asymptotic_;
    };

    /// id of the pre-asymptotic solution strategy
    unsigned preasymptotic_id_;

    /// id of the asymptotic solution strategy
    unsigned asymptotic_id_;

    /// which phase is currently active?
    bool is_asymptotic_;

    MinAWGap minawgap_;
  };
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_COMBO_STRATEGY_H_ */
