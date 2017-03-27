/*---------------------------------------------------------------------*/
/*!
\file contact_aug_combo_strategy.cpp

\brief This strategy allows the combination of an arbitrary number of
       augmented contact solving strategies.

\level 3

\maintainer Michael Hiermeier

\date Mar 20, 2017

*/
/*---------------------------------------------------------------------*/

#include "contact_aug_combo_strategy.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_strategy_factory.H"
#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_inpar/inpar_contact.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoAbstractStrategy> CONTACT::AUG::ComboStrategy::Create(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data,
    const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map,
    const Teuchos::ParameterList& params,
    const plain_interface_set& ref_interfaces,
    const int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm,
    const int maxdof )
{
  const Teuchos::ParameterList& p_combo =
      params.sublist( "AUGMENTED" ).sublist( "COMBO" );

  plain_strategy_set strategies(0);
  plain_interface_set strat_interfaces(0);

  unsigned count = 0;
  std::ostringstream strat_count;
  strat_count << "STRATEGY_" << count;
  while ( p_combo.isParameter( strat_count.str() ) )
  {
    const enum INPAR::CONTACT::SolvingStrategy strat_type =
        DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(
            p_combo, strat_count.str() );

    CreateStrategyInterfaces( strat_type, ref_interfaces, strat_interfaces );

    strategies.push_back( STRATEGY::Factory::BuildStrategy(
        strat_type, params, false, false, maxdof, strat_interfaces, dof_row_map,
        node_row_map, dim, comm, data ) );

    /// clear and increase strategy count string
    strat_count.str("");
    strat_count << "STRATEGY_" << ++count;
  }

  return Teuchos::rcp( new ComboStrategy( data, dof_row_map, node_row_map,
      params, strategies, dim, comm, maxdof ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::CreateStrategyInterfaces(
    const enum INPAR::CONTACT::SolvingStrategy strat_type,
    const plain_interface_set& ref_interfaces,
    plain_interface_set& strat_interfaces )
{
  strat_interfaces.clear();
  strat_interfaces.reserve( ref_interfaces.size() );

  Teuchos::RCP<CONTACT::CoInterface> newinterface = Teuchos::null;

  for ( plain_interface_set::const_iterator cit = ref_interfaces.begin();
        cit != ref_interfaces.end(); ++cit )
  {
    const CONTACT::AUG::Interface& ref_interface =
       dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    const Teuchos::RCP<CONTACT::AUG::IDataContainer> idata_ptr =
        ref_interface.SharedInterfaceDataPtr();
    CONTACT::AUG::IDataContainer& idata = *idata_ptr;

    /* create a new interface by copying the data pointer from the reference
     * interface */
    newinterface = STRATEGY::Factory::CreateInterface( strat_type, idata.Id(),
        idata.Comm(), idata.Dim(), idata.IMortar(), idata.IsSelfContact(),
        idata.RedundantStorage(), Teuchos::null, idata_ptr );

    strat_interfaces.push_back( newinterface );
    // reset pointer
    newinterface = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::ComboStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data,
    const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map,
    const Teuchos::ParameterList& params,
    const plain_strategy_set& strategies,
    const int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm,
    const int maxdof )
    : CONTACT::CoAbstractStrategy( data, dof_row_map, node_row_map,
        params, dim, comm, 0.0, maxdof ),
        strategies_( strategies ),
        interface_sets_( 0 ),
        data_( dynamic_cast<CONTACT::AUG::DataContainer&>( *data ) ),
        eval_action_type_( MORTAR::eval_none ),
        is_assemblecontactrhs_( true ),
        no_dbc_(),
        switch_( Switching::Create( *this ) )
{
  for ( plain_strategy_set::const_iterator cit = strategies_.begin();
        cit != strategies_.end(); ++cit )
  {
    CONTACT::CoAbstractStrategy& s = **cit;
    const plain_interface_set& sinterfaces = s.ContactInterfaces();
    const unsigned num_interfaces = sinterfaces.size();

    interface_sets_.push_back( plain_interface_set( num_interfaces, Teuchos::null ) );
    plain_interface_set& interfaces = interface_sets_.back();
    std::copy( sinterfaces.begin(), sinterfaces.end(), interfaces.begin() );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvaluate(
    CONTACT::ParamsInterface& cparams )
{
  switch ( eval_action_type_ )
  {
    case MORTAR::eval_force:
    {
      RunPostEvalForce( cparams );
      break;
    }
    case MORTAR::eval_force_stiff:
    {
      RunPostEvalForceStiff( cparams );
      break;
    }
    default:
      dserror( "Unexpected MORTAR::ActionType! ( actiontype = %s | %d )",
          MORTAR::ActionType2String( eval_action_type_ ).c_str(),
          eval_action_type_ );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Reset(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& dispnp,
    const Epetra_Vector& xnew )
{
  eval_action_type_ = MORTAR::eval_none;

  Get().Reset( cparams, dispnp, xnew );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::IsSaddlePointSystem() const
{
  return Get().IsSaddlePointSystem();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::ResetActiveSet()
{
  Get().ResetActiveSet();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::SaveReferenceState(
    Teuchos::RCP<const Epetra_Vector> dis )
{
  Get().SaveReferenceState( dis );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::ConstraintNorm() const
{
  return Get().ConstraintNorm();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::ActiveSetConverged()
{
  return Get().ActiveSetConverged();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::ComboStrategy::ActiveSetSteps()
{
  return Get().ActiveSetSteps();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::UpdateActiveSet()
{
  return Get().UpdateActiveSet();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvaluateRelMovPredict()
{
  Get().EvaluateRelMovPredict();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::ActiveSetSemiSmoothConverged() const
{
  return Get().ActiveSetSemiSmoothConverged();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::ComboStrategy::GetOldActiveRowNodes() const
{
  return Get().GetOldActiveRowNodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::ComboStrategy::GetOldSlipRowNodes() const
{
  return Get().GetOldSlipRowNodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::ComboStrategy::SlNormalDoFRowMapPtr(const bool& redist) const
{
  return Get().SlNormalDoFRowMapPtr( redist );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map&
CONTACT::AUG::ComboStrategy::SlNormalDoFRowMap(const bool& redist) const
{
  return Get().SlNormalDoFRowMap( redist );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::ComboStrategy::SlTangentialDoFRowMapPtr(const bool& redist) const
{
  return Get().SlTangentialDoFRowMapPtr( redist );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map&
CONTACT::AUG::ComboStrategy::SlTangentialDoFRowMap(const bool& redist) const
{
  return Get().SlTangentialDoFRowMap( redist );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  if ( not is_assemblecontactrhs_ )
    return Teuchos::null;

  return Get().GetRhsBlockPtr( bt );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::ComboStrategy::GetCondensedRhsPtr() const
{
  if ( not is_assemblecontactrhs_ )
    return Teuchos::null;

  return Get().GetCondensedRhsPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::ComboStrategy::GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const
{
  return Get().GetMatrixBlockPtr( bt );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
CONTACT::AUG::ComboStrategy::GetCondensedMatrixBlockPtr() const
{
  return Get().GetCondensedMatrixBlockPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::ComboStrategy::ConstrRhs()
{
  return Get().ConstrRhs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Initialize()
{
  dserror( "Unnecessary in this Strategy." );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalConstrRHS()
{
  dserror( "Unnecessary in this Strategy." );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::UpdateActiveSetSemiSmooth()
{
  dserror( "Unnecessary in this Strategy." );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::WasInContactLastIter() const
{
  return Get().WasInContactLastIter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::OutputStresses()
{
  Get().OutputStresses();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PostSetup( bool redistributed, bool init )
{
  if ( redistributed )
    no_dbc_.Redistribute( data_ );

  Get().PostSetup( redistributed, init );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PostStoreDirichletStatus(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmaps )
{
  no_dbc_.Assemble( *dbcmaps->CondMap(), data_ );

  Get().PostStoreDirichletStatus( dbcmaps );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Assemble(
    const Epetra_Map& dbcmap, const CONTACT::AUG::DataContainer& data )
{
  const Epetra_Map& gSlMaDofRowMap = *data.GSlMaDofRowMapPtr();

  Teuchos::RCP<Epetra_Map> gSlMaDbcDofRowMap =
      LINALG::IntersectMap( gSlMaDofRowMap, dbcmap );

  slMaMap_ = LINALG::SplitMap( gSlMaDofRowMap, *gSlMaDbcDofRowMap );

  Reset( *slMaMap_, data );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Redistribute(
    const CONTACT::AUG::DataContainer& data )
{
  RedistributeRowMap( *data.GSlMaDofRowMapPtr() , *slMaMap_ );

  Reset( *slMaMap_, data );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Reset(
    const Epetra_Map& slMaMap, const CONTACT::AUG::DataContainer& data )
{
  slMap_ = LINALG::IntersectMap( slMaMap, *data.GSlDofRowMapPtr() );
  maMap_ = LINALG::IntersectMap( slMaMap, *data.GMaDofRowMapPtr() );

  slMaForce_ = Teuchos::rcp( new Epetra_Vector( slMaMap, true ) );
  slForce_ = Teuchos::rcp( new Epetra_Vector( *slMap_, true ) );
  maForce_ = Teuchos::rcp( new Epetra_Vector( *maMap_, true ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  eval_action_type_ = cparams.GetActionType();

  Get().PreEvalForce( cparams );

  Get().AssembleGap();

  Get().EvalAugmentedForces();

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << INPAR::CONTACT::SolvingStrategy2String( Get().Type() ) << "\n";
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvalForce(
    CONTACT::ParamsInterface& cparams )
{
  switch_->Update( cparams );

  Get().AssembleContactRHS();
  Get().EvalStrContactRHS();
  Get().EvalConstrRHS();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForceStiff(
    CONTACT::ParamsInterface& cparams)
{
  EvalForce( cparams );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvalForceStiff(
    CONTACT::ParamsInterface& cparams )
{
  RunPostEvalForce( cparams );

  Get().AssembleContactStiff();
  Get().PostEvalForceStiff( cparams );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostIterate(
    const CONTACT::ParamsInterface& cparams )
{
  Get().RunPostIterate( cparams );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreComputeX(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    Epetra_Vector& dir_mutable )
{
  Get().RunPreComputeX( cparams, xold, dir_mutable );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RecoverState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  Get().RecoverState( cparams, xold, dir, xnew );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xnew)
{
  Get().ResetLagrangeMultipliers( cparams, xnew );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CONTACT::CoInterface> >&
CONTACT::AUG::ComboStrategy::Interfaces()
{
  return interface_sets_[ switch_->Id() ];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<CONTACT::CoInterface> >&
CONTACT::AUG::ComboStrategy::Interfaces() const
{
  return interface_sets_[ switch_->Id() ];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::Get()
{
  return dynamic_cast<CONTACT::AUG::Strategy&>( *strategies_[ switch_->Id() ] );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::Get() const
{
  return dynamic_cast<const CONTACT::AUG::Strategy&>( *strategies_[ switch_->Id() ] );
}
