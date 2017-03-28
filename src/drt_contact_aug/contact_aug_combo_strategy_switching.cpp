/*---------------------------------------------------------------------*/
/*!
\file contact_aug_combo_strategy_switching.cpp

\brief This file contains the implementation of the switching strategy
       for the AUG::ComboStrategy.

\level 3

\maintainer Michael Hiermeier

\date Mar 24, 2017

*/
/*---------------------------------------------------------------------*/


#include "contact_aug_combo_strategy.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_contact/contact_strategy_factory.H"

#include "../drt_io/io_pstream.H"

#include "../linalg/linalg_utils.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AUG::ComboStrategy::Switching>
CONTACT::AUG::ComboStrategy::Switching::Create( ComboStrategy& combo )
{
  const Teuchos::ParameterList& p_combo =
      combo.Params().sublist( "AUGMENTED" ).sublist( "COMBO" );

  const enum INPAR::CONTACT::SwitchingStrategy switch_type =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SwitchingStrategy>(
          p_combo, "SWITCHING_STRATEGY" );

  switch ( switch_type )
  {
    case INPAR::CONTACT::switch_preasymptotic:
      return Teuchos::rcp( new Switching( combo, p_combo ) );
    default:
      dserror( "Unknown switching strategy! (switch_type = %d)",
          switch_type );
      exit( EXIT_FAILURE );
  }

  dserror( "Impossible to reach this point!" );
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::Switching::Switching(
    ComboStrategy& combo,
    const Teuchos::ParameterList& p_combo )
    : combo_( combo ),
      preasymptotic_id_( 0 ),
      asymptotic_id_( 0 ),
      is_asymptotic_( false ),
      strat_types_( 0 )
{
  if ( combo_.strategies_.size() > 2 )
    dserror( "This basic switching strategy supports maximal a number of "
        "two strategies. Feel free to add a new switching strategy, if you "
        "need more." );

  GetStrategyTypes( combo_.strategies_, strat_types_ );

  const enum INPAR::CONTACT::SolvingStrategy preasymptotic =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(
          p_combo, "STRATEGY_0" );
  preasymptotic_id_ = FindId( preasymptotic );

  const enum INPAR::CONTACT::SolvingStrategy asymptotic =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(
          p_combo, "STRATEGY_1" );
  asymptotic_id_ = FindId( asymptotic );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Switching::GetStrategyTypes(
    const plain_strategy_set& strategies,
    plain_strattype_set& strat_types ) const
{
  for ( plain_strategy_set::const_iterator cit = strategies.begin();
        cit != strategies.end(); ++cit )
  {
    const CONTACT::CoAbstractStrategy& s = (**cit);

    switch ( s.Type() )
    {
      case INPAR::CONTACT::solution_augmented:
      case INPAR::CONTACT::solution_steepest_ascent:
      case INPAR::CONTACT::solution_std_lagrange:
        strat_types.push_back( s.Type() );
        break;
      default:
        dserror( "The strategy is of a non-supported type! ( type = "
            "%s | %d )", INPAR::CONTACT::SolvingStrategy2String( s.Type() ).c_str(),
            s.Type() );
        exit( EXIT_FAILURE );
    }
  }

  if ( strategies.size() != strat_types.size() )
    dserror("Size mismatch! Something went wrong." );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::Switching::FindId(
    INPAR::CONTACT::SolvingStrategy sol_type ) const
{
  unsigned id = 0;
  for ( plain_strattype_set::const_iterator cit = strat_types_.begin();
        cit != strat_types_.end(); ++cit )
  {
    if ( *cit == sol_type )
      return id;
    ++id;
  }

  dserror( "Couldn't find the given SolvingStrategy! (sol_type = %s | %d)",
      INPAR::CONTACT::SolvingStrategy2String( sol_type ).c_str(), sol_type );
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::Switching::Id() const
{
  if ( is_asymptotic_ )
    return asymptotic_id_;

  return preasymptotic_id_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::Switching::Id(
    enum INPAR::CONTACT::SolvingStrategy sol_type ) const
{
  return FindId( sol_type );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Switching::Update(
    CONTACT::ParamsInterface& cparams )
{
  bool is_asymptotic = ( CheckPenetration() and CheckResidual( cparams ) );

  // if the status is the same as before, do nothing
  if ( is_asymptotic == is_asymptotic_ )
    return;

  // --------------------------------------------------------------------------
  // switch back to the pre-asymptotic phase
  // --------------------------------------------------------------------------
  if ( not is_asymptotic and
       not combo_.Get().ActiveSetSemiSmoothConverged() )
  {
    is_asymptotic_ = false;
  }
  // --------------------------------------------------------------------------
  // switch to the asymptotic phase
  // --------------------------------------------------------------------------
  else if ( is_asymptotic )
  {
    is_asymptotic_ = true;
  }

  STRATEGY::Factory::PrintStrategyBanner( strat_types_[ Id() ] );

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "is_asymptotic_ = " << ( is_asymptotic_ ? "TRUE" : "FALSE" ) << "\n";
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::Switching::CheckResidual(
    CONTACT::ParamsInterface& cparams )
{
  Teuchos::RCP<Epetra_Vector> force_no_dbc_ptr =
      GetStructuralForceWithoutDbcDofs( cparams );

  Teuchos::RCP<Epetra_Vector> str_slmaforce = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> constr_slmaforce = Teuchos::null;
  GetActiveSlMaForces( *force_no_dbc_ptr, str_slmaforce, constr_slmaforce );

  const bool angle_check = CheckAngleBetweenStrForceAndContactForce( *str_slmaforce,
      *constr_slmaforce );
  IO::cout << "angle_check = " << ( angle_check ? "PASSED" : "FAILED" ) << "\n";

  const bool res_check = CheckContactResidualNorm( *str_slmaforce,
      *constr_slmaforce );
  IO::cout << "res_check   = " << ( res_check ? "PASSED" : "FAILED" ) << "\n";

  return ( res_check and angle_check );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::Switching::CheckContactResidualNorm(
    const Epetra_Vector& str_slmaforce,
    const Epetra_Vector& constr_slmaforce ) const
{
  double str_nrm = 0.0;
  double res_nrm = 0.0;

#ifdef DEBUG_COMBO_STRATEGY
  force.Norm2( &res_nrm );
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "res_nrm = " << res_nrm << "\n";
#endif

  Epetra_Vector res_slma_active( str_slmaforce.Map(), false );
  res_slma_active.Update( 1.0, str_slmaforce, -1.0, constr_slmaforce, 0.0 );

  res_slma_active.Norm2( &res_nrm );

  str_slmaforce.Norm2( &str_nrm );
//  IO::cout << "str_nrm = " << str_nrm << "\n";
//
//  IO::cout << "ABSOLUTE: res_nrm = " << res_nrm << "\n";
//  if ( str_nrm > 0.0 )
//    IO::cout << "RELATIVE: res_nrm/str_nrm = " << res_nrm/str_nrm << "\n";

  return ( res_nrm <= 1.0e-3 * str_nrm );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::Switching::CheckAngleBetweenStrForceAndContactForce(
    const Epetra_Vector& str_slmaforce,
    const Epetra_Vector& constr_slmaforce ) const
{
  double constr_slmaforce_nrm = 0.0;
  constr_slmaforce.Norm2( &constr_slmaforce_nrm );

  double str_slmaforce_nrm = 0.0;
  str_slmaforce.Norm2( &str_slmaforce_nrm );

  const double nrm_prod = constr_slmaforce_nrm * str_slmaforce_nrm;

  double inner_prod = 0.0;
  constr_slmaforce.Dot( str_slmaforce, &inner_prod );

  const double angle_tol = 1.0e-6;

  return ( nrm_prod - inner_prod <= nrm_prod * angle_tol );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Switching::GetActiveSlMaForces(
    const Epetra_Vector& str_force,
    Teuchos::RCP<Epetra_Vector>& str_slmaforce,
    Teuchos::RCP<Epetra_Vector>& constr_slmaforce ) const
{
  Epetra_Vector& slforce = *combo_.no_dbc_.slForce_;
  Epetra_Vector& maforce = *combo_.no_dbc_.maForce_;

  slforce.PutScalar( 0.0 );
  maforce.PutScalar( 0.0 );

  LINALG::ExtractMyVector( *combo_.data_.SlForceLmPtr(), slforce );
  LINALG::ExtractMyVector( *combo_.data_.MaForceLmPtr(), maforce );

  Teuchos::RCP<Epetra_Map> gSlActiveForceMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> gMaActiveForceMap = Teuchos::null;
  GetGlobalSlMaActiveForceMaps( slforce, maforce, gSlActiveForceMap, gMaActiveForceMap );
  Teuchos::RCP<Epetra_Map> gSlMaActiveForceMap =
      LINALG::MergeMap( *gSlActiveForceMap, *gMaActiveForceMap );

  Epetra_Vector& slmaforce = *combo_.no_dbc_.slMaForce_;
  slmaforce.PutScalar( 0.0 );

  LINALG::AssembleMyVector( 1.0, slmaforce, 1.0, slforce );
  LINALG::AssembleMyVector( 1.0, slmaforce, 1.0, maforce );

  constr_slmaforce = LINALG::CreateVector( *gSlMaActiveForceMap, true );
  LINALG::ExtractMyVector( slmaforce, *constr_slmaforce );

  str_slmaforce = LINALG::CreateVector( *gSlMaActiveForceMap, true );
  LINALG::ExtractMyVector( str_force, *str_slmaforce );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void
CONTACT::AUG::ComboStrategy::Switching::GetGlobalSlMaActiveForceMaps(
    const Epetra_Vector& slforce,
    const Epetra_Vector& maforce,
    Teuchos::RCP<Epetra_Map>& gSlActiveForceMap,
    Teuchos::RCP<Epetra_Map>& gMaActiveForceMap ) const
{
  // initialize the map pointers
  gSlActiveForceMap =
      Teuchos::rcp( new Epetra_Map( 0, 0, *combo_.data_.CommPtr() ) );
  gMaActiveForceMap =
      Teuchos::rcp( new Epetra_Map( 0, 0, *combo_.data_.CommPtr() ) );

  Teuchos::RCP<Epetra_Map> imap = Teuchos::null;

  const plain_interface_set::const_iterator cit_begin =
      combo_.Get().Interfaces().begin();
  const plain_interface_set::const_iterator cit_end =
      combo_.Get().Interfaces().end();
  for ( plain_interface_set::const_iterator cit = cit_begin; cit != cit_end; ++cit )
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<AUG::Interface&>( **cit );

    imap = Teuchos::null;
    imap = interface.BuildActiveForceMap( slforce );
    if ( not imap.is_null() )
      gSlActiveForceMap = LINALG::MergeMap( *gSlActiveForceMap, *imap );

    imap = Teuchos::null;
    imap = interface.BuildActiveForceMap( maforce );
    if ( not imap.is_null() )
      gMaActiveForceMap = LINALG::MergeMap( *gMaActiveForceMap, *imap );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::AUG::ComboStrategy::Switching::GetStructuralForceWithoutDbcDofs(
    const CONTACT::ParamsInterface& cparams )
{
  Teuchos::RCP<Epetra_Vector> force_no_dbc_ptr =
      Teuchos::rcp( new Epetra_Vector( *combo_.no_dbc_.slMaMap_ ) );

  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(
          cparams.GetModelEvaluator() );

  combo_.is_assemblecontactrhs_ = false;
  Teuchos::RCP<Epetra_Vector> force_ptr = cmodel.AssembleForceOfAllModels();
  combo_.is_assemblecontactrhs_ = true;

  LINALG::Export( *force_ptr, *force_no_dbc_ptr );

  return force_no_dbc_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::Switching::CheckPenetration() const
{
  // get the overall largest penetration value
  double min_awgap = 0.0;
  combo_.data_.AWGap().MinValue( &min_awgap );

  const double penbound = GetPenetrationBound();

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "min_awgap = " << min_awgap << std::endl;
#endif

  const bool pen_check = ( min_awgap > penbound );

  IO::cout << "pen_check   = " << ( pen_check ? "PASSED" : "FAILED" ) << "\n";

  return pen_check;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::Switching::GetPenetrationBound() const
{
  double penbound = 1.0e12;
  for ( plain_interface_set::const_iterator cit = combo_.Interfaces().begin();
        cit != combo_.Interfaces().end(); ++cit )
  {
    const CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );
    penbound = std::min( penbound, interface.PenetrationBound() );
  }

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "penbound = " << penbound << std::endl;
#endif

  return ( penbound * -1.0 );
}
