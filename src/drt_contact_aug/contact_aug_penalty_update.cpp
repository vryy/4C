/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_penalty_update.cpp

\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

\maintainer Michael Hiermeier
\date Jul 28, 2017

*/
/*----------------------------------------------------------------------------*/


#include "contact_aug_penalty_update.H"
#include "contact_aug_potential.H"
#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_lagrange_multiplier_function.H"
#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_multiply.H"
#include "../drt_lib/epetra_utils.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::PenaltyUpdate* CONTACT::AUG::PenaltyUpdate::Create(
    const Teuchos::ParameterList& sa_params )
{
  const INPAR::CONTACT::PenaltyUpdate update_type =
      Teuchos::getIntegralValue<INPAR::CONTACT::PenaltyUpdate>( sa_params,
      "PENALTY_UPDATE" );

  return Create( update_type );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::PenaltyUpdate* CONTACT::AUG::PenaltyUpdate::Create(
    const INPAR::CONTACT::PenaltyUpdate update_type,
    const PenaltyUpdate* pu_src )
{
  switch ( update_type )
  {
    case INPAR::CONTACT::PenaltyUpdate::sufficient_lin_reduction:
    {
      // call copy constructor
      if ( pu_src )
        return new PenaltyUpdate_SufficientLinReduction( *pu_src );
      return new PenaltyUpdate_SufficientLinReduction();
    }
    case INPAR::CONTACT::PenaltyUpdate::lm_gap_ratio:
    {
      // call copy constructor
      if ( pu_src )
        return new PenaltyUpdate_LagrMultiplierGapRatio( *pu_src );
      return new PenaltyUpdate_LagrMultiplierGapRatio();
    }
    case INPAR::CONTACT::PenaltyUpdate::complementarity:
    {
      // call copy constructor
      if ( pu_src )
        return new PenaltyUpdate_Complementarity( *pu_src );
      return new PenaltyUpdate_Complementarity();
    }
    case INPAR::CONTACT::PenaltyUpdate::sufficient_angle:
    {
      // call copy constructor
      if ( pu_src )
        return new PenaltyUpdate_SufficientAngle( *pu_src );
      return new PenaltyUpdate_SufficientAngle();
    }
    case INPAR::CONTACT::PenaltyUpdate::none:
    {
      // call copy constructor
      if ( pu_src )
        return new PenaltyUpdate_Empty( *pu_src );
      return new PenaltyUpdate_Empty();
    }
    case INPAR::CONTACT::PenaltyUpdate::vague:
    {
      dserror( "You specified no PENALTY_UPDATE routine. Fix you input file!\n"
          "enum = %d | \"%s\"", update_type,
          INPAR::CONTACT::PenaltyUpdate2String( update_type ).c_str() );
      exit( EXIT_FAILURE );
    }
    default:
    {
      dserror( "Unknown penalty update type! (enum = %d | \"%s\")", update_type,
          INPAR::CONTACT::PenaltyUpdate2String( update_type ).c_str() );
      exit( EXIT_FAILURE );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Init(
    CONTACT::AUG::Strategy* const strategy,
    CONTACT::AUG::DataContainer* const data )
{
  strategy_ptr_ = strategy;
  data_ptr_ = data;
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::ThrowIfNotInitialized() const
{
  if ( not isinit_ )
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::SetState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir )
{
  ThrowIfNotInitialized();

  state_.Set( xold, dir, Data() );

  double dir_nrm2 = 0.0;
  dir.Norm2( &dir_nrm2 );
  dir_norm2_ = dir_nrm2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PostUpdate()
{
  PrintUpdate( IO::cout.os(IO::standard) );
  state_.Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PrintInfo( std::ostream& os ) const
{
  os << "\nCONTACT::AUG::PenaltyUpdate\n";
  os << "Type = " << INPAR::CONTACT::PenaltyUpdate2String( Type() ) << "\n";
  os << "isinit_ = " << ( isinit_ ? "TRUE" : "FALSE" ) << "\n";
  os << "strategy_ptr_ = " << strategy_ptr_ << "\n";
  os << "data_ptr_ = " << data_ptr_ << "\n";
  os << "dir_norm2_ = " << dir_norm2_ << "\n";
  os << "ratio_ = " << ratio_ << "\n";
  state_.Print( os );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::State::Set(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const CONTACT::AUG::DataContainer& data )
{
  xold_ = Teuchos::rcp( new Epetra_Vector( xold ) );
  full_direction_ = Teuchos::rcp( new Epetra_Vector( dir ) );
  wgap_ = Teuchos::rcp( new Epetra_Vector( *data.WGapPtr() ) );
  tributary_area_active_ = Teuchos::rcp( new Epetra_Vector( *data.KappaVecPtr() ) );
  tributary_area_inactive_ = Teuchos::rcp( new Epetra_Vector( *data.AVecPtr() ) );

  CONTACT::AUG::Potential& pot = pu_.Data().Potential();
  pot.Compute();
  gn_gn_ = pot.Get( POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active );
  pot.ComputeLin(dir);
  gn_dgn_ = pot.GetLin( POTENTIAL::Type::infeasibility_measure,
      POTENTIAL::SetType::active, POTENTIAL::LinTerm::wrt_d );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::State::Reset()
{
  xold_ = Teuchos::null;
  full_direction_ = Teuchos::null;
  wgap_ = Teuchos::null;
  tributary_area_active_ = Teuchos::null;
  tributary_area_inactive_ = Teuchos::null;

  gn_gn_ = 0.0;
  gn_dgn_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetDirection() const
{
  if ( full_direction_.is_null() )
    dserror( "Set the state first!" );

  return *full_direction_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  CONTACT::AUG::PenaltyUpdate::State::Print( std::ostream& os ) const
{
  os << "--- CONTACT::AUG::PenaltyUpdate::State object\n";
  os << "    <RCP> fulldirection_:" << full_direction_ << ")\n";
  os << "    <RCP> xold_: " << xold_ << "\n";
  os << "    <RCP> wgap_: " << wgap_ << "\n";
  os << "    <RCP> tributary_area_active_: " << tributary_area_active_ << "\n";
  os << "    <RCP> tributary_area_inactive_: " << tributary_area_inactive_ << "\n";
  os << "    <double> gn_gn_: " << gn_gn_ << "\n";
  os << "    <double> gn_dgn_: " << gn_dgn_ << "\n" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector&
CONTACT::AUG::PenaltyUpdate::State::GetPreviouslyAcceptedState() const
{
  if ( xold_.is_null() )
    dserror( "Set the state first!" );

  return *xold_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetWGap() const
{
  if ( wgap_.is_null() )
    dserror( "Set the state first!" );

  return *wgap_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector&
CONTACT::AUG::PenaltyUpdate::State::GetActiveTributaryArea() const
{
  if ( tributary_area_active_.is_null() )
    dserror( "Set the state first!" );

  return *tributary_area_active_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector&
CONTACT::AUG::PenaltyUpdate::State::GetInactiveTributaryArea() const
{
  if ( tributary_area_inactive_.is_null() )
    dserror( "Set the state first!" );

  return *tributary_area_inactive_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Update( const CONTACT::ParamsInterface& cparams )
{
  PreUpdate();

  Execute( cparams );

  PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::PenaltyUpdate::Get_DGapN( const Epetra_Vector& dincr_slma ) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp( new Epetra_Vector( Data().DLmNWGapLinMatrixPtr()->RangeMap(), true ) );

  CATCH_EPETRA_ERROR( Data().DLmNWGapLinMatrixPtr()->Multiply( false,
      dincr_slma, *dgapn_ptr ) );

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::PenaltyUpdate::Get_inconsistent_DGapN(
    const Epetra_Vector& dincr_slma ) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp( new Epetra_Vector( Data().BMatrixPtr()->RangeMap(), true ) );

  CATCH_EPETRA_ERROR( Data().BMatrixPtr()->Multiply( false, dincr_slma,
      *dgapn_ptr ) );

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
CONTACT::AUG::PenaltyUpdate::GetProblemRhs(
    const CONTACT::ParamsInterface& cparams,
    const std::vector<INPAR::STR::ModelType>* without_these_models ) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>( model );

  return cmodel.AssembleForceOfModels(without_these_models,true).getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix>
    CONTACT::AUG::PenaltyUpdate::GetStructuralStiffnessMatrix(
        const CONTACT::ParamsInterface& cparams ) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>( model );

  // access the full stiffness matrix
  Teuchos::RCP<const LINALG::SparseMatrix> full_stiff_ptr =
      cmodel.GetJacobianBlock( DRT::UTILS::block_displ_displ );
  return full_stiff_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PrintUpdate( std::ostream& os ) const
{
  double cnmax = 0.0;
  Data().Cn().MaxValue( &cnmax );

  os << "\n=== Update of the regularization parameter ===\n";
  os << "Type   = " << INPAR::CONTACT::PenaltyUpdate2String( Type() ) << "\n";
  os << "Ratio (cN_new / cn_old) = " << std::setw(10) << std::setprecision(4)
     << std::scientific << Ratio() << "\n";
  os << "New cN = " << std::setw(10) << std::setprecision(4)
     << std::scientific << cnmax << "\n";
  os << "==============================================\n" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate
CONTACT::AUG::PenaltyUpdate_Empty::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::none;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate
CONTACT::AUG::PenaltyUpdate_LagrMultiplierGapRatio::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::lm_gap_ratio;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_LagrMultiplierGapRatio::Execute(
    const CONTACT::ParamsInterface& cparams )
{
  ThrowIfNotInitialized();

  Data().Potential().Compute();

  const double infeasibility_measure =
      Data().Potential().Get(
          AUG::POTENTIAL::Type::infeasibility_measure,
          AUG::POTENTIAL::SetType::all );

  const double old_infeasibility_measure = Data().SaData().GetOldInfeasibilityMeasure();

  const double constr_violation     = std::sqrt( infeasibility_measure );
  const double old_constr_violation = std::sqrt( old_infeasibility_measure );

  /* do nothing if one of the following criteria is fulfilled
   *
   * ( 1 ) current constraint violation is already very low
   * ( 2 ) new constraint violation is four times smaller than the previous one
   * ( 3 ) constraint violation became larger over the last step,
   *       in this case we do not want to make things even worse by increasing
   *       cn.
   *
   * hiermeier 04/17 */
  if ( constr_violation < 1.0e-5 or
       constr_violation < 0.25 * old_constr_violation or
       constr_violation > old_constr_violation )
    return;

  // get the largest current cn value
  Epetra_Vector& cn_vec = Data().Cn();
  double cn_max = 0.0;
  cn_vec.MaxValue( &cn_max );

  const double cn_upper_bound = Data().SaData().GetCnUpperBound();
  if ( cn_max >= cn_upper_bound )
    return;

  const Epetra_Vector& zn_active = Data().Potential().GetZnActive();
  const Epetra_Vector& aWGap     = Data().AWGap();
  Teuchos::RCP<const Epetra_Vector> zold_ptr = Data().OldLmPtr();

  Epetra_Vector diff_zn_active( zn_active.Map() );
  LINALG::Export( *zold_ptr, diff_zn_active );

  CATCH_EPETRA_ERROR(diff_zn_active.Update( -1.0, zn_active, 1.0 ));

  // get norms
  double diff_zn_active_nrm = 0.0;
  diff_zn_active.Norm2( &diff_zn_active_nrm );

  double aWGap_nrm = 0.0;
  aWGap.Norm2( &aWGap_nrm );

  const double cn_new =
      std::min( std::max( 2.0 * diff_zn_active_nrm/aWGap_nrm, cn_max),
          cn_upper_bound );

  cn_vec.PutScalar( cn_new );

  Print2Screen( constr_violation, old_constr_violation, cn_max, cn_new );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_LagrMultiplierGapRatio::Print2Screen(
    const double constr_violation,
    const double old_constr_violation,
    const double cn_max,
    const double cn_new ) const
{
  IO::cout << "====================================================================\n";
  IO::cout << "                      INCREASE OF CN                                \n";
  IO::cout << "       |constr_violation^(k)| > 0.25 * |constr_violation^(k-1)|     \n"
           << "                  " << std::setprecision(5) << constr_violation << " > "
           << std::setprecision(5) << 0.25 * old_constr_violation << "\n";
  IO::cout << "    ->->-> cN^(k): " << cn_max << ", cN^(k+1): " << cn_new << " <-<-<-" << IO::endl;
  IO::cout << "====================================================================" << IO::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate
CONTACT::AUG::PenaltyUpdate_Complementarity::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::complementarity;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_Complementarity::Execute(
    const CONTACT::ParamsInterface& cparams )
{
  if ( Ratio() <= 1.0 )
    return;

  Data().Cn().Scale( Ratio() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_Complementarity::SetState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir )
{
  PenaltyUpdate::SetState( cparams, xold, dir );

  // split the previously accepted state into its different parts
  Teuchos::RCP<Epetra_Vector> displ_slma, zn_active, zn_inactive;
  Strategy().SplitStateVector( GetState().GetPreviouslyAcceptedState(),
      displ_slma, zn_active, zn_inactive );

  // split the direction into its different parts
  Teuchos::RCP<Epetra_Vector> dincr_slma, znincr_active, znincr_inactive;
  Strategy().SplitStateVector( GetState().GetDirection(), dincr_slma, znincr_active,
      znincr_inactive );

  Teuchos::RCP<const Epetra_Vector> dgapn_ptr = Get_DGapN( *dincr_slma );

  Epetra_Vector tmp( dincr_slma->Map(), true );
  Data().DGLmLinMatrix().Multiply( false, *dincr_slma, tmp );
  double d_ddglm_d = 0.0;
  tmp.Dot( *dincr_slma, &d_ddglm_d );

  double gapn_zn = 0.0;
  zn_active->Dot( GetState().GetWGap(), &gapn_zn );
//
//  double dgapn_dzn = 0.0;
//  dgapn_ptr->Dot( *znincr_active, &dgapn_dzn );

  double dgapn_zn = 0.0;
  dgapn_ptr->Dot( *zn_active, &dgapn_zn );
//
//  double gapn_dzn = 0.0;
//  znincr_active->Dot( GetState().GetWGap(), &gapn_dzn );

  // --------------------------------------------------------------------------
  Epetra_Vector awgapn( GetState().GetWGap() );
  MultiplyElementwise( GetState().GetActiveTributaryArea(),
      Data().GActiveNodeRowMap(), awgapn, true );
  double awgapn_dgapn = 0.0;
  awgapn.Dot( *dgapn_ptr, &awgapn_dgapn );

  Epetra_Vector sc_dgapn( *dgapn_ptr );
  sc_dgapn.Scale(0.5);

  MultiplyElementwise( GetState().GetActiveTributaryArea(),
      Data().GActiveNodeRowMap(), sc_dgapn, true );

  double sc_dgapn_dgapn = 0.0;
  sc_dgapn.Dot( *dgapn_ptr, &sc_dgapn_dgapn );

  if ( std::abs( awgapn_dgapn + sc_dgapn_dgapn ) > 1.0e-14 )
  {
//    IO::cout << "gapn_zn             = " << gapn_zn << "\n";
//    IO::cout << "d_ddglm_d           = " << d_ddglm_d << "\n";
//    IO::cout << "dgapn_zn            = " << dgapn_zn << "\n";
//    IO::cout << "sc_gapn_dgapn_dgapn = " << sc_gapn_dgapn_dgapn << "\n";

    double cn_new = ( gapn_zn + dgapn_zn ) / // + 0.5* d_ddglm_d
        ( std::min( awgapn_dgapn, - 2.0 * sc_dgapn_dgapn ) + sc_dgapn_dgapn );

    Ratio() = cn_new / Data().Cn()[0];
  }
  else
    Ratio() = 1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_Complementarity::ScaleDirection(
    Epetra_Vector& dir )
{
  return Ratio();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate
CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::sufficient_lin_reduction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::SetState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir )
{
  PenaltyUpdate::SetState( cparams, xold, dir );
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::Execute(
    const CONTACT::ParamsInterface& cparams )
{
  const State& state = GetState();

  if ( std::abs( state.gn_gn_ ) < 1.0e-30 )
  {
    Ratio() = 1.0;
    return;
  }

  const double gamma_theta = GammaTheta();
  Ratio() = 1.0 / (1.0 - gamma_theta) * ( 1.0 + state.gn_dgn_/state.gn_gn_ );

  if ( Ratio() > 1.0 )
    Data().Cn().Scale( Ratio() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::GammaTheta() const
{
  return Data().SaData().GetPenaltyCorrectionParameter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate
CONTACT::AUG::PenaltyUpdate_SufficientAngle::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::sufficient_angle;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientAngle::SetState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir )
{
  PenaltyUpdate::SetState( cparams, xold, dir );

  if ( std::abs( GetState().gn_gn_ ) < 1.0e-30 )
  {
    Ratio() = 1.0;
    return;
  }

  // split the previously accepted state into its different parts
  Teuchos::RCP<Epetra_Vector> displ_slma, zn_active, zn_inactive;
  Strategy().SplitStateVector( GetState().GetPreviouslyAcceptedState(),
      displ_slma, zn_active, zn_inactive );

  // split the direction into its different parts
  Teuchos::RCP<Epetra_Vector> dincr_slma, znincr_active, znincr_inactive;
  Strategy().SplitStateVector( GetState().GetDirection(), dincr_slma, znincr_active,
      znincr_inactive );

  Epetra_Vector tmp( dincr_slma->Map(), true );
  Data().DGLmLinMatrix().Multiply( false, *dincr_slma, tmp );
  double d_ddglm_d = 0.0;
  tmp.Dot( *dincr_slma, &d_ddglm_d );

  // directional derivative of the structural gradient
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(
          cparams.GetModelEvaluator() );

  const std::vector<INPAR::STR::ModelType> without_contact( 1, cmodel.Type() );
  Teuchos::RCP<Epetra_Vector> str_gradient =
      cmodel.AssembleForceOfModels( &without_contact, true );

  Epetra_Vector dincr_str( str_gradient->Map() );
  LINALG::ExtractMyVector( dir, dincr_str );

  double dstr_grad = 0.0;
  CATCH_EPETRA_ERROR( str_gradient->Dot( dincr_str, &dstr_grad ) );

  // directional derivative of the active constraint gradient
  Teuchos::RCP<const Epetra_Vector> dgapn_ptr =
      Get_inconsistent_DGapN( *dincr_slma );
  Epetra_Vector dgapn_active( *Data().GActiveNDofRowMapPtr() );
  LINALG::ExtractMyVector( *dgapn_ptr, dgapn_active );

  double dgapn_zn = 0.0;
  CATCH_EPETRA_ERROR( dgapn_active.Dot( *zn_active, &dgapn_zn ) );

  Epetra_Vector sc_dgapn( dgapn_active );
  MultiplyElementwise( GetState().GetActiveTributaryArea(),
      Data().GActiveNodeRowMap(), sc_dgapn, true );

  double sc_dgapn_dgapn = 0.0;
  CATCH_EPETRA_ERROR( sc_dgapn.Dot( dgapn_active, &sc_dgapn_dgapn ) );

  Epetra_Vector awgapn( *GetState().wgap_ );
  MultiplyElementwise( GetState().GetActiveTributaryArea(),
      Data().GActiveNodeRowMap(), awgapn, true );

  double awgapn_nrm2 = 0.0;
  awgapn.Norm2( &awgapn_nrm2 );

  double dgapn_nrm2 = 0.0;
  dgapn_active.Norm2( &dgapn_nrm2 );

  double cn_old = 0.0;
  Data().CnPtr()->MaxValue( &cn_old );

  const double gamma_phi = GammaAngle();

  const double cn_new = ( cn_old*gamma_phi*dgapn_nrm2*awgapn_nrm2 +
      d_ddglm_d - dstr_grad + dgapn_zn ) / sc_dgapn_dgapn;

  Ratio() = cn_new / cn_old;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientAngle::Execute(
    const CONTACT::ParamsInterface& cparams )
{
  if ( Ratio() > 1.0 )
    Data().Cn().Scale( Ratio() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_SufficientAngle::GammaAngle() const
{
  return Data().SaData().GetPenaltyCorrectionParameter();
}
