/*---------------------------------------------------------------------*/
/*!
\file contact_aug_potential.cpp

\brief Class for the evaluation of the contact potential and its
       linearization.

\level 2

\maintainer Michael Hiermeier

\date Mar 14, 2017

*/
/*---------------------------------------------------------------------*/

#include "contact_aug_potential.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Potential::Potential( const CONTACT::AUG::Strategy& strategy,
    const CONTACT::AUG::DataContainer& data )
    : isvalidPotential_( false ),
      isvalidLinearization_( false ),
      isvalidState_( false ),
      isvalidDirection_( false ),
      issetup_( false ),
      strategy_( strategy ),
      data_( data ),
      zn_active_( Teuchos::null ),
      zn_inactive_( Teuchos::null ),
      zt_active_( Teuchos::null ),
      zt_inactive_( Teuchos::null ),
      dincrSlMa_( Teuchos::null ),
      znincr_active_( Teuchos::null ),
      znincr_inactive_( Teuchos::null ),
      potdata_(),
      lindata_()
{

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::ResetIsValid()
{
  isvalidPotential_     = false;
  isvalidLinearization_ = false;
  isvalidState_         = false;
  isvalidDirection_     = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::Setup( bool active_inactive_only )
{
  if ( not active_inactive_only )
  {
    dincrSlMa_ = Teuchos::rcp( new Epetra_Vector( *data_.GSlMaDofRowMapPtr() ) );
    isvalidDirection_ = false;
  }

  zn_active_ = Teuchos::rcp( new Epetra_Vector( *data_.GActiveNDofRowMapPtr() ) );
  znincr_active_ = Teuchos::rcp( new Epetra_Vector( *zn_active_ ) );

  Teuchos::RCP<Epetra_Map> ginactivendofs = LINALG::SplitMap(
      *data_.GSlNormalDofRowMapPtr(), *data_.GActiveNDofRowMapPtr() );
  zn_inactive_ = Teuchos::rcp( new Epetra_Vector( *ginactivendofs ) );
  znincr_inactive_ = Teuchos::rcp( new Epetra_Vector( *zn_inactive_ ) );

  isvalidState_ = false;
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::SetActiveInactiveState()
{
  LINALG::ExtractMyVector( *data_.LmPtr(), *zn_active_ );
  LINALG::ExtractMyVector( *data_.LmPtr(), *zn_inactive_ );

  LINALG::ExtractMyVector( *data_.LmIncrPtr(), *znincr_active_ );
  LINALG::ExtractMyVector( *data_.LmIncrPtr(), *znincr_inactive_ );

  isvalidState_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::SetDirection(
    const Epetra_Vector& direction )
{
  ResetIsValid();

  // extract slave/master part of the displacement increment
  Epetra_Vector dincr_exp( *data_.GDispDofRowMapPtr() );
  LINALG::Export( direction, dincr_exp );
  LINALG::ExtractMyVector( dincr_exp, *dincrSlMa_ );

  isvalidDirection_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::Compute()
{
  if ( isvalidPotential_ )
    return;

  if ( not issetup_ )
    dserror( "Call Setup() first!" );

  const Epetra_Vector& cn = data_.Cn();

  double lterms[4] = { 0.0, 0.0, 0.0, 0.0 };

  const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& co_interfaces =
      strategy_.ContactInterfaces();

  for ( std::vector<Teuchos::RCP<CONTACT::CoInterface> >::const_iterator cit =
      co_interfaces.begin(); cit != co_interfaces.end(); ++cit )
  {
    const CONTACT::AUG::Interface& interface =
        static_cast<const CONTACT::AUG::Interface&>( **cit );

    interface.AssembleContactPotentialTerms( cn, lterms[0], lterms[1], lterms[2], lterms[3] );
  }

  double gterms[4] = { 0.0, 0.0, 0.0, 0.0 };

  strategy_.Comm().SumAll( &lterms[0], &gterms[0], 4 );

  // copy results into the container
  potdata_.zn_gn_ = gterms[0];
  potdata_.gn_gn_ = gterms[1];
  potdata_.zn_zn_ = gterms[2];
  potdata_.zt_zt_ = gterms[3];

  // infeasibility values
  data_.WGapPtr()->Dot( *data_.WGapPtr(), &potdata_.inf_gn_gn_ );
  data_.InactiveRhsPtr()->Dot( *data_.InactiveRhsPtr(), &potdata_.inf_zn_zn_ );

  isvalidPotential_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Potential::ComputeLin()
{
  if ( isvalidLinearization_ )
    return;

  if ( not isvalidPotential_ )
    Compute();

  if ( not isvalidState_ )
    dserror( "Call SetState() first!" );

  if ( not isvalidDirection_ )
    dserror( "Call SetDirection() first!" );

  Epetra_Vector gradWGapD( *data_.GActiveNDofRowMapPtr() );

  int err = data_.DLmNWGapLinMatrixPtr()->Multiply( false, *dincrSlMa_, gradWGapD );
  if ( err )
    dserror( "Vector-matrix multiplication failed! (err=%d)", err );

  // --------------------------------------------------------------------------
  // Potential: Active contributions
  // --------------------------------------------------------------------------
  {
    // zn_k^T * gradWG(x_k)^T * dincr
    zn_active_->Dot( gradWGapD, &lindata_.zn_dgn_ );

    // wgn(x_k)^T * zincr
    data_.WGapPtr()->Dot( *znincr_active_, &lindata_.gn_dzn_ );

    // znincr^T * gradWG(x_k)^T * dincr
    znincr_active_->Dot( gradWGapD, &lindata_.dzn_dgn_ );

    // cn * awgn(x_k)^T * gradWG(x_k)^T * dincr
    Epetra_Vector scAWGap = Epetra_Vector( *data_.AWGapPtr() );
    const Epetra_Vector& cn = *data_.CnPtr();
    MultiplyElementwise( cn, *data_.GActiveNodeRowMapPtr(),
        scAWGap, false );
    scAWGap.Dot( gradWGapD, &lindata_.gn_dgn_ );
  }
  // --------------------------------------------------------------------------
  // Potential: Inactive contributions
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Map> ginactiveslnodes =
      LINALG::SplitMap( *data_.GSlNodeRowMapPtr(), *data_.GActiveNodeRowMapPtr() );
  Epetra_Vector scZnincr_inactive( *znincr_inactive_ );

  {
    // 1.0/cn * zn_k * A(x_k) * znincr
    MultiplyElementwise( *(data_.AVecPtr()), *ginactiveslnodes,
        scZnincr_inactive, false );
    MultiplyElementwise( *(data_.CnPtr()), *ginactiveslnodes,
        scZnincr_inactive, true );

    zn_inactive_->Dot( scZnincr_inactive, &lindata_.zn_dzn_ );

    znincr_inactive_->Dot( scZnincr_inactive, &lindata_.dzn_dzn_ );
  }
  // --------------------------------------------------------------------------
  // Infeasibility measure
  // --------------------------------------------------------------------------
  {
    // wGap^T * grad(wGap)^T * dincr [Active]
    data_.WGapPtr()->Dot( gradWGapD, &lindata_.inf_gn_dgn_ );

    // cn_k^(-2) * zn_k^T * A(x_k) * A(x_k) * znincr
    MultiplyElementwise( *(data_.AVecPtr()), *ginactiveslnodes,
        scZnincr_inactive, false );
    MultiplyElementwise( *(data_.CnPtr()), *ginactiveslnodes,
        scZnincr_inactive, true );

    zn_inactive_->Dot( scZnincr_inactive, &lindata_.inf_zn_dzn_ );
  }

  isvalidLinearization_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::Get( enum PotentialType pot_type,
    enum PotentialTerm pot_term ) const
{
  if ( not isvalidPotential_ )
    dserror( "Call Compute() first!" );

  switch( pot_type )
  {
    case type_lagrangian:
    {
      return GetLagrangian( pot_term );
    }
    case type_augmented_lagrangian:
    {
      return GetAugmentedLagrangian( pot_term );
    }
    case type_infeasibility_measure:
    {
      return GetInfeasibilityMeasure( pot_term );
    }
    default:
      dserror( "Unknown PotentialType enumerator ( enum = %d )", pot_type );
      exit( EXIT_FAILURE );
  }

  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetLin(
    enum PotentialType pot_type,
    enum LinearizationTerm lin_term ) const
{
  if ( not isvalidLinearization_ )
    dserror( "Call ComputeLin() first!" );

  switch( pot_type )
  {
    case type_lagrangian:
    {
      return GetLagrangianLin( lin_term );
    }
    case type_augmented_lagrangian:
    {
      return GetAugmentedLagrangianLin( lin_term );
    }
    case type_infeasibility_measure:
    {
      return GetInfeasibilityMeasureLin( lin_term );
    }
    default:
      dserror( "Unknown PotentialType enumerator ( enum = %d )", pot_type );
      exit( EXIT_FAILURE );
  }

  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetLagrangian( enum PotentialTerm pot_term ) const
{
  double val = 0.0;

  switch( pot_term )
  {
    case term_active:
    {
      val = - potdata_.zn_gn_;

      break;
    }
    case term_augmented:
    {
      val = 0.0;

      break;
    }
    case term_inactive:
    {
      val = - potdata_.zn_zn_;

      break;
    }
    case term_all:
    {
      val = -potdata_.zn_gn_ - potdata_.zn_zn_;

      break;
    }
    default:
      dserror( "Unknown PotentialTerm enumerator ( enum = %d )", pot_term );
      exit( EXIT_FAILURE );
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetLagrangianLin( enum LinearizationTerm lin_term ) const
{
  double val = 0.0;

  switch( lin_term )
  {
    case lin_active_wrt_d:
    {
      val = - lindata_.zn_dgn_;

      break;
    }
    case lin_active_wrt_z:
    {
      val = - lindata_.gn_dzn_;

      break;
    }
    case lin_augmented_wrt_d:
    {
      val = 0.0;

      break;
    }
    case lin_inactive_wrt_z:
    {
      val = -2.0 * lindata_.zn_dzn_;

      break;
    }
    case lin_active_wrt_d_and_z:
    {
      val = - lindata_.dzn_dgn_;

      break;
    }
    case lin_inactive_wrt_z_and_z:
    {
      val = - 2.0 * lindata_.dzn_dzn_;

      break;
    }
    default:
      dserror( "Unknown LinearizationTerm enumerator ( enum = %d )", lin_term );
      exit( EXIT_FAILURE );
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetAugmentedLagrangian(
    enum PotentialTerm pot_term ) const
{
  double val = 0.0;

  switch( pot_term )
  {
    case term_active:
    {
      val = - potdata_.zn_gn_;

      break;
    }
    case term_augmented:
    {
      val = 0.5 * potdata_.gn_gn_;

      break;
    }
    case term_inactive:
    {
      val = - 0.5 * potdata_.zn_zn_;

      break;
    }
    case term_all:
    {
      val = - potdata_.zn_gn_ + 0.5 * potdata_.gn_gn_ - 0.5 * potdata_.zn_zn_;

      break;
    }
    default:
      dserror( "Unknown PotentialTerm enumerator ( enum = %d )", pot_term );
      exit( EXIT_FAILURE );
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetAugmentedLagrangianLin(
    enum LinearizationTerm lin_term ) const
{
  double val = 0.0;

  switch( lin_term )
  {
    case lin_active_wrt_d:
    {
      val = - lindata_.zn_dgn_;

      break;
    }
    case lin_active_wrt_z:
    {
      val = - lindata_.gn_dzn_;

      break;
    }
    case lin_augmented_wrt_d:
    {
      val = lindata_.gn_dgn_;

      break;
    }
    case lin_inactive_wrt_z:
    {
      val = - lindata_.zn_dzn_;

      break;
    }
    case lin_active_wrt_d_and_z:
    {
      val = - lindata_.dzn_dgn_;

      break;
    }
    case lin_inactive_wrt_z_and_z:
    {
      val = - lindata_.dzn_dzn_;

      break;
    }
    default:
      dserror( "Unknown LinearizationTerm enumerator ( enum = %d )", lin_term );
      exit( EXIT_FAILURE );
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetInfeasibilityMeasure(
    enum PotentialTerm pot_term ) const
{
  double val = 0.0;

  switch( pot_term )
  {
    case term_active:
    {
      val = potdata_.inf_gn_gn_;

      break;
    }
    case term_augmented:
    {
      val = 0.0;

      break;
    }
    case term_inactive:
    {
      val = potdata_.inf_zn_zn_;

      break;
    }
    case term_all:
    {
      val = potdata_.inf_gn_gn_ + potdata_.inf_zn_zn_;

      break;
    }
    default:
      dserror( "Unknown PotentialTerm enumerator ( enum = %d )", pot_term );
      exit( EXIT_FAILURE );
  }

  return val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Potential::GetInfeasibilityMeasureLin(
    enum LinearizationTerm lin_term ) const
{
  double val = 0.0;

  switch( lin_term )
  {
    case lin_active_wrt_d:
    {
      val = - lindata_.inf_gn_dgn_;

      break;
    }
    case lin_active_wrt_z:
    {
      val = 0.0;

      break;
    }
    case lin_augmented_wrt_d:
    {
      val = 0.0;

      break;
    }
    case lin_inactive_wrt_z:
    {
      val = lindata_.inf_zn_dzn_;

      break;
    }
    case lin_active_wrt_d_and_z:
    {
      val = 0.0;

      break;
    }
    default:
      dserror( "Unknown LinearizationTerm enumerator ( enum = %d )", lin_term );
      exit( EXIT_FAILURE );
  }

  return val;
}
