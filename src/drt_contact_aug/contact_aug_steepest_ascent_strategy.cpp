/*---------------------------------------------------------------------*/
/*!
\file contact_aug_steepest_ascent_strategy.cpp

\brief Steepest ascent solution strategy based on the augmented contact
       formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 7, 2017

*/
/*---------------------------------------------------------------------*/

#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_steepest_ascent_interface.H"
#include "contact_aug_potential.H"

#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_mortar/mortar_utils.H"

#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::DataContainer::DataContainer()
    : AUG::DataContainer(),
      old_infeasibility_measure_( 0.0 ),
      cn_upper_bound_( 0.0 )
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params,
    const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
    int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof )
    : CONTACT::AUG::Strategy( data_ptr, DofRowMap, NodeRowMap, params, interfaces,
        dim, comm, maxdof )
{
  saDataPtr_ = Teuchos::rcp_dynamic_cast<CONTACT::AUG::STEEPESTASCENT::DataContainer>(
      data_ptr, true );

  const Teuchos::ParameterList& sa_params =
      Params().sublist("AUGMENTED",true).sublist("STEEPESTASCENT",true);
  Data().SetCnUpperBound( sa_params.get<double>( "CN_UPPER_BOUND" ) );

  // cast to steepest ascent interfaces
  for ( std::vector<Teuchos::RCP<CONTACT::CoInterface> >::const_iterator cit =
        interfaces.begin(); cit != interfaces.end(); ++cit )
  {
    const Teuchos::RCP<CONTACT::CoInterface> & interface = *cit;
    interface_.push_back( Teuchos::rcp_dynamic_cast<CONTACT::AUG::STEEPESTASCENT::Interface>(
        interface, true ) );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::STEEPESTASCENT::Strategy::AssembleContactRHS()
{
  if ( not AUG::Strategy::AssembleContactRHS() )
    return false;

  for ( InterfaceVector::const_iterator cit=interface_.begin();
        cit!=interface_.end(); ++cit )
  {
    const STEEPESTASCENT::Interface& interface = **cit;
    interface.AssembleAugAVector( Data().AVec(), Data().KappaVec() );
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::AddContributionsToConstrRHS(
    Epetra_Vector& augConstrRhs ) const
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
CONTACT::AUG::STEEPESTASCENT::Strategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new LINALG::SparseMatrix(SlMaDoFRowMap(true),100,false,true));

      // build matrix kdd
      AddContributionsToMatrixBlockDisplDispl( *mat_ptr );
      mat_ptr->Complete(SlMaDoFRowMap(true),SlMaDoFRowMap(true));

      // transform parallel row/column distribution of matrix kdd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            SlMaDoFRowMapPtr(false),SlMaDoFRowMapPtr(false));

      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      // do nothing

      break;
    }
    case DRT::UTILS::block_lm_displ:
    {
      // do nothing

      break;
    }
    case DRT::UTILS::block_lm_lm:
    {
      // do nothing

      break;
    }
    default:
    {
      dserror("Unknown STR::MatBlockType!");
      break;
    }
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::RunPreComputeX(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    Epetra_Vector& dir_mutable )
{
  AugmentDirection( cparams, xold, dir_mutable );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::RunPostIterate(
    const CONTACT::ParamsInterface& cparams )
{
  UpdateCn( cparams );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::AugmentDirection(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    Epetra_Vector& dir_mutable )
{
  // if there are no contact contributions, do a direct return
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;

  // just do it for the default step case
  if ( not cparams.IsDefaultStep() )
    return;

  // extract displ. increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_ptr =
      LINALG::ExtractMyVector( dir_mutable, *ProblemDofs() );
  Teuchos::RCP<const Epetra_Vector> displ_incr_redistributed_ptr = Teuchos::null;
  if ( ParRedist() )
  {
    Teuchos::RCP<Epetra_Vector> tmp_exp = Teuchos::rcp( new Epetra_Vector(
        *gdisprowmap_ ) );
    LINALG::Export( *displ_incr_ptr, *tmp_exp );
    displ_incr_redistributed_ptr = tmp_exp;
  }
  else
    displ_incr_redistributed_ptr = displ_incr_ptr;

  const Epetra_Vector& displ_incr = *displ_incr_redistributed_ptr;

  // --------------------------------------------------------------------------
  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active =
      ComputeActiveLagrangeIncrInNormalDirection( displ_incr );

  // --------------------------------------------------------------------------
  // inactive lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> zold_ptr =
      LINALG::ExtractMyVector( xold, LMDoFRowMap( false ) );
  zold_ptr->ReplaceMap( SlDoFRowMap( false ) );

  Teuchos::RCP<Epetra_Vector> zold_redistributed_ptr = Teuchos::null;
  if ( ParRedist() )
  {
    zold_redistributed_ptr =
        Teuchos::rcp( new Epetra_Vector( SlDoFRowMap( true ) ) );
    // export the zold vector to the zold redistributed vector
    LINALG::Export( *zold_ptr, *zold_redistributed_ptr );
  }
  else
    zold_redistributed_ptr = zold_ptr;

  const Epetra_Vector& zold = *zold_redistributed_ptr;

  Teuchos::RCP<Epetra_Vector> zincr_inactive =
      ComputeInactiveLagrangeIncrInNormalDirection( displ_incr, zold );

  // --------------------------------------------------------------------------
  // assemble the Lagrange multiplier contributions
  Epetra_Vector zincr_redistributed( SlDoFRowMap( true ) );
  LINALG::AssembleMyVector( 0.0, zincr_redistributed, 1.0, *znincr_active );
  LINALG::AssembleMyVector( 0.0, zincr_redistributed, 1.0, *zincr_inactive );
  zincr_redistributed.ReplaceMap( LMDoFRowMap( true ) );

  Epetra_Vector zincr_full( LMDoFRowMap( false ) );
  LINALG::Export( zincr_redistributed, zincr_full );

  LINALG::AssembleMyVector( 0.0, dir_mutable, 1.0, zincr_full );

  // run at the very end...
  PostAugmentDirection();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::STEEPESTASCENT::Strategy::
ComputeActiveLagrangeIncrInNormalDirection( const Epetra_Vector& displ_incr )
{

  // active lagrange multiplier increment in normal direction
  Teuchos::RCP<Epetra_Vector> znincr_active_ptr =
      Teuchos::rcp( new Epetra_Vector( Data().GActiveNDofRowMap(), true ) );
  Epetra_Vector& znincr_active = *znincr_active_ptr;

  // nothing to do, if there are no active contributions
  if ( not IsInContact() )
    return znincr_active_ptr;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx12,tempmtx21,tempmtx22;

  Teuchos::RCP<LINALG::SparseMatrix> gradWGapUpdate;
  // Split dLmNWGapLinMatrix_
  LINALG::SplitMatrix2x2( Data().DLmNWGapLinMatrixPtr(),
      Data().GActiveNDofRowMapPtr(), emptymap, gdisprowmap_, emptymap,
      gradWGapUpdate, tempmtx12, tempmtx21, tempmtx22 );

  // calculate the Uzawa Update increment
  // *** Attention: zincr_Update has the wrong sign here! ***
  int err = gradWGapUpdate->Multiply( false, displ_incr, znincr_active );
  if ( err )
    dserror("Multiply error! (err=%d)", err);

  znincr_active.Update( 1.0, Data().WGap(), 1.0 );

  // Scaling of the Lagrange multiplier increment
  // --> inverse area scaling
  MultiplyElementwise( Data().KappaVec(), Data().GActiveNodeRowMap(),
      znincr_active, true );

  // Update the final Lagrange multiplier increment.
  // These values will also be used to update the nodal quantities during
  // the recover routine.
  MultiplyElementwise( Data().Cn(), Data().GActiveNodeRowMap(), znincr_active, false );

  /* Step length parameter: Variations are possible, see for example:
   * "Constrained optimization and Lagrange multiplier methods",
   * Dimitri P. Bertsekas, 1996, pp.125-133
   * (interpolation strategy on the pages 132 and 133). -- hiermeier 03/17 */
  const double beta_zn = 1.0;
  // We correct the increment sign.
  znincr_active.Scale( -beta_zn );

  return znincr_active_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::STEEPESTASCENT::Strategy::
ComputeInactiveLagrangeIncrInNormalDirection(
    const Epetra_Vector& displ_incr,
    const Epetra_Vector& zold )
{
  Teuchos::RCP<Epetra_Map> ginactivedofs = LINALG::SplitMap(SlDoFRowMap(true),
      Data().GActiveDofRowMap());

  // inactive lagrange multiplier increment in normal and tangential direction
  Teuchos::RCP<Epetra_Vector> zincr_inactive_ptr =
      Teuchos::rcp( new Epetra_Vector( *ginactivedofs ) );
  Epetra_Vector& zincr_inactive = *zincr_inactive_ptr;

  // extract old lagrange multipliers in normal and tangential direction
  Teuchos::RCP<const Epetra_Vector> zold_inactive_ptr =
      LINALG::ExtractMyVector( zold, *ginactivedofs );
  const Epetra_Vector& zold_inactive = *zold_inactive_ptr;

  // extract displ increment
  Teuchos::RCP<const Epetra_Vector> displ_incr_sl_ptr =
        LINALG::ExtractMyVector( displ_incr, SlDoFRowMap(true) );

  int err = Data().InactiveLinMatrix().Multiply( false, *displ_incr_sl_ptr,
      zincr_inactive );
  if ( err )
    dserror("Multiply error (err=%d)!", err);

  zincr_inactive.ReciprocalMultiply( -1.0, Data().InactiveDiagMatrix(), zincr_inactive, 0.0 );

  zincr_inactive.Update( -1.0, zold_inactive, 1.0 );

  return zincr_inactive_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::PostAugmentDirection()
{
  Data().Potential().Compute();

  const double infeasibility_measure =
      Data().Potential().Get(
          AUG::Potential::type_infeasibility_measure,
          AUG::Potential::term_all );

  Data().SetOldInfeasibilityMeasure( infeasibility_measure );

//  Data().Potential().ComputeLin();
//  const double gndzn = Data().Potential().GetLin(
//      AUG::Potential::type_augmented_lagrangian, AUG::Potential::lin_active_wrt_z );
//  const double zndzn = Data().Potential().GetLin(
//      AUG::Potential::type_augmented_lagrangian, AUG::Potential::lin_inactive_wrt_z );
//  const double dzndzn = Data().Potential().GetLin(
//      AUG::Potential::type_augmented_lagrangian, AUG::Potential::lin_inactive_wrt_z_and_z );
//
//  std::cout << "\n\n=============================================\n";
//  std::cout << "( " << gndzn << " + " << zndzn << " ) / " << dzndzn << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Strategy::UpdateCn(
    const CONTACT::ParamsInterface& cparams )
{
  Data().Potential().Compute();

  const double infeasibility_measure =
      Data().Potential().Get(
          AUG::Potential::type_infeasibility_measure,
          AUG::Potential::term_all );

  const double old_infeasibility_measure = Data().GetOldInfeasibilityMeasure();

  const double constr_violation     = std::sqrt( infeasibility_measure );
  const double old_constr_violation = std::sqrt( old_infeasibility_measure );

  // do nothing if one of the following criteria is fulfilled
  if ( constr_violation < 1.0e-5 or
       constr_violation < 0.25 * old_constr_violation )
    return;

  // get the largest current cn value
  Epetra_Vector& cn_vec = Data().Cn();
  double cn_max = 0.0;
  cn_vec.MaxValue( &cn_max );

  const double cn_upper_bound = Data().GetCnUpperBound();
  if ( cn_max >= cn_upper_bound )
    return;

  const Epetra_Vector& zn_active = Data().Potential().GetZnActive();
  const Epetra_Vector& aWGap     = Data().AWGap();

  // get norms
  double zn_active_nrm = 0.0;
  zn_active.Norm2( &zn_active_nrm );

  double aWGap_nrm = 0.0;
  aWGap.Norm2( &aWGap_nrm );

  const double cn_new = std::min( std::max( 2.0 * zn_active_nrm/aWGap_nrm, cn_max), cn_upper_bound );
  cn_vec.PutScalar( cn_new );

  IO::cout << "====================================================================\n";
  IO::cout << "                      INCREASE OF CN                                \n";
  IO::cout << "       |constr_violation^(k)| > 0.25 * |constr_violation^(k-1)|     \n"
           << "                  " << std::setprecision(5) << constr_violation << " > "
           << std::setprecision(5) << 0.25 * old_constr_violation << "\n";
  IO::cout << "    ->->-> cN^(k): " << cn_max << ", cN^(k+1): " << cn_new << " <-<-<-" << IO::endl;
  IO::cout << "====================================================================" << IO::endl;
}
