/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_beaminteraction.cpp

\brief Evaluation of all beam interaction terms

\maintainer Jonas Eichinger

\level 3
*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_beaminteraction.H"
#include "beaminteraction_submodel_evaluator_factory.H"

#include "str_timint_base.H"
#include "str_utils.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "beaminteraction_submodel_evaluator_generic.H"
#include "beaminteraction_submodel_evaluator_crosslinking.H"
#include "beaminteraction_submodel_evaluator_beamcontact.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_crosslinking.H"

#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_particle/particle_handler.H"
#include "../drt_beam3/beam3_base.H"

#include "../drt_biopolynet/biopolynet_calc_utils.H"
#include "../drt_biopolynet/crosslinker_node.H"
#include "../drt_biopolynet/periodic_boundingbox.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FEVector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteraction::BeamInteraction():
    submodeltypes_ ( Teuchos::null ),
    me_map_ptr_ ( Teuchos::null ),
    me_vec_ptr_ ( Teuchos::null ),
    myrank_ ( -1 ),
    coupsia_ ( Teuchos::null ),
    siatransform_ ( Teuchos::null ),
    ia_discret_ ( Teuchos::null ),
    ia_state_ptr_ ( Teuchos::null ),
    ia_force_beaminteraction_ ( Teuchos::null ),
    force_beaminteraction_ ( Teuchos::null ),
    stiff_beaminteraction_ ( Teuchos::null ),
    particlehandler_ ( Teuchos::null ),
    bindis_ ( Teuchos::null ),
    rowbins_ ( Teuchos::null ),
    bin_beamcontent_ ( INPAR::BINSTRATEGY::Beam )
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Setup()
{
  CheckInit();

  // -------------------------------------------------------------------------
  // setup variables
  // -------------------------------------------------------------------------
  // stiff
  stiff_beaminteraction_ = Teuchos::rcp(new
      LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));
  // force
  force_beaminteraction_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));
  // get myrank
  myrank_ = DiscretPtr()->Comm().MyPID();

  // print logo
  Logo();

  // get submodel types
  SetSubModelTypes();

  // -------------------------------------------------------------------------
  // clone problem discretization, the idea is simple: we redistribute only
  // the new discretization to enable all interactions (including the required
  // search), calculate the resulting force and stiffness contributions, export
  // them no our initial discretization where all evaluation, assembly and
  // solving is done. Therefore the maps of our initial discretization don't
  // change, i.e. there is no need to rebuild the global state.
  // -------------------------------------------------------------------------
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase>  discloner =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  ia_discret_ = discloner->CreateMatchingDiscretization(DiscretPtr(),"intacdis");

  // init data container
  ia_state_ptr_ =  Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteractionDataState());
  ia_state_ptr_->Init();
  ia_state_ptr_->Setup(ia_discret_);

  ia_state_ptr_->GetMutableDisNp() = Teuchos::rcp( new Epetra_Vector(
      *GStatePtr()->GetMutableDisNp() ) );

  // -------------------------------------------------------------------------
  // initialize coupling adapter to transform matrices between the two discrets
  // (with distinct parallel distribution)
  // -------------------------------------------------------------------------
  coupsia_ = Teuchos::rcp(new ADAPTER::Coupling());
  siatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);

  // -------------------------------------------------------------------------
  // initialize and setup binning strategy and particle handler
  // -------------------------------------------------------------------------
  CreateBinDiscretization();

  // construct, init and setup particle handler and binning strategy
  particlehandler_ = Teuchos::rcp( new PARTICLE::ParticleHandler(myrank_) );
  particlehandler_->BinStrategy()->Init( bindis_, ia_discret_, ia_state_ptr_->GetMutableDisNp() );
  particlehandler_->BinStrategy()->Setup();

  // distribute problem according to bin distribution to procs
  PartitionProblem();

  // initialize submodel evaluators
  InitAndSetupSubModeEvaluators();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::SetSubModelTypes()
{
  CheckInit();

  submodeltypes_ = Teuchos::rcp(new std::set<enum INPAR::BEAMINTERACTION::SubModelType>());

  // ---------------------------------------------------------------------------
  // check for crosslinking in biopolymer networks
  // ---------------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CrosslinkingParams(), "CROSSLINKER"))
    submodeltypes_->insert(INPAR::BEAMINTERACTION::submodel_crosslinking);

  // ---------------------------------------------------------------------------
  // check for beam contact
  // ---------------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(
      DRT::Problem::Instance()->BeamContactParams(),"BEAMS_STRATEGY") != INPAR::BEAMCONTACT::bstr_none)
    submodeltypes_->insert(INPAR::BEAMINTERACTION::submodel_beamcontact);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::InitAndSetupSubModeEvaluators()
{
  CheckInit();

  // model map
  me_map_ptr_ =
      BEAMINTERACTION::SUBMODELEVALUATOR::BuildModelEvaluators(*submodeltypes_);

  // model vector
  me_vec_ptr_ = Teuchos::rcp(new STR::MODELEVALUATOR::BeamInteraction::Vector(0));
  STR::MODELEVALUATOR::BeamInteraction::Map::iterator miter;
  for ( miter = (*me_map_ptr_).begin(); miter != (*me_map_ptr_).end(); ++miter )
    me_vec_ptr_->push_back(miter->second);

  Vector::iterator me_iter;
  for ( me_iter = (*me_vec_ptr_).begin(); me_iter != (*me_vec_ptr_).end(); ++me_iter )
  {
    (*me_iter)->Init( ia_discret_, bindis_, GStatePtr(), ia_state_ptr_, particlehandler_,
        TimInt().GetDataSDynPtr()->GetPeriodicBoundingBox() );
    (*me_iter)->Setup();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::HaveModelType(
    INPAR::BEAMINTERACTION::SubModelType const& submodeltype) const
{
  CheckInit();
  return ( submodeltypes_->find(submodeltype) != submodeltypes_->end() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::PartitionProblem()
{
  CheckInit();

  // store structure discretization in vector
  std::vector< Teuchos::RCP<DRT::Discretization> > discret_vec(1);
  discret_vec[0] = ia_discret_;

  // displacement vector according to periodic boundary conditions
  std::vector< Teuchos::RCP<Epetra_Vector> > disnp(1);
  disnp[0] = Teuchos::rcp(new Epetra_Vector( *ia_discret_->DofColMap() ) );
  LINALG::Export( *ia_state_ptr_->GetMutableDisNp(), *disnp[0] );

  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector< std::map<int, std::vector<int> > > nodesinbin(1);

  // weight (this is experimental, choose what gives you best results)
  double const weight = 30.0;
  // get optimal row distribution of bins to procs
  rowbins_ = particlehandler_->BinStrategy()->WeightedDistributionOfBinsToProcs(
      discret_vec, disnp, nodesinbin, weight );

  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*bindis_->NodeRowMap()));
  // delete old bins
  bindis_->DeleteElements();
  particlehandler_->BinStrategy()->FillBinsIntoBinDiscretization( rowbins_ );

  // now node (=crosslinker) to bin (=element) relation needs to be
  // established in binning discretization. Therefore some nodes need to
  // change their owner according to the bins owner they are in
  particlehandler_->DistributeParticlesToBins(noderowmap);

  // redistribute, extend ghosting and assign element to bins
  RedistributeExtendGhostingAndAssign();

  // update maps of state vectors and matrices
  UpdateMaps();

  // reset transformation
  UpdateCouplingAdapterAndMatrixTransformation();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::RedistributeExtendGhostingAndAssign()
{
  CheckInit();

  std::map<int, std::set<int> > bintoelemap;
  RedistributeAndExtendGhostingOfProblemDiscret( bintoelemap );

  // build element to bin map according to extended ghosting
  // extbintoelemap_ than contains for each col element the bins that are
  // somehow touched by an axis aligned bounding box around the element
  bool cltoclinteraction = false; //fixme
  particlehandler_->BinStrategy()->ExtendBinGhosting( rowbins_,
      bintoelemap, cltoclinteraction );
  BuildEleToBinMap();

  // assign elements to bins
  particlehandler_->BinStrategy()->RemoveElesFromBins( bin_beamcontent_ );
  particlehandler_->BinStrategy()->AssignElesToBins( ia_discret_,
      ia_state_ptr_->GetMutableExtBinToEleMap(), bin_beamcontent_ );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::RedistributeAndExtendGhostingOfProblemDiscret(
    std::map<int, std::set<int> >& bintoelemap)
{
  CheckInit();

  // -------------------------------------------------------------------------
  // standard ghosting
  // -------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Map> dummy1;
  Teuchos::RCP<Epetra_Map> dummy2;
  std::map<int, std::vector<int> > nodesinbin;

  particlehandler_->BinStrategy()->StandardDiscretizationGhosting(ia_discret_,
      rowbins_, ia_state_ptr_->GetMutableDisNp(), dummy1, dummy2, nodesinbin);

  // ----------------------------------------------------------------------
  // extended ghosting
  // ----------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector( *ia_discret_->DofColMap() ) );
  LINALG::Export( *ia_state_ptr_->GetMutableDisNp(), *iadiscolnp );

  particlehandler_->BinStrategy()->ExtendDiscretizationGhosting(ia_discret_,
      rowbins_, iadiscolnp, bintoelemap, ia_state_ptr_->GetMutableExtBinToEleMap());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // get current displacement state and export to interaction discretization dofmap
  UpdateDofMapOfVector(ia_discret_, ia_state_ptr_->GetMutableDisNp(), GState().GetMutableDisNp());
  // update colume vector
  ia_state_ptr_->GetMutableDisColNp() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_state_ptr_->GetDisNp(), *ia_state_ptr_->GetMutableDisColNp());

  // submodel loop
  Vector::iterator me_iter;
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end() ; ++me_iter )
    (*me_iter)->Reset();

  // Zero out force and stiffness contributions
  force_beaminteraction_->PutScalar(0.0);
  ia_force_beaminteraction_->PutScalar(0.0);
  ia_state_ptr_->GetMutableForceNp()->PutScalar(0.0);
  stiff_beaminteraction_->Zero();
  ia_state_ptr_->GetMutableStiff()->Zero();

  // update gidmap_ and exporter in matrix transform object
  // Note: we need this in every evaluation call (i.e. every iteration) because a change
  //       in the active set of element pairs changes the entries of the used coarse
  //       system stiffness matrix (because we only assemble non-zero values).
  //       Therefore, the graph of the matrix changes and also the required gidmap
  //       (even in computation with one processor)
  UpdateCouplingAdapterAndMatrixTransformation();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateForce()
{
  CheckInitSetup();

  Vector::iterator me_iter;
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->EvaluateForce();

  // do communication
  if( ia_state_ptr_->GetMutableForceNp()->GlobalAssemble( Add, false ) != 0 )
    dserror("GlobalAssemble failed");
  // add to non fe vector
  if ( ia_force_beaminteraction_->Update( 1., *ia_state_ptr_->GetMutableForceNp(), 1. ) )
    dserror("update went wrong");

  // transformation from ia_discret to problem discret
  TransformForce();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateStiff()
{
  CheckInitSetup();

  ia_state_ptr_->GetMutableStiff()->UnComplete();

  Vector::iterator me_iter;
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->EvaluateStiff();

  if ( not ia_state_ptr_->GetMutableStiff()->Filled() )
    ia_state_ptr_->GetMutableStiff()->Complete();

  TransformStiff();

  if ( not stiff_beaminteraction_->Filled() )
    stiff_beaminteraction_->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::EvaluateForceStiff()
{
  CheckInitSetup();

  ia_state_ptr_->GetMutableStiff()->UnComplete();

  Vector::iterator me_iter;
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->EvaluateForceStiff();

  // do communication
  if( ia_state_ptr_->GetMutableForceNp()->GlobalAssemble( Add, false ) != 0 )
    dserror("GlobalAssemble failed");
  // add to non fe vector
  if ( ia_force_beaminteraction_->Update( 1., *ia_state_ptr_->GetMutableForceNp(), 1. ) )
    dserror("update went wrong");

  if (not ia_state_ptr_->GetMutableStiff()->Filled())
    ia_state_ptr_->GetMutableStiff()->Complete();

  TransformForceStiff();

  if (not stiff_beaminteraction_->Filled())
    stiff_beaminteraction_->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  CheckInitSetup();

  STR::AssembleVector( 1.0, f, timefac_np, *force_beaminteraction_ );

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::BeamInteraction::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  CheckInitSetup();

  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add( *stiff_beaminteraction_, false, timefac_np, 1.0 );

  // no need to keep it
  stiff_beaminteraction_->Zero();
  ia_state_ptr_->GetMutableStiff()->Zero();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
 // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetMutableFstructureOld();

  fstructold_ptr->Update( -timefac_n, *force_beaminteraction_, 1.0 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateStepElement()
{
  CheckInitSetup();

  // submodel loop
  Vector::iterator me_iter;
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->PreUpdateStepElement();

  // transfer beams and particles to new bins
  TransferBinContent();

  // update maps of state vectors and matrices
  UpdateMaps();

  // submodel loop
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->UpdateStepElement();

  // submodel loop
  for ( me_iter = me_vec_ptr_->begin(); me_iter != me_vec_ptr_->end(); ++me_iter )
    (*me_iter)->PostUpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::DetermineStressStrain()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::DetermineEnergy()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::BeamInteraction::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::BeamInteraction::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::BeamInteraction::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::PostOutput()
{
  // mesh is not written to disc, only maximum node id is important for output
  bindis_->Writer()->ParticleOutput(GState().GetStepN(), GState().GetTimeN(), false);
  bindis_->Writer()->NewStep(GState().GetStepN(), GState().GetTimeN());
  Teuchos::RCP<Epetra_Vector> dis = LINALG::CreateVector(*bindis_->DofRowMap(),true);
  Teuchos::RCP<Epetra_Vector> numbond = LINALG::CreateVector(*bindis_->NodeRowMap(),true);
  Teuchos::RCP<Epetra_Vector> owner = LINALG::CreateVector(*bindis_->NodeRowMap(),true);

  //todo: this is of course not nice, this needs to be done somewhere else
  for(int i=0;i<bindis_->NumMyRowNodes();++i)
  {
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(i));
    // std::vector holding gids of dofs
    std::vector<int> dofnode  = bindis_->Dof(crosslinker_i);

    // loop over all dofs
    for(int dim=0;dim<3;++dim)
    {
      int doflid = dis->Map().LID(dofnode[dim]);
      (*dis)[doflid] = crosslinker_i->X()[dim];
    }

    (*numbond)[i] = crosslinker_i->ClData()->GetNumberOfBonds();
    (*owner)[i] = crosslinker_i->Owner();
  }
  bindis_->Writer()->WriteVector("displacement", dis);
  bindis_->Writer()->WriteVector("numbond", numbond, bindis_->Writer()->nodevector);
  bindis_->Writer()->WriteVector("owner", owner, bindis_->Writer()->nodevector);


  // as we know that our maps have changed every time we write output, we can empty
  // the map cache as we can't get any advantage saving the maps anyway
  bindis_->Writer()->ClearMapCache();

//  particlehandler_->BinStrategy()->WriteBinOutput(GState().GetStepN(),GState().GetTimeN());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::ResetStepState()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::TransferBinContent()
{
  CheckInitSetup();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::TransferAndRedistribute");

  // transfer crosslinker to their new bins
  particlehandler_->TransferParticles();

  // redistribute, ghost and assign eles to bins
  RedistributeExtendGhostingAndAssign();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateCouplingAdapterAndMatrixTransformation()
{
  CheckInit();

  // reset transformation member variables (eg. exporter) by rebuilding
  // and provide new maps for coupling adapter
  siatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  coupsia_->SetupCoupling(*ia_discret_, Discret());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::CreateNewBins(bool newxaabb,
                                                         bool newcutoff)
{
  CheckInitSetup();

  // todo: unshifted current positions are needed here
  if(newcutoff)
    dserror("unshifted configuration is needed (not yet here) for calculation of cutoff.");

  // store structure discretization in vector
  std::vector<Teuchos::RCP<DRT::Discretization> > discret_vec(1);
  discret_vec[0] = ia_discret_;
  // displacement vector according to periodic boundary conditions
  std::vector<Teuchos::RCP<Epetra_Vector> > disnp(1);
  disnp[0] = ia_state_ptr_->GetMutableDisNp();

  // create XAABB and optionally set cutoff radius
  if(newxaabb)
    particlehandler_->BinStrategy()->CreateXAABB( discret_vec, disnp, newcutoff );
  // just set cutoff radius
  else if (newcutoff)
    particlehandler_->BinStrategy()->ComputeMinCutoff( discret_vec, disnp );

  particlehandler_->BinStrategy()->CreateBins();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::BuildEleToBinMap()
{
  CheckInit();

  // delete old map
  ia_state_ptr_->GetMutableExtEleToBinMap().clear();
  // loop over bins
  std::map<int, std::set<int> >::const_iterator biniter;
  for( biniter = ia_state_ptr_->GetExtBinToEleMap().begin() ; biniter != ia_state_ptr_->GetExtBinToEleMap().end(); ++biniter )
  {
    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for( eleiter = biniter->second.begin(); eleiter != biniter->second.end(); ++eleiter )
    {
      int elegid = *eleiter;
      int bingid = biniter->first;
      // assign bins to elements
      ia_state_ptr_->GetMutableExtEleToBinMap()[elegid].insert(bingid);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::CreateBinDiscretization()
{
  CheckInit();

  // clone communicator
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(DiscretPtr()->Comm().Clone());
  bindis_ = Teuchos::rcp(new DRT::Discretization("particle" ,com));
  // create discretization writer - in constructor set into and owned by
  // corresponding discret
  bindis_->SetWriter( Teuchos::rcp(new IO::DiscretizationWriter(bindis_)) );

  if(HaveModelType(INPAR::BEAMINTERACTION::submodel_crosslinking))
    AddParticlesToBinDiscret();

  // set row map of newly created particle discretization
  bindis_->FillComplete( false, false, false );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::AddParticlesToBinDiscret()
{
  CheckInit();

  const int numcrosslinker =
      DRT::Problem::Instance()->CrosslinkingParams().get<int> ("NUMCROSSLINK");

  // -------------------------------------------------------------------------
  // set range for uniform random number generator
  // -------------------------------------------------------------------------
  DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);
  std::vector<double> randpos;
  DRT::Problem::Instance()->Random()->Uni(randpos,3*numcrosslinker);

  // -------------------------------------------------------------------------
  // initialize crosslinker, i.e. add nodes (according to number of crosslinker
  // you want) to bin discretization and set their random reference position
  // -------------------------------------------------------------------------
  // only proc 0 is doing this (as the number of crosslinker is manageable)
  if(myrank_==0)
  {
    for (int i=0; i<numcrosslinker; i++)
    {
      // random reference position of crosslinker in bounding box
      std::vector<double> X(3);
      for (int dim=0; dim<3; dim++)
      {
        double edgelength = TimInt().GetDataSDynPtr()->GetPeriodicBoundingBox()->EdgeLength(dim);
        double min = TimInt().GetDataSDynPtr()->GetPeriodicBoundingBox()->min(dim);
        X[dim] = min + ( edgelength * randpos[i+dim]);
      }

      Teuchos::RCP<DRT::Node> newcrosslinker =
          Teuchos::rcp(new CROSSLINKING::CrosslinkerNode(i,&X[0],myrank_));

      // todo: put next two calls in constructor of CrosslinkerNode?
      // init crosslinker data container
      Teuchos::RCP<CROSSLINKING::CrosslinkerNode> clnode =
          Teuchos::rcp_dynamic_cast<CROSSLINKING::CrosslinkerNode>(newcrosslinker);
      clnode->InitializeDataContainer();
      // set material
      // fixme: assign matnum to crosslinker type in crosslinker section in input file
      //       for now, only one linker type with matnum 2
      clnode->SetMaterial(2);

      // add crosslinker to bin discretization
      bindis_->AddNode(newcrosslinker);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateMaps()
{
  CheckInit();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export. See in particle_timint
  // how this can be done.
  //todo: check if update is necessary (->SameAs())

  // beam displacement
  UpdateDofMapOfVector( ia_discret_, ia_state_ptr_->GetMutableDisNp() );

  // get current displacement state and export to interaction discretization dofmap
  UpdateDofMapOfVector( ia_discret_, ia_state_ptr_->GetMutableDisNp(), GState().GetMutableDisNp() );
  // update colume vector
  ia_state_ptr_->GetMutableDisColNp() = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_state_ptr_->GetDisNp(), *ia_state_ptr_->GetMutableDisColNp());

  // force
  ia_force_beaminteraction_ = Teuchos::rcp(new Epetra_Vector( *ia_discret_->DofRowMap() ,true) );
  ia_state_ptr_->GetMutableForceNp() = Teuchos::rcp(new Epetra_FEVector( *ia_discret_->DofRowMap() ,true) );

  // stiff
  ia_state_ptr_->GetMutableStiff() = Teuchos::rcp(new
      LINALG::SparseMatrix(*ia_discret_->DofRowMap(), 81, true, true,LINALG::SparseMatrix::FE_MATRIX));

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::UpdateDofMapOfVector(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Epetra_Vector>&      dofmapvec,
    Teuchos::RCP<Epetra_Vector>       old)
{
  CheckInit();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export. See in particle_timint
  // how this can be done.

  if (dofmapvec != Teuchos::null)
  {
    if(old==Teuchos::null)
      old = dofmapvec;
    dofmapvec = LINALG::CreateVector(*discret->DofRowMap(),true);
    LINALG::Export(*old, *dofmapvec);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::TransformForce()
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::BeamInteraction::TransformForce");

  // transform force vector to problem discret layout/distribution
  force_beaminteraction_ = coupsia_->MasterToSlave(ia_force_beaminteraction_);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::TransformStiff()
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::BeamInteraction::TransformStiff");

  stiff_beaminteraction_->UnComplete();
  // transform stiffness matrix to problem discret layout/distribution
  (*siatransform_)( *ia_state_ptr_->GetMutableStiff(), 1.0,
      ADAPTER::CouplingMasterConverter(*coupsia_), *stiff_beaminteraction_, false );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::TransformForceStiff()
{
  CheckInit();

  TransformForce();
  TransformStiff();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteraction::Logo() const
{
  CheckInit();

  if( myrank_ == 0 )
  {
    IO::cout << "\n****************************************************************" << IO::endl;
    IO::cout << "*                                                              *" << IO::endl;
    IO::cout << "*          Welcome to a Biopolymer Network Simulation          *" << IO::endl;
    IO::cout << "*                                                              *" << IO::endl;
    IO::cout << "****************************************************************" << IO::endl;
    IO::cout << "                                                                  " << IO::endl;
    IO::cout << "                                                                  " << IO::endl;
    IO::cout << "                      0=========================0                 " << IO::endl;
    IO::cout << "                    //|   \\            /       /||                " << IO::endl;
    IO::cout << "                   // |    \\ |       |/       //||                " << IO::endl;
    IO::cout << "                  //  |  /  \\|       /       // ||                " << IO::endl;
    IO::cout << "                 //   |  \\   \\   /  /|\\     //  ||                " << IO::endl;
    IO::cout << "                //    |  /   |\\ /  / | \\   //   ||                " << IO::endl;
    IO::cout << "               //     |  \\   | \\     |  \\ //  / ||                " << IO::endl;
    IO::cout << "              //  \\  /|  /   |/      |   //  /  ||                " << IO::endl;
    IO::cout << "              0=========================0 \\ /   ||                " << IO::endl;
    IO::cout << "             ||    /\\ |____          |  || \\    ||                " << IO::endl;
    IO::cout << "             ||   /  \\|    \\   ------   ||/ \\   ||                " << IO::endl;
    IO::cout << "             ||  /    |                 ||      ||                " << IO::endl;
    IO::cout << "             || /     0----------/------||------0-                " << IO::endl;
    IO::cout << "             ||      /   /       \\      ||     //                 " << IO::endl;
    IO::cout << "             ||     /___/  \\     /    / ||    //                  " << IO::endl;
    IO::cout << "             ||    /        \\    \\   /  ||   //                   " << IO::endl;
    IO::cout << "             ||   /  \\/\\/\\/  \\   /  /   ||  //                    " << IO::endl;
    IO::cout << "             ||  /      /     \\  \\ /    || //                     " << IO::endl;
    IO::cout << "             || /      /         /      ||//                      " << IO::endl;
    IO::cout << "             ||/                       /||/                       " << IO::endl;
    IO::cout << "              0=========================0                         " << IO::endl;
    IO::cout << "                                                                     " << IO::endl;
    IO::cout << "                                                                     " << IO::endl;
  }

}

///*-----------------------------------------------------------------------------*
// | nodes are checked and transferred if necessary              eichinger 09/16 |
// *-----------------------------------------------------------------------------*/
//void PARTICLE::Algorithm::TransferNodes(
//    Teuchos::RCP<DRT::Discretization> discret,
//    Teuchos::RCP<Epetra_Vector>       disnp)
//{
//  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferNodes");
//
//  // list of homeless nodes
//  std::list<Teuchos::RCP<DRT::Node> > homelessnodes;
//
//  // check in each bin whether nodes have moved out
//  // first run over nodes and then process whole bin in which node is located
//  // until all particles have been checked
//  std::map<int, std::set<int> > rownodetobin;
//  std::map<int, std::set<int> >::const_iterator rownodeiter;
//
//  // loop over all row nodes
//  for(rownodeiter=rownodetobin.begin(); rownodeiter!=rownodetobin.end(); ++rownodeiter)
//  {
//    // get current node
//    DRT::Node *currnode = discret->gNode(rownodeiter->first);
//    // as checked above, there is only one element in currele array
//    const int binId = rownodeiter->second;
//
//    // get current node position
//    double currpos[3] = {0.0,0.0,0.0};
//    GetCurrentNodePos(discret,currnode,disnp,currpos);
//
//    // get current bin of node
//    const int gidofbin = ConvertPosToGid(currpos);
//    if(gidofbin != binId) // node has left current bin
//    {
//      // gather all node Ids that will be removed and remove them afterwards
//      // (looping over nodes and deleting at the same time is detrimental)
//      tobemoved.push_back(currnode->Id());
//      // find new bin for particle
//      /*bool placed = */PlaceNodeCorrectly(Teuchos::rcp(currnode,false), currpos, homelessparticles);
//    }
//
//
//
//    // finally remove nodes from their old bin
//    for(size_t iter=0; iter<tobemoved.size(); iter++)
//    {
//      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currbin)->DeleteNode(tobemoved[iter]);
//    }
//
//  } // end for ibin
//
//#ifdef DEBUG
//  if(homelessparticles.size())
//    std::cout << "There are " << homelessparticles.size() << " homeless particles on proc" << myrank_ << std::endl;
//#endif
//
//  // homeless particles are sent to their new processors where they are inserted into their correct bin
//  FillParticlesIntoBinsRemoteIdList(homelessparticles);
//
//  // check whether all procs have a filled bindis_,
//  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
//  bindis_->CheckFilledGlobally();
//
//  // new ghosting if necessary
//  if (ghosting)
//    bindis_->ExtendedGhosting(*bincolmap_,true,false,true,false);
//  else
//    bindis_->FillComplete(true, false, true);
//
//  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
//  bool rebuildwallpointer = true;
//  if(moving_walls_)
//    rebuildwallpointer = false;
//  BuildElementToBinPointers(rebuildwallpointer);
//
//  // update state vectors in time integrator to the new layout
//  if(updatestates)
//  {
//    TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles::UpdateStates");
//    particles_->UpdateStatesAfterParticleTransfer();
//    UpdateStates();
//  }
//
//  return;
//}
