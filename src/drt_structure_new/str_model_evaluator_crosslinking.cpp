/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_crosslinking.cpp

\brief model evaluator for crosslinking in biopolymer networks

\maintainer Jonas Eichinger

\date May, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_crosslinking.H"

#include "str_model_evaluator_data.H"
#include "str_timint_databiopolynetdyn.H"
#include "str_timint_base.H"
#include "str_utils.H"
#include "str_integrator.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_particle/particle_algorithm.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_biopolynet/biopolynet_calc_utils.H"
#include "../drt_biopolynet/crosslinker_node.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Crosslinking::Crosslinking():
  eval_statmech_ptr_(Teuchos::null),
  myrank_(-1),
  numproc_(-1),
  intactdis_(Teuchos::null),
  particlealgo_(Teuchos::null),
  bindis_(Teuchos::null),
  rowbins_(Teuchos::null),
  bin_beamcontent_(INPAR::BINSTRATEGY::Beam)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Setup()
{
  CheckInit();
  // -------------------------------------------------------------------------
  // get pointer to statmech data
  // -------------------------------------------------------------------------
  eval_statmech_ptr_ = EvalData().StatMechPtr();
  // -------------------------------------------------------------------------
  // get myrank
  // -------------------------------------------------------------------------
  myrank_ = DiscretPtr()->Comm().MyPID();
  // -------------------------------------------------------------------------
  // get number of procs
  // -------------------------------------------------------------------------
  numproc_ = DiscretPtr()->Comm().NumProc();
  // -------------------------------------------------------------------------
  // adapt displacement vector so that node positions are consistent with
  // periodic boundary condition. Note: Input file contains the unshifted
  // configuration so that we do not have to set up the elements twice.
  // Shifting to periodic boundary configuration is done here. From now on the
  // global displacement vector always contains the shifted configuration.
  // -------------------------------------------------------------------------
  // todo: here on in brownian dyn (or both, but in which reset??)
  STATMECH::UTILS::PeriodicBoundaryConsistentDis(
      GStatePtr()->GetMutableDisN(),                            // disn
      eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength(),
      DiscretPtr());
  STATMECH::UTILS::PeriodicBoundaryConsistentDis(
      GStatePtr()->GetMutableDisNp(),                           // disnp
      eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength(),
      DiscretPtr());

  // -------------------------------------------------------------------------
  // clone problem discretization, the idea is simple: we redistribute only
  // the new discretiztion to enable all interactions (including the required
  // search), calculate the resulting force and stiffnes contributions, export
  // them no our initial discretization where all evaluation, assembly and
  // solving is done. Therefore the maps of our initial discretization don't
  // change, i.e. there is no need to rebuild the global state.
  // -------------------------------------------------------------------------
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase>  discloner =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  intactdis_ = discloner->CreateMatchingDiscretization(DiscretPtr(),"intacdis");

  // -------------------------------------------------------------------------
  // initialize crosslinker, i.e. add nodes (according to number of crosslinker
  // you want) to bin discretization and set their random reference position
  // -------------------------------------------------------------------------
  InitializeBinDiscret();
  // -------------------------------------------------------------------------
  // initialize particle algorithm, although we don't need/use any of the
  // actual particle "algorithm" or integrator, we are just using some nice
  // methods here.
  // -------------------------------------------------------------------------
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParams();
  /// algorithm is created here
  particlealgo_ =  Teuchos::rcp(new PARTICLE::Algorithm(DiscretPtr()->Comm(),params));

  // -------------------------------------------------------------------------
  // build periodic boundary conditions in binning strategy
  // -------------------------------------------------------------------------
  if(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(0)>0.0)
    particlealgo_->BuildPeriodicBC(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength());

  // -------------------------------------------------------------------------
  // a complete partitioning including the creation of bins, their weighted
  // distribution to procs as well a new distribution of the beam discret is
  // done here
  // -------------------------------------------------------------------------
  // get initial (shifted) displacement vector, from now on this vector is based
  // on the maps of the interaction discretization
  ia_disnp_ = Teuchos::rcp(new Epetra_Vector(*GStatePtr()->GetMutableDisNp()));
  UpdateBinStrategy(false,true);

  // gather data for all column crosslinker initially
  const int numcolcl = bindis_->NumMyColNodes();
  crosslinker_data_.resize(numcolcl);
  PreComputeCrosslinkerData(numcolcl);

  // gather data for all column beams
  const int numcolbeams = intactdis_->NumMyColElements();
  beam_data_.resize(numcolbeams);
  PreComputeBeamData(numcolbeams);

  // set flag
  issetup_ = true;

  // that's it
  return;

} // Setup()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  return;
} // Reset()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateForce()
{
  CheckInitSetup();
  bool ok = true;

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateStiff()
{
  CheckInitSetup();
  bool ok = true;

  // that's it
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateForceStiff()
{
  CheckInitSetup();
  bool ok = true;

  // that's it
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  CheckInitSetup();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  CheckInitSetup();
  // nothing to do
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{
  CheckInitSetup();

  return;
} // WriteRestart()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();

  return;
} // ReadRestart()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
 // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();


  return;
} // UpdateStepState()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateStepElement()
{
  // -------------------------------------------------------------------------
  // update all binding states and redistribute bindis and intactdis
  // -------------------------------------------------------------------------
  UpdateCrosslinking();

  return;

} // UpdateStepElement()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DetermineStressStrain()
{

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DetermineEnergy()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

//  // mesh is not written to disc, only maximum node id is important for output
//  bindis_->Writer()->ParticleOutput(GState().GetStepN(), GState().GetTimeN(), false);
//  bindis_->Writer()->NewStep(GState().GetStepN(), GState().GetTimeN());
//  Teuchos::RCP<Epetra_Vector> dis = LINALG::CreateVector(*bindis_->DofRowMap(),true);
//  Teuchos::RCP<Epetra_Vector> numbond = LINALG::CreateVector(*bindis_->NodeRowMap(),true);
//
//  //todo: this is of course not nice, this needs to be done somewhere else
//  for(int i=0;i<bindis_->NodeRowMap()->NumMyElements();++i)
//  {
//    CROSSLINKING::CrosslinkerNode *crosslinker_i =
//        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(i));
//    // std::vector holding gids of dofs
//    std::vector<int> dofnode  = bindis_->Dof(crosslinker_i);
//
//    // loop over all dofs
//    for(int dim=0;dim<3;++dim)
//    {
//      int doflid = dis->Map().LID(dofnode[dim]);
//      (*dis)[doflid] = crosslinker_i->X()[dim];
//    }
//
//    (*numbond)[i] = crosslinker_i->ClData()->GetNumberOfBonds();
//  }
//  bindis_->Writer()->WriteVector("displacement", dis);
//  bindis_->Writer()->WriteVector("numbond", numbond, bindis_->Writer()->nodevector);
//  // as we know that our maps have changed every time we write output, we can empty
//  // the map cache as we can't get any advantage saving the maps anyway
//  bindis_->Writer()->ClearMapCache();

  return;
} //OutputStepState()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Crosslinking::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Crosslinking::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Crosslinking::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PostOutput()
{
  CheckInitSetup();

  return;
} // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ResetStepState()
{
  CheckInitSetup();

  return;
} //ResetStepState()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateBinStrategy(bool transfer,
                                                          bool partition,
                                                          bool repartition,
                                                          bool createxaabb,
                                                          bool setcutoff)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::UpdateBinStrategy");

  if(repartition)
    dserror("repartioning not yet functional, will follow soon");

  // some safety checks
  if((transfer && partition) || (transfer && repartition) || (partition && repartition))
    dserror("You can either do a transfer, partitioning or repartitioning.");
  if((repartition || transfer) && (createxaabb || setcutoff))
    dserror("If you want to reset XAABB and/or cutoff you need to do a full partitioning");

  // store structure discretization in vector
  std::vector<Teuchos::RCP<DRT::Discretization> > discret_vec(1);
  discret_vec[0] = intactdis_;
  // displacement vector according to periodic boundary conditions
  std::vector<Teuchos::RCP<Epetra_Vector> > disnp(1);
  disnp[0] = ia_disnp_;
  // todo: unshifted current positions are needed here
  if(setcutoff)
    dserror("unshifted configuration is needed (not yet here) for calculation of cutoff.");

  // create XAABB and optionally set cutoff radius
  if(createxaabb)
    particlealgo_->CreateXAABB(discret_vec,disnp,setcutoff);
  // just set cutoff radius
  else if (setcutoff)
    particlealgo_->ComputeMaxCutoff(discret_vec,disnp);

  // -------------------------------------------------------------------------
  // Create bins according to XAABB and cutoff set in constructor (read
  // of input file)
  // -------------------------------------------------------------------------
  if(!repartition && !transfer)
    particlealgo_->CreateBins(Teuchos::null);

  // -------------------------------------------------------------------------
  // assign bins to procs according to a weighted partitioning, i.e. bins
  // are weighted with respect to the number of nodes they contain and then
  // distributed to the procs. Then an optimal distribution of bins to
  // procs can be obtained
  // -------------------------------------------------------------------------
  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int> > > nodesinbin(1);
  if(partition || repartition)
  {
    // get optimal row distribution of bins to procs
    rowbins_ = particlealgo_->WeightedDistributionOfBinsToProcs(
        discret_vec,disnp,nodesinbin,repartition);
  }

  // -------------------------------------------------------------------------
  // build element rowmap of binning discretization according to rowbins
  // -------------------------------------------------------------------------
  // completely new bins are created
  if(partition)
  {
    // extract noderowmap because it will be called Reset() after adding elements
    Teuchos::RCP<Epetra_Map> noderowmap =
        Teuchos::rcp(new Epetra_Map(*bindis_->NodeRowMap()));

    // delete old bins
    bindis_->DeleteElements();
    // loop over all bins on this proc, bindis than contains bins as elements
    for(int i=0; i<rowbins_->NumMyElements(); i++)
    {
      const int gid = rowbins_->GID(i);
      Teuchos::RCP<DRT::Element> bin =
          DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy",gid,myrank_);
      bindis_->AddElement(bin);
    }

    // -----------------------------------------------------------------------
    // now node (=crosslinker) to bin (=element) relation needs to be
    // established in binning discretization. Therefore some nodes need to
    // change their owner according to the bins owner they are in
    // -----------------------------------------------------------------------
    // i) create a set of homeless particles (owned by this proc) that do not
    //    reside in a bin owned by this proc
    std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less> homelesscrosslinker;
    for (int lid = 0; lid < noderowmap->NumMyElements(); ++lid)
    {
      DRT::Node* node = bindis_->gNode(noderowmap->GID(lid));
      const double* currpos = node->X();
      particlealgo_->PlaceNodeCorrectly(
          Teuchos::rcp(node,false), currpos, homelesscrosslinker);
    }
    // ii) round robin loop to assign homeless crosslinker to correct bin and owner
    particlealgo_->FillParticlesIntoBinsRoundRobin(homelesscrosslinker);

  }
  // bin gids are the same
  else if(repartition)
  {
    // -----------------------------------------------------------------------
    // export row elements/bins to new layout
    // -----------------------------------------------------------------------
    bindis_->ExportRowElements(*rowbins_);
    // -----------------------------------------------------------------------
    // export row nodes to new layout
    // -----------------------------------------------------------------------
    // create a set of row crosslinker IDs for each proc
    std::set<int> crosslinker;
    for (int lid=0; lid<rowbins_->NumMyElements(); ++lid)
    {
      DRT::Element* bin = bindis_->gElement(rowbins_->GID(lid));
      const int* crosslinkerids = bin->NodeIds();
      for(int icl=0; icl<bin->NumNode(); ++icl)
        crosslinker.insert(crosslinkerids[icl]);
    }
    // copy crosslinkergids to a vector and create crosslinkerrowmap
    std::vector<int> rowcrosslinker(crosslinker.begin(),crosslinker.end());
    Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(
        -1,(int)rowcrosslinker.size(),&rowcrosslinker[0],0,bindis_->Comm()));

    // place all nodes on the correct processor
    bindis_->ExportRowNodes(*noderowmap);
  }
  else if (transfer)
  {
    // transfer crosslinker to their new bins
    particlealgo_->TransferParticles();
  }
  else
  {
    dserror("You should not be here");
  }

  // -------------------------------------------------------------------------
  // call fill complete to build new node row map
  // -------------------------------------------------------------------------
  bindis_->FillComplete(false,false,false);

  // -------------------------------------------------------------------------
  // each owner of a bin gets owner of the nodes (of the structure discret) this
  // bin contains. All other nodes of elements, of which proc is owner of at
  // least one node, are ghosted.
  // -------------------------------------------------------------------------
  // as extended ghosting is applied to discret_vec, colmaps of standard ghosting
  // are given back separately if needed
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmap;
  particlealgo_->StandardGhosting(intactdis_,rowbins_,ia_disnp_,stdelecolmap,
      stdnodecolmap,nodesinbin[0]);

#ifdef DEBUG
  // print distribution after standard ghosting
  if(myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"   parallel distribution with standard ghosting   " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*intactdis_);
#endif

  // ----------------------------------------------------------------------
  // extended ghosting means the following here: Each proc ghosts
  // all elements whose XAABB cuts a bin that is next to a bin that is
  // owned by a proc an not empty. All associated nodes are ghosted as well
  // ----------------------------------------------------------------------
  // to distribute eles to bins correctly, we need a column map for disnp here
  // as the owner of a element does not need to be the owner of all nodes,
  // but we need the information of all nodes for correct distribution
  // (note this vector is col format before extenden ghosting, later we
  // a col vector based on maps of extended ghosting)
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector(*intactdis_->DofColMap()));
  LINALG::Export(*ia_disnp_, *iadiscolnp);

  // extend ghosting
  particlealgo_->ExtendBinGhosting(intactdis_,rowbins_,iadiscolnp,
      extbintoelemap_,true,true);

  // build element to bin map according to extended ghosting
  BuildEleToBinMap();

#ifdef DEBUG
  // print distribution after extended ghosting
  if(myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"   parallel distribution with extended ghosting   " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*intactdis_);
#endif

  // -------------------------------------------------------------------------
  // one layer ghosting of bins and particles, i.e. each proc ghosts the
  // bins (multibin elements, elemap) and particles (nodes, nodemap) of
  // its 26 neighbored bins -> final FillComplete() for maps included
  // -------------------------------------------------------------------------
  particlealgo_->BuildBinColMap(rowbins_,extbintoelemap_);

#ifdef DEBUG
  // print distribution after extended ghosting
  if (myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"           particles after ghosting     " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*bindis_);
#endif

  //--------------------------------------------------------------------------
  // update vectors and matrices as maps have changed during redistribution
  //--------------------------------------------------------------------------
  UpdateMaps();

  //--------------------------------------------------------------------------
  // assign beam elements to bins to enable crosslinker to beam interaction
  //--------------------------------------------------------------------------
  // in case have not build new bins we need to clear the content
  if(!partition)
    particlealgo_->RemoveElesFromBins(bin_beamcontent_);
  // note: extbintoelemap contains all necessary information for this step as
  // there is no assigning to do for empty bins (not in extbintoelemap) on a proc
  // bins that just have crosslinker are contained as well but do no damage
  // loop over all bins and remove assigned wall elements
  particlealgo_->AssignElesToBins(intactdis_, extbintoelemap_, bin_beamcontent_);

  // that's it
  return;

} // UpdateBinStrategy()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BuildEleToBinMap()
{
  CheckInit();

  // delete old map
  exteletobinmap_.clear();
  // loop over bins
  std::map<int, std::set<int> >::const_iterator biniter;
  for(biniter = extbintoelemap_.begin(); biniter!=extbintoelemap_.end(); ++biniter)
  {
    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for(eleiter = biniter->second.begin(); eleiter!=biniter->second.end(); ++eleiter)
    {
      int elegid = *eleiter;
      int bingid = biniter->first;
      // assign bins to elements
      exteletobinmap_[elegid].insert(bingid);
    }
  }

  // that is it
  return;
} //BuildEleToBinMap()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::InitializeBinDiscret()
{
  CheckInit();
  // -------------------------------------------------------------------------
  // add bin discretization manually (as it is not part of the input file)
  // to the global problem to enable output and other stuff (e.g. that particle
  // methods work here)
  // -------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(DiscretPtr()->Comm().Clone());
  bindis_ = Teuchos::rcp(new DRT::Discretization("particle" ,com));
  // create discretization writer - in constructor set into and owned by
  // corresponding discret
  bindis_->SetWriter(
      Teuchos::rcp(new IO::DiscretizationWriter(bindis_)));
  DRT::Problem::Instance()->AddDis("particle", bindis_);
  // -------------------------------------------------------------------------
  // set range for uniform random number generator
  // -------------------------------------------------------------------------
  DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);
  // -------------------------------------------------------------------------
  // add nodes to bin discretization
  // only proc 0 is doing this (as the number of crosslinker is manageable)
  // -------------------------------------------------------------------------
  if(!myrank_)
  {
    for (int i=0; i<eval_statmech_ptr_->GetDataSMDynPtr()->NumCrosslink(); i++)
    {
      // random reference position of crosslinker in bounding box
      std::vector<double> X(3);
      for (int dim=0; dim<3; dim++)
        X[dim] = eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(dim)
                 * DRT::Problem::Instance()->Random()->Uni();

      Teuchos::RCP<DRT::Node> newcrosslinker =
          Teuchos::rcp(new CROSSLINKING::CrosslinkerNode(i,&X[0],myrank_));

      // todo: put next two calls in constructor of CrosslinkerNode?
      // init crosslinker data container
      Teuchos::RCP<CROSSLINKING::CrosslinkerNode> clnode =
          Teuchos::rcp_dynamic_cast<CROSSLINKING::CrosslinkerNode>(newcrosslinker);
      clnode->InitializeDataContainer();
      // set material
      // todo: assign matnum to crosslinker type in crosslinker section in input file
      //       for now, only one linker type with matnum 2
      clnode->SetMaterial(2);

      // add crosslinker to bin discretization
      bindis_->AddNode(newcrosslinker);
    }
  }
  // -------------------------------------------------------------------------
  // set row map of newly created particle discretization
  // -------------------------------------------------------------------------
  bindis_->FillComplete(false,false,false);

  return;
} // InitializeBinDiscret()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DiffuseCrosslinker()
{
  CheckInit();

  // get standard deviation and mean value for crosslinker that are free to
  // diffuse
  double standarddev = sqrt(eval_statmech_ptr_->GetDataSMDynPtr()->KT() /
                       (2*M_PI * eval_statmech_ptr_->GetDataSMDynPtr()->Eta()
                       * eval_statmech_ptr_->GetDataSMDynPtr()->RLink())
                       * (*GState().GetDeltaTime())[0]);
  double meanvalue = 0.0;
  // Set mean value and standard deviation of normal distribution
  DRT::Problem::Instance()->Random()->SetMeanVariance(meanvalue,standarddev);

  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector(*intactdis_->DofColMap()));
  LINALG::Export(*ia_disnp_, *iadiscolnp);

  // loop over all row crosslinker (beam binding status not touched here)
  const int numrowcl = bindis_->NumMyRowNodes();
  for(int rowcli=0; rowcli<numrowcl; ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(rowcli));
    const int clgid = crosslinker_i->Id();
    const int clcollid = bindis_->NodeColMap()->LID(clgid);
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];
    LINALG::Matrix<3,1>& clpos_i = cldata_i.clpos;
    std::vector<std::pair<int, int> >& clbspots_i = cldata_i.clbspots;
    int& clnumbond_i = cldata_i.clnumbond;

    // different treatment according to number of bonds a crosslinker has
    switch(clnumbond_i)
    {
      // crosslinker has zero bonds, i.e. is free to diffuse according to
      // brownian dynamics
      case 0:
      {
        // -----------------------------------------------------------------
        // compute crosslinker diffusion for this time step
        // -----------------------------------------------------------------
        std::vector<double> randvec;
        int count = 3;
        DRT::Problem::Instance()->Random()->Normal(randvec,count);
        // note: check for compliance with periodic boundary conditions is
        // done during crosslinker transfer in UpdateBinStrategy()
        crosslinker_i->ChangePos(randvec);

        // that is it
        break;
      } // clnumbond_i = 0

      // crosslinker has one bond (cl gets current position of filament bspot it is
      // attached to)
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = 0;
        // check, which clbspot is free
        if(clbspots_i[0].first < 0)
          occbspotid  = 1;

#ifdef DEBUG
          // safety check
        if(clbspots_i[occbspotid].first < 0)
          dserror("clnumbond_i doesn't match clbspots_i");
#endif

        // get current position of binding spot of filament partner
        // note: we can not use our beam data container, as bspot position is not current position (as this
        // is the result of a sum, you can not have a reference to that)
        const int elegid = clbspots_i[occbspotid].first;

#ifdef DEBUG
        // safety check
        const int colelelid = intactdis_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",clnumbond_i,elegid,myrank_);
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->gElement(elegid));

        LINALG::Matrix<3,1> bbspotpos;
        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        ele->LocationVector(*intactdis_,lm,lmowner,lmstride);
        // get current displacements
        std::vector<double> eledisp(lm.size());
        DRT::UTILS::ExtractMyValues(*iadiscolnp,eledisp,lm);
        // extract bins from discretization after checking on existence

        // get current position of filament binding spot
        ele->GetPosOfBindingSpot(bbspotpos,eledisp,clbspots_i[occbspotid].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        // update position of single bond crosslinker to new position
        std::vector<double> newpos(3,0.0);
        for(int dim=0; dim<3; ++dim)
          newpos[dim] = bbspotpos(dim);
        crosslinker_i->SetPos(newpos);

        // that is it
        break;
      } // clnumbond_i = 1

      // crosslinker has two bonds (cl gets current mid position between the filament
      // binding spot it is attached to)
      case 2:
      {
        // -----------------------------------------------------------------
        // partner one
        // -----------------------------------------------------------------
        int elegid = clbspots_i[0].first;

#ifdef DEBUG
        // safety check
        int colelelid = intactdis_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is not \nghosted/owned on proc %i (owner of crosslinker)",clnumbond_i,elegid,myrank_);
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->gElement(elegid));
        LINALG::Matrix<3,1> bbspotposone;

        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        ele->LocationVector(*intactdis_,lm,lmowner,lmstride);
        // get current displacements
        std::vector<double> eledispone(lm.size());
        DRT::UTILS::ExtractMyValues(*iadiscolnp,eledispone,lm);

        // get current position of filament binding spot
        ele->GetPosOfBindingSpot(bbspotposone,eledispone,clbspots_i[0].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        // -----------------------------------------------------------------
        // partner two
        // -----------------------------------------------------------------
        elegid = clbspots_i[1].first;

#ifdef DEBUG
        // safety check
        colelelid = intactdis_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",clnumbond_i,elegid,myrank_);
#endif

        ele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->gElement(elegid));
        LINALG::Matrix<3,1> bbspotpostwo;

        ele->LocationVector(*intactdis_,lm,lmowner,lmstride);
        // get current displacements
        std::vector<double> eledisptwo(lm.size());
        DRT::UTILS::ExtractMyValues(*iadiscolnp,eledisptwo,lm);

        // get current position of filament binding spot
        ele->GetPosOfBindingSpot(bbspotpostwo,eledisptwo,clbspots_i[1].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        //cl gets current mid position between the filament binding spot it is attached to
        clpos_i.Update(0.5,bbspotposone,0.5,bbspotpostwo);

        // update position of single bond crosslinker to new position
        std::vector<double> newpos(3,0.0);
        for(int dim=0; dim<3; ++dim)
          newpos[dim] = clpos_i(dim);
        crosslinker_i->SetPos(newpos);

        // that is it
        break;
      }// clnumbond_i = 2

      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", clnumbond_i);
        break;
      }

    } // switch numbonds
  } // loop over all row linker

  // that is it
  return;

} //DiffuseCrosslinker()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateCrosslinking()
{
  CheckInit();

  // measure time for linking and redistribution
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::UpdateCrosslinking");

  // -------------------------------------------------------------------------
  // 1) crosslinker diffusion:
  //    - according to browninan dyn for free cl
  //    - according to beams for single and double bonded
  // -------------------------------------------------------------------------
  DiffuseCrosslinker();

  // -------------------------------------------------------------------------
  // 2) update parallel distribution for bindis and intacdis as crosslinker
  //    and beams have moved during the timestep. As we call export in the
  //    following as well, ghosted crosslinker and beams get current information
  //    regarding their binding status, position and so on.
  // -------------------------------------------------------------------------
  // get current beam positions based on maps of intactdis
  UpdateDofMapOfVector(intactdis_, ia_disnp_, GState().GetMutableDisNp());

  // todo: do we want this ? (how to choose between (re-) and partitioning?)
  // do a complete new weighted partitioning every 100th step
  int dlbevery = eval_statmech_ptr_->GetDataSMDynPtr()->DynLoadBalanceEvery();
  if(dlbevery)
  {
    if(GState().GetStepN()%dlbevery == 0 and numproc_ != 1)
      UpdateBinStrategy(false,true);
  }
  else
  {
    // transfer beams and crosslinker
    UpdateBinStrategy();
  }

  // erase temporary data container as both distributions as well the col data
  // might have changed
  crosslinker_data_.clear();
  beam_data_.clear();

  // -------------------------------------------------------------------------
  // 3) now we manage binding events, this includes:
  //    - find potential binding events on each myrank
  //    - make a decision by asking other procs
  //    - set bonds and adapt states accordingly
  // -------------------------------------------------------------------------
  // intended bonds row cl to row ele
  // note: this map is also not const, as another proc may get the permission
  // to set my row crosslink to a row element of his, so we need to erase
  // this entry out of the map
  std::map<int, BindEventData > mybonds; // key is clgid
  // intended bond col crosslinker to row element
  std::map<int, std::vector<BindEventData> > undecidedbonds; // key is owner!=myrank

  // fill binding event maps
  LookForBindingEvents(mybonds,undecidedbonds);

  // bind events where myrank only owns the elements, cl are taken care of by
  // their owner
  std::map<int, BindEventData > myelebonds; // key is clgid
  // now each row owner of a linker gets requests, makes a random decision and
  // informs back its requesters
  ManageBindingInParallel(mybonds,undecidedbonds,myelebonds);

  // actual update of binding states is done here
  BindCrosslinker(mybonds,myelebonds);

  // -------------------------------------------------------------------------
  // 4) unbinding events if probability check is passed
  // -------------------------------------------------------------------------
  UnBindCrosslinker();

  // that is it
  return;

} // UpdateCrosslinking()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::LookForBindingEvents(
  std::map<int, BindEventData >&              mybonds,
  std::map<int, std::vector<BindEventData> >& undecidedbonds)
{
  CheckInit();

  // measure time of interaction search crosslinker to beam
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::InteractionSearchAndBinding");

  // gather data for all column crosslinker initially
  const int numcolcl = bindis_->NumMyColNodes();
  crosslinker_data_.resize(numcolcl);
  PreComputeCrosslinkerData(numcolcl);

  // gather data for all column beams
  const int numcolbeams = intactdis_->NumMyColElements();
  beam_data_.resize(numcolbeams);
  PreComputeBeamData(numcolbeams);

  /* -------------------------------------------------------------------------
   note: only the owner of a beam element is allowed to change the status of
   of a binding spot. Therefore we utilize the one layer ghosting around bins
   containing a crosslinker and the ghosting around bins that are touched
   by a row element (this can lead to two layer ghosting) of a proc. Thus we
   exclude the binding of two crosslinker on different procs on the same
   binding spot without loosing any potential interaction.
   To ensure that no crosslinker is bonded to often but still totally random over
   all procs, each binding event of a col crosslinker to a row element needs to
   be communicated to the crosslinker owner, he randomly decides who is allowed
   to bind, sets the according stuff for the cl, and then informs back the
   requesting procs so they can set the stuff for the elements.
   As no proc on his own can decide whether a crosslink should be set, two
   binding events for one crosslinker in one time step are excluded (for this
   the proc must be sure that a crosslink is set as the binding range is
   different for a single bonded crosslinker compared to a free one)
  *  \author J. Eichinger
   -------------------------------------------------------------------------*/
  // store bins, which have already been examined
  std::set<int> examinedbins;
  // loop over all column crosslinker in random order
  // create random order of indices
  std::vector<int> rordercolcl = STATMECH::UTILS::Permutation(numcolcl);
  std::vector<int>::const_iterator icl;
  for(icl=rordercolcl.begin(); icl!=rordercolcl.end(); ++icl)
  {
    DRT::Node *currcrosslinker = bindis_->lColNode(*icl);

    // get bin that contains this crosslinker (can only be one)
    DRT::Element* CurrentBin = currcrosslinker->Elements()[0];
    const int currbinId = CurrentBin->Id();

    // if a bin has already been examined --> continue with next crosslinker
    if(examinedbins.find(currbinId) != examinedbins.end())
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins.insert(currbinId);

    // get neighboring bins
    // note: interaction distance cl to beam needs to be smaller than bin size
    std::vector<int> neighboring_binIds;
    neighboring_binIds.reserve(27);
    // do not check on existence here -> shifted to GetBinContent
    particlealgo_->GetNeighbouringBinGids(currbinId,neighboring_binIds);

    // get set of neighbouring beam elements (i.e. elements that somehow touch nb bins)
    // as explained above, we only need row elements
    std::set<DRT::Element*> neighboring_beams;
    particlealgo_->GetBinContent(neighboring_beams,bin_beamcontent_,neighboring_binIds,true);

    // get all crosslinker in current bin
    DRT::Node **ClInCurrentBin = CurrentBin->Nodes();
    const int numcrosslinker = CurrentBin->NumNode();

    // obtain random order in which crosslinker are addressed
    std::vector<int> randorder = STATMECH::UTILS::Permutation(numcrosslinker);

    // loop over all crosslinker in CurrentBin in random order
    std::vector<int>::const_iterator randcliter;
    for(randcliter=randorder.begin(); randcliter!=randorder.end(); ++randcliter)
    {
      // get random crosslinker in current bin
      DRT::Node *crosslinker_i = ClInCurrentBin[*randcliter];
      // get all potential binding events on myrank
      PrepareBinding(crosslinker_i,neighboring_beams,mybonds,undecidedbonds);

    } // loop over all crosslinker of current bin

  } // loop over all col crosslinker

  // that is it
  return;

} //LoodForBindingEvents()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PrepareBinding(
  DRT::Node*                                  crosslinker_i,
  const std::set<DRT::Element*>&              neighboring_beams,
  std::map<int, BindEventData >&              mybonds,
  std::map<int, std::vector<BindEventData> >& undecidedbonds)
{
  CheckInit();

  // -------------------------------------------------------------------------
  // get precomputed data of crosslinker i
  // -------------------------------------------------------------------------
  // note: at this stage we are not yet setting any bond or changing any status,
  // therefore everything is constant
  const int clcollid = crosslinker_i->LID();
  const CrosslinkerData& cldata_i = crosslinker_data_[clcollid];
  const LINALG::Matrix<3,1>& clpos_i = cldata_i.clpos;
  const int& clnumbond_i = cldata_i.clnumbond;
//  const Teuchos::RCP<MAT::CrosslinkerMat>& clmat_i = cldata_i.clmat;
  const int& clowner_i = cldata_i.clowner;

  // -------------------------------------------------------------------------
  // 1. criterion: in case crosslinker is double bonded, we can leave here
  // -------------------------------------------------------------------------
  if(clnumbond_i==2) return;

  // spherical shell around crosslinker center where binding spots have to lie
  // for potential binding event
  // todo: this needs to go crosslinker material
  double rmin = (eval_statmech_ptr_->GetDataSMDynPtr()->RLink()
               - eval_statmech_ptr_->GetDataSMDynPtr()->DeltaRLink()) / 2.0;
  double rmax = (eval_statmech_ptr_->GetDataSMDynPtr()->RLink()
               + eval_statmech_ptr_->GetDataSMDynPtr()->DeltaRLink()) / 2.0;
  // probability with which a crosslinker is established between crosslink
  // molecule and neighbor binding spot
  double plink = 1.0 - exp( -(*GState().GetDeltaTime())[0]*eval_statmech_ptr_->GetDataSMDynPtr()->KOnStart() );

  // -------------------------------------------------------------------------
  // look for potential interaction of crosslinker i and a binding spot of an
  // element, i.e. distance \Delta = R +/- \Delta R
  // -------------------------------------------------------------------------
  // loop over all neighboring beam elements in random order (keep in mind
  // we are only looping over row elements)
  std::vector<DRT::Element*> beamvec(neighboring_beams.begin(),neighboring_beams.end());
  const int numbeams = beamvec.size();
  std::vector<int> randorder = STATMECH::UTILS::Permutation(numbeams);
  std::vector<int> ::const_iterator randiter;
  for(randiter=randorder.begin(); randiter!=randorder.end();  ++randiter)
  {
    // get neighboring (nb) beam element
    DRT::ELEMENTS::Beam3Base* nbbeam =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(beamvec.at(*randiter));

#ifdef DEBUG
      if(nbbeam == NULL)
        dserror("Dynamic cast to beam3base failed");
#endif

    // -----------------------------------------------------------------------
    // get pre computed data of current nbbeam
    // -----------------------------------------------------------------------
    const BeamData& beamdata_i = beam_data_[nbbeam->LID()];
    const std::map<int, std::vector<LINALG::Matrix<3,1> > >& bbspotdofs_i = beamdata_i.bbspotdofs;
    const std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;
//    const int& bowner_i = beamdata_i.bowner;

    // loop over all binding spots of current element in random order
    std::vector<int> randbspot = STATMECH::UTILS::Permutation(bbspotstatus_i.size());
    std::vector<int> ::const_iterator rbspotiter;
    for(rbspotiter=randbspot.begin(); rbspotiter!=randbspot.end(); ++rbspotiter)
    {
      // get local number of binding spot in element
      const int locnbspot = *rbspotiter;

      // -----------------------------------------------------------------------
      // we are now doing some checks if a binding event could happen
      // -----------------------------------------------------------------------
      {
        // 2. criterion:
        // first check if binding spot is free, if not, check next bspot on curr ele
        // note: bspotstatus in bonded case holds cl gid, otherwise -1 (meaning free)
        if(bbspotstatus_i.at(locnbspot) != -1)
          continue;

        // get current position of free binding spot
        const LINALG::Matrix<3,1> currbbspos = bbspotdofs_i.at(locnbspot)[0];

        // compute distance between current binding spot and center (if free) or one
        // end (if single bonded) of crosslinker i
        static LINALG::Matrix<3,1> dist_vec;
        dist_vec.Update(1.0, currbbspos, -1.0, clpos_i);
        const double distance = dist_vec.Norm2();

        // 3. criterion:
        // if binding spot not in binding range, continue with next binding spot
        if(distance>rmax || distance<rmin)
          continue;

        // 4. criterion:
        // a crosslink is set if and only if it passes a probability check
        if(DRT::Problem::Instance()->Random()->Uni() > plink)
          continue;
      }

      // ---------------------------------------------------------------------
      // if we came this far, we can add this potential binding event to its
      // corresponding map
      // ---------------------------------------------------------------------
      BindEventData bindeventdata;
      bindeventdata.clgid = crosslinker_i->Id();
      bindeventdata.elegid = nbbeam->Id();
      bindeventdata.bspotlocn = locnbspot;
      bindeventdata.requestproc = myrank_;
      // this is default, is changed if owner of cl has something against it
      bindeventdata.permission = true;

      // in case myrank is owner, we add it to the mybonds map
      if(clowner_i == myrank_)
      {
        mybonds[bindeventdata.clgid] = bindeventdata;
      }
      // myrank is not owner, we add it to the map of events that need to be
      // communicated to make a decision
      else
      {
        undecidedbonds[clowner_i].push_back(bindeventdata);
      }

      // as we allow only one binding event for each cl in one time step,
      // we are done here, if we made it so far (i.e met 1., 2., 3. and 4.) criterion)
      return;

    } // loop over all binding spots

  }// loop over neighboring row eles

  // that is it
  return;

} //PrepareBinding()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BindCrosslinker(
  std::map<int, BindEventData >& mybonds,
  std::map<int, BindEventData >& myelebonds)
{
  CheckInit();

  // time needed to realize binding events
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::DoBinding");

  // store data of newly double bonded cl on proc = map key
  std::map<int, std::vector<DoubleBondedCl> > dbondcltosend;
  // this map holds all data for newly double bonded crosslinker on myrank
  // map key is crosslinker gid as we need to be able to uniquely address
  // one entry and each "element" that is set on a proc needs to have the same
  // unique key so we can delete it on all procs in case of a crosslinker
  // unbinding event
  std::map<int, DoubleBondedCl> mynewdbondcl;

  /* -------------------------------------------------------------------------
   now have two distinct maps of binding events on each proc, depending
   on ownership of crosslinker and elements myrank has different tasks:
   - mybonds: myrank takes care of crosslinker and (most) elements
   - myelebonds: myrank takes care of elements
   within those maps, different treatment is necessary for free and single
   bonded linker
   -------------------------------------------------------------------------*/

  // -------------------------------------------------------------------------
  // 1) mybonds
  // -------------------------------------------------------------------------
  {
    std::map<int, BindEventData>::const_iterator cliter;
    for(cliter=mybonds.begin(); cliter!=mybonds.end(); ++cliter)
    {
      // get binding event data
      BindEventData binevdata = cliter->second;

#ifdef DEBUG
      // safety check
      if(cliter->first!=binevdata.clgid)
        dserror("Map key does not match crosslinker gid of current binding event.");
#endif

      // get current linker
      const int clcollid = bindis_->NodeColMap()->LID(cliter->first);
      CROSSLINKING::CrosslinkerNode *crosslinker_i =
          dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lColNode(clcollid));

      // get crosslinker data
      CrosslinkerData& cldata_i = crosslinker_data_[clcollid];
      LINALG::Matrix<3,1>& clpos_i = cldata_i.clpos;
      std::vector<std::pair<int, int> >& clbspots_i = cldata_i.clbspots;
      int& clnumbond_i = cldata_i.clnumbond;

#ifdef DEBUG
      // safety check
      const int& clowner_i = cldata_i.clowner;
      if(clowner_i!=myrank_)
        dserror("Only row owner of crosslinker is changing its status");
#endif

      // get beam data
      const int colelelid = intactdis_->ElementColMap()->LID(binevdata.elegid);

#ifdef DEBUG
      // safety check
      if(colelelid<0)
        dserror("Binding element partner of current row crosslinker is not ghosted, "
                "this must be the case though.");
#endif

      // beam element i
      DRT::ELEMENTS::Beam3Base* ele_i =
          dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->lColElement(colelelid));
      BeamData& beamdata_i = beam_data_[colelelid];
      const std::map<int, std::vector<LINALG::Matrix<3,1> > >& bbspotdofs_i = beamdata_i.bbspotdofs;
      std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;
      const int& bowner_i = beamdata_i.bowner;

      // -------------------------------------------------------------------------
      // different treatment according to number of bonds crosslinker had before
      // this binding event
      // -------------------------------------------------------------------------
      switch(clnumbond_i)
      {
        // crosslinker with zero bonds before this binding event
        case 0:
        {
          // -----------------------------------------------------------------
          // update crosslinker status
          // -----------------------------------------------------------------
          // store gid and bspot local number of this element, first binding spot
          // always bonded first
          clbspots_i[0].first  = binevdata.elegid;
          clbspots_i[0].second = binevdata.bspotlocn;
          crosslinker_i->ClData()->SetClBSpotStatus(clbspots_i);

          // update number of bonds
          clnumbond_i = 1;
          crosslinker_i->ClData()->SetNumberOfBonds(clnumbond_i);

          // update position
          clpos_i = bbspotdofs_i.at(binevdata.bspotlocn)[0];
          std::vector<double> newpos(3,0.0);
          for(int dim=0; dim<3; ++dim)
            newpos[dim] = clpos_i(dim);
          crosslinker_i->SetPos(newpos);

#ifdef DEBUG
          // safety check
          if(not (clbspots_i[1].first < 0))
            dserror("Numbond does not fit to clbspot vector.");
#endif

          // -----------------------------------------------------------------
          // update beam status
          // -----------------------------------------------------------------
          // store crosslinker gid in status of beam binding spot if myrank
          // is owner of beam
          if(bowner_i == myrank_)
          {
            bbspotstatus_i.at(binevdata.bspotlocn) = binevdata.clgid;
            ele_i->SetBindingSpotStatus(bbspotstatus_i);
          }

          break;
        } // numbond = 0

        // crosslinker with one bond before this binding event
        case 1:
        {
          // -----------------------------------------------------------------
          // check which of the two clbspot is free
          // -----------------------------------------------------------------
          int freebspotid = 1;
          int occbspotid = 0;
          // check, which clbspot is free
          if(clbspots_i[0].first < 0)
          {
            freebspotid = 0;
            occbspotid  = 1;
          }

#ifdef DEBUG
          // safety check
          if(clbspots_i[occbspotid].first < 0)
            dserror("clnumbond_i doesn't match clbspots_i");
#endif

          // get owner of element of already existing bond (that was set before
          // this time step)
          const int oldbondpartnergid = clbspots_i[occbspotid].first;
          const int oldbondpartnerowner = intactdis_->gElement(oldbondpartnergid)->Owner();

#ifdef DEBUG
          // safety check
          const int colelelid = intactdis_->ElementColMap()->LID(oldbondpartnergid);
          if(colelelid<0)
            dserror("Binding element partner of current row crosslinker is not ghosted, but this must be the case");
#endif

          // -----------------------------------------------------------------
          // update crosslinker status
          // -----------------------------------------------------------------
          // store gid and bspot local number of this element
          clbspots_i[freebspotid].first = binevdata.elegid;
          clbspots_i[freebspotid].second = binevdata.bspotlocn;
          crosslinker_i->ClData()->SetClBSpotStatus(clbspots_i);

          // update number of bonds
          clnumbond_i = 2;
          crosslinker_i->ClData()->SetNumberOfBonds(clnumbond_i);

          // update position
          std::vector<double> newpos(3,0.0);
          for(int dim=0; dim<3; ++dim)
          {
            // newpos = 0.5*( pos new bbspot partner + old bbspot partner position(=cl pos))
            newpos[dim] = (bbspotdofs_i.at(clbspots_i[freebspotid].second)[0](dim) + clpos_i(dim))/2.0;
            clpos_i(dim ) = newpos[dim];
          }
          crosslinker_i->SetPos(newpos);

          // create double bond cl data
          DoubleBondedCl dbondcl;
          dbondcl.id = binevdata.clgid;
          dbondcl.eleone = clbspots_i[freebspotid];
          dbondcl.eletwo = clbspots_i[occbspotid];

          // first check if myrank is owner of element of current binding event
          // (additionally to being owner of cl)
          if(bowner_i == myrank_)
          {
            // -----------------------------------------------------------------
            // update beam data
            // -----------------------------------------------------------------
            // store crosslinker gid in status of beam binding spot
            bbspotstatus_i.at(binevdata.bspotlocn) = binevdata.clgid;
            ele_i->SetBindingSpotStatus(bbspotstatus_i);

            // insert pair in mypairs
            mynewdbondcl[dbondcl.id] = dbondcl;

            // -----------------------------------------------------------------
            // in this case we need to communicate again, as each proc that is
            // part of a binding event that leads to a new element needs to
            // either own or ghost this element
            // -----------------------------------------------------------------
            if(oldbondpartnerowner != myrank_)
            {
              /*
               * 1| 1__|2
               * 1|    |2
               * legend: | = beam; __= cl; 2 = owner; 1=myrank
               */
              dbondcltosend[oldbondpartnerowner].push_back(dbondcl);
            }
            // else: we do not need to communicate as myrank does all that is
            // necessary for this binding event
            /*
             * 1| 1__|1
             * 1|    |1
             * legend: | = beam; __= cl; 1=myrank
             */

          } // binevdata.requestproc==myrank_

          // -------------------------------------------------------------------
          // if myrank is just owner of second clbspot partner (i.e. bond that
          // was set before this time step)
          // -------------------------------------------------------------------
          else if (oldbondpartnerowner == myrank_)
          {
            /*
             * 1| 2__|2
             * 1|    |2
             * legend: | = beam; __= cl; 2 = owner; 1=myrank
             */
            // insert pair in mypairs
            mynewdbondcl[dbondcl.id] = dbondcl;

            // note: in this case we do not need to communicate this element, the
            // respective owner of the new bond partner will add this pair on his own

          } // ownerotherbspot == myrank_

          // -------------------------------------------------------------------
          // in this case myrank is owner of crosslinker but does not own none of
          // the two involved elements
          // -------------------------------------------------------------------
          else
          {
            // todo: do we need to add anything on myrank too in this case

            if(bowner_i != oldbondpartnerowner)
            {
              /*
               * 3| 1__|2
               * 3|    |2
               * legend: | = beam; __= cl; 2,3 = owner; 1=myrank
               */
              // in this case we need to communicate and inform the proc with the
              // old bond (not part of current binding event) that it needs to
              // compute the interaction of this now double bond crosslinker
              dbondcltosend[oldbondpartnerowner].push_back(dbondcl);;
            }

            // if we have the same owner in both partner eles, we neither need to
            // communicate as myelebonds list on other proc adds ele for us
            /*
             * 2| 1__|2
             * 2|    |2
             * legend: | = beam; __= cl; 2 = owner; 1=myrank
             */

          } // binevdata.requestproc != myrank_ and ownerotherbspot != myrank_

          break;
        } // numbond = 1

        default:
        {
          dserror("You should not be here, crosslinker has unrealistic number "
                  "%i of bonds.", clnumbond_i);
          break;
        }

      } // switch number of clnumbond
    }// loop through mybonds

  } // mybonds


  // -------------------------------------------------------------------------
  // 2) myelebonds (myrank only owner of current binding partner ele)
  // -------------------------------------------------------------------------
  {
    /*
      * 1| 2__|2  or  1| 3__|2  or  1| 2__|1
      * 1|    |2      1|    |2      1|    |1
      * legend: | = beam; __= cl; 2,3 = owner; 1=myrank
      */
    // loop through all binding events
    std::map<int, BindEventData>::const_iterator cliter;
    for(cliter=myelebonds.begin(); cliter!=myelebonds.end(); ++cliter)
    {
      // get binding event data
      BindEventData binevdata = cliter->second;

      // ---------------------------------------------------------------------
      // get linker data
      // ---------------------------------------------------------------------
      const int clcollid = bindis_->NodeColMap()->LID(cliter->first);

#ifdef DEBUG
      // safety check
      if(clcollid<0)
        dserror("Crosslinker needs to be ghosted, but this isn't the case.");
#endif

      CrosslinkerData& cldata_i = crosslinker_data_[clcollid];
//      LINALG::Matrix<3,1>& clpos_i = cldata_i.clpos;
      std::vector<std::pair<int, int> >& clbspots_i = cldata_i.clbspots;
      int& clnumbond_i = cldata_i.clnumbond;

      // ---------------------------------------------------------------------
      // get beam data
      // ---------------------------------------------------------------------
      const int colelelid = intactdis_->ElementColMap()->LID(binevdata.elegid);

#ifdef DEBUG
      // safety check
      if(colelelid<0)
        dserror("element with gid %i not ghosted on proc %i",binevdata.elegid,myrank_);
#endif

      // get beam element of current binding event
      DRT::ELEMENTS::Beam3Base* ele_i =
          dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->lColElement(colelelid));
      BeamData& beamdata_i = beam_data_[colelelid];
//      const std::map<int, std::vector<LINALG::Matrix<3,1> > >& bbspotdofs_i = beamdata_i.bbspotdofs;
      std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;

#ifdef DEBUG
      // safety check
      const int& bowner_i = beamdata_i.bowner;
      if(bowner_i!=myrank_)
        dserror("Only row owner of element is allowed to change its status");
      // safety checks
      const int& clowner_i = cldata_i.clowner;
      if(clowner_i==myrank_)
        dserror("myrank should not be owner of this crosslinker");
#endif

      // different treatment according to number of bonds crosslinker has before
      // this binding event
      switch(clnumbond_i)
      {
        case 0:
        {
          // -----------------------------------------------------------------
          // update beam data
          // -----------------------------------------------------------------
          // store crosslinker gid in status of beam binding spot
          bbspotstatus_i.at(binevdata.bspotlocn) = binevdata.clgid;
          ele_i->SetBindingSpotStatus(bbspotstatus_i);

          break;
        } // clnumbond = 0

        case 1:
        {
          // -----------------------------------------------------------------
          // update beam data
          // -----------------------------------------------------------------
          // store crosslinker gid in status of beam binding spot
          bbspotstatus_i.at(binevdata.bspotlocn) = binevdata.clgid;
          ele_i->SetBindingSpotStatus(bbspotstatus_i);

          // create double bond cl data
          DoubleBondedCl dbondcl;
          dbondcl.id = binevdata.clgid;
          dbondcl.eleone = clbspots_i[0];
          dbondcl.eletwo = clbspots_i[1];
          mynewdbondcl[dbondcl.id] = dbondcl;

          break;
        } // clnumbond_i = 1

        default:
        {
          dserror("You should not be here, crosslinker has to many bonds.");
          break;
        }

      } // switch number of bonds

    } // loop over myelebonds

  } // myelebonds

  // tell each proc which of the crosslinker that are bonded to one of its row
  // elements are now double bonded without him knowing (yet)
  CommunicateDoubleBondedCl(mynewdbondcl,dbondcltosend);

  //todo: HOW TO INSERT NEW DOUBLE BONDS IN NORMAL LIST

  // that is it
  return;

} // BindCrosslinker()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnBindCrosslinker()
{
  CheckInit();

  // time needed to communicate new crosslinker elements
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::CrosslinkerUnBinding");

  // get current off-rate for crosslinkers
  // todo: this needs to go somewhere else
  double koff = 0;
  double starttime = eval_statmech_ptr_->GetDataSMDynPtr()->ActionTime()->at(1);
  const double timenp = GState().GetTimeNp();
  const double dt = (*GState().GetDeltaTime())[0];

  if (timenp <= starttime || (timenp>starttime && fabs(starttime)<dt/1e4))
    koff = eval_statmech_ptr_->GetDataSMDynPtr()->KOffStart();
  else
    koff = eval_statmech_ptr_->GetDataSMDynPtr()->KOffEnd();

  // probability with which a crosslink breaks up in the current time step
  double p_unlink = 1.0 - exp(-dt * koff);

  // --------------------------------------------------------------------------
  // vectors that need to be filled
  // --------------------------------------------------------------------------
  // data containing information about elements that need to be updated and
  // crosslinker that are not double bonded anymore
  std::map<int, std::vector<UnBindEventData> > sendunbindevent;
  // crosslinker "eles" that need to be deleted on myrank
  std::vector<UnBindEventData> myrankunbindevent;

  // loop over all row linker (in random order) and dissolve bond if probability
  // criterion is met
  /* note: we loop over all row crosslinker, i.e. myrank needs to update all
   * crosslinker information. As it possible that a row crosslinker is linked
   * to col element, we potentially need to communicate if such an element
   * needs to be updated and/or if a crosslinker element needs to deleted
   */
  const int numrowcl = bindis_->NumMyRowNodes();
  std::vector<int> rorderrowcl = STATMECH::UTILS::Permutation(numrowcl);
  std::vector<int>::const_iterator rowcli;
  for(rowcli=rorderrowcl.begin(); rowcli!=rorderrowcl.end(); ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(*rowcli));

    // get linker data
    const int clgid = crosslinker_i->Id();
    const int clcollid = bindis_->NodeColMap()->LID(crosslinker_i->Id());
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];
    LINALG::Matrix<3,1>& clpos_i = cldata_i.clpos;
    std::vector<std::pair<int, int> >& clbspots_i = cldata_i.clbspots;
    int& clnumbond_i = cldata_i.clnumbond;
//    const int& clowner_i = cldata_i.clowner;

    // different treatment according to number of bonds a crosslinker already has
    switch(clnumbond_i)
    {
      // crosslinker has zero bonds
      case 0:
      {
        // nothing to do here
        break;
      } // clnumbond_i = 0

      // crosslinker has one bond
      case 1:
      {
        // -------------------------------------------------------------------
        // if probability criterion is not met, we are done here
        // -------------------------------------------------------------------
        if (DRT::Problem::Instance()->Random()->Uni() > p_unlink)
          break;

        // -------------------------------------------------------------------
        // dissolve bond
        // -------------------------------------------------------------------
        // check, which clbspot is occupied
        int occbspotid = 0;
        if(clbspots_i[0].first < 0)
          occbspotid  = 1;

#ifdef DEBUG
        // safety check
        if(clbspots_i[occbspotid].first < 0)
          dserror("clnumbond_i doesn't match clbspots_i");
#endif

        // store unbinding event data
        UnBindEventData unbindevent;
        // crosslinker has not been double bonded
        unbindevent.id = -1;
        unbindevent.eletoupdate = clbspots_i[occbspotid];

        // owner of beam
        const int beamowner =
            intactdis_->gElement(unbindevent.eletoupdate.first)->Owner();
        // check who needs to update the element status
        if(beamowner == myrank_)
          myrankunbindevent.push_back(unbindevent);
        else
          sendunbindevent[beamowner].push_back(unbindevent);

        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // update binding status of linker
        clbspots_i[occbspotid].first = -1;
        clbspots_i[occbspotid].second = -1;
        crosslinker_i->ClData()->SetClBSpotStatus(clbspots_i);

        // update number of bonds
        clnumbond_i = 0;
        crosslinker_i->ClData()->SetNumberOfBonds(clnumbond_i);

        // update position of crosslinker: generate vector in random direction
        // of length R_LINK to "reset" crosslink molecule position: it may now
        // reenter or leave the bonding proximity
        // todo: does this make sense?
        LINALG::Matrix<3,1> cldeltapos_i;
        std::vector<double> randunivec(3);
        int count = 3;
        DRT::Problem::Instance()->Random()->Uni(randunivec, count);
        for (int dim=0; dim<3; ++dim)
          cldeltapos_i(dim) = randunivec[dim];
        cldeltapos_i.Scale(eval_statmech_ptr_->GetDataSMDynPtr()->RLink() / cldeltapos_i.Norm2());
        // update position in cldata
        clpos_i.Update(1.0,cldeltapos_i,1.0);
        // real update of position
        std::vector<double> newpos(3,0.0);
        for(int dim=0; dim<3; ++dim)
          newpos[dim] = clpos_i(dim);
        crosslinker_i->SetPos(newpos);

        break;
      } // clnumbond_i = 1

      // crosslinker has two bonds
      case 2:
      {
        // get id of freed and still occupied bspot
        int freedbspotid = -1;
        int stayocc = -1;
        // loop through crosslinker bonds in random order
        std::vector<int> ro = STATMECH::UTILS::Permutation(clnumbond_i);
        std::vector<int>::const_iterator clbspotiter;
        for(clbspotiter=ro.begin(); clbspotiter!=ro.end(); ++clbspotiter)
        {
          // if probability criterion isn't met, go to next spot
          if (DRT::Problem::Instance()->Random()->Uni() > p_unlink)
            continue;

          // get id of freed and still occupied bspot
          freedbspotid = *clbspotiter;
          stayocc = 0;
          if(freedbspotid == 0)
            stayocc = 1;

          // -----------------------------------------------------------------
          // dissolve bond
          // -----------------------------------------------------------------
          // owner of unbonded beam element
          const int freedbeamowner =
              intactdis_->gElement(clbspots_i[freedbspotid].first)->Owner();
          // owner of beam that stays connected to crosslinker
          const int occbeamowner =
              intactdis_->gElement(clbspots_i[stayocc].first)->Owner();
          // initialize two versions of unbinding events
          // update and deletion
          UnBindEventData unbindeventupdel;
          unbindeventupdel.id = clgid;
          unbindeventupdel.eletoupdate = clbspots_i[freedbspotid];
          // just deletion
          UnBindEventData unbindeventdel;
          unbindeventdel.id = clgid;
          unbindeventdel.eletoupdate.first = -1;

          // check different configurations
          if(freedbeamowner == myrank_ && occbeamowner == myrank_)
          {
           /*
           * 1| 1__|1
           * 1|    |1
           * legend: | = beam; __= cl; 1 = myrank
           */
            // myrank updates beam element and deletes crosslinker element
            myrankunbindevent.push_back(unbindeventupdel);
          }
          else if(freedbeamowner == myrank_)
          {
            /*
            * 1| 1__|2
            * 1|    |2
            * legend: | = beam; __= cl; 2 = owner; 1=myrank
            */
            // myrank updates beam element and deletes crosslinker element
            myrankunbindevent.push_back(unbindeventupdel);

            // in this case the proc only needs to delete the crosslinker element
            sendunbindevent[occbeamowner].push_back(unbindeventdel);
          }
          else if(occbeamowner == myrank_)
          {
            /*
            * 2| 1__|1
            * 2|    |1
            * legend: | = beam; __= cl; 2 = owner; 1=myrank
            */
            // myrank only deletes crosslinker element
            myrankunbindevent.push_back(unbindeventdel);

            // freedbeamowner updates beam element and deletes crosslinker element
            sendunbindevent[freedbeamowner].push_back(unbindeventupdel);
          }
          else if (freedbeamowner == occbeamowner)
          {
            /*
            * 2| 2__|2
            * 2|    |2
            * legend: | = beam; __= cl; 2 = owner, 1=myrank
            */
            // freedbeamowner updates beam element and deletes crosslinker element
            sendunbindevent[freedbeamowner].push_back(unbindeventupdel);
          }
          else
          {
            /*
            * 2| 2__|3
            * 2|    |3
            * legend: | = beam; __= cl; 2,3 = owner, 1=myrank
            */
            // freedbeamowner updates beam element and deletes crosslinker element
            sendunbindevent[freedbeamowner].push_back(unbindeventupdel);

            // otherbeamowner only needs to delete the crosslinker element
            sendunbindevent[occbeamowner].push_back(unbindeventdel);
          }

          // -----------------------------------------------------------------
          // update crosslinker status
          // -----------------------------------------------------------------
          // reset binding status of freed crosslinker binding spot
          clbspots_i[freedbspotid].first = -1;
          clbspots_i[freedbspotid].second = -1;
          crosslinker_i->ClData()->SetClBSpotStatus(clbspots_i);

          // update number of bonds
          clnumbond_i = 1;
          crosslinker_i->ClData()->SetNumberOfBonds(clnumbond_i);

          // update postion
          const int collidoccbeam =
              intactdis_->ElementColMap()->LID(clbspots_i[stayocc].first);

#ifdef DEBUG
          // safety check
          if(collidoccbeam<0)
            dserror("element with gid %i not ghosted on proc %i",clbspots_i[stayocc].first,myrank_);
#endif

          BeamData& beamdata_i = beam_data_[collidoccbeam];
          std::map<int, std::vector<LINALG::Matrix<3,1> > > bspotdofs_i = beamdata_i.bbspotdofs;
          clpos_i = bspotdofs_i.at(clbspots_i[stayocc].second)[0];
          std::vector<double> newpos(3,0.0);
          for(int dim=0; dim<3; ++dim)
            newpos[dim] = clpos_i(dim);
          crosslinker_i->SetPos(newpos);

          // we only want to dissolve one bond per timestep, therefore we go to
          // next crosslinker if we made it so far (i.e. a bond got dissolved)
          break;

        } // loop over binding spots on crosslinker
        break;
      }// clnumbond_i = 2

      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", clnumbond_i);
        break;
      }
    } // switch numbonds
  } // loop over all crosslinker

  // communicate elements that need to be updated and elements that need to be
  // deleted
  CommunicateCrosslinkerUnbinding(sendunbindevent,myrankunbindevent);

  // loop through all unbinding events on myrank
  std::vector<UnBindEventData>::const_iterator iter;
  for(iter=myrankunbindevent.begin(); iter!=myrankunbindevent.end(); ++ iter)
  {
    // get data
    const int elegidtoupdate = iter->eletoupdate.first;
    const int bspotlocn = iter->eletoupdate.second;
    const int idtoerase = iter->id;

    // update beam data
    if(elegidtoupdate > -1)
    {
      const int colelelid = intactdis_->ElementColMap()->LID(elegidtoupdate);

#ifdef DEBUG
      // safety check
      if(colelelid<0)
        dserror("element with gid %i not ghosted on proc %i",elegidtoupdate,myrank_);
#endif

      // get beam element of current binding event
      DRT::ELEMENTS::Beam3Base* ele_i =
          dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->lColElement(colelelid));

#ifdef DEBUG
      // safety check
      if(ele_i->Owner()!=myrank_)
        dserror("Only row owner of elements change binding status of elements");
#endif

      BeamData& beamdata_i = beam_data_[colelelid];
      std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;

      // update beam data
      bbspotstatus_i.at(bspotlocn) = -1;
      ele_i->SetBindingSpotStatus(bbspotstatus_i);
    }

#ifdef DEBUG
        // safety check
        if(not doublebondcl_.count(idtoerase))
          dserror("willing to delete not existing entry, something went wrong");
#endif

    // delete double bonded crosslinker
    doublebondcl_.erase(idtoerase);

  }

  // that is it
  return;

} //UnbindCrosslinker()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PreComputeCrosslinkerData(
  const int numcolcl)
{
  CheckInit();

  // sanity check
  if((int)crosslinker_data_.size()!= bindis_->NumMyColNodes())
    dserror("temporary crosslinker data container has wrong size");

  // -------------------------------------------------------------------------
  // loop over all column crosslinker and pre compute their data that is needed
  // multiple times in the following, therefore this gives faster access
  // note: we get references here, i.e. changes in the crosslinker status will
  // be done with the variables created here (only row owner change something,
  // this is ensured in the actual algorithm)
  // -------------------------------------------------------------------------
  for (int i=0; i<numcolcl; ++i)
  {
    // crosslinker i for which data will be collected
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lColNode(i));

#ifdef DEBUG
      if(crosslinker_i == NULL)
        dserror("Dynamic cast to CrosslinkerNode failed");
#endif

    // col lid of this crosslinker
    const int collid = crosslinker_i->LID();
    // store date of crosslinker i according to column lid
    CrosslinkerData& cldata = crosslinker_data_[collid];

    // store positions
    for(int dim=0; dim<3; ++dim)
      cldata.clpos(dim) = crosslinker_i->X()[dim];
    // get current binding spot status of crosslinker
    cldata.clbspots = crosslinker_i->ClData()->GetClBSpotStatus();
    // get number of bonds
    cldata.clnumbond = crosslinker_i->ClData()->GetNumberOfBonds();
    // get type of crosslinker (i.e. its material)
    cldata.clmat = crosslinker_i->GetMaterial();
    // get owner
    cldata.clowner = crosslinker_i->Owner();

#ifdef DEBUG
      if(crosslinker_i->NumElement() != 1)
        dserror("More than one element for this crosslinker");
#endif

  } // loop over all column crosslinker

  // that is it
  return;

} //PreComputeCrosslinkerData()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PreComputeBeamData(
  const int numcolbeams)
{
  CheckInit();

  // sanity check
  if((int)beam_data_.size()!= intactdis_->NumMyColElements())
    dserror("temporary beam data container has wrong size");

  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector(*intactdis_->DofColMap()));
  LINALG::Export(*ia_disnp_, *iadiscolnp);

  // loop over all column beams elements
  for (int i=0; i<numcolbeams; ++i)
  {
    // beam element i for which data will be collected
    DRT::ELEMENTS::Beam3Base* ele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(intactdis_->lColElement(i));

#ifdef DEBUG
      if(ele_i == NULL)
        dserror("Dynamic cast to Beam3Base failed");
#endif

    // get elelid
    const int elelid = ele_i->LID();
    // store data
    BeamData& bdata = beam_data_[elelid];

    // get element location vector and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele_i->LocationVector(*intactdis_,lm,lmowner,lmstride);
    // get current displacements
    std::vector<double> eledisp(lm.size());
    // todo: discolnp_ with correct values here?
    DRT::UTILS::ExtractMyValues(*iadiscolnp,eledisp,lm);

    // loop over all binding spots of current element
    const int numbbspot = (int)ele_i->GetBindingSpotStatus().size();
    for(int j=0; j<numbbspot; ++j)
    {
      // get current position at binding spot xi
      LINALG::Matrix<3,1> bbspotpos_j;
      ele_i->GetPosOfBindingSpot(bbspotpos_j,eledisp,j,
          *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));
      bdata.bbspotdofs[j].push_back(bbspotpos_j);
    }

    // get status of beam binding spots
    bdata.bbspotstatus = ele_i->GetBindingSpotStatus();
    // get owner
    bdata.bowner = ele_i->Owner();

  } // loop over all column beams

  // that is it
  return;

} //PreComputeBeamData()


/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ManageBindingInParallel(
  std::map<int, BindEventData >&              mybonds,
  std::map<int, std::vector<BindEventData> >& undecidedbonds,
  std::map<int, BindEventData >&              myelebonds)
{
  CheckInit();

  // time needed to make binding decision in parallel
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::BindingDecision");

  // variable for safety check
  int numrecrequest;

  // -------------------------------------------------------------------------
  // 1) each procs makes his requests and receives the request of other procs
  // -------------------------------------------------------------------------
  // store requested cl and its data
  std::map<int, std::vector<BindEventData> > requestedcl;
  {
    // -----------------------------------------------------------------------
    // send
    // -----------------------------------------------------------------------
    DRT::Exporter exporter(bindis_->Comm());
    const int length = undecidedbonds.size();
    std::vector<MPI_Request> request(length);
    int tag = 0;
    std::map<int, std::vector<BindEventData> >::iterator p;
    for(p=undecidedbonds.begin(); p!=undecidedbonds.end(); ++p)
    {
      // ---- pack data for sending -----
      std::vector<char> sdata;
      std::vector<BindEventData>::const_iterator be;
      DRT::PackBuffer data;
      for(be=p->second.begin(); be!=p->second.end(); ++be)
      {
        Pack(data,*be);
      }
      data.StartPacking();
      for(be=p->second.begin(); be!=p->second.end(); ++be)
      {
        Pack(data,*be);
      }
      swap(sdata,data());

      // unblocking send
      exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
      ++tag;
    }
    if (tag != length) dserror("Number of messages is mixed up");

    // -----------------------------------------------------------------------
    // receive
    // -----------------------------------------------------------------------

    // ---- prepare receiving procs -----
    // get number of procs from which myrank receives data
    std::vector<int> targetprocs(numproc_,0);
    std::map<int, std::vector<BindEventData> >::iterator prociter;
    for(prociter=undecidedbonds.begin(); prociter!=undecidedbonds.end(); ++prociter)
      targetprocs[prociter->first] = 1;
    // store number of messages myrank receives
    std::vector<int> summedtargets(numproc_,0);
    bindis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc_);

    numrecrequest = summedtargets[myrank_];
    for(int rec=0; rec<numrecrequest; ++rec)
    {
      std::vector<char> rdata;
      int length = 0;
      int tag = -1;
      int from = -1;
      exporter.ReceiveAny(from,tag,rdata,length);
      if (tag != 1234)
        dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

      // store received data
      std::vector<char>::size_type position = 0;
      while (position < rdata.size())
      {
        // ---- extract received data -----
        BindEventData reccldata;
        UnPack(position,rdata,reccldata);

        // create map holding all requests
        requestedcl[reccldata.clgid].push_back(reccldata);

      } // loop over request of same proc to myrank

      if (position != rdata.size())
        dserror("Mismatch in size of data %d <-> %d",(int)rdata.size(),position);

    } // loop over all requesting procs

    // -----------------------------------------------------------------------
    // wait for all communication to finish
    // -----------------------------------------------------------------------
    for (int i=0; i<length; ++i)
      exporter.Wait(request[i]);

    // note: if we have done everything correct, this should be a no time operation
    bindis_->Comm().Barrier(); // I feel better this way ;-)

  } // sending and receiving requests

  // -------------------------------------------------------------------------
  // 2) now myrank needs to decide which proc is allowed to set the requested
  // link, add it to his own list as row owner of cl sets stuff for cls, send
  // back the answers to the row ele owner and receive the decisions made for
  // its own requests:
  // - if only one proc is requesting, the link can be set
  // - if two procs are requesting or the current proc wants to set a link with
  //   a requested crosslinker, a random decision who is allowed to set the link
  //   has to be made
  // -------------------------------------------------------------------------
  std::map<int, std::vector<BindEventData> > decidedbonds;
  {
    std::map<int, std::vector<BindEventData> >::iterator cliter;
    // loop over all requested cl (note myrank is owner of these)
    for(cliter=requestedcl.begin(); cliter!=requestedcl.end(); ++cliter)
    {
      // check if myrank wants to bind this crosslinker
      const bool myrankbond = (mybonds.find(cliter->first) != mybonds.end());

      // ---------------------------------------------------------------------
      // if only one request and myrank does not want to bind this cl,
      // requesting proc gets the permission to do so
      // ---------------------------------------------------------------------
      if((int)cliter->second.size() == 1 and not myrankbond)
      {
        // we send back the permission to the relevant proc, because myrank as row
        // owner of bspot needs to set the respective stuff for the element of this
        // binding event
        // note: permission = true was send as default, so this can be sent back
        // without changes
        decidedbonds[cliter->second[0].requestproc].push_back(cliter->second[0]);

        // insert this new binding event in map of myrank, because as row owner of
        // this cl he is responsible to set the respective stuff for the crosslinker
        // of this binding event
        mybonds[cliter->first] = cliter->second[0];

        // go to next crosslinker
        continue;
      }

      // ---------------------------------------------------------------------
      // in case number of requesting procs >1 for this cl or myrank wants to
      // set it itself
      // ---------------------------------------------------------------------
      int numrequprocs = cliter->second.size();
      if(myrankbond)
        numrequprocs += 1;

      // get random proc out of affected ones
      DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);
      // todo: what if random number exactly = 1?
      int randowner = (int)floor(numrequprocs* DRT::Problem::Instance()->Random()->Uni());

      // myrank is allowed to set link
      if(myrankbond and randowner==numrequprocs-1)
      {
        // note: this means link is set between row cl and row ele on myrank,
        // all relevant information for myrank is stored in mybonds
        // loop over all requesters and store their veto
        std::vector<BindEventData>::iterator iter;
        for(iter=cliter->second.begin(); iter!=cliter->second.end(); ++iter)
        {
          iter->permission = false;
          decidedbonds[iter->requestproc].push_back(*iter);
        }
      }
      // certain requester is allowed to set the link
      else
      {
        // loop over all requesters and store veto for all requester except for one
        std::vector<BindEventData>::iterator iter;

        int counter = 0;
        for(iter=cliter->second.begin(); iter!=cliter->second.end(); ++iter)
        {
          if(randowner==counter)
          {
            // permission for this random proc
            decidedbonds[iter->requestproc].push_back(*iter);

            // erase old binding event
            if(myrankbond)
              mybonds.erase(cliter->first);

            // insert new binding event
            mybonds[cliter->first] = *iter;
          }
          else
          {
            iter->permission = false;
            decidedbonds[iter->requestproc].push_back(*iter);
          }
          counter++;
        }
      }

    } // loop over all requested crosslinks
  } // making a decision for all request

  // -------------------------------------------------------------------------
  // 3) communicate the binding decisions made on myrank, receive decisions
  //     made for its own requests and create colbondmap accordingly
  // -------------------------------------------------------------------------
  {
    // -----------------------------------------------------------------------
    // send back decisions for all requests that were made
    // -----------------------------------------------------------------------
    DRT::Exporter exporter(bindis_->Comm());
    const int length = decidedbonds.size();
    std::vector<MPI_Request> request(length);
    int tag = 0;

    // safety check
    // todo: if tested, do it only in debug mode
    if(length != numrecrequest)
      dserror("Number of received requests %i unequal to number of answers %i",numrecrequest,length);

    std::map<int, std::vector<BindEventData> >::iterator p;
    for(p=decidedbonds.begin(); p!=decidedbonds.end(); ++p)
    {
      // ---- pack data for sending -----
      std::vector<char> sdata;
      std::vector<BindEventData>::const_iterator be;
      DRT::PackBuffer data;
      for(be=p->second.begin(); be!=p->second.end(); ++be)
      {
        Pack(data,*be);
      }
      data.StartPacking();
      for(be=p->second.begin(); be!=p->second.end(); ++be)
      {
        Pack(data,*be);
      }
      swap(sdata,data());
      // unblocking send
      exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
      ++tag;
    }
    if (tag != length) dserror("Number of messages is mixed up");

    // -----------------------------------------------------------------------
    // receive
    // -----------------------------------------------------------------------
    // store requested cl and its data
    int answersize = (int)undecidedbonds.size();
    for(int rec=0; rec<answersize; ++rec)
    {
      std::vector<char> rdata;
      int length = 0;
      int tag = -1;
      int from = -1;
      exporter.ReceiveAny(from,tag,rdata,length);
      if (tag != 1234)
        dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

      // store received data
      std::vector<char>::size_type position = 0;
      while (position < rdata.size())
      {
        // ---- extract received data -----
        BindEventData reccldata;
        UnPack(position,rdata,reccldata);

        // add binding events to new colbond map
        if(reccldata.permission)
          myelebonds[reccldata.clgid] = reccldata;

      } // loop over answers by specific proc

      if (position != rdata.size())
        dserror("Mismatch in size of data %d <-> %d",(int)rdata.size(),position);

    } // loop over all answering procs

    // -----------------------------------------------------------------------
    // wait for all communication to finish
    // -----------------------------------------------------------------------
    for (int i=0; i<length; ++i)
      exporter.Wait(request[i]);

    // if we have done everything correct, this should be a no time operation
    bindis_->Comm().Barrier(); // I feel better this way ;-)

  } // communicate decisions, receive and colbond map

  // that is it
  return;

} // MakeBindingDecision()

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateDoubleBondedCl(
  std::map<int, DoubleBondedCl>&               mynewdbondcl,
  std::map<int, std::vector<DoubleBondedCl> >& dbondcltosend)
{
  CheckInit();

  // time needed to make binding decision in parallel
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::CommunicateDoubleBondedCl");

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  DRT::Exporter exporter(bindis_->Comm());
  const int length = dbondcltosend.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  std::map<int, std::vector<DoubleBondedCl> >::iterator p;
  for(p=dbondcltosend.begin(); p!=dbondcltosend.end(); ++p)
  {
    // ---- pack data for sending -----
    std::vector<char> sdata;
    std::vector<DoubleBondedCl>::const_iterator dbcl;
    DRT::PackBuffer data;
    for(dbcl=p->second.begin(); dbcl!=p->second.end(); ++dbcl)
    {
      Pack(data,*dbcl);
    }
    data.StartPacking();
    for(dbcl=p->second.begin(); dbcl!=p->second.end(); ++dbcl)
    {
      Pack(data,*dbcl);
    }
    swap(sdata,data());

    // unblocking send
    exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  // get number of procs from which myrank receives data
  std::vector<int> targetprocs(numproc_,0);
  std::map<int, std::vector<DoubleBondedCl> >::iterator prociter;
  for(prociter=dbondcltosend.begin(); prociter!=dbondcltosend.end(); ++prociter)
    targetprocs[prociter->first] = 1;
  // store number of messages myrank receives
  std::vector<int> summedtargets(numproc_,0);
  bindis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc_);

  // myrank receive all packs that are sent to him
  for(int rec=0; rec<summedtargets[myrank_]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      // ---- extract received data -----
      DoubleBondedCl recdata;
      UnPack(position,rdata,recdata);

      // add received data to list of now double bonded cl on myrank
      mynewdbondcl[recdata.id] = recdata;

    } // loop over packs from one proc

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d",(int)rdata.size(),position);

  } // loop over all procs from which myrank receives data

  // -----------------------------------------------------------------------
  // wait for all communication to finish
  // -----------------------------------------------------------------------
  for (int i=0; i<length; ++i)
    exporter.Wait(request[i]);

  // note: if we have done everything correct, this should be a no time operation
  bindis_->Comm().Barrier(); // I feel better this way ;-)

  // that is it
  return;

} // CommunicateDoubleBondedCl()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateCrosslinkerUnbinding(
  std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
  std::vector<UnBindEventData>&                 myrankunbindevent)
{
  CheckInit();

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  DRT::Exporter exporter(bindis_->Comm());
  const int length = sendunbindevent.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  std::map<int, std::vector<UnBindEventData> >::iterator p;
  for(p=sendunbindevent.begin(); p!=sendunbindevent.end(); ++p)
  {
    // ---- pack data for sending -----
    std::vector<char> sdata;
    std::vector<UnBindEventData>::const_iterator ube;
    DRT::PackBuffer data;
    for(ube=p->second.begin(); ube!=p->second.end(); ++ube)
    {
      Pack(data,*ube);
    }
    data.StartPacking();
    for(ube=p->second.begin(); ube!=p->second.end(); ++ube)
    {
      Pack(data,*ube);
    }
    swap(sdata,data());

    // unblocking send
    exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  // get number of procs from which myrank receives data
  std::vector<int> targetprocs(numproc_,0);
  std::map<int, std::vector<UnBindEventData> >::iterator prociter;
  for(prociter=sendunbindevent.begin(); prociter!=sendunbindevent.end(); ++prociter)
    targetprocs[prociter->first] = 1;
  // store number of messages myrank receives
  std::vector<int> summedtargets(numproc_,0);
  bindis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc_);

  // myrank receive all packs that are sent to him
  for(int rec=0; rec<summedtargets[myrank_]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      // ---- extract received data -----
      UnBindEventData recdata;
      UnPack(position,rdata,recdata);

      // add received data to list of unbindevents on myrank
      myrankunbindevent.push_back(recdata);

    } // loop over packs from one proc

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d",(int)rdata.size(),position);

  } // loop over all procs from which myrank receives data

  // -----------------------------------------------------------------------
  // wait for all communication to finish
  // -----------------------------------------------------------------------
  for (int i=0; i<length; ++i)
    exporter.Wait(request[i]);

  // note: if we have done everything correct, this should be a no time operation
  bindis_->Comm().Barrier(); // I feel better this way ;-)

  // that is it
  return;

} // CommunicateCrosslinkerUnbinding()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetNeighboringEles()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::TimeForBeamToBeamInteraction");

#ifdef DEBUG
  if((int)exteletobinmap_.size() != intactdis_->ElementColMap()->NumMyElements())
    dserror("std::map does not equal elecolmap (check e.g. if extended ghosting contains "
            " standard ghosting). Therefore not every contact can be detected");
#endif

  // loop over all column elements
  std::map<int, std::set<int> >::const_iterator coleleiter;
  for(coleleiter=exteletobinmap_.begin(); coleleiter!=exteletobinmap_.end(); ++coleleiter)
  {
    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all row bins
    std::set<int>::const_iterator colbiniter;
    for(colbiniter=coleleiter->second.begin(); colbiniter!=coleleiter->second.end(); ++colbiniter)
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      particlealgo_->GetNeighbouringBinGids(*colbiniter,loc_neighboring_binIds);

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert(loc_neighboring_binIds.begin(), loc_neighboring_binIds.end());
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds(neighboring_binIds.begin(), neighboring_binIds.end());

    // set of elements that lie in neighboring bins
    std::set<DRT::Element*> neighboring_elements;

    particlealgo_->GetBinContent(neighboring_elements,bin_beamcontent_,glob_neighboring_binIds);

    // evaluate contact and potential
//      Evaluate();

  }
  // that is it
  return;

} //GetneighboringEles()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateMaps()
{
  CheckInit();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export. See in particle_timint
  // how this can be done.

  // beam displacement
  UpdateDofMapOfVector(intactdis_, ia_disnp_);

  // update stiff matrix
  // todo:
  // COUPLING ADAPTER, drt_fsi/fsi??????????????

  // that is it
  return;

} //TransformStates()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateDofMapOfVector(
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

  // that is it
  return;

} //UpdateVectorDofMap()

///*-----------------------------------------------------------------------------*
// *-----------------------------------------------------------------------------*/
//template<typename T>
//void STR::MODELEVALUATOR::Crosslinking::ISendRecv(
//  std::map<int, std::vector<T> >& send,
//  std::map<int, T>&               recv)
//{
//  CheckInit();
//
//   // time needed to make binding decision in parallel
//   TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::CommunicateDoubleBondedCl");
//
//   // -----------------------------------------------------------------------
//   // send
//   // -----------------------------------------------------------------------
//   DRT::Exporter exporter(bindis_->Comm());
//   const int length = send.size();
//   std::vector<MPI_Request> request(length);
//   int tag = 0;
//   typename std::map<int, std::vector<T> >::iterator p;
//   for(p=send.begin(); p!=send.end(); ++p)
//   {
//     // ---- pack data for sending -----
//     std::vector<char> sdata;
//     typename std::vector<T>::const_iterator iter;
//     DRT::PackBuffer data;
//     for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
//     {
//       Pack(data,*iter);
//     }
//     data.StartPacking();
//     for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
//     {
//       PackStruct(data,*iter);
//     }
//     swap(sdata,data());
//
//     // unblocking send
//     exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
//     ++tag;
//   }
//   if (tag != length) dserror("Number of messages is mixed up");
//
//   // -----------------------------------------------------------------------
//   // receive
//   // -----------------------------------------------------------------------
//   // ---- prepare receiving procs -----
//   // get number of procs from which myrank receives data
//   std::vector<int> targetprocs(numproc_,0);
//   typename std::map<int, std::vector<T> >::iterator prociter;
//   for(prociter=send.begin(); prociter!=send.end(); ++prociter)
//     targetprocs[prociter->first] = 1;
//   // store number of messages myrank receives
//   std::vector<int> summedtargets(numproc_,0);
//   bindis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc_);
//
//   // myrank receive all packs that are sent to him
//   for(int rec=0; rec<summedtargets[myrank_]; ++rec)
//   {
//     std::vector<char> rdata;
//     int length = 0;
//     int tag = -1;
//     int from = -1;
//     exporter.ReceiveAny(from,tag,rdata,length);
//     if (tag != 1234)
//       dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);
//
//     // store received data
//     std::vector<char>::size_type position = 0;
//     while (position < rdata.size())
//     {
//       // ---- extract received data -----
//       T recdata;
//       UnPack(position,rdata,recdata);
//
//       // add received data to list of now double bonded cl on myrank
////       mydbondcl[recdata.id] = recdata;
//
//     } // loop over packs from one proc
//
//     if (position != rdata.size())
//       dserror("Mismatch in size of data %d <-> %d",(int)rdata.size(),position);
//
//   } // loop over all procs from which myrank receives data
//
//   // -----------------------------------------------------------------------
//   // wait for all communication to finish
//   // -----------------------------------------------------------------------
//   for (int i=0; i<length; ++i)
//     exporter.Wait(request[i]);
//
//   // note: if we have done everything correct, this should be a no time operation
//   bindis_->Comm().Barrier(); // I feel better this way ;-)
//
//   // that is it
//   return;
//}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&     data,
  const BindEventData& bindeventdata)
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,bindeventdata.clgid);
  DRT::ParObject::AddtoPack(data,bindeventdata.elegid);
  DRT::ParObject::AddtoPack(data,bindeventdata.bspotlocn);
  DRT::ParObject::AddtoPack(data,bindeventdata.requestproc);
  DRT::ParObject::AddtoPack(data,bindeventdata.permission);

  // that is it
  return;

} //PackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&       data,
  const UnBindEventData& unbindeventdata)
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,unbindeventdata.id);
  DRT::ParObject::AddtoPack(data,unbindeventdata.eletoupdate);

  // that is it
  return;

} //PackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&      data,
  const DoubleBondedCl& dbondedcldata)
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,dbondedcldata.id);
  DRT::ParObject::AddtoPack(data,dbondedcldata.eleone);
  DRT::ParObject::AddtoPack(data,dbondedcldata.eletwo);

  // that is it
  return;

} //PackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  BindEventData&                bindeventdata)
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.clgid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.elegid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.bspotlocn);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.requestproc);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.permission);

  // that is it
  return;

} //UnPackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  DoubleBondedCl&               dbondedcldata)
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,dbondedcldata.id);
  DRT::ParObject::ExtractfromPack(position,data,dbondedcldata.eleone);
  DRT::ParObject::ExtractfromPack(position,data,dbondedcldata.eletwo);

  // that is it
  return;

} //UnPackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  UnBindEventData&              unbindeventdata)
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,unbindeventdata.id);
  DRT::ParObject::ExtractfromPack(position,data,unbindeventdata.eletoupdate);

  // that is it
  return;

} //UnPackStruct()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Logo()
{
  CheckInit();

//  if(myrank_==0)
//  {
//    IO::cout << "\n****************************************************************" << IO::endl;
//    IO::cout << "*                                                              *" << IO::endl;
//    IO::cout << "*          Welcome to a Biopolymer Network Simulation          *" << IO::endl;
//    IO::cout << "*                                                              *" << IO::endl;
//    IO::cout << "****************************************************************" << IO::endl;
//    IO::cout << "                                                                  " << IO::endl;
//    IO::cout << "                                                                  " << IO::endl;
//    IO::cout << "                      0=========================0                 " << IO::endl;
//    IO::cout << "                    //|   \\            /       /||                " << IO::endl;
//    IO::cout << "                   // |    \\ |       |/       //||                " << IO::endl;
//    IO::cout << "                  //  |  /  \\|       /       // ||                " << IO::endl;
//    IO::cout << "                 //   |  \\   \\   /  /|\\     //  ||                " << IO::endl;
//    IO::cout << "                //    |  /   |\\ /  / | \\   //   ||                " << IO::endl;
//    IO::cout << "               //     |  \\   | \\     |  \\ //  / ||                " << IO::endl;
//    IO::cout << "              //  \\  /|  /   |/      |   //  /  ||                " << IO::endl;
//    IO::cout << "              0=========================0 \\ /   ||                " << IO::endl;
//    IO::cout << "             ||    /\\ |____          |  || \\    ||                " << IO::endl;
//    IO::cout << "             ||   /  \\|    \\   ------   ||/ \\   ||                " << IO::endl;
//    IO::cout << "             ||  /    |                 ||      ||                " << IO::endl;
//    IO::cout << "             || /     0----------/------||------0-                " << IO::endl;
//    IO::cout << "             ||      /   /       \\      ||     //                 " << IO::endl;
//    IO::cout << "             ||     /___/  \\     /    / ||    //                  " << IO::endl;
//    IO::cout << "             ||    /        \\    \\   /  ||   //                   " << IO::endl;
//    IO::cout << "             ||   /  \\/\\/\\/  \\   /  /   ||  //                    " << IO::endl;
//    IO::cout << "             ||  /      /     \\  \\ /    || //                     " << IO::endl;
//    IO::cout << "             || /      /         /      ||//                      " << IO::endl;
//    IO::cout << "             ||/                       /||/                       " << IO::endl;
//    IO::cout << "              0=========================0                         " << IO::endl;
//    IO::cout << "                                                                     " << IO::endl;
//    IO::cout << "                                                                     " << IO::endl;
//  }

  // that is it
  return;

} //Logo()

///*----------------------------------------------------------------------------*
// *----------------------------------------------------------------------------*/
//void STR::MODELEVALUATOR::Crosslinking::CommunicateElementCreation(
//  std::map<int, std::set<Teuchos::RCP<DRT::Element>, STR::MODELEVALUATOR::Less > > elestosend,
//  std::set<Teuchos::RCP<DRT::Element>, STR::MODELEVALUATOR::Less >&                addelesonmyrank)
//{
//  CheckInit();
//
//  // time needed to communicate new crosslinker elements
//  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::CommunicateNewElements");
//
//  // -----------------------------------------------------------------------
//  // send
//  // -----------------------------------------------------------------------
//  DRT::Exporter exporter(bindis_->Comm());
//  const int length = elestosend.size();
//  std::vector<MPI_Request> request(length);
//  int tag = 0;
//  std::map<int, std::set<Teuchos::RCP<DRT::Element>, STR::MODELEVALUATOR::Less > >::const_iterator p;
//  for(p=elestosend.begin(); p!=elestosend.end(); ++p)
//  {
//    // ---- pack data for sending -----
//    std::vector<char> sdata;
//    std::set<Teuchos::RCP<DRT::Element>, STR::MODELEVALUATOR::Less >::const_iterator eleiter;
//    for(eleiter=p->second.begin(); eleiter!=p->second.end(); ++eleiter)
//    {
//      DRT::PackBuffer data;
//      Teuchos::RCP<DRT::Element> e = *eleiter;
//      e->Pack(data);
//      data.StartPacking();
//      e->Pack(data);
//      sdata.insert(sdata.end(),data().begin(),data().end());
//    }
//    // unblocking send
//    exporter.ISend(myrank_, p->first, &(sdata[0]), (int)(sdata.size()), 1234, request[tag]);
//    ++tag;
//  }
//  if (tag != length) dserror("Number of messages is mixed up");
//
//  // -----------------------------------------------------------------------
//  // receive
//  // -----------------------------------------------------------------------
//  // ---- prepare receiving procs -----
//
//  // ---- prepare receiving procs -----
//  std::vector<int> targetprocs(numproc_,0);
//  std::map<int, std::set<Teuchos::RCP<DRT::Element>, STR::MODELEVALUATOR::Less > >::const_iterator piter;
//  for(piter=elestosend.begin(); piter!=elestosend.end(); ++piter)
//    targetprocs[piter->first] = 1;
//  std::vector<int> summedtargets(numproc_,0);
//  DiscretPtr()->Comm().SumAll(targetprocs.data(),summedtargets.data(), numproc_);
//
//  // -----------------------------------------------------------------------
//  // receive
//  // -----------------------------------------------------------------------
//  for(int rec=0; rec<summedtargets[myrank_]; ++rec)
//  {
//    std::vector<char> rdata;
//    int length = 0;
//    int tag = -1;
//    int from = -1;
//    exporter.ReceiveAny(from,tag,rdata,length);
//    if (tag != 1234)
//      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);
//
//    // put add received element to discretization later on
//    std::vector<char>::size_type index = 0;
//    while (index < rdata.size())
//    {
//      std::vector<char> data;
//      DRT::ParObject::ExtractfromPack(index,rdata,data);
//      // this Teuchos::rcp holds the memory of the node
//      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);
//      Teuchos::RCP<DRT::Element> ele = Teuchos::rcp_dynamic_cast<DRT::Element>(object);
//
//      // add element to elements that will be added later on
//      addelesonmyrank.insert(ele);
//
//    } // loop over all eles of one proc myrank needs to add
//  } // loop over all procs myrank gets eles from
//
//  // -----------------------------------------------------------------------
//  // wait for all communication to finish
//  // -----------------------------------------------------------------------
//  for (int i=0; i<length; ++i)
//    exporter.Wait(request[i]);
//
//  // if we have done everything correct, this should be a no time operation
//  bindis_->Comm().Barrier(); // I feel better this way ;-)
//
//  // that is it
//  return;
//
//} //CommunicateElementCreation()

///*-----------------------------------------------------------------------------*
// | nodes are checked and transferred if necessary              eichinger 09/16 |
// *-----------------------------------------------------------------------------*/
//void PARTICLE::Algorithm::TransferNodes(
//    Teuchos::RCP<DRT::Discretization> discret,
//    Teuchos::RCP<Epetra_Vector>       disnp)
//{
//  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferNodes");
//
//  // set of homeless nodes
//  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less> homelessnodes;
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
