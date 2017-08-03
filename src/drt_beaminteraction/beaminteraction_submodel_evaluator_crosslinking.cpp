/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_crosslinking.cpp

\brief class for submodel crosslinking

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#include "beaminteraction_submodel_evaluator_crosslinking.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/runtime_vtp_writer.H"
#include "../drt_io/discretization_runtime_vtp_writer.H"
#include "../drt_structure_new/str_timint_basedataio_runtime_vtp_output.H"
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_new/str_timint_basedataio.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beaminteraction/drt_utils_parallel_proctoproc.H"
#include "../drt_geometry/intersection_math.H"

#include "../drt_beaminteraction/crosslinker_node.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_particle/particle_handler.H"
#include "../drt_beaminteraction/crosslinking_params.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"

#include "../drt_mat/crosslinkermat.H"
#include "../drt_beam3/beam3_base.H"

#include "beam_link.H"
#include "beam_link_beam3r_lin2_rigidjointed.H"
#include "beam_link_beam3r_lin2_pinjointed.H"
#include "../drt_binstrategy/drt_meshfree_multibin.H"

#include "beaminteraction_calc_utils.H"

#include <unordered_set>


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Crosslinking() :
    crosslinking_params_ptr_( Teuchos::null ),
    vtp_writer_ptr_( Teuchos::null )
{
  crosslinker_data_.clear();
  beam_data_.clear();
  doublebondcl_.clear();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Setup()
{
  CheckInit();

  // construct, init and setup data container for crosslinking
  crosslinking_params_ptr_ = Teuchos::rcp( new BEAMINTERACTION::CrosslinkingParams() );
  crosslinking_params_ptr_->Init( GState() );
  crosslinking_params_ptr_->Setup();

  // set binding spot positions on filament elements according input file specifications
  SetFilamentTypes();
  SetFilamentBindingSpotPositions();

  // add crosslinker to bin discretization
  AddCrosslinkerToBinDiscretization();

  // build runtime vtp writer
  if ( GInOutput().GetRuntimeVtpOutputParams() != Teuchos::null )
    InitOutputRuntimeVtpStructure();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostPartitionProblem()
{
  CheckInitSetup();
  bool repartition = DRT::Problem::Instance()->Restart() ? false : SetAllPossibleInitialDoubleBondedCrosslinker();

  // we need to partition our problem again
  return repartition;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostSetup()
{
  CheckInitSetup();

  // in case of initially set crosslinker, we
  if ( crosslinking_params_ptr_->TotalNumInitCrosslinker() > 0 and not DRT::Problem::Instance()->Restart() )
    UpdateMyDoubleBondsRemoteIdList();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::InitSubmodelDependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap )
{
  CheckInitSetup();
  // no active influence on other submodels
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetFilamentTypes()
{
  CheckInit();

  std::set<int> examined_fils;

  // loop over all col nodes
  for ( int coln = 0; coln < Discret().NumMyColNodes(); ++coln )
  {
    DRT::Node* currnode = Discret().lColNode(coln);
    // get filament number of current node ( requirement: node belongs to only one filament)
    DRT::Condition* cond = currnode->GetCondition("BeamLineFilamentCondition");

    // in case node (e.g. node of rigid sphere element) does not belong to a filament, go to next node
    if ( cond == NULL )
      continue;

    // get filament type
    INPAR::BEAMINTERACTION::FilamentType filtype =
        INPAR::BEAMINTERACTION::String2FilamentType( *(cond->Get<std::string>("Type") ) );

    for ( int i = 0; i < currnode->NumElement(); ++i )
    {
      DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>( currnode->Elements()[i] );

#ifdef DEBUG
      if ( beamele == NULL )
        dserror(" DESIGN LINE BEAM FILAMENT CONDITIONS only suitable for beam elements.");
#endif

      beamele->SetFilamentType(filtype);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetFilamentBindingSpotPositions()
{
  CheckInit();

  // get pointers to all filament number conditions set
  std::vector<DRT::Condition*> filamentconditions(0);
  DiscretPtr()->GetCondition( "BeamLineFilamentCondition", filamentconditions );

  // get distance between two binding spots on a filament
  double const filamentbspotinterval = crosslinking_params_ptr_->FilamentBspotInterval();
  // get arc parameter range on which binding spot positions are set on filament
  std::pair<double, double> const filamentbspotrange =
      std::make_pair( ( crosslinking_params_ptr_->FilamentBspotRange()[0] < 0.0 ) ? 0.0 : crosslinking_params_ptr_->FilamentBspotRange()[0],
          crosslinking_params_ptr_->FilamentBspotRange()[1] );
  // todo: set somewhere else
  double const tol = GEO::TOL7;

  // temporarily extend ghosting, should be reverted somewhere later on as it is needed
  // in this fashion only here
  std::set<int> relevantfilaments;
  ExtendGhostingForFilamentBspotSetup( relevantfilaments );

  // loop over all relevant (on myrank) filaments
  std::set<int>::const_iterator filiter;
  for ( filiter = relevantfilaments.begin(); filiter != relevantfilaments.end(); ++filiter )
  {
    // loop over all nodes of current filament, sort elements and calculate total filament length
    std::vector<int> const* nodeids = filamentconditions[*filiter]->Nodes();
    std::vector< DRT::Element* > sortedfilamenteles(0);
    double filreflength = 0.0;
    ComputeFilamentLengthAndSortItsElements( sortedfilamenteles, nodeids, filreflength );

    // get arc parameter range for binding spot positions for current filament
    double const start = ( filreflength < filamentbspotrange.first ) ?
        ( filreflength + 1.0 ) : (filamentbspotrange.first);
    double const end = ( filreflength < filamentbspotrange.second or filamentbspotrange.second < 0.0 ) ?
        filreflength : filamentbspotrange.second;
    // get number of binding spots for current filament in bonds
    int numbspot = ( start < filreflength and ( (end - start) > 0.0 ) ) ?
        ( std::floor( (end + tol - start) / filamentbspotinterval ) + 1 ) : 0;

    // in case current filament has no binding spots, we go to next filament
    if ( numbspot == 0 ) continue;

    // print number of binding spots for current filament
    IO::cout(IO::verbose) << "\n---------------------------------------------------------------"<< IO::endl;
    IO::cout(IO::verbose) << numbspot << " binding spots on filament " << *filiter << " (consists of " <<
        static_cast<int>(sortedfilamenteles.size()) << " elements)" << " with: " << IO::endl;

    // set xis on element level
    SetBindingSpotsPositionsOnFilament( sortedfilamenteles, start, numbspot,
        filamentbspotinterval, tol );

    IO::cout(IO::verbose) << "---------------------------------------------------------------\n"<< IO::endl;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ExtendGhostingForFilamentBspotSetup(
    std::set<int>& relevantfilaments)
{
  CheckInit();

  std::set<int> setofnodegidswithrequiredelecloud;
  DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
      relevantfilaments, setofnodegidswithrequiredelecloud );

  // collect element gids that need to be ghosted so that each proc ghosts all
  // elements of a filament containing at least one node of myrank
  // do communication to gather all elements to temporarily extend ghosting
  std::set<int> coleleset;
  for ( int iproc = 0; iproc < DiscretPtr()->Comm().NumProc(); ++iproc )
  {
    // myrank == iproc: copy set to vector in order to broadcast data
    std::vector<int> requirednodes(0);
    if ( iproc == GState().GetMyRank() )
      requirednodes.insert( requirednodes.begin(), setofnodegidswithrequiredelecloud.begin(),
          setofnodegidswithrequiredelecloud.end() );

    // proc i tells all procs how many nodegids it has
    int numnodes = requirednodes.size();
    DiscretPtr()->Comm().Broadcast( &numnodes, 1, iproc );

    // proc i actually sends nodegids
    requirednodes.resize(numnodes);
    DiscretPtr()->Comm().Broadcast( &requirednodes[0], numnodes, iproc );

    std::set<int> sdata;
    std::set<int> rdata;

    // each proc looks in node row map and adds element clouds of owned nodes to sdata
    for ( int i = 0; i < numnodes; ++i )
    {
      // only if myrank is owner of requested node
      if ( DiscretPtr()->NodeRowMap()->LID(requirednodes[i]) < 0 )
        continue;

      // insert element cloud of current node
      DRT::Node* currnode = DiscretPtr()->gNode(requirednodes[i]);
      for ( int j = 0; j < currnode->NumElement(); ++j )
        sdata.insert( currnode->Elements()[j]->Id() );
    }

    // gather and store information on iproc
    LINALG::Gather<int>( sdata, rdata, 1, &iproc, DiscretPtr()->Comm() );
    if ( iproc == GState().GetMyRank() )
      coleleset = rdata;
  }

  // insert previous ghosting
  for ( int lid = 0; lid < DiscretPtr()->NumMyColElements(); ++lid )
    coleleset.insert( DiscretPtr()->ElementColMap()->GID(lid) );
  std::vector<int> colgids( coleleset.begin(), coleleset.end() );

  // create new ele col map
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp( new Epetra_Map( -1,
      static_cast<int>( colgids.size() ), &colgids[0], 0, DiscretPtr()->Comm() ) );

  // temporarily extend ghosting
  BINSTRATEGY::UTILS::ExtendDiscretizationGhosting( DiscretPtr(), newelecolmap, true , false , true );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
        std::set<int>& relevantfilaments,
        std::set<int>& setofrequirednodes) const
{
  CheckInit();

  // loop over all row nodes
  for ( int rown = 0; rown < Discret().NumMyRowNodes(); ++rown )
  {
    // get filament number of current node ( requirement: node belongs to only one filament)
    DRT::Condition* cond = Discret().lRowNode(rown)->GetCondition("BeamLineFilamentCondition");

    // in case node (e.g. node of rigid sphere element) does not belong to a filament, go to next node
    if ( cond == NULL )
      continue;

    // get filament number
    int const currfilnum = cond->GetInt("FilamentId");

    // if a filament has already been examined --> continue with next node
    if ( relevantfilaments.find(currfilnum) != relevantfilaments.end() )
      continue;
    // filament is examined for the first time --> new entry for relevant filaments on myrank
    else
      relevantfilaments.insert(currfilnum);

    // loop over all nodes of current filament and store gids of nodes not owned by myrank
    std::vector<int> const* nodeids = cond->Nodes();
    for ( int i = 0; i < static_cast<int>( nodeids->size() ) ; ++i )
      if ( Discret().NodeRowMap()->LID((*nodeids)[i]) < 0 )
        setofrequirednodes.insert( (*nodeids)[i] );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ComputeFilamentLengthAndSortItsElements(
    std::vector< DRT::Element* >& sortedfilamenteles,
    std::vector<int> const* nodeids,
    double & filreflength) const
{
  CheckInit();

  // loop over all nodes associated with current filament
  for ( int nodei = 0; nodei < static_cast<int>( nodeids->size() ); ++nodei )
  {
    // insert element cloud of current node
    DRT::Node* node = Discret().gNode( (*nodeids)[nodei] );
    for ( int j = 0; j < node->NumElement(); ++j )
    {
      // only if element has not yet been added to the filaments elements
      if ( std::find( sortedfilamenteles.begin(), sortedfilamenteles.end(), node->Elements()[j] )
                != sortedfilamenteles.end() )
        continue;

#ifdef DEBUG
      if ( dynamic_cast<DRT::ELEMENTS::Beam3Base*>(node->Elements()[j]) == NULL )
        dserror("DESIGN LINE BEAM FILAMENT CONDITIONS only applicable to beam elements.");
#endif

      // add element reference length of new element to filament reference length
      filreflength += dynamic_cast<DRT::ELEMENTS::Beam3Base*>(node->Elements()[j])->RefLength();
      // add element
      sortedfilamenteles.push_back( node->Elements()[j] );
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetBindingSpotsPositionsOnFilament(
    std::vector< DRT::Element* >& sortedfilamenteles,
    double const start,
    int const    numbspot,
    double const filamentbspotinterval,
    double const tol)
{
  CheckInit();

  // default no occupied binding spots
  std::pair<double, double > elearcinterval = std::make_pair( 0.0, 0.0 );
  int bspotcounter = 0;
  double currbspotarcparam = start;
  // loop over all elements in sorted order
  for ( int se_iter = 0; se_iter < static_cast<int>( sortedfilamenteles.size() ); ++se_iter )
  {
    DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(sortedfilamenteles[se_iter]);

    // init variables set in beam eles
    std::map<int, int> bspotstatus;
    std::vector<double> bspotposxi;

    // get arc parameter range of current element
    elearcinterval = std::make_pair( elearcinterval.second, elearcinterval.second + beamele->RefLength() );

    while ( ( currbspotarcparam > ( elearcinterval.first - tol) and currbspotarcparam <= ( elearcinterval.second + tol ) )
        and bspotcounter < numbspot )
    {
      // default unoccupied binding spots
      bspotstatus[bspotposxi.size()] = -1;

      // linear mapping from arc pos to xi (local element system)
      double xi = ( 2.0 * currbspotarcparam - elearcinterval.first - elearcinterval.second ) / beamele->RefLength();
      xi = ( abs( round(xi) - xi ) < tol ) ? round(xi) : xi;
      bspotposxi.push_back(xi);

      // print to screen
      IO::cout(IO::verbose) << bspotcounter + 1 << ". binding spot: "<< "xi = " << xi << " on element "
          << se_iter << " (gid " << beamele->Id() << ")"<< IO::endl;

      ++bspotcounter;
      currbspotarcparam += filamentbspotinterval;
    }

    // set element variables for binding spot position and status
    beamele->SetBindingSpotPositionsAndStatus( bspotposxi, bspotstatus );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetAllPossibleInitialDoubleBondedCrosslinker()
{
  CheckInit();

  // in case no initial crosslinker are set, return without repartitioning the problem
  if ( crosslinking_params_ptr_->TotalNumInitCrosslinker() == 0 )
    return false;

  // id is unique pairing id of elegid and locbspot id of binding spot
  // built up two std::vectors instead of one std::vector< std::pair< int, int > > to ease
  // communication later on
  std::vector<int> myfirstbspotids, mysecondbspotids, mylinkertype;
  std::vector<BeamData> beam_data;
  PreComputeBeamData(beam_data);

  // loop over all row beam elements
  int unsigned const numbeams = EleTypeMapExtractorPtr()->BeamMap()->NumMyElements();
  for ( unsigned int rowbeam_i = 0; rowbeam_i < numbeams; ++rowbeam_i )
  {
    DRT::ELEMENTS::Beam3Base * beamele = dynamic_cast< DRT::ELEMENTS::Beam3Base* >(
        Discret().gElement( EleTypeMapExtractorPtr()->BeamMap()->GID(rowbeam_i) ) );

    BeamData const & beamdata_i = beam_data[ beamele->LID() ];

    // loop over all binding spots of current element
    unsigned int numbspots = beamele->GetNumberOfBindingSpots();
    for ( unsigned int locbspot_i = 0; locbspot_i < numbspots; ++locbspot_i )
    {
      // id uniquely identifying elegid and loc bspot number
      int const bspotgid = BEAMINTERACTION::UTILS::CantorPairing(
          std::make_pair( beamele->Id(), locbspot_i ), true );

      // get bin of current bspot
      int const bingid = BinStrategy().ConvertPosToGid( beamdata_i.bbspotpos.at(locbspot_i) );

      // get neighboring bins
      // note: interaction distance cl to beam needs to be smaller than the bin size
      std::vector<int> neighboring_binIds;
      neighboring_binIds.reserve(27);
      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds( bingid, neighboring_binIds );

      // get set of neighboring beam elements (i.e. elements that somehow touch nb bins)
      // we also need col elements (flag = false) here (in contrast to "normal" crosslinking)
      std::set< DRT::Element * > neighboring_beams;
      std::vector< BINSTRATEGY::UTILS::BinContentType > bc( 1, BINSTRATEGY::UTILS::Beam );
      BinStrategyPtr()->GetBinContent( neighboring_beams, bc, neighboring_binIds, false );

      // in case there are no neighbors, go to next binding spot
      if ( neighboring_beams.empty() )
        continue;

      // loop over all neighboring beam elements
      for ( auto const & nb_ele_i : neighboring_beams )
      {
        DRT::ELEMENTS::Beam3Base * nb_beamele =
            dynamic_cast< DRT::ELEMENTS::Beam3Base* >( nb_ele_i );

        BeamData const & nb_beamdata_i = beam_data[ nb_beamele->LID() ];

        // ------------------------------------------------------------
        // exclude linking of touching elements
        // -------------------------------------------------------------
        if ( BEAMINTERACTION::UTILS::DoBeamElementsShareNodes( beamele, nb_beamele ) )
          continue;

        // loop over binding spots of neighboring element
        for ( unsigned int nb_locbspot_i = 0; nb_locbspot_i < nb_beamele->GetNumberOfBindingSpots(); ++nb_locbspot_i )
        {
          // unique id of neighboring binding spot
          int const nb_bspotgid = BEAMINTERACTION::UTILS::CantorPairing(
              std::make_pair( nb_beamele->Id(), nb_locbspot_i ), true );

          // loop over different linker types and check for feasibility of link
          std::vector<int> const & matcrosslinkerpertype = crosslinking_params_ptr_->MatCrosslinkerPerType();
          for ( unsigned int type_i = 0; type_i < static_cast<unsigned int>(matcrosslinkerpertype.size()); ++type_i )
          {
            // check criteria for feasible bond
            Teuchos::RCP<MAT::CrosslinkerMat> mat =
                Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>( MAT::Material::Factory(matcrosslinkerpertype[type_i]) );

            // --------------------------------------
            // distance
            // --------------------------------------
            // minimum and maximum distance at which a double-bond crosslink can be established
            double const linkdistmin = mat->LinkingLength() - mat->LinkingLengthTolerance();
            double const linkdistmax = mat->LinkingLength() + mat->LinkingLengthTolerance();

#ifdef DEBUG
            if ( linkdistmax > BinStrategy().CutoffRadius() )
              dserror("The allowed binding distance of your linker material is greater than the "
                  "cutoff radius, this can lead to missed binding events" );
#endif

            if ( BEAMINTERACTION::UTILS::IsDistanceOutOfRange( beamdata_i.bbspotpos.at(locbspot_i),
                nb_beamdata_i.bbspotpos.at(nb_locbspot_i), linkdistmin, linkdistmax ) )
              continue;

            // ------------------------------------------------------------
            // orientation of centerline tangent vectors at binding spots
            // -------------------------------------------------------------
            // a crosslink (double-bonded crosslinker) will only be established if the
            // enclosed angle is in the specified range
            double const linkanglemin = mat->LinkingAngle() - mat->LinkingAngleTolerance();
            double const linkanglemax = mat->LinkingAngle() + mat->LinkingAngleTolerance();

            // get tangent of binding spot on beamele
            LINALG::Matrix< 3, 1 > bindingspot_beam_tangent(true);
            for ( unsigned int idim = 0; idim < 3; ++idim )
              bindingspot_beam_tangent(idim) = beamdata_i.bbspottriad.at(locbspot_i)( idim, 0 );

            // get tangent of binding spot on nb_beamele
            LINALG::Matrix< 3, 1 > nb_bindingspot_beam_tangent(true);
            for ( unsigned int idim = 0; idim < 3; ++idim )
              nb_bindingspot_beam_tangent(idim) = nb_beamdata_i.bbspottriad.at(nb_locbspot_i)( idim, 0 );

            if ( BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
                bindingspot_beam_tangent, nb_bindingspot_beam_tangent, linkanglemin, linkanglemax ) )
              continue;

            // if all criteria were met for a certain linker, add this bond to linker specific bond list
            myfirstbspotids.push_back(bspotgid);
            mysecondbspotids.push_back(nb_bspotgid);
            mylinkertype.push_back(matcrosslinkerpertype[type_i]);

            // go to next binding spot, if current one will poi
            break;
          }
        }
      }
    }
  }

  // gather all data over all procs
  Epetra_Comm const & com = Discret().Comm();
  const int numproc = com.NumProc();
  std::vector<int> firstbspotids, secondbspotids, linkertype;
  for ( int iproc = 0; iproc < numproc; ++iproc )
  {
    std::vector<int> firstbspotids_i, secondbspotids_i, linkertype_i;
    if ( iproc == GState().GetMyRank() )
    {
      firstbspotids_i = myfirstbspotids;
      secondbspotids_i = mysecondbspotids;
      linkertype_i = mylinkertype;
    }

    // first: proc i tells all procs how many pairs it has
    int numpairs = myfirstbspotids.size();
    com.Broadcast( &numpairs, 1, iproc );
    // second: proc i tells all procs which bins it has
    firstbspotids_i.resize(numpairs);
    secondbspotids_i.resize(numpairs);
    linkertype_i.resize(numpairs);

    com.Broadcast( &firstbspotids_i[0], numpairs, iproc );
    com.Broadcast( &secondbspotids_i[0], numpairs, iproc );
    com.Broadcast( &linkertype_i[0], numpairs, iproc );

    firstbspotids.insert( firstbspotids.end(), firstbspotids_i.begin(), firstbspotids_i.end() );
    secondbspotids.insert( secondbspotids.end(), secondbspotids_i.begin(), secondbspotids_i.end() );
    linkertype.insert( linkertype.end(), linkertype_i.begin(), linkertype_i.end() );
  }

  // safety check
  if ( firstbspotids.size() != secondbspotids.size() or firstbspotids.size() != linkertype.size() )
    dserror(" something went wrong during communication");

  // build map
  std::map< int, std::set< std::pair< int, int >, BEAMINTERACTION::UTILS::SetOfPairComparator > > bspotpairs;
  for ( unsigned int i = 0; i < static_cast<unsigned int>( firstbspotids.size() ); ++i )
    bspotpairs[firstbspotids[i]].insert( std::make_pair( secondbspotids[i], linkertype[i] ) );

  // store already occupied binding spots
  std::unordered_set<int> occupied;
  // vector of crosslinker added to bin discretization later on
  std::vector< CrosslinkerData > newlinker;
  // note: size can be bigger if more than one binding spot per element
  newlinker.reserve(Discret().ElementRowMap()->NumMyElements() );

  // loop over bspotpairs and make decision
  std::map< int, int > numpertype;
  int elegid = -1;
  int locbspotid = -1;
  int nb_elegid = -1;
  int nb_locbspotid = -1;
  for ( auto const & iter : bspotpairs )
  {
    int const bspotgid = iter.first;
    std::pair< int, int > bspotids = BEAMINTERACTION::UTILS::CantorDePairing(bspotgid);
    elegid = bspotids.first;
    locbspotid = bspotids.second;

    // check if binding spot has already been occupied
    if ( occupied.find(bspotgid) != occupied.end() )
      continue;
    else
      occupied.insert(bspotgid);

    // loop over neighboring bspot pairs
    for ( auto const & nbiter : iter.second )
    {
      int const nb_bspotgid = nbiter.first;
      std::pair< int, int > nb_bspotids = BEAMINTERACTION::UTILS::CantorDePairing(nb_bspotgid);
      nb_elegid = nb_bspotids.first;
      nb_locbspotid = nb_bspotids.second;

      // check if maximal number of crosslinker for current type is already exceeded
      ++numpertype[nbiter.second];
      if ( numpertype[nbiter.second] > crosslinking_params_ptr_->NumInitCrosslinkerPerCrosslinkerMatId(nbiter.second) )
        continue;

      // check if binding spot has already been occupied
      if ( occupied.find( nb_bspotgid ) != occupied.end() )
        continue;
      else
        occupied.insert( nb_bspotgid );

      // check ownership
      if ( Discret().ElementRowMap()->LID( elegid ) == -1 )
        break;

       // store data of new crosslinker
       CrosslinkerData cldata;

       // set positions
       SetPositionOfDoubleBondedCrosslinkerPBCconsistent( cldata.clpos, beam_data[DiscretPtr()->gElement(elegid)->LID()].bbspotpos.at(locbspotid),
           beam_data[DiscretPtr()->gElement(nb_elegid)->LID()].bbspotpos.at(nb_locbspotid) );
       // set current binding spot status of crosslinker
       cldata.clbspots.push_back( std::make_pair( elegid, locbspotid ) );
       cldata.clbspots.push_back( std::make_pair( nb_elegid, nb_locbspotid ) );
       // set number of bonds
       cldata.clnumbond = 2;
       // material pointer
       cldata.clmat = Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>( MAT::Material::Factory(nbiter.second) );
       // owner
       cldata.clowner = GState().GetMyRank();

       newlinker.push_back(cldata);

       // if current bspot has found a fitting neighbor, we go to next binding spot
       break;
    }
  }

  // determine unique gids on each proc (ascending over all procs)
  // gather numbers of new linker on each proc
  std::vector<int> nummynewlinks(1);
  nummynewlinks[0] = static_cast<int>( newlinker.size() );
  // initialize std::vector for communication
  std::vector<int> numnewlinks( com.NumProc(), 0 );
  // communicate
  com.GatherAll( &nummynewlinks[0], &numnewlinks[0], nummynewlinks.size() );
  com.Barrier();

  // calculate starting index on myrank
  int mystartgid = 0;
  for( int i = 0; i < GState().GetMyRank(); ++i )
    mystartgid += numnewlinks[i];

  // loop over new linker on myrank
  int gid = BinDiscretPtr()->NodeRowMap()->MaxAllGID() + 1 + mystartgid;
  std::vector< double > X(3);
  // map key is crosslinker gid to be able to uniquely address one entry over all procs
  std::map< int, NewDoubleBonds > mynewdbondcl;
  for( auto const & iter : newlinker )
  {
    for ( unsigned int dim = 0; dim < 3; ++dim )
      X[dim] = iter.clpos(dim);

    Teuchos::RCP<CROSSLINKING::CrosslinkerNode> newcrosslinker =
        Teuchos::rcp(new CROSSLINKING::CrosslinkerNode( gid, &X[0], GState().GetMyRank() ) );
    newcrosslinker->InitializeDataContainer();
    newcrosslinker->SetMaterial( iter.clmat );
    BinDiscretPtr()->AddNode(newcrosslinker);

    NewDoubleBonds dbondcl;
    dbondcl.id = gid;
    dbondcl.eleids.push_back( iter.clbspots[0] );
    dbondcl.eleids.push_back( iter.clbspots[1] );

    int colelelid = DiscretPtr()->ElementColMap()->LID( iter.clbspots[0].first );
    int nb_colelelid = DiscretPtr()->ElementColMap()->LID( iter.clbspots[1].first );
    dbondcl.bspotposs.push_back( beam_data[colelelid].bbspotpos.at( iter.clbspots[0].second) );
    dbondcl.bspotposs.push_back( beam_data[nb_colelelid].bbspotpos.at( iter.clbspots[1].second) );
    dbondcl.bspottriads.push_back( beam_data[colelelid].bbspottriad.at( iter.clbspots[0].second) );
    dbondcl.bspottriads.push_back( beam_data[nb_colelelid].bbspottriad.at( iter.clbspots[1].second) );

    mynewdbondcl[dbondcl.id] = dbondcl;

    // set correct states for linker
    newcrosslinker->ClData()->SetClBSpotStatus( iter.clbspots );
    newcrosslinker->ClData()->SetNumberOfBonds( iter.clnumbond );
    SetCrosslinkerPosition( newcrosslinker.get(), iter.clpos );

    DRT::ELEMENTS::Beam3Base* ele =
        dynamic_cast< DRT::ELEMENTS::Beam3Base * >( DiscretPtr()->gElement( iter.clbspots[0].first ) );
    ele->SetCertainBindingSpotStatus( std::make_pair( iter.clbspots[0].second, gid ) );

    if ( Discret().HaveGlobalElement( iter.clbspots[1].first ) )
    {
      DRT::ELEMENTS::Beam3Base* nb_ele =
          dynamic_cast< DRT::ELEMENTS::Beam3Base * >( DiscretPtr()->gElement( iter.clbspots[1].first ) );
      nb_ele->SetCertainBindingSpotStatus( std::make_pair( iter.clbspots[1].second, gid ) );
    }

    // update gid
    ++gid;
  }

  // setup new double bonds and insert them in doublebondcl_
  CreateNewDoubleBondedCrosslinkerElementPairs( mynewdbondcl );

  // set row map of newly created particle discretization
  BinDiscretPtr()->FillComplete( false, false, false );

  // flag indicating if a repartition needs to be done (true if initial crosslinker on some proc)
  return ( static_cast<int>( bspotpairs.size() ) > 0 ) ? true : false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::AddCrosslinkerToBinDiscretization()
{
  CheckInit();

  // loop over all linker types that should be added to simulation volume
  // (only proc 0 is doing this (as the number of crosslinker is manageable)
  std::vector<int> const & numcrosslinkerpertype = crosslinking_params_ptr_->NumCrosslinkerPerType();
  std::vector<int> const & matcrosslinkerpertype = crosslinking_params_ptr_->MatCrosslinkerPerType();
  int gid = 0;
  if ( GState().GetMyRank() == 0 )
  {
    for ( int cltype_i = 0; cltype_i < static_cast<int>( numcrosslinkerpertype.size() ); ++cltype_i )
    {
      for ( int cltype_i_cl_j = 0; cltype_i_cl_j < numcrosslinkerpertype[cltype_i]; ++cltype_i_cl_j )
      {
        // random reference position of crosslinker in bounding box
        LINALG::Matrix< 3, 1 > Xmat;
        PeriodicBoundingBoxPtr()->RandomPosWithin(Xmat);

        std::vector<double> X(3);
        for ( int dim = 0; dim < 3; ++dim )
          X[dim] = Xmat(dim);

        // construct node, init data container, set material and add to bin discret
        Teuchos::RCP<CROSSLINKING::CrosslinkerNode> newcrosslinker =
            Teuchos::rcp(new CROSSLINKING::CrosslinkerNode( gid++, &X[0], GState().GetMyRank() ) );
        newcrosslinker->InitializeDataContainer();
        newcrosslinker->SetMaterial( matcrosslinkerpertype[cltype_i] );
        BinDiscretPtr()->AddNode(newcrosslinker);
      }
    }
  }

  // set row map of newly created particle discretization
  BinDiscretPtr()->FillComplete( false, false, false );
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Reset()
{
  CheckInitSetup();

  // reset crosslinker pairs
  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

    CROSSLINKING::CrosslinkerNode *cl_i =
        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( BinDiscretPtr()->gNode( elepairptr->Id() ) );

    if ( cl_i->ClData()->GetNumberOfBonds() != 2 )
      dserror("Cl with gid %i Owner %i on myrank %i and numbonds %i",elepairptr->Id(),cl_i->Owner(),
          GStatePtr()->GetMyRank(), cl_i->ClData()->GetNumberOfBonds() );

    // init positions and triads
    std::vector<LINALG::Matrix< 3, 1 > > pos(2);
    std::vector<LINALG::Matrix< 3, 3 > > triad(2);

    for ( int i = 0; i < 2; ++i )
    {
      int elegid = elepairptr->GetEleGid(i);

      // safety check
      if ( DiscretPtr()->ElementColMap()->LID(elegid) < 0 )
      {
        elepairptr->Print(std::cout);
        dserror("Reset(): elegid %i not there on proc %i ", elegid, GState().GetMyRank() );
      }

      int locbspotnum = elepairptr->GetLocBSpotNum(i);
      DRT::Element* ele = DiscretPtr()->gElement(elegid);

      BEAMINTERACTION::UTILS::GetPosAndTriadOfBindingSpot( Discret(), ele,
          BeamInteractionDataStatePtr()->GetMutableDisColNp(),
          PeriodicBoundingBoxPtr(), locbspotnum, pos[i], triad[i] );
    }

    // unshift one of the positions if both are separated by a periodic boundary
    // condition, i.e. have been shifted before
    PeriodicBoundingBoxPtr()->UnShift3D( pos[1], pos[0] );

    // safety check until code is better tested for potential problems with periodic boundary conditions
    // **************************** DEBUG ****************************************
    LINALG::Matrix<3,1> dist ( true );
    dist.Update( 1.0, pos[0], -1.0, pos[1] );
    for ( unsigned int i = 0; i < 3; ++i )
    {
      if ( std::abs(dist(i)) > 0.5 * PeriodicBoundingBoxPtr()->EdgeLength(i) )
      {
        dserror("You are trying to set the binding spot positions of this crosslinker "
            "in at least one direction at a distance larger than %f, which is "
            " half of the period length in the respective direction)", 0.5 * PeriodicBoundingBoxPtr()->EdgeLength(i));
      }
    }
    // ********************** END DEBUG ****************************************

    // finally reset state
    elepairptr->ResetState( pos, triad );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateForce()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );

  // resulting discrete element force vectors of the two parent elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  std::vector< std::vector<LINALG::SerialDenseMatrix> > dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

    for ( int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();
    }

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForce( bspotforce[0], bspotforce[1] );

    // apply forces on binding spots to parent elements
    // and get their discrete element force vectors
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceToParentElements< DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base >(
        Discret(), PeriodicBoundingBoxPtr(), BeamInteractionDataStatePtr()->GetMutableDisColNp(),
        elepairptr, bspotforce, eleforce );

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
        eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateStiff()
{
  CheckInitSetup();

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector< std::vector<LINALG::SerialDenseMatrix> > bspotstiff( 2,
      std::vector<LINALG::SerialDenseMatrix>( 2, LINALG::SerialDenseMatrix(6,6) ) );

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector< std::vector<LINALG::SerialDenseMatrix> > elestiff( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  std::vector< LINALG::SerialDenseVector > dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

    for ( int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);

      for ( int j = 0; j< 2; ++j )
        bspotstiff[i][j].Zero();
    }

     // evaluate beam linkage object to get linearizations of forces and moments on binding spots
    elepairptr->EvaluateStiff( bspotstiff[0][0], bspotstiff[0][1], bspotstiff[1][0],
        bspotstiff[1][1] );

    // apply linearizations to parent elements and get their discrete element stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotStiffToParentElements< DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base >(
        Discret(),PeriodicBoundingBoxPtr(), BeamInteractionDataStatePtr()->GetMutableDisColNp(),
        elepairptr, bspotstiff, elestiff);

    // assemble the contributions into stiffness matrix class variable
    // stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(),
        elegids, dummyforce, elestiff, Teuchos::null,
        BeamInteractionDataStatePtr()->GetMutableStiff());
   }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateForceStiff()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector< std::vector<LINALG::SerialDenseMatrix> > bspotstiff( 2,
      std::vector<LINALG::SerialDenseMatrix>( 2, LINALG::SerialDenseMatrix(6,6) ) );

  // resulting discrete element force vectors of the two parent elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector< std::vector<LINALG::SerialDenseMatrix> > elestiff( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;
    for ( int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();

      for ( int j = 0; j< 2; ++j )
        bspotstiff[i][j].Zero();
    }

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForceStiff( bspotforce[0], bspotforce[1], bspotstiff[0][0],
        bspotstiff[0][1], bspotstiff[1][0], bspotstiff[1][1] );

    // apply forces on binding spots and corresponding linearizations to parent elements
    // and get their discrete element force vectors and stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceStiffToParentElements< DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base >(
        Discret(), PeriodicBoundingBoxPtr(), BeamInteractionDataStatePtr()->GetMutableDisColNp(),
        elepairptr, bspotforce, bspotstiff, eleforce, elestiff );

    // assemble the contributions into force and stiffness class variables
    // f_crosslink_np_ptr_, stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
        eleforce, elestiff, BeamInteractionDataStatePtr()->GetMutableForceNp(),
        BeamInteractionDataStatePtr()->GetMutableStiff() );
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreUpdateStepElement()
{
  CheckInitSetup();

  // crosslinker diffusion: - according to browninan dyn for free cl
  //                        - according to beams for single and double bonded
  DiffuseCrosslinker();

  // erase temporary data container as both distributions as well the col data
  // might have changed
  crosslinker_data_.clear();
  beam_data_.clear();

  // transfer crosslinker to new bins
  Teuchos::RCP<std::list<int> > lostcl = ParticleHandlerPtr()->TransferParticles();
  if ( not lostcl->empty() )
    dserror( "Crosslinker got lost during transfer, something went wrong");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateStepElement()
{
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::BindAndUnbindCrosslinker");

  CheckInitSetup();

  // -------------------------------------------------------------------------
  // update double bonded linker
  // -------------------------------------------------------------------------
  UpdateMyDoubleBondsAfterRedistribution();

  // -------------------------------------------------------------------------
  // gather data for all column crosslinker and column beams initially
  // -------------------------------------------------------------------------
  PreComputeCrosslinkerAndBeamData( crosslinker_data_, beam_data_ );

  // -------------------------------------------------------------------------
  // manage binding events
  // -------------------------------------------------------------------------
  BindCrosslinker();

  // -------------------------------------------------------------------------
  // manage unbinding events
  // -------------------------------------------------------------------------
  UnBindCrosslinker();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostUpdateStepElement()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();
  OutputStepState( false );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::OutputStepState(
    bool const restart ) const
{
  DRT::Discretization const& bindis = BinDiscret();

  // mesh is not written to disc, only maximum node id is important for output
  bindis.Writer()->ParticleOutput( GState().GetStepN(), GState().GetTimeN(), restart );
  bindis.Writer()->NewStep( GState().GetStepN(), GState().GetTimeN() );
  Teuchos::RCP<Epetra_Vector> dis         = LINALG::CreateVector( *bindis.DofRowMap(), true );
  Teuchos::RCP<Epetra_Vector> orientation = LINALG::CreateVector( *bindis.DofRowMap(), true );
  Teuchos::RCP<Epetra_Vector> numbond     = LINALG::CreateVector( *bindis.NodeRowMap(), true );
  Teuchos::RCP<Epetra_Vector> owner       = LINALG::CreateVector( *bindis.NodeRowMap(), true );

  FillStateDataVectorsForOutput( dis, orientation, numbond, owner);

  bindis.Writer()->WriteVector( "displacement", dis );
  bindis.Writer()->WriteVector( "orientation", orientation );
  bindis.Writer()->WriteVector( "numbond", numbond, IO::nodevector );
  bindis.Writer()->WriteVector( "owner", owner, IO::nodevector);

  // as we know that our maps have changed every time we write output, we can empty
  // the map cache as we can't get any advantage saving the maps anyway
  bindis.Writer()->ClearMapCache();

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RuntimeOutputStepState() const
{
  CheckInitSetup();

  if ( vtp_writer_ptr_ != Teuchos::null )
    WriteOutputRuntimeVtpStructure();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::InitOutputRuntimeVtpStructure()
{
  CheckInit();

  vtp_writer_ptr_ = Teuchos::rcp( new DiscretizationRuntimeVtpWriter() );

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  // initialize the writer object
  vtp_writer_ptr_->Initialize(
      BinDiscretPtr(),
      num_timesteps_in_simulation_upper_bound,
      GState().GetTimeN(),
      GInOutput().GetRuntimeVtpOutputParams()->WriteBinaryOutput() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::WriteOutputRuntimeVtpStructure() const
{
  CheckInitSetup();

  // ************** BEGIN RUNTIME VTP OUTPUT *** OPTION 1: DISCRETIZATION *********
  // this section compiles and seems to do the job correctly :-)

   // initialize the writer object
   vtp_writer_ptr_->SetGeometryFromParticleDiscretization();

   // reset time and time step of the writer object
   vtp_writer_ptr_->ResetTimeAndTimeStep( GState().GetTimeN(), GState().GetStepN() );

   // append all desired output data to the writer object's storage
   DRT::Discretization const& bindis = BinDiscret();
   Teuchos::RCP<Epetra_Vector> dis         = LINALG::CreateVector( *bindis.DofRowMap(), true );
   Teuchos::RCP<Epetra_Vector> orientation = LINALG::CreateVector( *bindis.DofRowMap(), true );
   Teuchos::RCP<Epetra_Vector> numbond     = LINALG::CreateVector( *bindis.NodeRowMap(), true );
   Teuchos::RCP<Epetra_Vector> owner       = LINALG::CreateVector( *bindis.NodeRowMap(), true );

   FillStateDataVectorsForOutput( dis, orientation, numbond, owner);

   // append displacement vector if desired
   // append displacement if desired
//   if ( GInOutput().GetRuntimeVtpOutputParams()->OutputDisplacementState() )
//     vtp_writer_ptr_->AppendDofBasedResultDataVector( dis, 3, "displacement" );

   // append owner if desired
   if ( GInOutput().GetRuntimeVtpOutputParams()->OutputOwner() )
     vtp_writer_ptr_->AppendNodeBasedResultDataVector( owner, 1, "owner" );

   // append orientation vector if desired
   if ( GInOutput().GetRuntimeVtpOutputParams()->OutputOrientation() )
     vtp_writer_ptr_->AppendDofBasedResultDataVector( orientation, 3, "orientation" );

   // append number of bonds if desired
   if ( GInOutput().GetRuntimeVtpOutputParams()->OutputNumberOfBonds() )
     vtp_writer_ptr_->AppendNodeBasedResultDataVector( numbond, 1, "numberofbonds" );

   // finalize everything and write all required VTU files to file system
   vtp_writer_ptr_->WriteFiles();

   // write collection files
   vtp_writer_ptr_->WriteCollectionFileOfAllWrittenFiles();



   // ************** BEGIN RUNTIME VTP OUTPUT *** OPTION 2: DIRECTLY *********
   // this section is just to get the idea and needs some minor modifications (indicated by Fixme)
   // this may serve as a template for visualization of stochastic&viscous forces, contact forces, ...

 //  Teuchos::RCP<RuntimeVtpWriter> vtp_writer_ptr =
 //      Teuchos::rcp( new RuntimeVtpWriter() );
 //
 //  // Todo: we need a better upper bound for total number of time steps here
 //  // however, this 'only' affects the number of leading zeros in the vtk file names
 //  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;
 //
 //
 //  // initialize the writer object
 //  vtp_writer_ptr->Initialize();   // Fixme
 //
 //  // set geometry manually
 //  const unsigned int num_spatial_dimensions = 3;
 //  unsigned int num_row_points = 2000;
 //
 //  // get and prepare storage for point coordinate values
 //  std::vector<double>& point_coordinates = vtp_writer_ptr->GetMutablePointCoordinateVector();
 //  point_coordinates.clear();
 //  point_coordinates.reserve( num_spatial_dimensions * num_row_points );
 //
 //
 //  // loop over my points and collect the geometry/grid data, i.e. reference positions of nodes
 //  for (unsigned int inode=0; inode < num_row_points; ++inode)
 //  {
 //    const DRT::Node* node = BinDiscretPtr()->lRowNode(inode);
 //
 //    for (unsigned int idim=0; idim<num_spatial_dimensions; ++idim)
 //      point_coordinates.push_back( node->X()[idim] );
 //  }
 //
 //
 //  // reset time and time step and geometry name in the writer object
 //  vtp_writer_ptr->SetupForNewTimeStepAndGeometry(
 //      GState().GetTimeN(), GState().GetStepN(), BinDiscretPtr()->Name() );
 //
 //
 //
 //  // append all desired output data to the writer object's storage
 //
 //  // number of bonds: collect data and append to visualization results if desired
 //  std::vector<double> num_bonds( num_row_points );
 //
 //  for ( unsigned int i = 0; i < num_row_points; ++i )
 //  {
 //    CROSSLINKING::CrosslinkerNode *crosslinker_i =
 //        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( bindis.lRowNode(i) );
 //
 //    num_bonds[i] = crosslinker_i->ClData()->GetNumberOfBonds();
 //  }
 //
 //  vtp_writer_ptr->AppendVisualizationPointDataVector( num_bonds, 1, "num_bonds" );
 //
 //
 //
 //  // finalize everything and write all required VTU files to filesystem
 //  vtp_writer_ptr->WriteFiles();
 //
 //
 //  // write a collection file summarizing all previously written files
 //  vtp_writer_ptr->WriteCollectionFileOfAllWrittenFiles();   // Fixme
   // ************** END RUNTIME VTP OUTPUT ***************************************
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FillStateDataVectorsForOutput(
    Teuchos::RCP<Epetra_Vector> displacement,
    Teuchos::RCP<Epetra_Vector> orientation,
    Teuchos::RCP<Epetra_Vector> numberofbonds,
    Teuchos::RCP<Epetra_Vector> owner
) const
{
  CheckInitSetup();
  DRT::Discretization const& bindis = BinDiscret();

  std::map< int, Teuchos::RCP< BEAMINTERACTION::BeamLink > >::const_iterator it;
  //todo: this is of course not nice, this needs to be done somewhere else
  for ( int i = 0; i < bindis.NumMyRowNodes(); ++i )
  {
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( bindis.lRowNode(i) );
    // std::vector holding gids of dofs
    std::vector<int> dofnode  = bindis.Dof(crosslinker_i);

    // loop over all dofs
    for ( int dim = 0; dim < 3; ++dim )
    {
      int doflid = displacement->Map().LID(dofnode[dim]);
      (*displacement)[doflid] = crosslinker_i->X()[dim];

      if (crosslinker_i->ClData()->GetNumberOfBonds() == 2)
      {
        it = doublebondcl_.find(crosslinker_i->Id());
        (*orientation)[doflid] = it->second->GetBindSpotPos1()(dim) - it->second->GetBindSpotPos2()(dim);
      }
      else
      {
        (*orientation)[doflid] = 0.0;
      }
    }

    (*numberofbonds)[i] = crosslinker_i->ClData()->GetNumberOfBonds();
    (*owner)[i] = crosslinker_i->Owner();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();
  // -------------------------------------------------------------------------
  // 1) write particle (restart) output (i.e. crosslinker nodes with corresponding data)
  // -------------------------------------------------------------------------
  OutputStepState( true );

  // -------------------------------------------------------------------------
  // 2) write list of double bonded crosslinker on each proc
  // -------------------------------------------------------------------------
  DRT::PackBuffer buffer;
  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP< BEAMINTERACTION::BeamLink > btbl = iter.second;
    btbl->Pack(buffer);
  }
  buffer.StartPacking();
  for ( auto const & iter : doublebondcl_ )
  {
    Teuchos::RCP< BEAMINTERACTION::BeamLink > btbl = iter.second;
    btbl->Pack(buffer);
  }

  Teuchos::RCP<std::vector<char> > block = Teuchos::rcp( new std::vector<char> );
  std::swap( *block, buffer() );

  // write data
  BinDiscret().Writer()->WriteCharVector( "BeamLink", block );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  // -------------------------------------------------------------------------
  // 1) rebuild bin discret correctly in case crosslinker were present
  // -------------------------------------------------------------------------
  ParticleHandlerPtr()->RemoveAllParticles();
  // read correct nodes
  IO::DiscretizationReader reader( BinDiscretPtr(), GState().GetStepN() );
  reader.ReadNodesOnly( GState().GetStepN() );
  // next, read step (as it was written next, do safety check)
  if ( GState().GetStepN() != reader.ReadInt("step") )
    dserror("Restart step of bin dis not consistent with that of problem discret ");

  BinDiscretPtr()->FillComplete( false, false, false );

  // -------------------------------------------------------------------------
  // 2) read list of double bonded crosslinker on each proc
  // -------------------------------------------------------------------------
  Teuchos::RCP< std::vector< char > > charvec;
  reader.ReadCharVector( charvec, "BeamLink" );

  std::vector<char>::size_type index = 0;
  while ( index < charvec->size() )
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack( index, *charvec, data );
    Teuchos::RCP< DRT::ParObject > object = Teuchos::rcp( DRT::UTILS::Factory(data), true );
    Teuchos::RCP< BEAMINTERACTION::BeamLink > beamtobeamlink =
        Teuchos::rcp_dynamic_cast< BEAMINTERACTION::BeamLink >(object);
    if ( beamtobeamlink == Teuchos::null )
      dserror("Failed to build a node from the node data");

    // insert in my list of double bonded crosslinker
    doublebondcl_[beamtobeamlink->Id()] = beamtobeamlink;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostReadRestart()
{
  CheckInitSetup();

  // bring each object in doublebondcl_ map to its correct owner
  UpdateMyDoubleBondsRemoteIdList();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::AddBinsToBinColMap(
    std::set< int >& colbins)
{
  // nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins ) const
{
  CheckInitSetup();

  std::set< int > colbinsext;
  std::set<int>::const_iterator iter;
  for ( iter = colbins.begin(); iter != colbins.end(); ++iter )
  {
    // get current bin
    DRT::Element* currbin = BinDiscret().gElement( *iter );
    // get all crosslinker in current bin
    DRT::Node **clincurrentbin = currbin->Nodes();

    // loop over all crosslinker in CurrentBin
    for ( int i = 0; i < currbin->NumNode(); ++i )
    {
      CROSSLINKING::CrosslinkerNode* crosslinker_i =
          dynamic_cast< CROSSLINKING::CrosslinkerNode* >( clincurrentbin[i] );

      if ( crosslinker_i->ClData()->GetNumberOfBonds() > 0 )
      {
        std::vector<int> binvec;
        ParticleHandler().BinStrategy().GetNeighborBinIds( *iter, binvec );
        colbinsext.insert( binvec.begin(), binvec.end() );
        // go to next bin
        break;
      }
    }
  }

  colbins.insert( colbinsext.begin(), colbinsext.end() );

//  // the following would be a default additional one ghost layer
//
//  std::set< int > colbinsext(colbins.begin(),colbins.end());
//  std::set< int >::const_iterator biniter;
//  for ( biniter = colbinsext.begin(); biniter != colbinsext.end() ; ++biniter )
//  {
//    std::vector<int> binvec;
//    BinStrategy().GetNeighborAndOwnBinIds( *biniter, binvec );
//    colbinsext.insert( binvec.begin(), binvec.end() );
//  }
//
//  colbins.insert( colbinsext.begin(), colbinsext.end() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DiffuseCrosslinker()
{
  CheckInitSetup();

  // loop over all row crosslinker (beam binding status not touched here)
  const int numrowcl = BinDiscretPtr()->NumMyRowNodes();
  for ( int rowcli = 0; rowcli < numrowcl; ++rowcli )
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode* crosslinker_i =
        dynamic_cast< CROSSLINKING::CrosslinkerNode* >(BinDiscretPtr()->lRowNode(rowcli));

#ifdef DEBUG
      if ( crosslinker_i == NULL )
        dserror("Dynamic cast to CrosslinkerNode failed");
      if ( crosslinker_i->NumElement() != 1 )
        dserror("More than one element for this crosslinker");
#endif

    Teuchos::RCP< CROSSLINKING::CrosslinkerNodeDataContainer > cldata_i = crosslinker_i->ClData();

    // different treatment according to number of bonds a crosslinker has
    switch( cldata_i->GetNumberOfBonds() )
    {
      case 0:
      {
        // crosslinker has zero bonds, i.e. is free to diffuse according to
        // brownian dynamics
        DiffuseUnboundCrosslinker( crosslinker_i );
        break;
      }
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = GetSingleOccupiedClBspot( cldata_i->GetClBSpotStatus() );

        // get current position of binding spot of filament partner
        // note: we can not use our beam data container, as bspot position is not current position (as this
        // is the result of a sum, you can not have a reference to that)
        const int elegid = cldata_i->GetClBSpotStatus()[occbspotid].first;

#ifdef DEBUG
        // safety check
        const int colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if ( colelelid < 0 )
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",
                  cldata_i->GetNumberOfBonds(),elegid,GState().GetMyRank());
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast< DRT::ELEMENTS::Beam3Base* >( DiscretPtr()->gElement(elegid) );

#ifdef DEBUG
        // safety check
        if ( ele == NULL)
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotpos;
        std::vector<double> eledisp;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp );
        ele->GetPosOfBindingSpot( bbspotpos, eledisp, cldata_i->GetClBSpotStatus()[occbspotid].second,
            PeriodicBoundingBox() );

        // note: a crosslinker can not leave the computational domain here, as no beam binding
        // spot can be outside the periodic box at this point
        SetCrosslinkerPosition( crosslinker_i, bbspotpos );

        break;
      }
      case 2:
      {
        // crosslinker has two bonds (cl gets current mid position between the filament
        // binding spot it is attached to)
        // -----------------------------------------------------------------
        // partner one
        // -----------------------------------------------------------------
        int elegid = cldata_i->GetClBSpotStatus()[0].first;

#ifdef DEBUG
        if ( elegid < 0 or cldata_i->GetClBSpotStatus()[0].second < 0 )
          dserror(" double bonded crosslinker has stored beam partner gid or loc bsponum of -1, "
                  " something went wrong");
        // safety check
        int colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if ( colelelid < 0 )
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is not \nghosted/owned on proc %i (owner of crosslinker)",cldata_i->GetNumberOfBonds(),elegid,GState().GetMyRank());
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

#ifdef DEBUG
        // safety check
        if ( ele == NULL)
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotposone;
        std::vector<double> eledisp;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp );
        ele->GetPosOfBindingSpot( bbspotposone, eledisp, cldata_i->GetClBSpotStatus()[0].second,
            PeriodicBoundingBox() );

        // -----------------------------------------------------------------
        // partner two
        // -----------------------------------------------------------------
        elegid = cldata_i->GetClBSpotStatus()[1].first;

#ifdef DEBUG
        // safety check
        if ( elegid < 0 or cldata_i->GetClBSpotStatus()[1].second < 0 )
          dserror(" double bonded crosslinker has stored beam partner gid or loc bsponum of -1, "
                  " something went wrong");
        colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if ( colelelid < 0 )
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",cldata_i->GetNumberOfBonds(),elegid,GState().GetMyRank());
#endif

        ele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

#ifdef DEBUG
        // safety check
        if ( ele == NULL )
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotpostwo;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);
        ele->GetPosOfBindingSpot( bbspotpostwo, eledisp, cldata_i->GetClBSpotStatus()[1].second,
            PeriodicBoundingBox() );

        LINALG::Matrix< 3, 1 > clpos( true );
        SetPositionOfDoubleBondedCrosslinkerPBCconsistent( clpos, bbspotposone, bbspotpostwo );

        SetCrosslinkerPosition( crosslinker_i, clpos );

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i->GetNumberOfBonds() );
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DiffuseUnboundCrosslinker(
    CROSSLINKING::CrosslinkerNode* crosslinker) const
{
  CheckInit();

  // get standard deviation and mean value for crosslinker that are free to
  // diffuse
  double standarddev = std::sqrt( 2.0 * crosslinking_params_ptr_->KT() /
                       ( 3.0*M_PI * crosslinking_params_ptr_->Viscosity()
                       * crosslinker->GetMaterial()->LinkingLength() )
                       * crosslinking_params_ptr_->DeltaTime() );
  double meanvalue = 0.0;
  // Set mean value and standard deviation of normal distribution
  // FixMe standard deviation = sqrt(variance) check this for potential error !!!
  DRT::Problem::Instance()->Random()->SetMeanVariance( meanvalue, standarddev );

  // diffuse crosslinker according to brownian dynamics
  LINALG::Matrix<3,1> newclpos ( true );
  std::vector<double> randvec;
  int count = 3;
  // maximal diffusion given by cutoff radius (sqrt(3) = 1.73..)
  double const maxmov = BinStrategy().CutoffRadius() / 1.74;
  DRT::Problem::Instance()->Random()->Normal( randvec, count );
  for ( int dim = 0; dim < 3; ++dim )
  {
    if ( abs(randvec[dim]) > maxmov )
    {
      double old = randvec[dim];
      randvec[dim] = ( abs(randvec[dim]) / randvec[dim] ) * maxmov;
      std::cout << "Movement of free crosslinker " << crosslinker->Id() << " was restricted by cutoff radius"
          " in " << dim << " direction. " << old << " to " << randvec[dim] << "\nThis should not happen to often "
          "to stay physical. Increase cutoff or reduce movement" << std::endl;
    }
    newclpos(dim) = crosslinker->X()[dim] + randvec[dim];
  }

  // check compliance with periodic boundary conditions
  PeriodicBoundingBox().Shift3D( newclpos );
  SetCrosslinkerPosition( crosslinker, newclpos );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
int BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetSingleOccupiedClBspot(
    std::vector<std::pair<int, int> > const& clbspots) const
{
  CheckInit();

  if ( clbspots[0].first > -1)
    return 0;
  else if ( clbspots[1].first > -1 )
    return 1;
  else
    dserror("numbond = 1 but both binding spots store invalid element GIDs!");

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    SetPositionOfDoubleBondedCrosslinkerPBCconsistent(
        LINALG::Matrix<3,1>& clpos,
        LINALG::Matrix<3,1> const& bspot1pos,
        LINALG::Matrix<3,1> const& bspot2pos ) const
{
  /* the position of (the center) of a double-bonded crosslinker is defined as
   * midpoint between the two given binding spot positions. (imagine a linker
   * being a slender body with a binding domain at each of both ends) */

  /* if the two binding spots are separated by a periodic boundary, we need to
   * shift one position back to get the interpolation right */
  clpos = bspot2pos;
  PeriodicBoundingBox().UnShift3D( clpos, bspot1pos );

  // fixme: to avoid senseless dserror in debug mode
  LINALG::Matrix<3,1> dummy(clpos);
  clpos.Update( 0.5, bspot1pos, 0.5, dummy);

  // shift the interpolated position back in the periodic box if necessary
  PeriodicBoundingBox().Shift3D(clpos);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetCrosslinkerPosition(
    DRT::Node* crosslinker,
    LINALG::Matrix<3,1> const & newclpos) const
{
  CheckInit();

//  double const& cutoff = BinStrategy().CutoffRadius();
//  double dist = 0.0;
//  for ( int dim = 0; dim < 3; ++dim )
//  {
//    dist = abs( newclpos(dim) - crosslinker->X()[dim] );
//    if ( dist > cutoff and (PeriodicBoundingBox().EdgeLength(dim) - dist) > cutoff )
//      dserror(" Crosslinker %i moved further than the bin length. %f > %f",
//          crosslinker->Id(), (dist > PeriodicBoundingBox().EdgeLength(dim)/2.0) ?
//              (PeriodicBoundingBox().EdgeLength(dim) - dist) : dist, cutoff );
//  }

  std::vector<double> newpos(3,0.0);
  for ( int dim = 0; dim < 3; ++dim )
    newpos[dim] = newclpos(dim);
  crosslinker->SetPos(newpos);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetPositionOfNewlyFreeCrosslinker(
    CROSSLINKING::CrosslinkerNode* crosslinker,
    LINALG::Matrix<3,1>& clpos) const
{
  CheckInit();

  //generate vector in random direction
  // of length half the linking length to "reset" crosslink molecule position: it may now
  // reenter or leave the bonding proximity
  // todo: does this make sense?
  LINALG::Matrix<3,1> cldeltapos_i;
  std::vector<double> randunivec(3);
  int count = 3;
  DRT::Problem::Instance()->Random()->Uni( randunivec, count );
  for ( unsigned int dim = 0 ; dim < 3; ++dim )
    cldeltapos_i(dim) = randunivec[dim];

  cldeltapos_i.Scale( crosslinker->GetMaterial()->LinkingLength() / cldeltapos_i.Norm2() );

  clpos.Update( 1.0, cldeltapos_i, 1.0 );

  PeriodicBoundingBox().Shift3D(clpos);
  SetCrosslinkerPosition( crosslinker, clpos );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetPositionOfNewlySingleBondedCrosslinker(
    DRT::Node* crosslinker,
    CrosslinkerData& cldata,
    int const stayoccpotid)
{
  CheckInit();

  // update postion
  const int collidoccbeam =
      DiscretPtr()->ElementColMap()->LID(cldata.clbspots[stayoccpotid].first);
#ifdef DEBUG
  // safety check
  if ( collidoccbeam < 0 )
    dserror("element with gid %i not ghosted on proc %i",cldata.clbspots[stayoccpotid].first, GState().GetMyRank() );
#endif
  BeamData& beamdata_i = beam_data_[collidoccbeam];
  cldata.clpos = beamdata_i.bbspotpos[cldata.clbspots[stayoccpotid].second];
  SetCrosslinkerPosition( crosslinker, cldata.clpos );
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreComputeCrosslinkerAndBeamData(
    std::vector<CrosslinkerData> & crosslinker_data,
    std::vector<BeamData> & beam_data)
{
  CheckInit();
  PreComputeCrosslinkerData(crosslinker_data);
  PreComputeBeamData(beam_data);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreComputeCrosslinkerData(
    std::vector<CrosslinkerData> & crosslinker_data )
{
  CheckInit();

  crosslinker_data.clear();
  const int numcolcl = BinDiscretPtr()->NumMyColNodes();
  crosslinker_data.resize(numcolcl);

  for ( int i = 0; i < numcolcl; ++i )
  {
    // crosslinker i for which data will be collected
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lColNode(i));

#ifdef DEBUG
      if (crosslinker_i == NULL)
        dserror("Dynamic cast to CrosslinkerNode failed");
      if (crosslinker_i->NumElement() != 1)
        dserror("More than one element for this crosslinker");
#endif

    // store data of crosslinker i according to column lid
    CrosslinkerData& cldata = crosslinker_data[i];

    // store positions
    for ( int dim = 0; dim < 3; ++dim )
      cldata.clpos(dim) = crosslinker_i->X()[dim];
    // get current binding spot status of crosslinker
    cldata.clbspots = crosslinker_i->ClData()->GetClBSpotStatus();
    // get number of bonds
    cldata.clnumbond = crosslinker_i->ClData()->GetNumberOfBonds();
    // get type of crosslinker (i.e. its material)
    cldata.clmat = crosslinker_i->GetMaterial();
    // get owner
    cldata.clowner = crosslinker_i->Owner();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreComputeBeamData(
    std::vector<BeamData> & beam_data)
{
  CheckInit();

  beam_data.clear();
  const int numcoleles = DiscretPtr()->NumMyColElements();
  beam_data.resize(numcoleles);

  // loop over all column beams elements
  for ( int i = 0; i < numcoleles; ++i )
  {
    // beam element i for which data will be collected
    DRT::ELEMENTS::Beam3Base* beamele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->lColElement(i));

    // go to next element in case the current one is not a beam element
    if ( beamele_i == NULL )
      continue;

    // store data
    BeamData& bdata = beam_data[i];

    std::vector<double> eledisp;
    BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), beamele_i,
        BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp );

    // loop over all binding spots of current element
    const int numbbspot = static_cast<int>(beamele_i->GetBindingSpotStatus().size());
    for (int j=0; j<numbbspot; ++j)
      BEAMINTERACTION::UTILS::GetPosAndTriadOfBindingSpot( beamele_i,
          PeriodicBoundingBoxPtr(), j, bdata.bbspotpos[j], bdata.bbspottriad[j], eledisp );

    // get status of beam binding spots
    bdata.bbspotstatus = beamele_i->GetBindingSpotStatus();
    // get owner
    bdata.bowner = beamele_i->Owner();
    // get filament type
    bdata.filamenttype = beamele_i->GetFilamentType();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::BindCrosslinker()
{
  CheckInit();

  DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );

  // intended bonds of row crosslinker on myrank (key is clgid)
  std::map< int, BindEventData > mybonds;
  // intended bond col crosslinker to row element (key is owner of crosslinker != myrank)
  std::map< int, std::vector<BindEventData> > undecidedbonds;

  // fill binding event maps
  FindPotentialBindingEvents( mybonds, undecidedbonds );

  // bind events where myrank only owns the elements, cl are taken care
  // of by their owner (key is clgid)
  std::map< int, BindEventData > myelebonds;

  // now each row owner of a linker gets requests, makes a random decision and
  // informs back its requesters
  ManageBindingInParallel( mybonds, undecidedbonds, myelebonds );

  // actual update of binding states is done here
  UpdateMyCrosslinkerAndElementBindingStates( mybonds, myelebonds );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FindPotentialBindingEvents(
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FindPotentialBindingEvents");

  // this variable is used to check if a beam binding spot is linked twice on
  // myrank during a time step
  std::vector< std::set <int> > intendedbeambonds( DiscretPtr()->NumMyRowElements() );

  // store bins that have already been examined
  std::unordered_set<int> examinedbins;
  // loop over all column crosslinker in random order
  // create random order of indices
  std::vector<int> rordercolcl =
      BEAMINTERACTION::UTILS::Permutation( BinDiscretPtr()->NumMyColNodes() );
  std::vector<int>::const_iterator icl;
  for ( icl = rordercolcl.begin(); icl != rordercolcl.end(); ++icl )
  {
    DRT::Node *currcrosslinker = BinDiscretPtr()->lColNode( *icl );

#ifdef DEBUG
    if ( currcrosslinker == NULL )
      dserror("Dynamic cast to CrosslinkerNode failed");
    if ( currcrosslinker->NumElement() != 1 )
      dserror("More than one element for this crosslinker");
#endif

    // get bin that contains this crosslinker (can only be one)
    DRT::Element* currentbin = currcrosslinker->Elements()[0];
    const int currbingid = currentbin->Id();

#ifdef DEBUG
    if ( currbingid < 0 )
      dserror(" negative bin id number %i ", currbingid );
#endif

    // if a bin has already been examined --> continue with next crosslinker
    if ( examinedbins.find(currbingid) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins.insert(currbingid);

    FindPotentialBindingEventsInBinAndNeighborhood( currentbin, mybonds,
        undecidedbonds, intendedbeambonds, true );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    FindPotentialBindingEventsInBinAndNeighborhood(
    DRT::Element* bin,
    std::map<int, BindEventData >& mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds,
    std::vector< std::set <int> >& intendedbeambonds,
    bool const checklinkingprop)
{
  CheckInit();

  // get neighboring bins
  // note: interaction distance cl to beam needs to be smaller than the bin size
  std::vector<int> neighboring_binIds;
  neighboring_binIds.reserve(27);
  // do not check on existence here -> shifted to GetBinContent
  BinStrategyPtr()->GetNeighborAndOwnBinIds( bin->Id(), neighboring_binIds );

  // get set of neighboring beam elements (i.e. elements that somehow touch nb bins)
  // as explained above, we only need row elements
  std::set<DRT::Element*> neighboring_beams;
  std::vector< BINSTRATEGY::UTILS::BinContentType > bc( 1, BINSTRATEGY::UTILS::Beam );
  BinStrategyPtr()->GetBinContent( neighboring_beams, bc, neighboring_binIds, true );

  // in case there are no neighbors, go to next crosslinker (an therefore bin)
  if ( neighboring_beams.empty() )
    return;

  // get all crosslinker in current bin
  DRT::Node **clincurrentbin = bin->Nodes();
  const int numcrosslinker = bin->NumNode();

  // obtain random order in which crosslinker are addressed
  std::vector<int> randorder = BEAMINTERACTION::UTILS::Permutation( numcrosslinker );

  // loop over all crosslinker in CurrentBin in random order
  std::vector<int>::const_iterator randcliter;
  for ( randcliter = randorder.begin(); randcliter != randorder.end(); ++randcliter )
  {
    // get random crosslinker in current bin
    DRT::Node *crosslinker_i = clincurrentbin[*randcliter];
    // get all potential binding events on myrank
    PrepareBinding( crosslinker_i, neighboring_beams, mybonds, undecidedbonds,
        intendedbeambonds, checklinkingprop );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareBinding(
    DRT::Node *                                  crosslinker_i,
    std::set<DRT::Element*> const &              neighboring_beams,
    std::map<int, BindEventData > &              mybonds,
    std::map<int, std::vector<BindEventData> > & undecidedbonds,
    std::vector< std::set <int> > &              intendedbeambonds,
    bool const checklinkingprop)
{
  CheckInit();

  // get precomputed data of crosslinker i
  CrosslinkerData const& cldata_i = crosslinker_data_[crosslinker_i->LID()];

  // -------------------------------------------------------------------------
  // We now check all criteria that need to be passed for a binding event one
  // after the other
  // -------------------------------------------------------------------------
  // 1. criterion: in case crosslinker is double bonded, we can leave here
  if ( cldata_i.clnumbond == 2 )
    return;

  // loop over all neighboring beam elements in random order (keep in mind
  // we are only looping over row elements)
  std::vector< DRT::Element* > beamvec( neighboring_beams.begin(), neighboring_beams.end() );

  // -------------------------------------------------------------------------
  // NOTE: This is crucial for reproducibility to ensure that computation does
  // not depend on pointer addresses (see also comment of class Less)
  // -------------------------------------------------------------------------
  std::sort( beamvec.begin(), beamvec.end(), BEAMINTERACTION::UTILS::Less() );

  std::vector<int> randorder = BEAMINTERACTION::UTILS::Permutation( static_cast<int>( beamvec.size() ) );
  std::vector<int> ::const_iterator randiter;
  for ( randiter = randorder.begin(); randiter != randorder.end();  ++randiter )
  {
    // get neighboring (nb) beam element
    DRT::Element* nbbeam = beamvec[*randiter];

#ifdef DEBUG
    // some safety checks
    if ( nbbeam->LID() < 0 )
      dserror("Beam lid < 0");
#endif

    // get pre computed data of current nbbeam
    BeamData const& beamdata_i = beam_data_[ nbbeam->LID() ];

    // 2. criterion:
    // exclude binding of a single bonded crosslinker in close proximity on the
    // same filament (i.e. element cloud of old element binding partner is excluded)
    if ( cldata_i.clnumbond == 1 and
        BEAMINTERACTION::UTILS::DoBeamElementsShareNodes( Discret().gElement(
            cldata_i.clbspots[GetSingleOccupiedClBspot( cldata_i.clbspots )].first ), nbbeam ) )
      continue;

    // loop over all binding spots of current element in random order
    std::vector<int> randbspot = BEAMINTERACTION::UTILS::Permutation( static_cast<int>( beamdata_i.bbspotstatus.size() ) );
    std::vector<int> ::const_iterator rbspotiter;
    for ( rbspotiter = randbspot.begin(); rbspotiter != randbspot.end(); ++rbspotiter )
    {
      // get local number of binding spot in element
      const int locnbspot = *rbspotiter;

      // we are now doing some additional checks if a binding event is feasible
      if ( not CheckBindEventCriteria( dynamic_cast< CROSSLINKING::CrosslinkerNode* >( crosslinker_i ), nbbeam,
                                       cldata_i, beamdata_i, locnbspot, intendedbeambonds, checklinkingprop )
         )
        continue;

      int const beamrowlid = Discret().ElementRowMap()->LID( nbbeam->Id() );
      // insert current event
      intendedbeambonds[beamrowlid].insert(locnbspot);

      // ---------------------------------------------------------------------
      // if we made it this far, we can add this potential binding event to its
      // corresponding map
      // ---------------------------------------------------------------------
      BindEventData bindeventdata;
      bindeventdata.clgid = crosslinker_i->Id();
      bindeventdata.elegid = nbbeam->Id();
      bindeventdata.bspotlocn = locnbspot;
      bindeventdata.requestproc = GState().GetMyRank();
      // this is default, is changed if owner of cl has something against it
      bindeventdata.permission = 1;

      // in case myrank is owner, we add it to the mybonds map
      if ( cldata_i.clowner == GState().GetMyRank() )
      {
        mybonds[bindeventdata.clgid] = bindeventdata;
      }
      else
      {
        // myrank is not owner, we add it to the map of events that need to be
        // communicated to make a decision
        undecidedbonds[cldata_i.clowner].push_back(bindeventdata);
      }

      // as we allow only one binding event for each cl in one time step,
      // we are done here, if we made it so far (i.e met criteria 1. - 7.)
      return;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CheckBindEventCriteria(
    CROSSLINKING::CrosslinkerNode const * const crosslinker_i,
    DRT::Element const * const potbeampartner,
    CrosslinkerData const& cldata_i,
    BeamData const& beamdata_i,
    int locnbspot,
    std::vector< std::set <int> >& intendedbeambonds,
    bool checklinkingprop) const
{
  CheckInit();

  // 3. criterion:
  // check compatibility of crosslinker type and filament type (some linker can only
  // bind to certain filament types)
  if ( not CheckLinkerAndFilamentTypeCompatibility( crosslinker_i->GetMaterial()->LinkingType(), beamdata_i.filamenttype ) )
    return false;

  // 4. criterion:
  // a crosslink is set if and only if it passes the probability check
  // for a binding event to happen
  double plink = 1.0 - exp( (-1.0) * crosslinking_params_ptr_->DeltaTime() * crosslinker_i->GetMaterial()->KOn() );

  if ( checklinkingprop and ( DRT::Problem::Instance()->Random()->Uni() > plink ) )
    return false;

  // 5. criterion:
  // first check if binding spot is free, if not, check next bspot on curr ele
  // note: bspotstatus in bonded case holds cl gid, otherwise -1 (meaning free)
  if ( beamdata_i.bbspotstatus.at(locnbspot) != -1 )
    return false;

  /* 6. criterion: check RELEVANT distance criterion
   * if free:
   *   distance between crosslinker center and current beam binding spot
   * if singly bound:
   *   distance between already bound bspot of crosslinker and current beam binding spot
   * note: as we set the crosslinker position to coincide with beam bspot position if singly bound,
   *       we can also use cldata_i.clpos in the second case*/

  // get current position and tangent vector of filament axis at free binding spot
  LINALG::Matrix<3,1> const& currbbspos = beamdata_i.bbspotpos.at(locnbspot);

  // minimum and maximum distance at which a double-bond crosslink can be established
  double const linkdistmin = crosslinker_i->GetMaterial()->LinkingLength()
                             - crosslinker_i->GetMaterial()->LinkingLengthTolerance();
  double const linkdistmax = crosslinker_i->GetMaterial()->LinkingLength()
                             + crosslinker_i->GetMaterial()->LinkingLengthTolerance();

#ifdef DEBUG
  //safety check
  if ( linkdistmax > BinStrategy().CutoffRadius() )
    dserror("The allowed binding distance of linker %i (in case it is single bonded) is"
            "\ngreater than the cutoff radius, this could lead to missing a binding event", crosslinker_i->Id() );
#endif

  if ( ( cldata_i.clnumbond == 0 and BEAMINTERACTION::UTILS::IsDistanceOutOfRange( cldata_i.clpos,
      currbbspos, 0.5 * linkdistmin, 0.5 * linkdistmax ) )
      or
      ( cldata_i.clnumbond == 1 and BEAMINTERACTION::UTILS::IsDistanceOutOfRange( cldata_i.clpos,
      currbbspos, linkdistmin, linkdistmax ) )
    )
    return false;

  // 7. criterion: orientation of centerline tangent vectors at binding spots
  // a crosslink (double-bonded crosslinker) will only be established if the
  // enclosed angle is in the specified range
  double const linkanglemin = crosslinker_i->GetMaterial()->LinkingAngle()
                              - crosslinker_i->GetMaterial()->LinkingAngleTolerance();
  double const linkanglemax = crosslinker_i->GetMaterial()->LinkingAngle()
                              + crosslinker_i->GetMaterial()->LinkingAngleTolerance();

  // if crosslinker is singly bound, we fetch the orientation vector
  LINALG::Matrix<3,1> occ_bindingspot_beam_tangent(true);
  if ( cldata_i.clnumbond == 1 )
    GetOccupiedClBSpotBeamTangent( cldata_i, occ_bindingspot_beam_tangent , crosslinker_i->Id() );

  // note: we use first base vector instead of tangent vector here
  LINALG::Matrix<3,1> curr_bindingspot_beam_tangent(true);
  for ( unsigned int idim = 0; idim < 3; ++idim )
    curr_bindingspot_beam_tangent(idim) = beamdata_i.bbspottriad.at(locnbspot)(idim,0);

  if ( cldata_i.clnumbond == 1 and BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
      occ_bindingspot_beam_tangent, curr_bindingspot_beam_tangent, linkanglemin, linkanglemax ) )
    return false;

  // 8. criterion
  // check if current beam binding spot yet intended to bind this timestep
  // by a crosslinker that came before in this random order
  int const beamrowlid = Discret().ElementRowMap()->LID( potbeampartner->Id() );
  if ( intendedbeambonds[beamrowlid].find( locnbspot ) != intendedbeambonds[beamrowlid].end() )
  {
    /* note: it is possible that the binding event that rejects the current one is rejected itself
     * later during communication with other procs and therefore the current one could be
     * valid. Just neglecting this here is a slight inconsistency, but should be ok as such an
     * coincidence is extremely rare in a simulation with realistic proportion of crosslinker
     * to beam binding spots. Additionally missing one event would not change any physics.
     * (Could be cured with additional communication)
     */
    if ( Discret().Comm().NumProc() > 1 )
      std::cout << " Warning: There is a minimal chance of missing a regular binding event on "
          "rank " << GState().GetMyRank() << std::endl;
    return false;
  }

  // bind event can happen
  // (rejection afterwards only possible in case of parallel simulation)
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CheckLinkerAndFilamentTypeCompatibility(
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype,
    INPAR::BEAMINTERACTION::FilamentType filamenttype) const
{
  switch ( linkertype )
  {
    case INPAR::BEAMINTERACTION::linkertype_arbitrary:
    {
      // no check of filament type necessary
      return true;
      break;
    }
    case INPAR::BEAMINTERACTION::linkertype_actin:
    {
      if ( filamenttype == INPAR::BEAMINTERACTION::filtype_actin )
        return true;
      else
        return false;
      break;
    }
    case INPAR::BEAMINTERACTION::linkertype_collagen:
    {
      if ( filamenttype == INPAR::BEAMINTERACTION::filtype_collagen )
        return true;
      else
        return false;
      break;
    }
    default:
    {
      dserror("Unknown linker type.");
      exit(EXIT_FAILURE);
    }
  }

  // default false
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetOccupiedClBSpotBeamTangent(
    CrosslinkerData const& cldata_i,
    LINALG::Matrix<3,1>& occ_bindingspot_beam_tangent,
    int const clgid ) const
{
  CheckInitSetup();

  int occbspotid = GetSingleOccupiedClBspot( cldata_i.clbspots );

  const int locbspotnum = cldata_i.clbspots[occbspotid].second;
  const int elegid = cldata_i.clbspots[occbspotid].first;
  const int elecollid = Discret().ElementColMap()->LID( elegid );

#ifdef DEBUG
  if ( elecollid < 0 )
    dserror (" Element with gid %i bonded to cl %i on rank %i not even ghosted",
             elegid, clgid, GState().GetMyRank() );
#endif

  // note: we use first base vector instead of tangent vector here
  for ( unsigned int idim = 0; idim < 3; ++idim )
    occ_bindingspot_beam_tangent(idim) =
        beam_data_[elecollid].bbspottriad.at(locbspotnum)(idim,0);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ManageBindingInParallel(
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds,
    std::map<int, BindEventData >&              myelebonds) const
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ManageBindingInParallel");

  // variable for safety check
  int numrecrequest;
  // exporter
  DRT::Exporter exporter( BinDiscret().Comm() );

  // -------------------------------------------------------------------------
  // 1) each procs makes his requests and receives the request of other procs
  // -------------------------------------------------------------------------
  // store requested cl and its data
  std::map<int, std::vector<BindEventData> > requestedcl;
  CommunicateUndecidedBonds( exporter, undecidedbonds, numrecrequest, requestedcl );

  // -------------------------------------------------------------------------
  // 2) now myrank needs to decide which proc is allowed to set the requested
  //    link
  // -------------------------------------------------------------------------
  std::map<int, std::vector<BindEventData> > decidedbonds;
  DecideBindingInParallel( requestedcl, mybonds, decidedbonds );

  // -------------------------------------------------------------------------
  // 3) communicate the binding decisions made on myrank, receive decisions
  //    made for its own requests and create colbondmap accordingly
  // -------------------------------------------------------------------------
  int answersize = static_cast<int>(undecidedbonds.size());
  CommunicateDecidedBonds( exporter, decidedbonds, myelebonds, numrecrequest, answersize );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyCrosslinkerAndElementBindingStates(
    std::map<int, BindEventData >& mybonds,
    std::map<int, BindEventData >& myelebonds)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::BindCrosslinker");

  // map key is crosslinker gid to be able to uniquely address one entry over all procs
  std::map< int, NewDoubleBonds > mynewdbondcl;

  // myrank owner of crosslinker and most elements
  UpdateMyCrosslinkerBindingStates( mybonds, mynewdbondcl );

  // myrank only owner of current binding partner ele
  UpdateMyElementBindingStates( myelebonds );

  // setup new double bonds and insert them in doublebondcl_
  CreateNewDoubleBondedCrosslinkerElementPairs( mynewdbondcl );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyCrosslinkerBindingStates(
    std::map<int, BindEventData > const& mybonds,
    std::map<int, NewDoubleBonds>&       mynewdbondcl)
{
  CheckInit();

  std::map<int, BindEventData>::const_iterator cliter;
  for ( cliter = mybonds.begin(); cliter != mybonds.end(); ++cliter )
  {
    // get binding event data
    BindEventData binevdata = cliter->second;

#ifdef DEBUG
    if ( binevdata.permission != 1)
      dserror(" Rank %i wants to bind crosslinker %i without permission, "
              " something went wrong", GState().GetMyRank(), cliter->first);
#endif

    // get current linker and beam data
    const int clcollid = BinDiscretPtr()->NodeColMap()->LID(cliter->first);
    const int colelelid = DiscretPtr()->ElementColMap()->LID(binevdata.elegid);

#ifdef DEBUG
    // safety checks
    if ( clcollid < 0 )
      dserror("Crosslinker not even ghosted, should be owned here.");
    if ( colelelid < 0 )
      dserror("Element with gid %i not ghosted.", binevdata.elegid );
#endif

    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lColNode(clcollid));
    // get crosslinker data
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    DRT::ELEMENTS::Beam3Base* beamele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->lColElement(colelelid));
    BeamData& beamdata_i = beam_data_[colelelid];

#ifdef DEBUG
    // safety checks
    if ( cliter->first != binevdata.clgid )
      dserror("Map key does not match crosslinker gid of current binding event.");

    if ( cldata_i.clowner != GState().GetMyRank() )
      dserror("Only row owner of crosslinker is changing its status");

    if ( colelelid < 0 )
      dserror("Binding element partner of current row crosslinker is not ghosted, "
              "this must be the case though.");
#endif

    // -------------------------------------------------------------------------
    // different treatment according to number of bonds crosslinker had before
    // this binding event
    // -------------------------------------------------------------------------
    switch( cldata_i.clnumbond )
    {
      case 0:
      {
        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element, first binding spot
        // always bonded first
        cldata_i.clbspots[0].first  = binevdata.elegid;
        cldata_i.clbspots[0].second = binevdata.bspotlocn;
        crosslinker_i->ClData()->SetClBSpotStatus( cldata_i.clbspots );

        // update number of bonds
        cldata_i.clnumbond = 1;
        crosslinker_i->ClData()->SetNumberOfBonds( cldata_i.clnumbond );

        // update position
        cldata_i.clpos = beamdata_i.bbspotpos[ binevdata.bspotlocn ];
        SetCrosslinkerPosition( crosslinker_i, cldata_i.clpos );

        // -----------------------------------------------------------------
        // update beam status
        // -----------------------------------------------------------------
        // store crosslinker gid in status of beam binding spot if myrank
        // is owner of beam
        if ( beamdata_i.bowner == GState().GetMyRank() )
        {
          beamdata_i.bbspotstatus[ binevdata.bspotlocn ] = binevdata.clgid;
          beamele_i->SetBindingSpotStatus( beamdata_i.bbspotstatus );
        }

#ifdef DEBUG
        // safety check
        if (not (cldata_i.clbspots[1].first < 0) )
          dserror("Numbond does not fit to clbspot vector.");
#endif

        break;
      }
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = GetSingleOccupiedClBspot( cldata_i.clbspots );
        int freebspotid = 1;
        if ( occbspotid == 1 )
          freebspotid = 0;

        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element
        cldata_i.clbspots[freebspotid].first = binevdata.elegid;
        cldata_i.clbspots[freebspotid].second = binevdata.bspotlocn;
        crosslinker_i->ClData()->SetClBSpotStatus( cldata_i.clbspots );

        // update number of bonds
        cldata_i.clnumbond = 2;
        crosslinker_i->ClData()->SetNumberOfBonds(cldata_i.clnumbond);

        // update position
        const LINALG::Matrix<3,1> occbspotpos_copy = cldata_i.clpos;
        SetPositionOfDoubleBondedCrosslinkerPBCconsistent(
            cldata_i.clpos,
            beamdata_i.bbspotpos[ cldata_i.clbspots[freebspotid].second ],
            occbspotpos_copy);

        SetCrosslinkerPosition( crosslinker_i, cldata_i.clpos );

        // create double bond cl data
        int occ_colelelid = DiscretPtr()->ElementColMap()->LID(cldata_i.clbspots[occbspotid].first);
        NewDoubleBonds dbondcl;
        dbondcl.id = binevdata.clgid;
        if ( cldata_i.clbspots[freebspotid].first > cldata_i.clbspots[occbspotid].first)
        {
          dbondcl.eleids.push_back( cldata_i.clbspots[freebspotid] );
          dbondcl.eleids.push_back( cldata_i.clbspots[occbspotid] );
          dbondcl.bspotposs.push_back( beam_data_[ colelelid ].bbspotpos.at(cldata_i.clbspots[freebspotid].second) );
          dbondcl.bspotposs.push_back( beam_data_[ occ_colelelid ].bbspotpos.at(cldata_i.clbspots[occbspotid].second) );
          dbondcl.bspottriads.push_back( beam_data_[ colelelid ].bbspottriad.at(cldata_i.clbspots[freebspotid].second) );
          dbondcl.bspottriads.push_back( beam_data_[ occ_colelelid ].bbspottriad.at(cldata_i.clbspots[occbspotid].second) );
        }
        else
        {
          dbondcl.eleids.push_back( cldata_i.clbspots[occbspotid] );
          dbondcl.eleids.push_back( cldata_i.clbspots[freebspotid] );
          dbondcl.bspotposs.push_back( beam_data_[ occ_colelelid ].bbspotpos.at(cldata_i.clbspots[occbspotid].second) );
          dbondcl.bspotposs.push_back( beam_data_[ colelelid ].bbspotpos.at(cldata_i.clbspots[freebspotid].second) );
          dbondcl.bspottriads.push_back( beam_data_[ occ_colelelid ].bbspottriad.at(cldata_i.clbspots[occbspotid].second) );
          dbondcl.bspottriads.push_back( beam_data_[ colelelid ].bbspottriad.at(cldata_i.clbspots[freebspotid].second) );
        }

        // insert pair in mypairs
        mynewdbondcl[dbondcl.id] = dbondcl;

        // first check if myrank is owner of element of current binding event
        // (additionally to being owner of cl)
        if ( beamdata_i.bowner == GState().GetMyRank() )
        {
          // update beam data
          // store crosslinker gid in status of beam binding spot
          beamdata_i.bbspotstatus[ binevdata.bspotlocn ] = binevdata.clgid;
          beamele_i->SetBindingSpotStatus( beamdata_i.bbspotstatus );
        }

        break;
      }
      default:
      {
        dserror("You should not be here, crosslinker has unrealistic number "
                "%i of bonds.", cldata_i.clnumbond);
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyElementBindingStates(
    std::map<int, BindEventData > const& myelebonds)
{
  CheckInit();

  /*
   * 1| 2__|2  or  1| 3__|2  or  1| 2__|1
   * 1|    |2      1|    |2      1|    |1
   * legend: | = beam; __= cl; 2,3 = owner; 1 = myrank
   */
  // loop through all binding events
  std::map<int, BindEventData>::const_iterator cliter;
  for ( cliter = myelebonds.begin(); cliter != myelebonds.end(); ++cliter )
  {
    // get binding event data
    BindEventData binevdata = cliter->second;

    // get linker data and beam data
    const int clcollid = BinDiscretPtr()->NodeColMap()->LID(cliter->first);
    const int colelelid = DiscretPtr()->ElementColMap()->LID(binevdata.elegid);

#ifdef DEBUG
    // safety checks
    if ( clcollid < 0 )
     dserror("Crosslinker needs to be ghosted, but this isn't the case.");
    if ( colelelid < 0 )
     dserror("element with gid %i not ghosted on proc %i", binevdata.elegid, GState().GetMyRank() );
#endif

    // linker
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // beam binding partner
    DRT::ELEMENTS::Beam3Base* ele_i =
       dynamic_cast<DRT::ELEMENTS::Beam3Base*>( DiscretPtr()->lColElement(colelelid) );
    BeamData& beamdata_i = beam_data_[colelelid];

#ifdef DEBUG
    // safety checks
    if ( beamdata_i.bowner != GState().GetMyRank() )
     dserror("Only row owner of element is allowed to change its status");
    if ( cldata_i.clowner == GState().GetMyRank() )
     dserror("myrank should not be owner of this crosslinker");
#endif

    // different treatment according to number of bonds crosslinker had before
    // this binding event
    switch( cldata_i.clnumbond )
    {
      case 0:
      {
        // update beam data
        // store crosslinker gid in status of beam binding spot
        beamdata_i.bbspotstatus[ binevdata.bspotlocn ] = binevdata.clgid;
        ele_i->SetBindingSpotStatus( beamdata_i.bbspotstatus );
        break;
      }
      case 1:
      {
        // update beam data
        // store crosslinker gid in status of beam binding spot
        beamdata_i.bbspotstatus[ binevdata.bspotlocn ] = binevdata.clgid;
        ele_i->SetBindingSpotStatus(beamdata_i.bbspotstatus );

        break;
      }
      default:
      {
        dserror("You should not be here, crosslinker has unrealistic number "
                "%i of bonds.", cldata_i.clnumbond );
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CreateNewDoubleBondedCrosslinkerElementPairs(
    std::map< int, NewDoubleBonds > const & mynewdbondcl)
{
  CheckInit();

  std::map< int, NewDoubleBonds >::const_iterator iter;
  for ( iter = mynewdbondcl.begin(); iter != mynewdbondcl.end(); ++iter )
  {
    NewDoubleBonds const& newdoublebond_i = iter->second;
    CROSSLINKING::CrosslinkerNode* cl_node = dynamic_cast<CROSSLINKING::CrosslinkerNode*>(
            BinDiscretPtr()->gNode(iter->first) );

    // create and initialize objects of beam-to-beam connections
    // Todo move this inside the create routines (or one create routine in BeamLink class)
    Teuchos::RCP<BEAMINTERACTION::BeamLink> linkelepairptr;
    if ( cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::beam3r_lin2_rigid )
      linkelepairptr = BEAMINTERACTION::BeamLinkRigidJointed::Create();
    else if ( cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::beam3r_lin2_pin
        or cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::truss )
      linkelepairptr = BEAMINTERACTION::BeamLinkPinJointed::Create( cl_node->GetMaterial()->JointType() );

    // finally initialize and setup object
    linkelepairptr->Init( iter->first, newdoublebond_i.eleids, newdoublebond_i.bspotposs, newdoublebond_i.bspottriads );
    linkelepairptr->Setup( cl_node->GetMaterial()->BeamElastHyperMatNum() );

    // add to my double bonds
    doublebondcl_[linkelepairptr->Id()] = linkelepairptr;

#ifdef DEBUG
    // safety check
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( BinDiscretPtr()->gNode( linkelepairptr->Id() ) );

    if ( crosslinker_i->ClData()->GetNumberOfBonds() != 2 )
     dserror("Setup: Cl with gid %i Owner %i on myrank %i and numbonds %i", linkelepairptr->Id(),
         crosslinker_i->Owner(), GStatePtr()->GetMyRank(), crosslinker_i->ClData()->GetNumberOfBonds() );
#endif
  }

  // print number of double bonded linker on each processor ( verbose)
  if ( mynewdbondcl.size() > 0 or doublebondcl_.size() > 0 )
  {
    IO::cout(IO::verbose) <<"************************************************"<<IO::endl;
    IO::cout(IO::verbose) << "PID " << GState().GetMyRank() << ": added " << mynewdbondcl.size()
        << " new db crosslinkers. Now have " << doublebondcl_.size() <<IO::endl;
  }

  // print total number of linker across all procs
  int num_local_linker = doublebondcl_.size();
  int num_global_linker = 0;
  MPI_Reduce( &num_local_linker, &num_global_linker, 1, MPI_INT, MPI_SUM, 0,
      dynamic_cast<const Epetra_MpiComm*>( &(Discret().Comm() ) )->Comm() );
  if ( GState().GetMyRank() == 0 )
  {
    IO::cout(IO::standard) << "\n************************************************\n" <<IO::endl;
    IO::cout(IO::standard) << "total number of double-bonded crosslinker: " << num_global_linker <<IO::endl;
    IO::cout(IO::standard) << "\n************************************************\n" <<IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    DoubleBindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids)
{
  CheckInit();

  for ( int nbonds = 0; nbonds < 2; ++nbonds )
  {
    // intended bonds of row crosslinker on myrank (key is clgid)
    std::map< int, BindEventData > mybonds;
    // intended bond col crosslinker to row element (key is owner of crosslinker != myrank)
    std::map< int, std::vector<BindEventData> > undecidedbonds;

    // store bins that have already been examined
    std::unordered_set<int> examinedbins;
    // this variable is used to check if a beam binding spot is linked twice on
    // myrank during a time step
    std::vector< std::set <int> > intendedbeambonds( DiscretPtr()->NumMyRowElements() );

    std::set<int>::const_iterator b_iter;
    for ( b_iter = bingids.begin(); b_iter != bingids.end(); ++b_iter )
    {
      // get neighboring bins whose bonds get dissolved
      std::vector<int> nb_binIds;
      nb_binIds.reserve(27);
      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds( *b_iter, nb_binIds );

      std::vector<int>::const_iterator nb_iter;
      for ( nb_iter = nb_binIds.begin(); nb_iter != nb_binIds.end(); ++nb_iter )
      {
        // check on existence of bin on this proc
        if ( not BinDiscretPtr()->HaveGlobalElement( *nb_iter ))
          continue;

        // if a bin has already been examined --> continue with next crosslinker
        if ( examinedbins.find( *nb_iter ) != examinedbins.end() )
          continue;
        //else: bin is examined for the first time --> new entry in examinedbins_
        else
          examinedbins.insert( *nb_iter );

        DRT::Element* currentbin =
            BinDiscretPtr()->gElement(*nb_iter);

        FindPotentialBindingEventsInBinAndNeighborhood( currentbin, mybonds,
            undecidedbonds, intendedbeambonds, false );
      }
    }

    // bind events where myrank only owns the elements, cl are taken care
    // of by their owner (key is clgid)
    std::map< int, BindEventData > myelebonds;

    // now each row owner of a linker gets requests, makes a random decision and
    // informs back its requesters
    ManageBindingInParallel( mybonds, undecidedbonds, myelebonds );

    // actual update of binding states is done here
    UpdateMyCrosslinkerAndElementBindingStates( mybonds, myelebonds );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnBindCrosslinker()
{
  CheckInit();

  DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );

  // data containing information about elements that need to be updated on
  // procs != myrank
  std::map< int, std::vector< UnBindEventData > > sendunbindevents;
  // elements that need to be updated on myrank
  std::vector< UnBindEventData > myrankunbindevents;

  // loop over all row linker (in random order) and dissolve bond if probability
  // criterion is met
  /* note: we loop over all row crosslinker, i.e. myrank needs to update all
   * crosslinker information. As it possible that a row crosslinker is linked
   * to col element, we potentially need to communicate if such an element
   * needs to be updated*/
  const int numrowcl = BinDiscretPtr()->NumMyRowNodes();
  std::vector<int> rorderrowcl = BEAMINTERACTION::UTILS::Permutation(numrowcl);
  std::vector<int>::const_iterator rowcli;
  for ( rowcli = rorderrowcl.begin(); rowcli != rorderrowcl.end(); ++rowcli )
  {
    CROSSLINKING::CrosslinkerNode* linker =
        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( BinDiscretPtr()->lRowNode( *rowcli ) );
    const int clcollid = linker->LID();
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // probability with which a crosslink breaks up in the current time step
    double p_unlink = 1.0 - exp( (-1.0) * crosslinking_params_ptr_->DeltaTime() * linker->GetMaterial()->KOff() );

    // different treatment according to number of bonds of a crosslinker
    switch( cldata_i.clnumbond )
    {
      case 0:
      {
        // nothing to do here
        break;
      }
      case 1:
      {
        // if probability criterion is not met, we are done here
        if ( DRT::Problem::Instance()->Random()->Uni() > p_unlink )
          break;

        // dissolve bond and update states
        DissolveBond( linker, GetSingleOccupiedClBspot( cldata_i.clbspots ), 1,
            sendunbindevents, myrankunbindevents );

        break;
      }
      case 2:
      {
        // calc unbind probability in case of force dependent off rate
        std::vector< double > p_unlink_db( 2, 0.0 );
        if ( abs( linker->GetMaterial()->DeltaBellEq() ) > 1.0e-8 )
          CalcBellsForceDependentUnbindProbability( linker, doublebondcl_[linker->Id()], p_unlink_db );
        else
          p_unlink_db[0] = p_unlink_db[1] = p_unlink;

        // loop through crosslinker bonds in random order
        std::vector<int> ro = BEAMINTERACTION::UTILS::Permutation( cldata_i.clnumbond );
        std::vector<int>::const_iterator clbspotiter;
        for ( clbspotiter = ro.begin(); clbspotiter != ro.end(); ++clbspotiter )
        {
          // if probability criterion isn't met, go to next spot
          if ( DRT::Problem::Instance()->Random()->Uni() > p_unlink_db[*clbspotiter] )
            continue;

          // dissolve bond and update states
          DissolveBond( linker, *clbspotiter, 2, sendunbindevents, myrankunbindevents );

          // we only want to dissolve one bond per timestep, therefore we go to
          // next crosslinker if we made it so far (i.e. a bond got dissolved)
          break;
        }

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i.clnumbond );
        exit(EXIT_FAILURE);
      }
    }
  }

  // communicate which elements need to be updated on rank != myrank
  CommunicateCrosslinkerUnbinding( sendunbindevents, myrankunbindevents );

  // update binding status of beam binding partners on myrank
  UpdateBeamBindingStatusAfterUnbinding( myrankunbindevents );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CalcBellsForceDependentUnbindProbability(
    CROSSLINKING::CrosslinkerNode* linker,
    Teuchos::RCP< BEAMINTERACTION::BeamLink > const& elepairptr,
    std::vector< double >& punlinkforcedependent) const
{
  CheckInitSetup();

  /* characteristic vond length (=nodereldis)
   * from B. Gui and W. Guilford: Mechanics of actomyosin bonds in different nucleotide states are tuned to muscle contraction
   * Fig 2: slip pathway, ADP, delta = 0.0004;
   * Note: delta < 0 -> catch bond, delta > 0 -> bond-weakening
   * see Kai Müller Dis p. 67/68 */
  double const delta = linker->GetMaterial()->DeltaBellEq();
  double const kt = crosslinking_params_ptr_->KT();
  double const koff = linker->GetMaterial()->KOff();

  // safety check
  if ( kt < 1e-08 )
    dserror(" Thermal energy (KT) set to zero, although you are about to divide by it. ");

  // force and moment exerted on the two binding sites of crosslinker with clgid
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );
  elepairptr->EvaluateForce( bspotforce[0], bspotforce[1] );

  // check if crosslinker is streched -> sgn+ or compressed -> sgn- by checking orientation of forcevector
  // note: this works only if there are no other forces (like inertia, stochastic, damping) acting on the cl
  LINALG::Matrix<3,1> dist_vec(true);
  LINALG::Matrix<3,1> bspotforceone(true);
  dist_vec.Update( 1.0, elepairptr->GetBindSpotPos1(), -1.0, elepairptr->GetBindSpotPos2() );
  for ( int j = 0; j < 3; ++j )
    bspotforceone(j) = bspotforce[0](j);
  double sgn = (dist_vec.Dot( bspotforceone ) < 0.0 ) ? -1.0 : 1.0;

  /* note: you have a sign criterion that is dependent on the axial strain, but the force you are
   * using also contains shear parts. This means in cases of axial strains near 0 and (large) shear
   * forces you are doing something strange (e.g. jumping between behaviour). Think about a two or
   * three dimensional calculation of the new koff, considering shear and axial forces independently */
//  if ( GState().GetMyRank() == 0 )
//    std::cout << "Warning: in cases of high shear forces and low axial strain you might "
//                 "be doing something strange using the force dependent off rate ..." << std::endl;

  // calculate new off rate
  std::vector< double > clbspotforcenorm ( 2, 0.0 );
  std::vector< double > forcedependentkoff ( 2, 0.0 );
  for ( unsigned int i = 0; i < 2; ++i )
  {
    // currently, only forces (not moments) considered
    bspotforce[i].Reshape(3,1);
    clbspotforcenorm[i] = bspotforce[i].Norm2();

    // adjusted off-rate according to Bell's equation (Howard, eq 5.10, p.89)
    if ( kt > 1e-16 )
      forcedependentkoff[i] = koff * exp( sgn * clbspotforcenorm[i] * delta / kt);

    // get respective force dependent unbind probability for each cl binding spot
    punlinkforcedependent[i] = 1.0 - exp( (-1.0) * crosslinking_params_ptr_->DeltaTime() * forcedependentkoff[i] );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateBeamBindingStatusAfterUnbinding(
    std::vector<UnBindEventData> const& unbindevent )
{
  CheckInit();

  // loop through all unbinding events on myrank
  std::vector<UnBindEventData>::const_iterator iter;
  for ( iter = unbindevent.begin(); iter != unbindevent.end(); ++ iter )
  {
    // get data
    const int elegidtoupdate = iter->eletoupdate.first;
    const int bspotlocn = iter->eletoupdate.second;
    const int colelelid = Discret().ElementColMap()->LID(elegidtoupdate);

#ifdef DEBUG
    // safety check
    if ( Discret().ElementRowMap()->LID(elegidtoupdate) < 0 )
      dserror("element with gid %i not owned by proc %i",elegidtoupdate,GState().GetMyRank());
#endif

    // get beam element of current binding event
    DRT::ELEMENTS::Beam3Base* ele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(Discret().lColElement(colelelid));

    BeamData& beamdata_i = beam_data_[colelelid];
    std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;

    // update beam data
    bbspotstatus_i[bspotlocn] = -1;
    ele_i->SetBindingSpotStatus(bbspotstatus_i);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyDoubleBondsAfterRedistribution()
{
  CheckInit();

  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink> > > dbcltosend;
  std::set<int> dbtoerase;

  // loop over all double bonds on myrank
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamLink> >::iterator iter;
  for ( iter = doublebondcl_.begin(); iter != doublebondcl_.end(); ++iter )
  {
    const int clgid = iter->first;

    // safety check
    if ( BinDiscretPtr()->NodeColMap()->LID(clgid) < 0 )
      dserror("Crosslinker %i moved further than the bin length in one time step on rank %i, "
              "this is not allowed (maybe increase cutoff radius). ", clgid, BeamInteractionDataStatePtr()->GetMyRank());

    DRT::Node* doublebondedcl_i = BinDiscretPtr()->gNode(clgid);

    // check ownership
    int owner = doublebondedcl_i->Owner();
    if ( owner != BeamInteractionDataStatePtr()->GetMyRank() )
    {
      if ( not doublebondcl_.count(clgid) )
        dserror("willing to delete double bond %i which is not existing", clgid);
      dbcltosend[owner].push_back( iter->second );
      dbtoerase.insert(clgid);
    }
  }

  std::set<int>::const_iterator i;
  for ( i = dbtoerase.begin(); i != dbtoerase.end(); ++i )
    doublebondcl_.erase(*i);

  // add new double bonds
  CommunicateBeamLinkAfterRedistribution(dbcltosend);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyDoubleBondsRemoteIdList()
{
  CheckInit();

  // loop over all double bonded crosslinker that were read during restart and
  // get the ones that do not belong to myrank
  std::set<int> notonmyrank;
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamLink> >::iterator iter;
  for ( iter = doublebondcl_.begin(); iter != doublebondcl_.end(); ++iter )
  {
    // double bonded crosslinker gid
    const int clgid = iter->first;

    // i) not even in myranks discret
    if ( BinDiscretPtr()->NodeColMap()->LID(clgid) < 0 )
    {
      notonmyrank.insert(clgid);
      continue;
    }

    // ii) myrank not owner
    if ( BinDiscretPtr()->gNode(clgid)->Owner() != BeamInteractionDataStatePtr()->GetMyRank() )
      notonmyrank.insert(clgid);
  }

  int const size = static_cast<int>( notonmyrank.size() );
  std::vector<int> unique_clgidlist( notonmyrank.begin(), notonmyrank.end() );
  std::vector<int> unique_pidlist(size);

  // find new host procs for double bonded crosslinker by communication
  int err = BinDiscretPtr()->NodeRowMap()->RemoteIDList( size, unique_clgidlist.data(), unique_pidlist.data(), NULL );
  if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d",err);

  std::map< int, std::vector< Teuchos::RCP<BEAMINTERACTION::BeamLink > > > dbcltosend;
  for ( unsigned int i = 0; i < static_cast<unsigned int>( unique_clgidlist.size() ); ++i )
    dbcltosend[ unique_pidlist[i] ].push_back( doublebondcl_[ unique_clgidlist[i] ] );

  // update myrank's map
  std::set< int >::const_iterator i;
  for ( i = notonmyrank.begin(); i != notonmyrank.end(); ++i )
    doublebondcl_.erase(*i);

  // send and receive double bonds
  CommunicateBeamLinkAfterRedistribution(dbcltosend);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    UnbindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids)
{
  CheckInit();

  UnbindCrosslinkerInBinsAndNeighborhood( bingids, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    UnbindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids,
        bool const doubleunbind)
{
  CheckInitSetup();

  std::set< int > binsonmyrank;
  DetermineResponsilbeProcsForForcedCrosslinkerUnbinding( bingids, binsonmyrank );

  // data containing information about elements that need to be updated on
  // procs != myrank
  std::map< int, std::vector< UnBindEventData > > sendunbindevents;
  // elements that need to be updated on myrank
  std::vector< UnBindEventData > myrankunbindevents;

  std::set<int>::const_iterator nb_iter;
  for ( nb_iter = binsonmyrank.begin(); nb_iter != binsonmyrank.end(); ++nb_iter )
  {
    // get all crosslinker in current bin
    DRT::Element* currentbin = BinDiscretPtr()->gElement( *nb_iter);
    DRT::Node **clincurrentbin = currentbin->Nodes();
    const int numcrosslinker = currentbin->NumNode();

    // loop over all crosslinker in current bin
    for ( int i = 0; i < numcrosslinker; ++i )
    {
      // get crosslinker in current bin
      DRT::Node *crosslinker_i = clincurrentbin[i];
      CrosslinkerData& cldata_i = crosslinker_data_[ crosslinker_i->LID() ];

#ifdef DEBUG
      // safety checks
      if ( crosslinker_i->Owner() != GState().GetMyRank() )
        dserror(" Only row owner of crosslinker changes its state, rank %i is not owner "
                "of linker with gid %i, but rank %i", GState().GetMyRank(), crosslinker_i->Id(), crosslinker_i->Owner() );
      if ( cldata_i.clnumbond != dynamic_cast< CROSSLINKING::CrosslinkerNode* >(crosslinker_i)->ClData()->GetNumberOfBonds() )
        dserror("Your crosslinker data container is not up to date, something went wrong ... ");
#endif

      switch( cldata_i.clnumbond )
      {
        case 0:
        {
          break;
        }
        case 1:
        {
          // dissolve bond and update states
          DissolveBond( crosslinker_i, GetSingleOccupiedClBspot( cldata_i.clbspots ),
              cldata_i.clnumbond, sendunbindevents, myrankunbindevents );
          break;
        }
        case 2:
        {
          // dissolve random bond and update states
          DissolveBond( crosslinker_i, BEAMINTERACTION::UTILS::Permutation( cldata_i.clnumbond )[0],
              cldata_i.clnumbond, sendunbindevents, myrankunbindevents );

          // in case we want to allow transition from double bonded to free, take same linker
          // again
          if ( doubleunbind )
            --i;
          break;
        }
        default:
        {
          dserror(" Unrealistic number %i of bonds for a crosslinker.", cldata_i.clnumbond );
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // communicate which elements need to be updated on rank != myrank
  CommunicateCrosslinkerUnbinding( sendunbindevents, myrankunbindevents );

  // update binding status of beam binding partners on myrank
  UpdateBeamBindingStatusAfterUnbinding( myrankunbindevents );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    DetermineResponsilbeProcsForForcedCrosslinkerUnbinding(
        std::set<int> const& bingids,
        std::set< int >& binsonmyrank) const
{
  CheckInitSetup();

  std::set< int > checkedbins;
  std::map< int, std::vector< int > > binstosend;

  // determine all bins that need to dissolve all bonds on myrank and
  // on other procs (these need to be informed by myrank)
  std::set<int>::const_iterator b_iter;
  for ( b_iter = bingids.begin(); b_iter != bingids.end(); ++b_iter )
  {
    std::vector<int> nb_binIds;
    nb_binIds.reserve(27);
    // do not check on existence here
    BinStrategy().GetNeighborAndOwnBinIds( *b_iter, nb_binIds );

    std::vector<int>::const_iterator nb_iter;
    for ( nb_iter = nb_binIds.begin(); nb_iter != nb_binIds.end(); ++nb_iter )
    {
      // safety check
      if ( not BinDiscret().HaveGlobalElement( *nb_iter ) )
        dserror("Not entire neighborhood ghosted, this is a problem in the following ");

      // if a bin has already been examined --> continue with next bin
      // like this we get a unique vector that myrank sends
      if ( checkedbins.find( *nb_iter ) != checkedbins.end() )
        continue;
      //else: bin is examined for the first time --> new entry in examinedbins_
      else
        checkedbins.insert( *nb_iter );

      // decide who needs to dissolve bonds
      const int owner = BinDiscret().gElement( *nb_iter )->Owner();
      if ( owner == GState().GetMyRank() )
        binsonmyrank.insert(*nb_iter);
      else
        binstosend[owner].push_back(*nb_iter);
    }
  }

  CommunicateBinIds( binstosend, binsonmyrank );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateBinIds(
    std::map< int, std::vector< int > > const& binstosend,
    std::set< int >& binsonmyrank) const
{
  CheckInitSetup();

  // build exporter
  DRT::Exporter exporter( Discret().Comm() );
  int const numproc = Discret().Comm().NumProc();
  int const myrank = GState().GetMyRank();

  // ---- send ---- ( we do not need to pack anything)
  int const length = binstosend.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  std::map<int, std::vector<int> >::const_iterator p;
  std::vector<int> targetprocs( numproc, 0 );
  for ( p = binstosend.begin(); p != binstosend.end(); ++p )
  {
    targetprocs[p->first] = 1;
    exporter.ISend( myrank, p->first, &((p->second)[0]), static_cast<int>( (p->second).size() ), 1234, request[tag] );
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // ---- prepare receiving procs -----
  std::vector<int> summedtargets( numproc, 0) ;
  Discret().Comm().SumAll( targetprocs.data(), summedtargets.data(), numproc );

  // ---- receive ----- (we do not need to unpack anything)
  for ( int rec = 0; rec < summedtargets[myrank]; ++rec )
  {
    std::vector<int> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny( from, tag ,rdata, length );
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank, from);

    // insert in binsonmyrank
    binsonmyrank.insert( rdata.begin(), rdata.end() );
  }

  // wait for all communication to finish
  Wait( exporter, request, static_cast<int>( binstosend.size() ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DissolveBond(
    DRT::Node* linker,
    int const freedbspotid,
    int const numbondsold,
    std::map< int, std::vector< UnBindEventData > >& sendunbindevents,
    std::vector< UnBindEventData >& myrankunbindevents)
{
  CheckInitSetup();

#ifdef DEBUG
    // safety check
    if ( numbondsold < 1 )
      dserror("dissolution of free crosslinker does not make any sense" );
#endif

  CROSSLINKING::CrosslinkerNode* crosslinker =
      dynamic_cast< CROSSLINKING::CrosslinkerNode* >(linker);

  // get linker data
  const int clcollid = crosslinker->LID();
  CrosslinkerData& cldata = crosslinker_data_[clcollid];

  // store unbinding event data
  UnBindEventData unbindevent;
  unbindevent.eletoupdate = cldata.clbspots[freedbspotid];

  // owner of beam
  const int beamowner =
      DiscretPtr()->gElement(unbindevent.eletoupdate.first)->Owner();

  // check who needs to update the element status
  if ( beamowner == GState().GetMyRank() )
    myrankunbindevents.push_back(unbindevent);
  else
    sendunbindevents[beamowner].push_back(unbindevent);

  // -----------------------------------------------------------------
  // update crosslinker status
  // -----------------------------------------------------------------
  // update binding status of linker
  cldata.clbspots[freedbspotid].first = -1;
  cldata.clbspots[freedbspotid].second = -1;
  crosslinker->ClData()->SetClBSpotStatus(cldata.clbspots);

  // update number of bonds
  cldata.clnumbond = numbondsold - 1;
  crosslinker->ClData()->SetNumberOfBonds(cldata.clnumbond);

  if ( numbondsold == 1 )
  {
    SetPositionOfNewlyFreeCrosslinker( crosslinker, cldata.clpos );
  }
  else if ( numbondsold == 2 )
  {
    int stayoccpotid = 0;
    if ( freedbspotid == 0 )
      stayoccpotid = 1;

    SetPositionOfNewlySingleBondedCrosslinker( linker, cldata, stayoccpotid );

    // safety check
    if ( not doublebondcl_.count( linker->Id() ) )
      dserror("crosslinker %i with %i bonds is not in double bonded map of rank %i",
             linker->Id(), crosslinker->ClData()->GetNumberOfBonds() + 1, GStatePtr()->GetMyRank() );

    // erase crosslinker from double bonded crosslinker list
    doublebondcl_.erase( linker->Id() );
  }
  else
  {
    dserror("dissolution of free linker does not make any sense, something went wrong." );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateUndecidedBonds(
    DRT::Exporter& exporter,
    std::map<int, std::vector<BindEventData> >& undecidedbonds,
    int& numrecrequest,
    std::map<int, std::vector<BindEventData> >& requestedcl) const
{
  CheckInit();

  // -----------------------------------------------------------------------
  // unblocking send
  // -----------------------------------------------------------------------
  std::vector<MPI_Request> request;
  ISend( exporter, request, undecidedbonds );

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs( undecidedbonds, summedtargets );

  numrecrequest = summedtargets[GState().GetMyRank()];
  for ( int rec = 0; rec < numrecrequest; ++rec )
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny( from, tag, rdata, length );
    if (tag != 1234)
     dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from );

    // store received data
    std::vector<char>::size_type position = 0;
    while ( position < rdata.size() )
    {
      // ---- extract received data -----
      BindEventData reccldata;
      UnPack( position, rdata, reccldata );

      // create map holding all requests
      requestedcl[reccldata.clgid].push_back(reccldata);
    }

    if ( position != rdata.size() )
      dserror("Mismatch in size of data %d <-> %d", static_cast<int>( rdata.size() ), position );
  }

  // wait for all communication to finish
  Wait( exporter, request, static_cast<int>( undecidedbonds.size() ) );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateDecidedBonds(
    DRT::Exporter& exporter,
    std::map<int, std::vector<BindEventData> >& decidedbonds,
    std::map<int, BindEventData >&              myelebonds,
    int const numrecrequest,
    int const answersize) const
{
  CheckInit();

  // -----------------------------------------------------------------------
  // send back decisions for all requests that were made
  // -----------------------------------------------------------------------
  // store requested cl and its data
  std::vector<MPI_Request> request;

#ifdef DEBUG
  // safety check
  if (static_cast<int>(decidedbonds.size()) != numrecrequest)
    dserror("Number of received requests %i unequal to number of answers %i",
        numrecrequest,decidedbonds.size());
#endif

  // unblocking send
  ISend( exporter, request, decidedbonds );

#ifdef DEBUG
  std::vector<int> summedtargets;
  PrepareReceivingProcs(decidedbonds,summedtargets);
  if ( answersize != summedtargets[GState().GetMyRank()] )
    dserror(" proc %i did not get an answer to all its questions, that it not fair.");
#endif

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // store requested cl and its data
  for ( int rec = 0; rec < answersize; ++rec )
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from);

    // store received data
    std::vector<char>::size_type position = 0;
    while ( position < rdata.size() )
    {
      // ---- extract received data -----
      BindEventData reccldata;
      UnPack(position,rdata,reccldata);

      // add binding events to new colbond map
      if (reccldata.permission)
        myelebonds[reccldata.clgid] = reccldata;

    }

    if ( position != rdata.size() )
      dserror("Mismatch in size of data %d <-> %d", static_cast<int>( rdata.size() ), position );
  }

  // wait for all communication to finish
  Wait( exporter, request, static_cast< int >( decidedbonds.size() ) );

}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DecideBindingInParallel(
    std::map<int, std::vector<BindEventData> >& requestedcl,
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& decidedbonds) const
{
  CheckInit();

  std::map<int, std::vector<BindEventData> >::iterator cliter;
  // loop over all requested cl (note myrank is owner of these)
  for ( cliter = requestedcl.begin(); cliter != requestedcl.end(); ++cliter )
  {
    // check if myrank wants to bind this crosslinker
    bool myrankbond = false;
    if ( mybonds.find(cliter->first) != mybonds.end() )
      myrankbond = true;

    // ---------------------------------------------------------------------
    // if only one request and myrank does not want to bind this cl,
    // requesting proc gets the permission to do so
    // ---------------------------------------------------------------------
    if ( static_cast<int>( cliter->second.size() ) == 1 and not myrankbond )
    {
      // we send back the permission to the relevant proc, because myrank as row
      // owner of bspot needs to set the respective stuff for the element of this
      // binding event
      // note: permission = true was send as default, so this can be sent back
      // without changes
      decidedbonds[cliter->second[0].requestproc].push_back(cliter->second[0]);

#ifdef DEBUG
      if ( cliter->second[0].permission != 1 )
        dserror(" something during communication went wrong, default true permission "
                " not received");
#endif

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
    int numrequprocs = static_cast<int>( cliter->second.size() );
    if (myrankbond)
      numrequprocs += 1;

    // get random proc out of affected ones
    DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );
    // fixme: what if random number exactly = 1?
    int rankwithpermission = std::floor( numrequprocs * DRT::Problem::Instance()->Random()->Uni() );

    // myrank is allowed to set link
    if ( myrankbond and rankwithpermission == (numrequprocs - 1) )
    {
      // note: this means link is set between row cl and row ele on myrank,
      // all relevant information for myrank is stored in mybonds
      // loop over all requesters and store their veto
      std::vector<BindEventData>::iterator iter;
      for ( iter = cliter->second.begin(); iter != cliter->second.end(); ++iter )
      {
        iter->permission = 0;
        decidedbonds[iter->requestproc].push_back(*iter);
      }
    }
    // certain requester is allowed to set the link
    else
    {
      // loop over all requesters and store veto for all requester except for one
      std::vector<BindEventData>::iterator iter;

      int counter = 0;
      for ( iter = cliter->second.begin(); iter != cliter->second.end(); ++iter )
      {
        if ( rankwithpermission == counter )
        {
          // permission for this random proc
          decidedbonds[iter->requestproc].push_back(*iter);

#ifdef DEBUG
        if ( iter->permission != 1 )
          dserror(" something during communication went wrong, default true permission "
                  " not received");
#endif

          // erase old binding event
          if ( myrankbond )
            mybonds.erase( cliter->first );

          // insert new binding event
          mybonds[cliter->first] = *iter;
        }
        else
        {
          iter->permission = 0;
          decidedbonds[iter->requestproc].push_back(*iter);
        }
        counter++;
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateBeamLinkAfterRedistribution(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink> > >& dbondcltosend)
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter( DiscretPtr()->Comm() );
  int const numproc = DiscretPtr()->Comm().NumProc();

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char> > sdata;
  std::vector<int> targetprocs( numproc, 0 );
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink> > >::const_iterator p;
  for ( p = dbondcltosend.begin(); p != dbondcltosend.end(); ++p )
  {
    std::vector< Teuchos::RCP< BEAMINTERACTION::BeamLink > >::const_iterator iter;
    for ( iter = p->second.begin(); iter != p->second.end(); ++iter )
    {
     DRT::PackBuffer data;
     (*iter)->Pack(data);
     data.StartPacking();
     (*iter)->Pack(data);
     sdata[p->first].insert( sdata[p->first].end(), data().begin(), data().end() );
    }
    targetprocs[p->first] = 1;
  }

  // ---- send ----
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for (std::map<int, std::vector<char> >::const_iterator p = sdata.begin(); p != sdata.end(); ++p )
  {
    exporter.ISend( GState().GetMyRank(), p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets( numproc, 0) ;
  DiscretPtr()->Comm().SumAll( targetprocs.data(), summedtargets.data(), numproc );

  // myrank receive all packs that are sent to him
  for ( int rec = 0; rec < summedtargets[GState().GetMyRank()]; ++rec )
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny( from, tag ,rdata, length );
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from);

    // store received data
    std::vector<char>::size_type position = 0;
    while ( position < rdata.size() )
    {
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack( position, rdata, data );
      // this Teuchos::rcp holds the memory
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp( DRT::UTILS::Factory(data), true );
      Teuchos::RCP<BEAMINTERACTION::BeamLink> beamtobeamlink =
          Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLink>(object);
      if ( beamtobeamlink == Teuchos::null )
        dserror("Received object is not a beam to beam linkage");

#ifdef DEBUG
      // some safety checks
      if ( BinDiscretPtr()->gNode( beamtobeamlink->Id() )->Owner() != GStatePtr()->GetMyRank() )
        dserror(" A double bond was sent to rank %i, although it is not the owner of "
                "the cl with gid %i ", GStatePtr()->GetMyRank(), beamtobeamlink->Id() );
      if ( doublebondcl_.count(beamtobeamlink->Id()) )
        dserror(" Rank %i got sent double bonded crosslinker %i which it already has ",  GStatePtr()->GetMyRank(), beamtobeamlink->Id() );
#endif

      // insert new double bonds in my list
      doublebondcl_[beamtobeamlink->Id()] = beamtobeamlink;
    }

    if ( position != rdata.size() )
      dserror("Mismatch in size of data %d <-> %d",static_cast<int>( rdata.size() ),position );
  }

  // wait for all communication to finish
  Wait( exporter, request, static_cast<int>( sdata.size() ) );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateCrosslinkerUnbinding(
    std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
    std::vector<UnBindEventData>&                 myrankunbindevent) const
{
  CheckInit();

  ISendRecvAny( sendunbindevent, myrankunbindevent );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter& exporter,
    std::vector<MPI_Request>& request,
    std::map<int, std::vector<T> > const& send) const
{
  CheckInit();

  // ---- pack data for sending -----
  std::map<int, std::vector<char> > sdata;
  typename std::map<int, std::vector<T> >::const_iterator p;
  for ( p = send.begin(); p != send.end(); ++p )
  {
    typename std::vector<T>::const_iterator iter;
    DRT::PackBuffer data;
    for ( iter = p->second.begin(); iter != p->second.end(); ++iter )
    {
     Pack(data,*iter);
    }
    data.StartPacking();
    for ( iter = p->second.begin(); iter != p->second.end(); ++iter )
    {
     Pack(data,*iter);
    }
    sdata[p->first].insert( sdata[p->first].end(), data().begin(), data().end() );
  }

  // ---- send ----
  const int length = sdata.size();
  request.resize(length);
  int tag = 0;
  for (std::map<int, std::vector<char> >::const_iterator p = sdata.begin(); p != sdata.end(); ++p )
  {
    exporter.ISend( GState().GetMyRank(), p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<T> > const& datasenttorank,
    std::vector<int>& summedtargets) const
{
  CheckInit();

  const int numproc = Discret().Comm().NumProc();

  // get number of procs from which myrank receives data
  std::vector<int> targetprocs(numproc,0);
  typename std::map<int, std::vector<T> >::const_iterator prociter;
  for (prociter=datasenttorank.begin(); prociter!=datasenttorank.end(); ++prociter)
    targetprocs[prociter->first] = 1;
  // store number of messages myrank receives
  summedtargets.resize(numproc,0);
  BinDiscret().Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&  exporter,
    int const receivesize,
    std::vector<T>& recv) const
{
  CheckInit();

  // myrank receive all packs that are sent to him
  for ( int rec = 0; rec < receivesize; ++rec )
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny( from, tag, rdata, length );
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from);

    // store received data
    std::vector<char>::size_type position = 0;
    while ( position < rdata.size() )
    {
      // ---- extract received data -----
      T recdata;
      UnPack( position, rdata, recdata );

      // add received data to list of unbindevents on myrank
      recv.push_back(recdata);
    }

    if ( position != rdata.size() )
      dserror( "Mismatch in size of data %d <-> %d", static_cast<int>( rdata.size() ), position );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
  const std::map<int, std::vector<T> >& send,
  std::vector<T>&                       recv) const
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter( BinDiscret().Comm() );

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // unblocking send
  std::vector<MPI_Request> request;
  ISend( exporter, request,send );

  // -----------------------------------------------------------------------
  // prepare receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs( send, summedtargets );

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  int receivesize = summedtargets[GState().GetMyRank()];
  RecvAny( exporter, receivesize, recv );

  // wait for all communication to finish
  Wait( exporter, request, static_cast<int>(send.size() ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Wait(
    DRT::Exporter& exporter,
    std::vector<MPI_Request>& request,
    int const length) const
{
  CheckInit();

  // wait for all communication to finish
  for ( int i = 0; i < length; ++i )
    exporter.Wait( request[i] );

  // note: if we have done everything correct, this should be a no time operation
  BinDiscret().Comm().Barrier(); // I feel better this way ;-)
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&     data,
  BindEventData const& bindeventdata) const
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,bindeventdata.clgid);
  DRT::ParObject::AddtoPack(data,bindeventdata.elegid);
  DRT::ParObject::AddtoPack(data,bindeventdata.bspotlocn);
  DRT::ParObject::AddtoPack(data,bindeventdata.requestproc);
  DRT::ParObject::AddtoPack(data,bindeventdata.permission);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&       data,
  UnBindEventData const& unbindeventdata) const
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,unbindeventdata.eletoupdate);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  BindEventData&                bindeventdata) const
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.clgid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.elegid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.bspotlocn);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.requestproc);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.permission);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  UnBindEventData&              unbindeventdata) const
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,unbindeventdata.eletoupdate);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrintAndCheckBindEventData(
  BindEventData const& bindeventdata) const
{
  CheckInit();

  // extract data
  std::cout << "\n Rank: " << GState().GetMyRank() << std::endl;
  std::cout << " crosslinker gid " << bindeventdata.clgid << std::endl;
  std::cout << " element gid " << bindeventdata.elegid << std::endl;
  std::cout << " bspot local number " << bindeventdata.bspotlocn << std::endl;
  std::cout << " requesting proc " << bindeventdata.requestproc << std::endl;
  std::cout << " permission " << bindeventdata.permission << std::endl;

  if ( bindeventdata.clgid < 0 or bindeventdata.elegid < 0 or bindeventdata.bspotlocn < 0
      or bindeventdata.requestproc < 0 or not (bindeventdata.permission == 0 or bindeventdata.permission == 1 ) )
    dserror(" your bindevent does not make sense.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------
void STR::MODELEVALUATOR::BeamInteraction::AddCrosslinkerToBinDiscret()
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
  if ( myrank_ == 0 )
  {
    for ( int i = 0; i < numcrosslinker; ++i )
    {
      // random reference position of crosslinker in bounding box
      std::vector<double> X(3);
      for ( int dim = 0; dim < 3; ++dim )
      {
        double min = 0;
        double edgelength = 15;
//        if (dim == 0)
//        {
//          min = 5.0;
//          edgelength = 15;
//        }
        if (dim == 2)
        {
          min = 7.0;
          edgelength = 1;
        }
        X[dim] = min + ( edgelength * randpos[ i + dim ] );
      }

      Teuchos::RCP<DRT::Node> newcrosslinker =
          Teuchos::rcp(new CROSSLINKING::CrosslinkerNode( i, &X[0], myrank_ ) );

      // todo: put next two calls in constructor of CrosslinkerNode?
      // init crosslinker data container
      Teuchos::RCP< CROSSLINKING::CrosslinkerNode > clnode =
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
}*/


//-----------------------------------------------------------------------------
// explicit template instantiation (to please every compiler)
//-----------------------------------------------------------------------------
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter&,std::vector<MPI_Request>&, std::map<int, std::vector<BindEventData> > const& ) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter&,std::vector<MPI_Request>&, std::map<int, std::vector<UnBindEventData> > const& ) const;

template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<BindEventData> > const&, std::vector<int>& ) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<UnBindEventData> > const&, std::vector<int>& ) const;

template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&, int const, std::vector<BindEventData>& ) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&, int const, std::vector<UnBindEventData>& ) const;

template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<BindEventData> > const&, std::vector<BindEventData>& ) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<UnBindEventData> > const&,std::vector<UnBindEventData>&) const;
