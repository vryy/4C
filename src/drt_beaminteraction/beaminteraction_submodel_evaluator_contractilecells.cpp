/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_contractilecells.cpp

\brief class for managing contractile cell to network interaction

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_contractilecells.H"
#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_crosslinking.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"
#include "../drt_beaminteraction/contractilecells_params.H"
#include "../drt_beaminteraction/crosslinker_node.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_beaminteraction/beam_link_pinjointed.H"
#include "../drt_beaminteraction/beam_link_beam3r_lin2_pinjointed.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_particle/particle_handler.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include <Teuchos_TimeMonitor.hpp>



/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::ContractileCells() :
    sm_crosslinkink_ptr ( Teuchos::null ),
    contractilecells_params_ptr_ ( Teuchos::null )
{
  cell_beam_element_pairs_.clear();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::Setup()
{
  CheckInit();

  // construct, init and setup data container for crosslinking
  contractilecells_params_ptr_ = Teuchos::rcp( new BEAMINTERACTION::ContractileCellsParams() );
  contractilecells_params_ptr_->Init();
  contractilecells_params_ptr_->Setup();

  // set flag
  issetup_ = true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::PostSetup()
{
  CheckInitSetup();
 // nothing to do (yet)
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::InitSubmodelDependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap)
{
  CheckInitSetup();

  // init pointer to crosslinker submodel
  STR::MODELEVALUATOR::BeamInteraction::Map::const_iterator miter;
  for ( miter = (*submodelmap).begin(); miter != (*submodelmap).end(); ++miter )
    if ( miter->first == INPAR::BEAMINTERACTION::submodel_crosslinking )
      sm_crosslinkink_ptr = Teuchos::rcp_dynamic_cast<BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking>(miter->second);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::Reset()
{
  CheckInitSetup();

  // reset crosslinker pairs
  for ( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = iter.second;

    // get elements
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elepairptr->GetEleGid(0) ) );
    DRT::ELEMENTS::Beam3Base const * beamele =
        dynamic_cast< DRT::ELEMENTS::Beam3Base const* >( Discret().gElement( elepairptr->GetEleGid(1) ) );

    // init position of linker nodes
    std::vector< LINALG::Matrix< 3, 1 > > pos( 2, LINALG::Matrix< 3, 1>(true) );

    // sphere current position
    std::vector<double> sphereeledisp;
    BEAMINTERACTION::UTILS::GetCurrentElementDis( Discret(), sphere,
        BeamInteractionDataState().GetDisColNp(), sphereeledisp);

    // note: sphere has just one node (with three translational dofs)
    for ( unsigned int dim = 0; dim < 3; ++dim )
      pos[0](dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

    // beam bspot pos
    std::vector<double> beameledisp;
    BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), beamele,
        BeamInteractionDataState().GetDisColNp(), PeriodicBoundingBox(), beameledisp );
    beamele->GetPosOfBindingSpot( pos[1], beameledisp, elepairptr->GetLocBSpotNum(1), PeriodicBoundingBox() );

    // unshift one of the positions if both are separated by a periodic boundary
    // condition, i.e. have been shifted before
    PeriodicBoundingBoxPtr()->UnShift3D( pos[1], pos[0] );

    // finally reset state
    elepairptr->ResetState( pos );
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateForce()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );

  // resulting discrete element force vectors of the two parent elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  std::vector< std::vector<LINALG::SerialDenseMatrix> > dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for ( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > elepairptr = iter.second;

    for ( unsigned int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();
    }

    // evaluate beam linkage object to get forces of binding spots
    elepairptr->EvaluateForce( bspotforce[0], bspotforce[1] );

    // apply forces on binding spots to parent elements
    // and get their discrete element force vectors
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceToParentElements< DRT::ELEMENTS::Rigidsphere, DRT::ELEMENTS::Beam3Base >(
        Discret(), PeriodicBoundingBoxPtr(), BeamInteractionDataStatePtr()->GetMutableDisColNp(),
        elepairptr, bspotforce, eleforce );

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
        eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null);
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateStiff()
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

  for ( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = iter.second;

    for ( unsigned int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);

      for ( unsigned int j = 0; j< 2; ++j )
        bspotstiff[i][j].Zero();
    }

     // evaluate beam linkage object to get linearizations of forces on binding spots
    elepairptr->EvaluateStiff( bspotstiff[0][0], bspotstiff[0][1], bspotstiff[1][0],
        bspotstiff[1][1] );

    // apply linearizations to parent elements and get their discrete element stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotStiffToParentElements< DRT::ELEMENTS::Rigidsphere, DRT::ELEMENTS::Beam3Base >(
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

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateForceStiff()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector< std::vector<LINALG::SerialDenseMatrix> > bspotstiff( 2,
      std::vector<LINALG::SerialDenseMatrix>( 2, LINALG::SerialDenseMatrix( 6, 6 ) ) );

  // resulting discrete element force vectors of the two parent elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector< std::vector<LINALG::SerialDenseMatrix> > elestiff( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for ( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = iter.second;
    for ( unsigned int i = 0; i < 2; ++i )
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();

      for ( int j = 0; j< 2; ++j )
        bspotstiff[i][j].Zero();
    }

    // evaluate beam linkage object to get forces on binding spots
    elepairptr->EvaluateForceStiff( bspotforce[0], bspotforce[1], bspotstiff[0][0],
        bspotstiff[0][1], bspotstiff[1][0], bspotstiff[1][1] );

    // apply forces on binding spots and corresponding linearizations to parent elements
    // and get their discrete element force vectors and stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceStiffToParentElements< DRT::ELEMENTS::Rigidsphere, DRT::ELEMENTS::Beam3Base >(
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

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

}
/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::PreUpdateStepElement()
{
  CheckInitSetup();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::UpdateStepElement()
{
  CheckInitSetup();
  AddNewElementPairsAfterFeasibilityCheck();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::PostUpdateStepElement()
{
  CheckInitSetup();

  // print total number of linker across all procs
  int num_local_linker = cell_beam_element_pairs_.size();
  int num_global_linker = 0;
  MPI_Reduce( &num_local_linker, &num_global_linker, 1, MPI_INT, MPI_SUM, 0,
      dynamic_cast<const Epetra_MpiComm*>( &(Discret().Comm() ) )->Comm() );
  if ( GState().GetMyRank() == 0 )
  {
    IO::cout(IO::standard) << "\n************************************************\n" <<IO::endl;
    IO::cout(IO::standard) << "total number of Integrins: " << num_global_linker <<IO::endl;
    IO::cout(IO::standard) << "\n************************************************\n" <<IO::endl;
  }

//  std::set<int> spherebingids;
//  int spheregid = -1;
//  for( int i = 0; i < EleTypeMapExtractorPtr()->SphereMap()->NumMyElements(); ++i )
//  {
//    spheregid =  EleTypeMapExtractorPtr()->SphereMap()->GID(i);
//    spherebingids.insert( BeamInteractionDataStatePtr()->GetRowEleToBinSet(spheregid).begin(),
//                          BeamInteractionDataStatePtr()->GetRowEleToBinSet(spheregid).end() );
//  }
//
//  sm_crosslinkink_ptr->UnbindCrosslinkerInBinsAndNeighborhood( spherebingids, false );

//  int const updateevery = 200;
//
//  if( (GState().GetStepN() + 1) % updateevery == 0 and GState().GetStepN() > updateevery)
//    sm_crosslinkink_ptr->DoubleBindCrosslinkerInBinsAndNeighborhood( spherebingids );
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::RuntimeOutputStepState() const
{
  CheckInitSetup();

}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();

  // -------------------------------------------------------------------------
  // write list of cell beam pinjointed links on each proc
  // -------------------------------------------------------------------------
  DRT::PackBuffer buffer;
  for( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > btbl = iter.second;
    btbl->Pack(buffer);
  }
  buffer.StartPacking();
  for( auto const & iter : cell_beam_element_pairs_ )
  {
    Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > btbl = iter.second;
    btbl->Pack(buffer);
  }

  Teuchos::RCP<std::vector<char> > block = Teuchos::rcp( new std::vector<char> );
  std::swap( *block, buffer() );

  // write data
  BinDiscret().Writer()->WriteCharVector( "Integrins", block );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  IO::DiscretizationReader reader( BinDiscretPtr(), GState().GetStepN() );

  // -------------------------------------------------------------------------
  // read list of cell beam pinjointed links on each proc
  // -------------------------------------------------------------------------
  Teuchos::RCP< std::vector< char > > charvec;
  reader.ReadCharVector( charvec, "Itegrins" );

  std::vector<char>::size_type index = 0;
  while ( index < charvec->size() )
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack( index, *charvec, data );
    Teuchos::RCP< DRT::ParObject > object = Teuchos::rcp( DRT::UTILS::Factory(data), true );
    Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > beamlinkpinjointed =
        Teuchos::rcp_dynamic_cast< BEAMINTERACTION::BeamLinkPinJointed >(object);
    if ( beamlinkpinjointed == Teuchos::null )
      dserror("Failed of dynamic cast from ParObject to BeamLinkPinJointed");

    // insert in my list cell to beam joints
    cell_beam_element_pairs_[beamlinkpinjointed->Id()] = beamlinkpinjointed;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::PostReadRestart()
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::AddBinsToBinColMap(
    std::set< int >& colbins)
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::
    AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins ) const
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::AddNewElementPairsAfterFeasibilityCheck()
{
  CheckInitSetup();

  std::map< int, std::vector< std::pair< int, int > > > newlinks;
  FindAndStoreNeighboringElements( newlinks );
  CreateBeamToSphereJoint( newlinks );
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::FindAndStoreNeighboringElements(
    std::map< int, std::vector< std::pair< int, int > > > & newlinks )
{
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::"
                            "ContractileCells::FindAndStoreNeighboringElements");

  CheckInitSetup();

  std::unordered_set< int > tobebonded;

  // loop over all sphere elements
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::Element* currsphere = DiscretPtr()->gElement(elegid);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all bins touched by currele
    for ( auto const & biniter : BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid) )
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds( biniter, loc_neighboring_binIds );

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert( loc_neighboring_binIds.begin(), loc_neighboring_binIds.end() );
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds( neighboring_binIds.begin(), neighboring_binIds.end() );

    // set of beam elements that reside in neighboring bins
    std::set< DRT::Element * > neighboring_elements;
    std::vector< BINSTRATEGY::UTILS::BinContentType > bc( 1, BINSTRATEGY::UTILS::Beam );
    BinStrategyPtr()->GetBinContent( neighboring_elements, bc, glob_neighboring_binIds );

    // sort out elements that do not meet bind event criteria
    CheckFeasibilityOfNewLink( currsphere, neighboring_elements, tobebonded, newlinks );

  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::CheckFeasibilityOfNewLink(
    DRT::Element const * currele,
    std::set< DRT::Element * > const & neighbors,
    std::unordered_set< int > & tobebonded,
    std::map< int, std::vector< std::pair< int, int > > > & newlinks ) const
{
  CheckInitSetup();

  DRT::ELEMENTS::Rigidsphere const * sphere =
      dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >(currele);

  // current sphere position
  std::vector<double> sphereeledisp;
  BEAMINTERACTION::UTILS::GetCurrentElementDis( Discret(), currele,
      BeamInteractionDataState().GetDisColNp(), sphereeledisp);

  // note: sphere has just one node (with three translational dofs)
  LINALG::Matrix< 3, 1 > spherepos(true);
  for ( unsigned int dim = 0; dim < 3; ++dim )
    spherepos(dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

  // loop over all neighboring elements
  for( auto const & eiter : neighbors )
  {
    DRT::ELEMENTS::Beam3Base const * beamele =
        dynamic_cast< DRT::ELEMENTS::Beam3Base const* >(eiter);

#ifdef DEBUG
    if ( sphere == NULL or beamele == NULL )
      dserror(" First element should be a sphere, second element a beam.");
#endif

    std::vector<double> beameledisp;
    BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), beamele,
        BeamInteractionDataState().GetDisColNp(), PeriodicBoundingBox(), beameledisp );

    LINALG::Matrix< 3, 1 > bspotpos(true);
    LINALG::Matrix< 3, 3 > bspottriad(true);

    // loop over binding spots of neighboring element
    for ( unsigned int bspot_i = 0; bspot_i < beamele->GetNumberOfBindingSpots(); ++bspot_i )
    {

      std::pair< int, int > bspotpair = std::make_pair( beamele->Id(), bspot_i );

      // 1. criterion: bond status of sphere
      if ( sphere->DoesBondExist( bspotpair ) )
        continue;

      // 2. criterion: is beam binding spot free (covered by criterion 1? Only in case of separate linker
      // for cell to beam and beam to beam (todo)


      // 3. criterion: distance
      // get current position at binding spot xi
      bspotpos.Clear();
      beamele->GetPosOfBindingSpot( bspotpos, beameledisp, bspot_i, PeriodicBoundingBox() );

      double const linkdistmin = sphere->Radius() - ( sphere->Radius() / 10 );
      double const linkdistmax = sphere->Radius() + ( sphere->Radius() / 10 );


      if ( BEAMINTERACTION::UTILS::IsDistanceOutOfRange( spherepos, bspotpos, linkdistmin, linkdistmax ) )
        continue;

//        //DEBUG
//      std::cout << "sphere gid " << sphere->Id() << std::endl;
//      std::cout << " beam gid " << beamele->Id() << std::endl;
//      std::cout << " loc " << bspot_i << std::endl;
//
//      spherepos.Print(std::cout);
//      bspotpos.Print(std::cout);

//      // 4. criterion: orientation
//      // get current triad at binding spot xi
//      bspottriad.Clear();
//      beamele->GetTriadOfBindingSpot( bspottriad, beameledisp, bspot_i );
//
//      // note: we use first base vector instead of tangent vector here
//      LINALG::Matrix< 3, 1> curr_bindingspot_beam_tangent(true);
//      for ( unsigned int idim = 0; idim < 3; ++idim )
//        curr_bindingspot_beam_tangent(idim) = bspottriad( idim, 0 );
//
//      if ( BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
//          occ_bindingspot_beam_tangent, curr_bindingspot_beam_tangent, linkanglemin, linkanglemax ) )
//      {
//        neighbors.erase( eiter++ );
//        continue;
//      }

      // 5. criterion:

      // check if bspot will be bonded this step
      int bpspotgid = BEAMINTERACTION::UTILS::CantorPairing(bspotpair);

      if ( tobebonded.find(bpspotgid) != tobebonded.end() )
        continue;
      else
        tobebonded.insert(bpspotgid);

      // add new link
      newlinks[sphere->Id()].push_back( bspotpair );


    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::CreateBeamToSphereJoint(
    std::map< int, std::vector< std::pair< int, int > > > const & newlinks )
{
  CheckInitSetup();

  // loop over cells
  for ( auto const & newlinkiter : newlinks )
  {
    // init position of linker nodes
    std::vector< LINALG::Matrix< 3, 1 > > pos( 2, LINALG::Matrix< 3, 1>(true) );

    int const spheregid = newlinkiter.first;
    // get elements
    DRT::ELEMENTS::Rigidsphere * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere * >( Discret().gElement(spheregid) );

    std::vector< std::pair< int, int > > eleids (2);
    // todo: for now, sphere has one binding spot, so local binding spot id 0
    eleids[0] = std::make_pair( spheregid, 0 );

    // sphere current position
    std::vector<double> sphereeledisp;
    BEAMINTERACTION::UTILS::GetCurrentElementDis( Discret(), sphere,
        BeamInteractionDataState().GetDisColNp(), sphereeledisp);

    // note: sphere has just one node (with three translational dofs)
    for ( unsigned int dim = 0; dim < 3; ++dim )
      pos[0](dim) = sphere->Nodes()[0]->X()[dim] + sphereeledisp[dim];

    // loop over all integrins that are about to be bonded
    for ( auto const & bspotiter : newlinkiter.second )
    {
      int const beamgid = bspotiter.first;
      eleids[1] = std::make_pair( beamgid, bspotiter.second );

      // get neighboring element
      DRT::ELEMENTS::Beam3Base * beamele =
          dynamic_cast< DRT::ELEMENTS::Beam3Base * >( Discret().gElement(beamgid) );

      // beam bspot pos
      std::vector<double> beameledisp;
      BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis( Discret(), beamele,
          BeamInteractionDataState().GetDisColNp(), PeriodicBoundingBox(), beameledisp );
      beamele->GetPosOfBindingSpot( pos[1], beameledisp, bspotiter.second, PeriodicBoundingBox() );

      // create and initialize objects of beam-to-beam connections
      // Todo introduce enum for type of linkage (only linear Beam3r element possible so far)
      //      and introduce corresponding input parameter
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr =
        BEAMINTERACTION::BeamLinkPinJointed::Create();

      // unique linker id is bspot elegid and locspot id paired
      int id = BEAMINTERACTION::UTILS::CantorPairing( eleids[1] );

      // finally initialize and setup object
      linkelepairptr->Init( id, eleids, pos );
      // fixme
      linkelepairptr->Setup( 4 );

      // set on status on element level
      sphere->AddBond(eleids[1]);

      // add to my double bonds
      cell_beam_element_pairs_[ id ] = linkelepairptr;
    }
  }

}

///*-------------------------------------------------------------------------------*
// *-------------------------------------------------------------------------------*/
//void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::UpdateCellsPositionRandomly()
//{
//  CheckInit();
//
//  DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );
//  int const numcells = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
//  std::vector<double> randpos;
//  DRT::Problem::Instance()->Random()->Uni( randpos, 3 * numcells );
//
//  //todo: this is of course not nice, this needs to be done somewhere else
//  for( int i = 0; i < numcells; ++i )
//  {
//    DRT::Element* eleptr = Discret().gElement(EleTypeMapExtractorPtr()->SphereMap()->GID(i) );
//
//    std::vector<int> dofnode  = Discret().Dof(eleptr->Nodes()[0]);
//
//    // random position of cell inside box
//    static std::vector< std::vector< double > > Xnew(numcells, std::vector< double >(3) );
//    if( GState().GetStepN() % 20 == 0 and GState().GetStepN() != 0)
//    {
//      for ( int dim = 0; dim < 3; ++dim )
//      {
//        double edgelength = 5.5;
//        double min = 0.5;
//        Xnew[i][dim] = min + ( edgelength * randpos[ i + dim ] );
//      }
//    }
//
//    // loop over all dofs
//    for( int dim = 0; dim < 3; ++dim )
//    {
//      int doflid = BeamInteractionDataStatePtr()->GetMutableDisNp()->Map().LID(dofnode[dim]);
//      (*BeamInteractionDataStatePtr()->GetMutableDisNp() )[doflid] = Xnew[i][dim] - eleptr->Nodes()[0]->X()[dim];
//    }
// /*   int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(i);
//    DiscretPtr()->gElement(elegid)->Nodes()[0]->SetPos(Xnew);*/
//   }
//}


