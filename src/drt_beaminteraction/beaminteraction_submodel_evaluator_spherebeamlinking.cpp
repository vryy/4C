/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_spherebeamlinking.cpp

\brief class for managing rigid sphere to beam crosslinking

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_spherebeamlinking.H"
#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_crosslinking.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"
#include "../drt_beaminteraction/crosslinker_node.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_beaminteraction/beam_link_pinjointed.H"
#include "../drt_beaminteraction/beam_link_beam3r_lin2_pinjointed.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"
#include "../drt_beaminteraction/spherebeamlinking_params.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_beaminteraction.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/runtime_vtp_writer.H"
#include "../drt_structure_new/str_timint_basedataio_runtime_vtp_output.H"
#include "../drt_io/io_control.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_particle/particle_handler.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_new/str_timint_basedataio.H"

#include "../drt_mat/crosslinkermat.H"

#include <Teuchos_TimeMonitor.hpp>





/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::SphereBeamLinking() :
    sm_crosslinkink_ptr ( Teuchos::null ),
    spherebeamlinking_params_ptr_ ( Teuchos::null )
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::Setup()
{
  CheckInit();

  // construct, init and setup data container for crosslinking
  spherebeamlinking_params_ptr_ = Teuchos::rcp( new BEAMINTERACTION::SphereBeamLinkingParams() );
  spherebeamlinking_params_ptr_->Init( GState() );
  spherebeamlinking_params_ptr_->Setup();

  // build runtime vtp writer
  if ( GInOutput().GetRuntimeVtpOutputParams() != Teuchos::null )
    InitOutputRuntimeVtp();

  // set flag
  issetup_ = true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PostSetup()
{
  CheckInitSetup();
 // nothing to do (yet)
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::InitSubmodelDependencies(
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
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::Reset()
{
  CheckInitSetup();

  // reset crosslinker pairs
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      // get elements
#ifdef DEBUG
      if ( sphere != dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elepairptr->GetEleGid(0) ) ) )
        dserror(" Rigid Sphere element has stored wrong linker. ");
#endif

      DRT::ELEMENTS::Beam3Base const * beamele =
          dynamic_cast< DRT::ELEMENTS::Beam3Base const* >( Discret().gElement( elepairptr->GetEleGid(1) ) );

      // init position of linker nodes
      std::vector< LINALG::Matrix< 3, 1 > > pos( 2, LINALG::Matrix< 3, 1 >(true) );

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

      // dummy triad
      std::vector< LINALG::Matrix< 3, 3 > > dummy_triad( 2, LINALG::Matrix< 3, 3 >(true) );

      // finally reset state
      elepairptr->ResetState( pos, dummy_triad );
    }
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::EvaluateForce()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector< LINALG::SerialDenseVector > bspotforce( 2, LINALG::SerialDenseVector(6) );

  // resulting discrete element force vectors of the two parent elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  std::vector< std::vector<LINALG::SerialDenseMatrix> > dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      for ( unsigned int i = 0; i < 2; ++i )
      {
        elegids[i] = elepairptr->GetEleGid(i);
        bspotforce[i].Zero();
      }

      // evaluate beam linkage object to get forces of binding spots
      elepairptr->EvaluateForce( bspotforce[0], bspotforce[1] );

      // apply forces on binding spots to parent elements
      // and get their discrete element force vectors
      BEAMINTERACTION::UTILS::ApplyBindingSpotForceToParentElements< DRT::ELEMENTS::Rigidsphere,
          DRT::ELEMENTS::Beam3Base >( Discret(), PeriodicBoundingBoxPtr(),
          BeamInteractionDataStatePtr()->GetMutableDisColNp(), elepairptr, bspotforce, eleforce );

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
          eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null );
    }
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::EvaluateStiff()
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

  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

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
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::EvaluateForceStiff()
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
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

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
  }

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

}
/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PreUpdateStepElement(
    bool beam_redist )
{
  CheckInitSetup();
  // not repartition of binning discretization necessary
  return false;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateStepElement(
    bool repartition_was_done )
{
  CheckInitSetup();

  // update linker
  if ( spherebeamlinking_params_ptr_->MaxNumberOfIntegrinsPerCell() > 0 )
  {
    // consider new bonds
    std::map< int, std::vector< std::pair< int, int > > > newlinks;
    FindAndStoreNeighboringElements( newlinks );
    CreateBeamToSphereJoint( newlinks );

    // consider possible unbinding
    UnbindSphereBeamBonds();

    // consider sphere linker contraction
    UpdateLinkerLength();
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PostUpdateStepElement()
{
  CheckInitSetup();

  // print total number of linker across all procs
  int num_local_linker = 0;
  int num_global_linker = 0;
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    num_local_linker += sphere->GetNumberOfBonds();
  }
  // build sum over all procs
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
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::RuntimeOutputStepState() const
{
  CheckInitSetup();

  if ( vtp_writer_ptr_ != Teuchos::null )
    WriteOutputRuntimeVtp();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::WriteRestart(
    IO::DiscretizationWriter & ia_writer,
    IO::DiscretizationWriter & bin_writer) const
{
  CheckInitSetup();

  // as bonds are stored in rigid sphere element, nothing to do here
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::ReadRestart(
    IO::DiscretizationReader & ia_reader,
    IO::DiscretizationReader & bin_reader)
{
  CheckInitSetup();

  // as bonds are stored in rigid sphere element, nothing to do here
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::PostReadRestart()
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::AddBinsToBinColMap(
    std::set< int >& colbins)
{
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::
    AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins ) const
{
  CheckInitSetup();
  // empty
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::
    GetHalfInteractionDistance( double & half_interaction_distance )
{
  double spherebeamlinking_half_interaction_distance = 0.0;

  // loop over all sphere elements (needed in case interaction distance should be
  // radius dependent in the future)
  double curr_ia_dist = 0.0;
  int unsigned const numrowsphereeles = EleTypeMapExtractor().SphereMap()->NumMyElements();
  for( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
//    int const elegid = EleTypeMapExtractor().SphereMap()->GID(rowele_i);
//    DRT::ELEMENTS::Rigidsphere * sphere =
//        dynamic_cast< DRT::ELEMENTS::Rigidsphere * >( Discret().gElement(elegid) );

    curr_ia_dist = 0.5 * spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLengthTolerance();

    // update distance
    spherebeamlinking_half_interaction_distance = ( curr_ia_dist > spherebeamlinking_half_interaction_distance ) ?
        curr_ia_dist : spherebeamlinking_half_interaction_distance;
  }

  // get global maximum
  double spherebeamlinking_half_interaction_distance_global = 0.0;
  // build sum over all procs
  MPI_Allreduce( &spherebeamlinking_half_interaction_distance, &spherebeamlinking_half_interaction_distance_global,
      1, MPI_DOUBLE, MPI_MAX, dynamic_cast<const Epetra_MpiComm*>( &(Discret().Comm() ) )->Comm() );

  if ( GState().GetMyRank() == 0 )
    std::cout << " spherebeamlinking half interaction distance " << spherebeamlinking_half_interaction_distance_global << std::endl;


  half_interaction_distance = ( spherebeamlinking_half_interaction_distance_global > half_interaction_distance ) ?
      spherebeamlinking_half_interaction_distance_global : half_interaction_distance;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::InitOutputRuntimeVtp()
{
  CheckInit();

  vtp_writer_ptr_ = Teuchos::rcp( new RuntimeVtpWriter() );

  // determine path of output directory
  const std::string outputfilename( DRT::Problem::Instance()->OutputControlFile()->FileName() );

  size_t pos = outputfilename.find_last_of("/");

  if ( pos == outputfilename.npos )
    pos = 0ul;
  else
    pos++;

  const std::string output_directory_path( outputfilename.substr(0ul, pos) );

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  vtp_writer_ptr_->Initialize(
      BinDiscretPtr()->Comm().MyPID(),
      BinDiscretPtr()->Comm().NumProc(),
      num_timesteps_in_simulation_upper_bound,
      output_directory_path,
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      "spherebeamlinker",
      DRT::Problem::Instance()->OutputControlFile()->RestartName(),
      GState().GetTimeN(),
      GInOutput().GetRuntimeVtpOutputParams()->WriteBinaryOutput()
      );
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::WriteOutputRuntimeVtp() const
{
  CheckInitSetup();

  // get number of linker on current proc
  unsigned int num_row_points = 0;
  int unsigned const numrowsphereeles = EleTypeMapExtractor().SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractor().SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    num_row_points += sphere->GetNumberOfBonds();
  }

  // set geometry manually
  const unsigned int num_spatial_dimensions = 3;

  // get and prepare storage for point coordinate values
  std::vector<double>& point_coordinates = vtp_writer_ptr_->GetMutablePointCoordinateVector();
  point_coordinates.clear();
  point_coordinates.reserve( num_spatial_dimensions * num_row_points );

  // init output desired output vectors
  std::vector<double> currlength ( num_row_points, 0.0 );
  std::vector<double> orientation ( num_row_points * num_spatial_dimensions, 0.0 );

  // set position of spherebeamlinks
  unsigned int bond_i = 0;
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractor().SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere const * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere const* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

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
      PeriodicBoundingBox().UnShift3D( pos[1], pos[0] );

      // set point coordiates value
      for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
        point_coordinates.push_back( pos[0](idim) + 0.5 * ( pos[1](idim) - pos[0](idim) ) );

      // set desired output vectors
      for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
        orientation[bond_i * num_spatial_dimensions + idim] = ( pos[1](idim) - pos[0](idim) );

      ++bond_i;
    }
  }

  // reset time and time step and geometry name in the writer object
  vtp_writer_ptr_->SetupForNewTimeStepAndGeometry(
      GState().GetTimeN(), GState().GetStepN(), "spherebeamlinker" );

  // append all desired output data to the writer object's storage
  // i) number of bonds: collect data and append to visualization results if desired
  vtp_writer_ptr_->AppendVisualizationPointDataVector( orientation, 3, "orientation" );

  // finalize everything and write all required VTU files to filesystem
  vtp_writer_ptr_->WriteFiles();

  // write a collection file summarizing all previously written files
  vtp_writer_ptr_->WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() +
      "-" + "spherebeamlinker" );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::FindAndStoreNeighboringElements(
    std::map< int, std::vector< std::pair< int, int > > > & newlinks )
{
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::"
                            "SphereBeamLinking::FindAndStoreNeighboringElements");

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
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::CheckFeasibilityOfNewLink(
    DRT::Element const * currele,
    std::set< DRT::Element * > const & neighbors,
    std::unordered_set< int > & tobebonded,
    std::map< int, std::vector< std::pair< int, int > > > & newlinks ) const
{
  CheckInitSetup();

  int numnewbondsthisstep = 0;

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
  for ( auto const & eiter : neighbors )
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
      // build unique linker id from elegid and local binding spot id
      std::pair< int, int > bspotpair = std::make_pair( beamele->Id(), bspot_i );
      int bpspotgid = BEAMINTERACTION::UTILS::CantorPairing( bspotpair);

      // 0. criterion: has sphere reached maximum number of admissible bonds?
      if ( ( sphere->GetNumberOfBonds() + numnewbondsthisstep ) ==
             spherebeamlinking_params_ptr_->MaxNumberOfIntegrinsPerCell() )
        continue;

#ifdef DEBUG
      // todo: do only in debug mode as soon tested enough
      // safety check
      if ( ( sphere->GetNumberOfBonds() + numnewbondsthisstep ) >
             spherebeamlinking_params_ptr_->MaxNumberOfIntegrinsPerCell() )
        dserror(" sphere has more bonds than allowed. Something went wrong.");
#endif

      // 1. criterion: does identical bond already exist?
      if ( sphere->DoesBondExist( bpspotgid ) )
        continue;

      // 2. criterion: is beam binding spot free (covered by criterion 1? Only in case of separate linker
      // for cell to beam and beam to beam
      // fixme: maybe replace beam to beam in this case
      // search literature if binding spots are the same


      // 3. criterion: distance
      // get current position at binding spot xi
      bspotpos.Clear();
      beamele->GetPosOfBindingSpot( bspotpos, beameledisp, bspot_i, PeriodicBoundingBox() );

      double const linkdistmin = spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLength()
                                 - spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLengthTolerance();
      double const linkdistmax = spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLength()
                                 + spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingLengthTolerance();


      if ( BEAMINTERACTION::UTILS::IsDistanceOutOfRange( spherepos, bspotpos, linkdistmin, linkdistmax ) )
        continue;

      // 4. criterion: orientation
      /* NOTE: this works slightly different to crosslinking of two beams: here we check the angle
       * between the beams first base vector and the direction vector from the spheres mid point to
       * the binding spot, i.e. LINKINGANGLE in the input file means something slightly different in
       * this case */

      // get current triad at binding spot xi
      bspottriad.Clear();
      beamele->GetTriadOfBindingSpot( bspottriad, beameledisp, bspot_i );

      // note: we use first base vector instead of tangent vector here
      LINALG::Matrix< 3, 1> curr_bindingspot_beam_tangent(true);
      for ( unsigned int idim = 0; idim < 3; ++idim )
        curr_bindingspot_beam_tangent(idim) = bspottriad( idim, 0 );

      LINALG::Matrix< 3, 1 > dist_vec(true);
      dist_vec.Update( 1.0, bspotpos, -1.0, spherepos );

      double const linkanglemin = spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngle()
                                  - spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngleTolerance();
      double const linkanglemax = spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngle()
                                  + spherebeamlinking_params_ptr_->GetLinkerMaterial()->LinkingAngleTolerance();

      if ( BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
          dist_vec, curr_bindingspot_beam_tangent, linkanglemin, linkanglemax ) )


      // 5. criterion: check if bspot will already be bonded this step
      if ( tobebonded.find(bpspotgid) != tobebonded.end() )
        continue;

      // update control variables
      tobebonded.insert(bpspotgid);
      ++numnewbondsthisstep;

      // add new link
      newlinks[sphere->Id()].push_back( bspotpair );

    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::CreateBeamToSphereJoint(
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
      Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > linkelepairptr =
        BEAMINTERACTION::BeamLinkPinJointed::Create( INPAR::BEAMINTERACTION::truss );

      // unique linker id is bspot elegid and locspot id paired
      int id = BEAMINTERACTION::UTILS::CantorPairing( eleids[1] );

      // dummy triad
      std::vector< LINALG::Matrix< 3, 3 > > dummy_triad( 2, LINALG::Matrix< 3, 3 >(true) );

      // finally initialize and setup object
      linkelepairptr->Init( id, eleids, pos, dummy_triad, GState().GetTimeNp() );
      // material id
      linkelepairptr->Setup( spherebeamlinking_params_ptr_->GetLinkerMaterial()->BeamElastHyperMatNum() );

      // set on status on element level
      sphere->AddBond( id, linkelepairptr );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UnbindSphereBeamBonds()
{
  CheckInitSetup();

  // init variables
  double p_unbind = 0.0;

  // adapt force/strain in linker
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere* >( Discret().gElement( elegid ) );

    // loop over bonds of current sphere
    std::vector< int > to_dissolve;
    to_dissolve.reserve( sphere->GetBondMap().size() );
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> elepairptr = ipair.second;

      // consider catch-slip bond behavior of integrin linkers
      if ( spherebeamlinking_params_ptr_->DeltaTime() > 0 )
      {
        CalcForceDependentCatchSlipBondUnbindProbability( elepairptr, p_unbind );
      }
      else
      {
        p_unbind = 0.0;
      }

      // if probability criterion is not met, we are done here
      if ( DRT::Problem::Instance()->Random()->Uni() > p_unbind )
        continue;

      to_dissolve.push_back( ipair.first );
    }

    std::cout << " unbound " << to_dissolve.size() << std::endl;

    // dissolve all bonds
    for ( unsigned int i_td = 0; i_td < to_dissolve.size(); ++i_td )
      sphere->DissolveBond( to_dissolve[i_td] );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::CalcForceDependentCatchSlipBondUnbindProbability(
    Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed> linkelepairptr,
    double & p_unbind) const
{
  CheckInitSetup();

  // check if link was set in current time step
  if ( linkelepairptr->GetTimeLinkWasSet() == GState().GetTimeNp() )
  {
    p_unbind = 0.0;
    return;
  }

  // Note: this needs to come after a contraction was equilibrated by the network. Doing this directly
  // after changing the linker reference length does not make sense.

  // todo: maybe add these to linker material input line
  double const & phi_FA_s = 7.78;
  double const & phi_FA_c = 4.02;
  double const & f_FA = 5.38;
  double const & k_FA_0 = 1;

  double dt = spherebeamlinking_params_ptr_->DeltaTime();

  // Fixme: is force 1 the correct one, do we need to take the force on the end that is connected to the beam?
  // note: as we only check unbinding in links that were set before the current time step, we do not need to
  // caclulate the forces again.
  LINALG::SerialDenseVector bspotforce_one(6);
  linkelepairptr->GetBindingSpotForce( 1, bspotforce_one );
  double f = bspotforce_one.Norm2();


  // check if linker is stretched -> sgn+ or compressed -> sgn- by checking orientation of force vector
  // note: this works only if there are no other forces (like inertia, stochastic, damping) acting on the linker
  LINALG::Matrix<3,1> dist_vec(true);
  LINALG::Matrix<3,1> bspotforceone(true);
  dist_vec.Update( -1.0, linkelepairptr->GetBindSpotPos1(), 1.0, linkelepairptr->GetBindSpotPos2() );
  for ( unsigned int j = 0; j < 3; ++j )
    bspotforceone(j) = bspotforce_one(j);
  double sgn = ( dist_vec.Dot( bspotforceone ) < 0.0 ) ? -1.0 : 1.0;

  /* alternative for linear centerline interpolation would be to compare
   * reference length and current length to see if element is stretched or compressed
   */

  // fixme: is that correct / does that make sense for compressive forces?
  f *= sgn;

  // calculate force dependent off rate for catch slip bond
  // ( see Wang et al. 2016 "Mechanosensitive subcellular rheostasis drives emergent
  //   single-cell mechanical homeostasis", nature materials. See supplementary information)
  double k_FA_off = k_FA_0 * ( exp( f / f_FA - phi_FA_s) + exp( ( (-1.0) * f / f_FA ) + phi_FA_c ) );

  // get respective force dependent unbind probability for each cl binding spot
  if ( std::isfinite( k_FA_off ) )
  {
    p_unbind = 1.0 - exp( (-1.0) * dt * k_FA_off );
  }
  else
  {
    p_unbind = 1.0;
    std::cout << "WARNING: You have some very high forces acting on your integrins. Are you "
        "sure this is what you want? " << std::endl;
  }

  std::cout << " force " << f << std::endl;
  std::cout << " koff: " << k_FA_off << std::endl;
  std::cout << " p_unbind " << p_unbind << std::endl;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateLinkerLength()
{
  CheckInitSetup();

  // adapt force/strain in linker
  // note: problem time step is used here
  double contraction_per_dt =
      spherebeamlinking_params_ptr_->ContractionRate() * (*GState().GetDeltaTime())[0];
  double scalefac = 0.0;
  int unsigned const numrowsphereeles = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  for ( unsigned int rowele_i = 0; rowele_i < numrowsphereeles; ++rowele_i )
  {
    int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(rowele_i);
    DRT::ELEMENTS::Rigidsphere * sphere =
        dynamic_cast< DRT::ELEMENTS::Rigidsphere * >( Discret().gElement( elegid ) );

    // todo: do we want this
    // change radius of contracting cell
//    if ( GState().GetStepN() % 10 == 0 )
//      sphere->ScaleRadius(0.95);

    // loop over bonds of current sphere
    for ( auto const & ipair : sphere->GetBondMap() )
    {
      // get pair object
      Teuchos::RCP< BEAMINTERACTION::BeamLinkPinJointed > elepairptr = ipair.second;

      // only contract if linker size > sphere radius
      if ( ( elepairptr->GetCurrentLinkerLength() <= sphere->Radius() ) and ( GState().GetStepN() == 0 ) )
        continue;

      // compute scaling factor for linker length
      scalefac = 1.0 - ( contraction_per_dt / elepairptr->GetReferenceLength() );

      // some safety checks
      if ( scalefac < 0.0 )
        dserror( "Contraction %f of a linker of more than its reference length %f in one time does not make sense.",
            contraction_per_dt, elepairptr->GetReferenceLength() );
      if ( contraction_per_dt > elepairptr->GetCurrentLinkerLength() )
        dserror( "Contraction of %f for a linker with current length %f does not make sense.",
            contraction_per_dt, elepairptr->GetCurrentLinkerLength() );

      // scale linker length / contract linker
      elepairptr->ScaleLinkerReferenceLength( scalefac );
    }
  }
}

///*-------------------------------------------------------------------------------*
// *-------------------------------------------------------------------------------*/
//void BEAMINTERACTION::SUBMODELEVALUATOR::SphereBeamLinking::UpdateCellsPositionRandomly()
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


