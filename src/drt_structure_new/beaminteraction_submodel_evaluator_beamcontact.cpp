/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact.cpp

\brief class for submodel beam contact

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_TimeMonitor.hpp>

#include "beaminteraction_submodel_evaluator_beamcontact.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "str_timint_basedataglobalstate.H"
#include "str_utils.H"
#include "../drt_biopolynet/biopolynet_calc_utils.H"
#include "../drt_particle/particle_handler.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_beamcontact/beam_contact_params.H"
#include "../drt_beamcontact/beam_to_beam_interaction.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::BeamContact()
    : beam_contact_params_ptr_(Teuchos::null),
      bin_beamcontent_(INPAR::BINSTRATEGY::Beam),
      BTB_contact_elepairs_(Teuchos::null)
{
  // clear stl stuff
  nearby_elements_map_.clear();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::Setup()
{
  CheckInit();

  // init and setup beamt to beam contact data container
  beam_contact_params_ptr_ =Teuchos::rcp(new BEAMINTERACTION::BeamContactParams() );
  beam_contact_params_ptr_->Init();
  beam_contact_params_ptr_->Setup();

  // Todo really needed here? maybe find better place
  // ensure that contact is evaluated correctly at beginning of first time step (initial overlap)
  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamToBeamContactElementPairs();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::Reset()
{
  CheckInitSetup();

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    std::vector<const DRT::Element*> element_ptr(2);

    element_ptr[0] = elepairptr->Element1();
    element_ptr[1] = elepairptr->Element2();

    // element Dof values relevant for centerline interpolation
    std::vector< std::vector<double> > element_centerline_dofvec_absolutevalues(2);

    for (unsigned int ielement=0; ielement<2; ++ielement)
    {
      std::vector<double> eledispvec;

      // extract the Dof values of this element from displacement vector
      BIOPOLYNET::UTILS::GetCurrentElementDis(
          Discret(),
          element_ptr[ielement],
          BeamInteractionDataStatePtr()->GetMutableDisColNp(),
          eledispvec);

      const DRT::ELEMENTS::Beam3Base* beam_element_ptr =
          dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element_ptr[ielement]);

      // get the current absolute values for those Dofs relevant for centerline interpolation
      // initial values are added by element itself
      beam_element_ptr->ExtractCenterlineDofValuesFromElementStateVector(
          eledispvec,
          element_centerline_dofvec_absolutevalues[ielement],
          true);
    }

    // update the Dof values in the interaction element pair object
    elepairptr->ResetState(
        element_centerline_dofvec_absolutevalues[0],
        element_centerline_dofvec_absolutevalues[1]);
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateForce()
{
  CheckInitSetup();

  // resulting discrete element force vectors of the two interacting elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector< LINALG::SerialDenseVector > eleforce_centerlineDOFs(2);

  std::vector< std::vector<LINALG::SerialDenseMatrix> > dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero force values returned which need assembly?
  bool pair_is_active = false;


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamToBeamInteraction > >::const_iterator iter;
  for ( iter = BTB_contact_elepairs_.begin(); iter != BTB_contact_elepairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    pair_is_active = elepairptr->Evaluate(
        &eleforce_centerlineDOFs[0],
        &eleforce_centerlineDOFs[1],
        NULL,
        NULL,
        NULL,
        NULL);

    if (pair_is_active)
    {
      elegids[0] = elepairptr->Element1()->Id();
      elegids[1] = elepairptr->Element2()->Id();

      // assemble force vector affecting the centerline DoFs only
      // into element force vector ('all DoFs' format, as usual)
      BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
          Discret(),
          elegids,
          eleforce_centerlineDOFs,
          dummystiff,
          &eleforce,
          NULL);

      // Fixme
      eleforce[0].Scale(-1.0);
      eleforce[1].Scale(-1.0);

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
          eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null);
    }

  }
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateStiff()
{
  CheckInitSetup();

  // linearizations
  std::vector< std::vector<LINALG::SerialDenseMatrix> > elestiff( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // linearizations (centerline DOFs only!)
  std::vector< std::vector< LINALG::SerialDenseMatrix > > elestiff_centerlineDOFs ( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  std::vector< LINALG::SerialDenseVector > dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamToBeamInteraction > >::const_iterator iter;
  for ( iter = BTB_contact_elepairs_.begin(); iter != BTB_contact_elepairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    pair_is_active = elepairptr->Evaluate(
        NULL,
        NULL,
        &elestiff_centerlineDOFs[0][0],
        &elestiff_centerlineDOFs[0][1],
        &elestiff_centerlineDOFs[1][0],
        &elestiff_centerlineDOFs[1][1]);

    if (pair_is_active)
    {
      elegids[0] = elepairptr->Element1()->Id();
      elegids[1] = elepairptr->Element2()->Id();

      // assemble stiffness matrix affecting the centerline DoFs only
      // into element stiffness matrix ('all DoFs' format, as usual)
      BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
          Discret(),
          elegids,
          dummyforce,
          elestiff_centerlineDOFs,
          NULL,
          &elestiff);

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(),
          elegids, dummyforce, elestiff, Teuchos::null, BeamInteractionDataStatePtr()->GetMutableStiff());
    }

  }
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateForceStiff()
{
  CheckInitSetup();

  // resulting discrete element force vectors of the two interacting elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector< LINALG::SerialDenseVector > eleforce_centerlineDOFs(2);

  // linearizations
  std::vector< std::vector<LINALG::SerialDenseMatrix> > elestiff( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // linearizations (centerline DOFs only!)
  std::vector< std::vector< LINALG::SerialDenseMatrix > > elestiff_centerlineDOFs ( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamToBeamInteraction > >::const_iterator iter;
  for ( iter = BTB_contact_elepairs_.begin(); iter != BTB_contact_elepairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    elegids[0] = elepairptr->Element1()->Id();
    elegids[1] = elepairptr->Element2()->Id();

    pair_is_active = elepairptr->Evaluate(
        &eleforce_centerlineDOFs[0],
        &eleforce_centerlineDOFs[1],
        &elestiff_centerlineDOFs[0][0],
        &elestiff_centerlineDOFs[0][1],
        &elestiff_centerlineDOFs[1][0],
        &elestiff_centerlineDOFs[1][1] );

    if (pair_is_active)
    {
      elegids[0] = elepairptr->Element1()->Id();
      elegids[1] = elepairptr->Element2()->Id();

      // assemble force vector and stiffness matrix affecting the centerline DoFs only
      // into element force vector and stiffness matrix ('all DoFs' format, as usual)
      BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
          Discret(),
          elegids,
          eleforce_centerlineDOFs,
          elestiff_centerlineDOFs,
          &eleforce,
          &elestiff);

      // Fixme
      eleforce[0].Scale(-1.0);
      eleforce[1].Scale(-1.0);

      // assemble the contributions into force vector class variable
      // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
      BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(),
          elegids, eleforce, elestiff, BeamInteractionDataStatePtr()->GetMutableForceNp(),
          BeamInteractionDataStatePtr()->GetMutableStiff() );
    }

  }

//  PrintActiveBeamToBeamContactSet(std::cout);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::PreUpdateStepElement()
{
  CheckInitSetup();

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::UpdateStepElement()
{
  CheckInitSetup();

//  PrintActiveBeamToBeamContactSet(std::cout);

  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamToBeamContactElementPairs();

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::PostUpdateStepElement()
{
  CheckInitSetup();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{


}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::FindAndStoreNeighboringElements()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::GetNeighboringEles");

#ifdef DEBUG
  if(static_cast<int>(BeamInteractionDataStatePtr()->GetExtEleToBinMap().size()) != DiscretPtr()->ElementColMap()->NumMyElements())
    dserror("std::map does not equal elecolmap (check e.g. if extended ghosting contains "
            " standard ghosting). Therefore not every contact can be detected");
#endif

  // loop over all row elements
  const int numroweles = DiscretPtr()->NumMyRowElements();
  for(int rowele_i=0; rowele_i<numroweles; ++rowele_i)
  {
    const int elegid = DiscretPtr()->ElementRowMap()->GID(rowele_i);
    DRT::Element* currele = DiscretPtr()->gElement(elegid);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all bind touched by currele
    std::set<int>::const_iterator biniter;
    for(biniter=BeamInteractionDataStatePtr()->GetExtEleToBinSet(elegid).begin();
        biniter!=BeamInteractionDataStatePtr()->GetExtEleToBinSet(elegid).end(); ++biniter)
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      ParticleHandlerPtr()->BinStrategy()->GetNeighborAndOwnBinIds(
          *biniter,loc_neighboring_binIds);

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert(loc_neighboring_binIds.begin(),
          loc_neighboring_binIds.end());
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds(neighboring_binIds.begin(),
        neighboring_binIds.end());

    // set of elements that lie in neighboring bins
    std::set<DRT::Element*> neighboring_elements;
    ParticleHandlerPtr()->GetBinContent(neighboring_elements,
        bin_beamcontent_,glob_neighboring_binIds);

    // sort out elements that should not be considered in contact evaluation
    SelectElesToBeConsideredForContactEvaluation(currele, neighboring_elements);

    nearby_elements_map_[elegid] = neighboring_elements;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::SelectElesToBeConsideredForContactEvaluation(
    DRT::Element*            currele,
    std::set<DRT::Element*>& neighbors) const
{
  CheckInit();

  // sort out elements that should not be considered in contact evaluation
  std::set<DRT::Element*>::iterator eiter;
  for( eiter = neighbors.begin(); eiter != neighbors.end();)
  {
    bool toerase = false;
    // 1) ensure each contact only evaluated once (keep in mind that we are
    //    using FEMatrices and FEvectors -> || (*eiter)->Owner() != myrank not necessary)
    if( not ( currele->Id() < (*eiter)->Id() ) )
      toerase = true;

    // 2) ensure that two elements sharing the same node do not get into contact
    for ( int i = 0; i < 2; i++ )
      for ( int j = 0; j < 2 ; j++ )
        if( (*eiter)->NodeIds()[i] == currele->NodeIds()[j] )
          toerase = true;

    if( toerase )
      neighbors.erase(eiter++);
    else
      ++eiter;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::
    CreateBeamToBeamContactElementPairs()
{
  // Todo maybe keep existing pairs and reuse them ?
  BTB_contact_elepairs_.clear();

  std::map<int, std::set<DRT::Element*> >::const_iterator nearbyeleiter;

  for (nearbyeleiter=nearby_elements_map_.begin(); nearbyeleiter!=nearby_elements_map_.end(); ++nearbyeleiter)
  {
    const int elegid = nearbyeleiter->first;
    const DRT::Element* firsteleptr = DiscretPtr()->gElement(elegid);

    const DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(firsteleptr);

    // at the moment, both elements of a beam contact pair must be of same type Todo
    const unsigned int numnodes_centerline = beamele1->NumCenterlineNodes();
    const unsigned int numnodalvalues = beamele1->HermiteCenterlineInterpolation() ? 2 : 1;

    std::set<DRT::Element*>::const_iterator secondeleiter;
    for (secondeleiter=nearbyeleiter->second.begin(); secondeleiter!=nearbyeleiter->second.end(); ++secondeleiter)
    {
      const DRT::Element* secondeleptr = *secondeleiter;

      Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> newbeaminteractionpair =
          BEAMINTERACTION::BeamToBeamInteraction::Create(numnodes_centerline,numnodalvalues);

      newbeaminteractionpair->Init(
          beam_contact_params_ptr_,
          firsteleptr,
          secondeleptr);

      newbeaminteractionpair->Setup();

      BTB_contact_elepairs_.push_back(newbeaminteractionpair);
    }
  }

  IO::cout(IO::standard) <<"\n\n************************************************"<<IO::endl;
  IO::cout(IO::standard) << "PID " << GState().GetMyRank() << " currently monitors " <<
      BTB_contact_elepairs_.size() << " beam contact pairs" << IO::endl;
  IO::cout(IO::standard) <<"************************************************"<<IO::endl;

}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::
    PrintAllBeamToBeamContactElementPairs( std::ostream& out ) const
{
  out << "\n\nCurrent BeamToBeamContactElementPairs: ";
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
    (*iter)->Print(out);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::
    PrintActiveBeamToBeamContactSet( std::ostream& out ) const
{
  out << "\n    Active BeamToBeam Contact Set (PID " << GState().GetMyRank() << "):-----------------------------------------\n";
  out << "    ID1            ID2              T xi       eta      angle    gap         force\n";

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
    (*iter)->PrintSummaryOneLinePerActiveSegmentPair(out);

  out << std::endl;
}
