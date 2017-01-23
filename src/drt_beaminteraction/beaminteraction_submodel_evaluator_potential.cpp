/*-----------------------------------------------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_potential.cpp

\brief submodel for potential-based beam interactions

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_potential.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_new/str_utils.H"

#include "../drt_particle/particle_handler.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

//#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_beaminteraction/beam3contact_utils.H"
#include "../drt_beaminteraction/biopolynet_calc_utils.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"
#include "beam_potential_pair.H"
#include "beam_potential_params.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::BeamPotential()
    : beam_potential_params_ptr_(Teuchos::null),
      bin_beamcontent_(BINSTRATEGY::UTILS::Beam),
      beam_potential_element_pairs_(Teuchos::null)
{
  // clear stl stuff
  nearby_elements_map_.clear();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::Setup()
{
  CheckInit();

  // init and setup beam to beam contact data container
  beam_potential_params_ptr_ =Teuchos::rcp(new BEAMINTERACTION::BeamPotentialParams() );
  BeamPotentialParams().Init();
  BeamPotentialParams().Setup();

  PrintConsoleWelcomeMessage(std::cout);

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PostSetup()
{
  CheckInitSetup();

  // Todo really needed here? maybe find better place
  // ensure that potential-based interactions are evaluated correctly at beginning of first time step (initially)
  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamPotentialElementPairs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::InitSubmodelDependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map>  const submodelmap)
{
  CheckInitSetup();
  // no active influence on other submodels
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::Reset()
{
  CheckInitSetup();

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> >::const_iterator iter;
  for (iter=beam_potential_element_pairs_.begin(); iter!=beam_potential_element_pairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

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

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateForce()
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


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamPotentialPair > >::const_iterator iter;
  for ( iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr,conditions_element1,conditions_element2);

    for (unsigned int k=0; k<conditions_element1.size(); ++k)
    {
      int npotlaw1 = conditions_element1[k]->GetInt("potlaw");

      for (unsigned int j=0; j<conditions_element2.size(); ++j)
      {
        int npotlaw2 = conditions_element2[j]->GetInt("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if(npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror("number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active = elepairptr->Evaluate(
              &eleforce_centerlineDOFs[0],
              &eleforce_centerlineDOFs[1],
              NULL,
              NULL,
              NULL,
              NULL,
              currconds,
              BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1-1),
              BeamPotentialParams().PotentialLawExponents().at(npotlaw1-1)
              );

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
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

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(), elegids,
                eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null);
          }

        }
      }
    }

  }
  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateStiff()
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


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamPotentialPair > >::const_iterator iter;
  for ( iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr,conditions_element1,conditions_element2);

    for (unsigned int k=0; k<conditions_element1.size(); ++k)
    {
      int npotlaw1 = conditions_element1[k]->GetInt("potlaw");

      for (unsigned int j=0; j<conditions_element2.size(); ++j)
      {
        int npotlaw2 = conditions_element2[j]->GetInt("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if(npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror("number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active = elepairptr->Evaluate(
              NULL,
              NULL,
              &elestiff_centerlineDOFs[0][0],
              &elestiff_centerlineDOFs[0][1],
              &elestiff_centerlineDOFs[1][0],
              &elestiff_centerlineDOFs[1][1],
              currconds,
              BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1-1),
              BeamPotentialParams().PotentialLawExponents().at(npotlaw1-1)
              );

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
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
      }
    }

  }
  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateForceStiff()
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


  std::vector< Teuchos::RCP< BEAMINTERACTION::BeamPotentialPair > >::const_iterator iter;
  for ( iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end(); ++iter )
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    elegids[0] = elepairptr->Element1()->Id();
    elegids[1] = elepairptr->Element2()->Id();

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr,conditions_element1,conditions_element2);

    for (unsigned int k=0; k<conditions_element1.size(); ++k)
    {
      int npotlaw1 = conditions_element1[k]->GetInt("potlaw");

      for (unsigned int j=0; j<conditions_element2.size(); ++j)
      {
        int npotlaw2 = conditions_element2[j]->GetInt("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if(npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror("number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active = elepairptr->Evaluate(
              &eleforce_centerlineDOFs[0],
              &eleforce_centerlineDOFs[1],
              &elestiff_centerlineDOFs[0][0],
              &elestiff_centerlineDOFs[0][1],
              &elestiff_centerlineDOFs[1][0],
              &elestiff_centerlineDOFs[1][1],
              currconds,
              BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1-1),
              BeamPotentialParams().PotentialLawExponents().at(npotlaw1-1)
              );

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
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

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix( Discret(),
                elegids, eleforce, elestiff, BeamInteractionDataStatePtr()->GetMutableForceNp(),
                BeamInteractionDataStatePtr()->GetMutableStiff() );
          }

        }
      }
    }


  }

//  PrintActiveBeamPotentialSet(std::cout);

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PreUpdateStepElement()
{
  CheckInitSetup();

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::UpdateStepElement()
{
  CheckInitSetup();

  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamPotentialElementPairs();

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PostUpdateStepElement()
{
  CheckInitSetup();

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();
  // nothing to do (so far)

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::AddBinsToBinColMap(
    std::set< int >& colbins)
{
  CheckInitSetup();
  // nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins ) const
{
  CheckInitSetup();
  // nothing to do
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::FindAndStoreNeighboringElements()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::FindAndStoreNeighboringElements");

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
    for(biniter=BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).begin();
        biniter!=BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).end(); ++biniter)
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
    std::vector< BINSTRATEGY::UTILS::BinContentType > bc(1);
    bc[0] = bin_beamcontent_;
    ParticleHandlerPtr()->BinStrategy()->GetBinContent(neighboring_elements,
        bc,glob_neighboring_binIds);

    // sort out elements that should not be considered in contact evaluation
    SelectElesToBeConsideredForPotentialEvaluation(currele, neighboring_elements);

    nearby_elements_map_[elegid] = neighboring_elements;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    SelectElesToBeConsideredForPotentialEvaluation(
    DRT::Element*            currele,
    std::set<DRT::Element*>& neighbors) const
{
  CheckInit();

  // sort out elements that should not be considered in potential evaluation
  std::set<DRT::Element*>::iterator eiter;
  for( eiter = neighbors.begin(); eiter != neighbors.end();)
  {
    bool toerase = false;

    DRT::Element* currneighborele = *eiter;

    // 1) ensure each interaction is only evaluated once (keep in mind that we are
    //    using FEMatrices and FEvectors -> || (*eiter)->Owner() != myrank not necessary)
    if ( not ( currele->Id() < currneighborele->Id() ) )
    {
      toerase = true;
    }

    // 2) exclude "self-interaction", i.e. a pair of elements on the same physical beam
    //    TODO introduce flag for self-interaction in input file
    else
    {
    // get the conditions applied to both elements of the pair and decide whether they need to be evaluated
    std::vector<DRT::Condition*> conds1, conds2;

    // since only the nodes know about their conditions, we need this workaround
    // we assume that a linecharge condition is always applied to the entire physical beam, i.e. it is sufficient to check only one node
    DRT::Node** nodes1;
    DRT::Node** nodes2;
    nodes1 = currele->Nodes();
    nodes2 = currneighborele->Nodes();

    dsassert(nodes1 != NULL and nodes2 != NULL, "pointer to nodes is NULL!");
    dsassert(nodes1[0] != NULL and nodes2[0] != NULL, "pointer to nodes is NULL!");

    nodes1[0]->GetCondition("BeamPotentialLineCharge", conds1);

    // get correct condition for beam or rigid sphere element
    if ( BEAMCONTACT::BeamElement(*currneighborele) )
      nodes2[0]->GetCondition("BeamPotentialLineCharge", conds2);
    else
      dserror("Only beam-to-beampotential-based interaction is implemented yet. "
          "No other types of elements allowed!");
//    else if(BEAMCONTACT::RigidsphereElement(*currneighborele)
//      nodes2[0]->GetCondition("RigidspherePotentialPointCharge", conds2);

    // validinteraction == true includes: both eles "loaded" by a charge condition of same potential law
    bool validinteraction = false;

    for (unsigned int i=0; i<conds1.size(); ++i)
    {
      int npotlaw1 = conds1[i]->GetInt("potlaw");

      for (unsigned int j=0; j<conds2.size(); ++j)
      {
        int npotlaw2 = conds2[j]->GetInt("potlaw");

        // here, we also exclude "self-interaction", i.e. a pair of elements on the same physical beam
        // TODO introduce flag for self-interaction in input file
        if (conds1[i] != conds2[j] and npotlaw1 == npotlaw2)
          validinteraction = true;
      }
    }

    if (not validinteraction)
      toerase = true;
    }


    if( toerase )
      neighbors.erase(eiter++);
    else
      ++eiter;
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    CreateBeamPotentialElementPairs()
{
  // Todo maybe keep existing pairs and reuse them ?
  beam_potential_element_pairs_.clear();

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

      Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> newbeaminteractionpair =
          BEAMINTERACTION::BeamPotentialPair::Create(numnodes_centerline,numnodalvalues);

      newbeaminteractionpair->Init(
          BeamPotentialParamsPtr(),
          firsteleptr,
          secondeleptr);

      newbeaminteractionpair->Setup();

      beam_potential_element_pairs_.push_back(newbeaminteractionpair);
    }
  }

  IO::cout(IO::standard) <<"\n\n************************************************"<<IO::endl;
  IO::cout(IO::standard) << "PID " << GState().GetMyRank() << " currently monitors " <<
      beam_potential_element_pairs_.size() << " beam potential pairs" << IO::endl;
  IO::cout(IO::standard) <<"************************************************"<<IO::endl;

}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    PrintAllBeamPotentialElementPairs( std::ostream& out ) const
{
  out << "\n\nCurrent BeamPotentialElementPairs: ";
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> >::const_iterator iter;
  for (iter=beam_potential_element_pairs_.begin(); iter!=beam_potential_element_pairs_.end(); ++iter)
    (*iter)->Print(out);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    PrintActiveBeamPotentialSet( std::ostream& out ) const
{
  // Todo
  out << "\n    Active BeamToBeam Potential Set (PID " << GState().GetMyRank() <<
      "):-----------------------------------------\n";
  out << "    ID1            ID2              T xi       eta      angle    gap         force\n";

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> >::const_iterator iter;
  for (iter=beam_potential_element_pairs_.begin(); iter!=beam_potential_element_pairs_.end(); ++iter)
    (*iter)->PrintSummaryOneLinePerActiveSegmentPair(out);

  out << std::endl;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    GetBeamPotentialConditionsAppliedToThisElementPair(
    BEAMINTERACTION::BeamPotentialPair const& elementpair,
    std::vector<DRT::Condition*>& conditions_element1,
    std::vector<DRT::Condition*>& conditions_element2) const
{
  // since only the nodes know about their conditions, we need this workaround
  // we assume that a linecharge condition is always applied to the entire physical beam, i.e. it is sufficient to check only one node
  const DRT::Element* ele1 = elementpair.Element1();
  const DRT::Element* ele2 = elementpair.Element2();

  const DRT::Node* const* nodes1;
  const DRT::Node* const* nodes2;
  nodes1 = ele1->Nodes();
  nodes2 = ele2->Nodes();

  dsassert(nodes1 != NULL and nodes2 != NULL, "pointer to nodes is NULL!");
  dsassert(nodes1[0] != NULL and nodes2[0] != NULL, "pointer to nodes is NULL!");

  nodes1[0]->GetCondition("BeamPotentialLineCharge", conditions_element1);
  nodes2[0]->GetCondition("BeamPotentialLineCharge", conditions_element2);

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PrintConsoleWelcomeMessage(
    std::ostream& out) const
{
  // console welcome message
  if (GState().GetMyRank() == 0)
  {
    std::cout << "=============== Beam Potential-Based Interaction ===============" << std::endl;

    switch (BeamPotentialParams().PotentialType())
    {
      case INPAR::BEAMPOTENTIAL::beampot_surf:
      {
        std::cout << "Potential Type:      Surface" << std::endl;
        break;
      }
      case INPAR::BEAMPOTENTIAL::beampot_vol:
      {
        std::cout << "Potential Type:      Volume" << std::endl;
        break;
      }
      default:
        dserror("Potential type not supported!");
    }

    std::cout << "Potential Law:       Phi(r) = ";
    for (unsigned int isummand=0; isummand<BeamPotentialParams().PotentialLawPrefactors().size(); ++isummand)
    {
      if (isummand>0)
        std::cout << " + ";

      std::cout << "(" << BeamPotentialParams().PotentialLawPrefactors().at(isummand) << ") * r^(-"
          << BeamPotentialParams().PotentialLawExponents().at(isummand) << ")";
    }
    std::cout << std::endl;

    std::cout << "================================================================\n" << std::endl;
  }
}
