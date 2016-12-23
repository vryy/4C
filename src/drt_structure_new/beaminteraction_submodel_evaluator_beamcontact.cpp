/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact.cpp

\brief class for submodel beam contact

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/


#include "beaminteraction_submodel_evaluator_beamcontact.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "str_timint_basedataglobalstate.H"
#include "str_utils.H"
#include "../drt_biopolynet/biopolynet_calc_utils.H"
#include "../drt_particle/particle_handler.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include <Epetra_Comm.h>
#include <fenv.h>
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_inpar/inpar_beamcontact.H"

#include "../drt_beamcontact/beam_contact_params.H"
#include "../drt_beamcontact/beam_to_beam_interaction.H"


#include <Epetra_FEVector.h>

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

    std::vector<double> ele1disp;
    BIOPOLYNET::UTILS::GetCurrentElementDis(Discret(), elepairptr->Element1(),
        BeamInteractionDataStatePtr()->GetMutableDisColNp(),ele1disp);

    const DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(elepairptr->Element1());

    unsigned int numcenterlinenodes = beamele1->NumCenterlineNodes();
    unsigned int numnodalvalues =
        beamele1->HermiteCenterlineInterpolation() ? 2 : 1;

    // Todo use LINALG::SerialDenseVector
    Epetra_SerialDenseVector ele1_centerline_dofvec(3*numcenterlinenodes*numnodalvalues);

    // ToDo safety check
    if (numcenterlinenodes!= 2 or numnodalvalues != 1)
      dserror("stop here and do your homework! generalize this method for all types of beams");

    LINALG::Matrix<3,1> nodalpos;

    for (unsigned int inode=0; inode<numcenterlinenodes; ++inode)
    {
      nodalpos.Clear();
      beamele1->GetPosAtXi(nodalpos,-1.0+2*inode/(numcenterlinenodes-1),ele1disp);

      for (unsigned int idim=0; idim<3; ++idim)
        ele1_centerline_dofvec(3*numnodalvalues*inode+idim) = nodalpos(idim);
    }

    std::vector<double> ele2disp;
    BIOPOLYNET::UTILS::GetCurrentElementDis(Discret(), elepairptr->Element2(),
        BeamInteractionDataStatePtr()->GetMutableDisColNp(),ele2disp);

    const DRT::ELEMENTS::Beam3Base* beamele2 =
        dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(elepairptr->Element2());

    Epetra_SerialDenseVector ele2_centerline_dofvec(3*numcenterlinenodes*numnodalvalues);

    for (unsigned int inode=0; inode<numcenterlinenodes; ++inode)
    {
      nodalpos.Clear();
      beamele2->GetPosAtXi(nodalpos,-1.0+2*inode/(numcenterlinenodes-1),ele2disp);

      for (unsigned int idim=0; idim<3; ++idim)
        ele2_centerline_dofvec(3*numnodalvalues*inode+idim) = nodalpos(idim);
    }

    elepairptr->ResetState(ele1_centerline_dofvec,ele2_centerline_dofvec);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateForce()
{
  CheckInitSetup();

  // resulting discrete element force vectors of the two interacting elements
  LINALG::SerialDenseVector ele1force;
  LINALG::SerialDenseVector ele2force;

  // resulting discrete force vectors (centerline DOFs only!) of the two interacting elements
  LINALG::SerialDenseVector ele1force_centerlineDOFs;
  LINALG::SerialDenseVector ele2force_centerlineDOFs;

  Epetra_SerialDenseMatrix dummystiff(0,0);

  // Todo generalize this
  const unsigned int numdof_centerline_ele1 = 6;
  const unsigned int numdof_centerline_ele2 = 6;

  const unsigned int numdof_ele1 = 12;
  const unsigned int numdof_ele2 = 12;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    ele1force_centerlineDOFs.Size(numdof_centerline_ele1);
    ele2force_centerlineDOFs.Size(numdof_centerline_ele2);

    elepairptr->Evaluate(
        &ele1force_centerlineDOFs,
        &ele2force_centerlineDOFs,
        NULL,
        NULL,
        NULL,
        NULL);

    // resize and clear values
    ele1force.Size(numdof_ele1);
    ele2force.Size(numdof_ele2);

    // assemble force vector and stiffness matrix affecting the centerline DoFs only
    // into element force vector and stiffness matrix ('all DoFs' format, as usual)
    BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        ele1force_centerlineDOFs,
        ele2force_centerlineDOFs,
        dummystiff,
        dummystiff,
        dummystiff,
        dummystiff,
        &ele1force,
        &ele2force,
        NULL,
        NULL,
        NULL,
        NULL);

    // Fixme
    ele1force.Scale(-1.0);
    ele2force.Scale(-1.0);

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        ele1force,
        ele2force,
        dummystiff,
        dummystiff,
        dummystiff,
        dummystiff,
        BeamInteractionDataStatePtr()->GetMutableForceNp(),
        Teuchos::null);
  }
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateStiff()
{
  CheckInitSetup();

  // linearizations
  LINALG::SerialDenseMatrix ele11stiff;
  LINALG::SerialDenseMatrix ele12stiff;
  LINALG::SerialDenseMatrix ele21stiff;
  LINALG::SerialDenseMatrix ele22stiff;

  Epetra_SerialDenseVector dummyforce(0);

  // linearizations (centerline DOFs only!)
  LINALG::SerialDenseMatrix ele11stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele12stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele21stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele22stiff_centerlineDOFs;

  // Todo generalize this
  const unsigned int numdof_centerline_ele1 = 6;
  const unsigned int numdof_centerline_ele2 = 6;

  const unsigned int numdof_ele1 = 12;
  const unsigned int numdof_ele2 = 12;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    ele11stiff_centerlineDOFs.Shape(numdof_centerline_ele1,numdof_centerline_ele1);
    ele12stiff_centerlineDOFs.Shape(numdof_centerline_ele1,numdof_centerline_ele2);
    ele21stiff_centerlineDOFs.Shape(numdof_centerline_ele2,numdof_centerline_ele1);
    ele22stiff_centerlineDOFs.Shape(numdof_centerline_ele2,numdof_centerline_ele2);

    elepairptr->Evaluate(
        NULL,
        NULL,
        &ele11stiff_centerlineDOFs,
        &ele12stiff_centerlineDOFs,
        &ele21stiff_centerlineDOFs,
        &ele22stiff_centerlineDOFs);

    // resize and clear values
    ele11stiff.Shape(numdof_ele1,numdof_ele1);
    ele12stiff.Shape(numdof_ele1,numdof_ele2);
    ele21stiff.Shape(numdof_ele2,numdof_ele1);
    ele22stiff.Shape(numdof_ele2,numdof_ele2);

    // assemble force vector and stiffness matrix affecting the centerline DoFs only
    // into element force vector and stiffness matrix ('all DoFs' format, as usual)
    BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        dummyforce,
        dummyforce,
        ele11stiff_centerlineDOFs,
        ele12stiff_centerlineDOFs,
        ele21stiff_centerlineDOFs,
        ele22stiff_centerlineDOFs,
        NULL,
        NULL,
        &ele11stiff,
        &ele12stiff,
        &ele21stiff,
        &ele22stiff);

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        dummyforce,
        dummyforce,
        ele11stiff,
        ele12stiff,
        ele21stiff,
        ele22stiff,
        Teuchos::null,
        BeamInteractionDataStatePtr()->GetMutableStiff());
  }
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact::EvaluateForceStiff()
{
  CheckInitSetup();

  // resulting discrete element force vectors of the two interacting elements
  LINALG::SerialDenseVector ele1force;
  LINALG::SerialDenseVector ele2force;

  // linearizations
  LINALG::SerialDenseMatrix ele11stiff;
  LINALG::SerialDenseMatrix ele12stiff;
  LINALG::SerialDenseMatrix ele21stiff;
  LINALG::SerialDenseMatrix ele22stiff;


  // resulting discrete force vectors (centerline DOFs only!) of the two interacting elements
  LINALG::SerialDenseVector ele1force_centerlineDOFs;
  LINALG::SerialDenseVector ele2force_centerlineDOFs;

  // linearizations (centerline DOFs only!)
  LINALG::SerialDenseMatrix ele11stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele12stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele21stiff_centerlineDOFs;
  LINALG::SerialDenseMatrix ele22stiff_centerlineDOFs;

  // Todo generalize this
  const unsigned int numdof_centerline_ele1 = 6;
  const unsigned int numdof_centerline_ele2 = 6;

  const unsigned int numdof_ele1 = 12;
  const unsigned int numdof_ele2 = 12;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> elepairptr = *iter;

    ele1force_centerlineDOFs.Size(numdof_centerline_ele1);
    ele2force_centerlineDOFs.Size(numdof_centerline_ele2);

    ele11stiff_centerlineDOFs.Shape(numdof_centerline_ele1,numdof_centerline_ele1);
    ele12stiff_centerlineDOFs.Shape(numdof_centerline_ele1,numdof_centerline_ele2);
    ele21stiff_centerlineDOFs.Shape(numdof_centerline_ele2,numdof_centerline_ele1);
    ele22stiff_centerlineDOFs.Shape(numdof_centerline_ele2,numdof_centerline_ele2);

    elepairptr->Evaluate(
        &ele1force_centerlineDOFs,
        &ele2force_centerlineDOFs,
        &ele11stiff_centerlineDOFs,
        &ele12stiff_centerlineDOFs,
        &ele21stiff_centerlineDOFs,
        &ele22stiff_centerlineDOFs);

    // resize and clear values
    ele1force.Size(numdof_ele1);
    ele2force.Size(numdof_ele2);

    ele11stiff.Shape(numdof_ele1,numdof_ele1);
    ele12stiff.Shape(numdof_ele1,numdof_ele2);
    ele21stiff.Shape(numdof_ele2,numdof_ele1);
    ele22stiff.Shape(numdof_ele2,numdof_ele2);

    // assemble force vector and stiffness matrix affecting the centerline DoFs only
    // into element force vector and stiffness matrix ('all DoFs' format, as usual)
    BIOPOLYNET::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        ele1force_centerlineDOFs,
        ele2force_centerlineDOFs,
        ele11stiff_centerlineDOFs,
        ele12stiff_centerlineDOFs,
        ele21stiff_centerlineDOFs,
        ele22stiff_centerlineDOFs,
        &ele1force,
        &ele2force,
        &ele11stiff,
        &ele12stiff,
        &ele21stiff,
        &ele22stiff);

    // Fixme
    ele1force.Scale(-1.0);
    ele2force.Scale(-1.0);

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BIOPOLYNET::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(
        Discret(),
        elepairptr->Element1()->Id(),
        elepairptr->Element2()->Id(),
        ele1force,
        ele2force,
        ele11stiff,
        ele12stiff,
        ele21stiff,
        ele22stiff,
        BeamInteractionDataStatePtr()->GetMutableForceNp(),
        BeamInteractionDataStatePtr()->GetMutableStiff());
  }
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

    std::set<DRT::Element*>::const_iterator secondeleiter;
    for (secondeleiter=nearbyeleiter->second.begin(); secondeleiter!=nearbyeleiter->second.end(); ++secondeleiter)
    {
      const DRT::Element* secondeleptr = *secondeleiter;

      const unsigned int numnodes = 2;
      const unsigned int numnodalvalues =1;

      Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> newbeaminteractionpair =
          BEAMINTERACTION::BeamToBeamInteraction::Create(numnodes,numnodalvalues);

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
  out << "\n    Active BeamToBeam Contact Set:------------------------------------------------\n";
  out << "    ID1            ID2              T xi       eta      angle    gap         force\n";

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> >::const_iterator iter;
  for (iter=BTB_contact_elepairs_.begin(); iter!=BTB_contact_elepairs_.end(); ++iter)
    (*iter)->PrintSummaryOneLinePerActiveSegmentPair(out);

  out << std::endl;
}
