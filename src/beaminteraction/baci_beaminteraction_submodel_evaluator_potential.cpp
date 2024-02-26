/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief submodel for potential-based beam interactions


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beaminteraction_submodel_evaluator_potential.hpp"

#include "baci_beam3_base.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "baci_beaminteraction_calc_utils.hpp"
#include "baci_beaminteraction_crosslinker_handler.hpp"
#include "baci_beaminteraction_potential_pair.hpp"
#include "baci_beaminteraction_potential_params.hpp"
#include "baci_beaminteraction_potential_runtime_visualization_output_params.hpp"
#include "baci_beaminteraction_str_model_evaluator_datastate.hpp"
#include "baci_io.hpp"
#include "baci_io_pstream.hpp"
#include "baci_io_visualization_manager.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_structure_new_timint_basedataglobalstate.hpp"
#include "baci_utils_exceptions.hpp"

#include <NOX_Solver_Generic.H>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::BeamPotential()
    : beam_potential_params_ptr_(Teuchos::null), beam_potential_element_pairs_(Teuchos::null)
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
  beam_potential_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::BeamPotentialParams());
  BeamPotentialParams().Init(GState().GetTimeN());
  BeamPotentialParams().Setup();

  PrintConsoleWelcomeMessage(std::cout);

  // build runtime visualization writer if desired
  if (BeamPotentialParams().RuntimeOutput()) InitOutputRuntimeBeamPotential();

  // set flag
  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PostSetup()
{
  CheckInitSetup();

  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamPotentialElementPairs();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::InitSubmodelDependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap)
{
  CheckInitSetup();
  // no active influence on other submodels
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::Reset()
{
  CheckInitSetup();

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    std::vector<const DRT::Element*> element_ptr(2);

    element_ptr[0] = elepairptr->Element1();
    element_ptr[1] = elepairptr->Element2();

    // element Dof values relevant for centerline interpolation
    std::vector<std::vector<double>> element_posdofvec_absolutevalues(2);

    for (unsigned int ielement = 0; ielement < 2; ++ielement)
    {
      // extract the Dof values of this element from displacement vector
      BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(Discret(), element_ptr[ielement],
          BeamInteractionDataStatePtr()->GetDisColNp(), element_posdofvec_absolutevalues[ielement]);
    }

    // update the Dof values in the interaction element pair object
    elepairptr->ResetState(GState().GetTimeNp(), element_posdofvec_absolutevalues[0],
        element_posdofvec_absolutevalues[1]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateForce()
{
  CheckInitSetup();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateForce");

  // resulting discrete element force vectors of the two interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce_centerlineDOFs(2);

  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero force values returned which need assembly?
  bool pair_is_active = false;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr, conditions_element1, conditions_element2);

    for (auto& k : conditions_element1)
    {
      int npotlaw1 = *k->Get<int>("potlaw");

      for (auto& j : conditions_element2)
      {
        int npotlaw2 = *j->Get<int>("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(k);
          currconds.push_back(j);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");

          pair_is_active = elepairptr->Evaluate(&(eleforce_centerlineDOFs[0]),
              &(eleforce_centerlineDOFs[1]), nullptr, nullptr, nullptr, nullptr, currconds,
              BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1 - 1),
              BeamPotentialParams().PotentialLawExponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->Element1()->Id();
            elegids[1] = elepairptr->Element2()->Id();

            // assemble force vector affecting the centerline DoFs only
            // into element force vector ('all DoFs' format, as usual)
            BEAMINTERACTION::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
                Discret(), elegids, eleforce_centerlineDOFs, dummystiff, &eleforce, nullptr);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(),
                elegids, eleforce, dummystiff, BeamInteractionDataStatePtr()->GetForceNp(),
                Teuchos::null);
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

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateStiff");

  // linearizations
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  std::vector<CORE::LINALG::SerialDenseVector> dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr, conditions_element1, conditions_element2);

    for (unsigned int k = 0; k < conditions_element1.size(); ++k)
    {
      int npotlaw1 = *conditions_element1[k]->Get<int>("potlaw");

      for (unsigned int j = 0; j < conditions_element2.size(); ++j)
      {
        int npotlaw2 = *conditions_element2[j]->Get<int>("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active = elepairptr->Evaluate(nullptr, nullptr, &(elestiff_centerlineDOFs[0][0]),
              &(elestiff_centerlineDOFs[0][1]), &(elestiff_centerlineDOFs[1][0]),
              &(elestiff_centerlineDOFs[1][1]), currconds,
              BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1 - 1),
              BeamPotentialParams().PotentialLawExponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->Element1()->Id();
            elegids[1] = elepairptr->Element2()->Id();

            // assemble stiffness matrix affecting the centerline DoFs only
            // into element stiffness matrix ('all DoFs' format, as usual)
            BEAMINTERACTION::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(
                Discret(), elegids, dummyforce, elestiff_centerlineDOFs, nullptr, &elestiff);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(),
                elegids, dummyforce, elestiff, Teuchos::null,
                BeamInteractionDataStatePtr()->GetStiff());
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

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::EvaluateForceStiff");

  // resulting discrete element force vectors of the two interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce_centerlineDOFs(2);

  // linearizations
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    elegids[0] = elepairptr->Element1()->Id();
    elegids[1] = elepairptr->Element2()->Id();

    // conditions applied to the elements of this pair
    std::vector<DRT::Condition*> conditions_element1;
    std::vector<DRT::Condition*> conditions_element2;
    GetBeamPotentialConditionsAppliedToThisElementPair(
        *elepairptr, conditions_element1, conditions_element2);

    for (unsigned int k = 0; k < conditions_element1.size(); ++k)
    {
      int npotlaw1 = *conditions_element1[k]->Get<int>("potlaw");

      for (unsigned int j = 0; j < conditions_element2.size(); ++j)
      {
        int npotlaw2 = *conditions_element2[j]->Get<int>("potlaw");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<DRT::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)BeamPotentialParams().PotentialLawPrefactors().size())
            dserror(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active =
              elepairptr->Evaluate(&(eleforce_centerlineDOFs[0]), &(eleforce_centerlineDOFs[1]),
                  &(elestiff_centerlineDOFs[0][0]), &(elestiff_centerlineDOFs[0][1]),
                  &(elestiff_centerlineDOFs[1][0]), &(elestiff_centerlineDOFs[1][1]), currconds,
                  BeamPotentialParams().PotentialLawPrefactors().at(npotlaw1 - 1),
                  BeamPotentialParams().PotentialLawExponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->Element1()->Id();
            elegids[1] = elepairptr->Element2()->Id();

            // assemble force vector and stiffness matrix affecting the centerline DoFs only
            // into element force vector and stiffness matrix ('all DoFs' format, as usual)
            BEAMINTERACTION::UTILS::AssembleCenterlineDofForceStiffIntoElementForceStiff(Discret(),
                elegids, eleforce_centerlineDOFs, elestiff_centerlineDOFs, &eleforce, &elestiff);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(),
                elegids, eleforce, elestiff, BeamInteractionDataStatePtr()->GetForceNp(),
                BeamInteractionDataStatePtr()->GetStiff());
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
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::UpdateStepState(const double& timefac_n)
{
  CheckInitSetup();

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PreUpdateStepElement(bool beam_redist)
{
  CheckInitSetup();
  // not repartition of binning discretization necessary
  return false;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::UpdateStepElement(bool repartition_was_done)
{
  CheckInitSetup();

  /* Fixme
   * writing vtk output needs to be done BEFORE updating (and thus clearing
   * element pairs)
   * move this to RuntimeOutputStepState as soon as we keep element pairs
   * from previous time step */
  if (visualization_manager_ != Teuchos::null and
      GState().GetStepNp() % BeamPotentialParams()
                                 .GetBeamPotentialVisualizationOutputParams()
                                 ->OutputIntervalInSteps() ==
          0)
  {
    WriteTimeStepOutputRuntimeBeamPotential();
  }

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
std::map<STR::EnergyType, double> BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::GetEnergy()
    const
{
  CheckInitSetup();

  std::map<STR::EnergyType, double> beam_interaction_potential;

  for (auto& elepairptr : beam_potential_element_pairs_)
  {
    beam_interaction_potential[STR::beam_interaction_potential] += elepairptr->GetEnergy();
  }

  return beam_interaction_potential;
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
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::RuntimeOutputStepState() const
{
  CheckInitSetup();
  // nothing to do (so far)
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::ResetStepState() { CheckInitSetup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::WriteRestart(
    IO::DiscretizationWriter& ia_writer, IO::DiscretizationWriter& bin_writer) const
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PreReadRestart()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::ReadRestart(
    IO::DiscretizationReader& ia_reader, IO::DiscretizationReader& bin_reader)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PostReadRestart()
{
  CheckInitSetup();
  nearby_elements_map_.clear();
  FindAndStoreNeighboringElements();
  CreateBeamPotentialElementPairs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::RunPostIterate(
    const ::NOX::Solver::Generic& solver)
{
  CheckInitSetup();

  if (visualization_manager_ != Teuchos::null and
      BeamPotentialParams().GetBeamPotentialVisualizationOutputParams()->OutputEveryIteration())
  {
    WriteIterationOutputRuntimeBeamPotential(solver.getNumIterations());
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::AddBinsToBinColMap(std::set<int>& colbins)
{
  CheckInitSetup();
  // nothing to do
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    AddBinsWithRelevantContentForIaDiscretColMap(std::set<int>& colbins) const
{
  CheckInitSetup();
  // nothing to do
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::GetHalfInteractionDistance(
    double& half_interaction_distance)
{
  CheckInitSetup();

  if (BeamPotentialParams().CutoffRadius() > 0.0)
  {
    half_interaction_distance = 0.5 * BeamPotentialParams().CutoffRadius();

    if (GState().GetMyRank() == 0)
      IO::cout(IO::verbose) << " beam potential half interaction distance "
                            << half_interaction_distance << IO::endl;
  }
  else
  {
    dserror(
        "You have to set a cutoff radius for beam-to-? potential-based interactions in order "
        "to use REPARTITIONSTRATEGY = Adaptive!");
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::FindAndStoreNeighboringElements()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::FindAndStoreNeighboringElements");

  // loop over all row elements
  int const numroweles = EleTypeMapExtractorPtr()->BeamMap()->NumMyElements();
  for (int rowele_i = 0; rowele_i < numroweles; ++rowele_i)
  {
    int const elegid = EleTypeMapExtractorPtr()->BeamMap()->GID(rowele_i);
    DRT::Element* currele = DiscretPtr()->gElement(elegid);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all bind touched by currele
    std::set<int>::const_iterator biniter;
    for (biniter = BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).begin();
         biniter != BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).end(); ++biniter)
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds(*biniter, loc_neighboring_binIds);

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert(loc_neighboring_binIds.begin(), loc_neighboring_binIds.end());
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds(neighboring_binIds.begin(), neighboring_binIds.end());

    // set of elements that lie in neighboring bins
    std::set<DRT::Element*> neighboring_elements;
    std::vector<BINSTRATEGY::UTILS::BinContentType> bc(2);
    bc[0] = BINSTRATEGY::UTILS::Beam;
    bc[1] = BINSTRATEGY::UTILS::RigidSphere;
    BinStrategyPtr()->GetBinContent(neighboring_elements, bc, glob_neighboring_binIds);

    // sort out elements that should not be considered in contact evaluation
    SelectElesToBeConsideredForPotentialEvaluation(currele, neighboring_elements);

    nearby_elements_map_[elegid] = neighboring_elements;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    SelectElesToBeConsideredForPotentialEvaluation(
        DRT::Element* currele, std::set<DRT::Element*>& neighbors) const
{
  CheckInit();

  // sort out elements that should not be considered in potential evaluation
  std::set<DRT::Element*>::iterator eiter;
  for (eiter = neighbors.begin(); eiter != neighbors.end();)
  {
    bool toerase = false;

    DRT::Element* currneighborele = *eiter;

    // 1) ensure each interaction is only evaluated once (keep in mind that we are
    //    using FEMatrices and FEvectors -> || (*eiter)->Owner() != myrank not necessary)
    if (not(currele->Id() < currneighborele->Id()))
    {
      toerase = true;
    }

    // 2) exclude "self-interaction", i.e. a pair of elements on the same physical beam
    //    TODO introduce flag for self-interaction in input file
    else
    {
      // get the conditions applied to both elements of the pair and decide whether they need to be
      // evaluated
      std::vector<DRT::Condition*> conds1, conds2;

      // since only the nodes know about their conditions, we need this workaround
      // we assume that a linecharge condition is always applied to the entire physical beam, i.e.
      // it is sufficient to check only one node
      DRT::Node** nodes1;
      DRT::Node** nodes2;
      nodes1 = currele->Nodes();
      nodes2 = currneighborele->Nodes();

      dsassert(nodes1 != nullptr and nodes2 != nullptr, "pointer to nodes is nullptr!");
      dsassert(nodes1[0] != nullptr and nodes2[0] != nullptr, "pointer to nodes is nullptr!");

      nodes1[0]->GetCondition("BeamPotentialLineCharge", conds1);

      // get correct condition for beam or rigid sphere element
      if (BEAMINTERACTION::BeamElement(*currneighborele))
        nodes2[0]->GetCondition("BeamPotentialLineCharge", conds2);
      else if (BEAMINTERACTION::RigidsphereElement(*currneighborele))
        nodes2[0]->GetCondition("RigidspherePotentialPointCharge", conds2);
      else
        dserror(
            "Only beam-to-beampotential or beam-to-sphere -based interaction is implemented yet. "
            "No other types of elements allowed!");

      // validinteraction == true includes: both eles "loaded" by a charge condition of same
      // potential law
      bool validinteraction = false;

      for (unsigned int i = 0; i < conds1.size(); ++i)
      {
        int npotlaw1 = *conds1[i]->Get<int>("potlaw");

        for (unsigned int j = 0; j < conds2.size(); ++j)
        {
          int npotlaw2 = *conds2[j]->Get<int>("potlaw");

          // here, we also exclude "self-interaction", i.e. a pair of elements on the same physical
          // beam
          // TODO introduce flag for self-interaction in input file
          if (conds1[i] != conds2[j] and npotlaw1 == npotlaw2) validinteraction = true;
        }
      }

      if (not validinteraction) toerase = true;
    }


    if (toerase)
      neighbors.erase(eiter++);
    else
      ++eiter;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::CreateBeamPotentialElementPairs()
{
  // Todo maybe keep existing pairs and reuse them ?
  beam_potential_element_pairs_.clear();

  std::map<int, std::set<DRT::Element*>>::const_iterator nearbyeleiter;

  for (nearbyeleiter = nearby_elements_map_.begin(); nearbyeleiter != nearby_elements_map_.end();
       ++nearbyeleiter)
  {
    const int elegid = nearbyeleiter->first;
    std::vector<DRT::Element const*> ele_ptrs(2);
    ele_ptrs[0] = DiscretPtr()->gElement(elegid);

    std::set<DRT::Element*>::const_iterator secondeleiter;
    for (secondeleiter = nearbyeleiter->second.begin();
         secondeleiter != nearbyeleiter->second.end(); ++secondeleiter)
    {
      ele_ptrs[1] = *secondeleiter;

      Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> newbeaminteractionpair =
          BEAMINTERACTION::BeamPotentialPair::Create(ele_ptrs, BeamPotentialParams());

      newbeaminteractionpair->Init(BeamPotentialParamsPtr(), ele_ptrs[0], ele_ptrs[1]);

      newbeaminteractionpair->Setup();

      beam_potential_element_pairs_.push_back(newbeaminteractionpair);
    }
  }

  if (static_cast<int>(beam_potential_element_pairs_.size()) > 0)
  {
    IO::cout(IO::standard) << "PID " << std::setw(2) << std::right << GState().GetMyRank()
                           << " currently monitors " << std::setw(5) << std::right
                           << beam_potential_element_pairs_.size() << " beam potential pairs"
                           << IO::endl;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PrintAllBeamPotentialElementPairs(
    std::ostream& out) const
{
  out << "\n\nCurrent BeamPotentialElementPairs: ";
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
    (*iter)->Print(out);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::PrintActiveBeamPotentialSet(
    std::ostream& out) const
{
  // Todo
  out << "\n    Active BeamToBeam Potential Set (PID " << GState().GetMyRank()
      << "):-----------------------------------------\n";
  out << "    ID1            ID2              T xi       eta      angle    gap         force\n";

  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
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
  // we assume that a linecharge condition is always applied to the entire physical beam, i.e. it is
  // sufficient to check only one node
  const DRT::Element* ele1 = elementpair.Element1();
  const DRT::Element* ele2 = elementpair.Element2();

  const DRT::Node* const* nodes1;
  const DRT::Node* const* nodes2;
  nodes1 = ele1->Nodes();
  nodes2 = ele2->Nodes();

  dsassert(nodes1 != nullptr and nodes2 != nullptr, "pointer to nodes is nullptr!");
  dsassert(nodes1[0] != nullptr and nodes2[0] != nullptr, "pointer to nodes is nullptr!");

  nodes1[0]->GetCondition("BeamPotentialLineCharge", conditions_element1);

  // get correct condition for beam or rigid sphere element
  if (BEAMINTERACTION::BeamElement(*ele2))
    nodes2[0]->GetCondition("BeamPotentialLineCharge", conditions_element2);
  else if (BEAMINTERACTION::RigidsphereElement(*ele2))
    nodes2[0]->GetCondition("RigidspherePotentialPointCharge", conditions_element2);
  else
    dserror(
        "Only beam-to-beam or beam-to-sphere potential-based interaction is implemented yet. "
        "No other types of elements allowed!");
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
    for (unsigned int isummand = 0;
         isummand < BeamPotentialParams().PotentialLawPrefactors().size(); ++isummand)
    {
      if (isummand > 0) std::cout << " + ";

      std::cout << "(" << BeamPotentialParams().PotentialLawPrefactors().at(isummand) << ") * r^(-"
                << BeamPotentialParams().PotentialLawExponents().at(isummand) << ")";
    }
    std::cout << std::endl;

    std::cout << "================================================================\n" << std::endl;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::InitOutputRuntimeBeamPotential()
{
  CheckInit();

  visualization_manager_ =
      Teuchos::rcp(new IO::VisualizationManager(BeamPotentialParams()
                                                    .GetBeamPotentialVisualizationOutputParams()
                                                    ->GetVisualizationParameters(),
          Discret().Comm(), "beam-potential"));
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::WriteTimeStepOutputRuntimeBeamPotential()
    const
{
  CheckInitSetup();

  auto [output_time, output_step] =
      IO::GetTimeAndTimeStepIndexForOutput(BeamPotentialParams()
                                               .GetBeamPotentialVisualizationOutputParams()
                                               ->GetVisualizationParameters(),
          GState().GetTimeN(), GState().GetStepN());
  WriteOutputRuntimeBeamPotential(output_step, output_time);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::WriteIterationOutputRuntimeBeamPotential(
    int iteration_number) const
{
  CheckInitSetup();

  auto [output_time, output_step] =
      IO::GetTimeAndTimeStepIndexForOutput(BeamPotentialParams()
                                               .GetBeamPotentialVisualizationOutputParams()
                                               ->GetVisualizationParameters(),
          GState().GetTimeN(), GState().GetStepN(), iteration_number);
  WriteOutputRuntimeBeamPotential(output_step, output_time);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::WriteOutputRuntimeBeamPotential(
    int timestep_number, double time) const
{
  CheckInitSetup();

  const unsigned int num_spatial_dimensions = 3;

  // estimate for number of interacting Gauss points = number of row points for writer object
  unsigned int num_row_points = 0;

  if (BeamPotentialParams()
          .GetBeamPotentialVisualizationOutputParams()
          ->IsWriteForcesMomentsPerElementPair())
  {
    num_row_points = 2 * beam_potential_element_pairs_.size() *
                     BeamPotentialParams().NumberIntegrationSegments() *
                     BeamPotentialParams().NumberGaussPoints();
  }
  else
  {
    // Todo: this won't perfectly work in parallel yet since some communication would be required
    //    if ( GState().GetMyRank() != 0 )
    //      dserror("visualization of resulting forces not implemented in parallel yet!");

    num_row_points = Discret().NumGlobalElements() *
                     BeamPotentialParams().NumberIntegrationSegments() *
                     BeamPotentialParams().NumberGaussPoints();
  }

  // get and prepare storage for point coordinate values
  auto& visualization_data = visualization_manager_->GetVisualizationData();
  std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates();
  point_coordinates.clear();
  point_coordinates.reserve(num_spatial_dimensions * num_row_points);

  // force values: collect data and append to visualization results if desired
  std::vector<double> potential_force_vector(0);
  potential_force_vector.reserve(num_spatial_dimensions * num_row_points);

  // moment values: collect data and append to visualization results if desired
  std::vector<double> potential_moment_vector(0);
  potential_moment_vector.reserve(num_spatial_dimensions * num_row_points);


  // loop over my points and collect the geometry/grid data, i.e. interacting points
  std::vector<CORE::LINALG::Matrix<3, 1, double>> coordinates_ele1_this_pair;
  std::vector<CORE::LINALG::Matrix<3, 1, double>> coordinates_ele2_this_pair;

  std::vector<CORE::LINALG::Matrix<3, 1, double>> potential_forces_ele1_this_pair;
  std::vector<CORE::LINALG::Matrix<3, 1, double>> potential_forces_ele2_this_pair;

  std::vector<CORE::LINALG::Matrix<3, 1, double>> potential_moments_ele1_this_pair;
  std::vector<CORE::LINALG::Matrix<3, 1, double>> potential_moments_ele2_this_pair;

  // loop over contact pairs and retrieve all active contact point coordinates
  std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>>::const_iterator pair_iter;
  for (pair_iter = beam_potential_element_pairs_.begin();
       pair_iter != beam_potential_element_pairs_.end(); ++pair_iter)
  {
    // retrieve data for interacting points of element 1 and element 2
    (*pair_iter)->GetAllInteractingPointCoordsElement1(coordinates_ele1_this_pair);
    (*pair_iter)->GetAllInteractingPointCoordsElement2(coordinates_ele2_this_pair);
    (*pair_iter)->GetForcesAtAllInteractingPointsElement1(potential_forces_ele1_this_pair);
    (*pair_iter)->GetForcesAtAllInteractingPointsElement2(potential_forces_ele2_this_pair);
    (*pair_iter)->GetMomentsAtAllInteractingPointsElement1(potential_moments_ele1_this_pair);
    (*pair_iter)->GetMomentsAtAllInteractingPointsElement2(potential_moments_ele2_this_pair);

    const unsigned int num_interacting_points_per_element =
        (unsigned int)coordinates_ele1_this_pair.size();

    dsassert(num_interacting_points_per_element == (unsigned int)coordinates_ele2_this_pair.size(),
        "number of interacting points on element 1 does not match number of interacting points "
        "on element 2!");

    dsassert(
        num_interacting_points_per_element == (unsigned int)potential_forces_ele1_this_pair.size(),
        "number of interacting points on element 1 does not match number of potential forces!");

    dsassert(
        num_interacting_points_per_element == (unsigned int)potential_forces_ele2_this_pair.size(),
        "number of interacting points on element 2 does not match number of potential forces!");


    for (unsigned int ipoint = 0; ipoint < num_interacting_points_per_element; ++ipoint)
    {
      // ignore point pairs with zero forces
      /* (e.g. if no valid point-to-curve projection in master-slave approach or
       * contribution is neglected on element pair level due to cutoff value) */
      if (potential_forces_ele1_this_pair[ipoint].Norm2() < 1e-16 and
          potential_forces_ele2_this_pair[ipoint].Norm2() < 1e-16 and
          potential_moments_ele1_this_pair[ipoint].Norm2() < 1e-16 and
          potential_moments_ele2_this_pair[ipoint].Norm2() < 1e-16)
      {
        continue;
      }


      // this is easier, since data is computed and stored in this 'element-pairwise' format
      if (BeamPotentialParams()
              .GetBeamPotentialVisualizationOutputParams()
              ->IsWriteForcesMomentsPerElementPair())
      {
        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele1_this_pair[ipoint](idim));

          potential_force_vector.push_back(potential_forces_ele1_this_pair[ipoint](idim));
          potential_moment_vector.push_back(potential_moments_ele1_this_pair[ipoint](idim));
        }

        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele2_this_pair[ipoint](idim));

          potential_force_vector.push_back(potential_forces_ele2_this_pair[ipoint](idim));
          potential_moment_vector.push_back(potential_moments_ele2_this_pair[ipoint](idim));
        }
      }
      /* in this case, we need to identify unique Gauss points based on their coordinate values and
       * compute resulting force/moment at this point by summation of contributions from all element
       * pairs */
      else
      {
        // interacting point on first element
        std::vector<double>::iterator xcoord_iter = point_coordinates.begin();

        // try to find data point with identical coordinates
        while (point_coordinates.size() >= 3 and xcoord_iter != point_coordinates.end() - 2)
        {
          // find identical x-coordinate value
          xcoord_iter = std::find(
              xcoord_iter, point_coordinates.end() - 2, coordinates_ele1_this_pair[ipoint](0));

          // check whether we've reached the end -> no match
          if (xcoord_iter == point_coordinates.end() - 2)
          {
            break;
          }
          // we have a match -> check whether also y- and z-coordinate value are identical
          else if ((xcoord_iter - point_coordinates.begin()) % 3 == 0 and
                   *(xcoord_iter + 1) == coordinates_ele1_this_pair[ipoint](1) and
                   *(xcoord_iter + 2) == coordinates_ele1_this_pair[ipoint](2))
          {
            int offset = xcoord_iter - point_coordinates.begin();

            for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
            {
              *(potential_force_vector.begin() + offset + idim) +=
                  potential_forces_ele1_this_pair[ipoint](idim);
              *(potential_moment_vector.begin() + offset + idim) +=
                  potential_moments_ele1_this_pair[ipoint](idim);
            }

            break;
          }
          // we have a matching value but it's not a point with the identical (x,y,z) coordinates
          else
          {
            xcoord_iter++;
          }
        }

        // add as a new point if not found above
        if (xcoord_iter == point_coordinates.end() - 2 or point_coordinates.empty())
        {
          for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
          {
            point_coordinates.push_back(coordinates_ele1_this_pair[ipoint](idim));

            potential_force_vector.push_back(potential_forces_ele1_this_pair[ipoint](idim));
            potential_moment_vector.push_back(potential_moments_ele1_this_pair[ipoint](idim));
          }
        }


        // interacting point on second element
        xcoord_iter = point_coordinates.begin();

        // try to find data point with identical coordinates
        while (xcoord_iter != point_coordinates.end() - 2)
        {
          // find identical x-coordinate value
          xcoord_iter = std::find(
              xcoord_iter, point_coordinates.end() - 2, coordinates_ele2_this_pair[ipoint](0));

          // check whether we've reached the end -> no match
          if (xcoord_iter == point_coordinates.end() - 2)
          {
            break;
          }
          // we have a match -> check whether also y- and z-coordinate value are identical
          else if ((xcoord_iter - point_coordinates.begin()) % 3 == 0 and
                   *(xcoord_iter + 1) == coordinates_ele2_this_pair[ipoint](1) and
                   *(xcoord_iter + 2) == coordinates_ele2_this_pair[ipoint](2))
          {
            int offset = xcoord_iter - point_coordinates.begin();

            for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
            {
              *(potential_force_vector.begin() + offset + idim) +=
                  potential_forces_ele2_this_pair[ipoint](idim);
              *(potential_moment_vector.begin() + offset + idim) +=
                  potential_moments_ele2_this_pair[ipoint](idim);
            }

            break;
          }
          // we have a matching value but it's not a point with the identical (x,y,z) coordinates
          else
          {
            xcoord_iter++;
          }
        }

        // add as a new point if not found above
        if (xcoord_iter == point_coordinates.end() - 2)
        {
          for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
          {
            point_coordinates.push_back(coordinates_ele2_this_pair[ipoint](idim));

            potential_force_vector.push_back(potential_forces_ele2_this_pair[ipoint](idim));
            potential_moment_vector.push_back(potential_moments_ele2_this_pair[ipoint](idim));
          }
        }
      }
    }
  }


  // append all desired output data to the writer object's storage
  if (BeamPotentialParams().GetBeamPotentialVisualizationOutputParams()->IsWriteForces())
  {
    visualization_manager_->GetVisualizationData().SetPointDataVector(
        "force", potential_force_vector, num_spatial_dimensions);
  }

  if (BeamPotentialParams().GetBeamPotentialVisualizationOutputParams()->IsWriteMoments())
  {
    visualization_manager_->GetVisualizationData().SetPointDataVector(
        "moment", potential_moment_vector, num_spatial_dimensions);
  }

  // finalize everything and write all required vtk files to filesystem
  visualization_manager_->WriteToDisk(time, timestep_number);
}

BACI_NAMESPACE_CLOSE
