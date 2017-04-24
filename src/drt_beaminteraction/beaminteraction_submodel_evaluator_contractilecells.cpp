/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_contractilecells.cpp

\brief class for managing contractile cell to network interaction

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#include "beaminteraction_submodel_evaluator_contractilecells.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

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
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_beaminteraction/contractilecells_params.H"
#include "../drt_beaminteraction/crosslinker_node.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "beaminteraction_calc_utils.H"
#include "beaminteraction_submodel_evaluator_crosslinking.H"
#include "str_model_evaluator_beaminteraction_datastate.H"

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::ContractileCells() :
    sm_crosslinkink_ptr(Teuchos::null),
    contractilecells_params_ptr_( Teuchos::null ),
    bin_beamcontent_( BINSTRATEGY::UTILS::Beam )
{
  // clear stl stuff
  nearby_elements_map_.clear();
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

}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateForce()
{
  CheckInitSetup();

  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateStiff()
{
  CheckInitSetup();


  return true;
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::EvaluateForceStiff()
{
  CheckInitSetup();

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

}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::PostUpdateStepElement()
{
  CheckInitSetup();

  std::set<int> spherebingids;
  int spheregid = -1;
  for( int i = 0; i < EleTypeMapExtractorPtr()->SphereMap()->NumMyElements(); ++i )
  {
    spheregid =  EleTypeMapExtractorPtr()->SphereMap()->GID(i);
    spherebingids.insert( BeamInteractionDataStatePtr()->GetRowEleToBinSet(spheregid).begin(),
                          BeamInteractionDataStatePtr()->GetRowEleToBinSet(spheregid).end() );
  }

  sm_crosslinkink_ptr->UnbindCrosslinkerInBinsAndNeighborhood( spherebingids, false );

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
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  // empty
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
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::FindAndStoreNeighboringElements()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::FindAndStoreNeighboringElements");

  // loop over all row elements
  const int numroweles = DiscretPtr()->NumMyRowElements();
  for( int rowele_i = 0; rowele_i < numroweles; ++rowele_i )
  {
    const int elegid = DiscretPtr()->ElementRowMap()->GID(rowele_i);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all bind touched by currele
    std::set<int>::const_iterator biniter;
    for( biniter = BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).begin();
        biniter != BeamInteractionDataStatePtr()->GetRowEleToBinSet(elegid).end(); ++biniter )
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds(
          *biniter, loc_neighboring_binIds );

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert( loc_neighboring_binIds.begin(),
          loc_neighboring_binIds.end() );
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds( neighboring_binIds.begin(),
        neighboring_binIds.end() );

    // set of elements that lie in neighboring bins
    std::vector< BINSTRATEGY::UTILS::BinContentType > bc( 1, bin_beamcontent_);
    std::set<DRT::Element*> neighboring_elements;
    BinStrategyPtr()->GetBinContent( neighboring_elements,
        bc, glob_neighboring_binIds );

    nearby_elements_map_[elegid] = neighboring_elements;
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::ContractileCells::UpdateCellsPositionRandomly()
{
  CheckInit();

  DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );
  int const numcells = EleTypeMapExtractorPtr()->SphereMap()->NumMyElements();
  std::vector<double> randpos;
  DRT::Problem::Instance()->Random()->Uni( randpos, 3 * numcells );

  //todo: this is of course not nice, this needs to be done somewhere else
  for( int i = 0; i < numcells; ++i )
  {
    DRT::Element* eleptr = Discret().gElement(EleTypeMapExtractorPtr()->SphereMap()->GID(i) );

    std::vector<int> dofnode  = Discret().Dof(eleptr->Nodes()[0]);

    // random position of cell inside box
    static std::vector< std::vector< double > > Xnew(numcells, std::vector< double >(3) );
    if( GState().GetStepN() % 20 == 0 and GState().GetStepN() != 0)
    {
      for ( int dim = 0; dim < 3; ++dim )
      {
        double edgelength = 5.5;
        double min = 0.5;
        Xnew[i][dim] = min + ( edgelength * randpos[ i + dim ] );
      }
    }

    // loop over all dofs
    for( int dim = 0; dim < 3; ++dim )
    {
      int doflid = BeamInteractionDataStatePtr()->GetMutableDisNp()->Map().LID(dofnode[dim]);
      (*BeamInteractionDataStatePtr()->GetMutableDisNp() )[doflid] = Xnew[i][dim] - eleptr->Nodes()[0]->X()[dim];
    }
 /*   int const elegid = EleTypeMapExtractorPtr()->SphereMap()->GID(i);
    DiscretPtr()->gElement(elegid)->Nodes()[0]->SetPos(Xnew);*/
   }
}


