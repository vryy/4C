/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_contractilecells.H

\brief class for managing contractile cell to network interaction

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_CONTRACTILECELLS_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_CONTRACTILECELLS_H_

#include "../drt_particle/binning_strategy_utils.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include <Epetra_MpiComm.h>
#include "beaminteraction_submodel_evaluator_generic.H"


// forward declaration
namespace DRT {
  class Element;
  class Node;
}
namespace BEAMINTERACTION {
class ContractileCellsParams;

namespace SUBMODELEVALUATOR {

  class ContractileCells : public Generic
  {
  public:

    //! constructor
    ContractileCells();

    //! destructor
    virtual ~ContractileCells() {};

    //! setup class variables
    virtual void Setup();

    //! derived
    virtual void PostSetup();

    //! Returns the type of the current model evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_contractilecells; }

    //! derived
    virtual void Reset();

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();

    //! derived
    virtual void UpdateStepState(const double& timefac_n);

    //! derived
    virtual void PreUpdateStepElement();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    virtual void PostUpdateStepElement();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual void ResetStepState();

    //! derived
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart( IO::DiscretizationReader& ioreader );

    //! derived
    virtual void PostReadRestart();

    //! derived
    virtual void InitSubmodelDependencies(Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap);

    //! derived
    virtual void AddBinsToBinColMap( std::set< int >& colbins);

    //! derived
    virtual void AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins) const;

    //! @}

private:

    //! @name submodel specific, not derived member functions
    //! @{

    /// get neighbouring eles in discret
    virtual void FindAndStoreNeighboringElements();

    /// update reference position of cells randomly inside a volume
    virtual void UpdateCellsPositionRandomly();



    //! @}

  private:

    //! cell discretization
    Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking> sm_crosslinkink_ptr;

    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::ContractileCellsParams> contractilecells_params_ptr_;

    //! type of eles in bins
    BINSTRATEGY::UTILS::BinContentType bin_beamcontent_;

    //! mapping ele (elegid) to set of spatially proximal eles (pointer to elements)
    std::map<int, std::set<DRT::Element*> > nearby_elements_map_;

  };

}
}

#endif
