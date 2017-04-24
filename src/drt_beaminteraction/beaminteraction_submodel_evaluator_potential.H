/*-----------------------------------------------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_potential.H

\brief submodel for potential-based beam interactions

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_POTENTIAL_H_


#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_generic.H"
#include "../drt_particle/binning_strategy_utils.H"

//forward declaration ...
namespace DRT {
  class Element;
  class Condition;
}
namespace BEAMINTERACTION {
class BeamPotentialParams;
class BeamPotentialPair;

namespace SUBMODELEVALUATOR {

  class BeamPotential : public Generic
  {
  public:
    //! constructor
    BeamPotential();

    //! destructor
    virtual ~BeamPotential() {};

    //! setup class variables
    virtual void Setup();

    //! derived
    virtual void PostSetup();

    //! Returns the type of the current submodel evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_potential; }

    //! @name Derived public BEAMINTERACTION::SUBMODELEVALUATOR::Generic methods
    //! @{
    //! \brief reset submodel specific variables
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
    virtual void RuntimeOutputStepState() const;

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

  protected:
    //!@name routines that are not derived and handle beam potential-based interactions
    //! @{
    /// print
    void PrintAllBeamPotentialElementPairs(std::ostream& out) const;

    /// print
    void PrintActiveBeamPotentialSet(std::ostream& out) const;

    //! @}

  private:

    inline BEAMINTERACTION::BeamPotentialParams const& BeamPotentialParams() const
    {
      CheckInit();
      return *beam_potential_params_ptr_;
    }

    inline BEAMINTERACTION::BeamPotentialParams& BeamPotentialParams()
    {
      CheckInit();
      return *beam_potential_params_ptr_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> BeamPotentialParamsPtr() const
    {
      CheckInit();
      return beam_potential_params_ptr_;
    }

    //!@name routines that are not derived and handle beam potential-based interactions
    //! @{
    /// get neighbouring eles in discret
    void FindAndStoreNeighboringElements();

    /// exclude certain neighbors from interaction evaluation
    void SelectElesToBeConsideredForPotentialEvaluation(
        DRT::Element* currele,
        std::set<DRT::Element*>& neighbors
    ) const;

    /// create instances of class BeamContactPair that will be evaluated
    //  to get force and stiffness contributions from beam interactions
    void CreateBeamPotentialElementPairs();

    void GetBeamPotentialConditionsAppliedToThisElementPair(
        BEAMINTERACTION::BeamPotentialPair const& elementpair,
        std::vector<DRT::Condition*>& conditions_element1,
        std::vector<DRT::Condition*>& conditions_element2) const;

    //! @}

    /** \brief print this beam potential-based element pair to screen
     *
     *  \author grill */
    void PrintConsoleWelcomeMessage(std::ostream& out) const;

  private:
    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> beam_potential_params_ptr_;

    //! type of eles in bins  // Todo kept line for future improvement
//    BINSTRATEGY::UTILS::BinContentType bin_beamcontent_;

     //! interacting pairs of beam elements that might exert forces on each other
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> > beam_potential_element_pairs_;

    //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
    std::map<int, std::set<DRT::Element*> > nearby_elements_map_;

  };

}
}

#endif
