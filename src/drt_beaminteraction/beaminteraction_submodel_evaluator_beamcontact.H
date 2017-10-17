/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact.H

\brief class for submodel beam contact

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_H_


#include "../drt_binstrategy/binning_strategy_utils.H"
#include "beaminteraction_submodel_evaluator_generic.H"

//forward declaration ...
class RuntimeVtpWriter;

namespace DRT {
  class Element;
}
namespace BEAMINTERACTION {
class BeamContactParams;
class BeamContactPair;

namespace SUBMODELEVALUATOR {

  class BeamContact : public Generic
  {
  public:
    //! constructor
    BeamContact();

    //! destructor
    virtual ~BeamContact() {};

    //! setup class variables
    virtual void Setup();

    //! derived
    virtual void PostSetup();

    //! Returns the type of the current submodel evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_beamcontact; }

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
    //!@name routines that are not derived and handle beam to beam contact
    //! @{
    /// print
    void PrintAllBeamContactElementPairs(std::ostream& out) const;

    /// print
    void PrintActiveBeamContactSet(std::ostream& out) const;

    //! @}

  private:
    //!@name routines that handle visualization output for contact forces
    //! @{

    //! writes output for contact forces in VTP format
    void WriteOutputRuntimeVtpBeamContact() const;

    //! init output for contact forces in VTP format
    void InitOutputRuntimeVtpBeamContact();

    //! @}

    //!@name routines that are not derived and handle beam to beam contact
    //! @{

    // init element types considered for beam to ? contact
    void InitElementTypesConsideredForContact();

    /// get neighbouring eles in discret
    void FindAndStoreNeighboringElements();

    /// exclude certain neighbors from interaction evaluation
    void SelectElesToBeConsideredForContactEvaluation(
        DRT::Element* currele,
        std::set<DRT::Element*>& neighbors
    ) const;

    /// create instances of class BeamContactPair that will be evaluated
    //  to get force and stiffness contributions from beam interactions
    void CreateBeamContactElementPairs();

    //! @}

  private:
    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_ptr_;

    //! element types considered for beam to ? contact
    std::vector< BINSTRATEGY::UTILS::BinContentType > contactelementtypes_;

     //! interacting pairs of beam elements that might exert forces on each other
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair> > contact_elepairs_;

    //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
    std::map<int, std::set<DRT::Element*> > nearby_elements_map_;

    //! runtime vtp writer for visualization of contact forces
    Teuchos::RCP<RuntimeVtpWriter> vtp_writer_ptr_;

  };

}
}

#endif
