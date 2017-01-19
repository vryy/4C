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


#include "../drt_structure_new/beaminteraction_submodel_evaluator_generic.H"

#include "../drt_inpar/inpar_binstrategy.H"

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

    //! Returns the type of the current submodel evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_potential; }

    //! @name Derived public BEAMINTERACTION::SUBMODELEVALUATOR::Generic methods
    //! @{
    //! \brief reset submodel specific variables
    virtual void Reset();

    //! \brief Evaluate the current right-hand-side at \f$t_{n+1}\f$
    virtual bool EvaluateForce();

    //! \brief Evaluate the current tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateStiff();

    //! \brief Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateForceStiff();

    virtual void UpdateStepState(const double& timefac_n);

    virtual void PreUpdateStepElement();

    virtual void UpdateStepElement();

    virtual void PostUpdateStepElement();

    //! output routine for model evlaluator
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! reset routine for model evlaluator
    virtual void ResetStepState();

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

    //! type of eles in bins
    INPAR::BINSTRATEGY::BinContent bin_beamcontent_;

     //! interacting pairs of beam elements that might exert forces on each other
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> > beam_potential_element_pairs_;

    //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
    std::map<int, std::set<DRT::Element*> > nearby_elements_map_;

  };

}
}

#endif
