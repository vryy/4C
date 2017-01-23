/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_generic.H

\brief Generic class for all beaminteraction submodel evaluators.

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_GENERIC_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_GENERIC_H_

#include "../drt_inpar/inpar_beaminteraction.H"
#include "str_model_evaluator_beaminteraction.H"

// forward declaration
namespace PARTICLE {
  class ParticleHandler;
}
namespace IO {
  class DiscretizationWriter;
  class DiscretizationReader;
}
namespace DRT {
  class Discretization;
}
namespace GEO{
namespace MESHFREE{
  class BoundingBox;
}
}
namespace STR {
namespace TIMINT{
  class BaseDataGlobalState;
}
namespace MODELEVALUATOR {
  class BeamInteractionDataState;
}
}
namespace BEAMINTERACTION {
namespace SUBMODELEVALUATOR {
  class Crosslinking;

  /*! \brief This is the abstract base class of all submodel evaluators for a beaminteraction problem
   *
   *  This class summarizes the functionality which all submodel evaluators share
   *  and/or have to implement. Look in the derived classes for examples. A minimal
   *  example can be found at \ref BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking.
   */
  class Generic
  {
  public:
    //! constructor
    Generic();

    //! destructor
    virtual ~Generic() {};

    //! initialize the class variables
    virtual void Init(
        Teuchos::RCP<DRT::Discretization> const& ia_discret,
        Teuchos::RCP<DRT::Discretization> const& bindis,
        Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> const& gstate,
        Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> const& ia_gstate_ptr,
        Teuchos::RCP<PARTICLE::ParticleHandler> const& particlehandler,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& periodic_boundingbox
    );

    //! setup class variables
    virtual void Setup() = 0;

  protected:
    //! Returns true, if Init() has been called
    inline const bool& IsInit() const { return isinit_; };

    //! Returns true, if Setup() has been called
    inline const bool& IsSetup() const { return issetup_; };

    //! Checks, if Init() and Setup() have been called
    virtual void CheckInitSetup() const;

    virtual void CheckInit() const;

  public:
    //! Returns the type of the current model evaluator
    virtual INPAR::BEAMINTERACTION::SubModelType Type() const = 0;

    //! \brief reset model specific variables (without jacobian)
    virtual void Reset() = 0;

    //! \brief Evaluate the current right-hand-side at \f$t_{n+1}\f$
    virtual bool EvaluateForce() = 0;

    //! \brief Evaluate the current tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateStiff() = 0;

    //! \brief Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateForceStiff() = 0;

    virtual void UpdateStepState(const double& timefac_n) = 0;

    virtual void PreUpdateStepElement() = 0;

    virtual void UpdateStepElement() = 0;

    virtual void PostUpdateStepElement() = 0;

    //! output routine for model evlaluator
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const = 0;

    //! reset routine for model evlaluator
    virtual void ResetStepState() = 0;

    //! reset routine for model evlaluator
    virtual void InitSubmodelDependencies(Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelvector) = 0;

    //! \brief add subproblem specific contributions to bin col map
    virtual void AddBinsToBinColMap( std::set< int >& colbins) = 0;

    //! \brief add subproblem specific contributions to bin col map
    virtual void AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins) const = 0;

    //! \brief do submodel specific stuff after partitioning
    virtual void PostSetup() = 0;

    //! @name internal accessors
    //! @{
    //! Returns the (structural) discretization
    DRT::Discretization& Discret();
    Teuchos::RCP<DRT::Discretization>& DiscretPtr();
    DRT::Discretization const& Discret() const;

    DRT::Discretization& BinDiscret();
    Teuchos::RCP<DRT::Discretization>& BinDiscretPtr();
    DRT::Discretization const& BinDiscret() const;

    //! Returns the global state data container
    STR::TIMINT::BaseDataGlobalState& GState();
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GStatePtr();
    STR::TIMINT::BaseDataGlobalState const& GState() const;

    //! Returns the global state data container
    STR::MODELEVALUATOR::BeamInteractionDataState& BeamInteractionDataState();
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState>& BeamInteractionDataStatePtr();
    STR::MODELEVALUATOR::BeamInteractionDataState const& BeamInteractionDataState() const;

    PARTICLE::ParticleHandler& ParticleHandler();
    Teuchos::RCP<PARTICLE::ParticleHandler>& ParticleHandlerPtr();
    PARTICLE::ParticleHandler const& ParticleHandler() const;

    GEO::MESHFREE::BoundingBox& PeriodicBoundingBox();
    Teuchos::RCP<GEO::MESHFREE::BoundingBox>& PeriodicBoundingBoxPtr();
    GEO::MESHFREE::BoundingBox const& PeriodicBoundingBox() const;

    //! @}
  protected:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;
  private:

    //! pointer to the interaction discretization
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

    //! pointer to the interaction discretization
    Teuchos::RCP<DRT::Discretization> bindis_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> beaminteractiondatastate_;

    //! particle handler and binning strategy
    Teuchos::RCP<PARTICLE::ParticleHandler> particlehandler_;

    //! periodic bounding box
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> periodic_boundingbox_;

  }; // class Generic

}
}

#endif
