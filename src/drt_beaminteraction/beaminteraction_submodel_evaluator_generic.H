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
namespace BINSTRATEGY{
  class BinningStrategy;
}
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
namespace NOX {
namespace Solver {
  class Generic;
}
}
namespace STR {
namespace TIMINT{
  class BaseDataGlobalState;
  class BaseDataIO;
}
namespace MODELEVALUATOR {
  class BeamInteractionDataState;
}
}
namespace BEAMINTERACTION {
namespace UTILS {
  class MapExtractor;
}
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
        Teuchos::RCP<STR::TIMINT::BaseDataIO> const& gio_ptr,
        Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> const& ia_gstate_ptr,
        Teuchos::RCP<PARTICLE::ParticleHandler> const& particlehandler,
        Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& periodic_boundingbox,
        Teuchos::RCP<BEAMINTERACTION::UTILS::MapExtractor> const& eletypeextractor
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

    //! update state
    virtual void UpdateStepState(const double& timefac_n) = 0;

    //! pre update step element
    virtual bool PreUpdateStepElement( bool beam_redist ) = 0;

    //! update step element
    virtual void UpdateStepElement( bool repartition_was_done ) = 0;

    //! post update step element
    virtual void PostUpdateStepElement() = 0;

    //! get contributions to system energy
    virtual std::map< STR::EnergyType, double > GetEnergy() const = 0;

    //! write submodel specific output
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const = 0;

    //! write submodel specific output during runtime
    virtual void RuntimeOutputStepState() const = 0;

    //! reset routine for model evlaluator
    virtual void ResetStepState() = 0;

    //! \brief write model specific restart
    virtual void WriteRestart(
        IO::DiscretizationWriter & ia_writer,
        IO::DiscretizationWriter & bin_writer) const = 0;

    /*! \brief read model specific restart information
     *
     *  \param ioreader (in) : input reader*/
    virtual void ReadRestart(
        IO::DiscretizationReader & ia_writer,
        IO::DiscretizationReader & bin_writer) = 0;

    //! \brief do stuff post read of model specific restart information
    virtual void PostReadRestart() = 0;

    /*! \brief Executed at the end of the NOX::Solver::Generic::Step() (f.k.a. Iterate()) method
     *
     *  \param solver (in) : reference to the non-linear nox solver object (read-only)
     *
     *  \author grill, hiermeier \date 10/17 */
    virtual void RunPostIterate( const NOX::Solver::Generic& solver ) = 0;

    //! reset routine for model evlaluator
    virtual void InitSubmodelDependencies(Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelvector) = 0;

    //! \brief add subproblem specific contributions to bin col map
    virtual void AddBinsToBinColMap( std::set< int >& colbins) = 0;

    //! \brief add subproblem specific contributions to bin col map
    virtual void AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins) const = 0;

    //! \brief add subproblem specific contributions to bin col map
    virtual void GetHalfInteractionDistance( double & half_interaction_distance ) = 0;

    //! \brief do submodel specific stuff after partitioning
    virtual bool PostPartitionProblem() {return false;};

    //! \brief do submodel specific stuff after setup
    virtual void PostSetup() = 0;

    //! @name internal accessors
    //! @{
    //! Returns the (structural) discretization
    DRT::Discretization& Discret();
    Teuchos::RCP<DRT::Discretization>& DiscretPtr();
    DRT::Discretization const& Discret() const;

    DRT::Discretization& BinDiscret();
    Teuchos::RCP<DRT::Discretization>& BinDiscretPtr();
    Teuchos::RCP<const DRT::Discretization> BinDiscretPtr() const;
    DRT::Discretization const& BinDiscret() const;

    //! Returns the global state data container
    STR::TIMINT::BaseDataGlobalState& GState();
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GStatePtr();
    STR::TIMINT::BaseDataGlobalState const& GState() const;

    //! Returns the global input/output data container
    STR::TIMINT::BaseDataIO& GInOutput();
    STR::TIMINT::BaseDataIO const& GInOutput() const;

    //! Returns the global state data container
    STR::MODELEVALUATOR::BeamInteractionDataState& BeamInteractionDataState();
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState>& BeamInteractionDataStatePtr();
    STR::MODELEVALUATOR::BeamInteractionDataState const& BeamInteractionDataState() const;

    PARTICLE::ParticleHandler& ParticleHandler();
    Teuchos::RCP<PARTICLE::ParticleHandler>& ParticleHandlerPtr();
    PARTICLE::ParticleHandler const& ParticleHandler() const;

    BINSTRATEGY::BinningStrategy& BinStrategy();
    Teuchos::RCP<BINSTRATEGY::BinningStrategy>& BinStrategyPtr();
    BINSTRATEGY::BinningStrategy const& BinStrategy() const;

    GEO::MESHFREE::BoundingBox& PeriodicBoundingBox();
    Teuchos::RCP<GEO::MESHFREE::BoundingBox>& PeriodicBoundingBoxPtr();
    GEO::MESHFREE::BoundingBox const& PeriodicBoundingBox() const;

    BEAMINTERACTION::UTILS::MapExtractor& EleTypeMapExtractor();
    Teuchos::RCP<BEAMINTERACTION::UTILS::MapExtractor>& EleTypeMapExtractorPtr();
    BEAMINTERACTION::UTILS::MapExtractor const& EleTypeMapExtractor() const;

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

    //! pointer to input/ouput data container
    Teuchos::RCP<STR::TIMINT::BaseDataIO> gio_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> beaminteractiondatastate_;

    //! particle handler and binning strategy
    Teuchos::RCP<PARTICLE::ParticleHandler> particlehandler_;

    //! binning strategy
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

    //! periodic bounding box
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> periodic_boundingbox_;

    /// map extractor for split of different element types
    Teuchos::RCP<BEAMINTERACTION::UTILS::MapExtractor> eletypeextractor_;

  }; // class Generic

}
}

#endif
