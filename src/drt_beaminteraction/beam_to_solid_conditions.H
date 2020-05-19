/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_CONDITIONS_H_
#define BEAM_TO_SOLID_CONDITIONS_H_

#include <set>

#include "beaminteraction_conditions.H"

#include "../drt_lib/drt_dserror.H"

#include <unordered_map>
#include <vector>
#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace DRT
{
  class FaceElement;
}
namespace BEAMINTERACTION
{
  class BeamToSolidParamsBase;
}


namespace BEAMINTERACTION
{
  /**
   * \brief This base class represents a single beam-to-solid interaction condition.
   */
  class BeamToSolidCondition : public BeamInteractionConditionBase
  {
   public:
    /**
     * \brief Constructor.
     *
     * @param condition_line (in) The line condition containing the beam elements.
     * @param condition_other (in) The surface / volume condition containing the solid elements
     * interacting with the beam.
     * @param beam_to_solid_params (in) Pointer to the beam-to-solid parameters.
     */
    BeamToSolidCondition(const Teuchos::RCP<const DRT::Condition>& condition_line,
        const Teuchos::RCP<const DRT::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidCondition(){};

    /**
     * \brief Check if a combination of beam and solid id is in this condition.
     */
    bool IdsInCondition(const int id_line, const int id_other) const override;

    /**
     * \brief Clear not reusable data (derived).
     */
    void Clear() override;

    /**
     * \brief Create the beam to solid pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<DRT::Element const*>& ele_ptrs) override;

    /**
     * \brief Return a pointer to the condition of the other geometry (volume or surface).
     * @return
     */
    const Teuchos::RCP<const DRT::Condition> GetOtherCondition() const { return condition_other_; }

    /**
     * \brief Create the indirect assembly manager for this condition.
     * @param discret (in) Discretization.
     */
    Teuchos::RCP<SUBMODELEVALUATOR::BeamContactAssemblyManager> CreateIndirectAssemblyManager(
        const Teuchos::RCP<const DRT::Discretization>& discret) override;

   protected:
    /**
     * \brief Check if a solid ID is in this condition.
     */
    virtual bool IdInOther(const int id_other) const = 0;

    /**
     * \brief Return the created beam contact pair for this condition.
     *
     * This function is called by CreateContactPair where the geometry pair in the created contact
     * pair is initialized.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @return Pointer to the created pair.
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPairInternal(
        const std::vector<DRT::Element const*>& ele_ptrs) = 0;

   protected:
    //! Pointer to the solid condition.
    Teuchos::RCP<const DRT::Condition> condition_other_;

    //! Vector containing all beam contact pairs created by this condition.
    std::vector<Teuchos::RCP<BeamContactPair>> condition_contact_pairs_;

    //! Pointer to the beam-to-solid parameters.
    Teuchos::RCP<const BeamToSolidParamsBase> beam_to_solid_params_;
  };

  /**
   * \brief This base class represents a single beam-to-solid volume mesh tying interaction
   * condition.
   */
  class BeamToSolidConditionVolumeMeshtying : public BeamToSolidCondition
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    BeamToSolidConditionVolumeMeshtying(const Teuchos::RCP<const DRT::Condition>& condition_line,
        const Teuchos::RCP<const DRT::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params);

    /**
     * \brief Destructor.
     */
    ~BeamToSolidConditionVolumeMeshtying(){};

    /**
     * \brief Build the volume ID sets for this condition.
     *
     * The BuildIdSets method from the base class is called to build the beam IDs.
     */
    void BuildIdSets() override;

   protected:
    /**
     * \brief Return the created beam contact pair for this condition. (derived)
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPairInternal(
        const std::vector<DRT::Element const*>& ele_ptrs) override;

    /**
     * \brief Check if a solid ID is in this condition.
     */
    inline bool IdInOther(const int id_other) const override
    {
      if (volume_ids_.find(id_other) != volume_ids_.end()) return true;
      return false;
    }

   private:
    //! Set containing the volume element IDs.
    std::set<int> volume_ids_;
  };

  /**
   * \brief This base class represents a single beam-to-solid surface mesh tying interaction
   * condition.
   */
  class BeamToSolidConditionSurfaceMeshtying : public BeamToSolidCondition
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    BeamToSolidConditionSurfaceMeshtying(const Teuchos::RCP<const DRT::Condition>& condition_line,
        const Teuchos::RCP<const DRT::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params);

    /**
     * \brief Destructor.
     */
    ~BeamToSolidConditionSurfaceMeshtying(){};

    /**
     * \brief Build the surface ID sets for this condition. The BuildIdSets method from the base
     * class is called to build the beam IDs.
     */
    void BuildIdSets() override;

    /**
     * \brief Here we get all face elements that are needed for the created pairs. This includes
     * elements which are not part of any pair, but share a node with a surface of a pair.
     *
     * @param discret (in) Discretization.
     */
    void Setup(const Teuchos::RCP<const DRT::Discretization>& discret) override;

    /**
     * \brief Set the displacement state (derived).
     */
    void SetState(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
            beaminteraction_data_state) override;

   protected:
    /**
     * \brief Return the created beam contact pair for this condition. (derived)
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPairInternal(
        const std::vector<DRT::Element const*>& ele_ptrs) override;

    /**
     * \brief Check if a solid ID is in this condition.
     */
    inline bool IdInOther(const int id_other) const override
    {
      if (surface_ids_.find(id_other) != surface_ids_.end()) return true;
      return false;
    }

   private:
    //! Map containing the global volume element IDs for each face element of the surface in this
    //! condition.
    std::unordered_map<int, Teuchos::RCP<const DRT::FaceElement>> surface_ids_;
  };  // namespace BEAMINTERACTION
}  // namespace BEAMINTERACTION

#endif
