/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_CONDITIONS_H_
#define BEAM_TO_SOLID_CONDITIONS_H_


#include "beaminteraction_conditions.H"

#include "../drt_lib/drt_dserror.H"

#include <map>
#include <vector>
#include <Teuchos_RCP.hpp>


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
     */
    BeamToSolidCondition(const Teuchos::RCP<const DRT::Condition>& condition_line,
        const Teuchos::RCP<const DRT::Condition>& condition_other);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidCondition(){};

    /**
     * \brief Check if a combination of beam and solid id is in this condition.
     */
    bool IdsInCondition(const int id_line, const int id_other) const override;

   protected:
    /**
     * \brief Check if a solid ID is in this condition.
     */
    virtual bool IdInOther(const int id_other) const = 0;

   protected:
    //! Pointer to the solid condition.
    Teuchos::RCP<const DRT::Condition> condition_other_;
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
        const Teuchos::RCP<const DRT::Condition>& condition_other);

    /**
     * \brief Destructor.
     */
    ~BeamToSolidConditionVolumeMeshtying(){};

    /**
     * \brief Create the beam to solid pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr) override;

    /**
     * \brief Build the volume ID sets for this condition. The BuildIdSets method from the base
     * class is called to build the beam IDs.
     */
    void BuildIdSets() override;

   protected:
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
        const Teuchos::RCP<const DRT::Condition>& condition_other);

    /**
     * \brief Destructor.
     */
    ~BeamToSolidConditionSurfaceMeshtying(){};

    /**
     * \brief Create the beam to solid pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr) override;

    /**
     * \brief Build the surface ID sets for this condition. The BuildIdSets method from the base
     * class is called to build the beam IDs.
     */
    void BuildIdSets() override { dserror("Not yet implemented!"); }

   protected:
    /**
     * \brief Check if a solid ID is in this condition.
     */
    inline bool IdInOther(const int id_other) const override
    {
      dserror("Not yet implemented!");
      return false;
    }
  };
}  // namespace BEAMINTERACTION

#endif
