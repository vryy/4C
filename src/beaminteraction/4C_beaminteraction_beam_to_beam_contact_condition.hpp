/*----------------------------------------------------------------------*/
/*! \file

\brief Object to hold and manage a beam-to-beam condition.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_CONDITION_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_CONDITION_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_conditions.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

#include <set>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace BEAMINTERACTION
{
  /**
   * \brief This base class represents a single beam-to-beam contact condition.
   */
  class BeamToBeamContactCondition : public BeamInteractionConditionBase
  {
   public:
    /**
     * \brief Constructor.
     *
     * @param condition_line_1 (in) The first line condition containing the beam elements.
     * @param condition_line_2 (in) The other line condition containing the beam elements.
     */
    BeamToBeamContactCondition(
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_line_1,
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_line_2);

    /**
     * \brief Build the ID sets for this condition.
     *
     * The BuildIdSets method from the base class is called to build the beam IDs.
     */
    void BuildIdSets(const Teuchos::RCP<const Discret::Discretization>& discretization) override;

    /**
     * \brief Check if a combination of beam and beam id is in this condition.
     */
    bool IdsInCondition(const int id_line, const int id_other) const override;

    /**
     * \brief Clear not reusable data (derived).
     */
    void Clear() override;

    /**
     * \brief Create the beam to beam pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) override;

   protected:
    /**
     * \brief Check if a ID is in a condition.
     */
    inline static bool id_is_in_condition(const std::set<int>& id_set, const int id)
    {
      return id_set.find(id) != id_set.end();
    }

   private:
    //! Pointer to the other line condition.
    Teuchos::RCP<const Core::Conditions::Condition> condition_other_;

    //! Vector containing all beam contact pairs created by this condition.
    std::vector<Teuchos::RCP<BeamContactPair>> condition_contact_pairs_;

    //! Set containing the other line element IDs.
    std::set<int> other_line_ids_;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
