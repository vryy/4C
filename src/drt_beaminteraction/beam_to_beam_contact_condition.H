/*----------------------------------------------------------------------*/
/*! \file

\brief Object to hold and manage a beam-to-beam condition.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_BEAM_CONTACT_CONDITION_H_
#define BEAM_TO_BEAM_CONTACT_CONDITION_H_

#include <set>

#include "beaminteraction_conditions.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_dserror.H"

#include <unordered_map>
#include <vector>
#include <Teuchos_RCP.hpp>


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
     * @param condition_line (in) The first line condition containing the beam elements.
     * @param condition_other (in) The other line condition containing the beam elements.
     */
    BeamToBeamContactCondition(const Teuchos::RCP<const DRT::Condition>& condition_line_1,
        const Teuchos::RCP<const DRT::Condition>& condition_line_2);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToBeamContactCondition(){};

    /**
     * \brief Build the ID sets for this condition.
     *
     * The BuildIdSets method from the base class is called to build the beam IDs.
     */
    void BuildIdSets() override;

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
        const std::vector<DRT::Element const*>& ele_ptrs) override;

   protected:
    /**
     * \brief Check if a ID is in a condition.
     */
    inline static bool IdIsInCondition(const std::set<int>& id_set, const int id)
    {
      if (id_set.find(id) != id_set.end()) return true;
      return false;
    }

   protected:
    //! Pointer to the other line condition.
    Teuchos::RCP<const DRT::Condition> condition_other_;

    //! Vector containing all beam contact pairs created by this condition.
    std::vector<Teuchos::RCP<BeamContactPair>> condition_contact_pairs_;

    //! Set containing the other line element IDs.
    std::set<int> other_line_ids_;
  };
}  // namespace BEAMINTERACTION

#endif
