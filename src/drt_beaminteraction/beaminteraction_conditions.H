/*----------------------------------------------------------------------*/
/*! \file

\brief Base classes to manage the beam interactions defined by conditions.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_INTERACTION_CONDITIONS_H_
#define BEAM_INTERACTION_CONDITIONS_H_


#include "../drt_inpar/inpar_beaminteraction.H"

#include "../drt_lib/drt_dserror.H"

#include <map>
#include <set>
#include <vector>
#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace DRT
{
  class Discretization;
  class Condition;
  class Element;
}  // namespace DRT
namespace BEAMINTERACTION
{
  class BeamContactPair;
  class BeamContactParams;
}  // namespace BEAMINTERACTION
namespace GEOMETRYPAIR
{
  class GeometryEvaluationDataBase;
}
namespace STR
{
  namespace MODELEVALUATOR
  {
    class BeamInteractionDataState;
  }
}  // namespace STR


namespace BEAMINTERACTION
{
  /**
   * \brief This abstract base class represents a single beam interaction condition.
   */
  class BeamInteractionConditionBase
  {
   public:
    /**
     * \brief Constructor.
     *
     * @param condition_line (in) The line condition containing the beam elements.
     */
    BeamInteractionConditionBase(const Teuchos::RCP<const DRT::Condition>& condition_line);

    /**
     * \brief Destructor.
     */
    virtual ~BeamInteractionConditionBase(){};

    /**
     * \brief Create the beam contact pair needed for this condition.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @param params_ptr (in) Pointer beam contact parameters.
     * @return Pointer to the created pair.
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr) = 0;

    /**
     * \brief Build the ID sets for this condition. The ID sets will be used to check if an element
     * is in this condition.
     */
    virtual void BuildIdSets();

    /**
     * \brief Set the displacement state.
     *
     * @param discret (in) Discretization.
     * @param beaminteraction_data_state (in) Datastate of the beaminteraction model evaluator.
     */
    virtual void SetState(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
            beaminteraction_data_state)
    {
    }

    /**
     * \brief Setup geometry data.
     * @param discret (in) Discretization.
     */
    virtual void Setup(const Teuchos::RCP<const DRT::Discretization>& discret);

    /**
     * \brief Clear not reusable data.
     */
    virtual void Clear();

    /**
     * \brief Check if a combination of a beam element ID and another element (beam, solid, ...)
     * ID is in this condition.
     */
    virtual bool IdsInCondition(const int id_line, const int id_other) const = 0;

    /**
     * \brief Return a pointer to the geometry evaluation data in this condition.
     */
    const Teuchos::RCP<const GEOMETRYPAIR::GeometryEvaluationDataBase> GetGeometryEvaluationData()
        const
    {
      return geometry_evaluation_data_;
    }

   protected:
    //! Pointer to the beam condition.
    Teuchos::RCP<const DRT::Condition> condition_line_;

    //! Set containing the beam element IDs.
    std::set<int> line_ids_;

    //! Pointer to the geometry evaluation data for this condition.
    Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase> geometry_evaluation_data_;
  };

  /**
   * \brief This class manages all beam interaction conditions.
   */
  class BeamInteractionConditions
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamInteractionConditions();

    /**
     * \brief Destructor.
     */
    virtual ~BeamInteractionConditions(){};

    /**
     * \brief Get all beam interaction conditions from the discretization.
     *
     * This method searches the discretization for input beam interaction conditions, finds the
     * correct line-to- line / surface / volume pairings and adds them to the class variable \ref
     * condition_map_.
     */
    void SetBeamInteractionConditions(const Teuchos::RCP<const DRT::Discretization>& discret);

    /**
     * \brief Build the ID sets on all contained beam interaction conditions.
     */
    void BuildIdSets();

    /**
     * \brief Set the displacement state.
     *
     * @param discret (in) Discretization.
     * @param beaminteraction_data_state (in) Datastate of the beaminteraction model evaluator.
     */
    virtual void SetState(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
            beaminteraction_data_state);

    /**
     * \brief Setup data in the conditions.
     * @param discret (in) Discretization.
     */
    virtual void Setup(const Teuchos::RCP<const DRT::Discretization>& discret);

    /**
     * \brief Clear not reusable data in the conditions.
     */
    virtual void Clear();

    /**
     * \brief Create the correct pair for the given element pointers.
     *
     * We assume, that each beam interaction pair can only be in one beam interaction condition.
     * This function checks which interaction condition contains both elements of this pair and
     * creates the correct pair.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @param params_ptr (in) Pointer beam contact parameters.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateContactPair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr);

    /**
     * \brief Return a const reference to the condition map.
     */
    inline const std::map<INPAR::BEAMINTERACTION::BeamInteractionConditions,
        std::vector<Teuchos::RCP<BeamInteractionConditionBase>>>&
    GetConditionMap() const
    {
      return condition_map_;
    }

    /**
     * \brief Return a mutable reference to the condition map.
     */
    inline std::map<INPAR::BEAMINTERACTION::BeamInteractionConditions,
        std::vector<Teuchos::RCP<BeamInteractionConditionBase>>>&
    GetConditionMapMutable()
    {
      return condition_map_;
    }

    /**
     * \brief Get the total number of beam interaction conditions.
     */
    inline unsigned int GetTotalNumberOfConditions() const
    {
      unsigned int count = 0;
      for (const auto& map_pair : condition_map_) count += map_pair.second.size();
      return count;
    }

    /**
     * \brief Check if a combination of beam element and other element id is in any beam interaction
     * condition.
     */
    bool IdsInConditions(const int id_line, const int id_other) const
    {
      for (auto& map_pair : condition_map_)
        for (auto& condition : map_pair.second)
          if (condition->IdsInCondition(id_line, id_other)) return true;
      return false;
    };

   private:
    //! A map containing all types of beam interaction conditions. The map keys are the beam
    //! interaction type, the values are vectors with conditions (since we can have multiple
    //! conditions of the same interaction type).
    std::map<INPAR::BEAMINTERACTION::BeamInteractionConditions,
        std::vector<Teuchos::RCP<BeamInteractionConditionBase>>>
        condition_map_;
  };

  /**
   * \brief Get the global element IDs of all elements in a condition.
   *
   * @param condition (in) A pointer to the condition.
   * @param element_ids (out) A vector with all global element IDs.
   */
  void ConditionToElementIds(
      const Teuchos::RCP<const DRT::Condition>& condition, std::vector<int>& element_ids);

}  // namespace BEAMINTERACTION

#endif
