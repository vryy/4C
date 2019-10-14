/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_CONDITIONS_H_
#define BEAM_TO_SOLID_CONDITIONS_H_


#include "../drt_lib/drt_dserror.H"

#include <map>
#include <vector>
#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace DRT
{
  class Discretization;
  class Condition;
  class Element;
}  // namespace DRT
namespace INPAR
{
  namespace BEAMTOSOLID
  {
    enum class BeamToSolidInteraction;
  }  // namespace BEAMTOSOLID
}  // namespace INPAR
namespace BEAMINTERACTION
{
  class BeamContactPair;
  class BeamContactParams;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  /**
   * \brief This base class represents a single beam-to-solid interaction condition.
   */
  class BeamToSolidCondition
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
     * \brief Create the beam to solid pairs needed for this condition.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @param params_ptr (in) Pointer beam contact parameters.
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidPair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr) = 0;

    /**
     * \brief Build the ID sets for this condition. The ID sets will be used to check if an element
     * is in this condition.
     */
    virtual void BuildIdSets();

    /**
     * \brief Check if a combination of beam and solid id is in this condition.
     */
    bool IdsInCondition(const int id_line, const int id_other) const;

   protected:
    /**
     * \brief Check if a solid ID is in this condition.
     */
    virtual bool IdInOther(const int id_other) const = 0;

   protected:
    //! Pointer to the beam condition.
    Teuchos::RCP<const DRT::Condition> condition_line_;

    //! Pointer to the solid condition.
    Teuchos::RCP<const DRT::Condition> condition_other_;

    //! Set containing the beam element IDs.
    std::set<int> line_ids_;
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
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidPair(
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
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidPair(
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


  /**
   * \brief This class manages all beam-to-solid interaction conditions.
   */
  class BeamToSolidConditions
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidConditions();

    /**
     * \brief Destructor.
     */
    ~BeamToSolidConditions(){};

    /**
     * \brief Get all beam-to-solid interaction conditions from the discretization.
     *
     * This method searches the discretization for input beam-to-solid conditions, finds the correct
     * line-to-solid / volume pairings and adds them to the class variable \ref condition_map_.
     */
    void SetBeamToSolidConditions(const Teuchos::RCP<const DRT::Discretization>& discret);

    /**
     * \brief Build the ID sets on all contained beam-to-solid interaction conditions.
     */
    void BuildIdSets();

    /**
     * \brief Create the correct pair for the given element pointers.
     *
     * We assume, that each beam-to-solid pair can only be in one beam-to-solid interaction. This
     * function checks which interaction condition contains both elements of this pair and creates
     * the correct pair.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @param params_ptr (in) Pointer beam contact parameters.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreatePair(
        const std::vector<DRT::Element const*>& ele_ptrs,
        const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr);

    /**
     * \brief Return a const reference to the condition map.
     */
    inline const std::map<INPAR::BEAMTOSOLID::BeamToSolidInteraction,
        std::vector<Teuchos::RCP<BeamToSolidCondition>>>&
    GetConditionMap() const
    {
      return condition_map_;
    }

    /**
     * \brief Return a mutable reference to the condition map.
     */
    inline std::map<INPAR::BEAMTOSOLID::BeamToSolidInteraction,
        std::vector<Teuchos::RCP<BeamToSolidCondition>>>&
    GetConditionMapMutable()
    {
      return condition_map_;
    }

    /**
     * \brief Get the total number of beam-to-solid interaction conditions.
     */
    inline unsigned int GetTotalNumberOfConditions() const
    {
      unsigned int count = 0;
      for (const auto& map_pair : condition_map_) count += map_pair.second.size();
      return count;
    }

    /**
     * \brief Check if a combination of beam and solid id is in any beam-to-solid condition.
     */
    bool IdsInConditions(const int id_line, const int id_other) const
    {
      for (auto& map_pair : condition_map_)
        for (auto& condition : map_pair.second)
          if (condition->IdsInCondition(id_line, id_other)) return true;
      return false;
    };

   private:
    //! A map containing all types of beam-to-solid conditions. The map keys are the beam-to-solid
    //! interaction type, the values are vectors with conditions (since we can have multiple
    //! conditions of the same interaction type).
    std::map<INPAR::BEAMTOSOLID::BeamToSolidInteraction,
        std::vector<Teuchos::RCP<BeamToSolidCondition>>>
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
