/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_CONDITIONS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_CONDITIONS_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_conditions.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

#include <set>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN


// Forward declarations.
namespace Core::Elements
{
  class FaceElement;
}
namespace BEAMINTERACTION
{
  class BeamToSolidParamsBase;
}
namespace Inpar
{
  namespace BeamToSolid
  {
    enum class BeamToSolidMortarShapefunctions;
  }
}  // namespace Inpar


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
    BeamToSolidCondition(const Teuchos::RCP<const Core::Conditions::Condition>& condition_line,
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params);


    /**
     * \brief Check if a combination of beam and solid id is in this condition.
     */
    bool ids_in_condition(const int id_line, const int id_other) const override;

    /**
     * \brief Clear not reusable data (derived).
     */
    void clear() override;

    /**
     * \brief Create the beam to solid pairs needed for this condition (derived).
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> create_contact_pair(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) override;

    /**
     * \brief Return a pointer to the condition of the other geometry (volume or surface).
     * @return
     */
    Teuchos::RCP<const Core::Conditions::Condition> get_other_condition() const
    {
      return condition_other_;
    }

    /**
     * \brief Create the indirect assembly manager for this condition.
     * @param discret (in) discretization.
     */
    Teuchos::RCP<SUBMODELEVALUATOR::BeamContactAssemblyManager> create_indirect_assembly_manager(
        const Teuchos::RCP<const Core::FE::Discretization>& discret) override;

    /**
     * \brief Return a pointer to the geometry evaluation data in this condition.
     */
    Teuchos::RCP<const GEOMETRYPAIR::GeometryEvaluationDataBase> get_geometry_evaluation_data()
        const
    {
      return geometry_evaluation_data_;
    }

   protected:
    /**
     * \brief Check if a solid ID is in this condition.
     */
    virtual bool id_in_other(const int id_other) const = 0;

    /**
     * \brief Return the created beam contact pair for this condition.
     *
     * This function is called by CreateContactPair where the geometry pair in the created contact
     * pair is initialized.
     *
     * @param ele_ptrs (in) Pointer to the two elements contained in the pair.
     * @return Pointer to the created pair.
     */
    virtual Teuchos::RCP<BEAMINTERACTION::BeamContactPair> create_contact_pair_internal(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) = 0;

   protected:
    //! Pointer to the geometry evaluation data for this condition.
    Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase> geometry_evaluation_data_;

    //! Pointer to the solid condition.
    Teuchos::RCP<const Core::Conditions::Condition> condition_other_;

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
    BeamToSolidConditionVolumeMeshtying(
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_line,
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params);


    /**
     * \brief Build the volume ID sets for this condition.
     *
     * The BuildIdSets method from the base class is called to build the beam IDs.
     */
    void build_id_sets(const Teuchos::RCP<const Core::FE::Discretization>& discretization) override;

   protected:
    /**
     * \brief Return the created beam contact pair for this condition. (derived)
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> create_contact_pair_internal(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) override;

    /**
     * \brief Check if a solid ID is in this condition.
     */
    inline bool id_in_other(const int id_other) const override
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
  class BeamToSolidConditionSurface : public BeamToSolidCondition
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    BeamToSolidConditionSurface(
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_line,
        const Teuchos::RCP<const Core::Conditions::Condition>& condition_other,
        const Teuchos::RCP<const BeamToSolidParamsBase>& beam_to_solid_params,
        const bool is_mesh_tying_in);


    /**
     * \brief Build the surface ID sets for this condition. The BuildIdSets method from the base
     * class is called to build the beam IDs.
     */
    void build_id_sets(const Teuchos::RCP<const Core::FE::Discretization>& discretization) override;

    /**
     * \brief Here we get all face elements that are needed for the created pairs. This includes
     * elements which are not part of any pair, but share a node with a surface of a pair.
     *
     * @param discret (in) discretization.
     */
    void setup(const Teuchos::RCP<const Core::FE::Discretization>& discret) override;

    /**
     * \brief Set the displacement state (derived).
     */
    void set_state(const Teuchos::RCP<const Core::FE::Discretization>& discret,
        const Teuchos::RCP<const Solid::MODELEVALUATOR::BeamInteractionDataState>&
            beaminteraction_data_state) override;

   protected:
    /**
     * \brief Return the created beam contact pair for this condition. (derived)
     */
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> create_contact_pair_internal(
        const std::vector<Core::Elements::Element const*>& ele_ptrs) override;

    /**
     * \brief Check if a solid ID is in this condition.
     */
    inline bool id_in_other(const int id_other) const override
    {
      if (surface_ids_.find(id_other) != surface_ids_.end()) return true;
      return false;
    }

    /**
     * \brief If this condition is a mesh tying condition.
     */
    inline bool is_mesh_tying() const { return is_mesh_tying_; }

    /**
     * \brief If this condition is a contact condition.
     */
    inline bool is_contact() const { return !is_mesh_tying_; }

   private:
    //! If the condition is mesh tying or contact.
    bool is_mesh_tying_;

    //! Map containing the global volume element IDs for each face element of the surface in this
    //! condition.
    std::unordered_map<int, Teuchos::RCP<const Core::Elements::FaceElement>> surface_ids_;
  };


  /**
   * \brief Create a beam-to-solid volume pair depending on the solid volume shape.
   * @tparam BtsClass Beam-to-solid class to create
   * @tparam bts_template_arguments Template arguments when creating the class
   * @param shape (in) Shape of the solid volume
   * @return The created beam-to-solid pair
   */
  template <template <typename...> class BtsClass, typename... BtsTemplateArguments>
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidVolumePairShape(
      const Core::FE::CellType shape);

  /**
   * \brief Create a beam-to-solid volume pair depending on the solid volume shape, without nurbs.
   *
   * This is for pairs which are not compatible with nurbs discretizations.
   *
   * @tparam BtsClass Beam-to-solid class to create
   * @tparam bts_template_arguments Template arguments when creating the class
   * @param shape
   * @return The created beam-to-solid pair
   */
  template <template <typename...> class BtsClass, typename... BtsTemplateArguments>
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidVolumePairShapeNoNurbs(
      const Core::FE::CellType shape);

  /**
   * \brief Create a beam-to-solid volume mortar pair depending on the solid volume shape and mortar
   * shape function(s).
   * @tparam BtsClass Beam-to-solid class to create
   * @tparam bts_mortar_template_arguments
   * @tparam bts_mortar_shape Template collection of type of optional given (e.g. rotational) mortar
   * shape functions
   * @param shape (in) Shape of the solid volume
   * @param mortar_shape_function (in) Mortar shape function to be appended to the
   * bts_template_arguments in this call
   * @param other_mortar_shape_function (in) Other mortar shape functions
   * @return The created beam-to-solid pair
   */
  template <template <typename...> class BtsClass, typename... BtsMortarTemplateArguments,
      typename... BtsMortarShape>
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidVolumePairMortar(
      const Core::FE::CellType shape,
      const Inpar::BeamToSolid::BeamToSolidMortarShapefunctions mortar_shape_function,
      BtsMortarShape... other_mortar_shape_function);

  /**
   * \brief Call the pair creation based on the shape once all mortar shape functions have been
   * converted to template arguments.
   * @tparam BtsClass Beam-to-solid class to create
   * @tparam bts_mortar_template_arguments Template arguments when creating the class
   * @param shape (in) Shape of the solid volume
   * @return The created beam-to-solid pair
   */
  template <template <typename...> class BtsClass, typename... BtsMortarTemplateArguments>
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreateBeamToSolidVolumePairMortar(
      const Core::FE::CellType shape);
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
