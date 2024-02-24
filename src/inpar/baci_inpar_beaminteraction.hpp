/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for beaminteraction


\level 2

*/
/*-----------------------------------------------------------*/
#ifndef BACI_INPAR_BEAMINTERACTION_HPP
#define BACI_INPAR_BEAMINTERACTION_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace BEAMINTERACTION
  {
    enum RepartitionStrategy
    {
      repstr_adaptive,  ///< only do repartitioning in case physically necessary
      repstr_everydt    ///< do repartitioning every time step
    };

    enum class SearchStrategy
    {
      bruteforce_with_binning,   ///< coupling pair search based on binning
      bounding_volume_hierarchy  ///< coupling pair search based on a bounding volume hierarchy
    };

    /// type of the used submodel for beaminteraction
    enum SubModelType
    {
      submodel_crosslinking,    ///< evaluate the structural model
      submodel_beamcontact,     ///< evaluate the contact model
      submodel_potential,       ///< evaluate the model for potential-based interactions
      submodel_spherebeamlink,  ///< evaluate model for cell filament interactions
      submodel_vague            ///< undefined model type
    };

    /// type of employed solving strategy for contact
    /// (this enum represents the input file parameter STRATEGY)
    enum Strategy
    {
      bstr_none,    ///< no beam contact
      bstr_penalty  ///< penalty method
    };

    /// type of linker
    enum JointType
    {
      beam3r_line2_rigid,  ///< rigid joint
      beam3r_line2_pin,    ///< pin joint
      truss                ///< truss
    };

    /// type of filament
    enum FilamentType
    {
      filtype_arbitrary,  ///< no special type
      filtype_actin,      ///< actin type
      filtype_collagen    ///< collagen type
    };

    /// type of crosslinker
    enum CrosslinkerType
    {
      linkertype_arbitrary,  ///< binds to all filament
      linkertype_actin,      ///< only binds to actin filaments
      linkertype_collagen,   ///< only binds to collagen filaments
      linkertype_integrin    ///< sphere to beam linker
    };

    /**
     * \brief Types of beam interaction conditions
     */
    enum class BeamInteractionConditions
    {
      //! Default value.
      none,
      //! Beam-to-beam contact.
      beam_to_beam_contact,
      //! Beam-to-solid volume mesh tying.
      beam_to_solid_volume_meshtying,
      //! Beam-to-solid surface mesh tying.
      beam_to_solid_surface_meshtying,
      //! beam-to-beam penalty point coupling.
      beam_to_beam_point_coupling,
      //! Beam-to-solid surface contact.
      beam_to_solid_surface_contact
    };

    //! Map type std::string to enum
    inline enum JointType String2JointType(const std::string& name)
    {
      JointType type = beam3r_line2_rigid;
      if (name == "beam3rline2rigid")
        type = beam3r_line2_rigid;
      else if (name == "beam3rline2pin")
        type = beam3r_line2_pin;
      else if (name == "truss")
        type = truss;
      else
        dserror("invalid filament type std::string ");

      return type;
    };

    //! Map type std::string to enum
    inline enum FilamentType String2FilamentType(const std::string& name)
    {
      FilamentType type = filtype_arbitrary;
      if (name == "arbitrary")
        type = filtype_arbitrary;
      else if (name == "actin")
        type = filtype_actin;
      else if (name == "collagen")
        type = filtype_collagen;
      else
        dserror("invalid filament type std::string ");

      return type;
    };

    //! Map type std::string to enum
    inline enum CrosslinkerType String2CrosslinkerType(const std::string& name)
    {
      CrosslinkerType type = linkertype_arbitrary;
      if (name == "arbitrary")
        type = linkertype_arbitrary;
      else if (name == "actin")
        type = linkertype_actin;
      else if (name == "collagen")
        type = linkertype_collagen;
      else if (name == "integrin")
        type = linkertype_integrin;
      else
      {
        dserror(
            "invalid crosslinker type %s. Possible values are arbitrary, actin, collagen and "
            "integrin.",
            name.c_str());
      }

      return type;
    };

    //! Map action type enum to std::string
    static inline std::string CrosslinkerType2String(const enum CrosslinkerType type)
    {
      switch (type)
      {
        case linkertype_arbitrary:
          return "arbitrary";
        case linkertype_actin:
          return "actin";
        case linkertype_collagen:
          return "collagen";
        case linkertype_integrin:
          return "integrin";
        default:
          return "unknown";
      }
      return "";
    };

    /**
     * \brief Get all available beam interaction conditions, excluding the default value.
     */
    void BeamInteractionConditionsGetAll(
        std::vector<INPAR::BEAMINTERACTION::BeamInteractionConditions>& interactions);

    /// set the beam interaction parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam interaction specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace BEAMINTERACTION

}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
