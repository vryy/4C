/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Contains ONLY lists of enumerators and is supposed to be included
       in the header files if required


\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_ENUM_LISTS_H_
#define SRC_DRT_STRUCTURE_NEW_STR_ENUM_LISTS_H_

#include <string>

namespace STR
{
  //! Supported types of energy contributions
  enum EnergyType : int
  {
    internal_energy,                 ///< internal, i.e. strain energy
    kinetic_energy,                  ///< kinetic energy
    beam_contact_penalty_potential,  ///< penalty potential for beam-to-? contact
    beam_interaction_potential,      ///< interaction potential for beam-to-? molecular interactions
    beam_to_beam_link_internal_energy,    ///< internal energy of beam to beam links
    beam_to_beam_link_kinetic_energy,     ///< kinetic energy of beam to beam links
    beam_to_sphere_link_internal_energy,  ///< internal energy of beam to sphere links
    beam_to_sphere_link_kinetic_energy    ///< kinetic energy of beam to shpere links
  };

  //! Map energy type to std::string
  inline std::string EnergyType2String(const enum EnergyType type)
  {
    switch (type)
    {
      case internal_energy:
        return "internal_energy";
      case kinetic_energy:
        return "kinetic_energy";
      case beam_contact_penalty_potential:
        return "beam_contact_penalty_potential";
      case beam_interaction_potential:
        return "beam_interaction_potential";
      case beam_to_beam_link_internal_energy:
        return "beam_to_beam_link_internal_energy";
      case beam_to_beam_link_kinetic_energy:
        return "beam_to_beam_link_kinetic_energy";
      case beam_to_sphere_link_internal_energy:
        return "beam_to_sphere_link_internal_energy";
      case beam_to_sphere_link_kinetic_energy:
        return "beam_to_sphere_link_kinetic_energy";
      default:
        return "unknown_type_of_energy";
    }
    exit(EXIT_FAILURE);
  };

  //! Map std::string to energy type
  inline EnergyType String2EnergyType(const std::string type)
  {
    if (type == "internal_energy")
      return internal_energy;
    else if (type == "kinetic_energy")
      return kinetic_energy;
    else if (type == "beam_contact_penalty_potential")
      return beam_contact_penalty_potential;
    else if (type == "beam_interaction_potential")
      return beam_interaction_potential;
    else if (type == "beam_to_beam_link_internal_energy")
      return beam_to_beam_link_internal_energy;
    else if (type == "beam_to_beam_link_kinetic_energy")
      return beam_to_beam_link_kinetic_energy;
    else if (type == "beam_to_sphere_link_internal_energy")
      return beam_to_sphere_link_internal_energy;
    else if (type == "beam_to_sphere_link_kinetic_energy")
      return beam_to_sphere_link_kinetic_energy;
    else
      dserror("Unknown type of energy %s", type.c_str());
    exit(EXIT_FAILURE);
  };


  enum DifferentiationType
  {
    none,
    elch,
    temp
  };

}  // namespace STR

#endif
