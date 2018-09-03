/*-----------------------------------------------------------------------------------------------*/
/*!
\file str_enum_lists.H

\brief Contains ONLY lists of enumerators and is supposed to be included
       in the header files if required

\maintainer Michael Hiermeier

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
        return "beam_contact";
      case beam_interaction_potential:
        return "beam_potential";
      case beam_to_beam_link_internal_energy:
        return "beam_to_beam_internal";
      case beam_to_beam_link_kinetic_energy:
        return "beam_to_beam_kinetic";
      case beam_to_sphere_link_internal_energy:
        return "beam_to_shpere_internal";
      case beam_to_sphere_link_kinetic_energy:
        return "beam_to_sphere_kinetic";
      default:
        return "unknown_type_of_energy";
    }
    exit(EXIT_FAILURE);
  };

}  // namespace STR

#endif
