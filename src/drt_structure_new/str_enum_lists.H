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

namespace STR {

  //! Supported types of energy contributions
  enum EnergyType : int
  {
    internal_energy,                ///< internal, i.e. strain energy
    kinetic_energy,                 ///< kinetic energy
    beam_contact_penalty_potential, ///< penalty potential for beam-to-? contact
    beam_interaction_potential      ///< interaction potential for beam-to-? molecular interactions
  };

  //! Map energy type to std::string
  inline std::string EnergyType2String( const enum EnergyType type )
  {
    switch ( type )
    {
      case internal_energy:
        return "internal_energy";
      case kinetic_energy:
        return "kinetic_energy";
      case beam_contact_penalty_potential:
        return "beam_contact";
      case beam_interaction_potential:
        return "beam_potential";
      default:
        return "unknown_type_of_energy";
    }
    exit(EXIT_FAILURE);
  };

}

#endif
