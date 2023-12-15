/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Contains ONLY lists of enumerators and is supposed to be included
       in the header files if required


\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_ENUM_LISTS_H
#define BACI_STRUCTURE_NEW_ENUM_LISTS_H

#include "baci_utils_exceptions.H"

#include <string>

BACI_NAMESPACE_OPEN

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


  //! for coupled, monolithic problems: linearization w.r.t. other primary variable
  enum class DifferentiationType : int
  {
    none,
    elch,
    temp
  };

  enum class MatBlockType
  {
    displ_displ,  ///< Kdd block (structural block)
    displ_lm,     ///< Kdz block (of the corresponding model evaluator)
    lm_displ,     ///< Kzd block (of the corresponding model evaluator)
    lm_lm,        ///< Kzz block (of the corresponding model evaluator)
  };

  inline std::string MatBlockType2String(const enum MatBlockType type)
  {
    switch (type)
    {
      case MatBlockType::displ_displ:
        return "block_displ_displ";
      case MatBlockType::displ_lm:
        return "block_displ_lm";
      case MatBlockType::lm_displ:
        return "block_lm_displ";
      case MatBlockType::lm_lm:
        return "block_lm_lm";
      default:
        return "unknown matrix block type";
    }
  }

}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif
