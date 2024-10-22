// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GLOBAL_DATA_READ_HPP
#define FOUR_C_GLOBAL_DATA_READ_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  /// setup the discretizations
  void read_fields(
      Global::Problem& problem, Core::IO::DatFileReader& reader, const bool read_mesh = true);

  void read_micro_fields(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// set up supporting processors for micro-scale discretizations
  void read_microfields_np_support(Global::Problem& problem);

  /// read global parameters
  void read_parameter(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of contact constitutive laws
  void read_contact_constitutive_laws(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of materials
  void read_materials(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// setup map between materials of original and cloned elements
  void read_cloning_material_map(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of conditions
  void read_conditions(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of result tests
  void read_result(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void read_knots(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of particles
  void read_particles(Global::Problem& problem, Core::IO::DatFileReader& reader);
}  // namespace Global

FOUR_C_NAMESPACE_CLOSE

#endif