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
#include "4C_io_input_file.hpp"
#include "4C_io_meshreader.hpp"

#include <filesystem>

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  /**
   * This legacy function prepares an InputFile object with all globally known sections. This object
   * can be used to read a file, convert file formats or emit the known sections as metadata.
   */
  [[nodiscard]] FOUR_C_API(FOUR_C_CORE) Core::IO::InputFile set_up_input_file(MPI_Comm comm);

  /**
   * Write the metadata that is not tied to any specific input of a physical problem and thus
   * not modified by developers of physics modules. This includes version information or general
   * definitions of cell geometries.
   */
  FOUR_C_API(FOUR_C_CORE) void emit_general_metadata(Core::IO::YamlNodeRef node);

  /// Read the discretization. This mainly means the mesh. The MeshReader is returned to the
  /// caller in case it needs to be used for further processing.
  FOUR_C_API(FOUR_C_CORE)
  std::unique_ptr<Core::IO::MeshReader> read_discretization(
      Global::Problem& problem, Core::IO::InputFile& input, const bool read_mesh = true);

  FOUR_C_API(FOUR_C_CORE)
  void read_micro_fields(Global::Problem& problem, const std::filesystem::path& input_path);

  /// set up supporting processors for micro-scale discretizations
  FOUR_C_API(FOUR_C_CORE) void read_microfields_np_support(Global::Problem& problem);

  /// read global parameters
  FOUR_C_API(FOUR_C_CORE) void read_parameter(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of contact constitutive laws
  FOUR_C_API(FOUR_C_CORE)
  void read_contact_constitutive_laws(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of materials
  FOUR_C_API(FOUR_C_CORE) void read_materials(Global::Problem& problem, Core::IO::InputFile& input);

  /// setup map between materials of original and cloned elements
  FOUR_C_API(FOUR_C_CORE)
  void read_cloning_material_map(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of conditions
  FOUR_C_API(FOUR_C_CORE)
  void read_conditions(Global::Problem& problem, Core::IO::InputFile& input,
      const Core::IO::MeshReader& mesh_reader);

  /// input of result tests
  FOUR_C_API(FOUR_C_CORE) void read_result(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of knots for isogeometric analysis
  FOUR_C_API(FOUR_C_CORE) void read_knots(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of particles
  FOUR_C_API(FOUR_C_CORE) void read_particles(Global::Problem& problem, Core::IO::InputFile& input);

  /// input of spatial fields
  FOUR_C_API(FOUR_C_CORE)
  void read_fields(Global::Problem& problem, Core::IO::InputFile& input,
      const Core::IO::MeshReader& mesh_reader);
}  // namespace Global

FOUR_C_NAMESPACE_CLOSE

#endif