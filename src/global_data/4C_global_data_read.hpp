/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities to read and fill global data

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_GLOBAL_DATA_READ_HPP
#define FOUR_C_GLOBAL_DATA_READ_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"
#include "4C_io_inputreader.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  /// setup the discretizations
  void ReadFields(
      Global::Problem& problem, Core::IO::DatFileReader& reader, const bool readmesh = true);

  void ReadMicroFields(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// set up supporting processors for micro-scale discretizations
  void ReadMicrofieldsNPsupport(Global::Problem& problem);

  /// read global parameters
  void ReadParameter(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of contact constitutive laws
  void ReadContactConstitutiveLaws(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of materials
  void ReadMaterials(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// setup map between materials of original and cloned elements
  void ReadCloningMaterialMap(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of conditions
  void ReadConditions(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of result tests
  void ReadResult(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void ReadKnots(Global::Problem& problem, Core::IO::DatFileReader& reader);

  /// input of particles
  void ReadParticles(Global::Problem& problem, Core::IO::DatFileReader& reader);
}  // namespace Global

FOUR_C_NAMESPACE_CLOSE

#endif