/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities to read and fill global data

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_GLOBAL_DATA_READ_HPP
#define FOUR_C_GLOBAL_DATA_READ_HPP

#include "baci_config.hpp"

#include "baci_global_data.hpp"
#include "baci_io_inputreader.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GLOBAL
{
  /// setup the discretizations
  void ReadFields(
      GLOBAL::Problem& problem, INPUT::DatFileReader& reader, const bool readmesh = true);

  void ReadMicroFields(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// set up supporting processors for micro-scale discretizations
  void ReadMicrofieldsNPsupport(GLOBAL::Problem& problem);

  /// read global parameters
  void ReadParameter(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of contact constitutive laws
  void ReadContactConstitutiveLaws(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of materials
  void ReadMaterials(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// setup map between materials of original and cloned elements
  void ReadCloningMaterialMap(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of conditions
  void ReadConditions(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of result tests
  void ReadResult(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void ReadKnots(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);

  /// input of particles
  void ReadParticles(GLOBAL::Problem& problem, INPUT::DatFileReader& reader);
}  // namespace GLOBAL

FOUR_C_NAMESPACE_CLOSE

#endif