/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities to read and fill global data

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_GLOBAL_DATA_READ_H
#define BACI_GLOBAL_DATA_READ_H

#include "baci_config.H"

#include "baci_global_data.H"
#include "baci_io_inputreader.H"

BACI_NAMESPACE_OPEN

namespace GLOBAL
{
  /// setup the discretizations
  void ReadFields(DRT::Problem& problem, INPUT::DatFileReader& reader, const bool readmesh = true);

  void ReadMicroFields(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// set up supporting processors for micro-scale discretizations
  void ReadMicrofieldsNPsupport(DRT::Problem& problem);

  /// read global parameters
  void ReadParameter(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of contact constitutive laws
  void ReadContactConstitutiveLaws(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of materials
  void ReadMaterials(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// setup map between materials of original and cloned elements
  void ReadCloningMaterialMap(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of conditions
  void ReadConditions(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of result tests
  void ReadResult(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void ReadKnots(DRT::Problem& problem, INPUT::DatFileReader& reader);

  /// input of particles
  void ReadParticles(DRT::Problem& problem, INPUT::DatFileReader& reader);
}  // namespace GLOBAL

BACI_NAMESPACE_CLOSE

#endif