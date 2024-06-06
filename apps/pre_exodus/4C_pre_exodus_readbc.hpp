/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus bc-file reader


\level 1

Here is everything related with reading a bc file
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_PRE_EXODUS_READBC_HPP
#define FOUR_C_PRE_EXODUS_READBC_HPP

#include "4C_config.hpp"

#include "4C_pre_exodus_reader.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  //! differentiate between underlying mesh entity
  enum MeshEntity
  {
    bceb,
    bcns,
    bcss
  };

  //! differentiate between corresponding condition_type
  enum CondType
  {
    element,
    dvol,
    dsurf,
    dline,
    dpoint,
    empty,
    invalid
  };

  //! this is what fully defines a 4C element
  struct ElemDef
  {
    int id;             ///< refering to mesh_entity id of eb,ns,ss
    MeshEntity me;      ///< refering to underlying mesh entity
    std::string sec;    ///< FLUID,STRUCTURE,ALE,etc.
    std::string desc;   ///< like "MAT 1 EAS full"
    std::string ename;  ///< FLUID,SOLIDSH8,etc
  };

  //! this is what fully defines a 4C condition
  struct CondDef
  {
    int id;            ///< refering to mesh_entity id of eb,ns,ss
    MeshEntity me;     ///< refering to underlying mesh entity
    std::string sec;   ///< see valid_condition 'sectionname'
    std::string desc;  ///< see valid_condition 'description'
    int e_id;          ///< refers to datfile 'E num -'
    Core::Conditions::GeometryType gtype;
  };

  void ReadBCFile(const std::string& bcfile, std::vector<EXODUS::ElemDef>& eledefs,
      std::vector<EXODUS::CondDef>& condefs);

  EXODUS::ElemDef ReadEdef(
      const std::string& mesh_entity, const int id, const std::string& actcond);

  EXODUS::CondDef ReadCdef(
      const std::string& mesh_entity, const int id, const std::string& actcond);

  //! Read bc_entity specifications
  std::vector<std::string> ReadBCEntity(const std::string actcond);

  //! Check condition type against valid types
  inline EXODUS::CondType CheckCondType(const std::string buffer);

  //! Conversion
  inline std::string CondTypeToString(const EXODUS::CondType);

  //! Print bc_entity
  void PrintBCDef(std::ostream& os, const EXODUS::ElemDef& def);
  void PrintBCDef(std::ostream& os, const EXODUS::CondDef& def);

  // ! Check if periodic boundary conditions are defined
  bool PeriodicBoundaryConditionsFound(std::vector<EXODUS::CondDef> condefs);

  // ! Correct nodal coordinates for periodic boundary conditions
  void CorrectNodalCoordinatesForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, std::vector<EXODUS::CondDef> condefs);

  // ! Correct nodal coordinates in the YZ plane for periodic boundary conditions
  void CorrectYZPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::CondDef>& condefs);

  // ! Correct nodal coordinates in the XZ plane for periodic boundary conditions
  void CorrectXZPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::CondDef>& condefs);

  // ! Correct nodal coordinates in the XY plane for periodic boundary conditions
  void CorrectXYPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::CondDef>& condefs);

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
