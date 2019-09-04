/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus bc-file reader

\maintainer Martin Kronbichler

\level 1

Here is everything related with reading a bc file
*/
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_READBC_H_
#define PRE_EXODUS_READBC_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"



namespace EXODUS
{
  //! differentiate between underlying mesh entity
  enum mesh_entity
  {
    bceb,
    bcns,
    bcss
  };

  //! differentiate between corresponding condition_type
  enum cond_type
  {
    element,
    dvol,
    dsurf,
    dline,
    dpoint,
    empty,
    invalid
  };

  //! this is what fully defines a baci element
  struct elem_def
  {
    int id;             ///< refering to mesh_entity id of eb,ns,ss
    mesh_entity me;     ///< refering to underlying mesh entity
    std::string sec;    ///< FLUID,STRUCTURE,ALE,etc.
    std::string desc;   ///< like "MAT 1 EAS full"
    std::string ename;  ///< FLUID,SOLIDSH8,etc
  };

  //! this is what fully defines a baci condition
  struct cond_def
  {
    int id;            ///< refering to mesh_entity id of eb,ns,ss
    mesh_entity me;    ///< refering to underlying mesh entity
    std::string sec;   ///< see valid_condition 'sectionname'
    std::string desc;  ///< see valid_condition 'description'
    int e_id;          ///< refers to datfile 'E num -'
    DRT::Condition::GeometryType gtype;
  };

  void ReadBCFile(const std::string& bcfile, std::vector<EXODUS::elem_def>& eledefs,
      std::vector<EXODUS::cond_def>& condefs);

  EXODUS::elem_def ReadEdef(
      const std::string& mesh_entity, const int id, const std::string& actcond);

  EXODUS::cond_def ReadCdef(
      const std::string& mesh_entity, const int id, const std::string& actcond);

  //! Read bc_entity specifications
  std::vector<std::string> ReadBCEntity(const std::string actcond);

  //! Check condition type against valid types
  inline EXODUS::cond_type CheckCondType(const std::string buffer);

  //! Conversion
  inline std::string CondTypeToString(const EXODUS::cond_type);

  //! Print bc_entity
  void PrintBCDef(std::ostream& os, const EXODUS::elem_def& def);
  void PrintBCDef(std::ostream& os, const EXODUS::cond_def& def);

  // ! Check if periodic boundary conditions are defined
  bool PeriodicBoundaryConditionsFound(std::vector<EXODUS::cond_def> condefs);

  // ! Correct nodal coordinates for periodic boundary conditions
  void CorrectNodalCoordinatesForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, std::vector<EXODUS::cond_def> condefs);

  // ! Correct nodal coordinates in the YZ plane for periodic boundary conditions
  void CorrectYZPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs);

  // ! Correct nodal coordinates in the XZ plane for periodic boundary conditions
  void CorrectXZPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs);

  // ! Correct nodal coordinates in the XY plane for periodic boundary conditions
  void CorrectXYPlaneForPeriodicBoundaryConditions(
      EXODUS::Mesh& mesh, const std::vector<EXODUS::cond_def>& condefs);

}  // namespace EXODUS


#endif /*PRE_EXODUS_READBC_H_*/
