/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus .dat-file writer

\level 1


Here is everything related with writing a dat-file
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_PRE_EXODUS_WRITEDAT_HPP
#define FOUR_C_PRE_EXODUS_WRITEDAT_HPP

#include "4C_config.hpp"

#include "4C_pre_exodus_readbc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  //! write datfile
  int WriteDatFile(const std::string& datfile, const EXODUS::Mesh& mymesh,
      const std::string& headfile, const std::vector<EXODUS::ElemDef>& eledefs,
      const std::vector<EXODUS::CondDef>& condefs);

  //! datfile Intro
  void WriteDatIntro(const std::string& headfile, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! headfile content into datfile
  void WriteDatHead(const std::string& headfile, std::ostream& dat);

  //! remove a specific input file section
  void RemoveDatSection(const std::string& secname, std::string& headstring);

  //! conditions
  void WriteDatConditions(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! special condition SurfLocsys: calculates Normal and Tangent
  std::vector<double> CalcNormalSurfLocsys(const int ns_id, const EXODUS::Mesh& m);

  //! DesignNode - Node Topology
  void WriteDatDesignTopology(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! return set of nodes in bc_entity
  std::set<int> GetNsFromBCEntity(const EXODUS::CondDef& e, const EXODUS::Mesh& m);

  //! Nodes into datfile
  void WriteDatNodes(const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements into datfile
  void WriteDatEles(
      const std::vector<EXODUS::ElemDef>& eledefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements from eblock
  void DatEles(Teuchos::RCP<const EXODUS::ElementBlock> eb, const EXODUS::ElemDef& acte,
      int& startele, std::ostream& dat, const int eb_id);

  inline std::string CondGeomTypeToString(const EXODUS::CondDef& def)
  {
    switch (def.gtype)
    {
      case CORE::Conditions::geometry_type_volume:
        return "DVOL  ";
      case CORE::Conditions::geometry_type_surface:
        return "DSURF ";
      case CORE::Conditions::geometry_type_line:
        return "DLINE ";
      case CORE::Conditions::geometry_type_point:
        return "DPOINT";
      case CORE::Conditions::geometry_type_no_geom:
        return "";
      default:
        FOUR_C_THROW("Unknown Condition GeometryType");
    }
    return "";
  }

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
