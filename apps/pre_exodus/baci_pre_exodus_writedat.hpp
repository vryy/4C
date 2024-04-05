/*----------------------------------------------------------------------*/
/*! \file

\brief pre_exodus .dat-file writer

\level 1


Here is everything related with writing a dat-file
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_PRE_EXODUS_WRITEDAT_HPP
#define FOUR_C_PRE_EXODUS_WRITEDAT_HPP

#include "baci_config.hpp"

#include "baci_pre_exodus_readbc.hpp"

BACI_NAMESPACE_OPEN

namespace EXODUS
{
  //! write datfile
  int WriteDatFile(const std::string& datfile, const EXODUS::Mesh& mymesh,
      const std::string& headfile, const std::vector<EXODUS::elem_def>& eledefs,
      const std::vector<EXODUS::cond_def>& condefs,
      const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo);

  //! datfile Intro
  void WriteDatIntro(const std::string& headfile, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! headfile content into datfile
  void WriteDatHead(const std::string& headfile, std::ostream& dat);

  //! remove a specific input file section
  void RemoveDatSection(const std::string& secname, std::string& headstring);

  //! conditions
  void WriteDatConditions(
      const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! special condition SurfLocsys: calculates Normal and Tangent
  std::vector<double> CalcNormalSurfLocsys(const int ns_id, const EXODUS::Mesh& m);

  //! DesignNode - Node Topology
  void WriteDatDesignTopology(
      const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! return set of nodes in bc_entity
  std::set<int> GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m);

  //! Nodes into datfile
  void WriteDatNodes(const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements into datfile
  void WriteDatEles(const std::vector<EXODUS::elem_def>& eledefs, const EXODUS::Mesh& mymesh,
      std::ostream& dat,
      const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo);

  //! Elements from eblock
  void DatEles(Teuchos::RCP<const EXODUS::ElementBlock> eb, const EXODUS::elem_def& acte,
      int& startele, std::ostream& dat,
      const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo,
      const int eb_id);

  inline std::string CondGeomTypeToString(const EXODUS::cond_def& def)
  {
    switch (def.gtype)
    {
      case DRT::Condition::Volume:
        return "DVOL  ";
      case DRT::Condition::Surface:
        return "DSURF ";
      case DRT::Condition::Line:
        return "DLINE ";
      case DRT::Condition::Point:
        return "DPOINT";
      case DRT::Condition::NoGeom:
        return "";
      default:
        dserror("Unknown Condition GeometryType");
    }
    return "";
  }

}  // namespace EXODUS

BACI_NAMESPACE_CLOSE

#endif
