/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_writedat.H

\brief pre_exodus .dat-file writer

\level 1

\maintainer Martin Kronbichler

Here is everything related with writing a dat-file
 */
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_WRITEDAT_H_
#define PRE_EXODUS_WRITEDAT_H_

#include "pre_exodus_readbc.H"


/*!
\brief EXODUS: namespace of the ExodusII-to-BACI preprocessing tool

*/
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
  void RemoveDatSection(const std::string secname, std::string& headstring);

  //! design description
  void WriteDatDesign(const std::vector<EXODUS::cond_def>& condefs, std::ostream& dat);

  //! conditions
  void WriteDatConditions(
      const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! special condition SurfLocsys: calculates Normal and Tangent
  std::vector<double> CalcNormalSurfLocsys(const int ns_id, const EXODUS::Mesh& m);

  //! DesignNode - Node Topology
  void WriteDatDesignTopology(
      const std::vector<EXODUS::cond_def>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! return set of nodes in bc_entity
  const std::set<int> GetNsFromBCEntity(const EXODUS::cond_def& e, const EXODUS::Mesh& m);

  //! Nodes into datfile
  void WriteDatNodes(const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements into datfile
  void WriteDatEles(const std::vector<EXODUS::elem_def>& eledefs, const EXODUS::Mesh& mymesh,
      std::ostream& dat,
      const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo);

  //! Elements from eblock
  void DatEles(Teuchos::RCP<const EXODUS::ElementBlock> eb, const EXODUS::elem_def& acte,
      int& struele, std::ostream& dat,
      const std::map<int, std::map<int, std::vector<std::vector<double>>>>& elecenterlineinfo,
      const int eb_id);

  inline std::string CondGeomTypeToString(const EXODUS::cond_def& def)
  {
    switch (def.gtype)
    {
      case DRT::Condition::Volume:
        return "DVOL  ";
        break;
      case DRT::Condition::Surface:
        return "DSURF ";
        break;
      case DRT::Condition::Line:
        return "DLINE ";
        break;
      case DRT::Condition::Point:
        return "DPOINT";
        break;
      case DRT::Condition::NoGeom:
        return "";
        break;
      default:
        dserror("Unknown Condition GeometryType");
    }
    return "";
  }

}  // namespace EXODUS


#endif /*PRE_EXODUS_WRITEDAT_H_*/
