#ifndef FOUR_C_PRE_EXODUS_WRITEDAT_HPP
#define FOUR_C_PRE_EXODUS_WRITEDAT_HPP

#include "4C_config.hpp"

#include "4C_pre_exodus_readbc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  //! write datfile
  int write_dat_file(const std::string& datfile, const EXODUS::Mesh& mymesh,
      const std::string& headfile, const std::vector<EXODUS::ElemDef>& eledefs,
      const std::vector<EXODUS::CondDef>& condefs);

  //! datfile Intro
  void write_dat_intro(const std::string& headfile, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! headfile content into datfile
  void write_dat_head(const std::string& headfile, std::ostream& dat);

  //! remove a specific input file section
  void remove_dat_section(const std::string& secname, std::string& headstring);

  //! conditions
  void write_dat_conditions(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! special condition SurfLocsys: calculates Normal and Tangent
  std::vector<double> calc_normal_surf_locsys(const int ns_id, const EXODUS::Mesh& m);

  //! DesignNode - Node Topology
  void write_dat_design_topology(
      const std::vector<EXODUS::CondDef>& condefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! return set of nodes in bc_entity
  std::set<int> get_ns_from_bc_entity(const EXODUS::CondDef& e, const EXODUS::Mesh& m);

  //! Nodes into datfile
  void write_dat_nodes(const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements into datfile
  void write_dat_eles(
      const std::vector<EXODUS::ElemDef>& eledefs, const EXODUS::Mesh& mymesh, std::ostream& dat);

  //! Elements from eblock
  void dat_eles(const EXODUS::ElementBlock& eb, const EXODUS::ElemDef& acte, int& startele,
      std::ostream& dat, const int eb_id);

  inline std::string cond_geom_type_to_string(const EXODUS::CondDef& def)
  {
    switch (def.gtype)
    {
      case Core::Conditions::geometry_type_volume:
        return "DVOL  ";
      case Core::Conditions::geometry_type_surface:
        return "DSURF ";
      case Core::Conditions::geometry_type_line:
        return "DLINE ";
      case Core::Conditions::geometry_type_point:
        return "DPOINT";
      case Core::Conditions::geometry_type_no_geom:
        return "";
      default:
        FOUR_C_THROW("Unknown Condition GeometryType");
    }
    return "";
  }

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
