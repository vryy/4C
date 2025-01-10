// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_coupling_adapter.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_geometric_search_matchingoctree.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"

#include <Epetra_Export.h>

#include <algorithm>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Coupling::Adapter::Coupling::Coupling()
    : masterdofmap_(nullptr),
      permmasterdofmap_(nullptr),
      slavedofmap_(nullptr),
      permslavedofmap_(nullptr),
      masterexport_(nullptr),
      slaveexport_(nullptr),
      matmm_(nullptr),
      matsm_(nullptr),
      matmm_trans_(nullptr),
      matsm_trans_(nullptr)
{
  // empty
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_condition_coupling(
    const Core::FE::Discretization& masterdis, std::shared_ptr<const Epetra_Map> mastercondmap,
    const Core::FE::Discretization& slavedis, std::shared_ptr<const Epetra_Map> slavecondmap,
    const std::string& condname, const std::vector<int>& masterdofs,
    const std::vector<int>& slavedofs, bool matchall, const int nds_master, const int nds_slave)
{
  const int numdof = masterdofs.size();
  const int numdof_slave = slavedofs.size();
  if (numdof != numdof_slave)
    FOUR_C_THROW("Received %d master DOFs, but %d slave DOFs", numdof, numdof_slave);

  std::vector<int> masternodes;
  Core::Conditions::find_conditioned_nodes(masterdis, condname, masternodes);
  std::vector<int> slavenodes;
  Core::Conditions::find_conditioned_nodes(slavedis, condname, slavenodes);

  int localmastercount = static_cast<int>(masternodes.size());
  int mastercount;
  int localslavecount = static_cast<int>(slavenodes.size());
  int slavecount;

  Core::Communication::sum_all(&localmastercount, &mastercount, 1, masterdis.get_comm());
  Core::Communication::sum_all(&localslavecount, &slavecount, 1, slavedis.get_comm());

  if (mastercount != slavecount)
    FOUR_C_THROW("got %d master nodes but %d slave nodes for coupling", mastercount, slavecount);

  setup_coupling(masterdis, slavedis, masternodes, slavenodes, masterdofs, slavedofs, matchall,
      1.0e-3, nds_master, nds_slave);

  // test for completeness
  if (static_cast<int>(masternodes.size()) * numdof != masterdofmap_->NumMyElements())
    FOUR_C_THROW("failed to setup master nodes properly");
  if (static_cast<int>(slavenodes.size()) * numdof != slavedofmap_->NumMyElements())
    FOUR_C_THROW("failed to setup slave nodes properly");

  // Now swap in the maps we already had.
  // So we did a little more work than required. But there are cases
  // where we have to do that work (fluid-ale coupling) and we want to
  // use just one setup implementation.
  //
  // The point is to make sure there is only one map for each
  // interface.

  if (not masterdofmap_->PointSameAs(*mastercondmap)) FOUR_C_THROW("master dof map mismatch");

  if (not slavedofmap_->PointSameAs(*slavecondmap))
  {
    FOUR_C_THROW("slave dof map mismatch");
  }

  masterdofmap_ = mastercondmap;
  masterexport_ = std::make_shared<Epetra_Export>(*permmasterdofmap_, *masterdofmap_);

  slavedofmap_ = slavecondmap;
  slaveexport_ = std::make_shared<Epetra_Export>(*permslavedofmap_, *slavedofmap_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_condition_coupling(
    const Core::FE::Discretization& masterdis, std::shared_ptr<const Epetra_Map> mastercondmap,
    const Core::FE::Discretization& slavedis, std::shared_ptr<const Epetra_Map> slavecondmap,
    const std::string& condname, const int numdof, bool matchall, const int nds_master,
    const int nds_slave)
{
  setup_condition_coupling(masterdis, mastercondmap, slavedis, slavecondmap, condname,
      build_dof_vector_from_num_dof(numdof), build_dof_vector_from_num_dof(numdof), matchall,
      nds_master, nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const std::vector<int>& masternodes,
    const std::vector<int>& slavenodes, const std::vector<int>& masterdofs,
    const std::vector<int>& slavedofs, const bool matchall, const double tolerance,
    const int nds_master, const int nds_slave)
{
  std::vector<int> patchedmasternodes(masternodes);
  std::vector<int> permslavenodes;
  match_nodes(
      masterdis, slavedis, patchedmasternodes, permslavenodes, slavenodes, matchall, tolerance);

  // Epetra maps in original distribution

  std::shared_ptr<Epetra_Map> masternodemap =
      std::make_shared<Epetra_Map>(-1, patchedmasternodes.size(), patchedmasternodes.data(), 0,
          Core::Communication::as_epetra_comm(masterdis.get_comm()));

  std::shared_ptr<Epetra_Map> slavenodemap = std::make_shared<Epetra_Map>(-1, slavenodes.size(),
      slavenodes.data(), 0, Core::Communication::as_epetra_comm(slavedis.get_comm()));

  std::shared_ptr<Epetra_Map> permslavenodemap =
      std::make_shared<Epetra_Map>(-1, permslavenodes.size(), permslavenodes.data(), 0,
          Core::Communication::as_epetra_comm(slavedis.get_comm()));

  finish_coupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap, masterdofs,
      slavedofs, nds_master, nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const std::vector<int>& masternodes,
    const std::vector<int>& slavenodes, const int numdof, const bool matchall,
    const double tolerance, const int nds_master, const int nds_slave)
{
  setup_coupling(masterdis, slavedis, masternodes, slavenodes,
      build_dof_vector_from_num_dof(numdof), build_dof_vector_from_num_dof(numdof), matchall,
      tolerance, nds_master, nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_constrained_condition_coupling(
    const Core::FE::Discretization& masterdis, std::shared_ptr<const Epetra_Map> mastercondmap,
    const Core::FE::Discretization& slavedis, std::shared_ptr<const Epetra_Map> slavecondmap,
    const std::string& condname1, const std::string& condname2, const int numdof, bool matchall)
{
  std::vector<int> masternodes1;
  Core::Conditions::find_conditioned_nodes(masterdis, condname1, masternodes1);
  std::vector<int> slavenodes1;
  Core::Conditions::find_conditioned_nodes(slavedis, condname1, slavenodes1);

  std::set<int> masternodes2;
  Core::Conditions::find_conditioned_nodes(masterdis, condname2, masternodes2);
  std::set<int> slavenodes2;
  Core::Conditions::find_conditioned_nodes(slavedis, condname2, slavenodes2);

  // now find all those elements of slavenodes1 and masternodes1 that
  // do not belong to slavenodes2 and masternodes2 at the same time

  std::vector<int> masternodes;
  std::vector<int> slavenodes;

  for (int& masternode1 : masternodes1)
  {
    if (masternodes2.find(masternode1) == masternodes2.end()) masternodes.push_back(masternode1);
  }

  for (int& slavenode1 : slavenodes1)
  {
    if (slavenodes2.find(slavenode1) == slavenodes2.end()) slavenodes.push_back(slavenode1);
  }

  int localmastercount = static_cast<int>(masternodes.size());
  int mastercount;
  int localslavecount = static_cast<int>(slavenodes.size());
  int slavecount;

  Core::Communication::sum_all(&localmastercount, &mastercount, 1, masterdis.get_comm());
  Core::Communication::sum_all(&localslavecount, &slavecount, 1, slavedis.get_comm());

  if (mastercount != slavecount and matchall)
    FOUR_C_THROW("got %d master nodes but %d slave nodes for coupling", mastercount, slavecount);

  setup_coupling(masterdis, slavedis, masternodes, slavenodes, numdof, matchall);

  // test for completeness
  if (static_cast<int>(masternodes.size()) * numdof != masterdofmap_->NumMyElements())
    FOUR_C_THROW("failed to setup master nodes properly");
  if (static_cast<int>(slavenodes.size()) * numdof != slavedofmap_->NumMyElements())
    FOUR_C_THROW("failed to setup slave nodes properly");

  // Now swap in the maps we already had.
  // So we did a little more work than required. But there are cases
  // where we have to do that work (fluid-ale coupling) and we want to
  // use just one setup implementation.
  //
  // The point is to make sure there is only one map for each
  // interface.

  if (not masterdofmap_->PointSameAs(*mastercondmap)) FOUR_C_THROW("master dof map mismatch");

  if (not slavedofmap_->PointSameAs(*slavecondmap)) FOUR_C_THROW("slave dof map mismatch");

  masterdofmap_ = mastercondmap;
  masterexport_ = std::make_shared<Epetra_Export>(*permmasterdofmap_, *masterdofmap_);

  slavedofmap_ = slavecondmap;
  slaveexport_ = std::make_shared<Epetra_Export>(*permslavedofmap_, *slavedofmap_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(std::shared_ptr<const Epetra_Map> slavedofmap,
    std::shared_ptr<const Epetra_Map> permslavedofmap,
    std::shared_ptr<const Epetra_Map> masterdofmap,
    std::shared_ptr<const Epetra_Map> permmasterdofmap)
{
  masterdofmap_ = masterdofmap;
  slavedofmap_ = slavedofmap;
  permmasterdofmap_ = permmasterdofmap;
  permslavedofmap_ = permslavedofmap;

  masterexport_ = std::make_shared<Epetra_Export>(*permmasterdofmap_, *masterdofmap_);
  slaveexport_ = std::make_shared<Epetra_Export>(*permslavedofmap_, *slavedofmap_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const Epetra_Map& masternodes,
    const Epetra_Map& slavenodes, const int numdof, const bool matchall, const double tolerance,
    const int nds_master, const int nds_slave)
{
  if (masternodes.NumGlobalElements() != slavenodes.NumGlobalElements() and matchall)
    FOUR_C_THROW("got %d master nodes but %d slave nodes for coupling",
        masternodes.NumGlobalElements(), slavenodes.NumGlobalElements());

  std::vector<int> mastervect(
      masternodes.MyGlobalElements(), masternodes.MyGlobalElements() + masternodes.NumMyElements());
  std::vector<int> slavevect(
      slavenodes.MyGlobalElements(), slavenodes.MyGlobalElements() + slavenodes.NumMyElements());
  std::vector<int> permslavenodes;

  match_nodes(masterdis, slavedis, mastervect, permslavenodes, slavevect, matchall, tolerance);

  // Epetra maps in original distribution

  std::shared_ptr<Epetra_Map> masternodemap = std::make_shared<Epetra_Map>(-1, mastervect.size(),
      mastervect.data(), 0, Core::Communication::as_epetra_comm(masterdis.get_comm()));

  std::shared_ptr<Epetra_Map> slavenodemap = std::make_shared<Epetra_Map>(slavenodes);

  std::shared_ptr<Epetra_Map> permslavenodemap =
      std::make_shared<Epetra_Map>(-1, permslavenodes.size(), permslavenodes.data(), 0,
          Core::Communication::as_epetra_comm(slavedis.get_comm()));

  finish_coupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap,
      build_dof_vector_from_num_dof(numdof), build_dof_vector_from_num_dof(numdof), nds_master,
      nds_slave);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const Epetra_Map& masternodemap,
    const Epetra_Map& slavenodemap, const Epetra_Map& permslavenodemap, const int numdof)
{
  if (masternodemap.NumGlobalElements() != slavenodemap.NumGlobalElements())
    FOUR_C_THROW("got %d master nodes but %d slave nodes for coupling",
        masternodemap.NumGlobalElements(), slavenodemap.NumGlobalElements());

  // just copy Epetra maps

  std::shared_ptr<Epetra_Map> mymasternodemap = std::make_shared<Epetra_Map>(masternodemap);

  std::shared_ptr<Epetra_Map> myslavenodemap = std::make_shared<Epetra_Map>(slavenodemap);

  std::shared_ptr<Epetra_Map> mypermslavenodemap = std::make_shared<Epetra_Map>(permslavenodemap);

  // build slave to master permutation and dof all maps
  finish_coupling(masterdis, slavedis, mymasternodemap, myslavenodemap, mypermslavenodemap,
      build_dof_vector_from_num_dof(numdof), build_dof_vector_from_num_dof(numdof));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(
    const Core::FE::Discretization& masterdis, const Core::FE::Discretization& slavedis)
{
  // safety check
  if (masterdis.dof_row_map()->NumGlobalElements() != slavedis.dof_row_map()->NumGlobalElements())
    FOUR_C_THROW("got %d master nodes but %d slave nodes for coupling",
        masterdis.dof_row_map()->NumGlobalElements(), slavedis.dof_row_map()->NumGlobalElements());

  // get master dof maps and build exporter
  permmasterdofmap_ = std::make_shared<Epetra_Map>(*slavedis.dof_row_map());
  masterdofmap_ = std::make_shared<Epetra_Map>(*masterdis.dof_row_map());
  masterexport_ = std::make_shared<Epetra_Export>(*permmasterdofmap_, *masterdofmap_);

  // get slave dof maps and build exporter
  permslavedofmap_ = std::make_shared<Epetra_Map>(*masterdis.dof_row_map());
  slavedofmap_ = std::make_shared<Epetra_Map>(*slavedis.dof_row_map());
  slaveexport_ = std::make_shared<Epetra_Export>(*permslavedofmap_, *slavedofmap_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const std::vector<std::vector<int>>& masternodes_vec,
    const std::vector<std::vector<int>>& slavenodes_vec, const int numdof, const bool matchall,
    const double tolerance, const int nds_master, const int nds_slave)
{
  // vectors with master and slave node maps (from input) for every coupling condition
  // Permuted slave node map for each coupling conditions from match_nodes()
  std::vector<std::shared_ptr<const Epetra_Map>> masternodemap_cond;
  std::vector<std::shared_ptr<const Epetra_Map>> slavenodemap_cond;
  std::vector<std::shared_ptr<const Epetra_Map>> permslavenodemap_cond;

  for (unsigned i = 0; i < masternodes_vec.size(); ++i)
  {
    std::vector<int> masternodes = masternodes_vec.at(i);
    std::vector<int> slavenodes = slavenodes_vec.at(i);

    std::vector<int> permslavenodes;

    match_nodes(masterdis, slavedis, masternodes, permslavenodes, slavenodes, matchall, tolerance);

    masternodemap_cond.push_back(std::make_shared<Epetra_Map>(-1, masternodes.size(),
        masternodes.data(), 0, Core::Communication::as_epetra_comm(masterdis.get_comm())));
    slavenodemap_cond.push_back(std::make_shared<Epetra_Map>(-1, slavenodes.size(),
        slavenodes.data(), 0, Core::Communication::as_epetra_comm(slavedis.get_comm())));
    permslavenodemap_cond.push_back(std::make_shared<Epetra_Map>(-1, permslavenodes.size(),
        permslavenodes.data(), 0, Core::Communication::as_epetra_comm(slavedis.get_comm())));
  }

  // merge maps for all conditions, but keep order (= keep assignment of permuted slave node map and
  // master map)
  auto masternodemap = Core::LinAlg::MultiMapExtractor::merge_maps_keep_order(masternodemap_cond);
  auto slavenodemap = Core::LinAlg::MultiMapExtractor::merge_maps_keep_order(slavenodemap_cond);
  auto permslavenodemap =
      Core::LinAlg::MultiMapExtractor::merge_maps_keep_order(permslavenodemap_cond);

  finish_coupling(masterdis, slavedis, masternodemap, slavenodemap, permslavenodemap,
      build_dof_vector_from_num_dof(numdof), build_dof_vector_from_num_dof(numdof), nds_master,
      nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::match_nodes(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, std::vector<int>& masternodes,
    std::vector<int>& permslavenodes, const std::vector<int>& slavenodes, const bool matchall,
    const double tolerance)
{
  // match master and slave nodes using octree
  auto tree = Core::GeometricSearch::NodeMatchingOctree();
  tree.init(masterdis, masternodes, 150, tolerance);
  tree.setup();

  std::map<int, std::pair<int, double>> coupling;
  tree.find_match(slavedis, slavenodes, coupling);

  if (masternodes.size() != coupling.size() and matchall)
    FOUR_C_THROW(
        "Did not get 1:1 correspondence. \nmasternodes.size()=%d (%s), coupling.size()=%d (%s)",
        masternodes.size(), masterdis.name().c_str(), coupling.size(), slavedis.name().c_str());

  // extract permutation

  std::vector<int> patchedmasternodes;
  patchedmasternodes.reserve(coupling.size());
  permslavenodes.reserve(slavenodes.size());

  for (int gid : masternodes)
  {
    // We allow to hand in master nodes that do not take part in the
    // coupling. If this is undesired behaviour the user has to make
    // sure all nodes were used.
    if (coupling.find(gid) != coupling.end())
    {
      std::pair<int, double>& coupled = coupling[gid];
      patchedmasternodes.push_back(gid);
      permslavenodes.push_back(coupled.first);
    }
  }

  // return new list of master nodes via reference
  swap(masternodes, patchedmasternodes);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::finish_coupling(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, std::shared_ptr<Epetra_Map> masternodemap,
    std::shared_ptr<Epetra_Map> slavenodemap, std::shared_ptr<Epetra_Map> permslavenodemap,
    const std::vector<int>& masterdofs, const std::vector<int>& slavedofs, const int nds_master,
    const int nds_slave)
{
  // we expect to get maps of exactly the same shape
  if (not masternodemap->PointSameAs(*permslavenodemap))
    FOUR_C_THROW("master and permuted slave node maps do not match");

  // export master nodes to slave node distribution

  // To do so we create vectors that contain the values of the master
  // maps, assigned to the slave maps. On the master side we actually
  // create just a view on the map! This vector must not be changed!
  std::shared_ptr<Core::LinAlg::Vector<int>> masternodevec =
      std::make_shared<Core::LinAlg::Vector<int>>(
          *permslavenodemap, masternodemap->MyGlobalElements());

  std::shared_ptr<Core::LinAlg::Vector<int>> permmasternodevec =
      std::make_shared<Core::LinAlg::Vector<int>>(*slavenodemap);

  Epetra_Export masternodeexport(*permslavenodemap, *slavenodemap);
  const int err = permmasternodevec->Export(*masternodevec, masternodeexport, Insert);
  if (err) FOUR_C_THROW("failed to export master nodes");

  std::shared_ptr<const Epetra_Map> permmasternodemap =
      std::make_shared<Epetra_Map>(-1, permmasternodevec->MyLength(), permmasternodevec->Values(),
          0, Core::Communication::as_epetra_comm(masterdis.get_comm()));

  if (not slavenodemap->PointSameAs(*permmasternodemap))
    FOUR_C_THROW("slave and permuted master node maps do not match");

  masternodevec = nullptr;
  permmasternodevec = nullptr;

  build_dof_maps(masterdis, slavedis, masternodemap, slavenodemap, permmasternodemap,
      permslavenodemap, masterdofs, slavedofs, nds_master, nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::build_dof_maps(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis,
    const std::shared_ptr<const Epetra_Map>& masternodemap,
    const std::shared_ptr<const Epetra_Map>& slavenodemap,
    const std::shared_ptr<const Epetra_Map>& permmasternodemap,
    const std::shared_ptr<const Epetra_Map>& permslavenodemap, const std::vector<int>& masterdofs,
    const std::vector<int>& slavedofs, const int nds_master, const int nds_slave)
{
  build_dof_maps(masterdis, *masternodemap, *permmasternodemap, masterdofmap_, permmasterdofmap_,
      masterexport_, masterdofs, nds_master);
  build_dof_maps(slavedis, *slavenodemap, *permslavenodemap, slavedofmap_, permslavedofmap_,
      slaveexport_, slavedofs, nds_slave);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> Coupling::Adapter::Coupling::build_dof_vector_from_num_dof(const int numdof)
{
  std::vector<int> dofvec;
  if (numdof > 0)
  {
    dofvec.resize(numdof);
    std::iota(dofvec.begin(), dofvec.end(), 0);
  }
  else
    dofvec = {-1};

  return dofvec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::build_dof_maps(const Core::FE::Discretization& dis,
    const Epetra_Map& nodemap, const Epetra_Map& permnodemap,
    std::shared_ptr<const Epetra_Map>& dofmap, std::shared_ptr<const Epetra_Map>& permdofmap,
    std::shared_ptr<Epetra_Export>& exporter, const std::vector<int>& coupled_dofs,
    const int nds) const
{
  // communicate dofs

  std::vector<int> dofmapvec;
  std::map<int, std::vector<int>> dofs;

  const int* nodes = nodemap.MyGlobalElements();
  const int numnode = nodemap.NumMyElements();

  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* actnode = dis.g_node(nodes[i]);

    // ----------------------------------------------------------------
    // get all periodic boundary conditions on this node
    // slave nodes do not contribute dofs, we skip them
    // ----------------------------------------------------------------
    std::vector<Core::Conditions::Condition*> thiscond;
    actnode->get_condition("SurfacePeriodic", thiscond);

    if (thiscond.empty())
    {
      actnode->get_condition("LinePeriodic", thiscond);
    }

    if (!thiscond.empty())
    {
      // loop them and check, whether this is a pbc pure master node
      // for all previous conditions
      unsigned ntimesmaster = 0;
      for (auto& cond : thiscond)
      {
        const auto& mymasterslavetoggle = cond->parameters().get<std::string>("MASTER_OR_SLAVE");

        if (mymasterslavetoggle == "Master")
        {
          ++ntimesmaster;
        }
      }

      if (ntimesmaster < thiscond.size())
      {
        // this node is not a master and does not own its own dofs
        continue;
      }
    }

    const std::vector<int> dof = dis.dof(nds, actnode);
    const int numdof = coupled_dofs.size();
    if (numdof > static_cast<int>(dof.size()))
      FOUR_C_THROW(
          "got just %d dofs at node %d (lid=%d) but expected %d", dof.size(), nodes[i], i, numdof);
    for (int idof = 0; idof < numdof; idof++)
    {
      copy(dof.data() + coupled_dofs[idof], dof.data() + coupled_dofs[idof] + 1,
          back_inserter(dofs[nodes[i]]));
      copy(dof.data() + coupled_dofs[idof], dof.data() + coupled_dofs[idof] + 1,
          back_inserter(dofmapvec));
    }
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(), dofmapvec.end());
  if (pos != dofmapvec.end() and *pos < 0) FOUR_C_THROW("illegal dof number %d", *pos);

  // dof map is the original, unpermuted distribution of dofs
  dofmap = std::make_shared<Epetra_Map>(-1, dofmapvec.size(), dofmapvec.data(), 0,
      Core::Communication::as_epetra_comm(dis.get_comm()));

  dofmapvec.clear();

  Core::Communication::Exporter exportdofs(nodemap, permnodemap, dis.get_comm());
  exportdofs.do_export(dofs);

  const int* permnodes = permnodemap.MyGlobalElements();
  const int permnumnode = permnodemap.NumMyElements();

  for (int i = 0; i < permnumnode; ++i)
  {
    const std::vector<int>& dof = dofs[permnodes[i]];
    copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  dofs.clear();

  // permuted dof map according to a given permuted node map
  permdofmap = std::make_shared<Epetra_Map>(-1, dofmapvec.size(), dofmapvec.data(), 0,
      Core::Communication::as_epetra_comm(dis.get_comm()));

  // prepare communication plan to create a dofmap out of a permuted
  // dof map
  exporter = std::make_shared<Epetra_Export>(*permdofmap, *dofmap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Coupling::Adapter::Coupling::master_to_slave(
    const Core::LinAlg::Vector<double>& mv) const
{
  std::shared_ptr<Core::LinAlg::Vector<double>> sv =
      std::make_shared<Core::LinAlg::Vector<double>>(*slavedofmap_);

  master_to_slave(mv, *sv);

  return sv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Coupling::Adapter::Coupling::slave_to_master(
    const Core::LinAlg::Vector<double>& sv) const
{
  std::shared_ptr<Core::LinAlg::Vector<double>> mv =
      std::make_shared<Core::LinAlg::Vector<double>>(*masterdofmap_);

  slave_to_master(sv, *mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_FEVector> Coupling::Adapter::Coupling::master_to_slave(
    const Epetra_FEVector& mv) const
{
  std::shared_ptr<Epetra_FEVector> sv =
      std::make_shared<Epetra_FEVector>(*slavedofmap_, mv.NumVectors());

  Core::LinAlg::VectorView sv_view(*sv);
  Core::LinAlg::VectorView mv_view(mv);
  master_to_slave(mv_view, sv_view);

  return sv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_FEVector> Coupling::Adapter::Coupling::slave_to_master(
    const Epetra_FEVector& sv) const
{
  std::shared_ptr<Epetra_FEVector> mv =
      std::make_shared<Epetra_FEVector>(*masterdofmap_, sv.NumVectors());

  Core::LinAlg::VectorView sv_view(sv);
  Core::LinAlg::VectorView mv_view(*mv);
  slave_to_master(sv_view, mv_view);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Coupling::Adapter::Coupling::master_to_slave(
    const Core::LinAlg::MultiVector<double>& mv) const
{
  std::shared_ptr<Core::LinAlg::MultiVector<double>> sv =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*slavedofmap_, mv.NumVectors());

  master_to_slave(mv, *sv);

  return sv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MultiVector<double>> Coupling::Adapter::Coupling::slave_to_master(
    const Core::LinAlg::MultiVector<double>& sv) const
{
  std::shared_ptr<Core::LinAlg::MultiVector<double>> mv =
      std::make_shared<Core::LinAlg::MultiVector<double>>(*masterdofmap_, sv.NumVectors());

  slave_to_master(sv, *mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::master_to_slave(
    const Core::LinAlg::MultiVector<double>& mv, Core::LinAlg::MultiVector<double>& sv) const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not mv.Map().PointSameAs(*masterdofmap_)) FOUR_C_THROW("master dof map vector expected");
  if (not sv.Map().PointSameAs(*slavedofmap_)) FOUR_C_THROW("slave dof map vector expected");
  if (sv.NumVectors() != mv.NumVectors())
    FOUR_C_THROW("column number mismatch %d!=%d", sv.NumVectors(), mv.NumVectors());
#endif

  Core::LinAlg::MultiVector<double> perm(*permslavedofmap_, mv.NumVectors());
  std::copy(mv.Values(), mv.Values() + (mv.MyLength() * mv.NumVectors()), perm.Values());

  const int err = sv.Export(perm, *slaveexport_, Insert);
  if (err) FOUR_C_THROW("Export to slave distribution returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::master_to_slave(
    const Core::LinAlg::Vector<int>& mv, Core::LinAlg::Vector<int>& sv) const
{
  Core::LinAlg::Vector<int> perm(*permslavedofmap_);
  std::copy(mv.Values(), mv.Values() + (mv.MyLength()), perm.Values());

  const int err = sv.Export(perm, *slaveexport_, Insert);
  if (err) FOUR_C_THROW("Export to slave distribution returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::slave_to_master(
    const Core::LinAlg::MultiVector<double>& sv, Core::LinAlg::MultiVector<double>& mv) const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not mv.Map().PointSameAs(*masterdofmap_)) FOUR_C_THROW("master dof map vector expected");
  if (not sv.Map().PointSameAs(*slavedofmap_))
  {
    std::cout << "slavedofmap_" << std::endl;
    std::cout << *slavedofmap_ << std::endl;
    std::cout << "sv" << std::endl;
    std::cout << sv.Map() << std::endl;
    FOUR_C_THROW("slave dof map vector expected");
  }
  if (sv.NumVectors() != mv.NumVectors())
    FOUR_C_THROW("column number mismatch %d!=%d", sv.NumVectors(), mv.NumVectors());
#endif

  Core::LinAlg::MultiVector<double> perm(*permmasterdofmap_, sv.NumVectors());
  std::copy(sv.Values(), sv.Values() + (sv.MyLength() * sv.NumVectors()), perm.Values());

  const int err = mv.Export(perm, *masterexport_, Insert);
  if (err) FOUR_C_THROW("Export to master distribution returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::slave_to_master(
    const Core::LinAlg::Vector<int>& sv, Core::LinAlg::Vector<int>& mv) const
{
  Core::LinAlg::Vector<int> perm(*permmasterdofmap_);
  std::copy(sv.Values(), sv.Values() + (sv.MyLength()), perm.Values());

  const int err = mv.Export(perm, *masterexport_, Insert);
  if (err) FOUR_C_THROW("Export to master distribution returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::fill_master_to_slave_map(std::map<int, int>& rowmap) const
{
  for (int i = 0; i < masterdofmap_->NumMyElements(); ++i)
  {
    rowmap[masterdofmap_->GID(i)] = permslavedofmap_->GID(i);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::fill_slave_to_master_map(std::map<int, int>& rowmap) const
{
  for (int i = 0; i < slavedofmap_->NumMyElements(); ++i)
  {
    rowmap[slavedofmap_->GID(i)] = permmasterdofmap_->GID(i);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Coupling::Adapter::Coupling::slave_to_master_map(Epetra_Map& slave)
{
  int nummyele = 0;
  std::vector<int> globalelements;
  const std::shared_ptr<Epetra_Map> slavemap = Core::LinAlg::allreduce_e_map(slave);
  for (int i = 0; i < slavemap->NumMyElements(); ++i)
  {
    int lid = permslavedofmap_->LID(slavemap->GID(i));
    if (lid != -1)
    {
      globalelements.push_back(masterdofmap_->GID(lid));
      nummyele++;
    }
  }

  return std::make_shared<Epetra_Map>(-1, nummyele, globalelements.data(), 0, slave.Comm());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Map> Coupling::Adapter::Coupling::master_to_slave_map(Epetra_Map& master)
{
  int nummyele = 0;
  std::vector<int> globalelements;
  const std::shared_ptr<Epetra_Map> mastermap = Core::LinAlg::allreduce_e_map(master);
  for (int i = 0; i < mastermap->NumMyElements(); ++i)
  {
    int lid = permmasterdofmap_->LID(mastermap->GID(i));
    if (lid != -1)
    {
      globalelements.push_back(slavedofmap_->GID(lid));
      nummyele++;
    }
  }

  return std::make_shared<Epetra_Map>(-1, nummyele, globalelements.data(), 0, master.Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::SparseMatrix> Coupling::Adapter::Coupling::master_to_perm_master(
    const Core::LinAlg::SparseMatrix& sm) const
{
  std::shared_ptr<Epetra_CrsMatrix> permsm =
      std::make_shared<Epetra_CrsMatrix>(Copy, *permmasterdofmap_, sm.max_num_entries());

  // OK. You cannot use the same exporter for different matrices. So we
  // recreate one all the time... This has to be optimized later on.
  Epetra_Export exporter(*permmasterdofmap_, *masterdofmap_);
  int err = permsm->Import(*sm.epetra_matrix(), exporter, Insert);

  if (err) FOUR_C_THROW("Import failed with err=%d", err);

  permsm->FillComplete(sm.domain_map(), *permmasterdofmap_);

  // create a SparseMatrix that wraps the new CrsMatrix.
  return std::make_shared<Core::LinAlg::SparseMatrix>(
      permsm, Core::LinAlg::View, sm.explicit_dirichlet(), sm.save_graph());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::SparseMatrix> Coupling::Adapter::Coupling::slave_to_perm_slave(
    const Core::LinAlg::SparseMatrix& sm) const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not sm.row_map().PointSameAs(*slavedofmap_)) FOUR_C_THROW("slave dof map vector expected");
  if (not sm.filled()) FOUR_C_THROW("matrix must be filled");
#endif

  std::shared_ptr<Epetra_CrsMatrix> permsm =
      std::make_shared<Epetra_CrsMatrix>(Copy, *permslavedofmap_, sm.max_num_entries());

  // OK. You cannot use the same exporter for different matrices. So we
  // recreate one all the time... This has to be optimized later on.

  Epetra_Export exporter(*permslavedofmap_, *slavedofmap_);
  int err = permsm->Import(*sm.epetra_matrix(), exporter, Insert);

  if (err) FOUR_C_THROW("Import failed with err=%d", err);

  permsm->FillComplete(sm.domain_map(), *permslavedofmap_);

  // create a SparseMatrix that wraps the new CrsMatrix.
  return std::make_shared<Core::LinAlg::SparseMatrix>(
      permsm, Core::LinAlg::View, sm.explicit_dirichlet(), sm.save_graph());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::Coupling::setup_coupling_matrices(const Epetra_Map& shiftedmastermap,
    const Epetra_Map& masterdomainmap, const Epetra_Map& slavedomainmap)
{
  // we always use the masterdofmap for the domain
  matmm_ = std::make_shared<Epetra_CrsMatrix>(Copy, shiftedmastermap, 1, true);
  matsm_ = std::make_shared<Epetra_CrsMatrix>(Copy, shiftedmastermap, 1, true);

  matmm_trans_ = std::make_shared<Epetra_CrsMatrix>(Copy, masterdomainmap, 1, true);
  matsm_trans_ = std::make_shared<Epetra_CrsMatrix>(Copy, *perm_slave_dof_map(), 1, true);

  int length = shiftedmastermap.NumMyElements();
  double one = 1.;
  for (int i = 0; i < length; ++i)
  {
    int sgid = perm_slave_dof_map()->GID(i);
    int mgid = master_dof_map()->GID(i);
    int shiftedmgid = shiftedmastermap.GID(i);

    int err = matmm_->InsertGlobalValues(shiftedmgid, 1, &one, &mgid);
    if (err != 0)
      FOUR_C_THROW(
          "InsertGlobalValues for entry (%d,%d) failed with err=%d", shiftedmgid, mgid, err);

    err = matsm_->InsertGlobalValues(shiftedmgid, 1, &one, &sgid);
    if (err != 0)
      FOUR_C_THROW(
          "InsertGlobalValues for entry (%d,%d) failed with err=%d", shiftedmgid, sgid, err);

    err = matmm_trans_->InsertGlobalValues(mgid, 1, &one, &shiftedmgid);
    if (err != 0)
      FOUR_C_THROW(
          "InsertGlobalValues for entry (%d,%d) failed with err=%d", mgid, shiftedmgid, err);

    err = matsm_trans_->InsertGlobalValues(sgid, 1, &one, &shiftedmgid);
    if (err != 0)
      FOUR_C_THROW(
          "InsertGlobalValues for entry (%d,%d) failed with err=%d", sgid, shiftedmgid, err);
  }

  matmm_->FillComplete(masterdomainmap, shiftedmastermap);
  matsm_->FillComplete(slavedomainmap, shiftedmastermap);

  matmm_trans_->FillComplete(shiftedmastermap, masterdomainmap);
  matsm_trans_->FillComplete(shiftedmastermap, *perm_slave_dof_map());

  // communicate slave to master matrix

  std::shared_ptr<Epetra_CrsMatrix> tmp =
      std::make_shared<Epetra_CrsMatrix>(Copy, slavedomainmap, 1);

  Epetra_Import exporter(slavedomainmap, *perm_slave_dof_map());
  int err = tmp->Import(*matsm_trans_, exporter, Insert);
  if (err) FOUR_C_THROW("Import failed with err=%d", err);

  tmp->FillComplete(shiftedmastermap, slavedomainmap);
  matsm_trans_ = tmp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map>& Coupling::Adapter::Coupling::ma_dof_map_ptr()
{
  return masterdofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& Coupling::Adapter::Coupling::ma_dof_map() const
{
  if (!masterdofmap_) FOUR_C_THROW("The masterdofmap_ has not been initialized correctly!");
  return *masterdofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map>& Coupling::Adapter::Coupling::permuted_ma_dof_map_ptr()
{
  return permmasterdofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& Coupling::Adapter::Coupling::permuted_ma_dof_map() const
{
  if (!permmasterdofmap_) FOUR_C_THROW("The permmasterdofmap_ has not been initialized correctly!");
  return *permmasterdofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map>& Coupling::Adapter::Coupling::sl_dof_map_ptr()
{
  return slavedofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& Coupling::Adapter::Coupling::sl_dof_map() const
{
  if (!slavedofmap_) FOUR_C_THROW("The slavedofmap_ has not been initialized correctly!");
  return *slavedofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map>& Coupling::Adapter::Coupling::permuted_sl_dof_map_ptr()
{
  return permslavedofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Map& Coupling::Adapter::Coupling::permuted_sl_dof_map() const
{
  if (!permslavedofmap_) FOUR_C_THROW("The permslavedofmap_ has not been initialized correctly!");
  return *permslavedofmap_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Export>& Coupling::Adapter::Coupling::ma_exporter_ptr()
{
  return masterexport_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Export& Coupling::Adapter::Coupling::ma_exporter() const
{
  if (!masterexport_) FOUR_C_THROW("The masterexport_ has not been initialized correctly!");
  return *masterexport_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Epetra_Export>& Coupling::Adapter::Coupling::sl_exporter_ptr()
{
  return slaveexport_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Epetra_Export& Coupling::Adapter::Coupling::sl_exporter() const
{
  if (!slaveexport_) FOUR_C_THROW("The slaveexport_ has not been initialized correctly!");
  return *slaveexport_;
}

FOUR_C_NAMESPACE_CLOSE
