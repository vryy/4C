/*----------------------------------------------------------------------------*/
/** \file

\brief Special coupling routines to handle the possible changing number of DoFs
       from node to node during XFEM simulations


\level 3

*/
/*----------------------------------------------------------------------------*/


#include "4C_xfem_xfield_field_coupling.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_coupling_adapter.hpp"

#include <Epetra_Export.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::XFieldField::Coupling::Coupling()
    : Core::Adapter::Coupling(), isinit_(false), min_dof_dis_(min_dof_unknown)
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::init(const enum MinDofDiscretization& min_dof_dis)
{
  min_dof_dis_ = min_dof_dis;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::XFieldField::Coupling::master_to_slave(
    const Teuchos::RCP<const Epetra_Vector>& mv, const enum XFEM::MapType& map_type) const
{
  Teuchos::RCP<Epetra_Vector> sv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return Core::Adapter::Coupling::master_to_slave(mv);
      break;
    case XFEM::map_nodes:
      sv = Teuchos::rcp(new Epetra_Vector(*slavenodemap_));
      break;
  }

  master_to_slave(mv, map_type, sv);
  return sv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::XFieldField::Coupling::slave_to_master(
    const Teuchos::RCP<const Epetra_Vector>& sv, const enum XFEM::MapType& map_type) const
{
  Teuchos::RCP<Epetra_Vector> mv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return Core::Adapter::Coupling::slave_to_master(sv);
      break;
    case XFEM::map_nodes:
      mv = Teuchos::rcp(new Epetra_Vector(*masternodemap_));
      break;
  }

  slave_to_master(sv, map_type, mv);
  return mv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::XFieldField::Coupling::master_to_slave(
    const Teuchos::RCP<const Epetra_MultiVector>& mv, const enum XFEM::MapType& map_type) const
{
  Teuchos::RCP<Epetra_MultiVector> sv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return Core::Adapter::Coupling::master_to_slave(mv);
      break;
    case XFEM::map_nodes:
      sv = Teuchos::rcp(new Epetra_MultiVector(*slavenodemap_, mv->NumVectors()));
      break;
  }

  master_to_slave(mv, map_type, sv);
  return sv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::XFieldField::Coupling::slave_to_master(
    const Teuchos::RCP<const Epetra_MultiVector>& sv, const enum XFEM::MapType& map_type) const
{
  Teuchos::RCP<Epetra_MultiVector> mv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return Core::Adapter::Coupling::slave_to_master(sv);
      break;
    case XFEM::map_nodes:
      mv = Teuchos::rcp(new Epetra_MultiVector(*masternodemap_, sv->NumVectors()));
      break;
  }

  slave_to_master(sv, map_type, mv);
  return mv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::master_to_slave(const Teuchos::RCP<const Epetra_MultiVector>& mv,
    const enum XFEM::MapType& map_type, Teuchos::RCP<Epetra_MultiVector> sv) const
{
  switch (map_type)
  {
    case XFEM::map_dofs:
    {
      return Core::Adapter::Coupling::master_to_slave(mv, sv);
      break;
    }
    case XFEM::map_nodes:
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not mv->Map().PointSameAs(*masternodemap_))
        FOUR_C_THROW("master node map vector expected");
      if (not sv->Map().PointSameAs(*slavenodemap_)) FOUR_C_THROW("slave node map vector expected");
      if (sv->NumVectors() != mv->NumVectors())
        FOUR_C_THROW("column number mismatch %d!=%d", sv->NumVectors(), mv->NumVectors());
#endif

      Epetra_MultiVector perm(*permslavenodemap_, mv->NumVectors());
      std::copy(mv->Values(), mv->Values() + (mv->MyLength() * mv->NumVectors()), perm.Values());

      const int err = sv->Export(perm, *nodal_slaveexport_, Insert);
      if (err) FOUR_C_THROW("Export to nodal slave distribution returned err=%d", err);
    }  // end: case XFEM::MultiFieldMapExtractor::map_nodes
  }    // end: switch (map_type)
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::slave_to_master(const Teuchos::RCP<const Epetra_MultiVector>& sv,
    const enum XFEM::MapType& map_type, Teuchos::RCP<Epetra_MultiVector> mv) const
{
  switch (map_type)
  {
    case XFEM::map_dofs:
    {
      return Core::Adapter::Coupling::slave_to_master(sv, mv);
      break;
    }
    case XFEM::map_nodes:
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not mv->Map().PointSameAs(*masternodemap_))
        FOUR_C_THROW("master node map vector expected");
      if (not sv->Map().PointSameAs(*slavenodemap_)) FOUR_C_THROW("slave node map vector expected");
      if (sv->NumVectors() != mv->NumVectors())
        FOUR_C_THROW("column number mismatch %d!=%d", sv->NumVectors(), mv->NumVectors());
#endif

      Epetra_MultiVector perm(*permmasternodemap_, sv->NumVectors());
      std::copy(sv->Values(), sv->Values() + (sv->MyLength() * sv->NumVectors()), perm.Values());

      const int err = mv->Export(perm, *nodal_masterexport_, Insert);
      if (err) FOUR_C_THROW("Export to nodal master distribution returned err=%d", err);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::build_dof_maps(const Core::FE::Discretization& masterdis,
    const Core::FE::Discretization& slavedis, const Teuchos::RCP<const Epetra_Map>& masternodemap,
    const Teuchos::RCP<const Epetra_Map>& slavenodemap,
    const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
    const Teuchos::RCP<const Epetra_Map>& permslavenodemap, const std::vector<int>& masterdofs,
    const std::vector<int>& slavedofs, const int nds_master, const int nds_slave)
{
  save_node_maps(masternodemap, slavenodemap, permmasternodemap, permslavenodemap);

  // call base class implementation
  if (masterdofs[0] != -1)
  {
    Core::Adapter::Coupling::build_dof_maps(masterdis, slavedis, masternodemap, slavenodemap,
        permmasternodemap, permslavenodemap, masterdofs, slavedofs, nds_master, nds_slave);
    return;
  }

  check_init();
  /* This map contains the information how many DoF's per node have to be
   * considered, i.e. the number of DoF's per node of the min-dof discretization.
   * The map-key is the corresponding max-dof discretization nodal coupling GID. */
  std::map<int, unsigned> my_mindofpernode;
  switch (min_dof_dis())
  {
    case min_dof_slave:
    {
      build_min_dof_maps(slavedis, *slavenodemap, *permslavenodemap, sl_dof_map_ptr(),
          permuted_sl_dof_map_ptr(), sl_exporter_ptr(), *masternodemap, my_mindofpernode);
      build_max_dof_maps(masterdis, *masternodemap, *permmasternodemap, ma_dof_map_ptr(),
          permuted_ma_dof_map_ptr(), ma_exporter_ptr(), my_mindofpernode);
      break;
    }
    case min_dof_master:
    {
      build_min_dof_maps(masterdis, *masternodemap, *permmasternodemap, ma_dof_map_ptr(),
          permuted_ma_dof_map_ptr(), ma_exporter_ptr(), *slavenodemap, my_mindofpernode);
      build_max_dof_maps(slavedis, *slavenodemap, *permslavenodemap, sl_dof_map_ptr(),
          permuted_sl_dof_map_ptr(), sl_exporter_ptr(), my_mindofpernode);
      break;
    }
    case min_dof_unknown:
    {
      FOUR_C_THROW(
          "The discretization with the minimum number of DoF's \n"
          "per node is unknown or cannot be identified, since it \n"
          "changes from node to node. This case needs extra \n"
          "communication effort and is currently unsupported.");
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::save_node_maps(
    const Teuchos::RCP<const Epetra_Map>& masternodemap,
    const Teuchos::RCP<const Epetra_Map>& slavenodemap,
    const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
    const Teuchos::RCP<const Epetra_Map>& permslavenodemap)
{
  masternodemap_ = masternodemap;
  slavenodemap_ = slavenodemap;
  permmasternodemap_ = permmasternodemap;
  permslavenodemap_ = permslavenodemap;

  nodal_masterexport_ =
      Teuchos::rcp<Epetra_Export>(new Epetra_Export(*permmasternodemap, *masternodemap));
  nodal_slaveexport_ =
      Teuchos::rcp<Epetra_Export>(new Epetra_Export(*permslavenodemap, *slavenodemap));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::build_min_dof_maps(const Core::FE::Discretization& min_dis,
    const Epetra_Map& min_nodemap, const Epetra_Map& min_permnodemap,
    Teuchos::RCP<const Epetra_Map>& min_dofmap, Teuchos::RCP<const Epetra_Map>& min_permdofmap,
    Teuchos::RCP<Epetra_Export>& min_exporter, const Epetra_Map& max_nodemap,
    std::map<int, unsigned>& my_mindofpernode) const
{
  std::vector<int> dofmapvec;
  std::map<int, std::vector<int>> dofs;

  const int* ngids = min_nodemap.MyGlobalElements();
  const int numnode = min_nodemap.NumMyElements();

  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* actnode = min_dis.g_node(ngids[i]);

    const int numdof = min_dis.num_dof(actnode);
    const std::vector<int> dof = min_dis.dof(0, actnode);
    std::copy(dof.data(), dof.data() + numdof, back_inserter(dofs[ngids[i]]));
    std::copy(dof.data(), dof.data() + numdof, back_inserter(dofmapvec));
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(), dofmapvec.end());
  if (pos != dofmapvec.end() and *pos < 0) FOUR_C_THROW("Illegal DoF number %d", *pos);

  // dof map is the original, unpermuted distribution of dofs
  min_dofmap =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), dofmapvec.data(), 0, min_dis.get_comm()));

  dofmapvec.clear();

  Core::Communication::Exporter exportdofs(min_nodemap, min_permnodemap, min_dis.get_comm());
  exportdofs.do_export(dofs);

  const int* permngids = min_permnodemap.MyGlobalElements();
  const int permnumnode = min_permnodemap.NumMyElements();

  for (int i = 0; i < permnumnode; ++i)
  {
    const std::vector<int>& dof = dofs[permngids[i]];
    std::copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  /* -------------------------------------------------------------------------
   * Get the number of dofs per node which have to be considered at the
   * coupling GID in the max-dof-discretization
   * -------------------------------------------------------------------------*/
  std::map<int, std::vector<int>>::const_iterator cit;
  for (cit = dofs.begin(); cit != dofs.end(); ++cit)
  {
    const int min_permlid = min_permnodemap.LID(cit->first);
    const int max_gid = max_nodemap.GID(min_permlid);
    my_mindofpernode[max_gid] = cit->second.size();
  }
  dofs.clear();

  // permuted dof map according to a given permuted node map
  min_permdofmap =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), dofmapvec.data(), 0, min_dis.get_comm()));

  /* prepare communication plan to create a dofmap out of a permuted
   * dof map */
  min_exporter = Teuchos::rcp<Epetra_Export>(new Epetra_Export(*min_permdofmap, *min_dofmap));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::build_max_dof_maps(const Core::FE::Discretization& max_dis,
    const Epetra_Map& max_nodemap, const Epetra_Map& max_permnodemap,
    Teuchos::RCP<const Epetra_Map>& max_dofmap, Teuchos::RCP<const Epetra_Map>& max_permdofmap,
    Teuchos::RCP<Epetra_Export>& max_exporter,
    const std::map<int, unsigned>& my_mindofpernode) const
{
  std::vector<int> dofmapvec;
  std::map<int, std::vector<int>> dofs;

  const int* ngids = max_nodemap.MyGlobalElements();
  const int numnode = max_nodemap.NumMyElements();

  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* actnode = max_dis.g_node(ngids[i]);

    // check if the nodal GID is part of the mindofmap
    std::map<int, unsigned>::const_iterator pos = my_mindofpernode.find(ngids[i]);
    if (pos == my_mindofpernode.end())
      FOUR_C_THROW("The GID %d could not be found in the my_mindofpernode map!", ngids[i]);

    // get the number of dofs to copy
    const unsigned numdof = pos->second;
    const std::vector<int> dof = max_dis.dof(0, actnode);
    if (numdof > dof.size())
      FOUR_C_THROW("Got just %d DoF's at node %d (LID=%d) but expected at least %d", dof.size(),
          ngids[i], i, numdof);

    // copy the first numdof dofs
    std::copy(dof.data(), dof.data() + numdof, back_inserter(dofs[ngids[i]]));
    std::copy(dof.data(), dof.data() + numdof, back_inserter(dofmapvec));
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(), dofmapvec.end());
  if (pos != dofmapvec.end() and *pos < 0) FOUR_C_THROW("Illegal DoF number %d", *pos);

  // dof map is the original, unpermuted distribution of dofs
  max_dofmap =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), dofmapvec.data(), 0, max_dis.get_comm()));

  dofmapvec.clear();

  Core::Communication::Exporter exportdofs(max_nodemap, max_permnodemap, max_dis.get_comm());
  exportdofs.do_export(dofs);

  const int* permngids = max_permnodemap.MyGlobalElements();
  const int permnumnode = max_permnodemap.NumMyElements();

  for (int i = 0; i < permnumnode; ++i)
  {
    const std::vector<int>& dof = dofs[permngids[i]];
    std::copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  dofs.clear();

  // permuted dof map according to a given permuted node map
  max_permdofmap =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), dofmapvec.data(), 0, max_dis.get_comm()));

  /* prepare communication plan to create a dofmap out of a permuted
   * dof map */
  max_exporter = Teuchos::rcp<Epetra_Export>(new Epetra_Export(*max_permdofmap, *max_dofmap));
}

FOUR_C_NAMESPACE_CLOSE
