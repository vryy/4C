/*----------------------------------------------------------------------------*/
/**
\file xfield_field_coupling.cpp

\brief Special coupling routines to handle the possible changing number of DoFs
       from node to node during XFEM simulations

\maintainer Michael Hiermeier

\date Sep 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xfield_field_coupling.H"

#include "../drt_lib/drt_exporter.H"

#include "Epetra_Export.h"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::XFieldField::Coupling::Coupling()
    : ADAPTER::Coupling(),
      isinit_(false),
      min_dof_dis_(min_dof_unknown)
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::Init(
    const enum MinDofDiscretization & min_dof_dis)
{
  min_dof_dis_ = min_dof_dis;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::XFieldField::Coupling::MasterToSlave(
    const Teuchos::RCP<const Epetra_Vector> & mv,
    const enum ::XFEM::MapType & map_type) const
{
  Teuchos::RCP<Epetra_Vector> sv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return ADAPTER::Coupling::MasterToSlave(mv);
      break;
    case XFEM::map_nodes:
      sv = Teuchos::rcp(new Epetra_Vector(*slavenodemap_));
      break;
  }

  MasterToSlave(mv,map_type,sv);
  return sv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::XFieldField::Coupling::SlaveToMaster(
    const Teuchos::RCP<const Epetra_Vector> & sv,
    const enum ::XFEM::MapType & map_type) const
{
  Teuchos::RCP<Epetra_Vector> mv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return ADAPTER::Coupling::SlaveToMaster(sv);
      break;
    case XFEM::map_nodes:
      mv = Teuchos::rcp(new Epetra_Vector(*masternodemap_));
      break;
  }

  SlaveToMaster(sv,map_type,mv);
  return mv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::XFieldField::Coupling::MasterToSlave(
    const Teuchos::RCP<const Epetra_MultiVector> & mv,
    const enum ::XFEM::MapType & map_type) const
{
  Teuchos::RCP<Epetra_MultiVector> sv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return ADAPTER::Coupling::MasterToSlave(mv);
      break;
    case XFEM::map_nodes:
      sv = Teuchos::rcp(new Epetra_MultiVector(*slavenodemap_,mv->NumVectors()));
      break;
  }

  MasterToSlave(mv,map_type,sv);
  return sv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::XFieldField::Coupling::SlaveToMaster(
    const Teuchos::RCP<const Epetra_MultiVector> & sv,
    const enum ::XFEM::MapType & map_type) const
{
  Teuchos::RCP<Epetra_MultiVector> mv = Teuchos::null;
  switch (map_type)
  {
    case XFEM::map_dofs:
      return ADAPTER::Coupling::SlaveToMaster(sv);
      break;
    case XFEM::map_nodes:
      mv = Teuchos::rcp(new Epetra_MultiVector(*masternodemap_,sv->NumVectors()));
      break;
  }

  SlaveToMaster(sv,map_type,mv);
  return mv;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::MasterToSlave(
    const Teuchos::RCP<const Epetra_MultiVector> & mv,
    const enum ::XFEM::MapType & map_type,
    Teuchos::RCP<Epetra_MultiVector> sv) const
{
  switch (map_type)
  {
    case XFEM::map_dofs:
    {
      return ADAPTER::Coupling::MasterToSlave(mv,sv);
      break;
    }
    case XFEM::map_nodes:
    {
#ifdef DEBUG
      if (not mv->Map().PointSameAs(*masternodemap_))
        dserror("master node map vector expected");
      if (not sv->Map().PointSameAs(*slavenodemap_))
        dserror("slave node map vector expected");
      if (sv->NumVectors()!=mv->NumVectors())
        dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

      Epetra_MultiVector perm(*permslavenodemap_,mv->NumVectors());
      std::copy(mv->Values(), mv->Values()+(mv->MyLength()*mv->NumVectors()), perm.Values());

      const int err = sv->Export(perm,*nodal_slaveexport_,Insert);
      if (err)
        dserror("Export to nodal slave distribution returned err=%d",err);
    } // end: case XFEM::MultiFieldMapExtractor::map_nodes
  } // end: switch (map_type)
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::SlaveToMaster(
    const Teuchos::RCP<const Epetra_MultiVector>& sv,
    const enum ::XFEM::MapType & map_type,
    Teuchos::RCP<Epetra_MultiVector> mv) const
{
  switch (map_type)
  {
    case XFEM::map_dofs:
    {
      return ADAPTER::Coupling::SlaveToMaster(sv,mv);
      break;
    }
    case XFEM::map_nodes:
    {
#ifdef DEBUG
      if (not mv->Map().PointSameAs(*masternodemap_))
        dserror("master node map vector expected");
      if (not sv->Map().PointSameAs(*slavenodemap_))
        dserror("slave node map vector expected");
      if (sv->NumVectors()!=mv->NumVectors())
        dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

      Epetra_MultiVector perm(*permmasternodemap_,sv->NumVectors());
      std::copy(sv->Values(), sv->Values()+(sv->MyLength()*sv->NumVectors()), perm.Values());

      const int err = mv->Export(perm,*nodal_masterexport_,Insert);
      if (err)
        dserror("Export to nodal master distribution returned err=%d",err);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::BuildDofMaps(
    const DRT::DiscretizationInterface& masterdis,
    const DRT::DiscretizationInterface& slavedis,
    const Teuchos::RCP<const Epetra_Map>& masternodemap,
    const Teuchos::RCP<const Epetra_Map>& slavenodemap,
    const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
    const Teuchos::RCP<const Epetra_Map>& permslavenodemap,
    const int& numdof)
{
  SaveNodeMaps(masternodemap,slavenodemap,permmasternodemap,permslavenodemap);

  // call base class implementation
  if (numdof!=-1)
  {
    ADAPTER::Coupling::BuildDofMaps(masterdis,slavedis,masternodemap,
        slavenodemap,permmasternodemap,permslavenodemap,numdof);
    return;
  }

  CheckInit();
  /* This map contains the information how many DoF's per node have to be
   * considered, i.e. the number of DoF's per node of the min-dof discretization.
   * The map-key is the corresponding max-dof discretization nodal coupling GID. */
  std::map<int,unsigned> my_mindofpernode;
  switch (MinDofDis())
  {
    case min_dof_slave:
    {
      BuildMinDofMaps(slavedis,*slavenodemap,*permslavenodemap,SlDofMapPtr(),
          PermutedSlDofMapPtr(),SlExporterPtr(),*masternodemap,my_mindofpernode);
      BuildMaxDofMaps(masterdis,*masternodemap,*permmasternodemap,MaDofMapPtr(),
          PermutedMaDofMapPtr(),MaExporterPtr(),my_mindofpernode);
      break;
    }
    case min_dof_master:
    {
      BuildMinDofMaps(masterdis,*masternodemap,*permmasternodemap,MaDofMapPtr(),
          PermutedMaDofMapPtr(),MaExporterPtr(),*slavenodemap,my_mindofpernode);
      BuildMaxDofMaps(slavedis,*slavenodemap,*permslavenodemap,SlDofMapPtr(),
          PermutedSlDofMapPtr(),SlExporterPtr(),my_mindofpernode);
      break;
    }
    case min_dof_unknown:
    {
      dserror("The discretization with the minimum number of DoF's \n"
          "per node is unknown or cannot be identified, since it \n"
          "changes from node to node. This case needs extra \n"
          "communication effort and is currently unsupported.");
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::SaveNodeMaps(
    const Teuchos::RCP<const Epetra_Map>& masternodemap,
    const Teuchos::RCP<const Epetra_Map>& slavenodemap,
    const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
    const Teuchos::RCP<const Epetra_Map>& permslavenodemap)
{
  masternodemap_     = masternodemap;
  slavenodemap_      = slavenodemap;
  permmasternodemap_ = permmasternodemap;
  permslavenodemap_  = permslavenodemap;

  nodal_masterexport_ = Teuchos::rcp<Epetra_Export>(
      new Epetra_Export(*permmasternodemap, *masternodemap));
  nodal_slaveexport_ = Teuchos::rcp<Epetra_Export>(
      new Epetra_Export(*permslavenodemap, *slavenodemap));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::BuildMinDofMaps(
    const DRT::DiscretizationInterface& min_dis,
    const Epetra_Map& min_nodemap,
    const Epetra_Map& min_permnodemap,
    Teuchos::RCP<const Epetra_Map>& min_dofmap,
    Teuchos::RCP<const Epetra_Map>& min_permdofmap,
    Teuchos::RCP<Epetra_Export>& min_exporter,
    const Epetra_Map& max_nodemap,
    std::map<int,unsigned>& my_mindofpernode) const
{
  std::vector<int> dofmapvec;
  std::map<int, std::vector<int> > dofs;

  const int* ngids = min_nodemap.MyGlobalElements();
  const int numnode = min_nodemap.NumMyElements();

  for (int i=0; i<numnode; ++i)
  {
    const DRT::Node* actnode = min_dis.gNode(ngids[i]);

    const int numdof = min_dis.NumDof(actnode);
    const std::vector<int> dof = min_dis.Dof(0,actnode);
    std::copy(&dof[0], &dof[0]+numdof, back_inserter(dofs[ngids[i]]));
    std::copy(&dof[0], &dof[0]+numdof, back_inserter(dofmapvec));
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(),
      dofmapvec.end());
  if (pos!=dofmapvec.end() and *pos < 0)
    dserror("Illegal DoF number %d", *pos);

  // dof map is the original, unpermuted distribution of dofs
  min_dofmap = Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0],
      0, min_dis.Comm()));

  dofmapvec.clear();

  DRT::Exporter exportdofs(min_nodemap,min_permnodemap,min_dis.Comm());
  exportdofs.Export(dofs);

  const int* permngids = min_permnodemap.MyGlobalElements();
  const int permnumnode = min_permnodemap.NumMyElements();

  for (int i=0; i<permnumnode; ++i)
  {
    const std::vector<int>& dof = dofs[permngids[i]];
    std::copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  /* -------------------------------------------------------------------------
   * Get the number of dofs per node which have to be considered at the
   * coupling GID in the max-dof-discretization
   * -------------------------------------------------------------------------*/
  std::map<int, std::vector<int> >::const_iterator cit;
  for (cit=dofs.begin();cit!=dofs.end();++cit)
  {
    const int min_permlid = min_permnodemap.LID(cit->first);
    const int max_gid     = max_nodemap.GID(min_permlid);
    my_mindofpernode[max_gid] = cit->second.size();
  }
  dofs.clear();

  // permuted dof map according to a given permuted node map
  min_permdofmap = Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0],
      0, min_dis.Comm()));

  /* prepare communication plan to create a dofmap out of a permuted
   * dof map */
  min_exporter = Teuchos::rcp<Epetra_Export>(new Epetra_Export(*min_permdofmap,
      *min_dofmap));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::Coupling::BuildMaxDofMaps(
    const DRT::DiscretizationInterface& max_dis,
    const Epetra_Map& max_nodemap,
    const Epetra_Map& max_permnodemap,
    Teuchos::RCP<const Epetra_Map>& max_dofmap,
    Teuchos::RCP<const Epetra_Map>& max_permdofmap,
    Teuchos::RCP<Epetra_Export>& max_exporter,
    const std::map<int,unsigned>& my_mindofpernode) const
{
  std::vector<int> dofmapvec;
  std::map<int, std::vector<int> > dofs;

  const int* ngids = max_nodemap.MyGlobalElements();
  const int numnode = max_nodemap.NumMyElements();

  for (int i=0; i<numnode; ++i)
  {
    const DRT::Node* actnode = max_dis.gNode(ngids[i]);

    // check if the nodal GID is part of the mindofmap
    std::map<int,unsigned>::const_iterator pos = my_mindofpernode.find(ngids[i]);
    if (pos == my_mindofpernode.end())
      dserror("The GID %d could not be found in the my_mindofpernode map!",ngids[i]);

    // get the number of dofs to copy
    const unsigned numdof = pos->second;
    const std::vector<int> dof = max_dis.Dof(0,actnode);
    if (numdof > dof.size())
      dserror("Got just %d DoF's at node %d (LID=%d) but expected at least %d",
          dof.size(),ngids[i],i,numdof);

    // copy the first numdof dofs
    std::copy(&dof[0], &dof[0]+numdof, back_inserter(dofs[ngids[i]]));
    std::copy(&dof[0], &dof[0]+numdof, back_inserter(dofmapvec));
  }

  std::vector<int>::const_iterator pos = std::min_element(dofmapvec.begin(),
      dofmapvec.end());
  if (pos!=dofmapvec.end() and *pos < 0)
    dserror("Illegal DoF number %d", *pos);

  // dof map is the original, unpermuted distribution of dofs
  max_dofmap = Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(), &dofmapvec[0],
      0, max_dis.Comm()));

  dofmapvec.clear();

  DRT::Exporter exportdofs(max_nodemap,max_permnodemap,max_dis.Comm());
  exportdofs.Export(dofs);

  const int * permngids = max_permnodemap.MyGlobalElements();
  const int permnumnode = max_permnodemap.NumMyElements();

  for (int i=0; i<permnumnode; ++i)
  {
    const std::vector<int>& dof = dofs[permngids[i]];
    std::copy(dof.begin(), dof.end(), back_inserter(dofmapvec));
  }

  dofs.clear();

  // permuted dof map according to a given permuted node map
  max_permdofmap = Teuchos::rcp(new Epetra_Map(-1, dofmapvec.size(),
      &dofmapvec[0], 0, max_dis.Comm()));

  /* prepare communication plan to create a dofmap out of a permuted
   * dof map */
  max_exporter = Teuchos::rcp<Epetra_Export>(new Epetra_Export(*max_permdofmap,
      *max_dofmap));
}
