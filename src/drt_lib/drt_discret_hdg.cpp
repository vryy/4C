/*---------------------------------------------------------------------*/
/*!

\brief a class to manage an enhanced discretization including all faces for HDG

\level 2

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_discret_hdg.H"

#include "drt_exporter.H"
#include "drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_fluid_ele/fluid_ele_hdg.H"
#include "../drt_fluid_ele/fluid_ele_hdg_weak_comp.H"
#include "../drt_fluid_ele/fluid_ele_calc_hdg.H"
#include "../drt_fluid_ele/fluid_ele_calc_hdg_weak_comp.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_acou/acou_ele.H"
#include "../drt_acou/acou_ele_action.H"
#include "../drt_elemag/elemag_ele_action.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_mat/matpar_parameter.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/fluid_murnaghantait.H"
#include "../drt_mat/fluid_weakly_compressible.H"


DRT::DiscretizationHDG::DiscretizationHDG(const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DiscretizationFaces(name, comm)
{
  this->doboundaryfaces_ = true;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                     kronbichler 12/13|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationHDG::FillComplete(
    bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions)
{
  // call FillComleteFaces of base class with create_faces set to true
  this->FillCompleteFaces(assigndegreesoffreedom, initelements, doboundaryconditions, true);

  // get the correct face orientation from the owner. since the elements in general do not allow
  // packing, extract the node ids, communicate them, and change the node ids in the element
  Exporter nodeexporter(*facerowmap_, *facecolmap_, Comm());
  std::map<int, std::vector<int>> nodeIds, trafoMap;
  for (std::map<int, Teuchos::RCP<DRT::FaceElement>>::const_iterator f = faces_.begin();
       f != faces_.end(); ++f)
  {
    std::vector<int> ids(f->second->NumNode());
    for (int i = 0; i < f->second->NumNode(); ++i) ids[i] = f->second->NodeIds()[i];
    nodeIds[f->first] = ids;
    trafoMap[f->first] = f->second->GetLocalTrafoMap();
  }

  nodeexporter.Export(nodeIds);
  nodeexporter.Export(trafoMap);

  for (std::map<int, Teuchos::RCP<DRT::FaceElement>>::iterator f = faces_.begin();
       f != faces_.end(); ++f)
  {
    if (f->second->Owner() == Comm().MyPID()) continue;
    std::vector<int>& ids = nodeIds[f->first];
    dsassert(ids.size() > 0, "Lost a face during communication");
    f->second->SetNodeIds(ids.size(), &ids[0]);
    f->second->SetLocalTrafoMap(trafoMap[f->first]);

    // refresh node pointers if they have been set up
    DRT::Node** oldnodes = f->second->Nodes();
    if (oldnodes != 0)
    {
      std::vector<DRT::Node*> nodes(ids.size(), 0);

      for (unsigned int i = 0; i < ids.size(); ++i)
      {
        for (unsigned int j = 0; j < ids.size(); ++j)
          if (oldnodes[j]->Id() == ids[i])
          {
            nodes[i] = oldnodes[j];
          }
        dsassert(nodes[i] != 0, "Could not find node.");
      }
      f->second->BuildNodalPointers(&nodes[0]);
    }

    // check master/slave relation of current face in terms of the local trafo map
    dsassert(
        f->second->ParentMasterElement() != NULL, "Unexpected topology between face and parent");
    const int* nodeIdsMaster = f->second->ParentMasterElement()->NodeIds();
    const int* nodeIds = f->second->NodeIds();

    std::vector<std::vector<int>> faceNodeOrder =
        DRT::UTILS::getEleNodeNumberingFaces(f->second->ParentMasterElement()->Shape());

    bool exchangeMasterAndSlave = false;
    for (int i = 0; i < f->second->NumNode(); ++i)
    {
      // TODO (MK): check that this is enough also on periodic B.C. where the
      // node ids are different in any case...
      if (nodeIdsMaster[faceNodeOrder[f->second->FaceMasterNumber()][i]] != nodeIds[i])
        exchangeMasterAndSlave = true;
    }
    if (exchangeMasterAndSlave)
    {
      DRT::Element* faceMaster = f->second->ParentMasterElement();
      const int faceMasterNo = f->second->FaceMasterNumber();
      // new master element might be NULL on MPI computations
      f->second->SetParentMasterElement(f->second->ParentSlaveElement(),
          f->second->ParentSlaveElement() != NULL ? f->second->FaceSlaveNumber() : -1);
      f->second->SetParentSlaveElement(faceMaster, faceMasterNo);
    }
  }

  // add dofsets:
  // nds = 0 used for trace values
  // nds = 1 used for interior values
  // nds = 2 used for nodal ALE values
  if (this->lRowElement(0)->ElementType().Name() == "FluidHDGWeakCompType")
  {
    // add nds 1
    if (this->NumDofSets() == 1)
    {
      int ndof_ele = this->NumMyRowElements() > 0
                         ? dynamic_cast<DRT::ELEMENTS::FluidHDGWeakComp*>(this->lRowElement(0))
                               ->NumDofPerElementAuxiliary()
                         : 0;
      Teuchos::RCP<DRT::DofSetInterface> dofset_ele =
          Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0, ndof_ele, 0, false));

      this->AddDofSet(dofset_ele);
    }

    // add nds 2
    if (this->NumDofSets() == 2)
    {
      int ndof_node = this->NumMyRowElements() > 0
                          ? dynamic_cast<DRT::ELEMENTS::FluidHDGWeakComp*>(this->lRowElement(0))
                                ->NumDofPerNodeAuxiliary()
                          : 0;
      Teuchos::RCP<DRT::DofSetInterface> dofset_node =
          Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(ndof_node, 0, 0, false));

      this->AddDofSet(dofset_node);
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | AssignGlobalIDs                                        schoeder 06/14|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationHDG::AssignGlobalIDs(const Epetra_Comm& comm,
    const std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>& elementmap,
    std::map<int, Teuchos::RCP<DRT::Element>>& finalelements)
{
  // The point here is to make sure the element gid are the same on any
  // parallel distribution of the elements. Thus we allreduce thing to
  // processor 0 and sort the element descriptions (vectors of nodal ids)
  // there. We also communicate the element degree! This is the difference
  // the base class function!
  //
  // This routine has not been optimized for efficiency. I don't think that is
  // needed.
  //
  // pack elements on all processors

  int size = 0;
  std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>::const_iterator elemsiter;
  for (elemsiter = elementmap.begin(); elemsiter != elementmap.end(); ++elemsiter)
  {
    size += elemsiter->first.size() + 2;
  }

  std::vector<int> sendblock;
  sendblock.reserve(size);
  for (elemsiter = elementmap.begin(); elemsiter != elementmap.end(); ++elemsiter)
  {
    sendblock.push_back(elemsiter->first.size());
    sendblock.push_back(elemsiter->second->Degree());
    std::copy(elemsiter->first.begin(), elemsiter->first.end(), std::back_inserter(sendblock));
  }

  // communicate elements to processor 0

  int mysize = sendblock.size();
  comm.SumAll(&mysize, &size, 1);
  int mypos = LINALG::FindMyPos(sendblock.size(), comm);

  std::vector<int> send(size);
  std::fill(send.begin(), send.end(), 0);
  std::copy(sendblock.begin(), sendblock.end(), &send[mypos]);
  sendblock.clear();
  std::vector<int> recv(size);
  comm.SumAll(&send[0], &recv[0], size);
  send.clear();

  // unpack, unify and sort elements on processor 0

  if (comm.MyPID() == 0)
  {
    std::map<std::vector<int>, int> elementsanddegree;
    int index = 0;
    while (index < static_cast<int>(recv.size()))
    {
      int esize = recv[index];
      int degree = recv[index + 1];
      index += 2;
      std::vector<int> element;
      element.reserve(esize);
      std::copy(&recv[index], &recv[index + esize], std::back_inserter(element));
      index += esize;

      // check if we already have this and if so, check for max degree
      std::map<std::vector<int>, int>::const_iterator iter = elementsanddegree.find(element);
      if (iter != elementsanddegree.end())
      {
        degree = iter->second > degree ? iter->second : degree;
        elementsanddegree.erase(
            element);  // is only inserted in the next line, if the entry does not exist
      }
      elementsanddegree.insert(std::pair<std::vector<int>, int>(element, degree));
    }
    recv.clear();

    // pack again to distribute pack to all processors
    std::map<std::vector<int>, int>::const_iterator iter;
    send.reserve(index);

    for (iter = elementsanddegree.begin(); iter != elementsanddegree.end(); ++iter)
    {
      send.push_back(iter->first.size());
      send.push_back(iter->second);
      std::copy(iter->first.begin(), iter->first.end(), std::back_inserter(send));
    }
    size = send.size();
  }
  else
  {
    recv.clear();
  }

  // broadcast sorted elements to all processors

  comm.Broadcast(&size, 1, 0);
  send.resize(size);
  comm.Broadcast(&send[0], send.size(), 0);

  // Unpack sorted elements. Take element position for gid.

  int index = 0;
  int gid = 0;
  while (index < static_cast<int>(send.size()))
  {
    int esize = send[index];
    int degree = send[index + 1];
    index += 2;
    std::vector<int> element;
    element.reserve(esize);
    std::copy(&send[index], &send[index + esize], std::back_inserter(element));
    index += esize;

    // set gid to my elements
    std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>::const_iterator iter =
        elementmap.find(element);
    if (iter != elementmap.end())
    {
      iter->second->SetId(gid);
      // TODO visc eles, fluid hdg eles
      Teuchos::RCP<DRT::ELEMENTS::AcouIntFace> acouele =
          Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::AcouIntFace>(iter->second);
      if (acouele != Teuchos::null) acouele->SetDegree(degree);

      finalelements[gid] = iter->second;
    }

    gid += 1;
  }
}



/*----------------------------------------------------------------------*
 |  AddElementGhostLayer                               kronbichler 04/15|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationHDG::AddElementGhostLayer()
{
  if (!Filled()) dserror("Discretization must be filled upon entry of AddElementGhostLayer");

  const int mypid = Comm().MyPID();

  // step 1: identify all nodes on row elements
  std::vector<std::pair<int, int>> nodeids;
  for (int ele = 0; ele < NumMyRowElements(); ++ele)
    for (int n = 0; n < lRowElement(ele)->NumNode(); ++n)
      nodeids.push_back(std::make_pair(
          lRowElement(ele)->Nodes()[n]->Id(), lRowElement(ele)->Nodes()[n]->Owner()));

  // sort and compress out duplicates
  std::sort(nodeids.begin(), nodeids.end());
  nodeids.resize(std::unique(nodeids.begin(), nodeids.end()) - nodeids.begin());

  // step 2: Get patch of elements around each node. this info must be provided by the owner
  // of the node.

  // step 2a: create a map with data we want to import
  std::vector<int> indices;
  for (unsigned int i = 0; i < nodeids.size(); ++i)
    if (nodeids[i].second != mypid) indices.push_back(nodeids[i].first);
  const int numindices = indices.size();
  if (numindices == 0) indices.push_back(-1);
  Epetra_Map targetmap(-1, numindices, &indices[0], 0, Comm());

  // step 2b: create a copy of the element topology for owned nodes
  std::map<int, std::vector<int>> nodetoelement;
  for (int n = 0; n < NumMyRowNodes(); ++n)
  {
    std::vector<int>& elements = nodetoelement[lRowNode(n)->Id()];
    elements.resize(lRowNode(n)->NumElement());
    for (unsigned int e = 0; e < elements.size(); ++e)
      elements[e] = lRowNode(n)->Elements()[e]->Id();
  }

  // step 3: do the communication
  {
    Exporter exporter(*NodeRowMap(), targetmap, Comm());
    exporter.Export(nodetoelement);
  }

  // step 4: collect the ids of the new set of col elements
  std::vector<int> newcolelements;
  for (int i = 0; i < NumMyColElements(); ++i) newcolelements.push_back(lColElement(i)->Id());

  for (std::map<int, std::vector<int>>::const_iterator it = nodetoelement.begin();
       it != nodetoelement.end(); ++it)
    for (unsigned int e = 0; e < it->second.size(); ++e) newcolelements.push_back(it->second[e]);

  nodetoelement.clear();

  // sort and compress out duplicates
  std::sort(newcolelements.begin(), newcolelements.end());
  newcolelements.resize(
      std::unique(newcolelements.begin(), newcolelements.end()) - newcolelements.begin());
  const int numcolelements = newcolelements.size();
  if (numcolelements == 0) newcolelements.push_back(-1);
  Epetra_Map elecolmap(-1, numcolelements, &newcolelements[0], 0, Comm());

  // step 5: find node col map that matches the selected elements, similarly to elements
  std::map<int, std::vector<int>> elementtonode;
  for (int e = 0; e < NumMyRowElements(); ++e)
  {
    std::vector<int>& nodes = elementtonode[lRowElement(e)->Id()];
    nodes.resize(lRowElement(e)->NumNode());
    for (unsigned int n = 0; n < nodes.size(); ++n) nodes[n] = lRowElement(e)->Nodes()[n]->Id();
  }
  {
    Exporter exporter(*ElementRowMap(), elecolmap, Comm());
    exporter.Export(elementtonode);
  }
  std::vector<int> newcolnodes;
  for (std::map<int, std::vector<int>>::const_iterator it = elementtonode.begin();
       it != elementtonode.end(); ++it)
    for (unsigned int n = 0; n < it->second.size(); ++n) newcolnodes.push_back(it->second[n]);

  elementtonode.clear();

  // sort and compress out duplicates
  std::sort(newcolnodes.begin(), newcolnodes.end());
  newcolnodes.resize(std::unique(newcolnodes.begin(), newcolnodes.end()) - newcolnodes.begin());
  const int numcolnodes = newcolnodes.size();
  if (numcolnodes == 0) newcolnodes.push_back(-1);
  Epetra_Map nodecolmap(-1, numcolnodes, &newcolnodes[0], 0, Comm());

  // step 6: pass information about new col elements to the discretization (i.e., myself)
  ExportColumnNodes(nodecolmap);
  ExportColumnElements(elecolmap);
  FillComplete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationHDG& dis)
{
  // print standard discretization info
  dis.Print(os);
  // print additional info about internal faces
  dis.PrintFaces(os);

  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
    const DRT::Condition& cond, Epetra_Vector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  ReadDirichletCondition(face_discret, cond, toggle, dbcgids);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::ReadDirichletCondition(const DRT::DiscretizationFaces& discret,
    const DRT::Condition& cond, Epetra_Vector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const

{
  // call to corresponding method in base class; safety checks inside
  DRT::UTILS::Dbc::ReadDirichletCondition(discret, cond, toggle, dbcgids);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == NULL) return;

  // get onoff toggles
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");

  if (discret.NumMyRowFaces() > 0)
  {
    // initialize with true on each proc except proc 0
    bool pressureDone = discret.Comm().MyPID() != 0;

    // loop over all faces
    for (int i = 0; i < discret.NumMyRowFaces(); ++i)
    {
      const DRT::FaceElement* faceele = dynamic_cast<const DRT::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->NumDofPerFace(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff->size() <= component || (*onoff)[component] == 0) pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.DofRowMap(0)->LID(gid);

          // set toggle vector
          toggle[lid] = 1.0;
          // amend vector of DOF-IDs which are Dirichlet BCs
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).insert(gid);
          pressureDone = true;
        }
      }

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      int nummynodes = discret.lRowFace(i)->NumNode();
      const int* mynodes = discret.lRowFace(i)->NodeIds();
      for (int j = 0; j < nummynodes; ++j)
        if (!cond.ContainsNode(mynodes[j]))
        {
          faceRelevant = false;
          break;
        }
      if (!faceRelevant) continue;

      // get dofs of current face element
      std::vector<int> dofs = discret.Dof(0, discret.lRowFace(i));

      // loop over dofs
      for (unsigned int j = 0; j < dofperface; ++j)
      {
        // get global id
        const int gid = dofs[j];
        // get corresponding local id
        const int lid = toggle.Map().LID(gid);
        if (lid < 0)
          dserror(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        if ((*onoff)[onesetj] == 0)
        {
          // no DBC on this dof, set toggle zero
          toggle[lid] = 0.0;
          // get rid of entry in DBC map - if it exists
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).erase(gid);
          continue;
        }
        else  // if ((*onoff)[onesetj]==1)
        {
          // dof has DBC, set toggle vector one
          toggle[lid] = 1.0;
          // amend vector of DOF-IDs which are dirichlet BCs
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).insert(gid);
        }

      }  // loop over DOFs of face
    }    // loop over all faces
  }      // if there are faces

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::DoDirichletCondition(const DRT::DiscretizationInterface& discret,
    const DRT::Condition& cond, const double& time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  DoDirichletCondition(face_discret, cond, time, systemvectors, toggle);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::DoDirichletCondition(const DRT::DiscretizationFaces& discret,
    const DRT::Condition& cond, const double& time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle) const
{
  // call corresponding method from base class; safety checks inside
  DRT::UTILS::Dbc::DoDirichletCondition(discret, cond, time, systemvectors, toggle, NULL);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == NULL) return;

  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");

  // get curves, functs, vals, and onoff toggles from the condition
  const std::vector<int>* funct = cond.Get<std::vector<int>>("funct");
  const std::vector<double>* val = cond.Get<std::vector<double>>("val");
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvectors[0] != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvectors[0];
  }
  if (systemvectors[1] != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[1];
  }
  if (systemvectors[2] != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[2];
  }

  // do we have faces?
  if (discret.NumMyRowFaces() > 0)
  {
    Epetra_SerialDenseVector elevec1, elevec2, elevec3;
    Epetra_SerialDenseMatrix elemat1, elemat2;
    DRT::Element::LocationArray dummy(1);
    Teuchos::ParameterList initParams;
    if (DRT::Problem::Instance(0)->ProblemType() == prb_acou)
      initParams.set<int>("action", ACOU::project_dirich_field);
    else if (DRT::Problem::Instance(0)->ProblemType() == prb_elemag)
      initParams.set<int>("action", ELEMAG::project_dirich_field);
    else if (DRT::Problem::Instance(0)->ProblemType() == prb_scatra)
      initParams.set<int>("action", SCATRA::project_dirich_field);
    else
      initParams.set<int>(
          "action", FLD::project_fluid_field);  // TODO: Introduce a general action type that is
                                                // valid for all problems
    if (funct != NULL)
    {
      Teuchos::Array<int> functarray(*funct);
      initParams.set("funct", functarray);
    }
    Teuchos::Array<int> onoffarray(*onoff);
    initParams.set("onoff", onoffarray);
    initParams.set("time", time);

    // initialize with true if proc is not proc 0
    bool pressureDone = discret.Comm().MyPID() != 0;

    // loop over all faces
    for (int i = 0; i < discret.NumMyRowFaces(); ++i)
    {
      const DRT::FaceElement* faceele = dynamic_cast<const DRT::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->NumDofPerFace(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff->size() <= component || (*onoff)[component] == 0) pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.DofRowMap(0)->LID(gid);

          // amend vector of DOF-IDs which are Dirichlet BCs
          if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = 0.0;
          if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = 0.0;
          if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = 0.0;

          // --------------------------------------------------------------------------------------
          // get parameters
          Teuchos::ParameterList params = DRT::Problem::Instance()->FluidDynamicParams();

          // check whether the imposition of the average pressure is requested
          const int dopressavgbc =
              DRT::INPUT::IntegralValue<INPAR::FLUID::PressAvgBc>(params, "PRESSAVGBC");

          if (dopressavgbc == INPAR::FLUID::yes_pressure_average_bc)
          {
            double pressureavgBC = 0.0;

            // get 1st element
            DRT::Element* ele = discret.lRowElement(0);
            DRT::ELEMENTS::Fluid* fluidele = dynamic_cast<DRT::ELEMENTS::Fluid*>(ele);

            // get material
            Teuchos::RCP<MAT::Material> mat = ele->Material();

            // get discretization type
            const DRT::Element::DiscretizationType distype = ele->Shape();

            // evaluate pressure average     //TODO als make it valid for every discretization type
            Epetra_SerialDenseVector elevec = Epetra_SerialDenseVector(1);
            if (distype == DRT::Element::DiscretizationType::quad4)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad4>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::quad8)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad8>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::quad9)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad9>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::tri3)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::tri3>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::tri6)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::tri6>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else
              dserror("Given distype currently not implemented.");
            pressureavgBC = elevec[0];

            (*systemvectors[0])[lid] = pressureavgBC;

            std::cout << "\n-----------------------------------------------------------------------"
                         "-------------------"
                      << std::endl;
            std::cout << "| Warning: Imposing the analytical average pressure in the first element "
                         "as Dirichlet BC |"
                      << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "-----------------\n"
                      << std::endl;
          }
          // --------------------------------------------------------------------------------------
          pressureDone = true;
        }
      }
      int nummynodes = discret.lRowFace(i)->NumNode();
      const int* mynodes = discret.lRowFace(i)->NodeIds();

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      for (int j = 0; j < nummynodes; ++j)
        if (!cond.ContainsNode(mynodes[j]))
        {
          faceRelevant = false;
          break;
        }
      if (!faceRelevant) continue;

      initParams.set<unsigned int>(
          "faceconsider", static_cast<unsigned int>(faceele->FaceMasterNumber()));
      if (static_cast<unsigned int>(elevec1.M()) != dofperface) elevec1.Shape(dofperface, 1);
      std::vector<int> dofs = discret.Dof(0, discret.lRowFace(i));

      bool do_evaluate = false;
      if (funct != NULL)
        for (unsigned int i = 0; i < component; ++i)
          if ((*funct)[i] > 0) do_evaluate = true;

      if (do_evaluate)
      {
        // cast the const qualifier away, thus the Evaluate routine can be called.
        DRT::DiscretizationFaces& non_const_dis = const_cast<DRT::DiscretizationFaces&>(discret);
        faceele->ParentMasterElement()->Evaluate(
            initParams, non_const_dis, dummy, elemat1, elemat2, elevec1, elevec2, elevec3);
      }
      else
        for (unsigned int i = 0; i < dofperface; ++i) elevec1(i) = 1.;

      // loop over face dofs
      for (unsigned int j = 0; j < dofperface; ++j)
      {
        // get global id
        const int gid = dofs[j];
        // get corresponding local id
        const int lid = toggle.Map().LID(gid);
        if (lid < 0)
          dserror(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        // check whether dof gid is a dbc gid
        if (std::abs(toggle[lid] - 1.0) > 1e-13) continue;

        std::vector<double> value(deg + 1, (*val)[onesetj]);

        // assign value
        if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = value[0] * elevec1(j);
        if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = value[1] * elevec1(j);
        if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = value[2] * elevec1(j);

      }  // loop over all DOFs
    }    // loop over all faces

  }  // if there are faces

  return;
}
