/*---------------------------------------------------------------------*/
/*! \file

\brief a class to manage an enhanced discretization including all faces for HDG

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_lib_discret_hdg.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_dg_element.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_utils_discret.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::DiscretizationHDG::DiscretizationHDG(const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DiscretizationFaces(name, comm)
{
  this->doboundaryfaces_ = true;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                     kronbichler 12/13|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationHDG::fill_complete(
    bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions)
{
  // call FillComleteFaces of base class with create_faces set to true
  this->FillCompleteFaces(assigndegreesoffreedom, initelements, doboundaryconditions, true);

  // get the correct face orientation from the owner. since the elements in general do not allow
  // packing, extract the node ids, communicate them, and change the node ids in the element
  CORE::COMM::Exporter nodeexporter(*facerowmap_, *facecolmap_, Comm());
  std::map<int, std::vector<int>> nodeIds, trafoMap;
  for (std::map<int, Teuchos::RCP<CORE::Elements::FaceElement>>::const_iterator f = faces_.begin();
       f != faces_.end(); ++f)
  {
    std::vector<int> ids(f->second->num_node());
    for (int i = 0; i < f->second->num_node(); ++i) ids[i] = f->second->NodeIds()[i];
    nodeIds[f->first] = ids;
    trafoMap[f->first] = f->second->GetLocalTrafoMap();
  }

  nodeexporter.Export(nodeIds);
  nodeexporter.Export(trafoMap);

  for (std::map<int, Teuchos::RCP<CORE::Elements::FaceElement>>::iterator f = faces_.begin();
       f != faces_.end(); ++f)
  {
    if (f->second->Owner() == Comm().MyPID()) continue;
    std::vector<int>& ids = nodeIds[f->first];
    FOUR_C_ASSERT(ids.size() > 0, "Lost a face during communication");
    f->second->SetNodeIds(ids.size(), ids.data());
    f->second->set_local_trafo_map(trafoMap[f->first]);

    // refresh node pointers if they have been set up
    DRT::Node** oldnodes = f->second->Nodes();
    if (oldnodes != nullptr)
    {
      std::vector<DRT::Node*> nodes(ids.size(), nullptr);

      for (unsigned int i = 0; i < ids.size(); ++i)
      {
        for (unsigned int j = 0; j < ids.size(); ++j)
          if (oldnodes[j]->Id() == ids[i])
          {
            nodes[i] = oldnodes[j];
          }
        FOUR_C_ASSERT(nodes[i] != 0, "Could not find node.");
      }
      f->second->BuildNodalPointers(nodes.data());
    }

    // check master/slave relation of current face in terms of the local trafo map
    FOUR_C_ASSERT(
        f->second->ParentMasterElement() != nullptr, "Unexpected topology between face and parent");
    const int* nodeIdsMaster = f->second->ParentMasterElement()->NodeIds();
    const int* nodeIds = f->second->NodeIds();

    std::vector<std::vector<int>> faceNodeOrder =
        CORE::FE::getEleNodeNumberingFaces(f->second->ParentMasterElement()->Shape());

    bool exchangeMasterAndSlave = false;
    for (int i = 0; i < f->second->num_node(); ++i)
    {
      // TODO (MK): check that this is enough also on periodic B.C. where the
      // node ids are different in any case...
      if (nodeIdsMaster[faceNodeOrder[f->second->FaceMasterNumber()][i]] != nodeIds[i])
        exchangeMasterAndSlave = true;
    }
    if (exchangeMasterAndSlave)
    {
      CORE::Elements::Element* faceMaster = f->second->ParentMasterElement();
      const int faceMasterNo = f->second->FaceMasterNumber();
      // new master element might be nullptr on MPI computations
      f->second->set_parent_master_element(f->second->ParentSlaveElement(),
          f->second->ParentSlaveElement() != nullptr ? f->second->FaceSlaveNumber() : -1);
      f->second->set_parent_slave_element(faceMaster, faceMasterNo);
    }
  }

  // add dofsets:
  // nds = 0 used for trace values
  // nds = 1 used for interior values
  // nds = 2 used for nodal ALE values
  if (this->NumMyRowElements())
    if (this->lRowElement(0)->ElementType().Name() == "FluidHDGWeakCompType")
    {
      // add nds 1
      if (this->NumDofSets() == 1)
      {
        int ndof_ele = this->NumMyRowElements() > 0
                           ? dynamic_cast<CORE::Elements::DgElement*>(this->lRowElement(0))
                                 ->num_dof_per_element_auxiliary()
                           : 0;
        Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofset_ele =
            Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(0, ndof_ele, 0, false));

        this->AddDofSet(dofset_ele);
      }

      // add nds 2
      if (this->NumDofSets() == 2)
      {
        int ndof_node = this->NumMyRowElements() > 0
                            ? dynamic_cast<CORE::Elements::DgElement*>(this->lRowElement(0))
                                  ->num_dof_per_node_auxiliary()
                            : 0;
        Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofset_node =
            Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(ndof_node, 0, 0, false));

        this->AddDofSet(dofset_node);
      }
    }

  return 0;
}


/*----------------------------------------------------------------------*
 | assign_global_i_ds                                        schoeder 06/14|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationHDG::assign_global_i_ds(const Epetra_Comm& comm,
    const std::map<std::vector<int>, Teuchos::RCP<CORE::Elements::Element>>& elementmap,
    std::map<int, Teuchos::RCP<CORE::Elements::Element>>& finalelements)
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
  std::map<std::vector<int>, Teuchos::RCP<CORE::Elements::Element>>::const_iterator elemsiter;
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
  int mypos = CORE::LINALG::FindMyPos(sendblock.size(), comm);

  std::vector<int> send(size);
  std::fill(send.begin(), send.end(), 0);
  std::copy(sendblock.begin(), sendblock.end(), &send[mypos]);
  sendblock.clear();
  std::vector<int> recv(size);
  comm.SumAll(send.data(), recv.data(), size);
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
  comm.Broadcast(send.data(), send.size(), 0);

  // Unpack sorted elements. Take element position for gid.

  int index = 0;
  int gid = 0;
  while (index < static_cast<int>(send.size()))
  {
    int esize = send[index];
    index += 2;
    std::vector<int> element;
    element.reserve(esize);
    std::copy(&send[index], &send[index + esize], std::back_inserter(element));
    index += esize;

    // set gid to my elements
    std::map<std::vector<int>, Teuchos::RCP<CORE::Elements::Element>>::const_iterator iter =
        elementmap.find(element);
    if (iter != elementmap.end())
    {
      iter->second->SetId(gid);

      finalelements[gid] = iter->second;
    }

    gid += 1;
  }
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
void DRT::UTILS::DbcHDG::read_dirichlet_condition(const DRT::Discretization& discret,
    const CORE::Conditions::Condition& cond, double time, DRT::UTILS::Dbc::DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  read_dirichlet_condition(face_discret, cond, time, info, dbcgids, hierarchical_order);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::read_dirichlet_condition(const DRT::DiscretizationFaces& discret,
    const CORE::Conditions::Condition& cond, double time, DRT::UTILS::Dbc::DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const

{
  // call to corresponding method in base class; safety checks inside
  DRT::UTILS::Dbc::read_dirichlet_condition(discret, cond, time, info, dbcgids, hierarchical_order);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == nullptr) return;

  // get onoff toggles
  const auto& onoff = cond.parameters().Get<std::vector<int>>("onoff");

  if (discret.NumMyRowFaces() > 0)
  {
    // initialize with true on each proc except proc 0
    bool pressureDone = discret.Comm().MyPID() != 0;

    // loop over all faces
    for (int i = 0; i < discret.NumMyRowFaces(); ++i)
    {
      const CORE::Elements::FaceElement* faceele =
          dynamic_cast<const CORE::Elements::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->num_dof_per_face(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff.size() <= component || onoff[component] == 0 ||
          GLOBAL::Problem::Instance(0)->GetProblemType() != GLOBAL::ProblemType::fluid)
        pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.dof_row_map(0)->LID(gid);

          // set toggle vector
          info.toggle[lid] = 1;
          // amend vector of DOF-IDs which are Dirichlet BCs
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).insert(gid);
          pressureDone = true;
        }
      }

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      int nummynodes = discret.lRowFace(i)->num_node();
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
        const int lid = info.toggle.Map().LID(gid);
        if (lid < 0)
          FOUR_C_THROW(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        if (onoff[onesetj] == 0)
        {
          // no DBC on this dof, set toggle zero
          info.toggle[lid] = 0;
          // get rid of entry in DBC map - if it exists
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).erase(gid);
          continue;
        }
        else  // if ((*onoff)[onesetj]==1)
        {
          // dof has DBC, set toggle vector one
          info.toggle[lid] = 1;
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
void DRT::UTILS::DbcHDG::do_dirichlet_condition(const DRT::Discretization& discret,
    const CORE::Conditions::Condition& cond, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  do_dirichlet_condition(face_discret, cond, time, systemvectors, toggle);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::DbcHDG::do_dirichlet_condition(const DRT::DiscretizationFaces& discret,
    const CORE::Conditions::Condition& cond, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle) const
{
  // call corresponding method from base class; safety checks inside
  DRT::UTILS::Dbc::do_dirichlet_condition(discret, cond, time, systemvectors, toggle, nullptr);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == nullptr) return;

  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.GetNodes();
  if (!nodeids) FOUR_C_THROW("Dirichlet condition does not have nodal cloud");

  // get curves, functs, vals, and onoff toggles from the condition
  const auto* funct = &cond.parameters().Get<std::vector<int>>("funct");
  const auto* val = &cond.parameters().Get<std::vector<double>>("val");
  const auto* onoff = &cond.parameters().Get<std::vector<int>>("onoff");

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
    CORE::LINALG::SerialDenseVector elevec1, elevec2, elevec3;
    CORE::LINALG::SerialDenseMatrix elemat1, elemat2;
    CORE::Elements::Element::LocationArray dummy(1);
    Teuchos::ParameterList initParams;
    if (GLOBAL::Problem::Instance(0)->GetProblemType() == GLOBAL::ProblemType::elemag or
        GLOBAL::Problem::Instance(0)->GetProblemType() == GLOBAL::ProblemType::scatra)
    {
      initParams.set("hdg_action", true);
      CORE::UTILS::AddEnumClassToParameterList<DRT::HDGAction>(
          "action", DRT::HDGAction::project_dirich_field, initParams);
    }

    // TODO: Introduce a general action type that is
    // valid for all problems
    if (funct != nullptr)
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
      const CORE::Elements::FaceElement* faceele =
          dynamic_cast<const CORE::Elements::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->num_dof_per_face(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff->size() <= component || (*onoff)[component] == 0 ||
          GLOBAL::Problem::Instance(0)->GetProblemType() != GLOBAL::ProblemType::fluid)
        pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.dof_row_map(0)->LID(gid);

          // amend vector of DOF-IDs which are Dirichlet BCs
          if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = 0.0;
          if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = 0.0;
          if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = 0.0;

          // --------------------------------------------------------------------------------------
          pressureDone = true;
        }
      }
      int nummynodes = discret.lRowFace(i)->num_node();
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
      if (static_cast<unsigned int>(elevec1.numRows()) != dofperface) elevec1.shape(dofperface, 1);
      std::vector<int> dofs = discret.Dof(0, discret.lRowFace(i));

      bool do_evaluate = false;
      if (funct != nullptr)
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
          FOUR_C_THROW(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        // check whether dof gid is a dbc gid
        if (toggle[lid] == 0) continue;

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

FOUR_C_NAMESPACE_CLOSE
