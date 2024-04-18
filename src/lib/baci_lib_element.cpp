/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of base element class in BACI with basic operations

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_element.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_discretization_geometric_search_bounding_volume.hpp"
#include "baci_discretization_geometric_search_params.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_condition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element_append_visualization.hpp"
#include "baci_lib_node.hpp"
#include "baci_mat_material.hpp"
#include "baci_utils_exceptions.hpp"

#include <Shards_BasicTopologies.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ShardsKeyToDisType(const unsigned& key)
{
  CORE::FE::CellType distype = CORE::FE::CellType::dis_none;
  switch (key)
  {
    case shards::Particle::key:
    {
      distype = CORE::FE::CellType::point1;
      break;
    }
    case shards::Line<2>::key:
    {
      distype = CORE::FE::CellType::line2;
      break;
    }
    case shards::Line<3>::key:
    {
      distype = CORE::FE::CellType::line3;
      break;
    }
    case shards::Quadrilateral<4>::key:
    {
      distype = CORE::FE::CellType::quad4;
      break;
    }
    case shards::Quadrilateral<8>::key:
    {
      distype = CORE::FE::CellType::quad8;
      break;
    }
    case shards::Quadrilateral<9>::key:
    {
      distype = CORE::FE::CellType::quad9;
      break;
    }
    case shards::Triangle<3>::key:
    {
      distype = CORE::FE::CellType::tri3;
      break;
    }
    case shards::Triangle<6>::key:
    {
      distype = CORE::FE::CellType::tri6;
      break;
    }
    case shards::Hexahedron<8>::key:
    {
      distype = CORE::FE::CellType::hex8;
      break;
    }
    case shards::Hexahedron<20>::key:
    {
      distype = CORE::FE::CellType::hex20;
      break;
    }
    case shards::Hexahedron<27>::key:
    {
      distype = CORE::FE::CellType::hex27;
      break;
    }
    case shards::Tetrahedron<4>::key:
    {
      distype = CORE::FE::CellType::tet4;
      break;
    }
    case shards::Tetrahedron<10>::key:
    {
      distype = CORE::FE::CellType::tet10;
      break;
    }
    case shards::Wedge<6>::key:
    {
      distype = CORE::FE::CellType::wedge6;
      break;
    }
    case shards::Wedge<15>::key:
    {
      distype = CORE::FE::CellType::wedge15;
      break;
    }
    case shards::Pyramid<5>::key:
    {
      distype = CORE::FE::CellType::pyramid5;
      break;
    }
    default:
      FOUR_C_THROW("Unknown conversion from Shards::key to disType!");
      break;
  }
  return distype;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(int id, int owner)
    : ParObject(), id_(id), lid_(-1), owner_(owner), mat_(1, Teuchos::null), is_nurbs_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(const DRT::Element& old)
    : ParObject(old),
      id_(old.id_),
      lid_(old.lid_),
      owner_(old.owner_),
      nodeid_(old.nodeid_),
      node_(old.node_),
      face_(old.face_),
      mat_(1, Teuchos::null),
      is_nurbs_(old.is_nurbs_)
{
  // we do NOT want a deep copy of the condition_ as the condition
  // is only a reference in the elements anyway
  std::map<std::string, Teuchos::RCP<Condition>>::const_iterator fool;
  for (fool = old.condition_.begin(); fool != old.condition_.end(); ++fool)
    SetCondition(fool->first, fool->second);

  if (!old.mat_.empty())
  {
    mat_.resize(old.mat_.size());
    for (unsigned iter = 0; iter < old.mat_.size(); ++iter)
      if (old.mat_[iter] != Teuchos::null) mat_[iter] = (old.mat_[iter]->Clone());
  }
  else
    mat_[0] = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::Element& element)
{
  element.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::Print(std::ostream& os) const
{
  os << std::setw(12) << Id() << " Owner " << std::setw(5) << Owner() << " ";
  const int nnode = NumNode();
  const int* nodeids = NodeIds();
  if (nnode > 0)
  {
    os << " Nodes ";
    for (int i = 0; i < nnode; ++i) os << std::setw(10) << nodeids[i] << " ";
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::Element::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  FOUR_C_THROW("subclass implementations missing");
  return false;
}


/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::SetNodeIds(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i = 0; i < nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Element::SetNodeIds(const std::string& distype, INPUT::LineDefinition* linedef)
{
  linedef->ExtractIntVector(distype, nodeid_);
  for (int& i : nodeid_) i -= 1;
  node_.resize(0);
}


/*----------------------------------------------------------------------*
 |  create material class (public)                                 05/07|
 *----------------------------------------------------------------------*/
void DRT::Element::SetMaterial(int matnum)
{
  Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matnum);
  if (mat == Teuchos::null)
    FOUR_C_THROW(
        "Invalid material given to the element. \n"
        "Invalid are Summands of the Elasthyper-Toolbox and single Growth-Materials. \n"
        "If you like to use a Summand of the Elasthyper-Material define it via MAT_ElastHyper. \n"
        "If you like to use a Growth-Material define it via the according base material.");
  mat_[0] = mat;
}

void DRT::Element::SetMaterial(int index, Teuchos::RCP<MAT::Material> mat)
{
  if (NumMaterial() > index)
    mat_[index] = mat;
  else if (NumMaterial() == index)
    AddMaterial(mat);
  else
    FOUR_C_THROW(
        "Setting material at index %d not possible (neither overwrite nor append) since currently  "
        "only %d materials are stored",
        index, NumMaterial());
}

/*----------------------------------------------------------------------*
 |  add material to element (public)                          vuong 02/14|
 *----------------------------------------------------------------------*/
int DRT::Element::AddMaterial(const Teuchos::RCP<MAT::Material>& mat)
{
  mat_.push_back(mat);

  return mat_.size();
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add id
  AddtoPack(data, id_);
  // add owner
  AddtoPack(data, owner_);
  // add vector nodeid_
  AddtoPack(data, nodeid_);
  // add material
  if (mat_[0] != Teuchos::null)
  {
    // pack only first material
    mat_[0]->Pack(data);
  }
  else
  {
    int size = 0;
    AddtoPack(data, size);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // id_
  ExtractfromPack(position, data, id_);
  // owner_
  ExtractfromPack(position, data, owner_);
  // nodeid_
  ExtractfromPack(position, data, nodeid_);
  // mat_
  std::vector<char> tmp;
  ExtractfromPack(position, data, tmp);
  if (!tmp.empty())
  {
    CORE::COMM::ParObject* o = CORE::COMM::Factory(tmp);
    auto* mat = dynamic_cast<MAT::Material*>(o);
    if (mat == nullptr) FOUR_C_THROW("failed to unpack material");
    // unpack only first material
    mat_[0] = Teuchos::rcp(mat);
  }
  else
  {
    mat_[0] = Teuchos::null;
  }

  // node_, face_, parent_master_, parent_slave_ are NOT communicated
  node_.resize(0);
  if (!face_.empty())
  {
    std::vector<Teuchos::RCP<DRT::FaceElement>> empty;
    std::swap(face_, empty);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(std::map<int, Teuchos::RCP<DRT::Node>>& nodes)
{
  int nnode = NumNode();
  const int* nodeids = NodeIds();
  node_.resize(nnode);
  for (int i = 0; i < nnode; ++i)
  {
    std::map<int, Teuchos::RCP<DRT::Node>>::const_iterator curr = nodes.find(nodeids[i]);
    // this node is not on this proc
    if (curr == nodes.end())
      FOUR_C_THROW("Element %d cannot find node %d", Id(), nodeids[i]);
    else
      node_[i] = curr->second.get();
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(DRT::Node** nodes)
{
  node_.resize(NumNode());
  for (int i = 0; i < NumNode(); ++i) node_[i] = nodes[i];
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal connectivity and weight nodes and edges        (public) |
 |                                                          ghamm 09/13 |
 *----------------------------------------------------------------------*/
void DRT::Element::NodalConnectivity(
    CORE::LINALG::SerialDenseMatrix& edgeweights, CORE::LINALG::SerialDenseVector& nodeweights)
{
  // weight for this element
  double weight = EvaluationCost();

  int numnode = NumNode();
  nodeweights.size(numnode);
  edgeweights.shape(numnode, numnode);

  // initialize weights
  for (int n = 0; n < numnode; ++n)
  {
    nodeweights[n] = weight;
    for (int k = 0; k < numnode; ++k)
    {
      edgeweights(n, k) = 1.0;
    }
  }

  // put squared weight on edges
  weight *= weight;

  std::vector<std::vector<int>> lines = CORE::FE::getEleNodeNumberingLines(Shape());
  size_t nodesperline = lines[0].size();
  if (nodesperline == 2)
  {
    for (auto& line : lines)
    {
      edgeweights(line[0], line[1]) = weight;
      edgeweights(line[1], line[0]) = weight;
    }
  }
  else if (nodesperline == 3)
  {
    for (auto& line : lines)
    {
      edgeweights(line[0], line[1]) = weight;
      edgeweights(line[1], line[0]) = weight;

      edgeweights(line[1], line[2]) = weight;
      edgeweights(line[2], line[1]) = weight;
    }
  }
  else
    FOUR_C_THROW("implementation is missing for this distype (%s)",
        CORE::FE::CellTypeToString(Shape()).c_str());

  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  auto startit = condition_.lower_bound(name);
  auto endit = condition_.upper_bound(name);
  int count = 0;
  std::multimap<std::string, Teuchos::RCP<Condition>>::const_iterator curr;
  for (curr = startit; curr != endit; ++curr) out[count++] = curr->second.get();
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Element::GetCondition(const std::string& name) const
{
  auto curr = condition_.find(name);
  if (curr == condition_.end()) return nullptr;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const DRT::Discretization& dis, const std::vector<int>& nds,
    DRT::Element::LocationArray& la, bool doDirichlet) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  if (numnode != static_cast<int>(nds.size()))
  {
    FOUR_C_THROW("wrong number of nodes");
  }

  la.Clear();

  // we need to look at all DofSets of our Discretization
  for (int dofset = 0; dofset < la.Size(); ++dofset)
  {
    std::vector<int>& lm = la[dofset].lm_;
    std::vector<int>& lmdirich = la[dofset].lmdirich_;
    std::vector<int>& lmowner = la[dofset].lmowner_;
    std::vector<int>& lmstride = la[dofset].stride_;

    // fill the vector with nodal dofs
    if (nodes)
    {
      for (int i = 0; i < numnode; ++i)
      {
        const DRT::Node* node = nodes[i];

        const int owner = node->Owner();
        std::vector<int> dof;
        dis.Dof(dof, node, dofset, nds[i]);
        const int size = dof.size();
        if (size) lmstride.push_back(size);

        for (int j = 0; j < size; ++j)
        {
          lmowner.push_back(owner);
          lm.push_back(dof[j]);
        }

        if (doDirichlet)
        {
          const std::vector<int>* flag = nullptr;
          DRT::Condition* dirich = node->GetCondition("Dirichlet");
          if (dirich)
          {
            if (dirich->Type() != DRT::Condition::PointDirichlet &&
                dirich->Type() != DRT::Condition::LineDirichlet &&
                dirich->Type() != DRT::Condition::SurfaceDirichlet &&
                dirich->Type() != DRT::Condition::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = dirich->Get<std::vector<int>>("onoff");
          }
          for (unsigned j = 0; j < dof.size(); ++j)
          {
            if (flag && (*flag)[j])
              lmdirich.push_back(1);
            else
              lmdirich.push_back(0);
          }
        }
      }
    }

    // fill the vector with element dofs
    const int owner = Owner();
    std::vector<int> dof = dis.Dof(dofset, this);
    if (!dof.empty()) lmstride.push_back(dof.size());
    for (int j : dof)
    {
      lmowner.push_back(owner);
      lm.push_back(j);
    }

    // fill the vector with face dofs
    if (this->NumDofPerFace(0) > 0)
    {
      for (int i = 0; i < NumFace(); ++i)
      {
        const int owner = face_[i]->Owner();
        std::vector<int> dof = dis.Dof(dofset, face_[i].getRawPtr());
        if (!dof.empty()) lmstride.push_back(dof.size());
        for (int j : dof)
        {
          lmowner.push_back(owner);
          lm.push_back(j);
        }
      }
    }

    if (doDirichlet)
    {
      const std::vector<int>* flag = nullptr;
      DRT::Condition* dirich = GetCondition("Dirichlet");
      if (dirich)
      {
        if (dirich->Type() != DRT::Condition::PointDirichlet &&
            dirich->Type() != DRT::Condition::LineDirichlet &&
            dirich->Type() != DRT::Condition::SurfaceDirichlet &&
            dirich->Type() != DRT::Condition::VolumeDirichlet)
          FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        flag = dirich->Get<std::vector<int>>("onoff");
      }
      for (unsigned j = 0; j < dof.size(); ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(
    const Discretization& dis, LocationArray& la, bool doDirichlet) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  la.Clear();

  // we need to look at all DofSets of our Discretization
  for (int dofset = 0; dofset < la.Size(); ++dofset)
  {
    std::vector<int>& lm = la[dofset].lm_;
    std::vector<int>& lmdirich = la[dofset].lmdirich_;
    std::vector<int>& lmowner = la[dofset].lmowner_;
    std::vector<int>& lmstride = la[dofset].stride_;

    // fill the vector with nodal dofs
    if (nodes)
    {
      for (int i = 0; i < numnode; ++i)
      {
        const DRT::Node* node = nodes[i];

        const int owner = node->Owner();
        std::vector<int> dof;
        dis.Dof(dof, node, dofset, 0, this);

        // if there are more dofs on the node than the element can handle, this cannot work
        FOUR_C_ASSERT(NumDofPerNode(*node) <= (int)dof.size() or dofset != 0,
            "More dofs on node than element can handle! Internal error!");

        // assume that the first dofs are the relevant ones
        const int size = dofset == 0 ? NumDofPerNode(*node) : dof.size();

        if (size) lmstride.push_back(size);
        for (int j = 0; j < size; ++j)
        {
          lmowner.push_back(owner);
          lm.push_back(dof[j]);
        }

        if (doDirichlet)
        {
          const std::vector<int>* flag = nullptr;
          DRT::Condition* dirich = node->GetCondition("Dirichlet");
          if (dirich)
          {
            if (dirich->Type() != DRT::Condition::PointDirichlet &&
                dirich->Type() != DRT::Condition::LineDirichlet &&
                dirich->Type() != DRT::Condition::SurfaceDirichlet &&
                dirich->Type() != DRT::Condition::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = dirich->Get<std::vector<int>>("onoff");
          }
          for (int j = 0; j < size; ++j)
          {
            if (flag && (*flag)[j])
              lmdirich.push_back(1);
            else
              lmdirich.push_back(0);
          }
        }
      }
    }

    // fill the vector with element dofs
    const int owner = Owner();
    std::vector<int> dof = dis.Dof(dofset, this);
    if (!dof.empty()) lmstride.push_back(dof.size());
    for (int j : dof)
    {
      lmowner.push_back(owner);
      lm.push_back(j);
    }

    // fill the vector with face dofs
    if (this->NumDofPerFace(0) > 0)
    {
      for (int i = 0; i < NumFace(); ++i)
      {
        const int owner = face_[i]->Owner();
        std::vector<int> dof = dis.Dof(dofset, face_[i].getRawPtr());
        if (!dof.empty()) lmstride.push_back(dof.size());
        for (int j : dof)
        {
          lmowner.push_back(owner);
          lm.push_back(j);
        }

        if (doDirichlet)
        {
          std::vector<DRT::Condition*> dirich_vec;
          dis.GetCondition("Dirichlet", dirich_vec);
          DRT::Condition* dirich;
          bool dirichRelevant = false;
          // Check if there exist a dirichlet condition
          if (!dirich_vec.empty())
          {
            // do only faces where all nodes are present in the node list
            const int nummynodes = face_[i]->NumNode();
            const int* mynodes = face_[i]->NodeIds();
            // Check if the face belongs to any condition
            for (auto& iter : dirich_vec)
            {
              bool faceRelevant = true;
              dirich = iter;
              for (int j = 0; j < nummynodes; ++j)
              {
                if (!dirich->ContainsNode(mynodes[j]))
                {
                  faceRelevant = false;
                  break;
                }
              }
              // If the face is not relevant the dirichlet flag is always zero
              if (!faceRelevant)
              {
                continue;  // This is related to the dirichlet conditions loop
              }
              else
              {
                dirichRelevant = true;
                break;  // We found the right dirichlet
              }
            }

            if (!dirichRelevant)
            {
              for (int j = 0; j < this->NumDofPerFace(i); ++j) lmdirich.push_back(0);
              continue;
            }

            const std::vector<int>* flag = nullptr;
            if (dirich->Type() != DRT::Condition::PointDirichlet &&
                dirich->Type() != DRT::Condition::LineDirichlet &&
                dirich->Type() != DRT::Condition::SurfaceDirichlet &&
                dirich->Type() != DRT::Condition::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = dirich->Get<std::vector<int>>("onoff");

            // Every component gets NumDofPerComponent ones or zeros
            for (unsigned j = 0; j < flag->size(); ++j)
              for (int k = 0; k < NumDofPerComponent(i); ++k)
              {
                if (flag && (*flag)[j])
                  lmdirich.push_back(1);
                else
                  lmdirich.push_back(0);
              }
          }
        }
      }
    }

    if (doDirichlet)
    {
      const std::vector<int>* flag = nullptr;
      DRT::Condition* dirich = GetCondition("Dirichlet");
      if (dirich)
      {
        if (dirich->Type() != DRT::Condition::PointDirichlet &&
            dirich->Type() != DRT::Condition::LineDirichlet &&
            dirich->Type() != DRT::Condition::SurfaceDirichlet &&
            dirich->Type() != DRT::Condition::VolumeDirichlet)
          FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        flag = dirich->Get<std::vector<int>>("onoff");
      }
      for (unsigned j = 0; j < dof.size(); ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
      }
    }

  }  // for (int dofset=0; dofset<la.Size(); ++dofset)

  return;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
    const std::string& condstring, Teuchos::ParameterList& params) const
{
  /* This method is intended to fill the LocationArray with the dofs
   * the element will assemble into. In the standard case implemented here
   * these dofs are the dofs of the element itself. For some special conditions (e.g.
   * the weak dirichlet boundary condtion) a surface element will assemble
   * into the dofs of a volume element. These elements need to overwrite this
   * method.
   */
  LocationVector(dis, la, doDirichlet);
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const Discretization& dis, std::vector<int>& lm,
    std::vector<int>& lmdirich, std::vector<int>& lmowner, std::vector<int>& lmstride) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  lm.clear();
  lmdirich.clear();
  lmowner.clear();
  lmstride.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i = 0; i < numnode; ++i)
    {
      DRT::Condition* dirich = nodes[i]->GetCondition("Dirichlet");
      const std::vector<int>* flag = nullptr;
      if (dirich)
      {
        if (dirich->Type() != DRT::Condition::PointDirichlet &&
            dirich->Type() != DRT::Condition::LineDirichlet &&
            dirich->Type() != DRT::Condition::SurfaceDirichlet &&
            dirich->Type() != DRT::Condition::VolumeDirichlet)
          FOUR_C_THROW("condition with name dirichlet is not of type Dirichlet");
        flag = dirich->Get<std::vector<int>>("onoff");
      }
      const int owner = nodes[i]->Owner();
      std::vector<int> dof;
      dis.Dof(dof, nodes[i], 0, 0);
      const int size = dof.size();
      lmstride.push_back(size);
      for (int j = 0; j < size; ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
        lmowner.push_back(owner);
        lm.push_back(dof[j]);
      }
    }
  }

  // fill the vectors with element dofs
  unsigned bef = lm.size();
  dis.Dof(0, this, lm);
  unsigned aft = lm.size();
  if (aft - bef) lmstride.push_back((int)(aft - bef));
  lmowner.resize(lm.size(), Owner());

  // fill the vector with face dofs
  if (this->NumDofPerFace(0) > 0)
  {
    for (int i = 0; i < NumFace(); ++i)
    {
      const int owner = face_[i]->Owner();
      std::vector<int> dof = dis.Dof(0, face_[i].getRawPtr());
      if (!dof.empty()) lmstride.push_back(dof.size());
      for (int j : dof)
      {
        lmowner.push_back(owner);
        lm.push_back(j);
      }
    }
  }

  // do dirichlet BCs
  const std::vector<int>* flag = nullptr;
  DRT::Condition* dirich = GetCondition("Dirichlet");
  if (dirich)
  {
    if (dirich->Type() != DRT::Condition::PointDirichlet &&
        dirich->Type() != DRT::Condition::LineDirichlet &&
        dirich->Type() != DRT::Condition::SurfaceDirichlet &&
        dirich->Type() != DRT::Condition::VolumeDirichlet)
      FOUR_C_THROW("condition with name dirichlet is not of type Dirichlet");
    flag = dirich->Get<std::vector<int>>("onoff");
  }
  const int owner = Owner();
  std::vector<int> dof = dis.Dof(this);
  if (!dof.empty()) lmstride.push_back((int)dof.size());
  for (unsigned j = 0; j < dof.size(); ++j)
  {
    if (flag && (*flag)[j])
      lmdirich.push_back(1);
    else
      lmdirich.push_back(0);
    lmowner.push_back(owner);
    lm.push_back(dof[j]);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(const Discretization& dis, std::vector<int>& lm,
    std::vector<int>& lmowner, std::vector<int>& lmstride) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  lm.clear();
  lmowner.clear();
  lmstride.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i = 0; i < numnode; ++i)
    {
      const Node* node = nodes[i];
      unsigned bef = lm.size();
      dis.Dof(0, this, node, lm);
      unsigned aft = lm.size();
      if (aft - bef) lmstride.push_back((int)(aft - bef));
      lmowner.resize(lm.size(), node->Owner());
    }
  }

  // fill the vector with element dofs
  unsigned bef = lm.size();
  dis.Dof(0, this, lm);
  unsigned aft = lm.size();
  if (aft - bef) lmstride.push_back((int)(aft - bef));
  lmowner.resize(lm.size(), Owner());

  // fill the vector with face dofs
  if (NumDofPerFace(0) > 0)
  {
    for (int i = 0; i < NumFace(); ++i)
    {
      const int owner = face_[i]->Owner();
      std::vector<int> dof = dis.Dof(0, face_[i].getRawPtr());
      if (!dof.empty()) lmstride.push_back(dof.size());
      for (int j : dof)
      {
        lmowner.push_back(owner);
        lm.push_back(j);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  return number of faces (public)                    kronbichler 05/13|
 *----------------------------------------------------------------------*/
int DRT::Element::NumFace() const
{
  switch (CORE::FE::getDimension(this->Shape()))
  {
    case 2:
      return NumLine();
    case 3:
      return NumSurface();
    default:
      FOUR_C_THROW("faces for discretization type %s not yet implemented",
          (CORE::FE::CellTypeToString(Shape())).c_str());
      return 0;
  }
}

/*----------------------------------------------------------------------*
 |  returns neighbor of element (public)               kronbichler 05/13|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Element::Neighbor(const int face) const
{
  if (face_.empty()) return nullptr;
  FOUR_C_ASSERT(face < NumFace(), "there is no face with the given index");
  DRT::FaceElement* faceelement = face_[face].getRawPtr();
  if (faceelement->ParentMasterElement() == this)
    return faceelement->ParentSlaveElement();
  else if (faceelement->ParentSlaveElement() == this)
    return faceelement->ParentMasterElement();
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  set faces (public)                                 kronbichler 05/13|
 *----------------------------------------------------------------------*/
void DRT::Element::SetFace(const int faceindex, DRT::FaceElement* faceelement)
{
  const int nface = NumFace();
  if (face_.empty()) face_.resize(nface, Teuchos::null);
  FOUR_C_ASSERT(faceindex < NumFace(), "there is no face with the given index");
  face_[faceindex] = Teuchos::rcpFromRef<DRT::FaceElement>(*faceelement);
}

/*----------------------------------------------------------------------*
 |  set faces (public)                                       seitz 12/16|
 *----------------------------------------------------------------------*/
void DRT::Element::SetFace(const int faceindex, Teuchos::RCP<DRT::FaceElement> faceelement)
{
  const int nface = NumFace();
  if (face_.empty()) face_.resize(nface, Teuchos::null);
  FOUR_C_ASSERT(faceindex < NumFace(), "there is no face with the given index");
  face_[faceindex] = std::move(faceelement);
}


/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Element::Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
    LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3)
{
  return Evaluate(params, discretization, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
}

/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Element::Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
    std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3)
{
  std::cout << "DRT::Element::Evaluate:\n"
            << "Base class dummy routine DRT::Element::Evaluate(...) called\n"
            << __FILE__ << ":" << __LINE__ << std::endl;
  return -1;
}

int DRT::Element::Degree() const { return CORE::FE::getDegree(Shape()); }

/*----------------------------------------------------------------------*
 |  check if the element has only ghost nodes (public)       vuong 09/14|
 *----------------------------------------------------------------------*/
bool DRT::Element::HasOnlyGhostNodes(const int mypid) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  // check for 'purely ghosted' element, i.e. only ghost nodes
  bool allghostnodes = true;
  for (int i = 0; i < numnode; ++i)
  {
    if (nodes[i]->Owner() == mypid)
    {
      // one node is not ghosted ->leave
      allghostnodes = false;
      break;
    }
  }
  return allghostnodes;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
unsigned int DRT::Element::AppendVisualizationGeometry(const DRT::Discretization& discret,
    std::vector<uint8_t>& cell_types, std::vector<double>& point_coordinates) const
{
  if (IsNurbsElement())
    return DRT::ELEMENTS::AppendVisualizationGeometryNURBSEle(
        *this, discret, cell_types, point_coordinates);
  else
    return DRT::ELEMENTS::AppendVisualizationGeometryLagrangeEle(
        *this, discret, cell_types, point_coordinates);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
unsigned int DRT::Element::AppendVisualizationDofBasedResultDataVector(
    const DRT::Discretization& discret, const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
    unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
    std::vector<double>& vtu_point_result_data) const
{
  if (IsNurbsElement())
    return DRT::ELEMENTS::AppendVisualizationDofBasedResultDataVectorNURBSEle(*this, discret,
        result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
        vtu_point_result_data);
  else
    return DRT::ELEMENTS::AppendVisualizationDofBasedResultDataVectorLagrangeEle(*this, discret,
        result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
        vtu_point_result_data);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::GEOMETRICSEARCH::BoundingVolume DRT::Element::GetBoundingVolume(
    const DRT::Discretization& discret, const Epetra_Vector& result_data_dofbased,
    const CORE::GEOMETRICSEARCH::GeometricSearchParams& params) const
{
  CORE::GEOMETRICSEARCH::BoundingVolume bounding_box;
  CORE::LINALG::Matrix<3, 1, double> point;

  // The default bounding box is simply the bounding box of all element nodes.
  for (unsigned int i_node = 0; i_node < (unsigned int)this->NumNode(); ++i_node)
  {
    std::vector<int> nodedofs;
    nodedofs.clear();

    // local storage position of desired dof gid
    const DRT::Node* node = this->Nodes()[i_node];
    discret.Dof(node, nodedofs);

    for (unsigned int i_dir = 0; i_dir < 3; ++i_dir)
    {
      const int lid = result_data_dofbased.Map().LID(nodedofs[i_dir]);

      if (lid > -1)
        point(i_dir) = node->X()[i_dir] + result_data_dofbased[lid];
      else
        FOUR_C_THROW("received illegal dof local id: %d", lid);
    }
    bounding_box.AddPoint(point);
  }

  return bounding_box;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                               kronbichler 03/15|
 *----------------------------------------------------------------------*/
DRT::FaceElement::FaceElement(const int id, const int owner)
    : DRT::Element(id, owner),
      parent_master_(nullptr),
      parent_slave_(nullptr),
      lface_master_(-1),
      lface_slave_(-1),
      parent_id_(-1)
{
}



/*----------------------------------------------------------------------*
 |  Copy constructor (public)                          kronbichler 03/15|
 *----------------------------------------------------------------------*/
DRT::FaceElement::FaceElement(const DRT::FaceElement& old)
    : DRT::Element(old),
      parent_master_(old.parent_master_),
      parent_slave_(old.parent_slave_),
      lface_master_(old.lface_master_),
      lface_slave_(old.lface_slave_),
      localtrafomap_(old.localtrafomap_),
      parent_id_(old.parent_id_)
{
}
/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           ager 06/15 |
 *----------------------------------------------------------------------*/
void DRT::FaceElement::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class DRT::Elememt
  DRT::Element::Pack(data);
  // add lface_master_
  AddtoPack(data, lface_master_);
  // Pack Parent Id, used to set parent_master_ after parallel communication!
  AddtoPack(data, parent_id_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           ager 06/15 |
 *----------------------------------------------------------------------*/
void DRT::FaceElement::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class DRT::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  // lface_master_
  lface_master_ = ExtractInt(position, data);
  // Parent Id
  parent_id_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  set the local trafo map (protected)                kronbichler 03/15|
 *----------------------------------------------------------------------*/
void DRT::FaceElement::SetLocalTrafoMap(const std::vector<int>& trafo) { localtrafomap_ = trafo; }

FOUR_C_NAMESPACE_CLOSE
