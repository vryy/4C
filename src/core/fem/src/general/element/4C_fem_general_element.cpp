/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of base element class in 4C with basic operations

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_fem_general_element.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_geometric_search_bounding_volume.hpp"
#include "4C_fem_geometric_search_params.hpp"
#include "4C_io_element_append_visualization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_material_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Shards_BasicTopologies.hpp>

#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::CellType Core::Elements::ShardsKeyToDisType(const unsigned& key)
{
  Core::FE::CellType distype = Core::FE::CellType::dis_none;
  switch (key)
  {
    case shards::Particle::key:
    {
      distype = Core::FE::CellType::point1;
      break;
    }
    case shards::Line<2>::key:
    {
      distype = Core::FE::CellType::line2;
      break;
    }
    case shards::Line<3>::key:
    {
      distype = Core::FE::CellType::line3;
      break;
    }
    case shards::Quadrilateral<4>::key:
    {
      distype = Core::FE::CellType::quad4;
      break;
    }
    case shards::Quadrilateral<8>::key:
    {
      distype = Core::FE::CellType::quad8;
      break;
    }
    case shards::Quadrilateral<9>::key:
    {
      distype = Core::FE::CellType::quad9;
      break;
    }
    case shards::Triangle<3>::key:
    {
      distype = Core::FE::CellType::tri3;
      break;
    }
    case shards::Triangle<6>::key:
    {
      distype = Core::FE::CellType::tri6;
      break;
    }
    case shards::Hexahedron<8>::key:
    {
      distype = Core::FE::CellType::hex8;
      break;
    }
    case shards::Hexahedron<20>::key:
    {
      distype = Core::FE::CellType::hex20;
      break;
    }
    case shards::Hexahedron<27>::key:
    {
      distype = Core::FE::CellType::hex27;
      break;
    }
    case shards::Tetrahedron<4>::key:
    {
      distype = Core::FE::CellType::tet4;
      break;
    }
    case shards::Tetrahedron<10>::key:
    {
      distype = Core::FE::CellType::tet10;
      break;
    }
    case shards::Wedge<6>::key:
    {
      distype = Core::FE::CellType::wedge6;
      break;
    }
    case shards::Wedge<15>::key:
    {
      distype = Core::FE::CellType::wedge15;
      break;
    }
    case shards::Pyramid<5>::key:
    {
      distype = Core::FE::CellType::pyramid5;
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
Core::Elements::Element::Element(int id, int owner)
    : ParObject(), id_(id), lid_(-1), owner_(owner), mat_(1, Teuchos::null), is_nurbs_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Core::Elements::Element::Element(const Element& old)
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
  std::map<std::string, Teuchos::RCP<Core::Conditions::Condition>>::const_iterator fool;
  for (fool = old.condition_.begin(); fool != old.condition_.end(); ++fool)
    set_condition(fool->first, fool->second);

  if (!old.mat_.empty())
  {
    mat_.resize(old.mat_.size());
    for (unsigned iter = 0; iter < old.mat_.size(); ++iter)
      if (old.mat_[iter] != Teuchos::null) mat_[iter] = (old.mat_[iter]->clone());
  }
  else
    mat_[0] = Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Core::Elements::Element& element)
{
  element.print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void Core::Elements::Element::print(std::ostream& os) const
{
  os << std::setw(12) << id() << " Owner " << std::setw(5) << owner() << " ";
  const int nnode = num_node();
  const int* nodeids = node_ids();
  if (nnode > 0)
  {
    os << " Nodes ";
    for (int i = 0; i < nnode; ++i) os << std::setw(10) << nodeids[i] << " ";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Core::Elements::Element::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  FOUR_C_THROW("subclass implementations missing");
  return false;
}


/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void Core::Elements::Element::set_node_ids(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i = 0; i < nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::Element::set_node_ids(
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  nodeid_ = container.get<std::vector<int>>(distype);
  for (int& i : nodeid_) i -= 1;
  node_.resize(0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Elements::Element::set_material(const int index, Teuchos::RCP<Core::Mat::Material> mat)
{
  FOUR_C_THROW_UNLESS(mat != Teuchos::null,
      "Invalid material given to the element. \n"
      "Invalid are Summands of the Elasthyper-Toolbox and single Growth-Materials. \n"
      "If you like to use a Summand of the Elasthyper-Material define it via MAT_ElastHyper. \n"
      "If you like to use a Growth-Material define it via the according base material.");

  if (num_material() > index)
    mat_[index] = mat;
  else if (num_material() == index)
    add_material(mat);
  else
    FOUR_C_THROW(
        "Setting material at index %d not possible (neither overwrite nor append) since currently  "
        "only %d materials are stored",
        index, num_material());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::Elements::Element::add_material(Teuchos::RCP<Core::Mat::Material> mat)
{
  mat_.push_back(mat);

  return mat_.size();
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add id
  add_to_pack(data, id_);
  // add owner
  add_to_pack(data, owner_);
  // add vector nodeid_
  add_to_pack(data, nodeid_);
  // add material
  if (mat_[0] != Teuchos::null)
  {
    // pack only first material
    mat_[0]->pack(data);
  }
  else
  {
    int size = 0;
    add_to_pack(data, size);
  }
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // id_
  extract_from_pack(position, data, id_);
  // owner_
  extract_from_pack(position, data, owner_);
  // nodeid_
  extract_from_pack(position, data, nodeid_);
  // mat_
  std::vector<char> tmp;
  extract_from_pack(position, data, tmp);
  if (!tmp.empty())
  {
    Core::Communication::ParObject* o = Core::Communication::Factory(tmp);
    auto* mat = dynamic_cast<Core::Mat::Material*>(o);
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
    std::vector<Teuchos::RCP<Core::Elements::FaceElement>> empty;
    std::swap(face_, empty);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool Core::Elements::Element::build_nodal_pointers(
    std::map<int, Teuchos::RCP<Core::Nodes::Node>>& nodes)
{
  int nnode = num_node();
  const int* nodeids = node_ids();
  node_.resize(nnode);
  for (int i = 0; i < nnode; ++i)
  {
    std::map<int, Teuchos::RCP<Core::Nodes::Node>>::const_iterator curr = nodes.find(nodeids[i]);
    // this node is not on this proc
    if (curr == nodes.end())
      FOUR_C_THROW("Element %d cannot find node %d", id(), nodeids[i]);
    else
      node_[i] = curr->second.get();
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool Core::Elements::Element::build_nodal_pointers(Core::Nodes::Node** nodes)
{
  node_.resize(num_node());
  for (int i = 0; i < num_node(); ++i) node_[i] = nodes[i];
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal connectivity and weight nodes and edges        (public) |
 |                                                          ghamm 09/13 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::nodal_connectivity(
    Core::LinAlg::SerialDenseMatrix& edgeweights, Core::LinAlg::SerialDenseVector& nodeweights)
{
  // weight for this element
  double weight = evaluation_cost();

  int numnode = num_node();
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

  std::vector<std::vector<int>> lines = Core::FE::getEleNodeNumberingLines(shape());
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
        Core::FE::CellTypeToString(shape()).c_str());
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::get_condition(
    const std::string& name, std::vector<Core::Conditions::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  auto startit = condition_.lower_bound(name);
  auto endit = condition_.upper_bound(name);
  int count = 0;
  std::multimap<std::string, Teuchos::RCP<Core::Conditions::Condition>>::const_iterator curr;
  for (curr = startit; curr != endit; ++curr) out[count++] = curr->second.get();
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
Core::Conditions::Condition* Core::Elements::Element::get_condition(const std::string& name) const
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
void Core::Elements::Element::location_vector(const Core::FE::Discretization& dis,
    const std::vector<int>& nds, Core::Elements::Element::LocationArray& la, bool doDirichlet) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Element::nodes();

  if (numnode != static_cast<int>(nds.size()))
  {
    FOUR_C_THROW("wrong number of nodes");
  }

  la.clear();

  // we need to look at all DofSets of our discretization
  for (int dofset = 0; dofset < la.size(); ++dofset)
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
        const Core::Nodes::Node* node = nodes[i];

        const int owner = node->owner();
        std::vector<int> dof;
        dis.dof(dof, node, dofset, nds[i]);
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
          Core::Conditions::Condition* dirich = node->get_condition("Dirichlet");
          if (dirich)
          {
            if (dirich->type() != Core::Conditions::PointDirichlet &&
                dirich->type() != Core::Conditions::LineDirichlet &&
                dirich->type() != Core::Conditions::SurfaceDirichlet &&
                dirich->type() != Core::Conditions::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = &dirich->parameters().get<std::vector<int>>("onoff");
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
    const int owner = Element::owner();
    std::vector<int> dof = dis.dof(dofset, this);
    if (!dof.empty()) lmstride.push_back(dof.size());
    for (int j : dof)
    {
      lmowner.push_back(owner);
      lm.push_back(j);
    }

    // fill the vector with face dofs
    if (this->num_dof_per_face(0) > 0)
    {
      for (int i = 0; i < num_face(); ++i)
      {
        const int owner = face_[i]->owner();
        std::vector<int> dof = dis.dof(dofset, face_[i].getRawPtr());
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
      Core::Conditions::Condition* dirich = get_condition("Dirichlet");
      if (dirich)
      {
        if (dirich->type() != Core::Conditions::PointDirichlet &&
            dirich->type() != Core::Conditions::LineDirichlet &&
            dirich->type() != Core::Conditions::SurfaceDirichlet &&
            dirich->type() != Core::Conditions::VolumeDirichlet)
          FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        flag = &dirich->parameters().get<std::vector<int>>("onoff");
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
void Core::Elements::Element::location_vector(
    const Core::FE::Discretization& dis, LocationArray& la, bool doDirichlet) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Element::nodes();

  la.clear();

  // we need to look at all DofSets of our discretization
  for (int dofset = 0; dofset < la.size(); ++dofset)
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
        const Core::Nodes::Node* node = nodes[i];

        const int owner = node->owner();
        std::vector<int> dof;
        dis.dof(dof, node, dofset, 0, this);

        // if there are more dofs on the node than the element can handle, this cannot work
        FOUR_C_ASSERT(num_dof_per_node(*node) <= (int)dof.size() or dofset != 0,
            "More dofs on node than element can handle! Internal error!");

        // assume that the first dofs are the relevant ones
        const int size = dofset == 0 ? num_dof_per_node(*node) : dof.size();

        if (size) lmstride.push_back(size);
        for (int j = 0; j < size; ++j)
        {
          lmowner.push_back(owner);
          lm.push_back(dof[j]);
        }

        if (doDirichlet)
        {
          const std::vector<int>* flag = nullptr;
          Core::Conditions::Condition* dirich = node->get_condition("Dirichlet");
          if (dirich)
          {
            if (dirich->type() != Core::Conditions::PointDirichlet &&
                dirich->type() != Core::Conditions::LineDirichlet &&
                dirich->type() != Core::Conditions::SurfaceDirichlet &&
                dirich->type() != Core::Conditions::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = &dirich->parameters().get<std::vector<int>>("onoff");
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
    const int owner = Element::owner();
    std::vector<int> dof = dis.dof(dofset, this);
    if (!dof.empty()) lmstride.push_back(dof.size());
    for (int j : dof)
    {
      lmowner.push_back(owner);
      lm.push_back(j);
    }

    // fill the vector with face dofs
    if (this->num_dof_per_face(0) > 0)
    {
      for (int i = 0; i < num_face(); ++i)
      {
        const int owner = face_[i]->owner();
        std::vector<int> dof = dis.dof(dofset, face_[i].getRawPtr());
        if (!dof.empty()) lmstride.push_back(dof.size());
        for (int j : dof)
        {
          lmowner.push_back(owner);
          lm.push_back(j);
        }

        if (doDirichlet)
        {
          std::vector<Core::Conditions::Condition*> dirich_vec;
          dis.get_condition("Dirichlet", dirich_vec);
          Core::Conditions::Condition* dirich;
          bool dirichRelevant = false;
          // Check if there exist a dirichlet condition
          if (!dirich_vec.empty())
          {
            // do only faces where all nodes are present in the node list
            const int nummynodes = face_[i]->num_node();
            const int* mynodes = face_[i]->node_ids();
            // Check if the face belongs to any condition
            for (auto& iter : dirich_vec)
            {
              bool faceRelevant = true;
              dirich = iter;
              for (int j = 0; j < nummynodes; ++j)
              {
                if (!dirich->contains_node(mynodes[j]))
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
              for (int j = 0; j < this->num_dof_per_face(i); ++j) lmdirich.push_back(0);
              continue;
            }

            const std::vector<int>* flag = nullptr;
            if (dirich->type() != Core::Conditions::PointDirichlet &&
                dirich->type() != Core::Conditions::LineDirichlet &&
                dirich->type() != Core::Conditions::SurfaceDirichlet &&
                dirich->type() != Core::Conditions::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = &dirich->parameters().get<std::vector<int>>("onoff");

            // Every component gets NumDofPerComponent ones or zeros
            for (unsigned j = 0; j < flag->size(); ++j)
              for (int k = 0; k < num_dof_per_component(i); ++k)
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
      Core::Conditions::Condition* dirich = get_condition("Dirichlet");
      if (dirich)
      {
        if (dirich->type() != Core::Conditions::PointDirichlet &&
            dirich->type() != Core::Conditions::LineDirichlet &&
            dirich->type() != Core::Conditions::SurfaceDirichlet &&
            dirich->type() != Core::Conditions::VolumeDirichlet)
          FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        flag = &dirich->parameters().get<std::vector<int>>("onoff");
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
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::location_vector(const Core::FE::Discretization& dis,
    LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  /* This method is intended to fill the LocationArray with the dofs
   * the element will assemble into. In the standard case implemented here
   * these dofs are the dofs of the element itself. For some special conditions (e.g.
   * the weak dirichlet boundary condtion) a surface element will assemble
   * into the dofs of a volume element. These elements need to overwrite this
   * method.
   */
  location_vector(dis, la, doDirichlet);
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::location_vector(const Core::FE::Discretization& dis,
    std::vector<int>& lm, std::vector<int>& lmdirich, std::vector<int>& lmowner,
    std::vector<int>& lmstride) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Element::nodes();

  lm.clear();
  lmdirich.clear();
  lmowner.clear();
  lmstride.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i = 0; i < numnode; ++i)
    {
      Core::Conditions::Condition* dirich = nodes[i]->get_condition("Dirichlet");
      const std::vector<int>* flag = nullptr;
      if (dirich)
      {
        if (dirich->type() != Core::Conditions::PointDirichlet &&
            dirich->type() != Core::Conditions::LineDirichlet &&
            dirich->type() != Core::Conditions::SurfaceDirichlet &&
            dirich->type() != Core::Conditions::VolumeDirichlet)
          FOUR_C_THROW("condition with name dirichlet is not of type Dirichlet");
        flag = &dirich->parameters().get<std::vector<int>>("onoff");
      }
      const int owner = nodes[i]->owner();
      std::vector<int> dof;
      dis.dof(dof, nodes[i], 0, 0);
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
  dis.dof(0, this, lm);
  unsigned aft = lm.size();
  if (aft - bef) lmstride.push_back((int)(aft - bef));
  lmowner.resize(lm.size(), owner());

  // fill the vector with face dofs
  if (this->num_dof_per_face(0) > 0)
  {
    for (int i = 0; i < num_face(); ++i)
    {
      const int owner = face_[i]->owner();
      std::vector<int> dof = dis.dof(0, face_[i].getRawPtr());
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
  Core::Conditions::Condition* dirich = get_condition("Dirichlet");
  if (dirich)
  {
    if (dirich->type() != Core::Conditions::PointDirichlet &&
        dirich->type() != Core::Conditions::LineDirichlet &&
        dirich->type() != Core::Conditions::SurfaceDirichlet &&
        dirich->type() != Core::Conditions::VolumeDirichlet)
      FOUR_C_THROW("condition with name dirichlet is not of type Dirichlet");
    flag = &dirich->parameters().get<std::vector<int>>("onoff");
  }
  const int owner = Element::owner();
  std::vector<int> dof = dis.dof(this);
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
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Core::Elements::Element::location_vector(const Core::FE::Discretization& dis,
    std::vector<int>& lm, std::vector<int>& lmowner, std::vector<int>& lmstride) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Element::nodes();

  lm.clear();
  lmowner.clear();
  lmstride.clear();

  // fill the vector with nodal dofs
  if (nodes)
  {
    for (int i = 0; i < numnode; ++i)
    {
      const Core::Nodes::Node* node = nodes[i];
      unsigned bef = lm.size();
      dis.dof(0, this, node, lm);
      unsigned aft = lm.size();
      if (aft - bef) lmstride.push_back((int)(aft - bef));
      lmowner.resize(lm.size(), node->owner());
    }
  }

  // fill the vector with element dofs
  unsigned bef = lm.size();
  dis.dof(0, this, lm);
  unsigned aft = lm.size();
  if (aft - bef) lmstride.push_back((int)(aft - bef));
  lmowner.resize(lm.size(), owner());

  // fill the vector with face dofs
  if (num_dof_per_face(0) > 0)
  {
    for (int i = 0; i < num_face(); ++i)
    {
      const int owner = face_[i]->owner();
      std::vector<int> dof = dis.dof(0, face_[i].getRawPtr());
      if (!dof.empty()) lmstride.push_back(dof.size());
      for (int j : dof)
      {
        lmowner.push_back(owner);
        lm.push_back(j);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  return number of faces (public)                    kronbichler 05/13|
 *----------------------------------------------------------------------*/
int Core::Elements::Element::num_face() const
{
  switch (Core::FE::getDimension(this->shape()))
  {
    case 2:
      return num_line();
    case 3:
      return num_surface();
    default:
      FOUR_C_THROW("faces for discretization type %s not yet implemented",
          (Core::FE::CellTypeToString(shape())).c_str());
      return 0;
  }
}

/*----------------------------------------------------------------------*
 |  returns neighbor of element (public)               kronbichler 05/13|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Core::Elements::Element::neighbor(const int face) const
{
  if (face_.empty()) return nullptr;
  FOUR_C_ASSERT(face < num_face(), "there is no face with the given index");
  Core::Elements::FaceElement* faceelement = face_[face].getRawPtr();
  if (faceelement->parent_master_element() == this)
    return faceelement->parent_slave_element();
  else if (faceelement->parent_slave_element() == this)
    return faceelement->parent_master_element();
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  set faces (public)                                 kronbichler 05/13|
 *----------------------------------------------------------------------*/
void Core::Elements::Element::set_face(
    const int faceindex, Core::Elements::FaceElement* faceelement)
{
  const int nface = num_face();
  if (face_.empty()) face_.resize(nface, Teuchos::null);
  FOUR_C_ASSERT(faceindex < num_face(), "there is no face with the given index");
  face_[faceindex] = Teuchos::rcpFromRef<Core::Elements::FaceElement>(*faceelement);
}

/*----------------------------------------------------------------------*
 |  set faces (public)                                       seitz 12/16|
 *----------------------------------------------------------------------*/
void Core::Elements::Element::set_face(
    const int faceindex, Teuchos::RCP<Core::Elements::FaceElement> faceelement)
{
  const int nface = num_face();
  if (face_.empty()) face_.resize(nface, Teuchos::null);
  FOUR_C_ASSERT(faceindex < num_face(), "there is no face with the given index");
  face_[faceindex] = std::move(faceelement);
}


/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int Core::Elements::Element::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  return evaluate(params, discretization, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
}

/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int Core::Elements::Element::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  std::cout << "Core::Elements::Element::Evaluate:\n"
            << "Base class dummy routine Core::Elements::Element::evaluate(...) called\n"
            << __FILE__ << ":" << __LINE__ << std::endl;
  return -1;
}

int Core::Elements::Element::degree() const { return Core::FE::getDegree(shape()); }

/*----------------------------------------------------------------------*
 |  check if the element has only ghost nodes (public)       vuong 09/14|
 *----------------------------------------------------------------------*/
bool Core::Elements::Element::has_only_ghost_nodes(const int mypid) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Element::nodes();

  // check for 'purely ghosted' element, i.e. only ghost nodes
  bool allghostnodes = true;
  for (int i = 0; i < numnode; ++i)
  {
    if (nodes[i]->owner() == mypid)
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
unsigned int Core::Elements::Element::append_visualization_geometry(
    const Core::FE::Discretization& discret, std::vector<uint8_t>& cell_types,
    std::vector<double>& point_coordinates) const
{
  if (is_nurbs_element())
    return IO::AppendVisualizationGeometryNURBSEle(*this, discret, cell_types, point_coordinates);
  else
    return IO::AppendVisualizationGeometryLagrangeEle(
        *this, discret, cell_types, point_coordinates);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
unsigned int Core::Elements::Element::append_visualization_dof_based_result_data_vector(
    const Core::FE::Discretization& discret,
    const Teuchos::RCP<Epetra_Vector>& result_data_dofbased, unsigned int& result_num_dofs_per_node,
    const unsigned int read_result_data_from_dofindex,
    std::vector<double>& vtu_point_result_data) const
{
  if (is_nurbs_element())
    return IO::AppendVisualizationDofBasedResultDataVectorNURBSEle(*this, discret,
        result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
        vtu_point_result_data);
  else
    return IO::AppendVisualizationDofBasedResultDataVectorLagrangeEle(*this, discret,
        result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
        vtu_point_result_data);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::GeometricSearch::BoundingVolume Core::Elements::Element::get_bounding_volume(
    const Core::FE::Discretization& discret, const Epetra_Vector& result_data_dofbased,
    const Core::GeometricSearch::GeometricSearchParams& params) const
{
  Core::GeometricSearch::BoundingVolume bounding_box;
  Core::LinAlg::Matrix<3, 1, double> point;

  // The default bounding box is simply the bounding box of all element nodes.
  for (unsigned int i_node = 0; i_node < (unsigned int)this->num_node(); ++i_node)
  {
    std::vector<int> nodedofs;
    nodedofs.clear();

    // local storage position of desired dof gid
    const Core::Nodes::Node* node = this->nodes()[i_node];
    discret.dof(node, nodedofs);

    for (unsigned int i_dir = 0; i_dir < 3; ++i_dir)
    {
      const int lid = result_data_dofbased.Map().LID(nodedofs[i_dir]);

      if (lid > -1)
        point(i_dir) = node->x()[i_dir] + result_data_dofbased[lid];
      else
        FOUR_C_THROW("received illegal dof local id: %d", lid);
    }
    bounding_box.add_point(point);
  }

  return bounding_box;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                               kronbichler 03/15|
 *----------------------------------------------------------------------*/
Core::Elements::FaceElement::FaceElement(const int id, const int owner)
    : Element(id, owner),
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
Core::Elements::FaceElement::FaceElement(const Core::Elements::FaceElement& old)
    : Element(old),
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
void Core::Elements::FaceElement::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Discret::Elememt
  Core::Elements::Element::pack(data);
  // add lface_master_
  add_to_pack(data, lface_master_);
  // Pack Parent Id, used to set parent_master_ after parallel communication!
  add_to_pack(data, parent_id_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           ager 06/15 |
 *----------------------------------------------------------------------*/
void Core::Elements::FaceElement::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::unpack(basedata);

  // lface_master_
  lface_master_ = extract_int(position, data);
  // Parent Id
  parent_id_ = extract_int(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  set the local trafo map (protected)                kronbichler 03/15|
 *----------------------------------------------------------------------*/
void Core::Elements::FaceElement::set_local_trafo_map(const std::vector<int>& trafo)
{
  localtrafomap_ = trafo;
}

FOUR_C_NAMESPACE_CLOSE
