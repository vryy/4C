/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 3


*----------------------------------------------------------------------*/

#include "4C_so3_nstet5.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::NStet5Type Discret::ELEMENTS::NStet5Type::instance_;

Discret::ELEMENTS::NStet5Type& Discret::ELEMENTS::NStet5Type::Instance() { return instance_; }


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Core::Communication::ParObject* Discret::ELEMENTS::NStet5Type::Create(const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::NStet5(-1, -1);
  object->Unpack(data);
  return object;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::NStet5Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::NStet5(id, owner));
    return ele;
  }
  return Teuchos::null;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::NStet5Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::NStet5(id, owner));
  return ele;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void Discret::ELEMENTS::NStet5Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::NStet5Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  // TODO: switch to correct data container!
  // do nullspace for element degrees of freedom
  /*
  const Epetra_Map* rowmap = dis.dof_row_map(0);
  const int lrows = rowmap->NumMyElements();

   double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowElements(); ++i)
  {
    Core::Elements::Element* ele = dis.lRowElement(i);
    auto* nstet = dynamic_cast<Discret::ELEMENTS::NStet5*>(ele);
    if (!nstet) continue;
    const double* x = nstet->MidX();
    std::vector<int> dofs = dis.Dof(0, ele);
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (dofs.size() != 3) FOUR_C_THROW("Wrong number of dofs");
#endif
    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) FOUR_C_THROW("Cannot find element dof in dofrowmap");
      switch (j)
      {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
          break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
          break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
          break;
        default:
          FOUR_C_THROW("Only dofs 0 - 5 supported");
          break;
      }  // switch (j)
    }
  }
  */

  if (numdof != 3)
    FOUR_C_THROW(
        "The computation of the solid nullspace in three dimensions requires three DOFs"
        "per solid node, however the current node carries %d DOFs.",
        numdof);

  if (dimnsp != 6)
    FOUR_C_THROW(
        "The computation of the solid nullspace in three dimensions requires six nullspace"
        "vectors per node, however the current node carries %d vectors.",
        dimnsp);

  Discret::ELEMENTS::NStet5* nstet = dynamic_cast<Discret::ELEMENTS::NStet5*>(node.Elements()[0]);
  if (!nstet) FOUR_C_THROW("Cannot cast to NStet5");
  const double* x = nstet->mid_x();

  Core::LinAlg::SerialDenseMatrix nullspace(numdof, dimnsp);
  // x-modes
  nullspace(0, 0) = 1.0;
  nullspace(0, 1) = 0.0;
  nullspace(0, 2) = 0.0;
  nullspace(0, 3) = 0.0;
  nullspace(0, 4) = x[2] - x0[2];
  nullspace(0, 5) = -x[1] + x0[1];
  // y-modes
  nullspace(1, 0) = 0.0;
  nullspace(1, 1) = 1.0;
  nullspace(1, 2) = 0.0;
  nullspace(1, 3) = -x[2] + x0[2];
  nullspace(1, 4) = 0.0;
  nullspace(1, 5) = x[0] - x0[0];
  // z-modes
  nullspace(2, 0) = 0.0;
  nullspace(2, 1) = 0.0;
  nullspace(2, 2) = 1.0;
  nullspace(2, 3) = x[1] - x0[1];
  nullspace(2, 4) = -x[0] + x0[0];
  nullspace(2, 5) = 0.0;

  return nullspace;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void Discret::ELEMENTS::NStet5Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = Input::LineDefinition::Builder()
                     .AddIntVector("TET4", 4)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("GROWTHTRIG")
                     .Build();
}


/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 03/12|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::NStet5::NStet5(int id, int owner)
    : Core::Elements::Element(id, owner),
      material_(0),
      V_(-1.0),
      pstype_(Inpar::STR::PreStress::none),
      pstime_(0.0),
      time_(0.0)
{
  sublm_[0] = 0;
  sublm_[1] = 1;
  sublm_[2] = 2;
  sublm_[3] = 4;
  sublm_[4] = 1;
  sublm_[5] = 3;
  sublm_[6] = 2;
  sublm_[7] = 4;
  sublm_[8] = 0;
  sublm_[9] = 3;
  sublm_[10] = 1;
  sublm_[11] = 4;
  sublm_[12] = 0;
  sublm_[13] = 2;
  sublm_[14] = 3;
  sublm_[15] = 4;

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = Prestress::GetType();
    pstime_ = Prestress::GetPrestressTime();

    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }
  if (Prestress::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(4, 4, true));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 03/128|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::NStet5::NStet5(const Discret::ELEMENTS::NStet5& old)
    : Core::Elements::Element(old),
      material_(old.material_),
      V_(old.V_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  for (int i = 0; i < 16; ++i) sublm_[i] = old.sublm_[i];

  if (Prestress::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(*(old.prestress_)));
}



/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);
  // material_
  add_to_pack(data, material_);
  // stresstype_
  add_to_pack(data, stresstype_);
  // V_
  add_to_pack(data, V_);

  // Pack prestress
  add_to_pack(data, static_cast<int>(pstype_));
  add_to_pack(data, pstime_);
  add_to_pack(data, time_);
  if (Prestress::IsMulf(pstype_))
  {
    Core::Communication::ParObject::add_to_pack(data, *prestress_);
  }
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);
  // material_
  extract_from_pack(position, data, material_);
  // stresstype_
  stresstype_ = static_cast<StressType>(ExtractInt(position, data));
  // V_
  extract_from_pack(position, data, V_);

  // Extract prestress
  pstype_ = static_cast<Inpar::STR::PreStress>(ExtractInt(position, data));
  extract_from_pack(position, data, pstime_);
  extract_from_pack(position, data, time_);
  if (Prestress::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    extract_from_pack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(4, 4, true));
    prestress_->Unpack(tmpprestress);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5::so_nstet5_expol(
    Core::LinAlg::Matrix<1, 6>& stresses, Core::LinAlg::Matrix<4, 6>& nodalstresses)
{
  Core::LinAlg::Matrix<4, 1> expol;
  expol(0, 0) = 1.0;
  expol(1, 0) = 1.0;
  expol(2, 0) = 1.0;
  expol(3, 0) = 1.0;
  nodalstresses.Multiply(expol, stresses);
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 03/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5::Print(std::ostream& os) const
{
  os << "NStet5 ";
  Element::Print(os);
  return;
}


/*====================================================================*/
/* 4-node tetrahedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
 * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
 *  4-node hexahedron: node 0,1,...,3
 *
 * -----------------------
 *- this is the numbering used in GiD & EXODUS!!
 *      3-
 *      |\ ---
 *      |  \    ---
 *      |    \      ---
 *      |      \        -2
 *      |        \       /\
 *      |          \   /   \
 *      |            X      \
 *      |          /   \     \
 *      |        /       \    \
 *      |      /           \   \
 *      |    /               \  \
 *      |  /                   \ \
 *      |/                       \\
 *      0--------------------------1
 */
/*====================================================================*/


/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             gee 03/12|
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::NStet5::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 03/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::NStet5::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5Type::init_elementsand_maps(
    std::map<int, Discret::ELEMENTS::NStet5*>& elecids, std::map<int, Core::Nodes::Node*>& noderids,
    const int myrank, const int numproc, Core::FE::Discretization& dis)
{
  const int numele = dis.NumMyColElements();

  for (int i = 0; i < numele; ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::NStet5*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to NStet5* failed");

    // init the element
    actele->init_element();

    // register element in list of column nstet elements
    elecids[actele->Id()] = actele;

    // compute a map of all row nodes adjacent to a NStet5 element
    for (int j = 0; j < actele->num_node(); ++j)
    {
      Core::Nodes::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner()) noderids[node->Id()] = node;
    }
  }  // i

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5Type::init_adjacency(
    std::map<int, Discret::ELEMENTS::NStet5*>& elecids, std::map<int, Core::Nodes::Node*>& noderids,
    std::map<int, std::vector<Discret::ELEMENTS::NStet5*>>& adjele,
    std::map<int, std::map<int, Core::Nodes::Node*>>& adjnode,
    std::map<int, std::vector<int>>& adjlm,
    std::map<int, std::map<int, std::vector<int>>>& adjsubele,
    std::map<int, std::vector<std::vector<std::vector<int>>>>& adjlmlm,
    Core::FE::Discretization& dis)
{
  std::map<int, Core::Nodes::Node*>::iterator node;
  for (node = noderids.begin(); node != noderids.end(); ++node)
  {
    Core::Nodes::Node* nodeL = node->second;
    const int nodeidL = nodeL->Id();

    //-----------------------------------------------------------------
    // list of adjacent elements
    std::vector<Discret::ELEMENTS::NStet5*> myadjele(0);
    for (int j = 0; j < nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      auto ele = elecids_.find(eleid);
      if (ele == elecids_.end()) continue;
      myadjele.push_back(ele->second);
    }
    adjele[nodeidL] = myadjele;

    //-----------------------------------------------------------------
    // patch of all nodes adjacent to adjacent elements
    std::map<int, Core::Nodes::Node*> nodepatch;
    for (auto& j : myadjele)
    {
      Core::Nodes::Node** nodes = j->Nodes();
      for (int k = 0; k < j->num_node(); ++k) nodepatch[nodes[k]->Id()] = nodes[k];
    }
    adjnode[nodeidL] = nodepatch;

    //-----------------------------------------------------------------
    // lm array
    const int ndofperpatch = ((int)nodepatch.size() + (int)myadjele.size()) * 3;

    // location and ownership vector of nodal patch
    std::vector<int> lm(ndofperpatch);
    std::map<int, Core::Nodes::Node*>::iterator pnode;
    int count = 0;
    // add dofs of nodes
    for (pnode = nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      const std::vector<int>& dofs = dis.Dof(pnode->second);
      for (int dof : dofs) lm[count++] = dof;
    }

    // add dofs of center nodes from elements. These appear as element dofs
    for (auto& j : myadjele)
    {
      const std::vector<int>& dofs = dis.Dof(j);
      for (int dof : dofs) lm[count++] = dof;
    }

    adjlm[nodeidL] = lm;

    //-----------------------------------------------------------------
    // for each adjele, find out which subelements I participate in
    std::map<int, std::vector<int>> masterele;
    for (auto ele : myadjele)
    {
      bool foundit = false;
      for (int i = 0; i < ele->num_node(); ++i)
      {
        if (ele->Nodes()[i]->Id() == nodeL->Id())
        {
          // found the center node on the element
          // local to the element, its node i
          foundit = true;
          // determine subelements node i is attached to
          // its attached to definitely 3 out of 4 subelements
          std::vector<int> subele;
          for (int k = 0; k < 4; ++k)
          {
            const int* sublm = ele->sub_lm(k);  // subelement k
            for (int l = 0; l < 3; ++l)         // the first 3 nodes of the subelement
              if (sublm[l] == i)
              {
                subele.push_back(k);
                break;
              }
          }
          if ((int)subele.size() != 3) FOUR_C_THROW("Node not attached to exactly 3 subelements");

          masterele[ele->Id()] = subele;

          // no longer need to look at this element
          break;
        }
      }
      if (!foundit) FOUR_C_THROW("Weired, this adjele seems not attached to me");
    }  // for (unsigned j=0; j<myadjele.size(); ++j)
    if (masterele.size() != myadjele.size()) FOUR_C_THROW("subelement connectivity wrong");

    adjsubele[nodeidL] = masterele;

    //-----------------------------------------------------------------
    // for each adjele and its subele, build local connectivity
    std::vector<std::vector<std::vector<int>>> lmlm((int)myadjele.size());
    for (unsigned j = 0; j < myadjele.size(); ++j)
    {
      Discret::ELEMENTS::NStet5* ele = myadjele[j];
      std::vector<int>& subele = masterele[ele->Id()];
      lmlm[j].resize((int)subele.size());
      for (unsigned k = 0; k < subele.size(); ++k)
      {
        const int subeleid = subele[k];
        const int* sublm = ele->sub_lm(subeleid);
        std::vector<int> elelm;
        for (int l = 0; l < 4; ++l)  // loop nodes of subelement and collect dofs
        {
          if (sublm[l] != 4)  // node 4 is center node owned by the element
          {
            std::vector<int> dofs = dis.Dof(ele->Nodes()[sublm[l]]);
            for (int dof : dofs) elelm.push_back(dof);
          }
          else
          {
            std::vector<int> dofs = dis.Dof(ele);
            for (int dof : dofs) elelm.push_back(dof);
          }
        }
        if ((int)elelm.size() != 12) FOUR_C_THROW("Subelement does not have 12 dofs");
        lmlm[j][k].resize(12);
        for (int l = 0; l < 12; ++l)
        {
          auto fool = find(lm.begin(), lm.end(), elelm[l]);
          lmlm[j][k][l] = fool - lm.begin();
        }
      }
    }
    adjlmlm[nodeidL] = lmlm;



  }  // for (node=noderids.begin(); node != noderids.end(); ++node)
  return;
}



/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 03/12|
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::NStet5Type::Initialize(Core::FE::Discretization& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("Discret::ELEMENTS::NStet5Type::Initialize");

  const int myrank = dis.Comm().MyPID();
  const int numproc = dis.Comm().NumProc();

  //----------------------------------------------------------------------
  // init elements, make maps of column elements and row nodes
  init_elementsand_maps(elecids_, noderids_, myrank, numproc, dis);

  //----------------------------------------------------------------------
  // compute adjacency for each row node
  // make patch of adjacent elements
  // make patch of adjacent nodes (including center node itself)
  // make lm for nodal patch
  init_adjacency(elecids_, noderids_, adjele_, adjnode_, adjlm_, adjsubele_, lmlm_, dis);


  return 0;
}

FOUR_C_NAMESPACE_CLOSE
