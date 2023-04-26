/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 3


*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include "lib_discret.H"
#include "lib_utils_factory.H"
#include "lib_dserror.H"
#include "lib_linedefinition.H"
#include "lib_globalproblem.H"
#include "lib_prestress_service.H"

#include "linalg_utils_nullspace.H"

#include "so3_nstet5.H"
#include "so3_surface.H"
#include "so3_line.H"

#include "so3_prestress.H"
#include "so3_utils.H"

DRT::ELEMENTS::NStet5Type DRT::ELEMENTS::NStet5Type::instance_;

DRT::ELEMENTS::NStet5Type& DRT::ELEMENTS::NStet5Type::Instance() { return instance_; }


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
DRT::ParObject* DRT::ELEMENTS::NStet5Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::NStet5(-1, -1);
  object->Unpack(data);
  return object;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NStet5(id, owner));
    return ele;
  }
  return Teuchos::null;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NStet5(id, owner));
  return ele;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void DRT::ELEMENTS::NStet5Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::SerialDenseMatrix<int, double> DRT::ELEMENTS::NStet5Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  // TODO: switch to correct data container!
  // do nullspace for element degrees of freedom
  /*
  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();

   double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowElements(); ++i)
  {
    DRT::Element* ele = dis.lRowElement(i);
    auto* nstet = dynamic_cast<DRT::ELEMENTS::NStet5*>(ele);
    if (!nstet) continue;
    const double* x = nstet->MidX();
    std::vector<int> dofs = dis.Dof(0, ele);
#ifdef DEBUG
    if (dofs.size() != 3) dserror("Wrong number of dofs");
#endif
    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find element dof in dofrowmap");
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
          dserror("Only dofs 0 - 5 supported");
          break;
      }  // switch (j)
    }
  }
  */

  if (numdof != 3)
    dserror(
        "The computation of the solid nullspace in three dimensions requires three DOFs"
        "per solid node, however the current node carries %d DOFs.",
        numdof);

  if (dimnsp != 6)
    dserror(
        "The computation of the solid nullspace in three dimensions requires six nullspace"
        "vectors per node, however the current node carries %d vectors.",
        dimnsp);

  DRT::ELEMENTS::NStet5* nstet = dynamic_cast<DRT::ELEMENTS::NStet5*>(node.Elements()[0]);
  if (!nstet) dserror("Cannot cast to NStet5");
  const double* x = nstet->MidX();

  Teuchos::SerialDenseMatrix<int, double> nullspace(numdof, dimnsp);
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
void DRT::ELEMENTS::NStet5Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET4"]
      .AddIntVector("TET4", 4)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3)
      .AddOptionalNamedDouble("GROWTHTRIG");
}


/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet5::NStet5(int id, int owner)
    : DRT::Element(id, owner),
      material_(0),
      V_(-1.0),
      pstype_(INPAR::STR::PreStress::none),
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

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = ::UTILS::PRESTRESS::GetType();
    pstime_ = ::UTILS::PRESTRESS::GetPrestressTime();

    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        DRT::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(4, 4, true));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 03/128|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet5::NStet5(const DRT::ELEMENTS::NStet5& old)
    : DRT::Element(old),
      material_(old.material_),
      V_(old.V_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  for (int i = 0; i < 16; ++i) sublm_[i] = old.sublm_[i];

  if (::UTILS::PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet5::~NStet5() { return; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // material_
  AddtoPack(data, material_);
  // stresstype_
  AddtoPack(data, stresstype_);
  // V_
  AddtoPack(data, V_);

  // Pack prestress
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    DRT::ParObject::AddtoPack(data, *prestress_);
  }
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // material_
  ExtractfromPack(position, data, material_);
  // stresstype_
  stresstype_ = static_cast<StressType>(ExtractInt(position, data));
  // V_
  ExtractfromPack(position, data, V_);

  // Extract prestress
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(4, 4, true));
    prestress_->Unpack(tmpprestress);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::so_nstet5_expol(
    LINALG::Matrix<1, 6>& stresses, LINALG::Matrix<4, 6>& nodalstresses)
{
  LINALG::Matrix<4, 1> expol;
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
void DRT::ELEMENTS::NStet5::Print(std::ostream& os) const
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
 |  get vector of volumes (length 1) (public)                  gee 03/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet5::Volumes()
{
  dserror("volume not impl. yet");
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}


/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             gee 03/12|
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet5::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 03/12|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet5::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      DRT::UTILS::buildLines, this);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::InitElementsandMaps(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
    std::map<int, DRT::Node*>& noderids, const int myrank, const int numproc,
    DRT::Discretization& dis)
{
  const int numele = dis.NumMyColElements();

  for (int i = 0; i < numele; ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::NStet5*>(dis.lColElement(i));
    if (!actele) dserror("cast to NStet5* failed");

    // init the element
    actele->InitElement();

    // register element in list of column nstet elements
    elecids[actele->Id()] = actele;

    // compute a map of all row nodes adjacent to a NStet5 element
    for (int j = 0; j < actele->NumNode(); ++j)
    {
      DRT::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner()) noderids[node->Id()] = node;
    }
  }  // i

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5Type::InitAdjacency(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
    std::map<int, DRT::Node*>& noderids, std::map<int, std::vector<DRT::ELEMENTS::NStet5*>>& adjele,
    std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& adjlm,
    std::map<int, std::map<int, std::vector<int>>>& adjsubele,
    std::map<int, std::vector<std::vector<std::vector<int>>>>& adjlmlm, DRT::Discretization& dis)
{
  std::map<int, DRT::Node*>::iterator node;
  for (node = noderids.begin(); node != noderids.end(); ++node)
  {
    DRT::Node* nodeL = node->second;
    const int nodeidL = nodeL->Id();

    //-----------------------------------------------------------------
    // list of adjacent elements
    std::vector<DRT::ELEMENTS::NStet5*> myadjele(0);
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
    std::map<int, DRT::Node*> nodepatch;
    for (auto& j : myadjele)
    {
      DRT::Node** nodes = j->Nodes();
      for (int k = 0; k < j->NumNode(); ++k) nodepatch[nodes[k]->Id()] = nodes[k];
    }
    adjnode[nodeidL] = nodepatch;

    //-----------------------------------------------------------------
    // lm array
    const int ndofperpatch = ((int)nodepatch.size() + (int)myadjele.size()) * 3;

    // location and ownership vector of nodal patch
    std::vector<int> lm(ndofperpatch);
    std::map<int, DRT::Node*>::iterator pnode;
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
      for (int i = 0; i < ele->NumNode(); ++i)
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
            const int* sublm = ele->SubLM(k);  // subelement k
            for (int l = 0; l < 3; ++l)        // the first 3 nodes of the subelement
              if (sublm[l] == i)
              {
                subele.push_back(k);
                break;
              }
          }
          if ((int)subele.size() != 3) dserror("Node not attached to exactly 3 subelements");

          masterele[ele->Id()] = subele;

          // no longer need to look at this element
          break;
        }
      }
      if (!foundit) dserror("Weired, this adjele seems not attached to me");
    }  // for (unsigned j=0; j<myadjele.size(); ++j)
    if (masterele.size() != myadjele.size()) dserror("subelement connectivity wrong");

    adjsubele[nodeidL] = masterele;

    //-----------------------------------------------------------------
    // for each adjele and its subele, build local connectivity
    std::vector<std::vector<std::vector<int>>> lmlm((int)myadjele.size());
    for (unsigned j = 0; j < myadjele.size(); ++j)
    {
      DRT::ELEMENTS::NStet5* ele = myadjele[j];
      std::vector<int>& subele = masterele[ele->Id()];
      lmlm[j].resize((int)subele.size());
      for (unsigned k = 0; k < subele.size(); ++k)
      {
        const int subeleid = subele[k];
        const int* sublm = ele->SubLM(subeleid);
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
        if ((int)elelm.size() != 12) dserror("Subelement does not have 12 dofs");
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
int DRT::ELEMENTS::NStet5Type::Initialize(DRT::Discretization& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStet5Type::Initialize");

  const int myrank = dis.Comm().MyPID();
  const int numproc = dis.Comm().NumProc();

  //----------------------------------------------------------------------
  // init elements, make maps of column elements and row nodes
  InitElementsandMaps(elecids_, noderids_, myrank, numproc, dis);

  //----------------------------------------------------------------------
  // compute adjacency for each row node
  // make patch of adjacent elements
  // make patch of adjacent nodes (including center node itself)
  // make lm for nodal patch
  InitAdjacency(elecids_, noderids_, adjele_, adjnode_, adjlm_, adjsubele_, lmlm_, dis);


  return 0;
}
