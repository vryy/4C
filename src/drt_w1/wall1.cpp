/*----------------------------------------------------------------------*/
/*!
\file wall1.cpp

\brief ToDo Add meaningful comment.

\level 1

\maintainer Michael Hiermeier

*/
/*----------------------------------------------------------------------*/

#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Wall1Type DRT::ELEMENTS::Wall1Type::instance_;

DRT::ELEMENTS::Wall1Type& DRT::ELEMENTS::Wall1Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Wall1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1* object = new DRT::ELEMENTS::Wall1(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALL")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::Wall1(id, owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1Type::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Wall1(id, owner));
}


void DRT::ELEMENTS::Wall1Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

void DRT::ELEMENTS::Wall1Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure2DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Wall1Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALL"];

  defs["QUAD4"]
      .AddIntVector("QUAD4", 4)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);

  defs["QUAD8"]
      .AddIntVector("QUAD8", 8)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);

  defs["QUAD9"]
      .AddIntVector("QUAD9", 9)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);

  defs["TRI3"]
      .AddIntVector("TRI3", 3)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2)
      //.AddNamedString("STRESSES")
      ;

  defs["TRI6"]
      .AddIntVector("TRI6", 6)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);

  defs["NURBS4"]
      .AddIntVector("NURBS4", 4)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);

  defs["NURBS9"]
      .AddIntVector("NURBS9", 9)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("EAS")
      .AddNamedDouble("THICK")
      .AddNamedString("STRESS_STRAIN")
      .AddNamedIntVector("GP", 2);
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1LineType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Wall1Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/08/|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(int id, int owner)
    : So_base(id, owner),
      data_(),
      material_(0),
      thickness_(0.0),
      old_step_length_(0.0),
      gaussrule_(DRT::UTILS::intrule2D_undefined),
      wtype_(plane_none),
      stresstype_(w1_none),
      iseas_(false),
      eastype_(eas_vague),
      structale_(false),
      distype_(dis_none)
{
  if (DRT::Problem::Instance()->ProblemType() == prb_struct_ale) structale_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(const DRT::ELEMENTS::Wall1& old)
    : So_base(old),
      data_(old.data_),
      material_(old.material_),
      thickness_(old.thickness_),
      old_step_length_(old.old_step_length_),
      gaussrule_(old.gaussrule_),
      wtype_(old.wtype_),
      stresstype_(old.stresstype_),
      iseas_(old.iseas_),
      eastype_(old.eas_vague),
      structale_(old.structale_),
      distype_(old.distype_)
// tsi_couptyp_(old.tsi_couptyp_)

{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1::Clone() const
{
  DRT::ELEMENTS::Wall1* newelement = new DRT::ELEMENTS::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Wall1::Shape() const { return distype_; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);
  // material_
  AddtoPack(data, material_);
  // thickness
  AddtoPack(data, thickness_);
  // plane strain or plane stress information
  AddtoPack(data, wtype_);
  // gaussrule_
  AddtoPack(data, gaussrule_);  // implicit conversion from enum to integer
  // stresstype
  AddtoPack(data, stresstype_);
  // eas
  AddtoPack(data, iseas_);
  // eas type
  AddtoPack(data, eastype_);
  // structale
  AddtoPack(data, structale_);
  // distype
  AddtoPack(data, distype_);
  // data
  AddtoPack(data, data_);
  // line search
  AddtoPack(data, old_step_length_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  So_base::Unpack(basedata);
  // material_
  ExtractfromPack(position, data, material_);
  // thickness_
  ExtractfromPack(position, data, thickness_);
  // plane strain or plane stress information_
  wtype_ = static_cast<DimensionalReduction>(ExtractInt(position, data));
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position, data, gausrule_integer);
  gaussrule_ =
      DRT::UTILS::GaussRule2D(gausrule_integer);  // explicit conversion from integer to enum
  // stresstype_
  stresstype_ = static_cast<StressType>(ExtractInt(position, data));
  // iseas_
  iseas_ = ExtractInt(position, data);
  // eastype_
  eastype_ = static_cast<EasType>(ExtractInt(position, data));
  // structale_
  structale_ = ExtractInt(position, data);
  // distype_
  distype_ = static_cast<DiscretizationType>(ExtractInt(position, data));
  // data
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);
  // line search
  ExtractfromPack(position, data, old_step_length_);
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::~Wall1() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Print(std::ostream& os) const
{
  os << "Wall1 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 07/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line, Wall1>(DRT::UTILS::buildLines, this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1::Surfaces()
{
  std::vector<Teuchos::RCP<Element>> surfaces(1);
  surfaces[0] = Teuchos::rcp(this, false);
  return surfaces;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      popp 08/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_expol(
    Epetra_SerialDenseMatrix& stresses, Epetra_MultiVector& expolstress)
{
  // get gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);
  const int numgp = intpoints.nquad;
  const int numnode = NumNode();
  Epetra_SerialDenseMatrix expol(numnode, numgp);

  const DiscretizationType dt = Shape();
  Epetra_SerialDenseVector funct(numgp);

  // quad4, quad8 and quad9
  if (dt == quad4 or dt == quad8 or dt == quad9)
  {
    // loop over all nodes
    for (int ip = 0; ip < numnode; ++ip)
    {
      // gaussian coordinates
      const double e1 = intpoints.qxg[ip][0];
      const double e2 = intpoints.qxg[ip][1];

      // coordinates of node in the fictitious GP element
      double e1expol;
      double e2expol;

      if (e1 != 0)
        e1expol = 1 / e1;
      else
        e1expol = 0;
      if (e2 != 0)
        e2expol = 1 / e2;
      else
        e2expol = 0;

      // shape functions for the extrapolated coordinates
      switch (numgp)
      {
        case 4:
        {
          DRT::UTILS::shape_function_2D(funct, e1expol, e2expol, quad4);
          break;
        }
        case 9:
        {
          DRT::UTILS::shape_function_2D(funct, e1expol, e2expol, quad9);
          break;
        }
        default:
        {
          dserror("ERROR: Quad4/8/9 nodal stress output only for 4 or 9 Gauss points!");
          break;
        }
      }

      // extrapolation matrix
      for (int i = 0; i < numgp; ++i) expol(ip, i) = funct(i);
    }
  }

  // tri3
  else if (dt == tri3)
  {
    // extrapolation matrix
    for (int ip = 0; ip < numnode; ++ip)
      for (int i = 0; i < numgp; ++i) expol(ip, i) = 1.0 / numgp;
  }

  // tri6
  else if (dt == tri6)
  {
    dserror("ERROR: Nodal stress output is not implemented for TRI6");
    // the following extrapolation is wrong
    /*    // loop over all nodes
        for (int ip=0; ip<numnode; ++ip)
        {
          // gaussian coordinates
          const double e1 = intpoints.qxg[ip][0];
          const double e2 = intpoints.qxg[ip][1];

          // coordinates of node in the fictitious GP element
          double e1expol = 2*e1 - 1.0/3.0;
          double e2expol = 2*e2 - 1.0/3.0;

          // shape functions for the extrapolated coordinates
          switch(numgp)
          {
            case 3:
            {
              DRT::UTILS::shape_function_2D(funct,e1expol,e2expol,tri3);
              break;
            }
            default:
            {
              dserror("ERROR: Tri6 nodal stresses only implemented for 3 Gauss points!");
              break;
            }
          }

          // extrapolation matrix
          for(int i=0;i<numgp;++i) expol(ip,i) = funct(i);
        }*/
  }

  // else
  else
    dserror("extrapolation not implemented for this element type");

  Epetra_SerialDenseMatrix nodalstresses(numnode, Wall1::numstr_);
  nodalstresses.Multiply('N', 'N', 1.0, expol, stresses, 0.0);

  // distribute nodal stresses to expolstress for assembling
  for (int i = 0; i < numnode; ++i)
  {
    int gid = NodeIds()[i];
    if (expolstress.Map().MyGID(NodeIds()[i]))  // rownode
    {
      int adjele = Nodes()[i]->NumElement();
      int lid = expolstress.Map().LID(gid);
      (*(expolstress(0)))[lid] += nodalstresses(i, 0) / adjele;
      (*(expolstress(1)))[lid] += nodalstresses(i, 1) / adjele;
      (*(expolstress(2)))[lid] += nodalstresses(i, 2) / adjele;
      (*(expolstress(3)))[lid] += nodalstresses(i, 3) / adjele;
    }
  }
}

/*-----------------------------------------------------------------------------*
| Map plane Green-Lagrange strains to 3d                       mayr.mt 05/2014 |
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::GreenLagrangePlane3d(
    const Epetra_SerialDenseVector& glplane, LINALG::Matrix<6, 1>& gl3d)
{
  gl3d(0) = glplane(0);               // E_{11}
  gl3d(1) = glplane(1);               // E_{22}
  gl3d(2) = 0.0;                      // E_{33}
  gl3d(3) = glplane(2) + glplane(3);  // 2*E_{12}=E_{12}+E_{21}
  gl3d(4) = 0.0;                      // 2*E_{23}
  gl3d(5) = 0.0;                      // 2*E_{31}

  return;
}
