/*----------------------------------------------------------------------*/
/*!
\brief discsh3 element

\level 3

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------*/
#include "discsh3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"


DRT::ELEMENTS::DiscSh3Type DRT::ELEMENTS::DiscSh3Type::instance_;

DRT::ELEMENTS::DiscSh3Type& DRT::ELEMENTS::DiscSh3Type::Instance() { return instance_; }


DRT::ParObject* DRT::ELEMENTS::DiscSh3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::DiscSh3* object = new DRT::ELEMENTS::DiscSh3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::DiscSh3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "DISCSH3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::DiscSh3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::DiscSh3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::DiscSh3(id, owner));
  return ele;
}


void DRT::ELEMENTS::DiscSh3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // Number of DOFs per node
  nv = 3;     // Number of velocity DOFs, i.e. time derivative of DOFs per node
  dimns = 6;  // Number of dimensions. 3 Translations, 3 rotations, Global dofs
}

void DRT::ELEMENTS::DiscSh3Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::DiscSh3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["DISCSH3"];

  defs["TRI3"].AddIntVector("TRI3", 3).AddNamedInt("MAT").AddNamedDouble("THICK");
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::DiscSh3LineType::Create(const int id, const int owner)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                       mukherjee 08/15 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3::DiscSh3(int id, int owner)
    : DRT::Element(id, owner),
      thickness_(0.0),
      ngptri_(0),
      nhyb_(0),
      ans_(0),
      sdc_(1.0),
      material_(0),
      x_n_(LINALG::Matrix<1, 9>(true)),
      x_n_1_(LINALG::Matrix<1, 9>(true)),
      data_()
{
  //  gaussrule_ = DRT::UTILS::intrule_tri_1point;
  //  gaussrule_ = DRT::UTILS::intrule_tri_3point;
  gaussrule_ = DRT::UTILS::intrule_tri_3point_gauss_radau;
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                  mukherjee 08/15 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3::DiscSh3(const DRT::ELEMENTS::DiscSh3& old)
    : DRT::Element(old),
      forcetype_(old.forcetype_),
      thickness_(old.thickness_),
      ngptri_(old.ngptri_),
      nhyb_(old.nhyb_),
      ans_(old.ans_),
      sdc_(old.sdc_),
      material_(old.material_),
      x_n_(old.x_n_),
      x_n_1_(old.x_n_1_),
      gaussrule_(old.gaussrule_),
      data_(old.data_)
{
  for (int i = 0; i < 3; ++i) ngp_[i] = old.ngp_[i];
  for (int i = 0; i < 5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*-----------------------------------------------------------------------*
 |  Deep copy this instance of DiscSh3 and return pointer to it (public) |
 |                                                       mukherjee 08/15 |
 *-----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::DiscSh3::Clone() const
{
  DRT::ELEMENTS::DiscSh3* newelement = new DRT::ELEMENTS::DiscSh3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                       mukherjee 08/15|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::DiscSh3::Shape() const
{
  switch (NumNode())
  {
    case 3:
      return tri3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // forcetype_
  AddtoPack(data, forcetype_);
  // thickness_
  AddtoPack(data, thickness_);
  // ngp_
  AddtoPack(data, ngp_, 3 * sizeof(int));
  // ngptri_
  AddtoPack(data, ngptri_);
  // nhyb_
  AddtoPack(data, nhyb_);
  // eas_
  AddtoPack(data, eas_, 5 * sizeof(int));
  // ans_
  AddtoPack(data, ans_);
  // sdc_
  AddtoPack(data, sdc_);
  // material_
  AddtoPack(data, material_);

  AddtoPack<1, 9>(data, x_n_);

  AddtoPack<1, 9>(data, x_n_1_);

  // data_
  AddtoPack(data, data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::Unpack(const std::vector<char>& data)
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
  // forcetype_
  forcetype_ = static_cast<ForceType>(ExtractInt(position, data));
  // thickness_
  ExtractfromPack(position, data, thickness_);
  // ngp_
  ExtractfromPack(position, data, ngp_, 3 * sizeof(int));
  // ngptri_
  ExtractfromPack(position, data, ngptri_);
  // nhyb_
  ExtractfromPack(position, data, nhyb_);
  // eas_
  ExtractfromPack(position, data, eas_, 5 * sizeof(int));
  // ans_
  ExtractfromPack(position, data, ans_);
  // sdc_
  ExtractfromPack(position, data, sdc_);
  //   material_
  ExtractfromPack(position, data, material_);

  ExtractfromPack<1, 9>(position, data, x_n_);

  ExtractfromPack<1, 9>(position, data, x_n_1_);

  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                        mukherjee 08/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3::~DiscSh3() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                          mukherjee 08/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::Print(std::ostream& os) const
{
  os << "DiscSh3 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         mukherjee 08/15|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::DiscSh3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<DiscSh3Line, DiscSh3>(DRT::UTILS::buildLines, this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                      mukherjee 08/15|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::DiscSh3::Surfaces()
{
  std::vector<Teuchos::RCP<Element>> surfaces(1);
  surfaces[0] = Teuchos::rcp(this, false);
  return surfaces;
}


/*----------------------------------------------------------------------------*
 |  Calculate spatial configuration of an element (private)    mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 9> DRT::ELEMENTS::DiscSh3::SpatialConfiguration(
    const std::vector<double>& disp) const
{
  LINALG::Matrix<1, 9> x(true);
  for (int i = 0; i < NumNode(); ++i)
  {
    x(3 * i) = Nodes()[i]->X()[0] + disp[i * 3 + 0];
    x(3 * i + 1) = Nodes()[i]->X()[1] + disp[i * 3 + 1];
    x(3 * i + 2) = Nodes()[i]->X()[2] + disp[i * 3 + 2];
  }
  return x;
}


/*----------------------------------------------------------------------------*
 |  Calculate spatial configuration of an element (private)    mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 9> DRT::ELEMENTS::DiscSh3::SpatialConfiguration(DRT::Discretization& dis) const
{
  Teuchos::RCP<const Epetra_Vector> discol = dis.GetState("displacement");

  if (discol == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
  LINALG::Matrix<1, 9> coord(true);

  // compute current nodal positions
  for (int dim = 0; dim < 3; ++dim)
  {
    for (int node = 0; node < NumNode(); ++node)
    {
      double referenceposition = ((Nodes())[node])->X()[dim];
      std::vector<int> dofnode = dis.Dof((Nodes())[node]);
      double displacement = (double)(*discol)[dis.DofColMap()->LID(dofnode[dim])];
      coord(3 * node + dim) = referenceposition + displacement;
    }
  }
  return coord;
}


/*----------------------------------------------------------------------------*
 |  Calculate spatial configuration of an element (private)    mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 9> DRT::ELEMENTS::DiscSh3::SpatialConfiguration(
    DRT::Discretization& dis, const Epetra_Vector& discol) const
{
  LINALG::Matrix<1, 9> coord(true);

  // compute current nodal positions
  for (int dim = 0; dim < 3; ++dim)
  {
    for (int node = 0; node < NumNode(); ++node)
    {
      double referenceposition = ((Nodes())[node])->X()[dim];
      std::vector<int> dofnode = dis.Dof((Nodes())[node]);
      double displacement = (double)(discol)[dis.DofColMap()->LID(dofnode[dim])];
      coord(3 * node + dim) = referenceposition + displacement;
    }
  }
  return coord;
}


/*----------------------------------------------------------------------------*
 |  Calculate current velocity of an element (private)         mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 9> DRT::ELEMENTS::DiscSh3::GetVel(DRT::Discretization& dis) const
{
  Teuchos::RCP<const Epetra_Vector> velcol = dis.GetState("velocity");

  LINALG::Matrix<1, 9> vel(true);

  // compute current nodal positions
  for (int dim = 0; dim < 3; ++dim)
  {
    for (int node = 0; node < NumNode(); ++node)
    {
      std::vector<int> dofnode = dis.Dof((Nodes())[node]);
      double velocity = (double)(*velcol)[dis.DofColMap()->LID(dofnode[dim])];
      vel(3 * node + dim) = velocity;
    }
  }
  return vel;
}


/*--------------------------------------------------------------------------------*
 | Add primary DOFs to master element                              mukherjee 09/15|
 *--------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::AddPrimaryDOFsMaster(DRT::Element& master, DRT::Element& neighbour,
    std::vector<int>& connectivity, std::vector<FAD>& x_FAD, DRT::Discretization& dis,
    bool refconfig)
{
  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(dis);
    int count = 0;
    for (int j = 0; j < master.NumNode(); j++)
    {
      if ((master.NodeIds()[j]) != (neighbour.NodeIds()[0]) &&
          (master.NodeIds()[j]) != (neighbour.NodeIds()[1]) &&
          (master.NodeIds()[j]) != (neighbour.NodeIds()[2]))
      {
        for (int dof = 0; dof < 3; dof++)
        {
          x_FAD[dof] = x(3 * j + dof);
          x_FAD[dof].diff(dof, NUMDOF_DISCSH3 + 3);
        }
        const int node0 = 0;
        connectivity.push_back(node0);  // 1
      }
      else if ((master.NodeIds()[j]) == (neighbour.NodeIds()[0]) ||
               (master.NodeIds()[j]) == (neighbour.NodeIds()[1]) ||
               (master.NodeIds()[j]) == (neighbour.NodeIds()[2]))
      {
        if (count == 0)
        {
          for (int dof = 0; dof < 3; dof++)
          {
            x_FAD[dof + 6] = x(3 * j + dof);
            x_FAD[dof + 6].diff(dof + 6, NUMDOF_DISCSH3 + 3);
          }
          const int node3 = 2;
          connectivity.push_back(node3);  // 3
        }
        else if (count == 1)
        {
          for (int dof = 0; dof < 3; dof++)
          {
            x_FAD[dof + 9] = x(3 * j + dof);
            x_FAD[dof + 9].diff(dof + 9, NUMDOF_DISCSH3 + 3);
          }
          const int node4 = 3;
          connectivity.push_back(node4);  // 4
        }
        else if (count == 2)
          dserror("invalid value of count!");

        count++;  //
      }
      else
      {
        dserror("Node not found");
      }
    }
  }

  return;
}

/*-------------------------- ----------------------------------------------------*
 | Add primary DOFs to slave element                              mukherjee 09/15|
 *-------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::AddPrimaryDOFsNeighbour(DRT::Element& master, DRT::Element& neighbour,
    std::vector<FAD>& x_FAD, DRT::Discretization& dis, bool refconfig)
{
  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(dis);

    for (int j = 0; j < neighbour.NumNode(); j++)
    {
      if ((neighbour.NodeIds()[j]) != (master.NodeIds()[0]) &&
          (neighbour.NodeIds()[j]) != (master.NodeIds()[1]) &&
          (neighbour.NodeIds()[j]) != (master.NodeIds()[2]))
      {
        for (int dof = 0; dof < 3; dof++)
        {
          x_FAD[dof + 3] = x(3 * j + dof);
          x_FAD[dof + 3].diff(dof + 3, NUMDOF_DISCSH3 + 3);
        }
      }
    }
  }

  return;
}

/*---------------------------------------------------------------------------------------*
 | Sort primary DOFs of master & slave elements                           mukherjee 09/15|
 *---------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::SortPrimaryDOFs(DRT::Element& master, DRT::Element& neighbour,
    std::vector<FAD>& x_FAD, std::vector<FAD>& x_FAD_master, std::vector<FAD>& x_FAD_neighbour)
{
  // It is important that Triangle 134 belongs to the master element
  // where 3-4 is the shared edge. What is also important that for neigbouring element
  // we change connectivity to Triangle 243
  for (int dof = 0; dof < 3; dof++)
  {
    // Triangle 134
    x_FAD_master[dof] = x_FAD[dof];
    x_FAD_master[dof + 3] = x_FAD[dof + 6];
    x_FAD_master[dof + 6] = x_FAD[dof + 9];

    // Triangle 243
    x_FAD_neighbour[dof] = x_FAD[dof + 3];
    x_FAD_neighbour[dof + 3] = x_FAD[dof + 9];
    x_FAD_neighbour[dof + 6] = x_FAD[dof + 6];
  }

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Check if surface normal of the element is poining outwords (private)    mukherjee 04/07|
 *-----------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CheckIfOutwardsNormal(
    Teuchos::ParameterList& params, const int NumGElements)
{
  LINALG::TMatrix<FAD, 1, 3> normal(true);
  LINALG::Matrix<1, 3> normal_val(true);

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  LINALG::Matrix<1, 9> X = MaterialConfiguration();
  for (int j = 0; j < 3; j++)
  {
    side1(j) = X(j + 3) - X(j);
    side2(j) = X(j + 6) - X(j);
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = pow(crossprod.Dot(crossprod), 0.5);

  for (int i = 0; i < 3; i++)
  {
    normal(i) = crossprod(i) / norm_crossprod;
    normal_val(i) = normal(i).val();
  }

  // Check if the normal vector is outwards
  Teuchos::RCP<Epetra_SerialDenseVector> CG_ref_rcp =
      params.get<Teuchos::RCP<Epetra_SerialDenseVector>>("reference CG", Teuchos::null);
  LINALG::Matrix<1, 3> CG(true);
  LINALG::Matrix<1, 3> Barycenter(true);
  for (int i = 0; i < 3; i++)
  {
    CG(i) = (*CG_ref_rcp)(i) / NumGElements;
  }


  // Calc barycenter
  LINALG::Matrix<1, 3> bary_center;
  for (int dim = 0; dim < 3; dim++)
  {
    bary_center(dim) = (X(dim) + X(3 + dim) + X(6 + dim)) / 3;
  }

  // Vector connecting center of gravity and barycenter
  LINALG::Matrix<1, 3> vector_aux;
  for (int dim = 0; dim < 3; dim++)
  {
    vector_aux(dim) = bary_center(dim) - CG(dim);
  }
  //  std::cout<<"vector_aux="<<vector_aux<<std::endl;

  double dot_product = vector_aux.Dot(normal_val);

  if (dot_product < 0)
  {
    dserror("Normal vector pointing inwards!");
  }



  return;
}

/*----------------------------------------------------------------------*
 * Compute surface area and its first and second derivatives    lw 05/08*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::ComputeAreaDeriv(const LINALG::SerialDenseMatrix& x, const int numnode,
    const int ndof, double& A, Teuchos::RCP<Epetra_SerialDenseVector> Adiff,
    Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2)
{
  // initialization
  A = 0.;
  Adiff->Size(ndof);

  if (Adiff2 != Teuchos::null) Adiff2->Shape(ndof, ndof);

  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);

  int ngp = intpoints.nquad;

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseMatrix deriv(2, numnode);
  LINALG::SerialDenseMatrix dxyzdrs(2, 3);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid = 0; gpid < ngp; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, Shape());

    std::vector<double> normal(3);
    double detA;
    SurfaceIntegration(detA, normal, x, deriv);
    A += detA * intpoints.qwgt[gpid];

    LINALG::SerialDenseMatrix ddet(3, ndof, true);
    LINALG::SerialDenseMatrix ddet2(3 * ndof, ndof, true);
    LINALG::SerialDenseVector jacobi_deriv(ndof, true);

    dxyzdrs.Multiply('N', 'N', 1.0, deriv, x, 0.0);

    /*--------------- derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int i = 0; i < numnode; ++i)
    {
      ddet(0, 3 * i) = 0.;
      ddet(0, 3 * i + 1) = deriv(0, i) * dxyzdrs(1, 2) - deriv(1, i) * dxyzdrs(0, 2);
      ddet(0, 3 * i + 2) = deriv(1, i) * dxyzdrs(0, 1) - deriv(0, i) * dxyzdrs(1, 1);

      ddet(1, 3 * i) = deriv(1, i) * dxyzdrs(0, 2) - deriv(0, i) * dxyzdrs(1, 2);
      ddet(1, 3 * i + 1) = 0.;
      ddet(1, 3 * i + 2) = deriv(0, i) * dxyzdrs(1, 0) - deriv(1, i) * dxyzdrs(0, 0);

      ddet(2, 3 * i) = deriv(0, i) * dxyzdrs(1, 1) - deriv(1, i) * dxyzdrs(0, 1);
      ddet(2, 3 * i + 1) = deriv(1, i) * dxyzdrs(0, 0) - deriv(0, i) * dxyzdrs(1, 0);
      ddet(2, 3 * i + 2) = 0.;

      jacobi_deriv(i * 3) = 1 / detA * (normal[2] * ddet(2, 3 * i) + normal[1] * ddet(1, 3 * i));
      jacobi_deriv(i * 3 + 1) =
          1 / detA * (normal[2] * ddet(2, 3 * i + 1) + normal[0] * ddet(0, 3 * i + 1));
      jacobi_deriv(i * 3 + 2) =
          1 / detA * (normal[0] * ddet(0, 3 * i + 2) + normal[1] * ddet(1, 3 * i + 2));
    }

    /*--- calculation of first derivatives of current interfacial area
     *----------------------------- with respect to the displacements */
    for (int i = 0; i < ndof; ++i)
    {
      (*Adiff)[i] += jacobi_deriv(i) * intpoints.qwgt[gpid];
    }

    if (Adiff2 != Teuchos::null)
    {
      /*--------- second derivates of minor determiants of the Jacobian
       *----------------------------- with respect to the displacements */
      for (int n = 0; n < numnode; ++n)
      {
        for (int o = 0; o < numnode; ++o)
        {
          ddet2(n * 3 + 1, o * 3 + 2) = deriv(0, n) * deriv(1, o) - deriv(1, n) * deriv(0, o);
          ddet2(n * 3 + 2, o * 3 + 1) = -ddet2(n * 3 + 1, o * 3 + 2);

          ddet2(ndof + n * 3, o * 3 + 2) = deriv(1, n) * deriv(0, o) - deriv(0, n) * deriv(1, o);
          ddet2(ndof + n * 3 + 2, o * 3) = -ddet2(ndof + n * 3, o * 3 + 2);

          ddet2(2 * ndof + n * 3, o * 3 + 1) = ddet2(n * 3 + 1, o * 3 + 2);
          ddet2(2 * ndof + n * 3 + 1, o * 3) = -ddet2(2 * ndof + n * 3, o * 3 + 1);
        }
      }

      /*- calculation of second derivatives of current interfacial areas
       *----------------------------- with respect to the displacements */
      for (int i = 0; i < ndof; ++i)
      {
        int var1, var2;

        if (i % 3 == 0)  // displacement in x-direction
        {
          var1 = 1;
          var2 = 2;
        }
        else if ((i - 1) % 3 == 0)  // displacement in y-direction
        {
          var1 = 0;
          var2 = 2;
        }
        else if ((i - 2) % 3 == 0)  // displacement in z-direction
        {
          var1 = 0;
          var2 = 1;
        }
        else
        {
          dserror("calculation of second derivatives of interfacial area failed");
          exit(1);
        }

        for (int j = 0; j < ndof; ++j)
        {
          (*Adiff2)(i, j) +=
              (-1 / detA * jacobi_deriv(j) * jacobi_deriv(i) +
                  1 / detA *
                      (ddet(var1, i) * ddet(var1, j) + normal[var1] * ddet2(var1 * ndof + i, j) +
                          ddet(var2, i) * ddet(var2, j) +
                          normal[var2] * ddet2(var2 * ndof + i, j))) *
              intpoints.qwgt[gpid];
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * Compute surface area at ref                                  lw 05/08*
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::ComputeAreaRef(
    const LINALG::SerialDenseMatrix& x0, const int numnode, const int ndof, double& A)
{
  // initialization
  A = 0.;

  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);

  int ngp = intpoints.nquad;

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseMatrix deriv(2, numnode);
  LINALG::SerialDenseMatrix dxyzdrs(2, 3);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid = 0; gpid < ngp; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, Shape());

    std::vector<double> normal(3);
    double detA;
    SurfaceIntegration(detA, normal, x0, deriv);
    A += detA * intpoints.qwgt[gpid];
  }


  return;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      gee 04/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::SurfaceIntegration(double& detA, std::vector<double>& normal,
    const Epetra_SerialDenseMatrix& x, const Epetra_SerialDenseMatrix& deriv)
{
  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2, 3);
  dxyzdrs.Multiply('N', 'N', 1.0, deriv, x, 0.0);

  /* compute covariant metric tensor G for surface element
  **                        | g11   g12 |
  **                    G = |           |
  **                        | g12   g22 |
  ** where (o denotes the inner product, xyz a vector)
  **
  **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
  ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
  **        dr     dr            dr     ds            ds     ds
  */
  LINALG::SerialDenseMatrix metrictensor(2, 2);
  metrictensor.Multiply('N', 'T', 1.0, dxyzdrs, dxyzdrs, 0.0);
  detA = sqrt(metrictensor(0, 0) * metrictensor(1, 1) - metrictensor(0, 1) * metrictensor(1, 0));
  normal[0] = dxyzdrs(0, 1) * dxyzdrs(1, 2) - dxyzdrs(0, 2) * dxyzdrs(1, 1);
  normal[1] = dxyzdrs(0, 2) * dxyzdrs(1, 0) - dxyzdrs(0, 0) * dxyzdrs(1, 2);
  normal[2] = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(0, 1) * dxyzdrs(1, 0);

  return;
}

/*-----------------------------------------------------------------------*
 |  get internal face  element (public)                  mukherjee 06 /15|
 *-----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::DiscSh3::CreateFaceElement(
    DRT::Element* parent_slave,            //!< parent slave fluid3 element
    int nnode,                             //!< number of surface nodes
    const int* nodeids,                    //!< node ids of surface element
    DRT::Node** nodes,                     //!< nodes of surface element
    const int lsurface_master,             //!< local surface number w.r.t master parent element
    const int lsurface_slave,              //!< local surface number w.r.t slave parent element
    const std::vector<int>& localtrafomap  //! local trafo map
)
{
  // dynamic cast for slave parent element
  DRT::ELEMENTS::DiscSh3* slave_pele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(parent_slave);

  // insert both parent elements
  return DRT::UTILS::ElementIntFaceFactory<DiscSh3Line, DiscSh3>(-1,  //!< internal face element id
      -1,               //!< owner of internal face element
      nnode,            //!< number of surface nodes
      nodeids,          //!< node ids of surface element
      nodes,            //!< nodes of surface element
      this,             //!< master parent element
      slave_pele,       //!< slave parent element
      lsurface_master,  //!< local surface number w.r.t master parent element
      lsurface_slave,   //!< local surface number w.r.t slave parent element
      localtrafomap     //!< local trafo map
  );
}
