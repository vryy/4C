/*-----------------------------------------------------------*/
/*! \file
\brief three dimensional spring element


\level 3
*/
/*-----------------------------------------------------------*/

#include "spring3.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"

DRT::ELEMENTS::Spring3Type DRT::ELEMENTS::Spring3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Spring3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Spring3* object = new DRT::ELEMENTS::Spring3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Spring3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SPRING3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Spring3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Spring3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Spring3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Spring3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::Spring3Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Spring3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SPRING3"];

  defs["LINE2"].AddIntVector("LINE2", 2).AddNamedInt("MAT");
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Spring3::Spring3(int id, int owner)
    : DRT::Element(id, owner),
      data_(),
      isinit_(false),
      diff_disp_ref_(LINALG::Matrix<1, 3>(true)),
      deltatheta_(LINALG::Matrix<1, 3>(true)),
      material_(0),
      lrefe_(0),
      lcurr_(0),
      crosssec_(0),
      NormMoment(0),
      NormForce(0),
      RatioNormForceMoment(0),
      Theta0_(LINALG::Matrix<3, 1>(true)),
      Theta_(LINALG::Matrix<3, 1>(true))
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Spring3::Spring3(const DRT::ELEMENTS::Spring3& old)
    : DRT::Element(old),
      data_(old.data_),
      isinit_(old.isinit_),
      X_(old.X_),
      trefNode_(old.trefNode_),
      ThetaRef_(old.ThetaRef_),
      diff_disp_ref_(old.diff_disp_ref_),
      deltatheta_(old.deltatheta_),
      Qnew_(old.Qnew_),
      Qold_(old.Qold_),
      dispthetanew_(old.dispthetanew_),
      dispthetaold_(old.dispthetaold_),
      material_(old.material_),
      lrefe_(old.lrefe_),
      lcurr_(old.lcurr_),
      jacobimass_(old.jacobimass_),
      jacobinode_(old.jacobinode_),
      crosssec_(old.crosssec_),
      NormMoment(old.NormMoment),
      NormForce(old.NormForce),
      RatioNormForceMoment(old.RatioNormForceMoment),
      Theta0_(old.Theta0_),
      Theta_(old.Theta_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Spring3 and return pointer to it (public)|
 |                                                       mukherjee 04/15|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Spring3::Clone() const
{
  DRT::ELEMENTS::Spring3* newelement = new DRT::ELEMENTS::Spring3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                       mukherjee 04/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Spring3::~Spring3() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                         mukherjee 04/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::Print(std::ostream& os) const
{
  os << "Spring3 ";
  Element::Print(os);
  return;
}



/*----------------------------------------------------------------------*
 | Print the change in angle of this element            mukherjee 04/15 |
 *----------------------------------------------------------------------*/
LINALG::Matrix<1, 3> DRT::ELEMENTS::Spring3::DeltaTheta() const { return deltatheta_; }

/*----------------------------------------------------------------------*
 |(public)                                               mukherjee 04/15|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Spring3::Shape() const { return line2; }

void DRT::ELEMENTS::Spring3::GetCurrTangents(
    std::vector<double>& disp, std::vector<LINALG::Matrix<3, 1>>& Tcurr)
{
  // rotational displacement at a certain node between this and last iteration step
  LINALG::Matrix<3, 1> deltatheta(true);
  // rotational displacement at a certain node between this and last iteration step in quaternion
  // form
  LINALG::Matrix<4, 1> deltaQ(true);
  // Compute current nodal triads
  for (int node = 0; node < 2; node++)
  {
    /*rotation increment relative to configuration in last iteration step is difference between
     *current rotation entry in displacement vector minus rotation entry in displacement vector in
     *last iteration step*/

    // This shift is necessary, since our beam formulation and the corresponding linearization is
    // based on multiplicative increments, while in BACI (Newtonfull()) the displacement of the last
    // iteration and the rotation increment of the current iteration are added in an additive
    // manner. This step is reversed in the next two lines in order to recover the multiplicative
    // rotation increment between the last and the current Newton step. This also means that
    // dispthetanew_ has no physical meaning since it is the additive sum of non-additive rotation
    // increments(meier, 03.2014)
    for (int i = 0; i < 3; i++) dispthetanew_[node](i) = disp[6 * node + 3 + i];

    deltatheta = dispthetanew_[node];
    deltatheta -= dispthetaold_[node];

    // compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta, deltaQ);
    // multiply relative rotation with rotation in last configuration to get rotation in new
    // configuration
    LARGEROTATIONS::quaternionproduct(Qold_[node], deltaQ, Qnew_[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and
    // intricate calculations
    Qnew_[node].Scale(1 / Qnew_[node].Norm2());

    LINALG::Matrix<3, 3> Triad(true);
    LARGEROTATIONS::quaterniontotriad(Qnew_[node], Triad);
    for (int i = 0; i < 3; i++)
      Tcurr[node](i) = Triad(i, 0);  // reference CheckOrientation StatMechManager

    Tcurr[node].Scale(1 / Tcurr[node].Norm2());
  }

  return;
}  // DRT::ELEMENTS::Beam3::UpdateNodalTriad

// brief! Return current tangent of beam3r elements connected to beam3 element
void DRT::ELEMENTS::Spring3::TcurrBeam3r(LINALG::Matrix<3, 1>& Tcurr1, LINALG::Matrix<3, 1>& Tcurr2)
{
  DRT::Node* node1 = this->Nodes()[0];
  DRT::Element* Element1 = node1->Elements()[0];
  const DRT::ElementType& eot_el1 = Element1->ElementType();
  if (eot_el1 == DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil1 = dynamic_cast<DRT::ELEMENTS::Beam3r*>(Element1);
    if (fil1 == NULL) return;
    int nodenumber = 0;
    if (node1->Id() != fil1->NodeIds()[0]) nodenumber = 1;
    Tcurr1 = fil1->Tcurr((int)fil1->NodeIds()[nodenumber]);
  }
  DRT::Node* node2 = this->Nodes()[1];
  DRT::Element* Element2 = node2->Elements()[0];
  if (Element2->ElementType() == DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil2 = dynamic_cast<DRT::ELEMENTS::Beam3r*>(Element2);
    if (fil2 == NULL) return;
    int nodenumber = 0;
    if (node1->Id() != fil2->NodeIds()[0]) nodenumber = 1;
    Tcurr2 = fil2->Tcurr((int)fil2->NodeIds()[nodenumber]);
  }
  return;
}

// brief! Return current tangent of beam3r elements connected to beam3 element
void DRT::ELEMENTS::Spring3::TrefBeam3r(LINALG::Matrix<3, 1>& Tref1, LINALG::Matrix<3, 1>& Tref2)
{
  DRT::Node* node1 = this->Nodes()[0];
  DRT::Element* Element1 = node1->Elements()[0];
  const DRT::ElementType& eot_el1 = Element1->ElementType();
  if (eot_el1 == DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil1 = dynamic_cast<DRT::ELEMENTS::Beam3r*>(Element1);
    if (fil1 == NULL) return;
    Tref1 = fil1->Treffirst();
  }
  DRT::Node* node2 = this->Nodes()[1];
  DRT::Element* Element2 = node2->Elements()[0];
  if (Element2->ElementType() == DRT::ELEMENTS::Beam3rType::Instance())
  {
    DRT::ELEMENTS::Beam3r* fil2 = dynamic_cast<DRT::ELEMENTS::Beam3r*>(Element2);
    if (fil2 == NULL) return;
    Tref2 = fil2->Treffirst();
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                       mukherjee 04/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  AddtoPack(data, isinit_);
  AddtoPack<6, 1>(data, X_);
  AddtoPack(data, trefNode_);
  AddtoPack(data, ThetaRef_);
  AddtoPack<1, 3>(data, diff_disp_ref_);
  AddtoPack<1, 3>(data, deltatheta_);
  AddtoPack<4, 1>(data, Qnew_);
  AddtoPack<4, 1>(data, Qold_);
  AddtoPack<3, 1>(data, dispthetanew_);
  AddtoPack<3, 1>(data, dispthetaold_);
  AddtoPack(data, material_);
  AddtoPack(data, lrefe_);
  AddtoPack(data, lcurr_);
  AddtoPack(data, jacobimass_);
  AddtoPack(data, jacobinode_);
  AddtoPack(data, crosssec_);
  AddtoPack(data, NormMoment);
  AddtoPack(data, NormForce);
  AddtoPack(data, RatioNormForceMoment);
  AddtoPack<3, 1>(data, Theta0_);
  AddtoPack<3, 1>(data, Theta_);
  AddtoPack(data, data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                       mukherjee 04/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::Unpack(const std::vector<char>& data)
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
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<6, 1>(position, data, X_);
  ExtractfromPack(position, data, trefNode_);
  ExtractfromPack(position, data, ThetaRef_);
  ExtractfromPack<1, 3>(position, data, diff_disp_ref_);
  ExtractfromPack<1, 3>(position, data, deltatheta_);
  ExtractfromPack<4, 1>(position, data, Qnew_);
  ExtractfromPack<4, 1>(position, data, Qold_);
  ExtractfromPack<3, 1>(position, data, dispthetanew_);
  ExtractfromPack<3, 1>(position, data, dispthetaold_);
  ExtractfromPack(position, data, material_);
  ExtractfromPack(position, data, lrefe_);
  ExtractfromPack(position, data, lcurr_);
  ExtractfromPack(position, data, jacobimass_);
  ExtractfromPack(position, data, jacobinode_);
  ExtractfromPack(position, data, crosssec_);
  ExtractfromPack(position, data, NormMoment);
  ExtractfromPack(position, data, NormForce);
  ExtractfromPack(position, data, RatioNormForceMoment);
  ExtractfromPack<3, 1>(position, data, Theta0_);
  ExtractfromPack<3, 1>(position, data, Theta_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position, data, gausrule_integer);
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                         mukherjee 04/15|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Spring3::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*-----------------------------------------------------------------------------*
 |  Initialize reference Tangents (public)                        mueller 10/12|
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::SetInitialTangents(std::vector<LINALG::Matrix<4, 1>>& initquaternions)
{
  //  if(initquaternions.Norm2()!=0.0)
  {
    Qnew_.resize(NumNode());
    Qold_.resize(NumNode());
    dispthetaold_.resize(NumNode());
    dispthetanew_.resize(NumNode());
    trefNode_.resize(2);

    LINALG::Matrix<3, 3> Triad(true);
    //     initquatern
    for (int node = 0; node < 2; node++)
    {
      for (int j = 0; j < 4; j++)
      {
        Qnew_[node](j) = initquaternions[node](j);
        Qold_[node](j) = initquaternions[node](j);
      }

      LARGEROTATIONS::quaterniontotriad(initquaternions[node], Triad);
      trefNode_[node].Clear();
      for (int i = 0; i < 3; i++)
      {
        trefNode_[node](i) = Triad(i, 0);
      }

      trefNode_[node].Scale(1 / trefNode_[node].Norm2());
    }
  }
  return;
}

void DRT::ELEMENTS::Spring3::SetUpReferenceGeometry(const std::vector<double>& xrefe,
    const std::vector<double>& rotrefe, const bool secondinit, const bool reissner)
{
  if (reissner)
    FilamentIsReissner_ = true;
  else
    FilamentIsReissner_ = false;

  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe; note: the isinit_
   *flag is important for avoiding reinitialization upon restart. However, it should be possible to
   *conduct a second initilization in principle (e.g. for periodic boundary conditions*/
  if (!isinit_ || secondinit)
  {
    isinit_ = true;

    // setting reference coordinates
    for (int i = 0; i < 6; i++) X_(i) = xrefe[i];

    // length in reference configuration
    lrefe_ = std::pow(pow(X_(3) - X_(0), 2) + pow(X_(4) - X_(1), 2) + pow(X_(5) - X_(2), 2), 0.5);

    // set jacobi determinants for integration of mass matrix and at nodes
    jacobimass_.resize(2);
    jacobimass_[0] = lrefe_ / 2.0;
    jacobimass_[1] = lrefe_ / 2.0;
    jacobinode_.resize(2);
    jacobinode_[0] = lrefe_ / 2.0;
    jacobinode_[1] = lrefe_ / 2.0;

    //    if (abs_rotrefe!=0)
    //    {
    // assign size to the vector
    ThetaRef_.resize(3);

    // The reference Tangent for Reissner type of filament
    //      // is set in statmech_manager
    if (!reissner)
    {
      trefNode_.resize(2);
      for (int node = 0; node < 2; node++)
      {
        trefNode_[node].Clear();
        for (int dof = 0; dof < 3; dof++) trefNode_[node](dof) = rotrefe[3 * node + dof];
      }
    }
    diff_disp_ref_.Clear();

    // Calculate reference directional vector of the truss element
    for (int j = 0; j < 3; ++j)
    {
      diff_disp_ref_(j) = Nodes()[1]->X()[j] - Nodes()[0]->X()[j];
    }

    for (int location = 0; location < 3;
         location++)  // Location of torsional spring. There are three locations
    {
      double dotprod = 0.0;
      LINALG::Matrix<1, 3> crossprod(true);
      double CosTheta = 0.0;
      double SinTheta = 0.0;

      if (location == 0)
      {
        double norm_v_ref = diff_disp_ref_.Norm2();
        double norm_t1_ref = trefNode_[location].Norm2();
        for (int j = 0; j < 3; ++j) dotprod += trefNode_[location](j) * diff_disp_ref_(j);

        CosTheta = dotprod / (norm_v_ref * norm_t1_ref);

        // Cross Product
        crossprod(0) =
            trefNode_[location](1) * diff_disp_ref_(2) - trefNode_[location](2) * diff_disp_ref_(1);
        crossprod(1) =
            trefNode_[location](2) * diff_disp_ref_(0) - trefNode_[location](0) * diff_disp_ref_(2);
        crossprod(2) =
            trefNode_[location](0) * diff_disp_ref_(1) - trefNode_[location](1) * diff_disp_ref_(0);

        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_v_ref * norm_t1_ref);
      }

      else if (location == 1)
      {
        double norm_v_ref = diff_disp_ref_.Norm2();
        double norm_t2_ref = trefNode_[location].Norm2();
        for (int j = 0; j < 3; ++j)
          dotprod +=
              trefNode_[location](j) * diff_disp_ref_(j);  // From the opposite direction v_2 =-v_1

        CosTheta = dotprod / (norm_v_ref * norm_t2_ref);

        // cross product
        crossprod(0) =
            trefNode_[location](1) * diff_disp_ref_(2) - trefNode_[location](2) * diff_disp_ref_(1);
        crossprod(1) =
            trefNode_[location](2) * diff_disp_ref_(0) - trefNode_[location](0) * diff_disp_ref_(2);
        crossprod(2) =
            trefNode_[location](0) * diff_disp_ref_(1) - trefNode_[location](1) * diff_disp_ref_(0);
        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_v_ref * norm_t2_ref);
      }

      else  // i.e. for calculation of reference angle between t1 & t2
      {
        double norm_t1_ref = trefNode_[location - 2].Norm2();
        double norm_t2_ref = trefNode_[location - 1].Norm2();
        for (int j = 0; j < 3; ++j)
          dotprod += trefNode_[location - 1](j) * trefNode_[location - 2](j);

        CosTheta = dotprod / (norm_t1_ref * norm_t2_ref);

        // cross product
        crossprod(0) = trefNode_[location - 2](1) * trefNode_[location - 1](2) -
                       trefNode_[location - 2](2) * trefNode_[location - 1](1);
        crossprod(1) = trefNode_[location - 2](2) * trefNode_[location - 1](0) -
                       trefNode_[location - 2](0) * trefNode_[location - 1](2);
        crossprod(2) = trefNode_[location - 2](0) * trefNode_[location - 1](1) -
                       trefNode_[location - 2](1) * trefNode_[location - 1](0);

        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_t1_ref * norm_t2_ref);
      }

      double ThetaBoundary1 = M_PI / 4;
      double ThetaBoundary2 = 3 * M_PI / 4;

      ThetaRef_[location] = 0;
      if (SinTheta >= 0)
      {
        if (CosTheta >= cos(ThetaBoundary1))
          ThetaRef_[location] = asin(SinTheta);
        else if (CosTheta <= cos(ThetaBoundary2))
          ThetaRef_[location] = M_PI - asin(SinTheta);
        else
          ThetaRef_[location] = acos(CosTheta);
      }
      else
        dserror("Angle more than 180 degrees!");

      Theta0_(location) = ThetaRef_[location];
    }
    //    }
  }

  return;
}


int DRT::ELEMENTS::Spring3Type::Initialize(DRT::Discretization& dis)
{
  // reference node positions
  std::vector<double> xrefe;

  // reference nodal tangent positions
  std::vector<double> rotrefe;
  LINALG::Matrix<3, 1> trefNodeAux(true);
  // resize vectors for the number of coordinates we need to store
  xrefe.resize(3 * 2);
  rotrefe.resize(3 * 2);
  for (int i = 0; i < 6; i++) rotrefe[i] = 0;

  // setting beam reference director correctly
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // in case that current element is not a truss3 element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    // if we get so far current element is a truss3 element and  we get a pointer at it
    DRT::ELEMENTS::Spring3* currele = dynamic_cast<DRT::ELEMENTS::Spring3*>(dis.lColElement(i));
    if (!currele) dserror("cast to Spring3* failed");

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++) xrefe[k * 3 + l] = currele->Nodes()[k]->X()[l];
    }

    // ask the spring element about the first element the first node is connected to
    DRT::Element* Element = currele->Nodes()[0]->Elements()[0];
    // Check via dynamic cast, if it's a beam3eb element
    DRT::ELEMENTS::Beam3eb* BeamElement = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(Element);
    if (BeamElement != NULL)
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++)
        {
          trefNodeAux = BeamElement->Tref()[k];
          rotrefe[k * 3 + l] = trefNodeAux(l);
        }
    }

    currele->SetUpReferenceGeometry(xrefe, rotrefe);


  }  // for (int i=0; i<dis_.NumMyColElements(); ++i)


  return 0;
}
