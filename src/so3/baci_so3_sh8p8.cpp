/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_so3_sh8p8.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_so3_nullspace.hpp"
#include "baci_so3_utils.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN


DRT::ELEMENTS::So_sh8p8Type DRT::ELEMENTS::So_sh8p8Type::instance_;

DRT::ELEMENTS::So_sh8p8Type& DRT::ELEMENTS::So_sh8p8Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::So_sh8p8Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So_sh8p8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8p8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8p8(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8p8Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8p8(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_sh8p8Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
  np = 1;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::So_sh8p8Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  dserror("method ComputeNullSpace not implemented!");
  return nullspace;
}

void DRT::ELEMENTS::So_sh8p8Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("STAB")
                     .AddNamedString("ANS")
                     .AddNamedString("LIN")
                     .AddNamedString("THICKDIR")
                     .AddNamedString("EAS")
                     .AddNamedString("ISO")
                     .AddOptionalNamedDoubleVector("RAD", 3)
                     .AddOptionalNamedDoubleVector("AXI", 3)
                     .AddOptionalNamedDoubleVector("CIR", 3)
                     .AddOptionalNamedDouble("STRENGTH")
                     .Build();
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 03/09|
 *----------------------------------------------------------------------*/
// 9-Voigt C-index                                         0 1 2  3 4 5  6 7 8
// TODO this notation does not fit our usual definitions
const int DRT::ELEMENTS::So_sh8p8::VOIGT9ROW_INCONSISTENT_[NUMDFGR_] = {0, 1, 2, 0, 1, 2, 0, 2, 1};
const int DRT::ELEMENTS::So_sh8p8::VOIGT9COL_INCONSISTENT_[NUMDFGR_] = {0, 1, 2, 1, 2, 0, 2, 1, 0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 9-Voigt C-indices    0   3   6   8   1   4   5   7   2
// TODO this notation does not fit our usual definitions
const int DRT::ELEMENTS::So_sh8p8::VOIGT3X3NONSYM_INCONSISTENT_[NUMDFGR_] = {
    0, 3, 6, 8, 1, 4, 5, 7, 2};

// 24 displacement and 8 pressure DOFs into 32 total element DOFs
const int DRT::ELEMENTS::So_sh8p8::DISPTODISPPRES_[NUMDISP_] = {
    0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 29, 30};
const int DRT::ELEMENTS::So_sh8p8::PRESTODISPPRES_[NUMPRES_] = {3, 7, 11, 15, 19, 23, 27, 31};


/*----------------------------------------------------------------------*
 |  ctor (public)                                            bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8p8::So_sh8p8(int id, int owner) : DRT::ELEMENTS::So_sh8(id, owner)
{
  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8p8::So_sh8p8(const DRT::ELEMENTS::So_sh8p8& old) : DRT::ELEMENTS::So_sh8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh8p8::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So_sh8p8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_sh8 Element
  DRT::ELEMENTS::So_sh8::Pack(data);
  // techniques
  AddtoPack(data, stab_);
  AddtoPack(data, ans_);
  AddtoPack(data, lin_);
  AddtoPack(data, iso_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_sh8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::So_sh8::Unpack(basedata);
  // techniques
  stab_ = static_cast<StabilisationType>(ExtractInt(position, data));
  ans_ = static_cast<AnsType>(ExtractInt(position, data));
  lin_ = static_cast<LinearizationType>(ExtractInt(position, data));
  iso_ = static_cast<IsochoricType>(ExtractInt(position, data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              bborn 03/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Print(std::ostream& os) const
{
  os << "So_sh8p8 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      tk 04/09   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::sosh8p8_expol(
    CORE::LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  // static variables, that are the same for every element
  static CORE::LINALG::Matrix<NUMNOD_, NUMGPT_> expol;
  static bool isfilled;


  if (isfilled == false)
  {
    double sq3 = sqrt(3.0);
    expol(0, 0) = 1.25 + 0.75 * sq3;
    expol(0, 1) = -0.25 - 0.25 * sq3;
    expol(0, 2) = -0.25 + 0.25 * sq3;
    expol(0, 3) = -0.25 - 0.25 * sq3;
    expol(0, 4) = -0.25 - 0.25 * sq3;
    expol(0, 5) = -0.25 + 0.25 * sq3;
    expol(0, 6) = 1.25 - 0.75 * sq3;
    expol(0, 7) = -0.25 + 0.25 * sq3;
    expol(1, 1) = 1.25 + 0.75 * sq3;
    expol(1, 2) = -0.25 - 0.25 * sq3;
    expol(1, 3) = -0.25 + 0.25 * sq3;
    expol(1, 4) = -0.25 + 0.25 * sq3;
    expol(1, 5) = -0.25 - 0.25 * sq3;
    expol(1, 6) = -0.25 + 0.25 * sq3;
    expol(1, 7) = 1.25 - 0.75 * sq3;
    expol(2, 2) = 1.25 + 0.75 * sq3;
    expol(2, 3) = -0.25 - 0.25 * sq3;
    expol(2, 4) = 1.25 - 0.75 * sq3;
    expol(2, 5) = -0.25 + 0.25 * sq3;
    expol(2, 6) = -0.25 - 0.25 * sq3;
    expol(2, 7) = -0.25 + 0.25 * sq3;
    expol(3, 3) = 1.25 + 0.75 * sq3;
    expol(3, 4) = -0.25 + 0.25 * sq3;
    expol(3, 5) = 1.25 - 0.75 * sq3;
    expol(3, 6) = -0.25 + 0.25 * sq3;
    expol(3, 7) = -0.25 - 0.25 * sq3;
    expol(4, 4) = 1.25 + 0.75 * sq3;
    expol(4, 5) = -0.25 - 0.25 * sq3;
    expol(4, 6) = -0.25 + 0.25 * sq3;
    expol(4, 7) = -0.25 - 0.25 * sq3;
    expol(5, 5) = 1.25 + 0.75 * sq3;
    expol(5, 6) = -0.25 - 0.25 * sq3;
    expol(5, 7) = -0.25 + 0.25 * sq3;
    expol(6, 6) = 1.25 + 0.75 * sq3;
    expol(6, 7) = -0.25 - 0.25 * sq3;
    expol(7, 7) = 1.25 + 0.75 * sq3;

    for (int i = 0; i < NUMNOD_; ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        expol(i, j) = expol(j, i);
      }
    }
    isfilled = true;
  }

  CORE::LINALG::Matrix<NUMNOD_, MAT::NUM_STRESS_3D> nodalstresses;

  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_; ++i)
  {
    const int lid = expolstresses.Map().LID(NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = 1.0 / Nodes()[i]->NumElement();
      for (int j = 0; j < MAT::NUM_STRESS_3D; ++j)
        (*(expolstresses(j)))[lid] += nodalstresses(i, j) * invmyadjele;
    }
  }
  return;
}

BACI_NAMESPACE_CLOSE
