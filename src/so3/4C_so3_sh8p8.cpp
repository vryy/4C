/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_so3_sh8p8.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::SoSh8p8Type Discret::ELEMENTS::SoSh8p8Type::instance_;

Discret::ELEMENTS::SoSh8p8Type& Discret::ELEMENTS::SoSh8p8Type::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::SoSh8p8Type::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoSh8p8(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoSh8p8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoSh8p8(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoSh8p8Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoSh8p8(id, owner));
  return ele;
}


void Discret::ELEMENTS::SoSh8p8Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
  np = 1;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoSh8p8Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}

void Discret::ELEMENTS::SoSh8p8Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder()
                     .add_int_vector("HEX8", 8)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_named_string("STAB")
                     .add_named_string("ANS")
                     .add_named_string("LIN")
                     .add_named_string("THICKDIR")
                     .add_named_string("EAS")
                     .add_named_string("ISO")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double("STRENGTH")
                     .build();
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 03/09|
 *----------------------------------------------------------------------*/
// 9-Voigt C-index                                         0 1 2  3 4 5  6 7 8
// TODO this notation does not fit our usual definitions
const int Discret::ELEMENTS::SoSh8p8::VOIGT9ROW_INCONSISTENT_[NUMDFGR_] = {
    0, 1, 2, 0, 1, 2, 0, 2, 1};
const int Discret::ELEMENTS::SoSh8p8::VOIGT9COL_INCONSISTENT_[NUMDFGR_] = {
    0, 1, 2, 1, 2, 0, 2, 1, 0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 9-Voigt C-indices    0   3   6   8   1   4   5   7   2
// TODO this notation does not fit our usual definitions
const int Discret::ELEMENTS::SoSh8p8::VOIGT3X3NONSYM_INCONSISTENT_[NUMDFGR_] = {
    0, 3, 6, 8, 1, 4, 5, 7, 2};

// 24 displacement and 8 pressure DOFs into 32 total element DOFs
const int Discret::ELEMENTS::SoSh8p8::DISPTODISPPRES_[NUMDISP_] = {
    0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 29, 30};
const int Discret::ELEMENTS::SoSh8p8::PRESTODISPPRES_[NUMPRES_] = {3, 7, 11, 15, 19, 23, 27, 31};


/*----------------------------------------------------------------------*
 |  ctor (public)                                            bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoSh8p8::SoSh8p8(int id, int owner) : Discret::ELEMENTS::SoSh8(id, owner)
{
  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoSh8p8::SoSh8p8(const Discret::ELEMENTS::SoSh8p8& old)
    : Discret::ELEMENTS::SoSh8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoSh8p8::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoSh8p8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh8p8::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class So_sh8 Element
  Discret::ELEMENTS::SoSh8::pack(data);
  // techniques
  add_to_pack(data, stab_);
  add_to_pack(data, ans_);
  add_to_pack(data, lin_);
  add_to_pack(data, iso_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh8p8::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_sh8 Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Discret::ELEMENTS::SoSh8::unpack(basedata);
  // techniques
  stab_ = static_cast<StabilisationType>(extract_int(position, data));
  ans_ = static_cast<AnsType>(extract_int(position, data));
  lin_ = static_cast<LinearizationType>(extract_int(position, data));
  iso_ = static_cast<IsochoricType>(extract_int(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              bborn 03/09|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh8p8::Print(std::ostream& os) const
{
  os << "So_sh8p8 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      tk 04/09   |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh8p8::sosh8p8_expol(
    Core::LinAlg::Matrix<NUMGPT_, Mat::NUM_STRESS_3D>& stresses, Epetra_MultiVector& expolstresses)
{
  // static variables, that are the same for every element
  static Core::LinAlg::Matrix<NUMNOD_, NUMGPT_> expol;
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

  Core::LinAlg::Matrix<NUMNOD_, Mat::NUM_STRESS_3D> nodalstresses;

  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i = 0; i < NUMNOD_; ++i)
  {
    const int lid = expolstresses.Map().LID(NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = 1.0 / Nodes()[i]->NumElement();
      for (int j = 0; j < Mat::NUM_STRESS_3D; ++j)
        (*(expolstresses(j)))[lid] += nodalstresses(i, j) * invmyadjele;
    }
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
