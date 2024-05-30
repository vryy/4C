/*----------------------------------------------------------------------*/
/*! \file

\brief pyramid shaped solid element

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_pyramid5.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_pyramid5fbar.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SoPyramid5Type DRT::ELEMENTS::SoPyramid5Type::instance_;

DRT::ELEMENTS::SoPyramid5Type& DRT::ELEMENTS::SoPyramid5Type::Instance() { return instance_; }


CORE::COMM::ParObject* DRT::ELEMENTS::SoPyramid5Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoPyramid5(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoPyramid5Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::SoPyramid5(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoPyramid5Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::SoPyramid5(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoPyramid5Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoPyramid5Type::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::SoPyramid5Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["PYRAMID5"] = INPUT::LineDefinition::Builder()
                         .AddIntVector("PYRAMID5", 5)
                         .AddNamedInt("MAT")
                         .AddNamedString("KINEM")
                         .add_optional_named_double_vector("RAD", 3)
                         .add_optional_named_double_vector("AXI", 3)
                         .add_optional_named_double_vector("CIR", 3)
                         .add_optional_named_double_vector("FIBER1", 3)
                         .add_optional_named_double_vector("FIBER2", 3)
                         .add_optional_named_double_vector("FIBER3", 3)
                         .add_optional_named_double("STRENGTH")
                         .add_optional_named_double("GROWTHTRIG")
                         .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoPyramid5::SoPyramid5(int id, int owner)
    : SoBase(id, owner), pstype_(INPAR::STR::PreStress::none), pstime_(0.0), time_(0.0)
{
  kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  invJ_.resize(NUMGPT_SOP5, CORE::LINALG::Matrix<NUMDIM_SOP5, NUMDIM_SOP5>(true));
  detJ_.resize(NUMGPT_SOP5, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = PRESTRESS::GetType();
    pstime_ = PRESTRESS::GetPrestressTime();

    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }
  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOP5, NUMGPT_SOP5));

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoPyramid5::SoPyramid5(const DRT::ELEMENTS::SoPyramid5& old)
    : SoBase(old),
      kintype_(old.kintype_),
      detJ_(old.detJ_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOP5 x NUMDIM_SOP5?
    invJ_[i] = old.invJ_[i];
  }

  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::SoPyramid5::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoPyramid5(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::SoPyramid5::Shape() const { return CORE::FE::CellType::pyramid5; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoPyramid5::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data, kintype_);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  // Pack prestress_
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    CORE::COMM::ParObject::AddtoPack(data, *prestress_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoPyramid5::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<NUMDIM_SOP5, NUMDIM_SOP5>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // Extract prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
    {
      int numgpt = NUMGPT_SOP5;
      // see whether I am actually a So_pyramid5fbar element
      auto* me = dynamic_cast<DRT::ELEMENTS::SoPyramid5fbar*>(this);
      if (me) numgpt += 1;  // one more history entry for centroid data in pyramid5fbar
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOP5, numgpt));
    }
    prestress_->Unpack(tmpprestress);
    // end
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoPyramid5::Print(std::ostream& os) const
{
  os << "So_pyramid5 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*====================================================================*/
/* 5-node pyramid node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (r,s,t) of nodes
 * of biunit pyramid [-1,1]x[-1,1]x[0,1]
 * 5-node pyramid: node 1,2,3,4,5


 *                /(5)\
 *              / //\\ \
 *            // //  \\ \\
 *          //  //    \\  \\
 *        //   //  t   \\   \\
 *      //    //   |    \\    \\
 *    //     //    |     \\     \\
 *  (4)-----//------------\\-----(3)
 *  ||     //      |       \\     ||
 *  ||    //       |        \\    ||
 *  ||   //        o---------\\---------s
 *  ||  //         |          \\  ||
 *  || //          |           \\ ||
 *  ||//           |            \\||
 *  (1)============|=============(2)
 *                 |
 *                 r
 *
 */
/*====================================================================*/

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::SoPyramid5::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface>(CORE::COMM::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::SoPyramid5::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, CORE::Elements::Element>(
      CORE::COMM::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoPyramid5::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);
  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                                  |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoPyramid5::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (CORE::Elements::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOP5, this->Id());
}

FOUR_C_NAMESPACE_CLOSE
