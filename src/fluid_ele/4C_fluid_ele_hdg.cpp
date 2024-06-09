/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid element based on the HDG method

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_hdg.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_factory.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret_faces.hpp"

FOUR_C_NAMESPACE_OPEN


// initialize static variable
Discret::ELEMENTS::FluidHDGType Discret::ELEMENTS::FluidHDGType::instance_;

Discret::ELEMENTS::FluidHDGType& Discret::ELEMENTS::FluidHDGType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::FluidHDGType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::FluidHDG* object = new Discret::ELEMENTS::FluidHDG(-1, -1);
  object->Unpack(data);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDHDG")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::FluidHDG(id, owner));
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::FluidHDG(id, owner));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = Core::FE::getDimension(dwele->Shape()) + 1;
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGType::ComputeNullSpace(
    Discret::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  if (Discret::DiscretizationFaces* facedis = dynamic_cast<Discret::DiscretizationFaces*>(&dis))
  {
    const Epetra_Map* rowmap = dis.dof_row_map();
    const int lrows = rowmap->NumMyElements();
    double* mode[6];
    for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

    const Epetra_Map* frowmap = facedis->FaceRowMap();
    for (int i = 0; i < frowmap->NumMyElements(); ++i)
    {
      std::vector<int> dofs = facedis->Dof(0, facedis->lRowFace(i));
      const unsigned int dim = Core::FE::getDimension(facedis->lRowFace(i)->Shape()) + 1;
      FOUR_C_ASSERT(dofs.size() % dim == 0, "Could not match face dofs");
      const unsigned int ndofs = dofs.size() / dim;
      for (unsigned int i = 0; i < dofs.size(); ++i)
      {
        const unsigned int lid = rowmap->LID(dofs[i]);
        for (unsigned int d = 0; d < dim + 1; ++d) mode[d][lid] = 0.;
        mode[i / ndofs][lid] = 1.;
      }
    }
    const Epetra_Map* erowmap = dis.ElementRowMap();
    for (int i = 0; i < erowmap->NumMyElements(); ++i)
    {
      std::vector<int> dofs = dis.Dof(0, dis.lRowElement(i));
      FOUR_C_ASSERT(dofs.size() == 1, "Expect a single pressure dof per element for fluid HDG");
      const unsigned int lid = rowmap->LID(dofs[0]);
      const unsigned int dim = Core::FE::getDimension(dis.lRowElement(i)->Shape());
      for (unsigned int d = 0; d < dim; ++d) mode[d][lid] = 0.;
      mode[dim][lid] = 1.;
    }
  }
  else
    FOUR_C_THROW("Faces not initialized");
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGType ::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  // Get the the fluid line definitions and amend them with data for HDG elements
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_fluid;
  FluidType::setup_element_definition(definitions_fluid);

  const std::map<std::string, Input::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, Input::LineDefinition>& defs_hdg = definitions["FLUIDHDG"];

  for (const auto& [key, fluid_line_def] : defs_fluid)
  {
    defs_hdg[key] = Input::LineDefinition::Builder(fluid_line_def)
                        .AddNamedInt("DEG")
                        .AddOptionalNamedInt("SPC")
                        .Build();
  }
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                      kronbichler 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidHDG::FluidHDG(int id, int owner)
    : Fluid(id, owner), degree_(1), completepol_(true)
{
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                 kronbichler 05/13|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidHDG::FluidHDG(const Discret::ELEMENTS::FluidHDG& old)
    : Fluid(old), degree_(old.degree_), completepol_(old.completepol_)
{
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public)  |
 |                                                    kronbichler 05/13 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::FluidHDG::Clone() const
{
  Discret::ELEMENTS::FluidHDG* newelement = new Discret::ELEMENTS::FluidHDG(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  Pack data (public)                                kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDG::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  Fluid::Pack(data);

  int degree = degree_;
  add_to_pack(data, degree);
  degree = completepol_;
  add_to_pack(data, degree);
}



/*----------------------------------------------------------------------*
 |  Unpack data (public)                              kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDG::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::extract_from_pack(position, data, basedata);
  Fluid::Unpack(basedata);

  int val = 0;
  extract_from_pack(position, data, val);
  FOUR_C_ASSERT(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  extract_from_pack(position, data, val);
  completepol_ = val;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 |  Read element from input (public)                  kronbichler 06/14 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::FluidHDG::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  bool success = Fluid::ReadElement(eletype, distype, linedef);
  int degree;
  linedef->ExtractInt("DEG", degree);
  degree_ = degree;

  if (linedef->HaveNamed("SPC"))
  {
    linedef->ExtractInt("SPC", degree);
    completepol_ = degree;
  }
  else
    completepol_ = false;

  return success;
}



/*---------------------------------------------------------------------*
|  evaluate the element (public)                      kronbichler 05/13|
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::FluidHDG::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = Core::UTILS::GetAsEnum<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = Material();

  // switch between different physical types as used below
  std::string impltype = "hdg";

  switch (act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      return Discret::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->Evaluate(
              this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    break;

    case FLD::calc_div_u:
    case FLD::calc_mass_matrix:
    case FLD::calc_fluid_error:
    case FLD::calc_dissipation:
    case FLD::integrate_shape:
    case FLD::calc_divop:
    case FLD::interpolate_hdg_to_node:
    case FLD::interpolate_hdg_for_hit:
    case FLD::project_hdg_force_on_dof_vec_for_hit:
    case FLD::project_hdg_initial_field_for_hit:
    case FLD::project_fluid_field:
    case FLD::calc_pressure_average:
    {
      return Discret::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->EvaluateService(
              this, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }

    case FLD::set_general_fluid_parameter:
    case FLD::set_time_parameter:
    case FLD::set_turbulence_parameter:
    case FLD::set_loma_parameter:
      break;

    default:
      FOUR_C_THROW("Unknown type of action '%i' for FluidHDG", act);
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                        kronbichler 05/13|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDG::Print(std::ostream& os) const
{
  os << "FluidHDG ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
