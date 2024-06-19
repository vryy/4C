/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly Compressible fluid element based on the HDG method

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_fluid_ele_hdg_weak_comp.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_factory.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN


// initialize static variable
Discret::ELEMENTS::FluidHDGWeakCompType Discret::ELEMENTS::FluidHDGWeakCompType::instance_;

Discret::ELEMENTS::FluidHDGWeakCompType& Discret::ELEMENTS::FluidHDGWeakCompType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::FluidHDGWeakCompType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::FluidHDGWeakComp* object = new Discret::ELEMENTS::FluidHDGWeakComp(-1, -1);
  object->unpack(data);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGWeakCompType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDHDGWEAKCOMP")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::FluidHDGWeakComp(id, owner));
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGWeakCompType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::FluidHDGWeakComp(id, owner));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakCompType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakCompType::ComputeNullSpace(
    Core::FE::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakCompType ::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_fluid;
  FluidType::setup_element_definition(definitions_fluid);

  std::map<std::string, Input::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, Input::LineDefinition>& defs_hdg = definitions["FLUIDHDGWEAKCOMP"];

  for (const auto& [key, fluid_line_def] : defs_fluid)
  {
    defs_hdg[key] = Input::LineDefinition::Builder(fluid_line_def)
                        .add_named_int("DEG")
                        .add_optional_named_int("SPC")
                        .build();
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidHDGWeakComp::FluidHDGWeakComp(int id, int owner)
    : Fluid(id, owner), degree_(1), completepol_(true)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidHDGWeakComp::FluidHDGWeakComp(
    const Discret::ELEMENTS::FluidHDGWeakComp& old)
    : Fluid(old), degree_(old.degree_), completepol_(old.completepol_)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::FluidHDGWeakComp::Clone() const
{
  Discret::ELEMENTS::FluidHDGWeakComp* newelement = new Discret::ELEMENTS::FluidHDGWeakComp(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakComp::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  Fluid::pack(data);

  int degree = degree_;
  add_to_pack(data, degree);
  degree = completepol_;
  add_to_pack(data, degree);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakComp::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::extract_from_pack(position, data, basedata);
  Fluid::unpack(basedata);

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
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::FluidHDGWeakComp::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  bool success = Fluid::ReadElement(eletype, distype, linedef);
  int degree;
  linedef->extract_int("DEG", degree);
  degree_ = degree;

  if (linedef->has_named("SPC"))
  {
    linedef->extract_int("SPC", degree);
    completepol_ = degree;
  }
  else
    completepol_ = false;

  return success;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::FluidHDGWeakComp::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = Core::UTILS::GetAsEnum<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = Material();

  // switch between different physical types as used below
  std::string impltype = "hdgweakcomp";

  switch (act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      return Discret::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->evaluate(
              this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    break;

    case FLD::calc_mass_matrix:
    case FLD::calc_fluid_error:
    case FLD::integrate_shape:
    case FLD::interpolate_hdg_to_node:
    case FLD::project_fluid_field:
    case FLD::update_local_solution:
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
      FOUR_C_THROW("Unknown type of action '%i' for FluidHDGWeakComp", act);
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGWeakComp::Print(std::ostream& os) const
{
  os << "FluidHDGWeakComp ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
