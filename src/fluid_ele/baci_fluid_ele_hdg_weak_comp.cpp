/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly Compressible fluid element based on the HDG method

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "baci_fluid_ele_hdg_weak_comp.hpp"

#include "baci_fluid_ele_action.hpp"
#include "baci_fluid_ele_factory.hpp"
#include "baci_fluid_ele_interface.hpp"
#include "baci_inpar_fluid.hpp"
#include "baci_io_linedefinition.hpp"

BACI_NAMESPACE_OPEN


// initialize static variable
DRT::ELEMENTS::FluidHDGWeakCompType DRT::ELEMENTS::FluidHDGWeakCompType::instance_;

DRT::ELEMENTS::FluidHDGWeakCompType& DRT::ELEMENTS::FluidHDGWeakCompType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::FluidHDGWeakCompType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::FluidHDGWeakComp* object = new DRT::ELEMENTS::FluidHDGWeakComp(-1, -1);
  object->Unpack(data);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidHDGWeakCompType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDHDGWEAKCOMP")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::FluidHDGWeakComp(id, owner));
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidHDGWeakCompType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::FluidHDGWeakComp(id, owner));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakCompType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakCompType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakCompType ::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_fluid;
  FluidType::SetupElementDefinition(definitions_fluid);

  std::map<std::string, INPUT::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, INPUT::LineDefinition>& defs_hdg = definitions["FLUIDHDGWEAKCOMP"];

  for (const auto& [key, fluid_line_def] : defs_fluid)
  {
    defs_hdg[key] = INPUT::LineDefinition::Builder(fluid_line_def)
                        .AddNamedInt("DEG")
                        .AddOptionalNamedInt("SPC")
                        .Build();
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidHDGWeakComp::FluidHDGWeakComp(int id, int owner)
    : Fluid(id, owner), degree_(1), completepol_(true)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidHDGWeakComp::FluidHDGWeakComp(const DRT::ELEMENTS::FluidHDGWeakComp& old)
    : Fluid(old), degree_(old.degree_), completepol_(old.completepol_)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidHDGWeakComp::Clone() const
{
  DRT::ELEMENTS::FluidHDGWeakComp* newelement = new DRT::ELEMENTS::FluidHDGWeakComp(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakComp::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Fluid::Pack(data);

  int degree = degree_;
  AddtoPack(data, degree);
  degree = completepol_;
  AddtoPack(data, degree);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakComp::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::ExtractfromPack(position, data, basedata);
  Fluid::Unpack(basedata);

  int val = 0;
  ExtractfromPack(position, data, val);
  dsassert(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  ExtractfromPack(position, data, val);
  completepol_ = val;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::FluidHDGWeakComp::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
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



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidHDGWeakComp::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = INPUT::get<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();

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
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->Evaluate(
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
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
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
      dserror("Unknown type of action '%i' for FluidHDGWeakComp", act);
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakComp::Print(std::ostream& os) const
{
  os << "FluidHDGWeakComp ";
  Element::Print(os);
}

BACI_NAMESPACE_CLOSE
