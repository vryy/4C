/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly Compressible fluid element based on the HDG method

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "fluid_ele_hdg_weak_comp.H"
#include "fluid_ele_action.H"
#include "fluid_ele_factory.H"
#include "fluid_ele_interface.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret_faces.H"
#include "../drt_fem_general/drt_utils_polynomial.H"


// initialize static variable
DRT::ELEMENTS::FluidHDGWeakCompType DRT::ELEMENTS::FluidHDGWeakCompType::instance_;

DRT::ELEMENTS::FluidHDGWeakCompType& DRT::ELEMENTS::FluidHDGWeakCompType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::FluidHDGWeakCompType::Create(const std::vector<char>& data)
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
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_fluid;
  FluidType::SetupElementDefinition(definitions_fluid);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["FLUIDHDGWEAKCOMP"];

  // 3D
  defs["HEX8"] = defs_fluid["HEX8"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["HEX20"] = defs_fluid["HEX20"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["HEX27"] = defs_fluid["HEX27"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["TET4"] = defs_fluid["TET4"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["TET10"] = defs_fluid["TET10"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["WEDGE6"] = defs_fluid["WEDGE6"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["WEDGE15"] = defs_fluid["WEDGE15"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["PYRAMID5"] = defs_fluid["PYRAMID5"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["NURBS8"] = defs_fluid["NURBS8"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["NURBS27"] = defs_fluid["NURBS27"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");

  // 2D
  defs["QUAD4"] = defs_fluid["QUAD4"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["QUAD8"] = defs_fluid["QUAD8"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["QUAD9"] = defs_fluid["QUAD9"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["TRI3"] = defs_fluid["TRI3"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["TRI6"] = defs_fluid["TRI6"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["NURBS4"] = defs_fluid["NURBS4"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
  defs["NURBS9"] = defs_fluid["NURBS9"].AddNamedInt("DEG").AddOptionalNamedInt("SPC");
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
DRT::ELEMENTS::FluidHDGWeakComp::~FluidHDGWeakComp() {}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGWeakComp::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

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
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
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
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params, "action");

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
    case FLD::set_topopt_parameter:
    case FLD::set_general_adjoint_parameter:
    case FLD::set_adjoint_time_parameter:
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
