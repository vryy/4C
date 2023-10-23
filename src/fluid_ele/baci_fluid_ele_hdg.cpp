/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid element based on the HDG method

\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_hdg.H"

#include "baci_fluid_ele_action.H"
#include "baci_fluid_ele_factory.H"
#include "baci_fluid_ele_interface.H"
#include "baci_inpar_fluid.H"
#include "baci_lib_discret_faces.H"
#include "baci_lib_linedefinition.H"


// initialize static variable
DRT::ELEMENTS::FluidHDGType DRT::ELEMENTS::FluidHDGType::instance_;

DRT::ELEMENTS::FluidHDGType& DRT::ELEMENTS::FluidHDGType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::FluidHDGType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::FluidHDG* object = new DRT::ELEMENTS::FluidHDG(-1, -1);
  object->Unpack(data);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidHDGType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDHDG")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::FluidHDG(id, owner));
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidHDGType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::FluidHDG(id, owner));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = CORE::DRT::UTILS::getDimension(dwele->Shape()) + 1;
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  if (DRT::DiscretizationFaces* facedis = dynamic_cast<DRT::DiscretizationFaces*>(&dis))
  {
    const Epetra_Map* rowmap = dis.DofRowMap();
    const int lrows = rowmap->NumMyElements();
    double* mode[6];
    for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

    const Epetra_Map* frowmap = facedis->FaceRowMap();
    for (int i = 0; i < frowmap->NumMyElements(); ++i)
    {
      std::vector<int> dofs = facedis->Dof(0, facedis->lRowFace(i));
      const unsigned int dim = CORE::DRT::UTILS::getDimension(facedis->lRowFace(i)->Shape()) + 1;
      dsassert(dofs.size() % dim == 0, "Could not match face dofs");
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
      dsassert(dofs.size() == 1, "Expect a single pressure dof per element for fluid HDG");
      const unsigned int lid = rowmap->LID(dofs[0]);
      const unsigned int dim = CORE::DRT::UTILS::getDimension(dis.lRowElement(i)->Shape());
      for (unsigned int d = 0; d < dim; ++d) mode[d][lid] = 0.;
      mode[dim][lid] = 1.;
    }
  }
  else
    dserror("Faces not initialized");
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDGType ::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  // Get the the fluid line definitions and amend them with data for HDG elements
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_fluid;
  FluidType::SetupElementDefinition(definitions_fluid);

  const std::map<std::string, DRT::INPUT::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hdg = definitions["FLUIDHDG"];

  for (const auto& [key, fluid_line_def] : defs_fluid)
  {
    defs_hdg[key] = INPUT::LineDefinition::Builder(fluid_line_def)
                        .AddNamedInt("DEG")
                        .AddOptionalNamedInt("SPC")
                        .Build();
  }
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                      kronbichler 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidHDG::FluidHDG(int id, int owner)
    : Fluid(id, owner), degree_(1), completepol_(true)
{
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                 kronbichler 05/13|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidHDG::FluidHDG(const DRT::ELEMENTS::FluidHDG& old)
    : Fluid(old), degree_(old.degree_), completepol_(old.completepol_)
{
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public)  |
 |                                                    kronbichler 05/13 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidHDG::Clone() const
{
  DRT::ELEMENTS::FluidHDG* newelement = new DRT::ELEMENTS::FluidHDG(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  Pack data (public)                                kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDG::Pack(DRT::PackBuffer& data) const
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
 |  Unpack data (public)                              kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDG::Unpack(const std::vector<char>& data)
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
 |  Read element from input (public)                  kronbichler 06/14 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::FluidHDG::ReadElement(
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



/*---------------------------------------------------------------------*
|  evaluate the element (public)                      kronbichler 05/13|
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidHDG::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();

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
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
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
      dserror("Unknown type of action '%i' for FluidHDG", act);
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                        kronbichler 05/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidHDG::Print(std::ostream& os) const
{
  os << "FluidHDG ";
  Element::Print(os);
}
