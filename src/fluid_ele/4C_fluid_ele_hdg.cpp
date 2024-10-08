/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid element based on the HDG method

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_hdg.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_factory.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_io_linedefinition.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


// initialize static variable
Discret::ELEMENTS::FluidHDGType Discret::ELEMENTS::FluidHDGType::instance_;

Discret::ELEMENTS::FluidHDGType& Discret::ELEMENTS::FluidHDGType::instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::FluidHDGType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::ELEMENTS::FluidHDG* object = new Discret::ELEMENTS::FluidHDG(-1, -1);
  object->unpack(buffer);
  return object;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDHDG")
  {
    return Teuchos::make_rcp<Discret::ELEMENTS::FluidHDG>(id, owner);
  }
  return Teuchos::null;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidHDGType::create(
    const int id, const int owner)
{
  return Teuchos::make_rcp<Discret::ELEMENTS::FluidHDG>(id, owner);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = Core::FE::get_dimension(dwele->shape()) + 1;
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDGType::compute_null_space(
    Core::FE::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  if (Core::FE::DiscretizationFaces* facedis = dynamic_cast<Core::FE::DiscretizationFaces*>(&dis))
  {
    const Epetra_Map* rowmap = dis.dof_row_map();
    const int lrows = rowmap->NumMyElements();
    double* mode[6];
    for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

    const Epetra_Map* frowmap = facedis->face_row_map();
    for (int i = 0; i < frowmap->NumMyElements(); ++i)
    {
      std::vector<int> dofs = facedis->dof(0, facedis->l_row_face(i));
      const unsigned int dim = Core::FE::get_dimension(facedis->l_row_face(i)->shape()) + 1;
      FOUR_C_ASSERT(dofs.size() % dim == 0, "Could not match face dofs");
      const unsigned int ndofs = dofs.size() / dim;
      for (unsigned int i = 0; i < dofs.size(); ++i)
      {
        const unsigned int lid = rowmap->LID(dofs[i]);
        for (unsigned int d = 0; d < dim + 1; ++d) mode[d][lid] = 0.;
        mode[i / ndofs][lid] = 1.;
      }
    }
    const Epetra_Map* erowmap = dis.element_row_map();
    for (int i = 0; i < erowmap->NumMyElements(); ++i)
    {
      std::vector<int> dofs = dis.dof(0, dis.l_row_element(i));
      FOUR_C_ASSERT(dofs.size() == 1, "Expect a single pressure dof per element for fluid HDG");
      const unsigned int lid = rowmap->LID(dofs[0]);
      const unsigned int dim = Core::FE::get_dimension(dis.l_row_element(i)->shape());
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
                        .add_named_int("DEG")
                        .add_optional_named_int("SPC")
                        .build();
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
Core::Elements::Element* Discret::ELEMENTS::FluidHDG::clone() const
{
  Discret::ELEMENTS::FluidHDG* newelement = new Discret::ELEMENTS::FluidHDG(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  Pack data (public)                                kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDG::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Fluid::pack(data);

  int degree = degree_;
  add_to_pack(data, degree);
  degree = completepol_;
  add_to_pack(data, degree);
}



/*----------------------------------------------------------------------*
 |  Unpack data (public)                              kronbichler 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidHDG::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Fluid::unpack(basedata_buffer);

  int val = 0;
  extract_from_pack(buffer, val);
  FOUR_C_ASSERT(val >= 0 && val < 255, "Degree out of range");
  degree_ = val;
  extract_from_pack(buffer, val);
  completepol_ = val;

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}


/*----------------------------------------------------------------------*
 |  Read element from input (public)                  kronbichler 06/14 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::FluidHDG::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  bool success = Fluid::read_element(eletype, distype, container);
  degree_ = container.get<int>("DEG");

  completepol_ = container.get_or<int>("SPC", false);

  return success;
}



/*---------------------------------------------------------------------*
|  evaluate the element (public)                      kronbichler 05/13|
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::FluidHDG::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const auto act = Teuchos::getIntegralValue<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = material();

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
      return Discret::ELEMENTS::FluidFactory::provide_impl(shape(), impltype)
          ->evaluate(
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
      return Discret::ELEMENTS::FluidFactory::provide_impl(shape(), impltype)
          ->evaluate_service(
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
void Discret::ELEMENTS::FluidHDG::print(std::ostream& os) const
{
  os << "FluidHDG ";
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
