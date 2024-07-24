/*----------------------------------------------------------------------*/
/*! \file
\brief 18-node hexahedral (bi-quadratic linear)
\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_hex18.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::SoHex18Type Discret::ELEMENTS::SoHex18Type::instance_;

Discret::ELEMENTS::SoHex18Type& Discret::ELEMENTS::SoHex18Type::instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::SoHex18Type::create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoHex18(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex18Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoHex18(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex18Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoHex18(id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex18Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoHex18Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::SoHex18Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX18"] = Input::LineDefinition::Builder()
                      .add_int_vector("HEX18", 18)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex18::SoHex18(int id, int owner) : SoBase(id, owner)
{
  invJ_.resize(NUMGPT_SOH18, Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  detJ_.resize(NUMGPT_SOH18, 0.0);
  init_gp();

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::instance()->get_parameter_list();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex18::SoHex18(const Discret::ELEMENTS::SoHex18& old)
    : SoBase(old), detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  // can this size be anything but NUMDIM_SOH27 x NUMDIM_SOH27?
  for (int i = 0; i < (int)invJ_.size(); ++i) invJ_[i] = old.invJ_[i];
  init_gp();

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoHex18::clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoHex18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  SoBase::pack(data);

  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  SoBase::unpack(basedata);

  // detJ_
  extract_from_pack(position, data, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(position, data, invJ_[i]);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::print(std::ostream& os) const
{
  os << "So_hex18 ";
  Element::print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                          seitz 11/14 |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SoHex18::surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
|  get vector of lines (public)                            seitz 11/14 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SoHex18::lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
|  Return names of visualization data (public)             seitz 11/14 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::vis_names(std::map<std::string, int>& names)
{
  solid_material()->vis_names(names);

  return;
}

/*----------------------------------------------------------------------*
|  Return visualization data (public)                      seitz 11/14 |
*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex18::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_material()->vis_data(name, data, NUMGPT_SOH18, this->id());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex18::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");

  set_material(0, Mat::Factory(material_id));

  // set up of materials with GP data (e.g., history variables)

  Teuchos::RCP<Core::Mat::Material> mat = material();

  solid_material()->setup(NUMGPT_SOH18, container);


  // read kinematic flag
  std::string kinem = container.get<std::string>("KINEM");
  if (kinem == "linear")
  {
    // kintype_ = soh8_linear;
    FOUR_C_THROW("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (kinem == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX18 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (solid_material()->uses_extended_update())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

void Discret::ELEMENTS::SoHex18::init_gp()
{
  xsi_.resize(NUMGPT_SOH18, Core::LinAlg::Matrix<NUMDIM_SOH18, 1>(true));
  wgt_.resize(NUMGPT_SOH18, 0.);
  Core::FE::IntPointsAndWeights<NUMDIM_SOH18> intpoints(Core::FE::GaussRule3D::hex_18point);
  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    wgt_.at(gp) = (intpoints.ip().qwgt)[gp];
    for (int idim = 0; idim < NUMDIM_SOH18; idim++)
      xsi_.at(gp)(idim) = (intpoints.ip().qxg)[gp][idim];
  }
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           seitz 11/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex18::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material post_setup() routine has already been called and call it if
  // not
  ensure_material_post_setup(params);

  Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18> elemat1(elemat1_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18> elemat2(elemat2_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_SOH18, 1> elevec1(elevec1_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_SOH18, 1> elevec2(elevec2_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_SOH18, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  Discret::ELEMENTS::SoHex18::ActionType act = SoHex18::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = SoHex18::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = SoHex18::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = SoHex18::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = SoHex18::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = SoHex18::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = SoHex18::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = SoHex18::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = SoHex18::calc_struct_eleload;
  else if (action == "calc_struct_update_istep")
    act = SoHex18::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = SoHex18::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = SoHex18::calc_struct_reset_all;
  else if (action == "calc_struct_recover")
    act = SoHex18::calc_recover;
  else if (action == "calc_struct_predict")
    return 0;
  else
    FOUR_C_THROW("Unknown type of action for So_hex8: %s", action.c_str());

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.get_state("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.get_state("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* matptr = nullptr;
      if (elemat1.is_initialized()) matptr = &elemat1;

      nlnstiffmass(lm, mydisp, myres, matptr, nullptr, &elevec1, nullptr, nullptr, params,
          Inpar::Solid::stress_none, Inpar::Solid::strain_none);
      break;
    }


    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.get_state("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.get_state("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18> myemat(true);

      nlnstiffmass(lm, mydisp, myres, &myemat, nullptr, &elevec1, nullptr, nullptr, params,
          Inpar::Solid::stress_none, Inpar::Solid::strain_none);

      break;
    }

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    case calc_struct_linstiffmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.get_state("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.get_state("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");

      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);

      nlnstiffmass(lm, mydisp, myres, &elemat1, &elemat2, &elevec1, nullptr, nullptr, params,
          Inpar::Solid::stress_none, Inpar::Solid::strain_none);

      if (act == calc_struct_nlnstifflmass) lumpmass(&elemat2);

      break;
    }

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.get_comm().MyPID() == owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.get_state("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.get_state("residual displacement");
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
        if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        Core::FE::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        Core::FE::ExtractMyValues(*res, myres, lm);
        Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D> stress;
        Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D> strain;
        auto iostress = Core::UTILS::GetAsEnum<Inpar::Solid::StressType>(
            params, "iostress", Inpar::Solid::stress_none);
        auto iostrain = Core::UTILS::GetAsEnum<Inpar::Solid::StrainType>(
            params, "iostrain", Inpar::Solid::strain_none);

        nlnstiffmass(lm, mydisp, myres, nullptr, nullptr, nullptr, &stress, &strain, params,
            iostress, iostrain);

        {
          Core::Communication::PackBuffer data;

          add_to_pack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }

        {
          Core::Communication::PackBuffer data;

          add_to_pack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      FOUR_C_THROW("this method is not supposed to evaluate a load, use evaluate_neumann(...)");
      break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      solid_material()->update();
      update();
    }
    break;

    //==================================================================================
    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      solid_material()->reset_step();
    }
    break;

    case calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.get_state("residual displacement");
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      recover(myres);
    }
    break;

    //==================================================================================
    default:
      FOUR_C_THROW("Unknown type of action for So_hex18");
      break;
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 11/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex18::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = &condition.parameters().get<std::vector<int>>("onoff");
  const auto* val = &condition.parameters().get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  const double time = std::invoke(
      [&]()
      {
        if (is_params_interface())
          return str_params_interface().get_total_time();
        else
          return params.get("total time", -1.0);
      });

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOH18)
    FOUR_C_THROW("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOH18; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = &condition.parameters().get<std::vector<int>>("funct");
  Core::LinAlg::Matrix<NUMDIM_SOH18, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOH18; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  /* ============================================================================*/

  // update element geometry
  Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH18; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    // shape function and derivatives
    Core::LinAlg::Matrix<NUMNOD_SOH18, 1> shapefunct;
    Core::FE::shape_function<Core::FE::CellType::hex18>(xsi_[gp], shapefunct);
    Core::LinAlg::Matrix<NUMDIM_SOH18, NUMNOD_SOH18> deriv;
    Core::FE::shape_function_deriv1<Core::FE::CellType::hex18>(xsi_[gp], deriv);

    // compute the Jacobian matrix
    Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> jac;
    jac.multiply(deriv, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.determinant();
    if (detJ == 0.0)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOH18; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOH18; ++nodid)
          xrefegp(dim) += shapefunct(nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = wgt_[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOH18; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0) ? Global::Problem::instance()
                                 ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .evaluate(xrefegp.data(), time, dim)
                           : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SOH18; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOH18 + dim] += shapefunct(nodid) * dim_fac;
        }
      }
    }


  } /* ==================================================== end of Loop over GP */

  return 0;
}  // Discret::ELEMENTS::So_hex18::evaluate_neumann

int Discret::ELEMENTS::SoHex18::init_jacobian_mapping()
{
  Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18> xrefe;
  for (int i = 0; i < NUMNOD_SOH18; ++i)
  {
    xrefe(i, 0) = nodes()[i]->x()[0];
    xrefe(i, 1) = nodes()[i]->x()[1];
    xrefe(i, 2) = nodes()[i]->x()[2];
  }
  invJ_.reserve(NUMGPT_SOH18);
  detJ_.reserve(NUMGPT_SOH18);


  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    // reset
    invJ_[gp].clear();
    detJ_[gp] = 0.;

    Core::LinAlg::Matrix<NUMDIM_SOH18, NUMNOD_SOH18> deriv;
    Core::FE::shape_function_deriv1<Core::FE::CellType::hex18>(xsi_[gp], deriv);

    invJ_[gp].multiply(deriv, xrefe);
    detJ_[gp] = invJ_[gp].invert();
    if (detJ_[gp] < 0.) return 1;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::nlnstiffmass(std::vector<int>& lm,  ///< location matrix
    std::vector<double>& disp,                                       ///< current displacements
    std::vector<double>& residual,                                   ///< current residual displ
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix,   ///< element stiffness matrix
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,    ///< element mass matrix
    Core::LinAlg::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
    Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,           ///< algorithmic parameters e.g. time
    const Inpar::Solid::StressType iostress,  ///< stress output option
    const Inpar::Solid::StrainType iostrain   ///< strain output option
)
{
  Core::LinAlg::Matrix<NUMNOD_SOH18, 3> xrefe(false);  // X, material coord. of element
  Core::LinAlg::Matrix<NUMNOD_SOH18, 3> xcurr(false);  // x, current  coord. of element


  for (int i = 0; i < NUMNOD_SOH18; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOH18 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOH18 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOH18 + 2];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  Core::LinAlg::Matrix<NUMDIM_SOH18, NUMNOD_SOH18> N_XYZ;
  // build deformation gradient wrt to material configuration
  Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> defgrd(false);
  // shape functions and their first derivatives
  Core::LinAlg::Matrix<NUMNOD_SOH18, 1> shapefunct;
  Core::LinAlg::Matrix<NUMDIM_SOH18, NUMNOD_SOH18> deriv;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::FE::shape_function<Core::FE::CellType::hex18>(xsi_[gp], shapefunct);
    Core::FE::shape_function_deriv1<Core::FE::CellType::hex18>(xsi_[gp], deriv);

    // by N_XYZ = J^-1 * N_rst
    N_XYZ.multiply(invJ_[gp], deriv);  // (6.21)
    double detJ = detJ_[gp];           // (6.22)

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.multiply_tt(xcurr, N_XYZ);

    // calcualte total rcg
    Core::LinAlg::Matrix<3, 3> cauchygreen(false);
    cauchygreen.multiply_tn(defgrd, defgrd);
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> glstrain(false);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    // B-operator
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_SOH18> bop(false);
    for (int i = 0; i < NUMNOD_SOH18; ++i)
    {
      bop(0, NUMDIM_SOH18 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NUMDIM_SOH18 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NUMDIM_SOH18 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NUMDIM_SOH18 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NUMDIM_SOH18 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NUMDIM_SOH18 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NUMDIM_SOH18 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NUMDIM_SOH18 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NUMDIM_SOH18 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NUMDIM_SOH18 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NUMDIM_SOH18 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NUMDIM_SOH18 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NUMDIM_SOH18 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NUMDIM_SOH18 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NUMDIM_SOH18 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NUMDIM_SOH18 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NUMDIM_SOH18 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NUMDIM_SOH18 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> cmat(true);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> stress(true);
    solid_material()->evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double detJ_w = detJ * wgt_[gp];
    // update internal force vector
    if (force) force->multiply_tn(detJ_w, bop, stress, 1.);

    // update stiffness matrix
    if (stiffmatrix)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      Core::LinAlg::Matrix<6, NUMDOF_SOH18> cb;
      cb.multiply(cmat, bop);
      stiffmatrix->multiply_tn(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      Core::LinAlg::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.scale(detJ_w);                       // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);             // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < NUMNOD_SOH18; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < NUMNOD_SOH18; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < NUMDIM_SOH18; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************
    }

    if (massmatrix)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      double density = material()->density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOH18; ++inod)
      {
        ifactor = shapefunct(inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH18; ++jnod)
        {
          massfactor = shapefunct(jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH18 * inod + 0, NUMDIM_SOH18 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH18 * inod + 1, NUMDIM_SOH18 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH18 * inod + 2, NUMDIM_SOH18 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
  }    // gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                              seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::lumpmass(Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).num_cols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).num_rows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                               seitz 11/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex18Type::initialize(Core::FE::Discretization& dis)
{
  // here we order the nodes such that we have a positive definite jacobian
  //       maybe the python script generating the hex18 elements would be a better place for this.
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoHex18*>(dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18* failed");
    if (actele->init_jacobian_mapping() == 1) actele->flip_t();
  }
  dis.fill_complete(false, false, false);

  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoHex18*>(dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18* failed");
    if (actele->init_jacobian_mapping() == 1) FOUR_C_THROW("why");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  revert the 3rd parameter direction                      seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18::flip_t()
{
  if (node_ids() == nullptr) FOUR_C_THROW("couldn't get node ids");
  // reorder nodes
  int new_nodeids[NUMNOD_SOH18];
  new_nodeids[0] = node_ids()[9];
  new_nodeids[1] = node_ids()[10];
  new_nodeids[2] = node_ids()[11];
  new_nodeids[3] = node_ids()[12];
  new_nodeids[4] = node_ids()[13];
  new_nodeids[5] = node_ids()[14];
  new_nodeids[6] = node_ids()[15];
  new_nodeids[7] = node_ids()[16];
  new_nodeids[8] = node_ids()[17];

  new_nodeids[9] = node_ids()[0];
  new_nodeids[10] = node_ids()[1];
  new_nodeids[11] = node_ids()[2];
  new_nodeids[12] = node_ids()[3];
  new_nodeids[13] = node_ids()[4];
  new_nodeids[14] = node_ids()[5];
  new_nodeids[15] = node_ids()[6];
  new_nodeids[16] = node_ids()[7];
  new_nodeids[17] = node_ids()[8];

  set_node_ids(NUMNOD_SOH18, new_nodeids);
  return;
}

Core::LinAlg::Matrix<18, 3> Discret::ELEMENTS::SoHex18::node_param_coord()
{
  Core::LinAlg::Matrix<18, 3> coord;
  for (int node = 0; node < NUMNOD_SOH18; ++node)
  {
    Core::LinAlg::Matrix<3, 1> nodeCoord = node_param_coord(node);
    for (int i = 0; i < 3; ++i) coord(node, i) = nodeCoord(i);
  }
  return coord;
}

Core::LinAlg::Matrix<3, 1> Discret::ELEMENTS::SoHex18::node_param_coord(const int node)
{
  Core::LinAlg::Matrix<3, 1> coord;

  switch (node % 9)
  {
    case 0:
    case 3:
    case 7:
      coord(0) = -1.;
      break;
    case 4:
    case 6:
    case 8:
      coord(0) = +0.;
      break;
    case 1:
    case 2:
    case 5:
      coord(0) = +1.;
      break;
  }
  switch (node % 9)
  {
    case 0:
    case 1:
    case 4:
      coord(1) = -1.;
      break;
    case 5:
    case 7:
    case 8:
      coord(1) = +0.;
      break;
    case 2:
    case 3:
    case 6:
      coord(1) = +1.;
      break;
  }

  if (node < 9)
    coord(2) = -1.;
  else
    coord(2) = +1.;

  return coord;
}

FOUR_C_NAMESPACE_CLOSE
