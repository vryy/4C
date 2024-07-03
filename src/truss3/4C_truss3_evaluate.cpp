/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element (can be connected to beam3 elements and
adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_lin_elast_1D.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_truss3.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public) cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int Discret::ELEMENTS::Truss3::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  set_params_interface_ptr(params);

  Core::Elements::ActionType act = Core::Elements::none;

  if (IsParamsInterface())
  {
    act = params_interface().get_action_type();
  }
  else  // Todo remove as soon as old structural time integration is gone
  {
    // get the action required
    std::string action = params.get<std::string>("action", "calc_none");
    if (action == "calc_none")
      FOUR_C_THROW("No action supplied");
    else if (action == "calc_struct_linstiff")
      act = Core::Elements::struct_calc_linstiff;
    else if (action == "calc_struct_nlnstiff")
      act = Core::Elements::struct_calc_nlnstiff;
    else if (action == "calc_struct_internalforce")
      act = Core::Elements::struct_calc_internalforce;
    else if (action == "calc_struct_linstiffmass")
      act = Core::Elements::struct_calc_linstiffmass;
    else if (action == "calc_struct_nlnstiffmass")
      act = Core::Elements::struct_calc_nlnstiffmass;
    else if (action == "calc_struct_nlnstifflmass")
      act = Core::Elements::struct_calc_nlnstifflmass;
    else if (action == "calc_struct_stress")
      act = Core::Elements::struct_calc_stress;
    else if (action == "calc_struct_update_istep")
      act = Core::Elements::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = Core::Elements::struct_calc_reset_istep;
    else if (action == "calc_struct_ptcstiff")
      act = Core::Elements::struct_calc_ptcstiff;
    else
    {
      std::cout << action << std::endl;
      FOUR_C_THROW("Unknown type of action for Truss3");
    }
  }

  switch (act)
  {
    case Core::Elements::struct_calc_ptcstiff:
    {
      FOUR_C_THROW("EvaluatePTC not implemented");

      break;
    }
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero
       displacement and residual values*/
    case Core::Elements::struct_calc_linstiff:
    {
      FOUR_C_THROW("linear stiffness matrix called, but not implemented");

      break;
    }
    // calculate internal energy
    case Core::Elements::struct_calc_energy:
    {
      std::map<std::string, std::vector<double>> ele_state;
      extract_elemental_variables(la, discretization, params, ele_state);

      energy(ele_state, params, elevec1);

      break;
    }
    // nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is
    // required
    case Core::Elements::struct_calc_nlnstiffmass:
    case Core::Elements::struct_calc_nlnstifflmass:
    case Core::Elements::struct_calc_nlnstiff:
    case Core::Elements::struct_calc_internalforce:
    {
      std::map<std::string, std::vector<double>> ele_state;
      extract_elemental_variables(la, discretization, params, ele_state);

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Core::Elements::struct_calc_nlnstiffmass)
        nln_stiff_mass(ele_state, &elemat1, &elemat2, &elevec1);
      else if (act == Core::Elements::struct_calc_nlnstifflmass)
      {
        nln_stiff_mass(ele_state, &elemat1, &elemat2, &elevec1);
        lump_mass(&elemat2);
      }
      else if (act == Core::Elements::struct_calc_nlnstiff)
        nln_stiff_mass(ele_state, &elemat1, nullptr, &elevec1);
      else if (act == Core::Elements::struct_calc_internalforce)
        nln_stiff_mass(ele_state, nullptr, nullptr, &elevec1);

      break;
    }
    case Core::Elements::struct_calc_stress:
    {
      std::map<std::string, std::vector<double>> ele_state;
      extract_elemental_variables(la, discretization, params, ele_state);

      CalcGPStresses(params, ele_state);
      break;
    }
    case Core::Elements::struct_calc_update_istep:
    case Core::Elements::struct_calc_reset_istep:
    case Core::Elements::struct_calc_recover:
    case Core::Elements::struct_calc_predict:
    {
      // do nothing here
      break;
    }
    default:
    {
      std::cout << "\ncalled element with action type " << ActionType2String(act);
      FOUR_C_THROW("Unknown type of action for Truss3");
      break;
    }
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public) cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int Discret::ELEMENTS::Truss3::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("This method needs to be modified for bio-polymer networks!");

  return 0;
}

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                             cyron 12/10|
 *--------------------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::energy(const std::map<std::string, std::vector<double>>& ele_state,
    Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector& intenergy)
{
  if (Material()->MaterialType() != Core::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic material supported for truss element");

  const std::vector<double>& disp_ele = ele_state.at("disp");

  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  Core::LinAlg::Matrix<6, 1> xcurr;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  Core::LinAlg::Matrix<3, 1> aux;

  const int ndof = 6;

  // current nodal position (first
  for (int j = 0; j < ndof; ++j) xcurr(j) = x_(j) + disp_ele[j];

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));

  double lcurr = std::sqrt(aux(0) * aux(0) + aux(1) * aux(1) + aux(2) * aux(2));

  // calculate deformation gradient
  const double def_grad = lcurr / lrefe_;

  //  strain for total Lagrange or engineering strain
  const double epsilon = kintype_ == KinematicType::tr3_totlag ? 0.5 * (def_grad * def_grad - 1.0)
                                                               : (lcurr - lrefe_) / lrefe_;

  // W_int = 1/2*PK2*A*lrefe*\epsilon
  const auto* mat = static_cast<const Mat::LinElast1D*>(Material().get());
  const double intenergy_calc = mat->evaluate_elastic_energy(epsilon) * crosssec_ * lrefe_;

  if (IsParamsInterface())  // new structural time integration
    params_interface().add_contribution_to_energy_type(intenergy_calc, Solid::internal_energy);
  else  // old structural time integration
  {
    // check length of elevec1
    if (intenergy.length() < 1) FOUR_C_THROW("The given result vector is too short.");

    intenergy(0) = intenergy_calc;
  }

  eint_ = intenergy_calc;
}

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      tk 11/08|
 *--------------------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::nln_stiff_mass(
    const std::map<std::string, std::vector<double>>& ele_state,
    Core::LinAlg::SerialDenseMatrix* stiffmatrix, Core::LinAlg::SerialDenseMatrix* massmatrix,
    Core::LinAlg::SerialDenseVector* force)
{
  /*
   * It is observed that for mixed problems, such is the case for biopolymer network simulations
   * (), the method "Evaluate" hands in the larger matrices and vectors of size of element described
   * in the input file. For example, if the computational volume contains both Beam and Truss
   * elements. The Evaluate hand into the method a 12x12 matrix. However, for truss element we need
   * only 6x6. Therefore, an appropriate mapping needs to be established to ensure proper assemblies
   * for corresponding DOFs. The algorithm implemented here is valid only for linear elements i.e.
   * element containing two nodes.
   */
  // 6x6 Stiffness Matrix of the Truss
  Core::LinAlg::SerialDenseMatrix DummyStiffMatrix;
  DummyStiffMatrix.shape(6, 6);
  DummyStiffMatrix.scale(0);
  // 6x6 force vector of the Truss
  Core::LinAlg::SerialDenseVector DummyForce;
  DummyForce.size(6);
  DummyForce.scale(0);

  switch (kintype_)
  {
    case KinematicType::tr3_totlag:
      nln_stiff_mass_tot_lag(ele_state, DummyStiffMatrix, massmatrix, DummyForce);
      break;
    case KinematicType::tr3_engstrain:
      nln_stiff_mass_eng_str(ele_state, DummyStiffMatrix, massmatrix, DummyForce);
      break;
    default:
      FOUR_C_THROW("Unknown type kintype_ for Truss3");
      break;
  }

  // Map element level into global 12 by 12 element
  if (force->length() > 12)
    FOUR_C_THROW("Vector is larger than 12. Please use different mapping strategy!");
  else if (force->length() == 6)
  {
    for (int i = 0; i < 6; i++) (*force)(i) += DummyForce(i);
  }
  else if (force->length() == 12)
  {
    for (int i = 0; i < 3; i++)
    {
      (*force)(i) += DummyForce(i);
      (*force)(i + 6) += DummyForce(i + 3);
    }
  }

  if (stiffmatrix != nullptr)
  {
    // Map element level into global 12 by 12 element
    if (stiffmatrix->numRows() > 12)
      FOUR_C_THROW("Matrix is larger than 12. Please use different mapping strategy!");
    else if (stiffmatrix->numRows() == 6)
    {
      for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) (*stiffmatrix)(i, j) += DummyStiffMatrix(i, j);
    }
    else if (stiffmatrix->numRows() == 12)
    {
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          (*stiffmatrix)(i, j) += DummyStiffMatrix(i, j);
          (*stiffmatrix)(i, j + 6) += DummyStiffMatrix(i, j + 3);
          (*stiffmatrix)(i + 6, j + 6) += DummyStiffMatrix(i + 3, j + 3);
          (*stiffmatrix)(i + 6, j) += DummyStiffMatrix(i + 3, j);
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::nln_stiff_mass_tot_lag(
    const std::map<std::string, std::vector<double>>& ele_state,
    Core::LinAlg::SerialDenseMatrix& DummyStiffMatrix, Core::LinAlg::SerialDenseMatrix* massmatrix,
    Core::LinAlg::SerialDenseVector& DummyForce)
{
  // calculate force vector and stiffness matrix
  calc_internal_force_stiff_tot_lag(ele_state, DummyForce, DummyStiffMatrix);

  const int ndof = 6;
  const int ndof_per_node = ndof / 2;

  const double density = Material()->Density();

  // calculating mass matrix
  if (massmatrix != nullptr)
  {
    for (int i = 0; i < ndof_per_node; ++i)
    {
      (*massmatrix)(i, i) = density * lrefe_ * crosssec_ / 3.0;
      (*massmatrix)(i + ndof_per_node, i + ndof_per_node) = density * lrefe_ * crosssec_ / 3.0;
      (*massmatrix)(i, i + ndof_per_node) = density * lrefe_ * crosssec_ / 6.0;
      (*massmatrix)(i + ndof_per_node, i) = density * lrefe_ * crosssec_ / 6.0;
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 | linear stiffness and mass matrix (private) | engineering strain measure, small displacements and
 rotations |
 *-----------------------------------------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::nln_stiff_mass_eng_str(
    const std::map<std::string, std::vector<double>>& ele_state,
    Core::LinAlg::SerialDenseMatrix& DummyStiffMatrix, Core::LinAlg::SerialDenseMatrix* massmatrix,
    Core::LinAlg::SerialDenseVector& DummyForce)
{
  if (Material()->MaterialType() != Core::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic material supported for truss element");

  const std::vector<double>& disp_ele = ele_state.at("disp");

  // auxiliary vector for both internal force and stiffness matrix
  Core::LinAlg::Matrix<6, 1> aux;

  const int ndof = 6;
  const int ndof_per_node = ndof / 2;

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xref
  aux(0) = x_(0) - x_(3);
  aux(1) = x_(1) - x_(4);
  aux(2) = x_(2) - x_(5);
  aux(3) = x_(3) - x_(0);
  aux(4) = x_(4) - x_(1);
  aux(5) = x_(5) - x_(2);

  // resulting force scaled by current length
  const auto* mat = static_cast<const Mat::LinElast1D*>(Material().get());

  // displacement vector
  Core::LinAlg::SerialDenseVector disp(ndof);

  // computing linear stiffness matrix
  for (int i = 0; i < ndof; ++i)
  {
    disp(i) = disp_ele[i];
    for (int j = 0; j < ndof; ++j)
      DummyStiffMatrix(i, j) =
          (mat->EvaluateStiffness() * crosssec_ / (lrefe_ * lrefe_ * lrefe_)) * aux(i) * aux(j);
  }

  // computing internal forces
  Core::LinAlg::multiply(DummyForce, DummyStiffMatrix, disp);

  // calculating mass matrix.
  if (massmatrix != nullptr)
  {
    const double density = Material()->Density();
    for (int i = 0; i < ndof_per_node; ++i)
    {
      (*massmatrix)(i, i) = density * lrefe_ * crosssec_ / 3.0;
      (*massmatrix)(i + ndof_per_node, i + ndof_per_node) = density * lrefe_ * crosssec_ / 3.0;
      (*massmatrix)(i, i + ndof_per_node) = density * lrefe_ * crosssec_ / 6.0;
      (*massmatrix)(i + ndof_per_node, i) = density * lrefe_ * crosssec_ / 6.0;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::prep_calc_internal_force_stiff_tot_lag(
    const std::map<std::string, std::vector<double>>& ele_state,
    Core::LinAlg::Matrix<6, 1>& curr_nodal_coords,
    Core::LinAlg::Matrix<6, 6>& dcurr_nodal_coords_du, Core::LinAlg::Matrix<6, 1>& dN_dx)
{
  const std::vector<double>& disp_ele = ele_state.at("disp");

  const int ndof = 6;
  static Core::LinAlg::Matrix<6, 1> xcurr;
  // current nodal position
  for (int j = 0; j < ndof; ++j) xcurr(j) = x_(j) + disp_ele[j];

  /* current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node)
   * compared to reference configuration; note: in general this is not equal to the values in disp
   * since the latter one refers to a nodal displacement compared to a reference configuration
   * before the first time step whereas the following variable refers to the displacement with
   * respect to a reference configuration which may have been set up at any point of time during the
   * simulation (usually this
   * is only important if an element changes its reference position during simulation)*/

  // current length of truss element (node1 - node2)
  curr_nodal_coords(0) = xcurr(0) - xcurr(3);
  curr_nodal_coords(1) = xcurr(1) - xcurr(4);
  curr_nodal_coords(2) = xcurr(2) - xcurr(5);
  curr_nodal_coords(3) = curr_nodal_coords(0);
  curr_nodal_coords(4) = curr_nodal_coords(1);
  curr_nodal_coords(5) = curr_nodal_coords(2);

  // derivative of current length w.r.t. nodal displacements
  dcurr_nodal_coords_du.put_scalar(0.0);
  dcurr_nodal_coords_du(0, 0) = dcurr_nodal_coords_du(1, 1) = dcurr_nodal_coords_du(2, 2) =
      dcurr_nodal_coords_du(3, 0) = dcurr_nodal_coords_du(4, 1) = dcurr_nodal_coords_du(5, 2) = 1.0;
  dcurr_nodal_coords_du(0, 3) = dcurr_nodal_coords_du(1, 4) = dcurr_nodal_coords_du(2, 5) =
      dcurr_nodal_coords_du(3, 3) = dcurr_nodal_coords_du(4, 4) = dcurr_nodal_coords_du(5, 5) =
          -1.0;

  // spatial derivative of shape functions
  dN_dx(0) = dN_dx(1) = dN_dx(2) = 1.0 / lrefe_;
  dN_dx(3) = dN_dx(4) = dN_dx(5) = -1.0 / lrefe_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::calc_internal_force_stiff_tot_lag(
    const std::map<std::string, std::vector<double>>& ele_state,
    Core::LinAlg::SerialDenseVector& forcevec, Core::LinAlg::SerialDenseMatrix& stiffmat)
{
  // safety check
  if (Material()->MaterialType() != Core::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic material supported for truss element");

  static Core::LinAlg::Matrix<6, 1> truss_disp;
  static Core::LinAlg::Matrix<6, 6> dtruss_disp_du;
  static Core::LinAlg::Matrix<6, 1> dN_dx;

  prep_calc_internal_force_stiff_tot_lag(ele_state, truss_disp, dtruss_disp_du, dN_dx);

  const int ndof = 6;

  // Green-Lagrange strain ( 1D truss: epsilon = 0.5 (l^2 - L^2)/L^2)
  const double epsilon_GL = 0.5 * (CurrLength2(truss_disp) - lrefe_ * lrefe_) / (lrefe_ * lrefe_);

  // 2nd Piola-Kirchhoff stress

  const auto* mat = static_cast<const Mat::LinElast1D*>(Material().get());
  const double PK2 = mat->EvaluatePK2(epsilon_GL);
  double stiffness = mat->EvaluateStiffness();
  // domain integration factor for linear shape functions -> constant strains and stresses ->
  // constant factor
  const double int_fac = crosssec_ * lrefe_;

  for (int row = 0; row < ndof; ++row)
  {
    const double def_grad = truss_disp(row) / lrefe_;

    forcevec(row) = dN_dx(row) * def_grad * PK2 * int_fac;
    for (int col = 0; col < ndof; ++col)
    {
      const double d_epsilon_dtruss_disp = truss_disp(col) / (lrefe_ * lrefe_);
      // derivative of deformation gradient w.r.t. nodal displacement
      const double ddef_grad_du = dtruss_disp_du(row, col) / lrefe_;
      // derivative of 2nd Piola Kirchhoff stress w.r.t. nodal displacement
      const double sign = (col < 3 ? 1.0 : -1.0);

      const double d_epsilon_du = d_epsilon_dtruss_disp * sign;
      const double dPK2_du = stiffness * d_epsilon_du;
      // product rule for derivative of forcevec w.r.t. nodal displacement = stiffmat
      const double first_part = dN_dx(row) * ddef_grad_du * PK2;
      const double second_part = dN_dx(row) * def_grad * dPK2_du;
      stiffmat(row, col) = (first_part + second_part) * int_fac;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::CalcGPStresses(
    Teuchos::ParameterList& params, const std::map<std::string, std::vector<double>>& ele_state)
{
  // safety check
  if (Material()->MaterialType() != Core::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic material supported for truss element");

  Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
  Inpar::Solid::StressType iostress;
  if (IsParamsInterface())
  {
    stressdata = params_interface().stress_data_ptr();
    iostress = params_interface().get_stress_output_type();
  }
  else
  {
    stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
    iostress = Core::UTILS::GetAsEnum<Inpar::Solid::StressType>(
        params, "iostress", Inpar::Solid::stress_none);
  }

  const Core::FE::IntegrationPoints1D intpoints(gaussrule_);

  Core::LinAlg::SerialDenseMatrix stress(intpoints.nquad, 1);

  static Core::LinAlg::Matrix<6, 1> curr_nodal_coords;
  static Core::LinAlg::Matrix<6, 6> dtruss_disp_du;
  static Core::LinAlg::Matrix<6, 1> dN_dx;

  prep_calc_internal_force_stiff_tot_lag(ele_state, curr_nodal_coords, dtruss_disp_du, dN_dx);

  // Green-Lagrange strain ( 1D truss: epsilon = 0.5 (l^2 - L^2)/L^2)
  const double epsilon_GL =
      0.5 * (CurrLength2(curr_nodal_coords) - lrefe_ * lrefe_) / (lrefe_ * lrefe_);

  // 2nd Piola-Kirchhoff stress
  const auto* mat = static_cast<const Mat::LinElast1D*>(Material().get());
  const double PK2 = mat->EvaluatePK2(epsilon_GL);

  for (int gp = 0; gp < intpoints.nquad; ++gp)
  {
    switch (iostress)
    {
      case Inpar::Solid::stress_2pk:
      {
        stress(gp, 0) = PK2;
        break;
      }
      case Inpar::Solid::stress_cauchy:
      {
        const double def_grad = CurrLength(curr_nodal_coords) / lrefe_;
        stress(gp, 0) = PK2 * def_grad;
        break;
      }

      case Inpar::Solid::stress_none:
        break;
      default:
        FOUR_C_THROW("Requested stress type not available");
        break;
    }
  }

  {
    Core::Communication::PackBuffer data;
    add_to_pack(data, stress);
    std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::lump_mass(Core::LinAlg::SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (int c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < (*emass).numRows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::extract_elemental_variables(LocationArray& la,
    const Core::FE::Discretization& discretization, const Teuchos::ParameterList& params,
    std::map<std::string, std::vector<double>>& ele_state)
{
  std::vector<double> disp_ele(la[0].lm_.size());

  auto disp = discretization.GetState("displacement");
  if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
  Core::FE::ExtractMyValues(*disp, disp_ele, la[0].lm_);

  if (ele_state.find("disp") == ele_state.end())
    ele_state.emplace(std::make_pair("disp", disp_ele));
  else
    ele_state["disp"] = disp_ele;
}
FOUR_C_NAMESPACE_CLOSE
