/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/

#include "4C_truss3_scatra.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_lin_elast_1D.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType DRT::ELEMENTS::Truss3ScatraType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType& DRT::ELEMENTS::Truss3ScatraType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::Truss3ScatraType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Truss3Scatra(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRUSS3SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3Scatra(id, owner));
    return ele;
  }
  // return base class
  else
    return DRT::ELEMENTS::Truss3Type::Create(eletype, eledistype, id, owner);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3Scatra(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["TRUSS3SCATRA"];

  // get definitions from standard truss element
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_truss;
  Truss3Type::setup_element_definition(definitions_truss);
  std::map<std::string, INPUT::LineDefinition>& defs_truss = definitions_truss["TRUSS3"];

  // copy definitions of standard truss element to truss element for scalar transport coupling
  defs["LINE2"] =
      INPUT::LineDefinition::Builder(defs_truss["LINE2"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3Scatra::Truss3Scatra(int id, int owner)
    : Truss3(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3Scatra::Truss3Scatra(const DRT::ELEMENTS::Truss3Scatra& old)
    : Truss3(static_cast<Truss3>(old)), impltype_(old.impltype_)
{
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss3Scatra::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Truss3Scatra(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Truss3::Pack(data);
  AddtoPack(data, impltype_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Truss3::Unpack(basedata);

  ExtractfromPack(position, data, impltype_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3Scatra::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read base element
  Truss3::ReadElement(eletype, distype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "ElchDiffCond")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond;
  else if (impltype == "ElchDiffCondMultiScale")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond_multiscale;
  else if (impltype == "ElchElectrode")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode;
  else
    FOUR_C_THROW("Invalid implementation type for Truss3Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::calc_internal_force_stiff_tot_lag(
    const std::map<std::string, std::vector<double>>& ele_state,
    CORE::LINALG::SerialDenseVector& forcevec, CORE::LINALG::SerialDenseMatrix& stiffmat)
{
  // safety check
  if (Material()->MaterialType() != CORE::Materials::m_linelast1D_growth and
      Material()->MaterialType() != CORE::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case CORE::Materials::m_linelast1D:
    {
      Truss3::calc_internal_force_stiff_tot_lag(ele_state, forcevec, stiffmat);
      break;
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      CORE::LINALG::Matrix<6, 1> curr_nodal_coords;
      CORE::LINALG::Matrix<6, 6> dtruss_disp_du;
      CORE::LINALG::Matrix<6, 1> dN_dx;
      CORE::LINALG::Matrix<2, 1> nodal_concentration;
      const int ndof = 6;

      prep_calc_internal_force_stiff_tot_lag_sca_tra(
          ele_state, curr_nodal_coords, dtruss_disp_du, dN_dx, nodal_concentration);

      // get data from input
      const auto* growth_mat = static_cast<const MAT::LinElast1DGrowth*>(Material().get());

      // get Gauss rule
      auto intpoints = CORE::FE::IntegrationPoints1D(gaussrule_);

      // computing forcevec and stiffmat
      forcevec.putScalar(0.0);
      stiffmat.putScalar(0.0);
      for (int gp = 0; gp < intpoints.nquad; ++gp)
      {
        const double dx_dxi = lrefe_ / 2.0;
        const double int_fac = dx_dxi * intpoints.qwgt[gp] * crosssec_;

        // get concentration at Gauss point
        const double c_GP =
            project_scalar_to_gauss_point(intpoints.qxg[gp][0], nodal_concentration);

        // calculate stress
        const double PK2_1D = growth_mat->EvaluatePK2(CurrLength(curr_nodal_coords) / lrefe_, c_GP);
        const double stiffness =
            growth_mat->EvaluateStiffness(CurrLength(curr_nodal_coords) / lrefe_, c_GP);

        // calculate residual (force.vec) and linearisation (stiffmat)
        for (int row = 0; row < ndof; ++row)
        {
          const double def_grad = curr_nodal_coords(row) / lrefe_;
          const double scalar_R = int_fac * def_grad * PK2_1D;
          forcevec(row) += dN_dx(row) * scalar_R;
          for (int col = 0; col < ndof; ++col)
          {
            const double ddef_grad_du = dtruss_disp_du(row, col) / lrefe_;
            const double sign = (col < 3 ? 1.0 : -1.0);
            const double dPK2_1D_du =
                2.0 * stiffness * dCurrLengthdu(curr_nodal_coords, col) / lrefe_ * sign;
            const double first_part = dN_dx(row) * ddef_grad_du * PK2_1D;
            const double second_part = dN_dx(row) * def_grad * dPK2_1D_du;
            stiffmat(row, col) += (first_part + second_part) * int_fac;
          }
        }
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("Material type is not supported");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::CalcGPStresses(
    Teuchos::ParameterList& params, const std::map<std::string, std::vector<double>>& ele_state)
{
  // safety check
  if (Material()->MaterialType() != CORE::Materials::m_linelast1D_growth and
      Material()->MaterialType() != CORE::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case CORE::Materials::m_linelast1D:
    {
      Truss3::CalcGPStresses(params, ele_state);
      break;
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
      INPAR::STR::StressType iostress;
      if (IsParamsInterface())
      {
        stressdata = params_interface().StressDataPtr();
        iostress = params_interface().GetStressOutputType();
      }
      else
      {
        stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
            params, "iostress", INPAR::STR::stress_none);
      }

      const CORE::FE::IntegrationPoints1D intpoints(gaussrule_);

      CORE::LINALG::SerialDenseMatrix stress(intpoints.nquad, 1);

      CORE::LINALG::Matrix<6, 1> curr_nodal_coords;
      CORE::LINALG::Matrix<6, 6> dtruss_disp_du;
      CORE::LINALG::Matrix<6, 1> dN_dx;
      CORE::LINALG::Matrix<2, 1> nodal_concentration;

      prep_calc_internal_force_stiff_tot_lag_sca_tra(
          ele_state, curr_nodal_coords, dtruss_disp_du, dN_dx, nodal_concentration);

      // get data from input
      const auto* growth_mat = static_cast<const MAT::LinElast1DGrowth*>(Material().get());

      const double def_grad = CurrLength(curr_nodal_coords) / lrefe_;
      for (int gp = 0; gp < intpoints.nquad; ++gp)
      {
        // get concentration at Gauss point
        const double c_GP =
            project_scalar_to_gauss_point(intpoints.qxg[gp][0], nodal_concentration);

        const double PK2 = growth_mat->EvaluatePK2(def_grad, c_GP);

        switch (iostress)
        {
          case INPAR::STR::stress_2pk:
          {
            stress(gp, 0) = PK2;
            break;
          }
          case INPAR::STR::stress_cauchy:
          {
            stress(gp, 0) = PK2 * def_grad;
            break;
          }
          case INPAR::STR::stress_none:
            break;
          default:
            FOUR_C_THROW("Requested stress type not available");
            break;
        }
      }
      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, stress);
        data.StartPacking();
        AddtoPack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }
    }
    break;
    default:
    {
      FOUR_C_THROW("Material type is not supported");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Truss3Scatra::project_scalar_to_gauss_point(
    const double xi, const CORE::LINALG::Matrix<2, 1>& c) const
{
  return (c(1) - c(0)) / 2.0 * xi + (c(1) + c(0)) / 2.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::extract_elemental_variables(LocationArray& la,
    const DRT::Discretization& discretization, const Teuchos::ParameterList& params,
    std::map<std::string, std::vector<double>>& ele_state)
{
  // add displacements
  Truss3::extract_elemental_variables(la, discretization, params, ele_state);

  // first: check, if micro state is set; if not -> take macro state
  // get nodal phi from micro state
  std::vector<double> phi_ele;
  if (discretization.NumDofSets() == 3 and discretization.HasState(2, "MicroCon"))
  {
    phi_ele.resize(la[2].lm_.size());
    phi_ele.clear();
    auto phi = discretization.GetState(2, "MicroCon");
    if (phi == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'MicroCon'");
    CORE::FE::ExtractMyValues(*phi, phi_ele, la[2].lm_);
  }
  // get nodal phi from micro state
  else if (discretization.HasState(1, "scalarfield"))
  {
    phi_ele.resize(la[1].lm_.size());
    phi_ele.clear();
    auto phi = discretization.GetState(1, "scalarfield");
    if (phi == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'scalar'");
    CORE::FE::ExtractMyValues(*phi, phi_ele, la[1].lm_);
  }
  else
    FOUR_C_THROW("Cannot find state vector");

  if (ele_state.find("phi") == ele_state.end())
    ele_state.emplace(std::make_pair("phi", phi_ele));
  else
    ele_state["phi"] = phi_ele;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::prep_calc_internal_force_stiff_tot_lag_sca_tra(
    const std::map<std::string, std::vector<double>>& ele_state,
    CORE::LINALG::Matrix<6, 1>& curr_nodal_coords,
    CORE::LINALG::Matrix<6, 6>& dcurr_nodal_coords_du, CORE::LINALG::Matrix<6, 1>& dN_dx,
    CORE::LINALG::Matrix<2, 1>& nodal_concentration)
{
  prep_calc_internal_force_stiff_tot_lag(
      ele_state, curr_nodal_coords, dcurr_nodal_coords_du, dN_dx);

  const std::vector<double>& phi_ele = ele_state.at("phi");

  nodal_concentration(0) = phi_ele[0];
  switch (phi_ele.size())
  {
    case 2:
      nodal_concentration(1) = phi_ele[1];
      break;
    case 4:
      nodal_concentration(1) = phi_ele[2];
      break;
    case 6:
      nodal_concentration(1) = phi_ele[3];
      break;
    default:
      FOUR_C_THROW("Vector has size other than 2,4, or 6. Please use different mapping strategy!");
      break;
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::energy(
    const std::map<std::string, std::vector<double>>& ele_state, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& intenergy)
{
  // safety check
  if (Material()->MaterialType() != CORE::Materials::m_linelast1D_growth and
      Material()->MaterialType() != CORE::Materials::m_linelast1D)
    FOUR_C_THROW("only linear elastic growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case CORE::Materials::m_linelast1D:
    {
      Truss3::energy(ele_state, params, intenergy);
      break;
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      CORE::LINALG::Matrix<6, 1> curr_nodal_coords;
      CORE::LINALG::Matrix<6, 6> dtruss_disp_du;
      CORE::LINALG::Matrix<6, 1> dN_dx;
      CORE::LINALG::Matrix<2, 1> nodal_concentration;

      prep_calc_internal_force_stiff_tot_lag_sca_tra(
          ele_state, curr_nodal_coords, dtruss_disp_du, dN_dx, nodal_concentration);

      // get data from input
      const auto* growth_mat = static_cast<const MAT::LinElast1DGrowth*>(Material().get());

      // get Gauss rule
      auto gauss_points = CORE::FE::IntegrationPoints1D(my_gauss_rule(2, gaussexactintegration));

      // internal energy
      for (int j = 0; j < gauss_points.nquad; ++j)
      {
        const double dx_dxi = lrefe_ / 2.0;
        const double int_fac = dx_dxi * gauss_points.qwgt[j] * crosssec_;

        const double c_GP =
            project_scalar_to_gauss_point(gauss_points.qxg[j][0], nodal_concentration);

        eint_ = growth_mat->evaluate_elastic_energy(CurrLength(curr_nodal_coords) / lrefe_, c_GP) *
                int_fac;
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("Material type is not supported");
      break;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
