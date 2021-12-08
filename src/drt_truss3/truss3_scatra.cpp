/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/

#include "truss3_scatra.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType DRT::ELEMENTS::Truss3ScatraType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType& DRT::ELEMENTS::Truss3ScatraType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Truss3ScatraType::Create(const std::vector<char>& data)
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
void DRT::ELEMENTS::Truss3ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["TRUSS3SCATRA"];

  // get definitions from standard truss element
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_truss;
  Truss3Type::SetupElementDefinition(definitions_truss);
  std::map<std::string, DRT::INPUT::LineDefinition>& defs_truss = definitions_truss["TRUSS3"];

  // copy definitions of standard truss element to truss element for scalar transport coupling
  defs["LINE2"] = defs_truss["LINE2"];

  // add scalar transport implementation type
  defs["LINE2"].AddNamedString("TYPE");
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
void DRT::ELEMENTS::Truss3Scatra::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Truss3::Unpack(basedata);

  ExtractfromPack(position, data, impltype_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3Scatra::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
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
    dserror("Invalid implementation type for Truss3Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::CalcInternalForceStiffTotLag(
    const std::map<std::string, std::vector<double>>& ele_state, Epetra_SerialDenseVector& forcevec,
    Epetra_SerialDenseMatrix& stiffmat)
{
  // safety check
  if (Material()->MaterialType() != INPAR::MAT::m_stvenant_growth and
      Material()->MaterialType() != INPAR::MAT::m_stvenant)
    dserror("only St. Venant Kirchhoff growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case INPAR::MAT::m_stvenant:
    {
      Truss3::CalcInternalForceStiffTotLag(ele_state, forcevec, stiffmat);
      break;
    }
    case INPAR::MAT::m_stvenant_growth:
    {
      LINALG::Matrix<6, 1> truss_disp;
      LINALG::Matrix<6, 6> dtruss_disp_du;
      LINALG::Matrix<6, 1> dN_dx;
      LINALG::Matrix<2, 1> nodal_concentration;
      const int ndof = 6;

      PrepCalcInternalForceStiffTotLagScaTra(
          truss_disp, dtruss_disp_du, dN_dx, nodal_concentration, ele_state);

      // get data from input
      const auto* stvk_growth_mat = static_cast<const MAT::StVKGrowth*>(Material().get());
      const double c_0 = stvk_growth_mat->C0();
      const std::vector<double> poly_params = stvk_growth_mat->PolyParams();
      const bool amount_prop_growth = stvk_growth_mat->AmountPropGrowth();
      const double youngs_modulus = stvk_growth_mat->Youngs();

      // get Gauss rule
      auto intpoints = DRT::UTILS::IntegrationPoints1D(gaussrule_);

      // computing forcevec and stiffmat
      forcevec.Scale(0.0);
      stiffmat.Scale(0.0);
      for (int i = 0; i < intpoints.nquad; ++i)
      {
        const double dx_dxi = lrefe_ / 2.0;
        const double int_fac = dx_dxi * intpoints.qwgt[i] * crosssec_;

        // get concentration at Gauss point
        const double c_GP = ProjectScalarToGaussPoint(intpoints.qxg[i][0], nodal_concentration);

        // growth propotional to amount of substance of proportional to concentration
        const double growth_factor =
            amount_prop_growth ? GetGrowthFactorAoSProp(c_GP, c_0, poly_params, truss_disp)
                               : GetGrowthFactorConcProp(c_GP, c_0, poly_params);

        // calculate stress
        const double E_el_1D =
            0.5 * (Lcurr2(truss_disp) / (lrefe_ * lrefe_ * std::pow(growth_factor, 2)) - 1.0);
        const double PK2_1D = 2.0 * youngs_modulus * E_el_1D / growth_factor;

        // calculate residual (force.vec) and linearisation (stiffmat)
        for (int row = 0; row < ndof; ++row)
        {
          const double def_grad = truss_disp(row) / lrefe_;
          const double scalar_R = int_fac * def_grad * PK2_1D;
          forcevec(row) += dN_dx(row) * scalar_R;
          for (int col = 0; col < ndof; ++col)
          {
            const double ddef_grad_du = dtruss_disp_du(row, col) / lrefe_;
            const double sign = (col < 3 ? 1.0 : -1.0);
            const double dPK2_1D_du = 2.0 * youngs_modulus / growth_factor * 1.0 /
                                      (lrefe_ * lrefe_ * std::pow(growth_factor, 2)) * sign *
                                      truss_disp(col);
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
      dserror("Material type is not supported");
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
  if (Material()->MaterialType() != INPAR::MAT::m_stvenant_growth and
      Material()->MaterialType() != INPAR::MAT::m_stvenant)
    dserror("only St. Venant Kirchhoff growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case INPAR::MAT::m_stvenant:
    {
      Truss3::CalcGPStresses(params, ele_state);
      break;
    }
    case INPAR::MAT::m_stvenant_growth:
    {
      Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
      INPAR::STR::StressType iostress;
      if (IsParamsInterface())
      {
        stressdata = ParamsInterface().MutableStressDataPtr();
        iostress = ParamsInterface().GetStressOutputType();
      }
      else
      {
        stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
      }

      const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);

      Epetra_SerialDenseMatrix stress(intpoints.nquad, MAT::NUM_STRESS_3D);

      switch (iostress)
      {
        case INPAR::STR::stress_2pk:
        {
          LINALG::Matrix<6, 1> truss_disp;
          LINALG::Matrix<6, 6> dtruss_disp_du;
          LINALG::Matrix<6, 1> dN_dx;
          LINALG::Matrix<2, 1> nodal_concentration;

          PrepCalcInternalForceStiffTotLagScaTra(
              truss_disp, dtruss_disp_du, dN_dx, nodal_concentration, ele_state);

          // get data from input
          const auto* stvk_growth_mat = static_cast<const MAT::StVKGrowth*>(Material().get());
          const double c_0 = stvk_growth_mat->C0();
          const std::vector<double> poly_params = stvk_growth_mat->PolyParams();
          const bool amount_prop_growth = stvk_growth_mat->AmountPropGrowth();
          const double youngs_modulus = stvk_growth_mat->Youngs();

          for (int i = 0; i < intpoints.nquad; ++i)
          {
            // get concentration at Gauss point
            const double c_GP = ProjectScalarToGaussPoint(intpoints.qxg[i][0], nodal_concentration);

            // growth propotional to amount of substance of proportional to concentration
            const double growth_factor =
                amount_prop_growth ? GetGrowthFactorAoSProp(c_GP, c_0, poly_params, truss_disp)
                                   : GetGrowthFactorConcProp(c_GP, c_0, poly_params);

            // calculate stress
            const double E_el_1D =
                0.5 * (Lcurr2(truss_disp) / (lrefe_ * lrefe_ * std::pow(growth_factor, 2)) - 1.0);
            const double PK2_1D = 2.0 * youngs_modulus * E_el_1D / growth_factor;

            stress(i, 0) = PK2_1D;
          }

          break;
        }
        case INPAR::STR::stress_cauchy:
        {
          dserror("Cauchy stress not supported for truss 3");
          break;
        }

        case INPAR::STR::stress_none:
          break;
        default:
          dserror("Requested stress type not available");
          break;
      }

      {
        DRT::PackBuffer data;
        AddtoPack(data, stress);
        data.StartPacking();
        AddtoPack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }
    }
    break;
    default:
    {
      dserror("Material type is not supported");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Truss3Scatra::ProjectScalarToGaussPoint(
    const double xi, const LINALG::Matrix<2, 1>& c) const
{
  return (c(1) - c(0)) / 2.0 * xi + (c(1) + c(0)) / 2.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Truss3Scatra::GetGrowthFactorConcProp(
    const double c_GP, const double c_0, const std::vector<double>& poly_params) const
{
  return DRT::UTILS::Polynomial(poly_params).Evaluate(c_GP - c_0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
// Get growth factor with Amount of substance instead of concentration
double DRT::ELEMENTS::Truss3Scatra::GetGrowthFactorAoSProp(const double c_GP, const double c_0,
    const std::vector<double>& poly_params, const LINALG::Matrix<6, 1>& truss_disp) const
{
  const double def_grad = Lcurr(truss_disp) / lrefe_;
  return DRT::UTILS::Polynomial(poly_params).Evaluate(c_GP * def_grad - c_0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::ExtractElementalVariables(LocationArray& la,
    const DRT::Discretization& discretization, const Teuchos::ParameterList& params,
    std::map<std::string, std::vector<double>>& ele_state)
{
  // add displacements
  Truss3::ExtractElementalVariables(la, discretization, params, ele_state);

  // first: check, if micro state is set; if not -> take macro state
  // get nodal phi from micro state
  std::vector<double> phi_ele;
  if (discretization.NumDofSets() == 3 and discretization.HasState(2, "MicroCon"))
  {
    phi_ele.resize(la[2].lm_.size());
    phi_ele.clear();
    auto phi = discretization.GetState(2, "MicroCon");
    if (phi == Teuchos::null) dserror("Cannot get state vector 'MicroCon'");
    DRT::UTILS::ExtractMyValues(*phi, phi_ele, la[2].lm_);
  }
  // get nodal phi from micro state
  else if (discretization.HasState(1, "scalarfield"))
  {
    phi_ele.resize(la[1].lm_.size());
    phi_ele.clear();
    auto phi = discretization.GetState(1, "scalarfield");
    if (phi == Teuchos::null) dserror("Cannot get state vectors 'scalar'");
    DRT::UTILS::ExtractMyValues(*phi, phi_ele, la[1].lm_);
  }
  else
    dserror("Cannot find state vector");

  if (ele_state.find("phi") == ele_state.end())
    ele_state.emplace(std::make_pair("phi", phi_ele));
  else
    ele_state["phi"] = phi_ele;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::PrepCalcInternalForceStiffTotLagScaTra(
    LINALG::Matrix<6, 1>& truss_disp, LINALG::Matrix<6, 6>& dtruss_disp_du,
    LINALG::Matrix<6, 1>& dN_dx, LINALG::Matrix<2, 1>& nodal_concentration,
    const std::map<std::string, std::vector<double>>& ele_state)
{
  PrepCalcInternalForceStiffTotLag(ele_state, truss_disp, dtruss_disp_du, dN_dx);

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
      dserror("Vector has size other than 2,4, or 6. Please use different mapping strategy!");
      break;
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3Scatra::Energy(
    const std::map<std::string, std::vector<double>>& ele_state, Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& intenergy)
{
  // safety check
  if (Material()->MaterialType() != INPAR::MAT::m_stvenant_growth and
      Material()->MaterialType() != INPAR::MAT::m_stvenant)
    dserror("only St. Venant Kirchhoff growth material supported for truss element");

  switch (Material()->MaterialType())
  {
    case INPAR::MAT::m_stvenant:
    {
      Truss3::Energy(ele_state, params, intenergy);
      break;
    }
    case INPAR::MAT::m_stvenant_growth:
    {
      LINALG::Matrix<6, 1> truss_disp;
      LINALG::Matrix<6, 6> dtruss_disp_du;
      LINALG::Matrix<6, 1> dN_dx;
      LINALG::Matrix<2, 1> nodal_concentration;

      PrepCalcInternalForceStiffTotLagScaTra(
          truss_disp, dtruss_disp_du, dN_dx, nodal_concentration, ele_state);

      // get data from input
      const auto* stvk_growth_mat = static_cast<const MAT::StVKGrowth*>(Material().get());
      const double c_0 = static_cast<const MAT::StVKGrowth*>(Material().get())->C0();
      const std::vector<double> poly_params = stvk_growth_mat->PolyParams();
      const bool amount_prop_growth = stvk_growth_mat->AmountPropGrowth();
      const double youngs_modulus = stvk_growth_mat->Youngs();

      // get Gauss rule
      auto gauss_points = DRT::UTILS::IntegrationPoints1D(MyGaussRule(2, gaussexactintegration));

      // internal energy
      for (int j = 0; j < gauss_points.nquad; ++j)
      {
        const double dx_dxi = lrefe_ / 2.0;
        const double int_fac = dx_dxi * gauss_points.qwgt[j] * crosssec_;

        const double c_GP = ProjectScalarToGaussPoint(gauss_points.qxg[j][0], nodal_concentration);

        const double growth_factor =
            amount_prop_growth ? GetGrowthFactorAoSProp(c_GP, c_0, poly_params, truss_disp)
                               : GetGrowthFactorConcProp(c_GP, c_0, poly_params);

        const double E_el_1D =
            0.5 * (Lcurr2(truss_disp) / (lrefe_ * lrefe_ * std::pow(growth_factor, 2)) - 1.0);
        const double PK2_1D = 2.0 * youngs_modulus * E_el_1D / growth_factor;
        eint_ = 0.5 * PK2_1D * E_el_1D * int_fac;
      }
      break;
    }
    default:
    {
      dserror("Material type is not supported");
      break;
    }
  }
}