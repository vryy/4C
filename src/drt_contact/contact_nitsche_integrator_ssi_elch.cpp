/*----------------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the ssi contact case including
electrochemistry

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "contact_nitsche_integrator_ssi_elch.H"
#include "contact_nitsche_utils.H"

#include "../drt_mat/electrode.H"

#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../drt_scatra_ele/scatra_ele_parameter_boundary.H"
#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_utils.H"

struct CONTACT::CoIntegratorNitscheSsiElch::ElementDataBundle
{
  MORTAR::MortarElement* element;
  const LINALG::SerialDenseVector* shape_funct;
  const LINALG::SerialDenseMatrix* shape_deriv;
  const std::vector<GEN::pairedvector<int, double>>* d_xi_dd;
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoIntegratorNitscheSsiElch::CoIntegratorNitscheSsiElch(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm)
    : CoIntegratorNitscheSsi(params, eletype, comm)
{
  if (std::abs(theta_) > 1.0e-16) dserror("SSI Elch Contact just implemented Adjoint free ...");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheSsiElch::IntegrateGP_3D(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
    GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
    GEN::pairedvector<int, double>& derivjac, double* normal,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
    GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<GEN::pairedvector<int, double>>& derivsxi,
    std::vector<GEN::pairedvector<int, double>>& derivmxi)
{
  GPTSForces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt, gap,
      deriv_gap, normal, dnmap_unit, sxi, mxi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::GPTSForces(MORTAR::MortarElement& slave_ele,
    MORTAR::MortarElement& master_ele, const LINALG::SerialDenseVector& slave_shape,
    const LINALG::SerialDenseMatrix& slave_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_slave_xi_dd,
    const LINALG::SerialDenseVector& master_shape,
    const LINALG::SerialDenseMatrix& master_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_master_xi_dd, const double jac,
    const GEN::pairedvector<int, double>& d_jac_dd, const double gp_wgt, const double gap,
    const GEN::pairedvector<int, double>& d_gap_dd, const double* gp_normal,
    const std::vector<GEN::pairedvector<int, double>>& d_gp_normal_dd, double* slave_xi,
    double* master_xi)
{
  if (slave_ele.Owner() != Comm_.MyPID()) return;

  static const bool do_fast_checks = true;
  // first rough check
  if (do_fast_checks)
  {
    if ((std::abs(theta_) < 1.0e-16) and
        (gap > std::max(slave_ele.MaxEdgeSize(), master_ele.MaxEdgeSize())))
      return;
  }

  dsassert(dim == Dim(), "dimension inconsistency");

  // calculate normals and derivatives
  const LINALG::Matrix<dim, 1> normal(gp_normal, true);
  LINALG::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<GEN::pairedvector<int, double>> d_slave_normal_dd(0, 0);
  std::vector<GEN::pairedvector<int, double>> d_master_normal_dd(0, 0);
  slave_ele.ComputeUnitNormalAtXi(slave_xi, slave_normal.A());
  master_ele.ComputeUnitNormalAtXi(master_xi, master_normal.A());
  slave_ele.DerivUnitNormalAtXi(slave_xi, d_slave_normal_dd);
  master_ele.DerivUnitNormalAtXi(master_xi, d_master_normal_dd);

  double pen = ppn_;
  double pet = ppt_;
  double nitsche_wgt_slave(0.0), nitsche_wgt_master(0.0);

  CONTACT::UTILS::NitscheWeightsAndScaling(
      slave_ele, master_ele, nit_wgt_, dt_, nitsche_wgt_slave, nitsche_wgt_master, pen, pet);

  double cauchy_nn_weighted_average(0.0);
  GEN::pairedvector<int, double> d_cauchy_nn_weighted_average_dd(
      slave_ele.NumNode() * 3 * 12 + slave_ele.MoData().ParentDisp().size() +
      master_ele.MoData().ParentDisp().size());
  GEN::pairedvector<int, double> d_cauchy_nn_weighted_average_ds(
      slave_ele.MoData().ParentScalarDof().size() + master_ele.MoData().ParentScalarDof().size());

  // evaluate cauchy stress components and derivatives
  SoEleCauchy<dim>(slave_ele, slave_xi, d_slave_xi_dd, gp_wgt, slave_normal, d_slave_normal_dd,
      normal, d_gp_normal_dd, nitsche_wgt_slave, cauchy_nn_weighted_average,
      d_cauchy_nn_weighted_average_dd, d_cauchy_nn_weighted_average_ds);
  SoEleCauchy<dim>(master_ele, master_xi, d_master_xi_dd, gp_wgt, master_normal, d_master_normal_dd,
      normal, d_gp_normal_dd, -nitsche_wgt_master, cauchy_nn_weighted_average,
      d_cauchy_nn_weighted_average_dd, d_cauchy_nn_weighted_average_ds);

  const double cauchy_nn_average_pen_gap = cauchy_nn_weighted_average + pen * gap;
  GEN::pairedvector<int, double> d_cauchy_nn_average_pen_gap_dd(
      d_cauchy_nn_weighted_average_dd.size() + d_gap_dd.size());
  for (const auto& p : d_cauchy_nn_weighted_average_dd)
    d_cauchy_nn_average_pen_gap_dd[p.first] += p.second;
  for (const auto& p : d_gap_dd) d_cauchy_nn_average_pen_gap_dd[p.first] += pen * p.second;

  if (cauchy_nn_average_pen_gap < 0.0)
  {
    // test in normal contact direction
    IntegrateTest<dim>(-1.0, slave_ele, slave_shape, slave_shape_deriv, d_slave_xi_dd, jac,
        d_jac_dd, gp_wgt, cauchy_nn_average_pen_gap, d_cauchy_nn_average_pen_gap_dd,
        d_cauchy_nn_weighted_average_ds, normal, d_gp_normal_dd);
    if (!two_half_pass_)
    {
      IntegrateTest<dim>(+1.0, master_ele, master_shape, master_shape_deriv, d_master_xi_dd, jac,
          d_jac_dd, gp_wgt, cauchy_nn_average_pen_gap, d_cauchy_nn_average_pen_gap_dd,
          d_cauchy_nn_weighted_average_ds, normal, d_gp_normal_dd);
    }

    // integrate the scatra-scatra interface condition
    IntegrateSSIInterfaceCondition<dim>(slave_ele, slave_shape, slave_shape_deriv, d_slave_xi_dd,
        master_ele, master_shape, master_shape_deriv, d_master_xi_dd, jac, d_jac_dd, gp_wgt);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::IntegrateTest(const double fac,
    MORTAR::MortarElement& ele, const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_xi_dd, const double jac,
    const GEN::pairedvector<int, double>& d_jac_dd, const double wgt, const double test_val,
    const GEN::pairedvector<int, double>& d_test_val_dd,
    const GEN::pairedvector<int, double>& d_test_val_ds, const LINALG::Matrix<dim, 1>& normal,
    const std::vector<GEN::pairedvector<int, double>>& d_normal_dd)
{
  if (std::abs(fac) < 1.0e-16) return;

  CONTACT::CoIntegratorNitsche::IntegrateTest<dim>(fac, ele, shape, shape_deriv, d_xi_dd, jac,
      d_jac_dd, wgt, test_val, d_test_val_dd, normal, d_normal_dd);

  for (const auto& d_testval_ds : d_test_val_ds)
  {
    double* row = ele.GetNitscheContainer().Kde(d_testval_ds.first);
    for (int s = 0; s < ele.NumNode(); ++s)
    {
      for (int d = 0; d < dim; ++d)
      {
        row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
                ele.ParentElement()->Shape(), ele.FaceParentNumber(), s) *
                dim +
            d] -= fac * jac * wgt * d_testval_ds.second * normal(d) * shape(s);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::IntegrateSSIInterfaceCondition(
    MORTAR::MortarElement& slave_ele, const LINALG::SerialDenseVector& slave_shape,
    const LINALG::SerialDenseMatrix& slave_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_slave_xi_dd,
    MORTAR::MortarElement& master_ele, const LINALG::SerialDenseVector& master_shape,
    const LINALG::SerialDenseMatrix& master_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_master_xi_dd, const double jac,
    const GEN::pairedvector<int, double>& d_jac_dd, const double wgt)
{
  // do only integrate if there is something to integrate!
  if (slave_ele.MoData().ParentScalarDof().empty()) return;
  if (master_ele.MoData().ParentScalarDof().empty()) dserror("This is not allowed!");

  ElementDataBundle electrode_quantities, electrolyte_quantities;
  bool slave_is_electrode(true);
  AssignElectrodeAndElectrolyteQuantities(slave_ele, slave_shape, slave_shape_deriv, d_slave_xi_dd,
      master_ele, master_shape, master_shape_deriv, d_master_xi_dd, slave_is_electrode,
      electrode_quantities, electrolyte_quantities);

  // get the scatra-scatra interface kinetic model
  const int kinetic_model = GetScaTraEleParameterBoundary()->KineticModel();

  // perform integration dependent on scatra-scatra interface kinetic model
  switch (kinetic_model)
  {
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      // access material of parent element for ELCH simulations
      Teuchos::RCP<const MAT::Electrode> electrode_material =
          Teuchos::rcp_dynamic_cast<const MAT::Electrode>(
              electrode_quantities.element->ParentElement()->Material(1));

      // get the relevant parameter
      const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
      const double frt = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT();
      const double kr = GetScaTraEleParameterBoundary()->ChargeTransferConstant();
      const double alphaa = GetScaTraEleParameterBoundary()->AlphaA();
      const double alphac = GetScaTraEleParameterBoundary()->AlphaC();

      // calculate the electrode side concentration, potential and their derivatives at the current
      // gauss point
      double electrode_conc(0.0), electrode_pot(0.0);
      GEN::pairedvector<int, double> d_electrode_conc_dc(0), d_electrode_pot_dpot(0),
          d_electrode_conc_dd(0), d_electrode_pot_dd(0);
      SetupGpElchProperties<dim>(electrode_quantities, electrode_conc, electrode_pot,
          d_electrode_conc_dc, d_electrode_conc_dd, d_electrode_pot_dpot, d_electrode_pot_dd);

      // calculate the electrolyte side concentration, potential and their derivatives at the
      // current gauss point
      double electrolyte_conc(0.0), electrolyte_pot(0.0);
      GEN::pairedvector<int, double> d_electrolyte_conc_dc(0), d_electrolyte_pot_dpot(0),
          d_electrolyte_conc_dd(0), d_electrolyte_pot_dd(0);
      SetupGpElchProperties<dim>(electrolyte_quantities, electrolyte_conc, electrolyte_pot,
          d_electrolyte_conc_dc, d_electrolyte_conc_dd, d_electrolyte_pot_dpot,
          d_electrolyte_pot_dd);

      // extract saturation value of intercalated lithium concentration from electrode material
      const double cmax = electrode_material->CMax();
      if (cmax < 1.0e-12)
        dserror("Saturation value c_max of intercalated lithium concentration is too small!");

      // equilibrium electric potential difference and derivative w.r.t. concentration at electrode
      // surface
      const double epd =
          electrode_material->ComputeOpenCircuitPotential(electrode_conc, faraday, frt);
      const double d_epd_dc = electrode_material->ComputeFirstDerivOpenCircuitPotentialConc(
          electrode_conc, faraday, frt);

      // Butler-Volmer exchange mass flux density
      const double j0(kinetic_model == INPAR::S2I::kinetics_butlervolmerreduced
                          ? kr
                          : kr * std::pow(electrolyte_conc, alphaa) *
                                std::pow(cmax - electrode_conc, alphaa) *
                                std::pow(electrode_conc, alphac));

      // skip further computation in case equilibrium electric potential difference is outside
      // physically meaningful range
      if (not std::isinf(epd))
      {
        // electrode-electrolyte overpotential at integration point
        const double eta = electrode_pot - electrolyte_pot - epd;

        // exponential Butler-Volmer terms
        const double expterm1 = std::exp(alphaa * frt * eta);
        const double expterm2 = std::exp(-alphac * frt * eta);
        const double expterm = expterm1 - expterm2;

        // safety check
        if (std::abs(expterm) > 1.0e5)
          dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
              expterm);

        // calculate Butler-Volmer mass flux density
        const double j = j0 * expterm;

        // initialize a dummy resistance as the method below requires a resistance which is not
        // relevant in this case
        const double dummyresistance(0.0);
        // define flux linearization terms
        double dj_dc_electrode(0.0), dj_dc_electrolyte(0.0), dj_dpot_electrode(0.0),
            dj_dpot_electrolyte(0.0);
        // calculate flux linearizations
        DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils::
            CalculateButlerVolmerElchLinearizations(kinetic_model, j0, frt, d_epd_dc, alphaa,
                alphac, dummyresistance, expterm1, expterm2, kr, faraday, electrolyte_conc,
                electrode_conc, cmax, dj_dc_electrode, dj_dc_electrolyte, dj_dpot_electrode,
                dj_dpot_electrolyte);

        // initialize derivatives of flux w.r.t. electrochemistry dofs
        GEN::pairedvector<int, double> dj_delch(
            d_electrode_conc_dc.size() + d_electrolyte_conc_dc.size() +
            d_electrode_pot_dpot.size() + d_electrolyte_pot_dpot.size());
        for (const auto& d_electrodeconc_dc : d_electrode_conc_dc)
          dj_delch[d_electrodeconc_dc.first] += dj_dc_electrode * d_electrodeconc_dc.second;
        for (const auto& d_electrolyteconc_dc : d_electrolyte_conc_dc)
          dj_delch[d_electrolyteconc_dc.first] += dj_dc_electrolyte * d_electrolyteconc_dc.second;
        for (const auto& d_electrodepot_dpot : d_electrode_pot_dpot)
          dj_delch[d_electrodepot_dpot.first] += dj_dpot_electrode * d_electrodepot_dpot.second;
        for (const auto& d_electrolytepot_dpot : d_electrolyte_pot_dpot)
          dj_delch[d_electrolytepot_dpot.first] +=
              dj_dpot_electrolyte * d_electrolytepot_dpot.second;

        // initialize derivatives of flux w.r.t. displacements
        GEN::pairedvector<int, double> dj_dd(
            d_electrode_conc_dd.size() + d_electrode_pot_dd.size() + d_electrolyte_conc_dd.size() +
            d_electrolyte_pot_dd.size());

        for (const auto& d_electrodeconc_dd : d_electrode_conc_dd)
          dj_dd[d_electrodeconc_dd.first] += dj_dc_electrode * d_electrodeconc_dd.second;
        for (const auto& d_electrodepot_dd : d_electrode_pot_dd)
          dj_dd[d_electrodepot_dd.first] += dj_dpot_electrode * d_electrodepot_dd.second;
        for (const auto& d_electrolyteconc_dd : d_electrolyte_conc_dd)
          dj_dd[d_electrolyteconc_dd.first] += dj_dc_electrolyte * d_electrolyteconc_dd.second;
        for (const auto& d_electrolytepot_dd : d_electrolyte_pot_dd)
          dj_dd[d_electrolytepot_dd.first] += dj_dpot_electrolyte * d_electrolytepot_dd.second;

        if (!two_half_pass_ or slave_is_electrode)
          IntegrateElchTest<dim>(1.0, electrode_quantities, jac, d_jac_dd, wgt, j, dj_dd, dj_delch);
        if (!two_half_pass_ or !slave_is_electrode)
          IntegrateElchTest<dim>(
              -1.0, electrolyte_quantities, jac, d_jac_dd, wgt, j, dj_dd, dj_delch);
      }

      break;
    }
    default:
    {
      dserror("Evaluation is not implemented for this scatra-scatra interface kinetic model: %i",
          kinetic_model);
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::IntegrateElchTest(double fac,
    ElementDataBundle& ele_data_bundle, double jac, const GEN::pairedvector<int, double>& d_jac_dd,
    double wgt, double test_val, const GEN::pairedvector<int, double>& d_test_val_dd,
    const GEN::pairedvector<int, double>& d_test_val_ds)
{
  MORTAR::MortarElement& ele = *ele_data_bundle.element;
  const LINALG::SerialDenseVector& shape_func = *ele_data_bundle.shape_funct;
  const LINALG::SerialDenseMatrix& shape_deriv = *ele_data_bundle.shape_deriv;
  const std::vector<GEN::pairedvector<int, double>>& d_xi_dd = *ele_data_bundle.d_xi_dd;

  // get time integration factors
  const double time_fac = GetScaTraEleParameterTimInt()->TimeFac();
  const double time_fac_rhs = GetScaTraEleParameterTimInt()->TimeFacRhs();

  // get number of electrons per charge transfer reaction
  const int num_electrons = GetScaTraEleParameterBoundary()->NumElectrons();

  // prepare the RHS integration value
  const double val = fac * jac * wgt * test_val;

  // setup derivative of RHS integration value w.r.t. the displacement dofs
  GEN::pairedvector<int, double> d_val_dd(d_jac_dd.size() + d_test_val_dd.size());
  for (const auto& djac_dd : d_jac_dd)
    d_val_dd[djac_dd.first] += fac * djac_dd.second * wgt * test_val;
  for (const auto& d_testval_dd : d_test_val_dd)
    d_val_dd[d_testval_dd.first] += fac * jac * wgt * d_testval_dd.second;

  for (int s = 0; s < ele.NumNode(); ++s)
  {
    const int slave_parent = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), s);
    const int slave_parent_conc = slave_parent * numdofpernode_;
    const int slave_parent_pot = slave_parent_conc + 1;

    *ele.GetNitscheContainer().RhsE(slave_parent_conc) -= shape_func(s) * val * time_fac_rhs;
    *ele.GetNitscheContainer().RhsE(slave_parent_pot) -=
        num_electrons * shape_func(s) * val * time_fac_rhs;

    for (const auto& d_testval_ds : d_test_val_ds)
    {
      double* row = ele.GetNitscheContainer().Kee(d_testval_ds.first);

      row[slave_parent_conc] += shape_func(s) * fac * jac * wgt * d_testval_ds.second * time_fac;
      row[slave_parent_pot] +=
          num_electrons * shape_func(s) * fac * jac * wgt * d_testval_ds.second * time_fac;
    }

    for (const auto& dval_dd : d_val_dd)
    {
      double* row = ele.GetNitscheContainer().Ked(dval_dd.first);

      row[slave_parent_conc] += shape_func(s) * dval_dd.second * time_fac;
      row[slave_parent_pot] += num_electrons * shape_func(s) * dval_dd.second * time_fac;
    }

    for (int e = 0; e < dim - 1; ++e)
    {
      for (const auto& d_xi_dd_e : d_xi_dd[e])
      {
        double* row = ele.GetNitscheContainer().Ked(d_xi_dd_e.first);

        row[slave_parent_conc] += shape_deriv(s, e) * d_xi_dd_e.second * val * time_fac;
        row[slave_parent_pot] +=
            num_electrons * shape_deriv(s, e) * d_xi_dd_e.second * val * time_fac;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::SetupGpElchProperties(ElementDataBundle& ele_data_bundle,
    double& gp_conc, double& gp_pot, GEN::pairedvector<int, double>& d_conc_dc,
    GEN::pairedvector<int, double>& d_conc_dd, GEN::pairedvector<int, double>& d_pot_dpot,
    GEN::pairedvector<int, double>& d_pot_dd)
{
  MORTAR::MortarElement& ele = *ele_data_bundle.element;
  const LINALG::SerialDenseVector& shape_func = *ele_data_bundle.shape_funct;
  const LINALG::SerialDenseMatrix& shape_deriv = *ele_data_bundle.shape_deriv;
  const std::vector<GEN::pairedvector<int, double>>& d_xi_dd = *ele_data_bundle.d_xi_dd;

  // resize and clear derivative vectors
  d_conc_dc.resize(shape_func.Length());
  d_pot_dpot.resize(shape_func.Length());
  d_conc_dc.clear();
  d_pot_dpot.clear();
  std::size_t deriv_size = 0;
  for (int i = 0; i < dim - 1; ++i) deriv_size += d_xi_dd.at(i).size();
  d_conc_dd.resize(deriv_size);
  d_pot_dd.resize(deriv_size);
  d_conc_dd.clear();
  d_pot_dd.clear();

  // calculate the nodal concentrations, potentials and derivatives w.r.t electrochemistry dofs
  LINALG::SerialDenseVector ele_conc(shape_func.Length());
  LINALG::SerialDenseVector ele_pot(shape_func.Length());
  for (int i = 0; i < ele.NumNode(); ++i)
  {
    const int iparent = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), i);
    const int iparent_conc = iparent * numdofpernode_;
    const int iparent_pot = iparent_conc + 1;

    ele_conc(i) = ele.MoData().ParentScalar().at(iparent_conc);
    ele_pot(i) = ele.MoData().ParentScalar().at(iparent_pot);

    d_conc_dc[ele.MoData().ParentScalarDof().at(iparent_conc)] = shape_func(i);
    d_pot_dpot[ele.MoData().ParentScalarDof().at(iparent_pot)] = shape_func(i);
  }

  // calculate the Gauss point concentration and potential
  gp_conc = shape_func.Dot(ele_conc);
  gp_pot = shape_func.Dot(ele_pot);

  // calculate the nodal concentrations, potentials and derivatives w.r.t displacement dofs
  for (int i = 0; i < dim - 1; ++i)
  {
    for (const auto& d_xi_dd_i : d_xi_dd.at(i))
    {
      double& dc_dd = d_conc_dd[d_xi_dd_i.first];
      double& dpot_dd = d_pot_dd[d_xi_dd_i.first];
      for (int n = 0; n < ele.NumNode(); ++n)
      {
        dc_dd += ele_conc(n) * shape_deriv(n, i) * d_xi_dd_i.second;
        dpot_dd += ele_pot(n) * shape_deriv(n, i) * d_xi_dd_i.second;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheSsiElch::SoEleCauchy(MORTAR::MortarElement& mortar_ele,
    double* gp_coord, const std::vector<GEN::pairedvector<int, double>>& d_gp_coord_dd,
    const double gp_wgt, const LINALG::Matrix<dim, 1>& gp_normal,
    const std::vector<GEN::pairedvector<int, double>>& d_gp_normal_dd,
    const LINALG::Matrix<dim, 1>& test_dir,
    const std::vector<GEN::pairedvector<int, double>>& d_test_dir_dd, const double nitsche_wgt,
    double& cauchy_nt_wgt, GEN::pairedvector<int, double>& d_cauchy_nt_dd,
    GEN::pairedvector<int, double>& d_cauchy_nt_de)
{
  Epetra_SerialDenseMatrix d_sigma_nt_de;

  SoEleCauchyStruct<dim>(mortar_ele, gp_coord, d_gp_coord_dd, gp_wgt, gp_normal, d_gp_normal_dd,
      test_dir, d_test_dir_dd, nitsche_wgt, cauchy_nt_wgt, d_cauchy_nt_dd, &d_sigma_nt_de);

  if (!mortar_ele.MoData().ParentScalar().empty())
  {
    for (int i = 0; i < mortar_ele.ParentElement()->NumNode(); ++i)
      d_cauchy_nt_de[mortar_ele.MoData().ParentScalarDof().at(i * numdofpernode_)] +=
          nitsche_wgt * d_sigma_nt_de(i, 0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheSsiElch::AssignElectrodeAndElectrolyteQuantities(
    MORTAR::MortarElement& slave_ele, const LINALG::SerialDenseVector& slave_shape,
    const LINALG::SerialDenseMatrix& slave_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_slave_xi_dd,
    MORTAR::MortarElement& master_ele, const LINALG::SerialDenseVector& master_shape,
    const LINALG::SerialDenseMatrix& master_shape_deriv,
    const std::vector<GEN::pairedvector<int, double>>& d_master_xi_dd, bool& slave_is_electrode,
    ElementDataBundle& electrode_quantitites, ElementDataBundle& electrolyte_quantities)
{
  Teuchos::RCP<const MAT::Electrode> electrode_material =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(slave_ele.ParentElement()->Material(1));
  if (electrode_material == Teuchos::null)
  {
    slave_is_electrode = false;

    electrode_material =
        Teuchos::rcp_dynamic_cast<const MAT::Electrode>(master_ele.ParentElement()->Material(1));

    // safety check
    if (electrode_material == Teuchos::null)
    {
      dserror(
          "Something went wrong, neither slave nor master side is electrode material. This is a "
          "fatal error!");
    }
  }

  if (slave_is_electrode)
  {
    electrode_quantitites.element = &slave_ele;
    electrode_quantitites.shape_funct = &slave_shape;
    electrode_quantitites.shape_deriv = &slave_shape_deriv;
    electrode_quantitites.d_xi_dd = &d_slave_xi_dd;

    electrolyte_quantities.element = &master_ele;
    electrolyte_quantities.shape_funct = &master_shape;
    electrolyte_quantities.shape_deriv = &master_shape_deriv;
    electrolyte_quantities.d_xi_dd = &d_master_xi_dd;
  }
  else
  {
    electrolyte_quantities.element = &slave_ele;
    electrolyte_quantities.shape_funct = &slave_shape;
    electrolyte_quantities.shape_deriv = &slave_shape_deriv;
    electrolyte_quantities.d_xi_dd = &d_slave_xi_dd;

    electrode_quantitites.element = &master_ele;
    electrode_quantitites.shape_funct = &master_shape;
    electrode_quantitites.shape_deriv = &master_shape_deriv;
    electrode_quantitites.d_xi_dd = &d_master_xi_dd;
  }
}