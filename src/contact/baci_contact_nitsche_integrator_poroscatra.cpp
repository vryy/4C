/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poro contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#include "baci_contact_nitsche_integrator_poroscatra.H"

#include "baci_contact_element.H"
#include "baci_contact_integrator.H"
#include "baci_contact_nitsche_utils.H"
#include "baci_contact_node.H"
#include "baci_contact_paramsinterface.H"
#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_mat_elasthyper.H"
#include "baci_mat_structporo.H"
#include "baci_scatra_ele_parameter_boundary.H"
#include "baci_so3_base.H"
#include "baci_so3_hex8.H"
#include "baci_so3_poro.H"

#include <Epetra_FEVector.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::IntegratorNitschePoroscatra::IntegratorNitschePoroscatra(
    Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
    : IntegratorNitschePoro(params, eletype, comm)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitschePoroscatra::SetupGpConcentrations(MORTAR::Element& ele,
    const CORE::LINALG::SerialDenseVector& shape_func,
    const CORE::LINALG::SerialDenseMatrix& shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_xi_dd, double& gp_conc,
    CORE::GEN::pairedvector<int, double>& d_conc_dc,
    CORE::GEN::pairedvector<int, double>& d_conc_dd)
{
  CORE::LINALG::SerialDenseVector ele_conc(shape_func.length());
  for (int i = 0; i < ele.NumNode(); ++i)
    ele_conc(i) = ele.MoData().ParentScalar().at(CORE::FE::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), i));

  // calculate gp concentration
  gp_conc = shape_func.dot(ele_conc);

  // calculate derivative of concentration w.r.t. concentration
  d_conc_dc.resize(shape_func.length());
  d_conc_dc.clear();
  for (int i = 0; i < ele.NumNode(); ++i)
    d_conc_dc[ele.MoData().ParentScalarDof().at(CORE::FE::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), i))] = shape_func(i);

  // calculate derivative of concentration w.r.t. displacements
  std::size_t deriv_size = 0;
  for (int i = 0; i < dim - 1; ++i) deriv_size += d_xi_dd.at(i).size();
  d_conc_dd.resize(deriv_size);
  d_conc_dd.clear();
  for (int i = 0; i < dim - 1; ++i)
  {
    for (const auto& d_xi_dd_i : d_xi_dd.at(i))
    {
      for (int n = 0; n < ele.NumNode(); ++n)
        d_conc_dd[d_xi_dd_i.first] += ele_conc(n) * shape_deriv(n, i) * d_xi_dd_i.second;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitschePoroscatra::IntegrateScaTraTest(const double fac,
    MORTAR::Element& ele, const CORE::LINALG::SerialDenseVector& shape_func,
    const CORE::LINALG::SerialDenseMatrix& shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_xi_dd, const double jac,
    const CORE::GEN::pairedvector<int, double>& d_jac_dd, const double wgt, const double test_val,
    const CORE::GEN::pairedvector<int, double>& d_test_val_dd,
    const CORE::GEN::pairedvector<int, double>& d_test_val_ds,
    const CORE::GEN::pairedvector<int, double>& d_test_val_dp)
{
  // get time integration factors
  // const double time_fac = GetScaTraEleParameterTimInt()->TimeFac();
  // const double time_fac_rhs = GetScaTraEleParameterTimInt()->TimeFacRhs();

  const double time_fac = 1.0;      // TODO parameterize this
  const double time_fac_rhs = 1.0;  // TODO parameterize this

  const double val = fac * jac * wgt * test_val;

  for (int s = 0; s < ele.NumNode(); ++s)
  {
    *(ele.GetNitscheContainer().RhsS(CORE::FE::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), s))) +=
        time_fac_rhs * val * shape_func(s);
  }

  for (const auto& d_testval_dp : d_test_val_dp)
  {
    double* row = ele.GetNitscheContainer().Kpp(d_testval_dp.first);
    for (int s = 0; s < ele.NumNode(); ++s)
    {
      row[CORE::FE::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), s)] -=
          time_fac * fac * jac * wgt * d_testval_dp.second * shape_func(s);
    }
  }

  for (const auto& d_testval_ds : d_test_val_ds)
  {
    double* row = ele.GetNitscheContainer().Kss(d_testval_ds.first);
    for (int s = 0; s < ele.NumNode(); ++s)
    {
      row[CORE::FE::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), s)] -=
          time_fac * fac * jac * wgt * d_testval_ds.second * shape_func(s);
    }
  }

  CORE::GEN::pairedvector<int, double> d_val_dd(d_jac_dd.size() + d_test_val_dd.size());
  for (const auto& djac_dd : d_jac_dd)
    d_val_dd[djac_dd.first] += fac * djac_dd.second * wgt * test_val;
  for (const auto& d_testval_dd : d_test_val_dd)
    d_val_dd[d_testval_dd.first] += fac * jac * wgt * d_testval_dd.second;

  for (const auto& dval_dd : d_val_dd)
  {
    double* row = ele.GetNitscheContainer().Ksd(dval_dd.first);
    for (int s = 0; s < ele.NumNode(); ++s)
      row[CORE::FE::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
          ele.FaceParentNumber(), s)] -= time_fac * dval_dd.second * shape_func(s);
  }

  for (int e = 0; e < dim - 1; ++e)
  {
    for (const auto& d_xi_dd_e : d_xi_dd[e])
    {
      double* row = ele.GetNitscheContainer().Ksd(d_xi_dd_e.first);
      for (int s = 0; s < ele.NumNode(); ++s)
        row[CORE::FE::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
            ele.FaceParentNumber(), s)] -= time_fac * val * shape_deriv(s, e) * d_xi_dd_e.second;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitschePoroscatra::IntegrateScatraInterface(MORTAR::Element& slave_ele,
    const CORE::LINALG::SerialDenseVector& slave_shape,
    const CORE::LINALG::SerialDenseMatrix& slave_shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_slave_xi_dd,
    MORTAR::Element& master_ele, const CORE::LINALG::SerialDenseVector& master_shape,
    const CORE::LINALG::SerialDenseMatrix& master_shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_master_xi_dd, const double jac,
    const CORE::GEN::pairedvector<int, double>& d_jac_dd, const double wgt, const double* gpn,
    std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi,
    double contact_stress, const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_d,
    const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_p)
{
  // do only integrate if there is something to integrate!
  if (slave_ele.MoData().ParentScalarDof().empty()) return;
  if (master_ele.MoData().ParentScalarDof().empty()) dserror("This is not allowed!");

  // prepare the slave and master side gauss point concentrations and derivatives w.r.t. the
  // concentration and the displacement
  double slave_conc(0.0), master_conc(0.0);
  CORE::GEN::pairedvector<int, double> d_slave_conc_dc(0), d_master_conc_dc(0), d_slave_conc_dd(0),
      d_master_conc_dd(0);
  SetupGpConcentrations<dim>(slave_ele, slave_shape, slave_shape_deriv, d_slave_xi_dd, slave_conc,
      d_slave_conc_dc, d_slave_conc_dd);
  SetupGpConcentrations<dim>(master_ele, master_shape, master_shape_deriv, d_master_xi_dd,
      master_conc, d_master_conc_dc, d_master_conc_dd);

  // get the scatra-scatra interface condition kinetic model
  // const int kinetic_model = GetScaTraEleParameterBoundary()->KineticModel();

  const int kinetic_model = INPAR::S2I::kinetics_constperm;  // TODO parameterize this
  const double permeability = 0.1;                           // TODO parameterize this

  double flux;

  CORE::GEN::pairedvector<int, double> dflux_dp(d_slave_conc_dc.size() + d_master_conc_dc.size());
  CORE::GEN::pairedvector<int, double> dflux_dc(d_slave_conc_dc.size() + d_master_conc_dc.size());
  CORE::GEN::pairedvector<int, double> dflux_dd(d_slave_conc_dd.size() + d_master_conc_dd.size());

  // perform integration according to kinetic model
  switch (kinetic_model)
  {
    case INPAR::S2I::kinetics_constperm:
    {
      /* buih: here the exchange terms are computed; only constant permeability works */
      // const double permeability = (*GetScaTraEleParameterBoundary()->Permeabilities())[0];

      // calculate the interface flux
      flux = permeability * (slave_conc - master_conc);

      // evaluate derivatives of flux w.r.t. pressures // TODO
      // do we have interaction between pressure and scatra?
      // for (const auto& p : d_slave_conc_dp) dflux_dp[p.first] += permeability * p.second;
      // for (const auto& p : d_master_conc_dp) dflux_dp[p.first] -= permeability * p.second;
      // assuming no interaction between poro pressure and scatra

      // evaluate derivatives of flux w.r.t. concentrations
      for (const auto& p : d_slave_conc_dc) dflux_dc[p.first] += permeability * p.second;
      for (const auto& p : d_master_conc_dc) dflux_dc[p.first] -= permeability * p.second;

      // evaluate derivatives of flux w.r.t. displacements
      for (const auto& p : d_slave_conc_dd) dflux_dd[p.first] += permeability * p.second;
      for (const auto& p : d_master_conc_dd) dflux_dd[p.first] -= permeability * p.second;

      break;
    }
    case INPAR::S2I::kinetics_linearperm:
    {
      // // obtain the positive value of contact stress
      // contact_stress = -contact_stress;

      // calculate the interface flux
      flux = permeability * contact_stress * (slave_conc - master_conc);

      // evaluate derivatives of flux w.r.t. pressures
      // do we have interaction between pressure and scatra?
      // assuming no interaction between poro pressure and scatra
      for (const auto& p : contact_stress_deriv_p)
        dflux_dp[p.first] += permeability * (slave_conc - master_conc) * p.second;

      // evaluate derivatives of flux w.r.t. concentrations
      for (const auto& p : d_slave_conc_dc)
        dflux_dc[p.first] += permeability * contact_stress * p.second;
      for (const auto& p : d_master_conc_dc)
        dflux_dc[p.first] -= permeability * contact_stress * p.second;

      // for (const auto& p : contact_stress_deriv_c) // TODO
      //   dflux_dc[p.first] += permeability * (slave_conc - master_conc) * p.second;
      // assuming the contact stress is independent with the scatra field

      // evaluate derivatives of flux w.r.t. displacements
      for (const auto& p : d_slave_conc_dd)
        dflux_dd[p.first] += permeability * contact_stress * p.second;
      for (const auto& p : d_master_conc_dd)
        dflux_dd[p.first] -= permeability * contact_stress * p.second;

      for (const auto& p : contact_stress_deriv_d)
        dflux_dd[p.first] += permeability * (slave_conc - master_conc) * p.second;

      break;
    }
    default:
    {
      dserror(
          "Integration can not be performed as kinetic model of scatra-scatra interface condition "
          "is not recognized: %i",
          kinetic_model);

      break;
    }
  }

  // contribute to the residual and linearization
  IntegrateScaTraTest<dim>(-1.0, slave_ele, slave_shape, slave_shape_deriv, d_slave_xi_dd, jac,
      d_jac_dd, wgt, flux, dflux_dd, dflux_dc, dflux_dp);
  if (!two_half_pass_)
  {
    IntegrateScaTraTest<dim>(1.0, master_ele, master_shape, master_shape_deriv, d_master_xi_dd, jac,
        d_jac_dd, wgt, flux, dflux_dd, dflux_dc, dflux_dp);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitschePoroscatra::GPTSForces(MORTAR::Element& sele, MORTAR::Element& mele,
    const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
    const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
    std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi)
{
  BaseType::GPTSForces<dim>(sele, mele, sval, sderiv, dsxi, mval, mderiv, dmxi, jac, jacintcellmap,
      wgt, gap, dgapgp, gpn, dnmap_unit, sxi, mxi);

  /* evaluation of contact-induced heat exchange */

  const CORE::LINALG::Matrix<dim, 1> normal(gpn, true);

  // evaluate slave and master normal
  CORE::LINALG::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<CORE::GEN::pairedvector<int, double>> deriv_slave_normal(0, 0);
  std::vector<CORE::GEN::pairedvector<int, double>> deriv_master_normal(0, 0);
  sele.ComputeUnitNormalAtXi(sxi, slave_normal.A());
  mele.ComputeUnitNormalAtXi(mxi, master_normal.A());
  sele.DerivUnitNormalAtXi(sxi, deriv_slave_normal);
  mele.DerivUnitNormalAtXi(mxi, deriv_master_normal);

  // obtain Nitsche parameters
  double pen = ppn_;
  double pet = ppt_;

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);

  // evaluate cauchy stress components and derivatives
  double cauchy_nn_weighted_average = 0.;
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_d(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_p(
      sele.MoData().ParentPFPres().size() + mele.MoData().ParentPFPres().size());

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, normal, dnmap_unit, ws,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);
  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, normal, dnmap_unit,
      -wm, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);

  // TODO check if jacintcellmap can be used for d_jac_dd

  const double cauchy_nn_average_pen_gap = cauchy_nn_weighted_average + pen * gap;

  if (cauchy_nn_average_pen_gap < 0.0)
  {
    CORE::GEN::pairedvector<int, double> cauchy_nn_average_pen_gap_deriv_d(
        cauchy_nn_weighted_average_deriv_d.size() + dgapgp.size());
    for (const auto& p : cauchy_nn_weighted_average_deriv_d)
      cauchy_nn_average_pen_gap_deriv_d[p.first] += p.second;
    for (const auto& p : dgapgp) cauchy_nn_average_pen_gap_deriv_d[p.first] += pen * p.second;

    CORE::GEN::pairedvector<int, double> cauchy_nn_average_pen_gap_deriv_p(
        cauchy_nn_weighted_average_deriv_p.size());
    for (const auto& p : cauchy_nn_weighted_average_deriv_p)
      cauchy_nn_average_pen_gap_deriv_p[p.first] += p.second;

    IntegrateScatraInterface<dim>(sele, sval, sderiv, dsxi, mele, mval, mderiv, dmxi, jac,
        jacintcellmap, wgt, gpn, dnmap_unit, sxi, mxi, cauchy_nn_average_pen_gap,
        cauchy_nn_average_pen_gap_deriv_d, cauchy_nn_weighted_average_deriv_p);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorNitschePoroscatra::IntegrateGP_3D(MORTAR::Element& sele,
    MORTAR::Element& mele, CORE::LINALG::SerialDenseVector& sval,
    CORE::LINALG::SerialDenseVector& lmval, CORE::LINALG::SerialDenseVector& mval,
    CORE::LINALG::SerialDenseMatrix& sderiv, CORE::LINALG::SerialDenseMatrix& mderiv,
    CORE::LINALG::SerialDenseMatrix& lmderiv,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
    CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
    std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi)
{
  // TEUCHOS_FUNC_TIME_MONITOR("CONTACT::IntegratorNitsche::IntegrateGP_3D");
  // We use the consistent element normal for poro contact!
  // if (nit_normal_==INPAR::CONTACT::NitNor_ele)
  {
    double n[3];
    sele.ComputeUnitNormalAtXi(sxi, n);
    std::vector<CORE::GEN::pairedvector<int, double>> dn(3, sele.NumNode() * 3);
    dynamic_cast<CONTACT::Element&>(sele).DerivUnitNormalAtXi(sxi, dn);

    GPTSForces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
        gap, deriv_gap, n, dn, sxi, mxi);
  }
  //  else if (nit_normal_==INPAR::CONTACT::NitNor_sm)
  //    dserror("Want to use the element normal!");
}

// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
// template <int dim>
// void CONTACT::IntegratorNitschePoro::SoEleCauchy(MORTAR::Element& moEle,
//     double* boundary_gpcoord,
//     std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
//     const CORE::LINALG::Matrix<dim, 1>& normal,
//     std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
//     const CORE::LINALG::Matrix<dim, 1>& direction,
//     std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, const double w,
//     double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_d,
//     CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_p,
//     CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_s)
// {

// }

/* template instantiation */

template void CONTACT::IntegratorNitschePoroscatra::IntegrateScatraInterface<3>(
    MORTAR::Element& slave_ele, const CORE::LINALG::SerialDenseVector& slave_shape,
    const CORE::LINALG::SerialDenseMatrix& slave_shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_slave_xi_dd,
    MORTAR::Element& master_ele, const CORE::LINALG::SerialDenseVector& master_shape,
    const CORE::LINALG::SerialDenseMatrix& master_shape_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& d_master_xi_dd, double jac,
    const CORE::GEN::pairedvector<int, double>& d_jac_dd, double wgt, const double* gpn,
    std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi,
    double contact_stress, const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_d,
    const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_p);

BACI_NAMESPACE_CLOSE
