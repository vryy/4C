/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_nitsche_integrator_tsi.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_nitsche_integrator.hpp"
#include "4C_contact_nitsche_utils.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_so3_base.hpp"
#include "4C_so3_plast_ssn.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorNitscheTsi::integrate_gp_3_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  gpts_forces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
      gap, deriv_gap, normal, dnmap_unit, sxi, mxi);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorNitscheTsi::integrate_gp_2_d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  gpts_forces<2>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
      gap, deriv_gap, normal, dnmap_unit, sxi, mxi);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitscheTsi::gpts_forces(Mortar::Element& sele, Mortar::Element& mele,
    const Core::LinAlg::SerialDenseVector& sval, const Core::LinAlg::SerialDenseMatrix& sderiv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dsxi,
    const Core::LinAlg::SerialDenseVector& mval, const Core::LinAlg::SerialDenseMatrix& mderiv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dmxi, const double jac,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const Core::Gen::Pairedvector<int, double>& dgapgp, const double* gpn,
    std::vector<Core::Gen::Pairedvector<int, double>>& deriv_contact_normal, double* sxi,
    double* mxi)
{
  if (sele.Owner() != Comm_.MyPID()) return;

  if (dim != Dim()) FOUR_C_THROW("dimension inconsistency");

  const Core::Gen::Pairedvector<int, double> empty(0);

  double s_gp_temp, m_gp_temp;
  Core::Gen::Pairedvector<int, double> d_s_gp_temp_dT(0), d_m_gp_temp_dT(0), d_s_gp_temp_dd(0),
      d_m_gp_temp_dd(0);
  setup_gp_temp<dim>(sele, sval, sderiv, dsxi, s_gp_temp, d_s_gp_temp_dT, d_s_gp_temp_dd);
  setup_gp_temp<dim>(mele, mval, mderiv, dmxi, m_gp_temp, d_m_gp_temp_dT, d_m_gp_temp_dd);

  Core::LinAlg::Matrix<dim, 1> xgp;
  for (int n = 0; n < sele.num_node(); ++n)
    for (int d = 0; d < dim; ++d)
      xgp(d) += sval(n) * dynamic_cast<Mortar::Node*>(sele.Nodes()[n])->xspatial()[d];

  if (frtype_ != Inpar::CONTACT::friction_none && dim != 3) FOUR_C_THROW("only 3D friction");
  if (frtype_ != Inpar::CONTACT::friction_none && frtype_ != Inpar::CONTACT::friction_coulomb &&
      frtype_ != Inpar::CONTACT::friction_tresca)
    FOUR_C_THROW("only coulomb or tresca friction");
  if (frtype_ == Inpar::CONTACT::friction_coulomb && frcoeff_ < 0.)
    FOUR_C_THROW("negative coulomb friction coefficient");
  if (frtype_ == Inpar::CONTACT::friction_tresca && frbound_ < 0.)
    FOUR_C_THROW("negative tresca friction bound");

  Core::LinAlg::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<Core::Gen::Pairedvector<int, double>> deriv_slave_normal(0, 0);
  std::vector<Core::Gen::Pairedvector<int, double>> deriv_master_normal(0, 0);
  sele.compute_unit_normal_at_xi(sxi, slave_normal.A());
  mele.compute_unit_normal_at_xi(mxi, master_normal.A());
  sele.DerivUnitNormalAtXi(sxi, deriv_slave_normal);
  mele.DerivUnitNormalAtXi(mxi, deriv_master_normal);

  double pen = ppn_;
  double pet = ppt_;

  const Core::LinAlg::Matrix<dim, 1> contact_normal(gpn, true);
  double cauchy_nn_weighted_average = 0.;
  Core::Gen::Pairedvector<int, double> cauchy_nn_weighted_average_deriv(
      sele.num_node() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  Core::Gen::Pairedvector<int, double> cauchy_nn_weighted_average_deriv_T(
      sele.parent_element()->num_node() + mele.parent_element()->num_node());

  Core::LinAlg::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_normal_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_normal_adjoint_test_slave_T(
      sele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));

  Core::LinAlg::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_normal_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_normal_adjoint_test_master_T(
      mele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);

  // variables for friction (declaration only)
  Core::LinAlg::Matrix<dim, 1> t1, t2;
  std::vector<Core::Gen::Pairedvector<int, double>> dt1, dt2;
  Core::LinAlg::Matrix<dim, 1> relVel;
  std::vector<Core::Gen::Pairedvector<int, double>> relVel_deriv(
      dim, sele.num_node() * dim + mele.num_node() * dim + dsxi[0].size() + dmxi[0].size());
  double vt1, vt2;
  Core::Gen::Pairedvector<int, double> dvt1(0);
  Core::Gen::Pairedvector<int, double> dvt2(0);
  double cauchy_nt1_weighted_average = 0.;
  Core::Gen::Pairedvector<int, double> cauchy_nt1_weighted_average_deriv(
      sele.num_node() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  Core::Gen::Pairedvector<int, double> cauchy_nt1_weighted_average_deriv_T(
      sele.parent_element()->num_node() + mele.parent_element()->num_node());
  Core::LinAlg::SerialDenseVector t1_adjoint_test_slave(sele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t1_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t1_adjoint_test_slave_T(
      sele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  Core::LinAlg::SerialDenseVector t1_adjoint_test_master(mele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t1_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t1_adjoint_test_master_T(
      mele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double cauchy_nt2_weighted_average = 0.;
  Core::Gen::Pairedvector<int, double> cauchy_nt2_weighted_average_deriv(
      sele.num_node() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  Core::Gen::Pairedvector<int, double> cauchy_nt2_weighted_average_deriv_T(
      sele.parent_element()->num_node() + mele.parent_element()->num_node());
  Core::LinAlg::SerialDenseVector t2_adjoint_test_slave(sele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t2_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t2_adjoint_test_slave_T(
      sele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  Core::LinAlg::SerialDenseVector t2_adjoint_test_master(mele.MoData().ParentDof().size());
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t2_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector> deriv_t2_adjoint_test_master_T(
      mele.parent_element()->num_node(), -1,
      Core::LinAlg::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  double sigma_nt1_pen_vt1 = 0;
  double sigma_nt2_pen_vt2 = 0;
  Core::Gen::Pairedvector<int, double> d_sigma_nt1_pen_vt1(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt1_weighted_average_deriv.capacity() + dvt1.capacity(),
      0, 0);
  Core::Gen::Pairedvector<int, double> d_sigma_nt2_pen_vt2(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt2_weighted_average_deriv.capacity() + dvt2.capacity(),
      0, 0);
  Core::Gen::Pairedvector<int, double> d_sigma_nt1_pen_vt1_T(
      sele.parent_element()->num_node() + mele.parent_element()->num_node());
  Core::Gen::Pairedvector<int, double> d_sigma_nt2_pen_vt2_T(
      sele.parent_element()->num_node() + mele.parent_element()->num_node());
  // variables for friction (end)

  so_ele_cauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, contact_normal,
      deriv_contact_normal, ws, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
      cauchy_nn_weighted_average_deriv_T, normal_adjoint_test_slave,
      deriv_normal_adjoint_test_slave, deriv_normal_adjoint_test_slave_T);

  so_ele_cauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, contact_normal,
      deriv_contact_normal, -wm, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
      cauchy_nn_weighted_average_deriv_T, normal_adjoint_test_master,
      deriv_normal_adjoint_test_master, deriv_normal_adjoint_test_master_T);

  const double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  Core::Gen::Pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv.size() + dgapgp.size());
  for (const auto& p : cauchy_nn_weighted_average_deriv) d_snn_av_pen_gap[p.first] += p.second;
  for (const auto& p : dgapgp) d_snn_av_pen_gap[p.first] += pen * p.second;

  // evaluation of tangential stuff
  if (frtype_)
  {
    CONTACT::UTILS::BuildTangentVectors<dim>(
        contact_normal.A(), deriv_contact_normal, t1.A(), dt1, t2.A(), dt2);
    CONTACT::UTILS::RelVelInvariant<dim>(sele, sxi, dsxi, sval, sderiv, mele, mxi, dmxi, mval,
        mderiv, gap, dgapgp, relVel, relVel_deriv);
    CONTACT::UTILS::VectorScalarProduct<dim>(t1, dt1, relVel, relVel_deriv, vt1, dvt1);
    CONTACT::UTILS::VectorScalarProduct<dim>(t2, dt2, relVel, relVel_deriv, vt2, dvt2);

    so_ele_cauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, t1, dt1, ws,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1_adjoint_test_slave, deriv_t1_adjoint_test_slave,
        deriv_t1_adjoint_test_slave_T);
    so_ele_cauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, t1, dt1, -wm,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1_adjoint_test_master, deriv_t1_adjoint_test_master,
        deriv_t1_adjoint_test_master_T);

    so_ele_cauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, t2, dt2, ws,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2_adjoint_test_slave, deriv_t2_adjoint_test_slave,
        deriv_t2_adjoint_test_slave_T);
    so_ele_cauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, t2, dt2, -wm,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2_adjoint_test_master, deriv_t2_adjoint_test_master,
        deriv_t2_adjoint_test_master_T);
  }  // evaluation of tangential stuff


  if (frtype_)
  {
    integrate_test<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1, dt1);
    integrate_test<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2, dt2);
    integrate_test<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1, dt1);
    integrate_test<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2, dt2);

    integrate_adjoint_test<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt1_weighted_average,
        cauchy_nt1_weighted_average_deriv, cauchy_nt1_weighted_average_deriv_T, sele,
        t1_adjoint_test_slave, deriv_t1_adjoint_test_slave, deriv_t1_adjoint_test_slave_T);
    integrate_adjoint_test<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt2_weighted_average,
        cauchy_nt2_weighted_average_deriv, cauchy_nt2_weighted_average_deriv_T, sele,
        t2_adjoint_test_slave, deriv_t2_adjoint_test_slave, deriv_t2_adjoint_test_slave_T);

    integrate_adjoint_test<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt1_weighted_average,
        cauchy_nt1_weighted_average_deriv, cauchy_nt1_weighted_average_deriv_T, mele,
        t1_adjoint_test_master, deriv_t1_adjoint_test_master, deriv_t1_adjoint_test_master_T);
    integrate_adjoint_test<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt2_weighted_average,
        cauchy_nt2_weighted_average_deriv, cauchy_nt2_weighted_average_deriv_T, mele,
        t2_adjoint_test_master, deriv_t2_adjoint_test_master, deriv_t2_adjoint_test_master_T);
  }

  if (snn_av_pen_gap >= 0.)
  {
    integrate_test<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);
    integrate_test<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);

    integrate_adjoint_test<dim>(-theta_ / pen, jac, jacintcellmap, wgt, cauchy_nn_weighted_average,
        cauchy_nn_weighted_average_deriv, cauchy_nn_weighted_average_deriv_T, sele,
        normal_adjoint_test_slave, deriv_normal_adjoint_test_slave,
        deriv_normal_adjoint_test_slave_T);
    integrate_adjoint_test<dim>(-theta_ / pen, jac, jacintcellmap, wgt, cauchy_nn_weighted_average,
        cauchy_nn_weighted_average_deriv, cauchy_nn_weighted_average_deriv_T, mele,
        normal_adjoint_test_master, deriv_normal_adjoint_test_master,
        deriv_normal_adjoint_test_master_T);
  }
  else
  {
    // test in normal contact direction
    integrate_test<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);
    integrate_test<dim>(+1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);

    integrate_test<dim>(-theta_2_ * pen, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, gap,
        dgapgp, empty, contact_normal, deriv_contact_normal);
    integrate_test<dim>(+theta_2_ * pen, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt, gap,
        dgapgp, empty, contact_normal, deriv_contact_normal);

    integrate_adjoint_test<dim>(theta_, jac, jacintcellmap, wgt, gap, dgapgp, empty, sele,
        normal_adjoint_test_slave, deriv_normal_adjoint_test_slave,
        deriv_normal_adjoint_test_slave_T);
    integrate_adjoint_test<dim>(theta_, jac, jacintcellmap, wgt, gap, dgapgp, empty, mele,
        normal_adjoint_test_master, deriv_normal_adjoint_test_master,
        deriv_normal_adjoint_test_master_T);

    if (frtype_)
    {
      double fr = 0.;
      Core::Gen::Pairedvector<int, double> d_fr_d(
          d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
      Core::Gen::Pairedvector<int, double> d_fr_T(cauchy_nn_weighted_average_deriv_T.size() +
                                                  d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
      switch (frtype_)
      {
        case Inpar::CONTACT::friction_coulomb:
        {
          double fr_temp_fac = std::pow((std::max(s_gp_temp, m_gp_temp) - temp_damage_), 2.) /
                               std::pow((temp_damage_ - temp_ref_), 2.);
          fr = frcoeff_ * (-1.) * (snn_av_pen_gap)*fr_temp_fac;
          for (const auto& p : d_snn_av_pen_gap)
            d_fr_d[p.first] += frcoeff_ * (-1.) * p.second * fr_temp_fac;
          for (const auto& p : cauchy_nn_weighted_average_deriv_T)
            d_fr_T[p.first] += frcoeff_ * (-1.) * p.second * fr_temp_fac;
          if (s_gp_temp >= m_gp_temp)
          {
            for (const auto& p : d_s_gp_temp_dd)
            {
              d_fr_d[p.first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                 (s_gp_temp - temp_damage_) /
                                 std::pow((temp_damage_ - temp_ref_), 2.) * p.second;
            }
            for (const auto& p : d_s_gp_temp_dT)
            {
              d_fr_T[p.first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                 (s_gp_temp - temp_damage_) /
                                 std::pow((temp_damage_ - temp_ref_), 2.) * p.second;
            }
          }
          else
          {
            for (const auto& p : d_m_gp_temp_dd)
            {
              d_fr_d[p.first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                 (m_gp_temp - temp_damage_) /
                                 std::pow((temp_damage_ - temp_ref_), 2.) * p.second;
            }
            for (const auto& p : d_m_gp_temp_dT)
            {
              d_fr_T[p.first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                 (m_gp_temp - temp_damage_) /
                                 std::pow((temp_damage_ - temp_ref_), 2.) * p.second;
            }
          }
          break;
        }
        case Inpar::CONTACT::friction_tresca:
          fr = frbound_;
          break;
        default:
          fr = 0.;
          FOUR_C_THROW("why are you here???");
          break;
      }

      double tan_tr = sqrt(
          (cauchy_nt1_weighted_average + pet * vt1) * (cauchy_nt1_weighted_average + pet * vt1) +
          (cauchy_nt2_weighted_average + pet * vt2) * (cauchy_nt2_weighted_average + pet * vt2));

      // stick
      if (tan_tr < fr)
      {
        sigma_nt1_pen_vt1 = cauchy_nt1_weighted_average + pet * vt1;
        for (const auto& p : dvt1) d_sigma_nt1_pen_vt1[p.first] += pet * p.second;
        for (const auto& p : cauchy_nt1_weighted_average_deriv)
          d_sigma_nt1_pen_vt1[p.first] += p.second;
        for (const auto& p : cauchy_nt1_weighted_average_deriv_T)
          d_sigma_nt1_pen_vt1_T[p.first] += p.second;

        sigma_nt2_pen_vt2 = cauchy_nt2_weighted_average + pet * vt2;
        for (const auto& p : dvt2) d_sigma_nt2_pen_vt2[p.first] += pet * p.second;
        for (const auto& p : cauchy_nt2_weighted_average_deriv)
          d_sigma_nt2_pen_vt2[p.first] += p.second;
        for (const auto& p : cauchy_nt2_weighted_average_deriv_T)
          d_sigma_nt2_pen_vt2_T[p.first] += p.second;
      }
      // slip
      else
      {
        Core::Gen::Pairedvector<int, double> tmp_d(
            dgapgp.size() + cauchy_nn_weighted_average_deriv.size() +
                cauchy_nt1_weighted_average_deriv.size() + dvt1.size(),
            0, 0);
        Core::Gen::Pairedvector<int, double> tmp_T(
            sele.parent_element()->num_node() + mele.parent_element()->num_node());
        if (frtype_ == Inpar::CONTACT::friction_coulomb)
        {
          for (const auto& p : d_fr_d) tmp_d[p.first] += p.second / tan_tr;
          for (const auto& p : d_fr_T) tmp_T[p.first] += p.second / tan_tr;
        }
        for (const auto& p : cauchy_nt1_weighted_average_deriv)
          tmp_d[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt1_weighted_average + pet * vt1) * p.second;
        for (const auto& p : dvt1)
          tmp_d[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt1_weighted_average + pet * vt1) * (+pet) * p.second;
        for (const auto& p : cauchy_nt1_weighted_average_deriv_T)
          tmp_T[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt1_weighted_average + pet * vt1) * p.second;

        for (const auto& p : cauchy_nt2_weighted_average_deriv)
          tmp_d[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt2_weighted_average + pet * vt2) * p.second;
        for (const auto& p : dvt2)
          tmp_d[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt2_weighted_average + pet * vt2) * (+pet) * p.second;
        for (const auto& p : cauchy_nt2_weighted_average_deriv_T)
          tmp_T[p.first] += -fr / (tan_tr * tan_tr * tan_tr) *
                            (cauchy_nt2_weighted_average + pet * vt2) * p.second;

        sigma_nt1_pen_vt1 = fr / tan_tr * (cauchy_nt1_weighted_average + pet * vt1);
        for (const auto& p : tmp_d)
          d_sigma_nt1_pen_vt1[p.first] += p.second * (cauchy_nt1_weighted_average + pet * vt1);
        for (const auto& p : cauchy_nt1_weighted_average_deriv)
          d_sigma_nt1_pen_vt1[p.first] += fr / tan_tr * p.second;
        for (const auto& p : dvt1) d_sigma_nt1_pen_vt1[p.first] += fr / tan_tr * pet * p.second;

        for (const auto& p : tmp_T)
          d_sigma_nt1_pen_vt1_T[p.first] += p.second * (cauchy_nt1_weighted_average + pet * vt1);
        for (const auto& p : cauchy_nt1_weighted_average_deriv_T)
          d_sigma_nt1_pen_vt1_T[p.first] += fr / tan_tr * p.second;

        sigma_nt2_pen_vt2 = fr / tan_tr * (cauchy_nt2_weighted_average + pet * vt2);
        for (const auto& p : tmp_d)
          d_sigma_nt2_pen_vt2[p.first] += p.second * (cauchy_nt2_weighted_average + pet * vt2);
        for (const auto& p : cauchy_nt2_weighted_average_deriv)
          d_sigma_nt2_pen_vt2[p.first] += fr / tan_tr * p.second;
        for (const auto& p : dvt2) d_sigma_nt2_pen_vt2[p.first] += fr / tan_tr * pet * p.second;

        for (const auto& p : tmp_T)
          d_sigma_nt2_pen_vt2_T[p.first] += p.second * (cauchy_nt2_weighted_average + pet * vt2);
        for (const auto& p : cauchy_nt2_weighted_average_deriv_T)
          d_sigma_nt2_pen_vt2_T[p.first] += fr / tan_tr * p.second;
      }
      integrate_test<dim>(-theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
          sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, t1, dt1);
      integrate_test<dim>(+theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
          sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, t1, dt1);
      integrate_test<dim>(-theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
          sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, t2, dt2);
      integrate_test<dim>(+theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
          sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, t2, dt2);


      integrate_adjoint_test<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt1_pen_vt1,
          d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, sele, t1_adjoint_test_slave,
          deriv_t1_adjoint_test_slave, deriv_t1_adjoint_test_slave_T);
      integrate_adjoint_test<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt1_pen_vt1,
          d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, mele, t1_adjoint_test_master,
          deriv_t1_adjoint_test_master, deriv_t1_adjoint_test_master_T);
      integrate_adjoint_test<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt2_pen_vt2,
          d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, sele, t2_adjoint_test_slave,
          deriv_t2_adjoint_test_slave, deriv_t2_adjoint_test_slave_T);
      integrate_adjoint_test<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt2_pen_vt2,
          d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, mele, t2_adjoint_test_master,
          deriv_t2_adjoint_test_master, deriv_t2_adjoint_test_master_T);
    }

    // ----------------------------------------------
    // thermo-stuff
    // ----------------------------------------------
    const double beta = gamma_slave_ * gamma_master_ / (gamma_slave_ + gamma_master_);
    const double delta_c = gamma_slave_ / (gamma_slave_ + gamma_master_);
    double diss = 0.;
    Core::Gen::Pairedvector<int, double> d_diss_d(
        d_sigma_nt1_pen_vt1.size() + d_sigma_nt2_pen_vt2.size() + dvt1.size() + dvt2.size());
    Core::Gen::Pairedvector<int, double> d_diss_T(
        d_sigma_nt1_pen_vt1_T.size() + d_sigma_nt2_pen_vt2_T.size());
    if (frtype_)
    {
      diss = (sigma_nt1_pen_vt1 * vt1 + sigma_nt2_pen_vt2 * vt2) / dt_;

      for (const auto& p : d_sigma_nt1_pen_vt1) d_diss_d[p.first] += vt1 * p.second / dt_;
      for (const auto& p : d_sigma_nt2_pen_vt2) d_diss_d[p.first] += vt2 * p.second / dt_;
      for (const auto& p : dvt1) d_diss_d[p.first] += sigma_nt1_pen_vt1 * p.second / dt_;
      for (const auto& p : dvt2) d_diss_d[p.first] += sigma_nt2_pen_vt2 * p.second / dt_;
      for (const auto& p : d_sigma_nt1_pen_vt1_T) d_diss_T[p.first] += vt1 * p.second / dt_;
      for (const auto& p : d_sigma_nt2_pen_vt2_T) d_diss_T[p.first] += vt2 * p.second / dt_;
    }

    switch (nit_thr_)
    {
      case Inpar::CONTACT::NitThr_substitution:
      {
        const double beta_bar = beta * (-snn_av_pen_gap);
        const double q1 = beta_bar * (s_gp_temp - m_gp_temp);

        Core::Gen::Pairedvector<int, double> d_q1_d(
            d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
        for (const auto& p : d_snn_av_pen_gap)
          d_q1_d[p.first] += beta * (-p.second) * (s_gp_temp - m_gp_temp);
        for (const auto& p : d_s_gp_temp_dd) d_q1_d[p.first] += beta_bar * p.second;
        for (const auto& p : d_m_gp_temp_dd) d_q1_d[p.first] += beta_bar * (-p.second);

        Core::Gen::Pairedvector<int, double> d_q1_T(cauchy_nn_weighted_average_deriv_T.size() +
                                                    d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
        for (const auto& p : cauchy_nn_weighted_average_deriv_T)
          d_q1_T[p.first] += beta * (-p.second) * (s_gp_temp - m_gp_temp);
        for (const auto& p : d_s_gp_temp_dT) d_q1_T[p.first] += beta_bar * p.second;
        for (const auto& p : d_m_gp_temp_dT) d_q1_T[p.first] += beta_bar * (-p.second);

        integrate_thermal_test<dim>(
            +1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, q1, d_q1_d, d_q1_T);
        integrate_thermal_test<dim>(
            -1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt, q1, d_q1_d, d_q1_T);

        if (frtype_)  // account for frictional contact
        {
          integrate_thermal_test<dim>(-delta_c, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              diss, d_diss_d, d_diss_T);
          integrate_thermal_test<dim>(-(1 - delta_c), mele, mval, mderiv, dmxi, jac, jacintcellmap,
              wgt, diss, d_diss_d, d_diss_T);
        }
        break;
      }
      case Inpar::CONTACT::NitThr_nitsche:
      {
        double pen_thermo = pp_thermo_;
        double ws_thermo = 0.;
        double wm_thermo = 0.;

        switch (nit_wgt_)
        {
          case Inpar::CONTACT::NitWgt_slave:
            ws_thermo = 1.;
            wm_thermo = 0.;
            pen_thermo /= dynamic_cast<CONTACT::Element&>(sele).TraceHCond();
            break;
          case Inpar::CONTACT::NitWgt_master:
            ws_thermo = 0.;
            wm_thermo = 1.;
            pen_thermo /= dynamic_cast<CONTACT::Element&>(mele).TraceHCond();
            break;
          case Inpar::CONTACT::NitWgt_harmonic:
            ws_thermo = 1. / dynamic_cast<CONTACT::Element&>(mele).TraceHCond();
            ws_thermo /= (ws_thermo + wm_thermo);
            wm_thermo = 1. - ws_thermo;
            pen_thermo =
                ws_thermo * pen_thermo / dynamic_cast<CONTACT::Element&>(sele).TraceHCond() +
                wm_thermo * pen_thermo / dynamic_cast<CONTACT::Element&>(mele).TraceHCond();
            break;
          case Inpar::CONTACT::NitWgt_phyiscal:
            ws_thermo = 1. - delta_c;
            wm_thermo = delta_c;
            pen_thermo =
                ws_thermo * pen_thermo / dynamic_cast<CONTACT::Element&>(sele).TraceHCond() +
                wm_thermo * pen_thermo / dynamic_cast<CONTACT::Element&>(mele).TraceHCond();
            break;
          default:
            FOUR_C_THROW("unknown Nitsche weighting");
            break;
        }

        double qn_weighted_average = 0.;
        Core::Gen::Pairedvector<int, double> deriv_qn_weighted_average_d(
            sele.parent_element()->num_node() * dim + mele.parent_element()->num_node() * dim +
            deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
        Core::Gen::Pairedvector<int, double> deriv_qn_weighted_average_T(
            sele.parent_element()->num_node() + mele.parent_element()->num_node());
        Core::LinAlg::SerialDenseVector thermo_adjoint_test_slave(
            sele.parent_element()->num_node());
        Core::LinAlg::SerialDenseVector thermo_adjoint_test_master(
            mele.parent_element()->num_node());
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>
            deriv_thermo_adjoint_test_slave_d(
                sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(),
                -1, Core::LinAlg::SerialDenseVector(sele.parent_element()->num_node(), true));
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>
            deriv_thermo_adjoint_test_master_d(
                mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(),
                -1, Core::LinAlg::SerialDenseVector(mele.parent_element()->num_node(), true));
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>
            deriv_thermo_adjoint_test_slave_T(
                1, -1, Core::LinAlg::SerialDenseVector(sele.parent_element()->num_node(), true));
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>
            deriv_thermo_adjoint_test_master_T(
                1, -1, Core::LinAlg::SerialDenseVector(mele.parent_element()->num_node(), true));

        so_ele_cauchy_heatflux<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal,
            ws_thermo, qn_weighted_average, deriv_qn_weighted_average_d,
            deriv_qn_weighted_average_T, thermo_adjoint_test_slave,
            deriv_thermo_adjoint_test_slave_d, deriv_thermo_adjoint_test_slave_T);
        so_ele_cauchy_heatflux<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal,
            -wm_thermo, qn_weighted_average, deriv_qn_weighted_average_d,
            deriv_qn_weighted_average_T, thermo_adjoint_test_master,
            deriv_thermo_adjoint_test_master_d, deriv_thermo_adjoint_test_master_T);

        {
          double test_val = 0.;
          Core::Gen::Pairedvector<int, double> deriv_test_val_d(
              sele.parent_element()->num_node() * dim + mele.parent_element()->num_node() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          Core::Gen::Pairedvector<int, double> deriv_test_val_T(
              sele.parent_element()->num_node() + mele.parent_element()->num_node());

          const double beta_bar = beta * (-snn_av_pen_gap);
          test_val += +(beta_bar) / (pen_thermo + beta_bar) * qn_weighted_average;
          test_val += +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (s_gp_temp - m_gp_temp);
          test_val += -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * diss;

          for (const auto& p : deriv_qn_weighted_average_d)
            deriv_test_val_d[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : deriv_qn_weighted_average_T)
            deriv_test_val_T[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * p.second;

          for (const auto& p : d_s_gp_temp_dd)
            deriv_test_val_d[p.first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : d_m_gp_temp_dd)
            deriv_test_val_d[p.first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (-p.second);
          for (const auto& p : d_s_gp_temp_dT)
            deriv_test_val_T[p.first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : d_m_gp_temp_dT)
            deriv_test_val_T[p.first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (-p.second);

          for (const auto& p : d_diss_d)
            deriv_test_val_d[p.first] +=
                -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p.second;
          for (const auto& p : d_diss_T)
            deriv_test_val_T[p.first] +=
                -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p.second;

          for (const auto& p : d_snn_av_pen_gap)
          {
            deriv_test_val_d[p.first] += beta * (-p.second) /
                                         (std::pow(pen_thermo + beta_bar, 2.)) *
                                         (+pen_thermo * qn_weighted_average +
                                             pen_thermo * pen_thermo * (s_gp_temp - m_gp_temp) -
                                             pen_thermo * (1. - delta_c - ws_thermo) * diss);
          }
          for (const auto& p : cauchy_nn_weighted_average_deriv_T)
          {
            deriv_test_val_T[p.first] += beta * (-p.second) /
                                         (std::pow(pen_thermo + beta_bar, 2.)) *
                                         (+pen_thermo * qn_weighted_average +
                                             pen_thermo * pen_thermo * (s_gp_temp - m_gp_temp) -
                                             pen_thermo * (1. - delta_c - ws_thermo) * diss);
          }

          integrate_thermal_test<dim>(+1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              test_val, deriv_test_val_d, deriv_test_val_T);
          integrate_thermal_test<dim>(-1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
              test_val, deriv_test_val_d, deriv_test_val_T);

          integrate_thermal_test<dim>(-delta_c, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              diss, d_diss_d, d_diss_T);
          integrate_thermal_test<dim>(-(1 - delta_c), mele, mval, mderiv, dmxi, jac, jacintcellmap,
              wgt, diss, d_diss_d, d_diss_T);
        }

        if (abs(theta_thermo_) > 1.e-12)
        {
          double test_val = 0.;
          Core::Gen::Pairedvector<int, double> deriv_test_val_d(
              sele.parent_element()->num_node() * dim + mele.parent_element()->num_node() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          Core::Gen::Pairedvector<int, double> deriv_test_val_T(
              sele.parent_element()->num_node() + mele.parent_element()->num_node());

          const double beta_bar = beta * (-snn_av_pen_gap);
          test_val += -(1.) / (pen_thermo + beta_bar) * qn_weighted_average;
          test_val += +(beta_bar) / (pen_thermo + beta_bar) * (s_gp_temp - m_gp_temp);
          test_val += +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * diss;

          for (const auto& p : deriv_qn_weighted_average_d)
            deriv_test_val_d[p.first] += -(1.) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : deriv_qn_weighted_average_T)
            deriv_test_val_T[p.first] += -(1.) / (pen_thermo + beta_bar) * p.second;

          for (const auto& p : d_s_gp_temp_dd)
            deriv_test_val_d[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : d_s_gp_temp_dT)
            deriv_test_val_T[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * p.second;
          for (const auto& p : d_m_gp_temp_dd)
            deriv_test_val_d[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * (-p.second);
          for (const auto& p : d_m_gp_temp_dT)
            deriv_test_val_T[p.first] += +(beta_bar) / (pen_thermo + beta_bar) * (-p.second);

          for (const auto& p : d_diss_d)
            deriv_test_val_d[p.first] +=
                +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p.second;
          for (const auto& p : d_diss_T)
            deriv_test_val_T[p.first] +=
                +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p.second;

          for (const auto& p : d_snn_av_pen_gap)
          {
            deriv_test_val_d[p.first] +=
                (+1. * qn_weighted_average + pen_thermo * (s_gp_temp - m_gp_temp) -
                    1. * (1. - delta_c - ws_thermo) * diss) /
                std::pow(pen_thermo + beta_bar, 2) * beta * (-p.second);
          }

          for (const auto& p : cauchy_nn_weighted_average_deriv_T)
          {
            deriv_test_val_T[p.first] +=
                (+1. * qn_weighted_average + pen_thermo * (s_gp_temp - m_gp_temp) -
                    1. * (1. - delta_c - ws_thermo) * diss) /
                std::pow(pen_thermo + beta_bar, 2) * beta * (-p.second);
          }

          integrate_thermal_adjoint_test<dim>(theta_thermo_, jac, jacintcellmap, wgt, test_val,
              deriv_test_val_d, deriv_test_val_T, sele, thermo_adjoint_test_slave,
              deriv_thermo_adjoint_test_slave_d, deriv_thermo_adjoint_test_slave_T);
          integrate_thermal_adjoint_test<dim>(theta_thermo_, jac, jacintcellmap, wgt, test_val,
              deriv_test_val_d, deriv_test_val_T, mele, thermo_adjoint_test_master,
              deriv_thermo_adjoint_test_master_d, deriv_thermo_adjoint_test_master_T);
        }
        break;
      }
      default:
        FOUR_C_THROW(
            "unknown method for thermal constraint enforcement in Nitsche contact integrator");
        break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType parentdistype, int dim>
void inline CONTACT::IntegratorNitscheTsi::so_ele_gp(Mortar::Element& sele, const double wgt,
    const double* gpcoord, Core::LinAlg::Matrix<dim, 1>& pxsi,
    Core::LinAlg::Matrix<dim, dim>& derivtrafo)
{
  Core::FE::CollectedGaussPoints intpoints =
      Core::FE::CollectedGaussPoints(1);  // reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1], 0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  Core::LinAlg::SerialDenseMatrix pqxg(1, dim);
  derivtrafo.Clear();

  Core::FE::BoundaryGPToParentGP<dim>(pqxg, derivtrafo, intpoints, sele.parent_element()->Shape(),
      sele.Shape(), sele.FaceParentNumber());

  // coordinates of the current integration point in parent coordinate system
  for (int idim = 0; idim < dim; idim++) pxsi(idim) = pqxg(0, idim);
}


template <int dim>
void CONTACT::IntegratorNitscheTsi::so_ele_cauchy(Mortar::Element& moEle, double* boundary_gpcoord,
    std::vector<Core::Gen::Pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
    const Core::LinAlg::Matrix<dim, 1>& normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv,
    const Core::LinAlg::Matrix<dim, 1>& direction,
    std::vector<Core::Gen::Pairedvector<int, double>>& direction_deriv, const double w,
    double& cauchy_nt, Core::Gen::Pairedvector<int, double>& deriv_sigma_nt_d,
    Core::Gen::Pairedvector<int, double>& deriv_sigma_nt_T,
    Core::LinAlg::SerialDenseVector& adjoint_test,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  Core::LinAlg::Matrix<dim, 1> pxsi(true);
  Core::LinAlg::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(moEle, boundary_gpcoord, gp_wgt, pxsi, derivtravo_slave);

  double sigma_nt;
  Core::LinAlg::SerialDenseMatrix dsntdd, d2sntdd2, d2sntDdDn, d2sntDdDt, d2sntDdDpxi, d2sntDdDT,
      dsntdT;
  Core::LinAlg::Matrix<dim, 1> dsntdn, dsntdt, dsntdpxi;
  dynamic_cast<Discret::ELEMENTS::SoBase*>(moEle.parent_element())
      ->get_cauchy_n_dir_and_derivatives_at_xi(pxsi, moEle.MoData().ParentDisp(), normal, direction,
          sigma_nt, &dsntdd, &d2sntdd2, &d2sntDdDn, &d2sntDdDt, &d2sntDdDpxi, &dsntdn, &dsntdt,
          &dsntdpxi, &moEle.MoData().ParentTemp(), &dsntdT, &d2sntDdDT, nullptr, nullptr);

  cauchy_nt += w * sigma_nt;

  for (int i = 0; i < moEle.parent_element()->num_node() * dim; ++i)
    deriv_sigma_nt_d[moEle.MoData().ParentDof().at(i)] += w * dsntdd(i, 0);

  for (int i = 0; i < dim - 1; ++i)
  {
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = boundary_gpcoord_lin[i].begin();
         p != boundary_gpcoord_lin[i].end(); ++p)
    {
      double& ref = deriv_sigma_nt_d[p->first];
      for (int k = 0; k < dim; ++k) ref += dsntdpxi(k) * derivtravo_slave(k, i) * p->second * w;
    }
  }

  for (int d = 0; d < dim; ++d)
  {
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdn(d) * p->second * w;
  }

  for (int d = 0; d < dim; ++d)
  {
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = direction_deriv[d].begin();
         p != direction_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdt(d) * p->second * w;
  }

  if (moEle.MoData().ParentTempDof().size() != 0)
    for (int i = 0; i < moEle.parent_element()->num_node(); ++i)
      deriv_sigma_nt_T[moEle.MoData().ParentTempDof().at(i)] += dsntdT(i, 0) * w;

  if (abs(theta_) > 1.e-12)
  {
    build_adjoint_test<dim>(moEle, w, dsntdd, d2sntdd2, d2sntDdDn, d2sntDdDt, d2sntDdDpxi,
        boundary_gpcoord_lin, derivtravo_slave, normal_deriv, direction_deriv, adjoint_test,
        deriv_adjoint_test_d);
    build_adjoint_test_tsi<dim>(moEle, w, d2sntDdDT, deriv_adjoint_test_T);
  }
}
template <int dim>
void CONTACT::IntegratorNitscheTsi::integrate_test(const double fac, Mortar::Element& ele,
    const Core::LinAlg::SerialDenseVector& shape, const Core::LinAlg::SerialDenseMatrix& deriv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, const double jac,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test_val, const Core::Gen::Pairedvector<int, double>& test_deriv_d,
    const Core::Gen::Pairedvector<int, double>& test_deriv_T,
    const Core::LinAlg::Matrix<dim, 1>& test_dir,
    const std::vector<Core::Gen::Pairedvector<int, double>>& test_dir_deriv)
{
  CONTACT::IntegratorNitsche::integrate_test<dim>(fac, ele, shape, deriv, dxi, jac, jacintcellmap,
      wgt, test_val, test_deriv_d, test_dir, test_dir_deriv);

  for (const auto& p : test_deriv_T)
  {
    double* row = ele.GetNitscheContainer().Kdt(p.first);
    for (int s = 0; s < ele.num_node(); ++s)
    {
      for (int d = 0; d < Dim(); ++d)
      {
        row[Core::FE::getParentNodeNumberFromFaceNodeNumber(
                ele.parent_element()->Shape(), ele.FaceParentNumber(), s) *
                dim +
            d] += fac * jac * wgt * test_dir(d) * p.second * shape(s);
      }
    }
  }
}

template <int dim>
void CONTACT::IntegratorNitscheTsi::integrate_adjoint_test(const double fac, const double jac,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt, const double test,
    const Core::Gen::Pairedvector<int, double>& deriv_test_d,
    const Core::Gen::Pairedvector<int, double>& deriv_test_T, Mortar::Element& moEle,
    Core::LinAlg::SerialDenseVector& adjoint_test,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  CONTACT::IntegratorNitsche::integrate_adjoint_test<dim>(
      fac, jac, jacintcellmap, wgt, test, deriv_test_d, moEle, adjoint_test, deriv_adjoint_test_d);

  for (const auto& p : deriv_test_T)
  {
    Core::LinAlg::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Kdt(p.first), moEle.MoData().ParentDof().size());
    Core::LinAlg::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_adjoint_test_T)
  {
    Core::LinAlg::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Kdt(p.first), moEle.MoData().ParentDof().size());
    Core::LinAlg::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }
}


template <int dim>
void CONTACT::IntegratorNitscheTsi::integrate_thermal_test(const double fac, Mortar::Element& ele,
    const Core::LinAlg::SerialDenseVector& shape, const Core::LinAlg::SerialDenseMatrix& deriv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, const double jac,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test_val, const Core::Gen::Pairedvector<int, double>& test_deriv_d,
    const Core::Gen::Pairedvector<int, double>& test_deriv_T)
{
  double val = fac * jac * wgt * test_val;

  for (int s = 0; s < ele.num_node(); ++s)
    *(ele.GetNitscheContainer().RhsT(Core::FE::getParentNodeNumberFromFaceNodeNumber(
        ele.parent_element()->Shape(), ele.FaceParentNumber(), s))) += val * shape(s);

  Core::Gen::Pairedvector<int, double> val_deriv_d(jacintcellmap.size() + test_deriv_d.size());
  for (const auto& p : jacintcellmap) val_deriv_d[p.first] += fac * p.second * wgt * test_val;
  for (const auto& p : test_deriv_d) val_deriv_d[p.first] += fac * jac * wgt * p.second;

  for (const auto& p : val_deriv_d)
  {
    double* row = ele.GetNitscheContainer().Ktd(p.first);
    for (int s = 0; s < ele.num_node(); ++s)
      row[Core::FE::getParentNodeNumberFromFaceNodeNumber(
          ele.parent_element()->Shape(), ele.FaceParentNumber(), s)] += p.second * shape(s);
  }

  for (int e = 0; e < Dim() - 1; ++e)
  {
    for (auto p = dxi[e].begin(); p != dxi[e].end(); ++p)
    {
      double* row = ele.GetNitscheContainer().Ktd(p->first);
      for (int s = 0; s < ele.num_node(); ++s)
        row[Core::FE::getParentNodeNumberFromFaceNodeNumber(ele.parent_element()->Shape(),
            ele.FaceParentNumber(), s)] += val * deriv(s, e) * p->second;
    }
  }

  for (const auto& p : test_deriv_T)
  {
    double* row = ele.GetNitscheContainer().Ktt(p.first);
    for (int s = 0; s < ele.num_node(); ++s)
      row[Core::FE::getParentNodeNumberFromFaceNodeNumber(ele.parent_element()->Shape(),
          ele.FaceParentNumber(), s)] += fac * jac * wgt * p.second * shape(s);
  }
}

template <int dim>
void CONTACT::IntegratorNitscheTsi::integrate_thermal_adjoint_test(const double fac,
    const double jac, const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test, const Core::Gen::Pairedvector<int, double>& deriv_test_d,
    const Core::Gen::Pairedvector<int, double>& deriv_test_T, Mortar::Element& moEle,
    Core::LinAlg::SerialDenseVector& adjoint_test,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  Core::LinAlg::SerialDenseVector Tmp(
      Teuchos::View, moEle.GetNitscheContainer().RhsT(), moEle.parent_element()->num_node());
  Core::LinAlg::Update(fac * jac * wgt * test, adjoint_test, 1., Tmp);

  for (const auto& p : deriv_adjoint_test_d)
  {
    Core::LinAlg::SerialDenseVector Tmp(Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first),
        moEle.parent_element()->num_node());
    Core::LinAlg::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }

  for (const auto& p : jacintcellmap)
  {
    Core::LinAlg::SerialDenseVector Tmp(Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first),
        moEle.parent_element()->num_node());
    Core::LinAlg::Update(fac * p.second * wgt * test, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_test_d)
  {
    Core::LinAlg::SerialDenseVector Tmp(Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first),
        moEle.parent_element()->num_node());
    Core::LinAlg::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_adjoint_test_T)
  {
    Core::LinAlg::SerialDenseVector Tmp(Teuchos::View, moEle.GetNitscheContainer().Ktt(p.first),
        moEle.parent_element()->num_node());
    Core::LinAlg::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }

  for (const auto& p : deriv_test_T)
  {
    Core::LinAlg::SerialDenseVector Tmp(Teuchos::View, moEle.GetNitscheContainer().Ktt(p.first),
        moEle.parent_element()->num_node());
    Core::LinAlg::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }
}

template <int dim>
void CONTACT::IntegratorNitscheTsi::build_adjoint_test_tsi(Mortar::Element& moEle, const double fac,
    const Core::LinAlg::SerialDenseMatrix& d2sntDdDT,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (moEle.MoData().ParentTempDof().size())
  {
    for (int p = 0; p < moEle.parent_element()->num_node(); ++p)
    {
      Core::LinAlg::SerialDenseVector& at = deriv_adjoint_test_T[moEle.MoData().ParentTempDof()[p]];
      for (int i = 0; i < moEle.parent_element()->num_node() * dim; ++i)
        at(i) += fac * d2sntDdDT(i, p);
    }
  }
}

template <int dim>
void CONTACT::IntegratorNitscheTsi::setup_gp_temp(Mortar::Element& moEle,
    const Core::LinAlg::SerialDenseVector& val, const Core::LinAlg::SerialDenseMatrix& deriv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, double& temp,
    Core::Gen::Pairedvector<int, double>& d_temp_dT,
    Core::Gen::Pairedvector<int, double>& d_temp_dd)
{
  if (moEle.MoData().ParentTemp().size() == 0)
  {
    temp = 0.;
    return;
  }

  Core::LinAlg::SerialDenseVector moele_temp(val.length());
  for (int i = 0; i < moEle.num_node(); ++i)
  {
    moele_temp(i) = moEle.MoData().ParentTemp().at(Core::FE::getParentNodeNumberFromFaceNodeNumber(
        moEle.parent_element()->Shape(), moEle.FaceParentNumber(), i));
  }
  temp = val.dot(moele_temp);

  d_temp_dT.resize(val.length());
  d_temp_dT.clear();
  for (int i = 0; i < moEle.num_node(); ++i)
    d_temp_dT[moEle.MoData().ParentTempDof().at(Core::FE::getParentNodeNumberFromFaceNodeNumber(
        moEle.parent_element()->Shape(), moEle.FaceParentNumber(), i))] = val(i);

  int deriv_size = 0.;
  for (int i = 0; i < dim - 1; ++i) deriv_size += dxi.at(i).size();
  d_temp_dd.resize(deriv_size);
  d_temp_dd.clear();
  for (int i = 0; i < dim - 1; ++i)
  {
    for (auto p = dxi.at(i).begin(); p != dxi.at(i).end(); ++p)
    {
      double& a = d_temp_dd[p->first];
      for (int n = 0; n < moEle.num_node(); ++n) a += moele_temp(n) * deriv(n, i) * p->second;
    }
  }
}


template <int dim>
void CONTACT::IntegratorNitscheTsi::so_ele_cauchy_heatflux(Mortar::Element& moEle,
    double* boundary_gpcoord,
    const std::vector<Core::Gen::Pairedvector<int, double>>& boundary_gpcoord_lin,
    const double gp_wgt, const Core::LinAlg::Matrix<dim, 1>& normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv, const double w,
    double& heatflux, Core::Gen::Pairedvector<int, double>& dq_dd,
    Core::Gen::Pairedvector<int, double>& dq_dT, Core::LinAlg::SerialDenseVector& adjoint_test,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (moEle.MoData().ParentTemp().size() == 0)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    std::cout << "***************************************************************\n"
                 "WARNING: we are skipping the evaluation of the cauchy heat flux\n"
                 "         because parent temperatures are not set. This can\n"
                 "         happen in the constructor phase but may not happen\n"
                 "         during a running simulation\n"
                 "***************************************************************\n"
              << std::endl;
#endif
    return;
  }

  Core::LinAlg::Matrix<dim, 1> pxsi(true);
  Core::LinAlg::Matrix<dim, dim> derivtravo_slave;
  Core::FE::CellType distype = moEle.parent_element()->Shape();

  double q = 0;
  Core::LinAlg::SerialDenseMatrix dq_dT_ele, dq_dd_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi;
  Core::LinAlg::Matrix<dim, 1> dq_dn, dq_dpxi;

  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      auto* parent_ele = dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>*>(
          moEle.parent_element());
      if (!parent_ele) FOUR_C_THROW("thermo-mechanical Nitsche contact only for So3Plast for now.");

      so_ele_gp<Core::FE::CellType::hex8, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case Core::FE::CellType::hex27:
    {
      auto* parent_ele = dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>*>(
          moEle.parent_element());
      if (!parent_ele) FOUR_C_THROW("thermo-mechanical Nitsche contact only for So3Plast for now.");

      so_ele_gp<Core::FE::CellType::hex27, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case Core::FE::CellType::tet4:
    {
      auto* parent_ele = dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>*>(
          moEle.parent_element());
      if (!parent_ele) FOUR_C_THROW("thermo-mechanical Nitsche contact only for So3Plast for now.");

      so_ele_gp<Core::FE::CellType::tet4, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    default:
      FOUR_C_THROW("Nitsche contact not implemented for used (bulk) elements");
      break;
  }

  heatflux += w * q;

  for (int i = 0; i < moEle.parent_element()->num_node() * dim; ++i)
    dq_dd[moEle.MoData().ParentDof().at(i)] += w * dq_dd_ele(i, 0);

  for (int i = 0; i < dim - 1; ++i)
  {
    for (auto p = boundary_gpcoord_lin[i].begin(); p != boundary_gpcoord_lin[i].end(); ++p)
    {
      double& ref = dq_dd[p->first];
      for (int k = 0; k < dim; ++k) ref += dq_dpxi(k) * derivtravo_slave(k, i) * p->second * w;
    }
  }

  for (int d = 0; d < dim; ++d)
  {
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      dq_dd[p->first] += dq_dn(d) * p->second * w;
  }

  for (int i = 0; i < moEle.parent_element()->num_node(); ++i)
    dq_dT[moEle.MoData().ParentTempDof().at(i)] += w * dq_dT_ele(i, 0);

  if (abs(theta_thermo_) > 1.e-12)
  {
    build_adjoint_test_thermo<dim>(moEle, w, dq_dT_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi,
        normal_deriv, boundary_gpcoord_lin, derivtravo_slave, adjoint_test, deriv_adjoint_test_d,
        deriv_adjoint_test_T);
  }
}


template <int dim>
void CONTACT::IntegratorNitscheTsi::build_adjoint_test_thermo(Mortar::Element& moEle,
    const double fac, const Core::LinAlg::SerialDenseMatrix& dq_dT_ele,
    const Core::LinAlg::SerialDenseMatrix& d2q_dT_dd,
    const Core::LinAlg::SerialDenseMatrix& d2q_dT_dn,
    const Core::LinAlg::SerialDenseMatrix& d2q_dT_dpxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv,
    const std::vector<Core::Gen::Pairedvector<int, double>>& boundary_gpcoord_lin,
    Core::LinAlg::Matrix<dim, dim>& derivtravo_slave, Core::LinAlg::SerialDenseVector& adjoint_test,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T)
{
  for (int i = 0; i < moEle.parent_element()->num_node(); ++i)
    adjoint_test(i) = fac * dq_dT_ele(i, 0);

  for (int i = 0; i < moEle.parent_element()->num_node() * dim; ++i)
  {
    Core::LinAlg::SerialDenseVector& at = deriv_adjoint_test_d[moEle.MoData().ParentDof().at(i)];
    for (int j = 0; j < moEle.num_node(); ++j) at(j) += fac * d2q_dT_dd(j, i);
  }

  for (int d = 0; d < dim; ++d)
  {
    for (Core::Gen::Pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
    {
      Core::LinAlg::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.parent_element()->num_node(); ++i)
        at(i) += fac * d2q_dT_dn(i, d) * p->second;
    }
  }

  Core::LinAlg::SerialDenseMatrix tmp(moEle.parent_element()->num_node(), dim, false);
  Core::LinAlg::SerialDenseMatrix deriv_trafo(Teuchos::View, derivtravo_slave.A(),
      derivtravo_slave.numRows(), derivtravo_slave.numRows(), derivtravo_slave.numCols());
  if (Core::LinAlg::multiply(tmp, d2q_dT_dpxi, deriv_trafo)) FOUR_C_THROW("multiply failed");
  for (int d = 0; d < dim - 1; ++d)
  {
    for (auto p = boundary_gpcoord_lin[d].begin(); p != boundary_gpcoord_lin[d].end(); ++p)
    {
      Core::LinAlg::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.parent_element()->num_node(); ++i)
        at(i) += fac * tmp(i, d) * p->second;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
