/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#include "baci_contact_nitsche_integrator_tsi.H"

#include "baci_contact_element.H"
#include "baci_contact_nitsche_integrator.H"
#include "baci_contact_nitsche_utils.H"
#include "baci_contact_node.H"
#include "baci_contact_paramsinterface.H"
#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_linalg_utils_densematrix_multiply.H"
#include "baci_mat_elasthyper.H"
#include "baci_so3_base.H"
#include "baci_so3_plast_ssn.H"

#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheTsi::IntegrateGP_3D(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, CORE::LINALG::SerialDenseVector& sval,
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
  GPTSForces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt, gap,
      deriv_gap, normal, dnmap_unit, sxi, mxi);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheTsi::IntegrateGP_2D(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, CORE::LINALG::SerialDenseVector& sval,
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
  GPTSForces<2>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt, gap,
      deriv_gap, normal, dnmap_unit, sxi, mxi);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheTsi::GPTSForces(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const CORE::LINALG::SerialDenseVector& sval,
    const CORE::LINALG::SerialDenseMatrix& sderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
    const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
    std::vector<CORE::GEN::pairedvector<int, double>>& deriv_contact_normal, double* sxi,
    double* mxi)
{
  if (sele.Owner() != Comm_.MyPID()) return;

  if (dim != Dim()) dserror("dimension inconsistency");

  const CORE::GEN::pairedvector<int, double> empty(0);

  double s_gp_temp, m_gp_temp;
  CORE::GEN::pairedvector<int, double> d_s_gp_temp_dT(0), d_m_gp_temp_dT(0), d_s_gp_temp_dd(0),
      d_m_gp_temp_dd(0);
  SetupGpTemp<dim>(sele, sval, sderiv, dsxi, s_gp_temp, d_s_gp_temp_dT, d_s_gp_temp_dd);
  SetupGpTemp<dim>(mele, mval, mderiv, dmxi, m_gp_temp, d_m_gp_temp_dT, d_m_gp_temp_dd);

  CORE::LINALG::Matrix<dim, 1> xgp;
  for (int n = 0; n < sele.NumNode(); ++n)
    for (int d = 0; d < dim; ++d)
      xgp(d) += sval(n) * dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[n])->xspatial()[d];

  if (frtype_ != INPAR::CONTACT::friction_none && dim != 3) dserror("only 3D friction");
  if (frtype_ != INPAR::CONTACT::friction_none && frtype_ != INPAR::CONTACT::friction_coulomb &&
      frtype_ != INPAR::CONTACT::friction_tresca)
    dserror("only coulomb or tresca friction");
  if (frtype_ == INPAR::CONTACT::friction_coulomb && frcoeff_ < 0.)
    dserror("negative coulomb friction coefficient");
  if (frtype_ == INPAR::CONTACT::friction_tresca && frbound_ < 0.)
    dserror("negative tresca friction bound");

  CORE::LINALG::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<CORE::GEN::pairedvector<int, double>> deriv_slave_normal(0, 0);
  std::vector<CORE::GEN::pairedvector<int, double>> deriv_master_normal(0, 0);
  sele.ComputeUnitNormalAtXi(sxi, slave_normal.A());
  mele.ComputeUnitNormalAtXi(mxi, master_normal.A());
  sele.DerivUnitNormalAtXi(sxi, deriv_slave_normal);
  mele.DerivUnitNormalAtXi(mxi, deriv_master_normal);

  double pen = ppn_;
  double pet = ppt_;

  const CORE::LINALG::Matrix<dim, 1> contact_normal(gpn, true);
  double cauchy_nn_weighted_average = 0.;
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

  CORE::LINALG::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));

  CORE::LINALG::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);

  // variables for friction (declaration only)
  CORE::LINALG::Matrix<dim, 1> t1, t2;
  std::vector<CORE::GEN::pairedvector<int, double>> dt1, dt2;
  CORE::LINALG::Matrix<dim, 1> relVel;
  std::vector<CORE::GEN::pairedvector<int, double>> relVel_deriv(
      dim, sele.NumNode() * dim + mele.NumNode() * dim + dsxi[0].size() + dmxi[0].size());
  double vt1, vt2;
  CORE::GEN::pairedvector<int, double> dvt1(0);
  CORE::GEN::pairedvector<int, double> dvt2(0);
  double cauchy_nt1_weighted_average = 0.;
  CORE::GEN::pairedvector<int, double> cauchy_nt1_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  CORE::GEN::pairedvector<int, double> cauchy_nt1_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  CORE::LINALG::SerialDenseVector t1_adjoint_test_slave(sele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t1_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t1_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  CORE::LINALG::SerialDenseVector t1_adjoint_test_master(mele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t1_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t1_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double cauchy_nt2_weighted_average = 0.;
  CORE::GEN::pairedvector<int, double> cauchy_nt2_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  CORE::GEN::pairedvector<int, double> cauchy_nt2_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  CORE::LINALG::SerialDenseVector t2_adjoint_test_slave(sele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t2_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t2_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  CORE::LINALG::SerialDenseVector t2_adjoint_test_master(mele.MoData().ParentDof().size());
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t2_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_t2_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  double sigma_nt1_pen_vt1 = 0;
  double sigma_nt2_pen_vt2 = 0;
  CORE::GEN::pairedvector<int, double> d_sigma_nt1_pen_vt1(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt1_weighted_average_deriv.capacity() + dvt1.capacity(),
      0, 0);
  CORE::GEN::pairedvector<int, double> d_sigma_nt2_pen_vt2(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt2_weighted_average_deriv.capacity() + dvt2.capacity(),
      0, 0);
  CORE::GEN::pairedvector<int, double> d_sigma_nt1_pen_vt1_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  CORE::GEN::pairedvector<int, double> d_sigma_nt2_pen_vt2_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  // variables for friction (end)

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, contact_normal,
      deriv_contact_normal, ws, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
      cauchy_nn_weighted_average_deriv_T, normal_adjoint_test_slave,
      deriv_normal_adjoint_test_slave, deriv_normal_adjoint_test_slave_T);

  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, contact_normal,
      deriv_contact_normal, -wm, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
      cauchy_nn_weighted_average_deriv_T, normal_adjoint_test_master,
      deriv_normal_adjoint_test_master, deriv_normal_adjoint_test_master_T);

  const double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  CORE::GEN::pairedvector<int, double> d_snn_av_pen_gap(
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

    SoEleCauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, t1, dt1, ws,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1_adjoint_test_slave, deriv_t1_adjoint_test_slave,
        deriv_t1_adjoint_test_slave_T);
    SoEleCauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, t1, dt1, -wm,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1_adjoint_test_master, deriv_t1_adjoint_test_master,
        deriv_t1_adjoint_test_master_T);

    SoEleCauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, t2, dt2, ws,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2_adjoint_test_slave, deriv_t2_adjoint_test_slave,
        deriv_t2_adjoint_test_slave_T);
    SoEleCauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, t2, dt2, -wm,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2_adjoint_test_master, deriv_t2_adjoint_test_master,
        deriv_t2_adjoint_test_master_T);
  }  // evaluation of tangential stuff


  if (frtype_)
  {
    IntegrateTest<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1, dt1);
    IntegrateTest<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2, dt2);
    IntegrateTest<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nt1_weighted_average, cauchy_nt1_weighted_average_deriv,
        cauchy_nt1_weighted_average_deriv_T, t1, dt1);
    IntegrateTest<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nt2_weighted_average, cauchy_nt2_weighted_average_deriv,
        cauchy_nt2_weighted_average_deriv_T, t2, dt2);

    IntegrateAdjointTest<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt1_weighted_average,
        cauchy_nt1_weighted_average_deriv, cauchy_nt1_weighted_average_deriv_T, sele,
        t1_adjoint_test_slave, deriv_t1_adjoint_test_slave, deriv_t1_adjoint_test_slave_T);
    IntegrateAdjointTest<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt2_weighted_average,
        cauchy_nt2_weighted_average_deriv, cauchy_nt2_weighted_average_deriv_T, sele,
        t2_adjoint_test_slave, deriv_t2_adjoint_test_slave, deriv_t2_adjoint_test_slave_T);

    IntegrateAdjointTest<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt1_weighted_average,
        cauchy_nt1_weighted_average_deriv, cauchy_nt1_weighted_average_deriv_T, mele,
        t1_adjoint_test_master, deriv_t1_adjoint_test_master, deriv_t1_adjoint_test_master_T);
    IntegrateAdjointTest<dim>(-theta_ / pet, jac, jacintcellmap, wgt, cauchy_nt2_weighted_average,
        cauchy_nt2_weighted_average_deriv, cauchy_nt2_weighted_average_deriv_T, mele,
        t2_adjoint_test_master, deriv_t2_adjoint_test_master, deriv_t2_adjoint_test_master_T);
  }

  if (snn_av_pen_gap >= 0.)
  {
    IntegrateTest<dim>(-1. + theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);
    IntegrateTest<dim>(+1. - theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);

    IntegrateAdjointTest<dim>(-theta_ / pen, jac, jacintcellmap, wgt, cauchy_nn_weighted_average,
        cauchy_nn_weighted_average_deriv, cauchy_nn_weighted_average_deriv_T, sele,
        normal_adjoint_test_slave, deriv_normal_adjoint_test_slave,
        deriv_normal_adjoint_test_slave_T);
    IntegrateAdjointTest<dim>(-theta_ / pen, jac, jacintcellmap, wgt, cauchy_nn_weighted_average,
        cauchy_nn_weighted_average_deriv, cauchy_nn_weighted_average_deriv_T, mele,
        normal_adjoint_test_master, deriv_normal_adjoint_test_master,
        deriv_normal_adjoint_test_master_T);
  }
  else
  {
    // test in normal contact direction
    IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);
    IntegrateTest<dim>(+1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
        cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv,
        cauchy_nn_weighted_average_deriv_T, contact_normal, deriv_contact_normal);

    IntegrateTest<dim>(-theta_2_ * pen, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, gap,
        dgapgp, empty, contact_normal, deriv_contact_normal);
    IntegrateTest<dim>(+theta_2_ * pen, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt, gap,
        dgapgp, empty, contact_normal, deriv_contact_normal);

    IntegrateAdjointTest<dim>(theta_, jac, jacintcellmap, wgt, gap, dgapgp, empty, sele,
        normal_adjoint_test_slave, deriv_normal_adjoint_test_slave,
        deriv_normal_adjoint_test_slave_T);
    IntegrateAdjointTest<dim>(theta_, jac, jacintcellmap, wgt, gap, dgapgp, empty, mele,
        normal_adjoint_test_master, deriv_normal_adjoint_test_master,
        deriv_normal_adjoint_test_master_T);

    if (frtype_)
    {
      double fr = 0.;
      CORE::GEN::pairedvector<int, double> d_fr_d(
          d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
      CORE::GEN::pairedvector<int, double> d_fr_T(cauchy_nn_weighted_average_deriv_T.size() +
                                                  d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
      switch (frtype_)
      {
        case INPAR::CONTACT::friction_coulomb:
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
        case INPAR::CONTACT::friction_tresca:
          fr = frbound_;
          break;
        default:
          fr = 0.;
          dserror("why are you here???");
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
        CORE::GEN::pairedvector<int, double> tmp_d(
            dgapgp.size() + cauchy_nn_weighted_average_deriv.size() +
                cauchy_nt1_weighted_average_deriv.size() + dvt1.size(),
            0, 0);
        CORE::GEN::pairedvector<int, double> tmp_T(
            sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
        if (frtype_ == INPAR::CONTACT::friction_coulomb)
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
      IntegrateTest<dim>(-theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
          sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, t1, dt1);
      IntegrateTest<dim>(+theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
          sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, t1, dt1);
      IntegrateTest<dim>(-theta_2_, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
          sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, t2, dt2);
      IntegrateTest<dim>(+theta_2_, mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
          sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, t2, dt2);


      IntegrateAdjointTest<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt1_pen_vt1,
          d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, sele, t1_adjoint_test_slave,
          deriv_t1_adjoint_test_slave, deriv_t1_adjoint_test_slave_T);
      IntegrateAdjointTest<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt1_pen_vt1,
          d_sigma_nt1_pen_vt1, d_sigma_nt1_pen_vt1_T, mele, t1_adjoint_test_master,
          deriv_t1_adjoint_test_master, deriv_t1_adjoint_test_master_T);
      IntegrateAdjointTest<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt2_pen_vt2,
          d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, sele, t2_adjoint_test_slave,
          deriv_t2_adjoint_test_slave, deriv_t2_adjoint_test_slave_T);
      IntegrateAdjointTest<dim>(theta_ / pet, jac, jacintcellmap, wgt, sigma_nt2_pen_vt2,
          d_sigma_nt2_pen_vt2, d_sigma_nt2_pen_vt2_T, mele, t2_adjoint_test_master,
          deriv_t2_adjoint_test_master, deriv_t2_adjoint_test_master_T);
    }

    // ----------------------------------------------
    // thermo-stuff
    // ----------------------------------------------
    const double beta = gamma_slave_ * gamma_master_ / (gamma_slave_ + gamma_master_);
    const double delta_c = gamma_slave_ / (gamma_slave_ + gamma_master_);
    double diss = 0.;
    CORE::GEN::pairedvector<int, double> d_diss_d(
        d_sigma_nt1_pen_vt1.size() + d_sigma_nt2_pen_vt2.size() + dvt1.size() + dvt2.size());
    CORE::GEN::pairedvector<int, double> d_diss_T(
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
      case INPAR::CONTACT::NitThr_substitution:
      {
        const double beta_bar = beta * (-snn_av_pen_gap);
        const double q1 = beta_bar * (s_gp_temp - m_gp_temp);

        CORE::GEN::pairedvector<int, double> d_q1_d(
            d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
        for (const auto& p : d_snn_av_pen_gap)
          d_q1_d[p.first] += beta * (-p.second) * (s_gp_temp - m_gp_temp);
        for (const auto& p : d_s_gp_temp_dd) d_q1_d[p.first] += beta_bar * p.second;
        for (const auto& p : d_m_gp_temp_dd) d_q1_d[p.first] += beta_bar * (-p.second);

        CORE::GEN::pairedvector<int, double> d_q1_T(cauchy_nn_weighted_average_deriv_T.size() +
                                                    d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
        for (const auto& p : cauchy_nn_weighted_average_deriv_T)
          d_q1_T[p.first] += beta * (-p.second) * (s_gp_temp - m_gp_temp);
        for (const auto& p : d_s_gp_temp_dT) d_q1_T[p.first] += beta_bar * p.second;
        for (const auto& p : d_m_gp_temp_dT) d_q1_T[p.first] += beta_bar * (-p.second);

        IntegrateThermalTest<dim>(
            +1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, q1, d_q1_d, d_q1_T);
        IntegrateThermalTest<dim>(
            -1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt, q1, d_q1_d, d_q1_T);

        if (frtype_)
        {
          IntegrateThermalTest<dim>(-delta_c, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              diss, d_diss_d, d_diss_T);
          IntegrateThermalTest<dim>(-(1 - delta_c), mele, mval, mderiv, dmxi, jac, jacintcellmap,
              wgt, diss, d_diss_d, d_diss_T);
        }
        break;
      }
      case INPAR::CONTACT::NitThr_nitsche:
      {
        double pen_thermo = pp_thermo_;
        double ws_thermo = 0.;
        double wm_thermo = 0.;

        switch (nit_wgt_)
        {
          case INPAR::CONTACT::NitWgt_slave:
            ws_thermo = 1.;
            wm_thermo = 0.;
            pen_thermo /= dynamic_cast<CONTACT::CoElement&>(sele).TraceHCond();
            break;
          case INPAR::CONTACT::NitWgt_master:
            ws_thermo = 0.;
            wm_thermo = 1.;
            pen_thermo /= dynamic_cast<CONTACT::CoElement&>(mele).TraceHCond();
            break;
          case INPAR::CONTACT::NitWgt_harmonic:
            ws_thermo = 1. / dynamic_cast<CONTACT::CoElement&>(mele).TraceHCond();
            ws_thermo /= (ws_thermo + wm_thermo);
            wm_thermo = 1. - ws_thermo;
            pen_thermo =
                ws_thermo * pen_thermo / dynamic_cast<CONTACT::CoElement&>(sele).TraceHCond() +
                wm_thermo * pen_thermo / dynamic_cast<CONTACT::CoElement&>(mele).TraceHCond();
            break;
          case INPAR::CONTACT::NitWgt_phyiscal:
            ws_thermo = 1. - delta_c;
            wm_thermo = delta_c;
            pen_thermo =
                ws_thermo * pen_thermo / dynamic_cast<CONTACT::CoElement&>(sele).TraceHCond() +
                wm_thermo * pen_thermo / dynamic_cast<CONTACT::CoElement&>(mele).TraceHCond();
            break;
          default:
            dserror("unknown Nitsche weighting");
            break;
        }

        double qn_weighted_average = 0.;
        CORE::GEN::pairedvector<int, double> deriv_qn_weighted_average_d(
            sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
            deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
        CORE::GEN::pairedvector<int, double> deriv_qn_weighted_average_T(
            sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
        CORE::LINALG::SerialDenseVector thermo_adjoint_test_slave(sele.ParentElement()->NumNode());
        CORE::LINALG::SerialDenseVector thermo_adjoint_test_master(mele.ParentElement()->NumNode());
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>
            deriv_thermo_adjoint_test_slave_d(
                sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(),
                -1, CORE::LINALG::SerialDenseVector(sele.ParentElement()->NumNode(), true));
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>
            deriv_thermo_adjoint_test_master_d(
                mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(),
                -1, CORE::LINALG::SerialDenseVector(mele.ParentElement()->NumNode(), true));
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>
            deriv_thermo_adjoint_test_slave_T(
                1, -1, CORE::LINALG::SerialDenseVector(sele.ParentElement()->NumNode(), true));
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>
            deriv_thermo_adjoint_test_master_T(
                1, -1, CORE::LINALG::SerialDenseVector(mele.ParentElement()->NumNode(), true));

        SoEleCauchyHeatflux<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, ws_thermo,
            qn_weighted_average, deriv_qn_weighted_average_d, deriv_qn_weighted_average_T,
            thermo_adjoint_test_slave, deriv_thermo_adjoint_test_slave_d,
            deriv_thermo_adjoint_test_slave_T);
        SoEleCauchyHeatflux<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal,
            -wm_thermo, qn_weighted_average, deriv_qn_weighted_average_d,
            deriv_qn_weighted_average_T, thermo_adjoint_test_master,
            deriv_thermo_adjoint_test_master_d, deriv_thermo_adjoint_test_master_T);

        {
          double test_val = 0.;
          CORE::GEN::pairedvector<int, double> deriv_test_val_d(
              sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          CORE::GEN::pairedvector<int, double> deriv_test_val_T(
              sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

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

          IntegrateThermalTest<dim>(+1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              test_val, deriv_test_val_d, deriv_test_val_T);
          IntegrateThermalTest<dim>(-1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt,
              test_val, deriv_test_val_d, deriv_test_val_T);

          IntegrateThermalTest<dim>(-delta_c, sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
              diss, d_diss_d, d_diss_T);
          IntegrateThermalTest<dim>(-(1 - delta_c), mele, mval, mderiv, dmxi, jac, jacintcellmap,
              wgt, diss, d_diss_d, d_diss_T);
        }

        if (abs(theta_thermo_) > 1.e-12)
        {
          double test_val = 0.;
          CORE::GEN::pairedvector<int, double> deriv_test_val_d(
              sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          CORE::GEN::pairedvector<int, double> deriv_test_val_T(
              sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

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

          IntegrateThermalAdjointTest<dim>(theta_thermo_, jac, jacintcellmap, wgt, test_val,
              deriv_test_val_d, deriv_test_val_T, sele, thermo_adjoint_test_slave,
              deriv_thermo_adjoint_test_slave_d, deriv_thermo_adjoint_test_slave_T);
          IntegrateThermalAdjointTest<dim>(theta_thermo_, jac, jacintcellmap, wgt, test_val,
              deriv_test_val_d, deriv_test_val_T, mele, thermo_adjoint_test_master,
              deriv_thermo_adjoint_test_master_d, deriv_thermo_adjoint_test_master_T);
        }
        break;
      }
      default:
        dserror("unknown method for thermal constraint enforcement in Nitsche contact integrator");
        break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType parentdistype, int dim>
void inline CONTACT::CoIntegratorNitscheTsi::SoEleGP(MORTAR::MortarElement& sele, const double wgt,
    const double* gpcoord, CORE::LINALG::Matrix<dim, 1>& pxsi,
    CORE::LINALG::Matrix<dim, dim>& derivtrafo)
{
  CORE::DRT::UTILS::CollectedGaussPoints intpoints =
      CORE::DRT::UTILS::CollectedGaussPoints(1);  // reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1], 0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  CORE::LINALG::SerialDenseMatrix pqxg(1, dim);
  derivtrafo.Clear();

  CORE::DRT::UTILS::BoundaryGPToParentGP<dim>(pqxg, derivtrafo, intpoints,
      sele.ParentElement()->Shape(), sele.Shape(), sele.FaceParentNumber());

  // coordinates of the current integration point in parent coordinate system
  for (int idim = 0; idim < dim; idim++) pxsi(idim) = pqxg(0, idim);
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SoEleCauchy(MORTAR::MortarElement& moEle,
    double* boundary_gpcoord,
    std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
    const CORE::LINALG::Matrix<dim, 1>& normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
    const CORE::LINALG::Matrix<dim, 1>& direction,
    std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, const double w,
    double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_d,
    CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_T,
    CORE::LINALG::SerialDenseVector& adjoint_test,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  CORE::LINALG::Matrix<dim, 1> pxsi(true);
  CORE::LINALG::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(moEle, boundary_gpcoord, gp_wgt, pxsi, derivtravo_slave);

  double sigma_nt;
  CORE::LINALG::SerialDenseMatrix dsntdd, d2sntdd2, d2sntDdDn, d2sntDdDt, d2sntDdDpxi, d2sntDdDT,
      dsntdT;
  CORE::LINALG::Matrix<dim, 1> dsntdn, dsntdt, dsntdpxi;
  dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())
      ->GetCauchyNDirAndDerivativesAtXi(pxsi, moEle.MoData().ParentDisp(), normal, direction,
          sigma_nt, &dsntdd, &d2sntdd2, &d2sntDdDn, &d2sntDdDt, &d2sntDdDpxi, &dsntdn, &dsntdt,
          &dsntdpxi, &moEle.MoData().ParentTemp(), &dsntdT, &d2sntDdDT, nullptr, nullptr);

  cauchy_nt += w * sigma_nt;

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
    deriv_sigma_nt_d[moEle.MoData().ParentDof().at(i)] += w * dsntdd(i, 0);

  for (int i = 0; i < dim - 1; ++i)
  {
    for (CORE::GEN::pairedvector<int, double>::const_iterator p = boundary_gpcoord_lin[i].begin();
         p != boundary_gpcoord_lin[i].end(); ++p)
    {
      double& ref = deriv_sigma_nt_d[p->first];
      for (int k = 0; k < dim; ++k) ref += dsntdpxi(k) * derivtravo_slave(k, i) * p->second * w;
    }
  }

  for (int d = 0; d < dim; ++d)
  {
    for (CORE::GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdn(d) * p->second * w;
  }

  for (int d = 0; d < dim; ++d)
  {
    for (CORE::GEN::pairedvector<int, double>::const_iterator p = direction_deriv[d].begin();
         p != direction_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdt(d) * p->second * w;
  }

  if (moEle.MoData().ParentTempDof().size() != 0)
    for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
      deriv_sigma_nt_T[moEle.MoData().ParentTempDof().at(i)] += dsntdT(i, 0) * w;

  if (abs(theta_) > 1.e-12)
  {
    BuildAdjointTest<dim>(moEle, w, dsntdd, d2sntdd2, d2sntDdDn, d2sntDdDt, d2sntDdDpxi,
        boundary_gpcoord_lin, derivtravo_slave, normal_deriv, direction_deriv, adjoint_test,
        deriv_adjoint_test_d);
    BuildAdjointTestTsi<dim>(moEle, w, d2sntDdDT, deriv_adjoint_test_T);
  }
}
template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateTest(const double fac, MORTAR::MortarElement& ele,
    const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test_val, const CORE::GEN::pairedvector<int, double>& test_deriv_d,
    const CORE::GEN::pairedvector<int, double>& test_deriv_T,
    const CORE::LINALG::Matrix<dim, 1>& test_dir,
    const std::vector<CORE::GEN::pairedvector<int, double>>& test_dir_deriv)
{
  CONTACT::CoIntegratorNitsche::IntegrateTest<dim>(fac, ele, shape, deriv, dxi, jac, jacintcellmap,
      wgt, test_val, test_deriv_d, test_dir, test_dir_deriv);

  for (const auto& p : test_deriv_T)
  {
    double* row = ele.GetNitscheContainer().Kdt(p.first);
    for (int s = 0; s < ele.NumNode(); ++s)
    {
      for (int d = 0; d < Dim(); ++d)
      {
        row[CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
                ele.ParentElement()->Shape(), ele.FaceParentNumber(), s) *
                dim +
            d] += fac * jac * wgt * test_dir(d) * p.second * shape(s);
      }
    }
  }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateAdjointTest(const double fac, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test,
    const CORE::GEN::pairedvector<int, double>& deriv_test_d,
    const CORE::GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
    CORE::LINALG::SerialDenseVector& adjoint_test,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  CONTACT::CoIntegratorNitsche::IntegrateAdjointTest<dim>(
      fac, jac, jacintcellmap, wgt, test, deriv_test_d, moEle, adjoint_test, deriv_adjoint_test_d);

  for (const auto& p : deriv_test_T)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Kdt(p.first), moEle.MoData().ParentDof().size());
    CORE::LINALG::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_adjoint_test_T)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Kdt(p.first), moEle.MoData().ParentDof().size());
    CORE::LINALG::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateThermalTest(const double fac,
    MORTAR::MortarElement& ele, const CORE::LINALG::SerialDenseVector& shape,
    const CORE::LINALG::SerialDenseMatrix& deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test_val, const CORE::GEN::pairedvector<int, double>& test_deriv_d,
    const CORE::GEN::pairedvector<int, double>& test_deriv_T)
{
  double val = fac * jac * wgt * test_val;

  for (int s = 0; s < ele.NumNode(); ++s)
    *(ele.GetNitscheContainer().RhsT(CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), s))) += val * shape(s);

  CORE::GEN::pairedvector<int, double> val_deriv_d(jacintcellmap.size() + test_deriv_d.size());
  for (const auto& p : jacintcellmap) val_deriv_d[p.first] += fac * p.second * wgt * test_val;
  for (const auto& p : test_deriv_d) val_deriv_d[p.first] += fac * jac * wgt * p.second;

  for (const auto& p : val_deriv_d)
  {
    double* row = ele.GetNitscheContainer().Ktd(p.first);
    for (int s = 0; s < ele.NumNode(); ++s)
      row[CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), s)] += p.second * shape(s);
  }

  for (int e = 0; e < Dim() - 1; ++e)
  {
    for (auto p = dxi[e].begin(); p != dxi[e].end(); ++p)
    {
      double* row = ele.GetNitscheContainer().Ktd(p->first);
      for (int s = 0; s < ele.NumNode(); ++s)
        row[CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
            ele.FaceParentNumber(), s)] += val * deriv(s, e) * p->second;
    }
  }

  for (const auto& p : test_deriv_T)
  {
    double* row = ele.GetNitscheContainer().Ktt(p.first);
    for (int s = 0; s < ele.NumNode(); ++s)
      row[CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
          ele.FaceParentNumber(), s)] += fac * jac * wgt * p.second * shape(s);
  }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateThermalAdjointTest(const double fac,
    const double jac, const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test, const CORE::GEN::pairedvector<int, double>& deriv_test_d,
    const CORE::GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
    CORE::LINALG::SerialDenseVector& adjoint_test,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  CORE::LINALG::SerialDenseVector Tmp(
      Teuchos::View, moEle.GetNitscheContainer().RhsT(), moEle.ParentElement()->NumNode());
  CORE::LINALG::Update(fac * jac * wgt * test, adjoint_test, 1., Tmp);

  for (const auto& p : deriv_adjoint_test_d)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first), moEle.ParentElement()->NumNode());
    CORE::LINALG::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }

  for (const auto& p : jacintcellmap)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first), moEle.ParentElement()->NumNode());
    CORE::LINALG::Update(fac * p.second * wgt * test, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_test_d)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Ktd(p.first), moEle.ParentElement()->NumNode());
    CORE::LINALG::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }

  for (const auto& p : deriv_adjoint_test_T)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Ktt(p.first), moEle.ParentElement()->NumNode());
    CORE::LINALG::Update(fac * jac * wgt * test, p.second, 1., Tmp);
  }

  for (const auto& p : deriv_test_T)
  {
    CORE::LINALG::SerialDenseVector Tmp(
        Teuchos::View, moEle.GetNitscheContainer().Ktt(p.first), moEle.ParentElement()->NumNode());
    CORE::LINALG::Update(fac * jac * wgt * p.second, adjoint_test, 1., Tmp);
  }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::BuildAdjointTestTsi(MORTAR::MortarElement& moEle,
    const double fac, const CORE::LINALG::SerialDenseMatrix& d2sntDdDT,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (moEle.MoData().ParentTempDof().size())
  {
    for (int p = 0; p < moEle.ParentElement()->NumNode(); ++p)
    {
      CORE::LINALG::SerialDenseVector& at = deriv_adjoint_test_T[moEle.MoData().ParentTempDof()[p]];
      for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
        at(i) += fac * d2sntDdDT(i, p);
    }
  }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SetupGpTemp(MORTAR::MortarElement& moEle,
    const CORE::LINALG::SerialDenseVector& val, const CORE::LINALG::SerialDenseMatrix& deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, double& temp,
    CORE::GEN::pairedvector<int, double>& d_temp_dT,
    CORE::GEN::pairedvector<int, double>& d_temp_dd)
{
  if (moEle.MoData().ParentTemp().size() == 0)
  {
    temp = 0.;
    return;
  }

  CORE::LINALG::SerialDenseVector moele_temp(val.length());
  for (int i = 0; i < moEle.NumNode(); ++i)
  {
    moele_temp(i) =
        moEle.MoData().ParentTemp().at(CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
            moEle.ParentElement()->Shape(), moEle.FaceParentNumber(), i));
  }
  temp = val.dot(moele_temp);

  d_temp_dT.resize(val.length());
  d_temp_dT.clear();
  for (int i = 0; i < moEle.NumNode(); ++i)
    d_temp_dT[moEle.MoData().ParentTempDof().at(
        CORE::DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
            moEle.ParentElement()->Shape(), moEle.FaceParentNumber(), i))] = val(i);

  int deriv_size = 0.;
  for (int i = 0; i < dim - 1; ++i) deriv_size += dxi.at(i).size();
  d_temp_dd.resize(deriv_size);
  d_temp_dd.clear();
  for (int i = 0; i < dim - 1; ++i)
  {
    for (auto p = dxi.at(i).begin(); p != dxi.at(i).end(); ++p)
    {
      double& a = d_temp_dd[p->first];
      for (int n = 0; n < moEle.NumNode(); ++n) a += moele_temp(n) * deriv(n, i) * p->second;
    }
  }
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SoEleCauchyHeatflux(MORTAR::MortarElement& moEle,
    double* boundary_gpcoord,
    const std::vector<CORE::GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
    const double gp_wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv, const double w,
    double& heatflux, CORE::GEN::pairedvector<int, double>& dq_dd,
    CORE::GEN::pairedvector<int, double>& dq_dT, CORE::LINALG::SerialDenseVector& adjoint_test,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (moEle.MoData().ParentTemp().size() == 0)
  {
#ifdef DEBUG
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

  CORE::LINALG::Matrix<dim, 1> pxsi(true);
  CORE::LINALG::Matrix<dim, dim> derivtravo_slave;
  CORE::FE::CellType distype = moEle.ParentElement()->Shape();

  double q = 0;
  CORE::LINALG::SerialDenseMatrix dq_dT_ele, dq_dd_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi;
  CORE::LINALG::Matrix<dim, 1> dq_dn, dq_dpxi;

  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      auto* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<CORE::FE::CellType::hex8, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      auto* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<CORE::FE::CellType::hex27, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      auto* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<CORE::FE::CellType::tet4, dim>(
          moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    default:
      dserror("Nitsche contact not implemented for used (bulk) elements");
      break;
  }

  heatflux += w * q;

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
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
    for (CORE::GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      dq_dd[p->first] += dq_dn(d) * p->second * w;
  }

  for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
    dq_dT[moEle.MoData().ParentTempDof().at(i)] += w * dq_dT_ele(i, 0);

  if (abs(theta_thermo_) > 1.e-12)
  {
    BuildAdjointTestThermo<dim>(moEle, w, dq_dT_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi,
        normal_deriv, boundary_gpcoord_lin, derivtravo_slave, adjoint_test, deriv_adjoint_test_d,
        deriv_adjoint_test_T);
  }
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::BuildAdjointTestThermo(MORTAR::MortarElement& moEle,
    const double fac, const CORE::LINALG::SerialDenseMatrix& dq_dT_ele,
    const CORE::LINALG::SerialDenseMatrix& d2q_dT_dd,
    const CORE::LINALG::SerialDenseMatrix& d2q_dT_dn,
    const CORE::LINALG::SerialDenseMatrix& d2q_dT_dpxi,
    std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
    CORE::LINALG::Matrix<dim, dim>& derivtravo_slave, CORE::LINALG::SerialDenseVector& adjoint_test,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
    adjoint_test(i) = fac * dq_dT_ele(i, 0);

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
  {
    CORE::LINALG::SerialDenseVector& at = deriv_adjoint_test_d[moEle.MoData().ParentDof().at(i)];
    for (int j = 0; j < moEle.NumNode(); ++j) at(j) += fac * d2q_dT_dd(j, i);
  }

  for (int d = 0; d < dim; ++d)
  {
    for (CORE::GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
    {
      CORE::LINALG::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
        at(i) += fac * d2q_dT_dn(i, d) * p->second;
    }
  }

  CORE::LINALG::SerialDenseMatrix tmp(moEle.ParentElement()->NumNode(), dim, false);
  CORE::LINALG::SerialDenseMatrix deriv_trafo(Teuchos::View, derivtravo_slave.A(),
      derivtravo_slave.numRows(), derivtravo_slave.numRows(), derivtravo_slave.numCols());
  if (CORE::LINALG::multiply(tmp, d2q_dT_dpxi, deriv_trafo)) dserror("multiply failed");
  for (int d = 0; d < dim - 1; ++d)
  {
    for (auto p = boundary_gpcoord_lin[d].begin(); p != boundary_gpcoord_lin[d].end(); ++p)
    {
      CORE::LINALG::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
        at(i) += fac * tmp(i, d) * p->second;
    }
  }
}

BACI_NAMESPACE_CLOSE
