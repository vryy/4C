/*---------------------------------------------------------------------*/
/*!
\brief A class to perform integrations of nitsche related terms

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#include "contact_nitsche_integrator_tsi.H"
#include "contact_nitsche_integrator.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_paramsinterface.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_so3/so_base.H"
#include "../drt_so3/so3_plast/so3_ssn_plast.H"

#include "../drt_mat/elasthyper.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include "../linalg/linalg_utils.H"
#include "contact_nitsche_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheTsi::IntegrateGP_3D(MORTAR::MortarElement& sele,
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
  GPTS_forces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
      gap, deriv_gap, normal, dnmap_unit, sxi, mxi);

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheTsi::IntegrateGP_2D(MORTAR::MortarElement& sele,
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
  GPTS_forces<2>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
      gap, deriv_gap, normal, dnmap_unit, sxi, mxi);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheTsi::GPTS_forces(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const LINALG::SerialDenseVector& sval,
    const LINALG::SerialDenseMatrix& sderiv,
    const std::vector<GEN::pairedvector<int, double>>& dsxi, const LINALG::SerialDenseVector& mval,
    const LINALG::SerialDenseMatrix& mderiv,
    const std::vector<GEN::pairedvector<int, double>>& dmxi, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const GEN::pairedvector<int, double>& dgapgp, double* gpn,
    std::vector<GEN::pairedvector<int, double>>& deriv_contact_normal, double* sxi, double* mxi)
{
  if (sele.Owner() != Comm_.MyPID()) return;

  if (dim != Dim()) dserror("dimension inconsistency");

  const GEN::pairedvector<int, double> empty(0);

  double s_gp_temp, m_gp_temp;
  GEN::pairedvector<int, double> d_s_gp_temp_dT(0), d_m_gp_temp_dT(0), d_s_gp_temp_dd(0),
      d_m_gp_temp_dd(0);
  SetupGpTemp<dim>(sele, sval, sderiv, dsxi, s_gp_temp, d_s_gp_temp_dT, d_s_gp_temp_dd);
  SetupGpTemp<dim>(mele, mval, mderiv, dmxi, m_gp_temp, d_m_gp_temp_dT, d_m_gp_temp_dd);

  LINALG::Matrix<dim, 1> xgp;
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

  LINALG::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<GEN::pairedvector<int, double>> deriv_slave_normal(0, 0);
  std::vector<GEN::pairedvector<int, double>> deriv_master_normal(0, 0);
  sele.ComputeUnitNormalAtXi(sxi, slave_normal.A());
  mele.ComputeUnitNormalAtXi(mxi, master_normal.A());
  sele.DerivUnitNormalAtXi(sxi, deriv_slave_normal);
  mele.DerivUnitNormalAtXi(mxi, deriv_master_normal);

  double pen = ppn_;
  double pet = ppt_;
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  const LINALG::Matrix<dim, 1> contact_normal(gpn, true);
  double cauchy_nn_weighted_average = 0.;
  GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

  LINALG::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));

  LINALG::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);

  // variables for friction (declaration only)
  LINALG::Matrix<dim, 1> t1, t2;
  std::vector<GEN::pairedvector<int, double>> dt1, dt2;
  LINALG::Matrix<dim, 1> relVel;
  std::vector<GEN::pairedvector<int, double>> relVel_deriv(
      dim, sele.NumNode() * dim + mele.NumNode() * dim + dsxi[0].size() + dmxi[0].size());
  double vt1, vt2;
  GEN::pairedvector<int, double> dvt1(0);
  GEN::pairedvector<int, double> dvt2(0);
  double cauchy_nt1_weighted_average = 0.;
  GEN::pairedvector<int, double> cauchy_nt1_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  GEN::pairedvector<int, double> cauchy_nt1_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  LINALG::SerialDenseVector t1_adjoint_test_slave(sele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t1_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t1_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  LINALG::SerialDenseVector t1_adjoint_test_master(mele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t1_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t1_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  double cauchy_nt2_weighted_average = 0.;
  GEN::pairedvector<int, double> cauchy_nt2_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  GEN::pairedvector<int, double> cauchy_nt2_weighted_average_deriv_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  LINALG::SerialDenseVector t2_adjoint_test_slave(sele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t2_adjoint_test_slave(
      sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t2_adjoint_test_slave_T(
      sele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));
  LINALG::SerialDenseVector t2_adjoint_test_master(mele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t2_adjoint_test_master(
      mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_t2_adjoint_test_master_T(
      mele.ParentElement()->NumNode(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));
  double sigma_nt1_pen_vt1 = 0;
  double sigma_nt2_pen_vt2 = 0;
  GEN::pairedvector<int, double> d_sigma_nt1_pen_vt1(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt1_weighted_average_deriv.capacity() + dvt1.capacity(),
      0, 0);
  GEN::pairedvector<int, double> d_sigma_nt2_pen_vt2(
      dgapgp.capacity() + cauchy_nn_weighted_average_deriv.capacity() +
          cauchy_nt2_weighted_average_deriv.capacity() + dvt2.capacity(),
      0, 0);
  GEN::pairedvector<int, double> d_sigma_nt1_pen_vt1_T(
      sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
  GEN::pairedvector<int, double> d_sigma_nt2_pen_vt2_T(
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
  GEN::pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv.size() + dgapgp.size());
  for (_CI p = cauchy_nn_weighted_average_deriv.begin();
       p != cauchy_nn_weighted_average_deriv.end(); ++p)
    d_snn_av_pen_gap[p->first] += p->second;
  for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p)
    d_snn_av_pen_gap[p->first] += pen * p->second;

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
      GEN::pairedvector<int, double> d_fr_d(
          d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
      GEN::pairedvector<int, double> d_fr_T(cauchy_nn_weighted_average_deriv_T.size() +
                                            d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
      switch (frtype_)
      {
        case INPAR::CONTACT::friction_coulomb:
        {
          double fr_temp_fac = std::pow((std::max(s_gp_temp, m_gp_temp) - temp_damage_), 2.) /
                               std::pow((temp_damage_ - temp_ref_), 2.);
          fr = frcoeff_ * (-1.) * (snn_av_pen_gap)*fr_temp_fac;
          for (_CI p = d_snn_av_pen_gap.begin(); p != d_snn_av_pen_gap.end(); ++p)
            d_fr_d[p->first] += frcoeff_ * (-1.) * p->second * fr_temp_fac;
          for (_CI p = cauchy_nn_weighted_average_deriv_T.begin();
               p != cauchy_nn_weighted_average_deriv_T.end(); ++p)
            d_fr_T[p->first] += frcoeff_ * (-1.) * p->second * fr_temp_fac;
          if (s_gp_temp >= m_gp_temp)
          {
            for (_CI p = d_s_gp_temp_dd.begin(); p != d_s_gp_temp_dd.end(); ++p)
              d_fr_d[p->first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                  (s_gp_temp - temp_damage_) /
                                  std::pow((temp_damage_ - temp_ref_), 2.) * p->second;
            for (_CI p = d_s_gp_temp_dT.begin(); p != d_s_gp_temp_dT.end(); ++p)
              d_fr_T[p->first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                  (s_gp_temp - temp_damage_) /
                                  std::pow((temp_damage_ - temp_ref_), 2.) * p->second;
          }
          else
          {
            for (_CI p = d_m_gp_temp_dd.begin(); p != d_m_gp_temp_dd.end(); ++p)
              d_fr_d[p->first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                  (m_gp_temp - temp_damage_) /
                                  std::pow((temp_damage_ - temp_ref_), 2.) * p->second;
            for (_CI p = d_m_gp_temp_dT.begin(); p != d_m_gp_temp_dT.end(); ++p)
              d_fr_T[p->first] += frcoeff_ * (-1.) * (snn_av_pen_gap)*2. *
                                  (m_gp_temp - temp_damage_) /
                                  std::pow((temp_damage_ - temp_ref_), 2.) * p->second;
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
        for (_CI p = dvt1.begin(); p != dvt1.end(); ++p)
          d_sigma_nt1_pen_vt1[p->first] += pet * p->second;
        for (_CI p = cauchy_nt1_weighted_average_deriv.begin();
             p != cauchy_nt1_weighted_average_deriv.end(); ++p)
          d_sigma_nt1_pen_vt1[p->first] += p->second;
        for (_CI p = cauchy_nt1_weighted_average_deriv_T.begin();
             p != cauchy_nt1_weighted_average_deriv_T.end(); ++p)
          d_sigma_nt1_pen_vt1_T[p->first] += p->second;

        sigma_nt2_pen_vt2 = cauchy_nt2_weighted_average + pet * vt2;
        for (_CI p = dvt2.begin(); p != dvt2.end(); ++p)
          d_sigma_nt2_pen_vt2[p->first] += pet * p->second;
        for (_CI p = cauchy_nt2_weighted_average_deriv.begin();
             p != cauchy_nt2_weighted_average_deriv.end(); ++p)
          d_sigma_nt2_pen_vt2[p->first] += p->second;
        for (_CI p = cauchy_nt2_weighted_average_deriv_T.begin();
             p != cauchy_nt2_weighted_average_deriv_T.end(); ++p)
          d_sigma_nt2_pen_vt2_T[p->first] += p->second;
      }
      // slip
      else
      {
        GEN::pairedvector<int, double> tmp_d(
            dgapgp.size() + cauchy_nn_weighted_average_deriv.size() +
                cauchy_nt1_weighted_average_deriv.size() + dvt1.size(),
            0, 0);
        GEN::pairedvector<int, double> tmp_T(
            sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
        if (frtype_ == INPAR::CONTACT::friction_coulomb)
        {
          for (_CI p = d_fr_d.begin(); p != d_fr_d.end(); ++p)
            tmp_d[p->first] += p->second / tan_tr;
          for (_CI p = d_fr_T.begin(); p != d_fr_T.end(); ++p)
            tmp_T[p->first] += p->second / tan_tr;
        }
        for (_CI p = cauchy_nt1_weighted_average_deriv.begin();
             p != cauchy_nt1_weighted_average_deriv.end(); ++p)
          tmp_d[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt1_weighted_average + pet * vt1) * p->second;
        for (_CI p = dvt1.begin(); p != dvt1.end(); ++p)
          tmp_d[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt1_weighted_average + pet * vt1) * (+pet) * p->second;
        for (_CI p = cauchy_nt1_weighted_average_deriv_T.begin();
             p != cauchy_nt1_weighted_average_deriv_T.end(); ++p)
          tmp_T[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt1_weighted_average + pet * vt1) * p->second;

        for (_CI p = cauchy_nt2_weighted_average_deriv.begin();
             p != cauchy_nt2_weighted_average_deriv.end(); ++p)
          tmp_d[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt2_weighted_average + pet * vt2) * p->second;
        for (_CI p = dvt2.begin(); p != dvt2.end(); ++p)
          tmp_d[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt2_weighted_average + pet * vt2) * (+pet) * p->second;
        for (_CI p = cauchy_nt2_weighted_average_deriv_T.begin();
             p != cauchy_nt2_weighted_average_deriv_T.end(); ++p)
          tmp_T[p->first] += -fr / (tan_tr * tan_tr * tan_tr) *
                             (cauchy_nt2_weighted_average + pet * vt2) * p->second;

        sigma_nt1_pen_vt1 = fr / tan_tr * (cauchy_nt1_weighted_average + pet * vt1);
        for (_CI p = tmp_d.begin(); p != tmp_d.end(); ++p)
          d_sigma_nt1_pen_vt1[p->first] += p->second * (cauchy_nt1_weighted_average + pet * vt1);
        for (_CI p = cauchy_nt1_weighted_average_deriv.begin();
             p != cauchy_nt1_weighted_average_deriv.end(); ++p)
          d_sigma_nt1_pen_vt1[p->first] += fr / tan_tr * p->second;
        for (_CI p = dvt1.begin(); p != dvt1.end(); ++p)
          d_sigma_nt1_pen_vt1[p->first] += fr / tan_tr * pet * p->second;

        for (_CI p = tmp_T.begin(); p != tmp_T.end(); ++p)
          d_sigma_nt1_pen_vt1_T[p->first] += p->second * (cauchy_nt1_weighted_average + pet * vt1);
        for (_CI p = cauchy_nt1_weighted_average_deriv_T.begin();
             p != cauchy_nt1_weighted_average_deriv_T.end(); ++p)
          d_sigma_nt1_pen_vt1_T[p->first] += fr / tan_tr * p->second;

        sigma_nt2_pen_vt2 = fr / tan_tr * (cauchy_nt2_weighted_average + pet * vt2);
        for (_CI p = tmp_d.begin(); p != tmp_d.end(); ++p)
          d_sigma_nt2_pen_vt2[p->first] += p->second * (cauchy_nt2_weighted_average + pet * vt2);
        for (_CI p = cauchy_nt2_weighted_average_deriv.begin();
             p != cauchy_nt2_weighted_average_deriv.end(); ++p)
          d_sigma_nt2_pen_vt2[p->first] += fr / tan_tr * p->second;
        for (_CI p = dvt2.begin(); p != dvt2.end(); ++p)
          d_sigma_nt2_pen_vt2[p->first] += fr / tan_tr * pet * p->second;

        for (_CI p = tmp_T.begin(); p != tmp_T.end(); ++p)
          d_sigma_nt2_pen_vt2_T[p->first] += p->second * (cauchy_nt2_weighted_average + pet * vt2);
        for (_CI p = cauchy_nt2_weighted_average_deriv_T.begin();
             p != cauchy_nt2_weighted_average_deriv_T.end(); ++p)
          d_sigma_nt2_pen_vt2_T[p->first] += fr / tan_tr * p->second;
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
    GEN::pairedvector<int, double> d_diss_d(
        d_sigma_nt1_pen_vt1.size() + d_sigma_nt2_pen_vt2.size() + dvt1.size() + dvt2.size());
    GEN::pairedvector<int, double> d_diss_T(
        d_sigma_nt1_pen_vt1_T.size() + d_sigma_nt2_pen_vt2_T.size());
    if (frtype_)
    {
      diss = (sigma_nt1_pen_vt1 * vt1 + sigma_nt2_pen_vt2 * vt2) / dt_;

      for (_CI p = d_sigma_nt1_pen_vt1.begin(); p != d_sigma_nt1_pen_vt1.end(); ++p)
        d_diss_d[p->first] += vt1 * p->second / dt_;
      for (_CI p = d_sigma_nt2_pen_vt2.begin(); p != d_sigma_nt2_pen_vt2.end(); ++p)
        d_diss_d[p->first] += vt2 * p->second / dt_;
      for (_CI p = dvt1.begin(); p != dvt1.end(); ++p)
        d_diss_d[p->first] += sigma_nt1_pen_vt1 * p->second / dt_;
      for (_CI p = dvt2.begin(); p != dvt2.end(); ++p)
        d_diss_d[p->first] += sigma_nt2_pen_vt2 * p->second / dt_;
      for (_CI p = d_sigma_nt1_pen_vt1_T.begin(); p != d_sigma_nt1_pen_vt1_T.end(); ++p)
        d_diss_T[p->first] += vt1 * p->second / dt_;
      for (_CI p = d_sigma_nt2_pen_vt2_T.begin(); p != d_sigma_nt2_pen_vt2_T.end(); ++p)
        d_diss_T[p->first] += vt2 * p->second / dt_;
    }

    switch (nit_thr_)
    {
      case INPAR::CONTACT::NitThr_substitution:
      {
        const double beta_bar = beta * (-snn_av_pen_gap);
        const double q1 = beta_bar * (s_gp_temp - m_gp_temp);

        GEN::pairedvector<int, double> d_q1_d(
            d_snn_av_pen_gap.size() + d_s_gp_temp_dd.size() + d_m_gp_temp_dd.size());
        for (_CI p = d_snn_av_pen_gap.begin(); p != d_snn_av_pen_gap.end(); ++p)
          d_q1_d[p->first] += beta * (-p->second) * (s_gp_temp - m_gp_temp);
        for (_CI p = d_s_gp_temp_dd.begin(); p != d_s_gp_temp_dd.end(); ++p)
          d_q1_d[p->first] += beta_bar * p->second;
        for (_CI p = d_m_gp_temp_dd.begin(); p != d_m_gp_temp_dd.end(); ++p)
          d_q1_d[p->first] += beta_bar * (-p->second);

        GEN::pairedvector<int, double> d_q1_T(cauchy_nn_weighted_average_deriv_T.size() +
                                              d_s_gp_temp_dT.size() + d_m_gp_temp_dT.size());
        for (_CI p = cauchy_nn_weighted_average_deriv_T.begin();
             p != cauchy_nn_weighted_average_deriv_T.end(); ++p)
          d_q1_T[p->first] += beta * (-p->second) * (s_gp_temp - m_gp_temp);
        for (_CI p = d_s_gp_temp_dT.begin(); p != d_s_gp_temp_dT.end(); ++p)
          d_q1_T[p->first] += beta_bar * p->second;
        for (_CI p = d_m_gp_temp_dT.begin(); p != d_m_gp_temp_dT.end(); ++p)
          d_q1_T[p->first] += beta_bar * (-p->second);

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
        GEN::pairedvector<int, double> deriv_qn_weighted_average_d(
            sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
            deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
        GEN::pairedvector<int, double> deriv_qn_weighted_average_T(
            sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());
        LINALG::SerialDenseVector thermo_adjoint_test_slave(sele.ParentElement()->NumNode());
        LINALG::SerialDenseVector thermo_adjoint_test_master(mele.ParentElement()->NumNode());
        GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_thermo_adjoint_test_slave_d(
            sele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dsxi[0].size(), -1,
            LINALG::SerialDenseVector(sele.ParentElement()->NumNode(), true));
        GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_thermo_adjoint_test_master_d(
            mele.MoData().ParentDof().size() + deriv_contact_normal[0].size() + dmxi[0].size(), -1,
            LINALG::SerialDenseVector(mele.ParentElement()->NumNode(), true));
        GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_thermo_adjoint_test_slave_T(
            1, -1, LINALG::SerialDenseVector(sele.ParentElement()->NumNode(), true));
        GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_thermo_adjoint_test_master_T(
            1, -1, LINALG::SerialDenseVector(mele.ParentElement()->NumNode(), true));

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
          GEN::pairedvector<int, double> deriv_test_val_d(
              sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          GEN::pairedvector<int, double> deriv_test_val_T(
              sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

          const double beta_bar = beta * (-snn_av_pen_gap);
          test_val += +(beta_bar) / (pen_thermo + beta_bar) * qn_weighted_average;
          test_val += +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (s_gp_temp - m_gp_temp);
          test_val += -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * diss;

          for (_CI p = deriv_qn_weighted_average_d.begin(); p != deriv_qn_weighted_average_d.end();
               ++p)
            deriv_test_val_d[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = deriv_qn_weighted_average_T.begin(); p != deriv_qn_weighted_average_T.end();
               ++p)
            deriv_test_val_T[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * p->second;

          for (_CI p = d_s_gp_temp_dd.begin(); p != d_s_gp_temp_dd.end(); ++p)
            deriv_test_val_d[p->first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = d_m_gp_temp_dd.begin(); p != d_m_gp_temp_dd.end(); ++p)
            deriv_test_val_d[p->first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (-p->second);
          for (_CI p = d_s_gp_temp_dT.begin(); p != d_s_gp_temp_dT.end(); ++p)
            deriv_test_val_T[p->first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = d_m_gp_temp_dT.begin(); p != d_m_gp_temp_dT.end(); ++p)
            deriv_test_val_T[p->first] +=
                +(beta_bar * pen_thermo) / (pen_thermo + beta_bar) * (-p->second);

          for (_CI p = d_diss_d.begin(); p != d_diss_d.end(); ++p)
            deriv_test_val_d[p->first] +=
                -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p->second;
          for (_CI p = d_diss_T.begin(); p != d_diss_T.end(); ++p)
            deriv_test_val_T[p->first] +=
                -(beta_bar) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p->second;

          for (_CI p = d_snn_av_pen_gap.begin(); p != d_snn_av_pen_gap.end(); ++p)
            deriv_test_val_d[p->first] += beta * (-p->second) /
                                          (std::pow(pen_thermo + beta_bar, 2.)) *
                                          (+pen_thermo * qn_weighted_average +
                                              pen_thermo * pen_thermo * (s_gp_temp - m_gp_temp) -
                                              pen_thermo * (1. - delta_c - ws_thermo) * diss);
          for (_CI p = cauchy_nn_weighted_average_deriv_T.begin();
               p != cauchy_nn_weighted_average_deriv_T.end(); ++p)
            deriv_test_val_T[p->first] += beta * (-p->second) /
                                          (std::pow(pen_thermo + beta_bar, 2.)) *
                                          (+pen_thermo * qn_weighted_average +
                                              pen_thermo * pen_thermo * (s_gp_temp - m_gp_temp) -
                                              pen_thermo * (1. - delta_c - ws_thermo) * diss);

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
          GEN::pairedvector<int, double> deriv_test_val_d(
              sele.ParentElement()->NumNode() * dim + mele.ParentElement()->NumNode() * dim +
              deriv_contact_normal[0].size() + dsxi[0].size() + dmxi[0].size());
          GEN::pairedvector<int, double> deriv_test_val_T(
              sele.ParentElement()->NumNode() + mele.ParentElement()->NumNode());

          const double beta_bar = beta * (-snn_av_pen_gap);
          test_val += -(1.) / (pen_thermo + beta_bar) * qn_weighted_average;
          test_val += +(beta_bar) / (pen_thermo + beta_bar) * (s_gp_temp - m_gp_temp);
          test_val += +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * diss;

          for (_CI p = deriv_qn_weighted_average_d.begin(); p != deriv_qn_weighted_average_d.end();
               ++p)
            deriv_test_val_d[p->first] += -(1.) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = deriv_qn_weighted_average_T.begin(); p != deriv_qn_weighted_average_T.end();
               ++p)
            deriv_test_val_T[p->first] += -(1.) / (pen_thermo + beta_bar) * p->second;

          for (_CI p = d_s_gp_temp_dd.begin(); p != d_s_gp_temp_dd.end(); ++p)
            deriv_test_val_d[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = d_s_gp_temp_dT.begin(); p != d_s_gp_temp_dT.end(); ++p)
            deriv_test_val_T[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * p->second;
          for (_CI p = d_m_gp_temp_dd.begin(); p != d_m_gp_temp_dd.end(); ++p)
            deriv_test_val_d[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * (-p->second);
          for (_CI p = d_m_gp_temp_dT.begin(); p != d_m_gp_temp_dT.end(); ++p)
            deriv_test_val_T[p->first] += +(beta_bar) / (pen_thermo + beta_bar) * (-p->second);

          for (_CI p = d_diss_d.begin(); p != d_diss_d.end(); ++p)
            deriv_test_val_d[p->first] +=
                +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p->second;
          for (_CI p = d_diss_T.begin(); p != d_diss_T.end(); ++p)
            deriv_test_val_T[p->first] +=
                +(1.) / (pen_thermo + beta_bar) * (1. - delta_c - ws_thermo) * p->second;

          for (_CI p = d_snn_av_pen_gap.begin(); p != d_snn_av_pen_gap.end(); ++p)
            deriv_test_val_d[p->first] +=
                (+1. * qn_weighted_average + pen_thermo * (s_gp_temp - m_gp_temp) -
                    1. * (1. - delta_c - ws_thermo) * diss) /
                std::pow(pen_thermo + beta_bar, 2) * beta * (-p->second);

          for (_CI p = cauchy_nn_weighted_average_deriv_T.begin();
               p != cauchy_nn_weighted_average_deriv_T.end(); ++p)
            deriv_test_val_T[p->first] +=
                (+1. * qn_weighted_average + pen_thermo * (s_gp_temp - m_gp_temp) -
                    1. * (1. - delta_c - ws_thermo) * diss) /
                std::pow(pen_thermo + beta_bar, 2) * beta * (-p->second);

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

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType parentdistype, int dim>
void inline CONTACT::CoIntegratorNitscheTsi::SoEleGP(MORTAR::MortarElement& sele, const double wgt,
    const double* gpcoord, LINALG::Matrix<dim, 1>& pxsi, LINALG::Matrix<dim, dim>& derivtrafo)
{
  DRT::UTILS::CollectedGaussPoints intpoints =
      DRT::UTILS::CollectedGaussPoints(1);  // reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1], 0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(1, dim);
  derivtrafo.Clear();

  DRT::UTILS::BoundaryGPToParentGP<dim>(pqxg, derivtrafo, intpoints, sele.ParentElement()->Shape(),
      sele.Shape(), sele.FaceParentNumber());

  // coordinates of the current integration point in parent coordinate system
  for (int idim = 0; idim < dim; idim++) pxsi(idim) = pqxg(0, idim);
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SoEleCauchy(MORTAR::MortarElement& moEle,
    double* boundary_gpcoord, std::vector<GEN::pairedvector<int, double>> boundary_gpcoord_lin,
    const double gp_wgt, const LINALG::Matrix<dim, 1>& normal,
    std::vector<GEN::pairedvector<int, double>>& normal_deriv,
    const LINALG::Matrix<dim, 1>& direction,
    std::vector<GEN::pairedvector<int, double>>& direction_deriv, const double w, double& cauchy_nt,
    GEN::pairedvector<int, double>& deriv_sigma_nt_d,
    GEN::pairedvector<int, double>& deriv_sigma_nt_T, LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  LINALG::Matrix<dim, 1> pxsi(true);
  LINALG::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(moEle, boundary_gpcoord, gp_wgt, pxsi, derivtravo_slave);

  double sigma_nt;
  Epetra_SerialDenseMatrix dsntdd, d2sntdd2, d2sntDdDn, d2sntDdDt, d2sntDdDpxi, d2sntDdDT, dsntdT;
  LINALG::Matrix<dim, 1> dsntdn, dsntdt, dsntdpxi;
  dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())
      ->GetCauchyAtXi(pxsi, moEle.MoData().ParentDisp(), normal, direction, sigma_nt, &dsntdd,
          &d2sntdd2, &d2sntDdDn, &d2sntDdDt, &d2sntDdDpxi, &dsntdn, &dsntdt, &dsntdpxi,
          &moEle.MoData().ParentTemp(), &dsntdT, &d2sntDdDT);

  cauchy_nt += w * sigma_nt;

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
    deriv_sigma_nt_d[moEle.MoData().ParentDof().at(i)] += w * dsntdd(i, 0);

  for (int i = 0; i < dim - 1; ++i)
    for (GEN::pairedvector<int, double>::const_iterator p = boundary_gpcoord_lin[i].begin();
         p != boundary_gpcoord_lin[i].end(); ++p)
    {
      double& ref = deriv_sigma_nt_d[p->first];
      for (int k = 0; k < dim; ++k) ref += dsntdpxi(k) * derivtravo_slave(k, i) * p->second * w;
    }

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdn(d) * p->second * w;

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = direction_deriv[d].begin();
         p != direction_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdt(d) * p->second * w;

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

  return;
}
template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateTest(const double fac, MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int, double>>& dxi, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test_val,
    const GEN::pairedvector<int, double>& test_deriv_d,
    const GEN::pairedvector<int, double>& test_deriv_T, const LINALG::Matrix<dim, 1>& test_dir,
    const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv)
{
  CONTACT::CoIntegratorNitsche::IntegrateTest<dim>(fac, ele, shape, deriv, dxi, jac, jacintcellmap,
      wgt, test_val, test_deriv_d, test_dir, test_dir_deriv);

  for (GEN::pairedvector<int, double>::const_iterator p = test_deriv_T.begin();
       p != test_deriv_T.end(); ++p)
  {
    double* row = ele.GetNitscheContainer().k_dt(p->first);
    for (int s = 0; s < ele.NumNode(); ++s)
      for (int d = 0; d < Dim(); ++d)
        row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
                ele.ParentElement()->Shape(), ele.FaceParentNumber(), s) *
                dim +
            d] += fac * jac * wgt * test_dir(d) * p->second * shape(s);
  }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateAdjointTest(const double fac, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test,
    const GEN::pairedvector<int, double>& deriv_test_d,
    const GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  CONTACT::CoIntegratorNitsche::IntegrateAdjointTest<dim>(
      fac, jac, jacintcellmap, wgt, test, deriv_test_d, moEle, adjoint_test, deriv_adjoint_test_d);

  for (GEN::pairedvector<int, double>::const_iterator p = deriv_test_T.begin();
       p != deriv_test_T.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_dt(p->first), moEle.MoData().ParentDof().size())
        .Update(fac * jac * wgt * p->second, adjoint_test, 1.);

  for (GEN::pairedvector<int, LINALG::SerialDenseVector>::const_iterator p =
           deriv_adjoint_test_T.begin();
       p != deriv_adjoint_test_T.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_dt(p->first), moEle.MoData().ParentDof().size())
        .Update(fac * jac * wgt * test, p->second, 1.);
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateThermalTest(const double fac,
    MORTAR::MortarElement& ele, const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& deriv, const std::vector<GEN::pairedvector<int, double>>& dxi,
    const double jac, const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test_val, const GEN::pairedvector<int, double>& test_deriv_d,
    const GEN::pairedvector<int, double>& test_deriv_T)
{
  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  double val = fac * jac * wgt * test_val;

  for (int s = 0; s < ele.NumNode(); ++s)
    *(ele.GetNitscheContainer().rhs_t(DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), s))) += val * shape(s);

  GEN::pairedvector<int, double> val_deriv_d(jacintcellmap.size() + test_deriv_d.size());
  for (_CI p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
    val_deriv_d[p->first] += fac * p->second * wgt * test_val;
  for (_CI p = test_deriv_d.begin(); p != test_deriv_d.end(); ++p)
    val_deriv_d[p->first] += fac * jac * wgt * p->second;

  for (_CI p = val_deriv_d.begin(); p != val_deriv_d.end(); ++p)
  {
    double* row = ele.GetNitscheContainer().k_td(p->first);
    for (int s = 0; s < ele.NumNode(); ++s)
      row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), s)] += p->second * shape(s);
  }

  for (int e = 0; e < Dim() - 1; ++e)
    for (_CI p = dxi[e].begin(); p != dxi[e].end(); ++p)
    {
      double* row = ele.GetNitscheContainer().k_td(p->first);
      for (int s = 0; s < ele.NumNode(); ++s)
        row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
            ele.FaceParentNumber(), s)] += val * deriv(s, e) * p->second;
    }

  for (_CI p = test_deriv_T.begin(); p != test_deriv_T.end(); ++p)
  {
    double* row = ele.GetNitscheContainer().k_tt(p->first);
    for (int s = 0; s < ele.NumNode(); ++s)
      row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),
          ele.FaceParentNumber(), s)] += fac * jac * wgt * p->second * shape(s);
  }

  return;
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::IntegrateThermalAdjointTest(const double fac,
    const double jac, const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const double test, const GEN::pairedvector<int, double>& deriv_test_d,
    const GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (abs(fac) < 1.e-16) return;

  LINALG::SerialDenseVector(
      View, moEle.GetNitscheContainer().rhs_t(), moEle.ParentElement()->NumNode())
      .Update(fac * jac * wgt * test, adjoint_test, 1.);

  for (GEN::pairedvector<int, LINALG::SerialDenseVector>::const_iterator p =
           deriv_adjoint_test_d.begin();
       p != deriv_adjoint_test_d.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_td(p->first), moEle.ParentElement()->NumNode())
        .Update(fac * jac * wgt * test, p->second, 1.);

  for (GEN::pairedvector<int, double>::const_iterator p = jacintcellmap.begin();
       p != jacintcellmap.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_td(p->first), moEle.ParentElement()->NumNode())
        .Update(fac * p->second * wgt * test, adjoint_test, 1.);

  for (GEN::pairedvector<int, double>::const_iterator p = deriv_test_d.begin();
       p != deriv_test_d.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_td(p->first), moEle.ParentElement()->NumNode())
        .Update(fac * jac * wgt * p->second, adjoint_test, 1.);

  for (GEN::pairedvector<int, LINALG::SerialDenseVector>::const_iterator p =
           deriv_adjoint_test_T.begin();
       p != deriv_adjoint_test_T.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_tt(p->first), moEle.ParentElement()->NumNode())
        .Update(fac * jac * wgt * test, p->second, 1.);

  for (GEN::pairedvector<int, double>::const_iterator p = deriv_test_T.begin();
       p != deriv_test_T.end(); ++p)
    LINALG::SerialDenseVector(
        View, moEle.GetNitscheContainer().k_tt(p->first), moEle.ParentElement()->NumNode())
        .Update(fac * jac * wgt * p->second, adjoint_test, 1.);
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::BuildAdjointTestTsi(MORTAR::MortarElement& moEle,
    const double fac, const Epetra_SerialDenseMatrix& d2sntDdDT,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  if (moEle.MoData().ParentTempDof().size())
    for (int p = 0; p < moEle.ParentElement()->NumNode(); ++p)
    {
      LINALG::SerialDenseVector& at = deriv_adjoint_test_T[moEle.MoData().ParentTempDof()[p]];
      for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
        at(i) += fac * d2sntDdDT(i, p);
    }
}

template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SetupGpTemp(MORTAR::MortarElement& moEle,
    const LINALG::SerialDenseVector& val, const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int, double>> dxi, double& temp,
    GEN::pairedvector<int, double>& d_temp_dT, GEN::pairedvector<int, double>& d_temp_dd)
{
  if (moEle.MoData().ParentTemp().size() == 0)
  {
    temp = 0.;
    return;
  }

  LINALG::SerialDenseVector moele_temp(val.Length());
  for (int i = 0; i < moEle.NumNode(); ++i)
    moele_temp(i) =
        moEle.MoData().ParentTemp().at(DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
            moEle.ParentElement()->Shape(), moEle.FaceParentNumber(), i));
  temp = val.Dot(moele_temp);

  d_temp_dT.resize(val.Length());
  d_temp_dT.clear();
  for (int i = 0; i < moEle.NumNode(); ++i)
    d_temp_dT[moEle.MoData().ParentTempDof().at(DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        moEle.ParentElement()->Shape(), moEle.FaceParentNumber(), i))] = val(i);

  int deriv_size = 0.;
  for (int i = 0; i < dim - 1; ++i) deriv_size += dxi.at(i).size();
  d_temp_dd.resize(deriv_size);
  d_temp_dd.clear();
  for (int i = 0; i < dim - 1; ++i)
    for (GEN::pairedvector<int, double>::const_iterator p = dxi.at(i).begin(); p != dxi.at(i).end();
         ++p)
    {
      double& a = d_temp_dd[p->first];
      for (int n = 0; n < moEle.NumNode(); ++n) a += moele_temp(n) * deriv(n, i) * p->second;
    }
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::SoEleCauchyHeatflux(MORTAR::MortarElement& moEle,
    double* boundary_gpcoord,
    const std::vector<GEN::pairedvector<int, double>>& boundary_gpcoord_lin, const double gp_wgt,
    const LINALG::Matrix<dim, 1>& normal, std::vector<GEN::pairedvector<int, double>>& normal_deriv,
    const double w, double& heatflux, GEN::pairedvector<int, double>& dq_dd,
    GEN::pairedvector<int, double>& dq_dT, LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
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

  LINALG::Matrix<dim, 1> pxsi(true);
  LINALG::Matrix<dim, dim> derivtravo_slave;
  DRT::Element::DiscretizationType distype = moEle.ParentElement()->Shape();

  double q = 0;
  Epetra_SerialDenseMatrix dq_dT_ele, dq_dd_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi;
  LINALG::Matrix<dim, 1> dq_dn, dq_dpxi;

  switch (distype)
  {
    case DRT::Element::hex8:
    {
      DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<DRT::Element::hex8, dim>(moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case DRT::Element::hex27:
    {
      DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<DRT::Element::hex27, dim>(moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

      parent_ele->HeatFlux(moEle.MoData().ParentTemp(), moEle.MoData().ParentDisp(), pxsi, normal,
          q, &dq_dT_ele, &dq_dd_ele, &dq_dn, &dq_dpxi, &d2q_dT_dd, &d2q_dT_dn, &d2q_dT_dpxi);
      break;
    }
    case DRT::Element::tet4:
    {
      DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>* parent_ele =
          dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>*>(moEle.ParentElement());
      if (!parent_ele) dserror("thermo-mechanical Nitsche contact only for So3_Plast for now.");

      SoEleGP<DRT::Element::tet4, dim>(moEle, gp_wgt, boundary_gpcoord, pxsi, derivtravo_slave);

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
    for (GEN::pairedvector<int, double>::const_iterator p = boundary_gpcoord_lin[i].begin();
         p != boundary_gpcoord_lin[i].end(); ++p)
    {
      double& ref = dq_dd[p->first];
      for (int k = 0; k < dim; ++k) ref += dq_dpxi(k) * derivtravo_slave(k, i) * p->second * w;
    }

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      dq_dd[p->first] += dq_dn(d) * p->second * w;

  for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
    dq_dT[moEle.MoData().ParentTempDof().at(i)] += w * dq_dT_ele(i, 0);

  if (abs(theta_thermo_) > 1.e-12)
    BuildAdjointTestThermo<dim>(moEle, w, dq_dT_ele, d2q_dT_dd, d2q_dT_dn, d2q_dT_dpxi,
        normal_deriv, boundary_gpcoord_lin, derivtravo_slave, adjoint_test, deriv_adjoint_test_d,
        deriv_adjoint_test_T);
}


template <int dim>
void CONTACT::CoIntegratorNitscheTsi::BuildAdjointTestThermo(MORTAR::MortarElement& moEle,
    const double fac, const Epetra_SerialDenseMatrix& dq_dT_ele,
    const Epetra_SerialDenseMatrix& d2q_dT_dd, const Epetra_SerialDenseMatrix& d2q_dT_dn,
    const Epetra_SerialDenseMatrix& d2q_dT_dpxi,
    std::vector<GEN::pairedvector<int, double>>& normal_deriv,
    const std::vector<GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
    LINALG::Matrix<dim, dim>& derivtravo_slave, LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
    GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T)
{
  for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
    adjoint_test(i) = fac * dq_dT_ele(i, 0);

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
  {
    LINALG::SerialDenseVector& at = deriv_adjoint_test_d[moEle.MoData().ParentDof().at(i)];
    for (int j = 0; j < moEle.NumNode(); ++j) at(j) += fac * d2q_dT_dd(j, i);
  }

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
    {
      LINALG::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
        at(i) += fac * d2q_dT_dn(i, d) * p->second;
    }

  Epetra_SerialDenseMatrix tmp(moEle.ParentElement()->NumNode(), dim, false);
  Epetra_SerialDenseMatrix deriv_trafo(::View, derivtravo_slave.A(), derivtravo_slave.Rows(),
      derivtravo_slave.Rows(), derivtravo_slave.Columns());
  if (tmp.Multiply('N', 'N', 1., d2q_dT_dpxi, deriv_trafo, 0.)) dserror("multiply failed");
  for (int d = 0; d < dim - 1; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = boundary_gpcoord_lin[d].begin();
         p != boundary_gpcoord_lin[d].end(); ++p)
    {
      LINALG::SerialDenseVector& at = deriv_adjoint_test_d[p->first];
      for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
        at(i) += fac * tmp(i, d) * p->second;
    }
}
