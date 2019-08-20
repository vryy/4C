/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poro contact case

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#include "contact_nitsche_integrator_poro.H"
#include "contact_integrator.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_paramsinterface.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_so3/so_base.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so3_poro.H"

#include "../drt_mat/elasthyper.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include "../linalg/linalg_utils.H"
#include "contact_nitsche_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_mat/structporo.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoIntegratorNitschePoro::CoIntegratorNitschePoro(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm)
    : CoIntegratorNitsche(params, eletype, comm),
      no_penetration_(params.get<bool>("CONTACTNOPEN")),
      dv_dd_(params.get<double>("porotimefac"))
{
  if (fabs(theta_) > 1e-16) dserror("Poro Contact just implemented Adjoint free ...");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitschePoro::IntegrateGP_3D(MORTAR::MortarElement& sele,
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
  // TEUCHOS_FUNC_TIME_MONITOR("CONTACT::CoIntegratorNitsche::IntegrateGP_3D");
  // We use the consistent element normal for poro contact!
  // if (nit_normal_==INPAR::CONTACT::NitNor_ele)
  {
    double n[3];
    sele.ComputeUnitNormalAtXi(sxi, n);
    std::vector<GEN::pairedvector<int, double>> dn(3, sele.NumNode() * 3);
    dynamic_cast<CONTACT::CoElement&>(sele).DerivUnitNormalAtXi(sxi, dn);

    GPTS_forces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
        gap, deriv_gap, n, dn, sxi, mxi);
  }
  //  else if (nit_normal_==INPAR::CONTACT::NitNor_sm)
  //    dserror("Want to use the element normal!");

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitschePoro::IntegrateGP_2D(MORTAR::MortarElement& sele,
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
  dserror("2D is not implemented!");
  //    GPTS_forces<2>(sele,mele,sval,sderiv,derivsxi,mval,mderiv,derivmxi,
  //        jac,derivjac,wgt,gap,deriv_gap,normal,dnmap_unit,sxi,mxi);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitschePoro::GPTS_forces(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const LINALG::SerialDenseVector& sval,
    const LINALG::SerialDenseMatrix& sderiv,
    const std::vector<GEN::pairedvector<int, double>>& dsxi, const LINALG::SerialDenseVector& mval,
    const LINALG::SerialDenseMatrix& mderiv,
    const std::vector<GEN::pairedvector<int, double>>& dmxi, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const GEN::pairedvector<int, double>& dgapgp, double* gpn,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi)
{
  if (sele.Owner() != Comm_.MyPID()) return;

  static const bool do_fast_checks = true;
  // first rough check
  if (do_fast_checks)
    if (abs(theta_) < 1.e-12 && gap > std::max(sele.MaxEdgeSize(), mele.MaxEdgeSize())) return;

  const LINALG::Matrix<dim, 1> normal(gpn, true);

  if (dim != Dim()) dserror("dimension inconsistency");

  LINALG::Matrix<dim, 1> slave_normal, master_normal;
  std::vector<GEN::pairedvector<int, double>> deriv_slave_normal(0, 0);
  std::vector<GEN::pairedvector<int, double>> deriv_master_normal(0, 0);
  sele.ComputeUnitNormalAtXi(sxi, slave_normal.A());
  mele.ComputeUnitNormalAtXi(mxi, master_normal.A());
  sele.DerivUnitNormalAtXi(sxi, deriv_slave_normal);
  mele.DerivUnitNormalAtXi(mxi, deriv_master_normal);

  double pen = ppn_;
  double pet = ppt_;

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);


  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  double cauchy_nn_weighted_average = 0.;
  GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_d(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_p(
      sele.MoData().ParentPFPres().size() + mele.MoData().ParentPFPres().size());

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, slave_normal, deriv_slave_normal, normal, dnmap_unit, ws,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);
  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, master_normal, deriv_master_normal, normal, dnmap_unit,
      -wm, cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);

  const double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  GEN::pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv_d.size() + dgapgp.size());
  for (_CI p = cauchy_nn_weighted_average_deriv_d.begin();
       p != cauchy_nn_weighted_average_deriv_d.end(); ++p)
    d_snn_av_pen_gap[p->first] += p->second;
  for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p)
    d_snn_av_pen_gap[p->first] += pen * p->second;

  if (snn_av_pen_gap < 0.)
  {
    // test in normal contact direction
    IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, snn_av_pen_gap,
        d_snn_av_pen_gap, cauchy_nn_weighted_average_deriv_p, normal, dnmap_unit);
    IntegrateTest<dim>(+1., mele, mval, mderiv, dmxi, jac, jacintcellmap, wgt, snn_av_pen_gap,
        d_snn_av_pen_gap, cauchy_nn_weighted_average_deriv_p, normal, dnmap_unit);

    IntegratePoroNoOutFlow<dim>(
        -1, sele, sxi, sval, sderiv, jac, jacintcellmap, wgt, normal, dnmap_unit, mele, mval);
    IntegratePoroNoOutFlow<dim>(
        +1, mele, mxi, mval, mderiv, jac, jacintcellmap, wgt, normal, dnmap_unit, sele, sval);
  }
  return;
}


template <int dim>
void CONTACT::CoIntegratorNitschePoro::SoEleCauchy(MORTAR::MortarElement& moEle,
    double* boundary_gpcoord, std::vector<GEN::pairedvector<int, double>> boundary_gpcoord_lin,
    const double gp_wgt, const LINALG::Matrix<dim, 1>& normal,
    std::vector<GEN::pairedvector<int, double>>& normal_deriv,
    const LINALG::Matrix<dim, 1>& direction,
    std::vector<GEN::pairedvector<int, double>>& direction_deriv, const double w, double& cauchy_nt,
    GEN::pairedvector<int, double>& deriv_sigma_nt_d,
    GEN::pairedvector<int, double>& deriv_sigma_nt_p)
{
  LINALG::Matrix<dim, 1> pxsi(true);
  LINALG::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(moEle, boundary_gpcoord, gp_wgt, pxsi, derivtravo_slave);

  double sigma_nt;
  Epetra_SerialDenseMatrix dsdd;
  Epetra_SerialDenseMatrix dsntdd, dsntdp;
  LINALG::Matrix<dim, 1> dsntdn, dsntdt, dsntdpxi;

  if (!moEle.MoData().ParentPFPres().size())
  {
    dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())
        ->GetCauchyAtXi(pxsi, moEle.MoData().ParentDisp(), normal, direction, sigma_nt, &dsntdd,
            NULL, NULL, NULL, NULL, &dsntdn, &dsntdt, &dsntdpxi);
  }
  else
  {
    dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
        moEle.ParentElement())
        ->GetCauchyAtXi(pxsi, moEle.MoData().ParentDisp(), moEle.MoData().ParentPFPres(), normal,
            direction, sigma_nt, &dsntdd, &dsntdp, &dsntdn, &dsntdt, &dsntdpxi);
  }

  cauchy_nt += w * sigma_nt;

  for (int i = 0; i < moEle.ParentElement()->NumNode() * dim; ++i)
    deriv_sigma_nt_d[moEle.MoData().ParentDof().at(i)] += w * dsntdd(i, 0);

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
         p != normal_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdn(d) * p->second * w;

  for (int d = 0; d < dim; ++d)
    for (GEN::pairedvector<int, double>::const_iterator p = direction_deriv[d].begin();
         p != direction_deriv[d].end(); ++p)
      deriv_sigma_nt_d[p->first] += dsntdt(d) * p->second * w;

  if (moEle.MoData().ParentPFPres().size())
  {
    for (int i = 0; i < moEle.ParentElement()->NumNode(); ++i)
      deriv_sigma_nt_p[moEle.MoData().ParentPFDof()[i * (dim + 1) + 3]] += dsntdp(i, 0) * w;
  }

  return;
}

template <int dim>
void CONTACT::CoIntegratorNitschePoro::IntegrateTest(const double fac, MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int, double>>& dxi, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test_val,
    const GEN::pairedvector<int, double>& test_deriv_d,
    const GEN::pairedvector<int, double>& test_deriv_p, const LINALG::Matrix<dim, 1>& test_dir,
    const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv)
{
  if (abs(fac) < 1.e-16) return;

  CONTACT::CoIntegratorNitsche::IntegrateTest<dim>(fac, ele, shape, deriv, dxi, jac, jacintcellmap,
      wgt, test_val, test_deriv_d, test_dir, test_dir_deriv);

  for (GEN::pairedvector<int, double>::const_iterator p = test_deriv_p.begin();
       p != test_deriv_p.end(); ++p)
  {
    double* row = ele.GetNitscheContainer().k_dp(p->first);
    for (int s = 0; s < ele.NumNode(); ++s)
      for (int d = 0; d < Dim(); ++d)
        row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
                ele.ParentElement()->Shape(), ele.FaceParentNumber(), s) *
                dim +
            d] -= fac * jac * wgt * test_dir(d) * p->second * shape(s);
  }
}

template <int dim>
void CONTACT::CoIntegratorNitschePoro::IntegratePoroNoOutFlow(const double fac,
    MORTAR::MortarElement& ele, double* xi, const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& deriv, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
    const LINALG::Matrix<dim, 1>& normal,
    const std::vector<GEN::pairedvector<int, double>>& normal_deriv,
    MORTAR::MortarElement& otherele, const LINALG::SerialDenseVector& othershape)
{
  if (abs(fac) < 1e-16) return;
  if (!no_penetration_) return;

  if (!ele.MoData().ParentPFDof().size()) return;

  // weighting for poro pressure depenent if two or onesided porocontact
  double sweight = 1;
  double oweight = 0;

  double spresgp = 0;
  double srelveln = 0;
  for (int j = 0; j < ele.NumNode(); ++j)
  {
    int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), j);
    spresgp += ele.MoData().ParentPFPres()[pj] * shape(j);
    for (int d = 0; d < dim; ++d)
    {
      srelveln +=
          (ele.MoData().ParentPFVel()[pj * dim + d] - ele.MoData().ParentVel()[pj * dim + d]) *
          shape(j) * normal(d);
    }
  }

  double sporosity = -1;
  double sJ = -1;
  std::map<int, double> sJLin;
  double sdphi_dp;
  double sdphi_dJ;
  GetPoroQuantitiesatGP(ele, xi, spresgp, sJ, sJLin, sporosity, sdphi_dp, sdphi_dJ);

  double opresgp = 0;
  if (otherele.MoData().ParentPFDof().size())  // two sided poro contact case
  {
    sweight = 0.5;
    oweight = 1 - sweight;
    for (int j = 0; j < otherele.NumNode(); ++j)
    {
      int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          otherele.ParentElement()->Shape(), otherele.FaceParentNumber(), j);
      opresgp += otherele.MoData().ParentPFPres()[pj] * othershape(j);
    }
  }

  double val = fac * jac * wgt / dv_dd_;  //*1./dv_dd_;
  for (int i = 0; i < ele.NumNode(); ++i)
  {
    int pi = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
        ele.ParentElement()->Shape(), ele.FaceParentNumber(), i);
    for (int j = 0; j < ele.NumNode(); ++j)
    {
      int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), j);
      for (int d = 0; d < dim; ++d)
      {
        (*ele.GetNitscheContainer().rhs_p(pi * (dim + 1) + d)) +=
            shape(i) * ele.MoData().ParentPFPres()[pj] * shape(j) * normal(d) * val *
            sweight;  //(v,k1 p1 n)
        ele.GetNitscheContainer().k_pp(
            ele.MoData().ParentPFDof()[pj * (dim + 1) + dim])[(pi * (dim + 1) + d)] -=
            shape(i) * shape(j) * normal(d) * val * sweight;  //(v,k1 dp1/dp n)
        for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
             p != normal_deriv[d].end(); ++p)
        {
          ele.GetNitscheContainer().k_pd(p->first)[(pi * (dim + 1) + d)] -=
              shape(i) * shape(j) * ele.MoData().ParentPFPres()[pj] * p->second * val *
              sweight;  //(v,k1 p1 dn/dd)

          ele.GetNitscheContainer().k_pd(p->first)[(pi * (dim + 1) + dim)] +=
              shape(i) * sporosity * p->second *
              (ele.MoData().ParentPFVel()[pj * dim + d] - ele.MoData().ParentVel()[pj * dim + d]) *
              shape(j) * val * sweight;  // (k1 q1, phi (vF-vS) dn/dd)
        }
        ele.GetNitscheContainer().k_pd(ele.MoData().ParentDof()[pj * dim + d])[(
            pi * (dim + 1) + dim)] -=  // (k1 q1, phi (vF-dvS/dd) n)
            shape(i) * sporosity * normal(d) * shape(j) * val * dv_dd_ * sweight;

        ele.GetNitscheContainer().k_pp(ele.MoData().ParentPFDof()[pj * (dim + 1) + d])[(
            pi * (dim + 1) + dim)] +=  // (k1 q1, phi (dvF/dvF-vS) n)
            shape(i) * sporosity * normal(d) * shape(j) * val * sweight;

        (*ele.GetNitscheContainer().rhs_p(pi * (dim + 1) + (dim))) -=
            shape(i) * sporosity * normal(d) *
            (ele.MoData().ParentPFVel()[pj * dim + d] - ele.MoData().ParentVel()[pj * dim + d]) *
            shape(j) * val * sweight;  // (k1 q1, phi (vF-vS) n)
      }
      ele.GetNitscheContainer().k_pp(
          ele.MoData().ParentPFDof()[pj * (dim + 1) + dim])[(pi * (dim + 1) + dim)] +=
          val * srelveln * shape(i) * sdphi_dp * shape(j) * sweight;  // (k1 q1, phi/dp (vF-vS) n)
    }

    for (std::map<int, double>::iterator dJit = sJLin.begin(); dJit != sJLin.end(); ++dJit)
      ele.GetNitscheContainer().k_pd(dJit->first)[(pi * (dim + 1) + dim)] +=
          val * srelveln * shape(i) * sdphi_dJ * dJit->second *
          sweight;  // (k1 q1, phi/dd (vF-vS) n)
  }

  if (oweight > 1e-16)
  {
    for (int i = 0; i < ele.NumNode(); ++i)
    {
      int pi = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), i);
      for (int j = 0; j < otherele.NumNode(); ++j)
      {
        int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
            otherele.ParentElement()->Shape(), otherele.FaceParentNumber(), j);
        for (int d = 0; d < dim; ++d)
        {
          (*ele.GetNitscheContainer().rhs_p(pi * (dim + 1) + d)) +=
              shape(i) * otherele.MoData().ParentPFPres()[pj] * othershape(j) * normal(d) * val *
              oweight;  //(v,k2 p2 n)
          ele.GetNitscheContainer().k_pp(
              otherele.MoData().ParentPFDof()[pj * (dim + 1) + dim])[(pi * (dim + 1) + d)] -=
              shape(i) * othershape(j) * normal(d) * val * oweight;  //(v,k2 dp2/dp n)
          for (GEN::pairedvector<int, double>::const_iterator p = normal_deriv[d].begin();
               p != normal_deriv[d].end(); ++p)
          {
            ele.GetNitscheContainer().k_pd(p->first)[(pi * (dim + 1) + d)] -=
                shape(i) * othershape(j) * otherele.MoData().ParentPFPres()[pj] * p->second * val *
                oweight;  //(v,k2 p2 dn/dd)

            otherele.GetNitscheContainer().k_pd(p->first)[(pj * (dim + 1) + dim)] +=
                othershape(j) * sporosity * p->second *
                (ele.MoData().ParentPFVel()[pi * dim + d] -
                    ele.MoData().ParentVel()[pi * dim + d]) *
                shape(i) * val * oweight;  // (k2 q2, phi (vF-vS) dn/dd)
          }

          otherele.GetNitscheContainer().k_pd(ele.MoData().ParentDof()[pi * dim + d])[(
              pj * (dim + 1) + dim)] -=  // (k2 q2, phi (vF-dvS/dd) n)
              othershape(j) * sporosity * normal(d) * shape(i) * val * dv_dd_ * oweight;

          otherele.GetNitscheContainer().k_pp(ele.MoData().ParentPFDof()[pi * (dim + 1) + d])[(
              pj * (dim + 1) + dim)] +=  // (k2 q2, phi (dvF/dvF-vS) n)
              othershape(j) * sporosity * normal(d) * shape(i) * val * oweight;

          (*otherele.GetNitscheContainer().rhs_p(pj * (dim + 1) + (dim))) -=
              othershape(j) * sporosity * normal(d) *
              (ele.MoData().ParentPFVel()[pi * dim + d] - ele.MoData().ParentVel()[pi * dim + d]) *
              shape(i) * val * oweight;  // (k2 q2, phi (vF-vS) n)
        }
        otherele.GetNitscheContainer().k_pp(
            ele.MoData().ParentPFDof()[pi * (dim + 1) + dim])[(pj * (dim + 1) + dim)] +=
            val * srelveln * othershape(j) * sdphi_dp * shape(i) *
            oweight;  // (k2 q2, phi/dp (vF-vS) n)
      }
    }
    for (int j = 0; j < otherele.NumNode(); ++j)
    {
      int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          otherele.ParentElement()->Shape(), otherele.FaceParentNumber(), j);
      for (std::map<int, double>::iterator dJit = sJLin.begin(); dJit != sJLin.end(); ++dJit)
        otherele.GetNitscheContainer().k_pd(dJit->first)[(pj * (dim + 1) + dim)] +=
            val * srelveln * othershape(j) * sdphi_dJ * dJit->second *
            oweight;  // (k2 q2, phi/dd (vF-vS) n)
    }
  }
}

bool CONTACT::CoIntegratorNitschePoro::GetPoroPressure(MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape, MORTAR::MortarElement& otherele,
    const LINALG::SerialDenseVector& othershape, double& poropressure)
{
  if (!ele.MoData().ParentPFDof().size() && !otherele.MoData().ParentPFDof().size()) return false;

  double w1 = 1;
  if (ele.MoData().ParentPFDof().size() && otherele.MoData().ParentPFDof().size())
    w1 = 0.5;
  else if (ele.MoData().ParentPFDof().size())
    w1 = 1.0;
  else if (otherele.MoData().ParentPFDof().size())
    w1 = 0.0;
  else
    dserror("Thats not exptected...!");
  double w2 = 1.0 - w1;

  poropressure = 0.0;
  if (ele.MoData().ParentPFDof().size())
  {
    for (int j = 0; j < ele.NumNode(); ++j)
    {
      int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          ele.ParentElement()->Shape(), ele.FaceParentNumber(), j);
      poropressure += w1 * ele.MoData().ParentPFPres()[pj] * shape(j);
    }
  }

  if (otherele.MoData().ParentPFDof().size())
  {
    for (int j = 0; j < otherele.NumNode(); ++j)
    {
      int pj = DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(
          otherele.ParentElement()->Shape(), otherele.FaceParentNumber(), j);
      poropressure += w2 * otherele.MoData().ParentPFPres()[pj] * othershape(j);
    }
  }
  return true;
}


void CONTACT::CoIntegratorNitschePoro::GetPoroQuantitiesatGP(MORTAR::MortarElement& ele, double* xi,
    double& spresgp,  //(in)
    double& sJ, std::map<int, double>& sJLin, double& sporosity, double& sdphi_dp,
    double& sdphi_dJ)  // out
{
  static double dummy = 1.0;
  sJ = DetDeformationGradient(ele, dummy, xi, sJLin);
  Teuchos::ParameterList sparams;  // empty parameter list;

  Teuchos::RCP<MAT::StructPoro> sstructmat =
      Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele.ParentElement()->Material(0));
  if (sstructmat == Teuchos::null)
    sstructmat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele.ParentElement()->Material(1));
  if (sstructmat == Teuchos::null) dserror("Cast to StructPoro failed!");
  sstructmat->ComputeSurfPorosity(sparams, spresgp, sJ, ele.FaceParentNumber(), 1, sporosity,
      &sdphi_dp, &sdphi_dJ, NULL, NULL, NULL, false);
}

template void CONTACT::CoIntegratorNitschePoro::IntegrateTest<2>(const double,
    MORTAR::MortarElement&, const LINALG::SerialDenseVector&, const LINALG::SerialDenseMatrix&,
    const std::vector<GEN::pairedvector<int, double>>& i, const double,
    const GEN::pairedvector<int, double>&, const double, const double,
    const GEN::pairedvector<int, double>&, const GEN::pairedvector<int, double>&,
    const LINALG::Matrix<2, 1>& test_dir,
    const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv);
template void CONTACT::CoIntegratorNitschePoro::IntegrateTest<3>(const double,
    MORTAR::MortarElement&, const LINALG::SerialDenseVector&, const LINALG::SerialDenseMatrix&,
    const std::vector<GEN::pairedvector<int, double>>& i, const double,
    const GEN::pairedvector<int, double>&, const double, const double,
    const GEN::pairedvector<int, double>&, const GEN::pairedvector<int, double>&,
    const LINALG::Matrix<3, 1>& test_dir,
    const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv);
