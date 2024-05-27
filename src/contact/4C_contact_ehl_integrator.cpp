/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of ehl related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_ehl_integrator.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_nitsche_integrator.hpp"  // for CONTACT::UTILS:: functions
#include "4C_contact_node.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorEhl::integrate_gp_3_d(MORTAR::Element& sele, MORTAR::Element& mele,
    CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
    CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
    CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
    CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, CORE::GEN::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double& gap,
    CORE::GEN::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
    std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi)
{
  // check bound
  bool bound = false;
  for (int i = 0; i < sele.num_node(); ++i)
    if (dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->IsOnBoundorCE())
      FOUR_C_THROW("no boundary modification for EHL implemented");

  // is quadratic case?
  bool quad = sele.IsQuad();

  // weighted gap
  gp_3_d_w_gap(sele, sval, lmval, &gap, jac, wgt, quad);
  for (int j = 0; j < sele.num_node(); ++j)
    gp_g_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt, deriv_gap,
        derivjac, derivsxi, dualmap);

  // integrate D and M matrix
  gp_dm(sele, mele, lmval, sval, mval, jac, wgt, bound);
  gp_3_d_dm_lin(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, wgt, jac, derivsxi,
      derivmxi, derivjac, dualmap);

  // get second derivative of shape function
  CORE::LINALG::SerialDenseMatrix ssecderiv(sele.num_node(), 3);
  sele.evaluate2nd_deriv_shape(sxi, ssecderiv, sele.num_node());

  // weighted surface gradient
  gp_weighted_surf_grad_and_deriv(
      sele, sxi, derivsxi, lmval, lmderiv, dualmap, sval, sderiv, ssecderiv, wgt, jac, derivjac);

  //  // weighted tangential velocity (average and relative)
  gp_weighted_av_rel_vel(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
      derivjac, normal, dnmap_unit, gap, deriv_gap, sxi, mxi, derivsxi, derivmxi);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorEhl::integrate_gp_2_d(MORTAR::Element& sele, MORTAR::Element& mele,
    CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
    CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
    CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
    CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, CORE::GEN::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double& gap,
    CORE::GEN::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
    std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi)
{
  FOUR_C_THROW("2D EHL integration not supported");
}


void CONTACT::IntegratorEhl::gp_weighted_surf_grad_and_deriv(MORTAR::Element& sele,
    const double* xi, const std::vector<CORE::GEN::Pairedvector<int, double>>& dsxigp,
    const CORE::LINALG::SerialDenseVector& lmval, const CORE::LINALG::SerialDenseMatrix& lmderiv,
    const CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
    const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
    const CORE::LINALG::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
    const CORE::GEN::Pairedvector<int, double>& jacintcellmap)
{
  // empty local basis vectors
  std::vector<std::vector<double>> gxi(2, std::vector<double>(3, 0));

  // metrics routine gives local basis vectors
  sele.Metrics(xi, gxi.at(0).data(), gxi.at(1).data());

  if (Dim() == 2)
  {
    gxi.at(1).at(0) = gxi.at(1).at(1) = 0.;
    gxi.at(1).at(2) = 1.;
  }

  CORE::LINALG::Matrix<3, 1> t1(gxi.at(0).data(), true);
  CORE::LINALG::Matrix<3, 1> t2(gxi.at(1).data(), true);
  CORE::LINALG::Matrix<3, 1> n;
  n.CrossProduct(t1, t2);
  n.Scale(1. / n.Norm2());
  CORE::LINALG::Matrix<3, 3> covariant_metric;
  for (int i = 0; i < 3; ++i)
  {
    covariant_metric(i, 0) = t1(i);
    covariant_metric(i, 1) = t2(i);
    covariant_metric(i, 2) = n(i);
  }
  CORE::LINALG::Matrix<3, 3> contravariant_metric;
  contravariant_metric.Invert(covariant_metric);

  std::vector<std::vector<double>> gxi_contra(2, std::vector<double>(3, 0));
  for (int i = 0; i < Dim() - 1; ++i)
    for (int d = 0; d < Dim(); ++d) gxi_contra.at(i).at(d) = contravariant_metric(i, d);

  for (int a = 0; a < sele.num_node(); ++a)
  {
    DRT::Node* node = sele.Nodes()[a];
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("this is not a contact node");

    for (int c = 0; c < sele.num_node(); ++c)
    {
      CORE::LINALG::Matrix<3, 1>& tmp =
          cnode->EhlData().GetSurfGrad()[dynamic_cast<CONTACT::Node*>(sele.Nodes()[c])->Dofs()[0]];
      for (int d = 0; d < Dim(); ++d)
        for (int al = 0; al < Dim() - 1; ++al)
          tmp(d) += wgt * jac * lmval(a) * sderiv(c, al) * gxi_contra.at(al).at(d);
    }

    for (int c = 0; c < sele.num_node(); ++c)
    {
      for (auto p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      {
        CORE::LINALG::Matrix<3, 1>& tmp =
            cnode->EhlData()
                .GetSurfGradDeriv()[p->first]
                                   [dynamic_cast<CONTACT::Node*>(sele.Nodes()[c])->Dofs()[0]];
        for (int d = 0; d < Dim(); ++d)
          for (int al = 0; al < Dim() - 1; ++al)
            tmp(d) += wgt * p->second * lmval(a) * sderiv(c, al) * gxi_contra.at(al).at(d);
      }
      for (int e = 0; e < Dim() - 1; ++e)
        for (auto p = dsxigp.at(e).begin(); p != dsxigp.at(e).end(); ++p)
        {
          CORE::LINALG::Matrix<3, 1>& tmp =
              cnode->EhlData()
                  .GetSurfGradDeriv()[p->first]
                                     [dynamic_cast<CONTACT::Node*>(sele.Nodes()[c])->Dofs()[0]];
          for (int d = 0; d < Dim(); ++d)
            for (int al = 0; al < Dim() - 1; ++al)
              tmp(d) +=
                  wgt * jac * lmderiv(a, e) * p->second * sderiv(c, al) * gxi_contra.at(al).at(d);
        }

      for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
      {
        CORE::LINALG::Matrix<3, 1>& tmp =
            cnode->EhlData()
                .GetSurfGradDeriv()[p->first]
                                   [dynamic_cast<CONTACT::Node*>(sele.Nodes()[c])->Dofs()[0]];
        for (int d = 0; d < Dim(); ++d)
          for (int al = 0; al < Dim() - 1; ++al)
            for (int m = 0; m < sele.num_node(); ++m)
              tmp(d) +=
                  wgt * jac * p->second(a, m) * sval(m) * sderiv(c, al) * gxi_contra.at(al).at(d);
      }
    }
  }

  return;
}

void CONTACT::IntegratorEhl::gp_weighted_av_rel_vel(MORTAR::Element& sele, MORTAR::Element& mele,
    const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseVector& lmval,
    const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& sderiv,
    const CORE::LINALG::SerialDenseMatrix& mderiv, const CORE::LINALG::SerialDenseMatrix& lmderiv,
    const CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, const double& wgt,
    const double& jac, const CORE::GEN::Pairedvector<int, double>& derivjac, const double* normal,
    const std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, const double& gap,
    const CORE::GEN::Pairedvector<int, double>& deriv_gap, const double* sxi, const double* mxi,
    const std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
    const std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi)
{
  const int dim = 3;
  if (Dim() != dim)
    FOUR_C_THROW("dimension inconsistency, or is this not implemented for all spatial dimensions?");

  CORE::LINALG::Matrix<dim, 1> t1, t2;
  std::vector<CORE::GEN::Pairedvector<int, double>> dt1, dt2;
  CORE::LINALG::Matrix<dim, 1> relVel;
  std::vector<CORE::GEN::Pairedvector<int, double>> relVel_deriv(
      dim, sele.num_node() * dim + mele.num_node() * dim + derivsxi[0].size() + derivmxi[0].size());
  double vt1, vt2;
  CORE::GEN::Pairedvector<int, double> dvt1(0);
  CORE::GEN::Pairedvector<int, double> dvt2(0);

  CONTACT::UTILS::BuildTangentVectors<dim>(normal, dnmap_unit, t1.A(), dt1, t2.A(), dt2);
  CONTACT::UTILS::RelVelInvariant<dim>(sele, sxi, derivsxi, sval, sderiv, mele, mxi, derivmxi, mval,
      mderiv, gap, deriv_gap, relVel, relVel_deriv, -.5);

  CONTACT::UTILS::VectorScalarProduct<dim>(t1, dt1, relVel, relVel_deriv, vt1, dvt1);
  CONTACT::UTILS::VectorScalarProduct<dim>(t2, dt2, relVel, relVel_deriv, vt2, dvt2);

  for (int i = 0; i < sele.num_node(); ++i)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(sele.Nodes()[i]);
    for (int d = 0; d < dim; ++d)
      cnode->EhlData().get_weighted_rel_tang_vel()(d) +=
          jac * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));

    for (auto p = derivjac.begin(); p != derivjac.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp =
          cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        tmp(d) += p->second * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (int e = 0; e < dim - 1; ++e)
      for (auto p = derivsxi.at(e).begin(); p != derivsxi.at(e).end(); ++p)
      {
        CORE::LINALG::Matrix<3, 1>& tmp =
            cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first];
        for (int d = 0; d < dim; ++d)
          tmp(d) += jac * wgt * lmderiv(i, e) * p->second * (vt1 * t1(d) + vt2 * t2(d));
      }

    for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp =
          cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        for (int m = 0; m < sele.num_node(); ++m)
          tmp(d) += jac * wgt * p->second(i, m) * sval(m) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (auto p = dvt1.begin(); p != dvt1.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp =
          cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) += jac * wgt * lmval(i) * p->second * t1(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt1.at(d).begin(); p != dt1.at(d).end(); ++p)
        cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first](d) +=
            jac * wgt * lmval(i) * vt1 * p->second;

    for (auto p = dvt2.begin(); p != dvt2.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp =
          cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) += jac * wgt * lmval(i) * p->second * t2(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt2.at(d).begin(); p != dt2.at(d).end(); ++p)
        cnode->EhlData().get_weighted_rel_tang_vel_deriv()[p->first](d) +=
            jac * wgt * lmval(i) * vt2 * p->second;
  }

  CONTACT::UTILS::RelVel<dim>(sele, sval, sderiv, derivsxi, -1., relVel, relVel_deriv);
  CONTACT::UTILS::VectorScalarProduct<dim>(t1, dt1, relVel, relVel_deriv, vt1, dvt1);
  CONTACT::UTILS::VectorScalarProduct<dim>(t2, dt2, relVel, relVel_deriv, vt2, dvt2);

  for (int i = 0; i < sele.num_node(); ++i)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(sele.Nodes()[i]);
    for (int d = 0; d < dim; ++d)
      cnode->EhlData().get_weighted_av_tang_vel()(d) -=
          jac * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));

    for (auto p = derivjac.begin(); p != derivjac.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp = cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        tmp(d) -= p->second * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (int e = 0; e < dim - 1; ++e)
      for (auto p = derivsxi.at(e).begin(); p != derivsxi.at(e).end(); ++p)
      {
        CORE::LINALG::Matrix<3, 1>& tmp =
            cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first];
        for (int d = 0; d < dim; ++d)
          tmp(d) -= jac * wgt * lmderiv(i, e) * p->second * (vt1 * t1(d) + vt2 * t2(d));
      }

    for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp = cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        for (int m = 0; m < sele.num_node(); ++m)
          tmp(d) -= jac * wgt * p->second(i, m) * sval(m) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (auto p = dvt1.begin(); p != dvt1.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp = cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) -= jac * wgt * lmval(i) * p->second * t1(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt1.at(d).begin(); p != dt1.at(d).end(); ++p)
        cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first](d) -=
            jac * wgt * lmval(i) * vt1 * p->second;

    for (auto p = dvt2.begin(); p != dvt2.end(); ++p)
    {
      CORE::LINALG::Matrix<3, 1>& tmp = cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) -= jac * wgt * lmval(i) * p->second * t2(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt2.at(d).begin(); p != dt2.at(d).end(); ++p)
        cnode->EhlData().get_weighted_av_tang_vel_deriv()[p->first](d) -=
            jac * wgt * lmval(i) * vt2 * p->second;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
