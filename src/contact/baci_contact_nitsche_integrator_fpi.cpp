/*---------------------------------------------------------------------*/
/*! \file

\brief A class to perform integrations of nitsche related terms for the fpi contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#include "baci_contact_nitsche_integrator_fpi.H"

#include "baci_contact_element.H"
#include "baci_contact_nitsche_integrator_fsi.H"
#include "baci_contact_node.H"
#include "baci_xfem_xfluid_contact_communicator.H"
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoIntegratorNitscheFpi::CoIntegratorNitscheFpi(
    Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
    : CoIntegratorNitschePoro(params, eletype, comm), ele_contact_state_(-2)
{
  if (imortar_.isParameter("XFluid_Contact_Comm"))
    xf_c_comm_ = imortar_.get<Teuchos::RCP<XFEM::XFluid_Contact_Comm>>("XFluid_Contact_Comm");
  else
    dserror("Couldn't find XFluid_Contact_Comm!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheFpi::IntegrateDerivEle3D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  auto* csele = dynamic_cast<CONTACT::CoElement*>(&sele);
  if (!csele) dserror("Could cast to Contact Element!");

  // do quick orientation check
  CORE::LINALG::Matrix<3, 1> sn, mn;
  double center[2] = {0., 0.};
  sele.ComputeUnitNormalAtXi(center, sn.A());
  for (auto mit = meles.begin(); mit != meles.end(); ++mit)
  {
    (*mit)->ComputeUnitNormalAtXi(center, mn.A());
    if (sn.Dot(mn) > -1e-1)
    {
      meles.erase(mit);
      --mit;
    }
  }

  if (!meles.size()) return;

  if (xf_c_comm_->HigherIntegrationforContactElement(sele.Id()))
    xf_c_comm_->GetCutSideIntegrationPoints(sele.Id(), coords_, weights_, ngp_);

  // Call Base Contact Integratederiv with potentially increased number of GPs!
  CONTACT::CoIntegrator::IntegrateDerivEle3D(sele, meles, boundary_ele, proj_, comm, cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheFpi::IntegrateGP_3D(MORTAR::MortarElement& sele,
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
  // Here the consistent element normal is use to allow for a continous transition between FSI and
  // Contact
  double n[3];
  sele.ComputeUnitNormalAtXi(sxi, n);
  std::vector<CORE::GEN::pairedvector<int, double>> dn(3, sele.NumNode() * 3);
  dynamic_cast<CONTACT::CoElement&>(sele).DerivUnitNormalAtXi(sxi, dn);

  GPTSForces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt, gap,
      deriv_gap, n, dn, sxi, mxi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheFpi::GPTSForces(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const CORE::LINALG::SerialDenseVector& sval,
    const CORE::LINALG::SerialDenseMatrix& sderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
    const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
    const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, const double jac,
    const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
    std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi)
{
  // first rough check
  if (gap > 10 * std::max(sele.MaxEdgeSize(), mele.MaxEdgeSize())) return;

  const CORE::LINALG::Matrix<dim, 1> normal(gpn, true);

  if (dim != Dim()) dserror("dimension inconsistency");


  double pen = ppn_;
  double pet = ppt_;

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);

  bool FSI_integrated = true;  // bool indicates if fsi condition is already evaluated ... --> if
                               // true no contribution here ...

  CORE::LINALG::Matrix<dim, 1> pxsi(true);
  CORE::LINALG::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(sele, sxi, wgt, pxsi, derivtravo_slave);

  bool gp_on_this_proc;
  double normal_contact_transition = GetNormalContactTransition<dim>(
      sele, mele, sval, mval, sxi, pxsi, normal, FSI_integrated, gp_on_this_proc);

#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
    xf_c_comm_->Gmsh_Write(sgp_x, gp_on_this_proc, 7);
  }
#endif

  if (!gp_on_this_proc) return;

  static int processed_gps = 0;
  ++processed_gps;
  if (processed_gps == 100000)
  {
    std::cout << "==| Processed again 100000 C-Gps! (" << Comm_.MyPID() << ") |==" << std::endl;
    processed_gps = 0;
  }


  // fast check
  const double snn_pengap =
      ws * CONTACT::UTILS::SolidCauchyAtXi(dynamic_cast<CONTACT::CoElement*>(&sele),
               CORE::LINALG::Matrix<dim - 1, 1>(sxi, true), normal, normal) +
      wm * CONTACT::UTILS::SolidCauchyAtXi(dynamic_cast<CONTACT::CoElement*>(&mele),
               CORE::LINALG::Matrix<dim - 1, 1>(mxi, true), normal, normal) +
      pen * gap;

#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
    xf_c_comm_->Gmsh_Write(sgp_x, snn_pengap, 4);
    xf_c_comm_->Gmsh_Write(sgp_x, normal_contact_transition, 5);
  }
#endif


  if (!FSI_integrated || (gap < (1 + xf_c_comm_->Get_fpi_pcontact_fullfraction()) *
                                     xf_c_comm_->Get_fpi_pcontact_exchange_dist() &&
                             xf_c_comm_->Get_fpi_pcontact_exchange_dist() > 1e-16))
  {
    double ffac = 1;
    if (gap < (1 + xf_c_comm_->Get_fpi_pcontact_fullfraction()) *
                  xf_c_comm_->Get_fpi_pcontact_exchange_dist() &&
        FSI_integrated &&
        gap > xf_c_comm_->Get_fpi_pcontact_fullfraction() *
                  xf_c_comm_->Get_fpi_pcontact_exchange_dist())
      ffac = 1. - (gap / (xf_c_comm_->Get_fpi_pcontact_exchange_dist()) -
                      xf_c_comm_->Get_fpi_pcontact_fullfraction());

    IntegratePoroNoOutFlow<dim>(
        -ffac, sele, sxi, sval, sderiv, jac, jacintcellmap, wgt, normal, dnmap_unit, mele, mval);
#ifdef WRITE_GMSH
    {
      CORE::LINALG::Matrix<3, 1> sgp_x;
      for (int i = 0; i < sele.NumNode(); ++i)
        for (int d = 0; d < dim; ++d)
          sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
      xf_c_comm_->Gmsh_Write(sgp_x, ffac, 8);
    }
#endif
  }
  else
  {
#ifdef WRITE_GMSH
    {
      CORE::LINALG::Matrix<3, 1> sgp_x;
      for (int i = 0; i < sele.NumNode(); ++i)
        for (int d = 0; d < dim; ++d)
          sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
      xf_c_comm_->Gmsh_Write(sgp_x, gap, 9);
    }
#endif
  }


  if (snn_pengap >= normal_contact_transition && !FSI_integrated)
  {
    CORE::GEN::pairedvector<int, double> lin_fluid_traction(0);
    IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        normal_contact_transition, lin_fluid_traction, lin_fluid_traction, normal, dnmap_unit);
#ifdef WRITE_GMSH
    {
      CORE::LINALG::Matrix<3, 1> sgp_x;
      for (int i = 0; i < sele.NumNode(); ++i)
        for (int d = 0; d < dim; ++d)
          sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
      xf_c_comm_->Gmsh_Write(sgp_x, normal_contact_transition, 0);
      xf_c_comm_->Gmsh_Write(sgp_x, 2.0, 2);
    }
#endif
    UpdateEleContactState(sele, 0);
  }

  if (snn_pengap >= normal_contact_transition)
  {
    UpdateEleContactState(sele, -1);
    if (!FSI_integrated)
      xf_c_comm_->Inc_GP(1);
    else
      xf_c_comm_->Inc_GP(2);
    return;
  }

  double cauchy_nn_weighted_average = 0.;
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_d(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());
  CORE::GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv_p(
      sele.MoData().ParentPFPres().size() + mele.MoData().ParentPFPres().size());

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, normal, dnmap_unit, normal, dnmap_unit, ws,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);
  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, normal, dnmap_unit, normal, dnmap_unit, wm,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv_d,
      cauchy_nn_weighted_average_deriv_p);

  const double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  CORE::GEN::pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv_d.size() + dgapgp.size());
  for (const auto& p : cauchy_nn_weighted_average_deriv_d) d_snn_av_pen_gap[p.first] += p.second;
  for (const auto& p : dgapgp) d_snn_av_pen_gap[p.first] += pen * p.second;

  // test in normal contact direction
  IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, snn_av_pen_gap,
      d_snn_av_pen_gap, cauchy_nn_weighted_average_deriv_p, normal, dnmap_unit);

  UpdateEleContactState(sele, 1);
#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];

    xf_c_comm_->Gmsh_Write(sgp_x, snn_av_pen_gap, 0);
    xf_c_comm_->Gmsh_Write(sgp_x, 1.0, 2);
  }
#endif
  xf_c_comm_->Inc_GP(0);
}

void CONTACT::CoIntegratorNitscheFpi::UpdateEleContactState(MORTAR::MortarElement& sele, int state)
{
  if (!state && ele_contact_state_)
  {
    ele_contact_state_ = state;
    xf_c_comm_->RegisterContactElementforHigherIntegration(sele.Id());
  }
  else if (ele_contact_state_ == -2)
    ele_contact_state_ = state;
  else if (ele_contact_state_ == -state)  // switch between contact and no contact
  {
    ele_contact_state_ = 0;
    xf_c_comm_->RegisterContactElementforHigherIntegration(sele.Id());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
double CONTACT::CoIntegratorNitscheFpi::GetNormalContactTransition(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const CORE::LINALG::SerialDenseVector& sval,
    const CORE::LINALG::SerialDenseVector& mval, const double* sxi,
    const CORE::LINALG::Matrix<dim, 1>& pxsi, const CORE::LINALG::Matrix<dim, 1>& normal,
    bool& FSI_integrated, bool& gp_on_this_proc)
{
  double poropressure(0.0);
  if (GetPoroPressure(sele, sval, mele, mval, poropressure))
  {
    return xf_c_comm_->Get_FSI_Traction(&sele, pxsi, CORE::LINALG::Matrix<dim - 1, 1>(sxi, false),
        normal, FSI_integrated, gp_on_this_proc, &poropressure);
  }
  else
    return xf_c_comm_->Get_FSI_Traction(&sele, pxsi, CORE::LINALG::Matrix<dim - 1, 1>(sxi, false),
        normal, FSI_integrated, gp_on_this_proc);
}
