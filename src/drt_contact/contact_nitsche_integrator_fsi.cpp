/*---------------------------------------------------------------------*/
/*! \file

\brief A class to perform integrations of nitsche related terms for the fsi contact case

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#include "contact_nitsche_integrator_fsi.H"

#include "contact_element.H"

#include "contact_node.H"

#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so3_poro.H"

#include "../drt_xfem/xfem_xfluid_contact_communicator.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoIntegratorNitscheFsi::CoIntegratorNitscheFsi(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm)
    : CoIntegratorNitsche(params, eletype, comm), ele_contact_state_(-2)
{
  if (fabs(theta_) > 1e-12)
    dserror("No Adjoint Consistency term for Nitsche Contact FSI implemented!");

  if (imortar_.isParameter("XFluid_Contact_Comm"))
    xf_c_comm_ = imortar_.get<Teuchos::RCP<XFEM::XFluid_Contact_Comm>>("XFluid_Contact_Comm");
  else
    dserror("Couldn't find XFluid_Contact_Comm!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitscheFsi::IntegrateDerivEle3D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  CONTACT::CoElement* csele = dynamic_cast<CONTACT::CoElement*>(&sele);
  if (!csele) dserror("Could cast to Contact Element!");

  // do quick orientation check
  LINALG::Matrix<3, 1> sn, mn;
  double center[2] = {0., 0.};
  sele.ComputeUnitNormalAtXi(center, sn.A());
  for (std::vector<MORTAR::MortarElement*>::iterator mit = meles.begin(); mit != meles.end(); ++mit)
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
void CONTACT::CoIntegratorNitscheFsi::IntegrateGP_3D(MORTAR::MortarElement& sele,
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
  // Here the consistent element normal is use to allow for a continous transition between FSI and
  // Contact
  double n[3];
  sele.ComputeUnitNormalAtXi(sxi, n);
  std::vector<GEN::pairedvector<int, double>> dn(3, sele.NumNode() * 3);
  dynamic_cast<CONTACT::CoElement&>(sele).DerivUnitNormalAtXi(sxi, dn);

  GPTS_forces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt,
      gap, deriv_gap, n, dn, sxi, mxi);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitscheFsi::GPTS_forces(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const LINALG::SerialDenseVector& sval,
    const LINALG::SerialDenseMatrix& sderiv,
    const std::vector<GEN::pairedvector<int, double>>& dsxi, const LINALG::SerialDenseVector& mval,
    const LINALG::SerialDenseMatrix& mderiv,
    const std::vector<GEN::pairedvector<int, double>>& dmxi, const double jac,
    const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const GEN::pairedvector<int, double>& dgapgp, double* gpn,
    std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi)
{
  // first rough check
  if (gap > 10 * std::max(sele.MaxEdgeSize(), mele.MaxEdgeSize())) return;

  const LINALG::Matrix<dim, 1> normal(gpn, true);

  if (dim != Dim()) dserror("dimension inconsistency");


  double pen = ppn_;
  double pet = ppt_;

  double ws = 0.;
  double wm = 0.;
  CONTACT::UTILS::NitscheWeightsAndScaling(sele, mele, nit_wgt_, dt_, ws, wm, pen, pet);


  double normal_contact_transition = 0.;
  bool FSI_integrated = true;  // bool indicates if fsi condition is already evaluated ... --> if
                               // true no contribution here ...

  LINALG::Matrix<dim, 1> pxsi(true);
  LINALG::Matrix<dim, dim> derivtravo_slave;
  CONTACT::UTILS::MapGPtoParent<dim>(sele, sxi, wgt, pxsi, derivtravo_slave);

  bool gp_on_this_proc;

  normal_contact_transition = xf_c_comm_->Get_FSI_Traction(
      &sele, pxsi, LINALG::Matrix<dim - 1, 1>(sxi, false), normal, FSI_integrated, gp_on_this_proc);

  {
    LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
    xf_c_comm_->Gmsh_Write(sgp_x, gp_on_this_proc, 7);
  }

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
               LINALG::Matrix<dim - 1, 1>(sxi, true), normal, normal) +
      wm * CONTACT::UTILS::SolidCauchyAtXi(dynamic_cast<CONTACT::CoElement*>(&mele),
               LINALG::Matrix<dim - 1, 1>(mxi, true), normal, normal) +
      pen * gap;

  {
    LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
    xf_c_comm_->Gmsh_Write(sgp_x, snn_pengap, 4);
    xf_c_comm_->Gmsh_Write(sgp_x, normal_contact_transition, 5);
  }


  if (snn_pengap >= normal_contact_transition && !FSI_integrated)
  {
    GEN::pairedvector<int, double> lin_fluid_traction(0);
    IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        normal_contact_transition, lin_fluid_traction, normal, dnmap_unit);

    {
      LINALG::Matrix<3, 1> sgp_x;
      for (int i = 0; i < sele.NumNode(); ++i)
        for (int d = 0; d < dim; ++d)
          sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];
      xf_c_comm_->Gmsh_Write(sgp_x, normal_contact_transition, 0);
      xf_c_comm_->Gmsh_Write(sgp_x, 2.0, 2);
    }
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


  typedef GEN::pairedvector<int, double>::const_iterator _CI;

  double cauchy_nn_weighted_average = 0.;
  GEN::pairedvector<int, double> cauchy_nn_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());

  LINALG::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave(
      sele.MoData().ParentDof().size() + dnmap_unit[0].size() + dsxi[0].size(), -1,
      LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));

  LINALG::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
  GEN::pairedvector<int, LINALG::SerialDenseVector> deriv_normal_adjoint_test_master(
      mele.MoData().ParentDof().size() + dnmap_unit[0].size() + dmxi[0].size(), -1,
      LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, normal, dnmap_unit, normal, dnmap_unit, ws,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv, normal_adjoint_test_slave,
      deriv_normal_adjoint_test_slave);
  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, normal, dnmap_unit, normal, dnmap_unit, wm,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv, normal_adjoint_test_master,
      deriv_normal_adjoint_test_master);

  double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  GEN::pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv.size() + dgapgp.size());
  for (_CI p = cauchy_nn_weighted_average_deriv.begin();
       p != cauchy_nn_weighted_average_deriv.end(); ++p)
    d_snn_av_pen_gap[p->first] += p->second;
  for (_CI p = dgapgp.begin(); p != dgapgp.end(); ++p)
    d_snn_av_pen_gap[p->first] += pen * p->second;

  // test in normal contact direction
  IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, snn_av_pen_gap,
      d_snn_av_pen_gap, normal, dnmap_unit);

  UpdateEleContactState(sele, 1);

  {
    LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[i])->xspatial()[d];

    xf_c_comm_->Gmsh_Write(sgp_x, snn_av_pen_gap, 0);
    xf_c_comm_->Gmsh_Write(sgp_x, 1.0, 2);
    xf_c_comm_->Inc_GP(0);
  }

  return;
}

void CONTACT::CoIntegratorNitscheFsi::UpdateEleContactState(MORTAR::MortarElement& sele, int state)
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

double CONTACT::UTILS::SolidCauchyAtXi(CONTACT::CoElement* cele, const LINALG::Matrix<2, 1>& xsi,
    const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& dir)
{
  if (cele->ParentElement()->Shape() != DRT::Element::hex8)
    dserror("This Element shape is not implemented for CONTACT::UTILS::CauchyStressatXi");

  LINALG::Matrix<3, 1> pxsi(true);
  LINALG::Matrix<3, 3> trafo;
  CONTACT::UTILS::SoEleGP<DRT::Element::hex8, 3>(*cele, 1., xsi.A(), pxsi, trafo);

  double sigma_nt;

  if (!cele->MoData().ParentPFPres().size())
    dynamic_cast<DRT::ELEMENTS::So_base*>(cele->ParentElement())
        ->GetCauchyAtXi(pxsi, cele->MoData().ParentDisp(), n, dir, sigma_nt, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL);
  else
    dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
        cele->ParentElement())
        ->GetCauchyAtXi(
            pxsi, cele->MoData().ParentDisp(), cele->MoData().ParentPFPres(), n, dir, sigma_nt);
  return sigma_nt;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType parentdistype, int dim>
void inline CONTACT::UTILS::SoEleGP(MORTAR::MortarElement& sele, const double wgt,
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
