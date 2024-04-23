/*---------------------------------------------------------------------*/
/*! \file

\brief A class to perform integrations of nitsche related terms for the fsi contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_nitsche_integrator_fsi.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_node.hpp"
#include "4C_discretization_fem_general_utils_boundary_integration.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_poro.hpp"
#include "4C_xfem_xfluid_contact_communicator.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::IntegratorNitscheFsi::IntegratorNitscheFsi(
    Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
    : IntegratorNitsche(params, eletype, comm), ele_contact_state_(-2)
{
  if (fabs(theta_) > 1e-12)
    FOUR_C_THROW("No Adjoint Consistency term for Nitsche Contact FSI implemented!");

  if (imortar_.isParameter("XFluidContactComm"))
    xf_c_comm_ = imortar_.get<Teuchos::RCP<XFEM::XFluidContactComm>>("XFluidContactComm");
  else
    FOUR_C_THROW("Couldn't find XFluidContactComm!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorNitscheFsi::IntegrateDerivEle3D(MORTAR::Element& sele,
    std::vector<MORTAR::Element*> meles, bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  auto* csele = dynamic_cast<CONTACT::Element*>(&sele);
  if (!csele) FOUR_C_THROW("Could cast to Contact Element!");

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
  CONTACT::Integrator::IntegrateDerivEle3D(sele, meles, boundary_ele, proj_, comm, cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorNitscheFsi::IntegrateGP_3D(MORTAR::Element& sele, MORTAR::Element& mele,
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
  // Here the consistent element normal is use to allow for a continous transition between FSI and
  // Contact
  double n[3];
  sele.ComputeUnitNormalAtXi(sxi, n);
  std::vector<CORE::GEN::Pairedvector<int, double>> dn(3, sele.NumNode() * 3);
  dynamic_cast<CONTACT::Element&>(sele).DerivUnitNormalAtXi(sxi, dn);

  GPTSForces<3>(sele, mele, sval, sderiv, derivsxi, mval, mderiv, derivmxi, jac, derivjac, wgt, gap,
      deriv_gap, n, dn, sxi, mxi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::IntegratorNitscheFsi::GPTSForces(MORTAR::Element& sele, MORTAR::Element& mele,
    const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
    const std::vector<CORE::GEN::Pairedvector<int, double>>& dsxi,
    const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
    const std::vector<CORE::GEN::Pairedvector<int, double>>& dmxi, const double jac,
    const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
    const CORE::GEN::Pairedvector<int, double>& dgapgp, const double* gpn,
    std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi)
{
  // first rough check
  if (gap > 10 * std::max(sele.MaxEdgeSize(), mele.MaxEdgeSize())) return;

  const CORE::LINALG::Matrix<dim, 1> normal(gpn, true);

  if (dim != Dim()) FOUR_C_THROW("dimension inconsistency");


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

  double normal_contact_transition = xf_c_comm_->Get_FSI_Traction(&sele, pxsi,
      CORE::LINALG::Matrix<dim - 1, 1>(sxi, false), normal, FSI_integrated, gp_on_this_proc);
#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->xspatial()[d];
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
      ws * CONTACT::UTILS::SolidCauchyAtXi(dynamic_cast<CONTACT::Element*>(&sele),
               CORE::LINALG::Matrix<dim - 1, 1>(sxi, true), normal, normal) +
      wm * CONTACT::UTILS::SolidCauchyAtXi(dynamic_cast<CONTACT::Element*>(&mele),
               CORE::LINALG::Matrix<dim - 1, 1>(mxi, true), normal, normal) +
      pen * gap;
#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->xspatial()[d];
    xf_c_comm_->Gmsh_Write(sgp_x, snn_pengap, 4);
    xf_c_comm_->Gmsh_Write(sgp_x, normal_contact_transition, 5);
  }
#endif


  if (snn_pengap >= normal_contact_transition && !FSI_integrated)
  {
    CORE::GEN::Pairedvector<int, double> lin_fluid_traction(0);
    IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt,
        normal_contact_transition, lin_fluid_traction, normal, dnmap_unit);
#ifdef WRITE_GMSH
    {
      CORE::LINALG::Matrix<3, 1> sgp_x;
      for (int i = 0; i < sele.NumNode(); ++i)
        for (int d = 0; d < dim; ++d)
          sgp_x(d) += sval(i) * dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->xspatial()[d];
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
  CORE::GEN::Pairedvector<int, double> cauchy_nn_weighted_average_deriv(
      sele.NumNode() * 3 * 12 + sele.MoData().ParentDisp().size() +
      mele.MoData().ParentDisp().size());

  CORE::LINALG::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
  CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave(
      sele.MoData().ParentDof().size() + dnmap_unit[0].size() + dsxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(sele.MoData().ParentDof().size(), true));

  CORE::LINALG::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
  CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector> deriv_normal_adjoint_test_master(
      mele.MoData().ParentDof().size() + dnmap_unit[0].size() + dmxi[0].size(), -1,
      CORE::LINALG::SerialDenseVector(mele.MoData().ParentDof().size(), true));

  SoEleCauchy<dim>(sele, sxi, dsxi, wgt, normal, dnmap_unit, normal, dnmap_unit, ws,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv, normal_adjoint_test_slave,
      deriv_normal_adjoint_test_slave);
  SoEleCauchy<dim>(mele, mxi, dmxi, wgt, normal, dnmap_unit, normal, dnmap_unit, wm,
      cauchy_nn_weighted_average, cauchy_nn_weighted_average_deriv, normal_adjoint_test_master,
      deriv_normal_adjoint_test_master);

  double snn_av_pen_gap = cauchy_nn_weighted_average + pen * gap;
  CORE::GEN::Pairedvector<int, double> d_snn_av_pen_gap(
      cauchy_nn_weighted_average_deriv.size() + dgapgp.size());
  for (const auto& p : cauchy_nn_weighted_average_deriv) d_snn_av_pen_gap[p.first] += p.second;
  for (const auto& p : dgapgp) d_snn_av_pen_gap[p.first] += pen * p.second;

  // test in normal contact direction
  IntegrateTest<dim>(-1., sele, sval, sderiv, dsxi, jac, jacintcellmap, wgt, snn_av_pen_gap,
      d_snn_av_pen_gap, normal, dnmap_unit);

  UpdateEleContactState(sele, 1);
#ifdef WRITE_GMSH
  {
    CORE::LINALG::Matrix<3, 1> sgp_x;
    for (int i = 0; i < sele.NumNode(); ++i)
      for (int d = 0; d < dim; ++d)
        sgp_x(d) += sval(i) * dynamic_cast<CONTACT::Node*>(sele.Nodes()[i])->xspatial()[d];

    xf_c_comm_->Gmsh_Write(sgp_x, snn_av_pen_gap, 0);
    xf_c_comm_->Gmsh_Write(sgp_x, 1.0, 2);
  }
#endif
  xf_c_comm_->Inc_GP(0);
}

void CONTACT::IntegratorNitscheFsi::UpdateEleContactState(MORTAR::Element& sele, int state)
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

double CONTACT::UTILS::SolidCauchyAtXi(CONTACT::Element* cele,
    const CORE::LINALG::Matrix<2, 1>& xsi, const CORE::LINALG::Matrix<3, 1>& n,
    const CORE::LINALG::Matrix<3, 1>& dir)
{
  if (cele->ParentElement()->Shape() != CORE::FE::CellType::hex8)
    FOUR_C_THROW("This Element shape is not implemented for CONTACT::UTILS::CauchyStressatXi");

  CORE::LINALG::Matrix<3, 1> pxsi(true);
  CORE::LINALG::Matrix<3, 3> trafo;
  CONTACT::UTILS::SoEleGP<CORE::FE::CellType::hex8, 3>(*cele, 1., xsi.A(), pxsi, trafo);

  double sigma_nt;

  if (!cele->MoData().ParentPFPres().size())
  {
    dynamic_cast<DRT::ELEMENTS::SoBase*>(cele->ParentElement())
        ->GetCauchyNDirAndDerivativesAtXi(pxsi, cele->MoData().ParentDisp(), n, dir, sigma_nt,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr);
  }
  else
  {
    dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
        cele->ParentElement())
        ->GetCauchyNDirAndDerivativesAtXi(pxsi, cele->MoData().ParentDisp(),
            cele->MoData().ParentPFPres(), n, dir, sigma_nt, nullptr, nullptr, nullptr, nullptr,
            nullptr);
  }
  return sigma_nt;
}


FOUR_C_NAMESPACE_CLOSE
