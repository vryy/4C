/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.cpp

\brief Classes for interface coupling in the XFEM

<pre>
Maintainer: Raffaela Kruse /Benedikt Schott
            kruse@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include <Teuchos_TimeMonitor.hpp>

namespace DRT
{
namespace ELEMENTS
{
namespace XFLUID
{

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::AddSlaveEleDisp(
  const DRT::Discretization &  slavedis,     ///< coupling slave discretization
  const std::vector<int>&      lm            ///< local map
)
{
  // leave, if displacements are not set
  if (!slavedis.HasState(disp_statename_))
    return;
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(disp_statename_);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", disp_statename_.c_str());

  // extract local values of the global vector
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (unsigned inode=0; inode<slave_nen_; ++inode)  // number of nodes
  {
    for (unsigned idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (slave_disp_)(idim,inode) = mymatrix[idim+(inode*slave_numdof)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
    }
  }

  // add the displacement of the interface
  for (unsigned inode = 0; inode < slave_nen_; ++inode)
  {
    slave_xyze_(0,inode) += slave_disp_(0, inode);
    slave_xyze_(1,inode) += slave_disp_(1, inode);
    slave_xyze_(2,inode) += slave_disp_(2, inode);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::SetSlaveState(
  const DRT::Discretization &  slavedis,       ///< coupling slave discretization
  const std::vector<int>&      lm              ///< local map
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(vel_statename_);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", vel_statename_.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (unsigned inode=0; inode<slave_nen_; ++inode)  // number of nodes
  {
    for(unsigned idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (slave_vel_)(idim,inode) = mymatrix[idim+(inode*slave_numdof)];  // state vector includes velocity and pressure
    }
    if (slave_numdof == nsd_+1) (slave_pres_)(inode,0) = mymatrix[nsd_+(inode*slave_numdof)];
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::SetSlaveStaten(
  const DRT::Discretization &  slavedis,       ///< coupling slave discretization
  const std::vector<int>&      lm              ///< local map
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(veln_statename_);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", veln_statename_.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (unsigned inode=0; inode<slave_nen_; ++inode)  // number of nodes
  {
    for(unsigned idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (slave_veln_)(idim,inode) = mymatrix[idim+(inode*slave_numdof)];  // state vector includes velocity and pressure
    }
    if (slave_numdof == nsd_+1) (slave_presn_)(inode,0) = mymatrix[nsd_+(inode*slave_numdof)];
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceVelnp(
  LINALG::Matrix<nsd_,1>& ivelint  ///< interface velocity at coupling slave side
) const
{
  ivelint.Multiply(slave_vel_,slave_funct_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceVeln(
  LINALG::Matrix<nsd_,1>& ivelintn  ///< interface velocity at coupling slave side
) const
{
  ivelintn.Multiply(slave_veln_,slave_funct_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfacePresnp(
  double & ipres ///< interface pressure at coupling slave side
) const
{
  // pressure at current gauss-point
  ipres = slave_funct_.Dot(slave_pres_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfacePresn(
  double & ipresn ///< interface pressure at coupling slave side
) const
{
  // pressure at current gauss-point
  ipresn = slave_funct_.Dot(slave_presn_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceVelGradnp(
  LINALG::Matrix<nsd_,nsd_>& velgradint  ///< interface velocity gradients at coupling slave side
) const
{
  velgradint = slave_vderxy_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceVelGradn(
  LINALG::Matrix<nsd_,nsd_>& velgradintn  ///< interface velocity gradients at coupling slave side
) const
{
  velgradintn = slave_vderxyn_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetSlaveFunct(
  LINALG::Matrix<slave_nen_,1> & slave_funct ///< coupling slave shape functions
  ) const
{
  slave_funct = slave_funct_;
}

 /*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::SetInterfaceJumpStatenp(
  const DRT::Discretization &  cutterdis,      ///< cutter discretization
  const std::string            state,          ///< state
  const std::vector<int>&      lm              ///< local map
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = cutterdis.GetState(state);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (unsigned inode=0; inode<slave_nen_; ++inode)  // number of nodes
  {
    for(unsigned idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (interface_velnp_jump_)(idim,inode) = mymatrix[idim+(inode*slave_numdof)];  // state vector includes velocity and pressure
      // no pressure jump required
    }
  }

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::SetInterfaceJumpStaten(
  const DRT::Discretization &  cutterdis,      ///< cutter discretization
  const std::string            state,          ///< state
  const std::vector<int>&      lm              ///< local map
)
{
// get state of the global vector
Teuchos::RCP<const Epetra_Vector> matrix_state = cutterdis.GetState(state);
if (matrix_state == Teuchos::null)
  dserror("Cannot get state vector %s", state.c_str());

// extract local values of the global vectors
std::vector<double> mymatrix(lm.size());
DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

for (unsigned inode=0; inode<slave_nen_; ++inode)  // number of nodes
{
  for(unsigned idim=0; idim<nsd_; ++idim) // number of dimensions
  {
    (interface_veln_jump_)(idim,inode) = mymatrix[idim+(inode*slave_numdof)];  // state vector includes velocity and pressure
    // no pressure jump required
  }
}

return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceJumpVelnp(
  LINALG::Matrix<nsd_,1> & ivelint_jump ///< cutter element interface velocity jump or prescribed DBC at Gaussian point
  ) const
{
  ivelint_jump.Multiply(interface_velnp_jump_, slave_funct_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetInterfaceJumpVeln(
  LINALG::Matrix<nsd_,1> & ivelintn_jump ///< cutter element interface velocity jump or prescribed DBC at Gaussian point
  ) const
{
  ivelintn_jump.Multiply(interface_veln_jump_, slave_funct_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::Evaluate( LINALG::Matrix<nsd_,1> & xslave )
{
  // coupling with a 2D element
  if (slave_nsd_ == nsd_ - 1)
  {
    // evaluate shape function at solution
    DRT::UTILS::shape_function_2D( slave_funct_, xslave( 0 ), xslave( 1 ), slave_distype );
//    dserror("You called 3D evaluation routine when coupling with a 2D element.");
    return;
  }

  // find element local position of gauss point
  GEO::CUT::Position<slave_distype> pos( slave_xyze_, xslave );
  pos.Compute();

  const LINALG::Matrix<nsd_,1> rst_slave = pos.LocalCoordinates();

  DRT::UTILS::shape_function_3D( slave_funct_, rst_slave( 0 ), rst_slave( 1 ), rst_slave( 2 ), slave_distype );
  DRT::UTILS::shape_function_3D_deriv1( slave_deriv_, rst_slave( 0 ), rst_slave( 1 ), rst_slave( 2 ), slave_distype );

  LINALG::Matrix<nsd_,nsd_> slave_xjm(true);
  LINALG::Matrix<nsd_,nsd_> slave_xji(true);

  slave_xjm.MultiplyNT(slave_deriv_,slave_xyze_);
  slave_xji.Invert(slave_xjm);

  // compute global first derivates
  slave_derxy_.Multiply(slave_xji,slave_deriv_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  slave_vderxy_.MultiplyNT(slave_vel_,slave_derxy_);
  // previous time step
  slave_vderxyn_.MultiplyNT(slave_veln_,slave_derxy_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::ComputeCharElementLength( double & h_k )
{
  // only for 3d elements
  if (slave_nsd_ != nsd_)
    dserror("Computation of characteristic element length only possible for a 3d slave element!");

  //const double vol = VolumeViaNumIntegration<DISTYPE>(emb_xyze_);
  const unsigned numnode = DRT::UTILS::DisTypeToNumNodePerEle<slave_distype>::numNodePerElement;

  // use one point integration rule to calculate hk at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili = DRT::UTILS::intrule3D_undefined;

  switch (slave_distype)
  {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
        integrationrule_stabili = DRT::UTILS::intrule_hex_1point;
      break;
    case DRT::Element::tet4:
    case DRT::Element::tet10:
        integrationrule_stabili = DRT::UTILS::intrule_tet_1point;
      break;
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
        integrationrule_stabili = DRT::UTILS::intrule_wedge_1point;
      break;
    case DRT::Element::pyramid5:
        integrationrule_stabili = DRT::UTILS::intrule_pyramid_1point;
      break;
    default:
      dserror("invalid discretization type for fluid3");
      integrationrule_stabili = DRT::UTILS::intrule3D_undefined;
      break;
  }

  // integration points
  const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  LINALG::Matrix<3,1> e;
  e(0) = intpoints.qxg[0][0];
  e(1) = intpoints.qxg[0][1];
  e(2) = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  LINALG::Matrix<3,numnode> deriv;
  DRT::UTILS::shape_function_3D_deriv1(deriv, e(0), e(1), e(2), slave_distype);

  // get Jacobian matrix and determinant
  // xjm_ = deriv_(i,k)*xyze(j,k);
  LINALG::Matrix<3,3> xjm;
  xjm.MultiplyNT(deriv,slave_xyze_);

  const double vol = wquad * xjm.Determinant();

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  h_k = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::ComputeInterfaceForce(
  Epetra_SerialDenseVector &   iforce,     ///< interface force vector
  LINALG::Matrix<nsd_,1> &     traction,   ///< traction vector at gaussian point
  const double &               fac         ///< integration factor
)
{
  for (unsigned inode = 0; inode < slave_nen_; ++inode)
  {
    double tmp = slave_funct_(inode) * fac;

    for(unsigned idim=0; idim<nsd_; ++idim )
    {
      // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
      iforce[idim+inode*nsd_] += tmp * traction(idim);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::ProjectOnSide(
  LINALG::Matrix<nsd_,1> & x_gp_lin,       ///< global coordinates of gaussian point w.r.t linearized interface
  LINALG::Matrix<nsd_,1> & x_side,         ///< projected gaussian point on side
  LINALG::Matrix<nsd_,1> & xi_side         ///< local coordinates of projected gaussian point w.r.t side
)
{
//  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );

#if(0)
  // REMARK: the current cut_kernel implementation is a factor of 3 slower than this implementation

  GEO::CUT::Position2d<slave_distype> pos( slave_xyze_, x_gp_lin);
  pos.Compute(true);
  pos.LocalCoordinates(xi_side);

  // evaluate shape function at solution
  DRT::UTILS::shape_function_2D( slave_funct_, xi_side( 0 ), xi_side( 1 ), slave_distype );

  // get projected gauss point
  x_side.Multiply(slave_xyze_, slave_funct_);
#endif

  // check, if called on a 3D-element
#ifdef DEBUG
  if (slave_nsd_ == nsd_)
  {
    dserror("You can't project onto a 3D coupling slave element directly. You need an associated boundary element!");
  }
#endif

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );


  if(slave_distype == DRT::Element::tri3 or
     slave_distype == DRT::Element::tri6)
  {
    proj_sol_(0) = 0.333333333333333;
    proj_sol_(1) = 0.333333333333333;
  }
  else if( slave_distype == DRT::Element::quad4 or
           slave_distype == DRT::Element::quad8 or
           slave_distype == DRT::Element::quad9)
  {
    proj_sol_(0) = 0.0;
    proj_sol_(1) = 0.0;
  }
  else
  {
    dserror("Define start side xi-coordinates for unsupported cell type");
  }
  proj_sol_(2) = 0.0;

  // we only use absolute tolerances, since we compute local coordinates
  const double absTolIncr = 1.0e-9;   // abs tolerance for the local coordinates increment
  const double absTolRes  = 1.0e-9;   // abs tolerance for the whole residual
  const double absTOLdist = 1.0e-9;   // abs tolerance for distance

  unsigned iter=0;
  const unsigned maxiter = 7;

  bool converged = false;

  while(iter < maxiter && !converged)
  {
    iter++;

    // get current values
    DRT::UTILS::shape_function_2D( proj_funct_, proj_sol_( 0 ), proj_sol_( 1 ), slave_distype );
    DRT::UTILS::shape_function_2D_deriv1( proj_deriv_, proj_sol_( 0 ), proj_sol_( 1 ), slave_distype );
    DRT::UTILS::shape_function_2D_deriv2( proj_deriv2_, proj_sol_( 0 ), proj_sol_( 1 ), slave_distype );

    proj_x_.Multiply(slave_xyze_, proj_funct_);
    proj_derxy_.MultiplyNT(slave_xyze_, proj_deriv_);
    proj_derxy2_.MultiplyNT(slave_xyze_, proj_deriv2_);

    proj_dx_drdr_times_dx_ds_(0) = proj_derxy2_(1,0)*proj_derxy_(2,1)-proj_derxy_(1,1)*proj_derxy2_(2,0);
    proj_dx_drdr_times_dx_ds_(1) = proj_derxy2_(2,0)*proj_derxy_(0,1)-proj_derxy_(2,1)*proj_derxy2_(0,0);
    proj_dx_drdr_times_dx_ds_(2) = proj_derxy2_(0,0)*proj_derxy_(1,1)-proj_derxy_(0,1)*proj_derxy2_(1,0);

    proj_dx_dr_times_dx_drds_(0) = proj_derxy_(1,0)*proj_derxy2_(2,2)-proj_derxy2_(1,2)*proj_derxy_(2,0);
    proj_dx_dr_times_dx_drds_(1) = proj_derxy_(2,0)*proj_derxy2_(0,2)-proj_derxy2_(2,2)*proj_derxy_(0,0);
    proj_dx_dr_times_dx_drds_(2) = proj_derxy_(0,0)*proj_derxy2_(1,2)-proj_derxy2_(0,2)*proj_derxy_(1,0);

    proj_dx_drds_times_dx_ds_(0) = proj_derxy2_(1,2)*proj_derxy_(2,1)-proj_derxy_(1,1)*proj_derxy2_(2,2);
    proj_dx_drds_times_dx_ds_(1) = proj_derxy2_(2,2)*proj_derxy_(0,1)-proj_derxy_(2,1)*proj_derxy2_(0,2);
    proj_dx_drds_times_dx_ds_(2) = proj_derxy2_(0,2)*proj_derxy_(1,1)-proj_derxy_(0,1)*proj_derxy2_(1,2);

    proj_dx_dr_times_dx_dsds_(0) = proj_derxy_(1,0)*proj_derxy2_(2,1)-proj_derxy2_(1,1)*proj_derxy_(2,0);
    proj_dx_dr_times_dx_dsds_(1) = proj_derxy_(2,0)*proj_derxy2_(0,1)-proj_derxy2_(2,1)*proj_derxy_(0,0);
    proj_dx_dr_times_dx_dsds_(2) = proj_derxy_(0,0)*proj_derxy2_(1,1)-proj_derxy2_(0,1)*proj_derxy_(1,0);

    proj_dx_dr_times_dx_ds_(0) = proj_derxy_(1,0)*proj_derxy_(2,1)-proj_derxy_(1,1)*proj_derxy_(2,0);
    proj_dx_dr_times_dx_ds_(1) = proj_derxy_(2,0)*proj_derxy_(0,1)-proj_derxy_(2,1)*proj_derxy_(0,0);
    proj_dx_dr_times_dx_ds_(2) = proj_derxy_(0,0)*proj_derxy_(1,1)-proj_derxy_(0,1)*proj_derxy_(1,0);

    // define sysmat
    for(unsigned i=0; i< 3; i++)
    {
      // d/dr
      proj_sysmat_(i,0) = proj_derxy_(i,0) - proj_sol_(2) * (proj_dx_drdr_times_dx_ds_(i) + proj_dx_dr_times_dx_drds_(i));

      // d/ds
      proj_sysmat_(i,1) = proj_derxy_(i,1) - proj_sol_(2) * (proj_dx_drds_times_dx_ds_(i) + proj_dx_dr_times_dx_dsds_(i));

      // d/d(dist)
      proj_sysmat_(i,2) = - proj_dx_dr_times_dx_ds_(i);


      // residual
      proj_residuum_(i) = proj_x_(i) - proj_sol_(2) * proj_dx_dr_times_dx_ds_(i) - x_gp_lin(i);

    }

    proj_sysmat_.Invert();

    //solve Newton iteration
    proj_incr_.Multiply(-1.0,proj_sysmat_,proj_residuum_); // incr = -Systemmatrix^-1 * residuum

    // update solution
    proj_sol_.Update(1.0, proj_incr_, 1.0);

    // check ° absolute criterion for local coordinates (between [-1,1]^2)
    //       ° absolute criterion for distance (-> 0)
    //       ° absolute criterion for whole residuum
    if( proj_incr_(2) < absTOLdist
        && sqrt(proj_incr_(0)*proj_incr_(0)+proj_incr_(1)*proj_incr_(1)) <  absTolIncr
        && proj_residuum_.Norm2() < absTolRes)
    {
      converged = true;
    }

  }

  if(!converged)
  {
    std::cout.precision(15);

    std::cout << "increment criterion loc coord "
        //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
        << sqrt(proj_incr_(0)*proj_incr_(0)+proj_incr_(1)*proj_incr_(1))
        << " \tabsTOL: " << absTolIncr
        << std::endl;
    std::cout << "absolute criterion for distance "
        << proj_incr_(2)
        << " \tabsTOL: " << absTOLdist
        << std::endl;
    std::cout << "relative criterion whole residuum "
        << proj_residuum_.Norm2()
        << " \tabsTOL: " << absTolRes
        << std::endl;


    std::cout << "sysmat.Invert" << proj_sysmat_ << std::endl;
    std::cout << "sol-norm " << proj_sol_.Norm2() << std::endl;
    std::cout << "sol " << proj_sol_ << std::endl;
    std::cout << "x_gp_lin" << x_gp_lin << std::endl;
    std::cout << "side " << slave_xyze_ << std::endl;

    dserror( "Newton scheme in ProjectOnSide not converged! " );
  }

  // evaluate shape function at solution
  DRT::UTILS::shape_function_2D( slave_funct_, proj_sol_( 0 ), proj_sol_( 1 ), slave_distype );

  // get projected gauss point
  x_side.Multiply(slave_xyze_, slave_funct_);

  // set local coordinates w.r.t side
  xi_side(0) = proj_sol_(0);
  xi_side(1) = proj_sol_(1);
  xi_side(2) = 0.0; // actually 2D coordinates
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
double SlaveElementRepresentation<distype,slave_distype,slave_numdof>::EvalShapeFuncAndDerivsAtEleCenter()
{
  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<slave_nsd_> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<slave_distype>::rule);

  return EvalShapeFuncAndDerivsAtIntPoint((intpoints_stab.IP().qxg)[0],intpoints_stab.IP().qwgt[0]);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
double SlaveElementRepresentation<distype,slave_distype,slave_numdof>::EvalShapeFuncAndDerivsAtIntPoint(
  const double* gpcoord,  // actual integration point (coords)
  double gpweight// actual integration point (weight)
  )
{
  LINALG::Matrix<nsd_,1>    slave_xsi;     ///< local coordinates of Gaussian point
  LINALG::Matrix<nsd_,nsd_> slave_xjm;     ///< jacobi matrix
  LINALG::Matrix<nsd_,nsd_> slave_xji;     ///< inverse of jacobi matrix

  for (unsigned idim=0;idim<nsd_;idim++)
  {
    slave_xsi(idim) = gpcoord[idim];
  }


  // shape functions and their first derivatives
  DRT::UTILS::shape_function<slave_distype>(slave_xsi,slave_funct_);
  DRT::UTILS::shape_function_deriv1<slave_distype>(slave_xsi,slave_deriv_);

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  slave_xjm.MultiplyNT(slave_deriv_,slave_xyze_);
  double det = slave_xji.Invert(slave_xjm);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", det);

  // compute global first derivates
  slave_derxy_.Multiply(slave_xji,slave_deriv_);

  // compute integration factor
  return gpweight*det;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::GetSlaveFunctDeriv(LINALG::Matrix<nsd_,slave_nen_>& slave_derxy) const
{
  slave_derxy = slave_derxy_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::ApplyConvStabTerms(
  const Teuchos::RCP<SlaveElementInterface<distype> > & slave_ele,              ///< associated slave element coupling object
  const LINALG::Matrix<nen_,1> &                        funct_m,                ///< master shape functions
  const LINALG::Matrix<nsd_,1> &                        velint_m,               ///< vector of slave shape functions
  const LINALG::Matrix<nsd_,1> &                        normal,                 ///< normal vector n^b
  const double &                                        density_m,              ///< fluid density (master)
  const double &                                        NIT_stab_fac_conv,      ///< full Nitsche's penalty term scaling (viscous+convective part)
  const double &                                        timefacfac,             ///< theta*dt
  const LINALG::Matrix<nsd_,1> &                        ivelint_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
  const INPAR::XFEM::EleCouplingCondType &              cond_type,              ///< condition type
  INPAR::XFEM::XFF_ConvStabScaling                      xff_conv_stab           ///< type of convective stabilization in XFF-problems
)
{
  if (cond_type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID &&
      xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_none)
    dserror("Cannot apply convective stabilization terms for XFF_ConvStabScaling_none!");

  // funct_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_timefacfac(funct_m);
  funct_m_timefacfac.Scale(timefacfac);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  LINALG::Matrix<nen_,nen_> funct_m_m_dyad_timefacfac(true);
  funct_m_m_dyad_timefacfac.MultiplyNT(funct_m_timefacfac, funct_m);

  // velint_s
  velint_s_.Clear();

  if(eval_coupling_)
    slave_ele->GetInterfaceVelnp(velint_s_);

  // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump height for coupled problems
  velint_s_.Update(1.0, ivelint_jump, 1.0);

  velint_diff_.Update(1.0, velint_m, -1.0, velint_s_, 0.0);

  velint_diff_normal_ = velint_diff_.Dot(normal);


  // REMARK:
  // the (additional) convective stabilization is included in NIT_full_stab_fac
  // (in case of mixed/hybrid LM approaches, we don't compute the penalty term explicitly -
  // it 'evolves'); in that case we therefore don't chose the maximum, but add the penalty term scaled
  // with conv_stab_fac to the viscous counterpart;
  // this happens by calling NIT_Stab_Penalty


  switch (cond_type)
  {
    case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
    case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
    case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
    case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
    {
      velint_diff_stabfac_.Update(NIT_stab_fac_conv, velint_diff_, 0.0);

      NIT_Stab_Penalty_MasterTerms(
        funct_m_timefacfac,
        funct_m_m_dyad_timefacfac,
        NIT_stab_fac_conv
      );
      break;
    }
    case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
    {
      // funct_s
      Teuchos::RCP<SlaveElementRepresentation<distype,slave_distype,slave_numdof> > ser =
        Teuchos::rcp_dynamic_cast<SlaveElementRepresentation<distype,slave_distype,slave_numdof> >(slave_ele);
      if (ser == Teuchos::null)
        dserror("Failed to cast slave_ele to SlaveElementRepresentation!");
      LINALG::Matrix<slave_nen_,1> funct_s;
      ser->GetSlaveFunct(funct_s);

      // funct_s * timefac * fac
      LINALG::Matrix<slave_nen_,1> funct_s_timefacfac(funct_s);
      funct_s_timefacfac.Scale(timefacfac);

      // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
      LINALG::Matrix<slave_nen_,slave_nen_> funct_s_s_dyad_timefacfac(true);
      funct_s_s_dyad_timefacfac.MultiplyNT(funct_s_timefacfac, funct_s);

      LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac(true);
      funct_s_m_dyad_timefacfac.MultiplyNT(funct_s_timefacfac, funct_m);

      if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_upwinding)
      {
        NIT_Stab_Penalty(
          funct_m_timefacfac,
          funct_m_m_dyad_timefacfac,
          funct_s_timefacfac,
          funct_s_m_dyad_timefacfac,
          funct_s_s_dyad_timefacfac,
          NIT_stab_fac_conv
        );
      }

      // prevent instabilities due to convective mass transport across the fluid-fluid interface
      if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_only_averaged)
      {
        NIT_Stab_Inflow_AveragedTerm(
          velint_m,
          velint_s_,
          funct_m_timefacfac,
          funct_m_m_dyad_timefacfac,
          funct_s_timefacfac,
          funct_s_m_dyad_timefacfac,
          funct_s_s_dyad_timefacfac,
          normal,
          density_m
        );
      }
      break;
    }
    case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
    {
      dserror("Convective stabilization in monolithic XFSI is not yet available!");
      break;
    }
    default:
      dserror("Unsupported coupling condition type. Cannot apply convective stabilization terms."); break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_evaluateCoupling(
  const LINALG::Matrix<nsd_,1> &    normal,                 ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
  const double &                    timefacfac,             ///< theta*dt
  const double &                    visceff_m,              ///< viscosity in coupling master fluid
  const double &                    visceff_s,              ///< viscosity in coupling slave fluid
  const double &                    kappa_m,                ///< mortaring weight for coupling master
  const double &                    kappa_s,                ///< mortaring weight for coupling slave
  const double &                    density_m,              ///< fluid density (master) USED IN XFF
  const double &                    NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
  const LINALG::Matrix<nen_,1> &    funct_m,                ///< coupling master shape functions
  const LINALG::Matrix<nsd_,nen_> & derxy_m,                ///< spatial derivatives of coupling master shape functions
  const LINALG::Matrix<nsd_,nsd_> & vderxy_m,               ///< coupling master spatial velocity derivatives
  const double &                    pres_m,                 ///< coupling master pressure
  const LINALG::Matrix<nsd_,1> &    velint_m,               ///< coupling master interface velocity
  const LINALG::Matrix<nsd_,1> &    ivelint_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
  const LINALG::Matrix<nsd_,1> &    itraction_jump,         ///< prescribed interface traction, jump height for coupled problems
  INPAR::XFEM::XFF_ConvStabScaling  xff_conv_stab           ///< type of convective stabilization in XFF-problems
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_evaluateCoupling");

  //--------------------------------------------

  // define the coupling between two not matching grids
  // for fluidfluidcoupling
  // domain Omega^m := Coupling master (XFluid)
  // domain Omega^s := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

  // [| v |] := vm - vs
  //  { v }  := kappa_m * vm + kappas * vs = kappam * vm (for Dirichlet coupling km=1.0, ks = 0.0)
  //  < v >  := kappa_s * vm + kappam * vs = kappas * vm (for Dirichlet coupling km=1.0, ks = 0.0)

  //--------------------------------------------

  //--------------------------------------------

  // plausibility check
  if (kappa_s * kappa_m < 0)
    dserror("Please choose reasonable, non-negative weights");

  // check if this is a one-sided background fluid-master approach
  const bool full_mastersided(fabs(kappa_m-1.0) < 1.e-14);
  // check if this is a slave-sided (embedded-sided approach=
  const bool full_slavesided(fabs(kappa_s-1.0) < 1.e-14);

  // check the weights
  if (fabs(kappa_m + kappa_s -1.0) > 1.e-14)
    dserror("Invalid weights kappa_m=%d and kappa_s=%d. The have to sum up to 1.0!", kappa_m, kappa_s);

  //--------------------------------------------

  half_normal_.Update(0.5,normal,0.0);

  half_normal_deriv_m_.MultiplyTN(derxy_m, half_normal_); // half_normal(k)*derxy_m(k,ic);

  vderxy_m_normal_.Multiply(vderxy_m, half_normal_);

  vderxy_m_normal_transposed_.MultiplyTN(vderxy_m, half_normal_);
  vderxy_m_normal_transposed_.Update(1.0,vderxy_m_normal_,1.0);


  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  // interface velocity vector in gausspoint
  velint_s_.Clear();

  if(eval_coupling_)
    this->GetInterfaceVelnp(velint_s_);

  // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump height for coupled problems
  velint_s_.Update(1.0, ivelint_jump, 1.0);

  velint_diff_.Update(1.0, velint_m, -1.0, velint_s_, 0.0);

  velint_diff_stabfac_.Update(NIT_full_stab_fac, velint_diff_, 0.0);

  velint_diff_normal_ = velint_diff_.Dot(normal);

  // funct_m * timefac * fac
  funct_m_timefacfac_.Update(timefacfac,funct_m, 0.0 );

  // funct_s * timefac * fac
  funct_s_.Clear();
  if (slave_distype != DRT::Element::dis_none)
    this->GetSlaveFunct(funct_s_);

  funct_s_timefacfac_.Update(timefacfac, funct_s_, 0.0);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  funct_m_m_dyad_timefacfac_.MultiplyNT(funct_m_timefacfac_, funct_m);

  // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
  funct_s_s_dyad_timefacfac_.MultiplyNT(funct_s_timefacfac_, funct_s_);

  // funct_s * timefac * fac * funct_m (dyadic product)
  funct_s_m_dyad_timefacfac_.MultiplyNT(funct_s_timefacfac_, funct_m);

  // 2 * mu_m * timefacfac
  const double viscm_fac = 2.0 * timefacfac * visceff_m;

  // compute normal velocity components
  const double velint_normal_m = velint_m.Dot(normal);
  const double velint_normal_s = velint_s_.Dot(normal);


  //-----------------------------------------------------------------
  // viscous stability term
  // REMARK: this term includes also inflow coercivity in case of XFSI
  // with modified stabfac (see NIT_ComputeStabfac)
  if ( !eval_coupling_ )
  {
    NIT_Stab_Penalty_MasterTerms(
      funct_m_timefacfac_,
      funct_m_m_dyad_timefacfac_,
      NIT_full_stab_fac
    );
  }
  else
  {
    NIT_Stab_Penalty(
      funct_m_timefacfac_,
      funct_m_m_dyad_timefacfac_,
      funct_s_timefacfac_,
      funct_s_m_dyad_timefacfac_,
      funct_s_s_dyad_timefacfac_,
      NIT_full_stab_fac
    );
  }

  // add averaged term
  if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
      xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_only_averaged )
  {
    NIT_Stab_Inflow_AveragedTerm(
      velint_m,
      velint_s_,
      funct_m_timefacfac_,
      funct_m_m_dyad_timefacfac_,
      funct_s_timefacfac_,
      funct_s_m_dyad_timefacfac_,
      funct_s_s_dyad_timefacfac_,
      normal,
      density_m
    );
  }

  // funct_s * timefac * fac * kappa_m
  funct_s_timefacfac_km_.Update(kappa_m, funct_s_timefacfac_, 0.0);

  // funct_m * timefac * fac * kappa_s
  funct_m_timefacfac_ks_.Update(kappa_s, funct_m_timefacfac_, 0.0);

   //-----------------------------------------------------------------
   // standard consistency traction jump term
   if( eval_coupling_ )
   {
     NIT_Traction_Consistency_Term(
         funct_m_timefacfac_ks_,
         funct_s_timefacfac_km_,
         itraction_jump
         );
   }

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - standard Dirichlet case/pure xfluid-sided case

  if (!eval_coupling_ || full_mastersided)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterTerms(
        pres_m,
        funct_m_timefacfac_,
        funct_s_timefacfac_,
        normal,
        funct_s_m_dyad_timefacfac_,
        funct_m_m_dyad_timefacfac_
        );

    //-----------------------------------------------------------------
    // pressure adjoint consistency term

    NIT_p_AdjointConsistency_MasterTerms(
        velint_normal_m,
        velint_normal_s,
        funct_m_timefacfac_,
        funct_s_m_dyad_timefacfac_,
        funct_m_m_dyad_timefacfac_,
        normal
        );

    //-----------------------------------------------------------------
    // viscous consistency term

    // funct_m * 2 * mu_m * kappa_m * timefac * fac
    funct_m_viscm_timefacfac_.Update(viscm_fac,funct_m,0.0);

    // funct_s * 2 * mu_m * kappa_m * timefac * fac
    funct_s_viscm_timefacfac_.Update(viscm_fac,funct_s_,0.0);

    NIT_visc_Consistency_MasterTerms(
        derxy_m,
        vderxy_m,
        funct_m_viscm_timefacfac_,
        funct_s_viscm_timefacfac_,
        normal);

    //-----------------------------------------------------------------
    // viscous adjoint consistency term

    double tmp_fac = adj_visc_scale_*viscm_fac;
    derxy_m_viscm_timefacfac_.Update(tmp_fac,derxy_m);
    NIT_visc_AdjointConsistency_MasterTerms(
        velint_m,
        velint_s_,
        derxy_m_viscm_timefacfac_,
        funct_m,
        funct_s_,
        normal );

    // we are done here!
    return;
  }

  //TODO: introduce member for the subsequent matrices or introduce the scalings otherwise

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - two-sided or xfluid-sided:

  // funct_m * timefac * fac * kappa_m
  LINALG::Matrix<nen_,1> funct_m_timefacfac_km(funct_m_timefacfac_);
  funct_m_timefacfac_km.Scale(kappa_m);



  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  LINALG::Matrix<nen_,nen_> funct_m_m_dyad_timefacfac_km(funct_m_m_dyad_timefacfac_);
  funct_m_m_dyad_timefacfac_km.Scale(kappa_m);

  // funct_s * timefac * fac * funct_m * kappa_m
  LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac_km(funct_s_m_dyad_timefacfac_);
  funct_s_m_dyad_timefacfac_km.Scale(kappa_m);

  // 2 * mu_m * kappa_m * timefac * fac
  const double km_viscm_fac = viscm_fac * kappa_m;

  // funct_m * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_viscm_timefacfac_km(funct_m);
  funct_m_viscm_timefacfac_km.Scale(km_viscm_fac);

  // funct_s * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<slave_nen_,1> funct_s_viscm_timefacfac_km(funct_s_);
  funct_s_viscm_timefacfac_km.Scale(km_viscm_fac);


  if (! full_slavesided)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterTerms(
        pres_m,
        funct_m_timefacfac_km,
        funct_s_timefacfac_km_,
        normal,
        funct_s_m_dyad_timefacfac_km,
        funct_m_m_dyad_timefacfac_km
        );

    //-----------------------------------------------------------------
    // pressure adjoint consistency term

    NIT_p_AdjointConsistency_MasterTerms(
        velint_normal_m,
        velint_normal_s,
        funct_m_timefacfac_km,
        funct_s_m_dyad_timefacfac_km,
        funct_m_m_dyad_timefacfac_km,
        normal
        );

    //-----------------------------------------------------------------
    // viscous consistency term
    NIT_visc_Consistency_MasterTerms(
        derxy_m,
        vderxy_m,
        funct_m_viscm_timefacfac_km,
        funct_s_viscm_timefacfac_km,
        normal);

    //-----------------------------------------------------------------
    // viscous adjoint consistency term

    LINALG::Matrix<nsd_,nen_> derxy_m_viscm_timefacfac_km(derxy_m);
    derxy_m_viscm_timefacfac_km.Scale(adj_visc_scale_*km_viscm_fac);
    NIT_visc_AdjointConsistency_MasterTerms(
        velint_m,
        velint_s_,
        derxy_m_viscm_timefacfac_km,
        funct_m,
        funct_s_,
        normal
        );

  }

  //-----------------------------------------------------------------
  // the following quantities are only required for two-sided coupling
  // kappa_s > 0.0

  // 2 * mu_s * kappa_s * timefac * fac
  const double ks_viscs_fac = 2.0 * kappa_s * visceff_s * timefacfac;

  // funct_s * timefac * fac * funct_m * kappa_s
  LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac_ks(funct_s_m_dyad_timefacfac_);
  funct_s_m_dyad_timefacfac_ks.Scale(kappa_s);

  // funct_s * timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_ks(funct_s_timefacfac_);
  funct_s_timefacfac_ks.Scale(kappa_s);

  // funct_s * timefac * fac * funct_s * kappa_s
  LINALG::Matrix<slave_nen_,slave_nen_> funct_s_s_dyad_timefacfac_ks(funct_s_s_dyad_timefacfac_);
  funct_s_s_dyad_timefacfac_ks.Scale(kappa_s);

  // funct_m * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<nen_,1> funct_m_viscs_timefacfac_ks(funct_m);
  funct_m_viscs_timefacfac_ks.Scale(ks_viscs_fac);

  // funct_s * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_viscs_timefacfac_ks(funct_s_);
  funct_s_viscs_timefacfac_ks.Scale(ks_viscs_fac);

  //-----------------------------------------------------------------
  // pressure consistency term

  double pres_s = 0.0;
  // must use this-pointer because of two-stage lookup!
  this->GetInterfacePresnp(pres_s);

  NIT_p_Consistency_SlaveTerms(
      pres_s,
      funct_m_timefacfac_ks_,
      funct_s_timefacfac_ks,
      funct_s_m_dyad_timefacfac_ks,
      funct_s_s_dyad_timefacfac_ks,
      normal);

  //-----------------------------------------------------------------
  // pressure adjoint consistency term

  NIT_p_AdjointConsistency_SlaveTerms(
      velint_normal_m,
      velint_normal_s,
      funct_s_timefacfac_ks,
      funct_s_m_dyad_timefacfac_ks,
      funct_s_s_dyad_timefacfac_ks,
      normal);

  //-----------------------------------------------------------------
  // viscous consistency term

  // Shape function derivatives for slave side
  LINALG::Matrix<nsd_,slave_nen_> derxy_s;
  this->GetSlaveFunctDeriv(derxy_s);

  // Spatial velocity gradient for slave side
  LINALG::Matrix<nsd_,nsd_> vderxy_s;
  this->GetInterfaceVelGradnp(vderxy_s);

  half_normal_deriv_s_.MultiplyTN(derxy_s, half_normal_); // half_normal(k)*derxy_s(k,ic);
  vderxy_s_normal_.Multiply(vderxy_s, half_normal_);
  vderxy_s_normal_transposed_.MultiplyTN(vderxy_s, half_normal_);
  vderxy_s_normal_transposed_.Update(1.0,vderxy_s_normal_,1.0);

  NIT_visc_Consistency_SlaveTerms(
      derxy_s,
      funct_m_viscs_timefacfac_ks,
      funct_s_viscs_timefacfac_ks,
      normal);

  //-----------------------------------------------------------------
  // viscous adjoint consistency term

  LINALG::Matrix<nsd_,slave_nen_> derxy_s_viscs_timefacfac_ks(derxy_s);
  derxy_s_viscs_timefacfac_ks.Scale(adj_visc_scale_*ks_viscs_fac);

  NIT_visc_AdjointConsistency_SlaveTerms(
      velint_m,
      velint_s_,
      funct_m,
      funct_s_,
      derxy_s_viscs_timefacfac_ks,
      normal
  );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_evaluateCouplingOldState(
  const LINALG::Matrix<nsd_,1> &           normal,                  ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
  const double &                           timefacfacn,             ///< dt*(1-theta)*fac
  const double &                           visceff_m,               ///< viscosity in coupling master fluid
  const double &                           visceff_s,               ///< viscosity in coupling slave fluid
  const double &                           kappa_m,                 ///< mortaring weight for coupling master
  const double &                           kappa_s,                 ///< mortaring weight for coupling slave
  const double &                           density_m,              ///< fluid density (master) USED IN XFF
  const LINALG::Matrix<nen_,1> &           funct_m,                 ///< coupling master shape functions
  const LINALG::Matrix<nsd_,nen_> &        derxy_m,                 ///< spatial derivatives of coupling master shape functions
  const LINALG::Matrix<nsd_,nsd_> &        vderxyn_m,               ///< coupling master spatial velocity derivatives
  const double &                           presn_m,                 ///< coupling master pressure
  const LINALG::Matrix<nsd_,1> &           velintn_m,               ///< coupling master interface velocity
  const LINALG::Matrix<nsd_,1> &           ivelintn_jump,           ///< prescribed interface velocity, Dirichlet values or jump height for coupled problems
  const LINALG::Matrix<nsd_,1> &           itractionn_jump,         ///< prescribed interface traction, jump height for coupled problems
  INPAR::XFEM::InterfaceTermsPreviousState prev_state,              ///< evaluate all terms or only consistency term at previous state
  const double                             NIT_full_stab_facn,      ///< full Nitsche's penalty term scaling (viscous+convective part)
  INPAR::XFEM::XFF_ConvStabScaling         xff_conv_stab            ///< type of convective stabilization in XFF-problems
)
{
  //--------------------------------------------

  // define the coupling between two not matching grids
  // for fluidfluidcoupling
  // domain Omega^m := Coupling master (XFluid)
  // domain Omega^s := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

  // [| v |] := vm - vs
  //  { v }  := kappa_m * vm + kappas * vs = kappam * vm (for Dirichlet coupling km=1.0, ks = 0.0)
  //  < v >  := kappa_s * vm + kappam * vs = kappas * vm (for Dirichlet coupling km=1.0, ks = 0.0)

  //--------------------------------------------

  //--------------------------------------------

  // plausibility check
  if (kappa_s * kappa_m < 0)
    dserror("Please choose reasonable, non-negative weights");

  // check if this is a one-sided background fluid-master approach
  const bool full_mastersided(fabs(kappa_m-1.0) < 1.e-14);
  // check if this is a slave-sided (embedded-sided approach=
  const bool full_slavesided(fabs(kappa_s-1.0) < 1.e-14);

  // check the weights
  if (fabs(kappa_m + kappa_s -1.0) > 1.e-14)
    dserror("Invalid weights kappa_m=%d and kappa_s=%d. The have to sum up to 1.0!", kappa_m, kappa_s);

  half_normal_.Update(0.5,normal,0.0);
  half_normal_deriv_m_.MultiplyTN(derxy_m, half_normal_); // half_normal(k)*derxy_m(k,ic);

  vderxy_m_normal_.Multiply(vderxyn_m, half_normal_);

  vderxy_m_normal_transposed_.MultiplyTN(vderxyn_m, half_normal_);
  vderxy_m_normal_transposed_.Update(1.0,vderxy_m_normal_,1.0);

  //--------------------------------------------

  // get velocity at integration point
  // (values at n)
  // interface velocity vector in gausspoint
  velint_s_.Clear();
  if (eval_coupling_)
    this->GetInterfaceVeln(velint_s_);

  // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump height for coupled problems
  velint_s_.Update(1.0, ivelintn_jump, 1.0);

  velint_diff_.Update(1.0, velintn_m, -1.0, velint_s_, 0.0);

  velint_diff_normal_ = velint_diff_.Dot(normal);

  // funct_m * timefac * fac
  // funct_m * timefac * fac
  funct_m_timefacfac_.Update(timefacfacn,funct_m, 0.0 );

  // funct_s * timefac * fac
  funct_s_.Clear();
  if (slave_distype != DRT::Element::dis_none)
    this->GetSlaveFunct(funct_s_);
  funct_s_timefacfac_.Update(timefacfacn, funct_s_, 0.0);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  funct_m_m_dyad_timefacfac_.MultiplyNT(funct_m_timefacfac_, funct_m);

  // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
  funct_s_s_dyad_timefacfac_.MultiplyNT(funct_s_timefacfac_, funct_s_);

  // funct_s * timefac * fac * funct_m (dyadic product)
  funct_s_m_dyad_timefacfac_.MultiplyNT(funct_s_timefacfac_, funct_m);

  // 2 * mu_m * timefacfacn
  const double viscm_facn = 2.0 * timefacfacn * visceff_m;

  // compute normal velocity components
  const double velintn_normal_m = velintn_m.Dot(normal);
  const double velintn_normal_s = velint_s_.Dot(normal);

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - standard Dirichlet case/pure xfluid-sided case

  if (!eval_coupling_ || full_mastersided)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterRHS(
        presn_m,
        funct_m_timefacfac_,
        funct_s_timefacfac_,
        normal
    );

    //-----------------------------------------------------------------
    // viscous consistency term

    // funct_m * 2 * mu_m * kappa_m * timefac * fac
    funct_m_viscm_timefacfac_.Update(viscm_facn,funct_m,0.0);

    // funct_s * 2 * mu_m * kappa_m * timefac * fac
    funct_s_viscm_timefacfac_.Update(viscm_facn,funct_s_,0.0);

    NIT_visc_Consistency_MasterRHS(
        funct_m_viscm_timefacfac_,
        funct_s_viscm_timefacfac_,
        normal);

    // done, if only consistency terms are evaluated for previos time step
    if (prev_state != INPAR::XFEM::PreviousState_full)
      return;

    //-----------------------------------------------------------------
    // pressure adjoint consistency term

    NIT_p_AdjointConsistency_MasterRHS(
        velintn_normal_m,
        velintn_normal_s,
        funct_m_timefacfac_,
        normal
    );
    //-----------------------------------------------------------------
    // viscous adjoint consistency term

    LINALG::Matrix<nsd_,nen_> derxy_m_viscm_timefacfacn(derxy_m);
    derxy_m_viscm_timefacfacn.Scale(adj_visc_scale_*viscm_facn);
    NIT_visc_AdjointConsistency_MasterRHS(
        velintn_m,
        velint_s_,
        derxy_m_viscm_timefacfacn,
        normal );

    //-----------------------------------------------------------------
    // penalty term

    NIT_Stab_Penalty_MasterRHS(
      velintn_m,
      velint_s_,
      funct_m_timefacfac_,
      NIT_full_stab_facn
    );

    // we are done here!
    return;
  }

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - two-sided or xfluid-sided:

  // funct_m * timefac * fac * kappa_m
  LINALG::Matrix<nen_,1> funct_m_timefacfacn_km(funct_m_timefacfac_);
  funct_m_timefacfacn_km.Scale(kappa_m);

  // funct_s * timefac * fac * kappa_m
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfacn_km(funct_s_timefacfac_);
  funct_s_timefacfacn_km.Scale(kappa_m);

  // 2 * mu_m * kappa_m * timefac * fac
  const double km_viscm_facn = viscm_facn * kappa_m;

  // funct_m * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_viscm_timefacfacn_km(funct_m);
  funct_m_viscm_timefacfacn_km.Scale(km_viscm_facn);

  // funct_s * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<slave_nen_,1> funct_s_viscm_timefacfacn_km(funct_s_);
  funct_s_viscm_timefacfacn_km.Scale(km_viscm_facn);

  if (! full_slavesided)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterRHS(
        presn_m,
        funct_m_timefacfacn_km,
        funct_s_timefacfacn_km,
        normal
    );

    //-----------------------------------------------------------------
    // viscous consistency term
    NIT_visc_Consistency_MasterRHS(
        funct_m_viscm_timefacfacn_km,
        funct_s_viscm_timefacfacn_km,
        normal);

    if (prev_state == INPAR::XFEM::PreviousState_full)
    {
      //-----------------------------------------------------------------
      // pressure adjoint consistency term
      NIT_p_AdjointConsistency_MasterRHS(
          velintn_normal_m,
          velintn_normal_s,
          funct_m_timefacfacn_km,
          normal
      );
      //-----------------------------------------------------------------
      // viscous adjoint consistency term

      LINALG::Matrix<nsd_,nen_> derxy_m_viscm_timefacfacn_km(derxy_m);
      derxy_m_viscm_timefacfacn_km.Scale(adj_visc_scale_*km_viscm_facn);
      NIT_visc_AdjointConsistency_MasterRHS(
          velintn_m,
          velint_s_,
          derxy_m_viscm_timefacfacn_km,
          normal
      );
    }
  }

  //-----------------------------------------------------------------
  // the following quantities are only required for two-sided coupling
  // kappa_s > 0.0

  // 2 * mu_s * kappa_s * timefac * fac
  const double ks_viscs_facn = 2.0 * kappa_s * visceff_s * timefacfacn;

  // funct_m * timefac * fac * kappa_s
  LINALG::Matrix<nen_,1> funct_m_timefacfacn_ks(funct_m_timefacfac_);
  funct_m_timefacfacn_ks.Scale(kappa_s);

  // funct_s * timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfacn_ks(funct_s_timefacfac_);
  funct_s_timefacfacn_ks.Scale(kappa_s);


  // funct_m * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<nen_,1> funct_m_viscs_timefacfacn_ks(funct_m);
  funct_m_viscs_timefacfacn_ks.Scale(ks_viscs_facn);

  // funct_s * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_viscs_timefacfacn_ks(funct_s_);
  funct_s_viscs_timefacfacn_ks.Scale(ks_viscs_facn);

  //-----------------------------------------------------------------
  // standard consistency traction jump term
  if( eval_coupling_ )
  {
    NIT_Traction_Consistency_Term(
        funct_m_timefacfacn_ks,
        funct_s_timefacfacn_km,
        itractionn_jump
        );
  }

 //-----------------------------------------------------------------


  //-----------------------------------------------------------------
  // pressure consistency term

  double presn_s = 0.0;
  // must use this-pointer because of two-stage lookup!
  this->GetInterfacePresn(presn_s);

  NIT_p_Consistency_SlaveRHS(
      presn_s,
      funct_m_timefacfacn_ks,
      funct_s_timefacfacn_ks,
      normal);

  //-----------------------------------------------------------------
  // viscous consistency term

  // Spatial velocity gradient for slave side
  LINALG::Matrix<nsd_,nsd_> vderxyn_s;
  this->GetInterfaceVelGradn(vderxyn_s);

  vderxy_s_normal_.Multiply(vderxyn_s, half_normal_);
  vderxy_s_normal_transposed_.MultiplyTN(vderxyn_s, half_normal_);
  vderxy_s_normal_transposed_.Update(1.0,vderxy_s_normal_,1.0);

  NIT_visc_Consistency_SlaveRHS(
      vderxyn_s,
      funct_m_viscs_timefacfacn_ks,
      funct_s_viscs_timefacfacn_ks,
      normal);

  if (prev_state != INPAR::XFEM::PreviousState_full)
    return;

  //-----------------------------------------------------------------
  // pressure adjoint consistency term

  NIT_p_AdjointConsistency_SlaveRHS(
      velintn_normal_m,
      velintn_normal_s,
      funct_s_timefacfacn_ks);
  //-----------------------------------------------------------------
  // viscous adjoint consistency term
  // Shape function derivatives for slave side
  LINALG::Matrix<nsd_,slave_nen_> derxyn_s;
  this->GetSlaveFunctDeriv(derxyn_s);

  LINALG::Matrix<nsd_,slave_nen_> derxy_s_viscs_timefacfacn_ks(derxyn_s);
  derxy_s_viscs_timefacfacn_ks.Scale(adj_visc_scale_*ks_viscs_facn);

  NIT_visc_AdjointConsistency_SlaveRHS(
      velintn_m,
      velint_s_,
      derxy_s_viscs_timefacfacn_ks,
      normal
  );

  //-----------------------------------------------------------------
  // penalty term

  NIT_Stab_Penalty_SlaveRHS(
    velintn_m,
    velint_s_,
    funct_s_timefacfac_,
    NIT_full_stab_facn
  );

  // add averaged term
  if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_upwinding ||
      xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_only_averaged)
  {
    NIT_Stab_Inflow_AveragedTermRHS(
      velintn_m,
      velint_s_,
      funct_m_timefacfac_,
      funct_s_timefacfac_,
      normal,
      density_m
    );
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Traction_Consistency_Term(
  const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_s
  const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &            itraction_jump                ///< prescribed interface traction, jump height for coupled problems
)
{
        /*            \
     - |  < v >,   t   |   with t = [sigma * n]
        \             /     */

  // Two-Phase Flow:
  //
  //     t_{n+1}          [| sigma*n |] =      gamma * curv * n
  //                                                             with curv = div(grad(phi)/||grad(phi)||)
  //
  //      t_{n}           [| sigma*n |] = [|  -pI + \mu*[\nabla u + (\nabla u)^T]  |] * n
  //

  // Combustion:        TO BE IMPLEMENTED.
  //

  // All else:            [| sigma*n |] = 0


  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    //-----------------------------------------------
    //    - (vm, ks * t)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const double funct_m_ks_timefacfac_traction = funct_m_timefacfac_ks(ir)*itraction_jump(ivel);

      const unsigned row = mIndex(ir,ivel);
      rhC_um_(row,0) += funct_m_ks_timefacfac_traction;
    }

    if (!eval_coupling_ ) continue;

    //-----------------------------------------------
    //    + (vs, km * t)
    //-----------------------------------------------
    for(unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const double funct_s_km_timefacfac_traction = funct_s_timefacfac_km(ir)*itraction_jump(ivel);

      const unsigned row = sIndex(ir,ivel);
      rhC_us_(row,0) += funct_s_km_timefacfac_traction;
    }
  } // end loop over velocity components
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_MasterTerms(
  const double &                            pres_m,                       ///< master pressure
  const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &            normal,                       ///< normal vector
  const LINALG::Matrix<slave_nen_,nen_> &   funct_s_m_dyad_timefacfac_km, ///< (funct_s^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nen_,nen_> &         funct_m_m_dyad_timefacfac_km  ///< (funct^T * funct) * timefacfac *kappa_m
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_p_Consistency_MasterTerms");


  /*                  \       /          i      \
+ |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
  \                   /       \                */

  //-----------------------------------------------
  //    + (vm, km *(Dpm)*n)
  //-----------------------------------------------
  for (unsigned ic =0; ic<nen_; ic++)
  {
    const unsigned col = mPres(ic);

    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const double tmp = funct_m_m_dyad_timefacfac_km(ir,ic);

      // loop over velocity components
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // (v,Dp*n)
        C_umum_(mIndex(ir,ivel), col) += tmp*normal(ivel);
      }
    }
  }

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    const double funct_m_km_timefacfac_press = funct_m_timefacfac_km(ir)*pres_m;

    // loop over velocity components
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // -(v,p*n)
      rhC_um_( mIndex(ir,ivel),0) -= funct_m_km_timefacfac_press*normal(ivel);
    }
  }


  if (!eval_coupling_) return;


  for(unsigned ic =0; ic<nen_; ic++)
  {
    const unsigned col = mPres(ic);

    for(unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const double tmp = funct_s_m_dyad_timefacfac_km(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        //-----------------------------------------------
        //    - (vs, km *(Dpm)*n)
        //-----------------------------------------------

        // (v,Dp*n)
        C_usum_(sIndex(ir,ivel), col) -= tmp*normal(ivel);
      }
    }
  }

  for(unsigned ir = 0; ir<slave_nen_; ir++)
  {
    const double funct_s_km_timefacfac_pres = funct_s_timefacfac_km(ir)*pres_m;

    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // -(v,p*n)
      rhC_us_(sIndex(ir,ivel),0) += funct_s_km_timefacfac_pres*normal(ivel);
    }
  } // end loop over velocity components

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_SlaveTerms(
  const double &                                pres_s,                       ///< slave pressure
  const LINALG::Matrix<nen_,1> &                funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,nen_> &       funct_s_m_dyad_timefacfac_ks, ///< (funct_s^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac_ks, ///< (funct^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
)
{
  for (unsigned ic =0; ic<slave_nen_; ic++)
  {
    unsigned col = sPres(ic);

    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const double tmp = funct_s_m_dyad_timefacfac_ks(ic,ir);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        //-----------------------------------------------
        //    + (vm, ks *(Dps)*n)
        //-----------------------------------------------

        // (vm, ks * Dps*n)
        C_umus_(mIndex(ir,ivel),col) += tmp*normal(ivel);
      }
    }

    //-----------------------------------------------
    //    - (vs, ks *(Dps)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const double tmp = funct_s_s_dyad_timefacfac_ks(ir,ic);
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        C_usus_(sIndex(ir,ivel),col) -= tmp*normal(ivel);
      }
    }
  }

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    const double tmp = pres_s*funct_m_timefacfac_ks(ir);
    // loop over velocity components
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // -(vm, ks * ps*n)
      rhC_um_(mIndex(ir,ivel),0) -= tmp*normal(ivel);
    }
  }

  //-----------------------------------------------
  //    - (vs, ks *(Dps)*n)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    const double tmp = funct_s_timefacfac_ks(ir)*pres_s;
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // +(vs,ks * ps*n)
      rhC_us_(sIndex(ir,ivel),0) += tmp*normal(ivel);
    }
  } // end loop over velocity components
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_AdjointConsistency_MasterTerms(
  const double &                          velint_normal_m,              ///< velocity in normal direction
  const double &                          velint_normal_s,              ///< interface velocity in normal direction
  const LINALG::Matrix<nen_,1> &          funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,nen_> & funct_s_m_dyad_timefacfac_km, ///< (funct^T * funct_s) * timefacfac *kappa_m
  const LINALG::Matrix<nen_,nen_> &       funct_m_m_dyad_timefacfac_km, ///< (funct^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_p_AdjointConsistency_MasterTerms");


  /*                   \     /              i   \
- |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
  \                    /     \                 */

  // REMARK:
  // the sign of the pressure adjoint consistency term is opposite to the sign of the pressure consistency term (interface),
  // as a non-symmetric pressure formulation is chosen in the standard fluid
  // the sign of the standard volumetric pressure consistency term is opposite to the (chosen) sign of the pressure-weighted continuity residual;
  // think about the Schur-complement for the Stokes problem: S_pp = A_pp - A_pu A_uu^-1 A_up
  // (--> A_pu == -A_up^T; sgn(A_pp) == sgn(- A_pu A_uu^-1 Aup), where A_pp comes from pressure-stabilizing terms)
  // a symmetric adjoint pressure consistency term would also affect the sign of the pressure stabilizing terms
  // for Stokes' problem, this sign choice leads to a symmetric, positive definite Schur-complement matrix S
  // (v, p*n)--> A_up; -(q,u*n)--> -A_up^T; S_pp = A_pp + A_up^T A_uu A_up

  //-----------------------------------------------
  //    - (qm*n, km *(Dum))
  //-----------------------------------------------
  for (unsigned ic =0; ic<nen_; ic++)
  {
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const unsigned row = mPres(ir);

      const double tmp = funct_m_m_dyad_timefacfac_km(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // - (qm*n, km *(Dum))
        C_umum_(row, mIndex(ic,ivel)) -= tmp*normal(ivel);
      }
    }
  }

  //TODO: used more often?
  const double velint_normal_diff = velint_normal_m - velint_normal_s;

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    // (qm*n, km * um)
    // -(qm*n,km * u_DBC) for weak DBC or
    // -(qm*n,km * us)
    rhC_um_(mPres(ir),0) += funct_m_timefacfac_km(ir)*velint_normal_diff;
  }


  if (!eval_coupling_) return;


  //-----------------------------------------------
  //    + (qm*n, km *(Dus))
  //-----------------------------------------------
  for (unsigned ic =0; ic<slave_nen_; ic++)
  {
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const unsigned row = mPres(ir);
      const double tmp = funct_s_m_dyad_timefacfac_km(ic,ir);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // -(qm*n, km * Dus)
        C_umus_(row, sIndex(ic,ivel)) += tmp*normal(ivel);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_AdjointConsistency_SlaveTerms(
  const double &                                velint_normal_m,              ///< velocity in normal direction
  const double &                                velint_normal_s,              ///< interface velocity in normal direction
  const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac * kappas
  const LINALG::Matrix<slave_nen_,nen_> &       funct_s_m_dyad_timefacfac_ks, ///< (funct_s^T * funct) * timefacfac *kappas
  const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac_ks, ///< (funct_s^T * funct_s) * timefacfac *kappas
  const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
)
{

  //-----------------------------------------------
  //    - (qs*n, ks *(Dum))
  //-----------------------------------------------
  for (unsigned ic =0; ic<nen_; ic++)
  {
    for (unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const unsigned row = sPres(ir);

      const double tmp = funct_s_m_dyad_timefacfac_ks(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // -(qs*n, ks* Dum)
        C_usum_(row,mIndex(ic,ivel)) -= tmp*normal(ivel);
      }
    }
  }

  //-----------------------------------------------
  //    + (qs*n, ks *(Dus))
  //-----------------------------------------------
  for (unsigned ic = 0; ic<slave_nen_; ic++)
  {
    for (unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const unsigned row = sPres(ir);

      const double tmp = funct_s_s_dyad_timefacfac_ks(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // +(qs*n, ks* Dus)
        C_usus_(row,sIndex(ic,ivel)) += tmp*normal(ivel);
      }
    }
  }

  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    // (qs*n,ks* um)
    rhC_us_(sPres(ir),0) += funct_s_timefacfac_ks(ir)*velint_diff_normal_;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_MasterTerms(
  const LINALG::Matrix<nsd_,nen_> &     derxy_m,                      ///< master deriv
  const LINALG::Matrix<nsd_,nsd_> &     vderxy_m,                     ///< master deriv^n
  const LINALG::Matrix<nen_,1> &        funct_m_viscm_timefacfac_km,  ///< funct_m*mu_m*timefacfac
  const LINALG::Matrix<slave_nen_,1> &  funct_s_viscm_timefacfac_km,  ///< funct_s*mu_m*timefacfac
  const LINALG::Matrix<nsd_,1> &        normal                        ///< normal vector
)
{
  // viscous consistency term

  /*                           \       /                   i      \
- |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
  \                            /       \                         */


  // here we use a non-optimal order to assemble the values into C_umum
  // however for this term we have to save operations
  for (unsigned ic =0; ic<nen_; ic++)
  {
    const double normal_deriv_tmp = half_normal_deriv_m_(ic);

    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      const double tmp_derxy_m = derxy_m(ivel,ic);
      for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
      {
        const unsigned col = mIndex(ic,jvel);

        double tmp = half_normal_(jvel)*tmp_derxy_m;
        if(ivel==jvel) tmp+=normal_deriv_tmp;

        for (unsigned ir = 0; ir<nen_; ir++)
        {
          C_umum_(mIndex(ir,ivel), col) -= funct_m_viscm_timefacfac_km(ir) * tmp;
        }


        if (!eval_coupling_) continue;


        for (unsigned ir = 0; ir<slave_nen_; ir++)
        {
          C_usum_(sIndex(ir,ivel), col) += funct_s_viscm_timefacfac_km(ir) * tmp;
        }
      }
    }
  }


  for (unsigned ir = 0; ir<nen_; ir++)
  {
    const double tmp_val = funct_m_viscm_timefacfac_km(ir);

    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      //-----------------------------------------------
      //    - (vm, (2*km*mum) *eps(Dum)*n)
      //-----------------------------------------------
      rhC_um_(mIndex(ir,ivel),0) += tmp_val * vderxy_m_normal_transposed_(ivel);
    }
  }


  if (!eval_coupling_) return;


  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    const double tmp_val = funct_s_viscm_timefacfac_km(ir);

    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      //-----------------------------------------------
      //    + (vs, (2*km*mum) *eps(Dum)*n)
      //-----------------------------------------------

      // rhs
      rhC_us_(sIndex(ir,ivel),0) -= tmp_val * vderxy_m_normal_transposed_(ivel);
    }
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_SlaveTerms(
  const LINALG::Matrix<nsd_,slave_nen_> & derxy_s,                      ///< slave shape function derivatives
  const LINALG::Matrix<nen_,1> &          funct_m_viscs_timefacfac_ks,  ///< funct_m*mu_m*timefacfac
  const LINALG::Matrix<slave_nen_,1> &    funct_s_viscs_timefacfac_ks,  ///< funct_s*mu_m*timefacfac
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
  // diagonal block (i,i): +/-2*ks*mus * ...
  /*       nsd_
   *       *---*
   *       \    dN                     dN
   *  N  *  *   --  * 0.5 * n_j + N  * --  * n_i * 0.5
   *       /    dxj                    dxi
   *       *---*
   *       j = 1
   */

  // off-diagonal block (i,j) : +/-2*ks*mus * ...
  /*       dN
   *  N  * -- * n_j * 0.5
   *       dxi
   */

  for (unsigned ic =0; ic<slave_nen_; ic++)
  {
    const double normal_deriv_tmp = half_normal_deriv_s_(ic);

    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      const double tmp_derxy_s = derxy_s(ivel,ic);

      for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
      {
        const unsigned col = sIndex(ic,jvel);

        double tmp = half_normal_(jvel)*tmp_derxy_s;

        if(ivel==jvel) tmp+=normal_deriv_tmp;

        for (unsigned ir = 0; ir<nen_; ir++)
        {
          //-----------------------------------------------
          //    - (vm, (2*ks*mus) *eps(Dus)*n)
          //-----------------------------------------------
          C_umus_(mIndex(ir,ivel), col) -= funct_m_viscs_timefacfac_ks(ir)*tmp;
        }

        //-----------------------------------------------
        //    + (vs, (2*ks*mus) *eps(Dus)*n)
        //-----------------------------------------------
        for (unsigned ir = 0; ir<slave_nen_; ir++)
        {
          // diagonal block
          C_usus_(sIndex(ir,ivel),col) += funct_s_viscs_timefacfac_ks(ir) * tmp;
        }
      }
    }
  }

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    const double tmp_val = funct_m_viscs_timefacfac_ks(ir);
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // rhs
      rhC_um_(mIndex(ir,ivel),0) += tmp_val * vderxy_s_normal_transposed_(ivel);
    }
  }

  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    const double tmp_val = funct_s_viscs_timefacfac_ks(ir);
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      // rhs
      rhC_us_(sIndex(ir,ivel),0) -= tmp_val * vderxy_s_normal_transposed_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_AdjointConsistency_MasterTerms(
  const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
  const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
  const LINALG::Matrix<nsd_,nen_> &        derxy_m_viscm_timefacfac_km,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
  const LINALG::Matrix<nen_,1> &           funct_m,                           ///< embedded element funct *mu*timefacfac
  const LINALG::Matrix<slave_nen_,1> &     funct_s,                           ///< embedded element funct *mu*timefacfac
  const LINALG::Matrix<nsd_,1> &           normal                             ///< normal vector
)
{
    /*                                \       /                             i   \
  - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
    \                                 /       \                                */
  // (see Burman, Fernandez 2009)
  // +1.0 symmetric
  // -1.0 antisymmetric

  normal_deriv_m_viscm_km_.MultiplyTN(derxy_m_viscm_timefacfac_km, half_normal_); //half_normal(k)*derxy_m(k,ic)*viscm*km

  // here we use a non-optimal order to assemble the values into C_umum
  // however for this term we have to save operations
  for (unsigned ir =0; ir<nen_; ir++)
  {
    const double normal_deriv_tmp = normal_deriv_m_viscm_km_(ir);

    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      const double tmp_derxy_m = derxy_m_viscm_timefacfac_km(jvel,ir);
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        const unsigned row = mIndex(ir,ivel);

        double tmp = half_normal_(ivel)*tmp_derxy_m;
        if(ivel==jvel) tmp+=normal_deriv_tmp;

        for (unsigned ic = 0; ic<nen_; ic++)
        {
          C_umum_(row, mIndex(ic,jvel)) -= funct_m(ic) * tmp;
        }


        if (!eval_coupling_) continue;


        for (unsigned ic = 0; ic<slave_nen_; ic++)
        {
          C_umus_(row, sIndex(ic,jvel)) += funct_s(ic) * tmp;
        }
      }
    }
  }


  static LINALG::Matrix<nsd_,nsd_> velint_diff_dyad_normal, velint_diff_dyad_normal_symm;
  velint_diff_dyad_normal.MultiplyNT(velint_diff_, normal);

  for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      velint_diff_dyad_normal_symm(ivel,jvel) = velint_diff_dyad_normal(ivel,jvel) + velint_diff_dyad_normal(jvel,ivel);
    }
  }

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      const double derxy_m_viscm_timefacfac_km_half_tmp = derxy_m_viscm_timefacfac_km(jvel,ir)*0.5;

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // rhs
        rhC_um_(mIndex(ir,ivel),0) += derxy_m_viscm_timefacfac_km_half_tmp * velint_diff_dyad_normal_symm(ivel,jvel);
      }
    }
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_AdjointConsistency_SlaveTerms(
  const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
  const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
  const LINALG::Matrix<nen_,1> &           funct_m,                           ///< embedded element funct *mu*timefacfac
  const LINALG::Matrix<slave_nen_,1> &     funct_s,                           ///< embedded element funct *mu*timefacfac
  const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s_viscs_timefacfac_ks,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
  const LINALG::Matrix<nsd_,1> &           normal
)
{
    /*                                \       /                             i   \
  - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
    \                                 /       \                                */
  // (see Burman, Fernandez 2009)
  // +1.0 symmetric
  // -1.0 antisymmetric

  // diagonal block (i,i): +/-2*km*mum * alpha * ...
  /*
   *       nsd_
   *       *---*
   *       \    dN                    dN
   *        *   --  * 0.5 * n_j *N +  --  * n_i * 0.5 *N
   *       /    dxj                   dxi
   *       *---*
   *       j = 1
   */

  // off-diagonal block (i,j) : +/-2*km*mum * alpha * ...
  /*   dN
   *   -- * n_i * 0.5 * N
   *   dxj
   */
  normal_deriv_s_viscs_ks_.MultiplyTN(derxy_s_viscs_timefacfac_ks, half_normal_); //half_normal(k)*derxy_m(k,ic)*viscm*km

  for (unsigned ir =0; ir<slave_nen_; ir++)
  {
    const double normal_deriv_tmp = normal_deriv_s_viscs_ks_(ir);

    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      const double tmp_derxy_s = derxy_s_viscs_timefacfac_ks(jvel,ir);
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        const unsigned row = sIndex(ir,ivel);

        double tmp = half_normal_(ivel)*tmp_derxy_s;
        if(ivel==jvel) tmp+=normal_deriv_tmp;

        for (unsigned ic = 0; ic<nen_; ic++)
        {
          C_usum_(row, mIndex(ic,jvel)) -= funct_m(ic) * tmp;
        }


        if (!eval_coupling_) continue;


        for (unsigned ic = 0; ic<slave_nen_; ic++)
        {
          C_usus_(row, sIndex(ic,jvel)) += funct_s(ic) * tmp;
        }
      }
    }
  }


  static LINALG::Matrix<nsd_,nsd_> velint_diff_dyad_normal, velint_diff_dyad_normal_symm;
  velint_diff_dyad_normal.MultiplyNT(velint_diff_, normal);

  for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      velint_diff_dyad_normal_symm(ivel,jvel) = velint_diff_dyad_normal(ivel,jvel) + velint_diff_dyad_normal(jvel,ivel);
    }
  }

  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      const double derxy_s_viscs_timefacfac_ks_half_tmp = derxy_s_viscs_timefacfac_ks(jvel,ir)*0.5;

      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // rhs
        rhC_us_(sIndex(ir,ivel),0) += derxy_s_viscs_timefacfac_ks_half_tmp * velint_diff_dyad_normal_symm(ivel,jvel);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_MasterRHS(
  const double &                            pres_m,                      ///< master pressure
  const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &            normal                       ///< normal vector
)
{
             /*                  \       /          i      \
          + |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
             \                   /       \                */

  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    const double tmp_val = pres_m*normal(ivel);
    //-----------------------------------------------
    //    + (vm, km *(Dpm)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      // -(v,p*n)
      rhC_um_(mIndex(ir,ivel),0) -= funct_m_timefacfac_km(ir)*tmp_val;
    }

    if (!eval_coupling_) continue;

    //-----------------------------------------------
    //    - (vs, km *(Dpm)*n)
    //-----------------------------------------------
    for(unsigned ir = 0; ir<slave_nen_; ir++)
    {
      // -(v,p*n)
      rhC_us_(sIndex(ir,ivel),0) += funct_s_timefacfac_km(ir)*tmp_val;
    }
  } // end loop over velocity components
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_SlaveRHS(
  const double &                                pres_s,                       ///< slave pressure
  const LINALG::Matrix<nen_,1> &                funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
)
{
  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    const double tmp_val = pres_s*normal(ivel);
    //-----------------------------------------------
    //    + (vm, ks *(Dps)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      // -(vm, ks * ps*n)
      rhC_um_(mIndex(ir,ivel),0) -= funct_m_timefacfac_ks(ir)*tmp_val;
    }

    //-----------------------------------------------
    //    - (vs, ks *(Dps)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<slave_nen_; ir++)
    {
      // +(vs,ks * ps*n)
      rhC_us_(sIndex(ir,ivel),0) += funct_s_timefacfac_ks(ir)*tmp_val;
    }
  } // end loop over velocity components
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_AdjointConsistency_MasterRHS(
  const double &                          velint_normal_m,              ///< velocity in normal direction
  const double &                          velint_normal_s,              ///< interface velocity in normal direction
  const LINALG::Matrix<nen_,1> &          funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
      /*                   \     /              i   \
   - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
      \                    /     \                 */

  //-----------------------------------------------
  //    - (qm*n, km *(Dum))
  //-----------------------------------------------
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    // (qm*n, km * um)
    // -(qm*n,km * u_DBC) for weak DBC or
    // -(qm*n,km * us)
    rhC_um_(mPres(ir),0) += funct_m_timefacfac_km(ir)*velint_diff_normal_;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_AdjointConsistency_SlaveRHS(
  const double &                                velint_normal_m,              ///< velocity in normal direction
  const double &                                velint_normal_s,              ///< interface velocity in normal direction
  const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks         ///< funct_s * timefacfac * kappas
)
{
  //-----------------------------------------------
  //    - (qs*n, ks *(Dum))
  //    + (qs*n, ks *(Dus))
  //-----------------------------------------------
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    // (qs*n,ks* um)
    rhC_us_(sPres(ir),0) += funct_s_timefacfac_ks(ir)*velint_diff_normal_;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_MasterRHS(
  const LINALG::Matrix<nen_,1> &        funct_m_viscm_timefacfac_km,  ///< funct_m*mu_m*timefacfac
  const LINALG::Matrix<slave_nen_,1> &  funct_s_viscm_timefacfac_km,  ///< funct_s*mu_m*timefacfac
  const LINALG::Matrix<nsd_,1> &        normal                        ///< normal vector
)
{
  // viscous consistency term

  /*                           \       /                   i      \
- |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
  \                            /       \                         */



  // diagonal terms (i,i): +/-2*km*mum * ...
  /*
   *        nsd_
   *       *---*
   *       \    dN                     dN
   *  N  *  *   --  * 0.5 * n_j + N  * --  * n_i * 0.5
   *       /    dxj                    dxi
   *       *---*
   *       j = 1
   *
   * REMARK: 2nd summand == 'off-diagonal'-term for i=j
   */

  // off-diagonal terms (i,j) : +/-2*km*mum * ...
  /*       dN
   *  N  * -- * n_j * 0.5
   *       dxi
   */

  //-----------------------------------------------
  //    - (vm, (2*km*mum) *eps(Dum)*n)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    const double tmp_val = funct_m_viscm_timefacfac_km(ir);
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      rhC_um_(mIndex(ir,ivel),0) += tmp_val * vderxy_m_normal_transposed_(ivel);
    }
  }
  if (!eval_coupling_) return;

  //-----------------------------------------------
  //    + (vs, (2*km*mum) *eps(Dum)*n)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    const double tmp_val = funct_s_viscm_timefacfac_km(ir);
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      rhC_us_(sIndex(ir,ivel),0) -= tmp_val * vderxy_m_normal_transposed_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_SlaveRHS(
  const LINALG::Matrix<nsd_,nsd_> &       vderxy_s,                     ///< slave velocity gradient
  const LINALG::Matrix<nen_,1> &          funct_m_viscs_timefacfac_ks,  ///< funct_m*mu_m*timefacfac
  const LINALG::Matrix<slave_nen_,1> &    funct_s_viscs_timefacfac_ks,  ///< funct_s*mu_m*timefacfac
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
  // diagonal block (i,i): +/-2*ks*mus * ...
  /*       nsd_
   *       *---*
   *       \    dN                     dN
   *  N  *  *   --  * 0.5 * n_j + N  * --  * n_i * 0.5
   *       /    dxj                    dxi
   *       *---*
   *       j = 1
   */

  // off-diagonal block (i,j) : +/-2*ks*mus * ...
  /*       dN
   *  N  * -- * n_j * 0.5
   *       dxi
   */

  //-----------------------------------------------
  //    - (vm, (2*ks*mus) *eps(Dus)*n)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      rhC_um_(mIndex(ir,ivel),0) += funct_m_viscs_timefacfac_ks(ir) * vderxy_s_normal_transposed_(ivel);
    }
  }

  //-----------------------------------------------
  //    + (vs, (2*ks*mus) *eps(Dus)*n)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      rhC_us_(sIndex(ir,ivel),0) -= funct_s_viscs_timefacfac_ks(ir) * vderxy_s_normal_transposed_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_AdjointConsistency_MasterRHS(
  const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
  const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
  const LINALG::Matrix<nsd_,nen_> &        derxy_m_viscm_timefacfac_km,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
  const LINALG::Matrix<nsd_,1> &           normal                             ///< normal vector
)
{
    /*                                \       /                             i   \
  - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
    \                                 /       \                                */
  // (see Burman, Fernandez 2009)
  // +1.0 symmetric
  // -1.0 antisymmetric


  // diagonal block (i,i): +/-2*km*mum * alpha * ...
  /*       nsd_
   *       *---*
   *       \    dN                    dN
   *        *   --  * 0.5 * n_j *N +  --  * n_i * 0.5 * N
   *       /    dxj                   dxi
   *       *---*
   *       j = 1
   */

  // off-diagonal block (i,j) : +/-2*km*mum * alpha * ...
  /*       dN
   *  N  * -- * n_i * 0.5
   *       dxj
   */

  //-----------------------------------------------
  //    - alpha * ((2*km*mum) *eps(vm)*n , um)
  //-----------------------------------------------
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      const double tmp = derxy_m_viscm_timefacfac_km(ivel,ir) * 0.5;
      for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
      {

        rhC_um_(mIndex(ir,ivel),0) += tmp * (normal(jvel) * velint_diff_(ivel) + normal(ivel)*velint_diff_(jvel));
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_AdjointConsistency_SlaveRHS(
  const LINALG::Matrix<nsd_,1>&            velint_m,                          ///< velocity
  const LINALG::Matrix<nsd_,1>&            velint_s,                          ///< interface velocity
  const LINALG::Matrix<nsd_,slave_nen_> &  derxy_s_viscs_timefacfac_ks,       ///< master shape function derivatives * timefacfac * 2 * mu_m * kappa_m
  const LINALG::Matrix<nsd_,1> &           normal
)
{
    /*                                \       /                             i   \
  - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
    \                                 /       \                                */
  // (see Burman, Fernandez 2009)
  // +1.0 symmetric
  // -1.0 antisymmetric
  // diagonal block (i,i): +/-2*km*mum * alpha * ...
  /*
   *       nsd_
   *       *---*
   *       \    dN                    dN
   *        *   --  * 0.5 * n_j *N +  --  * n_i * 0.5 *N
   *       /    dxj                   dxi
   *       *---*
   *       j = 1
   */

  // off-diagonal block (i,j) : +/-2*km*mum * alpha * ...
  /*   dN
   *   -- * n_i * 0.5 * N
   *   dxj
   */
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      const unsigned row = sIndex(ir,ivel);
      for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
      {
        rhC_us_(row,0) += derxy_s_viscs_timefacfac_ks(jvel,ir) * 0.5 * (normal(jvel) * velint_diff_(ivel) + normal(ivel)*velint_diff_(jvel));
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty(
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,slave_nen_>&  funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
  const double &                                stabfac                       ///< stabilization factor
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_Stab_Penalty");

  // viscous stability term

  // combined viscous and inflow stabilization for one-sided problems (XFSI)
  // gamma_combined = max(alpha*mu/hk, |u*n| )
   /*                      _        \        /                     i   _   \
  |  gamma_combined *  v , u - u     | =  - |   gamma/h_K *  v , (u  - u)   |
   \                                /        \                            */

  // just viscous stabilization for two-sided problems (XFF, XFFSI)
  /*                                  \        /                           i   \
 |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
  \                                   /        \                              */

  funct_m_m_dyad_timefacfac_stabfac_.Update(stabfac, funct_m_m_dyad_timefacfac, 0.0);

  // + gamma*mu/h_K (vm, um))
  for (unsigned ic=0; ic<nen_; ic++)
  {
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const double tmp_val = funct_m_m_dyad_timefacfac_stabfac_(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ivel ++)
      {
        C_umum_(mIndex(ir,ivel), mIndex(ic,ivel)) += tmp_val;
      }
    }
  }

  for (unsigned ir=0; ir<nen_; ir++)
  {
    const double tmp = funct_m_timefacfac(ir);

    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      rhC_um_(mIndex(ir,ivel),0) -= tmp*velint_diff_stabfac_(ivel);

      // +(stab * vm, u_DBC) (weak dirichlet case) or from
      // +(stab * vm, u_s)
    }
  }


  if (!eval_coupling_) return;


  funct_s_m_dyad_timefacfac_stabfac_.Update(stabfac, funct_s_m_dyad_timefacfac, 0.0);

  // - gamma*mu/h_K (vm, us))
  // - gamma*mu/h_K (vs, um))

  for (unsigned ic=0; ic<slave_nen_; ic++)
  {
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const double tmp = funct_s_m_dyad_timefacfac_stabfac_(ic,ir);

      for (unsigned ivel = 0; ivel < nsd_; ivel ++)
      {
        C_umus_(mIndex(ir,ivel),sIndex(ic,ivel)) -= tmp;
      }
    }
  }

  for (unsigned ic=0; ic<nen_; ic++)
  {
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const double tmp = funct_s_m_dyad_timefacfac_stabfac_(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ivel ++)
      {
        C_usum_(sIndex(ir,ivel),mIndex(ic,ivel)) -= tmp;
      }
    }
  }

  funct_s_s_dyad_timefacfac_stabfac_.Update(stabfac, funct_s_s_dyad_timefacfac, 0.0);

  for (unsigned ic=0; ic<slave_nen_; ic++)
  {
    // + gamma*mu/h_K (vs, us))
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const double tmp = funct_s_s_dyad_timefacfac_stabfac_(ir,ic);

      for (unsigned ivel = 0; ivel < nsd_; ivel ++)
      {
        C_usus_(sIndex(ir,ivel),sIndex(ic,ivel)) += tmp;
      }
    }
  }

  for (unsigned ir=0; ir<slave_nen_; ir++)
  {
    const double tmp = funct_s_timefacfac(ir);

    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      // +(stab * vs, um)
      // -(stab * vs, us)
      rhC_us_(sIndex(ir,ivel),0) += tmp*velint_diff_stabfac_(ivel);
    }
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty_MasterTerms(
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const double &                                stabfac
)
{
  for (unsigned ic=0; ic<nen_; ic++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      const unsigned col = mIndex(ic,ivel);
      // + gamma*mu/h_K (vm, um))
      for (unsigned ir=0; ir<nen_; ir++)
      {
        C_umum_(mIndex(ir,ivel),col) += funct_m_m_dyad_timefacfac(ir,ic)*stabfac;
      }
    }
  }

  // + gamma*mu/h_K (vm, um))
  for (unsigned ir=0; ir<nen_; ir++)
  {
    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      // +(stab * vm, u_DBC) (weak dirichlet case) or from
      // +(stab * vm, u_s)
       rhC_um_(mIndex(ir,ivel),0) -= funct_m_timefacfac(ir)*velint_diff_stabfac_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty_MasterRHS(
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const double &                                stabfac
)
{
  // + gamma*mu/h_K (vm, um))
  for (unsigned ir=0; ir<nen_; ir++)
  {
    const double tmp_val = funct_m_timefacfac(ir)*stabfac;
    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      rhC_um_(mIndex(ir,ivel),0) -= tmp_val*velint_diff_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty_SlaveRHS(
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const double &                                stabfac                      ///< stabilization factor
)
{
  // - gamma*mu/h_K (vm, us))
  // - gamma*mu/h_K (vs, um))
  for (unsigned ir=0; ir<slave_nen_; ir++)
  {
    const double tmp_val = funct_s_timefacfac(ir)*stabfac;
    for (unsigned ivel = 0; ivel < nsd_; ivel ++)
    {
      rhC_us_(sIndex(ir,ivel),0) += tmp_val * velint_diff_(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 * add averaged term to balance instabilities due to convective
 * mass transport across the fluid-fluid interface
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Inflow_AveragedTerm(
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
  const LINALG::Matrix<nsd_,1>&                 normal,                       ///< normal vector n^m
  const double &                                density                       ///< fluid density
)
{
  //
  /*                                        \
 -|  [rho * (beta * n)] *  { v }_m , [   u ] |
  \ ----stab_avg-----                       / */

  // { v }_m = 0.5* (v^b + v^e) leads to the scaling with 0.5;
  // beta: convective velocity, currently beta=u^b_Gamma;
  // n:= n^b
  const double stabfac_avg_scaled = 0.5 * velint_m.Dot(normal) * density;

  for (unsigned ivel = 0; ivel<nsd_; ivel ++)
  {
    //  [rho * (beta * n^b)] (0.5*vb,ub)
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const unsigned mrow = mIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        C_umum_(mrow,mIndex(ic,ivel)) -= funct_m_m_dyad_timefacfac(ir,ic)*stabfac_avg_scaled;
      }

      const double tmp = funct_m_timefacfac(ir)*stabfac_avg_scaled;

      rhC_um_(mrow,0) += tmp*velint_m(ivel);

    //  -[rho * (beta * n^b)] (0.5*vb,ue)
      for (unsigned ic=0; ic<slave_nen_; ic++)
      {
        C_umus_(mrow,sIndex(ic,ivel)) += funct_s_m_dyad_timefacfac(ic,ir)*stabfac_avg_scaled;
      }

      rhC_um_(mrow,0) -= tmp*velint_s(ivel);
    }

    //  [rho * (beta * n^b)] (0.5*ve,ub)
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const unsigned srow = sIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        C_usum_(srow,mIndex(ic,ivel)) -= funct_s_m_dyad_timefacfac(ir,ic)*stabfac_avg_scaled;
      }

      const double tmp = funct_s_timefacfac(ir)*stabfac_avg_scaled;

      rhC_us_(srow,0) += tmp*velint_m(ivel);

    //-[rho * (beta * n^b)] (0.5*ve,ue)

      for (unsigned ic=0; ic<slave_nen_; ic++)
      {
        C_usus_(srow,sIndex(ic,ivel)) += funct_s_s_dyad_timefacfac(ir,ic)*stabfac_avg_scaled;
      }

      rhC_us_(srow,0) -= tmp*velint_s(ivel);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 * add averaged term to balance instabilities due to convective
 * mass transport across the fluid-fluid interface (only RHS)
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Inflow_AveragedTermRHS(
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const LINALG::Matrix<nsd_,1>&                 normal,                       ///< normal vector n^m
  const double &                                density                       ///< fluid density
)
{
  //
  /*                                        \
 -|  [rho * (beta * n)] *  { v }_m , [   u ] |
  \ ----stab_avg-----                      */

  // { v }_m = 0.5* (v^b + v^e) leads to the scaling with 0.5;
  // beta: convective velocity, currently beta=u^b_Gamma;
  // n:= n^b
  const double stabfac_avg_scaled = 0.5 * velint_m.Dot(normal) * density;

  for (unsigned ivel = 0; ivel<nsd_; ivel ++)
  {
    //  [rho * (beta * n^b)] (0.5*vb,ub)
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const unsigned mrow = mIndex(ir,ivel);

      const double tmp = funct_m_timefacfac(ir)*stabfac_avg_scaled;

      //  -[rho * (beta * n^b)] (0.5*vb,ue)
      rhC_um_(mrow,0) += tmp*velint_diff_(ivel);
    }

    //  [rho * (beta * n^b)] (0.5*ve,ub)
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const unsigned srow = sIndex(ir,ivel);

      const double tmp = funct_s_timefacfac(ir)*stabfac_avg_scaled;

      //-[rho * (beta * n^b)] (0.5*ve,ue)
      rhC_us_(srow,0) += tmp*velint_diff_(ivel);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void HybridLMCoupling<distype,slave_distype,slave_numdof>::MHCS_buildCouplingMatrices(
  const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
  const double &                                                fac,        ///< integration factor
  const LINALG::Matrix<nen_,1> &                                funct,      ///< shape function
  LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s,      ///< block rhs vector \f$ rhs_{\sigma} \f$
  const LINALG::Matrix<nsd_,1> &                                ivelint_jump,  ///< prescribed interface velocity or interface jump height
  const LINALG::Matrix<nsd_,1> &                                itraction_jump ///< prescribed interface traction or interface jump height
)
{
  LINALG::Matrix<nen_,slave_nen_> bK_ms;
  LINALG::Matrix<slave_nen_,nen_> bK_sm;

  // interface velocity at gauss-point (current gauss-point in calling method)
  LINALG::Matrix<nsd_,1> velint_s(true);
  this->GetInterfaceVelnp(velint_s);

  // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump height for coupled problems
  velint_s.Update(1.0, ivelint_jump, 1.0);

  // get nodal shape function vector
  LINALG::Matrix<slave_nen_,1> slave_funct(true);
  this->GetSlaveFunct(slave_funct);

  bK_ms.MultiplyNT(funct,slave_funct);
  bK_sm.UpdateT(bK_ms);

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      const double tmp = fac*normal(jvel);

      const unsigned sigma = stressIndex(ivel,jvel);
      // G_sus
      BG_sus_( sigma, ivel )->Update( tmp, bK_ms, 1.0 );
      rhs_s( sigma, 0 )->Update( -tmp*velint_s(ivel), funct, 1.0 );

      // G_uss
      BG_uss_( ivel, sigma )->Update( tmp, bK_sm, 1.0 );
    }
  }

  const double km=1.0; // only master-sided weighting
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_km;        ///< funct_s * timefacfac *kappa_m
  funct_s_timefacfac_km.Update(km, slave_funct, 0.0);

  // Traction Standard Consistency term
  MH_Traction_Consistency_Term( funct_s_timefacfac_km, itraction_jump );

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void HybridLMCoupling<distype,slave_distype,slave_numdof>::MHVS_buildCouplingMatrices(
  const LINALG::Matrix<nsd_,1> &                                  normal,     ///< normal vector
  const double &                                                  fac,        ///< integration factor
  const LINALG::Matrix<nen_,1> &                                  funct,      ///< background element shape functions
  LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> &   rhs_s,      ///< block rhs vector \f$ rhs_{\sigma}\f$
  const double &                                                  press,      ///< background element pressure
  LINALG::Matrix<nen_,1> &                                        rhs_pmus,   ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
  const LINALG::Matrix<nsd_,1> &                                  ivelint_jump,  ///< prescribed interface velocity or interface jump height
  const LINALG::Matrix<nsd_,1> &                                  itraction_jump ///< prescribed interface traction or interface jump height
)
{
  // interface velocity at gauss-point (current gauss-point in calling method)
  LINALG::Matrix<nsd_,1> velint_s(true);
  this->GetInterfaceVelnp(velint_s);

  // add the prescribed interface velocity for weak Dirichlet boundary conditions or the jump height for coupled problems
  velint_s.Update(1.0, ivelint_jump, 1.0);

  // block submatrices for interface coupling; s: slave side, m: master side (always background here)
  LINALG::Matrix<nen_, slave_nen_> bG_ms(true);
  LINALG::Matrix<slave_nen_, nen_> bG_sm(true);

  LINALG::Matrix<slave_nen_,1> slave_funct;
  this->GetSlaveFunct(slave_funct);

  bG_ms.MultiplyNT(funct, slave_funct);
  bG_sm.MultiplyNT(slave_funct, funct);

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      /*
       * G_sus
       *
       * from:
       *   /              \
       *  | \tau^m n, u^s  |
       *   \              /
       */
      BG_sus_(stressIndex(ivel,jvel), ivel)->Update(fac * normal(jvel), bG_ms, 1.0);
      rhs_s(stressIndex(ivel,jvel), 0)->Update(-fac * normal(jvel) * velint_s(ivel), funct, 1.0);

      /*
       *  G_uss
       *
       *  from:
       *   /                 \
       *  | v^s, \sigma^m n   |
       *   \                 /
       *
       */
      BG_uss_(ivel, stressIndex(ivel,jvel))->Update(fac * normal(jvel), bG_sm, 1.0);
    }

    //Build cross-interface pressure-velocity coupling matrices G_uip, G_pui!

    //G_pmus - from adjoint pressure consistency term
    /*
     *  /            \
     *  | q^m, u^s n  |
     *  \            /
     *
     */

    BG_pmus_(0,ivel)->Update(fac * normal(ivel), bG_ms, 1.0);


    //G_uspm - from pressure consistency term
    /*
     *  /           \
     *  | -v^s, p n  |
     *  \           /
     *
     */

    BG_uspm_(ivel,0)->Update(-fac * normal(ivel), bG_sm, 1.0);
  }

  // add normal interface velocity to rhs vector (pressure row)
  const double svelnormal = velint_s.Dot(normal);

  // (Viscous) stress/Velocities

  // rhs_pm_us, residual from:
  /*
   *  /            \
   *  | q^m, u^s n  |
   *  \            /
   *
   */
  rhs_pmus.Update(-fac * svelnormal, funct, 1.0);

  // rhs_us_pm, residual from:
  /*
   *  /            \
   *  | -v^s, p^m n |
   *  \            /
   *
   */
  // belongs to the side and therefore contributes to rhC_ui_!
  for (unsigned ir = 0; ir < slave_nen_; ++ir)
  {
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      rhC_us_(sIndex(ir,ivel), 0) += press * fac * normal(ivel) * slave_funct(ir);
    }

    rhC_us_(sIndex(ir,Pres), 0) = 0.0;
  }

  const double km=1.0; // only master-sided weighting
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_km;        ///< funct_s * timefacfac *kappa_m
  funct_s_timefacfac_km.Update(km, slave_funct, 0.0);

  // Traction Standard Consistency term
  MH_Traction_Consistency_Term( funct_s_timefacfac_km, itraction_jump );


  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void HybridLMCoupling<distype,slave_distype,slave_numdof>::MH_Traction_Consistency_Term(
  const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &            itraction_jump                ///< prescribed interface traction, jump height for coupled problems
)
{
        /*            \
     - |  < v >,   t   |   with t = [sigma * n]
        \             /     */

  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    /*
    //-----------------------------------------------
    //    - (vm, ks * t) = 0 as ks=0
    //-----------------------------------------------
    */

    //-----------------------------------------------
    //    + (vs, km * t)
    //-----------------------------------------------
    for(unsigned ir = 0; ir<slave_nen_; ir++)
    {
      const double funct_s_km_timefacfac_traction = funct_s_timefacfac_km(ir)*itraction_jump(ivel);

      const unsigned row = sIndex(ir,ivel);
      rhC_us_(row,0) += funct_s_km_timefacfac_traction;
    }
  } // end loop over velocity components
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void HybridLMCoupling<distype,slave_distype,slave_numdof>::HybridLM_buildFinalCouplingMatrices(
  LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,numstressdof_> &     BinvK_ss,    ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
  LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,master_numdof_,numstressdof_> &    BKumsInvKss, ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
  LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,master_numdof_> &    BK_sum,      ///< block matrix \f$ K_{\sigma u} \f$
  LINALG::BlockMatrix<LINALG::Matrix<nen_,1>, numstressdof_,1>      &              rhs_s        ///< block rhs vector \f$ rhs_{\sigma}\f$
)
{
  // final coupling matrices in block form
  LINALG::BlockMatrix<LINALG::Matrix<nen_,slave_nen_>,master_numdof_,nsd_> BCumus;
  LINALG::BlockMatrix<LINALG::Matrix<slave_nen_,nen_>,nsd_,master_numdof_> BCusum;

  // auxiliary matrices for intermediate calculation results of the matrix product
  LINALG::BlockMatrix<LINALG::Matrix<slave_nen_,nen_>,nsd_,numstressdof_>  BGussInvKss;
  LINALG::BlockMatrix<LINALG::Matrix<slave_nen_,1>,nsd_,1>                 BGussinvKssrhs_s;

  // BKusInvKss:
  // (K_ums + G_ums) * K_ss^-1 (MHVS) or
  // G_ums * K_ss^-1 (MHCS)

  // (K_ums + G_ums) * K_ss^-1 G_sus (MHVS) or
  // G_ums * K_ss^-1 G_sus (MHCS)
  BCumus.Multiply(BKumsInvKss,BG_sus_);

  // G_uss * K_ss^-1
  BGussInvKss.Multiply(BG_uss_,BinvK_ss);

  // G_uss * K_ss^-1 * (K_sum + G_sum) (MHVS) or
  // G_uss * K_ss^-1 * (K_sum + G_sum + K_spm) (MHCS)
  BCusum.Multiply(BGussInvKss,BK_sum);

  // G_uss K_ss^-1 rhs_s
  BGussinvKssrhs_s.Multiply(BGussInvKss,rhs_s);

  // transfer the entries from Cumus,Cusum, rhCus (in case of MHCS, coupling term us-p is
  // included in Cusum!)

  // loop over slave velocity dof
  for (unsigned isvel = 0; isvel < nsd_; ++ isvel)
  {
    // loop over background element dof (velocity & pressure)
    for (unsigned imdof = 0; imdof < master_numdof_; ++ imdof)
    {
      // (um-us)
      if (BCumus.IsUsed(imdof,isvel))
      {
        const LINALG::Matrix<nen_,slave_nen_> & bCumus = *BCumus(imdof,isvel);
        // loop over slave element nodes
        for (unsigned isn = 0; isn < slave_nen_; ++ isn)
        {
          // loop over background element nodes
          for (unsigned imn = 0; imn < nen_; ++ imn)
          {
            C_umus_(mIndex(imn,imdof),sIndex(isn,isvel)) -= bCumus(imn,isn);
          }
        }
      } // (um-us)

      // (us-um), MHCS: (us-pm)
      if (BCusum.IsUsed(isvel,imdof))
      {
        const LINALG::Matrix<slave_nen_,nen_> & bCusum = *BCusum(isvel,imdof);
        // loop over slave element nodes
        for (unsigned isn = 0; isn < slave_nen_; ++ isn)
        {
          // loop over background element nodes
          for (unsigned imn = 0; imn < nen_; ++ imn)
          {
            C_usum_(sIndex(isn,isvel),mIndex(imn,imdof)) -= bCusum(isn,imn);
          }
        }
      } // (us-um), MHCS: (us-pm)
    }

    // add surface-based pressure coupling terms (only MHVS)

    // (us-pm)
    if (BG_uspm_.IsUsed(isvel,0))
    {
      const LINALG::Matrix<slave_nen_,nen_> & bGuspm = *BG_uspm_(isvel,0);
      // loop over slave element nodes
      for (unsigned isn = 0; isn < slave_nen_; ++ isn)
      {
        // loop over background element nodes
        for (unsigned imn = 0; imn < nen_; ++ imn)
        {
          C_usum_(sIndex(isn,isvel),mIndex(imn,Pres)) = bGuspm(isn,imn);
        }
      }
    } // (us-pm)

    // (pm-us)
    if (BG_pmus_.IsUsed(0,isvel))
    {
      const LINALG::Matrix<nen_,slave_nen_> & bGpmus = *BG_pmus_(0,isvel);
      // loop over slave element nodes
      for (unsigned isn = 0; isn < slave_nen_; ++ isn)
      {
        // loop over background element nodes
        for (unsigned imn = 0; imn < nen_; ++ imn)
        {
          C_umus_(mIndex(imn,Pres),sIndex(isn,isvel)) = bGpmus(imn,isn);
        }
      }
    } // (pm-us)

    // rhs - us
    if (BGussinvKssrhs_s.IsUsed(isvel,0))
    {
      const LINALG::Matrix<slave_nen_,1> & bGussinvKssrhs_s = *BGussinvKssrhs_s(isvel,0);

      // loop over slave element nodes
      for (unsigned isn = 0; isn < slave_nen_; ++ isn)
      {
        rhC_us_(sIndex(isn,isvel),0) -= bGussinvKssrhs_s(isn,0);
      }
    } // rhs - us
  } // end loop over slave velocity dof

  // finally, build G_uss & G_sus for C_usus

  // G_sus
  // loop over block rows
  for (unsigned ibr = 0; ibr < numstressdof_; ++ibr)
  {
    // loop over block columns (interface velocity)
    for (unsigned ibc = 0; ibc < nsd_; ++ibc)
    {
      // extract the stress-velocity coupling submatrix
      if (BG_sus_.IsUsed(ibr,ibc))
      {
        LINALG::Matrix<nen_,slave_nen_> & bGsus = * BG_sus_(ibr,ibc);

        // transfer the entries
        for (unsigned ir =0; ir<nen_; ++ir)
        {
          unsigned stressrow = ibr + ir * numstressdof_;

          for (unsigned ic = 0; ic < slave_nen_; ++ic)
          {
            unsigned slavevelcol = ibc + ic * slave_numdof;

            G_sus_(stressrow,slavevelcol) = bGsus(ir,ic);
          }
        }
      }
    }
  } // G_sus filled

  // fill G_uss_ from BG_uss_

  // loop over block columns (sigmaxx, sigmaxy, ...)
  for (unsigned ibc = 0; ibc < numstressdof_; ++ibc)
  {
    // loop over block rows (interface velocity)
    for (unsigned ibr = 0; ibr < nsd_; ++ibr)
    {
      if (BG_uss_.IsUsed(ibr,ibc))
      {
        LINALG::Matrix< slave_nen_, nen_> & bGuss = * BG_uss_(ibr,ibc);

        // transfer the entries
        for (unsigned ic=0; ic<nen_; ++ic)
        {
          unsigned stresscol = ibc + ic * numstressdof_;

          for (unsigned ir=0; ir<slave_nen_; ++ir)
          {
            unsigned slavevelrow = ibr + ir * slave_numdof;

            G_uss_(slavevelrow,stresscol) = bGuss(ir,ic);
          }
        }
      }
    }
  } // G_uss filled

}

} // namespace XFLUID
} // namespace ELEMENTS
} // namespace DRT

// pairs with numdof=3
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad9,3>;

template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6,DRT::Element::dis_none,3>;


// pairs with numdof=4
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::quad9,4>;
//

//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tet4, 4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::hex27,4>;

template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex8,  DRT::Element::dis_none,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex20, DRT::Element::dis_none,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::hex27, DRT::Element::dis_none,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet4,  DRT::Element::dis_none,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::tet10, DRT::Element::dis_none,4>;
template class DRT::ELEMENTS::XFLUID::NitscheCoupling<DRT::Element::wedge6, DRT::Element::dis_none,4>;


// pairs with numdof=3
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad9,3>;

template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::dis_none,3>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6,DRT::Element::dis_none,3>;

// pairs with numdof=4
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex8,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex20, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::hex27, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet4,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::tet10, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::HybridLMCoupling<DRT::Element::wedge6, DRT::Element::quad9,4>;

// pairs with numdof=3
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad9,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,  DRT::Element::tri3,3>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6,  DRT::Element::tri6,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad9,3>;

// pairs with numdof=4
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::quad9,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tri3,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tri6,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::quad9,4>;

//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tet4, 4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex8,  DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex20, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::hex27, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet4,  DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::tet10, DRT::Element::hex27,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tet4,4>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::tet10,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::hex8,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::hex20,4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DRT::Element::wedge6, DRT::Element::hex27,4>;
