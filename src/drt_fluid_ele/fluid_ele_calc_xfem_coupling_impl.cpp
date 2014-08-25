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

#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"

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
  const std::string            state,        ///< state
  const std::vector<int>&      lm            ///< local map
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(state);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

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
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::SetSlaveVel(
  const DRT::Discretization &  slavedis,       ///< coupling slave discretization
  const std::string            state,          ///< state
  const std::vector<int>&      lm              ///< local map
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = slavedis.GetState(state);
  if (matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

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
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::Evaluate(
  const LINALG::Matrix<nsd_-1,1>  & eta,     ///< local coordinates w.r.t slave element
  LINALG::Matrix<nsd_,1>          & x,       ///< global coordinates of gaussian point
  LINALG::Matrix<nsd_,1>          & normal,  ///< normal vector
  double                          & drs      ///< transformation factor
)
{
  // Todo: remove this routine - it goes together with old boundary transformation
  if (slave_nsd_ == nsd_)
  {
    dserror("You called 2D evaluation routine when coupling with a 3D element.");
    return;
  }

  LINALG::Matrix<slave_nsd_,slave_nsd_> metrictensor;
  LINALG::Matrix<slave_nsd_+1,slave_nen_> bele_xyze(true);
  LINALG::Matrix<slave_nsd_,slave_nen_> bele_deriv(true);
  // auxiliary variable to match the interface for metric
  // tensor computation
  LINALG::Matrix<slave_nsd_+1,1> aux(true);
  for (unsigned i = 0; i < slave_nsd_+1; ++ i)
  {
    if (i < nsd_)
    {
        aux(i) = normal(i);
      for (unsigned j = 0; j < slave_nen_; ++ j)
        bele_xyze(i,j) = slave_xyze_(i,j);
    }
  }

  DRT::UTILS::shape_function_2D( slave_funct_, eta( 0 ), eta( 1 ), slave_distype );
  DRT::UTILS::shape_function_2D_deriv1( bele_deriv, eta( 0 ), eta( 1 ), slave_distype );
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<slave_distype>( bele_xyze, bele_deriv, metrictensor, drs, &aux);
  aux.Clear();
  aux.Multiply(bele_xyze, slave_funct_);

  for (unsigned i = 0; i < slave_nsd_+1; ++ i)
  {
    if (i < nsd_)
      x(i) = aux(i);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void SlaveElementRepresentation<distype,slave_distype,slave_numdof>::Evaluate( LINALG::Matrix<nsd_,1> & xslave )
{
  if (slave_nsd_ != nsd_)
  {
    dserror("You called 3D evaluation routine when coupling with a 2D element.");
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
    for(unsigned idim=0; idim<nsd_; ++idim )
    {
      // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
      iforce[idim+inode*nsd_] += slave_funct_(inode) * traction(idim) * fac;
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
  LINALG::Matrix<nsd_-1,1> & xi_side       ///< local coordinates of projected gaussian point w.r.t side
)
{
  // check, if called on a 3D-element
  if (slave_nsd_ == nsd_)
  {
    dserror("You can't project onto a 3D coupling slave element directly. You need an associated boundary element!");
  }

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );
  // Initialization
  LINALG::Matrix<slave_nen_,1> funct(true);      // shape functions
  LINALG::Matrix<2,slave_nen_> deriv(true);      // derivatives dr, ds
  LINALG::Matrix<3,slave_nen_> deriv2(true);     // 2nd derivatives drdr, dsds, drds

  LINALG::Matrix<3,1> x(true);

  LINALG::Matrix<3,2> derxy (true);
  LINALG::Matrix<3,1> dx_dr (true);
  LINALG::Matrix<3,1> dx_ds (true);

  LINALG::Matrix<3,3> derxy2 (true);
  LINALG::Matrix<3,1> dx_drdr (true);
  LINALG::Matrix<3,1> dx_dsds (true);
  LINALG::Matrix<3,1> dx_drds (true);

  LINALG::Matrix<3,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<3,3> sysmat(true);               // matrix for the newton system
  LINALG::Matrix<3,1> incr(true);                 // increment of the newton system

  LINALG::Matrix<3,1> sol(true); // sol carries xi_1, xi_2, d (distance)

  if(slave_distype == DRT::Element::tri3 or
     slave_distype == DRT::Element::tri6)
  {
    sol(0) = 0.333333333333333;
    sol(1) = 0.333333333333333;
  }
  else if( slave_distype == DRT::Element::quad4 or
           slave_distype == DRT::Element::quad8 or
           slave_distype == DRT::Element::quad9)
  {

    sol(0) = 0.0;
    sol(1) = 0.0;
  }
  else
  {
    dserror("Define start side xi-coordinates for unsupported cell type");
  }


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
    DRT::UTILS::shape_function_2D( funct, sol( 0 ), sol( 1 ), slave_distype );
    DRT::UTILS::shape_function_2D_deriv1( deriv, sol( 0 ), sol( 1 ), slave_distype );
    DRT::UTILS::shape_function_2D_deriv2( deriv2, sol( 0 ), sol( 1 ), slave_distype );

    x.Multiply(slave_xyze_, funct);

    derxy.MultiplyNT(slave_xyze_, deriv);

    derxy2.MultiplyNT(slave_xyze_, deriv2);

    // set dx_dr and dx_ds
    for (unsigned i=0; i< 3; i++)
    {
      dx_dr(i) = derxy(i,0);
      dx_ds(i) = derxy(i,1);

      dx_drdr(i) = derxy2(i,0);
      dx_dsds(i) = derxy2(i,1);
      dx_drds(i) = derxy2(i,2);
    }

    // get vector products
    LINALG::Matrix<3,1> dx_drdr_times_dx_ds(true);
    LINALG::Matrix<3,1> dx_dr_times_dx_drds(true);
    LINALG::Matrix<3,1> dx_drds_times_dx_ds(true);
    LINALG::Matrix<3,1> dx_dr_times_dx_dsds(true);
    LINALG::Matrix<3,1> dx_dr_times_dx_ds(true);

    dx_drdr_times_dx_ds(0) = dx_drdr(1)*dx_ds(2)-dx_ds(1)*dx_drdr(2);
    dx_drdr_times_dx_ds(1) = dx_drdr(2)*dx_ds(0)-dx_ds(2)*dx_drdr(0);
    dx_drdr_times_dx_ds(2) = dx_drdr(0)*dx_ds(1)-dx_ds(0)*dx_drdr(1);

    dx_dr_times_dx_drds(0) = dx_dr(1)*dx_drds(2)-dx_drds(1)*dx_dr(2);
    dx_dr_times_dx_drds(1) = dx_dr(2)*dx_drds(0)-dx_drds(2)*dx_dr(0);
    dx_dr_times_dx_drds(2) = dx_dr(0)*dx_drds(1)-dx_drds(0)*dx_dr(1);

    dx_drds_times_dx_ds(0) = dx_drds(1)*dx_ds(2)-dx_ds(1)*dx_drds(2);
    dx_drds_times_dx_ds(1) = dx_drds(2)*dx_ds(0)-dx_ds(2)*dx_drds(0);
    dx_drds_times_dx_ds(2) = dx_drds(0)*dx_ds(1)-dx_ds(0)*dx_drds(1);

    dx_dr_times_dx_dsds(0) = dx_dr(1)*dx_dsds(2)-dx_dsds(1)*dx_dr(2);
    dx_dr_times_dx_dsds(1) = dx_dr(2)*dx_dsds(0)-dx_dsds(2)*dx_dr(0);
    dx_dr_times_dx_dsds(2) = dx_dr(0)*dx_dsds(1)-dx_dsds(0)*dx_dr(1);

    dx_dr_times_dx_ds(0) = dx_dr(1)*dx_ds(2)-dx_ds(1)*dx_dr(2);
    dx_dr_times_dx_ds(1) = dx_dr(2)*dx_ds(0)-dx_ds(2)*dx_dr(0);
    dx_dr_times_dx_ds(2) = dx_dr(0)*dx_ds(1)-dx_ds(0)*dx_dr(1);

    // define sysmat
    for(unsigned i=0; i< 3; i++)
    {
      // d/dr
      sysmat(i,0) = dx_dr(i) - sol(2) * (dx_drdr_times_dx_ds(i) + dx_dr_times_dx_drds(i));

      // d/ds
      sysmat(i,1) = dx_ds(i) - sol(2) * (dx_drds_times_dx_ds(i) + dx_dr_times_dx_dsds(i));

      // d/d(dist)
      sysmat(i,2) = - dx_dr_times_dx_ds(i);


      // residual
      residuum(i) = x(i) - sol(2) * dx_dr_times_dx_ds(i) - x_gp_lin(i);

    }

    sysmat.Invert();

    //solve Newton iteration
    incr.Clear();
    incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

    // update solution
    sol.Update(1.0, incr, 1.0);

    // check ° absolute criterion for local coordinates (between [-1,1]^2)
    //       ° absolute criterion for distance (-> 0)
    //       ° absolute criterion for whole residuum
    if( sqrt(incr(0)*incr(0)+incr(1)*incr(1)) <  absTolIncr
        && incr(2) < absTOLdist
        && residuum.Norm2() < absTolRes)
    {
      converged = true;
    }

  }

  if(!converged)
  {
    std::cout.precision(15);

    std::cout << "increment criterion loc coord "
        //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
        << sqrt(incr(0)*incr(0)+incr(1)*incr(1))
        << " \tabsTOL: " << absTolIncr
        << std::endl;
    std::cout << "absolute criterion for distance "
        << incr(2)
        << " \tabsTOL: " << absTOLdist
        << std::endl;
    std::cout << "relative criterion whole residuum "
        << residuum.Norm2()
        << " \tabsTOL: " << absTolRes
        << std::endl;


    std::cout << "sysmat.Invert" << sysmat << std::endl;
    std::cout << "sol-norm " << sol.Norm2() << std::endl;
    std::cout << "sol " << sol << std::endl;
    std::cout << "x_gp_lin" << x_gp_lin << std::endl;
    std::cout << "side " << slave_xyze_ << std::endl;

    dserror( "Newton scheme in ProjectOnSide not converged! " );
  }

  // evaluate shape function at solution
  DRT::UTILS::shape_function_2D( slave_funct_, sol( 0 ), sol( 1 ), slave_distype );

  // get projected gauss point
  x_side.Multiply(slave_xyze_, slave_funct_);

  // set local coordinates w.r.t side
  xi_side(0) = sol(0);
  xi_side(1) = sol(1);
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
void NitscheCoupling<distype,slave_distype,slave_numdof>::ApplyConvStabTerms(
  const Teuchos::RCP<SlaveElementInterface<distype> > & slave_ele,              ///< associated slave element coupling object
  Epetra_SerialDenseMatrix &                            C_umum,                 ///< standard master-master-matrix
  Epetra_SerialDenseVector &                            rhs_Cum,                ///< master rhs
  const LINALG::Matrix<nen_,1> &                        funct_m,                ///< master shape functions
  const LINALG::Matrix<nsd_,1> &                        velint_m,               ///< vector of slave shape functions
  const double &                                        NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
  const double &                                        avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
  const double &                                        timefacfac,             ///< theta*dt
  INPAR::XFEM::XFF_ConvStabScaling                      xff_conv_stab           ///< type of convective stabilization in XFF-problems
)
{
  if (applicationType_ == SlaveElementInterface<distype>::XFluidFluid &&
      xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_none)
    dserror("Cannot apply convective stabilization terms for XFF_ConvStabScaling_none!");

  // funct_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_timefacfac(funct_m);
  funct_m_timefacfac.Scale(timefacfac);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  LINALG::Matrix<nen_,nen_> funct_m_m_dyad_timefacfac(true);
  funct_m_m_dyad_timefacfac.MultiplyNT(funct_m_timefacfac, funct_m);

  // velint_s
  LINALG::Matrix<nsd_,1> velint_s;
  slave_ele->GetInterfaceVel(velint_s);

  // REMARK:
  // the (additional) convective stabilization is included in NIT_full_stab_fac
  // (in case of mixed/hybrid LM approaches, we don't compute the penalty term explicitly -
  // it 'evolves'); in that case we therefore don't chose the maximum, but add the penalty term scaled
  // with conv_stab_fac to the viscous counterpart;
  // this happens by calling NIT_Stab_Penalty

  switch (applicationType_)
  {
    case SlaveElementInterface<distype>::XFluidWDBC:
    {
      NIT_Stab_Penalty_MasterTerms(
        C_umum,
        rhs_Cum,
        velint_m,
        velint_s,
        funct_m_timefacfac,
        funct_m_m_dyad_timefacfac,
        NIT_full_stab_fac
      );
      break;
    }
    case SlaveElementInterface<distype>::XFluidFluid:
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

      if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty)
      {
        NIT_Stab_Penalty(
          C_umum,
          rhs_Cum,
          velint_m,
          velint_s,
          funct_m_timefacfac,
          funct_m_m_dyad_timefacfac,
          funct_s_timefacfac,
          funct_s_m_dyad_timefacfac,
          funct_s_s_dyad_timefacfac,
          NIT_full_stab_fac
        );
      }

      // prevent instabilities due to convective mass transport across the fluid-fluid interface
      if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_averaged)
      {
        NIT_Stab_ConvAveraged(
          C_umum,
          rhs_Cum,
          velint_m,
          velint_s,
          funct_m_timefacfac,
          funct_m_m_dyad_timefacfac,
          funct_s_timefacfac,
          funct_s_m_dyad_timefacfac,
          funct_s_s_dyad_timefacfac,
          avg_conv_stab_fac
        );
      }
      break;
    }
    case SlaveElementInterface<distype>::MonolithicXFSI:
    {
      dserror("Convective stabilization in monolithic XFSI is not yet available!");
      break;
    }
    default:
      dserror("Unsupported application type. Cannot apply convective stabilization terms."); break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_buildCouplingMatrices(
  Epetra_SerialDenseMatrix &        C_umum,                 ///< standard master-master-matrix
  Epetra_SerialDenseVector &        rhs_Cum,                ///< standard master-rhs
  const LINALG::Matrix<nsd_,1> &    normal,                 ///< outward pointing normal (defined by the coupling partner, that determines the interface traction)
  const double &                    timefacfac,             ///< theta*dt
  const double &                    visceff_m,              ///< viscosity in coupling master fluid
  const double &                    visceff_s,              ///< viscosity in coupling slave fluid
  const double &                    kappa_m,                ///< mortaring weight for coupling master
  const double &                    kappa_s,                ///< mortaring weight for coupling slave
  const double &                    NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
  const double &                    avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
  const LINALG::Matrix<nen_,1> &    funct_m,                ///< coupling master shape functions
  const LINALG::Matrix<nsd_,nen_> & derxy_m,                ///< spatial derivatives of coupling master shape functions
  const LINALG::Matrix<nsd_,nsd_> & vderxy_m,               ///< coupling master spatial velocity derivatives
  const double &                    press_m,                ///< coupling master pressure
  const LINALG::Matrix<nsd_,1> &    velint_m,               ///< coupling master interface velocity
  INPAR::XFEM::XFF_ConvStabScaling  xff_conv_stab           ///< type of convective stabilization in XFF-problems
)
{
//  TEUCHOS_FUNC_TIME_MONITOR("FLD::NIT_buildCouplingMatricess");

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

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  // interface velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> velint_s;
  this->GetInterfaceVel(velint_s);

  // funct_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_timefacfac(funct_m);
  funct_m_timefacfac.Scale(timefacfac);

  // funct_s * timefac * fac
  LINALG::Matrix<slave_nen_,1> funct_s;
  this->GetSlaveFunct(funct_s);
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac(funct_s);
  funct_s_timefacfac.Scale(timefacfac);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  LINALG::Matrix<nen_,nen_> funct_m_m_dyad_timefacfac;
  funct_m_m_dyad_timefacfac.MultiplyNT(funct_m_timefacfac, funct_m);

  // funct_s * timefac * fac * funct_s * kappa_s (dyadic product)
  LINALG::Matrix<slave_nen_,slave_nen_> funct_s_s_dyad_timefacfac;
  funct_s_s_dyad_timefacfac.MultiplyNT(funct_s_timefacfac, funct_s);

  // funct_s * timefac * fac * funct_m (dyadic product)
  LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac;
  funct_s_m_dyad_timefacfac.MultiplyNT(funct_s_timefacfac, funct_m);

  // 2 * mu_m * timefacfac
  const double viscm_fac = 2.0 * timefacfac * visceff_m;

  // compute normal velocity components
  const double velint_normal_m = velint_m.Dot(normal);
  const double velint_normal_s = velint_s.Dot(normal);

  // Todo: no convective stabilization terms for monolithic XFSI available yet

  //-----------------------------------------------------------------
  // viscous stability term
  // REMARK: this term includes also inflow coercivity in case of XFSI
  // with modified stabfac (see NIT_ComputeStabfac)

  switch (applicationType_)
  {
    case SlaveElementInterface<distype>::XFluidWDBC:
    {
      NIT_Stab_Penalty_MasterTerms(
        C_umum,
        rhs_Cum,
        velint_m,
        velint_s,
        funct_m_timefacfac,
        funct_m_m_dyad_timefacfac,
        NIT_full_stab_fac
      );
      break;
    }
    case SlaveElementInterface<distype>::XFluidFluid:
    {
      NIT_Stab_Penalty(
        C_umum,
        rhs_Cum,
        velint_m,
        velint_s,
        funct_m_timefacfac,
        funct_m_m_dyad_timefacfac,
        funct_s_timefacfac,
        funct_s_m_dyad_timefacfac,
        funct_s_s_dyad_timefacfac,
        NIT_full_stab_fac
      );

      if (xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_averaged ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty ||
          xff_conv_stab == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty)
      {
        NIT_Stab_ConvAveraged(
          C_umum,
          rhs_Cum,
          velint_m,
          velint_s,
          funct_m_timefacfac,
          funct_m_m_dyad_timefacfac,
          funct_s_timefacfac,
          funct_s_m_dyad_timefacfac,
          funct_s_s_dyad_timefacfac,
          avg_conv_stab_fac
        );
      }
      break;
    }
    default:
      break;
  }

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - standard Dirichlet case

  if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterTerms(
        C_umum,
        rhs_Cum,
        press_m,
        funct_m_timefacfac,
        funct_s_timefacfac,
        funct_s_m_dyad_timefacfac,
        funct_m_m_dyad_timefacfac,
        normal
        );

    //-----------------------------------------------------------------
    // pressure adjoint consistency term

    NIT_p_AdjointConsistency_MasterTerms(
        C_umum,
        rhs_Cum,
        velint_normal_m,
        velint_normal_s,
        funct_m_timefacfac,
        funct_s_m_dyad_timefacfac,
        funct_m_m_dyad_timefacfac,
        normal
        );

    //-----------------------------------------------------------------
    // viscous consistency term

    // funct_m * 2 * mu_m * kappa_m * timefac * fac
    LINALG::Matrix<nen_,1> funct_m_viscm_timefacfac(funct_m);
    funct_m_viscm_timefacfac.Scale(viscm_fac);

    // funct_s * 2 * mu_m * kappa_m * timefac * fac
    LINALG::Matrix<slave_nen_,1> funct_s_viscm_timefacfac(funct_s);
    funct_s_viscm_timefacfac.Scale(viscm_fac);

    NIT_visc_Consistency_MasterTerms(
        C_umum,
        rhs_Cum,
        derxy_m,
        vderxy_m,
        funct_m_viscm_timefacfac,
        funct_s_viscm_timefacfac,
        normal);

    //-----------------------------------------------------------------
    // viscous adjoint consistency term

    LINALG::Matrix<nsd_,nen_> derxy_m_viscm_timefacfac(derxy_m);
    derxy_m_viscm_timefacfac.Scale(adj_visc_scale_*viscm_fac);
    NIT_visc_AdjointConsistency_MasterTerms(
        C_umum,
        rhs_Cum,
        velint_m,
        velint_s,
        derxy_m_viscm_timefacfac,
        funct_m,
        funct_s,
        normal );

    // we are done here!
    return;
  }

  //-----------------------------------------------------------------

  // evaluate the terms, that contribute to the background fluid
  // system - two-sided or xfluid-sided:

  // funct_m * timefac * fac * kappa_m
  LINALG::Matrix<nen_,1> funct_m_timefacfac_km(funct_m_timefacfac);
  funct_m_timefacfac_km.Scale(kappa_m);

  // funct_s * timefac * fac * kappa_m
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_km(funct_s_timefacfac);
  funct_s_timefacfac_km.Scale(kappa_m);

  // funct_m * timefac * fac * funct_m  * kappa_m (dyadic product)
  LINALG::Matrix<nen_,nen_> funct_m_m_dyad_timefacfac_km(funct_m_m_dyad_timefacfac);
  funct_m_m_dyad_timefacfac_km.Scale(kappa_m);

  // funct_s * timefac * fac * funct_m * kappa_m
  LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac_km(funct_s_m_dyad_timefacfac);
  funct_s_m_dyad_timefacfac_km.Scale(kappa_m);

  // 2 * mu_m * kappa_m * timefac * fac
  const double km_viscm_fac = viscm_fac * kappa_m;

  // funct_m * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<nen_,1> funct_m_viscm_timefacfac_km(funct_m);
  funct_m_viscm_timefacfac_km.Scale(km_viscm_fac);

  // funct_s * 2 * mu_m * kappa_m * timefac * fac
  LINALG::Matrix<slave_nen_,1> funct_s_viscm_timefacfac_km(funct_s);
  funct_s_viscm_timefacfac_km.Scale(km_viscm_fac);

  if (! full_slavesided)
  {
    //-----------------------------------------------------------------
    // pressure consistency term
    NIT_p_Consistency_MasterTerms(
        C_umum,
        rhs_Cum,
        press_m,
        funct_m_timefacfac_km,
        funct_s_timefacfac_km,
        funct_s_m_dyad_timefacfac_km,
        funct_m_m_dyad_timefacfac_km,
        normal
        );

    //-----------------------------------------------------------------
    // pressure adjoint consistency term

    NIT_p_AdjointConsistency_MasterTerms(
        C_umum,
        rhs_Cum,
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
        C_umum,
        rhs_Cum,
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
        C_umum,
        rhs_Cum,
        velint_m,
        velint_s,
        derxy_m_viscm_timefacfac_km,
        funct_m,
        funct_s,
        normal
        );
  }

  // in case of a purely xfluid-sided evaluation, we are done
  if (full_mastersided) return;

  //-----------------------------------------------------------------
  // the following quantities are only required for two-sided coupling
  // kappa_s > 0.0

  // 2 * mu_s * kappa_s * timefac * fac
  const double ks_viscs_fac = 2.0 * kappa_s * visceff_s * timefacfac;

  // funct_m * timefac * fac * kappa_s
  LINALG::Matrix<nen_,1> funct_m_timefacfac_ks(funct_m_timefacfac);
  funct_m_timefacfac_ks.Scale(kappa_s);

  // funct_s * timefac * fac * funct_m * kappa_s
  LINALG::Matrix<slave_nen_,nen_> funct_s_m_dyad_timefacfac_ks(funct_s_m_dyad_timefacfac);
  funct_s_m_dyad_timefacfac_ks.Scale(kappa_s);

  // funct_s * timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_timefacfac_ks(funct_s_timefacfac);
  funct_s_timefacfac_ks.Scale(kappa_s);

  // funct_s * timefac * fac * funct_s * kappa_s
  LINALG::Matrix<slave_nen_,slave_nen_> funct_s_s_dyad_timefacfac_ks(funct_s_s_dyad_timefacfac);
  funct_s_s_dyad_timefacfac_ks.Scale(kappa_s);

  // funct_m * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<nen_,1> funct_m_viscs_timefacfac_ks(funct_m);
  funct_m_viscs_timefacfac_ks.Scale(ks_viscs_fac);

  // funct_s * 2* mu_s  timefac * fac * kappa_s
  LINALG::Matrix<slave_nen_,1> funct_s_viscs_timefacfac_ks(funct_s);
  funct_s_viscs_timefacfac_ks.Scale(ks_viscs_fac);

  //-----------------------------------------------------------------
  // pressure consistency term

  double press_s = 0.0;
  // must use this-pointer because of two-stage lookup!
  this->GetInterfacePres(press_s);

  NIT_p_Consistency_SlaveTerms(
      rhs_Cum,
      press_s,
      funct_m_timefacfac_ks,
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
  this->GetInterfaceVelGrad(vderxy_s);

  NIT_visc_Consistency_SlaveTerms(
      rhs_Cum,
      derxy_s,
      vderxy_s,
      funct_m_viscs_timefacfac_ks,
      funct_s_viscs_timefacfac_ks,
      normal);

  //-----------------------------------------------------------------
  // viscous adjoint consistency term

  LINALG::Matrix<nsd_,slave_nen_> derxy_s_viscs_timefacfac_ks(derxy_s);
  derxy_s_viscs_timefacfac_ks.Scale(adj_visc_scale_*ks_viscs_fac);

  NIT_visc_AdjointConsistency_SlaveTerms(
      velint_m,
      velint_s,
      funct_m,
      funct_s,
      derxy_s_viscs_timefacfac_ks,
      normal
  );

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_applyAdditionalInterfacePenaltyTerms(
  Epetra_SerialDenseMatrix &          C_umum,                        ///< standard master-master couplingmatrix
  Epetra_SerialDenseVector &          rhs_Cum,                      ///< standard master coupling rhs
  const LINALG::Matrix<nsd_,1> &      normal,                       ///< normal vector
  const double &                      timefacfac,                   ///< theta*dt
  LINALG::Matrix<nen_,1> &            funct_m,                      ///< master shape function
  LINALG::Matrix<nsd_,nen_> &         derxy_m,                      ///< master shape function derivatives
  LINALG::Matrix<nsd_,nsd_> &         vderxy_m,                     ///< master spatial velocity derivatives
  const double &                      press_m,                      ///< master pressure
  const bool &                        velgrad_interface_stab,       ///< penalty term for velocity gradients at the interface
  const double &                      velgrad_interface_fac,        ///< stabilization fac for velocity gradients at the interface
  const bool &                        presscoupling_interface_stab, ///< penalty term for pressure coupling at the interface
  const double &                      presscoupling_interface_fac   ///< stabilization fac for pressure coupling at the interface
)
{
  //-----------------------------------------------------------------
  // penalty term for velocity gradients at the interface
  if (velgrad_interface_stab)
  {
    const double tau_timefacfac = velgrad_interface_fac * timefacfac;

    // get shape function gradient in normal direction
    LINALG::Matrix<nsd_,slave_nen_> deriv_s(true);
    this->GetSlaveFunctDeriv(deriv_s);
    LINALG::Matrix<slave_nen_,1> deriv_s_normal(true);
    deriv_s_normal.MultiplyTN(deriv_s,normal);

    LINALG::Matrix<nen_,1> deriv_m_normal(true);
    deriv_m_normal.MultiplyTN(derxy_m,normal);

    LINALG::Matrix<nsd_,nsd_> vderxy_s(true);
    this->GetInterfaceVelGrad(vderxy_s);
    LINALG::Matrix<nsd_,1> vderxy_s_normal(true);
    vderxy_s_normal.Multiply(vderxy_s,normal);

    LINALG::Matrix<nsd_,1> vderxy_m_normal(true);
    vderxy_m_normal.Multiply(vderxy_m,normal);

    // additional stability of gradients
    for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
    {
      for (unsigned ic=0; ic<slave_nen_; ++ic)
      {
        for (unsigned ir=0; ir<slave_nen_; ++ir)
        {
          const double tmp = tau_timefacfac*deriv_s_normal(ic)*deriv_s_normal(ir);
          C_usus_(sIndex(ir,ivel), sIndex(ic,ivel)) += tmp;
        }

        // neighbor row
        for (unsigned ir=0; ir<nen_; ++ir)
        {
          const double tmp = tau_timefacfac*deriv_s_normal(ic)*deriv_m_normal(ir);
          C_umus_(mIndex(ir,ivel), sIndex(ic,ivel)) -= tmp;
          C_usum_(sIndex(ic,ivel), mIndex(ir,ivel)) -= tmp;
        }
      }

      for (unsigned ir=0; ir<nen_; ++ir)
      {
        for (unsigned ic=0; ic<nen_; ++ic)
        {
          const double tmp = tau_timefacfac*deriv_m_normal(ir)*deriv_m_normal(ic);
          C_umum(mIndex(ir,ivel), mIndex(ic,ivel)) += tmp;
        }
      }

      // v_parent (u_neighbor-u_parent)
      for (unsigned ir=0; ir<slave_nen_; ++ir)
      {
        rhC_us_(sIndex(ir,ivel),0) += tau_timefacfac * deriv_s_normal(ir) * (vderxy_m_normal(ivel)-vderxy_s_normal(ivel));
      }

      // v_neighbor (u_neighbor-u_parent)
      for (unsigned ir=0; ir<nen_; ++ir)
      {
        rhs_Cum(mIndex(ir,ivel),0) +=  tau_timefacfac * deriv_m_normal(ir) * (vderxy_s_normal(ivel)-vderxy_m_normal(ivel));
      }
    }
  }// velgrad_interface_fac
  // --------------------------------------------------------
  // pressure coupling penalty term at the interface

  if (presscoupling_interface_stab)
  {
    double press_s = 0.0;
    this->GetInterfacePres(press_s);

    LINALG::Matrix<slave_nen_,1> funct_s;
    this->GetSlaveFunct(funct_s);

    const double tau_timefacfac_press = presscoupling_interface_fac * timefacfac;

    const double press_jump = press_m - press_s;

    for (unsigned ir=0; ir<slave_nen_; ++ir)
    {
      for (unsigned ic=0; ic<slave_nen_; ++ic)
      {
        const double tmp = tau_timefacfac_press*funct_s(ir)*funct_s(ic);
        C_usus_(sPres(ir),sPres(ic)) += tmp;

      }

      for (unsigned ic=0; ic<nen_; ++ic)
      {
        const double tmp = tau_timefacfac_press*funct_s(ir)*funct_m(ic);
        C_umus_(mPres(ic),sPres(ir)) -= tmp;

        C_usum_(sPres(ir),mPres(ic)) -= tmp;
      }

      rhC_us_(sPres(ir),0) += funct_s(ir)*press_jump*tau_timefacfac_press;
    }

    for (unsigned ir=0; ir<nen_; ++ir)
    {
      for (unsigned ic=0; ic<nen_; ++ic)
      {
        const double tmp = tau_timefacfac_press*funct_m(ir)*funct_m(ic);
        C_umum(mPres(ir),mPres(ic)) += tmp;
      }

      const double tmp = funct_m(ir)*press_jump*tau_timefacfac_press;
      rhs_Cum(mPres(ir),0) -= tmp;
    }
  } // presscoupling_interface_stab
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_MasterTerms(
  Epetra_SerialDenseMatrix &                C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &                rhs_Cum,                      ///< standard master-rhs
  const double &                            press_m,                      ///< master pressure
  const LINALG::Matrix<nen_,1> &            funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &      funct_s_timefacfac_km,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,nen_> &   funct_s_m_dyad_timefacfac_km, ///< (funct_s^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nen_,nen_> &         funct_m_m_dyad_timefacfac_km, ///< (funct^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &            normal                        ///< normal vector
)
{
             /*                  \       /          i      \
          + |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
             \                   /       \                */

  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    //-----------------------------------------------
    //    + (vm, km *(Dpm)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      const unsigned row = mIndex(ir,ivel);
      for (unsigned ic =0; ic<nen_; ic++)
      {
        // (v,Dp*n)
        C_umum(row, mPres(ic)) += funct_m_m_dyad_timefacfac_km(ir,ic)*normal(ivel);
      }

      // -(v,p*n)
      const double funct_m_km_timefacfac_press = funct_m_timefacfac_km(ir)*press_m;
      rhs_Cum(row,0) -= funct_m_km_timefacfac_press*normal(ivel);
    }

    if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC) continue;

    //-----------------------------------------------
    //    - (vs, km *(Dpm)*n)
    //-----------------------------------------------
    for(unsigned ir = 0; ir<slave_nen_; ir++)
    {
      // (v,Dp*n)
      for(unsigned ic =0; ic<nen_; ic++)
      {
        C_usum_(sIndex(ir,ivel),mPres(ic)) -= funct_s_m_dyad_timefacfac_km(ir,ic)*normal(ivel);
      }

      // -(v,p*n)
      const double funct_s_km_timefacfac_press = funct_s_timefacfac_km(ir)*press_m;
      rhC_us_(sIndex(ir,ivel),0) += funct_s_km_timefacfac_press*normal(ivel);
    }
  } // end loop over velocity components
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_Consistency_SlaveTerms(
  Epetra_SerialDenseVector &                    rhs_Cum,                      ///< standard master-rhs
  const double &                                pres_s,                       ///< slave pressure
  const LINALG::Matrix<nen_,1> &                funct_m_timefacfac_ks,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,1> &          funct_s_timefacfac_ks,        ///< funct_s * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,nen_> &       funct_s_m_dyad_timefacfac_ks, ///< (funct_s^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac_ks, ///< (funct^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &                normal                        ///< normal vector
)
{
  // loop over velocity components
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    //-----------------------------------------------
    //    + (vm, ks *(Dps)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<nen_; ir++)
    {
      for (unsigned ic =0; ic<slave_nen_; ic++)
      {
        // (vm, ks * Dps*n)
        C_umus_(mIndex(ir,ivel),sPres(ic)) += funct_s_m_dyad_timefacfac_ks(ic,ir)*normal(ivel);
      }

      // -(vm, ks * ps*n)
      const double funct_m_timefacfac_ks_press = funct_m_timefacfac_ks(ir)*pres_s;
      rhs_Cum(mIndex(ir,ivel),0) -= funct_m_timefacfac_ks_press*normal(ivel);
    }

    //-----------------------------------------------
    //    - (vs, ks *(Dps)*n)
    //-----------------------------------------------
    for (unsigned ir = 0; ir<slave_nen_; ir++)
    {
      // -(vs,ks *Dps*n)
      for (unsigned ic =0; ic<slave_nen_; ic++)
      {
        C_usus_(sIndex(ir,ivel),sPres(ic)) -= funct_s_s_dyad_timefacfac_ks(ir,ic)*normal(ivel);
      }

      // +(vs,ks * ps*n)
      const double funct_s_timefacfac_ks_press = funct_s_timefacfac_ks(ir)*pres_s;
      rhC_us_(sIndex(ir,ivel),0) += funct_s_timefacfac_ks_press*normal(ivel);
    }
  } // end loop over velocity components
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_p_AdjointConsistency_MasterTerms(
  Epetra_SerialDenseMatrix &              C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &              rhs_Cum,                      ///< standard master-rhs
  const double &                          velint_normal_m,              ///< velocity in normal direction
  const double &                          velint_normal_s,              ///< interface velocity in normal direction
  const LINALG::Matrix<nen_,1> &          funct_m_timefacfac_km,        ///< funct * timefacfac *kappa_m
  const LINALG::Matrix<slave_nen_,nen_> & funct_s_m_dyad_timefacfac_km, ///< (funct^T * funct_s) * timefacfac *kappa_m
  const LINALG::Matrix<nen_,nen_> &       funct_m_m_dyad_timefacfac_km, ///< (funct^T * funct) * timefacfac *kappa_m
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
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
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned ic =0; ic<nen_; ic++)
    {
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // - (qm*n, km *(Dum))
        C_umum(mPres(ir),mIndex(ic,ivel)) -= funct_m_m_dyad_timefacfac_km(ir,ic)*normal(ivel);
      }
    }

    // (qm*n, km * um)
    rhs_Cum(mPres(ir),0) += funct_m_timefacfac_km(ir)*velint_normal_m;

    // -(qm*n,km * u_DBC) for weak DBC or
    // -(qm*n,km * us)
    rhs_Cum(mPres(ir),0) -= funct_m_timefacfac_km(ir)*velint_normal_s;
  }

  if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC) return;

  //-----------------------------------------------
  //    + (qm*n, km *(Dus))
  //-----------------------------------------------
  for (unsigned ir = 0; ir<nen_; ir++)
  {
    for (unsigned ic =0; ic<slave_nen_; ic++)
    {
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // -(qm*n, km * Dus)
        C_umus_(mPres(ir),sIndex(ic,ivel)) += funct_s_m_dyad_timefacfac_km(ic,ir)*normal(ivel);
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
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    for (unsigned ic =0; ic<nen_; ic++)
    {
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // -(qs*n, ks* Dum)
        C_usum_(sPres(ir),mIndex(ic,ivel)) -= funct_s_m_dyad_timefacfac_ks(ir,ic)*normal(ivel);
      }
    }

    // (qs*n,ks* um)
    rhC_us_(sPres(ir),0) += funct_s_timefacfac_ks(ir)*velint_normal_m;
  }


  //-----------------------------------------------
  //    + (qs*n, ks *(Dus))
  //-----------------------------------------------
  for (unsigned ir = 0; ir<slave_nen_; ir++)
  {
    for (unsigned ic =0; ic<slave_nen_; ic++)
    {
      for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
      {
        // +(qs*n, ks* Dus)
        C_usus_(sPres(ir),sIndex(ic,ivel)) += funct_s_s_dyad_timefacfac_ks(ir,ic)*normal(ivel);
      }
    }

    // -(qs*n,ks *us)
    rhC_us_(sPres(ir),0) -= funct_s_timefacfac_ks(ir)*velint_normal_s;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_MasterTerms(
  Epetra_SerialDenseMatrix &            C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &            rhs_Cum,                      ///< standard master-rhs
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

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {

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
        const unsigned row = mIndex(ir,ivel);
        for (unsigned ic =0; ic<nen_; ic++)
        {
          // diagonal block
          C_umum(row, mIndex(ic,ivel)) -= funct_m_viscm_timefacfac_km(ir) * 0.5 * normal(jvel)*derxy_m(jvel,ic);

          // off diagonal block
          C_umum(row, mIndex(ic,jvel)) -= funct_m_viscm_timefacfac_km(ir) * 0.5 * normal(jvel)*derxy_m(ivel,ic);
        }

        // rhs
        rhs_Cum(row,0) += 0.5 * funct_m_viscm_timefacfac_km(ir) * ( vderxy_m(ivel,jvel) + vderxy_m(jvel,ivel) )*normal(jvel);
      }

      if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC) continue;

      //-----------------------------------------------
      //    + (vs, (2*km*mum) *eps(Dum)*n)
      //-----------------------------------------------
      for (unsigned ir = 0; ir<slave_nen_; ir++)
      {
        const unsigned row = sIndex(ir,ivel);
        for (unsigned ic =0; ic<nen_; ic++)
        {
          // diagonal block
          C_usum_(row,mIndex(ic,ivel)) += funct_s_viscm_timefacfac_km(ir) * 0.5 * normal(jvel)*derxy_m(jvel,ic);

          // off-diagonal block
          C_usum_(row,mIndex(ic,jvel)) += funct_s_viscm_timefacfac_km(ir) * 0.5 * normal(jvel)*derxy_m(ivel,ic);
        }

        // rhs
        rhC_us_(row,0) -= 0.5 * funct_s_viscm_timefacfac_km(ir) * ( vderxy_m(ivel,jvel) + vderxy_m(jvel,ivel) )*normal(jvel);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_Consistency_SlaveTerms(
  Epetra_SerialDenseVector &              rhs_Cum,                      ///< standard master-rhs
  const LINALG::Matrix<nsd_,slave_nen_> & derxy_s,                      ///< slave shape function derivatives
  const LINALG::Matrix<nsd_,nsd_> &       vderxy_s,                     ///< slave velocity gradient
  const LINALG::Matrix<nen_,1> &          funct_m_viscs_timefacfac_ks,  ///< funct_m*mu_m*timefacfac
  const LINALG::Matrix<slave_nen_,1> &    funct_s_viscs_timefacfac_ks,  ///< funct_s*mu_m*timefacfac
  const LINALG::Matrix<nsd_,1> &          normal                        ///< normal vector
)
{
  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
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
        const unsigned row = mIndex(ir,ivel);
        for (unsigned ic =0; ic<slave_nen_; ic++)
        {
          // diagonal block
          C_umus_(row, sIndex(ic,ivel)) -= funct_m_viscs_timefacfac_ks(ir)* 0.5 * normal(jvel)*derxy_s(jvel,ic);

          // off diagonal block
          C_umus_(row, sIndex(ic,jvel)) -= funct_m_viscs_timefacfac_ks(ir)* 0.5 * normal(jvel)*derxy_s(ivel,ic);
        }

        // rhs
        rhs_Cum(row,0) += 0.5 * funct_m_viscs_timefacfac_ks(ir) * ( vderxy_s(ivel,jvel) + vderxy_s(jvel,ivel) )*normal(jvel);
      }

      //-----------------------------------------------
      //    + (vs, (2*ks*mus) *eps(Dus)*n)
      //-----------------------------------------------
      for (unsigned ir = 0; ir<slave_nen_; ir++)
      {
        const unsigned row = sIndex(ir,ivel);
        for (unsigned ic =0; ic<slave_nen_; ic++)
        {
          // diagonal block
          C_usus_(row,sIndex(ic,ivel)) += funct_s_viscs_timefacfac_ks(ir) * 0.5 * normal(jvel)*derxy_s(jvel,ic);

          // off diagonal block
          C_usus_(row,sIndex(ic,jvel)) += funct_s_viscs_timefacfac_ks(ir) * 0.5 * normal(jvel)*derxy_s(ivel,ic);
        }

        // rhs
        rhC_us_(row,0) -= 0.5 * funct_s_viscs_timefacfac_ks(ir) * ( vderxy_s(ivel,jvel) + vderxy_s(jvel,ivel) )*normal(jvel);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_visc_AdjointConsistency_MasterTerms(
  Epetra_SerialDenseMatrix &               C_umum,                            ///< standard master-master-matrix
  Epetra_SerialDenseVector &               rhs_Cum,                           ///< standard master-rhs
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

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
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
        const unsigned mrow = mIndex(ir,ivel);
        for (unsigned ic =0; ic<nen_; ic++)
        {
          const double tmp = funct_m(ic)* 0.5 * derxy_m_viscm_timefacfac_km(jvel,ir);
          // diagonal block
          C_umum(mrow, mIndex(ic,ivel)) -= normal(jvel) * tmp;

          // off diagonal block
          C_umum(mrow, mIndex(ic,jvel)) -= normal(ivel) * tmp;
        }

        // rhs
        const double tmp = derxy_m_viscm_timefacfac_km(jvel,ir) * 0.5;
        rhs_Cum(mrow,0) += tmp * (normal(jvel) * velint_m(ivel) + normal(ivel)*velint_m(jvel));
        rhs_Cum(mrow,0) -= tmp * (normal(jvel) * velint_s(ivel) + normal(ivel)*velint_s(jvel));
      }

      if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC) continue;

      //-----------------------------------------------
      //    + alpha * ((2*km*mum) *eps(vm)*n , us)
      //-----------------------------------------------
      for (unsigned ir = 0; ir<nen_; ir++)
      {
        const unsigned mrow = mIndex(ir,ivel);
        for (unsigned ic =0; ic<slave_nen_; ic++)
        {
          const double tmp = funct_s(ic)* 0.5 * derxy_m_viscm_timefacfac_km(jvel,ir);
          // diagonal block
          C_umus_(mrow,sIndex(ic,ivel)) += normal(jvel) * tmp;

          // off diagonal block
          C_umus_(mrow,sIndex(ic,jvel)) += normal(ivel) * tmp;
        }
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

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
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

      //-----------------------------------------------
      //    - alpha * ((2*ks*mus) *eps(vs)*n , um)
      //-----------------------------------------------
      for (unsigned ir = 0; ir<nen_; ir++)
      {
        const unsigned row = sIndex(ir,ivel);
        for (unsigned ic =0; ic<nen_; ic++)
        {
          // diagonal block
          C_usum_(row, mIndex(ic,ivel)) -= funct_m(ic)* 0.5 * normal(jvel)* derxy_s_viscs_timefacfac_ks(jvel,ir);

          // off diagonal block
          C_usum_(row, mIndex(ic,jvel)) -= funct_m(ic)* 0.5 * normal(ivel)* derxy_s_viscs_timefacfac_ks(jvel,ir);
        }

        // rhs
        rhC_us_(row,0) += derxy_s_viscs_timefacfac_ks(jvel,ir) * 0.5* (normal(jvel) * velint_m(ivel) + normal(ivel)*velint_m(jvel));
      }

      //-----------------------------------------------
      //    + alpha * ((2*km*mum) *eps(vm)*n , us)
      //-----------------------------------------------
      for (unsigned ir = 0; ir<slave_nen_; ir++)
      {
        const unsigned row = sIndex(ir,ivel);
        for (unsigned ic =0; ic<slave_nen_; ic++)
        {
          // diagonal block
          C_usus_(row,sIndex(ic,ivel)) += funct_s(ic)* 0.5 * normal(jvel)* derxy_s_viscs_timefacfac_ks(jvel,ir);

          // off diagonal block
          C_usus_(row,sIndex(ic,jvel)) += funct_s(ic)* 0.5 * normal(ivel)* derxy_s_viscs_timefacfac_ks(jvel,ir);
        }

        // rhs
        rhC_us_(row,0) -= derxy_s_viscs_timefacfac_ks(jvel,ir) * 0.5* (normal(jvel) * velint_s(ivel) + normal(ivel)*velint_s(jvel));
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty(
  Epetra_SerialDenseMatrix &                    C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &                    rhs_Cum,                      ///< standard master-rhs
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,slave_nen_>&  funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
  const double &                                stabfac                      ///< stabilization factor
)
{
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

  for (unsigned ivel = 0; ivel < nsd_; ivel ++)
  {
    // + gamma*mu/h_K (vm, um))
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const unsigned m = mIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        const double tmp = funct_m_m_dyad_timefacfac(ir,ic)*stabfac;

        C_umum(m,mIndex(ic,ivel)) += tmp;
      }

      const double tmp = funct_m_timefacfac(ir)*stabfac;

      rhs_Cum(m,0) -= tmp*velint_m(ivel);

      // +(stab * vm, u_DBC) (weak dirichlet case) or from
      // +(stab * vm, u_s)
      rhs_Cum(m,0) += tmp*velint_s(ivel);
    }

    if (applicationType_ == SlaveElementInterface<distype>::XFluidWDBC) continue;

    // - gamma*mu/h_K (vm, us))
    // - gamma*mu/h_K (vs, um))
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const unsigned s = sIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        const double tmp = funct_s_m_dyad_timefacfac(ir,ic)*stabfac;

        C_usum_(s,mIndex(ic,ivel)) -= tmp;

        C_umus_(mIndex(ic,ivel),s) -= tmp;
      }

      const double tmp = funct_s_timefacfac(ir)*stabfac;

      // +(stab * vs, um)
      rhC_us_(s,0) += tmp*velint_m(ivel);

      // + gamma*mu/h_K (vs, us))
      for (unsigned ic=0; ic<slave_nen_; ic++)
      {
        C_usus_(s,sIndex(ic,ivel)) += funct_s_s_dyad_timefacfac(ir,ic)*stabfac;
      }

      // -(stab * vs, us)
      rhC_us_(s,0) -= tmp*velint_s(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_Penalty_MasterTerms(
  Epetra_SerialDenseMatrix &                    C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &                    rhs_Cum,                      ///< standard master-rhs
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const double &                                stabfac
)
{
  for (unsigned ivel = 0; ivel < nsd_; ivel ++)
  {
    // + gamma*mu/h_K (vm, um))
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const unsigned mrow = mIndex(ir,ivel);

      for (unsigned ic=0; ic<nen_; ic++)
      {
        C_umum(mrow,mIndex(ic,ivel)) += funct_m_m_dyad_timefacfac(ir,ic)*stabfac;
      }

      const double tmp = funct_m_timefacfac(ir)*stabfac;

      rhs_Cum(mrow,0) -= tmp*velint_m(ivel);

      // +(stab * vm, u_DBC) (weak dirichlet case) or from
      // +(stab * vm, u_s)
      rhs_Cum(mrow,0) += tmp*velint_s(ivel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType slave_distype, unsigned int slave_numdof>
void NitscheCoupling<distype,slave_distype,slave_numdof>::NIT_Stab_ConvAveraged(
  Epetra_SerialDenseMatrix &                    C_umum,                       ///< standard master-master-matrix
  Epetra_SerialDenseVector &                    rhs_Cum,                      ///< master rhs
  const LINALG::Matrix<nsd_,1>&                 velint_m,                     ///< velocity at integration point
  const LINALG::Matrix<nsd_,1>&                 velint_s,                     ///< interface velocity at integration point
  const LINALG::Matrix<nen_,1>&                 funct_m_timefacfac,           ///< funct * timefacfac
  const LINALG::Matrix<nen_,nen_>&              funct_m_m_dyad_timefacfac,    ///< (funct^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,1>&           funct_s_timefacfac,           ///< funct_s^T * timefacfac
  const LINALG::Matrix<slave_nen_,nen_>&        funct_s_m_dyad_timefacfac,    ///< (funct_s^T * funct) * timefacfac
  const LINALG::Matrix<slave_nen_,slave_nen_> & funct_s_s_dyad_timefacfac,    ///< (funct_s^T * funct_s) * timefacfac
  const double &                                stabfac_avg                  ///< stabilization factor
)
{
  /*                                        \        /                                       i \
 -|  [rho * (beta * n)] *  { v }_m , [ Du ] | =  + |  [rho * (beta * n)] * { v }_m,  [ u ]      |
  \ ----stab_avg-----                       /         \ ----stab_avg-----                      */

  for (unsigned ivel = 0; ivel<nsd_; ivel ++)
  {
    //  [rho * (beta * n^b)] (0.5*vb,ub)
    for (unsigned ir=0; ir<nen_; ir++)
    {
      const unsigned mrow = mIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        C_umum(mrow,mIndex(ic,ivel)) -= 0.5*funct_m_m_dyad_timefacfac(ir,ic)*stabfac_avg;
      }

      rhs_Cum(mrow,0) += 0.5*funct_m_timefacfac(ir)*stabfac_avg*velint_m(ivel);

    //  -[rho * (beta * n^b)] (0.5*vb,ue)
      for (unsigned ic=0; ic<slave_nen_; ic++)
      {
        C_umus_(mrow,sIndex(ic,ivel)) += 0.5*funct_s_m_dyad_timefacfac(ic,ir)*stabfac_avg;
      }

      rhs_Cum(mrow,0) -= 0.5*funct_m_timefacfac(ir)*stabfac_avg*velint_s(ivel);
    }

    //  [rho * (beta * n^b)] (0.5*ve,ub)
    for (unsigned ir=0; ir<slave_nen_; ir++)
    {
      const unsigned srow = sIndex(ir,ivel);
      for (unsigned ic=0; ic<nen_; ic++)
      {
        C_usum_(srow,mIndex(ic,ivel)) -= 0.5*funct_s_m_dyad_timefacfac(ir,ic)*stabfac_avg;
      }
      rhC_us_(srow,0) += 0.5*funct_s_timefacfac(ir)*stabfac_avg*velint_m(ivel);

    //-[rho * (beta * n^b)] (0.5*ve,ue)

      for (unsigned ic=0; ic<slave_nen_; ic++)
      {
        C_usus_(srow,sIndex(ic,ivel)) += 0.5*funct_s_s_dyad_timefacfac(ir,ic)*stabfac_avg;
      }

      rhC_us_(srow,0) -= 0.5*funct_s_timefacfac(ir)*stabfac_avg*velint_s(ivel);
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
  LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s       ///< block rhs vector \f$ rhs_{\sigma} \f$
)
{
  LINALG::Matrix<nen_,slave_nen_> bK_ms;
  LINALG::Matrix<slave_nen_,nen_> bK_sm;

  // interface velocity at gauss-point (current gauss-point in calling method)
  LINALG::Matrix<nsd_,1> velint_s(true);
  this->GetInterfaceVel(velint_s);

  // get nodal shape function vector
  LINALG::Matrix<slave_nen_,1> slave_funct(true);
  this->GetSlaveFunct(slave_funct);

  bK_ms.MultiplyNT(funct,slave_funct);
  bK_sm.UpdateT(bK_ms);

  for (unsigned ivel = 0; ivel < nsd_; ++ ivel)
  {
    for (unsigned jvel = 0; jvel < nsd_; ++ jvel)
    {
      // G_sus
      BG_sus_( stressIndex(ivel,jvel), ivel )->Update( fac*normal(jvel), bK_ms, 1.0 );
      rhs_s( stressIndex(ivel,jvel), 0 )->Update( -fac*normal(jvel)*velint_s(ivel), funct, 1.0 );

      // G_uss
      BG_uss_( ivel, stressIndex(ivel,jvel) )->Update( fac*normal(jvel), bK_sm, 1.0 );
    }
  }

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
  LINALG::Matrix<nen_,1> &                                        rhs_pmus    ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
)
{
  // interface velocity at gauss-point (current gauss-point in calling method)
  LINALG::Matrix<nsd_,1> velint_s(true);
  this->GetInterfaceVel(velint_s);

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

  return;
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
