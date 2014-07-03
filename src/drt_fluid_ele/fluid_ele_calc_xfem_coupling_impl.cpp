/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling_impl.cpp

\brief Implementation class for coupling of two different meshes using Stress/Hybrid method of Nitsche's method to enforce
       interface conditions weakly
<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

#include "../drt_lib/drt_utils.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"

#include "../linalg/linalg_utils.H"

#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"


namespace DRT
{
namespace ELEMENTS
{
namespace XFLUID
{

/*--------------------------------------------------------------------------------
 * add side's interface displacements and set current side node coordinates
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::addeidisp(
    const DRT::Discretization &  cutdis,       ///< cut discretization
    const std::string            state,        ///< state
    const std::vector<int>&      lm            ///< local map
    )
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vector
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<side_nen_; ++inode)  // number of nodes
  {
    for(int idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (eidisp_)(idim,inode) = mymatrix[idim+(inode*side_numdof)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
    }
  }

  // add the displacement of the interface
  for (int inode = 0; inode < side_nen_; ++inode)
  {
    xyze_(0,inode) += eidisp_(0, inode);
    xyze_(1,inode) += eidisp_(1, inode);
    xyze_(2,inode) += eidisp_(2, inode);
  }

  return;
} // addeidisp


/*--------------------------------------------------------------------------------
 * extract/set side's interface velocity
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::eivel(
    const DRT::Discretization &  cutdis,  ///< cut discretization
    const std::string            state,   ///< state
    const std::vector<int>&      lm       ///< local map
    )
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<side_nen_; ++inode)  // number of nodes
  {
    for(int idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      (eivel_)(idim,inode) = mymatrix[idim+(inode*side_numdof)];  // state vector includes velocity and pressure
    }
    if (side_numdof == nsd_+1) (eipres_)(inode,0) = mymatrix[nsd_+(inode*side_numdof)];
  }

  return;
} // eVel

/*--------------------------------------------------------------------------------
 * extract/set side's interface velocity
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::getivelint(
    LINALG::Matrix<nsd_,1>& ivelint  ///< interface velocity at embedded side
    )
{
  // interface velocity vector in gausspoint
  ivelint.Multiply(eivel_,side_funct_);

  return;
}

/*--------------------------------------------------------------------------------
 * evaluate shape function, derivatives, normal and transformation w.r.t
 * side element at gaussian point
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::Evaluate(
    const LINALG::Matrix<nsd_-1,1>  & eta,     ///< local coordinates w.r.t side element
    LINALG::Matrix<nsd_,1>          & x,       ///< global coordinates of gaussian point
    LINALG::Matrix<nsd_,1>          & normal,  ///< normal vector
    double                          & drs      ///< transformation factor
)
{

  LINALG::Matrix<nsd_-1,2> metrictensor;
  DRT::UTILS::shape_function_2D( side_funct_, eta( 0 ), eta( 1 ), side_distype );
  DRT::UTILS::shape_function_2D_deriv1( side_deriv_, eta( 0 ), eta( 1 ), side_distype );
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<side_distype>( xyze_,side_deriv_, metrictensor, drs, &normal );
  x.Multiply( xyze_,side_funct_ );

  return;
} // Evaluate


/*--------------------------------------------------------------------------------
 * compute interface force for side nodes
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::InterfaceForce(
    Epetra_SerialDenseVector &   iforce,     ///< interface force vector
    LINALG::Matrix<nsd_,1> &     traction,   ///< traction vector at gaussian point
    const double &               fac         ///< integration factor
    )
{

  if(side_numdof != nsd_) dserror(" pay attention in buildInterfaceForce: side_numdof != nsd_");

  for (int inode = 0; inode < side_nen_; ++inode)
  {
    for(int idim=0; idim<nsd_; ++idim )
    {
      // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
      iforce[idim+(inode*side_numdof)] += side_funct_(inode) * traction(idim) * fac;
    }
  }

  return;
}


/*--------------------------------------------------------------------------------
 * project gaussian point from linearized interfac in normal direction onto
 * corresponding side
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::ProjectOnSide(
    LINALG::Matrix<3,1> & x_gp_lin,  ///< global coordinates of gaussian point w.r.t linearized interface
    LINALG::Matrix<3,1> & x_side,    ///< projected gaussian point on side
    LINALG::Matrix<2,1> & xi_side    ///< local coordinates of projected gaussian point w.r.t side
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::ProjectOnSide" );
  // Initialization
  LINALG::Matrix<side_nen_,1> funct(true);      // shape functions
  LINALG::Matrix<2,side_nen_> deriv(true);      // derivatives dr, ds
  LINALG::Matrix<3,side_nen_> deriv2(true);     // 2nd derivatives drdr, dsds, drds


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

  if(side_distype == DRT::Element::tri3 or
     side_distype == DRT::Element::tri6)
  {
    sol(0) = 0.333333333333333;
    sol(1) = 0.333333333333333;
  }
  else if( side_distype == DRT::Element::quad4 or
           side_distype == DRT::Element::quad8 or
           side_distype == DRT::Element::quad9)
  {

    sol(0) = 0.0;
    sol(1) = 0.0;
  }
  else
  {
    dserror("define start side xi-coordinates for unsupported cell type");
  }


  // we only use absolute tolerances, since we compute local coordinates
  const double absTolIncr = 1.0e-9;   // abs tolerance for the local coordinates increment
  const double absTolRes  = 1.0e-9;   // abs tolerance for the whole residual
  const double absTOLdist = 1.0e-9;   // abs tolerance for distance

  int iter=0;
  const int maxiter = 7;

  bool converged = false;

  while(iter < maxiter && !converged)
  {

    iter++;


    // get current values
    DRT::UTILS::shape_function_2D( funct, sol( 0 ), sol( 1 ), side_distype );
    DRT::UTILS::shape_function_2D_deriv1( deriv, sol( 0 ), sol( 1 ), side_distype );
    DRT::UTILS::shape_function_2D_deriv2( deriv2, sol( 0 ), sol( 1 ), side_distype );

    x.Multiply(xyze_, funct);

    derxy.MultiplyNT(xyze_, deriv);

    derxy2.MultiplyNT(xyze_, deriv2);

    // set dx_dr and dx_ds
    for (int i=0; i< 3; i++)
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
    for(int i=0; i< 3; i++)
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
    std::cout << "side " << xyze_ << std::endl;

    dserror( "newton scheme in ProjectOnSide not converged! " );
  }



  // evaluate shape function at solution
  DRT::UTILS::shape_function_2D( side_funct_, sol( 0 ), sol( 1 ), side_distype );

  // get projected gauss point
  x_side.Multiply(xyze_, side_funct_);

  // set local coordinates w.r.t side
  xi_side(0) = sol(0);
  xi_side(1) = sol(1);

  return;
} // ProjectOnSide











/*--------------------------------------------------------------------------------
 * build side coupling matrices for mixed/hybrid Cauchy stress-based LM
 * (MHCS) coupling
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::MHCS_buildCouplingMatrices(
    const LINALG::Matrix<nsd_,1> &                                normal,     ///< normal vector
    const double &                                                fac,        ///< integration factor
    const LINALG::Matrix<nen_,1> &                                funct,      ///< shape function
    LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> & rhs_s       ///< rhs block matrix
)
{
  LINALG::Matrix<nen_,side_nen_> bKi_ss;
  LINALG::Matrix<side_nen_,nen_> bKiT_ss;

  // interface velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> ivelint;
  getivelint(ivelint);

  const unsigned Sigmaxx = 0;
  const unsigned Sigmaxy = 1;
  const unsigned Sigmaxz = 2;
  const unsigned Sigmayx = 1;
  const unsigned Sigmayy = 3;
  const unsigned Sigmayz = 4;
  const unsigned Sigmazx = 2;
  const unsigned Sigmazy = 4;
  const unsigned Sigmazz = 5;
  const unsigned Velxi = 0;
  const unsigned Velyi = 1;
  const unsigned Velzi = 2;

  for(int i=0;i<nen_;++i)
  {
    for(int j=0;j<side_nen_;++j)
    {
      bKi_ss(i,j) = funct(i)*side_funct_(j);
    }
  }

  // G_sui (coupling only for sigma and u (not for sigma and p)
  BG_sui_( Sigmaxx, Velxi )->Update( fac*normal(0), bKi_ss, 1.0 );
  BG_sui_( Sigmaxy, Velxi )->Update( fac*normal(1), bKi_ss, 1.0 );
  BG_sui_( Sigmaxz, Velxi )->Update( fac*normal(2), bKi_ss, 1.0 );
  BG_sui_( Sigmayx, Velyi )->Update( fac*normal(0), bKi_ss, 1.0 );
  BG_sui_( Sigmayy, Velyi )->Update( fac*normal(1), bKi_ss, 1.0 );
  BG_sui_( Sigmayz, Velyi )->Update( fac*normal(2), bKi_ss, 1.0 );
  BG_sui_( Sigmazx, Velzi )->Update( fac*normal(0), bKi_ss, 1.0 );
  BG_sui_( Sigmazy, Velzi )->Update( fac*normal(1), bKi_ss, 1.0 );
  BG_sui_( Sigmazz, Velzi )->Update( fac*normal(2), bKi_ss, 1.0 );


  rhs_s( Sigmaxx, 0 )->Update( -fac*normal(0)*ivelint(0), funct, 1.0 );
  rhs_s( Sigmaxy, 0 )->Update( -fac*normal(1)*ivelint(0), funct, 1.0 );
  rhs_s( Sigmaxz, 0 )->Update( -fac*normal(2)*ivelint(0), funct, 1.0 );
  rhs_s( Sigmayx, 0 )->Update( -fac*normal(0)*ivelint(1), funct, 1.0 );
  rhs_s( Sigmayy, 0 )->Update( -fac*normal(1)*ivelint(1), funct, 1.0 );
  rhs_s( Sigmayz, 0 )->Update( -fac*normal(2)*ivelint(1), funct, 1.0 );
  rhs_s( Sigmazx, 0 )->Update( -fac*normal(0)*ivelint(2), funct, 1.0 );
  rhs_s( Sigmazy, 0 )->Update( -fac*normal(1)*ivelint(2), funct, 1.0 );
  rhs_s( Sigmazz, 0 )->Update( -fac*normal(2)*ivelint(2), funct, 1.0 );

  bKiT_ss.UpdateT(bKi_ss);

  // G_uis
  BG_uis_( Velxi, Sigmaxx )->Update( fac*normal(0), bKiT_ss, 1.0 );
  BG_uis_( Velxi, Sigmaxy )->Update( fac*normal(1), bKiT_ss, 1.0 );
  BG_uis_( Velxi, Sigmaxz )->Update( fac*normal(2), bKiT_ss, 1.0 );
  BG_uis_( Velyi, Sigmayx )->Update( fac*normal(0), bKiT_ss, 1.0 );
  BG_uis_( Velyi, Sigmayy )->Update( fac*normal(1), bKiT_ss, 1.0 );
  BG_uis_( Velyi, Sigmayz )->Update( fac*normal(2), bKiT_ss, 1.0 );
  BG_uis_( Velzi, Sigmazx )->Update( fac*normal(0), bKiT_ss, 1.0 );
  BG_uis_( Velzi, Sigmazy )->Update( fac*normal(1), bKiT_ss, 1.0 );
  BG_uis_( Velzi, Sigmazz )->Update( fac*normal(2), bKiT_ss, 1.0 );

  return;
}// MHCS_buildCouplingMatrices

/*--------------------------------------------------------------------------------
 * build side coupling matrices for mixed/hybrid viscous stress-based LM
 * (MHCS) coupling
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::MHVS_buildCouplingMatrices(
    const LINALG::Matrix<nsd_,1> &                                    normal,     ///< normal vector
    const double &                                                    fac,        ///< integration factor
    const LINALG::Matrix<nen_,1> &                                    funct,      ///< background element shape functions
    LINALG::BlockMatrix<LINALG::Matrix<nen_,1>,numstressdof_,1> &     rhs_s,      ///< block rhs vector \f$ rhs_{\sigma}\f$
    const double &                                                    press,      ///< background element pressure
    LINALG::Matrix<nen_,1> &                                          rhs_pui     ///< part of block rhs vector \f$rhs_p\f$ including interface velocity terms
)
{
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;
  const unsigned Pres = 3;

  const unsigned Sigmaxx = 0;
  const unsigned Sigmaxy = 1;
  const unsigned Sigmaxz = 2;
  const unsigned Sigmayx = 1;
  const unsigned Sigmayy = 3;
  const unsigned Sigmayz = 4;
  const unsigned Sigmazx = 2;
  const unsigned Sigmazy = 4;
  const unsigned Sigmazz = 5;

  // interface velocity at gauss-point (current gauss-point in calling method)
  LINALG::Matrix<nsd_,1> ivelint(true);
  getivelint(ivelint);

  // block submatrices for side-background coupling; i: interface, b: background
  LINALG::Matrix<nen_, side_nen_> bG_bi(true);
  LINALG::Matrix<side_nen_, nen_> bG_ib(true);

  bG_bi.MultiplyNT(funct, side_funct_);
  bG_ib.MultiplyNT(side_funct_, funct);

  /*
   * G_sui
   *
   * from:
   *   /            \
   *  | \tau n, u^i  |
   *   \            /
   */
  BG_sui_(Sigmaxx, Velx)->Update(fac * normal(Velx), bG_bi, 1.0);
  BG_sui_(Sigmaxy, Velx)->Update(fac * normal(Vely), bG_bi, 1.0);
  BG_sui_(Sigmaxz, Velx)->Update(fac * normal(Velz), bG_bi, 1.0);
  BG_sui_(Sigmayx, Vely)->Update(fac * normal(Velx), bG_bi, 1.0);
  BG_sui_(Sigmayy, Vely)->Update(fac * normal(Vely), bG_bi, 1.0);
  BG_sui_(Sigmayz, Vely)->Update(fac * normal(Velz), bG_bi, 1.0);
  BG_sui_(Sigmazx, Velz)->Update(fac * normal(Velx), bG_bi, 1.0);
  BG_sui_(Sigmazy, Velz)->Update(fac * normal(Vely), bG_bi, 1.0);
  BG_sui_(Sigmazz, Velz)->Update(fac * normal(Velz), bG_bi, 1.0);

  /*
   *  G_uis
   *
   *  from:
   *   /               \
   *  | v^i, \sigma n   |
   *   \               /
   *
   */
  BG_uis_(Velx, Sigmaxx)->Update(fac * normal(Velx), bG_ib, 1.0);
  BG_uis_(Velx, Sigmaxy)->Update(fac * normal(Vely), bG_ib, 1.0);
  BG_uis_(Velx, Sigmaxz)->Update(fac * normal(Velz), bG_ib, 1.0);
  BG_uis_(Vely, Sigmayx)->Update(fac * normal(Velx), bG_ib, 1.0);
  BG_uis_(Vely, Sigmayy)->Update(fac * normal(Vely), bG_ib, 1.0);
  BG_uis_(Vely, Sigmayz)->Update(fac * normal(Velz), bG_ib, 1.0);
  BG_uis_(Velz, Sigmazx)->Update(fac * normal(Velx), bG_ib, 1.0);
  BG_uis_(Velz, Sigmazy)->Update(fac * normal(Vely), bG_ib, 1.0);
  BG_uis_(Velz, Sigmazz)->Update(fac * normal(Velz), bG_ib, 1.0);

  rhs_s(Sigmaxx, 0)->Update(-fac * normal(Velx) * ivelint(Velx), funct, 1.0);
  rhs_s(Sigmayy, 0)->Update(-fac * normal(Vely) * ivelint(Vely), funct, 1.0);
  rhs_s(Sigmazz, 0)->Update(-fac * normal(Velz) * ivelint(Velz), funct, 1.0);

  rhs_s(Sigmaxy, 0)->Update(-fac * (normal(Vely) * ivelint(Velx) + normal(Velx) * ivelint(Vely)), funct, 1.0);
  rhs_s(Sigmaxz, 0)->Update(-fac * (normal(Velz) * ivelint(Velx) + normal(Velx) * ivelint(Velz)), funct, 1.0);
  rhs_s(Sigmayz, 0)->Update(-fac * (normal(Velz) * ivelint(Vely) + normal(Vely) * ivelint(Velz)), funct, 1.0);


  //Build cross-interface pressure-velocity coupling matrices G_uip, G_pui!

  //G_pui - from adjoint pressure consistency term
  /*
   *  /          \
   *  | q, u^i n  |
   *  \          /
   *
   */

  BG_pui_(0,Velx)->Update(fac * normal(Velx), bG_bi, 1.0);
  BG_pui_(0,Vely)->Update(fac * normal(Vely), bG_bi, 1.0);
  BG_pui_(0,Velz)->Update(fac * normal(Velz), bG_bi, 1.0);


  //G_uip - from pressure consistency term
  /*
   *  /           \
   *  | -v^i, p n  |
   *  \           /
   *
   */

  BG_uip_(Velx,0)->Update(-fac * normal(Velx), bG_ib, 1.0);
  BG_uip_(Vely,0)->Update(-fac * normal(Vely), bG_ib, 1.0);
  BG_uip_(Velz,0)->Update(-fac * normal(Velz), bG_ib, 1.0);

  // add normal interface velocity to rhs vector (pressure row)
  const double ivelnormal = ivelint.Dot(normal);

  // (Viscous) stress/Velocities

  // rhs_p_u, residual from:
  /*
   *  /         \
   *  | q, u^i n |
   *  \         /
   *
   */
  rhs_pui.Update(-fac * ivelnormal, funct, 1.0);

  // rhs_ui_p, residual from:
  /*
   *  /          \
   *  | -v^i, p n |
   *  \          /
   *
   */
  // belongs to the side and therefore contributes to rhC_ui_!
  for (int ir = 0; ir < side_nen_; ++ir)
  {
    unsigned velxrow = ir * side_numdof + Velx;
    unsigned velyrow = ir * side_numdof + Vely;
    unsigned velzrow = ir * side_numdof + Velz;
    unsigned presrow = ir * side_numdof + Pres;

    rhC_ui_(velxrow, 0) += press * fac * normal(Velx) * side_funct_(ir);
    rhC_ui_(velyrow, 0) += press * fac * normal(Vely) * side_funct_(ir);
    rhC_ui_(velzrow, 0) += press * fac * normal(Velz) * side_funct_(ir);
    rhC_ui_(presrow, 0) = 0.0;
  }

  return;
}// MHVS_buildCouplingMatrices

/*--------------------------------------------------------------------------------
 * build final coupling matrices for Mixed/Stress/Hybrid (MSH) method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::HybridLM_buildFinalCouplingMatrices(
    LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,numstressdof_> &     BinvK_ss,   ///< block inverse \f$ K^{-1}_{\sigma\sigma} \f$
    LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,bg_numdof_,numstressdof_> &        BKusInvKss, ///< block matrix \f$ K_{u\sigma} \cdot K^{-1}_{\sigma\sigma} \f$
    LINALG::BlockMatrix<LINALG::Matrix<nen_,nen_>,numstressdof_,bg_numdof_> &        BK_su,      ///< block matrix \f$ K_{\sigma u} \f$
    LINALG::BlockMatrix<LINALG::Matrix<nen_,1>, numstressdof_,1>      &              rhs_s       ///< block rhs vector \f$ rhs_{\sigma}\f$
)
{
  // final coupling matrices in block form
  LINALG::BlockMatrix<LINALG::Matrix<nen_,side_nen_>,bg_numdof_,nsd_>     BCuui;
  LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,nsd_,bg_numdof_>     BCuiu;

  // auxiliary matrices for intermediate calculation results of the matrix product
  LINALG::BlockMatrix<LINALG::Matrix<side_nen_,nen_>,nsd_,numstressdof_>  BGuisInvKss;
  LINALG::BlockMatrix<LINALG::Matrix<side_nen_,1>,nsd_,1>                 BGuisinvKssrhs_s;

  // BKusInvKss:
  // (K_us + G_us) * K_ss^-1 (MHVS) or
  // G_us * K_ss^-1 (MHCS)

  // (K_us + G_us) * K_ss^-1 G_sui (MHVS) or
  // G_us * K_ss^-1 G_sui (MHCS)
  BCuui.Multiply(BKusInvKss,BG_sui_);

  // G_uis * K_ss^-1
  BGuisInvKss.Multiply(BG_uis_,BinvK_ss);

  // G_uis * K_ss^-1 * (K_su + G_su) (MHVS) or
  // G_uis * K_ss^-1 * (K_su + G_su + K_sp) (MHCS)
  BCuiu.Multiply(BGuisInvKss,BK_su);

  // G_uis K_ss^-1 rhs_s
  BGuisinvKssrhs_s.Multiply(BGuisInvKss,rhs_s);

  // transfer the entries from Cuui,Cuiu, rhCui (in case of MHCS, coupling term ui-p is
  // included in Cuiu!)

  // loop over side velocity dof
  for (int isivel = 0; isivel < nsd_; ++ isivel)
  {
    // loop over background element dof (velocity & pressure)
    for (int ibgdof = 0; ibgdof < bg_numdof_; ++ ibgdof)
    {
      // (u-ui)
      if (BCuui.IsUsed(ibgdof,isivel))
      {
        const LINALG::Matrix<nen_,side_nen_> & bCuui = *BCuui(ibgdof,isivel);
        // loop over side element nodes
        for (int isin = 0; isin < side_nen_; ++ isin)
        {
          unsigned col = isivel + side_numdof * isin;
          // loop over background element nodes
          for (int ibgn = 0; ibgn < nen_; ++ ibgn)
          {
            unsigned row = ibgdof + bg_numdof_ * ibgn;
            C_uui_(row,col) -= bCuui(ibgn,isin);
          }
        }
      } // (u-ui)

      // (ui-u), MHCS: (ui-p)
      if (BCuiu.IsUsed(isivel,ibgdof))
      {
        const LINALG::Matrix<side_nen_,nen_> & bCuiu = *BCuiu(isivel,ibgdof);
        // loop over side element nodes
        for (int isin = 0; isin < side_nen_; ++ isin)
        {
          unsigned row = isivel + side_numdof * isin;
          // loop over background element nodes
          for (int ibgn = 0; ibgn < nen_; ++ ibgn)
          {
            unsigned col = ibgdof + bg_numdof_ * ibgn;
            C_uiu_(row,col) -= bCuiu(isin,ibgn);
          }
        }
      } // (ui-u), MHCS: (ui-p)
    }

    // add surface-based pressure coupling terms (only MHVS)

    // (ui-p)
    if (BG_uip_.IsUsed(isivel,0))
    {
      const LINALG::Matrix<side_nen_, nen_> & bGuip = *BG_uip_(isivel,0);
      // loop over side element nodes
      for (int isin = 0; isin < side_nen_; ++ isin)
      {
        unsigned row = isivel + side_numdof * isin;
        // loop over background element nodes
        for (int ibgn = 0; ibgn < nen_; ++ ibgn)
        {
          unsigned col = nsd_ + bg_numdof_ * ibgn;
          C_uiu_(row,col) = bGuip(isin,ibgn);
        }
      }
    } // (ui-p)

    // (p-ui)
    if (BG_pui_.IsUsed(0,isivel))
    {
      const LINALG::Matrix<nen_, side_nen_> & bGpui = *BG_pui_(0,isivel);
      // loop over side element nodes
      for (int isin = 0; isin < side_nen_; ++ isin)
      {
        unsigned col = isivel + side_numdof * isin;
        // loop over background element nodes
        for (int ibgn = 0; ibgn < nen_; ++ ibgn)
        {
          unsigned row = nsd_ + bg_numdof_ * ibgn;
          C_uui_(row,col) = bGpui(ibgn,isin);
        }
      }
    } // (p-ui)

    // rhs - ui
    if (BGuisinvKssrhs_s.IsUsed(isivel,0))
    {
      const LINALG::Matrix<side_nen_,1> & bGuisinvKssrhs_s = *BGuisinvKssrhs_s(isivel,0);

      // loop over side element nodes
      for (int isin = 0; isin < side_nen_; ++ isin)
      {
        unsigned row = isivel + side_numdof * isin;
        rhC_ui_(row,0) -= bGuisinvKssrhs_s(isin,0);
      }
    } // rhs - ui
  } // end loop over side velocity dof

  // finally, build G_uis & G_uis for C_uiui

  // G_sui
  // loop over block rows
  for (int ibr = 0; ibr < numstressdof_; ++ibr)
  {
    // loop over block columns (interface velocity)
    for (int ibc = 0; ibc < nsd_; ++ibc)
    {
      // extract the stress-velocity coupling submatrix
      if (BG_sui_.IsUsed(ibr,ibc))
      {
        LINALG::Matrix<nen_,side_nen_> & bGsui = * BG_sui_(ibr,ibc);

        // transfer the entries
        for (int ir =0; ir<nen_; ++ir)
        {
          unsigned stressrow = ibr + ir * numstressdof_;

          for (int ic = 0; ic < side_nen_; ++ic)
          {
            unsigned sidevelcol = ibc + ic * side_numdof;

            G_sui_(stressrow,sidevelcol) = bGsui(ir,ic);
          }
        }
      }
    }
  } // G_sui filled

  // fill G_uis_ from BG_uis_

  // loop over block columns (sigmaxx, sigmaxy, ...)
  for (int ibc = 0; ibc < numstressdof_; ++ibc)
  {
    // loop over block rows (interface velocity)
    for (int ibr = 0; ibr < nsd_; ++ibr)
    {
      if (BG_uis_.IsUsed(ibr,ibc))
      {
        LINALG::Matrix< side_nen_, nen_> & bGuis = * BG_uis_(ibr,ibc);

        // transfer the entries
        for (int ic=0; ic<nen_; ++ic)
        {
          unsigned stresscol = ibc + ic * numstressdof_;

          for (int ir=0; ir<side_nen_; ++ir)
          {
            unsigned sidevelrow = ibr + ir * side_numdof;

            G_uis_(sidevelrow,stresscol) = bGuis(ir,ic);
          }
        }
      }
    }
  } // G_uis filled

  return;
} // Hybrid_LM_buildFinalCouplingMatrices



/*--------------------------------------------------------------------------------
 * build convective/inflow stabilization for hybrid LM approach
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::HybridLM_Stab_InflowCoercivity(
    Epetra_SerialDenseMatrix &          C_uu_,                ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &          rhs_Cu_,              ///< standard bg-rhs
    const bool &                        coupling,             ///< indicates the need to assemble coupling terms for the side (yes/no)
    const bool &                        bg_mortaring,         ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
    const LINALG::Matrix<nsd_,1> &      normal,               ///< normal vector
    const double &                      timefacfac,           ///< theta*dt
    const double &                      visceff_1,            ///< viscosity in background fluid
    const double &                      visceff_2,            ///< viscosity in embedded fluid
    const double &                      kappa1,               ///< mortaring weighting
    const double &                      kappa2,               ///< mortaring weighting
    const double &                      stabfac,              ///< Nitsche penalty
    const double &                      stabfac_avg,          ///< Nitsche convective non-dimensionless stabilization factor
    const LINALG::Matrix<nen_,1> &      funct_,               ///< bg shape functions
    const LINALG::Matrix<nsd_,nen_> &   derxy_,               ///< bg deriv
    const LINALG::Matrix<nsd_,nsd_> &   vderxy_,              ///< bg deriv^n
    const LINALG::Matrix<nsd_,1> &      velint,               ///< bg u^n
    const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,    ///< Dirichlet velocity vector or prescribed jump vector
    INPAR::XFEM::ConvStabScaling        conv_stab_scaling,    ///< Inflow term strategies xfluid
    INPAR::XFEM::XFF_ConvStabScaling    xff_conv_stab_scaling ///< Inflow term strategies xfluidfluid
  )
{
  // if we don't have a reasonable convective stabilization factor, leave
  if(stabfac_avg <= 0.0) return;
  //--------------------------------------------

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  // interface velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> ivelint;
  ivelint.Multiply(eivel_,side_funct_);

  // funct_ * timefac * fac * funct_ (dyadic product)
  LINALG::Matrix<nen_,1> funct_timefacfac(true);
  funct_timefacfac.Update(timefacfac,funct_,0.0);

  //  LINALG::Matrix<nen_,1> funct_timefacfac_k1(true);
  //  funct_timefacfac_k1.Update(kappa1,funct_timefacfac,0.0);

  LINALG::Matrix<side_nen_,1> side_funct_timefacfac(true);
  side_funct_timefacfac.Update(timefacfac,side_funct_,0.0);

  LINALG::Matrix<side_nen_,1> side_funct_timefacfac_k1(true);
  side_funct_timefacfac_k1.Update(kappa1,side_funct_timefacfac,0.0);

  LINALG::Matrix<nen_,nen_> funct_dyad_timefacfac(true);
  funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct_);

  LINALG::Matrix<side_nen_,nen_> side_funct_dyad_timefacfac(true);
  side_funct_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, funct_);

  LINALG::Matrix<nen_,side_nen_> funct_side_dyad_timefacfac(true);
  funct_side_dyad_timefacfac.MultiplyNT(funct_timefacfac, side_funct_);

  LINALG::Matrix<side_nen_,side_nen_> side_side_dyad_timefacfac(true);
  side_side_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, side_funct_);

  if(!coupling)
  {
    //-----------------------------------------------------------------
    // same as Nitsche's viscous stability term
    // REMARK: this term includes inflow coercivity in case of XFSI

    NIT_Stab_ViscCoercivity(
        C_uu_,          // standard bg-bg-matrix
        rhs_Cu_,        // standard bg-rhs
        velint,
        ivelint,
        ivelint_WDBC_JUMP,
        funct_timefacfac,
        funct_dyad_timefacfac,
        side_funct_timefacfac,
        side_funct_dyad_timefacfac,
        stabfac,       //Nitsche parameter
        false,         // assemble coupling terms (Cuiu, Cuui, etc.) (yes/no)
        normal         // normal vector
    );
  }
  else
  {
    if(conv_stab_scaling == INPAR::XFEM::ConvStabScaling_inflow or
       xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
       xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty)
    {
      NIT_Stab_ViscCoercivity(C_uu_,          // standard bg-bg-matrix
                              rhs_Cu_,        // standard bg-rhs
                              velint,
                              ivelint,
                              ivelint_WDBC_JUMP,
                              funct_timefacfac,
                              funct_dyad_timefacfac,
                              side_funct_timefacfac,
                              side_funct_dyad_timefacfac,
                              stabfac,
                              true,       // assemble coupling terms (yes/no)
                              normal         // normal vector
                              );
    }

    if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_averaged or
        conv_stab_scaling == INPAR::XFEM::ConvStabScaling_inflow or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty)
    {
      NIT_Stab_ConvAveraged(C_uu_,          // standard bg-bg-matrix
                            rhs_Cu_,        // standard bg-rhs
                            velint,
                            ivelint,
                            ivelint_WDBC_JUMP,
                            funct_timefacfac,
                            funct_dyad_timefacfac,
                            side_funct_timefacfac,
                            side_funct_dyad_timefacfac,
                            stabfac_avg,
                            coupling,       // assemble coupling terms (yes/no)
                            normal         // normal vector
                          );
   }
  }
  return;
}

/*--------------------------------------------------------------------------------
 * build coupling matrices and assemble terms for Nitsche's (NIT) method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_buildCouplingMatrices(
    Epetra_SerialDenseMatrix &          C_uu_,                  ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &          rhs_Cu_,                ///< standard bg-rhs
    const bool &                        coupling,               ///< assemble coupling terms (yes/no)
    const bool &                        bg_mortaring,           ///< yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring, two-sided impl)
    const LINALG::Matrix<nsd_,1> &      normal,                 ///< normal vector
    const double &                      timefacfac,             ///< theta*dt
    const double &                      visceff_1,              ///< viscosity in background fluid
    const double &                      visceff_2,              ///< viscosity in embedded fluid
    const double &                      kappa1,                 ///< mortaring weighting
    const double &                      kappa2,                 ///< mortaring weighting
    const double &                      NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
    const double &                      avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
    const LINALG::Matrix<nen_,1> &      funct_,                 ///< bg shape functions
    const LINALG::Matrix<nsd_,nen_> &   derxy_,                 ///< bg deriv
    const LINALG::Matrix<nsd_,nsd_> &   vderxy_,                ///< bg deriv^n
    const double &                      press,                  ///< bg p^n
    const LINALG::Matrix<nsd_,1> &      velint,                 ///< bg u^n
    const LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP,      ///< Dirichlet velocity vector or prescribed jump vector
    INPAR::XFEM::ConvStabScaling        conv_stab_scaling,      ///< Inflow term strategies xfluid
    INPAR::XFEM::XFF_ConvStabScaling    xff_conv_stab_scaling   ///< Inflow term strategies xfluidfluid
  )
{

  //--------------------------------------------

  // define the coupling between two not matching grids
  // for fluidfluidcoupling
  // domain Omega^1 := Xfluid
  // domain Omega^2 := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

  // [| v |] := v1 - v2
  //  { v }  := kappa1 * v1 + kappa2 * v2 = kappa1 * v1 (for Dirichlet coupling k1=1.0, k2 = 0.0)
  //  < v >  := kappa2 * v1 + kappa1 * v2 = kappa1 * v2 (for Dirichlet coupling k1=1.0, k2 = 0.0)

  //--------------------------------------------

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  // interface velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> ivelint;
  ivelint.Multiply(eivel_,side_funct_);

  // funct_ * timefac * fac * funct_ (dyadic product)
  LINALG::Matrix<nen_,1> funct_timefacfac(true);
  funct_timefacfac.Update(timefacfac,funct_,0.0);

  LINALG::Matrix<nen_,1> funct_timefacfac_k1(true);
  funct_timefacfac_k1.Update(kappa1,funct_timefacfac,0.0);

  LINALG::Matrix<side_nen_,1> side_funct_timefacfac(true);
  side_funct_timefacfac.Update(timefacfac,side_funct_,0.0);

  LINALG::Matrix<side_nen_,1> side_funct_timefacfac_k1(true);
  side_funct_timefacfac_k1.Update(kappa1,side_funct_timefacfac,0.0);

  LINALG::Matrix<nen_,nen_> funct_dyad_timefacfac(true);
  LINALG::Matrix<nen_,nen_> funct_dyad_k1_timefacfac(true);
  funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct_);
  funct_dyad_k1_timefacfac.Update(kappa1,funct_dyad_timefacfac,0.0);

  LINALG::Matrix<side_nen_,nen_> side_funct_dyad_timefacfac(true);
  LINALG::Matrix<side_nen_,nen_> side_funct_dyad_k1_timefacfac(true);
  side_funct_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, funct_);
  side_funct_dyad_k1_timefacfac.Update(kappa1, side_funct_dyad_timefacfac, 0.0);


  double k1mu1_fac = 2.0 * timefacfac * kappa1 * visceff_1;

  LINALG::Matrix<nen_,1> e_funct_visc1_timefacfac(true);
  e_funct_visc1_timefacfac.Update(k1mu1_fac, funct_, 0.0);

  LINALG::Matrix<side_nen_,1> s_funct_visc1_timefacfac(true);
  s_funct_visc1_timefacfac.Update(k1mu1_fac, side_funct_, 0.0);

  //-----------------------------------------------------------------
  // pressure consistency term

  NIT_p_Consistency(  C_uu_,          // standard bg-bg-matrix
                      rhs_Cu_,        // standard bg-rhs
                      press,
                      funct_timefacfac_k1,
                      side_funct_timefacfac_k1,
                      side_funct_dyad_k1_timefacfac,
                      funct_dyad_k1_timefacfac,
                      coupling,       // assemble coupling terms (yes/no)
                      normal         // normal vector
  );



  double velint_normal = velint.Dot(normal);
  double ivelint_normal = ivelint.Dot(normal);
  double ivelint_WDBC_JUMP_normal = ivelint_WDBC_JUMP.Dot(normal);

  //-----------------------------------------------------------------
  // pressure adjoint consistency term

  NIT_p_AdjointConsistency(  C_uu_,          // standard bg-bg-matrix
                             rhs_Cu_,        // standard bg-rhs
                             velint_normal,
                             ivelint_normal,
                             ivelint_WDBC_JUMP_normal,
                             funct_timefacfac_k1,
                             side_funct_dyad_k1_timefacfac,
                             funct_dyad_k1_timefacfac,
                             coupling,       // assemble coupling terms (yes/no)
                             normal         // normal vector
  );


  //-----------------------------------------------------------------
  // viscous consistency term

  NIT_visc_Consistency(  C_uu_,          // standard bg-bg-matrix
                         rhs_Cu_,        // standard bg-rhs
                         derxy_,         // bg deriv
                         vderxy_,        // bg deriv^n
                         e_funct_visc1_timefacfac,
                         s_funct_visc1_timefacfac,
                         coupling,       // assemble coupling terms (yes/no)
                         normal         // normal vector
  );


  //-----------------------------------------------------------------
  // viscous adjoint consistency term

  NIT_visc_AdjointConsistency(  C_uu_,          // standard bg-bg-matrix
                                rhs_Cu_,        // standard bg-rhs
                                velint,
                                ivelint,
                                ivelint_WDBC_JUMP,
                                derxy_,         // bg deriv
                                visceff_1,
                                timefacfac,
                                e_funct_visc1_timefacfac,
                                s_funct_visc1_timefacfac,
                                coupling,       // assemble coupling terms (yes/no)
                                normal         // normal vector);
  );


  //-----------------------------------------------------------------
  // viscous stability term
  // REMARK: this term includes also inflow coercivity in case of XFSI with modified stabfac (see NIT_ComputeStabfac)

  //if (coupling == false)
    NIT_Stab_ViscCoercivity(C_uu_,          // standard bg-bg-matrix
                            rhs_Cu_,        // standard bg-rhs
                            velint,
                            ivelint,
                            ivelint_WDBC_JUMP,
                            funct_timefacfac,
                            funct_dyad_timefacfac,
                            side_funct_timefacfac,
                            side_funct_dyad_timefacfac,
                            NIT_full_stab_fac,
                            coupling,       // assemble coupling terms (yes/no)
                            normal         // normal vector
      );



    if(coupling)
    {
      if (xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
          xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged or
          xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty or
          xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty)
      {
        NIT_Stab_ConvAveraged(C_uu_,          // standard bg-bg-matrix
                              rhs_Cu_,        // standard bg-rhs
                              velint,
                              ivelint,
                              ivelint_WDBC_JUMP,
                              funct_timefacfac,
                              funct_dyad_timefacfac,
                              side_funct_timefacfac,
                              side_funct_dyad_timefacfac,
                              avg_conv_stab_fac,
                              coupling,       // assemble coupling terms (yes/no)
                              normal          // normal vector
        );

      }
    }


  // no penetration stabilization for xfsi
    /*                            \       /                       i   _     \
   |  gamma/h_K *  v*n , Du*n     | =  - |   gamma/h_K *  v*n , (u  - u)*n   |
    \                             /       \                                */
//
//  for(int ir=0; ir<nen_; ir++)
//  {
//    int idVelx = ir*(nsd_+1) + 0;
//    int idVely = ir*(nsd_+1) + 1;
//    int idVelz = ir*(nsd_+1) + 2;
//
//    // (stab * v, Du)
//    for(int ic=0; ic<nen_; ic++)
//    {
//      int iVelx = ic*(nsd_+1)+0;
//      int iVely = ic*(nsd_+1)+1;
//      int iVelz = ic*(nsd_+1)+2;
//
//      double tmp_x = funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(iVelx);
//      double tmp_y = funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(iVely);
//      double tmp_z = funct_dyad_timefacfac(ir,ic)*stabfac_conv*normal(iVelz);
//
//      C_uu_(idVelx, iVelx) += tmp_x*normal(Velx);
//      C_uu_(idVelx, iVely) += tmp_x*normal(Vely);
//      C_uu_(idVelx, iVelz) += tmp_x*normal(Velz);
//
//      C_uu_(idVely, iVelx) += tmp_y*normal(Velx);
//      C_uu_(idVely, iVely) += tmp_y*normal(Vely);
//      C_uu_(idVely, iVelz) += tmp_y*normal(Velz);
//
//      C_uu_(idVelz, iVelx) += tmp_z*normal(Velx);
//      C_uu_(idVelz, iVely) += tmp_z*normal(Vely);
//      C_uu_(idVelz, iVelz) += tmp_z*normal(Velz);
//    }
//
//    double velint_normal = velint.Dot(normal);
//    double tmp1 = funct_timefacfac(ir)*stabfac_conv*velint_normal;
//
//    // -(stab * v*n, u*n)
//    rhs_Cu_(idVelx) -= tmp1*normal(Velx);
//    rhs_Cu_(idVely) -= tmp1*normal(Vely);
//    rhs_Cu_(idVelz) -= tmp1*normal(Velz);
//
//
//    double velint_WDBC_normal = ivelint_WDBC_JUMP.Dot(normal);
//    double tmp2 = funct_timefacfac(ir)*stabfac_conv*velint_WDBC_normal;
//
//    // +(stab * v*n, u_DBC*n)
//    rhs_Cu_(idVelx) += tmp2*normal(Velx);
//    rhs_Cu_(idVely) += tmp2*normal(Vely);
//    rhs_Cu_(idVelz) += tmp2*normal(Velz);
//  }





  return;
}// NIT_buildCouplingMatrices


/*--------------------------------------------------------------------------------
 * evaluate pressure-consistency term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_p_Consistency(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const double &                           press,                            ///< pressure
    const LINALG::Matrix<nen_,1> &           funct_timefacfac_k1,              ///< funct * timefacfac *kappa1
    const LINALG::Matrix<side_nen_,1> &      side_funct_timefacfac_k1,         ///< sidefunct * timefacfac *kappa1
    const LINALG::Matrix<side_nen_,nen_> &   side_funct_dyad_k1_timefacfac,    ///< (sidefunct^T * funct) * timefacfac *kappa1
    const LINALG::Matrix<nen_,nen_> &        funct_dyad_k1_timefacfac,         ///< (funct^T * funct) * timefacfac *kappa1
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  // NIT_p_Consistency
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;

             /*                  \       /          i      \
          + |  [ v ],   {Dp}*n    | = - | [ v ], { p }* n   |
             \                   /       \                */

  //-----------------------------------------------
  //    + (v1, k1 *(Dp1)*n)
  //-----------------------------------------------


  for(int ir = 0; ir<nen_; ir++)
  {
    int idVelx = ir*bg_numdof_ + 0;
    int idVely = ir*bg_numdof_ + 1;
    int idVelz = ir*bg_numdof_ + 2;

    // (v,Dp*n)
    for(int ic =0; ic<nen_; ic++)
    {
      int iPres = ic*bg_numdof_ + 3;

      C_uu_(idVelx, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
      C_uu_(idVely, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
      C_uu_(idVelz, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
    }

    // -(v,p*n)
    double funct_k1_timefacfac_press = funct_timefacfac_k1(ir)*press;
    rhs_Cu_(idVelx,0) -= funct_k1_timefacfac_press*normal(Velx);
    rhs_Cu_(idVely,0) -= funct_k1_timefacfac_press*normal(Vely);
    rhs_Cu_(idVelz,0) -= funct_k1_timefacfac_press*normal(Velz);
  }




  if(coupling)
  {
    //-----------------------------------------------
    //    - (v2, k1 *(Dp1)*n)
    //-----------------------------------------------
    for(int ir = 0; ir<side_nen_; ir++)
    {

      int idVelx = ir*side_numdof + 0;
      int idVely = ir*side_numdof + 1;
      int idVelz = ir*side_numdof + 2;

      // (v,Dp*n)
      for(int ic =0; ic<nen_; ic++)
      {
        int iPres = ic*bg_numdof_ + 3;

        C_uiu_(idVelx, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
        C_uiu_(idVely, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
        C_uiu_(idVelz, iPres) -= side_funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
      }

      // -(v,p*n)
      double side_funct_k1_timefacfac_press = side_funct_timefacfac_k1(ir)*press;
      rhC_ui_(idVelx,0) += side_funct_k1_timefacfac_press*normal(Velx);
      rhC_ui_(idVely,0) += side_funct_k1_timefacfac_press*normal(Vely);
      rhC_ui_(idVelz,0) += side_funct_k1_timefacfac_press*normal(Velz);
    }
  }// end coupling


  return;
}// NIT_p_Consistency


/*--------------------------------------------------------------------------------
 * evaluate pressure-adjoint-consistency term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_p_AdjointConsistency(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    double &                                 velint_normal,                    ///< velocity in normal direction
    double &                                 ivelint_normal,                   ///< interface velocity in normal direction
    double &                                 ivelint_WDBC_JUMP_normal,         ///< prescribed interface velocity in normal direction
    const LINALG::Matrix<nen_,1> &           funct_timefacfac_k1,              ///< funct * timefacfac *kappa1
    const LINALG::Matrix<side_nen_,nen_> &   side_funct_dyad_k1_timefacfac,    ///< (sidefunct^T * funct) * timefacfac *kappa1
    const LINALG::Matrix<nen_,nen_> &        funct_dyad_k1_timefacfac,         ///< (funct^T * funct) * timefacfac *kappa1
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  // NIT_p_AdjointConsistency
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;

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
  //    - (q1*n, k1 *(Du1))
  //-----------------------------------------------
  for(int ir = 0; ir<nen_; ir++)
  {
    int idPres = ir*bg_numdof_ + 3;

    // -(q*n,Du)
    for(int ic =0; ic<nen_; ic++)
    {
      int iVelx = ic*bg_numdof_ + 0;
      int iVely = ic*bg_numdof_ + 1;
      int iVelz = ic*bg_numdof_ + 2;

      C_uu_(idPres, iVelx) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
      C_uu_(idPres, iVely) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
      C_uu_(idPres, iVelz) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
    }

    // (q*n,u)
    rhs_Cu_(idPres,0) += funct_timefacfac_k1(ir)*velint_normal;

    if(!coupling) // weak Dirichlet case
    {
      // -(q*n,u_DBC)
      rhs_Cu_(idPres,0) -= funct_timefacfac_k1(ir)*ivelint_WDBC_JUMP_normal;
    }
  }


  if(coupling)
  {
    //-----------------------------------------------
    //    + (q1*n, k1 *(Du2))
    //-----------------------------------------------
    for(int ir = 0; ir<nen_; ir++)
    {
      int idPres = ir*bg_numdof_ + 3;

      // -(q*n,Du)
      for(int ic =0; ic<side_nen_; ic++)
      {
        int iVelx = ic*side_numdof+0;
        int iVely = ic*side_numdof+1;
        int iVelz = ic*side_numdof+2;

        C_uui_(idPres, iVelx) += side_funct_dyad_k1_timefacfac(ic,ir)*normal(Velx);
        C_uui_(idPres, iVely) += side_funct_dyad_k1_timefacfac(ic,ir)*normal(Vely);
        C_uui_(idPres, iVelz) += side_funct_dyad_k1_timefacfac(ic,ir)*normal(Velz);
      }

      // -(q*n,u)
      rhs_Cu_(idPres,0) -= funct_timefacfac_k1(ir)*ivelint_normal;

    }
  }// end coupling

  return;
}// NIT_p_AdjointConsistency


/*--------------------------------------------------------------------------------
 * evaluate viscous-consistency term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_visc_Consistency(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const LINALG::Matrix<nsd_,nen_> &        derxy,                            ///< bg deriv
    const LINALG::Matrix<nsd_,nsd_> &        vderxy,                           ///< bg deriv^n
    const LINALG::Matrix<nen_,1> &           e_funct_visc1_timefacfac,         ///< embedded element funct *mu*timefacfac
    const LINALG::Matrix<side_nen_,1> &      s_funct_visc1_timefacfac,         ///< side element funct *mu*timefacfac
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  // viscous consistency term
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;


  /*                           \       /                   i      \
- |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
  \                            /       \                         */

  //-----------------------------------------------
  //    - (v1, (2*k1*mu1) *eps(Du1)*n)
  //-----------------------------------------------

  for(int ir = 0; ir<nen_; ir++)
  {
    int idVelx = ir*bg_numdof_ + 0;
    int idVely = ir*bg_numdof_ + 1;
    int idVelz = ir*bg_numdof_ + 2;

    for(int ic =0; ic<nen_; ic++)
    {
      int iVelx = ic*bg_numdof_ + 0;
      int iVely = ic*bg_numdof_ + 1;
      int iVelz = ic*bg_numdof_ + 2;

      // - (v1, (2*k1*mu1) *eps(Du1)*n)

      //(x,x)
      C_uu_(idVelx, iVelx) -= e_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy(Velx,ic)
                                                             + 0.5 * normal(Vely)*derxy(Vely,ic)
                                                             + 0.5 * normal(Velz)*derxy(Velz,ic)  );
      //(x,y)
      C_uu_(idVelx, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velx,ic);
      //(x,z)
      C_uu_(idVelx, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Velx,ic);

      //(y,x)
      C_uu_(idVely, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Vely,ic);
      //(y,y)
      C_uu_(idVely, iVely) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
                                                             +       normal(Vely)*derxy(Vely,ic)
                                                             + 0.5 * normal(Velz)*derxy(Velz,ic)  );
      //(y,z)
      C_uu_(idVely, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Vely,ic);

      //(z,x)
      C_uu_(idVelz, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Velz,ic);
      //(z,y)
      C_uu_(idVelz, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velz,ic);
      //(z,z)
      C_uu_(idVelz, iVelz) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
                                                             + 0.5 * normal(Vely)*derxy(Vely,ic)
                                                             +       normal(Velz)*derxy(Velz,ic)  );
    }

    // - (v1, (2*k1*mu1) *eps(Du1)*n)
    rhs_Cu_(idVelx,0) += e_funct_visc1_timefacfac(ir)*(            vderxy(Velx,Velx)                     *normal(Velx)
                                                         + 0.5 * ( vderxy(Velx,Vely) + vderxy(Vely,Velx))*normal(Vely)
                                                         + 0.5 * ( vderxy(Velx,Velz) + vderxy(Velz,Velx))*normal(Velz)  );
    rhs_Cu_(idVely,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Vely,Velx) + vderxy(Velx,Vely))*normal(Velx)
                                                         +         vderxy(Vely,Vely)                     *normal(Vely)
                                                         + 0.5 * ( vderxy(Vely,Velz) + vderxy(Velz,Vely))*normal(Velz)  );
    rhs_Cu_(idVelz,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Velz,Velx) + vderxy(Velx,Velz))*normal(Velx)
                                                         + 0.5 * ( vderxy(Velz,Vely) + vderxy(Vely,Velz))*normal(Vely)
                                                         +         vderxy(Velz,Velz)                     *normal(Velz)  );

  }


  if(coupling)
  {
  //-----------------------------------------------
  //    + (v2, (2*k1*mu1) *eps(Du1)*n)
  //-----------------------------------------------
  for(int ir = 0; ir<side_nen_; ir++)
  {
    int idVelx = ir*side_numdof + 0;
    int idVely = ir*side_numdof + 1;
    int idVelz = ir*side_numdof + 2;

    for(int ic =0; ic<nen_; ic++)
    {
      int iVelx = ic*bg_numdof_ + 0;
      int iVely = ic*bg_numdof_ + 1;
      int iVelz = ic*bg_numdof_ + 2;

      // + (v2, (2*k1*mu1) *eps(Du1)*n)

      //(x,x)
      C_uiu_(idVelx, iVelx) += s_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy(Velx,ic)
                                                              + 0.5 * normal(Vely)*derxy(Vely,ic)
                                                              + 0.5 * normal(Velz)*derxy(Velz,ic)  );
      //(x,y)
      C_uiu_(idVelx, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velx,ic);
      //(x,z)
      C_uiu_(idVelx, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Velx,ic);

      //(y,x)
      C_uiu_(idVely, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Vely,ic);
      //(y,y)
      C_uiu_(idVely, iVely) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
                                                              +       normal(Vely)*derxy(Vely,ic)
                                                              + 0.5 * normal(Velz)*derxy(Velz,ic)  );
      //(y,z)
      C_uiu_(idVely, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Vely,ic);

      //(z,x)
      C_uiu_(idVelz, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Velz,ic);
      //(z,y)
      C_uiu_(idVelz, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velz,ic);
      //(z,z)
      C_uiu_(idVelz, iVelz) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
                                                              + 0.5 * normal(Vely)*derxy(Vely,ic)
                                                              +       normal(Velz)*derxy(Velz,ic)  );
    }

    // - (v2, (2*k1*mu1) *eps(Du1)*n)
    rhC_ui_(idVelx,0) -= s_funct_visc1_timefacfac(ir)*(            vderxy(Velx,Velx)                     *normal(Velx)
                                                         + 0.5 * ( vderxy(Velx,Vely) + vderxy(Vely,Velx))*normal(Vely)
                                                         + 0.5 * ( vderxy(Velx,Velz) + vderxy(Velz,Velx))*normal(Velz)  );
    rhC_ui_(idVely,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Vely,Velx) + vderxy(Velx,Vely))*normal(Velx)
                                                         +         vderxy(Vely,Vely)                     *normal(Vely)
                                                         + 0.5 * ( vderxy(Vely,Velz) + vderxy(Velz,Vely))*normal(Velz)  );
    rhC_ui_(idVelz,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Velz,Velx) + vderxy(Velx,Velz))*normal(Velx)
                                                         + 0.5 * ( vderxy(Velz,Vely) + vderxy(Vely,Velz))*normal(Vely)
                                                         +         vderxy(Velz,Velz)                     *normal(Velz)  );

  }
  }// end coupling

  return;
}// NIT_visc_Consistency


/*--------------------------------------------------------------------------------
 * evaluate viscous-adjoint-consistency term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_visc_AdjointConsistency(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity
    const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity
    const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
    const LINALG::Matrix<nsd_,nen_> &        derxy,                            ///< bg deriv
    const double &                           visceff_1,                        ///< viscosity in background fluid element
    const double &                           timefacfac,                       ///< timefacfac
    const LINALG::Matrix<nen_,1> &           e_funct_visc1_timefacfac,         ///< embedded element funct *mu*timefacfac
    const LINALG::Matrix<side_nen_,1> &      s_funct_visc1_timefacfac,         ///< side element funct *mu*timefacfac
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  // viscous adjoint consistency term
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;


     /*                                \       /                             i   \
  - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
     \                                 /       \                                */
  // (see Burman, Fernandez 2009)
  // +1.0 symmetric
  // -1.0 antisymmetric
  const double alpha = adj_visc_scale_;

  //-----------------------------------------------
  //    - ((2*k1*mu1) *eps(v1)*n , u1)
  //-----------------------------------------------
  for(int ir = 0; ir<nen_; ir++)
  {
    int idVelx = ir*bg_numdof_ + 0;
    int idVely = ir*bg_numdof_ + 1;
    int idVelz = ir*bg_numdof_ + 2;

    // -(2mu*eps(v)*n, Du)
    for(int ic =0; ic<nen_; ic++)
    {
      int iVelx = ic*bg_numdof_ + 0;
      int iVely = ic*bg_numdof_ + 1;
      int iVelz = ic*bg_numdof_ + 2;

      //(x,x)
      C_uu_(idVelx, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy(Velx,ir)
                                                                   + 0.5 * normal(Vely)*derxy(Vely,ir)
                                                                   + 0.5 * normal(Velz)*derxy(Velz,ir)  );
      //(y,x)
      C_uu_(idVely, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velx,ir);
      //(z,x)
      C_uu_(idVelz, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Velx,ir);

      //(x,y)
      C_uu_(idVelx, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Vely,ir);
      //(y,y)
      C_uu_(idVely, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
                                                                   +       normal(Vely)*derxy(Vely,ir)
                                                                   + 0.5 * normal(Velz)*derxy(Velz,ir)  );
      //(z,y)
      C_uu_(idVelz, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Vely,ir);

      //(x,z)
      C_uu_(idVelx, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Velz,ir);
      //(y,z)
      C_uu_(idVely, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velz,ir);
      //(z,z)
      C_uu_(idVelz, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
                                                                   + 0.5 * normal(Vely)*derxy(Vely,ir)
                                                                   +       normal(Velz)*derxy(Velz,ir)  );
    }
    //  (2mu1*eps(v1)*n, u1)
    double timefacfac_visc = alpha*timefacfac*2.0*visceff_1;
    rhs_Cu_(idVelx,0) += timefacfac_visc* (     derxy(Velx,ir) *       normal(Velx) * velint(Velx)
                                              + derxy(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                              + derxy(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

    rhs_Cu_(idVely,0) += timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
                                              + derxy(Vely,ir) *       normal(Vely) * velint(Vely)
                                              + derxy(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

    rhs_Cu_(idVelz,0) += timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
                                              + derxy(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
                                              + derxy(Velz,ir) *       normal(Velz) * velint(Velz));

    if(!coupling) // weak Dirichlet case
    {
      // -(2mu*eps(v)*n, u_DBC)
      rhs_Cu_(idVelx,0) -= timefacfac_visc* (  derxy(Velx,ir) *       normal(Velx) * ivelint_WDBC_JUMP(Velx)
                                             + derxy(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
          + derxy(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Velz)));

      rhs_Cu_(idVely,0) -= timefacfac_visc* (  derxy(Velx,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
                                             + derxy(Vely,ir) *       normal(Vely) * ivelint_WDBC_JUMP(Vely)
                                             + derxy(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Vely) + normal(Vely)*ivelint_WDBC_JUMP(Velz)));

      rhs_Cu_(idVelz,0) -= timefacfac_visc* (  derxy(Velx,ir) * 0.5* (normal(Velx) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Velx))
                                             + derxy(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Vely))
                                             + derxy(Velz,ir) *       normal(Velz) * ivelint_WDBC_JUMP(Velz));

    }
  }

  if(coupling)
  {
    //-----------------------------------------------
    //    + ((2*k1*mu1) *eps(v1)*n , u2)
    //-----------------------------------------------
    for(int ir = 0; ir<nen_; ir++)
    {
      int idVelx = ir*bg_numdof_ + 0;
      int idVely = ir*bg_numdof_ + 1;
      int idVelz = ir*bg_numdof_ + 2;

      // -(2mu*eps(v1)*n, Du2)
      for(int ic =0; ic<side_nen_; ic++)
      {
        int iVelx = ic*side_numdof + 0;
        int iVely = ic*side_numdof + 1;
        int iVelz = ic*side_numdof + 2;

        //(x,x)
        C_uui_(idVelx, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy(Velx,ir)
                                                                      + 0.5 * normal(Vely)*derxy(Vely,ir)
                                                                      + 0.5 * normal(Velz)*derxy(Velz,ir)  );
        //(y,x)
        C_uui_(idVely, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velx,ir);
        //(z,x)
        C_uui_(idVelz, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Velx,ir);

        //(x,y)
        C_uui_(idVelx, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Vely,ir);
        //(y,y)
        C_uui_(idVely, iVely) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
                                                                      +       normal(Vely)*derxy(Vely,ir)
                                                                      + 0.5 * normal(Velz)*derxy(Velz,ir)  );
        //(z,y)
        C_uui_(idVelz, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Vely,ir);

        //(x,z)
        C_uui_(idVelx, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Velz,ir);
        //(y,z)
        C_uui_(idVely, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velz,ir);
        //(z,z)
        C_uui_(idVelz, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
                                                                      + 0.5 * normal(Vely)*derxy(Vely,ir)
                                                                      +       normal(Velz)*derxy(Velz,ir)  );
      }
      //  (2mu1*eps(v1)*n, u2)
      double timefacfac_visc = alpha*timefacfac*2.0*visceff_1;
      rhs_Cu_(idVelx,0) -= timefacfac_visc* (     derxy(Velx,ir) *       normal(Velx) * ivelint(Velx)
                                                + derxy(Vely,ir) * 0.5* (normal(Vely) * ivelint(Velx) + normal(Velx)*ivelint(Vely))
                                                + derxy(Velz,ir) * 0.5* (normal(Velz) * ivelint(Velx) + normal(Velx)*ivelint(Velz)));

      rhs_Cu_(idVely,0) -= timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Vely) * ivelint(Velx) + normal(Velx)*ivelint(Vely))
                                                + derxy(Vely,ir) *       normal(Vely) * ivelint(Vely)
                                                + derxy(Velz,ir) * 0.5* (normal(Velz) * ivelint(Vely) + normal(Vely)*ivelint(Velz)));

      rhs_Cu_(idVelz,0) -= timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Velx) * ivelint(Velz) + normal(Velz)*ivelint(Velx))
                                                + derxy(Vely,ir) * 0.5* (normal(Vely) * ivelint(Velz) + normal(Velz)*ivelint(Vely))
                                                + derxy(Velz,ir) *       normal(Velz) * ivelint(Velz));


    }
  }// end coupling

  return;
}// NIT_visc_AdjointConsistency


/*--------------------------------------------------------------------------------
 * evaluate stabilizing viscous term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_Stab_ViscCoercivity(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
    const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
    const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
    const LINALG::Matrix<side_nen_,1>&       side_funct_timefacfac,            ///< sidefunct^T * timefacfac
    const LINALG::Matrix<side_nen_,nen_>&    side_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
    const double &                           stabfac,                          ///< stabilization factor
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  // viscous stability term
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;

  // combined viscous and inflow stabilization for one-sided problems (XFSI)
  // gamma_combined = max(alpha*mu/hk, |u*n| )
   /*                      _        \        /                     i   _   \
  |  gamma_combined *  v , u - u     | =  - |   gamma/h_K *  v , (u  - u)   |
   \                                /        \                            */



  // for fluidfluidcoupling
  // just viscous stabilization for two-sided problems (XFF, XFFSI)
  /*                                  \        /                           i   \
 |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
  \                                   /        \                              */


  // + gamma*mu/h_K (v1, u1))
  for(int ir=0; ir<nen_; ir++)
  {
    int idVelx = ir*bg_numdof_ + 0;
    int idVely = ir*bg_numdof_ + 1;
    int idVelz = ir*bg_numdof_ + 2;

    for(int ic=0; ic<nen_; ic++)
    {
      int Velx = ic*bg_numdof_ + 0;
      int Vely = ic*bg_numdof_ + 1;
      int Velz = ic*bg_numdof_ + 2;

      double tmp = funct_dyad_timefacfac(ir,ic)*stabfac;

      C_uu_(idVelx, Velx) += tmp;
      C_uu_(idVely, Vely) += tmp;
      C_uu_(idVelz, Velz) += tmp;
    }

    double tmp = funct_timefacfac(ir)*stabfac;

    // -(stab * v, u)
    rhs_Cu_(idVelx,0) -= tmp*velint(Velx);
    rhs_Cu_(idVely,0) -= tmp*velint(Vely);
    rhs_Cu_(idVelz,0) -= tmp*velint(Velz);

    if(!coupling) // weak Dirichlet case
    {
      // +(stab * v, u_DBC)
      rhs_Cu_(idVelx,0) += tmp*ivelint_WDBC_JUMP(Velx);
      rhs_Cu_(idVely,0) += tmp*ivelint_WDBC_JUMP(Vely);
      rhs_Cu_(idVelz,0) += tmp*ivelint_WDBC_JUMP(Velz);
    }

  }

  // fluidfluid Coupling
  if(coupling)
  {
    // - gamma*mu/h_K (v1, u2))
    for(int ir=0; ir<nen_; ir++)
    {
      int idVelx = ir*bg_numdof_ + 0;
      int idVely = ir*bg_numdof_ + 1;
      int idVelz = ir*bg_numdof_ + 2;

      for(int ic=0; ic<side_nen_; ic++)
      {
        int iVelx = ic*side_numdof+0;
        int iVely = ic*side_numdof+1;
        int iVelz = ic*side_numdof+2;

        double tmp = side_funct_dyad_timefacfac(ic,ir)*stabfac;

        C_uui_(idVelx, iVelx) -= tmp;
        C_uui_(idVely, iVely) -= tmp;
        C_uui_(idVelz, iVelz) -= tmp;
      }

      double tmp = funct_timefacfac(ir)*stabfac;

      // -(stab * v, u)
      rhs_Cu_(idVelx,0) += tmp*ivelint(Velx);
      rhs_Cu_(idVely,0) += tmp*ivelint(Vely);
      rhs_Cu_(idVelz,0) += tmp*ivelint(Velz);


    }


    // - gamma*mu/h_K (v2, u1))
    for(int ir=0; ir<side_nen_; ir++)
    {
      int idVelx = ir*side_numdof + 0;
      int idVely = ir*side_numdof + 1;
      int idVelz = ir*side_numdof + 2;

      for(int ic=0; ic<nen_; ic++)
      {
        int iVelx = ic*bg_numdof_ + 0;
        int iVely = ic*bg_numdof_ + 1;
        int iVelz = ic*bg_numdof_ + 2;

        double tmp = side_funct_dyad_timefacfac(ir,ic)*stabfac;

        C_uiu_(idVelx, iVelx) -= tmp;
        C_uiu_(idVely, iVely) -= tmp;
        C_uiu_(idVelz, iVelz) -= tmp;
      }

      double tmp = side_funct_timefacfac(ir)*stabfac;

      // +(stab * v2, u1)
      rhC_ui_(idVelx,0) += tmp*velint(Velx);
      rhC_ui_(idVely,0) += tmp*velint(Vely);
      rhC_ui_(idVelz,0) += tmp*velint(Velz);

    }

    LINALG::Matrix<side_nen_,side_nen_> side_side_dyad_timefacfac(true);
    side_side_dyad_timefacfac.MultiplyNT(side_funct_timefacfac, side_funct_);

    // + gamma*mu/h_K (v2, u2))
    for(int ir=0; ir<side_nen_; ir++)
    {
      int idVelx = ir*side_numdof + 0;
      int idVely = ir*side_numdof + 1;
      int idVelz = ir*side_numdof + 2;

      for(int ic=0; ic<side_nen_; ic++)
      {
        int iVelx = ic*side_numdof+0;
        int iVely = ic*side_numdof+1;
        int iVelz = ic*side_numdof+2;

        double tmp = side_side_dyad_timefacfac(ir,ic)*stabfac;

        C_uiui_(idVelx, iVelx) += tmp;
        C_uiui_(idVely, iVely) += tmp;
        C_uiui_(idVelz, iVelz) += tmp;
      }

      double tmp = side_funct_timefacfac(ir)*stabfac;

      // -(stab * v2, u2)
      rhC_ui_(idVelx,0) -= tmp*ivelint(Velx);
      rhC_ui_(idVely,0) -= tmp*ivelint(Vely);
      rhC_ui_(idVelz,0) -= tmp*ivelint(Velz);

    }
  } // end coupling

  return;
}// NIT_Stab_ViscCoercivity

/*--------------------------------------------------------------------------------
 * evaluate stabilizing viscous term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType side_distype, const int side_numdof>
void SideImpl<distype, side_distype, side_numdof>::NIT_Stab_ConvAveraged(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
    const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
    const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
    const LINALG::Matrix<side_nen_,1>&       side_funct_timefacfac,            ///< sidefunct^T * timefacfac
    const LINALG::Matrix<side_nen_,nen_>&    side_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
    const double &                           stabfac_avg,                      ///< stabilization factor
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;

/*                                        \        /                                       i \
|  [rho * (beta * n^e)] *  { v }_m , [ Du ] | =  - |  [rho * (beta * n^e)] * { v }_m,  [ u ]  |
\ ----stab_avg-----                      /         \ ----stab_avg-----                      */


  //  [rho * (beta * n^e)] (k1*vb,ub)
  for(int ir=0; ir<nen_; ir++)
  {
    int idVelx = ir*(bg_numdof_) + 0;
    int idVely = ir*(bg_numdof_) + 1;
    int idVelz = ir*(bg_numdof_) + 2;

    for(int ic=0; ic<nen_; ic++)
    {
      int Velx = ic*(bg_numdof_)+0;
      int Vely = ic*(bg_numdof_)+1;
      int Velz = ic*(bg_numdof_)+2;

      //minus?
      double tmp = 0.5*funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uu_(idVelx, Velx) += tmp;
      C_uu_(idVely, Vely) += tmp;
      C_uu_(idVelz, Velz) += tmp;
    }

    double tmp = 0.5*funct_timefacfac(ir)*stabfac_avg;

    rhs_Cu_(idVelx,0) -= tmp*velint(Velx);
    rhs_Cu_(idVely,0) -= tmp*velint(Vely);
    rhs_Cu_(idVelz,0) -= tmp*velint(Velz);
  }

  //  -[rho * (beta * n^e)] (k1*vb,ue)
  for(int ir=0; ir<nen_; ir++)
  {
    int idVelx = ir*(bg_numdof_) + 0;
    int idVely = ir*(bg_numdof_) + 1;
    int idVelz = ir*(bg_numdof_) + 2;

    for(int ic=0; ic<side_nen_; ic++)
    {
      int Velx = ic*(side_numdof)+0;
      int Vely = ic*(side_numdof)+1;
      int Velz = ic*(side_numdof)+2;

      double tmp = -0.5*side_funct_dyad_timefacfac(ic,ir)*stabfac_avg;

      C_uui_(idVelx, Velx) += tmp;
      C_uui_(idVely, Vely) += tmp;
      C_uui_(idVelz, Velz) += tmp;
    }

    double tmp = -0.5*funct_timefacfac(ir)*stabfac_avg;

    rhs_Cu_(idVelx,0) -= tmp*ivelint(Velx);
    rhs_Cu_(idVely,0) -= tmp*ivelint(Vely);
    rhs_Cu_(idVelz,0) -= tmp*ivelint(Velz);
  }

  //  [rho * (beta * n^e)] (k1*ve,ub)
  for(int ir=0; ir<side_nen_; ir++)
  {
    int idVelx = ir*(side_numdof) + 0;
    int idVely = ir*(side_numdof) + 1;
    int idVelz = ir*(side_numdof) + 2;

    for(int ic=0; ic<nen_; ic++)
    {
      int iVelx = ic*(bg_numdof_)+0;
      int iVely = ic*(bg_numdof_)+1;
      int iVelz = ic*(bg_numdof_)+2;

      double tmp = 0.5*side_funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uiu_(idVelx, iVelx) += tmp;
      C_uiu_(idVely, iVely) += tmp;
      C_uiu_(idVelz, iVelz) += tmp;
    }

    double tmp = 0.5*side_funct_timefacfac(ir)*stabfac_avg;

    rhC_ui_(idVelx,0) -= tmp*velint(Velx);
    rhC_ui_(idVely,0) -= tmp*velint(Vely);
    rhC_ui_(idVelz,0) -= tmp*velint(Velz);
  }

  //-[rho * (beta * n^e)] (k1*ve,ue)
  for(int ir=0; ir<side_nen_; ir++)
  {
    int idVelx = ir*(side_numdof) + 0;
    int idVely = ir*(side_numdof) + 1;
    int idVelz = ir*(side_numdof) + 2;

    for(int ic=0; ic<side_nen_; ic++)
    {
      int iVelx = ic*(side_numdof)+0;
      int iVely = ic*(side_numdof)+1;
      int iVelz = ic*(side_numdof)+2;

      double tmp = -0.5*side_funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uiui_(idVelx, iVelx) += tmp;
      C_uiui_(idVely, iVely) += tmp;
      C_uiui_(idVelz, iVelz) += tmp;
    }

    double tmp = -0.5*side_funct_timefacfac(ir)*stabfac_avg;

    rhC_ui_(idVelx,0) -= tmp*ivelint(Velx);
    rhC_ui_(idVely,0) -= tmp*ivelint(Vely);
    rhC_ui_(idVelz,0) -= tmp*ivelint(Velz);
  }

  return;
}


/*--------------------------------------------------------------------------------
 * add embedded element displacements and set current element node coordinates
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::addembdisp(
    const DRT::Discretization &  embdis,       ///< cut discretization
    const std::string            state,        ///< state
    const std::vector<int>&      lm            ///< local map
    )
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = embdis.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<emb_nen_; ++inode)  // number of nodes
  {
    for(int idim=0; idim<3; ++idim) // number of dimensions
    {
      (emb_disp_)(idim,inode) = mymatrix[idim+(inode*4)]; // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
    }
  }

  // add the displacement of the interface
  for (int inode = 0; inode < emb_nen_; ++inode)
  {
    emb_xyze_(0,inode) += emb_disp_(0, inode);
    emb_xyze_(1,inode) += emb_disp_(1, inode);
    emb_xyze_(2,inode) += emb_disp_(2, inode);
  }

  return;
}


/*--------------------------------------------------------------------------------
 * extract embedded  pressure
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::getembpress(
  double& press  ///< interface velocity at embedded side
  )
{
  // interface velocity vector in gausspoint
  press = emb_funct_.Dot(emb_pres_);

  return;
}

/*--------------------------------------------------------------------------------
 * extract/set side's interface velocity gradients
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::getembvelgradint(
  LINALG::Matrix<nsd_,nsd_>& velgradint    ///< interface velocity gradients at embedded side
  )
{
  // get velocity derivatives at integration point
  velgradint.MultiplyNT(emb_vel_,emb_derxy_);

  return;
}

/*--------------------------------------------------------------------------------
 * extract/set embedded element velocity
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::emb_vel(
    const DRT::Discretization &  embdis,       ///< embedded discretization
    const std::string            state,        ///< state
    const std::vector<int>&      lm            ///< local map
    )
{

  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = embdis.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<emb_nen_; ++inode)  // number of nodes
  {
    for(int idim=0; idim<3; ++idim) // number of dimensions
    {
      (emb_vel_)(idim,inode) = mymatrix[idim+(inode*4)];  // state vector includes velocity and pressure
    }
    (emb_pres_)(inode,0) = mymatrix[3+(inode*4)];
  }

  return;
}


/*--------------------------------------------------------------------------------
 * compute embedded element's element length
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::element_length( double & hk_emb )
{
  //const double vol = VolumeViaNumIntegration<DISTYPE>(emb_xyze_);
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<emb_distype>::numNodePerElement;

  // use one point integration rule to calculate hk at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili = DRT::UTILS::intrule3D_undefined;

  switch (emb_distype)
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
  DRT::UTILS::shape_function_3D_deriv1(deriv, e(0), e(1), e(2), emb_distype);

  // get Jacobian matrix and determinant
  // xjm_ = deriv_(i,k)*xyze(j,k);
  LINALG::Matrix<3,3> xjm;
  xjm.MultiplyNT(deriv,emb_xyze_);

  const double vol = wquad * xjm.Determinant();

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  hk_emb = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  return;
}// element_length


/*--------------------------------------------------------------------------------
 * evaluate shape function, derivatives and transformation w.r.t
 * embedded element at gaussian point
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::EvaluateEmb( LINALG::Matrix<nsd_,1> & xside )
{

  // find element local position of gauss point
  GEO::CUT::Position<emb_distype> pos( emb_xyze_, xside );
  pos.Compute();

  const LINALG::Matrix<3,1> & rst_emb = pos.LocalCoordinates();

  DRT::UTILS::shape_function_3D( emb_funct_, rst_emb( 0 ), rst_emb( 1 ), rst_emb( 2 ), emb_distype );
  DRT::UTILS::shape_function_3D_deriv1( emb_deriv_, rst_emb( 0 ), rst_emb( 1 ), rst_emb( 2 ), emb_distype );


  LINALG::Matrix<nsd_,nsd_> emb_xjm(true);
  LINALG::Matrix<nsd_,nsd_> emb_xji(true);


  emb_xjm.MultiplyNT(emb_deriv_,emb_xyze_);
  emb_xji.Invert(emb_xjm);

  // compute global first derivates
  emb_derxy_.Multiply(emb_xji,emb_deriv_);

  // get pressure derivatives at integration point
  emb_prederxy_.MultiplyNN(emb_derxy_, emb_pres_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  emb_vderxy_.MultiplyNT(emb_vel_,emb_derxy_);

  return;
}// EvaluateEmb


/*--------------------------------------------------------------------------------
 * build coupling matrices and assemble terms for two-sided mortaring Nitsche's (NIT2) method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::NIT2_buildCouplingMatrices(
    Epetra_SerialDenseMatrix &    C_uu,          // standard bg-bg-matrix
    Epetra_SerialDenseVector &    rhs_Cu,        // standard bg-rhs
    bool &                        coupling,       // assemble coupling terms (yes/no)
    bool &                        bg_mortaring,   // yes: background-sided mortaring, no: coupling between two meshes (mixed mortaring)
    LINALG::Matrix<nsd_,1> &      normal,         // normal vector
    const double                  timefacfac,     // theta*dt
    const double                  visceff_1,      // viscosity in background fluid
    const double                  visceff_2,      // viscosity in embedded fluid
    double &                      kappa1,         // mortaring weighting
    double &                      kappa2,         // mortaring weighting
    const double &                NIT_full_stab_fac,      ///< full Nitsche's penalty term scaling (viscous+convective part)
    const double &                avg_conv_stab_fac,      ///< scaling of the convective average coupling term for fluidfluid problems
    bool &                        velgrad_interface_stab,// penalty term for velocity gradients at the interface
    double &                      velgrad_interface_fac, //stabilization fac for velocity gradients at the interface
    bool &                        presscoupling_interface_stab,// penalty term for pressure coupling at the interface
    double &                      presscoupling_interface_fac,//stabilization fac for pressure coupling at the interface
    LINALG::Matrix<nen_,1> &      funct,         // bg shape functions
    LINALG::Matrix<nsd_,nen_> &   derxy,         // bg deriv
    LINALG::Matrix<nsd_,nsd_> &   vderxy,        // bg deriv^n
    double &                      press,          // bg pressure at integration point
    LINALG::Matrix<nsd_,1> &      velint,         // bg u^n
    LINALG::Matrix<nsd_,1> &      ivelint_WDBC_JUMP, // Dirichlet velocity vector or prescribed jump vector
    INPAR::XFEM::XFF_ConvStabScaling  xff_conv_stab_scaling // Inflow term strategies
  )
{
  const unsigned Velx = 0;
  const unsigned Vely = 1;
  const unsigned Velz = 2;
  //const unsigned Pres = 3;

  //TODO: cleanup

  //--------------------------------------------
  // define the coupling between two not matching grids
  // for fluidfluidcoupling
  // domain Omega^1 := Xfluid
  // domain Omega^2 := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

  // [| v |] := v1 - v2
  //  { v }  := kappa1 * v1 + kappa2 * v2 = kappa1 * v1 (for Dirichlet coupling k1=1.0, k2 = 0.0)
  //  < v >  := kappa2 * v1 + kappa1 * v2 = kappa1 * v2 (for Dirichlet coupling k1=1.0, k2 = 0.0)


  double k1mu1_fac = 2.0 * timefacfac * kappa1 * visceff_1;
  double k2mu2_fac = 2.0 * timefacfac * kappa2 * visceff_2;


  //--------------------------------------------
  // get fields at interface


  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  // interface velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> emb_velint(true);
  emb_velint.Multiply(emb_vel_,emb_funct_);

  double emb_press = emb_funct_.Dot(emb_pres_);

  // funct_ * timefac * fac * funct_ (dyadic product)
  LINALG::Matrix<nen_,1> funct_timefacfac(true);
  funct_timefacfac.Update(timefacfac,funct,0.0);

  LINALG::Matrix<emb_nen_,1> emb_funct_timefacfac(true);
  emb_funct_timefacfac.Update(timefacfac,emb_funct_,0.0);

  LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_timefacfac(true);
  emb_emb_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, emb_funct_);

  LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_k2_timefacfac(true);
  emb_emb_dyad_k2_timefacfac.Update(kappa2, emb_emb_dyad_timefacfac, 0.0);


  LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_timefacfac(true);
  LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k1_timefacfac(true);
  emb_funct_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, funct);
  emb_funct_dyad_k1_timefacfac.Update(kappa1, emb_funct_dyad_timefacfac, 0.0);


  LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k2_timefacfac(true);
  emb_funct_dyad_k2_timefacfac.Update(kappa2, emb_funct_dyad_timefacfac, 0.0);

	          /*                  \       /          i      \
	        + |  [ v ], + {Dp}*n    | = - | [ v ], { p }* n   |
	          \                   /       \                */

  //-----------------------------------------------
  //    + (v1, k1 *(Dp1)*n)
  //-----------------------------------------------
  LINALG::Matrix<nen_,nen_> funct_dyad_timefacfac(true);
  LINALG::Matrix<nen_,nen_> funct_dyad_k1_timefacfac(true);
  funct_dyad_timefacfac.MultiplyNT(funct_timefacfac, funct);
  funct_dyad_k1_timefacfac.Update(kappa1,funct_dyad_timefacfac,0.0);
  for(int ir = 0; ir<nen_; ir++)
  {
    int idVelx = ir*(nsd_+1) + 0;
	int idVely = ir*(nsd_+1) + 1;
	int idVelz = ir*(nsd_+1) + 2;

	// (v,Dp*n)
	for(int ic =0; ic<nen_; ic++)
	{
	  int iPres = ic*(nsd_+1)+3;

	  C_uu(idVelx, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
	  C_uu(idVely, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
	  C_uu(idVelz, iPres) += funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
	 }

	 // -(v,p*n)
	 double funct_k1_timefacfac_press = funct_timefacfac(ir)*press*kappa1;
	 rhs_Cu(idVelx,0) -= funct_k1_timefacfac_press*normal(Velx);
	 rhs_Cu(idVely,0) -= funct_k1_timefacfac_press*normal(Vely);
	 rhs_Cu(idVelz,0) -= funct_k1_timefacfac_press*normal(Velz);
	}



	if(coupling)
	{
	  //-----------------------------------------------
	  //    - (v2, k1 *(Dp1)*n)
	  //-----------------------------------------------
	  for(int ir = 0; ir<emb_nen_; ir++)
	  {
	    int idVelx = ir*(nsd_+1) + 0;
	    int idVely = ir*(nsd_+1) + 1;
	    int idVelz = ir*(nsd_+1) + 2;

	    // (v,Dp*n)
	    for(int ic =0; ic<nen_; ic++)
	    {
	      int iPres = ic*(nsd_+1)+3;

	      C_uiu_(idVelx, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
	      C_uiu_(idVely, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
	      C_uiu_(idVelz, iPres) -= emb_funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
	     }

	    // -(v,p*n)
	    double emb_funct_k1_timefacfac_press = emb_funct_timefacfac(ir)*press*kappa1;
	    rhC_ui_(idVelx,0) += emb_funct_k1_timefacfac_press*normal(Velx);
	    rhC_ui_(idVely,0) += emb_funct_k1_timefacfac_press*normal(Vely);
	    rhC_ui_(idVelz,0) += emb_funct_k1_timefacfac_press*normal(Velz);
	   }


	   if(!bg_mortaring)
	   {
	     LINALG::Matrix<emb_nen_,nen_> emb_funct_dyad_k2_timefacfac(true);
	     emb_funct_dyad_k2_timefacfac.Update(kappa2, emb_funct_dyad_timefacfac, 0.0);

	     //-----------------------------------------------
	     //    + (v1, k2 *(Dp2)*n)
	     //-----------------------------------------------

	     for(int ir = 0; ir<nen_; ir++)
	     {
	       int idVelx = ir*(nsd_+1) + 0;
	       int idVely = ir*(nsd_+1) + 1;
	       int idVelz = ir*(nsd_+1) + 2;

	       // (v,Dp*n)
	       for(int ic =0; ic<emb_nen_; ic++)
	       {
	         int iPres = ic*(nsd_+1)+3;

	         C_uui_(idVelx, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Velx);
	         C_uui_(idVely, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Vely);
	         C_uui_(idVelz, iPres) += emb_funct_dyad_k2_timefacfac(ic,ir)*normal(Velz);
	       }

	       // -(v,p*n)
	       double funct_k2_timefacfac_press = funct_timefacfac(ir)*emb_press*kappa2;
	       rhs_Cu(idVelx,0) -= funct_k2_timefacfac_press*normal(Velx);
	       rhs_Cu(idVely,0) -= funct_k2_timefacfac_press*normal(Vely);
	       rhs_Cu(idVelz,0) -= funct_k2_timefacfac_press*normal(Velz);
	     }




	     //-----------------------------------------------
	     //    - (v2, k1 *(Dp2)*n)
	     //-----------------------------------------------
	     for(int ir = 0; ir<emb_nen_; ir++)
	     {
	       int idVelx = ir*(nsd_+1) + 0;
	       int idVely = ir*(nsd_+1) + 1;
	       int idVelz = ir*(nsd_+1) + 2;

	       // (v,Dp*n)
	       for(int ic =0; ic<nen_; ic++)
	       {
	         int iPres = ic*(nsd_+1)+3;

	          C_uiui_(idVelx, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Velx);
	          C_uiui_(idVely, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Vely);
	          C_uiui_(idVelz, iPres) -= emb_emb_dyad_k2_timefacfac(ir,ic)*normal(Velz);
	       }

	        // -(v,p*n)
	        double emb_funct_k2_timefacfac_press = emb_funct_timefacfac(ir)*emb_press*kappa2;
	        rhC_ui_(idVelx,0) += emb_funct_k2_timefacfac_press*normal(Velx);
	        rhC_ui_(idVely,0) += emb_funct_k2_timefacfac_press*normal(Vely);
	        rhC_ui_(idVelz,0) += emb_funct_k2_timefacfac_press*normal(Velz);
	      }

	    } // end !bgmortaring



	  }// end coupling



	    /*                   \     /              i   \
	 - |  { q }*n ,[ Du ]     | = |  { q }*n  ,[ u ]   |
	    \                    /     \                 */
    // the sign of the pressure adjoint consistency term is opposite to the sign of the pressure consistency term (interface)

	  //-----------------------------------------------
	  //    - (q1*n, k1 *(Du1))
	  //-----------------------------------------------
	  for(int ir = 0; ir<nen_; ir++)
	  {
	    int idPres = ir*(nsd_+1) + 3;

	    // -(q*n,Du)
	    for(int ic =0; ic<nen_; ic++)
	    {
	      int iVelx = ic*(nsd_+1)+0;
	      int iVely = ic*(nsd_+1)+1;
	      int iVelz = ic*(nsd_+1)+2;

	      C_uu(idPres, iVelx) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Velx);
	      C_uu(idPres, iVely) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Vely);
	      C_uu(idPres, iVelz) -= funct_dyad_k1_timefacfac(ir,ic)*normal(Velz);
	    }

	    // (q*n,u)
	    double velint_normal = velint.Dot(normal);
	    rhs_Cu(idPres,0) += funct_timefacfac(ir)*kappa1*velint_normal;

	    //        if(!coupling) // weak Dirichlet case
	    //        {
	    //          // -(q*n,u_DBC)
	    //          double ivelint_WDBC_JUMP_normal = ivelint_WDBC_JUMP.Dot(normal);
	    //          rhs_Cu(idPres,0) -= funct_timefacfac(ir)*ivelint_WDBC_JUMP_normal;
	    //        }
	  }

	  if(coupling)
	  {
	    //-----------------------------------------------
	    //    + (q1*n, k1 *(Du2))
	    //-----------------------------------------------
	    for(int ir = 0; ir<nen_; ir++)
	    {
	      int idPres = ir*(nsd_+1) + 3;

	      // -(q*n,Du)
	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uui_(idPres, iVelx) += emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Velx);
	        C_uui_(idPres, iVely) += emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Vely);
	        C_uui_(idPres, iVelz) += emb_funct_dyad_k1_timefacfac(ic,ir)*normal(Velz);
	      }

	      // -(q*n,u)
	      double emb_velint_normal = emb_velint.Dot(normal);
	      rhs_Cu(idPres,0) -= funct_timefacfac(ir)*kappa1*emb_velint_normal;

	    }
	  }// end coupling

	  if(!bg_mortaring)
	  {
	    //-----------------------------------------------
	    //    - (q2*n, k2 *(Du1))
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idPres = ir*(nsd_+1) + 3;

	      // -(q*n,Du)
	      for(int ic =0; ic<nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uiu_(idPres, iVelx) -= emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Velx);
	        C_uiu_(idPres, iVely) -= emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Vely);
	        C_uiu_(idPres, iVelz) -= emb_funct_dyad_k2_timefacfac(ir,ic)*normal(Velz);
	      }

	      // (q*n,u)
	      double velint_normal = velint.Dot(normal);
	      rhC_ui_(idPres,0) += emb_funct_timefacfac(ir)*kappa2*velint_normal;


	    }


	    //-----------------------------------------------
	    //    + (q2*n, k2 *(Du2))
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idPres = ir*(nsd_+1) + 3;

	      // -(q*n,Du)
	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uiui_(idPres, iVelx) += emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Velx);
	        C_uiui_(idPres, iVely) += emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Vely);
	        C_uiui_(idPres, iVelz) += emb_emb_dyad_k2_timefacfac(ic,ir)*normal(Velz);
	      }

	      // -(q*n,u)
	      double emb_velint_normal = emb_velint.Dot(normal);
	      rhC_ui_(idPres,0) -= emb_funct_timefacfac(ir)*kappa2*emb_velint_normal;

	    }

	  }


	  /*                           \       /                   i       \
	- |  [ v ],  { 2mu eps(u) }*n    | = + | [ v ],  { 2mu eps(u ) }*n  |
	  \                            /       \                          */

	  //-----------------------------------------------
	  //    - (v1, (2*k1*mu1) *eps(Du1)*n)
	  //-----------------------------------------------


	  LINALG::Matrix<nen_,1> e_funct_visc1_timefacfac(true);
	  e_funct_visc1_timefacfac.Update(k1mu1_fac, funct, 0.0);

	  LINALG::Matrix<nen_,1> e_funct_visc2_timefacfac(true);
	  e_funct_visc2_timefacfac.Update(k2mu2_fac, funct, 0.0);

	  //            LINALG::Matrix<emb_nen_,1> s_funct_visc_timefacfac(true);
	  //            s_funct_visc_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);

	  for(int ir = 0; ir<nen_; ir++)
	  {
	    int idVelx = ir*(nsd_+1) + 0;
	    int idVely = ir*(nsd_+1) + 1;
	    int idVelz = ir*(nsd_+1) + 2;

	    for(int ic =0; ic<nen_; ic++)
	    {
	      int iVelx = ic*(nsd_+1)+0;
	      int iVely = ic*(nsd_+1)+1;
	      int iVelz = ic*(nsd_+1)+2;


	      // - (v1, (2*k1*mu1) *eps(Du1)*n)

	      //(x,x)
	      C_uu(idVelx, iVelx) -= e_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy(Velx,ic)
	          + 0.5 * normal(Vely)*derxy(Vely,ic)
	          + 0.5 * normal(Velz)*derxy(Velz,ic)  );
	      //(x,y)
	      C_uu(idVelx, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velx,ic);
	      //(x,z)
	      C_uu(idVelx, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Velx,ic);

	      //(y,x)
	      C_uu(idVely, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Vely,ic);
	      //(y,y)
	      C_uu(idVely, iVely) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
	          +       normal(Vely)*derxy(Vely,ic)
	          + 0.5 * normal(Velz)*derxy(Velz,ic)  );
	      //(y,z)
	      C_uu(idVely, iVelz) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Vely,ic);

	      //(z,x)
	      C_uu(idVelz, iVelx) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Velz,ic);
	      //(z,y)
	      C_uu(idVelz, iVely) -= e_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velz,ic);
	      //(z,z)
	      C_uu(idVelz, iVelz) -= e_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
	          + 0.5 * normal(Vely)*derxy(Vely,ic)
	          +       normal(Velz)*derxy(Velz,ic)  );
	    }

	    // - (v1, (2*k1*mu1) *eps(Du1)*n)
	    rhs_Cu(idVelx,0) += e_funct_visc1_timefacfac(ir)*(            vderxy(Velx,Velx)                      *normal(Velx)
	        + 0.5 * ( vderxy(Velx,Vely) + vderxy(Vely,Velx))*normal(Vely)
	        + 0.5 * ( vderxy(Velx,Velz) + vderxy(Velz,Velx))*normal(Velz)  );
	    rhs_Cu(idVely,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Vely,Velx) + vderxy(Velx,Vely))*normal(Velx)
	        +         vderxy(Vely,Vely)                      *normal(Vely)
	        + 0.5 * ( vderxy(Vely,Velz) + vderxy(Velz,Vely))*normal(Velz)  );
	    rhs_Cu(idVelz,0) += e_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Velz,Velx) + vderxy(Velx,Velz))*normal(Velx)
	        + 0.5 * ( vderxy(Velz,Vely) + vderxy(Vely,Velz))*normal(Vely)
	        +         vderxy(Velz,Velz)                      *normal(Velz)  );

	  }


	  LINALG::Matrix<emb_nen_,1> s_funct_visc1_timefacfac(true);
	  s_funct_visc1_timefacfac.Update(k1mu1_fac, emb_funct_, 0.0);
	  LINALG::Matrix<emb_nen_,1> s_funct_visc2_timefacfac(true);
	  s_funct_visc2_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);
	  if(coupling)
	  {
	    //-----------------------------------------------
	    //    + (v2, (2*k1*mu1) *eps(Du1)*n)
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;
	      //int idPres = ir*(nsd_+1) + 3;


	      for(int ic =0; ic<nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;
	        //int iPres = ic*(nsd_+1)+3;


	        // + (v2, (2*k1*mu1) *eps(Du1)*n)

	        //(x,x)
	        C_uiu_(idVelx, iVelx) += s_funct_visc1_timefacfac(ir)*(         normal(Velx)*derxy(Velx,ic)
	            + 0.5 * normal(Vely)*derxy(Vely,ic)
	            + 0.5 * normal(Velz)*derxy(Velz,ic)  );
	        //(x,y)
	        C_uiu_(idVelx, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velx,ic);
	        //(x,z)
	        C_uiu_(idVelx, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Velx,ic);

	        //(y,x)
	        C_uiu_(idVely, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Vely,ic);
	        //(y,y)
	        C_uiu_(idVely, iVely) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
	            +       normal(Vely)*derxy(Vely,ic)
	            + 0.5 * normal(Velz)*derxy(Velz,ic)  );
	        //(y,z)
	        C_uiu_(idVely, iVelz) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velz)*derxy(Vely,ic);

	        //(z,x)
	        C_uiu_(idVelz, iVelx) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Velx)*derxy(Velz,ic);
	        //(z,y)
	        C_uiu_(idVelz, iVely) += s_funct_visc1_timefacfac(ir)*    0.5 * normal(Vely)*derxy(Velz,ic);
	        //(z,z)
	        C_uiu_(idVelz, iVelz) += s_funct_visc1_timefacfac(ir)*(   0.5 * normal(Velx)*derxy(Velx,ic)
	            + 0.5 * normal(Vely)*derxy(Vely,ic)
	            +       normal(Velz)*derxy(Velz,ic)  );
	      }

	      // - (v2, (2*k1*mu1) *eps(Du1)*n)
	      rhC_ui_(idVelx,0) -= s_funct_visc1_timefacfac(ir)*(            vderxy(Velx,Velx)                      *normal(Velx)
	          + 0.5 * ( vderxy(Velx,Vely) + vderxy(Vely,Velx))*normal(Vely)
	          + 0.5 * ( vderxy(Velx,Velz) + vderxy(Velz,Velx))*normal(Velz)  );
	      rhC_ui_(idVely,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Vely,Velx) + vderxy(Velx,Vely))*normal(Velx)
	          +         vderxy(Vely,Vely)                      *normal(Vely)
	          + 0.5 * ( vderxy(Vely,Velz) + vderxy(Velz,Vely))*normal(Velz)  );
	      rhC_ui_(idVelz,0) -= s_funct_visc1_timefacfac(ir)*(    0.5 * ( vderxy(Velz,Velx) + vderxy(Velx,Velz))*normal(Velx)
	          + 0.5 * ( vderxy(Velz,Vely) + vderxy(Vely,Velz))*normal(Vely)
	          +         vderxy(Velz,Velz)                      *normal(Velz)  );

	    }
	  }// end coupling

	  if(!bg_mortaring)
	  {
	    //-----------------------------------------------
	    //    - (v1, (2*k2*mu2) *eps(Du2)*n)
	    //-----------------------------------------------


	    LINALG::Matrix<nen_,1> e_funct_visc2_timefacfac(true);
	    e_funct_visc2_timefacfac.Update(k2mu2_fac, funct, 0.0);


	    for(int ir = 0; ir<nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;


	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        // - (v1, (2*k2*mu2) *eps(Du2)*n)

	        //(x,x)
	        C_uui_(idVelx, iVelx) -= e_funct_visc2_timefacfac(ir)*(         normal(Velx)*emb_derxy_(Velx,ic)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
	        //(x,y)
	        C_uui_(idVelx, iVely) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velx,ic);
	        //(x,z)
	        C_uui_(idVelx, iVelz) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Velx,ic);

	        //(y,x)
	        C_uui_(idVely, iVelx) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Vely,ic);
	        //(y,y)
	        C_uui_(idVely, iVely) -= e_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
	            +       normal(Vely)*emb_derxy_(Vely,ic)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
	        //(y,z)
	        C_uui_(idVely, iVelz) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Vely,ic);

	        //(z,x)
	        C_uui_(idVelz, iVelx) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Velz,ic);
	        //(z,y)
	        C_uui_(idVelz, iVely) -= e_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velz,ic);
	        //(z,z)
	        C_uui_(idVelz, iVelz) -= e_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
	            +       normal(Velz)*emb_derxy_(Velz,ic)  );
	      }

	      // - (v1, (2*k2*mu2) *eps(u2)*n)
	      rhs_Cu(idVelx,0) += e_funct_visc2_timefacfac(ir)*(            emb_vderxy_(Velx,Velx)                          *normal(Velx)
	          + 0.5 * ( emb_vderxy_(Velx,Vely) + emb_vderxy_(Vely,Velx))*normal(Vely)
	          + 0.5 * ( emb_vderxy_(Velx,Velz) + emb_vderxy_(Velz,Velx))*normal(Velz)  );
	      rhs_Cu(idVely,0) += e_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Vely,Velx) + emb_vderxy_(Velx,Vely))*normal(Velx)
	          +         emb_vderxy_(Vely,Vely)                          *normal(Vely)
	          + 0.5 * ( emb_vderxy_(Vely,Velz) + emb_vderxy_(Velz,Vely))*normal(Velz)  );
	      rhs_Cu(idVelz,0) += e_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Velz,Velx) + emb_vderxy_(Velx,Velz))*normal(Velx)
	          + 0.5 * ( emb_vderxy_(Velz,Vely) + emb_vderxy_(Vely,Velz))*normal(Vely)
	          +         emb_vderxy_(Velz,Velz)                          *normal(Velz)  );

	    }


	    LINALG::Matrix<emb_nen_,1> s_funct_visc2_timefacfac(true);
	    s_funct_visc2_timefacfac.Update(k2mu2_fac, emb_funct_, 0.0);

	    //-----------------------------------------------
	    //    + (v2, (2*k2*mu2) *eps(Du2)*n)
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;
	      //int idPres = ir*(nsd_+1) + 3;


	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;
	        //int iPres = ic*(nsd_+1)+3;


	        // + (v2, (2*k2*mu2) *eps(Du2)*n)

	        //(x,x)
	        C_uiui_(idVelx, iVelx) += s_funct_visc2_timefacfac(ir)*(         normal(Velx)*emb_derxy_(Velx,ic)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
	        //(x,y)
	        C_uiui_(idVelx, iVely) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velx,ic);
	        //(x,z)
	        C_uiui_(idVelx, iVelz) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Velx,ic);

	        //(y,x)
	        C_uiui_(idVely, iVelx) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Vely,ic);
	        //(y,y)
	        C_uiui_(idVely, iVely) += s_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
	            +       normal(Vely)*emb_derxy_(Vely,ic)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ic)  );
	        //(y,z)
	        C_uiui_(idVely, iVelz) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velz)*emb_derxy_(Vely,ic);

	        //(z,x)
	        C_uiui_(idVelz, iVelx) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Velx)*emb_derxy_(Velz,ic);
	        //(z,y)
	        C_uiui_(idVelz, iVely) += s_funct_visc2_timefacfac(ir)*    0.5 * normal(Vely)*emb_derxy_(Velz,ic);
	        //(z,z)
	        C_uiui_(idVelz, iVelz) += s_funct_visc2_timefacfac(ir)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ic)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ic)
	            +       normal(Velz)*emb_derxy_(Velz,ic)  );
	      }

	      // - (v2, (2*k2*mu2) *eps(Du2)*n)
	      rhC_ui_(idVelx,0) -= s_funct_visc2_timefacfac(ir)*(            emb_vderxy_(Velx,Velx)                          *normal(Velx)
	          + 0.5 * ( emb_vderxy_(Velx,Vely) + emb_vderxy_(Vely,Velx))*normal(Vely)
	          + 0.5 * ( emb_vderxy_(Velx,Velz) + emb_vderxy_(Velz,Velx))*normal(Velz)  );
	      rhC_ui_(idVely,0) -= s_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Vely,Velx) + emb_vderxy_(Velx,Vely))*normal(Velx)
	          +         emb_vderxy_(Vely,Vely)                          *normal(Vely)
	          + 0.5 * ( emb_vderxy_(Vely,Velz) + emb_vderxy_(Velz,Vely))*normal(Velz)  );
	      rhC_ui_(idVelz,0) -= s_funct_visc2_timefacfac(ir)*(    0.5 * ( emb_vderxy_(Velz,Velx) + emb_vderxy_(Velx,Velz))*normal(Velx)
	          + 0.5 * ( emb_vderxy_(Velz,Vely) + emb_vderxy_(Vely,Velz))*normal(Vely)
	          +         emb_vderxy_(Velz,Velz)                          *normal(Velz)  );

	    }

	  }


      /*                                \       /                             i   \
    - |  alpha* { 2mu*eps(v) }*n , [ Du ] |  =  |  alpha* { 2mu eps(v) }*n ,[ u ]   |
      \                                 /       \                                */
    // (see Burman, Fernandez 2009)
    // +1.0 symmetric
    // -1.0 antisymmetric
	  const double alpha = adj_visc_scale_;


	  //-----------------------------------------------
	  //    - ((2*k1*mu1) *eps(v1)*n , u1)
	  //-----------------------------------------------
	  for(int ir = 0; ir<nen_; ir++)
	  {
	    int idVelx = ir*(nsd_+1) + 0;
	    int idVely = ir*(nsd_+1) + 1;
	    int idVelz = ir*(nsd_+1) + 2;

	    // -(2mu*eps(v)*n, Du)
	    for(int ic =0; ic<nen_; ic++)
	    {
	      int iVelx = ic*(nsd_+1)+0;
	      int iVely = ic*(nsd_+1)+1;
	      int iVelz = ic*(nsd_+1)+2;

	      //(x,x)
	      C_uu(idVelx, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy(Velx,ir)
	          + 0.5 * normal(Vely)*derxy(Vely,ir)
	          + 0.5 * normal(Velz)*derxy(Velz,ir)  );
	      //(y,x)
	      C_uu(idVely, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velx,ir);
	      //(z,x)
	      C_uu(idVelz, iVelx) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Velx,ir);

	      //(x,y)
	      C_uu(idVelx, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Vely,ir);
	      //(y,y)
	      C_uu(idVely, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
	          +       normal(Vely)*derxy(Vely,ir)
	          + 0.5 * normal(Velz)*derxy(Velz,ir)  );
	      //(z,y)
	      C_uu(idVelz, iVely) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Vely,ir);

	      //(x,z)
	      C_uu(idVelx, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Velz,ir);
	      //(y,z)
	      C_uu(idVely, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velz,ir);
	      //(z,z)
	      C_uu(idVelz, iVelz) -= alpha*e_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
	          + 0.5 * normal(Vely)*derxy(Vely,ir)
	          +       normal(Velz)*derxy(Velz,ir)  );
	    }
	    //  (2mu1k1*eps(v1)*n, u1)
	    double timefacfac_visc = alpha*timefacfac*2.0*visceff_1*kappa1;
	    rhs_Cu(idVelx,0) += timefacfac_visc* (     derxy(Velx,ir) *       normal(Velx) * velint(Velx)
	        + derxy(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
	        + derxy(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

	    rhs_Cu(idVely,0) += timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
	        + derxy(Vely,ir) *       normal(Vely) * velint(Vely)
	        + derxy(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

	    rhs_Cu(idVelz,0) += timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
	        + derxy(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
	        + derxy(Velz,ir) *       normal(Velz) * velint(Velz));

	    //       if(!coupling) // weak Dirichlet case
	        //       {
	      //         // -(2mu*eps(v)*n, u_DBC)
	    //         rhs_Cu(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * ivelint_WDBC_JUMP(Velx)
	    //                                                + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
	    //                                                + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Velz)));
	    //
	    //         rhs_Cu_(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velx) + normal(Velx)*ivelint_WDBC_JUMP(Vely))
	    //                                                + derxy_(Vely,ir) *       normal(Vely) * ivelint_WDBC_JUMP(Vely)
	    //                                                + derxy_(Velz,ir) * 0.5* (normal(Velz) * ivelint_WDBC_JUMP(Vely) + normal(Vely)*ivelint_WDBC_JUMP(Velz)));
	    //
	    //         rhs_Cu_(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Velx))
	    //                                                + derxy_(Vely,ir) * 0.5* (normal(Vely) * ivelint_WDBC_JUMP(Velz) + normal(Velz)*ivelint_WDBC_JUMP(Vely))
	    //                                                + derxy_(Velz,ir) *       normal(Velz) * ivelint_WDBC_JUMP(Velz));
	    //
	    //       }
	  }

	  if(coupling)
	  {
	    //-----------------------------------------------
	    //    + ((2*k1*mu1) *eps(v1)*n , u2)
	    //-----------------------------------------------
	    for(int ir = 0; ir<nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      // -(2mu*eps(v1)*n, Du2)
	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        //(x,x)
	        C_uui_(idVelx, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*(         normal(Velx)*derxy(Velx,ir)
	            + 0.5 * normal(Vely)*derxy(Vely,ir)
	            + 0.5 * normal(Velz)*derxy(Velz,ir)  );
	        //(y,x)
	        C_uui_(idVely, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velx,ir);
	        //(z,x)
	        C_uui_(idVelz, iVelx) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Velx,ir);

	        //(x,y)
	        C_uui_(idVelx, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Vely,ir);
	        //(y,y)
	        C_uui_(idVely, iVely) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
	            +       normal(Vely)*derxy(Vely,ir)
	            + 0.5 * normal(Velz)*derxy(Velz,ir)  );
	        //(z,y)
	        C_uui_(idVelz, iVely) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velz)*derxy(Vely,ir);

	        //(x,z)
	        C_uui_(idVelx, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Velx)*derxy(Velz,ir);
	        //(y,z)
	        C_uui_(idVely, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*    0.5 * normal(Vely)*derxy(Velz,ir);
	        //(z,z)
	        C_uui_(idVelz, iVelz) += alpha*s_funct_visc1_timefacfac(ic)*(   0.5 * normal(Velx)*derxy(Velx,ir)
	            + 0.5 * normal(Vely)*derxy(Vely,ir)
	            +       normal(Velz)*derxy(Velz,ir)  );
	      }
	      //  (2k1mu1*eps(v1)*n, u2)
	      double timefacfac_visc = alpha*timefacfac*2.0*visceff_1*kappa1;
	      rhs_Cu(idVelx,0) -= timefacfac_visc* (     derxy(Velx,ir) *       normal(Velx) * emb_velint(Velx)
	          + derxy(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
	          + derxy(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Velx) + normal(Velx)*emb_velint(Velz)));

	      rhs_Cu(idVely,0) -= timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
	          + derxy(Vely,ir) *       normal(Vely) * emb_velint(Vely)
	          + derxy(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Vely) + normal(Vely)*emb_velint(Velz)));

	      rhs_Cu(idVelz,0) -= timefacfac_visc* (     derxy(Velx,ir) * 0.5* (normal(Velx) * emb_velint(Velz) + normal(Velz)*emb_velint(Velx))
	          + derxy(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velz) + normal(Velz)*emb_velint(Vely))
	          + derxy(Velz,ir) *       normal(Velz) * emb_velint(Velz));

	      //         // -(2mu*eps(v)*n, u_DBC)
	    //         elevec1_epetra(idVelx,0) -= timefacfac_visc* (  derxy_(Velx,ir) *       normal(Velx) * velint_WDBC(Velx)
	        //                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
	      //                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Velz)));
	      //
	      //         elevec1_epetra(idVely,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Vely) * velint_WDBC(Velx) + normal(Velx)*velint_WDBC(Vely))
	      //                                                   + derxy_(Vely,ir) *       normal(Vely) * velint_WDBC(Vely)
	      //                                                   + derxy_(Velz,ir) * 0.5* (normal(Velz) * velint_WDBC(Vely) + normal(Vely)*velint_WDBC(Velz)));
	      //
	      //         elevec1_epetra(idVelz,0) -= timefacfac_visc* (  derxy_(Velx,ir) * 0.5* (normal(Velx) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Velx))
	      //                                                   + derxy_(Vely,ir) * 0.5* (normal(Vely) * velint_WDBC(Velz) + normal(Velz)*velint_WDBC(Vely))
	      //                                                   + derxy_(Velz,ir) *       normal(Velz) * velint_WDBC(Velz));
	    }
	  }// end coupling
	  if(!bg_mortaring)
	  {
	    //-----------------------------------------------
	    //    - ((2*k2*mu2) *eps(v2)*n , u1)
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      // -(2mu*eps(v2)*n, Du1)
	      for(int ic =0; ic<nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        //(x,x)
	        C_uiu_(idVelx, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*(         normal(Velx)*emb_derxy_(Velx,ir)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
	        //(y,x)
	        C_uiu_(idVely, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velx,ir);
	        //(z,x)
	        C_uiu_(idVelz, iVelx) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Velx,ir);

	        //(x,y)
	        C_uiu_(idVelx, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Vely,ir);
	        //(y,y)
	        C_uiu_(idVely, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
	            +       normal(Vely)*emb_derxy_(Vely,ir)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
	        //(z,y)
	        C_uiu_(idVelz, iVely) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Vely,ir);

	        //(x,z)
	        C_uiu_(idVelx, iVelz) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Velz,ir);
	        //(y,z)
	        C_uiu_(idVely, iVelz) -= alpha*e_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velz,ir);
	        //(z,z)
	        C_uiu_(idVelz, Velz) -= alpha*e_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
	            +       normal(Velz)*emb_derxy_(Velz,ir)  );
	      }
	      //  (2k2mu2*eps(v2)*n, u1)
	      double timefacfac_visc = alpha*timefacfac*2.0*visceff_2*kappa2;
	      rhC_ui_(idVelx,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) *       normal(Velx) * velint(Velx)
	          + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
	          + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Velx) + normal(Velx)*velint(Velz)));

	      rhC_ui_(idVely,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Vely) * velint(Velx) + normal(Velx)*velint(Vely))
	          + emb_derxy_(Vely,ir) *       normal(Vely) * velint(Vely)
	          + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * velint(Vely) + normal(Vely)*velint(Velz)));

	      rhC_ui_(idVelz,0) += timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Velx) * velint(Velz) + normal(Velz)*velint(Velx))
	          + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * velint(Velz) + normal(Velz)*velint(Vely))
	          + emb_derxy_(Velz,ir) *       normal(Velz) * velint(Velz));
	    }

	    //-----------------------------------------------
	    //    + ((2*k2*mu2) *eps(v2)*n , u2)
	    //-----------------------------------------------
	    for(int ir = 0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      // -(2mu2*eps(v2)*n, Du2)
	      for(int ic =0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        //(x,x)
	        C_uiui_(idVelx, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*(         normal(Velx)*emb_derxy_(Velx,ir)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
	        //(y,x)
	        C_uiui_(idVely, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velx,ir);
	        //(z,x)
	        C_uiui_(idVelz, iVelx) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Velx,ir);

	        //(x,y)
	        C_uiui_(idVelx, iVely) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Vely,ir);
	        //(y,y)
	        C_uiui_(idVely, iVely) += alpha*s_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
	            +       normal(Vely)*emb_derxy_(Vely,ir)
	            + 0.5 * normal(Velz)*emb_derxy_(Velz,ir)  );
	        //(z,y)
	        C_uiui_(idVelz, iVely) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velz)*emb_derxy_(Vely,ir);

	        //(x,z)
	        C_uiui_(idVelx, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Velx)*emb_derxy_(Velz,ir);
	        //(y,z)
	        C_uiui_(idVely, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*    0.5 * normal(Vely)*emb_derxy_(Velz,ir);
	        //(z,z)
	        C_uiui_(idVelz, iVelz) += alpha*s_funct_visc2_timefacfac(ic)*(   0.5 * normal(Velx)*emb_derxy_(Velx,ir)
	            + 0.5 * normal(Vely)*emb_derxy_(Vely,ir)
	            +       normal(Velz)*emb_derxy_(Velz,ir)  );
	      }
	      //  (2mu2*eps(v2)*n, u2)
	      double timefacfac_visc = alpha*timefacfac*2.0*visceff_2*kappa2;
	      rhC_ui_(idVelx,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) *       normal(Velx) * emb_velint(Velx)
	          + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
	          + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Velx) + normal(Velx)*emb_velint(Velz)));

	      rhC_ui_(idVely,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Vely) * emb_velint(Velx) + normal(Velx)*emb_velint(Vely))
	          + emb_derxy_(Vely,ir) *       normal(Vely) * emb_velint(Vely)
	          + emb_derxy_(Velz,ir) * 0.5* (normal(Velz) * emb_velint(Vely) + normal(Vely)*emb_velint(Velz)));

	      rhC_ui_(idVelz,0) -= timefacfac_visc* (     emb_derxy_(Velx,ir) * 0.5* (normal(Velx) * emb_velint(Velz) + normal(Velz)*emb_velint(Velx))
	          + emb_derxy_(Vely,ir) * 0.5* (normal(Vely) * emb_velint(Velz) + normal(Velz)*emb_velint(Vely))
	          + emb_derxy_(Velz,ir) *       normal(Velz) * emb_velint(Velz));

	    }

	  }


	  /*                                  \        /                           i   \
	 |  gamma*mu/h_K *  [ v ] , [ Du ]     | =  - |   gamma*mu/h_K * [ v ], [ u ]   |
	  \                                   /        \                              */

/////////////////NITSCHE Penalty
	  // + gamma*mu/h_K (v1, u1))
	  for(int ir=0; ir<nen_; ir++)
	  {
	    int idVelx = ir*(nsd_+1) + 0;
	    int idVely = ir*(nsd_+1) + 1;
	    int idVelz = ir*(nsd_+1) + 2;

	    for(int ic=0; ic<nen_; ic++)
	    {
	      int iVelx = ic*(nsd_+1)+0;
	      int iVely = ic*(nsd_+1)+1;
	      int iVelz = ic*(nsd_+1)+2;

	      C_uu(idVelx, iVelx) += funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	      C_uu(idVely, iVely) += funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	      C_uu(idVelz, iVelz) += funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	    }

	    // -(stab * v, u)
	    rhs_Cu(idVelx,0) -= funct_timefacfac(ir)*NIT_full_stab_fac*velint(Velx);
	    rhs_Cu(idVely,0) -= funct_timefacfac(ir)*NIT_full_stab_fac*velint(Vely);
	    rhs_Cu(idVelz,0) -= funct_timefacfac(ir)*NIT_full_stab_fac*velint(Velz);

	    if(!coupling) // weak Dirichlet case
	    {
	      // +(stab * v, u_DBC)
	      rhs_Cu(idVelx,0) += funct_timefacfac(ir)*NIT_full_stab_fac*ivelint_WDBC_JUMP(Velx);
	      rhs_Cu(idVely,0) += funct_timefacfac(ir)*NIT_full_stab_fac*ivelint_WDBC_JUMP(Vely);
	      rhs_Cu(idVelz,0) += funct_timefacfac(ir)*NIT_full_stab_fac*ivelint_WDBC_JUMP(Velz);
	    }

	  }

	  if(coupling)
	  {
	    // - gamma*mu/h_K (v1, u2))
	    for(int ir=0; ir<nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      for(int ic=0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uui_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ic,ir)*NIT_full_stab_fac;
	        C_uui_(idVely, iVely) -= emb_funct_dyad_timefacfac(ic,ir)*NIT_full_stab_fac;
	        C_uui_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ic,ir)*NIT_full_stab_fac;
	      }

	      // -(stab * v, u)
	      rhs_Cu(idVelx,0) += funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Velx);
	      rhs_Cu(idVely,0) += funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Vely);
	      rhs_Cu(idVelz,0) += funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Velz);


	    }

	    // - gamma*mu/h_K (v2, u1))
	    for(int ir=0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      for(int ic=0; ic<nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uiu_(idVelx, iVelx) -= emb_funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	        C_uiu_(idVely, iVely) -= emb_funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	        C_uiu_(idVelz, iVelz) -= emb_funct_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	      }

	      // +(stab * v2, u1)
	      rhC_ui_(idVelx,0) += emb_funct_timefacfac(ir)*NIT_full_stab_fac*velint(Velx);
	      rhC_ui_(idVely,0) += emb_funct_timefacfac(ir)*NIT_full_stab_fac*velint(Vely);
	      rhC_ui_(idVelz,0) += emb_funct_timefacfac(ir)*NIT_full_stab_fac*velint(Velz);

	    }

	    //       LINALG::Matrix<emb_nen_,emb_nen_> emb_emb_dyad_timefacfac(true);
	    //       emb_emb_dyad_timefacfac.MultiplyNT(emb_funct_timefacfac, emb_funct_);
	    // + gamma*mu/h_K (v2, u2))
	    for(int ir=0; ir<emb_nen_; ir++)
	    {
	      int idVelx = ir*(nsd_+1) + 0;
	      int idVely = ir*(nsd_+1) + 1;
	      int idVelz = ir*(nsd_+1) + 2;

	      for(int ic=0; ic<emb_nen_; ic++)
	      {
	        int iVelx = ic*(nsd_+1)+0;
	        int iVely = ic*(nsd_+1)+1;
	        int iVelz = ic*(nsd_+1)+2;

	        C_uiui_(idVelx, iVelx) += emb_emb_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	        C_uiui_(idVely, iVely) += emb_emb_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	        C_uiui_(idVelz, iVelz) += emb_emb_dyad_timefacfac(ir,ic)*NIT_full_stab_fac;
	      }

	      // -(stab * v2, u2)
	      rhC_ui_(idVelx,0) -= emb_funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Velx);
	      rhC_ui_(idVely,0) -= emb_funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Vely);
	      rhC_ui_(idVelz,0) -= emb_funct_timefacfac(ir)*NIT_full_stab_fac*emb_velint(Velz);

	    }
	  } // end coupling




	  //-----------------------------------------------------------------
	   // penalty term for velocity gradients at the interface
	   if (velgrad_interface_stab)
	   {
	     double tau_timefacfac = velgrad_interface_fac * timefacfac;

	     // get grad(u)*n
	     LINALG::Matrix<emb_nen_,1> emb_normal_deriv(true);
	     emb_normal_deriv.MultiplyTN(emb_derxy_,normal);

	     LINALG::Matrix<nen_,1> e_normal_deriv(true);
	     e_normal_deriv.MultiplyTN(derxy,normal);

	     // additional stability of gradients
	     // parent column
	     for (int ui=0; ui<emb_nen_; ++ui)
	     {
	       for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
	       {
	         int col = ui*4+ijdim;

	         // v_parent * u_parent
	         //parent row
	         for (int vi=0; vi<emb_nen_; ++vi)
	         {
	           C_uiui_(vi*4+ijdim, col) += tau_timefacfac*emb_normal_deriv(vi)*emb_normal_deriv(ui);
	         }

	        // neighbor row
	         for (int vi=0; vi<nen_; ++vi)
	         {
	           C_uui_(vi*4+ijdim, col) -= tau_timefacfac*e_normal_deriv(vi)*emb_normal_deriv(ui);
	         }
	       }
	     }

	     for (int ui=0; ui<nen_; ++ui)
	     {
	       for (int ijdim = 0; ijdim <3; ++ijdim) // combined components of u and v
	       {
	        int col = ui*4+ijdim;

	         // v_parent * u_parent
	         //parent row
	         for (int vi=0; vi<emb_nen_; ++vi)
	         {
	           C_uiu_(vi*4+ijdim, col) -= tau_timefacfac*emb_normal_deriv(vi)*e_normal_deriv(ui);
	         }

	         //neighbor row
	         for (int vi=0; vi<nen_; ++vi)
	         {
	           C_uu(vi*4+ijdim, col) += tau_timefacfac*e_normal_deriv(vi)*e_normal_deriv(ui);
	         }
	       }
	     }

	     LINALG::Matrix<3,1> emb_grad_u_n(true);
	     emb_grad_u_n.Multiply(emb_vderxy_,normal);
	     LINALG::Matrix<3,1> e_grad_u_n(true);
	     e_grad_u_n.Multiply(vderxy,normal);

	     for(int idim = 0; idim <3; ++idim)
	     {
	       double diff_grad_u_n = tau_timefacfac * (e_grad_u_n(idim)-emb_grad_u_n(idim));

	       // v_parent (u_neighbor-u_parent)
	       for (int vi=0; vi<emb_nen_; ++vi)
	       {
	         rhC_ui_(vi*4+idim,0) +=  emb_normal_deriv(vi)*diff_grad_u_n;
	       }

	       // v_neighbor (u_neighbor-u_parent)
	       for (int vi=0; vi<nen_; ++vi)
	       {
	         rhs_Cu(vi*4+idim,0) -=  e_normal_deriv(vi)*diff_grad_u_n;
	       }
	     }
	   }// velgrad_interface_fac
	   // --------------------------------------------------------
	    // pressure coupling penalty term at the interface
	    if (presscoupling_interface_stab)
	    {
	      double tau_timefacfac_press = presscoupling_interface_fac * timefacfac;
	      for (int ui=0; ui<emb_nen_; ++ui)
	      {
	        int col = ui*4+3;
	        // v_parent * u_parent
	        //parent row
	        for (int vi=0; vi<emb_nen_; ++vi)
	        {
	          C_uiui_(vi*4+3, col) += tau_timefacfac_press*emb_funct_(vi)*emb_funct_(ui);
	        }

	        // neighbor row
	        for (int vi=0; vi<nen_; ++vi)
	        {
	          C_uui_(vi*4+3, col) -= tau_timefacfac_press*funct(vi)*emb_funct_(ui);
	        }
	      }

	      for (int ui=0; ui<nen_; ++ui)
	      {
	         int col = ui*4+3;
	         // v_parent * u_parent
	         //parent row
	         for (int vi=0; vi<emb_nen_; ++vi)
	         {
	           C_uiu_(vi*4+3, col) -= tau_timefacfac_press*emb_funct_(vi)*funct(ui);
	         }

	         //neighbor row
	         for (int vi=0; vi<nen_; ++vi)
	         {
	           C_uu(vi*4+3, col) += tau_timefacfac_press*funct(vi)*funct(ui);
	         }
	      }

	      double press_jump = press - emb_press;

	      for (int vi=0; vi<emb_nen_; ++vi)
	      {
	        int row = vi*4+3;
	        rhC_ui_(row,0) += emb_funct_(vi)*press_jump*tau_timefacfac_press;
	      }

	      for (int vi=0; vi<nen_; ++vi)
	      {
	        int row = vi*4+3;
	        rhs_Cu(row,0) -= funct(vi)*press_jump*tau_timefacfac_press;
	      }

	    } // presscoupling_interface_stab

//   if(coupling)
	 if(  xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty or
        xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty)
   {

   NIT2_Stab_ConvAveraged(C_uu,          // standard bg-bg-matrix
	       rhs_Cu,        // standard bg-rhs
	       velint,
	       emb_velint,
	       ivelint_WDBC_JUMP,
	       funct_timefacfac,
	       funct_dyad_timefacfac,
	       emb_funct_timefacfac,
	       emb_funct_dyad_timefacfac,
	       avg_conv_stab_fac,
	       coupling,       // assemble coupling terms (yes/no)
	       normal         // normal vector
	   );

   }

   return;
}

/*--------------------------------------------------------------------------------
 * evaluate stabilizing viscous term for Nitsche's method
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, DRT::Element::DiscretizationType emb_distype>
void EmbImpl<distype, emb_distype>::NIT2_Stab_ConvAveraged(
    Epetra_SerialDenseMatrix &               C_uu_,                            ///< standard bg-bg-matrix
    Epetra_SerialDenseVector &               rhs_Cu_,                          ///< standard bg-rhs
    const LINALG::Matrix<nsd_,1>&            velint,                           ///< velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint,                          ///< interface velocity at integration point
    const LINALG::Matrix<nsd_,1>&            ivelint_WDBC_JUMP,                ///< prescribed interface velocity or jump vector
    const LINALG::Matrix<nen_,1>&            funct_timefacfac,                 ///< funct * timefacfac
    const LINALG::Matrix<nen_,nen_>&         funct_dyad_timefacfac,            ///< (funct^T * funct) * timefacfac
    const LINALG::Matrix<emb_nen_,1>&       emb_funct_timefacfac,            ///< sidefunct^T * timefacfac
    const LINALG::Matrix<emb_nen_,nen_>&    emb_funct_dyad_timefacfac,       ///< (sidefunct^T * funct) * timefacfac
    const double &                           stabfac_avg,                          ///< stabilization factor
    const bool &                             coupling,                         ///< assemble coupling terms (yes/no)
    const LINALG::Matrix<nsd_,1> &           normal                            ///< normal vector
)
{
	  const unsigned Velx = 0;
	  const unsigned Vely = 1;
	  const unsigned Velz = 2;

/*                                        \        /                                       i \
|  [rho * (beta * n^e)] *  { v }_m , [ Du ] | =  - |  [rho * (beta * n^e)] * { v }_m,  [ u ]   |
\ ----stab_avg-----                      /         \ ----stab_avg-----                    */


  //  [rho * (beta * n^e)] (k1*vb,ub)
  for(int ir=0; ir<nen_; ir++)
  {
    int idVelx = ir*(nsd_+1) + 0;
    int idVely = ir*(nsd_+1) + 1;
    int idVelz = ir*(nsd_+1) + 2;

    for(int ic=0; ic<nen_; ic++)
    {
      int iVelx = ic*(nsd_+1)+0;
      int iVely = ic*(nsd_+1)+1;
      int iVelz = ic*(nsd_+1)+2;

      //minus?
      double tmp = 0.5*funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uu_(idVelx, iVelx) += tmp;
      C_uu_(idVely, iVely) += tmp;
      C_uu_(idVelz, iVelz) += tmp;
    }

    double tmp = 0.5*funct_timefacfac(ir)*stabfac_avg;

    rhs_Cu_(idVelx,0) -= tmp*velint(Velx);
    rhs_Cu_(idVely,0) -= tmp*velint(Vely);
    rhs_Cu_(idVelz,0) -= tmp*velint(Velz);
  }

  //  -[rho * (beta * n^e)] (k1*vb,ue)
  for(int ir=0; ir<nen_; ir++)
  {
    int idVelx = ir*(nsd_+1) + 0;
    int idVely = ir*(nsd_+1) + 1;
    int idVelz = ir*(nsd_+1) + 2;

    for(int ic=0; ic<emb_nen_; ic++)
    {
      int iVelx = ic*(nsd_+1)+0;
      int iVely = ic*(nsd_+1)+1;
      int iVelz = ic*(nsd_+1)+2;

      double tmp = -0.5*emb_funct_dyad_timefacfac(ic,ir)*stabfac_avg;

      C_uui_(idVelx, iVelx) += tmp;
      C_uui_(idVely, iVely) += tmp;
      C_uui_(idVelz, iVelz) += tmp;
    }

    double tmp = -0.5*funct_timefacfac(ir)*stabfac_avg;

    rhs_Cu_(idVelx,0) -= tmp*ivelint(Velx);
    rhs_Cu_(idVely,0) -= tmp*ivelint(Vely);
    rhs_Cu_(idVelz,0) -= tmp*ivelint(Velz);
  }

  //  [rho * (beta * n^e)] (k1*ve,ub)
  for(int ir=0; ir<emb_nen_; ir++)
  {
    int idVelx = ir*(nsd_+1) + 0;
    int idVely = ir*(nsd_+1) + 1;
    int idVelz = ir*(nsd_+1) + 2;

    for(int ic=0; ic<nen_; ic++)
    {
      int iVelx = ic*(nsd_+1)+0;
      int iVely = ic*(nsd_+1)+1;
      int iVelz = ic*(nsd_+1)+2;

      double tmp = 0.5*emb_funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uiu_(idVelx, iVelx) += tmp;
      C_uiu_(idVely, iVely) += tmp;
      C_uiu_(idVelz, iVelz) += tmp;
    }

    double tmp = 0.5*emb_funct_timefacfac(ir)*stabfac_avg;

    rhC_ui_(idVelx,0) -= tmp*velint(Velx);
    rhC_ui_(idVely,0) -= tmp*velint(Vely);
    rhC_ui_(idVelz,0) -= tmp*velint(Velz);
  }

  //-[rho * (beta * n^e)] (k1*ve,ue)
  for(int ir=0; ir<emb_nen_; ir++)
  {
    int idVelx = ir*(nsd_+1) + 0;
    int idVely = ir*(nsd_+1) + 1;
    int idVelz = ir*(nsd_+1) + 2;

    for(int ic=0; ic<emb_nen_; ic++)
    {
      int iVelx = ic*(nsd_+1)+0;
      int iVely = ic*(nsd_+1)+1;
      int iVelz = ic*(nsd_+1)+2;

      double tmp = -0.5*emb_funct_dyad_timefacfac(ir,ic)*stabfac_avg;

      C_uiui_(idVelx, iVelx) += tmp;
      C_uiui_(idVely, iVely) += tmp;
      C_uiui_(idVelz, iVelz) += tmp;
    }

    double tmp = -0.5*emb_funct_timefacfac(ir)*stabfac_avg;

    rhC_ui_(idVelx,0) -= tmp*ivelint(Velx);
    rhC_ui_(idVely,0) -= tmp*ivelint(Vely);
    rhC_ui_(idVelz,0) -= tmp*ivelint(Velz);
  }
  return;
}

} // namespace XFLUID
} // namespace ELEMENTS
} // namespace DRT

// pairs with numdof=3
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex8,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex8,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex20, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex20, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex27, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex27, DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet4,  DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet4,  DRT::Element::quad8,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet10, DRT::Element::quad4,3>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet10, DRT::Element::quad8,3>;


// pairs with numdof=4
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex8,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex8,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex20, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex20, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex27, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::hex27, DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet4,  DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet4,  DRT::Element::quad8,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet10, DRT::Element::quad4,4>;
template class DRT::ELEMENTS::XFLUID::SideImpl<DRT::Element::tet10, DRT::Element::quad8,4>;



template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex8,  DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex8,  DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex20, DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex20, DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex27, DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::hex27, DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::tet4,  DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::tet4,  DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::tet10, DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbImpl<DRT::Element::tet10, DRT::Element::hex20>;





