/*----------------------------------------------------------------------*/
/*!
\file xfluid3_impl.cpp

\brief Internal implementation of XFluid3 element

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3_impl.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"

#include <Epetra_SerialDenseSolver.h>


DRT::Elements::XFluid3Impl::XFluid3Impl(
		const int iel,
		const int numparamvelx,
		const int numparamvely,
		const int numparamvelz,
		const int numparampres)
  : iel_(iel),
    vart_(),
    xyze_(3,iel_,blitz::ColumnMajorArray<2>()),
    edeadng_(3,numparamvelx,blitz::ColumnMajorArray<2>()),
    funct_(iel_),
    deriv_(3,iel_,blitz::ColumnMajorArray<2>()),
    deriv2_(6,iel_,blitz::ColumnMajorArray<2>()),
    xjm_(3,3,blitz::ColumnMajorArray<2>()),
    xji_(3,3,blitz::ColumnMajorArray<2>()),
    vderxy_(3,3,blitz::ColumnMajorArray<2>()),
    pderxy_(3),
    vderxy2_(3,6,blitz::ColumnMajorArray<2>()),
    derxy_(3,iel_,blitz::ColumnMajorArray<2>()),
    derxy2_(6,iel_,blitz::ColumnMajorArray<2>()),
    bodyforce_(3),
    histvec_(3),
    velino_(3),
    velint_(3),
    gridvelint_(3),
    gradp_(3),
    tau_(3),
    viscs2_(3,3,numparamvelx,blitz::ColumnMajorArray<3>()),
    conv_c_(numparamvelx),
    conv_g_(numparamvelx),
    conv_r_(3,3,numparamvelx,blitz::ColumnMajorArray<3>()),
    rhsint_(3),
    conv_old_(3),
    visc_old_(3),
    res_old_(3),
    xder2_(6,3,blitz::ColumnMajorArray<2>())
{
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                g.bau 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Impl::Sysmat(XFluid3* ele,
									   const XFEM::DomainIntCells&   domainIntCells,
									   const XFEM::BoundaryIntCells& boundaryIntCells,
                                       const blitz::Array<double,2>&     evelnp,
                                       const blitz::Array<double,1>&     eprenp,
                                       const blitz::Array<double,2>&     evhist,
                                       const blitz::Array<double,2>&     edispnp,
                                       const blitz::Array<double,2>&     egridv,
                                       blitz::Array<double,2>&           estif,
                                       blitz::Array<double,2>&           esv,
                                       blitz::Array<double,1>&           eforce,
                                       struct _MATERIAL*       material,
                                       double                  time,
                                       double                  timefac,
                                       bool                    newton ,
                                       int                     fssgv ,
                                       bool                    pstab  ,
                                       bool                    supg   ,
                                       bool                    vstab  ,
                                       bool                    cstab
  )
{

  // set element data
  const DRT::Element::DiscretizationType distype = ele->Shape();

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }

  // add displacement, when fluid nodes move in the ALE case
  if (ele->is_ale_)
  {
    xyze_ += edispnp;
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
//   blitz::fourthIndex l;   // Placeholder for the fourth index

  blitz::Range _  = blitz::Range::all();
//   blitz::Range ux = blitz::Range(0, 4*iel_-4, 4);
//   blitz::Range uy = blitz::Range(1, 4*iel_-3, 4);
//   blitz::Range uz = blitz::Range(2, 4*iel_-2, 4);
//   blitz::Range p  = blitz::Range(3, 4*iel_-1, 4);

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  Caltau(ele,evelnp,distype,visc,timefac,fssgv);

  // flag for higher order elements
  const bool higher_order_ele = ele->isHigherOrderElement(distype);

  // integrate of integration cell
  dsassert(domainIntCells.empty() == false, "this is a bug!");
  
  DRT::Utils::GaussRule3D gaussrule;
  bool standard_integration = true;
  if (domainIntCells.size() == 1){
      gaussrule = ele->gaussrule_;
      standard_integration = true;
  }
  else {
      gaussrule = DRT::Utils::intrule_tet_10point;
      standard_integration = false;
  }
  
  for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell < domainIntCells.end(); ++cell)
  {
  
  // gaussian points
  const DRT::Utils::IntegrationPoints3D intpoints(gaussrule);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point in cell coordinates \eta
    const double cell_e0 = intpoints.qxg[iquad][0];
    const double cell_e1 = intpoints.qxg[iquad][1];
    const double cell_e2 = intpoints.qxg[iquad][2];
    
    const vector<double> e = cell->modifyGaussRule3D(standard_integration,cell_e0,cell_e1,cell_e2);

    // coordinates of the current integration point in element coordinates \xi
    const double e1 = e[0];
    const double e2 = e[1];
    const double e3 = e[2];
    const double detcell = e[3];

    // shape functions and their derivatives
    DRT::Utils::shape_function_3D(funct_,e1,e2,e3,distype);
    DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

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
    xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                       xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                       xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                       xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                       xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                       xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
    const double fac = intpoints.qwgt[iquad]*det*detcell;

    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %lf", ele->Id(), det);
    }

    // inverse of jacobian
    xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
    xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
    xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
    xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
    xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
    xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
    xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
    xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
    xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;

    // compute global derivates
    derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

    // compute second global derivative
    if (higher_order_ele)
    {
      DRT::Utils::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      gder2(ele);

      // calculate 2nd velocity derivatives at integration point
      vderxy2_ = blitz::sum(derxy2_(j,k)*evelnp(i,k),k);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get velocities (n+g,i) at integration point
    velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

    // get history data (n,i) at integration point
    histvec_ = blitz::sum(funct_(j)*evhist(i,j),j);

    // get velocity (np,i) derivatives at integration point
    vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

    // get grid velocity at integration point
    if (ele->is_ale_)
    {
      gridvelint_ = blitz::sum(funct_(j)*egridv(i,j),j);
    }
    else
    {
      gridvelint_ = 0.;
    }

    // get pressure gradients
    gradp_ = blitz::sum(derxy_(i,j)*eprenp(j),j);

    double press = blitz::sum(funct_*eprenp);

    // get bodyforce in gausspoint
    bodyforce_ = blitz::sum(edeadng_(i,j)*funct_(j),j);

    // perform integration for entire matrix and rhs

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    // integration factors and coefficients of single terms
    const double timetauM   = timefac * tau_M;
    const double timetauMp  = timefac * tau_Mp;

    const double ttimetauM  = timefac * timetauM;
    const double ttimetauMp = timefac * timetauMp;
    const double timefacfac = timefac * fac;

    const double vartfac = vart_*timefacfac;

    /*------------------------- evaluate rhs vector at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsint_ = histvec_(i) + bodyforce_(i)*timefac;

    /*----------------- get numerical representation of single operators ---*/

    /* Convective term  u_old * grad u_old: */
    conv_old_ = blitz::sum(vderxy_(i, j)*velint_(j), j);

    /* Viscous term  div epsilon(u_old) */
    visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,3) + vderxy2_(0,2) + vderxy2_(2,4));
    visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,3) + vderxy2_(1,2) + vderxy2_(2,5));
    visc_old_(2) = vderxy2_(2,2) + 0.5 * (vderxy2_(2,0) + vderxy2_(0,4) + vderxy2_(2,1) + vderxy2_(1,5));

    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_c_ = blitz::sum(derxy_(j,i)*velint_(j), j);

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    if (ele->is_ale_)
    {
      conv_g_ = - blitz::sum(derxy_(j,i) * gridvelint_(j), j);
    }
    else
    {
      conv_g_ = 0.0;
    }

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                                     \
       |  u_old_x,x   u_old_x,y   u_old x,z  |
       |                                     |
       |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
       |                                     |
       |  u_old_z,x   u_old_z,y   u_old_z,z  |
       \                                     /
       with  N .. form function matrix                                   */
    conv_r_ = vderxy_(i, j)*funct_(k);

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                                                \
         |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
       1 |                                                |
     - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
       2 |                                                |
         |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
         \                                                /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */

    viscs2_(0,0,_) = - 0.5 * (2.0 * derxy2_(0,_) + derxy2_(1,_) + derxy2_(2,_));
    viscs2_(0,1,_) = - 0.5 *  derxy2_(3,_);
    viscs2_(0,2,_) = - 0.5 *  derxy2_(4,_);
    viscs2_(1,0,_) = - 0.5 *  derxy2_(3,_);
    viscs2_(1,1,_) = - 0.5 * (derxy2_(0,_) + 2.0 * derxy2_(1,_) + derxy2_(2,_));
    viscs2_(1,2,_) = - 0.5 *  derxy2_(5,_);
    viscs2_(2,0,_) = - 0.5 *  derxy2_(4,_);
    viscs2_(2,1,_) = - 0.5 *  derxy2_(5,_);
    viscs2_(2,2,_) = - 0.5 * (derxy2_(0,_) + derxy2_(1,_) + 2.0 * derxy2_(2,_));

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--------------------------------- now build single stiffness terms ---*/
#if 0
    if (ele->is_ale_)
    {
      for (int vi=0; vi<iel_; ++vi)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          /* Konvektionsterm */
          estif(vi*4, ui*4)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
          estif(vi*4, ui*4 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
          estif(vi*4, ui*4 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
          estif(vi*4 + 1, ui*4)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
          estif(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
          estif(vi*4 + 1, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
          estif(vi*4 + 2, ui*4)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
          estif(vi*4 + 2, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
          estif(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;

          /* Stabilisierung der Konvektion ( L_conv_u) */
          estif(vi*4, ui*4)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + conv_g_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxy_(0, vi)*conv_r_(0, 2, ui) - gridvelint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) - gridvelint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) - gridvelint_(2)*derxy_(0, vi)*conv_r_(0, 2, ui)) ;
          estif(vi*4, ui*4 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + conv_g_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxy_(1, vi)*conv_r_(0, 2, ui) - gridvelint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) - gridvelint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) - gridvelint_(2)*derxy_(1, vi)*conv_r_(0, 2, ui)) ;
          estif(vi*4, ui*4 + 2)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + conv_g_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxy_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxy_(2, vi)*conv_r_(0, 2, ui) - gridvelint_(0)*derxy_(2, vi)*conv_r_(0, 0, ui) - gridvelint_(1)*derxy_(2, vi)*conv_r_(0, 1, ui) - gridvelint_(2)*derxy_(2, vi)*conv_r_(0, 2, ui)) ;
          estif(vi*4 + 1, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + conv_g_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxy_(0, vi)*conv_r_(1, 2, ui) - gridvelint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) - gridvelint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) - gridvelint_(2)*derxy_(0, vi)*conv_r_(1, 2, ui)) ;
          estif(vi*4 + 1, ui*4 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + conv_g_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxy_(1, vi)*conv_r_(1, 2, ui) - gridvelint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) - gridvelint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) - gridvelint_(2)*derxy_(1, vi)*conv_r_(1, 2, ui)) ;
          estif(vi*4 + 1, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + conv_g_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxy_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxy_(2, vi)*conv_r_(1, 2, ui) - gridvelint_(0)*derxy_(2, vi)*conv_r_(1, 0, ui) - gridvelint_(1)*derxy_(2, vi)*conv_r_(1, 1, ui) - gridvelint_(2)*derxy_(2, vi)*conv_r_(1, 2, ui)) ;
          estif(vi*4 + 2, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + conv_g_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxy_(0, vi)*conv_r_(2, 2, ui) - gridvelint_(0)*derxy_(0, vi)*conv_r_(2, 0, ui) - gridvelint_(1)*derxy_(0, vi)*conv_r_(2, 1, ui) - gridvelint_(2)*derxy_(0, vi)*conv_r_(2, 2, ui)) ;
          estif(vi*4 + 2, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + conv_g_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxy_(1, vi)*conv_r_(2, 2, ui) - gridvelint_(0)*derxy_(1, vi)*conv_r_(2, 0, ui) - gridvelint_(1)*derxy_(1, vi)*conv_r_(2, 1, ui) - gridvelint_(2)*derxy_(1, vi)*conv_r_(2, 2, ui)) ;
          estif(vi*4 + 2, ui*4 + 2) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + conv_g_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxy_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxy_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxy_(2, vi)*conv_r_(2, 2, ui) - gridvelint_(0)*derxy_(2, vi)*conv_r_(2, 0, ui) - gridvelint_(1)*derxy_(2, vi)*conv_r_(2, 1, ui) - gridvelint_(2)*derxy_(2, vi)*conv_r_(2, 2, ui)) ;

          /* Stabilisierung der Konvektion (-L_visc_u) */
          estif(vi*4, ui*4)         += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxy_(0, vi)) ;
          estif(vi*4, ui*4 + 1)     += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxy_(1, vi)) ;
          estif(vi*4, ui*4 + 2)     += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4)     += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxy_(0, vi)) ;
          estif(vi*4 + 1, ui*4 + 1) += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxy_(1, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxy_(2, vi)) ;
          estif(vi*4 + 2, ui*4)     += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxy_(0, vi)) ;
          estif(vi*4 + 2, ui*4 + 1) += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxy_(1, vi)) ;
          estif(vi*4 + 2, ui*4 + 2) += -2.0*visc*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxy_(2, vi)) ;

          /* Stabilisierung der Konvektion ( L_pres_p) */
          estif(vi*4, ui*4)         += ttimetauM*funct_(ui)*gradp_(0)*derxy_(0, vi) ;
          estif(vi*4, ui*4 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxy_(1, vi) ;
          estif(vi*4, ui*4 + 2)     += ttimetauM*funct_(ui)*gradp_(0)*derxy_(2, vi) ;
          estif(vi*4, ui*4 + 3)     += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxy_(0, ui) ;
          estif(vi*4 + 1, ui*4)     += ttimetauM*funct_(ui)*gradp_(1)*derxy_(0, vi) ;
          estif(vi*4 + 1, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxy_(1, vi) ;
          estif(vi*4 + 1, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(1)*derxy_(2, vi) ;
          estif(vi*4 + 1, ui*4 + 3) += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxy_(1, ui) ;
          estif(vi*4 + 2, ui*4)     += ttimetauM*funct_(ui)*gradp_(2)*derxy_(0, vi) ;
          estif(vi*4 + 2, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(2)*derxy_(1, vi) ;
          estif(vi*4 + 2, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(2)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4 + 3) += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxy_(2, ui) ;

          /* Viskosit�tsterm */
          estif(vi*4, ui*4)         += visc*timefacfac*(2.0*derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4, ui*4 + 1)     += visc*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*4, ui*4 + 2)     += visc*timefacfac*derxy_(0, ui)*derxy_(2, vi) ;
          estif(vi*4 + 1, ui*4)     += visc*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*4 + 1, ui*4 + 1) += visc*timefacfac*(derxy_(0, ui)*derxy_(0, vi) + 2.0*derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += visc*timefacfac*derxy_(1, ui)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4)     += visc*timefacfac*derxy_(0, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 1) += visc*timefacfac*derxy_(1, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 2) += visc*timefacfac*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + 2.0*derxy_(2, ui)*derxy_(2, vi)) ;

          /* Stabilisierung der Viskosit�t ( L_conv_u) */
          estif(vi*4, ui*4)         += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
          estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
          estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
          estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
          estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
          estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
          estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
          estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
          estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;

          /* Stabilisierung der Viskosit�t (-L_visc_u) */
          estif(vi*4, ui*4)         += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
          estif(vi*4, ui*4 + 1)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
          estif(vi*4, ui*4 + 2)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
          estif(vi*4 + 1, ui*4)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
          estif(vi*4 + 1, ui*4 + 1) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
          estif(vi*4 + 2, ui*4)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
          estif(vi*4 + 2, ui*4 + 1) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
          estif(vi*4 + 2, ui*4 + 2) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

          /* Stabilisierung der Viskosit�t ( L_pres_p) */
          estif(vi*4, ui*4 + 3)     += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 0, vi) + derxy_(1, ui)*viscs2_(0, 1, vi) + derxy_(2, ui)*viscs2_(0, 2, vi)) ;
          estif(vi*4 + 1, ui*4 + 3) += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 1, vi) + derxy_(1, ui)*viscs2_(1, 1, vi) + derxy_(2, ui)*viscs2_(1, 2, vi)) ;
          estif(vi*4 + 2, ui*4 + 3) += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 2, vi) + derxy_(1, ui)*viscs2_(1, 2, vi) + derxy_(2, ui)*viscs2_(2, 2, vi)) ;

          /* Druckterm */
          estif(vi*4, ui*4 + 3)     += -(timefacfac*funct_(ui)*derxy_(0, vi)) ;
          estif(vi*4 + 1, ui*4 + 3) += -(timefacfac*funct_(ui)*derxy_(1, vi)) ;
          estif(vi*4 + 2, ui*4 + 3) += -(timefacfac*funct_(ui)*derxy_(2, vi)) ;

          /* Stabilisierung des Drucks ( L_conv_u) */
          estif(vi*4 + 3, ui*4)     += ttimetauMp*(conv_c_(ui)*derxy_(0, vi) + conv_g_(ui)*derxy_(0, vi) + derxy_(0, vi)*conv_r_(0, 0, ui) + derxy_(1, vi)*conv_r_(1, 0, ui) + derxy_(2, vi)*conv_r_(2, 0, ui)) ;
          estif(vi*4 + 3, ui*4 + 1) += ttimetauMp*(conv_c_(ui)*derxy_(1, vi) + conv_g_(ui)*derxy_(1, vi) + derxy_(0, vi)*conv_r_(0, 1, ui) + derxy_(1, vi)*conv_r_(1, 1, ui) + derxy_(2, vi)*conv_r_(2, 1, ui)) ;
          estif(vi*4 + 3, ui*4 + 2) += ttimetauMp*(conv_c_(ui)*derxy_(2, vi) + conv_g_(ui)*derxy_(2, vi) + derxy_(0, vi)*conv_r_(0, 2, ui) + derxy_(1, vi)*conv_r_(1, 2, ui) + derxy_(2, vi)*conv_r_(2, 2, ui)) ;

          /* Stabilisierung des Drucks (-L_visc_u) */
          estif(vi*4 + 3, ui*4)     += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 0, ui) + derxy_(1, vi)*viscs2_(0, 1, ui) + derxy_(2, vi)*viscs2_(0, 2, ui)) ;
          estif(vi*4 + 3, ui*4 + 1) += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 1, ui) + derxy_(1, vi)*viscs2_(1, 1, ui) + derxy_(2, vi)*viscs2_(1, 2, ui)) ;
          estif(vi*4 + 3, ui*4 + 2) += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 2, ui) + derxy_(1, vi)*viscs2_(1, 2, ui) + derxy_(2, vi)*viscs2_(2, 2, ui)) ;

          /* Stabilisierung des Drucks ( L_pres_p) */
          estif(vi*4 + 3, ui*4 + 3) += ttimetauMp*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi)) ;

          /* Divergenzfreiheit */
          estif(vi*4 + 3, ui*4)     += timefacfac*funct_(vi)*derxy_(0, ui) ;
          estif(vi*4 + 3, ui*4 + 1) += timefacfac*funct_(vi)*derxy_(1, ui) ;
          estif(vi*4 + 3, ui*4 + 2) += timefacfac*funct_(vi)*derxy_(2, ui) ;

          /* Kontinuit�tsstabilisierung */
          estif(vi*4, ui*4)         += (timefac*timefac)*tau_C*derxy_(0, ui)*derxy_(0, vi) ;
          estif(vi*4, ui*4 + 1)     += (timefac*timefac)*tau_C*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*4, ui*4 + 2)     += (timefac*timefac)*tau_C*derxy_(0, vi)*derxy_(2, ui) ;
          estif(vi*4 + 1, ui*4)     += (timefac*timefac)*tau_C*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*4 + 1, ui*4 + 1) += (timefac*timefac)*tau_C*derxy_(1, ui)*derxy_(1, vi) ;
          estif(vi*4 + 1, ui*4 + 2) += (timefac*timefac)*tau_C*derxy_(1, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4)     += (timefac*timefac)*tau_C*derxy_(0, ui)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4 + 1) += (timefac*timefac)*tau_C*derxy_(1, ui)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4 + 2) += (timefac*timefac)*tau_C*derxy_(2, ui)*derxy_(2, vi) ;

          /* Massenterm */
          estif(vi*4, ui*4)         += fac*funct_(ui)*funct_(vi) ;
          estif(vi*4 + 1, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
          estif(vi*4 + 2, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;

          /* Konvektionsstabilisierung */
          estif(vi*4, ui*4)         += timetauM*funct_(ui)*(conv_g_(vi) + 2.0*velint_(0)*derxy_(0, vi) + velint_(1)*derxy_(1, vi) + velint_(2)*derxy_(2, vi)) ;
          estif(vi*4, ui*4 + 1)     += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
          estif(vi*4, ui*4 + 2)     += timetauM*funct_(ui)*velint_(0)*derxy_(2, vi) ;
          estif(vi*4 + 1, ui*4)     += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
          estif(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxy_(0, vi) + 2.0*velint_(1)*derxy_(1, vi) + velint_(2)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += timetauM*funct_(ui)*velint_(1)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4)     += timetauM*funct_(ui)*velint_(2)*derxy_(0, vi) ;
          estif(vi*4 + 2, ui*4 + 1) += timetauM*funct_(ui)*velint_(2)*derxy_(1, vi) ;
          estif(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxy_(0, vi) + velint_(1)*derxy_(1, vi) + 2.0*velint_(2)*derxy_(2, vi)) ;

          /* Viskosit�tsstabilisierung */
          estif(vi*4, ui*4)         += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
          estif(vi*4, ui*4 + 1)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
          estif(vi*4, ui*4 + 2)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
          estif(vi*4 + 1, ui*4)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
          estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
          estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
          estif(vi*4 + 2, ui*4)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
          estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
          estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;

          /* Stabilisierung der Druckgleichung */
          estif(vi*4 + 3, ui*4)     += timetauMp*funct_(ui)*derxy_(0, vi) ;
          estif(vi*4 + 3, ui*4 + 1) += timetauMp*funct_(ui)*derxy_(1, vi) ;
          estif(vi*4 + 3, ui*4 + 2) += timetauMp*funct_(ui)*derxy_(2, vi) ;

          /* Quellterm der rechten Seite */

          /* Konvektionsstabilisierung */
          estif(vi*4, ui*4)         += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
          estif(vi*4, ui*4 + 1)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
          estif(vi*4, ui*4 + 2)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4)     += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
          estif(vi*4 + 1, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(2, vi)) ;
          estif(vi*4 + 2, ui*4)     += -(timetauM*funct_(ui)*rhsint_(2)*derxy_(0, vi)) ;
          estif(vi*4 + 2, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(2)*derxy_(1, vi)) ;
          estif(vi*4 + 2, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(2)*derxy_(2, vi)) ;

          /* Viskosit�tsstabilisierung */

          /* Stabilisierung der Druckgleichung */

        }
      }
      for (int vi=0; vi<iel_; ++vi)
      {
        /* Konvektionsterm */
        eforce(vi*4)     += timefacfac*(-(velint_(0)*conv_r_(0, 0, vi)) - velint_(1)*conv_r_(0, 1, vi) - velint_(2)*conv_r_(0, 2, vi) + gridvelint_(0)*conv_r_(0, 0, vi) + gridvelint_(1)*conv_r_(0, 1, vi) + gridvelint_(2)*conv_r_(0, 2, vi)) ;
        eforce(vi*4 + 1) += timefacfac*(-(velint_(0)*conv_r_(1, 0, vi)) - velint_(1)*conv_r_(1, 1, vi) - velint_(2)*conv_r_(1, 2, vi) + gridvelint_(0)*conv_r_(1, 0, vi) + gridvelint_(1)*conv_r_(1, 1, vi) + gridvelint_(2)*conv_r_(1, 2, vi)) ;
        eforce(vi*4 + 2) += -(timefacfac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi) - gridvelint_(0)*conv_r_(2, 0, vi) - gridvelint_(1)*conv_r_(2, 1, vi) - gridvelint_(2)*conv_r_(2, 2, vi))) ;

        /* Stabilisierung der Konvektion ( L_conv_u) */
        eforce(vi*4)     += ttimetauM*(-(conv_c_(vi)*conv_old_(0)) + conv_g_(vi)*gridvelint_(0)*vderxy_(0, 0) - (gridvelint_(1)*gridvelint_(1))*derxy_(1, vi)*vderxy_(0, 1) - (gridvelint_(2)*gridvelint_(2))*derxy_(2, vi)*vderxy_(0, 2) + 2.0*velint_(0)*gridvelint_(0)*derxy_(0, vi)*vderxy_(0, 0) + velint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(0, 1) + velint_(0)*gridvelint_(1)*derxy_(1, vi)*vderxy_(0, 0) + velint_(1)*gridvelint_(0)*derxy_(0, vi)*vderxy_(0, 1) + velint_(1)*gridvelint_(0)*derxy_(1, vi)*vderxy_(0, 0) + velint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(0, 2) + velint_(0)*gridvelint_(2)*derxy_(2, vi)*vderxy_(0, 0) + velint_(2)*gridvelint_(0)*derxy_(0, vi)*vderxy_(0, 2) + velint_(2)*gridvelint_(0)*derxy_(2, vi)*vderxy_(0, 0) + 2.0*velint_(1)*gridvelint_(1)*derxy_(1, vi)*vderxy_(0, 1) + velint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(0, 2) + velint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(0, 1) + velint_(2)*gridvelint_(1)*derxy_(1, vi)*vderxy_(0, 2) + velint_(2)*gridvelint_(1)*derxy_(2, vi)*vderxy_(0, 1) + 2.0*velint_(2)*gridvelint_(2)*derxy_(2, vi)*vderxy_(0, 2) - gridvelint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(0, 1) - gridvelint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(0, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(0, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(0, 1)) ;
        eforce(vi*4 + 1) += ttimetauM*(-(conv_c_(vi)*conv_old_(1)) + conv_g_(vi)*gridvelint_(0)*vderxy_(1, 0) - (gridvelint_(1)*gridvelint_(1))*derxy_(1, vi)*vderxy_(1, 1) - (gridvelint_(2)*gridvelint_(2))*derxy_(2, vi)*vderxy_(1, 2) + 2.0*velint_(0)*gridvelint_(0)*derxy_(0, vi)*vderxy_(1, 0) + velint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(1, 1) + velint_(0)*gridvelint_(1)*derxy_(1, vi)*vderxy_(1, 0) + velint_(1)*gridvelint_(0)*derxy_(0, vi)*vderxy_(1, 1) + velint_(1)*gridvelint_(0)*derxy_(1, vi)*vderxy_(1, 0) + velint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(1, 2) + velint_(0)*gridvelint_(2)*derxy_(2, vi)*vderxy_(1, 0) + velint_(2)*gridvelint_(0)*derxy_(0, vi)*vderxy_(1, 2) + velint_(2)*gridvelint_(0)*derxy_(2, vi)*vderxy_(1, 0) + 2.0*velint_(1)*gridvelint_(1)*derxy_(1, vi)*vderxy_(1, 1) + velint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(1, 2) + velint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(1, 1) + velint_(2)*gridvelint_(1)*derxy_(1, vi)*vderxy_(1, 2) + velint_(2)*gridvelint_(1)*derxy_(2, vi)*vderxy_(1, 1) + 2.0*velint_(2)*gridvelint_(2)*derxy_(2, vi)*vderxy_(1, 2) - gridvelint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(1, 1) - gridvelint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(1, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(1, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(1, 1)) ;
        eforce(vi*4 + 2) += ttimetauM*(-(conv_c_(vi)*conv_old_(2)) + conv_g_(vi)*gridvelint_(0)*vderxy_(2, 0) - (gridvelint_(1)*gridvelint_(1))*derxy_(1, vi)*vderxy_(2, 1) - (gridvelint_(2)*gridvelint_(2))*derxy_(2, vi)*vderxy_(2, 2) + 2.0*velint_(0)*gridvelint_(0)*derxy_(0, vi)*vderxy_(2, 0) + velint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(2, 1) + velint_(0)*gridvelint_(1)*derxy_(1, vi)*vderxy_(2, 0) + velint_(1)*gridvelint_(0)*derxy_(0, vi)*vderxy_(2, 1) + velint_(1)*gridvelint_(0)*derxy_(1, vi)*vderxy_(2, 0) + velint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(2, 2) + velint_(0)*gridvelint_(2)*derxy_(2, vi)*vderxy_(2, 0) + velint_(2)*gridvelint_(0)*derxy_(0, vi)*vderxy_(2, 2) + velint_(2)*gridvelint_(0)*derxy_(2, vi)*vderxy_(2, 0) + 2.0*velint_(1)*gridvelint_(1)*derxy_(1, vi)*vderxy_(2, 1) + velint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(2, 2) + velint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(2, 1) + velint_(2)*gridvelint_(1)*derxy_(1, vi)*vderxy_(2, 2) + velint_(2)*gridvelint_(1)*derxy_(2, vi)*vderxy_(2, 1) + 2.0*velint_(2)*gridvelint_(2)*derxy_(2, vi)*vderxy_(2, 2) - gridvelint_(0)*gridvelint_(1)*derxy_(0, vi)*vderxy_(2, 1) - gridvelint_(0)*gridvelint_(2)*derxy_(0, vi)*vderxy_(2, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(1, vi)*vderxy_(2, 2) - gridvelint_(1)*gridvelint_(2)*derxy_(2, vi)*vderxy_(2, 1)) ;

        /* Stabilisierung der Konvektion (-L_visc_u) */
        eforce(vi*4)     += 2.0*visc*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(0) ;
        eforce(vi*4 + 1) += 2.0*visc*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(1) ;
        eforce(vi*4 + 2) += 2.0*visc*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(2) ;

        /* Stabilisierung der Konvektion ( L_pres_p) */
        eforce(vi*4)     += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(0)) ;
        eforce(vi*4 + 1) += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(1)) ;
        eforce(vi*4 + 2) += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(2)) ;

        /* Viskosit�tsterm */
        eforce(vi*4)     += -(visc*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0) + derxy_(2, vi)*vderxy_(0, 2) + derxy_(2, vi)*vderxy_(2, 0))) ;
        eforce(vi*4 + 1) += -(visc*timefacfac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1) + derxy_(2, vi)*vderxy_(1, 2) + derxy_(2, vi)*vderxy_(2, 1))) ;
        eforce(vi*4 + 2) += -(visc*timefacfac*(derxy_(0, vi)*vderxy_(0, 2) + derxy_(0, vi)*vderxy_(2, 0) + derxy_(1, vi)*vderxy_(1, 2) + derxy_(1, vi)*vderxy_(2, 1) + 2.0*derxy_(2, vi)*vderxy_(2, 2))) ;

        /* Stabilisierung der Viskosit�t ( L_conv_u) */
        eforce(vi*4)     += -2.0*visc*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi) + conv_old_(2)*viscs2_(0, 2, vi) - gridvelint_(0)*vderxy_(0, 0)*viscs2_(0, 0, vi) - gridvelint_(0)*vderxy_(1, 0)*viscs2_(0, 1, vi) - gridvelint_(1)*vderxy_(0, 1)*viscs2_(0, 0, vi) - gridvelint_(0)*vderxy_(2, 0)*viscs2_(0, 2, vi) - gridvelint_(2)*vderxy_(0, 2)*viscs2_(0, 0, vi) - gridvelint_(1)*vderxy_(1, 1)*viscs2_(0, 1, vi) - gridvelint_(1)*vderxy_(2, 1)*viscs2_(0, 2, vi) - gridvelint_(2)*vderxy_(1, 2)*viscs2_(0, 1, vi) - gridvelint_(2)*vderxy_(2, 2)*viscs2_(0, 2, vi)) ;
        eforce(vi*4 + 1) += 2.0*visc*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 1, vi)) - conv_old_(1)*viscs2_(1, 1, vi) - conv_old_(2)*viscs2_(1, 2, vi) + gridvelint_(0)*vderxy_(0, 0)*viscs2_(0, 1, vi) + gridvelint_(0)*vderxy_(1, 0)*viscs2_(1, 1, vi) + gridvelint_(1)*vderxy_(0, 1)*viscs2_(0, 1, vi) + gridvelint_(0)*vderxy_(2, 0)*viscs2_(1, 2, vi) + gridvelint_(2)*vderxy_(0, 2)*viscs2_(0, 1, vi) + gridvelint_(1)*vderxy_(1, 1)*viscs2_(1, 1, vi) + gridvelint_(1)*vderxy_(2, 1)*viscs2_(1, 2, vi) + gridvelint_(2)*vderxy_(1, 2)*viscs2_(1, 1, vi) + gridvelint_(2)*vderxy_(2, 2)*viscs2_(1, 2, vi)) ;
        eforce(vi*4 + 2) += 2.0*visc*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 2, vi)) - conv_old_(1)*viscs2_(1, 2, vi) - conv_old_(2)*viscs2_(2, 2, vi) + gridvelint_(0)*vderxy_(0, 0)*viscs2_(0, 2, vi) + gridvelint_(0)*vderxy_(1, 0)*viscs2_(1, 2, vi) + gridvelint_(1)*vderxy_(0, 1)*viscs2_(0, 2, vi) + gridvelint_(0)*vderxy_(2, 0)*viscs2_(2, 2, vi) + gridvelint_(2)*vderxy_(0, 2)*viscs2_(0, 2, vi) + gridvelint_(1)*vderxy_(1, 1)*viscs2_(1, 2, vi) + gridvelint_(1)*vderxy_(2, 1)*viscs2_(2, 2, vi) + gridvelint_(2)*vderxy_(1, 2)*viscs2_(1, 2, vi) + gridvelint_(2)*vderxy_(2, 2)*viscs2_(2, 2, vi)) ;

        /* Stabilisierung der Viskosit�t (-L_visc_u) */
        eforce(vi*4)     += 4.0*(visc*visc)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi) + visc_old_(2)*viscs2_(0, 2, vi)) ;
        eforce(vi*4 + 1) += 4.0*(visc*visc)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi) + visc_old_(2)*viscs2_(1, 2, vi)) ;
        eforce(vi*4 + 2) += 4.0*(visc*visc)*ttimetauMp*(visc_old_(0)*viscs2_(0, 2, vi) + visc_old_(1)*viscs2_(1, 2, vi) + visc_old_(2)*viscs2_(2, 2, vi)) ;

        /* Stabilisierung der Viskosit�t ( L_pres_p) */
        eforce(vi*4)     += -2.0*visc*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi) + gradp_(2)*viscs2_(0, 2, vi)) ;
        eforce(vi*4 + 1) += -2.0*visc*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi) + gradp_(2)*viscs2_(1, 2, vi)) ;
        eforce(vi*4 + 2) += -2.0*visc*ttimetauMp*(gradp_(0)*viscs2_(0, 2, vi) + gradp_(1)*viscs2_(1, 2, vi) + gradp_(2)*viscs2_(2, 2, vi)) ;

        /* Druckterm */
        eforce(vi*4)     += press*timefacfac*derxy_(0, vi) ;
        eforce(vi*4 + 1) += press*timefacfac*derxy_(1, vi) ;
        eforce(vi*4 + 2) += press*timefacfac*derxy_(2, vi) ;

        /* Stabilisierung des Drucks ( L_conv_u) */
        eforce(vi*4 + 3) += ttimetauMp*(-(conv_old_(0)*derxy_(0, vi)) -
                                        conv_old_(1)*derxy_(1, vi) -
                                        conv_old_(2)*derxy_(2, vi) +
                                        gridvelint_(0)*derxy_(0, vi)*vderxy_(0, 0) +
                                        gridvelint_(0)*derxy_(1, vi)*vderxy_(1, 0) +
                                        gridvelint_(1)*derxy_(0, vi)*vderxy_(0, 1) +
                                        gridvelint_(0)*derxy_(2, vi)*vderxy_(2, 0) +
                                        gridvelint_(2)*derxy_(0, vi)*vderxy_(0, 2) +
                                        gridvelint_(1)*derxy_(1, vi)*vderxy_(1, 1) +
                                        gridvelint_(1)*derxy_(2, vi)*vderxy_(2, 1) +
                                        gridvelint_(2)*derxy_(1, vi)*vderxy_(1, 2) +
                                        gridvelint_(2)*derxy_(2, vi)*vderxy_(2, 2)) ;

        /* Stabilisierung des Drucks (-L_visc_u) */
        eforce(vi*4 + 3) += 2.0*visc*ttimetauMp*(visc_old_(0)*derxy_(0, vi) + visc_old_(1)*derxy_(1, vi) + visc_old_(2)*derxy_(2, vi)) ;

        /* Stabilisierung des Drucks ( L_pres_p) */
        eforce(vi*4 + 3) += -(ttimetauMp*(gradp_(0)*derxy_(0, vi) + gradp_(1)*derxy_(1, vi) + gradp_(2)*derxy_(2, vi))) ;

        /* Divergenzfreiheit */
        eforce(vi*4 + 3) += -(timefacfac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi) + conv_r_(2, 2, vi))) ;

        /* Kontinuit�tsstabilisierung */
        eforce(vi*4)     += -((timefac*timefac)*tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2))) ;
        eforce(vi*4 + 1) += -((timefac*timefac)*tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2))) ;
        eforce(vi*4 + 2) += -((timefac*timefac)*tau_C*derxy_(2, vi)*(vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2))) ;

        /* Massenterm */
        eforce(vi*4)     += -(fac*funct_(vi)*velint_(0)) ;
        eforce(vi*4 + 1) += -(fac*funct_(vi)*velint_(1)) ;
        eforce(vi*4 + 2) += -(fac*funct_(vi)*velint_(2)) ;

        /* Konvektionsstabilisierung */
        eforce(vi*4)     += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(0)) ;
        eforce(vi*4 + 1) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(1)) ;
        eforce(vi*4 + 2) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(2)) ;

        /* Viskosit�tsstabilisierung */
        eforce(vi*4)     += -2.0*visc*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi) + velint_(2)*viscs2_(0, 2, vi)) ;
        eforce(vi*4 + 1) += -2.0*visc*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi) + velint_(2)*viscs2_(1, 2, vi)) ;
        eforce(vi*4 + 2) += -2.0*visc*timetauMp*(velint_(0)*viscs2_(0, 2, vi) + velint_(1)*viscs2_(1, 2, vi) + velint_(2)*viscs2_(2, 2, vi)) ;

        /* Stabilisierung der Druckgleichung */
        eforce(vi*4 + 3) += -(timetauMp*conv_c_(vi)) ;

        /* Quellterm der rechten Seite */
        eforce(vi*4)     += fac*funct_(vi)*rhsint_(0) ;
        eforce(vi*4 + 1) += fac*funct_(vi)*rhsint_(1) ;
        eforce(vi*4 + 2) += fac*funct_(vi)*rhsint_(2) ;

        /* Konvektionsstabilisierung */
        eforce(vi*4)     += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(0) ;
        eforce(vi*4 + 1) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(1) ;
        eforce(vi*4 + 2) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(2) ;

        /* Viskosit�tsstabilisierung */
        eforce(vi*4)     += 2.0*visc*timetauMp*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi) + rhsint_(2)*viscs2_(0, 2, vi)) ;
        eforce(vi*4 + 1) += 2.0*visc*timetauMp*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi) + rhsint_(2)*viscs2_(1, 2, vi)) ;
        eforce(vi*4 + 2) += 2.0*visc*timetauMp*(rhsint_(0)*viscs2_(0, 2, vi) + rhsint_(1)*viscs2_(1, 2, vi) + rhsint_(2)*viscs2_(2, 2, vi)) ;

        /* Stabilisierung der Druckgleichung */
        eforce(vi*4 + 3) += timetauMp*(rhsint_(0)*derxy_(0, vi) + rhsint_(1)*derxy_(1, vi) + rhsint_(2)*derxy_(2, vi)) ;

      }
    }
    else
#endif
    {
      // evaluate residual once for all stabilisation right hand sides
      res_old_ = velint_-rhsint_+timefac*(conv_old_+gradp_-2*visc*visc_old_);

      if (ele->is_ale_)
      {
        // correct convection with grid velocity
        res_old_ -= timefac*blitz::sum(vderxy_(i, j)*gridvelint_(j), j);
      }

      //----------------------------------------------------------------------
      //                            GALERKIN PART

      for (int ui=0; ui<iel_; ++ui)
      {
        for (int vi=0; vi<iel_; ++vi)
        {

          /* inertia (contribution to mass matrix) */
          /*

                               /        \
                              |          |
                              |  Du , v  |
                              |          |
                               \        /
          */
          estif(vi*4, ui*4)         += fac*funct_(ui)*funct_(vi) ;
          estif(vi*4 + 1, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
          estif(vi*4 + 2, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;

          /* convection, convective part */
          /*

                       /                       \
                      |  / n+1       \          |
                      | | u   o nabla | Du , v  |
                      |  \ (i)       /          |
                       \                       /

          */
          estif(vi*4, ui*4)         += timefacfac*funct_(vi)*conv_c_(ui) ;
          estif(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*conv_c_(ui) ;
          estif(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*conv_c_(ui) ;

          /* Viskosit�tsterm */
          /*
                        /                        \
                       |       /  \         / \   |
                       |  eps | Du | , eps | v |  |
                       |       \  /         \ /   |
                        \                        /
          */
          estif(vi*4, ui*4)         += visc*timefacfac*(2.0*derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4, ui*4 + 1)     += visc*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*4, ui*4 + 2)     += visc*timefacfac*derxy_(0, ui)*derxy_(2, vi) ;
          estif(vi*4 + 1, ui*4)     += visc*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*4 + 1, ui*4 + 1) += visc*timefacfac*(derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        2.0*derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += visc*timefacfac*derxy_(1, ui)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4)     += visc*timefacfac*derxy_(0, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 1) += visc*timefacfac*derxy_(1, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 2) += visc*timefacfac*(derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        2.0*derxy_(2, ui)*derxy_(2, vi)) ;

          /* Druckterm */
          /*

                          /                \
                         |                  |
                         |  Dp , nabla o v  |
                         |                  |
                          \                /
          */

          estif(vi*4, ui*4 + 3)     += -(timefacfac*funct_(ui)*derxy_(0, vi)) ;
          estif(vi*4 + 1, ui*4 + 3) += -(timefacfac*funct_(ui)*derxy_(1, vi)) ;
          estif(vi*4 + 2, ui*4 + 3) += -(timefacfac*funct_(ui)*derxy_(2, vi)) ;


          /* Divergenzfreiheit */
          /*
                         /                \
                        |                  |
                        | nabla o Du  , q  |
                        |                  |
                         \                /
          */
          estif(vi*4 + 3, ui*4)     += timefacfac*funct_(vi)*derxy_(0, ui) ;
          estif(vi*4 + 3, ui*4 + 1) += timefacfac*funct_(vi)*derxy_(1, ui) ;
          estif(vi*4 + 3, ui*4 + 2) += timefacfac*funct_(vi)*derxy_(2, ui) ;

        }
      }

      if (ele->is_ale_)
      {
        for (int vi=0; vi<iel_; ++vi)
        {
          for (int ui=0; ui<iel_; ++ui)
          {

            /*  reduced convection through grid motion


                       /                      \
                      |  /          \          |
                    - | | u  o nabla | Du , v  |
                      |  \ G        /          |
                       \                      /

            */
            estif(vi*4, ui*4)         += timefacfac*funct_(vi)*conv_g_(ui) ;
            estif(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*conv_g_(ui) ;
            estif(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*conv_g_(ui)  ;
          }

        }
      } // end if is_ale

      if (newton)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {

            /*  convection, reactive part

                   /                         \
                  |  /          \   n+1       |
                  | | Du o nabla | u     , v  |
                  |  \          /   (i)       |
                   \                         /
            */
            estif(vi*4, ui*4)         += timefacfac*funct_(vi)*conv_r_(0, 0, ui) ;
            estif(vi*4, ui*4 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
            estif(vi*4, ui*4 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
            estif(vi*4 + 1, ui*4)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
            estif(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(1, 1, ui) ;
            estif(vi*4 + 1, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
            estif(vi*4 + 2, ui*4)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
            estif(vi*4 + 2, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
            estif(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(2, 2, ui) ;
          }
        }
      }

      for (int vi=0; vi<iel_; ++vi)
      {
        /* inertia */
        eforce(vi*4)     += -(fac*funct_(vi)*velint_(0)) ;
        eforce(vi*4 + 1) += -(fac*funct_(vi)*velint_(1)) ;
        eforce(vi*4 + 2) += -(fac*funct_(vi)*velint_(2)) ;

        /* convection */
        eforce(vi*4)     += -(timefacfac*(velint_(0)*conv_r_(0, 0, vi)
                                          +
                                          velint_(1)*conv_r_(0, 1, vi)
                                          +
                                          velint_(2)*conv_r_(0, 2, vi))) ;
        eforce(vi*4 + 1) += -(timefacfac*(velint_(0)*conv_r_(1, 0, vi)
                                          +
                                          velint_(1)*conv_r_(1, 1, vi)
                                          +
                                          velint_(2)*conv_r_(1, 2, vi))) ;
        eforce(vi*4 + 2) += -(timefacfac*(velint_(0)*conv_r_(2, 0, vi)
                                          +
                                          velint_(1)*conv_r_(2, 1, vi)
                                          +
                                          velint_(2)*conv_r_(2, 2, vi))) ;

        /* pressure */
        eforce(vi*4)     += press*timefacfac*derxy_(0, vi) ;
        eforce(vi*4 + 1) += press*timefacfac*derxy_(1, vi) ;
        eforce(vi*4 + 2) += press*timefacfac*derxy_(2, vi) ;

        /* viscosity */
        eforce(vi*4)     += -(visc*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                                               +
                                               derxy_(1, vi)*vderxy_(0, 1)
                                               +
                                               derxy_(1, vi)*vderxy_(1, 0)
                                               +
                                               derxy_(2, vi)*vderxy_(0, 2)
                                               +
                                               derxy_(2, vi)*vderxy_(2, 0))) ;
        eforce(vi*4 + 1) += -(visc*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
                                               +
                                               derxy_(0, vi)*vderxy_(1, 0)
                                               +
                                               2.0*derxy_(1, vi)*vderxy_(1, 1)
                                               +
                                               derxy_(2, vi)*vderxy_(1, 2)
                                               +
                                               derxy_(2, vi)*vderxy_(2, 1))) ;
        eforce(vi*4 + 2) += -(visc*timefacfac*(derxy_(0, vi)*vderxy_(0, 2)
                                               +
                                               derxy_(0, vi)*vderxy_(2, 0)
                                               +
                                               derxy_(1, vi)*vderxy_(1, 2)
                                               +
                                               derxy_(1, vi)*vderxy_(2, 1)
                                               +
                                               2.0*derxy_(2, vi)*vderxy_(2, 2))) ;

        // source term of the right hand side
        eforce(vi*4)     += fac*funct_(vi)*rhsint_(0) ;
        eforce(vi*4 + 1) += fac*funct_(vi)*rhsint_(1) ;
        eforce(vi*4 + 2) += fac*funct_(vi)*rhsint_(2) ;

        // continuity equation
        eforce(vi*4 + 3) += -(timefacfac*(conv_r_(0, 0, vi)
                                          +
                                          conv_r_(1, 1, vi)
                                          +
                                          conv_r_(2, 2, vi))) ;
      } // vi

      if (ele->is_ale_)
      {
        for (int vi=0; vi<iel_; ++vi)
        {
          eforce(vi*4)     += timefacfac*(gridvelint_(0)*conv_r_(0, 0, vi)
                                          +
                                          gridvelint_(1)*conv_r_(0, 1, vi)
                                          +
                                          gridvelint_(2)*conv_r_(0, 2, vi)) ;
          eforce(vi*4 + 1) += timefacfac*(gridvelint_(0)*conv_r_(1, 0, vi)
                                          +
                                          gridvelint_(1)*conv_r_(1, 1, vi)
                                          +
                                          gridvelint_(2)*conv_r_(1, 2, vi)) ;
          eforce(vi*4 + 2) += timefacfac*(gridvelint_(0)*conv_r_(2, 0, vi)
                                          +
                                          gridvelint_(1)*conv_r_(2, 1, vi)
                                          +
                                          gridvelint_(2)*conv_r_(2, 2, vi)) ;
        } // vi
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if(pstab)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {

            /* pressure stabilisation: inertia */
            /*
                        /              \
                       |                |
                       |  Du , nabla q  |
                       |                |
                        \              /
            */

            estif(vi*4 + 3, ui*4)     += timetauMp*funct_(ui)*derxy_(0, vi) ;
            estif(vi*4 + 3, ui*4 + 1) += timetauMp*funct_(ui)*derxy_(1, vi) ;
            estif(vi*4 + 3, ui*4 + 2) += timetauMp*funct_(ui)*derxy_(2, vi) ;

            /* pressure stabilisation: convection, convective part */
            /*

                      /                            \
                     |  / n+1       \               |
                     | | u   o nabla | Du , nabla q |
                     |  \ (i)       /               |
                      \                            /

            */
            estif(vi*4 + 3, ui*4)     += ttimetauMp*conv_c_(ui)*derxy_(0, vi) ;
            estif(vi*4 + 3, ui*4 + 1) += ttimetauMp*conv_c_(ui)*derxy_(1, vi) ;
            estif(vi*4 + 3, ui*4 + 2) += ttimetauMp*conv_c_(ui)*derxy_(2, vi) ;


            /* pressure stabilisation: viscosity (-L_visc_u) */
            /*
                     /                              \
                    |               /  \             |
                    |  nabla o eps | Du | , nabla q  |
                    |               \  /             |
                     \                              /
            */
            estif(vi*4 + 3, ui*4)     += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 0, ui)
                                                              +
                                                              derxy_(1, vi)*viscs2_(0, 1, ui)
                                                              +
                                                              derxy_(2, vi)*viscs2_(0, 2, ui)) ;
            estif(vi*4 + 3, ui*4 + 1) += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 1, ui)
                                                              +
                                                              derxy_(1, vi)*viscs2_(1, 1, ui)
                                                              +
                                                              derxy_(2, vi)*viscs2_(1, 2, ui)) ;
            estif(vi*4 + 3, ui*4 + 2) += 2.0*visc*ttimetauMp*(derxy_(0, vi)*viscs2_(0, 2, ui)
                                                              +
                                                              derxy_(1, vi)*viscs2_(1, 2, ui)
                                                              +
                                                              derxy_(2, vi)*viscs2_(2, 2, ui)) ;

            /* pressure stabilisation: pressure( L_pres_p) */
            /*
                        /                    \
                       |                      |
                       |  nabla Dp , nabla q  |
                       |                      |
                        \                    /
            */
            estif(vi*4 + 3, ui*4 + 3) += ttimetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                     +
                                                     derxy_(1, ui)*derxy_(1, vi)
                                                     +
                                                     derxy_(2, ui)*derxy_(2, vi)) ;


          } // vi
        } // ui

        if (ele->is_ale_)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            for (int ui=0; ui<iel_; ++ui)
            {
              /*  reduced convection through grid motion


                       /                          \
                      |  /          \              |
                    - | | u  o nabla | Du , grad q |
                      |  \ G        /              |
                       \                          /

              */
              estif(vi*4 + 3, ui*4)     += ttimetauMp*conv_g_(ui)*derxy_(0, vi) ;
              estif(vi*4 + 3, ui*4 + 1) += ttimetauMp*conv_g_(ui)*derxy_(1, vi) ;
              estif(vi*4 + 3, ui*4 + 2) += ttimetauMp*conv_g_(ui)*derxy_(2, vi) ;
            }// ui
          }//vi
        } // end if is_ale

        if (newton)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            for (int vi=0; vi<iel_; ++vi)
            {
              /*  pressure stabilisation: convection, reactive part

                  /                             \
                 |  /          \   n+1           |
                 | | Du o nabla | u     , grad q |
                 |  \          /   (i)           |
                  \                             /

              */
              estif(vi*4 + 3, ui*4)     += ttimetauMp*(derxy_(0, vi)*conv_r_(0, 0, ui)
                                                       +
                                                       derxy_(1, vi)*conv_r_(1, 0, ui)
                                                       +
                                                       derxy_(2, vi)*conv_r_(2, 0, ui)) ;
              estif(vi*4 + 3, ui*4 + 1) += ttimetauMp*(derxy_(0, vi)*conv_r_(0, 1, ui)
                                                       +
                                                       derxy_(1, vi)*conv_r_(1, 1, ui)
                                                       +
                                                       derxy_(2, vi)*conv_r_(2, 1, ui)) ;
              estif(vi*4 + 3, ui*4 + 2) += ttimetauMp*(derxy_(0, vi)*conv_r_(0, 2, ui)
                                                       +
                                                       derxy_(1, vi)*conv_r_(1, 2, ui)
                                                       +
                                                       derxy_(2, vi)*conv_r_(2, 2, ui)) ;

            } // vi
          } // ui
        } // if newton


        for (int vi=0; vi<iel_; ++vi)
        {
          // pressure stabilisation
          eforce(vi*4 + 3) -= timetauMp*(res_old_(0)*derxy_(0, vi)
                                         +
                                         res_old_(1)*derxy_(1, vi)
                                         +
                                         res_old_(2)*derxy_(2, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            /* supg stabilisation: inertia  */
            /*
                      /                        \
                     |        / n+1       \     |
                     |  Du , | u   o nabla | v  |
                     |        \ (i)       /     |
                      \                        /
            */
            estif(vi*4, ui*4)         += timetauM*funct_(ui)*conv_c_(vi);
            estif(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*conv_c_(vi);
            estif(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*conv_c_(vi);

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

                 /                                           \
                |    / n+1        \        / n+1        \     |
                |   | u    o nabla | Du , | u    o nabla | v  |
                |    \ (i)        /        \ (i)        /     |
                 \                                           /

            */

            estif(vi*4, ui*4)         += ttimetauM*conv_c_(ui)*conv_c_(vi) ;
            estif(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_c_(ui)*conv_c_(vi) ;
            estif(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_c_(ui)*conv_c_(vi) ;

            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
                      /                              \
                     |              / n+1       \     |
                     |  nabla Dp , | u   o nabla | v  |
                     |              \ (i)       /     |
                      \                              /
            */
            estif(vi*4, ui*4 + 3)     += ttimetauM*conv_c_(vi)*derxy_(0, ui) ;
            estif(vi*4 + 1, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxy_(1, ui) ;
            estif(vi*4 + 2, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxy_(2, ui) ;

            /* supg stabilisation: viscous part  (-L_visc_u) */
            /*
                  /                                        \
                 |               /  \    / n+1        \     |
                 |  nabla o eps | Du |, | u    o nabla | v  |
                 |               \  /    \ (i)        /     |
                  \                                        /
            */
            estif(vi*4, ui*4)         += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(0, 0, ui) ;
            estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(0, 1, ui) ;
            estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(0, 2, ui) ;
            estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(0, 1, ui) ;
            estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(1, 1, ui) ;
            estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(1, 2, ui) ;
            estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(0, 2, ui) ;
            estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(1, 2, ui) ;
            estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauM*conv_c_(vi)*viscs2_(2, 2, ui) ;


          } // vi
        } // ui

        if (ele->is_ale_)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            for (int ui=0; ui<iel_; ++ui)
            {


              estif(vi*4, ui*4)         += ttimetauM*conv_g_(ui)*conv_c_(vi) ;
              estif(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_g_(ui)*conv_c_(vi) ;
              estif(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_g_(ui)*conv_c_(vi) ;


              /*
                      /                       \
                     |        /          \     |
                     |  Du , | u  o nabla | v  |
                     |        \ G        /     |
                      \                       /
              */

              estif(vi*4, ui*4)         += timetauM*funct_(ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*conv_g_(vi) ;

              /*

                 /                                         \
                |    / n+1        \        /          \     |
                |   | u    o nabla | Du , | u  o nabla | v  |
                |    \ (i)        /        \ G        /     |
                 \                                         /

              */
              estif(vi*4, ui*4)         += ttimetauM*conv_c_(ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_c_(ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_c_(ui)*conv_g_(vi) ;

              /*
                      /                             \
                     |              /          \     |
                     |  nabla Dp , | u  o nabla | v  |
                     |              \ G        /     |
                      \                             /
              */
              estif(vi*4, ui*4 + 3)     += ttimetauM*derxy_(0, ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 3) += ttimetauM*derxy_(1, ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 3) += ttimetauM*derxy_(2, ui)*conv_g_(vi) ;

              /*
                  /                                      \
                 |               /  \    /          \     |
                 |  nabla o eps | Du |, | u  o nabla | v  |
                 |               \  /    \ G        /     |
                  \                                      /
              */
              estif(vi*4, ui*4)         += 2.0*visc*ttimetauM*viscs2_(0, 0, ui)*conv_g_(vi) ;
              estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauM*viscs2_(0, 1, ui)*conv_g_(vi) ;
              estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauM*viscs2_(0, 2, ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauM*viscs2_(0, 1, ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauM*viscs2_(1, 1, ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauM*viscs2_(1, 2, ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauM*viscs2_(0, 2, ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauM*viscs2_(1, 2, ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauM*viscs2_(2, 2, ui)*conv_g_(vi) ;

              /*

                 /                                       \
                |    /          \        /          \     |
                |   | u  o nabla | Du , | u  o nabla | v  |
                |    \ G        /        \ G        /     |
                 \                                       /

              */
              estif(vi*4, ui*4)         += ttimetauM*conv_g_(ui)*conv_g_(vi) ;
              estif(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_g_(ui)*conv_g_(vi) ;
              estif(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_g_(ui)*conv_g_(vi) ;
            }// ui
          }//vi


          if (newton)
          {
            for (int vi=0; vi<iel_; ++vi)
            {
              for (int ui=0; ui<iel_; ++ui)
              {
                /*
                       /                                         \
                      |    /          \   n+1    /          \     |
                      |   | u  o nabla | u    , | Du o nabla | v  |
                      |    \ G        /   (i)    \          /     |
                       \                                         /

                */
                estif(vi*4, ui*4)         += ttimetauM*( - gridvelint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui)
                                                         - gridvelint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui)
                                                         - gridvelint_(2)*derxy_(0, vi)*conv_r_(0, 2, ui)) ;
                estif(vi*4, ui*4 + 1)     += ttimetauM*( - gridvelint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui)
                                                         - gridvelint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui)
                                                         - gridvelint_(2)*derxy_(1, vi)*conv_r_(0, 2, ui)) ;
                estif(vi*4, ui*4 + 2)     += ttimetauM*( - gridvelint_(0)*derxy_(2, vi)*conv_r_(0, 0, ui)
                                                         - gridvelint_(1)*derxy_(2, vi)*conv_r_(0, 1, ui)
                                                         - gridvelint_(2)*derxy_(2, vi)*conv_r_(0, 2, ui)) ;
                estif(vi*4 + 1, ui*4)     += ttimetauM*( - gridvelint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui)
                                                         - gridvelint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui)
                                                         - gridvelint_(2)*derxy_(0, vi)*conv_r_(1, 2, ui)) ;
                estif(vi*4 + 1, ui*4 + 1) += ttimetauM*( - gridvelint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui)
                                                         - gridvelint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui)
                                                         - gridvelint_(2)*derxy_(1, vi)*conv_r_(1, 2, ui)) ;
                estif(vi*4 + 1, ui*4 + 2) += ttimetauM*( - gridvelint_(0)*derxy_(2, vi)*conv_r_(1, 0, ui)
                                                         - gridvelint_(1)*derxy_(2, vi)*conv_r_(1, 1, ui)
                                                         - gridvelint_(2)*derxy_(2, vi)*conv_r_(1, 2, ui)) ;
                estif(vi*4 + 2, ui*4)     += ttimetauM*( - gridvelint_(0)*derxy_(0, vi)*conv_r_(2, 0, ui)
                                                         - gridvelint_(1)*derxy_(0, vi)*conv_r_(2, 1, ui)
                                                         - gridvelint_(2)*derxy_(0, vi)*conv_r_(2, 2, ui)) ;
                estif(vi*4 + 2, ui*4 + 1) += ttimetauM*( - gridvelint_(0)*derxy_(1, vi)*conv_r_(2, 0, ui)
                                                         - gridvelint_(1)*derxy_(1, vi)*conv_r_(2, 1, ui)
                                                         - gridvelint_(2)*derxy_(1, vi)*conv_r_(2, 2, ui)) ;
                estif(vi*4 + 2, ui*4 + 2) += ttimetauM*( - gridvelint_(0)*derxy_(2, vi)*conv_r_(2, 0, ui)
                                                         - gridvelint_(1)*derxy_(2, vi)*conv_r_(2, 1, ui)
                                                         - gridvelint_(2)*derxy_(2, vi)*conv_r_(2, 2, ui)) ;

                /*
                       /                                         \
                      |    /          \   n+1    /          \     |
                      |   | Du o nabla | u    , | u  o nabla | v  |
                      |    \          /   (i)    \ G        /     |
                       \                                         /
                */
                estif(vi*4, ui*4)         += ttimetauM*conv_r_(0, 0, ui)*conv_g_(vi) ;
                estif(vi*4, ui*4 + 1)     += ttimetauM*conv_r_(0, 1, ui)*conv_g_(vi) ;
                estif(vi*4, ui*4 + 2)     += ttimetauM*conv_r_(0, 2, ui)*conv_g_(vi) ;
                estif(vi*4 + 1, ui*4)     += ttimetauM*conv_r_(1, 0, ui)*conv_g_(vi) ;
                estif(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_r_(1, 1, ui)*conv_g_(vi) ;
                estif(vi*4 + 1, ui*4 + 2) += ttimetauM*conv_r_(1, 2, ui)*conv_g_(vi) ;
                estif(vi*4 + 2, ui*4)     += ttimetauM*conv_r_(2, 0, ui)*conv_g_(vi) ;
                estif(vi*4 + 2, ui*4 + 1) += ttimetauM*conv_r_(2, 1, ui)*conv_g_(vi) ;
                estif(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_r_(2, 2, ui)*conv_g_(vi) ;
              }// ui
            }//vi
          } // end if newton
        } // end if is_ale




        if (newton)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            for (int vi=0; vi<iel_; ++vi)
            {
              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                         /                           \
                        |   n+1      /          \     |
                        |  u      , | Du o nabla | v  |
                        |   (i)      \          /     |
                         \                           /

              */
              estif(vi*4, ui*4)         += timetauM*funct_(ui)*velint_(0)*derxy_(0, vi) ;
              estif(vi*4, ui*4 + 1)     += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
              estif(vi*4, ui*4 + 2)     += timetauM*funct_(ui)*velint_(0)*derxy_(2, vi) ;
              estif(vi*4 + 1, ui*4)     += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
              estif(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*velint_(1)*derxy_(1, vi) ;
              estif(vi*4 + 1, ui*4 + 2) += timetauM*funct_(ui)*velint_(1)*derxy_(2, vi) ;
              estif(vi*4 + 2, ui*4)     += timetauM*funct_(ui)*velint_(2)*derxy_(0, vi) ;
              estif(vi*4 + 2, ui*4 + 1) += timetauM*funct_(ui)*velint_(2)*derxy_(1, vi) ;
              estif(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*velint_(2)*derxy_(2, vi) ;


              /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
              /*
                       /                                           \
                      |    / n+1        \   n+1    /          \     |
                      |   | u    o nabla | u    , | Du o nabla | v  |
                      |    \ (i)        /   (i)    \          /     |
                       \                                           /

                       /                                           \
                      |    /          \   n+1    / n+1        \     |
                      |   | Du o nabla | u    , | u    o nabla | v  |
                      |    \          /   (i)    \ (i)        /     |
                       \                                           /
              */
              estif(vi*4, ui*4)         += ttimetauM*(conv_c_(vi)*conv_r_(0, 0, ui)
                                                      +
                                                      velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(0, vi)*conv_r_(0, 2, ui)) ;
              estif(vi*4, ui*4 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui)
                                                      +
                                                      velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(1, vi)*conv_r_(0, 2, ui)) ;
              estif(vi*4, ui*4 + 2)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui)
                                                      +
                                                      velint_(0)*derxy_(2, vi)*conv_r_(0, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(2, vi)*conv_r_(0, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(2, vi)*conv_r_(0, 2, ui)) ;
              estif(vi*4 + 1, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui)
                                                      +
                                                      velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(0, vi)*conv_r_(1, 2, ui)) ;
              estif(vi*4 + 1, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(1, 1, ui)
                                                      +
                                                      velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(1, vi)*conv_r_(1, 2, ui)) ;
              estif(vi*4 + 1, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui)
                                                      +
                                                      velint_(0)*derxy_(2, vi)*conv_r_(1, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(2, vi)*conv_r_(1, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(2, vi)*conv_r_(1, 2, ui)) ;
              estif(vi*4 + 2, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui)
                                                      +
                                                      velint_(0)*derxy_(0, vi)*conv_r_(2, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(0, vi)*conv_r_(2, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(0, vi)*conv_r_(2, 2, ui)) ;
              estif(vi*4 + 2, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui)
                                                      +
                                                      velint_(0)*derxy_(1, vi)*conv_r_(2, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(1, vi)*conv_r_(2, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(1, vi)*conv_r_(2, 2, ui)) ;
              estif(vi*4 + 2, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(2, 2, ui)
                                                      +
                                                      velint_(0)*derxy_(2, vi)*conv_r_(2, 0, ui)
                                                      +
                                                      velint_(1)*derxy_(2, vi)*conv_r_(2, 1, ui)
                                                      +
                                                      velint_(2)*derxy_(2, vi)*conv_r_(2, 2, ui)) ;


              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                            /                               \
                           |         n+1    /          \     |
                           |  nabla p    , | Du o nabla | v  |
                           |         (i)    \          /     |
                            \                               /
              */
              estif(vi*4, ui*4)         += ttimetauM*funct_(ui)*gradp_(0)*derxy_(0, vi) ;
              estif(vi*4, ui*4 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxy_(1, vi) ;
              estif(vi*4, ui*4 + 2)     += ttimetauM*funct_(ui)*gradp_(0)*derxy_(2, vi) ;
              estif(vi*4 + 1, ui*4)     += ttimetauM*funct_(ui)*gradp_(1)*derxy_(0, vi) ;
              estif(vi*4 + 1, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxy_(1, vi) ;
              estif(vi*4 + 1, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(1)*derxy_(2, vi) ;
              estif(vi*4 + 2, ui*4)     += ttimetauM*funct_(ui)*gradp_(2)*derxy_(0, vi) ;
              estif(vi*4 + 2, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(2)*derxy_(1, vi) ;
              estif(vi*4 + 2, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(2)*derxy_(2, vi) ;


              /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
              /*
                      /                                         \
                     |               / n+1 \    /          \     |
                     |  nabla o eps | u     |, | Du o nabla | v  |
                     |               \ (i) /    \          /     |
                      \                                         /
              */
              estif(vi*4, ui*4)         += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(0)*derxy_(0, vi) ;
              estif(vi*4, ui*4 + 1)     += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(0)*derxy_(1, vi) ;
              estif(vi*4, ui*4 + 2)     += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(0)*derxy_(2, vi) ;
              estif(vi*4 + 1, ui*4)     += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(1)*derxy_(0, vi) ;
              estif(vi*4 + 1, ui*4 + 1) += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(1)*derxy_(1, vi) ;
              estif(vi*4 + 1, ui*4 + 2) += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(1)*derxy_(2, vi) ;
              estif(vi*4 + 2, ui*4)     += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(2)*derxy_(0, vi) ;
              estif(vi*4 + 2, ui*4 + 1) += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(2)*derxy_(1, vi) ;
              estif(vi*4 + 2, ui*4 + 2) += -2.0*visc*ttimetauM*funct_(ui)*visc_old_(2)*derxy_(2, vi) ;


              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                          /                             \
                         |              /          \     |
                         |  rhsint   , | Du o nabla | v  |
                         |              \          /     |
                          \                             /

              */
              estif(vi*4    , ui*4)     += -(timetauM*funct_(ui)*derxy_(0, vi)*rhsint_(0)) ;
              estif(vi*4    , ui*4 + 1) += -(timetauM*funct_(ui)*derxy_(1, vi)*rhsint_(0)) ;
              estif(vi*4    , ui*4 + 2) += -(timetauM*funct_(ui)*derxy_(2, vi)*rhsint_(0)) ;
              estif(vi*4 + 1, ui*4)     += -(timetauM*funct_(ui)*derxy_(0, vi)*rhsint_(1)) ;
              estif(vi*4 + 1, ui*4 + 1) += -(timetauM*funct_(ui)*derxy_(1, vi)*rhsint_(1)) ;
              estif(vi*4 + 1, ui*4 + 2) += -(timetauM*funct_(ui)*derxy_(2, vi)*rhsint_(1)) ;
              estif(vi*4 + 2, ui*4)     += -(timetauM*funct_(ui)*derxy_(0, vi)*rhsint_(2)) ;
              estif(vi*4 + 2, ui*4 + 1) += -(timetauM*funct_(ui)*derxy_(1, vi)*rhsint_(2)) ;
              estif(vi*4 + 2, ui*4 + 2) += -(timetauM*funct_(ui)*derxy_(2, vi)*rhsint_(2)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<iel_; ++vi)
        {
          // supg stabilisation
          eforce(vi*4)     += -(timetauM*conv_c_(vi)*res_old_(0)) ;
          eforce(vi*4 + 1) += -(timetauM*conv_c_(vi)*res_old_(1)) ;
          eforce(vi*4 + 2) += -(timetauM*conv_c_(vi)*res_old_(2)) ;
        }
      }


      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART
      if(vstab)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            /* viscous stabilisation, inertia part */
            /*
                        /                  \
                       |                    |
                       |  Du , div eps (v)  |
                       |                    |
                        \                  /
            */
            estif(vi*4, ui*4)         += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
            estif(vi*4, ui*4 + 1)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
            estif(vi*4, ui*4 + 2)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
            estif(vi*4 + 1, ui*4)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
            estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
            estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
            estif(vi*4 + 2, ui*4)     += 2.0*visc*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
            estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
            estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;

            /* viscous stabilisation, convective part */
            /*
                 /                                \
                |  / n+1       \                   |
                | | u   o nabla | Du , div eps (v) |
                |  \ (i)       /                   |
                 \                                /
            */
            estif(vi*4, ui*4)         += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(0, 0, vi) ;
            estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(0, 1, vi) ;
            estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(0, 2, vi) ;
            estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(0, 1, vi) ;
            estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(1, 1, vi) ;
            estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(1, 2, vi) ;
            estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(0, 2, vi) ;
            estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(1, 2, vi) ;
            estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauMp*conv_c_(ui)*viscs2_(2, 2, vi) ;


            /* viscous stabilisation, pressure part ( L_pres_p) */
            /*
                     /                        \
                    |                          |
                    |  nabla Dp , div eps (v)  |
                    |                          |
                     \                        /
            */
            estif(vi*4, ui*4 + 3)     += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 0, vi)
                                                              +
                                                              derxy_(1, ui)*viscs2_(0, 1, vi)
                                                              +
                                                              derxy_(2, ui)*viscs2_(0, 2, vi)) ;
            estif(vi*4 + 1, ui*4 + 3) += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 1, vi)
                                                              +
                                                              derxy_(1, ui)*viscs2_(1, 1, vi)
                                                              +
                                                              derxy_(2, ui)*viscs2_(1, 2, vi)) ;
            estif(vi*4 + 2, ui*4 + 3) += 2.0*visc*ttimetauMp*(derxy_(0, ui)*viscs2_(0, 2, vi)
                                                              +
                                                              derxy_(1, ui)*viscs2_(1, 2, vi)
                                                              +
                                                              derxy_(2, ui)*viscs2_(2, 2, vi)) ;

            /* viscous stabilisation, viscous part (-L_visc_u) */
            /*
               /                                 \
              |               /  \                |
              |  nabla o eps | Du | , div eps (v) |
              |               \  /                |
               \                                 /
            */
            estif(vi*4, ui*4)         += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi)
                                                                     +
                                                                     viscs2_(0, 1, ui)*viscs2_(0, 1, vi)
                                                                     +
                                                                     viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
            estif(vi*4, ui*4 + 1)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui)
                                                                     +
                                                                     viscs2_(0, 1, vi)*viscs2_(1, 1, ui)
                                                                     +
                                                                     viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
            estif(vi*4, ui*4 + 2)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui)
                                                                     +
                                                                     viscs2_(0, 1, vi)*viscs2_(1, 2, ui)
                                                                     +
                                                                     viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
            estif(vi*4 + 1, ui*4)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi)
                                                                     +
                                                                     viscs2_(0, 1, ui)*viscs2_(1, 1, vi)
                                                                     +
                                                                     viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
            estif(vi*4 + 1, ui*4 + 1) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi)
                                                                     +
                                                                     viscs2_(1, 1, ui)*viscs2_(1, 1, vi)
                                                                     +
                                                                     viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
            estif(vi*4 + 1, ui*4 + 2) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui)
                                                                     +
                                                                     viscs2_(1, 1, vi)*viscs2_(1, 2, ui)
                                                                     +
                                                                     viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
            estif(vi*4 + 2, ui*4)     += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi)
                                                                     +
                                                                     viscs2_(0, 1, ui)*viscs2_(1, 2, vi)
                                                                     +
                                                                     viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
            estif(vi*4 + 2, ui*4 + 1) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi)
                                                                     +
                                                                     viscs2_(1, 1, ui)*viscs2_(1, 2, vi)
                                                                     +
                                                                     viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
            estif(vi*4 + 2, ui*4 + 2) += 4.0*(visc*visc)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi)
                                                                     +
                                                                     viscs2_(1, 2, ui)*viscs2_(1, 2, vi)
                                                                     +
                                                                     viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;
          } // vi
        } // ui

        if (ele->is_ale_)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            for (int ui=0; ui<iel_; ++ui)
            {

              /*  reduced convection through grid motion


                       /                                   \
                      |  /          \                       |
                      | | u  o nabla | Du ,  nabla o eps (v)|
                      |  \ G        /                       |
                       \                                   /

              */

              estif(vi*4, ui*4)         += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(0, 0, vi) ;
              estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(0, 1, vi) ;
              estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(0, 2, vi) ;
              estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(0, 1, vi) ;
              estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(1, 1, vi) ;
              estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(1, 2, vi) ;
              estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(0, 2, vi) ;
              estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(1, 2, vi) ;
              estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauMp*conv_g_(ui)*viscs2_(2, 2, vi) ;
            }// ui
          }//vi
        } // end if is_ale


        if (newton)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            for (int vi=0; vi<iel_; ++vi)
            {
              /* viscous stabilisation, reactive part of convection */
              /*
                   /                                 \
                  |  /          \   n+1               |
                  | | Du o nabla | u    , div eps (v) |
                  |  \          /   (i)               |
                   \                                 /
              */
              estif(vi*4, ui*4)         += 2.0*visc*ttimetauMp*(viscs2_(0, 0, vi)*conv_r_(0, 0, ui)
                                                                +
                                                                viscs2_(0, 1, vi)*conv_r_(1, 0, ui)
                                                                +
                                                                viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
              estif(vi*4, ui*4 + 1)     += 2.0*visc*ttimetauMp*(viscs2_(0, 0, vi)*conv_r_(0, 1, ui)
                                                                +
                                                                viscs2_(0, 1, vi)*conv_r_(1, 1, ui)
                                                                +
                                                                viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
              estif(vi*4, ui*4 + 2)     += 2.0*visc*ttimetauMp*(viscs2_(0, 0, vi)*conv_r_(0, 2, ui)
                                                                +
                                                                viscs2_(0, 1, vi)*conv_r_(1, 2, ui)
                                                                +
                                                                viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
              estif(vi*4 + 1, ui*4)     += 2.0*visc*ttimetauMp*(viscs2_(0, 1, vi)*conv_r_(0, 0, ui)
                                                                +
                                                                viscs2_(1, 1, vi)*conv_r_(1, 0, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
              estif(vi*4 + 1, ui*4 + 1) += 2.0*visc*ttimetauMp*(viscs2_(0, 1, vi)*conv_r_(0, 1, ui)
                                                                +
                                                                viscs2_(1, 1, vi)*conv_r_(1, 1, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
              estif(vi*4 + 1, ui*4 + 2) += 2.0*visc*ttimetauMp*(viscs2_(0, 1, vi)*conv_r_(0, 2, ui)
                                                                +
                                                                viscs2_(1, 1, vi)*conv_r_(1, 2, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
              estif(vi*4 + 2, ui*4)     += 2.0*visc*ttimetauMp*(viscs2_(0, 2, vi)*conv_r_(0, 0, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(1, 0, ui)
                                                                +
                                                                viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
              estif(vi*4 + 2, ui*4 + 1) += 2.0*visc*ttimetauMp*(viscs2_(0, 2, vi)*conv_r_(0, 1, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(1, 1, ui)
                                                                +
                                                                viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
              estif(vi*4 + 2, ui*4 + 2) += 2.0*visc*ttimetauMp*(viscs2_(0, 2, vi)*conv_r_(0, 2, ui)
                                                                +
                                                                viscs2_(1, 2, vi)*conv_r_(1, 2, ui)
                                                                +
                                                                viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;
            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<iel_; ++vi)
        {

          /* viscous stabilisation */
          eforce(vi*4)     += -2.0*visc*timetauMp*(res_old_(0)*viscs2_(0, 0, vi)
                                                   +
                                                   res_old_(1)*viscs2_(0, 1, vi)
                                                   +
                                                   res_old_(2)*viscs2_(0, 2, vi)) ;
          eforce(vi*4 + 1) += -2.0*visc*timetauMp*(res_old_(0)*viscs2_(0, 1, vi)
                                                   +
                                                   res_old_(1)*viscs2_(1, 1, vi)
                                                   +
                                                   res_old_(2)*viscs2_(1, 2, vi)) ;
          eforce(vi*4 + 2) += -2.0*visc*timetauMp*(res_old_(0)*viscs2_(0, 2, vi)
                                                   +
                                                   res_old_(1)*viscs2_(1, 2, vi)
                                                   +
                                                   res_old_(2)*viscs2_(2, 2, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if(cstab)
      {
        const double timefac_timefac_tau_C=timefac*timefac*tau_C;
        const double timefac_timefac_tau_C_divunp=timefac_timefac_tau_C*(vderxy_(0, 0)+vderxy_(1, 1)+vderxy_(2, 2));

        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
            /* continuity stabilisation */
            /*
                       /                        \
                      |                          |
                      | nabla o Du  , nabla o v  |
                      |                          |
                       \                        /
            */
            estif(vi*4, ui*4)         += timefac_timefac_tau_C*derxy_(0, ui)*derxy_(0, vi) ;
            estif(vi*4, ui*4 + 1)     += timefac_timefac_tau_C*derxy_(0, vi)*derxy_(1, ui) ;
            estif(vi*4, ui*4 + 2)     += timefac_timefac_tau_C*derxy_(0, vi)*derxy_(2, ui) ;
            estif(vi*4 + 1, ui*4)     += timefac_timefac_tau_C*derxy_(0, ui)*derxy_(1, vi) ;
            estif(vi*4 + 1, ui*4 + 1) += timefac_timefac_tau_C*derxy_(1, ui)*derxy_(1, vi) ;
            estif(vi*4 + 1, ui*4 + 2) += timefac_timefac_tau_C*derxy_(1, vi)*derxy_(2, ui) ;
            estif(vi*4 + 2, ui*4)     += timefac_timefac_tau_C*derxy_(0, ui)*derxy_(2, vi) ;
            estif(vi*4 + 2, ui*4 + 1) += timefac_timefac_tau_C*derxy_(1, ui)*derxy_(2, vi) ;
            estif(vi*4 + 2, ui*4 + 2) += timefac_timefac_tau_C*derxy_(2, ui)*derxy_(2, vi) ;
          }
        }

        for (int vi=0; vi<iel_; ++vi)
        {
          /* Kontinuit�tsstabilisierung */
          eforce(vi*4)     += -timefac_timefac_tau_C_divunp*derxy_(0, vi) ;
          eforce(vi*4 + 1) += -timefac_timefac_tau_C_divunp*derxy_(1, vi) ;
          eforce(vi*4 + 2) += -timefac_timefac_tau_C_divunp*derxy_(2, vi) ;
        }

      }

      if(fssgv > 0)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          for (int vi=0; vi<iel_; ++vi)
          {
          /* subgrid-viscosity term */
          /*
                        /                        \
                       |       /  \         / \   |
              nu_art * |  eps | Du | , eps | v |  |
                       |       \  /         \ /   |
                        \                        /
          */
          esv(vi*4, ui*4)         += vartfac*(2.0*derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        derxy_(2, ui)*derxy_(2, vi)) ;
          esv(vi*4, ui*4 + 1)     += vartfac*derxy_(0, ui)*derxy_(1, vi) ;
          esv(vi*4, ui*4 + 2)     += vartfac*derxy_(0, ui)*derxy_(2, vi) ;
          esv(vi*4 + 1, ui*4)     += vartfac*derxy_(0, vi)*derxy_(1, ui) ;
          esv(vi*4 + 1, ui*4 + 1) += vartfac*(derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        2.0*derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        derxy_(2, ui)*derxy_(2, vi)) ;
          esv(vi*4 + 1, ui*4 + 2) += vartfac*derxy_(1, ui)*derxy_(2, vi) ;
          esv(vi*4 + 2, ui*4)     += vartfac*derxy_(0, vi)*derxy_(2, ui) ;
          esv(vi*4 + 2, ui*4 + 1) += vartfac*derxy_(1, vi)*derxy_(2, ui) ;
          esv(vi*4 + 2, ui*4 + 2) += vartfac*(derxy_(0, ui)*derxy_(0, vi)
                                                        +
                                                        derxy_(1, ui)*derxy_(1, vi)
                                                        +
                                                        2.0*derxy_(2, ui)*derxy_(2, vi)) ;
          }
        }
      }
    }
  }
  } // end loop over integration cells
  return;
}



//
// calculate stabilization parameter
//
void DRT::Elements::XFluid3Impl::Caltau(
  XFluid3* ele,
  const blitz::Array<double,2>&           evelnp,
  const DRT::Element::DiscretizationType  distype,
  const double                            visc,
  const double                            timefac,
  const int                               fssgv
  )
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // use one point gauss rule to calculate tau at element center
  DRT::Utils::GaussRule3D integrationrule_stabili=DRT::Utils::intrule3D_undefined;
  switch (distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    integrationrule_stabili = DRT::Utils::intrule_hex_1point;
    break;
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    integrationrule_stabili = DRT::Utils::intrule_tet_1point;
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    integrationrule_stabili = DRT::Utils::intrule_wedge_1point;
    break;
  case DRT::Element::pyramid5:
    integrationrule_stabili = DRT::Utils::intrule_pyramid_1point;
    break;
  default:
    dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::Utils::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  DRT::Utils::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
  case DRT::Element::tet4:
  case DRT::Element::pyramid5:
  case DRT::Element::hex8:
  case DRT::Element::wedge6:
    mk = 0.333333333333333333333;
    break;
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
    mk = 0.083333333333333333333;
    break;
  default:
    dserror("type unknown!\n");
  }

  // get velocities at element center
  velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

  // get Jacobian matrix and determinant
  xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  const double vol = wquad*det;

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);

  // inverse of jacobian
  xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
  xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
  xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
  xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
  xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
  xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
  xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
  xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
  xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;

  // compute global derivates
  derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

  // get velocity norm
  const double vel_norm = sqrt(blitz::sum(velint_*velint_));

  // normed velocity at element centre
  if (vel_norm>=1e-6)
  {
    velino_ = velint_/vel_norm;
  }
  else
  {
    velino_ = 0.;
    velino_(0) = 1;
  }

  // get streamlength
  const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  const double strle = 2.0/val;

  // calculate tau

        /*----------------------------------------------------- compute tau_Mu ---*/
        /* stability parameter definition according to

                  Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
                  element method for a generalized Stokes problem. Numerische
                  Mathematik, Vol. 92, pp. 652-677, 2002.
                  http://www.lncc.br/~valentin/publication.htm
        and:
                  Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
                  Finite Element Method for the Advective-Reactive-Diffusive
                  Equation. Computer Methods in Applied Mechanics and Enginnering,
                  Vol. 190, pp. 1785-1800, 2000.
                  http://www.lncc.br/~valentin/publication.htm                   */


    const double re1 =/* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
    const double re2 = mk * vel_norm * strle / /* *1.0 */(2.0 * visc);    /* convective : viscous forces */

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

    // compute tau_Mp
    //    stability parameter definition according to Franca and Valentin (2000)
    //                                       and Barrenechea and Valentin (2002)
    const double re_viscous = /* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
    const double re_convect = mk * vel_norm * hk / /* *1.0 */(2.0 * visc);     /* convective : viscous forces */

    const double xi_viscous = DMAX(re_viscous,1.0);
    const double xi_convect = DMAX(re_convect,1.0);

    /*
                  xi1,xi2 ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re1,re2
                              1
    */
    tau_(1) = DSQR(hk) / (DSQR(hk) * xi_viscous + (/* 2.0*/ 4.0 * timefac * visc/mk) * xi_convect);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    //tau[2] = sqrt(DSQR(visc)+DSQR(0.5*vel_norm*hk));

    // Wall Diss. 99
    /*
                      xi2 ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> Re2
                              1
    */
    const double xi_tau_c = DMIN(re2,1.0);
    tau_(2) = vel_norm * hk * 0.5 * xi_tau_c /timefac;

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv == 1)
  {
    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = 2.0 * re_convect;  /* convective : viscous forces */
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(vel_norm))/(2.0*visc*xi);

  }
  else if (fssgv == 2)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                               +-                                 -+ 1
    //                           2   |          / h \           / h \    | -
    //    visc          = (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent              |          \   / ij        \   / ij |
    //                               +-                                 -+
    //                               |                                   |
    //                               +-----------------------------------+
    //                                    'resolved' rate of strain
    //

    double rateofstrain = 0.0;
    {
      // get velocity (np,i) derivatives at element center
      vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

      blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());
      epsilon = 0.5 * ( vderxy_(i,j) + vderxy_(j,i) );

      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<3;mm++)
        {
          rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
        }
      }
      rateofstrain *= 2.0;
      rateofstrain = sqrt(rateofstrain);
    }
    //
    // Choices of the Smagorinsky constant Cs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)

    const double Cs = 0.17;

    vart_ = Cs * Cs * hk * hk * rateofstrain;
  }
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Impl::BodyForce(XFluid3* ele, const double time)
{
  vector<DRT::Condition*> myneumcond;
  DRT::Node** nodes = ele->Nodes();

  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<iel_;inode++)
  {
    nodes[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  if (nodecount == iel_)
  {
    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
        //curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;isd++)
      {
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_ = 0.;
  }
}


/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y,z at point r,s,t
 |                                            (private)      gammi 07/07
 |
 | From the six equations
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ------ = -- | --*-- + --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy   ds dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dt^2     dt | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | dt dr     dt | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dt     ds | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 | the matrix (jacobian-bar matrix) system
 |
 | +-                                                                                         -+   +-    -+
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
 | |                                                                                           | * |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
 | +-                                                                                         -+   +-    -+
 |
 |                  +-    -+     +-                           -+
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
 |              =   |      |  -  |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drds |     | drds dx   drds dy   drds dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drdt |     | drdt dx   drdt dy   drdt dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dtds |     | dtds dx   dtds dy   dtds dz |
 |                  +-    -+     +-                           -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Impl::gder2(XFluid3* ele)
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // initialize and zero out everything
  static Epetra_SerialDenseMatrix bm(6,6);

  // calculate elements of jacobian_bar matrix
  bm(0,0) = xjm_(0,0)*xjm_(0,0);
  bm(1,0) = xjm_(1,0)*xjm_(1,0);
  bm(2,0) = xjm_(2,0)*xjm_(2,0);
  bm(3,0) = xjm_(0,0)*xjm_(1,0);
  bm(4,0) = xjm_(0,0)*xjm_(2,0);
  bm(5,0) = xjm_(2,0)*xjm_(1,0);

  bm(0,1) = xjm_(0,1)*xjm_(0,1);
  bm(1,1) = xjm_(1,1)*xjm_(1,1);
  bm(2,1) = xjm_(2,1)*xjm_(2,1);
  bm(3,1) = xjm_(0,1)*xjm_(1,1);
  bm(4,1) = xjm_(0,1)*xjm_(2,1);
  bm(5,1) = xjm_(2,1)*xjm_(1,1);

  bm(0,2) = xjm_(0,2)*xjm_(0,2);
  bm(1,2) = xjm_(1,2)*xjm_(1,2);
  bm(2,2) = xjm_(2,2)*xjm_(2,2);
  bm(3,2) = xjm_(0,2)*xjm_(1,2);
  bm(4,2) = xjm_(0,2)*xjm_(2,2);
  bm(5,2) = xjm_(2,2)*xjm_(1,2);

  bm(0,3) = 2.*xjm_(0,0)*xjm_(0,1);
  bm(1,3) = 2.*xjm_(1,0)*xjm_(1,1);
  bm(2,3) = 2.*xjm_(2,0)*xjm_(2,1);
  bm(3,3) = xjm_(0,0)*xjm_(1,1)+xjm_(1,0)*xjm_(0,1);
  bm(4,3) = xjm_(0,0)*xjm_(2,1)+xjm_(2,0)*xjm_(0,1);
  bm(5,3) = xjm_(1,0)*xjm_(2,1)+xjm_(2,0)*xjm_(1,1);

  bm(0,4) = 2.*xjm_(0,0)*xjm_(0,2);
  bm(1,4) = 2.*xjm_(1,0)*xjm_(1,2);
  bm(2,4) = 2.*xjm_(2,0)*xjm_(2,2);
  bm(3,4) = xjm_(0,0)*xjm_(1,2)+xjm_(1,0)*xjm_(0,2);
  bm(4,4) = xjm_(0,0)*xjm_(2,2)+xjm_(2,0)*xjm_(0,2);
  bm(5,4) = xjm_(1,0)*xjm_(2,2)+xjm_(2,0)*xjm_(1,2);

  bm(0,5) = 2.*xjm_(0,1)*xjm_(0,2);
  bm(1,5) = 2.*xjm_(1,1)*xjm_(1,2);
  bm(2,5) = 2.*xjm_(2,1)*xjm_(2,2);
  bm(3,5) = xjm_(0,1)*xjm_(1,2)+xjm_(1,1)*xjm_(0,2);
  bm(4,5) = xjm_(0,1)*xjm_(2,2)+xjm_(2,1)*xjm_(0,2);
  bm(5,5) = xjm_(1,1)*xjm_(2,2)+xjm_(2,1)*xjm_(1,2);

  /*------------------ determine 2nd derivatives of coord.-functions */

  /*
  |
  |         0 1 2              0...iel-1
  |        +-+-+-+             +-+-+-+-+        0 1 2
  |        | | | | 0           | | | | | 0     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | 0
  |        | | | | 1           | | | | | 1   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 2           | | | | | 2     +-+-+-+
  |        +-+-+-+       =     +-+-+-+-+       | | | | .
  |        | | | | 3           | | | | | 3     +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 4           | | | | | 4   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 5           | | | | | 5     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | iel-1
  |		     	      	     	       +-+-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-                  -+
  |  	   	    	    	        | d^2x   d^2y   d^2z |
  |  	   	    	    	        | ----   ----   ---- |
  | 	   	   	   	        | dr^2   dr^2   dr^2 |
  | 	   	   	   	        |                    |
  | 	   	   	   	        | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  | 	   	   	   	        | ds^2   ds^2   ds^2 |
  | 	   	   	   	        |                    |
  | 	   	   	   	        | d^2x   d^2y   d^2z |
  | 	   	   	   	        | ----   ----   ---- |
  | 	   	   	   	        | dt^2   dt^2   dt^2 |
  |               yields    xder2  =    |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drds   drds   drds |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drdt   drdt   drdt |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dsdt   dsdt   dsdt |
  | 	   	   	   	        +-                  -+
  |
  |
  */

  xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);

  /*
  |        0...iel-1             0 1 2
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 0          | | | | 0
  |        +-+-+-+-+            +-+-+-+            0...iel-1
  |        | | | | | 1          | | | | 1         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 0
  |        | | | | | 2          | | | | 2         +-+-+-+-+
  |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
  |        | | | | | 3          | | | | 3         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 2
  |        | | | | | 4          | | | | 4         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 5          | | | | 5          derxy
  |        +-+-+-+-+            +-+-+-+
  |
  |       chainrulerhs          xder2
  */

  derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);

  /*
  |        0...iel-1            0...iel-1         0...iel-1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 0          | | | | | 0       | | | | | 0
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 1          | | | | | 1       | | | | | 1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 2          | | | | | 2       | | | | | 2
  |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
  |        | | | | | 3          | | | | | 3       | | | | | 3
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 4          | | | | | 4       | | | | | 4
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 5          | | | | | 5       | | | | | 5
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |
  |       chainrulerhs         chainrulerhs        deriv2
  */

  derxy2_ += deriv2_;

  /* make LR decomposition and solve system for all right hand sides
   * (i.e. the components of chainrulerhs)
  |
  |          0  1  2  3  4  5         i        i
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 0     | | 0    | | 0
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 1     | | 1    | | 1
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 2     | | 2    | | 2
  | 	   +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
  |        |  |  |  |  |  |  | 3     | | 3    | | 3
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 4     | | 4    | | 4
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 5     | | 5    | | 5
  |        +--+--+--+--+--+--+       +-+      +-+
  |                                   |        |
  |                                   |        |
  |                                   derxy2[i]|
  |		                               |
  |		                               chainrulerhs[i]
  |
  |	  yields
  |
  |                      0...iel-1
  |                      +-+-+-+-+
  |                      | | | | | 0 = drdr
  |                      +-+-+-+-+
  |                      | | | | | 1 = dsds
  |                      +-+-+-+-+
  |                      | | | | | 2 = dtdt
  |            derxy2 =  +-+-+-+-+
  |                      | | | | | 3 = drds
  |                      +-+-+-+-+
  |                      | | | | | 4 = drdt
  |                      +-+-+-+-+
  |                      | | | | | 5 = dsdt
  |    	          	 +-+-+-+-+
  */

  Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,iel_);

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(ederxy2,ederxy2);
  solver.Solve();

  return;
}


#endif
#endif
