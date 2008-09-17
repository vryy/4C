/*----------------------------------------------------------------------*/
/*!
\file fluid2_impl.cpp

\brief Internal implementation of Fluid2 element (one-step-theta/BDF2)

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_LAPACK.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Impl* DRT::ELEMENTS::Fluid2Impl::Impl(DRT::ELEMENTS::Fluid2* f2)
{
  switch (f2->NumNode())
  {
  case 4:
  {
    static Fluid2Impl* f4;
    if (f4==NULL)
      f4 = new Fluid2Impl(4);
    return f4;
  }
  case 8:
  {
    static Fluid2Impl* f8;
    if (f8==NULL)
      f8 = new Fluid2Impl(8);
    return f8;
  }
  case 9:
  {
    static Fluid2Impl* f9;
    if (f9==NULL)
      f9 = new Fluid2Impl(9);
    return f9;
  }
  case 3:
  {
    static Fluid2Impl* f3;
    if (f3==NULL)
      f3 = new Fluid2Impl(3);
    return f3;
  }
  case 6:
  {
    static Fluid2Impl* f6;
    if (f6==NULL)
      f6 = new Fluid2Impl(6);
    return f6;
  }
  default:
    dserror("node number %d not supported", f2->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Impl::Fluid2Impl(int iel)
  : iel_(iel),
    vart_(),
    xyze_(2,iel_,blitz::ColumnMajorArray<2>()),
    edeadng_(2,iel_,blitz::ColumnMajorArray<2>()),
    funct_(iel_),
    densfunct_(iel_),
    functdens_(iel_),
    deriv_(2,iel_,blitz::ColumnMajorArray<2>()),
    deriv2_(3,iel_,blitz::ColumnMajorArray<2>()),
    xjm_(2,2,blitz::ColumnMajorArray<2>()),
    xji_(2,2,blitz::ColumnMajorArray<2>()),
    vderxy_(2,2,blitz::ColumnMajorArray<2>()),
    mderxy_(2,2,blitz::ColumnMajorArray<2>()),
    fsvderxy_(2,2,blitz::ColumnMajorArray<2>()),
    vderxy2_(2,3,blitz::ColumnMajorArray<2>()),
    derxy_(2,iel_,blitz::ColumnMajorArray<2>()),
    densderxy_(2,iel_,blitz::ColumnMajorArray<2>()),
    derxy2_(3,iel_,blitz::ColumnMajorArray<2>()),
    bodyforce_(2),
    histmom_(2),
    histcon_(),
    //velino_(2),
    velint_(2),
    fsvelint_(2),
    convvelint_(2),
    gradp_(2),
    tau_(3),
    viscs2_(2,2,iel_,blitz::ColumnMajorArray<3>()),
    conv_c_(iel_),
    mdiv_(),
    vdiv_(),
    rhsmom_(2),
    rhscon_(),
    conv_old_(2),
    visc_old_(2),
    res_old_(2),
    conv_resM_(iel_),
    xder2_(3,2,blitz::ColumnMajorArray<2>()),
    vderiv_(2,2,blitz::ColumnMajorArray<2>())
{
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                  vg 08/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Impl::Sysmat(
  Fluid2*                                 ele,
  const blitz::Array<double,2>&           evelnp,
  const blitz::Array<double,2>&           fsevelnp,
  const blitz::Array<double,1>&           eprenp,
  const blitz::Array<double,1>&           edensnp,
  const blitz::Array<double,2>&           emhist,
  const blitz::Array<double,1>&           echist,
  const blitz::Array<double,2>&           edispnp,
  const blitz::Array<double,2>&           egridv,
  blitz::Array<double,2>&                 estif,
  blitz::Array<double,2>&                 emesh,
  blitz::Array<double,1>&                 eforce,
  struct _MATERIAL*                       material,
  double                                  time,
  double                                  dt,
  double                                  timefac,
  bool                                    newton,
  bool                                    loma,
  const bool                              higher_order_ele,
  const enum Fluid2::StabilisationAction  fssgv,
  const enum Fluid2::StabilisationAction  pspg,
  const enum Fluid2::StabilisationAction  supg,
  const enum Fluid2::StabilisationAction  vstab,
  const enum Fluid2::StabilisationAction  cstab,
  const enum Fluid2::StabilisationAction  cross,
  const enum Fluid2::StabilisationAction  reynolds,
  const enum Fluid2::TauType              whichtau,
  const enum Fluid2::TurbModelAction      turb_mod_action,
  double&                                 Cs,
  double&                                 visceff
  )
{
  // set element data
  const DRT::Element::DiscretizationType distype = ele->Shape();
  const int numnode = iel_;

  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // add displacement, when fluid nodes move in the ALE case
  if (ele->is_ale_) xyze_ += edispnp;

  // dead load in element nodes
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->mattyp != m_fluid
   && material->mattyp != m_carreauyasuda
   && material->mattyp != m_modpowerlaw) dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(material->mattyp == m_fluid) visc = material->m.fluid->viscosity;

  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally
  Caltau(ele,
         evelnp,
         fsevelnp,
         edensnp,
         distype,
         whichtau,
         material,
         visc,
         timefac,
         dt,
         turb_mod_action,
         Cs,
         visceff,
         fssgv);

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == Fluid2::viscous_stab_usfem || vstab == Fluid2::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordiantes of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
      +-       -+ T      +-       -+
      | dx   dx |        | dx   dy |
      | --   -- |        | --   -- |
      | dr   ds |        | dr   dr |
      |         |   =    |         |
      | dy   dy |        | dx   dy |
      | --   -- |        | --   -- |
      | dr   ds |        | ds   ds |
      +-       -+        +-       -+
    */
    xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);

    // The determinant is computed using Sarrus's rule:
    const double det = xjm_(0,0)*xjm_(1,1)-xjm_(0,1)*xjm_(1,0);

    if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

    const double fac = intpoints.qwgt[iquad]*det;

    //--------------------------------------------------------------
    //             compute global first derivates
    //--------------------------------------------------------------
    //
    /*
      Use the Jacobian and the known derivatives in element coordinate
      directions on the right hand side to compute the derivatives in
      global coordinate directions

          +-          -+     +-    -+      +-    -+
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  dr    dr  |     |  dx  |      |  dr  |
          |            |  *  |      |   =  |      | for all k
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  ds    ds  |     |  dy  |      |  ds  |
          +-          -+     +-    -+      +-    -+

          Matrix is inverted analytically
    */
    // inverse of jacobian
    xji_(0,0) = ( xjm_(1,1))/det;
    xji_(0,1) = (-xjm_(0,1))/det;
    xji_(1,0) = (-xjm_(1,0))/det;
    xji_(1,1) = ( xjm_(0,0))/det;

    // compute global derivates
    derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<numnode; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);
      functdens_(inode) = funct_(inode)/edensnp(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    //--------------------------------------------------------------
    //             compute second global derivative
    //--------------------------------------------------------------

    /*----------------------------------------------------------------------*
     |  calculate second global derivatives w.r.t. x,y at point r,s
     |                                            (private)      gammi 02/08
     |
     | From the three equations
     |
     |              +-             -+
     |  d^2N     d  | dx dN   dy dN |
     |  ----   = -- | --*-- + --*-- |
     |  dr^2     dr | dr dx   dr dy |
     |              +-             -+
     |
     |              +-             -+
     |  d^2N     d  | dx dN   dy dN |
     |  ------ = -- | --*-- + --*-- |
     |  ds^2     ds | ds dx   ds dy |
     |              +-             -+
     |
     |              +-             -+
     |  d^2N     d  | dx dN   dy dN |
     | -----   = -- | --*-- + --*-- |
     | ds dr     ds | dr dx   dr dy |
     |              +-             -+
     |
     | the matrix (jacobian-bar matrix) system
     |
     | +-                                          -+   +-    -+
     | |   /dx\^2        /dy\^2         dy dx       |   | d^2N |
     | |  | -- |        | ---|        2*--*--       |   | ---- |
     | |   \dr/          \dr/           dr dr       |   | dx^2 |
     | |                                            |   |      |
     | |   /dx\^2        /dy\^2         dy dx       |   | d^2N |
     | |  | -- |        | ---|        2*--*--       |   | ---- |
     | |   \ds/          \ds/           ds ds       |   | dy^2 |
     | |                                            | * |      |
     | |   dx dx         dy dy      dx dy   dx dy   |   | d^2N |
     | |   --*--         --*--      --*-- + --*--   |   | ---- |
     | |   dr ds         dr ds      dr ds   ds dr   |   | dxdy |
     | +-                                          -+   +-    -+
     |
     |                  +-    -+     +-                 -+
     |                  | d^2N |     | d^2x dN   d^2y dN |
     |                  | ---- |     | ----*-- + ----*-- |
     |                  | dr^2 |     | dr^2 dx   dr^2 dy |
     |                  |      |     |                   |
     |                  | d^2N |     | d^2x dN   d^2y dN |
     |              =   | ---- |  -  | ----*-- + ----*-- |
     |                  | ds^2 |     | ds^2 dx   ds^2 dy |
     |                  |      |     |                   |
     |                  | d^2N |     | d^2x dN   d^2y dN |
     |                  | ---- |     | ----*-- + ----*-- |
     |                  | drds |     | drds dx   drds dy |
     |                  +-    -+     +-                 -+
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
    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      // initialize and zero out everything
      blitz::Array<double,2> bm(3,3,blitz::ColumnMajorArray<2>());

      // calculate elements of jacobian_bar matrix
      bm(0,0) =                     xjm_(0,0)*xjm_(0,0);
      bm(0,1) =                     xjm_(0,1)*xjm_(0,1);
      bm(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

      bm(1,0) =                     xjm_(1,0)*xjm_(1,0);
      bm(1,1) =                     xjm_(1,1)*xjm_(1,1);
      bm(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

      bm(2,0) =                     xjm_(0,0)*xjm_(1,0);
      bm(2,1) =                     xjm_(0,1)*xjm_(1,1);
      bm(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);


      /*------------------ determine 2nd derivatives of coord.-functions */
      /*
       |                                             0 1
       |         0 1              0...iel-1         +-+-+
       |        +-+-+             +-+-+-+-+         | | | 0
       |        | | | 0           | | | | | 0       +-+-+
       |        +-+-+             +-+-+-+-+         | | | .
       |        | | | 1     =     | | | | | 1     * +-+-+ .
       |        +-+-+             +-+-+-+-+         | | | .
       |        | | | 2           | | | | | 2       +-+-+
       |        +-+-+             +-+-+-+-+         | | | iel-1
       |                                            +-+-+
       |
       |        xder2               deriv2          xyze^T
       |
       |
       |                                        +-           -+
       |  	   	    	    	        | d^2x   d^2y |
       |  	   	    	    	        | ----   ---- |
       | 	   	   	   	        | dr^2   dr^2 |
       | 	   	   	   	        |             |
       | 	   	   	   	        | d^2x   d^2y |
       |                    yields    xder2  =  | ----   ---- |
       | 	   	   	   	        | ds^2   ds^2 |
       | 	   	   	   	        |             |
       | 	   	   	   	        | d^2x   d^2y |
       | 	   	   	   	        | ----   ---- |
       | 	   	   	   	        | drds   drds |
       | 	   	   	   	        +-           -+
      */
      xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);


      /*
       |        0...iel-1             0 1
       |        +-+-+-+-+            +-+-+               0...iel-1
       |        | | | | | 0          | | | 0             +-+-+-+-+
       |        +-+-+-+-+            +-+-+               | | | | | 0
       |        | | | | | 1     =    | | | 1     *       +-+-+-+-+   * (-1)
       |        +-+-+-+-+            +-+-+               | | | | | 1
       |        | | | | | 2          | | | 2             +-+-+-+-+
       |        +-+-+-+-+            +-+-+
       |
       |       chainrulerhs          xder2                 derxy
      */
      derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);

      /*
       |        0...iel-1             0...iel-1             0...iel-1
       |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
       |        | | | | | 0           | | | | | 0           | | | | | 0
       |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
       |        | | | | | 1     =     | | | | | 1     +     | | | | | 1
       |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
       |        | | | | | 2           | | | | | 2           | | | | | 2
       |        +-+-+-+-+             +-+-+-+-+             +-+-+-+-+
       |
       |       chainrulerhs          chainrulerhs             deriv2
      */
      derxy2_ += deriv2_;

       /* make LU decomposition and solve system for all right hand sides
       * (i.e. the components of chainrulerhs)

       |
       |            0  1  2          i        i
       | 	   +--+--+--+       +-+      +-+
       | 	   |  |  |  | 0     | | 0    | | 0
       | 	   +--+--+--+       +-+	     +-+
       | 	   |  |  |  | 1  *  | | 1 =  | | 1  for i=0...iel-1
       | 	   +--+--+--+       +-+	     +-+
       | 	   |  |  |  | 2     | | 2    | | 2
       | 	   +--+--+--+       +-+	     +-+
       |                             |        |
       |                             |        |
       |                           derxy2[i]  |
       |                                      |
       |                                chainrulerhs[i]
       |
       |
       |
       |                      0...iel-1
       |		     +-+-+-+-+
       |		     | | | | | 0
       |		     +-+-+-+-+
       |	  yields     | | | | | 1
       |		     +-+-+-+-+
       |                     | | | | | 2
       | 		     +-+-+-+-+
       |
       |                      derxy2
       |
       */

      // Use LAPACK
      Epetra_LAPACK          solver;

      // a vector specifying the pivots (reordering)
      int pivot[3];

      // error code
      int ierr = 0;

      // Perform LU factorisation --- this call replaces bm with its factorisation
      solver.GETRF(3,3,bm.data(),3,&(pivot[0]),&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform LU factorisation during computation of derxy2");
      }

      // backward substitution. GETRS replaces the input (chainrulerhs, currently
      // stored on derxy2) with the result
      solver.GETRS('N',3,iel_,bm.data(),3,&(pivot[0]),derxy2_.data(),3,&ierr);

      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");

      // calculate 2nd velocity derivatives at integration point
      vderxy2_ = blitz::sum(derxy2_(j,k)*evelnp(i,k),k);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get momentum (n+g,i) at integration point
    velint_ = blitz::sum(densfunct_(j)*evelnp(i,j),j);

    // get history data (n,i) at integration point
    histmom_ = blitz::sum(densfunct_(j)*emhist(i,j),j);
    histcon_ = blitz::sum(funct_*echist);

    // get velocity (np,i) derivatives at integration point
    vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

    // get momentum (np,i) derivatives at integration point
    mderxy_ = blitz::sum(densderxy_(j,k)*evelnp(i,k),k);

    // get fine-scale velocity (np,i) derivatives at integration point
    if (fssgv != Fluid2::fssgv_no) fsvderxy_ = blitz::sum(derxy_(j,k)*fsevelnp(i,k),k);
    else fsvderxy_ = 0.;

    // get convective velocity at integration point
    // We handle the ale case very implicitely here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    if (ele->is_ale_) convvelint_ = velint_ - blitz::sum(densfunct_(j)*egridv(i,j),j);
    else              convvelint_ = velint_;

    // get pressure gradient at integration point
    gradp_ = blitz::sum(derxy_(i,j)*eprenp(j),j);

    // get pressure at integration point
    double press = blitz::sum(funct_*eprenp);

    // get density at integration point
    double dens = blitz::sum(funct_*edensnp);

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_ = blitz::sum(edeadng_(i,j)*densfunct_(j),j);

    //--------------------------------------------------------------
    // perform integration for entire matrix and rhs
    //--------------------------------------------------------------

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

    // subgrid-viscosity factor
    const double vartfac = vart_*timefacfac;

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsmom_ = histmom_(i) + bodyforce_(i)*timefac;
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_ = blitz::sum(vderxy_(i, j)*convvelint_(j), j);

    /* Viscous term  div epsilon(u_old) */
    if (higher_order_ele)
    {
      visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,2));
      visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,2));
    }
    else
    {
      visc_old_ = 0;
    }

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_c_ = blitz::sum(derxy_(j,i)*convvelint_(j), j);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
           1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
           - |                              |         N_y .. y-line of N
           2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
             \                              /                            */

      for (int _=0; _<numnode; ++_) viscs2_(0,0,_) = 0.5 * (2.0 * derxy2_(0,_) + derxy2_(1,_));
      for (int _=0; _<numnode; ++_) viscs2_(0,1,_) = 0.5 *  derxy2_(2,_);
      for (int _=0; _<numnode; ++_) viscs2_(1,0,_) = 0.5 *  derxy2_(2,_);
      for (int _=0; _<numnode; ++_) viscs2_(1,1,_) = 0.5 * (derxy2_(0,_) + 2.0 * derxy2_(1,_));

      if (loma)
      {
        /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
        /*    /                   \
            1 |  N_x,xx + N_y,yx  |
         -  - |                   |
            3 |  N_x,xy + N_y,yy  |
              \                   /

                 with N_x .. x-line of N
                 N_y .. y-line of N                                      */

        for (int _=0; _<numnode; ++_) viscs2_(0,0,_) -=  derxy2_(0,_)/3.0;
        for (int _=0; _<numnode; ++_) viscs2_(0,1,_) -=  derxy2_(2,_)/3.0;
        for (int _=0; _<numnode; ++_) viscs2_(1,0,_) -=  derxy2_(2,_)/3.0;
        for (int _=0; _<numnode; ++_) viscs2_(1,1,_) -=  derxy2_(1,_)/3.0;
      }
    }
    else
    {
      viscs2_ = 0;
    }

    /* momentum and velocity divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    // evaluate residual once for all stabilization right hand sides
    res_old_ = velint_-rhsmom_+timefac*(conv_old_+gradp_-2*visceff*visc_old_);

    /*
      This is the operator

                  /               \
                 | resM    o nabla |
                  \    (i)        /

                  required for (lhs) cross- and (rhs) Reynolds-stress calculation

    */

    if (cross    == Fluid2::cross_stress_stab ||
        reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      conv_resM_ =  blitz::sum(res_old_(j)*densderxy_(j,i),j);

    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      for (int ui=0; ui<numnode; ++ui)
      {
        double v = fac*densfunct_(ui)
#if 1
                   + timefacfac*conv_c_(ui)
#endif
                   ;
        for (int vi=0; vi<numnode; ++vi)
        {

          /* inertia (contribution to mass matrix) */
          /*

          /          \
          |          |
          |  Du , v  |
          |          |
          \          /
          */

          /* convection, convective part */
          /*

          /                         \
          |  / n+1       \          |
          | | u   o nabla | Du , v  |
          |  \ (i)       /          |
          \                         /

          */
          double v2 = v*funct_(vi) ;
          estif(vi*3,     ui*3    ) += v2;
          estif(vi*3 + 1, ui*3 + 1) += v2;
        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        for (int vi=0; vi<numnode; ++vi)
        {

          /* viscosity term */
          /*

                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(vi*3,     ui*3    ) += visceff*timefacfac*(2.0*derxy_(0, ui)*derxy_(0, vi)
                                                           +
                                                           derxy_(1, ui)*derxy_(1, vi)) ;
          estif(vi*3,     ui*3 + 1) += visceff*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*3 + 1, ui*3    ) += visceff*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*3 + 1, ui*3 + 1) += visceff*timefacfac*(derxy_(0, ui)*derxy_(0, vi)
                                                           +
                                                           2.0*derxy_(1, ui)*derxy_(1, vi)) ;

        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          /* Druckterm */
          /*

          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */

          estif(vi*3,     ui*3 + 2) += v*derxy_(0, vi) ;
          estif(vi*3 + 1, ui*3 + 2) += v*derxy_(1, vi) ;

        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        double v = timefacfac*functdens_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {

          /* Divergenzfreiheit */
          /*
            /                  \
            |                  |
            | nabla o Du  , q  |
            |                  |
            \                  /
          */
          estif(vi*3 + 2, ui*3    ) += v*densderxy_(0, ui) ;
          estif(vi*3 + 2, ui*3 + 1) += v*densderxy_(1, ui) ;
        }
      }

      if (newton)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          double v = timefacfac*funct_(vi);
          for (int ui=0; ui<numnode; ++ui)
          {

            /*  convection, reactive part

            /                           \
            |  /          \   n+1       |
            | | Du o nabla | u     , v  |
            |  \          /   (i)       |
            \                           /
            */
            estif(vi*3,     ui*3    ) += v*vderxy_(0, 0)*densfunct_(ui) ;
            estif(vi*3,     ui*3 + 1) += v*vderxy_(0, 1)*densfunct_(ui) ;
            estif(vi*3 + 1, ui*3    ) += v*vderxy_(1, 0)*densfunct_(ui) ;
            estif(vi*3 + 1, ui*3 + 1) += v*vderxy_(1, 1)*densfunct_(ui) ;
          }
        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        /* inertia */
        double v = -fac*funct_(vi);
        eforce(vi*3    ) += v*velint_(0) ;
        eforce(vi*3 + 1) += v*velint_(1) ;
      }

#if 1
      for (int vi=0; vi<numnode; ++vi)
      {
        /* convection */
        double v = -timefacfac*funct_(vi);
        eforce(vi*3    ) += v*(convvelint_(0)*vderxy_(0, 0)
                               +
                               convvelint_(1)*vderxy_(0, 1)) ;
        eforce(vi*3 + 1) += v*(convvelint_(0)*vderxy_(1, 0)
                               +
                               convvelint_(1)*vderxy_(1, 1)) ;
      }
#endif

      for (int vi=0; vi<numnode; ++vi)
      {
        /* pressure */
        double v = press*timefacfac;
        eforce(vi*3    ) += v*derxy_(0, vi) ;
        eforce(vi*3 + 1) += v*derxy_(1, vi) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        /* viscosity */
        eforce(vi*3    ) -= visceff*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                                                +
                                                derxy_(1, vi)*vderxy_(0, 1)
                                                +
                                                derxy_(1, vi)*vderxy_(1, 0)) ;
        eforce(vi*3 + 1) -= visceff*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
                                                +
                                                derxy_(0, vi)*vderxy_(1, 0)
                                                +
                                                2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        // source term of the right hand side
        double v = fac*funct_(vi);
        eforce(vi*3    ) += v*rhsmom_(0) ;
        eforce(vi*3 + 1) += v*rhsmom_(1) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        // continuity equation
        eforce(vi*3 + 2) -= timefacfac*functdens_(vi)*mdiv_ ;
      }

      if (loma)
      {
        double v = -(2.0/3.0)*visceff*timefacfac ;
        for (int ui=0; ui<numnode; ++ui)
        {
          for (int vi=0; vi<numnode; ++vi)
          {

            /* viscosity term - subtraction for low-Mach-number flow */
            /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
            */
            estif(vi*3,     ui*3    ) += v*derxy_(0, vi)*derxy_(0, ui) ;
            estif(vi*3,     ui*3 + 1) += v*derxy_(0, vi)*derxy_(1, ui) ;
            estif(vi*3 + 1, ui*3    ) += v*derxy_(1, vi)*derxy_(0, ui) ;
            estif(vi*3 + 1, ui*3 + 1) += v*derxy_(1, vi)*derxy_(1, ui) ;

          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          /* viscosity term - subtraction for low-Mach-number flow */
          eforce(vi*3    ) -= derxy_(0, vi)*v*vdiv_ ;
          eforce(vi*3 + 1) -= derxy_(1, vi)*v*vdiv_ ;

          /* rhs term of continuity equation */
          eforce(vi*3 + 2) += fac*functdens_(vi)*rhscon_ ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg == Fluid2::pstab_use_pspg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          double v = timetauMp*densfunct_(ui)
#if 1
                     + ttimetauMp*conv_c_(ui)
#endif
                     ;
          for (int vi=0; vi<numnode; ++vi)
          {

            /* pressure stabilisation: inertia */
            /*
              /              \
              |                |
              |  Du , nabla q  |
              |                |
              \              /
            */
            /* pressure stabilisation: convection, convective part */
            /*

            /                            \
            |  / n+1       \               |
            | | u   o nabla | Du , nabla q |
            |  \ (i)       /               |
            \                            /

            */

            estif(vi*3 + 2, ui*3    ) += v*derxy_(0, vi) ;
            estif(vi*3 + 2, ui*3 + 1) += v*derxy_(1, vi) ;
          }
        }

        if (higher_order_ele)
        {
          double v = -2.0*visceff*ttimetauMp;
          for (int ui=0; ui<numnode; ++ui)
          {
            for (int vi=0; vi<numnode; ++vi)
            {

              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(vi*3 + 2, ui*3    ) += v*(derxy_(0, vi)*viscs2_(0, 0, ui)
                                              +
                                              derxy_(1, vi)*viscs2_(0, 1, ui)) ;
              estif(vi*3 + 2, ui*3 + 1) += v*(derxy_(0, vi)*viscs2_(0, 1, ui)
                                              +
                                              derxy_(1, vi)*viscs2_(1, 1, ui)) ;
            }
          }
        }

        for (int ui=0; ui<numnode; ++ui)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            /* pressure stabilisation: pressure( L_pres_p) */
            /*
              /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
              \                    /
            */
            estif(vi*3 + 2, ui*3 + 2) += ttimetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                     +
                                                     derxy_(1, ui)*derxy_(1, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            double v = ttimetauMp*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              /*  pressure stabilisation: convection, reactive part

              /                             \
              |  /          \   n+1           |
              | | Du o nabla | u     , grad q |
              |  \          /   (i)           |
              \                             /

              */
              estif(vi*3 + 2, ui*3    ) += v*(derxy_(0, vi)*vderxy_(0, 0)
                                              +
                                              derxy_(1, vi)*vderxy_(1, 0)) ;
              estif(vi*3 + 2, ui*3 + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                                              +
                                              derxy_(1, vi)*vderxy_(1, 1)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<numnode; ++vi)
        {
          // pressure stabilisation
          eforce(vi*3 + 2) -= timetauMp*(res_old_(0)*derxy_(0, vi)
                                         +
                                         res_old_(1)*derxy_(1, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg == Fluid2::convective_stab_supg)
      {
#if 1
        for (int ui=0; ui<numnode; ++ui)
        {
          double v = timetauM*densfunct_(ui) + ttimetauM*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            /* supg stabilisation: inertia  */
            /*
              /                        \
              |        / n+1       \     |
              |  Du , | u   o nabla | v  |
              |        \ (i)       /     |
              \                        /
            */

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                           \
            |    / n+1        \        / n+1        \     |
            |   | u    o nabla | Du , | u    o nabla | v  |
            |    \ (i)        /        \ (i)        /     |
            \                                           /

            */

            estif(vi*3,     ui*3    ) += v*conv_c_(vi);
            estif(vi*3 + 1, ui*3 + 1) += v*conv_c_(vi);
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          double v = ttimetauM*conv_c_(vi);
          for (int ui=0; ui<numnode; ++ui)
          {

            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                              \
              |              / n+1       \     |
              |  nabla Dp , | u   o nabla | v  |
              |              \ (i)       /     |
              \                              /
            */
            estif(vi*3,     ui*3 + 2) += v*derxy_(0, ui) ;
            estif(vi*3 + 1, ui*3 + 2) += v*derxy_(1, ui) ;
          }
        }

        if (higher_order_ele)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            double v = -2.0*visceff*ttimetauM*conv_c_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {

              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                        \
                |               /  \    / n+1        \     |
                |  nabla o eps | Du |, | u    o nabla | v  |
                |               \  /    \ (i)        /     |
                \                                        /
              */
              estif(vi*3,     ui*3    ) += v*viscs2_(0, 0, ui) ;
              estif(vi*3 + 1, ui*3    ) += v*viscs2_(0, 1, ui) ;

              estif(vi*3,     ui*3 + 1) += v*viscs2_(0, 1, ui) ;
              estif(vi*3 + 1, ui*3 + 1) += v*viscs2_(1, 1, ui) ;
            }
          }
        }
#endif

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            double v = timetauM*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                /                           \
                |   n+1      /          \     |
                |  u      , | Du o nabla | v  |
                |   (i)      \          /     |
                \                           /

              */
              estif(vi*3,     ui*3    ) += v*velint_(0)*derxy_(0, vi) ;
              estif(vi*3 + 1, ui*3    ) += v*velint_(1)*derxy_(0, vi) ;

              estif(vi*3,     ui*3 + 1) += v*velint_(0)*derxy_(1, vi) ;
              estif(vi*3 + 1, ui*3 + 1) += v*velint_(1)*derxy_(1, vi) ;
            }
          }

#if 1
          {
            double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1);
            double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1);

            for (int ui=0; ui<numnode; ++ui)
            {
              double v = ttimetauM*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {

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
                estif(vi*3,     ui*3    ) += (conv_c_(vi)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(vi*3 + 1, ui*3    ) += (conv_c_(vi)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;

                estif(vi*3,     ui*3 + 1) += (conv_c_(vi)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(vi*3 + 1, ui*3 + 1) += (conv_c_(vi)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
              }
            }
          }
#endif

          for (int ui=0; ui<numnode; ++ui)
          {
            double v = ttimetauM*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {

              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                               \
                |         n+1    /          \     |
                |  nabla p    , | Du o nabla | v  |
                |         (i)    \          /     |
                \                               /
              */
              estif(vi*3,     ui*3    ) += v*gradp_(0)*derxy_(0, vi) ;
              estif(vi*3 + 1, ui*3    ) += v*gradp_(1)*derxy_(0, vi) ;

              estif(vi*3,     ui*3 + 1) += v*gradp_(0)*derxy_(1, vi) ;
              estif(vi*3 + 1, ui*3 + 1) += v*gradp_(1)*derxy_(1, vi) ;

            }
          }

          if (higher_order_ele)
          {
            for (int ui=0; ui<numnode; ++ui)
            {
              double v = -2.0*visceff*ttimetauM*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                           \
                  |               / n+1 \    /          \     |
                  |  nabla o eps | u     |, | Du o nabla | v  |
                  |               \ (i) /    \          /     |
                  \                                           /
                */
                estif(vi*3,     ui*3    ) += v*visc_old_(0)*derxy_(0, vi) ;
                estif(vi*3 + 1, ui*3    ) += v*visc_old_(1)*derxy_(0, vi) ;

                estif(vi*3,     ui*3 + 1) += v*visc_old_(0)*derxy_(1, vi) ;
                estif(vi*3 + 1, ui*3 + 1) += v*visc_old_(1)*derxy_(1, vi) ;

              }
            }
          }

          for (int ui=0; ui<numnode; ++ui)
          {
            double v = -timetauM*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {

              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                             \
                |              /          \     |
                |  rhsint   , | Du o nabla | v  |
                |              \          /     |
                \                             /

              */
              estif(vi*3    , ui*3    ) += v*derxy_(0, vi)*rhsmom_(0) ;
              estif(vi*3 + 1, ui*3    ) += v*derxy_(0, vi)*rhsmom_(1) ;

              estif(vi*3    , ui*3 + 1) += v*derxy_(1, vi)*rhsmom_(0) ;
              estif(vi*3 + 1, ui*3 + 1) += v*derxy_(1, vi)*rhsmom_(1) ;

            } // vi
          } // ui
        } // if newton

#if 1
        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<numnode; ++vi)
        {
          // supg stabilisation
          double v = -timetauM*conv_c_(vi);
          eforce(vi*3    ) += v*res_old_(0) ;
          eforce(vi*3 + 1) += v*res_old_(1) ;
        }
#endif
      }

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART

      if (higher_order_ele)
      {
        if(vstab != Fluid2::viscous_stab_none)
        {
          const double two_visc_timefac = vstabfac*2.0*visc*timetauMp;

          // viscous stabilization either on left hand side or on right hand side
          if (vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_usfem)
          {
            const double two_visc_ttimefac = vstabfac*2.0*visc*ttimetauMp;
            const double four_visc2_ttimefac = vstabfac*4.0*visceff*visc*ttimetauMp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<numnode; ++ui)
            {
              double v = two_visc_timefac*densfunct_(ui)
#if 1
                         + two_visc_ttimefac*conv_c_(ui)
#endif
                         ;
              for (int vi=0; vi<numnode; ++vi)
              {
                /* viscous stabilisation, inertia part */
                /*
                  /                    \
                  |                    |
              +/- |  Du , div eps (v)  |
                  |                    |
                  \                    /
                */
                /* viscous stabilisation, convective part */
                /*
                  /                                  \
                  |  / n+1       \                   |
              +/- | | u   o nabla | Du , div eps (v) |
                  |  \ (i)       /                   |
                  \                                  /
                */
                estif(vi*3,     ui*3    ) += v*viscs2_(0, 0, vi) ;
                estif(vi*3 + 1, ui*3    ) += v*viscs2_(0, 1, vi) ;

                estif(vi*3,     ui*3 + 1) += v*viscs2_(0, 1, vi) ;
                estif(vi*3 + 1, ui*3 + 1) += v*viscs2_(1, 1, vi) ;
              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              for (int vi=0; vi<numnode; ++vi)
              {

                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(vi*3,     ui*3 + 2) += two_visc_ttimefac*(derxy_(0, ui)*viscs2_(0, 0, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(0, 1, vi)) ;
                estif(vi*3 + 1, ui*3 + 2) += two_visc_ttimefac*(derxy_(0, ui)*viscs2_(0, 1, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(1, 1, vi)) ;

              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              for (int vi=0; vi<numnode; ++vi)
              {
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(vi*3,     ui*3    ) -= four_visc2_ttimefac*(viscs2_(0,0,ui)*viscs2_(0,0,vi)+viscs2_(0,1,ui)*viscs2_(0,1,vi)) ;
                estif(vi*3 + 1, ui*3    ) -= four_visc2_ttimefac*(viscs2_(0,0,ui)*viscs2_(0,1,vi)+viscs2_(0,1,ui)*viscs2_(1,1,vi)) ;

                estif(vi*3,     ui*3 + 1) -= four_visc2_ttimefac*(viscs2_(0,0,vi)*viscs2_(0,1,ui)+viscs2_(0,1,vi)*viscs2_(1,1,ui)) ;
                estif(vi*3 + 1, ui*3 + 1) -= four_visc2_ttimefac*(viscs2_(0,1,ui)*viscs2_(0,1,vi)+viscs2_(1,1,ui)*viscs2_(1,1,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<numnode; ++ui)
              {
                double v = two_visc_ttimefac*densfunct_(ui);
                for (int vi=0; vi<numnode; ++vi)
                {
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                   \
                    |  /          \   n+1               |
                +/- | | Du o nabla | u    , div eps (v) |
                    |  \          /   (i)               |
                    \                                   /
                  */
                  estif(vi*3,     ui*3    ) += v*(viscs2_(0,0,vi)*vderxy_(0,0)+viscs2_(0,1,vi)*vderxy_(1,0)) ;
                  estif(vi*3 + 1, ui*3    ) += v*(viscs2_(0,1,vi)*vderxy_(0,0)+viscs2_(1,1,vi)*vderxy_(1,0)) ;

                  estif(vi*3,     ui*3 + 1) += v*(viscs2_(0,0,vi)*vderxy_(0,1)+viscs2_(0,1,vi)*vderxy_(1,1)) ;
                  estif(vi*3 + 1, ui*3 + 1) += v*(viscs2_(0,1,vi)*vderxy_(0,1)+viscs2_(1,1,vi)*vderxy_(1,1)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<numnode; ++vi)
          {
            /* viscous stabilisation */
            eforce(vi*3    ) -= two_visc_timefac*(res_old_(0)*viscs2_(0, 0, vi)+res_old_(1)*viscs2_(0, 1, vi)) ;
            eforce(vi*3 + 1) -= two_visc_timefac*(res_old_(0)*viscs2_(0, 1, vi)+res_old_(1)*viscs2_(1, 1, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab == Fluid2::continuity_stab_yes)
      {
        const double timefac_tau_C=timefac*tau_C;
        const double timefac_timefac_tau_C=timefac*timefac*tau_C;
        const double timefac_timefac_tau_C_divunp=timefac_timefac_tau_C*mdiv_;

        for (int ui=0; ui<numnode; ++ui)
        {
          double v0 = timefac_timefac_tau_C*densderxy_(0, ui);
          double v1 = timefac_timefac_tau_C*densderxy_(1, ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            /* continuity stabilisation on left hand side */
            /*
              /                          \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
              \                          /
            */
            estif(vi*3,     ui*3    ) += v0*densderxy_(0, vi) ;
            estif(vi*3 + 1, ui*3    ) += v0*densderxy_(1, vi) ;

            estif(vi*3,     ui*3 + 1) += v1*densderxy_(0, vi) ;
            estif(vi*3 + 1, ui*3 + 1) += v1*densderxy_(1, vi) ;
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          /* continuity stabilisation on right hand side */
          eforce(vi*3    ) -= timefac_timefac_tau_C_divunp*densderxy_(0, vi) ;
          eforce(vi*3 + 1) -= timefac_timefac_tau_C_divunp*densderxy_(1, vi) ;

          if (loma)
          {
            /* continuity stabilisation of rhs term of continuity equation */
            eforce(vi*3    ) -= timefac_tau_C*densderxy_(0, vi)*rhscon_ ;
            eforce(vi*3 + 1) -= timefac_tau_C*densderxy_(1, vi)*rhscon_ ;
          }
        }
      }

      //----------------------------------------------------------------------
      //     STABILIZATION, CROSS-STRESS PART (RESIDUAL-BASED VMM)

      if (cross == Fluid2::cross_stress_stab_only_rhs || cross == Fluid2::cross_stress_stab)
      {
        if (cross == Fluid2::cross_stress_stab)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            double v = ttimetauM*conv_resM_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              /* cross-stress part on lhs */
              /*

                          /                        \
                         |  /            \          |
                      -  | | resM o nabla | Du , v  |
                         |  \            /          |
                          \                        /
              */
              double v2 = v*funct_(vi);
              estif(vi*3    , ui*3    ) -= v2 ;
              estif(vi*3 + 1, ui*3 + 1) -= v2 ;
            }
          }
        } // end cross-stress part on left hand side

        for (int vi=0; vi<numnode; ++vi)
        {
          /* cross-stress part on rhs */
          /*

                          /                         \
                         |  /            \           |
                         | | resM o nabla | u   , v  |
                         |  \            /  (i)      |
                          \                         /
          */
          double v = ttimetauM*funct_(vi);
          eforce(vi*3    ) += v*(res_old_(0)*mderxy_(0,0)+res_old_(1)*mderxy_(0,1));
          eforce(vi*3 + 1) += v*(res_old_(0)*mderxy_(1,0)+res_old_(1)*mderxy_(1,1));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if (reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      {
        const double ttimetauMtauM = ttimetauM*tau_M/fac;
        for (int vi=0; vi<numnode; ++vi)
        {
          /* Reynolds-stress part on rhs */
          /*

                  /                             \
                 |                               |
                 |  resM   , ( resM o nabla ) v  |
                 |                               |
                  \                             /
          */
          double v = ttimetauMtauM*conv_resM_(vi);
          eforce(vi*3    ) += v*res_old_(0);
          eforce(vi*3 + 1) += v*res_old_(1);
        }
      } // end Reynolds-stress part on right hand side

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv != Fluid2::fssgv_no)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(vi*3    ) -= vartfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                      +    derxy_(1, vi)*fsvderxy_(0, 1)
                                      +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
          eforce(vi*3 + 1) -= vartfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                      +    derxy_(0, vi)*fsvderxy_(1, 0)
                                      +2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
        }
      }
    }

    // linearization with respect to mesh motion
    if (emesh.shape()[0]==estif.shape()[0] and
        emesh.shape()[1]==estif.shape()[1])
    {

      // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
      // xGderiv_ == xjm_

      // mass + rhs
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = fac*funct_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*3    , ui*3    ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
          emesh(vi*3    , ui*3 + 1) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);

          emesh(vi*3 + 1, ui*3    ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
          emesh(vi*3 + 1, ui*3 + 1) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);

          emesh(vi*3 + 2, ui*3    ) += v*(velint_(2)-rhsmom_(2))*derxy_(0, ui);
          emesh(vi*3 + 2, ui*3 + 1) += v*(velint_(2)-rhsmom_(2))*derxy_(1, ui);
        }
      }

      vderiv_  = sum(evelnp(i,k) * deriv_(j,k), k);

//#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

//#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
//#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
//#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
//#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))

      for (int vi=0; vi<numnode; ++vi)
      {
        double v = timefacfac/det*funct_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {

          emesh(vi*3 + 0, ui*3 + 0) += v*(
          + convvelint_(1)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(vi*3 + 0, ui*3 + 1) += v*(
          + convvelint_(0)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(vi*3 + 1, ui*3 + 0) += v*(
          + convvelint_(1)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );

          emesh(vi*3 + 1, ui*3 + 1) += v*(
          + convvelint_(0)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );
        }
      }

      // pressure
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = press*timefacfac/det;
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*3 + 0, ui*3 + 1) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
          emesh(vi*3 + 1, ui*3 + 0) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
        }
      }

      // div u
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = timefacfac/det*functdens_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*3 + 2, ui*3 + 0) += v*(
          deriv_(0,ui)*vderiv_(1,1) - deriv_(1,ui)*vderiv_(1,0)
          ) ;

          emesh(vi*3 + 2, ui*3 + 1) += v*(
          deriv_(0,ui)*vderiv_(0,1) - deriv_(1,ui)*vderiv_(0,0)
          ) ;
        }
      }
    }
  } // loop gausspoints
}



//
// calculate stabilization parameter
//
void DRT::ELEMENTS::Fluid2Impl::Caltau(
  Fluid2*                                 ele,
  const blitz::Array<double,2>&           evelnp,
  const blitz::Array<double,2>&           fsevelnp,
  const blitz::Array<double,1>&           edensnp,
  const DRT::Element::DiscretizationType  distype,
  const enum Fluid2::TauType              whichtau,
  struct _MATERIAL*                       material,
  double&                           	  visc,
  const double                            timefac,
  const double                            dt,
  const enum Fluid2::TurbModelAction      turb_mod_action,
  double&                                 Cs,
  double&                                 visceff,
  const enum Fluid2::StabilisationAction  fssgv
  )
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule2D integrationrule_stabili=DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      integrationrule_stabili = DRT::UTILS::intrule_quad_1point;
      break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      integrationrule_stabili = DRT::UTILS::intrule_tri_1point;
      break;
    default:
      dserror("invalid discretization type for fluid2");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];

  DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
  DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
    case DRT::Element::tri3:
    case DRT::Element::quad4:
      mk = 0.333333333333333333333;
      break;
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::tri6:
      mk = 0.083333333333333333333;
      break;
    default:
      dserror("type unknown!\n");
  }

  // get velocities at element center
  velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

  // get density at element center
  const double dens = blitz::sum(funct_*edensnp);

  // get Jacobian matrix and determinant
  xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  const double det = xjm_(0,0)*xjm_(1,1) - xjm_(0,1)*xjm_(1,0);

  // check for degenerated elements
  if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

  // get characteristic element length: square root of element area
  double area=0;
  double a,b,c;

  switch (distype)
  {
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      a = (xyze_(0,0)-xyze_(0,1))*(xyze_(0,0)-xyze_(0,1))
          +(xyze_(1,0)-xyze_(1,1))*(xyze_(1,0)-xyze_(1,1)); /* line 0-1 squared */
      b = (xyze_(0,1)-xyze_(0,2))*(xyze_(0,1)-xyze_(0,2))
          +(xyze_(1,1)-xyze_(1,2))*(xyze_(1,1)-xyze_(1,2)); /* line 1-2 squared */
      c = (xyze_(0,2)-xyze_(0,0))*(xyze_(0,2)-xyze_(0,0))
          +(xyze_(1,2)-xyze_(1,0))*(xyze_(1,2)-xyze_(1,0)); /* diag 2-0 squared */
      area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      break;
    }
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    {
      a = (xyze_(0,0)-xyze_(0,1))*(xyze_(0,0)-xyze_(0,1))
          +(xyze_(1,0)-xyze_(1,1))*(xyze_(1,0)-xyze_(1,1)); /* line 0-1 squared */
      b = (xyze_(0,1)-xyze_(0,2))*(xyze_(0,1)-xyze_(0,2))
          +(xyze_(1,1)-xyze_(1,2))*(xyze_(1,1)-xyze_(1,2)); /* line 1-2 squared */
      c = (xyze_(0,2)-xyze_(0,0))*(xyze_(0,2)-xyze_(0,0))
          +(xyze_(1,2)-xyze_(1,0))*(xyze_(1,2)-xyze_(1,0)); /* diag 2-0 squared */
      area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      a = (xyze_(0,2)-xyze_(0,3))*(xyze_(0,2)-xyze_(0,3))
          +(xyze_(1,2)-xyze_(1,3))*(xyze_(1,2)-xyze_(1,3)); /* line 2-3 squared */
      b = (xyze_(0,3)-xyze_(0,0))*(xyze_(0,3)-xyze_(0,0))
          +(xyze_(1,3)-xyze_(1,0))*(xyze_(1,3)-xyze_(1,0)); /* line 3-0 squared */
      area += 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      break;
    }
    default: dserror("type unknown!\n");
  }

  const double hk = sqrt(area);

  //             compute global first derivates
  //
  // this is necessary only for the calculation of the
  // streamlength (required by the quasistatic formulation)
  //
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-          -+     +-    -+      +-    -+
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  dr    dr  |     |  dx  |      |  dr  |
          |            |  *  |      |   =  |      | for all k
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  ds    ds  |     |  dy  |      |  ds  |
          +-          -+     +-    -+      +-    -+

          Matrix is inverted analytically
  */
  // inverse of jacobian
  xji_(0,0) = ( xjm_(1,1))/det;
  xji_(0,1) = (-xjm_(0,1))/det;
  xji_(1,0) = (-xjm_(1,0))/det;
  xji_(1,1) = ( xjm_(0,0))/det;

  // compute global derivates
  derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

  // get velocity (np,i) derivatives at integration point
  vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

  // get velocity norm
  const double vel_norm = sqrt(blitz::sum(velint_*velint_));

  // normed velocity at element centre (currently not used)
  //if (vel_norm>=1e-6)
  //{
  //  velino_ = velint_/vel_norm;
  //}
  //else
  //{
  //  velino_ = 0.;
  //  velino_(0) = 1;
  //}

  // get streamlength (currently not used)
  //const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  //const double strle = 2.0/val;

  /*------------------------------------------------------------------*/
  /*                                                                  */
  /*                 GET EFFECTIVE VISCOSITY IN GAUSSPOINT            */
  /*                                                                  */
  /* This part is used to specify an effective viscosity. This eff.   */
  /* viscosity may be caused by a Smagorinsky model                   */
  /*                                                                  */
  /*          visc    = visc + visc                                   */
  /*              eff              turbulent                          */
  /*                                                                  */
  /* here, the latter turbulent viscosity is not a material thing,    */
  /* but a flow feature!                                              */
  /*                                                                  */
  /* Another cause for the necessity of an effective viscosity might  */
  /* be the use of a shear thinning Non-Newtonian fluid               */
  /*                                                                  */
  /*                            /         \                           */
  /*            visc    = visc | shearrate |                          */
  /*                eff         \         /                           */
  /*                                                                  */
  /*                                                                  */
  /* Mind that at the moment all stabilization (tau and viscous test  */
  /* functions if applied) are based on the material viscosity not    */
  /* the effective viscosity. We do this since we do not evaluate the */
  /* stabilisation parameter in the gausspoints but just once in the  */
  /* middle of the element.                                           */
  /*------------------------------------------------------------------*/

  // compute nonlinear viscosity according to the Carreau-Yasuda model
  if( material->mattyp != m_fluid ) CalVisc( material, visc);


  if (turb_mod_action == Fluid2::smagorinsky)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                                   +-                                 -+ 1
    //                               2   |          / h \           / h \    | -
    //    visc          = dens * lmix  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent           |      |          \   / ij        \   / ij |
    //                            |      +-                                 -+
    //                            |
    //                            |      |                                   |
    //                            |      +-----------------------------------+
    //                            |           'resolved' rate of strain
    //                         mixing length
    //

    double rateofstrain = 0;
    {
      blitz::Array<double,2> epsilon(2,2,blitz::ColumnMajorArray<2>());
      epsilon = 0.5 * ( vderxy_(i,j) + vderxy_(j,i) );

      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<2;mm++)
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
    //
    // mixing length set proportional to grid witdh
    //
    //                     lmix = Cs * hk

    double lmix = Cs * hk;

    //
    //          visc    = visc + visc
    //              eff              turbulent

    visceff = visc + dens * lmix * lmix * rateofstrain;
  }
  else
  {
    visceff = visc;
  }

  // calculate tau

  if (whichtau == Fluid2::franca_barrenechea_valentin_wall)
  {
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


    /* viscous : reactive forces */
    const double re1 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));

    /* convective : viscous forces */
    const double re2 = mk * dens * vel_norm * hk / (2.0 * visceff);

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
    tau_(1) = tau_(0);

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
    tau_(2) = vel_norm * hk * 0.5 * xi_tau_c / ( timefac * dens );
  }
  else if(whichtau == Fluid2::bazilevs)
  {
    /* INSTATIONARY FLOW PROBLEM, ONE-STEP-THETA, BDF2

    tau_M: Bazilevs et al.
                                                               1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+1       n+1          2         |
          tau  = | --- + u     * G u     + C * nu  * G : G |
             M   |   2           -          I        -   - |
                 | dt            -                   -   - |
                 +-                                       -+

   tau_C: Bazilevs et al., derived from the fine scale complement Shur
          operator of the pressure equation


                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */

    /*            +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
    */
    blitz::Array<double,2> G(2,2,blitz::ColumnMajorArray<2>());

    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
        G(nn,rr) = xji_(nn,0)*xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
        }
      }
    }

    /*            +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
    */
    double normG = 0;
    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
        normG+=G(nn,rr)*G(nn,rr);
      }
    }

    /*                      +----
           n+1       n+1     \     n+1          n+1
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
    */
    double Gnormu = 0;
    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
        Gnormu+=dens*velint_(nn)*G(nn,rr)*dens*velint_(rr);
      }
    }

    // definition of constant
    // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
    //  brought 144.0 from Austin...)
    const double CI = 12.0/mk;

    /*                                                         1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+1       n+1          2         |
          tau  = | --- + u     * G u     + C * nu  * G : G |
             M   |   2           -          I        -   - |
                 | dt            -                   -   - |
                 +-                                       -+
    */
    tau_(0) = 1.0/(timefac*sqrt((4.0*dens*dens)/(dt*dt)+Gnormu+CI*visceff*visceff*normG));
    tau_(1) = tau_(0);

    /*           +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
    */
    blitz::Array<double,1> g(2);

    for (int rr=0;rr<2;++rr)
    {
      g(rr) = xji_(rr,0);
      for (int mm=1;mm<2;++mm)
      {
        g(rr) += xji_(rr,mm);
      }
    }

    /*           +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
    */
    const double normgsq = g(0)*g(0)+g(1)*g(1);

    /*
                                1.0
                  tau  = -----------------
                     C            /     \
                          tau  * | g * g |
                             M    \-   -/
    */
    tau_(2) = 1./(tau_(0)*normgsq*timefac*timefac*dens*dens);

  }
  else if(whichtau == Fluid2::codina)
  {
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


    /* viscous : reactive forces */
    const double re1 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));

    /* convective : viscous forces */
    const double re2 = mk * dens * vel_norm * hk / (2.0 * visceff);

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff)+DSQR(0.5*dens*vel_norm*hk)) / ( timefac*dens*dens );

  }
  else
  {
    dserror("unknown definition of tau\n");
  }

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv == Fluid2::fssgv_artificial_all || fssgv == Fluid2::fssgv_artificial_small)
  {
    double fsvel_norm = 0.0;
    if (fssgv == Fluid2::fssgv_artificial_small)
    {
      // get fine-scale velocities at element center
      fsvelint_ = blitz::sum(funct_(j)*fsevelnp(i,j),j);

      // get fine-scale velocity norm
      fsvel_norm = sqrt(blitz::sum(fsvelint_*fsvelint_));
    }
    // get all-scale velocity norm
    else fsvel_norm = vel_norm;

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = mk * dens * fsvel_norm * hk / visc; /* convective : viscous forces */
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(dens)*DSQR(fsvel_norm))/(2.0*visc*xi);

  }
  else if (fssgv == Fluid2::fssgv_Smagorinsky_all or
           fssgv == Fluid2::fssgv_Smagorinsky_small)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                                      +-                                 -+ 1
    //                                  2   |          / h \           / h \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    double rateofstrain = 0.0;
    {
      // get fine-scale or all-scale velocity (np,i) derivatives at element center
      if (fssgv == Fluid2::fssgv_Smagorinsky_small)
           fsvderxy_ = blitz::sum(derxy_(j,k)*fsevelnp(i,k),k);
      else fsvderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

      blitz::Array<double,2> epsilon(2,2,blitz::ColumnMajorArray<2>());
      epsilon = 0.5 * ( fsvderxy_(i,j) + fsvderxy_(j,i) );

      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<2;mm++)
        {
          rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
        }
      }
      rateofstrain *= 2.0;
      rateofstrain = sqrt(rateofstrain);
    }
    //
    // Choices of the fine-scale Smagorinsky constant Cs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)

    vart_ = dens * Cs * Cs * hk * hk * rateofstrain;
  }
}



//
// calculate material viscosity    u.may 05/08
//
void DRT::ELEMENTS::Fluid2Impl::CalVisc(
  const struct _MATERIAL*                 material,
  double&                                 visc)
{

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index

  // compute shear rate
  double rateofshear = 0.0;
  blitz::Array<double,2> epsilon(2,2,blitz::ColumnMajorArray<2>());   // strain rate tensor
  epsilon = 0.5 * ( vderxy_(i,j) + vderxy_(j,i) );

  for(int rr=0;rr<2;rr++)
    for(int mm=0;mm<2;mm++)
      rateofshear += epsilon(rr,mm)*epsilon(rr,mm);

  rateofshear = sqrt(2.0*rateofshear);

  if(material->mattyp == m_carreauyasuda)
  {
    double nu_0   = material->m.carreauyasuda->nu_0;    // parameter for zero-shear viscosity
    double nu_inf = material->m.carreauyasuda->nu_inf;  // parameter for infinite-shear viscosity
    double lambda = material->m.carreauyasuda->lambda;  // parameter for characteristic time
    double a 	  = material->m.carreauyasuda->a_param; // constant parameter
    double b      = material->m.carreauyasuda->b_param; // constant parameter

    // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    const double tmp = pow(lambda*rateofshear,b);
    visc = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
  }
  else if(material->mattyp == m_modpowerlaw)
  {
    // get material parameters
    double m     = material->m.modpowerlaw->m_cons;     // consistency constant
    double delta = material->m.modpowerlaw->delta;      // safety factor
    double a     = material->m.modpowerlaw->a_exp;      // exponent

    // compute viscosity according to a modified power law model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    visc = m * pow((delta + rateofshear), (-1)*a);
  }
  else
    dserror("material type is not yet implemented");
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Impl::BodyForce(Fluid2*      ele,
                                          const double time)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique surface Neumann condition
  DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);

  if (myneumcond.size()>1) 
    dserror("more than one SurfaceNeumann cond on one node");

  if (myneumcond.size()==1)
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
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
        // A negative time value indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );
    const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      for(int isd=0;isd<2;isd++)
      {
        // get factor given by spatial function
        if (functions)
          functnum = (*functions)[isd];
        else 
          functnum = -1;

        if (functnum>0)
        {
          // evaluate function at the position of the current node
          functionfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X());
        }
        else
          functionfac = 1.0;

        // compute and store the (normalized) bodyforce value
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac*functionfac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_ = 0.;
  }
}



#endif
#endif
