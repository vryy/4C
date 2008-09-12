/*----------------------------------------------------------------------*/
/*!
\file fluid2_stationary.cpp

\brief Internal implementation of Fluid2 element (stationary formulation)

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

#include "fluid2_stationary.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_LAPACK.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Stationary* DRT::ELEMENTS::Fluid2Stationary::StationaryImpl(DRT::ELEMENTS::Fluid2* f2)
{
  switch (f2->NumNode())
  {
  case 4:
  {
    static Fluid2Stationary* f4;
    if (f4==NULL)
      f4 = new Fluid2Stationary(4);
    return f4;
  }
  case 8:
  {
    static Fluid2Stationary* f8;
    if (f8==NULL)
      f8 = new Fluid2Stationary(8);
    return f8;
  }
  case 9:
  {
    static Fluid2Stationary* f9;
    if (f9==NULL)
      f9 = new Fluid2Stationary(9);
    return f9;
  }
  case 3:
  {
    static Fluid2Stationary* f3;
    if (f3==NULL)
      f3 = new Fluid2Stationary(3);
    return f3;
  }
  case 6:
  {
    static Fluid2Stationary* f6;
    if (f6==NULL)
      f6 = new Fluid2Stationary(6);
    return f6;
  }

  default:
    dserror("node number %d not supported", f2->NumNode());
  }
  return NULL;
}


DRT::ELEMENTS::Fluid2Stationary::Fluid2Stationary(int iel)
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
    //velino_(2),
    velint_(2),
    fsvelint_(2),
    gradp_(2),
    tau_(3),
    viscs2_(2,2,iel_,blitz::ColumnMajorArray<3>()),
    conv_c_(iel_),
    mdiv_(),
    vdiv_(),
    rhsmom_(2),
    conv_old_(2),
    visc_old_(2),
    res_old_(2),
    conv_resM_(iel_),
    xder2_(3,2,blitz::ColumnMajorArray<2>())
{
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                  vg 08/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Stationary::Sysmat(
  Fluid2*                                 ele,
  const blitz::Array<double,2>&           evelnp,
  const blitz::Array<double,2>&           fsevelnp,
  const blitz::Array<double,1>&           eprenp,
  const blitz::Array<double,1>&           edensnp,
  blitz::Array<double,2>&                 estif,
  blitz::Array<double,1>&                 eforce,
  struct _MATERIAL*                       material,
  double                                  pseudotime,
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
  const double                            Cs
  )
{

  // set element data
  const DRT::Element::DiscretizationType distype = ele->Shape();
  const int numnode = iel_;

  // get node coordinates and number of elements per node
  DRT::Node** const nodes = ele->Nodes();
  for (int inode=0; inode<iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // dead load in element nodes
  BodyForce(ele,pseudotime);

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

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  CalTauStationary(ele,
                   evelnp,
                   fsevelnp,
                   edensnp,
                   distype,
                   visc,
                   fssgv,
                   Cs);

  // in case of viscous stabilisation decide whether to use GLS or usfemM
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
    // coordinates of the current integration point
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

    // get velocity (np,i) derivatives at integration point
    vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

    // get momentum (np,i) derivatives at integration point
    mderxy_ = blitz::sum(densderxy_(j,k)*evelnp(i,k),k);

    // get fine-scale velocity (np,i) derivatives at integration point
    if (fssgv != Fluid2::fssgv_no) fsvderxy_ = blitz::sum(derxy_(j,k)*fsevelnp(i,k),k);
    else fsvderxy_ = 0.;

    // get pressure gradients
    gradp_ = blitz::sum(derxy_(i,j)*eprenp(j),j);

    // get pressure at integration point
    double press = blitz::sum(funct_*eprenp);

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_ = blitz::sum(edeadng_(i,j)*densfunct_(j),j);

    // perform integration for entire matrix and rhs

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;

    /*------------------------- evaluate rhs vector at integration point ---*/
    // histvectors are always zero in stationary case (!):
    rhsmom_ = bodyforce_(i);

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_ = blitz::sum(vderxy_(i, j)*velint_(j), j);

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
    conv_c_ = blitz::sum(derxy_(j,i)*velint_(j), j);

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

    /*--------------------------------- now build single stiffness terms ---*/
    // evaluate residual once for all stabilisation right hand sides
    res_old_ = -rhsmom_+(conv_old_+gradp_-2*visc*visc_old_);

    /*
      This is the operator

                /               \
               | resM    o nabla |
                \    (i)        /

      required for (lhs) cross- and (rhs) Reynolds-stress calculation

    */

    if(cross    == Fluid2::cross_stress_stab ||
       reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      conv_resM_ =  blitz::sum(res_old_(j)*densderxy_(j,i),j);

   {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      for (int ui=0; ui<numnode; ++ui)
      {
        double v = fac*conv_c_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          /* convection, convective part */
          /*

                       /                                \
                      |  /                   \          |
                      | | dens   u   o nabla | Du , v   |
                      |  \   (i)  (i)        /          |
                       \                                /

          */
          double v2 = v*funct_(vi) ;
          estif(vi*3,     ui*3    ) += v2 ;
          estif(vi*3 + 1, ui*3 + 1) += v2 ;
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
          estif(vi*3,     ui*3    ) += visc*fac*(2.0*derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)) ;
          estif(vi*3,     ui*3 + 1) += visc*fac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*3 + 1, ui*3    ) += visc*fac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*3 + 1, ui*3 + 1) += visc*fac*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 2.0*derxy_(1, ui)*derxy_(1, vi)) ;
        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        double v = -fac*funct_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          /*  prssure term */
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
        double v = fac*functdens_(vi);
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
          double v = fac*funct_(vi);
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

      /* right hand side */
      for (int vi=0; vi<numnode; ++vi)
      {
        /* convection */
        double v = -fac*funct_(vi);
        eforce(vi*3    ) += v*(velint_(0)*vderxy_(0, 0)+velint_(1)*vderxy_(0, 1)) ;
        eforce(vi*3 + 1) += v*(velint_(0)*vderxy_(1, 0)+velint_(1)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        /* pressure */
        double v = press*fac;
        eforce(vi*3    ) += v*derxy_(0, vi) ;
        eforce(vi*3 + 1) += v*derxy_(1, vi) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        /* viscosity */
        eforce(vi*3    ) -= visc*fac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                                      +
                                      derxy_(1, vi)*vderxy_(0, 1)
                                      +
                                      derxy_(1, vi)*vderxy_(1, 0)) ;
        eforce(vi*3 + 1) -= visc*fac*(derxy_(0, vi)*vderxy_(0, 1)
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
        eforce(vi*3 + 2) -= fac*functdens_(vi)*mdiv_ ;
      }

      if (loma)
      {
        double v = -(2.0/3.0)*visc*fac ;
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
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if(pspg == Fluid2::pstab_use_pspg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          double v = tau_Mp*conv_c_(ui) ;
          for (int vi=0; vi<numnode; ++vi)
          {

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
          double v = -2.0*visc*tau_Mp;
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
            estif(vi*3 + 2, ui*3 + 2) += tau_Mp*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            double v = tau_Mp*densfunct_(ui);
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
          eforce(vi*3 + 2) -= tau_Mp*(res_old_(0)*derxy_(0, vi)
                                      +
                                      res_old_(1)*derxy_(1, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg == Fluid2::convective_stab_supg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          double v = tau_M*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
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
          double v = tau_M*conv_c_(vi);
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
            double v = -2.0*visc*tau_M*conv_c_(vi);
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

        if (newton)
        {
          {
            double v0 = velint_(0)*vderxy_(0, 0) + velint_(1)*vderxy_(0, 1);
            double v1 = velint_(0)*vderxy_(1, 0) + velint_(1)*vderxy_(1, 1);

            for (int ui=0; ui<numnode; ++ui)
            {
              double v = tau_M*densfunct_(ui);
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

          for (int ui=0; ui<numnode; ++ui)
          {
            double v = tau_M*densfunct_(ui);
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
              double v = -2.0*visc*tau_M*densfunct_(ui);
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
            double v = -tau_M*densfunct_(ui);
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

        /* right hand side */
        for (int vi=0; vi<numnode; ++vi)
        {
          // supg stabilisation
          double v = -tau_M*conv_c_(vi);
          eforce(vi*3    ) += v*res_old_(0) ;
          eforce(vi*3 + 1) += v*res_old_(1) ;
        }
      }


      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART
      if (higher_order_ele)
      {
        if(vstab != Fluid2::viscous_stab_none)
        {
          const double two_visc_tauMp = vstabfac*2.0*visc*tau_Mp;

          // viscous stabilization either on left hand side or on right hand side
          if (vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_usfem)
          {
            const double four_visc2_tauMp = vstabfac*4.0*visc*visc*tau_Mp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<numnode; ++ui)
            {
              double v = two_visc_tauMp*conv_c_(ui)
                         ;
              for (int vi=0; vi<numnode; ++vi)
              {
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
                estif(vi*3,     ui*3 + 2) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(0, 0, vi)
                                                              +
                                                              derxy_(1, ui)*viscs2_(0, 1, vi)) ;
                estif(vi*3 + 1, ui*3 + 2) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(0, 1, vi)
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
                estif(vi*3,     ui*3    ) -= four_visc2_tauMp*(viscs2_(0,0,ui)*viscs2_(0,0,vi)+viscs2_(0,1,ui)*viscs2_(0,1,vi)) ;
                estif(vi*3 + 1, ui*3    ) -= four_visc2_tauMp*(viscs2_(0,0,ui)*viscs2_(0,1,vi)+viscs2_(0,1,ui)*viscs2_(1,1,vi)) ;

                estif(vi*3,     ui*3 + 1) -= four_visc2_tauMp*(viscs2_(0,0,vi)*viscs2_(0,1,ui)+viscs2_(0,1,vi)*viscs2_(1,1,ui)) ;
                estif(vi*3 + 1, ui*3 + 1) -= four_visc2_tauMp*(viscs2_(0,1,ui)*viscs2_(0,1,vi)+viscs2_(1,1,ui)*viscs2_(1,1,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<numnode; ++ui)
              {
                double v = two_visc_tauMp*densfunct_(ui);
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
            eforce(vi*3    ) -= two_visc_tauMp*(res_old_(0)*viscs2_(0, 0, vi)+res_old_(1)*viscs2_(0, 1, vi)) ;
            eforce(vi*3 + 1) -= two_visc_tauMp*(res_old_(0)*viscs2_(0, 1, vi)+res_old_(1)*viscs2_(1, 1, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab == Fluid2::continuity_stab_yes)
      {
        const double tau_C_divunp=tau_C*mdiv_;

        for (int ui=0; ui<numnode; ++ui)
        {
          double v0 = tau_C*densderxy_(0, ui);
          double v1 = tau_C*densderxy_(1, ui);
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
          eforce(vi*3    ) -= tau_C_divunp*densderxy_(0, vi) ;
          eforce(vi*3 + 1) -= tau_C_divunp*densderxy_(1, vi) ;
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
            double v = tau_M*conv_resM_(ui);
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
          double v = tau_M*funct_(vi);
          eforce(vi*3    ) += v*(res_old_(0)*mderxy_(0,0)+res_old_(1)*mderxy_(0,1));
          eforce(vi*3 + 1) += v*(res_old_(0)*mderxy_(1,0)+res_old_(1)*mderxy_(1,1));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if (reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      {
        const double tauMtauM = tau_M*tau_M/fac;
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
          double v = tauMtauM*conv_resM_(vi);
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
  }
  return;
}



//
// calculate stabilization parameter
//
void DRT::ELEMENTS::Fluid2Stationary::CalTauStationary(
  Fluid2*                                 ele,
  const blitz::Array<double,2>&           evelnp,
  const blitz::Array<double,2>&           fsevelnp,
  const blitz::Array<double,1>&           edensnp,
  const DRT::Element::DiscretizationType  distype,
  const double                            visc,
  const enum Fluid2::StabilisationAction  fssgv,
  const double                            Cs
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

  // calculate stabilization parameters for stationary case

  // compute tau_Mu and tau_Mp
  /* convective : viscous forces */
  const double re = mk * dens * vel_norm * hk / (2.0 * visc);
  const double xi = DMAX(re, 1.0);
  tau_(0) = (DSQR(hk)*mk)/(4.0*visc*xi);
  tau_(1) = tau_(0);

  // compute tau_C
  const double xi_tau_c = DMIN(re, 1.0);
  tau_(2) = 0.5*vel_norm*hk*xi_tau_c/dens;

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
    const double re_sv = mk * dens * fsvel_norm * hk / visc; /* convective : viscous forces */
    const double xi_sv = DMAX(re_sv,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(dens)*DSQR(fsvel_norm))/(2.0*visc*xi_sv);

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
      // get fine-scale or all-scale velocity derivatives at element center
      if (fssgv == Fluid2::fssgv_Smagorinsky_small)
            fsvderxy_ = blitz::sum(derxy_(j,k)*fsevelnp(i,k),k);
      else  fsvderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

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



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Stationary::BodyForce(Fluid2*      ele,
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

        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*functionfac;
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
