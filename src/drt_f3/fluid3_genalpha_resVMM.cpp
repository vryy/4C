/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.cpp

\brief Internal implementation of Fluid3 element with a generalised alpha
       time integration.

       This element is designed for the solution of the Navier-Stokes
       equations using a residual based stabilised method. The
       stabilisation terms are derived in a variational multiscale sense.

       Subscales are either treated as quasi-static or time dependent.

       There is no ALE-ability of the element up to now.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_genalpha_resVMM.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_LAPACK.h>


DRT::ELEMENTS::Fluid3GenalphaResVMM* DRT::ELEMENTS::Fluid3GenalphaResVMM::Impl(DRT::ELEMENTS::Fluid3* ele)
{
  switch (ele->NumNode())
  {
  case 8:
  {
    static Fluid3GenalphaResVMM* f8;
    if (f8==NULL)
      f8 = new Fluid3GenalphaResVMM(8);
    return f8;
  }
  case 20:
  {
    static Fluid3GenalphaResVMM* f20;
    if (f20==NULL)
      f20 = new Fluid3GenalphaResVMM(20);
    return f20;
  }
  case 27:
  {
    static Fluid3GenalphaResVMM* f27;
    if (f27==NULL)
      f27 = new Fluid3GenalphaResVMM(27);
    return f27;
  }
  case 4:
  {
    static Fluid3GenalphaResVMM* f4;
    if (f4==NULL)
      f4 = new Fluid3GenalphaResVMM(4);
    return f4;
  }
  case 10:
  {
    static Fluid3GenalphaResVMM* f10;
    if (f10==NULL)
      f10 = new Fluid3GenalphaResVMM(10);
    return f10;
  }
  case 6:
  {
    static Fluid3GenalphaResVMM* f6;
    if (f6==NULL)
      f6 = new Fluid3GenalphaResVMM(6);
    return f6;
  }
  case 15:
  {
    static Fluid3GenalphaResVMM* f15;
    if (f15==NULL)
      f15 = new Fluid3GenalphaResVMM(15);
    return f15;
  }
  default:
    dserror("node number %d not supported", ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
  |  constructor allocating arrays whose sizes may depend on the number |
  | of nodes of the element                                             |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3GenalphaResVMM::Fluid3GenalphaResVMM(int iel)
  : iel_        (iel),
// fine-scale subgrid viscosity
    vart_(),
// nodal data
//-----------------------+------------+------------------------------------
//                  dim  | derivative | node
    xyze_         (  3   ,              iel_,blitz::ColumnMajorArray<2>()),
    edeadaf_      (  3   ,              iel_,blitz::ColumnMajorArray<2>()),
//-----------------------+------------+------------------------------------
// gausspoint data
//------------------------------------------------------------------------
//                  dim  | derivative | node
//-----------------------+------------+------------------------------------
    funct_        (                     iel_                             ),
    deriv_        (            3      , iel_,blitz::ColumnMajorArray<2>()),
    deriv2_       (            6      , iel_,blitz::ColumnMajorArray<2>()),
    derxy_        (            3      , iel_,blitz::ColumnMajorArray<2>()),
    derxy2_       (            6      , iel_,blitz::ColumnMajorArray<2>()),
    viscs2_       (  3   ,     3      , iel_,blitz::ColumnMajorArray<3>()),
    xjm_          (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    xji_          (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    xder2_        (  6   ,     3            ,blitz::ColumnMajorArray<2>()),
    accintam_     (  3                                                   ),
    velintnp_     (  3                                                   ),
    velintaf_     (  3                                                   ),
    csvelintaf_   (  3                                                   ),
    fsvelintaf_   (  3                                                   ),
    csconvintaf_  (  3                                                   ),
    aleconvintaf_ (  3                                                   ),
    pderxynp_     (  3                                                   ),
    vderxynp_     (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    vderxyaf_     (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    csvderxyaf_   (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    fsvderxyaf_   (  3   ,     3            ,blitz::ColumnMajorArray<2>()),
    vderxy2af_    (  3   ,     6            ,blitz::ColumnMajorArray<2>()),
    bodyforceaf_  (  3                                                   ),
    conv_c_af_    (                     iel_                             ),
    conv_r_af_    (  3   ,     3      , iel_,blitz::ColumnMajorArray<3>()),
    conv_g_af_    (                     iel_                             ),
//----------------------+------------+------------------------------------
// element data
//------------------------------------------------------------------------
    tau_          (3),
    svelaf_       (3),
    convaf_old_   (3),
    convaf_s_     (3),
    convsubaf_old_(3),
    viscaf_old_   (3),
    resM_         (3),
    conv_resM_    (                      iel_),
    conv_subaf_   (                      iel_)
{
}


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3GenalphaResVMM::Sysmat(
  Fluid3*                                               ele,
  Epetra_SerialDenseMatrix&                             elemat,
  Epetra_SerialDenseVector&                             elevec,
  const blitz::Array<double,2>&                         edispnp,
  const blitz::Array<double,2>&                         egridvaf,
  const blitz::Array<double,2>&                         evelnp,
  const blitz::Array<double,1>&                         eprenp,
  const blitz::Array<double,2>&                         eaccam,
  const blitz::Array<double,2>&                         evelaf,
  const blitz::Array<double,2>&                         csevelaf,
  const blitz::Array<double,2>&                         fsevelaf,
  const blitz::Array<double,2>&                         cseconvaf,
  const struct _MATERIAL*                               material,
  const double                                          alphaM,
  const double                                          alphaF,
  const double                                          gamma,
  const double                                          dt,
  const double                                          time,
  const bool                                            newton,
  const bool                                            higher_order_ele,
  const enum Fluid3::StabilisationAction                fssgv,
  const enum Fluid3::StabilisationAction                tds,
  const enum Fluid3::StabilisationAction                inertia,
  const enum Fluid3::StabilisationAction                pspg,
  const enum Fluid3::StabilisationAction                supg,
  const enum Fluid3::StabilisationAction                vstab,
  const enum Fluid3::StabilisationAction                cstab,
  const enum Fluid3::StabilisationAction                cross,
  const enum Fluid3::StabilisationAction                reynolds,
  const enum Fluid3::TauType                            whichtau,
  const enum Fluid3::TurbModelAction                    turb_mod_action,
  double&                                               Cs,
  double&                                               Cs_delta_sq,
  double&                                               visceff,
  const double                                          l_tau,
#ifdef PERF
  RefCountPtr<Time>                                     timeelederxy2    ,
  RefCountPtr<Time>                                     timeelederxy     ,
  RefCountPtr<Time>                                     timeeletau       ,
  RefCountPtr<Time>                                     timeelegalerkin  ,
  RefCountPtr<Time>                                     timeelepspg      ,
  RefCountPtr<Time>                                     timeelesupg      ,
  RefCountPtr<Time>                                     timeelecstab     ,
  RefCountPtr<Time>                                     timeelevstab     ,
  RefCountPtr<Time>                                     timeelecrossrey  ,
  RefCountPtr<Time>                                     timeeleintertogp ,
  RefCountPtr<Time>                                     timeeleseteledata,
  RefCountPtr<Time>                                     timeeletdextras  ,
#endif
  const bool                                            compute_elemat
  )
{
#ifdef PERF
    RefCountPtr<TimeMonitor> timeeleseteledata_ref = rcp(new TimeMonitor(*timeeleseteledata));
#endif
  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  blitz::Range       _ = blitz::Range::all();

  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt

  const double timealphaF = time-(1-alphaF)*dt;

  //------------------------------------------------------------------
  //                      SET MATERIAL DATA
  //-------------------------------------------------
  // check here, if we really have a fluid !!
    if( material->mattyp != m_carreauyasuda
    	&&	material->mattyp != m_modpowerlaw
        && material->mattyp != m_fluid)
    	  dserror("Material law is not a fluid");

    // get viscosity
    double visc = 0.0;
    double viscturb = 0.0;
    if(material->mattyp == m_fluid)
  	  visc = material->m.fluid->viscosity;

  //------------------------------------------------------------------
  //                      SET ELEMENT DATA
  //------------------------------------------------------------------
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
    for (int inode=0; inode<iel_; inode++)
    {
      xyze_(0,inode) += edispnp(0,inode);
      xyze_(1,inode) += edispnp(1,inode);
      xyze_(2,inode) += edispnp(2,inode);
    }
  }

  // dead load in element nodes
  GetNodalBodyForce(ele,timealphaF);

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }
#ifdef PERF
  timeeleseteledata_ref = null;
#endif

  //----------------------------------------------------------------------------
  //            STABILIZATION PARAMETER, SMAGORINSKY MODEL
  //      and everything else that is evaluated in the element center
  //
  // This has to be done before anything else is calculated because we use
  // the same arrays internally.
  //----------------------------------------------------------------------------
#ifdef PERF
    RefCountPtr<TimeMonitor> timeeletau_ref = rcp(new TimeMonitor(*timeeletau));
#endif

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili=DRT::UTILS::intrule3D_undefined;
  switch (distype)
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
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  const double wquad = intpoints_onepoint.qwgt[0];

  DRT::UTILS::shape_function_3D       (funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // get transposed Jacobian matrix and determinant
  //
  //        +-            -+ T      +-            -+
  //        | dx   dx   dx |        | dx   dy   dz |
  //        | --   --   -- |        | --   --   -- |
  //        | dr   ds   dt |        | dr   dr   dr |
  //        |              |        |              |
  //        | dy   dy   dy |        | dx   dy   dz |
  //        | --   --   -- |   =    | --   --   -- |
  //        | dr   ds   dt |        | ds   ds   ds |
  //        |              |        |              |
  //        | dz   dz   dz |        | dx   dy   dz |
  //        | --   --   -- |        | --   --   -- |
  //        | dr   ds   dt |        | dt   dt   dt |
  //        +-            -+        +-            -+
  //
  // The Jacobian is computed using the formula
  //
  //            +-----
  //   dx_j(r)   \      dN_k(r)
  //   -------  = +     ------- * (x_j)_k
  //    dr_i     /       dr_i       |
  //            +-----    |         |
  //            node k    |         |
  //                  derivative    |
  //                   of shape     |
  //                   function     |
  //                           component of
  //                          node coordinate
  //
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        xjm_(nn,rr)+=deriv_(nn,mm)*xyze_(rr,mm);
      }
    }
  }

  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  vol_ = wquad*det;

  // get element length for tau_M and tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol_/PI),(1.0/3.0))/sqrt(3.0);

  //
  //             compute global first derivates
  //
  // this is necessary only for the calculation of the
  // streamlength (required by the quasistatic formulation) and
  // the Smagorinsky model.
  //
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

  */

  // inverse of jacobian (transposed)
  /*
          +-                 -+     +-                 -+ -1
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dx    dx    dx   |     |  dr    dr    dr   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |  =  |  --    --    --   |
          |  dy    dy    dy   |     |  ds    ds    ds   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dz    dz    dz   |     |  dt    dt    dt   |
          +-                 -+     +-                 -+

  */
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
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<iel_;++rr)
    {
      derxy_(nn,rr)=deriv_(0,rr)*xji_(nn,0);
      for (int mm=1;mm<3;++mm)
      {
        derxy_(nn,rr)+=deriv_(mm,rr)*xji_(nn,mm);
      }
    }
  }

  // get velocities (n+alpha_F,i) at integration point
  //
  //                 +-----
  //       n+af       \                  n+af
  //    vel    (x) =   +      N (x) * vel
  //                  /        j         j
  //                 +-----
  //                 node j
  //
  for (int rr=0;rr<3;++rr)
  {
    velintaf_(rr)=funct_(0)*evelaf(rr,0);
    for (int mm=1;mm<iel_;++mm)
    {
      velintaf_(rr)+=funct_(mm)*evelaf(rr,mm);
    }
  }

  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //       n+af      +-----  dN (x)
  //   dvel    (x)    \        k         n+af
  //   ----------- =   +     ------ * vel
  //       dx         /        dx        k
  //         j       +-----      j
  //                 node k
  //
  // j : direction of derivative x/y/z
  //
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      vderxyaf_(nn,rr)=derxy_(rr,0)*evelaf(nn,0);
      for (int mm=1;mm<iel_;++mm)
      {
        vderxyaf_(nn,rr)+=derxy_(rr,mm)*evelaf(nn,mm);
      }
    }
  }


  // get velocities (n+1,i)  at integration point
  //
  //                +-----
  //       n+1       \                  n+1
  //    vel   (x) =   +      N (x) * vel
  //                 /        j         j
  //                +-----
  //                node j
  //
  for (int rr=0;rr<3;++rr)
  {
    velintnp_(rr)=funct_(0)*evelnp(rr,0);
    for (int mm=1;mm<iel_;++mm)
    {
      velintnp_(rr)+=funct_(mm)*evelnp(rr,mm);
    }
  }

  // get velocity norms
  const double vel_normaf = sqrt(velintaf_(0)*velintaf_(0)
                                 +
                                 velintaf_(1)*velintaf_(1)
                                 +
                                 velintaf_(2)*velintaf_(2));
  const double vel_normnp = sqrt(velintnp_(0)*velintnp_(0)
                                 +
                                 velintnp_(1)*velintnp_(1)
                                 +
                                 velintnp_(2)*velintnp_(2));

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
  if( material->mattyp != m_fluid )
	  CalVisc_AtGaussianPoint( material, visc);


  if (turb_mod_action == Fluid3::smagorinsky_with_wall_damping
      ||
      turb_mod_action == Fluid3::smagorinsky)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                            +-                                 -+ 1
    //                        2   |          / h \           / h \    | -
    //    visc          = lmix  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent    |      |          \   / ij        \   / ij |
    //                     |      +-                                 -+
    //                     |
    //                     |      |                                   |
    //                     |      +-----------------------------------+
    //                     |           'resolved' rate of strain
    //                    mixing length
    //

    double rateofstrain = 0;
    {
      blitz::Array<double,2> epsilon(3,3);

      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<3;mm++)
        {
          epsilon(rr,mm) = 0.5 * ( vderxyaf_(rr,mm) + vderxyaf_(mm,rr) );
        }
      }

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
    //
    //             Cs dynamic  (Germano model. Use several filter
    //                          resolutions to determine Cs)

    if (turb_mod_action == Fluid3::smagorinsky_with_wall_damping)
    {
      // since the Smagorinsky constant is only valid if hk is in the inertial
      // subrange of turbulent flows, the mixing length is damped in the
      // viscous near wall region using the van Driest damping function
      /*
                                       /         /   y+ \ \
                     lmix = Cs * hk * | 1 - exp | - ---- | |
                                       \         \   A+ / /
      */
      // A+ is a constant parameter, y+ the distance from the wall in wall
      // units
      const double A_plus = 26.0;
      double y_plus;

      // the integration point coordinate is defined by the isometric approach
      /*
                  +-----
                   \
              x =   +      N (x) * x
                   /        j       j
                  +-----
                  node j
      */
      blitz::Array<double,1> centernodecoord(3);

      for(int rr=0;rr<3;++rr)
      {
        centernodecoord(rr)=funct_(0)*xyze_(rr,0);

        for(int mm=1;mm<iel_;++mm)
        {
          centernodecoord(rr)+=funct_(mm)*xyze_(rr,mm);
        }
      }

      if(centernodecoord(1)>0)
      {
        y_plus=(1.0-centernodecoord(1))/l_tau;
      }
      else
      {
        y_plus=(1.0+centernodecoord(1))/l_tau;
      }

      // multiply with van Driest damping function
      Cs *= (1.0-exp(-y_plus/A_plus));
    }

    const double hk = pow((vol_),(1.0/3.0));

    //
    // mixing length set proportional to grid witdh
    //
    //                     lmix = Cs * hk

    double lmix = Cs * hk;

    Cs_delta_sq = lmix * lmix;

    //
    //          visc    = visc + visc
    //              eff              turbulent

    viscturb = Cs_delta_sq * rateofstrain;
    visceff = visc + viscturb;
  }
  else if(turb_mod_action == Fluid3::dynamic_smagorinsky)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                            +-                                 -+ 1
    //                        2   |          / h \           / h \    | -
    //    visc          = lmix  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent    |      |          \   / ij        \   / ij |
    //                     |      +-                                 -+
    //                     |
    //                     |      |                                   |
    //                     |      +-----------------------------------+
    //                     |           'resolved' rate of strain
    //                    mixing length
    //               provided by the dynamic model
    //            procedure and stored in Cs_delta_sq
    //

    double rateofstrain = 0;
    {
      blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());

      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<3;mm++)
        {
          epsilon(rr,mm) = 0.5 * ( vderxyaf_(rr,mm) + vderxyaf_(mm,rr) );
        }
      }

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

    viscturb = Cs_delta_sq * rateofstrain;
    visceff = visc + viscturb;

    // for evaluation of statistics: remember the 'real' Cs
    Cs=sqrt(Cs_delta_sq)/pow((vol_),(1.0/3.0));
  }
  else
  {
    visceff = visc;
  }

  if(tds == Fluid3::subscales_time_dependent)
  {
    //-------------------------------------------------------
    //          TAUS FOR TIME DEPENDENT SUBSCALES
    //-------------------------------------------------------

    if(whichtau == Fluid3::bazilevs)
    {
      /* INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA

         tau_M: Bazilevs et al. + ideas from Codina
                                                         1.0
                 +-                                 -+ - ---
                 |                                   |   2.0
             td  |  n+af      n+af         2         |
          tau  = | u     * G u     + C * nu  * G : G |
             M   |         -          I        -   - |
                 |         -                   -   - |
                 +-                                 -+

         tau_C: Bazilevs et al., derived from the fine scale complement Shur
                                 operator of the pressure equation


                       td         1.0
                    tau  = -----------------
                       C       td   /     \
                            tau  * | g * g |
                               M    \-   -/
      */

      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      blitz::Array<double,2> G(3,3,blitz::ColumnMajorArray<2>());

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
          }
        }
      }

      /*          +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
      */
      double normG = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          normG+=G(nn,rr)*G(nn,rr);
        }
      }

      /*                    +----
           n+af      n+af    \     n+af         n+af
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
      */
      double Gnormu = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          Gnormu+=velintaf_(nn)*G(nn,rr)*velintaf_(rr);
        }
      }

      // definition of constant
      // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
      //  brought 144.0 from Austin...)
      const double CI = 12.0/mk;

      /*                                                 1.0
                 +-                                 -+ - ---
                 |                                   |   2.0
                 |  n+af      n+af         2         |
          tau  = | u     * G u     + C * nu  * G : G |
             M   |         -          I        -   - |
                 |         -                   -   - |
                 +-                                 -+
      */
      tau_(0) = 1.0/sqrt(Gnormu+CI*visceff*visceff*normG);
      tau_(1) = tau_(0);

      /*         +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
      */
      blitz::Array<double,1> g(3);

      for (int rr=0;rr<3;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          g(rr) += xji_(rr,mm);
        }
      }

      /*         +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
      */
      const double normgsq = g(0)*g(0)+g(1)*g(1)+g(2)*g(2);

      /*
                                1.0
                  tau  = -----------------
                     C           /      \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);

    }
    else if(whichtau == Fluid3::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      //
      // tau_M: modification of
      //
      //    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
      //    Finite Element Method for the Advective-Reactive-Diffusive
      //    Equation. Computer Methods in Applied Mechanics and Enginnering,
      //    Vol. 190, pp. 1785-1800, 2000.
      //    http://www.lncc.br/~valentin/publication.htm                   */
      //
      // tau_Mp: modification of Barrenechea, G.R. and Valentin, F.
      //
      //    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
      //    element method for a generalized Stokes problem. Numerische
      //    Mathematik, Vol. 92, pp. 652-677, 2002.
      //    http://www.lncc.br/~valentin/publication.htm
      //
      //
      // tau_C: kept Wall definition
      //
      // for the modifications see Codina, Principe, Guasch, Badia
      //    "Time dependent subscales in the stabilized finite  element
      //     approximation of incompressible flow problems"
      //
      //
      // see also: Codina, R. and Soto, O.: Approximation of the incompressible
      //    Navier-Stokes equations using orthogonal subscale stabilisation
      //    and pressure segregation on anisotropic finite element meshes.
      //    Computer methods in Applied Mechanics and Engineering,
      //    Vol 193, pp. 1403-1419, 2004.

      //---------------------------------------------- compute tau_Mu = tau_Mp
      /* convective : viscous forces (element reynolds number)*/
      const double re_convectaf = (vel_normaf * hk / visceff ) * (mk/2.0);

      const double xi_convectaf = DMAX(re_convectaf,1.0);

      /*
               xi_convect ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re_convect
                              1
      */

      /* the 4.0 instead of the Franca's definition 2.0 results from the viscous
       * term in the Navier-Stokes-equations, which is scaled by 2.0*nu         */

      tau_(0) = DSQR(hk) / (4.0 * visceff / mk + ( 4.0 * visceff/mk) * xi_convectaf);

      /*------------------------------------------------------ compute tau_C ---*/

      //-- stability parameter definition according to Wall Diss. 99
      /*
               xi_convect ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> Re_convect
                              1
      */
      const double re_convectnp = (vel_normnp * hk / visceff ) * (mk/2.0);

      const double xi_tau_c = DMIN(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else if(whichtau == Fluid3::codina)
    {
      // Parameter from Codina, Badia (Constants are chosen according to
      // the values in the standard definition above)

      const double CI  = 4.0/mk;
      const double CII = 2.0/mk;

      // in contrast to the original definition, we neglect the influence of
      // the subscale velocity on velnormaf
      tau_(0)=1.0/(CI*visceff/(hk*hk)+CII*vel_normaf/hk);

      tau_(1)=tau_(0);

      tau_(2)=(hk*hk)/(CI*tau_(0));
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter\n");
    }

  } // end Fluid3::subscales_time_dependent
  else
  {
    //-------------------------------------------------------
    //        TAUS FOR THE QUASISTATIC FORMULATION
    //-------------------------------------------------------


    if(whichtau == Fluid3::bazilevs)
    {
      /* INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA

         tau_M: Bazilevs et al.
                                                               1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+af      n+af         2         |
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

      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      blitz::Array<double,2> G(3,3,blitz::ColumnMajorArray<2>());

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
          }
        }
      }

      /*          +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
      */
      double normG = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          normG+=G(nn,rr)*G(nn,rr);
        }
      }

      /*                    +----
           n+af      n+af    \     n+af         n+af
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
      */
      double Gnormu = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          Gnormu+=velintaf_(nn)*G(nn,rr)*velintaf_(rr);
        }
      }

      // definition of constant
      // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
      //  brought 144.0 from Austin...)
      const double CI = 12.0/mk;

      /*                                                       1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+af      n+af         2         |
          tau  = | --- + u     * G u     + C * nu  * G : G |
             M   |   2           -          I        -   - |
                 | dt            -                   -   - |
                 +-                                       -+
      */
      tau_(0) = 1.0/sqrt(4.0/(dt*dt)+Gnormu+CI*visceff*visceff*normG);
      tau_(1) = tau_(0);

      /*         +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
      */
      blitz::Array<double,1> g(3);

      for (int rr=0;rr<3;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          g(rr) += xji_(rr,mm);
        }
      }

      /*         +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
      */
      const double normgsq = g(0)*g(0)+g(1)*g(1)+g(2)*g(2);

      /*
                                1.0
                  tau  = -----------------
                     C           /      \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);
    }
    else if (whichtau == Fluid3::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Wall


      // this copy of velintaf_ will be used to store the normed velocity
      blitz::Array<double,1> normed_velintaf(3);
      normed_velintaf=velintaf_.copy();

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        normed_velintaf = velintaf_/vel_normaf;
      }
      else
      {
        normed_velintaf    = 0.;
        normed_velintaf(0) = 1.;
      }

      // get streamlength
      const double val = blitz::sum(blitz::abs(blitz::sum(normed_velintaf(j)*derxy_(j,i),j)));
      const double strle = 2.0/val;

      // time factor
      const double timefac = gamma*dt;

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


      const double re1 = 4.0 * timefac * visceff / (mk * DSQR(strle));   /* viscous : reactive forces   */
      const double re2 = mk * vel_normaf * strle / (2.0 * visceff);      /* convective : viscous forces */

      const double xi1 = DMAX(re1,1.0);
      const double xi2 = DMAX(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visceff);    /* convective : viscous forces */

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
      tau_(1) = timefac * DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visceff/mk) * xi_convect);

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
      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;

    }
    else if(whichtau == Fluid3::codina)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Codina


      // this copy of velintaf_ will be used to store the normed velocity
      blitz::Array<double,1> normed_velintaf(3);
      normed_velintaf=velintaf_.copy();

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        normed_velintaf = velintaf_/vel_normaf;
      }
      else
      {
        normed_velintaf    = 0.;
        normed_velintaf(0) = 1.;
      }

      // get streamlength
      const double val = blitz::sum(blitz::abs(blitz::sum(normed_velintaf(j)*derxy_(j,i),j)));
      const double strle = 2.0/val;

      // time factor
      const double timefac = gamma*dt;

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


      const double re1 = 4.0 * timefac * visceff / (mk * DSQR(strle));   /* viscous : reactive forces   */
      const double re2 = mk * vel_normaf * strle / (2.0 * visceff);      /* convective : viscous forces */

      const double xi1 = DMAX(re1,1.0);
      const double xi2 = DMAX(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visceff);    /* convective : viscous forces */

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
      tau_(1) = timefac * DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visceff/mk) * xi_convect);

      /*------------------------------------------------------ compute tau_C ---*/
      /*-- stability parameter definition according to Codina (2002), CMAME 191
       *
       * Analysis of a stabilized finite element approximation of the transient
       * convection-diffusion-reaction equation using orthogonal subscales.
       * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
       *
       * */
      tau_(2) = sqrt(DSQR(visceff)+DSQR(0.5*vel_normnp*hk));
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter\n");
    }
  }

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv == Fluid3::fssgv_artificial_all or
      fssgv == Fluid3::fssgv_artificial_small)
  {
    double fsvel_normaf = 0.0;
    if (fssgv == Fluid3::fssgv_artificial_small)
    {
      // get fine-scale velocities at element center
      fsvelintaf_ = blitz::sum(funct_(j)*fsevelaf(i,j),j);

      // get fine-scale velocity norm
      fsvel_normaf = sqrt(blitz::sum(fsvelintaf_*fsvelintaf_));
    }
    // get all-scale velocity norm
    else fsvel_normaf = vel_normaf;

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = mk * fsvel_normaf * hk / visc; /* convective : viscous forces */
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(fsvel_normaf))/(2.0*visc*xi);

  }
  else if (fssgv == Fluid3::fssgv_Smagorinsky_all or
           fssgv == Fluid3::fssgv_Smagorinsky_small or
           fssgv == Fluid3::fssgv_mixed_Smagorinsky_all or
           fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
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
      // get fine-scale or all-scale velocity (np,i) derivatives at element center
      if (fssgv == Fluid3::fssgv_Smagorinsky_small or fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
           fsvderxyaf_ = blitz::sum(derxy_(j,k)*fsevelaf(i,k),k);
      else fsvderxyaf_ = vderxyaf_;

      blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());
      epsilon = 0.5 * ( fsvderxyaf_(i,j) + fsvderxyaf_(j,i) );

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
    // Choices of the fine-scale Smagorinsky constant Cs_fs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)

    vart_ = Cs * Cs * hk * hk * rateofstrain;
  }

#ifdef PERF
  timeeletau_ref=null;
#endif

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

#ifdef PERF
  RefCountPtr<TimeMonitor> timeeletdextras_ref = rcp(new TimeMonitor(*timeeletdextras));
#endif

  // remember whether the subscale quantities have been allocated an set to zero.
  if(tds == Fluid3::subscales_time_dependent)
  {
    // if not available, the arrays for the subscale quantities have to
    // be resized and initialised to zero
    if(ele->sub_acc_old_.extent(blitz::firstDim) != 3 || ele->sub_acc_old_.extent(blitz::secondDim) != intpoints.nquad)
    {
      ele->sub_acc_old_ .resize(3,intpoints.nquad);
      ele->sub_acc_old_  = 0.;
    }
    if(ele->sub_vel_old_.extent(blitz::firstDim) != 3 || ele->sub_vel_old_.extent(blitz::secondDim) != intpoints.nquad)
    {
      ele->sub_vel_old_ .resize(3,intpoints.nquad);
      ele->sub_vel_old_  = 0.;

      ele->sub_vel_.resize(3,intpoints.nquad);
      ele->sub_vel_ = 0.;
    }
    if(ele->sub_pre_old_ .extent(blitz::firstDim) != intpoints.nquad)
    {
      ele->sub_pre_old_ .resize(intpoints.nquad);
      ele->sub_pre_old_ = 0.;

      ele->sub_pre_.resize(intpoints.nquad);
      ele->sub_pre_ = 0.;
    }
  }

  // get subscale information from element --- this is just a reference
  // to the element data
  blitz::Array<double,2> saccn (ele->sub_acc_old_);
  blitz::Array<double,2> sveln (ele->sub_vel_old_);
  blitz::Array<double,2> svelnp(ele->sub_vel_    );
  blitz::Array<double,1> spren (ele->sub_pre_old_);
  blitz::Array<double,1> sprenp(ele->sub_pre_    );


#ifdef PERF
  timeeletdextras_ref = null;
#endif


  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {

    // set gauss point coordinates
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

#ifdef PERF
    RefCountPtr<TimeMonitor> timeelederxy_ref = rcp(new TimeMonitor(*timeelederxy));
#endif

    // get values of shape functions and derivatives in the gausspoint
    DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
    }

    // get transposed Jacobian matrix and determinant
    //
    //        +-            -+ T      +-            -+
    //        | dx   dx   dx |        | dx   dy   dz |
    //        | --   --   -- |        | --   --   -- |
    //        | dr   ds   dt |        | dr   dr   dr |
    //        |              |        |              |
    //        | dy   dy   dy |        | dx   dy   dz |
    //        | --   --   -- |   =    | --   --   -- |
    //        | dr   ds   dt |        | ds   ds   ds |
    //        |              |        |              |
    //        | dz   dz   dz |        | dx   dy   dz |
    //        | --   --   -- |        | --   --   -- |
    //        | dr   ds   dt |        | dt   dt   dt |
    //        +-            -+        +-            -+
    //
    // The Jacobian is computed using the formula
    //
    //            +-----
    //   dx_j(r)   \      dN_k(r)
    //   -------  = +     ------- * (x_j)_k
    //    dr_i     /       dr_i       |
    //            +-----    |         |
    //            node k    |         |
    //                  derivative    |
    //                   of shape     |
    //                   function     |
    //                           component of
    //                          node coordinate
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
        for (int mm=1;mm<iel_;++mm)
        {
          xjm_(nn,rr)+=deriv_(nn,mm)*xyze_(rr,mm);
        }
      }
    }
    // The determinant ist computed using Sarrus's rule
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                       xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                       xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                       xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                       xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                       xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

    // check for degenerated elements
    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
    }

    // set total integration factor
    const double fac = intpoints.qwgt[iquad]*det;

    //--------------------------------------------------------------
    //             compute global first derivates
    //--------------------------------------------------------------
    /*
      Use the Jacobian and the known derivatives in element coordinate
      directions on the right hand side to compute the derivatives in
      global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+
    */

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
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<iel_;++rr)
      {
        derxy_(nn,rr)=deriv_(0,rr)*xji_(nn,0);
        for (int mm=1;mm<3;++mm)
        {
          derxy_(nn,rr)+=deriv_(mm,rr)*xji_(nn,mm);
        }
      }
    }

#ifdef PERF
    timeelederxy_ref=null;
#endif

    //--------------------------------------------------------------
    //             compute second global derivative
    //--------------------------------------------------------------

#ifdef PERF
    RefCountPtr<TimeMonitor> timeelederxy2_ref = rcp(new TimeMonitor(*timeelederxy2));
#endif

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
    if (higher_order_ele)
    {
      // initialize and zero out everything
      blitz::Array<double,2> bm(6,6,blitz::ColumnMajorArray<2>());

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
       |        | | | | 1           | | | | | 1     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 2           | | | | | 2     +-+-+-+
       |        +-+-+-+       =     +-+-+-+-+    *  | | | | .
       |        | | | | 3           | | | | | 3     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 4           | | | | | 4     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 5           | | | | | 5     +-+-+-+
       |        +-+-+-+             +-+-+-+-+       | | | | iel-1
       |                                            +-+-+-+
       |
       |        xder2               deriv2          xyze^T
       |
       |
       |                                     +-                  -+
       |  	   	    	    	     | d^2x   d^2y   d^2z |
       |  	   	    	    	     | ----   ----   ---- |
       | 	   	   	   	     | dr^2   dr^2   dr^2 |
       | 	   	   	   	     |                    |
       | 	   	   	   	     | d^2x   d^2y   d^2z |
       |                                     | ----   ----   ---- |
       | 	   	   	   	     | ds^2   ds^2   ds^2 |
       | 	   	   	   	     |                    |
       | 	   	   	   	     | d^2x   d^2y   d^2z |
       | 	   	   	   	     | ----   ----   ---- |
       | 	   	   	   	     | dt^2   dt^2   dt^2 |
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
       | 	   	   	   	     +-                  -+
       |
       |
      */
      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          xder2_(nn,rr)= deriv2_(nn,0)*xyze_(rr,0);

          for (int mm=1;mm<iel_;++mm)
          {
            xder2_(nn,rr)+= deriv2_(nn,mm)*xyze_(rr,mm);
          }
        }
      }

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
      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<iel_;++rr)
        {
          derxy2_(nn,rr)= (-1)*xder2_(nn,0)*derxy_(0,rr);

          for (int mm=1;mm<3;++mm)
          {
            derxy2_(nn,rr)-= xder2_(nn,mm)*derxy_(mm,rr);
          }
        }
      }

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

      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<iel_;++rr)
        {
          derxy2_(nn,rr) += deriv2_(nn,rr);
        }
      }

      /* make LU decomposition and solve system for all right hand sides
       * (i.e. the components of chainrulerhs)
       |
       |             0  1  2  3  4  5         i        i
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 0     | | 0    | | 0
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 1     | | 1    | | 1
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 2     | | 2    | | 2
       | 	   +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
       |           |  |  |  |  |  |  | 3     | | 3    | | 3
       |           +--+--+--+--+--+--+       +-+      +-+
       |           |  |  |  |  |  |  | 4     | | 4    | | 4
       |           +--+--+--+--+--+--+       +-+      +-+
       |           |  |  |  |  |  |  | 5     | | 5    | | 5
       |           +--+--+--+--+--+--+       +-+      +-+
       |                                      |        |
       |                                      |        |
       |                                      derxy2[i]|
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
       |    	       	      +-+-+-+-+
      */
      // Use LAPACK
      Epetra_LAPACK          solver;

      // a vector specifying the pivots (reordering)
      int pivot[6];

      // error code
      int ierr = 0;

      // Perform LU factorisation --- this call replaces bm with its factorisation
      solver.GETRF(6,6,bm.data(),6,&(pivot[0]),&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform LU factorisation during computation of derxy2");
      }

      // backward substitution. GETRS replaces the input (chainrulerhs, currently
      // stored on derxy2) with the result
      solver.GETRS('N',6,iel_,bm.data(),6,&(pivot[0]),derxy2_.data(),6,&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform backward substitution after factorisation of jacobian");
      }
    }

#ifdef PERF
    timeelederxy2_ref = null;
#endif

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------

#ifdef PERF
    RefCountPtr<TimeMonitor> timeeleintertogp_ref = rcp(new TimeMonitor(*timeeleintertogp));
#endif

    // get intermediate accelerations (n+alpha_M,i) at integration point
    //
    //                 +-----
    //       n+am       \                  n+am
    //    acc    (x) =   +      N (x) * acc
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    // i         : space dimension u/v/w
    //
    for (int rr=0;rr<3;++rr)
    {
      accintam_(rr)=funct_(0)*eaccam(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        accintam_(rr)+=funct_(mm)*eaccam(rr,mm);
      }
    }


    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    for (int rr=0;rr<3;++rr)
    {
      velintaf_(rr)=funct_(0)*evelaf(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        velintaf_(rr)+=funct_(mm)*evelaf(rr,mm);
      }
    }

    // get bodyforce in gausspoint, time (n+alpha_F)
    //
    //                 +-----
    //       n+af       \                n+af
    //      f    (x) =   +      N (x) * f
    //                  /        j       j
    //                 +-----
    //                 node j
    //
    for (int rr=0;rr<3;++rr)
    {
      bodyforceaf_(rr)=funct_(0)*edeadaf_(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        bodyforceaf_(rr)+=funct_(mm)*edeadaf_(rr,mm);
      }
    }

    // get velocities (n+1,i)  at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    vel   (x) =   +      N (x) * vel
    //                 /        j         j
    //                +-----
    //                node j
    //
    for (int rr=0;rr<3;++rr)
    {
      velintnp_(rr)=funct_(0)*evelnp(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        velintnp_(rr)+=funct_(mm)*evelnp(rr,mm);
      }
    }

    // get pressure (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    pre   (x) =   +      N (x) * pre
    //                 /        i         i
    //                +-----
    //                node i
    //
    prenp_=funct_(0)*eprenp(0);
    for (int mm=1;mm<iel_;++mm)
    {
      prenp_+=funct_(mm)*eprenp(mm);
    }

    // get pressure gradient (n+1,i) at integration point
    //
    //       n+1      +-----  dN (x)
    //   dpre   (x)    \        j         n+1
    //   ---------- =   +     ------ * pre
    //       dx        /        dx        j
    //         i      +-----      i
    //                node j
    //
    // i : direction of derivative
    //
    for (int rr=0;rr<3;++rr)
    {
      pderxynp_(rr)=derxy_(rr,0)*eprenp(0);
      for (int mm=1;mm<iel_;++mm)
      {
        pderxynp_(rr)+=derxy_(rr,mm)*eprenp(mm);
      }
    }

    // get velocity (n+alpha_F,i) derivatives at integration point
    //
    //       n+af      +-----  dN (x)
    //   dvel    (x)    \        k         n+af
    //   ----------- =   +     ------ * vel
    //       dx         /        dx        k
    //         j       +-----      j
    //                 node k
    //
    // j : direction of derivative x/y/z
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        vderxyaf_(nn,rr)=derxy_(rr,0)*evelaf(nn,0);
        for (int mm=1;mm<iel_;++mm)
        {
          vderxyaf_(nn,rr)+=derxy_(rr,mm)*evelaf(nn,mm);
        }
      }
    }

    //TODO
    // compute material and effective viscosity at every Gaussian point
    if( material->mattyp != m_fluid )
    {
    	CalVisc_AtGaussianPoint(material, visc);
    	visceff = visc + viscturb;
    }
    // get fine-scale velocity (n+alpha_F,i) derivatives at integration point
    //
    //       n+af      +-----  dN (x)
    //   dfsvel  (x)    \        k           n+af
    //   ----------- =   +     ------ * fsvel
    //       dx         /        dx          k
    //         j       +-----      j
    //                 node k
    //
    // j : direction of derivative x/y/z
    //
    // get fine-scale velocity (np,i) derivatives at integration point
    if (fssgv != Fluid3::fssgv_no  && fssgv != Fluid3::fssgv_scale_similarity)
         fsvderxyaf_ = blitz::sum(derxy_(j,k)*fsevelaf(i,k),k);
    else fsvderxyaf_ = 0.;

    // get values at integration point required for scale-similarity model
    if(fssgv == Fluid3::fssgv_scale_similarity ||
       fssgv == Fluid3::fssgv_mixed_Smagorinsky_all ||
       fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
    {
      // get coarse-scale velocities at integration point
      csvelintaf_ = blitz::sum(funct_(j)*csevelaf(i,j),j);

      // get coarse-scale velocity (np,i) derivatives at integration point
      csvderxyaf_ = blitz::sum(derxy_(j,k)*csevelaf(i,k),k);

      // PR(u) * grad PR(u): */
      convaf_s_ = blitz::sum(csvderxyaf_(j,i)*csvelintaf_(j), j);

      // get coarse-scale convective stresses at integration point
      csconvintaf_ = blitz::sum(funct_(j)*cseconvaf(i,j),j);
    }

    // get velocity (n+1,i) derivatives at integration point
    //
    //       n+1      +-----  dN (x)
    //   dvel   (x)    \        k         n+1
    //   ---------- =   +     ------ * vel
    //       dx        /        dx        k
    //         j      +-----      j
    //                node k
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        vderxynp_(nn,rr)=derxy_(rr,0)*evelnp(nn,0);
        for (int mm=1;mm<iel_;++mm)
        {
          vderxynp_(nn,rr)+=derxy_(rr,mm)*evelnp(nn,mm);
        }
      }
    }

    // get ale convective velocity (n+alpha_F,i) at integration point

    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----           +-                   -+
    //       n+af       \               |   n+af      n+alphaF|
    //      c    (x) =   +      N (x) * |vel     - u_G        |
    //                  /        j      |   j         j       |
    //                 +-----           +-                   -+
    //                 node j
    //
    //
    // where u_G is the grid velocity at the integration point,
    // time (n+alpha_F,i)
    //
    //                 +-----
    //       n+af       \                  n+af
    //    u_G    (x) =   +      N (x) * u_G
    //                  /        j         j
    //                 +-----
    //                 node j
    //

    for (int rr=0;rr<3;++rr)
    {
      aleconvintaf_(rr)=velintaf_(rr);
    }

    if (ele->is_ale_)
    {
      for (int rr=0;rr<3;++rr)
      {
        for (int mm=0;mm<iel_;++mm)
        {
          aleconvintaf_(rr)-=funct_(mm)*egridvaf(rr,mm);
        }
      }
    }


    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    for (int rr=0;rr<iel_;++rr)
    {
      conv_c_af_(rr)=aleconvintaf_(0)*derxy_(0,rr);
      for (int mm=1;mm<3;++mm)
      {
        conv_c_af_(rr)+=aleconvintaf_(mm)*derxy_(mm,rr);
      }
    }

    // calculate 2nd velocity derivatives at integration point, time(n+alpha_F)
    //
    //    2   n+af       +-----   dN (x)
    //   d vel    (x)     \         k          n+af
    //   ------------  =   +     -------- * vel
    //    dx  dx          /      dx  dx        k
    //      j1  j2       +-----    j1  j2
    //                   node k
    //
    // j=(j1,j2) : direction of derivative x/y/z
    if(higher_order_ele)
    {
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<6;++rr)
        {
          vderxy2af_(nn,rr)=derxy2_(rr,0)*evelaf(nn,0);
          for (int mm=1;mm<iel_;++mm)
          {
            vderxy2af_(nn,rr)+=derxy2_(rr,mm)*evelaf(nn,mm);
          }
        }
      }

    }


    //TODO



    if (higher_order_ele)
    {
      /*--- viscous term  grad * epsilon(u): ------------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
         - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
         2 |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                             */
      for (int mm=0;mm<iel_;++mm)
      {
        viscs2_(0,1,mm) = 0.5 *  derxy2_(3,mm);
        viscs2_(1,0,mm) = viscs2_(0,1,mm)     ;
      }
      for (int mm=0;mm<iel_;++mm)
      {
        viscs2_(0,2,mm) = 0.5 *  derxy2_(4,mm);
        viscs2_(2,0,mm) = viscs2_(0,2,mm);
      }
      for (int mm=0;mm<iel_;++mm)
      {
        viscs2_(1,2,mm) = 0.5 *  derxy2_(5,mm);
        viscs2_(2,1,mm) = viscs2_(1,2,mm);
      }
      for (int mm=0;mm<iel_;++mm)
      {
        const double  sum = 0.5 * (derxy2_(0,mm) + derxy2_(1,mm) + derxy2_(2,mm));

        viscs2_(0,0,mm) = sum + 0.5 * derxy2_(0,mm);
        viscs2_(1,1,mm) = sum + 0.5 * derxy2_(1,mm);
        viscs2_(2,2,mm) = sum + 0.5 * derxy2_(2,mm);
      }
    }

    /* divergence new time step n+1 */
    const double divunp = (vderxynp_(0,0)+vderxynp_(1,1)+vderxynp_(2,2));

    /* Convective term  u_old * grad u_old: */
    for (int rr=0;rr<3;++rr)
    {
      convaf_old_(rr)=vderxyaf_(rr, 0)*aleconvintaf_(0);
      for (int mm=1;mm<3;++mm)
      {
        convaf_old_(rr)+=vderxyaf_(rr,mm)*aleconvintaf_(mm);
      }
    }

    if (higher_order_ele)
    {
      /* Viscous term  div epsilon(u_old) */
      viscaf_old_(0) = vderxy2af_(0,0) + 0.5 * (vderxy2af_(0,1) + vderxy2af_(1,3) + vderxy2af_(0,2) + vderxy2af_(2,4));
      viscaf_old_(1) = vderxy2af_(1,1) + 0.5 * (vderxy2af_(1,0) + vderxy2af_(0,3) + vderxy2af_(1,2) + vderxy2af_(2,5));
      viscaf_old_(2) = vderxy2af_(2,2) + 0.5 * (vderxy2af_(2,0) + vderxy2af_(0,4) + vderxy2af_(2,1) + vderxy2af_(1,5));
    }

    /* compute residual in gausspoint --- second derivatives only
     * exist for higher order elements,  convaf_old_ is based on
     * ale-convective velocity */
    for (int rr=0;rr<3;++rr)
    {
      resM_(rr) = accintam_(rr) + convaf_old_(rr) + pderxynp_(rr) - bodyforceaf_(rr);
    }

    /* the residual is based on the effective viscosity! */

    if (higher_order_ele)
    {
      const double two_visceff=2.0*visceff;

      for (int rr=0;rr<3;++rr)
      {
        resM_(rr) -= two_visceff*viscaf_old_(rr);
      }
    }

    /*
      This is the operator

                  /               \
                 | resM    o nabla |
                  \    (i)        /

      required for the cross and reynolds stress calculation

    */
    for (int rr=0;rr<iel_;++rr)
    {
      conv_resM_(rr)=resM_(0)*derxy_(0,rr);
      for (int mm=1;mm<3;++mm)
      {
        conv_resM_(rr)+=resM_(mm)*derxy_(mm,rr);
      }
    }

#ifdef PERF
    timeeleintertogp_ref = null;
#endif

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //
    //    ELEMENT FORMULATION BASED ON TIME DEPENDENT SUBSCALES
    //
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    if(tds == Fluid3::subscales_time_dependent)
    {
      const double tauM   = tau_(0);
      const double tauC   = tau_(2);

#ifdef PERF
      timeeletdextras_ref = rcp(new TimeMonitor(*timeeletdextras));
#endif

      // update estimates for the subscale quantities

      const double factauC                  = tauC/(tauC+dt);
      const double facMtau                  = 1./(alphaM*tauM+afgdt);

      /*-------------------------------------------------------------------*
       *                                                                   *
       *                  update of SUBSCALE PRESSURE                      *
       *                                                                   *
       *-------------------------------------------------------------------*/

      /*
        ~n+1      tauC     ~n   tauC * dt            n+1
        p    = --------- * p  - --------- * nabla o u
         (i)   tauC + dt        tauC + dt            (i)
      */
      sprenp(iquad)=(spren(iquad)-dt*divunp)*factauC;

      /*-------------------------------------------------------------------*
       *                                                                   *
       *                  update of SUBSCALE VELOCITY                      *
       *                                                                   *
       *-------------------------------------------------------------------*/

      /*
        ~n+1                1.0
        u    = ----------------------------- *
         (i)   alpha_M*tauM+alpha_F*gamma*dt

                +-
                | +-                                  -+   ~n
               *| |alpha_M*tauM +gamma*dt*(alpha_F-1.0)| * u +
                | +-                                  -+
                +-


                    +-                      -+    ~ n
                  + | dt*tauM*(alphaM-gamma) | * acc -
                    +-                      -+

                                           -+
                                       n+1  |
                  - gamma*dt*tauM * res     |
                                       (i)  |
                                           -+
      */
      {
        const double fac1=(alphaM*tauM+gamma*dt*(alphaF-1.0))*facMtau;
        const double fac2=(dt*tauM*(alphaM-gamma))*facMtau;
        const double fac3=(gamma*dt*tauM)*facMtau;


        for (int rr=0;rr<3;++rr)
        {
          svelnp(rr,iquad)= fac1*sveln(rr,iquad)+fac2*saccn(rr,iquad)-fac3*resM_(rr);
        }
      }

      /*-------------------------------------------------------------------*
       *                                                                   *
       *               update of intermediate quantities                   *
       *                                                                   *
       *-------------------------------------------------------------------*/

      /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

      */
      for (int rr=0;rr<3;++rr)
      {
        svelaf_(rr)=alphaF*svelnp(rr,iquad)+(1.0-alphaF)*sveln(rr,iquad) ;
      }

      /* the intermediate value of subscale acceleration is not needed to be
       * computed anymore --- we use the governing ODE to replace it ....

             ~ n+am    alphaM     / ~n+1   ~n \    gamma - alphaM    ~ n
            acc     = -------- * |  u    - u   | + -------------- * acc
               (i)    gamma*dt    \  (i)      /         gamma

      */

      /*
        This is the operator

                  /~n+af         \
                 | u      o nabla |
                  \   (i)        /

                  required for the cross and reynolds stress calculation

      */

      for (int rr=0;rr<iel_;++rr)
      {
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);

        for (int mm=1;mm<3;++mm)
        {
               conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
        }
      }

      /* Most recent value for subgrid velocity convective term

                  /~n+af         \   n+af
                 | u      o nabla | u
                  \   (i)        /   (i)
      */
      for (int rr=0;rr<3;++rr)
      {
        convsubaf_old_(rr) = vderxyaf_(rr, 0)*svelaf_(0);

        for (int mm=1;mm<3;++mm)
        {
          convsubaf_old_(rr) += vderxyaf_(rr, mm)*svelaf_(mm);
        }
      }

#ifdef PERF
      timeeletdextras_ref = null;
#endif


      //--------------------------------------------------------------
      //--------------------------------------------------------------
      //
      //                       SYSTEM MATRIX
      //
      //--------------------------------------------------------------
      //--------------------------------------------------------------
      if(compute_elemat)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelegalerkin_ref = rcp(new TimeMonitor(*timeelegalerkin));
#endif

        //---------------------------------------------------------------
        //
        //              TIME-DEPENDENT SUBGRID-SCALES
        //   GALERKIN PART 1 AND SUBSCALE ACCELERATION STABILISATION
        //
        //---------------------------------------------------------------
        if(inertia == Fluid3::inertia_stab_keep)
        {
          const double fac_alphaM_tauM_facMtau    = fac*alphaM*tauM*facMtau;
          const double fac_afgdt_afgdt_facMtau    = fac*afgdt*afgdt*facMtau;

          {
            const double fac_alphaM_afgdt_facMtau
              =
              fac*alphaM*afgdt*facMtau;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double inertia_and_conv_ui = fac_alphaM_afgdt_facMtau*funct_(ui)+fac_afgdt_afgdt_facMtau*conv_c_af_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                const double inertia_and_conv_ui_funct_vi=inertia_and_conv_ui*funct_(vi);

                /*
                  inertia term (intermediate)
                  +convection (intermediate)


                                                          /          \
                                alphaF*gamma*dt          |            |
                    alphaM*--------------------------- * |  Dacc , v  | +
                           alphaM*tauM+alphaF*gamma*dt   |            |
                                                          \          /


                                                                  /                          \
                                        alphaF*gamma*dt          |  / n+af       \            |
                  +alphaF*gamma*dt*--------------------------- * | | c    o nabla | Dacc , v  |
                                   alphaM*tauM+alphaF*gamma*dt   |  \            /            |
                                                                  \                          /
                */

                elemat(vi*4    , ui*4    ) += inertia_and_conv_ui_funct_vi;
                elemat(vi*4 + 1, ui*4 + 1) += inertia_and_conv_ui_funct_vi;
                elemat(vi*4 + 2, ui*4 + 2) += inertia_and_conv_ui_funct_vi;

              } // end loop rows (test functions for matrix)
            } // end loop rows (solution for matrix, test function for vector)
          }


          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {
            const double fac_alphaM_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*funct_(vi);

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              /* pressure (implicit) */

              /*  factor:
                             alphaM*tauM
                    ---------------------------
                    alphaM*tauM+alphaF*gamma*dt

                 /               \
                |                 |
                |  nabla Dp ,  v  |
                |                 |
                 \               /
              */
              elemat(vi*4    , ui*4 + 3) -= fac_alphaM_tauM_facMtau_funct_vi*derxy_(0,ui);
              elemat(vi*4 + 1, ui*4 + 3) -= fac_alphaM_tauM_facMtau_funct_vi*derxy_(1,ui);
              elemat(vi*4 + 2, ui*4 + 3) -= fac_alphaM_tauM_facMtau_funct_vi*derxy_(2,ui);
            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

          if(higher_order_ele)
          {
            const double fac_two_visceff_afgdt_alphaM_tauM_facMtau
              =
              fac*2.0*visceff*afgdt*alphaM*tauM*facMtau;

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
                =
                fac_two_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                /* viscous term (intermediate) */
                /*  factor:
                                                 alphaM*tauM
                      2*nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                  /                         \
                 |               /    \      |
                 |  nabla o eps | Dacc | , v |
                 |               \    /      |
                  \                         /

                */
                elemat(vi*4    , ui*4    ) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,0,ui);
                elemat(vi*4    , ui*4 + 1) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,1,ui);
                elemat(vi*4    , ui*4 + 2) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,2,ui);
                elemat(vi*4 + 1, ui*4    ) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,1,ui);
                elemat(vi*4 + 1, ui*4 + 1) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,1,ui);
                elemat(vi*4 + 1, ui*4 + 2) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,2,ui);
                elemat(vi*4 + 2, ui*4    ) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,2,ui);
                elemat(vi*4 + 2, ui*4 + 1) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,2,ui);
                elemat(vi*4 + 2, ui*4 + 2) += fac_two_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(2,2,ui);
              } // end loop rows (test functions for matrix)
            } // end loop rows (solution for matrix, test function for vector)
          }

          if (newton) // if inertia and newton
          {
            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double fac_afgdt_afgdt_facMtau_funct_ui
                =
                fac_afgdt_afgdt_facMtau*funct_(ui);
              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                const double fac_afgdt_afgdt_facMtau_funct_ui_funct_vi
                  =
                  fac_afgdt_afgdt_facMtau_funct_ui*funct_(vi);


                /* convection (intermediate)

                factor:

                                     alphaF*gamma*dt
               +alphaF*gamma*dt*---------------------------
                                alphaM*tauM+alphaF*gamma*dt

                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
                */


                elemat(vi*4    , ui*4    ) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(0,0);
                elemat(vi*4    , ui*4 + 1) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(0,1);
                elemat(vi*4    , ui*4 + 2) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(0,2);
                elemat(vi*4 + 1, ui*4    ) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(1,0);
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(1,1);
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(1,2);
                elemat(vi*4 + 2, ui*4    ) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(2,0);
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(2,1);
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_afgdt_facMtau_funct_ui_funct_vi*vderxyaf_(2,2);
              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)
          } // end if inertia and newton
        } //   end if inertia stabilisation
        else
        {
          //---------------------------------------------------------------
          //
          //              TIME-DEPENDENT SUBGRID-SCALES
          //    GALERKIN PART I WITH NEGLECTED SUBSCALE INERTIA TERM
          //
          //---------------------------------------------------------------

          // if no inertia stabilisation
          const double fac_alphaM = fac*alphaM;
          const double fac_afgdt  = fac*afgdt;


          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double inertia_and_conv_ui_ = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double inertia_and_conv_ui_funct_vi = inertia_and_conv_ui_*funct_(vi);
              /*

                inertia term (intermediate)

                                   /          \
                                  |            |
                          +alphaM*|  Dacc , v  |
                                  |            |
                                   \          /

                convective term (intermediate)


                                   /                          \
                                  |  / n+af       \            |
                 +alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                                  |  \            /            |
                                   \                          /

              */
              elemat(vi*4    , ui*4    ) += inertia_and_conv_ui_funct_vi;
              elemat(vi*4 + 1, ui*4 + 1) += inertia_and_conv_ui_funct_vi;
              elemat(vi*4 + 2, ui*4 + 2) += inertia_and_conv_ui_funct_vi;

            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton) // if no inertia and newton
          {

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double fac_afgdt_funct_ui = fac_afgdt*funct_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                const double fac_afgdt_funct_ui_funct_vi = fac_afgdt_funct_ui*funct_(vi);

                /*  factor:

                +alphaF*gamma*dt

                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
                */
                elemat(vi*4    , ui*4    ) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(0,0);
                elemat(vi*4    , ui*4 + 1) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(0,1);
                elemat(vi*4    , ui*4 + 2) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(0,2);
                elemat(vi*4 + 1, ui*4    ) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(1,0);
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(1,1);
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(1,2);
                elemat(vi*4 + 2, ui*4    ) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(2,0);
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(2,1);
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_funct_ui_funct_vi*vderxyaf_(2,2);

              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)
          } // end if no inertia and newton
        } // end if no inertia stabilisation


        const double fac_afgdt_visceff        = fac*visceff*afgdt;
        const double fac_gamma_dt             = fac*gamma*dt;

        //---------------------------------------------------------------
        //
        //              TIME-DEPENDENT SUBGRID-SCALES
        //        GALERKIN PART 2 (REMAINING EXPRESSIONS)
        //
        //---------------------------------------------------------------
        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const double fac_funct_ui=fac*funct_(ui);

          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {

            /* pressure (implicit) */

            /*  factor: -1

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
            */

            elemat(vi*4    , ui*4 + 3) -= fac_funct_ui*derxy_(0, vi) ;
            elemat(vi*4 + 1, ui*4 + 3) -= fac_funct_ui*derxy_(1, vi) ;
            elemat(vi*4 + 2, ui*4 + 3) -= fac_funct_ui*derxy_(2, vi) ;
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {
            const double sum = derxy_(0,ui)*derxy_(0,vi)+derxy_(1,ui)*derxy_(1,vi)+derxy_(2,ui)*derxy_(2,vi);

            /* viscous term (intermediate) */

            /*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
            */

            elemat(vi*4    , ui*4    ) += fac_afgdt_visceff*(sum+derxy_(0,ui)*derxy_(0,vi));
            elemat(vi*4    , ui*4 + 1) += fac_afgdt_visceff*derxy_(0,ui)*derxy_(1,vi) ;
            elemat(vi*4    , ui*4 + 2) += fac_afgdt_visceff*derxy_(0,ui)*derxy_(2,vi) ;
            elemat(vi*4 + 1, ui*4    ) += fac_afgdt_visceff*derxy_(1,ui)*derxy_(0,vi) ;
            elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_visceff*(sum+derxy_(1,ui)*derxy_(1,vi));
            elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_visceff*derxy_(1,ui)*derxy_(2,vi) ;
            elemat(vi*4 + 2, ui*4    ) += fac_afgdt_visceff*derxy_(2,ui)*derxy_(0,vi) ;
            elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_visceff*derxy_(2,ui)*derxy_(1,vi) ;
            elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_visceff*(sum+derxy_(2,ui)*derxy_(2,vi));

          } // end loop rows (test functions for matrix)
        }

        for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
        {
          const double fac_gamma_dt_funct_vi= fac_gamma_dt*funct_(vi);

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            /* continuity equation (implicit) */

            /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
            */

            elemat(vi*4 + 3, ui*4    ) += fac_gamma_dt_funct_vi*derxy_(0,ui);
            elemat(vi*4 + 3, ui*4 + 1) += fac_gamma_dt_funct_vi*derxy_(1,ui);
            elemat(vi*4 + 3, ui*4 + 2) += fac_gamma_dt_funct_vi*derxy_(2,ui);


          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)
        // end remaining Galerkin terms

        //---------------------------------------------------------------
        //
        //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
        //                    PRESSURE STABILISATION
        //
        //---------------------------------------------------------------
        if(pspg == Fluid3::pstab_use_pspg)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelepspg_ref = rcp(new TimeMonitor(*timeelepspg));
#endif
          const double fac_gamma_dt_tauM_facMtau                   = fac*gamma*dt*tauM*facMtau;

          const double fac_afgdt_gamma_dt_tauM_facMtau             = fac*afgdt*gamma*dt*tauM*facMtau;

          const double fac_alphaM_gamma_dt_tauM_facMtau            = fac*alphaM*gamma*dt*tauM*facMtau;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)+fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /* pressure stabilisation --- inertia    */

              /*
                           gamma*dt*tau_M
                     ------------------------------ * alpha_M *
                     alpha_M*tau_M+alpha_F*gamma*dt


                                /                \
                               |                  |
                             * |  Dacc , nabla q  | +
                               |                  |
                                \                /

               pressure stabilisation --- convection


                             gamma*dt*tau_M
                   + ------------------------------ * alpha_F*gamma*dt *
                     alpha_M*tau_M+alpha_F*gamma*dt


                        /                                \
                       |  / n+af       \                  |
                     * | | c    o nabla | Dacc , nabla q  |
                       |  \            /                  |
                        \                                /
              */
              elemat(vi*4 + 3, ui*4    ) += inertia_and_conv_ui*derxy_(0,vi) ;
              elemat(vi*4 + 3, ui*4 + 1) += inertia_and_conv_ui*derxy_(1,vi) ;
              elemat(vi*4 + 3, ui*4 + 2) += inertia_and_conv_ui*derxy_(2,vi) ;

            }
          }
          if(higher_order_ele)
          {
            const double fac_two_visceff_afgdt_gamma_dt_tauM_facMtau
              =
              fac*2.0*visceff*afgdt*gamma*dt*tauM*facMtau;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                /* pressure stabilisation --- diffusion  */


            /*
                           gamma*dt*tau_M
            factor:  ------------------------------ * alpha_F*gamma*dt * 2 * nu
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                                \
                   |               /    \             |
                   |  nabla o eps | Dacc | , nabla q  |
                   |               \    /             |
                    \                                /
            */

                elemat(vi*4 + 3, ui*4    ) -= fac_two_visceff_afgdt_gamma_dt_tauM_facMtau*
                                              (derxy_(0,vi)*viscs2_(0,0,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(0,1,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(0,2,ui)) ;
                elemat(vi*4 + 3, ui*4 + 1) -= fac_two_visceff_afgdt_gamma_dt_tauM_facMtau*
                                              (derxy_(0,vi)*viscs2_(0,1,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(1,1,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(1,2,ui)) ;
                elemat(vi*4 + 3, ui*4 + 2) -= fac_two_visceff_afgdt_gamma_dt_tauM_facMtau*
                                              (derxy_(0,vi)*viscs2_(0,2,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(1,2,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(2,2,ui)) ;
              }
            }
          }
          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
                /* pressure stabilisation --- pressure   */

              /*
                          gamma*dt*tau_M
            factor:  ------------------------------
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
              */

              elemat(vi*4 + 3, ui*4 + 3) += fac_gamma_dt_tauM_facMtau*
                                            (derxy_(0,ui)*derxy_(0,vi)
                                             +
                                             derxy_(1,ui)*derxy_(1,vi)
                                             +
                                             derxy_(2,ui)*derxy_(2,vi)) ;

            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton) // if pspg and newton
          {

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double fac_afgdt_gamma_dt_tauM_facMtau_funct_ui
                =
                fac_afgdt_gamma_dt_tauM_facMtau*funct_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                /* pressure stabilisation --- convection */

                /*
                                gamma*dt*tau_M
                factor:  ------------------------------ * alpha_F*gamma*dt
                         alpha_M*tau_M+alpha_F*gamma*dt

                       /                                  \
                      |  /            \   n+af             |
                      | | Dacc o nabla | u      , nabla q  |
                      |  \            /                    |
                       \                                  /

                */

                elemat(vi*4 + 3, ui*4    ) += fac_afgdt_gamma_dt_tauM_facMtau_funct_ui*
                                              (derxy_(0,vi)*vderxyaf_(0,0)
                                               +
                                               derxy_(1,vi)*vderxyaf_(1,0)
                                               +
                                               derxy_(2,vi)*vderxyaf_(2,0)) ;
                elemat(vi*4 + 3, ui*4 + 1) += fac_afgdt_gamma_dt_tauM_facMtau_funct_ui*
                                              (derxy_(0,vi)*vderxyaf_(0,1)
                                               +
                                               derxy_(1,vi)*vderxyaf_(1,1)
                                               +
                                               derxy_(2,vi)*vderxyaf_(2,1)) ;
                elemat(vi*4 + 3, ui*4 + 2) += fac_afgdt_gamma_dt_tauM_facMtau_funct_ui*
                                              (derxy_(0,vi)*vderxyaf_(0,2)
                                               +
                                               derxy_(1,vi)*vderxyaf_(1,2)
                                               +
                                               derxy_(2,vi)*vderxyaf_(2,2)) ;
              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)
          }// end if pspg and newton

#ifdef PERF
          timeelepspg_ref = null;
#endif

        } // end pressure stabilisation


        //---------------------------------------------------------------
        //
        //         STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
        //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
        //
        //---------------------------------------------------------------
        if(supg == Fluid3::convective_stab_supg)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelesupg_ref = rcp(new TimeMonitor(*timeelesupg));
#endif

          const double fac_alphaM_afgdt_tauM_facMtau            = fac*alphaM*afgdt*facMtau*tauM;
          const double fac_afgdt_tauM_afgdt_facMtau             = fac*afgdt*afgdt*facMtau*tauM;
          const double fac_afgdt_tauM_facMtau                   = fac*afgdt*tauM*facMtau;
          const double fac_two_visceff_afgdt_afgdt_tauM_facMtau = fac*2.0*visceff*afgdt*afgdt*tauM*facMtau;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double inertia_and_conv_ui
              =
              fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_tauM_afgdt_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double inertia_and_conv_ui_conv_c_af_vi
                = inertia_and_conv_ui*conv_c_af_(vi);

              /* SUPG stabilisation --- inertia

               factor:
                           alphaF*gamma*dt*tauM
                        --------------------------- * alphaM
                        alphaM*tauM+alphaF*gamma*dt


                    /                           \
                   |          / n+af       \     |
                   |  Dacc , | c    o nabla | v  |
                   |          \            /     |
                    \                           /
              */

              /* SUPG stabilisation --- convection


               factor:
                           alphaF*gamma*dt*tauM
                        --------------------------- * alphaF * gamma * dt
                        alphaM*tauM+alphaF*gamma*dt

                    /                                               \
                   |    / n+af        \          / n+af        \     |
                   |   | c     o nabla | Dacc , | c     o nabla | v  |
                   |    \             /          \             /     |
                    \                                               /
              */
              elemat(vi*4    , ui*4    ) += inertia_and_conv_ui_conv_c_af_vi ;
              elemat(vi*4 + 1, ui*4 + 1) += inertia_and_conv_ui_conv_c_af_vi ;
              elemat(vi*4 + 2, ui*4 + 2) += inertia_and_conv_ui_conv_c_af_vi ;
            }
          }

          if(higher_order_ele)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi=
                fac_two_visceff_afgdt_afgdt_tauM_facMtau*conv_c_af_(vi);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                /* SUPG stabilisation --- diffusion

                factor:
                               alphaF*gamma*tauM*dt
                  - 2 * nu  --------------------------- * alphaF * gamma * dt
                            alphaM*tauM+alphaF*gamma*dt


                    /                                            \
                   |               /     \    / n+af        \     |
                   |  nabla o eps | Dacc  |, | c     o nabla | v  |
                   |               \     /    \             /     |
                    \                                            /
                */
                elemat(vi*4, ui*4)         -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(0, 0, ui);
                elemat(vi*4, ui*4 + 1)     -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(0, 1, ui);
                elemat(vi*4, ui*4 + 2)     -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(0, 2, ui);
                elemat(vi*4 + 1, ui*4)     -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(0, 1, ui);
                elemat(vi*4 + 1, ui*4 + 1) -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(1, 1, ui);
                elemat(vi*4 + 1, ui*4 + 2) -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(1, 2, ui);
                elemat(vi*4 + 2, ui*4)     -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(0, 2, ui);
                elemat(vi*4 + 2, ui*4 + 1) -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(1, 2, ui);
                elemat(vi*4 + 2, ui*4 + 2) -= fac_two_visceff_afgdt_afgdt_tauM_facMtau_conv_c_af_vi*viscs2_(2, 2, ui);
              }
            }
          }

          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {
            const double fac_afgdt_tauM_facMtau_conv_c_af_vi
              =
              fac_afgdt_tauM_facMtau*conv_c_af_(vi);

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
                /* SUPG stabilisation --- pressure

               factor:
                               alphaF*gamma*tauM*dt
                            ---------------------------
                            alphaM*tauM+alphaF*gamma*dt


                    /                               \
                   |              / n+af       \     |
                   |  nabla Dp , | c    o nabla | v  |
                   |              \            /     |
                    \                               /
              */
              elemat(vi*4    , ui*4 + 3) += fac_afgdt_tauM_facMtau_conv_c_af_vi*derxy_(0,ui);
              elemat(vi*4 + 1, ui*4 + 3) += fac_afgdt_tauM_facMtau_conv_c_af_vi*derxy_(1,ui);
              elemat(vi*4 + 2, ui*4 + 3) += fac_afgdt_tauM_facMtau_conv_c_af_vi*derxy_(2,ui);

            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton)
          {
            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double fac_afgdt_tauM_afgdt_facMtau_funct_ui
                =
                fac_afgdt_tauM_afgdt_facMtau*funct_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                const double fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi
                  =
                  fac_afgdt_tauM_afgdt_facMtau_funct_ui*conv_c_af_(vi);

                /* SUPG stabilisation --- convection


               factor:
                           alphaF*gamma*dt*tauM
                        --------------------------- * alphaF * gamma * dt
                        alphaM*tauM+alphaF*gamma*dt

                    /                                               \
                   |    /            \   n+af    / n+af        \     |
                   |   | Dacc o nabla | u     , | c     o nabla | v  |
                   |    \            /           \             /     |
                    \                                               /

                */

                elemat(vi*4, ui*4)         += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(0,0);
                elemat(vi*4, ui*4 + 1)     += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(0,1);
                elemat(vi*4, ui*4 + 2)     += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(0,2);
                elemat(vi*4 + 1, ui*4)     += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(1,0);
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(1,1);
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(1,2);
                elemat(vi*4 + 2, ui*4)     += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(2,0);
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(2,1);
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_tauM_afgdt_facMtau_funct_ui_conv_c_af_vi*vderxyaf_(2,2);
              }
            }

            const double fac_afgdt_svelaf_x = fac*afgdt*svelaf_(0);
            const double fac_afgdt_svelaf_y = fac*afgdt*svelaf_(1);
            const double fac_afgdt_svelaf_z = fac*afgdt*svelaf_(2);


            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {

              const double fac_afgdt_svelaf_x_funct_ui = fac_afgdt_svelaf_x*funct_(ui);
              const double fac_afgdt_svelaf_y_funct_ui = fac_afgdt_svelaf_y*funct_(ui);
              const double fac_afgdt_svelaf_z_funct_ui = fac_afgdt_svelaf_z*funct_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {

                /* SUPG stabilisation --- subscale velocity, nonlinear part from testfunction

                factor:
                          alphaF * gamma * dt


                    /                            \
                   |  ~n+af    /            \     |
                   |  u     , | Dacc o nabla | v  |
                   |   (i)     \            /     |
                    \                            /

                */

                elemat(vi*4    , ui*4)     -= fac_afgdt_svelaf_x_funct_ui*derxy_(0,vi) ;
                elemat(vi*4    , ui*4 + 1) -= fac_afgdt_svelaf_x_funct_ui*derxy_(1,vi) ;
                elemat(vi*4    , ui*4 + 2) -= fac_afgdt_svelaf_x_funct_ui*derxy_(2,vi) ;
                elemat(vi*4 + 1, ui*4)     -= fac_afgdt_svelaf_y_funct_ui*derxy_(0,vi) ;
                elemat(vi*4 + 1, ui*4 + 1) -= fac_afgdt_svelaf_y_funct_ui*derxy_(1,vi) ;
                elemat(vi*4 + 1, ui*4 + 2) -= fac_afgdt_svelaf_y_funct_ui*derxy_(2,vi) ;
                elemat(vi*4 + 2, ui*4)     -= fac_afgdt_svelaf_z_funct_ui*derxy_(0,vi) ;
                elemat(vi*4 + 2, ui*4 + 1) -= fac_afgdt_svelaf_z_funct_ui*derxy_(1,vi) ;
                elemat(vi*4 + 2, ui*4 + 2) -= fac_afgdt_svelaf_z_funct_ui*derxy_(2,vi) ;
              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)
          } // end if newton and supg
#ifdef PERF
          timeelesupg_ref = null;
#endif

        } // end supg stabilisation



        //---------------------------------------------------------------
        //
        //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
        //            VISCOUS STABILISATION TERMS FOR (A)GLS
        //
        //---------------------------------------------------------------
        if (higher_order_ele)
        {
          if(vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_gls)
          {
#ifdef PERF
            RefCountPtr<TimeMonitor> timeelevstab_ref = rcp(new TimeMonitor(*timeelevstab));
#endif
            const double fac_alphaM_two_visc_afgdt_tauM_facMtau         = vstabfac*fac*alphaM*2.0*visc*afgdt*tauM*facMtau;
            const double fac_afgdt_two_visc_afgdt_tauM_facMtau          = vstabfac*fac*afgdt*2.0*visc*afgdt*tauM*facMtau;
            const double fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau = vstabfac*fac*afgdt*4.0*visceff*visc*afgdt*tauM*facMtau;
            const double fac_two_visc_afgdt_tauM_facMtau                = vstabfac*fac*2.0*visc*afgdt*tauM*facMtau;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const double fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui
                =
                fac_alphaM_two_visc_afgdt_tauM_facMtau*funct_(ui);

              const double fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui
                =
                fac_afgdt_two_visc_afgdt_tauM_facMtau*conv_c_af_(ui);

              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                /* viscous stabilisation --- inertia     */

                /* factor:

                                        alphaF*gamma*tauM*dt
                   +(-)alphaM*2*nu* ---------------------------
                                    alphaM*tauM+alphaF*gamma*dt

                     /                    \
                    |                      |
                    |  Dacc , div eps (v)  |
                    |                      |
                     \                    /
                */
                elemat(vi*4    , ui*4    ) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(0,0,vi);
                elemat(vi*4    , ui*4 + 1) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(1,0,vi);
                elemat(vi*4    , ui*4 + 2) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(2,0,vi);
                elemat(vi*4 + 1, ui*4    ) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(0,1,vi);
                elemat(vi*4 + 1, ui*4 + 1) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(1,1,vi);
                elemat(vi*4 + 1, ui*4 + 2) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(2,1,vi);
                elemat(vi*4 + 2, ui*4    ) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(0,2,vi);
                elemat(vi*4 + 2, ui*4 + 1) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(1,2,vi);
                elemat(vi*4 + 2, ui*4 + 2) += fac_alphaM_two_visc_afgdt_tauM_facMtau_funct_ui*viscs2_(2,2,vi);

                /* viscous stabilisation --- convection */
                /*  factor:
                                         alphaF*gamma*dt*tauM
            +(-)alphaF*gamma*dt*2*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                       /                                  \
                      |  / n+af       \                    |
                      | | c    o nabla | Dacc, div eps (v) |
                      |  \            /                    |
                       \                                  /

                */

                elemat(vi*4    , ui*4    ) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(0, 0, vi) ;
                elemat(vi*4    , ui*4 + 1) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(1, 0, vi) ;
                elemat(vi*4    , ui*4 + 2) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(2, 0, vi) ;
                elemat(vi*4 + 1, ui*4    ) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(0, 1, vi) ;
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(1, 1, vi) ;
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(2, 1, vi) ;
                elemat(vi*4 + 2, ui*4    ) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(0, 2, vi) ;
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(1, 2, vi) ;
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_two_visc_afgdt_tauM_facMtau_conv_c_af_ui*viscs2_(2, 2, vi) ;

                /* viscous stabilisation --- diffusion  */

              /* factor:

                                           alphaF*gamma*tauM*dt
            -(+)alphaF*gamma*dt*4*nu*nu ---------------------------
                                        alphaM*tauM+alphaF*gamma*dt

                    /                                   \
                   |               /    \                |
                   |  nabla o eps | Dacc | , div eps (v) |
                   |               \    /                |
                    \                                   /
              */
                elemat(vi*4    , ui*4    ) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 0, vi)
                                               +
                                               viscs2_(1, 0, ui)*viscs2_(1, 0, vi)
                                               +
                                               viscs2_(2, 0, ui)*viscs2_(2, 0, vi)) ;

                elemat(vi*4    , ui*4 + 1) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 0, vi)*viscs2_(0, 1, ui)
                                               +
                                               viscs2_(1, 0, vi)*viscs2_(1, 1, ui)
                                               +
                                               viscs2_(2, 0, vi)*viscs2_(2, 1, ui)) ;

                elemat(vi*4    , ui*4 + 2) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 0, vi)*viscs2_(0, 2, ui)
                                               +
                                               viscs2_(1, 0, vi)*viscs2_(1, 2, ui)
                                               +
                                               viscs2_(2, 0, vi)*viscs2_(2, 2, ui)) ;

                elemat(vi*4 + 1, ui*4    ) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 1, vi)
                                               +
                                               viscs2_(1, 0, ui)*viscs2_(1, 1, vi)
                                               +
                                               viscs2_(2, 0, ui)*viscs2_(2, 1, vi)) ;

                elemat(vi*4 + 1, ui*4 + 1) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 1, ui)*viscs2_(0, 1, vi)
                                               +
                                               viscs2_(1, 1, ui)*viscs2_(1, 1, vi)
                                               +
                                               viscs2_(2, 1, ui)*viscs2_(2, 1, vi)) ;

                elemat(vi*4 + 1, ui*4 + 2) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 1, vi)*viscs2_(0, 2, ui)
                                               +
                                               viscs2_(1, 1, vi)*viscs2_(1, 2, ui)
                                               +
                                               viscs2_(2, 1, vi)*viscs2_(2, 2, ui)) ;

                elemat(vi*4 + 2, ui*4    ) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(1, 0, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(2, 0, ui)*viscs2_(2, 2, vi)) ;
                elemat(vi*4 + 2, ui*4 + 1) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 1, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(1, 1, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(2, 1, ui)*viscs2_(2, 2, vi)) ;

                elemat(vi*4 + 2, ui*4 + 2) -= fac_afgdt_four_visceff_visc_afgdt_tauM_facMtau*
                                              (viscs2_(0, 2, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(1, 2, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

                /* viscous stabilisation --- pressure   */

                /* factor:

                                    alphaF*gamma*tauM*dt
                     +(-)2*nu * ---------------------------
                                alphaM*tauM+alphaF*gamma*dt


                    /                        \
                   |                          |
                   |  nabla Dp , div eps (v)  |
                   |                          |
                    \                        /
                */
                elemat(vi*4    , ui*4 + 3) += fac_two_visc_afgdt_tauM_facMtau*
                                              (derxy_(0,ui)*viscs2_(0,0,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(1,0,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(2,0,vi)) ;
                elemat(vi*4 + 1, ui*4 + 3) += fac_two_visc_afgdt_tauM_facMtau*
                                              (derxy_(0,ui)*viscs2_(0,1,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(1,1,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(2,1,vi)) ;
                elemat(vi*4 + 2, ui*4 + 3) += fac_two_visc_afgdt_tauM_facMtau*
                                              (derxy_(0,ui)*viscs2_(0,2,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(1,2,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(2,2,vi)) ;
              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)

            if (newton)
            {
              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui
                  =
                  fac_afgdt_two_visc_afgdt_tauM_facMtau*funct_(ui);

                for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
                {

                  /* viscous stabilisation --- convection
                     factor:
                                         alphaF*gamma*dt*tauM
            +(-)alphaF*gamma*dt*2*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                     /                                     \
                    |   /            \   n+af               |
                    |  | Dacc o nabla | u     , div eps (v) |
                    |   \            /                      |
                     \                                     /


                  */

                  elemat(vi*4     , ui*4    )+= fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(1, 0, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(2, 0, vi)*vderxyaf_(2, 0));
                  elemat(vi*4     , ui*4 + 1)+= fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(1, 0, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(2, 0, vi)*vderxyaf_(2, 1));
                  elemat(vi*4     , ui*4 + 2)+= fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(1, 0, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(2, 0, vi)*vderxyaf_(2, 2));
                  elemat(vi*4 + 1, ui*4     )+= fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(2, 1, vi)*vderxyaf_(2, 0));
                  elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(2, 1, vi)*vderxyaf_(2, 1)) ;
                  elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(2, 1, vi)*vderxyaf_(2, 2));
                  elemat(vi*4 + 2, ui*4    ) += fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 0)) ;
                  elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 1)) ;
                  elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_two_visc_afgdt_tauM_facMtau_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 2));
                } // end loop rows (test functions for matrix)
              } // end loop columns (solution for matrix, test function for vector)

            } // end if (a)gls and newton

#ifdef PERF
            timeelevstab_ref = null;
#endif
          } // end (a)gls stabilisation
        } // end higher_order_element

        if(cstab == Fluid3::continuity_stab_yes)
        {
          //---------------------------------------------------------------
          //
          //  STABILISATION PART, QUASISTATIC SUBGRID-PRESSURE VERSION
          //                  CONTINUITY STABILISATION
          //
          //---------------------------------------------------------------

#ifdef PERF
          RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif

          const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_gamma_dt_tauC_derxy_x_ui=fac_gamma_dt_tauC*derxy_(0,ui);
            const double fac_gamma_dt_tauC_derxy_y_ui=fac_gamma_dt_tauC*derxy_(1,ui);
            const double fac_gamma_dt_tauC_derxy_z_ui=fac_gamma_dt_tauC*derxy_(2,ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /*  factor: +gamma*dt*tauC

                    /                          \
                   |                            |
                   | nabla o Dacc  , nabla o v  |
                   |                            |
                    \                          /
              */

              elemat(vi*4    , ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(0,vi) ;
              elemat(vi*4    , ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi) ;
              elemat(vi*4    , ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(0,vi) ;
              elemat(vi*4 + 1, ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi) ;
              elemat(vi*4 + 1, ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(1,vi) ;
              elemat(vi*4 + 1, ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(1,vi) ;
              elemat(vi*4 + 2, ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(2,vi) ;
              elemat(vi*4 + 2, ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(2,vi) ;
              elemat(vi*4 + 2, ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(2,vi) ;
            } // end loop rows vi (test functions for matrix)
          } // end loop columns ui (solution for matrix, test function for vector)
#ifdef PERF
          timeelecstab_ref = null;
#endif
        } // end cstab
        else if(cstab == Fluid3::continuity_stab_td)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif
          //---------------------------------------------------------------
          //
          //  STABILISATION PART, TIME-DEPENDENT SUBGRID-PRESSURE VERSION
          //                  CONTINUITY STABILISATION
          //
          //---------------------------------------------------------------

          const double fac_gamma_dt_dt_factauC = fac*gamma*dt*dt*factauC;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_gamma_dt_dt_factauC_derxy_x_ui
              =
              fac_gamma_dt_dt_factauC*derxy_(0,ui);

            const double fac_gamma_dt_dt_factauC_derxy_y_ui
              =
              fac_gamma_dt_dt_factauC*derxy_(1,ui);

            const double fac_gamma_dt_dt_factauC_derxy_z_ui
              =
              fac_gamma_dt_dt_factauC*derxy_(2,ui);


              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /*
                                 tauC * dt
            factor: +gamma* dt * ---------
                                 tauC + dt
                    /                          \
                   |                            |
                   | nabla o Dacc  , nabla o v  |
                   |                            |
                    \                          /
              */

              elemat(vi*4    , ui*4    ) += fac_gamma_dt_dt_factauC_derxy_x_ui*derxy_(0,vi);
              elemat(vi*4    , ui*4 + 1) += fac_gamma_dt_dt_factauC_derxy_y_ui*derxy_(0,vi);
              elemat(vi*4    , ui*4 + 2) += fac_gamma_dt_dt_factauC_derxy_z_ui*derxy_(0,vi);
              elemat(vi*4 + 1, ui*4    ) += fac_gamma_dt_dt_factauC_derxy_x_ui*derxy_(1,vi);
              elemat(vi*4 + 1, ui*4 + 1) += fac_gamma_dt_dt_factauC_derxy_y_ui*derxy_(1,vi);
              elemat(vi*4 + 1, ui*4 + 2) += fac_gamma_dt_dt_factauC_derxy_z_ui*derxy_(1,vi);
              elemat(vi*4 + 2, ui*4    ) += fac_gamma_dt_dt_factauC_derxy_x_ui*derxy_(2,vi);
              elemat(vi*4 + 2, ui*4 + 1) += fac_gamma_dt_dt_factauC_derxy_y_ui*derxy_(2,vi);
              elemat(vi*4 + 2, ui*4 + 2) += fac_gamma_dt_dt_factauC_derxy_z_ui*derxy_(2,vi);

            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

#ifdef PERF
          timeelecstab_ref = null;
#endif
        }
        // end continuity stabilisation

        //---------------------------------------------------------------
        //
        //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
        //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
        //
        //---------------------------------------------------------------
        if(cross == Fluid3::cross_stress_stab)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif

          const double fac_afgdt=fac*afgdt;
          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_afgdt_conv_subaf_ui=fac_afgdt*conv_subaf_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {

              /*  factor:

               +alphaF*gamma*dt

                          /                          \
                         |  /~n+af       \            |
                         | | u    o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
              */
              const double fac_afgdt_conv_subaf_ui_funct_vi=fac_afgdt_conv_subaf_ui*funct_(vi);

              elemat(vi*4    , ui*4    ) += fac_afgdt_conv_subaf_ui_funct_vi;
              elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_conv_subaf_ui_funct_vi;
              elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_conv_subaf_ui_funct_vi;
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
#ifdef PERF
          timeelecrossrey_ref = null;
#endif
        } // end cross
      } // end if compute_elemat

      //---------------------------------------------------------------
      //---------------------------------------------------------------
      //
      //                       RIGHT HAND SIDE
      //
      //---------------------------------------------------------------
      //---------------------------------------------------------------
#ifdef PERF
      RefCountPtr<TimeMonitor> timeelegalerkin_ref = rcp(new TimeMonitor(*timeelegalerkin));
#endif
      //---------------------------------------------------------------
      //
      // (MODIFIED) GALERKIN PART, SUBSCALE ACCELERATION STABILISATION
      //
      //---------------------------------------------------------------
      if(inertia == Fluid3::inertia_stab_keep)
      {

        double aux_x =(-svelaf_(0)/tauM-pderxynp_(0)) ;
        double aux_y =(-svelaf_(1)/tauM-pderxynp_(1)) ;
        double aux_z =(-svelaf_(2)/tauM-pderxynp_(2)) ;

        if(higher_order_ele)
        {
          const double fact =2.0*visceff;

          aux_x += fact*viscaf_old_(0);
          aux_y += fact*viscaf_old_(1);
          aux_z += fact*viscaf_old_(2);
        }

        const double fac_sacc_plus_resM_not_partially_integrated_x =fac*aux_x ;
        const double fac_sacc_plus_resM_not_partially_integrated_y =fac*aux_y ;
        const double fac_sacc_plus_resM_not_partially_integrated_z =fac*aux_z ;


        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          //---------------------------------------------------------------
          //
          //     GALERKIN PART I AND SUBSCALE ACCELERATION STABILISATION
          //
          //---------------------------------------------------------------
          /*  factor: +1

               /             \     /
              |   ~ n+am      |   |     n+am    / n+af        \   n+af
              |  acc     , v  | + |  acc     + | c     o nabla | u     +
              |     (i)       |   |     (i)     \ (i)         /   (i)
               \             /     \

                                                   \
                                        n+af        |
                                     - f       , v  |
                                                    |
                                                   /

             using
                                                        /
                        ~ n+am        1.0      ~n+af   |    n+am
                       acc     = - --------- * u     - | acc     +
                          (i)           n+af    (i)    |    (i)
                                   tau_M                \

                                    / n+af        \   n+af            n+1
                                 + | c     o nabla | u     + nabla o p    -
                                    \ (i)         /   (i)             (i)

                                                            / n+af \
                                 - 2 * nu * grad o epsilon | u      | -
                                                            \ (i)  /
                                         \
                                    n+af  |
                                 - f      |
                                          |
                                         /

          */

          elevec[ui*4    ] -= fac_sacc_plus_resM_not_partially_integrated_x*funct_(ui) ;
          elevec[ui*4 + 1] -= fac_sacc_plus_resM_not_partially_integrated_y*funct_(ui) ;
          elevec[ui*4 + 2] -= fac_sacc_plus_resM_not_partially_integrated_z*funct_(ui) ;
        }

        //---------------------------------------------------------------
        //
        //   GALERKIN PART 2 (REMAINING EXPRESSIONS)
        //
        //---------------------------------------------------------------
        {
          const double fac_prenp_  = fac*prenp_ ;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            /* pressure */

            /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
            */

            elevec[ui*4    ] += fac_prenp_*derxy_(0,ui) ;
            elevec[ui*4 + 1] += fac_prenp_*derxy_(1,ui) ;
            elevec[ui*4 + 2] += fac_prenp_*derxy_(2,ui) ;
          }

          {
            const double fac_visceff = fac*visceff;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              /* viscous term */

            /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
            */

              elevec[ui*4    ] -= fac_visceff*
                                  (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                                   +
                                   derxy_(1,ui)*vderxyaf_(0,1)
                                   +
                                   derxy_(1,ui)*vderxyaf_(1,0)
                                   +
                                   derxy_(2,ui)*vderxyaf_(0,2)
                                   +
                                   derxy_(2,ui)*vderxyaf_(2,0)) ;
              elevec[ui*4 + 1] -= fac_visceff*
                                  (derxy_(0,ui)*vderxyaf_(0,1)
                                   +
                                   derxy_(0,ui)*vderxyaf_(1,0)
                                   +
                                   derxy_(1,ui)*vderxyaf_(1,1)*2.0
                                   +
                                   derxy_(2,ui)*vderxyaf_(1,2)
                                   +
                                   derxy_(2,ui)*vderxyaf_(2,1)) ;
              elevec[ui*4 + 2] -= fac_visceff*
                                  (derxy_(0,ui)*vderxyaf_(0,2)
                                   +
                                   derxy_(0,ui)*vderxyaf_(2,0)
                                   +
                                   derxy_(1,ui)*vderxyaf_(1,2)
                                   +
                                   derxy_(1,ui)*vderxyaf_(2,1)
                                   +
                                   derxy_(2,ui)*vderxyaf_(2,2)*2.0) ;

            }
          }
          {
            const double fac_divunp  = fac*divunp;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              /* continuity equation */

              /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
              */

              elevec[ui*4 + 3] -= fac_divunp*funct_(ui);

            } // end loop rows (solution for matrix, test function for vector)
          }
        }
      }
      else
      {
        //---------------------------------------------------------------
        //
        //        GALERKIN PART, NEGLECTING SUBSCALE ACCLERATIONS
        //
        //---------------------------------------------------------------
        {
          const double fac_inertia_convection_dead_load_x
            =
            fac*(accintam_(0)+convaf_old_(0)-bodyforceaf_(0));

          const double fac_inertia_convection_dead_load_y
            =
            fac*(accintam_(1)+convaf_old_(1)-bodyforceaf_(1));

          const double fac_inertia_convection_dead_load_z
            =
            fac*(accintam_(2)+convaf_old_(2)-bodyforceaf_(2));

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            /* inertia terms */

            /*  factor: +1

               /             \
              |     n+am      |
              |  acc     , v  |
              |               |
               \             /
            */

            /* convection */

            /*  factor: +1

               /                             \
              |  / n+af       \    n+af       |
              | | c    o nabla |  u      , v  |
              |  \            /               |
               \                             /
            */

            /* body force (dead load...) */

            /*  factor: -1

               /           \
              |   n+af      |
              |  f     , v  |
              |             |
               \           /
            */

            elevec[ui*4    ] -= funct_(ui)*fac_inertia_convection_dead_load_x;
            elevec[ui*4 + 1] -= funct_(ui)*fac_inertia_convection_dead_load_y;
            elevec[ui*4 + 2] -= funct_(ui)*fac_inertia_convection_dead_load_z;
          }
        }

        {
          const double fac_prenp=fac*prenp_;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {

            /* pressure */

            /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
            */

            elevec[ui*4    ] += fac_prenp*derxy_(0,ui);
            elevec[ui*4 + 1] += fac_prenp*derxy_(1,ui);
            elevec[ui*4 + 2] += fac_prenp*derxy_(2,ui);
          }
        }

        {
          const double visceff_fac=visceff*fac;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {

            /* viscous term */

            /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
            */

            elevec[ui*4    ] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                                 +
                                 derxy_(1,ui)*vderxyaf_(0,1)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,0)
                                 +
                                 derxy_(2,ui)*vderxyaf_(0,2)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,0)) ;
            elevec[ui*4 + 1] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,1)
                                 +
                                 derxy_(0,ui)*vderxyaf_(1,0)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,1)*2.0
                                 +
                                 derxy_(2,ui)*vderxyaf_(1,2)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,1)) ;
            elevec[ui*4 + 2] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,2)
                                 +
                                 derxy_(0,ui)*vderxyaf_(2,0)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,2)
                                 +
                                 derxy_(1,ui)*vderxyaf_(2,1)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,2)*2.0) ;
          }
        }

        {
          const double fac_divunp=fac*divunp;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {

            /* continuity equation */

            /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
            */

            elevec[ui*4 + 3] -= fac_divunp*funct_(ui);
          } // end loop rows (solution for matrix, test function for vector)
        }
      }
#ifdef PERF
      timeelegalerkin_ref = null;
#endif

      //---------------------------------------------------------------
      //
      //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //                    PRESSURE STABILISATION
      //
      //---------------------------------------------------------------
      if(pspg == Fluid3::pstab_use_pspg)
      {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelepspg_ref = rcp(new TimeMonitor(*timeelepspg));
#endif

        const double fac_svelnpx                      = fac*svelnp(0,iquad);
        const double fac_svelnpy                      = fac*svelnp(1,iquad);
        const double fac_svelnpz                      = fac*svelnp(2,iquad);

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          /* factor: -1

                       /                 \
                      |  ~n+1             |
                      |  u    , nabla  q  |
                      |   (i)             |
                       \                 /
          */

          elevec[ui*4 + 3] += fac_svelnpx*derxy_(0,ui)+fac_svelnpy*derxy_(1,ui)+fac_svelnpz*derxy_(2,ui);

        } // end loop rows (solution for matrix, test function for vector)

#ifdef PERF
          timeelepspg_ref = null;
#endif
      }

      //---------------------------------------------------------------
      //
      //         STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
      //
      //---------------------------------------------------------------
      if(supg == Fluid3::convective_stab_supg)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelesupg_ref = rcp(new TimeMonitor(*timeelesupg));
#endif

        const double fac_svelaf_x=fac*svelaf_(0);
        const double fac_svelaf_y=fac*svelaf_(1);
        const double fac_svelaf_z=fac*svelaf_(2);

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          /*
                  /                             \
                 |  ~n+af    / n+af        \     |
                 |  u     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /

          */

          elevec[ui*4    ] += fac_svelaf_x*conv_c_af_(ui);
          elevec[ui*4 + 1] += fac_svelaf_y*conv_c_af_(ui);
          elevec[ui*4 + 2] += fac_svelaf_z*conv_c_af_(ui);

        } // end loop rows (solution for matrix, test function for vector)
#ifdef PERF
        timeelesupg_ref = null;
#endif

      }

      //---------------------------------------------------------------
      //
      //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //             VISCOUS STABILISATION (FOR (A)GLS)
      //
      //---------------------------------------------------------------
      if (higher_order_ele)
      {
        if (vstab != Fluid3::viscous_stab_none)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelevstab_ref = rcp(new TimeMonitor(*timeelevstab));
#endif

          const double fac_two_visc_svelaf_x = vstabfac*fac*2.0*visc*svelaf_(0);
          const double fac_two_visc_svelaf_y = vstabfac*fac*2.0*visc*svelaf_(1);
          const double fac_two_visc_svelaf_z = vstabfac*fac*2.0*visc*svelaf_(2);

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {

            /*
                 /                      \
                |  ~n+af                 |
                |  u      , div eps (v)  |
                |                        |
                 \                      /

            */
            elevec[ui*4    ] += fac_two_visc_svelaf_x*viscs2_(0, 0, ui)
                                +
                                fac_two_visc_svelaf_y*viscs2_(1, 0, ui)
                                +
                                fac_two_visc_svelaf_z*viscs2_(2, 0, ui);

            elevec[ui*4 + 1] += fac_two_visc_svelaf_x*viscs2_(0, 1, ui)
                                +
                                fac_two_visc_svelaf_y*viscs2_(1, 1, ui)
                                +
                                fac_two_visc_svelaf_z*viscs2_(2, 1, ui);

            elevec[ui*4 + 2] += fac_two_visc_svelaf_x*viscs2_(0, 2, ui)
                                +
                                fac_two_visc_svelaf_y*viscs2_(1, 2, ui)
                                +
                                fac_two_visc_svelaf_z*viscs2_(2, 2, ui);

          } // end loop rows (solution for matrix, test function for vector)
#ifdef PERF
          timeelevstab_ref = null;
#endif
        } // endif (a)gls
      }// end if higher order ele

      if(cstab == Fluid3::continuity_stab_yes)
      {
        //---------------------------------------------------------------
        //
        //            QUASISTATIC VERSION OF RESIDUAL BASED
        //                  CONTINUITY STABILISATION
        //          (the original version proposed by Codina)
        //
        //---------------------------------------------------------------
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif
        const double fac_tauC_divunp = fac*tauC*divunp;
        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
          /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
          */

          elevec[ui*4    ] -= fac_tauC_divunp*derxy_(0,ui) ;
          elevec[ui*4 + 1] -= fac_tauC_divunp*derxy_(1,ui) ;
          elevec[ui*4 + 2] -= fac_tauC_divunp*derxy_(2,ui) ;
        } // end loop rows
#ifdef PERF
        timeelecstab_ref = null;
#endif
      }
      else if (cstab == Fluid3::continuity_stab_td)
      {
        //---------------------------------------------------------------
        //
        //           TIME DEPENDENT VERSION OF RESIDUAL BASED
        //                  CONTINUITY STABILISATION
        //
        //---------------------------------------------------------------

#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif
        const double fac_sprenp = fac*sprenp(iquad);

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {


          /* factor: -1

                       /                  \
                      |  ~n+1              |
                      |  p    , nabla o v  |
                      |   (i)              |
                       \                  /
          */
          elevec[ui*4    ] += fac_sprenp*derxy_(0,ui) ;
          elevec[ui*4 + 1] += fac_sprenp*derxy_(1,ui) ;
          elevec[ui*4 + 2] += fac_sprenp*derxy_(2,ui) ;
        } // end loop rows (solution for matrix, test function for vector)
#ifdef PERF
        timeelecstab_ref = null;
#endif
      }


      //---------------------------------------------------------------
      //
      //        TIME-DEPENDENT SUBGRID-SCALE STABILISATION
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif
        const double fac_convsubaf_old_x=fac*convsubaf_old_(0);
        const double fac_convsubaf_old_y=fac*convsubaf_old_(1);
        const double fac_convsubaf_old_z=fac*convsubaf_old_(2);

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          /* factor:

                  /                           \
                 |   ~n+af           n+af      |
                 | ( u    o nabla ) u     , v  |
                 |    (i)            (i)       |
                  \                           /
          */
          elevec[ui*4    ] -= fac_convsubaf_old_x*funct_(ui);
          elevec[ui*4 + 1] -= fac_convsubaf_old_y*funct_(ui);
          elevec[ui*4 + 2] -= fac_convsubaf_old_z*funct_(ui);
        }
#ifdef PERF
          timeelecrossrey_ref = null;
#endif
      }

      //---------------------------------------------------------------
      //
      //       TIME DEPENDENT SUBGRID-SCALE STABILISATION PART
      //     RESIDUAL BASED VMM STABILISATION --- REYNOLDS STRESS
      //
      //---------------------------------------------------------------
      if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif
        const double fac_svelaf_x=fac*svelaf_(0);
        const double fac_svelaf_y=fac*svelaf_(1);
        const double fac_svelaf_z=fac*svelaf_(2);

        for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
        {

          /* factor:

                  /                             \
                 |  ~n+af      ~n+af             |
                 |  u      , ( u    o nabla ) v  |
                 |                               |
                  \                             /
          */
          elevec[ui*4    ] += fac_svelaf_x*(svelaf_(0)*derxy_(0,ui)
                                            +
                                            svelaf_(1)*derxy_(1,ui)
                                            +
                                            svelaf_(2)*derxy_(2,ui));
          elevec[ui*4 + 1] += fac_svelaf_y*(svelaf_(0)*derxy_(0,ui)
                                            +
                                            svelaf_(1)*derxy_(1,ui)
                                            +
                                            svelaf_(2)*derxy_(2,ui));
          elevec[ui*4 + 2] += fac_svelaf_z*(svelaf_(0)*derxy_(0,ui)
                                            +
                                            svelaf_(1)*derxy_(1,ui)
                                            +
                                            svelaf_(2)*derxy_(2,ui));

        } // end loop rows (solution for matrix, test function for vector)
#ifdef PERF
      timeelecrossrey_ref = null;
#endif
      }
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//     ELEMENT FORMULATION BASED ON QUASISTATIC SUBSCALES
//
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
    else
    {
      // stabilisation parameters
      const double tauM   = tau_(0);
      const double tauMp  = tau_(1);
      const double tauC   = tau_(2);

      // subgrid-viscosity factor
      const double vartfac = vart_*fac;

      //--------------------------------------------------------------
      //--------------------------------------------------------------
      //
      //              SYSTEM MATRIX, QUASISTATIC FORMULATION
      //
      //--------------------------------------------------------------
      //--------------------------------------------------------------
      if(compute_elemat)
      {

        //---------------------------------------------------------------
        //
        //             GALERKIN PART, QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelegalerkin_ref = rcp(new TimeMonitor(*timeelegalerkin));
#endif

          const double fac_visceff_afgdt = fac*visceff*afgdt;
          const double fac_gamma_dt      = fac*gamma*dt;

          const double fac_alphaM        = fac*alphaM;
          const double fac_afgdt         = fac*afgdt;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double inertia_and_conv_ui
              = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double inertia_and_conv_ui_funct_vi
                = inertia_and_conv_ui*funct_(vi);


              /*
                inertia term (intermediate) + convection (intermediate)

                factor:

                          /          \                   /                          \
                         |            |                 |  / n+af       \            |
                +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                         |            |                 |  \            /            |
                          \          /                   \                          /
              */
              elemat(vi*4    , ui*4    ) += inertia_and_conv_ui_funct_vi;
              elemat(vi*4 + 1, ui*4 + 1) += inertia_and_conv_ui_funct_vi;
              elemat(vi*4 + 2, ui*4 + 2) += inertia_and_conv_ui_funct_vi;

            }
          }


          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_funct_ui = fac*funct_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /* pressure (implicit) */

              /*  factor: -1

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
              */

              elemat(vi*4    , ui*4 + 3) -= fac_funct_ui*derxy_(0, vi) ;
              elemat(vi*4 + 1, ui*4 + 3) -= fac_funct_ui*derxy_(1, vi) ;
              elemat(vi*4 + 2, ui*4 + 3) -= fac_funct_ui*derxy_(2, vi) ;
            }
          }

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double sum = derxy_(0,ui)*derxy_(0,vi)+derxy_(1,ui)*derxy_(1,vi)+derxy_(2,ui)*derxy_(2,vi);

              /* viscous term (intermediate) */

              /*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
              */

              elemat(vi*4    , ui*4    ) += fac_visceff_afgdt*(sum+derxy_(0,ui)*derxy_(0,vi)) ;
              elemat(vi*4    , ui*4 + 1) += fac_visceff_afgdt*     derxy_(0,ui)*derxy_(1,vi) ;
              elemat(vi*4    , ui*4 + 2) += fac_visceff_afgdt*     derxy_(0,ui)*derxy_(2,vi) ;
              elemat(vi*4 + 1, ui*4    ) += fac_visceff_afgdt*     derxy_(1,ui)*derxy_(0,vi) ;
              elemat(vi*4 + 1, ui*4 + 1) += fac_visceff_afgdt*(sum+derxy_(1,ui)*derxy_(1,vi)) ;
              elemat(vi*4 + 1, ui*4 + 2) += fac_visceff_afgdt*     derxy_(1,ui)*derxy_(2,vi) ;
              elemat(vi*4 + 2, ui*4    ) += fac_visceff_afgdt*     derxy_(2,ui)*derxy_(0,vi) ;
              elemat(vi*4 + 2, ui*4 + 1) += fac_visceff_afgdt*     derxy_(2,ui)*derxy_(1,vi) ;
              elemat(vi*4 + 2, ui*4 + 2) += fac_visceff_afgdt*(sum+derxy_(2,ui)*derxy_(2,vi)) ;
            }
          }

          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {
            const double fac_gamma_dt_funct_vi =fac_gamma_dt*funct_(vi);

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {

              /* continuity equation (implicit) */

              /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
              */

              elemat(vi*4 + 3, ui*4    ) += fac_gamma_dt_funct_vi*derxy_(0,ui) ;
              elemat(vi*4 + 3, ui*4 + 1) += fac_gamma_dt_funct_vi*derxy_(1,ui) ;
              elemat(vi*4 + 3, ui*4 + 2) += fac_gamma_dt_funct_vi*derxy_(2,ui) ;

            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

          if (newton)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {

              const double fac_afgdt_funct_vi = fac_afgdt*funct_(vi);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double  fac_afgdt_funct_vi_funct_ui=fac_afgdt_funct_vi*funct_(ui);

                /* convection (intermediate)

                 factor:

                        +alphaF*gamma*dt


                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
                */

                elemat(vi*4    , ui*4    ) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(0,0);
                elemat(vi*4    , ui*4 + 1) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(0,1);
                elemat(vi*4    , ui*4 + 2) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(0,2);
                elemat(vi*4 + 1, ui*4    ) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(1,0);
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(1,1);
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(1,2);
                elemat(vi*4 + 2, ui*4    ) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(2,0);
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(2,1);
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_funct_vi_funct_ui*vderxyaf_(2,2);
              } // end loop rows (test functions for matrix)
            } // end loop rows (solution for matrix, test function for vector)
          }

#ifdef PERF
          timeelegalerkin_ref = null;
#endif
        }

        //---------------------------------------------------------------
        //
        //             PSPG PART, QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        if(pspg == Fluid3::pstab_use_pspg)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelepspg_ref = rcp(new TimeMonitor(*timeelepspg));
#endif
          const double fac_afgdt_tauMp             = fac*afgdt*tauMp;
          const double fac_tauMp                   = fac*tauMp;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {


              /* pressure stabilisation --- inertia+convection    */

              /* factor:

                             /                \
                            |                  |
              +alphaM*tauMp*|  Dacc , nabla q  |+
                            |                  |
                             \                /
                                          /                                \
                                         |  / n+af       \                  |
                  +alphaF*gamma*dt*tauMp*| | c    o nabla | Dacc , nabla q  |
                                         |  \            /                  |
                                          \                                /
              */

              elemat(vi*4 + 3, ui*4    ) += fac_tauMp_inertia_and_conv*derxy_(0,vi) ;
              elemat(vi*4 + 3, ui*4 + 1) += fac_tauMp_inertia_and_conv*derxy_(1,vi) ;
              elemat(vi*4 + 3, ui*4 + 2) += fac_tauMp_inertia_and_conv*derxy_(2,vi) ;
            }
          }


          if (higher_order_ele)
          {
            const double fac_two_visceff_afgdt_tauMp
              =
              fac*2.0*visceff*afgdt*tauMp;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                /* pressure stabilisation --- diffusion  */

                /* factor: -2*nu*alphaF*gamma*dt*tauMp

                    /                                \
                   |               /    \             |
                   |  nabla o eps | Dacc | , nabla q  |
                   |               \    /             |
                    \                                /
                */

                elemat(vi*4 + 3, ui*4    ) -= fac_two_visceff_afgdt_tauMp*
                                              (derxy_(0,vi)*viscs2_(0,0,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(0,1,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(0,2,ui)) ;
                elemat(vi*4 + 3, ui*4 + 1) -= fac_two_visceff_afgdt_tauMp*
                                              (derxy_(0,vi)*viscs2_(0,1,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(1,1,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(1,2,ui)) ;
                elemat(vi*4 + 3, ui*4 + 2) -= fac_two_visceff_afgdt_tauMp*
                                              (derxy_(0,vi)*viscs2_(0,2,ui)
                                               +
                                               derxy_(1,vi)*viscs2_(1,2,ui)
                                               +
                                               derxy_(2,vi)*viscs2_(2,2,ui)) ;
              }
            }
          }


          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /* pressure stabilisation --- pressure   */

              /* factor: +tauMp

                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
              */

              elemat(vi*4 + 3, ui*4 + 3) += fac_tauMp*
                                            (derxy_(0,ui)*derxy_(0,vi)
                                             +
                                             derxy_(1,ui)*derxy_(1,vi)
                                             +
                                             derxy_(2,ui)*derxy_(2,vi)) ;

            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

          if (newton)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              int vidx = vi*4 + 3;
              double v1 = derxy_(0,vi)*vderxyaf_(0,0) + derxy_(1,vi)*vderxyaf_(1,0) + derxy_(2,vi)*vderxyaf_(2,0);
              double v2 = derxy_(0,vi)*vderxyaf_(0,1) + derxy_(1,vi)*vderxyaf_(1,1) + derxy_(2,vi)*vderxyaf_(2,1);
              double v3 = derxy_(0,vi)*vderxyaf_(0,2) + derxy_(1,vi)*vderxyaf_(1,2) + derxy_(2,vi)*vderxyaf_(2,2);
              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double fac_afgdt_tauMp_funct_ui = fac_afgdt_tauMp*funct_(ui);
                int uidx = ui*4;

                /* pressure stabilisation --- convection */

                /*  factor: +alphaF*gamma*dt*tauMp

                       /                                  \
                      |  /            \   n+af             |
                      | | Dacc o nabla | u      , nabla q  |
                      |  \            /                    |
                       \                                  /
                */

                elemat(vidx, uidx    ) += fac_afgdt_tauMp_funct_ui*v1;
                elemat(vidx, uidx + 1) += fac_afgdt_tauMp_funct_ui*v2;
                elemat(vidx, uidx + 2) += fac_afgdt_tauMp_funct_ui*v3;
              } // end loop rows (test functions for matrix)
            } // end loop rows (solution for matrix, test function for vector)
          }
#ifdef PERF
          timeelepspg_ref = null;
#endif
        }

        //---------------------------------------------------------------
        //
        //             SUPG PART, QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        if(supg == Fluid3::convective_stab_supg)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelesupg_ref = rcp(new TimeMonitor(*timeelesupg));
#endif
          const double fac_afgdt_tauM             = fac*tauM*afgdt;
          const double fac_two_visceff_afgdt_tauM = fac*tauM*afgdt*2.0*visceff;
          const double fac_tauM                   = fac*tauM;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_tauM_inertia_and_conv
              = fac_tauM*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_tauM_inertia_and_conv_conv_c_af_vi
                =
                fac_tauM_inertia_and_conv*conv_c_af_(vi);

              /* SUPG stabilisation --- inertia and convection */

              /* factor: +alphaM*tauM

                           /                           \
                          |          / n+af       \     |
             +alphaM*tauM*|  Dacc , | c    o nabla | v  |+
                          |          \            /     |
                           \                           /


                                     /                                               \
                                    |    / n+af        \          / n+af        \     |
              +alphaF*gamma*dt*tauM*|   | c     o nabla | Dacc , | c     o nabla | v  |
                                    |    \             /          \             /     |
                                     \                                               /

              */

              elemat(vi*4    , ui*4    ) += fac_tauM_inertia_and_conv_conv_c_af_vi;
              elemat(vi*4 + 1, ui*4 + 1) += fac_tauM_inertia_and_conv_conv_c_af_vi;
              elemat(vi*4 + 2, ui*4 + 2) += fac_tauM_inertia_and_conv_conv_c_af_vi;
            }
          }

          if (higher_order_ele)
          {
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_two_visceff_afgdt_tauM_conv_c_af_vi = fac_two_visceff_afgdt_tauM*conv_c_af_(vi);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                /* SUPG stabilisation --- diffusion   */

                /* factor: -2*nu*alphaF*gamma*dt*tauM


                    /                                            \
                   |               /     \    / n+af        \     |
                   |  nabla o eps | Dacc  |, | c     o nabla | v  |
                   |               \     /    \             /     |
                    \                                            /

                */

                elemat(vi*4, ui*4)         -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(0, 0, ui) ;
                elemat(vi*4, ui*4 + 1)     -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(0, 1, ui) ;
                elemat(vi*4, ui*4 + 2)     -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(0, 2, ui) ;
                elemat(vi*4 + 1, ui*4)     -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(0, 1, ui) ;
                elemat(vi*4 + 1, ui*4 + 1) -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(1, 1, ui) ;
                elemat(vi*4 + 1, ui*4 + 2) -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(1, 2, ui) ;
                elemat(vi*4 + 2, ui*4)     -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(0, 2, ui) ;
                elemat(vi*4 + 2, ui*4 + 1) -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(1, 2, ui) ;
                elemat(vi*4 + 2, ui*4 + 2) -= fac_two_visceff_afgdt_tauM_conv_c_af_vi*viscs2_(2, 2, ui) ;
              } //end ui
            } // end vi
          }// end higher_order_ele

          for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
          {
            const double fac_tauM_conv_c_af_vi = fac_tauM*conv_c_af_(vi);

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              /* SUPG stabilisation --- pressure    */

              /* factor: +tauM

                    /                               \
                   |              / n+af       \     |
                   |  nabla Dp , | c    o nabla | v  |
                   |              \            /     |
                    \                               /
              */

              elemat(vi*4    , ui*4 + 3) += fac_tauM_conv_c_af_vi*derxy_(0,ui) ;
              elemat(vi*4 + 1, ui*4 + 3) += fac_tauM_conv_c_af_vi*derxy_(1,ui) ;
              elemat(vi*4 + 2, ui*4 + 3) += fac_tauM_conv_c_af_vi*derxy_(2,ui) ;

            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

          if (newton)
          {
            {
              const double fac_afgdt_tauM_resM_x=fac_afgdt_tauM*resM_(0);
              const double fac_afgdt_tauM_resM_y=fac_afgdt_tauM*resM_(1);
              const double fac_afgdt_tauM_resM_z=fac_afgdt_tauM*resM_(2);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double fac_afgdt_tauM_resM_x_funct_ui = fac_afgdt_tauM_resM_x*funct_(ui);
                const double fac_afgdt_tauM_resM_y_funct_ui = fac_afgdt_tauM_resM_y*funct_(ui);
                const double fac_afgdt_tauM_resM_z_funct_ui = fac_afgdt_tauM_resM_z*funct_(ui);

                for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
                {
                  /* SUPG stabilisation:

                                         /                               \
                                        |     n+am     /            \     |
                  +alphaF*gamma*dt*tauM*|  acc      , | Dacc o nabla | v  | -
                                        |              \            /     |
                                         \                               /

                                         /                                               \
                                        |    / n+af        \   n+af    /            \     |
                  +alphaF*gamma*dt*tauM*|   | c     o nabla | u     , | Dacc o nabla | v  |
                                        |    \             /           \            /     |
                                         \                                               /

                                         /                                 \
                                        |         n+1    /            \     |
                  +alphaF*gamma*dt*tauM*|  nabla p    , | Dacc o nabla | v  |
                                        |                \            /     |
                                         \                                 /

                                              /                                            \
                                             |               / n+af \    /            \     |
                  -2*nu*alphaF*gamma*dt*tauM*|  nabla o eps | u      |, | Dacc o nabla | v  |
                                             |               \      /    \            /     |
                                              \                                            /

                                         /                            \
                                        |   n+af    /            \     |
                  -tauM*alphaF*gamma*dt*|  f     , | Dacc o nabla | v  |
                                        |           \            /     |
                                         \                            /

                  */

                  elemat(vi*4    , ui*4    ) += fac_afgdt_tauM_resM_x_funct_ui*derxy_(0,vi) ;
                  elemat(vi*4    , ui*4 + 1) += fac_afgdt_tauM_resM_x_funct_ui*derxy_(1,vi) ;
                  elemat(vi*4    , ui*4 + 2) += fac_afgdt_tauM_resM_x_funct_ui*derxy_(2,vi) ;
                  elemat(vi*4 + 1, ui*4    ) += fac_afgdt_tauM_resM_y_funct_ui*derxy_(0,vi) ;
                  elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_tauM_resM_y_funct_ui*derxy_(1,vi) ;
                  elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_tauM_resM_y_funct_ui*derxy_(2,vi) ;
                  elemat(vi*4 + 2, ui*4    ) += fac_afgdt_tauM_resM_z_funct_ui*derxy_(0,vi) ;
                  elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_tauM_resM_z_funct_ui*derxy_(1,vi) ;
                  elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_tauM_resM_z_funct_ui*derxy_(2,vi) ;
                }
              }
            }
            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_afgdt_tauM_conv_c_af_vi =fac_afgdt_tauM*conv_c_af_(vi);

              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double fac_afgdt_tauM_conv_c_af_vi_funct_ui =fac_afgdt_tauM_conv_c_af_vi*funct_(ui);

                /* SUPG stabilisation --- convection  */

                /* factor: +alphaF*gamma*dt*tauM

                    /                                               \
                   |    /            \   n+af    / n+af        \     |
                   |   | Dacc o nabla | u     , | c     o nabla | v  |
                   |    \            /           \             /     |
                    \                                               /
                */

                elemat(vi*4    , ui*4    ) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(0,0);
                elemat(vi*4    , ui*4 + 1) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(0,1);
                elemat(vi*4    , ui*4 + 2) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(0,2);
                elemat(vi*4 + 1, ui*4    ) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(1,0);
                elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(1,1);
                elemat(vi*4 + 1, ui*4 + 2) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(1,2);
                elemat(vi*4 + 2, ui*4    ) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(2,0);
                elemat(vi*4 + 2, ui*4 + 1) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(2,1);
                elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_tauM_conv_c_af_vi_funct_ui*vderxyaf_(2,2);
              }
            }
          }
#ifdef PERF
          timeelesupg_ref = null;
#endif
        }

        //---------------------------------------------------------------
        //
        //      VISCOUS STABILISATION PART, QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        if (higher_order_ele)
        {
          if((vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_usfem)&&higher_order_ele)
          {
#ifdef PERF
            RefCountPtr<TimeMonitor> timeelevstab_ref = rcp(new TimeMonitor(*timeelevstab));
#endif
            const double fac_two_visc_tauMp             = vstabfac*fac*2.0*visc*tauMp;
            const double fac_two_visc_afgdt_tauMp       = vstabfac*fac*2.0*visc*afgdt*tauMp;
            const double fac_two_visc_alphaM_tauMp      = vstabfac*fac*2.0*visc*alphaM*tauMp;
            const double fac_four_visceff_visc_afgdt_tauMp = vstabfac*fac*4.0*visceff*visc*afgdt*tauMp;

            for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
              {
                /* viscous stabilisation --- inertia     */

                /* factor: +(-)alphaM*tauMp*2*nu

                    /                    \
                   |                      |
                   |  Dacc , div eps (v)  |
                   |                      |
                    \                    /
                */
                elemat(vi*4    , ui*4    ) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(0,0,vi);
                elemat(vi*4    , ui*4 + 1) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(0,1,vi);
                elemat(vi*4    , ui*4 + 2) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(0,2,vi);
                elemat(vi*4 + 1, ui*4    ) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(0,1,vi);
                elemat(vi*4 + 1, ui*4 + 1) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(1,1,vi);
                elemat(vi*4 + 1, ui*4 + 2) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(1,2,vi);
                elemat(vi*4 + 2, ui*4    ) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(0,2,vi);
                elemat(vi*4 + 2, ui*4 + 1) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(1,2,vi);
                elemat(vi*4 + 2, ui*4 + 2) += fac_two_visc_alphaM_tauMp*funct_(ui)*viscs2_(2,2,vi);


                /* viscous stabilisation --- convection */

                /*  factor: +(-)2*nu*alphaF*gamma*dt*tauMp

                       /                                  \
                      |  / n+af       \                    |
                      | | c    o nabla | Dacc, div eps (v) |
                      |  \            /                    |
                       \                                  /

                */
                elemat(vi*4     , ui*4    )+= fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(0, 0, vi) ;
                elemat(vi*4     , ui*4 + 1)+= fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(0, 1, vi) ;
                elemat(vi*4     , ui*4 + 2)+= fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(0, 2, vi) ;
                elemat(vi*4 + 1, ui*4     )+= fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(0, 1, vi) ;
                elemat(vi*4 + 1, ui*4 + 1) += fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(1, 1, vi) ;
                elemat(vi*4 + 1, ui*4 + 2) += fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(1, 2, vi) ;
                elemat(vi*4 + 2, ui*4    ) += fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(0, 2, vi) ;
                elemat(vi*4 + 2, ui*4 + 1) += fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(1, 2, vi) ;
                elemat(vi*4 + 2, ui*4 + 2) += fac_two_visc_afgdt_tauMp*conv_c_af_(ui)*viscs2_(2, 2, vi) ;

                /* viscous stabilisation --- diffusion  */

                /* factor: -(+)4*nu*nu*alphaF*gamma*dt*tauMp

                    /                                   \
                   |               /    \                |
                   |  nabla o eps | Dacc | , div eps (v) |
                   |               \    /                |
                    \                                   /
                */

                elemat(vi*4    , ui*4    ) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 0, vi)
                                               +
                                               viscs2_(0, 1, ui)*viscs2_(0, 1, vi)
                                               +
                                               viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
                elemat(vi*4    , ui*4 + 1) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 0, vi)*viscs2_(0, 1, ui)
                                               +
                                               viscs2_(0, 1, vi)*viscs2_(1, 1, ui)
                                               +
                                               viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
                elemat(vi*4    , ui*4 + 2) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 0, vi)*viscs2_(0, 2, ui)
                                               +
                                               viscs2_(0, 1, vi)*viscs2_(1, 2, ui)
                                               +
                                               viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
                elemat(vi*4 + 1, ui*4    ) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 1, vi)
                                               +
                                               viscs2_(0, 1, ui)*viscs2_(1, 1, vi)
                                               +
                                               viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
                elemat(vi*4 + 1, ui*4 + 1) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 1, ui)*viscs2_(0, 1, vi)
                                               +
                                               viscs2_(1, 1, ui)*viscs2_(1, 1, vi)
                                               +
                                               viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
                elemat(vi*4 + 1, ui*4 + 2) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 1, vi)*viscs2_(0, 2, ui)
                                               +
                                               viscs2_(1, 1, vi)*viscs2_(1, 2, ui)
                                               +
                                               viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
                elemat(vi*4 + 2, ui*4    ) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 0, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(0, 1, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
                elemat(vi*4 + 2, ui*4 + 1) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 1, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(1, 1, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
                elemat(vi*4 + 2, ui*4 + 2) -= fac_four_visceff_visc_afgdt_tauMp*
                                              (viscs2_(0, 2, ui)*viscs2_(0, 2, vi)
                                               +
                                               viscs2_(1, 2, ui)*viscs2_(1, 2, vi)
                                               +
                                               viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

                /* viscous stabilisation --- pressure   */

                /* factor: +(-)tauMp*2*nu

                    /                        \
                   |                          |
                   |  nabla Dp , div eps (v)  |
                   |                          |
                    \                        /
                */
                elemat(vi*4    , ui*4 + 3) += fac_two_visc_tauMp*
                                              (derxy_(0,ui)*viscs2_(0,0,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(0,1,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(0,2,vi)) ;
                elemat(vi*4 + 1, ui*4 + 3) += fac_two_visc_tauMp*
                                              (derxy_(0,ui)*viscs2_(0,1,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(1,1,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(1,2,vi)) ;
                elemat(vi*4 + 2, ui*4 + 3) += fac_two_visc_tauMp*
                                              (derxy_(0,ui)*viscs2_(0,2,vi)
                                               +
                                               derxy_(1,ui)*viscs2_(1,2,vi)
                                               +
                                               derxy_(2,ui)*viscs2_(2,2,vi)) ;

              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)

            if (newton)
            {
              for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
              {
                const double fac_two_visc_afgdt_tauMp_funct_ui=fac_two_visc_afgdt_tauMp*funct_(ui);

                for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
                {
                  /* viscous stabilisation --- convection */

                  /*  factor: +(-)2*nu*alphaF*gamma*dt*tauMp

                     /                                     \
                    |   /            \   n+af               |
                    |  | Dacc o nabla | u     , div eps (v) |
                    |   \            /                      |
                     \                                     /


                  */
                  elemat(vi*4     , ui*4    )+= fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(0, 1, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(0, 2, vi)*vderxyaf_(2, 0)) ;
                  elemat(vi*4     , ui*4 + 1)+= fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(0, 1, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(0, 2, vi)*vderxyaf_(2, 1)) ;
                  elemat(vi*4     , ui*4 + 2)+= fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 0, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(0, 1, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(0, 2, vi)*vderxyaf_(2, 2)) ;
                  elemat(vi*4 + 1, ui*4     )+= fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(2, 0)) ;
                  elemat(vi*4 + 1, ui*4 + 1) += fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(2, 1)) ;
                  elemat(vi*4 + 1, ui*4 + 2) += fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 1, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(1, 1, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(2, 2)) ;
                  elemat(vi*4 + 2, ui*4    ) += fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 0)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 0)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 0)) ;
                  elemat(vi*4 + 2, ui*4 + 1) += fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 1)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 1)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 1)) ;
                  elemat(vi*4 + 2, ui*4 + 2) += fac_two_visc_afgdt_tauMp_funct_ui*
                                                (viscs2_(0, 2, vi)*vderxyaf_(0, 2)
                                                 +
                                                 viscs2_(1, 2, vi)*vderxyaf_(1, 2)
                                                 +
                                                 viscs2_(2, 2, vi)*vderxyaf_(2, 2)) ;
                } // end loop rows (test functions for matrix)
              } // end loop columns (solution for matrix, test function for vector)
            }
#ifdef PERF
            timeelevstab_ref = null;
#endif
          } // endif (a)gls
        }

        //---------------------------------------------------------------
        //
        //         LEAST SQUARES CONTINUITY STABILISATION PART,
        //                  QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        if(cstab == Fluid3::continuity_stab_yes)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif

          const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
            const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);
            const double fac_gamma_dt_tauC_derxy_z_ui = fac_gamma_dt_tauC*derxy_(2,ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              /*  factor: +gamma*dt*tauC

                    /                          \
                   |                            |
                   | nabla o Dacc  , nabla o v  |
                   |                            |
                    \                          /
              */

              elemat(vi*4    , ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(0,vi) ;
              elemat(vi*4    , ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi) ;
              elemat(vi*4    , ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(0,vi) ;
              elemat(vi*4 + 1, ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi) ;
              elemat(vi*4 + 1, ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(1,vi) ;
              elemat(vi*4 + 1, ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(1,vi) ;
              elemat(vi*4 + 2, ui*4    ) += fac_gamma_dt_tauC_derxy_x_ui*derxy_(2,vi) ;
              elemat(vi*4 + 2, ui*4 + 1) += fac_gamma_dt_tauC_derxy_y_ui*derxy_(2,vi) ;
              elemat(vi*4 + 2, ui*4 + 2) += fac_gamma_dt_tauC_derxy_z_ui*derxy_(2,vi) ;
            } // end loop rows vi (test functions for matrix)
          } // end loop columns ui (solution for matrix, test function for vector)
#ifdef PERF
          timeelecstab_ref = null;
#endif
        } // end cstab


        //---------------------------------------------------------------
        //
        //               QUASISTATIC STABILISATION PART
        //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
        //
        //---------------------------------------------------------------
        if(cross == Fluid3::cross_stress_stab)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif
          const double fac_afgdt_tauM = fac*afgdt*tauM;
          for (int ui=0; ui<iel_; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const double fac_afgdt_tauM_conv_resM_ui = fac_afgdt_tauM*conv_resM_(ui);

            for (int vi=0; vi<iel_; ++vi)  // loop rows (test functions for matrix)
            {
              const double fac_afgdt_tauM_conv_resM_ui_funct_vi = fac_afgdt_tauM_conv_resM_ui*funct_(vi);

              /*  factor:

               -alphaF*gamma*dt*tauM

                          /                          \
                         |  /            \            |
                         | | resM o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
              */
              elemat(vi*4    , ui*4    ) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
              elemat(vi*4 + 1, ui*4 + 1) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
              elemat(vi*4 + 2, ui*4 + 2) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
#ifdef PERF
          timeelecrossrey_ref = null;
#endif
        } // end cross
      } // end if compute_elemat

      //---------------------------------------------------------------
      //---------------------------------------------------------------
      //
      //          RIGHT HAND SIDE, QUASISTATIC SUBGRID SCALES
      //
      //---------------------------------------------------------------
      //---------------------------------------------------------------

      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelegalerkin_ref = rcp(new TimeMonitor(*timeelegalerkin));
#endif

        {
          const double fac_inertia_conv_and_bodyforce_x = fac*(accintam_(0)+convaf_old_(0)-edeadaf_(0));
          const double fac_inertia_conv_and_bodyforce_y = fac*(accintam_(1)+convaf_old_(1)-edeadaf_(1));
          const double fac_inertia_conv_and_bodyforce_z = fac*(accintam_(2)+convaf_old_(2)-edeadaf_(2));

          for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
          {

            /* inertia, convective and dead load terms */

            /*

                /             \
               |     n+am      |
              -|  acc     , v  |
               |               |
                \             /


                /                             \
               |  / n+af       \    n+af       |
              -| | c    o nabla |  u      , v  |
               |  \            /               |
                \                             /

                /           \
               |   n+af      |
              +|  f     , v  |
               |             |
                \           /
            */

            elevec[ui*4    ] -= fac_inertia_conv_and_bodyforce_x*funct_(ui);
            elevec[ui*4 + 1] -= fac_inertia_conv_and_bodyforce_y*funct_(ui);
            elevec[ui*4 + 2] -= fac_inertia_conv_and_bodyforce_z*funct_(ui);
          }
        }

        {
          const double fac_divunp = fac*divunp;

          for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
          {

            /* continuity equation */

            /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
            */

            elevec[ui*4 + 3] -= fac_divunp*funct_(ui);
          }
        }

        {
          const double fac_prenp = fac*prenp_;

          for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
          {
            /* pressure */

            /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
            */

            elevec[ui*4    ] += fac_prenp*derxy_(0,ui) ;
            elevec[ui*4 + 1] += fac_prenp*derxy_(1,ui) ;
            elevec[ui*4 + 2] += fac_prenp*derxy_(2,ui) ;
          }
        }


        {
          const double visceff_fac = visceff*fac;
          for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
          {
            /* viscous term */

            /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
            */

            elevec[ui*4    ] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                                 +
                                 derxy_(1,ui)*vderxyaf_(0,1)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,0)
                                 +
                                 derxy_(2,ui)*vderxyaf_(0,2)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,0)) ;
            elevec[ui*4 + 1] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,1)
                                 +
                                 derxy_(0,ui)*vderxyaf_(1,0)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,1)*2.0
                                 +
                                 derxy_(2,ui)*vderxyaf_(1,2)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,1)) ;
            elevec[ui*4 + 2] -= visceff_fac*
                                (derxy_(0,ui)*vderxyaf_(0,2)
                                 +
                                 derxy_(0,ui)*vderxyaf_(2,0)
                                 +
                                 derxy_(1,ui)*vderxyaf_(1,2)
                                 +
                                 derxy_(1,ui)*vderxyaf_(2,1)
                                 +
                                 derxy_(2,ui)*vderxyaf_(2,2)*2.0) ;
          } // end loop rows
        }

#ifdef PERF
        timeelegalerkin_ref = null;
#endif

      }

      if(pspg == Fluid3::pstab_use_pspg)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelepspg_ref = rcp(new TimeMonitor(*timeelepspg));
#endif

        const double fac_tauMp = fac*tauMp;

        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
          /*
            factor: +tauMp

            pressure stabilisation --- inertia


                  /                  \
                 |     n+am           |
                 |  acc    , nabla q  |
                 |                    |
                  \                  /

            pressure stabilisation --- convection


                  /                                   \
                 |  / n+af       \    n+af             |
                 | | c    o nabla |  u      , nabla q  |
                 |  \            /                     |
                  \                                   /


            pressure stabilisation --- diffusion

                  /                                  \
                 |               / n+af \             |
                 |  nabla o eps | u      | , nabla q  |
                 |               \      /             |
                  \                                  /

            pressure stabilisation --- pressure

                  /                      \
                 |         n+1            |
                 |  nabla p    , nabla q  |
                 |                        |
                  \                      /


            pressure stabilisation --- bodyforce
                  /                 \
                 |    n+af           |
                 |  f     , nabla q  |
                 |                   |
                  \                 /
          */
          elevec[ui*4 + 3] -= fac_tauMp*
                              (derxy_(0,ui)*resM_(0)
                               +
                               derxy_(1,ui)*resM_(1)
                               +
                               derxy_(2,ui)*resM_(2));
        } // end loop rows

#ifdef PERF
        timeelepspg_ref = null;
#endif
      }

      if(supg == Fluid3::convective_stab_supg)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelesupg_ref = rcp(new TimeMonitor(*timeelesupg));
#endif

        const double fac_tauM = fac*tauM;

        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
          const double fac_tauM_conv_c_af_ui = fac_tauM*conv_c_af_(ui);
          /*
            factor: +tauM

            SUPG stabilisation --- inertia

                  /                              \
                 |     n+am   / n+af        \     |
                 |  acc    , | c     o nabla | v  |
                 |            \             /     |
                  \                              /

           SUPG stabilisation --- convection


                  /                                                \
                 |    / n+af        \   n+af    / n+af        \     |
                 |   | c     o nabla | u     , | c     o nabla | v  |
                 |    \             /           \             /     |
                  \                                                /


           SUPG stabilisation --- diffusion

                  /                                               \
                 |               / n+af \      / n+af        \     |
                 |  nabla o eps | u      |  , | c     o nabla | v  |
                 |               \      /      \             /     |
                  \                                               /

           SUPG stabilisation --- pressure

                  /                                  \
                 |         n+1    / n+af        \     |
                 |  nabla p    , | c     o nabla | v  |
                 |                \             /     |
                  \                                  /

           SUPG stabilisation --- bodyforce


                  /                             \
                 |   n+af    / n+af        \     |
                 |  f     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /
          */

          elevec[ui*4    ] -= fac_tauM_conv_c_af_ui*resM_(0) ;
          elevec[ui*4 + 1] -= fac_tauM_conv_c_af_ui*resM_(1) ;
          elevec[ui*4 + 2] -= fac_tauM_conv_c_af_ui*resM_(2) ;

        } // end loop rows
#ifdef PERF
        timeelesupg_ref = null;
#endif
      }

      if (higher_order_ele)
      {
        if(vstab != Fluid3::viscous_stab_none && higher_order_ele)
        {
#ifdef PERF
          RefCountPtr<TimeMonitor> timeelevstab_ref = rcp(new TimeMonitor(*timeelevstab));
#endif

          const double fac_two_visc_tauMp = vstabfac * fac*2.0*visc*tauMp;

          for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
          {
            /*
              factor: -(+)tauMp*2*nu


              viscous stabilisation --- inertia


                 /                         \
                |      n+am                 |
                |  Dacc      , div eps (v)  |
                |                           |
                 \                         /

            viscous stabilisation --- convection

            /                                     \
           |  / n+af       \    n+af               |
           | | c    o nabla |  u     , div eps (v) |
           |  \            /                       |
            \                                     /

            viscous stabilisation --- diffusion

               /                                      \
              |               /  n+af \                |
              |  nabla o eps |  u      | , div eps (v) |
              |               \       /                |
               \                                      /

            viscous stabilisation --- pressure

                 /                           \
                |                             |
                |  nabla p , nabla o eps (v)  |
                |                             |
                 \                           /

           viscous stabilisation --- bodyforce

                  /                         \
                 |    n+af                   |
                 |  f     ,  nabla o eps (v) |
                 |                           |
                  \                         /
            */
            elevec[ui*4    ] -= fac_two_visc_tauMp*
                                (resM_(0)*viscs2_(0, 0, ui)
                                 +
                                 resM_(1)*viscs2_(0, 1, ui)
                                 +
                                 resM_(2)*viscs2_(0, 2, ui)) ;
            elevec[ui*4 + 1] -= fac_two_visc_tauMp*
                                (resM_(0)*viscs2_(0, 1, ui)
                                 +
                                 resM_(1)*viscs2_(1, 1, ui)
                                 +
                                 resM_(2)*viscs2_(1, 2, ui)) ;
            elevec[ui*4 + 2] -= fac_two_visc_tauMp*
                                (resM_(0)*viscs2_(0, 2, ui)
                                 +
                                 resM_(1)*viscs2_(1, 2, ui)
                                 +
                                 resM_(2)*viscs2_(2, 2, ui)) ;
          } // end loop rows ui
#ifdef PERF
          timeelevstab_ref = null;
#endif
        } // endif (a)gls
      }

      if(cstab == Fluid3::continuity_stab_yes)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecstab_ref = rcp(new TimeMonitor(*timeelecstab));
#endif

        const double fac_tauC_divunp = fac*tauC*divunp;
        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
          /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
          */

          elevec[ui*4    ] -= fac_tauC_divunp*derxy_(0,ui) ;
          elevec[ui*4 + 1] -= fac_tauC_divunp*derxy_(1,ui) ;
          elevec[ui*4 + 2] -= fac_tauC_divunp*derxy_(2,ui) ;
        } // end loop rows

#ifdef PERF
        timeelecstab_ref = null;
#endif

      } // end cstab

      if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif

        const double fac_tauM = fac*tauM;
        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
        const double fac_tauM_funct_ui = fac_tauM*funct_(ui);

          /* factor: +tauM

                  /                            \
                 |                    n+af      |
                 |  ( resM o nabla ) u    ,  v  |
                 |                    (i)       |
                  \                            /
          */

          elevec[ui*4    ] += fac_tauM_funct_ui*(resM_(0)*vderxyaf_(0,0)
                                                 +
                                                 resM_(1)*vderxyaf_(0,1)
                                                 +
                                                 resM_(2)*vderxyaf_(0,2));
          elevec[ui*4 + 1] += fac_tauM_funct_ui*(resM_(0)*vderxyaf_(1,0)
                                                 +
                                                 resM_(1)*vderxyaf_(1,1)
                                                 +
                                        resM_(2)*vderxyaf_(1,2));
          elevec[ui*4 + 2] += fac_tauM_funct_ui*(resM_(0)*vderxyaf_(2,0)
                                                 +
                                                 resM_(1)*vderxyaf_(2,1)
                                                 +
                                                 resM_(2)*vderxyaf_(2,2));

        } // end loop rows
#ifdef PERF
        timeelecrossrey_ref = null;
#endif
      }

      if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
      {
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelecrossrey_ref = rcp(new TimeMonitor(*timeelecrossrey));
#endif

        const double fac_tauM_tauM = fac*tauM*tauM;
        for (int ui=0; ui<iel_; ++ui) // loop rows  (test functions)
        {
          const double fac_tauM_tauM_conv_resM_ui = fac_tauM_tauM*conv_resM_(ui);

          /* factor: -tauM*tauM

                  /                             \
                 |                               |
                 |  resM   , ( resM o nabla ) v  |
                 |                               |
                  \                             /
          */
          elevec[ui*4    ] += fac_tauM_tauM_conv_resM_ui*resM_(0);
          elevec[ui*4 + 1] += fac_tauM_tauM_conv_resM_ui*resM_(1);
          elevec[ui*4 + 2] += fac_tauM_tauM_conv_resM_ui*resM_(2);
        } // end loop rows
#ifdef PERF
        timeelecrossrey_ref = null;
#endif
      }

      if(fssgv == Fluid3::fssgv_scale_similarity ||
         fssgv == Fluid3::fssgv_mixed_Smagorinsky_all ||
         fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
      {
        //----------------------------------------------------------------------
        //     SCALE-SIMILARITY TERM (ON RIGHT HAND SIDE)

        for (int ui=0; ui<iel_; ++ui)
        {
          elevec[ui*4    ] -= fac*(csconvintaf_(0) - convaf_s_(0))*funct_(ui);
          elevec[ui*4 + 1] -= fac*(csconvintaf_(1) - convaf_s_(1))*funct_(ui);
          elevec[ui*4 + 2] -= fac*(csconvintaf_(2) - convaf_s_(2))*funct_(ui);
        }
      }

      if(fssgv != Fluid3::fssgv_no && fssgv != Fluid3::fssgv_scale_similarity)
      {
        //----------------------------------------------------------------------
        //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

        for (int ui=0; ui<iel_; ++ui)
        {
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                                  /                              \
                         n+af    |       /    n+af\         / \   |
             - nu_art(fsu    ) * |  eps | Dfsu     | , eps | v |  |
                                 |       \        /         \ /   |
                                  \                              /
          */
          elevec[ui*4    ] -= vartfac*(2.0*derxy_(0, ui)*fsvderxyaf_(0, 0)
                                      +    derxy_(1, ui)*fsvderxyaf_(0, 1)
                                      +    derxy_(1, ui)*fsvderxyaf_(1, 0)
                                      +    derxy_(2, ui)*fsvderxyaf_(0, 2)
                                      +    derxy_(2, ui)*fsvderxyaf_(2, 0)) ;
          elevec[ui*4 + 1] -= vartfac*(    derxy_(0, ui)*fsvderxyaf_(0, 1)
                                      +    derxy_(0, ui)*fsvderxyaf_(1, 0)
                                      +2.0*derxy_(1, ui)*fsvderxyaf_(1, 1)
                                      +    derxy_(2, ui)*fsvderxyaf_(1, 2)
                                      +    derxy_(2, ui)*fsvderxyaf_(2, 1)) ;
          elevec[ui*4 + 2] -= vartfac*(    derxy_(0, ui)*fsvderxyaf_(0, 2)
                                      +    derxy_(0, ui)*fsvderxyaf_(2, 0)
                                      +    derxy_(1, ui)*fsvderxyaf_(1, 2)
                                      +    derxy_(1, ui)*fsvderxyaf_(2, 1)
                                      +2.0*derxy_(2, ui)*fsvderxyaf_(2, 2)) ;
        }
      }
    }
  } // end loop iquad
  return;
}

// this is just for comparison of dynamic/quasistatic subscales --- NOT for
// the comparison with physical turbulence models (Smagorinsky etc.)

void DRT::ELEMENTS::Fluid3GenalphaResVMM::CalcRes(
  Fluid3*                                               ele,
  const blitz::Array<double,2>&                         edispnp,
  const blitz::Array<double,2>&                         egridvaf,
  const blitz::Array<double,2>&                         evelnp,
  const blitz::Array<double,1>&                         eprenp,
  const blitz::Array<double,2>&                         eaccam,
  const blitz::Array<double,2>&                         evelaf,
  const struct _MATERIAL*                               material,
  const double                                          alphaM,
  const double                                          alphaF,
  const double                                          gamma,
  const double                                          dt,
  const double                                          time,
  const enum Fluid3::StabilisationAction                tds,
  const enum Fluid3::TauType                            whichtau,
  const bool                                            higher_order_ele,
  blitz::Array<double,1>&                               mean_res,
  blitz::Array<double,1>&                               mean_sacc,
  blitz::Array<double,1>&                               mean_svelaf,
  blitz::Array<double,1>&                               mean_res_sq,
  blitz::Array<double,1>&                               mean_sacc_sq,
  blitz::Array<double,1>&                               mean_svelaf_sq,
  double&                                               mean_resC,
  double&                                               mean_spreacc,
  double&                                               mean_sprenp,
  double&                                               mean_resC_sq,
  double&                                               mean_spreacc_sq,
  double&                                               mean_sprenp_sq
  )
{
  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  blitz::Range       _ = blitz::Range::all();

  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt

  const double timealphaF = time-(1-alphaF)*dt;

  //------------------------------------------------------------------
  //                      SET MATERIAL DATA
  //------------------------------------------------------------------
  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  //------------------------------------------------------------------
  //                      SET ELEMENT DATA
  //------------------------------------------------------------------
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
    for (int inode=0; inode<iel_; inode++)
    {
      xyze_(0,inode) += edispnp(0,inode);
      xyze_(1,inode) += edispnp(1,inode);
      xyze_(2,inode) += edispnp(2,inode);
    }
  }

  // dead load in element nodes
  GetNodalBodyForce(ele,timealphaF);

  //----------------------------------------------------------------------------
  //            STABILIZATION PARAMETER, SMAGORINSKY MODEL
  //      and everything else that is evaluated in the element center
  //
  // This has to be done before anything else is calculated because we use
  // the same arrays internally.
  //----------------------------------------------------------------------------

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili=DRT::UTILS::intrule3D_undefined;
  switch (distype)
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
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  const double wquad = intpoints_onepoint.qwgt[0];

  DRT::UTILS::shape_function_3D       (funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // get Jacobian matrix and determinant
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        xjm_(nn,rr)+=deriv_(nn,mm)*xyze_(rr,mm);
      }
    }
  }

  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  vol_ = wquad*det;

  // get element length for tau_M and tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol_/PI),(1.0/3.0))/sqrt(3.0);

  //
  //             compute global first derivates
  //
  // this is necessary only for the calculation of the
  // streamlength (required by the quasistatic formulation) and
  // the Smagorinsky model.
  //
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

  */

  // inverse of jacobian (transposed)
  /*
          +-                 -+     +-                 -+ -1
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dx    dx    dx   |     |  dr    dr    dr   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |  =  |  --    --    --   |
          |  dy    dy    dy   |     |  ds    ds    ds   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dz    dz    dz   |     |  dt    dt    dt   |
          +-                 -+     +-                 -+

  */
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
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<iel_;++rr)
    {
      derxy_(nn,rr)=deriv_(0,rr)*xji_(nn,0);
      for (int mm=1;mm<3;++mm)
      {
        derxy_(nn,rr)+=deriv_(mm,rr)*xji_(nn,mm);
      }
    }
  }

  // get velocities (n+alpha_F,i) at integration point
  //
  //                 +-----
  //       n+af       \                  n+af
  //    vel    (x) =   +      N (x) * vel
  //                  /        j         j
  //                 +-----
  //                 node j
  //
  for (int rr=0;rr<3;++rr)
  {
    velintaf_(rr)=funct_(0)*evelaf(rr,0);
    for (int mm=1;mm<iel_;++mm)
    {
      velintaf_(rr)+=funct_(mm)*evelaf(rr,mm);
    }
  }

  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //       n+af      +-----  dN (x)
  //   dvel    (x)    \        k         n+af
  //   ----------- =   +     ------ * vel
  //       dx         /        dx        k
  //         j       +-----      j
  //                 node k
  //
  // j : direction of derivative x/y/z
  //
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      vderxyaf_(nn,rr)=derxy_(rr,0)*evelaf(nn,0);
      for (int mm=1;mm<iel_;++mm)
      {
        vderxyaf_(nn,rr)+=derxy_(rr,mm)*evelaf(nn,mm);
      }
    }
  }


  // get velocities (n+1,i)  at integration point
  //
  //                +-----
  //       n+1       \                  n+1
  //    vel   (x) =   +      N (x) * vel
  //                 /        j         j
  //                +-----
  //                node j
  //
  for (int rr=0;rr<3;++rr)
  {
    velintnp_(rr)=funct_(0)*evelnp(rr,0);
    for (int mm=1;mm<iel_;++mm)
    {
      velintnp_(rr)+=funct_(mm)*evelnp(rr,mm);
    }
  }

  // get velocity norms
  const double vel_normaf = sqrt(velintaf_(0)*velintaf_(0)
                                 +
                                 velintaf_(1)*velintaf_(1)
                                 +
                                 velintaf_(2)*velintaf_(2));
  const double vel_normnp = sqrt(velintnp_(0)*velintnp_(0)
                                 +
                                 velintnp_(1)*velintnp_(1)
                                 +
                                 velintnp_(2)*velintnp_(2));

  if(tds == Fluid3::subscales_time_dependent)
  {
    //-------------------------------------------------------
    //          TAUS FOR TIME DEPENDENT SUBSCALES
    //-------------------------------------------------------

    if(whichtau == Fluid3::bazilevs)
    {
      /* INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA

         tau_M: Bazilevs et al. + ideas from Codina
                                                         1.0
                 +-                                 -+ - ---
                 |                                   |   2.0
             td  |  n+af      n+af         2         |
          tau  = | u     * G u     + C * nu  * G : G |
             M   |         -          I        -   - |
                 |         -                   -   - |
                 +-                                 -+

         tau_C: Bazilevs et al., derived from the fine scale complement Shur
                                 operator of the pressure equation


                       td         1.0
                    tau  = -----------------
                       C       td   /     \
                            tau  * | g * g |
                               M    \-   -/
      */

      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      blitz::Array<double,2> G(3,3,blitz::ColumnMajorArray<2>());

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
          }
        }
      }

      /*          +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
      */
      double normG = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          normG+=G(nn,rr)*G(nn,rr);
        }
      }

      /*                    +----
           n+af      n+af    \     n+af         n+af
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
      */
      double Gnormu = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          Gnormu+=velintaf_(nn)*G(nn,rr)*velintaf_(rr);
        }
      }

      // definition of constant from Akkerman et al. (2008)
      const double CI = 36.0;

      /*                                                 1.0
                 +-                                 -+ - ---
                 |                                   |   2.0
                 |  n+af      n+af         2         |
          tau  = | u     * G u     + C * nu  * G : G |
             M   |         -          I        -   - |
                 |         -                   -   - |
                 +-                                 -+
      */
      tau_(0) = 1.0/sqrt(Gnormu+CI*visc*visc*normG);
      tau_(1) = tau_(0);

      /*         +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
      */
      blitz::Array<double,1> g(3);

      for (int rr=0;rr<3;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          g(rr) += xji_(rr,mm);
        }
      }

      /*         +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
      */
      const double normgsq = g(0)*g(0)+g(1)*g(1)+g(2)*g(2);

      /*
                                1.0
                  tau  = -----------------
                     C           /      \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);

    }
    else if(whichtau == Fluid3::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      //
      // tau_M: modification of
      //
      //    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
      //    Finite Element Method for the Advective-Reactive-Diffusive
      //    Equation. Computer Methods in Applied Mechanics and Enginnering,
      //    Vol. 190, pp. 1785-1800, 2000.
      //    http://www.lncc.br/~valentin/publication.htm                   */
      //
      // tau_Mp: modification of Barrenechea, G.R. and Valentin, F.
      //
      //    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
      //    element method for a generalized Stokes problem. Numerische
      //    Mathematik, Vol. 92, pp. 652-677, 2002.
      //    http://www.lncc.br/~valentin/publication.htm
      //
      //
      // tau_C: kept Wall definition
      //
      // for the modifications see Codina, Principe, Guasch, Badia
      //    "Time dependent subscales in the stabilized finite  element
      //     approximation of incompressible flow problems"
      //
      //
      // see also: Codina, R. and Soto, O.: Approximation of the incompressible
      //    Navier-Stokes equations using orthogonal subscale stabilisation
      //    and pressure segregation on anisotropic finite element meshes.
      //    Computer methods in Applied Mechanics and Engineering,
      //    Vol 193, pp. 1403-1419, 2004.

      //---------------------------------------------- compute tau_Mu = tau_Mp
      /* convective : viscous forces (element reynolds number)*/
      const double re_convectaf = (vel_normaf * hk / visc ) * (mk/2.0);

      const double xi_convectaf = DMAX(re_convectaf,1.0);

      /*
               xi_convect ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re_convect
                              1
      */

      /* the 4.0 instead of the Franca's definition 2.0 results from the viscous
       * term in the Navier-Stokes-equations, which is scaled by 2.0*nu         */

      tau_(0) = DSQR(hk) / (4.0 * visc / mk + ( 4.0 * visc/mk) * xi_convectaf);

      /*------------------------------------------------------ compute tau_C ---*/

      //-- stability parameter definition according to Wall Diss. 99
      /*
               xi_convect ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> Re_convect
                              1
      */
      const double re_convectnp = (vel_normnp * hk / visc ) * (mk/2.0);

      const double xi_tau_c = DMIN(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else if(whichtau == Fluid3::codina)
    {
      // Parameter from Codina, Badia (Constants are chosen according to
      // the values in the standard definition above)

      const double CI  = 4.0/mk;
      const double CII = 2.0/mk;

      // in contrast to the original definition, we neglect the influence of
      // the subscale velocity on velnormaf
      tau_(0)=1.0/(CI*visc/(hk*hk)+CII*vel_normaf/hk);

      tau_(1)=tau_(0);

      tau_(2)=(hk*hk)/(CI*tau_(0));
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter\n");
    }
  }
  else
  {
    //-------------------------------------------------------
    //        TAUS FOR THE QUASISTATIC FORMULATION
    //-------------------------------------------------------

    if(whichtau == Fluid3::bazilevs)
    {
      /* INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA

         tau_M: Bazilevs et al.
                                                               1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+af      n+af         2         |
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

      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      blitz::Array<double,2> G(3,3,blitz::ColumnMajorArray<2>());

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
          }
        }
      }

      /*          +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
      */
      double normG = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          normG+=G(nn,rr)*G(nn,rr);
        }
      }

      /*                    +----
           n+af      n+af    \     n+af         n+af
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
      */
      double Gnormu = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          Gnormu+=velintaf_(nn)*G(nn,rr)*velintaf_(rr);
        }
      }

      // definition of constant from Akkerman et al. (2008)
      const double CI = 36.0;

      /*                                                       1.0
                 +-                                       -+ - ---
                 |                                         |   2.0
                 | 4.0    n+af      n+af         2         |
          tau  = | --- + u     * G u     + C * nu  * G : G |
             M   |   2           -          I        -   - |
                 | dt            -                   -   - |
                 +-                                       -+
      */
      tau_(0) = 1.0/sqrt(4.0/(dt*dt)+Gnormu+CI*visc*visc*normG);
      tau_(1) = tau_(0);

      /*         +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
      */
      blitz::Array<double,1> g(3);

      for (int rr=0;rr<3;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          g(rr) += xji_(rr,mm);
        }
      }

      /*         +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
      */
      const double normgsq = g(0)*g(0)+g(1)*g(1)+g(2)*g(2);

      /*
                                1.0
                  tau  = -----------------
                     C           /      \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);
    }
    else if (whichtau == Fluid3::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Wall


      // this copy of velintaf_ will be used to store the normed velocity
      blitz::Array<double,1> normed_velintaf(3);
      normed_velintaf=velintaf_.copy();

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        normed_velintaf = velintaf_/vel_normaf;
      }
      else
      {
        normed_velintaf    = 0.;
        normed_velintaf(0) = 1.;
      }

      // get streamlength
      const double val = blitz::sum(blitz::abs(blitz::sum(normed_velintaf(j)*derxy_(j,i),j)));
      const double strle = 2.0/val;

      // time factor
      const double timefac = gamma*dt;

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


      const double re1 = 4.0 * timefac * visc / (mk * DSQR(strle));   /* viscous : reactive forces   */
      const double re2 = mk * vel_normaf * strle / (2.0 * visc);      /* convective : viscous forces */

      const double xi1 = DMAX(re1,1.0);
      const double xi2 = DMAX(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visc/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visc);    /* convective : viscous forces */

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
      tau_(1) = timefac * DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visc/mk) * xi_convect);

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
      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;

    }
    else if(whichtau == Fluid3::codina)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Codina


      // this copy of velintaf_ will be used to store the normed velocity
      blitz::Array<double,1> normed_velintaf(3);
      normed_velintaf=velintaf_.copy();

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        normed_velintaf = velintaf_/vel_normaf;
      }
      else
      {
        normed_velintaf    = 0.;
        normed_velintaf(0) = 1.;
      }

      // get streamlength
      const double val = blitz::sum(blitz::abs(blitz::sum(normed_velintaf(j)*derxy_(j,i),j)));
      const double strle = 2.0/val;

      // time factor
      const double timefac = gamma*dt;

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


      const double re1 = 4.0 * timefac * visc / (mk * DSQR(strle));   /* viscous : reactive forces   */
      const double re2 = mk * vel_normaf * strle / (2.0 * visc);      /* convective : viscous forces */

      const double xi1 = DMAX(re1,1.0);
      const double xi2 = DMAX(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visc/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visc);    /* convective : viscous forces */

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
      tau_(1) = timefac * DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visc/mk) * xi_convect);

      /*------------------------------------------------------ compute tau_C ---*/
      /*-- stability parameter definition according to Codina (2002), CMAME 191
       *
       * Analysis of a stabilized finite element approximation of the transient
       * convection-diffusion-reaction equation using orthogonal subscales.
       * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
       *
       * */
      tau_(2) = sqrt(DSQR(visc)+DSQR(0.5*vel_normnp*hk));
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter\n");
    }
  }

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // get subscale information from element --- this is just a reference
  // to the element data
  blitz::Array<double,2> saccn (ele->sub_acc_old_);
  blitz::Array<double,2> sveln (ele->sub_vel_old_);
  blitz::Array<double,2> svelnp(ele->sub_vel_    );
  blitz::Array<double,1> spren (ele->sub_pre_old_);
  blitz::Array<double,1> sprenp(ele->sub_pre_    );

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {

    // set gauss point coordinates
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

    // get values of shape functions and derivatives in the gausspoint
    DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
    }

    // get transposed Jacobian matrix and determinant
    //
    //        +-            -+ T      +-            -+
    //        | dx   dx   dx |        | dx   dy   dz |
    //        | --   --   -- |        | --   --   -- |
    //        | dr   ds   dt |        | dr   dr   dr |
    //        |              |        |              |
    //        | dy   dy   dy |        | dx   dy   dz |
    //        | --   --   -- |   =    | --   --   -- |
    //        | dr   ds   dt |        | ds   ds   ds |
    //        |              |        |              |
    //        | dz   dz   dz |        | dx   dy   dz |
    //        | --   --   -- |        | --   --   -- |
    //        | dr   ds   dt |        | dt   dt   dt |
    //        +-            -+        +-            -+
    //
    // The Jacobian is computed using the formula
    //
    //            +-----
    //   dx_j(r)   \      dN_k(r)
    //   -------  = +     ------- * (x_j)_k
    //    dr_i     /       dr_i       |
    //            +-----    |         |
    //            node k    |         |
    //                  derivative    |
    //                   of shape     |
    //                   function     |
    //                           component of
    //                          node coordinate
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
        for (int mm=1;mm<iel_;++mm)
        {
          xjm_(nn,rr)+=deriv_(nn,mm)*xyze_(rr,mm);
        }
      }
    }
    // The determinant ist computed using Sarrus's rule
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                       xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                       xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                       xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                       xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                       xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

    // check for degenerated elements
    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
    }

    //--------------------------------------------------------------
    //             compute global first derivates
    //--------------------------------------------------------------
    /*
      Use the Jacobian and the known derivatives in element coordinate
      directions on the right hand side to compute the derivatives in
      global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+
    */

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
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<iel_;++rr)
      {
        derxy_(nn,rr)=deriv_(0,rr)*xji_(nn,0);
        for (int mm=1;mm<3;++mm)
        {
          derxy_(nn,rr)+=deriv_(mm,rr)*xji_(nn,mm);
        }
      }
    }

    //--------------------------------------------------------------
    //             compute second global derivative
    //--------------------------------------------------------------

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
    if (higher_order_ele)
    {
      // initialize and zero out everything
      blitz::Array<double,2> bm(6,6,blitz::ColumnMajorArray<2>());

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
       |        | | | | 1           | | | | | 1     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 2           | | | | | 2     +-+-+-+
       |        +-+-+-+       =     +-+-+-+-+    *  | | | | .
       |        | | | | 3           | | | | | 3     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 4           | | | | | 4     +-+-+-+ .
       |        +-+-+-+             +-+-+-+-+       | | | | .
       |        | | | | 5           | | | | | 5     +-+-+-+
       |        +-+-+-+             +-+-+-+-+       | | | | iel-1
       |                                            +-+-+-+
       |
       |        xder2               deriv2          xyze^T
       |
       |
       |                                     +-                  -+
       |  	   	    	    	     | d^2x   d^2y   d^2z |
       |  	   	    	    	     | ----   ----   ---- |
       | 	   	   	   	     | dr^2   dr^2   dr^2 |
       | 	   	   	   	     |                    |
       | 	   	   	   	     | d^2x   d^2y   d^2z |
       |                                     | ----   ----   ---- |
       | 	   	   	   	     | ds^2   ds^2   ds^2 |
       | 	   	   	   	     |                    |
       | 	   	   	   	     | d^2x   d^2y   d^2z |
       | 	   	   	   	     | ----   ----   ---- |
       | 	   	   	   	     | dt^2   dt^2   dt^2 |
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
       | 	   	   	   	     +-                  -+
       |
       |
      */
      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          xder2_(nn,rr)= deriv2_(nn,0)*xyze_(rr,0);

          for (int mm=1;mm<iel_;++mm)
          {
            xder2_(nn,rr)+= deriv2_(nn,mm)*xyze_(rr,mm);
          }
        }
      }

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
      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<iel_;++rr)
        {
          derxy2_(nn,rr)= (-1)*xder2_(nn,0)*derxy_(0,rr);

          for (int mm=1;mm<3;++mm)
          {
            derxy2_(nn,rr)-= xder2_(nn,mm)*derxy_(mm,rr);
          }
        }
      }

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

      for (int nn=0;nn<6;++nn)
      {
        for (int rr=0;rr<iel_;++rr)
        {
          derxy2_(nn,rr) += deriv2_(nn,rr);
        }
      }

      /* make LU decomposition and solve system for all right hand sides
       * (i.e. the components of chainrulerhs)
       |
       |             0  1  2  3  4  5         i        i
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 0     | | 0    | | 0
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 1     | | 1    | | 1
       | 	   +--+--+--+--+--+--+       +-+      +-+
       | 	   |  |  |  |  |  |  | 2     | | 2    | | 2
       | 	   +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
       |           |  |  |  |  |  |  | 3     | | 3    | | 3
       |           +--+--+--+--+--+--+       +-+      +-+
       |           |  |  |  |  |  |  | 4     | | 4    | | 4
       |           +--+--+--+--+--+--+       +-+      +-+
       |           |  |  |  |  |  |  | 5     | | 5    | | 5
       |           +--+--+--+--+--+--+       +-+      +-+
       |                                      |        |
       |                                      |        |
       |                                      derxy2[i]|
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
       |    	       	      +-+-+-+-+
      */
      // Use LAPACK
      Epetra_LAPACK          solver;

      // a vector specifying the pivots (reordering)
      int pivot[6];

      // error code
      int ierr = 0;

      // Perform LU factorisation --- this call replaces bm with its factorisation
      solver.GETRF(6,6,bm.data(),6,&(pivot[0]),&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform LU factorisation during computation of derxy2");
      }

      // backward substitution. GETRS replaces the input (chainrulerhs, currently
      // stored on derxy2) with the result
      solver.GETRS('N',6,iel_,bm.data(),6,&(pivot[0]),derxy2_.data(),6,&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform backward substitution after factorisation of jacobian");
      }
    }

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------

    // get intermediate accelerations (n+alpha_M,i) at integration point
    //
    //                 +-----
    //       n+am       \                  n+am
    //    acc    (x) =   +      N (x) * acc
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    // i         : space dimension u/v/w
    //
    for (int rr=0;rr<3;++rr)
    {
      accintam_(rr)=funct_(0)*eaccam(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        accintam_(rr)+=funct_(mm)*eaccam(rr,mm);
      }
    }

    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    for (int rr=0;rr<3;++rr)
    {
      velintaf_(rr)=funct_(0)*evelaf(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        velintaf_(rr)+=funct_(mm)*evelaf(rr,mm);
      }
    }

    // get bodyforce in gausspoint, time (n+alpha_F)
    //
    //                 +-----
    //       n+af       \                n+af
    //      f    (x) =   +      N (x) * f
    //                  /        j       j
    //                 +-----
    //                 node j
    //
    for (int rr=0;rr<3;++rr)
    {
      bodyforceaf_(rr)=funct_(0)*edeadaf_(rr,0);
      for (int mm=1;mm<iel_;++mm)
      {
        bodyforceaf_(rr)+=funct_(mm)*edeadaf_(rr,mm);
      }
    }

    // get pressure gradient (n+1,i) at integration point
    //
    //       n+1      +-----  dN (x)
    //   dpre   (x)    \        j         n+1
    //   ---------- =   +     ------ * pre
    //       dx        /        dx        j
    //         i      +-----      i
    //                node j
    //
    // i : direction of derivative
    //
    for (int rr=0;rr<3;++rr)
    {
      pderxynp_(rr)=derxy_(rr,0)*eprenp(0);
      for (int mm=1;mm<iel_;++mm)
      {
        pderxynp_(rr)+=derxy_(rr,mm)*eprenp(mm);
      }
    }

    // get velocity (n+alpha_F,i) derivatives at integration point
    //
    //       n+af      +-----  dN (x)
    //   dvel    (x)    \        k         n+af
    //   ----------- =   +     ------ * vel
    //       dx         /        dx        k
    //         j       +-----      j
    //                 node k
    //
    // j : direction of derivative x/y/z
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        vderxyaf_(nn,rr)=derxy_(rr,0)*evelaf(nn,0);
        for (int mm=1;mm<iel_;++mm)
        {
          vderxyaf_(nn,rr)+=derxy_(rr,mm)*evelaf(nn,mm);
        }
      }
    }

    // get velocity (n+1,i) derivatives at integration point
    //
    //       n+1      +-----  dN (x)
    //   dvel   (x)    \        k         n+1
    //   ---------- =   +     ------ * vel
    //       dx        /        dx        k
    //         j      +-----      j
    //                node k
    //
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        vderxynp_(nn,rr)=derxy_(rr,0)*evelnp(nn,0);
        for (int mm=1;mm<iel_;++mm)
        {
          vderxynp_(nn,rr)+=derxy_(rr,mm)*evelnp(nn,mm);
        }
      }
    }

    // get ale convective velocity (n+alpha_F,i) at integration point

    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----           +-                   -+
    //       n+af       \               |   n+af      n+alphaF|
    //      c    (x) =   +      N (x) * |vel     - u_G        |
    //                  /        j      |   j         j       |
    //                 +-----           +-                   -+
    //                 node j
    //
    //
    // where u_G is the grid velocity at the integration point,
    // time (n+alpha_F,i)
    //
    //                 +-----
    //       n+af       \                  n+af
    //    u_G    (x) =   +      N (x) * u_G
    //                  /        j         j
    //                 +-----
    //                 node j
    //

    for (int rr=0;rr<3;++rr)
    {
      aleconvintaf_(rr)=velintaf_(rr);
    }

    if (ele->is_ale_)
    {
      for (int rr=0;rr<3;++rr)
      {
        for (int mm=0;mm<iel_;++mm)
        {
          aleconvintaf_(rr)-=funct_(mm)*egridvaf(rr,mm);
        }
      }
    }

    // calculate 2nd velocity derivatives at integration point, time(n+alpha_F)
    //
    //    2   n+af       +-----   dN (x)
    //   d vel    (x)     \         k          n+af
    //   ------------  =   +     -------- * vel
    //    dx  dx          /      dx  dx        k
    //      j1  j2       +-----    j1  j2
    //                   node k
    //
    // j=(j1,j2) : direction of derivative x/y/z
    if(higher_order_ele)
    {
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<6;++rr)
        {
          vderxy2af_(nn,rr)=derxy2_(rr,0)*evelaf(nn,0);
          for (int mm=1;mm<iel_;++mm)
          {
            vderxy2af_(nn,rr)+=derxy2_(rr,mm)*evelaf(nn,mm);
          }
        }
      }
    }

    /* divergence new time step n+1 */
    const double divunp = (vderxynp_(0,0)+vderxynp_(1,1)+vderxynp_(2,2));

    /* Convective term  u_old * grad u_old: */
    for (int rr=0;rr<3;++rr)
    {
      convaf_old_(rr)=vderxyaf_(rr, 0)*aleconvintaf_(0);
      for (int mm=1;mm<3;++mm)
      {
        convaf_old_(rr)+=vderxyaf_(rr,mm)*aleconvintaf_(mm);
      }
    }

    if (higher_order_ele)
    {
      /* Viscous term  div epsilon(u_old) */
      viscaf_old_(0) = vderxy2af_(0,0) + 0.5 * (vderxy2af_(0,1) + vderxy2af_(1,3) + vderxy2af_(0,2) + vderxy2af_(2,4));
      viscaf_old_(1) = vderxy2af_(1,1) + 0.5 * (vderxy2af_(1,0) + vderxy2af_(0,3) + vderxy2af_(1,2) + vderxy2af_(2,5));
      viscaf_old_(2) = vderxy2af_(2,2) + 0.5 * (vderxy2af_(2,0) + vderxy2af_(0,4) + vderxy2af_(2,1) + vderxy2af_(1,5));
    }

    /* compute residual in gausspoint --- second derivatives only
     * exist for higher order elements,  convaf_old_ is based on
     * ale-convective velocity */
    for (int rr=0;rr<3;++rr)
    {
      resM_(rr) = accintam_(rr) + convaf_old_(rr) + pderxynp_(rr) - bodyforceaf_(rr);
    }

    /* the residual is based on the effective viscosity! */

    if (higher_order_ele)
    {
      const double two_visc=2.0*visc;

      for (int rr=0;rr<3;++rr)
      {
        resM_(rr) -= two_visc*viscaf_old_(rr);
      }
    }

    if(tds == Fluid3::subscales_time_dependent)
    {
      /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

      */
      svelaf_(_) = alphaF*svelnp(_,iquad)+(1.0-alphaF)*sveln(_,iquad);
    }

    for(int rr=0;rr<3;++rr)
    {
      mean_res    (rr) += resM_(rr);
      mean_res_sq (rr) += resM_(rr)*resM_(rr);
    }

    if(tds == Fluid3::subscales_time_dependent)
    {
      for(int rr=0;rr<3;++rr)
      {
        const double aux = -1.0/tau_(0)*svelaf_(rr) -resM_(rr);

        mean_sacc   (rr) += aux;
        mean_sacc_sq(rr) += aux*aux;
      }
    }

    if(tds == Fluid3::subscales_time_dependent)
    {
      for(int rr=0;rr<3;++rr)
      {
        mean_svelaf   (rr) += svelaf_(rr);
        mean_svelaf_sq(rr) += svelaf_(rr)*svelaf_(rr);
      }
    }
    else
    {
      for(int rr=0;rr<3;++rr)
      {
        const double aux = tau_(0)*resM_(rr);

        mean_svelaf   (rr) -= aux;
        mean_svelaf_sq(rr) += aux*aux;
      }
    }

    mean_resC    += divunp;
    mean_resC_sq += divunp*divunp;

    if(tds == Fluid3::subscales_time_dependent)
    {
      const double preaccinc=(sprenp(iquad)-spren(iquad))/dt;

      mean_spreacc   += preaccinc;
      mean_spreacc_sq+= preaccinc*preaccinc;

      mean_sprenp    += sprenp(iquad);

      mean_sprenp_sq += sprenp(iquad)*sprenp(iquad);
    }
    else
    {
      const double aux = tau_(2)*divunp;

      mean_sprenp     -= aux;
      mean_sprenp_sq  += aux*aux;
    }
  }

  mean_res        /= intpoints.nquad;
  mean_res_sq     /= intpoints.nquad;
  mean_sacc       /= intpoints.nquad;
  mean_sacc_sq    /= intpoints.nquad;
  mean_svelaf     /= intpoints.nquad;
  mean_svelaf_sq  /= intpoints.nquad;
  mean_resC       /= intpoints.nquad;
  mean_resC_sq    /= intpoints.nquad;
  mean_spreacc    /= intpoints.nquad;
  mean_spreacc_sq /= intpoints.nquad;
  mean_sprenp     /= intpoints.nquad;
  mean_sprenp_sq  /= intpoints.nquad;


  return;
}




/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3GenalphaResVMM::GetNodalBodyForce(Fluid3* ele, const double time)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

  if (myneumcond.size()>1)
  {
    dserror("more than one VolumeNeumann cond on one node");
  }

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
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
        //curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      for(int isd=0;isd<3;isd++)
      {
        edeadaf_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadaf_ = 0.;
  }
}


//
// calculate material viscosity at gaussian point  u.may 05/08
//
void DRT::ELEMENTS::Fluid3GenalphaResVMM::CalVisc_AtGaussianPoint(
  const struct _MATERIAL*                 material,
  double&                           	  visc)
{

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index

  // compute shear rate
  double rateofshear = 0.0;
  blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());   // strain rate tensor
  epsilon = 0.5 * ( vderxyaf_(i,j) + vderxyaf_(j,i) );

  for(int rr=0;rr<3;rr++)
  	for(int mm=0;mm<3;mm++)
  		rateofshear += epsilon(rr,mm)*epsilon(rr,mm);

  rateofshear = sqrt(2.0*rateofshear);

  if(material->mattyp == m_carreauyasuda)
  {
	double nu_0 	= material->m.carreauyasuda->nu_0;          // parameter for zero-shear viscosity
	double nu_inf   = material->m.carreauyasuda->nu_inf;      	// parameter for infinite-shear viscosity
	double lambda   = material->m.carreauyasuda->lambda;      	// parameter for characteristic time
	double a 		= material->m.carreauyasuda->a_param;  		// constant parameter
	double b 		= material->m.carreauyasuda->b_param;  		// constant parameter

	// compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
	// see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
	const double tmp = pow(lambda*rateofshear,b);
	visc = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
  }
  else if(material->mattyp == m_modpowerlaw)
  {
	// get material parameters
	double m  	  = material->m.modpowerlaw->m_cons;    // consistency constant
	double delta  = material->m.modpowerlaw->delta;       // safety factor
	double a      = material->m.modpowerlaw->a_exp;       // exponent

    // compute viscosity according to a modified power law model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    visc = m * pow((delta + rateofshear), (-1)*a);
  }
  else
	dserror("material type not yet implemented");
}



#endif
#endif
