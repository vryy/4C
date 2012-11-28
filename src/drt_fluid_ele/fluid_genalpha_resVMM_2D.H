/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.cpp

\brief Internal implementation of Fluid2 element with a generalised alpha
       time integration.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


#include "fluid_ele_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  advective version based on quasistatic subgrid scales              |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_adv_qs_2D(
  Fluid*                                ele,
  std::vector<Epetra_SerialDenseVector>& myknots,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
  LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
  const LINALG::Matrix<nsd_,iel>&           edispnp,
  const LINALG::Matrix<nsd_,iel>&           egridvaf,
  const LINALG::Matrix<nsd_,iel>&           evelnp,
  const LINALG::Matrix<iel,1>&           eprenp,
  const LINALG::Matrix<nsd_,iel>&           eaccam,
  const LINALG::Matrix<nsd_,iel>&           evelaf,
  Teuchos::RCP<const MAT::Material>      material,
  const double                           alphaM,
  const double                           alphaF,
  const double                           gamma,
  const double                           dt,
  const double                           time,
  const enum INPAR::FLUID::LinearisationAction newton,
  const bool                             higher_order_ele,
  const enum INPAR::FLUID::Transient         inertia         ,
  const enum INPAR::FLUID::PSPG              pspg            ,
  const enum INPAR::FLUID::SUPG              supg            ,
  const enum INPAR::FLUID::VStab             vstab           ,
  const enum INPAR::FLUID::CStab             cstab           ,
  const enum INPAR::FLUID::CrossStress       cross           ,
  const enum INPAR::FLUID::ReynoldsStress    reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha             whichtau,
  double&                                visceff,
  const bool                             compute_elemat
  )
{
  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt
  //
  const double timealphaF = time-(1-alphaF)*dt;

  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }


  //------------------------------------------------------------------
  //                    SET ALL ELEMENT DATA
  // o including element geometry (node coordinates)
  // o including dead loads in nodes
  // o including hk, mk, element area
  // o including material viscosity, effective viscosity by
  //   Non-Newtonian fluids
  //------------------------------------------------------------------

  double hk   = 0.0;
  double mk   = 0.0;
  double visc = 0.0;

  SetElementData_2D(ele            ,
                 edispnp        ,
                 evelaf         ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 visceff        );

#if 1

  {
    /*
    // use one point gauss rule to calculate volume at element center
    DRT::UTILS::GaussRule2D integrationrule_stabili=DRT::UTILS::intrule2D_undefined;
    switch (distype)
    {
    case DRT::Element::quad4:
    case DRT::Element::nurbs4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs9:
      integrationrule_stabili = DRT::UTILS::intrule_quad_1point;
      break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      integrationrule_stabili = DRT::UTILS::intrule_tri_1point;
      break;
    default:
      dserror("invalid discretization type for fluid2");
    }
*/
    // gaussian points
    //const DRT::UTILS::IntegrationPoints2D intpoints_onepoint(integrationrule_stabili);
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_onepoint(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);


    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      0           ,
      intpoints_onepoint       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_quasistatic,gamma,dt,hk,mk,visceff);
  }
#endif


  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  //const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {

    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);


    /*
         This is the operator

                  /               \
                 | resM    o nabla |
                  \    (i)        /

         required for the cross stress linearisation
    */
    //
    //                    +-----
    //          n+af       \         n+af      dN
    // conv_resM    (x) =   +    resM    (x) * --- (x)
    //                     /         j         dx
    //                    +-----                 j
    //                     dim j
    if(cross == INPAR::FLUID::cross_stress_stab)
    {
      for(int nn=0;nn<iel;++nn)
      {
        conv_resM_(nn)=resM_(0)*derxy_(0,nn);

        for(int rr=1;rr<2;++rr)
        {
          conv_resM_(nn)+=resM_(rr)*derxy_(rr,nn);
        }
      }
    }

#if 0
    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_quasistatic,gamma,dt,hk,mk,visceff);
#endif
    // stabilisation parameters
    const double tauM   = tau_(0);
    const double tauMp  = tau_(1);

    if(cstab == INPAR::FLUID::continuity_stab_none)
    {
      tau_(2)=0.0;
    }

    const double tauC   = tau_(2);

    double supg_active_tauM;
    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      supg_active_tauM=tauM;
    }
    else
    {
      supg_active_tauM=0.0;
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


    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //
    //              SYSTEM MATRIX, QUASISTATIC FORMULATION
    //
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    if(compute_elemat)
    {

      /* get combined convective linearisation (n+alpha_F,i) at
         integration point
         takes care of half of the linearisation of reynolds part
         (if necessary)


                         n+af
         conv_c_plus_svel_   (x) =


                   +-----  /                   \
                    \     |  n+af      ~n+af    |   dN
            = tauM * +    | c    (x) + u    (x) | * --- (x)
                    /     |  j          j       |   dx
                   +-----  \                   /      j
                    dim j
                           +-------+  +-------+
                              if         if
                             supg      reynolds

      */
      for(int nn=0;nn<iel;++nn)
      {
        conv_c_plus_svel_af_(nn)=supg_active_tauM*conv_c_af_(nn);
      }

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {

        /* half of the reynolds linearisation is done by modifiing
           the supg testfunction, see above */

        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(0)*derxy_(0,nn);

          for(int rr=1;rr<2;++rr)
          {
            conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(rr)*derxy_(rr,nn);
          }
        }

        /*
                  /                           \
                 |                             |
                 |  resM , ( resM o nabla ) v  |
                 |                             |
                  \                           /
                            +----+
                              ^
                              |
                              linearisation of this expression
        */
        const double fac_alphaM_tauM_tauM=fac*alphaM*tauM*tauM;

        const double fac_alphaM_tauM_tauM_resM_x=fac_alphaM_tauM_tauM*resM_(0);
        const double fac_alphaM_tauM_tauM_resM_y=fac_alphaM_tauM_tauM*resM_(1);

        const double fac_afgdt_tauM_tauM=fac*afgdt*tauM*tauM;

        double fac_afgdt_tauM_tauM_resM[2];
        fac_afgdt_tauM_tauM_resM[0]=fac_afgdt_tauM_tauM*resM_(0);
        fac_afgdt_tauM_tauM_resM[1]=fac_afgdt_tauM_tauM*resM_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double inertia_and_conv[2];

          inertia_and_conv[0]=fac_afgdt_tauM_tauM_resM[0]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_tauM_tauM_resM[1]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_y*funct_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: -alphaM * tauM * tauM

                  /                           \
                 |                             |
                 |  resM , ( Dacc o nabla ) v  |
                 |                             |
                  \                           /

            */

            /*
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / / n+af        \       \         \     |
             |  resM , | | | u     o nabla | Dacc  | o nabla | v  |
             |          \ \ \             /       /         /     |
              \                                                  /

            */

            elemat(tvi  ,tui  ) -= inertia_and_conv[0]*derxy_(0,vi);
            elemat(tvi  ,tuip ) -= inertia_and_conv[0]*derxy_(1,vi);

            elemat(tvip ,tui  ) -= inertia_and_conv[1]*derxy_(0,vi);
            elemat(tvip ,tuip ) -= inertia_and_conv[1]*derxy_(1,vi);
          } // vi
        } // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          double temp[2];
          temp[0]=fac_afgdt_tauM_tauM*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi));
          temp[1]=fac_afgdt_tauM_tauM*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi));

          double rowtemp[2][2];

          rowtemp[0][0]=resM_(0)*temp[0];
          rowtemp[0][1]=resM_(0)*temp[1];

          rowtemp[1][0]=resM_(1)*temp[0];
          rowtemp[1][1]=resM_(1)*temp[1];

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            /*
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / /            \   n+af \         \     |
             |  resM , | | | Dacc o nabla | u      | o nabla | v  |
             |          \ \ \            /        /         /     |
              \                                                  /

            */

            elemat(tvi  ,tui  ) -= funct_(ui)*rowtemp[0][0];
            elemat(tvi  ,tuip ) -= funct_(ui)*rowtemp[0][1];

            elemat(tvip ,tui  ) -= funct_(ui)*rowtemp[1][0];
            elemat(tvip ,tuip ) -= funct_(ui)*rowtemp[1][1];
          } // ui
        } // vi


        const double fac_gdt_tauM_tauM       =fac*gamma*dt*tauM*tauM;
        const double fac_gdt_tauM_tauM_resM_x=fac_gdt_tauM_tauM*resM_(0);
        const double fac_gdt_tauM_tauM_resM_y=fac_gdt_tauM_tauM*resM_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp =3*ui+2;

          double coltemp[2][2];

          coltemp[0][0]=fac_gdt_tauM_tauM_resM_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_tauM_tauM_resM_x*derxy_(1,ui);
          coltemp[1][0]=fac_gdt_tauM_tauM_resM_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_tauM_tauM_resM_y*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: - gamma * dt * tauM * tauM (rescaled)

              /                               \
             |          /                \     |
             |  resM , | nabla Dp o nabla | v  |
             |          \                /     |
              \                               /

            */

            elemat(tvi  ,tuipp) -= coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi);
            elemat(tvip ,tuipp) -= coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi);

          } // vi
        } // ui


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_nu_afgdt_tauM_tauM=fac*visceff*afgdt*tauM*tauM;

          double temp[2];

          temp[0]=fac_nu_afgdt_tauM_tauM*resM_(0);
          temp[1]=fac_nu_afgdt_tauM_tauM*resM_(1);


          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            double rowtemp[2][2];

            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);

            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);


            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui   =3*ui;
              const int tuip  =tui+1;

              /*
                   factor: + 2.0 * visc * alphaF * gamma * dt * tauM * tauM

                    /                                                \
                   |          / /             /    \  \         \     |
                   |  resM , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /
              */

              elemat(tvi  ,tui  ) += viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
              elemat(tvi  ,tuip ) += derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

              elemat(tvip ,tui  ) += viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
              elemat(tvip ,tuip ) += derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
            } // ui
          } //vi
        } // end hoel
      } // end if reynolds stab


      //---------------------------------------------------------------
      /*
             GALERKIN PART, INERTIA, CONVECTION AND VISCOUS TERMS
	                  QUASISTATIC FORMULATION

        ---------------------------------------------------------------

          inertia term (intermediate) + convection (intermediate)

                /          \                   /                          \
               |            |                 |  / n+af       \            |
      +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
               |            |                 |  \            /            |
                \          /                   \                          /


	  viscous term (intermediate), factor: +2*nu*alphaF*gamma*dt

                                    /                          \
                                   |       /    \         / \   |
             +2*nu*alphaF*gamma*dt |  eps | Dacc | , eps | v |  |
                                   |       \    /         \ /   |
                                    \                          /

|	  convection (intermediate)
|
N                                /                            \
E                               |  /            \   n+af       |
W              +alphaF*gamma*dt | | Dacc o nabla | u      , v  |
T                               |  \            /              |
O                                \                            /
N
      */
      //---------------------------------------------------------------


      /*---------------------------------------------------------------

                     SUPG PART, INERTIA AND CONVECTION TERMS
                  REYNOLDS PART, SUPG-TESTFUNCTION TYPE TERMS
                      QUASISTATIC FORMULATION (IF ACTIVE)

        ---------------------------------------------------------------

          inertia and convection, factor: +alphaM*tauM

                             /                                        \                    -+
                            |          / / n+af  ~n+af \         \     |                    |     c
               +alphaM*tauM*|  Dacc , | | c    + u      | o nabla | v  |                    |     o
                            |          \ \             /         /     |                    | i   n
                             \                                        /                     | n   v
                                                                                            | e a e
                                                                                            | r n c
                             /                                                          \   | t d t
                            |   / n+af        \          / / n+af  ~n+af \         \     |  | i   i
      +alphaF*gamma*dt*tauM*|  | c     o nabla | Dacc , | | c    + u      | o nabla | v  |  | a   o
                            |   \             /          \ \             /         /     |  |     n
                             \                                                          /  -+


                                                                                              p
                             /                                            \                -+ r
                            |              / / n+af  ~n+af \         \     |                | e
            +tauM*gamma*dt* |  nabla Dp , | | c    + u      | o nabla | v  |                | s
                            |              \ \             /         /     |                | s
                             \                                            /                -+ u
                                                                                              r
                                                                                              e

                                                                                              d
                                                                                              i
                             /                                                           \ -+ f
                            |                 /     \    /  / n+af  ~n+af \         \     | | f
   -nu*alphaF*gamma*dt*tauM*|  2*nabla o eps | Dacc  |, |  | c    + u      | o nabla | v  | | u
                            |                 \     /    \  \             /         /     | | s
                             \                                                           / -+ i
                                                                                              o
                                                                                              n


|         linearised convective term in residual
|
N                            /                                                           \
E                           |    /            \   n+af    / / n+af  ~n+af \         \     |
W     +alphaF*gamma*dt*tauM |   | Dacc o nabla | u     , | | c    + u      | o nabla | v  |
T                           |    \            /           \ \             /         /     |
O                            \                                                           /
N

|	  linearisation of testfunction
|
N                            /                            \
E                           |   n+af    /            \     |
W     +alphaF*gamma*dt*tauM*|  r     , | Dacc o nabla | v  |
T                           |   M       \            /     |
O                            \                            /
N

      */
      //---------------------------------------------------------------


      //---------------------------------------------------------------
      /*
	           LEAST SQUARES CONTINUITY STABILISATION PART,
	              QUASISTATIC FORMULATION (IF ACTIVE)

        ---------------------------------------------------------------

          factor: +gamma*dt*tauC

                         /                          \
                        |                            |
                        | nabla o Dacc  , nabla o v  |
                        |                            |
                         \                          /
      */


      const double fac_afgdt         = fac*afgdt;
      const double fac_visceff_afgdt = fac_afgdt*visceff;
      const double fac_gamma_dt      = fac*gamma*dt;
      const double fac_alphaM        = fac*alphaM;

      const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tui  =3*ui;
        const int tuip =tui+1;

        /* GALERKIN inertia term (intermediate) + convection (intermediate) */
        const double inertia_and_conv_ui
          = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);

        /* viscous term (intermediate), diagonal parts */
        const double fac_visceff_afgdt_derxy0_ui=fac_visceff_afgdt*derxy_(0,ui);
        const double fac_visceff_afgdt_derxy1_ui=fac_visceff_afgdt*derxy_(1,ui);

        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          /* add:                                                             */
          /* GALERKIN inertia term (intermediate) + convection (intermediate) */
          /* SUPG stabilisation --- inertia and convection                    */
          /* viscous term (intermediate), diagonal parts                      */
          const double sum =
            inertia_and_conv_ui*(funct_(vi)+conv_c_plus_svel_af_(vi))
            +
            fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
            +
            fac_visceff_afgdt_derxy1_ui*derxy_(1,vi);

          elemat(tvi  ,tui  ) += sum+(fac_visceff_afgdt_derxy0_ui+fac_gamma_dt_tauC_derxy_x_ui)*derxy_(0,vi);
          elemat(tvi  ,tuip ) +=      fac_visceff_afgdt_derxy0_ui*derxy_(1,vi)+fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi);
          elemat(tvip ,tui  ) +=      fac_visceff_afgdt_derxy1_ui*derxy_(0,vi)+fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi);
          elemat(tvip ,tuip ) += sum+(fac_visceff_afgdt_derxy1_ui+fac_gamma_dt_tauC_derxy_y_ui)*derxy_(1,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp=3*ui+2;

        const double fac_gamma_dt_derxy_0_ui = fac_gamma_dt*derxy_(0,ui);
        const double fac_gamma_dt_derxy_1_ui = fac_gamma_dt*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          int tvi =vi*3;

          /* SUPG stabilisation --- pressure     */
          /* factor: +tauM, rescaled by gamma*dt */

          elemat(tvi++,tuipp) += fac_gamma_dt_derxy_0_ui*conv_c_plus_svel_af_(vi);
          elemat(tvi  ,tuipp) += fac_gamma_dt_derxy_1_ui*conv_c_plus_svel_af_(vi);

        } // vi
      } // ui

      if (higher_order_ele && newton!=INPAR::FLUID::minimal)
      {
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui  =ui*3 ;
          const int tuip =tui+1;

          /* SUPG stabilisation --- diffusion   */
          /* factor: -nu*alphaF*gamma*dt*tauM   */
          const double fac_visceff_afgdt_viscs2_0_ui=fac_visceff_afgdt*viscs2_(0,ui);
          const double fac_visceff_afgdt_viscs2_1_ui=fac_visceff_afgdt*viscs2_(1,ui);
          const double fac_visceff_afgdt_derxy2_2_ui=fac_visceff_afgdt*derxy2_(2,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int tvi  =vi*3 ;
            int tvip =tvi+1;

            elemat(tvi  ,tui  ) -= fac_visceff_afgdt_viscs2_0_ui*conv_c_plus_svel_af_(vi);
            elemat(tvi  ,tuip ) -= fac_visceff_afgdt_derxy2_2_ui*conv_c_plus_svel_af_(vi);
            elemat(tvip ,tui  ) -= fac_visceff_afgdt_derxy2_2_ui*conv_c_plus_svel_af_(vi);
            elemat(tvip ,tuip ) -= fac_visceff_afgdt_viscs2_1_ui*conv_c_plus_svel_af_(vi);
          } //end vi
        } // end ui
      }// end higher_order_ele and linearisation of viscous term

      //---------------------------------------------------------------
      //
      //                  GALERKIN AND SUPG PART
      //    REACTIVE TYPE LINEARISATIONS, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if (newton==INPAR::FLUID::Newton)
      {
        double temp[2][2];
        double supg_active_tauMresM[2];

        /* for linearisation of supg testfunction */
        supg_active_tauMresM[0]=supg_active_tauM*resM_(0);
        supg_active_tauMresM[1]=supg_active_tauM*resM_(1);

        // loop rows
        for (int vi=0; vi<iel; ++vi)
        {
          int tvi   =vi*3;
          int tvip  =tvi+1;

          /*  add linearised convective term in residual (supg),
              linearisation of testfunction (supg)
              and linearised Galerkin term                */
          temp[0][0]=fac_afgdt*(supg_active_tauMresM[0]*derxy_(0,vi)+vderxyaf_(0,0)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[0][1]=fac_afgdt*(supg_active_tauMresM[0]*derxy_(1,vi)+vderxyaf_(0,1)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[1][0]=fac_afgdt*(supg_active_tauMresM[1]*derxy_(0,vi)+vderxyaf_(1,0)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[1][1]=fac_afgdt*(supg_active_tauMresM[1]*derxy_(1,vi)+vderxyaf_(1,1)*(conv_c_plus_svel_af_(vi)+funct_(vi)));

          // loop columns
          for (int ui=0; ui<iel; ++ui)
          {
            int tui=3*ui;

            elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
            elemat(tvip ,tui++) += temp[1][0]*funct_(ui);
            elemat(tvi  ,tui  ) += temp[0][1]*funct_(ui);
            elemat(tvip ,tui  ) += temp[1][1]*funct_(ui);
          } // ui
        } // vi
      } // end newton

      //---------------------------------------------------------------
      //
      //      GALERKIN PART, CONTINUITY AND PRESSURE PART
      //                QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------

      for (int vi=0; vi<iel; ++vi)  // loop rows
      {
        const int tvi    =3*vi;
        const int tvip   =tvi+1;

        const double fac_gamma_dt_derxy_0_vi = fac_gamma_dt*derxy_(0,vi);
        const double fac_gamma_dt_derxy_1_vi = fac_gamma_dt*derxy_(1,vi);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp=3*ui+2;

          /* GALERKIN pressure (implicit, rescaled to keep symmetry) */

          /*  factor: -1, rescaled by gamma*dt

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
          */

          elemat(tvi  ,tuipp) -= fac_gamma_dt_derxy_0_vi*funct_(ui);
          elemat(tvip ,tuipp) -= fac_gamma_dt_derxy_1_vi*funct_(ui);

          /* GALERKIN continuity equation (implicit, transposed of above equation) */

          /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
          */

          elemat(tuipp,tvi  ) += fac_gamma_dt_derxy_0_vi*funct_(ui);
          elemat(tuipp,tvip ) += fac_gamma_dt_derxy_1_vi*funct_(ui);
        } // ui
      } // vi

      //---------------------------------------------------------------
      //
      //             PSPG PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if(pspg == INPAR::FLUID::pstab_use_pspg)
      {
        const double fac_tauMp                   = fac*tauMp;
        const double fac_alphaM_tauMp            = fac_tauMp*alphaM;
        const double fac_gamma_dt_tauMp          = fac_tauMp*gamma*dt;
        const double fac_afgdt_tauMp             = fac_tauMp*afgdt;

        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_tauMp
            =
            fac*visceff*afgdt*tauMp;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =ui*3;
            const int tuip =tui+1;

            /* pressure stabilisation --- diffusion  */

            /* factor: -nu*alphaF*gamma*dt*tauMp

                    /                                  \
                   |                 /    \             |
                   |  2*nabla o eps | Dacc | , nabla q  |
                   |                 \    /             |
                    \                                  /
            */

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
            const double fac_tauMp_inertia_and_conv
              =
              fac_alphaM_tauMp*funct_(ui)+fac_afgdt_tauMp*conv_c_af_(ui);

            const double pspg_diffusion_inertia_convect_0_ui
              =
              fac_visceff_afgdt_tauMp*viscs2_(0,ui)-fac_tauMp_inertia_and_conv;
            const double pspg_diffusion_inertia_convect_1_ui
              =
              fac_visceff_afgdt_tauMp*viscs2_(1,ui)-fac_tauMp_inertia_and_conv;

            const double fac_visceff_afgdt_tauMp_derxy2_2_ui=fac_visceff_afgdt_tauMp*derxy2_(2,ui);


            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvipp=vi*3+2;

              elemat(tvipp,tui  ) -=
                pspg_diffusion_inertia_convect_0_ui*derxy_(0,vi)
                +
                fac_visceff_afgdt_tauMp_derxy2_2_ui*derxy_(1,vi);
              elemat(tvipp,tuip ) -=
                fac_visceff_afgdt_tauMp_derxy2_2_ui*derxy_(0,vi)
                +
                pspg_diffusion_inertia_convect_1_ui*derxy_(1,vi);
            } // vi
          } // ui
        } // this is a higher order element and linearisation is not minimal
        else
        { // either this ain't a higher order element or a
          // linearisation of the viscous term is not necessary
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =ui*3 ;
            const int tuip =tui+1;

            const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvipp=vi*3+2;

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

              elemat(tvipp,tui ) += fac_tauMp_inertia_and_conv*derxy_(0,vi) ;
              elemat(tvipp,tuip) += fac_tauMp_inertia_and_conv*derxy_(1,vi) ;
            } // vi
          } // ui
        } // no linearisation of viscous part of residual is
          // performed for pspg stabilisation cause either this
          // ain't a higher order element or a linearisation of
          // the viscous term is not necessary

        if (newton==INPAR::FLUID::Newton)
        {
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int tvipp = vi*3+2;
            double v1 = derxy_(0,vi)*vderxyaf_(0,0) + derxy_(1,vi)*vderxyaf_(1,0);
            double v2 = derxy_(0,vi)*vderxyaf_(0,1) + derxy_(1,vi)*vderxyaf_(1,1)
;
            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const double fac_afgdt_tauMp_funct_ui = fac_afgdt_tauMp*funct_(ui);
              int tui = ui*3;

              /* pressure stabilisation --- convection */

              /*  factor: +alphaF*gamma*dt*tauMp

                       /                                  \
                      |  /            \   n+af             |
                      | | Dacc o nabla | u      , nabla q  |
                      |  \            /                    |
                       \                                  /
              */

              elemat(tvipp, tui  ) += fac_afgdt_tauMp_funct_ui*v1;
              elemat(tvipp, tui+1) += fac_afgdt_tauMp_funct_ui*v2;
            } // ui
          } // vi
        } // end newton

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp=ui*3+2;

          const double fac_gamma_dt_tauMp_derxy_0_ui=fac_gamma_dt_tauMp*derxy_(0,ui);
          const double fac_gamma_dt_tauMp_derxy_1_ui=fac_gamma_dt_tauMp*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            /* pressure stabilisation --- rescaled pressure   */

            /* factor: +tauMp, rescaled by gamma*dt

                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
            */

            elemat(vi*3+2,tuipp) +=
              fac_gamma_dt_tauMp_derxy_0_ui*derxy_(0,vi)
              +
              fac_gamma_dt_tauMp_derxy_1_ui*derxy_(1,vi);

          } // vi
        } // ui
      } // end pspg

      //---------------------------------------------------------------
      //
      //      VISCOUS STABILISATION PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if (higher_order_ele)
      {
        if(vstab == INPAR::FLUID::viscous_stab_gls
           ||
           vstab == INPAR::FLUID::viscous_stab_usfem)
        {
          const double fac_visc_tauMp_gamma_dt      = vstabfac*fac*visc*tauMp*gamma*dt;
          const double fac_visc_afgdt_tauMp         = vstabfac*fac*visc*afgdt*tauMp;
          const double fac_visc_alphaM_tauMp        = vstabfac*fac*visc*alphaM*tauMp;
          const double fac_visceff_visc_afgdt_tauMp = vstabfac*fac*visceff*visc*afgdt*tauMp;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui    = ui*3;
            const int tuip   = tui+1;
            const int tuipp  = tui+2;

            const double acc_conv=(fac_visc_alphaM_tauMp*funct_(ui)
                                   +
                                   fac_visc_afgdt_tauMp*conv_c_af_(ui));

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvi   = vi*3;
              const int tvip  = tvi+1;

              /* viscous stabilisation --- inertia     */

              /* factor: +(-)alphaM*tauMp*nu

                    /                      \
                   |                        |
                   |  Dacc , 2*div eps (v)  |
                   |                        |
                    \                      /
              */
              /* viscous stabilisation --- convection */

              /*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                       /                                    \
                      |  / n+af       \                      |
                      | | c    o nabla | Dacc, 2*div eps (v) |
                      |  \            /                      |
                       \                                    /

              */

              elemat(tvi  ,tui  ) += acc_conv*viscs2_(0,vi);
              elemat(tvi  ,tuip ) += acc_conv*derxy2_(2,vi);
              elemat(tvip ,tui  ) += acc_conv*derxy2_(2,vi);
              elemat(tvip ,tuip ) += acc_conv*viscs2_(1,vi);

              /* viscous stabilisation --- diffusion  */

              /* factor: -(+)nu*nu*alphaF*gamma*dt*tauMp

                    /                                       \
                   |                 /    \                  |
                   |  2*nabla o eps | Dacc | , 2*div eps (v) |
                   |                 \    /                  |
                    \                                       /
              */
              elemat(tvi  ,tui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy2_(2,ui)*derxy2_(2,vi));
              elemat(tvi  ,tuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(2,ui)
                                      +
                                      derxy2_(2,vi)*viscs2_(1,ui));
              elemat(tvip ,tui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(2,vi)
                                      +
                                      derxy2_(2,ui)*viscs2_(1,vi));
              elemat(tvip ,tuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(2,ui)*derxy2_(2,vi)
                                      +
                                      viscs2_(1,ui)*viscs2_(1,vi));

              /* viscous stabilisation --- pressure   */

              /* factor: +(-)tauMp*nu, rescaled by gamma*dt

                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
                */
              elemat(tvi  ,tuipp ) += fac_visc_tauMp_gamma_dt*
                                   (derxy_(0,ui)*viscs2_(0,vi)
                                    +
                                    derxy_(1,ui)*derxy2_(2,vi));
              elemat(tvip ,tuipp ) += fac_visc_tauMp_gamma_dt*
                                   (derxy_(0,ui)*derxy2_(2,vi)
                                    +
                                    derxy_(1,ui)*viscs2_(1,vi));

            } // vi
          } // ui

          if (newton==INPAR::FLUID::Newton)
          {
            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui    = ui*3;
              const int tuip   = tui+1;

              const double fac_visc_afgdt_tauMp_funct_ui=fac_visc_afgdt_tauMp*funct_(ui);

              for (int vi=0; vi<iel; ++vi)  // loop rows
              {
                const int tvi   = vi*3;
                const int tvip  = tvi+1;

                /* viscous stabilisation --- convection */

                /*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /


                */
                elemat(tvi  ,tui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,0)
                                        +
                                        derxy2_(2,vi)*vderxyaf_(1,0));
                elemat(tvi  ,tuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,1)
                                        +
                                        derxy2_(2,vi)*vderxyaf_(1,1));
                elemat(tvip ,tui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(2,vi)*vderxyaf_(0,0)
                                        +
                                        viscs2_(1,vi)*vderxyaf_(1,0));
                elemat(tvip ,tuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(2,vi)*vderxyaf_(0,1)
                                        +
                                        viscs2_(1,vi)*vderxyaf_(1,1));
              } // vi
            } // ui
          } // end newton
        } // endif (a)gls
      } // end hoel

      //---------------------------------------------------------------
      //
      //               QUASISTATIC STABILISATION PART
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == INPAR::FLUID::cross_stress_stab)
      {
        const double fac_afgdt_tauM = fac*afgdt*tauM;

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double fac_afgdt_tauM_conv_resM_ui = fac_afgdt_tauM*conv_resM_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int tvi =3*vi;
            const double fac_afgdt_tauM_conv_resM_ui_funct_vi = fac_afgdt_tauM_conv_resM_ui*funct_(vi);

            /*  factor:

                -alphaF*gamma*dt*tauM

                          /                          \
                         |  /            \            |
                         | | resM o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
            */
            elemat(tvi++,tui  ) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
            elemat(tvi  ,tuip ) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
          } // vi
        } // ui

        const double fac_alphaM_tauM = fac*alphaM*tauM;

        double  am_nabla_u_afgdt_nabla_u_nabla_u[2][2];

        am_nabla_u_afgdt_nabla_u_nabla_u[0][0]=
          fac_alphaM_tauM*vderxyaf_(0,0)
          +
          fac_afgdt_tauM*(vderxyaf_(0,0)*vderxyaf_(0,0)
                          +
                          vderxyaf_(0,1)*vderxyaf_(1,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][1]=
          fac_alphaM_tauM*vderxyaf_(0,1)
          +
          fac_afgdt_tauM*(vderxyaf_(0,0)*vderxyaf_(0,1)
                          +
                          vderxyaf_(0,1)*vderxyaf_(1,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][0]=
          fac_alphaM_tauM*vderxyaf_(1,0)
          +
          fac_afgdt_tauM*(vderxyaf_(1,0)*vderxyaf_(0,0)
                          +
                          vderxyaf_(1,1)*vderxyaf_(1,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][1]=
          fac_alphaM_tauM*vderxyaf_(1,1)
          +
          fac_afgdt_tauM*(vderxyaf_(1,0)*vderxyaf_(0,1)
                          +
                          vderxyaf_(1,1)*vderxyaf_(1,1));

        double nabla_u[2][2];
        nabla_u[0][0]=fac_afgdt_tauM*vderxyaf_(0,0);
        nabla_u[0][1]=fac_afgdt_tauM*vderxyaf_(0,1);

        nabla_u[1][0]=fac_afgdt_tauM*vderxyaf_(1,0);
        nabla_u[1][1]=fac_afgdt_tauM*vderxyaf_(1,1);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_nabla_ui = velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double coltemp[2][2];

          coltemp[0][0]=am_nabla_u_afgdt_nabla_u_nabla_u[0][0]*funct_(ui)+nabla_u[0][0]*u_nabla_ui;
          coltemp[0][1]=am_nabla_u_afgdt_nabla_u_nabla_u[0][1]*funct_(ui)+nabla_u[0][1]*u_nabla_ui;

          coltemp[1][0]=am_nabla_u_afgdt_nabla_u_nabla_u[1][0]*funct_(ui)+nabla_u[1][0]*u_nabla_ui;
          coltemp[1][1]=am_nabla_u_afgdt_nabla_u_nabla_u[1][1]*funct_(ui)+nabla_u[1][1]*u_nabla_ui;

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;


            /*  factor:

                -alphaM*tauM

                          /                           \
                         |  /            \   n+af      |
                         | | Dacc o nabla | u     , v  |
                         |  \            /             |
                          \                           /
            */

            /*  factor:

                -alphaF*gamma*dt*tauM

                          /                                                \
                         |  / / /            \   n+af \         \   n+af    |
                         | | | | Dacc o nabla | u      | o nabla | u   , v  |
                         |  \ \ \            /        /         /           |
                          \                                                /
            */

            /*  factor:

                -alphaF*gamma*dt*tauM

                          /                                                 \
                         |  / / / n+af        \       \         \   n+af     |
                         | | | | u     o nabla | Dacc  | o nabla | u    , v  |
                         |  \ \ \             /       /         /            |
                          \                                                 /
            */

            elemat(tvi  ,tui  ) -= funct_(vi)*coltemp[0][0];
            elemat(tvi  ,tuip ) -= funct_(vi)*coltemp[0][1];

            elemat(tvip ,tui  ) -= funct_(vi)*coltemp[1][0];
            elemat(tvip ,tuip ) -= funct_(vi)*coltemp[1][1];
          } // vi
        } // ui

        const double fac_gdt_tauM = fac*gamma*dt*tauM;
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp =3*ui+2;

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*  factor:

               -gamma*dt*tauM (rescaled for consistency)

                          /                               \
                         |  /                \   n+af      |
                         | | nabla Dp o nabla | u     , v  |
                         |  \                /             |
                          \                               /
            */
            elemat(tvi  ,tuipp) -= fac_gdt_tauM*funct_(vi)*(vderxyaf_(0,0)*derxy_(0,ui)+vderxyaf_(0,1)*derxy_(1,ui));
            elemat(tvip ,tuipp) -= fac_gdt_tauM*funct_(vi)*(vderxyaf_(1,0)*derxy_(0,ui)+vderxyaf_(1,1)*derxy_(1,ui));
          } // vi
        } // ui

        if (higher_order_ele)
        {
          const double fac_visceff_afgdt_tauM = fac_afgdt_tauM*visceff;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            double coltemp[2][2];

            coltemp[0][0]=fac_visceff_afgdt_tauM*(viscs2_(0,ui)+derxy2_(2,ui)*vderxyaf_(0,1));
            coltemp[0][1]=fac_visceff_afgdt_tauM*(viscs2_(1,ui)+derxy2_(2,ui)*vderxyaf_(0,0));

            coltemp[1][0]=fac_visceff_afgdt_tauM*(viscs2_(0,ui)+derxy2_(2,ui)*vderxyaf_(1,1));
            coltemp[1][1]=fac_visceff_afgdt_tauM*(viscs2_(1,ui)+derxy2_(2,ui)*vderxyaf_(1,0));

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /*  factor:

                  +alphaF*gamma*dt*tauM

                          /                                               \
                         |  / /             /    \ \         \   n+af      |
                         | | | nabla o eps | Dacc | | o nabla | u     , v  |
                         |  \ \             \    / /         /             |
                          \                                               /
              */
              elemat(tvi  ,tui  ) += coltemp[0][0]*funct_(vi);
              elemat(tvi  ,tuip ) += coltemp[0][1]*funct_(vi);

              elemat(tvip ,tui  ) += coltemp[1][0]*funct_(vi);
              elemat(tvip ,tuip ) += coltemp[1][1]*funct_(vi);

            } // vi
          } // ui
        } // hoel
      } // end cross

    } // end if compute_elemat

    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //
    //          RIGHT HAND SIDE, QUASISTATIC SUBGRID SCALES
    //
    //---------------------------------------------------------------
    //---------------------------------------------------------------

    /* inertia, convective and dead load terms -- all tested
       against shapefunctions, as well as cross terms            */
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

    double fac_inertia_conv_and_bodyforce[2];
    fac_inertia_conv_and_bodyforce[0] = fac*(accintam_(0)+convaf_old_(0)-bodyforceaf_(0));
    fac_inertia_conv_and_bodyforce[1] = fac*(accintam_(1)+convaf_old_(1)-bodyforceaf_(1));

    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs
       ||
       cross == INPAR::FLUID::cross_stress_stab)
    {
      const double fac_tauM = fac*tauM;

      /* factor: +tauM

                  /                            \
                 |                    n+af      |
                 |  ( resM o nabla ) u    ,  v  |
                 |                    (i)       |
                  \                            /
      */

      fac_inertia_conv_and_bodyforce[0]-=
        fac_tauM*
        (resM_(0)*vderxyaf_(0,0)+
         resM_(1)*vderxyaf_(0,1));
      fac_inertia_conv_and_bodyforce[1]-=
        fac_tauM*
        (resM_(0)*vderxyaf_(1,0)+
         resM_(1)*vderxyaf_(1,1));
    }

    /*
      pressure and viscous term combined in viscous_and_pres
      cross and reynolds stabilisation are combined with the
      same testfunctions (of derivative type)
    */

    // continuity stabilisation adds a small-scale pressure

    /*
      factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /

    */

    /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
    */
    const double fac_prenp   = fac*prenp_-fac*tauC*divunp_;

    /*
      factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
    */

    const double visceff_fac = visceff*fac;

    double viscous_and_pres[4];
    viscous_and_pres[0]=visceff_fac*vderxyaf_(0,0)*2.0-fac_prenp;
    viscous_and_pres[1]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_and_pres[2]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_and_pres[3]=visceff_fac*vderxyaf_(1,1)*2.0-fac_prenp;

    if(reynolds == INPAR::FLUID::reynolds_stress_stab_only_rhs
       ||
       reynolds == INPAR::FLUID::reynolds_stress_stab)
    {

      /* factor: -tauM*tauM

                  /                             \
                 |                               |
                 |  resM   , ( resM o nabla ) v  |
                 |                               |
                  \                             /
      */
      const double fac_tauM_tauM        = fac*tauM*tauM;
      const double fac_tauM_tauM_resM_0 = fac_tauM_tauM*resM_(0);
      const double fac_tauM_tauM_resM_1 = fac_tauM_tauM*resM_(1);

      viscous_and_pres[0]-=fac_tauM_tauM_resM_0*resM_(0);
      viscous_and_pres[1]-=fac_tauM_tauM_resM_0*resM_(1);
      viscous_and_pres[2]-=fac_tauM_tauM_resM_0*resM_(1);
      viscous_and_pres[3]-=fac_tauM_tauM_resM_1*resM_(1);
    }

    /* continuity equation, factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
    */
    const double fac_divunp  = fac*divunp_;

    for (int vi=0; vi<iel; ++vi) // loop rows
    {
      const int tvi  =3*vi;
      const int tvip =tvi+1;
      const int tvipp=tvi+2;
      /* inertia, convective and dead load, cross terms fith
         funct                                              */
      /* viscous, pressure, reynolds, cstab terms with
         derxy                                              */

      elevec(tvi  ) -=
        fac_inertia_conv_and_bodyforce[0]*funct_(vi)
        +
        derxy_(0,vi)*viscous_and_pres[0]
        +
        derxy_(1,vi)*viscous_and_pres[1];
      elevec(tvip ) -=
        fac_inertia_conv_and_bodyforce[1]*funct_(vi)
        +
        derxy_(0,vi)*viscous_and_pres[2]
        +
        derxy_(1,vi)*viscous_and_pres[3];

      /* continuity equation */
      elevec(tvipp) -= fac_divunp*funct_(vi);
    }

    if(pspg == INPAR::FLUID::pstab_use_pspg)
    {
      /*
            pressure stabilisation

            factor: +tauMp

                  /                 \
                 |    n+af           |
                 |  r     , nabla q  |
                 |   M               |
                  \                 /

      */
      const double fac_tauMp = fac*tauMp;

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        elevec(3*vi+2)-=fac_tauMp*(resM_(0)*derxy_(0,vi)+resM_(1)*derxy_(1,vi));
      }
    } // end pspg

    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      const double fac_tauM = fac*supg_active_tauM;

      for (int vi=0; vi<iel; ++vi) // loop rows  (test functions)
      {
        int tvi=3*vi;

        const double fac_tauM_conv_c_af_vi = fac_tauM*conv_c_af_(vi);
        /*
          factor: +tauM

          SUPG stabilisation


                  /                             \
                 |   n+af    / n+af        \     |
                 |  r     , | c     o nabla | v  |
                 |   M       \             /     |
                  \                             /
        */

        elevec(tvi++) -= fac_tauM_conv_c_af_vi*resM_(0);
        elevec(tvi  ) -= fac_tauM_conv_c_af_vi*resM_(1);
      } // end loop rows
    } // end supg

    if (higher_order_ele)
    {
      if(vstab != INPAR::FLUID::viscous_stab_none && higher_order_ele)
      {
        const double fac_visc_tauMp = vstabfac * fac*visc*tauMp;

        for (int vi=0; vi<iel; ++vi) // loop rows
        {
          int tvi=3*vi;
          /*
              factor: -(+)tauMp*nu

              viscous stabilisation --- inertia


                 /                      \
                |   n+af                 |
                |  r    , 2*div eps (v)  |
                |   M                    |
                 \                      /

          */
          elevec(tvi++) -= fac_visc_tauMp*
                           (resM_(0)*viscs2_(0,vi)
                            +
                            resM_(1)*derxy2_(2,vi)) ;
          elevec(tvi  ) -= fac_visc_tauMp*
                           (resM_(0)*derxy2_(2,vi)
                            +
                            resM_(1)*viscs2_(1,vi)) ;
        } // end loop rows vi
      } // endif (a)gls
    } // hoel

  } // end loop iquad

  return;
} // Sysmat_adv_qs


/*----------------------------------------------------------------------*
 |  calculate system matrix for a generalised alpha time integration    |
 |  advective version using time dependent subgrid scales               |
 |                            (public)                      gammi 06/07 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_adv_td_2D(
  Fluid*                                       ele,
  std::vector<Epetra_SerialDenseVector>&        myknots,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&    elemat,
  LINALG::Matrix<(nsd_+1)*iel,    1>&     elevec,
  const LINALG::Matrix<nsd_,iel>&         edispnp,
  const LINALG::Matrix<nsd_,iel>&         egridvaf,
  const LINALG::Matrix<nsd_,iel>&         evelnp,
  const LINALG::Matrix<iel,1>&            eprenp,
  const LINALG::Matrix<nsd_,iel>&         eaccam,
  const LINALG::Matrix<nsd_,iel>&         evelaf,
  Teuchos::RCP<const MAT::Material>       material,
  const double                            alphaM,
  const double                            alphaF,
  const double                            gamma,
  const double                            dt,
  const double                            time,
  const enum INPAR::FLUID::LinearisationAction  newton,
  const bool                              higher_order_ele,
  const enum INPAR::FLUID::Transient      inertia         ,
  const enum INPAR::FLUID::PSPG           pspg            ,
  const enum INPAR::FLUID::SUPG           supg            ,
  const enum INPAR::FLUID::VStab          vstab           ,
  const enum INPAR::FLUID::CStab          cstab           ,
  const enum INPAR::FLUID::CrossStress    cross           ,
  const enum INPAR::FLUID::ReynoldsStress reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha              whichtau,
  double&                                 visceff,
  const bool                              compute_elemat
  )
{

  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt
  //
  const double timealphaF = time-(1-alphaF)*dt;

  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }


  //------------------------------------------------------------------
  //                    SET ALL ELEMENT DATA
  // o including element geometry (node coordinates)
  // o including dead loads in nodes
  // o including hk, mk, element area
  // o including material viscosity, effective viscosity by
  //   Non-Newtonian fluids
  //------------------------------------------------------------------

  double hk   = 0.0;
  double mk   = 0.0;
  double visc = 0.0;

  SetElementData_2D(ele            ,
                 edispnp        ,
                 evelaf         ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 visceff        );

#if 1

  {
    // gaussian points
    //const DRT::UTILS::IntegrationPoints2D intpoints_onepoint(integrationrule_stabili);
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_onepoint(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);


    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      0           ,
      intpoints_onepoint       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_time_dependent,gamma,dt,hk,mk,visceff);
  }
#endif

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  //const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // remember whether the subscale quantities have been allocated an set to zero.
  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero

//  ele->ActivateTDS(intpoints.IP().nquad, nsd_);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {

    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);

#if 0
    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_time_dependent,gamma,dt,hk,mk,visceff);
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

    const double tauM   = tau_(0);

    if(cstab == INPAR::FLUID::continuity_stab_none)
    {
      tau_(2)=0.0;
    }

    const double tauC   = tau_(2);

    double supg_active;
    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      supg_active=1.0;
    }
    else
    {
      supg_active=0.0;
    }

    // update estimates for the subscale quantities
    const double facMtau = 1./(alphaM*tauM+afgdt);
//    const double fac1=(alphaM*tauM+gamma*dt*(alphaF-1.0))*facMtau;
//    const double fac2=(dt*tauM*(alphaM-gamma))*facMtau;
//    const double fac3=(gamma*dt*tauM)*facMtau;

    /*-------------------------------------------------------------------*
     *                                                                   *
     *                  update of SUBSCALE VELOCITY                      *
     *             and update of intermediate quantities                 *
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

     compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

  */

//    for(int rr=0;rr<nsd_;++rr)
//    {
//      ele->UpdateSvelnpInOneDirection(
//          fac1       ,
//          fac2       ,
//          fac3       ,
//          resM_(rr)  ,
//          alphaF     ,
//          rr         ,
//          iquad      ,
//          svelnp_(rr),
//          svelaf_(rr));
//    }

    /* the intermediate value of subscale acceleration is not needed to be
     * computed anymore --- we use the governing ODE to replace it ....

             ~ n+am    alphaM     / ~n+1   ~n \    gamma - alphaM    ~ n
            acc     = -------- * |  u    - u   | + -------------- * acc
               (i)    gamma*dt    \  (i)      /         gamma

    */

    // prepare possible modification of convective linearisation for
    // combined reynolds/supg test function
    for(int nn=0;nn<iel;++nn)
    {
      conv_c_plus_svel_af_(nn)=conv_c_af_(nn)*supg_active;
    }

    /*
        This is the operator

                  /~n+af         \
                 | u      o nabla |
                  \   (i)        /

        required for the cross/reynolds stress linearisation

    */
    if(cross    == INPAR::FLUID::cross_stress_stab
       ||
       reynolds == INPAR::FLUID::reynolds_stress_stab)
    {
      for (int rr=0;rr<iel;++rr)
      {
        conv_subaf_(rr)=svelaf_(0)*derxy_(0,rr)+svelaf_(1)*derxy_(1,rr);
      }

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {
        /* get modified convective linearisation (n+alpha_F,i) at
           integration point takes care of half of the linearisation

                                   +-----  /                   \
                         n+af       \     |  n+af      ~n+af    |   dN
         conv_c_plus_svel_   (x) =   +    | c    (x) + u    (x) | * --- (x)
                                    /     |  j          j       |   dx
                                   +-----  \                   /      j
                                   dim j    +------+   +------+
                                               if         if
                                              supg     reynolds

        */
        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)+=conv_subaf_(nn);
        }
      }
    }

    /* Most recent value for subgrid velocity convective term

                  /~n+af         \   n+af
                 | u      o nabla | u
                  \   (i)        /   (i)
    */
    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs
       ||
       cross == INPAR::FLUID::cross_stress_stab
      )
    {
      for (int rr=0;rr<2;++rr)
      {
        convsubaf_old_(rr)=vderxyaf_(rr,0)*svelaf_(0)+vderxyaf_(rr,1)*svelaf_(1);
      }
    }


    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //
    //                       SYSTEM MATRIX
    //
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    if(compute_elemat)
    {

      // scaling factors for Galerkin 1 terms
      double fac_inertia   =fac*alphaM;
      double fac_convection=fac*afgdt ;

      // select continuity stabilisation
      const double cstabfac=fac*gamma*dt*tauC;

      const double fac_gamma_dt      = fac*gamma*dt;
      const double fac_afgdt_visceff = fac*visceff*afgdt;

      //---------------------------------------------------------------
      //
      //              SUBSCALE ACCELERATION PART
      //        RESCALING FACTORS FOR GALERKIN 1 TERMS AND
      //              COMPUTATION OF EXTRA TERMS
      //
      //---------------------------------------------------------------

      if(inertia == INPAR::FLUID::inertia_stab_keep
         ||
         inertia == INPAR::FLUID::inertia_stab_keep_complete)
      {
        // rescale time factors terms affected by inertia stabilisation
        fac_inertia   *=afgdt*facMtau;
        fac_convection*=afgdt*facMtau;

        // do inertia stabilisation terms which are not scaled
        // Galerkin terms since they are not partially integrated

        const double fac_alphaM_tauM_facMtau = fac*alphaM*tauM*facMtau;

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          const double fac_alphaM_gamma_dt_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*gamma*dt*funct_(vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tuipp  =3*ui+2;
            /* pressure (implicit) */

            /*  factor:
                             alphaM*tauM
                    ---------------------------, rescaled by gamma*dt
                    alphaM*tauM+alphaF*gamma*dt

                 /               \
                |                 |
                |  nabla Dp ,  v  |
                |                 |
                 \               /
            */
            /* pressure (implicit) */

            elemat(tvi  ,tuipp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(0,ui);
            elemat(tvip ,tuipp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(1,ui);
          } // ui
        } // vi

        if(higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_alphaM_tauM_facMtau
            =
            fac*visceff*afgdt*alphaM*tauM*facMtau;

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi    =3*vi;
            const int tvip   =tvi+1;

            const double fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
              =
              fac_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui    =3*ui;
              const int tuip   =tui+1;

              /* viscous term (intermediate) */
              /*  factor:
                                                 alphaM*tauM
                        nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                  /                           \
                 |                 /    \      |
                 |  2*nabla o eps | Dacc | , v |
                 |                 \    /      |
                  \                           /

              */
              elemat(tvi  ,tui  ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,ui);
              elemat(tvi  ,tuip ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(2,ui);
              elemat(tvip ,tui  ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(2,ui);
              elemat(tvip ,tuip ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,ui);
            } // ui
          } // vi
        } // end higher order element and  linearisation of linear terms not supressed

        if(inertia == INPAR::FLUID::inertia_stab_keep_complete)
        {

          /*
                                  immediately enters the matrix
                                  |
                                  v
                               +--------------+
                               |              |
                                /            \
                      1.0      |  ~n+af       |
                 - --------- * |  u     ,  v  |
                        n+af   |   (i)        |
                   tau_M        \            /

                   |       |
                   +-------+
                       ^
                       |
                       consider linearisation of this expression

          */
          const double norm = velintaf_.Norm2();

          // normed velocity at element center (we use the copy for safety reasons!)
          if (norm>=1e-6)
          {
            for (int rr=0;rr<2;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=velintaf_(rr)/norm;
            }
          }
          else
          {
            normed_velintaf_(0) = 0.0;
            for (int rr=1;rr<2;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=0.0;
            }
          }

          double temp=0.0;
          if(whichtau==INPAR::FLUID::codina)
          {
            /*
                                                  || n+af||
                       1.0           visc         ||u    ||
                    --------- = CI * ---- + CII * ---------
                         n+af           2
                    tau_M             hk             hk


                    where CII=2.0/mk
            */

            temp=fac*afgdt/hk*2.0/mk;
          }
          else if(whichtau==INPAR::FLUID::smoothed_franca_barrenechea_valentin_wall)
          {
            /*
                                  -x   '       -x
                    using f(x)=x+e  , f (x)=1-e


                                                +-                                -+
                                                |          / || n+af||          \  |
                       1.0      4.0 * visceff   |         |  ||u    || * hk * mk | |
                    --------- = ------------- * | 1.0 + f |  ------------------- | |
                         n+af           2       |         |                      | |
                    tau_M         mk* hk        |          \    2.0 * visceff   /  |
                                                +-                                -+

            */

            temp=fac*afgdt/hk*2.0*(1-exp(-1.0*(norm*hk/visceff)*(mk/2.0)));


          }
          else if(whichtau==INPAR::FLUID::franca_barrenechea_valentin_wall)
          {

            /*
                                             +-                                  -+
                                             |            / || n+af||          \  |
                       1.0      4.0 * visc   |           |  ||u    || * hk * mk | |
                    --------- = ---------- * | 1.0 + max |  ------------------- | |
                         n+af           2    |           |                      | |
                    tau_M         mk* hk     |            \    2.0 * visceff   /  |
                                             +-                                  -+

            */

            if((norm*hk/visceff)*(mk/2.0)>1)
            {
              temp=fac*afgdt/hk*2.0;
            }
          }
          else
          {
            dserror("There's no linearisation of 1/tau available for this tau definition\n");
          }

          /*
                        || n+af||             n+af
                      d ||u    ||            u    * Dacc
                      ----------- = afgdt *  -----------
                                              || n+af||
                        d Dacc                ||u    ||

          */

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi    =3*vi;
            const int tvip   =tvi+1;


            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui  =3*ui;
              const int tuip =tui+1;

              elemat(tvi  ,tui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(tvi  ,tuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(0);

              elemat(tvip ,tui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(tvip ,tuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(1);
            } // ui
          } // vi
        } // end linearisation of 1/tauM
      } // extra terms for inertia stab

      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //
      //      GALERKIN PART 1 (INERTIA, CONVECTION, VISCOUS)
      // GALERKIN PART 2 (REMAINING PRESSURE AND CONTINUITY EXPRESSIONS)
      //
      //               CONTINUITY STABILISATION
      //
      //---------------------------------------------------------------

      /*
        inertia term (intermediate)

                                                 /          \
                         alphaF*gamma*dt        |            |
             alphaM*---------------------------*|  Dacc , v  |
                    alphaM*tauM+alphaF*gamma*dt |            |
                                                 \          /
             |                                 |
             +---------------------------------+
               	            alphaM
           	without inertia stabilisation

       convection (intermediate)

                                                 /                          \
                         alphaF*gamma*dt        |  / n+af       \            |
    alphaF*gamma*dt*---------------------------*| | c    o nabla | Dacc , v  |
                    alphaM*tauM+alphaF*gamma*dt |  \            /            |
                                                 \                          /
    |                                          |
    +------------------------------------------+
		  +alphaF*gamma*dt
          without inertia stabilisation


      convection (intermediate)
|
|                                                /                            \
N                         alphaF*gamma*dt       |  /            \   n+af       |
E  +alphaF*gamma*dt*---------------------------*| | Dacc o nabla | u      , v  |
W                   alphaM*tauM+alphaF*gamma*dt |  \            /              |
T                                                \                            /
O  |                                          |
N  +------------------------------------------+
              +alphaF*gamma*dt
        without inertia stabilisation


      pressure (implicit)

                                                 /                \
                                                |                  |
                                      -gamma*dt |  Dp , nabla o v  |
                                                |                  |
                                                 \                /

     viscous term (intermediate)


                                                 /                          \
		                                |       /    \         / \   |
                          +2*nu*alphaF*gamma*dt*|  eps | Dacc | , eps | v |  |
                                                |       \    /         \ /   |
                                                 \                          /


     continuity equation (implicit)



                                                 /                  \
                                                |                    |
                                     +gamma*dt* | nabla o Dacc  , q  |
                                                |                    |
                                                 \                  /


      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //               CONTINUITY STABILISATION
      //
      //---------------------------------------------------------------

                                                 /                          \
                                                |                            |
                                +gamma*dt*tauC* | nabla o Dacc  , nabla o v  |
                                                |                            |
                                                 \                          /
                                +-------------+
                               zero for no cstab


      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //
      //                   SUPG STABILISATION
      //            SUPG TYPE REYNOLDS LINEARISATIONS
      //
      //---------------------------------------------------------------
         SUPG stabilisation --- subscale velocity, nonlinear part from testfunction
|
|
N                                       /                            \
E                                      |  ~n+af    /            \     |
W                 alphaF * gamma * dt* |  u     , | Dacc o nabla | v  |
T                                      |   (i)     \            /     |
O                                       \                            /
N

         SUPG stabilisation --- inertia

                              alphaF*gamma*dt
                         --------------------------- * alphaM * tauM *
                         alphaM*tauM+alphaF*gamma*dt


                     /                                        \
                    |          / / n+af  ~n+af \         \     |
                    |  Dacc , | | c    + u      | o nabla | v  |
                    |          \ \             /         /     |
                     \                                        /

        SUPG stabilisation --- convection

                               alphaF*gamma*dt
                         --------------------------- * alphaF * gamma * dt * tauM
                         alphaM*tauM+alphaF*gamma*dt

                     /                                                           \
                    |    / n+af        \          / / n+af  ~n+af \         \     |
                    |   | c     o nabla | Dacc , | | c    + u      | o nabla | v  |
                    |    \             /          \ \             /         /     |
                     \                                                           /

        SUPG stabilisation --- convection

                              alphaF*gamma*dt
|                       --------------------------- * alphaF * gamma * dt * tauM
|                       alphaM*tauM+alphaF*gamma*dt
N
E                   /                                                           \
W                  |    /            \   n+af    / / n+af  ~n+af \         \     |
T                  |   | Dacc o nabla | u     , | | c    + u      | o nabla | v  |
O                  |    \            /           \ \             /         /     |
N                   \                                                           /

        SUPG stabilisation --- pressure

                               alphaF*gamma*dt*tauM
                            ---------------------------, rescaled by gamma*dt
                            alphaM*tauM+alphaF*gamma*dt


                    /                                            \
                   |              / / n+af  ~n+af \         \     |
                   |  nabla Dp , | | c    + u      | o nabla | v  |
                   |              \ \             /         /     |
                    \                                            /

        SUPG stabilisation --- diffusion

                                              alphaF*gamma*dt*tauM
                        nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt

                    /                                                          \
                   |  /             /      \     / / n+af  ~n+af \         \    |
                   | | nabla o eps |  Dacc  | , | | c    + u      | o nabla | v |
                   |  \             \      /     \ \             /         /    |
                    \                                                          /
      */

      const double fac_afgdt_afgdt_tauM_facMtau  = fac*afgdt   *afgdt*tauM*facMtau;
      const double fac_gdt_afgdt_tauM_facMtau    = fac*gamma*dt*afgdt*tauM*facMtau;
      const double fac_alphaM_afgdt_tauM_facMtau = fac*alphaM  *afgdt*tauM*facMtau;

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tui  =3*ui;
        const int tuip =tui+1;

        /* GALERKIN inertia term (intermediate) + convection (intermediate) */
        const double inertia_and_conv_ui
          =
          fac_inertia*funct_(ui)
          +
          fac_convection*conv_c_af_(ui);

        /* viscous term (intermediate), 'diagonal' parts */
        const double visc_0=fac_afgdt_visceff*derxy_(0,ui);
        const double visc_1=fac_afgdt_visceff*derxy_(1,ui);

        /* SUPG stabilisation --- inertia and convection */
        const double supg_inertia_and_conv_ui
          =
          fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_afgdt_tauM_facMtau*conv_c_af_(ui);

        /* CSTAB entries */
        const double cstab_0 = cstabfac*derxy_(0,ui);
        const double cstab_1 = cstabfac*derxy_(1,ui);

        /* combined CSTAB/viscous entires */
        const double visc_and_cstab_0 = visc_0+cstab_0;
        const double visc_and_cstab_1 = visc_1+cstab_1;

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          /* inertia term (intermediate)                */
          /* convection   (intermediate)                */
          /* supg inertia and convection                */
          /* viscous term (intermediate, diagonal part) */
          const double sum =
            inertia_and_conv_ui*funct_(vi)
            +
            supg_inertia_and_conv_ui*conv_c_plus_svel_af_(vi)
            +
            visc_0*derxy_(0,vi)
            +
            visc_1*derxy_(1,vi);

          /* CONTINUITY stabilisation                     */
          /* viscous term (intermediate, remaining parts) */

          const double a=visc_0*derxy_(1,vi)+cstab_1*derxy_(0,vi);

          elemat(tvi  ,tui  ) += sum+visc_and_cstab_0*derxy_(0,vi);
          elemat(tvi  ,tuip ) += a;
          elemat(tuip ,tvi  ) += a;
          elemat(tvip ,tuip ) += sum+visc_and_cstab_1*derxy_(1,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp  =3*ui+2;

        const double fac_gamma_dt_funct_ui=fac_gamma_dt*funct_(ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi  =3*vi;
          const int tvip =tvi+1;

          /* GALERKIN pressure   (implicit), rescaled by gamma*dt */
          /* continuity equation (implicit)                       */

          elemat(tvi  ,tuipp) -= fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(tvip ,tuipp) -= fac_gamma_dt_funct_ui*derxy_(1,vi);

          elemat(tuipp,tvi  ) += fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(tuipp,tvip ) += fac_gamma_dt_funct_ui*derxy_(1,vi);
        } // vi
      } // ui

      if (newton==INPAR::FLUID::Newton) // if newton and supg
      {
        const double fac_afgdt_afgdt_tauM_facMtau = fac*afgdt*afgdt*facMtau*tauM;

        // linearisation of SUPG testfunction
        double temp[2][2];

        const double fac_afgdt_svelaf_0 = fac*afgdt*supg_active*svelaf_(0);
        const double fac_afgdt_svelaf_1 = fac*afgdt*supg_active*svelaf_(1);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi  =3*vi;
          const int tvip =tvi+1;

          // linearisations of reactive Galerkin part (remaining after inertia_stab)
          // and SUPG part (reactive part from residual)
          const double scaled_inertia_and_conv_vi
            =
            fac_convection*funct_(vi)
            +
            fac_afgdt_afgdt_tauM_facMtau*conv_c_plus_svel_af_(vi);

          temp[0][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,0)-fac_afgdt_svelaf_0*derxy_(0,vi);
          temp[1][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,1)-fac_afgdt_svelaf_0*derxy_(1,vi);
          temp[0][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,0)-fac_afgdt_svelaf_1*derxy_(0,vi);
          temp[1][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,1)-fac_afgdt_svelaf_1*derxy_(1,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
            elemat(tvi  ,tuip ) += temp[1][0]*funct_(ui);
            elemat(tvip ,tui  ) += temp[0][1]*funct_(ui);
            elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);
          } // ui
        } // vi
      } // end if newton

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp=3*ui+2;

        const double scaled_gradp_0 = fac_gdt_afgdt_tauM_facMtau*derxy_(0,ui);
        const double scaled_gradp_1 = fac_gdt_afgdt_tauM_facMtau*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi=3*vi;

          /* SUPG stabilisation --- pressure, rescaled by gamma*dt */
          elemat(tvi  ,tuipp) += scaled_gradp_0*conv_c_plus_svel_af_(vi);
          elemat(tvi+1,tuipp) += scaled_gradp_1*conv_c_plus_svel_af_(vi);
        } // vi
      } // ui

      if(higher_order_ele && newton!=INPAR::FLUID::minimal)
      {
        const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          double coltemp[2][2];

          coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(0,ui);
          coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(2,ui);

          coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(2,ui);
          coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(1,ui);


          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*  SUPG stabilisation, diffusion */
            elemat(tvi  ,tui  ) -= coltemp[0][0]*conv_c_plus_svel_af_(vi);
            elemat(tvi  ,tuip ) -= coltemp[0][1]*conv_c_plus_svel_af_(vi);

            elemat(tvip ,tui  ) -= coltemp[1][0]*conv_c_plus_svel_af_(vi);
            elemat(tvip ,tuip ) -= coltemp[1][1]*conv_c_plus_svel_af_(vi);
          } // vi
        } // ui
      } // hoel

      //---------------------------------------------------------------
      //
      //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //
      //                    PRESSURE STABILISATION
      //
      //---------------------------------------------------------------
      if(pspg == INPAR::FLUID::pstab_use_pspg)
      {
        const double fac_afgdt_gamma_dt_tauM_facMtau  = fac*afgdt*gamma*dt*tauM*facMtau;
        const double fac_gdt_gdt_tauM_facMtau         = fac*gamma*dt*tauM*facMtau*gamma*dt;
        const double fac_alphaM_gamma_dt_tauM_facMtau = fac*alphaM*gamma*dt*tauM*facMtau;

        if(higher_order_ele  && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_gamma_dt_tauM_facMtau
            =
            fac*visceff*afgdt*gamma*dt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);


            const double pspg_diffusion_inertia_convect_0_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(0,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_1_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(1,ui)-inertia_and_conv_ui;

            const double scaled_derxy2_2_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(2,ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvipp =3*vi+2;

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

              /* pressure stabilisation --- diffusion  */


              /*
                           gamma*dt*tau_M
            factor:  ------------------------------ * alpha_F*gamma*dt * nu
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                                  \
                   |                 /    \             |
                   |  2*nabla o eps | Dacc | , nabla q  |
                   |                 \    /             |
                    \                                  /
              */

              elemat(tvipp,tui  ) -= derxy_(0,vi)*pspg_diffusion_inertia_convect_0_ui+derxy_(1,vi)*scaled_derxy2_2_ui;
              elemat(tvipp,tuip ) -= derxy_(1,vi)*pspg_diffusion_inertia_convect_1_ui+derxy_(0,vi)*scaled_derxy2_2_ui ;
            }
          }
        }
        else
        {
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi) // loop rows
            {
              const int tvipp =3*vi+2;

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

              elemat(tvipp,tui  ) +=derxy_(0,vi)*inertia_and_conv_ui;
              elemat(tvipp,tuip ) +=derxy_(1,vi)*inertia_and_conv_ui;
            }
          }
        } // neglect viscous linearisations, do just inertia and convective

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp=3*ui+2;
          const double scaled_derxy_0=fac_gdt_gdt_tauM_facMtau*derxy_(0,ui);
          const double scaled_derxy_1=fac_gdt_gdt_tauM_facMtau*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            /* pressure stabilisation --- pressure   */

            /*
                          gamma*dt*tau_M
            factor:  ------------------------------, rescaled by gamma*dt
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
            */

            elemat(vi*3+2,tuipp) +=
              (scaled_derxy_0*derxy_(0,vi)
               +
               scaled_derxy_1*derxy_(1,vi)) ;

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (newton==INPAR::FLUID::Newton) // if pspg and newton
        {

          for (int vi=0; vi<iel; ++vi) // loop columns
          {
            const int tvipp=3*vi+2;

            const double a=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0));
            const double b=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1));

            for (int ui=0; ui<iel; ++ui)  // loop rows
            {
              const int tui=3*ui;
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

              elemat(tvipp,tui  ) += a*funct_(ui);
              elemat(tvipp,tui+1) += b*funct_(ui);
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
        }// end if pspg and newton
      } // end pressure stabilisation

      //---------------------------------------------------------------
      //
      //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //            VISCOUS STABILISATION TERMS FOR (A)GLS
      //
      //---------------------------------------------------------------
      if (higher_order_ele)
      {
        if(vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_gls)
        {
          const double tauMqs = afgdt*tauM*facMtau;

          const double fac_visc_tauMqs_alphaM        = vstabfac*fac*visc*tauMqs*alphaM;
          const double fac_visc_tauMqs_afgdt         = vstabfac*fac*visc*tauMqs*afgdt;
          const double fac_visc_tauMqs_afgdt_visceff = vstabfac*fac*visc*tauMqs*afgdt*visceff;
          const double fac_visc_tauMqs_gamma_dt      = vstabfac*fac*visc*tauMqs*gamma*dt;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;

            const double inertia_and_conv
              =
              fac_visc_tauMqs_alphaM*funct_(ui)+fac_visc_tauMqs_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;
              /* viscous stabilisation --- inertia     */

              /* factor:

                                        alphaF*gamma*tauM*dt
                     +(-)alphaM*nu* ---------------------------
                                    alphaM*tauM+alphaF*gamma*dt

                     /                      \
                    |                        |
                    |  Dacc , 2*div eps (v)  |
                    |                        |
                     \                      /
              */

              /* viscous stabilisation --- convection */
              /*  factor:
                                         alphaF*gamma*dt*tauM
              +(-)alphaF*gamma*dt*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                       /                                    \
                      |  / n+af       \                      |
                      | | c    o nabla | Dacc, 2*div eps (v) |
                      |  \            /                      |
                       \                                    /

              */


              const double a = inertia_and_conv*derxy2_(2,vi);

              elemat(tvi  ,tui  ) += inertia_and_conv*viscs2_(0,vi);
              elemat(tvi  ,tuip ) += a;
              elemat(tvip ,tui  ) += a;
              elemat(tvip ,tuip ) += inertia_and_conv*viscs2_(1,vi);
            }
          }

          for (int ui=0;ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;

            for (int vi=0;vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /* viscous stabilisation --- diffusion  */

              /* factor:

                                             alphaF*gamma*tauM*dt
                -(+)alphaF*gamma*dt*nu*nu ---------------------------
                                          alphaM*tauM+alphaF*gamma*dt

                    /                                        \
                   |                  /    \                  |
                   |  2* nabla o eps | Dacc | , 2*div eps (v) |
                   |                  \    /                  |
                    \                                        /
              */

              const double a = fac_visc_tauMqs_afgdt_visceff*
                               (viscs2_(0,vi)*derxy2_(2,ui)
                                +
                                derxy2_(2,vi)*viscs2_(1,ui));

              elemat(tvi  ,tuip ) -= a;
              elemat(tuip ,tvi  ) -= a;


              elemat(tvi   ,tui ) -= fac_visc_tauMqs_afgdt_visceff*
		                     (viscs2_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy2_(2,ui)*derxy2_(2,vi));

              elemat(tvip ,tuip ) -= fac_visc_tauMqs_afgdt_visceff*
                                     (derxy2_(2,ui)*derxy2_(2,vi)
                                      +
                                      viscs2_(1,ui)*viscs2_(1,vi));

            }
          }

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;
            const int tuipp  =tuip+1;

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /* viscous stabilisation --- pressure   */

              /* factor:

                                    alphaF*gamma*tauM*dt
                       +(-)nu * ---------------------------, rescaled by gamma*dt
                                alphaM*tauM+alphaF*gamma*dt


                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
              */
              elemat(tvi  ,tuipp) += fac_visc_tauMqs_gamma_dt*
                                     (derxy_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy_(1,ui)*derxy2_(2,vi)) ;
              elemat(tvip ,tuipp) += fac_visc_tauMqs_gamma_dt*
                                     (derxy_(0,ui)*derxy2_(2,vi)
                                      +
                                      derxy_(1,ui)*viscs2_(1,vi)) ;
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton==INPAR::FLUID::Newton)
          {

            double temp[2][2];
            for (int vi=0; vi<iel; ++vi) // loop columns (solution for matrix, test function for vector)
            {
              const int tvi    =3*vi;
              const int tvip   =tvi+1;

              temp[0][0]=(viscs2_(0,vi)*vderxyaf_(0,0)
			  +
                          derxy2_(2,vi)*vderxyaf_(1,0))*fac_visc_tauMqs_afgdt;
              temp[1][0]=(viscs2_(0,vi)*vderxyaf_(0,1)
			  +
                          derxy2_(2,vi)*vderxyaf_(1,1))*fac_visc_tauMqs_afgdt;
              temp[0][1]=(derxy2_(2,vi)*vderxyaf_(0,0)
			  +
                          viscs2_(1,vi)*vderxyaf_(1,0))*fac_visc_tauMqs_afgdt;
              temp[1][1]=(derxy2_(2,vi)*vderxyaf_(0,1)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,1))*fac_visc_tauMqs_afgdt;

              for (int ui=0; ui<iel; ++ui)  // loop rows (test functions for matrix)
              {
                const int tui    =3*ui;
                const int tuip   =tui+1;

                /* viscous stabilisation --- convection
                     factor:
                                         alphaF*gamma*dt*tauM
              +(-)alphaF*gamma*dt*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /


                */
                elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
                elemat(tvi  ,tuip ) += temp[1][0]*funct_(ui);
                elemat(tvip ,tui  ) += temp[0][1]*funct_(ui);
                elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);

              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)

          } // end if (a)gls and newton
        } // end (a)gls stabilisation
      } // end higher_order_element

      //---------------------------------------------------------------
      //
      //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == INPAR::FLUID::cross_stress_stab)
      {
        const double fac_afgdt=fac*afgdt;

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const double fac_afgdt_conv_subaf_ui=fac_afgdt*conv_subaf_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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

            elemat(vi*3    , ui*3    ) += fac_afgdt_conv_subaf_ui_funct_vi;
            elemat(vi*3 + 1, ui*3 + 1) += fac_afgdt_conv_subaf_ui_funct_vi;
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

          /*
                                                  alphaM*tauM
                          -alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt

          */
        const double fac_afgdt_alphaM_tauM_facMtau = fac*afgdt*alphaM*tauM*facMtau;
        /*

                                              alphaF*gamma*dt*tauM
                          -alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt
        */
        const double fac_afgdt_afgdt_tauM_facMtau  = fac*afgdt*afgdt*tauM*facMtau;

        double  am_nabla_u_afgdt_nabla_u_nabla_u[2][2];

        am_nabla_u_afgdt_nabla_u_nabla_u[0][0]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(0,0)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*vderxyaf_(0,0)
                                        +
                                        vderxyaf_(0,1)*vderxyaf_(1,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][1]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(0,1)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*vderxyaf_(0,1)
                                        +
                                        vderxyaf_(0,1)*vderxyaf_(1,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][0]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(1,0)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(1,0)*vderxyaf_(0,0)
                                        +
                                        vderxyaf_(1,1)*vderxyaf_(1,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][1]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(1,1)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(1,0)*vderxyaf_(0,1)
                                        +
                                        vderxyaf_(1,1)*vderxyaf_(1,1));

        double nabla_u[2][2];
        nabla_u[0][0]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(0,0);
        nabla_u[0][1]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(0,1);

        nabla_u[1][0]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(1,0);
        nabla_u[1][1]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(1,1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_nabla_ui = velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double coltemp[2][2];

          coltemp[0][0]=am_nabla_u_afgdt_nabla_u_nabla_u[0][0]*funct_(ui)+nabla_u[0][0]*u_nabla_ui;
          coltemp[0][1]=am_nabla_u_afgdt_nabla_u_nabla_u[0][1]*funct_(ui)+nabla_u[0][1]*u_nabla_ui;

          coltemp[1][0]=am_nabla_u_afgdt_nabla_u_nabla_u[1][0]*funct_(ui)+nabla_u[1][0]*u_nabla_ui;
          coltemp[1][1]=am_nabla_u_afgdt_nabla_u_nabla_u[1][1]*funct_(ui)+nabla_u[1][1]*u_nabla_ui;

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*  factor:

                                                  alphaM*tauM
                          -alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                -alphaM*tauM

                          /                           \
                         |  /            \   n+af      |
                         | | Dacc o nabla | u     , v  |
                         |  \            /             |
                          \                           /
            */

            /*  factor:

                                              alphaF*gamma*dt*tauM
                          -alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt



                          /                                                \
                         |  / / /            \   n+af \         \   n+af    |
                         | | | | Dacc o nabla | u      | o nabla | u   , v  |
                         |  \ \ \            /        /         /           |
                          \                                                /
            */

            /*  factor:

                                              alphaF*gamma*dt*tauM
                          -alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt

                          /                                                 \
                         |  / / / n+af        \       \         \   n+af     |
                         | | | | u     o nabla | Dacc  | o nabla | u    , v  |
                         |  \ \ \             /       /         /            |
                          \                                                 /
            */

            elemat(tvi  ,tui  ) -= funct_(vi)*coltemp[0][0];
            elemat(tvi  ,tuip ) -= funct_(vi)*coltemp[0][1];

            elemat(tvip ,tui  ) -= funct_(vi)*coltemp[1][0];
            elemat(tvip ,tuip ) -= funct_(vi)*coltemp[1][1];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        const double fac_afgdt_tauM_facMtau_gdt = fac*alphaF*gamma*dt*tauM*facMtau*gamma*dt;

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tuipp =3*ui+2;

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*

               factor:
                                alpha_F*gamma*dt*tau_M
                            ------------------------------, rescaled by gamma*dt
                            alpha_M*tau_M+alpha_F*gamma*dt


                          /                               \
                         |  /                \   n+af      |
                         | | nabla Dp o nabla | u     , v  |
                         |  \                /             |
                          \                               /
            */
            elemat(tvi  ,tuipp) -= fac_afgdt_tauM_facMtau_gdt*funct_(vi)*(vderxyaf_(0,0)*derxy_(0,ui)+vderxyaf_(0,1)*derxy_(1,ui));
            elemat(tvip ,tuipp) -= fac_afgdt_tauM_facMtau_gdt*funct_(vi)*(vderxyaf_(1,0)*derxy_(0,ui)+vderxyaf_(1,1)*derxy_(1,ui));
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            double coltemp[2][2];

            coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(0,ui)+derxy2_(2,ui)*vderxyaf_(0,1));
            coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(1,ui)+derxy2_(2,ui)*vderxyaf_(0,0));

            coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(0,ui)+derxy2_(2,ui)*vderxyaf_(1,1));
            coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(1,ui)+derxy2_(2,ui)*vderxyaf_(1,0));

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /*  factor:

                                              alphaF*gamma*dt*tauM
                        nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt

                          /                                              \
                         |  / /             /      \         \   n+af     |
                         | | | nabla o eps |  Dacc  | o nabla | u    , v  |
                         |  \ \             \      /         /            |
                          \                                              /
              */
              elemat(tvi  ,tui  ) += coltemp[0][0]*funct_(vi);
              elemat(tvi  ,tuip ) += coltemp[0][1]*funct_(vi);

              elemat(tvip ,tui  ) += coltemp[1][0]*funct_(vi);
              elemat(tvip ,tuip ) += coltemp[1][1]*funct_(vi);

            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
        } // end if higher_order_element
      } // end cross

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {
        /*
                  /                            \
                 |  ~n+af    ~n+af              |
               - |  u    , ( u     o nabla ) v  |
                 |                              |
                  \                            /
                             +----+
                               ^
                               |
                               linearisation of this expression
        */
        const double fac_alphaM_afgdt_tauM_facMtau=fac*alphaM*afgdt*tauM*facMtau;

        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_x = fac_alphaM_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_y = fac_alphaM_afgdt_tauM_facMtau*svelaf_(1);

        const double fac_afgdt_afgdt_tauM_facMtau =fac*afgdt*afgdt*tauM*facMtau;

        double fac_afgdt_afgdt_tauM_facMtau_svelaf[2];
        fac_afgdt_afgdt_tauM_facMtau_svelaf[0]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[1]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double inertia_and_conv[2];

          inertia_and_conv[0]=fac_afgdt_afgdt_tauM_facMtau_svelaf[0]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_afgdt_tauM_facMtau_svelaf[1]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_y*funct_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
               factor: +alphaM * alphaF * gamma * dt * tauM * facMtau

                  /                            \
                 |  ~n+af                       |
                 |  u     , ( Dacc o nabla ) v  |
                 |                              |
                  \                            /

            */

            /*
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / / n+af        \       \         \     |
             |  u     , | | | u     o nabla | Dacc  | o nabla | v  |
             |           \ \ \             /       /         /     |
              \                                                   /

            */

            elemat(tvi  ,tui  ) += inertia_and_conv[0]*derxy_(0,vi);
            elemat(tvi  ,tuip ) += inertia_and_conv[0]*derxy_(1,vi);

            elemat(tvip ,tui  ) += inertia_and_conv[1]*derxy_(0,vi);
            elemat(tvip ,tuip ) += inertia_and_conv[1]*derxy_(1,vi);
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          double temp[2];
          temp[0]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi));
          temp[1]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi));

          double rowtemp[2][2];

          rowtemp[0][0]=svelaf_(0)*temp[0];
          rowtemp[0][1]=svelaf_(0)*temp[1];

          rowtemp[1][0]=svelaf_(1)*temp[0];
          rowtemp[1][1]=svelaf_(1)*temp[1];

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            /*
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / /            \   n+af \         \     |
             |  u     , | | | Dacc o nabla | u      | o nabla | v  |
             |           \ \ \            /        /         /     |
              \                                                   /

            */

            elemat(tvi  ,tui  ) += funct_(ui)*rowtemp[0][0];
            elemat(tvi  ,tuip ) += funct_(ui)*rowtemp[0][1];

            elemat(tvip ,tui  ) += funct_(ui)*rowtemp[1][0];
            elemat(tvip ,tuip ) += funct_(ui)*rowtemp[1][1];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        const double fac_gdt_afgdt_tauM_facMtau         =fac*gamma*dt*afgdt*tauM*facMtau;
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_x=fac_gdt_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_y=fac_gdt_afgdt_tauM_facMtau*svelaf_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tuipp =3*ui+2;

          double coltemp[2][2];

          coltemp[0][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(1,ui);

          coltemp[1][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: + gamma * dt * alphaF * gamma * dt * tauM *facMtau (rescaled)

              /                                \
             |  ~n+af    /                \     |
             |  u     , | nabla Dp o nabla | v  |
             |           \                /     |
              \                                /

            */

            elemat(tvi  ,tuipp) += coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi);
            elemat(tvip ,tuipp) += coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi);

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_nu_afgdt_afgdt_tauM_facMtau =fac*visceff*afgdt*afgdt*tauM*facMtau;

          double temp[2];

          temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(0);
          temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(1);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            double rowtemp[2][2];

            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);

            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const int tui   =3*ui;
              const int tuip  =tui+1;

              /*
                   factor: - 2.0 * visc * alphaF * gamma * dt * alphaF * gamma * dt * tauM * facMtauM

                    /                                                 \
                   |  ~n+af    / /             /    \  \         \     |
                   |  u     , | | nabla o eps | Dacc |  | o nabla | v  |
                   |           \ \             \    /  /         /     |
                    \                                                 /
              */

              elemat(tvi  ,tui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
              elemat(tvi  ,tuip ) -= derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

              elemat(tvip ,tui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
              elemat(tvip ,tuip ) -= derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
            }
          }
        }// end higher order ele
      } // end if reynolds stab

    } // end if compute_elemat

    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //
    //                       RIGHT HAND SIDE
    //
    //---------------------------------------------------------------
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    //
    // (MODIFIED) GALERKIN PART, SUBSCALE ACCELERATION STABILISATION
    //
    //---------------------------------------------------------------
    if(inertia == INPAR::FLUID::inertia_stab_keep
       ||
       inertia == INPAR::FLUID::inertia_stab_keep_complete)
    {

      double aux_x =(-svelaf_(0)/tauM-pderxynp_(0));
      double aux_y =(-svelaf_(1)/tauM-pderxynp_(1));

      if(higher_order_ele)
      {
        const double fact =visceff;

        aux_x += fact*viscaf_old_(0);
        aux_y += fact*viscaf_old_(1);
      }

      const double fac_sacc_plus_resM_not_partially_integrated_x =fac*aux_x ;
      const double fac_sacc_plus_resM_not_partially_integrated_y =fac*aux_y ;

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        const int tvi=3*vi;
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

        elevec(tvi  ) -= fac_sacc_plus_resM_not_partially_integrated_x*funct_(vi);
        elevec(tvi+1) -= fac_sacc_plus_resM_not_partially_integrated_y*funct_(vi);
      } // vi
    }
    else
    {
      //---------------------------------------------------------------
      //
      //        GALERKIN PART, NEGLECTING SUBSCALE ACCLERATIONS
      //
      //---------------------------------------------------------------
      const double fac_inertia_convection_dead_load_x
        =
        fac*(accintam_(0)+convaf_old_(0)-bodyforceaf_(0));

      const double fac_inertia_convection_dead_load_y
        =
        fac*(accintam_(1)+convaf_old_(1)-bodyforceaf_(1));

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        const int tvi=3*vi;
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

        elevec(tvi  ) -= funct_(vi)*fac_inertia_convection_dead_load_x;
        elevec(tvi+1) -= funct_(vi)*fac_inertia_convection_dead_load_y;
      } // vi
    }
    //---------------------------------------------------------------
    //
    //            GALERKIN PART 2, REMAINING EXPRESSIONS
    //
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    //
    //         RESIDUAL BASED CONTINUITY STABILISATION
    //          (the original version proposed by Codina)
    //
    //---------------------------------------------------------------

    const double fac_prenp_      = fac*prenp_-fac*tauC*divunp_;

    for (int vi=0; vi<iel; ++vi) // loop rows
    {
      const int tvi =3*vi;
      /* pressure */

      /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
      */

      /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
      */

      elevec(tvi  ) += fac_prenp_*derxy_(0,vi) ;
      elevec(tvi+1) += fac_prenp_*derxy_(1,vi) ;
    } // vi

    const double visceff_fac=visceff*fac;

    for (int vi=0; vi<iel; ++vi) // loop rows
    {
      const int tvi=3*vi;

      /* viscous term */

      /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
      */

      elevec(tvi   ) -= visceff_fac*
                        (derxy_(0,vi)*vderxyaf_(0,0)*2.0
                         +
                         derxy_(1,vi)*(vderxyaf_(0,1)+vderxyaf_(1,0)));
      elevec(tvi+1) -= visceff_fac*
	               (derxy_(0,vi)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                        +
                        derxy_(1,vi)*vderxyaf_(1,1)*2.0);
    } // vi

    const double fac_divunp  = fac*divunp_;

    for (int vi=0; vi<iel; ++vi) // loop rows
    {
      /* continuity equation */

      /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
      */

      elevec(vi*3 + 2) -= fac_divunp*funct_(vi);
    } // vi

    //---------------------------------------------------------------
    //
    //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //                    PRESSURE STABILISATION
    //
    //---------------------------------------------------------------
//    if(pspg == INPAR::FLUID::pstab_use_pspg)
//    {

//      const double fac_svelnpx = fac*(ele->Svelnp())(0,iquad);
//      const double fac_svelnpy = fac*(ele->Svelnp())(1,iquad);
//
//      for (int vi=0; vi<iel; ++vi) // loop rows
//      {
        /* factor: -1

                       /                 \
                      |  ~n+1             |
                      |  u    , nabla  q  |
                      |   (i)             |
                       \                 /
        */
//
//        elevec(vi*3 + 2) += fac_svelnpx*derxy_(0,vi)+fac_svelnpy*derxy_(1,vi);
//
//      } // vi
//    }

    //---------------------------------------------------------------
    //
    //         STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
    //
    //---------------------------------------------------------------
    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      const double fac_svelaf_x=fac*svelaf_(0);
      const double fac_svelaf_y=fac*svelaf_(1);

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        const int tvi=3*vi;
        /*
                  /                             \
                 |  ~n+af    / n+af        \     |
                 |  u     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /

        */

        elevec(tvi  ) += fac_svelaf_x*conv_c_af_(vi);
        elevec(tvi+1) += fac_svelaf_y*conv_c_af_(vi);
      } // vi
    }

    //---------------------------------------------------------------
    //
    //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //             VISCOUS STABILISATION (FOR (A)GLS)
    //
    //---------------------------------------------------------------
    if (higher_order_ele)
    {
      if (vstab != INPAR::FLUID::viscous_stab_none)
      {
        const double fac_visc_svelaf_x = vstabfac*fac*visc*svelaf_(0);
        const double fac_visc_svelaf_y = vstabfac*fac*visc*svelaf_(1);

        for (int vi=0; vi<iel; ++vi) // loop rows
        {
          const int tvi=3*vi;
          /*
                 /                        \
                |  ~n+af                   |
                |  u      , 2*div eps (v)  |
                |                          |
                 \                        /

          */
          elevec(tvi  ) += fac_visc_svelaf_x*viscs2_(0,vi)
	                   +
                           fac_visc_svelaf_y*derxy2_(2,vi);

          elevec(tvi+1) += fac_visc_svelaf_x*derxy2_(2,vi)
                           +
                           fac_visc_svelaf_y*viscs2_(1,vi);

        } // end vi
      } // endif (a)gls
    }// end if higher order ele

    //---------------------------------------------------------------
    //
    //        TIME-DEPENDENT SUBGRID-SCALE STABILISATION
    //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
    //
    //---------------------------------------------------------------
    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs
       ||
       cross == INPAR::FLUID::cross_stress_stab)
    {
      const double fac_convsubaf_old_x=fac*convsubaf_old_(0);
      const double fac_convsubaf_old_y=fac*convsubaf_old_(1);

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        const int tvi=3*vi;

        /* factor:

                  /                           \
                 |   ~n+af           n+af      |
                 | ( u    o nabla ) u     , v  |
                 |    (i)            (i)       |
                  \                           /
        */
        elevec(tvi  ) -= fac_convsubaf_old_x*funct_(vi);
        elevec(tvi+1) -= fac_convsubaf_old_y*funct_(vi);
      } // vi
    } // end cross

    //---------------------------------------------------------------
    //
    //       TIME DEPENDENT SUBGRID-SCALE STABILISATION PART
    //     RESIDUAL BASED VMM STABILISATION --- REYNOLDS STRESS
    //
    //---------------------------------------------------------------
    if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
    {
      const double fac_svelaf_x=fac*svelaf_(0);
      const double fac_svelaf_y=fac*svelaf_(1);

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
        const int tvi=3*vi;

        /* factor:

                  /                             \
                 |  ~n+af      ~n+af             |
                 |  u      , ( u    o nabla ) v  |
                 |                               |
                  \                             /
        */
        elevec(tvi  ) += fac_svelaf_x*(svelaf_(0)*derxy_(0,vi)
                                       +
                                       svelaf_(1)*derxy_(1,vi));
        elevec(tvi+1) += fac_svelaf_y*(svelaf_(0)*derxy_(0,vi)
                                       +
                                       svelaf_(1)*derxy_(1,vi));
      } // vi
    } // end reynolds

  } // end loop iquad

  return;
} // Sysmat_adv_td

/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  conservative, quasistatic version                                  |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_cons_qs_2D(
  Fluid*                                ele,
  std::vector<Epetra_SerialDenseVector>& myknots,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
  LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
  const LINALG::Matrix<nsd_,iel>&           edispnp,
  const LINALG::Matrix<nsd_,iel>&           egridvaf,
  const LINALG::Matrix<nsd_,iel>&           evelnp,
  const LINALG::Matrix<iel,1>&                eprenp,
  const LINALG::Matrix<nsd_,iel>&             eaccam,
  const LINALG::Matrix<nsd_,iel>&             evelaf,
  Teuchos::RCP<const MAT::Material>           material,
  const double                                alphaM,
  const double                                alphaF,
  const double                                gamma,
  const double                                dt,
  const double                                time,
  const enum INPAR::FLUID::LinearisationAction      newton,
  const bool                                  higher_order_ele,
  const enum INPAR::FLUID::PSPG               pspg            ,
  const enum INPAR::FLUID::SUPG               supg            ,
  const enum INPAR::FLUID::VStab              vstab           ,
  const enum INPAR::FLUID::CStab              cstab           ,
  const enum INPAR::FLUID::CrossStress        cross           ,
  const enum INPAR::FLUID::ReynoldsStress     reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha             whichtau,
  double&                                visceff,
  const bool                             compute_elemat
  )
{
  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt
  //
  const double timealphaF = time-(1-alphaF)*dt;

  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }

  //------------------------------------------------------------------
  //                    SET ALL ELEMENT DATA
  // o including element geometry (node coordinates)
  // o including dead loads in nodes
  // o including hk, mk, element volume
  // o including material viscosity, effective viscosity by
  //   Non-Newtonian fluids or fine/large scale Smagorinsky models
  //------------------------------------------------------------------

  double hk   = 0.0;
  double mk   = 0.0;
  double visc = 0.0;

  SetElementData_2D(ele            ,
                 edispnp        ,
                 evelaf         ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 visceff        );

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  //const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);

    /*
      This is the operator

                  /               \
                 | resM    o nabla |
                  \    (i)        /

                  required for the cross stress linearisation
    */
    //
    //                    +-----
    //          n+af       \         n+af      dN
    // conv_resM    (x) =   +    resM    (x) * --- (x)
    //                     /         j         dx
    //                    +-----                 j
    //                     dim j
    if(cross == INPAR::FLUID::cross_stress_stab)
    {
      for(int nn=0;nn<iel;++nn)
      {
        conv_resM_(nn)=resM_(0)*derxy_(0,nn);

        for(int rr=1;rr<2;++rr)
        {
          conv_resM_(nn)+=resM_(rr)*derxy_(rr,nn);
        }
      }
    }

    // get convective linearisation (n+alpha_F,i) at integration point
    // (convection by grid velocity)
    //
    //                    +-----
    //         n+af        \      n+af      dN
    // conv_u_G_    (x) =   +    u    (x) * --- (x)
    //                     /      G,j       dx
    //                    +-----              j
    //                    dim j
    //
    if(ele->IsAle())
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=u_G_af_(0)*derxy_(0,nn);

	for(int rr=1;rr<2;++rr)
	{
	  conv_u_G_af_(nn)+=u_G_af_(rr)*derxy_(rr,nn);
	}
      }
    }
    else
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=0.0;
      }
    }

    /* Convective term  u_G_old * grad u_old: */
    /*
    //     /    n+af        \   n+af
    //    |  u_G     o nabla | u
    //     \                /
    */
    for(int rr=0;rr<2;++rr)
    {
      convu_G_af_old_(rr)=u_G_af_(0)*vderxyaf_(rr,0);
      for(int mm=1;mm<2;++mm)
      {
        convu_G_af_old_(rr)+=u_G_af_(mm)*vderxyaf_(rr,mm);
      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_quasistatic,gamma,dt,hk,mk,visceff);

    // stabilisation parameters
    const double tauM   = tau_(0);
    const double tauMp  = tau_(1);

    if(cstab == INPAR::FLUID::continuity_stab_none)
    {
      tau_(2)=0.0;
    }
    const double tauC    = tau_(2);

    double supg_active_tauM;
    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      supg_active_tauM=tauM;
    }
    else
    {
      supg_active_tauM=0.0;
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

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //
    //              SYSTEM MATRIX, QUASISTATIC FORMULATION
    //
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    if(compute_elemat)
    {
      /* get combined convective linearisation (n+alpha_F,i) at
         integration point
         takes care of half of the linearisation of reynolds part
         (if necessary)


                         n+af
         conv_c_plus_svel_   (x) =


                   +-----  /                   \
                    \     |  n+af      ~n+af    |   dN
            = tauM * +    | c    (x) + u    (x) | * --- (x)
                    /     |  j          j       |   dx
                   +-----  \                   /      j
                    dim j
                           +-------+  +-------+
                              if         if
                             supg      reynolds

      */
      for(int nn=0;nn<iel;++nn)
      {
        conv_c_plus_svel_af_(nn)=supg_active_tauM*conv_c_af_(nn);
      }

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {

        /* half of the reynolds linearisation is done by modifiing
           the supg testfunction, see above */

        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(0)*derxy_(0,nn);

          for(int rr=1;rr<2;++rr)
          {
            conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(rr)*derxy_(rr,nn);
          }
        }

        /*
                  /                           \
                 |                             |
                 |  resM , ( resM o nabla ) v  |
                 |                             |
                  \                           /
                            +----+
                              ^
                              |
                              linearisation of this expression
        */
        const double fac_alphaM_tauM_tauM=fac*alphaM*tauM*tauM;

        const double fac_alphaM_tauM_tauM_resM_x=fac_alphaM_tauM_tauM*resM_(0);
        const double fac_alphaM_tauM_tauM_resM_y=fac_alphaM_tauM_tauM*resM_(1);

        const double fac_afgdt_tauM_tauM=fac*afgdt*tauM*tauM;

        double fac_afgdt_tauM_tauM_resM[2];
        fac_afgdt_tauM_tauM_resM[0]=fac_afgdt_tauM_tauM*resM_(0);
        fac_afgdt_tauM_tauM_resM[1]=fac_afgdt_tauM_tauM*resM_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double inertia_and_conv[2];

          inertia_and_conv[0]=fac_afgdt_tauM_tauM_resM[0]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_tauM_tauM_resM[1]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_y*funct_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: -alphaM * tauM * tauM

                  /                           \
                 |                             |
                 |  resM , ( Dacc o nabla ) v  |
                 |                             |
                  \                           /

            */

            /*
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / / n+af        \       \         \     |
             |  resM , | | | u     o nabla | Dacc  | o nabla | v  |
             |          \ \ \             /       /         /     |
              \                                                  /

            */

            elemat(tvi  ,tui  ) -= inertia_and_conv[0]*derxy_(0,vi);
            elemat(tvi  ,tuip ) -= inertia_and_conv[0]*derxy_(1,vi);

            elemat(tvip ,tui  ) -= inertia_and_conv[1]*derxy_(0,vi);
            elemat(tvip ,tuip ) -= inertia_and_conv[1]*derxy_(1,vi);
          } // vi
        } // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          double temp[2];
          temp[0]=fac_afgdt_tauM_tauM*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi));
          temp[1]=fac_afgdt_tauM_tauM*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi));

          double rowtemp[2][2];

          rowtemp[0][0]=resM_(0)*temp[0];
          rowtemp[0][1]=resM_(0)*temp[1];

          rowtemp[1][0]=resM_(1)*temp[0];
          rowtemp[1][1]=resM_(1)*temp[1];

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            /*
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / /            \   n+af \         \     |
             |  resM , | | | Dacc o nabla | u      | o nabla | v  |
             |          \ \ \            /        /         /     |
              \                                                  /

            */

            elemat(tvi  ,tui  ) -= funct_(ui)*rowtemp[0][0];
            elemat(tvi  ,tuip ) -= funct_(ui)*rowtemp[0][1];

            elemat(tvip ,tui  ) -= funct_(ui)*rowtemp[1][0];
            elemat(tvip ,tuip ) -= funct_(ui)*rowtemp[1][1];
          } // ui
        } // vi


        const double fac_gdt_tauM_tauM       =fac*gamma*dt*tauM*tauM;
        const double fac_gdt_tauM_tauM_resM_x=fac_gdt_tauM_tauM*resM_(0);
        const double fac_gdt_tauM_tauM_resM_y=fac_gdt_tauM_tauM*resM_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp =3*ui+2;

          double coltemp[2][2];

          coltemp[0][0]=fac_gdt_tauM_tauM_resM_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_tauM_tauM_resM_x*derxy_(1,ui);
          coltemp[1][0]=fac_gdt_tauM_tauM_resM_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_tauM_tauM_resM_y*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: - gamma * dt * tauM * tauM (rescaled)

              /                               \
             |          /                \     |
             |  resM , | nabla Dp o nabla | v  |
             |          \                /     |
              \                               /

            */

            elemat(tvi  ,tuipp) -= coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi);
            elemat(tvip ,tuipp) -= coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi);

          } // vi
        } // ui


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_nu_afgdt_tauM_tauM=fac*visceff*afgdt*tauM*tauM;

          double temp[2];

          temp[0]=fac_nu_afgdt_tauM_tauM*resM_(0);
          temp[1]=fac_nu_afgdt_tauM_tauM*resM_(1);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            double rowtemp[2][2];

            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);

            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui   =3*ui;
              const int tuip  =tui+1;

              /*
                   factor: + 2.0 * visc * alphaF * gamma * dt * tauM * tauM

                    /                                                \
                   |          / /             /    \  \         \     |
                   |  resM , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /
              */

              elemat(tvi  ,tui  ) += viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
              elemat(tvi  ,tuip ) += derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

              elemat(tvip ,tui  ) += viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
              elemat(tvip ,tuip ) += derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
            } // ui
          } //vi
        } // end hoel
      } // end if reynolds stab


      //---------------------------------------------------------------
      /*
             GALERKIN PART, INERTIA, CONVECTION AND VISCOUS TERMS
	                  QUASISTATIC FORMULATION

          ---------------------------------------------------------------

          inertia term (intermediate) + convection (intermediate)

                /          \                   /                          \
               |            |                 |  / n+af       \            |
      +alphaM *|  Dacc , v  |-alphaF*gamma*dt*| | u    o nabla | Dacc , v  |
               |            |                 |  \ G          /            |
                \          /                   \                          /


                                   /                            \
                                  |          / n+af        \     |
                 -alphaF*gamma*dt |  Dacc , | u     o nabla | v  |
                                  |          \             /     |
                                   \                            /



|	  convection (intermediate)
|
N                                  /                            \
E                                 |   n+af    /            \     |
W                -alphaF*gamma*dt |  u     , | Dacc o nabla | v  |
T                                 |           \            /     |
O                                  \                            /
N


	  viscous term (intermediate), factor: +2*nu*alphaF*gamma*dt

                                   /                          \
                                  |       /    \         / \   |
            +2*nu*alphaF*gamma*dt |  eps | Dacc | , eps | v |  |
                                  |       \    /         \ /   |
                                   \                          /

          pressure

                                   /                \
	                          |                  |
                        -gamma*dt*|  Dp , nabla o v  |
                                  |                  |
                                   \                /

          continuity
                                   /                  \
                                  |                    |
                        gamma*dt* | nabla o Dacc  , q  |
                                  |                    |
                                   \                  /
      */
      //---------------------------------------------------------------


      /*---------------------------------------------------------------

                     SUPG PART, INERTIA AND CONVECTION TERMS
                REYNOLDS SUPG TYPE LINEARISATIONS, IF NECESSARY
                       QUASISTATIC FORMULATION (IF ACTIVE)

        ---------------------------------------------------------------

          inertia and convection, factor: +alphaM*tauM

                               /                                        \
                              |          / / n+af  ~n+af \         \     |
                 +alphaM*tauM*|  Dacc , | | c    + u      | o nabla | v  |+
                              |          \ \             /         /     |
                               \                                        /


                               /                                                           \
                              |    / n+af        \          / / n+af  ~n+af \         \     |
        +alphaF*gamma*dt*tauM*|   | c     o nabla | Dacc , | | c    + u      | o nabla | v  |
                              |    \             /          \ \             /         /     |
                               \                                                           /


                               /                                            \
                              |              / / n+af  ~n+af \         \     |
               +tauM*gamma*dt*|  nabla Dp , | | c    + u      | o nabla | v  |
                              |              \ \             /         /     |
                               \                                            /


                               /                                                          \
                              |                 /     \    / / n+af  ~n+af \         \     |
     -nu*alphaF*gamma*dt*tauM*|  2*nabla o eps | Dacc  |, | | c    + u      | o nabla | v  |
                              |                 \     /    \ \             /         /     |
                               \                                                          /



|         linearised convective term in residual
|
N                              /                                                           \
E                             |    /            \   n+af    / / n+af  ~n+af \         \     |
W       +alphaF*gamma*dt*tauM |   | Dacc o nabla | u     , | | c    + u      | o nabla | v  |
T                             |    \            /           \ \             /         /     |
O                              \                                                           /
N


|	  linearisation of testfunction
|
N                              /                            \
E                             |   n+af    /            \     |
W       +alphaF*gamma*dt*tauM*|  r     , | Dacc o nabla | v  |
T                             |   M       \            /     |
O                              \                            /
N
      */
      //---------------------------------------------------------------


      //---------------------------------------------------------------
      /*
	           LEAST SQUARES CONTINUITY STABILISATION PART,
	              QUASISTATIC FORMULATION (IF ACTIVE)

        ---------------------------------------------------------------

          factor: +gamma*dt*tauC

                         /                          \
                        |                            |
                        | nabla o Dacc  , nabla o v  |
                        |                            |
                         \                          /
      */


      const double fac_afgdt         = fac*afgdt;
      const double fac_visceff_afgdt = fac_afgdt*visceff;
      const double fac_gamma_dt      = fac*gamma*dt;
      const double fac_alphaM        = fac*alphaM;


      const double fac_afgdt_velintaf_x=fac_afgdt*velintaf_(0);
      const double fac_afgdt_velintaf_y=fac_afgdt*velintaf_(1);

      // supg and cstab conservative
      const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tui  =3*ui;
        const int tuip =tui+1;

        /* GALERKIN inertia term (intermediate) + convection, mesh velocity (intermediate) */
        const double inertia_and_gridconv_ui
          = fac_alphaM*funct_(ui)-fac_afgdt*conv_u_G_af_(ui);

        /* SUPG stabilisation --- inertia and convection */
        const double inertia_and_conv
          = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);

        // convection GALERKIN and diagonal parts of viscous term (intermediate)
        const double convection_and_viscous_x=fac_visceff_afgdt*derxy_(0,ui)-fac_afgdt_velintaf_x*funct_(ui);
        const double convection_and_viscous_y=fac_visceff_afgdt*derxy_(1,ui)-fac_afgdt_velintaf_y*funct_(ui);

        // viscous GALERKIN term
        const double viscous_x=fac_visceff_afgdt*derxy_(0,ui);
        const double viscous_y=fac_visceff_afgdt*derxy_(1,ui);

        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          const double sum =
            inertia_and_gridconv_ui*funct_(vi)
            +
            inertia_and_conv*conv_c_plus_svel_af_(vi)
            +
            convection_and_viscous_x*derxy_(0,vi)
            +
            convection_and_viscous_y*derxy_(1,vi);

          /* adding GALERKIN convection, convective linearisation (intermediate), viscous and cstab */

          elemat(tvi ,tui ) += sum+(fac_gamma_dt_tauC_derxy_x_ui             + viscous_x)*derxy_(0,vi);
          elemat(tvi ,tuip) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi)+(viscous_x)*derxy_(1,vi);
          elemat(tvip,tui ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi)+(viscous_y)*derxy_(0,vi);
          elemat(tvip,tuip) += sum+(fac_gamma_dt_tauC_derxy_y_ui             + viscous_y)*derxy_(1,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp=3*ui+2;

        const double fac_gamma_dt_derxy_0_ui = fac_gamma_dt*derxy_(0,ui);
        const double fac_gamma_dt_derxy_1_ui = fac_gamma_dt*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          int tvi =vi*3;

          /* SUPG stabilisation --- pressure    */
          /* factor: +tauM, rescaled by gamma*dt*/

          elemat(tvi++,tuipp) += fac_gamma_dt_derxy_0_ui*conv_c_plus_svel_af_(vi);
          elemat(tvi++,tuipp) += fac_gamma_dt_derxy_1_ui*conv_c_plus_svel_af_(vi);

        } // vi
      } // ui

      if (higher_order_ele && newton!=INPAR::FLUID::minimal)
      {
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui  =ui*3 ;
          const int tuip =tui+1;

          /* SUPG stabilisation --- diffusion   */
          /* factor: -nu*alphaF*gamma*dt*tauM   */

          const double fac_visceff_afgdt_viscs2_0_ui=fac_visceff_afgdt*viscs2_(0,ui);
          const double fac_visceff_afgdt_viscs2_1_ui=fac_visceff_afgdt*viscs2_(1,ui);
          const double fac_visceff_afgdt_derxy2_2_ui=fac_visceff_afgdt*derxy2_(3,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int tvi  =vi*3 ;
            int tvip =tvi+1;

            elemat(tvi ,tui ) -= fac_visceff_afgdt_viscs2_0_ui*conv_c_plus_svel_af_(vi);
            elemat(tvi ,tuip) -= fac_visceff_afgdt_derxy2_2_ui*conv_c_plus_svel_af_(vi);
            elemat(tvip,tui ) -= fac_visceff_afgdt_derxy2_2_ui*conv_c_plus_svel_af_(vi);
            elemat(tvip,tuip) -= fac_visceff_afgdt_viscs2_1_ui*conv_c_plus_svel_af_(vi);
          } // vi
        } // ui
      }// end higher_order_ele and linearisation of viscous term

      //---------------------------------------------------------------
      //
      //                  GALERKIN AND SUPG PART
      //        REYNOLDS LINEARISATIONS VIA SUPG TESTFUNCTION
      //    REACTIVE TYPE LINEARISATIONS, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------

      if (newton==INPAR::FLUID::Newton)
      {
        double temp[2][2];
        double testlin[2];

        /* for linearisation of testfunction (SUPG) and reactive GALERKIN part */
        testlin[0]=supg_active_tauM*resM_(0)-velintaf_(0);
        testlin[1]=supg_active_tauM*resM_(1)-velintaf_(1);

        // loop rows
        for (int vi=0; vi<iel; ++vi)
        {
          int tvi   =vi*3;
          int tvip  =tvi+1;

          /*  add linearised convective term in residual (SUPG), reactive
              GALERKIN part and linearisation of testfunction (SUPG) */
          temp[0][0]=fac_afgdt*(testlin[0]*derxy_(0,vi)+vderxyaf_(0,0)*conv_c_plus_svel_af_(vi));
          temp[0][1]=fac_afgdt*(testlin[0]*derxy_(1,vi)+vderxyaf_(0,1)*conv_c_plus_svel_af_(vi));
          temp[1][0]=fac_afgdt*(testlin[1]*derxy_(0,vi)+vderxyaf_(1,0)*conv_c_plus_svel_af_(vi));
          temp[1][1]=fac_afgdt*(testlin[1]*derxy_(1,vi)+vderxyaf_(1,1)*conv_c_plus_svel_af_(vi));

          // loop columns
          for (int ui=0; ui<iel; ++ui)
          {
            int tui=3*ui;

            elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
            elemat(tvip ,tui++) += temp[1][0]*funct_(ui);
            elemat(tvi  ,tui  ) += temp[0][1]*funct_(ui);
            elemat(tvip ,tui  ) += temp[1][1]*funct_(ui);
          } // ui
        } // vi
      } // end newton

      //---------------------------------------------------------------
      //
      //      GALERKIN PART, CONTINUITY AND PRESSURE PART
      //                QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------

      for (int vi=0; vi<iel; ++vi)  // loop rows
      {
	const int tvi    =3*vi;
	const int tvip   =tvi+1;

	const double fac_gamma_dt_derxy_0_vi = fac_gamma_dt*derxy_(0,vi);
	const double fac_gamma_dt_derxy_1_vi = fac_gamma_dt*derxy_(1,vi);

	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int tuipp=3*ui+2;

	  /* GALERKIN pressure (implicit, rescaled to keep symmetry) */
	  /*  factor: -1, rescaled by gamma*dt */

	  elemat(tvi  ,tuipp) -= fac_gamma_dt_derxy_0_vi*funct_(ui);
	  elemat(tvip ,tuipp) -= fac_gamma_dt_derxy_1_vi*funct_(ui);

	  /* continuity equation (implicit, transposed of above equation) */
	  /*  factor: +gamma*dt */

	  elemat(tuipp,tvi  ) += fac_gamma_dt_derxy_0_vi*funct_(ui);
	  elemat(tuipp,tvip ) += fac_gamma_dt_derxy_1_vi*funct_(ui);
	} // ui
      } // vi

      //---------------------------------------------------------------
      //
      //             PSPG PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if(pspg == INPAR::FLUID::pstab_use_pspg)
      {
	const double fac_tauMp                   = fac*tauMp;
	const double fac_alphaM_tauMp            = fac_tauMp*alphaM;
	const double fac_gamma_dt_tauMp          = fac_tauMp*gamma*dt;
	const double fac_afgdt_tauMp             = fac_tauMp*afgdt;

	if (higher_order_ele && newton!=INPAR::FLUID::minimal)
	{
	  const double fac_visceff_afgdt_tauMp
	    =
	    fac*visceff*afgdt*tauMp;

	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int tui  =ui*3;
	    const int tuip =tui+1;

	    /* pressure stabilisation --- diffusion  */

	    /* factor: -nu*alphaF*gamma*dt*tauMp

                    /                                  \
                   |                 /    \             |
                   |  2*nabla o eps | Dacc | , nabla q  |
                   |                 \    /             |
                    \                                  /
	    */

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
	    const double fac_tauMp_inertia_and_conv
	      =
	      fac_alphaM_tauMp*funct_(ui)+fac_afgdt_tauMp*conv_c_af_(ui);

	    const double pspg_diffusion_inertia_convect_0_ui
	      =
	      fac_visceff_afgdt_tauMp*viscs2_(0,ui)-fac_tauMp_inertia_and_conv;
	    const double pspg_diffusion_inertia_convect_1_ui
	      =
	      fac_visceff_afgdt_tauMp*viscs2_(1,ui)-fac_tauMp_inertia_and_conv;

	    const double fac_visceff_afgdt_tauMp_derxy2_2_ui=fac_visceff_afgdt_tauMp*derxy2_(2,ui);

	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int tvipp=vi*3+2;

	      elemat(tvipp,tui  ) -=
		pspg_diffusion_inertia_convect_0_ui*derxy_(0,vi)
		+
		fac_visceff_afgdt_tauMp_derxy2_2_ui*derxy_(1,vi);
	      elemat(tvipp,tuip ) -=
		fac_visceff_afgdt_tauMp_derxy2_2_ui*derxy_(0,vi)
		+
		pspg_diffusion_inertia_convect_1_ui*derxy_(1,vi);
	    } // vi
	  } // ui
	} // this is a higher order element and linearisation is not minimal
	else
	{ // either this ain't a higher order element or a
	  // linearisation of the viscous term is not necessary
	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int tui  =ui*3 ;
	    const int tuip =tui+1;

	    const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));

	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int tvipp=vi*3+2;

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

	      elemat(tvipp,tui ) += fac_tauMp_inertia_and_conv*derxy_(0,vi) ;
	      elemat(tvipp,tuip) += fac_tauMp_inertia_and_conv*derxy_(1,vi) ;
	    } // vi
	  } // ui
	} // no linearisation of viscous part of residual is
	  // performed for pspg stabilisation cause either this
  	  // ain't a higher order element or a linearisation of
	  // the viscous term is not necessary

	if (newton==INPAR::FLUID::Newton)
	{
	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    int vidx = vi*3 + 2;
	    double v1 = derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0);
	    double v2 = derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1);
	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const double fac_afgdt_tauMp_funct_ui = fac_afgdt_tauMp*funct_(ui);
	      int uidx = ui*3;

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
	    } // ui
	  } // vi
	} // end newton

	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int tuipp=ui*3+2;

	  const double fac_gamma_dt_tauMp_derxy_0_ui=fac_gamma_dt_tauMp*derxy_(0,ui);
	  const double fac_gamma_dt_tauMp_derxy_1_ui=fac_gamma_dt_tauMp*derxy_(1,ui);

	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    /* pressure stabilisation --- rescaled pressure   */

	    /* factor: +tauMp, rescaled by gamma*dt

                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
	    */

	    elemat(vi*3+2,tuipp) +=
	      fac_gamma_dt_tauMp_derxy_0_ui*derxy_(0,vi)
	      +
	      fac_gamma_dt_tauMp_derxy_1_ui*derxy_(1,vi);

	  } // vi
	} // ui
      } // end pspg


      //---------------------------------------------------------------
      //
      //      VISCOUS STABILISATION PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if (higher_order_ele)
      {
	if((vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_usfem)&&higher_order_ele)
	{
	  const double fac_visc_tauMp_gamma_dt      = vstabfac*fac*visc*tauMp*gamma*dt;
	  const double fac_visc_afgdt_tauMp         = vstabfac*fac*visc*afgdt*tauMp;
	  const double fac_visc_alphaM_tauMp        = vstabfac*fac*visc*alphaM*tauMp;
	  const double fac_visceff_visc_afgdt_tauMp = vstabfac*fac*visceff*visc*afgdt*tauMp;

	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int tui    = ui*3;
	    const int tuip   = tui+1;
	    const int tuipp  = tui+2;

	    const double acc_conv=(fac_visc_alphaM_tauMp*funct_(ui)
				   +
				   fac_visc_afgdt_tauMp*conv_c_af_(ui));

	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int tvi   = vi*3;
	      const int tvip  = tvi+1;

	      /* viscous stabilisation --- inertia     */

	      /* factor: +(-)alphaM*tauMp*nu

                    /                      \
                   |                        |
                   |  Dacc , 2*div eps (v)  |
                   |                        |
                    \                      /
	      */
	      /* viscous stabilisation --- convection */

	      /*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                       /                                    \
                      |  / n+af       \                      |
                      | | c    o nabla | Dacc, 2*div eps (v) |
                      |  \            /                      |
                       \                                    /

	      */

	      elemat(tvi  ,tui  ) += acc_conv*viscs2_(0,vi);
	      elemat(tvi  ,tuip ) += acc_conv*derxy2_(2,vi);
	      elemat(tvip ,tui  ) += acc_conv*derxy2_(2,vi);
	      elemat(tvip ,tuip ) += acc_conv*viscs2_(1,vi);

	      /* viscous stabilisation --- diffusion  */

	      /* factor: -(+)nu*nu*alphaF*gamma*dt*tauMp

                    /                                       \
                   |                 /    \                  |
                   |  2*nabla o eps | Dacc | , 2*div eps (v) |
                   |                 \    /                  |
                    \                                       /
	      */
	      elemat(tvi  ,tui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*viscs2_(0,vi)
				      +
				      derxy2_(2,ui)*derxy2_(2,vi));
	      elemat(tvi  ,tuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(2,ui)
				      +
				      derxy2_(2,vi)*viscs2_(1,ui));
	      elemat(tvip ,tui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(2,vi)
				      +
				      derxy2_(2,ui)*viscs2_(1,vi));
	      elemat(tvip ,tuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(2,ui)*derxy2_(2,vi)
				      +
				      viscs2_(1,ui)*viscs2_(1,vi));

	      /* viscous stabilisation --- pressure   */

	      /* factor: +(-)tauMp*nu, rescaled by gamma*dt

                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
	      */
	      elemat(tvi  ,tuipp) += fac_visc_tauMp_gamma_dt*
                                     (derxy_(0,ui)*viscs2_(0,vi)
				      +
				      derxy_(1,ui)*derxy2_(2,vi));
	      elemat(tvip ,tuipp) += fac_visc_tauMp_gamma_dt*
                                     (derxy_(0,ui)*derxy2_(2,vi)
				      +
				      derxy_(1,ui)*viscs2_(1,vi));
	    } // vi
	  } // ui

	  if (newton==INPAR::FLUID::Newton)
	  {
	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int tui    = ui*3;
	      const int tuip   = tui+1;

	      const double fac_visc_afgdt_tauMp_funct_ui=fac_visc_afgdt_tauMp*funct_(ui);

	      for (int vi=0; vi<iel; ++vi)  // loop rows
	      {
		const int tvi   = vi*3;
		const int tvip  = tvi+1;

		/* viscous stabilisation --- convection */

		/*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /

		*/
		elemat(tvi  ,tui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,0)
					+
					derxy2_(2,vi)*vderxyaf_(1,0));
		elemat(tvi  ,tuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,1)
					+
					derxy2_(2,vi)*vderxyaf_(1,1));
		elemat(tvip ,tui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(2,vi)*vderxyaf_(0,0)
					+
					viscs2_(1,vi)*vderxyaf_(1,0));
		elemat(tvip ,tuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(2,vi)*vderxyaf_(0,1)
					+
					viscs2_(1,vi)*vderxyaf_(1,1));
	      } // end vi
	    } // end ui
	  } // end newton
	} // endif (a)gls
      } // end hoel

      //---------------------------------------------------------------
      //
      //               QUASISTATIC STABILISATION PART
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == INPAR::FLUID::cross_stress_stab)
      {
        const double fac_tauM         =fac*tauM;
        const double fac_tauM_alphaM  =fac_tauM*alphaM;
        const double fac_tauM_afgdt   =fac_tauM*afgdt;
        const double fac_tauM_gdt     =fac_tauM*gamma*dt;

        double fac_tauM_alphaM_velintaf[2];
        fac_tauM_alphaM_velintaf[0]=fac_tauM_alphaM*velintaf_(0);
        fac_tauM_alphaM_velintaf[1]=fac_tauM_alphaM*velintaf_(1);

        double fac_tauM_afgdt_velintaf[2];
        fac_tauM_afgdt_velintaf[0]=fac_tauM_afgdt*velintaf_(0);
        fac_tauM_afgdt_velintaf[1]=fac_tauM_afgdt*velintaf_(1);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;


          /*
                  /                         \
                 |    n+af                   |
                 |   u     , resM o nabla v  |
                 |                           |
                  \                         /
                    +----+
                      ^
                      |
                      +------ linearisation of this part
          */


          /* factor: tauM*afgdt

                  /                         \
                 |                           |
                 |   Dacc  , resM o nabla v  |
                 |                           |
                  \                         /
          */
          const double fac_tauM_afgdt_conv_resM_vi=fac_tauM_afgdt*conv_resM_(vi);

          double aux[2];

          /*
                  /                         \
                 |    n+af                   |
                 |   u     , resM o nabla v  |
                 |                           |
                  \                         /
                            +----+
                               ^
                               |
                               +------ linearisation of second part
          */

          /* factor: tauM*afgdt

                  /                                               \
                 |    n+af    / /            \   n+af \            |
                 |   u     , | | Dacc o nabla | u      | o nabla v |
                 |            \ \            /        /            |
                  \                                               /
          */
          aux[0]=vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi);
          aux[1]=vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi);

          double temp[2][2];

          /* factor: tauM*alpha_M

                  /                         \
                 |    n+af                   |
                 |   u     , Dacc o nabla v  |
                 |                           |
                  \                         /
          */
          temp[0][0]=fac_tauM_alphaM_velintaf[0]*derxy_(0,vi)+fac_tauM_afgdt_velintaf[0]*aux[0]+fac_tauM_afgdt_conv_resM_vi;
          temp[0][1]=fac_tauM_alphaM_velintaf[0]*derxy_(1,vi)+fac_tauM_afgdt_velintaf[0]*aux[1];

          temp[1][0]=fac_tauM_alphaM_velintaf[1]*derxy_(0,vi)+fac_tauM_afgdt_velintaf[1]*aux[0];
          temp[1][1]=fac_tauM_alphaM_velintaf[1]*derxy_(1,vi)+fac_tauM_afgdt_velintaf[1]*aux[1]+fac_tauM_afgdt_conv_resM_vi;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

	    elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
	    elemat(tvi  ,tuip ) += temp[0][1]*funct_(ui);

	    elemat(tvip ,tui  ) += temp[1][0]*funct_(ui);
	    elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);
	  } // ui
	} // vi

        double fac_tauM_gdt_velintaf[2];
        fac_tauM_gdt_velintaf[0]=fac_tauM_gdt*velintaf_(0);
        fac_tauM_gdt_velintaf[1]=fac_tauM_gdt*velintaf_(1);

	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int tuipp =3*ui+2;

	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int tvi   =3*vi;
	    const int tvip  =tvi+1;

            /* factor: tauM, rescaled by gamma*dt

                         /                                  \
                        |    n+af    /          \            |
                        |   u     , |  nabla Dp  | o nabla v |
                        |            \          /            |
                         \                                  /
            */
	    const double aux=derxy_(0,vi)*derxy_(0,ui)+derxy_(1,vi)*derxy_(1,ui);

	    elemat(tvi  ,tuipp) += fac_tauM_gdt_velintaf[0]*aux;
	    elemat(tvip ,tuipp) += fac_tauM_gdt_velintaf[1]*aux;
	  } // vi
	} // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          /* factor: tauM*afgdt

                  /                                               \
                 |    n+af    / /  n+af       \       \            |
                 |   u     , | |  u    o nabla | Dacc  | o nabla v |
                 |            \ \             /       /            |
                  \                                               /
          */
          double temp[2][2];

          temp[0][0]=fac_tauM_afgdt_velintaf[0]*derxy_(0,vi);
          temp[0][1]=fac_tauM_afgdt_velintaf[0]*derxy_(1,vi);

          temp[1][0]=fac_tauM_afgdt_velintaf[1]*derxy_(0,vi);
          temp[1][1]=fac_tauM_afgdt_velintaf[1]*derxy_(1,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

	    elemat(tvi ,tui ) += temp[0][0]*conv_c_af_(ui);
	    elemat(tvi ,tuip) += temp[0][1]*conv_c_af_(ui);

	    elemat(tvip,tui ) += temp[1][0]*conv_c_af_(ui);
	    elemat(tvip,tuip) += temp[1][1]*conv_c_af_(ui);
	  } // ui
	} // vi

	if (higher_order_ele && newton!=INPAR::FLUID::minimal)
	{
	  const double fac_nu_afgdt_tauM=fac*visceff*afgdt*tauM;

	  double temp[2];

	  temp[0]=fac_nu_afgdt_tauM*velintaf_(0);
	  temp[1]=fac_nu_afgdt_tauM*velintaf_(1);

	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int tvi   =3*vi;
	    const int tvip  =tvi+1;

	    double rowtemp[2][2];

	    rowtemp[0][0]=temp[0]*derxy_(0,vi);
	    rowtemp[0][1]=temp[0]*derxy_(1,vi);

	    rowtemp[1][0]=temp[1]*derxy_(0,vi);
	    rowtemp[1][1]=temp[1]*derxy_(1,vi);

	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int tui   =3*ui;
	      const int tuip  =tui+1;

	      /*
		 factor: - 2.0 * visc * alphaF * gamma * dt * tauM

                    /                                                \
                   |   n+af   / /             /    \  \         \     |
                   |  u    , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /
	      */

	      elemat(tvi ,tui ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
	      elemat(tvi ,tuip) -= derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

	      elemat(tvip,tui ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
	      elemat(tvip,tuip) -= derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
	    } //  ui
	  } // vi
	} // hoel
      } // cross

    } // compute_elemat

    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //
    //          RIGHT HAND SIDE, QUASISTATIC SUBGRID SCALES
    //
    //---------------------------------------------------------------
    //---------------------------------------------------------------

   /* inertia, convective and dead load terms -- all tested
       against shapefunctions, as well as cross terms            */
    /*

                /             \
               |     n+am      |
              -|  acc     , v  |
               |               |
                \             /


                /                             \
               |  / n+af       \    n+af       |
              +| | u    o nabla |  u      , v  |
               |  \ G          /               |
                \                             /

                /           \
               |   n+af      |
              +|  f     , v  |
               |             |
                \           /

    */

    double fac_inertia_gridconv_and_bodyforce[2];
    fac_inertia_gridconv_and_bodyforce[0] = fac*(accintam_(0)-convu_G_af_old_(0)-bodyforceaf_(0));
    fac_inertia_gridconv_and_bodyforce[1] = fac*(accintam_(1)-convu_G_af_old_(1)-bodyforceaf_(1));

    /*
      pressure, partially integrated convective term and viscous
      term combined in viscous_conv_and_pres
      cross and reynolds stabilisation are combined with the
      same testfunctions (of derivative type).
      continuity stabilisation adds a small-scale pressure
    */

    /*
       factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /

    */
    /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
    */

    const double fac_prenp   = fac*prenp_-fac*tauC*divunp_;


    /*
      factor: +2*nu
        /                            \     /                              \
       |       / n+af \         / \   |   |       / n+af \           / \   |
       |  eps | u      | , eps | v |  | = |  eps | u      | , nabla | v |  |
       |       \      /         \ /   |   |       \      /           \ /   |
        \                            /     \                              /

    */

    const double visceff_fac = visceff*fac;

    double viscous_conv_and_pres[4];
    viscous_conv_and_pres[0]=visceff_fac*vderxyaf_(0,0)*2.0-fac_prenp;
    viscous_conv_and_pres[1]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_conv_and_pres[2]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_conv_and_pres[3]=visceff_fac*vderxyaf_(1,1)*2.0-fac_prenp;

    /*
          factor: -1.0

               /                                       \
              |   / n+af \     / n+af \           / \   |
              |  | u      | X | u      | , nabla | v |  |
              |   \      /     \      /           \ /   |
               \                                       /
    */

    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs || cross == INPAR::FLUID::cross_stress_stab)
    {
      const double fac_tauM = fac*tauM;

      /* factor: -tauM

                  /                             \
                 |     n+af                      |
                 |  ( u     x resM ) ,  nabla v  |
                 |     (i)                       |
                  \                             /
      */
      viscous_conv_and_pres[0]-=velintaf_(0)*(-fac_tauM*resM_(0)+fac*velintaf_(0));
      viscous_conv_and_pres[1]-=velintaf_(0)*(-fac_tauM*resM_(1)+fac*velintaf_(1));
      viscous_conv_and_pres[2]-=velintaf_(1)*(-fac_tauM*resM_(0)+fac*velintaf_(0));
      viscous_conv_and_pres[3]-=velintaf_(1)*(-fac_tauM*resM_(1)+fac*velintaf_(1));
    }
    else
    {
      viscous_conv_and_pres[0]-=velintaf_(0)*velintaf_(0)*fac;
      viscous_conv_and_pres[1]-=velintaf_(0)*velintaf_(1)*fac;
      viscous_conv_and_pres[2]-=velintaf_(1)*velintaf_(0)*fac;
      viscous_conv_and_pres[3]-=velintaf_(1)*velintaf_(1)*fac;
    }

    if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
    {

      /* factor: -tauM*tauM

                  /                             \
                 |                               |
                 |  resM   , ( resM o nabla ) v  |
                 |                               |
                  \                             /
      */
      const double fac_tauM_tauM        = fac*tauM*tauM;
      const double fac_tauM_tauM_resM_0 = fac_tauM_tauM*resM_(0);
      const double fac_tauM_tauM_resM_1 = fac_tauM_tauM*resM_(1);

      viscous_conv_and_pres[0]-=fac_tauM_tauM_resM_0*resM_(0);
      viscous_conv_and_pres[1]-=fac_tauM_tauM_resM_0*resM_(1);
      viscous_conv_and_pres[2]-=fac_tauM_tauM_resM_0*resM_(1);
      viscous_conv_and_pres[3]-=fac_tauM_tauM_resM_1*resM_(1);
    }

    /* continuity equation, factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
    */
    const double fac_divunp  = fac*divunp_;

    for (int vi=0; vi<iel; ++vi) // loop rows  (test functions)
    {
      int tvi=3*vi;
      /* inertia, convective and dead load, cross terms fith
	 funct                                              */
      /* viscous, pressure, reynolds, cstab terms with
	 derxy                                              */

      elevec(tvi++) -=
	fac_inertia_gridconv_and_bodyforce[0]*funct_(vi)
	+
	derxy_(0,vi)*viscous_conv_and_pres[0]
	+
	derxy_(1,vi)*viscous_conv_and_pres[1];
      elevec(tvi++) -=
	fac_inertia_gridconv_and_bodyforce[1]*funct_(vi)
	+
	derxy_(0,vi)*viscous_conv_and_pres[2]
	+
	derxy_(1,vi)*viscous_conv_and_pres[3];

      /* continuity equation */
      elevec(tvi  ) -= fac_divunp*funct_(vi);
    }

    if(pspg == INPAR::FLUID::pstab_use_pspg)
    {
      /*
      pressure stabilisation

      factor: +tauMp

                  /                 \
                 |    n+af           |
                 |  r     , nabla q  |
                 |   M               |
                  \                 /

      */
      const double fac_tauMp = fac*tauMp;

      for (int vi=0; vi<iel; ++vi) // loop rows
      {
	elevec(3*vi+2)-=fac_tauMp*(resM_(0)*derxy_(0,vi)+resM_(1)*derxy_(1,vi));
      }
    } // end pspg

    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      const double fac_tauM = fac*tauM;

      for (int vi=0; vi<iel; ++vi) // loop rows  (test functions)
      {
	int tvi=3*vi;

	const double fac_tauM_conv_c_af_vi = fac_tauM*conv_c_af_(vi);
	/*
	  factor: +tauM

	  SUPG stabilisation


                  /                             \
                 |   n+af    / n+af        \     |
                 |  r     , | c     o nabla | v  |
                 |   M       \             /     |
                  \                             /
	*/

	elevec(tvi++) -= fac_tauM_conv_c_af_vi*resM_(0);
	elevec(tvi  ) -= fac_tauM_conv_c_af_vi*resM_(1);

      } // end loop rows
    } // end supg

    if (higher_order_ele)
    {
      if(vstab != INPAR::FLUID::viscous_stab_none && higher_order_ele)
      {
	const double fac_visc_tauMp = vstabfac * fac*visc*tauMp;

	for (int vi=0; vi<iel; ++vi) // loop rows
	{
	  int tvi=3*vi;
	  /*
	    factor: -(+)tauMp*nu

	    viscous stabilisation --- inertia


                 /                      \
                |   n+af                 |
                |  r    , 2*div eps (v)  |
                |   M                    |
                 \                      /

	  */
	  elevec(tvi++) -= fac_visc_tauMp*
	                   (resM_(0)*viscs2_(0,vi)
			    +
			    resM_(1)*derxy2_(2,vi));
	  elevec(tvi  ) -= fac_visc_tauMp*
                           (resM_(0)*derxy2_(2,vi)
			    +
			    resM_(1)*viscs2_(1,vi));
	} // end loop rows ui
      } // endif (a)gls
    } // end hoel

  } // end loop iquad

  return;
} // Sysmat_cons_qs


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  conservative, time-dependent version                               |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_cons_td_2D(
  Fluid*                                ele,
  std::vector<Epetra_SerialDenseVector>& myknots,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&           elemat,
  LINALG::Matrix<(nsd_+1)*iel,    1>&           elevec,
  const LINALG::Matrix<nsd_,iel>&           edispnp,
  const LINALG::Matrix<nsd_,iel>&           egridvaf,
  const LINALG::Matrix<nsd_,iel>&           evelnp,
  const LINALG::Matrix<iel,1>&           eprenp,
  const LINALG::Matrix<nsd_,iel>&           eaccam,
  const LINALG::Matrix<nsd_,iel>&           evelaf,
  Teuchos::RCP<const MAT::Material>      material,
  const double                           alphaM,
  const double                           alphaF,
  const double                           gamma,
  const double                           dt,
  const double                           time,
  const enum INPAR::FLUID::LinearisationAction newton,
  const bool                             higher_order_ele,
  const enum INPAR::FLUID::Transient         inertia         ,
  const enum INPAR::FLUID::PSPG              pspg            ,
  const enum INPAR::FLUID::SUPG              supg            ,
  const enum INPAR::FLUID::VStab             vstab           ,
  const enum INPAR::FLUID::CStab             cstab           ,
  const enum INPAR::FLUID::CrossStress       cross           ,
  const enum INPAR::FLUID::ReynoldsStress    reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha             whichtau,
  double&                                visceff,
  const bool                             compute_elemat
  )
{

  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt
  //
  const double timealphaF = time-(1-alphaF)*dt;

  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_gls_only_rhs)
  {
    vstabfac = -1.0;
  }

  //------------------------------------------------------------------
  //                    SET ALL ELEMENT DATA
  // o including element geometry (node coordinates)
  // o including dead loads in nodes
  // o including hk, mk, element area
  // o including material viscosity, effective viscosity by
  //   Non-Newtonian fluids
  //------------------------------------------------------------------

  double hk   = 0.0;
  double mk   = 0.0;
  double visc = 0.0;

  SetElementData_2D(ele            ,
                 edispnp        ,
                 evelaf         ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 visceff        );

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  //const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // remember whether the subscale quantities have been allocated an set to zero.
  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero

//  ele->ActivateTDS(intpoints.IP().nquad, nsd_);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {

    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second
    // derivatives in a gausspoint and integration weight including
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives_2D(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint_2D(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            visceff         ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau_2D(whichtau,INPAR::FLUID::subscales_time_dependent,gamma,dt,hk,mk,visceff);

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

    const double tauM   = tau_(0);

    if(cstab == INPAR::FLUID::continuity_stab_none)
    {
      tau_(2)=0.0;
    }

    const double tauC   = tau_(2);

    double supg_active;
    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      supg_active=1.0;
    }
    else
    {
      supg_active=0.0;
    }

    // update estimates for the subscale quantities
    const double facMtau = 1./(alphaM*tauM+afgdt);
//    const double fac1=(alphaM*tauM+gamma*dt*(alphaF-1.0))*facMtau;
//    const double fac2=(dt*tauM*(alphaM-gamma))*facMtau;
//    const double fac3=(gamma*dt*tauM)*facMtau;

    /*-------------------------------------------------------------------*
     *                                                                   *
     *                  update of SUBSCALE VELOCITY                      *
     *             and update of intermediate quantities                 *
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

     compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

  */

//    for(int rr=0;rr<nsd_;++rr)
//    {
//      ele->UpdateSvelnpInOneDirection(
//          fac1       ,
//          fac2       ,
//          fac3       ,
//          resM_(rr)  ,
//          alphaF     ,
//          rr         ,
//          iquad      ,
//          svelnp_(rr),
//          svelaf_(rr));
//    }

    /* the intermediate value of subscale acceleration is not needed to be
     * computed anymore --- we use the governing ODE to replace it ....

             ~ n+am    alphaM     / ~n+1   ~n \    gamma - alphaM    ~ n
            acc     = -------- * |  u    - u   | + -------------- * acc
               (i)    gamma*dt    \  (i)      /         gamma

    */

    // prepare possible modification of convective linearisation for
    // combined reynolds/supg test function
    for(int nn=0;nn<iel;++nn)
    {
      conv_c_plus_svel_af_(nn)=conv_c_af_(nn)*supg_active;
    }

    /*
        This is the operator

                  /~n+af         \
                 | u      o nabla |
                  \   (i)        /

        required for the cross/reynolds stress linearisation

    */
    if(cross    == INPAR::FLUID::cross_stress_stab
       ||
       reynolds == INPAR::FLUID::reynolds_stress_stab)
    {
      for (int rr=0;rr<iel;++rr)
      {
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr)+svelaf_(1)*derxy_(1,rr);
      }

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {
        /* get modified convective linearisation (n+alpha_F,i) at
           integration point takes care of half of the linearisation

                                   +-----  /                   \
                         n+af       \     |  n+af      ~n+af    |   dN
         conv_c_plus_svel_   (x) =   +    | c    (x) + u    (x) | * --- (x)
                                    /     |  j          j       |   dx
                                   +-----  \                   /      j
                                   dim j    +------+   +------+
                                               if         if
                                              supg     reynolds

        */
        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)+=conv_subaf_(nn);
        }
      }
    }

    /* Most recent value for subgrid velocity convective term

                  /~n+af         \   n+af
                 | u      o nabla | u
                  \   (i)        /   (i)
    */
    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs
       ||
       cross == INPAR::FLUID::cross_stress_stab
      )
    {
      for (int rr=0;rr<2;++rr)
      {
        convsubaf_old_(rr) = vderxyaf_(rr,0)*svelaf_(0)+vderxyaf_(rr,1)*svelaf_(1);
      }
    }

    // get convective linearisation (n+alpha_F,i) at integration point
    // (convection by grid velocity)
    //
    //                    +-----
    //         n+af        \      n+af      dN
    // conv_u_G_    (x) =   +    u    (x) * --- (x)
    //                     /      G,j       dx
    //                    +-----              j
    //                    dim j
    //
    if(ele->IsAle())
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=u_G_af_(0)*derxy_(0,nn)+u_G_af_(1)*derxy_(1,nn);
      }
    }
    else
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=0.0;
      }
    }

    /* Convective term  u_G_old * grad u_old: */
    /*
    //     /    n+af        \   n+af
    //    |  u_G     o nabla | u
    //     \                /
    */
    for(int rr=0;rr<2;++rr)
    {
      convu_G_af_old_(rr)=u_G_af_(0)*vderxyaf_(rr,0)+u_G_af_(1)*vderxyaf_(rr,1);
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

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //
    //          SYSTEM MATRIX, TIME DEPENDENT FORMULATION
    //
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    if(compute_elemat)
    {


      // scaling factors for Galerkin 1 terms
      double fac_inertia   =fac*alphaM;

      const double fac_gamma_dt      = fac*gamma*dt;

      //---------------------------------------------------------------
      //
      //              SUBSCALE ACCELERATION PART
      //        RESCALING FACTORS FOR GALERKIN 1 TERMS AND
      //              COMPUTATION OF EXTRA TERMS
      //
      //---------------------------------------------------------------

      if(inertia == INPAR::FLUID::inertia_stab_keep
         ||
         inertia == INPAR::FLUID::inertia_stab_keep_complete)
      {
        // rescale time factors terms affected by inertia stabilisation
        fac_inertia*=afgdt*facMtau;

        // do inertia stabilisation terms which are not scaled
        // Galerkin terms since they are not partially integrated

        const double fac_alphaM_tauM_facMtau = fac*alphaM*tauM*facMtau;

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          const double fac_alphaM_gamma_dt_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*gamma*dt*funct_(vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tuipp  =3*ui+2;
            /* pressure (implicit) */

            /*  factor:
                             alphaM*tauM
                  - ---------------------------, rescaled by gamma*dt
                    alphaM*tauM+alphaF*gamma*dt

                 /               \
                |                 |
                |  nabla Dp ,  v  |
                |                 |
                 \               /
            */
            /* pressure (implicit) */

            elemat(tvi  ,tuipp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(0,ui);
            elemat(tvip ,tuipp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(1,ui);
          } // ui
        } // vi

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;

            /* convective term (intermediate), convective linearisation */
            /*  factor:
                                                 alphaM*tauM
                           alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                  /                          \
                 |   / n+af       \           |
               - |  | c    o nabla | Dacc , v |
                 |   \            /           |
                  \                          /

            */

              elemat(tvi ,tui ) -= afgdt*fac_alphaM_tauM_facMtau*conv_c_af_(ui)*funct_(vi);
              elemat(tvip,tuip) -= afgdt*fac_alphaM_tauM_facMtau*conv_c_af_(ui)*funct_(vi);
          }
        }
        if(newton==INPAR::FLUID::Newton)
        {
          double temp[2][2];

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi    =3*vi;
            const int tvip   =tvi+1;

            const double aux=afgdt*fac_alphaM_tauM_facMtau*funct_(vi);

            temp[0][0]=aux*vderxyaf_(0,0);
            temp[1][0]=aux*vderxyaf_(0,1);
            temp[0][1]=aux*vderxyaf_(1,0);
            temp[1][1]=aux*vderxyaf_(1,1);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui    =3*ui;
              const int tuip   =tui+1;

              /* convective term (intermediate), reactive part from linearisation */
              /*  factor:
                                                 alphaM*tauM
                           alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                  /                          \
                 |   /            \   n+af    |
               - |  | Dacc o nabla | u    , v |
                 |   \            /           |
                  \                          /

              */


              elemat(tvi ,tui ) -= temp[0][0]*funct_(ui);
              elemat(tvi ,tuip) -= temp[1][0]*funct_(ui);
              elemat(tvip,tui ) -= temp[0][1]*funct_(ui);
              elemat(tvip,tuip) -= temp[1][1]*funct_(ui);
            }
          }
        }

        if(higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_alphaM_tauM_facMtau
            =
            fac*visceff*afgdt*alphaM*tauM*facMtau;

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi    =3*vi;
            const int tvip   =tvi+1;

            const double fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
              =
              fac_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui    =3*ui;
              const int tuip   =tui+1;

              /* viscous term (intermediate) */
              /*  factor:
                                                 alphaM*tauM
                        nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt


                  /                           \
                 |                 /    \      |
                 |  2*nabla o eps | Dacc | , v |
                 |                 \    /      |
                  \                           /

              */
              const double a = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(2,ui);

              elemat(tvi ,tui ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,ui);
              elemat(tvi ,tuip) += a;
              elemat(tvip,tui ) += a;
              elemat(tvip,tuip) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,ui);
            } // ui
          } // vi
        } // end higher order element and linearisation of linear terms not supressed

        if(inertia == INPAR::FLUID::inertia_stab_keep_complete)
        {

          /*
                                  immediately enters the matrix
                                  |
                                  v
                               +--------------+
                               |              |
                                /            \
                      1.0      |  ~n+af       |
                 - --------- * |  u     ,  v  |
                        n+af   |   (i)        |
                   tau_M        \            /

                   |       |
                   +-------+
                       ^
                       |
                       consider linearisation of this expression

          */
          const double norm = sqrt(velintaf_(0)*velintaf_(0)+velintaf_(1)*velintaf_(1));

          // normed velocity at element center (we use the copy for safety reasons!)
          if (norm>=1e-6)
          {
            for (int rr=0;rr<2;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=velintaf_(rr)/norm;
            }
          }
          else
          {
            normed_velintaf_(0) = 0.0;
            normed_velintaf_(1) = 0.0;
          }

          double temp=0.0;
          if(whichtau==INPAR::FLUID::codina)
          {
            /*
                                                  || n+af||
                       1.0           visc         ||u    ||
                    --------- = CI * ---- + CII * ---------
                         n+af           2
                    tau_M             hk             hk


                    where CII=2.0/mk
            */

            temp=fac*afgdt/hk*2.0/mk;
          }
          else if(whichtau==INPAR::FLUID::smoothed_franca_barrenechea_valentin_wall)
          {
            /*
                                  -x   '       -x
                    using f(x)=x+e  , f (x)=1-e


                                                +-                                -+
                                                |          / || n+af||          \  |
                       1.0      4.0 * visceff   |         |  ||u    || * hk * mk | |
                    --------- = ------------- * | 1.0 + f |  ------------------- | |
                         n+af           2       |         |                      | |
                    tau_M         mk* hk        |          \    2.0 * visceff   /  |
                                                +-                                -+

            */

            temp=fac*afgdt/hk*2.0*(1-exp(-1.0*(norm*hk/visceff)*(mk/2.0)));


          }
          else if(whichtau==INPAR::FLUID::franca_barrenechea_valentin_wall)
          {

            /*
                                             +-                                  -+
                                             |            / || n+af||          \  |
                       1.0      4.0 * visc   |           |  ||u    || * hk * mk | |
                    --------- = ---------- * | 1.0 + max |  ------------------- | |
                         n+af           2    |           |                      | |
                    tau_M         mk* hk     |            \    2.0 * visceff   /  |
                                             +-                                  -+

            */

            if((norm*hk/visceff)*(mk/2.0)>1)
            {
              temp=fac*afgdt/hk*2.0;
            }
          }
          else
          {
            dserror("There's no linearisation of 1/tau available for this tau definition\n");
          }

          /*
                        || n+af||             n+af
                      d ||u    ||            u    * Dacc
                      ----------- = afgdt *  -----------
                                              || n+af||
                        d Dacc                ||u    ||

          */

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi    =3*vi;
            const int tvip   =tvi+1;

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int tui  =3*ui;
              const int tuip =tui+1;

              elemat(tvi  ,tui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(tvi  ,tuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(0);

              elemat(tvip ,tui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(tvip ,tuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(1);
            } // ui
          } // vi
        } // end linearisation of 1/tauM
      } // extra terms for inertia stab

      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //
      //      GALERKIN PART 1 (INERTIA, CONVECTION, VISCOUS)
      // GALERKIN PART 2 (REMAINING PRESSURE AND CONTINUITY EXPRESSIONS)
      //
      //               CONTINUITY STABILISATION
      //
      //---------------------------------------------------------------

      /*
        inertia term (intermediate)

                                                 /          \
                         alphaF*gamma*dt        |            |
             alphaM*---------------------------*|  Dacc , v  |
                    alphaM*tauM+alphaF*gamma*dt |            |
                                                 \          /
             |                                 |
             +---------------------------------+
               	            alphaM
           	without inertia stabilisation



                                   /                            \
                                  |          / n+af        \     |
                 -alphaF*gamma*dt |  Dacc , | u     o nabla | v  |
                                  |          \             /     |
                                   \                            /



|	  convection (intermediate)
|
N                                  /                            \
E                                 |   n+af    /            \     |
W                -alphaF*gamma*dt |  u     , | Dacc o nabla | v  |
T                                 |           \            /     |
O                                  \                            /
N


      pressure (implicit)

                                                 /                \
                                                |                  |
                                      -gamma*dt |  Dp , nabla o v  |
                                                |                  |
                                                 \                /

     viscous term (intermediate)


                                                 /                          \
		                                |       /    \         / \   |
                          +2*nu*alphaF*gamma*dt*|  eps | Dacc | , eps | v |  |
                                                |       \    /         \ /   |
                                                 \                          /


     continuity equation (implicit)



                                                 /                  \
                                                |                    |
                                     +gamma*dt* | nabla o Dacc  , q  |
                                                |                    |
                                                 \                  /


      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //               CONTINUITY STABILISATION
      //
      //---------------------------------------------------------------

                                                 /                          \
                                                |                            |
                                +gamma*dt*tauC* | nabla o Dacc  , nabla o v  |
                                                |                            |
                                                 \                          /
                                +-------------+
                               zero for no cstab


      //---------------------------------------------------------------
      //
      //              TIME-DEPENDENT SUBGRID-SCALES
      //
      //                   SUPG STABILISATION
      //            SUPG TYPE REYNOLDS LINEARISATIONS
      //
      //---------------------------------------------------------------
         SUPG stabilisation --- subscale velocity, nonlinear part from testfunction
|
|
N                                       /                            \
E                                      |  ~n+af    /            \     |
W                 alphaF * gamma * dt* |  u     , | Dacc o nabla | v  |
T                                      |   (i)     \            /     |
O                                       \                            /
N

         SUPG stabilisation --- inertia

                              alphaF*gamma*dt
                         --------------------------- * alphaM * tauM *
                         alphaM*tauM+alphaF*gamma*dt


                     /                                        \
                    |          / / n+af  ~n+af \         \     |
                    |  Dacc , | | c    + u      | o nabla | v  |
                    |          \ \             /         /     |
                     \                                        /

        SUPG stabilisation --- convection

                               alphaF*gamma*dt
                         --------------------------- * alphaF * gamma * dt * tauM
                         alphaM*tauM+alphaF*gamma*dt

                     /                                                           \
                    |    / n+af        \          / / n+af  ~n+af \         \     |
                    |   | c     o nabla | Dacc , | | c    + u      | o nabla | v  |
                    |    \             /          \ \             /         /     |
                     \                                                           /

        SUPG stabilisation --- convection

                              alphaF*gamma*dt
|                       --------------------------- * alphaF * gamma * dt * tauM
|                       alphaM*tauM+alphaF*gamma*dt
N
E                   /                                                           \
W                  |    /            \   n+af    / / n+af  ~n+af \         \     |
T                  |   | Dacc o nabla | u     , | | c    + u      | o nabla | v  |
O                  |    \            /           \ \             /         /     |
N                   \                                                           /

        SUPG stabilisation --- pressure

                               alphaF*gamma*dt*tauM
                            ---------------------------, rescaled by gamma*dt
                            alphaM*tauM+alphaF*gamma*dt


                    /                                            \
                   |              / / n+af  ~n+af \         \     |
                   |  nabla Dp , | | c    + u      | o nabla | v  |
                   |              \ \             /         /     |
                    \                                            /

        SUPG stabilisation --- diffusion

                                              alphaF*gamma*dt*tauM
                        nu*alphaF*gamma*dt*---------------------------
                                           alphaM*tauM+alphaF*gamma*dt

                    /                                                          \
                   |  /             /      \     / / n+af  ~n+af \         \    |
                   | | nabla o eps |  Dacc  | , | | c    + u      | o nabla | v |
                   |  \             \      /     \ \             /         /    |
                    \                                                          /
      */

      const double fac_afgdt_afgdt_tauM_facMtau  = fac*afgdt   *afgdt*tauM*facMtau;
      const double fac_gdt_afgdt_tauM_facMtau    = fac*gamma*dt*afgdt*tauM*facMtau;
      const double fac_alphaM_afgdt_tauM_facMtau = fac*alphaM  *afgdt*tauM*facMtau;

      const double fac_afgdt         = fac*afgdt;
      const double fac_visceff_afgdt = fac_afgdt*visceff;

      const double fac_afgdt_velintaf_x=fac_afgdt*velintaf_(0);
      const double fac_afgdt_velintaf_y=fac_afgdt*velintaf_(1);

      // supg and cstab conservative
      const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tui  =3*ui;
        const int tuip =tui+1;

        /* GALERKIN inertia term (intermediate) + convection, mesh velocity (intermediate) */
        const double inertia_and_gridconv_ui = fac_inertia*funct_(ui)-fac_afgdt*conv_u_G_af_(ui);

        /* SUPG stabilisation --- inertia and convection */
        const double supg_inertia_and_conv_ui
          = fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_afgdt_tauM_facMtau*conv_c_af_(ui);

        // convection GALERKIN and diagonal parts of viscous term (intermediate)
        const double convection_and_viscous_x=fac_visceff_afgdt*derxy_(0,ui)-fac_afgdt_velintaf_x*funct_(ui);
        const double convection_and_viscous_y=fac_visceff_afgdt*derxy_(1,ui)-fac_afgdt_velintaf_y*funct_(ui);

        // viscous GALERKIN term
        const double viscous_x=fac_visceff_afgdt*derxy_(0,ui);
        const double viscous_y=fac_visceff_afgdt*derxy_(1,ui);

        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int tvi    =3*vi;
          const int tvip   =tvi+1;

          const double sum =
            inertia_and_gridconv_ui*funct_(vi)
            +
            supg_inertia_and_conv_ui*conv_c_plus_svel_af_(vi)
            +
            convection_and_viscous_x*derxy_(0,vi)
            +
            convection_and_viscous_y*derxy_(1,vi);

          /* adding GALERKIN convection, convective linearisation (intermediate), viscous and cstab */

          elemat(tvi ,tui ) += sum+(fac_gamma_dt_tauC_derxy_x_ui             + viscous_x)*derxy_(0,vi);
          elemat(tvi ,tuip) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi)+(viscous_x)*derxy_(1,vi);
          elemat(tvip,tui ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi)+(viscous_y)*derxy_(0,vi);
          elemat(tvip,tuip) += sum+(fac_gamma_dt_tauC_derxy_y_ui             + viscous_y)*derxy_(1,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp  =3*ui+2;

        const double fac_gamma_dt_funct_ui=fac_gamma_dt*funct_(ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi  =3*vi;
          const int tvip =tvi+1;

          /* GALERKIN pressure   (implicit), rescaled by gamma*dt */
          /* continuity equation (implicit)                       */

          elemat(tvi  ,tuipp) -= fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(tvip ,tuipp) -= fac_gamma_dt_funct_ui*derxy_(1,vi);

          elemat(tuipp,tvi  ) += fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(tuipp,tvip ) += fac_gamma_dt_funct_ui*derxy_(1,vi);
        } // vi
      } // ui

      if (newton==INPAR::FLUID::Newton) // if newton and supg
      {
        const double fac_afgdt_afgdt_tauM_facMtau = fac*afgdt*afgdt*facMtau*tauM;

        // linearisation of SUPG testfunction and GALERKIN reactive part of convection
        double temp[2][2];

        const double fac_afgdt_svelaf_0 = fac*afgdt*supg_active*svelaf_(0)+fac*afgdt*velintaf_(0);
        const double fac_afgdt_svelaf_1 = fac*afgdt*supg_active*svelaf_(1)+fac*afgdt*velintaf_(1);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi  =3*vi;
          const int tvip =tvi+1;

          // SUPG part (reactive part from residual)
          const double scaled_inertia_and_conv_vi
            =
            fac_afgdt_afgdt_tauM_facMtau*conv_c_plus_svel_af_(vi);

          temp[0][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,0)-fac_afgdt_svelaf_0*derxy_(0,vi);
          temp[1][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,1)-fac_afgdt_svelaf_0*derxy_(1,vi);
          temp[0][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,0)-fac_afgdt_svelaf_1*derxy_(0,vi);
          temp[1][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,1)-fac_afgdt_svelaf_1*derxy_(1,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
            elemat(tvi  ,tuip ) += temp[1][0]*funct_(ui);
            elemat(tvip ,tui  ) += temp[0][1]*funct_(ui);
            elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);
          } // ui
        } // vi
      } // end if newton

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int tuipp=3*ui+2;

        const double scaled_gradp_0 = fac_gdt_afgdt_tauM_facMtau*derxy_(0,ui);
        const double scaled_gradp_1 = fac_gdt_afgdt_tauM_facMtau*derxy_(1,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi=3*vi;

          /* SUPG stabilisation --- pressure, rescaled by gamma*dt */
          elemat(tvi  ,tuipp) += scaled_gradp_0*conv_c_plus_svel_af_(vi);
          elemat(tvi+1,tuipp) += scaled_gradp_1*conv_c_plus_svel_af_(vi);
        } // vi
      } // ui

      if(higher_order_ele && newton!=INPAR::FLUID::minimal)
      {
        const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          double coltemp[2][2];

          coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(0,ui);
          coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(2,ui);

          coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(2,ui);
          coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*  SUPG stabilisation, diffusion */
            elemat(tvi  ,tui  ) -= coltemp[0][0]*conv_c_plus_svel_af_(vi);
            elemat(tvi  ,tuip ) -= coltemp[0][1]*conv_c_plus_svel_af_(vi);

            elemat(tvip ,tui  ) -= coltemp[1][0]*conv_c_plus_svel_af_(vi);
            elemat(tvip ,tuip ) -= coltemp[1][1]*conv_c_plus_svel_af_(vi);
          } // vi
        } // ui
      } // hoel

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
      {
        /*
                  /                            \
                 |  ~n+af    ~n+af              |
               - |  u    , ( u     o nabla ) v  |
                 |                              |
                  \                            /
                             +----+
                               ^
                               |
                               linearisation of this expression
        */
        const double fac_alphaM_afgdt_tauM_facMtau=fac*alphaM*afgdt*tauM*facMtau;

        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_x = fac_alphaM_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_y = fac_alphaM_afgdt_tauM_facMtau*svelaf_(1);

        const double fac_afgdt_afgdt_tauM_facMtau =fac*afgdt*afgdt*tauM*facMtau;

        double fac_afgdt_afgdt_tauM_facMtau_svelaf[2];
        fac_afgdt_afgdt_tauM_facMtau_svelaf[0]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[1]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tui   =3*ui;
          const int tuip  =tui+1;

          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui);

          double inertia_and_conv[2];

          inertia_and_conv[0]=fac_afgdt_afgdt_tauM_facMtau_svelaf[0]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_afgdt_tauM_facMtau_svelaf[1]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_y*funct_(ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
               factor: +alphaM * alphaF * gamma * dt * tauM * facMtau

                  /                            \
                 |  ~n+af                       |
                 |  u     , ( Dacc o nabla ) v  |
                 |                              |
                  \                            /

            */

            /*
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / / n+af        \       \         \     |
             |  u     , | | | u     o nabla | Dacc  | o nabla | v  |
             |           \ \ \             /       /         /     |
              \                                                   /

            */

            elemat(tvi  ,tui  ) += inertia_and_conv[0]*derxy_(0,vi);
            elemat(tvi  ,tuip ) += inertia_and_conv[0]*derxy_(1,vi);

            elemat(tvip ,tui  ) += inertia_and_conv[1]*derxy_(0,vi);
            elemat(tvip ,tuip ) += inertia_and_conv[1]*derxy_(1,vi);
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          double temp[3];
          temp[0]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi));
          temp[1]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi));

          double rowtemp[2][2];

          rowtemp[0][0]=svelaf_(0)*temp[0];
          rowtemp[0][1]=svelaf_(0)*temp[1];

          rowtemp[1][0]=svelaf_(1)*temp[0];
          rowtemp[1][1]=svelaf_(1)*temp[1];

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

            /*
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / /            \   n+af \         \     |
             |  u     , | | | Dacc o nabla | u      | o nabla | v  |
             |           \ \ \            /        /         /     |
              \                                                   /

            */

            elemat(tvi  ,tui  ) += funct_(ui)*rowtemp[0][0];
            elemat(tvi  ,tuip ) += funct_(ui)*rowtemp[0][1];

            elemat(tvip ,tui  ) += funct_(ui)*rowtemp[1][0];
            elemat(tvip ,tuip ) += funct_(ui)*rowtemp[1][1];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        const double fac_gdt_afgdt_tauM_facMtau         =fac*gamma*dt*afgdt*tauM*facMtau;

        const double fac_gdt_afgdt_tauM_facMtau_svelaf_x=fac_gdt_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_y=fac_gdt_afgdt_tauM_facMtau*svelaf_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tuipp =3*ui+2;

          double coltemp[2][2];

          coltemp[0][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(1,ui);

          coltemp[1][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            /*
                 factor: + gamma * dt * alphaF * gamma * dt * tauM *facMtau (rescaled)

              /                                \
             |  ~n+af    /                \     |
             |  u     , | nabla Dp o nabla | v  |
             |           \                /     |
              \                                /

            */

            elemat(tvi  ,tuipp) += coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi);
            elemat(tvip ,tuipp) += coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi);

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
        {
          const double fac_nu_afgdt_afgdt_tauM_facMtau =fac*visceff*afgdt*afgdt*tauM*facMtau;

          double temp[2];

          temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(0);
          temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(1);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int tvi   =3*vi;
            const int tvip  =tvi+1;

            double rowtemp[2][2];

            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);

            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const int tui   =3*ui;
              const int tuip  =tui+1;

              /*
                   factor: - 2.0 * visc * alphaF * gamma * dt * alphaF * gamma * dt * tauM * facMtauM

                    /                                                 \
                   |  ~n+af    / /             /    \  \         \     |
                   |  u     , | | nabla o eps | Dacc |  | o nabla | v  |
                   |           \ \             \    /  /         /     |
                    \                                                 /
              */

              elemat(tvi  ,tui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
              elemat(tvi  ,tuip ) -= derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

              elemat(tvip ,tui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
              elemat(tvip ,tuip ) -= derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
            }
          }
        }// end higher order ele
      } // end if reynolds stab

      //---------------------------------------------------------------
      //
      //               TIME DEPENDENT STABILISATION PART
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == INPAR::FLUID::cross_stress_stab)
      {
        const double fac_afgdt_afgdt_tauM_facMtau  = fac*afgdt   *afgdt*tauM*facMtau;
        const double fac_gdt_afgdt_tauM_facMtau    = fac*gamma*dt*afgdt*tauM*facMtau;
        const double fac_alphaM_afgdt_tauM_facMtau = fac*alphaM  *afgdt*tauM*facMtau;

        double fac_alphaM_afgdt_tauM_velintaf[2];
        fac_alphaM_afgdt_tauM_velintaf[0]=fac_alphaM_afgdt_tauM_facMtau*velintaf_(0);
        fac_alphaM_afgdt_tauM_velintaf[1]=fac_alphaM_afgdt_tauM_facMtau*velintaf_(1);

        double fac_afgdt_afgdt_tauM_facMtau_velintaf[2];
        fac_afgdt_afgdt_tauM_facMtau_velintaf[0]=fac_afgdt_afgdt_tauM_facMtau*velintaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_velintaf[1]=fac_afgdt_afgdt_tauM_facMtau*velintaf_(1);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;


          /*
                  /                              \
                 |    n+af    / ~n+af        \    |
               - |   u     , |  u     o nabla | v |
                 |            \              /    |
                  \                              /
                    +----+
                      ^
                      |
                      +------ linearisation of this part
          */

          /* factor:

                  /                              \
                 |            / ~n+af        \    |
               - |   Dacc  , |  u     o nabla | v |
                 |            \              /    |
                  \                              /
          */
          const double fac_afgdt_conv_subaf_vi=fac_afgdt*conv_subaf_(vi);

          double aux[2];

          /*
                  /                          \
                 |    n+af   ~n+af            |
               - |   u     , u     o nabla v  |
                 |                            |
                  \                          /
                            +----+
                               ^
                               |
                               +------ linearisation of second part
          */

          /* factor:

                  /                                                   \
                 |    n+af    / / /            \   n+af \         \    |
               - |   u     , | | | Dacc o nabla | u      | o nabla | v |
                 |            \ \ \            /        /         /    |
                  \                                                   /
          */
          aux[0]=vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi);
          aux[1]=vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi);

          double temp[2][2];

          /* factor:

                  /                            \
                 |    n+af    /            \    |
                 |   u     , | Dacc o nabla | v |
                 |            \            /    |
                  \                            /
          */
          temp[0][0]=fac_alphaM_afgdt_tauM_velintaf[0]*derxy_(0,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*aux[0]-fac_afgdt_conv_subaf_vi;
          temp[0][1]=fac_alphaM_afgdt_tauM_velintaf[0]*derxy_(1,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*aux[1];

          temp[1][0]=fac_alphaM_afgdt_tauM_velintaf[1]*derxy_(0,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*aux[0];
          temp[1][1]=fac_alphaM_afgdt_tauM_velintaf[1]*derxy_(1,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*aux[1]-fac_afgdt_conv_subaf_vi;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

	    elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
	    elemat(tvi  ,tuip ) += temp[0][1]*funct_(ui);

	    elemat(tvip ,tui  ) += temp[1][0]*funct_(ui);
	    elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);

	  } // ui
	} // vi

        double fac_gdt_afgdt_tauM_facMtau_velintaf[2];
        fac_gdt_afgdt_tauM_facMtau_velintaf[0]=fac_gdt_afgdt_tauM_facMtau*velintaf_(0);
        fac_gdt_afgdt_tauM_facMtau_velintaf[1]=fac_gdt_afgdt_tauM_facMtau*velintaf_(1);

	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int tuipp =3*ui+2;

	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int tvi   =3*vi;
	    const int tvip  =tvi+1;

            /* factor: tauM, rescaled by gamma*dt

                         /                                      \
                        |    n+af    / /          \         \    |
                        |   u     , | |  nabla Dp  | o nabla | v |
                        |            \ \          /         /    |
                         \                                      /
            */
	    const double aux=derxy_(0,vi)*derxy_(0,ui)+derxy_(1,vi)*derxy_(1,ui);

	    elemat(tvi  ,tuipp) += fac_gdt_afgdt_tauM_facMtau_velintaf[0]*aux;
	    elemat(tvip ,tuipp) += fac_gdt_afgdt_tauM_facMtau_velintaf[1]*aux;
	  } // vi
	} // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int tvi   =3*vi;
          const int tvip  =tvi+1;

          /* factor: tauM*afgdt

                  /                                                   \
                 |    n+af    / / /  n+af       \       \         \    |
                 |   u     , | | |  u    o nabla | Dacc  | o nabla | v |
                 |            \ \ \             /       /         /    |
                  \                                                   /
          */
          double temp[2][2];

          temp[0][0]=fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*derxy_(0,vi);
          temp[0][1]=fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*derxy_(1,vi);

          temp[1][0]=fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*derxy_(0,vi);
          temp[1][1]=fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*derxy_(1,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui   =3*ui;
            const int tuip  =tui+1;

	    elemat(tvi  ,tui  ) += temp[0][0]*conv_c_af_(ui);
	    elemat(tvi  ,tuip ) += temp[0][1]*conv_c_af_(ui);

	    elemat(tvip ,tui  ) += temp[1][0]*conv_c_af_(ui);
	    elemat(tvip ,tuip ) += temp[1][1]*conv_c_af_(ui);
	  } // ui
	} // vi

	if (higher_order_ele && newton!=INPAR::FLUID::minimal)
	{
	  const double fac_nu_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

	  double temp[2];

	  temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*velintaf_(0);
	  temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*velintaf_(1);


	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int tvi   =3*vi;
	    const int tvip  =tvi+1;

	    double rowtemp[2][2];

	    rowtemp[0][0]=temp[0]*derxy_(0,vi);
	    rowtemp[0][1]=temp[0]*derxy_(1,vi);

	    rowtemp[1][0]=temp[1]*derxy_(0,vi);
	    rowtemp[1][1]=temp[1]*derxy_(1,vi);

	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int tui   =3*ui;
	      const int tuip  =tui+1;

	      /*
		 factor: 2.0 * visc * alphaF * gamma * dt * tauM

                    /                                                \
                   |   n+af   / /             /    \  \         \     |
                 - |  u    , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /
	      */

	      elemat(tvi  ,tui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(2,ui)*rowtemp[0][1];
	      elemat(tvi  ,tuip ) -= derxy2_(2,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1];

	      elemat(tvip ,tui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(2,ui)*rowtemp[1][1];
	      elemat(tvip ,tuip ) -= derxy2_(2,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1];
	    } //  ui
	  } // vi
	} // hoel
      } // cross

      //---------------------------------------------------------------
      //
      //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //
      //                    PRESSURE STABILISATION
      //
      //---------------------------------------------------------------
      if(pspg == INPAR::FLUID::pstab_use_pspg)
      {
        const double fac_afgdt_gamma_dt_tauM_facMtau  = fac*afgdt*gamma*dt*tauM*facMtau;
        const double fac_gdt_gdt_tauM_facMtau         = fac*gamma*dt*tauM*facMtau*gamma*dt;
        const double fac_alphaM_gamma_dt_tauM_facMtau = fac*alphaM*gamma*dt*tauM*facMtau;

        if(higher_order_ele  && newton!=INPAR::FLUID::minimal)
        {
          const double fac_visceff_afgdt_gamma_dt_tauM_facMtau
            =
            fac*visceff*afgdt*gamma*dt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);


            const double pspg_diffusion_inertia_convect_0_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(0,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_1_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(1,ui)-inertia_and_conv_ui;

            const double scaled_derxy2_2_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(2,ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int tvipp =3*vi+2;

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

              /* pressure stabilisation --- diffusion  */


              /*
                           gamma*dt*tau_M
            factor:  ------------------------------ * alpha_F*gamma*dt * nu
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                                  \
                   |                 /    \             |
                   |  2*nabla o eps | Dacc | , nabla q  |
                   |                 \    /             |
                    \                                  /
              */

              elemat(tvipp,tui ) -=
                derxy_(0,vi)*pspg_diffusion_inertia_convect_0_ui
                +
                derxy_(1,vi)*scaled_derxy2_2_ui;
              elemat(tvipp,tuip) -=
                derxy_(0,vi)*scaled_derxy2_2_ui
                +
                derxy_(1,vi)*pspg_diffusion_inertia_convect_1_ui;
            }
          }
        }
        else
        {
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int tui  =3*ui;
            const int tuip =tui+1;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi) // loop rows
            {
              const int tvipp =3*vi+2;

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

              elemat(tvipp,tui  ) +=derxy_(0,vi)*inertia_and_conv_ui;
              elemat(tvipp,tuip ) +=derxy_(1,vi)*inertia_and_conv_ui;
            }
          }
        } // neglect viscous linearisations, do just inertia and convective

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int tuipp=3*ui+2;

          const double scaled_derxy_0=fac_gdt_gdt_tauM_facMtau*derxy_(0,ui);
          const double scaled_derxy_1=fac_gdt_gdt_tauM_facMtau*derxy_(1,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            /* pressure stabilisation --- pressure   */

            /*
                          gamma*dt*tau_M
            factor:  ------------------------------, rescaled by gamma*dt
                     alpha_M*tau_M+alpha_F*gamma*dt


                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
            */

            elemat(vi*3+2,tuipp) +=
              (scaled_derxy_0*derxy_(0,vi)
               +
               scaled_derxy_1*derxy_(1,vi));

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (newton==INPAR::FLUID::Newton) // if pspg and newton
        {

          for (int vi=0; vi<iel; ++vi) // loop columns
          {
            const int tvipp=3*vi+2;

            const double a=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0));
            const double b=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1));

            for (int ui=0; ui<iel; ++ui)  // loop rows
            {
              const int tui=3*ui;
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

              elemat(tvipp,tui  ) += a*funct_(ui);
              elemat(tvipp,tui+1) += b*funct_(ui);
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
        }// end if pspg and newton
      } // end pressure stabilisation

      //---------------------------------------------------------------
      //
      //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
      //            VISCOUS STABILISATION TERMS FOR (A)GLS
      //
      //---------------------------------------------------------------
      if (higher_order_ele)
      {
        if(vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_gls)
        {
          const double tauMqs = afgdt*tauM*facMtau;

          const double fac_visc_tauMqs_alphaM        = vstabfac*fac*visc*tauMqs*alphaM;
          const double fac_visc_tauMqs_afgdt         = vstabfac*fac*visc*tauMqs*afgdt;
          const double fac_visc_tauMqs_afgdt_visceff = vstabfac*fac*visc*tauMqs*afgdt*visceff;
          const double fac_visc_tauMqs_gamma_dt      = vstabfac*fac*visc*tauMqs*gamma*dt;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;

            const double inertia_and_conv
              =
              fac_visc_tauMqs_alphaM*funct_(ui)+fac_visc_tauMqs_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;
              /* viscous stabilisation --- inertia     */

              /* factor:

                                        alphaF*gamma*tauM*dt
                     +(-)alphaM*nu* ---------------------------
                                    alphaM*tauM+alphaF*gamma*dt

                     /                      \
                    |                        |
                    |  Dacc , 2*div eps (v)  |
                    |                        |
                     \                      /
              */

              /* viscous stabilisation --- convection */
              /*  factor:
                                         alphaF*gamma*dt*tauM
              +(-)alphaF*gamma*dt*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                       /                                    \
                      |  / n+af       \                      |
                      | | c    o nabla | Dacc, 2*div eps (v) |
                      |  \            /                      |
                       \                                    /

              */


              const double a = inertia_and_conv*derxy2_(2,vi);

              elemat(tvi  ,tui  ) += inertia_and_conv*viscs2_(0,vi);
              elemat(tvi  ,tuip ) += a;
              elemat(tvip ,tui  ) += a;
              elemat(tvip ,tuip ) += inertia_and_conv*viscs2_(1,vi);
            }
          }

          for (int ui=0;ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;

            for (int vi=0;vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /* viscous stabilisation --- diffusion  */

              /* factor:

                                             alphaF*gamma*tauM*dt
                -(+)alphaF*gamma*dt*nu*nu ---------------------------
                                          alphaM*tauM+alphaF*gamma*dt

                    /                                        \
                   |                  /    \                  |
                   |  2* nabla o eps | Dacc | , 2*div eps (v) |
                   |                  \    /                  |
                    \                                        /
              */

              const double a = fac_visc_tauMqs_afgdt_visceff*
                               (viscs2_(0,vi)*derxy2_(2,ui)
                                +
                                derxy2_(2,vi)*viscs2_(1,ui));

              elemat(tvi  ,tuip) -= a;
              elemat(tuip ,tvi ) -= a;

              elemat(tvi  ,tui ) -= fac_visc_tauMqs_afgdt_visceff*
		                    (viscs2_(0,ui)*viscs2_(0,vi)
                                     +
                                     derxy2_(2,ui)*derxy2_(2,vi));

              elemat(tvip,tuip ) -= fac_visc_tauMqs_afgdt_visceff*
                                    (derxy2_(2,ui)*derxy2_(2,vi)
                                     +
                                     viscs2_(1,ui)*viscs2_(1,vi));


            }
          }

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int tui    =3*ui;
            const int tuip   =tui+1;
            const int tuipp  =tuip+1;

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int tvi   =3*vi;
              const int tvip  =tvi+1;

              /* viscous stabilisation --- pressure   */

              /* factor:

                                    alphaF*gamma*tauM*dt
                       +(-)nu * ---------------------------, rescaled by gamma*dt
                                alphaM*tauM+alphaF*gamma*dt


                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
              */
              elemat(tvi  ,tuipp) += fac_visc_tauMqs_gamma_dt*
                                     (derxy_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy_(1,ui)*derxy2_(2,vi));
              elemat(tvip ,tuipp) += fac_visc_tauMqs_gamma_dt*
                                     (derxy_(0,ui)*derxy2_(2,vi)
                                      +
                                      derxy_(1,ui)*viscs2_(1,vi));
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton==INPAR::FLUID::Newton)
          {

            double temp[2][2];
            for (int vi=0; vi<iel; ++vi) // loop columns (solution for matrix, test function for vector)
            {
              const int tvi    =3*vi;
              const int tvip   =tvi+1;

              temp[0][0]=(viscs2_(0,vi)*vderxyaf_(0,0)
			  +
                          derxy2_(2,vi)*vderxyaf_(1,0))*fac_visc_tauMqs_afgdt;
              temp[1][0]=(viscs2_(0,vi)*vderxyaf_(0,1)
			  +
                          derxy2_(2,vi)*vderxyaf_(1,1))*fac_visc_tauMqs_afgdt;
              temp[0][1]=(derxy2_(2,vi)*vderxyaf_(0,0)
			  +
                          viscs2_(1,vi)*vderxyaf_(1,0))*fac_visc_tauMqs_afgdt;
              temp[1][1]=(derxy2_(2,vi)*vderxyaf_(0,1)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,1))*fac_visc_tauMqs_afgdt;

              for (int ui=0; ui<iel; ++ui)  // loop rows (test functions for matrix)
              {
                const int tui    =3*ui;
                const int tuip   =tui+1;

                /* viscous stabilisation --- convection
                     factor:
                                         alphaF*gamma*dt*tauM
              +(-)alphaF*gamma*dt*nu* ---------------------------
                                      alphaM*tauM+alphaF*gamma*dt

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /


                */
                elemat(tvi  ,tui  ) += temp[0][0]*funct_(ui);
                elemat(tvi  ,tuip ) += temp[1][0]*funct_(ui);
                elemat(tvip ,tui  ) += temp[0][1]*funct_(ui);
                elemat(tvip ,tuip ) += temp[1][1]*funct_(ui);

              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)

          } // end if (a)gls and newton
        } // end (a)gls stabilisation
      } // end higher_order_element

    } // compute_elemat

    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //
    //         RIGHT HAND SIDE, TIME-DEPENDENT SUBGRID SCALES
    //
    //---------------------------------------------------------------
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    //
    // (MODIFIED) GALERKIN PART, SUBSCALE ACCELERATION STABILISATION
    //
    //---------------------------------------------------------------
    if(inertia == INPAR::FLUID::inertia_stab_keep
       ||
       inertia == INPAR::FLUID::inertia_stab_keep_complete)
    {

      double aux_x =(-svelaf_(0)/tauM-pderxynp_(0)-convaf_old_(0)) ;
      double aux_y =(-svelaf_(1)/tauM-pderxynp_(1)-convaf_old_(1)) ;

      if(higher_order_ele)
      {
        const double fact =visceff;

        aux_x += fact*viscaf_old_(0);
        aux_y += fact*viscaf_old_(1);
      }

      const double fac_sacc_plus_resM_not_partially_integrated_x =fac*aux_x;
      const double fac_sacc_plus_resM_not_partially_integrated_y =fac*aux_y;

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int tui=3*ui;
        //---------------------------------------------------------------
        //
        //     GALERKIN PART I AND SUBSCALE ACCELERATION STABILISATION
        //
        //---------------------------------------------------------------
        /*  factor: +1

               /             \     /                     \
              |   ~ n+am      |   |     n+am    n+af      |
              |  acc     , v  | + |  acc     - f     , v  |
              |     (i)       |   |     (i)               |
               \             /     \	                 /


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

        elevec(tui  ) -= fac_sacc_plus_resM_not_partially_integrated_x*funct_(ui) ;
        elevec(tui+1) -= fac_sacc_plus_resM_not_partially_integrated_y*funct_(ui) ;
      }
    }
    else
    {
      //---------------------------------------------------------------
      //
      //        GALERKIN PART, NEGLECTING SUBSCALE ACCLERATIONS
      //
      //---------------------------------------------------------------
      const double fac_inertia_dead_load_x
        =
        fac*(accintam_(0)-bodyforceaf_(0));

      const double fac_inertia_dead_load_y
        =
        fac*(accintam_(1)-bodyforceaf_(1));

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int tui=3*ui;
        /* inertia terms */

        /*  factor: +1

               /             \
              |     n+am      |
              |  acc     , v  |
              |               |
               \             /
        */

        /* body force (dead load...) */

        /*  factor: -1

               /           \
              |   n+af      |
              |  f     , v  |
              |             |
               \           /
        */

        elevec(tui  ) -= funct_(ui)*fac_inertia_dead_load_x;
        elevec(tui+1) -= funct_(ui)*fac_inertia_dead_load_y;
      }
    }
    //---------------------------------------------------------------
    //
    //            GALERKIN PART 2, REMAINING EXPRESSIONS
    //
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    //
    //         RESIDUAL BASED CONTINUITY STABILISATION
    //          (the original version proposed by Codina)
    //
    //---------------------------------------------------------------

    const double fac_prenp_      = fac*prenp_-fac*tauC*divunp_;

    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
    {
      const int tui =3*ui;
      /* pressure */

      /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
      */

      /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
      */

      elevec(tui  ) += fac_prenp_*derxy_(0,ui) ;
      elevec(tui+1) += fac_prenp_*derxy_(1,ui) ;
    }

    const double visceff_fac=visceff*fac;

    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
    {
      const int tui=3*ui;

      /* viscous term */

      /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
      */

      elevec(tui  ) -= visceff_fac*
                       (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                        +
                        derxy_(1,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0)));
      elevec(tui+1) -= visceff_fac*
	               (derxy_(0,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                        +
                        derxy_(1,ui)*vderxyaf_(1,1)*2.0);
    }

    const double fac_divunp  = fac*divunp_;

    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
    {
      /* continuity equation */

      /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
      */

      elevec(ui*3+2) -= fac_divunp*funct_(ui);
    } // end loop rows (solution for matrix, test function for vector)

    /*
                /                             \
               |  / n+af       \    n+af       |
              +| | u    o nabla |  u      , v  |
               |  \ G          /               |
                \                             /
    */

    double fac_gridconv[2];
    fac_gridconv[0] = -fac*convu_G_af_old_(0);
    fac_gridconv[1] = -fac*convu_G_af_old_(1);

    //---------------------------------------------------------------
    //
    //         STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //
    //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
    //        REYNOLDS CONTRIBUTION FOR CONVECTION DOMINATED FLOWS
    //         CROSS CONTRIBUTION FOR CONVECTION DOMINATED FLOWS
    //
    //---------------------------------------------------------------
    /*
          factor: -1.0

               /                     \
              |                 / \   |
              |  u X u , nabla | v |  |
              |                 \ /   |
               \                     /
    */
    double conv_and_cross_and_re[4];

    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs || cross == INPAR::FLUID::cross_stress_stab)
    {
      /*
                  /                             \
                 |     n+af    n+af              |
               - |  ( u     x u    ) ,  nabla v  |
                 |                               |
                  \                             /

                  /                             \
                 |     n+af   ~n+af              |
               - |  ( u     x u    ) ,  nabla v  |
                 |                               |
                  \                             /
      */
      conv_and_cross_and_re[0]=-velintaf_(0)*fac*(svelaf_(0)+velintaf_(0));
      conv_and_cross_and_re[1]=-velintaf_(0)*fac*(svelaf_(1)+velintaf_(1));
      conv_and_cross_and_re[2]=-velintaf_(1)*fac*(svelaf_(0)+velintaf_(0));
      conv_and_cross_and_re[3]=-velintaf_(1)*fac*(svelaf_(1)+velintaf_(1));
    }
    else
    {
      /*
                  /                             \
                 |     n+af    n+af              |
               - |  ( u     x u    ) ,  nabla v  |
                 |                               |
                  \                             /
      */
      conv_and_cross_and_re[0]=-velintaf_(0)*velintaf_(0)*fac;
      conv_and_cross_and_re[1]=-velintaf_(0)*velintaf_(1)*fac;
      conv_and_cross_and_re[2]=-velintaf_(1)*velintaf_(0)*fac;
      conv_and_cross_and_re[3]=-velintaf_(1)*velintaf_(1)*fac;
    }

    if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
    {
      /*
                  /                             \
                 |    ~n+af   ~n+af              |
               - |  ( u     x u    ) ,  nabla v  |
                 |                               |
                  \                             /
      */

      conv_and_cross_and_re[0]-=fac*svelaf_(0)*svelaf_(0);
      conv_and_cross_and_re[1]-=fac*svelaf_(0)*svelaf_(1);
      conv_and_cross_and_re[2]-=fac*svelaf_(1)*svelaf_(0);
      conv_and_cross_and_re[3]-=fac*svelaf_(1)*svelaf_(1);
    }

    for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
    {
      int tui=3*ui;
      /* gridconv with funct                                */
      /* conv, cross, reynolds with derxy                   */

      elevec(tui++) -=
	fac_gridconv[0]*funct_(ui)
	+
	derxy_(0,ui)*conv_and_cross_and_re[0]
	+
	derxy_(1,ui)*conv_and_cross_and_re[1];
      elevec(tui  ) -=
	fac_gridconv[1]*funct_(ui)
	+
	derxy_(0,ui)*conv_and_cross_and_re[2]
	+
	derxy_(1,ui)*conv_and_cross_and_re[3];
    }

    if(supg == INPAR::FLUID::convective_stab_supg)
    {
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
	int tui=3*ui;

	const double fac_conv_c_af_ui = fac*conv_c_af_(ui);
	/*
	  SUPG stabilisation


                  /                             \
                 |  ~n+af    / n+af        \     |
               - |  u     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /
	*/

	elevec(tui++) += fac_conv_c_af_ui*svelaf_(0);
	elevec(tui  ) += fac_conv_c_af_ui*svelaf_(1);

      } // end loop rows
    } // end supg

    //---------------------------------------------------------------
    //
    //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //                    PRESSURE STABILISATION
    //
    //---------------------------------------------------------------
//    if(pspg == INPAR::FLUID::pstab_use_pspg)
//    {
//
//      const double fac_svelnpx                      = fac*(ele->Svelnp())(0,iquad);
//      const double fac_svelnpy                      = fac*(ele->Svelnp())(1,iquad);
//
//      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
//      {
        /* factor: -1

                       /                 \
                      |  ~n+1             |
                      |  u    , nabla  q  |
                      |   (i)             |
                       \                 /
        */
//
//        elevec(ui*3 + 2) += fac_svelnpx*derxy_(0,ui)+fac_svelnpy*derxy_(1,ui);
//
//      } // end loop rows (solution for matrix, test function for vector)
//    }

    //---------------------------------------------------------------
    //
    //       STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //             VISCOUS STABILISATION (FOR (A)GLS)
    //
    //---------------------------------------------------------------
    if (higher_order_ele)
    {
      if (vstab != INPAR::FLUID::viscous_stab_none)
      {
        const double fac_visc_svelaf_x = vstabfac*fac*visc*svelaf_(0);
        const double fac_visc_svelaf_y = vstabfac*fac*visc*svelaf_(1);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int tui=3*ui;
          /*
                 /                        \
                |  ~n+af                   |
                |  u      , 2*div eps (v)  |
                |                          |
                 \                        /

          */
          elevec(tui  ) += fac_visc_svelaf_x*viscs2_(0,ui)
	                   +
                           fac_visc_svelaf_y*derxy2_(2,ui);

          elevec(tui+1) += fac_visc_svelaf_x*derxy2_(2,ui)
                           +
                           fac_visc_svelaf_y*viscs2_(1,ui);

        } // end loop rows (solution for matrix, test function for vector)
      } // endif (a)gls
    }// end if higher order ele

  } // end loop iquad

  return;
} // Sysmat_cons_td


/*----------------------------------------------------------------------*
 |  extract velocities, pressure and accelerations from the global      |
 |  distributed vectors                           (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::ExtractValuesFromGlobalVectors_2D(
  const bool                        is_ale        ,
  const DRT::Discretization&        discretization,
  const std::vector<int>&           lm            ,
  LINALG::Matrix<iel,1>&            eprenp        ,
  LINALG::Matrix<nsd_,iel>&            evelnp        ,
  LINALG::Matrix<nsd_,iel>&            evelaf        ,
  LINALG::Matrix<nsd_,iel>&            eaccam        ,
  LINALG::Matrix<nsd_,iel>&            edispnp       ,
  LINALG::Matrix<nsd_,iel>&            egridvelaf
  )
{
  // velocity and pressure values (current iterate, n+1)
  RCP<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

  // velocities    (intermediate time step, n+alpha_F)
  RCP<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

  // accelerations (intermediate time step, n+alpha_M)
  RCP<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");

  if (velnp==Teuchos::null || velaf==Teuchos::null || accam==Teuchos::null)
  {
    dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
  }

  // extract local values from the global vectors
  std::vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  std::vector<double> myvelaf(lm.size());
  DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);

  std::vector<double> myaccam(lm.size());
  DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);

  // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
  // Additionally only the 'velocity' components of my_velaf
  // and my_accam are important!
  for (int i=0;i<iel;++i)
  {
    int ti    =3*i;

    eaccam(0,i) = myaccam[ti  ];
    evelnp(0,i) = myvelnp[ti  ];
    evelaf(0,i) = myvelaf[ti++];
    evelnp(1,i) = myvelnp[ti  ];
    eaccam(1,i) = myaccam[ti  ];
    evelaf(1,i) = myvelaf[ti++];
    eprenp(i)   = myvelnp[ti  ];
  }

  if(is_ale)
  {
    // get most recent displacements
    RCP<const Epetra_Vector> dispnp
      =
      discretization.GetState("dispnp");

    // get intermediate grid velocities
    RCP<const Epetra_Vector> gridvelaf
      =
      discretization.GetState("gridvelaf");

    if (dispnp==Teuchos::null || gridvelaf==Teuchos::null)
    {
      dserror("Cannot get state vectors 'dispnp' and/or 'gridvelaf'");
    }

    std::vector<double> mydispnp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

    std::vector<double> mygridvelaf(lm.size());
    DRT::UTILS::ExtractMyValues(*gridvelaf,mygridvelaf,lm);

    // extract velocity part from "mygridvelaf" and get
    // set element displacements
    for (int i=0;i<iel;++i)
    {
      int ti    =3*i;
      int tip   =ti+1;

      egridvelaf(0,i) = mygridvelaf[ti ];
      egridvelaf(1,i) = mygridvelaf[tip];

      edispnp(0,i)    = mydispnp   [ti ];
      edispnp(1,i)    = mydispnp   [tip];
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a SurfaceNeumann condition     |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::GetNodalBodyForce_2D(
  const Fluid* ele ,
  const double  time)
{
  constant_bodyforce_=false;

  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique SurfaceNeumann condition
  DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);

  if (myneumcond.size()>1)
  {
    dserror("more than one VolumeNeumann cond on one node");
  }

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<std::vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
        //curvefac = DRT::Problem::Instance()->Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const std::vector<int>*    onoff = myneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myneumcond[0]->Get<std::vector<double> >("val"  );

    // set this condition to the edeadng array
    for(int isd=0;isd<2;isd++)
    {
      const double value=(*onoff)[isd]*(*val)[isd]*curvefac;

      for (int jnode=0; jnode<iel; jnode++)
      {
        edeadaf_(isd,jnode) = value;
      }
    }

    // this is a constant bodyforce
    constant_bodyforce_=true;
  }
  else
  {
    // we have no dead load
    edeadaf_ = 0.;

    // this is a constant bodyforce
    constant_bodyforce_=true;
  }

  return;
}



/*----------------------------------------------------------------------*
 |                                                                      |
 | Get all global shape functions, first and eventually second          |
 | derivatives in a gausspoint                                          |
 |                                                         gammi 12/08  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::FluidGenalphaResVMM<distype>::ShapeFunctionsFirstAndSecondDerivatives_2D(
  const Fluid*                                 ele      ,
  const int                                  &  iquad    ,
  const DRT::UTILS::IntPointsAndWeights<nsd_>&        intpoints,
  const std::vector<Epetra_SerialDenseVector>&  myknots  ,
  const bool                                    hoel
  )
{
  // set gauss point coordinates
  LINALG::Matrix<2,1> gp;

  gp(0)=intpoints.IP().qxg[iquad][0];
  gp(1)=intpoints.IP().qxg[iquad][1];


  if(not (distype == DRT::Element::nurbs4
       ||
       distype == DRT::Element::nurbs9))
  {
    // get values of shape functions and derivatives in the gausspoint
    DRT::UTILS::shape_function_2D(funct_,gp(0),gp(1),distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv_,gp(0),gp(1),distype);

    if (hoel)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,gp(0),gp(1),distype);
    }
  }
  else if (ele->Shape()==Fluid::nurbs8 || ele->Shape()==Fluid::nurbs27)
    dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());
  else
  {
    if (hoel)
    {
      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv_deriv2
        (funct_  ,
         deriv_  ,
         deriv2_ ,
         gp      ,
         myknots ,
         weights_,
         distype );
    }
    else
    {
      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct_  ,
         deriv_  ,
         gp      ,
         myknots ,
         weights_,
         distype );
    }
  }

  // get transposed Jacobian matrix and determinant
  //
  //        +-       -+ T      +-       -+
  //        | dx   dx |        | dx   dy |
  //        | --   -- |        | --   -- |
  //        | dr   ds |        | dr   dr |
  //        |         |    =   |         |
  //        | dy   dy |        | dx   dy |
  //        | --   -- |        | --   -- |
  //        | dr   ds |        | ds   ds |
  //        +-       -+        +-       -+
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
  for(int rr=0;rr<2;++rr)
  {
    for(int mm=0;mm<2;++mm)
    {
      xjm_(rr,mm)=0.0;
      for(int i=0;i<iel;++i)
      {
        xjm_(rr,mm)+=deriv_(rr,i)*xyze_(mm,i);
      }
    }
  }

  // The determinant ist computed using Sarrus's rule
  const double det = xjm_(0,0)*xjm_(1,1)-xjm_(0,1)*xjm_(1,0);

  // check for degenerated elements
  if (det < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  }

  // set total integration factor
  double fac = intpoints.IP().qwgt[iquad]*det;

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

  // compute global derivates at integration point
  //
  //   dN    +-----  dN (xi)    dxi
  //     i    \        i           k
  //   --- =   +     ------- * -----
  //   dx     /        dxi      dx
  //     j   +-----       k       j
  //         node k
  //
  // j : direction of derivative x/y
  //
  for(int rr=0;rr<2;++rr)
  {
    for(int i=0;i<iel;++i)
    {
      derxy_(rr,i)=0.0;
      for(int mm=0;mm<2;++mm)
      {
        derxy_(rr,i)+=xji_(rr,mm)*deriv_(mm,i);
      }
    }
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
  if (hoel)
  {
    // calculate elements of jacobian_bar matrix
    bm_(0,0) =                     xjm_(0,0)*xjm_(0,0);
    bm_(0,1) =                     xjm_(0,1)*xjm_(0,1);
    bm_(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

    bm_(1,0) =                     xjm_(1,0)*xjm_(1,0);
    bm_(1,1) =                     xjm_(1,1)*xjm_(1,1);
    bm_(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

    bm_(2,0) =                     xjm_(0,0)*xjm_(1,0);
    bm_(2,1) =                     xjm_(0,1)*xjm_(1,1);
    bm_(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);


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
      |  	   	    	   	       | d^2x   d^2y |
      |  	   	    	    	       | ----   ---- |
      | 	   	   	   	       | dr^2   dr^2 |
      | 	   	   	   	       |             |
      | 	   	   	   	       | d^2x   d^2y |
      |                    yields    xder2  =  | ----   ---- |
      | 	   	   	   	       | ds^2   ds^2 |
      | 	   	   	   	       |             |
      | 	   	   	   	       | d^2x   d^2y |
      | 	   	   	   	       | ----   ---- |
      | 	   	   	   	       | drds   drds |
      | 	   	   	   	       +-           -+
    */
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<2;++mm)
      {
        xder2_(rr,mm)=deriv2_(rr,0)*xyze_(mm,0);
        for(int i=1;i<iel;++i)
        {
          xder2_(rr,mm)+=deriv2_(rr,i)*xyze_(mm,i);
        }
      }
    }

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
    for(int rr=0;rr<3;++rr)
    {
      for(int i=0;i<iel;++i)
      {
        derxy2_(rr,i)=deriv2_(rr,i);
        for(int mm=0;mm<2;++mm)
        {
          derxy2_(rr,i)-=xder2_(rr,mm)*derxy_(mm,i);
        }
      }
    }

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
    solver.GETRF(3,3,bm_.A(),3,&(pivot[0]),&ierr);

    if (ierr!=0)
    {
      dserror("Unable to perform LU factorisation during computation of derxy2");
    }

    // backward substitution. GETRS replaces the input (chainrulerhs, currently
    // stored on derxy2) with the result
    solver.GETRS('N',3,iel,bm_.A(),3,&(pivot[0]),derxy2_.A(),3,&ierr);

    if (ierr!=0)
    {
      dserror("Unable to perform backward substitution after factorisation of jacobian");
    }
  }
  else
  {
    for(int rr=0;rr<2;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        derxy2_(rr,mm) = 0.0;
      }
    }
  }

  return(fac);
}

/*----------------------------------------------------------------------*
 |                                                                      |
 | calculates material viscosity                            u.may 05/08 |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalVisc_2D(
  Teuchos::RCP<const MAT::Material> material,
  double&                           visc,
  const double &                    rateofshear)
{
  if(material->MaterialType() == INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(material.get());

    double nu_0 = actmat->Nu0();      // parameter for zero-shear viscosity
    double nu_inf = actmat->NuInf();  // parameter for infinite-shear viscosity
    double lambda = actmat->Lambda();  // parameter for characteristic time
    double a = actmat->AParam();      // constant parameter
    double b = actmat->BParam();  	// constant parameter

    // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    const double tmp = std::pow(lambda*rateofshear,b);
    visc = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
  }
  else if(material->MaterialType() == INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(material.get());

    // get material parameters
    double m  	  = actmat->MCons();    // consistency constant
    double delta  = actmat->Delta();     // safety factor
    double a      = actmat->AExp();     // exponent

    // compute viscosity according to a modified power law model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    visc = m * pow((delta + rateofshear), (-1)*a);
  }
  else
    dserror("material type not yet implemented");

  return;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 |  calculation of stabilisation parameter          gammi 12/08         |
 |    (element center or Gaussian point)                                |
 |                                                                      |
 *----------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalcTau_2D(
  const enum INPAR::FLUID::TauType_genalpha             whichtau ,
  const enum INPAR::FLUID::SubscalesTD           tds      ,
  const double &                         gamma    ,
  const double &                         dt       ,
  const double &                         hk       ,
  const double &                         mk       ,
  const double &                         visceff  )
{
  const int dim=2;

  // get velocity norms
  const double vel_normaf = velintaf_.Norm2();
  const double vel_normnp = velintnp_.Norm2();

  if(tds == INPAR::FLUID::subscales_time_dependent)
  {
    //-------------------------------------------------------
    //          TAUS FOR TIME DEPENDENT SUBSCALES
    //-------------------------------------------------------

    if(whichtau == INPAR::FLUID::taylor_hughes_zarins_whiting_jansen)
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
      LINALG::Matrix<dim,dim> G;

      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<dim;++mm)
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
      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
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
      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
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
      LINALG::Matrix<dim,1> g;

      for (int rr=0;rr<dim;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<dim;++mm)
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
      double normgsq = 0.0;

      for (int rr=0;rr<dim;++rr)
      {
        normgsq+=g(rr)*g(rr);
      }

      /*
                                1.0
                  tau  = -----------------
                     C            /      \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);

    }
    else if(whichtau == INPAR::FLUID::franca_barrenechea_valentin_wall
            ||
            whichtau == INPAR::FLUID::fbvw_wo_dt                      )
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
      const double xi_convectaf = std::max(re_convectaf,1.0);

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

      tau_(1) = tau_(0);

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

      const double xi_tau_c = std::min(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else if(whichtau == INPAR::FLUID::smoothed_franca_barrenechea_valentin_wall)
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

      const double xi_convectaf = re_convectaf+exp(-1.0*re_convectaf);

      /*
                                           -x
                                   f(x)=x+e
               xi_convect ^       -
                          |      -
                          |     -
                          |   --
                        1 +---/
                          |  /
                          | /
                          |/
                          +--------------> re_convect

      */

      /* the 4.0 instead of the Franca's definition 2.0 results from the viscous
       * term in the Navier-Stokes-equations, which is scaled by 2.0*nu         */

      tau_(0) = DSQR(hk) / (4.0 * visceff / mk + ( 4.0 * visceff/mk) * xi_convectaf);

      tau_(1) = tau_(0);

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

      const double xi_tau_c = std::min(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else if(whichtau == INPAR::FLUID::codina)
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
    else if(whichtau == INPAR::FLUID::fbvw_gradient_based_hk)
    {
      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<dim,1> normed_velgrad;

      for (int rr=0;rr<dim;++rr)
      {
        normed_velgrad(rr)=0.0;

        for (int mm=0;mm<dim;++mm)
        {
          normed_velgrad(rr)+=vderxyaf_(mm,rr)*vderxyaf_(mm,rr);
        }
        normed_velgrad(rr)=sqrt(normed_velgrad(rr));
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<dim;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<dim;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;

      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        double temp = 0.0;

        for (int mm=0;mm<dim;++mm)
        {
          temp+=normed_velgrad(mm)*derxy_(mm,rr);
        }

        val += fabs(temp);
      } /* end of loop over element nodes */

      const double gradle = 2.0/val;


      //---------------------------------------------- compute tau_Mu = tau_Mp
      /* convective : viscous forces (element reynolds number)*/
      const double re_convectaf = (vel_normaf * gradle / visceff ) * (mk/2.0);
      const double xi_convectaf = std::max(re_convectaf,1.0);

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

      tau_(0) = DSQR(gradle) / (4.0 * visceff / mk + ( 4.0 * visceff/mk) * xi_convectaf);

      tau_(1) = tau_(0);

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
      const double re_convectnp = (vel_normnp * gradle / visceff ) * (mk/2.0);

      const double xi_tau_c = std::min(re_convectnp,1.0);

      tau_(2) = vel_normnp * gradle * 0.5 * xi_tau_c;
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter for time-dependent formulation\n");
    }

  } // end Fluid::subscales_time_dependent
  else
  {
    //-------------------------------------------------------
    //        TAUS FOR THE QUASISTATIC FORMULATION
    //-------------------------------------------------------
    if(whichtau == INPAR::FLUID::taylor_hughes_zarins_whiting_jansen)
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
      LINALG::Matrix<dim,dim> G;
      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<dim;++mm)
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
      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
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
      for (int nn=0;nn<dim;++nn)
      {
        for (int rr=0;rr<dim;++rr)
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
      LINALG::Matrix<dim,1> g;

      for (int rr=0;rr<dim;++rr)
      {
        g(rr) = xji_(rr,0);
        for (int mm=1;mm<dim;++mm)
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
      double normgsq = 0.0;

      for (int rr=0;rr<dim;++rr)
      {
        normgsq+=g(rr)*g(rr);
      }

      /*
                                1.0
                  tau  = -----------------
                     C            /     \
                          tau  * | g * g |
                             M    \-   -/
      */
      tau_(2) = 1./(tau_(0)*normgsq);
    }
    else if (whichtau == INPAR::FLUID::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Wall


      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<dim,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<dim;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=velintaf_(rr)/vel_normaf;
        }
      }
      else
      {
        normed_velintaf(0) = 1.;
        for (int rr=1;rr<dim;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=0.0;
        }
      }

      // get streamlength

      double val = 0.0;

      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {

        double temp=0.0;

        for(int mm=0;mm<dim;++mm)
        {
          temp+=normed_velintaf(mm)*derxy_(mm,rr);
        }

        val += fabs(temp);
      } /* end of loop over element nodes */

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

      const double xi1 = std::max(re1,1.0);
      const double xi2 = std::max(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visceff);    /* convective : viscous forces */

      const double xi_viscous = std::max(re_viscous,1.0);
      const double xi_convect = std::max(re_convect,1.0);

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
      const double xi_tau_c = std::min(re_convect,1.0);
      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;

    }
    else if(whichtau == INPAR::FLUID::franca_barrenechea_valentin_codina)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Codina


      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<dim,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<dim;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=velintaf_(rr)/vel_normaf;
        }
      }
      else
      {
        normed_velintaf(0) = 1.;
        for (int rr=1;rr<dim;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=0.0;
        }
      }

      // get streamlength
      double val = 0.0;

      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        double temp=0.0;

        for(int mm=0;mm<dim;++mm)
        {
          temp+=normed_velintaf(mm)*derxy_(mm,rr);
        }

        val += fabs(temp);
      } /* end of loop over element nodes */

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

      const double xi1 = std::max(re1,1.0);
      const double xi2 = std::max(re2,1.0);

      tau_(0) = timefac * DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)
      const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk)); /* viscous : reactive forces   */
      const double re_convect = mk * vel_normaf * hk / (2.0 * visceff);    /* convective : viscous forces */

      const double xi_viscous = std::max(re_viscous,1.0);
      const double xi_convect = std::max(re_convect,1.0);

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
    else if (whichtau == INPAR::FLUID::codina)
    {

      // time factor
      const double timefac = gamma*dt;

      // Parameter from Codina, Badia (Constants are chosen according to
      // the values in the standard definition above)

      const double CI  = 4.0/mk;
      const double CII = 2.0/mk;

      // in contrast to the original definition, we neglect the influence of
      // the subscale velocity on velnormaf
      tau_(0)=1.0/(1./timefac+CI*visceff/(hk*hk)+CII*vel_normaf/hk);

      tau_(1)=tau_(0);

      tau_(2)=(hk*hk)/(CI*tau_(0));
    }
    else if(whichtau == INPAR::FLUID::fbvw_wo_dt)
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
      // tau_C: kept Wall definition
      //
      // for the modifications see Codina, Principe, Guasch, Badia
      //    "Time dependent subscales in the stabilized finite  element
      //     approximation of incompressible flow problems"

      //---------------------------------------------- compute tau_Mu = tau_Mp
      /* convective : viscous forces (element reynolds number)*/
      const double re_convectaf = (vel_normaf * hk / visceff ) * (mk/2.0);
      const double xi_convectaf = std::max(re_convectaf,1.0);

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
      tau_(1) = tau_(0);

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

      const double xi_tau_c = std::min(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else if(whichtau == INPAR::FLUID::fbvw_gradient_based_hk)
    {
      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<dim,1> normed_velgrad;

      for (int rr=0;rr<dim;++rr)
      {
        double temp=0.0;

        for(int mm=0;mm<dim;++mm)
        {
          temp+=vderxyaf_(mm,rr)*vderxyaf_(mm,rr);
        }

        normed_velgrad(rr)=sqrt(temp);
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<dim;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<dim;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        double temp=0.0;
        for (int mm=0;mm<dim;++mm) /* loop element nodes */
        {
          temp+=normed_velgrad(mm)*derxy_(mm,rr);
        }

        val += fabs(temp);
      } /* end of loop over element nodes */

      const double gradle = 2.0/val;

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

      // time factor
      const double timefac = gamma*dt;

      const double re1 = 4.0 * timefac * visceff / (mk * DSQR(gradle));   /* viscous : reactive forces   */
      const double re2 = mk * vel_normaf * gradle / (2.0 * visceff);      /* convective : viscous forces */

      const double xi1 = std::max(re1,1.0);
      const double xi2 = std::max(re2,1.0);

      tau_(0) = timefac * DSQR(gradle) / (DSQR(gradle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);
      tau_(1) = tau_(0);

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
      const double xi_tau_c = std::min(re2,1.0);
      tau_(2) = vel_normnp * gradle * 0.5 * xi_tau_c;

    }
    else
    {
      dserror("Unknown definition of stabilisation parameter for quasistatic formulation\n");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  calculates all quantities which are defined at the element center   |
 |  or for the whole element                                            |
 |     o element geometry (xyze_ etc)                                   |
 |     o element volume area_                                           |
 |     o element size hk, constant mk from inverse estimate             |
 |     o dead load                                                      |
 |     o viscosity, effective viscosity                                 |
 |                                                         gammi 01/09  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::SetElementData_2D
(
  Fluid*                                      ele            ,
  const LINALG::Matrix<nsd_,iel>                & edispnp        ,
  const LINALG::Matrix<nsd_,iel>                & evelaf         ,
  const std::vector<Epetra_SerialDenseVector>& myknots        ,
  const double                               & timealphaF     ,
  double                                     & hk             ,
  double                                     & mk             ,
  Teuchos::RCP<const MAT::Material>            material       ,
  double                                     & visc           ,
  double                                     & visceff        )
{
  const int dim=2;

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();

    for(int rr=0;rr<dim;++rr)
    {
      xyze_(rr,inode) = x[rr];
    }
  }

  // get node weights for nurbs elements
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights_(inode) = cp->W();
    }
  }
  else if (ele->Shape()==Fluid::nurbs8 || ele->Shape()==Fluid::nurbs27)
    dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());


  // add displacement, when fluid nodes move in the ALE case
  if (ele->IsAle())
  {
    for (int inode=0; inode<iel; inode++)
    {
      for(int rr=0;rr<dim;++rr)
      {
        xyze_(rr,inode) += edispnp(rr,inode);
      }
    }
  }

  //----------------------------------------------------------------------------
  //                  GET DEAD LOAD IN ELEMENT NODES
  //----------------------------------------------------------------------------
  GetNodalBodyForce_2D(ele,timealphaF);

  //------------------------------------------------------------------
  //                      SET MATERIAL DATA
  //------------------------------------------------------------------
  // check here, if we really have a fluid !!
  if( material->MaterialType() != INPAR::MAT::m_carreauyasuda
      &&
      material->MaterialType() != INPAR::MAT::m_modpowerlaw
      &&
      material->MaterialType() != INPAR::MAT::m_fluid)
  dserror("Material law is not a fluid");

  // get material viscosity
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    // kinematic viscosity
    visc = actmat->Viscosity()/actmat->Density();
  }
  // initialise visceff to visc
  visceff=visc;

  // ---------------------------------------------------------------------------
  // Initialisation of tau computation: mk and hk

  // get element type constant mk for tau and the fssgv_artificial approach
  switch (distype)
  {
      case DRT::Element::tri3:
      case DRT::Element::quad4:
      case DRT::Element::nurbs4:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::tri6:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      case DRT::Element::nurbs9:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // use one point gauss rule to calculate volume at element center
  DRT::UTILS::GaussRule2D integrationrule_stabili=DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
      case DRT::Element::quad4:
      case DRT::Element::nurbs4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      case DRT::Element::nurbs9:
        integrationrule_stabili = DRT::UTILS::intrule_quad_1point;
        break;
      case DRT::Element::tri3:
      case DRT::Element::tri6:
        integrationrule_stabili = DRT::UTILS::intrule_tri_1point;
        break;
      default:
        dserror("invalid discretization type for fluid 2D");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints_onepoint(integrationrule_stabili);

  // shape functions and derivs at element center
  const double wquad = intpoints_onepoint.qwgt[0];

  LINALG::Matrix<dim,1> gp;
  gp(0)=intpoints_onepoint.qxg[0][0];
  gp(1)=intpoints_onepoint.qxg[0][1];

  if(distype == DRT::Element::nurbs4
     ||
     distype == DRT::Element::nurbs9)
  {
    DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
      (funct_  ,
       deriv_  ,
       gp      ,
       myknots ,
       weights_,
       distype );
  }
  else
  {
    DRT::UTILS::shape_function_2D       (funct_,gp(0),gp(1),distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv_,gp(0),gp(1),distype);
  }

  // get transposed Jacobian matrix and determinant
  //
  //        +-       -+ T      +-       -+
  //        | dx   dx |        | dx   dy |
  //        | --   -- |        | --   -- |
  //        | dr   ds |        | dr   dr |
  //        |         |    =   |         |
  //        | dy   dy |        | dx   dy |
  //        | --   -- |        | --   -- |
  //        | dr   ds |        | ds   ds |
  //        +-       -+        +-       -+
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
  for(int rr=0;rr<dim;++rr)
  {
    for(int mm=0;mm<dim;++mm)
    {
      xjm_(rr,mm)=0.0;
      for(int i=0;i<iel;++i)
      {
        xjm_(rr,mm)+=deriv_(rr,i)*xyze_(mm,i);
      }
    }
  }

  // The determinant ist computed using Sarrus's rule
  const double det = xjm_(0,0)*xjm_(1,1)-xjm_(0,1)*xjm_(1,0);

  // check for degenerated elements
  if (det < 0.0)
  {
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  }

  // area is used for some stabilisation parameters
  const double area = wquad*det;

  // get element length hk for tau_M and tau_C: volume-equival. sqrt(area)
  // the same hk is used for the fssgv_artificial approach
  hk = sqrt(area);

  /*------------------------------------------------------------------*/
  /*                                                                  */
  /*                 GET EFFECTIVE VISCOSITY IN ELEMENT               */
  /*                                                                  */
  /* This part is used to specify an effective viscosity.             */
  /*                                                                  */
  /* A cause for the necessity of an effective viscosity might be the */
  /* use of a shear thinning Non-Newtonian fluid                      */
  /*                                                                  */
  /*                            /         \                           */
  /*            visc    = visc | shearrate |                          */
  /*                eff         \         /                           */
  /*                                                                  */
  /*                                                                  */
  /* Mind that at the moment all stabilization (tau and viscous test  */
  /* functions if applied) are based on the material viscosity not    */
  /* the effective viscosity.                                         */
  /* This has to be done before anything else is calculated because   */
  /* we use the same arrays internally. We need hk, mk as well as the */
  /* element data computed above!                                     */
  /*------------------------------------------------------------------*/

  // -------------------------------------------------------------------
  // strain rate based models

  if(material->MaterialType()!= INPAR::MAT::m_fluid)
  {
    //
    //             compute global first derivates
    //
    //
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
    for(int rr=0;rr<dim;++rr)
    {
      for(int i=0;i<iel;++i)
      {
        derxy_(rr,i)=0.0;
        for(int mm=0;mm<dim;++mm)
        {
          derxy_(rr,i)+=xji_(rr,mm)*deriv_(mm,i);
        }
      }
    }

    // ---------------------------------------------------------------
    // compute nonlinear viscosity according to Carreau-Yasuda
    // ---------------------------------------------------------------
    if(material->MaterialType() != INPAR::MAT::m_fluid)
    {
      // compute the rate of strain
      double rateofstrain=GetStrainRate(evelaf  ,derxy_,vderxyaf_  );

      CalVisc_2D(material,visceff,rateofstrain);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 | Interpolates standard quantities to gausspoint, including            |
 |        o accintam_                                                   |
 |        o velintaf_                                                   |
 |        o velintnp_                                                   |
 |        o bodyforce_                                                  |
 |        o prenp_                                                      |
 |        o pderxynp_                                                   |
 |        o vderxyaf_                                                   |
 |        o vderxynp_                                                   |
 |        o divunp_                                                     |
 |        o u_G_af                                                      |
 |        o aleconvintaf_                                               |
 |        o conv_af_old_                                                |
 |        o resM_                                                       |
 |        o conv_c_af_                                                  |
 |        o viscs2_ (if hoel)                                           |
 |        o viscaf_old_ (if hoel)                                       |
 |                                                         gammi 01/09  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::InterpolateToGausspoint_2D(
  Fluid*                                     ele             ,
  const LINALG::Matrix<nsd_,iel>               & egridvaf        ,
  const LINALG::Matrix<nsd_,iel>               & evelnp          ,
  const LINALG::Matrix<iel,1>               & eprenp          ,
  const LINALG::Matrix<nsd_,iel>               & eaccam          ,
  const LINALG::Matrix<nsd_,iel>               & evelaf          ,
  const double                              & visceff         ,
  const bool                                  higher_order_ele)
{
  const double dim=2;

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
  for(int rr=0;rr<dim;++rr)
  {
    accintam_(rr)=funct_(0)*eaccam(rr,0);
    for(int nn=1;nn<iel;++nn)
    {
      accintam_(rr)+=funct_(nn)*eaccam(rr,nn);
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
  for(int rr=0;rr<dim;++rr)
  {
    velintaf_(rr)=funct_(0)*evelaf(rr,0);
    for(int nn=1;nn<iel;++nn)
    {
      velintaf_(rr)+=funct_(nn)*evelaf(rr,nn);
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
  // required for computation of tauC
  for(int rr=0;rr<dim;++rr)
  {
    velintnp_(rr)=funct_(0)*evelnp(rr,0);
    for(int nn=1;nn<iel;++nn)
    {
      velintnp_(rr)+=funct_(nn)*evelnp(rr,nn);
    }
  }

  if(!constant_bodyforce_)
  {
    // get bodyforce in gausspoint, time (n+alpha_F)
    //
    //                 +-----
    //       n+af       \                n+af
    //      f    (x) =   +      N (x) * f
    //                  /        j       j
    //                 +-----
    //                 node j
    //
    for(int rr=0;rr<dim;++rr)
    {
      bodyforceaf_(rr)=funct_(0)*edeadaf_(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        bodyforceaf_(rr)+=funct_(nn)*edeadaf_(rr,nn);
      }
    }
  }
  else
  {
    // a constant bodyforce doesn't require
    // interpolation to gausspoint
    //
    //
    //       n+af       n+af
    //      f    (x) = f     = onst.
    //
    for(int rr=0;rr<dim;++rr)
    {
      bodyforceaf_(rr)=edeadaf_(rr,0);
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
  prenp_=0;
  for(int mm=0;mm<iel;++mm)
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
  for(int rr=0;rr<dim;++rr)
  {
    pderxynp_(rr)=derxy_(rr,0)*eprenp(0);
    for(int nn=1;nn<iel;++nn)
    {
      pderxynp_(rr)+=derxy_(rr,nn)*eprenp(nn);
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
  for(int rr=0;rr<dim;++rr)
  {
    for(int mm=0;mm<dim;++mm)
    {
      vderxyaf_(rr,mm)=derxy_(mm,0)*evelaf(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        vderxyaf_(rr,mm)+=derxy_(mm,nn)*evelaf(rr,nn);
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
  for(int rr=0;rr<dim;++rr)
  {
    for(int mm=0;mm<dim;++mm)
    {
      vderxynp_(rr,mm)=derxy_(mm,0)*evelnp(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        vderxynp_(rr,mm)+=derxy_(mm,nn)*evelnp(rr,nn);
      }
    }
  }

  /* divergence new time step n+1 */
  //
  //                   +-----     n+1
  //          n+1       \     dvel   (x)
  //   div vel   (x) =   +    ----------
  //                    /         dx
  //                   +-----       j
  //                    dim j
  //

  divunp_ = vderxynp_(0,0);
  for(int rr=1;rr<dim;++rr)
  {
    divunp_+=vderxynp_(rr,rr);
  }

  // get ale convective velocity (n+alpha_F,i) at integration point
  for(int rr=0;rr<dim;++rr)
  {
    aleconvintaf_(rr)=velintaf_(rr);
  }

  if (ele->IsAle())
  {
    // u_G is the grid velocity at the integration point,
    // time (n+alpha_F,i)
    //
    //                 +-----
    //       n+af       \                  n+af
    //    u_G    (x) =   +      N (x) * u_G
    //                  /        j         j
    //                 +-----
    //                 node j
    //

    for(int rr=0;rr<dim;++rr)
    {
      u_G_af_(rr)=funct_(0)*egridvaf(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        u_G_af_(rr)+=funct_(nn)*egridvaf(rr,nn);
      }
    }
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

    for(int rr=0;rr<dim;++rr)
    {
      aleconvintaf_(rr)-=u_G_af_(rr);
    }
  }
  else
  {
    for(int rr=0;rr<dim;++rr)
    {
      u_G_af_(rr)=0.0;
    }
  }

  /* Convective term  u_old * grad u_old: */
  /*
  //     /  n+af        \   n+af
  //    |  c     o nabla | u
  //     \              /
  */
  for(int rr=0;rr<dim;++rr)
  {
    convaf_old_(rr)=aleconvintaf_(0)*vderxyaf_(rr,0);
    for(int mm=1;mm<dim;++mm)
    {
      convaf_old_(rr)+=aleconvintaf_(mm)*vderxyaf_(rr,mm);
    }
  }

  // compute residual in gausspoint --- second derivatives only
  // exist for higher order elements, so we subtract them later
  // convaf_old_ is based on ale-convective velocity
  //
  //   n+af         n+am       /   n+af           \     n+af
  //  r    (x) = acc    (x) + | vel    (x) o nabla | vel    (x) +
  //   M                       \                  /
  //                      n+1    n+af
  //             + nabla p    - f             (not higher order)
  //
  for (int rr=0;rr<dim;++rr)
  {
    resM_(rr) = accintam_(rr) + convaf_old_(rr) + pderxynp_(rr) - bodyforceaf_(rr);
  }

  // get convective linearisation (n+alpha_F,i) at integration point
  //
  //                 +-----
  //       n+af       \      n+af      dN
  // conv_c    (x) =   +    c    (x) * --- (x)
  //                  /      j         dx
  //                 +-----              j
  //                  dim j
  //
  for(int nn=0;nn<iel;++nn)
  {
    conv_c_af_(nn)=aleconvintaf_(0)*derxy_(0,nn);

    for(int rr=1;rr<dim;++rr)
    {
      conv_c_af_(nn)+=aleconvintaf_(rr)*derxy_(rr,nn);
    }
  }

  if (higher_order_ele)
  {
    /*--- viscous term  2* grad * epsilon(u): --------------------------*/
    /*   /                                                \
        |   2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz   |
        |                                                  |
        |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz    |
        |                                                  |
        |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz    |
         \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                           */


    /* Viscous term  div epsilon(u_old)
    //
    //              /             \
    //             |     / n+af \  |
    //     nabla o | eps| u      | | =
    //             |     \      /  |
    //              \             / j
    //
    //              / +-----  / +----- dN (x)             +----- dN (x)           \ \
    //             |   \     |   \       k         n+af    \       k         n+af  | |
    //       1.0   |    +    |    +    ------ * vel     +   +   ------- * vel      | |
    //     = --- * |   /     |   /     dx dx       k,i     /     dx dx       k,j   | |
    //       2.0   |  +----- |  +-----   i  j             +-----   i  i            | |
    //              \ node k  \  dim i                     dim i                  / /
    */
    double sum = derxy2_(0,0)+derxy2_(1,0);

    viscs2_(0,0) = sum + derxy2_(0,0);
    viscs2_(1,0) = sum + derxy2_(1,0);

    viscaf_old_(0) = (viscs2_(0,0)*evelaf(0,0)
                      +
                      derxy2_(2,0)*evelaf(1,0));
    viscaf_old_(1) = (derxy2_(2,0)*evelaf(0,0)
                      +
                      viscs2_(1,0)*evelaf(1,0));

    for (int mm=1;mm<iel;++mm)
    {
      sum = derxy2_(0,mm)+derxy2_(1,mm);

      viscs2_(0,mm) = sum + derxy2_(0,mm);
      viscs2_(1,mm) = sum + derxy2_(1,mm);

      viscaf_old_(0) += (viscs2_(0,mm)*evelaf(0,mm)
                         +
                         derxy2_(2,mm)*evelaf(1,mm));
      viscaf_old_(1) += (derxy2_(2,mm)*evelaf(0,mm)
                         +
                         viscs2_(1,mm)*evelaf(1,mm));
    }

    /* the residual is based on the effective viscosity!
    //
    //   n+af         n+am       /   n+af           \     n+af
    //  r    (x) = acc    (x) + | vel    (x) o nabla | vel    (x) +
    //   M                       \                  /
    //                      n+1                     /   n+af \    n+af
    //             + nabla p    - 2*nu nabla o eps | vel      |- f
    //                                              \        /
    */
    for (int rr=0;rr<dim;++rr)
    {
      resM_(rr) -= visceff*viscaf_old_(rr);
    }
  } // end if higher order

  return;
}


