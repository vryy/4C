/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.cpp

\brief Internal implementation of Fluid3 element with a generalised alpha
       time integration.

       This element is designed for the solution of the Navier-Stokes
       equations using a residual based stabilised method. The
       stabilisation terms are derived in a variational multiscale sense.

       Subscales are either treated as quasi-static or time dependent.

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
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_LAPACK.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3GenalphaResVMMInterface* DRT::ELEMENTS::Fluid3GenalphaResVMMInterface::Impl(DRT::ELEMENTS::Fluid3* f3)
{
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Fluid3GenalphaResVMM<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new Fluid3GenalphaResVMM<DRT::Element::hex8>();
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static Fluid3GenalphaResVMM<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new Fluid3GenalphaResVMM<DRT::Element::hex20>();
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static Fluid3GenalphaResVMM<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new Fluid3GenalphaResVMM<DRT::Element::hex27>();
    return fh27;
  }
  case DRT::Element::nurbs8:
  {
    static Fluid3GenalphaResVMM<DRT::Element::nurbs8>* fn8;
    if (fn8==NULL)
      fn8 = new Fluid3GenalphaResVMM<DRT::Element::nurbs8>();
    return fn8;
  }
  case DRT::Element::nurbs27:
  {
    static Fluid3GenalphaResVMM<DRT::Element::nurbs27>* fn27;
    if (fn27==NULL)
      fn27 = new Fluid3GenalphaResVMM<DRT::Element::nurbs27>();
    return fn27;
  }
  case DRT::Element::tet4:
  {
    static Fluid3GenalphaResVMM<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new Fluid3GenalphaResVMM<DRT::Element::tet4>();
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static Fluid3GenalphaResVMM<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new Fluid3GenalphaResVMM<DRT::Element::tet10>();
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static Fluid3GenalphaResVMM<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new Fluid3GenalphaResVMM<DRT::Element::wedge6>();
    return fw6;
  }
  case DRT::Element::wedge15:
  {
    static Fluid3GenalphaResVMM<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new Fluid3GenalphaResVMM<DRT::Element::wedge15>();
    return fw15;
  }
  case DRT::Element::pyramid5:
  {
    static Fluid3GenalphaResVMM<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new Fluid3GenalphaResVMM<DRT::Element::pyramid5>();
    return fp5;
  }
  default:
    dserror("shape %d (%d nodes) not supported", f3->Shape(), f3->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
  |  constructor allocating arrays whose sizes may depend on the number |
  | of nodes of the element                                             |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Fluid3GenalphaResVMM()
// fine-scale subgrid viscosity
  : vart_(0.0)
{
  return;
}


template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Evaluate(
  Fluid3*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RefCountPtr<MAT::Material> mat,
  MATERIAL*                  actmat)
{

  // --------------------------------------------------
  // construct views
  LINALG::Matrix<4*iel,4*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<4*iel,    1> elevec1(elevec1_epetra.A(),true);

  // --------------------------------------------------
  // create matrix objects for nodal values
  LINALG::Matrix<iel,1> eprenp    ;
  LINALG::Matrix<3,iel> evelnp    ;
  LINALG::Matrix<3,iel> evelaf    ;
  LINALG::Matrix<3,iel> eaccam    ;
  LINALG::Matrix<3,iel> edispnp   ;
  LINALG::Matrix<3,iel> egridvelaf;
  LINALG::Matrix<3,iel> fsevelaf  ;

  // --------------------------------------------------
  // set parameters for time integration
  ParameterList& timelist = params.sublist("time integration parameters");

  const double alphaM = timelist.get<double>("alpha_M");
  const double alphaF = timelist.get<double>("alpha_F");
  const double gamma  = timelist.get<double>("gamma");
  const double dt     = timelist.get<double>("dt");
  const double time   = timelist.get<double>("time");

  // --------------------------------------------------
  // set parameters for nonlinear treatment
  string newtonstr=params.get<string>("Linearisation");

  Fluid3::LinearisationAction newton=Fluid3::no_linearisation;
  if(newtonstr=="Newton")
  {
    newton=Fluid3::Newton;
  }
  else if (newtonstr=="fixed_point_like")
  {
    newton=Fluid3::fixed_point_like;
  }
  else if (newtonstr=="minimal")
  {
    newton=Fluid3::minimal;
  }

  // --------------------------------------------------
  // get flag for fine-scale subgrid viscosity
  Fluid3::StabilisationAction fssgv =
    ele->ConvertStringToStabAction(params.get<string>("fs subgrid viscosity","No"));

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");
  
  // specify which residual based stabilisation terms
  // will be used
  Fluid3::StabilisationAction tds      = ele->ConvertStringToStabAction(stablist.get<string>("TDS"));
  Fluid3::StabilisationAction inertia  = ele->ConvertStringToStabAction(stablist.get<string>("TRANSIENT"));
  Fluid3::StabilisationAction pspg     = ele->ConvertStringToStabAction(stablist.get<string>("PSPG"));
  Fluid3::StabilisationAction supg     = ele->ConvertStringToStabAction(stablist.get<string>("SUPG"));
  Fluid3::StabilisationAction vstab    = ele->ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  Fluid3::StabilisationAction cstab    = ele->ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  Fluid3::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  Fluid3::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // select tau definition
  Fluid3::TauType whichtau = Fluid3::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid3::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid3::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid3::codina;
    }
    else if(taudef == "FBVW_without_dt")
    {
      whichtau = Fluid3::fbvw_wo_dt;
    }
    else if(taudef == "Franca_Barrenechea_Valentin_Codina")
    {
      whichtau = Fluid3::franca_barrenechea_valentin_codina;
    }
    else if(taudef == "Smoothed_FBVW")
    {
      whichtau = Fluid3::smoothed_franca_barrenechea_valentin_wall;
    }
    else
    {
      dserror("unknown tau definition\n");
    }
  }

  // flag for higher order elements
  bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());
  
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if(stablist.get<string>("STABTYPE") == "inconsistent")
  {
    higher_order_ele = false;
  }

  // flag conservative form on/off 
  string conservativestr =params.get<string>("CONVFORM");

  bool conservative =false;
  if(conservativestr=="conservative")
  {
    conservative =true;
  }
  else if(conservativestr=="convective")
  {
    conservative =false;
  }

  // --------------------------------------------------
  // set parameters for turbulence model
  ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

  // the default action is no model
  Fluid3::TurbModelAction turb_mod_action = Fluid3::no_model;

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  double Cs            = 0.0;
  double Cs_delta_sq   = 0.0;
  double l_tau         = 0.0;

  // number of the layer in a turbulent plane channel flow --- used 
  // to compute averaged viscosity etc
  int    nlayer        = 0;

  SetParametersForTurbulenceModel(
    ele            ,
    turbmodelparams,
    fssgv          ,
    turb_mod_action, 
    Cs             ,
    Cs_delta_sq    ,
    l_tau          ,
    nlayer         );

  // --------------------------------------------------
  // specify whether to compute the element matrix or not
  const bool compute_elemat = params.get<bool>("compute element matrix");

  // --------------------------------------------------
  // extract velocities, pressure and accelerations from the
  // global distributed vectors
  ExtractValuesFromGlobalVectors(
        fssgv         ,
        ele->is_ale_  ,
        discretization,
        lm            ,
        eprenp        ,
        evelnp        ,
        evelaf        ,
        eaccam        ,
        edispnp       ,
        egridvelaf    ,
        fsevelaf
    );

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<blitz::Array<double,1> > myknots(3);

  // for isogeometric elements
  if(ele->Shape()==Fluid3::nurbs8 || ele->Shape()==Fluid3::nurbs27)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());
  }

  // on output of Sysmat, visceff will contain the computed effective viscosity
  double visceff       = 0.0;

  // --------------------------------------------------
  // calculate element coefficient matrix
  if(!conservative)
  {
    if(tds!=Fluid3::subscales_time_dependent)
    {
      Sysmat_adv_qs(
        ele,
        myknots,
        elemat1,
        elevec1,
        edispnp,
        egridvelaf,
        evelnp,
        eprenp,
        eaccam,
        evelaf,
        fsevelaf,
        actmat,
        alphaM,
        alphaF,
        gamma,
        dt,
        time,
        newton,
        higher_order_ele,
        fssgv,
        inertia,
        pspg,
        supg,
        vstab,
        cstab,
        cross,
        reynolds,
        whichtau,
        turb_mod_action,
        Cs,
        Cs_delta_sq,
        visceff,
        l_tau,
        compute_elemat
        );
    }
    else
    {
      Sysmat_adv_td(
        ele,
        myknots,
        elemat1,
        elevec1,
        edispnp,
        egridvelaf,
        evelnp,
        eprenp,
        eaccam,
        evelaf,
        fsevelaf,
        actmat,
        alphaM,
        alphaF,
        gamma,
        dt,
        time,
        newton,
        higher_order_ele,
        fssgv,
        inertia,
        pspg,
        supg,
        vstab,
        cstab,
        cross,
        reynolds,
        whichtau,
        turb_mod_action,
        Cs,
        Cs_delta_sq,
        visceff,
        l_tau,
        compute_elemat
        );
    }
  }
  else
  {
    if(tds!=Fluid3::subscales_time_dependent)
    {
      Sysmat_cons_qs(
	ele,
	myknots,
	elemat1,
	elevec1,
	edispnp,
	egridvelaf,
	evelnp,
	eprenp,
	eaccam,
	evelaf,
	fsevelaf,
	actmat,
	alphaM,
	alphaF,
	gamma,
	dt,
	time,
	newton,
	higher_order_ele,
	fssgv,
	pspg,
	supg,
	vstab,
	cstab,
	cross,
	reynolds,
	whichtau,
	turb_mod_action,
	Cs,
	Cs_delta_sq,
	visceff,
	l_tau,
	compute_elemat
	);
    }
    else
    {
      Sysmat_cons_td(
        ele             ,
        myknots         ,
        elemat1         ,
        elevec1         ,
        edispnp         ,
        egridvelaf      ,
        evelnp          ,
        eprenp          ,
        eaccam          ,
        evelaf          ,
        fsevelaf        ,
        actmat          ,
        alphaM          ,
        alphaF          ,
        gamma           ,
        dt              ,
        time            ,
        newton          ,
        higher_order_ele,
        fssgv           ,
        inertia         ,
        pspg            ,
        supg            ,
        vstab           ,
        cstab           ,
        cross           ,
        reynolds        ,
        whichtau        ,
        turb_mod_action ,
        Cs              ,
        Cs_delta_sq     ,
        visceff         ,
        l_tau           ,
        compute_elemat
        );
    }
  }


  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    // this output is only interesting for a plane channel flow
    if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
        ==
        "channel_flow_of_height_2")
    {
      if (physical_turbulence_model == "Dynamic_Smagorinsky"
          ||
          physical_turbulence_model ==  "Smagorinsky_with_van_Driest_damping"
        )
      {
        // Cs was changed in Sysmat (Cs->sqrt(Cs/hk)) to compare it with the standard
        // Smagorinsky Cs

        if(ele->Owner() == discretization.Comm().MyPID())
        {
          (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_sum")))         [nlayer]+=Cs;
          (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_delta_sq_sum")))[nlayer]+=Cs_delta_sq;
          (*(turbmodelparams.get<RCP<vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff;
        }
      }
    }
  }

  {
    // This is a very poor way to transport the density to the
    // outside world. Is there a better one?
    double dens = 0.0;
    if(mat->MaterialType()== m_fluid)
      dens = actmat->m.fluid->density;
    else if(mat->MaterialType()== m_carreauyasuda)
      dens = actmat->m.carreauyasuda->density;
    else if(mat->MaterialType()== m_modpowerlaw)
      dens = actmat->m.modpowerlaw->density;
    else
      dserror("no fluid material found");
    
    params.set("density", dens);
  }

  return 0;
}

/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  advective version based on quasistatic subgrid scales              |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Sysmat_adv_qs(
  Fluid3*                                 ele             ,
  std::vector<blitz::Array<double,1> >&   myknots         ,
  LINALG::Matrix<4*iel,4*iel>&            elemat          ,
  LINALG::Matrix<4*iel,1>&                elevec          ,
  const LINALG::Matrix<3,iel>&            edispnp         ,
  const LINALG::Matrix<3,iel>&            egridvaf        ,
  const LINALG::Matrix<3,iel>&            evelnp          ,
  const LINALG::Matrix<iel,1>&            eprenp          ,
  const LINALG::Matrix<3,iel>&            eaccam          ,
  const LINALG::Matrix<3,iel>&            evelaf          ,
  const LINALG::Matrix<3,iel>&            fsevelaf        ,
  const struct _MATERIAL*                 material        ,
  const double                            alphaM          ,
  const double                            alphaF          ,
  const double                            gamma           ,
  const double                            dt              ,
  const double                            time            ,
  const enum Fluid3::LinearisationAction  newton          ,
  const bool                              higher_order_ele,
  const enum Fluid3::StabilisationAction  fssgv           ,
  const enum Fluid3::StabilisationAction  inertia         ,
  const enum Fluid3::StabilisationAction  pspg            ,
  const enum Fluid3::StabilisationAction  supg            ,
  const enum Fluid3::StabilisationAction  vstab           ,
  const enum Fluid3::StabilisationAction  cstab           ,
  const enum Fluid3::StabilisationAction  cross           ,
  const enum Fluid3::StabilisationAction  reynolds        ,
  const enum Fluid3::TauType              whichtau        ,
  const enum Fluid3::TurbModelAction      turb_mod_action ,
  double&                                 Cs              ,
  double&                                 Cs_delta_sq     ,
  double&                                 visceff         ,
  const double                            l_tau           ,
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
  if (vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_gls_only_rhs)
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

  SetElementData(ele            ,
                 edispnp        ,
                 evelaf         ,
                 fsevelaf       ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 fssgv          ,
                 turb_mod_action,
                 l_tau          ,
                 Cs             ,
                 Cs_delta_sq    ,
                 visceff        );
 
 
  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second 
    // derivatives in a gausspoint and integration weight including 
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
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
    if(cross == Fluid3::cross_stress_stab)
    {
      for(int nn=0;nn<iel;++nn)
      {
        conv_resM_(nn)=resM_(0)*derxy_(0,nn);
          
        for(int rr=1;rr<3;++rr)
        {
          conv_resM_(nn)+=resM_(rr)*derxy_(rr,nn);
        }
      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,Fluid3::subscales_quasistatic,gamma,dt,hk,mk,visceff);


    // stabilisation parameters
    const double tauM   = tau_(0);
    const double tauMp  = tau_(1);

    if(cstab == Fluid3::continuity_stab_none)
    {
      tau_(2)=0.0;
    }

    const double tauC   = tau_(2);

    double supg_active_tauM;
    if(supg == Fluid3::convective_stab_supg)
    {
      supg_active_tauM=tauM;
    }
    else
    {
      supg_active_tauM=0.0;
    }

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;

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

      if(reynolds == Fluid3::reynolds_stress_stab)
      {

        /* half of the reynolds linearisation is done by modifiing 
           the supg testfunction, see above */

        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(0)*derxy_(0,nn);
            
          for(int rr=1;rr<3;++rr)
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
        const double fac_alphaM_tauM_tauM_resM_z=fac_alphaM_tauM_tauM*resM_(2);

        const double fac_afgdt_tauM_tauM=fac*afgdt*tauM*tauM;

        double fac_afgdt_tauM_tauM_resM[3];
        fac_afgdt_tauM_tauM_resM[0]=fac_afgdt_tauM_tauM*resM_(0);
        fac_afgdt_tauM_tauM_resM[1]=fac_afgdt_tauM_tauM*resM_(1);
        fac_afgdt_tauM_tauM_resM[2]=fac_afgdt_tauM_tauM*resM_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
            
          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
            
          double inertia_and_conv[3];

          inertia_and_conv[0]=fac_afgdt_tauM_tauM_resM[0]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_tauM_tauM_resM[1]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_y*funct_(ui);
          inertia_and_conv[2]=fac_afgdt_tauM_tauM_resM[2]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_z*funct_(ui);
              
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

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

            elemat(fvi  ,fui  ) -= inertia_and_conv[0]*derxy_(0,vi);
            elemat(fvi  ,fuip ) -= inertia_and_conv[0]*derxy_(1,vi);
            elemat(fvi  ,fuipp) -= inertia_and_conv[0]*derxy_(2,vi);
                                     
            elemat(fvip ,fui  ) -= inertia_and_conv[1]*derxy_(0,vi);
            elemat(fvip ,fuip ) -= inertia_and_conv[1]*derxy_(1,vi);
            elemat(fvip ,fuipp) -= inertia_and_conv[1]*derxy_(2,vi);
                                     
            elemat(fvipp,fui  ) -= inertia_and_conv[2]*derxy_(0,vi);
            elemat(fvipp,fuip ) -= inertia_and_conv[2]*derxy_(1,vi);
            elemat(fvipp,fuipp) -= inertia_and_conv[2]*derxy_(2,vi);
          } // vi
        } // ui

              
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
            
          double temp[3];
          temp[0]=fac_afgdt_tauM_tauM*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi));
          temp[1]=fac_afgdt_tauM_tauM*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi));
          temp[2]=fac_afgdt_tauM_tauM*(vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi));

          double rowtemp[3][3];
            
          rowtemp[0][0]=resM_(0)*temp[0];
          rowtemp[0][1]=resM_(0)*temp[1];
          rowtemp[0][2]=resM_(0)*temp[2];
                                   
          rowtemp[1][0]=resM_(1)*temp[0];
          rowtemp[1][1]=resM_(1)*temp[1];
          rowtemp[1][2]=resM_(1)*temp[2];
                                   
          rowtemp[2][0]=resM_(2)*temp[0];
          rowtemp[2][1]=resM_(2)*temp[1];
          rowtemp[2][2]=resM_(2)*temp[2];
              
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
            
            /* 
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / /            \   n+af \         \     |
             |  resM , | | | Dacc o nabla | u      | o nabla | v  |
             |          \ \ \            /        /         /     |
              \                                                  /      

            */

            elemat(fvi  ,fui  ) -= funct_(ui)*rowtemp[0][0];
            elemat(fvi  ,fuip ) -= funct_(ui)*rowtemp[0][1];
            elemat(fvi  ,fuipp) -= funct_(ui)*rowtemp[0][2];
                                                             
            elemat(fvip ,fui  ) -= funct_(ui)*rowtemp[1][0];
            elemat(fvip ,fuip ) -= funct_(ui)*rowtemp[1][1];
            elemat(fvip ,fuipp) -= funct_(ui)*rowtemp[1][2];
                                                             
            elemat(fvipp,fui  ) -= funct_(ui)*rowtemp[2][0];
            elemat(fvipp,fuip ) -= funct_(ui)*rowtemp[2][1];
            elemat(fvipp,fuipp) -= funct_(ui)*rowtemp[2][2];
          } // ui
        } // vi


        const double fac_gdt_tauM_tauM       =fac*gamma*dt*tauM*tauM;
        const double fac_gdt_tauM_tauM_resM_x=fac_gdt_tauM_tauM*resM_(0);
        const double fac_gdt_tauM_tauM_resM_y=fac_gdt_tauM_tauM*resM_(1);
        const double fac_gdt_tauM_tauM_resM_z=fac_gdt_tauM_tauM*resM_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fuippp =4*ui+3;
            
          double coltemp[3][3];

          coltemp[0][0]=fac_gdt_tauM_tauM_resM_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_tauM_tauM_resM_x*derxy_(1,ui);
          coltemp[0][2]=fac_gdt_tauM_tauM_resM_x*derxy_(2,ui);
          coltemp[1][0]=fac_gdt_tauM_tauM_resM_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_tauM_tauM_resM_y*derxy_(1,ui);
          coltemp[1][2]=fac_gdt_tauM_tauM_resM_y*derxy_(2,ui);
          coltemp[2][0]=fac_gdt_tauM_tauM_resM_z*derxy_(0,ui);
          coltemp[2][1]=fac_gdt_tauM_tauM_resM_z*derxy_(1,ui);
          coltemp[2][2]=fac_gdt_tauM_tauM_resM_z*derxy_(2,ui);
  
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
            
            /* 
                 factor: - gamma * dt * tauM * tauM (rescaled)
             
              /                               \
             |          /                \     |
             |  resM , | nabla Dp o nabla | v  |
             |          \                /     |
              \                               /      

            */

            elemat(fvi  ,fuippp) -= coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi)+coltemp[0][2]*derxy_(2,vi);
            elemat(fvip ,fuippp) -= coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi)+coltemp[1][2]*derxy_(2,vi);
            elemat(fvipp,fuippp) -= coltemp[2][0]*derxy_(0,vi)+coltemp[2][1]*derxy_(1,vi)+coltemp[2][2]*derxy_(2,vi);

          } // vi
        } // ui


        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_nu_afgdt_tauM_tauM=fac*visceff*afgdt*tauM*tauM;
          
          double temp[3];

          temp[0]=fac_nu_afgdt_tauM_tauM*resM_(0);
          temp[1]=fac_nu_afgdt_tauM_tauM*resM_(1);
          temp[2]=fac_nu_afgdt_tauM_tauM*resM_(2);


          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
              
            double rowtemp[3][3];
              
            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);
            rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);
            rowtemp[1][2]=temp[1]*derxy_(2,vi);
              
            rowtemp[2][0]=temp[2]*derxy_(0,vi);
            rowtemp[2][1]=temp[2]*derxy_(1,vi);
            rowtemp[2][2]=temp[2]*derxy_(2,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int fui   =4*ui;
              const int fuip  =fui+1;
              const int fuipp =fui+2;
                
              /* 
                   factor: + 2.0 * visc * alphaF * gamma * dt * tauM * tauM 

                    /                                                \
                   |          / /             /    \  \         \     |
                   |  resM , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /      
              */

              elemat(fvi  ,fui  ) += viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
              elemat(fvi  ,fuip ) += derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
              elemat(fvi  ,fuipp) += derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
              elemat(fvip ,fui  ) += viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
              elemat(fvip ,fuip ) += derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
              elemat(fvip ,fuipp) += derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
              elemat(fvipp,fui  ) += viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
              elemat(fvipp,fuip ) += derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
              elemat(fvipp,fuipp) += derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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
        const int fui  =4*ui;
        const int fuip =fui+1;
        const int fuipp=fui+2;
                
        /* GALERKIN inertia term (intermediate) + convection (intermediate) */
        const double inertia_and_conv_ui
          = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);
	      
        /* viscous term (intermediate), diagonal parts */
        const double fac_visceff_afgdt_derxy0_ui=fac_visceff_afgdt*derxy_(0,ui);
        const double fac_visceff_afgdt_derxy1_ui=fac_visceff_afgdt*derxy_(1,ui);
        const double fac_visceff_afgdt_derxy2_ui=fac_visceff_afgdt*derxy_(2,ui);
	      
        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);
        const double fac_gamma_dt_tauC_derxy_z_ui = fac_gamma_dt_tauC*derxy_(2,ui);
	      
        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
		
          /* add:                                                             */
          /* GALERKIN inertia term (intermediate) + convection (intermediate) */
          /* SUPG stabilisation --- inertia and convection                    */
          /* viscous term (intermediate), diagonal parts                      */
          const double sum =
            inertia_and_conv_ui*(funct_(vi)+conv_c_plus_svel_af_(vi))
            +
            fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
            +
            fac_visceff_afgdt_derxy1_ui*derxy_(1,vi)
            +
            fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
          
          elemat(fvi  ,fui  ) += sum+(fac_visceff_afgdt_derxy0_ui+fac_gamma_dt_tauC_derxy_x_ui)*derxy_(0,vi);
          elemat(fvi  ,fuip ) +=      fac_visceff_afgdt_derxy0_ui*derxy_(1,vi)+fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi);
          elemat(fvi  ,fuipp) +=      fac_visceff_afgdt_derxy0_ui*derxy_(2,vi)+fac_gamma_dt_tauC_derxy_z_ui*derxy_(0,vi);
          elemat(fvip ,fui  ) +=      fac_visceff_afgdt_derxy1_ui*derxy_(0,vi)+fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi);
          elemat(fvip ,fuip ) += sum+(fac_visceff_afgdt_derxy1_ui+fac_gamma_dt_tauC_derxy_y_ui)*derxy_(1,vi);
          elemat(fvip ,fuipp) +=      fac_visceff_afgdt_derxy1_ui*derxy_(2,vi)+fac_gamma_dt_tauC_derxy_z_ui*derxy_(1,vi);
          elemat(fvipp,fui  ) +=      fac_visceff_afgdt_derxy2_ui*derxy_(0,vi)+fac_gamma_dt_tauC_derxy_x_ui*derxy_(2,vi);
          elemat(fvipp,fuip ) +=      fac_visceff_afgdt_derxy2_ui*derxy_(1,vi)+fac_gamma_dt_tauC_derxy_y_ui*derxy_(2,vi);
          elemat(fvipp,fuipp) += sum+(fac_visceff_afgdt_derxy2_ui+fac_gamma_dt_tauC_derxy_z_ui)*derxy_(2,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int fuippp=4*ui+3;
	    
        const double fac_gamma_dt_derxy_0_ui = fac_gamma_dt*derxy_(0,ui);
        const double fac_gamma_dt_derxy_1_ui = fac_gamma_dt*derxy_(1,ui);
        const double fac_gamma_dt_derxy_2_ui = fac_gamma_dt*derxy_(2,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          int fvi =vi*4;
          
          /* SUPG stabilisation --- pressure     */
          /* factor: +tauM, rescaled by gamma*dt */

          elemat(fvi++,fuippp) += fac_gamma_dt_derxy_0_ui*conv_c_plus_svel_af_(vi);
          elemat(fvi++,fuippp) += fac_gamma_dt_derxy_1_ui*conv_c_plus_svel_af_(vi);
          elemat(fvi  ,fuippp) += fac_gamma_dt_derxy_2_ui*conv_c_plus_svel_af_(vi);
          
        } // vi
      } // ui

      if (higher_order_ele && newton!=Fluid3::minimal)
      {
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui  =ui*4 ;
          const int fuip =fui+1;
          const int fuipp=fui+2;
          
          /* SUPG stabilisation --- diffusion   */
          /* factor: -nu*alphaF*gamma*dt*tauM   */
          const double fac_visceff_afgdt_viscs2_0_ui=fac_visceff_afgdt*viscs2_(0,ui);
          const double fac_visceff_afgdt_viscs2_1_ui=fac_visceff_afgdt*viscs2_(1,ui);
          const double fac_visceff_afgdt_viscs2_2_ui=fac_visceff_afgdt*viscs2_(2,ui);
          const double fac_visceff_afgdt_derxy2_3_ui=fac_visceff_afgdt*derxy2_(3,ui);
          const double fac_visceff_afgdt_derxy2_4_ui=fac_visceff_afgdt*derxy2_(4,ui);
          const double fac_visceff_afgdt_derxy2_5_ui=fac_visceff_afgdt*derxy2_(5,ui);
          
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            int fvi  =vi*4 ;
            int fvip =fvi+1;
            int fvipp=fvi+2;
            
            elemat(fvi  ,fui  ) -= fac_visceff_afgdt_viscs2_0_ui*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuip ) -= fac_visceff_afgdt_derxy2_3_ui*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuipp) -= fac_visceff_afgdt_derxy2_4_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fui  ) -= fac_visceff_afgdt_derxy2_3_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuip ) -= fac_visceff_afgdt_viscs2_1_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuipp) -= fac_visceff_afgdt_derxy2_5_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fui  ) -= fac_visceff_afgdt_derxy2_4_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuip ) -= fac_visceff_afgdt_derxy2_5_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuipp) -= fac_visceff_afgdt_viscs2_2_ui*conv_c_plus_svel_af_(vi);
          } //end vi
        } // end ui
      }// end higher_order_ele and linearisation of viscous term

      //---------------------------------------------------------------
      //
      //                  GALERKIN AND SUPG PART
      //    REACTIVE TYPE LINEARISATIONS, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if (newton==Fluid3::Newton)
      {
        double temp[3][3];
        double supg_active_tauMresM[3];

        /* for linearisation of supg testfunction */
        supg_active_tauMresM[0]=supg_active_tauM*resM_(0);
        supg_active_tauMresM[1]=supg_active_tauM*resM_(1);
        supg_active_tauMresM[2]=supg_active_tauM*resM_(2);

        // loop rows (test functions for matrix)
        for (int vi=0; vi<iel; ++vi)
        {
          int fvi   =vi*4;
          int fvip  =fvi+1;
          int fvipp =fvi+2;

          /*  add linearised convective term in residual (supg), 
              linearisation of testfunction (supg)
              and linearised Galerkin term                */
          temp[0][0]=fac_afgdt*(supg_active_tauMresM[0]*derxy_(0,vi)+vderxyaf_(0,0)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[0][1]=fac_afgdt*(supg_active_tauMresM[0]*derxy_(1,vi)+vderxyaf_(0,1)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[0][2]=fac_afgdt*(supg_active_tauMresM[0]*derxy_(2,vi)+vderxyaf_(0,2)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[1][0]=fac_afgdt*(supg_active_tauMresM[1]*derxy_(0,vi)+vderxyaf_(1,0)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[1][1]=fac_afgdt*(supg_active_tauMresM[1]*derxy_(1,vi)+vderxyaf_(1,1)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[1][2]=fac_afgdt*(supg_active_tauMresM[1]*derxy_(2,vi)+vderxyaf_(1,2)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[2][0]=fac_afgdt*(supg_active_tauMresM[2]*derxy_(0,vi)+vderxyaf_(2,0)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[2][1]=fac_afgdt*(supg_active_tauMresM[2]*derxy_(1,vi)+vderxyaf_(2,1)*(conv_c_plus_svel_af_(vi)+funct_(vi)));
          temp[2][2]=fac_afgdt*(supg_active_tauMresM[2]*derxy_(2,vi)+vderxyaf_(2,2)*(conv_c_plus_svel_af_(vi)+funct_(vi)));

          // loop columns (solution for matrix, test function for vector)
          for (int ui=0; ui<iel; ++ui)
          {
            int fui=4*ui;

            elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][0]*funct_(ui);
            elemat(fvipp,fui++) += temp[2][0]*funct_(ui);
            elemat(fvi  ,fui  ) += temp[0][1]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][1]*funct_(ui);
            elemat(fvipp,fui++) += temp[2][1]*funct_(ui);
            elemat(fvi  ,fui  ) += temp[0][2]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][2]*funct_(ui);
            elemat(fvipp,fui  ) += temp[2][2]*funct_(ui);
          } // ui
        } // vi
      } // end newton

      //---------------------------------------------------------------
      //
      //      GALERKIN PART, CONTINUITY AND PRESSURE PART
      //                QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------

      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
      {
        const int fvi    =4*vi;
        const int fvip   =fvi+1;
        const int fvipp  =fvi+2;
	  
        const double fac_gamma_dt_derxy_0_vi = fac_gamma_dt*derxy_(0,vi);
        const double fac_gamma_dt_derxy_1_vi = fac_gamma_dt*derxy_(1,vi);
        const double fac_gamma_dt_derxy_2_vi = fac_gamma_dt*derxy_(2,vi);
	  
        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fuippp=4*ui+3;
	    
          /* GALERKIN pressure (implicit, rescaled to keep symmetry) */
	    
          /*  factor: -1, rescaled by gamma*dt

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
          */

          elemat(fvi  ,fuippp) -= fac_gamma_dt_derxy_0_vi*funct_(ui);
          elemat(fvip ,fuippp) -= fac_gamma_dt_derxy_1_vi*funct_(ui);
          elemat(fvipp,fuippp) -= fac_gamma_dt_derxy_2_vi*funct_(ui);

          /* GALERKIN continuity equation (implicit, transposed of above equation) */

          /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
          */

          elemat(fuippp,fvi  ) += fac_gamma_dt_derxy_0_vi*funct_(ui);
          elemat(fuippp,fvip ) += fac_gamma_dt_derxy_1_vi*funct_(ui);
          elemat(fuippp,fvipp) += fac_gamma_dt_derxy_2_vi*funct_(ui);
        } // ui
      } // vi

      //---------------------------------------------------------------
      //
      //             PSPG PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if(pspg == Fluid3::pstab_use_pspg)
      {
        const double fac_tauMp                   = fac*tauMp;
        const double fac_alphaM_tauMp            = fac_tauMp*alphaM;
        const double fac_gamma_dt_tauMp          = fac_tauMp*gamma*dt;
        const double fac_afgdt_tauMp             = fac_tauMp*afgdt;

        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_tauMp
            =
            fac*visceff*afgdt*tauMp;

          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =ui*4;
            const int fuip =fui+1;
            const int fuipp=fui+2;
	      
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
            const double pspg_diffusion_inertia_convect_2_ui
              =
              fac_visceff_afgdt_tauMp*viscs2_(2,ui)-fac_tauMp_inertia_and_conv;

            const double fac_visceff_afgdt_tauMp_derxy2_3_ui=fac_visceff_afgdt_tauMp*derxy2_(3,ui);
            const double fac_visceff_afgdt_tauMp_derxy2_4_ui=fac_visceff_afgdt_tauMp*derxy2_(4,ui);
            const double fac_visceff_afgdt_tauMp_derxy2_5_ui=fac_visceff_afgdt_tauMp*derxy2_(5,ui);
                           
	      
            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int fvippp=vi*4+3;

              elemat(fvippp,fui  ) -= 
                pspg_diffusion_inertia_convect_0_ui*derxy_(0,vi)
                +				                  
                fac_visceff_afgdt_tauMp_derxy2_3_ui*derxy_(1,vi)
                +				                  
                fac_visceff_afgdt_tauMp_derxy2_4_ui*derxy_(2,vi);
              elemat(fvippp,fuip ) -= 	                  
                fac_visceff_afgdt_tauMp_derxy2_3_ui*derxy_(0,vi)
                +				                  
                pspg_diffusion_inertia_convect_1_ui*derxy_(1,vi)
                +				                  
                fac_visceff_afgdt_tauMp_derxy2_5_ui*derxy_(2,vi);
              elemat(fvippp,fuipp) -=		                  
                fac_visceff_afgdt_tauMp_derxy2_4_ui*derxy_(0,vi)
                +				                  
                fac_visceff_afgdt_tauMp_derxy2_5_ui*derxy_(1,vi)
                +				                  
                pspg_diffusion_inertia_convect_2_ui*derxy_(2,vi);
            } // vi
          } // ui
        } // this is a higher order element and linearisation is not minimal
        else
        { // either this ain't a higher order element or a
          // linearisation of the viscous term is not necessary
          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =ui*4 ;
            const int fuip =fui+1;
            const int fuipp=fui+2;
	      
            const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));
	      
            for (int vi=0; vi<iel; ++vi)  // loop rows 
            {
              const int fvippp=vi*4+3;
		
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

              elemat(fvippp,fui  ) += fac_tauMp_inertia_and_conv*derxy_(0,vi) ;
              elemat(fvippp,fuip ) += fac_tauMp_inertia_and_conv*derxy_(1,vi) ;
              elemat(fvippp,fuipp) += fac_tauMp_inertia_and_conv*derxy_(2,vi) ;
            } // vi
          } // ui
        } // no linearisation of viscous part of residual is 
          // performed for pspg stabilisation cause either this 
          // ain't a higher order element or a linearisation of 
          // the viscous term is not necessary

        if (newton==Fluid3::Newton)
        {
          for (int vi=0; vi<iel; ++vi)  // loop rows 
          {
            int vidx = vi*4 + 3;
            double v1 = derxy_(0,vi)*vderxyaf_(0,0) + derxy_(1,vi)*vderxyaf_(1,0) + derxy_(2,vi)*vderxyaf_(2,0);
            double v2 = derxy_(0,vi)*vderxyaf_(0,1) + derxy_(1,vi)*vderxyaf_(1,1) + derxy_(2,vi)*vderxyaf_(2,1);
            double v3 = derxy_(0,vi)*vderxyaf_(0,2) + derxy_(1,vi)*vderxyaf_(1,2) + derxy_(2,vi)*vderxyaf_(2,2);
            for (int ui=0; ui<iel; ++ui) // loop columns
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
            } // ui
          } // vi
        } // end newton

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fuippp=ui*4+3;

          const double fac_gamma_dt_tauMp_derxy_0_ui=fac_gamma_dt_tauMp*derxy_(0,ui);
          const double fac_gamma_dt_tauMp_derxy_1_ui=fac_gamma_dt_tauMp*derxy_(1,ui);
          const double fac_gamma_dt_tauMp_derxy_2_ui=fac_gamma_dt_tauMp*derxy_(2,ui);

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

            elemat(vi*4+3,fuippp) += 
              fac_gamma_dt_tauMp_derxy_0_ui*derxy_(0,vi)
              +
              fac_gamma_dt_tauMp_derxy_1_ui*derxy_(1,vi)
              +
              fac_gamma_dt_tauMp_derxy_2_ui*derxy_(2,vi);
		
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
        if((vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_usfem)&&higher_order_ele)
        {
          const double fac_visc_tauMp_gamma_dt      = vstabfac*fac*visc*tauMp*gamma*dt;
          const double fac_visc_afgdt_tauMp         = vstabfac*fac*visc*afgdt*tauMp;
          const double fac_visc_alphaM_tauMp        = vstabfac*fac*visc*alphaM*tauMp;
          const double fac_visceff_visc_afgdt_tauMp = vstabfac*fac*visceff*visc*afgdt*tauMp;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    = ui*4;
            const int fuip   = fui+1;
            const int fuipp  = fui+2;
            const int fuippp = fui+3;

            const double acc_conv=(fac_visc_alphaM_tauMp*funct_(ui)
                                   +
                                   fac_visc_afgdt_tauMp*conv_c_af_(ui));

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   = vi*4;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;

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

              elemat(fvi  ,fui  ) += acc_conv*viscs2_(0,vi);
              elemat(fvi  ,fuip ) += acc_conv*derxy2_(3,vi);
              elemat(fvi  ,fuipp) += acc_conv*derxy2_(4,vi);
              elemat(fvip ,fui  ) += acc_conv*derxy2_(3,vi);
              elemat(fvip ,fuip ) += acc_conv*viscs2_(1,vi);
              elemat(fvip ,fuipp) += acc_conv*derxy2_(5,vi);
              elemat(fvipp,fui  ) += acc_conv*derxy2_(4,vi);
              elemat(fvipp,fuip ) += acc_conv*derxy2_(5,vi);
              elemat(fvipp,fuipp) += acc_conv*viscs2_(2,vi);

              /* viscous stabilisation --- diffusion  */

              /* factor: -(+)nu*nu*alphaF*gamma*dt*tauMp

                    /                                       \
                   |                 /    \                  |
                   |  2*nabla o eps | Dacc | , 2*div eps (v) |
                   |                 \    /                  |
                    \                                       /
              */
              elemat(fvi  ,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      derxy2_(4,ui)*derxy2_(4,vi));
              elemat(fvi  ,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(3,ui)
                                      +
                                      derxy2_(3,vi)*viscs2_(1,ui)
                                      +
                                      derxy2_(4,vi)*derxy2_(5,ui));
              elemat(fvi  ,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(4,ui)
                                      +
                                      derxy2_(3,vi)*derxy2_(5,ui)
                                      +
                                      derxy2_(4,vi)*viscs2_(2,ui));
              elemat(fvip ,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(3,vi)
                                      +
                                      derxy2_(3,ui)*viscs2_(1,vi)
                                      +
                                      derxy2_(4,ui)*derxy2_(5,vi));
              elemat(fvip ,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      viscs2_(1,ui)*viscs2_(1,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi));
              elemat(fvip ,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,vi)*derxy2_(4,ui)
                                      +
                                      viscs2_(1,vi)*derxy2_(5,ui)
                                      +
                                      derxy2_(5,vi)*viscs2_(2,ui));
              elemat(fvipp,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(4,vi)
                                      +
                                      derxy2_(3,ui)*derxy2_(5,vi)
                                      +
                                      derxy2_(4,ui)*viscs2_(2,vi));
              elemat(fvipp,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,ui)*derxy2_(4,vi)
                                      +
                                      viscs2_(1,ui)*derxy2_(5,vi)
                                      +
                                      derxy2_(5,ui)*viscs2_(2,vi));
              elemat(fvipp,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(4,ui)*derxy2_(4,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi)
                                      +
                                      viscs2_(2,ui)*viscs2_(2,vi));

              /* viscous stabilisation --- pressure   */

              /* factor: +(-)tauMp*nu, rescaled by gamma*dt

                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
                */
              elemat(fvi  ,fuippp) += fac_visc_tauMp_gamma_dt*
                                   (derxy_(0,ui)*viscs2_(0,vi)
                                    +
                                    derxy_(1,ui)*derxy2_(3,vi)
                                    +
                                    derxy_(2,ui)*derxy2_(4,vi)) ;
              elemat(fvip ,fuippp) += fac_visc_tauMp_gamma_dt*
                                   (derxy_(0,ui)*derxy2_(3,vi)
                                    +
                                    derxy_(1,ui)*viscs2_(1,vi)
                                    +
                                    derxy_(2,ui)*derxy2_(5,vi)) ;
              elemat(fvipp,fuippp) += fac_visc_tauMp_gamma_dt*
                                   (derxy_(0,ui)*derxy2_(4,vi)
                                    +
                                    derxy_(1,ui)*derxy2_(5,vi)
                                    +
                                    derxy_(2,ui)*viscs2_(2,vi)) ;

            } // vi
          } // ui
          
          if (newton==Fluid3::Newton)
          {
            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const int fui    = ui*4;
              const int fuip   = fui+1;
              const int fuipp  = fui+2;

              const double fac_visc_afgdt_tauMp_funct_ui=fac_visc_afgdt_tauMp*funct_(ui);

              for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
              {
                const int fvi   = vi*4;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;

                /* viscous stabilisation --- convection */

                /*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /


                */
                elemat(fvi  ,fui  ) += fac_visc_afgdt_tauMp_funct_ui*        
                                       (viscs2_(0,vi)*vderxyaf_(0,0)
                                        +
                                        derxy2_(3,vi)*vderxyaf_(1,0)
                                        +
                                        derxy2_(4,vi)*vderxyaf_(2,0));
                elemat(fvi  ,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,1)
                                        +
                                        derxy2_(3,vi)*vderxyaf_(1,1)
                                        +
                                        derxy2_(4,vi)*vderxyaf_(2,1));
                elemat(fvi  ,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,2)
                                        +
                                        derxy2_(3,vi)*vderxyaf_(1,2)
                                        +
                                        derxy2_(4,vi)*vderxyaf_(2,2));
                elemat(fvip ,fui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,0)
                                        +
                                        viscs2_(1,vi)*vderxyaf_(1,0)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(2,0));
                elemat(fvip ,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,1)
                                        +
                                        viscs2_(1,vi)*vderxyaf_(1,1)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(2,1));
                elemat(fvip ,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,2)
                                        +
                                        viscs2_(1,vi)*vderxyaf_(1,2)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(2,2));
                elemat(fvipp,fui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,0)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(1,0)
                                        +
                                        viscs2_(2,vi)*vderxyaf_(2,0));
                elemat(fvipp,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,1)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(1,1)
                                        +
                                        viscs2_(2,vi)*vderxyaf_(2,1));
                elemat(fvipp,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,2)
                                        +
                                        derxy2_(5,vi)*vderxyaf_(1,2)
                                        +
                                        viscs2_(2,vi)*vderxyaf_(2,2));
              } // vi
            } // ui
          } // end newton
        } // endif (a)gls
      }

      //---------------------------------------------------------------
      //
      //               QUASISTATIC STABILISATION PART
      //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
      //
      //---------------------------------------------------------------
      if(cross == Fluid3::cross_stress_stab)
      {
        const double fac_afgdt_tauM = fac*afgdt*tauM;
          
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
              
          const double fac_afgdt_tauM_conv_resM_ui = fac_afgdt_tauM*conv_resM_(ui);
              
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int fvi =4*vi;
            const double fac_afgdt_tauM_conv_resM_ui_funct_vi = fac_afgdt_tauM_conv_resM_ui*funct_(vi);
                
            /*  factor:
                    
                -alphaF*gamma*dt*tauM
                
                          /                          \
                         |  /            \            |
                         | | resM o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
            */
            elemat(fvi++,fui  ) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
            elemat(fvi++,fuip ) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
            elemat(fvi  ,fuipp) -= fac_afgdt_tauM_conv_resM_ui_funct_vi;
          } // vi
        } // ui
          
        const double fac_alphaM_tauM = fac*alphaM*tauM;

        double  am_nabla_u_afgdt_nabla_u_nabla_u[3][3];
          
        am_nabla_u_afgdt_nabla_u_nabla_u[0][0]=
          fac_alphaM_tauM*vderxyaf_(0,0)
          +
          fac_afgdt_tauM*(vderxyaf_(0,0)*vderxyaf_(0,0)
                          +
                          vderxyaf_(0,1)*vderxyaf_(1,0)
                          +
                          vderxyaf_(0,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][1]=
          fac_alphaM_tauM*vderxyaf_(0,1)
          +
          fac_afgdt_tauM*(vderxyaf_(0,0)*vderxyaf_(0,1)
                          +
                          vderxyaf_(0,1)*vderxyaf_(1,1)
                          +
                          vderxyaf_(0,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][2]=
          fac_alphaM_tauM*vderxyaf_(0,2)
          +
          fac_afgdt_tauM*(vderxyaf_(0,0)*vderxyaf_(0,2)
                          +
                          vderxyaf_(0,1)*vderxyaf_(1,2)
                          +
                          vderxyaf_(0,2)*vderxyaf_(2,2));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][0]=
          fac_alphaM_tauM*vderxyaf_(1,0)
          +
          fac_afgdt_tauM*(vderxyaf_(1,0)*vderxyaf_(0,0)
                          +
                          vderxyaf_(1,1)*vderxyaf_(1,0)
                          +
                          vderxyaf_(1,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][1]=
          fac_alphaM_tauM*vderxyaf_(1,1)
          +
          fac_afgdt_tauM*(vderxyaf_(1,0)*vderxyaf_(0,1)
                          +
                          vderxyaf_(1,1)*vderxyaf_(1,1)
                          +
                          vderxyaf_(1,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][2]=
          fac_alphaM_tauM*vderxyaf_(1,2)
          +
          fac_afgdt_tauM*(vderxyaf_(1,0)*vderxyaf_(0,2)
                          +
                          vderxyaf_(1,1)*vderxyaf_(1,2)
                          +
                          vderxyaf_(1,2)*vderxyaf_(2,2));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][0]=
          fac_alphaM_tauM*vderxyaf_(2,0)
          +
          fac_afgdt_tauM*(vderxyaf_(2,0)*vderxyaf_(0,0)
                          +
                          vderxyaf_(2,1)*vderxyaf_(1,0)
                          +
                          vderxyaf_(2,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][1]=
          fac_alphaM_tauM*vderxyaf_(2,1)
          +
          fac_afgdt_tauM*(vderxyaf_(2,0)*vderxyaf_(0,1)
                          +
                          vderxyaf_(2,1)*vderxyaf_(1,1)
                          +
                          vderxyaf_(2,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][2]=
          fac_alphaM_tauM*vderxyaf_(2,2)
          +
          fac_afgdt_tauM*(vderxyaf_(2,0)*vderxyaf_(0,2)
                          +
                          vderxyaf_(2,1)*vderxyaf_(1,2)
                          +
                          vderxyaf_(2,2)*vderxyaf_(2,2));
            
        double nabla_u[3][3];
        nabla_u[0][0]=fac_afgdt_tauM*vderxyaf_(0,0);
        nabla_u[0][1]=fac_afgdt_tauM*vderxyaf_(0,1);
        nabla_u[0][2]=fac_afgdt_tauM*vderxyaf_(0,2);
            
        nabla_u[1][0]=fac_afgdt_tauM*vderxyaf_(1,0);
        nabla_u[1][1]=fac_afgdt_tauM*vderxyaf_(1,1);
        nabla_u[1][2]=fac_afgdt_tauM*vderxyaf_(1,2);
            
        nabla_u[2][0]=fac_afgdt_tauM*vderxyaf_(2,0);
        nabla_u[2][1]=fac_afgdt_tauM*vderxyaf_(2,1);
        nabla_u[2][2]=fac_afgdt_tauM*vderxyaf_(2,2);
          
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
              
          const double u_nabla_ui = velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
              
          double coltemp[3][3];
            
          coltemp[0][0]=am_nabla_u_afgdt_nabla_u_nabla_u[0][0]*funct_(ui)+nabla_u[0][0]*u_nabla_ui;
          coltemp[0][1]=am_nabla_u_afgdt_nabla_u_nabla_u[0][1]*funct_(ui)+nabla_u[0][1]*u_nabla_ui;
          coltemp[0][2]=am_nabla_u_afgdt_nabla_u_nabla_u[0][2]*funct_(ui)+nabla_u[0][2]*u_nabla_ui;
              
          coltemp[1][0]=am_nabla_u_afgdt_nabla_u_nabla_u[1][0]*funct_(ui)+nabla_u[1][0]*u_nabla_ui;
          coltemp[1][1]=am_nabla_u_afgdt_nabla_u_nabla_u[1][1]*funct_(ui)+nabla_u[1][1]*u_nabla_ui;
          coltemp[1][2]=am_nabla_u_afgdt_nabla_u_nabla_u[1][2]*funct_(ui)+nabla_u[1][2]*u_nabla_ui;

          coltemp[2][0]=am_nabla_u_afgdt_nabla_u_nabla_u[2][0]*funct_(ui)+nabla_u[2][0]*u_nabla_ui;
          coltemp[2][1]=am_nabla_u_afgdt_nabla_u_nabla_u[2][1]*funct_(ui)+nabla_u[2][1]*u_nabla_ui;
          coltemp[2][2]=am_nabla_u_afgdt_nabla_u_nabla_u[2][2]*funct_(ui)+nabla_u[2][2]*u_nabla_ui;
            
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
                

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
              
            elemat(fvi  ,fui  ) -= funct_(vi)*coltemp[0][0];
            elemat(fvi  ,fuip ) -= funct_(vi)*coltemp[0][1];
            elemat(fvi  ,fuipp) -= funct_(vi)*coltemp[0][2];
            
            elemat(fvip ,fui  ) -= funct_(vi)*coltemp[1][0];
            elemat(fvip ,fuip ) -= funct_(vi)*coltemp[1][1];
            elemat(fvip ,fuipp) -= funct_(vi)*coltemp[1][2];
                
            elemat(fvipp,fui  ) -= funct_(vi)*coltemp[2][0];
            elemat(fvipp,fuip ) -= funct_(vi)*coltemp[2][1];
            elemat(fvipp,fuipp) -= funct_(vi)*coltemp[2][2];
          } // vi
        } // ui
            
        const double fac_gdt_tauM = fac*gamma*dt*tauM;
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fuippp =4*ui+3;
              
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

            /*  factor:

               -gamma*dt*tauM (rescaled for consistency)

                          /                               \
                         |  /                \   n+af      |
                         | | nabla Dp o nabla | u     , v  |
                         |  \                /             |
                          \                               /
            */
            elemat(fvi  ,fuippp) -= fac_gdt_tauM*funct_(vi)*(vderxyaf_(0,0)*derxy_(0,ui)+vderxyaf_(0,1)*derxy_(1,ui)+vderxyaf_(0,2)*derxy_(2,ui));
            elemat(fvip ,fuippp) -= fac_gdt_tauM*funct_(vi)*(vderxyaf_(1,0)*derxy_(0,ui)+vderxyaf_(1,1)*derxy_(1,ui)+vderxyaf_(1,2)*derxy_(2,ui));
            elemat(fvipp,fuippp) -= fac_gdt_tauM*funct_(vi)*(vderxyaf_(2,0)*derxy_(0,ui)+vderxyaf_(2,1)*derxy_(1,ui)+vderxyaf_(2,2)*derxy_(2,ui));
          } // vi
        } // ui
            
        if (higher_order_ele)
        {
          const double fac_visceff_afgdt_tauM = fac_afgdt_tauM*visceff;
              
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
                
            double coltemp[3][3];
                
            coltemp[0][0]=fac_visceff_afgdt_tauM*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(0,1)+derxy2_(4,ui)*vderxyaf_(0,2));
            coltemp[0][1]=fac_visceff_afgdt_tauM*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(0,0)+derxy2_(5,ui)*vderxyaf_(0,2));
            coltemp[0][2]=fac_visceff_afgdt_tauM*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(0,0)+derxy2_(5,ui)*vderxyaf_(0,1));

            coltemp[1][0]=fac_visceff_afgdt_tauM*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(1,1)+derxy2_(4,ui)*vderxyaf_(1,2));
            coltemp[1][1]=fac_visceff_afgdt_tauM*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(1,0)+derxy2_(5,ui)*vderxyaf_(1,2));
            coltemp[1][2]=fac_visceff_afgdt_tauM*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(1,0)+derxy2_(5,ui)*vderxyaf_(1,1));

            coltemp[2][0]=fac_visceff_afgdt_tauM*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(2,1)+derxy2_(4,ui)*vderxyaf_(2,2));
            coltemp[2][1]=fac_visceff_afgdt_tauM*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(2,0)+derxy2_(5,ui)*vderxyaf_(2,2));
            coltemp[2][2]=fac_visceff_afgdt_tauM*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(2,0)+derxy2_(5,ui)*vderxyaf_(2,1));

            for (int vi=0; vi<iel; ++vi)  // loop rows
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvi+2;
                
              /*  factor:
                  
                  +alphaF*gamma*dt*tauM

                          /                                               \
                         |  / /             /    \ \         \   n+af      |
                         | | | nabla o eps | Dacc | | o nabla | u     , v  |
                         |  \ \             \    / /         /             |
                          \                                               /
              */
              elemat(fvi  ,fui  ) += coltemp[0][0]*funct_(vi);
              elemat(fvi  ,fuip ) += coltemp[0][1]*funct_(vi);
              elemat(fvi  ,fuipp) += coltemp[0][2]*funct_(vi);

              elemat(fvip ,fui  ) += coltemp[1][0]*funct_(vi);
              elemat(fvip ,fuip ) += coltemp[1][1]*funct_(vi);
              elemat(fvip ,fuipp) += coltemp[1][2]*funct_(vi);

              elemat(fvipp,fui  ) += coltemp[2][0]*funct_(vi);
              elemat(fvipp,fuip ) += coltemp[2][1]*funct_(vi);
              elemat(fvipp,fuipp) += coltemp[2][2]*funct_(vi);
                  
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

    double fac_inertia_conv_and_bodyforce[3];
    fac_inertia_conv_and_bodyforce[0] = fac*(accintam_(0)+convaf_old_(0)-bodyforceaf_(0));
    fac_inertia_conv_and_bodyforce[1] = fac*(accintam_(1)+convaf_old_(1)-bodyforceaf_(1));
    fac_inertia_conv_and_bodyforce[2] = fac*(accintam_(2)+convaf_old_(2)-bodyforceaf_(2));
    
    if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
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
         resM_(1)*vderxyaf_(0,1)+
         resM_(2)*vderxyaf_(0,2));
      fac_inertia_conv_and_bodyforce[1]-=
        fac_tauM*
        (resM_(0)*vderxyaf_(1,0)+
         resM_(1)*vderxyaf_(1,1)+
         resM_(2)*vderxyaf_(1,2));
      fac_inertia_conv_and_bodyforce[2]-=
        fac_tauM*
        (resM_(0)*vderxyaf_(2,0)+
         resM_(1)*vderxyaf_(2,1)+
         resM_(2)*vderxyaf_(2,2));
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
	
    double viscous_and_pres[9];
    viscous_and_pres[0]=visceff_fac*vderxyaf_(0,0)*2.0-fac_prenp;
    viscous_and_pres[1]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_and_pres[2]=visceff_fac*(vderxyaf_(0,2)+vderxyaf_(2,0));
    viscous_and_pres[3]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_and_pres[4]=visceff_fac*vderxyaf_(1,1)*2.0-fac_prenp;
    viscous_and_pres[5]=visceff_fac*(vderxyaf_(1,2)+vderxyaf_(2,1));
    viscous_and_pres[6]=visceff_fac*(vderxyaf_(0,2)+vderxyaf_(2,0));
    viscous_and_pres[7]=visceff_fac*(vderxyaf_(1,2)+vderxyaf_(2,1));
    viscous_and_pres[8]=visceff_fac*vderxyaf_(2,2)*2.0-fac_prenp;
	
    if(reynolds == Fluid3::reynolds_stress_stab_only_rhs
       ||
       reynolds == Fluid3::reynolds_stress_stab)
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
      
      viscous_and_pres[0]-=fac_tauM_tauM_resM_0  *resM_(0);
      viscous_and_pres[1]-=fac_tauM_tauM_resM_0  *resM_(1);
      viscous_and_pres[2]-=fac_tauM_tauM_resM_0  *resM_(2);
      viscous_and_pres[3]-=fac_tauM_tauM_resM_0  *resM_(1);
      viscous_and_pres[4]-=fac_tauM_tauM_resM_1  *resM_(1);
      viscous_and_pres[5]-=fac_tauM_tauM_resM_1  *resM_(2);
      viscous_and_pres[6]-=fac_tauM_tauM_resM_0  *resM_(2);
      viscous_and_pres[7]-=fac_tauM_tauM_resM_1  *resM_(2);
      viscous_and_pres[8]-=fac_tauM_tauM*resM_(2)*resM_(2);
    }

    /* continuity equation, factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
    */
    const double fac_divunp  = fac*divunp_;
    
    for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
    {
      int fui=4*ui;
      /* inertia, convective and dead load, cross terms fith 
         funct                                              */
      /* viscous, pressure, reynolds, cstab terms with 
         derxy                                              */
      
      elevec(fui++) -= 
        fac_inertia_conv_and_bodyforce[0]*funct_(ui)
        +
        derxy_(0,ui)*viscous_and_pres[0]
        +
        derxy_(1,ui)*viscous_and_pres[1]
        +
        derxy_(2,ui)*viscous_and_pres[2];
      elevec(fui++) -= 
        fac_inertia_conv_and_bodyforce[1]*funct_(ui)
        +
        derxy_(0,ui)*viscous_and_pres[3]
        +
        derxy_(1,ui)*viscous_and_pres[4]
        +
        derxy_(2,ui)*viscous_and_pres[5];
      elevec(fui++) -= 
        fac_inertia_conv_and_bodyforce[2]*funct_(ui)
        +
        derxy_(0,ui)*viscous_and_pres[6]
        +
        derxy_(1,ui)*viscous_and_pres[7]
        +
        derxy_(2,ui)*viscous_and_pres[8];
	  
      /* continuity equation */
      elevec(fui  ) -= fac_divunp*funct_(ui);
    }

    if(pspg == Fluid3::pstab_use_pspg)
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
      
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
        elevec(4*ui+3)-=fac_tauMp*(resM_(0)*derxy_(0,ui)+resM_(1)*derxy_(1,ui)+resM_(2)*derxy_(2,ui));
      }
    } // end pspg

    if(supg == Fluid3::convective_stab_supg)
    {
      const double fac_tauM = fac*supg_active_tauM;
      
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
        int fui=4*ui;
        
        const double fac_tauM_conv_c_af_ui = fac_tauM*conv_c_af_(ui);
        /*
          factor: +tauM
          
          SUPG stabilisation


                  /                             \
                 |   n+af    / n+af        \     |
                 |  r     , | c     o nabla | v  |
                 |   M       \             /     |
                  \                             /
        */

        elevec(fui++) -= fac_tauM_conv_c_af_ui*resM_(0);
        elevec(fui++) -= fac_tauM_conv_c_af_ui*resM_(1);
        elevec(fui  ) -= fac_tauM_conv_c_af_ui*resM_(2);
      } // end loop rows
    } // end supg
    
    if (higher_order_ele)
    {
      if(vstab != Fluid3::viscous_stab_none && higher_order_ele)
      {
        const double fac_visc_tauMp = vstabfac * fac*visc*tauMp;
        
        for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
        {
          int fui=4*ui;
          /*
              factor: -(+)tauMp*nu

              viscous stabilisation --- inertia


                 /                      \
                |   n+af                 |
                |  r    , 2*div eps (v)  |
                |   M                    |
                 \                      /

          */
          elevec(fui++) -= fac_visc_tauMp*
                           (resM_(0)*viscs2_(0,ui)
                            +
                            resM_(1)*derxy2_(3,ui)
                            +
                            resM_(2)*derxy2_(4,ui)) ;
          elevec(fui++) -= fac_visc_tauMp*
                           (resM_(0)*derxy2_(3,ui)
                            +
                            resM_(1)*viscs2_(1,ui)
                            +
                            resM_(2)*derxy2_(5,ui)) ;
          elevec(fui  ) -= fac_visc_tauMp*
                           (resM_(0)*derxy2_(4,ui)
                            +
                            resM_(1)*derxy2_(5,ui)
                            +
                            resM_(2)*viscs2_(2,ui)) ;
        } // end loop rows ui
      } // endif (a)gls
    } // hoel

    if(fssgv != Fluid3::fssgv_no)
    {
      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)
      
      for (int ui=0; ui<iel; ++ui)
      {
        /* fine-scale subgrid-viscosity term on right hand side */
        /*
                                  /                              \
                         n+af    |       /    n+af\         / \   |
             - nu_art(fsu    ) * |  eps | Dfsu     | , eps | v |  |
                                 |       \        /         \ /   |
                                  \                              /
        */
        elevec(ui*4    ) -= vartfac*( 2.0*derxy_(0, ui)*fsvderxyaf_(0, 0)
                                     +    derxy_(1, ui)*fsvderxyaf_(0, 1)
                                     +    derxy_(1, ui)*fsvderxyaf_(1, 0)
                                     +    derxy_(2, ui)*fsvderxyaf_(0, 2)
                                     +    derxy_(2, ui)*fsvderxyaf_(2, 0));
        elevec(ui*4 + 1) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 1)
                                     +    derxy_(0, ui)*fsvderxyaf_(1, 0)
                                     +2.0*derxy_(1, ui)*fsvderxyaf_(1, 1)
                                     +    derxy_(2, ui)*fsvderxyaf_(1, 2)
                                     +    derxy_(2, ui)*fsvderxyaf_(2, 1));
        elevec(ui*4 + 2) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 2)
                                     +    derxy_(0, ui)*fsvderxyaf_(2, 0)
                                     +    derxy_(1, ui)*fsvderxyaf_(1, 2)
                                     +    derxy_(1, ui)*fsvderxyaf_(2, 1)
                                     +2.0*derxy_(2, ui)*fsvderxyaf_(2, 2));
      } // end loop ui
    } // end not fssgv_no
  } // end loop iquad
  return;
} //Sysmat_adv_qs


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  advective version using time dependent subgrid scales              |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Sysmat_adv_td(
  Fluid3*                                 ele             ,
  std::vector<blitz::Array<double,1> >&   myknots         ,
  LINALG::Matrix<4*iel,4*iel>&            elemat          ,
  LINALG::Matrix<4*iel,1>&                elevec          ,
  const LINALG::Matrix<3,iel>&            edispnp         ,
  const LINALG::Matrix<3,iel>&            egridvaf        ,
  const LINALG::Matrix<3,iel>&            evelnp          ,
  const LINALG::Matrix<iel,1>&            eprenp          ,
  const LINALG::Matrix<3,iel>&            eaccam          ,
  const LINALG::Matrix<3,iel>&            evelaf          ,
  const LINALG::Matrix<3,iel>&            fsevelaf        ,
  const struct _MATERIAL*                 material        ,
  const double                            alphaM          ,
  const double                            alphaF          ,
  const double                            gamma           ,
  const double                            dt              ,
  const double                            time            ,
  const enum Fluid3::LinearisationAction  newton          ,
  const bool                              higher_order_ele,
  const enum Fluid3::StabilisationAction  fssgv           ,
  const enum Fluid3::StabilisationAction  inertia         ,
  const enum Fluid3::StabilisationAction  pspg            ,
  const enum Fluid3::StabilisationAction  supg            ,
  const enum Fluid3::StabilisationAction  vstab           ,
  const enum Fluid3::StabilisationAction  cstab           ,
  const enum Fluid3::StabilisationAction  cross           ,
  const enum Fluid3::StabilisationAction  reynolds        ,
  const enum Fluid3::TauType              whichtau        ,
  const enum Fluid3::TurbModelAction      turb_mod_action ,
  double&                                 Cs              ,
  double&                                 Cs_delta_sq     ,
  double&                                 visceff         ,
  const double                            l_tau           ,
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
  if (vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_gls_only_rhs)
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

  SetElementData(ele            ,
                 edispnp        ,
                 evelaf         ,
                 fsevelaf       ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 fssgv          ,
                 turb_mod_action,
                 l_tau          ,
                 Cs             ,
                 Cs_delta_sq    ,
                 visceff        );
 
  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero
  if(ele->sub_acc_old_.extent(blitz::firstDim) != 3 
     ||
     ele->sub_acc_old_.extent(blitz::secondDim) != intpoints.nquad)
  {
    ele->sub_acc_old_ .resize(3,intpoints.nquad);
    ele->sub_acc_old_  = 0.;
  }
  if(ele->sub_vel_old_.extent(blitz::firstDim) != 3 
     || 
     ele->sub_vel_old_.extent(blitz::secondDim) != intpoints.nquad)
  {
    ele->sub_vel_old_ .resize(3,intpoints.nquad);
    ele->sub_vel_old_  = 0.;
    
    ele->sub_vel_.resize(3,intpoints.nquad);
    ele->sub_vel_ = 0.;
  }
  
  // get subscale information from element --- this is just a reference
  // to the element data
  blitz::Array<double,2> saccn (ele->sub_acc_old_);
  blitz::Array<double,2> sveln (ele->sub_vel_old_);
  blitz::Array<double,2> svelnp(ele->sub_vel_    );

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second 
    // derivatives in a gausspoint and integration weight including 
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);

    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,Fluid3::subscales_time_dependent,gamma,dt,hk,mk,visceff);

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

    if(cstab == Fluid3::continuity_stab_none)
    {
      tau_(2)=0.0;
    }

    const double tauC   = tau_(2);

    double supg_active;
    if(supg == Fluid3::convective_stab_supg)
    {
      supg_active=1.0;
    }
    else
    {
      supg_active=0.0;
    }

    // update estimates for the subscale quantities
    const double facMtau = 1./(alphaM*tauM+afgdt);

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
      svelaf_(rr)=alphaF*svelnp(rr,iquad)+(1.0-alphaF)*sveln(rr,iquad);
    }

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
    if(cross    == Fluid3::cross_stress_stab
       ||
       reynolds == Fluid3::reynolds_stress_stab)
    {
      for (int rr=0;rr<iel;++rr)
      {
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);
	
        for (int mm=1;mm<3;++mm)
        {
          conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
        }
      }

      if(reynolds == Fluid3::reynolds_stress_stab)
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
    if(cross == Fluid3::cross_stress_stab_only_rhs
       ||
       cross == Fluid3::cross_stress_stab
      )
    {
      for (int rr=0;rr<3;++rr)
      {
        convsubaf_old_(rr) = vderxyaf_(rr, 0)*svelaf_(0);
	
        for (int mm=1;mm<3;++mm)
        {
          convsubaf_old_(rr) += vderxyaf_(rr, mm)*svelaf_(mm);
        }
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
	  
      if(inertia == Fluid3::inertia_stab_keep
         || 
         inertia == Fluid3::inertia_stab_keep_complete)
      {
        // rescale time factors terms affected by inertia stabilisation
        fac_inertia   *=afgdt*facMtau;
        fac_convection*=afgdt*facMtau;
	  
        // do inertia stabilisation terms which are not scaled
        // Galerkin terms since they are not partially integrated

        const double fac_alphaM_tauM_facMtau = fac*alphaM*tauM*facMtau;
	  
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
	    
          const double fac_alphaM_gamma_dt_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*gamma*dt*funct_(vi);
	    
          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fuippp  =4*ui+3;
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
	      
            elemat(fvi  ,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(0,ui);
            elemat(fvip ,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(1,ui);
            elemat(fvipp,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(2,ui);
          } // ui
        } // vi

        if(higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_alphaM_tauM_facMtau
            =
            fac*visceff*afgdt*alphaM*tauM*facMtau;
	    
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi    =4*vi;
            const int fvip   =fvi+1;
            const int fvipp  =fvi+2;
		
            const double fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
              =
              fac_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);
		  
            for (int ui=0; ui<iel; ++ui) // loop columns
            {	  
              const int fui    =4*ui;
              const int fuip   =fui+1;
              const int fuipp  =fui+2;
	  
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
              const double a = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(3,ui);
              const double b = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(4,ui);
              const double c = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(5,ui);

              elemat(fvi  ,fui  ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,ui);
              elemat(fvi  ,fuip ) += a;
              elemat(fvi  ,fuipp) += b;
              elemat(fvip ,fui  ) += a;
              elemat(fvip ,fuip ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,ui);
              elemat(fvip ,fuipp) += c;
              elemat(fvipp,fui  ) += b;
              elemat(fvipp,fuip ) += c;
              elemat(fvipp,fuipp) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(2,ui);
            } // ui
          } // vi
        } // end higher order element and  linearisation of linear terms not supressed

        if(inertia == Fluid3::inertia_stab_keep_complete)
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
            for (int rr=0;rr<3;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=velintaf_(rr)/norm;
            }
          }
          else
          {
            normed_velintaf_(0) = 0.0;
            for (int rr=1;rr<3;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=0.0;
            }      
          }
	      
          double temp=0.0;
          if(whichtau==Fluid3::codina)
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
          else if(whichtau==Fluid3::smoothed_franca_barrenechea_valentin_wall)
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
          else if(whichtau==Fluid3::franca_barrenechea_valentin_wall)
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
            const int fvi    =4*vi;
            const int fvip   =fvi+1;
            const int fvipp  =fvi+2;

              
            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int fui  =4*ui;
              const int fuip =fui+1;
              const int fuipp=fui+2;

              elemat(fvi  ,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(fvi  ,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(fvi  ,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(0);

              elemat(fvip ,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(fvip ,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(fvip ,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(1);

              elemat(fvipp,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(2);
              elemat(fvipp,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(2);
              elemat(fvipp,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(2);
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
        const int fui  =4*ui;
        const int fuip =fui+1;
        const int fuipp=fui+2;
        
        /* GALERKIN inertia term (intermediate) + convection (intermediate) */
        const double inertia_and_conv_ui 
          = 
          fac_inertia*funct_(ui)
          +
          fac_convection*conv_c_af_(ui);
	      
        /* viscous term (intermediate), 'diagonal' parts */
        const double visc_0=fac_afgdt_visceff*derxy_(0,ui);
        const double visc_1=fac_afgdt_visceff*derxy_(1,ui);
        const double visc_2=fac_afgdt_visceff*derxy_(2,ui);
	    
        /* SUPG stabilisation --- inertia and convection */
        const double supg_inertia_and_conv_ui
          =
          fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_afgdt_tauM_facMtau*conv_c_af_(ui);
        
        /* CSTAB entries */
        const double cstab_0 = cstabfac*derxy_(0,ui);
        const double cstab_1 = cstabfac*derxy_(1,ui);
        const double cstab_2 = cstabfac*derxy_(2,ui);
        
        /* combined CSTAB/viscous entires */
        const double visc_and_cstab_0 = visc_0+cstab_0;
        const double visc_and_cstab_1 = visc_1+cstab_1;
        const double visc_and_cstab_2 = visc_2+cstab_2;
	      
        for (int vi=0; vi<iel; ++vi)  // loop rows 
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
	  
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
            visc_1*derxy_(1,vi)
            +
            visc_2*derxy_(2,vi);
		  
          /* CONTINUITY stabilisation                     */
          /* viscous term (intermediate, remaining parts) */

          const double a=visc_0*derxy_(1,vi)+cstab_1*derxy_(0,vi);
          const double b=visc_0*derxy_(2,vi)+cstab_2*derxy_(0,vi);
          const double c=visc_2*derxy_(1,vi)+cstab_1*derxy_(2,vi);

          elemat(fvi  ,fui  ) += sum+visc_and_cstab_0*derxy_(0,vi);
          elemat(fvip ,fuip ) += sum+visc_and_cstab_1*derxy_(1,vi);
          elemat(fvipp,fuipp) += sum+visc_and_cstab_2*derxy_(2,vi);
	      
          elemat(fvi  ,fuip ) += a;
          elemat(fvi  ,fuipp) += b;
          elemat(fvipp,fuip ) += c;
          elemat(fuip ,fvipp) += c;
          elemat(fuipp,fvi  ) += b;
          elemat(fuip ,fvi  ) += a;
        } // vi
      } // ui
      
      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int fuippp  =4*ui+3;

        const double fac_gamma_dt_funct_ui=fac_gamma_dt*funct_(ui);
	  
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi  =4*vi;
          const int fvip =fvi+1;
          const int fvipp=fvi+2;

          /* GALERKIN pressure   (implicit), rescaled by gamma*dt */
          /* continuity equation (implicit)                       */

          elemat(fvi   ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(fvip  ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(1,vi);
          elemat(fvipp ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(2,vi);

          elemat(fuippp,fvi   ) += fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(fuippp,fvip  ) += fac_gamma_dt_funct_ui*derxy_(1,vi);
          elemat(fuippp,fvipp ) += fac_gamma_dt_funct_ui*derxy_(2,vi);
        } // vi
      } // ui

      if (newton==Fluid3::Newton) // if newton and supg
      {
        const double fac_afgdt_afgdt_tauM_facMtau = fac*afgdt*afgdt*facMtau*tauM;

        // linearisation of SUPG testfunction
        double temp[3][3];

        const double fac_afgdt_svelaf_0 = fac*afgdt*supg_active*svelaf_(0);
        const double fac_afgdt_svelaf_1 = fac*afgdt*supg_active*svelaf_(1);
        const double fac_afgdt_svelaf_2 = fac*afgdt*supg_active*svelaf_(2);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi  =4*vi;
          const int fvip =fvi+1;
          const int fvipp=fvi+2;
            
          // linearisations of reactive Galerkin part (remaining after inertia_stab)
          // and SUPG part (reactive part from residual)
          const double scaled_inertia_and_conv_vi 
            = 
            fac_convection*funct_(vi)
            +
            fac_afgdt_afgdt_tauM_facMtau*conv_c_plus_svel_af_(vi);
	      
          temp[0][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,0)-fac_afgdt_svelaf_0*derxy_(0,vi);
          temp[1][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,1)-fac_afgdt_svelaf_0*derxy_(1,vi);
          temp[2][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,2)-fac_afgdt_svelaf_0*derxy_(2,vi);
          temp[0][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,0)-fac_afgdt_svelaf_1*derxy_(0,vi);
          temp[1][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,1)-fac_afgdt_svelaf_1*derxy_(1,vi);
          temp[2][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,2)-fac_afgdt_svelaf_1*derxy_(2,vi);
          temp[0][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,0)-fac_afgdt_svelaf_2*derxy_(0,vi);
          temp[1][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,1)-fac_afgdt_svelaf_2*derxy_(1,vi);
          temp[2][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,2)-fac_afgdt_svelaf_2*derxy_(2,vi);
	      
          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
            elemat(fvi  ,fuip ) += temp[1][0]*funct_(ui);
            elemat(fvi  ,fuipp) += temp[2][0]*funct_(ui);
            elemat(fvip ,fui  ) += temp[0][1]*funct_(ui);
            elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);
            elemat(fvip ,fuipp) += temp[2][1]*funct_(ui);
            elemat(fvipp,fui  ) += temp[0][2]*funct_(ui);
            elemat(fvipp,fuip ) += temp[1][2]*funct_(ui);
            elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);
          } // ui
        } // vi
      } // end if newton

      for (int ui=0; ui<iel; ++ui) // loop columns 
      {
        const int fuippp=4*ui+3;

        const double scaled_gradp_0 = fac_gdt_afgdt_tauM_facMtau*derxy_(0,ui);
        const double scaled_gradp_1 = fac_gdt_afgdt_tauM_facMtau*derxy_(1,ui);
        const double scaled_gradp_2 = fac_gdt_afgdt_tauM_facMtau*derxy_(2,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows 
        {
          const int fvi=4*vi;
 
          /* SUPG stabilisation --- pressure, rescaled by gamma*dt */
          elemat(fvi  ,fuippp) += scaled_gradp_0*conv_c_plus_svel_af_(vi);
          elemat(fvi+1,fuippp) += scaled_gradp_1*conv_c_plus_svel_af_(vi);
          elemat(fvi+2,fuippp) += scaled_gradp_2*conv_c_plus_svel_af_(vi);
        } // vi
      } // ui

      if(higher_order_ele && newton!=Fluid3::minimal)
      {
        const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
    
          double coltemp[3][3];
              
          coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(0,ui);
          coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(3,ui);
          coltemp[0][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(4,ui);
                                                                                                                     
          coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(3,ui);
          coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(1,ui);
          coltemp[1][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(5,ui);
                                                                                                                     
          coltemp[2][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(4,ui);
          coltemp[2][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(5,ui);
          coltemp[2][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(2,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

            /*  SUPG stabilisation, diffusion */
            elemat(fvi  ,fui  ) -= coltemp[0][0]*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuip ) -= coltemp[0][1]*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuipp) -= coltemp[0][2]*conv_c_plus_svel_af_(vi);
            
            elemat(fvip ,fui  ) -= coltemp[1][0]*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuip ) -= coltemp[1][1]*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuipp) -= coltemp[1][2]*conv_c_plus_svel_af_(vi);
            
            elemat(fvipp,fui  ) -= coltemp[2][0]*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuip ) -= coltemp[2][1]*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuipp) -= coltemp[2][2]*conv_c_plus_svel_af_(vi);
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
      if(pspg == Fluid3::pstab_use_pspg)
      {
        const double fac_afgdt_gamma_dt_tauM_facMtau  = fac*afgdt*gamma*dt*tauM*facMtau;
        const double fac_gdt_gdt_tauM_facMtau         = fac*gamma*dt*tauM*facMtau*gamma*dt;
        const double fac_alphaM_gamma_dt_tauM_facMtau = fac*alphaM*gamma*dt*tauM*facMtau;

        if(higher_order_ele  && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_gamma_dt_tauM_facMtau
            =
            fac*visceff*afgdt*gamma*dt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);


            const double pspg_diffusion_inertia_convect_0_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(0,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_1_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(1,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_2_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(2,ui)-inertia_and_conv_ui;

            const double scaled_derxy2_3_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(3,ui);
            const double scaled_derxy2_4_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(4,ui);
            const double scaled_derxy2_5_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(5,ui);
                 
            for (int vi=0; vi<iel; ++vi)  // loop rows 
            {
              const int fvippp =4*vi+3;
	      
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

              elemat(fvippp,fui  ) -= 
                derxy_(0,vi)*pspg_diffusion_inertia_convect_0_ui
                +
                derxy_(1,vi)*scaled_derxy2_3_ui
                +
                derxy_(2,vi)*scaled_derxy2_4_ui;
              elemat(fvippp,fuip ) -= 
                derxy_(0,vi)*scaled_derxy2_3_ui
                +
                derxy_(1,vi)*pspg_diffusion_inertia_convect_1_ui
                +
                derxy_(2,vi)*scaled_derxy2_5_ui;
              elemat(fvippp,fuipp) -= 
                derxy_(0,vi)*scaled_derxy2_4_ui
                +
                derxy_(1,vi)*scaled_derxy2_5_ui
                +
                derxy_(2,vi)*pspg_diffusion_inertia_convect_2_ui;
            }
          }
        }
        else
        {
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi) // loop rows
            {
              const int fvippp =4*vi+3;
	      
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

              elemat(fvippp,fui  ) +=derxy_(0,vi)*inertia_and_conv_ui;
              elemat(fvippp,fuip ) +=derxy_(1,vi)*inertia_and_conv_ui;
              elemat(fvippp,fuipp) +=derxy_(2,vi)*inertia_and_conv_ui;
            }
          }
        } // neglect viscous linearisations, do just inertia and convective

        for (int ui=0; ui<iel; ++ui) // loop columns 
        {
          const int fuippp=4*ui+3;
          const double scaled_derxy_0=fac_gdt_gdt_tauM_facMtau*derxy_(0,ui);
          const double scaled_derxy_1=fac_gdt_gdt_tauM_facMtau*derxy_(1,ui);
          const double scaled_derxy_2=fac_gdt_gdt_tauM_facMtau*derxy_(2,ui);

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

            elemat(vi*4+3,fuippp) += 
              (scaled_derxy_0*derxy_(0,vi)
               +
               scaled_derxy_1*derxy_(1,vi)
               +
               scaled_derxy_2*derxy_(2,vi)) ;

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (newton==Fluid3::Newton) // if pspg and newton
        {

          for (int vi=0; vi<iel; ++vi) // loop columns
          {
            const int fvippp=4*vi+3;

            const double a=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0)+derxy_(2,vi)*vderxyaf_(2,0));
            const double b=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1)+derxy_(2,vi)*vderxyaf_(2,1));
            const double c=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,2)+derxy_(1,vi)*vderxyaf_(1,2)+derxy_(2,vi)*vderxyaf_(2,2));

            for (int ui=0; ui<iel; ++ui)  // loop rows 
            {
              const int fui=4*ui;
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

              elemat(fvippp,fui  ) += a*funct_(ui);
              elemat(fvippp,fui+1) += b*funct_(ui);
              elemat(fvippp,fui+2) += c*funct_(ui);
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
        if(vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_gls)
        {
          const double tauMqs = afgdt*tauM*facMtau;

          const double fac_visc_tauMqs_alphaM        = vstabfac*fac*visc*tauMqs*alphaM;
          const double fac_visc_tauMqs_afgdt         = vstabfac*fac*visc*tauMqs*afgdt;
          const double fac_visc_tauMqs_afgdt_visceff = vstabfac*fac*visc*tauMqs*afgdt*visceff;
          const double fac_visc_tauMqs_gamma_dt      = vstabfac*fac*visc*tauMqs*gamma*dt;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;

            const double inertia_and_conv
              =
              fac_visc_tauMqs_alphaM*funct_(ui)+fac_visc_tauMqs_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;
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


              const double a = inertia_and_conv*derxy2_(3,vi);
              const double b = inertia_and_conv*derxy2_(4,vi);
              const double c = inertia_and_conv*derxy2_(5,vi);

              elemat(fvi  ,fui  ) += inertia_and_conv*viscs2_(0,vi);
              elemat(fvi  ,fuip ) += a;
              elemat(fvi  ,fuipp) += b;
              elemat(fvip ,fui  ) += a;
              elemat(fvip ,fuip ) += inertia_and_conv*viscs2_(1,vi);
              elemat(fvip ,fuipp) += c;
              elemat(fvipp,fui  ) += b;
              elemat(fvipp,fuip ) += c;
              elemat(fvipp,fuipp) += inertia_and_conv*viscs2_(2,vi);
            }
          }

          for (int ui=0;ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;

            for (int vi=0;vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;

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
                               (viscs2_(0,vi)*derxy2_(3,ui)
                                +
                                derxy2_(3,vi)*viscs2_(1,ui)
                                +
                                derxy2_(4,vi)*derxy2_(5,ui));
              
              elemat(fvi  ,fuip ) -= a; 
              elemat(fuip ,fvi  ) -= a;
		
              const double b = fac_visc_tauMqs_afgdt_visceff*
                               (viscs2_(0,ui)*derxy2_(4,vi)
                                +
                                derxy2_(3,ui)*derxy2_(5,vi)
                                +
                                derxy2_(4,ui)*viscs2_(2,vi));
              
              elemat(fvipp,fui  ) -= b;
              elemat(fui  ,fvipp) -= b;
              
              const double c = fac_visc_tauMqs_afgdt_visceff*
                               (derxy2_(3,ui)*derxy2_(4,vi)
                                +
                                viscs2_(1,ui)*derxy2_(5,vi)
                                +
                                derxy2_(5,ui)*viscs2_(2,vi));
              
              elemat(fvipp,fuip ) -= c;
              elemat(fuip ,fvipp) -= c;

              elemat(fvi   ,fui ) -= fac_visc_tauMqs_afgdt_visceff*
		                     (viscs2_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      derxy2_(4,ui)*derxy2_(4,vi));

              elemat(fvip ,fuip ) -= fac_visc_tauMqs_afgdt_visceff*
                                     (derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      viscs2_(1,ui)*viscs2_(1,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi));

              elemat(fvipp,fuipp) -= fac_visc_tauMqs_afgdt_visceff*
                                     (derxy2_(4,ui)*derxy2_(4,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi)
                                      +
                                      viscs2_(2,ui)*viscs2_(2,vi));
              
            }
          }

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;
            const int fuippp =fuipp+1;
            
            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;

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
              elemat(fvi  ,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*viscs2_(0,vi)
                                       +
                                       derxy_(1,ui)*derxy2_(3,vi)
                                       +
                                       derxy_(2,ui)*derxy2_(4,vi)) ;
              elemat(fvip ,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(3,vi)
                                       +
                                       derxy_(1,ui)*viscs2_(1,vi)
                                       +
                                       derxy_(2,ui)*derxy2_(5,vi)) ;
              elemat(fvipp,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(4,vi)
                                       +
                                       derxy_(1,ui)*derxy2_(5,vi)
                                       +
                                       derxy_(2,ui)*viscs2_(2,vi)) ;
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton==Fluid3::Newton)
          {
            
            double temp[3][3];
            for (int vi=0; vi<iel; ++vi) // loop columns (solution for matrix, test function for vector)
            {
              const int fvi    =4*vi;
              const int fvip   =fvi+1;
              const int fvipp  =fvip+1;

              temp[0][0]=(viscs2_(0,vi)*vderxyaf_(0,0)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,0)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][0]=(viscs2_(0,vi)*vderxyaf_(0,1)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,1)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][0]=(viscs2_(0,vi)*vderxyaf_(0,2)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,2)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;
              temp[0][1]=(derxy2_(3,vi)*vderxyaf_(0,0)
			  +
                          viscs2_(1,vi)*vderxyaf_(1,0)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][1]=(derxy2_(3,vi)*vderxyaf_(0,1)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,1)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][1]=(derxy2_(3,vi)*vderxyaf_(0,2)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,2)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;
              temp[0][2]=(derxy2_(4,vi)*vderxyaf_(0,0)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,0)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][2]=(derxy2_(4,vi)*vderxyaf_(0,1)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,1)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][2]=(derxy2_(4,vi)*vderxyaf_(0,2)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,2)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;

              for (int ui=0; ui<iel; ++ui)  // loop rows (test functions for matrix)
              {
                const int fui    =4*ui;
                const int fuip   =fui+1;
                const int fuipp  =fuip+1;

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
                elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
                elemat(fvi  ,fuip ) += temp[1][0]*funct_(ui);                                     
                elemat(fvi  ,fuipp) += temp[2][0]*funct_(ui);                                     
                elemat(fvip ,fui  ) += temp[0][1]*funct_(ui);                                       
                elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);                                   
                elemat(fvip ,fuipp) += temp[2][1]*funct_(ui);                                   
                elemat(fvipp,fui  ) += temp[0][2]*funct_(ui);                                   
                elemat(fvipp,fuip ) += temp[1][2]*funct_(ui);		                         
                elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);
                                         
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
      if(cross == Fluid3::cross_stress_stab)
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

            elemat(vi*4    , ui*4    ) += fac_afgdt_conv_subaf_ui_funct_vi;
            elemat(vi*4 + 1, ui*4 + 1) += fac_afgdt_conv_subaf_ui_funct_vi;
            elemat(vi*4 + 2, ui*4 + 2) += fac_afgdt_conv_subaf_ui_funct_vi;
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

        double  am_nabla_u_afgdt_nabla_u_nabla_u[3][3];
          
        am_nabla_u_afgdt_nabla_u_nabla_u[0][0]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(0,0)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*vderxyaf_(0,0)
                                        +
                                        vderxyaf_(0,1)*vderxyaf_(1,0)
                                        +
                                        vderxyaf_(0,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][1]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(0,1)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*vderxyaf_(0,1)
                                        +
                                        vderxyaf_(0,1)*vderxyaf_(1,1)
                                        +
                                        vderxyaf_(0,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[0][2]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(0,2)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*vderxyaf_(0,2)
                                        +
                                        vderxyaf_(0,1)*vderxyaf_(1,2)
                                        +
                                        vderxyaf_(0,2)*vderxyaf_(2,2));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][0]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(1,0)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(1,0)*vderxyaf_(0,0)
                                        +
                                        vderxyaf_(1,1)*vderxyaf_(1,0)
                                        +
                                        vderxyaf_(1,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][1]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(1,1)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(1,0)*vderxyaf_(0,1)
                                        +
                                        vderxyaf_(1,1)*vderxyaf_(1,1)
                                        +
                                        vderxyaf_(1,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[1][2]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(1,2)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(1,0)*vderxyaf_(0,2)
                                        +
                                        vderxyaf_(1,1)*vderxyaf_(1,2)
                                        +
                                        vderxyaf_(1,2)*vderxyaf_(2,2));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][0]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(2,0)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(2,0)*vderxyaf_(0,0)
                                        +
                                        vderxyaf_(2,1)*vderxyaf_(1,0)
                                        +
                                        vderxyaf_(2,2)*vderxyaf_(2,0));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][1]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(2,1)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(2,0)*vderxyaf_(0,1)
                                        +
                                        vderxyaf_(2,1)*vderxyaf_(1,1)
                                        +
                                        vderxyaf_(2,2)*vderxyaf_(2,1));
        am_nabla_u_afgdt_nabla_u_nabla_u[2][2]=
          fac_afgdt_alphaM_tauM_facMtau*vderxyaf_(2,2)
          +
          fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(2,0)*vderxyaf_(0,2)
                                        +
                                        vderxyaf_(2,1)*vderxyaf_(1,2)
                                        +
                                        vderxyaf_(2,2)*vderxyaf_(2,2));
          
        double nabla_u[3][3];
        nabla_u[0][0]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(0,0);
        nabla_u[0][1]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(0,1);
        nabla_u[0][2]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(0,2);
          
        nabla_u[1][0]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(1,0);
        nabla_u[1][1]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(1,1);
        nabla_u[1][2]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(1,2);
          
        nabla_u[2][0]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(2,0);
        nabla_u[2][1]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(2,1);
        nabla_u[2][2]=fac_afgdt_afgdt_tauM_facMtau*vderxyaf_(2,2);
          
        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
            
          const double u_nabla_ui = velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
            
          double coltemp[3][3];
            
          coltemp[0][0]=am_nabla_u_afgdt_nabla_u_nabla_u[0][0]*funct_(ui)+nabla_u[0][0]*u_nabla_ui;
          coltemp[0][1]=am_nabla_u_afgdt_nabla_u_nabla_u[0][1]*funct_(ui)+nabla_u[0][1]*u_nabla_ui;
          coltemp[0][2]=am_nabla_u_afgdt_nabla_u_nabla_u[0][2]*funct_(ui)+nabla_u[0][2]*u_nabla_ui;
            
          coltemp[1][0]=am_nabla_u_afgdt_nabla_u_nabla_u[1][0]*funct_(ui)+nabla_u[1][0]*u_nabla_ui;
          coltemp[1][1]=am_nabla_u_afgdt_nabla_u_nabla_u[1][1]*funct_(ui)+nabla_u[1][1]*u_nabla_ui;
          coltemp[1][2]=am_nabla_u_afgdt_nabla_u_nabla_u[1][2]*funct_(ui)+nabla_u[1][2]*u_nabla_ui;

          coltemp[2][0]=am_nabla_u_afgdt_nabla_u_nabla_u[2][0]*funct_(ui)+nabla_u[2][0]*u_nabla_ui;
          coltemp[2][1]=am_nabla_u_afgdt_nabla_u_nabla_u[2][1]*funct_(ui)+nabla_u[2][1]*u_nabla_ui;
          coltemp[2][2]=am_nabla_u_afgdt_nabla_u_nabla_u[2][2]*funct_(ui)+nabla_u[2][2]*u_nabla_ui;
            
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
               
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
              
            elemat(fvi  ,fui  ) -= funct_(vi)*coltemp[0][0];
            elemat(fvi  ,fuip ) -= funct_(vi)*coltemp[0][1];
            elemat(fvi  ,fuipp) -= funct_(vi)*coltemp[0][2];
              
            elemat(fvip ,fui  ) -= funct_(vi)*coltemp[1][0];
            elemat(fvip ,fuip ) -= funct_(vi)*coltemp[1][1];
            elemat(fvip ,fuipp) -= funct_(vi)*coltemp[1][2];
              
            elemat(fvipp,fui  ) -= funct_(vi)*coltemp[2][0];
            elemat(fvipp,fuip ) -= funct_(vi)*coltemp[2][1];
            elemat(fvipp,fuipp) -= funct_(vi)*coltemp[2][2];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)
        

        const double fac_afgdt_tauM_facMtau_gdt = fac*alphaF*gamma*dt*tauM*facMtau*gamma*dt;

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fuippp =4*ui+3;

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

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
            elemat(fvi  ,fuippp) -= fac_afgdt_tauM_facMtau_gdt*funct_(vi)*(vderxyaf_(0,0)*derxy_(0,ui)+vderxyaf_(0,1)*derxy_(1,ui)+vderxyaf_(0,2)*derxy_(2,ui));
            elemat(fvip ,fuippp) -= fac_afgdt_tauM_facMtau_gdt*funct_(vi)*(vderxyaf_(1,0)*derxy_(0,ui)+vderxyaf_(1,1)*derxy_(1,ui)+vderxyaf_(1,2)*derxy_(2,ui));
            elemat(fvipp,fuippp) -= fac_afgdt_tauM_facMtau_gdt*funct_(vi)*(vderxyaf_(2,0)*derxy_(0,ui)+vderxyaf_(2,1)*derxy_(1,ui)+vderxyaf_(2,2)*derxy_(2,ui));
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
    
            double coltemp[3][3];
              
            coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(0,1)+derxy2_(4,ui)*vderxyaf_(0,2));
            coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(0,0)+derxy2_(5,ui)*vderxyaf_(0,2));
            coltemp[0][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(0,0)+derxy2_(5,ui)*vderxyaf_(0,1));

            coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(1,1)+derxy2_(4,ui)*vderxyaf_(1,2));
            coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(1,0)+derxy2_(5,ui)*vderxyaf_(1,2));
            coltemp[1][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(1,0)+derxy2_(5,ui)*vderxyaf_(1,1));

            coltemp[2][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(0,ui)+derxy2_(3,ui)*vderxyaf_(2,1)+derxy2_(4,ui)*vderxyaf_(2,2));
            coltemp[2][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(1,ui)+derxy2_(3,ui)*vderxyaf_(2,0)+derxy2_(5,ui)*vderxyaf_(2,2));
            coltemp[2][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*(viscs2_(2,ui)+derxy2_(4,ui)*vderxyaf_(2,0)+derxy2_(5,ui)*vderxyaf_(2,1));

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvi+2;

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
              elemat(fvi  ,fui  ) += coltemp[0][0]*funct_(vi);
              elemat(fvi  ,fuip ) += coltemp[0][1]*funct_(vi);
              elemat(fvi  ,fuipp) += coltemp[0][2]*funct_(vi);

              elemat(fvip ,fui  ) += coltemp[1][0]*funct_(vi);
              elemat(fvip ,fuip ) += coltemp[1][1]*funct_(vi);
              elemat(fvip ,fuipp) += coltemp[1][2]*funct_(vi);

              elemat(fvipp,fui  ) += coltemp[2][0]*funct_(vi);
              elemat(fvipp,fuip ) += coltemp[2][1]*funct_(vi);
              elemat(fvipp,fuipp) += coltemp[2][2]*funct_(vi);

            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
        } // end if higher_order_element
      } // end cross


      if(reynolds == Fluid3::reynolds_stress_stab)
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
        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_z = fac_alphaM_afgdt_tauM_facMtau*svelaf_(2);

        const double fac_afgdt_afgdt_tauM_facMtau =fac*afgdt*afgdt*tauM*facMtau;

        double fac_afgdt_afgdt_tauM_facMtau_svelaf[3];
        fac_afgdt_afgdt_tauM_facMtau_svelaf[0]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[1]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(1);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[2]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
            
          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
            
          double inertia_and_conv[3];

          inertia_and_conv[0]=fac_afgdt_afgdt_tauM_facMtau_svelaf[0]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_afgdt_tauM_facMtau_svelaf[1]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_y*funct_(ui);
          inertia_and_conv[2]=fac_afgdt_afgdt_tauM_facMtau_svelaf[2]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_z*funct_(ui);
              
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

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

            elemat(fvi  ,fui  ) += inertia_and_conv[0]*derxy_(0,vi);
            elemat(fvi  ,fuip ) += inertia_and_conv[0]*derxy_(1,vi);
            elemat(fvi  ,fuipp) += inertia_and_conv[0]*derxy_(2,vi);
                                     
            elemat(fvip ,fui  ) += inertia_and_conv[1]*derxy_(0,vi);
            elemat(fvip ,fuip ) += inertia_and_conv[1]*derxy_(1,vi);
            elemat(fvip ,fuipp) += inertia_and_conv[1]*derxy_(2,vi);
                                     
            elemat(fvipp,fui  ) += inertia_and_conv[2]*derxy_(0,vi);
            elemat(fvipp,fuip ) += inertia_and_conv[2]*derxy_(1,vi);
            elemat(fvipp,fuipp) += inertia_and_conv[2]*derxy_(2,vi);
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
            
          double temp[3];
          temp[0]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi));
          temp[1]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi));
          temp[2]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi));

          double rowtemp[3][3];
            
          rowtemp[0][0]=svelaf_(0)*temp[0];
          rowtemp[0][1]=svelaf_(0)*temp[1];
          rowtemp[0][2]=svelaf_(0)*temp[2];
                               
          rowtemp[1][0]=svelaf_(1)*temp[0];
          rowtemp[1][1]=svelaf_(1)*temp[1];
          rowtemp[1][2]=svelaf_(1)*temp[2];
                               
          rowtemp[2][0]=svelaf_(2)*temp[0];
          rowtemp[2][1]=svelaf_(2)*temp[1];
          rowtemp[2][2]=svelaf_(2)*temp[2];
              
          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
            
            /* 
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / /            \   n+af \         \     |
             |  u     , | | | Dacc o nabla | u      | o nabla | v  |
             |           \ \ \            /        /         /     |
              \                                                   /      

            */

            elemat(fvi  ,fui  ) += funct_(ui)*rowtemp[0][0];
            elemat(fvi  ,fuip ) += funct_(ui)*rowtemp[0][1];
            elemat(fvi  ,fuipp) += funct_(ui)*rowtemp[0][2];
                                                             
            elemat(fvip ,fui  ) += funct_(ui)*rowtemp[1][0];
            elemat(fvip ,fuip ) += funct_(ui)*rowtemp[1][1];
            elemat(fvip ,fuipp) += funct_(ui)*rowtemp[1][2];
                                                             
            elemat(fvipp,fui  ) += funct_(ui)*rowtemp[2][0];
            elemat(fvipp,fuip ) += funct_(ui)*rowtemp[2][1];
            elemat(fvipp,fuipp) += funct_(ui)*rowtemp[2][2];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        const double fac_gdt_afgdt_tauM_facMtau         =fac*gamma*dt*afgdt*tauM*facMtau;
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_x=fac_gdt_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_y=fac_gdt_afgdt_tauM_facMtau*svelaf_(1);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_z=fac_gdt_afgdt_tauM_facMtau*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fuippp =4*ui+3;
            
          double coltemp[3][3];

          coltemp[0][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(1,ui);
          coltemp[0][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(2,ui);
          coltemp[1][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(1,ui);
          coltemp[1][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(2,ui);
          coltemp[2][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(0,ui);
          coltemp[2][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(1,ui);
          coltemp[2][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(2,ui);
  
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
            
            /* 
                 factor: + gamma * dt * alphaF * gamma * dt * tauM *facMtau (rescaled)
             
              /                                \
             |  ~n+af    /                \     |
             |  u     , | nabla Dp o nabla | v  |
             |           \                /     |
              \                                /      

            */

            elemat(fvi  ,fuippp) += coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi)+coltemp[0][2]*derxy_(2,vi);
            elemat(fvip ,fuippp) += coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi)+coltemp[1][2]*derxy_(2,vi);
            elemat(fvipp,fuippp) += coltemp[2][0]*derxy_(0,vi)+coltemp[2][1]*derxy_(1,vi)+coltemp[2][2]*derxy_(2,vi);

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_nu_afgdt_afgdt_tauM_facMtau =fac*visceff*afgdt*afgdt*tauM*facMtau;
          
          double temp[3];

          temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(0);
          temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(1);
          temp[2]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(2);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
              
            double rowtemp[3][3];
              
            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);
            rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);
            rowtemp[1][2]=temp[1]*derxy_(2,vi);
              
            rowtemp[2][0]=temp[2]*derxy_(0,vi);
            rowtemp[2][1]=temp[2]*derxy_(1,vi);
            rowtemp[2][2]=temp[2]*derxy_(2,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const int fui   =4*ui;
              const int fuip  =fui+1;
              const int fuipp =fui+2;
                
              /* 
                   factor: - 2.0 * visc * alphaF * gamma * dt * alphaF * gamma * dt * tauM * facMtauM

                    /                                                 \
                   |  ~n+af    / /             /    \  \         \     |
                   |  u     , | | nabla o eps | Dacc |  | o nabla | v  |
                   |           \ \             \    /  /         /     |
                    \                                                 /      
              */

              elemat(fvi  ,fui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
              elemat(fvi  ,fuip ) -= derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
              elemat(fvi  ,fuipp) -= derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
              elemat(fvip ,fui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
              elemat(fvip ,fuip ) -= derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
              elemat(fvip ,fuipp) -= derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
              elemat(fvipp,fui  ) -= viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
              elemat(fvipp,fuip ) -= derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
              elemat(fvipp,fuipp) -= derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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
    if(inertia == Fluid3::inertia_stab_keep
       ||
       inertia == Fluid3::inertia_stab_keep_complete)
    {
      
      double aux_x =(-svelaf_(0)/tauM-pderxynp_(0)) ;
      double aux_y =(-svelaf_(1)/tauM-pderxynp_(1)) ;
      double aux_z =(-svelaf_(2)/tauM-pderxynp_(2)) ;
      
      if(higher_order_ele)
      {
        const double fact =visceff;
        
        aux_x += fact*viscaf_old_(0);
        aux_y += fact*viscaf_old_(1);
        aux_z += fact*viscaf_old_(2);
      }
      
      const double fac_sacc_plus_resM_not_partially_integrated_x =fac*aux_x ;
      const double fac_sacc_plus_resM_not_partially_integrated_y =fac*aux_y ;
      const double fac_sacc_plus_resM_not_partially_integrated_z =fac*aux_z ;
      
      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;
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
        
        elevec(fui  ) -= fac_sacc_plus_resM_not_partially_integrated_x*funct_(ui) ;
        elevec(fui+1) -= fac_sacc_plus_resM_not_partially_integrated_y*funct_(ui) ;
        elevec(fui+2) -= fac_sacc_plus_resM_not_partially_integrated_z*funct_(ui) ;
      }
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
      
      const double fac_inertia_convection_dead_load_z
        =
        fac*(accintam_(2)+convaf_old_(2)-bodyforceaf_(2));
      
      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;
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
        
        elevec(fui  ) -= funct_(ui)*fac_inertia_convection_dead_load_x;
        elevec(fui+1) -= funct_(ui)*fac_inertia_convection_dead_load_y;
        elevec(fui+2) -= funct_(ui)*fac_inertia_convection_dead_load_z;
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
      const int fui =4*ui;
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
    
      elevec(fui  ) += fac_prenp_*derxy_(0,ui) ;
      elevec(fui+1) += fac_prenp_*derxy_(1,ui) ;
      elevec(fui+2) += fac_prenp_*derxy_(2,ui) ;
    }

    const double visceff_fac=visceff*fac;
      
    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
    {
      const int fui=4*ui;
          
      /* viscous term */
      
      /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
      */
          
      elevec(fui   ) -= visceff_fac*
                        (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                         +
                         derxy_(1,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                         +
                         derxy_(2,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0)));
      elevec(fui+1) -= visceff_fac*
	               (derxy_(0,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                        +
                        derxy_(1,ui)*vderxyaf_(1,1)*2.0
                        +
                        derxy_(2,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1)));
      elevec(fui+2) -= visceff_fac*
	               (derxy_(0,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0))
                        +
                        derxy_(1,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1))
                        +
                        derxy_(2,ui)*vderxyaf_(2,2)*2.0);
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
      
      elevec(ui*4 + 3) -= fac_divunp*funct_(ui);
    } // end loop rows (solution for matrix, test function for vector)

    //---------------------------------------------------------------
    //
    //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //                    PRESSURE STABILISATION
    //
    //---------------------------------------------------------------
    if(pspg == Fluid3::pstab_use_pspg)
    {
      
      const double fac_svelnpx                      = fac*svelnp(0,iquad);
      const double fac_svelnpy                      = fac*svelnp(1,iquad);
      const double fac_svelnpz                      = fac*svelnp(2,iquad);

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        /* factor: -1

                       /                 \
                      |  ~n+1             |
                      |  u    , nabla  q  |
                      |   (i)             |
                       \                 /
        */
        
        elevec(ui*4 + 3) += fac_svelnpx*derxy_(0,ui)+fac_svelnpy*derxy_(1,ui)+fac_svelnpz*derxy_(2,ui);
        
      } // end loop rows (solution for matrix, test function for vector)
    }
    
    //---------------------------------------------------------------
    //
    //         STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
    //
    //---------------------------------------------------------------
    if(supg == Fluid3::convective_stab_supg)
    {
      const double fac_svelaf_x=fac*svelaf_(0);
      const double fac_svelaf_y=fac*svelaf_(1);
      const double fac_svelaf_z=fac*svelaf_(2);

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;
        /*
                  /                             \
                 |  ~n+af    / n+af        \     |
                 |  u     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /

        */

        elevec(fui  ) += fac_svelaf_x*conv_c_af_(ui);
        elevec(fui+1) += fac_svelaf_y*conv_c_af_(ui);
        elevec(fui+2) += fac_svelaf_z*conv_c_af_(ui);

      } // end loop rows (solution for matrix, test function for vector)
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
        const double fac_visc_svelaf_x = vstabfac*fac*visc*svelaf_(0);
        const double fac_visc_svelaf_y = vstabfac*fac*visc*svelaf_(1);
        const double fac_visc_svelaf_z = vstabfac*fac*visc*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fui=4*ui;
          /*
                 /                        \
                |  ~n+af                   |
                |  u      , 2*div eps (v)  |
                |                          |
                 \                        /

          */
          elevec(fui  ) += fac_visc_svelaf_x*viscs2_(0,ui)
	                   +
                           fac_visc_svelaf_y*derxy2_(3,ui)
                           +
                           fac_visc_svelaf_z*derxy2_(4,ui);

          elevec(fui+1) += fac_visc_svelaf_x*derxy2_(3,ui)
                           +
                           fac_visc_svelaf_y*viscs2_(1,ui)
                           +
                           fac_visc_svelaf_z*derxy2_(5,ui);

          elevec(fui+2) += fac_visc_svelaf_x*derxy2_(4,ui)
                           +
                           fac_visc_svelaf_y*derxy2_(5,ui)
                           +
                           fac_visc_svelaf_z*viscs2_(2,ui);

        } // end loop rows (solution for matrix, test function for vector)
      } // endif (a)gls
    }// end if higher order ele

    //---------------------------------------------------------------
    //
    //        TIME-DEPENDENT SUBGRID-SCALE STABILISATION
    //       RESIDUAL BASED VMM STABILISATION --- CROSS STRESS
    //
    //---------------------------------------------------------------
    if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
    {
      const double fac_convsubaf_old_x=fac*convsubaf_old_(0);
      const double fac_convsubaf_old_y=fac*convsubaf_old_(1);
      const double fac_convsubaf_old_z=fac*convsubaf_old_(2);

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;

        /* factor:

                  /                           \
                 |   ~n+af           n+af      |
                 | ( u    o nabla ) u     , v  |
                 |    (i)            (i)       |
                  \                           /
        */
        elevec(fui  ) -= fac_convsubaf_old_x*funct_(ui);
        elevec(fui+1) -= fac_convsubaf_old_y*funct_(ui);
        elevec(fui+2) -= fac_convsubaf_old_z*funct_(ui);
      } // ui
    } // end cross

    //---------------------------------------------------------------
    //
    //       TIME DEPENDENT SUBGRID-SCALE STABILISATION PART
    //     RESIDUAL BASED VMM STABILISATION --- REYNOLDS STRESS
    //
    //---------------------------------------------------------------
    if(reynolds != Fluid3::reynolds_stress_stab_none)
    {
      const double fac_svelaf_x=fac*svelaf_(0);
      const double fac_svelaf_y=fac*svelaf_(1);
      const double fac_svelaf_z=fac*svelaf_(2);

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;

        /* factor:

                  /                             \
                 |  ~n+af      ~n+af             |
                 |  u      , ( u    o nabla ) v  |
                 |                               |
                  \                             /
        */
        elevec(fui  ) += fac_svelaf_x*(svelaf_(0)*derxy_(0,ui)
                                       +
                                       svelaf_(1)*derxy_(1,ui)
                                       +
                                       svelaf_(2)*derxy_(2,ui));
        elevec(fui+1) += fac_svelaf_y*(svelaf_(0)*derxy_(0,ui)
                                       +
                                       svelaf_(1)*derxy_(1,ui)
                                       +
                                       svelaf_(2)*derxy_(2,ui));
        elevec(fui+2) += fac_svelaf_z*(svelaf_(0)*derxy_(0,ui)
                                       +
                                       svelaf_(1)*derxy_(1,ui)
                                       +
                                       svelaf_(2)*derxy_(2,ui));

      } // end ui
    } // end reynolds
  } // end loop iquad
  return;
} //Sysmat_adv_td

/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  conservative, quasistatic version                                  |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Sysmat_cons_qs(
  Fluid3*                                          ele             ,
  std::vector<blitz::Array<double,1> >&            myknots         ,
  LINALG::Matrix<4*iel,4*iel>&                     elemat          ,
  LINALG::Matrix<4*iel,1>&                         elevec          ,
  const LINALG::Matrix<3,iel>&                     edispnp         ,
  const LINALG::Matrix<3,iel>&                     egridvaf        ,
  const LINALG::Matrix<3,iel>&                     evelnp          ,
  const LINALG::Matrix<iel,1>&                     eprenp          ,
  const LINALG::Matrix<3,iel>&                     eaccam          ,
  const LINALG::Matrix<3,iel>&                     evelaf          ,
  const LINALG::Matrix<3,iel>&                     fsevelaf        ,
  const struct _MATERIAL*                          material        ,
  const double                                     alphaM          ,
  const double                                     alphaF          ,
  const double                                     gamma           ,
  const double                                     dt              ,
  const double                                     time            ,
  const enum Fluid3::LinearisationAction           newton          ,
  const bool                                       higher_order_ele,
  const enum Fluid3::StabilisationAction           fssgv           ,
  const enum Fluid3::StabilisationAction           pspg            ,
  const enum Fluid3::StabilisationAction           supg            ,
  const enum Fluid3::StabilisationAction           vstab           ,
  const enum Fluid3::StabilisationAction           cstab           ,
  const enum Fluid3::StabilisationAction           cross           ,
  const enum Fluid3::StabilisationAction           reynolds        ,
  const enum Fluid3::TauType                       whichtau        ,
  const enum Fluid3::TurbModelAction               turb_mod_action ,
  double&                                          Cs              ,
  double&                                          Cs_delta_sq     ,
  double&                                          visceff         ,
  const double                                     l_tau           ,
  const bool                                       compute_elemat
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
  if (vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_gls_only_rhs)
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

  SetElementData(ele            ,
                 edispnp        ,
                 evelaf         ,
                 fsevelaf       ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 fssgv          ,
                 turb_mod_action,
                 l_tau          ,
                 Cs             ,
                 Cs_delta_sq    ,
                 visceff        );

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second 
    // derivatives in a gausspoint and integration weight including 
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);
    
    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
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
    if(cross == Fluid3::cross_stress_stab)
    {
      for(int nn=0;nn<iel;++nn)
      {
        conv_resM_(nn)=resM_(0)*derxy_(0,nn);
        
        for(int rr=1;rr<3;++rr)
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
    if(ele->is_ale_)
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=u_G_af_(0)*derxy_(0,nn);
        
	for(int rr=1;rr<3;++rr)
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
    for(int rr=0;rr<3;++rr)
    {
      convu_G_af_old_(rr)=u_G_af_(0)*vderxyaf_(rr,0);
      for(int mm=1;mm<3;++mm)
      {
        convu_G_af_old_(rr)+=u_G_af_(mm)*vderxyaf_(rr,mm);
      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,Fluid3::subscales_quasistatic,gamma,dt,hk,mk,visceff);

    // stabilisation parameters
    const double tauM   = tau_(0);
    const double tauMp  = tau_(1);

    if(cstab == Fluid3::continuity_stab_none)
    {
      tau_(2)=0.0;
    }
    const double tauC    = tau_(2);

    double supg_active_tauM;
    if(supg == Fluid3::convective_stab_supg)
    {
      supg_active_tauM=tauM;
    }
    else
    {
      supg_active_tauM=0.0;
    }

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;

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

      if(reynolds == Fluid3::reynolds_stress_stab)
      {

        /* half of the reynolds linearisation is done by modifiing 
           the supg testfunction, see above */

        for(int nn=0;nn<iel;++nn)
        {
          conv_c_plus_svel_af_(nn)-=tauM*tauM*resM_(0)*derxy_(0,nn);
            
          for(int rr=1;rr<3;++rr)
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
        const double fac_alphaM_tauM_tauM_resM_z=fac_alphaM_tauM_tauM*resM_(2);

        const double fac_afgdt_tauM_tauM=fac*afgdt*tauM*tauM;

        double fac_afgdt_tauM_tauM_resM[3];
        fac_afgdt_tauM_tauM_resM[0]=fac_afgdt_tauM_tauM*resM_(0);
        fac_afgdt_tauM_tauM_resM[1]=fac_afgdt_tauM_tauM*resM_(1);
        fac_afgdt_tauM_tauM_resM[2]=fac_afgdt_tauM_tauM*resM_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
            
          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
            
          double inertia_and_conv[3];

          inertia_and_conv[0]=fac_afgdt_tauM_tauM_resM[0]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_tauM_tauM_resM[1]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_y*funct_(ui);
          inertia_and_conv[2]=fac_afgdt_tauM_tauM_resM[2]*u_o_nabla_ui+fac_alphaM_tauM_tauM_resM_z*funct_(ui);
              
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

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

            elemat(fvi  ,fui  ) -= inertia_and_conv[0]*derxy_(0,vi);
            elemat(fvi  ,fuip ) -= inertia_and_conv[0]*derxy_(1,vi);
            elemat(fvi  ,fuipp) -= inertia_and_conv[0]*derxy_(2,vi);
                                     
            elemat(fvip ,fui  ) -= inertia_and_conv[1]*derxy_(0,vi);
            elemat(fvip ,fuip ) -= inertia_and_conv[1]*derxy_(1,vi);
            elemat(fvip ,fuipp) -= inertia_and_conv[1]*derxy_(2,vi);
                                     
            elemat(fvipp,fui  ) -= inertia_and_conv[2]*derxy_(0,vi);
            elemat(fvipp,fuip ) -= inertia_and_conv[2]*derxy_(1,vi);
            elemat(fvipp,fuipp) -= inertia_and_conv[2]*derxy_(2,vi);
          } // vi
        } // ui

              
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
            
          double temp[3];
          temp[0]=fac_afgdt_tauM_tauM*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi));
          temp[1]=fac_afgdt_tauM_tauM*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi));
          temp[2]=fac_afgdt_tauM_tauM*(vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi));

          double rowtemp[3][3];
            
          rowtemp[0][0]=resM_(0)*temp[0];
          rowtemp[0][1]=resM_(0)*temp[1];
          rowtemp[0][2]=resM_(0)*temp[2];
                                   
          rowtemp[1][0]=resM_(1)*temp[0];
          rowtemp[1][1]=resM_(1)*temp[1];
          rowtemp[1][2]=resM_(1)*temp[2];
                                   
          rowtemp[2][0]=resM_(2)*temp[0];
          rowtemp[2][1]=resM_(2)*temp[1];
          rowtemp[2][2]=resM_(2)*temp[2];
              
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
            
            /* 
                 factor: -alphaF * gamma * dt * tauM * tauM

              /                                                  \
             |          / / /            \   n+af \         \     |
             |  resM , | | | Dacc o nabla | u      | o nabla | v  |
             |          \ \ \            /        /         /     |
              \                                                  /      

            */

            elemat(fvi  ,fui  ) -= funct_(ui)*rowtemp[0][0];
            elemat(fvi  ,fuip ) -= funct_(ui)*rowtemp[0][1];
            elemat(fvi  ,fuipp) -= funct_(ui)*rowtemp[0][2];
                                                             
            elemat(fvip ,fui  ) -= funct_(ui)*rowtemp[1][0];
            elemat(fvip ,fuip ) -= funct_(ui)*rowtemp[1][1];
            elemat(fvip ,fuipp) -= funct_(ui)*rowtemp[1][2];
                                                             
            elemat(fvipp,fui  ) -= funct_(ui)*rowtemp[2][0];
            elemat(fvipp,fuip ) -= funct_(ui)*rowtemp[2][1];
            elemat(fvipp,fuipp) -= funct_(ui)*rowtemp[2][2];
          } // ui
        } // vi


        const double fac_gdt_tauM_tauM       =fac*gamma*dt*tauM*tauM;
        const double fac_gdt_tauM_tauM_resM_x=fac_gdt_tauM_tauM*resM_(0);
        const double fac_gdt_tauM_tauM_resM_y=fac_gdt_tauM_tauM*resM_(1);
        const double fac_gdt_tauM_tauM_resM_z=fac_gdt_tauM_tauM*resM_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fuippp =4*ui+3;
            
          double coltemp[3][3];

          coltemp[0][0]=fac_gdt_tauM_tauM_resM_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_tauM_tauM_resM_x*derxy_(1,ui);
          coltemp[0][2]=fac_gdt_tauM_tauM_resM_x*derxy_(2,ui);
          coltemp[1][0]=fac_gdt_tauM_tauM_resM_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_tauM_tauM_resM_y*derxy_(1,ui);
          coltemp[1][2]=fac_gdt_tauM_tauM_resM_y*derxy_(2,ui);
          coltemp[2][0]=fac_gdt_tauM_tauM_resM_z*derxy_(0,ui);
          coltemp[2][1]=fac_gdt_tauM_tauM_resM_z*derxy_(1,ui);
          coltemp[2][2]=fac_gdt_tauM_tauM_resM_z*derxy_(2,ui);
  
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
            
            /* 
                 factor: - gamma * dt * tauM * tauM (rescaled)
             
              /                               \
             |          /                \     |
             |  resM , | nabla Dp o nabla | v  |
             |          \                /     |
              \                               /      

            */

            elemat(fvi  ,fuippp) -= coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi)+coltemp[0][2]*derxy_(2,vi);
            elemat(fvip ,fuippp) -= coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi)+coltemp[1][2]*derxy_(2,vi);
            elemat(fvipp,fuippp) -= coltemp[2][0]*derxy_(0,vi)+coltemp[2][1]*derxy_(1,vi)+coltemp[2][2]*derxy_(2,vi);

          } // vi
        } // ui


        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_nu_afgdt_tauM_tauM=fac*visceff*afgdt*tauM*tauM;
          
          double temp[3];

          temp[0]=fac_nu_afgdt_tauM_tauM*resM_(0);
          temp[1]=fac_nu_afgdt_tauM_tauM*resM_(1);
          temp[2]=fac_nu_afgdt_tauM_tauM*resM_(2);


          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
              
            double rowtemp[3][3];
              
            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);
            rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);
            rowtemp[1][2]=temp[1]*derxy_(2,vi);
              
            rowtemp[2][0]=temp[2]*derxy_(0,vi);
            rowtemp[2][1]=temp[2]*derxy_(1,vi);
            rowtemp[2][2]=temp[2]*derxy_(2,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int fui   =4*ui;
              const int fuip  =fui+1;
              const int fuipp =fui+2;
                
              /* 
                   factor: + 2.0 * visc * alphaF * gamma * dt * tauM * tauM 

                    /                                                \
                   |          / /             /    \  \         \     |
                   |  resM , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /      
              */

              elemat(fvi  ,fui  ) += viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
              elemat(fvi  ,fuip ) += derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
              elemat(fvi  ,fuipp) += derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
              elemat(fvip ,fui  ) += viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
              elemat(fvip ,fuip ) += derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
              elemat(fvip ,fuipp) += derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
              elemat(fvipp,fui  ) += viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
              elemat(fvipp,fuip ) += derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
              elemat(fvipp,fuipp) += derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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
      const double fac_afgdt_velintaf_z=fac_afgdt*velintaf_(2);

      // supg and cstab conservative
      const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;
            
      for (int ui=0; ui<iel; ++ui) // loop columns 
      {
        const int fui  =4*ui;
        const int fuip =fui+1;
        const int fuipp=fui+2;
            
        /* GALERKIN inertia term (intermediate) + convection, mesh velocity (intermediate) */
        const double inertia_and_gridconv_ui
          = fac_alphaM*funct_(ui)-fac_afgdt*conv_u_G_af_(ui);
        
        /* SUPG stabilisation --- inertia and convection */
        const double inertia_and_conv
          = fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);
          
        // convection GALERKIN and diagonal parts of viscous term (intermediate)
        const double convection_and_viscous_x=fac_visceff_afgdt*derxy_(0,ui)-fac_afgdt_velintaf_x*funct_(ui);
        const double convection_and_viscous_y=fac_visceff_afgdt*derxy_(1,ui)-fac_afgdt_velintaf_y*funct_(ui);
        const double convection_and_viscous_z=fac_visceff_afgdt*derxy_(2,ui)-fac_afgdt_velintaf_z*funct_(ui);
        
        // viscous GALERKIN term
        const double viscous_x=fac_visceff_afgdt*derxy_(0,ui);
        const double viscous_y=fac_visceff_afgdt*derxy_(1,ui);
        const double viscous_z=fac_visceff_afgdt*derxy_(2,ui);
        
        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);
        const double fac_gamma_dt_tauC_derxy_z_ui = fac_gamma_dt_tauC*derxy_(2,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
              
          const double sum =
            inertia_and_gridconv_ui*funct_(vi)
            +
            inertia_and_conv*conv_c_plus_svel_af_(vi)
            +
            convection_and_viscous_x*derxy_(0,vi)
            +
            convection_and_viscous_y*derxy_(1,vi)
            +
            convection_and_viscous_z*derxy_(2,vi);
          
          /* adding GALERKIN convection, convective linearisation (intermediate), viscous and cstab */
          
          elemat(fvi  ,fui  ) += sum+(fac_gamma_dt_tauC_derxy_x_ui             + viscous_x)*derxy_(0,vi);
          elemat(fvi  ,fuip ) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi)+(viscous_x)*derxy_(1,vi);
          elemat(fvi  ,fuipp) +=      fac_gamma_dt_tauC_derxy_z_ui*derxy_(0,vi)+(viscous_x)*derxy_(2,vi);
          elemat(fvip ,fui  ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi)+(viscous_y)*derxy_(0,vi);
          elemat(fvip ,fuip ) += sum+(fac_gamma_dt_tauC_derxy_y_ui             + viscous_y)*derxy_(1,vi);
          elemat(fvip ,fuipp) +=      fac_gamma_dt_tauC_derxy_z_ui*derxy_(1,vi)+(viscous_y)*derxy_(2,vi);
          elemat(fvipp,fui  ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(2,vi)+(viscous_z)*derxy_(0,vi);
          elemat(fvipp,fuip ) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(2,vi)+(viscous_z)*derxy_(1,vi);
          elemat(fvipp,fuipp) += sum+(fac_gamma_dt_tauC_derxy_z_ui             + viscous_z)*derxy_(2,vi);
        } // vi
      } // ui

      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int fuippp=4*ui+3;
	
        const double fac_gamma_dt_derxy_0_ui = fac_gamma_dt*derxy_(0,ui);
        const double fac_gamma_dt_derxy_1_ui = fac_gamma_dt*derxy_(1,ui);
        const double fac_gamma_dt_derxy_2_ui = fac_gamma_dt*derxy_(2,ui);
	  
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          int fvi =vi*4;
	  
          /* SUPG stabilisation --- pressure    */
          /* factor: +tauM, rescaled by gamma*dt*/
   
          elemat(fvi++,fuippp) += fac_gamma_dt_derxy_0_ui*conv_c_plus_svel_af_(vi);
          elemat(fvi++,fuippp) += fac_gamma_dt_derxy_1_ui*conv_c_plus_svel_af_(vi);
          elemat(fvi  ,fuippp) += fac_gamma_dt_derxy_2_ui*conv_c_plus_svel_af_(vi);
          
        } // vi
      } // ui

      if (higher_order_ele && newton!=Fluid3::minimal)
      {
        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui  =ui*4 ;
          const int fuip =fui+1;
          const int fuipp=fui+2;
	  
          /* SUPG stabilisation --- diffusion   */
          /* factor: -nu*alphaF*gamma*dt*tauM   */

          const double fac_visceff_afgdt_viscs2_0_ui=fac_visceff_afgdt*viscs2_(0,ui);
          const double fac_visceff_afgdt_viscs2_1_ui=fac_visceff_afgdt*viscs2_(1,ui);
          const double fac_visceff_afgdt_viscs2_2_ui=fac_visceff_afgdt*viscs2_(2,ui);
          const double fac_visceff_afgdt_derxy2_3_ui=fac_visceff_afgdt*derxy2_(3,ui);
          const double fac_visceff_afgdt_derxy2_4_ui=fac_visceff_afgdt*derxy2_(4,ui);
          const double fac_visceff_afgdt_derxy2_5_ui=fac_visceff_afgdt*derxy2_(5,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            int fvi  =vi*4 ;
            int fvip =fvi+1;
            int fvipp=fvi+2;
	      
            elemat(fvi  ,fui  ) -= fac_visceff_afgdt_viscs2_0_ui*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuip ) -= fac_visceff_afgdt_derxy2_3_ui*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuipp) -= fac_visceff_afgdt_derxy2_4_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fui  ) -= fac_visceff_afgdt_derxy2_3_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuip ) -= fac_visceff_afgdt_viscs2_1_ui*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuipp) -= fac_visceff_afgdt_derxy2_5_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fui  ) -= fac_visceff_afgdt_derxy2_4_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuip ) -= fac_visceff_afgdt_derxy2_5_ui*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuipp) -= fac_visceff_afgdt_viscs2_2_ui*conv_c_plus_svel_af_(vi);
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
          
      if (newton==Fluid3::Newton)
      {
        double temp[3][3];
        double testlin[3];
        
        /* for linearisation of testfunction (SUPG) and reactive GALERKIN part */
        testlin[0]=supg_active_tauM*resM_(0)-velintaf_(0);
        testlin[1]=supg_active_tauM*resM_(1)-velintaf_(1);
        testlin[2]=supg_active_tauM*resM_(2)-velintaf_(2);
            
        // loop rows
        for (int vi=0; vi<iel; ++vi)
        {
          int fvi   =vi*4;
          int fvip  =fvi+1;
          int fvipp =fvi+2;
              
          /*  add linearised convective term in residual (SUPG), reactive 
              GALERKIN part and linearisation of testfunction (SUPG) */
          temp[0][0]=fac_afgdt*(testlin[0]*derxy_(0,vi)+vderxyaf_(0,0)*conv_c_plus_svel_af_(vi));
          temp[0][1]=fac_afgdt*(testlin[0]*derxy_(1,vi)+vderxyaf_(0,1)*conv_c_plus_svel_af_(vi));
          temp[0][2]=fac_afgdt*(testlin[0]*derxy_(2,vi)+vderxyaf_(0,2)*conv_c_plus_svel_af_(vi));
          temp[1][0]=fac_afgdt*(testlin[1]*derxy_(0,vi)+vderxyaf_(1,0)*conv_c_plus_svel_af_(vi));
          temp[1][1]=fac_afgdt*(testlin[1]*derxy_(1,vi)+vderxyaf_(1,1)*conv_c_plus_svel_af_(vi));
          temp[1][2]=fac_afgdt*(testlin[1]*derxy_(2,vi)+vderxyaf_(1,2)*conv_c_plus_svel_af_(vi));
          temp[2][0]=fac_afgdt*(testlin[2]*derxy_(0,vi)+vderxyaf_(2,0)*conv_c_plus_svel_af_(vi));
          temp[2][1]=fac_afgdt*(testlin[2]*derxy_(1,vi)+vderxyaf_(2,1)*conv_c_plus_svel_af_(vi));
          temp[2][2]=fac_afgdt*(testlin[2]*derxy_(2,vi)+vderxyaf_(2,2)*conv_c_plus_svel_af_(vi));
              
          // loop columns
          for (int ui=0; ui<iel; ++ui)
          {
            int fui=4*ui;
            
            elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][0]*funct_(ui);
            elemat(fvipp,fui++) += temp[2][0]*funct_(ui);
            elemat(fvi  ,fui  ) += temp[0][1]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][1]*funct_(ui);
            elemat(fvipp,fui++) += temp[2][1]*funct_(ui);
            elemat(fvi  ,fui  ) += temp[0][2]*funct_(ui);
            elemat(fvip ,fui  ) += temp[1][2]*funct_(ui);
            elemat(fvipp,fui  ) += temp[2][2]*funct_(ui);
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
	const int fvi    =4*vi;
	const int fvip   =fvi+1;
	const int fvipp  =fvi+2;
	  
	const double fac_gamma_dt_derxy_0_vi = fac_gamma_dt*derxy_(0,vi);
	const double fac_gamma_dt_derxy_1_vi = fac_gamma_dt*derxy_(1,vi);
	const double fac_gamma_dt_derxy_2_vi = fac_gamma_dt*derxy_(2,vi);
	  
	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int fuippp=4*ui+3;
	    
	  /* GALERKIN pressure (implicit, rescaled to keep symmetry) */
	  /*  factor: -1, rescaled by gamma*dt */

	  elemat(fvi  ,fuippp) -= fac_gamma_dt_derxy_0_vi*funct_(ui);
	  elemat(fvip ,fuippp) -= fac_gamma_dt_derxy_1_vi*funct_(ui);
	  elemat(fvipp,fuippp) -= fac_gamma_dt_derxy_2_vi*funct_(ui);

	  /* continuity equation (implicit, transposed of above equation) */
	  /*  factor: +gamma*dt */

	  elemat(fuippp,fvi  ) += fac_gamma_dt_derxy_0_vi*funct_(ui);
	  elemat(fuippp,fvip ) += fac_gamma_dt_derxy_1_vi*funct_(ui);
	  elemat(fuippp,fvipp) += fac_gamma_dt_derxy_2_vi*funct_(ui);
	} // ui
      } // vi
      
      //---------------------------------------------------------------
      //
      //             PSPG PART, QUASISTATIC FORMULATION
      //
      //---------------------------------------------------------------
      if(pspg == Fluid3::pstab_use_pspg)
      {
	const double fac_tauMp                   = fac*tauMp;
	const double fac_alphaM_tauMp            = fac_tauMp*alphaM;
	const double fac_gamma_dt_tauMp          = fac_tauMp*gamma*dt;
	const double fac_afgdt_tauMp             = fac_tauMp*afgdt;
	
	if (higher_order_ele && newton!=Fluid3::minimal)
	{
	  const double fac_visceff_afgdt_tauMp
	    =
	    fac*visceff*afgdt*tauMp;
	  
	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int fui  =ui*4;
	    const int fuip =fui+1;
	    const int fuipp=fui+2;
	      
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
	    const double pspg_diffusion_inertia_convect_2_ui
	      =
	      fac_visceff_afgdt_tauMp*viscs2_(2,ui)-fac_tauMp_inertia_and_conv;

	    const double fac_visceff_afgdt_tauMp_derxy2_3_ui=fac_visceff_afgdt_tauMp*derxy2_(3,ui);
	    const double fac_visceff_afgdt_tauMp_derxy2_4_ui=fac_visceff_afgdt_tauMp*derxy2_(4,ui);
	    const double fac_visceff_afgdt_tauMp_derxy2_5_ui=fac_visceff_afgdt_tauMp*derxy2_(5,ui);
                           
	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int fvippp=vi*4+3;

	      elemat(fvippp,fui  ) -= 
		pspg_diffusion_inertia_convect_0_ui*derxy_(0,vi)
		+				                  
		fac_visceff_afgdt_tauMp_derxy2_3_ui*derxy_(1,vi)
		+				                  
		fac_visceff_afgdt_tauMp_derxy2_4_ui*derxy_(2,vi);
	      elemat(fvippp,fuip ) -= 	                  
		fac_visceff_afgdt_tauMp_derxy2_3_ui*derxy_(0,vi)
		+				                  
		pspg_diffusion_inertia_convect_1_ui*derxy_(1,vi)
		+				                  
		fac_visceff_afgdt_tauMp_derxy2_5_ui*derxy_(2,vi);
	      elemat(fvippp,fuipp) -=		                  
		fac_visceff_afgdt_tauMp_derxy2_4_ui*derxy_(0,vi)
		+				                  
		fac_visceff_afgdt_tauMp_derxy2_5_ui*derxy_(1,vi)
		+				                  
		pspg_diffusion_inertia_convect_2_ui*derxy_(2,vi);
	    } // vi
	  } // ui
	} // this is a higher order element and linearisation is not minimal
	else
	{ // either this ain't a higher order element or a
	  // linearisation of the viscous term is not necessary
	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int fui  =ui*4 ;
	    const int fuip =fui+1;
	    const int fuipp=fui+2;
	    
	    const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));
	    
	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int fvippp=vi*4+3;
	      
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
	      
	      elemat(fvippp,fui  ) += fac_tauMp_inertia_and_conv*derxy_(0,vi) ;
	      elemat(fvippp,fuip ) += fac_tauMp_inertia_and_conv*derxy_(1,vi) ;
	      elemat(fvippp,fuipp) += fac_tauMp_inertia_and_conv*derxy_(2,vi) ;
	    } // vi
	  } // ui
	} // no linearisation of viscous part of residual is 
	  // performed for pspg stabilisation cause either this 
  	  // ain't a higher order element or a linearisation of 
	  // the viscous term is not necessary

	if (newton==Fluid3::Newton)
	{
	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    int vidx = vi*4 + 3;
	    double v1 = derxy_(0,vi)*vderxyaf_(0,0) + derxy_(1,vi)*vderxyaf_(1,0) + derxy_(2,vi)*vderxyaf_(2,0);
	    double v2 = derxy_(0,vi)*vderxyaf_(0,1) + derxy_(1,vi)*vderxyaf_(1,1) + derxy_(2,vi)*vderxyaf_(2,1);
	    double v3 = derxy_(0,vi)*vderxyaf_(0,2) + derxy_(1,vi)*vderxyaf_(1,2) + derxy_(2,vi)*vderxyaf_(2,2);
	    for (int ui=0; ui<iel; ++ui) // loop columns
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
	    } // ui
	  } // vi
	} // end newton
	
	for (int ui=0; ui<iel; ++ui) // loop columns
	{
	  const int fuippp=ui*4+3;

	  const double fac_gamma_dt_tauMp_derxy_0_ui=fac_gamma_dt_tauMp*derxy_(0,ui);
	  const double fac_gamma_dt_tauMp_derxy_1_ui=fac_gamma_dt_tauMp*derxy_(1,ui);
	  const double fac_gamma_dt_tauMp_derxy_2_ui=fac_gamma_dt_tauMp*derxy_(2,ui);
	  
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

	    elemat(vi*4+3,fuippp) += 
	      fac_gamma_dt_tauMp_derxy_0_ui*derxy_(0,vi)
	      +
	      fac_gamma_dt_tauMp_derxy_1_ui*derxy_(1,vi)
	      +
	      fac_gamma_dt_tauMp_derxy_2_ui*derxy_(2,vi);
	    
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
	if((vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_usfem)&&higher_order_ele)
	{
	  const double fac_visc_tauMp_gamma_dt      = vstabfac*fac*visc*tauMp*gamma*dt;
	  const double fac_visc_afgdt_tauMp         = vstabfac*fac*visc*afgdt*tauMp;
	  const double fac_visc_alphaM_tauMp        = vstabfac*fac*visc*alphaM*tauMp;
	  const double fac_visceff_visc_afgdt_tauMp = vstabfac*fac*visceff*visc*afgdt*tauMp;

	  for (int ui=0; ui<iel; ++ui) // loop columns
	  {
	    const int fui    = ui*4;
	    const int fuip   = fui+1;
	    const int fuipp  = fui+2;
	    const int fuippp = fui+3;

	    const double acc_conv=(fac_visc_alphaM_tauMp*funct_(ui)
				   +
				   fac_visc_afgdt_tauMp*conv_c_af_(ui));

	    for (int vi=0; vi<iel; ++vi)  // loop rows
	    {
	      const int fvi   = vi*4;
	      const int fvip  = fvi+1;
	      const int fvipp = fvi+2;

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

	      elemat(fvi  ,fui  ) += acc_conv*viscs2_(0,vi);
	      elemat(fvi  ,fuip ) += acc_conv*derxy2_(3,vi);
	      elemat(fvi  ,fuipp) += acc_conv*derxy2_(4,vi);
	      elemat(fvip ,fui  ) += acc_conv*derxy2_(3,vi);
	      elemat(fvip ,fuip ) += acc_conv*viscs2_(1,vi);
	      elemat(fvip ,fuipp) += acc_conv*derxy2_(5,vi);
	      elemat(fvipp,fui  ) += acc_conv*derxy2_(4,vi);
	      elemat(fvipp,fuip ) += acc_conv*derxy2_(5,vi);
	      elemat(fvipp,fuipp) += acc_conv*viscs2_(2,vi);

	      /* viscous stabilisation --- diffusion  */

	      /* factor: -(+)nu*nu*alphaF*gamma*dt*tauMp

                    /                                       \
                   |                 /    \                  |
                   |  2*nabla o eps | Dacc | , 2*div eps (v) |
                   |                 \    /                  |
                    \                                       /
	      */
	      elemat(fvi  ,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*viscs2_(0,vi)
				      +
				      derxy2_(3,ui)*derxy2_(3,vi)
				      +
				      derxy2_(4,ui)*derxy2_(4,vi)) ;
	      elemat(fvi  ,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(3,ui)
				      +
				      derxy2_(3,vi)*viscs2_(1,ui)
				      +
				      derxy2_(4,vi)*derxy2_(5,ui)) ;
	      elemat(fvi  ,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,vi)*derxy2_(4,ui)
				      +
				      derxy2_(3,vi)*derxy2_(5,ui)
				      +
				      derxy2_(4,vi)*viscs2_(2,ui)) ;
	      elemat(fvip ,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(3,vi)
				      +
				      derxy2_(3,ui)*viscs2_(1,vi)
				      +
				      derxy2_(4,ui)*derxy2_(5,vi)) ;
	      elemat(fvip ,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,ui)*derxy2_(3,vi)
				      +
				      viscs2_(1,ui)*viscs2_(1,vi)
				      +
				      derxy2_(5,ui)*derxy2_(5,vi)) ;
	      elemat(fvip ,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,vi)*derxy2_(4,ui)
				      +
				      viscs2_(1,vi)*derxy2_(5,ui)
				      +
				      derxy2_(5,vi)*viscs2_(2,ui)) ;
	      elemat(fvipp,fui  ) -= fac_visceff_visc_afgdt_tauMp*
                                     (viscs2_(0,ui)*derxy2_(4,vi)
				      +
				      derxy2_(3,ui)*derxy2_(5,vi)
				      +
				      derxy2_(4,ui)*viscs2_(2,vi)) ;
	      elemat(fvipp,fuip ) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(3,ui)*derxy2_(4,vi)
				      +
				      viscs2_(1,ui)*derxy2_(5,vi)
				      +
				      derxy2_(5,ui)*viscs2_(2,vi)) ;
	      elemat(fvipp,fuipp) -= fac_visceff_visc_afgdt_tauMp*
                                     (derxy2_(4,ui)*derxy2_(4,vi)
				      +
				      derxy2_(5,ui)*derxy2_(5,vi)
				      +
				      viscs2_(2,ui)*viscs2_(2,vi)) ;

	      /* viscous stabilisation --- pressure   */

	      /* factor: +(-)tauMp*nu, rescaled by gamma*dt

                    /                          \
                   |                            |
                   |  nabla Dp , 2*div eps (v)  |
                   |                            |
                    \                          /
	      */
	      elemat(fvi  ,fuippp) += fac_visc_tauMp_gamma_dt*
                                      (derxy_(0,ui)*viscs2_(0,vi)
				       +
				       derxy_(1,ui)*derxy2_(3,vi)
				       +
				       derxy_(2,ui)*derxy2_(4,vi)) ;
	      elemat(fvip ,fuippp) += fac_visc_tauMp_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(3,vi)
				       +
				       derxy_(1,ui)*viscs2_(1,vi)
				       +
				       derxy_(2,ui)*derxy2_(5,vi)) ;
	      elemat(fvipp,fuippp) += fac_visc_tauMp_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(4,vi)
				       +
				       derxy_(1,ui)*derxy2_(5,vi)
				       +
				       derxy_(2,ui)*viscs2_(2,vi)) ;

	    } // vi
	  } // ui

	  if (newton==Fluid3::Newton)
	  {
	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int fui    = ui*4;
	      const int fuip   = fui+1;
	      const int fuipp  = fui+2;
	      
	      const double fac_visc_afgdt_tauMp_funct_ui=fac_visc_afgdt_tauMp*funct_(ui);
	      
	      for (int vi=0; vi<iel; ++vi)  // loop rows
	      {
		const int fvi   = vi*4;
		const int fvip  = fvi+1;
		const int fvipp = fvi+2;

		/* viscous stabilisation --- convection */

		/*  factor: +(-)nu*alphaF*gamma*dt*tauMp

                     /                                       \
                    |   /            \   n+af                 |
                    |  | Dacc o nabla | u     , 2*div eps (v) |
                    |   \            /                        |
                     \                                       /
		     
		*/
		elemat(fvi  ,fui  ) += fac_visc_afgdt_tauMp_funct_ui*        
                                       (viscs2_(0,vi)*vderxyaf_(0,0)
					+
					derxy2_(3,vi)*vderxyaf_(1,0)
					+
					derxy2_(4,vi)*vderxyaf_(2,0));
		elemat(fvi  ,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,1)
					+
					derxy2_(3,vi)*vderxyaf_(1,1)
					+
					derxy2_(4,vi)*vderxyaf_(2,1));
		elemat(fvi  ,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (viscs2_(0,vi)*vderxyaf_(0,2)
					+
					derxy2_(3,vi)*vderxyaf_(1,2)
					+
					derxy2_(4,vi)*vderxyaf_(2,2));
		elemat(fvip ,fui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,0)
					+
					viscs2_(1,vi)*vderxyaf_(1,0)
					+
					derxy2_(5,vi)*vderxyaf_(2,0));
		elemat(fvip ,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,1)
					+
					viscs2_(1,vi)*vderxyaf_(1,1)
					+
					derxy2_(5,vi)*vderxyaf_(2,1));
		elemat(fvip ,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(3,vi)*vderxyaf_(0,2)
					+
					viscs2_(1,vi)*vderxyaf_(1,2)
					+
					derxy2_(5,vi)*vderxyaf_(2,2));
		elemat(fvipp,fui  ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,0)
					+
					derxy2_(5,vi)*vderxyaf_(1,0)
					+
					viscs2_(2,vi)*vderxyaf_(2,0));
		elemat(fvipp,fuip ) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,1)
					+
					derxy2_(5,vi)*vderxyaf_(1,1)
					+
					viscs2_(2,vi)*vderxyaf_(2,1));
		elemat(fvipp,fuipp) += fac_visc_afgdt_tauMp_funct_ui*
                                       (derxy2_(4,vi)*vderxyaf_(0,2)
					+
					derxy2_(5,vi)*vderxyaf_(1,2)
					+
					viscs2_(2,vi)*vderxyaf_(2,2));
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
      if(cross == Fluid3::cross_stress_stab)
      {
        const double fac_tauM         =fac*tauM;
        const double fac_tauM_alphaM  =fac_tauM*alphaM;
        const double fac_tauM_afgdt   =fac_tauM*afgdt;
        const double fac_tauM_gdt     =fac_tauM*gamma*dt;

        double fac_tauM_alphaM_velintaf[3];
        fac_tauM_alphaM_velintaf[0]=fac_tauM_alphaM*velintaf_(0);
        fac_tauM_alphaM_velintaf[1]=fac_tauM_alphaM*velintaf_(1);
        fac_tauM_alphaM_velintaf[2]=fac_tauM_alphaM*velintaf_(2);

        double fac_tauM_afgdt_velintaf[3];
        fac_tauM_afgdt_velintaf[0]=fac_tauM_afgdt*velintaf_(0);
        fac_tauM_afgdt_velintaf[1]=fac_tauM_afgdt*velintaf_(1);
        fac_tauM_afgdt_velintaf[2]=fac_tauM_afgdt*velintaf_(2);

        for (int vi=0; vi<iel; ++vi)  // loop rows 
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;

          
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
         
          double aux[3];

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
          aux[0]=vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi);
          aux[1]=vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi);
          aux[2]=vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi);

          double temp[3][3];
                          
          /* factor: tauM*alpha_M
	       
                  /                         \ 
                 |    n+af                   |
                 |   u     , Dacc o nabla v  |
                 |                           |
                  \                         / 
          */
          temp[0][0]=fac_tauM_alphaM_velintaf[0]*derxy_(0,vi)+fac_tauM_afgdt_velintaf[0]*aux[0]+fac_tauM_afgdt_conv_resM_vi;
          temp[0][1]=fac_tauM_alphaM_velintaf[0]*derxy_(1,vi)+fac_tauM_afgdt_velintaf[0]*aux[1];
          temp[0][2]=fac_tauM_alphaM_velintaf[0]*derxy_(2,vi)+fac_tauM_afgdt_velintaf[0]*aux[2];
                                                                                     
          temp[1][0]=fac_tauM_alphaM_velintaf[1]*derxy_(0,vi)+fac_tauM_afgdt_velintaf[1]*aux[0];
          temp[1][1]=fac_tauM_alphaM_velintaf[1]*derxy_(1,vi)+fac_tauM_afgdt_velintaf[1]*aux[1]+fac_tauM_afgdt_conv_resM_vi;
          temp[1][2]=fac_tauM_alphaM_velintaf[1]*derxy_(2,vi)+fac_tauM_afgdt_velintaf[1]*aux[2];
                                                                                     
          temp[2][0]=fac_tauM_alphaM_velintaf[2]*derxy_(0,vi)+fac_tauM_afgdt_velintaf[2]*aux[0];
          temp[2][1]=fac_tauM_alphaM_velintaf[2]*derxy_(1,vi)+fac_tauM_afgdt_velintaf[2]*aux[1];
          temp[2][2]=fac_tauM_alphaM_velintaf[2]*derxy_(2,vi)+fac_tauM_afgdt_velintaf[2]*aux[2]+fac_tauM_afgdt_conv_resM_vi;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;

	    elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
	    elemat(fvi  ,fuip ) += temp[0][1]*funct_(ui);
	    elemat(fvi  ,fuipp) += temp[0][2]*funct_(ui);
	                                     
	    elemat(fvip ,fui  ) += temp[1][0]*funct_(ui);
	    elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);
	    elemat(fvip ,fuipp) += temp[1][2]*funct_(ui);
	                                     
	    elemat(fvipp,fui  ) += temp[2][0]*funct_(ui);
	    elemat(fvipp,fuip ) += temp[2][1]*funct_(ui);
	    elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);

	  } // ui
	} // vi
	
        double fac_tauM_gdt_velintaf[3];
        fac_tauM_gdt_velintaf[0]=fac_tauM_gdt*velintaf_(0);
        fac_tauM_gdt_velintaf[1]=fac_tauM_gdt*velintaf_(1);
        fac_tauM_gdt_velintaf[2]=fac_tauM_gdt*velintaf_(2);

	for (int ui=0; ui<iel; ++ui) // loop columns 
	{
	  const int fuippp =4*ui+3;
          
	  for (int vi=0; vi<iel; ++vi)  // loop rows 
	  {
	    const int fvi   =4*vi;
	    const int fvip  =fvi+1;
	    const int fvipp =fvi+2;
            
            /* factor: tauM, rescaled by gamma*dt
	       
                         /                                  \
                        |    n+af    /          \            |
                        |   u     , |  nabla Dp  | o nabla v |
                        |            \          /            |
                         \                                  /
            */
	    const double aux=derxy_(0,vi)*derxy_(0,ui)+derxy_(1,vi)*derxy_(1,ui)+derxy_(2,vi)*derxy_(2,ui);
	    
	    elemat(fvi  ,fuippp) += fac_tauM_gdt_velintaf[0]*aux;
	    elemat(fvip ,fuippp) += fac_tauM_gdt_velintaf[1]*aux;
	    elemat(fvipp,fuippp) += fac_tauM_gdt_velintaf[2]*aux;
	  } // vi
	} // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
          
          /* factor: tauM*afgdt
	       
                  /                                               \
                 |    n+af    / /  n+af       \       \            |
                 |   u     , | |  u    o nabla | Dacc  | o nabla v |
                 |            \ \             /       /            |
                  \                                               /
          */
          double temp[3][3];

          temp[0][0]=fac_tauM_afgdt_velintaf[0]*derxy_(0,vi);
          temp[0][1]=fac_tauM_afgdt_velintaf[0]*derxy_(1,vi);
          temp[0][2]=fac_tauM_afgdt_velintaf[0]*derxy_(2,vi);
                                                                                                                       
          temp[1][0]=fac_tauM_afgdt_velintaf[1]*derxy_(0,vi);
          temp[1][1]=fac_tauM_afgdt_velintaf[1]*derxy_(1,vi);
          temp[1][2]=fac_tauM_afgdt_velintaf[1]*derxy_(2,vi);
                                                                                                                       
          temp[2][0]=fac_tauM_afgdt_velintaf[2]*derxy_(0,vi);
          temp[2][1]=fac_tauM_afgdt_velintaf[2]*derxy_(1,vi);
          temp[2][2]=fac_tauM_afgdt_velintaf[2]*derxy_(2,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
	    
	    elemat(fvi  ,fui  ) += temp[0][0]*conv_c_af_(ui);
	    elemat(fvi  ,fuip ) += temp[0][1]*conv_c_af_(ui);
	    elemat(fvi  ,fuipp) += temp[0][2]*conv_c_af_(ui);
                                             
	    elemat(fvip ,fui  ) += temp[1][0]*conv_c_af_(ui);
	    elemat(fvip ,fuip ) += temp[1][1]*conv_c_af_(ui);
	    elemat(fvip ,fuipp) += temp[1][2]*conv_c_af_(ui);
                                             
	    elemat(fvipp,fui  ) += temp[2][0]*conv_c_af_(ui);
	    elemat(fvipp,fuip ) += temp[2][1]*conv_c_af_(ui);
	    elemat(fvipp,fuipp) += temp[2][2]*conv_c_af_(ui);
	  } // ui
	} // vi

	if (higher_order_ele && newton!=Fluid3::minimal)
	{
	  const double fac_nu_afgdt_tauM=fac*visceff*afgdt*tauM;
          
	  double temp[3];

	  temp[0]=fac_nu_afgdt_tauM*velintaf_(0);
	  temp[1]=fac_nu_afgdt_tauM*velintaf_(1);
	  temp[2]=fac_nu_afgdt_tauM*velintaf_(2);


	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int fvi   =4*vi;
	    const int fvip  =fvi+1;
	    const int fvipp =fvi+2;
              
	    double rowtemp[3][3];
            
	    rowtemp[0][0]=temp[0]*derxy_(0,vi);
	    rowtemp[0][1]=temp[0]*derxy_(1,vi);
	    rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
	    rowtemp[1][0]=temp[1]*derxy_(0,vi);
	    rowtemp[1][1]=temp[1]*derxy_(1,vi);
	    rowtemp[1][2]=temp[1]*derxy_(2,vi);
            
	    rowtemp[2][0]=temp[2]*derxy_(0,vi);
	    rowtemp[2][1]=temp[2]*derxy_(1,vi);
	    rowtemp[2][2]=temp[2]*derxy_(2,vi);

	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int fui   =4*ui;
	      const int fuip  =fui+1;
	      const int fuipp =fui+2;
                
	      /* 
		 factor: - 2.0 * visc * alphaF * gamma * dt * tauM 

                    /                                                \
                   |   n+af   / /             /    \  \         \     |
                   |  u    , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /      
	      */
	      
	      elemat(fvi  ,fui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
	      elemat(fvi  ,fuip ) -= derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
	      elemat(fvi  ,fuipp) -= derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
	      elemat(fvip ,fui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
	      elemat(fvip ,fuip ) -= derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
	      elemat(fvip ,fuipp) -= derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
	      elemat(fvipp,fui  ) -= viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
	      elemat(fvipp,fuip ) -= derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
	      elemat(fvipp,fuipp) -= derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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

    double fac_inertia_gridconv_and_bodyforce[3];
    fac_inertia_gridconv_and_bodyforce[0] = fac*(accintam_(0)-convu_G_af_old_(0)-bodyforceaf_(0));
    fac_inertia_gridconv_and_bodyforce[1] = fac*(accintam_(1)-convu_G_af_old_(1)-bodyforceaf_(1));
    fac_inertia_gridconv_and_bodyforce[2] = fac*(accintam_(2)-convu_G_af_old_(2)-bodyforceaf_(2));
      
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
	
    double viscous_conv_and_pres[9];
    viscous_conv_and_pres[0]=visceff_fac*vderxyaf_(0,0)*2.0-fac_prenp;
    viscous_conv_and_pres[1]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_conv_and_pres[2]=visceff_fac*(vderxyaf_(0,2)+vderxyaf_(2,0));
    viscous_conv_and_pres[3]=visceff_fac*(vderxyaf_(0,1)+vderxyaf_(1,0));
    viscous_conv_and_pres[4]=visceff_fac*vderxyaf_(1,1)*2.0-fac_prenp;
    viscous_conv_and_pres[5]=visceff_fac*(vderxyaf_(1,2)+vderxyaf_(2,1));
    viscous_conv_and_pres[6]=visceff_fac*(vderxyaf_(0,2)+vderxyaf_(2,0));
    viscous_conv_and_pres[7]=visceff_fac*(vderxyaf_(1,2)+vderxyaf_(2,1));
    viscous_conv_and_pres[8]=visceff_fac*vderxyaf_(2,2)*2.0-fac_prenp;
     
    /*
          factor: -1.0

               /                                       \
              |   / n+af \     / n+af \           / \   |
              |  | u      | X | u      | , nabla | v |  |
              |   \      /     \      /           \ /   |
               \                                       /
    */
	
    if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
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
      viscous_conv_and_pres[2]-=velintaf_(0)*(-fac_tauM*resM_(2)+fac*velintaf_(2));
      viscous_conv_and_pres[3]-=velintaf_(1)*(-fac_tauM*resM_(0)+fac*velintaf_(0));
      viscous_conv_and_pres[4]-=velintaf_(1)*(-fac_tauM*resM_(1)+fac*velintaf_(1));
      viscous_conv_and_pres[5]-=velintaf_(1)*(-fac_tauM*resM_(2)+fac*velintaf_(2));
      viscous_conv_and_pres[6]-=velintaf_(2)*(-fac_tauM*resM_(0)+fac*velintaf_(0));
      viscous_conv_and_pres[7]-=velintaf_(2)*(-fac_tauM*resM_(1)+fac*velintaf_(1));
      viscous_conv_and_pres[8]-=velintaf_(2)*(-fac_tauM*resM_(2)+fac*velintaf_(2));
    }
    else
    {
      viscous_conv_and_pres[0]-=velintaf_(0)*velintaf_(0)*fac;
      viscous_conv_and_pres[1]-=velintaf_(0)*velintaf_(1)*fac;
      viscous_conv_and_pres[2]-=velintaf_(0)*velintaf_(2)*fac;
      viscous_conv_and_pres[3]-=velintaf_(1)*velintaf_(0)*fac;
      viscous_conv_and_pres[4]-=velintaf_(1)*velintaf_(1)*fac;
      viscous_conv_and_pres[5]-=velintaf_(1)*velintaf_(2)*fac;
      viscous_conv_and_pres[6]-=velintaf_(2)*velintaf_(0)*fac;
      viscous_conv_and_pres[7]-=velintaf_(2)*velintaf_(1)*fac;
      viscous_conv_and_pres[8]-=velintaf_(2)*velintaf_(2)*fac;
    }

    if(reynolds != Fluid3::reynolds_stress_stab_none)
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
	
      viscous_conv_and_pres[0]-=fac_tauM_tauM_resM_0  *resM_(0);
      viscous_conv_and_pres[1]-=fac_tauM_tauM_resM_0  *resM_(1);
      viscous_conv_and_pres[2]-=fac_tauM_tauM_resM_0  *resM_(2);
      viscous_conv_and_pres[3]-=fac_tauM_tauM_resM_0  *resM_(1);
      viscous_conv_and_pres[4]-=fac_tauM_tauM_resM_1  *resM_(1);
      viscous_conv_and_pres[5]-=fac_tauM_tauM_resM_1  *resM_(2);
      viscous_conv_and_pres[6]-=fac_tauM_tauM_resM_0  *resM_(2);
      viscous_conv_and_pres[7]-=fac_tauM_tauM_resM_1  *resM_(2);
      viscous_conv_and_pres[8]-=fac_tauM_tauM*resM_(2)*resM_(2);
    }
    
    /* continuity equation, factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
    */
    const double fac_divunp  = fac*divunp_;
    
    for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
    {
      int fui=4*ui;
      /* inertia, convective and dead load, cross terms fith 
	 funct                                              */
      /* viscous, pressure, reynolds, cstab terms with 
	 derxy                                              */

      elevec(fui++) -= 
	fac_inertia_gridconv_and_bodyforce[0]*funct_(ui)
	+
	derxy_(0,ui)*viscous_conv_and_pres[0]
	+
	derxy_(1,ui)*viscous_conv_and_pres[1]
	+
	derxy_(2,ui)*viscous_conv_and_pres[2];
      elevec(fui++) -= 
	fac_inertia_gridconv_and_bodyforce[1]*funct_(ui)
	+
	derxy_(0,ui)*viscous_conv_and_pres[3]
	+
	derxy_(1,ui)*viscous_conv_and_pres[4]
	+
	derxy_(2,ui)*viscous_conv_and_pres[5];
      elevec(fui++) -= 
	fac_inertia_gridconv_and_bodyforce[2]*funct_(ui)
	+
	derxy_(0,ui)*viscous_conv_and_pres[6]
	+
	derxy_(1,ui)*viscous_conv_and_pres[7]
	+
	derxy_(2,ui)*viscous_conv_and_pres[8];
	  
      /* continuity equation */
      elevec(fui  ) -= fac_divunp*funct_(ui);
    }

    if(pspg == Fluid3::pstab_use_pspg)
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
      
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
	elevec(4*ui+3)-=fac_tauMp*(resM_(0)*derxy_(0,ui)+resM_(1)*derxy_(1,ui)+resM_(2)*derxy_(2,ui));
      }
    } // end pspg
    
    if(supg == Fluid3::convective_stab_supg)
    {
      const double fac_tauM = fac*tauM;
      
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
	int fui=4*ui;
	
	const double fac_tauM_conv_c_af_ui = fac_tauM*conv_c_af_(ui);
	/*
	  factor: +tauM
	  
	  SUPG stabilisation


                  /                             \
                 |   n+af    / n+af        \     |
                 |  r     , | c     o nabla | v  |
                 |   M       \             /     |
                  \                             /
	*/
	
	elevec(fui++) -= fac_tauM_conv_c_af_ui*resM_(0);
	elevec(fui++) -= fac_tauM_conv_c_af_ui*resM_(1);
	elevec(fui  ) -= fac_tauM_conv_c_af_ui*resM_(2);
	
      } // end loop rows
    } // end supg

    if (higher_order_ele)
    {
      if(vstab != Fluid3::viscous_stab_none && higher_order_ele)
      {
	const double fac_visc_tauMp = vstabfac * fac*visc*tauMp;
	
	for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
	{
	  int fui=4*ui;
	  /*
	    factor: -(+)tauMp*nu

	    viscous stabilisation --- inertia


                 /                      \
                |   n+af                 |
                |  r    , 2*div eps (v)  |
                |   M                    |
                 \                      /

	  */
	  elevec(fui++) -= fac_visc_tauMp*
	                   (resM_(0)*viscs2_(0,ui)
			    +
			    resM_(1)*derxy2_(3,ui)
			    +
			    resM_(2)*derxy2_(4,ui));
	  elevec(fui++) -= fac_visc_tauMp*
                           (resM_(0)*derxy2_(3,ui)
			    +
			    resM_(1)*viscs2_(1,ui)
			    +
			    resM_(2)*derxy2_(5,ui));
	  elevec(fui  ) -= fac_visc_tauMp*
	                   (resM_(0)*derxy2_(4,ui)
			    +
			    resM_(1)*derxy2_(5,ui)
			    +
			    resM_(2)*viscs2_(2,ui));
	} // end loop rows ui
      } // endif (a)gls
    } // end hoel

    if(fssgv != Fluid3::fssgv_no)
    {
      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)
      
      for (int ui=0; ui<iel; ++ui)
      {
	/* fine-scale subgrid-viscosity term on right hand side */
	/*
                                  /                              \
                         n+af    |       /    n+af\         / \   |
             - nu_art(fsu    ) * |  eps | Dfsu     | , eps | v |  |
                                 |       \        /         \ /   |
                                  \                              /
	*/
	elevec(ui*4    ) -= vartfac*(2.0*derxy_(0, ui)*fsvderxyaf_(0, 0)
				     +    derxy_(1, ui)*fsvderxyaf_(0, 1)
				     +    derxy_(1, ui)*fsvderxyaf_(1, 0)
				     +    derxy_(2, ui)*fsvderxyaf_(0, 2)
				     +    derxy_(2, ui)*fsvderxyaf_(2, 0));
	elevec(ui*4 + 1) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 1)
				     +    derxy_(0, ui)*fsvderxyaf_(1, 0)
				     +2.0*derxy_(1, ui)*fsvderxyaf_(1, 1)
				     +    derxy_(2, ui)*fsvderxyaf_(1, 2)
				     +    derxy_(2, ui)*fsvderxyaf_(2, 1));
        elevec(ui*4 + 2) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 2)
				     +    derxy_(0, ui)*fsvderxyaf_(2, 0)
				     +    derxy_(1, ui)*fsvderxyaf_(1, 2)
				     +    derxy_(1, ui)*fsvderxyaf_(2, 1)
				     +2.0*derxy_(2, ui)*fsvderxyaf_(2, 2));
      } // end loop ui
    } // end not fssgv_no
  } // end loop iquad
  return;
} // end Sysmat_cons_qs


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  conservative, time-dependent version                               |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Sysmat_cons_td(
  Fluid3*                                          ele             ,
  std::vector<blitz::Array<double,1> >&            myknots         ,
  LINALG::Matrix<4*iel,4*iel>&                     elemat          ,
  LINALG::Matrix<4*iel,1>&                         elevec          ,
  const LINALG::Matrix<3,iel>&                     edispnp         ,
  const LINALG::Matrix<3,iel>&                     egridvaf        ,
  const LINALG::Matrix<3,iel>&                     evelnp          ,
  const LINALG::Matrix<iel,1>&                     eprenp          ,
  const LINALG::Matrix<3,iel>&                     eaccam          ,
  const LINALG::Matrix<3,iel>&                     evelaf          ,
  const LINALG::Matrix<3,iel>&                     fsevelaf        ,
  const struct _MATERIAL*                          material        ,
  const double                                     alphaM          ,
  const double                                     alphaF          ,
  const double                                     gamma           ,
  const double                                     dt              ,
  const double                                     time            ,
  const enum Fluid3::LinearisationAction           newton          ,
  const bool                                       higher_order_ele,
  const enum Fluid3::StabilisationAction           fssgv           ,
  const enum Fluid3::StabilisationAction           inertia         ,
  const enum Fluid3::StabilisationAction           pspg            ,
  const enum Fluid3::StabilisationAction           supg            ,
  const enum Fluid3::StabilisationAction           vstab           ,
  const enum Fluid3::StabilisationAction           cstab           ,
  const enum Fluid3::StabilisationAction           cross           ,
  const enum Fluid3::StabilisationAction           reynolds        ,
  const enum Fluid3::TauType                       whichtau        ,
  const enum Fluid3::TurbModelAction               turb_mod_action ,
  double&                                          Cs              ,
  double&                                          Cs_delta_sq     ,
  double&                                          visceff         ,
  const double                                     l_tau           ,
  const bool                                       compute_elemat
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
  if (vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_usfem_only_rhs)
  {
    vstabfac =  1.0;
  }
  else if(vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_gls_only_rhs)
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

  SetElementData(ele            ,
                 edispnp        ,
                 evelaf         ,
                 fsevelaf       ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 fssgv          ,
                 turb_mod_action,
                 l_tau          ,
                 Cs             ,
                 Cs_delta_sq    ,
                 visceff        );

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero
  if(ele->sub_acc_old_.extent(blitz::firstDim) != 3 
     ||
     ele->sub_acc_old_.extent(blitz::secondDim) != intpoints.nquad)
  {
    ele->sub_acc_old_ .resize(3,intpoints.nquad);
    ele->sub_acc_old_  = 0.;
  }
  if(ele->sub_vel_old_.extent(blitz::firstDim) != 3 
     || 
     ele->sub_vel_old_.extent(blitz::secondDim) != intpoints.nquad)
  {
    ele->sub_vel_old_ .resize(3,intpoints.nquad);
    ele->sub_vel_old_  = 0.;
    
    ele->sub_vel_.resize(3,intpoints.nquad);
    ele->sub_vel_ = 0.;
  }
  
  // get subscale information from element --- this is just a reference
  // to the element data
  blitz::Array<double,2> saccn (ele->sub_acc_old_);
  blitz::Array<double,2> sveln (ele->sub_vel_old_);
  blitz::Array<double,2> svelnp(ele->sub_vel_    );

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second 
    // derivatives in a gausspoint and integration weight including 
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);
    
    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint(ele             ,
                            egridvaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,Fluid3::subscales_time_dependent,gamma,dt,hk,mk,visceff);

    // stabilisation parameters
    const double tauM   = tau_(0);

    if(cstab == Fluid3::continuity_stab_none)
    {
      tau_(2)=0.0;
    }
    const double tauC    = tau_(2);

    double supg_active;
    if(supg == Fluid3::convective_stab_supg)
    {
      supg_active=1.0;
    }
    else
    {
      supg_active=0.0;
    }

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;
    // update estimates for the subscale quantities
    const double facMtau = 1./(alphaM*tauM+afgdt);

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
    const double fac1=(alphaM*tauM+gamma*dt*(alphaF-1.0))*facMtau;
    const double fac2=(dt*tauM*(alphaM-gamma))*facMtau;
    const double fac3=(gamma*dt*tauM)*facMtau;
    
    for (int rr=0;rr<3;++rr)
    {
      svelnp(rr,iquad)= fac1*sveln(rr,iquad)+fac2*saccn(rr,iquad)-fac3*resM_(rr);
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
      svelaf_(rr)=alphaF*svelnp(rr,iquad)+(1.0-alphaF)*sveln(rr,iquad);
    }

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
    if(cross    == Fluid3::cross_stress_stab
       ||
       reynolds == Fluid3::reynolds_stress_stab)
    {
      for (int rr=0;rr<iel;++rr)
      {
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);
	
        for (int mm=1;mm<3;++mm)
        {
          conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
        }
      }

      if(reynolds == Fluid3::reynolds_stress_stab)
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
    if(cross == Fluid3::cross_stress_stab_only_rhs
       ||
       cross == Fluid3::cross_stress_stab
      )
    {
      for (int rr=0;rr<3;++rr)
      {
        convsubaf_old_(rr) = vderxyaf_(rr, 0)*svelaf_(0);
	
        for (int mm=1;mm<3;++mm)
        {
          convsubaf_old_(rr) += vderxyaf_(rr, mm)*svelaf_(mm);
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
    if(ele->is_ale_)
    {
      for(int nn=0;nn<iel;++nn)
      {
	conv_u_G_af_(nn)=u_G_af_(0)*derxy_(0,nn);
        
	for(int rr=1;rr<3;++rr)
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
    for(int rr=0;rr<3;++rr)
    {
      convu_G_af_old_(rr)=u_G_af_(0)*vderxyaf_(rr,0);
      for(int mm=1;mm<3;++mm)
      {
        convu_G_af_old_(rr)+=u_G_af_(mm)*vderxyaf_(rr,mm);
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
	  
      if(inertia == Fluid3::inertia_stab_keep
         || 
         inertia == Fluid3::inertia_stab_keep_complete)
      {
        // rescale time factors terms affected by inertia stabilisation
        fac_inertia*=afgdt*facMtau;
	  
        // do inertia stabilisation terms which are not scaled
        // Galerkin terms since they are not partially integrated

        const double fac_alphaM_tauM_facMtau = fac*alphaM*tauM*facMtau;
	  
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
	    
          const double fac_alphaM_gamma_dt_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*gamma*dt*funct_(vi);
	    
          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fuippp  =4*ui+3;
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
	      
            elemat(fvi  ,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(0,ui);
            elemat(fvip ,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(1,ui);
            elemat(fvipp,fuippp) -= fac_alphaM_gamma_dt_tauM_facMtau_funct_vi*derxy_(2,ui);
          } // ui
        } // vi

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
            
          for (int ui=0; ui<iel; ++ui) // loop columns
          {	  
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fui+2;

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

              elemat(fvi  ,fui  ) -= afgdt*fac_alphaM_tauM_facMtau*conv_c_af_(ui)*funct_(vi);
              elemat(fvip ,fuip ) -= afgdt*fac_alphaM_tauM_facMtau*conv_c_af_(ui)*funct_(vi);
              elemat(fvipp,fuipp) -= afgdt*fac_alphaM_tauM_facMtau*conv_c_af_(ui)*funct_(vi);
          }
        }
        if(newton==Fluid3::Newton)
        {
          double temp[3][3];

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi    =4*vi;
            const int fvip   =fvi+1;
            const int fvipp  =fvi+2;
          
            const double aux=afgdt*fac_alphaM_tauM_facMtau*funct_(vi);

            temp[0][0]=aux*vderxyaf_(0,0);
            temp[1][0]=aux*vderxyaf_(0,1);
            temp[2][0]=aux*vderxyaf_(0,2);
            temp[0][1]=aux*vderxyaf_(1,0);
            temp[1][1]=aux*vderxyaf_(1,1);
            temp[2][1]=aux*vderxyaf_(1,2);
            temp[0][2]=aux*vderxyaf_(2,0);
            temp[1][2]=aux*vderxyaf_(2,1);
            temp[2][2]=aux*vderxyaf_(2,2);

  
            for (int ui=0; ui<iel; ++ui) // loop columns
            {	  
              const int fui    =4*ui;
              const int fuip   =fui+1;
              const int fuipp  =fui+2;
	  
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


              elemat(fvi  ,fui  ) -= temp[0][0]*funct_(ui);
              elemat(fvi  ,fuip ) -= temp[1][0]*funct_(ui);
              elemat(fvi  ,fuipp) -= temp[2][0]*funct_(ui);
              elemat(fvip ,fui  ) -= temp[0][1]*funct_(ui);
              elemat(fvip ,fuip ) -= temp[1][1]*funct_(ui);
              elemat(fvip ,fuipp) -= temp[2][1]*funct_(ui);
              elemat(fvipp,fui  ) -= temp[0][2]*funct_(ui);
              elemat(fvipp,fuip ) -= temp[1][2]*funct_(ui);
              elemat(fvipp,fuipp) -= temp[2][2]*funct_(ui);

            }
          }
        }

        if(higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_alphaM_tauM_facMtau
            =
            fac*visceff*afgdt*alphaM*tauM*facMtau;
	    
          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi    =4*vi;
            const int fvip   =fvi+1;
            const int fvipp  =fvi+2;
		
            const double fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
              =
              fac_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);
		  
            for (int ui=0; ui<iel; ++ui) // loop columns
            {	  
              const int fui    =4*ui;
              const int fuip   =fui+1;
              const int fuipp  =fui+2;
	  
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
              const double a = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(3,ui);
              const double b = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(4,ui);
              const double c = fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*derxy2_(5,ui);

              elemat(fvi  ,fui  ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(0,ui);
              elemat(fvi  ,fuip ) += a;
              elemat(fvi  ,fuipp) += b;
              elemat(fvip ,fui  ) += a;
              elemat(fvip ,fuip ) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(1,ui);
              elemat(fvip ,fuipp) += c;
              elemat(fvipp,fui  ) += b;
              elemat(fvipp,fuip ) += c;
              elemat(fvipp,fuipp) += fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi*viscs2_(2,ui);
            } // ui
          } // vi
        } // end higher order element and  linearisation of linear terms not supressed

        if(inertia == Fluid3::inertia_stab_keep_complete)
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
          const double norm = sqrt(velintaf_(0)*velintaf_(0)+velintaf_(1)*velintaf_(1)+velintaf_(2)*velintaf_(2));

          // normed velocity at element center (we use the copy for safety reasons!)
          if (norm>=1e-6)
          {
            for (int rr=0;rr<3;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=velintaf_(rr)/norm;
            }
          }
          else
          {
            normed_velintaf_(0) = 0.0;
            for (int rr=1;rr<3;++rr) /* loop element nodes */
            {
              normed_velintaf_(rr)=0.0;
            }      
          }
	      
          double temp=0.0;
          if(whichtau==Fluid3::codina)
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
          else if(whichtau==Fluid3::smoothed_franca_barrenechea_valentin_wall)
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
          else if(whichtau==Fluid3::franca_barrenechea_valentin_wall)
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
            const int fvi    =4*vi;
            const int fvip   =fvi+1;
            const int fvipp  =fvi+2;

              
            for (int ui=0; ui<iel; ++ui) // loop columns
            {
              const int fui  =4*ui;
              const int fuip =fui+1;
              const int fuipp=fui+2;

              elemat(fvi  ,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(fvi  ,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(0);
              elemat(fvi  ,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(0);

              elemat(fvip ,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(fvip ,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(1);
              elemat(fvip ,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(1);

              elemat(fvipp,fui  ) -= temp*normed_velintaf_(0)*funct_(ui)*funct_(vi)*svelaf_(2);
              elemat(fvipp,fuip ) -= temp*normed_velintaf_(1)*funct_(ui)*funct_(vi)*svelaf_(2);
              elemat(fvipp,fuipp) -= temp*normed_velintaf_(2)*funct_(ui)*funct_(vi)*svelaf_(2);
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
      const double fac_afgdt_velintaf_z=fac_afgdt*velintaf_(2);

      // supg and cstab conservative
      const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;
        
      for (int ui=0; ui<iel; ++ui) // loop columns 
      {
        const int fui  =4*ui;
        const int fuip =fui+1;
        const int fuipp=fui+2;
        
        /* GALERKIN inertia term (intermediate) + convection, mesh velocity (intermediate) */
        const double inertia_and_gridconv_ui = fac_inertia*funct_(ui)-fac_afgdt*conv_u_G_af_(ui);

        /* SUPG stabilisation --- inertia and convection */
        const double supg_inertia_and_conv_ui
          = fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_afgdt_tauM_facMtau*conv_c_af_(ui);
          
        // convection GALERKIN and diagonal parts of viscous term (intermediate)
        const double convection_and_viscous_x=fac_visceff_afgdt*derxy_(0,ui)-fac_afgdt_velintaf_x*funct_(ui);
        const double convection_and_viscous_y=fac_visceff_afgdt*derxy_(1,ui)-fac_afgdt_velintaf_y*funct_(ui);
        const double convection_and_viscous_z=fac_visceff_afgdt*derxy_(2,ui)-fac_afgdt_velintaf_z*funct_(ui);
        
        // viscous GALERKIN term
        const double viscous_x=fac_visceff_afgdt*derxy_(0,ui);
        const double viscous_y=fac_visceff_afgdt*derxy_(1,ui);
        const double viscous_z=fac_visceff_afgdt*derxy_(2,ui);
        
        /* CSTAB entries */
        const double fac_gamma_dt_tauC_derxy_x_ui = fac_gamma_dt_tauC*derxy_(0,ui);
        const double fac_gamma_dt_tauC_derxy_y_ui = fac_gamma_dt_tauC*derxy_(1,ui);
        const double fac_gamma_dt_tauC_derxy_z_ui = fac_gamma_dt_tauC*derxy_(2,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int fvi    =4*vi;
          const int fvip   =fvi+1;
          const int fvipp  =fvi+2;
              
          const double sum =
            inertia_and_gridconv_ui*funct_(vi)
            +
            supg_inertia_and_conv_ui*conv_c_plus_svel_af_(vi)
            +
            convection_and_viscous_x*derxy_(0,vi)
            +
            convection_and_viscous_y*derxy_(1,vi)
            +
            convection_and_viscous_z*derxy_(2,vi);
          
          /* adding GALERKIN convection, convective linearisation (intermediate), viscous and cstab */
          
          elemat(fvi  ,fui  ) += sum+(fac_gamma_dt_tauC_derxy_x_ui             + viscous_x)*derxy_(0,vi);
          elemat(fvi  ,fuip ) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(0,vi)+(viscous_x)*derxy_(1,vi);
          elemat(fvi  ,fuipp) +=      fac_gamma_dt_tauC_derxy_z_ui*derxy_(0,vi)+(viscous_x)*derxy_(2,vi);
          elemat(fvip ,fui  ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(1,vi)+(viscous_y)*derxy_(0,vi);
          elemat(fvip ,fuip ) += sum+(fac_gamma_dt_tauC_derxy_y_ui             + viscous_y)*derxy_(1,vi);
          elemat(fvip ,fuipp) +=      fac_gamma_dt_tauC_derxy_z_ui*derxy_(1,vi)+(viscous_y)*derxy_(2,vi);
          elemat(fvipp,fui  ) +=      fac_gamma_dt_tauC_derxy_x_ui*derxy_(2,vi)+(viscous_z)*derxy_(0,vi);
          elemat(fvipp,fuip ) +=      fac_gamma_dt_tauC_derxy_y_ui*derxy_(2,vi)+(viscous_z)*derxy_(1,vi);
          elemat(fvipp,fuipp) += sum+(fac_gamma_dt_tauC_derxy_z_ui             + viscous_z)*derxy_(2,vi);
        } // vi
      } // ui
      
      for (int ui=0; ui<iel; ++ui) // loop columns
      {
        const int fuippp  =4*ui+3;

        const double fac_gamma_dt_funct_ui=fac_gamma_dt*funct_(ui);
	  
        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi  =4*vi;
          const int fvip =fvi+1;
          const int fvipp=fvi+2;

          /* GALERKIN pressure   (implicit), rescaled by gamma*dt */
          /* continuity equation (implicit)                       */

          elemat(fvi   ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(fvip  ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(1,vi);
          elemat(fvipp ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(2,vi);

          elemat(fuippp,fvi   ) += fac_gamma_dt_funct_ui*derxy_(0,vi);
          elemat(fuippp,fvip  ) += fac_gamma_dt_funct_ui*derxy_(1,vi);
          elemat(fuippp,fvipp ) += fac_gamma_dt_funct_ui*derxy_(2,vi);
        } // vi
      } // ui

      if (newton==Fluid3::Newton) // if newton and supg
      {
        const double fac_afgdt_afgdt_tauM_facMtau = fac*afgdt*afgdt*facMtau*tauM;

        // linearisation of SUPG testfunction and GALERKIN reactive part of convection
        double temp[3][3];

        const double fac_afgdt_svelaf_0 = fac*afgdt*supg_active*svelaf_(0)+fac*afgdt*velintaf_(0);
        const double fac_afgdt_svelaf_1 = fac*afgdt*supg_active*svelaf_(1)+fac*afgdt*velintaf_(1);
        const double fac_afgdt_svelaf_2 = fac*afgdt*supg_active*svelaf_(2)+fac*afgdt*velintaf_(2);

        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi  =4*vi;
          const int fvip =fvi+1;
          const int fvipp=fvi+2;
            
          // SUPG part (reactive part from residual)
          const double scaled_inertia_and_conv_vi 
            = 
            fac_afgdt_afgdt_tauM_facMtau*conv_c_plus_svel_af_(vi);
	      
          temp[0][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,0)-fac_afgdt_svelaf_0*derxy_(0,vi);
          temp[1][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,1)-fac_afgdt_svelaf_0*derxy_(1,vi);
          temp[2][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,2)-fac_afgdt_svelaf_0*derxy_(2,vi);
          temp[0][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,0)-fac_afgdt_svelaf_1*derxy_(0,vi);
          temp[1][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,1)-fac_afgdt_svelaf_1*derxy_(1,vi);
          temp[2][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,2)-fac_afgdt_svelaf_1*derxy_(2,vi);
          temp[0][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,0)-fac_afgdt_svelaf_2*derxy_(0,vi);
          temp[1][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,1)-fac_afgdt_svelaf_2*derxy_(1,vi);
          temp[2][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,2)-fac_afgdt_svelaf_2*derxy_(2,vi);
	      
          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
            elemat(fvi  ,fuip ) += temp[1][0]*funct_(ui);
            elemat(fvi  ,fuipp) += temp[2][0]*funct_(ui);
            elemat(fvip ,fui  ) += temp[0][1]*funct_(ui);
            elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);
            elemat(fvip ,fuipp) += temp[2][1]*funct_(ui);
            elemat(fvipp,fui  ) += temp[0][2]*funct_(ui);
            elemat(fvipp,fuip ) += temp[1][2]*funct_(ui);
            elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);
          } // ui
        } // vi
      } // end if newton

      for (int ui=0; ui<iel; ++ui) // loop columns 
      {
        const int fuippp=4*ui+3;

        const double scaled_gradp_0 = fac_gdt_afgdt_tauM_facMtau*derxy_(0,ui);
        const double scaled_gradp_1 = fac_gdt_afgdt_tauM_facMtau*derxy_(1,ui);
        const double scaled_gradp_2 = fac_gdt_afgdt_tauM_facMtau*derxy_(2,ui);

        for (int vi=0; vi<iel; ++vi)  // loop rows 
        {
          const int fvi=4*vi;
 
          /* SUPG stabilisation --- pressure, rescaled by gamma*dt */
          elemat(fvi  ,fuippp) += scaled_gradp_0*conv_c_plus_svel_af_(vi);
          elemat(fvi+1,fuippp) += scaled_gradp_1*conv_c_plus_svel_af_(vi);
          elemat(fvi+2,fuippp) += scaled_gradp_2*conv_c_plus_svel_af_(vi);
        } // vi
      } // ui

      if(higher_order_ele && newton!=Fluid3::minimal)
      {
        const double fac_visceff_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;

        for (int ui=0; ui<iel; ++ui) // loop columns
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
    
          double coltemp[3][3];
              
          coltemp[0][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(0,ui);
          coltemp[0][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(3,ui);
          coltemp[0][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(4,ui);
                                                                                                                     
          coltemp[1][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(3,ui);
          coltemp[1][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(1,ui);
          coltemp[1][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(5,ui);
                                                                                                                     
          coltemp[2][0]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(4,ui);
          coltemp[2][1]=fac_visceff_afgdt_afgdt_tauM_facMtau*derxy2_(5,ui);
          coltemp[2][2]=fac_visceff_afgdt_afgdt_tauM_facMtau*viscs2_(2,ui);

          for (int vi=0; vi<iel; ++vi)  // loop rows
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

            /*  SUPG stabilisation, diffusion */
            elemat(fvi  ,fui  ) -= coltemp[0][0]*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuip ) -= coltemp[0][1]*conv_c_plus_svel_af_(vi);
            elemat(fvi  ,fuipp) -= coltemp[0][2]*conv_c_plus_svel_af_(vi);
            
            elemat(fvip ,fui  ) -= coltemp[1][0]*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuip ) -= coltemp[1][1]*conv_c_plus_svel_af_(vi);
            elemat(fvip ,fuipp) -= coltemp[1][2]*conv_c_plus_svel_af_(vi);
            
            elemat(fvipp,fui  ) -= coltemp[2][0]*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuip ) -= coltemp[2][1]*conv_c_plus_svel_af_(vi);
            elemat(fvipp,fuipp) -= coltemp[2][2]*conv_c_plus_svel_af_(vi);
          } // vi
        } // ui
      } // hoel

      if(reynolds == Fluid3::reynolds_stress_stab)
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
        const double fac_alphaM_afgdt_tauM_facMtau_svelaf_z = fac_alphaM_afgdt_tauM_facMtau*svelaf_(2);

        const double fac_afgdt_afgdt_tauM_facMtau =fac*afgdt*afgdt*tauM*facMtau;

        double fac_afgdt_afgdt_tauM_facMtau_svelaf[3];
        fac_afgdt_afgdt_tauM_facMtau_svelaf[0]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[1]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(1);
        fac_afgdt_afgdt_tauM_facMtau_svelaf[2]=fac_afgdt_afgdt_tauM_facMtau*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fui   =4*ui;
          const int fuip  =fui+1;
          const int fuipp =fui+2;
            
          const double u_o_nabla_ui=velintaf_(0)*derxy_(0,ui)+velintaf_(1)*derxy_(1,ui)+velintaf_(2)*derxy_(2,ui);
            
          double inertia_and_conv[3];

          inertia_and_conv[0]=fac_afgdt_afgdt_tauM_facMtau_svelaf[0]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_x*funct_(ui);
          inertia_and_conv[1]=fac_afgdt_afgdt_tauM_facMtau_svelaf[1]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_y*funct_(ui);
          inertia_and_conv[2]=fac_afgdt_afgdt_tauM_facMtau_svelaf[2]*u_o_nabla_ui+fac_alphaM_afgdt_tauM_facMtau_svelaf_z*funct_(ui);
              
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;

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

            elemat(fvi  ,fui  ) += inertia_and_conv[0]*derxy_(0,vi);
            elemat(fvi  ,fuip ) += inertia_and_conv[0]*derxy_(1,vi);
            elemat(fvi  ,fuipp) += inertia_and_conv[0]*derxy_(2,vi);
                                     
            elemat(fvip ,fui  ) += inertia_and_conv[1]*derxy_(0,vi);
            elemat(fvip ,fuip ) += inertia_and_conv[1]*derxy_(1,vi);
            elemat(fvip ,fuipp) += inertia_and_conv[1]*derxy_(2,vi);
                                     
            elemat(fvipp,fui  ) += inertia_and_conv[2]*derxy_(0,vi);
            elemat(fvipp,fuip ) += inertia_and_conv[2]*derxy_(1,vi);
            elemat(fvipp,fuipp) += inertia_and_conv[2]*derxy_(2,vi);
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
            
          double temp[3];
          temp[0]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi));
          temp[1]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi));
          temp[2]=fac_afgdt_afgdt_tauM_facMtau*(vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi));

          double rowtemp[3][3];
            
          rowtemp[0][0]=svelaf_(0)*temp[0];
          rowtemp[0][1]=svelaf_(0)*temp[1];
          rowtemp[0][2]=svelaf_(0)*temp[2];
                               
          rowtemp[1][0]=svelaf_(1)*temp[0];
          rowtemp[1][1]=svelaf_(1)*temp[1];
          rowtemp[1][2]=svelaf_(1)*temp[2];
                               
          rowtemp[2][0]=svelaf_(2)*temp[0];
          rowtemp[2][1]=svelaf_(2)*temp[1];
          rowtemp[2][2]=svelaf_(2)*temp[2];
              
          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
            
            /* 
                 factor: + alphaF * gamma * dt * alphaF * gamma * dt * tauM *facMtau

              /                                                   \
             |  ~n+af    / / /            \   n+af \         \     |
             |  u     , | | | Dacc o nabla | u      | o nabla | v  |
             |           \ \ \            /        /         /     |
              \                                                   /      

            */

            elemat(fvi  ,fui  ) += funct_(ui)*rowtemp[0][0];
            elemat(fvi  ,fuip ) += funct_(ui)*rowtemp[0][1];
            elemat(fvi  ,fuipp) += funct_(ui)*rowtemp[0][2];
                                                             
            elemat(fvip ,fui  ) += funct_(ui)*rowtemp[1][0];
            elemat(fvip ,fuip ) += funct_(ui)*rowtemp[1][1];
            elemat(fvip ,fuipp) += funct_(ui)*rowtemp[1][2];
                                                             
            elemat(fvipp,fui  ) += funct_(ui)*rowtemp[2][0];
            elemat(fvipp,fuip ) += funct_(ui)*rowtemp[2][1];
            elemat(fvipp,fuipp) += funct_(ui)*rowtemp[2][2];
          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        const double fac_gdt_afgdt_tauM_facMtau         =fac*gamma*dt*afgdt*tauM*facMtau;
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_x=fac_gdt_afgdt_tauM_facMtau*svelaf_(0);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_y=fac_gdt_afgdt_tauM_facMtau*svelaf_(1);
        const double fac_gdt_afgdt_tauM_facMtau_svelaf_z=fac_gdt_afgdt_tauM_facMtau*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fuippp =4*ui+3;
            
          double coltemp[3][3];

          coltemp[0][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(0,ui);
          coltemp[0][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(1,ui);
          coltemp[0][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_x*derxy_(2,ui);
          coltemp[1][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(0,ui);
          coltemp[1][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(1,ui);
          coltemp[1][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_y*derxy_(2,ui);
          coltemp[2][0]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(0,ui);
          coltemp[2][1]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(1,ui);
          coltemp[2][2]=fac_gdt_afgdt_tauM_facMtau_svelaf_z*derxy_(2,ui);
  
          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
            
            /* 
                 factor: + gamma * dt * alphaF * gamma * dt * tauM *facMtau (rescaled)
             
              /                                \
             |  ~n+af    /                \     |
             |  u     , | nabla Dp o nabla | v  |
             |           \                /     |
              \                                /      

            */

            elemat(fvi  ,fuippp) += coltemp[0][0]*derxy_(0,vi)+coltemp[0][1]*derxy_(1,vi)+coltemp[0][2]*derxy_(2,vi);
            elemat(fvip ,fuippp) += coltemp[1][0]*derxy_(0,vi)+coltemp[1][1]*derxy_(1,vi)+coltemp[1][2]*derxy_(2,vi);
            elemat(fvipp,fuippp) += coltemp[2][0]*derxy_(0,vi)+coltemp[2][1]*derxy_(1,vi)+coltemp[2][2]*derxy_(2,vi);

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)


        if (higher_order_ele && newton!=Fluid3::minimal)
        {
          const double fac_nu_afgdt_afgdt_tauM_facMtau =fac*visceff*afgdt*afgdt*tauM*facMtau;
          
          double temp[3];

          temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(0);
          temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(1);
          temp[2]=fac_nu_afgdt_afgdt_tauM_facMtau*svelaf_(2);

          for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
          {
            const int fvi   =4*vi;
            const int fvip  =fvi+1;
            const int fvipp =fvi+2;
              
            double rowtemp[3][3];
              
            rowtemp[0][0]=temp[0]*derxy_(0,vi);
            rowtemp[0][1]=temp[0]*derxy_(1,vi);
            rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
            rowtemp[1][0]=temp[1]*derxy_(0,vi);
            rowtemp[1][1]=temp[1]*derxy_(1,vi);
            rowtemp[1][2]=temp[1]*derxy_(2,vi);
              
            rowtemp[2][0]=temp[2]*derxy_(0,vi);
            rowtemp[2][1]=temp[2]*derxy_(1,vi);
            rowtemp[2][2]=temp[2]*derxy_(2,vi);

            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
              const int fui   =4*ui;
              const int fuip  =fui+1;
              const int fuipp =fui+2;
                
              /* 
                   factor: - 2.0 * visc * alphaF * gamma * dt * alphaF * gamma * dt * tauM * facMtauM

                    /                                                 \
                   |  ~n+af    / /             /    \  \         \     |
                   |  u     , | | nabla o eps | Dacc |  | o nabla | v  |
                   |           \ \             \    /  /         /     |
                    \                                                 /      
              */

              elemat(fvi  ,fui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
              elemat(fvi  ,fuip ) -= derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
              elemat(fvi  ,fuipp) -= derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
              elemat(fvip ,fui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
              elemat(fvip ,fuip ) -= derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
              elemat(fvip ,fuipp) -= derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
              elemat(fvipp,fui  ) -= viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
              elemat(fvipp,fuip ) -= derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
              elemat(fvipp,fuipp) -= derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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
      if(cross == Fluid3::cross_stress_stab)
      {
        const double fac_afgdt_afgdt_tauM_facMtau  = fac*afgdt   *afgdt*tauM*facMtau;
        const double fac_gdt_afgdt_tauM_facMtau    = fac*gamma*dt*afgdt*tauM*facMtau;
        const double fac_alphaM_afgdt_tauM_facMtau = fac*alphaM  *afgdt*tauM*facMtau;

        double fac_alphaM_afgdt_tauM_velintaf[3];
        fac_alphaM_afgdt_tauM_velintaf[0]=fac_alphaM_afgdt_tauM_facMtau*velintaf_(0);
        fac_alphaM_afgdt_tauM_velintaf[1]=fac_alphaM_afgdt_tauM_facMtau*velintaf_(1);
        fac_alphaM_afgdt_tauM_velintaf[2]=fac_alphaM_afgdt_tauM_facMtau*velintaf_(2);

        double fac_afgdt_afgdt_tauM_facMtau_velintaf[3];
        fac_afgdt_afgdt_tauM_facMtau_velintaf[0]=fac_afgdt_afgdt_tauM_facMtau*velintaf_(0);
        fac_afgdt_afgdt_tauM_facMtau_velintaf[1]=fac_afgdt_afgdt_tauM_facMtau*velintaf_(1);
        fac_afgdt_afgdt_tauM_facMtau_velintaf[2]=fac_afgdt_afgdt_tauM_facMtau*velintaf_(2);

        for (int vi=0; vi<iel; ++vi)  // loop rows 
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;

          
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
         
          double aux[3];

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
          aux[0]=vderxyaf_(0,0)*derxy_(0,vi)+vderxyaf_(1,0)*derxy_(1,vi)+vderxyaf_(2,0)*derxy_(2,vi);
          aux[1]=vderxyaf_(0,1)*derxy_(0,vi)+vderxyaf_(1,1)*derxy_(1,vi)+vderxyaf_(2,1)*derxy_(2,vi);
          aux[2]=vderxyaf_(0,2)*derxy_(0,vi)+vderxyaf_(1,2)*derxy_(1,vi)+vderxyaf_(2,2)*derxy_(2,vi);

          double temp[3][3];
                          
          /* factor: 
	       
                  /                            \ 
                 |    n+af    /            \    |
                 |   u     , | Dacc o nabla | v |
                 |            \            /    |
                  \                            / 
          */
          temp[0][0]=fac_alphaM_afgdt_tauM_velintaf[0]*derxy_(0,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*aux[0]-fac_afgdt_conv_subaf_vi;
          temp[0][1]=fac_alphaM_afgdt_tauM_velintaf[0]*derxy_(1,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*aux[1];
          temp[0][2]=fac_alphaM_afgdt_tauM_velintaf[0]*derxy_(2,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*aux[2];
                                                                      
          temp[1][0]=fac_alphaM_afgdt_tauM_velintaf[1]*derxy_(0,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*aux[0];
          temp[1][1]=fac_alphaM_afgdt_tauM_velintaf[1]*derxy_(1,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*aux[1]-fac_afgdt_conv_subaf_vi;
          temp[1][2]=fac_alphaM_afgdt_tauM_velintaf[1]*derxy_(2,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*aux[2];
                                                                      
          temp[2][0]=fac_alphaM_afgdt_tauM_velintaf[2]*derxy_(0,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*aux[0];
          temp[2][1]=fac_alphaM_afgdt_tauM_velintaf[2]*derxy_(1,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*aux[1];
          temp[2][2]=fac_alphaM_afgdt_tauM_velintaf[2]*derxy_(2,vi)+fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*aux[2]-fac_afgdt_conv_subaf_vi;

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;

	    elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
	    elemat(fvi  ,fuip ) += temp[0][1]*funct_(ui);
	    elemat(fvi  ,fuipp) += temp[0][2]*funct_(ui);
	                                     
	    elemat(fvip ,fui  ) += temp[1][0]*funct_(ui);
	    elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);
	    elemat(fvip ,fuipp) += temp[1][2]*funct_(ui);
	                                     
	    elemat(fvipp,fui  ) += temp[2][0]*funct_(ui);
	    elemat(fvipp,fuip ) += temp[2][1]*funct_(ui);
	    elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);

	  } // ui
	} // vi
	
        double fac_gdt_afgdt_tauM_facMtau_velintaf[3];
        fac_gdt_afgdt_tauM_facMtau_velintaf[0]=fac_gdt_afgdt_tauM_facMtau*velintaf_(0);
        fac_gdt_afgdt_tauM_facMtau_velintaf[1]=fac_gdt_afgdt_tauM_facMtau*velintaf_(1);
        fac_gdt_afgdt_tauM_facMtau_velintaf[2]=fac_gdt_afgdt_tauM_facMtau*velintaf_(2);

	for (int ui=0; ui<iel; ++ui) // loop columns 
	{
	  const int fuippp =4*ui+3;
          
	  for (int vi=0; vi<iel; ++vi)  // loop rows 
	  {
	    const int fvi   =4*vi;
	    const int fvip  =fvi+1;
	    const int fvipp =fvi+2;
            
            /* factor: tauM, rescaled by gamma*dt
	       
                         /                                      \
                        |    n+af    / /          \         \    |
                        |   u     , | |  nabla Dp  | o nabla | v |
                        |            \ \          /         /    |
                         \                                      /
            */
	    const double aux=derxy_(0,vi)*derxy_(0,ui)+derxy_(1,vi)*derxy_(1,ui)+derxy_(2,vi)*derxy_(2,ui);
	    
	    elemat(fvi  ,fuippp) += fac_gdt_afgdt_tauM_facMtau_velintaf[0]*aux;
	    elemat(fvip ,fuippp) += fac_gdt_afgdt_tauM_facMtau_velintaf[1]*aux;
	    elemat(fvipp,fuippp) += fac_gdt_afgdt_tauM_facMtau_velintaf[2]*aux;
	  } // vi
	} // ui


        for (int vi=0; vi<iel; ++vi)  // loop rows
        {
          const int fvi   =4*vi;
          const int fvip  =fvi+1;
          const int fvipp =fvi+2;
          
          /* factor: tauM*afgdt
	       
                  /                                                   \
                 |    n+af    / / /  n+af       \       \         \    |
                 |   u     , | | |  u    o nabla | Dacc  | o nabla | v |
                 |            \ \ \             /       /         /    |
                  \                                                   /
          */
          double temp[3][3];

          temp[0][0]=fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*derxy_(0,vi);
          temp[0][1]=fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*derxy_(1,vi);
          temp[0][2]=fac_afgdt_afgdt_tauM_facMtau_velintaf[0]*derxy_(2,vi);
                                                                                                        
          temp[1][0]=fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*derxy_(0,vi);
          temp[1][1]=fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*derxy_(1,vi);
          temp[1][2]=fac_afgdt_afgdt_tauM_facMtau_velintaf[1]*derxy_(2,vi);
                                                                                                        
          temp[2][0]=fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*derxy_(0,vi);
          temp[2][1]=fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*derxy_(1,vi);
          temp[2][2]=fac_afgdt_afgdt_tauM_facMtau_velintaf[2]*derxy_(2,vi);

          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui   =4*ui;
            const int fuip  =fui+1;
            const int fuipp =fui+2;
	    
	    elemat(fvi  ,fui  ) += temp[0][0]*conv_c_af_(ui);
	    elemat(fvi  ,fuip ) += temp[0][1]*conv_c_af_(ui);
	    elemat(fvi  ,fuipp) += temp[0][2]*conv_c_af_(ui);
                                             
	    elemat(fvip ,fui  ) += temp[1][0]*conv_c_af_(ui);
	    elemat(fvip ,fuip ) += temp[1][1]*conv_c_af_(ui);
	    elemat(fvip ,fuipp) += temp[1][2]*conv_c_af_(ui);
                                             
	    elemat(fvipp,fui  ) += temp[2][0]*conv_c_af_(ui);
	    elemat(fvipp,fuip ) += temp[2][1]*conv_c_af_(ui);
	    elemat(fvipp,fuipp) += temp[2][2]*conv_c_af_(ui);
	  } // ui
	} // vi

	if (higher_order_ele && newton!=Fluid3::minimal)
	{
	  const double fac_nu_afgdt_afgdt_tauM_facMtau=fac*visceff*afgdt*afgdt*tauM*facMtau;
          
	  double temp[3];

	  temp[0]=fac_nu_afgdt_afgdt_tauM_facMtau*velintaf_(0);
	  temp[1]=fac_nu_afgdt_afgdt_tauM_facMtau*velintaf_(1);
	  temp[2]=fac_nu_afgdt_afgdt_tauM_facMtau*velintaf_(2);


	  for (int vi=0; vi<iel; ++vi)  // loop rows
	  {
	    const int fvi   =4*vi;
	    const int fvip  =fvi+1;
	    const int fvipp =fvi+2;
              
	    double rowtemp[3][3];
            
	    rowtemp[0][0]=temp[0]*derxy_(0,vi);
	    rowtemp[0][1]=temp[0]*derxy_(1,vi);
	    rowtemp[0][2]=temp[0]*derxy_(2,vi);
              
	    rowtemp[1][0]=temp[1]*derxy_(0,vi);
	    rowtemp[1][1]=temp[1]*derxy_(1,vi);
	    rowtemp[1][2]=temp[1]*derxy_(2,vi);
            
	    rowtemp[2][0]=temp[2]*derxy_(0,vi);
	    rowtemp[2][1]=temp[2]*derxy_(1,vi);
	    rowtemp[2][2]=temp[2]*derxy_(2,vi);

	    for (int ui=0; ui<iel; ++ui) // loop columns
	    {
	      const int fui   =4*ui;
	      const int fuip  =fui+1;
	      const int fuipp =fui+2;
                
	      /* 
		 factor: 2.0 * visc * alphaF * gamma * dt * tauM 

                    /                                                \
                   |   n+af   / /             /    \  \         \     |
                 - |  u    , | | nabla o eps | Dacc |  | o nabla | v  |
                   |          \ \             \    /  /         /     |
                    \                                                /      
	      */
	      
	      elemat(fvi  ,fui  ) -= viscs2_(0,ui)*rowtemp[0][0]+derxy2_(3,ui)*rowtemp[0][1]+derxy2_(4,ui)*rowtemp[0][2];
	      elemat(fvi  ,fuip ) -= derxy2_(3,ui)*rowtemp[0][0]+viscs2_(1,ui)*rowtemp[0][1]+derxy2_(5,ui)*rowtemp[0][2];
	      elemat(fvi  ,fuipp) -= derxy2_(4,ui)*rowtemp[0][0]+derxy2_(5,ui)*rowtemp[0][1]+viscs2_(2,ui)*rowtemp[0][2];
                
	      elemat(fvip ,fui  ) -= viscs2_(0,ui)*rowtemp[1][0]+derxy2_(3,ui)*rowtemp[1][1]+derxy2_(4,ui)*rowtemp[1][2];
	      elemat(fvip ,fuip ) -= derxy2_(3,ui)*rowtemp[1][0]+viscs2_(1,ui)*rowtemp[1][1]+derxy2_(5,ui)*rowtemp[1][2];
	      elemat(fvip ,fuipp) -= derxy2_(4,ui)*rowtemp[1][0]+derxy2_(5,ui)*rowtemp[1][1]+viscs2_(2,ui)*rowtemp[1][2];
                
	      elemat(fvipp,fui  ) -= viscs2_(0,ui)*rowtemp[2][0]+derxy2_(3,ui)*rowtemp[2][1]+derxy2_(4,ui)*rowtemp[2][2];
	      elemat(fvipp,fuip ) -= derxy2_(3,ui)*rowtemp[2][0]+viscs2_(1,ui)*rowtemp[2][1]+derxy2_(5,ui)*rowtemp[2][2];
	      elemat(fvipp,fuipp) -= derxy2_(4,ui)*rowtemp[2][0]+derxy2_(5,ui)*rowtemp[2][1]+viscs2_(2,ui)*rowtemp[2][2];
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
      if(pspg == Fluid3::pstab_use_pspg)
      {
        const double fac_afgdt_gamma_dt_tauM_facMtau  = fac*afgdt*gamma*dt*tauM*facMtau;
        const double fac_gdt_gdt_tauM_facMtau         = fac*gamma*dt*tauM*facMtau*gamma*dt;
        const double fac_alphaM_gamma_dt_tauM_facMtau = fac*alphaM*gamma*dt*tauM*facMtau;

        if(higher_order_ele  && newton!=Fluid3::minimal)
        {
          const double fac_visceff_afgdt_gamma_dt_tauM_facMtau
            =
            fac*visceff*afgdt*gamma*dt*tauM*facMtau;

          for (int ui=0; ui<iel; ++ui) // loop columns 
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);


            const double pspg_diffusion_inertia_convect_0_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(0,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_1_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(1,ui)-inertia_and_conv_ui;
            const double pspg_diffusion_inertia_convect_2_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*viscs2_(2,ui)-inertia_and_conv_ui;

            const double scaled_derxy2_3_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(3,ui);
            const double scaled_derxy2_4_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(4,ui);
            const double scaled_derxy2_5_ui=fac_visceff_afgdt_gamma_dt_tauM_facMtau*derxy2_(5,ui);
                 
            for (int vi=0; vi<iel; ++vi)  // loop rows 
            {
              const int fvippp =4*vi+3;
	      
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

              elemat(fvippp,fui  ) -= 
                derxy_(0,vi)*pspg_diffusion_inertia_convect_0_ui
                +
                derxy_(1,vi)*scaled_derxy2_3_ui
                +
                derxy_(2,vi)*scaled_derxy2_4_ui;
              elemat(fvippp,fuip ) -= 
                derxy_(0,vi)*scaled_derxy2_3_ui
                +
                derxy_(1,vi)*pspg_diffusion_inertia_convect_1_ui
                +
                derxy_(2,vi)*scaled_derxy2_5_ui;
              elemat(fvippp,fuipp) -= 
                derxy_(0,vi)*scaled_derxy2_4_ui
                +
                derxy_(1,vi)*scaled_derxy2_5_ui
                +
                derxy_(2,vi)*pspg_diffusion_inertia_convect_2_ui;
            }
          }
        }
        else
        {
          for (int ui=0; ui<iel; ++ui) // loop columns
          {
            const int fui  =4*ui;
            const int fuip =fui+1;
            const int fuipp=fui+2;

            const double inertia_and_conv_ui
              =
              fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
              +
              fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi) // loop rows
            {
              const int fvippp =4*vi+3;
	      
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

              elemat(fvippp,fui  ) +=derxy_(0,vi)*inertia_and_conv_ui;
              elemat(fvippp,fuip ) +=derxy_(1,vi)*inertia_and_conv_ui;
              elemat(fvippp,fuipp) +=derxy_(2,vi)*inertia_and_conv_ui;
            }
          }
        } // neglect viscous linearisations, do just inertia and convective

        for (int ui=0; ui<iel; ++ui) // loop columns 
        {
          const int fuippp=4*ui+3;
          const double scaled_derxy_0=fac_gdt_gdt_tauM_facMtau*derxy_(0,ui);
          const double scaled_derxy_1=fac_gdt_gdt_tauM_facMtau*derxy_(1,ui);
          const double scaled_derxy_2=fac_gdt_gdt_tauM_facMtau*derxy_(2,ui);

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

            elemat(vi*4+3,fuippp) += 
              (scaled_derxy_0*derxy_(0,vi)
               +
               scaled_derxy_1*derxy_(1,vi)
               +
               scaled_derxy_2*derxy_(2,vi)) ;

          } // end loop rows (test functions for matrix)
        } // end loop columns (solution for matrix, test function for vector)

        if (newton==Fluid3::Newton) // if pspg and newton
        {

          for (int vi=0; vi<iel; ++vi) // loop columns
          {
            const int fvippp=4*vi+3;

            const double a=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0)+derxy_(2,vi)*vderxyaf_(2,0));
            const double b=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1)+derxy_(2,vi)*vderxyaf_(2,1));
            const double c=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,2)+derxy_(1,vi)*vderxyaf_(1,2)+derxy_(2,vi)*vderxyaf_(2,2));

            for (int ui=0; ui<iel; ++ui)  // loop rows 
            {
              const int fui=4*ui;
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

              elemat(fvippp,fui  ) += a*funct_(ui);
              elemat(fvippp,fui+1) += b*funct_(ui);
              elemat(fvippp,fui+2) += c*funct_(ui);
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
        if(vstab == Fluid3::viscous_stab_usfem || vstab == Fluid3::viscous_stab_gls)
        {
          const double tauMqs = afgdt*tauM*facMtau;

          const double fac_visc_tauMqs_alphaM        = vstabfac*fac*visc*tauMqs*alphaM;
          const double fac_visc_tauMqs_afgdt         = vstabfac*fac*visc*tauMqs*afgdt;
          const double fac_visc_tauMqs_afgdt_visceff = vstabfac*fac*visc*tauMqs*afgdt*visceff;
          const double fac_visc_tauMqs_gamma_dt      = vstabfac*fac*visc*tauMqs*gamma*dt;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;

            const double inertia_and_conv
              =
              fac_visc_tauMqs_alphaM*funct_(ui)+fac_visc_tauMqs_afgdt*conv_c_af_(ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;
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


              const double a = inertia_and_conv*derxy2_(3,vi);
              const double b = inertia_and_conv*derxy2_(4,vi);
              const double c = inertia_and_conv*derxy2_(5,vi);

              elemat(fvi  ,fui  ) += inertia_and_conv*viscs2_(0,vi);
              elemat(fvi  ,fuip ) += a;
              elemat(fvi  ,fuipp) += b;
              elemat(fvip ,fui  ) += a;
              elemat(fvip ,fuip ) += inertia_and_conv*viscs2_(1,vi);
              elemat(fvip ,fuipp) += c;
              elemat(fvipp,fui  ) += b;
              elemat(fvipp,fuip ) += c;
              elemat(fvipp,fuipp) += inertia_and_conv*viscs2_(2,vi);
            }
          }

          for (int ui=0;ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;

            for (int vi=0;vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;

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
                               (viscs2_(0,vi)*derxy2_(3,ui)
                                +
                                derxy2_(3,vi)*viscs2_(1,ui)
                                +
                                derxy2_(4,vi)*derxy2_(5,ui));
              
              elemat(fvi  ,fuip ) -= a; 
              elemat(fuip ,fvi  ) -= a;
		
              const double b = fac_visc_tauMqs_afgdt_visceff*
                               (viscs2_(0,ui)*derxy2_(4,vi)
                                +
                                derxy2_(3,ui)*derxy2_(5,vi)
                                +
                                derxy2_(4,ui)*viscs2_(2,vi));
              
              elemat(fvipp,fui  ) -= b;
              elemat(fui  ,fvipp) -= b;
              
              const double c = fac_visc_tauMqs_afgdt_visceff*
                               (derxy2_(3,ui)*derxy2_(4,vi)
                                +
                                viscs2_(1,ui)*derxy2_(5,vi)
                                +
                                derxy2_(5,ui)*viscs2_(2,vi));
              
              elemat(fvipp,fuip ) -= c;
              elemat(fuip ,fvipp) -= c;

              elemat(fvi   ,fui ) -= fac_visc_tauMqs_afgdt_visceff*
		                     (viscs2_(0,ui)*viscs2_(0,vi)
                                      +
                                      derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      derxy2_(4,ui)*derxy2_(4,vi));

              elemat(fvip ,fuip ) -= fac_visc_tauMqs_afgdt_visceff*
                                     (derxy2_(3,ui)*derxy2_(3,vi)
                                      +
                                      viscs2_(1,ui)*viscs2_(1,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi));

              elemat(fvipp,fuipp) -= fac_visc_tauMqs_afgdt_visceff*
                                     (derxy2_(4,ui)*derxy2_(4,vi)
                                      +
                                      derxy2_(5,ui)*derxy2_(5,vi)
                                      +
                                      viscs2_(2,ui)*viscs2_(2,vi));
              
            }
          }

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
            const int fui    =4*ui;
            const int fuip   =fui+1;
            const int fuipp  =fuip+1;
            const int fuippp =fuipp+1;
            
            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              const int fvi   =4*vi;
              const int fvip  =fvi+1;
              const int fvipp =fvip+1;

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
              elemat(fvi  ,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*viscs2_(0,vi)
                                       +
                                       derxy_(1,ui)*derxy2_(3,vi)
                                       +
                                       derxy_(2,ui)*derxy2_(4,vi)) ;
              elemat(fvip ,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(3,vi)
                                       +
                                       derxy_(1,ui)*viscs2_(1,vi)
                                       +
                                       derxy_(2,ui)*derxy2_(5,vi)) ;
              elemat(fvipp,fuippp) += fac_visc_tauMqs_gamma_dt*
                                      (derxy_(0,ui)*derxy2_(4,vi)
                                       +
                                       derxy_(1,ui)*derxy2_(5,vi)
                                       +
                                       derxy_(2,ui)*viscs2_(2,vi)) ;
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

          if (newton==Fluid3::Newton)
          {
            
            double temp[3][3];
            for (int vi=0; vi<iel; ++vi) // loop columns (solution for matrix, test function for vector)
            {
              const int fvi    =4*vi;
              const int fvip   =fvi+1;
              const int fvipp  =fvip+1;

              temp[0][0]=(viscs2_(0,vi)*vderxyaf_(0,0)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,0)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][0]=(viscs2_(0,vi)*vderxyaf_(0,1)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,1)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][0]=(viscs2_(0,vi)*vderxyaf_(0,2)
			  +
                          derxy2_(3,vi)*vderxyaf_(1,2)
                          +
                          derxy2_(4,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;
              temp[0][1]=(derxy2_(3,vi)*vderxyaf_(0,0)
			  +
                          viscs2_(1,vi)*vderxyaf_(1,0)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][1]=(derxy2_(3,vi)*vderxyaf_(0,1)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,1)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][1]=(derxy2_(3,vi)*vderxyaf_(0,2)
                          +
                          viscs2_(1,vi)*vderxyaf_(1,2)
                          +
                          derxy2_(5,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;
              temp[0][2]=(derxy2_(4,vi)*vderxyaf_(0,0)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,0)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,0))*fac_visc_tauMqs_afgdt;
              temp[1][2]=(derxy2_(4,vi)*vderxyaf_(0,1)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,1)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,1))*fac_visc_tauMqs_afgdt;
              temp[2][2]=(derxy2_(4,vi)*vderxyaf_(0,2)
                          +
                          derxy2_(5,vi)*vderxyaf_(1,2)
                          +
                          viscs2_(2,vi)*vderxyaf_(2,2))*fac_visc_tauMqs_afgdt;

              for (int ui=0; ui<iel; ++ui)  // loop rows (test functions for matrix)
              {
                const int fui    =4*ui;
                const int fuip   =fui+1;
                const int fuipp  =fuip+1;

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
                elemat(fvi  ,fui  ) += temp[0][0]*funct_(ui);
                elemat(fvi  ,fuip ) += temp[1][0]*funct_(ui);                                     
                elemat(fvi  ,fuipp) += temp[2][0]*funct_(ui);                                     
                elemat(fvip ,fui  ) += temp[0][1]*funct_(ui);                                       
                elemat(fvip ,fuip ) += temp[1][1]*funct_(ui);                                   
                elemat(fvip ,fuipp) += temp[2][1]*funct_(ui);                                   
                elemat(fvipp,fui  ) += temp[0][2]*funct_(ui);                                   
                elemat(fvipp,fuip ) += temp[1][2]*funct_(ui);		                         
                elemat(fvipp,fuipp) += temp[2][2]*funct_(ui);
                                         
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
    if(inertia == Fluid3::inertia_stab_keep
       ||
       inertia == Fluid3::inertia_stab_keep_complete)
    {
      
      double aux_x =(-svelaf_(0)/tauM-pderxynp_(0)-convaf_old_(0)) ;
      double aux_y =(-svelaf_(1)/tauM-pderxynp_(1)-convaf_old_(1)) ;
      double aux_z =(-svelaf_(2)/tauM-pderxynp_(2)-convaf_old_(2)) ;
      
      if(higher_order_ele)
      {
        const double fact =visceff;
        
        aux_x += fact*viscaf_old_(0);
        aux_y += fact*viscaf_old_(1);
        aux_z += fact*viscaf_old_(2);
      }
      
      const double fac_sacc_plus_resM_not_partially_integrated_x =fac*aux_x;
      const double fac_sacc_plus_resM_not_partially_integrated_y =fac*aux_y;
      const double fac_sacc_plus_resM_not_partially_integrated_z =fac*aux_z;
      
      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;
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
        
        elevec(fui  ) -= fac_sacc_plus_resM_not_partially_integrated_x*funct_(ui) ;
        elevec(fui+1) -= fac_sacc_plus_resM_not_partially_integrated_y*funct_(ui) ;
        elevec(fui+2) -= fac_sacc_plus_resM_not_partially_integrated_z*funct_(ui) ;
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
      
      const double fac_inertia_dead_load_z
        =
        fac*(accintam_(2)-bodyforceaf_(2));
      
      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        const int fui=4*ui;
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
        
        elevec(fui  ) -= funct_(ui)*fac_inertia_dead_load_x;
        elevec(fui+1) -= funct_(ui)*fac_inertia_dead_load_y;
        elevec(fui+2) -= funct_(ui)*fac_inertia_dead_load_z;
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
      const int fui =4*ui;
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
    
      elevec(fui  ) += fac_prenp_*derxy_(0,ui) ;
      elevec(fui+1) += fac_prenp_*derxy_(1,ui) ;
      elevec(fui+2) += fac_prenp_*derxy_(2,ui) ;
    }

    const double visceff_fac=visceff*fac;
      
    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
    {
      const int fui=4*ui;
          
      /* viscous term */
      
      /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
      */
          
      elevec(fui   ) -= visceff_fac*
                        (derxy_(0,ui)*vderxyaf_(0,0)*2.0
                         +
                         derxy_(1,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                         +
                         derxy_(2,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0)));
      elevec(fui+1) -= visceff_fac*
	               (derxy_(0,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
                        +
                        derxy_(1,ui)*vderxyaf_(1,1)*2.0
                        +
                        derxy_(2,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1)));
      elevec(fui+2) -= visceff_fac*
	               (derxy_(0,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0))
                        +
                        derxy_(1,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1))
                        +
                        derxy_(2,ui)*vderxyaf_(2,2)*2.0);
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
      
      elevec(ui*4 + 3) -= fac_divunp*funct_(ui);
    } // end loop rows (solution for matrix, test function for vector)

    /*
                /                             \
               |  / n+af       \    n+af       |
              +| | u    o nabla |  u      , v  |
               |  \ G          /               |
                \                             /
    */

    double fac_gridconv[3];
    fac_gridconv[0] = -fac*convu_G_af_old_(0);
    fac_gridconv[1] = -fac*convu_G_af_old_(1);
    fac_gridconv[2] = -fac*convu_G_af_old_(2);

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
    double conv_and_cross_and_re[9];
	
    if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
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
      conv_and_cross_and_re[2]=-velintaf_(0)*fac*(svelaf_(2)+velintaf_(2));
      conv_and_cross_and_re[3]=-velintaf_(1)*fac*(svelaf_(0)+velintaf_(0));
      conv_and_cross_and_re[4]=-velintaf_(1)*fac*(svelaf_(1)+velintaf_(1));
      conv_and_cross_and_re[5]=-velintaf_(1)*fac*(svelaf_(2)+velintaf_(2));
      conv_and_cross_and_re[6]=-velintaf_(2)*fac*(svelaf_(0)+velintaf_(0));
      conv_and_cross_and_re[7]=-velintaf_(2)*fac*(svelaf_(1)+velintaf_(1));
      conv_and_cross_and_re[8]=-velintaf_(2)*fac*(svelaf_(2)+velintaf_(2));
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
      conv_and_cross_and_re[2]=-velintaf_(0)*velintaf_(2)*fac;
      conv_and_cross_and_re[3]=-velintaf_(1)*velintaf_(0)*fac;
      conv_and_cross_and_re[4]=-velintaf_(1)*velintaf_(1)*fac;
      conv_and_cross_and_re[5]=-velintaf_(1)*velintaf_(2)*fac;
      conv_and_cross_and_re[6]=-velintaf_(2)*velintaf_(0)*fac;
      conv_and_cross_and_re[7]=-velintaf_(2)*velintaf_(1)*fac;
      conv_and_cross_and_re[8]=-velintaf_(2)*velintaf_(2)*fac;
    }

    if(reynolds != Fluid3::reynolds_stress_stab_none)
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
      conv_and_cross_and_re[2]-=fac*svelaf_(0)*svelaf_(2);
      conv_and_cross_and_re[3]-=fac*svelaf_(0)*svelaf_(1);
      conv_and_cross_and_re[4]-=fac*svelaf_(1)*svelaf_(1);
      conv_and_cross_and_re[5]-=fac*svelaf_(1)*svelaf_(2);
      conv_and_cross_and_re[6]-=fac*svelaf_(0)*svelaf_(2);
      conv_and_cross_and_re[7]-=fac*svelaf_(1)*svelaf_(2);
      conv_and_cross_and_re[8]-=fac*svelaf_(2)*svelaf_(2);
    }

    for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
    {
      int fui=4*ui;
      /* gridconv with funct                                */
      /* conv, cross, reynolds with derxy                   */

      elevec(fui++) -= 
	fac_gridconv[0]*funct_(ui)
	+
	derxy_(0,ui)*conv_and_cross_and_re[0]
	+
	derxy_(1,ui)*conv_and_cross_and_re[1]
	+
	derxy_(2,ui)*conv_and_cross_and_re[2];
      elevec(fui++) -= 
	fac_gridconv[1]*funct_(ui)
	+
	derxy_(0,ui)*conv_and_cross_and_re[3]
	+
	derxy_(1,ui)*conv_and_cross_and_re[4]
	+
	derxy_(2,ui)*conv_and_cross_and_re[5];
      elevec(fui++) -= 
	fac_gridconv[2]*funct_(ui)
	+
	derxy_(0,ui)*conv_and_cross_and_re[6]
	+
	derxy_(1,ui)*conv_and_cross_and_re[7]
	+
	derxy_(2,ui)*conv_and_cross_and_re[8];
    }

    if(supg == Fluid3::convective_stab_supg)
    {
      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
	int fui=4*ui;
	
	const double fac_conv_c_af_ui = fac*conv_c_af_(ui);
	/*
	  SUPG stabilisation


                  /                             \
                 |  ~n+af    / n+af        \     |
               - |  u     , | c     o nabla | v  |
                 |           \             /     |
                  \                             /
	*/
	
	elevec(fui++) += fac_conv_c_af_ui*svelaf_(0);
	elevec(fui++) += fac_conv_c_af_ui*svelaf_(1);
	elevec(fui  ) += fac_conv_c_af_ui*svelaf_(2);
	
      } // end loop rows
    } // end supg

    //---------------------------------------------------------------
    //
    //        STABILISATION PART, TIME-DEPENDENT SUBGRID-SCALES
    //                    PRESSURE STABILISATION
    //
    //---------------------------------------------------------------
    if(pspg == Fluid3::pstab_use_pspg)
    {
      
      const double fac_svelnpx                      = fac*svelnp(0,iquad);
      const double fac_svelnpy                      = fac*svelnp(1,iquad);
      const double fac_svelnpz                      = fac*svelnp(2,iquad);

      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
      {
        /* factor: -1

                       /                 \
                      |  ~n+1             |
                      |  u    , nabla  q  |
                      |   (i)             |
                       \                 /
        */
        
        elevec(ui*4 + 3) += fac_svelnpx*derxy_(0,ui)+fac_svelnpy*derxy_(1,ui)+fac_svelnpz*derxy_(2,ui);
        
      } // end loop rows (solution for matrix, test function for vector)
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
        const double fac_visc_svelaf_x = vstabfac*fac*visc*svelaf_(0);
        const double fac_visc_svelaf_y = vstabfac*fac*visc*svelaf_(1);
        const double fac_visc_svelaf_z = vstabfac*fac*visc*svelaf_(2);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
          const int fui=4*ui;
          /*
                 /                        \
                |  ~n+af                   |
                |  u      , 2*div eps (v)  |
                |                          |
                 \                        /

          */
          elevec(fui  ) += fac_visc_svelaf_x*viscs2_(0,ui)
	                   +
                           fac_visc_svelaf_y*derxy2_(3,ui)
                           +
                           fac_visc_svelaf_z*derxy2_(4,ui);

          elevec(fui+1) += fac_visc_svelaf_x*derxy2_(3,ui)
                           +
                           fac_visc_svelaf_y*viscs2_(1,ui)
                           +
                           fac_visc_svelaf_z*derxy2_(5,ui);

          elevec(fui+2) += fac_visc_svelaf_x*derxy2_(4,ui)
                           +
                           fac_visc_svelaf_y*derxy2_(5,ui)
                           +
                           fac_visc_svelaf_z*viscs2_(2,ui);

        } // end loop rows (solution for matrix, test function for vector)
      } // endif (a)gls
    }// end if higher order ele

    if(fssgv != Fluid3::fssgv_no)
    {
      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)
      
      for (int ui=0; ui<iel; ++ui)
      {
	/* fine-scale subgrid-viscosity term on right hand side */
	/*
                                  /                              \
                         n+af    |       /    n+af\         / \   |
             - nu_art(fsu    ) * |  eps | Dfsu     | , eps | v |  |
                                 |       \        /         \ /   |
                                  \                              /
	*/
	elevec(ui*4    ) -= vartfac*( 2.0*derxy_(0, ui)*fsvderxyaf_(0, 0)
				     +    derxy_(1, ui)*fsvderxyaf_(0, 1)
				     +    derxy_(1, ui)*fsvderxyaf_(1, 0)
				     +    derxy_(2, ui)*fsvderxyaf_(0, 2)
				     +    derxy_(2, ui)*fsvderxyaf_(2, 0));
	elevec(ui*4 + 1) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 1)
				     +    derxy_(0, ui)*fsvderxyaf_(1, 0)
				     +2.0*derxy_(1, ui)*fsvderxyaf_(1, 1)
				     +    derxy_(2, ui)*fsvderxyaf_(1, 2)
				     +    derxy_(2, ui)*fsvderxyaf_(2, 1));
        elevec(ui*4 + 2) -= vartfac*(     derxy_(0, ui)*fsvderxyaf_(0, 2)
				     +    derxy_(0, ui)*fsvderxyaf_(2, 0)
				     +    derxy_(1, ui)*fsvderxyaf_(1, 2)
				     +    derxy_(1, ui)*fsvderxyaf_(2, 1)
				     +2.0*derxy_(2, ui)*fsvderxyaf_(2, 2));
      } // end loop ui
    } // end not fssgv_no
  } // end loop iquad
  return;
} // end Sysmat_cons_td


/*----------------------------------------------------------------------*
 |  evaluate residuals resM and resC on the element and return          |
 |  averaged values together with averaged subgrid-scale quantities.    |
 |                                                (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalcResAvgs(
  Fluid3*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  RefCountPtr<MAT::Material> mat,
  _MATERIAL*                 material)
{

  // --------------------------------------------------
  // create matrix objects for nodal values
  LINALG::Matrix<iel,1> eprenp    ;
  LINALG::Matrix<3,iel> evelnp    ;
  LINALG::Matrix<3,iel> evelaf    ;
  LINALG::Matrix<3,iel> eaccam    ;
  LINALG::Matrix<3,iel> edispnp   ;
  LINALG::Matrix<3,iel> egridvelaf;
  LINALG::Matrix<3,iel> fsevelaf  ;

  // --------------------------------------------------
  // set parameters for time integration
  ParameterList& timelist = params.sublist("time integration parameters");

  const double alphaF = timelist.get<double>("alpha_F");
  const double gamma  = timelist.get<double>("gamma");
  const double dt     = timelist.get<double>("dt");
  const double time   = timelist.get<double>("time");

  // --------------------------------------------------
  // disable fine-scale subgrid viscosity
  Fluid3::StabilisationAction fssgv = ele->ConvertStringToStabAction("No");

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");
  
  // specify which residual based stabilisation terms
  // will be used
  Fluid3::StabilisationAction tds      = ele->ConvertStringToStabAction(stablist.get<string>("TDS"));
  Fluid3::StabilisationAction inertia  = ele->ConvertStringToStabAction(stablist.get<string>("TRANSIENT"));
  Fluid3::StabilisationAction pspg     = ele->ConvertStringToStabAction(stablist.get<string>("PSPG"));
  Fluid3::StabilisationAction supg     = ele->ConvertStringToStabAction(stablist.get<string>("SUPG"));
  Fluid3::StabilisationAction vstab    = ele->ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  Fluid3::StabilisationAction cstab    = ele->ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  Fluid3::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  Fluid3::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // select tau definition
  Fluid3::TauType whichtau = Fluid3::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid3::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid3::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid3::codina;
    }
    else if(taudef == "Smoothed_FBVW")
    {
      whichtau = Fluid3::smoothed_franca_barrenechea_valentin_wall;
    }
    else if(taudef == "FBVW_without_dt")
    {
      whichtau = Fluid3::fbvw_wo_dt;
    }
    else
    {
      dserror("unknown tau definition\n");
    }
  }

  // flag for higher order elements
  bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());
  
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if(stablist.get<string>("STABTYPE") == "inconsistent")
  {
    higher_order_ele = false;
  }

  // --------------------------------------------------
  // set parameters for turbulence model
  ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

  // the default action is no model
  Fluid3::TurbModelAction turb_mod_action = Fluid3::no_model;

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  double Cs            = 0.0;
  double Cs_delta_sq   = 0.0;
  double l_tau         = 0.0;

  // number of the layer in a turbulent plane channel flow --- used 
  // to compute averaged viscosity etc
  int    nlayer        = 0;

  SetParametersForTurbulenceModel(
    ele            ,
    turbmodelparams,
    fssgv          ,
    turb_mod_action, 
    Cs             ,
    Cs_delta_sq    ,
    l_tau          ,
    nlayer         );

  // --------------------------------------------------
  // extract velocities, pressure and accelerations from the
  // global distributed vectors
  ExtractValuesFromGlobalVectors(
        fssgv         ,
        ele->is_ale_  ,
        discretization,
        lm            ,
        eprenp        ,
        evelnp        ,
        evelaf        ,
        eaccam        ,
        edispnp       ,
        egridvelaf    ,
        fsevelaf
    );

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<blitz::Array<double,1> > myknots(3);

  // for isogeometric elements
  if(ele->Shape()==Fluid3::nurbs8 || ele->Shape()==Fluid3::nurbs27)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());
  }

  // visceff will contain the computed effective viscosity
  double visceff       = 0.0;

  // the coordinates of the element layers in the channel
  RefCountPtr<vector<double> > planecoords  = params.get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);
  
  if(planecoords==Teuchos::null)
    dserror("planecoords is null, but need channel_flow_of_height_2\n");

  // ---------------------------------------------------
  // working arrays for the quantities we want to compute
  LINALG::Matrix<3,1>  mean_res       ;
  LINALG::Matrix<3,1>  mean_sacc      ;
  LINALG::Matrix<3,1>  mean_svelaf    ;
  LINALG::Matrix<3,1>  mean_res_sq    ;
  LINALG::Matrix<3,1>  mean_sacc_sq   ;
  LINALG::Matrix<3,1>  mean_svelaf_sq ;
  LINALG::Matrix<3,1>  mean_tauinvsvel;

  double vol             = 0.0;

  double averaged_tauC   = 0.0;
  double averaged_tauM   = 0.0;

  double abs_res         = 0.0;
  double abs_svel        = 0.0;
  double abs_sacc        = 0.0;

  double mean_resC       = 0.0;
  double mean_resC_sq    = 0.0;
  double mean_sprenp     = 0.0;
  double mean_sprenp_sq  = 0.0;

  double eps_sacc        = 0.0;
  double eps_pspg        = 0.0;
  double eps_supg        = 0.0;
  double eps_cross       = 0.0;
  double eps_rey         = 0.0;
  double eps_cstab       = 0.0;
  double eps_vstab       = 0.0;
  double eps_eddyvisc    = 0.0;
  double eps_visc        = 0.0;
  double eps_conv        = 0.0;

  
  mean_res       .Clear();
  mean_sacc      .Clear();
  mean_svelaf    .Clear();
  mean_res_sq    .Clear();
  mean_sacc_sq   .Clear();
  mean_svelaf_sq .Clear();
  mean_tauinvsvel.Clear();

  //------------------------------------------------------------------
  //           SET TIME INTEGRATION SCHEME RELATED DATA
  //------------------------------------------------------------------

  //         n+alpha_F     n+1
  //        t          = t     - (1-alpha_F) * dt
  //
  const double timealphaF = time-(1-alphaF)*dt;

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

  SetElementData(ele            ,
                 edispnp        ,
                 evelaf         ,
                 fsevelaf       ,
                 myknots        ,
                 timealphaF     ,
                 hk             ,
                 mk             ,
                 material       ,
                 visc           ,
                 fssgv          ,
                 turb_mod_action,
                 l_tau          ,
                 Cs             ,
                 Cs_delta_sq    ,
                 visceff        );

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

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    //--------------------------------------------------------------
    // Get all global shape functions, first and eventually second 
    // derivatives in a gausspoint and integration weight including 
    //                   jacobi-determinant
    //--------------------------------------------------------------

    const double fac=ShapeFunctionsFirstAndSecondDerivatives(
      ele             ,
      iquad           ,
      intpoints       ,
      myknots         ,
      higher_order_ele);
    
    //--------------------------------------------------------------
    //            interpolate nodal values to gausspoint
    //--------------------------------------------------------------
    InterpolateToGausspoint(ele             ,
                            egridvelaf      ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
                            higher_order_ele);

    if(tds == Fluid3::subscales_time_dependent)
    {
      /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

      */
      for (int rr=0;rr<3;++rr)
      {
        svelaf_(rr) = alphaF*svelnp(rr,iquad)+(1.0-alphaF)*sveln(rr,iquad);
      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,tds,gamma,dt,hk,mk,visceff);

    // -----------------------------------------------------
    // Element volume
    vol               += fac;

    // ------------------------------------------------------
    // Element average:
    //     o residuals
    //     o subscale acceleration
    //     o subscale velocity
    //     o subscale pressure

    for(int rr=0;rr<3;++rr)
    {
      mean_res    (rr) += resM_(rr)*fac;
      mean_res_sq (rr) += resM_(rr)*resM_(rr)*fac;
    }
    abs_res    += sqrt(resM_(0)*resM_(0)+resM_(1)*resM_(1)+resM_(2)*resM_(2))*fac;

    if(tds == Fluid3::subscales_time_dependent
       &&
       (inertia == Fluid3::inertia_stab_keep
        ||
        inertia == Fluid3::inertia_stab_keep_complete)) 
    {
      for(int rr=0;rr<3;++rr)
      {
        const double aux = -1.0/tau_(0)*svelaf_(rr) -resM_(rr);

        mean_sacc   (rr) += aux*fac;
        mean_sacc_sq(rr) += aux*aux*fac;
      }
      const double aux0 = -1.0/tau_(0)*svelaf_(0)-resM_(0);
      const double aux1 = -1.0/tau_(0)*svelaf_(1)-resM_(1);
      const double aux2 = -1.0/tau_(0)*svelaf_(2)-resM_(2);

      abs_sacc += sqrt(aux0*aux0+aux1*aux1+aux2*aux2)*fac;
    }

    if(tds == Fluid3::subscales_time_dependent)
    {
      for(int rr=0;rr<3;++rr)
      {
        mean_svelaf   (rr) += svelaf_(rr)*fac;
        mean_svelaf_sq(rr) += svelaf_(rr)*svelaf_(rr)*fac;
      }

      abs_svel +=sqrt(svelaf_(0)*svelaf_(0)
                      +
                      svelaf_(1)*svelaf_(1)
                      +
                      svelaf_(2)*svelaf_(2))*fac;
    }
    else
    {
      for(int rr=0;rr<3;++rr)
      {
        const double aux = tau_(0)*resM_(rr);

        mean_svelaf   (rr) -= aux*fac;
        mean_svelaf_sq(rr) += aux*aux*fac;
      }

      abs_svel +=sqrt(resM_(0)*resM_(0)
                      +
                      resM_(1)*resM_(1)
                      +
                      resM_(2)*resM_(2))*tau_(0)*fac;
    }

    for(int rr=0;rr<3;++rr)
    {
      mean_tauinvsvel(rr)+=mean_svelaf(rr)/tau_(0);
    }
    

    {
      const double aux = tau_(2)*divunp_;

      mean_sprenp     -= aux*fac;
      mean_sprenp_sq  += aux*aux*fac;
    }

    // ------------------------------------------------------
    // Element average:
    //      o tauM
    //      o tauC

    averaged_tauM+=tau_(0)*fac;
    averaged_tauC+=tau_(2)*fac;

    mean_resC    += divunp_*fac;
    mean_resC_sq += divunp_*divunp_*fac;

    // ------------------------------------------------------
    // Element average dissipation/production rates:
    //      o inertia
    //      o pspg
    //      o supg
    //      o viscous
    //      o cross
    //      o reynolds
    //      o cstab
    //      o Smagorinsky
    //      o Galerkin visc
    //      o Galerkin convection

    if(tds      == Fluid3::subscales_time_dependent
       &&
       (inertia  == Fluid3::inertia_stab_keep
        ||
        inertia == Fluid3::inertia_stab_keep_complete))
    {
      double sacc[3];
      sacc[0]= -1.0/tau_(0)*svelaf_(0) -resM_(0);
      sacc[1]= -1.0/tau_(0)*svelaf_(1) -resM_(1);
      sacc[2]= -1.0/tau_(0)*svelaf_(2) -resM_(2);
        
      /*

               /                 \  
              |   ~ n+am    n+af  | 
              |  acc     , u      | 
              |     (i)     (i)   | 
               \                 /  

      */
      eps_sacc+=(sacc[0]*velintaf_(0)+sacc[1]*velintaf_(1)+sacc[2]*velintaf_(2))*fac;
    }

    if(pspg     == Fluid3::pstab_use_pspg)
    {
      // contribution of this gausspoint to energy dissipated by pspg
      /*
                        /                  \
                       |  ^ n+1        n+1  |
                     - |  u   , nabla p     |
                       |                    |
                        \                  /
      */
      if(tds == Fluid3::subscales_time_dependent)
      {
        eps_pspg-=
          (svelnp(0,iquad)*pderxynp_(0)
           +
           svelnp(1,iquad)*pderxynp_(1)
           +
           svelnp(2,iquad)*pderxynp_(2))*fac;
      }
      else
      {
        eps_pspg+=tau_(0)*fac*
          (resM_(0)*pderxynp_(0)
           +
           resM_(1)*pderxynp_(1)
           +
           resM_(2)*pderxynp_(2));
      }
    }

    if(supg     == Fluid3::convective_stab_supg)
    {
      // contribution of this gausspoint to energy dissipated by supg
      /*

                   /                                      \
                  |  ~n+af    / n+af \           / n+af\   |
          - rho * |  u     , | u      | o nabla | u     |  |
                  |           \      /           \     /   |
                   \                                      /
      */
      if(tds == Fluid3::subscales_time_dependent)
      {
        eps_supg-=fac*
          (svelaf_(0)*convaf_old_(0)
           +
           svelaf_(1)*convaf_old_(1)
           +
           svelaf_(2)*convaf_old_(2));
      }
      else
      {
        eps_supg+=tau_(0)*fac*
          (resM_(0)*convaf_old_(0)
           +
           resM_(1)*convaf_old_(1)
           +
           resM_(2)*convaf_old_(2));
      }
    }

    if(cross    != Fluid3::cross_stress_stab_none)
    {
      /*
               /                                       \
              |    n+af    /~n+af \           / n+af\   |
              |   u     , | u      | o nabla | u     |  |
              |            \      /           \     /   |
               \                                       /
        
      */
      if(tds == Fluid3::subscales_time_dependent)
      {
        eps_cross+=fac*(velintaf_(0)*(svelaf_(0)*vderxyaf_(0,0)+
                                      svelaf_(1)*vderxyaf_(0,1)+
                                      svelaf_(2)*vderxyaf_(0,2))
                        +
                        velintaf_(1)*(svelaf_(0)*vderxyaf_(1,0)+
                                      svelaf_(1)*vderxyaf_(1,1)+
                                      svelaf_(2)*vderxyaf_(1,2))
                        +
                        velintaf_(2)*(svelaf_(0)*vderxyaf_(2,0)+
                                      svelaf_(1)*vderxyaf_(2,1)+
                                      svelaf_(2)*vderxyaf_(2,2)));
      }
      else
      {
        eps_cross-=
          tau_(0)*fac*(velintaf_(0)*(resM_(0)*vderxyaf_(0,0)+
                                     resM_(1)*vderxyaf_(0,1)+
                                     resM_(2)*vderxyaf_(0,2))
                       +
                       velintaf_(1)*(resM_(0)*vderxyaf_(1,0)+
                                     resM_(1)*vderxyaf_(1,1)+
                                     resM_(2)*vderxyaf_(1,2))
                       +
                       velintaf_(2)*(resM_(0)*vderxyaf_(2,0)+
                                     resM_(1)*vderxyaf_(2,1)+
                                     resM_(2)*vderxyaf_(2,2)));
      }
    }

    if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
    {
      // contribution of this gausspoint to energy
      // dissipated/produced by reynolds 
      /*

                 /                                 \
                |  ~n+af    /~n+af        \   n+af  |
        - rho * |  u     , | u     o nabla | u      |
                |           \             /         |
                 \                                 /
      */
      if(tds == Fluid3::subscales_time_dependent)
      {
        eps_rey-=fac*
          (svelaf_(0)*(svelaf_(0)*vderxyaf_(0,0)+
                       svelaf_(1)*vderxyaf_(0,1)+
                       svelaf_(2)*vderxyaf_(0,2))
           +
           svelaf_(1)*(svelaf_(0)*vderxyaf_(1,0)+
                       svelaf_(1)*vderxyaf_(1,1)+
                       svelaf_(2)*vderxyaf_(1,2))
           +
           svelaf_(2)*(svelaf_(0)*vderxyaf_(2,0)+
                       svelaf_(1)*vderxyaf_(2,1)+
                       svelaf_(2)*vderxyaf_(2,2)));
      }
      else
      {
        eps_rey-=
          tau_(0)*tau_(0)*fac*
          (resM_(0)*(resM_(0)*vderxyaf_(0,0)+
                     resM_(1)*vderxyaf_(0,1)+
                     resM_(2)*vderxyaf_(0,2))
           +
           resM_(1)*(resM_(0)*vderxyaf_(1,0)+
                     resM_(1)*vderxyaf_(1,1)+
                     resM_(2)*vderxyaf_(1,2))
           +
           resM_(2)*(resM_(0)*vderxyaf_(2,0)+
                     resM_(1)*vderxyaf_(2,1)+
                     resM_(2)*vderxyaf_(2,2)));
      }
    }

    if(vstab    != Fluid3::viscous_stab_none)
    {

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

      /*
         /                        \
        |  ~n+af                   |
        |  u      , 2*div eps (v)  |
        |                          |
         \                        /
        
      */
      if(tds == Fluid3::subscales_time_dependent)
      {
        eps_vstab-=
          vstabfac*fac*
          (svelaf_(0)*visceff*viscaf_old_(0)
           +
           svelaf_(1)*visceff*viscaf_old_(1)
           +
           svelaf_(2)*visceff*viscaf_old_(2));
      }
      else
      {
        eps_vstab+=
          vstabfac*tau_(0)*fac*
          (resM_(0)*visceff*viscaf_old_(0)
           +
           resM_(1)*visceff*viscaf_old_(1)
           +
           resM_(2)*visceff*viscaf_old_(2));
      }
    }

    if(cstab    != Fluid3::continuity_stab_none)
    {
      // contribution of this gausspoint to energy dissipated/produced by pprime
      /*
                        /                             \
                       |           n+1            n+1  |
                 +tauC |  nabla o u    , nabla o u     |
                       |                               |
                        \                             /
      */
      eps_cstab+=tau_(2)*fac*divunp_*divunp_;
    }
    
    {
      // contribution of this gausspoint to convective
      // dissipation/production (Galerkin)  
      /*
                   /                                \
                  |   n+af   / n+af \     /  n+af\   |
                  |  u    , | u      | o | u      |  |
                  |          \      /     \      /   |
                   \                                /
      */
        eps_conv+=fac*
          (velintaf_(0)*(velintaf_(0)*vderxyaf_(0,0)+
                         velintaf_(1)*vderxyaf_(0,1)+
                         velintaf_(2)*vderxyaf_(0,2))
           +
           velintaf_(1)*(velintaf_(0)*vderxyaf_(1,0)+
                         velintaf_(1)*vderxyaf_(1,1)+
                         velintaf_(2)*vderxyaf_(1,2))
           +
           velintaf_(2)*(velintaf_(0)*vderxyaf_(2,0)+
                         velintaf_(1)*vderxyaf_(2,1)+
                         velintaf_(2)*vderxyaf_(2,2)));
    }

    {
      // contribution of this gausspoint to viscous energy
      // dissipation (Galerkin)  
      /*
                   /                                \
                  |       / n+af \         / n+af\   |
        2* visc * |  eps | u      | , eps | u     |  |
                  |       \      /         \     /   |
                   \                                /
      */
	
      double two_eps_uaf[9];
      two_eps_uaf[0]=(vderxyaf_(0,0))*2.0;
      two_eps_uaf[1]=(vderxyaf_(0,1)+vderxyaf_(1,0));
      two_eps_uaf[2]=(vderxyaf_(0,2)+vderxyaf_(2,0));
      two_eps_uaf[3]=(vderxyaf_(0,1)+vderxyaf_(1,0));
      two_eps_uaf[4]=(vderxyaf_(1,1))*2.0;
      two_eps_uaf[5]=(vderxyaf_(1,2)+vderxyaf_(2,1));
      two_eps_uaf[6]=(vderxyaf_(0,2)+vderxyaf_(2,0));
      two_eps_uaf[7]=(vderxyaf_(1,2)+vderxyaf_(2,1));
      two_eps_uaf[8]=(vderxyaf_(2,2))*2.0;

      eps_visc+=
        0.5*visceff*fac*(two_eps_uaf[0]*two_eps_uaf[0]+
                         two_eps_uaf[1]*two_eps_uaf[1]+
                         two_eps_uaf[2]*two_eps_uaf[2]+
                         two_eps_uaf[3]*two_eps_uaf[3]+
                         two_eps_uaf[4]*two_eps_uaf[4]+
                         two_eps_uaf[5]*two_eps_uaf[5]+
                         two_eps_uaf[6]*two_eps_uaf[6]+
                         two_eps_uaf[7]*two_eps_uaf[7]+
                         two_eps_uaf[8]*two_eps_uaf[8]);

      // contribution of this gausspoint to viscous energy
      // dissipation (Smagorinsky)  
      /*
                      /                                \
                     |       / n+af \         / n+af\   |
        2* visc    * |  eps | u      | , eps | u     |  |
               art   |       \      /         \     /   |
                      \                                /
      */
      if(turb_mod_action == Fluid3::dynamic_smagorinsky
         ||
         turb_mod_action == Fluid3::smagorinsky_with_wall_damping
         ||
         turb_mod_action == Fluid3::smagorinsky)
      {
        eps_eddyvisc+=
          0.5*(visceff-visc)*fac*(two_eps_uaf[0]*two_eps_uaf[0]+
                                  two_eps_uaf[1]*two_eps_uaf[1]+
                                  two_eps_uaf[2]*two_eps_uaf[2]+
                                  two_eps_uaf[3]*two_eps_uaf[3]+
                                  two_eps_uaf[4]*two_eps_uaf[4]+
                                  two_eps_uaf[5]*two_eps_uaf[5]+
                                  two_eps_uaf[6]*two_eps_uaf[6]+
                                  two_eps_uaf[7]*two_eps_uaf[7]+
                                  two_eps_uaf[8]*two_eps_uaf[8]);
      }
    }
  } // end integration loop


  for(int rr=0;rr<3;++rr)
  {
    mean_res        (rr)/= vol;
    mean_res_sq     (rr)/= vol;
    mean_sacc       (rr)/= vol;
    mean_sacc_sq    (rr)/= vol;
    mean_svelaf     (rr)/= vol;
    mean_svelaf_sq  (rr)/= vol;
    mean_tauinvsvel (rr)/= vol;
  }

  abs_res         /= vol;
  abs_sacc        /= vol;
  abs_svel        /= vol;

  mean_resC       /= vol;
  mean_resC_sq    /= vol;
  mean_sprenp     /= vol;
  mean_sprenp_sq  /= vol;

  averaged_tauC   /= vol;
  averaged_tauM   /= vol;

  eps_sacc        /= vol;
  eps_pspg        /= vol;
  eps_supg        /= vol;
  eps_cross       /= vol;
  eps_rey         /= vol;
  eps_cstab       /= vol;
  eps_vstab       /= vol;
  eps_eddyvisc    /= vol;
  eps_visc        /= vol;
  eps_conv        /= vol;

  // ---------------------------------------------------
  // the vectors containing the local sums over layers

  RefCountPtr<vector<double> > incrvol           = params.get<RefCountPtr<vector<double> > >("incrvol"          );

  RefCountPtr<vector<double> > incrres           = params.get<RefCountPtr<vector<double> > >("incrres"          );
  RefCountPtr<vector<double> > incrres_sq        = params.get<RefCountPtr<vector<double> > >("incrres_sq"       );
  RefCountPtr<vector<double> > incrabsres        = params.get<RefCountPtr<vector<double> > >("incrabsres"       );
  RefCountPtr<vector<double> > incrtauinvsvel    = params.get<RefCountPtr<vector<double> > >("incrtauinvsvel"   );
 
  RefCountPtr<vector<double> > incrsacc          = params.get<RefCountPtr<vector<double> > >("incrsacc"         );
  RefCountPtr<vector<double> > incrsacc_sq       = params.get<RefCountPtr<vector<double> > >("incrsacc_sq"      );
  RefCountPtr<vector<double> > incrabssacc       = params.get<RefCountPtr<vector<double> > >("incrabssacc"      );

  RefCountPtr<vector<double> > incrsvelaf        = params.get<RefCountPtr<vector<double> > >("incrsvelaf"       );
  RefCountPtr<vector<double> > incrsvelaf_sq     = params.get<RefCountPtr<vector<double> > >("incrsvelaf_sq"    );
  RefCountPtr<vector<double> > incrabssvelaf     = params.get<RefCountPtr<vector<double> > >("incrabssvelaf"    );
  
  RefCountPtr<vector<double> > incrresC          = params.get<RefCountPtr<vector<double> > >("incrresC"         );
  RefCountPtr<vector<double> > incrresC_sq       = params.get<RefCountPtr<vector<double> > >("incrresC_sq"      );
  RefCountPtr<vector<double> > spressnp          = params.get<RefCountPtr<vector<double> > >("incrspressnp"     );
  RefCountPtr<vector<double> > spressnp_sq       = params.get<RefCountPtr<vector<double> > >("incrspressnp_sq"  );

  RefCountPtr<vector<double> > incrtauC          = params.get<RefCountPtr<vector<double> > >("incrtauC"         );
  RefCountPtr<vector<double> > incrtauM          = params.get<RefCountPtr<vector<double> > >("incrtauM"         );

  RefCountPtr<vector<double> > incr_eps_sacc     = params.get<RefCountPtr<vector<double> > >("incr_eps_sacc"    );
  RefCountPtr<vector<double> > incr_eps_pspg     = params.get<RefCountPtr<vector<double> > >("incr_eps_pspg"    );
  RefCountPtr<vector<double> > incr_eps_supg     = params.get<RefCountPtr<vector<double> > >("incr_eps_supg"    );
  RefCountPtr<vector<double> > incr_eps_cross    = params.get<RefCountPtr<vector<double> > >("incr_eps_cross"   );
  RefCountPtr<vector<double> > incr_eps_rey      = params.get<RefCountPtr<vector<double> > >("incr_eps_rey"     );
  RefCountPtr<vector<double> > incr_eps_cstab    = params.get<RefCountPtr<vector<double> > >("incr_eps_cstab"   );
  RefCountPtr<vector<double> > incr_eps_vstab    = params.get<RefCountPtr<vector<double> > >("incr_eps_vstab"   );
  RefCountPtr<vector<double> > incr_eps_eddyvisc = params.get<RefCountPtr<vector<double> > >("incr_eps_eddyvisc");
  RefCountPtr<vector<double> > incr_eps_visc     = params.get<RefCountPtr<vector<double> > >("incr_eps_visc"    );
  RefCountPtr<vector<double> > incr_eps_conv     = params.get<RefCountPtr<vector<double> > >("incr_eps_conv"    );


  //this will be the y-coordinate of a point in the element interior
  double center = 0;

  // get node coordinates of element
  LINALG::Matrix<3,iel>  xyze;
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=ele->Nodes()[inode]->X()[0];
    xyze(1,inode)=ele->Nodes()[inode]->X()[1];
    xyze(2,inode)=ele->Nodes()[inode]->X()[2];
    
    center+=xyze(1,inode);
  }
  center/=iel;
  
  bool found = false;
  
  for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
  {
    if(center<(*planecoords)[nlayer+1])
    {
      found = true;
      break;
    }
    nlayer++;
  }
  if (found ==false)
  {
    dserror("could not determine element layer");
  }

  // collect layer volume
  (*incrvol  )[nlayer] += vol;

  // averages of stabilisation parameters
  (*incrtauC )[nlayer] += averaged_tauC;
  (*incrtauM )[nlayer] += averaged_tauM;
  
  // averages of momentum residuals, subscale velocity and accelerations
  for(int mm=0;mm<3;++mm)
  {
    (*incrres       )[3*nlayer+mm] += mean_res       (mm);
    (*incrres_sq    )[3*nlayer+mm] += mean_res_sq    (mm);
    (*incrsacc      )[3*nlayer+mm] += mean_sacc      (mm);
    (*incrsacc_sq   )[3*nlayer+mm] += mean_sacc_sq   (mm);
    
    (*incrsvelaf    )[3*nlayer+mm] += mean_svelaf    (mm);
    (*incrsvelaf_sq )[3*nlayer+mm] += mean_svelaf_sq (mm);

    (*incrtauinvsvel)[3*nlayer+mm] += mean_tauinvsvel(mm);
  }

  (*incrabsres       )[nlayer] += abs_res;
  (*incrabssacc      )[nlayer] += abs_sacc;
  (*incrabssvelaf    )[nlayer] += abs_svel;

  // averages of subscale pressure and continuity residuals
  (*incrresC         )[nlayer] += mean_resC      ;
  (*incrresC_sq      )[nlayer] += mean_resC_sq   ;

  (*spressnp         )[nlayer] += mean_sprenp    ;
  (*spressnp_sq      )[nlayer] += mean_sprenp_sq ;

  (*incr_eps_sacc    )[nlayer] += eps_sacc       ;
  (*incr_eps_pspg    )[nlayer] += eps_pspg       ;
  (*incr_eps_supg    )[nlayer] += eps_supg       ;
  (*incr_eps_cross   )[nlayer] += eps_cross      ;
  (*incr_eps_rey     )[nlayer] += eps_rey        ;
  (*incr_eps_cstab   )[nlayer] += eps_cstab      ;
  (*incr_eps_vstab   )[nlayer] += eps_vstab      ;
  (*incr_eps_eddyvisc)[nlayer] += eps_eddyvisc   ;
  (*incr_eps_visc    )[nlayer] += eps_visc       ;
  (*incr_eps_conv    )[nlayer] += eps_conv       ;

  return(0);
}

/*----------------------------------------------------------------------*
 |  extract velocities, pressure and accelerations from the global      |
 |  distributed vectors                           (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::ExtractValuesFromGlobalVectors(
        const Fluid3::StabilisationAction          fssgv         ,
        const bool                                 is_ale        ,
        const DRT::Discretization&                 discretization,
        const vector<int>&                         lm            ,
        LINALG::Matrix<iel,1>& eprenp        ,
        LINALG::Matrix<3,iel>& evelnp        ,
        LINALG::Matrix<3,iel>& evelaf        ,
        LINALG::Matrix<3,iel>& eaccam        ,
        LINALG::Matrix<3,iel>& edispnp       ,
        LINALG::Matrix<3,iel>& egridvelaf    ,
        LINALG::Matrix<3,iel>& fsevelaf
        )
{
    // velocity and pressure values (current iterate, n+1)
    RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

    // velocities    (intermediate time step, n+alpha_F)
    RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

    // accelerations (intermediate time step, n+alpha_M)
    RefCountPtr<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");

    if (velnp==null || velaf==null || accam==null)
    {
      dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
    }

    // extract local values from the global vectors
    vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
    
    vector<double> myvelaf(lm.size());
    DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
    
    vector<double> myaccam(lm.size());
    DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);

    // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
    // Additionally only the 'velocity' components of my_velaf
    // and my_accam are important!
    for (int i=0;i<iel;++i)
    {
      int fi    =4*i;
      
      eaccam(0,i) = myaccam[fi  ];
      evelnp(0,i) = myvelnp[fi  ];
      evelaf(0,i) = myvelaf[fi++];
      evelnp(1,i) = myvelnp[fi  ];
      eaccam(1,i) = myaccam[fi  ];
      evelaf(1,i) = myvelaf[fi++];
      eaccam(2,i) = myaccam[fi  ];
      evelnp(2,i) = myvelnp[fi  ];
      evelaf(2,i) = myvelaf[fi++];
      eprenp(i)   = myvelnp[fi  ];
    }
    
    if(is_ale)
    {
      // get most recent displacements
      RefCountPtr<const Epetra_Vector> dispnp
        =
        discretization.GetState("dispnp");
      
      // get intermediate grid velocities
      RefCountPtr<const Epetra_Vector> gridvelaf
        =
        discretization.GetState("gridvelaf");
      
      if (dispnp==null || gridvelaf==null)
      {
        dserror("Cannot get state vectors 'dispnp' and/or 'gridvelaf'");
      }

      vector<double> mydispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      
      vector<double> mygridvelaf(lm.size());
      DRT::UTILS::ExtractMyValues(*gridvelaf,mygridvelaf,lm);

      // extract velocity part from "mygridvelaf" and get
      // set element displacements
      for (int i=0;i<iel;++i)
      {
        int fi    =4*i;
        int fip   =fi+1;
        int fipp  =fip+1;

        egridvelaf(0,i) = mygridvelaf[fi  ];
        egridvelaf(1,i) = mygridvelaf[fip ];
        egridvelaf(2,i) = mygridvelaf[fipp];

        edispnp(0,i)    = mydispnp   [fi  ];
        edispnp(1,i)    = mydispnp   [fip ];
        edispnp(2,i)    = mydispnp   [fipp];
      }
    }

    // get fine-scale velocity
    if (fssgv != Fluid3::fssgv_no)
    {
      RCP<const Epetra_Vector> fsvelaf;

      fsvelaf = discretization.GetState("fsvelaf");
      if (fsvelaf==null) dserror("Cannot get state vector 'fsvelaf'");
      vector<double> myfsvelaf(lm.size());
      DRT::UTILS::ExtractMyValues(*fsvelaf,myfsvelaf,lm);

      // get fine-scale velocity and insert into element arrays
      for (int i=0;i<iel;++i)
      {
        fsevelaf(0,i) = myfsvelaf[0+(i*4)];
        fsevelaf(1,i) = myfsvelaf[1+(i*4)];
        fsevelaf(2,i) = myfsvelaf[2+(i*4)];
      }
    }
    else
    {
      for (int i=0;i<iel;++i)
      {
        fsevelaf(0,i) = 0.0;
        fsevelaf(1,i) = 0.0;
        fsevelaf(2,i) = 0.0;
      }
    }
  return;
} //ExtractValuesFromGlobalVectors


/*----------------------------------------------------------------------*
 |  Set parameters for turbulence models                                |
 |                                                (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::SetParametersForTurbulenceModel(
  const Fluid3*                       ele            ,
  ParameterList                     & turbmodelparams,
  const Fluid3::StabilisationAction & fssgv          ,
  Fluid3::TurbModelAction           & turb_mod_action, 
  double                            & Cs             ,
  double                            & Cs_delta_sq    ,
  double                            & l_tau          ,
  int                               & nlayer         
  )
{

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv != Fluid3::fssgv_no && turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
    dserror("No combination of a classical (all-scale) turbulence model and a fine-scale subgrid-viscosity approach currently possible!");
  if (fssgv != Fluid3::fssgv_no) Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action = Fluid3::smagorinsky;
      Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      RefCountPtr<vector<double> > planecoords      = turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);

      if(planecoords==Teuchos::null)
        dserror("planecoords is null, but need channel_flow_of_height_2\n");

      // for the Smagorinsky model with van Driest damping, we need a viscous length to determine
      // the y+ (heigth in wall units)
      turb_mod_action = Fluid3::smagorinsky_with_wall_damping;
      Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
      l_tau           = turbmodelparams.get<double>("CHANNEL_L_TAU");

      //this will be the y-coordinate of a point in the element interior
      double center = 0;
      for(int inode=0;inode<iel;inode++)
      {
        center+=ele->Nodes()[inode]->X()[1];
      }
      center/=iel;

      bool found = false;
      for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
      {
        if(center<(*planecoords)[nlayer+1])
        {
          found = true;
          break;
        }
        nlayer++;
      }
      if (found ==false)
      {
        dserror("could not determine element layer");
      }
    }
    else if (physical_turbulence_model == "Dynamic_Smagorinsky")
    {
      turb_mod_action = Fluid3::dynamic_smagorinsky;

      if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
          ==
          "channel_flow_of_height_2")
      {
        RCP<vector<double> > averaged_LijMij
          =
          turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
        RCP<vector<double> > averaged_MijMij
          =
          turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");

        RCP<vector<double> > planecoords
          =
          turbmodelparams.get<RCP<vector<double> > >("planecoords_");

        //this will be the y-coordinate of a point in the element interior
        double center = 0;
        for(int inode=0;inode<iel;inode++)
        {
          center+=ele->Nodes()[inode]->X()[1];
        }
        center/=iel;

        bool found = false;
        for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
        {
          if(center<(*planecoords)[nlayer+1])
          {
            found = true;
            break;
          }
          nlayer++;
        }
        if (found ==false)
        {
          dserror("could not determine element layer");
        }

        Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;

        // clipping to get algorithm stable
        if (Cs_delta_sq<0)
        {
          Cs_delta_sq=0;
        }
      }
      else
      {
        Cs_delta_sq = ele->Cs_delta_sq_;
      }
    }
    else
    {
      dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
    }
  }

  return;
} // end SetParametersForTurbulenceModel


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::GetNodalBodyForce(
  const Fluid3* ele, 
  const double  time)
{
  constant_bodyforce_=false;

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
    for(int isd=0;isd<3;isd++)
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
 | calculates material viscosity                            u.may 05/08 |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalVisc(
  const struct _MATERIAL*                 material   ,
  double&                           	  visc       ,
  const double &                          rateofshear
  )
{
  if(material->mattyp == m_carreauyasuda)
  {
    double nu_0 = material->m.carreauyasuda->nu_0;      // parameter for zero-shear viscosity
    double nu_inf = material->m.carreauyasuda->nu_inf;  // parameter for infinite-shear viscosity
    double lambda = material->m.carreauyasuda->lambda;  // parameter for characteristic time
    double a = material->m.carreauyasuda->a_param;      // constant parameter
    double b = material->m.carreauyasuda->b_param;  	// constant parameter

    // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    const double tmp = pow(lambda*rateofshear,b);
    visc = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
  }
  else if(material->mattyp == m_modpowerlaw)
  {
    // get material parameters
    double m  	  = material->m.modpowerlaw->m_cons;    // consistency constant
    double delta  = material->m.modpowerlaw->delta;     // safety factor
    double a      = material->m.modpowerlaw->a_exp;     // exponent

    // compute viscosity according to a modified power law model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    visc = m * pow((delta + rateofshear), (-1)*a);
  }
  else
    dserror("material type not yet implemented");
} // DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalVisc


/*----------------------------------------------------------------------*
 |                                                                      |
 |  calculation of stabilisation parameter          gammi 12/08         |
 |    (element center or Gaussian point)                                |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalcTau(
  const enum Fluid3::TauType             whichtau  ,
  const enum Fluid3::StabilisationAction tds       ,
  const double &                         gamma     ,
  const double &                         dt        ,
  const double &                         hk        ,
  const double &                         mk        ,
  const double &                         visceff   )
{
  // get velocity norms
  const double vel_normaf = velintaf_.Norm2();
  const double vel_normnp = velintnp_.Norm2();

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
      LINALG::Matrix<3,3> G;

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
      LINALG::Matrix<3,1> g;

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
    else if(whichtau == Fluid3::franca_barrenechea_valentin_wall
            ||
            whichtau == Fluid3::fbvw_wo_dt                      )
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
    else if(whichtau == Fluid3::smoothed_franca_barrenechea_valentin_wall)
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
      dserror("Unknown definition of stabilisation parameter for time-dependent formulation\n");
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
      LINALG::Matrix<3,3> G;
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
      LINALG::Matrix<3,1> g;

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
                     C            /     \
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
      LINALG::Matrix<3,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=velintaf_(rr)/vel_normaf;
        }      
      }
      else
      {
        normed_velintaf(0) = 1.;
        for (int rr=1;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=0.0;
        }      
      }

      // get streamlength
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        val += FABS( normed_velintaf(0)*derxy_(0,rr)         
                    +normed_velintaf(1)*derxy_(1,rr)        
                    +normed_velintaf(2)*derxy_(2,rr));
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
    else if(whichtau == Fluid3::franca_barrenechea_valentin_codina)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Codina


      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=velintaf_(rr)/vel_normaf;
        }      
      }
      else
      {
        normed_velintaf(0) = 1.;
        for (int rr=1;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=0.0;
        }      
      }

      // get streamlength
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        val += FABS( normed_velintaf(0)*derxy_(0,rr)         
                    +normed_velintaf(1)*derxy_(1,rr)        
                    +normed_velintaf(2)*derxy_(2,rr));
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
    else if (whichtau == Fluid3::codina)
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
    else if(whichtau == Fluid3::fbvw_wo_dt)
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

      const double xi_tau_c = DMIN(re_convectnp,1.0);

      tau_(2) = vel_normnp * hk * 0.5 * xi_tau_c;
    }
    else
    {
      dserror("Unknown definition of stabilisation parameter for quasistatic formulation\n");
    }
  }

  return;
} // DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalcTau


/*----------------------------------------------------------------------*
 |                                                                      |
 | Get all global shape functions, first and eventually second          |
 | derivatives in a gausspoint                                          |
 |                                                         gammi 12/08  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::ShapeFunctionsFirstAndSecondDerivatives(
  const Fluid3*                                ele             ,
  const int                                 &  iquad           ,
  const DRT::UTILS::IntegrationPoints3D     &  intpoints       ,
  const std::vector<blitz::Array<double,1> >&  myknots         ,
  const bool                                   higher_order_ele
  )
{
    // set gauss point coordinates
    LINALG::Matrix<3,1> gp;

    gp(0)=intpoints.qxg[iquad][0];
    gp(1)=intpoints.qxg[iquad][1];
    gp(2)=intpoints.qxg[iquad][2];

    if(!(distype == DRT::Element::nurbs8
         ||
         distype == DRT::Element::nurbs27))
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_3D(funct_,gp(0),gp(1),gp(2),distype);
      DRT::UTILS::shape_function_3D_deriv1(deriv_,gp(0),gp(1),gp(2),distype);

      if (higher_order_ele)
      {
        // get values of shape functions and derivatives in the gausspoint
        DRT::UTILS::shape_function_3D_deriv2(derxy2_,gp(0),gp(1),gp(2),distype);
      }
    }
    else
    {
      if (higher_order_ele)
      {
        DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv_deriv2
          (funct_  ,
           deriv_  ,
           derxy2_ ,
           gp      ,
           myknots ,
           weights_,
           distype );
      }
      else
      {
        DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
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
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        xjm_(rr,mm)=deriv_(rr,0)*xyze_(mm,0);
        for(int nn=1;nn<iel;++nn)
        {
          xjm_(rr,mm)+=deriv_(rr,nn)*xyze_(mm,nn);
        }
      }
    }

    // the bm as well as the other const doubles values 
    // here will be reused later for the linear system 
    // for the second derivatives
    bm_(3,3) = xjm_(0,0)*xjm_(1,1);
    bm_(3,4) = xjm_(0,0)*xjm_(1,2);
    bm_(3,5) = xjm_(0,1)*xjm_(1,2);

    bm_(4,3) = xjm_(0,0)*xjm_(2,1);
    bm_(4,4) = xjm_(0,0)*xjm_(2,2);
    bm_(4,5) = xjm_(0,1)*xjm_(2,2);

    bm_(5,3) = xjm_(1,0)*xjm_(2,1);
    bm_(5,4) = xjm_(1,0)*xjm_(2,2);
    bm_(5,5) = xjm_(1,1)*xjm_(2,2);
   
    const double xjm_2_1_xjm_1_2=xjm_(2,1)*xjm_(1,2);
    const double xjm_2_0_xjm_1_2=xjm_(2,0)*xjm_(1,2);
    const double xjm_2_0_xjm_1_1=xjm_(2,0)*xjm_(1,1);
    const double xjm_2_1_xjm_0_2=xjm_(2,1)*xjm_(0,2);
    const double xjm_2_0_xjm_0_2=xjm_(2,0)*xjm_(0,2);
    const double xjm_2_0_xjm_0_1=xjm_(2,0)*xjm_(0,1);
    const double xjm_1_1_xjm_0_2=xjm_(1,1)*xjm_(0,2);
    const double xjm_1_0_xjm_0_2=xjm_(1,0)*xjm_(0,2);
    const double xjm_1_0_xjm_0_1=xjm_(1,0)*xjm_(0,1);
   
    // The determinant ist computed using Sarrus's rule
    const double det = xjm_(0,0)*bm_(5,5)+
                       xjm_(2,0)*bm_(3,5)+
                       xjm_(0,2)*(bm_(5,3)-xjm_2_0_xjm_1_1)-
                       xjm_(2,1)*bm_(3,4)-
                       xjm_(0,1)*bm_(5,4);

    // check for degenerated elements
    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
    }

    // set total integration factor
    double fac = intpoints.qwgt[iquad]*det;

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
    xji_(0,0) = ( bm_(5,5)-xjm_2_1_xjm_1_2)/det;
    xji_(1,0) = (-bm_(5,4)+xjm_2_0_xjm_1_2)/det;
    xji_(2,0) = ( bm_(5,3)-xjm_2_0_xjm_1_1)/det;
    xji_(0,1) = (-bm_(4,5)+xjm_2_1_xjm_0_2)/det;
    xji_(1,1) = ( bm_(4,4)-xjm_2_0_xjm_0_2)/det;
    xji_(2,1) = (-bm_(4,3)+xjm_2_0_xjm_0_1)/det;
    xji_(0,2) = ( bm_(3,5)-xjm_1_1_xjm_0_2)/det;
    xji_(1,2) = (-bm_(3,4)+xjm_1_0_xjm_0_2)/det;
    xji_(2,2) = ( bm_(3,3)-xjm_1_0_xjm_0_1)/det;

    // compute global derivates at integration point
    //
    //   dN    +-----  dN (xi)    dxi
    //     i    \        i           k
    //   --- =   +     ------- * -----
    //   dx     /        dxi      dx
    //     j   +-----       k       j
    //         node k
    //
    // j : direction of derivative x/y/z
    //
    for(int nn=0;nn<iel;++nn)
    {
      for(int rr=0;rr<3;++rr)
      {
        derxy_(rr,nn)=xji_(rr,0)*deriv_(0,nn);
        for(int mm=1;mm<3;++mm)
        {
          derxy_(rr,nn)+=xji_(rr,mm)*deriv_(mm,nn);
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
      // calculate elements of jacobian_bar matrix
      bm_(0,0)  = xjm_(0,0)*xjm_(0,0);
      bm_(1,0)  = xjm_(1,0)*xjm_(1,0);
      bm_(2,0)  = xjm_(2,0)*xjm_(2,0);
      bm_(3,0)  = xjm_(0,0)*xjm_(1,0);
      bm_(4,0)  = xjm_(0,0)*xjm_(2,0);
      bm_(5,0)  = xjm_(2,0)*xjm_(1,0);

      bm_(0,1)  = xjm_(0,1)*xjm_(0,1);
      bm_(1,1)  = xjm_(1,1)*xjm_(1,1);
      bm_(2,1)  = xjm_(2,1)*xjm_(2,1);
      bm_(3,1)  = xjm_(0,1)*xjm_(1,1);
      bm_(4,1)  = xjm_(0,1)*xjm_(2,1);
      bm_(5,1)  = xjm_(2,1)*xjm_(1,1);

      bm_(0,2)  = xjm_(0,2)*xjm_(0,2);
      bm_(1,2)  = xjm_(1,2)*xjm_(1,2);
      bm_(2,2)  = xjm_(2,2)*xjm_(2,2);
      bm_(3,2)  = xjm_(0,2)*xjm_(1,2);
      bm_(4,2)  = xjm_(0,2)*xjm_(2,2);
      bm_(5,2)  = xjm_(2,2)*xjm_(1,2);

      bm_(0,3)  = 2.*xjm_(0,0)*xjm_(0,1);
      bm_(1,3)  = 2.*xjm_(1,0)*xjm_(1,1);
      bm_(2,3)  = 2.*xjm_(2,0)*xjm_(2,1);
      bm_(3,3) += xjm_1_0_xjm_0_1;
      bm_(4,3) += xjm_2_0_xjm_0_1;
      bm_(5,3) += xjm_2_0_xjm_1_1;

      bm_(0,4)  = 2.*xjm_(0,0)*xjm_(0,2);
      bm_(1,4)  = 2.*xjm_(1,0)*xjm_(1,2);
      bm_(2,4)  = 2.*xjm_(2,0)*xjm_(2,2);
      bm_(3,4) += xjm_1_0_xjm_0_2;
      bm_(4,4) += xjm_2_0_xjm_0_2;
      bm_(5,4) += xjm_2_0_xjm_1_2;

      bm_(0,5)  = 2.*xjm_(0,1)*xjm_(0,2);
      bm_(1,5)  = 2.*xjm_(1,1)*xjm_(1,2);
      bm_(2,5)  = 2.*xjm_(2,1)*xjm_(2,2);
      bm_(3,5) += xjm_1_1_xjm_0_2;
      bm_(4,5) += xjm_2_1_xjm_0_2;
      bm_(5,5) += xjm_2_1_xjm_1_2;

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
      for(int rr=0;rr<6;++rr)
      {
        for(int mm=0;mm<3;++mm)
        {
          xder2_(rr,mm)=derxy2_(rr,0)*xyze_(mm,0);
          for(int nn=1;nn<iel;++nn)
          {
            xder2_(rr,mm)+=derxy2_(rr,nn)*xyze_(mm,nn);
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
      for(int rr=0;rr<6;++rr)
      {
        for(int nn=0;nn<iel;++nn)
        {
          for(int mm=0;mm<3;++mm)
          {
            derxy2_(rr,nn)-=xder2_(rr,mm)*derxy_(mm,nn);
          }
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

      // a vector specifying the pivots (reordering)
      int pivot[6];

      // error code
      int ierr = 0;

      // Perform LU factorisation --- this call replaces bm with its factorisation
      solver_.GETRF(6,6,bm_.A(),6,&(pivot[0]),&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform LU factorisation during computation of derxy2");
      }

      // backward substitution. GETRS replaces the input (chainrulerhs, currently
      // stored on derxy2) with the result
      solver_.GETRS('N',6,iel,bm_.A(),6,&(pivot[0]),derxy2_.A(),6,&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform backward substitution after factorisation of jacobian");
      }
    }


  return(fac);
} // ShapeFunctionsFirstAndSecondDerivatives

/*----------------------------------------------------------------------*
 |  calculates all quantities which are defined at the element center   |
 |  or for the whole element                                            |
 |     o element geometry (xyze_ etc)                                   |
 |     o element volume vol_                                            |
 |     o element size hk, constant mk from inverse estimate             |
 |     o dead load                                                      |
 |     o viscosity, effective viscosity                                 |
 |                                                         gammi 01/09  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::SetElementData(
  Fluid3*                                     ele            ,
  const LINALG::Matrix<3,iel>               & edispnp        ,
  const LINALG::Matrix<3,iel>               & evelaf         ,
  const LINALG::Matrix<3,iel>               & fsevelaf       ,
  const std::vector<blitz::Array<double,1> >& myknots        ,
  const double                              & timealphaF     ,
  double                                    & hk             ,
  double                                    & mk             ,
  const struct _MATERIAL*                     material       ,
  double                                    & visc           ,
  const enum Fluid3::StabilisationAction      fssgv          ,
  const enum Fluid3::TurbModelAction          turb_mod_action,
  const double                                l_tau          ,
  double                                    & Cs             ,
  double                                    & Cs_delta_sq    ,
  double                                    & visceff        )
{

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }
  
  // get node weights for nurbs elements
  if(distype==DRT::Element::nurbs8 || distype==DRT::Element::nurbs27)
  {
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
      
      weights_(inode) = cp->W();
    }
  }
  
  // add displacement, when fluid nodes move in the ALE case
  if (ele->is_ale_)
  {
    for (int inode=0; inode<iel; inode++)
    {
      xyze_(0,inode) += edispnp(0,inode);
      xyze_(1,inode) += edispnp(1,inode);
      xyze_(2,inode) += edispnp(2,inode);
    }
  }
  
  //----------------------------------------------------------------------------
  //                  GET DEAD LOAD IN ELEMENT NODES
  //----------------------------------------------------------------------------
  GetNodalBodyForce(ele,timealphaF);

  //------------------------------------------------------------------
  //                      SET MATERIAL DATA
  //------------------------------------------------------------------
  // check here, if we really have a fluid !!
  if( material->mattyp != m_carreauyasuda
      && 
      material->mattyp != m_modpowerlaw
      && 
      material->mattyp != m_fluid)
  dserror("Material law is not a fluid");

  // get material viscosity
  if(material->mattyp == m_fluid)
  {
    visc = material->m.fluid->viscosity;
  }
  // initialise visceff to visc
  visceff=visc;

  // ---------------------------------------------------------------------------
  // Initialisation of tau computation: mk and hk

  // get element type constant mk for tau and the fssgv_artificial approach
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
      case DRT::Element::nurbs8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::nurbs27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // use one point gauss rule to calculate volume at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili=DRT::UTILS::intrule3D_undefined;
  switch (distype)
  {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::nurbs8:
      case DRT::Element::nurbs27:
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
  const double wquad = intpoints_onepoint.qwgt[0];
        
  LINALG::Matrix<3,1> gp;
  gp(0)=intpoints_onepoint.qxg[0][0];
  gp(1)=intpoints_onepoint.qxg[0][1];
  gp(2)=intpoints_onepoint.qxg[0][2];

  if(distype == DRT::Element::nurbs8
     ||
     distype == DRT::Element::nurbs27)
  {
    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
      (funct_  ,
       deriv_  ,
       gp      ,
       myknots ,
       weights_,
       distype );
  }
  else
  {
    DRT::UTILS::shape_function_3D       (funct_,gp(0),gp(1),gp(2),distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv_,gp(0),gp(1),gp(2),distype);
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
  for(int rr=0;rr<3;++rr)
  {
    for(int mm=0;mm<3;++mm)
    {
      xjm_(rr,mm)=deriv_(rr,0)*xyze_(mm,0);
      for(int nn=1;nn<iel;++nn)
      {
        xjm_(rr,mm)+=deriv_(rr,nn)*xyze_(mm,nn);
      }
    }
  }

  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

  // vol is used for some stabilisation parameters and the mixing length
  // definition of the Smagorinsky models
  vol_ = wquad*det;

  // get element length hk for tau_M and tau_C: volume-equival. diameter/sqrt(3)
  // the same hk is used for the fssgv_artificial approach
  hk = pow((6.*vol_/PI),(1.0/3.0))/sqrt(3.0);

  /*------------------------------------------------------------------*/
  /*                                                                  */
  /*                 GET EFFECTIVE VISCOSITY IN ELEMENT               */
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
  /* the effective viscosity.                                         */
  /* This has to be done before anything else is calculated because   */
  /* we use the same arrays internally. We need hk, mk as well as the */
  /* element data computed above!                                     */
  /*------------------------------------------------------------------*/

  // -------------------------------------------------------------------
  // strain rate based models

  if(material->mattyp!= m_fluid
     ||
     fssgv==Fluid3::fssgv_Smagorinsky_all
     ||
     fssgv==Fluid3::fssgv_Smagorinsky_small
     ||
     turb_mod_action !=Fluid3::no_model)
  {
    //
    //             compute global first derivates
    //
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
    
    // compute global derivates at integration point
    //
    //   dN    +-----  dN (xi)    dxi
    //     i    \        i           k
    //   --- =   +     ------- * -----
    //   dx     /        dxi      dx
    //     j   +-----       k       j
    //         node k
    //
    // j : direction of derivative x/y/z
    //
    for(int nn=0;nn<iel;++nn)
    {
      for(int rr=0;rr<3;++rr)
      {
        derxy_(rr,nn)=xji_(rr,0)*deriv_(0,nn);
        for(int mm=1;mm<3;++mm)
        {
          derxy_(rr,nn)+=xji_(rr,mm)*deriv_(mm,nn);
        }
      }
    }

    // compute the rate of strain
    double rateofstrain   = -1.0e30;
    double fsrateofstrain = -1.0e30;

    // large scale strains
    if(material->mattyp!=m_fluid
       ||
       fssgv==Fluid3::fssgv_Smagorinsky_all
       ||
       turb_mod_action!=Fluid3::no_model)
    {
      rateofstrain  =GetStrainRate(evelaf  ,derxy_,vderxyaf_  );
    }

    // fine scale strains
    if(fssgv==Fluid3::fssgv_Smagorinsky_small)
    {
      fsrateofstrain=GetStrainRate(fsevelaf,derxy_,fsvderxyaf_);
    }


    // ---------------------------------------------------------------
    // compute nonlinear viscosity according to Carreau-Yasuda
    // ---------------------------------------------------------------
    if(material->mattyp != m_fluid)
    {
      CalVisc(material,visceff,rateofstrain);
    }


    // ---------------------------------------------------------------
    // classical large scale turbulence models of Smagorinsky type
    // ---------------------------------------------------------------

    //
    // Choices of the Smagorinsky constant Cs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)
    //
    //             Cs rescaled using law of the wall (Van Driest)
    //
    //             Cs dynamic  (Germano model. Use several filter
    //                          resolutions to determine Cs)
    if (turb_mod_action == Fluid3::smagorinsky_with_wall_damping
        ||
        turb_mod_action == Fluid3::smagorinsky)
    {
      //
      // CONSTANT COEFFICIENT SMAGORINSKY MODEL
      // --------------------------------------
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

      if (turb_mod_action == Fluid3::smagorinsky_with_wall_damping)
      {
        //
        // VAN DRIEST DAMPING (option)
        // ------------------
        
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
        LINALG::Matrix<3,1> centernodecoord;
        
        for(int rr=0;rr<3;++rr)
        {
          centernodecoord(rr)=funct_(0)*xyze_(rr,0);
          
          for(int mm=1;mm<iel;++mm)
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
      
      const double h_grid = pow((vol_),(1.0/3.0));

      //
      // mixing length set proportional to grid witdh
      //
      //              lmix = Cs * h_grid
      
      double lmix = Cs * h_grid;
      
      Cs_delta_sq = lmix * lmix;

      //
      //          visc    = visc + visc
      //              eff              turbulent
      
      visceff = visc + Cs_delta_sq * rateofstrain;
    }
    else if(turb_mod_action == Fluid3::dynamic_smagorinsky)
    {
      //
      // DYNAMIC SMAGORINSKY MODEL
      // -------------------------
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
      
      visceff = visc + Cs_delta_sq * rateofstrain;

      // for evaluation of statistics: remember the 'real' Cs
      Cs=sqrt(Cs_delta_sq)/pow((vol_),(1.0/3.0));
    }


    // ---------------------------------------------------------------
    // fine scale turbulence models of Smagorinsky type
    // ---------------------------------------------------------------

    if (fssgv == Fluid3::fssgv_Smagorinsky_all)
    {
      //
      // SMALL SCALE LARGE SCALE SMAGORINSKY MODEL
      // -----------------------------------------
      //                               +-                                 -+ 1
      //                           2   |          / * \           / * \    | -
      //    visc          = (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
      //        turbulent              |          \   / ij        \   / ij |
      //                               +-                                 -+
      //                               |                                   |
      //                               +-----------------------------------+
      //                                  * 'resolved' rate of strain
      
      // computes the turbulent viscosity (acting on the fine scales) 
      // based on the resolved-scale rate of strain
      //
      // Choices of the fine-scale Smagorinsky constant Cs_fs:
      //
      //             Cs = 0.17   (Lilly --- Determined from filter
      //                          analysis of Kolmogorov spectrum of
      //                          isotropic turbulence)
      //
      //             0.1 < Cs < 0.24 (depending on the flow)
      //
      vart_ = Cs * Cs * hk * hk * rateofstrain;
    }
    else if(fssgv == Fluid3::fssgv_Smagorinsky_small)
    {
      //
      // SMALL SCALE SMALL SCALE SMAGORINSKY MODEL
      // -----------------------------------------
      //                               +-                                 -+ 1
      //                           2   |          /   \           /   \    | -
      //    visc          = (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
      //        turbulent              |          \   / ij        \   / ij |
      //                               +-                                 -+
      //                               |                                   |
      //                               +-----------------------------------+
      //                                     'fine-scale'  rate of strain
      
      // compute turbulent viscosity (acting on the fine scales) 
      // based on the fine-scale rate of strain
    
      //
      // Choices of the fine-scale Smagorinsky constant Cs_fs:
      //
      //             Cs = 0.17   (Lilly --- Determined from filter
      //                          analysis of Kolmogorov spectrum of
      //                          isotropic turbulence)
      //
      //             0.1 < Cs < 0.24 (depending on the flow)
      //
      vart_ = Cs * Cs * hk * hk * fsrateofstrain;
    }
  }

  // ------------------------------------------------------------------
  // fine scale models based on artificial subgrid viscosities
  // ------------------------------------------------------------------

  if (fssgv == Fluid3::fssgv_artificial_all   or
      fssgv == Fluid3::fssgv_artificial_small
    )
  {
    double fsvel_normaf = 0.0;
    if (fssgv == Fluid3::fssgv_artificial_small)
    {
      // SMALL SCALE ARTIFICIAL SUBGRID VISCOSITY MODEL
      // ----------------------------------------------

      // get fine-scale velocities at element center
      fsvelintaf_.Multiply(fsevelaf,funct_);

      // get fine-scale velocity norm
      fsvel_normaf = fsvelintaf_.Norm2();
    }
    // get all-scale velocity norm
    else 
    {
      // LARGE SCALE ARTIFICIAL SUBGRID VISCOSITY MODEL
      // ----------------------------------------------

      // get velocities (n+alpha_F,i) at integration point
      velintaf_.Multiply(evelaf,funct_);

      // get velocity norms
      const double vel_normaf = velintaf_.Norm2();

      fsvel_normaf = vel_normaf;
    }

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = mk * fsvel_normaf * hk / visc; /* convective : viscous forces */
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(fsvel_normaf))/(2.0*visc*xi);
  }

  return;
} // SetElementData

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
 |        o fsvderxyaf_ (if fssgv)                                      |
 |        o viscs2_ (if hoel)                                           |
 |        o viscaf_old_ (if hoel)                                       |
 |                                                         gammi 01/09  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::InterpolateToGausspoint(
  Fluid3*                                     ele             ,
  const LINALG::Matrix<3,iel>               & egridvaf        ,
  const LINALG::Matrix<3,iel>               & evelnp          ,
  const LINALG::Matrix<iel,1>               & eprenp          ,
  const LINALG::Matrix<3,iel>               & eaccam          ,
  const LINALG::Matrix<3,iel>               & evelaf          ,
  const LINALG::Matrix<3,iel>               & fsevelaf        ,
  const double                              & visceff         ,
  const enum Fluid3::StabilisationAction      fssgv           ,
  const bool                                  higher_order_ele)
{
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
  for(int rr=0;rr<3;++rr)
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
  for(int rr=0;rr<3;++rr)
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
  for(int rr=0;rr<3;++rr)
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
    for(int rr=0;rr<3;++rr)
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
    for(int rr=0;rr<3;++rr)
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
  for(int rr=0;rr<3;++rr)
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
  for(int rr=0;rr<3;++rr)
  {
    for(int mm=0;mm<3;++mm)
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
  for(int rr=0;rr<3;++rr)
  {
    for(int mm=0;mm<3;++mm)
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
  divunp_ = (vderxynp_(0,0)+vderxynp_(1,1)+vderxynp_(2,2));

  // get ale convective velocity (n+alpha_F,i) at integration point
  for(int rr=0;rr<3;++rr)
  {
    aleconvintaf_(rr)=velintaf_(rr);
  }    

  if (ele->is_ale_)
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
      
    for(int rr=0;rr<3;++rr)
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

    for(int rr=0;rr<3;++rr)
    {
      aleconvintaf_(rr)-=u_G_af_(rr);
    }
  }
  else
  {
    for(int rr=0;rr<3;++rr)
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
  for(int rr=0;rr<3;++rr)
  {
    convaf_old_(rr)=aleconvintaf_(0)*vderxyaf_(rr,0);
    for(int mm=1;mm<3;++mm)
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
  for (int rr=0;rr<3;++rr)
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

    for(int rr=1;rr<3;++rr)
    {
      conv_c_af_(nn)+=aleconvintaf_(rr)*derxy_(rr,nn);
    }
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
  if (fssgv != Fluid3::fssgv_no) fsvderxyaf_.MultiplyNT(fsevelaf,derxy_);
  else fsvderxyaf_.Clear();

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
    double sum = derxy2_(0,0)+derxy2_(1,0)+derxy2_(2,0);

    viscs2_(0,0) = sum + derxy2_(0,0);
    viscs2_(1,0) = sum + derxy2_(1,0);
    viscs2_(2,0) = sum + derxy2_(2,0);

    viscaf_old_(0) =
      (viscs2_(0,0)*evelaf(0,0)
       +
       derxy2_(3,0)*evelaf(1,0)
       +
       derxy2_(4,0)*evelaf(2,0));
    viscaf_old_(1) =
      (derxy2_(3,0)*evelaf(0,0)
       +
       viscs2_(1,0)*evelaf(1,0)
       +
       derxy2_(5,0)*evelaf(2,0));
    viscaf_old_(2) =
      (derxy2_(4,0)*evelaf(0,0)
       +
       derxy2_(5,0)*evelaf(1,0)
       +
       viscs2_(2,0)*evelaf(2,0));

    for (int mm=1;mm<iel;++mm)
    {
      sum = derxy2_(0,mm)+derxy2_(1,mm)+derxy2_(2,mm);

      viscs2_(0,mm) = sum + derxy2_(0,mm);
      viscs2_(1,mm) = sum + derxy2_(1,mm);
      viscs2_(2,mm) = sum + derxy2_(2,mm);

      viscaf_old_(0) += (viscs2_(0,mm)*evelaf(0,mm)
                         +
                         derxy2_(3,mm)*evelaf(1,mm)
                         +
                         derxy2_(4,mm)*evelaf(2,mm));
      viscaf_old_(1) += (derxy2_(3,mm)*evelaf(0,mm)
                         +
                         viscs2_(1,mm)*evelaf(1,mm)
                         +
                         derxy2_(5,mm)*evelaf(2,mm));
      viscaf_old_(2) += (derxy2_(4,mm)*evelaf(0,mm)
                         +
                         derxy2_(5,mm)*evelaf(1,mm)
                         +
                         viscs2_(2,mm)*evelaf(2,mm));
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
    for (int rr=0;rr<3;++rr)
    {
      resM_(rr) -= visceff*viscaf_old_(rr);
    }
  } // end if higher order

  return;
} // InterpolateToGausspoint

#endif
#endif
