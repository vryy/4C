/*----------------------------------------------------------------------*/
/*!
\file fluid3_genalpha_resVMM.cpp

\brief Internal implementation of Fluid element with a generalised alpha
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


#include "fluid_ele.H"
#include "fluid_genalpha_resVMM.H"
#include "fluid_genalpha_resVMM_2D.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidGenalphaResVMMInterface* DRT::ELEMENTS::FluidGenalphaResVMMInterface::Impl(DRT::ELEMENTS::Fluid* f3)
{
  const int numdofpernode = f3->NumDofPerNode(*(f3->Nodes()[0]));
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static FluidGenalphaResVMM<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new FluidGenalphaResVMM<DRT::Element::hex8>(numdofpernode);
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static FluidGenalphaResVMM<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new FluidGenalphaResVMM<DRT::Element::hex20>(numdofpernode);
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static FluidGenalphaResVMM<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new FluidGenalphaResVMM<DRT::Element::hex27>(numdofpernode);
    return fh27;
  }
  case DRT::Element::nurbs8:
  {
    static FluidGenalphaResVMM<DRT::Element::nurbs8>* fn8;
    if (fn8==NULL)
      fn8 = new FluidGenalphaResVMM<DRT::Element::nurbs8>(numdofpernode);
    return fn8;
  }
  case DRT::Element::nurbs27:
  {
    static FluidGenalphaResVMM<DRT::Element::nurbs27>* fn27;
    if (fn27==NULL)
      fn27 = new FluidGenalphaResVMM<DRT::Element::nurbs27>(numdofpernode);
    return fn27;
  }
  case DRT::Element::tet4:
  {
    static FluidGenalphaResVMM<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new FluidGenalphaResVMM<DRT::Element::tet4>(numdofpernode);
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static FluidGenalphaResVMM<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new FluidGenalphaResVMM<DRT::Element::tet10>(numdofpernode);
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static FluidGenalphaResVMM<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new FluidGenalphaResVMM<DRT::Element::wedge6>(numdofpernode);
    return fw6;
  }
  /*
  case DRT::Element::wedge15:
  {
    static FluidGenalphaResVMM<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new FluidGenalphaResVMM<DRT::Element::wedge15>(numdofpernode);
    return fw15;
  }*/
  case DRT::Element::pyramid5:
  {
    static FluidGenalphaResVMM<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new FluidGenalphaResVMM<DRT::Element::pyramid5>(numdofpernode);
    return fp5;
  }
  case DRT::Element::quad4:
  {
    static FluidGenalphaResVMM<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new FluidGenalphaResVMM<DRT::Element::quad4>(numdofpernode);
    return cp4;
  }
  case DRT::Element::quad8:
  {
    static FluidGenalphaResVMM<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new FluidGenalphaResVMM<DRT::Element::quad8>(numdofpernode);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static FluidGenalphaResVMM<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new FluidGenalphaResVMM<DRT::Element::quad9>(numdofpernode);
    return cp9;
  }
  case DRT::Element::tri3:
  {
    static FluidGenalphaResVMM<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new FluidGenalphaResVMM<DRT::Element::tri3>(numdofpernode);
    return cp3;
  }
  case DRT::Element::tri6:
  {
    static FluidGenalphaResVMM<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new FluidGenalphaResVMM<DRT::Element::tri6>(numdofpernode);
    return cp6;
  }
  case DRT::Element::nurbs4:
  {
    static FluidGenalphaResVMM<DRT::Element::nurbs4>* fn4;
    if (fn4==NULL)
      fn4 = new FluidGenalphaResVMM<DRT::Element::nurbs4>(numdofpernode);
    return fn4;
  }
  case DRT::Element::nurbs9:
  {
    static FluidGenalphaResVMM<DRT::Element::nurbs9>* fn9;
    if (fn9==NULL)
      fn9 = new FluidGenalphaResVMM<DRT::Element::nurbs9>(numdofpernode);
    return fn9;
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
DRT::ELEMENTS::FluidGenalphaResVMM<distype>::FluidGenalphaResVMM(int numdofpernode)
// fine-scale subgrid viscosity
  : numdofpernode_(numdofpernode),
    vart_(0.0)

{
  return;
}


template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Evaluate(
  Fluid*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RefCountPtr<MAT::Material> mat)
{
  // --------------------------------------------------
  // construct views

  // standard fluid linearisation
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel> elemat1(elemat1_epetra.A(),true);
  // linearisation of fluid with respect to mesh motion
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel> elemat2(elemat2_epetra.A(),true);
  // residual vector (right hand side of equation system)
  LINALG::Matrix<(nsd_+1)*iel,    1>        elevec1(elevec1_epetra.A(),true);


  // --------------------------------------------------
  // create matrix objects for nodal values
  LINALG::Matrix<iel,1> eprenp    ;
  LINALG::Matrix<nsd_,iel> evelnp    ;
  LINALG::Matrix<nsd_,iel> evelaf    ;
  LINALG::Matrix<nsd_,iel> eaccam    ;
  LINALG::Matrix<nsd_,iel> edispnp   ;
  LINALG::Matrix<nsd_,iel> egridvelaf;
  LINALG::Matrix<nsd_,iel> fsevelaf  ;

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
  INPAR::FLUID::LinearisationAction newton = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(params, "Linearisation");

  // --------------------------------------------------
  // get flag for fine-scale subgrid-viscosity approach
  INPAR::FLUID::FineSubgridVisc fssgv = INPAR::FLUID::no_fssgv;
  {
    const string fssgvdef = params.get<string>("fs subgrid viscosity","No");

    if (fssgvdef == "Smagorinsky_all")        fssgv = INPAR::FLUID::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv = INPAR::FLUID::smagorinsky_small;
  }

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");

  // specify which residual based stabilisation terms
  // will be used
  INPAR::FLUID::SubscalesTD             tds = DRT::INPUT::IntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
  INPAR::FLUID::Transient       inertia = DRT::INPUT::IntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
  INPAR::FLUID::PSPG            pspg = DRT::INPUT::IntegralValue<INPAR::FLUID::PSPG>(stablist,"PSPG");
  INPAR::FLUID::SUPG            supg = DRT::INPUT::IntegralValue<INPAR::FLUID::SUPG>(stablist,"SUPG");
  INPAR::FLUID::VStab           vstab = DRT::INPUT::IntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
  INPAR::FLUID::CStab           cstab = DRT::INPUT::IntegralValue<INPAR::FLUID::CStab>(stablist,"CSTAB");
  INPAR::FLUID::CrossStress     cross = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
  INPAR::FLUID::ReynoldsStress  reynolds = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");

  // select tau definition
  INPAR::FLUID::TauType_genalpha whichtau = INPAR::FLUID::tautype_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");

    if(taudef == "Franca_Barrenechea_Valentin_Frey_Wall")
    {
      whichtau = INPAR::FLUID::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Taylor_Hughes_Zarins_Whiting_Jansen")
    {
      whichtau = INPAR::FLUID::taylor_hughes_zarins_whiting_jansen;
    }
    else if(taudef == "Codina")
    {
      whichtau = INPAR::FLUID::codina;
    }
    else if(taudef == "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt")
    {
      whichtau = INPAR::FLUID::fbvw_wo_dt;
    }
    else if(taudef == "Franca_Barrenechea_Valentin_Codina")
    {
      whichtau = INPAR::FLUID::franca_barrenechea_valentin_codina;
    }
    else if(taudef == "Smoothed_FBVW")
    {
      whichtau = INPAR::FLUID::smoothed_franca_barrenechea_valentin_wall;
    }
    else if(taudef == "BFVW_gradient_based_hk")
    {
      whichtau = INPAR::FLUID::fbvw_gradient_based_hk;
    }
    else
    {
      dserror("unknown tau definition\n");
    }
  }

  // flag for higher order elements
  //bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());
  bool higher_order_ele = IsHigherOrder<distype>::ishigherorder;

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
  ParameterList& sgviscparams    = params.sublist("SUBGRID VISCOSITY");

  // the default action is no model
  INPAR::FLUID::TurbModelAction turb_mod_action = INPAR::FLUID::no_model;

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  double Cs            = 0.0;
  double Cs_delta_sq   = 0.0;
  double l_tau         = 0.0;

  // number of the layer in a turbulent plane channel flow --- used
  // to compute averaged viscosity etc
  int    nlayer        = 0;

  // turbulence model is only implemented for 3D flows
  if(nsd_ == 3)
    SetParametersForTurbulenceModel(
      ele            ,
      turbmodelparams,
      sgviscparams   ,
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
  if (nsd_ == 3)
    ExtractValuesFromGlobalVectors(
          fssgv         ,
          ele->IsAle()  ,
          discretization,
          lm,
          eprenp        ,
          evelnp        ,
          evelaf        ,
          eaccam        ,
          edispnp       ,
          egridvelaf    ,
          fsevelaf
      );
  else if (nsd_==2)
    ExtractValuesFromGlobalVectors_2D(
          ele->IsAle()  ,
          discretization,
          lm,
          eprenp        ,
          evelnp        ,
          evelaf        ,
          eaccam        ,
          edispnp       ,
          egridvelaf
      );
  else dserror("Only 2D and 3D fluid supported");

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(nsd_);

  // for isogeometric elements
  //if(ele->Shape()==Fluid::nurbs8 || ele->Shape()==Fluid::nurbs27)
  if(IsNurbs<distype>::isnurbs)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    bool zero_size = false;
    zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

    // if we have a zero sized element due to a interpolated
    // point --- exit here
    if(zero_size)
    {
      return(0);
    }
  }

  // on output of Sysmat, visceff will contain the computed effective viscosity
  double visceff       = 0.0;

  // --------------------------------------------------
  // calculate element coefficient matrix

  // call 3D Sysmat
  if (nsd_ == 3)
  {
  if(!conservative)
  {
    if(tds!=INPAR::FLUID::subscales_time_dependent)
    {
      Sysmat_adv_qs(
        ele,
        myknots,
        elemat1,
        elemat2,
        elevec1,
        edispnp,
        egridvelaf,
        evelnp,
        eprenp,
        eaccam,
        evelaf,
        fsevelaf,
        mat,
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
        mat,
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
    if(tds!=INPAR::FLUID::subscales_time_dependent)
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
	mat,
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
        mat             ,
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
  } //end 3D Sysmat
  // call 2D Sysmat
  else if (nsd_ == 2)
  {
    if(!conservative)
    {
      if(tds!=INPAR::FLUID::subscales_time_dependent)
      {
        Sysmat_adv_qs_2D(
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
          mat,
          alphaM,
          alphaF,
          gamma,
          dt,
          time,
          newton,
          higher_order_ele,
          inertia,
          pspg,
          supg,
          vstab,
          cstab,
          cross,
          reynolds,
          whichtau,
          visceff,
          compute_elemat
          );
      }
      else
      {
        Sysmat_adv_td_2D(
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
          mat,
          alphaM,
          alphaF,
          gamma,
          dt,
          time,
          newton,
          higher_order_ele,
          inertia,
          pspg,
          supg,
          vstab,
          cstab,
          cross,
          reynolds,
          whichtau,
          visceff,
          compute_elemat
          );
      }
    }
    else
    {
      if(tds!=INPAR::FLUID::subscales_time_dependent)
      {
        Sysmat_cons_qs_2D(
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
    mat,
    alphaM,
    alphaF,
    gamma,
    dt,
    time,
    newton,
    higher_order_ele,
    pspg,
    supg,
    vstab,
    cstab,
    cross,
    reynolds,
    whichtau,
    visceff,
    compute_elemat
    );
      }
      else
      {
        Sysmat_cons_td_2D(
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
          mat             ,
          alphaM          ,
          alphaF          ,
          gamma           ,
          dt              ,
          time            ,
          newton          ,
          higher_order_ele,
          inertia         ,
          pspg            ,
          supg            ,
          vstab           ,
          cstab           ,
          cross           ,
          reynolds        ,
          whichtau        ,
          visceff         ,
          compute_elemat
          );
      }
    }
  } //end 2D Sysmat
  else dserror("Genalpha solver does supported only Fluid2 and FLUID3 elements");

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
          ||
          physical_turbulence_model ==  "Smagorinsky"
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
    if(mat->MaterialType()== INPAR::MAT::m_fluid)
    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
      dens = actmat->Density();
    }
    else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
    {
      const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
      dens = actmat->Density();
    }
    else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
    {
      const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
      dens = actmat->Density();
    }
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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_adv_qs(
  Fluid*                                    ele             ,
  std::vector<Epetra_SerialDenseVector> &    myknots         ,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat          ,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& meshmat         ,
  LINALG::Matrix<(nsd_+1)*iel,    1>&        elevec          ,
  const LINALG::Matrix<nsd_,iel>&            edispnp         ,
  const LINALG::Matrix<nsd_,iel>&            egridvelaf      ,
  const LINALG::Matrix<nsd_,iel>&            evelnp          ,
  const LINALG::Matrix<iel,1>&               eprenp          ,
  const LINALG::Matrix<nsd_,iel>&            eaccam          ,
  const LINALG::Matrix<nsd_,iel>&            evelaf          ,
  const LINALG::Matrix<nsd_,iel>&            fsevelaf        ,
  Teuchos::RCP<const MAT::Material>          material        ,
  const double                               alphaM          ,
  const double                               alphaF          ,
  const double                               gamma           ,
  const double                               dt              ,
  const double                               time            ,
  const enum INPAR::FLUID::LinearisationAction     newton          ,
  const bool                                 higher_order_ele,
  const enum INPAR::FLUID::FineSubgridVisc         fssgv           ,
  const enum INPAR::FLUID::Transient         inertia         ,
  const enum INPAR::FLUID::PSPG              pspg            ,
  const enum INPAR::FLUID::SUPG              supg            ,
  const enum INPAR::FLUID::VStab             vstab           ,
  const enum INPAR::FLUID::CStab             cstab           ,
  const enum INPAR::FLUID::CrossStress       cross           ,
  const enum INPAR::FLUID::ReynoldsStress    reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha                 whichtau        ,
  const enum INPAR::FLUID::TurbModelAction         turb_mod_action ,
  double&                                    Cs              ,
  double&                                    Cs_delta_sq     ,
  double&                                    visceff         ,
  const double                               l_tau           ,
  const bool                                 compute_elemat
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
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
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

        for(int rr=1;rr<3;++rr)
        {
          conv_resM_(nn)+=resM_(rr)*derxy_(rr,nn);
        }
      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,INPAR::FLUID::subscales_quasistatic,gamma,dt,hk,mk,visceff);


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

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
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


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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
      if (newton==INPAR::FLUID::Newton)
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

        if (newton==INPAR::FLUID::Newton)
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
        if((vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_usfem)&&higher_order_ele)
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

          if (newton==INPAR::FLUID::Newton)
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
      if(cross == INPAR::FLUID::cross_stress_stab)
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

    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs || cross == INPAR::FLUID::cross_stress_stab)
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

      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
        elevec(4*ui+3)-=fac_tauMp*(resM_(0)*derxy_(0,ui)+resM_(1)*derxy_(1,ui)+resM_(2)*derxy_(2,ui));
      }
    } // end pspg

    if(supg == INPAR::FLUID::convective_stab_supg)
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
      if(vstab != INPAR::FLUID::viscous_stab_none && higher_order_ele)
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

    if(fssgv != INPAR::FLUID::no_fssgv)
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


    // linearisation with respect to mesh motion
    if (meshmat.IsInitialized())
    {
      const double fac_acc_minus_bf_x=fac*(accintam_(0)-bodyforceaf_(0));
      const double fac_acc_minus_bf_y=fac*(accintam_(1)-bodyforceaf_(1));
      const double fac_acc_minus_bf_z=fac*(accintam_(2)-bodyforceaf_(2));

      // inertia + dead load
      for (int vi=0; vi<iel; ++vi)
      {
	const int fvi  =4*vi;
	const int fvip =fvi+1;
	const int fvipp=fvip+1;

	const double fac_acc_minus_bf_x_funct_vi=fac_acc_minus_bf_x*funct_(vi);
	const double fac_acc_minus_bf_y_funct_vi=fac_acc_minus_bf_y*funct_(vi);
	const double fac_acc_minus_bf_z_funct_vi=fac_acc_minus_bf_z*funct_(vi);

        for (int ui=0; ui<iel; ++ui)
        {
	  const int fui  =4*ui;
	  const int fuip =fui+1;
	  const int fuipp=fuip+1;

          meshmat(fvi  ,fui  ) += fac_acc_minus_bf_x_funct_vi*derxy_(0,ui);
          meshmat(fvi  ,fuip ) += fac_acc_minus_bf_x_funct_vi*derxy_(1,ui);
          meshmat(fvi  ,fuipp) += fac_acc_minus_bf_x_funct_vi*derxy_(2,ui);

          meshmat(fvip ,fui  ) += fac_acc_minus_bf_y_funct_vi*derxy_(0,ui);
          meshmat(fvip ,fuip ) += fac_acc_minus_bf_y_funct_vi*derxy_(1,ui);
          meshmat(fvip ,fuipp) += fac_acc_minus_bf_y_funct_vi*derxy_(2,ui);

          meshmat(fvipp,fui  ) += fac_acc_minus_bf_z_funct_vi*derxy_(0,ui);
          meshmat(fvipp,fuip ) += fac_acc_minus_bf_z_funct_vi*derxy_(1,ui);
          meshmat(fvipp,fuipp) += fac_acc_minus_bf_z_funct_vi*derxy_(2,ui);
        }
      }


      //vderiv_  = sum(evelaf(i,k) * deriv_(j,k), k);
      vderiv_.MultiplyNT(evelaf,deriv_);

#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
#define derxjm_002(ui) (deriv_(1, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(1, 1))

#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
#define derxjm_102(ui) (deriv_(2, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(2, 0))

#define derxjm_200(ui) (deriv_(2, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(2, 1))
#define derxjm_201(ui) (deriv_(1, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(1, 0))

#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
#define derxjm_012(ui) (deriv_(2, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(2, 1))

#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))
#define derxjm_112(ui) (deriv_(0, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(0, 0))

#define derxjm_210(ui) (deriv_(0, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(0, 1))
#define derxjm_211(ui) (deriv_(2, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(2, 0))

#define derxjm_021(ui) (deriv_(1, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(1, 2))
#define derxjm_022(ui) (deriv_(0, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(0, 1))

#define derxjm_120(ui) (deriv_(0, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(0, 2))
#define derxjm_122(ui) (deriv_(1, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(1, 0))

#define derxjm_220(ui) (deriv_(1, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(1, 1))
#define derxjm_221(ui) (deriv_(0, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(0, 0))

      for (int ui=0; ui<iel; ++ui)
      {
        double v00 = + convaf_old_(1)*(vderiv_(0, 0)*derxjm_(0,0,1,ui)
				       +
				       vderiv_(0, 1)*derxjm_(0,1,1,ui)
				       +
				       vderiv_(0, 2)*derxjm_(0,2,1,ui))
                     + convaf_old_(2)*(vderiv_(0, 0)*derxjm_(0,0,2,ui)
				       +
				       vderiv_(0, 1)*derxjm_(0,1,2,ui)
				       +
				       vderiv_(0, 2)*derxjm_(0,2,2,ui));
        double v01 = + convaf_old_(0)*(vderiv_(0, 0)*derxjm_(1,0,0,ui)
				       +
				       vderiv_(0, 1)*derxjm_(1,1,0,ui)
				       +
				       vderiv_(0, 2)*derxjm_(1,2,0,ui))
                     + convaf_old_(2)*(vderiv_(0, 0)*derxjm_(1,0,2,ui)
				       +
				       vderiv_(0, 1)*derxjm_(1,1,2,ui)
				       +
				       vderiv_(0, 2)*derxjm_(1,2,2,ui));
        double v02 = + convaf_old_(0)*(vderiv_(0, 0)*derxjm_(2,0,0,ui)
				       +
				       vderiv_(0, 1)*derxjm_(2,1,0,ui)
				       +
				       vderiv_(0, 2)*derxjm_(2,2,0,ui))
                     + convaf_old_(1)*(vderiv_(0, 0)*derxjm_(2,0,1,ui)
				       +
				       vderiv_(0, 1)*derxjm_(2,1,1,ui)
				       +
				       vderiv_(0, 2)*derxjm_(2,2,1,ui));
        double v10 = + convaf_old_(1)*(vderiv_(1, 0)*derxjm_(0,0,1,ui)
				       +
				       vderiv_(1, 1)*derxjm_(0,1,1,ui)
				       +
				       vderiv_(1, 2)*derxjm_(0,2,1,ui))
                     + convaf_old_(2)*(vderiv_(1, 0)*derxjm_(0,0,2,ui)
				       +
				       vderiv_(1, 1)*derxjm_(0,1,2,ui)
				       +
				       vderiv_(1, 2)*derxjm_(0,2,2,ui));
        double v11 = + convaf_old_(0)*(vderiv_(1, 0)*derxjm_(1,0,0,ui)
				       +
				       vderiv_(1, 1)*derxjm_(1,1,0,ui)
				       +
				       vderiv_(1, 2)*derxjm_(1,2,0,ui))
                     + convaf_old_(2)*(vderiv_(1, 0)*derxjm_(1,0,2,ui)
				       +
				       vderiv_(1, 1)*derxjm_(1,1,2,ui)
				       +
				       vderiv_(1, 2)*derxjm_(1,2,2,ui));
        double v12 = + convaf_old_(0)*(vderiv_(1, 0)*derxjm_(2,0,0,ui)
				       +
				       vderiv_(1, 1)*derxjm_(2,1,0,ui)
				       +
				       vderiv_(1, 2)*derxjm_(2,2,0,ui))
                     + convaf_old_(1)*(vderiv_(1, 0)*derxjm_(2,0,1,ui)
				       +
				       vderiv_(1, 1)*derxjm_(2,1,1,ui)
				       +
				       vderiv_(1, 2)*derxjm_(2,2,1,ui));
        double v20 = + convaf_old_(1)*(vderiv_(2, 0)*derxjm_(0,0,1,ui)
				       +
				       vderiv_(2, 1)*derxjm_(0,1,1,ui)
				       +
				       vderiv_(2, 2)*derxjm_(0,2,1,ui))
                     + convaf_old_(2)*(vderiv_(2, 0)*derxjm_(0,0,2,ui)
				       +
				       vderiv_(2, 1)*derxjm_(0,1,2,ui)
				       +
				       vderiv_(2, 2)*derxjm_(0,2,2,ui));
        double v21 = + convaf_old_(0)*(vderiv_(2, 0)*derxjm_(1,0,0,ui)
				       +
				       vderiv_(2, 1)*derxjm_(1,1,0,ui)
				       +
				       vderiv_(2, 2)*derxjm_(1,2,0,ui))
                     + convaf_old_(2)*(vderiv_(2, 0)*derxjm_(1,0,2,ui)
				       +
				       vderiv_(2, 1)*derxjm_(1,1,2,ui)
				       +
				       vderiv_(2, 2)*derxjm_(1,2,2,ui));
        double v22 = + convaf_old_(0)*(vderiv_(2, 0)*derxjm_(2,0,0,ui)
				       +
				       vderiv_(2, 1)*derxjm_(2,1,0,ui)
				       +
				       vderiv_(2, 2)*derxjm_(2,2,0,ui))
                     + convaf_old_(1)*(vderiv_(2, 0)*derxjm_(2,0,1,ui)
				       +
				       vderiv_(2, 1)*derxjm_(2,1,1,ui)
				       +
				       vderiv_(2, 2)*derxjm_(2,2,1,ui));

        for (int vi=0; vi<iel; ++vi)
        {
          double v = fac/det_*funct_(vi);

          meshmat(vi*4    , ui*4    ) += v*v00;
          meshmat(vi*4    , ui*4 + 1) += v*v01;
          meshmat(vi*4    , ui*4 + 2) += v*v02;

          meshmat(vi*4 + 1, ui*4    ) += v*v10;
          meshmat(vi*4 + 1, ui*4 + 1) += v*v11;
          meshmat(vi*4 + 1, ui*4 + 2) += v*v12;

          meshmat(vi*4 + 2, ui*4    ) += v*v20;
          meshmat(vi*4 + 2, ui*4 + 1) += v*v21;
          meshmat(vi*4 + 2, ui*4 + 2) += v*v22;
        }
      }


      // viscosity

#define xji_00 xji_(0,0)
#define xji_01 xji_(0,1)
#define xji_02 xji_(0,2)
#define xji_10 xji_(1,0)
#define xji_11 xji_(1,1)
#define xji_12 xji_(1,2)
#define xji_20 xji_(2,0)
#define xji_21 xji_(2,1)
#define xji_22 xji_(2,2)

#define xjm(i,j) xjm_(i,j)

      // part 1: derivative of 1/det

      double v = visceff*fac;
      for (int ui=0; ui<iel; ++ui)
      {
        double derinvJ0 = -v*(deriv_(0,ui)*xji_00 + deriv_(1,ui)*xji_01 + deriv_(2,ui)*xji_02);
        double derinvJ1 = -v*(deriv_(0,ui)*xji_10 + deriv_(1,ui)*xji_11 + deriv_(2,ui)*xji_12);
        double derinvJ2 = -v*(deriv_(0,ui)*xji_20 + deriv_(1,ui)*xji_21 + deriv_(2,ui)*xji_22);
        for (int vi=0; vi<iel; ++vi)
        {
          double visres0 =   2.0*derxy_(0, vi)* vderxyaf_(0, 0)
                             +     derxy_(1, vi)*(vderxyaf_(0, 1) + vderxyaf_(1, 0))
                             +     derxy_(2, vi)*(vderxyaf_(0, 2) + vderxyaf_(2, 0)) ;
          double visres1 =         derxy_(0, vi)*(vderxyaf_(0, 1) + vderxyaf_(1, 0))
                             + 2.0*derxy_(1, vi)* vderxyaf_(1, 1)
                             +     derxy_(2, vi)*(vderxyaf_(1, 2) + vderxyaf_(2, 1)) ;
          double visres2 =         derxy_(0, vi)*(vderxyaf_(0, 2) + vderxyaf_(2, 0))
                             +     derxy_(1, vi)*(vderxyaf_(1, 2) + vderxyaf_(2, 1))
                             + 2.0*derxy_(2, vi)* vderxyaf_(2, 2) ;
          meshmat(vi*4 + 0, ui*4 + 0) += derinvJ0*visres0;
          meshmat(vi*4 + 1, ui*4 + 0) += derinvJ0*visres1;
          meshmat(vi*4 + 2, ui*4 + 0) += derinvJ0*visres2;

          meshmat(vi*4 + 0, ui*4 + 1) += derinvJ1*visres0;
          meshmat(vi*4 + 1, ui*4 + 1) += derinvJ1*visres1;
          meshmat(vi*4 + 2, ui*4 + 1) += derinvJ1*visres2;

          meshmat(vi*4 + 0, ui*4 + 2) += derinvJ2*visres0;
          meshmat(vi*4 + 1, ui*4 + 2) += derinvJ2*visres1;
          meshmat(vi*4 + 2, ui*4 + 2) += derinvJ2*visres2;
        }
      }

      // part 2: derivative of viscosity residual

      v = fac*visceff/det_;
      for (int ui=0; ui<iel; ++ui)
      {
        double v0 = - vderiv_(0,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_100(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_100(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_100(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_200(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_200(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_200(ui)*xji_02);
        double v1 = - vderiv_(0,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_110(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_110(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_110(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_210(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_210(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_210(ui)*xji_02);
        double v2 = - vderiv_(0,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_120(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_120(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_120(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_220(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_220(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_220(ui)*xji_02);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 0, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(2*derxjm_001(ui)*xji_00 + 2*derxjm_001(ui)*xji_00 + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_10)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_10)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_10)
             - vderiv_(2,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20);
        v1 = - vderiv_(0,0)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_01 + 2*derxjm_011(ui)*xji_01 + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_11)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_11)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_11)
             - vderiv_(2,0)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21);
        v2 = - vderiv_(0,0)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_021(ui)*xji_02 + 2*derxjm_021(ui)*xji_02 + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_12)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_12)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_12)
             - vderiv_(2,0)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 0, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(2*derxjm_002(ui)*xji_00 + 2*derxjm_002(ui)*xji_00 + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_20);
        v1 = - vderiv_(0,0)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_01 + 2*derxjm_012(ui)*xji_01 + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_21);
        v2 = - vderiv_(0,0)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_022(ui)*xji_02 + 2*derxjm_022(ui)*xji_02 + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_22);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 0, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_100(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_00)
             - vderiv_(1,0)*(2*xji_10*derxjm_100(ui) + 2*xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_100(ui) + 2*xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_100(ui) + 2*xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20);
        v1 = - vderiv_(0,0)*(derxjm_100(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_01)
             - vderiv_(1,0)*(2*xji_10*derxjm_110(ui) + 2*xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_110(ui) + 2*xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_110(ui) + 2*xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21);
        v2 = - vderiv_(0,0)*(derxjm_100(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_02)
             - vderiv_(1,0)*(2*xji_10*derxjm_120(ui) + 2*xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_120(ui) + 2*xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_120(ui) + 2*xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 1, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_001(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_001(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_001(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_201(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_201(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_201(ui)*xji_12);
        v1 = - vderiv_(0,0)*(derxjm_011(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_011(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_011(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + xji_20*derxjm_211(ui) + xji_21*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + xji_22*derxjm_211(ui) + xji_21*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_211(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_211(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_211(ui)*xji_12);
        v2 = - vderiv_(0,0)*(derxjm_021(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_021(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_021(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + xji_20*derxjm_221(ui) + xji_22*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_221(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_221(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_221(ui)*xji_12);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 1, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
             - vderiv_(1,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + 2*xji_10*derxjm_102(ui) + 2*xji_10*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + 2*xji_11*derxjm_102(ui) + 2*xji_10*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + 2*xji_12*derxjm_102(ui) + 2*xji_10*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_20);
        v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
             - vderiv_(1,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + 2*xji_10*derxjm_112(ui) + 2*xji_11*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + 2*xji_11*derxjm_112(ui) + 2*xji_11*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + 2*xji_12*derxjm_112(ui) + 2*xji_11*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_21);
        v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
             - vderiv_(1,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + 2*xji_10*derxjm_122(ui) + 2*xji_12*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + 2*xji_11*derxjm_122(ui) + 2*xji_12*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + 2*xji_12*derxjm_122(ui) + 2*xji_12*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_22);

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 1, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_200(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_00)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + 2*xji_20*derxjm_200(ui) + 2*xji_20*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + 2*xji_21*derxjm_200(ui) + 2*xji_20*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + 2*xji_22*derxjm_200(ui) + 2*xji_20*derxjm_220(ui));
        v1 = - vderiv_(0,0)*(derxjm_200(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_01)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + 2*xji_20*derxjm_210(ui) + 2*xji_21*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + 2*xji_21*derxjm_210(ui) + 2*xji_21*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + 2*xji_22*derxjm_210(ui) + 2*xji_21*derxjm_220(ui));
        v2 = - vderiv_(0,0)*(derxjm_200(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_02)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + 2*xji_20*derxjm_220(ui) + 2*xji_22*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + 2*xji_21*derxjm_220(ui) + 2*xji_22*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + 2*xji_22*derxjm_220(ui) + 2*xji_22*derxjm_220(ui));

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 2, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_10)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_10)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_10)
             - vderiv_(2,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + 2*xji_20*derxjm_201(ui) + 2*xji_20*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + 2*xji_21*derxjm_201(ui) + 2*xji_20*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + 2*xji_22*derxjm_201(ui) + 2*xji_20*derxjm_221(ui));
        v1 = - vderiv_(0,0)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_11)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_11)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_11)
             - vderiv_(2,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + 2*xji_20*derxjm_211(ui) + 2*xji_21*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + 2*xji_21*derxjm_211(ui) + 2*xji_21*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + 2*xji_22*derxjm_211(ui) + 2*xji_21*derxjm_221(ui));
        v2 = - vderiv_(0,0)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_12)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_12)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_12)
             - vderiv_(2,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + 2*xji_20*derxjm_221(ui) + 2*xji_22*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + 2*xji_21*derxjm_221(ui) + 2*xji_22*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + 2*xji_22*derxjm_221(ui) + 2*xji_22*derxjm_221(ui));

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 2, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_002(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_002(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_102(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_102(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_102(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui));
        v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_012(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_012(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_112(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_112(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_112(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + xji_10*derxjm_112(ui) + xji_11*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + xji_12*derxjm_112(ui) + xji_11*derxjm_122(ui));
        v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_022(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_022(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_122(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_122(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_122(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + xji_10*derxjm_122(ui) + xji_12*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui));

        for (int vi=0; vi<iel; ++vi)
        {
          meshmat(vi*4 + 2, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }
      }


      // pressure
      for (int vi=0; vi<iel; ++vi)
      {
	const int fvi  =4*vi;
	const int fvip =fvi+1;
	const int fvipp=fvip+1;

        double v = prenp_*fac/det_;
        for (int ui=0; ui<iel; ++ui)
        {
	  const int fui  =4*ui;
	  const int fuip =fui+1;
	  const int fuipp=fuip+1;

          meshmat(fvi  ,fuip ) += v*(deriv_(0, vi)*derxjm_(0,0,1,ui)
				     +
				     deriv_(1, vi)*derxjm_(0,1,1,ui)
				     +
				     deriv_(2, vi)*derxjm_(0,2,1,ui)) ;
          meshmat(fvi  ,fuipp) += v*(deriv_(0, vi)*derxjm_(0,0,2,ui)
				     +
				     deriv_(1, vi)*derxjm_(0,1,2,ui)
				     +
				     deriv_(2, vi)*derxjm_(0,2,2,ui)) ;

          meshmat(fvip ,fui  ) += v*(deriv_(0, vi)*derxjm_(1,0,0,ui)
				     +
				     deriv_(1, vi)*derxjm_(1,1,0,ui)
				     +
				     deriv_(2, vi)*derxjm_(1,2,0,ui)) ;
          meshmat(fvip ,fuipp) += v*(deriv_(0, vi)*derxjm_(1,0,2,ui)
				     +
				     deriv_(1, vi)*derxjm_(1,1,2,ui)
				     +
				     deriv_(2, vi)*derxjm_(1,2,2,ui)) ;

          meshmat(fvipp,fui  ) += v*(deriv_(0, vi)*derxjm_(2,0,0,ui)
				     +
				     deriv_(1, vi)*derxjm_(2,1,0,ui)
				     +
				     deriv_(2, vi)*derxjm_(2,2,0,ui)) ;
          meshmat(fvipp,fuip ) += v*(deriv_(0, vi)*derxjm_(2,0,1,ui)
				     +
				     deriv_(1, vi)*derxjm_(2,1,1,ui)
				     +
				     deriv_(2, vi)*derxjm_(2,2,1,ui)) ;
        }
      }

      // div u
      for (int vi=0; vi<iel; ++vi)
      {
	const int fvippp = 4*vi+3;

        double v = fac/det_*funct_(vi);
        for (int ui=0; ui<iel; ++ui)
        {
	  const int fui   = 4*ui;
	  const int fuip  = fui+1;
	  const int fuipp = fuip+1;

          meshmat(fvippp,fui  ) += v*(vderiv_(1, 0)*derxjm_(0,0,1,ui)
				      +
				      vderiv_(1, 1)*derxjm_(0,1,1,ui)
				      +
				      vderiv_(1, 2)*derxjm_(0,2,1,ui)
				      +
				      vderiv_(2, 0)*derxjm_(0,0,2,ui)
				      +
				      vderiv_(2, 1)*derxjm_(0,1,2,ui)
				      +
				      vderiv_(2, 2)*derxjm_(0,2,2,ui));

          meshmat(fvippp,fuip ) += v*(vderiv_(0, 0)*derxjm_(1,0,0,ui)
				      +
				      vderiv_(0, 1)*derxjm_(1,1,0,ui)
				      +
				      vderiv_(0, 2)*derxjm_(1,2,0,ui)
				      +
				      vderiv_(2, 0)*derxjm_(1,0,2,ui)
				      +
				      vderiv_(2, 1)*derxjm_(1,1,2,ui)
				      +
				      vderiv_(2, 2)*derxjm_(1,2,2,ui));

          meshmat(fvippp,fuipp) += v*(vderiv_(0, 0)*derxjm_(2,0,0,ui)
				      +
				      vderiv_(0, 1)*derxjm_(2,1,0,ui)
				      +
				      vderiv_(0, 2)*derxjm_(2,2,0,ui)
				      +
				      vderiv_(1, 0)*derxjm_(2,0,1,ui)
				      +
				      vderiv_(1, 1)*derxjm_(2,1,1,ui)
				      +
				      vderiv_(1, 2)*derxjm_(2,2,1,ui));
        }
      }


    } // end linearisation with respect to mesh motion

  } // end loop iquad
  return;
} //Sysmat_adv_qs


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  advective version using time dependent subgrid scales              |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_adv_td(
  Fluid*                                    ele             ,
  std::vector<Epetra_SerialDenseVector>&     myknots         ,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>& elemat          ,
  LINALG::Matrix<(nsd_+1)*iel,1>&            elevec          ,
  const LINALG::Matrix<nsd_,iel>&            edispnp         ,
  const LINALG::Matrix<nsd_,iel>&            egridvelaf        ,
  const LINALG::Matrix<nsd_,iel>&            evelnp          ,
  const LINALG::Matrix<iel,1>&               eprenp          ,
  const LINALG::Matrix<nsd_,iel>&            eaccam          ,
  const LINALG::Matrix<nsd_,iel>&            evelaf          ,
  const LINALG::Matrix<nsd_,iel>&            fsevelaf        ,
  Teuchos::RCP<const MAT::Material>          material        ,
  const double                               alphaM          ,
  const double                               alphaF          ,
  const double                               gamma           ,
  const double                               dt              ,
  const double                               time            ,
  const enum INPAR::FLUID::LinearisationAction     newton          ,
  const bool                                 higher_order_ele,
  const enum INPAR::FLUID::FineSubgridVisc         fssgv           ,
  const enum INPAR::FLUID::Transient         inertia         ,
  const enum INPAR::FLUID::PSPG              pspg            ,
  const enum INPAR::FLUID::SUPG              supg            ,
  const enum INPAR::FLUID::VStab             vstab           ,
  const enum INPAR::FLUID::CStab             cstab           ,
  const enum INPAR::FLUID::CrossStress       cross           ,
  const enum INPAR::FLUID::ReynoldsStress    reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha                 whichtau        ,
  const enum INPAR::FLUID::TurbModelAction         turb_mod_action ,
  double&                                    Cs              ,
  double&                                    Cs_delta_sq     ,
  double&                                    visceff         ,
  const double                               l_tau           ,
  const bool                                 compute_elemat
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
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero

//  ele->ActivateTDS(intpoints.IP().nquad,nsd_);

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
                            egridvelaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,INPAR::FLUID::subscales_time_dependent,gamma,dt,hk,mk,visceff);

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
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);

        for (int mm=1;mm<3;++mm)
        {
          conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
        }
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

        if(higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      if (newton==INPAR::FLUID::Newton) // if newton and supg
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

      if(higher_order_ele && newton!=INPAR::FLUID::minimal)
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

        if (newton==INPAR::FLUID::Newton) // if pspg and newton
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
        if(vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_gls)
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

          if (newton==INPAR::FLUID::Newton)
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

        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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
    if(inertia == INPAR::FLUID::inertia_stab_keep
       ||
       inertia == INPAR::FLUID::inertia_stab_keep_complete)
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
    if(pspg == INPAR::FLUID::pstab_use_pspg)
    {

      const double fac_svelnpx                      = fac*svelnp_(0);
      const double fac_svelnpy                      = fac*svelnp_(1);
      const double fac_svelnpz                      = fac*svelnp_(2);

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
    if(supg == INPAR::FLUID::convective_stab_supg)
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
      if (vstab != INPAR::FLUID::viscous_stab_none)
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
    if(cross == INPAR::FLUID::cross_stress_stab_only_rhs || cross == INPAR::FLUID::cross_stress_stab)
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
    if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_cons_qs(
  Fluid*                                             ele             ,
  std::vector<Epetra_SerialDenseVector>&              myknots         ,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&          elemat          ,
  LINALG::Matrix<(nsd_+1)*iel,1>&                     elevec          ,
  const LINALG::Matrix<nsd_,iel>&                     edispnp         ,
  const LINALG::Matrix<nsd_,iel>&                     egridvelaf      ,
  const LINALG::Matrix<nsd_,iel>&                     evelnp          ,
  const LINALG::Matrix<iel,1>&                        eprenp          ,
  const LINALG::Matrix<nsd_,iel>&                     eaccam          ,
  const LINALG::Matrix<nsd_,iel>&                     evelaf          ,
  const LINALG::Matrix<nsd_,iel>&                     fsevelaf        ,
  Teuchos::RCP<const MAT::Material>                   material        ,
  const double                                        alphaM          ,
  const double                                        alphaF          ,
  const double                                        gamma           ,
  const double                                        dt              ,
  const double                                        time            ,
  const enum INPAR::FLUID::LinearisationAction              newton          ,
  const bool                                          higher_order_ele,
  const enum INPAR::FLUID::FineSubgridVisc                  fssgv           ,
  const enum INPAR::FLUID::PSPG                       pspg            ,
  const enum INPAR::FLUID::SUPG                       supg            ,
  const enum INPAR::FLUID::VStab                      vstab           ,
  const enum INPAR::FLUID::CStab                      cstab           ,
  const enum INPAR::FLUID::CrossStress                cross           ,
  const enum INPAR::FLUID::ReynoldsStress             reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha                          whichtau        ,
  const enum INPAR::FLUID::TurbModelAction                  turb_mod_action ,
  double&                                             Cs              ,
  double&                                             Cs_delta_sq     ,
  double&                                             visceff         ,
  const double                                        l_tau           ,
  const bool                                          compute_elemat
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
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
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
                            egridvelaf        ,
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
    if(cross == INPAR::FLUID::cross_stress_stab)
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
    if(ele->IsAle())
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
    CalcTau(whichtau,INPAR::FLUID::subscales_quasistatic,gamma,dt,hk,mk,visceff);

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

      if(reynolds == INPAR::FLUID::reynolds_stress_stab)
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


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      if (newton==INPAR::FLUID::Newton)
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

	if (newton==INPAR::FLUID::Newton)
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
	if((vstab == INPAR::FLUID::viscous_stab_gls || vstab == INPAR::FLUID::viscous_stab_usfem)&&higher_order_ele)
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

	  if (newton==INPAR::FLUID::Newton)
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
      if(cross == INPAR::FLUID::cross_stress_stab)
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

	if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
      {
	elevec(4*ui+3)-=fac_tauMp*(resM_(0)*derxy_(0,ui)+resM_(1)*derxy_(1,ui)+resM_(2)*derxy_(2,ui));
      }
    } // end pspg

    if(supg == INPAR::FLUID::convective_stab_supg)
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
      if(vstab != INPAR::FLUID::viscous_stab_none && higher_order_ele)
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

    if(fssgv != INPAR::FLUID::no_fssgv)
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
    } // end not no_fssgv
  } // end loop iquad
  return;
} // end Sysmat_cons_qs


/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |  conservative, time-dependent version                               |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::Sysmat_cons_td(
  Fluid*                                          ele             ,
  std::vector<Epetra_SerialDenseVector>&           myknots         ,
  LINALG::Matrix<(nsd_+1)*iel,(nsd_+1)*iel>&       elemat          ,
  LINALG::Matrix<(nsd_+1)*iel,1>&                  elevec          ,
  const LINALG::Matrix<nsd_,iel>&                  edispnp         ,
  const LINALG::Matrix<nsd_,iel>&                  egridvelaf        ,
  const LINALG::Matrix<nsd_,iel>&                  evelnp          ,
  const LINALG::Matrix<iel,1>&                     eprenp          ,
  const LINALG::Matrix<nsd_,iel>&                  eaccam          ,
  const LINALG::Matrix<nsd_,iel>&                  evelaf          ,
  const LINALG::Matrix<nsd_,iel>&                  fsevelaf        ,
  Teuchos::RCP<const MAT::Material>                material        ,
  const double                                     alphaM          ,
  const double                                     alphaF          ,
  const double                                     gamma           ,
  const double                                     dt              ,
  const double                                     time            ,
  const enum INPAR::FLUID::LinearisationAction           newton          ,
  const bool                                       higher_order_ele,
  const enum INPAR::FLUID::FineSubgridVisc               fssgv           ,
  const enum INPAR::FLUID::Transient               inertia         ,
  const enum INPAR::FLUID::PSPG                    pspg            ,
  const enum INPAR::FLUID::SUPG                    supg            ,
  const enum INPAR::FLUID::VStab                   vstab           ,
  const enum INPAR::FLUID::CStab                   cstab           ,
  const enum INPAR::FLUID::CrossStress             cross           ,
  const enum INPAR::FLUID::ReynoldsStress          reynolds        ,
  const enum INPAR::FLUID::TauType_genalpha                       whichtau        ,
  const enum INPAR::FLUID::TurbModelAction               turb_mod_action ,
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
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero

//  ele->ActivateTDS(intpoints.IP().nquad,nsd_);

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
                            egridvelaf        ,
                            evelnp          ,
                            eprenp          ,
                            eaccam          ,
                            evelaf          ,
                            fsevelaf        ,
                            visceff         ,
                            fssgv           ,
                            higher_order_ele);

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,INPAR::FLUID::subscales_time_dependent,gamma,dt,hk,mk,visceff);

    // stabilisation parameters
    const double tauM   = tau_(0);

    if(cstab == INPAR::FLUID::continuity_stab_none)
    {
      tau_(2)=0.0;
    }
    const double tauC    = tau_(2);

    double supg_active;
    if(supg == INPAR::FLUID::convective_stab_supg)
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
        conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);

        for (int mm=1;mm<3;++mm)
        {
          conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
        }
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
    if(ele->IsAle())
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
    //     ELEMENT FORMULATION BASED ON DYNAMIC SUBSCALES
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
        if(newton==INPAR::FLUID::Newton)
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

        if(higher_order_ele && newton!=INPAR::FLUID::minimal)
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

      if (newton==INPAR::FLUID::Newton) // if newton and supg
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

      if(higher_order_ele && newton!=INPAR::FLUID::minimal)
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


        if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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
      if(cross == INPAR::FLUID::cross_stress_stab)
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

	if (higher_order_ele && newton!=INPAR::FLUID::minimal)
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

        if (newton==INPAR::FLUID::Newton) // if pspg and newton
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
        if(vstab == INPAR::FLUID::viscous_stab_usfem || vstab == INPAR::FLUID::viscous_stab_gls)
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

          if (newton==INPAR::FLUID::Newton)
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
    if(inertia == INPAR::FLUID::inertia_stab_keep
       ||
       inertia == INPAR::FLUID::inertia_stab_keep_complete)
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

    if(supg == INPAR::FLUID::convective_stab_supg)
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
    if(pspg == INPAR::FLUID::pstab_use_pspg)
    {

      const double fac_svelnpx                      = fac*svelnp_(0);
      const double fac_svelnpy                      = fac*svelnp_(1);
      const double fac_svelnpz                      = fac*svelnp_(2);

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
      if (vstab != INPAR::FLUID::viscous_stab_none)
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

    if(fssgv != INPAR::FLUID::no_fssgv)
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
    } // end not no_fssgv
  } // end loop iquad
  return;
} // end Sysmat_cons_td


/*----------------------------------------------------------------------*
 |  evaluate residuals resM and resC on the element and return          |
 |  averaged values together with averaged subgrid-scale quantities.    |
 |                                                (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalcResAvgs(
  Fluid*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  RefCountPtr<MAT::Material> mat)
{

  // --------------------------------------------------
  // create matrix objects for nodal values
  LINALG::Matrix<iel,1> eprenp    ;
  LINALG::Matrix<nsd_,iel> evelnp    ;
  LINALG::Matrix<nsd_,iel> evelaf    ;
  LINALG::Matrix<nsd_,iel> eaccam    ;
  LINALG::Matrix<nsd_,iel> edispnp   ;
  LINALG::Matrix<nsd_,iel> egridvelaf;
  LINALG::Matrix<nsd_,iel> fsevelaf  ;

  // --------------------------------------------------
  // set parameters for time integration
  ParameterList& timelist = params.sublist("time integration parameters");

  const double alphaF = timelist.get<double>("alpha_F");
  const double gamma  = timelist.get<double>("gamma");
  const double dt     = timelist.get<double>("dt");
  const double time   = timelist.get<double>("time");

  // --------------------------------------------------
  // disable fine-scale subgrid viscosity
  INPAR::FLUID::FineSubgridVisc fssgv = INPAR::FLUID::no_fssgv;

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");

  // specify which residual based stabilisation terms
  // will be used
  INPAR::FLUID::SubscalesTD             tds = DRT::INPUT::IntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
  INPAR::FLUID::Transient       inertia = DRT::INPUT::IntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
  INPAR::FLUID::PSPG            pspg = DRT::INPUT::IntegralValue<INPAR::FLUID::PSPG>(stablist,"PSPG");
  INPAR::FLUID::SUPG            supg = DRT::INPUT::IntegralValue<INPAR::FLUID::SUPG>(stablist,"SUPG");
  INPAR::FLUID::VStab           vstab = DRT::INPUT::IntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
  INPAR::FLUID::CStab           cstab = DRT::INPUT::IntegralValue<INPAR::FLUID::CStab>(stablist,"CSTAB");
  INPAR::FLUID::CrossStress     cross = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
  INPAR::FLUID::ReynoldsStress  reynolds = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");


  // flag conservative form on/off
  string conservativestr =params.get<string>("CONVFORM");

  double fac_conservative =1.0;
  if(conservativestr=="conservative")
  {
    fac_conservative =-1.0;
  }


  // select tau definition
  INPAR::FLUID::TauType_genalpha whichtau = INPAR::FLUID::tautype_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");

    if(taudef == "Franca_Barrenechea_Valentin_Frey_Wall")
    {
      whichtau = INPAR::FLUID::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Taylor_Hughes_Zarins_Whiting_Jansen")
    {
      whichtau = INPAR::FLUID::taylor_hughes_zarins_whiting_jansen;
    }
    else if(taudef == "Codina")
    {
      whichtau = INPAR::FLUID::codina;
    }
    else if(taudef == "Smoothed_FBVW")
    {
      whichtau = INPAR::FLUID::smoothed_franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt")
    {
      whichtau = INPAR::FLUID::fbvw_wo_dt;
    }
    else if(taudef == "BFVW_gradient_based_hk")
    {
      whichtau = INPAR::FLUID::fbvw_gradient_based_hk;
    }
    else
    {
      dserror("unknown tau definition\n");
    }
  }

  // flag for higher order elements
  //bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());
  bool higher_order_ele = IsHigherOrder<distype>::ishigherorder;

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
  ParameterList& sgviscparams    = params.sublist("SUBGRID VISCOSITY");

  // the default action is no model
  INPAR::FLUID::TurbModelAction turb_mod_action = INPAR::FLUID::no_model;

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
    sgviscparams   ,
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
        ele->IsAle()  ,
        discretization,
        lm,
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
  std::vector<Epetra_SerialDenseVector> myknots(nsd_);

  // for isogeometric elements
  if(ele->Shape()==Fluid::nurbs8 || ele->Shape()==Fluid::nurbs27)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    bool zero_size = false;
    zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

    // if we have a zero sized element due to a interpolated
    // point --- exit here
    if(zero_size)
    {
      return(0);
    }
  }
  else if (ele->Shape()==Fluid::nurbs4 || ele->Shape()==Fluid::nurbs9)
  {
    dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());
  }

  // visceff will contain the computed effective viscosity
  double visceff       = 0.0;

  // the coordinates of the element layers in the channel
  RefCountPtr<vector<double> > planecoords  = params.get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);

  if(planecoords==Teuchos::null)
    dserror("planecoords is null, but need channel_flow_of_height_2\n");

  // ---------------------------------------------------
  // working arrays for the quantities we want to compute
  LINALG::Matrix<nsd_,1>  mean_res        ;
  LINALG::Matrix<nsd_,1>  mean_sacc       ;
  LINALG::Matrix<nsd_,1>  mean_svelaf     ;
  LINALG::Matrix<nsd_,1>  mean_res_sq     ;
  LINALG::Matrix<nsd_,1>  mean_sacc_sq    ;
  LINALG::Matrix<nsd_,1>  mean_svelaf_sq  ;
  LINALG::Matrix<nsd_,1>  mean_tauinvsvel ;

  LINALG::Matrix<2*nsd_,1>  mean_crossstress;
  LINALG::Matrix<2*nsd_,1>  mean_reystress  ;

  double vol             = 0.0;

  double h               = 0.0;
  double h_bazilevs      = 0.0;
  double strle           = 0.0;
  double gradle          = 0.0;
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


  mean_res        .Clear();
  mean_sacc       .Clear();
  mean_svelaf     .Clear();
  mean_res_sq     .Clear();
  mean_sacc_sq    .Clear();
  mean_svelaf_sq  .Clear();
  mean_tauinvsvel .Clear();
  mean_crossstress.Clear();
  mean_reystress  .Clear();

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
                 mat            ,
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
  // const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // TODO
  // the integration loop is only a 3D function, since it may influence the performance

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

    if(tds == INPAR::FLUID::subscales_time_dependent)
    {
      /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

      */
//      for (int rr=0;rr<3;++rr)
//      {
//        svelaf_(rr) =
//          alphaF      *(ele->Svelnp())(rr,iquad)
//          +
//          (1.0-alphaF)*(ele->Sveln())(rr,iquad);
//      }
    }

    /*---------------------------- get stabilisation parameter ---*/
    CalcTau(whichtau,tds,gamma,dt,hk,mk,visceff);

    // -----------------------------------------------------
    // Element volume
    vol               += fac;

    // volume based element size
    h                 += fac*hk;

    // streamlength based element size
    {
      const double vel_normaf=velintaf_.Norm2();

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
        val += fabs( normed_velintaf(0)*derxy_(0,rr)
                    +normed_velintaf(1)*derxy_(1,rr)
                    +normed_velintaf(2)*derxy_(2,rr));
      } /* end of loop over element nodes */
      strle += 2.0/val*fac;
    }

    // element size in main gradient direction
    {
      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velgrad;

      for (int rr=0;rr<3;++rr)
      {
        normed_velgrad(rr)=sqrt(vderxyaf_(0,rr)*vderxyaf_(0,rr)
                                +
                                vderxyaf_(1,rr)*vderxyaf_(1,rr)
                                +
                                vderxyaf_(2,rr)*vderxyaf_(2,rr));
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<3;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<3;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        val += fabs( normed_velgrad(0)*derxy_(0,rr)
                    +normed_velgrad(1)*derxy_(1,rr)
                    +normed_velgrad(2)*derxy_(2,rr));
      } /* end of loop over element nodes */
      gradle += 2.0/val*fac;
    }

    // -----------------------------------------------------
    {
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

      h_bazilevs+=1./sqrt(sqrt(normG))*fac;
     }

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

    if(tds == INPAR::FLUID::subscales_time_dependent
       &&
       (inertia == INPAR::FLUID::inertia_stab_keep
        ||
        inertia == INPAR::FLUID::inertia_stab_keep_complete))
    {
      for(int rr=0;rr<3;++rr)
      {
        const double aux = -1.0/tau_(1)*svelaf_(rr) -resM_(rr);

        mean_sacc   (rr) += aux*fac;
        mean_sacc_sq(rr) += aux*aux*fac;
      }
      const double aux0 = -1.0/tau_(1)*svelaf_(0)-resM_(0);
      const double aux1 = -1.0/tau_(1)*svelaf_(1)-resM_(1);
      const double aux2 = -1.0/tau_(1)*svelaf_(2)-resM_(2);

      abs_sacc += sqrt(aux0*aux0+aux1*aux1+aux2*aux2)*fac;
    }

    if(tds == INPAR::FLUID::subscales_time_dependent)
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
        const double aux = tau_(1)*resM_(rr);

        mean_svelaf   (rr) -= aux*fac;
        mean_svelaf_sq(rr) += aux*aux*fac;
      }

      abs_svel +=sqrt(resM_(0)*resM_(0)
                      +
                      resM_(1)*resM_(1)
                      +
                      resM_(2)*resM_(2))*tau_(1)*fac;
    }

    for(int rr=0;rr<3;++rr)
    {
      mean_tauinvsvel(rr)+=mean_svelaf(rr)/tau_(1);
    }


    {
      const double aux = tau_(2)*divunp_;

      mean_sprenp     -= aux*fac;
      mean_sprenp_sq  += aux*aux*fac;
    }

    // ------------------------------------------------------
    // Element average:
    //      o cross stresses 11,22,33,12,23,31
    //      o reynolds stresses 11,22,33,12,23,31
    if(tds      == INPAR::FLUID::subscales_time_dependent)
    {
      if(cross    != INPAR::FLUID::cross_stress_stab_none)
      {
        mean_crossstress(0)+=fac*(svelaf_(0)*velintaf_(0)+svelaf_(0)*velintaf_(0));
        mean_crossstress(1)+=fac*(svelaf_(1)*velintaf_(1)+svelaf_(1)*velintaf_(1));
        mean_crossstress(2)+=fac*(svelaf_(2)*velintaf_(2)+svelaf_(2)*velintaf_(2));
        mean_crossstress(3)+=fac*(svelaf_(0)*velintaf_(1)+svelaf_(1)*velintaf_(0));
        mean_crossstress(4)+=fac*(svelaf_(1)*velintaf_(2)+svelaf_(2)*velintaf_(1));
        mean_crossstress(5)+=fac*(svelaf_(2)*velintaf_(0)+svelaf_(0)*velintaf_(2));
      }

      if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
      {
        mean_reystress(0)  +=fac*(svelaf_(0)*svelaf_(0)+svelaf_(0)*svelaf_(0));
        mean_reystress(1)  +=fac*(svelaf_(1)*svelaf_(1)+svelaf_(1)*svelaf_(1));
        mean_reystress(2)  +=fac*(svelaf_(2)*svelaf_(2)+svelaf_(2)*svelaf_(2));
        mean_reystress(3)  +=fac*(svelaf_(0)*svelaf_(1)+svelaf_(1)*svelaf_(0));
        mean_reystress(4)  +=fac*(svelaf_(1)*svelaf_(2)+svelaf_(2)*svelaf_(1));
        mean_reystress(5)  +=fac*(svelaf_(2)*svelaf_(0)+svelaf_(0)*svelaf_(2));
      }
    }
    else
    {
      if(cross    != INPAR::FLUID::cross_stress_stab_none)
      {
        mean_crossstress(0)-=fac*tau_(1)*(resM_(0)*velintaf_(0)+velintaf_(0)*resM_(0));
        mean_crossstress(1)-=fac*tau_(1)*(resM_(1)*velintaf_(1)+velintaf_(1)*resM_(1));
        mean_crossstress(2)-=fac*tau_(1)*(resM_(2)*velintaf_(2)+velintaf_(2)*resM_(2));
        mean_crossstress(3)-=fac*tau_(1)*(resM_(0)*velintaf_(1)+velintaf_(0)*resM_(1));
        mean_crossstress(4)-=fac*tau_(1)*(resM_(1)*velintaf_(2)+velintaf_(1)*resM_(2));
        mean_crossstress(5)-=fac*tau_(1)*(resM_(2)*velintaf_(0)+velintaf_(2)*resM_(0));
      }

      if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
      {
        mean_reystress(0)  +=fac*tau_(1)*tau_(1)*(resM_(0)*resM_(0)+resM_(0)*resM_(0));
        mean_reystress(1)  +=fac*tau_(1)*tau_(1)*(resM_(1)*resM_(1)+resM_(1)*resM_(1));
        mean_reystress(2)  +=fac*tau_(1)*tau_(1)*(resM_(2)*resM_(2)+resM_(2)*resM_(2));
        mean_reystress(3)  +=fac*tau_(1)*tau_(1)*(resM_(0)*resM_(1)+resM_(1)*resM_(0));
        mean_reystress(4)  +=fac*tau_(1)*tau_(1)*(resM_(1)*resM_(2)+resM_(2)*resM_(1));
        mean_reystress(5)  +=fac*tau_(1)*tau_(1)*(resM_(2)*resM_(0)+resM_(0)*resM_(2));
      }
    }


    // ------------------------------------------------------
    // Element average:
    //      o tauM
    //      o tauC

    averaged_tauM+=tau_(1)*fac;
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

    if(tds      == INPAR::FLUID::subscales_time_dependent
       &&
       (inertia  == INPAR::FLUID::inertia_stab_keep
        ||
        inertia == INPAR::FLUID::inertia_stab_keep_complete))
    {
      double sacc[3];
      sacc[0]= -1.0/tau_(1)*svelaf_(0) -resM_(0);
      sacc[1]= -1.0/tau_(1)*svelaf_(1) -resM_(1);
      sacc[2]= -1.0/tau_(1)*svelaf_(2) -resM_(2);

      /*

               /                 \
              |   ~ n+am    n+af  |
              |  acc     , u      |
              |     (i)     (i)   |
               \                 /

      */
      eps_sacc+=(sacc[0]*velintaf_(0)+sacc[1]*velintaf_(1)+sacc[2]*velintaf_(2))*fac;
    }

    if(pspg     == INPAR::FLUID::pstab_use_pspg)
    {
      // contribution of this gausspoint to energy dissipated by pspg
      /*
                        /                  \
                       |  ^ n+1        n+1  |
                     - |  u   , nabla p     |
                       |                    |
                        \                  /
      */
//      if(tds == INPAR::FLUID::subscales_time_dependent)
//      {
//        eps_pspg-=
//          ((ele->Svelnp())(0,iquad)*pderxynp_(0)
//           +
//           (ele->Svelnp())(1,iquad)*pderxynp_(1)
//           +
//           (ele->Svelnp())(2,iquad)*pderxynp_(2))*fac;
//      }
//      else
//      {
        eps_pspg+=tau_(1)*fac*
          (resM_(0)*pderxynp_(0)
           +
           resM_(1)*pderxynp_(1)
           +
           resM_(2)*pderxynp_(2));
//      }
    }

    if(supg     == INPAR::FLUID::convective_stab_supg)
    {
      // contribution of this gausspoint to energy dissipated by supg
      /*

                   /                                      \
                  |  ~n+af    / n+af \           / n+af\   |
          - rho * |  u     , | u      | o nabla | u     |  |
                  |           \      /           \     /   |
                   \                                      /
      */
      if(tds == INPAR::FLUID::subscales_time_dependent)
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
        eps_supg+=tau_(1)*fac*
          (resM_(0)*convaf_old_(0)
           +
           resM_(1)*convaf_old_(1)
           +
           resM_(2)*convaf_old_(2));
      }
    }

    if(cross    != INPAR::FLUID::cross_stress_stab_none)
    {
      /*
               /                                       \
              |    n+af    /~n+af \           / n+af\   |
              |   u     , | u      | o nabla | u     |  |
              |            \      /           \     /   |
               \                                       /

      */
      if(tds == INPAR::FLUID::subscales_time_dependent)
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
                                      svelaf_(2)*vderxyaf_(2,2)))*fac_conservative;
      }
      else
      {
        eps_cross-=
          tau_(1)*fac*(velintaf_(0)*(resM_(0)*vderxyaf_(0,0)+
                                     resM_(1)*vderxyaf_(0,1)+
                                     resM_(2)*vderxyaf_(0,2))
                       +
                       velintaf_(1)*(resM_(0)*vderxyaf_(1,0)+
                                     resM_(1)*vderxyaf_(1,1)+
                                     resM_(2)*vderxyaf_(1,2))
                       +
                       velintaf_(2)*(resM_(0)*vderxyaf_(2,0)+
                                     resM_(1)*vderxyaf_(2,1)+
                                     resM_(2)*vderxyaf_(2,2)))*fac_conservative;
      }
    }

    if(reynolds != INPAR::FLUID::reynolds_stress_stab_none)
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
      if(tds == INPAR::FLUID::subscales_time_dependent)
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
          tau_(1)*tau_(1)*fac*
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

    if(vstab    != INPAR::FLUID::viscous_stab_none)
    {

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

      /*
         /                        \
        |  ~n+af                   |
        |  u      , 2*div eps (v)  |
        |                          |
         \                        /

      */
      if(tds == INPAR::FLUID::subscales_time_dependent)
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
          vstabfac*tau_(1)*fac*
          (resM_(0)*visceff*viscaf_old_(0)
           +
           resM_(1)*visceff*viscaf_old_(1)
           +
           resM_(2)*visceff*viscaf_old_(2));
      }
    }

    if(cstab    != INPAR::FLUID::continuity_stab_none)
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
                         velintaf_(2)*vderxyaf_(2,2)))*fac_conservative;
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
      if(turb_mod_action == INPAR::FLUID::dynamic_smagorinsky
         ||
         turb_mod_action == INPAR::FLUID::smagorinsky_with_van_Driest_damping
         ||
         turb_mod_action == INPAR::FLUID::smagorinsky)
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


  for(int rr=0;rr<6;++rr)
  {
    mean_crossstress(rr)/=vol;
    mean_reystress  (rr)/=vol;
  }

  abs_res         /= vol;
  abs_sacc        /= vol;
  abs_svel        /= vol;

  mean_resC       /= vol;
  mean_resC_sq    /= vol;
  mean_sprenp     /= vol;
  mean_sprenp_sq  /= vol;

  h               /= vol;
  h_bazilevs      /= vol;
  strle           /= vol;
  gradle          /= vol;

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

  RefCountPtr<vector<double> > incrhk            = params.get<RefCountPtr<vector<double> > >("incrhk"           );
  RefCountPtr<vector<double> > incrhbazilevs     = params.get<RefCountPtr<vector<double> > >("incrhbazilevs"    );
  RefCountPtr<vector<double> > incrstrle         = params.get<RefCountPtr<vector<double> > >("incrstrle"        );
  RefCountPtr<vector<double> > incrgradle        = params.get<RefCountPtr<vector<double> > >("incrgradle"       );

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

  RefCountPtr<vector<double> > incrcrossstress   = params.get<RefCountPtr<vector<double> > >("incrcrossstress"  );
  RefCountPtr<vector<double> > incrreystress     = params.get<RefCountPtr<vector<double> > >("incrreystress"    );


  if(not (distype == DRT::Element::nurbs8
       ||
       distype == DRT::Element::nurbs27))
  {

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
  }
  else if((distype == DRT::Element::nurbs8
       ||
       distype == DRT::Element::nurbs27))
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    bool zero_size = false;
    zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

    // if we have a zero sized element due to a interpolated
    // point --- exit here
    if(zero_size)
    {
      return(0);
    }

    int patchid=0;
    int gid    =ele->Id();
    vector<int> ele_cart_id(3);

    (nurbsdis->GetKnotVector())->ConvertEleGidToKnotIds(gid,patchid,ele_cart_id);

    nlayer=ele_cart_id[1];
  }
  else dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());

  // collect layer volume
  (*incrvol      )[nlayer] += vol;

  // element length in stabilisation parameter
  (*incrhk       )[nlayer] += hk;

  // element length in viscous regime defined by the Bazilevs parameter
  (*incrhbazilevs)[nlayer] += h_bazilevs;

  // stream length
  (*incrstrle    )[nlayer] += strle;

  // gradient based element length
  (*incrgradle   )[nlayer] += gradle;

  // averages of stabilisation parameters
  (*incrtauC     )[nlayer] += averaged_tauC;
  (*incrtauM     )[nlayer] += averaged_tauM;

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

  // averages of subgrid stress tensors
  for(int mm=0;mm<6;++mm)
  {
    (*incrcrossstress)[6*nlayer+mm] += mean_crossstress(mm);
    (*incrreystress  )[6*nlayer+mm] += mean_reystress  (mm);
  }

  return(0);
}

/*----------------------------------------------------------------------*
 |  extract velocities, pressure and accelerations from the global      |
 |  distributed vectors                           (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::ExtractValuesFromGlobalVectors(
        const INPAR::FLUID::FineSubgridVisc              fssgv         ,
        const bool                                 is_ale        ,
        const DRT::Discretization&                 discretization,
        const vector<int>&                         lm,
        LINALG::Matrix<iel,1>& eprenp        ,
        LINALG::Matrix<nsd_,iel>& evelnp        ,
        LINALG::Matrix<nsd_,iel>& evelaf        ,
        LINALG::Matrix<nsd_,iel>& eaccam        ,
        LINALG::Matrix<nsd_,iel>& edispnp       ,
        LINALG::Matrix<nsd_,iel>& egridvelaf    ,
        LINALG::Matrix<nsd_,iel>& fsevelaf
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
      /*
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
      */
      for (int idim=0;idim<nsd_;++idim)
      {
        eaccam(idim,i) = myaccam[idim+(numdofpernode_*i)];
        evelnp(idim,i) = myvelnp[idim+(numdofpernode_*i)];
        evelaf(idim,i) = myvelaf[idim+(numdofpernode_*i)];
      }
      eprenp(i)   = myvelnp[nsd_+(numdofpernode_*i)];
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
        for (int idim=0;idim<nsd_;++idim)
        {
          egridvelaf(idim,i) = mygridvelaf[idim+(numdofpernode_*i)];
          edispnp(idim,i)    = mydispnp[idim+(numdofpernode_*i)];
        }
        /*
        int fi    =4*i;
        int fip   =fi+1;
        int fipp  =fip+1;

        egridvelaf(0,i) = mygridvelaf[fi  ];
        egridvelaf(1,i) = mygridvelaf[fip ];
        egridvelaf(2,i) = mygridvelaf[fipp];

        edispnp(0,i)    = mydispnp   [fi  ];
        edispnp(1,i)    = mydispnp   [fip ];
        edispnp(2,i)    = mydispnp   [fipp];
        */
      }
    }

    // get fine-scale velocity
    if (fssgv != INPAR::FLUID::no_fssgv)
    {
      RCP<const Epetra_Vector> fsvelaf;

      fsvelaf = discretization.GetState("fsvelaf");
      if (fsvelaf==null) dserror("Cannot get state vector 'fsvelaf'");
      vector<double> myfsvelaf(lm.size());
      DRT::UTILS::ExtractMyValues(*fsvelaf,myfsvelaf,lm);

      // get fine-scale velocity and insert into element arrays
      for (int i=0;i<iel;++i)
      {
        for (int idim=0;idim<nsd_;++idim)
        {
          fsevelaf(idim,i) = myfsvelaf[idim+(numdofpernode_*i)];
        }
        /*
        fsevelaf(0,i) = myfsvelaf[0+(i*4)];
        fsevelaf(1,i) = myfsvelaf[1+(i*4)];
        fsevelaf(2,i) = myfsvelaf[2+(i*4)];
        */
      }
    }
    else
    {
      for (int i=0;i<iel;++i)
      {
        for (int idim=0;idim<nsd_;++idim)
        {
          fsevelaf(idim,i) = 0.0;
        }
        /*
        fsevelaf(0,i) = 0.0;
        fsevelaf(1,i) = 0.0;
        fsevelaf(2,i) = 0.0;
        */
      }
    }
  return;
} //ExtractValuesFromGlobalVectors


/*----------------------------------------------------------------------*
 |  Set parameters for turbulence models                                |
 |                                                (private) gammi 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::SetParametersForTurbulenceModel(
  const Fluid*                       ele            ,
  ParameterList                     & turbmodelparams,
  ParameterList                     & sgviscparams,
  const INPAR::FLUID::FineSubgridVisc     & fssgv          ,
  INPAR::FLUID::TurbModelAction           & turb_mod_action,
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
  if (fssgv != INPAR::FLUID::no_fssgv && turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
    dserror("No combination of a classical (all-scale) turbulence model and a fine-scale subgrid-viscosity approach currently possible!");
  if (fssgv != INPAR::FLUID::no_fssgv) Cs = sgviscparams.get<double>("C_SMAGORINSKY",0.0);

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action = INPAR::FLUID::smagorinsky;
      Cs              = sgviscparams.get<double>("C_SMAGORINSKY");
    }
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      RefCountPtr<vector<double> > planecoords      = turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);

      if(planecoords==Teuchos::null)
        dserror("planecoords is null, but need channel_flow_of_height_2\n");

      // for the Smagorinsky model with van Driest damping, we need a viscous length to determine
      // the y+ (heigth in wall units)
      turb_mod_action = INPAR::FLUID::smagorinsky_with_van_Driest_damping;
      Cs              = sgviscparams.get<double>("C_SMAGORINSKY");
      l_tau           = sgviscparams.get<double>("CHANNEL_L_TAU");

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
      turb_mod_action = INPAR::FLUID::dynamic_smagorinsky;

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
        Cs_delta_sq = 0.0; //ele->CsDeltaSq();
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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::GetNodalBodyForce(
  const Fluid* ele,
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
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

    // set this condition to the edeadng array
    for(int isd=0;isd<nsd_;isd++)
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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalVisc(
  Teuchos::RCP<const MAT::Material>       material   ,
  double&                           	  visc       ,
  const double &                          rateofshear
  )
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
} // DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalVisc


/*----------------------------------------------------------------------*
 |                                                                      |
 |  calculation of stabilisation parameter          gammi 12/08         |
 |    (element center or Gaussian point)                                |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalcTau(
  const enum INPAR::FLUID::TauType_genalpha             whichtau  ,
  const enum INPAR::FLUID::SubscalesTD           tds       ,
  const double &                         gamma     ,
  const double &                         dt        ,
  const double &                         hk        ,
  const double &                         mk        ,
  const double &                         visceff   )
{
  // TODO
  // This function is only a 3D function, since it may influence the performance
  // get velocity norms
  // For the time being, the ALE convective velocity, which appears to be the
  // more appropriate velocity, is used here. Both at n+1 and n+alpha_F, the
  // the ALE convective velocity at n+alpha_F is used, since the ALE convective
  // at n+1 is not available.
  const double vel_normaf = aleconvintaf_.Norm2();
  const double vel_normnp = aleconvintaf_.Norm2();

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
          Gnormu+=aleconvintaf_(nn)*G(nn,rr)*aleconvintaf_(rr);
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
      // this copy of aleconvintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velgrad;

      for (int rr=0;rr<3;++rr)
      {
        normed_velgrad(rr)=sqrt(vderxyaf_(0,rr)*vderxyaf_(0,rr)
                                +
                                vderxyaf_(1,rr)*vderxyaf_(1,rr)
                                +
                                vderxyaf_(2,rr)*vderxyaf_(2,rr));
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<3;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<3;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        val += fabs( normed_velgrad(0)*derxy_(0,rr)
                    +normed_velgrad(1)*derxy_(1,rr)
                    +normed_velgrad(2)*derxy_(2,rr));
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
          Gnormu+=aleconvintaf_(nn)*G(nn,rr)*aleconvintaf_(rr);
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
    else if (whichtau == INPAR::FLUID::franca_barrenechea_valentin_wall)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Wall


      // this copy of aleconvintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=aleconvintaf_(rr)/vel_normaf;
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
        val += fabs( normed_velintaf(0)*derxy_(0,rr)
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


      // this copy of aleconvintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velintaf;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_normaf>=1e-6)
      {
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          normed_velintaf(rr)=aleconvintaf_(rr)/vel_normaf;
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
        val += fabs( normed_velintaf(0)*derxy_(0,rr)
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
      // this copy of aleconvintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velgrad;

      for (int rr=0;rr<3;++rr)
      {
        normed_velgrad(rr)=sqrt(vderxyaf_(0,rr)*vderxyaf_(0,rr)
                                +
                                vderxyaf_(1,rr)*vderxyaf_(1,rr)
                                +
                                vderxyaf_(2,rr)*vderxyaf_(2,rr));
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<3;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<3;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;
      for (int rr=0;rr<iel;++rr) /* loop element nodes */
      {
        val += fabs( normed_velgrad(0)*derxy_(0,rr)
                    +normed_velgrad(1)*derxy_(1,rr)
                    +normed_velgrad(2)*derxy_(2,rr));
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
} // DRT::ELEMENTS::FluidGenalphaResVMM<distype>::CalcTau


/*----------------------------------------------------------------------*
 |                                                                      |
 | Get all global shape functions, first and eventually second          |
 | derivatives in a gausspoint                                          |
 |                                                         gammi 12/08  |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::FluidGenalphaResVMM<distype>::ShapeFunctionsFirstAndSecondDerivatives(
  const Fluid*                                ele             ,
  const int                                  & iquad           ,
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints       ,
  const std::vector<Epetra_SerialDenseVector>& myknots         ,
  const bool                                   higher_order_ele
  )
{
    // set gauss point coordinates
    LINALG::Matrix<nsd_,1> gp;

    // coordinates of the current integration point
    const double* gpcoord = (intpoints.IP().qxg)[iquad];
    for (int idim=0;idim<nsd_;idim++)
    {
      gp(idim) = gpcoord[idim];
    }

    //gp(0)=intpoints.IP().qxg[iquad][0];
    //gp(1)=intpoints.IP().qxg[iquad][1];
    //gp(2)=intpoints.IP().qxg[iquad][2];

    // TODO: 2D & 3D
    if(not (distype == DRT::Element::nurbs8
         ||
         distype == DRT::Element::nurbs27
         ||
         distype == DRT::Element::nurbs4
         ||
         distype == DRT::Element::nurbs9))
    {
      // get values of shape functions and derivatives in the gausspoint
      //DRT::UTILS::shape_function_3D(funct_,gp(0),gp(1),gp(2),distype);
      //DRT::UTILS::shape_function_3D_deriv1(deriv_,gp(0),gp(1),gp(2),distype);
      DRT::UTILS::shape_function<distype>(gp,funct_);
      DRT::UTILS::shape_function_deriv1<distype>(gp,deriv_);

      if (higher_order_ele)
      {
        // get values of shape functions and derivatives in the gausspoint
        //DRT::UTILS::shape_function_3D_deriv2(derxy2_,gp(0),gp(1),gp(2),distype);
        DRT::UTILS::shape_function_deriv2<distype>(gp, derxy2_);
      }
    }
    else
    {
      if (distype == DRT::Element::nurbs4 or distype == DRT::Element::nurbs9)
        dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());
      else if (higher_order_ele)
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
    // TODO
    // Take a closer look to the calculation of Jacobian, ... (performance)

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
    det_ = xjm_(0,0)*bm_(5,5)+
           xjm_(2,0)*bm_(3,5)+
           xjm_(0,2)*(bm_(5,3)-xjm_2_0_xjm_1_1)-
           xjm_(2,1)*bm_(3,4)-
           xjm_(0,1)*bm_(5,4);

    // check for degenerated elements
    if (det_ < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det_);
    }

    // set total integration factor
    double fac = intpoints.IP().qwgt[iquad]*det_;

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
    xji_(0,0) = ( bm_(5,5)-xjm_2_1_xjm_1_2)/det_;
    xji_(1,0) = (-bm_(5,4)+xjm_2_0_xjm_1_2)/det_;
    xji_(2,0) = ( bm_(5,3)-xjm_2_0_xjm_1_1)/det_;
    xji_(0,1) = (-bm_(4,5)+xjm_2_1_xjm_0_2)/det_;
    xji_(1,1) = ( bm_(4,4)-xjm_2_0_xjm_0_2)/det_;
    xji_(2,1) = (-bm_(4,3)+xjm_2_0_xjm_0_1)/det_;
    xji_(0,2) = ( bm_(3,5)-xjm_1_1_xjm_0_2)/det_;
    xji_(1,2) = (-bm_(3,4)+xjm_1_0_xjm_0_2)/det_;
    xji_(2,2) = ( bm_(3,3)-xjm_1_0_xjm_0_1)/det_;

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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::SetElementData(
  Fluid*                                      ele            ,
  const LINALG::Matrix<nsd_,iel>             & edispnp        ,
  const LINALG::Matrix<nsd_,iel>             & evelaf         ,
  const LINALG::Matrix<nsd_,iel>             & fsevelaf       ,
  const std::vector<Epetra_SerialDenseVector>& myknots        ,
  const double                               & timealphaF     ,
  double                                     & hk             ,
  double                                     & mk             ,
  Teuchos::RCP<const MAT::Material>            material       ,
  double                                     & visc           ,
  const enum INPAR::FLUID::FineSubgridVisc           fssgv          ,
  const enum INPAR::FLUID::TurbModelAction           turb_mod_action,
  const double                                 l_tau          ,
  double                                     & Cs             ,
  double                                     & Cs_delta_sq    ,
  double                                     & visceff        )
{

  // TODO
  // only a 3D function

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
  else if (ele->Shape()==Fluid::nurbs4 || ele->Shape()==Fluid::nurbs9)
    dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());


  // add displacement, when fluid nodes move in the ALE case
  if (ele->IsAle())
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
    //kinematic viscosity
    visc = actmat->Viscosity()/actmat->Density();
  }
  // initialise visceff to visc
  visceff=visc;

  // ---------------------------------------------------------------------------
  // Initialisation of tau computation: mk and hk

  // get element type constant mk for tau
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
  else if (distype == DRT::Element::nurbs4
     ||
     distype == DRT::Element::nurbs9)
   {
     dserror("%s is not a 3D Nurbs element",(DRT::DistypeToString(ele->Shape())).c_str());
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

  det_ = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
         xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
         xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
         xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
         xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
         xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

  // vol is used for some stabilisation parameters and the mixing length
  // definition of the Smagorinsky models
  vol_ = wquad*det_;

  // get element length hk for tau_M and tau_C: volume-equival. diameter/sqrt(3)
  hk = std::pow((6.*vol_/PI),(1.0/3.0))/sqrt(3.0);

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

  if(material->MaterialType()!= INPAR::MAT::m_fluid
     ||
     fssgv==INPAR::FLUID::smagorinsky_all
     ||
     fssgv==INPAR::FLUID::smagorinsky_small
     ||
     turb_mod_action !=INPAR::FLUID::no_model)
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
    xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det_;
    xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det_;
    xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det_;
    xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det_;
    xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det_;
    xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det_;
    xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det_;
    xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det_;
    xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det_;

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
    if(material->MaterialType()!=INPAR::MAT::m_fluid
       ||
       fssgv==INPAR::FLUID::smagorinsky_all
       ||
       turb_mod_action!=INPAR::FLUID::no_model)
    {
      rateofstrain  =GetStrainRate(evelaf,derxy_,vderxyaf_  );
    }

    // fine scale strains
    if(fssgv==INPAR::FLUID::smagorinsky_small)
    {
      fsrateofstrain=GetStrainRate(fsevelaf,derxy_,fsvderxyaf_);
    }


    // ---------------------------------------------------------------
    // compute nonlinear viscosity according to Carreau-Yasuda
    // ---------------------------------------------------------------
    if(material->MaterialType() != INPAR::MAT::m_fluid)
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
    if (turb_mod_action == INPAR::FLUID::smagorinsky_with_van_Driest_damping
        ||
        turb_mod_action == INPAR::FLUID::smagorinsky)
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

      if (turb_mod_action == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
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

      const double h_grid = std::pow((vol_),(1.0/3.0));

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
    else if(turb_mod_action == INPAR::FLUID::dynamic_smagorinsky)
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

    if (fssgv == INPAR::FLUID::smagorinsky_all)
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
    else if(fssgv == INPAR::FLUID::smagorinsky_small)
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
void DRT::ELEMENTS::FluidGenalphaResVMM<distype>::InterpolateToGausspoint(
  Fluid*                                     ele             ,
  const LINALG::Matrix<nsd_,iel>            & egridvelaf        ,
  const LINALG::Matrix<nsd_,iel>            & evelnp          ,
  const LINALG::Matrix<iel,1>               & eprenp          ,
  const LINALG::Matrix<nsd_,iel>            & eaccam          ,
  const LINALG::Matrix<nsd_,iel>            & evelaf          ,
  const LINALG::Matrix<nsd_,iel>            & fsevelaf        ,
  const double                              & visceff         ,
  const enum INPAR::FLUID::FineSubgridVisc          fssgv           ,
  const bool                                  higher_order_ele)
{
  // TODO
  // only a 3D function -> Interpolating may influence the calculation speed

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
    for(int mm=0;mm<nsd_;++mm)
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

    for(int rr=0;rr<3;++rr)
    {
      u_G_af_(rr)=funct_(0)*egridvelaf(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        u_G_af_(rr)+=funct_(nn)*egridvelaf(rr,nn);
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
  if (fssgv != INPAR::FLUID::no_fssgv) fsvderxyaf_.MultiplyNT(fsevelaf,derxy_);
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


