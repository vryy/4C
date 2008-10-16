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
  : vart_(0.0),
// nodal data
//-----------------------+------------+------------------------------------
//                  dim  | derivative | node
    weights_      (                     iel                             )
//-----------------------+------------+-----------------------------------
// gausspoint data
//------------------------------------------------------------------------
//                  dim  | derivative | node
//-----------------------+------------+------------------------------------
//----------------------+------------+------------------------------------
// element data
//------------------------------------------------------------------------
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


  // construct views
  LINALG::FixedSizeSerialDenseMatrix<4*iel,4*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<4*iel,    1> elevec1(elevec1_epetra.A(),true);

  // --------------------------------------------------
  // extract velocities, pressure and accelerations from the
  // global distributed vectors

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

  // create blitz matrix objects
  LINALG::FixedSizeSerialDenseMatrix<iel,1> eprenp    ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> evelnp    ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> evelaf    ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> eaccam    ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> edispnp   ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> egridvelaf;


  // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
  // Additionally only the 'velocity' components of my_velaf
  // and my_accam are important!
  for (int i=0;i<iel;++i)
  {
    int fi    =4*i;
    int fip   =fi+1;
    int fipp  =fip+1;


    eprenp(i)   = myvelnp[fipp+1];
    
    evelnp(0,i) = myvelnp[fi  ];
    evelnp(1,i) = myvelnp[fip ];
    evelnp(2,i) = myvelnp[fipp];

    evelaf(0,i) = myvelaf[fi  ];
    evelaf(1,i) = myvelaf[fip ];
    evelaf(2,i) = myvelaf[fipp];

    eaccam(0,i) = myaccam[fi  ];
    eaccam(1,i) = myaccam[fip ];
    eaccam(2,i) = myaccam[fipp];
  }

  if(ele->is_ale_)
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

  // --------------------------------------------------
  // set parameters for time integration
  ParameterList& timelist = params.sublist("time integration parameters");

  const double alphaM = timelist.get<double>("alpha_M");
  const double alphaF = timelist.get<double>("alpha_F");
  const double gamma  = timelist.get<double>("gamma");
  const double dt     = timelist.get<double>("dt");
  const double time   = timelist.get<double>("time");

  // if not available, define map from string to action
  if(ele->stabstrtoact_.empty())
  {
    ele->stabstrtoact_["quasistatic"            ]=Fluid3::subscales_quasistatic;
    ele->stabstrtoact_["time_dependent"         ]=Fluid3::subscales_time_dependent;
    ele->stabstrtoact_["no_transient"           ]=Fluid3::inertia_stab_drop;
    ele->stabstrtoact_["yes_transient"          ]=Fluid3::inertia_stab_keep;
    ele->stabstrtoact_["no_pspg"                ]=Fluid3::pstab_assume_inf_sup_stable;
    ele->stabstrtoact_["yes_pspg"               ]=Fluid3::pstab_use_pspg;
    ele->stabstrtoact_["no_supg"                ]=Fluid3::convective_stab_none;
    ele->stabstrtoact_["yes_supg"               ]=Fluid3::convective_stab_supg;
    ele->stabstrtoact_["no_vstab"               ]=Fluid3::viscous_stab_none;
    ele->stabstrtoact_["vstab_gls"              ]=Fluid3::viscous_stab_gls;
    ele->stabstrtoact_["vstab_gls_rhs"          ]=Fluid3::viscous_stab_gls_only_rhs;
    ele->stabstrtoact_["vstab_usfem"            ]=Fluid3::viscous_stab_usfem;
    ele->stabstrtoact_["vstab_usfem_rhs"        ]=Fluid3::viscous_stab_usfem_only_rhs;
    ele->stabstrtoact_["no_cstab"               ]=Fluid3::continuity_stab_none;
    ele->stabstrtoact_["cstab_qs"               ]=Fluid3::continuity_stab_yes;
    ele->stabstrtoact_["cstab_td"               ]=Fluid3::continuity_stab_td;
    ele->stabstrtoact_["no_cross"               ]=Fluid3::cross_stress_stab_none;
    ele->stabstrtoact_["cross_complete"         ]=Fluid3::cross_stress_stab;
    ele->stabstrtoact_["cross_rhs"              ]=Fluid3::cross_stress_stab_only_rhs;
    ele->stabstrtoact_["no_reynolds"            ]=Fluid3::reynolds_stress_stab_none;
    ele->stabstrtoact_["reynolds_rhs"           ]=Fluid3::reynolds_stress_stab_only_rhs;
    ele->stabstrtoact_["No"                     ]=Fluid3::fssgv_no;
    ele->stabstrtoact_["artificial_all"         ]=Fluid3::fssgv_artificial_all;
    ele->stabstrtoact_["artificial_small"       ]=Fluid3::fssgv_artificial_small;
    ele->stabstrtoact_["Smagorinsky_all"        ]=Fluid3::fssgv_Smagorinsky_all;
    ele->stabstrtoact_["Smagorinsky_small"      ]=Fluid3::fssgv_Smagorinsky_small;
    ele->stabstrtoact_["mixed_Smagorinsky_all"  ]=Fluid3::fssgv_mixed_Smagorinsky_all;
    ele->stabstrtoact_["mixed_Smagorinsky_small"]=Fluid3::fssgv_mixed_Smagorinsky_small;
    ele->stabstrtoact_["scale_similarity"       ]=Fluid3::fssgv_scale_similarity;
  }

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

  // get flag for fine-scale subgrid viscosity
  Fluid3::StabilisationAction fssgv =
    ele->ConvertStringToStabAction(params.get<string>("fs subgrid viscosity","No"));

  // get coarse- and fine-scale velocity as well as coarse-scale convective stress
  RCP<const Epetra_Vector> csvelaf;
  RCP<const Epetra_Vector> fsvelaf;
  RCP<const Epetra_Vector> csconvaf;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> fsevelaf ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> csevelaf ;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> cseconvaf;

  if (fssgv != Fluid3::fssgv_no)
  {
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
    if (fssgv == Fluid3::fssgv_mixed_Smagorinsky_all ||
        fssgv == Fluid3::fssgv_mixed_Smagorinsky_small ||
        fssgv == Fluid3::fssgv_scale_similarity)
    {
      csvelaf = discretization.GetState("csvelaf");
      if (csvelaf==null) dserror("Cannot get state vector 'csvelaf'");
      vector<double> mycsvelaf(lm.size());
      DRT::UTILS::ExtractMyValues(*csvelaf,mycsvelaf,lm);
      
      // get coarse-scale velocity and insert into element arrays
      for (int i=0;i<iel;++i)
      {
        csevelaf(0,i) = mycsvelaf[0+(i*4)];
        csevelaf(1,i) = mycsvelaf[1+(i*4)];
        csevelaf(2,i) = mycsvelaf[2+(i*4)];
      }
      
      csconvaf = discretization.GetState("csconvaf");
      if (csconvaf==null) dserror("Cannot get state vector 'csconvaf'");
      vector<double> mycsconvaf(lm.size());
      DRT::UTILS::ExtractMyValues(*csconvaf,mycsconvaf,lm);
      
      // get coarse-scale velocity and insert into element arrays
      for (int i=0;i<iel;++i)
      {
        cseconvaf(0,i) = mycsconvaf[0+(i*4)];
        cseconvaf(1,i) = mycsconvaf[1+(i*4)];
        cseconvaf(2,i) = mycsconvaf[2+(i*4)];
      }
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

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  double Cs            = 0.0;
  double Cs_delta_sq   = 0.0;
  double l_tau         = 0.0;
  double visceff       = 0.0;

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv != Fluid3::fssgv_no && turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
    dserror("No combination of a classical (all-scale) turbulence model and a fine-scale subgrid-viscosity approach currently possible!");
  if (fssgv != Fluid3::fssgv_no) Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);
  
  // the default action is no model
  Fluid3::TurbModelAction turb_mod_action = Fluid3::no_model;

  int  nlayer=0;
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
      RefCountPtr<vector<double> > planecoords      = turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");

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

  // --------------------------------------------------
  // specify what to compute
  const bool compute_elemat = params.get<bool>("compute element matrix");

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

  // --------------------------------------------------
  // calculate element coefficient matrix
  Sysmat(
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
    csevelaf,
    fsevelaf,
    cseconvaf,
    actmat,
    alphaM,
    alphaF,
    gamma,
    dt,
    time,
    newton,
    higher_order_ele,
    fssgv,
    tds,
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

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

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

  return 0;
}

/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |                            (public)                      gammi 06/07|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::Sysmat(
  Fluid3*                                          ele             ,
  std::vector<blitz::Array<double,1> >&            myknots         ,
  LINALG::FixedSizeSerialDenseMatrix<4*iel,4*iel>& elemat          ,
  LINALG::FixedSizeSerialDenseMatrix<4*iel,1>&     elevec          ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& edispnp         ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& egridvaf        ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evelnp          ,
  const LINALG::FixedSizeSerialDenseMatrix<iel,1>& eprenp          ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& eaccam          ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evelaf          ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& csevelaf        ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& fsevelaf        ,
  const LINALG::FixedSizeSerialDenseMatrix<3,iel>& cseconvaf       ,
  struct _MATERIAL*                                material        ,
  double                                           alphaM          ,
  double                                           alphaF          ,
  double                                           gamma           ,
  double                                           dt              ,
  double                                           time            ,
  enum Fluid3::LinearisationAction                 newton          ,
  bool                                             higher_order_ele,
  enum Fluid3::StabilisationAction                 fssgv           ,
  enum Fluid3::StabilisationAction                 tds             ,
  enum Fluid3::StabilisationAction                 inertia         ,
  enum Fluid3::StabilisationAction                 pspg            ,
  enum Fluid3::StabilisationAction                 supg            ,
  enum Fluid3::StabilisationAction                 vstab           ,
  enum Fluid3::StabilisationAction                 cstab           ,
  enum Fluid3::StabilisationAction                 cross           ,
  enum Fluid3::StabilisationAction                 reynolds        ,
  enum Fluid3::TauType                             whichtau        ,
  enum Fluid3::TurbModelAction                     turb_mod_action ,
  double&                                          Cs              ,
  double&                                          Cs_delta_sq     ,
  double&                                          visceff         ,
  double                                           l_tau           ,
  bool                                             compute_elemat
  )
{

  //------------------------------------------------------------------
  //   We use LAPACK libraries --- get the epetra wrapper
  //------------------------------------------------------------------
  Epetra_LAPACK      solver;

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
  if(material->mattyp == m_fluid)
  {
    visc = material->m.fluid->viscosity;
  }

  //------------------------------------------------------------------
  //                      SET ELEMENT DATA
  //------------------------------------------------------------------
  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
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
        
  LINALG::FixedSizeSerialDenseMatrix<3,1> gp;
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

  // get element type constant for tau
  double mk=0.0;
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
  for(int rr=0;rr<3;++rr)
  {
    velintnp_(rr)=funct_(0)*evelnp(rr,0);
    for(int nn=1;nn<iel;++nn)
    {
      velintnp_(rr)+=funct_(nn)*evelnp(rr,nn);
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

  // get velocity norms
  const double vel_normaf = velintaf_.Norm2();
  const double vel_normnp = velintnp_.Norm2();

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
  {
    CalVisc( material, visc);
  }

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
      LINALG::FixedSizeSerialDenseMatrix<3,3> epsilon;

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
      LINALG::FixedSizeSerialDenseMatrix<3,1> centernodecoord;

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

    visceff = visc + Cs_delta_sq * rateofstrain;
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
      LINALG::FixedSizeSerialDenseMatrix<3,3> epsilon;

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

    visceff = visc + Cs_delta_sq * rateofstrain;

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
      LINALG::FixedSizeSerialDenseMatrix<3,3> G;

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
      LINALG::FixedSizeSerialDenseMatrix<3,1> g;

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
      LINALG::FixedSizeSerialDenseMatrix<3,3> G;
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
      LINALG::FixedSizeSerialDenseMatrix<3,1> g;

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
      LINALG::FixedSizeSerialDenseMatrix<3,1> normed_velintaf;

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
    else if(whichtau == Fluid3::codina)
    {
      // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
      // tau_M: Barrenechea, G.R. and Valentin, F.
      // tau_C: Codina


      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::FixedSizeSerialDenseMatrix<3,1> normed_velintaf;

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
      fsvelintaf_.Multiply(fsevelaf,funct_);

      // get fine-scale velocity norm
      fsvel_normaf = fsvelintaf_.Norm2();
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
      {
        fsvderxyaf_.MultiplyNT(fsevelaf,derxy_);
      }
      else
      {
        for(int rr=0;rr<3;rr++)
        {
          for(int mm=0;mm<3;mm++)
          {
            fsvderxyaf_(rr,mm) = vderxyaf_(rr,mm);
          }
        }
      }

      LINALG::FixedSizeSerialDenseMatrix<3,3> epsilon;

      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<3;mm++)
        {
          epsilon(rr,mm)=0.5 * ( fsvderxyaf_(rr,mm) + fsvderxyaf_(mm,rr) );
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
    // Choices of the fine-scale Smagorinsky constant Cs_fs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)

    vart_ = Cs * Cs * hk * hk * rateofstrain;
  }

  //----------------------------------------------------------------------------
  //
  //    From here onwards, we are working on the gausspoints of the element
  //            integration, not on the element center anymore!
  //
  //----------------------------------------------------------------------------

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  // remember whether the subscale quantities have been allocated an set to zero.
  if(tds == Fluid3::subscales_time_dependent)
  {
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

  // just define certain constants for conveniance
  const double afgdt  = alphaF * gamma * dt;

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    // set gauss point coordinates
    LINALG::FixedSizeSerialDenseMatrix<3,1> gp;

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
      solver.GETRF(6,6,bm_.A(),6,&(pivot[0]),&ierr);

      if (ierr!=0)
      {
        dserror("Unable to perform LU factorisation during computation of derxy2");
      }

      // backward substitution. GETRS replaces the input (chainrulerhs, currently
      // stored on derxy2) with the result
      solver.GETRS('N',6,iel,bm_.A(),6,&(pivot[0]),derxy2_.A(),6,&ierr);

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
    for(int rr=0;rr<3;++rr)
    {
      velintnp_(rr)=funct_(0)*evelnp(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        velintnp_(rr)+=funct_(nn)*evelnp(rr,nn);
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
    for(int rr=0;rr<3;++rr)
    {
      bodyforceaf_(rr)=funct_(0)*edeadaf_(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        bodyforceaf_(rr)+=funct_(nn)*edeadaf_(rr,nn);
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
    const double divunp = (vderxynp_(0,0)+vderxynp_(1,1)+vderxynp_(2,2));

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
    for(int rr=0;rr<3;++rr)
    {
      aleconvintaf_(rr)=velintaf_(rr);
    }    

    if (ele->is_ale_)
    {
      for(int rr=0;rr<3;++rr)
      {
        for(int nn=0;nn<iel;++nn)
        {
          aleconvintaf_(rr)-=funct_(nn)*egridvaf(rr,nn);
        }
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
    if (fssgv != Fluid3::fssgv_no  && fssgv != Fluid3::fssgv_scale_similarity)
    {
      fsvderxyaf_.MultiplyNT(fsevelaf,derxy_);
    }
    else 
    {
      fsvderxyaf_.Clear();
    }

    // get values at integration point required for scale-similarity model
    if(fssgv == Fluid3::fssgv_scale_similarity ||
       fssgv == Fluid3::fssgv_mixed_Smagorinsky_all ||
       fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
    {
      // get coarse-scale velocities at integration point
      csvelintaf_.Multiply(csevelaf,funct_);

      // get coarse-scale velocity (np,i) derivatives at integration point
      csvderxyaf_.MultiplyNT(csevelaf,derxy_);

      // PR(u) * grad PR(u): */
      convaf_s_.MultiplyTN(csvderxyaf_,csvelintaf_);

      // get coarse-scale convective stresses at integration point
      csconvintaf_.Multiply(cseconvaf,funct_);
    }

    if (higher_order_ele)
    {
      /*--- viscous term  2* grad * epsilon(u): --------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
           |                                                |
           |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
           |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                             */


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

	viscaf_old_(0) += 
	      (viscs2_(0,mm)*evelaf(0,mm)
	       +
	       derxy2_(3,mm)*evelaf(1,mm)
	       +
	       derxy2_(4,mm)*evelaf(2,mm));
	viscaf_old_(1) += 
	      (derxy2_(3,mm)*evelaf(0,mm)
	       +
	       viscs2_(1,mm)*evelaf(1,mm)
	       +
	       derxy2_(5,mm)*evelaf(2,mm));
	viscaf_old_(2) +=
	      (derxy2_(4,mm)*evelaf(0,mm)
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

      // update estimates for the subscale quantities

      const double factauC                  = tauC/(tauC+dt);
      const double facMtau                  = 1./(alphaM*tauM+afgdt);

      if (cstab == Fluid3::continuity_stab_td)
      {
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
      }

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

      /*
        This is the operator

                  /~n+af         \
                 | u      o nabla |
                  \   (i)        /

        required for the cross stress linearisation

      */
      if(cross == Fluid3::cross_stress_stab)
      {
	for (int rr=0;rr<iel;++rr)
	{
	  conv_subaf_(rr) = svelaf_(0)*derxy_(0,rr);
	  
	  for (int mm=1;mm<3;++mm)
	  {
	    conv_subaf_(rr) += svelaf_(mm)*derxy_(mm,rr);
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
	 cross == Fluid3::cross_stress_stab)
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
	double cstabfac;
	if(cstab == Fluid3::continuity_stab_yes)
	{ // cstab Codina style
	  cstabfac=fac*gamma*dt*tauC;
	} // cstab Codina style
	else if(cstab == Fluid3::continuity_stab_td)
	{ // time dependent version
	  cstabfac=fac*gamma*dt*dt*factauC;
	} // time dependent version
	else
	{ // no cstab
	  cstabfac=0;
	} // no cstab
	
        const double fac_gamma_dt      = fac*gamma*dt;
	const double fac_afgdt_visceff = fac*visceff*afgdt;

        //---------------------------------------------------------------
        //
	//          SUBSCALE ACCELERATION STABILISATION
	//        RESCALING FACTORS FOR GALERKIN 1 TERMS AND
	//              COMPUTATION OF EXTRA TERMS
        //
        //---------------------------------------------------------------
	  
	if(inertia == Fluid3::inertia_stab_keep)
	{
	  // rescale time factors terms affected by inertia stabilisation
	  fac_inertia   *=afgdt*facMtau;
	  fac_convection*=afgdt*facMtau;
	  
	  // do inertia stabilisation terms which are not scaled
	  // Galerkin terms since they are not partially integrated

	  const double fac_alphaM_tauM_facMtau          = fac*alphaM*tauM*facMtau;
	  
	  for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	  {
	    const int fvi    =4*vi;
	    const int fvip   =fvi+1;
	    const int fvipp  =fvi+2;
	    
	    const double fac_alphaM_gamma_dt_tauM_facMtau_funct_vi=fac_alphaM_tauM_facMtau*gamma*dt*funct_(vi);
	    
	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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
	    } // end loop rows (test functions for matrix)
	  } // end loop rows (solution for matrix, test function for vector)

	  if(higher_order_ele && newton!=Fluid3::minimal)
	  {
	    const double fac_visceff_afgdt_alphaM_tauM_facMtau
	      =
	      fac*visceff*afgdt*alphaM*tauM*facMtau;
	    
	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      const int fvi    =4*vi;
	      const int fvip   =fvi+1;
	      const int fvipp  =fvi+2;
		
	      const double fac_visceff_afgdt_alphaM_tauM_facMtau_funct_vi
		=
		fac_visceff_afgdt_alphaM_tauM_facMtau*funct_(vi);
		  
	      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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
	      } // end loop rows (test functions for matrix)
	    } // end loop rows (solution for matrix, test function for vector)
	  } // end higher order element and  linearisation of linear terms not supressed
	} // extra terms for inertia stab
	
	//---------------------------------------------------------------
        //
        //              TIME-DEPENDENT SUBGRID-SCALES
	//
        //       GALERKIN PART 1 (INERTIA, CONVECTION, VISCOUS)
	//                     SUPG STABILISATION
	//                  CONTINUITY STABILISATION
        //
        //---------------------------------------------------------------

	if(supg == Fluid3::convective_stab_supg)
	{
	  const double fac_afgdt_tauM_afgdt_facMtau  = fac*afgdt*afgdt*facMtau*tauM;

          const double fac_afgdt_tauM_facMtau        = fac*afgdt*tauM*facMtau;

	  const double fac_alphaM_afgdt_tauM_facMtau
	    =
	    fac*alphaM*afgdt*tauM*facMtau;

	  for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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
	      fac_alphaM_afgdt_tauM_facMtau*funct_(ui)+fac_afgdt_tauM_afgdt_facMtau*conv_c_af_(ui);
	    
	    /* CSTAB entries */
	    const double cstab_0 = cstabfac*derxy_(0,ui);
	    const double cstab_1 = cstabfac*derxy_(1,ui);
	    const double cstab_2 = cstabfac*derxy_(2,ui);

	    /* combined CSTAB/viscous entires */
	    const double visc_and_cstab_0 = visc_0+cstab_0;
	    const double visc_and_cstab_1 = visc_1+cstab_1;
	    const double visc_and_cstab_2 = visc_2+cstab_2;
	      
	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      const int fvi    =4*vi;
	      const int fvip   =fvi+1;
	      const int fvipp  =fvi+2;
	      
	      const double sum =
		inertia_and_conv_ui*funct_(vi)
		+
		supg_inertia_and_conv_ui*conv_c_af_(vi)
		+
		visc_0*derxy_(0,vi)
		+
		visc_1*derxy_(1,vi)
		+
		visc_2*derxy_(2,vi);
		  
	      /*
		inertia term (intermediate)
		
		factor with inertia stabilisation:

                                                       
                                alphaF*gamma*dt        
                    alphaM*--------------------------- 
                           alphaM*tauM+alphaF*gamma*dt 
                                                       

	        factor without inertia stabilisation:
		  
		    alphaM


                              /          \
                             |            |
                             |  Dacc , v  | +
                             |            |
                              \          /

	      */
	      /*

		  +convection (intermediate)

		  factor without inertia stabilisation: 

		  +alphaF*gamma*dt


		  factor with inertia stabilisation: 

                                           alphaF*gamma*dt
                    +alphaF*gamma*dt*---------------------------
                                     alphaM*tauM+alphaF*gamma*dt


                          /                          \
                         |  / n+af       \            |
                         | | c    o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
	      */

	      /* viscous term (intermediate) */
	      
	      /*  factor: +2*nu*alphaF*gamma*dt
		   
                          /                          \
		         |       /    \         / \   |
		         |  eps | Dacc | , eps | v |  |
	        	 |       \    /         \ /   |
                          \                          /
	      */
		  
	      /* CONTINUITY stabilisation */

	      /*
		  Codina style: (recommended)
		  
                  factor: +gamma*dt*tauC

		  Time dependent: (buggy)
 
		                       tauC * dt
		  factor: +gamma* dt * ---------
                                       tauC + dt  

		  No cstab:
				      
                  factor: 0
		  
                     /                          \
                    |                            |
                    | nabla o Dacc  , nabla o v  |
                    |                            |
                     \                          /
	      */

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

	  if (newton==Fluid3::Newton) // if newton and supg
	  {
	    const double fac_afgdt_tauM_afgdt_facMtau = fac*afgdt*afgdt*facMtau*tauM;

	    double temp[3][3];

	    const double fac_afgdt_svelaf_0 = fac*afgdt*svelaf_(0);
	    const double fac_afgdt_svelaf_1 = fac*afgdt*svelaf_(1);
	    const double fac_afgdt_svelaf_2 = fac*afgdt*svelaf_(2);

	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      const int fvi  =4*vi;
	      const int fvip =fvi+1;
	      const int fvipp=fvi+2;

	      const double scaled_inertia_and_conv_vi 
		= 
		fac_convection*funct_(vi)
		+
		fac_afgdt_tauM_afgdt_facMtau*conv_c_af_(vi);
	      
	      temp[0][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,0)-fac_afgdt_svelaf_0*derxy_(0,vi);
	      temp[1][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,1)-fac_afgdt_svelaf_0*derxy_(1,vi);
	      temp[2][0]=scaled_inertia_and_conv_vi*vderxyaf_(0,2)-fac_afgdt_svelaf_0*derxy_(2,vi);
	      temp[0][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,0)-fac_afgdt_svelaf_1*derxy_(0,vi);
	      temp[1][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,1)-fac_afgdt_svelaf_1*derxy_(1,vi);
	      temp[2][1]=scaled_inertia_and_conv_vi*vderxyaf_(1,2)-fac_afgdt_svelaf_1*derxy_(2,vi);
	      temp[0][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,0)-fac_afgdt_svelaf_2*derxy_(0,vi);
	      temp[1][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,1)-fac_afgdt_svelaf_2*derxy_(1,vi);
	      temp[2][2]=scaled_inertia_and_conv_vi*vderxyaf_(2,2)-fac_afgdt_svelaf_2*derxy_(2,vi);
	      
	      for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	      {
		const int fui  =4*ui;
		const int fuip =fui+1;
		const int fuipp=fui+2;

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

		/* SUPG stabilisation --- subscale velocity, nonlinear part from testfunction

		   factor:
                          alphaF * gamma * dt


                    /                            \
                   |  ~n+af    /            \     |
                   |  u     , | Dacc o nabla | v  |
                   |   (i)     \            /     |
                    \                            /

		*/

		/* convection (intermediate)

		   factor without inertia stabilisation: 

   		      +alphaF*gamma*dt


		   factor with inertia stabilisation: 

                                             alphaF*gamma*dt
                      +alphaF*gamma*dt*---------------------------
                                       alphaM*tauM+alphaF*gamma*dt

                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
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
	      }
	    }
	  } // end if supg and newton

	  for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	  {
	    const int fuippp=4*ui+3;

	    const double scaled_gradp_0 = fac_afgdt_tauM_facMtau*gamma*dt*derxy_(0,ui);
	    const double scaled_gradp_1 = fac_afgdt_tauM_facMtau*gamma*dt*derxy_(1,ui);
	    const double scaled_gradp_2 = fac_afgdt_tauM_facMtau*gamma*dt*derxy_(2,ui);

	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      const int fvi=4*vi;
 
	      /* SUPG stabilisation --- pressure

               factor:
                               alphaF*gamma*tauM*dt
                            ---------------------------, rescaled by gamma*dt
                            alphaM*tauM+alphaF*gamma*dt


                    /                               \
                   |              / n+af       \     |
                   |  nabla Dp , | c    o nabla | v  |
                   |              \            /     |
                    \                               /
              */
              elemat(fvi  ,fuippp) += scaled_gradp_0*conv_c_af_(vi);
              elemat(fvi+1,fuippp) += scaled_gradp_1*conv_c_af_(vi);
              elemat(fvi+2,fuippp) += scaled_gradp_2*conv_c_af_(vi);
	      
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)

	} // end supg
	else
	{ // no supg

	  for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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

	    /* CSTAB entries */
	    const double cstab_0 = cstabfac*derxy_(0,ui);
	    const double cstab_1 = cstabfac*derxy_(1,ui);
	    const double cstab_2 = cstabfac*derxy_(2,ui);

	    /* combined CSTAB/viscous entires */
	    const double visc_and_cstab_0 = visc_0+cstab_0;
	    const double visc_and_cstab_1 = visc_1+cstab_1;
	    const double visc_and_cstab_2 = visc_2+cstab_2;
	      
	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      const int fvi    =4*vi;
	      const int fvip   =fvi+1;
	      const int fvipp  =fvi+2;
	      
	      const double sum =
		inertia_and_conv_ui*funct_(vi)
		+
		visc_0*derxy_(0,vi)
		+
		visc_1*derxy_(1,vi)
		+
		visc_2*derxy_(2,vi);
		  
	      /*
		inertia term (intermediate)


                              /          \
                             |            |
                    alphaM * |  Dacc , v  | +
                             |            |
                              \          /

	      */
	      /*

		  +convection (intermediate)


		  factor without inertia stabilisation: 

   		      +alphaF*gamma*dt


		  factor with inertia stabilisation: 

                                             alphaF*gamma*dt
                      +alphaF*gamma*dt*---------------------------
                                       alphaM*tauM+alphaF*gamma*dt


                          /                          \
                         |  / n+af       \            |
                         | | c    o nabla | Dacc , v  |
                         |  \            /            |
                          \                          /
	      */

	      /* viscous term (intermediate) */
	      
	      /*  factor: +2*nu*alphaF*gamma*dt
		   
                          /                          \
		         |       /    \         / \   |
		         |  eps | Dacc | , eps | v |  |
	        	 |       \    /         \ /   |
                          \                          /
		  */
		  
	      /* CONTINUITY stabilisation */

	      /*
		Codina style: (recommended)
		  
                  factor: +gamma*dt*tauC

		Time dependent: (buggy)
 
		                       tauC * dt
		  factor: +gamma* dt * ---------
                                       tauC + dt  

		No cstab:
				      
                  factor: 0
		  
                     /                          \
                    |                            |
                    | nabla o Dacc  , nabla o v  |
                    |                            |
                     \                          /
	      */
	      
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

	  if (newton==Fluid3::Newton) // if no supg and newton
	  {
	      
	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	    {
	      const int fui    =4*ui;
	      const int fuip   =fui+1;
	      const int fuipp  =fui+2;
	      
	      const double fac_convection_funct_ui = fac_convection*funct_(ui);
	      
	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	      {
		const int fvi    =4*vi;
		const int fvip   =fvi+1;
		const int fvipp  =fvi+2;
	      
		const double fac_convection_funct_ui_funct_vi = fac_convection_funct_ui*funct_(vi);
		
		/*

		  factor without inertia stabilisation: 
		  
   		      +alphaF*gamma*dt


		  factor with inertia stabilisation: 

                                             alphaF*gamma*dt
                      +alphaF*gamma*dt*---------------------------
                                       alphaM*tauM+alphaF*gamma*dt

                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
		*/
		elemat(fvi  ,fui  ) += fac_convection_funct_ui_funct_vi*vderxyaf_(0,0);
		elemat(fvi  ,fuip ) += fac_convection_funct_ui_funct_vi*vderxyaf_(0,1);
		elemat(fvi  ,fuipp) += fac_convection_funct_ui_funct_vi*vderxyaf_(0,2);
		elemat(fvip ,fui  ) += fac_convection_funct_ui_funct_vi*vderxyaf_(1,0);
		elemat(fvip ,fuip ) += fac_convection_funct_ui_funct_vi*vderxyaf_(1,1);
		elemat(fvip ,fuipp) += fac_convection_funct_ui_funct_vi*vderxyaf_(1,2);
		elemat(fvipp,fui  ) += fac_convection_funct_ui_funct_vi*vderxyaf_(2,0);
		elemat(fvipp,fuip ) += fac_convection_funct_ui_funct_vi*vderxyaf_(2,1);
		elemat(fvipp,fuipp) += fac_convection_funct_ui_funct_vi*vderxyaf_(2,2);
		
	      } // end loop rows (test functions for matrix)
	    } // end loop columns (solution for matrix, test function for vector)
	  } // end if no inertia and newton
	} // no supg
       
        //---------------------------------------------------------------
        //
        //               TIME-DEPENDENT SUBGRID-SCALES
	//
        //                     GALERKIN PART 2 
	//     (REMAINING PRESSURE AND CONTINUITY EXPRESSIONS)
        //
        //---------------------------------------------------------------

	for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	{
	  const int fuippp  =4*ui+3;

	  const double fac_gamma_dt_funct_ui=fac_gamma_dt*funct_(ui);
	  
	  for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	  {
	    const int fvi  =4*vi;
	    const int fvip =fvi+1;
	    const int fvipp=fvi+2;


	    /* pressure (implicit) */

	    /*  factor: -1, rescaled by gamma*dt

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
	    */

            /* continuity equation (implicit) */
	    
            /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
            */

	    elemat(fvi   ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(0,vi);
	    elemat(fvip  ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(1,vi);
	    elemat(fvipp ,fuippp) -= fac_gamma_dt_funct_ui*derxy_(2,vi);

            elemat(fuippp,fvi   ) += fac_gamma_dt_funct_ui*derxy_(0,vi);
            elemat(fuippp,fvip  ) += fac_gamma_dt_funct_ui*derxy_(1,vi);
            elemat(fuippp,fvipp ) += fac_gamma_dt_funct_ui*derxy_(2,vi);

	  } // end loop rows (test functions for matrix)
	} // end loop columns (solution for matrix, test function for vector)
        // end remaining Galerkin terms

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

            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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
                 
              for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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
            for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
            {
	      const int fui  =4*ui;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;

	      const double inertia_and_conv_ui
		=
		fac_alphaM_gamma_dt_tauM_facMtau*funct_(ui)
		+
		fac_afgdt_gamma_dt_tauM_facMtau*conv_c_af_(ui);

              for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
	    const int fuippp=4*ui+3;
	    const double scaled_derxy_0=fac_gdt_gdt_tauM_facMtau*derxy_(0,ui);
	    const double scaled_derxy_1=fac_gdt_gdt_tauM_facMtau*derxy_(1,ui);
	    const double scaled_derxy_2=fac_gdt_gdt_tauM_facMtau*derxy_(2,ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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

            for (int vi=0; vi<iel; ++vi) // loop columns (solution for matrix, test function for vector)
            {
	      const int fvippp=4*vi+3;

	      const double a=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,0)+derxy_(1,vi)*vderxyaf_(1,0)+derxy_(2,vi)*vderxyaf_(2,0));
	      const double b=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,1)+derxy_(1,vi)*vderxyaf_(1,1)+derxy_(2,vi)*vderxyaf_(2,1));
	      const double c=fac_afgdt_gamma_dt_tauM_facMtau*(derxy_(0,vi)*vderxyaf_(0,2)+derxy_(1,vi)*vderxyaf_(1,2)+derxy_(2,vi)*vderxyaf_(2,2));

              for (int ui=0; ui<iel; ++ui)  // loop rows (test functions for matrix)
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
        } // end cross
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
      if(inertia == Fluid3::inertia_stab_keep)
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

        //---------------------------------------------------------------
        //
        //   GALERKIN PART 2 (REMAINING EXPRESSIONS)
        //
        //---------------------------------------------------------------
        {
          const double fac_prenp_  = fac*prenp_ ;

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

            elevec(fui  ) += fac_prenp_*derxy_(0,ui) ;
            elevec(fui+1) += fac_prenp_*derxy_(1,ui) ;
            elevec(fui+2) += fac_prenp_*derxy_(2,ui) ;
          }

          {
            const double fac_visceff = fac*visceff;

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

              elevec(fui  ) -= fac_visceff*
                               (derxy_(0,ui)*vderxyaf_(0,0)*2.0
				+
				derxy_(1,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
				+
				derxy_(2,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0)));
              elevec(fui+1) -= fac_visceff*
      		               (derxy_(0,ui)*(vderxyaf_(0,1)+vderxyaf_(1,0))
				+
				derxy_(1,ui)*vderxyaf_(1,1)*2.0
				+
				derxy_(2,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1)));
              elevec(fui+2) -= fac_visceff*
		               (derxy_(0,ui)*(vderxyaf_(0,2)+vderxyaf_(2,0))
				+
				derxy_(1,ui)*(vderxyaf_(1,2)+vderxyaf_(2,1))
				+
				derxy_(2,ui)*vderxyaf_(2,2)*2.0);
            }
          }

          {
            const double fac_divunp  = fac*divunp;

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

        {
          const double fac_prenp=fac*prenp_;

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
	    const int fui=4*ui;
            /* pressure */

            /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
            */

            elevec(fui  ) += fac_prenp*derxy_(0,ui);
            elevec(fui+1) += fac_prenp*derxy_(1,ui);
            elevec(fui+2) += fac_prenp*derxy_(2,ui);
          }
        }

        {
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
        }

        {
          const double fac_divunp=fac*divunp;

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
        }
      }

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

      if(cstab == Fluid3::continuity_stab_yes)
      {
        //---------------------------------------------------------------
        //
        //            QUASISTATIC VERSION OF RESIDUAL BASED
        //                  CONTINUITY STABILISATION
        //          (the original version proposed by Codina)
        //
        //---------------------------------------------------------------

        const double fac_tauC_divunp = fac*tauC*divunp;

        for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
        {
	  const int fui=4*ui;

          /* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
          */

          elevec(fui  ) -= fac_tauC_divunp*derxy_(0,ui) ;
          elevec(fui+1) -= fac_tauC_divunp*derxy_(1,ui) ;
          elevec(fui+2) -= fac_tauC_divunp*derxy_(2,ui) ;
        } // end loop rows
      } // cstab_qs
      else if (cstab == Fluid3::continuity_stab_td)
      {
        //---------------------------------------------------------------
        //
        //           TIME DEPENDENT VERSION OF RESIDUAL BASED
        //                  CONTINUITY STABILISATION
        //
        //---------------------------------------------------------------

        const double fac_sprenp = fac*sprenp(iquad);

        for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
        {
	  const int fui=4*ui;

          /* factor: -1

                       /                  \
                      |  ~n+1              |
                      |  p    , nabla o v  |
                      |   (i)              |
                       \                  /
          */
          elevec(fui  ) += fac_sprenp*derxy_(0,ui);
          elevec(fui+1) += fac_sprenp*derxy_(1,ui);
          elevec(fui+2) += fac_sprenp*derxy_(2,ui);
        } // ui
      } // cstab_td


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
      if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
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
    } // end time dependent subscales
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
        //   GALERKIN PART, INERTIA, CONVECTION AND VISCOUS TERMS
	//                QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------

	//---------------------------------------------------------------
	//
	//         LEAST SQUARES CONTINUITY STABILISATION PART,
	//            QUASISTATIC FORMULATION (IF ACTIVE)
	//
	//---------------------------------------------------------------

        //---------------------------------------------------------------
        //
        //           SUPG PART, INERTIA AND CONVECTION TERMS
        //             QUASISTATIC FORMULATION (IF ACTIVE)
        //
        //---------------------------------------------------------------

	const double fac_afgdt         = fac*afgdt;
	const double fac_visceff_afgdt = fac_afgdt*visceff;
	const double fac_gamma_dt      = fac*gamma*dt;
	const double fac_alphaM        = fac*alphaM;

	if(supg == Fluid3::convective_stab_supg)
	{
	  if(cstab == Fluid3::continuity_stab_yes)
	  { // end supg and cstab
	    const double fac_gamma_dt_tauC = fac*gamma*dt*tauC;

	    for (int ui=0; ui<iel; ++ui) // loop columns 
	    {
	      const int fui  =4*ui;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
	      
	      /* GALERKIN inertia term (intermediate) + convection (intermediate) */
	      
	      const double inertia_and_conv_ui
		= fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);
	      
	      /* SUPG stabilisation --- inertia and convection */
	      const double fac_tauM_inertia_and_conv
		= tauM*inertia_and_conv_ui;
	      
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
		
		const double sum =
		  inertia_and_conv_ui*funct_(vi)
		  +
		  fac_tauM_inertia_and_conv*conv_c_af_(vi)
		  +
		  fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
		  +
		  fac_visceff_afgdt_derxy1_ui*derxy_(1,vi)
		  +
		  fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
	  
		/*
   		  GALERKIN inertia term (intermediate) + convection (intermediate)

		  factor:

                          /          \                   /                          \
                         |            |                 |  / n+af       \            |
                +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                         |            |                 |  \            /            |
                          \          /                   \                          /
		*/

		/*  GALERKIN viscous term (intermediate) */
	    
		/*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
		*/
		
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

		/* CONTINUITY stabilisation */

		/* factor: +gamma*dt*tauC
		  
                    /                          \
                   |                            |
                   | nabla o Dacc  , nabla o v  |
                   |                            |
                    \                          /
		*/

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
	  } // end supg and cstab
	  else
	  { // supg without cstab
	    for (int ui=0; ui<iel; ++ui) // loop columns 
	    {
	      const int fui  =4*ui;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
	      
	      /* GALERKIN inertia term (intermediate) + convection (intermediate) */
	      
	      const double inertia_and_conv_ui
		= fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);
	      
	      /* SUPG stabilisation --- inertia and convection */
	      const double fac_tauM_inertia_and_conv
		= tauM*inertia_and_conv_ui;
	      
	      /* viscous term (intermediate), diagonal parts */
	      const double fac_visceff_afgdt_derxy0_ui=fac_visceff_afgdt*derxy_(0,ui);
	      const double fac_visceff_afgdt_derxy1_ui=fac_visceff_afgdt*derxy_(1,ui);
	      const double fac_visceff_afgdt_derxy2_ui=fac_visceff_afgdt*derxy_(2,ui);
	      
	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	      {
		const int fvi    =4*vi;
		const int fvip   =fvi+1;
		const int fvipp  =fvi+2;
		
		const double sum =
		  inertia_and_conv_ui*funct_(vi)
		  +
		  fac_tauM_inertia_and_conv*conv_c_af_(vi)
		  +
		  fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
		  +
		  fac_visceff_afgdt_derxy1_ui*derxy_(1,vi)
		  +
		  fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
		
		/*
		  GALERKIN inertia term (intermediate) + convection (intermediate)
		  
		  factor:
		  
                          /          \                   /                          \
                         |            |                 |  / n+af       \            |
                +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                         |            |                 |  \            /            |
                          \          /                   \                          /
		*/

		/*  GALERKIN viscous term (intermediate) */
	    
		/*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
		*/
		
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

		elemat(fvi  ,fui  ) += sum+fac_visceff_afgdt_derxy0_ui*derxy_(0,vi);
		elemat(fvi  ,fuip ) +=     fac_visceff_afgdt_derxy0_ui*derxy_(1,vi);
		elemat(fvi  ,fuipp) +=     fac_visceff_afgdt_derxy0_ui*derxy_(2,vi);
		elemat(fvip ,fui  ) +=     fac_visceff_afgdt_derxy1_ui*derxy_(0,vi);
		elemat(fvip ,fuip ) += sum+fac_visceff_afgdt_derxy1_ui*derxy_(1,vi);
		elemat(fvip ,fuipp) +=     fac_visceff_afgdt_derxy1_ui*derxy_(2,vi);
		elemat(fvipp,fui  ) +=     fac_visceff_afgdt_derxy2_ui*derxy_(0,vi);
		elemat(fvipp,fuip ) +=     fac_visceff_afgdt_derxy2_ui*derxy_(1,vi);
		elemat(fvipp,fuipp) += sum+fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
	      } // vi
	    } // ui
	  } // end supg without cstab


	  //---------------------------------------------------------------
	  //
	  //                  GALERKIN AND SUPG PART
	  //    REACTIVE TYPE LINEARISATIONS, QUASISTATIC FORMULATION
	  //
	  //---------------------------------------------------------------

          if (newton==Fluid3::Newton)
          {
	    double temp[3][3];
	    double tauMresM[3];

 	    tauMresM[0]=tauM*resM_(0);
 	    tauMresM[1]=tauM*resM_(1);
 	    tauMresM[2]=tauM*resM_(2);

	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      int fvi   =vi*4;
	      int fvip  =fvi+1;
	      int fvipp =fvi+2;

	      temp[0][0]=fac_afgdt*(tauMresM[0]*derxy_(0,vi)+vderxyaf_(0,0)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[0][1]=fac_afgdt*(tauMresM[0]*derxy_(1,vi)+vderxyaf_(0,1)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[0][2]=fac_afgdt*(tauMresM[0]*derxy_(2,vi)+vderxyaf_(0,2)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[1][0]=fac_afgdt*(tauMresM[1]*derxy_(0,vi)+vderxyaf_(1,0)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[1][1]=fac_afgdt*(tauMresM[1]*derxy_(1,vi)+vderxyaf_(1,1)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[1][2]=fac_afgdt*(tauMresM[1]*derxy_(2,vi)+vderxyaf_(1,2)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[2][0]=fac_afgdt*(tauMresM[2]*derxy_(0,vi)+vderxyaf_(2,0)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[2][1]=fac_afgdt*(tauMresM[2]*derxy_(1,vi)+vderxyaf_(2,1)*(tauM*conv_c_af_(vi)+funct_(vi)));
	      temp[2][2]=fac_afgdt*(tauMresM[2]*derxy_(2,vi)+vderxyaf_(2,2)*(tauM*conv_c_af_(vi)+funct_(vi)));

              for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
              {
		int fui=4*ui;

                /* SUPG stabilisation --- convection (reactive parts)  */

		/* linearisation of testfunction

                                         /                            \
                                        |   n+af    /            \     |
                  +alphaF*gamma*dt*tauM*|  r     , | Dacc o nabla | v  | 
                                        |   M       \            /     |
                                         \                            /
		*/

                /*  linearised convective term in residual
		   
		   factor: +alphaF*gamma*dt*tauM

                    /                                               \
                   |    /            \   n+af    / n+af        \     |
                   |   | Dacc o nabla | u     , | c     o nabla | v  |
                   |    \            /           \             /     |
                    \                                               /
                */

		/* convection (intermediate)

                 factor:

                        +alphaF*gamma*dt


                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
		*/

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
	}
	else
	{ // no supg
	  if(cstab == Fluid3::continuity_stab_yes)
	  {
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
		
		const double sum =
		  inertia_and_conv_ui*funct_(vi)
		  +
		  fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
		  +
		  fac_visceff_afgdt_derxy1_ui*derxy_(1,vi)
		  +
		  fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
	  
		/*
		  GALERKIN inertia term (intermediate) + convection (intermediate)
		
		  factor:

                          /          \                   /                          \
                         |            |                 |  / n+af       \            |
                +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                         |            |                 |  \            /            |
                          \          /                   \                          /
		*/

		/*  GALERKIN viscous term (intermediate) */
	    
		/*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
		*/

		/* CONTINUITY stabilisation */

		/* factor: +gamma*dt*tauC
		  
                    /                          \
                   |                            |
                   | nabla o Dacc  , nabla o v  |
                   |                            |
                    \                          /
		*/

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
	  } // no supg but cstab
	  else
	  { // no supg, no cstab, just galerkin

	    for (int ui=0; ui<iel; ++ui) // loop columns 
	    {
	      const int fui  =4*ui;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
	    	    
	      /* GALERKIN inertia term (intermediate) + convection (intermediate) */

	      const double inertia_and_conv_ui
		= fac_alphaM*funct_(ui)+fac_afgdt*conv_c_af_(ui);

	      /*  GALERKIN viscous term (intermediate), diagonal parts */
	      const double fac_visceff_afgdt_derxy0_ui=fac_visceff_afgdt*derxy_(0,ui);
	      const double fac_visceff_afgdt_derxy1_ui=fac_visceff_afgdt*derxy_(1,ui);
	      const double fac_visceff_afgdt_derxy2_ui=fac_visceff_afgdt*derxy_(2,ui);

	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	      {
		const int fvi    =4*vi;
		const int fvip   =fvi+1;
		const int fvipp  =fvi+2;
		
		const double sum =
		  inertia_and_conv_ui*funct_(vi)
		  +
		  fac_visceff_afgdt_derxy0_ui*derxy_(0,vi)
		  +
		  fac_visceff_afgdt_derxy1_ui*derxy_(1,vi)
		  +
		  fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
		
		/*
		  GALERKIN inertia term (intermediate) + convection (intermediate)

		  factor:

                          /          \                   /                          \
                         |            |                 |  / n+af       \            |
                +alphaM *|  Dacc , v  |+alphaF*gamma*dt*| | c    o nabla | Dacc , v  |
                         |            |                 |  \            /            |
                          \          /                   \                          /
		*/

		/* GALERKIN viscous term (intermediate) */
	    
		/*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
		*/

		elemat(fvi  ,fui  ) += sum+fac_visceff_afgdt_derxy0_ui*derxy_(0,vi);
		elemat(fvi  ,fuip ) +=     fac_visceff_afgdt_derxy0_ui*derxy_(1,vi);
		elemat(fvi  ,fuipp) +=     fac_visceff_afgdt_derxy0_ui*derxy_(2,vi);
		elemat(fvip ,fui  ) +=     fac_visceff_afgdt_derxy1_ui*derxy_(0,vi);
		elemat(fvip ,fuip ) += sum+fac_visceff_afgdt_derxy1_ui*derxy_(1,vi);
		elemat(fvip ,fuipp) +=     fac_visceff_afgdt_derxy1_ui*derxy_(2,vi);
		elemat(fvipp,fui  ) +=     fac_visceff_afgdt_derxy2_ui*derxy_(0,vi);
		elemat(fvipp,fuip ) +=     fac_visceff_afgdt_derxy2_ui*derxy_(1,vi);
		elemat(fvipp,fuipp) += sum+fac_visceff_afgdt_derxy2_ui*derxy_(2,vi);
	      } // vi
	    } // ui
	  } // no supg, no cstab, just galerkin
	  
	  //---------------------------------------------------------------
	  //
	  //                      GALERKIN PART
	  //    REACTIVE TYPE LINEARISATIONS, QUASISTATIC FORMULATION
	  //
	  //---------------------------------------------------------------
	  
	  if (newton==Fluid3::Newton)
	  {
	    double fac_afgdt_vderxyaf_0_0_funct_ui=fac_afgdt*vderxyaf_(0,0);
	    double fac_afgdt_vderxyaf_0_1_funct_ui=fac_afgdt*vderxyaf_(0,1);
	    double fac_afgdt_vderxyaf_0_2_funct_ui=fac_afgdt*vderxyaf_(0,2);
	    double fac_afgdt_vderxyaf_1_0_funct_ui=fac_afgdt*vderxyaf_(1,0);
	    double fac_afgdt_vderxyaf_1_1_funct_ui=fac_afgdt*vderxyaf_(1,1);
	    double fac_afgdt_vderxyaf_1_2_funct_ui=fac_afgdt*vderxyaf_(1,2);
	    double fac_afgdt_vderxyaf_2_0_funct_ui=fac_afgdt*vderxyaf_(2,0);
	    double fac_afgdt_vderxyaf_2_1_funct_ui=fac_afgdt*vderxyaf_(2,1);
	    double fac_afgdt_vderxyaf_2_2_funct_ui=fac_afgdt*vderxyaf_(2,2);
	    
	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	    {
	      const int fui  =4*ui;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
	      
	      fac_afgdt_vderxyaf_0_0_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_0_1_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_0_2_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_1_0_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_1_1_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_1_2_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_2_0_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_2_1_funct_ui*=funct_(ui);
	      fac_afgdt_vderxyaf_2_2_funct_ui*=funct_(ui);

	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	      {
		int fvi  =4*vi;
		
		/* convection (intermediate)

		   factor:

                        +alphaF*gamma*dt


                         /                            \
                        |  /            \   n+af       |
                        | | Dacc o nabla | u      , v  |
                        |  \            /              |
                         \                            /
		*/

		elemat(fvi  ,fui  ) += fac_afgdt_vderxyaf_0_0_funct_ui*funct_(vi);
		elemat(fvi  ,fuip ) += fac_afgdt_vderxyaf_0_1_funct_ui*funct_(vi);
		elemat(fvi++,fuipp) += fac_afgdt_vderxyaf_0_2_funct_ui*funct_(vi);
		elemat(fvi  ,fui  ) += fac_afgdt_vderxyaf_1_0_funct_ui*funct_(vi);
		elemat(fvi  ,fuip ) += fac_afgdt_vderxyaf_1_1_funct_ui*funct_(vi);
		elemat(fvi++,fuipp) += fac_afgdt_vderxyaf_1_2_funct_ui*funct_(vi);
		elemat(fvi  ,fui  ) += fac_afgdt_vderxyaf_2_0_funct_ui*funct_(vi);
		elemat(fvi  ,fuip ) += fac_afgdt_vderxyaf_2_1_funct_ui*funct_(vi);
		elemat(fvi  ,fuipp) += fac_afgdt_vderxyaf_2_2_funct_ui*funct_(vi);
	      } // end loop rows (test functions for matrix)
	    } // end loop rows (solution for matrix, test function for vector)
	  } // end newton
	} // no supg

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

	    /* continuity equation (implicit, transposed of above equation) */

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
	  } // vi rows (Dp) and columns (q)
	} // ui columns (Dp) and rows (q)

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

	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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
                           
	      
	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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
	      }
            }
          } // this is a higher order element and linearisation is not minimal
	  else
	  { // either this ain't a higher order element or a
	    // linearisation of the viscous term is not necessary
	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	    {
	      const int fui  =ui*4 ;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
	      
	      const double fac_tauMp_inertia_and_conv=fac_tauMp*(alphaM*funct_(ui)+afgdt*conv_c_af_(ui));
	      
	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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
            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
            {
              int vidx = vi*4 + 3;
              double v1 = derxy_(0,vi)*vderxyaf_(0,0) + derxy_(1,vi)*vderxyaf_(1,0) + derxy_(2,vi)*vderxyaf_(2,0);
              double v2 = derxy_(0,vi)*vderxyaf_(0,1) + derxy_(1,vi)*vderxyaf_(1,1) + derxy_(2,vi)*vderxyaf_(2,1);
              double v3 = derxy_(0,vi)*vderxyaf_(0,2) + derxy_(1,vi)*vderxyaf_(1,2) + derxy_(2,vi)*vderxyaf_(2,2);
              for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
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

	  for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	  {
	    const int fuippp=ui*4+3;

	    const double fac_gamma_dt_tauMp_derxy_0_ui=fac_gamma_dt_tauMp*derxy_(0,ui);
	    const double fac_gamma_dt_tauMp_derxy_1_ui=fac_gamma_dt_tauMp*derxy_(1,ui);
	    const double fac_gamma_dt_tauMp_derxy_2_ui=fac_gamma_dt_tauMp*derxy_(2,ui);

	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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
        //             SUPG PART, QUASISTATIC FORMULATION
        //
        //---------------------------------------------------------------
        if(supg == Fluid3::convective_stab_supg)
        {
          const double fac_tauM               = fac*tauM;
	  const double fac_tauM_gamma_dt      = fac*tauM*gamma*dt;
	  const double fac_tauM_afgdt         = fac_tauM*afgdt;
          const double fac_visceff_afgdt_tauM = fac_tauM_afgdt*visceff;


	  for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	  {
	    const int fuippp=4*ui+3;
	    
	    const double fac_tauM_gamma_dt_derxy_0_ui = fac_tauM_gamma_dt*derxy_(0,ui);
	    const double fac_tauM_gamma_dt_derxy_1_ui = fac_tauM_gamma_dt*derxy_(1,ui);
	    const double fac_tauM_gamma_dt_derxy_2_ui = fac_tauM_gamma_dt*derxy_(2,ui);

	    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	    {
	      int fvi =vi*4;
  
              /* SUPG stabilisation --- pressure    */

              /* factor: +tauM, rescaled by gamma*dt

                    /                               \
                   |              / n+af       \     |
                   |  nabla Dp , | c    o nabla | v  |
                   |              \            /     |
                    \                               /
              */

              elemat(fvi++,fuippp) += fac_tauM_gamma_dt_derxy_0_ui*conv_c_af_(vi);
              elemat(fvi++,fuippp) += fac_tauM_gamma_dt_derxy_1_ui*conv_c_af_(vi);
              elemat(fvi  ,fuippp) += fac_tauM_gamma_dt_derxy_2_ui*conv_c_af_(vi);

            } // end loop rows (test functions for matrix)
          } // end loop rows (solution for matrix, test function for vector)

          if (higher_order_ele && newton!=Fluid3::minimal)
          {
	    for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
	    {
	      const int fui  =ui*4 ;
	      const int fuip =fui+1;
	      const int fuipp=fui+2;
   
	      /* SUPG stabilisation --- diffusion   */
	      
	      /* factor: -nu*alphaF*gamma*dt*tauM


                    /                                              \
                   |                 /     \    / n+af        \     |
                   |  2*nabla o eps | Dacc  |, | c     o nabla | v  |
                   |                 \     /    \             /     |
                    \                                              /

	      */
	      const double fac_visceff_afgdt_tauM_viscs2_0_ui=fac_visceff_afgdt_tauM*viscs2_(0,ui);
	      const double fac_visceff_afgdt_tauM_viscs2_1_ui=fac_visceff_afgdt_tauM*viscs2_(1,ui);
	      const double fac_visceff_afgdt_tauM_viscs2_2_ui=fac_visceff_afgdt_tauM*viscs2_(2,ui);
	      const double fac_visceff_afgdt_tauM_derxy2_3_ui=fac_visceff_afgdt_tauM*derxy2_(3,ui);
	      const double fac_visceff_afgdt_tauM_derxy2_4_ui=fac_visceff_afgdt_tauM*derxy2_(4,ui);
	      const double fac_visceff_afgdt_tauM_derxy2_5_ui=fac_visceff_afgdt_tauM*derxy2_(5,ui);

	      for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
	      {
		int fvi  =vi*4 ;
		int fvip =fvi+1;
		int fvipp=fvi+2;

                elemat(fvi  ,fui  ) -= fac_visceff_afgdt_tauM_viscs2_0_ui*conv_c_af_(vi);
		elemat(fvi  ,fuip ) -= fac_visceff_afgdt_tauM_derxy2_3_ui*conv_c_af_(vi);
                elemat(fvi  ,fuipp) -= fac_visceff_afgdt_tauM_derxy2_4_ui*conv_c_af_(vi);
                elemat(fvip ,fui  ) -= fac_visceff_afgdt_tauM_derxy2_3_ui*conv_c_af_(vi);
                elemat(fvip ,fuip ) -= fac_visceff_afgdt_tauM_viscs2_1_ui*conv_c_af_(vi);
                elemat(fvip ,fuipp) -= fac_visceff_afgdt_tauM_derxy2_5_ui*conv_c_af_(vi);
                elemat(fvipp,fui  ) -= fac_visceff_afgdt_tauM_derxy2_4_ui*conv_c_af_(vi);
                elemat(fvipp,fuip ) -= fac_visceff_afgdt_tauM_derxy2_5_ui*conv_c_af_(vi);
                elemat(fvipp,fuipp) -= fac_visceff_afgdt_tauM_viscs2_2_ui*conv_c_af_(vi);
              } //end ui
            } // end vi
          }// end higher_order_ele and linearisation of viscous term
        } // end supg

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

              } // end loop rows (test functions for matrix)
            } // end loop columns (solution for matrix, test function for vector)

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
                } // end loop rows (test functions for matrix)
              } // end loop columns (solution for matrix, test function for vector)
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

          for (int ui=0; ui<iel; ++ui) // loop columns (solution for matrix, test function for vector)
          {
	    const int fui   =4*ui;
	    const int fuip  =fui+1;
	    const int fuipp =fui+2;

            const double fac_afgdt_tauM_conv_resM_ui = fac_afgdt_tauM*conv_resM_(ui);

            for (int vi=0; vi<iel; ++vi)  // loop rows (test functions for matrix)
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
            } // end loop rows (test functions for matrix)
          } // end loop columns (solution for matrix, test function for vector)
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
      
      /*
	factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /

      */
      double fac_prenp   = fac*prenp_;
      
      if(cstab == Fluid3::continuity_stab_yes)
      {
	// continuity stabilisation adds a small-scale pressure
	
	/* factor: +tauC

                  /                          \
                 |           n+1              |
                 |  nabla o u    , nabla o v  |
                 |                            |
                  \                          /
	*/

	fac_prenp -= fac*tauC*divunp;

      } // end cstab

      
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
	
      if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
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
      const double fac_divunp  = fac*divunp;

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
      }

      if(fssgv == Fluid3::fssgv_scale_similarity ||
         fssgv == Fluid3::fssgv_mixed_Smagorinsky_all ||
         fssgv == Fluid3::fssgv_mixed_Smagorinsky_small)
      {
        //----------------------------------------------------------------------
        //     SCALE-SIMILARITY TERM (ON RIGHT HAND SIDE)

        for (int ui=0; ui<iel; ++ui)
        {
          elevec(ui*4    ) -= fac*(csconvintaf_(0) - convaf_s_(0))*funct_(ui);
          elevec(ui*4 + 1) -= fac*(csconvintaf_(1) - convaf_s_(1))*funct_(ui);
          elevec(ui*4 + 2) -= fac*(csconvintaf_(2) - convaf_s_(2))*funct_(ui);
        }
      }

      if(fssgv != Fluid3::fssgv_no && fssgv != Fluid3::fssgv_scale_similarity)
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
                                      +    derxy_(2, ui)*fsvderxyaf_(2, 0)) ;
          elevec(ui*4 + 1) -= vartfac*(    derxy_(0, ui)*fsvderxyaf_(0, 1)
                                      +    derxy_(0, ui)*fsvderxyaf_(1, 0)
                                      +2.0*derxy_(1, ui)*fsvderxyaf_(1, 1)
                                      +    derxy_(2, ui)*fsvderxyaf_(1, 2)
                                      +    derxy_(2, ui)*fsvderxyaf_(2, 1)) ;
          elevec(ui*4 + 2) -= vartfac*(    derxy_(0, ui)*fsvderxyaf_(0, 2)
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalcRes(
     Fluid3*                                               ele,
     std::vector<blitz::Array<double,1> >&                 myknots,
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
     enum Fluid3::TurbModelAction                          turb_mod_action,
     double&                                               Cs,
     double&                                               Cs_delta_sq,
     double&                                               visceff,
     double                                                l_tau,
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
#if 0
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
  double visc = material->m.fluid->viscosity;

  //------------------------------------------------------------------
  //                      SET ELEMENT DATA
  //------------------------------------------------------------------
  // set element data
  //  const DRT::Element::DiscretizationType distype = ele->Shape();

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
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

  // set gauss point coordinates
  blitz::Array<double, 1> gp(3);
  
  gp(0)=intpoints_onepoint.qxg[0][0];
  gp(1)=intpoints_onepoint.qxg[0][1];
  gp(2)=intpoints_onepoint.qxg[0][2];
  
  const double wquad = intpoints_onepoint.qwgt[0];

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

  // get element type constant for tau
  double mk=0.0;
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

  // get Jacobian matrix and determinant
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
      for (int mm=1;mm<iel;++mm)
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
    for (int rr=0;rr<iel;++rr)
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
    for (int mm=1;mm<iel;++mm)
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
      for (int mm=1;mm<iel;++mm)
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
    for (int mm=1;mm<iel;++mm)
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
  {
    CalVisc( material, visc);
  }

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

    visceff = visc + Cs_delta_sq * rateofstrain;
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

    visceff = visc + Cs_delta_sq * rateofstrain;

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
    LINALG::FixedSizeSerialDenseMatrix<3,1>  gp;

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
        DRT::UTILS::shape_function_3D_deriv2(deriv2_,gp(0),gp(1),gp(2),distype);
      }
    }
    else
    {
      if (higher_order_ele)
      {
        DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv_deriv2
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
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        xjm_(nn,rr)=deriv_(nn,0)*xyze_(rr,0);
        for (int mm=1;mm<iel;++mm)
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
      for (int rr=0;rr<iel;++rr)
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

          for (int mm=1;mm<iel;++mm)
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
        for (int rr=0;rr<iel;++rr)
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
        for (int rr=0;rr<iel;++rr)
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
      solver.GETRS('N',6,iel,bm.data(),6,&(pivot[0]),derxy2_.data(),6,&ierr);

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
      for (int mm=1;mm<iel;++mm)
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
      for (int mm=1;mm<iel;++mm)
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
      for (int mm=1;mm<iel;++mm)
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
      for (int mm=1;mm<iel;++mm)
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
        for (int mm=1;mm<iel;++mm)
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
        for (int mm=1;mm<iel;++mm)
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
        for (int mm=0;mm<iel;++mm)
        {
          aleconvintaf_(rr)-=funct_(mm)*egridvaf(rr,mm);
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
      /*--- viscous term  2* grad * epsilon(u): --------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
           |                                                |
           |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
           |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                             */


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

	viscaf_old_(0) += 
	      (viscs2_(0,mm)*evelaf(0,mm)
	       +
	       derxy2_(3,mm)*evelaf(1,mm)
	       +
	       derxy2_(4,mm)*evelaf(2,mm));
	viscaf_old_(1) += 
	      (derxy2_(3,mm)*evelaf(0,mm)
	       +
	       viscs2_(1,mm)*evelaf(1,mm)
	       +
	       derxy2_(5,mm)*evelaf(2,mm));
	viscaf_old_(2) +=
	      (derxy2_(4,mm)*evelaf(0,mm)
	       +
	       derxy2_(5,mm)*evelaf(1,mm)
	       +
	       viscs2_(2,mm)*evelaf(2,mm));
      }
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
      for (int rr=0;rr<3;++rr)
      {
        resM_(rr) -= visceff*viscaf_old_(rr);
      }
    }

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

#endif
  return;
}




/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::GetNodalBodyForce(Fluid3* ele, const double time)
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
    for(int isd=0;isd<3;isd++)
    {
      const double value=(*onoff)[isd]*(*val)[isd]*curvefac;

      for (int jnode=0; jnode<iel; jnode++)
      {
        edeadaf_(isd,jnode) = value;
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
// calculates material viscosity   u.may 05/08
//
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3GenalphaResVMM<distype>::CalVisc(
  const struct _MATERIAL*                 material,
  double&                           	  visc)
{

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index

  // compute shear rate
  double rateofshear = 0.0;
  blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());   // strain rate tensor
  for(int rr=0;rr<3;rr++)
    for(int mm=0;mm<3;mm++)
      epsilon(rr,mm)=0.5*(vderxyaf_(rr,mm)+vderxyaf_(mm,rr));

  for(int rr=0;rr<3;rr++)
    for(int mm=0;mm<3;mm++)
      rateofshear += epsilon(rr,mm)*epsilon(rr,mm);

  rateofshear = sqrt(2.0*rateofshear);

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
}



#endif
#endif
