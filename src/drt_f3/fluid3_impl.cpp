/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl.cpp

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Epetra_SerialDenseSolver.h>

#ifdef DEBUG
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplInterface* DRT::ELEMENTS::Fluid3ImplInterface::Impl(DRT::ELEMENTS::Fluid3* f3)
{
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Fluid3Impl<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new Fluid3Impl<DRT::Element::hex8>();
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static Fluid3Impl<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new Fluid3Impl<DRT::Element::hex20>();
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static Fluid3Impl<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new Fluid3Impl<DRT::Element::hex27>();
    return fh27;
  }
  case DRT::Element::tet4:
  {
    static Fluid3Impl<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new Fluid3Impl<DRT::Element::tet4>();
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static Fluid3Impl<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new Fluid3Impl<DRT::Element::tet10>();
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static Fluid3Impl<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new Fluid3Impl<DRT::Element::wedge6>();
    return fw6;
  }
  case DRT::Element::wedge15:
  {
    static Fluid3Impl<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new Fluid3Impl<DRT::Element::wedge15>();
    return fw15;
  }
  case DRT::Element::pyramid5:
  {
    static Fluid3Impl<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new Fluid3Impl<DRT::Element::pyramid5>();
    return fp5;
  }

  default:
    dserror("shape %d (%d nodes) not supported", f3->Shape(), f3->NumNode());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3Impl<distype>::Fluid3Impl()
  : vart_(),
    xyze_(),
    edeadng_(),
    funct_(),
    densfunct_(),
    densamfunct_(),
    functdens_(),
    deriv_(),
    deriv2_(),
    xjm_(),
    xji_(),
    vderxy_(),
    mderxy_(),
    fsvderxy_(),
    derxy_(),
    densderxy_(),
    derxy2_(),
    bodyforce_(),
    histmom_(),
    histcon_(),
    velino_(),
    velint_(),
    fsvelint_(),
    convvelint_(),
    accintam_(),
    gradp_(),
    tau_(),
    viscs2_(),
    conv_c_(),
    mdiv_(),
    vdiv_(),
    rhsmom_(),
    rhscon_(),
    conv_old_(),
    visc_old_(),
    res_old_(true),  // initialize to zero
    conv_resM_(true),
    xder2_(),
    vderiv_()
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3Impl<distype>::Evaluate(
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
  // the number of nodes
  const int numnode = iel;

  // construct views
  LINALG::Matrix<4*iel,4*iel> elemat1(elemat1_epetra,true);
  LINALG::Matrix<4*iel,4*iel> elemat2(elemat2_epetra,true);
  LINALG::Matrix<4*iel,    1> elevec1(elevec1_epetra,true);
  // elevec2 and elevec3 are never used anyway

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // check whether we have a generalized-alpha time-integration scheme
  const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  const double time = params.get<double>("total time",-1.0);

  // get time-step size
  const double dt = params.get<double>("dt");

  // get timefactor for left hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");
  // initially set timefactor for right hand side (only required for
  // generalized-alpha scheme in low-Mach-number case)
  // -> will be specified below if required
  double timefacrhs=0.0;

  // ---------------------------------------------------------------------
  // get control parameters for linearization and low-Mach-number solver
  //----------------------------------------------------------------------
  string newtonstr=params.get<string>("Linearisation");
  string lomastr  =params.get<string>("low-Mach-number solver");
  bool newton = false;
  bool loma   = false;
  if(newtonstr=="Newton") newton=true;
  if(lomastr  =="Yes"   ) loma  =true;

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  ParameterList& stablist = params.sublist("STABILIZATION");

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
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  // ---------------------------------------------------------------------
  // get all general state vectors: vel./press., vel./dens. and hist
  // vel./press., vel./dens. values are at time n+alpha_F for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  RefCountPtr<const Epetra_Vector> velnp  = discretization.GetState("velnp");
  RefCountPtr<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  RefCountPtr<const Epetra_Vector> hist   = discretization.GetState("hist");
  if (velnp==null || vedenp==null || hist==null)
    dserror("Cannot get state vectors 'velnp', 'vedenp' and/or 'hist'");

  // extract local values from the global vectors
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  vector<double> myvedenp(lm.size());
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);
  vector<double> myhist(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eprenp;
  LINALG::Matrix<3,numnode> evelnp;
  LINALG::Matrix<numnode,1> edensnp;
  LINALG::Matrix<3,numnode> emhist;
  LINALG::Matrix<numnode,1> echist;

  for (int i=0;i<numnode;++i)
  {
    // split velocity and pressure/density, insert into element arrays
    evelnp(0,i) = myvelnp[0+(i*4)];
    evelnp(1,i) = myvelnp[1+(i*4)];
    evelnp(2,i) = myvelnp[2+(i*4)];

    eprenp(i) = myvelnp[3+(i*4)];

    edensnp(i) = myvedenp[3+(i*4)];

    // momentum equation part of history vector
    // (containing information of time step t_n (mass rhs!))
    emhist(0,i) = myhist[0+(i*4)];
    emhist(1,i) = myhist[1+(i*4)];
    emhist(2,i) = myhist[2+(i*4)];

    // continuity equation part (only non-trivial for low-Mach-number flow)
    echist(i) = myhist[3+(i*4)];
  }

  // ---------------------------------------------------------------------
  // get additional state vectors for generalized-alpha scheme:
  // vel./dens. and acceleration/density time derivative at time n+alpha_M
  // ---------------------------------------------------------------------
  RefCountPtr<const Epetra_Vector> vedeam;
  vector<double> myvedeam;
  RefCountPtr<const Epetra_Vector> accam;
  vector<double> myaccam;

  // create objects for element arrays
  LINALG::Matrix<3,numnode> eaccam;
  LINALG::Matrix<numnode,1> ededtam;
  LINALG::Matrix<numnode,1> edensam;

  if (is_genalpha)
  {
    vedeam = discretization.GetState("vedeam");
    accam  = discretization.GetState("accam");
    if (vedeam==null || accam==null)
      dserror("Cannot get state vectors 'vedeam' and/or 'accam'");

    // extract local values from the global vectors
    myvedeam.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*vedeam,myvedeam,lm);
    myaccam.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);

    for (int i=0;i<numnode;++i)
    {
      // split acceleration and density time derivative
      eaccam(0,i) = myaccam[0+(i*4)];
      eaccam(1,i) = myaccam[1+(i*4)];
      eaccam(2,i) = myaccam[2+(i*4)];

      ededtam(i)  = myaccam[3+(i*4)];

      // extract density
      edensam(i)  = myvedeam[3+(i*4)];
    }

    // set timefactor for right hand side (only required for
    // low-Mach-number case)
    timefacrhs = params.get<double>("timefacrhs",-1.0);
    if (timefacrhs < 0.0) dserror("No timefacrhs supplied");
  }
  else
  {
    // ensure that density at time n+alpha_M is correctly set for
    // all schemes other than generalized-alpha
    for (int i=0;i<numnode;++i)
    {
      edensam(i,0) = myvedenp[3+(i*4)];
    }
  }

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  RefCountPtr<const Epetra_Vector> dispnp;
  vector<double> mydispnp;
  RefCountPtr<const Epetra_Vector> gridv;
  vector<double> mygridv;

  // create objects for element arrays
  LINALG::Matrix<3,numnode> edispnp;
  LINALG::Matrix<3,numnode> egridv;

  if (ele->is_ale_)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vectors 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);

    for (int i=0;i<numnode;++i)
    {
      // set grid displacements
      edispnp(0,i) = mydispnp[0+(i*4)];
      edispnp(1,i) = mydispnp[1+(i*4)];
      edispnp(2,i) = mydispnp[2+(i*4)];

      // set grid velocities
      egridv(0,i) = mygridv[0+(i*4)];
      egridv(1,i) = mygridv[1+(i*4)];
      egridv(2,i) = mygridv[2+(i*4)];
    }
  }

  // ---------------------------------------------------------------------
  // get additional state vector for AVM3 case: fine-scale velocity
  // values are at time n+alpha_F for generalized-alpha scheme and at
  // time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // get flag for fine-scale subgrid-viscosity approach
  Fluid3::StabilisationAction fssgv =
    ele->ConvertStringToStabAction(params.get<string>("fs subgrid viscosity","No"));

  RCP<const Epetra_Vector> fsvelnp;
  vector<double> myfsvelnp;

  // create object for element array
  LINALG::Matrix<3,numnode> fsevelnp;

  if (fssgv != Fluid3::fssgv_no)
  {
    fsvelnp = discretization.GetState("fsvelnp");
    if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
    myfsvelnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

    for (int i=0;i<numnode;++i)
    {
      // get fine-scale velocity
      fsevelnp(0,i) = myfsvelnp[0+(i*4)];
      fsevelnp(1,i) = myfsvelnp[1+(i*4)];
      fsevelnp(2,i) = myfsvelnp[2+(i*4)];
    }
  }

  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
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

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  nlayer=0;

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action = Fluid3::smagorinsky;
      Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    // --------------------------------------------------
    // Smagorinsky model with van Driest damping
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      // that's only implemented for turbulent channel flow
      if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
          !=
          "channel_flow_of_height_2")
      {
        dserror("van_Driest_damping only for channel_flow_of_height_2\n");
      }

      // for the Smagorinsky model with van Driest damping, we need
      // a viscous length to determine the y+ (heigth in wall units)
      turb_mod_action = Fluid3::smagorinsky_with_wall_damping;

      // get parameters of model
      Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
      l_tau           = turbmodelparams.get<double>("CHANNEL_L_TAU");

      // this will be the y-coordinate of a point in the element interior
      // we will determine the element layer in which he is contained to
      // be able to do the output of visceff etc.
      double center = 0;

      DRT::Node** nodes = ele->Nodes();
      for(int inode=0;inode<numnode;inode++)
      {
        center+=nodes[inode]->X()[1];
      }
      center/=numnode;

      // node coordinates of plane to the element layer
      RefCountPtr<vector<double> > planecoords
        =
        turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");

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
    // --------------------------------------------------
    // Smagorinsky model with dynamic Computation of Cs
    else if (physical_turbulence_model == "Dynamic_Smagorinsky")
    {
      turb_mod_action = Fluid3::dynamic_smagorinsky;

      // for turbulent channel flow, use averaged quantities
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

        //this will be the y-coordinate of a point in the element interior
        // here, the layer is determined in order to get the correct
        // averaged value from the vector of averaged (M/L)ijMij
        double center = 0;
        DRT::Node** nodes = ele->Nodes();
        for(int inode=0;inode<numnode;inode++)
        {
          center+=nodes[inode]->X()[1];
        }
        center/=numnode;

        RCP<vector<double> > planecoords
          =
          turbmodelparams.get<RCP<vector<double> > >("planecoords_");

        bool found = false;
        for (nlayer=0;nlayer < static_cast<int>((*planecoords).size()-1);)
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

        // Cs_delta_sq is set by the averaged quantities
        Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;

        // clipping to get algorithm stable
        if (Cs_delta_sq<0)
        {
          Cs_delta_sq=0;
        }
      }
      else
      {
        // when no averaging was done, we just keep the calculated (clipped) value
        Cs_delta_sq = ele->Cs_delta_sq_;
      }
    }
    else
    {
      dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
    }
  }

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele,
         evelnp,
         fsevelnp,
         eprenp,
         edensnp,
         eaccam,
         ededtam,
         edensam,
         emhist,
         echist,
         edispnp,
         egridv,
         elemat1,
         elemat2,
         elevec1,
         actmat,
         time,
         dt,
         timefac,
         timefacrhs,
         newton,
         loma,
         is_genalpha,
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
         l_tau);

  // ---------------------------------------------------------------------
  // output values of Cs, visceff and Cs_delta_sq
  // ---------------------------------------------------------------------
  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    if (physical_turbulence_model == "Dynamic_Smagorinsky"
        ||
        physical_turbulence_model ==  "Smagorinsky_with_van_Driest_damping"
      )
    {
      if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
          ==
          "channel_flow_of_height_2")
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

  return 0;
}



/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Sysmat(
  Fluid3*                                 ele,
  const LINALG::Matrix<3,iel>&            evelnp,
  const LINALG::Matrix<3,iel>&            fsevelnp,
  const LINALG::Matrix<iel,1>&            eprenp,
  const LINALG::Matrix<iel,1>&            edensnp,
  const LINALG::Matrix<3,iel>&            eaccam,
  const LINALG::Matrix<iel,1>&            ededtam,
  const LINALG::Matrix<iel,1>&            edensam,
  const LINALG::Matrix<3,iel>&            emhist,
  const LINALG::Matrix<iel,1>&            echist,
  const LINALG::Matrix<3,iel>&            edispnp,
  const LINALG::Matrix<3,iel>&            egridv,
  LINALG::Matrix<4*iel,4*iel>&            estif,
  LINALG::Matrix<4*iel,4*iel>&            emesh,
  LINALG::Matrix<4*iel,1>&                eforce,
  struct _MATERIAL*                       material,
  double                                  time,
  double                                  dt,
  double                                  timefac,
  double                                  timefacrhs,
  bool                                    newton,
  bool                                    loma,
  const bool                              is_genalpha,
  const bool                              higher_order_ele,
  const enum Fluid3::StabilisationAction  fssgv,
  const enum Fluid3::StabilisationAction  pspg,
  const enum Fluid3::StabilisationAction  supg,
  const enum Fluid3::StabilisationAction  vstab,
  const enum Fluid3::StabilisationAction  cstab,
  const enum Fluid3::StabilisationAction  cross,
  const enum Fluid3::StabilisationAction  reynolds,
  const enum Fluid3::TauType              whichtau,
  const enum Fluid3::TurbModelAction      turb_mod_action,
  double&                                 Cs,
  double&                                 Cs_delta_sq,
  double&                                 visceff,
  double&                                 l_tau
  )
{
  // set element data
  const int numnode = iel;

  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }

  // add displacement when fluid nodes move in the ALE case
  if (ele->is_ale_) xyze_ += edispnp;

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->mattyp != m_fluid
   && material->mattyp != m_carreauyasuda
   && material->mattyp != m_modpowerlaw) dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(material->mattyp == m_fluid) visc = material->m.fluid->viscosity;

  // ---------------------------------------------------------------------
  // call routine for calculation of stabilization parameter
  // (needs to be done before anything else is calculated, since
  //  we use the same arrays internally)
  // ---------------------------------------------------------------------
  Caltau(ele,
         evelnp,
         fsevelnp,
         edensnp,
         whichtau,
         material,
         visc,
         timefac,
         dt,
         turb_mod_action,
         Cs,
         Cs_delta_sq,
         visceff,
         l_tau,
         fssgv);

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

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

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
    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);

    if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

    const double fac = intpoints.qwgt[iquad]*det;

    // compute global first derivates
    derxy_.Multiply(xji_,deriv_);

    // density-weighted shape functions at n+1/n+alpha_F and n+1/n+alpha_M
    densfunct_.EMultiply(funct_,edensnp);
    densamfunct_.EMultiply(funct_,edensam);

    // inverse-density-weighted shape functions and density-weighted
    // global first derivatives at n+1/n+alpha_F
    for (int inode=0; inode<numnode; inode++)
    {
      functdens_(inode) = funct_(inode)/edensnp(inode);

      densderxy_(0,inode) = edensnp(inode,0)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode,0)*derxy_(1,inode);
      densderxy_(2,inode) = edensnp(inode,0)*derxy_(2,inode);
    }

    //--------------------------------------------------------------
    //             compute global second derivatives
    //--------------------------------------------------------------
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      gder2(ele);
    }
    else derxy2_.Clear();

    // get momentum (i.e., density times velocity) at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelnp,densfunct_);

    // get (density-weighted) history data at integration point
    histmom_.Multiply(emhist,densfunct_);
    histcon_ = funct_.Dot(echist);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelnp,derxy_);

    // get momentum derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    mderxy_.MultiplyNT(evelnp,densderxy_);

    // get fine-scale velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (fssgv != Fluid3::fssgv_no) fsvderxy_.MultiplyNT(fsevelnp,derxy_);
    else                           fsvderxy_.Clear();

    // get density-weighted convective velocity at integration point
    // We handle the ale case very implicitly here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_.Update(velint_);
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, densfunct_, 1.0);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    gradp_.Multiply(derxy_,eprenp);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = funct_.Dot(eprenp);

    // get density-weighted bodyforce in gausspoint
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(edeadng_,densfunct_);

    //--------------------------------------------------------------------
    // get numerical representation of some single operators
    //--------------------------------------------------------------------
    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_,convvelint_);

    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_,convvelint_);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) --------------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
         - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
         2 |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                             */

      /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
      /*   /                            \
           |  N_x,xx + N_y,yx + N_z,zx  |
         1 |                            |
      -  - |  N_x,xy + N_y,yy + N_z,zy  |
         3 |                            |
           |  N_x,xz + N_y,yz + N_z,zz  |
           \                            /

             with N_x .. x-line of N
             N_y .. y-line of N                                             */

      double prefac;
      if (loma)
      {
        prefac = 1.0/3.0;
        derxy2_.Scale(prefac);
      }
      else prefac = 1.0;

      double sum = (derxy2_(0,0)+derxy2_(1,0)+derxy2_(2,0))/prefac;

      viscs2_(0,0) = 0.5 * (sum + derxy2_(0,0));
      viscs2_(1,0) = 0.5 *  derxy2_(3,0);
      viscs2_(2,0) = 0.5 *  derxy2_(4,0);
      viscs2_(3,0) = 0.5 *  derxy2_(3,0);
      viscs2_(4,0) = 0.5 * (sum + derxy2_(1,0));
      viscs2_(5,0) = 0.5 *  derxy2_(5,0);
      viscs2_(6,0) = 0.5 *  derxy2_(4,0);
      viscs2_(7,0) = 0.5 *  derxy2_(5,0);
      viscs2_(8,0) = 0.5 * (sum + derxy2_(2,0));

      visc_old_(0) = viscs2_(0,0)*evelnp(0,0)+viscs2_(1,0)*evelnp(1,0)+viscs2_(2,0)*evelnp(2,0);
      visc_old_(1) = viscs2_(3,0)*evelnp(0,0)+viscs2_(4,0)*evelnp(1,0)+viscs2_(5,0)*evelnp(2,0);
      visc_old_(2) = viscs2_(6,0)*evelnp(0,0)+viscs2_(7,0)*evelnp(1,0)+viscs2_(8,0)*evelnp(2,0);

      for (int i=1; i<numnode; ++i)
      {
        double sum = (derxy2_(0,i)+derxy2_(1,i)+derxy2_(2,i))/prefac;

        viscs2_(0,i) = 0.5 * (sum + derxy2_(0,i));
        viscs2_(1,i) = 0.5 *  derxy2_(3,i);
        viscs2_(2,i) = 0.5 *  derxy2_(4,i);
        viscs2_(3,i) = 0.5 *  derxy2_(3,i);
        viscs2_(4,i) = 0.5 * (sum + derxy2_(1,i));
        viscs2_(5,i) = 0.5 *  derxy2_(5,i);
        viscs2_(6,i) = 0.5 *  derxy2_(4,i);
        viscs2_(7,i) = 0.5 *  derxy2_(5,i);
        viscs2_(8,i) = 0.5 * (sum + derxy2_(2,i));

        visc_old_(0) += viscs2_(0,i)*evelnp(0,i)+viscs2_(1,i)*evelnp(1,i)+viscs2_(2,i)*evelnp(2,i);
        visc_old_(1) += viscs2_(3,i)*evelnp(0,i)+viscs2_(4,i)*evelnp(1,i)+viscs2_(5,i)*evelnp(2,i);
        visc_old_(2) += viscs2_(6,i)*evelnp(0,i)+viscs2_(7,i)*evelnp(1,i)+viscs2_(8,i)*evelnp(2,i);
      }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }

    // momentum and velocity divergence
    // (the later is only required for low-Mach-number flow)
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1) + mderxy_(2, 2);
    if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2);

    //--------------------------------------------------------------------
    // stabilization, time-integration and subgrid-viscosity factors
    //--------------------------------------------------------------------
    const double tau_M      = tau_(0)*fac;
    const double tau_Mp     = tau_(1)*fac;
    const double tau_C      = tau_(2)*fac;

    const double timefacfac = timefac * fac;
    const double timetauM   = timefac * tau_M;
    const double timetauMp  = timefac * tau_Mp;
    double       rhsfac     = fac;

    const double vartfac    = vart_*timefacfac;

    //--------------------------------------------------------------------
    // calculation of rhs for momentum/continuity equation and residual
    // (different for generalized-alpha and other schemes)
    //--------------------------------------------------------------------
    if(is_genalpha)
    {
      // rhs of momentum equation: only bodyforce at n+alpha_F
      rhsmom_.Update(1.0,bodyforce_,0.0);

      // get time derivative of density at integration point
      double densdt = funct_.Dot(ededtam);

      // rhs of continuity equation (only relevant for low-Mach-number flow)
      rhscon_ = - timefacrhs * densdt;

      // get acceleration at time n+alpha_M at integration point
      accintam_.Multiply(eaccam,densamfunct_);

      // evaluate residual once for all stabilization right hand sides
      for (int rr=0;rr<3;++rr)
      {
        res_old_(rr) = accintam_(rr)+conv_old_(rr)+gradp_(rr)-2*visceff*visc_old_(rr)-bodyforce_(rr);
      }
    }
    else
    {
      // rhs of momentum equation: timefac*bodyforce at n+1 + histmom
      rhsmom_.Update(1.0,histmom_,timefac,bodyforce_);

      // get density at integration point
      double dens = funct_.Dot(edensnp);

      // rhs of continuity equation (only relevant for low-Mach-number flow)
      rhscon_ = histcon_ - dens;

      // modify integration factor for Galerkin rhs and continuity stabilization
      rhsfac *= timefac;

      // evaluate residual once for all stabilization right hand sides
      for (int rr=0;rr<3;++rr)
      {
        res_old_(rr) = velint_(rr)-rhsmom_(rr)+timefac*(conv_old_(rr)+gradp_(rr)-2*visceff*visc_old_(rr));
      }
    }

    /*
      This is the operator

                  /               \
                 | resM    o nabla |
                  \    (i)        /

                  required for (lhs) cross- and (rhs) Reynolds-stress calculation

    */

    if (cross    == Fluid3::cross_stress_stab ||
        reynolds == Fluid3::reynolds_stress_stab_only_rhs)
      conv_resM_.MultiplyTN (densderxy_,res_old_);

    //--------------------------------------------------------------------
    // perform integration for element matrix and right hand side
    //--------------------------------------------------------------------
    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      for (int ui=0; ui<numnode; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;
        double v = fac*densamfunct_(ui)
#if 1
                   + timefacfac*conv_c_(ui)
#endif
                   ;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 4*vi;

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
          estif(fvi, fui)       += v2;
          estif(fvi + 1, fuip)  += v2;
          estif(fvi + 2, fuipp) += v2;
        }
      }

      const double visceff_timefacfac = visceff*timefacfac;
      for (int ui=0; ui<numnode; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
          const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);
          /* viscosity term */
          /*

                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(fvi, fui)     += visceff_timefacfac*(2.0*derxy_0ui_0vi
                                                     +
                                                     derxy_1ui_1vi
                                                     +
                                                     derxy_2ui_2vi) ;
          estif(fvi , fuip)   += visceff_timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(fvi , fuipp)  += visceff_timefacfac*derxy_(0, ui)*derxy_(2, vi) ;
          estif(fvip, fui)    += visceff_timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(fvip, fuip)   += visceff_timefacfac*(derxy_0ui_0vi
                                                     +
                                                     2.0*derxy_1ui_1vi
                                                     +
                                                     derxy_2ui_2vi) ;
          estif(fvip , fuipp) += visceff_timefacfac*derxy_(1, ui)*derxy_(2, vi) ;
          estif(fvipp, fui)   += visceff_timefacfac*derxy_(0, vi)*derxy_(2, ui) ;
          estif(fvipp, fuip)  += visceff_timefacfac*derxy_(1, vi)*derxy_(2, ui) ;
          estif(fvipp, fuipp) += visceff_timefacfac*(derxy_0ui_0vi
                                                     +
                                                     derxy_1ui_1vi
                                                     +
                                                     2.0*derxy_2ui_2vi) ;

        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        const int fuippp = 4*ui+3;

        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi   = 4*vi;

          /* Druckterm */
          /*

          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */

          estif(fvi,     fuippp) += (v*derxy_(0, vi)) ;
          estif(fvi + 1, fuippp) += (v*derxy_(1, vi)) ;
          estif(fvi + 2, fuippp) += (v*derxy_(2, vi)) ;

        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const double v = timefacfac*functdens_(vi);

        const int fvippp = 4*vi+3;
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui   = 4*ui;

          /* Divergenzfreiheit */
          /*
            /                  \
            |                  |
            | nabla o Du  , q  |
            |                  |
            \                  /
          */
          estif(fvippp, fui)     += v*densderxy_(0, ui) ;
          estif(fvippp, fui + 1) += v*densderxy_(1, ui) ;
          estif(fvippp, fui + 2) += v*densderxy_(2, ui) ;
        }
      }

      if (newton)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const double v = timefacfac*funct_(vi);
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;

          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;

            const double v2 = v*densfunct_(ui);
            /*  convection, reactive part

            /                           \
            |  /          \   n+1       |
            | | Du o nabla | u     , v  |
            |  \          /   (i)       |
            \                           /
            */
            estif(fvi,   fui)   += v2*vderxy_(0, 0) ;
            estif(fvi,   fuip)  += v2*vderxy_(0, 1) ;
            estif(fvi,   fuipp) += v2*vderxy_(0, 2) ;
            estif(fvip,  fui)   += v2*vderxy_(1, 0) ;
            estif(fvip,  fuip)  += v2*vderxy_(1, 1) ;
            estif(fvip,  fuipp) += v2*vderxy_(1, 2) ;
            estif(fvipp, fui)   += v2*vderxy_(2, 0) ;
            estif(fvipp, fuip)  += v2*vderxy_(2, 1) ;
            estif(fvipp, fuipp) += v2*vderxy_(2, 2) ;
          }
        }
      }

      if(is_genalpha)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 4*vi;
          /* inertia */
          const double v = -fac*funct_(vi);
          eforce(fvi)     += v*accintam_(0) ;
          eforce(fvi + 1) += v*accintam_(1) ;
          eforce(fvi + 2) += v*accintam_(2) ;
        }
      }
      else
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 4*vi;
          /* inertia */
          const double v = -fac*funct_(vi);
          eforce(fvi)     += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
          eforce(fvi + 2) += v*velint_(2) ;
        }
      }

#if 1
      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi = 4*vi;
        /* convection */
        double v = -rhsfac*funct_(vi);
        eforce(fvi    ) += v*conv_old_(0) ;
        eforce(fvi + 1) += v*conv_old_(1) ;
        eforce(fvi + 2) += v*conv_old_(2) ;
      }
#endif

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi = 4*vi;
        /* pressure */
        const double v = press*rhsfac;
        eforce(fvi    ) += v*derxy_(0, vi) ;
        eforce(fvi + 1) += v*derxy_(1, vi) ;
        eforce(fvi + 2) += v*derxy_(2, vi) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi = 4*vi;
        const double v = -visceff*rhsfac;
        /* viscosity */
        eforce(fvi)     += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                              +
                              derxy_(1, vi)*vderxy_(0, 1)
                              +
                              derxy_(1, vi)*vderxy_(1, 0)
                              +
                              derxy_(2, vi)*vderxy_(0, 2)
                              +
                              derxy_(2, vi)*vderxy_(2, 0)) ;
        eforce(fvi + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                              +
                              derxy_(0, vi)*vderxy_(1, 0)
                              +
                              2.0*derxy_(1, vi)*vderxy_(1, 1)
                              +
                              derxy_(2, vi)*vderxy_(1, 2)
                              +
                              derxy_(2, vi)*vderxy_(2, 1)) ;
        eforce(fvi + 2) += v*(derxy_(0, vi)*vderxy_(0, 2)
                              +
                              derxy_(0, vi)*vderxy_(2, 0)
                              +
                              derxy_(1, vi)*vderxy_(1, 2)
                              +
                              derxy_(1, vi)*vderxy_(2, 1)
                              +
                              2.0*derxy_(2, vi)*vderxy_(2, 2)) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi = 4*vi;
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eforce(fvi    ) += v*rhsmom_(0) ;
        eforce(fvi + 1) += v*rhsmom_(1) ;
        eforce(fvi + 2) += v*rhsmom_(2) ;
      }

      const double rhsfac_mdiv = -rhsfac * mdiv_;
      for (int vi=0; vi<numnode; ++vi)
      {
        // continuity equation
        eforce(vi*4 + 3) += rhsfac_mdiv*functdens_(vi) ;
      }

      if (loma)
      {
        const double v = -(2.0/3.0)*visceff*timefacfac ;
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* viscosity term - subtraction for low-Mach-number flow */
            /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
            */
            estif(fvi,   fui  ) += v*derxy_(0, vi)*derxy_(0, ui) ;
            estif(fvi,   fuip ) += v*derxy_(0, vi)*derxy_(1, ui) ;
            estif(fvi,   fuipp) += v*derxy_(0, vi)*derxy_(2, ui) ;
            estif(fvip,  fui  ) += v*derxy_(1, vi)*derxy_(0, ui) ;
            estif(fvip,  fuip ) += v*derxy_(1, vi)*derxy_(1, ui) ;
            estif(fvip,  fuipp) += v*derxy_(1, vi)*derxy_(2, ui) ;
            estif(fvipp, fui  ) += v*derxy_(2, vi)*derxy_(0, ui) ;
            estif(fvipp, fuip ) += v*derxy_(2, vi)*derxy_(1, ui) ;
            estif(fvipp, fuipp) += v*derxy_(2, vi)*derxy_(2, ui) ;

          }
        }

        const double v_div = (2.0/3.0)*visceff*rhsfac*vdiv_ ;
        const double fac_rhscon = fac*rhscon_;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 4*vi;
          /* viscosity term - subtraction for low-Mach-number flow */
          eforce(fvi    ) += derxy_(0, vi)*v_div ;
          eforce(fvi + 1) += derxy_(1, vi)*v_div ;
          eforce(fvi + 2) += derxy_(2, vi)*v_div ;

          /* additional rhs term of continuity equation */
          eforce(fvi + 3) += fac_rhscon*functdens_(vi) ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg == Fluid3::pstab_use_pspg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v = tau_Mp*densamfunct_(ui)
#if 1
                     + timetauMp*conv_c_(ui)
#endif
                     ;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvippp = 4*vi+3;

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

            estif(fvippp, fui)   += v*derxy_(0, vi) ;
            estif(fvippp, fuip)  += v*derxy_(1, vi) ;
            estif(fvippp, fuipp) += v*derxy_(2, vi) ;
          }
        }

        if (higher_order_ele)
        {
          const double v = -2.0*visceff*timetauMp;
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            for (int vi=0; vi<numnode; ++vi)
            {

              const int fvippp = 4*vi+3;
              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(fvippp, fui)   += v*(derxy_(0, vi)*viscs2_(0, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(1, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(2, ui)) ;
              estif(fvippp, fuip)  += v*(derxy_(0, vi)*viscs2_(1, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(4, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(5, ui)) ;
              estif(fvippp, fuipp) += v*(derxy_(0, vi)*viscs2_(2, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(5, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(8, ui)) ;
            }
          }
        }

        for (int ui=0; ui<numnode; ++ui)
        {
          const int fuippp = 4*ui+3;
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
            estif(vi*4 + 3, fuippp) += timetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                  +
                                                  derxy_(1, ui)*derxy_(1, vi)
                                                  +
                                                  derxy_(2, ui)*derxy_(2, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauMp*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvippp = 4*vi + 3;
              /*  pressure stabilisation: convection, reactive part

              /                             \
              |  /          \   n+1           |
              | | Du o nabla | u     , grad q |
              |  \          /   (i)           |
              \                             /

              */
              estif(fvippp, fui)   += v*(derxy_(0, vi)*vderxy_(0, 0)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 0)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 0)) ;
              estif(fvippp, fuip)  += v*(derxy_(0, vi)*vderxy_(0, 1)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 1)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 1)) ;
              estif(fvippp, fuipp) += v*(derxy_(0, vi)*vderxy_(0, 2)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 2)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 2)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<numnode; ++vi)
        {
          // pressure stabilisation
          eforce(vi*4 + 3) -= tau_Mp*(res_old_(0)*derxy_(0, vi)
                                      +
                                      res_old_(1)*derxy_(1, vi)
                                      +
                                      res_old_(2)*derxy_(2, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg == Fluid3::convective_stab_supg)
      {
#if 1
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v = tau_M*densamfunct_(ui) + timetauM*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
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

            estif(fvi,   fui)   += v*conv_c_(vi);
            estif(fvip,  fuip)  += v*conv_c_(vi);
            estif(fvipp, fuipp) += v*conv_c_(vi);
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          const double v = timetauM*conv_c_(vi,0);
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fuippp = 4*ui + 3;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                              \
              |              / n+1       \     |
              |  nabla Dp , | u   o nabla | v  |
              |              \ (i)       /     |
              \                              /
            */
            estif(fvi,   fuippp) += v*derxy_(0, ui) ;
            estif(fvip,  fuippp) += v*derxy_(1, ui) ;
            estif(fvipp, fuippp) += v*derxy_(2, ui) ;
          }
        }

        if (higher_order_ele)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v = -2.0*visceff*timetauM*conv_c_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;

              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                        \
                |               /  \    / n+1        \     |
                |  nabla o eps | Du |, | u    o nabla | v  |
                |               \  /    \ (i)        /     |
                \                                        /
              */
              estif(fvi, fui)     += v*viscs2_(0, ui) ;
              estif(fvip, fui)    += v*viscs2_(1, ui) ;
              estif(fvipp, fui)   += v*viscs2_(2, ui) ;

              estif(fvi, fuip)    += v*viscs2_(1, ui) ;
              estif(fvip, fuip)   += v*viscs2_(4, ui) ;
              estif(fvipp, fuip)  += v*viscs2_(5, ui) ;

              estif(fvi, fuipp)   += v*viscs2_(2, ui) ;
              estif(fvip, fuipp)  += v*viscs2_(5, ui) ;
              estif(fvipp, fuipp) += v*viscs2_(8, ui) ;
            }
          }
        }
#endif

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = tau_M*densamfunct_(ui);
            const double v0 = v*velint_(0);
            const double v1 = v*velint_(1);
            const double v2 = v*velint_(2);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;

              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                /                            \
                |   n+1      /          \     |
                |  u      , | Du o nabla | v  |
                |   (i)      \          /     |
                \                            /

              */
              estif(fvi,  fui)     += v0*derxy_(0, vi) ;
              estif(fvip,  fui)    += v1*derxy_(0, vi) ;
              estif(fvipp, fui)    += v2*derxy_(0, vi) ;

              estif(fvi,   fuip)   += v0*derxy_(1, vi) ;
              estif(fvip,  fuip)   += v1*derxy_(1, vi) ;
              estif(fvipp, fuip)   += v2*derxy_(1, vi) ;

              estif(fvi,   fuipp)  += v0*derxy_(2, vi) ;
              estif(fvip,  fuipp)  += v1*derxy_(2, vi) ;
              estif(fvipp, fuipp)  += v2*derxy_(2, vi) ;
            }
          }

#if 1
          {
            const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1) + convvelint_(2)*vderxy_(0, 2);
            const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1) + convvelint_(2)*vderxy_(1, 2);
            const double v2 = convvelint_(0)*vderxy_(2, 0) + convvelint_(1)*vderxy_(2, 1) + convvelint_(2)*vderxy_(2, 2);

            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = timetauM*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;

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
                estif(fvi, fui)     += (conv_c_(vi,0)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(fvip, fui)    += (conv_c_(vi,0)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;
                estif(fvipp, fui)   += (conv_c_(vi,0)*vderxy_(2, 0) + v2*derxy_(0, vi))*v;

                estif(fvi, fuip)    += (conv_c_(vi,0)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(fvip, fuip)   += (conv_c_(vi,0)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
                estif(fvipp, fuip)  += (conv_c_(vi,0)*vderxy_(2, 1) + v2*derxy_(1, vi))*v;

                estif(fvi, fuipp)   += (conv_c_(vi,0)*vderxy_(0, 2) + v0*derxy_(2, vi))*v;
                estif(fvip, fuipp)  += (conv_c_(vi,0)*vderxy_(1, 2) + v1*derxy_(2, vi))*v;
                estif(fvipp, fuipp) += (conv_c_(vi,0)*vderxy_(2, 2) + v2*derxy_(2, vi))*v;
              }
            }
          }
#endif

          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauM*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;

              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                               \
                |         n+1    /          \     |
                |  nabla p    , | Du o nabla | v  |
                |         (i)    \          /     |
                \                               /
              */
              estif(fvi, fui)     += v*gradp_(0,0)*derxy_(0, vi) ;
              estif(fvip, fui)    += v*gradp_(1,0)*derxy_(0, vi) ;
              estif(fvipp, fui)   += v*gradp_(2,0)*derxy_(0, vi) ;

              estif(fvi, fuip)    += v*gradp_(0,0)*derxy_(1, vi) ;
              estif(fvip, fuip)   += v*gradp_(1,0)*derxy_(1, vi) ;
              estif(fvipp, fuip)  += v*gradp_(2,0)*derxy_(1, vi) ;

              estif(fvi, fuipp)   += v*gradp_(0,0)*derxy_(2, vi) ;
              estif(fvip, fuipp)  += v*gradp_(1,0)*derxy_(2, vi) ;
              estif(fvipp, fuipp) += v*gradp_(2,0)*derxy_(2, vi) ;

            }
          }

          if (higher_order_ele)
          {
            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = -2.0*visceff*timetauM*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                const double v0 = v*visc_old_(0);
                const double v1 = v*visc_old_(1);
                const double v2 = v*visc_old_(2);

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                           \
                  |               / n+1 \    /          \     |
                  |  nabla o eps | u     |, | Du o nabla | v  |
                  |               \ (i) /    \          /     |
                  \                                           /
                */
                estif(fvi, fui)     += v0*derxy_(0, vi) ;
                estif(fvip, fui)    += v1*derxy_(0, vi) ;
                estif(fvipp, fui)   += v2*derxy_(0, vi) ;

                estif(fvi, fuip)    += v0*derxy_(1, vi) ;
                estif(fvip, fuip)   += v1*derxy_(1, vi) ;
                estif(fvipp, fuip)  += v2*derxy_(1, vi) ;

                estif(fvi, fuipp)   += v0*derxy_(2, vi) ;
                estif(fvip, fuipp)  += v1*derxy_(2, vi) ;
                estif(fvipp, fuipp) += v2*derxy_(2, vi) ;
              }
            }
          }

          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = -tau_M*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              const double v0 = v*rhsmom_(0);
              const double v1 = v*rhsmom_(1);
              const double v2 = v*rhsmom_(2);

              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                             \
                |              /          \     |
                |  rhsint   , | Du o nabla | v  |
                |              \          /     |
                \                             /

              */
              estif(fvi    , fui)   += (v0*derxy_(0, vi)) ;
              estif(fvip   , fui)   += (v1*derxy_(0, vi)) ;
              estif(fvipp  , fui)   += (v2*derxy_(0, vi)) ;

              estif(fvi    , fuip)  += (v0*derxy_(1, vi)) ;
              estif(fvip   , fuip)  += (v1*derxy_(1, vi)) ;
              estif(fvipp  , fuip)  += (v2*derxy_(1, vi)) ;

              estif(fvi    , fuipp) += (v0*derxy_(2, vi)) ;
              estif(fvip   , fuipp) += (v1*derxy_(2, vi)) ;
              estif(fvipp  , fuipp) += (v2*derxy_(2, vi)) ;

            } // vi
          } // ui
        } // if newton

#if 1
        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 4*vi;
          // supg stabilisation
          double v = -tau_M*conv_c_(vi);
          eforce(fvi)     += v*res_old_(0) ;
          eforce(fvi + 1) += v*res_old_(1) ;
          eforce(fvi + 2) += v*res_old_(2) ;
        }
#endif
      }

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART

      if (higher_order_ele)
      {
        if(vstab != Fluid3::viscous_stab_none)
        {
          const double two_visc_tauMp = vstabfac*2.0*visc*tau_Mp;
          // viscous stabilization either on left hand side or on right hand side
          if (vstab == Fluid3::viscous_stab_gls || vstab == Fluid3::viscous_stab_usfem)
          {
            const double two_visc_timetauMp   = vstabfac*2.0*visc*timetauMp;
            const double four_visc2_timetauMp = vstabfac*4.0*visceff*visc*timetauMp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<numnode; ++ui)
            {
              const double v = two_visc_tauMp*densamfunct_(ui)
#if 1
                         + two_visc_timetauMp*conv_c_(ui)
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
                estif(vi*4,     ui*4    ) += v*viscs2_(0, vi) ;
                estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, vi) ;
                estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, vi) ;

                estif(vi*4,     ui*4 + 1) += v*viscs2_(1, vi) ;
                estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, vi) ;
                estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, vi) ;

                estif(vi*4,     ui*4 + 2) += v*viscs2_(2, vi) ;
                estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, vi) ;
                estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, vi) ;
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
                estif(vi*4,     ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(1, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(2, vi)) ;
                estif(vi*4 + 1, ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(4, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(5, vi)) ;
                estif(vi*4 + 2, ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(2, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(5, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(8, vi)) ;

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
                estif(vi*4,     ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(2,ui)*viscs2_(2,vi)) ;
                estif(vi*4 + 1, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(4,vi)+viscs2_(2,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(2,vi)+viscs2_(1,ui)*viscs2_(5,vi)+viscs2_(2,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(4,ui)+viscs2_(2,vi)*viscs2_(5,ui)) ;
                estif(vi*4 + 1, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(4,ui)*viscs2_(4,vi)+viscs2_(5,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(2,vi)+viscs2_(4,ui)*viscs2_(5,vi)+viscs2_(5,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(2,ui)+viscs2_(1,vi)*viscs2_(5,ui)+viscs2_(2,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 1, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(1,vi)*viscs2_(2,ui)+viscs2_(4,vi)*viscs2_(5,ui)+viscs2_(5,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 2, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(2,ui)*viscs2_(2,vi)+viscs2_(5,ui)*viscs2_(5,vi)+viscs2_(8,ui)*viscs2_(8,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<numnode; ++ui)
              {
                double v = two_visc_timetauMp*densfunct_(ui,0);
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
                  estif(vi*4,     ui*4    ) += v*(viscs2_(0,vi)*vderxy_(0,0)+
                                                  viscs2_(1,vi)*vderxy_(1,0)+
                                                  viscs2_(2,vi)*vderxy_(2,0)) ;
                  estif(vi*4 + 1, ui*4    ) += v*(viscs2_(1,vi)*vderxy_(0,0)+
                                                  viscs2_(4,vi)*vderxy_(1,0)+
                                                  viscs2_(5,vi)*vderxy_(2,0)) ;
                  estif(vi*4 + 2, ui*4    ) += v*(viscs2_(2,vi)*vderxy_(0,0)+
                                                  viscs2_(5,vi)*vderxy_(1,0)+
                                                  viscs2_(8,vi)*vderxy_(2,0)) ;

                  estif(vi*4,     ui*4 + 1) += v*(viscs2_(0,vi)*vderxy_(0,1)+
                                                  viscs2_(1,vi)*vderxy_(1,1)+
                                                  viscs2_(2,vi)*vderxy_(2,1)) ;
                  estif(vi*4 + 1, ui*4 + 1) += v*(viscs2_(1,vi)*vderxy_(0,1)+
                                                  viscs2_(4,vi)*vderxy_(1,1)+
                                                  viscs2_(5,vi)*vderxy_(2,1)) ;
                  estif(vi*4 + 2, ui*4 + 1) += v*(viscs2_(2,vi)*vderxy_(0,1)+
                                                  viscs2_(5,vi)*vderxy_(1,1)+
                                                  viscs2_(8,vi)*vderxy_(2,1)) ;

                  estif(vi*4,     ui*4 + 2) += v*(viscs2_(0,vi)*vderxy_(0,2)+
                                                  viscs2_(1,vi)*vderxy_(1,2)+
                                                  viscs2_(2,vi)*vderxy_(2,2)) ;
                  estif(vi*4 + 1, ui*4 + 2) += v*(viscs2_(1,vi)*vderxy_(0,2)+
                                                  viscs2_(4,vi)*vderxy_(1,2)+
                                                  viscs2_(5,vi)*vderxy_(2,2)) ;
                  estif(vi*4 + 2, ui*4 + 2) += v*(viscs2_(2,vi)*vderxy_(0,2)+
                                                  viscs2_(5,vi)*vderxy_(1,2)+
                                                  viscs2_(8,vi)*vderxy_(2,2)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<numnode; ++vi)
          {

            /* viscous stabilisation */
            eforce(vi*4    ) -= two_visc_tauMp*(res_old_(0)*viscs2_(0, vi)+res_old_(1)*viscs2_(1, vi)+res_old_(2)*viscs2_(2, vi)) ;
            eforce(vi*4 + 1) -= two_visc_tauMp*(res_old_(0)*viscs2_(1, vi)+res_old_(1)*viscs2_(4, vi)+res_old_(2)*viscs2_(5, vi)) ;
            eforce(vi*4 + 2) -= two_visc_tauMp*(res_old_(0)*viscs2_(2, vi)+res_old_(1)*viscs2_(5, vi)+res_old_(2)*viscs2_(8, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab == Fluid3::continuity_stab_yes)
      {
        const double timetauC     = timefac*tau_C;
        const double rhs_tauC_div = rhsfac*tau_C*mdiv_/fac;

        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v0 = timetauC*densderxy_(0, ui);
          double v1 = timetauC*densderxy_(1, ui);
          double v2 = timetauC*densderxy_(2, ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* continuity stabilisation on left hand side */
            /*
              /                          \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
              \                          /
            */
            estif(fvi,  fui  ) += v0*densderxy_(0, vi) ;
            estif(fvip, fui  ) += v0*densderxy_(1, vi) ;
            estif(fvipp,fui  ) += v0*densderxy_(2, vi) ;

            estif(fvi,  fuip ) += v1*densderxy_(0, vi) ;
            estif(fvip, fuip ) += v1*densderxy_(1, vi) ;
            estif(fvipp,fuip ) += v1*densderxy_(2, vi) ;

            estif(fvi,  fuipp) += v2*densderxy_(0, vi) ;
            estif(fvip, fuipp) += v2*densderxy_(1, vi) ;
            estif(fvipp,fuipp) += v2*densderxy_(2, vi) ;
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* continuity stabilisation on right hand side */
          eforce(fvi  ) -= rhs_tauC_div*densderxy_(0, vi) ;
          eforce(fvip ) -= rhs_tauC_div*densderxy_(1, vi) ;
          eforce(fvipp) -= rhs_tauC_div*densderxy_(2, vi) ;
        }

        if (loma)
        {
          const double v = tau_C*rhscon_;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* continuity stabilisation of rhs term of continuity equation */
            eforce(fvi  ) += v*densderxy_(0, vi) ;
            eforce(fvip ) += v*densderxy_(1, vi) ;
            eforce(fvipp) += v*densderxy_(2, vi) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //     STABILIZATION, CROSS-STRESS PART (RESIDUAL-BASED VMM)

      if (cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
      {
        if (cross == Fluid3::cross_stress_stab)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            double v = timetauM*conv_resM_(ui,0);
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
              double v2 = v*funct_(vi,0);
              estif(vi*4    , ui*4    ) -= v2 ;
              estif(vi*4 + 1, ui*4 + 1) -= v2 ;
              estif(vi*4 + 2, ui*4 + 2) -= v2 ;
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
          double v = tau_M*funct_(vi,0);
          eforce(vi*4)     += v*(res_old_(0,0)*mderxy_(0,0) +
                                 res_old_(1,0)*mderxy_(0,1) +
                                 res_old_(2,0)*mderxy_(0,2));
          eforce(vi*4 + 1) += v*(res_old_(0,0)*mderxy_(1,0) +
                                 res_old_(1,0)*mderxy_(1,1) +
                                 res_old_(2,0)*mderxy_(1,2));
          eforce(vi*4 + 2) += v*(res_old_(0,0)*mderxy_(2,0) +
                                 res_old_(1,0)*mderxy_(2,1) +
                                 res_old_(2,0)*mderxy_(2,2));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if (reynolds == Fluid3::reynolds_stress_stab_only_rhs)
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
          eforce(vi*4    ) += v*res_old_(0);
          eforce(vi*4 + 1) += v*res_old_(1);
          eforce(vi*4 + 2) += v*res_old_(2);
        }
      } // end Reynolds-stress part on right hand side

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv != Fluid3::fssgv_no)
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
          eforce(vi*4    ) -= vartfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                      +    derxy_(1, vi)*fsvderxy_(0, 1)
                                      +    derxy_(1, vi)*fsvderxy_(1, 0)
                                      +    derxy_(2, vi)*fsvderxy_(0, 2)
                                      +    derxy_(2, vi)*fsvderxy_(2, 0)) ;
          eforce(vi*4 + 1) -= vartfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                      +    derxy_(0, vi)*fsvderxy_(1, 0)
                                      +2.0*derxy_(1, vi)*fsvderxy_(1, 1)
                                      +    derxy_(2, vi)*fsvderxy_(1, 2)
                                      +    derxy_(2, vi)*fsvderxy_(2, 1)) ;
          eforce(vi*4 + 2) -= vartfac*(    derxy_(0, vi)*fsvderxy_(0, 2)
                                      +    derxy_(0, vi)*fsvderxy_(2, 0)
                                      +    derxy_(1, vi)*fsvderxy_(1, 2)
                                      +    derxy_(1, vi)*fsvderxy_(2, 1)
                                      +2.0*derxy_(2, vi)*fsvderxy_(2, 2)) ;
        }
      }
    }

    // linearization with respect to mesh motion
    if (emesh.IsInitialized())
    {

      // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
      // xGderiv_ == xjm_

      // mass + rhs
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = fac*funct_(vi,0);
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*4    , ui*4    ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
          emesh(vi*4    , ui*4 + 1) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);
          emesh(vi*4    , ui*4 + 2) += v*(velint_(0)-rhsmom_(0))*derxy_(2, ui);

          emesh(vi*4 + 1, ui*4    ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
          emesh(vi*4 + 1, ui*4 + 1) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
          emesh(vi*4 + 1, ui*4 + 2) += v*(velint_(1)-rhsmom_(1))*derxy_(2, ui);

          emesh(vi*4 + 2, ui*4    ) += v*(velint_(2)-rhsmom_(2))*derxy_(0, ui);
          emesh(vi*4 + 2, ui*4 + 1) += v*(velint_(2)-rhsmom_(2))*derxy_(1, ui);
          emesh(vi*4 + 2, ui*4 + 2) += v*(velint_(2)-rhsmom_(2))*derxy_(2, ui);
        }
      }

      //vderiv_  = sum(evelnp(i,k) * deriv_(j,k), k);
      vderiv_.MultiplyNT(evelnp,deriv_);

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

      for (int ui=0; ui<numnode; ++ui)
      {
        double v00 = + convvelint_(1)*(vderiv_(0, 0)*derxjm_(0,0,1,ui) + vderiv_(0, 1)*derxjm_(0,1,1,ui) + vderiv_(0, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(0, 0)*derxjm_(0,0,2,ui) + vderiv_(0, 1)*derxjm_(0,1,2,ui) + vderiv_(0, 2)*derxjm_(0,2,2,ui));
        double v01 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(0, 0)*derxjm_(1,0,2,ui) + vderiv_(0, 1)*derxjm_(1,1,2,ui) + vderiv_(0, 2)*derxjm_(1,2,2,ui));
        double v02 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(0, 0)*derxjm_(2,0,1,ui) + vderiv_(0, 1)*derxjm_(2,1,1,ui) + vderiv_(0, 2)*derxjm_(2,2,1,ui));
        double v10 = + convvelint_(1)*(vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(1, 0)*derxjm_(0,0,2,ui) + vderiv_(1, 1)*derxjm_(0,1,2,ui) + vderiv_(1, 2)*derxjm_(0,2,2,ui));
        double v11 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(1,0,0,ui) + vderiv_(1, 1)*derxjm_(1,1,0,ui) + vderiv_(1, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(1, 0)*derxjm_(1,0,2,ui) + vderiv_(1, 1)*derxjm_(1,1,2,ui) + vderiv_(1, 2)*derxjm_(1,2,2,ui));
        double v12 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(2,0,0,ui) + vderiv_(1, 1)*derxjm_(2,1,0,ui) + vderiv_(1, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui));
        double v20 = + convvelint_(1)*(vderiv_(2, 0)*derxjm_(0,0,1,ui) + vderiv_(2, 1)*derxjm_(0,1,1,ui) + vderiv_(2, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui));
        double v21 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(1,0,0,ui) + vderiv_(2, 1)*derxjm_(1,1,0,ui) + vderiv_(2, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui));
        double v22 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(2,0,0,ui) + vderiv_(2, 1)*derxjm_(2,1,0,ui) + vderiv_(2, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(2, 0)*derxjm_(2,0,1,ui) + vderiv_(2, 1)*derxjm_(2,1,1,ui) + vderiv_(2, 2)*derxjm_(2,2,1,ui));

        for (int vi=0; vi<numnode; ++vi)
        {
          double v = timefacfac/det*funct_(vi);

          emesh(vi*4 + 0, ui*4 + 0) += v*v00;
          emesh(vi*4 + 0, ui*4 + 1) += v*v01;
          emesh(vi*4 + 0, ui*4 + 2) += v*v02;

          emesh(vi*4 + 1, ui*4 + 0) += v*v10;
          emesh(vi*4 + 1, ui*4 + 1) += v*v11;
          emesh(vi*4 + 1, ui*4 + 2) += v*v12;

          emesh(vi*4 + 2, ui*4 + 0) += v*v20;
          emesh(vi*4 + 2, ui*4 + 1) += v*v21;
          emesh(vi*4 + 2, ui*4 + 2) += v*v22;
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

      double v = visceff*timefac*fac;
      for (int ui=0; ui<numnode; ++ui)
      {
        double derinvJ0 = -v*(deriv_(0,ui)*xji_00 + deriv_(1,ui)*xji_01 + deriv_(2,ui)*xji_02);
        double derinvJ1 = -v*(deriv_(0,ui)*xji_10 + deriv_(1,ui)*xji_11 + deriv_(2,ui)*xji_12);
        double derinvJ2 = -v*(deriv_(0,ui)*xji_20 + deriv_(1,ui)*xji_21 + deriv_(2,ui)*xji_22);
        for (int vi=0; vi<numnode; ++vi)
        {
          double visres0 =   2.0*derxy_(0, vi)* vderxy_(0, 0)
                             +     derxy_(1, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                             +     derxy_(2, vi)*(vderxy_(0, 2) + vderxy_(2, 0)) ;
          double visres1 =         derxy_(0, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                             + 2.0*derxy_(1, vi)* vderxy_(1, 1)
                             +     derxy_(2, vi)*(vderxy_(1, 2) + vderxy_(2, 1)) ;
          double visres2 =         derxy_(0, vi)*(vderxy_(0, 2) + vderxy_(2, 0))
                             +     derxy_(1, vi)*(vderxy_(1, 2) + vderxy_(2, 1))
                             + 2.0*derxy_(2, vi)* vderxy_(2, 2) ;
          emesh(vi*4 + 0, ui*4 + 0) += derinvJ0*visres0;
          emesh(vi*4 + 1, ui*4 + 0) += derinvJ0*visres1;
          emesh(vi*4 + 2, ui*4 + 0) += derinvJ0*visres2;

          emesh(vi*4 + 0, ui*4 + 1) += derinvJ1*visres0;
          emesh(vi*4 + 1, ui*4 + 1) += derinvJ1*visres1;
          emesh(vi*4 + 2, ui*4 + 1) += derinvJ1*visres2;

          emesh(vi*4 + 0, ui*4 + 2) += derinvJ2*visres0;
          emesh(vi*4 + 1, ui*4 + 2) += derinvJ2*visres1;
          emesh(vi*4 + 2, ui*4 + 2) += derinvJ2*visres2;
        }
      }

      // part 2: derivative of viscosity residual

      v = timefacfac*visceff/det;
      for (int ui=0; ui<numnode; ++ui)
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
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

        for (int vi=0; vi<numnode; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }
      }


      // pressure
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = press*timefacfac/det;
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*4    , ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(0,0,1,ui) + deriv_(1, vi)*derxjm_(0,1,1,ui) + deriv_(2, vi)*derxjm_(0,2,1,ui)) ;
          emesh(vi*4    , ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(0,0,2,ui) + deriv_(1, vi)*derxjm_(0,1,2,ui) + deriv_(2, vi)*derxjm_(0,2,2,ui)) ;

          emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(1,0,0,ui) + deriv_(1, vi)*derxjm_(1,1,0,ui) + deriv_(2, vi)*derxjm_(1,2,0,ui)) ;
          emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(1,0,2,ui) + deriv_(1, vi)*derxjm_(1,1,2,ui) + deriv_(2, vi)*derxjm_(1,2,2,ui)) ;

          emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(2,0,0,ui) + deriv_(1, vi)*derxjm_(2,1,0,ui) + deriv_(2, vi)*derxjm_(2,2,0,ui)) ;
          emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(2,0,1,ui) + deriv_(1, vi)*derxjm_(2,1,1,ui) + deriv_(2, vi)*derxjm_(2,2,1,ui)) ;
        }
      }

      // div u
      for (int vi=0; vi<numnode; ++vi)
      {
        double v = timefacfac/det*funct_(vi,0);
        for (int ui=0; ui<numnode; ++ui)
        {
          emesh(vi*4 + 3, ui*4 + 0) += v*(
            + vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui)
            + vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui)
            ) ;

          emesh(vi*4 + 3, ui*4 + 1) += v*(
            + vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui)
            + vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui)
            ) ;

          emesh(vi*4 + 3, ui*4 + 2) += v*(
            + vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui)
            + vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui)
            ) ;
        }
      }

    }
  } // loop gausspoints
}



//
// calculate stabilization parameter
//
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Caltau(
  Fluid3*                                 ele,
  const LINALG::Matrix<3,iel>&            evelnp,
  const LINALG::Matrix<3,iel>&            fsevelnp,
  const LINALG::Matrix<iel,1>&            edensnp,
  const enum Fluid3::TauType              whichtau,
  struct _MATERIAL*                       material,
  double&                                 visc,
  const double                            timefac,
  const double                            dt,
  const enum Fluid3::TurbModelAction      turb_mod_action,
  double&                                 Cs,
  double&                                 Cs_delta_sq,
  double&                                 visceff,
  double&                                 l_tau,
  const enum Fluid3::StabilisationAction  fssgv
  )
{
  // use one-point Gauss rule to calculate tau at element center
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
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    integrationrule_stabili = DRT::UTILS::intrule_wedge_1point;
    break;
  case DRT::Element::pyramid5:
    integrationrule_stabili = DRT::UTILS::intrule_pyramid_1point;
    break;
  default:
    dserror("invalid discretization type for fluid3");
  }

  // Gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element-type constant for tau
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
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  velint_.Multiply(evelnp,funct_);

  // get density at element center
  const double dens = funct_.Dot(edensnp);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

  // compute element volume
  const double vol = wquad*det;

  // get characteristic element length for tau_Mp/tau_C:
  // volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  vderxy_.MultiplyNT(evelnp,derxy_);

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // normed velocity at element centre
  if (vel_norm>=1e-6)
  {
    velino_.Update(1.0/vel_norm,velint_);
  }
  else
  {
    velino_.Clear();
    velino_(0,0) = 1;
  }

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


  if (turb_mod_action == Fluid3::smagorinsky_with_wall_damping
      ||
      turb_mod_action == Fluid3::smagorinsky)
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
      double epsilon;
      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * vderxy_(rr,rr)*vderxy_(rr,rr);
      }
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
     //blitz::Array<double,1> centernodecoord(3);
      LINALG::Matrix<3,1> centernodecoord;
      //centernodecoord = blitz::sum(funct_(j)*xyze_(i,j),j);
      centernodecoord.Multiply(xyze_,funct_);

      if(centernodecoord(1,0)>0)
      {
        y_plus=(1.0-centernodecoord(1,0))/l_tau;
      }
      else
      {
        y_plus=(1.0+centernodecoord(1,0))/l_tau;
      }

//      lmix *= (1.0-exp(-y_plus/A_plus));
      // multiply with van Driest damping function
      Cs *= (1.0-exp(-y_plus/A_plus));
    }

    const double hk = pow((vol),(1.0/3.0));

    //
    // mixing length set proportional to grid witdh
    //
    //                     lmix = Cs * hk

    double lmix = Cs * hk;

    Cs_delta_sq = lmix * lmix;

    //
    //          visc    = visc + visc
    //              eff              turbulent

    visceff = visc + dens * Cs_delta_sq * rateofstrain;
  }
  else if(turb_mod_action == Fluid3::dynamic_smagorinsky)
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
    //                    mixing length
    //               provided by the dynamic model
    //            procedure and stored in Cs_delta_sq
    //

    double rateofstrain = 0;
    {
      double epsilon;
      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * vderxy_(rr,rr)*vderxy_(rr,rr);
      }
      rateofstrain = sqrt(rateofstrain);
    }

    visceff = visc + dens * Cs_delta_sq * rateofstrain;

    // for evaluation of statistics: remember the 'real' Cs
    Cs=sqrt(Cs_delta_sq)/pow((vol),(1.0/3.0));
  }
  else
  {
    visceff = visc;
  }

  // calculate tau

  if (whichtau == Fluid3::franca_barrenechea_valentin_wall)
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


    // get streamlength
    LINALG::Matrix<iel,1> tmp;
    tmp.MultiplyTN(derxy_,velino_);
    const double val = tmp.Norm1();
    const double strle = 2.0/val;

    /* viscous : reactive forces */
    const double re01 = 4.0 * timefac * visceff / (mk * dens * DSQR(strle));

    /* convective : viscous forces */
    const double re02 = mk * dens * vel_norm * strle / (2.0 * visceff);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*dens*xi01+(4.0*timefac*visceff/mk)*xi02);

    // compute tau_Mp
    //    stability parameter definition according to Franca and Valentin (2000)
    //                                       and Barrenechea and Valentin (2002)

    /* viscous : reactive forces */
    const double re11 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));

    /* convective : viscous forces */
    const double re12 = mk * dens * vel_norm * hk / (2.0 * visceff);

    const double xi11 = DMAX(re11,1.0);
    const double xi12 = DMAX(re12,1.0);

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
    tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*dens*xi11+(4.0*timefac*visceff/mk)*xi12);

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
    const double xi_tau_c = DMIN(re02,1.0);
    tau_(2) = vel_norm * hk * 0.5 * xi_tau_c / dens;
  }
  else if(whichtau == Fluid3::bazilevs)
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
    /*            +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
    */
    /*                      +----
           n+1       n+1     \     n+1          n+1
          u     * G u     =   +   u    * G   * u
                  -          /     i     -ij    j
                  -         +----        -
                             i,j
    */
    double G;
    double normG = 0;
    double Gnormu = 0;

    const double dens_sqr = dens*dens;
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0) + xji_(nn,1)*xji_(rr,1) + xji_(nn,2)*xji_(rr,2);
        normG+=G*G;
        Gnormu+=dens_sqr*velint_(nn,0)*G*velint_(rr,0);
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
    tau_(0) = 1.0/(sqrt((4.0*dens_sqr)/(dt*dt)+Gnormu+CI*visceff*visceff*normG));
    tau_(1) = tau_(0);

    /*           +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
    */
    /*           +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
    */
    double g;
    double normgsq = 0;
    for (int rr=0;rr<3;++rr)
    {
      g = xji_(rr,0) + xji_(rr,1) + xji_(rr,2);
      normgsq += g*g;
    }

    /*
                                1.0
                  tau  = -----------------
                     C            /     \
                          tau  * | g * g |
                             M    \-   -/
    */
    tau_(2) = 1./(tau_(0)*normgsq*dens_sqr);
  }
  else if(whichtau == Fluid3::codina)
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


    // get streamlength
    LINALG::Matrix<iel,1> tmp;
    tmp.MultiplyTN(derxy_,velino_);
    const double val = tmp.Norm1();
    const double strle = 2.0/val;

    /* viscous : reactive forces */
    const double re01 = 4.0 * timefac * visceff / (mk * dens * DSQR(strle));

    /* convective : viscous forces */
    const double re02 = mk * dens * vel_norm * strle / (2.0 * visceff);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*dens*xi01+(4.0*timefac*visceff/mk)*xi02);

    // compute tau_Mp
    //    stability parameter definition according to Franca and Valentin (2000)
    //                                       and Barrenechea and Valentin (2002)

    /* viscous : reactive forces */
    const double re11 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));

    /* convective : viscous forces */
    const double re12 = mk * dens * vel_norm * hk / (2.0 * visceff);

    const double xi11 = DMAX(re11,1.0);
    const double xi12 = DMAX(re12,1.0);

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
    tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*dens*xi11+(4.0*timefac*visceff/mk)*xi12);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff)+DSQR(0.5*dens*vel_norm*hk)) / ( dens*dens );
  }
  else
  {
    dserror("unknown definition of tau\n");
  }

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv == Fluid3::fssgv_artificial_all or
      fssgv == Fluid3::fssgv_artificial_small)
  {
    double fsvel_norm = 0.0;
    if (fssgv == Fluid3::fssgv_artificial_small)
    {
      // get fine-scale velocities at element center
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      fsvelint_.Multiply(fsevelnp,funct_);

      // get fine-scale velocity norm
      fsvel_norm = fsvelint_.Norm2();
    }
    // get all-scale velocity norm
    else fsvel_norm = vel_norm;

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = mk * dens * fsvel_norm * hk / visc; /* convective : viscous forces */
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(dens)*DSQR(fsvel_norm))/(2.0*visc*xi);

  }
  else if (fssgv == Fluid3::fssgv_Smagorinsky_all or
           fssgv == Fluid3::fssgv_Smagorinsky_small)
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
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      if (fssgv == Fluid3::fssgv_Smagorinsky_small)
           fsvderxy_.MultiplyNT(fsevelnp,derxy_);
      else fsvderxy_.MultiplyNT(evelnp,derxy_);

      double epsilon;
      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = fsvderxy_(rr,mm) + fsvderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * fsvderxy_(rr,rr)*fsvderxy_(rr,rr);
      }
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalVisc(
  const struct _MATERIAL* material,
  double&                 visc)
{

  // compute shear rate
  double rateofshear = 0.0;
  double tmp;
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<rr;mm++)
    {
      tmp = ( vderxy_(rr,mm) + vderxy_(mm,rr) );
      rateofshear += tmp*tmp;
    }
    rateofshear += 2.0 * vderxy_(rr,rr)*vderxy_(rr,rr);
  }
  rateofshear = sqrt(rateofshear);

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





/* If you make changes in this method please consider also changes in
   DRT::ELEMENTS::Fluid3lin_Impl::BodyForce() in fluid3_lin_impl.cpp */
/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::BodyForce(Fluid3* ele,
                                                   const double time)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

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
    const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the edeadng array
    for(int isd=0;isd<3;isd++)
    {
      // get factor given by spatial function
      if (functions) functnum = (*functions)[isd];
      else functnum = -1;

      double num = (*onoff)[isd]*(*val)[isd]*curvefac;

      for (int jnode=0; jnode<iel; jnode++)
      {
        if (functnum>0)
        {
          // evaluate function at the position of the current node
          functionfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X());
        }
        else functionfac = 1.0;

        edeadng_(isd,jnode) = num*functionfac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_.Clear();
  }
}


/* If you make changes in this method please consider also changes in
   DRT::ELEMENTS::Fluid3lin_Impl::gder2() in fluid3_lin_impl.cpp */
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::gder2(Fluid3* ele)
{
  // initialize and zero out everything
  static LINALG::Matrix<6,6> bm;

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

  //xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);
  xder2_.MultiplyNT(deriv2_,xyze_);

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

  //derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);
  derxy2_.Multiply(-1.0,xder2_,derxy_);

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

  //derxy2_ += deriv2_;
  derxy2_.Update(1.0,deriv2_,1.0);

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

  LINALG::FixedSizeSerialDenseSolver<6,6,iel> solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(derxy2_,derxy2_);
  solver.Solve();
}


#endif
#endif
