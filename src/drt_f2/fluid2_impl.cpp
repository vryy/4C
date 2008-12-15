/*----------------------------------------------------------------------*/
/*!
\file fluid2_impl.cpp

\brief Internal implementation of Fluid2 element

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
DRT::ELEMENTS::Fluid2ImplInterface* DRT::ELEMENTS::Fluid2ImplInterface::Impl(DRT::ELEMENTS::Fluid2* f2)
{
  switch (f2->Shape())
  {
  case DRT::Element::quad4:
  {
    static Fluid2Impl<DRT::Element::quad4>* fq4;
    if (fq4==NULL)
      fq4 = new Fluid2Impl<DRT::Element::quad4>;
    return fq4;
  }
  case DRT::Element::quad8:
  {
    static Fluid2Impl<DRT::Element::quad8>* fq8;
    if (fq8==NULL)
      fq8 = new Fluid2Impl<DRT::Element::quad8>;
    return fq8;
  }
  case DRT::Element::quad9:
  {
    static Fluid2Impl<DRT::Element::quad9>* fq9;
    if (fq9==NULL)
      fq9 = new Fluid2Impl<DRT::Element::quad9>;
    return fq9;
  }
  case DRT::Element::tri3:
  {
    static Fluid2Impl<DRT::Element::tri3>* ft3;
    if (ft3==NULL)
      ft3 = new Fluid2Impl<DRT::Element::tri3>;
    return ft3;
  }
  case DRT::Element::tri6:
  {
    static Fluid2Impl<DRT::Element::tri6>* ft6;
    if (ft6==NULL)
      ft6 = new Fluid2Impl<DRT::Element::tri6>;
    return ft6;
  }
  default:
    dserror("shape %d (%d nodes) not supported", f2->Shape(), f2->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid2Impl<distype>::Fluid2Impl()
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
    //velino_(2),
    velint_(),
    ndwvelint_(),
    fsvelint_(),
    convvelint_(),
    accintam_(),
    gradp_(),
    tau_(),
    viscs2_(),
    conv_c_(),
    ndwconv_c_(),
    mdiv_(),
    vdiv_(),
    rhsmom_(),
    rhscon_(),
    conv_old_(),
    visc_old_(),
    res_old_(),
    dwres_old_(),
    conv_resM_(),
    xder2_(),
    vderiv_()
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2Impl<distype>::Evaluate(
  Fluid2*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RefCountPtr<MAT::Material> mat,
  _MATERIAL*                 actmat)
{
  // the number of nodes
  const int numnode = iel;

  // construct views
  LINALG::Matrix<3*iel,3*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<3*iel,3*iel> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<3*iel,    1> elevec1(elevec1_epetra.A(),true);
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
  // get control parameters for linearization, low-Mach-number solver
  // and form of convective term
  //----------------------------------------------------------------------
  string newtonstr   = params.get<string>("Linearisation");
  string lomastr     = params.get<string>("low-Mach-number solver");
  string convformstr = params.get<string>("form of convective term");
  bool newton = false;
  bool loma   = false;
  bool conservative = false;
  if(newtonstr=="Newton")          newton       = true;
  if(lomastr  =="Yes")             loma         = true;
  if(convformstr =="conservative") conservative = true;

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid2::StabilisationAction pspg     = ele->ConvertStringToStabAction(stablist.get<string>("PSPG"));
  Fluid2::StabilisationAction supg     = ele->ConvertStringToStabAction(stablist.get<string>("SUPG"));
  Fluid2::StabilisationAction vstab    = ele->ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  Fluid2::StabilisationAction cstab    = ele->ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  Fluid2::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  Fluid2::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // select tau definition
  Fluid2::TauType whichtau = Fluid2::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");

    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid2::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid2::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid2::codina;
    }
  }

  // flag for higher order elements
  // this could be done better with XFEM::isHigherOrderElement, but
  // this is not the right place to include XFEM-stuff.
  bool higher_order_ele = ele->isHigherOrderElement(distype);

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
  LINALG::Matrix<2,numnode> evelnp;
  LINALG::Matrix<numnode,1> edensnp;
  LINALG::Matrix<2,numnode> emhist;
  LINALG::Matrix<numnode,1> echist;

  for (int i=0;i<numnode;++i)
  {
    // split velocity and pressure/density, insert into element arrays
    evelnp(0,i) = myvelnp[0+(i*3)];
    evelnp(1,i) = myvelnp[1+(i*3)];

    eprenp(i) = myvelnp[2+(i*3)];

    edensnp(i) = myvedenp[2+(i*3)];

    // momentum equation part of history vector
    // (containing information of time step t_n (mass rhs!))
    emhist(0,i) = myhist[0+(i*3)];
    emhist(1,i) = myhist[1+(i*3)];

    // continuity equation part (only non-trivial for low-Mach-number flow)
    echist(i) = myhist[2+(i*3)];
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
  LINALG::Matrix<2,numnode> eaccam;
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
      eaccam(0,i) = myaccam[0+(i*3)];
      eaccam(1,i) = myaccam[1+(i*3)];

      ededtam(i)  = myaccam[2+(i*3)];

      // extract density
      edensam(i)  = myvedeam[2+(i*3)];
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
      edensam(i,0) = myvedenp[2+(i*3)];
    }
  }

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  RCP<const Epetra_Vector> dispnp;
  vector<double> mydispnp;
  RCP<const Epetra_Vector> gridv;
  vector<double> mygridv;

  // create objects for element arrays
  LINALG::Matrix<2,numnode> edispnp;
  LINALG::Matrix<2,numnode> egridv;

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
      edispnp(0,i) = mydispnp[0+(i*3)];
      edispnp(1,i) = mydispnp[1+(i*3)];

      // set grid velocities
      egridv(0,i) = mygridv[0+(i*3)];
      egridv(1,i) = mygridv[1+(i*3)];
    }
  }

  // ---------------------------------------------------------------------
  // get additional state vector for AVM3 case: fine-scale velocity
  // values are at time n+alpha_F for generalized-alpha scheme and at
  // time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // get flag for fine-scale subgrid-viscosity approach
  Fluid2::StabilisationAction fssgv =
    ele->ConvertStringToStabAction(params.get<string>("fs subgrid viscosity","No"));

  RCP<const Epetra_Vector> fsvelnp;
  vector<double> myfsvelnp;

  // create object for element array
  LINALG::Matrix<2,numnode> fsevelnp;

  if (fssgv != Fluid2::fssgv_no)
  {
    fsvelnp = discretization.GetState("fsvelnp");
    if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
    myfsvelnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

    for (int i=0;i<numnode;++i)
    {
      // get fine-scale velocity
      fsevelnp(0,i) = myfsvelnp[0+(i*3)];
      fsevelnp(1,i) = myfsvelnp[1+(i*3)];
    }
  }

  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
  ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  // initialise the Smagorinsky constant Cs to zero
  double Cs            = 0.0;
  double visceff       = 0.0;

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv != Fluid2::fssgv_no && turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
    dserror("No combination of a classical (all-scale) turbulence model and a fine-scale subgrid-viscosity approach currently possible!");
  if (fssgv != Fluid2::fssgv_no) Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

  // the default action is no model
  Fluid2::TurbModelAction turb_mod_action = Fluid2::no_model;

  if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action = Fluid2::smagorinsky;
      Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    else
      dserror("For 2-D, up to now, only constant-coefficient Smagorinsky model is available");
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
         conservative,
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
         visceff);

  return 0;
}

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Impl<distype>::Sysmat(
  Fluid2*                                 ele,
  const LINALG::Matrix<2,iel>&            evelnp,
  const LINALG::Matrix<2,iel>&            fsevelnp,
  const LINALG::Matrix<iel,1>&            eprenp,
  const LINALG::Matrix<iel,1>&            edensnp,
  const LINALG::Matrix<2,iel>&            eaccam,
  const LINALG::Matrix<iel,1>&            ededtam,
  const LINALG::Matrix<iel,1>&            edensam,
  const LINALG::Matrix<2,iel>&            emhist,
  const LINALG::Matrix<iel,1>&            echist,
  const LINALG::Matrix<2,iel>&            edispnp,
  const LINALG::Matrix<2,iel>&            egridv,
  LINALG::Matrix<3*iel,3*iel>&            estif,
  LINALG::Matrix<3*iel,3*iel>&            emesh,
  LINALG::Matrix<3*iel,    1>&            eforce,
  struct _MATERIAL*                       material,
  double                                  time,
  double                                  dt,
  double                                  timefac,
  double                                  timefacrhs,
  const bool                              newton,
  const bool                              loma,
  const bool                              conservative,
  const bool                              is_genalpha,
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
  const int numnode = iel;

  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
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
    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);

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

    */

    // compute global first derivates
    derxy_.Multiply(xji_, deriv_);

    // density-weighted shape functions at n+1/n+alpha_F and n+1/n+alpha_M
    densfunct_.EMultiply(funct_,edensnp);
    densamfunct_.EMultiply(funct_,edensam);

    // inverse-density-weighted shape functions and density-weighted
    // global first derivatives at n+1/n+alpha_F
    for (int inode=0; inode<numnode; inode++)
    {
      functdens_(inode) = funct_(inode)/edensnp(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    //--------------------------------------------------------------
    //             compute global second derivatives
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
      // initialize everything
      LINALG::Matrix<3,3> bm;

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
      xder2_.MultiplyNT(deriv2_, xyze_);


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
      derxy2_.Multiply(-1.0, xder2_, derxy_);

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

      LINALG::FixedSizeSerialDenseSolver<3,3,iel> solver_bm;
      solver_bm.SetMatrix(bm);                // A = bm
      // X == B is not a problem, derxy2_ will contain the solution.
      solver_bm.SetVectors(derxy2_,derxy2_);  // X = B = derxy2_
      // error code
      int ierr = 0;

      // Factor. Calling this directly is not necessary, only to find
      // out where an error came from.
      ierr = solver_bm.Factor();
      if (ierr!=0)
        dserror("Unable to perform LU factorisation during computation of derxy2");

      solver_bm.Solve();                      // Solve A*X = B
      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");
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
    if (fssgv != Fluid2::fssgv_no) fsvderxy_.MultiplyNT(fsevelnp,derxy_);
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

    // For conservative form of convective term, we also need "pure"
    // (i.e., non-density-weighted) velocity (including potential ALE
    // velocity) and convective part.
    if (conservative)
    {
      ndwvelint_.Multiply(evelnp,funct_);
      if (ele->is_ale_) ndwvelint_.Multiply(-1.0, egridv, funct_, 1.0);
      ndwconv_c_.MultiplyTN(derxy_,ndwvelint_);
    }

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
           1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
           - |                              |         N_y .. y-line of N
           2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
             \                              /                            */

      /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
      /*    /                   \
          1 |  N_x,xx + N_y,yx  |
       -  - |                   |
          3 |  N_x,xy + N_y,yy  |
            \                   /

               with N_x .. x-line of N
               N_y .. y-line of N                                      */

      double prefac;
      if (loma)
      {
        prefac = 1.0/3.0;
        derxy2_.Scale(prefac);
      }
      else prefac = 1.0;

      double sum = (derxy2_(0,0)+derxy2_(1,0))/prefac;

      viscs2_(0,0) = 0.5 * (sum + derxy2_(0,0));
      viscs2_(1,0) = 0.5 * derxy2_(2,0);
      viscs2_(3,0) = 0.5 * (sum + derxy2_(1,0));

      /* viscous term  div epsilon(u_old) */
      visc_old_(0) = viscs2_(0,0)*evelnp(0,0)+viscs2_(1,0)*evelnp(1,0);
      visc_old_(1) = viscs2_(1,0)*evelnp(0,0)+viscs2_(3,0)*evelnp(1,0);

      for (int i=1; i<numnode; ++i)
      {
        double sum = (derxy2_(0,i)+derxy2_(1,i))/prefac;

        viscs2_(0,i) = 0.5 * (sum + derxy2_(0,i));
        viscs2_(1,i) = 0.5 * derxy2_(2,i);
        viscs2_(3,i) = 0.5 * (sum + derxy2_(1,i));

        /* viscous term  div epsilon(u_old) */
        visc_old_(0) += viscs2_(0,i)*evelnp(0,i)+viscs2_(1,i)*evelnp(1,i);
        visc_old_(1) += viscs2_(1,i)*evelnp(0,i)+viscs2_(3,i)*evelnp(1,i);
      }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }

    // momentum and velocity divergence
    // (the later is only required for low-Mach-number flow)
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

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
      for (int rr=0;rr<2;++rr)
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
      for (int rr=0;rr<2;++rr)
      {
        res_old_(rr) = velint_(rr)-rhsmom_(rr)+timefac*(conv_old_(rr)+gradp_(rr)-2*visceff*visc_old_(rr));
      }
    }

    //--------------------------------------------------------------------
    // calculation of additional terms when cross- and Reynolds-stress
    // are included
    //--------------------------------------------------------------------
    if (cross    == Fluid2::cross_stress_stab ||
        cross    == Fluid2::cross_stress_stab_only_rhs ||
        reynolds == Fluid2::reynolds_stress_stab_only_rhs)
    {
      // get density at element center
      const double dens = funct_.Dot(edensnp);

      // get density-weighted residual
      dwres_old_.Update(dens,res_old_,0.0);

      /*
        This is the operator

                  /                      \
                 | (rho*resM)    o nabla |
                  \         (i)          /

                  required for (lhs) cross- and (rhs) Reynolds-stress calculation

      */
      conv_resM_.MultiplyTN (derxy_,dwres_old_);
    }

    //--------------------------------------------------------------------
    // perform integration for element matrix and right hand side
    //--------------------------------------------------------------------
    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      // computation of inertia and convection (convective and reactive part)
      // for conservative form including right-hand-side contribution
      if (conservative)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = fac*densamfunct_(ui);

          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            /* inertia (contribution to mass matrix) */
            /*

            /               \
            |                |
            |  D(rho*u) , v  |
            |                |
            \               /
            */
            double v2 = v*funct_(vi) ;
            estif(tvi , tui ) += v2;
            estif(tvip, tuip) += v2;
          }
        }

#if 1
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui = 3*ui;
          const int tuip = tui+1;
          const double v0 = timefacfac*funct_(ui)*velint_(0);
          const double v1 = timefacfac*funct_(ui)*velint_(1);

          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi = 3*vi;
            const int tvip = tvi+1;
            /* convection, convective part */
            /*

            /                                   \
            |  /              n+1  \             |
            | | Du (x) (rho*u)     |  , nabla v  |
            |  \              (i)  /             |
            \                                   /

            */
            estif(tvi,  tui ) -= v0*derxy_(0, vi);
            estif(tvi,  tuip) -= v0*derxy_(1, vi);
            estif(tvip, tui ) -= v1*derxy_(0, vi);
            estif(tvip, tuip) -= v1*derxy_(1, vi);
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi  = 3*vi;
          const int tvip = tvi+1;
          /* convection */
          double v = rhsfac*ndwconv_c_(vi);
          eforce(tvi ) += v*velint_(0) ;
          eforce(tvip) += v*velint_(1) ;
        }
#endif

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui = 3*ui;
            const int tuip = tui+1;
            const double v = timefacfac*densfunct_(ui);

            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvi = 3*vi;
              const int tvip = tvi+1;
              /* convection, reactive part */
              /*

              /                                     \
              |  / n+1                \             |
              | | u    (x) D(rho*u)   |  , nabla v  |
              |  \ (i)               /              |
              \                                    /

              */
              double v2 = v*ndwconv_c_(vi) ;
              estif(tvi,  tui ) -= v2;
              estif(tvip, tuip) -= v2;
            }
          }
        }
      }
      // computation of inertia and convection (convective and reactive part)
      // for convective form including right-hand-side contribution
      else
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = fac*densamfunct_(ui)
#if 1
                           + timefacfac*conv_c_(ui)
#endif
                           ;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            /* inertia (contribution to mass matrix) */
            /*

            /                \
            |                |
            |  D(rho*u) , v  |
            |                |
            \                /
            */

            /* convection, convective part */
            /*

            /                               \
            |  /       n+1        \         |
            | | (rho*u)   o nabla | Du , v  |
            |  \      (i)        /          |
            \                              /

            */
            double v2 = v*funct_(vi) ;
            estif(tvi , tui ) += v2;
            estif(tvip, tuip) += v2;
          }
        }

        if (newton)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            const double v = timefacfac*funct_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v2 = v*densfunct_(ui);

              /*  convection, reactive part

              /                                 \
              |  /                \   n+1       |
              | | D(rho*u) o nabla | u     , v  |
              |  \                /   (i)       |
              \                                /
              */
              estif(tvi,  tui ) += v2*vderxy_(0, 0) ;
              estif(tvi,  tuip) += v2*vderxy_(0, 1) ;
              estif(tvip, tui ) += v2*vderxy_(1, 0) ;
              estif(tvip, tuip) += v2*vderxy_(1, 1) ;
            }
          }
        }

#if 1
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          const int tvip = tvi+1;
          /* convection */
          double v = -rhsfac*funct_(vi);
          eforce(tvi ) += v*conv_old_(0) ;
          eforce(tvip) += v*conv_old_(1) ;
        }
#endif
      }

      const double viscefftimefacfac = visceff*timefacfac;
      for (int ui=0; ui<numnode; ++ui)
      {
        const int tui  = 3*ui;
        const int tuip = tui+1;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi  = 3*vi;
          const int tvip = tvi+1;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
          /* viscosity term */
          /*

                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(tvi,  tui ) += viscefftimefacfac*(2.0*derxy_0ui_0vi
                                                  +
                                                  derxy_1ui_1vi) ;
          estif(tvi,  tuip) += viscefftimefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(tvip, tui ) += viscefftimefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(tvip, tuip) += viscefftimefacfac*(derxy_0ui_0vi
                                                  +
                                                  2.0*derxy_1ui_1vi) ;

        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        const int tuipp = 3*ui + 2;
        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* pressure */
          /*

          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */

          estif(tvi,     tuipp) += v*derxy_(0, vi) ;
          estif(tvi + 1, tuipp) += v*derxy_(1, vi) ;

        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvipp = 3*vi + 2;
        const double v = timefacfac*functdens_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui = 3*ui;
          /* divergence */
          /*
            /                              \
            |                              |
            | nabla o D(rho*u)  , (q/rho)  |
            |                              |
            \                              /
          */
          estif(tvipp, tui    ) += v*densderxy_(0, ui) ;
          estif(tvipp, tui + 1) += v*densderxy_(1, ui) ;
        }
      }

      if(is_genalpha)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* inertia */
          const double v = -fac*funct_(vi);
          eforce(tvi    ) += v*accintam_(0) ;
          eforce(tvi + 1) += v*accintam_(1) ;
        }
      }
      else
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* inertia */
          const double v = -fac*funct_(vi);
          eforce(tvi    ) += v*velint_(0) ;
          eforce(tvi + 1) += v*velint_(1) ;
        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvi = 3*vi;
        /* pressure */
        double v = press*rhsfac;
        eforce(tvi    ) += v*derxy_(0, vi) ;
        eforce(tvi + 1) += v*derxy_(1, vi) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvi = 3*vi;
        const double v = -visceff*rhsfac;
        /* viscosity */
        eforce(tvi    ) += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                              +
                              derxy_(1, vi)*vderxy_(0, 1)
                              +
                              derxy_(1, vi)*vderxy_(1, 0)) ;
        eforce(tvi + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                              +
                              derxy_(0, vi)*vderxy_(1, 0)
                              +
                              2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvi = 3*vi;
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eforce(tvi    ) += v*rhsmom_(0) ;
        eforce(tvi + 1) += v*rhsmom_(1) ;
      }

      const double rhsfac_mdiv = -rhsfac * mdiv_;
      for (int vi=0; vi<numnode; ++vi)
      {
        // continuity equation
        eforce(vi*3 + 2) += rhsfac_mdiv*functdens_(vi) ;
      }

      if (loma)
      {
        const double v = -(2.0/3.0)*visceff*timefacfac ;
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v0 = v*derxy_(0,ui);
          const double v1 = v*derxy_(1,ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;

            /* viscosity term - subtraction for low-Mach-number flow */
            /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
            */
            estif(tvi,  tui ) += v0*derxy_(0, vi) ;
            estif(tvi,  tuip) += v1*derxy_(0, vi) ;
            estif(tvip, tui ) += v0*derxy_(1, vi) ;
            estif(tvip, tuip) += v1*derxy_(1, vi) ;

          }
        }

        const double v_div = (2.0/3.0)*visceff*rhsfac*vdiv_ ;
        const double fac_rhscon = fac*rhscon_;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* viscosity term - subtraction for low-Mach-number flow */
          eforce(tvi    ) += derxy_(0, vi)*v_div ;
          eforce(tvi + 1) += derxy_(1, vi)*v_div ;

          /* additional rhs term of continuity equation */
          eforce(tvi + 2) += fac_rhscon*functdens_(vi) ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg == Fluid2::pstab_use_pspg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = tau_Mp*densamfunct_(ui)
#if 1
                           + timetauMp*conv_c_(ui)
#endif
                           ;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvipp  = 3*vi + 2;

            /* pressure stabilisation: inertia */
            /*
              /                      \
              |                      |
              |  D(rho*u) , nabla q  |
              |                      |
              \                     /
            */
            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | | (rho*u)   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */

            estif(tvipp, tui ) += v*derxy_(0, vi) ;
            estif(tvipp, tuip) += v*derxy_(1, vi) ;
          }
        }

        if (higher_order_ele)
        {
          const double v = -2.0*visceff*timetauMp;
          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvipp  = 3*vi + 2;

              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(tvipp, tui ) += v*(derxy_(0, vi)*viscs2_(0, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(1, ui)) ;
              estif(tvipp, tuip) += v*(derxy_(0, vi)*viscs2_(1, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(3, ui)) ;
            }
          }
        }

        for (int ui=0; ui<numnode; ++ui)
        {
          const int tuipp = 3*ui + 2;
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
            estif(vi*3 + 2, tuipp) += timetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = timetauMp*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvipp = 3*vi + 2;
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | | D(rho*u) o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

              */
              estif(tvipp, tui ) += v*(derxy_(0, vi)*vderxy_(0, 0)
                                       +
                                       derxy_(1, vi)*vderxy_(1, 0)) ;
              estif(tvipp, tuip) += v*(derxy_(0, vi)*vderxy_(0, 1)
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
#if 1
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = tau_M*densamfunct_(ui) + timetauM*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi = 3*vi;
            /* supg stabilisation: inertia  */
            /*
              /                                      \
              |             /        n+1       \     |
              |  D(rho*u) , | (rho*u)   o nabla | v  |
              |              \       (i)       /     |
              \                                      /
            */

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                         \
            |    /       n+1        \        /       n+1        \     |
            |   | (rho*u)    o nabla | Du , | (rho*u)    o nabla | v  |
            |    \       (i)        /        \       (i)        /     |
            \                                                         /

            */

            estif(tvi,     tui ) += v*conv_c_(vi);
            estif(tvi + 1, tuip) += v*conv_c_(vi);
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi  = 3*vi;
          const int tvip = tvi+1;
          const double v = timetauM*conv_c_(vi);
          for (int ui=0; ui<numnode; ++ui)
          {
            const int tuipp = 3*ui + 2;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , | (rho*u)   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(tvi,  tuipp) += v*derxy_(0, ui) ;
            estif(tvip, tuipp) += v*derxy_(1, ui) ;
          }
        }

        if (higher_order_ele)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            const double v = -2.0*visceff*timetauM*conv_c_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, | (rho*u)    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
              */
              estif(tvi,  tui ) += v*viscs2_(0, ui) ;
              estif(tvip, tui ) += v*viscs2_(1, ui) ;

              estif(tvi,  tuip) += v*viscs2_(1, ui) ;
              estif(tvip, tuip) += v*viscs2_(3, ui) ;
            }
          }
        }
#endif

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = tau_M*densamfunct_(ui);
            const double v0 = v*velint_(0);
            const double v1 = v*velint_(1);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;
              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                /                                         \
                |         n+1      /                \     |
                |  (rho*u)      , | D(rho*u) o nabla | v  |
                |         (i)      \                /     |
                \                                         /

              */
              estif(tvi,  tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi,  tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;
            }
          }

#if 1
          {
            const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1);
            const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1);

            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = timetauM*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;

                /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                /*
                  /                                                         \
                  |    /       n+1        \   n+1    /                \     |
                  |   | (rho*u)    o nabla | u    , | D(rho*u) o nabla | v  |
                  |    \       (i)        /   (i)    \                /     |
                  \                                                        /

                  /                                                         \
                  |    /                \   n+1    /       n+1        \     |
                  |   | D(rho*u) o nabla | u    , | (rho*u)    o nabla | v  |
                  |    \                /   (i)    \       (i)        /     |
                  \                                                        /
                */
                estif(tvi,  tui ) += (conv_c_(vi)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(tvip, tui ) += (conv_c_(vi)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;

                estif(tvi,  tuip) += (conv_c_(vi)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(tvip, tuip) += (conv_c_(vi)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
              }
            }
          }
#endif

          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = timetauM*densfunct_(ui);
            const double v0 = v*gradp_(0);
            const double v1 = v*gradp_(1);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , | D(rho*u) o nabla | v  |
                |         (i)    \                /     |
                \                                      /
              */
              estif(tvi,  tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi,  tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;

            }
          }

          if (higher_order_ele)
          {
            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = -2.0*visceff*timetauM*densfunct_(ui);
              const double v0 = v*visc_old_(0);
              const double v1 = v*visc_old_(1);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, | D(rho*u) o nabla | v  |
                  |               \ (i) /    \                /     |
                  \                                                 /
                */
                estif(tvi,  tui ) += v0*derxy_(0, vi) ;
                estif(tvip, tui ) += v1*derxy_(0, vi) ;

                estif(tvi,  tuip) += v0*derxy_(1, vi) ;
                estif(tvip, tuip) += v1*derxy_(1, vi) ;

              }
            }
          }

          for (int ui=0; ui<numnode; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = -tau_M*densfunct_(ui);
            const double v0 = v*rhsmom_(0);
            const double v1 = v*rhsmom_(1);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;

              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                                     \
                |              /                \     |
                |  rhsint   , | D(rho*u) o nabla | v  |
                |              \                /     |
                \                                     /

              */
              estif(tvi , tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi , tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;

            } // vi
          } // ui
        } // if newton

#if 1
        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          // supg stabilisation
          const double v = -tau_M*conv_c_(vi);
          eforce(tvi    ) += v*res_old_(0) ;
          eforce(tvi + 1) += v*res_old_(1) ;
        }
#endif
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
            const double two_visc_timetauMp   = vstabfac*2.0*visc*timetauMp;
            const double four_visc2_timetauMp = vstabfac*4.0*visceff*visc*timetauMp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = two_visc_tauMp*densamfunct_(ui)
#if 1
                         + two_visc_timetauMp*conv_c_(ui)
#endif
                               ;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;
                /* viscous stabilisation, inertia part */
                /*
                  /                          \
                  |                          |
              +/- |  D(rho*u) , div eps (v)  |
                  |                          |
                  \                          /
                */
                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | | (rho*u)   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
                */
                estif(tvi,  tui ) += v*viscs2_(0, vi) ;
                estif(tvip, tui ) += v*viscs2_(1, vi) ;

                estif(tvi,  tuip) += v*viscs2_(1, vi) ;
                estif(tvip, tuip) += v*viscs2_(3, vi) ;
              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              const int tuipp = 3*ui + 2;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int tvi = 3*vi;
                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(tvi,     tuipp) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                            +
                                                            derxy_(1, ui)*viscs2_(1, vi)) ;
                estif(tvi + 1, tuipp) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                            +
                                                            derxy_(1, ui)*viscs2_(3, vi)) ;

              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(tvi,  tui ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)) ;
                estif(tvip, tui ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(3,vi)) ;

                estif(tvi,  tuip) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(3,ui)) ;
                estif(tvip, tuip) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(3,ui)*viscs2_(3,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<numnode; ++ui)
              {
                const int tui  = 3*ui;
                const int tuip = tui+1;
                const double v = two_visc_timetauMp*densfunct_(ui);
                for (int vi=0; vi<numnode; ++vi)
                {
                  const int tvi  = 3*vi;
                  const int tvip = tvi+1;
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | | D(rho*u) o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
                  */
                  estif(tvi,  tui ) += v*(viscs2_(0,vi)*vderxy_(0,0)+viscs2_(1,vi)*vderxy_(1,0)) ;
                  estif(tvip, tui ) += v*(viscs2_(1,vi)*vderxy_(0,0)+viscs2_(3,vi)*vderxy_(1,0)) ;

                  estif(tvi,  tuip) += v*(viscs2_(0,vi)*vderxy_(0,1)+viscs2_(1,vi)*vderxy_(1,1)) ;
                  estif(tvip, tuip) += v*(viscs2_(1,vi)*vderxy_(0,1)+viscs2_(3,vi)*vderxy_(1,1)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi = 3*vi;
            /* viscous stabilisation */
            eforce(tvi    ) -= two_visc_tauMp*(res_old_(0)*viscs2_(0, vi)+res_old_(1)*viscs2_(1, vi)) ;
            eforce(tvi + 1) -= two_visc_tauMp*(res_old_(0)*viscs2_(1, vi)+res_old_(1)*viscs2_(3, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab == Fluid2::continuity_stab_yes)
      {
        const double timetauC     = timefac*tau_C;
        const double rhs_tauC_div = rhsfac*tau_C*mdiv_/fac;

        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v0 = timetauC*densderxy_(0, ui);
          const double v1 = timetauC*densderxy_(1, ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            /* continuity stabilisation on left hand side */
            /*
              /                                      \
              |                                      |
              | nabla o D(rho*u)  , nabla o (rho*v)  |
              |                                      |
              \                                      /
            */
            estif(tvi,  tui ) += v0*densderxy_(0, vi) ;
            estif(tvip, tui ) += v0*densderxy_(1, vi) ;

            estif(tvi,  tuip) += v1*densderxy_(0, vi) ;
            estif(tvip, tuip) += v1*densderxy_(1, vi) ;
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* continuity stabilisation on right hand side */
          eforce(tvi    ) -= rhs_tauC_div*densderxy_(0, vi) ;
          eforce(tvi + 1) -= rhs_tauC_div*densderxy_(1, vi) ;
        }

        if (loma)
        {
          const double v = tau_C*rhscon_;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int tvi = 3*vi;
            /* continuity stabilisation of rhs term of continuity equation */
            eforce(tvi    ) += v*densderxy_(0, vi) ;
            eforce(tvi + 1) += v*densderxy_(1, vi) ;
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
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = timetauM*conv_resM_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int tvi = 3*vi;
              /* cross-stress part on lhs */
              /*

                          /                               \
                         |  /                  \          |
                      -  | | (rho*resM) o nabla | Du , v  |
                         |  \                  /          |
                          \                               /
              */
              double v2 = v*funct_(vi);
              estif(tvi    , tui ) -= v2 ;
              estif(tvi + 1, tuip) -= v2 ;
            }
          }
        } // end cross-stress part on left hand side

        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* cross-stress part on rhs */
          /*

                          /                                \
                         |  /                  \           |
                         | | (rho*resM) o nabla | u   , v  |
                         |  \                  /  (i)      |
                          \                               /
          */
          const double v = tau_M*funct_(vi);
          eforce(tvi    ) += v*(dwres_old_(0)*vderxy_(0,0)+dwres_old_(1)*vderxy_(0,1));
          eforce(tvi + 1) += v*(dwres_old_(0)*vderxy_(1,0)+dwres_old_(1)*vderxy_(1,1));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if (reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      {
        const double tauMtauM = tau_M*tau_M/fac;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* Reynolds-stress part on rhs */
          /*

                  /                                    \
                 |                                     |
                 |  resM   ,   (rho*resM) o nabla ) v  |
                 |                                     |
                  \                                    /
          */
          double v = tauMtauM*conv_resM_(vi);
          eforce(tvi    ) += v*res_old_(0);
          eforce(tvi + 1) += v*res_old_(1);
        }
      } // end Reynolds-stress part on right hand side

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv != Fluid2::fssgv_no)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int tvi = 3*vi;
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(tvi    ) -= vartfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                      +    derxy_(1, vi)*fsvderxy_(0, 1)
                                      +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
          eforce(tvi + 1) -= vartfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                      +    derxy_(0, vi)*fsvderxy_(1, 0)
                                      +2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
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
        const int tvi   = 3*vi;
        const int tvip  = tvi + 1;

        const double v = fac*funct_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui   = 3*ui;
          const int tuip  = tui + 1;

          emesh(tvi,   tui ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
          emesh(tvi,   tuip) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);

          emesh(tvip,  tui ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
          emesh(tvip,  tuip) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
        }
      }

      vderiv_.MultiplyNT(evelnp, deriv_);

//#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

//#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
//#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
//#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
//#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))

      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvi  = 3*vi;
        const int tvip = tvi+1;
        const double v = timefacfac/det*funct_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;

          emesh(tvi , tui ) += v*(
          + convvelint_(1)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(tvi , tuip) += v*(
          + convvelint_(0)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(tvip, tui ) += v*(
          + convvelint_(1)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );

          emesh(tvip, tuip) += v*(
          + convvelint_(0)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );
        }
      }

      // pressure
      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvi  = 3*vi;
        const int tvip = tvi+1;
        const double v = press*timefacfac/det;
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui = 3*ui;
          emesh(tvi,  tui + 1) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
          emesh(tvip, tui    ) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
        }
      }

      // div u
      for (int vi=0; vi<numnode; ++vi)
      {
        const int tvipp = 3*vi + 2;
        const double v = timefacfac/det*functdens_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          const int tui = 3*ui;
          emesh(tvipp, tui) += v*(
          deriv_(0,ui)*vderiv_(1,1) - deriv_(1,ui)*vderiv_(1,0)
          ) ;

          emesh(tvipp, tui + 1) += v*(
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Impl<distype>::Caltau(
  Fluid2*                                 ele,
  const LINALG::Matrix<2,iel>&            evelnp,
  const LINALG::Matrix<2,iel>&            fsevelnp,
  const LINALG::Matrix<iel,1>&            edensnp,
  const enum Fluid2::TauType              whichtau,
  struct _MATERIAL*                       material,
  double&                                 visc,
  const double                            timefac,
  const double                            dt,
  const enum Fluid2::TurbModelAction      turb_mod_action,
  double&                                 Cs,
  double&                                 visceff,
  const enum Fluid2::StabilisationAction  fssgv
  )
{
  // use one-point Gauss rule to calculate tau at element center
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

  // Gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];

  DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
  DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);

  // get element-type constant for tau
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
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  velint_.Multiply(evelnp,funct_);

  // get density at element center
  const double dens = funct_.Dot(edensnp);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xji_.Invert(xjm_);

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

  */

  // compute global derivates
  derxy_.Multiply(xji_,deriv_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  vderxy_.MultiplyNT(evelnp,derxy_);

  // get velocity norm
  const double vel_norm = velint_.Norm2();

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
      double epsilon;
      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * vderxy_(rr,rr) * vderxy_(rr,rr);
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

    tau_(0) = timefac*DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
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
    tau_(2) = vel_norm * hk * 0.5 * xi_tau_c / dens;
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
    double g;
    double normG = 0;
    double Gnormu = 0;

    const double dens_sqr = dens*dens;
    for (int nn=0;nn<2;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*velint_(nn);
      for (int rr=0;rr<2;++rr)
      {
        g = xji_(nn,0)*xji_(rr,0);
        for (int mm=1;mm<2;++mm)
        {
          g += xji_(nn,mm)*xji_(rr,mm);
        }
        normG += g*g;
        Gnormu+=dens_sqr_velint_nn*g*velint_(rr);
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
    const double g0 = xji_(0,0) + xji_(0,1);
    const double g1 = xji_(1,0) + xji_(1,1);

    /*           +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
    */
    const double normgsq = g0*g0+g1*g1;

    /*
                                1.0
                  tau  = -----------------
                     C            /     \
                          tau  * | g * g |
                             M    \-   -/
    */
    tau_(2) = 1./(tau_(0)*normgsq*dens_sqr);
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

    tau_(0) = timefac*DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
    tau_(1) = tau_(0);

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
  if (fssgv == Fluid2::fssgv_artificial_all || fssgv == Fluid2::fssgv_artificial_small)
  {
    double fsvel_norm = 0.0;
    if (fssgv == Fluid2::fssgv_artificial_small)
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
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      if (fssgv == Fluid2::fssgv_Smagorinsky_small)
           fsvderxy_.MultiplyNT(fsevelnp,derxy_);
      else fsvderxy_.MultiplyNT(evelnp,derxy_);

      double epsilon;
      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = fsvderxy_(rr,mm) + fsvderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * fsvderxy_(rr,rr) * fsvderxy_(rr,rr);
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
void DRT::ELEMENTS::Fluid2Impl<distype>::CalVisc(
  const struct _MATERIAL*                 material,
  double&                                 visc)
{

  // compute shear rate
  double rateofshear = 0.0;

  double epsilon;
  for(int rr=0;rr<2;rr++) {
    for(int mm=0;mm<rr;mm++) {
      epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
      rateofshear += epsilon*epsilon;
    }
    rateofshear += 2.0 * vderxy_(rr,rr) * vderxy_(rr,rr);
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



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Impl<distype>::BodyForce(Fluid2*      ele,
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
    for (int jnode=0; jnode<iel; jnode++)
    {
      const double* x = (ele->Nodes()[jnode])->X();
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
          functionfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(isd,x);
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
    edeadng_.Clear();
  }
}



#endif
#endif
