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
#include "fluid3_impl_sysmat2D.H"
#include "fluid3_impl_sysmat2D3D.H"

#include "fluid3_stabilization.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Epetra_SerialDenseSolver.h>

#ifdef DEBUG
#endif

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplInterface* DRT::ELEMENTS::Fluid3ImplInterface::Impl(DRT::ELEMENTS::Fluid3* f3)
{
  const int numdofpernode = f3->NumDofPerNode(*(f3->Nodes()[0]));
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Fluid3Impl<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new Fluid3Impl<DRT::Element::hex8>(numdofpernode);
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static Fluid3Impl<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new Fluid3Impl<DRT::Element::hex20>(numdofpernode);
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static Fluid3Impl<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new Fluid3Impl<DRT::Element::hex27>(numdofpernode);
    return fh27;
  }
  case DRT::Element::tet4:
  {
    static Fluid3Impl<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new Fluid3Impl<DRT::Element::tet4>(numdofpernode);
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static Fluid3Impl<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new Fluid3Impl<DRT::Element::tet10>(numdofpernode);
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static Fluid3Impl<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new Fluid3Impl<DRT::Element::wedge6>(numdofpernode);
    return fw6;
  }
  /* wedge15 cannot be used since no mesh generator exists
  case DRT::Element::wedge15:
  {
    static Fluid3Impl<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new Fluid3Impl<DRT::Element::wedge15>(numdofpernode);
    return fw15;
  }
  */
  case DRT::Element::pyramid5:
  {
    static Fluid3Impl<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new Fluid3Impl<DRT::Element::pyramid5>(numdofpernode);
    return fp5;
  }
  case DRT::Element::quad4:
  {
    static Fluid3Impl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new Fluid3Impl<DRT::Element::quad4>(numdofpernode);
    return cp4;
  }
  case DRT::Element::quad8:
  {
    static Fluid3Impl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new Fluid3Impl<DRT::Element::quad8>(numdofpernode);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static Fluid3Impl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new Fluid3Impl<DRT::Element::quad9>(numdofpernode);
    return cp9;
  }
  case DRT::Element::tri3:
  {
    static Fluid3Impl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new Fluid3Impl<DRT::Element::tri3>(numdofpernode);
    return cp3;
  }
  case DRT::Element::tri6:
  {
    static Fluid3Impl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new Fluid3Impl<DRT::Element::tri6>(numdofpernode);
    return cp6;
  }
  // no 1D elements
  /*
  case DRT::Element::line2:
  {
    static Fluid3Impl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new Fluid3Impl<DRT::Element::line2>(numdofpernode);
    return cl2;
  }
  case DRT::Element::line3:
  {
    static Fluid3Impl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new Fluid3Impl<DRT::Element::line3>(numdofpernode);
    return cl3;
  }
  case DRT::Element::line4:
  {
    static Fluid3Impl<DRT::Element::line4>* cl4;
    if (cl4==NULL)
      cl4 = new Fluid3Impl<DRT::Element::line4>(numdofpernode);
    return cl4;
  }
  case DRT::Element::line5:
  {
    static Fluid3Impl<DRT::Element::line5>* cl5;
    if (cl5==NULL)
      cl5 = new Fluid3Impl<DRT::Element::line5>(numdofpernode);
    return cl5;
  }
  //nurbs are not integrated yet
  case DRT::Element::nurbs2:
  {
    static Fluid3Impl<DRT::Element::nurbs2>* nu2;
    if (nu2==NULL)
      nu2 = new Fluid3Impl<DRT::Element::nurbs2>(numdofpernode);
    return nu2;
  }
  case DRT::Element::nurbs3:
  {
    static Fluid3Impl<DRT::Element::nurbs3>* nu3;
    if (nu3==NULL)
      nu3 = new Fluid3Impl<DRT::Element::nurbs3>(numdofpernode);
    return nu3;
  }
  case DRT::Element::nurbs4:
  {
    static Fluid3Impl<DRT::Element::nurbs4>* nu4;
    if (nu4==NULL)
      nu4 = new Fluid3Impl<DRT::Element::nurbs4>(numdofpernode);
    return nu4;
  }
  case DRT::Element::nurbs9:
  {
    static Fluid3Impl<DRT::Element::nurbs9>* nu9;
    if (nu9==NULL)
      nu9 = new Fluid3Impl<DRT::Element::nurbs9>(numdofpernode);
    return nu9;
  }
  case DRT::Element::nurbs8:
  {
    static Fluid3Impl<DRT::Element::nurbs8>* nu8;
    if (nu8==NULL)
      nu8 = new Fluid3Impl<DRT::Element::nurbs8>(numdofpernode);
    return nu8;
  }
  case DRT::Element::nurbs27:
  {
    static Fluid3Impl<DRT::Element::nurbs27>* nu27;
    if (nu27==NULL)
      nu27 = new Fluid3Impl<DRT::Element::nurbs27>(numdofpernode);
    return nu27;
  }
  */
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(f3->Shape()).c_str());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3Impl<distype>::Fluid3Impl(int numdofpernode)
  : numdofpernode_(numdofpernode),
    xyze_(true),
    edeadaf_(true),
    funct_(true),
    deriv_(true),
    deriv2_(true),
    xjm_(true),
    xji_(true),
    vderxy_(true),
    fsvderxy_(true),
    derxy_(true),
    derxy2_(true),
    bodyforce_(true),
    histmom_(true),
    velino_(true),
    velint_(true),
    fsvelint_(true),
    sgvelint_(true),
    convvelint_(true),
    accint_(true),
    gradp_(true),
    tau_(true),
    viscs2_(true),
    conv_c_(true),
    sgconv_c_(true),  // initialize to zero
    vdiv_(0.0),
    rhsmom_(true),
    conv_old_(true),
    visc_old_(true),
    momres_old_(true),  // initialize to zero
    conres_old_(true),
    xder2_(true),
    vderiv_(true),
    mat_gp_(false),     // standard evaluation of the material at the element center
    tau_gp_(false),     // standard evaluation of tau at the element center
    det_(0.0),
    fac_(0.0),
    visc_(0.0),
    sgvisc_(0.0),
    visceff_(0.0),
    fssgvisc_(0.0),
    dt_(0.0),
    omtheta_(0.0),
    gamma_(0.0),
    alphaF_(0.0),
    alphaM_(0.0),
    afgdt_(0.0),
    rhscon_(true),
    densaf_(1.0),         // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    densam_(1.0),         // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    densn_(1.0),          // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    scadtfac_(0.0),       // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    scaconvfacaf_(0.0),   // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    scaconvfacn_(0.0),    // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    thermpressadd_(0.0),  // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    vderxyn_(true),
    grad_scaaf_(true),
    grad_scan_(true),
    conv_scaaf_(0.0),
    conv_scan_(0.0),
    //thermpressaf_(),
    //thermpressam_(),
    //thermpressdtam_(),
    rotsymmpbc_(NULL),
    //flags
    is_genalpha_(false),
    is_conservative_(false),
    is_stationary_(false),
    is_newton_(false),
    is_higher_order_ele_(false),
    physicaltype_(INPAR::FLUID::incompressible),
    pspg_(INPAR::FLUID::pstab_use_pspg),
    supg_(INPAR::FLUID::convective_stab_supg),
    vstab_(INPAR::FLUID::viscous_stab_none),
    cstab_(INPAR::FLUID::continuity_stab_yes),
    cross_(INPAR::FLUID::cross_stress_stab_none),
    reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
    whichtau_(INPAR::FLUID::tau_not_defined),
    fssgv_(Fluid3::no_fssgv),
    turb_mod_action_(Fluid3::no_model)
{
  rotsymmpbc_= new FLD::RotationallySymmetricPeriodicBC<distype>();
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
  RefCountPtr<MAT::Material> mat)
{
  // construct views
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat1(elemat1_epetra,true);
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat2(elemat2_epetra,true);
  LINALG::Matrix<(nsd_+1)*nen_,    1> elevec1(elevec1_epetra,true);
  // elevec2 and elevec3 are currently not use

  // rotationally symmetric periodic bc's: do setup for current element
  rotsymmpbc_->Setup(ele);

  // set all solution parameters and flag, which does not change anymore at element level
  const double timefac = SetSolutionParameter(ele, params);

  // set thermodynamic pressure at n+1/n+alpha_F and n+alpha_M/n and
  // its time derivative at n+alpha_M/n+1
  const double thermpressaf   = params.get<double>("thermpress at n+alpha_F/n+1");
  const double thermpressam   = params.get<double>("thermpress at n+alpha_M/n");
  const double thermpressdtam = params.get<double>("thermpressderiv at n+alpha_M/n+1");

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, scalar,
  // acceleration/scalar time derivative and history
  // velocity/pressure and scalar values are at time n+alpha_F/n+alpha_M
  // for generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration/scalar time derivative values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
/*
  //RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("velaf");
  RefCountPtr<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  RefCountPtr<const Epetra_Vector> accam = discretization.GetState("accam");
  RefCountPtr<const Epetra_Vector> scaam = discretization.GetState("scaam");
  RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("hist");
  //if (velaf==null || accam==null || scaaf==null || scaam==null || hist==null)
  //  dserror("Cannot get state vectors 'velaf', 'accam', 'scaaf', 'scaam' and/or 'hist'");

  // extract local values from the global vectors
  //vector<double> myvelaf(lm.size());
  //DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
  vector<double> myaccam(lm.size());
  DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);
  vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);
  vector<double> myscaam(lm.size());
  DRT::UTILS::ExtractMyValues(*scaam,myscaam,lm);
  vector<double> myhist(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // rotation of velocity components
  // for rotationally symmetric boundary conditions
  //rotsymmpbc_->RotateMyValuesIfNecessary(myvelaf);
  //rotsymmpbc_->RotateMyValuesIfNecessary(myhist);
*/

  // fill the local element vector/matrix with the global values
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1> epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, &evelaf, &epreaf,"velaf");

  LINALG::Matrix<nen_,1> escaaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, NULL, &escaaf,"scaaf");

  LINALG::Matrix<nsd_,nen_> emhist(true);
  ExtractValuesFromGlobalVector(discretization,lm, &emhist, NULL,"hist");

  LINALG::Matrix<nsd_,nen_> eaccam(true);
  LINALG::Matrix<nen_,1> escadtam(true);
  ExtractValuesFromGlobalVector(discretization,lm, &eaccam, &escadtam,"accam");

  LINALG::Matrix<nsd_,nen_> eveln(true);
  LINALG::Matrix<nen_,1> escaam(true);
  ExtractValuesFromGlobalVector(discretization,lm, &eveln, &escaam,"scaam");

  if (is_genalpha_)
    eveln.Clear();
  else
    eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------

  LINALG::Matrix<nsd_, nen_> edispnp(true);
  LINALG::Matrix<nsd_, nen_> egridv(true);

  if(ele-> IsAle())
  {
    ExtractValuesFromGlobalVector(discretization,lm, &edispnp, NULL,"dispnp");
    ExtractValuesFromGlobalVector(discretization,lm, &egridv, NULL,"gridv");
  }

  // ---------------------------------------------------------------------
  // get additional state vector for AVM3 case: fine-scale velocity
  // values are at time n+alpha_F for generalized-alpha scheme and at
  // time n+1 for all other schemes
  // ---------------------------------------------------------------------

  // fine-scale velocity at time n+alpha_F/n+1
  LINALG::Matrix<nsd_,nen_> fsevelaf(true);
  if (fssgv_ != Fluid3::no_fssgv)
  {
    ExtractValuesFromGlobalVector(discretization,lm, &fsevelaf, NULL,"fsvelaf");
  }


  /*
  for (int i=0;i<nen_;++i)
  {
    for (int idim=0; idim<nsd_;++idim)
    {
      // velocity at n+alpha_F or n+1
      //evelaf(idim,i) = myvelaf[idim+(i*numdofpernode_)];

      // momentum equation part of history vector
      // (containing information of time step t_n (mass rhs!))
      emhist(idim,i) = myhist[idim+(i*numdofpernode_)];
    }

    // pressure at n+alpha_F or n+1
    //epreaf(i) = myvelaf[nsd_+(i*numdofpernode_)];

    // scalar at n+alpha_F or n+1
    escaaf(i) = myscaaf[nsd_+(i*numdofpernode_)];

    // scalar at n+alpha_M or n
    escaam(i) = myscaam[nsd_+(i*numdofpernode_)];

    // scalar time derivative at n+alpha_M or n+1
    escadtam(i) = myaccam[nsd_+(i*numdofpernode_)];
  }
*/
  // ---------------------------------------------------------------------
  // get additional state vector for generalized-alpha scheme or OST/BDF2:
  // acceleration at time n+alpha_M or velocity at time n
  // ---------------------------------------------------------------------

  // create objects for element arrays
  //LINALG::Matrix<nsd_,nen_> eaccam(true);
  //LINALG::Matrix<nsd_,nen_> eveln(true);

  /*
  if (is_genalpha_)
  {
    for (int i=0;i<nen_;++i)
    {
      for (int idim=0; idim<nsd_;++idim)
      {
      // acceleration at n+alpha_M
      eaccam(idim,i) = myaccam[idim+(i*numdofpernode_)];
      }
    }
  }
  else
  {
    for (int i=0;i<nen_;++i)
    {
      for (int idim=0; idim<nsd_;++idim)
      {
      // acceleration at n+alpha_M
      eveln(idim,i) = myscaam[idim+(i*numdofpernode_)];
      }
    }
  }
*/
  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  //RCP<const Epetra_Vector> dispnp;
  //vector<double> mydispnp;
  //RCP<const Epetra_Vector> gridv;
  //vector<double> mygridv;

  // create objects for element arrays
  //LINALG::Matrix<nsd_, nen_> edispnp(true);
  //LINALG::Matrix<nsd_, nen_> egridv(true);

  /*
  if (ele->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vectors 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);

    for (int i=0;i<nen_;++i)
    {
      for (int idim=0; idim<nsd_;++idim)
      {
      // set grid displacements
      edispnp(idim,i) = mydispnp[idim+(i*numdofpernode_)];

      // set grid velocities
      egridv(idim,i) = mygridv[idim+(i*numdofpernode_)];
      }
    }
  }
*/
  // ---------------------------------------------------------------------
  // get additional state vector for AVM3 case: fine-scale velocity
  // values are at time n+alpha_F for generalized-alpha scheme and at
  // time n+1 for all other schemes
  // ---------------------------------------------------------------------
/*
  // fine-scale velocity at time n+alpha_F/n+1
  RCP<const Epetra_Vector> fsvelaf;
  vector<double> myfsvelaf;

  // create object for element array
  LINALG::Matrix<nsd_,nen_> fsevelaf(true);

  if (fssgv_ != Fluid3::no_fssgv)
  {
    fsvelaf = discretization.GetState("fsvelaf");
    if (fsvelaf==null) dserror("Cannot get state vector 'fsvelaf'");
    myfsvelaf.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelaf,myfsvelaf,lm);

    for (int i=0;i<nen_;++i)
    {
      for (int idim=0; idim<nsd_;++idim)
      {
        // get fine-scale velocity
        fsevelaf(idim,i) = myfsvelaf[idim+(i*numdofpernode_)];
      }
    }
  }
*/
  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
  ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

  // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
  double Cs            = 0.0;
  double Cs_delta_sq   = 0.0;
  double l_tau         = 0.0;
  visceff_  = 0.0;

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  nlayer=0;

  GetTurbulenceParams(ele,
                      turbmodelparams,
                      Cs,
                      Cs_delta_sq,
                      l_tau,
                      nlayer);

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------

  bool mixed_formulation = params.get<bool>("mixed_formulation",false);

  if (mixed_formulation == false)
  {
    if (nsd_ == 3)
    {
      Sysmat(ele,
             evelaf,
             eveln,
             fsevelaf,
             epreaf,
             eaccam,
             escaaf,
             escaam,
             escadtam,
             emhist,
             edispnp,
             egridv,
             elemat1,
             elemat2,
             elevec1,
             thermpressaf,
             thermpressam,
             thermpressdtam,
             mat,
             timefac,
             Cs,
             Cs_delta_sq,
             l_tau);

#if 0
      if(ele->Id()==100 && 1)
      {
        FDcheck(ele,
                evelaf,
                eveln,
                fsevelaf,
                epreaf,
                eaccam,
                escaaf,
                escaam,
                escadtam,
                emhist,
                edispnp,
                egridv,
                elemat1,
                elemat2,
                elevec1,
                thermpressaf,
                thermpressam,
                thermpressdtam,
                mat,
                timefac,
                Cs,
                Cs_delta_sq,
                l_tau);
      }
#endif
    }
    else if (nsd_==2)
    {
      Sysmat2D(ele,
               evelaf,
               eveln,
               fsevelaf,
               epreaf,
               eaccam,
               escaaf,
               escaam,
               escadtam,
               emhist,
               edispnp,
               egridv,
               elemat1,
               elemat2,
               elevec1,
               thermpressaf,
               thermpressam,
               thermpressdtam,
               mat,
               timefac,
               Cs);
    }
    else dserror("Dimension is not working");
  }
  else
  {
    Sysmat2D3D(ele,
               evelaf,
               eveln,
               fsevelaf,
               epreaf,
               eaccam,
               escaaf,
               escaam,
               escadtam,
               emhist,
               edispnp,
               egridv,
               elemat1,
               elemat2,
               elevec1,
               thermpressaf,
               thermpressam,
               thermpressdtam,
               mat,
               timefac,
               Cs,
               Cs_delta_sq,
               l_tau);

#if 0
    if(ele->Id()==100 && 1)
    {
      FDcheck(ele,
              evelaf,
              eveln,
              fsevelaf,
              epreaf,
              eaccam,
              escaaf,
              escaam,
              escadtam,
              emhist,
              edispnp,
              egridv,
              elemat1,
              elemat2,
              elevec1,
              thermpressaf,
              thermpressam,
              thermpressdtam,
              mat,
              timefac,
              Cs,
              Cs_delta_sq,
              l_tau);
    }
#endif
  }

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
          (*(turbmodelparams.get<RCP<vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff_;
        }
      }
    }
  }

  //rotate matrices and vectors if we have a rotationally symmetric problem
  rotsymmpbc_->RotateMatandVecIfNecessary(elemat1,elemat2,elevec1);

  return 0;
}

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Sysmat(
  Fluid3*                                     ele,
  const LINALG::Matrix<nsd_,nen_>&             evelaf,
  const LINALG::Matrix<nsd_,nen_>&             eveln,
  const LINALG::Matrix<nsd_,nen_>&             fsevelaf,
  const LINALG::Matrix<nen_,1>&                epreaf,
  const LINALG::Matrix<nsd_,nen_>&             eaccam,
  const LINALG::Matrix<nen_,1>&                escaaf,
  const LINALG::Matrix<nen_,1>&                escaam,
  const LINALG::Matrix<nen_,1>&                escadtam,
  const LINALG::Matrix<nsd_,nen_>&             emhist,
  const LINALG::Matrix<nsd_,nen_>&             edispnp,
  const LINALG::Matrix<nsd_,nen_>&             egridv,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
  LINALG::Matrix<(nsd_+1)*nen_,1>&             eforce,
  const double                                 thermpressaf,
  const double                                 thermpressam,
  const double                                 thermpressdtam,
  Teuchos::RCP<const MAT::Material>           material,
  const double                                timefac,
  double&                                     Cs,
  double&                                     Cs_delta_sq,
  double&                                     l_tau
  )
{
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<nen_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }

  // add displacement when fluid nodes move in the ALE case
  if (ele->IsAle()) xyze_ += edispnp;

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(ele->Id());

  // element aera or volume
  const double vol = fac_;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab_ == INPAR::FLUID::viscous_stab_usfem or
      vstab_ == INPAR::FLUID::viscous_stab_usfem_only_rhs)   vstabfac =  1.0;
  else if(vstab_ == INPAR::FLUID::viscous_stab_gls or
          vstab_ == INPAR::FLUID::viscous_stab_gls_only_rhs) vstabfac = -1.0;

  //----------------------------------------------------------------------
  // get material parameters at element center
  //----------------------------------------------------------------------
  if (not mat_gp_ or not tau_gp_)
    GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

  if (not tau_gp_)
  {
    // ---------------------------------------------------------------------
    // calculate all-scale or fine-scale subgrid viscosity at element center
    // ---------------------------------------------------------------------
    visceff_ = visc_;
    if (turb_mod_action_ != Fluid3::no_model)
    {
      CalcSubgrVisc(evelaf,vol,Cs,Cs_delta_sq,l_tau);

      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fssgv_ != Fluid3::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,Cs);

    // get velocity at element center
    velint_.Multiply(evelaf,funct_);

    // ---------------------------------------------------------------------
    // calculate stabilization parameter at element center
    // ---------------------------------------------------------------------
    CalcStabParameter(timefac,vol);
  }

  // Gaussian integration points
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    if (tau_gp_)
    {
      // ---------------------------------------------------------------------
      // calculate all-scale or fine-scale subgrid viscosity at element center
      // ---------------------------------------------------------------------
      visceff_ = visc_;
      if (turb_mod_action_ != Fluid3::no_model)
      {
        CalcSubgrVisc(evelaf,vol,Cs,Cs_delta_sq,l_tau);

        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fssgv_ != Fluid3::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,Cs);

      // ---------------------------------------------------------------------
      // calculate stabilization parameter at element center
      // ---------------------------------------------------------------------
      CalcStabParameter(timefac,vol);
    }

    // get momentum history data at integration point
    histmom_.Multiply(emhist,funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf,derxy_);

    // get fine-scale velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (fssgv_ != Fluid3::no_fssgv) fsvderxy_.MultiplyNT(fsevelaf,derxy_);
    else                           fsvderxy_.Clear();

    // get convective velocity at integration point
    // We handle the ale case very implicitly here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_.Update(velint_);
    if (ele->IsAle()) convvelint_.Multiply(-1.0, egridv, funct_, 1.0);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    gradp_.Multiply(derxy_,epreaf);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = funct_.Dot(epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(edeadaf_,funct_);

    //--------------------------------------------------------------------
    // get numerical representation of some single operators
    //--------------------------------------------------------------------
    if (is_higher_order_ele_)
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
      if(physicaltype_ == INPAR::FLUID::loma)
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

      visc_old_(0) = viscs2_(0,0)*evelaf(0,0)+viscs2_(1,0)*evelaf(1,0)+viscs2_(2,0)*evelaf(2,0);
      visc_old_(1) = viscs2_(3,0)*evelaf(0,0)+viscs2_(4,0)*evelaf(1,0)+viscs2_(5,0)*evelaf(2,0);
      visc_old_(2) = viscs2_(6,0)*evelaf(0,0)+viscs2_(7,0)*evelaf(1,0)+viscs2_(8,0)*evelaf(2,0);

      for (int i=1; i<nen_; ++i)
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

        visc_old_(0) += viscs2_(0,i)*evelaf(0,i)+viscs2_(1,i)*evelaf(1,i)+viscs2_(2,i)*evelaf(2,i);
        visc_old_(1) += viscs2_(3,i)*evelaf(0,i)+viscs2_(4,i)*evelaf(1,i)+viscs2_(5,i)*evelaf(2,i);
        visc_old_(2) += viscs2_(6,i)*evelaf(0,i)+viscs2_(7,i)*evelaf(1,i)+viscs2_(8,i)*evelaf(2,i);
      }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }

    // convective term from previous iteration
    conv_old_.Multiply(vderxy_,convvelint_);

    // compute convective operator
    conv_c_.MultiplyTN(derxy_,convvelint_);

    // velocity divergence from previous iteration
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2);

    //--------------------------------------------------------------------
    // factors for stabilization, time integration
    // and fine-scale subgrid-viscosity
    //--------------------------------------------------------------------
    const double tau_M       = tau_(0)*fac_;
    const double tau_Mp      = tau_(1)*fac_;
    const double tau_C       = tau_(2)*fac_;

    const double timefacfac  = timefac * fac_;
    const double timetauM    = timefac * tau_M;
    const double timetauMp   = timefac * tau_Mp;

    double rhsfac            = fac_;

    const double fssgviscfac = fssgvisc_*timefacfac;

    //--------------------------------------------------------------------
    // The following computations are performed depending on
    // time-integration, that is, whether it is generalized-alpha or not,
    // since several terms differ with respect to the scheme:
    //
    // 1) calculation of rhs for momentum equation and momentum residual
    // -> different for generalized-alpha and other schemes
    //
    // 2) calculation of additional subgrid-scale velocity when cross-
    //    and Reynolds-stress are included:
    // - Cross- and Reynolds-stress are always included simultaneously.
    // - They are included in a complete form on left- and right-hand side.
    // - For this purpose, a subgrid-scale convective term is computed.
    // - Within a Newton linearization, the present formulation is not
    //   consistent for the reactive terms.
    // - To turn them off, both flags must be "no".
    //
    // 3) calculation of convective scalar term, rhs for continuity
    //    equation and residual of continuity equations
    // -> only required for low-Mach-number flow
    // -> for incompressible flow, residual is velocity divergence only
    //--------------------------------------------------------------------
    if (is_genalpha_)
    {
      // rhs of momentum equation: density*bodyforce at n+alpha_F
      rhsmom_.Update(densaf_,bodyforce_,0.0);

      // get acceleration at time n+alpha_M at integration point
      accint_.Multiply(eaccam,funct_);

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<3;++rr)
      {
        momres_old_(rr) = densam_*accint_(rr)+densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr)-densaf_*bodyforce_(rr);
      }

      if (cross_    != INPAR::FLUID::cross_stress_stab_none or
          reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-tau_(1),momres_old_,0.0);

        // compute subgrid-scale convective operator
        sgconv_c_.MultiplyTN(derxy_,sgvelint_);

        // re-calculate convective term from previous iteration if cross-stress
        // is included
        convvelint_.Update(1.0,sgvelint_,1.0);
        conv_old_.Multiply(vderxy_,convvelint_);
      }
      else sgconv_c_.Clear();

      // "incompressible" part of continuity residual: velocity divergence
      conres_old_ = vdiv_;

      if(physicaltype_ == INPAR::FLUID::loma)
      {
        // time derivative of scalar at n+alpha_M
        const double tder_sca = funct_.Dot(escadtam);

        // gradient of scalar value at n+alpha_F
        grad_scaaf_.Multiply(derxy_,escaaf);

        // convective scalar term at n+alpha_F
        conv_scaaf_ = velint_.Dot(grad_scaaf_);

        // add subgrid-scale velocity part also to convective scalar term
        // -> currently not taken into account
        /*if (cross    != Fluid3::cross_stress_stab_none or
            reynolds != Fluid3::reynolds_stress_stab_none)
          conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);*/

        // rhs of continuity equation (only relevant for low-Mach-number flow)
        rhscon_ = scadtfac_*tder_sca + scaconvfacaf_*conv_scaaf_ + thermpressadd_;

        // residual of continuity equation
        conres_old_ -= rhscon_;
      }
    }
    else
    {
      // rhs of momentum equation:
      // density*timefac*bodyforce at n+1 + density*histmom at n

      // in the case of a Boussinesq approximation: f = (rho - rho_0)/rho_0 *g
      // else: 										f = rho * g

      if (physicaltype_ == INPAR::FLUID::boussinesq)
    	rhsmom_.Update(densn_,histmom_,deltadens_*timefac,bodyforce_);
      else rhsmom_.Update(densn_,histmom_,densaf_*timefac,bodyforce_);

      // modify integration factor for Galerkin rhs
      rhsfac *= timefac;

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<3;++rr)
      {
        momres_old_(rr) = densaf_*velint_(rr)+timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))-rhsmom_(rr);
      }

      if (cross_    != INPAR::FLUID::cross_stress_stab_none or
          reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-(tau_(1)/dt_),momres_old_,0.0);

        // compute subgrid-scale convective operator
        sgconv_c_.MultiplyTN(derxy_,sgvelint_);

        // re-calculate convective term from previous iteration if cross-stress
        // is included
        convvelint_.Update(1.0,sgvelint_,1.0);
        conv_old_.Multiply(vderxy_,convvelint_);
      }
      else sgconv_c_.Clear();

      // "incompressible" part of continuity residual: velocity divergence
      conres_old_ = timefac*vdiv_;

      if (physicaltype_ == INPAR::FLUID::loma or physicaltype_ == INPAR::FLUID::varying_density)
      {
        // get velocity derivatives at n
        vderxyn_.MultiplyNT(eveln,derxy_);

        // velocity divergence at n
        const double vdivn = vderxyn_(0, 0) + vderxyn_(1, 1) + vderxyn_(2, 2);

        // time derivative of scalar
        const double tder_sca = funct_.Dot(escadtam);

        // gradient of scalar value at n+1
        grad_scaaf_.Multiply(derxy_,escaaf);

        // convective scalar term at n+1
        conv_scaaf_ = velint_.Dot(grad_scaaf_);

        // gradient of scalar value at n
        grad_scan_.Multiply(derxy_,escaam);

        // convective scalar term at n
        conv_scan_ = velint_.Dot(grad_scan_);

        // add subgrid-scale velocity part also to convective scalar term
        // (subgrid-scale velocity at n+1 also approximately used at n)
        // -> currently not taken into account
        /*if (cross    != Fluid3::cross_stress_stab_none or
            reynolds != Fluid3::reynolds_stress_stab_none)
        {
          conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);
          conv_scan_  += sgvelint_.Dot(grad_scan_);
        }*/

        // rhs of continuity equation (only relevant for low-Mach-number flow)
        rhscon_ = scadtfac_*tder_sca + timefac*scaconvfacaf_*conv_scaaf_ + omtheta_*dt_*(scaconvfacn_*conv_scan_-vdivn) + thermpressadd_;

        // residual of continuity equation
        conres_old_ -= rhscon_;
      }
    }

    //------------------------------------------------------------------------
    // perform integration for element matrix and right hand side
    //------------------------------------------------------------------------
    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      //----------------------------------------------------------------------
      // computation of inertia term and convection term (convective and
      // reactive part) for convective form of convection term including
      // right-hand-side contribution and potential cross-stress term
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;
        const double v = fac_*densam_*funct_(ui)
#if 1
                         + timefacfac*densaf_*(conv_c_(ui)+sgconv_c_(ui))
#endif
                         ;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* inertia (contribution to mass matrix) */
          /*
          /                \
          |                |
          |    rho*Du , v  |
          |                |
          \                /
          */
          /* convection, convective part (convective form) */
          /*
          /                               \
          |  /       n+1        \         |
          | |   rho*u   o nabla | Du , v  |
          |  \      (i)        /          |
          \                              /
          */
          double v2 = v*funct_(vi) ;
          estif(fvi  , fui  ) += v2;
          estif(fvip , fuip ) += v2;
          estif(fvipp, fuipp) += v2;
        }
      }

      if (is_newton_)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          const double v = timefacfac*densaf_*funct_(vi);
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v2 = v*funct_(ui);
            /*  convection, reactive part (convective form)
            /                                 \
            |  /                \   n+1       |
            | |  rho*Du o nabla | u     , v   |
            |  \                /   (i)       |
            \                                /
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

      if (is_genalpha_)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* inertia term on right-hand side for generalized-alpha scheme */
          const double v = -fac_*densam_*funct_(vi);
          eforce(fvi    ) += v*accint_(0) ;
          eforce(fvi + 1) += v*accint_(1) ;
          eforce(fvi + 2) += v*accint_(2) ;
        }
      }
      else
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* inertia term on right-hand side for one-step-theta/BDF2 schem */
          const double v = -fac_*densaf_*funct_(vi);
          eforce(fvi    ) += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
          eforce(fvi + 2) += v*velint_(2) ;
        }
      }

#if 1
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = 4*vi;
        /* convection (convective form) on right-hand side */
        double v = -rhsfac*densaf_*funct_(vi);
        eforce(fvi    ) += v*conv_old_(0) ;
        eforce(fvi + 1) += v*conv_old_(1) ;
        eforce(fvi + 2) += v*conv_old_(2) ;
      }
#endif

      //----------------------------------------------------------------------
      // computation of additions to convection term (convective and
      // reactive part) for conservative form of convection term including
      // right-hand-side contribution
      //----------------------------------------------------------------------
      if (is_conservative_)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v = timefacfac*densaf_*funct_(ui)*vdiv_;
          if (physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
          // only with linear density-concentration correlation
          else if(physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* convection, convective part (conservative addition) */
            /*
            /                                                   \
            |      /              n+1    n+1           \      |
            |  Du | rho*nabla o u    +  u   *nabla rho | , v  |
            |      \             (i)     (i)          /       |
            \                                                 /
            */
            double v2 = v*funct_(vi) ;
            estif(fvi  , fui  ) += v2;
            estif(fvip , fuip ) += v2;
            estif(fvipp, fuipp) += v2;
          }
        }

        if (is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
            const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
            const double v2 = timefacfac*densaf_*velint_(2)*funct_(vi);
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              /*  convection, reactive part (conservative addition) */
              /*
              /                              \
              |  n+1  /               \      |
              | u    | rho*nabla o Du | , v  |
              |  (i)  \              /       |
              \                             /
              */
              estif(fvi,  fui  ) += v0*derxy_(0, ui) ;
              estif(fvi,  fuip ) += v0*derxy_(1, ui) ;
              estif(fvi,  fuipp) += v0*derxy_(2, ui) ;
              estif(fvip, fui  ) += v1*derxy_(0, ui) ;
              estif(fvip, fuip ) += v1*derxy_(1, ui) ;
              estif(fvip, fuipp) += v1*derxy_(2, ui) ;
              estif(fvipp,fui  ) += v2*derxy_(0, ui) ;
              estif(fvipp,fuip ) += v2*derxy_(1, ui) ;
              estif(fvipp,fuipp) += v2*derxy_(2, ui) ;
            }
          }

          if (physicaltype_ == INPAR::FLUID::loma)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
              const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
              const double v2 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(2)*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = 4*ui;
                const int fuip  = fui+1;
                const int fuipp = fui+2;
                /*  convection, reactive part (conservative addition) */
                /*
                /                           \
                |  n+1  /             \      |
                | u    | Du*nabla rho | , v  |
                |  (i)  \            /       |
                \                           /
                */
                estif(fvi,  fui  ) += v0*funct_(ui) ;
                estif(fvi,  fuip ) += v0*funct_(ui) ;
                estif(fvi,  fuipp) += v0*funct_(ui) ;
                estif(fvip, fui  ) += v1*funct_(ui) ;
                estif(fvip, fuip ) += v1*funct_(ui) ;
                estif(fvip, fuipp) += v1*funct_(ui) ;
                estif(fvipp,fui  ) += v2*funct_(ui) ;
                estif(fvipp,fuip ) += v2*funct_(ui) ;
                estif(fvipp,fuipp) += v2*funct_(ui) ;
              }
            }
          }
          if (physicaltype_ == INPAR::FLUID::varying_density)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
			  const int fvi   = 4*vi;
			  const int fvip  = fvi+1;
			  const int fvipp = fvi+2;
			  const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
			  const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
			  const double v2 = +timefacfac*grad_scaaf_(2)*velint_(2)*funct_(vi);
			  for (int ui=0; ui<nen_; ++ui)
			  {
				const int fui   = 4*ui;
				const int fuip  = fui+1;
				const int fuipp = fui+2;
				/*  convection, reactive part (conservative addition) */
				/*
				/                           \
				|  n+1  /             \      |
				| u    | Du*nabla rho | , v  |
				|  (i)  \            /       |
				\                           /
				*/
				estif(fvi,  fui  ) += v0*funct_(ui) ;
				estif(fvi,  fuip ) += v0*funct_(ui) ;
				estif(fvi,  fuipp) += v0*funct_(ui) ;
				estif(fvip, fui  ) += v1*funct_(ui) ;
				estif(fvip, fuip ) += v1*funct_(ui) ;
				estif(fvip, fuipp) += v1*funct_(ui) ;
				estif(fvipp,fui  ) += v2*funct_(ui) ;
				estif(fvipp,fuip ) += v2*funct_(ui) ;
				estif(fvipp,fuipp) += v2*funct_(ui) ;
			  }
			}
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          /* convection (conservative addition) on right-hand side */
          double v = -rhsfac*densaf_*funct_(vi)*vdiv_;
          eforce(fvi    ) += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
          eforce(fvi + 2) += v*velint_(2) ;
        }

        if (physicaltype_ == INPAR::FLUID::loma)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = rhsfac*densaf_*scaconvfacaf_*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
            eforce(fvi + 2) += v*velint_(2) ;
          }
        }
        if (physicaltype_ == INPAR::FLUID::varying_density)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = -rhsfac*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
            eforce(fvi + 2) += v*velint_(2) ;
          }
        }
      }

      //----------------------------------------------------------------------
      // computation of viscosity term including right-hand-side contribution
      //----------------------------------------------------------------------
      const double visceff_timefacfac = visceff_*timefacfac;
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;
        for (int vi=0; vi<nen_; ++vi)
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

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        const double v = -visceff_*rhsfac;
        /* viscosity term on right-hand side */
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

      //----------------------------------------------------------------------
      // computation of pressure term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fuippp = 4*ui+3;
        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          /* pressure term */
          /*
          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */
          estif(fvi,     fuippp) += v*derxy_(0, vi) ;
          estif(fvi + 1, fuippp) += v*derxy_(1, vi) ;
          estif(fvi + 2, fuippp) += v*derxy_(2, vi) ;
        }
      }

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        /* pressure term on right-hand side */
        const double v = press*rhsfac;
        eforce(fvi    ) += v*derxy_(0, vi) ;
        eforce(fvi + 1) += v*derxy_(1, vi) ;
        eforce(fvi + 2) += v*derxy_(2, vi) ;
      }

      //----------------------------------------------------------------------
      // computation of continuity term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvippp = 4*vi+3;
        const double v = timefacfac*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          /* continuity term */
          /*
            /                  \
            |                  |
            | nabla o Du  , q  |
            |                  |
            \                  /
          */
          estif(fvippp, fui)     += v*derxy_(0, ui) ;
          estif(fvippp, fui + 1) += v*derxy_(1, ui) ;
          estif(fvippp, fui + 2) += v*derxy_(2, ui) ;
        }
      }

      const double rhsfac_vdiv = -rhsfac * vdiv_;
      for (int vi=0; vi<nen_; ++vi)
      {
        // continuity term on right-hand side
        eforce(vi*4 + 3) += rhsfac_vdiv*funct_(vi) ;
      }

      //----------------------------------------------------------------------
      // computation of body-force term on right-hand side
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        const double v = fac_*funct_(vi);
        eforce(fvi    ) += v*rhsmom_(0) ;
        eforce(fvi + 1) += v*rhsmom_(1) ;
        eforce(fvi + 2) += v*rhsmom_(2) ;
      }

      //----------------------------------------------------------------------
      // computation of additional terms for low-Mach-number flow:
      // 1) subtracted viscosity term including right-hand-side contribution
      // 2) additional rhs term of continuity equation
      //----------------------------------------------------------------------
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        const double v = -(2.0/3.0)*visceff_*timefacfac ;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v0 = v*derxy_(0,ui);
          const double v1 = v*derxy_(1,ui);
          const double v2 = v*derxy_(2,ui);
          for (int vi=0; vi<nen_; ++vi)
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
            estif(fvi,   fui  ) += v0*derxy_(0, vi) ;
            estif(fvi,   fuip ) += v1*derxy_(0, vi) ;
            estif(fvi,   fuipp) += v2*derxy_(0, vi) ;
            estif(fvip,  fui  ) += v0*derxy_(1, vi) ;
            estif(fvip,  fuip ) += v1*derxy_(1, vi) ;
            estif(fvip,  fuipp) += v2*derxy_(1, vi) ;
            estif(fvipp, fui  ) += v0*derxy_(2, vi) ;
            estif(fvipp, fuip ) += v1*derxy_(2, vi) ;
            estif(fvipp, fuipp) += v2*derxy_(2, vi) ;
          }
        }


        if (is_newton_)
        {
          const double timefacfac_scaconvfacaf=timefacfac*scaconvfacaf_;

          LINALG::Matrix<nsd_,1> temp;

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp= numdofpernode_*vi+nsd_;
            const double timefacfac_scaconvfacaf_funct_vi=timefacfac_scaconvfacaf*funct_(vi);

            for(int jdim=0;jdim<nsd_;++jdim)
            {
              temp(jdim)=timefacfac_scaconvfacaf_funct_vi*grad_scaaf_(jdim);
            }

            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui=numdofpernode_*ui;

              for(int jdim=0;jdim<nsd_;++jdim)
              {
                /*
                  factor afgtd/am

                          /                    \
                    1    |       /         \    |
                   --- * |  q , | Du o grad | T |
                    T    |       \         /    |
                          \                    /
                */
                estif(fvippp,fui+jdim) -= temp(jdim)*funct_(ui);
              }
            }
          }
        } // end if (is_newton_)

        const double v_div = (2.0/3.0)*visceff_*rhsfac*vdiv_ ;
        const double fac_rhscon = fac_*rhscon_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* viscosity term on rhs - subtraction for low-Mach-number flow */
          eforce(fvi    ) += derxy_(0, vi)*v_div ;
          eforce(fvi + 1) += derxy_(1, vi)*v_div ;
          eforce(fvi + 2) += derxy_(2, vi)*v_div ;

          /* additional rhs term of continuity equation */
          eforce(fvi + 3) += fac_rhscon*funct_(vi) ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg_ == INPAR::FLUID::pstab_use_pspg)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v = tau_Mp*densam_*funct_(ui)
#if 1
                     + timetauMp*densaf_*conv_c_(ui)
#endif
                     ;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp = 4*vi+3;

            /* pressure stabilisation: inertia */
            /*
              /                    \
              |                    |
              |  rho*Du , nabla q  |
              |                    |
              \                   /
            */
            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | |   rho*u   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */

            estif(fvippp, fui)   += v*derxy_(0, vi) ;
            estif(fvippp, fuip)  += v*derxy_(1, vi) ;
            estif(fvippp, fuipp) += v*derxy_(2, vi) ;
          }
        }

        if (is_higher_order_ele_)
        {
          const double v = -2.0*visceff_*timetauMp;
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            for (int vi=0; vi<nen_; ++vi)
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

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fuippp = 4*ui+3;
          for (int vi=0; vi<nen_; ++vi)
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

        if (is_newton_)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauMp*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvippp = 4*vi + 3;
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | |   rho*Du o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

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

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          eforce(vi*4 + 3) -= tau_Mp*(momres_old_(0)*derxy_(0, vi)
                                      +
                                      momres_old_(1)*derxy_(1, vi)
                                      +
                                      momres_old_(2)*derxy_(2, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg_ == INPAR::FLUID::convective_stab_supg)
      {
#if 1
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v = densaf_*(tau_M*densam_*funct_(ui) + timetauM*densaf_*conv_c_(ui));
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* supg stabilisation: inertia  */
            /*
              /                                 \
              |           /      n+1        \    |
              |  rho*Du , | rho*u   o nabla | v  |
              |           \      (i)        /    |
              \                                 /
            */

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                      \
            |    /       n+1         \      /      n+1        \     |
            |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
            |    \       (i)         /      \      (i)        /     |
            \                                                       /

            */

            const double v2 = v*(conv_c_(vi)+sgconv_c_(vi));
            estif(fvi,   fui)   += v2;
            estif(fvip,  fuip)  += v2;
            estif(fvipp, fuipp) += v2;
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          const double v = timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fuippp = 4*ui + 3;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , |   rho*u   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(fvi,   fuippp) += v*derxy_(0, ui) ;
            estif(fvip,  fuippp) += v*derxy_(1, ui) ;
            estif(fvipp, fuippp) += v*derxy_(2, ui) ;
          }
        }

        if (is_higher_order_ele_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v = -2.0*visceff_*timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, |   rho*u    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
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

        if (is_newton_)
        {
          if(is_genalpha_)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = timetauM*densam_*densaf_*funct_(ui);
              const double v0 = v*accint_(0);
              const double v1 = v*accint_(1);
              const double v2 = v*accint_(2);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                /* supg stabilisation: inertia, linearisation of testfunction  */
                /*
                /                                         \
                |         n+1      /                \     |
                |    rho*u      ,  | rho*Du o nabla | v   |
                |         (i)      \                /     |
                \                                         /

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


            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = -timetauM*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                const double v0 = v*rhsmom_(0);
                const double v1 = v*rhsmom_(1);
                const double v2 = v*rhsmom_(2);

                /* supg stabilisation: bodyforce part, linearisation of test function */

                /*
                /                                       \
                |                 /                \     |
                |  rho*rhsint   , |  rho*Du o nabla | v  |
                |                 \                /     |
                \                                        /

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
          } // end is_genalpha
          else
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = tau_M*densam_*densaf_*funct_(ui);
              const double v0 = v*velint_(0);
              const double v1 = v*velint_(1);
              const double v2 = v*velint_(2);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                /* supg stabilisation: inertia, linearisation of testfunction  */
                /*
                /                                         \
                |         n+1      /                \     |
                |    rho*u      ,  | rho*Du o nabla | v   |
                |         (i)      \                /     |
                \                                         /

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


            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = -tau_M*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                const double v0 = v*rhsmom_(0);
                const double v1 = v*rhsmom_(1);
                const double v2 = v*rhsmom_(2);

                /* supg stabilisation: bodyforce part, linearisation of test function */

                /*
                /                                       \
                |                 /                \     |
                |  rho*rhsint   , |  rho*Du o nabla | v  |
                |                 \                /     |
                \                                        /

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
          } // end not genalpha
#if 1
          {
            const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1) + convvelint_(2)*vderxy_(0, 2);
            const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1) + convvelint_(2)*vderxy_(1, 2);
            const double v2 = convvelint_(0)*vderxy_(2, 0) + convvelint_(1)*vderxy_(2, 1) + convvelint_(2)*vderxy_(2, 2);

            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = timetauM*densaf_*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                /*
                  /                                                         \
                  |    /       n+1        \   n+1    /                \     |
                  |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                  |    \       (i)        /   (i)    \                /     |
                  \                                                        /

                  /                                                         \
                  |    /                \   n+1    /       n+1        \     |
                  |    | rho*Du o nabla | u    ,   | rho*u   o nabla  | v   |
                  |    \                /   (i)    \       (i)        /     |
                  \                                                        /
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

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauM*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , |   rho*Du o nabla | v  |
                |         (i)    \                /     |
                \                                      /
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

          if (is_higher_order_ele_)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = -2.0*visceff_*timetauM*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                const double v0 = v*visc_old_(0);
                const double v1 = v*visc_old_(1);
                const double v2 = v*visc_old_(2);

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \                /     |
                  \                                                 /
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
        } // if newton

#if 1
        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          // supg stabilisation
          const double v = -tau_M*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          eforce(fvi)     += v*momres_old_(0) ;
          eforce(fvi + 1) += v*momres_old_(1) ;
          eforce(fvi + 2) += v*momres_old_(2) ;
        }
#endif
      }

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART

      if (is_higher_order_ele_)
      {
        if(vstab_ != INPAR::FLUID::viscous_stab_none)
        {
          const double two_visc_tauMp = vstabfac*2.0*visc_*tau_Mp;
          // viscous stabilization either on left hand side or on right hand side
          if (vstab_ == INPAR::FLUID::viscous_stab_gls || vstab_ == INPAR::FLUID::viscous_stab_usfem)
          {
            const double two_visc_timetauMp   = vstabfac*2.0*visc_*timetauMp;
            const double four_visc2_timetauMp = vstabfac*4.0*visceff_*visc_*timetauMp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<nen_; ++ui)
            {
              const double v = two_visc_tauMp*densam_*funct_(ui)
#if 1
                         + two_visc_timetauMp*densaf_*conv_c_(ui)
#endif
                         ;
              for (int vi=0; vi<nen_; ++vi)
              {
                /* viscous stabilisation, inertia part */
                /*
                  /                          \
                  |                          |
              +/- |    rho*Du , div eps (v)  |
                  |                          |
                  \                          /
                */
                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | |   rho*u   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
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

            for (int ui=0; ui<nen_; ++ui)
            {
              for (int vi=0; vi<nen_; ++vi)
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

            for (int ui=0; ui<nen_; ++ui)
            {
              for (int vi=0; vi<nen_; ++vi)
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

            if (is_newton_)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                double v = two_visc_timetauMp*densaf_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | |   rho*Du o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
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

          for (int vi=0; vi<nen_; ++vi)
          {

            /* viscous stabilisation */
            eforce(vi*4    ) -= two_visc_tauMp*(momres_old_(0)*viscs2_(0, vi)+momres_old_(1)*viscs2_(1, vi)+momres_old_(2)*viscs2_(2, vi)) ;
            eforce(vi*4 + 1) -= two_visc_tauMp*(momres_old_(0)*viscs2_(1, vi)+momres_old_(1)*viscs2_(4, vi)+momres_old_(2)*viscs2_(5, vi)) ;
            eforce(vi*4 + 2) -= two_visc_tauMp*(momres_old_(0)*viscs2_(2, vi)+momres_old_(1)*viscs2_(5, vi)+momres_old_(2)*viscs2_(8, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab_ == INPAR::FLUID::continuity_stab_yes)
      {
        const double timetauC = timefac*tau_C;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v0 = timetauC*derxy_(0, ui);
          double v1 = timetauC*derxy_(1, ui);
          double v2 = timetauC*derxy_(2, ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* continuity stabilisation on left hand side */
            /*
              /                         \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
              \                         /
            */
            estif(fvi,  fui  ) += v0*derxy_(0, vi) ;
            estif(fvip, fui  ) += v0*derxy_(1, vi) ;
            estif(fvipp,fui  ) += v0*derxy_(2, vi) ;

            estif(fvi,  fuip ) += v1*derxy_(0, vi) ;
            estif(fvip, fuip ) += v1*derxy_(1, vi) ;
            estif(fvipp,fuip ) += v1*derxy_(2, vi) ;

            estif(fvi,  fuipp) += v2*derxy_(0, vi) ;
            estif(fvip, fuipp) += v2*derxy_(1, vi) ;
            estif(fvipp,fuipp) += v2*derxy_(2, vi) ;
          }
        }

        const double tauC_conres = tau_C*conres_old_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* continuity stabilisation on right hand side */
          eforce(fvi  ) -= tauC_conres*derxy_(0, vi) ;
          eforce(fvip ) -= tauC_conres*derxy_(1, vi) ;
          eforce(fvipp) -= tauC_conres*derxy_(2, vi) ;
        }
      }

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv_ != Fluid3::no_fssgv)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(fvi    ) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                         +    derxy_(1, vi)*fsvderxy_(0, 1)
                                         +    derxy_(1, vi)*fsvderxy_(1, 0)
                                         +    derxy_(2, vi)*fsvderxy_(0, 2)
                                         +    derxy_(2, vi)*fsvderxy_(2, 0)) ;
          eforce(fvi + 1) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                         +    derxy_(0, vi)*fsvderxy_(1, 0)
                                         +2.0*derxy_(1, vi)*fsvderxy_(1, 1)
                                         +    derxy_(2, vi)*fsvderxy_(1, 2)
                                         +    derxy_(2, vi)*fsvderxy_(2, 1)) ;
          eforce(fvi + 2) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 2)
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
      for (int vi=0; vi<nen_; ++vi)
      {
        double v = fac_*funct_(vi,0);
        for (int ui=0; ui<nen_; ++ui)
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

      for (int ui=0; ui<nen_; ++ui)
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

        for (int vi=0; vi<nen_; ++vi)
        {
          double v = timefacfac/det_*funct_(vi);

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

      double v = visceff_*timefac*fac_;
      for (int ui=0; ui<nen_; ++ui)
      {
        double derinvJ0 = -v*(deriv_(0,ui)*xji_00 + deriv_(1,ui)*xji_01 + deriv_(2,ui)*xji_02);
        double derinvJ1 = -v*(deriv_(0,ui)*xji_10 + deriv_(1,ui)*xji_11 + deriv_(2,ui)*xji_12);
        double derinvJ2 = -v*(deriv_(0,ui)*xji_20 + deriv_(1,ui)*xji_21 + deriv_(2,ui)*xji_22);
        for (int vi=0; vi<nen_; ++vi)
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

      v = timefacfac*visceff_/det_;
      for (int ui=0; ui<nen_; ++ui)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
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

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }
      }


      // pressure
      for (int vi=0; vi<nen_; ++vi)
      {
        double v = press*timefacfac/det_;
        for (int ui=0; ui<nen_; ++ui)
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
      for (int vi=0; vi<nen_; ++vi)
      {
        double v = timefacfac/det_*funct_(vi,0);
        for (int ui=0; ui<nen_; ++ui)
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

  return;
}


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadaf only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::BodyForce(Fluid3* ele)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  if(nsd_==3)
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
  else if (nsd_==2)
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
  else
    dserror("Body force for a 1D problem is not yet implemented");

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
      if(time_ >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time_);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time_);
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
    const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the edeadaf array
    for(int isd=0;isd<nsd_;isd++)
    {
      // get factor given by spatial function
      if (functions) functnum = (*functions)[isd];
      else functnum = -1;

      double num = (*onoff)[isd]*(*val)[isd]*curvefac;

      for (int jnode=0; jnode<nen_; jnode++)
      {
        if (functnum>0)
        {
          // evaluate function at the position of the current node
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X(),time_,NULL);
        }
        else functionfac = 1.0;

        edeadaf_(isd,jnode) = num*functionfac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadaf_.Clear();
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at element center  vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::EvalShapeFuncAndDerivsAtEleCenter(
  const int  eleid
)
{
  // use one-point Gauss rule
  //DRT::UTILS::GaussRule3D integrationrule_stabili=DRT::UTILS::intrule3D_undefined;
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  /*
  switch (distype)
  {
  // 3D
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    intpoints_stab
    //integrationrule_stabili = DRT::UTILS::intrule_hex_1point;
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
  */

  // Gaussian points
  // const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim=0;idim<nsd_;idim++)
  {
	  xsi_(idim) = gpcoord[idim];
  }
  const double wquad = intpoints_stab.IP().qwgt[0];

  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // shape functions and derivs at element center
  //const double e1    = intpoints.qxg[0][0];
  //const double e2    = intpoints.qxg[0][1];
  //const double e3    = intpoints.qxg[0][2];
  //const double wquad = intpoints.qwgt[0];

  //DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  //DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // compute Jacobian matrix and determinant
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

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = wquad*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point   vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,		// integration points
    const int                              iquad,				// actual integration point
    const int                              eleid				// element ID
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {
	  xsi_(idim) = gpcoord[idim];
  }

  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // coordinates of the current integration point
  //const double e1 = intpoints.IP()qxg[iquad][0];
  //const double e2 = intpoints.qxg[iquad][1];
  //const double e3 = intpoints.qxg[iquad][2];

  // shape functions and their derivatives
  //DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  //DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

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
  det_ = xji_.Invert(xjm_);

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = intpoints.IP().qwgt[iquad]*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::shape_function_deriv2<distype>(xsi_, deriv2_);
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}


/*----------------------------------------------------------------------*
 |  compute material parameters                                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material>  material,
  const LINALG::Matrix<nsd_,nen_>&       evelaf,
  const LINALG::Matrix<nen_,1>&       escaaf,
  const LINALG::Matrix<nen_,1>&       escaam,
  const double                        thermpressaf,
  const double                        thermpressam,
  const double                        thermpressdtam
)
{
// initially set density values and values with respect to continuity rhs
// -> will only be changed for low-Mach-number flow "material" below
densam_        = 1.0;
densaf_        = 1.0;
densn_         = 1.0;
scadtfac_      = 0.0;
scaconvfacaf_  = 0.0;
scaconvfacn_   = 0.0;
thermpressadd_ = 0.0;

if (material->MaterialType() == INPAR::MAT::m_fluid)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

  // get constant viscosity
  visc_ = actmat->Viscosity();

  // Varying Density
  if (physicaltype_ == INPAR::FLUID::varying_density)
  {
	  densaf_ = funct_.Dot(escaaf);
	  densam_ = densaf_;

	  densn_ = funct_.Dot(escaam);
  }
  // Boussinesq approximation: Calculation of delta rho
  else if (physicaltype_ == INPAR::FLUID::boussinesq)
  {
	  const double density_0 = actmat->Density();
	  deltadens_ =  (funct_.Dot(escaaf)- density_0)/ density_0;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_carreauyasuda)
{
  const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(material.get());

  double nu_0   = actmat->Nu0();    // parameter for zero-shear viscosity
  double nu_inf = actmat->NuInf();  // parameter for infinite-shear viscosity
  double lambda = actmat->Lambda();  // parameter for characteristic time
  double a      = actmat->AParam(); // constant parameter
  double b      = actmat->BParam(); // constant parameter

  // compute rate of strain at n+alpha_F or n+1
  double rateofstrain   = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  const double tmp = pow(lambda*rateofstrain,b);
  visc_ = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
}
else if (material->MaterialType() == INPAR::MAT::m_modpowerlaw)
{
  const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(material.get());

  // get material parameters
  double m     = actmat->MCons();     // consistency constant
  double delta = actmat->Delta();      // safety factor
  double a     = actmat->AExp();      // exponent

  // compute rate of strain at n+alpha_F or n+1
  double rateofstrain   = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  // compute viscosity according to a modified power law model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  visc_ = m * pow((delta + rateofstrain), (-1)*a);
}
else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
{
  const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

  // compute mixture fraction at n+alpha_F or n+1
  const double mixfracaf = funct_.Dot(escaaf);

  // compute dynamic viscosity at n+alpha_F or n+1 based on mixture fraction
  visc_ = actmat->ComputeViscosity(mixfracaf);

  // compute density at n+alpha_F or n+1 based on mixture fraction
  densaf_ = actmat->ComputeDensity(mixfracaf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->EosFacA()*densaf_;

  if (is_genalpha_)
  {
    // compute density at n+alpha_M based on mixture fraction
    const double mixfracam = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(mixfracam);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->EosFacA()*densam_;
  }
  else
  {
    // compute density at n based on mixture fraction
    const double mixfracn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(mixfracn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->EosFacA()*densn_;

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = dt_*actmat->EosFacA()*densam_;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_sutherland)
{
  const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

  // compute temperature at n+alpha_F or n+1
  const double tempaf = funct_.Dot(escaaf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on temperature
  // and thermodynamic pressure
  densaf_ = actmat->ComputeDensity(tempaf,thermpressaf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = 1.0/tempaf;

  if (is_genalpha_)
  {
    // compute temperature at n+alpha_M
    const double tempam = funct_.Dot(escaam);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = 1.0/tempam;

    // compute density at n+alpha_M based on temperature
    densam_ = actmat->ComputeDensity(tempam,thermpressam);

    // addition due to thermodynamic pressure at n+alpha_M
    thermpressadd_ = -thermpressdtam/thermpressam;
  }
  else
  {
    // compute temperature at n
    const double tempn = funct_.Dot(escaam);

    // compute density at n based on temperature at n and
    // (approximately) thermodynamic pressure at n+1
    densn_ = actmat->ComputeDensity(tempn,thermpressaf);

    // factor for convective scalar term at n
    scaconvfacn_ = 1.0/tempn;

    // factor for scalar time derivative
    scadtfac_ = dt_/tempaf;

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // addition due to thermodynamic pressure
    thermpressadd_ = -dt_*thermpressdtam/thermpressaf;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
{
  const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->ComputeFactor(provaraf);

  if (is_genalpha_)
  {
    // compute density at n+alpha_M based on progress variable
    const double provaram = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(provaram);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->ComputeFactor(provaram);
  }
  else
  {
    // compute density at n based on progress variable
    const double provarn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(provarn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->ComputeFactor(provarn);

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = dt_*scaconvfacaf_;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
{
  const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->ComputeFactor(provaraf);

  if (is_genalpha_)
  {
    // compute density at n+alpha_M based on progress variable
    const double provaram = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(provaram);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->ComputeFactor(provaram);
  }
  else
  {
    // compute density at n based on progress variable
    const double provarn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(provarn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->ComputeFactor(provarn);

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = dt_*scaconvfacaf_;
  }
}
else dserror("Material type is not supported");

// check whether there is zero or negative (physical) viscosity
if (visc_ < EPS15) dserror("zero or negative (physical) diffusivity");

return;
} // Fluid3Impl::GetMaterialParams


/*----------------------------------------------------------------------*
 |  compute turbulence parameters                            ehrl 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::Fluid3Impl<distype>::SetSolutionParameter(
                                  Fluid3*                    ele,
                                  ParameterList&             params)
{
  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // definition of a local variable to compute the timefac
  double timefac = 0.0;

  // set flag, if it is a stationary fluid system
  is_stationary_ = params.get<bool>("is stationary");
  //TODO: check, when is_stationary must be set

  // check whether we have a generalized-alpha time-integration scheme
  is_genalpha_ = params.get<bool>("using generalized-alpha time integration");

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  time_ = params.get<double>("total time",-1.0);

  if (is_stationary_ == false)
  {
    // get time-step length and time-integration parameters
    dt_                = params.get<double>("dt");
    const double theta = params.get<double>("theta",-1.0);
    omtheta_           = params.get<double>("omtheta",-1.0);

    // compute timefactor for left-hand side
    // One-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
    timefac = theta*dt_;

    if(is_genalpha_)
    {
      gamma_ =params.get<double>("gamma");
      alphaF_=params.get<double>("alphaF");
      alphaM_=params.get<double>("alphaM");
    }
    else
    {
      gamma_ =theta;
      alphaF_=1.0;
      alphaM_=1.0;
    }

    afgdt_=alphaF_*gamma_*dt_;
  }
  else
  {
    if (ele->IsAle()) dserror("No ALE support within stationary fluid solver.");

    // timefac stationary = 1.0
    timefac = 1.0;
  }

  if (timefac < 0.0 or time_ < 0.0)
    dserror("Negative time-integration parameter or time-step length supplied");

  // ---------------------------------------------------------------------
  // get control parameters for linearization, physical type of fluid flow for solver
  // and form of convective term
  //----------------------------------------------------------------------
  string newtonstr   = params.get<string>("Linearisation");
  string convformstr = params.get<string>("form of convective term");
  if (newtonstr=="Newton")          is_newton_       = true;
  if (convformstr =="conservative") is_conservative_ = true;

  // set flag physical type
  physicaltype_ = params.get<INPAR::FLUID::PhysicalType>("Physical Type");
  if (((physicaltype_ != INPAR::FLUID::boussinesq) and (physicaltype_ != INPAR::FLUID::incompressible))
      and (is_stationary_ == true))
    dserror("physical type not supported in stationary FLUID3 implementation.");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  ParameterList& stablist = params.sublist("STABILIZATION");

  // no safety check necessary since all options are used
  tds_      = Teuchos::getIntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
  transient_= Teuchos::getIntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
  pspg_     = Teuchos::getIntegralValue<INPAR::FLUID::PSPG>(stablist,"PSPG");
  supg_     = Teuchos::getIntegralValue<INPAR::FLUID::SUPG>(stablist,"SUPG");
  vstab_    = Teuchos::getIntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
  cstab_    = Teuchos::getIntegralValue<INPAR::FLUID::CStab>(stablist,"CSTAB");
  cross_    = Teuchos::getIntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
  reynolds_ = Teuchos::getIntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");

  //-------------------------------
  // get tau definition
  //-------------------------------

  whichtau_ =  Teuchos::getIntegralValue<INPAR::FLUID::TauType>(stablist,"DEFINITION_TAU");
  // check if tau can be handled
  if (not(whichtau_ == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall or
      INPAR::FLUID::tautype_bazilevs or
      INPAR::FLUID::tautype_franca_barrenechea_valentin_codina))
    dserror("Definition of Tau cannot handled by the element");

  // TODO: Adapt test cases, that it is necessary to use stationary tau definition in input file
  if (is_stationary_ == true)
    whichtau_ = INPAR::FLUID::tautype_stationary;

  // set flags for potential evaluation of tau and material law at int. point
  // default value: evaluation at element center
  const string tauloc = stablist.get<string>("EVALUATION_TAU");
  if (tauloc == "integration_point") tau_gp_ = true;
  else                               tau_gp_ = false;
  const string matloc = stablist.get<string>("EVALUATION_MAT");
  if (matloc == "integration_point") mat_gp_ = true;
  else                               mat_gp_ = false;

  // flag for higher order elements
  //bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  // TODO: stationary and inconsistent
  if (stablist.get<string>("STABTYPE") == "inconsistent") is_higher_order_ele_ = false;

  // get flag for fine-scale subgrid-viscosity approach
  {
    const string fssgvdef = params.get<string>("fs subgrid viscosity","No");

    if (fssgvdef == "Smagorinsky_all")        fssgv_ = Fluid3::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv_ = Fluid3::smagorinsky_small;
  }

  // the returned variable is saved to a const expression
  return timefac;
} // Fluid3Impl::SetSolutionParameter

/*----------------------------------------------------------------------*
 |  compute turbulence parameters                            ehrl 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetTurbulenceParams(
                               Fluid3*                    ele,
                               ParameterList&             turbmodelparams,
                               double&                    Cs,
                               double&                    Cs_delta_sq,
                               double&                    l_tau,
                               int&                       nlayer)
{
  if (is_stationary_ == true)
  {
    Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);
    Cs_delta_sq = 0.0;
    l_tau = 0.0;
    nlayer = 0;
  }
  else
  {
    //TODO: Check parameter for 2D and 3D case

    // get Smagorinsky model parameter for fine-scale subgrid viscosity
    // (Since either all-scale Smagorinsky model (i.e., classical LES model
    // as will be inititalized below) or fine-scale Smagorinsky model is
    // used (and never both), the same input parameter can be exploited.)
    if (fssgv_ != Fluid3::no_fssgv) Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

    if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
    {
      string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

      // --------------------------------------------------
      // standard constant coefficient Smagorinsky model
      if (physical_turbulence_model == "Smagorinsky")
      {
        // the classic Smagorinsky model only requires one constant parameter
        turb_mod_action_ = Fluid3::smagorinsky;
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
        turb_mod_action_ = Fluid3::smagorinsky_with_wall_damping;

        // get parameters of model
        Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
        l_tau           = turbmodelparams.get<double>("CHANNEL_L_TAU");

        // this will be the y-coordinate of a point in the element interior
        // we will determine the element layer in which he is contained to
        // be able to do the output of visceff etc.
        double center = 0;

        DRT::Node** nodes = ele->Nodes();
        for(int inode=0;inode<nen_;inode++)
        {
          center+=nodes[inode]->X()[1];
        }
        center/=nen_;

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
        turb_mod_action_ = Fluid3::dynamic_smagorinsky;

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
          for(int inode=0;inode<nen_;inode++)
          {
            center+=nodes[inode]->X()[1];
          }
          center/=nen_;

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
          Cs_delta_sq = ele->CsDeltaSq();
        }
      }
      else
      {
        dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
      }
    }
  } // end else(is_stationary_)
  return;
} // Fluid3Impl::GetturbulenceParams

/*----------------------------------------------------------------------*
 |  calculation of (all-scale) subgrid viscosity               vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const double                            vol,
  double&                                 Cs,
  double&                                 Cs_delta_sq,
  double&                                 l_tau
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);
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
  // -> either provided by dynamic modeling procedure and stored in Cs_delta_sq
  // -> or computed based on fixed Smagorinsky constant Cs:
  //             Cs = 0.17   (Lilly --- Determined from filter
  //                          analysis of Kolmogorov spectrum of
  //                          isotropic turbulence)
  //             0.1 < Cs < 0.24 (depending on the flow)
  //

  // compute (all-scale) rate of strain
  double rateofstrain = -1.0e30;
  // TODO: calculation of vderxy_ is realized two times!!
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  if (turb_mod_action_ == Fluid3::dynamic_smagorinsky)
  {
    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;

    // for evaluation of statistics: remember the 'real' Cs
    Cs = sqrt(Cs_delta_sq)/pow((vol),(1.0/3.0));
  }
  else if (turb_mod_action_ == Fluid3::dynamic_smagorinsky and nsd_==2)
    dserror("dynamic smagorinsky is not implemented in 2D");
  else
  {
    if (turb_mod_action_ == Fluid3::smagorinsky_with_wall_damping)
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

      LINALG::Matrix<nsd_,1> centernodecoord;
      centernodecoord.Multiply(xyze_,funct_);

      if (centernodecoord(1,0)>0) y_plus=(1.0-centernodecoord(1,0))/l_tau;
      else                        y_plus=(1.0+centernodecoord(1,0))/l_tau;

      //   lmix *= (1.0-exp(-y_plus/A_plus));
      // multiply with van Driest damping function
      Cs *= (1.0-exp(-y_plus/A_plus));
    }
    // TODO: Is it good to check it each time??
    else if (turb_mod_action_ == Fluid3::smagorinsky_with_wall_damping and nsd_==2)
      dserror("smagorinsky with wall damping is not implemented in 2D");

    // get characteristic element length for Smagorinsky model for 2D and 3D
    // 3D: hk = V^1/3
    // 2D: hk = A^1/2
    const double hk = pow(vol,(1.0/dim));

    // mixing length set proportional to grid witdh: lmix = Cs * hk
    double lmix = Cs * hk;

    Cs_delta_sq = lmix * lmix;

    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of fine-scale subgrid viscosity                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcFineScaleSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const LINALG::Matrix<nsd_,nen_>&        fsevelaf,
  const double                            vol,
  double&                                 Cs
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  //     // get characteristic element length for Smagorinsky model for 2D and 3D
  // 3D: hk = V^1/3
  // 2D: hk = A^1/2
  const double hk = pow(vol,(1.0/dim));

  if (fssgv_ == Fluid3::smagorinsky_all)
  {
    //
    // ALL-SCALE SMAGORINSKY MODEL
    // ---------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          / h \           / h \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * rateofstrain;
  }
  else if (fssgv_ == Fluid3::smagorinsky_small)
  {
    //
    // FINE-SCALE SMAGORINSKY MODEL
    // ----------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          /    \          /   \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // fine-scale rate of strain
    double fsrateofstrain = -1.0e30;
    fsrateofstrain = GetStrainRate(fsevelaf,derxy_,fsvderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * fsrateofstrain;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter                     vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcStabParameter(
  const double               timefac,
  const double               vol
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // ---------------------------------------------------------------
  // computation of stabilization parameter tau
  // ---------------------------------------------------------------
  //a)  Franca, Barrenechea, Valentin, Wall
  //b)  Franca, Barrenechea, Valentin, Wall stationary
  //c)  Bazilevs
  //d)  Codina

  // TODO: different element length calculations for 2D and 3D flow is
  // just a temporary work around (not to change the results of the test cases)

  //----------------------------------------------------------------
  // definiton of element size 'strle' for tau_Mu
  //----------------------------------------------------------------

  double strle = 0.0;
// a) streamlength (based on velocity vector at element centre) -> default

    if(nsd_==3)
    {
      // normed velocity at element centre
      if (vel_norm>=1e-6) velino_.Update(1.0/vel_norm,velint_);
      else
      {
        velino_.Clear();
        velino_(0,0) = 1;
      }

      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_, velino_);
      const double val = tmp.Norm1();
      strle = 2.0/val;
      // const double strle = 2.0/val;
    }

//  b) volume-equival. diameter -> not default
      // warning: 3D formula
      // const double strle = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

//  c) cubic/square root of the element volume/area  -> not default
    else if(nsd_==2)
      strle = pow(vol,1/dim);
      //const double strle = pow(vol,1/dim);

    else dserror("elment length calculation is not implemented for 1D flow");


  //----------------------------------------------------------------
  // definiton of element size 'hk' for tau_Mp and tau_C
  //----------------------------------------------------------------

    double hk = 0.0;
// a) streamlength (based on velocity vector at element centre) -> default

      // normed velocity at element centre
 /*     if (vel_norm>=1e-6) velino_.Update(1.0/vel_norm,velint_);
      else
      {
        velino_.Clear();
        velino_(0,0) = 1;
      }

      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_, velino_);
      const double val = tmp.Norm1();
      const double hk = 2.0/val;*/

//  b) volume-equival. diameter -> not default
      // warning: 3D formula
      if (nsd_==3)
        hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
        //const double hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

//  c) cubic/square root of the element volume/area  -> not default
      else if (nsd_==2)
        hk = pow(vol,1/dim);
        //const double hk = pow(vol,1/dim);

      else dserror("elment length calculation is not implemented for 1D flow");

  switch (whichtau_)
  {
  case INPAR::FLUID::tautype_franca_barrenechea_valentin_wall:
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

    // viscous : reactive forces
    const double re01 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(strle));

    // convective : viscous forces
    const double re02 = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*densaf_*xi01+(4.0*timefac*visceff_/mk)*xi02);

    if (nsd_ == 2)
      //TODO: in 2D hk and strlng was defined in the same way
      tau_(1)= tau_(0);
    else if (nsd_ == 3)
    {
      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)

      // viscous : reactive forces
      const double re11 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

      // convective : viscous forces
      const double re12 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

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

      tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi11+(4.0*timefac*visceff_/mk)*xi12);
    }
    else
      dserror("stabilization 'franca_barrenechea_valentin_wall' is not implemented for 1D flow");

    /*------------------------------------------------------ compute tau_C ---*/
    // PhD thesis Wall (1999)
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
    tau_(2) = densaf_ * vel_norm * hk * 0.5 * xi_tau_c;
  }
  break; // end franca_barrenechea_valentin_wall

  case INPAR::FLUID::tautype_stationary:
  {
    // TODO: density dependency is still missing

    // compute tau_Mu
    /* convective : viscous forces */
    const double re_tau_mu = mk * vel_norm * strle / (2.0 * visceff_);
    const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
    tau_(0) = (DSQR(strle)*mk)/(4.0*visceff_*xi_tau_mu);

    // compute tau_Mp
    /* convective : viscous forces */
    const double re_tau_mp = mk  * vel_norm * hk / (2.0 * visceff_);
    const double xi_tau_mp = DMAX(re_tau_mp,1.0);
    tau_(1) = (DSQR(hk)*mk)/(4.0*visceff_*xi_tau_mp);

    // compute tau_C
    const double re_tau_c = mk * vel_norm * hk / (2.0 * visceff_);
    const double xi_tau_c = DMIN(re_tau_c, 1.0);
    tau_(2) = 0.5*vel_norm*hk*xi_tau_c;

    break;
  }

  case INPAR::FLUID::tautype_bazilevs:
  {
    /*

    tau_M: Bazilevs et al. (2007)
    e.g.: Variational multiscale residual-based turbulence modeling for
          large eddy simulation of incompressible flows
                                                                              1.0
                 +-                                                      -+ - ---
                 |        2                                               |   2.0
                 | 4.0*rho         n+1             n+1          2         |
          tau  = | -------  + rho*u     * G * rho*u     + C * mu  * G : G |
             M   |     2                  -                I        -   - |
                 |   dt                   -                         -   - |
                 +-                                                      -+

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
    /*                                 +----
               n+1             n+1     \         n+1              n+1
          rho*u     * G * rho*u     =   +   rho*u    * G   * rho*u
                      -                /         i     -ij        j
                      -               +----        -
                                        i,j
    */

    double G;
    double normG = 0;
    double Gnormu = 0;

    const double dens_sqr = densaf_*densaf_;
    for (int nn=0;nn<nsd_;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*velint_(nn);
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG+=G*G;
        Gnormu+=dens_sqr_velint_nn*G*velint_(rr);
      }
    }

    // definition of constant:
    // 12.0/m_k = 36.0 for linear elements and 144.0 for quadratic elements
    // (differently defined, e.g., in Akkerman et al. (2008))
    const double CI = 12.0/mk;

    tau_(0) = 1.0/(sqrt((4.0*dens_sqr)/(dt_*dt_)+Gnormu+CI*visceff_*visceff_*normG));
    tau_(1) = tau_(0);

    /*
      tau_C: Bazilevs et al. (2007), derived from fine-scale complement Shur
                                     operator of the pressure equation

                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */
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
    for (int rr=0;rr<nsd_;++rr)
    {
      g = xji_(rr,0);
      for(int mm=1;mm<nsd_;++mm)
      {
        g += xji_(rr,mm);
      }
      normgsq += g*g;
    }

    tau_(2) = 1.0/(tau_(0)*normgsq);
  }
  break;  // end Bazilev

  case INPAR::FLUID::tautype_franca_barrenechea_valentin_codina:
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
    const double re01 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(strle));

    /* convective : viscous forces */
    const double re02 = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*densaf_*xi01+(4.0*timefac*visceff_/mk)*xi02);

    if(nsd_==2)
      //TODO: in 2D hk and strlng was defined in the same way
      tau_(1) = tau_(0);

    else if(nsd_==3)
    {

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)

      /* viscous : reactive forces */
      const double re11 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

      /* convective : viscous forces */
      const double re12 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

      const double xi11 = DMAX(re11,1.0);
      const double xi12 = DMAX(re12,1.0);

      /*
                    xi1,xi2 ^
                            |      /
                            |     /
                            |    /
                          1 +---+
                            |dens
                            |
                            |
                            +--------------> re1,re2
                                1
      */
      tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi11+(4.0*timefac*visceff_/mk)*xi12);
    }
    else
      dserror("stabilization 'codina' is not implemented for 1D flow");

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff_)+DSQR(0.5*densaf_*vel_norm*hk));
  }
  break; // end Codina


  case INPAR::FLUID::tautype_bazilevs_wo_dt:
  {
    /*
    tau_M by Bazilevs et al. (2007) adapted for the time-dependent subgrid-scale
    approach by Codina. Note that the density dependency is included in the ODE
    for the unresolved scales

                                                        1.0
                 +-                                -+ - ---
                 |                                  |   2.0
                 |  n+1      n+1          2         |
          tau  = | u     * Gu     + C * nu  * G : G |
             M   |         -         I        -   - |
                 |         -                  -   - |
                 +-                                -+

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
    /*                        +----
          n+af        n+af     \     n+af          n+af
         u     * G * u      =   +   u     * G   * u
                 -             /     i      -ij    j
                 -            +----         -
                               i,j
    */

    double G;
    double normG = 0;
    double Gnormu = 0;

    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG+=G*G;
        Gnormu+=velint_(nn)*G*velint_(rr);
      }
    }

    const double dens_sqr = densaf_*densaf_;

    // definition of constant:
    // 12.0/m_k = 36.0 for linear elements and 144.0 for quadratic elements
    // (differently defined, e.g., in Akkerman et al. (2008))
    const double CI = 12.0/mk;

    tau_(0) = 1.0/(sqrt(Gnormu+CI*(visceff_*visceff_/dens_sqr)*normG));
    tau_(1) = tau_(0);

    /*
      tau_C: Bazilevs et al. (2007), derived from fine-scale complement Shur
                                     operator of the pressure equation

                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */
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
    for (int rr=0;rr<nsd_;++rr)
    {
      g = xji_(rr,0);
      for(int mm=1;mm<nsd_;++mm)
      {
        g += xji_(rr,mm);
      }
      normgsq += g*g;
    }

    tau_(2) = 1.0/(tau_(0)*normgsq);
  }
  break;  // end Bazilevs for tds

  default: dserror("unknown definition of tau\n");
  }  // end switch

  return;
}

#endif
#endif

