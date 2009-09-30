/*----------------------------------------------------------------------*/
/*!
\file fluid3_stationary.cpp

\brief Internal implementation of Fluid3 element (stationary formulation)

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_stationary.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_condition_utils.H"

#ifdef DEBUG
//#define PRINTDEBUG
#endif

#ifdef PRINTDEBUG
#include <string>
#include <sstream>
#include <cstring>
template <class T>
void writeArray(const T& mat, std::string name = "unnamed")
{
  int M = mat.M();
  int N = mat.N();
  if (mat.A() == NULL) {
    M = N = 0;
  }
  std::stringstream header;
  header << 'M' << name << ':' << M << 'x' << N << ':';
  unsigned int s = header.str().size() + M*N*sizeof(double);
  std::cerr.write(reinterpret_cast<const char*>(&s),sizeof(unsigned int));
  std::cerr << header.str();
  if (M*N)
    std::cerr.write(reinterpret_cast<const char*>(mat.A()),M*N*sizeof(double));
  if (not std::cerr.good()) {
    std::cout << "Error: "<< std::cerr.good() << std::cerr.eof() << std::cerr.fail() << std::cerr.bad() << '\n';
    std::cerr.clear();
  }
}

extern void writeComment(const std::string v);

#endif // PB


DRT::ELEMENTS::Fluid3StationaryImplInterface* DRT::ELEMENTS::Fluid3StationaryImplInterface::Impl(Fluid3* f3)
{
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Fluid3StationaryImpl<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new Fluid3StationaryImpl<DRT::Element::hex8>();
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static Fluid3StationaryImpl<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new Fluid3StationaryImpl<DRT::Element::hex20>();
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static Fluid3StationaryImpl<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new Fluid3StationaryImpl<DRT::Element::hex27>();
    return fh27;
  }
  case DRT::Element::tet4:
  {
    static Fluid3StationaryImpl<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new Fluid3StationaryImpl<DRT::Element::tet4>();
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static Fluid3StationaryImpl<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new Fluid3StationaryImpl<DRT::Element::tet10>();
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static Fluid3StationaryImpl<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new Fluid3StationaryImpl<DRT::Element::wedge6>();
    return fw6;
  }
  case DRT::Element::wedge15:
  {
    static Fluid3StationaryImpl<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new Fluid3StationaryImpl<DRT::Element::wedge15>();
    return fw15;
  }
  case DRT::Element::pyramid5:
  {
    static Fluid3StationaryImpl<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new Fluid3StationaryImpl<DRT::Element::pyramid5>();
    return fp5;
  }

  default:
    dserror("shape %d (%d nodes) not supported", f3->Shape(), f3->NumNode());
  }
  return NULL;
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3StationaryImpl<distype>::Fluid3StationaryImpl()
  : Fluid3StationaryImplInterface(),
    vart_(),
    xyze_(),
    edeadng_(),
    funct_(),
    deriv_(),
    deriv2_(),
    xjm_(),
    xji_(),
    vderxy_(),
    fsvderxy_(),
    derxy_(),
    derxy2_(),
    bodyforce_(),
    velino_(),
    velint_(),
    fsvelint_(),
    gradp_(),
    tau_(),
    viscs2_(),
    conv_c_(),
    vdiv_(),
    rhsmom_(),
    conv_old_(),
    visc_old_(),
    res_old_(),
    conv_resM_(),
    xder2_()
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3StationaryImpl<distype>::Evaluate(
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

  // the number of nodes
  const int numnode = iel;
  // construct views
  LINALG::Matrix<4*iel,4*iel> elemat1(elemat1_epetra.A(),true);
  //LINALG::Matrix<4*iel,4*iel> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<4*iel,     1> elevec1(elevec1_epetra.A(),true);

  // need current velocity/pressure vector
  RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velaf");
  if (velnp==null) dserror("Cannot get state vector 'velaf'");

  // extract local values from the global vector
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  if (ele->is_ale_) dserror("No ALE support within stationary fluid solver.");

  // split velocity and pressure and set density
  LINALG::Matrix<numnode,1> eprenp;
  LINALG::Matrix<3,numnode> evelnp;

  for (int i=0;i<numnode;++i)
  {
    // split velocity and pressure, insert into element arrays
    evelnp(0,i) = myvelnp[0+(i*4)];
    evelnp(1,i) = myvelnp[1+(i*4)];
    evelnp(2,i) = myvelnp[2+(i*4)];

    eprenp(i) = myvelnp[3+(i*4)];
  }

  // get control parameter
  const double pseudotime = params.get<double>("total time",-1.0);
  if (pseudotime < 0.0)
    dserror("no value for total (pseudo-)time in the parameter list");

  // ---------------------------------------------------------------------
  // get control parameters for linearization, low-Mach-number solver
  // and form of convective term
  //----------------------------------------------------------------------
  string newtonstr   = params.get<string>("Linearisation");
  string convformstr = params.get<string>("form of convective term");
  bool newton = false;
  bool conservative = false;
  if(newtonstr=="Newton")          newton       = true;
  if(convformstr =="conservative") conservative = true;

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid3::StabilisationAction pspg     = ele->ConvertStringToStabAction(stablist.get<string>("PSPG"));
  Fluid3::StabilisationAction supg     = ele->ConvertStringToStabAction(stablist.get<string>("SUPG"));
  Fluid3::StabilisationAction vstab    = ele->ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  Fluid3::StabilisationAction cstab    = ele->ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  Fluid3::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  Fluid3::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // get flag for fine-scale subgrid-viscosity approach
  Fluid3::FineSubgridVisc fssgv = Fluid3::no_fssgv;
  {
    const string fssgvdef = params.get<string>("fs subgrid viscosity","No");

    if (fssgvdef == "artificial_all")         fssgv = Fluid3::artificial_all;
    else if (fssgvdef == "artificial_small")  fssgv = Fluid3::artificial_small;
    else if (fssgvdef == "Smagorinsky_all")   fssgv = Fluid3::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv = Fluid3::smagorinsky_small;
  }

  // flag for higher order elements
  bool higher_order_ele = ele->isHigherOrderElement(ele->Shape());

  // get fine-scale velocity
  RCP<const Epetra_Vector> fsvelnp;
  LINALG::Matrix<3,numnode> fsevelnp;

  if (fssgv != Fluid3::no_fssgv)
  {
    fsvelnp = discretization.GetState("fsvelnp");
    if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
    vector<double> myfsvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

    // get fine-scale velocity and insert into element arrays
    for (int i=0;i<numnode;++i)
    {
      fsevelnp(0,i) = myfsvelnp[0+(i*4)];
      fsevelnp(1,i) = myfsvelnp[1+(i*4)];
      fsevelnp(2,i) = myfsvelnp[2+(i*4)];
    }
  }
  else
  {
    for (int i=0;i<numnode;++i)
    {
      fsevelnp(0,i) = 0.0;
      fsevelnp(1,i) = 0.0;
      fsevelnp(2,i) = 0.0;
    }
  }

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");
  const double Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

  // calculate element coefficient matrix and rhs
  Sysmat(ele,
         evelnp,
         fsevelnp,
         eprenp,
         elemat1,
         elevec1,
         mat,
         pseudotime,
         newton,
         conservative,
         higher_order_ele,
         fssgv,
         pspg,
         supg,
         vstab,
         cstab,
         cross,
         reynolds,
         Cs);

  // This is a very poor way to transport the density to the
  // outside world. Is there a better one?
  /*double dens = 0.0;
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

  params.set("density", dens);*/
#ifdef PRINTDEBUG
  writeArray(elemat1,"elemat1");
  writeArray(elevec1,"elevec1");
#endif

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                  gjb 11/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3StationaryImpl<distype>::Sysmat(
  Fluid3*                                          ele,
  const LINALG::Matrix<3,iel>&                     evelnp,
  const LINALG::Matrix<3,iel>&                     fsevelnp,
  const LINALG::Matrix<iel,1>&                     eprenp,
  LINALG::Matrix<4*iel,4*iel>&                     estif,
  LINALG::Matrix<4*iel,1>&                         eforce,
  Teuchos::RCP<const MAT::Material>                material,
  double                                           pseudotime,
  const bool                                       newton,
  const bool                                       conservative,
  const bool                                       higher_order_ele,
  const enum Fluid3::FineSubgridVisc               fssgv,
  const enum Fluid3::StabilisationAction           pspg,
  const enum Fluid3::StabilisationAction           supg,
  const enum Fluid3::StabilisationAction           vstab,
  const enum Fluid3::StabilisationAction           cstab,
  const enum Fluid3::StabilisationAction           cross,
  const enum Fluid3::StabilisationAction           reynolds,
  const double                                     Cs
  )
{

  // get node coordinates and number of elements per node
  DRT::Node** const nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }

  // dead load in element nodes
  BodyForce(ele,pseudotime);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
  const double visc = actmat->Viscosity();

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  CalTauStationary(ele,evelnp,fsevelnp,visc,fssgv,Cs);

  // in case of viscous stabilisation decide whether to use GLS or usfemM
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
    //xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
#ifdef PRINTDEBUG
    writeArray(xjm_,"xjm");
    writeArray(xji_,"xji");
#endif
    if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

    const double fac = intpoints.qwgt[iquad]*det;

    // compute global derivates
    //derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);
    derxy_.Multiply(xji_,deriv_);

    //--------------------------------------------------------------
    //             compute second global derivative
    //--------------------------------------------------------------
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
    }
    else derxy2_.Clear();

    // get velocity at integration point
    velint_.Multiply(evelnp,funct_);

    // get velocity (np,i) derivatives at integration point
    vderxy_.MultiplyNT(evelnp,derxy_);

    // get fine-scale velocity (np,i) derivatives at integration point
    if (fssgv != Fluid3::no_fssgv) fsvderxy_.MultiplyNT(fsevelnp,derxy_);
    else                           fsvderxy_ = 0.;
#ifdef PRINTDEBUG
    writeArray(velint_,"velint");
    writeArray(vderxy_,"vderxy");
    writeArray(fsvderxy_,"fsvderxy");
#endif

    // get pressure gradients
    gradp_.Multiply(derxy_,eprenp);

    // get pressure at integration point
    const double press = funct_.Dot(eprenp);

    // get bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_,funct_);
#ifdef PRINTDEBUG
      writeArray(gradp_,"gradp");
      writeArray(bodyforce_,"bodyforce");
#endif

    //--------------------------------------------------------------
    // perform integration for entire matrix and rhs
    //--------------------------------------------------------------

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;

    /*------------------------- evaluate rhs vector at integration point ---*/
    //   rhsmom_ = histvec_(i) + bodyforce_(i);
    // histvec is always zero in stationary case (!):
    rhsmom_ = bodyforce_;  //copy
#ifdef PRINTDEBUG
      writeArray(rhsmom_,"rhsmom");
#endif

    /*----------------- get numerical representation of single operators ---*/
    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    //conv_c_ = blitz::sum(derxy_(j,i)*velint_(j), j);
    conv_c_.MultiplyTN(derxy_,velint_);

    /* Convective term  u_old * grad u_old: */
    //conv_old_ = blitz::sum(vderxy_(i, j)*velint_(j), j);
    conv_old_.Multiply(vderxy_,velint_);

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

      double prefac = 1.0;

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

      for (int i=1; i<iel; ++i)
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

    /* velocity divergence: */
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2);

    /*--------------------------------- now build single stiffness terms ---*/
    // evaluate residual once for all stabilisation right hand sides
    res_old_(0) = -rhsmom_(0)+(conv_old_(0)+gradp_(0)-2*visc*visc_old_(0));
    res_old_(1) = -rhsmom_(1)+(conv_old_(1)+gradp_(1)-2*visc*visc_old_(1));
    res_old_(2) = -rhsmom_(2)+(conv_old_(2)+gradp_(2)-2*visc*visc_old_(2));

#ifdef PRINTDEBUG
    writeArray(rhsmom_,"rhsmom");
    writeArray(conv_old_,"conv_old_");
    writeArray(conv_c_,"conv_c");
    writeArray(gradp_,"gradp_");
    writeArray(viscs2_,"viscs2");
    writeArray(visc_old_,"visc_old_");
    writeArray(res_old_,"res_old");
#endif

    /*
      This is the operator

                  /                      \
                 | (rho*resM)    o nabla |
                  \         (i)          /

                  required for (lhs) cross- and (rhs) Reynolds-stress calculation

    */

    if(cross    == Fluid3::cross_stress_stab ||
       reynolds == Fluid3::reynolds_stress_stab_only_rhs)
      conv_resM_.MultiplyTN(derxy_,res_old_);
#ifdef PRINTDEBUG
      writeArray(conv_resM_,"conv_resM");
#endif

    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      // computation of convection (convective and reactive part)
      // for conservative form including right-hand-side contribution
      if (conservative)
      {
        for (int ui=0; ui<iel; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v0 = fac*funct_(ui)*velint_(0);
          const double v1 = fac*funct_(ui)*velint_(1);
          const double v2 = fac*funct_(ui)*velint_(2);

          for (int vi=0; vi<iel; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* convection, convective part */
            /*

            /                                   \
            |  /              n+1  \             |
            | | Du (x) (rho*u)     |  , nabla v  |
            |  \              (i)  /             |
            \                                   /

            */
            estif(fvi  , fui  ) -= v0*derxy_(0, vi);
            estif(fvi  , fuip ) -= v0*derxy_(1, vi);
            estif(fvi  , fuipp) -= v0*derxy_(2, vi);
            estif(fvip , fui  ) -= v1*derxy_(0, vi);
            estif(fvip , fuip ) -= v1*derxy_(1, vi);
            estif(fvip , fuipp) -= v1*derxy_(2, vi);
            estif(fvipp, fui  ) -= v2*derxy_(0, vi);
            estif(fvipp, fuip ) -= v2*derxy_(1, vi);
            estif(fvipp, fuipp) -= v2*derxy_(2, vi);
          }
        }

        if (newton)
        {
          for (int ui=0; ui<iel; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = fac*funct_(ui);

            for (int vi=0; vi<iel; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              /* convection, reactive part */
              /*

              /                                     \
              |  / n+1                \             |
              | | u    (x) D(rho*u)   |  , nabla v  |
              |  \ (i)               /              |
              \                                    /

              */
              double v2 = v*conv_c_(vi) ;
              estif(fvi,   fui  ) -= v2;
              estif(fvip,  fuip ) -= v2;
              estif(fvipp, fuipp) -= v2;
            }
          }
        }

        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* convection */
          double v = fac*conv_c_(vi);
          eforce(fvi  ) += v*velint_(0) ;
          eforce(fvip ) += v*velint_(1) ;
          eforce(fvipp) += v*velint_(2) ;
        }
      }
      // computation of convection (convective and reactive part)
      // for convective form including right-hand-side contribution
      else
      {
        for (int ui=0; ui<iel; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v = fac*conv_c_(ui);
          for (int vi=0; vi<iel; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* convection, convective part */
            /*

            /                               \
            |  /       n+1        \         |
            | | (rho*u)   o nabla | Du , v  |
            |  \      (i)        /          |
            \                              /

            */
            double v2 = v*funct_(vi) ;
            estif(fvi  , fui  ) += v2;
            estif(fvip , fuip ) += v2;
            estif(fvipp, fuipp) += v2;
          }
        }

        if (newton)
        {
          for (int vi=0; vi<iel; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v = fac*funct_(vi);
            for (int ui=0; ui<iel; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v2 = v*funct_(ui);

              /*  convection, reactive part

              /                                 \
              |  /                \   n+1       |
              | | D(rho*u) o nabla | u     , v  |
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

        for (int vi=0; vi<iel; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* convection */
          double v = -fac*funct_(vi);
          eforce(fvi  ) += v*conv_old_(0) ;
          eforce(fvip ) += v*conv_old_(1) ;
          eforce(fvipp) += v*conv_old_(2) ;
        }
      }

      for (int ui=0; ui<iel; ++ui)
      {
        for (int vi=0; vi<iel; ++vi)
        {

          /* viscosity term */
          /*

                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(vi*4,     ui*4    ) += visc*fac*(2.0*derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)
                                                 +
                                                 derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4,     ui*4 + 1) += visc*fac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(vi*4,     ui*4 + 2) += visc*fac*derxy_(0, ui)*derxy_(2, vi) ;
          estif(vi*4 + 1, ui*4    ) += visc*fac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*4 + 1, ui*4 + 1) += visc*fac*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 2.0*derxy_(1, ui)*derxy_(1, vi)
                                                 +
                                                 derxy_(2, ui)*derxy_(2, vi)) ;
          estif(vi*4 + 1, ui*4 + 2) += visc*fac*derxy_(1, ui)*derxy_(2, vi) ;
          estif(vi*4 + 2, ui*4    ) += visc*fac*derxy_(0, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 1) += visc*fac*derxy_(1, vi)*derxy_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 2) += visc*fac*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)
                                                 +
                                                 2.0*derxy_(2, ui)*derxy_(2, vi)) ;

        }
      }

      for (int ui=0; ui<iel; ++ui)
      {
        double v = -fac*funct_(ui);
        for (int vi=0; vi<iel; ++vi)
        {
          /* pressure term */
          /*

          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */

          estif(vi*4,     ui*4 + 3) += v*derxy_(0, vi) ;
          estif(vi*4 + 1, ui*4 + 3) += v*derxy_(1, vi) ;
          estif(vi*4 + 2, ui*4 + 3) += v*derxy_(2, vi) ;

        }
      }

      for (int vi=0; vi<iel; ++vi)
      {
        double v = fac*funct_(vi);
        for (int ui=0; ui<iel; ++ui)
        {

          /* divergence */
          /*
            /                              \
            |                              |
            | nabla o D(rho*u)  , (q/rho)  |
            |                              |
            \                              /
          */
          estif(vi*4 + 3, ui*4    ) += v*derxy_(0, ui) ;
          estif(vi*4 + 3, ui*4 + 1) += v*derxy_(1, ui) ;
          estif(vi*4 + 3, ui*4 + 2) += v*derxy_(2, ui) ;
        }
      }

      for (int vi=0; vi<iel; ++vi)
      {
        /* pressure */
        double v = press*fac;
        eforce(vi*4    ) += v*derxy_(0, vi) ;
        eforce(vi*4 + 1) += v*derxy_(1, vi) ;
        eforce(vi*4 + 2) += v*derxy_(2, vi) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        /* viscosity */
        eforce(vi*4    ) -= visc*fac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                                       +
                                       derxy_(1, vi)*vderxy_(0, 1)
                                       +
                                       derxy_(1, vi)*vderxy_(1, 0)
                                       +
                                       derxy_(2, vi)*vderxy_(0, 2)
                                       +
                                       derxy_(2, vi)*vderxy_(2, 0)) ;
        eforce(vi*4 + 1) -= visc*fac*(derxy_(0, vi)*vderxy_(0, 1)
                                      +
                                      derxy_(0, vi)*vderxy_(1, 0)
                                      +
                                      2.0*derxy_(1, vi)*vderxy_(1, 1)
                                      +
                                      derxy_(2, vi)*vderxy_(1, 2)
                                      +
                                      derxy_(2, vi)*vderxy_(2, 1)) ;
        eforce(vi*4 + 2) -= visc*fac*(derxy_(0, vi)*vderxy_(0, 2)
                                      +
                                      derxy_(0, vi)*vderxy_(2, 0)
                                      +
                                      derxy_(1, vi)*vderxy_(1, 2)
                                      +
                                      derxy_(1, vi)*vderxy_(2, 1)
                                      +
                                      2.0*derxy_(2, vi)*vderxy_(2, 2)) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        // source term of the right hand side
        double v = fac*funct_(vi);
        eforce(vi*4    ) += v*rhsmom_(0) ;
        eforce(vi*4 + 1) += v*rhsmom_(1) ;
        eforce(vi*4 + 2) += v*rhsmom_(2) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        // continuity equation
        eforce(vi*4 + 3) -= fac*funct_(vi)*vdiv_ ;
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg == Fluid3::pstab_use_pspg)
      {
        for (int ui=0; ui<iel; ++ui)
        {
          double v = tau_Mp*conv_c_(ui) ;
          for (int vi=0; vi<iel; ++vi)
          {
            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | | (rho*u)   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */
            estif(vi*4 + 3, ui*4    ) += v*derxy_(0, vi) ;
            estif(vi*4 + 3, ui*4 + 1) += v*derxy_(1, vi) ;
            estif(vi*4 + 3, ui*4 + 2) += v*derxy_(2, vi) ;
          }
        }

        if (higher_order_ele)
        {
          double v = -2.0*visc*tau_Mp;
          for (int ui=0; ui<iel; ++ui)
          {
            for (int vi=0; vi<iel; ++vi)
            {

              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(vi*4 + 3, ui*4    ) += v*(derxy_(0, vi)*viscs2_(0, ui)
                                              +
                                              derxy_(1, vi)*viscs2_(1, ui)
                                              +
                                              derxy_(2, vi)*viscs2_(2, ui)) ;
              estif(vi*4 + 3, ui*4 + 1) += v*(derxy_(0, vi)*viscs2_(1, ui)
                                              +
                                              derxy_(1, vi)*viscs2_(4, ui)
                                              +
                                              derxy_(2, vi)*viscs2_(5, ui)) ;
              estif(vi*4 + 3, ui*4 + 2) += v*(derxy_(0, vi)*viscs2_(2, ui)
                                              +
                                              derxy_(1, vi)*viscs2_(5, ui)
                                              +
                                              derxy_(2, vi)*viscs2_(8, ui)) ;
            }
          }
        }

        for (int ui=0; ui<iel; ++ui)
        {
          for (int vi=0; vi<iel; ++vi)
          {
            /* pressure stabilisation: pressure( L_pres_p) */
            /*
              /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
              \                    /
            */
            estif(vi*4 + 3, ui*4 + 3) += tau_Mp*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)
                                                 +
                                                 derxy_(2, ui)*derxy_(2, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<iel; ++ui)
          {
            double v = tau_Mp*funct_(ui);
            for (int vi=0; vi<iel; ++vi)
            {
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | | D(rho*u) o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

              */
              estif(vi*4 + 3, ui*4    ) += v*(derxy_(0, vi)*vderxy_(0, 0)
                                              +
                                              derxy_(1, vi)*vderxy_(1, 0)
                                              +
                                              derxy_(2, vi)*vderxy_(2, 0)) ;
              estif(vi*4 + 3, ui*4 + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                                              +
                                              derxy_(1, vi)*vderxy_(1, 1)
                                              +
                                              derxy_(2, vi)*vderxy_(2, 1)) ;
              estif(vi*4 + 3, ui*4 + 2) += v*(derxy_(0, vi)*vderxy_(0, 2)
                                              +
                                              derxy_(1, vi)*vderxy_(1, 2)
                                              +
                                              derxy_(2, vi)*vderxy_(2, 2)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<iel; ++vi)
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
        for (int ui=0; ui<iel; ++ui)
        {
          double v = tau_M*conv_c_(ui);
          for (int vi=0; vi<iel; ++vi)
          {
            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                         \
            |    /       n+1        \        /       n+1        \     |
            |   | (rho*u)    o nabla | Du , | (rho*u)    o nabla | v  |
            |    \       (i)        /        \       (i)        /     |
            \                                                         /

            */
            estif(vi*4,     ui*4    ) += v*conv_c_(vi);
            estif(vi*4 + 1, ui*4 + 1) += v*conv_c_(vi);
            estif(vi*4 + 2, ui*4 + 2) += v*conv_c_(vi);
          }
        }

        for (int vi=0; vi<iel; ++vi)
        {
          double v = tau_M*conv_c_(vi);
          for (int ui=0; ui<iel; ++ui)
          {
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , | (rho*u)   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(vi*4,     ui*4 + 3) += v*derxy_(0, ui) ;
            estif(vi*4 + 1, ui*4 + 3) += v*derxy_(1, ui) ;
            estif(vi*4 + 2, ui*4 + 3) += v*derxy_(2, ui) ;
          }
        }

        if (higher_order_ele)
        {
          for (int vi=0; vi<iel; ++vi)
          {
            double v = -2.0*visc*tau_M*conv_c_(vi);
            for (int ui=0; ui<iel; ++ui)
            {
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, | (rho*u)    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
              */
              estif(vi*4,     ui*4    ) += v*viscs2_(0, ui) ;
              estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, ui) ;
              estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, ui) ;

              estif(vi*4,     ui*4 + 1) += v*viscs2_(1, ui) ;
              estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, ui) ;
              estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, ui) ;

              estif(vi*4,     ui*4 + 2) += v*viscs2_(2, ui) ;
              estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, ui) ;
              estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, ui) ;
            }
          }
        }

        if (newton)
        {
          {
            double v0 = velint_(0)*vderxy_(0, 0) + velint_(1)*vderxy_(0, 1) + velint_(2)*vderxy_(0, 2);
            double v1 = velint_(0)*vderxy_(1, 0) + velint_(1)*vderxy_(1, 1) + velint_(2)*vderxy_(1, 2);
            double v2 = velint_(0)*vderxy_(2, 0) + velint_(1)*vderxy_(2, 1) + velint_(2)*vderxy_(2, 2);

            for (int ui=0; ui<iel; ++ui)
            {
              double v = tau_M*funct_(ui);
              for (int vi=0; vi<iel; ++vi)
              {
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
                estif(vi*4,     ui*4    ) += (conv_c_(vi)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(vi*4 + 1, ui*4    ) += (conv_c_(vi)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;
                estif(vi*4 + 2, ui*4    ) += (conv_c_(vi)*vderxy_(2, 0) + v2*derxy_(0, vi))*v;

                estif(vi*4,     ui*4 + 1) += (conv_c_(vi)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(vi*4 + 1, ui*4 + 1) += (conv_c_(vi)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
                estif(vi*4 + 2, ui*4 + 1) += (conv_c_(vi)*vderxy_(2, 1) + v2*derxy_(1, vi))*v;

                estif(vi*4,     ui*4 + 2) += (conv_c_(vi)*vderxy_(0, 2) + v0*derxy_(2, vi))*v;
                estif(vi*4 + 1, ui*4 + 2) += (conv_c_(vi)*vderxy_(1, 2) + v1*derxy_(2, vi))*v;
                estif(vi*4 + 2, ui*4 + 2) += (conv_c_(vi)*vderxy_(2, 2) + v2*derxy_(2, vi))*v;
              }
            }
          }

          for (int ui=0; ui<iel; ++ui)
          {
            double v = tau_M*funct_(ui);
            for (int vi=0; vi<iel; ++vi)
            {
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , | D(rho*u) o nabla | v  |
                |         (i)    \                /     |
                \                                      /
              */
              estif(vi*4,     ui*4    ) += v*gradp_(0)*derxy_(0, vi) ;
              estif(vi*4 + 1, ui*4    ) += v*gradp_(1)*derxy_(0, vi) ;
              estif(vi*4 + 2, ui*4    ) += v*gradp_(2)*derxy_(0, vi) ;

              estif(vi*4,     ui*4 + 1) += v*gradp_(0)*derxy_(1, vi) ;
              estif(vi*4 + 1, ui*4 + 1) += v*gradp_(1)*derxy_(1, vi) ;
              estif(vi*4 + 2, ui*4 + 1) += v*gradp_(2)*derxy_(1, vi) ;

              estif(vi*4,     ui*4 + 2) += v*gradp_(0)*derxy_(2, vi) ;
              estif(vi*4 + 1, ui*4 + 2) += v*gradp_(1)*derxy_(2, vi) ;
              estif(vi*4 + 2, ui*4 + 2) += v*gradp_(2)*derxy_(2, vi) ;

            }
          }

          if (higher_order_ele)
          {
            for (int ui=0; ui<iel; ++ui)
            {
              double v = -2.0*visc*tau_M*funct_(ui);
              for (int vi=0; vi<iel; ++vi)
              {
                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, | D(rho*u) o nabla | v  |
                  |               \ (i) /    \                /     |
                  \                                                 /
                */
                estif(vi*4,     ui*4    ) += v*visc_old_(0)*derxy_(0, vi) ;
                estif(vi*4 + 1, ui*4    ) += v*visc_old_(1)*derxy_(0, vi) ;
                estif(vi*4 + 2, ui*4    ) += v*visc_old_(2)*derxy_(0, vi) ;

                estif(vi*4,     ui*4 + 1) += v*visc_old_(0)*derxy_(1, vi) ;
                estif(vi*4 + 1, ui*4 + 1) += v*visc_old_(1)*derxy_(1, vi) ;
                estif(vi*4 + 2, ui*4 + 1) += v*visc_old_(2)*derxy_(1, vi) ;

                estif(vi*4,     ui*4 + 2) += v*visc_old_(0)*derxy_(2, vi) ;
                estif(vi*4 + 1, ui*4 + 2) += v*visc_old_(1)*derxy_(2, vi) ;
                estif(vi*4 + 2, ui*4 + 2) += v*visc_old_(2)*derxy_(2, vi) ;
              }
            }
          }

          for (int ui=0; ui<iel; ++ui)
          {
            double v = -tau_M*funct_(ui);
            for (int vi=0; vi<iel; ++vi)
            {
              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                                     \
                |              /                \     |
                |  rhsint   , | D(rho*u) o nabla | v  |
                |              \                /     |
                \                                     /

              */
              estif(vi*4    , ui*4    ) += v*derxy_(0, vi)*rhsmom_(0) ;
              estif(vi*4 + 1, ui*4    ) += v*derxy_(0, vi)*rhsmom_(1) ;
              estif(vi*4 + 2, ui*4    ) += v*derxy_(0, vi)*rhsmom_(2) ;

              estif(vi*4    , ui*4 + 1) += v*derxy_(1, vi)*rhsmom_(0) ;
              estif(vi*4 + 1, ui*4 + 1) += v*derxy_(1, vi)*rhsmom_(1) ;
              estif(vi*4 + 2, ui*4 + 1) += v*derxy_(1, vi)*rhsmom_(2) ;

              estif(vi*4    , ui*4 + 2) += v*derxy_(2, vi)*rhsmom_(0) ;
              estif(vi*4 + 1, ui*4 + 2) += v*derxy_(2, vi)*rhsmom_(1) ;
              estif(vi*4 + 2, ui*4 + 2) += v*derxy_(2, vi)*rhsmom_(2) ;

            } // vi
          } // ui
        } // if newton

        /* right hand side */
        for (int vi=0; vi<iel; ++vi)
        {
          // supg stabilisation
          double v = -tau_M*conv_c_(vi);
          eforce(vi*4    ) += v*res_old_(0) ;
          eforce(vi*4 + 1) += v*res_old_(1) ;
          eforce(vi*4 + 2) += v*res_old_(2) ;
        }
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
            const double four_visc2_tauMp = vstabfac*4.0*visc*visc*tau_Mp;

            // viscous stabilization on left hand side

            for (int ui=0; ui<iel; ++ui)
            {
              double v = two_visc_tauMp*conv_c_(ui);
              for (int vi=0; vi<iel; ++vi)
              {
                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | | (rho*u)   o nabla | Du , div eps (v) |
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

            for (int ui=0; ui<iel; ++ui)
            {
              for (int vi=0; vi<iel; ++vi)
              {

                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(vi*4, ui*4 + 3    ) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                             +
                                                             derxy_(1, ui)*viscs2_(1, vi)
                                                             +
                                                             derxy_(2, ui)*viscs2_(2, vi)) ;
                estif(vi*4 + 1, ui*4 + 3) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                             +
                                                             derxy_(1, ui)*viscs2_(4, vi)
                                                             +
                                                             derxy_(2, ui)*viscs2_(5, vi)) ;
                estif(vi*4 + 2, ui*4 + 3) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(2, vi)
                                                             +
                                                             derxy_(1, ui)*viscs2_(5, vi)
                                                             +
                                                             derxy_(2, ui)*viscs2_(8, vi)) ;

              }
            }

            for (int ui=0; ui<iel; ++ui)
            {
              for (int vi=0; vi<iel; ++vi)
              {
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(vi*4,     ui*4    ) -= four_visc2_tauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(2,ui)*viscs2_(2,vi)) ;
                estif(vi*4 + 1, ui*4    ) -= four_visc2_tauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(4,vi)+viscs2_(2,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4    ) -= four_visc2_tauMp*(viscs2_(0,ui)*viscs2_(2,vi)+viscs2_(1,ui)*viscs2_(5,vi)+viscs2_(2,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 1) -= four_visc2_tauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(4,ui)+viscs2_(2,vi)*viscs2_(5,ui)) ;
                estif(vi*4 + 1, ui*4 + 1) -= four_visc2_tauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(4,ui)*viscs2_(4,vi)+viscs2_(5,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4 + 1) -= four_visc2_tauMp*(viscs2_(1,ui)*viscs2_(2,vi)+viscs2_(4,ui)*viscs2_(5,vi)+viscs2_(5,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 2) -= four_visc2_tauMp*(viscs2_(0,vi)*viscs2_(2,ui)+viscs2_(1,vi)*viscs2_(5,ui)+viscs2_(2,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 1, ui*4 + 2) -= four_visc2_tauMp*(viscs2_(1,vi)*viscs2_(2,ui)+viscs2_(4,vi)*viscs2_(5,ui)+viscs2_(5,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 2, ui*4 + 2) -= four_visc2_tauMp*(viscs2_(2,ui)*viscs2_(2,vi)+viscs2_(5,ui)*viscs2_(5,vi)+viscs2_(8,ui)*viscs2_(8,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<iel; ++ui)
              {
                double v = two_visc_tauMp*funct_(ui);
                for (int vi=0; vi<iel; ++vi)
                {
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | | D(rho*u) o nabla | u    , div eps (v) |
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

          for (int vi=0; vi<iel; ++vi)
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

      if(cstab == Fluid3::continuity_stab_yes)
      {
        const double tau_C_divunp=tau_C*vdiv_;

        for (int ui=0; ui<iel; ++ui)
        {
          double v0 = tau_C*derxy_(0, ui);
          double v1 = tau_C*derxy_(1, ui);
          double v2 = tau_C*derxy_(2, ui);
          for (int vi=0; vi<iel; ++vi)
          {
            /* continuity stabilisation on left hand side */
            /*
              /                                      \
              |                                      |
              | nabla o D(rho*u)  , nabla o (rho*v)  |
              |                                      |
              \                                      /
            */
            estif(vi*4,     ui*4    ) += v0*derxy_(0, vi) ;
            estif(vi*4 + 1, ui*4    ) += v0*derxy_(1, vi) ;
            estif(vi*4 + 2, ui*4    ) += v0*derxy_(2, vi) ;

            estif(vi*4,     ui*4 + 1) += v1*derxy_(0, vi) ;
            estif(vi*4 + 1, ui*4 + 1) += v1*derxy_(1, vi) ;
            estif(vi*4 + 2, ui*4 + 1) += v1*derxy_(2, vi) ;

            estif(vi*4,     ui*4 + 2) += v2*derxy_(0, vi) ;
            estif(vi*4 + 1, ui*4 + 2) += v2*derxy_(1, vi) ;
            estif(vi*4 + 2, ui*4 + 2) += v2*derxy_(2, vi) ;
          }
        }

        for (int vi=0; vi<iel; ++vi)
        {
          /* continuity stabilisation on right hand side */
          eforce(vi*4    ) -= tau_C_divunp*derxy_(0, vi) ;
          eforce(vi*4 + 1) -= tau_C_divunp*derxy_(1, vi) ;
          eforce(vi*4 + 2) -= tau_C_divunp*derxy_(2, vi) ;

        }
      }

      //----------------------------------------------------------------------
      //     STABILIZATION, CROSS-STRESS PART (RESIDUAL-BASED VMM)

      if(cross == Fluid3::cross_stress_stab_only_rhs || cross == Fluid3::cross_stress_stab)
      {
        if(cross == Fluid3::cross_stress_stab)
        {
          for (int ui=0; ui<iel; ++ui)
          {
            double v = tau_M*conv_resM_(ui);
            for (int vi=0; vi<iel; ++vi)
            {
              /* cross-stress part on lhs */
              /*

                          /                               \
                         |  /                  \          |
                      -  | | (rho*resM) o nabla | Du , v  |
                         |  \                  /          |
                          \                               /
              */
              double v2 = v*funct_(vi);
              estif(vi*4    , ui*4    ) -= v2 ;
              estif(vi*4 + 1, ui*4 + 1) -= v2 ;
              estif(vi*4 + 2, ui*4 + 2) -= v2 ;
            }
          }
        } // end cross-stress part on left hand side

        for (int vi=0; vi<iel; ++vi)
        {
          /* cross-stress part on rhs */
          /*

                          /                                \
                         |  /                  \           |
                         | | (rho*resM) o nabla | u   , v  |
                         |  \                  /  (i)      |
                          \                               /
          */
          double v = tau_M*funct_(vi);
          eforce(vi*4    ) += v*(res_old_(0)*vderxy_(0,0) +
                                 res_old_(1)*vderxy_(0,1) +
                                 res_old_(2)*vderxy_(0,2));
          eforce(vi*4 + 1) += v*(res_old_(0)*vderxy_(1,0) +
                                 res_old_(1)*vderxy_(1,1) +
                                 res_old_(2)*vderxy_(1,2));
          eforce(vi*4 + 2) += v*(res_old_(0)*vderxy_(2,0) +
                                 res_old_(1)*vderxy_(2,1) +
                                 res_old_(2)*vderxy_(2,2));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if(reynolds == Fluid3::reynolds_stress_stab_only_rhs)
      {
        const double tauMtauM = tau_M*tau_M/fac;
        for (int vi=0; vi<iel; ++vi)
        {
          /* Reynolds-stress part on rhs */
          /*

                  /                                    \
                 |                                     |
                 |  resM   ,   (rho*resM) o nabla ) v  |
                 |                                     |
                  \                                    /
          */
          double v = tauMtauM*conv_resM_(vi);
          eforce(vi*4    ) += v*res_old_(0);
          eforce(vi*4 + 1) += v*res_old_(1);
          eforce(vi*4 + 2) += v*res_old_(2);
        }
      } // end Reynolds-stress part on right hand side

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv != Fluid3::no_fssgv)
      {
        for (int vi=0; vi<iel; ++vi)
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
  }
  return;
}



//
// calculate stabilization parameter
//
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3StationaryImpl<distype>::CalTauStationary(
  Fluid3*                             ele,
  const LINALG::Matrix<3,iel>&        evelnp,
  const LINALG::Matrix<3,iel>&        fsevelnp,
  const double                        visc,
  const enum Fluid3::FineSubgridVisc  fssgv,
  const double                        Cs
  )
{
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

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element type constant for tau
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
  //velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);
  velint_.Multiply(evelnp,funct_);

  // get Jacobian matrix and determinant
  //xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xji_.Invert(xjm_);
#ifdef PRINTDEBUG
  writeArray(velint_,"velint");
  writeArray(xjm_,"xjm");
  writeArray(xji_,"xji");
#endif

  // check for degenerated elements
  if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

  const double vol = wquad*det;

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // compute global derivates
  //derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);
  derxy_.Multiply(xji_,deriv_);

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // normed velocity at element centre
  if (vel_norm>=1e-6)
  {
    //velino_ = velint_/vel_norm;
    velino_.Update(1.0/vel_norm,velint_);
  }
  else
  {
    velino_ = 0.;
    velino_(0) = 1;
  }

  // get streamlength
  LINALG::Matrix<iel,1> tmp;
  tmp.MultiplyTN(derxy_,velino_);
  const double val = tmp.Norm1();
  const double strle = 2.0/val;
#ifdef PRINTDEBUG
  writeArray(velint_,"velint");
  writeArray(derxy_,"derxy");
  writeArray(velino_,"velino");
  writeArray(xji_,"xji");
  LINALG::Matrix<3,1> vvs;
  vvs(0) = vel_norm;
  vvs(1) = val;
  vvs(2) = strle;
  writeArray(vvs,"velnorm,val,strle");
#endif

  // calculate tau
  // stabilization parameters for stationary case

  // compute tau_Mu
  /* convective : viscous forces */
  const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);
  const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
  tau_(0) = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

  // compute tau_Mp
  /* convective : viscous forces */
  const double re_tau_mp = mk  * vel_norm * hk / (2.0 * visc);
  const double xi_tau_mp = DMAX(re_tau_mp,1.0);
  tau_(1) = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

  // compute tau_C
  const double re_tau_c = mk * vel_norm * hk / (2.0 * visc);
  const double xi_tau_c = DMIN(re_tau_c, 1.0);
  tau_(2) = 0.5*vel_norm*hk*xi_tau_c;

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv == Fluid3::artificial_all or fssgv == Fluid3::artificial_small)
  {
    double fsvel_norm = 0.0;
    if (fssgv == Fluid3::artificial_small)
    {
      // get fine-scale velocities at element center
      //fsvelint_ = blitz::sum(funct_(j)*fsevelnp(i,j),j);
      fsvelint_.Multiply(fsevelnp,funct_);

      // get fine-scale velocity norm
      //fsvel_norm = sqrt(blitz::sum(fsvelint_*fsvelint_));
      fsvel_norm = fsvelint_.Norm2();
    }
    // get all-scale velocity norm
    else fsvel_norm = vel_norm;
#ifdef PRINTDEBUG
  writeArray(fsvelint_,"fsvelint");
#endif

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re = mk * fsvel_norm * hk / visc;
    const double xi = DMAX(re,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(fsvel_norm))/(2.0*visc*xi);

  }
  else if (fssgv == Fluid3::smagorinsky_all or
           fssgv == Fluid3::smagorinsky_small)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                                      +-                                 -+ 1
    //                                  2   |          / h \           / h \    | -
    //    visc          =  (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    double rateofstrain = 0.0;
    {
      // get fine-scale or all-scale velocity derivatives at element center
      if (fssgv == Fluid3::smagorinsky_small)
        //fsvderxy_ = blitz::sum(derxy_(j,k)*fsevelnp(i,k),k);
        fsvderxy_.MultiplyNT(fsevelnp,derxy_);
      else  //fsvderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);
        fsvderxy_.MultiplyNT(evelnp,derxy_);
#ifdef PRINTDEBUG
  writeArray(fsvderxy_,"fsvderxy");
#endif

      double tmp;
      for(int rr=0;rr<3;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          tmp = (fsvderxy_(rr,mm) + fsvderxy_(mm,rr));
          rateofstrain += tmp*tmp;
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

    vart_ = Cs * Cs * hk * hk * rateofstrain;
#ifdef PRINTDEBUG
    LINALG::Matrix<2,1> rv;
    rv(0) = rateofstrain;
    rv(1) = vart_;
    writeArray(rv,"rateofstrain,vart");
#endif
  }
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3StationaryImpl<distype>::BodyForce(Fluid3*      ele,
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
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X(),time,NULL);
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

  return;
}


#endif
#endif
