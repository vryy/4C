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
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_condition_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2StationaryInterface* DRT::ELEMENTS::Fluid2StationaryInterface::Impl(DRT::ELEMENTS::Fluid2* f2)
{
  switch (f2->Shape())
  {
  case DRT::Element::quad4:
  {
    static Fluid2Stationary<DRT::Element::quad4>* fq4;
    if (fq4==NULL)
      fq4 = new Fluid2Stationary<DRT::Element::quad4>;
    return fq4;
  }
  case DRT::Element::quad8:
  {
    static Fluid2Stationary<DRT::Element::quad8>* fq8;
    if (fq8==NULL)
      fq8 = new Fluid2Stationary<DRT::Element::quad8>;
    return fq8;
  }
  case DRT::Element::quad9:
  {
    static Fluid2Stationary<DRT::Element::quad9>* fq9;
    if (fq9==NULL)
      fq9 = new Fluid2Stationary<DRT::Element::quad9>;
    return fq9;
  }
  case DRT::Element::tri3:
  {
    static Fluid2Stationary<DRT::Element::tri3>* ft3;
    if (ft3==NULL)
      ft3 = new Fluid2Stationary<DRT::Element::tri3>;
    return ft3;
  }
  case DRT::Element::tri6:
  {
    static Fluid2Stationary<DRT::Element::tri6>* ft6;
    if (ft6==NULL)
      ft6 = new Fluid2Stationary<DRT::Element::tri6>;
    return ft6;
  }
  default:
    dserror("shape %d (%d nodes) not supported", f2->Shape(), f2->NumNode());
  }
  return NULL;
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid2Stationary<distype>::Fluid2Stationary()
  : vart_(),
    xyze_(),
    edeadng_(),
    funct_(),
    densfunct_(),
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
    //velino_(2),
    velint_(),
    ndwvelint_(),
    fsvelint_(),
    gradp_(),
    tau_(),
    viscs2_(),
    conv_c_(),
    ndwconv_c_(),
    mdiv_(),
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
int DRT::ELEMENTS::Fluid2Stationary<distype>::Evaluate(
  Fluid2*                    ele,
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
  const int numnode = iel;

  LINALG::Matrix<3*iel,3*iel> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<3*iel,3*iel> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<3*iel,    1> elevec1(elevec1_epetra.A(),true);

  // need current velocity/pressure and velocity/density vector
  RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
  RefCountPtr<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  if (velnp==null || vedenp==null)
    dserror("Cannot get state vectors 'velnp' and/or 'vedenp'");

  // extract local values from the global vectors
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  vector<double> myvedenp(lm.size());
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);

  if (ele->is_ale_) dserror("No ALE support within stationary fluid solver.");

  // create objects for element arrays
  LINALG::Matrix<numnode, 1> eprenp;
  LINALG::Matrix<2, numnode> evelnp;
  LINALG::Matrix<numnode, 1> edensnp;

  for (int i=0;i<numnode;++i)
  {
    // split velocity and pressure, insert into element arrays
    evelnp(0,i) = myvelnp[0+(i*3)];
    evelnp(1,i) = myvelnp[1+(i*3)];

    eprenp(i) = myvelnp[2+(i*3)];

    // insert density vector into element array
    edensnp(i) = myvedenp[2+(i*3)];
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
  string lomastr     = params.get<string>("low-Mach-number solver");
  string convformstr = params.get<string>("form of convective term");
  bool newton = false;
  bool loma   = false;
  bool conservative = false;
  if(newtonstr=="Newton")          newton       = true;
  if(lomastr  =="Yes")             loma         = true;
  if(convformstr =="conservative") conservative = true;

  // --------------------------------------------------
  // set parameters for stabilisation
  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid2::StabilisationAction pspg     = ele->ConvertStringToStabAction(stablist.get<string>("PSPG"));
  Fluid2::StabilisationAction supg     = ele->ConvertStringToStabAction(stablist.get<string>("SUPG"));
  Fluid2::StabilisationAction vstab    = ele->ConvertStringToStabAction(stablist.get<string>("VSTAB"));
  Fluid2::StabilisationAction cstab    = ele->ConvertStringToStabAction(stablist.get<string>("CSTAB"));
  Fluid2::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  Fluid2::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  // flag for higher order elements
  bool higher_order_ele = ele->isHigherOrderElement(distype);

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  // get fine-scale velocity
  RCP<const Epetra_Vector> fsvelnp;
  LINALG::Matrix<2,numnode> fsevelnp;

  // get flag for fine-scale subgrid-viscosity approach
  Fluid2::FineSubgridVisc fssgv = Fluid2::no_fssgv;
  {
    const string fssgvdef = params.get<string>("fs subgrid viscosity","No");

    if (fssgvdef == "artificial_all")         fssgv = Fluid2::artificial_all;
    else if (fssgvdef == "artificial_small")  fssgv = Fluid2::artificial_small;
    else if (fssgvdef == "Smagorinsky_all")   fssgv = Fluid2::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv = Fluid2::smagorinsky_small;
  }

  if (fssgv != Fluid2::no_fssgv)
  {
    fsvelnp = discretization.GetState("fsvelnp");
    if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
    vector<double> myfsvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

    // get fine-scale velocity and insert into element arrays
    for (int i=0;i<numnode;++i)
    {
      fsevelnp(0,i) = myfsvelnp[0+(i*3)];
      fsevelnp(1,i) = myfsvelnp[1+(i*3)];
    }
  }
  else
  {
    for (int i=0;i<numnode;++i)
    {
      fsevelnp.Clear();
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
         edensnp,
         elemat1,
         elevec1,
         mat,
         pseudotime,
         newton,
         loma,
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
      MAT::NewtonianFluid* actmat = static_cast<MAT::NewtonianFluid*>(mat.get());
      dens = actmat->Density();
    }
    else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
    {
      MAT::CarreauYasuda* actmat = static_cast<MAT::CarreauYasuda*>(mat.get());
      dens = actmat->Density();
    }
    else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
    {
      MAT::ModPowerLaw* actmat = static_cast<MAT::ModPowerLaw*>(mat.get());
      dens = actmat->Density();
    }
    else
      dserror("no fluid material found");

    params.set("density", dens);*/

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                  vg 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Stationary<distype>::Sysmat(
  Fluid2*                                 ele,
  const LINALG::Matrix<2,iel>&            evelnp,
  const LINALG::Matrix<2,iel>&            fsevelnp,
  const LINALG::Matrix<iel,1>&            eprenp,
  const LINALG::Matrix<iel,1>&            edensnp,
  LINALG::Matrix<3*iel,3*iel>&            estif,
  LINALG::Matrix<3*iel,    1>&            eforce,
  Teuchos::RCP<const MAT::Material>       material,
  double                                  pseudotime,
  const bool                              newton,
  const bool                              loma,
  const bool                              conservative,
  const bool                              higher_order_ele,
  const enum Fluid2::FineSubgridVisc      fssgv,
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
  const int numnode = iel;

  // get node coordinates and number of elements per node
  DRT::Node** const nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
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
  CalTauStationary(ele,
                   evelnp,
                   fsevelnp,
                   edensnp,
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
    // inverse of jacobian
    xjm_.MultiplyNT(deriv_, xyze_);

    // The determinant is computed using Sarrus's rule:
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

    // compute global derivates
    derxy_.Multiply(xji_, deriv_);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<numnode; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);
      functdens_(inode) = funct_(inode)/edensnp(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    //--------------------------------------------------------------
    //             compute global second derivatives
    //--------------------------------------------------------------
    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
    }
    else derxy2_.Clear();

    // get momentum (n+g,i) at integration point
    velint_.Multiply(evelnp, densfunct_);

    // get velocity (np,i) derivatives at integration point
    vderxy_.MultiplyNT(evelnp, derxy_);

    // get momentum (np,i) derivatives at integration point
    mderxy_.MultiplyNT(evelnp, densderxy_);

    // get fine-scale velocity (np,i) derivatives at integration point
    if (fssgv != Fluid2::no_fssgv) fsvderxy_.MultiplyNT(fsevelnp, derxy_);
    else fsvderxy_.Clear();

    // get pressure gradients
    gradp_.Multiply(derxy_, eprenp);

    // get pressure at integration point
    double press = funct_.Dot(eprenp);

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_, densfunct_);

    // perform integration for entire matrix and rhs

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    // subgrid-viscosity factor
    const double vartfac = vart_*fac;

    /*------------------------- evaluate rhs vector at integration point ---*/
    // histvectors are always zero in stationary case (!):
    rhsmom_.Update(bodyforce_);

    /*----------------- get numerical representation of single operators ---*/
    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_, velint_);

    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_, velint_);

    // For conservative form of convective term, we also need "pure"
    // (i.e., non-density-weighted) velocity (including potential ALE
    // velocity) and convective part.
    if (conservative)
    {
      ndwvelint_.Multiply(evelnp,funct_);
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

    /* momentum and velocity divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    /*--------------------------------- now build single stiffness terms ---*/
    // evaluate residual once for all stabilisation right hand sides
    res_old_(0) = -rhsmom_(0)+(conv_old_(0)+gradp_(0)-2*visc*visc_old_(0));
    res_old_(1) = -rhsmom_(1)+(conv_old_(1)+gradp_(1)-2*visc*visc_old_(1));

    /*
      This is the operator

                  /                      \
                 | (rho*resM)    o nabla |
                  \         (i)          /

                  required for (lhs) cross- and (rhs) Reynolds-stress calculation

    */

    if(cross    == Fluid2::cross_stress_stab ||
       reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      conv_resM_.MultiplyTN(densderxy_, res_old_);

   {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      // computation of convection (convective and reactive part)
      // for conservative form including right-hand-side contribution
      if (conservative)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui = 3*ui;
          const int fuip = fui+1;
          const double v0 = fac*funct_(ui)*velint_(0);
          const double v1 = fac*funct_(ui)*velint_(1);

          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;
            /* convection, convective part */
            /*

            /                                   \
            |  /              n+1  \             |
            | | Du (x) (rho*u)     |  , nabla v  |
            |  \              (i)  /             |
            \                                   /

            */
            estif(fvi,  fui ) -= v0*derxy_(0, vi);
            estif(fvi,  fuip) -= v0*derxy_(1, vi);
            estif(fvip, fui ) -= v1*derxy_(0, vi);
            estif(fvip, fuip) -= v1*derxy_(1, vi);
          }
        }

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui  = 3*ui;
            const int fuip = fui+1;
            const double v = fac*densfunct_(ui);

            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi  = 3*vi;
              const int fvip = fvi+1;
              /* convection, reactive part */
              /*

              /                                     \
              |  / n+1                \             |
              | | u    (x) D(rho*u)   |  , nabla v  |
              |  \ (i)               /              |
              \                                    /

              */
              double v2 = v*ndwconv_c_(vi) ;
              estif(fvi,  fui ) -= v2;
              estif(fvip, fuip) -= v2;
            }
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          const int fvip = fvi+1;
          /* convection */
          double v = fac*ndwconv_c_(vi);
          eforce(fvi ) += v*velint_(0) ;
          eforce(fvip) += v*velint_(1) ;
        }
      }
      // computation of convection (convective and reactive part)
      // for convective form including right-hand-side contribution
      else
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui  = 3*ui;
          const int fuip = fui+1;
          const double v = fac*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;
            /* convection, convective part */
            /*

            /                               \
            |  /       n+1        \         |
            | | (rho*u)   o nabla | Du , v  |
            |  \      (i)        /          |
            \                              /

            */
            double v2 = v*funct_(vi) ;
            estif(fvi , fui ) += v2;
            estif(fvip, fuip) += v2;
          }
        }

        if (newton)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;
            const double v = fac*funct_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui  = 3*ui;
              const int fuip = fui+1;
              const double v2 = v*densfunct_(ui);

              /*  convection, reactive part

              /                                 \
              |  /                \   n+1       |
              | | D(rho*u) o nabla | u     , v  |
              |  \                /   (i)       |
              \                                /
              */
              estif(fvi,  fui ) += v2*vderxy_(0, 0) ;
              estif(fvi,  fuip) += v2*vderxy_(0, 1) ;
              estif(fvip, fui ) += v2*vderxy_(1, 0) ;
              estif(fvip, fuip) += v2*vderxy_(1, 1) ;
            }
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi = 3*vi;
          const int fvip = fvi+1;
          /* convection */
          double v = -fac*funct_(vi);
          eforce(fvi ) += v*conv_old_(0) ;
          eforce(fvip) += v*conv_old_(1) ;
        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        const int fui  = 3*ui;
        const int fuip = fui+1;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          const int fvip = fvi+1;
          const double derxy_0ui_0vi = derxy_(0,ui)*derxy_(0,vi);
          const double derxy_1ui_1vi = derxy_(1,ui)*derxy_(1,vi);
          /* viscosity term */
          /*

                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(fvi,  fui ) += visc*fac*(2.0*derxy_0ui_0vi
                               +
                               derxy_1ui_1vi) ;
          estif(fvi,  fuip) += visc*fac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(fvip, fui ) += visc*fac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(fvip, fuip) += visc*fac*(derxy_0ui_0vi
                                         +
                                         2.0*derxy_1ui_1vi) ;
        }
      }

      for (int ui=0; ui<numnode; ++ui)
      {
        const int fuipp = 3*ui+2;
        const double v = -fac*funct_(ui);
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /*  pressure term */
          /*

          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */

          estif(fvi,     fuipp) += v*derxy_(0, vi) ;
          estif(fvi + 1, fuipp) += v*derxy_(1, vi) ;

        }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvipp = 3*vi+2;
        const double v = fac*functdens_(vi);
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui = 3*ui;

          /* divergence */
          /*
            /                              \
            |                              |
            | nabla o D(rho*u)  , (q/rho)  |
            |                              |
            \                              /
          */
          estif(fvipp, fui    ) += v*densderxy_(0, ui) ;
          estif(fvipp, fui + 1) += v*densderxy_(1, ui) ;
        }
      }

      {
      const double v = press*fac;
      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi  = 3*vi;
        /* pressure */
        eforce(fvi    ) += v*derxy_(0, vi) ;
        eforce(fvi + 1) += v*derxy_(1, vi) ;
      }
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi  = 3*vi;
        /* viscosity */
        eforce(fvi    ) -= visc*fac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                                     +
                                     derxy_(1, vi)*vderxy_(0, 1)
                                     +
                                     derxy_(1, vi)*vderxy_(1, 0)) ;
        eforce(fvi + 1) -= visc*fac*(derxy_(0, vi)*vderxy_(0, 1)
                                      +
                                     derxy_(0, vi)*vderxy_(1, 0)
                                      +
                                     2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<numnode; ++vi)
      {
        const int fvi  = 3*vi;
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eforce(fvi    ) += v*rhsmom_(0) ;
        eforce(fvi + 1) += v*rhsmom_(1) ;
      }

      {
      const double fac_mdiv = fac*mdiv_;
      for (int vi=0; vi<numnode; ++vi)
      {
        // continuity equation
        eforce(vi*3 + 2) -= fac_mdiv*functdens_(vi) ;
      }
      }

      if (loma)
      {
        const double v = -(2.0/3.0)*visc*fac ;
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui  = 3*ui;
          const int fuip = fui+1;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;

            /* viscosity term - subtraction for low-Mach-number flow */
            /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
            */
            estif(fvi,  fui ) += v*derxy_(0, vi)*derxy_(0, ui) ;
            estif(fvi,  fuip) += v*derxy_(0, vi)*derxy_(1, ui) ;
            estif(fvip, fui ) += v*derxy_(1, vi)*derxy_(0, ui) ;
            estif(fvip, fuip) += v*derxy_(1, vi)*derxy_(1, ui) ;

          }
        }

        {
        const double vvdiv = v*vdiv_;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /* viscosity term - subtraction for low-Mach-number flow */
          eforce(fvi    ) -= derxy_(0, vi)*vvdiv ;
          eforce(fvi + 1) -= derxy_(1, vi)*vvdiv ;
        }
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if(pspg == Fluid2::pstab_use_pspg)
      {
        for (int ui=0; ui<numnode; ++ui)
        {
          const int fui  = 3*ui;
          const int fuip = fui+1;
          const double v = tau_Mp*conv_c_(ui) ;
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvipp  = 3*vi+2;
            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | | (rho*u)   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */
            estif(fvipp, fui ) += v*derxy_(0, vi) ;
            estif(fvipp, fuip) += v*derxy_(1, vi) ;
          }
        }

        if (higher_order_ele)
        {
          const double v = -2.0*visc*tau_Mp;
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui  = 3*ui;
            const int fuip = fui+1;
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvipp = 3*vi+2;

              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(fvipp, fui ) += v*(derxy_(0, vi)*viscs2_(0, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(1, ui)) ;
              estif(fvipp, fuip) += v*(derxy_(0, vi)*viscs2_(1, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(3, ui)) ;
            }
          }
        }

        for (int ui=0; ui<numnode; ++ui)
        {
          const int fuipp  = 3*ui+2;
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
            estif(vi*3 + 2, fuipp) += tau_Mp*(derxy_(0, ui)*derxy_(0, vi)
                                              +
                                              derxy_(1, ui)*derxy_(1, vi)) ;

          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui  = 3*ui;
            const int fuip = fui+1;
            const double v = tau_Mp*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvipp  = 3*vi+2;
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | | D(rho*u) o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

              */
              estif(fvipp, fui ) += v*(derxy_(0, vi)*vderxy_(0, 0)
                                       +
                                       derxy_(1, vi)*vderxy_(1, 0)) ;
              estif(fvipp, fuip) += v*(derxy_(0, vi)*vderxy_(0, 1)
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
          const int fui  = 3*ui;
          const int fuip = fui+1;
          const double v = tau_M*conv_c_(ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                         \
            |    /       n+1        \        /       n+1        \     |
            |   | (rho*u)    o nabla | Du , | (rho*u)    o nabla | v  |
            |    \       (i)        /        \       (i)        /     |
            \                                                         /

            */
            const double v2 = v*conv_c_(vi);
            estif(fvi,     fui ) += v2;
            estif(fvi + 1, fuip) += v2;
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          const int fvip = fvi+1;
          const double v = tau_M*conv_c_(vi);
          for (int ui=0; ui<numnode; ++ui)
          {
            const int fuipp = 3*ui+2;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , | (rho*u)   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(fvi,  fuipp) += v*derxy_(0, ui) ;
            estif(fvip, fuipp) += v*derxy_(1, ui) ;
          }
        }

        if (higher_order_ele)
        {
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;
            const double v = -2.0*visc*tau_M*conv_c_(vi);
            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui  = 3*ui;
              const int fuip = fui+1;
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, | (rho*u)    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
              */
              estif(fvi,  fui ) += v*viscs2_(0, ui) ;
              estif(fvip, fui ) += v*viscs2_(1, ui) ;

              estif(fvi,  fuip) += v*viscs2_(1, ui) ;
              estif(fvip, fuip) += v*viscs2_(3, ui) ;
            }
          }
        }

        if (newton)
        {
          {
            const double v0 = velint_(0)*vderxy_(0, 0) + velint_(1)*vderxy_(0, 1);
            const double v1 = velint_(0)*vderxy_(1, 0) + velint_(1)*vderxy_(1, 1);

            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui  = 3*ui;
              const int fuip = fui+1;
              const double v = tau_M*densfunct_(ui);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi  = 3*vi;
                const int fvip = fvi+1;
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
                estif(fvi,  fui ) += (conv_c_(vi)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(fvip, fui ) += (conv_c_(vi)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;

                estif(fvi,  fuip) += (conv_c_(vi)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(fvip, fuip) += (conv_c_(vi)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
              }
            }
          }

          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui  = 3*ui;
            const int fuip = fui+1;
            const double v = tau_M*densfunct_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi  = 3*vi;
              const int fvip = fvi+1;
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , | D(rho*u) o nabla | v  |
                |         (i)    \                /     |
                \                                      /
              */
              estif(fvi,  fui ) += v*gradp_(0)*derxy_(0, vi) ;
              estif(fvip, fui ) += v*gradp_(1)*derxy_(0, vi) ;

              estif(fvi,  fuip) += v*gradp_(0)*derxy_(1, vi) ;
              estif(fvip, fuip) += v*gradp_(1)*derxy_(1, vi) ;

            }
          }

          if (higher_order_ele)
          {
            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui  = 3*ui;
              const int fuip = fui+1;
              const double v = -2.0*visc*tau_M*densfunct_(ui);
              const double v0 = v*visc_old_(0);
              const double v1 = v*visc_old_(1);
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi  = 3*vi;
                const int fvip = fvi+1;
                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, | D(rho*u) o nabla | v  |
                  |               \ (i) /    \                /     |
                  \                                                 /
                */
                estif(fvi,  fui ) += v0*derxy_(0, vi) ;
                estif(fvip, fui ) += v1*derxy_(0, vi) ;

                estif(fvi,  fuip) += v0*derxy_(1, vi) ;
                estif(fvip, fuip) += v1*derxy_(1, vi) ;

              }
            }
          }

          for (int ui=0; ui<numnode; ++ui)
          {
            const int fui  = 3*ui;
            const int fuip = fui+1;
            const double v = -tau_M*densfunct_(ui);
            const double v0 = v*rhsmom_(0);
            const double v1 = v*rhsmom_(1);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi  = 3*vi;
              const int fvip = fvi+1;
              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                                     \
                |              /                \     |
                |  rhsint   , | D(rho*u) o nabla | v  |
                |              \                /     |
                \                                     /

              */
              estif(fvi , fui ) += v0*derxy_(0, vi) ;
              estif(fvip, fui ) += v1*derxy_(0, vi) ;

              estif(fvi , fuip) += v0*derxy_(1, vi) ;
              estif(fvip, fuip) += v1*derxy_(1, vi) ;

            } // vi
          } // ui
        } // if newton

        /* right hand side */
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          // supg stabilisation
          const double v = -tau_M*conv_c_(vi);
          eforce(fvi    ) += v*res_old_(0) ;
          eforce(fvi + 1) += v*res_old_(1) ;
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
              const int fui  = 3*ui;
              const int fuip = fui+1;
              const double v = two_visc_tauMp*conv_c_(ui)
                         ;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi  = 3*vi;
                const int fvip = fvi+1;
                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | | (rho*u)   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
                */
                estif(fvi,  fui ) += v*viscs2_(0, vi) ;
                estif(fvip, fui ) += v*viscs2_(1, vi) ;

                estif(fvi,  fuip) += v*viscs2_(1, vi) ;
                estif(fvip, fuip) += v*viscs2_(3, vi) ;
              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              const int fuipp = 3*ui+2;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi  = 3*vi;

                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(fvi,     fuipp) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                         +
                                                         derxy_(1, ui)*viscs2_(1, vi)) ;
                estif(fvi + 1, fuipp) += two_visc_tauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                         +
                                                         derxy_(1, ui)*viscs2_(3, vi)) ;

              }
            }

            for (int ui=0; ui<numnode; ++ui)
            {
              const int fui  = 3*ui;
              const int fuip = fui+1;
              for (int vi=0; vi<numnode; ++vi)
              {
                const int fvi  = 3*vi;
                const int fvip = fvi+1;
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(fvi,  fui ) -= four_visc2_tauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)) ;
                estif(fvip, fui ) -= four_visc2_tauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(3,vi)) ;

                estif(fvi,  fuip) -= four_visc2_tauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(3,ui)) ;
                estif(fvip, fuip) -= four_visc2_tauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(3,ui)*viscs2_(3,vi)) ;
              } // vi
            } // ui

            if (newton)
            {
              for (int ui=0; ui<numnode; ++ui)
              {
                const int fui  = 3*ui;
                const int fuip = fui+1;
                const double v = two_visc_tauMp*densfunct_(ui);
                for (int vi=0; vi<numnode; ++vi)
                {
                  const int fvi  = 3*vi;
                  const int fvip = fvi+1;
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | | D(rho*u) o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
                  */
                  estif(fvi,  fui ) += v*(viscs2_(0,vi)*vderxy_(0,0)+viscs2_(1,vi)*vderxy_(1,0)) ;
                  estif(fvip, fui ) += v*(viscs2_(1,vi)*vderxy_(0,0)+viscs2_(3,vi)*vderxy_(1,0)) ;

                  estif(fvi,  fuip) += v*(viscs2_(0,vi)*vderxy_(0,1)+viscs2_(1,vi)*vderxy_(1,1)) ;
                  estif(fvip, fuip) += v*(viscs2_(1,vi)*vderxy_(0,1)+viscs2_(3,vi)*vderxy_(1,1)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            /* viscous stabilisation */
            eforce(fvi    ) -= two_visc_tauMp*(res_old_(0)*viscs2_(0, vi)+res_old_(1)*viscs2_(1, vi)) ;
            eforce(fvi + 1) -= two_visc_tauMp*(res_old_(0)*viscs2_(1, vi)+res_old_(1)*viscs2_(3, vi)) ;
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
          const int fui  = 3*ui;
          const int fuip = fui+1;
          const double v0 = tau_C*densderxy_(0, ui);
          const double v1 = tau_C*densderxy_(1, ui);
          for (int vi=0; vi<numnode; ++vi)
          {
            const int fvi  = 3*vi;
            const int fvip = fvi+1;
            /* continuity stabilisation on left hand side */
            /*
              /                                      \
              |                                      |
              | nabla o D(rho*u)  , nabla o (rho*v)  |
              |                                      |
              \                                      /
            */
            estif(fvi,  fui ) += v0*densderxy_(0, vi) ;
            estif(fvip, fui ) += v0*densderxy_(1, vi) ;

            estif(fvi,  fuip) += v1*densderxy_(0, vi) ;
            estif(fvip, fuip) += v1*densderxy_(1, vi) ;
          }
        }

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /* continuity stabilisation on right hand side */
          eforce(fvi    ) -= tau_C_divunp*densderxy_(0, vi) ;
          eforce(fvi + 1) -= tau_C_divunp*densderxy_(1, vi) ;
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
            const int fui  = 3*ui;
            const int fuip = fui+1;
            const double v = tau_M*conv_resM_(ui);
            for (int vi=0; vi<numnode; ++vi)
            {
              const int fvi  = 3*vi;
              const int fvip = fvi+1;
              /* cross-stress part on lhs */
              /*

                          /                               \
                         |  /                  \          |
                      -  | | (rho*resM) o nabla | Du , v  |
                         |  \                  /          |
                          \                               /
              */
              const double v2 = v*funct_(vi);
              estif(fvi , fui ) -= v2 ;
              estif(fvip, fuip) -= v2 ;
            }
          }
        } // end cross-stress part on left hand side

        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /* cross-stress part on rhs */
          /*

                          /                                \
                         |  /                  \           |
                         | | (rho*resM) o nabla | u   , v  |
                         |  \                  /  (i)      |
                          \                               /
          */
          double v = tau_M*funct_(vi);
          eforce(fvi    ) += v*(res_old_(0)*mderxy_(0,0)+res_old_(1)*mderxy_(0,1));
          eforce(fvi + 1) += v*(res_old_(0)*mderxy_(1,0)+res_old_(1)*mderxy_(1,1));
        }
      } // end cross-stress part on right hand side

      //----------------------------------------------------------------------
      //     STABILIZATION, REYNOLDS-STRESS PART (RESIDUAL-BASED VMM)

      if (reynolds == Fluid2::reynolds_stress_stab_only_rhs)
      {
        const double tauMtauM = tau_M*tau_M/fac;
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /* Reynolds-stress part on rhs */
          /*

                  /                                    \
                 |                                     |
                 |  resM   ,   (rho*resM) o nabla ) v  |
                 |                                     |
                  \                                    /
          */
          const double v = tauMtauM*conv_resM_(vi);
          eforce(fvi    ) += v*res_old_(0);
          eforce(fvi + 1) += v*res_old_(1);
        }
      } // end Reynolds-stress part on right hand side

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv != Fluid2::no_fssgv)
      {
        for (int vi=0; vi<numnode; ++vi)
        {
          const int fvi  = 3*vi;
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(fvi    ) -= vartfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                      +    derxy_(1, vi)*fsvderxy_(0, 1)
                                      +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
          eforce(fvi + 1) -= vartfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Stationary<distype>::CalTauStationary(
  Fluid2*                             ele,
  const LINALG::Matrix<2,iel>&        evelnp,
  const LINALG::Matrix<2,iel>&        fsevelnp,
  const LINALG::Matrix<iel,1>&        edensnp,
  const double                        visc,
  const enum Fluid2::FineSubgridVisc  fssgv,
  const double                        Cs
  )
{

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
  velint_.Multiply(evelnp, funct_);

  // get density at element center
  const double dens = funct_.Dot(edensnp);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_, xyze_);
  // inverse of jacobian
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
  derxy_.Multiply(xji_, deriv_);

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
  if (fssgv == Fluid2::artificial_all or fssgv == Fluid2::artificial_small)
  {
    double fsvel_norm = 0.0;
    if (fssgv == Fluid2::artificial_small)
    {
      // get fine-scale velocities at element center
      fsvelint_.Multiply(fsevelnp, funct_);

      // get fine-scale velocity norm
      fsvel_norm = fsvelint_.Norm2();
    }
    // get all-scale velocity norm
    else fsvel_norm = vel_norm;

    /*----------------------------- compute artificial subgrid viscosity ---*/
    const double re_sv = mk * dens * fsvel_norm * hk / visc; /* convective : viscous forces */
    const double xi_sv = DMAX(re_sv,1.0);

    vart_ = (DSQR(hk)*mk*DSQR(dens)*DSQR(fsvel_norm))/(2.0*visc*xi_sv);

  }
  else if (fssgv == Fluid2::smagorinsky_all or
           fssgv == Fluid2::smagorinsky_small)
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
      if (fssgv == Fluid2::smagorinsky_small)
            fsvderxy_.MultiplyNT(fsevelnp, derxy_);
      else  fsvderxy_.MultiplyNT(evelnp, derxy_);

      double epsilon;
      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = fsvderxy_(mm,rr) + fsvderxy_(rr,mm);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0*fsvderxy_(rr,rr)*fsvderxy_(rr,rr);
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



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2Stationary<distype>::BodyForce(Fluid2*      ele,
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
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
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
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,x,time,NULL);
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
    edeadng_.Clear();
  }
}


#endif
#endif
