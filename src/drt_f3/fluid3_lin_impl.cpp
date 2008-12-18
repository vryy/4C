/*----------------------------------------------------------------------*/
/*!
\file fluid3_lin_impl.cpp

\brief Internal implementation of linearised (fast) Fluid3 element

<pre>
Maintainer: Christiane Foerster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_lin_impl.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"

#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3lin_ImplInterface* DRT::ELEMENTS::Fluid3lin_ImplInterface::Impl(DRT::ELEMENTS::Fluid3* f3)
{
  switch (f3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Fluid3lin_Impl<DRT::Element::hex8>* fh8;
    if (fh8==NULL)
      fh8 = new Fluid3lin_Impl<DRT::Element::hex8>();
    return fh8;
  }
  case DRT::Element::hex20:
  {
    static Fluid3lin_Impl<DRT::Element::hex20>* fh20;
    if (fh20==NULL)
      fh20 = new Fluid3lin_Impl<DRT::Element::hex20>();
    return fh20;
  }
  case DRT::Element::hex27:
  {
    static Fluid3lin_Impl<DRT::Element::hex27>* fh27;
    if (fh27==NULL)
      fh27 = new Fluid3lin_Impl<DRT::Element::hex27>();
    return fh27;
  }
  case DRT::Element::tet4:
  {
    static Fluid3lin_Impl<DRT::Element::tet4>* ft4;
    if (ft4==NULL)
      ft4 = new Fluid3lin_Impl<DRT::Element::tet4>();
    return ft4;
  }
  case DRT::Element::tet10:
  {
    static Fluid3lin_Impl<DRT::Element::tet10>* ft10;
    if (ft10==NULL)
      ft10 = new Fluid3lin_Impl<DRT::Element::tet10>();
    return ft10;
  }
  case DRT::Element::wedge6:
  {
    static Fluid3lin_Impl<DRT::Element::wedge6>* fw6;
    if (fw6==NULL)
      fw6 = new Fluid3lin_Impl<DRT::Element::wedge6>();
    return fw6;
  }
  case DRT::Element::wedge15:
  {
    static Fluid3lin_Impl<DRT::Element::wedge15>* fw15;
    if (fw15==NULL)
      fw15 = new Fluid3lin_Impl<DRT::Element::wedge15>();
    return fw15;
  }
  case DRT::Element::pyramid5:
  {
    static Fluid3lin_Impl<DRT::Element::pyramid5>* fp5;
    if (fp5==NULL)
      fp5 = new Fluid3lin_Impl<DRT::Element::pyramid5>();
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
DRT::ELEMENTS::Fluid3lin_Impl<distype>::Fluid3lin_Impl()
  : Fluid3lin_ImplInterface(),
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
    pderxy_(),
    vderxy2_(),
    derxy_(),
    densderxy_(),
    derxy2_(),
    bodyforce_(),
    histmom_(),
    histcon_(),
    velino_(),
    velint_(),
    gradp_(),
    tau_(),
    viscs2_(),
    conv_(),
    mdiv_(),
    rhsmom_(),
    rhscon_(),
    xder2_(),
    numepn_(iel)
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3lin_Impl<distype>::Evaluate(
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
  const int numnode = iel;

  // construct views
  LINALG::Matrix<4*iel,4*iel> elemat1(elemat1_epetra.A(), true);
  //LINALG::Matrix<4*iel,4*iel> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<4*iel,1> elevec1(elevec1_epetra.A(), true);
  // elemat2, elevec2 and elevec3 are never used anyway

  // need current velocity and history vector
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

  RefCountPtr<const Epetra_Vector> dispnp;
  vector<double> mydispnp;

  // create objects for element arrays
  LINALG::Matrix<numnode,1> eprenp;
  LINALG::Matrix<3,numnode> evelnp;
  LINALG::Matrix<numnode,1> edensnp;
  LINALG::Matrix<3,numnode> emhist;
  LINALG::Matrix<numnode,1> echist;

  // split velocity and pressure, insert into element arrays
  for (int i=0;i<numnode;++i)
  {
    evelnp(0,i) = myvelnp[0+(i*4)];
    evelnp(1,i) = myvelnp[1+(i*4)];
    evelnp(2,i) = myvelnp[2+(i*4)];

    eprenp(i) = myvelnp[3+(i*4)];

    // insert density vector into element array
    edensnp(i,0) = myvedenp[3+(i*4)];

    // the history vectors contain information of time step t_n (mass rhs!)
    // momentum equation part
    emhist(0,i) = myhist[0+(i*4)];
    emhist(1,i) = myhist[1+(i*4)];
    emhist(2,i) = myhist[2+(i*4)];

    // continuity equation part (only non-trivial for low-Mach-number flow)
    echist(i,0) = myhist[3+(i*4)];
  }

  // set parameter for potential low-Mach-number solver
  string lomastr  =params.get<string>("low-Mach-number solver");
  bool loma   = false;
  if(lomastr  =="Yes") loma  =true;

  // get control parameter
  const double time = params.get<double>("total time",-1.0);

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // calculate element coefficient matrix and rhs
  Sysmat(ele,
         evelnp,
         eprenp,
         edensnp,
         emhist,
         echist,
         elemat1,
         elevec1,
         actmat,
         time,
         timefac,
         loma);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                chfoe 02/08|
 *----------------------------------------------------------------------*/
/*
 Note: This routine is made for total rather than incremental solution.

       THE PRESENT FLUID ELEMENT IS ENTIRELY LINEAR!!!!!!!!!!!!!!!
*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3lin_Impl<distype>::Sysmat(
  Fluid3*                                              ele,
  const LINALG::Matrix<3,iel>&     evelnp,
  const LINALG::Matrix<iel,1>&     eprenp,
  const LINALG::Matrix<iel,1>&     edensnp,
  const LINALG::Matrix<3,iel>&     emhist,
  const LINALG::Matrix<iel,1>&     echist,
  LINALG::Matrix<4*iel,4*iel>&     estif,
  LINALG::Matrix<4*iel,1>&         eforce,
  struct _MATERIAL*                                    material,
  double                                               time,
  double                                               timefac,
  bool                                                 loma
  )
{
  const double numnode = iel;

  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];

    numepn_[inode] = nodes[inode]->NumElement();
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  Caltau(ele,evelnp,edensnp,visc,timefac);

  // flag for higher order elements
  const bool higher_order_ele = ele->isHigherOrderElement(distype);

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordiantes of the current integration point
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
    const double fac = intpoints.qwgt[iquad]*det;

    if (det < 0.0)
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

    // compute global derivates
    //derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);
    derxy_.Multiply(xji_,deriv_);

    // (inverse-)density-weighted shape functions and global derivatives
    // remark: this loop is to be replaced soon!
    for (int inode=0; inode<iel; inode++)
    {
      densfunct_(inode,0) = edensnp(inode,0)*funct_(inode,0);
      functdens_(inode,0) = funct_(inode,0)/edensnp(inode,0);

      densderxy_(0,inode) = edensnp(inode,0)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode,0)*derxy_(1,inode);
      densderxy_(2,inode) = edensnp(inode,0)*derxy_(2,inode);
    }

    // compute second global derivative
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);

      // calculate 2nd velocity derivatives at integration point
      //vderxy2_ = blitz::sum(derxy2_(j,k)*evelnp(i,k),k);
      vderxy2_.MultiplyNT(evelnp,derxy2_);
    }
    else
    {
      derxy2_.Clear();
      vderxy2_.Clear();
    }

    // get density-weighted velocities (n+g,i) at integration point
    //velint_ = blitz::sum(densfunct_(j)*evelnp(i,j),j);
    velint_.Multiply(evelnp,densfunct_);

    // get history data (n,i) at integration point
    //histmom_ = blitz::sum(funct_(j)*emhist(i,j),j);
    //histcon_ = blitz::sum(funct_*echist);
    histmom_.Multiply(emhist,densfunct_);
    histcon_ = funct_.Dot(echist);

    // get velocity (np,i) derivatives at integration point
    //vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);
    vderxy_.MultiplyNT(evelnp,derxy_);

    // get momentum (np,i) derivatives at integration point
    //mderxy_ = blitz::sum(densderxy_(j,k)*evelnp(i,k),k);
    mderxy_.MultiplyNT(evelnp,densderxy_);

    // get pressure gradients
    //gradp_ = blitz::sum(derxy_(i,j)*eprenp(j),j);
    gradp_.Multiply(derxy_,eprenp);

    // get density at integration point
    //double dens = blitz::sum(funct_*edensnp);
    double dens = funct_.Dot(edensnp);

    // get (density-weighted) bodyforce in gausspoint
    //bodyforce_ = blitz::sum(edeadng_(i,j)*densfunct_(j),j);
    bodyforce_.Multiply(edeadng_,densfunct_);

    // perform integration for entire matrix and rhs

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    // integration factors and coefficients of single terms
    const double timetauM   = timefac * tau_M;
    const double timetauMp  = timefac * tau_Mp;
    const double timetau_C  = timefac * tau_C;

    const double ttimetauM  = timefac * timetauM;
    const double ttimetauMp = timefac * timetauMp;
    const double timefacfac = timefac * fac;

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    //rhsmom_ = histmom_(i) + bodyforce_(i)*timefac;
    rhsmom_.Update(1.0,histmom_,timefac,bodyforce_);
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    //conv_ = blitz::sum(derxy_(j,i)*velint_(j), j);
    conv_.MultiplyTN(derxy_,velint_);

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

      for (int i=0; i<numnode; ++i) {
        viscs2_(0,i) = 0.5 * (2.0 * derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
        viscs2_(1,i) = 0.5 *  derxy2_(3,i);
        viscs2_(2,i) = 0.5 *  derxy2_(4,i);
        viscs2_(3,i) = 0.5 *  derxy2_(3,i);
        viscs2_(4,i) = 0.5 * (derxy2_(0,i) + 2.0 * derxy2_(1,i) + derxy2_(2,i));
        viscs2_(5,i) = 0.5 *  derxy2_(5,i);
        viscs2_(6,i) = 0.5 *  derxy2_(4,i);
        viscs2_(7,i) = 0.5 *  derxy2_(5,i);
        viscs2_(8,i) = 0.5 * (derxy2_(0,i) + derxy2_(1,i) + 2.0 * derxy2_(2,i));
      }

      if (loma)
      {
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

        for (int i=0; i<numnode; ++i) {
          viscs2_(0,i) -=  derxy2_(0,i)/3.0;
          viscs2_(1,i) -=  derxy2_(3,i)/3.0;
          viscs2_(2,i) -=  derxy2_(4,i)/3.0;
          viscs2_(3,i) -=  derxy2_(3,i)/3.0;
          viscs2_(4,i) -=  derxy2_(1,i)/3.0;
          viscs2_(5,i) -=  derxy2_(5,i)/3.0;
          viscs2_(6,i) -=  derxy2_(4,i)/3.0;
          viscs2_(7,i) -=  derxy2_(5,i)/3.0;
          viscs2_(8,i) -=  derxy2_(2,i)/3.0;
        }
      }
    }
    else
    {
      viscs2_.Clear();
    }

    /* momentum divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1) + mderxy_(2, 2);

    /*--------------------------------- now build single stiffness terms ---*/

    //----------------------------------------------------------------------
    //                            GALERKIN PART

    for (int ui=0; ui<iel; ++ui)
    {
      for (int vi=0; vi<iel; ++vi)
      {
	/* an auxiliary variable for speedup */
        double aux;

	/* inertia (contribution to mass matrix) */
	/*
                               /       \
                              |  u , v  |
                               \       /
	*/
        aux = fac*densfunct_(ui)*funct_(vi);
	estif(vi*4, ui*4)         += aux ;
	estif(vi*4 + 1, ui*4 + 1) += aux ;
	estif(vi*4 + 2, ui*4 + 2) += aux ;

	/* convection, convective part */
	/*
                     /                      \
                    |  u_old  o  grad u , v  |
                     \                      /
	*/
        aux = timefacfac*funct_(vi)*conv_(ui);
	estif(vi*4, ui*4)         += aux;
	estif(vi*4 + 1, ui*4 + 1) += aux;
	estif(vi*4 + 2, ui*4 + 2) += aux;

	/* convection, correction part for linearisation*/
	/*
                   1  /                  \
                   - |  div u_old  u , v  |
                   2  \                  /
	*/
        aux = 0.5*timefacfac*mdiv_*funct_(ui)*funct_(vi);
	estif(vi*4, ui*4)         += aux ;
	estif(vi*4 + 1, ui*4 + 1) += aux ;
	estif(vi*4 + 2, ui*4 + 2) += aux ;

	/* Viskositaetsterm */
	/*
                        /                 \
                       |  eps(u) , eps(v)  |
                        \                 /
	*/
        aux = visc*timefacfac;
	estif(vi*4, ui*4)         += aux * ( 2.0*derxy_(0, ui)*derxy_(0, vi)
					    +derxy_(1, ui)*derxy_(1, vi)
					    +derxy_(2, ui)*derxy_(2, vi)) ;
	estif(vi*4, ui*4 + 1)     += aux * derxy_(0, ui)*derxy_(1, vi) ;
	estif(vi*4, ui*4 + 2)     += aux * derxy_(0, ui)*derxy_(2, vi) ;
	estif(vi*4 + 1, ui*4)     += aux * derxy_(0, vi)*derxy_(1, ui) ;
	estif(vi*4 + 1, ui*4 + 1) += aux * ( derxy_(0, ui)*derxy_(0, vi)
					    +2.0*derxy_(1, ui)*derxy_(1, vi)
					    +derxy_(2, ui)*derxy_(2, vi)) ;
	estif(vi*4 + 1, ui*4 + 2) += aux * derxy_(1, ui)*derxy_(2, vi) ;
	estif(vi*4 + 2, ui*4)     += aux * derxy_(0, vi)*derxy_(2, ui) ;
	estif(vi*4 + 2, ui*4 + 1) += aux * derxy_(1, vi)*derxy_(2, ui) ;
	estif(vi*4 + 2, ui*4 + 2) += aux * ( derxy_(0, ui)*derxy_(0, vi)
					    +derxy_(1, ui)*derxy_(1, vi)
					    +2.0*derxy_(2, ui)*derxy_(2, vi)) ;

        if (loma)
        {
          aux = -(2.0/3.0)*visc*timefacfac;
          /* viscosity term - subtraction for low-Mach-number flow */
          /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
          */
          estif(vi*4,   ui*4  ) += aux*derxy_(0, vi)*derxy_(0, ui) ;
          estif(vi*4,   ui*4+1) += aux*derxy_(0, vi)*derxy_(1, ui) ;
          estif(vi*4,   ui*4+2) += aux*derxy_(0, vi)*derxy_(2, ui) ;
          estif(vi*4+1, ui*4  ) += aux*derxy_(1, vi)*derxy_(0, ui) ;
          estif(vi*4+1, ui*4+1) += aux*derxy_(1, vi)*derxy_(1, ui) ;
          estif(vi*4+1, ui*4+2) += aux*derxy_(1, vi)*derxy_(2, ui) ;
          estif(vi*4+2, ui*4  ) += aux*derxy_(2, vi)*derxy_(0, ui) ;
          estif(vi*4+2, ui*4+1) += aux*derxy_(2, vi)*derxy_(1, ui) ;
          estif(vi*4+2, ui*4+2) += aux*derxy_(2, vi)*derxy_(2, ui) ;
        }

	/* Druckterm */
	/*
                          /           \
                         |  p , div v  |
                          \           /
	*/
        aux = timefacfac*funct_(ui);
	estif(vi*4, ui*4 + 3)     -= aux * derxy_(0, vi);
	estif(vi*4 + 1, ui*4 + 3) -= aux * derxy_(1, vi);
	estif(vi*4 + 2, ui*4 + 3) -= aux * derxy_(2, vi);

	/* Divergenzfreiheit */
	/*
                         /          \
                        | div u , q  |
                         \          /
	*/
        aux = timefacfac*functdens_(vi);
	estif(vi*4 + 3, ui*4)     += aux * densderxy_(0, ui) ;
	estif(vi*4 + 3, ui*4 + 1) += aux * densderxy_(1, ui) ;
	estif(vi*4 + 3, ui*4 + 2) += aux * densderxy_(2, ui) ;


	//-----------------------------------------------------------------
	//                           STABILISATION PART

	/* pressure stabilisation: inertia PLUS */
	/* pressure stabilisation: convection, convective part */
	/*
                        /            \     /                           \
                       |  u , grad q  | + |  u_old  o  grad u , grad q  |
                        \            /     \                           /
	*/
	aux = timetauMp*densfunct_(ui) + ttimetauMp*conv_(ui);
	estif(vi*4 + 3, ui*4)     += aux * derxy_(0, vi) ;
	estif(vi*4 + 3, ui*4 + 1) += aux * derxy_(1, vi) ;
	estif(vi*4 + 3, ui*4 + 2) += aux * derxy_(2, vi) ;

        if (higher_order_ele)
        {
          /* pressure stabilisation: viscosity (-L_visc_u) */
          /*
                     /                     \
                    |  div eps(u) , grad q  |
                     \                     /
          */
          aux = 2.0*visc*ttimetauMp;
          estif(vi*4 + 3, ui*4)     += aux*( derxy_(0, vi)*viscs2_(0, ui)
                                            +derxy_(1, vi)*viscs2_(1, ui)
                                            +derxy_(2, vi)*viscs2_(2, ui));
          estif(vi*4 + 3, ui*4 + 1) += aux*( derxy_(0, vi)*viscs2_(1, ui)
                                            +derxy_(1, vi)*viscs2_(4, ui)
                                            +derxy_(2, vi)*viscs2_(5, ui));
          estif(vi*4 + 3, ui*4 + 2) += aux*( derxy_(0, vi)*viscs2_(2, ui)
                                            +derxy_(1, vi)*viscs2_(5, ui)
                                            +derxy_(2, vi)*viscs2_(8, ui));
        }

	/* pressure stabilisation: pressure( L_pres_p) */
	/*
	                /	          \
                       |  grad p , grad q  |
                        \                 /
	*/
	estif(vi*4 + 3, ui*4 + 3) += ttimetauMp*( derxy_(0, ui)*derxy_(0, vi)
						 +derxy_(1, ui)*derxy_(1, vi)
						 +derxy_(2, ui)*derxy_(2, vi));

	//----------------------------------------------------------------

	/* supg stabilisation: inertia  PLUS */
	/* supg stabilisation: convective part ( L_conv_u) */

	/*
                      /                     \     /                                        \
                     |  u ,  u_old  grad  v  | + |   u_old  o  grad u ,  u_old  o  grad  v  |
                      \                     /     \                                        /
	*/
        aux = timetauM*densfunct_(ui)*conv_(vi) + ttimetauM*conv_(ui)*conv_(vi);
	estif(vi*4, ui*4)         += aux;
	estif(vi*4 + 1, ui*4 + 1) += aux;
	estif(vi*4 + 2, ui*4 + 2) += aux;

        if (higher_order_ele)
        {
          /* supg stabilisation: viscous part  (-L_visc_u) */
          /*
                  /                               \
                 |  div eps(u),  u_old  o  grad v  |
                  \                               /
          */
          aux = 2.0*visc*ttimetauM*conv_(vi);
          estif(vi*4, ui*4)         += aux*viscs2_(0, ui) ;
          estif(vi*4, ui*4 + 1)     += aux*viscs2_(1, ui) ;
          estif(vi*4, ui*4 + 2)     += aux*viscs2_(2, ui) ;
          estif(vi*4 + 1, ui*4)     += aux*viscs2_(1, ui) ;
          estif(vi*4 + 1, ui*4 + 1) += aux*viscs2_(4, ui) ;
          estif(vi*4 + 1, ui*4 + 2) += aux*viscs2_(5, ui) ;
          estif(vi*4 + 2, ui*4)     += aux*viscs2_(2, ui) ;
          estif(vi*4 + 2, ui*4 + 1) += aux*viscs2_(5, ui) ;
          estif(vi*4 + 2, ui*4 + 2) += aux*viscs2_(8, ui) ;
        }

	/* supg stabilisation: pressure part  ( L_pres_p) */
	/*
                      /                            \
                     |  grad p ,  u_old  o  grad v  |
                      \                            /
	*/
	aux = ttimetauM*conv_(vi);
	estif(vi*4, ui*4 + 3)     += aux*derxy_(0, ui) ;
	estif(vi*4 + 1, ui*4 + 3) += aux*derxy_(1, ui) ;
	estif(vi*4 + 2, ui*4 + 3) += aux*derxy_(2, ui) ;

	//----------------------------------------------------------------

	/* continuity stabilisation */
	/*
                       /              \
                      | div u , div v  |
                       \              /
	*/
	estif(vi*4, ui*4)         += timetau_C*densderxy_(0, ui)*densderxy_(0, vi);
	estif(vi*4, ui*4 + 1)     += timetau_C*densderxy_(0, vi)*densderxy_(1, ui);
	estif(vi*4, ui*4 + 2)     += timetau_C*densderxy_(0, vi)*densderxy_(2, ui);
	estif(vi*4 + 1, ui*4)     += timetau_C*densderxy_(0, ui)*densderxy_(1, vi);
	estif(vi*4 + 1, ui*4 + 1) += timetau_C*densderxy_(1, ui)*densderxy_(1, vi);
	estif(vi*4 + 1, ui*4 + 2) += timetau_C*densderxy_(1, vi)*densderxy_(2, ui);
	estif(vi*4 + 2, ui*4)     += timetau_C*densderxy_(0, ui)*densderxy_(2, vi);
	estif(vi*4 + 2, ui*4 + 1) += timetau_C*densderxy_(1, ui)*densderxy_(2, vi);
	estif(vi*4 + 2, ui*4 + 2) += timetau_C*densderxy_(2, ui)*densderxy_(2, vi);
        }
      }

    /*------------------------------------ now build single rhs terms ---*/
    for (int vi=0; vi<iel; ++vi)
    {
      double aux;
      //-------------------------------------------------------------------
      //                            GALERKIN PART

      // source term of the right hand side
      aux = fac*funct_(vi);
      eforce(vi*4)     += aux*rhsmom_(0) ;
      eforce(vi*4 + 1) += aux*rhsmom_(1) ;
      eforce(vi*4 + 2) += aux*rhsmom_(2) ;

      /* rhs term of continuity equation */
      if (loma) eforce(vi*4 + 3) += fac*functdens_(vi)*rhscon_ ;

      //-------------------------------------------------------------------
      //                           STABILISATION PART

      // pressure stabilisation
      eforce(vi*4 + 3) += timetauMp*( rhsmom_(0)*derxy_(0, vi)
				     +rhsmom_(1)*derxy_(1, vi)
				     +rhsmom_(2)*derxy_(2, vi)) ;
      // supg stabilisation
      aux =  timetauM*conv_(vi);
      eforce(vi*4)     += aux*rhsmom_(0) ;
      eforce(vi*4 + 1) += aux*rhsmom_(1) ;
      eforce(vi*4 + 2) += aux*rhsmom_(2) ;

      if (loma)
      {
        // continuity stabilisation of rhs term of continuity equation
        aux = timetau_C*rhscon_;
        eforce(vi*4    ) += densderxy_(0, vi)*aux ;
        eforce(vi*4 + 1) += densderxy_(1, vi)*aux ;
        eforce(vi*4 + 2) += densderxy_(2, vi)*aux ;
      }
    } // vi
  }
}




//
// calculate stabilization parameter
//
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3lin_Impl<distype>::Caltau(
  Fluid3*                                           ele,
  const LINALG::Matrix<3,iel>&  evelnp,
  const LINALG::Matrix<iel,1>&  edensnp,
  const double                                      visc,
  const double                                      timefac
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

  // get density at element center
  //const double dens = blitz::sum(funct_*edensnp);
  const double dens = funct_.Dot(edensnp);

  // get Jacobian matrix and determinant
  //xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xji_.Invert(xjm_);

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
    velino_.Update(1.0/vel_norm,velint_);
  }
  else
  {
    velino_ = 0.;
    velino_(0,0) = 1;
  }

  // get streamlength
  //const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  LINALG::Matrix<iel,1> tmp;
  tmp.MultiplyTN(derxy_,velino_);
  const double val = tmp.Norm1();
  const double strle = 2.0/val;

  // calculate tau

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
  const double re1 = 4.0 * timefac * visc / (mk * dens * DSQR(strle));

  /* convective : viscous forces */
  const double re2 = mk * dens * vel_norm * strle / (2.0 * visc);

  const double xi1 = DMAX(re1,1.0);
  const double xi2 = DMAX(re2,1.0);

  tau_(0,0) = DSQR(strle)/(DSQR(strle)*dens*xi1+(4.0*timefac*visc/mk)*xi2);

  // compute tau_Mp
  //    stability parameter definition according to Franca and Valentin (2000)
  //                                       and Barrenechea and Valentin (2002)

  /* viscous : reactive forces */
  const double re_viscous = 4.0 * timefac * visc / (mk * dens * DSQR(hk));
  /* convective : viscous forces */
  const double re_convect = mk * dens * vel_norm * hk / (2.0 * visc);

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
    tau_(1) = DSQR(hk)/(DSQR(hk)*dens*xi_viscous+(4.0*timefac*visc/mk)*xi_convect);

  /*------------------------------------------------------ compute tau_C ---*/
  /*-- stability parameter definition according to Codina (2002), CMAME 191
   *
   * Analysis of a stabilized finite element approximation of the transient
   * convection-diffusion-reaction equation using orthogonal subscales.
   * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
   *
   * */
  tau_(2) = sqrt(4.0*DSQR(visc)+DSQR(0.5*vel_norm*hk/dens));

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
  //const double xi_tau_c = DMIN(re2,1.0);
  //tau_(2) = vel_norm * hk * 0.5 * xi_tau_c / ( timefac * dens );
}



// NOTE: this is an exact copy of the routine in fluid3_impl.cpp
/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3lin_Impl<distype>::BodyForce( Fluid3* ele,
                                                        const double time)
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

  return;
}


#endif
#endif
