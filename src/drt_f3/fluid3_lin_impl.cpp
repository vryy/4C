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
#include "../drt_lib/drt_timecurve.H"

#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3lin_Impl* DRT::ELEMENTS::Fluid3lin_Impl::Impl(DRT::ELEMENTS::Fluid3* f3)
{
  switch (f3->NumNode())
  {
  case 8:
  {
    static Fluid3lin_Impl* f8;
    if (f8==NULL)
      f8 = new Fluid3lin_Impl(8);
    return f8;
  }
  case 20:
  {
    static Fluid3lin_Impl* f20;
    if (f20==NULL)
      f20 = new Fluid3lin_Impl(20);
    return f20;
  }
  case 27:
  {
    static Fluid3lin_Impl* f27;
    if (f27==NULL)
      f27 = new Fluid3lin_Impl(27);
    return f27;
  }
  case 4:
  {
    static Fluid3lin_Impl* f4;
    if (f4==NULL)
      f4 = new Fluid3lin_Impl(4);
    return f4;
  }
  case 10:
  {
    static Fluid3lin_Impl* f10;
    if (f10==NULL)
      f10 = new Fluid3lin_Impl(10);
    return f10;
  }
  case 6:
  {
    static Fluid3lin_Impl* f6;
    if (f6==NULL)
      f6 = new Fluid3lin_Impl(6);
    return f6;
  }
  case 15:
  {
    static Fluid3lin_Impl* f15;
    if (f15==NULL)
      f15 = new Fluid3lin_Impl(15);
    return f15;
  }
  case 5:
  {
    static Fluid3lin_Impl* f5;
    if (f5==NULL)
      f5 = new Fluid3lin_Impl(5);
    return f5;
  }

  default:
    dserror("node number %d not supported", f3->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3lin_Impl::Fluid3lin_Impl(int iel)
  : iel_(iel),
    xyze_(3,iel_,blitz::ColumnMajorArray<2>()),
    edeadng_(3,iel_,blitz::ColumnMajorArray<2>()),
    funct_(iel_),
    deriv_(3,iel_,blitz::ColumnMajorArray<2>()),
    deriv2_(6,iel_,blitz::ColumnMajorArray<2>()),
    xjm_(3,3,blitz::ColumnMajorArray<2>()),
    xji_(3,3,blitz::ColumnMajorArray<2>()),
    vderxy_(3,3,blitz::ColumnMajorArray<2>()),
    pderxy_(3),
    vderxy2_(3,6,blitz::ColumnMajorArray<2>()),
    derxy_(3,iel_,blitz::ColumnMajorArray<2>()),
    derxy2_(6,iel_,blitz::ColumnMajorArray<2>()),
    bodyforce_(3),
    histvec_(3),
    velino_(3),
    velint_(3),
    gradp_(3),
    tau_(3),
    viscs2_(3,3,iel_,blitz::ColumnMajorArray<3>()),
    conv_(iel_),
    rhsint_(3),
    xder2_(6,3,blitz::ColumnMajorArray<2>()),
    numepn_(iel_)
{
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                chfoe 02/08|
 *----------------------------------------------------------------------*/
/*
 Note: This routine is made for total rather than incremental soltution.

       THE PRESENT FLUID ELEMENT IS ENTIRELY LINEAR!!!!!!!!!!!!!!!
*/
void DRT::ELEMENTS::Fluid3lin_Impl::Sysmat(Fluid3* ele,
                                           const blitz::Array<double,2>&     evelnp,
                                           const blitz::Array<double,1>&     eprenp,
                                           const blitz::Array<double,2>&     evhist,
                                           blitz::Array<double,2>&           estif,
                                           blitz::Array<double,1>&           eforce,
                                           struct _MATERIAL*       material,
                                           double                  time,
                                           double                  timefac
  )
{

  // set element data
  const DRT::Element::DiscretizationType distype = ele->Shape();

  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];

    numepn_(inode) = nodes[inode]->NumElement();
  }


  // dead load in element nodes
  BodyForce(ele,time);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index

  blitz::Range _  = blitz::Range::all();
//   blitz::Range ux = blitz::Range(0, 4*iel_-4, 4);
//   blitz::Range uy = blitz::Range(1, 4*iel_-3, 4);
//   blitz::Range uz = blitz::Range(2, 4*iel_-2, 4);
//   blitz::Range p  = blitz::Range(3, 4*iel_-1, 4);

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  Caltau(ele,evelnp,distype,visc,timefac);

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
    xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                       xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                       xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                       xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                       xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                       xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
    const double fac = intpoints.qwgt[iquad]*det;

    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %lf", ele->Id(), det);
    }

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
    derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

    // compute second global derivative
    if (higher_order_ele)
    {
      DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      gder2(ele);

      // calculate 2nd velocity derivatives at integration point
      vderxy2_ = blitz::sum(derxy2_(j,k)*evelnp(i,k),k);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get velocities (n+g,i) at integration point
    velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

    // get history data (n,i) at integration point
    histvec_ = blitz::sum(funct_(j)*evhist(i,j),j);

    // get velocity (np,i) derivatives at integration point
    vderxy_ = blitz::sum(derxy_(j,k)*evelnp(i,k),k);

    // get pressure gradients
    gradp_ = blitz::sum(derxy_(i,j)*eprenp(j),j);

    // get bodyforce in gausspoint
    bodyforce_ = blitz::sum(edeadng_(i,j)*funct_(j),j);

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

    /*--------------------------------------- divergence of old velocity ---*/
    double divu_old = vderxy_(0,0) + vderxy_(1,1) + vderxy_(2,2);

    /*------------------------- evaluate rhs vector at integration point ---*/
    rhsint_ = histvec_(i) + bodyforce_(i)*timefac;

    /*----------------- get numerical representation of single operators ---*/

    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_ = blitz::sum(derxy_(j,i)*velint_(j), j);

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                                                \
         |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
       1 |                                                |
     - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
       2 |                                                |
         |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
         \                                                /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */

    viscs2_(0,0,_) = - 0.5 * (2.0 * derxy2_(0,_) + derxy2_(1,_) + derxy2_(2,_));
    viscs2_(0,1,_) = - 0.5 *  derxy2_(3,_);
    viscs2_(0,2,_) = - 0.5 *  derxy2_(4,_);
    viscs2_(1,0,_) = - 0.5 *  derxy2_(3,_);
    viscs2_(1,1,_) = - 0.5 * (derxy2_(0,_) + 2.0 * derxy2_(1,_) + derxy2_(2,_));
    viscs2_(1,2,_) = - 0.5 *  derxy2_(5,_);
    viscs2_(2,0,_) = - 0.5 *  derxy2_(4,_);
    viscs2_(2,1,_) = - 0.5 *  derxy2_(5,_);
    viscs2_(2,2,_) = - 0.5 * (derxy2_(0,_) + derxy2_(1,_) + 2.0 * derxy2_(2,_));

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--------------------------------- now build single stiffness terms ---*/

    //----------------------------------------------------------------------
    //                            GALERKIN PART

    for (int ui=0; ui<iel_; ++ui)
    {
      for (int vi=0; vi<iel_; ++vi)
      {
	/* an auxiliary variable for speedup */
        double aux;

	/* inertia (contribution to mass matrix) */
	/*
                               /       \
                              |  u , v  |
                               \       /
	*/
        aux = fac*funct_(ui)*funct_(vi);
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
        aux = 0.5*timefacfac*divu_old*funct_(ui)*funct_(vi);
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
        aux = timefacfac*funct_(vi);
	estif(vi*4 + 3, ui*4)     += aux * derxy_(0, ui) ;
	estif(vi*4 + 3, ui*4 + 1) += aux * derxy_(1, ui) ;
	estif(vi*4 + 3, ui*4 + 2) += aux * derxy_(2, ui) ;


	//-----------------------------------------------------------------
	//                           STABILISATION PART

	/* pressure stabilisation: inertia PLUS */
	/* pressure stabilisation: convection, convective part */
	/*
                        /            \     /                           \
                       |  u , grad q  | + |  u_old  o  grad u , grad q  |
                        \            /     \                           /
	*/
	aux = timetauMp*funct_(ui) + ttimetauMp*conv_(ui);
	estif(vi*4 + 3, ui*4)     += aux * derxy_(0, vi) ;
	estif(vi*4 + 3, ui*4 + 1) += aux * derxy_(1, vi) ;
	estif(vi*4 + 3, ui*4 + 2) += aux * derxy_(2, vi) ;

	/* pressure stabilisation: viscosity (-L_visc_u) */
	/*
                     /                     \
                    |  div eps(u) , grad q  |
                     \                     /
	*/
	aux = 2.0*visc*ttimetauMp;
	estif(vi*4 + 3, ui*4)     += aux*( derxy_(0, vi)*viscs2_(0, 0, ui)
					   +derxy_(1, vi)*viscs2_(0, 1, ui)
					   +derxy_(2, vi)*viscs2_(0, 2, ui));
	estif(vi*4 + 3, ui*4 + 1) += aux*( derxy_(0, vi)*viscs2_(0, 1, ui)
					   +derxy_(1, vi)*viscs2_(1, 1, ui)
					   +derxy_(2, vi)*viscs2_(1, 2, ui));
	estif(vi*4 + 3, ui*4 + 2) += aux*( derxy_(0, vi)*viscs2_(0, 2, ui)
					   +derxy_(1, vi)*viscs2_(1, 2, ui)
					   +derxy_(2, vi)*viscs2_(2, 2, ui));

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
        aux = timetauM*funct_(ui)*conv_(vi) + ttimetauM*conv_(ui)*conv_(vi);
	estif(vi*4, ui*4)         += aux;
	estif(vi*4 + 1, ui*4 + 1) += aux;
	estif(vi*4 + 2, ui*4 + 2) += aux;

	/* supg stabilisation: viscous part  (-L_visc_u) */
	/*
                  /                               \
                 |  div eps(u),  u_old  o  grad v  |
                  \                               /
	*/
	aux = 2.0*visc*ttimetauM*conv_(vi);
	estif(vi*4, ui*4)         += aux*viscs2_(0, 0, ui) ;
	estif(vi*4, ui*4 + 1)     += aux*viscs2_(0, 1, ui) ;
	estif(vi*4, ui*4 + 2)     += aux*viscs2_(0, 2, ui) ;
	estif(vi*4 + 1, ui*4)     += aux*viscs2_(0, 1, ui) ;
	estif(vi*4 + 1, ui*4 + 1) += aux*viscs2_(1, 1, ui) ;
	estif(vi*4 + 1, ui*4 + 2) += aux*viscs2_(1, 2, ui) ;
	estif(vi*4 + 2, ui*4)     += aux*viscs2_(0, 2, ui) ;
	estif(vi*4 + 2, ui*4 + 1) += aux*viscs2_(1, 2, ui) ;
	estif(vi*4 + 2, ui*4 + 2) += aux*viscs2_(2, 2, ui) ;

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
	estif(vi*4, ui*4)         += timetau_C*derxy_(0, ui)*derxy_(0, vi);
	estif(vi*4, ui*4 + 1)     += timetau_C*derxy_(0, vi)*derxy_(1, ui);
	estif(vi*4, ui*4 + 2)     += timetau_C*derxy_(0, vi)*derxy_(2, ui);
	estif(vi*4 + 1, ui*4)     += timetau_C*derxy_(0, ui)*derxy_(1, vi);
	estif(vi*4 + 1, ui*4 + 1) += timetau_C*derxy_(1, ui)*derxy_(1, vi);
	estif(vi*4 + 1, ui*4 + 2) += timetau_C*derxy_(1, vi)*derxy_(2, ui);
	estif(vi*4 + 2, ui*4)     += timetau_C*derxy_(0, ui)*derxy_(2, vi);
	estif(vi*4 + 2, ui*4 + 1) += timetau_C*derxy_(1, ui)*derxy_(2, vi);
	estif(vi*4 + 2, ui*4 + 2) += timetau_C*derxy_(2, ui)*derxy_(2, vi);
        }
      }

    /*------------------------------------ now build single rhs terms ---*/
    for (int vi=0; vi<iel_; ++vi)
    {
      double aux;
      //-------------------------------------------------------------------
      //                            GALERKIN PART

      // source term of the right hand side
      aux = fac*funct_(vi);
      eforce(vi*4)     += aux*rhsint_(0) ;
      eforce(vi*4 + 1) += aux*rhsint_(1) ;
      eforce(vi*4 + 2) += aux*rhsint_(2) ;

      //-------------------------------------------------------------------
      //                           STABILISATION PART

      // pressure stabilisation
      eforce(vi*4 + 3) += timetauMp*( rhsint_(0)*derxy_(0, vi)
				     +rhsint_(1)*derxy_(1, vi)
				     +rhsint_(2)*derxy_(2, vi)) ;
      // supg stabilisation
      aux =  timetauM*conv_(vi);
      eforce(vi*4)     += aux*rhsint_(0) ;
      eforce(vi*4 + 1) += aux*rhsint_(1) ;
      eforce(vi*4 + 2) += aux*rhsint_(2) ;
    } // vi
  }
}




//
// calculate stabilization parameter
//
void DRT::ELEMENTS::Fluid3lin_Impl::Caltau(
  Fluid3* ele,
  const blitz::Array<double,2>&           evelnp,
  const DRT::Element::DiscretizationType  distype,
  const double                            visc,
  const double                            timefac
  )
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

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
  velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

  // get Jacobian matrix and determinant
  xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  const double vol = wquad*det;

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);

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
  derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

  // get velocity norm
  const double vel_norm = sqrt(blitz::sum(velint_*velint_));

  // normed velocity at element centre
  if (vel_norm>=1e-6)
  {
    velino_ = velint_/vel_norm;
  }
  else
  {
    velino_ = 0.;
    velino_(0) = 1;
  }

  // get streamlength
  const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
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


  const double re1 =/* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
  const double re2 = mk * vel_norm * strle / /* *1.0 */(2.0 * visc);    /* convective : viscous forces */

  const double xi1 = DMAX(re1,1.0);
  const double xi2 = DMAX(re2,1.0);

  tau_(0) = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

  // compute tau_Mp
  //    stability parameter definition according to Franca and Valentin (2000)
  //                                       and Barrenechea and Valentin (2002)
  const double re_viscous = /* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
  const double re_convect = mk * vel_norm * hk / /* *1.0 */(2.0 * visc);     /* convective : viscous forces */

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
  tau_(1) = DSQR(hk) / (DSQR(hk) * xi_viscous + (/* 2.0*/ 4.0 * timefac * visc/mk) * xi_convect);

  /*------------------------------------------------------ compute tau_C ---*/
  /*-- stability parameter definition according to Codina (2002), CMAME 191
   *
   * Analysis of a stabilized finite element approximation of the transient
   * convection-diffusion-reaction equation using orthogonal subscales.
   * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
   *
   * */
  tau_(2) = sqrt(4.0*DSQR(visc)+DSQR(0.5*vel_norm*hk));

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
  //tau_(2) = vel_norm * hk * 0.5 * xi_tau_c;
}



// NOTE: this is an exact copy of the routine in fluid3_impl.cpp
/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3lin_Impl::BodyForce(Fluid3* ele, const double time)
{
  vector<DRT::Condition*> myneumcond;
  DRT::Node** nodes = ele->Nodes();

  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<iel_;inode++)
  {
    nodes[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  if (nodecount == iel_)
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

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;isd++)
      {
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_ = 0.;
  }
}

// NOTE: this is an exact copy of the routine in fluid3_impl.cpp
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
void DRT::ELEMENTS::Fluid3lin_Impl::gder2(Fluid3* ele)
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // initialize and zero out everything
  static Epetra_SerialDenseMatrix bm(6,6);

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

  xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);

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

  derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);

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

  derxy2_ += deriv2_;

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

  Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,iel_);

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(ederxy2,ederxy2);
  solver.Solve();

  return;
}


#endif
#endif
