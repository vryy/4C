/*!----------------------------------------------------------------------
\file condif2_evaluate.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "condif2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/convecdiffus.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 05/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Condif2::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::Elements::Condif2::ActionType act = Condif2::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_condif_systemmat_and_residual")
  	act = Condif2::calc_condif_systemmat_and_residual;
  else dserror("Unknown type of action for Condif2");

  // get the material
  RefCountPtr<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_condif)
    dserror("convection-diffusion material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::ConvecDiffus*>(mat.get())->MaterialData();

  switch(act)
  {
    case calc_condif_systemmat_and_residual:
    {
      // need current velocity and history vector
      RefCountPtr<const Epetra_Vector> hist = discretization.GetState("hist");
      if (hist==null) dserror("Cannot get state vector 'hist'");

      // extract local values from the global vector
      vector<double> myhist(lm.size());
      DRT::Utils::ExtractMyValues(*hist,myhist,lm);

      // get control parameter
      const bool is_stationary = params.get<bool>("using stationary formulation",false);
      const double time = params.get<double>("total time",-1.0);

      // One-step-Theta: timefac = theta*dt
      // BDF2:           timefac = 2/3 * dt
      double timefac = 0;
      if (not is_stationary)
      {
        timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");
      }

      // get type of velocity field
      const int vel_field = params.get<int>("condif velocity field",0);

      // get flag for discontinuity capturing (1=yes, 0=no)
      const int discap = params.get<int>("discontinuity capturing",0);

      // calculate element coefficient matrix and rhs
      condif2_sys_mat(lm,myhist,&elemat1,&elemat2,&elevec1,actmat,time,timefac,vel_field,discap,
                      is_stationary);

    }
    break;
    default:
      dserror("Unknown type of action for Condif2");
  } // end of switch(act)

  return 0;
} // end of DRT::Elements::Condif2::Evaluate



/*----------------------------------------------------------------------*
 |  do nothing (public)                                         vg 08/07|
 |                                                                      |
 |  The function is just a dummy. For the condif2 elements, the         |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilization terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::Condif2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs for convec.-diff. (private) vg 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::condif2_sys_mat(vector<int>&              lm,
                                         vector<double>&           ehist,
                                         Epetra_SerialDenseMatrix* sys_mat,
                                         Epetra_SerialDenseMatrix* sys_mat_dc,
                                         Epetra_SerialDenseVector* residual,
                                         struct _MATERIAL*         material,
                                         double                    time,
                                         double                    timefac,
                                         int                       vel_field,
                                         int                       discap,
                                         bool                      is_stationary)
{

  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();
  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // dead load in element nodes
  const Epetra_SerialDenseVector bodyforce = condif2_getbodyforce();

  /*---------------------------------------------- get diffusivity ---*/
  if(material->mattyp != m_condif) dserror("Material law is not of type m_condif.");
  const double  diffus = material->m.condif->diffusivity;

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector      funct(iel);
  Epetra_SerialDenseMatrix 	deriv(2,iel);
  Epetra_SerialDenseMatrix 	deriv2(3,iel);
  static Epetra_SerialDenseMatrix 	xjm(2,2);
  Epetra_SerialDenseMatrix 	derxy(2,iel);
  Epetra_SerialDenseMatrix 	derxy2(3,iel);
  static vector<double>         velint(2);
  double                        edeadng;
  double              		hist; /* history data at integration point      */
  double         		hk;  /* element length for calculation of tau  */
  double         		vel_norm, epe1, epe2, xi1, xi2;
  double         		mk=0.0;
  double                        tau; // stabilization parameter
  double                        kart; // artificial diffusivity

  /*----------------------------------------------------------------------*/
  // calculation of stabilization parameter
  /*----------------------------------------------------------------------*/
  /*------------------------------------------------------- initialize ---*/
    // use one point gauss rule to calculate tau at element center
  GaussRule2D integrationrule_stabili;
  switch(distype)
  {
  case quad4: case quad8: case quad9:
      integrationrule_stabili = intrule_quad_1point;
      break;
  case tri3: case tri6:
      integrationrule_stabili = intrule_tri_1point;
      break;
  default:
      dserror("invalid discretization type");
  }
  // gaussian points
  const IntegrationPoints2D  intpoints_tau = getIntegrationPoints2D(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints_tau.qxg[0][0];
  const double e2    = intpoints_tau.qxg[0][1];
  // shape functions and their derivatives
  DRT::Utils::shape_function_2D(funct,e1,e2,distype);
  DRT::Utils::shape_function_2D_deriv1(deriv,e1,e2,distype);

/*------------------------------- get element type constant for tau ---*/
  switch(iel)
  {
  case 3:
  case 4:
    mk = 0.333333333333333333333;
    break;
  case 6:
  case 8:
  case 9:
    mk = 0.083333333333333333333;
    break;
  default: dserror("type unknown!\n");
  }

/*--------------------------------- get velocities at element center ---*/
  switch(vel_field)
  {
  case 1:
    velint[0]=1.0;
    velint[1]=0.0;
    break;
  case 2:
    velint[0]=0.5*sqrt(3.0);
    velint[1]=0.5;
    break;
  case 3:
    velint[0]=0.5;
    velint[1]=0.5*sqrt(3.0);
    break;
  case 4:
    velint[0]=0.5;
    velint[1]=-0.5*sqrt(3.0);
    break;
  default: dserror("condif velocity field unknown!\n");
  }

/*------------------------------ get Jacobian matrix and determinant ---*/
  double  det;
  condif2_jaco(xyze,deriv,xjm,&det,iel);

/*----------------------------------------------- get element length ---*/
/*  the element length is chosen as the square root of the element area */
  {
    double area=0;
    double a,b,c;

    switch(iel)
    {
    case 3:
    case 6:
      a = (xyze(0,0)-xyze(0,1))*(xyze(0,0)-xyze(0,1))
          +(xyze(1,0)-xyze(1,1))*(xyze(1,0)-xyze(1,1)); /* line 0-1 squared */
      b = (xyze(0,1)-xyze(0,2))*(xyze(0,1)-xyze(0,2))
          +(xyze(1,1)-xyze(1,2))*(xyze(1,1)-xyze(1,2)); /* line 1-2 squared */
      c = (xyze(0,2)-xyze(0,0))*(xyze(0,2)-xyze(0,0))
          +(xyze(1,2)-xyze(1,0))*(xyze(1,2)-xyze(1,0)); /* diag 2-0 squared */
      area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      break;
    case 4:
    case 8:
    case 9:
    {
      a = (xyze(0,0)-xyze(0,1))*(xyze(0,0)-xyze(0,1))
          +(xyze(1,0)-xyze(1,1))*(xyze(1,0)-xyze(1,1)); /* line 0-1 squared */
      b = (xyze(0,1)-xyze(0,2))*(xyze(0,1)-xyze(0,2))
          +(xyze(1,1)-xyze(1,2))*(xyze(1,1)-xyze(1,2)); /* line 1-2 squared */
      c = (xyze(0,2)-xyze(0,0))*(xyze(0,2)-xyze(0,0))
          +(xyze(1,2)-xyze(1,0))*(xyze(1,2)-xyze(1,0)); /* diag 2-0 squared */
      area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      a = (xyze(0,2)-xyze(0,3))*(xyze(0,2)-xyze(0,3))
          +(xyze(1,2)-xyze(1,3))*(xyze(1,2)-xyze(1,3)); /* line 2-3 squared */
      b = (xyze(0,3)-xyze(0,0))*(xyze(0,3)-xyze(0,0))
          +(xyze(1,3)-xyze(1,0))*(xyze(1,3)-xyze(1,0)); /* line 3-0 squared */
      /*-------------------------------- evaluate element area ---*/
      area += 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
      break;
    }
    default: dserror("type unknown!\n");
    }

    hk = sqrt(area);
  }

  /*----------------------------------------------- get vel_norm ---*/
  vel_norm = sqrt(DSQR(velint[0]) + DSQR(velint[1]));

  if (is_stationary == false)
  {// stabilization parameters for instationary case (default)

    /* parameter relating diffusive : reactive forces */
    epe1 = 2.0 * timefac * diffus / (mk * DSQR(hk));
    /* parameter relating convective : diffusive forces */
    epe2 = mk * vel_norm * hk / diffus;
    xi1 = DMAX(epe1,1.0);
    xi2 = DMAX(epe2,1.0);

    /*--------------------------------------------------- compute tau ---*/
    tau = DSQR(hk)/((DSQR(hk)*xi1)/timefac + (2.0*diffus/mk)*xi2);

  }
  else
  {// stabilization parameters for stationary case

    /*------------------------------------------------------ compute tau ---*/
    /* stability parameter definition according to Franca and Valentin (2000) */
    epe1 = mk * vel_norm * hk / diffus;      /* convective : diffusive forces */
    xi1 = DMAX(epe1,1.0);

    tau = (DSQR(hk)*mk)/(2.0*diffus*xi1);

  }
  if (discap == 1)
  {
    /*-------------------------- compute artificial diffusivity kappa_art ---*/
    epe1 = mk * vel_norm * hk / diffus;     /* convective : diffusive forces */
    xi1 = DMAX(epe1,1.0);

    kart = (DSQR(hk)*mk*DSQR(vel_norm))/(2.0*diffus*xi1);
  }

  /*----------------------------------------------------------------------*/
  // integration loop for one condif2 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = is_higher_order_element(distype);

  // gaussian points
  const GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule);

  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
      const double e1 = intpoints.qxg[iquad][0];
      const double e2 = intpoints.qxg[iquad][1];

      // shape functions and their derivatives
      shape_function_2D(funct,e1,e2,distype);
      shape_function_2D_deriv1(deriv,e1,e2,distype);
      if (higher_order_ele)
      {
         shape_function_2D_deriv2(deriv2,e1,e2,distype);
      }

      /*----------------------------------------- compute Jacobian matrix */
      condif2_jaco(xyze,deriv,xjm,&det,iel);
      const double fac = intpoints.qwgt[iquad]*det;

      /*---------------------------------------- compute global derivates */
      condif2_gder(derxy,deriv,xjm,det,iel);

      /*--------------------------------- compute second global derivative */
      if (higher_order_ele) condif2_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

      /*---------------------- get velocity at integration point */
      switch(vel_field)
      {
      case 1:
        velint[0]=1.0;
        velint[1]=0.0;
        break;
      case 2:
        velint[0]=0.5*sqrt(3.0);
        velint[1]=0.5;
        break;
      case 3:
        velint[0]=0.5;
        velint[1]=0.5*sqrt(3.0);
        break;
      case 4:
        velint[0]=0.5;
        velint[1]=-0.5*sqrt(3.0);
        break;
      default: dserror("condif velocity field unknown!\n");
      }

      /*---------------- get history data (n,i) at integration point */
      hist=ZERO;
      for (int j=0;j<iel;j++)
      {
        hist += funct[j]*ehist[j];
      } /* end of loop over j */

      // get bodyforce in gausspoint
      edeadng = 0.0;
      for (int inode=0;inode<iel;inode++)
      {
        edeadng+= bodyforce[inode]*funct[inode];
      }

      /*-------------- perform integration for entire matrix and rhs ---*/
      if(is_stationary==false)
        condif2_calmat(*sys_mat,*sys_mat_dc,*residual,velint,hist,funct,derxy,derxy2,
                       edeadng,tau,kart,fac,diffus,iel,discap,timefac);
      else
        condif2_calmat_stat(*sys_mat,*sys_mat_dc,*residual,velint,hist,funct,derxy,
                            derxy2,edeadng,tau,kart,fac,diffus,iel,discap);

  } // end of loop over integration points

    return;
} // DRT::Elements::Condif2::condif2_sys_mat


/*----------------------------------------------------------------------*
 | calculate Jacobian matrix and it's determinant (private)     vg 05/07|
 |
 |
 |     +-        -+ T
 |     | dx    dx |
 |     | ---   -- |
 |     | dr    ds |
 |     |	  |
 |     | dy    dy |
 |     | ---   -- |
 |     | dr    ds |
 |     +-        -+
 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::condif2_jaco(const Epetra_SerialDenseMatrix& xyze,
				    const Epetra_SerialDenseMatrix& deriv,
                                    Epetra_SerialDenseMatrix& xjm,
				    double* det,
                                    const int iel
				    )
{
  double dum;

  /*-------------------------------- determine jacobian at point r,s ---*/
  for (int i=0; i<2; i++)
  {
     for (int j=0; j<2; j++)
     {
        dum=0.0;
        for (int l=0; l<iel; l++)
        {
           dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
     } /* end of loop j */
  } /* end of loop i */

  /*------------------------------------------ determinant of jacobian---*/
  *det = xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1);

  if(*det<0.0)
  {
     printf("\n");
     printf("GLOBAL ELEMENT NO.%i\n",Id());
     printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n",*det);
     dserror("Stopped not regulary!\n");
  }

} //end of DRT::Elements::Condif2::condif2_jaco


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private)    vg 08/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a surface Neumann condition    |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector DRT::Elements::Condif2::condif2_getbodyforce()
{
  const int iel = NumNode();
  Epetra_SerialDenseVector edeadng(iel);

  // we have only a constant dead load so far
  for(int inode=0;inode<iel;inode++)
  {
    edeadng(inode)=1.0;
  }

  return edeadng;
} // end of DRT:Elements:Condif2:condif2_getbodyforce


/*----------------------------------------------------------------------*
 |  calculate global derivatives w.r.t. x,y at point r,s (private)
 |                                                            vg 05/07
 |
 |      +-  -+      +-        -+     +-  -+
 |      | dN |      | dr    ds |     | dN |
 |      | -- |      | ---   -- |     | -- |
 |      | dx |      | dx    dx |     | dr |
 |      |    |  =   |	       |  *  |    |
 |      | dN |      | dr    ds |     | dN |
 |      | -- |      | ---   -- |     | -- |
 |      | dy |      | dy    dy |     | ds |
 |      +-  -+      +-	      -+     +-  -+
 |                        |
 |                        |
 |                        |
 |                      J^{-T}
 |
 | To calculate the derivatives, the actual Jacobian matrix is
 | inverted.
 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::condif2_gder(Epetra_SerialDenseMatrix& derxy,
				    const Epetra_SerialDenseMatrix& deriv,
                                    Epetra_SerialDenseMatrix& xjm,
				    double& det,
                                    const int iel
				    )
{
  static Epetra_SerialDenseMatrix       xji(2,2);   // inverse of jacobian matrix


  /*----------calculate global derivatives w.r.t. x,y at point r,s ---*/

  /*------------------------------------------------------- initialistion */
  for(int k=0;k<iel;k++)
  {
    derxy(0,k)=0.0;
    derxy(1,k)=0.0;
  } /* end of loop over k */


  /*------------------------------------------------- inverse of jacobian */
  xji(0,0) =  xjm(1,1)/det;
  xji(1,0) = -xjm(1,0)/det;
  xji(0,1) = -xjm(0,1)/det;
  xji(1,1) =  xjm(0,0)/det;

  /*---------------------------------------- calculate global derivatives */
  for (int k=0;k<iel;k++)
  {
    derxy(0,k) +=   xji(0,0) * deriv(0,k) + xji(0,1) * deriv(1,k) ;
    derxy(1,k) +=   xji(1,0) * deriv(0,k) + xji(1,1) * deriv(1,k) ;
  } /* end of loop over k */

  /*----------------------------------------------------------------------*/

  return;
} // end of DRT:Elements:Condif2:condif2_gder

/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y at point r,s (private)
 |                                                             vg 05/07
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
 | the matrix system
 |
 | +-                                        -+   +-    -+
 | |   /dx\^2        /dy\^2         dy dx     |	  | d^2N |
 | |  | -- |        | ---|        2*--*--     |	  | ---- |
 | |   \dr/	     \dr/ 	    dr dr     |	  | dx^2 |
 | |					      |	  |      |
 | |   /dx\^2        /dy\^2         dy dx     |	  | d^2N |
 | |  | -- |        | -- |        2*--*--     |	* | ---- |
 | |   \ds/	     \ds/ 	    ds ds     |   | dy^2 | =
 | |  					      |	  |      |
 | |   dx dx         dy dy      dx dy   dy dx |	  | d^2N |
 | |   --*--         --*--      --*-- + --*-- |   | ---- |
 | |   dr ds	     dr ds	dr ds   dr ds |	  | dxdy |
 | +-					     -+	  +-    -+
 |
 |             +-    -+   +-                 -+
 | 	       | d^2N |	  | d^2x dN   d^2y dN |
 | 	       | ---- |	  | ----*-- + ----*-- |
 |	       | dr^2 |	  | dr^2 dx   dr^2 dy |
 |	       |      |	  |                   |
 |	       | d^2N |	  | d^2x dN   d^2y dN |
 |          =  | ---- | - | ----*-- + ----*-- |
 |	       | ds^2 |	  | ds^2 dx   ds^2 dy |
 |	       |      |	  |                   |
 |	       | d^2N |	  | d^2x dN   d^2y dN |
 |	       | ---- |	  | ----*-- + ----*-- |
 |	       | drds |	  | drds dx   drds dy |
 |	       +-    -+	  +-                 -+
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

void DRT::Elements::Condif2::condif2_gder2(const Epetra_SerialDenseMatrix& xyze,
				     const Epetra_SerialDenseMatrix& xjm,
				     const Epetra_SerialDenseMatrix& derxy,
				     Epetra_SerialDenseMatrix& derxy2,
				     const Epetra_SerialDenseMatrix& deriv2,
				     const int iel
				    )
{
//--------------------------------------------initialize and zero out everything
    static Epetra_SerialDenseMatrix bm(3,3);
    static Epetra_SerialDenseMatrix xder2(3,2);
    Epetra_SerialDenseMatrix chainrulerhs(3,iel);

/*--------------------------- calculate elements of jacobian_bar matrix */
    bm(0,0) =                   xjm(0,0)*xjm(0,0);
    bm(0,1) =                   xjm(0,1)*xjm(0,1);
    bm(0,2) =               TWO*xjm(0,0)*xjm(0,1);

    bm(1,0) =                   xjm(1,0)*xjm(1,0);
    bm(1,1) =                   xjm(1,1)*xjm(1,1);
    bm(1,2) =               TWO*xjm(1,1)*xjm(1,0);

    bm(2,0) =                   xjm(0,0)*xjm(1,0);
    bm(2,1) =                   xjm(0,1)*xjm(1,1);
    bm(2,2) = xjm(0,0)*xjm(1,1)+xjm(0,1)*xjm(1,0);

    //init sol to zero
    memset(derxy2.A(),0,derxy2.M()*derxy2.N()*sizeof(double));


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
  |		     	      	     	       +-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-           -+
  |  	   	    	    	        | d^2x   d^2y |
  |  	   	    	    	        | ----   ---- |
  | 	   	   	   	        | dr^2   dr^2 |
  | 	   	   	   	        |             |
  | 	   	   	   	        | d^2x   d^2y |
  |                 yields    xder2  =  | ----   ---- |
  | 	   	   	   	        | ds^2   ds^2 |
  | 	   	   	   	        |             |
  | 	   	   	   	        | d^2x   d^2y |
  | 	   	   	   	        | ----   ---- |
  | 	   	   	   	        | drds   drds |
  | 	   	   	   	        +-           -+
  |
  |
  */
    xder2.Multiply('N','T',1.0,deriv2,xyze,0.0);

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
    xder2.Multiply(false,derxy,chainrulerhs);
    chainrulerhs.Scale(-1.0);

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
    chainrulerhs+=deriv2;


  /*
  |
  |          0  1  2         i        i
  | 	   +--+--+--+       +-+      +-+
  | 	   |  |  |  | 0     | | 0    | | 0
  | 	   +--+--+--+       +-+	     +-+
  | 	   |  |  |  | 1  *  | | 1 =  | | 1  for i=0...iel-1
  | 	   +--+--+--+       +-+	     +-+
  | 	   |  |  |  | 2     | | 2    | | 2
  | 	   +--+--+--+       +-+	     +-+
  |                          |        |
  |                          |        |
  |                        derxy2[i]  |
  |                                   |
  |                              chainrulerhs[i]
  |
  |
  |
  |                   0...iel-1
  |		     +-+-+-+-+
  |		     | | | | | 0
  |		     +-+-+-+-+
  |	  yields     | | | | | 1
  |		     +-+-+-+-+
  |                  | | | | | 2
  | 		     +-+-+-+-+
  |
  |                    derxy2
  |
  */

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix (bm);
    solver.SetVectors(derxy2,chainrulerhs);
    solver.Solve();
/*----------------------------------------------------------------------*/

    return;
} // end of DRT:Elements:Condif2:condif2_gder2


/*----------------------------------------------------------------------*
 |  evaluate instationary convection-diffusion matrix (private) vg 05/07|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif2 element are calculated for the instationary
case. The procedure is based on the Rothe method of first integrating in
time. Hence the resulting terms include coefficients containing time
integration variables such as theta or delta t which are represented by
'timefac'.

The stabilization is based on the residuum:

R = phi + timefac * u * grad(phi) - timefac * diffus * laplace(phi) - rhsint

The corresponding weighting operators are
L = timefac * u * grad(w) +/- timefac * diffus *  laplace(w)

'+': USFEM
'-': GLS


time-integration schemes:

one-step-theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old

BDF2:

generalised alpha:


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *hist       DOUBLE        (i)   rhs at INT point
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param   edeadng    DOUBLE        (i)   dead load
\param   tau        DOUBLE        (i)   stabilization parameter
\param   fac        DOUBLE        (i)   weighting factor
\param   diffus     DOUBLE        (i)   diffusivity
\param   iel        INT           (i)   number of nodes of act. ele
\return void
------------------------------------------------------------------------*/

void DRT::Elements::Condif2::condif2_calmat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseMatrix& edc,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    const double&             hist,
    Epetra_SerialDenseVector& funct,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const double&             edeadng,
    const double&             tau,
    const double&             kart,
    const double&             fac,
    const double&             diffus,
    const int&                iel,
    const int&                discap,
    double                    timefac
    )
{
/*========================= further variables =========================*/

vector<double>            conv(iel); 	    /* convective part       */
vector<double>            diff(iel); 	    /* diffusive part        */
static double             rhsint;           /* rhs at int. point     */

// stabilization parameter
const double taufac = tau*fac;

// integration factors and coefficients of single terms
const double timefacfac  = timefac * fac;
const double timetaufac  = timefac * taufac;

/*-------------------------------- evaluate rhs at integration point ---*/
rhsint = hist + edeadng*timefac;

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* convective part */
   /* u_x * N,x  +  u_y * N,y          with  N .. form function matrix */
   conv[i] = velint[0] * derxy(0,i) + velint[1] * derxy(1,i) ;

   /* diffusive part */
   /* diffus * ( N,xx  +  N,yy  ) */
   diff[i] = diffus * ( derxy2(0,i) + derxy2(1,i) );

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define derxy_(i,j)    derxy(i,j)
#define velint_(j)     velint[j]
#define conv_(j)       conv[j]
#define diff_(j)       diff[j]
#define rhsint_        rhsint
#define diffus_        diffus

#include "condif2_stiff.cpp"
#include "condif2_rhs.cpp"

#undef estif_
#undef eforce_
#undef funct_
#undef derxy_
#undef velint_
#undef conv_
#undef diff_
#undef rhsint_
#undef diffus_

if (discap == 1)
{
  // parameter for artificial diffusivity
  const double kartfac = kart*timefacfac;
  const double taumfac = tau*timefacfac;

  #define edc_(i,j)    edc(i,j)
  #define derxy_(i,j)  derxy(i,j)
  #define conv_(j)       conv[j]

  #include "condif2_kart.cpp"

  #undef edc_
  #undef derxy_
  #undef conv_
}

return;
} // end of DRT:Elements:Condif2:condif2_calmat


/*----------------------------------------------------------------------*
 |  evaluate stationary convection-diffusion matrix (private)   vg 05/07|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif2 element are calculated for the stationary
case.

The stabilization is based on the residuum:

R = u * grad(phi) - diffus *  laplace(phi) - rhsint

The corresponding weighting operators are
L = u * grad(w) +/- diffus *  laplace(w)

'+': USFEM
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *hist       DOUBLE        (i)   rhs at INT point
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param   edeadng    DOUBLE        (i)   dead load
\param   tau        DOUBLE        (i)   stabilization parameter
\param   fac        DOUBLE        (i)   weighting factor
\param   diffus     DOUBLE        (i)   diffusivity
\param   iel        INT           (i)   number of nodes of act. ele
\return void
------------------------------------------------------------------------*/

void DRT::Elements::Condif2::condif2_calmat_stat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseMatrix& edc,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    const double&             hist,
    Epetra_SerialDenseVector& funct,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const double&             edeadng,
    const double&             tau,
    const double&             kart,
    const double&             fac,
    const double&             diffus,
    const int&                iel,
    const int&                discap
    )
{
/*========================= further variables =========================*/

vector<double>            conv(iel); 	    /* convective part       */
vector<double>            diff(iel); 	    /* diffusive part        */
static double             rhsint;           /* rhs at int. point     */

// stabilization parameter
const double taufac = tau*fac;

/*-------------------------------- evaluate rhs at integration point ---*/
rhsint = hist + edeadng;

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* convective part */
   /* u_x * N,x  +  u_y * N,y          with  N .. form function matrix */
   conv[i] = velint[0] * derxy(0,i) + velint[1] * derxy(1,i) ;

   /* diffusive part */
   /* diffus * ( N,xx  +  N,yy  ) */
   diff[i] = diffus * ( derxy2(0,i) + derxy2(1,i) );

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define derxy_(i,j)    derxy(i,j)
#define velint_(j)     velint[j]
#define conv_(j)       conv[j]
#define diff_(j)       diff[j]
#define rhsint_        rhsint
#define diffus_        diffus

#include "condif2_stiff_stat.cpp"
#include "condif2_rhs_stat.cpp"

#undef estif_
#undef eforce_
#undef funct_
#undef derxy_
#undef velint_
#undef conv_
#undef diff_
#undef rhsint_
#undef diffus_

if (discap == 1)
{
  // parameter for artificial diffusivity
  const double kartfac = kart*fac;
  const double taumfac = tau*fac;

  #define edc_(i,j)    edc(i,j)
  #define derxy_(i,j)  derxy(i,j)
  #define conv_(j)       conv[j]

  #include "condif2_kart.cpp"

  #undef edc_
  #undef derxy_
  #undef conv_
}

return;
} // end of DRT:Elements:Condif2:condif2_calmat_stat

// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::Elements::Condif2::is_higher_order_element(
              const DRT::Element::DiscretizationType  distype) const
{
    bool hoel = true;
    switch (distype)
    {
    case quad4: case quad8: case quad9: case tri6:
        hoel = true;
        break;
    case tri3:
        hoel = false;
        break;
    default:
        dserror("distype unknown!");
    }
    return hoel;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Condif2Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
