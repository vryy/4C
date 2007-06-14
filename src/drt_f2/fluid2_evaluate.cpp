/*!----------------------------------------------------------------------
\file fluid2_evaluate.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "fluid2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"

extern "C"
{
#include "../headers/standardtypes.h"
//#include "../fluid2/fluid2.h"
}


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid2::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::Elements::Fluid2::ActionType act = Fluid2::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_fluid_systemmat_and_residual")
  	act = Fluid2::calc_fluid_systemmat_and_residual;
  else dserror("Unknown type of action for Fluid2");

  // get the material
  MATERIAL* actmat = &(mat[material_-1]);

  switch(act)
  {
    case calc_fluid_systemmat_and_residual:
    {
      // need current velocity and history vector
      RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (trial)");
      RefCountPtr<const Epetra_Vector> hist = discretization.GetState("old solution data for rhs");
      if (vel_pre_np==null || hist==null) dserror("Cannot get state vectors 'velnp' and/or 'hist'");

      // extract local values from the global vectors
      vector<double> my_vel_pre_np(lm.size());
      DRT::Utils::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);
      vector<double> myhist(lm.size());
      DRT::Utils::ExtractMyValues(*hist,myhist,lm);

      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;
      RefCountPtr<const Epetra_Vector> gridv;
      vector<double> mygridv;

      if (is_ale_)
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
        mydispnp.resize(lm.size());
        DRT::Utils::ExtractMyValues(*dispnp,mydispnp,lm);

        gridv = discretization.GetState("gridv");
        if (gridv==null) dserror("Cannot get state vectors 'gridv'");
        mygridv.resize(lm.size());
        DRT::Utils::ExtractMyValues(*gridv,mygridv,lm);
      }

      // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
      // Additionally only the velocity components of myhist are important!
      int numnode = NumNode();
      vector<double> myprenp(numnode);
      vector<double> myvelnp(2*numnode);
      vector<double> myvhist(2*numnode);

      for (int i=0;i<numnode;++i)
       {
       myvelnp[0+(i*2)]=my_vel_pre_np[0+(i*3)];
       myvelnp[1+(i*2)]=my_vel_pre_np[1+(i*3)];

       myprenp[i]=my_vel_pre_np[2+(i*3)];

       myvhist[0+(i*2)]=myhist[0+(i*3)];
       myvhist[1+(i*2)]=myhist[1+(i*3)];
       }

      // calculate element coefficient matrix and rhs
      f2_sys_mat(lm,myvelnp,myprenp,myvhist,mydispnp,mygridv,&elemat1,&elevec1,actmat,params);

      // This is a very poor way to transport the density to the
      // outside world. Is there a better one?
      params.set("density", actmat->m.fluid->density);

#if 0

	for (int i=0;i<elevec1.size();++i)
	    {
	    printf("eforce[%d]: %26.16e\n",i,elevec1[i]);
	    ;
	    }
	    printf("\n");
#endif
#if 0
        //if (Id()==0)
	    for (int i=0;i<elemat1.ColDim();++i)
	{
	    for (int j=0;j<elemat1.RowDim();++j)
	    {
		printf("%26.16e\n",elemat1(i,j));
//		printf("%3d res %26.19e\n",Id(),elevec1[i]);

	    }
	    printf("\n");
	}
#endif

#if 0
        for (unsigned int i=0;i<myvelnp.size();++i){
	    printf("vel %26.16e ",myvelnp[i]);
	    printf("\n");
	    }
#endif


    }
    break;
    default:
      dserror("Unknown type of action for Fluid2");
  } // end of switch(act)

  return 0;
} // end of DRT::Elements::Fluid2::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)   gammi 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  dserror("Volume Neumann conditions not yet implemented for Fluid2");
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                gammi 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid2::f2_sys_mat(vector<int>&              lm,
                                       vector<double>&           evelnp,
				       vector<double>&           eprenp,
                                       vector<double>&           evhist,
                                       vector<double>&           edispnp,
                                       vector<double>&           egridv,
                                       Epetra_SerialDenseMatrix* sys_mat,
                                       Epetra_SerialDenseVector* residual,
                                       struct _MATERIAL*         material,
				       ParameterList& 		 params
  )
{

  /*---------------------------------------------------- set element data */
  const int iel = NumNode();
  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  if (is_ale_)
  {
    for (int i=0;i<iel;i++)
    {
      xyze(0,i) += edispnp[3*i];
      xyze(1,i) += edispnp[3*i+1];
    }
  }

  /*---------------------------------------------- get viscosity ---*/
  // check here, if we really have a fluid !!
  if(material->mattyp != m_fluid) dserror("Material law is not of type m_fluid.");
  const double  visc = material->m.fluid->viscosity;

  /*--------------------------------------------- stab-parameter ---*/
  // USFEM stabilization is default. No switch here at the moment.

  /*----------------------------------------- declaration of variables ---*/
  vector<double> 		funct(iel);
  Epetra_SerialDenseMatrix 	deriv(2,iel);
  Epetra_SerialDenseMatrix 	deriv2(3,iel);
  Epetra_SerialDenseMatrix 	xjm(2,2);
  Epetra_SerialDenseMatrix 	vderxy(2,2);
  vector<double> 		pderxy(2);
  Epetra_SerialDenseMatrix 	vderxy2(2,3);
  Epetra_SerialDenseMatrix 	derxy(2,iel);
  Epetra_SerialDenseMatrix 	derxy2(3,iel);
  vector<double> 		ephin(iel);
  vector<double> 		ephing(iel);
  vector<double> 		iedgnod(iel);
  Epetra_SerialDenseMatrix 	wa1(100,100);  // working matrix used as dummy
  vector<double>    		histvec(2); /* history data at integration point      */
  double         		hk;         /* element length for calculation of tau  */
  double         		vel_norm, pe, re, xi1, xi2, xi;
  vector<double>         	velino(2); /* normed velocity at element centre */
  double         		det;
  FLUID_DATA            	data;
  double         		e1, e2;
  double         		facr=0.0, facs=0.0;
  double         		mk=0.0;
  vector<double>     		velint(2);
  double 			timefac;
  vector<double>               tau(3); // stab parameters

  /*------------------------------------------------------- initialise ---*/
  // gaussian points
  f2_integration_points(data);
  timefac=params.get<double>("time constant for integration",0.0);


  /*---------------------- shape functions and derivs at element center --*/
  switch(iel)
  {
  case 4: case 8: case 9:   /* --> quad - element */
    e1   = data.qxg [0][0];
    facr = data.qwgt[0][0];
    e2   = data.qxg [0][0];
    facs = data.qwgt[0][0];

    f2_shape_function(funct,deriv,wa1,e1,e2,iel,2); //wa1 as dummy for not wanted second derivatives
    break;
  case 3: case 6:   /* --> tri - element */
    e1   = data.txgr[0][0];
    e2   = data.txgs[0][0];

    facr = data.twgt[0][0];
    facs = ONE;

    f2_shape_function(funct,deriv,wa1,e1,e2,iel,2); //wa1 as dummy for not wanted second derivatives
    break;
  default:
    dserror("type unknown!\n");
  } /*end switch(iel) */

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
  for (int i=0;i<2;i++)
  {
    velint[i]=ZERO;
    for (int j=0;j<iel;j++)
    {
      velint[i] += funct[j]*evelnp[i+(2*j)];
    }
  } //end loop over i

/*------------------------------ get Jacobian matrix and determinant ---*/
  f2_jaco(xyze,deriv,xjm,&det,iel);

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

  // get control parameter
  bool is_stationary = params.get<bool>("using stationary formulation",false);

  if (is_stationary == false)
  {// stabilization parameters for instationary case (default)

    /* parameter relating viscous : reactive forces */
    pe = 4.0 * timefac * visc / (mk * DSQR(hk));
    /* parameter relating advective : viscous forces */
    re = mk * vel_norm * hk / (2.0 * visc);

    xi1 = DMAX(pe,1.0);
    xi2 = DMAX(re,1.0);


    /*-------------------------------------------- compute tau_M ---*/
    tau[0] = DSQR(hk)/(DSQR(hk)*xi1 + (4.0*timefac*visc/mk)*xi2);

    /*-------------------------------------------- compute tau_C ---*/
    xi2 = DMIN(re,1.0);
    tau[2] = vel_norm * hk * 0.5 * xi2 /timefac;
  }
  else
  {// stabilization parameters for stationary case

    /*------------------------------------------------------ compute tau_M ---*/
    /* stability parameter definition according to Franca and Valentin (2000) */
    re = mk * vel_norm * hk / (2.0 * visc);       /* convective : viscous forces */

    xi = DMAX(re,1.0);

    tau[0] = (DSQR(hk)*mk)/(4.0*visc*xi);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191 */
    xi = DMIN(re,1.0);
    tau[2] = 0.5*vel_norm*hk*xi;
  }



/*----------------------------------------------------------------------*/
// end of f2_caltau function
/*----------------------------------------------------------------------*/

// integration loop for one Fluid2 element using USFEM

  int       intc=0;      /* "integration case" for tri for further infos
                            see f2_inpele.c and f2_intg.c                 */
  int       nir=0;       /* number of integration nodesin r direction     */
  int       nis=0;       /* number of integration nodesin s direction     */
  int       ihoel=0;     /* flag for higher order elements                 */
  int       icode=2;     /* flag for eveluation of shape functions         */
  double    fac;         /* total integration factor */
  double    press;
  vector<double>    gridvelint(2); /* grid velocity                       */
  vector<double>    gradp(2);      /* pressure gradient at integration point         */

  switch (iel)
  {
  case 4: case 8: case 9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ngp_[0];
    nis = ngp_[1];
    break;
  case 6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case 3:
    /* initialise integration */
    nir  = ngp_[0];
    nis  = 1;
    intc = ngp_[1];
    break;
  default:
    dserror("typ unknown!");
  } // end switch (iel) //


/*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int lr=0;lr<nir;lr++)
  {
    for (int ls=0;ls<nis;ls++)
    {
/*------------- get values of  shape functions and their derivatives ---*/
      switch(iel)
      {
      case 4: case 8: case 9:   /* --> quad - element */
        e1   = data.qxg [lr][nir-1];
        facr = data.qwgt[lr][nir-1];
        e2   = data.qxg [ls][nis-1];
        facs = data.qwgt[ls][nis-1];
        f2_shape_function(funct,deriv,deriv2,e1,e2,iel,icode);
        break;
      case 3: case 6:   /* --> tri - element */
        e1   = data.txgr[lr][intc];
        facr = data.twgt[lr][intc];
        e2   = data.txgs[lr][intc];
        facs = ONE;
        f2_shape_function(funct,deriv,deriv2,e1,e2,iel,icode);
        break;
      default:
        facr = facs = 0.0;
        e1 = e2 = 0.0;
        dserror("typ unknown!");
      } /* end switch (iel) */

      /*----------------------------------------- compute Jacobian matrix */

      f2_jaco(xyze,deriv,xjm,&det,iel);
      fac = facr*facs*det;

      /*---------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

      /*--------------------------------- compute second global derivative */
      if (ihoel!=0)
      {
        f2_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

        /*------calculate 2nd velocity derivatives at integration point */
        // former f2_vder2(vderxy2,derxy2,evelnp,iel);
        for (int i=0;i<3;i++)
        {
          vderxy2(0,i)=ZERO;
          vderxy2(1,i)=ZERO;
          for (int j=0;j<iel;j++)
          {
            vderxy2(0,i) += derxy2(i,j)*evelnp[0+(2*j)];
            vderxy2(1,i) += derxy2(i,j)*evelnp[1+(2*j)];
          } /* end of loop over j */
        } /* end of loop over i */
      }

      /*---------------------- get velocities (n+g,i) at integration point */
      // expression for f2_veci(velint,funct,evelnp,iel);
      for (int i=0;i<2;i++)
      {
        velint[i]=ZERO;
        for (int j=0;j<iel;j++)
        {
          velint[i] += funct[j]*evelnp[i+(2*j)];
        }
      } //end loop over i

      /*---------------- get history data (n,i) at integration point ---*/
      //expression for f2_veci(histvec,funct,evhist,iel);
      for (int i=0;i<2;i++)
      {
        histvec[i]=ZERO;
        for (int j=0;j<iel;j++)
        {
          histvec[i] += funct[j]*evhist[i+(2*j)];
        } /* end of loop over j */
      } /* end of loop over i */

      /*----------- get velocity (np,i) derivatives at integration point */
      // expression for f2_vder(vderxy,derxy,evelnp,iel);
      for (int i=0;i<2;i++)
      {
        vderxy(0,i)=ZERO;
        vderxy(1,i)=ZERO;
        for (int j=0;j<iel;j++)
        {
          vderxy(0,i) += derxy(i,j)*evelnp[0+(2*j)];
          vderxy(1,i) += derxy(i,j)*evelnp[1+(2*j)];
        } /* end of loop over j */
      } /* end of loop over i */

      /*--------------------- get grid velocity at integration point ---*/
      if (is_ale_)
      {
        for (int i=0; i<2; i++)
        {
          gridvelint[i] = 0.;
          for (int j=0; j<iel; j++)
          {
            gridvelint[i] += derxy(i,j)*egridv[i+(3*j)];
          }
        }
      }
      else
      {
        gridvelint[0] = 0.0;
        gridvelint[1] = 0.0;
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp[0] = gradp[1] = 0.0;

      for (int i=0; i<iel; i++)
      {
        gradp[0] += derxy(0,i) * eprenp[i];
        gradp[1] += derxy(1,i) * eprenp[i];
      }

      press = 0;
      for (int i=0;i<iel;i++)
      {
        press += funct[i]*eprenp[i];
      }

      /*-------------- perform integration for entire matrix and rhs ---*/
      if(is_stationary==false)
        f2_calmat(*sys_mat,*residual,
                  velint,histvec,gridvelint,press,
                  vderxy,vderxy2,gradp,
                  funct,tau,
                  derxy,derxy2,
                  fac,visc,iel,
                  params);
      else
        f2_calmat_stationary(*sys_mat,*residual,
                             velint,histvec,gridvelint,press,
                             vderxy,vderxy2,gradp,
                             funct,tau,
                             derxy,derxy2,
                             fac,visc,iel,
                             params);


    } /* end of loop over integration points ls */
  } /* end of loop over integration points lr */


  return;
} // DRT::Elements::Fluid2::f2_sys_mat



/*----------------------------------------------------------------------*
 |  evaluate the element integration points (private)        gammi 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid2::f2_integration_points(struct _FLUID_DATA& data)
{
  const DOUBLE Q12 = ONE/TWO;
  const DOUBLE Q13 = ONE/THREE;
  const DOUBLE Q16 = ONE/SIX;
  const DOUBLE Q23 = TWO/THREE;

/*----------------------------------------------------------------------*/
  /* inintialise triX arrays */
  for (int i=0; i<MAXTINTP; i++) /* loop over all integration points    */
  {
      for (int k=0; k<MAXTINTC; k++) /* loop integration cases          */
      {
	  /* set coordinates (r,s) coordinates of integration point     */
	  data.txgr[i][k] = ZERO;
	  data.txgs[i][k] = ZERO;

	  /* innitialise the vector of gaussweights */
	  data.twgt[i][k] = ZERO;
      }
  }
  /* inintialise quadX arrays */
  for (int i=0; i<MAXQINTP; i++) /* loop over integration points       */
  {
      for (int k=0; k<MAXQINTC; k++) /* loop integration cases         */
      {
	  /* set one coordinate of integration points --- the rest is
	   * 'symmetric'                                               */
	  data.qxg [i][k] = ZERO;

	  /* innitialise the vector of gaussweights */
	  data.qwgt[i][k] = ZERO;
      }
  }

/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    R E C T A N G U L A R   ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 |    data.qxg[i][j]                                                    |
 |   data.qwgt[i][j]:  i+1 - actual number of gausspoint                |
 |                     j+1 - total number of gausspoints                |
 *----------------------------------------------------------------------*/
/* coordinates for two gauss points */
      data.qxg[0][1]  =  -0.5773502691896;
      data.qxg[1][1]  =  -data.qxg[0][1] ;
/* coordinates for three gauss points */
      data.qxg[0][2]  =  -0.7745966692415;
      data.qxg[2][2]  =  -data.qxg[0][2] ;
/* coordinates for four gauss points */
      data.qxg[0][3]  =  -0.8611363115941;
      data.qxg[1][3]  =  -0.3399810435849;
      data.qxg[2][3]  =  -data.qxg[1][3] ;
      data.qxg[3][3]  =  -data.qxg[0][3] ;
/* coordinates for five gauss points */
      data.qxg[0][4]  =  -0.9061798459387;
      data.qxg[1][4]  =  -0.5384693101057;
      data.qxg[3][4]  =  -data.qxg[1][4] ;
      data.qxg[4][4]  =  -data.qxg[0][4] ;
/* coordinates for six gauss points */
      data.qxg[0][5]  =  -0.9324695142032;
      data.qxg[1][5]  =  -0.6612093864663;
      data.qxg[2][5]  =  -0.2386191860832;
      data.qxg[3][5]  =  -data.qxg[2][5] ;
      data.qxg[4][5]  =  -data.qxg[1][5] ;
      data.qxg[5][5]  =  -data.qxg[0][5] ;

/* weights for one gauss points */
      data.qwgt[0][0] =  TWO             ;
/* weights for two gauss points */
      data.qwgt[0][1] =  ONE             ;
      data.qwgt[1][1] =  ONE             ;
/* weights for three gauss points */
      data.qwgt[0][2] =  0.5555555555556 ;
      data.qwgt[1][2] =  0.8888888888889 ;
      data.qwgt[2][2] =  data.qwgt[0][2] ;
/* weights for four gauss points */
      data.qwgt[0][3] =  0.3478548451375 ;
      data.qwgt[1][3] =  0.6521451548625 ;
      data.qwgt[2][3] =  data.qwgt[1][3] ;
      data.qwgt[3][3] =  data.qwgt[0][3] ;
/* weights for five gauss points */
      data.qwgt[0][4] =  0.2369268850562 ;
      data.qwgt[1][4] =  0.4786286704994 ;
      data.qwgt[2][4] =  0.5688888888889 ;
      data.qwgt[3][4] =  data.qwgt[1][4] ;
      data.qwgt[4][4] =  data.qwgt[0][4] ;
/* weights for six gauss points */
      data.qwgt[0][5] =  0.1713244923792 ;
      data.qwgt[1][5] =  0.3607615730481 ;
      data.qwgt[2][5] =  0.4679139345727 ;
      data.qwgt[3][5] =  data.qwgt[2][5] ;
      data.qwgt[4][5] =  data.qwgt[1][5] ;
      data.qwgt[5][5] =  data.qwgt[0][5] ;


/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    T R I A N G U L A R     ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 |   data.txgr[i][j]                                                    |
 |  data.twgts[i][j]:  i+1 - actual number of gausspoint                |
 |                     j+1 - number for integration case (from input)   |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
 |                              CASE 0                                  |
 *----------------------------------------------------------------------*/
      data.txgr[0][0]    =  Q13 ;
      data.txgs[0][0]    =  Q13 ;

      data.twgt[0][0]    =  Q12 ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        3 SAMPLING POINTS, DEG.OF PRECISION 2    |
 |                             CASE 1                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][1]    =  Q12  ;
      data.txgr[1][1]    =  Q12  ;
      data.txgr[2][1]    =  ZERO ;
      data.txgs[0][1]    =  ZERO ;
      data.txgs[1][1]    =  Q12  ;
      data.txgs[2][1]    =  Q12  ;

      data.twgt[0][1]    =  Q16  ;
      data.twgt[1][1]    =  Q16  ;
      data.twgt[2][1]    =  Q16  ;
/*----------------------------------------------------------------------*
 |    ALT.GAUSS INTEGRATION    3 SAMPLING POINTS, DEG.OF PRECISION 2    |
 |                             CASE 2                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][2]    =  Q16  ;
      data.txgr[1][2]    =  Q23  ;
      data.txgr[2][2]    =  Q16  ;
      data.txgs[0][2]    =  Q16  ;
      data.txgs[1][2]    =  Q16  ;
      data.txgs[2][2]    =  Q23  ;

      data.twgt[0][2]    =  Q16  ;
      data.twgt[1][2]    =  Q16  ;
      data.twgt[2][2]    =  Q16  ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 3    |
 |                             CASE 3                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][3]    =  0.2                ;
      data.txgr[1][3]    =  0.6                ;
      data.txgr[2][3]    =  0.2                ;
      data.txgr[3][3]    =  Q13                ;
      data.txgs[0][3]    =  0.2                ;
      data.txgs[1][3]    =  0.2                ;
      data.txgs[2][3]    =  0.6                ;
      data.txgs[3][3]    =  Q13                ;

      data.twgt[0][3]    =  0.2604166666667    ;
      data.twgt[1][3]    =  data.twgt[0][2]    ;
      data.twgt[2][3]    =  data.twgt[0][2]    ;
      data.twgt[3][3]    = -0.28125            ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        6 SAMPLING POINTS, DEG.OF PRECISION 4    |
 |                             CASE 4                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][4]    =  0.0915762135098	;
      data.txgr[1][4]    =  0.8168475729805	;
      data.txgr[2][4]    =  0.0915762135098 	;
      data.txgr[3][4]    =  0.4459484909160	;
      data.txgr[4][4]    =  0.4459484909160 	;
      data.txgr[5][4]    =  0.1081030181681	;
      data.txgs[0][4]    =  0.0915762135098 	;
      data.txgs[1][4]    =  0.0915762135098 	;
      data.txgs[2][4]    =  0.8168475729805 	;
      data.txgs[3][4]    =  0.1081030181681 	;
      data.txgs[4][4]    =  0.4459484909160 	;
      data.txgs[5][4]    =  0.4459484909160 	;

      data.twgt[0][4]   =  0.0549758718277	;
      data.twgt[1][4]   =  0.0549758718277	;
      data.twgt[2][4]   =  0.0549758718277	;
      data.twgt[3][4]   =  0.1116907948390	;
      data.twgt[4][4]   =  0.1116907948390	;
      data.twgt[5][4]   =  0.1116907948390	;
/*----------------------------------------------------------------------*
 |    ALT.GAUSS INTEGRATION    6 SAMPLING POINTS, DEG.OF PRECISION 3    |
 |                             CASE 5                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][5]    =  0.1090390090729	;
      data.txgr[1][5]    =  0.2319333685530	;
      data.txgr[2][5]    =  0.6590276223741	;
      data.txgr[3][5]    =  0.0915762135098 	;
      data.txgr[4][5]    =  0.8168475729805 	;
      data.txgr[5][5]    =  0.0915762135098 	;
      data.txgs[0][5]    =  0.8168475729805 	;
      data.txgs[1][5]    =  0.0915762135098 	;
      data.txgs[2][5]    =  0.0915762135098 	;
      data.txgs[3][5]    =  0.8168475729805 	;
      data.txgs[4][5]    =  0.0915762135098 	;
      data.txgs[5][5]    =  0.0915762135098 	;



      data.twgt[0][5]   =  0.0833333333333	;
      data.twgt[1][5]   =  0.0833333333333	;
      data.twgt[2][5]   =  0.0833333333333	;
      data.twgt[3][5]   =  0.0833333333333	;
      data.twgt[4][5]   =  0.0833333333333	;
      data.twgt[5][5]   =  0.0833333333333	;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        7 SAMPLING POINTS, DEG.OF PRECISION 5    |
 |                             CASE 6                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][6]    =  0.1012865073235 ;
      data.txgr[1][6]    =  0.4701420641051 ;
      data.txgr[2][6]    =  0.7974269853531 ;
      data.txgr[3][6]    =  data.txgr[1][4]       ;
      data.txgr[4][6]    =  data.txgr[0][4]       ;
      data.txgr[5][6]    =  0.0597158717898 ;
      data.txgr[6][6]    =  Q13	      ;
      data.txgs[0][6]    =  data.txgr[0][4]       ;
      data.txgs[1][6]    =  data.txgr[5][4]       ;
      data.txgs[2][6]    =  data.txgr[0][4]       ;
      data.txgs[3][6]    =  data.txgr[1][4]       ;
      data.txgs[4][6]    =  data.txgr[2][4]       ;
      data.txgs[5][6]    =  data.txgr[1][4]       ;
      data.txgs[6][6]    =  Q13	      ;

      data.twgt[0][6]    =  0.0629695902724 ;
      data.twgt[1][6]    =  0.0661970763943 ;
      data.twgt[2][6]    =  data.twgt[0][4]      ;
      data.twgt[3][6]    =  data.twgt[1][4]      ;
      data.twgt[4][6]    =  data.twgt[0][4]      ;
      data.twgt[5][6]    =  data.twgt[1][4]      ;
      data.twgt[6][6]    =  0.1125	      ;
/*----------------------------------------------------------------------*
 |    ALT.GAUSS INTEGRATION    7 SAMPLING POINTS, DEG.OF PRECISION 4    |
 |                             CASE 7                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][7]    =  0.2379323664724 ;
      data.txgr[1][7]    =  0.7367124989684 ;
      data.txgr[2][7]    =  data.txgr[1][4] ;
      data.txgr[3][7]    =  data.txgr[0][4] ;
      data.txgr[4][7]    =  0.0253551345591 ;
      data.txgr[5][7]    =  data.txgr[4][4] ;
      data.txgr[6][7]    =  Q13	            ;
      data.txgs[0][7]    =  data.txgr[4][4] ;
      data.txgs[1][7]    =  data.txgr[4][4] ;
      data.txgs[2][7]    =  data.txgr[0][4] ;
      data.txgs[3][7]    =  data.txgr[1][4] ;
      data.txgs[4][7]    =  data.txgr[1][4] ;
      data.txgs[5][7]    =  data.txgr[0][4] ;
      data.txgs[6][7]    =  Q13	            ;

      data.twgt[0][7]    =  0.0520833333333 ;
      data.twgt[1][7]    =  data.twgt[0][4] ;
      data.twgt[2][7]    =  data.twgt[0][4] ;
      data.twgt[3][7]    =  data.twgt[0][4] ;
      data.twgt[4][7]    =  data.twgt[0][4] ;
      data.twgt[5][7]    =  data.twgt[0][4] ;
      data.twgt[6][7]    =  0.1875	    ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        9 SAMPLING POINTS, DEG.OF PRECISION 5    |
 |                             CASE 8                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][8]    =  0.1654099273898     ;
      data.txgr[1][8]    =  0.4375252483834     ;
      data.txgr[2][8]    =  0.7971126518601     ;
      data.txgr[3][8]    =  data.txgr[2][5]     ;
      data.txgr[4][8]    =  data.txgr[1][5]     ;
      data.txgr[5][8]    =  data.txgr[0][5]     ;
      data.txgr[6][8]    =  0.0374774207501     ;
      data.txgr[7][8]    =  0.1249495032332     ;
      data.txgr[8][8]    =  data.txgr[6][5]     ;
      data.txgs[0][8]    =  data.txgr[6][5]     ;
      data.txgs[1][8]    =  data.txgr[7][5]     ;
      data.txgs[2][8]    =  data.txgr[6][5]     ;
      data.txgs[3][8]    =  data.txgr[0][5]     ;
      data.txgs[4][8]    =  data.txgr[1][5]     ;
      data.txgs[5][8]    =  data.txgr[2][5]     ;
      data.txgs[6][8]    =  data.txgr[2][5]     ;
      data.txgs[7][8]    =  data.txgr[1][5]     ;
      data.txgs[8][8]    =  data.txgr[0][5]     ;

      data.twgt[0][8]    =  0.0318457071431     ;
      data.twgt[1][8]    =  0.1029752523804     ;
      data.twgt[2][8]    =  data.twgt[0][5]     ;
      data.twgt[3][8]    =  data.twgt[0][5]     ;
      data.twgt[4][8]    =  data.twgt[1][5]     ;
      data.twgt[5][8]    =  data.twgt[0][5]     ;
      data.twgt[6][8]    =  data.twgt[0][5]     ;
      data.twgt[7][8]    =  data.twgt[1][5]     ;
      data.twgt[8][8]    =  data.twgt[0][5]     ;
 /*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION       12 SAMPLING POINTS, DEG.OF PRECISION 6    |
 |                            CASE 9                                    |
 *----------------------------------------------------------------------*/
      data.txgr[ 0][9]   =  0.0630890144915     ;
      data.txgr[ 1][9]   =  0.3103524510338     ;
      data.txgr[ 2][9]   =  0.6365024991214     ;
      data.txgr[ 3][9]   =  0.8738219710170     ;
      data.txgr[ 4][9]   =  data.txgr[ 2][6]    ;
      data.txgr[ 5][9]   =  data.txgr[ 1][6]    ;
      data.txgr[ 6][9]   =  data.txgr[ 0][6]    ;
      data.txgr[ 7][9]   =  0.0531450498448     ;
      data.txgr[ 8][9]   =  data.txgr[ 7][6]    ;
      data.txgr[ 9][9]   =  0.2492867451709     ;
      data.txgr[10][9]   =  0.5014265096582     ;
      data.txgr[11][9]   =  data.txgr[ 9][6]    ;
      data.txgs[ 0][9]   =  data.txgr[ 0][6]    ;
      data.txgs[ 1][9]   =  data.txgr[ 7][6]    ;
      data.txgs[ 2][9]   =  data.txgr[ 7][6]    ;
      data.txgs[ 3][9]   =  data.txgr[ 0][6]    ;
      data.txgs[ 4][9]   =  data.txgr[ 1][6]    ;
      data.txgs[ 5][9]   =  data.txgr[ 2][6]    ;
      data.txgs[ 6][9]   =  data.txgr[ 3][6]    ;
      data.txgs[ 7][9]   =  data.txgr[ 2][6]    ;
      data.txgs[ 8][9]   =  data.txgr[ 1][6]    ;
      data.txgs[ 9][9]   =  data.txgr[ 9][6]    ;
      data.txgs[10][9]   =  data.txgr[ 9][6]    ;
      data.txgs[11][9]   =  data.txgr[10][6]    ;
      data.twgt[ 0][9]   =  0.0254224531851     ;
      data.twgt[ 1][9]   =  0.0414255378092     ;
      data.twgt[ 2][9]   =  data.twgt[ 1][6]    ;
      data.twgt[ 3][9]   =  data.twgt[ 0][6]    ;
      data.twgt[ 4][9]   =  data.twgt[ 1][6]    ;
      data.twgt[ 5][9]   =  data.twgt[ 1][6]    ;
      data.twgt[ 6][9]   =  data.twgt[ 0][6]    ;
      data.twgt[ 7][9]   =  data.twgt[ 1][6]    ;
      data.twgt[ 8][9]   =  data.twgt[ 1][6]    ;
      data.twgt[ 9][9]   =  0.0583931378632     ;
      data.twgt[10][9]   =  data.twgt[ 9][6]    ;
      data.twgt[11][6]   =  data.twgt[ 9][6]    ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION       13 SAMPLING POINTS, DEG.OF PRECISION 7    |
 |                            CASE 10                                   |
 *----------------------------------------------------------------------*/
      data.txgr[ 0][10]  =  0.0651301029022     ;
      data.txgr[ 1][10]  =  0.3128654960049     ;
      data.txgr[ 2][10]  =  0.6384441885698     ;
      data.txgr[ 3][10]  =  0.8697397941956     ;
      data.txgr[ 4][10]  =  data.txgr[ 2][7]    ;
      data.txgr[ 5][10]  =  data.txgr[ 1][7]    ;
      data.txgr[ 6][10]  =  data.txgr[ 0][7]    ;
      data.txgr[ 7][10]  =  0.0486903154253     ;
      data.txgr[ 8][10]  =  data.txgr[ 7][7]    ;
      data.txgr[ 9][10]  =  0.2603459660790     ;
      data.txgr[10][10]  =  0.4793080678419     ;
      data.txgr[11][10]  =  data.txgr[ 9][7]    ;
      data.txgr[12][10]  =  Q13                 ;
      data.txgs[ 0][10]  =  data.txgr[ 0][7]    ;
      data.txgs[ 1][10]  =  data.txgr[ 7][7]    ;
      data.txgs[ 2][10]  =  data.txgr[ 7][7]    ;
      data.txgs[ 3][10]  =  data.txgr[ 0][7]    ;
      data.txgs[ 4][10]  =  data.txgr[ 1][7]    ;
      data.txgs[ 5][10]  =  data.txgr[ 2][7]    ;
      data.txgs[ 6][10]  =  data.txgr[ 3][7]    ;
      data.txgs[ 7][10]  =  data.txgr[ 2][7]    ;
      data.txgs[ 8][10]  =  data.txgr[ 1][7]    ;
      data.txgs[ 9][10]  =  data.txgr[ 9][7]    ;
      data.txgs[10][10]  =  data.txgr[ 9][7]    ;
      data.txgs[11][10]  =  data.txgr[10][7]    ;
      data.txgs[12][10]  =  Q13                 ;
      data.twgt[ 0][10]  =  0.0266736178044     ;
      data.twgt[ 1][10]  =  0.0385568804451     ;
      data.twgt[ 2][10]  =  data.twgt[ 1][7]    ;
      data.twgt[ 3][10]  =  data.twgt[ 0][7]    ;
      data.twgt[ 4][10]  =  data.twgt[ 1][7]    ;
      data.twgt[ 5][10]  =  data.twgt[ 1][7]    ;
      data.twgt[ 6][10]  =  data.twgt[ 0][7]    ;
      data.twgt[ 7][10]  =  data.twgt[ 1][7]    ;
      data.twgt[ 8][10]  =  data.twgt[ 1][7]    ;
      data.twgt[ 9][10]  =  0.0878076287166     ;
      data.twgt[10][10]  =  data.twgt[ 9][7]    ;
      data.twgt[11][10]  =  data.twgt[ 9][7]    ;
      data.twgt[12][10]  = -0.0747850222338     ;

  return;
} //end of DRT::Elements::Fluid2::f2_integration_points




/*----------------------------------------------------------------------*
 | calculate Jacobian matrix and it's determinant (private) gammi  04/07|
 | Well, I think we actually compute its transpose....
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
void DRT::Elements::Fluid2::f2_jaco(const Epetra_SerialDenseMatrix& xyze,
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
        dum=ZERO;
        for (int l=0; l<iel; l++)
        {
           dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
     } /* end of loop j */
  } /* end of loop i */

  /*------------------------------------------ determinant of jacobian---*/
  *det = xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1);

  if(*det<ZERO)
  {
     printf("\n");
     printf("GLOBAL ELEMENT NO.%i\n",Id());
     printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n",*det);
     dserror("Stopped not regulary!\n");
  }

} //end of DRT::Elements::Fluid2::f2_jaco





/*----------------------------------------------------------------------*
  shape functions and natural derivatives for quadrilaterals (private)
                                                             gammi 04/07

In this routine the shape functions (always) and their natural first
(icode==2) and second derivatives (icode==3) with respect to r/s are
evaluated for  R E C T A N G L E S or T R I A N G L E S

Numbering of the nodes:


                    ^ s
                    |
              1     |4    0
              o-----o-----o
              |           |
              |           |7
            5 o     o     o -----> r
              |     8     |
              |           |
              o-----o-----o
              2     6     3





                    ^ s
                    |
                    |
                   2|
                    o
                    |\
                    | \
                   5o  o4
                    |   \
                    |    \
                    o--o--o -----> r
	           0   3   1




 \param   funct    vector<double>&             (o)    shape functions
 \param   deriv    Epetra_SerialDenseMatrix&   (o)    1st natural deriv.
                                                      of shape funct.
 \param   deriv2   Epetra_SerialDenseMatrix&   (o)    2nd natural deriv.
                                                      of shape funct.
 \param   r        DOUBLE&                     (i)    coordinate
 \param   s        DOUBLE&                     (i)    coordinate
 \param   iel      const int&                  (i)    number of nodes
 \param   icode    int         	               (i)    flag on (higher
                                                      order) derivatives

 *----------------------------------------------------------------------*/

void DRT::Elements::Fluid2::f2_shape_function(
               vector<double>&    		funct,
               Epetra_SerialDenseMatrix& 	deriv,
               Epetra_SerialDenseMatrix& 	deriv2,
               const double&      r,
               const double&      s,
               const int&         iel,
               int         	  icode
            )
{
  const DOUBLE Q12 = ONE/TWO;
  const DOUBLE Q14 = ONE/FOUR;

/*------------------------------- selection of polynomial interpolation */
  switch (iel)
  {
    case 4: /* LINEAR shape functions for quad4 and their natural
	     *                                          derivatives ----*/
    {
      /*--------------------------------------------- form basic values */
	double rp=ONE+r;
	double rm=ONE-r;
	double sp=ONE+s;
	double sm=ONE-s;

	funct[0]=Q14*rp*sp;
	funct[1]=Q14*rm*sp;
	funct[2]=Q14*rm*sm;
	funct[3]=Q14*rp*sm;

	if(icode>1) /* --> first derivative evaluation */
	{
	    deriv(0,0)= Q14*sp;
	    deriv(1,0)= Q14*rp;

	    deriv(0,1)=-Q14*sp;
	    deriv(1,1)= Q14*rm;

	    deriv(0,2)=-Q14*sm;
	    deriv(1,2)=-Q14*rm;

	    deriv(0,3)= Q14*sm;
	    deriv(1,3)=-Q14*rp;
	} /* endif (icode>1) */

	if(icode==3) /* --> second derivative evaluation */
	{
	    deriv2(0,0)= ZERO;
	    deriv2(1,0)= ZERO;
	    deriv2(2,0)= Q14;

	    deriv2(0,1)= ZERO;
	    deriv2(1,1)= ZERO;
	    deriv2(2,1)=-Q14;

	    deriv2(0,2)= ZERO;
	    deriv2(1,2)= ZERO;
	    deriv2(2,2)= Q14;

	    deriv2(0,3)= ZERO;
	    deriv2(1,3)=ZERO;
	    deriv2(2,3)=-Q14;
	} /* endif (icode==3) */
	break;
    }
    case 8: /* QUADRATIC shape functions for quadrilaterals without
  	       central node and their natural derivatives (serendipity) */
    {
	double rp=ONE+r;
	double rm=ONE-r;
	double sp=ONE+s;
	double sm=ONE-s;
	double r2=ONE-r*r;
	double s2=ONE-s*s;

	funct[4]=Q12*r2*sp;
	funct[5]=Q12*rm*s2;
	funct[6]=Q12*r2*sm;
	funct[7]=Q12*rp*s2;
	funct[0]=Q14*rp*sp-Q12*(funct[4]+funct[7]);
	funct[1]=Q14*rm*sp-Q12*(funct[4]+funct[5]);
	funct[2]=Q14*rm*sm-Q12*(funct[5]+funct[6]);
	funct[3]=Q14*rp*sm-Q12*(funct[6]+funct[7]);

	if(icode>1) /* --> first derivative evaluation */
	{
	    deriv(0,0)= Q14*sp;
	    deriv(1,0)= Q14*rp;

	    deriv(0,1)=-Q14*sp;
	    deriv(1,1)= Q14*rm;

	    deriv(0,2)=-Q14*sm;
	    deriv(1,2)=-Q14*rm;

	    deriv(0,3)= Q14*sm;
	    deriv(1,3)=-Q14*rp;

	    deriv(0,4)=-ONE*r*sp;
	    deriv(1,4)= Q12*r2;

	    deriv(0,5)=-Q12*  s2;
	    deriv(1,5)=-ONE*rm*s;

	    deriv(0,6)=-ONE*r*sm;
	    deriv(1,6)=-Q12*r2;

	    deriv(0,7)= Q12*  s2;
	    deriv(1,7)=-ONE*rp*s;

	    deriv(0,0)-= Q12*(deriv(0,4)+deriv(0,7));
	    deriv(1,0)-= Q12*(deriv(1,4)+deriv(1,7));

	    for(int i=1;i<4;i++)
	    {
		int ii=i+3;
		deriv(0,i) -= Q12*(deriv(0,ii)+deriv(0,ii+1));
		deriv(1,i) -= Q12*(deriv(1,ii)+deriv(1,ii+1));
	    } /* end loop over i */
	} /* endif (icode>1) */

	if(icode==3) /* --> second derivative evaluation */
	{
	    deriv2(0,0)= ZERO;
	    deriv2(1,0)= ZERO;
	    deriv2(2,0)= Q14;

	    deriv2(0,1)= ZERO;
	    deriv2(1,1)= ZERO;
	    deriv2(2,1)=-Q14;

	    deriv2(0,2)= ZERO;
	    deriv2(1,2)= ZERO;
	    deriv2(2,2)= Q14;

	    deriv2(0,3)= ZERO;
	    deriv2(1,3)= ZERO;
	    deriv2(2,3)=-Q14;

	    deriv2(0,4)=-(ONE+s);
	    deriv2(1,4)= ZERO;
	    deriv2(2,4)=-r;

	    deriv2(0,5)= ZERO;
	    deriv2(1,5)=-(ONE-r);
	    deriv2(2,5)= s;

	    deriv2(0,6)=-(ONE-s);
	    deriv2(1,6)= ZERO;
	    deriv2(2,6)= r;

	    deriv2(0,7)= ZERO;
	    deriv2(1,7)=-(ONE+r);
	    deriv2(2,7)=-s;

	    deriv2(0,0) -= Q12*(deriv2(0,4)+deriv2(0,7));
	    deriv2(1,0) -= Q12*(deriv2(1,4)+deriv2(1,7));
	    deriv2(2,0) -= Q12*(deriv2(2,4)+deriv2(2,7));

	    for(int i=1;i<4;i++)
	    {
		int ii=i+3;
		deriv2(0,i) -= Q12*(deriv2(0,ii)+deriv2(0,ii+1));
		deriv2(1,i) -= Q12*(deriv2(1,ii)+deriv2(1,ii+1));
		deriv2(2,i) -= Q12*(deriv2(2,ii)+deriv2(2,ii+1));
	    } /* end loop over i */
	} /* endif (icode==3) */
	break;
    }
    case 9: /* full QUADRATIC shape functions for quadrilaterals with
	                     central node and their natural derivatives */
    {
/*--------------------------------------------------- form basic values */
	double rp=ONE+r;
	double rm=ONE-r;
	double sp=ONE+s;
	double sm=ONE-s;
	double r2=ONE-r*r;
	double s2=ONE-s*s;
	double rh=Q12*r;
	double sh=Q12*s;
	double rs=rh*sh;
	double rhp=r+Q12;
	double rhm=r-Q12;
	double shp=s+Q12;
	double shm=s-Q12;

	funct[0]= rs*rp*sp;
	funct[1]=-rs*rm*sp;
	funct[2]= rs*rm*sm;
	funct[3]=-rs*rp*sm;
	funct[4]= sh*sp*r2;
	funct[5]=-rh*rm*s2;
	funct[6]=-sh*sm*r2;
	funct[7]= rh*rp*s2;
	funct[8]= r2*s2;

	if(icode>1) /* --> first derivative evaluation */
	{
	    deriv(0,0)= rhp*sh*sp;
	    deriv(1,0)= shp*rh*rp;

	    deriv(0,1)= rhm*sh*sp;
	    deriv(1,1)=-shp*rh*rm;

	    deriv(0,2)=-rhm*sh*sm;
	    deriv(1,2)=-shm*rh*rm;

	    deriv(0,3)=-rhp*sh*sm;
	    deriv(1,3)= shm*rh*rp;

	    deriv(0,4)=-TWO*r*sh*sp;
	    deriv(1,4)= shp*r2;

	    deriv(0,5)= rhm*s2;
	    deriv(1,5)= TWO*s*rh*rm;

	    deriv(0,6)= TWO*r*sh*sm;
	    deriv(1,6)= shm*r2;

	    deriv(0,7)= rhp*s2;
	    deriv(1,7)=-TWO*s*rh*rp;

	    deriv(0,8)=-TWO*r*s2;
	    deriv(1,8)=-TWO*s*r2;
	} /* endif (icode>1) */

	if(icode==3) /* --> second derivative evaluation */
	{
	    deriv2(0,0)= sh*sp;
	    deriv2(1,0)= rh*rp;
	    deriv2(2,0)= shp*rhp;

	    deriv2(0,1)= sh*sp;
	    deriv2(1,1)=-rh*rm;
	    deriv2(2,1)= shp*rhm;

	    deriv2(0,2)=-sh*sm;
	    deriv2(1,2)=-rh*rm;
	    deriv2(2,2)= shm*rhm;

	    deriv2(0,3)=-sh*sm;
	    deriv2(1,3)= rh*rp;
	    deriv2(2,3)= shm*rhp;

	    deriv2(0,4)=-TWO*sh*sp;
	    deriv2(1,4)= r2;
	    deriv2(2,4)=-TWO*r*shp;

	    deriv2(0,5)= s2;
	    deriv2(1,5)= TWO*rh*rm;
	    deriv2(2,5)=-TWO*s*rhm;

	    deriv2(0,6)= TWO*sh*sm;
	    deriv2(1,6)= r2;
	    deriv2(2,6)=-TWO*r*shm;

	    deriv2(0,7)= s2;
	    deriv2(1,7)=-TWO*rh*rp;
	    deriv2(2,7)=-TWO*s*rhp;

	    deriv2(0,8)=-TWO*s2;
	    deriv2(1,8)=-TWO*r2;
	    deriv2(2,8)= TWO*s*TWO*r;
	} /* endif (icode==3) */
	break;
    }
    case 3: /* LINEAR shape functions for triangles and their natural
	     *                                         derivatives -----*/
    {
        /*------------------------------------------- form basic values */
	funct[0]=ONE-r-s;
	funct[1]=r;
	funct[2]=s;

	if(icode>1) /* --> first derivative evaluation */
	{
	    deriv(0,0)=-ONE;
	    deriv(1,0)=-ONE;
	    deriv(0,1)= ONE;
	    deriv(1,1)=ZERO;
	    deriv(0,2)=ZERO;
	    deriv(1,2)= ONE;
	} /* endif (icode>1) */
	break;
    }
    case 6: /* QUADRATIC shape functions for triangles and their natural
	     *                                             derivatives -*/
    {
        /*------------------------------------------- form basic values */
	double rr=r*r;
	double ss=s*s;
	double rs=r*s;

	funct[0]=(ONE-TWO*r-TWO*s)*(ONE-r-s);
	funct[1]=TWO*rr-r;
	funct[2]=TWO*ss-s;
	funct[3]=FOUR*(r-rr-rs);
	funct[4]=FOUR*rs;
	funct[5]=FOUR*(s-rs-ss);

	if(icode>1) /* --> first derivative evaluation */
	{
	    deriv(0,0)=-THREE+FOUR*(r+s);
	    deriv(1,0)= deriv(0,0);

	    deriv(0,1)= FOUR*r-ONE;
	    deriv(1,1)= ZERO;

	    deriv(0,2)= ZERO;
	    deriv(1,2)= FOUR*s-ONE;

	    deriv(0,3)= FOUR*(ONE-TWO*r-s);
	    deriv(1,3)=-FOUR*r;

	    deriv(0,4)= FOUR*s;
	    deriv(1,4)= FOUR*r;

	    deriv(0,5)=-FOUR*s;
	    deriv(1,5)= FOUR*(ONE-r-TWO*s);
	} /* endif (icode>1) */

	if(icode==3) /* --> second derivative evaluation */
	{
	    deriv2(0,0)= FOUR;
	    deriv2(1,0)= FOUR;
	    deriv2(2,0)= FOUR;

	    deriv2(0,1)= FOUR;
	    deriv2(1,1)= ZERO;
	    deriv2(2,1)= ZERO;

	    deriv2(0,2)= ZERO;
	    deriv2(1,2)= FOUR;
	    deriv2(2,2)= ZERO;

	    deriv2(0,3)=-EIGHT;
	    deriv2(1,3)= ZERO;
	    deriv2(2,3)=-FOUR;

	    deriv2(0,4)= ZERO;
	    deriv2(1,4)= ZERO;
	    deriv2(2,4)= FOUR;

	    deriv2(0,5)= ZERO;
	    deriv2(1,5)=-EIGHT;
	    deriv2(2,5)=-FOUR;
	} /* endif (icode==3) */
	break;
    }
    /*------------------------------------------------------------------*/
    default:
	dserror("distyp unknown\n");
  } /* end switch(iel) */

return;
} // end of DRT:Elements:Fluid2:f2_shape_function



/*----------------------------------------------------------------------*
 |  calculate global derivatives w.r.t. x,y at point r,s (private)
 |                                                        gammi 04/07
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
void DRT::Elements::Fluid2::f2_gder(Epetra_SerialDenseMatrix& derxy,
				    const Epetra_SerialDenseMatrix& deriv,
                                    Epetra_SerialDenseMatrix& xjm,
				    double& det,
                                    const int iel
				    )
{
  Epetra_SerialDenseMatrix 	xji(2,2);   // inverse of jacobian matrix


  /*----------calculate global derivatives w.r.t. x,y at point r,s ---*/

  /*------------------------------------------------------- initialistion */
  for(int k=0;k<iel;k++)
  {
    derxy(0,k)=ZERO;
    derxy(1,k)=ZERO;
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
} // end of DRT:Elements:Fluid2:f2_gder

/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y at point r,s (private)
 |                                                          gammi 04/07
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

void DRT::Elements::Fluid2::f2_gder2(const Epetra_SerialDenseMatrix& xyze,
				     const Epetra_SerialDenseMatrix& xjm,
				     const Epetra_SerialDenseMatrix& derxy,
				     Epetra_SerialDenseMatrix& derxy2,
				     const Epetra_SerialDenseMatrix& deriv2,
				     const int iel
				    )
{
//--------------------------------------------initialize and zero out everything
    Epetra_SerialDenseMatrix bm(3,3);
    Epetra_SerialDenseMatrix xder2(3,2);
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
} // end of DRT:Elements:Fluid2:f2_gder2



/*----------------------------------------------------------------------*
 |  evaluate fluid coefficient matrix (private)              chfoe 04/04|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilised fluid2 element are calculated. The procedure is
based on the Rothe method of first integrating in time. Hence the
resulting terms include coefficients containing time integration variables
such as theta or delta t which are represented by 'timefac'.

The routine was completed to contain ALE-terms also.         chfoe 11/04

The stabilisation is based on the residuum:

R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u)
    + timefac grad p - rhsint

R_C = div u

The corresponding weighting operators are
L_M = v + timefac u_old * grad v + timefac v * grad u_old
    - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q

L_C = div v

where alpha = -1
      beta  = -1
are sign regulating factors and rhsint differs for different time
These factores are worked in now and cannot be changed any more.

integration schemes:

One-step-Theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old

BDF2:

generalised alpha:


The stabilisation by means of the momentum residuum R_M is of the unusual
type:
   Galerkin parts MINUS sum over elements (stabilising parts)
The stabilisation by means of the continuity equation R_C is done in the
usual way:
   Galerkin parts PLUS sum over elements (stabilising parts)

The calculation proceeds as follows.
1) obtain single (linearised) operators of R_M, R_C, L_M and L_C
2) build Galerkin terms from them
3) build stabilising terms from them
4) build Galerkin and stabilising terms of RHS

NOTE: u_old represents the last iteration value. (The most recent one
      we've got!)

NOTE: Galerkin and stabilisation matrices are calculated within one
      routine.

NOTE: In order to increase the performance plenty of terms are concentrated
      and worked into each other. A lengthy version of the file is available
      from the author.


Notational remarks:

                   /              \
                  | u_x,x   u_x,y |
vderxy = grad u = |               |
                  | u_y,x   u_y,y |
                  \               /

           /                         \
          | u_x,xx   u_x,yy   u_x,xy |
vderxy2 = |                          |
          | u_y,xx   u_y,yy   u_y,xy |
          \                          /

for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *histvec    DOUBLE        (i)   rhs at INT point
\param  *gridvint   DOUBLE        (i)   gridvel at INT point
\param **vderxy     DOUBLE        (i)   global vel derivatives
\param  *vderxy2    DOUBLE        (i)   2nd global vel derivatives
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param  *edeadng    DOUBLE        (i)   dead load at time n+1
\param   fac        DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel        INT           (i)   number of nodes of act. ele
\param  *hasext     INT           (i)   flag, if element has volume load
\param   isale      INT           (i)   flag, if ALE or EULER
\return void
------------------------------------------------------------------------*/

void DRT::Elements::Fluid2::f2_calmat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    vector<double>&           histvec,
    vector<double>&           gridvint,
    double&   	              press,
    Epetra_SerialDenseMatrix& vderxy,
    Epetra_SerialDenseMatrix& vderxy2,
    vector<double>&           gradp,
    vector<double>&           funct,
    vector<double>&           tau,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    double&                   fac,
    const double&             visc,
    const int&                iel,
    ParameterList& 	      params
    )
{
/*========================= further variables =========================*/

  Epetra_SerialDenseMatrix  viscs2(2,2*iel);  /* viscous term incluiding 2nd derivatives         */
  vector<double>            conv_c(iel); 	    /* linearisation of convect, convective part       */
  vector<double>            conv_g(iel); 	    /* linearisation of convect, grid part             */
  Epetra_SerialDenseMatrix  conv_r(2,2*iel);  /* linearisation of convect, reactive part         */
  vector<double>            div(2*iel);  	    /* divergence of u or v                            */
  Epetra_SerialDenseMatrix  ugradv(iel,2*iel);/* linearisation of u * grad v                     */
  vector<double>            conv_old(2);      /* convective term evalaluated with old velocities */
  //vector<double>            conv_g_old(2);
  vector<double>            visc_old(2); 	    /* viscous term evaluated with old velocities      */
  vector<double>            rhsint(2);   	    /* total right hand side terms at int.-point       */

/*========================== initialisation ============================*/
// One-step-Theta: timefac = theta*dt
// BDF2:           timefac = 2/3 * dt
  double timefac = params.get<double>("time constant for integration",-1.0);
  if (timefac < 0.0) dserror("No time constant for integration supplied");

// stabilisation parameter
  double tau_M  = tau[0]*fac;
  double tau_Mp = tau[0]*fac;
  double tau_C  = tau[2]*fac;

// integration factors and coefficients of single terms
  double time2nue   = timefac * 2.0 * visc;
  double timetauM   = timefac * tau_M;
  double timetauMp  = timefac * tau_Mp;

  double ttimetauM  = timefac * timetauM;
  double ttimetauMp = timefac * timetauMp;
  double timefacfac = timefac * fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
  rhsint[0] = histvec[0];
  rhsint[1] = histvec[1];

/*----------------- get numerical representation of single operators ---*/
/* Convective term  u_old * grad u_old: */
  conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1];
  conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1];

/* Viscous term  div epsilon(u_old) */
  visc_old[0] = 0.5 * (2.0*vderxy2(0,0) + vderxy2(0,1) + vderxy2(1,2));
  visc_old[1] = 0.5 * (2.0*vderxy2(1,1) + vderxy2(1,0) + vderxy2(0,2));

  for (int i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    conv_c[i] = derxy(0,i) * velint[0] + derxy(1,i) * velint[1] ;

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    if (is_ale_)
    {
      conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1];
    }
    else
    {
      conv_g[i] = 0;
    }

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                          \
       |  u_old_x,x   u_old_x,y   |
       |                          | * N   with  N .. form function matrix
       |  u_old_y,x   u_old_y,y   |
       \                         /                                       */
    conv_r(0,2*i  ) = vderxy(0,0)*funct[i];
    conv_r(0,2*i+1) = vderxy(0,1)*funct[i];
    conv_r(1,2*i  ) = vderxy(1,0)*funct[i];
    conv_r(1,2*i+1) = vderxy(1,1)*funct[i];

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                              \
         1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
         - - |                              |         N_y .. y-line of N
         2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
         \                             /                                 */
    viscs2(0,2*i  ) = - 0.5 * ( 2.0 * derxy2(0,i) + derxy2(1,i) );
    viscs2(0,2*i+1) = - 0.5 * ( derxy2(2,i) );
    viscs2(1,2*i  ) = - 0.5 * ( derxy2(2,i) );
    viscs2(1,2*i+1) = - 0.5 * ( derxy2(0,i) + 2.0 * derxy2(1,i) );

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--- divergence u term ---------------------------------------------*/
    div[2*i]   = derxy(0,i);
    div[2*i+1] = derxy(1,i);

    /*--- ugradv-Term ---------------------------------------------------*/
    /*
      /                                                          \
      |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
      |                                                          |
      |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
      |                                                          |
      |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
      |                                           .              |
      |  . . .                                        .          |
      |                                                  Ni*Ni,y |
      \                                                          /       */
    /* remark: vgradu = ugradv^T */
    for (int j=0; j<iel; j++)
    {
      ugradv(i,2*j  ) = derxy(0,i) * funct[j];
      ugradv(i,2*j+1) = derxy(1,i) * funct[j];
    }

  } // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define vderxy_(i,j)   vderxy(i,j)
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r(i,2*(k)+j)
#define conv_old_(j)   conv_old[j]
//#define conv_g_old_(j) conv_g_old[j]
#define derxy_(i,j)    derxy(i,j)
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2(i,2*(k)+j)
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

  /* We keep two versions: with and without ale. The laster one is a
   * little faster. (more than 10%) */

  if (is_ale_)
  {
    int vi;
    int ui;
#include "fluid2_stiff_ale.cpp"
#include "fluid2_rhs_incr_ale.cpp"
  }
  else
  {
#include "fluid2_stiff.cpp"
#include "fluid2_rhs_incr.cpp"
  }

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef conv_old_
//#undef conv_g_old_
#undef derxy_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef funct_
#undef vderxy_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

  return;
} // end of DRT:Elements:Fluid2:f2_calmat


void DRT::Elements::Fluid2::f2_calmat_stationary(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    vector<double>&           histvec,
    vector<double>&           gridvint,
    double&   	              press,
    Epetra_SerialDenseMatrix& vderxy,
    Epetra_SerialDenseMatrix& vderxy2,
    vector<double>&           gradp,
    vector<double>&           funct,
    vector<double>&           tau,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    double&                   fac,
    const double&             visc,
    const int&                iel,
    ParameterList& 	      params
    )
{
/*========================= further variables =========================*/

Epetra_SerialDenseMatrix  viscs2(2,2*iel);  /* viscous term incluiding 2nd derivatives         */
vector<double>            conv_c(iel); 	    /* linearisation of convect, convective part       */
vector<double>            conv_g(iel); 	    /* linearisation of convect, grid part             */
Epetra_SerialDenseMatrix  conv_r(2,2*iel);  /* linearisation of convect, reactive part         */
vector<double>            div(2*iel);  	    /* divergence of u or v                            */
Epetra_SerialDenseMatrix  ugradv(iel,2*iel);/* linearisation of u * grad v                     */
vector<double>            conv_old(2);      /* convective term evalaluated with old velocities */
vector<double>            conv_g_old(2);
vector<double>            visc_old(2); 	    /* viscous term evaluated with old velocities      */
vector<double>            rhsint(2);   	    /* total right hand side terms at int.-point       */

/*========================== initialisation ============================*/

// stabilisation parameter
double tau_M  = tau[0]*fac;
double tau_Mp = tau[0]*fac;
double tau_C  = tau[2]*fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
rhsint[0] = histvec[0];
rhsint[1] = histvec[1];

/*----------------- get numerical representation of single operators ---*/
/* Convective term  u_old * grad u_old: */
conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1];
conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1];

/* Viscous term  div epsilon(u_old) */
visc_old[0] = 0.5 * (2.0*vderxy2(0,0) + vderxy2(0,1) + vderxy2(1,2));
visc_old[1] = 0.5 * (2.0*vderxy2(1,1) + vderxy2(1,0) + vderxy2(0,2));

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* Reactive term  u:  funct */
   /* linearise convective term */

   /*--- convective part u_old * grad (funct) --------------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    conv_c[i] = derxy(0,i) * velint[0] + derxy(1,i) * velint[1] ;

   /*--- convective grid part u_G * grad (funct) -----------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   if(is_ale_)
   {
       conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1];
   }
   else
   {
     conv_g[i] = 0;
   }

   /*--- reactive part funct * grad (u_old) ----------------------------*/
   /* /                          \
      |  u_old_x,x   u_old_x,y   |
      |                          | * N   with  N .. form function matrix
      |  u_old_y,x   u_old_y,y   |
      \                         /                                       */
   conv_r(0,2*i  ) = vderxy(0,0)*funct[i];
   conv_r(0,2*i+1) = vderxy(0,1)*funct[i];
   conv_r(1,2*i  ) = vderxy(1,0)*funct[i];
   conv_r(1,2*i+1) = vderxy(1,1)*funct[i];

   /*--- viscous term  - grad * epsilon(u): ----------------------------*/
   /*   /                              \
      1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
    - - |                              |         N_y .. y-line of N
      2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
        \                             /                                 */
   viscs2(0,2*i  ) = - 0.5 * ( 2.0 * derxy2(0,i) + derxy2(1,i) );
   viscs2(0,2*i+1) = - 0.5 * ( derxy2(2,i) );
   viscs2(1,2*i  ) = - 0.5 * ( derxy2(2,i) );
   viscs2(1,2*i+1) = - 0.5 * ( derxy2(0,i) + 2.0 * derxy2(1,i) );

   /* pressure gradient term derxy, funct without or with integration   *
    * by parts, respectively                                            */

   /*--- divergence u term ---------------------------------------------*/
   div[2*i]   = derxy(0,i);
   div[2*i+1] = derxy(1,i);

   /*--- ugradv-Term ---------------------------------------------------*/
   /*
     /                                                          \
     |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
     |                                                          |
     |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
     |                                                          |
     |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
     |                                           .              |
     |  . . .                                        .          |
     |                                                  Ni*Ni,y |
     \                                                          /       */
   /* remark: vgradu = ugradv^T */
   for (int j=0; j<iel; j++)
   {
       ugradv(i,2*j  ) = derxy(0,i) * funct[j];
       ugradv(i,2*j+1) = derxy(1,i) * funct[j];
   }

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define vderxy_(i,j)   vderxy(i,j)
#define conv_c_(j)     conv_c[j]
//#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r(i,2*(k)+j)
#define conv_old_(j)   conv_old[j]
//#define conv_g_old_(j) conv_g_old[j]
#define derxy_(i,j)    derxy(i,j)
//#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2(i,2*(k)+j)
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

  {
   #include "fluid2_stiff_stationary.cpp"
   #include "fluid2_rhs_incr_stationary.cpp"
  }

#undef estif_
#undef eforce_
#undef conv_c_
//#undef conv_g_
#undef conv_r_
#undef conv_old_
//#undef conv_g_old_
#undef derxy_
//#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef funct_
#undef vderxy_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

return;
} // end of DRT:Elements:Fluid2:f2_calmat_stationary


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid2Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
