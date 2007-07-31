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
#include "../drt_lib/drt_timecurve.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/newtonianfluid.H"

extern "C"
{
#include "../headers/standardtypes.h"
}



using namespace DRT::Utils;

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
  RefCountPtr<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

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

      // get control parameter
      const bool is_stationary = params.get<bool>("using stationary formulation",false);
      const double time = params.get<double>("total time",-1.0);

      // One-step-Theta: timefac = theta*dt
      // BDF2:           timefac = 2/3 * dt
      double timefac = 0;
      if (not is_stationary)
      {
        timefac = params.get<double>("time constant for integration",-1.0);
        if (timefac < 0.0) dserror("No time constant for integration supplied");
      }

      // calculate element coefficient matrix and rhs
      f2_sys_mat(lm,myvelnp,myprenp,myvhist,mydispnp,mygridv,&elemat1,&elevec1,actmat,
                 time,timefac,is_stationary);

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
 |  do nothing (public)                                      g.bau 07/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid2 elements, the          |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
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
                                       double                    time,
                                       double                    timefac,
                                       bool                      is_stationary
  )
{

  /*---------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();
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

  // dead load in element nodes
  const Epetra_SerialDenseMatrix bodyforce = f2_getbodyforce(time);

  /*---------------------------------------------- get viscosity ---*/
  // check here, if we really have a fluid !!
  if(material->mattyp != m_fluid) dserror("Material law is not of type m_fluid.");
  const double  visc = material->m.fluid->viscosity;

  /*--------------------------------------------- stab-parameter ---*/
  // USFEM stabilization is default. No switch here at the moment.

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector      funct(iel);
  Epetra_SerialDenseMatrix 	deriv(2,iel);
  Epetra_SerialDenseMatrix 	deriv2(3,iel);
  static Epetra_SerialDenseMatrix 	xjm(2,2);
  static Epetra_SerialDenseMatrix 	vderxy(2,2);
  static vector<double> 	pderxy(2);
  static Epetra_SerialDenseMatrix 	vderxy2(2,3);
  Epetra_SerialDenseMatrix 	derxy(2,iel);
  Epetra_SerialDenseMatrix 	derxy2(3,iel);
  vector<double>                edeadng(3);
  vector<double> 		ephin(iel);
  vector<double> 		ephing(iel);
  vector<double> 		iedgnod(iel);
  //Epetra_SerialDenseMatrix 	wa1(100,100);  // working matrix used as dummy
  static vector<double>    	histvec(2); /* history data at integration point      */
  double         		hk;         /* element length for calculation of tau  */
  double         		vel_norm, pe, re, xi1, xi2, xi;
  vector<double>         	velino(2); /* normed velocity at element centre */
  double         		mk=0.0;
  static vector<double>         velint(2);
  static vector<double>         tau(3); // stab parameters

  /*------------------------------------------------------- initialise ---*/
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
  for (int i=0;i<2;i++)
  {
    velint[i]=0.0;
    for (int j=0;j<iel;j++)
    {
      velint[i] += funct[j]*evelnp[i+(2*j)];
    }
  } //end loop over i

  // get Jacobian matrix and determinant
  double  det;
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

  double    press;
  static vector<double>    gridvelint(2); /* grid velocity                       */
  static vector<double>    gradp(2);      /* pressure gradient at integration point         */

  // flag for higher order elements
  const bool higher_order_ele = is_higher_order_element(distype);

  // gaussian points
  const GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule);


  // integration loops
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

      f2_jaco(xyze,deriv,xjm,&det,iel);
      const double fac = intpoints.qwgt[iquad]*det;

      /*---------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

      /*--------------------------------- compute second global derivative */
      if (higher_order_ele)
      {
        f2_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

        /*------calculate 2nd velocity derivatives at integration point */
        // former f2_vder2(vderxy2,derxy2,evelnp,iel);
        for (int i=0;i<3;i++)
        {
          vderxy2(0,i)=0.0;
          vderxy2(1,i)=0.0;
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
        velint[i]=0.0;
        for (int j=0;j<iel;j++)
        {
          velint[i] += funct[j]*evelnp[i+(2*j)];
        }
      } //end loop over i

      /*---------------- get history data (n,i) at integration point ---*/
      //expression for f2_veci(histvec,funct,evhist,iel);
      for (int i=0;i<2;i++)
      {
        histvec[i]=0.0;
        for (int j=0;j<iel;j++)
        {
          histvec[i] += funct[j]*evhist[i+(2*j)];
        } /* end of loop over j */
      } /* end of loop over i */

      /*----------- get velocity (np,i) derivatives at integration point */
      // expression for f2_vder(vderxy,derxy,evelnp,iel);
      for (int i=0;i<2;i++)
      {
        vderxy(0,i)=0.0;
        vderxy(1,i)=0.0;
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
            gridvelint[i] += funct(j)*egridv[i+(3*j)];
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

     // get bodyforce in gausspoint
     for (int isd=0;isd<3;isd++)
     {
       edeadng[isd] = 0.0;
       for (int inode=0;inode<iel;inode++)
       {
         edeadng[isd]+= bodyforce(isd,inode)*funct[inode];
       }
     }

      // perform integration for entire matrix and rhs
      if(is_stationary==false)
        f2_calmat(*sys_mat,*residual,
                  velint,histvec,gridvelint,press,
                  vderxy,vderxy2,gradp,
                  funct,tau,
                  derxy,derxy2,edeadng,
                  fac,visc,iel,
                  timefac);
      else
        f2_calmat_stationary(*sys_mat,*residual,
                             velint,histvec,gridvelint,press,
                             vderxy,vderxy2,gradp,
                             funct,tau,
                             derxy,derxy2,edeadng,
                             fac,visc,iel);


  } // end of loop over integration points


  return;
} // DRT::Elements::Fluid2::f2_sys_mat


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

} //end of DRT::Elements::Fluid2::f2_jaco


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) g.bau 07/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a surface Neumann condition   |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::Elements::Fluid2::f2_getbodyforce(
        const double          time
)
{
  const int iel = NumNode();
  const int nsd = 2; // number of space dimensions
  Epetra_SerialDenseMatrix edeadng(nsd,iel);

  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique surface Neumann condition
  int nodecount = 0;
  for(int inode=0;inode<iel;inode++)
  {
    Nodes()[inode]->GetCondition("SurfaceNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one SurfaceNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  if (nodecount == iel)
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
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
         dserror("Negative time value in body force calculation: time = %f",time);
	// curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for(int jnode=0;jnode<iel;jnode++)
    {
      Nodes()[jnode]->GetCondition("SurfaceNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<nsd;isd++)
      {
        edeadng(isd,jnode)=(*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    for(int inode=0;inode<iel;inode++)
    {
      for(int isd=0;isd<nsd;isd++)
      {
        edeadng(isd,inode)=0.0;
      }
    }
  }

  return edeadng;
} // end of DRT:Elements:Fluid2:f2_getbodyforce


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
  static Epetra_SerialDenseMatrix 	xji(2,2);   // inverse of jacobian matrix


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
    Epetra_SerialDenseVector& funct,
    vector<double>&           tau,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const vector<double>&     edeadng,
    const double&             fac,
    const double&             visc,
    const int&                iel,
    double                    timefac
    )
{
/*========================= further variables =========================*/

  Epetra_SerialDenseMatrix  viscs2(2,2*iel);  /* viscous term incluiding 2nd derivatives         */
  vector<double>            conv_c(iel); 	    /* linearisation of convect, convective part       */
  vector<double>            conv_g(iel); 	    /* linearisation of convect, grid part             */
  Epetra_SerialDenseMatrix  conv_r(2,2*iel);  /* linearisation of convect, reactive part         */
  vector<double>            div(2*iel);  	    /* divergence of u or v                            */
  Epetra_SerialDenseMatrix  ugradv(iel,2*iel);/* linearisation of u * grad v                     */
  static vector<double>            conv_old(2);      /* convective term evalaluated with old velocities */
  //vector<double>            conv_g_old(2);
  static vector<double>            visc_old(2); 	    /* viscous term evaluated with old velocities      */
  static vector<double>            rhsint(2);   	    /* total right hand side terms at int.-point       */

// stabilisation parameter
  const double tau_M  = tau[0]*fac;
  const double tau_Mp = tau[0]*fac;
  const double tau_C  = tau[2]*fac;

// integration factors and coefficients of single terms
  const double time2nue   = timefac * 2.0 * visc;
  const double timetauM   = timefac * tau_M;
  const double timetauMp  = timefac * tau_Mp;

  const double ttimetauM  = timefac * timetauM;
  const double ttimetauMp = timefac * timetauMp;
  const double timefacfac = timefac * fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
  rhsint[0] = histvec[0]+ edeadng[0]*timefac;;
  rhsint[1] = histvec[1]+ edeadng[1]*timefac;;

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
    Epetra_SerialDenseVector& funct,
    vector<double>&           tau,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const vector<double>&     edeadng,
    const double&             fac,
    const double&             visc,
    const int&                iel
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
double timefac = 1.0;

rhsint[0] = histvec[0] + edeadng[0]*timefac;
rhsint[1] = histvec[1] + edeadng[1]*timefac;

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

// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::Elements::Fluid2::is_higher_order_element(
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
int DRT::Elements::Fluid2Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
