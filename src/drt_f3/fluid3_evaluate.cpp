/*!----------------------------------------------------------------------
\file fluid3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "fluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::Elements::Fluid3::ActionType act = Fluid3::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_fluid_systemmat_and_residual")
  	act = Fluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_genalpha_sysmat")
  	act = Fluid3::calc_fluid_genalpha_sysmat;
  else if (action == "calc_fluid_genalpha_residual")
  	act = Fluid3::calc_fluid_genalpha_residual;
  else if (action == "calc_fluid_beltrami_error")
  	act = Fluid3::calc_fluid_beltrami_error;
  else dserror("Unknown type of action for Fluid3");

  // get the material
  MATERIAL* actmat = &(mat[material_-1]);

  switch(act)
  {
      case calc_fluid_systemmat_and_residual:
      {
        // need current velocity and history vector
        RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (trial)");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("old solution data for rhs");
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
        vector<double> myvelnp(3*numnode);
        vector<double> myvhist(3*numnode);

        for (int i=0;i<numnode;++i)
        {
          myvelnp[0+(i*3)]=my_vel_pre_np[0+(i*4)];
          myvelnp[1+(i*3)]=my_vel_pre_np[1+(i*4)];
          myvelnp[2+(i*3)]=my_vel_pre_np[2+(i*4)];

          myprenp[i]=my_vel_pre_np[3+(i*4)];

          myvhist[0+(i*3)]=myhist[0+(i*4)];
          myvhist[1+(i*3)]=myhist[1+(i*4)];
          myvhist[2+(i*3)]=myhist[2+(i*4)];
        }

        // calculate element coefficient matrix and rhs
        f3_sys_mat(lm,myvelnp,myprenp,myvhist,mydispnp,mygridv,&elemat1,&elevec1,actmat,params);

        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        params.set("density", actmat->m.fluid->density);

        /* the following has to be checked again !!! */
        // use local variables instead of directly write into elemat1, elevec1.
        // this speeds up computations by 3%-5%
        //Epetra_SerialDenseVector  eforce(4*numnode);      	// rhs vector
        //Epetra_SerialDenseMatrix 	estif(4*numnode,4*numnode); 	// element coefficient matrix

        // calculate element coefficient matrix and rhs
        //f3_sys_mat(lm,myvelnp,myprenp,myvhist,&estif,&eforce,actmat,params);

        // copy values
        //elemat1 = estif;
        //elevec1 = eforce;


// outputs for debugging

// if (Id()==10 || Id()==21)
        {
          //printf("Element %5d\n",Id());
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
        } // end of debug part

      }
      break;
      case calc_fluid_beltrami_error:
      {
        // add error only for elements which are not ghosted
        if(this->Owner() == discretization.Comm().MyPID())
        {

          // need current velocity and history vector
          RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
          if (vel_pre_np==null) dserror("Cannot get state vectors 'velnp'");

          // extract local values from the global vectors
          vector<double> my_vel_pre_np(lm.size());
          DRT::Utils::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

          // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
          int numnode = NumNode();
          vector<double> myprenp(numnode);
          vector<double> myvelnp(3*numnode);

          for (int i=0;i<numnode;++i)
          {
            myvelnp[0+(i*3)]=my_vel_pre_np[0+(i*4)];
            myvelnp[1+(i*3)]=my_vel_pre_np[1+(i*4)];
            myvelnp[2+(i*3)]=my_vel_pre_np[2+(i*4)];

            myprenp[i]=my_vel_pre_np[3+(i*4)];
          }

          // integrate beltrami error
          f3_int_beltrami_err(myvelnp,myprenp,actmat,params);
        }
      }
      break;
      case calc_fluid_genalpha_sysmat:
      {
        // --------------------------------------------------
        // extract velocities, pressure and accelerations from the
        // global distributed vectors

        // velocity and pressure values (current iterate, n+1)
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

        // velocities    (intermediate time step, n+alpha_F)
        RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

        // accelerations (intermediate time step, n+alpha_M)
        RefCountPtr<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");


        if (velnp==null || velaf==null || accam==null)
        {
          dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
        }

        // extract local values from the global vectors
        vector<double> my_velnp(lm.size());
        DRT::Utils::ExtractMyValues(*velnp,my_velnp,lm);

        vector<double> my_velaf(lm.size());
        DRT::Utils::ExtractMyValues(*velaf,my_velaf,lm);

        vector<double> my_accam(lm.size());
        DRT::Utils::ExtractMyValues(*accam,my_accam,lm);

        // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
        // Additionally only the 'velocity' components of my_velaf
        // and my_accam are important!
        int numnode = NumNode();
        vector<double> myprenp(numnode);
        vector<double> myvelnp(3*numnode);
        vector<double> myvelaf(3*numnode);
        vector<double> myaccam(3*numnode);

        for (int i=0;i<numnode;++i)
        {
          myvelnp[0+(i*3)] = my_velnp[0+(i*4)];
          myvelnp[1+(i*3)] = my_velnp[1+(i*4)];
          myvelnp[2+(i*3)] = my_velnp[2+(i*4)];

          myprenp[  (  i)] = my_velnp[3+(i*4)];

          myvelaf[0+(i*3)] = my_velaf[0+(i*4)];
          myvelaf[1+(i*3)] = my_velaf[1+(i*4)];
          myvelaf[2+(i*3)] = my_velaf[2+(i*4)];

          myaccam[0+(i*3)] = my_accam[0+(i*4)];
          myaccam[1+(i*3)] = my_accam[1+(i*4)];
          myaccam[2+(i*3)] = my_accam[2+(i*4)];
        }

        // --------------------------------------------------
        // calculate element coefficient matrix
        f3_genalpha_sys_mat(lm,
                            myvelnp,
                            myprenp,
                            myvelaf,
                            myaccam,
                            &elemat1,
                            actmat,
                            params);
      }
      break;
      case calc_fluid_genalpha_residual:
      {
        // --------------------------------------------------
        // extract velocities, pressure and accelerations from the
        // global distributed vectors

        // velocity and pressure values (current iterate, n+1)
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

        // velocities    (intermediate time step, n+alpha_F)
        RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

        // accelerations (intermediate time step, n+alpha_M)
        RefCountPtr<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");


        if (velnp==null || velaf==null || accam==null)
        {
          dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
        }

        // extract local values from the global vectors
        vector<double> my_velnp(lm.size());
        DRT::Utils::ExtractMyValues(*velnp,my_velnp,lm);

        vector<double> my_velaf(lm.size());
        DRT::Utils::ExtractMyValues(*velaf,my_velaf,lm);

        vector<double> my_accam(lm.size());
        DRT::Utils::ExtractMyValues(*accam,my_accam,lm);

        // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
        // Additionally only the 'velocity' components of my_velaf
        // and my_accam are important!
        int numnode = NumNode();
        vector<double> myprenp(numnode);
        vector<double> myvelnp(3*numnode);
        vector<double> myvelaf(3*numnode);
        vector<double> myaccam(3*numnode);

        for (int i=0;i<numnode;++i)
        {
          myvelnp[0+(i*3)] = my_velnp[0+(i*4)];
          myvelnp[1+(i*3)] = my_velnp[1+(i*4)];
          myvelnp[2+(i*3)] = my_velnp[2+(i*4)];

          myprenp[  (  i)] = my_velnp[3+(i*4)];

          myvelaf[0+(i*3)] = my_velaf[0+(i*4)];
          myvelaf[1+(i*3)] = my_velaf[1+(i*4)];
          myvelaf[2+(i*3)] = my_velaf[2+(i*4)];

          myaccam[0+(i*3)] = my_accam[0+(i*4)];
          myaccam[1+(i*3)] = my_accam[1+(i*4)];
          myaccam[2+(i*3)] = my_accam[2+(i*4)];
        }

        // --------------------------------------------------
        // calculate element right hand side
        f3_genalpha_rhs(lm,
                        myvelnp,
                        myprenp,
                        myvelaf,
                        myaccam,
                        elevec1,
                        actmat,
                        params);

      }
      break;
      default:
        dserror("Unknown type of action for Fluid3");
  } // end of switch(act)

  return 0;
} // end of DRT::Elements::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
  |  calculate system matrix and rhs (private)                g.bau 03/07|
  *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_sys_mat(vector<int>&              lm,
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
    Epetra_SerialDenseMatrix xyze(3,iel);

    // get node coordinates
    for (int i=0;i<iel;i++)
    {
      xyze(0,i)=Nodes()[i]->X()[0];
      xyze(1,i)=Nodes()[i]->X()[1];
      xyze(2,i)=Nodes()[i]->X()[2];
    }

    if (is_ale_)
    {
      for (int i=0;i<iel;i++)
      {
        xyze(0,i) += edispnp[4*i];
        xyze(1,i) += edispnp[4*i+1];
        xyze(2,i) += edispnp[4*i+2];
      }
    }

    // dead load in element nodes
    double time = params.get("total time",-1.0);
    Epetra_SerialDenseMatrix bodyforce(3,iel);
    this->f3_getbodyforce(bodyforce,time,iel,params);

    /*---------------------------------------------- get viscosity ---*/
    // check here, if we really have a fluid !!
    if(material->mattyp != m_fluid) dserror("Material law is not of type m_fluid.");
    const double  visc = material->m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    // USFEM stabilization is default. No switch here at the moment.

    /*----------------------------------------- declaration of variables ---*/
    vector<double> 		funct(iel);
    Epetra_SerialDenseMatrix 	deriv(3,iel);
    Epetra_SerialDenseMatrix 	deriv2(6,iel);
    Epetra_SerialDenseMatrix 	xjm(3,3);
    Epetra_SerialDenseMatrix 	vderxy(3,3);
    vector<double> 		pderxy(3);
    Epetra_SerialDenseMatrix 	vderxy2(3,6);
    Epetra_SerialDenseMatrix 	derxy(3,iel);
    Epetra_SerialDenseMatrix 	derxy2(6,iel);
    vector<double>              edeadng(3);
    Epetra_SerialDenseMatrix 	wa1(100,100);  // working matrix used as dummy
    vector<double>    		histvec(3); /* history data at integration point              */
    double         		hk;
    double         		val, strle;
    vector<double>         	velino(3); /* normed velocity at element centre */
    double         		det, vol;
    FLUID_DATA            	data;
    double         		e1, e2, e3;
    double         		facr=0.0, facs=0.0, fact=0.0;
    double         		mk=0.0;
    vector<double>     		velint(3);
    double 			timefac;
    vector<double>              tau(3); // stab parameters

    /*------------------------------------------------------- initialise ---*/
    // gaussian points
    f3_integration_points(data);
    timefac=params.get<double>("time constant for integration",0.0);

    // get control parameter
    bool is_stationary = params.get<bool>("using stationary formulation",false);

    /*---------------------- shape functions and derivs at element center --*/
    switch(iel)
    {
        case 8: case 20: case 27:   /* --> hex - element */
          e1   = data.qxg[0][0];
          facr = data.qwgt[0][0];
          e2   = data.qxg[0][0];
          facs = data.qwgt[0][0];
          e3   = data.qxg[0][0];
          fact = data.qwgt[0][0];

          f3_shape_function(funct,deriv,wa1,e1,e2,e3,iel,2); //wa1 as dummy for not wanted second derivatives
          break;
        case 4: case 10:   /* --> tet - element */
          e1   = data.txgr[0][0];
          facr = data.twgt[0][0];
          e2   = data.txgs[0][0];
          facs = ONE;
          e3   = data.txgs[0][0];
          fact = ONE;
          f3_shape_function(funct,deriv,wa1,e1,e2,e3,iel,2); //wa1 as dummy for not wanted second derivatives
          break;
        default:
          dserror("type unknown!\n");
    } /*end switch(iel) */

/*------------------------------- get element type constant for tau ---*/
    switch(iel)
    {
        case 4:
        case 8:
          mk = 0.333333333333333333333;
          break;
        case 20:
        case 27:
        case 10:
          mk = 0.083333333333333333333;
          break;
        default: dserror("type unknown!\n");
    }
/*--------------------------------- get velocities at element center ---*/
    for (int i=0;i<3;i++)
    {
      velint[i]=ZERO;
      for (int j=0;j<iel;j++)
      {
        velint[i] += funct[j]*evelnp[i+(3*j)];
      }
    } //end loop over i

    {
      double vel_norm, re1, re2, xi1, xi2, re, xi;

      /*------------------------------ get Jacobian matrix and determinant ---*/
      f3_jaco(xyze,deriv,xjm,&det,iel);
      vol=facr*facs*fact*det;

      /* get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)*/
      hk = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);

      /*------------------------------------------------- get streamlength ---*/
      f3_gder(derxy,deriv,xjm,det,iel);
      val = ZERO;


      /* get velocity norm */
      vel_norm=sqrt( velint[0]*velint[0]
                     + velint[1]*velint[1]
                     + velint[2]*velint[2]);
      if(vel_norm>=EPS6)
      {
        velino[0] = velint[0]/vel_norm;
        velino[1] = velint[1]/vel_norm;
        velino[2] = velint[2]/vel_norm;
      }
      else
      {
        velino[0] = ONE;
        velino[1] = ZERO;
        velino[2] = ZERO;
      }
      for (int i=0;i<iel;i++) /* loop element nodes */
      {
        val += FABS(velino[0]*derxy(0,i) \
                    +velino[1]*derxy(1,i) \
                    +velino[2]*derxy(2,i));
      } /* end of loop over elements */
      strle=TWO/val;

      if (is_stationary == false)
      {// stabilization parameters for instationary case (default)

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


      re1 =/* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
      re2 = mk * vel_norm * strle / /* *1.0 */(2.0 * visc);    /* convective : viscous forces */

      xi1 = DMAX(re1,1.0);
      xi2 = DMAX(re2,1.0);

      tau[0] = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

      /*------------------------------------------------------compute tau_Mp ---*/
      /* stability parameter definition according to Franca and Valentin (2000)
       *                                    and Barrenechea and Valentin (2002) */
      re1 = /* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
      re2 = mk * vel_norm * hk / /* *1.0 */(2.0 * visc);     /* convective : viscous forces */

      xi1 = DMAX(re1,1.0);
      xi2 = DMAX(re2,1.0);

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


      tau[1] = DSQR(hk) / (DSQR(hk) * xi1 + (/* 2.0*/ 4.0 * timefac * visc/mk) * xi2);

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
      xi2 = DMIN(re2,1.0);

      tau[2] = vel_norm * hk * 0.5 * xi2 /timefac;
      }
      else
      {// stabilization parameters for stationary case

	/*----------------------------------------------------- compute tau_Mu ---*/
	re = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
 	xi = DMAX(re,1.0);

	tau[0] = (DSQR(strle)*mk)/(4.0*visc*xi);

      	/*------------------------------------------------------compute tau_Mp ---*/
	re = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
      	xi = DMAX(re,1.0);

      	tau[1] = (DSQR(hk)*mk)/(4.0*visc*xi);

	/*------------------------------------------------------ compute tau_C ---*/
      	xi = DMIN(re,1.0);
	tau[2] = 0.5*vel_norm*hk*xi;
      }
    }
    /*----------------------------------------------------------------------*/
    // end of old f3_caltau function
    /*----------------------------------------------------------------------*/

    // integration loop for one Fluid3 element using USFEM

    int       intc=0;      /* "integration case" for tri for further infos
                              see f2_inpele.c and f2_intg.c                 */
    int       nir=0;       /* number of integration nodesin r direction     */
    int       nis=0;       /* number of integration nodesin s direction     */
    int       nit=0;       /* number of integration nodesin t direction      */
    int       ihoel=0;     /* flag for higher order elements                 */
    int       icode=2;     /* flag for eveluation of shape functions         */
    double    fac;         /* total integration factor */
    double    press;
    vector<double>    gridvelint(3); /* grid velocity                       */
    vector<double>    gradp(3);      /* pressure gradient at integration point         */


    switch (iel)
    {
        case 8: case 20: case 27:  /* --> hex - element */
          icode   = 3;
          ihoel   = 1;
          /* initialise integration */
          nir = ngp_[0];
          nis = ngp_[1];
          nit = ngp_[2];
          intc= 0;
          break;
        case 10: /* --> tet - element */
          icode   = 3;
          ihoel   = 1;
          /* do NOT break at this point!!! */
        case 4:    /* initialise integration */
          nir  = ngp_[0]; // for tets in ngp_[0] the number of gauss points is stored !
          nis  = 1;
          nit  = 1;
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
        for (int lt=0;lt<nit;lt++)
        {
          /*------------- get values of  shape functions and their derivatives ---*/
          switch(iel)
          {
              case 8: case 20: case 27:   /* --> hex - element */
                e1   = data.qxg[lr][nir-1];
                facr = data.qwgt[lr][nir-1];
                e2   = data.qxg[ls][nis-1];
                facs = data.qwgt[ls][nis-1];
                e3   = data.qxg[lt][nit-1];
                fact = data.qwgt[lt][nit-1];
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              case 4: case 10:   /* --> tet - element */
                e1   = data.txgr[lr][intc];
                facr = data.twgt[lr][intc];
                e2   = data.txgs[lr][intc];
                facs = ONE;
                e3   = data.txgt[lr][intc];
                fact = ONE;
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              default:
                facr = facs = fact = 0.0;
                e1 = e2 = e3 = 0.0;
                dserror("typ unknown!");
          } /* end switch (iel) */

          /*----------------------------------------- compute Jacobian matrix */
          f3_jaco(xyze,deriv,xjm,&det,iel);
          fac = facr*facs*fact*det;

          /*---------------------------------------- compute global derivates */
          f3_gder(derxy,deriv,xjm,det,iel);

          /*--------------------------------- compute second global derivative */
          if (ihoel!=0)
          {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

            /*------calculate 2nd velocity derivatives at integration point */
            // former f3_vder2(vderxy2,derxy2,evelnp,iel);
            for (int i=0;i<6;i++)
            {
              vderxy2(0,i)=ZERO;
              vderxy2(1,i)=ZERO;
              vderxy2(2,i)=ZERO;
              for (int j=0;j<iel;j++)
              {
                vderxy2(0,i) += derxy2(i,j)*evelnp[0+(3*j)];
                vderxy2(1,i) += derxy2(i,j)*evelnp[1+(3*j)];
                vderxy2(2,i) += derxy2(i,j)*evelnp[2+(3*j)];
              } /* end of loop over j */
            } /* end of loop over i */
          }

          /*---------------------- get velocities (n+g,i) at integration point */
          // expression for f3_veci(velint,funct,evelnp,iel);
          for (int i=0;i<3;i++)
          {
            velint[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              velint[i] += funct[j]*evelnp[i+(3*j)];
            }
          } //end loop over i

          /*---------------- get history data (n,i) at integration point ---*/
          //expression for f3_veci(histvec,funct,evhist,iel);
          for (int i=0;i<3;i++)
          {
            histvec[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              histvec[i] += funct[j]*evhist[i+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */

          /*----------- get velocity (np,i) derivatives at integration point */
          // expression for f3_vder(vderxy,derxy,evelnp,iel);
          for (int i=0;i<3;i++)
          {
            vderxy(0,i)=ZERO;
            vderxy(1,i)=ZERO;
            vderxy(2,i)=ZERO;
            for (int j=0;j<iel;j++)
            {
              vderxy(0,i) += derxy(i,j)*evelnp[0+(3*j)];
              vderxy(1,i) += derxy(i,j)*evelnp[1+(3*j)];
              vderxy(2,i) += derxy(i,j)*evelnp[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */

          /*--------------------- get grid velocity at integration point ---*/
          if (is_ale_)
          {
            for (int i=0; i<3; i++)
            {
              gridvelint[i] = 0.;
              for (int j=0; j<iel; j++)
              {
                gridvelint[i] += derxy(i,j)*egridv[i+(4*j)];
              }
            }
          }
          else
          {
            gridvelint[0] = 0.0;
            gridvelint[1] = 0.0;
            gridvelint[2] = 0.0;
          }

          /*------------------------------------- get pressure gradients ---*/
          gradp[0] = gradp[1] = gradp[2] = 0.0;

          for (int i=0; i<iel; i++)
          {
            gradp[0] += derxy(0,i) * eprenp[i];
            gradp[1] += derxy(1,i) * eprenp[i];
            gradp[2] += derxy(2,i) * eprenp[i];
          }

          press = 0;
          for (int i=0;i<iel;i++)
          {
            press += funct[i]*eprenp[i];
          }

          /*--------------------------------- get bodyforce in gausspoint---*/
          for (int dim=0;dim<3;dim++)
          {
            edeadng[dim] = 0;

            for (int i=0;i<iel;i++)
            {
              edeadng[dim]+= bodyforce(dim,i)*funct[i];
            }
          }

          /*-------------- perform integration for entire matrix and rhs ---*/
          if(is_stationary==false)
            f3_calmat(*sys_mat,*residual,velint,histvec,gridvelint,
                    press,vderxy,vderxy2,gradp,funct,tau,
                    derxy,derxy2,edeadng,fac,visc,iel,params);
	  else
	    f3_calmat_stationary(*sys_mat,*residual,velint,histvec,gridvelint,
                    press,vderxy,vderxy2,gradp,funct,tau,
                    derxy,derxy2,edeadng,fac,visc,iel,params);


        } /* end of loop over integration points lt*/
      } /* end of loop over integration points ls */
    } /* end of loop over integration points lr */

  return;
} // DRT::Elements::Fluid3::f3_sys_mat



/*----------------------------------------------------------------------*
  |  calculate system matrix for a generalised alpha time integration   |
  |                            (private)                     gammi 06/07|
  *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_genalpha_sys_mat(
  vector<int>&              lm,
  vector<double>&           myvelnp,
  vector<double>&           myprenp,
  vector<double>&           myvelaf,
  vector<double>&           myaccam,
  Epetra_SerialDenseMatrix* elemat,
  struct _MATERIAL*         material,
  ParameterList& 	    params)
{

  if(!is_ale_)
  {
    /*---------------------------------------------------- set element data */
    const int iel = NumNode();
    Epetra_SerialDenseMatrix xyze(3,iel);

    // get node coordinates
    for(int i=0;i<iel;i++)
    {
      xyze(0,i)=Nodes()[i]->X()[0];
      xyze(1,i)=Nodes()[i]->X()[1];
      xyze(2,i)=Nodes()[i]->X()[2];
    }

    /*-------- dead load in element nodes, evaluated at t_(n+alpha_F) */
    Epetra_SerialDenseMatrix bodyforce(3,iel);
    {
      double acttime    = params.get<double>("time");
      double dt         = params.get<double>("dt");
      double alphaF     = params.get<double>("alpha_F");

      //         n+alpha_F     n+1
      //        t          = t     - (1-alpha_F) * dt

      double timealphaF = acttime-(1-alphaF)*dt;

      this->f3_getbodyforce(bodyforce,timealphaF,iel,params);
    }

    /*---------------------------------------------- get viscosity ---*/
    // check here, if we really have a fluid !!
    if(material->mattyp != m_fluid)
    {
      dserror("Material law is not of type m_fluid.");
    }
    const double  visc = material->m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    vector<double>   tau(3); // stab parameters
    {
      int version =2; // evaluate stabilisation parameter for genalpha
                      // time integration

      f3_calc_stabpar(tau,iel,xyze,myvelaf,visc,params,version);
    }

    /*-------------------------------- initialise the integration ----*/
    FLUID_DATA            	data;
    f3_integration_points(data);


    int       intc=0;      /* "integration case" for tri for further infos
                              see f2_inpele.c and f2_intg.c                 */
    int       nir=0;       /* number of integration nodesin r direction     */
    int       nis=0;       /* number of integration nodesin s direction     */
    int       nit=0;       /* number of integration nodesin t direction      */
    int       ihoel=0;     /* flag for higher order elements                 */
    int       icode=2;     /* flag for eveluation of shape functions         */

    switch (iel)
    {
        case 8: case 20: case 27:  /* --> hex - element */
          icode   = 3;
          ihoel   = 1;
          /* initialise integration */
          nir = ngp_[0];
          nis = ngp_[1];
          nit = ngp_[2];
          intc= 0;
          break;
        case 10: /* --> tet - element */
          icode   = 3;
          ihoel   = 1;
          /* do NOT break at this point!!! */
        case 4:    /* initialise integration */
          nir  = ngp_[0]; // for tets in ngp_[0] the number of gauss points is stored !
          nis  = 1;
          nit  = 1;
          intc = ngp_[1];
          break;
        default:
          dserror("typ unknown!");
    } // end switch (iel) //


    /*------------------------------------------------------------------*
     |              start loop over integration points                  |
     *------------------------------------------------------------------*/
    for (int lr=0;lr<nir;lr++)
    {
      for (int ls=0;ls<nis;ls++)
      {
        for (int lt=0;lt<nit;lt++)
        {
          /*------------------- declaration of gauss point variables ---*/

          // shape functions and derivatives
          vector<double> 		funct (iel);
          Epetra_SerialDenseMatrix 	deriv (3,iel);
          Epetra_SerialDenseMatrix 	deriv2(6,iel);
          Epetra_SerialDenseMatrix 	derxy (3,iel);
          Epetra_SerialDenseMatrix 	derxy2(6,iel);
          Epetra_SerialDenseMatrix 	xjm   (3,3);

          // intermediate accelerations (n+alpha_M)
          vector<double>     		accintam (3);
          // intermediate velocities    (n+alpha_F) and its derivatives
          vector<double>     		velintaf (3);
          Epetra_SerialDenseMatrix 	vderxyaf (3,3);
          Epetra_SerialDenseMatrix 	vderxyaf2(3,6);
          // new velocities for continuity equation and its derivatives
          vector<double>     		velintnp (3);
          Epetra_SerialDenseMatrix 	vderxynp (3,3);
          // new pressure and its derivatives
          double                        prenp;
          vector<double> 		pderxynp(3);

          // dead load
          vector<double>                edeadaf(3);

          // working matrix used as dummy
          Epetra_SerialDenseMatrix 	wa1(100,100);

          double         		det;

          // integration data
          // -> gauss point coordinates
          double         		e1, e2, e3;
          // -> gauss weights
          double         		facr=0.0, facs=0.0, fact=0.0;
          // -> total integration factor
          double                        fac;

          /*---- get values of  shape functions and their derivatives ---*/
          switch(iel)
          {
              case 8: case 20: case 27:   /* --> hex - element */
                e1   = data.qxg [lr][nir-1];
                facr = data.qwgt[lr][nir-1];
                e2   = data.qxg [ls][nis-1];
                facs = data.qwgt[ls][nis-1];
                e3   = data.qxg [lt][nit-1];
                fact = data.qwgt[lt][nit-1];
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              case 4: case 10:   /* --> tet - element */
                e1   = data.txgr[lr][intc];
                facr = data.twgt[lr][intc];
                e2   = data.txgs[lr][intc];
                facs = ONE;
                e3   = data.txgt[lr][intc];
                fact = ONE;
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              default:
                dserror("typ unknown!");
          } /* end switch (iel) */

          /*----------------------------------- compute Jacobian matrix */
          f3_jaco(xyze,deriv,xjm,&det,iel);

          // set total integration factor
          fac = facr*facs*fact*det;

          /*---------------------------------- compute global derivates */
          f3_gder(derxy,deriv,xjm,det,iel);

          /*-------------------------- compute second global derivative */
          if (ihoel!=0)
          {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);
          }


          /*-- get intermediate accelerations (n+1,i)  at integration
                                                                  point */
          for (int i=0;i<3;i++)
          {
            accintam[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              accintam[i] += funct[j]*myaccam[i+(3*j)];
            }
          } //end loop over i


          /*-------------- get velocities (n+1,i)  at integration point */
          for (int i=0;i<3;i++)
          {
            velintnp[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              velintnp[i] += funct[j]*myvelnp[i+(3*j)];
            }
          } //end loop over i


          /*--------- get velocities (n+alpha_F,i) at integration point */
          for (int i=0;i<3;i++)
          {
            velintaf[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              velintaf[i] += funct[j]*myvelaf[i+(3*j)];
            }
          } //end loop over i


          /*----- get velocity (n+1,i) derivatives at integration point */
          for (int i=0;i<3;i++)
          {
            vderxynp(0,i)=ZERO;
            vderxynp(1,i)=ZERO;
            vderxynp(2,i)=ZERO;
            for (int j=0;j<iel;j++)
            {
              vderxynp(0,i) += derxy(i,j)*myvelnp[0+(3*j)];
              vderxynp(1,i) += derxy(i,j)*myvelnp[1+(3*j)];
              vderxynp(2,i) += derxy(i,j)*myvelnp[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*----------- get velocity (n+alpha_F,i) derivatives at
                                                      integration point */
          for (int i=0;i<3;i++)
          {
            vderxyaf(0,i)=ZERO;
            vderxyaf(1,i)=ZERO;
            vderxyaf(2,i)=ZERO;
            for (int j=0;j<iel;j++)
            {
              vderxyaf(0,i) += derxy(i,j)*myvelaf[0+(3*j)];
              vderxyaf(1,i) += derxy(i,j)*myvelaf[1+(3*j)];
              vderxyaf(2,i) += derxy(i,j)*myvelaf[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*------calculate 2nd velocity derivatives at integration
                                                     point (n+alpha_F,i)*/
          if (ihoel!=0)
          {
            for (int i=0;i<6;i++)
            {
              vderxyaf2(0,i)=ZERO;
              vderxyaf2(1,i)=ZERO;
              vderxyaf2(2,i)=ZERO;
              for (int j=0;j<iel;j++)
              {
                vderxyaf2(0,i) += derxy2(i,j)*myvelaf[0+(3*j)];
                vderxyaf2(1,i) += derxy2(i,j)*myvelaf[1+(3*j)];
                vderxyaf2(2,i) += derxy2(i,j)*myvelaf[2+(3*j)];
              } /* end of loop over j */
            } /* end of loop over i */
          }


          /*--------------------------------- get pressure at time (n+1) ---*/
          prenp = 0;
          for (int i=0;i<iel;i++)
          {
            prenp += funct[i]*myprenp[i];
          }

          /*------------------------ get pressure gradient at time (n+1) ---*/
          pderxynp[0] = pderxynp[1] = pderxynp[2] = 0.0;

          for (int i=0; i<iel; i++)
          {
            pderxynp[0] += derxy(0,i) * myprenp[i];
            pderxynp[1] += derxy(1,i) * myprenp[i];
            pderxynp[2] += derxy(2,i) * myprenp[i];
          }


          /*----------- get bodyforce in gausspoint, time (n+alpha_F)---*/
          for (int dim=0;dim<3;dim++)
          {
            edeadaf[dim] = 0;

            for (int i=0;i<iel;i++)
            {
              edeadaf[dim]+= bodyforce(dim,i)*funct[i];
            }
          }

          /*--- assemble all contributions into to the element matrix --*/
          f3_genalpha_calmat(*elemat,
                             accintam,
                             velintaf,
                             vderxyaf,
                             vderxyaf2,
                             velintnp,
                             vderxynp,
                             prenp,
                             pderxynp,
                             edeadaf,
                             funct,
                             derxy,
                             derxy2,
                             tau,
                             fac,
                             visc,
                             iel,
                             params);

        } /* end of loop over integration points lt*/
      } /* end of loop over integration points ls */
    } /* end of loop over integration points lr */

  }
  return;
} // end of DRT:Elements:Fluid3:f3_genalpha_sys_mat

/*----------------------------------------------------------------------*
  |  calculate residual (rhs) for a generalised alpha time integration  |
  |                            (private)                     gammi 06/07|
  *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_genalpha_rhs(
  vector<int>&              lm,
  vector<double>&           myvelnp,
  vector<double>&           myprenp,
  vector<double>&           myvelaf,
  vector<double>&           myaccam,
  Epetra_SerialDenseVector& elevec,
  struct _MATERIAL*         material,
  ParameterList& 	    params)
{

  if(!is_ale_)
  {
    /*---------------------------------------------------- set element data */
    const int iel = NumNode();
    Epetra_SerialDenseMatrix xyze(3,iel);

    // get node coordinates
    for(int i=0;i<iel;i++)
    {
      xyze(0,i)=Nodes()[i]->X()[0];
      xyze(1,i)=Nodes()[i]->X()[1];
      xyze(2,i)=Nodes()[i]->X()[2];
    }

    /*-------- dead load in element nodes, evaluated at t_(n+alpha_F) */
    Epetra_SerialDenseMatrix bodyforce(3,iel);
    {
      double acttime    = params.get<double>("time");
      double dt         = params.get<double>("dt");
      double alphaF     = params.get<double>("alpha_F");

      //         n+alpha_F     n+1
      //        t          = t     - (1-alpha_F) * dt

      double timealphaF = acttime-(1-alphaF)*dt;

      this->f3_getbodyforce(bodyforce,timealphaF,iel,params);
    }

    /*---------------------------------------------- get viscosity ---*/
    // check here, if we really have a fluid !!
    if(material->mattyp != m_fluid)
    {
      dserror("Material law is not of type m_fluid.");
    }
    const double  visc = material->m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    vector<double>   tau(3); // stab parameters
    {
      int version =2; // evaluate stabilisation parameter for genalpha
                      // time integration

      f3_calc_stabpar(tau,iel,xyze,myvelaf,visc,params,version);
    }

    /*-------------------------------- initialise the integration ----*/
    FLUID_DATA            	data;
    f3_integration_points(data);


    int       intc=0;      /* "integration case" for tri for further infos
                              see f2_inpele.c and f2_intg.c                 */
    int       nir=0;       /* number of integration nodesin r direction     */
    int       nis=0;       /* number of integration nodesin s direction     */
    int       nit=0;       /* number of integration nodesin t direction      */
    int       ihoel=0;     /* flag for higher order elements                 */
    int       icode=2;     /* flag for eveluation of shape functions         */

    switch (iel)
    {
        case 8: case 20: case 27:  /* --> hex - element */
          icode   = 3;
          ihoel   = 1;
          /* initialise integration */
          nir = ngp_[0];
          nis = ngp_[1];
          nit = ngp_[2];
          intc= 0;
          break;
        case 10: /* --> tet - element */
          icode   = 3;
          ihoel   = 1;
          /* do NOT break at this point!!! */
        case 4:    /* initialise integration */
          nir  = ngp_[0]; // for tets in ngp_[0] the number of gauss points is stored !
          nis  = 1;
          nit  = 1;
          intc = ngp_[1];
          break;
        default:
          dserror("typ unknown!");
    } // end switch (iel) //


    /*------------------------------------------------------------------*
     |              start loop over integration points                  |
     *------------------------------------------------------------------*/
    for (int lr=0;lr<nir;lr++)
    {
      for (int ls=0;ls<nis;ls++)
      {
        for (int lt=0;lt<nit;lt++)
        {
          /*------------------- declaration of gauss point variables ---*/

          // shape functions and derivatives
          vector<double> 		funct (iel);
          Epetra_SerialDenseMatrix 	deriv (3,iel);
          Epetra_SerialDenseMatrix 	deriv2(6,iel);
          Epetra_SerialDenseMatrix 	derxy (3,iel);
          Epetra_SerialDenseMatrix 	derxy2(6,iel);
          Epetra_SerialDenseMatrix 	xjm   (3,3);

          // intermediate accelerations (n+alpha_M)
          vector<double>     		accintam (3);
          // intermediate velocities    (n+alpha_F) and its derivatives
          vector<double>     		velintaf (3);
          Epetra_SerialDenseMatrix 	vderxyaf (3,3);
          Epetra_SerialDenseMatrix 	vderxyaf2(3,6);
          // new velocities for continuity equation and its derivatives
          vector<double>     		velintnp (3);
          Epetra_SerialDenseMatrix 	vderxynp (3,3);
          // new pressure and its derivatives
          double                        prenp;
          vector<double> 		pderxynp(3);

          // dead load
          vector<double>                edeadaf(3);

          // working matrix used as dummy
          Epetra_SerialDenseMatrix 	wa1(100,100);

          double         		det;

          // integration data
          // -> gauss point coordinates
          double         		e1, e2, e3;
          // -> gauss weights
          double         		facr=0.0, facs=0.0, fact=0.0;
          // -> total integration factor
          double                        fac;

          /*---- get values of  shape functions and their derivatives ---*/
          switch(iel)
          {
              case 8: case 20: case 27:   /* --> hex - element */
                e1   = data.qxg [lr][nir-1];
                facr = data.qwgt[lr][nir-1];
                e2   = data.qxg [ls][nis-1];
                facs = data.qwgt[ls][nis-1];
                e3   = data.qxg [lt][nit-1];
                fact = data.qwgt[lt][nit-1];
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              case 4: case 10:   /* --> tet - element */
                e1   = data.txgr[lr][intc];
                facr = data.twgt[lr][intc];
                e2   = data.txgs[lr][intc];
                facs = ONE;
                e3   = data.txgt[lr][intc];
                fact = ONE;
                f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
                break;
              default:
                dserror("typ unknown!");
          } /* end switch (iel) */

          /*----------------------------------- compute Jacobian matrix */
          f3_jaco(xyze,deriv,xjm,&det,iel);

          // set total integration factor
          fac = facr*facs*fact*det;

          /*---------------------------------- compute global derivates */
          f3_gder(derxy,deriv,xjm,det,iel);

          /*-------------------------- compute second global derivative */
          if (ihoel!=0)
          {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);
          }


          /*-- get intermediate accelerations (n+am,i)  at integration
                                                                  point */
          for (int i=0;i<3;i++)
          {
            accintam[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              accintam[i] += funct[j]*myaccam[i+(3*j)];
            }
          } //end loop over i


          /*-------------- get velocities (n+1,i)  at integration point */
          for (int i=0;i<3;i++)
          {
            velintnp[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              velintnp[i] += funct[j]*myvelnp[i+(3*j)];
            }
          } //end loop over i


          /*--------- get velocities (n+alpha_F,i) at integration point */
          for (int i=0;i<3;i++)
          {
            velintaf[i]=ZERO;
            for (int j=0;j<iel;j++)
            {
              velintaf[i] += funct[j]*myvelaf[i+(3*j)];
            }
          } //end loop over i


          /*----- get velocity (n+1,i) derivatives at integration point */
          for (int i=0;i<3;i++)
          {
            vderxynp(0,i)=ZERO;
            vderxynp(1,i)=ZERO;
            vderxynp(2,i)=ZERO;
            for (int j=0;j<iel;j++)
            {
              vderxynp(0,i) += derxy(i,j)*myvelnp[0+(3*j)];
              vderxynp(1,i) += derxy(i,j)*myvelnp[1+(3*j)];
              vderxynp(2,i) += derxy(i,j)*myvelnp[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*----------- get velocity (n+alpha_F,i) derivatives at
                                                      integration point */
          for (int i=0;i<3;i++)
          {
            vderxyaf(0,i)=ZERO;
            vderxyaf(1,i)=ZERO;
            vderxyaf(2,i)=ZERO;
            for (int j=0;j<iel;j++)
            {
              vderxyaf(0,i) += derxy(i,j)*myvelaf[0+(3*j)];
              vderxyaf(1,i) += derxy(i,j)*myvelaf[1+(3*j)];
              vderxyaf(2,i) += derxy(i,j)*myvelaf[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*------calculate 2nd velocity derivatives at integration
                                                     point (n+alpha_F,i)*/
          if (ihoel!=0)
          {
            for (int i=0;i<6;i++)
            {
              vderxyaf2(0,i)=ZERO;
              vderxyaf2(1,i)=ZERO;
              vderxyaf2(2,i)=ZERO;
              for (int j=0;j<iel;j++)
              {
                vderxyaf2(0,i) += derxy2(i,j)*myvelaf[0+(3*j)];
                vderxyaf2(1,i) += derxy2(i,j)*myvelaf[1+(3*j)];
                vderxyaf2(2,i) += derxy2(i,j)*myvelaf[2+(3*j)];
              } /* end of loop over j */
            } /* end of loop over i */
          }


          /*--------------------------------- get pressure at time (n+1) ---*/
          prenp = 0;
          for (int i=0;i<iel;i++)
          {
            prenp += funct[i]*myprenp[i];
          }

          /*------------------------ get pressure gradient at time (n+1) ---*/
          pderxynp[0] = pderxynp[1] = pderxynp[2] = 0.0;

          for (int i=0; i<iel; i++)
          {
            pderxynp[0] += derxy(0,i) * myprenp[i];
            pderxynp[1] += derxy(1,i) * myprenp[i];
            pderxynp[2] += derxy(2,i) * myprenp[i];
          }


          /*----------- get bodyforce in gausspoint, time (n+alpha_F)---*/
          for (int dim=0;dim<3;dim++)
          {
            edeadaf[dim] = 0;

            for (int i=0;i<iel;i++)
            {
              edeadaf[dim]+= bodyforce(dim,i)*funct[i];
            }
          }

          /*--- assemble all contributions into to the element rhs --*/
          f3_genalpha_calrhs(elevec,
                             accintam,
                             velintaf,
                             vderxyaf,
                             vderxyaf2,
                             velintnp,
                             vderxynp,
                             prenp,
                             pderxynp,
                             edeadaf,
                             funct,
                             derxy,
                             derxy2,
                             tau,
                             fac,
                             visc,
                             iel,
                             params);

        } /* end of loop over integration points lt*/
      } /* end of loop over integration points ls */
    } /* end of loop over integration points lr */

  }
  return;
} // end of DRT:Elements:Fluid3:f3_genalpha_rhs


/*----------------------------------------------------------------------*
 |  evaluate the stabilisation parameter at the element center.         |
 |                               (private)                   gammi 06/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_calc_stabpar(
  vector<double>           &tau,
  int                       iel,
  Epetra_SerialDenseMatrix &xyze,
  vector<double>           &myvelnp,
  double                    visc,
  ParameterList            &params,
  int                       version
)
{

  switch(version)
  {
      case 0:
      {
        // STATIONARY FLOW PROBLEM
        // tau_M: Barrenechea, G.R. and Valentin, F.
        // tau_C: Wall
        dserror("code for this version of tau is still inline!");

      }
      break;
      case 1:
      {
        // INSTATIONARY FLOW PROBLEM, ONE-STEP THETA AND BDF2
        // tau_M: Barrenechea, G.R. and Valentin, F.
        // tau_C: Wall
        dserror("code for this version of tau is still inline!");

      }
      break;
      case 2:
      {
        // INSTATIONARY FLOW PROBLEM, GENERALISED ALPHA
        // tau_M: Barrenechea, G.R. and Valentin, F.
        // tau_C: Wall

        /*----------------------------- declaration of variables ---*/
        vector<double> 	   	    funct(iel);
        Epetra_SerialDenseMatrix    deriv(3,iel);
        Epetra_SerialDenseMatrix    xjm(3,3);
        Epetra_SerialDenseMatrix    derxy(3,iel);
        Epetra_SerialDenseMatrix    wa1(100,100);  // working matrix used as dummy
        double                      hk;
        double                      val, strle;
        vector<double>              velino(3); /* normed velocity at element centre */
        double                      det, vol;
        FLUID_DATA                  data;
        double                      e1, e2, e3;
        double                      facr=0.0, facs=0.0, fact=0.0;
        double                      mk=0.0;
        vector<double>              velint(3);
        double                      timefac;

        /*------------------------------------------------- initialise ---*/
        // gaussian points
        f3_integration_points(data);
        timefac=params.get<double>("dt",0.0);

        /*---------------------- shape functions and derivs at element center --*/
        switch(iel)
        {
            case 8: case 20: case 27:   /* --> hex - element */
              e1   = data.qxg[0][0];
              facr = data.qwgt[0][0];
              e2   = data.qxg[0][0];
              facs = data.qwgt[0][0];
              e3   = data.qxg[0][0];
              fact = data.qwgt[0][0];

              f3_shape_function(funct,deriv,wa1,e1,e2,e3,iel,2); //wa1 as dummy for not wanted second derivatives
              break;
            case 4: case 10:   /* --> tet - element */
              e1   = data.txgr[0][0];
              facr = data.twgt[0][0];
              e2   = data.txgs[0][0];
              facs = ONE;
              e3   = data.txgs[0][0];
              fact = ONE;
              f3_shape_function(funct,deriv,wa1,e1,e2,e3,iel,2); //wa1 as dummy for not wanted second derivatives
              break;
            default:
              dserror("type unknown!\n");
        } /*end switch(iel) */

        /*------------------------------- get element type constant for tau ---*/
        switch(iel)
        {
            case 4:
            case 8:
              mk = 0.333333333333333333333;
              break;
            case 20:
            case 27:
            case 10:
              mk = 0.083333333333333333333;
              break;
            default: dserror("type unknown!\n");
        }
        /*--------------------------------- get velocities at element center ---*/
        for (int i=0;i<3;i++)
        {
          velint[i]=ZERO;
          for (int j=0;j<iel;j++)
          {
            velint[i] += funct[j]*myvelnp[i+(3*j)];
          }
        } //end loop over i

        {
          double vel_norm, re1, re2, xi1, xi2;

          /*------------------------------ get Jacobian matrix and determinant ---*/
          f3_jaco(xyze,deriv,xjm,&det,iel);
          vol=facr*facs*fact*det;

          /* get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)*/
          hk = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);

          /*------------------------------------------------- get streamlength ---*/
          f3_gder(derxy,deriv,xjm,det,iel);
          val = ZERO;


          /* get velocity norm */
          vel_norm=sqrt( velint[0]*velint[0]
                         + velint[1]*velint[1]
                         + velint[2]*velint[2]);
          if(vel_norm>=EPS6)
          {
            velino[0] = velint[0]/vel_norm;
            velino[1] = velint[1]/vel_norm;
            velino[2] = velint[2]/vel_norm;
          }
          else
          {
            velino[0] = ONE;
            velino[1] = ZERO;
            velino[2] = ZERO;
          }
          for (int i=0;i<iel;i++) /* loop element nodes */
          {
            val += FABS(velino[0]*derxy(0,i)    \
                        +velino[1]*derxy(1,i)   \
                        +velino[2]*derxy(2,i));
          } /* end of loop over elements */
          strle=TWO/val;

          {// stabilization parameters for instationary case (default)

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


            re1 = 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous    : reactive forces */
            re2 = mk * vel_norm * strle / (2.0 * visc);      /* convective : viscous forces */

            xi1 = DMAX(re1,1.0);
            xi2 = DMAX(re2,1.0);

            tau[0] = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

            /*------------------------------------------------------compute tau_Mp ---*/
            /* stability parameter definition according to Franca and Valentin (2000)
             *                                    and Barrenechea and Valentin (2002) */
            re1 = 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous    : reactive forces */
            re2 = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces  */

            xi1 = DMAX(re1,1.0);
            xi2 = DMAX(re2,1.0);

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


            tau[1] = DSQR(hk) / (DSQR(hk) * xi1 + ( 4.0 * timefac * visc/mk) * xi2);

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
            xi2 = DMIN(re2,1.0);

            tau[2] = vel_norm * hk * 0.5 * xi2 /timefac;
          }

        }

      }
      break;
      default:
        dserror("unknown version of tau!");

  }

  return;
} // end of DRT:Elements:Fluid3:f3_calc_stabpar

/*----------------------------------------------------------------------*
 |  evaluate the element integration points (private)        g.bau 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_integration_points(struct _FLUID_DATA& data)
{
const double Q12  = ONE/TWO;
const double Q14  = ONE/FOUR;
const double Q16  = ONE/SIX;
const double Q124 = ONE/SIX/FOUR;
const double Q430 = FOUR/FIVE/SIX;
const double Q9120= NINE/FOUR/FIVE/SIX;

double  palpha;
double  pbeta;

/*----------------------------------------------------------------------*/
                                               /* initialize arrays */
  for (int i=0; i<MAXTINTP; i++)
  {
    for (int k=0; k<MAXTINTC; k++)
    {
       data.txgr[i][k] = ZERO;
       data.txgs[i][k] = ZERO;
       data.txgt[i][k] = ZERO;
       data.twgt[i][k] = ZERO;
    }
  }
  for (int i=0; i<MAXQINTP; i++)
  {
    for (int k=0; k<MAXQINTC; k++)
    {
       data.qxg[i][k] = ZERO;
       data.qwgt[i][k] = ZERO;
    }
  }

  palpha = (FIVE+THREE*sqrt(FIVE))/20.0;
  pbeta  = (FIVE-sqrt(FIVE))/20.0;

/*----------------------------------------------------------------------*
 |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L    ELEMENTS    |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 |     xg[i][j]                                                         |
 |    wgt[i][j]:  i+1 - actual number of gausspoint                     |
 |                j+1 - total number of gausspoints                     |
 *----------------------------------------------------------------------*/
/* coordinates for two gauss points */
      data.qxg[0][1]  =  -0.5773502691896;
      data.qxg[1][1]  =  -data.qxg[0][1]       ;
/* coordinates for three gauss points */
      data.qxg[0][2]  =  -0.7745966692415;
      data.qxg[2][2]  =  -data.qxg[0][2]       ;
/* coordinates for four gauss points */
      data.qxg[0][3]  =  -0.8611363115941;
      data.qxg[1][3]  =  -0.3399810435849;
      data.qxg[2][3]  =  -data.qxg[1][3]       ;
      data.qxg[3][3]  =  -data.qxg[0][3]       ;
/* coordinates for five gauss points */
      data.qxg[0][4]  =  -0.9061798459387;
      data.qxg[1][4]  =  -0.5384693101057;
      data.qxg[3][4]  =  -data.qxg[1][4]       ;
      data.qxg[4][4]  =  -data.qxg[0][4]       ;
/* coordinates for six gauss points */
      data.qxg[0][5]  =  -0.9324695142032;
      data.qxg[1][5]  =  -0.6612093864663;
      data.qxg[2][5]  =  -0.2386191860832;
      data.qxg[3][5]  =  -data.qxg[2][5]       ;
      data.qxg[4][5]  =  -data.qxg[1][5]       ;
      data.qxg[5][5]  =  -data.qxg[0][5]       ;

/* weights for one gauss points */
      data.qwgt[0][0] =  TWO             ;
/* weights for two gauss points */
      data.qwgt[0][1] =  ONE             ;
      data.qwgt[1][1] =  ONE             ;
/* weights for three gauss points */
      data.qwgt[0][2] =  0.5555555555556 ;
      data.qwgt[1][2] =  0.8888888888889 ;
      data.qwgt[2][2] =  data.qwgt[0][2]       ;
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
 |     xgr[i][j]                                                        |
 |    wgts[i][j]:  i+1 - actual number of gausspoint                    |
 |                 j+1 - number for integration case (from input)       |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
 |                              CASE 0                                  |
 *----------------------------------------------------------------------*/
      data.txgr[0][0]    =  Q14 ;
      data.txgs[0][0]    =  Q14 ;
      data.txgt[0][0]    =  Q14 ;
      data.twgt[0][0]   =  Q16 ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 2    |
 |                             CASE 1                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][1]    =    pbeta ;
      data.txgr[1][1]    =    palpha;
      data.txgr[2][1]    =    pbeta ;
      data.txgr[3][1]    =    pbeta ;
      data.txgs[0][1]    =    pbeta ;
      data.txgs[1][1]    =    pbeta ;
      data.txgs[2][1]    =    palpha;
      data.txgs[3][1]    =    pbeta ;
      data.txgt[0][1]    =    pbeta ;
      data.txgt[1][1]    =    pbeta ;
      data.txgt[2][1]    =    pbeta ;
      data.txgt[3][1]    =    palpha;
      data.twgt[0][1]   =    Q124  ;
      data.twgt[1][1]   =    Q124  ;
      data.twgt[2][1]   =    Q124  ;
      data.twgt[3][1]   =    Q124  ;
/*----------------------------------------------------------------------*
 |    ALT.GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 1    |
 |                             CASE 2                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][2]    =     ZERO;
      data.txgr[1][2]    =     ONE ;
      data.txgr[2][2]    =     ZERO;
      data.txgr[3][2]    =     ZERO;
      data.txgs[0][2]    =     ZERO;
      data.txgs[1][2]    =     ZERO;
      data.txgs[2][2]    =     ONE ;
      data.txgs[3][2]    =     ZERO;
      data.txgt[0][2]    =     ZERO;
      data.txgt[1][2]    =     ZERO;
      data.txgt[2][2]    =     ZERO;
      data.txgt[3][2]    =     ONE ;
      data.twgt[0][2]   =     Q124;
      data.twgt[1][2]   =     Q124;
      data.twgt[2][2]   =     Q124;
      data.twgt[3][2]   =     Q124;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION        5 SAMPLING POINTS, DEG.OF PRECISION 3    |
 |                             CASE 3                                   |
 *----------------------------------------------------------------------*/
      data.txgr[0][3]    =     Q14  ;
      data.txgr[1][3]    =     Q12  ;
      data.txgr[2][3]    =     Q16  ;
      data.txgr[3][3]    =     Q16  ;
      data.txgr[4][3]    =     Q16  ;
      data.txgs[0][3]    =     Q14  ;
      data.txgs[1][3]    =     Q16  ;
      data.txgs[2][3]    =     Q16  ;
      data.txgs[3][3]    =     Q16  ;
      data.txgs[4][3]    =     Q12  ;
      data.txgt[0][3]    =     Q14  ;
      data.txgt[1][3]    =     Q16  ;
      data.txgt[2][3]    =     Q16  ;
      data.txgt[3][3]    =     Q12  ;
      data.txgt[4][3]    =     Q16  ;
      data.twgt[0][3]   =    -Q430 ;
      data.twgt[1][3]   =     Q9120;
      data.twgt[2][3]   =     Q9120;
      data.twgt[3][3]   =     Q9120;
      data.twgt[4][3]   =     Q9120;

  return;
} //end of DRT::Elements::Fluid3::f3_integration_points




/*----------------------------------------------------------------------*
 |  calculate Jacobian matrix and it's determinant (private) g.bau 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_jaco(const Epetra_SerialDenseMatrix& xyze,
				    const Epetra_SerialDenseMatrix& deriv,
                                    Epetra_SerialDenseMatrix& xjm,
				    double* det,
                                    const int iel
				    )
{
  double dum;

  /*-------------------------------- determine jacobian at point r,s,t---*/
  for (int i=0; i<3; i++)
  {
     for (int j=0; j<3; j++)
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
  *det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
         xjm(0,1)*xjm(1,2)*xjm(2,0)+
         xjm(0,2)*xjm(1,0)*xjm(2,1)-
         xjm(0,2)*xjm(1,1)*xjm(2,0)-
         xjm(0,0)*xjm(1,2)*xjm(2,1)-
         xjm(0,1)*xjm(1,0)*xjm(2,2);

  if(*det<ZERO)
  {
     printf("\n");
     printf("GLOBAL ELEMENT NO.%i\n",Id());
     printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n",*det);
     dserror("Stopped not regulary!\n");
  }

} //end of DRT::Elements::Fluid3::f3_jaco





/*----------------------------------------------------------------------*
 |  shape functions and natural derivatives for hexaeder (private) genk 02/04|
 *----------------------------------------------------------------------*/

void DRT::Elements::Fluid3::f3_shape_function(
               vector<double>&    		funct,
               Epetra_SerialDenseMatrix& 	deriv,
               Epetra_SerialDenseMatrix& 	deriv2,
               const double&      r,
               const double&      s,
               const double&      t,
               const int&         iel,
               int         	  icode
            )
{
/*
In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for
 H E X A H E D E R

   Numbering of the nodes:

                           ^ t
                           |
                           |
                           |
                    8      |  15        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /	|
              16o   |     o       14o   |
               /    o20       o    /	o19
              /     |             /     |
             /      |  13      6 /	|
          5 o---------o---------o	|
            |   o   |     o   	|   o   |  ---------->
            |       o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          17o    /    o         o18  /
            | 12o         o     |   o10
            |  /                |  /
            | /                 | /
            |/	                |/
            o---------o---------o
	    1	/     9         2
	       /
	      /
	     /
	    r

   GiD:

                           ^ t
                           |
                           |
                           |
                    8      |  19        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /	|
              20o   |   26o       18o   |
               /    o16     24o    /	o15
              /     |             /     |
             /      |  17      6 /	|
          5 o---------o---------o   23	|
            |   o   |   27o   	|   o   |  ---------->
            |  25   o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          13o    /  22o         o14  /
            | 12o         o     |   o10
            |  /         21     |  /
            | /                 | /
            |/	                |/
            o---------o---------o
	    1	/     9         2
	       /
	      /
	     /
	    r



   PROBLEM: GID has a different numbering of the element nodes than this one.
            So either the shape functions for hex20 and hex27 (see drawing)
	    has to be adapted or during the input phase the numbering has to
	    be adapted to the shape functions.
	    This is all in progress and should be done for fluid3 and
	    brick1 the same way!!!!

   There are no HEX27 Elements in brick1 so we just go ahead here and
   use the GiD numbering for HEX27.

\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   t        DOUBLE   (i)    coordinate

------------------------------------------------------------------------*/
// variables for use in hex elements
const double Q12 = ONE/TWO;
const double Q14 = ONE/FOUR;
const double Q18 = ONE/EIGHT;
double rp,rm,sp,sm,tp,tm;
double rrm,ssm,ttm;
// variables for use in tet elements
double t1,t2,t3,t4;

/*------------------------------- selection of polynomial interpolation */
switch (iel)
{
case 8: /* LINEAR shape functions and their natural derivatives ----*/

/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   tp=ONE+t;
   tm=ONE-t;

   funct[0]=Q18*rp*sm*tm;
   funct[1]=Q18*rp*sp*tm;
   funct[2]=Q18*rm*sp*tm;
   funct[3]=Q18*rm*sm*tm;
   funct[4]=Q18*rp*sm*tp;
   funct[5]=Q18*rp*sp*tp;
   funct[6]=Q18*rm*sp*tp;
   funct[7]=Q18*rm*sm*tp;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv(0,0)= Q18*sm*tm  ;
      deriv(0,1)= Q18*sp*tm  ;
      deriv(0,2)=-deriv(0,1);
      deriv(0,3)=-deriv(0,0);
      deriv(0,4)= Q18*sm*tp  ;
      deriv(0,5)= Q18*sp*tp  ;
      deriv(0,6)=-deriv(0,5);
      deriv(0,7)=-deriv(0,4);

      deriv(1,0)=-Q18*tm*rp  ;
      deriv(1,1)=-deriv(1,0);
      deriv(1,2)= Q18*tm*rm  ;
      deriv(1,3)=-deriv(1,2);
      deriv(1,4)=-Q18*tp*rp  ;
      deriv(1,5)=-deriv(1,4);
      deriv(1,6)= Q18*tp*rm  ;
      deriv(1,7)=-deriv(1,6);

      deriv(2,0)=-Q18*rp*sm  ;
      deriv(2,1)=-Q18*rp*sp  ;
      deriv(2,2)=-Q18*rm*sp  ;
      deriv(2,3)=-Q18*rm*sm  ;
      deriv(2,4)=-deriv(2,0);
      deriv(2,5)=-deriv(2,1);
      deriv(2,6)=-deriv(2,2);
      deriv(2,7)=-deriv(2,3);
   } /* endif (icode>1) */
   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2(0,0) =  ZERO;
      deriv2(1,0) =  ZERO;
      deriv2(2,0) =  ZERO;
      deriv2(3,0) = -Q18*tm;
      deriv2(4,0) = -Q18*sm;
      deriv2(5,0) =  Q18*rp;

      deriv2(0,1) =  ZERO;
      deriv2(1,1) =  ZERO;
      deriv2(2,1) =  ZERO;
      deriv2(3,1) = -deriv2(3,0);
      deriv2(4,1) = -Q18*sp;
      deriv2(5,1) = -deriv2(5,0);

      deriv2(0,2) =  ZERO;
      deriv2(1,2) =  ZERO;
      deriv2(2,2) =  ZERO;
      deriv2(3,2) =  deriv2(3,0);
      deriv2(4,2) = -deriv2(4,1);
      deriv2(5,2) = -Q18*rm;

      deriv2(0,3) =  ZERO;
      deriv2(1,3) =  ZERO;
      deriv2(2,3) =  ZERO;
      deriv2(3,3) = -deriv2(3,0);
      deriv2(4,3) = -deriv2(4,0);
      deriv2(5,3) = -deriv2(5,2);

      deriv2(0,4) =  ZERO;
      deriv2(1,4) =  ZERO;
      deriv2(2,4) =  ZERO;
      deriv2(3,4) = -Q18*tp;
      deriv2(4,4) = -deriv2(4,0);
      deriv2(5,4) = -deriv2(5,0);

      deriv2(0,5) =  ZERO;
      deriv2(1,5) =  ZERO;
      deriv2(2,5) =  ZERO;
      deriv2(3,5) = -deriv2(3,4);
      deriv2(4,5) = -deriv2(4,1);
      deriv2(5,5) =  deriv2(5,0);

      deriv2(0,6) =  ZERO;
      deriv2(1,6) =  ZERO;
      deriv2(2,6) =  ZERO;
      deriv2(3,6) =  deriv2(3,4);
      deriv2(4,6) =  deriv2(4,1);
      deriv2(5,6) = -deriv2(5,2);

      deriv2(0,7) =  ZERO;
      deriv2(1,7) =  ZERO;
      deriv2(2,7) =  ZERO;
      deriv2(3,7) = -deriv2(3,4);
      deriv2(4,7) =  deriv2(4,0);
      deriv2(5,7) =  deriv2(5,2);
   } /* endif icode==3) */
break;

case 20: /* QUADRATIC shape functions and their natural derivatives
                         without central nodes                      ----*/

   dserror("shape functions for hex20 are not correct!!! Compare with f3f_calfunctderiv.f \n");

/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   tp=ONE+t;
   tm=ONE-t;
   rrm=ONE-r*r;
   ssm=ONE-s*s;
   ttm=ONE-t*t;

   funct[0] =Q18*rp*sm*tm*(rp+sm+tm-FIVE);
   funct[1] =Q18*rp*sp*tm*(rp+sp+tm-FIVE);
   funct[2] =Q18*rm*sp*tm*(rm+sp+tm-FIVE);
   funct[3] =Q18*rm*sm*tm*(rm+sm+tm-FIVE);
   funct[4] =Q18*rp*sm*tp*(rp+sm+tp-FIVE);
   funct[5] =Q18*rp*sp*tp*(rp+sp+tp-FIVE);
   funct[6] =Q18*rm*sp*tp*(rm+sp+tp-FIVE);
   funct[7] =Q18*rm*sm*tp*(rm+sm+tp-FIVE);
   funct[8] =Q14*rp*ssm*tm;
   funct[9] =Q14*rrm*sp*tm;
   funct[10]=Q14*rm*ssm*tm;
   funct[11]=Q14*rrm*sm*tm;
   funct[12]=Q14*rp*ssm*tp;
   funct[13]=Q14*rrm*sp*tp;
   funct[14]=Q14*rm*ssm*tp;
   funct[15]=Q14*rrm*sm*tp;
   funct[16]=Q14*rp*sm*ttm;
   funct[17]=Q14*rp*sp*ttm;
   funct[18]=Q14*rm*sp*ttm;
   funct[19]=Q14*rm*sm*ttm;  /* analytisch gecheckt und fuer OK erklaert!!! */

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv(0,0) = Q18*sm*tm*(TWO*rp+sm+tm-FIVE);
      deriv(1,0) =-Q18*tm*rp*(TWO*sm+tm+rp-FIVE);
      deriv(2,0) =-Q18*rp*sm*(TWO*tm+rp+sm-FIVE);

      deriv(0,1) = Q18*sp*tm*(TWO*rp+sp+tm-FIVE);
      deriv(1,1) = Q18*tm*rp*(TWO*sp+tm+rp-FIVE);
      deriv(2,1) =-Q18*rp*sp*(TWO*tm+rp+sp-FIVE);

      deriv(0,2) =-Q18*sp*tm*(TWO*rm+sp+tm-FIVE);
      deriv(1,2) = Q18*tm*rm*(TWO*sp+tm+rm-FIVE);
      deriv(2,2) =-Q18*rm*sp*(TWO*tm+rm+sp-FIVE);

      deriv(0,3) =-Q18*sm*tm*(TWO*rm+sm+tm-FIVE);
      deriv(1,3) =-Q18*tm*rm*(TWO*sm+tm+rm-FIVE);
      deriv(2,3) =-Q18*rm*sm*(TWO*tm+rm+sm-FIVE);

      deriv(0,4) = Q18*sm*tp*(TWO*rp+sm+tp-FIVE);
      deriv(1,4) =-Q18*tp*rp*(TWO*sm+tp+rp-FIVE);
      deriv(2,4) = Q18*rp*sm*(TWO*tp+rp+sm-FIVE);

      deriv(0,5) = Q18*sp*tp*(TWO*rp+sp+tp-FIVE);
      deriv(1,5) = Q18*tp*rp*(TWO*sp+tp+rp-FIVE);
      deriv(2,5) = Q18*rp*sp*(TWO*tp+rp+sp-FIVE);

      deriv(0,6) =-Q18*sp*tp*(TWO*rm+sp+tp-FIVE);
      deriv(1,6) = Q18*tp*rm*(TWO*sp+tp+rm-FIVE);
      deriv(2,6) = Q18*rm*sp*(TWO*tp+rm+sp-FIVE);

      deriv(0,7) =-Q18*sm*tp*(TWO*rm+sm+tp-FIVE);
      deriv(1,7) =-Q18*tp*rm*(TWO*sm+tp+rm-FIVE);
      deriv(2,7) = Q18*rm*sm*(TWO*tp+rm+sm-FIVE);

      deriv(0,8) = Q14*ssm*tm;
      deriv(1,8) =-Q12*s*tm*rp;
      deriv(2,8) =-Q14*ssm*rp;

      deriv(0,9) =-Q12*r*sp*tm;
      deriv(1,9) = Q14*rrm*tm;
      deriv(2,9) =-Q14*rrm*sp;

      deriv(0,10)=-deriv(0,8);
      deriv(1,10)=-Q12*s*tm*rm;
      deriv(2,10)=-Q14*ssm*rm;

      deriv(0,11)=-Q12*r*sm*tm;
      deriv(1,11)=-deriv(1,9);
      deriv(2,11)=-Q14*rrm*sm;

      deriv(0,12)= Q14*ssm*tp;
      deriv(1,12)=-Q12*s*tp*rp;
      deriv(2,12)=-deriv(2,8);

      deriv(0,13)=-Q12*r*sp*tp;
      deriv(1,13)= Q14*rrm*tp;
      deriv(2,13)=-deriv(2,8);

      deriv(0,14)=-deriv(0,12);
      deriv(1,14)=-Q12*s*tp*rm;
      deriv(2,14)=-deriv(2,10);

      deriv(0,15)=-Q12*r*sm*tp;
      deriv(1,15)=-deriv(1,13);
      deriv(2,15)=-deriv(2,11);

      deriv(0,16)= Q14*sm*ttm;
      deriv(1,16)=-Q14*ttm*rp;
      deriv(2,16)=-Q12*t*rp*sm;

      deriv(0,17)= Q14*sp*ttm;
      deriv(1,17)=-deriv(1,16);
      deriv(2,17)=-Q12*t*rp*sp;

      deriv(0,18)=-deriv(0,17);
      deriv(1,18)= Q14*ttm*rm;
      deriv(2,18)=-Q12*t*rm*sp;

      deriv(0,19)=-deriv(0,16);
      deriv(1,19)=-deriv(1,18);
      deriv(2,19)=-Q12*t*rm*sm;
   } /* endif (icode>1) */
   if(icode==3) /* --> second derivative evaluation  */
   {
      deriv2(0,0) = Q14*sm*tm;
      deriv2(1,0) = Q14*tm*rp;
      deriv2(2,0) = Q14*rp*sm;
      deriv2(3,0) =-Q18*(tm*(2*rp+sm+tm-FIVE+sm*tm));
      deriv2(4,0) =-Q18*(sm*(2*rp+sm+tm-FIVE+sm*tm));
      deriv2(5,0) = Q18*(rp*(2*sm+tm+rp-FIVE+tm*rp));

      deriv2(0,1) = Q14*sp*tm;
      deriv2(1,1) = deriv2(2,1);
      deriv2(2,1) = Q14*rp*sp;
      deriv2(3,1) =-Q18*(tm*(2*rp+sp+tm-FIVE+sp*tm));
      deriv2(4,1) =-Q18*(sp*(2*rp+sp+tm-FIVE+sp*tm));
      deriv2(5,1) =-Q18*(rp*(2*sp+tm+rp-FIVE+tm*rp));

      deriv2(0,2) =-deriv2(1,2);
      deriv2(1,2) = Q14*tm*rm;
      deriv2(2,2) = Q14*rm*sp;
      deriv2(3,2) =-Q18*(tm*(2*rm+sp+tm-FIVE+sp*tm));
      deriv2(4,2) = Q18*(sp*(2*rm+sp+tm-FIVE+sp*tm));
      deriv2(5,2) =-Q18*(rm*(2*sp+tm+rm-FIVE+tm*rm));

      deriv2(0,3) =-deriv2(1,1);
      deriv2(1,3) = deriv2(2,3);
      deriv2(2,3) = Q14*rm*sm;
      deriv2(3,3) =-Q18*(tm*(2*rm+sm+tm-FIVE+sm*tm));
      deriv2(4,3) = Q18*(sm*(2*rm+sm+tm-FIVE+sm*tm));
      deriv2(5,3) = Q18*(rm*(2*sm+tm+rm-FIVE+tm*rm));

      deriv2(0,4) = Q14*sm*tp;
      deriv2(1,4) = Q14*tp*rp;
      deriv2(2,4) = deriv2(3,1);
      deriv2(3,4) =-Q18*(tp*(2*rp+sm+tp-FIVE+sm*tp));
      deriv2(4,4) = Q18*(sm*(2*rp+sm+tp-FIVE+sm*tp));
      deriv2(5,4) =-Q18*(rp*(2*sm+tp+rp-FIVE+tp*rp));

      deriv2(0,5) = Q14*sp*tp;
      deriv2(1,5) = deriv2(2,5);
      deriv2(2,5) = deriv2(3,2);
      deriv2(3,5) =-Q18*(tp*(2*rp+sp+tp-FIVE+sp*tp));
      deriv2(4,5) = Q18*(sp*(2*rp+sp+tp-FIVE+sp*tp));
      deriv2(5,5) = Q18*(rp*(2*sp+tp+rp-FIVE+tp*rp));

      deriv2(0,6) =-deriv2(1,6);
      deriv2(1,6) = Q14*tp*rm;
      deriv2(2,6) = deriv2(3,3);
      deriv2(3,6) =-Q18*(tp*(2*rm+sp+tp-FIVE+sp*tp));
      deriv2(4,6) =-Q18*(sp*(2*rm+sp+tp-FIVE+sp*tp));
      deriv2(5,6) = Q18*(rm*(2*sp+tp+rm-FIVE+tp*rm));

      deriv2(0,7) =-deriv2(1,5);
      deriv2(1,7) = deriv2(2,7);
      deriv2(2,7) = deriv2(3,4);
      deriv2(3,7) =-Q18*(tp*(2*rm+sm+tp-FIVE+sm*tp));
      deriv2(4,7) =-Q18*(sm*(2*rm+sm+tp-FIVE+sm*tp));
      deriv2(5,7) =-Q18*(rm*(2*sm+tp+rm-FIVE+tp*rm));

      deriv2(0,8) = ZERO;
      deriv2(1,8) = -Q12*tm*rp;
      deriv2(2,8) = ZERO;
      deriv2(3,8) =-Q12*s*tm;
      deriv2(4,8) =-Q14*ssm;
      deriv2(5,8) = Q12*s*rp;

      deriv2(0,9)=-Q12*sp*tm;
      deriv2(1,9)= ZERO;
      deriv2(2,9)= ZERO;
      deriv2(3,9)=-Q12*r*tm;
      deriv2(4,9)= Q12*r*sp;
      deriv2(5,9)=-Q14*rrm ;

      deriv2(0,10)= ZERO;
      deriv2(1,10)= -Q12*tm*rm;
      deriv2(2,10)= ZERO;
      deriv2(3,10)= Q12*s*tm;
      deriv2(4,10)=-deriv2(4,8);
      deriv2(5,10)= Q12*s*rm;

      deriv2(0,11)=-Q12*sm*tm;
      deriv2(1,11)= ZERO;
      deriv2(2,11)= ZERO;
      deriv2(3,11)= Q12*r*tm;
      deriv2(4,11)= Q12*r*sm;
      deriv2(5,11)=-deriv2(5,9);

      deriv2(0,12)= ZERO;
      deriv2(1,12)= -Q12*tp*rp;
      deriv2(2,12)= ZERO;
      deriv2(3,12)=-Q12*s*tp;
      deriv2(4,12)=-deriv2(4,8);
      deriv2(5,12)=-deriv2(5,8);

      deriv2(0,13)=-Q12*sp*tp;
      deriv2(1,13)= ZERO;
      deriv2(2,13)= ZERO;
      deriv2(3,13)=-Q12*r*tp;
      deriv2(4,13)=-deriv2(4,9);
      deriv2(5,13)=-deriv2(5,9);

      deriv2(0,14)= ZERO;
      deriv2(1,14)= -Q12*tp*rm;
      deriv2(2,14)= ZERO;
      deriv2(3,14)= Q12*s*tp;
      deriv2(4,14)= deriv2(4,8);
      deriv2(5,14)=-deriv2(5,10);

      deriv2(0,15)=-Q12*sm*tp;
      deriv2(1,15)= ZERO;
      deriv2(2,15)= ZERO;
      deriv2(3,15)= Q12*r*tp;
      deriv2(4,15)=-deriv2(4,11);
      deriv2(5,15)= deriv2(5,9);

      deriv2(0,16)= ZERO;
      deriv2(1,16)= ZERO;
      deriv2(2,16)= ZERO;
      deriv2(3,16)=-Q14*ttm;
      deriv2(4,16)=-Q12*t*sm;
      deriv2(5,16)= Q12*t*rp;

      deriv2(0,17)= ZERO;
      deriv2(1,17)= ZERO;
      deriv2(2,17)= ZERO;
      deriv2(3,17)= Q14*ttm;
      deriv2(4,17)=-Q12*t*sp;
      deriv2(5,17)=-deriv2(5,16);

      deriv2(0,18)= ZERO;
      deriv2(1,18)= ZERO;
      deriv2(2,18)= ZERO;
      deriv2(3,18)= deriv2(3,16);
      deriv2(4,18)= Q12*t*sp;
      deriv2(5,18)= Q12*t*rm;

      deriv2(0,19)= ZERO;
      deriv2(1,19)= ZERO;
      deriv2(2,19)= ZERO;
      deriv2(3,19)= deriv2(3,17);
      deriv2(4,19)= Q12*t*sm;
      deriv2(5,19)=-deriv2(5,18);
   } /* endif (icode==3) */
break;

case 27: /* QUADRATIC shape functions and their natural derivatives
               with central nodes                         ----*/
/*--------------------------------------------------- form basic values */
{
  double drm1,dr00,drp1,dsm1,ds00,dsp1,dtm1,dt00,dtp1;
  double rm1,r00,rp1,sm1,s00,sp1,tm1,t00,tp1;

  rm1=Q12*r*(r - ONE);
  r00=(ONE - r*r);
  rp1=Q12*r*(r + ONE);
  sm1=Q12*s*(s - ONE);
  s00=(ONE - s*s);
  sp1=Q12*s*(s + ONE);
  tm1=Q12*t*(t - ONE);
  t00=(ONE - t*t);
  tp1=Q12*t*(t + ONE);

  drm1 = r - Q12;
  dr00 = -TWO * r;
  drp1 = r + Q12;
  dsm1 = s - Q12;
  ds00 = -TWO * s;
  dsp1 = s + Q12;
  dtm1 = t - Q12;
  dt00 = -TWO * t;
  dtp1 = t + Q12;

  funct[0] = rp1*sp1*tp1;
  funct[1] = sm1*rp1*tp1;
  funct[2] = rm1*sm1*tp1;
  funct[3] = rm1*sp1*tp1;
  funct[4] = tm1*rp1*sp1;
  funct[5] = sm1*tm1*rp1;
  funct[6] = rm1*sm1*tm1;
  funct[7] = rm1*tm1*sp1;
  funct[8] = s00*rp1*tp1;
  funct[9] = r00*sm1*tp1;
  funct[10] = s00*rm1*tp1;
  funct[11] = r00*sp1*tp1;
  funct[12] = t00*rp1*sp1;
  funct[13] = t00*sm1*rp1;
  funct[14] = t00*rm1*sm1;
  funct[15] = t00*rm1*sp1;
  funct[16] = s00*tm1*rp1;
  funct[17] = r00*sm1*tm1;
  funct[18] = s00*rm1*tm1;
  funct[19] = r00*tm1*sp1;
  funct[20] = r00*s00*tp1;
  funct[21] = s00*t00*rp1;
  funct[22] = r00*t00*sm1;
  funct[23] = s00*t00*rm1;
  funct[24] = r00*t00*sp1;
  funct[25] = r00*s00*tm1;
  funct[26] = r00*s00*t00;

  if (icode>1) /* --> first derivative evaluation */
  {
    deriv(0,0) = sp1*tp1*drp1;
    deriv(0,1) = sm1*tp1*drp1;
    deriv(0,2) = sm1*tp1*drm1;
    deriv(0,3) = sp1*tp1*drm1;
    deriv(0,4) = tm1*sp1*drp1;
    deriv(0,5) = sm1*tm1*drp1;
    deriv(0,6) = sm1*tm1*drm1;
    deriv(0,7) = tm1*sp1*drm1;
    deriv(0,8) = s00*tp1*drp1;
    deriv(0,9) = sm1*tp1*dr00;
    deriv(0,10) = s00*tp1*drm1;
    deriv(0,11) = sp1*tp1*dr00;
    deriv(0,12) = t00*sp1*drp1;
    deriv(0,13) = t00*sm1*drp1;
    deriv(0,14) = t00*sm1*drm1;
    deriv(0,15) = t00*sp1*drm1;
    deriv(0,16) = s00*tm1*drp1;
    deriv(0,17) = sm1*tm1*dr00;
    deriv(0,18) = s00*tm1*drm1;
    deriv(0,19) = tm1*sp1*dr00;
    deriv(0,20) = s00*tp1*dr00;
    deriv(0,21) = s00*t00*drp1;
    deriv(0,22) = t00*sm1*dr00;
    deriv(0,23) = s00*t00*drm1;
    deriv(0,24) = t00*sp1*dr00;
    deriv(0,25) = s00*tm1*dr00;
    deriv(0,26) = s00*t00*dr00;

    deriv(1,0) = rp1*tp1*dsp1;
    deriv(1,1) = rp1*tp1*dsm1;
    deriv(1,2) = rm1*tp1*dsm1;
    deriv(1,3) = rm1*tp1*dsp1;
    deriv(1,4) = tm1*rp1*dsp1;
    deriv(1,5) = tm1*rp1*dsm1;
    deriv(1,6) = rm1*tm1*dsm1;
    deriv(1,7) = rm1*tm1*dsp1;
    deriv(1,8) = rp1*tp1*ds00;
    deriv(1,9) = r00*tp1*dsm1;
    deriv(1,10) = rm1*tp1*ds00;
    deriv(1,11) = r00*tp1*dsp1;
    deriv(1,12) = t00*rp1*dsp1;
    deriv(1,13) = t00*rp1*dsm1;
    deriv(1,14) = t00*rm1*dsm1;
    deriv(1,15) = t00*rm1*dsp1;
    deriv(1,16) = tm1*rp1*ds00;
    deriv(1,17) = r00*tm1*dsm1;
    deriv(1,18) = rm1*tm1*ds00;
    deriv(1,19) = r00*tm1*dsp1;
    deriv(1,20) = r00*tp1*ds00;
    deriv(1,21) = t00*rp1*ds00;
    deriv(1,22) = r00*t00*dsm1;
    deriv(1,23) = t00*rm1*ds00;
    deriv(1,24) = r00*t00*dsp1;
    deriv(1,25) = r00*tm1*ds00;
    deriv(1,26) = r00*t00*ds00;

    deriv(2,0) = rp1*sp1*dtp1;
    deriv(2,1) = sm1*rp1*dtp1;
    deriv(2,2) = rm1*sm1*dtp1;
    deriv(2,3) = rm1*sp1*dtp1;
    deriv(2,4) = rp1*sp1*dtm1;
    deriv(2,5) = sm1*rp1*dtm1;
    deriv(2,6) = rm1*sm1*dtm1;
    deriv(2,7) = rm1*sp1*dtm1;
    deriv(2,8) = s00*rp1*dtp1;
    deriv(2,9) = r00*sm1*dtp1;
    deriv(2,10) = s00*rm1*dtp1;
    deriv(2,11) = r00*sp1*dtp1;
    deriv(2,12) = rp1*sp1*dt00;
    deriv(2,13) = sm1*rp1*dt00;
    deriv(2,14) = rm1*sm1*dt00;
    deriv(2,15) = rm1*sp1*dt00;
    deriv(2,16) = s00*rp1*dtm1;
    deriv(2,17) = r00*sm1*dtm1;
    deriv(2,18) = s00*rm1*dtm1;
    deriv(2,19) = r00*sp1*dtm1;
    deriv(2,20) = r00*s00*dtp1;
    deriv(2,21) = s00*rp1*dt00;
    deriv(2,22) = r00*sm1*dt00;
    deriv(2,23) = s00*rm1*dt00;
    deriv(2,24) = r00*sp1*dt00;
    deriv(2,25) = r00*s00*dtm1;
    deriv(2,26) = r00*s00*dt00;
  }
  if (icode==3) /* --> second derivative evaluation */
  {
    deriv2(0,0) = sp1*tp1;
    deriv2(0,1) = sm1*tp1;
    deriv2(0,2) = sm1*tp1;
    deriv2(0,3) = sp1*tp1;
    deriv2(0,4) = tm1*sp1;
    deriv2(0,5) = sm1*tm1;
    deriv2(0,6) = sm1*tm1;
    deriv2(0,7) = tm1*sp1;
    deriv2(0,8) = s00*tp1;
    deriv2(0,9) = -2*sm1*tp1;
    deriv2(0,10) = s00*tp1;
    deriv2(0,11) = -2*sp1*tp1;
    deriv2(0,12) = t00*sp1;
    deriv2(0,13) = t00*sm1;
    deriv2(0,14) = t00*sm1;
    deriv2(0,15) = t00*sp1;
    deriv2(0,16) = s00*tm1;
    deriv2(0,17) = -2*sm1*tm1;
    deriv2(0,18) = s00*tm1;
    deriv2(0,19) = -2*tm1*sp1;
    deriv2(0,20) = -2*s00*tp1;
    deriv2(0,21) = s00*t00;
    deriv2(0,22) = -2*t00*sm1;
    deriv2(0,23) = s00*t00;
    deriv2(0,24) = -2*t00*sp1;
    deriv2(0,25) = -2*s00*tm1;
    deriv2(0,26) = -2*s00*t00;

    deriv2(1,0) = rp1*tp1;
    deriv2(1,1) = rp1*tp1;
    deriv2(1,2) = rm1*tp1;
    deriv2(1,3) = rm1*tp1;
    deriv2(1,4) = tm1*rp1;
    deriv2(1,5) = tm1*rp1;
    deriv2(1,6) = rm1*tm1;
    deriv2(1,7) = rm1*tm1;
    deriv2(1,8) = -2*rp1*tp1;
    deriv2(1,9) = r00*tp1;
    deriv2(1,10) = -2*rm1*tp1;
    deriv2(1,11) = r00*tp1;
    deriv2(1,12) = t00*rp1;
    deriv2(1,13) = t00*rp1;
    deriv2(1,14) = t00*rm1;
    deriv2(1,15) = t00*rm1;
    deriv2(1,16) = -2*tm1*rp1;
    deriv2(1,17) = r00*tm1;
    deriv2(1,18) = -2*rm1*tm1;
    deriv2(1,19) = r00*tm1;
    deriv2(1,20) = -2*r00*tp1;
    deriv2(1,21) = -2*t00*rp1;
    deriv2(1,22) = r00*t00;
    deriv2(1,23) = -2*t00*rm1;
    deriv2(1,24) = r00*t00;
    deriv2(1,25) = -2*r00*tm1;
    deriv2(1,26) = -2*r00*t00;

    deriv2(2,0) = rp1*sp1;
    deriv2(2,1) = sm1*rp1;
    deriv2(2,2) = rm1*sm1;
    deriv2(2,3) = rm1*sp1;
    deriv2(2,4) = rp1*sp1;
    deriv2(2,5) = sm1*rp1;
    deriv2(2,6) = rm1*sm1;
    deriv2(2,7) = rm1*sp1;
    deriv2(2,8) = s00*rp1;
    deriv2(2,9) = r00*sm1;
    deriv2(2,10) = s00*rm1;
    deriv2(2,11) = r00*sp1;
    deriv2(2,12) = -2*rp1*sp1;
    deriv2(2,13) = -2*sm1*rp1;
    deriv2(2,14) = -2*rm1*sm1;
    deriv2(2,15) = -2*rm1*sp1;
    deriv2(2,16) = s00*rp1;
    deriv2(2,17) = r00*sm1;
    deriv2(2,18) = s00*rm1;
    deriv2(2,19) = r00*sp1;
    deriv2(2,20) = r00*s00;
    deriv2(2,21) = -2*s00*rp1;
    deriv2(2,22) = -2*r00*sm1;
    deriv2(2,23) = -2*s00*rm1;
    deriv2(2,24) = -2*r00*sp1;
    deriv2(2,25) = r00*s00;
    deriv2(2,26) = -2*r00*s00;

    deriv2(3,0) = tp1*drp1*dsp1;
    deriv2(3,1) = tp1*dsm1*drp1;
    deriv2(3,2) = tp1*drm1*dsm1;
    deriv2(3,3) = tp1*drm1*dsp1;
    deriv2(3,4) = tm1*drp1*dsp1;
    deriv2(3,5) = tm1*dsm1*drp1;
    deriv2(3,6) = tm1*drm1*dsm1;
    deriv2(3,7) = tm1*drm1*dsp1;
    deriv2(3,8) = tp1*ds00*drp1;
    deriv2(3,9) = tp1*dr00*dsm1;
    deriv2(3,10) = tp1*ds00*drm1;
    deriv2(3,11) = tp1*dr00*dsp1;
    deriv2(3,12) = t00*drp1*dsp1;
    deriv2(3,13) = t00*dsm1*drp1;
    deriv2(3,14) = t00*drm1*dsm1;
    deriv2(3,15) = t00*drm1*dsp1;
    deriv2(3,16) = tm1*ds00*drp1;
    deriv2(3,17) = tm1*dr00*dsm1;
    deriv2(3,18) = tm1*ds00*drm1;
    deriv2(3,19) = tm1*dr00*dsp1;
    deriv2(3,20) = 4*r*s*tp1;
    deriv2(3,21) = t00*ds00*drp1;
    deriv2(3,22) = t00*dr00*dsm1;
    deriv2(3,23) = t00*ds00*drm1;
    deriv2(3,24) = t00*dr00*dsp1;
    deriv2(3,25) = 4*r*s*tm1;
    deriv2(3,26) = 4*r*s*t00;

    deriv2(4,0) = sp1*drp1*dtp1;
    deriv2(4,1) = sm1*drp1*dtp1;
    deriv2(4,2) = sm1*drm1*dtp1;
    deriv2(4,3) = sp1*drm1*dtp1;
    deriv2(4,4) = sp1*dtm1*drp1;
    deriv2(4,5) = sm1*dtm1*drp1;
    deriv2(4,6) = sm1*drm1*dtm1;
    deriv2(4,7) = sp1*drm1*dtm1;
    deriv2(4,8) = s00*drp1*dtp1;
    deriv2(4,9) = sm1*dr00*dtp1;
    deriv2(4,10) = s00*drm1*dtp1;
    deriv2(4,11) = sp1*dr00*dtp1;
    deriv2(4,12) = sp1*dt00*drp1;
    deriv2(4,13) = sm1*dt00*drp1;
    deriv2(4,14) = sm1*dt00*drm1;
    deriv2(4,15) = sp1*dt00*drm1;
    deriv2(4,16) = s00*dtm1*drp1;
    deriv2(4,17) = sm1*dr00*dtm1;
    deriv2(4,18) = s00*drm1*dtm1;
    deriv2(4,19) = sp1*dr00*dtm1;
    deriv2(4,20) = s00*dr00*dtp1;
    deriv2(4,21) = s00*dt00*drp1;
    deriv2(4,22) = 4*r*t*sm1;
    deriv2(4,23) = s00*dt00*drm1;
    deriv2(4,24) = 4*r*t*sp1;
    deriv2(4,25) = s00*dr00*dtm1;
    deriv2(4,26) = 4*r*t*s00;

    deriv2(5,0) = rp1*dsp1*dtp1;
    deriv2(5,1) = rp1*dsm1*dtp1;
    deriv2(5,2) = rm1*dsm1*dtp1;
    deriv2(5,3) = rm1*dsp1*dtp1;
    deriv2(5,4) = rp1*dtm1*dsp1;
    deriv2(5,5) = rp1*dsm1*dtm1;
    deriv2(5,6) = rm1*dsm1*dtm1;
    deriv2(5,7) = rm1*dtm1*dsp1;
    deriv2(5,8) = rp1*ds00*dtp1;
    deriv2(5,9) = r00*dsm1*dtp1;
    deriv2(5,10) = rm1*ds00*dtp1;
    deriv2(5,11) = r00*dsp1*dtp1;
    deriv2(5,12) = rp1*dt00*dsp1;
    deriv2(5,13) = rp1*dt00*dsm1;
    deriv2(5,14) = rm1*dt00*dsm1;
    deriv2(5,15) = rm1*dt00*dsp1;
    deriv2(5,16) = rp1*ds00*dtm1;
    deriv2(5,17) = r00*dsm1*dtm1;
    deriv2(5,18) = rm1*ds00*dtm1;
    deriv2(5,19) = r00*dtm1*dsp1;
    deriv2(5,20) = r00*ds00*dtp1;
    deriv2(5,21) = 4*s*t*rp1;
    deriv2(5,22) = r00*dt00*dsm1;
    deriv2(5,23) = 4*s*t*rm1;
    deriv2(5,24) = r00*dt00*dsp1;
    deriv2(5,25) = r00*ds00*dtm1;
    deriv2(5,26) = 4*s*t*r00;
  }
  break;
}

case 4:

/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for tetraeder

<pre>                                                         genk 08/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for
T E T R A E D E R

</pre>
\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   t        DOUBLE   (i)    coordinate
\return void
\warning shape functions for TET10 not implemented yet!!!
*/

/* LINEAR shape functions and their natural derivatives -----*/
/*--------------------------------------------------- form basic values */

  /*
   Numbering of the nodes:
   -----------------------
   - this is the numbering used in GiD!!


          4 o---
            |\  ---
            |  \   -o3
            |   \  / \
            |     \   \
            |    / \   \
            |   /    \  \
            |  /      \  \
            | /         \ \
            |/            \\
            o---------------o
           1                2
   */


   t1=ONE-r-s-t;
   t2=r;
   t3=s;
   t4=t;

   funct[0]= t1;
   funct[1]= t2;
   funct[2]= t3;
   funct[3]= t4;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv(0,0)=-ONE;
      deriv(0,1)= ONE;
      deriv(0,2)= ZERO;
      deriv(0,3)= ZERO;

      deriv(1,0)=-ONE;
      deriv(1,1)= ZERO;
      deriv(1,2)= ONE;
      deriv(1,3)= ZERO;

      deriv(2,0)=-ONE;
      deriv(2,1)= ZERO;
      deriv(2,2)= ZERO;
      deriv(2,3)= ONE;
   } /* endif (icode>1) */
break;

case 10: /*  QUADRATIC shape functions and their natural derivatives */

   dserror("shape functions for tet10 not implemented yet!\n");

/*--------------------------------------------------- form basic values */
#if 0
   t1=r;
   t2=s;
   t3=t;
   t4=ONE-r-s-t;

   /*These are the shape functions used by Bathe (p. 439) using the corresponding numbering (p.438).*/
   /*If these shape functions go together with the numbering used for the elements, was not checked. -> could be wrong!*/
   funct[0] =1-r-s-t-2*r*(1-r-s-t)-2*s*(1-r-s-t)-2*t*(1-r-s-t);
   funct[1] =r-2*r*(1-r-s-t)-2*r*s-2*r*t;
   funct[2] =(s-2*r*s-2*s*(1-r-s-t)-2*s*t;
   funct[3] =t-2*r*t-2*s*t-2*t*(1-r-s-t);
   funct[4] =4*r*(1-r-s-t);
   funct[5] =4*r*s;
   funct[6] =4*s*(1-r-s-t);
   funct[7] =4*r*t;
   funct[8] =4*s*t;
   funct[9] =4*t*(1-r-s-t);

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv(0,0) = ;
      deriv(1,0) = ;
      deriv(2,0) = ;

      deriv(0,1) = ;
      deriv(1,1) = ;
      deriv(2,1) = ;

      deriv(0,2) = ;
      deriv(1,2) = ;
      deriv(2,2) = ;

      deriv(0,3) = ;
      deriv(1,3) = ;
      deriv(2,3) = ;

      deriv(0,4) = ;
      deriv(1,4) = ;
      deriv(2,4) = ;

      deriv(0,5) = ;
      deriv(1,5) = ;
      deriv(2,5) = ;

      deriv(0,6) = ;
      deriv(1,6) = ;
      deriv(2,6) = ;

      deriv(0,7) = ;
      deriv(1,7) = ;
      deriv(2,7) = ;

      deriv(0,8) = ;
      deriv(1,8) = ;
      deriv(2,8) = ;

      deriv(0,9) = ;
      deriv(1,9) = ;
      deriv(2,9) = ;

   }
   if(icode==3) /* --> second derivative evaluation  */
   {
      deriv2(0,0) =  ;
      deriv2(1,0) =  ;
      deriv2(2,0) =  ;
      deriv2(3,0) = ;
      deriv2(4,0) = ;
      deriv2(5,0) = ;

      deriv2(0,1) =  ;
      deriv2(1,1) =  ;
      deriv2(2,1) =  ;
      deriv2(3,1) = ;
      deriv2(4,1) = ;
      deriv2(5,1) = ;

      deriv2(0,2) =  ;
      deriv2(1,2) = ;
      deriv2(2,2) =  ;
      deriv2(3,2) = ;
      deriv2(4,2) = ;
      deriv2(5,2) = ;

      deriv2(0,3) = ;
      deriv2(1,3) =  ;
      deriv2(2,3) =  ;
      deriv2(3,3) = ;
      deriv2(4,3) = ;
      deriv2(5,3) = ;

      deriv2(0,4) =  ;
      deriv2(1,4) =  ;
      deriv2(2,4) =  ;
      deriv2(3,4) = ;
      deriv2(4,4) = ;
      deriv2(5,4) = ;

      deriv2(0,5) =  ;
      deriv2(1,5) =  ;
      deriv2(2,5) =  ;
      deriv2(3,5) = ;
      deriv2(4,5) = ;
      deriv2(5,5) = ;

      deriv2(0,6) = ;
      deriv2(1,6) =  ;
      deriv2(2,6) =  ;
      deriv2(3,6) = ;
      deriv2(4,6) = ;
      deriv2(5,6) = ;

      deriv2(0,7) = ;
      deriv2(1,7) = ;
      deriv2(2,7) = ;
      deriv2(3,7) = ;
      deriv2(4,7) = ;
      deriv2(5,7) = ;

      deriv2(0,8) =  ;
      deriv2(1,8) =  ;
      deriv2(2,8) =  ;
      deriv2(3,8) = ;
      deriv2(4,8) = ;
      deriv2(5,8) =  ;

      deriv2(0,9)= ;
      deriv2(1,9)=  ;
      deriv2(2,9)=  ;
      deriv2(3,9)= ;
      deriv2(4,9)= ;
      deriv2(5,9)=  ;
   }
#endif
break;

/*----------------------------------------------------------------------*/
default:
   dserror("distyp unknown\n");
} /* end switch(iel) */

return;
} // end of DRT:Elements:Fluid3:f3_shape_function


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_getbodyforce
  (Epetra_SerialDenseMatrix& edeadng,
   double                    time   ,
   const int                 iel    ,
   ParameterList& 	     params
)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  int count=0;
  for(int nn=0;nn<iel;nn++)
  {
    Nodes()[nn]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      count++;
    }
  }

  if (count == iel)
  {

    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum>=0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0)
      {
        curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
        curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for(int nn=0;nn<iel;nn++)
    {
      Nodes()[nn]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int dim=0;dim<3;dim++)
      {
        edeadng(dim,nn)=(*onoff)[dim]*(*val)[dim]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    for(int nn=0;nn<iel;nn++)
    {
      for(int dim=0;dim<3;dim++)
      {
        edeadng(dim,nn)=0.0;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  calculate global derivatives w.r.t. x,y,z at point r,s,t (private)genk05/02
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_gder(Epetra_SerialDenseMatrix& derxy,
				    const Epetra_SerialDenseMatrix& deriv,
                                    Epetra_SerialDenseMatrix& xjm,
				    double& det,
                                    const int iel
				    )
{
  Epetra_SerialDenseMatrix 	xji(3,3);   // inverse of jacobian matrix


 /*----------calculate global derivatives w.r.t. x,y,z at point r,s,t ---*/

	/*------------------------------------------------------- initialistion */
for(int k=0;k<iel;k++)
{
   derxy(0,k)=ZERO;
   derxy(1,k)=ZERO;
   derxy(2,k)=ZERO;
} /* end of loop over k */

	/*------------------------------------------------- inverse of jacobian */
xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

	/*---------------------------------------- calculate global derivatives */
for (int k=0;k<iel;k++)
{
   derxy(0,k) +=   xji(0,0) * deriv(0,k) \
                  + xji(0,1) * deriv(1,k) \
                  + xji(0,2) * deriv(2,k) ;
   derxy(1,k) +=   xji(1,0) * deriv(0,k) \
                  + xji(1,1) * deriv(1,k) \
                  + xji(1,2) * deriv(2,k) ;
   derxy(2,k) +=   xji(2,0) * deriv(0,k) \
                  + xji(2,1) * deriv(1,k) \
                  + xji(2,2) * deriv(2,k) ;
} /* end of loop over k */

	/*----------------------------------------------------------------------*/

return;
} // end of DRT:Elements:Fluid3:f3_gder


void DRT::Elements::Fluid3::f3_gder2(const Epetra_SerialDenseMatrix& xyze,
				    const Epetra_SerialDenseMatrix& xjm,
                                    const Epetra_SerialDenseMatrix& derxy,
				    Epetra_SerialDenseMatrix& derxy2,
				    const Epetra_SerialDenseMatrix& deriv2,
                                    const int iel
				    )
{
double r0,r1,r2,r3,r4,r5;
//--------------------------------------------initialize and zero out everything
Epetra_SerialDenseMatrix bm(6,6);
Epetra_SerialDenseMatrix xder2(6,3);

/*--------------------------- calculate elements of jacobian_bar matrix */
bm(0,0) = xjm(0,0)*xjm(0,0);
bm(1,0) = xjm(1,0)*xjm(1,0);
bm(2,0) = xjm(2,0)*xjm(2,0);
bm(3,0) = xjm(0,0)*xjm(1,0);
bm(4,0) = xjm(0,0)*xjm(2,0);
bm(5,0) = xjm(1,0)*xjm(2,0);

bm(0,1) = xjm(0,1)*xjm(0,1);
bm(1,1) = xjm(1,1)*xjm(1,1);
bm(2,1) = xjm(2,1)*xjm(2,1);
bm(3,1) = xjm(0,1)*xjm(1,1);
bm(4,1) = xjm(0,1)*xjm(2,1);
bm(5,1) = xjm(1,1)*xjm(2,1);

bm(0,2) = xjm(0,2)*xjm(0,2);
bm(1,2) = xjm(1,2)*xjm(1,2);
bm(2,2) = xjm(2,2)*xjm(2,2);
bm(3,2) = xjm(0,2)*xjm(1,2);
bm(4,2) = xjm(0,2)*xjm(2,2);
bm(5,2) = xjm(1,2)*xjm(2,2);

bm(0,3) = TWO*xjm(0,0)*xjm(0,1);
bm(1,3) = TWO*xjm(1,0)*xjm(1,1);
bm(2,3) = TWO*xjm(2,0)*xjm(2,1);
bm(3,3) = xjm(0,0)*xjm(1,1)+xjm(1,0)*xjm(0,1);
bm(4,3) = xjm(0,0)*xjm(2,1)+xjm(2,0)*xjm(0,1);
bm(5,3) = xjm(1,0)*xjm(2,1)+xjm(2,0)*xjm(1,1);

bm(0,4) = TWO*xjm(0,0)*xjm(0,2);
bm(1,4) = TWO*xjm(1,0)*xjm(1,2);
bm(2,4) = TWO*xjm(2,0)*xjm(2,2);
bm(3,4) = xjm(0,0)*xjm(1,2)+xjm(1,0)*xjm(0,2);
bm(4,4) = xjm(0,0)*xjm(2,2)+xjm(2,0)*xjm(0,2);
bm(5,4) = xjm(1,0)*xjm(2,2)+xjm(2,0)*xjm(1,2);

bm(0,5) = TWO*xjm(0,1)*xjm(0,2);
bm(1,5) = TWO*xjm(1,1)*xjm(1,2);
bm(2,5) = TWO*xjm(2,1)*xjm(2,2);
bm(3,5) = xjm(0,1)*xjm(1,2)+xjm(1,1)*xjm(0,2);
bm(4,5) = xjm(0,1)*xjm(2,2)+xjm(2,1)*xjm(0,2);
bm(5,5) = xjm(1,1)*xjm(2,2)+xjm(2,1)*xjm(1,2);

/*-------------------------------------- inverse of jacobian_bar matrix */

LINALG::NonSymmetricInverse(bm,6);

// output for debug
// (for more details see the comments in definition of NonSymmetricInverse()
#if 0
for (int i = 0 ; i < 6; ++i)
{
for (int j = 0 ; j < 6; ++j)
{
       if (bm(i,j)!=0.0)
       printf("bm[%d][%d] %22.16e ",i,j,bm(i,j));
       else
        printf("bm[%d][%d] 0.000 ",i,j);
       printf("\n");
}
printf("\n");
}
#endif

/*----------------------------------------------------------- initialise*/
/*   already initialized by constructor of EpetraSerialDenseMeatrix
for (int i=0;i<3;i++)
{
   for (int j=0;j<6;j++) xder2(j,i)=ZERO;
}
*/

for (int i=0;i<iel;i++)
{
   for (int j=0;j<6;j++) derxy2(j,i)=ZERO;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (int i=0;i<iel;i++)
{
   xder2(0,0) += deriv2(0,i) * xyze(0,i);
   xder2(1,0) += deriv2(1,i) * xyze(0,i);
   xder2(2,0) += deriv2(2,i) * xyze(0,i);
   xder2(3,0) += deriv2(3,i) * xyze(0,i);
   xder2(4,0) += deriv2(4,i) * xyze(0,i);
   xder2(5,0) += deriv2(5,i) * xyze(0,i);

   xder2(0,1) += deriv2(0,i) * xyze(1,i);
   xder2(1,1) += deriv2(1,i) * xyze(1,i);
   xder2(2,1) += deriv2(2,i) * xyze(1,i);
   xder2(3,1) += deriv2(3,i) * xyze(1,i);
   xder2(4,1) += deriv2(4,i) * xyze(1,i);
   xder2(5,1) += deriv2(5,i) * xyze(1,i);

   xder2(0,2) += deriv2(0,i) * xyze(2,i);
   xder2(1,2) += deriv2(1,i) * xyze(2,i);
   xder2(2,2) += deriv2(2,i) * xyze(2,i);
   xder2(3,2) += deriv2(3,i) * xyze(2,i);
   xder2(4,2) += deriv2(4,i) * xyze(2,i);
   xder2(5,2) += deriv2(5,i) * xyze(2,i);
} /* end of loop over i */

/*--------------------------------- calculate second global derivatives */
for (int i=0;i<iel;i++)
{
   r0 = deriv2(0,i) - xder2(0,0)*derxy(0,i) - xder2(0,1)*derxy(1,i) \
                     - xder2(0,2)*derxy(2,i);
   r1 = deriv2(1,i) - xder2(1,0)*derxy(0,i) - xder2(1,1)*derxy(1,i) \
                     - xder2(1,2)*derxy(2,i);
   r2 = deriv2(2,i) - xder2(2,0)*derxy(0,i) - xder2(2,1)*derxy(1,i) \
                     - xder2(2,2)*derxy(2,i);
   r3 = deriv2(3,i) - xder2(3,0)*derxy(0,i) - xder2(3,1)*derxy(1,i) \
                     - xder2(3,2)*derxy(2,i);
   r4 = deriv2(4,i) - xder2(4,0)*derxy(0,i) - xder2(4,1)*derxy(1,i) \
                     - xder2(4,2)*derxy(2,i);
   r5 = deriv2(5,i) - xder2(5,0)*derxy(0,i) - xder2(5,1)*derxy(1,i) \
                     - xder2(5,2)*derxy(2,i);

   derxy2(0,i) += bm(0,0)*r0 + bm(0,1)*r1 + bm(0,2)*r2 \
                +  bm(0,3)*r3 + bm(0,4)*r4 + bm(0,5)*r5;
   derxy2(1,i) += bm(1,0)*r0 + bm(1,1)*r1 + bm(1,2)*r2 \
                +  bm(1,3)*r3 + bm(1,4)*r4 + bm(1,5)*r5;
   derxy2(2,i) += bm(2,0)*r0 + bm(2,1)*r1 + bm(2,2)*r2 \
                +  bm(2,3)*r3 + bm(2,4)*r4 + bm(2,5)*r5;
   derxy2(3,i) += bm(3,0)*r0 + bm(3,1)*r1 + bm(3,2)*r2 \
                +  bm(3,3)*r3 + bm(3,4)*r4 + bm(3,5)*r5;
   derxy2(4,i) += bm(4,0)*r0 + bm(4,1)*r1 + bm(4,2)*r2 \
                +  bm(4,3)*r3 + bm(4,4)*r4 + bm(4,5)*r5;
   derxy2(5,i) += bm(5,0)*r0 + bm(5,1)*r1 + bm(5,2)*r2 \
                +  bm(5,3)*r3 + bm(5,4)*r4 + bm(5,5)*r5;
} /* end of loop over i */

/*----------------------------------------------------------------------*/

return;
} // end of DRT:Elements:Fluid3:f3_gder2



/*----------------------------------------------------------------------*
 |  evaluate fluid coefficient matrix (private)              chfoe 04/04|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilised fluid3 element are calculated. The procedure is
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

void DRT::Elements::Fluid3::f3_calmat( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  eforce,
                vector<double>&            velint,
                vector<double>&            histvec,
                vector<double>&            gridvint,
		double&   	           press,
                Epetra_SerialDenseMatrix&  vderxy,
                Epetra_SerialDenseMatrix&  vderxy2,
                vector<double>&            gradp,
                vector<double>&            funct,
                vector<double>&            tau,
                Epetra_SerialDenseMatrix&  derxy,
                Epetra_SerialDenseMatrix&  derxy2,
                vector<double>&            edeadng,
                double&                    fac,
                const double&              visc,
                const int&                 iel,
		ParameterList& 	           params
                )
{
//DOUBLE  viscous[3][3][3*iel];	/* viscous term partially integrated */

/*========================= further variables =========================*/

Epetra_SerialDenseMatrix  viscs2(3,3*iel);   	/* viscous term incluiding 2nd derivatives */
vector<double>  conv_c(iel); 		/* linearisation of convect, convective part */
vector<double>  conv_g(iel);       	/* linearisation of convect, grid part */
Epetra_SerialDenseMatrix  conv_r(3,3*iel);	/* linearisation of convect, reactive part */
vector<double>  div(3*iel);          	/* divergence of u or v              */
Epetra_SerialDenseMatrix  ugradv(iel,3*iel);	/* linearisation of u * grad v   */
vector<double>  conv_old(3); 		/* convective term evalaluated with old velocities */
vector<double>  conv_g_old(3);
vector<double>  visc_old(3); 		/* viscous term evaluated with old velocities      */
vector<double>  rhsint(3);   		/* total right hand side terms at int.-point       */
Epetra_SerialDenseMatrix  vconv_r(3,iel);

/*========================== initialisation ============================*/
// One-step-Theta: timefac = theta*dt
// BDF2:           timefac = 2/3 * dt
double timefac = params.get<double>("time constant for integration",-1.0);
  if (timefac < 0.0) dserror("No time constant for integration supplied");

// time step size
//double dt = params.get<double>("delta time",-1.0);
//  if (dt == -1.0) dserror("No dta supplied");

// stabilisation parameter
double tau_M  = tau[0]*fac;
double tau_Mp = tau[1]*fac;
double tau_C  = tau[2]*fac;

// integration factors and coefficients of single terms
// double time2nue   = timefac * 2.0 * visc;
double timetauM   = timefac * tau_M;
double timetauMp  = timefac * tau_Mp;

double ttimetauM  = timefac * timetauM;
double ttimetauMp = timefac * timetauMp;
double timefacfac = timefac * fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
    rhsint[0] = histvec[0] + edeadng[0]*timefac;
    rhsint[1] = histvec[1] + edeadng[1]*timefac;
    rhsint[2] = histvec[2] + edeadng[2]*timefac;
/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1]
            + vderxy(0,2) * velint[2];
conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1]
            + vderxy(1,2) * velint[2];
conv_old[2] = vderxy(2,0) * velint[0] + vderxy(2,1) * velint[1]
            + vderxy(2,2) * velint[2];

/* new for incremental formulation: */
/* Convective term  u_G_old * grad u_old: */
conv_g_old[0] = (vderxy(0,0) * gridvint[0] +
		 vderxy(0,1) * gridvint[1] +
		 vderxy(0,2) * gridvint[2]);
conv_g_old[1] = (vderxy(1,0) * gridvint[0] +
		 vderxy(1,1) * gridvint[1] +
		 vderxy(1,2) * gridvint[2]);
conv_g_old[2] = (vderxy(2,0) * gridvint[0] +
		 vderxy(2,1) * gridvint[1] +
		 vderxy(2,2) * gridvint[2]);

/* Viscous term  div epsilon(u_old) */
visc_old[0] = vderxy2(0,0) + 0.5 * ( vderxy2(0,1) + vderxy2(1,3)
                                    + vderxy2(0,2) + vderxy2(2,4));
visc_old[1] = vderxy2(1,1) + 0.5 * ( vderxy2(1,0) + vderxy2(0,3)
                                    + vderxy2(1,2) + vderxy2(2,5));
visc_old[2] = vderxy2(2,2) + 0.5 * ( vderxy2(2,0) + vderxy2(0,4)
                                    + vderxy2(2,1) + vderxy2(1,5));

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* Reactive term  u:  funct */
   /* linearise convective term */

   /*--- convective part u_old * grad (funct) --------------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      with  N .. form function matrix                                   */
   conv_c[i] = derxy(0,i) * velint[0] + derxy(1,i) * velint[1]
             + derxy(2,i) * velint[2];

   /*--- convective grid part u_G * grad (funct) -----------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   if (is_ale_)
   {
     conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1]
                 - derxy(2,i) * gridvint[2];
   }
   else
   {
     conv_g[i] = 0.0;
   }

   /*--- reactive part funct * grad (u_old) ----------------------------*/
   /* /                                     \
      |  u_old_x,x   u_old_x,y   u_old x,z  |
      |                                     |
      |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
      |                                     |
      |  u_old_z,x   u_old_z,y   u_old_z,z  |
      \                                     /
      with  N .. form function matrix                                   */

   conv_r(0,3*i)   = vderxy(0,0)*funct[i];
   conv_r(0,3*i+1) = vderxy(0,1)*funct[i];
   conv_r(0,3*i+2) = vderxy(0,2)*funct[i];
   conv_r(1,3*i)   = vderxy(1,0)*funct[i];
   conv_r(1,3*i+1) = vderxy(1,1)*funct[i];
   conv_r(1,3*i+2) = vderxy(1,2)*funct[i];
   conv_r(2,3*i)   = vderxy(2,0)*funct[i];
   conv_r(2,3*i+1) = vderxy(2,1)*funct[i];
   conv_r(2,3*i+2) = vderxy(2,2)*funct[i];

   vconv_r(0,i) = conv_r(0,3*i)*velint[0] + conv_r(0,3*i+1)*velint[1] + conv_r(0,3*i+2)*velint[2];
   vconv_r(1,i) = conv_r(1,3*i)*velint[0] + conv_r(1,3*i+1)*velint[1] + conv_r(1,3*i+2)*velint[2];
   vconv_r(2,i) = conv_r(2,3*i)*velint[0] + conv_r(2,3*i+1)*velint[1] + conv_r(2,3*i+2)*velint[2];

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

   viscs2(0,3*i)   = - 0.5 * (2.0 * derxy2(0,i) + derxy2(1,i) + derxy2(2,i));
   viscs2(0,3*i+1) = - 0.5 *  derxy2(3,i);
   viscs2(0,3*i+2) = - 0.5 *  derxy2(4,i);
   viscs2(1,3*i)   = - 0.5 *  derxy2(3,i);
   viscs2(1,3*i+1) = - 0.5 * (derxy2(0,i) + 2.0 * derxy2(1,i) + derxy2(2,i));
   viscs2(1,3*i+2) = - 0.5 *  derxy2(5,i);
   viscs2(2,3*i)   = - 0.5 *  derxy2(4,i);
   viscs2(2,3*i+1) = - 0.5 *  derxy2(5,i);
   viscs2(2,3*i+2) = - 0.5 * (derxy2(0,i) + derxy2(1,i) + 2.0 * derxy2(2,i));

   /*--- viscous term (after integr. by parts) -------------------------*/
   /*   /                                             \
        |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
      1 |                                             |
      - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
      2 |                                             |
        |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
        \                                             /
   with N_x .. x-line of N
        N_y .. y-line of N
        N_z .. z-line of N                                              */


/* not needed for incremental solver routine     g.bau 03/07
   viscous[0][0][3*i]   = derxy(0,i);
   viscous[0][0][3*i+1] = 0.0;
   viscous[0][0][3*i+2] = 0.0;                // 1st index:
   viscous[0][1][3*i]   = 0.5 * derxy(1,i);  //   line of epsilon
   viscous[0][1][3*i+1] = 0.5 * derxy(0,i);  // 2nd index:
   viscous[0][1][3*i+2] = 0.0;                //   column of epsilon
   viscous[0][2][3*i]   = 0.5 * derxy(2,i);  // 3rd index:
   viscous[0][2][3*i+1] = 0.0;                //   elemental vel dof
   viscous[0][2][3*i+2] = 0.5 * derxy(0,i);
   viscous[1][0][3*i]   = 0.5 * derxy(1,i);
   viscous[1][0][3*i+1] = 0.5 * derxy(0,i);
   viscous[1][0][3*i+2] = 0.0;
   viscous[1][1][3*i]   = 0.0;
   viscous[1][1][3*i+1] = derxy(1,i);
   viscous[1][1][3*i+2] = 0.0;
   viscous[1][2][3*i]   = 0.0;
   viscous[1][2][3*i+1] = 0.5 * derxy(2,i);
   viscous[1][2][3*i+2] = 0.5 * derxy(1,i);
   viscous[2][0][3*i]   = 0.5 * derxy(2,i);
   viscous[2][0][3*i+1] = 0.0;
   viscous[2][0][3*i+2] = 0.5 * derxy(0,i);
   viscous[2][1][3*i]   = 0.0;
   viscous[2][1][3*i+1] = 0.5 * derxy(2,i);
   viscous[2][1][3*i+2] = 0.5 * derxy(1,i);
   viscous[2][2][3*i]   = 0.0;
   viscous[2][2][3*i+1] = 0.0;
   viscous[2][2][3*i+2] = derxy(2,i);
*/

   /* pressure gradient term derxy, funct without or with integration   *
    * by parts, respectively                                            */

   /*--- divergence u term ---------------------------------------------*/
   div[3*i]   = derxy(0,i);
   div[3*i+1] = derxy(1,i);
   div[3*i+2] = derxy(2,i);

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
      ugradv(i,3*j)   = derxy(0,i) * funct[j];
      ugradv(i,3*j+1) = derxy(1,i) * funct[j];
      ugradv(i,3*j+2) = derxy(2,i) * funct[j];
   }

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define vderxyz_(i,j)  vderxy(i,j)
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r(i,3*(k)+j)
#define vconv_r_(i,j)  vconv_r(i,j)
#define conv_old_(j)   conv_old[j]
#define conv_g_old_(j) conv_g_old[j]
#define derxyz_(i,j)   derxy(i,j)
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2(i,3*(k)+j)
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
#include "fluid3_stiff_ale.cpp"
#include "fluid3_rhs_incr_ale.cpp"
  }
  else
  {
#include "fluid3_stiff.cpp"
#include "fluid3_rhs_incr.cpp"
  }

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef vconv_r_
#undef conv_old_
#undef conv_g_old_
#undef derxyz_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef funct_
#undef vderxyz_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

return;
} // end of DRT:Elements:Fluid3:f3_calmat


void DRT::Elements::Fluid3::f3_calmat_stationary( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  eforce,
                vector<double>&            velint,
                vector<double>&            histvec,
                vector<double>&            gridvint,
		double&   	           press,
                Epetra_SerialDenseMatrix&  vderxy,
                Epetra_SerialDenseMatrix&  vderxy2,
                vector<double>&            gradp,
                vector<double>&            funct,
                vector<double>&            tau,
                Epetra_SerialDenseMatrix&  derxy,
                Epetra_SerialDenseMatrix&  derxy2,
                vector<double>&            edeadng,
                double&                    fac,
                const double&              visc,
                const int&                 iel,
		ParameterList& 	           params
                )
{

/*========================= further variables =========================*/

Epetra_SerialDenseMatrix  viscs2(3,3*iel);   	/* viscous term incluiding 2nd derivatives */
vector<double>  conv_c(iel); 		/* linearisation of convect, convective part */
vector<double>  conv_g(iel);       	/* linearisation of convect, grid part */
Epetra_SerialDenseMatrix  conv_r(3,3*iel);	/* linearisation of convect, reactive part */
vector<double>  div(3*iel);          	/* divergence of u or v              */
Epetra_SerialDenseMatrix  ugradv(iel,3*iel);	/* linearisation of u * grad v   */
vector<double>  conv_old(3); 		/* convective term evalaluated with old velocities */
vector<double>  conv_g_old(3);
vector<double>  visc_old(3); 		/* viscous term evaluated with old velocities      */
vector<double>  rhsint(3);   		/* total right hand side terms at int.-point       */
Epetra_SerialDenseMatrix  vconv_r(3,iel);

/*========================== initialisation ============================*/
// One-step-Theta: timefac = theta*dt
// BDF2:           timefac = 2/3 * dt
  double timefac = params.get<double>("time constant for integration",-1.0);
  if (timefac < 0.0) dserror("No time constant for integration supplied");

// stabilisation parameter
double tau_M  = tau[0]*fac;
double tau_Mp = tau[1]*fac;
double tau_C  = tau[2]*fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
    rhsint[0] = histvec[0] + edeadng[0]*timefac;
    rhsint[1] = histvec[1] + edeadng[1]*timefac;
    rhsint[2] = histvec[2] + edeadng[2]*timefac;
/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1]
            + vderxy(0,2) * velint[2];
conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1]
            + vderxy(1,2) * velint[2];
conv_old[2] = vderxy(2,0) * velint[0] + vderxy(2,1) * velint[1]
            + vderxy(2,2) * velint[2];

/* new for incremental formulation: */
/* Convective term  u_G_old * grad u_old: */
conv_g_old[0] = (vderxy(0,0) * gridvint[0] +
		 vderxy(0,1) * gridvint[1] +
		 vderxy(0,2) * gridvint[2]);
conv_g_old[1] = (vderxy(1,0) * gridvint[0] +
		 vderxy(1,1) * gridvint[1] +
		 vderxy(1,2) * gridvint[2]);
conv_g_old[2] = (vderxy(2,0) * gridvint[0] +
		 vderxy(2,1) * gridvint[1] +
		 vderxy(2,2) * gridvint[2]);

/* Viscous term  div epsilon(u_old) */
visc_old[0] = vderxy2(0,0) + 0.5 * ( vderxy2(0,1) + vderxy2(1,3)
                                    + vderxy2(0,2) + vderxy2(2,4));
visc_old[1] = vderxy2(1,1) + 0.5 * ( vderxy2(1,0) + vderxy2(0,3)
                                    + vderxy2(1,2) + vderxy2(2,5));
visc_old[2] = vderxy2(2,2) + 0.5 * ( vderxy2(2,0) + vderxy2(0,4)
                                    + vderxy2(2,1) + vderxy2(1,5));

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* Reactive term  u:  funct */
   /* linearise convective term */

   /*--- convective part u_old * grad (funct) --------------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      with  N .. form function matrix                                   */
   conv_c[i] = derxy(0,i) * velint[0] + derxy(1,i) * velint[1]
             + derxy(2,i) * velint[2];

   /*--- convective grid part u_G * grad (funct) -----------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   if(is_ale_)
   {
     dserror("No ALE supported by Fluid3 at the moment.");
      //    conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1]
      //           - derxy(2,i) * gridvint[2];
   }
   else
   {
     conv_g[i] = 0.0;
   }

   /*--- reactive part funct * grad (u_old) ----------------------------*/
   /* /                                     \
      |  u_old_x,x   u_old_x,y   u_old x,z  |
      |                                     |
      |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
      |                                     |
      |  u_old_z,x   u_old_z,y   u_old_z,z  |
      \                                     /
      with  N .. form function matrix                                   */

   conv_r(0,3*i)   = vderxy(0,0)*funct[i];
   conv_r(0,3*i+1) = vderxy(0,1)*funct[i];
   conv_r(0,3*i+2) = vderxy(0,2)*funct[i];
   conv_r(1,3*i)   = vderxy(1,0)*funct[i];
   conv_r(1,3*i+1) = vderxy(1,1)*funct[i];
   conv_r(1,3*i+2) = vderxy(1,2)*funct[i];
   conv_r(2,3*i)   = vderxy(2,0)*funct[i];
   conv_r(2,3*i+1) = vderxy(2,1)*funct[i];
   conv_r(2,3*i+2) = vderxy(2,2)*funct[i];

   vconv_r(0,i) = conv_r(0,3*i)*velint[0] + conv_r(0,3*i+1)*velint[1] + conv_r(0,3*i+2)*velint[2];
   vconv_r(1,i) = conv_r(1,3*i)*velint[0] + conv_r(1,3*i+1)*velint[1] + conv_r(1,3*i+2)*velint[2];
   vconv_r(2,i) = conv_r(2,3*i)*velint[0] + conv_r(2,3*i+1)*velint[1] + conv_r(2,3*i+2)*velint[2];

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

   viscs2(0,3*i)   = - 0.5 * (2.0 * derxy2(0,i) + derxy2(1,i) + derxy2(2,i));
   viscs2(0,3*i+1) = - 0.5 *  derxy2(3,i);
   viscs2(0,3*i+2) = - 0.5 *  derxy2(4,i);
   viscs2(1,3*i)   = - 0.5 *  derxy2(3,i);
   viscs2(1,3*i+1) = - 0.5 * (derxy2(0,i) + 2.0 * derxy2(1,i) + derxy2(2,i));
   viscs2(1,3*i+2) = - 0.5 *  derxy2(5,i);
   viscs2(2,3*i)   = - 0.5 *  derxy2(4,i);
   viscs2(2,3*i+1) = - 0.5 *  derxy2(5,i);
   viscs2(2,3*i+2) = - 0.5 * (derxy2(0,i) + derxy2(1,i) + 2.0 * derxy2(2,i));

   /*--- viscous term (after integr. by parts) -------------------------*/
   /*   /                                             \
        |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
      1 |                                             |
      - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
      2 |                                             |
        |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
        \                                             /
   with N_x .. x-line of N
        N_y .. y-line of N
        N_z .. z-line of N                                              */


   /* pressure gradient term derxy, funct without or with integration   *
    * by parts, respectively                                            */

   /*--- divergence u term ---------------------------------------------*/
   div[3*i]   = derxy(0,i);
   div[3*i+1] = derxy(1,i);
   div[3*i+2] = derxy(2,i);

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
      ugradv(i,3*j)   = derxy(0,i) * funct[j];
      ugradv(i,3*j+1) = derxy(1,i) * funct[j];
      ugradv(i,3*j+2) = derxy(2,i) * funct[j];
   }

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define vderxyz_(i,j)  vderxy(i,j)
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r(i,3*(k)+j)
#define vconv_r_(i,j)  vconv_r(i,j)
#define conv_old_(j)   conv_old[j]
#define conv_g_old_(j) conv_g_old[j]
#define derxyz_(i,j)   derxy(i,j)
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2(i,3*(k)+j)
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

  {
    #include "fluid3_stiff_stationary.cpp"
    #include "fluid3_rhs_incr_stationary.cpp"
  }

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef vconv_r_
#undef conv_old_
#undef conv_g_old_
#undef derxyz_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef funct_
#undef vderxyz_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

return;
} // end of DRT:Elements:Fluid3:f3_calmat_stationary


/*----------------------------------------------------------------------*
 |  evaluate fluid coefficient matrix for generalised alpha             |
 |                            (private)                      gammi 06/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Fluid3::f3_genalpha_calmat(
    Epetra_SerialDenseMatrix&  elemat,
    vector<double>&  	       accintam,
    vector<double>&  	       velintaf,
    Epetra_SerialDenseMatrix&  vderxyaf,
    Epetra_SerialDenseMatrix&  vderxyaf2,
    vector<double>&  	       velintnp,
    Epetra_SerialDenseMatrix&  vderxynp,
    double&                    prenp,
    vector<double>&  	       pderxynp,
    vector<double>&  	       edeadaf,
    vector<double>&  	       funct,
    Epetra_SerialDenseMatrix&  derxy,
    Epetra_SerialDenseMatrix&  derxy2,
    vector<double>&  	       tau,
    double&                    fac,
    const double&              visc,
    const int&                 iel,
    ParameterList& 	       params)
{

  // set parameters
  double alphaM = params.get<double>("alpha_M");
  double alphaF = params.get<double>("alpha_F");
  double gamma  = params.get<double>("gamma");
  double dt     = params.get<double>("dt");

  double afgdt  = alphaF * gamma * dt;

  double tauM   = tau[0];
  double tauMp  = tau[1];
  double tauC   = tau[2];

  // further variables
  vector<double>            conv_c(iel);     /* linearisation of convect, convective part */
  Epetra_SerialDenseMatrix  conv_r(3,3*iel); /* linearisation of convect, reactive part   */
  Epetra_SerialDenseMatrix  viscs2(3,3*iel); /* viscous term including 2nd derivatives    */
  vector<double>            viscaf_old(3);     /* viscous term evaluated with old velocities      */

  bool supg =true;
  bool pstab=true;
  bool cstab=true;


  for (int i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix

                    n+alpha_F
       and u_old = u
                    (i)

       */
    conv_c[i] = derxy(0,i) * velintaf[0] + derxy(1,i) * velintaf[1]
      + derxy(2,i) * velintaf[2];

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /*       /                                     \
             |  u_old_x,x   u_old_x,y   u_old x,z  |
             |                                     |
             |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
             |                                     |
             |  u_old_z,x   u_old_z,y   u_old_z,z  |
             \                                     /

       with  N .. form function matrix

                    n+alpha_F
       and u_old = u
                    (i)

       */

    conv_r(0,3*i)   = vderxyaf(0,0)*funct[i];
    conv_r(0,3*i+1) = vderxyaf(0,1)*funct[i];
    conv_r(0,3*i+2) = vderxyaf(0,2)*funct[i];
    conv_r(1,3*i)   = vderxyaf(1,0)*funct[i];
    conv_r(1,3*i+1) = vderxyaf(1,1)*funct[i];
    conv_r(1,3*i+2) = vderxyaf(1,2)*funct[i];
    conv_r(2,3*i)   = vderxyaf(2,0)*funct[i];
    conv_r(2,3*i+1) = vderxyaf(2,1)*funct[i];
    conv_r(2,3*i+2) = vderxyaf(2,2)*funct[i];

   /*--- viscous term  grad * epsilon(u): ------------------------------*/
   /*      /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
       + - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
         2 |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

       with N_x .. x-line of N
            N_y .. y-line of N                                          */

    viscs2(0,3*i  ) = 0.5 * (2.0 * derxy2(0,i) + derxy2(1,i) + derxy2(2,i));
    viscs2(0,3*i+1) = 0.5 *  derxy2(3,i);
    viscs2(0,3*i+2) = 0.5 *  derxy2(4,i);
    viscs2(1,3*i  ) = 0.5 *  derxy2(3,i);
    viscs2(1,3*i+1) = 0.5 * (derxy2(0,i) + 2.0 * derxy2(1,i) + derxy2(2,i));
    viscs2(1,3*i+2) = 0.5 *  derxy2(5,i);
    viscs2(2,3*i  ) = 0.5 *  derxy2(4,i);
    viscs2(2,3*i+1) = 0.5 *  derxy2(5,i);
    viscs2(2,3*i+2) = 0.5 * (derxy2(0,i) + derxy2(1,i) + 2.0 * derxy2(2,i));


    /* Viscous term  div epsilon(u_old)

                    n+alpha_F
      with u_old = u
                    (i)

    */
    viscaf_old[0] = vderxyaf2(0,0) + 0.5 * ( vderxyaf2(0,1) + vderxyaf2(1,3)
                                             + vderxyaf2(0,2) + vderxyaf2(2,4));
    viscaf_old[1] = vderxyaf2(1,1) + 0.5 * ( vderxyaf2(1,0) + vderxyaf2(0,3)
                                             + vderxyaf2(1,2) + vderxyaf2(2,5));
    viscaf_old[2] = vderxyaf2(2,2) + 0.5 * ( vderxyaf2(2,0) + vderxyaf2(0,4)
                                             + vderxyaf2(2,1) + vderxyaf2(1,5));



  } // end of loop over nodes of element

#define conv_c_(j)     conv_c[j]
#define conv_r_(i,j,k) conv_r(i,3*(k)+j)
#define viscs2_(i,j,k) viscs2(i,3*(k)+j)

  for (int ui=0; ui<iel; ++ui) // loop columns (solution)
  {
    for (int vi=0; vi<iel; ++vi)  // loop rows (test functions)
    {
      //---------------------------------------------------------------
      //
      //                       GALERKIN PART
      //
      //---------------------------------------------------------------

      /*
       inertia term (intermediate)

       factor: +alphaM

                 /          \
                |            |
                |  Dacc , v  |
                |            |
                 \          /
      */

      elemat(vi*4    , ui*4    ) += fac*alphaM*funct[ui]*funct[vi] ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*alphaM*funct[ui]*funct[vi] ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*alphaM*funct[ui]*funct[vi] ;

#if 1
      /* convection (intermediate) */

      /*  factor: +alphaF*gamma*dt

                 /                                                    \
                |  / n+af       \          /            \   n+af       |
                | | u    o nabla | Dacc + | Dacc o nabla | u      , v  |
                |  \            /          \            /              |
                 \                                                    /
      */

      elemat(vi*4    , ui*4    ) += fac*afgdt*funct[vi]*(conv_r_(0, 0, ui)+conv_c_(ui)) ;
      elemat(vi*4    , ui*4 + 1) += fac*afgdt*funct[vi]*(conv_r_(0, 1, ui)            ) ;
      elemat(vi*4    , ui*4 + 2) += fac*afgdt*funct[vi]*(conv_r_(0, 2, ui)            ) ;
      elemat(vi*4 + 1, ui*4    ) += fac*afgdt*funct[vi]*(conv_r_(1, 0, ui)            ) ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*afgdt*funct[vi]*(conv_r_(1, 1, ui)+conv_c_(ui)) ;
      elemat(vi*4 + 1, ui*4 + 2) += fac*afgdt*funct[vi]*(conv_r_(1, 2, ui)            ) ;
      elemat(vi*4 + 2, ui*4    ) += fac*afgdt*funct[vi]*(conv_r_(2, 0, ui)            ) ;
      elemat(vi*4 + 2, ui*4 + 1) += fac*afgdt*funct[vi]*(conv_r_(2, 1, ui)            ) ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*afgdt*funct[vi]*(conv_r_(2, 2, ui)+conv_c_(ui)) ;
#endif

      /* pressure (implicit) */

      /*  factor: -1

                 /                \
                |                  |
                |  Dp , nabla o v  |
                |                  |
                 \                /
      */

      elemat(vi*4    , ui*4 + 3) -= fac*funct[ui]*derxy(0, vi) ;
      elemat(vi*4 + 1, ui*4 + 3) -= fac*funct[ui]*derxy(1, vi) ;
      elemat(vi*4 + 2, ui*4 + 3) -= fac*funct[ui]*derxy(2, vi) ;

      /* viscous term (intermediate) */

      /*  factor: +2*nu*alphaF*gamma*dt

                 /                          \
                |       /    \         / \   |
                |  eps | Dacc | , eps | v |  |
                |       \    /         \ /   |
                 \                          /
      */

      elemat(vi*4    , ui*4    ) += visc*afgdt*fac*(2.0*derxy(0,ui)*derxy(0,vi)
                                                    +
                                                    derxy(1,ui)*derxy(1,vi)
                                                    +
                                                    derxy(2,ui)*derxy(2,vi)) ;
      elemat(vi*4    , ui*4 + 1) += visc*afgdt*fac*derxy(0,ui)*derxy(1,vi) ;
      elemat(vi*4    , ui*4 + 2) += visc*afgdt*fac*derxy(0,ui)*derxy(2,vi) ;
      elemat(vi*4 + 1, ui*4    ) += visc*afgdt*fac*derxy(0,vi)*derxy(1,ui) ;
      elemat(vi*4 + 1, ui*4 + 1) += visc*afgdt*fac*(derxy(0,ui)*derxy(0,vi)
                                                    +
                                                    2.0*derxy(1,ui)*derxy(1,vi)
                                                    +
                                                    derxy(2,ui)*derxy(2,vi)) ;
      elemat(vi*4 + 1, ui*4 + 2) += visc*afgdt*fac*derxy(1,ui)*derxy(2,vi) ;
      elemat(vi*4 + 2, ui*4    ) += visc*afgdt*fac*derxy(0,vi)*derxy(2,ui) ;
      elemat(vi*4 + 2, ui*4 + 1) += visc*afgdt*fac*derxy(1,vi)*derxy(2,ui) ;
      elemat(vi*4 + 2, ui*4 + 2) += visc*afgdt*fac*(derxy(0,ui)*derxy(0,vi)
                                                    +
                                                    derxy(1,ui)*derxy(1,vi)
                                                    +
                                                    2.0*derxy(2,ui)*derxy(2,vi)) ;

      /* continuity equation (implicit) */

      /*  factor: +gamma*dt

                 /                  \
                |                    |
                | nabla o Dacc  , q  |
                |                    |
                 \                  /
      */

      elemat(vi*4 + 3, ui*4    ) += fac*gamma*dt*funct[vi]*derxy(0,ui) ;
      elemat(vi*4 + 3, ui*4 + 1) += fac*gamma*dt*funct[vi]*derxy(1,ui) ;
      elemat(vi*4 + 3, ui*4 + 2) += fac*gamma*dt*funct[vi]*derxy(2,ui) ;


      //---------------------------------------------------------------
      //
      //                     STABILISATION PART
      //                    PRESSURE STABILISATION
      //
      //---------------------------------------------------------------
      if(pstab)
      {
      /* pressure stabilisation --- inertia    */

      /* factor: +alphaM*tauMp

                 /                \
                |                  |
                |  Dacc , nabla q  |
                |                  |
                 \                /
      */

      elemat(vi*4 + 3, ui*4    ) += fac*alphaM*tauMp*funct[ui]*derxy(0,vi) ;
      elemat(vi*4 + 3, ui*4 + 1) += fac*alphaM*tauMp*funct[ui]*derxy(1,vi) ;
      elemat(vi*4 + 3, ui*4 + 2) += fac*alphaM*tauMp*funct[ui]*derxy(2,vi) ;

#if 1
      /* pressure stabilisation --- convection */


      /*  factor: +alphaF*gamma*dt*tauMp

                 /                                                          \
                |  / n+af       \          /            \   n+af             |
                | | u    o nabla | Dacc + | Dacc o nabla | u      , nabla q  |
                |  \            /          \            /                    |
                 \                                                          /
      */

      elemat(vi*4 + 3, ui*4    ) += fac*afgdt*tauMp*
                                    (conv_c_(ui)*derxy(0,vi)
                                     +
                                     derxy(0,vi)*conv_r_(0,0,ui)
                                     +
                                     derxy(1,vi)*conv_r_(1,0,ui)
                                     +
                                     derxy(2,vi)*conv_r_(2,0,ui)) ;
      elemat(vi*4 + 3, ui*4 + 1) += fac*afgdt*tauMp*
                                    (conv_c_(ui)*derxy(1,vi)
                                     +
                                     derxy(0,vi)*conv_r_(0,1,ui)
                                     +
                                     derxy(1,vi)*conv_r_(1,1,ui)
                                     +
                                     derxy(2,vi)*conv_r_(2,1,ui)) ;
      elemat(vi*4 + 3, ui*4 + 2) += fac*afgdt*tauMp*
                                    (conv_c_(ui)*derxy(2,vi)
                                     +
                                     derxy(0,vi)*conv_r_(0,2,ui)
                                     +
                                     derxy(1,vi)*conv_r_(1,2,ui)
                                     +
                                     derxy(2,vi)*conv_r_(2,2,ui)) ;
#endif

      /* pressure stabilisation --- diffusion  */


      /* factor: -2*nu*alphaF*gamma*dt*tauMp

                 /                                \
                |               /    \             |
                |  nabla o eps | Dacc | , nabla q  |
                |               \    /             |
                 \                                /
      */

      elemat(vi*4 + 3, ui*4    ) -= 2.0*visc*fac*afgdt*tauMp*
                                    (derxy(0,vi)*viscs2_(0,0,ui)
                                     +
                                     derxy(1,vi)*viscs2_(0,1,ui)
                                     +
                                     derxy(2,vi)*viscs2_(0,2,ui)) ;
      elemat(vi*4 + 3, ui*4 + 1) -= 2.0*visc*fac*afgdt*tauMp*
                                    (derxy(0,vi)*viscs2_(0,1,ui)
                                     +
                                     derxy(1,vi)*viscs2_(1,1,ui)
                                     +
                                     derxy(2,vi)*viscs2_(1,2,ui)) ;
      elemat(vi*4 + 3, ui*4 + 2) -= 2.0*visc*fac*afgdt*tauMp*
                                    (derxy(0,vi)*viscs2_(0,2,ui)
                                     +
                                     derxy(1,vi)*viscs2_(1,2,ui)
                                     +
                                     derxy(2,vi)*viscs2_(2,2,ui)) ;

      /* pressure stabilisation --- pressure   */

      /* factor: +tauMp

                 /                   \
                |                     |
                |  nabla p , nabla q  |
                |                     |
                 \                   /
      */

      elemat(vi*4 + 3, ui*4 + 3) += fac*tauMp*
                                    (derxy(0,ui)*derxy(0,vi)
                                     +
                                     derxy(1,ui)*derxy(1,vi)
                                     +
                                     derxy(2,ui)*derxy(2,vi)) ;
      }
      //---------------------------------------------------------------
      //
      //                     STABILISATION PART
      //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
      //
      //---------------------------------------------------------------
      if(supg)
      {
      /* SUPG stabilisation --- inertia     */

      /* factor: +alphaF*gamma*dt*tauM

                 /                               \
                |     n+am     /            \     |
                |  acc      , | Dacc o nabla | v  |
                |              \            /     |
                 \                               /
      */

      elemat(vi*4    , ui*4    ) += fac*afgdt*tauM*funct[ui]*accintam[0]*derxy(0,vi) ;
      elemat(vi*4    , ui*4 + 1) += fac*afgdt*tauM*funct[ui]*accintam[0]*derxy(1,vi) ;
      elemat(vi*4    , ui*4 + 2) += fac*afgdt*tauM*funct[ui]*accintam[0]*derxy(2,vi) ;
      elemat(vi*4 + 1, ui*4    ) += fac*afgdt*tauM*funct[ui]*accintam[1]*derxy(0,vi) ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*afgdt*tauM*funct[ui]*accintam[1]*derxy(1,vi) ;
      elemat(vi*4 + 1, ui*4 + 2) += fac*afgdt*tauM*funct[ui]*accintam[1]*derxy(2,vi) ;
      elemat(vi*4 + 2, ui*4)     += fac*afgdt*tauM*funct[ui]*accintam[2]*derxy(0,vi) ;
      elemat(vi*4 + 2, ui*4 + 1) += fac*afgdt*tauM*funct[ui]*accintam[2]*derxy(1,vi) ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*afgdt*tauM*funct[ui]*accintam[2]*derxy(2,vi) ;

      /* factor: +alphaM*tauM

                 /                           \
                |          / n+af       \     |
                |  Dacc , | u    o nabla | v  |
                |          \            /     |
                 \                           /
      */

      elemat(vi*4    , ui*4 + 3) += fac*alphaM*tauM*conv_c_(vi)*funct[ui] ;
      elemat(vi*4 + 1, ui*4 + 3) += fac*alphaM*tauM*conv_c_(vi)*funct[ui] ;
      elemat(vi*4 + 2, ui*4 + 3) += fac*alphaM*tauM*conv_c_(vi)*funct[ui] ;

#if 1
      /* SUPG stabilisation --- convection  */

      /* factor: +alphaF*gamma*dt*tauM

                 /                                               \
                |    / n+af        \   n+af    /            \     |
                |   | u     o nabla | u     , | Dacc o nabla | v  |
                |    \             /           \            /     |
                 \                                               /

                 /                                               \
                |    / n+af        \          / n+af        \     |
                |   | u     o nabla | Dacc , | u     o nabla | v  |
                |    \             /          \             /     |
                 \                                               /

                 /                                               \
                |    /            \   n+af    / n+af        \     |
                |   | Dacc o nabla | u     , | u     o nabla | v  |
                |    \            /           \             /     |
                 \                                               /
      */

      elemat(vi*4, ui*4)         += fac*afgdt*tauM*
                                    (conv_c_(ui)*conv_c_(vi)
                                     +
                                     conv_c_(vi)*conv_r_(0, 0, ui)
                                     +
                                     velintaf[0]*derxy(0, vi)*conv_r_(0, 0, ui)
                                     +
                                     velintaf[1]*derxy(0, vi)*conv_r_(0, 1, ui)
                                     +
                                     velintaf[2]*derxy(0, vi)*conv_r_(0, 2, ui)) ;
      elemat(vi*4, ui*4 + 1)     += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(0, 1, ui)
                                     +
                                     velintaf[0]*derxy(1, vi)*conv_r_(0, 0, ui)
                                     +
                                     velintaf[1]*derxy(1, vi)*conv_r_(0, 1, ui)
                                     +
                                     velintaf[2]*derxy(1, vi)*conv_r_(0, 2, ui)) ;
      elemat(vi*4, ui*4 + 2)     += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(0, 2, ui)
                                     +
                                     velintaf[0]*derxy(2, vi)*conv_r_(0, 0, ui)
                                     +
                                     velintaf[1]*derxy(2, vi)*conv_r_(0, 1, ui)
                                     +
                                     velintaf[2]*derxy(2, vi)*conv_r_(0, 2, ui)) ;
      elemat(vi*4 + 1, ui*4)     += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(1, 0, ui)
                                     +
                                     velintaf[0]*derxy(0, vi)*conv_r_(1, 0, ui)
                                     +
                                     velintaf[1]*derxy(0, vi)*conv_r_(1, 1, ui)
                                     +
                                     velintaf[2]*derxy(0, vi)*conv_r_(1, 2, ui)) ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*afgdt*tauM*
                                    (conv_c_(ui)*conv_c_(vi)
                                     +
                                     conv_c_(vi)*conv_r_(1, 1, ui)
                                     +
                                     velintaf[0]*derxy(1, vi)*conv_r_(1, 0, ui)
                                     +
                                     velintaf[1]*derxy(1, vi)*conv_r_(1, 1, ui)
                                     +
                                     velintaf[2]*derxy(1, vi)*conv_r_(1, 2, ui)) ;
      elemat(vi*4 + 1, ui*4 + 2) += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(1, 2, ui)
                                     +
                                     velintaf[0]*derxy(2, vi)*conv_r_(1, 0, ui)
                                     +
                                     velintaf[1]*derxy(2, vi)*conv_r_(1, 1, ui)
                                     +
                                     velintaf[2]*derxy(2, vi)*conv_r_(1, 2, ui)) ;
      elemat(vi*4 + 2, ui*4)     += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(2, 0, ui)
                                     +
                                     velintaf[0]*derxy(0, vi)*conv_r_(2, 0, ui)
                                     +
                                     velintaf[1]*derxy(0, vi)*conv_r_(2, 1, ui)
                                     +
                                     velintaf[2]*derxy(0, vi)*conv_r_(2, 2, ui)) ;
      elemat(vi*4 + 2, ui*4 + 1) += fac*afgdt*tauM*
                                    (conv_c_(vi)*conv_r_(2, 1, ui)
                                     +
                                     velintaf[0]*derxy(1, vi)*conv_r_(2, 0, ui)
                                     +
                                     velintaf[1]*derxy(1, vi)*conv_r_(2, 1, ui)
                                     +
                                     velintaf[2]*derxy(1, vi)*conv_r_(2, 2, ui)) ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*afgdt*tauM*
                                    (conv_c_(ui)*conv_c_(vi)
                                     +
                                     conv_c_(vi)*conv_r_(2, 2, ui)
                                     +
                                     velintaf[0]*derxy(2, vi)*conv_r_(2, 0, ui)
                                     +
                                     velintaf[1]*derxy(2, vi)*conv_r_(2, 1, ui)
                                     +
                                     velintaf[2]*derxy(2, vi)*conv_r_(2, 2, ui)) ;
#endif

      /* SUPG stabilisation --- diffusion   */

      /* factor: -2*nu*alphaF*gamma*dt*tauM

                 /                                            \
                |               / n+af \    /            \     |
                |  nabla o eps | u      |, | Dacc o nabla | v  |
                |               \      /    \            /     |
                 \                                            /


                 /                                            \
                |               /     \    / n+af        \     |
                |  nabla o eps | Dacc  |, | u     o nabla | v  |
                |               \     /    \             /     |
                 \                                            /
      */

      elemat(vi*4, ui*4)         -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(0, 0, ui)
                                     +
                                     funct[ui]*viscaf_old[0]*derxy(0, vi)) ;
      elemat(vi*4, ui*4 + 1)     -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(0, 1, ui)
                                     +
                                     funct[ui]*viscaf_old[0]*derxy(1, vi)) ;
      elemat(vi*4, ui*4 + 2)     -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(0, 2, ui)
                                     +
                                     funct[ui]*viscaf_old[0]*derxy(2, vi)) ;
      elemat(vi*4 + 1, ui*4)     -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(0, 1, ui)
                                     +
                                     funct[ui]*viscaf_old[1]*derxy(0, vi)) ;
      elemat(vi*4 + 1, ui*4 + 1) -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(1, 1, ui)
                                     +
                                     funct[ui]*viscaf_old[1]*derxy(1, vi)) ;
      elemat(vi*4 + 1, ui*4 + 2) -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(1, 2, ui)
                                     +
                                     funct[ui]*viscaf_old[1]*derxy(2, vi)) ;
      elemat(vi*4 + 2, ui*4)     -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(0, 2, ui)
                                     +
                                     funct[ui]*viscaf_old[2]*derxy(0, vi)) ;
      elemat(vi*4 + 2, ui*4 + 1) -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(1, 2, ui)
                                     +
                                     funct[ui]*viscaf_old[2]*derxy(1, vi)) ;
      elemat(vi*4 + 2, ui*4 + 2) -= 2.0*visc*fac*afgdt*tauM*
                                    (conv_c_(vi)*viscs2_(2, 2, ui)
                                     +
                                     funct[ui]*viscaf_old[2]*derxy(2, vi)) ;

      /* SUPG stabilisation --- pressure    */

      /* factor: +alphaF*gamma*dt*tauM

                 /                                 \
                |         n+1    /            \     |
                |  nabla p    , | Dacc o nabla | v  |
                |                \            /     |
                 \                                 /
      */

      elemat(vi*4    , ui*4    ) += fac*afgdt*tauM*funct[ui]*pderxynp[0]*derxy(0,vi) ;
      elemat(vi*4    , ui*4 + 1) += fac*afgdt*tauM*funct[ui]*pderxynp[0]*derxy(1,vi) ;
      elemat(vi*4    , ui*4 + 2) += fac*afgdt*tauM*funct[ui]*pderxynp[0]*derxy(2,vi) ;
      elemat(vi*4 + 1, ui*4    ) += fac*afgdt*tauM*funct[ui]*pderxynp[1]*derxy(0,vi) ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*afgdt*tauM*funct[ui]*pderxynp[1]*derxy(1,vi) ;
      elemat(vi*4 + 1, ui*4 + 2) += fac*afgdt*tauM*funct[ui]*pderxynp[1]*derxy(2,vi) ;
      elemat(vi*4 + 2, ui*4    ) += fac*afgdt*tauM*funct[ui]*pderxynp[2]*derxy(0,vi) ;
      elemat(vi*4 + 2, ui*4 + 1) += fac*afgdt*tauM*funct[ui]*pderxynp[2]*derxy(1,vi) ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*afgdt*tauM*funct[ui]*pderxynp[2]*derxy(2,vi) ;

      /* factor: +tauM

                 /                               \
                |              / n+af       \     |
                |  nabla Dp , | u    o nabla | v  |
                |              \            /     |
                 \                               /
      */

      elemat(vi*4    , ui*4 + 3) += fac*tauM*conv_c_(vi)*derxy(0, ui) ;
      elemat(vi*4 + 1, ui*4 + 3) += fac*tauM*conv_c_(vi)*derxy(1, ui) ;
      elemat(vi*4 + 2, ui*4 + 3) += fac*tauM*conv_c_(vi)*derxy(2, ui) ;
      }



      //---------------------------------------------------------------
      //
      //                     STABILISATION PART
      //            VISCOUS STABILISATION TERMS FOR USFEM
      //
      //---------------------------------------------------------------




      //---------------------------------------------------------------
      //
      //                     STABILISATION PART
      //                  CONTINUITY STABILISATION
      //
      //---------------------------------------------------------------

      /*  factor: +gamma*dt*tauC

                 /                          \
                |                            |
                | nabla o Dacc  , nabla o v  |
                |                            |
                 \                          /
      */

      if(cstab)
      {
      elemat(vi*4    , ui*4    ) += fac*gamma*dt*tauC*derxy(0,ui)*derxy(0,vi) ;
      elemat(vi*4    , ui*4 + 1) += fac*gamma*dt*tauC*derxy(0,vi)*derxy(1,ui) ;
      elemat(vi*4    , ui*4 + 2) += fac*gamma*dt*tauC*derxy(0,vi)*derxy(2,ui) ;
      elemat(vi*4 + 1, ui*4    ) += fac*gamma*dt*tauC*derxy(0,ui)*derxy(1,vi) ;
      elemat(vi*4 + 1, ui*4 + 1) += fac*gamma*dt*tauC*derxy(1,ui)*derxy(1,vi) ;
      elemat(vi*4 + 1, ui*4 + 2) += fac*gamma*dt*tauC*derxy(1,vi)*derxy(2,ui) ;
      elemat(vi*4 + 2, ui*4    ) += fac*gamma*dt*tauC*derxy(0,ui)*derxy(2,vi) ;
      elemat(vi*4 + 2, ui*4 + 1) += fac*gamma*dt*tauC*derxy(1,ui)*derxy(2,vi) ;
      elemat(vi*4 + 2, ui*4 + 2) += fac*gamma*dt*tauC*derxy(2,ui)*derxy(2,vi) ;
      }
    } // end loop vi
  } // end loop ui

#undef conv_c_
#undef conv_r_
#undef viscs2_

  return;
} // end of DRT:Elements:Fluid3:f3_genalpha_calmat


/*----------------------------------------------------------------------*
 |  evaluate fluid rhs (residual) for generalised alpha                 |
 |                            (private)                      gammi 06/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Fluid3::f3_genalpha_calrhs(
    Epetra_SerialDenseVector&  elevec,
    vector<double>&  	       accintam,
    vector<double>&  	       velintaf,
    Epetra_SerialDenseMatrix&  vderxyaf,
    Epetra_SerialDenseMatrix&  vderxyaf2,
    vector<double>&  	       velintnp,
    Epetra_SerialDenseMatrix&  vderxynp,
    double&                    prenp,
    vector<double>&  	       pderxynp,
    vector<double>&  	       edeadaf,
    vector<double>&  	       funct,
    Epetra_SerialDenseMatrix&  derxy,
    Epetra_SerialDenseMatrix&  derxy2,
    vector<double>&  	       tau,
    double&                    fac,
    const double&              visc,
    const int&                 iel,
    ParameterList& 	       params)
{

  // set parameters
  double tauM   = tau[0];
  double tauMp  = tau[1];
  double tauC   = tau[2];

  // further variables
  vector<double>            conv_c(iel);     /* linearisation of convect, convective part */
  Epetra_SerialDenseMatrix  conv_r(3,3*iel); /* linearisation of convect, reactive part   */
  vector<double>            conv_old(3);     /* convective term evaluated with old velocities */
  vector<double>            viscaf_old(3);     /* viscous term evaluated with old velocities      */

  bool supg =true;
  bool pstab=true;
  bool cstab=true;

  for (int i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix

                    n+alpha_F
       and u_old = u
                    (i)

       */
    conv_c[i] = derxy(0,i) * velintaf[0] + derxy(1,i) * velintaf[1]
      + derxy(2,i) * velintaf[2];

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /*       /                                     \
             |  u_old_x,x   u_old_x,y   u_old x,z  |
             |                                     |
             |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
             |                                     |
             |  u_old_z,x   u_old_z,y   u_old_z,z  |
             \                                     /

       with  N .. form function matrix

                    n+alpha_F
       and u_old = u
                    (i)

       */

    conv_r(0,3*i)   = vderxyaf(0,0)*funct[i];
    conv_r(0,3*i+1) = vderxyaf(0,1)*funct[i];
    conv_r(0,3*i+2) = vderxyaf(0,2)*funct[i];
    conv_r(1,3*i)   = vderxyaf(1,0)*funct[i];
    conv_r(1,3*i+1) = vderxyaf(1,1)*funct[i];
    conv_r(1,3*i+2) = vderxyaf(1,2)*funct[i];
    conv_r(2,3*i)   = vderxyaf(2,0)*funct[i];
    conv_r(2,3*i+1) = vderxyaf(2,1)*funct[i];
    conv_r(2,3*i+2) = vderxyaf(2,2)*funct[i];

  } // end of loop over nodes of element


  /* Viscous term  div epsilon(u_old) */
  viscaf_old[0] = vderxyaf2(0,0) + 0.5 * ( vderxyaf2(0,1) + vderxyaf2(1,3)
                                       + vderxyaf2(0,2) + vderxyaf2(2,4));
  viscaf_old[1] = vderxyaf2(1,1) + 0.5 * ( vderxyaf2(1,0) + vderxyaf2(0,3)
                                       + vderxyaf2(1,2) + vderxyaf2(2,5));
  viscaf_old[2] = vderxyaf2(2,2) + 0.5 * ( vderxyaf2(2,0) + vderxyaf2(0,4)
                                       + vderxyaf2(2,1) + vderxyaf2(1,5));


  /* Convective term  u_old * grad u_old: */
  conv_old[0] = vderxyaf(0,0) * velintaf[0] + vderxyaf(0,1) * velintaf[1]
    + vderxyaf(0,2) * velintaf[2];
  conv_old[1] = vderxyaf(1,0) * velintaf[0] + vderxyaf(1,1) * velintaf[1]
    + vderxyaf(1,2) * velintaf[2];
  conv_old[2] = vderxyaf(2,0) * velintaf[0] + vderxyaf(2,1) * velintaf[1]
    + vderxyaf(2,2) * velintaf[2];


#define conv_c_(j)     conv_c[j]
#define conv_r_(i,j,k) conv_r(i,3*(k)+j)
#define conv_old_(j)   conv_old[j]

  for (int vi=0; vi<iel; ++vi)  // loop rows (test functions)
  {
    //---------------------------------------------------------------
    //
    //                       GALERKIN PART
    //
    //---------------------------------------------------------------

    /* inertia terms */

    /*  factor: +1

               /             \
              |     n+am      |
              |  acc     , v  |
              |               |
               \             /
    */

    elevec[vi*4    ] -= fac*funct[vi]*accintam[0] ;
    elevec[vi*4 + 1] -= fac*funct[vi]*accintam[1] ;
    elevec[vi*4 + 2] -= fac*funct[vi]*accintam[2] ;

#if 1
    /* convection */

    /*  factor: +1

               /                             \
              |  / n+af       \    n+af       |
              | | u    o nabla |  u      , v  |
              |  \            /               |
               \                             /
    */

    elevec[vi*4    ] -= fac*(velintaf[0]*conv_r_(0,0,vi)
                             +
                             velintaf[1]*conv_r_(0,1,vi)
                             +
                             velintaf[2]*conv_r_(0,2,vi)) ;
    elevec[vi*4 + 1] -= fac*(velintaf[0]*conv_r_(1,0,vi)
                             +
                             velintaf[1]*conv_r_(1,1,vi)
                             +
                             velintaf[2]*conv_r_(1,2,vi)) ;
    elevec[vi*4 + 2] -= fac*(velintaf[0]*conv_r_(2,0,vi)
                             +
                             velintaf[1]*conv_r_(2,1,vi)
                             +
                             velintaf[2]*conv_r_(2,2,vi)) ;
#endif

    /* pressure */

    /*  factor: -1

               /                  \
              |   n+1              |
              |  p    , nabla o v  |
              |                    |
               \                  /
    */

    elevec[vi*4    ] += fac*prenp*derxy(0,vi) ;
    elevec[vi*4 + 1] += fac*prenp*derxy(1,vi) ;
    elevec[vi*4 + 2] += fac*prenp*derxy(2,vi) ;

    /* viscous term */

    /*  factor: +2*nu

               /                            \
              |       / n+af \         / \   |
              |  eps | u      | , eps | v |  |
              |       \      /         \ /   |
               \                            /
    */

    elevec[vi*4    ] -= visc*fac*
                        (derxy(0,vi)*vderxyaf(0,0)*2.0
                         +
                         derxy(1,vi)*vderxyaf(0,1)
                         +
                         derxy(1,vi)*vderxyaf(1,0)
                         +
                         derxy(2,vi)*vderxyaf(0,2)
                         +
                         derxy(2,vi)*vderxyaf(2,0)) ;
    elevec[vi*4 + 1] -= visc*fac*
                        (derxy(0,vi)*vderxyaf(0,1)
                         +
                         derxy(0,vi)*vderxyaf(1,0)
                         +
                         derxy(1,vi)*vderxyaf(1,1)*2.0
                         +
                         derxy(2,vi)*vderxyaf(1,2)
                         +
                         derxy(2,vi)*vderxyaf(2,1)) ;
    elevec[vi*4 + 2] -= visc*fac*
                        (derxy(0,vi)*vderxyaf(0,2)
                         +
                         derxy(0,vi)*vderxyaf(2,0)
                         +
                         derxy(1,vi)*vderxyaf(1,2)
                         +
                         derxy(1,vi)*vderxyaf(2,1)
                         +
                         derxy(2,vi)*vderxyaf(2,2)*2.0) ;

    /* body force (dead load...) */

    /*  factor: -1

               /           \
              |   n+af      |
              |  f     , v  |
              |             |
               \           /
    */

    elevec[vi*4    ] += fac*edeadaf[0]*funct[vi];
    elevec[vi*4 + 1] += fac*edeadaf[1]*funct[vi];
    elevec[vi*4 + 2] += fac*edeadaf[2]*funct[vi];

    /* continuity equation */

    /*  factor: +1

               /                \
              |          n+1     |
              | nabla o u   , q  |
              |                  |
               \                /
    */

    elevec[vi*4 + 3] -= fac*(vderxynp(0,0)+vderxynp(1,1)+vderxynp(2,2))*funct[vi];

    //---------------------------------------------------------------
    //
    //                     STABILISATION PART
    //                    PRESSURE STABILISATION
    //
    //---------------------------------------------------------------
    if(pstab)
    {
    /* pressure stabilisation --- inertia    */

    /* factor: +tauMp

               /                  \
              |     n+am           |
              |  acc    , nabla q  |
              |                    |
               \                  /
    */

    elevec[vi*4 + 3] -= fac*tauMp*
                        (derxy(0,vi)*accintam[0]
                         +
                         derxy(1,vi)*accintam[1]
                         +
                         derxy(2,vi)*accintam[2]);

#if 1
    /* pressure stabilisation --- convection */

    /*  factor: +tauMp

               /                                   \
              |  / n+af       \    n+af             |
              | | u    o nabla |  u      , nabla q  |
              |  \            /                     |
               \                                   /
    */

    elevec[vi*4 + 3] -= fac*tauMp*
                        (conv_old_(0)*derxy(0,vi)
                         +
                         conv_old_(1)*derxy(1,vi)
                         +
                         conv_old_(2)*derxy(2,vi)) ;
#endif

    /* pressure stabilisation --- diffusion  */

    /* factor: -2*nu*tauMp

               /                                  \
              |               / n+af \             |
              |  nabla o eps | u      | , nabla q  |
              |               \      /             |
               \                                  /
    */

    elevec[vi*4 + 3] += 2.0*visc*fac*tauMp*
                        (viscaf_old[0]*derxy(0,vi)
                         +
                         viscaf_old[1]*derxy(1,vi)
                         +
                         viscaf_old[2]*derxy(2,vi)) ;

    /* pressure stabilisation --- pressure   */

    /* factor: +tauMp

               /                      \
              |         n+1            |
              |  nabla p    , nabla q  |
              |                        |
               \                      /
    */

    elevec[vi*4 + 3] -= fac*tauMp*
                        (pderxynp[0]*derxy(0,vi)
                         +
                         pderxynp[1]*derxy(1,vi)
                         +
                         pderxynp[2]*derxy(2,vi)) ;

    /* pressure stabilisation --- bodyforce   */

    /* factor: -tauMp

               /                 \
              |    n+af           |
              |  f     , nabla q  |
              |                   |
               \                 /
    */

    elevec[vi*4 + 3] += fac*tauMp*(edeadaf[0]*derxy(0,vi)
                                   +
                                   edeadaf[1]*derxy(1,vi)
                                   +
                                   edeadaf[2]*derxy(2,vi)) ;
    }
    //---------------------------------------------------------------
    //
    //                     STABILISATION PART
    //         SUPG STABILISATION FOR CONVECTION DOMINATED FLOWS
    //
    //---------------------------------------------------------------
    if(supg)
    {
    /* SUPG stabilisation --- inertia     */

    /* factor: +tauM

               /                              \
              |     n+am   / n+af        \     |
              |  acc    , | u     o nabla | v  |
              |            \             /     |
               \                              /
    */

    elevec[vi*4    ] -= fac*tauM*conv_c_(vi)*accintam[0] ;
    elevec[vi*4 + 1] -= fac*tauM*conv_c_(vi)*accintam[1] ;
    elevec[vi*4 + 2] -= fac*tauM*conv_c_(vi)*accintam[2] ;

#if 1
    /* SUPG stabilisation --- convection  */

    /* factor: +tauM

               /                                                \
              |    / n+af        \   n+af    / n+af        \     |
              |   | u     o nabla | u     , | u     o nabla | v  |
              |    \             /           \             /     |
               \                                                /
    */

    elevec[vi*4    ] -= fac*tauM*conv_c_(vi)*conv_old_(0) ;
    elevec[vi*4 + 1] -= fac*tauM*conv_c_(vi)*conv_old_(1) ;
    elevec[vi*4 + 2] -= fac*tauM*conv_c_(vi)*conv_old_(2) ;
#endif

    /* SUPG stabilisation --- diffusion   */


    /* factor: -2*nu*tauM

               /                                               \
              |               / n+af \      / n+af        \     |
              |  nabla o eps | u      |  , | u     o nabla | v  |
              |               \      /      \             /     |
               \                                               /
    */

    elevec[vi*4    ] += fac*tauM*2.0*visc*conv_c_(vi)*viscaf_old[0] ;
    elevec[vi*4 + 1] += fac*tauM*2.0*visc*conv_c_(vi)*viscaf_old[1] ;
    elevec[vi*4 + 2] += fac*tauM*2.0*visc*conv_c_(vi)*viscaf_old[2] ;

    /* SUPG stabilisation --- pressure    */

    /* factor: +tauM

               /                                  \
              |         n+1    / n+af        \     |
              |  nabla p    , | u     o nabla | v  |
              |                \             /     |
               \                                  /
    */

    elevec[vi*4    ] -= fac*tauM*conv_c_(vi)*pderxynp[0] ;
    elevec[vi*4 + 1] -= fac*tauM*conv_c_(vi)*pderxynp[1] ;
    elevec[vi*4 + 2] -= fac*tauM*conv_c_(vi)*pderxynp[2] ;

    /* SUPG stabilisation --- bodyforce   */

    /* factor: -tauM

               /                             \
              |   n+af    / n+af        \     |
              |  f     , | u     o nabla | v  |
              |           \             /     |
               \                             /
    */

    elevec[vi*4    ] += fac*tauM*conv_c_(vi)*edeadaf[0] ;
    elevec[vi*4 + 1] += fac*tauM*conv_c_(vi)*edeadaf[1] ;
    elevec[vi*4 + 2] += fac*tauM*conv_c_(vi)*edeadaf[2] ;
    }


    //---------------------------------------------------------------
    //
    //                     STABILISATION PART
    //                  CONTINUITY STABILISATION
    //
    //---------------------------------------------------------------


    /* factor: +tauC

               /                          \
              |           n+1              |
              |  nabla o u    , nabla o v  |
              |                            |
               \                          /
    */

    if(cstab)
    {
    elevec[vi*4    ] -= fac*tauC*derxy(0,vi)*(vderxynp(0,0)+vderxynp(1,1)+vderxynp(2,2)) ;
    elevec[vi*4 + 1] -= fac*tauC*derxy(1,vi)*(vderxynp(0,0)+vderxynp(1,1)+vderxynp(2,2)) ;
    elevec[vi*4 + 2] -= fac*tauC*derxy(2,vi)*(vderxynp(0,0)+vderxynp(1,1)+vderxynp(2,2)) ;
    }
  } // end loop vi

#undef conv_c_
#undef conv_r_
#undef conv_old_

  return;
} // end of DRT:Elements:Fluid3:f3_genalpha_calrhs



/*---------------------------------------------------------------------*
 |  calculate error for beltrami test problem (private)     gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_int_beltrami_err(
  vector<double>&           evelnp,
  vector<double>&           eprenp,
  struct _MATERIAL*         material,
  ParameterList& 	    params
  )
{

  /*-------------------------- add element error to "integrated" error */
  double velerr = params.get<double>("L2 integrated velocity error");
  double preerr = params.get<double>("L2 integrated pressure error");

  /*------------------------------------------------- set element data */
  const int iel = NumNode();

  int       intc=0;   /* "integration case" for tri for further infos
                         see f2_inpele.c and f2_intg.c                 */
  int       nir=0;    /* number of integration nodesin r direction     */
  int       nis=0;    /* number of integration nodesin s direction     */
  int       nit=0;    /* number of integration nodesin t direction     */
  int       ihoel=0;  /* flag for higher order elements                */
  int       icode=2;  /* flag for eveluation of shape functions        */
  double    fac;      /* total integration factor                      */


  vector<double>                funct(iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	deriv(3,iel);
  Epetra_SerialDenseMatrix 	deriv2(6,iel);

  double         		det;
  double         		e1, e2, e3;
  double         		facr=0.0, facs=0.0, fact=0.0;

  FLUID_DATA            	data;


  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  //------------------------------ set constants for analytical solution
  const double t = params.get("total time",-1.0);
  if (t<0)
  {
    dserror("beltrami: no total time for error calculation");
  }

  double a      = PI/4.0;
  double d      = PI/2.0;

  /* get viscosity ---*/
  const double  visc = material->m.fluid->viscosity;


  /*-------------------------------------------------- initialise ---*/
  // gaussian points
  f3_integration_points(data);



  switch (iel)
  {
      case 8: case 20: case 27:  /* --> hex - element */
        icode   = 3;
        ihoel   = 1;
        /* initialise integration */
        nir = ngp_[0];
        nis = ngp_[1];
        nit = ngp_[2];
        intc= 0;
        break;
      case 10: /* --> tet - element */
        icode   = 3;
        ihoel   = 1;
        /* do NOT break at this point!!! */
      case 4:    /* initialise integration */
        nir  = ngp_[0]; // for tets in ngp_[0] the number of gauss points is stored !
        nis  = 1;
        nit  = 1;
        intc = ngp_[1];
        break;
      default:
        dserror("typ unknown!");
  } // end switch (iel) //


  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  double         p;
  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);


  for (int lr=0;lr<nir;lr++)
  {
    for (int ls=0;ls<nis;ls++)
    {
      for (int lt=0;lt<nit;lt++)
      {
        /*------------- get values of  shape functions and their derivatives ---*/
        switch(iel)
        {
            case 8: case 20: case 27:   /* --> hex - element */
              e1   = data.qxg[lr][nir-1];
              facr = data.qwgt[lr][nir-1];
              e2   = data.qxg[ls][nis-1];
              facs = data.qwgt[ls][nis-1];
              e3   = data.qxg[lt][nit-1];
              fact = data.qwgt[lt][nit-1];
              f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
              break;
            case 4: case 10:   /* --> tet - element */
              e1   = data.txgr[lr][intc];
              facr = data.twgt[lr][intc];
              e2   = data.txgs[lr][intc];
              facs = ONE;
              e3   = data.txgt[lr][intc];
              fact = ONE;
              f3_shape_function(funct,deriv,deriv2,e1,e2,e3,iel,icode);
              break;
            default:
              facr = facs = fact = 0.0;
              e1 = e2 = e3 = 0.0;
              dserror("typ unknown!");
        } /* end switch (iel) */

          /*------------------------------------ compute Jacobian matrix */
        f3_jaco(xyze,deriv,xjm,&det,iel);
        fac = facr*facs*fact*det;

        /*---------------------- get velocity sol at integration point */
        for (int i=0;i<3;i++)
        {
          velint[i]=ZERO;
          for (int j=0;j<iel;j++)
          {
            velint[i] += funct[j]*evelnp[i+(3*j)];
          }
        } //end loop over i

          /*---------------------- get pressure sol at integration point */
        preint = 0;
        for (int i=0;i<iel;i++)
        {
          preint += funct[i]*eprenp[i];
        }

        /*---------------------- get velocity sol at integration point */
        for (int i=0;i<3;i++)
        {
          xint[i]=ZERO;
          for (int j=0;j<iel;j++)
          {
            xint[i] += funct[j]*xyze(i,j);
          }
        } //end loop over i


          // compute analytical pressure
        p = -a*a/2.0 *
          ( exp(2.0*a*xint[0])
            + exp(2.0*a*xint[1])
            + exp(2.0*a*xint[2])
            + 2.0 * sin(a*xint[0] + d*xint[1]) * cos(a*xint[2] + d*xint[0]) * exp(a*(xint[1]+xint[2]))
            + 2.0 * sin(a*xint[1] + d*xint[2]) * cos(a*xint[0] + d*xint[1]) * exp(a*(xint[2]+xint[0]))
            + 2.0 * sin(a*xint[2] + d*xint[0]) * cos(a*xint[1] + d*xint[2]) * exp(a*(xint[0]+xint[1]))
            )* exp(-2.0*visc*d*d*t);

        // compute analytical velocities
        u[0] = -a * ( exp(a*xint[0]) * sin(a*xint[1] + d*xint[2]) +
                      exp(a*xint[2]) * cos(a*xint[0] + d*xint[1]) ) * exp(-visc*d*d*t);
        u[1] = -a * ( exp(a*xint[1]) * sin(a*xint[2] + d*xint[0]) +
                      exp(a*xint[0]) * cos(a*xint[1] + d*xint[2]) ) * exp(-visc*d*d*t);
        u[2] = -a * ( exp(a*xint[2]) * sin(a*xint[0] + d*xint[1]) +
                      exp(a*xint[1]) * cos(a*xint[2] + d*xint[0]) ) * exp(-visc*d*d*t);

        // compute difference between analytical solution and numerical solution
        deltap = preint - p;

        for (int dim=0;dim<3;dim++)
        {
          deltavel[dim]=velint[dim]-u[dim];
        }

        // add square to L2 error
        for (int dim=0;dim<3;dim++)
        {
          velerr += deltavel[dim]*deltavel[dim]*fac;
        }
        preerr += deltap*deltap*fac;

      } /* end of loop over integration points lt*/
    } /* end of loop over integration points ls */
  } /* end of loop over integration points lr */


  // we use the parameterlist as a container to transport the calculated
  // errors from the elements to the dynamic routine

  params.set<double>("L2 integrated velocity error",velerr);
  params.set<double>("L2 integrated pressure error",preerr);

  return;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
