/*!
\file fluid3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
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
#include "../drt_mat/newtonianfluid.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


// converts a string into an Action for this element
DRT::Elements::Fluid3::ActionType DRT::Elements::Fluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");
  
  DRT::Elements::Fluid3::ActionType act = Fluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")      
    act = Fluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_genalpha_sysmat")
    act = Fluid3::calc_fluid_genalpha_sysmat;
  else if (action == "calc_fluid_genalpha_residual")
    act = Fluid3::calc_fluid_genalpha_residual;
  else if (action == "calc_fluid_beltrami_error")      
    act = Fluid3::calc_fluid_beltrami_error;
  else if (action == "calc_Shapefunction")
    act = Fluid3::calc_Shapefunction;
  else if (action == "calc_ShapeDeriv1")
    act = Fluid3::calc_ShapeDeriv1;
  else if (action == "calc_ShapeDeriv2")
    act = Fluid3::calc_ShapeDeriv2;
  else 
    dserror("Unknown type of action for Fluid3");
  return act;
}

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
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::Elements::Fluid3::ActionType act = convertStringToActionType(action);

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
        const int numnode = NumNode();
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
          const int numnode = NumNode();
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
      case calc_Shapefunction:
        shape_function_3D(elevec1,elevec2[0],elevec2[1],elevec2[2],this->Shape());
        break;
      case calc_ShapeDeriv1:
        shape_function_3D_deriv1(elemat1,elevec2[0],elevec2[1],elevec2[2],this->Shape());
        break;
      case calc_ShapeDeriv2:
        shape_function_3D_deriv2(elemat2,elevec2[0],elevec2[1],elevec2[2],this->Shape());
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
void DRT::Elements::Fluid3::f3_sys_mat(const vector<int>&        lm,
                                       const vector<double>&     evelnp,
                                       const vector<double>&     eprenp,
                                       const vector<double>&     evhist,
                                       const vector<double>&     edispnp,
                                       const vector<double>&     egridv,
                                       Epetra_SerialDenseMatrix* sys_mat,
                                       Epetra_SerialDenseVector* residual,
                                       struct _MATERIAL*         material,
                                       ParameterList&            params
  )
{
  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  // get node coordinates
  const Epetra_SerialDenseMatrix xyze = f3_getPositionArray(edispnp);
  
  // dead load in element nodes
  const double time = params.get("total time",-1.0);
  const Epetra_SerialDenseMatrix bodyforce = f3_getbodyforce(time,params);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  // declaration of variables
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	deriv(3,iel);
  Epetra_SerialDenseMatrix 	deriv2(6,iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	vderxy(3,3);
  vector<double>            pderxy(3);
  Epetra_SerialDenseMatrix 	vderxy2(3,6);
  Epetra_SerialDenseMatrix 	derxy(3,iel);
  Epetra_SerialDenseMatrix 	derxy2(6,iel);
  vector<double>            edeadng(3);
  vector<double>            histvec(3);   ///< history data at integration point
  vector<double>         	  velino(3);    ///< normed velocity at element centre
  vector<double>     	      velint(3);

  const double timefac=params.get<double>("time constant for integration",0.0);

  // get control parameter to switch between stationary and instationary problem
  const bool is_stationary = params.get<bool>("using stationary formulation",false);

  // stabilization parameter
  const vector<double> tau = f3_caltau(xyze,evelnp,distype,visc,iel,timefac,is_stationary);

  // flag for higher order elements
  const bool higher_order_ele = is_higher_order_element(distype);

  // gaussian points
//  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
//  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule_);

  // integration loop
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // coordiantes of the current integration point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    // shape functions and their derivatives
    Epetra_SerialDenseVector    funct(iel);
    Epetra_SerialDenseMatrix    deriv(3,iel);
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
          
    // get Jacobian matrix and determinant
    const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
    const double det = getDeterminante(xjm);
    const double fac = intpoints.qwgt[iquad]*det;

    // compute global derivates
    f3_gder(derxy,deriv,xjm,det,iel);

    // compute second global derivative
    if (higher_order_ele)
    {
      shape_function_3D_deriv2(deriv2,e1,e2,e3,distype);
      f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);

      // calculate 2nd velocity derivatives at integration point
      // former f3_vder2(vderxy2,derxy2,evelnp,iel);
      for (int i=0;i<6;i++)
      {
        vderxy2(0,i)=0.0;
        vderxy2(1,i)=0.0;
        vderxy2(2,i)=0.0;
        for (int inode=0;inode<iel;inode++)
        {
          vderxy2(0,i) += derxy2(i,inode)*evelnp[0+(3*inode)];
          vderxy2(1,i) += derxy2(i,inode)*evelnp[1+(3*inode)];
          vderxy2(2,i) += derxy2(i,inode)*evelnp[2+(3*inode)];
        }
      }
    }

    // get velocities (n+g,i) at integration point
    // expression for f3_veci(velint,funct,evelnp,iel);
    for (int isd=0;isd<NSD_;isd++)
    {
      velint[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        velint[isd] += funct[inode]*evelnp[isd+(3*inode)];
      }
    }

    // get history data (n,i) at integration point
    //expression for f3_veci(histvec,funct,evhist,iel);
    for (int isd=0;isd<NSD_;isd++)
    {
      histvec[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        histvec[isd] += funct[inode]*evhist[isd+(3*inode)];
      }
    } 

    // get velocity (np,i) derivatives at integration point
    // expression for f3_vder(vderxy,derxy,evelnp,iel);
    for (int isd=0;isd<NSD_;isd++)
    {
      vderxy(0,isd)=0.0;
      vderxy(1,isd)=0.0;
      vderxy(2,isd)=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        vderxy(0,isd) += derxy(isd,inode)*evelnp[0+(3*inode)];
        vderxy(1,isd) += derxy(isd,inode)*evelnp[1+(3*inode)];
        vderxy(2,isd) += derxy(isd,inode)*evelnp[2+(3*inode)];
      }
    }

    // get grid velocity at integration point
    vector<double>    gridvelint(NSD_);
    if (is_ale_)
    {
      for (int isd=0; isd<NSD_; isd++)
      {
        gridvelint[isd] = 0.;
        for (int inode=0; inode<iel; inode++)
        {
          gridvelint[isd] += derxy(isd,inode)*egridv[isd+(4*inode)];
        }
      }
    }
    else
    {
      gridvelint[0] = 0.0;
      gridvelint[1] = 0.0;
      gridvelint[2] = 0.0;
    }

    // get pressure gradients
    vector<double>    gradp(NSD_);
    gradp[0] = gradp[1] = gradp[2] = 0.0;
    for (int inode=0; inode<iel; inode++)
    {
      gradp[0] += derxy(0,inode) * eprenp[inode];
      gradp[1] += derxy(1,inode) * eprenp[inode];
      gradp[2] += derxy(2,inode) * eprenp[inode];
    }

    double press = 0.0;
    for (int inode=0;inode<iel;inode++)
    {
      press += funct[inode]*eprenp[inode];
    }

    // get bodyforce in gausspoint
    for (int isd=0;isd<NSD_;isd++)
    {
      edeadng[isd] = 0.0;
      for (int inode=0;inode<iel;inode++)
      {
        edeadng[isd]+= bodyforce(isd,inode)*funct[inode];
      }
    }

    // perform integration for entire matrix and rhs
    if(is_stationary==false)
        f3_calmat(*sys_mat,*residual,velint,histvec,gridvelint,
                press,vderxy,vderxy2,gradp,funct,tau,
                derxy,derxy2,edeadng,fac,visc,iel,params);
    else
        f3_calmat_stationary(*sys_mat,*residual,velint,histvec,gridvelint,
                press,vderxy,vderxy2,gradp,funct,tau,
                derxy,derxy2,edeadng,fac,visc,iel,params);

  } // end of loop over integration points

  return;
} // DRT::Elements::Fluid3::f3_sys_mat




Epetra_SerialDenseMatrix DRT::Elements::Fluid3::f3_getPositionArray(
          const vector<double>&   edispnp)
{
    
    const int iel = NumNode();
    Epetra_SerialDenseMatrix xyze(NSD_,iel);

    // get initial position
    for (int inode=0;inode<iel;inode++)
    {
      xyze(0,inode) = Nodes()[inode]->X()[0];
      xyze(1,inode) = Nodes()[inode]->X()[1];
      xyze(2,inode) = Nodes()[inode]->X()[2];
    }

    // add displacement, when fluid nodes move in the ALE case
    if (is_ale_)
    {
      for (int inode=0;inode<iel;inode++)
      {
        xyze(0,inode) += edispnp[4*inode];
        xyze(1,inode) += edispnp[4*inode+1];
        xyze(2,inode) += edispnp[4*inode+2];
      }
    }
    
    return xyze;
}    



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
  const DiscretizationType distype = this->Shape();
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
    const double acttime    = params.get<double>("time");
    const double dt         = params.get<double>("dt");
    const double alphaF     = params.get<double>("alpha_F");

    //         n+alpha_F     n+1
    //        t          = t     - (1-alpha_F) * dt

    const double timealphaF = acttime-(1-alphaF)*dt;

    const Epetra_SerialDenseMatrix bodyforce = f3_getbodyforce(timealphaF,params);

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

    const bool higher_order_ele = is_higher_order_element(distype);
    // gaussian points
    const GaussRule3D gaussrule = getOptimalGaussrule(distype);
    const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
          // gauss point coordinates
          const double e1 = intpoints.qxg[iquad][0];
          const double e2 = intpoints.qxg[iquad][1];
          const double e3 = intpoints.qxg[iquad][2];
          
          // shape functions and derivatives
          Epetra_SerialDenseVector  funct (iel);
          Epetra_SerialDenseMatrix 	deriv (3,iel);
          Epetra_SerialDenseMatrix 	deriv2(6,iel);
          Epetra_SerialDenseMatrix 	derxy (3,iel);
          Epetra_SerialDenseMatrix 	derxy2(6,iel);

          // intermediate accelerations (n+alpha_M)
          vector<double>     		    accintam (3);
          // intermediate velocities    (n+alpha_F) and its derivatives
          vector<double>     		    velintaf (3);
          Epetra_SerialDenseMatrix 	vderxyaf (3,3);
          Epetra_SerialDenseMatrix 	vderxyaf2(3,6);
          // new velocities for continuity equation and its derivatives
          vector<double>            velintnp (3);
          Epetra_SerialDenseMatrix 	vderxynp (3,3);
          // new pressure and its derivatives
          double                    prenp;
          vector<double>            pderxynp(3);

          // dead load
          vector<double>                edeadaf(3);

          shape_function_3D(funct,e1,e2,e3,distype);
          shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
          if (higher_order_ele)
          {
            shape_function_3D_deriv2(deriv2,e1,e2,e3,distype);
          }

          // get Jacobian matrix and determinant
          const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
          const double det = getDeterminante(xjm);

          // set total integration factor
          const double fac = intpoints.qwgt[iquad]*det;

          // compute global derivates
          f3_gder(derxy,deriv,xjm,det,iel);

          // compute second global derivative
          if (higher_order_ele)
          {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);
          }


          // get intermediate accelerations (n+1,i)  at integration point
          for (int i=0;i<3;i++)
          {
            accintam[i]=0.0;
            for (int j=0;j<iel;j++)
            {
              accintam[i] += funct[j]*myaccam[i+(3*j)];
            }
          }


          //get velocities (n+1,i)  at integration point
          for (int isd=0;isd<NSD_;isd++)
          {
            velintnp[isd]=0.0;
            for (int inode=0;inode<iel;inode++)
            {
              velintnp[isd] += funct[inode]*myvelnp[isd+(3*inode)];
            }
          }


          // get velocities (n+alpha_F,i) at integration point
          for (int isd=0;isd<3;isd++)
          {
            velintaf[isd]=0.0;
            for (int inode=0;inode<iel;inode++)
            {
              velintaf[isd] += funct[inode]*myvelaf[isd+(3*inode)];
            }
          }


          /*----- get velocity (n+1,i) derivatives at integration point */
          for (int isd=0;isd<3;isd++)
          {
            vderxynp(0,isd) = 0.0;
            vderxynp(1,isd) = 0.0;
            vderxynp(2,isd) = 0.0;
            for (int inode=0;inode<iel;inode++)
            {
              vderxynp(0,isd) += derxy(isd,inode)*myvelnp[0+(3*inode)];
              vderxynp(1,isd) += derxy(isd,inode)*myvelnp[1+(3*inode)];
              vderxynp(2,isd) += derxy(isd,inode)*myvelnp[2+(3*inode)];
            }
          }


          /*----------- get velocity (n+alpha_F,i) derivatives at
                                                      integration point */
          for (int i=0;i<3;i++)
          {
            vderxyaf(0,i)=0.0;
            vderxyaf(1,i)=0.0;
            vderxyaf(2,i)=0.0;
            for (int j=0;j<iel;j++)
            {
              vderxyaf(0,i) += derxy(i,j)*myvelaf[0+(3*j)];
              vderxyaf(1,i) += derxy(i,j)*myvelaf[1+(3*j)];
              vderxyaf(2,i) += derxy(i,j)*myvelaf[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*------calculate 2nd velocity derivatives at integration
                                                     point (n+alpha_F,i)*/
          if (higher_order_ele)
          {
            for (int i=0;i<6;i++)
            {
              vderxyaf2(0,i)=0.0;
              vderxyaf2(1,i)=0.0;
              vderxyaf2(2,i)=0.0;
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

    } // end of loop over integration points

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
  const DiscretizationType distype = this->Shape();
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
    const double acttime    = params.get<double>("time");
    const double dt         = params.get<double>("dt");
    const double alphaF     = params.get<double>("alpha_F");

    //         n+alpha_F     n+1
    //        t          = t     - (1-alpha_F) * dt

    const double timealphaF = acttime-(1-alphaF)*dt;
    const Epetra_SerialDenseMatrix bodyforce = f3_getbodyforce(timealphaF,params);

    /*---------------------------------------------- get viscosity ---*/
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double  visc = material->m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    vector<double>   tau(3); // stab parameters
    {
      int version =2; // evaluate stabilisation parameter for genalpha
                      // time integration

      f3_calc_stabpar(tau,iel,xyze,myvelaf,visc,params,version);
    }

    const bool higher_order_ele = is_higher_order_element(distype);
    // gaussian points
    const GaussRule3D gaussrule = getOptimalGaussrule(distype);
    const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
          /*------------------- declaration of gauss point variables ---*/
          const double e1 = intpoints.qxg[iquad][0];
          const double e2 = intpoints.qxg[iquad][1];
          const double e3 = intpoints.qxg[iquad][2];
          /*------------------- declaration of gauss point variables ---*/

          // shape functions and derivatives
          Epetra_SerialDenseVector  funct (iel);
          Epetra_SerialDenseMatrix 	deriv (3,iel);
          Epetra_SerialDenseMatrix 	deriv2(6,iel);
          Epetra_SerialDenseMatrix 	derxy (3,iel);
          Epetra_SerialDenseMatrix 	derxy2(6,iel);

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


          shape_function_3D(funct,e1,e2,e3,distype);
          shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
          if (higher_order_ele)
          {
              shape_function_3D_deriv2(deriv2,e1,e2,e3,distype);
          }

          // get Jacobian matrix and determinant
          const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
          const double det = getDeterminante(xjm);

          // integration constant
          const double fac = intpoints.qwgt[iquad]*det;

          /*---------------------------------- compute global derivates */
          f3_gder(derxy,deriv,xjm,det,iel);

          /*-------------------------- compute second global derivative */
          if (higher_order_ele)
          {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,iel);
          }


          /*-- get intermediate accelerations (n+am,i)  at integration
                                                                  point */
          for (int i=0;i<3;i++)
          {
            accintam[i]=0.0;
            for (int j=0;j<iel;j++)
            {
              accintam[i] += funct[j]*myaccam[i+(3*j)];
            }
          } //end loop over i


          /*-------------- get velocities (n+1,i)  at integration point */
          for (int i=0;i<3;i++)
          {
            velintnp[i]=0.0;
            for (int j=0;j<iel;j++)
            {
              velintnp[i] += funct[j]*myvelnp[i+(3*j)];
            }
          } //end loop over i


          /*--------- get velocities (n+alpha_F,i) at integration point */
          for (int i=0;i<3;i++)
          {
            velintaf[i]=0.0;
            for (int j=0;j<iel;j++)
            {
              velintaf[i] += funct[j]*myvelaf[i+(3*j)];
            }
          } //end loop over i


          /*----- get velocity (n+1,i) derivatives at integration point */
          for (int i=0;i<3;i++)
          {
            vderxynp(0,i)=0.0;
            vderxynp(1,i)=0.0;
            vderxynp(2,i)=0.0;
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
            vderxyaf(0,i)=0.0;
            vderxyaf(1,i)=0.0;
            vderxyaf(2,i)=0.0;
            for (int j=0;j<iel;j++)
            {
              vderxyaf(0,i) += derxy(i,j)*myvelaf[0+(3*j)];
              vderxyaf(1,i) += derxy(i,j)*myvelaf[1+(3*j)];
              vderxyaf(2,i) += derxy(i,j)*myvelaf[2+(3*j)];
            } /* end of loop over j */
          } /* end of loop over i */


          /*------calculate 2nd velocity derivatives at integration
                                                     point (n+alpha_F,i)*/
          if (higher_order_ele)
          {
            for (int i=0;i<6;i++)
            {
              vderxyaf2(0,i)=0.0;
              vderxyaf2(1,i)=0.0;
              vderxyaf2(2,i)=0.0;
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

    } // end of loop over integration points

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
  const DiscretizationType distype = this->Shape();
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
        Epetra_SerialDenseVector    funct(iel);
        Epetra_SerialDenseMatrix    deriv(3,iel);
        Epetra_SerialDenseMatrix    xjm(3,3);
        Epetra_SerialDenseMatrix    derxy(3,iel);
        double                      hk;
        double                      val, strle;
        vector<double>              velino(3); /* normed velocity at element centre */
        double                      mk=0.0;
        vector<double>              velint(3);

        /*------------------------------------------------- initialise ---*/
        // gaussian points
        // use one point gauss rule to calculate tau at element center
        GaussRule3D integrationrule_stabili;
        switch(distype)
        {
        case hex8: case hex20: case hex27:
            integrationrule_stabili = intrule_hex_1point;
            break;
        case tet4: case tet10:
            integrationrule_stabili = intrule_tet_1point;
            break;
        default:
            dserror("invalid discretization type for fluid3");
        }
        const double timefac = params.get<double>("dt",0.0) * params.get<double>("gamma",0.0);

        // gaussian points
        const IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_stabili);

        const double e1    = intpoints.qxg[0][0];
        const double e2    = intpoints.qxg[0][1];
        const double e3    = intpoints.qxg[0][2];
        const double wquad = intpoints.qwgt[0];
        DRT::Utils::shape_function_3D(funct,e1,e2,e3,distype);
        DRT::Utils::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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
          velint[i]=0.0;
          for (int j=0;j<iel;j++)
          {
            velint[i] += funct[j]*myvelnp[i+(3*j)];
          }
        } //end loop over i

        {
          double vel_norm, re1, re2, xi1, xi2;

          // get Jacobian matrix and determinant
          const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
          const double det = getDeterminante(xjm);
          const double vol=wquad*det;

          /* get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)*/
          hk = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);

          /*------------------------------------------------- get streamlength ---*/
          f3_gder(derxy,deriv,xjm,det,iel);
          val = 0.0;


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
            velino[1] = 0.0;
            velino[2] = 0.0;
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

            tau[0] = timefac* DSQR(strle) / (DSQR(strle)*xi1+(4.0 * timefac*visc/mk)*xi2);

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


            tau[1] = timefac * DSQR(hk) / (DSQR(hk) * xi1 + ( 4.0 * timefac * visc/mk) * xi2);

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

            tau[2] = vel_norm * hk * 0.5 * xi2;
          }

        }

      }
      break;
      default:
        dserror("unknown version of tau!");

  }

  return;
} // end of DRT:Elements:Fluid3:f3_calc_stabpar

//
// calculate stabilization parameter
//
vector<double> DRT::Elements::Fluid3::f3_caltau(
    const Epetra_SerialDenseMatrix&         xyze,
    const vector<double>&                   evelnp,
    const DRT::Element::DiscretizationType  distype,
    const double                            visc,
    const int                               numnode,
    const double                            timefac,
    const bool                              is_stationary
    )
{
    // use one point gauss rule to calculate tau at element center
    GaussRule3D integrationrule_stabili;
    switch(distype)
    {
    case hex8: case hex20: case hex27:
        integrationrule_stabili = intrule_hex_1point;
        break;
    case tet4: case tet10:
        integrationrule_stabili = intrule_tet_1point;
        break;
    default: 
        dserror("invalid discretization type for fluid3");
    }

    // gaussian points
    const IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_stabili);

    // shape functions and derivs at element center
    const double e1    = intpoints.qxg[0][0];
    const double e2    = intpoints.qxg[0][1];
    const double e3    = intpoints.qxg[0][2];
    const double wquad = intpoints.qwgt[0];
    
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(NSD_, numnode);
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
    
    // get element type constant for tau
    double mk=0.0;
    switch(distype)
    {
    case tet4: case hex8:
        mk = 0.333333333333333333333;
        break;
    case hex20: case hex27: case tet10:
        mk = 0.083333333333333333333;
        break;
    default: 
        dserror("type unknown!\n");
    }
    
    // get velocities at element center
    vector<double>  velint(NSD_);
    for (int isd=0;isd<NSD_;isd++)
    {
        velint[isd]=0.0;
        for (int inode=0;inode<numnode;inode++)
        {
            velint[isd] += funct[inode]*evelnp[isd+(NSD_*inode)];
        }
    }

    // get Jacobian matrix and determinant
    const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,numnode);
    const double det = getDeterminante(xjm);
    const double vol = wquad*det;

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = pow((SIX*vol/PI),(1.0/3.0))/sqrt(3.0);

    // get derivatives
    Epetra_SerialDenseMatrix    derxy(NSD_, numnode);
    f3_gder(derxy,deriv,xjm,det,numnode);

    // get velocity norm
    const double vel_norm=sqrt( velint[0]*velint[0]
                              + velint[1]*velint[1]
                              + velint[2]*velint[2]);
    
    // normed velocity at element centre
    vector<double>  velino(NSD_);     
    if(vel_norm>=EPS6)
    {
        velino[0] = velint[0]/vel_norm;
        velino[1] = velint[1]/vel_norm;
        velino[2] = velint[2]/vel_norm;
    }
    else
    {
        velino[0] = 1.0;
        velino[1] = 0.0;
        velino[2] = 0.0;
    }
    
    // get streamlength
    double val = 0.0;
    for (int inode=0;inode<numnode;inode++)
    {
        val += abs(velino[0]*derxy(0,inode) 
                  +velino[1]*derxy(1,inode) 
                  +velino[2]*derxy(2,inode));
    }
    const double strle = 2.0/val;

    // calculate tau
    vector<double>  tau(3); // stab parameters
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


        const double re1 =/* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
        const double re2 = mk * vel_norm * strle / /* *1.0 */(2.0 * visc);    /* convective : viscous forces */

        const double xi1 = DMAX(re1,1.0);
        const double xi2 = DMAX(re2,1.0);

        tau[0] = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

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
        tau[1] = DSQR(hk) / (DSQR(hk) * xi_viscous + (/* 2.0*/ 4.0 * timefac * visc/mk) * xi_convect);
    
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
        tau[2] = vel_norm * hk * 0.5 * xi_tau_c /timefac;
      
    }
    else
    {// stabilization parameters for stationary case
    
        // compute tau_Mu    
        const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
        const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
        tau[0] = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);
 
        // compute tau_Mp
        const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
        const double xi_tau_mp = DMAX(re_tau_mp,1.0);
        tau[1] = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);    

        // compute tau_C
        const double xi_tau_c = DMIN(re_tau_mp, 1.0);
        tau[2] = 0.5*vel_norm*hk*xi_tau_c;
    }
    return tau;
}

//
// calculate Jacobian matrix
//
Epetra_SerialDenseMatrix DRT::Elements::Fluid3::getJacobiMatrix(
                      const Epetra_SerialDenseMatrix& xyze,
                      const Epetra_SerialDenseMatrix& deriv,
                      const int                       iel) const
{
  Epetra_SerialDenseMatrix    xjm(NSD_,NSD_);
    
  // determine jacobian matrix at point r,s,t
  for (int isd=0; isd<NSD_; isd++)
  {
    for (int jsd=0; jsd<NSD_; jsd++)
      {
      double dum = 0.0;
      for (int inode=0; inode<iel; inode++)
      {
        dum += deriv(isd,inode)*xyze(jsd,inode);
      }
      xjm(isd,jsd) = dum;
    }
  }

  return xjm;
}


double DRT::Elements::Fluid3::getDeterminante(const Epetra_SerialDenseMatrix&  xjm) const
{
    // determinant of jacobian matrix
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);

    if(det < 0.0)
    {
        printf("\n");
        printf("GLOBAL ELEMENT NO.%i\n",Id());
        printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n", det);
        dserror("Stopped not regulary!\n");
    }

    return det;
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::Elements::Fluid3::f3_getbodyforce(
        const double          time,
        const ParameterList&  params
)
{
  const int iel = NumNode();
  Epetra_SerialDenseMatrix edeadng(NSD_,iel);
    
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for(int inode=0;inode<iel;inode++)
  {
    Nodes()[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
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
      if(time >= 0)
      {
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for(int jnode=0;jnode<iel;jnode++)
    {
      Nodes()[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<NSD_;isd++)
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
      for(int isd=0;isd<3;isd++)
      {
        edeadng(isd,inode)=0.0;
      }
    }
  }

  return edeadng;
}

//
//  calculate global derivatives w.r.t. x,y,z at point r,s,t
//
void DRT::Elements::Fluid3::f3_gder(
  Epetra_SerialDenseMatrix& derxy,
	const Epetra_SerialDenseMatrix& deriv,
  const Epetra_SerialDenseMatrix& xjm,
	const double& det,
  const int iel
	)
{
  Epetra_SerialDenseMatrix 	xji(3,3);  // inverse of jacobian matrix

  // initialistion
  for(int k=0;k<iel;k++)
  {
    derxy(0,k) = 0.0;
    derxy(1,k) = 0.0;
    derxy(2,k) = 0.0;
  }

	// inverse of jacobian
  xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
  xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
  xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
  xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
  xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
  xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
  xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
  xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
  xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

	// calculate global derivatives
  for (int inode=0;inode<iel;inode++)
  {
    derxy(0,inode) +=   xji(0,0) * deriv(0,inode)
                      + xji(0,1) * deriv(1,inode)
                      + xji(0,2) * deriv(2,inode) ;
    derxy(1,inode) +=   xji(1,0) * deriv(0,inode)
                      + xji(1,1) * deriv(1,inode)
                      + xji(1,2) * deriv(2,inode) ;
    derxy(2,inode) +=   xji(2,0) * deriv(0,inode)
                      + xji(2,1) * deriv(1,inode)
                      + xji(2,2) * deriv(2,inode) ;
  }
  return;
} // end of DRT:Elements:Fluid3:f3_gder


void DRT::Elements::Fluid3::f3_gder2(
  const Epetra_SerialDenseMatrix& xyze,
  const Epetra_SerialDenseMatrix& xjm,
  const Epetra_SerialDenseMatrix& derxy,
  Epetra_SerialDenseMatrix& derxy2,
  const Epetra_SerialDenseMatrix& deriv2,
  const int iel
	)
{
  // initialize and zero out everything
  Epetra_SerialDenseMatrix bm(6,6);
  Epetra_SerialDenseMatrix xder2(6,3);

  // calculate elements of jacobian_bar matrix
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

  // inverse of jacobian_bar matrix

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

  // initialise
  /*   already initialized by constructor of EpetraSerialDenseMeatrix
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<6;j++) 
      xder2(j,i)=0.0;
  }
  */

  for (int i=0;i<iel;i++)
  {
    for (int j=0;j<6;j++) derxy2(j,i)=0.0;
  }

  // determine 2nd derivatives of coord.-functions
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
  }

  // calculate second global derivatives
  for (int inode=0;inode<iel;inode++)
  {
    const double r0 = deriv2(0,inode) - xder2(0,0)*derxy(0,inode) 
                                      - xder2(0,1)*derxy(1,inode)
                                      - xder2(0,2)*derxy(2,inode);
    const double r1 = deriv2(1,inode) - xder2(1,0)*derxy(0,inode) 
                                      - xder2(1,1)*derxy(1,inode)
                                      - xder2(1,2)*derxy(2,inode);
    const double r2 = deriv2(2,inode) - xder2(2,0)*derxy(0,inode) 
                                      - xder2(2,1)*derxy(1,inode)
                                      - xder2(2,2)*derxy(2,inode);
    const double r3 = deriv2(3,inode) - xder2(3,0)*derxy(0,inode)
                                      - xder2(3,1)*derxy(1,inode)
                                      - xder2(3,2)*derxy(2,inode);
    const double r4 = deriv2(4,inode) - xder2(4,0)*derxy(0,inode)
                                      - xder2(4,1)*derxy(1,inode)
                                      - xder2(4,2)*derxy(2,inode);
    const double r5 = deriv2(5,inode) - xder2(5,0)*derxy(0,inode)
                                      - xder2(5,1)*derxy(1,inode)
                                      - xder2(5,2)*derxy(2,inode);

    derxy2(0,inode) += bm(0,0)*r0 + bm(0,1)*r1 + bm(0,2)*r2
                    +  bm(0,3)*r3 + bm(0,4)*r4 + bm(0,5)*r5;
    derxy2(1,inode) += bm(1,0)*r0 + bm(1,1)*r1 + bm(1,2)*r2
                    +  bm(1,3)*r3 + bm(1,4)*r4 + bm(1,5)*r5;
    derxy2(2,inode) += bm(2,0)*r0 + bm(2,1)*r1 + bm(2,2)*r2
                    +  bm(2,3)*r3 + bm(2,4)*r4 + bm(2,5)*r5;
    derxy2(3,inode) += bm(3,0)*r0 + bm(3,1)*r1 + bm(3,2)*r2
                    +  bm(3,3)*r3 + bm(3,4)*r4 + bm(3,5)*r5;
    derxy2(4,inode) += bm(4,0)*r0 + bm(4,1)*r1 + bm(4,2)*r2
                    +  bm(4,3)*r3 + bm(4,4)*r4 + bm(4,5)*r5;
    derxy2(5,inode) += bm(5,0)*r0 + bm(5,1)*r1 + bm(5,2)*r2
                    +  bm(5,3)*r3 + bm(5,4)*r4 + bm(5,5)*r5;
  }

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

void DRT::Elements::Fluid3::f3_calmat(
       Epetra_SerialDenseMatrix&        estif,
       Epetra_SerialDenseVector&        eforce,
       const vector<double>&            velint,
       const vector<double>&            histvec,
       const vector<double>&            gridvint,
       const double&   	                press,
       const Epetra_SerialDenseMatrix&  vderxy,
       const Epetra_SerialDenseMatrix&  vderxy2,
       const vector<double>&            gradp,
       const Epetra_SerialDenseVector&  funct,
       const vector<double>&            tau,
       const Epetra_SerialDenseMatrix&  derxy,
       const Epetra_SerialDenseMatrix&  derxy2,
       const vector<double>&            edeadng,
       const double&                    fac,
       const double&                    visc,
       const int&                       iel,
       ParameterList&                   params
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




Epetra_SerialDenseMatrix  vconv_r(3,iel);

/*========================== initialisation ============================*/
// One-step-Theta: timefac = theta*dt
// BDF2:           timefac = 2/3 * dt
const double timefac = params.get<double>("time constant for integration",-1.0);
  if (timefac < 0.0) dserror("No time constant for integration supplied");

// time step size
//double dt = params.get<double>("delta time",-1.0);
//  if (dt == -1.0) dserror("No dta supplied");

// stabilisation parameter
const double tau_M  = tau[0]*fac;
const double tau_Mp = tau[1]*fac;
const double tau_C  = tau[2]*fac;

// integration factors and coefficients of single terms
// double time2nue   = timefac * 2.0 * visc;
const double timetauM   = timefac * tau_M;
const double timetauMp  = timefac * tau_Mp;

const double ttimetauM  = timefac * timetauM;
const double ttimetauMp = timefac * timetauMp;
const double timefacfac = timefac * fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
vector<double>  rhsint(3);          /* total right hand side terms at int.-point       */
rhsint[0] = histvec[0] + edeadng[0]*timefac;
rhsint[1] = histvec[1] + edeadng[1]*timefac;
rhsint[2] = histvec[2] + edeadng[2]*timefac;
/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
vector<double>  conv_old(3);
conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1]
            + vderxy(0,2) * velint[2];
conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1]
            + vderxy(1,2) * velint[2];
conv_old[2] = vderxy(2,0) * velint[0] + vderxy(2,1) * velint[1]
            + vderxy(2,2) * velint[2];

/* new for incremental formulation: */
/* Convective term  u_G_old * grad u_old: */
vector<double>  conv_g_old(3);
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
vector<double>  visc_old(3);
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


void DRT::Elements::Fluid3::f3_calmat_stationary(
       Epetra_SerialDenseMatrix&        estif,
       Epetra_SerialDenseVector&        eforce,
       const vector<double>&            velint,
       const vector<double>&            histvec,
       const vector<double>&            gridvint,
       const double&   	                press,
       const Epetra_SerialDenseMatrix&  vderxy,
       const Epetra_SerialDenseMatrix&  vderxy2,
       const vector<double>&            gradp,
       const Epetra_SerialDenseVector&  funct,
       const vector<double>&            tau,
       const Epetra_SerialDenseMatrix&  derxy,
       const Epetra_SerialDenseMatrix&  derxy2,
       const vector<double>&            edeadng,
       const double&                    fac,
       const double&                    visc,
       const int&                       iel,
       ParameterList& 	                params
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
    Epetra_SerialDenseVector&  funct,
    Epetra_SerialDenseMatrix&  derxy,
    Epetra_SerialDenseMatrix&  derxy2,
    vector<double>&  	       tau,
    const double&                    fac,
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
    Epetra_SerialDenseVector&  funct,
    Epetra_SerialDenseMatrix&  derxy,
    Epetra_SerialDenseMatrix&  derxy2,
    vector<double>&  	       tau,
    const double&                    fac,
    const double&              visc,
    const int&                 iel,
    ParameterList& 	       params)
{

  // set parameters
  const double tauM   = tau[0];
  const double tauMp  = tau[1];
  const double tauC   = tau[2];

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


// get optimal gaussrule for discretization type
GaussRule3D DRT::Elements::Fluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_10point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


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
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	deriv(3,iel);

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

  const double a      = PI/4.0;
  const double d      = PI/2.0;

  /* get viscosity ---*/
  const double  visc = material->m.fluid->viscosity;

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  double         p;
  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
          /*------------------- declaration of gauss point variables ---*/
        const double e1 = intpoints.qxg[iquad][0];
        const double e2 = intpoints.qxg[iquad][1];
        const double e3 = intpoints.qxg[iquad][2];
        shape_function_3D(funct,e1,e2,e3,distype);
        shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

        // get Jacobian matrix and determinant
        const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
        const double det = getDeterminante(xjm);
        const double fac = intpoints.qwgt[iquad]*det;

        /*---------------------- get velocity sol at integration point */
        for (int i=0;i<3;i++)
        {
          velint[i]=0.0;
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
          xint[i]=0.0;
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

  } // end of loop over integration points


  // we use the parameterlist as a container to transport the calculated
  // errors from the elements to the dynamic routine

  params.set<double>("L2 integrated velocity error",velerr);
  params.set<double>("L2 integrated pressure error",preerr);

  return;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::Elements::Fluid3::is_higher_order_element(
              const DRT::Element::DiscretizationType  distype) const
{
    bool hoel = true;
    switch (distype)
    {
    case hex8: case hex20: case hex27: case tet10:
        hoel = true;
        break;
    case tet4:
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
int DRT::Elements::Fluid3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
