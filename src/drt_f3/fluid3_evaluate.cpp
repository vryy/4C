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
#include "fluid3_impl.H"
#include "fluid3_genalpha_resVMM.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"

#include <blitz/array.h>
#include <Epetra_SerialDenseSolver.h>

using namespace DRT::Utils;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


DRT::Elements::Fluid3Impl* DRT::Elements::Fluid3::Impl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3Impl* f;
    if (f==NULL)
      f = new Fluid3Impl(8);
    return f;
  }
  case 20:
  {
    static Fluid3Impl* f;
    if (f==NULL)
      f = new Fluid3Impl(20);
    return f;
  }
  case 27:
  {
    static Fluid3Impl* f;
    if (f==NULL)
      f = new Fluid3Impl(27);
    return f;
  }
  case 4:
  {
    static Fluid3Impl* f;
    if (f==NULL)
      f = new Fluid3Impl(4);
    return f;
  }
  case 10:
  {
    static Fluid3Impl* f;
    if (f==NULL)
      f = new Fluid3Impl(10);
    return f;
  }
  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

DRT::Elements::Fluid3GenalphaResVMM* DRT::Elements::Fluid3::GenalphaResVMM()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3GenalphaResVMM* f;
    if (f==NULL)
      f = new Fluid3GenalphaResVMM(8);
    return f;
  }
  case 20:
  {
    static Fluid3GenalphaResVMM* f;
    if (f==NULL)
      f = new Fluid3GenalphaResVMM(20);
    return f;
  }
  case 27:
  {
    static Fluid3GenalphaResVMM* f;
    if (f==NULL)
      f = new Fluid3GenalphaResVMM(27);
    return f;
  }
  case 4:
  {
    static Fluid3GenalphaResVMM* f;
    if (f==NULL)
      f = new Fluid3GenalphaResVMM(4);
    return f;
  }
  case 10:
  {
    static Fluid3GenalphaResVMM* f;
    if (f==NULL)
      f = new Fluid3GenalphaResVMM(10);
    return f;
  }
  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

// converts a string into an Action for this element
DRT::Elements::Fluid3::ActionType DRT::Elements::Fluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::Elements::Fluid3::ActionType act = Fluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = Fluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
    act = Fluid3::calc_fluid_genalpha_sysmat_and_residual;
  else if (action == "time update for subscales")
    act = Fluid3::calc_fluid_genalpha_update_for_subscales;
  else if (action == "calc_fluid_beltrami_error")
    act = Fluid3::calc_fluid_beltrami_error;
  else if (action == "calc_turbulence_statistics")
    act = Fluid3::calc_turbulence_statistics;
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
 // converts a string into an stabilisation action for this element
 //                                                          gammi 08/07
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3::StabilisationAction DRT::Elements::Fluid3::ConvertStringToStabAction(
  const string& action) const
{
  DRT::Elements::Fluid3::StabilisationAction act = stabaction_unspecified;
  
  map<string,StabilisationAction>::iterator iter=stabstrtoact_.find(action);

  if (iter != stabstrtoact_.end())
  {
    act = (*iter).second;
  }
  else
  {
    char errorout[200];
    sprintf(errorout,"looking for stab action (%s) not contained in map",action.c_str());
    dserror(errorout);
  }
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
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p at time n+1 (trial)");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("old solution data for rhs");
        if (velnp==null || hist==null)
          dserror("Cannot get state vectors 'velnp' and/or 'hist'");

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::Utils::ExtractMyValues(*velnp,myvelnp,lm);
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

        // split velocity and pressure
        // create blitz objects

        const int numnode = NumNode();
        blitz::Array<double, 1> eprenp(numnode);
        blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> edispnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> egridv(3,numnode,blitz::ColumnMajorArray<2>());

        for (int i=0;i<numnode;++i)
        {
          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          eprenp(i) = myvelnp[3+(i*4)];

          evhist(0,i) = myhist[0+(i*4)];
          evhist(1,i) = myhist[1+(i*4)];
          evhist(2,i) = myhist[2+(i*4)];
        }

        if (is_ale_)
        {
          for (int i=0;i<numnode;++i)
          {
            edispnp(0,i) = mydispnp[0+(i*4)];
            edispnp(1,i) = mydispnp[1+(i*4)];
            edispnp(2,i) = mydispnp[2+(i*4)];

            egridv(0,i) = mygridv[0+(i*4)];
            egridv(1,i) = mygridv[1+(i*4)];
            egridv(2,i) = mygridv[2+(i*4)];
          }
        }

        // get control parameter
        const bool is_stationary = params.get<bool>("using stationary formulation",false);
        const double time = params.get<double>("total time",-1.0);

        bool newton              = params.get<bool>("include reactive terms for linearisation",false);
        bool pstab  =true;
        bool supg   =true;
        bool vstab  =true;
        bool cstab  =true;
        
        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        double timefac = 0;
        if (not is_stationary)
        {
          timefac = params.get<double>("time constant for integration",-1.0);
          if (timefac < 0.0) dserror("No time constant for integration supplied");
        }

        // wrap epetra serial dense objects in blitz objects
        blitz::Array<double, 2> estif(elemat1.A(),
                                      blitz::shape(elemat1.M(),elemat1.N()),
                                      blitz::neverDeleteData,
                                      blitz::ColumnMajorArray<2>());
        blitz::Array<double, 1> eforce(elevec1.Values(),
                                       blitz::shape(elevec1.Length()),
                                       blitz::neverDeleteData);

        // calculate element coefficient matrix and rhs
        Impl()->Sysmat(this,
                       evelnp,
                       eprenp,
                       evhist,
                       edispnp,
                       egridv,
                       estif,
                       eforce,
                       actmat,
                       time,
                       timefac,
                       newton ,
                       pstab  ,
                       supg   ,
                       vstab  ,
                       cstab  ,
                       is_stationary);

        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        params.set("density", actmat->m.fluid->density);

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
      case calc_turbulence_statistics:
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocities and pressure from the global distributed vectors

          // velocity and pressure values (n+1)
          RefCountPtr<const Epetra_Vector> velnp
            = discretization.GetState("u and p (n+1,converged)");


          if (velnp==null)
          {
            dserror("Cannot get state vectors 'velnp'");
          }

          // extract local values from the global vectors
          vector<double> mysol  (lm.size());
          DRT::Utils::ExtractMyValues(*velnp,mysol,lm);

          // integrate mean values
          f3_calc_means(mysol,params);

        }
      }
      break;
      case calc_fluid_genalpha_sysmat_and_residual:
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
        vector<double> myvelnp(lm.size());
        DRT::Utils::ExtractMyValues(*velnp,myvelnp,lm);

        vector<double> myvelaf(lm.size());
        DRT::Utils::ExtractMyValues(*velaf,myvelaf,lm);

        vector<double> myaccam(lm.size());
        DRT::Utils::ExtractMyValues(*accam,myaccam,lm);

        // create blitz matrix objects

        const int numnode = NumNode();
        blitz::Array<double, 1> eprenp (  numnode);
        blitz::Array<double, 2> evelnp (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evelaf (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> eaccam (3,numnode,blitz::ColumnMajorArray<2>());

        
        // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
        // Additionally only the 'velocity' components of my_velaf
        // and my_accam are important!
        for (int i=0;i<numnode;++i)
        {
          eprenp(i)   = myvelnp[3+(i*4)];

          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          evelaf(0,i) = myvelaf[0+(i*4)];
          evelaf(1,i) = myvelaf[1+(i*4)];
          evelaf(2,i) = myvelaf[2+(i*4)];

          eaccam(0,i) = myaccam[0+(i*4)];
          eaccam(1,i) = myaccam[1+(i*4)];
          eaccam(2,i) = myaccam[2+(i*4)];
        }

        // --------------------------------------------------
        // set parameters for time integration
        ParameterList& timelist = params.sublist("time integration parameters");
        
        const double alphaM = timelist.get<double>("alpha_M");
        const double alphaF = timelist.get<double>("alpha_F");
        const double gamma  = timelist.get<double>("gamma");
        const double dt     = timelist.get<double>("dt");
        const double time   = timelist.get<double>("time");

        // --------------------------------------------------
        // set parameters for nonlinear treatment

        const bool newton = params.get<bool>("include reactive terms for linearisation");
   
        // --------------------------------------------------
        // set parameters for stabilisation
        ParameterList& stablist = params.sublist("stabilisation");

        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["subscales quasistatic"                        ]=subscales_quasistatic;
          stabstrtoact_["subscales time dependent"                     ]=subscales_time_dependent;
          stabstrtoact_["drop inertia stabilisation"                   ]=inertia_stab_drop;
          stabstrtoact_["keep inertia stabilisation"                   ]=inertia_stab_keep;
          stabstrtoact_["assume inf-sup stable"                        ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["use pspg stabilisation"                       ]=pstab_use_pspg;
          stabstrtoact_["no convective stabilisation"                  ]=convective_stab_none;
          stabstrtoact_["supg convective stabilisation"                ]=convective_stab_supg;
          stabstrtoact_["no viscous stabilisation"                     ]=viscous_stab_none;
          stabstrtoact_["viscous stabilisation of gls type"            ]=viscous_stab_gls;
          stabstrtoact_["viscous stabilisation of gls type (only rhs)" ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["viscous stabilisation of agls type"           ]=viscous_stab_agls;
          stabstrtoact_["viscous stabilisation of agls type (only rhs)"]=viscous_stab_agls_only_rhs;
          stabstrtoact_["use continuity stabilisation"                 ]=continuity_stab_yes;
          stabstrtoact_["no continuity stabilisation"                  ]=continuity_stab_none;
          stabstrtoact_["cross stress stabilisation (on rhs)"          ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no cross stress stabilisation"                ]=cross_stress_stab_none;
          stabstrtoact_["reynolds stress stabilisation (on rhs)"       ]=reynolds_stress_stab_only_rhs;
          stabstrtoact_["no reynolds stress stabilisation"             ]=reynolds_stress_stab_none;
        }
        
        StabilisationAction tds      = ConvertStringToStabAction(stablist.get<string>("time tracking of subscales"));
        StabilisationAction inertia  = ConvertStringToStabAction(stablist.get<string>("use subscale acceleration term in weak form"));
        StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("stabilisation (saddle point problem)"));
        StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("convective stabilisation"));
        StabilisationAction agls     = ConvertStringToStabAction(stablist.get<string>("viscous stabilisation"));
        StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("continuity stabilisation"));
        StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("cross stress stabilisation"));
        StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("reynolds stress stabilisation"));

        // --------------------------------------------------
        // specify what to compute
        const bool compute_elemat = params.get<bool>("compute element matrix");
        
        
        // --------------------------------------------------
        // calculate element coefficient matrix
        GenalphaResVMM()->Sysmat(this,
                                 elemat1,
                                 elevec1,
                                 evelnp, 
                                 eprenp, 
                                 eaccam, 
                                 evelaf, 
                                 actmat,
                                 alphaM, 
                                 alphaF, 
                                 gamma,  
                                 dt,     
                                 time,   
                                 newton, 
                                 tds,    
                                 inertia,
                                 pspg,   
                                 supg,   
                                 agls,   
                                 cstab,   
                                 cross,  
                                 reynolds,
                                 compute_elemat
          );
      }
      break;
      case calc_fluid_genalpha_update_for_subscales:
        // most recent subscale pressure becomes the old subscale pressure
        // for the next timestep
        //
        //  ~n   ~n+1
        //  p <- p
        //
        sub_pre_old_ = sub_pre_;
        
        // the old subscale acceleration for the next timestep is calculated
        // on the fly, not stored on the element
        /*
                     ~n+1   ~n
             ~ n     u    - u     ~ n   / gamma-1.0 \
            acc  <-  --------- - acc * |  ---------  |
                     gamma*dt           \   gamma   /
        */
        {
          const double dt     = params.get<double>("dt");
          const double gamma  = params.get<double>("gamma");

          sub_acc_old_ = (sub_vel_-sub_vel_old_)/(gamma*dt)
                         -
                         sub_acc_old_*(gamma-1.0)/gamma;
        }
        // most recent subscale velocity becomes the old subscale velocity
        // for the next timestep
        //
        //  ~n   ~n+1
        //  u <- u
        //
        sub_vel_old_=sub_vel_;
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
  const int NSD = 3;

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
  const bool higher_order_ele = isHigherOrderElement(distype);

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
    for (int isd=0;isd<NSD;isd++)
    {
      velint[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        velint[isd] += funct[inode]*evelnp[isd+(3*inode)];
      }
    }

    // get history data (n,i) at integration point
    //expression for f3_veci(histvec,funct,evhist,iel);
    for (int isd=0;isd<NSD;isd++)
    {
      histvec[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        histvec[isd] += funct[inode]*evhist[isd+(3*inode)];
      }
    }

    // get velocity (np,i) derivatives at integration point
    // expression for f3_vder(vderxy,derxy,evelnp,iel);
    for (int isd=0;isd<NSD;isd++)
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
    vector<double>    gridvelint(NSD);
    if (is_ale_)
    {
      for (int isd=0; isd<NSD; isd++)
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
    vector<double>    gradp(NSD);
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
    for (int isd=0;isd<NSD;isd++)
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
    const int NSD = 3;
    const int iel = NumNode();
    Epetra_SerialDenseMatrix xyze(NSD,iel);

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
    const int NSD = 3;
    // use one point gauss rule to calculate tau at element center
    GaussRule3D integrationrule_stabili = intrule_hex_1point;
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
    Epetra_SerialDenseMatrix    deriv(NSD, numnode);
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
    vector<double>  velint(NSD);
    for (int isd=0;isd<NSD;isd++)
    {
        velint[isd]=0.0;
        for (int inode=0;inode<numnode;inode++)
        {
            velint[isd] += funct[inode]*evelnp[isd+(NSD*inode)];
        }
    }

    // get Jacobian matrix and determinant
    const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,numnode);
    const double det = getDeterminante(xjm);
    const double vol = wquad*det;

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = pow((SIX*vol/PI),(1.0/3.0))/sqrt(3.0);

    // get derivatives
    Epetra_SerialDenseMatrix    derxy(NSD, numnode);
    f3_gder(derxy,deriv,xjm,det,numnode);

    // get velocity norm
    const double vel_norm=sqrt( velint[0]*velint[0]
                              + velint[1]*velint[1]
                              + velint[2]*velint[2]);

    // normed velocity at element centre
    vector<double>  velino(NSD);
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

/*----------------------------------------------------------------------*
 | calculate Jacobian matrix and it's determinant (private) gammi  07/07|
 | Well, I think we actually compute its transpose....
 |
 |     +-            -+ T      +-            -+
 |     | dx   dx   dx |        | dx   dy   dz |
 |     | --   --   -- |        | --   --   -- |
 |     | dr   ds   dt |        | dr   dr   dr |
 |     |              |        |              |
 |     | dy   dy   dy |        | dx   dy   dz |
 |     | --   --   -- |   =    | --   --   -- |
 |     | dr   ds   dt |        | ds   ds   ds |
 |     |              |        |              |
 |     | dz   dz   dz |        | dx   dy   dz |
 |     | --   --   -- |        | --   --   -- |
 |     | dr   ds   dt |        | dt   dt   dt |
 |     +-            -+        +-            -+
 |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::Elements::Fluid3::getJacobiMatrix(
                      const Epetra_SerialDenseMatrix& xyze,
                      const Epetra_SerialDenseMatrix& deriv,
                      const int                       iel) const
{
  const int NSD = 3;
  Epetra_SerialDenseMatrix    xjm(NSD,NSD);

  // determine jacobian matrix at point r,s,t
  for (int isd=0; isd<NSD; isd++)
  {
    for (int jsd=0; jsd<NSD; jsd++)
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
  const int NSD = 3;
  const int iel = NumNode();
  Epetra_SerialDenseMatrix edeadng(NSD,iel);

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
      if(time >= 0.0)
      {
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
         dserror("Negative time value in body force calculation: time = %f",time);
        //curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
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

      for(int isd=0;isd<NSD;isd++)
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
  Epetra_SerialDenseMatrix chainrulerhs(6,iel);

  // calculate elements of jacobian_bar matrix
  bm(0,0) = xjm(0,0)*xjm(0,0);
  bm(1,0) = xjm(1,0)*xjm(1,0);
  bm(2,0) = xjm(2,0)*xjm(2,0);
  bm(3,0) = xjm(0,0)*xjm(1,0);
  bm(4,0) = xjm(0,0)*xjm(2,0);
  bm(5,0) = xjm(2,0)*xjm(1,0);

  bm(0,1) = xjm(0,1)*xjm(0,1);
  bm(1,1) = xjm(1,1)*xjm(1,1);
  bm(2,1) = xjm(2,1)*xjm(2,1);
  bm(3,1) = xjm(0,1)*xjm(1,1);
  bm(4,1) = xjm(0,1)*xjm(2,1);
  bm(5,1) = xjm(2,1)*xjm(1,1);

  bm(0,2) = xjm(0,2)*xjm(0,2);
  bm(1,2) = xjm(1,2)*xjm(1,2);
  bm(2,2) = xjm(2,2)*xjm(2,2);
  bm(3,2) = xjm(0,2)*xjm(1,2);
  bm(4,2) = xjm(0,2)*xjm(2,2);
  bm(5,2) = xjm(2,2)*xjm(1,2);

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

  //init sol to zero
  memset(derxy2.A(),0,derxy2.M()*derxy2.N()*sizeof(double));


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
  xder2.Multiply('N','T',1.0,deriv2,xyze,0.0);

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

  xder2.Multiply(false,derxy,chainrulerhs);
  chainrulerhs.Scale(-1.0);

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

  chainrulerhs+=deriv2;

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

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix (bm);
  solver.SetVectors(derxy2,chainrulerhs);
  solver.Solve();

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

// get optimal gaussrule for discretization type
GaussRule3D DRT::Elements::Fluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule_hex_8point;
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
 |  calculate error for beltrami test problem               gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_int_beltrami_err(
  vector<double>&           evelnp,
  vector<double>&           eprenp,
  struct _MATERIAL*         material,
  ParameterList& 	    params
  )
{
  const int NSD = 3;

  // add element error to "integrated" error
  double velerr = params.get<double>("L2 integrated velocity error");
  double preerr = params.get<double>("L2 integrated pressure error");

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	deriv(3,iel);

  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // set constants for analytical solution
  const double t = params.get("total time",-1.0);
  dsassert (t >= 0.0, "beltrami: no total time for error calculation");

  const double a      = PI/4.0;
  const double d      = PI/2.0;

  // get viscosity
  const double  visc = material->m.fluid->viscosity;

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // declaration of gauss point variables
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    // get Jacobian matrix and determinant
    const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
    const double det = getDeterminante(xjm);
    const double fac = intpoints.qwgt[iquad]*det;

    // get velocity sol at integration point
    for (int i=0;i<3;i++)
    {
      velint[i]=0.0;
      for (int j=0;j<iel;j++)
      {
        velint[i] += funct[j]*evelnp[i+(3*j)];
      }
    }

    // get pressure sol at integration point
    preint = 0;
    for (int inode=0;inode<iel;inode++)
    {
      preint += funct[inode]*eprenp[inode];
    }

    // get velocity sol at integration point
    for (int isd=0;isd<3;isd++)
    {
      xint[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        xint[isd] += funct[inode]*xyze(isd,inode);
      }
    }

    // compute analytical pressure
    const double p = -a*a/2.0 *
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

    for (int isd=0;isd<NSD;isd++)
    {
      deltavel[isd] = velint[isd]-u[isd];
    }

    // add square to L2 error
    for (int isd=0;isd<NSD;isd++)
    {
      velerr += deltavel[isd]*deltavel[isd]*fac;
    }
    preerr += deltap*deltap*fac;

  } // end of loop over integration points


  // we use the parameterlist as a container to transport the calculated
  // errors from the elements to the dynamic routine

  params.set<double>("L2 integrated velocity error",velerr);
  params.set<double>("L2 integrated pressure error",preerr);

  return;
}

/*---------------------------------------------------------------------*
 | Calculate spatial mean values for channel flow (cartesian mesh)
 |                                                           gammi 07/07
 |
 | The necessary element integration is performed in here. The element
 | is cut into at least two (HEX8) or three (quadratic elements) planes,
 | the spatial functions (velocity, pressure etc.) are integrated over
 | this plane and this element contribution is added to a processor local
 | vector (see formulas below for a exact description of the output).
 | The method assumes, that all elements are of the same rectangular
 | shape in the "inplanedirection". In addition, it is assumed that
 | the sampling planes are distributed equidistant in the element.
 |
 |
 |                      ^ normdirect       integration plane
 |                      |                /
 |                      |               /
 |                      |
 |                +-----|-------------+
 |               /|     |            /|
 |              / |     |           / |
 |             /  |     |          /  |
 |            /   |     |         /   |
 |           /    +-----|--------/----+ ---- additional integration
 |          /    /|     |       /    /|      plane (for quadratic elements)
 |         /    / |     |      /    / |
 |        +-------------------+    /  |
 |        |   /   |     *-----|---+------------>
 |        |  /    +----/------|--/----+         inplanedirect[1]
 |        | /    /    /       | /    /
 |        |/    /    /        |/    /   \
 |        +---------+---------+    /     \
 |        |   /    /          |   /       integration plane
 |        |  /    /           |  /
 |        | /    /            | /
 |        |/    /             |/
 |        +----/--------------+
 |            /
 |           /   inplanedirect[0]
 |
 |
 |  Example for a mean value evaluation:
 |
 |         1.0       /
 |  _               |
 |  u = -------- *  | u(x,y,z) dx dy dz =
 |      +---        |
 |       \         / A
 |       / area
 |      +---
 |
 |
 |        1.0      /
 |                |            area
 |  =  -------- * | u(r,s,t) * ---- dr ds dt
 |     +---       |              4
 |      \        /  [-1:1]^2
 |      / area
 |     +---
 |
 |
 |
 |         1.0      /
 |                 |            1
 |  =   -------- * | u(r,s,t) * - dr ds dt
 |                 |            4
 |       numele   /  [-1:1]^2
 |
 |                |                        |
 |                +------------------------+
 |             this is the integral we compute!
 |
 | The factor 1/4 is necessary since we use a reference element of
 | size 2x2
 |
 | The method computes:
 |                      _             _             _             _
 |             numele * u  , numele * v  , numele * w  , numele * p
 |                      ___           ___           ___           ___
 |                       ^2            ^2            ^2            ^2
 | and         numele * u  , numele * v  , numele * w  , numele * p
 |
 | as well as numele.
 | All results are communicated vi the parameter list!
 |
 *---------------------------------------------------------------------*/
void DRT::Elements::Fluid3::f3_calc_means(
  vector<double>&           sol  ,
  ParameterList& 	    params
  )
{

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();


  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");


  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  RefCountPtr<vector<double> > planes = params.get<RefCountPtr<vector<double> > >("coordinate vector for hom. planes");

  // get the pointers to the solution vectors
  RefCountPtr<vector<double> > sumu   = params.get<RefCountPtr<vector<double> > >("mean velocities x direction");
  RefCountPtr<vector<double> > sumv   = params.get<RefCountPtr<vector<double> > >("mean velocities y direction");
  RefCountPtr<vector<double> > sumw   = params.get<RefCountPtr<vector<double> > >("mean velocities z direction");
  RefCountPtr<vector<double> > sump   = params.get<RefCountPtr<vector<double> > >("mean pressure");

  RefCountPtr<vector<double> > sumsqu = params.get<RefCountPtr<vector<double> > >("variance velocities x direction");
  RefCountPtr<vector<double> > sumsqv = params.get<RefCountPtr<vector<double> > >("variance velocities y direction");
  RefCountPtr<vector<double> > sumsqw = params.get<RefCountPtr<vector<double> > >("variance velocities z direction");
  RefCountPtr<vector<double> > sumsqp = params.get<RefCountPtr<vector<double> > >("variance pressure");


  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  double min = xyze(normdirect,0);
  double max = xyze(normdirect,0);

  // set maximum and minimum value in wall normal direction
  for(int inode=0;inode<iel;inode++)
  {
    if(min > xyze(normdirect,inode))
    {
      min=xyze(normdirect,inode);
    }
    if(max < xyze(normdirect,inode))
    {
      max=xyze(normdirect,inode);
    }
  }

  // determine the ids of the homogeneous planes intersecting this element
  set<int> planesinele;
  for(unsigned nplane=0;nplane<planes->size();++nplane)
  {
    // get all available wall normal coordinates
    for(int nn=0;nn<iel;++nn)
    {
      if (min-2e-9 < (*planes)[nplane] && max+2e-9 > (*planes)[nplane])
      {
        planesinele.insert(nplane);
      }
    }
  }

  // remove lowest layer from planesinele to avoid double calculations. This is not done
  // for the first level (index 0) --- if deleted, shift the first integration point in
  // wall normal direction
  // the shift depends on the number of sampling planes in the element
  double shift=0;

  // set the number of planes which cut the element
  const int numplanesinele = planesinele.size();

  if(*planesinele.begin() != 0)
  {
    // this is not an element of the lowest element layer
    planesinele.erase(planesinele.begin());

    shift=2.0/((double) numplanesinele - 1.0);
  }
  else
  {
    // this is an element of the lowest element layer. Increase the counter
    // in order to compute the total number of elements in one layer
    int* count = params.get<int*>("count processed elements");

    (*count)++;
  }

  // determine the orientation of the rst system compared to the xyz system
  int elenormdirect=-1;
  bool upsidedown =false;
  // the only thing of interest is how normdirect is oriented in the
  // element coordinate system
  if(xyze(normdirect,4)-xyze(normdirect,0)>2e-9)
  {
    // t aligned
    elenormdirect =2;
    cout << "upsidedown false" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)>2e-9)
  {
    // s aligned
    elenormdirect =1;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)>2e-9)
  {
    // r aligned
    elenormdirect =0;
  }
  else if(xyze(normdirect,4)-xyze(normdirect,0)<-2e-9)
  {
    cout << xyze(normdirect,4)-xyze(normdirect,0) << &endl;
    // -t aligned
    elenormdirect =2;
    upsidedown =true;
    cout << "upsidedown true" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)<-2e-9)
  {
    // -s aligned
    elenormdirect =1;
    upsidedown =true;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)<-2e-9)
  {
    // -r aligned
    elenormdirect =0;
    upsidedown =true;
  }
  else
  {
    dserror("cannot determine orientation of plane normal in local coordinate system of element");
  }
  vector<int> inplanedirect;
  {
    set <int> inplanedirectset;
    for(int i=0;i<3;++i)
    {
      inplanedirectset.insert(i);
    }
    inplanedirectset.erase(elenormdirect);

    for(set<int>::iterator id = inplanedirectset.begin();id!=inplanedirectset.end() ;++id)
    {
      inplanedirect.push_back(*id);
    }
  }

  // allocate vector for shapefunctions
  Epetra_SerialDenseVector  funct(iel);

  // get the quad9 gaussrule for the in plane integration
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(intrule_quad_9point);

  // a hex8 element has two levels, the hex20 and hex27 element have three layers to sample
  // (now we allow even more)
  double layershift=0;

  // loop all levels in element
  for(set<int>::iterator id = planesinele.begin();id!=planesinele.end() ;++id)
  {
    // reset temporary values
    double ubar=0;
    double vbar=0;
    double wbar=0;
    double pbar=0;

    double usqbar=0;
    double vsqbar=0;
    double wsqbar=0;
    double psqbar=0;

    // get the intgration point in wall normal direction
    double e[3];

    e[elenormdirect]=-1.0+shift+layershift;
    if(upsidedown)
    {
      e[elenormdirect]*=-1;
    }

    // start loop over integration points in layer
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
      // get the other gauss point coordinates
      for(int i=0;i<2;++i)
      {
        e[inplanedirect[i]]=intpoints.qxg[iquad][i];
      }

      // compute the shape function values
      shape_function_3D(funct,e[0],e[1],e[2],distype);

      // check whether this gausspoint is really inside the desired plane
      {
        double x[3];
        x[0]=0;
        x[1]=0;
        x[2]=0;
        for(int inode=0;inode<iel;inode++)
        {
          x[0]+=funct[inode]*xyze(0,inode);
          x[1]+=funct[inode]*xyze(1,inode);
          x[2]+=funct[inode]*xyze(2,inode);
        }

        if(abs(x[normdirect]-(*planes)[*id])>2e-9)
        {
          dserror("Mixing up element cut planes during integration");
        }
      }

      //interpolated values at gausspoints
      double ugp=0;
      double vgp=0;
      double wgp=0;
      double pgp=0;

      // we assume that every 2d element we are integrating here is of
      // rectangular shape and every element is of the same size.
      // 1/4 is necessary since we use a reference element of size 2x2
      // the factor fac is omitting the element area up to now
      double fac=0.25*intpoints.qwgt[iquad];

      for(int inode=0;inode<iel;inode++)
      {
        ugp += funct[inode]*sol[inode*4  ];
        vgp += funct[inode]*sol[inode*4+1];
        wgp += funct[inode]*sol[inode*4+2];
        pgp += funct[inode]*sol[inode*4+3];
      }

      // add contribution to integral
      ubar   += ugp*fac;
      vbar   += vgp*fac;
      wbar   += wgp*fac;
      pbar   += pgp*fac;

      usqbar += ugp*ugp*fac;
      vsqbar += vgp*vgp*fac;
      wsqbar += wgp*wgp*fac;
      psqbar += pgp*pgp*fac;
    } // end loop integration points

    // add increments from this layer to processor local vectors
    (*sumu  )[*id] += ubar;
    (*sumv  )[*id] += vbar;
    (*sumw  )[*id] += wbar;
    (*sump  )[*id] += pbar;

    (*sumsqu)[*id] += usqbar;
    (*sumsqv)[*id] += vsqbar;
    (*sumsqw)[*id] += wsqbar;
    (*sumsqp)[*id] += psqbar;

    // jump to the next layer in the element.
    // in case of an hex8 element, the two coordinates are -1 and 1(+2)
    // for quadratic elements with three sample planes, we have -1,0(+1),1(+2)

    layershift+=2.0/((double) numplanesinele - 1.0);
  }


  return;
}

//
// check for higher order derivatives for shape functions
//
bool DRT::Elements::Fluid3::isHigherOrderElement(
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
