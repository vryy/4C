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

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "fluid3.H"
#include "fluid3_impl.H"
#include "fluid3_genalpha_resVMM.H"
#include "fluid3_stationary.H"

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


/*
  Depending on the type of the algorithm (the implementation) and the
  element type (tet, hex etc.), the elements allocate common static
  arrays.

  That means that for example all hex8 fluid elements of the stationary
  implementation have a pointer f8 to the same 'implementation class'
  containing all the element arrays for eight noded elements, and all
  wedge15 fluid elements of the same problem have a pointer f15 to
  the 'implementation class' containing all the element arrays for the
  15 noded element.
  
  */

DRT::Elements::Fluid3Impl* DRT::Elements::Fluid3::Impl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3Impl* f8;
    if (f8==NULL)
      f8 = new Fluid3Impl(8);
    return f8;
  }
  case 20:
  {
    static Fluid3Impl* f20;
    if (f20==NULL)
      f20 = new Fluid3Impl(20);
    return f20;
  }
  case 27:
  {
    static Fluid3Impl* f27;
    if (f27==NULL)
      f27 = new Fluid3Impl(27);
    return f27;
  }
  case 4:
  {
    static Fluid3Impl* f4;
    if (f4==NULL)
      f4 = new Fluid3Impl(4);
    return f4;
  }
  case 10:
  {
    static Fluid3Impl* f10;
    if (f10==NULL)
      f10 = new Fluid3Impl(10);
    return f10;
  }
  case 6:
  {
    static Fluid3Impl* f6;
    if (f6==NULL)
      f6 = new Fluid3Impl(6);
    return f6;
  }
  case 15:
  {
    static Fluid3Impl* f15;
    if (f15==NULL)
      f15 = new Fluid3Impl(15);
    return f15;
  }
  case 5:
  {
    static Fluid3Impl* f5;
    if (f5==NULL)
      f5 = new Fluid3Impl(5);
    return f5;
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
    static Fluid3GenalphaResVMM* f8;
    if (f8==NULL)
      f8 = new Fluid3GenalphaResVMM(8);
    return f8;
  }
  case 20:
  {
    static Fluid3GenalphaResVMM* f20;
    if (f20==NULL)
      f20 = new Fluid3GenalphaResVMM(20);
    return f20;
  }
  case 27:
  {
    static Fluid3GenalphaResVMM* f27;
    if (f27==NULL)
      f27 = new Fluid3GenalphaResVMM(27);
    return f27;
  }
  case 4:
  {
    static Fluid3GenalphaResVMM* f4;
    if (f4==NULL)
      f4 = new Fluid3GenalphaResVMM(4);
    return f4;
  }
  case 10:
  {
    static Fluid3GenalphaResVMM* f10;
    if (f10==NULL)
      f10 = new Fluid3GenalphaResVMM(10);
    return f10;
  }
  case 6:
  {
    static Fluid3GenalphaResVMM* f6;
    if (f6==NULL)
      f6 = new Fluid3GenalphaResVMM(6);
    return f6;
  }
  case 15:
  {
    static Fluid3GenalphaResVMM* f15;
    if (f15==NULL)
      f15 = new Fluid3GenalphaResVMM(15);
    return f15;
  }
  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}


DRT::Elements::Fluid3Stationary* DRT::Elements::Fluid3::StationaryImpl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3Stationary* f8;
    if (f8==NULL)
      f8 = new Fluid3Stationary(8);
    return f8;
  }
  case 20:
  {
    static Fluid3Stationary* f20;
    if (f20==NULL)
      f20 = new Fluid3Stationary(20);
    return f20;
  }
  case 27:
  {
    static Fluid3Stationary* f27;
    if (f27==NULL)
      f27 = new Fluid3Stationary(27);
    return f27;
  }
  case 4:
  {
    static Fluid3Stationary* f4;
    if (f4==NULL)
      f4 = new Fluid3Stationary(4);
    return f4;
  }
  case 10:
  {
    static Fluid3Stationary* f10;
    if (f10==NULL)
      f10 = new Fluid3Stationary(10);
    return f10;
  }
  case 6:
  {
    static Fluid3Stationary* f6;
    if (f6==NULL)
      f6 = new Fluid3Stationary(6);
    return f6;
  }
  case 15:
  {
    static Fluid3Stationary* f15;
    if (f15==NULL)
      f15 = new Fluid3Stationary(15);
    return f15;
  }
  case 5:
  {
    static Fluid3Stationary* f5;
    if (f5==NULL)
      f5 = new Fluid3Stationary(5);
    return f5;
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
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Fluid3::calc_fluid_stationary_systemmat_and_residual;  
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
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("hist");
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

        // create blitz objects for element arrays
        const int numnode = NumNode();
        blitz::Array<double, 1> eprenp(numnode);
        blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> edispnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> egridv(3,numnode,blitz::ColumnMajorArray<2>());

        // split velocity and pressure, insert into element arrays
        for (int i=0;i<numnode;++i)
        {
          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          eprenp(i) = myvelnp[3+(i*4)];

          // the history vector contains the information of time step t_n (mass rhs!)
          evhist(0,i) = myhist[0+(i*4)];
          evhist(1,i) = myhist[1+(i*4)];
          evhist(2,i) = myhist[2+(i*4)];
        }

        if (is_ale_)
        {
          // assign grid velocity and grid displacement to element arrays
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
        const double time = params.get<double>("total time",-1.0);

        bool newton = params.get<bool>("include reactive terms for linearisation",false);

        // the stabilisation scheme is hardcoded up to now --- maybe it's worth taking
        // this into the input or to choose a standard implementation and drop all ifs
        // on the element level -- if so, I would recommend to drop vstab...
        
        bool pstab  = true;
        bool supg   = true;
        bool vstab  = true;
        bool cstab  = true;

        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        double timefac = 0.0;
        timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");

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
                       cstab  );

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
          stabstrtoact_["quasistatic subscales"                        ]=subscales_quasistatic;
          stabstrtoact_["time dependent subscales"                     ]=subscales_time_dependent;
          stabstrtoact_["drop inertia stabilisation"                   ]=inertia_stab_drop;
          stabstrtoact_["keep inertia stabilisation"                   ]=inertia_stab_keep;
          stabstrtoact_["inf-sup-stable (off)"                         ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["(svel,nabla q)"                               ]=pstab_use_pspg;
          stabstrtoact_["off"                                          ]=convective_stab_none;
          stabstrtoact_["(svel,(u o nabla)v)"                          ]=convective_stab_supg;
          stabstrtoact_["off"                                          ]=viscous_stab_none;
          stabstrtoact_["(svel,+2 visc nabla o eps(v))"                ]=viscous_stab_gls;
          stabstrtoact_["(svel,+2 visc nabla o eps(v)) [RHS]"          ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["(svel,-2 visc nabla o eps(v))"                ]=viscous_stab_agls;
          stabstrtoact_["(svel,-2 visc nabla o eps(v)) [RHS]"          ]=viscous_stab_agls_only_rhs;
          stabstrtoact_["(spres,nabla o v)"                            ]=continuity_stab_yes;
          stabstrtoact_["off"                                          ]=continuity_stab_none;
          stabstrtoact_["((svel o nabla)u,v)"                          ]=cross_stress_stab;
          stabstrtoact_["((svel o nabla)u,v) [RHS]"                    ]=cross_stress_stab_only_rhs;
          stabstrtoact_["off"                                          ]=cross_stress_stab_none;
          stabstrtoact_["(svel,(svel o nabla)v) [RHS]"                 ]=reynolds_stress_stab_only_rhs;
          stabstrtoact_["off"                                          ]=reynolds_stress_stab_none;
        }

        StabilisationAction tds      = ConvertStringToStabAction(stablist.get<string>("TDS"));
        StabilisationAction inertia  = ConvertStringToStabAction(stablist.get<string>("INERTIA"));
        StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
        StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
        StabilisationAction vstab    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
        StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
        StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
        StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

        // --------------------------------------------------
        // set parameters for turbulence model
        ParameterList& turbmodelparams    = params.sublist("turbulence model");

        // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
        double Cs    = 0.0;
        double l_tau = 0.0;

        // the default action is no model        
        TurbModelAction turb_mod_action = no_model;
        
        if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") != "none")
        {
        
          string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");

 
          if (physical_turbulence_model == "Smagorinsky")
          {
            // the classic Smagorinsky model only requires one constant parameter
            turb_mod_action = smagorinsky;
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
            
          }
          else if (physical_turbulence_model == "Smagorinsky with van Driest damping")
          {
            // for the Smagorinsky model with van Driest damping, we need a viscous length to determine
            // the y+ (heigth in wall units)
            turb_mod_action = smagorinsky_with_wall_damping;
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
            l_tau           = turbmodelparams.get<double>("L_TAU");
          }
          else
          {
            dserror("Up to now, only Smagorinsky with and without wall function is available");
          }
        }

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
                                 vstab,
                                 cstab,
                                 cross,
                                 reynolds,
                                 turb_mod_action,
                                 Cs,
                                 l_tau,
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
             ~ n     u    - u     ~ n   / 1.0-gamma \
            acc  <-  --------- - acc * |  ---------  |
                     gamma*dt           \   gamma   /
        */
        {
          const double dt     = params.get<double>("dt");
          const double gamma  = params.get<double>("gamma");

          sub_acc_old_ = (sub_vel_-sub_vel_old_)/(gamma*dt)
                         -
                         sub_acc_old_*(1.0-gamma)/gamma;
        }
        // most recent subscale velocity becomes the old subscale velocity
        // for the next timestep
        //
        //  ~n   ~n+1
        //  u <- u
        //
        sub_vel_old_=sub_vel_;
      break;
      case calc_fluid_stationary_systemmat_and_residual:
      {
          // need current velocity/pressure 
          RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
          if (velnp==null)
            dserror("Cannot get state vector 'velnp'");

          // extract local values from the global vector
          vector<double> myvelnp(lm.size());
          DRT::Utils::ExtractMyValues(*velnp,myvelnp,lm);


          if (is_ale_)
          {
        	  dserror("No ALE support within stationary fluid solver.");
          }

          // split velocity and pressure
          // create blitz objects

          const int numnode = NumNode();
          blitz::Array<double, 1> eprenp(numnode);
          blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());

          for (int i=0;i<numnode;++i)
          {
            evelnp(0,i) = myvelnp[0+(i*4)];
            evelnp(1,i) = myvelnp[1+(i*4)];
            evelnp(2,i) = myvelnp[2+(i*4)];

            eprenp(i) = myvelnp[3+(i*4)];
          }

          // get control parameter
          const double pseudotime = params.get<double>("total time",-1.0);
          if (pseudotime < 0.0)
        	  dserror("no value for total (pseudo-)time in the parameter list");

          bool newton = params.get<bool>("include reactive terms for linearisation",false);
          bool pstab  =true;
          bool supg   =true;
          bool vstab  =false;  // viscous stabilisation part switched off !!
          bool cstab  =true;        

          // wrap epetra serial dense objects in blitz objects
          blitz::Array<double, 2> estif(elemat1.A(),
                                        blitz::shape(elemat1.M(),elemat1.N()),
                                        blitz::neverDeleteData,
                                        blitz::ColumnMajorArray<2>());
          blitz::Array<double, 1> eforce(elevec1.Values(),
                                         blitz::shape(elevec1.Length()),
                                         blitz::neverDeleteData);

          // calculate element coefficient matrix and rhs         
          StationaryImpl()->Sysmat(this,
                         evelnp,
                         eprenp,
                         estif,
                         eforce,
                         actmat,
                         pseudotime,
                         newton ,
                         pstab  ,
                         supg   ,
                         vstab  ,
                         cstab  );

          // This is a very poor way to transport the density to the
          // outside world. Is there a better one?
          params.set("density", actmat->m.fluid->density);	  
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
        rule = intrule_tet_5point;
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
    Epetra_SerialDenseMatrix    xjm(NSD,NSD);

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
  case hex8: case hex20: case hex27: case tet10: case wedge15:
    hoel = true;
    break;
  case tet4: case wedge6: case pyramid5:
    hoel = false;
    break;
  default:
    dserror("distype unknown!");
  }
  return hoel;
}

//
// check for element rewinding based on Jacobian determinant
//
bool DRT::Elements::Fluid3::checkRewinding()
{
  const DiscretizationType distype = this->Shape();
  const int iel = NumNode();
  // use one point gauss rule to calculate tau at element center
  GaussRule3D integrationrule_1point = intrule_hex_1point;
  switch(distype)
  {
  case hex8: case hex20: case hex27:
      integrationrule_1point = intrule_hex_1point;
      break;
  case tet4: case tet10:
      integrationrule_1point = intrule_tet_1point;
      break;
  case wedge6: case wedge15:
      integrationrule_1point = intrule_wedge_1point;
      break;
  case pyramid5:
      integrationrule_1point = intrule_pyramid_1point;
      break;
  default:
      dserror("invalid discretization type for fluid3");
  }
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_1point);

  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);
  Epetra_SerialDenseMatrix    xyze(NSD,iel);
  DRT::Utils::shape_function_3D_deriv1(deriv,intpoints.qxg[0][0],intpoints.qxg[0][1],intpoints.qxg[0][2],distype);
  // get node coordinates
  DRT::Node** nodes = this->Nodes();
  for (int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

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
  Epetra_SerialDenseMatrix xjm(NSD,NSD);

  xjm.Multiply('N','T',1.0,deriv,xyze,0.0);


  
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);
  if (det < 0.0) return true;

  return false;
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
  bool dofillcompleteagain = false;
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_fluid3) continue;
    DRT::Elements::Fluid3* actele = dynamic_cast<DRT::Elements::Fluid3*>(dis.lColElement(i));
    if (!actele) dserror("cast to Fluid3* failed");
    
    const DRT::Element::DiscretizationType distype = actele->Shape();
    bool possiblytorewind = false;
    switch(distype)
    {
    case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        possiblytorewind = true;
        break;
    case DRT::Element::tet4: case DRT::Element::tet10:
        possiblytorewind = true;
        break;
    case DRT::Element::wedge6: case DRT::Element::wedge15:
        possiblytorewind = true;
        break;
    case DRT::Element::pyramid5:
        possiblytorewind = true;
        break;
    default:
        dserror("invalid discretization type for fluid3");
    }
    
    if ( (possiblytorewind) && (!actele->donerewinding_) ) {
      actele->rewind_ = actele->checkRewinding();

      if (actele->rewind_) {
        if (distype==DRT::Element::tet4){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[1] = old_nodeids[2];
          new_nodeids[2] = old_nodeids[1];
          new_nodeids[3] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::hex8){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[4];
          new_nodeids[1] = old_nodeids[5];
          new_nodeids[2] = old_nodeids[6];
          new_nodeids[3] = old_nodeids[7];
          new_nodeids[4] = old_nodeids[0];
          new_nodeids[5] = old_nodeids[1];
          new_nodeids[6] = old_nodeids[2];
          new_nodeids[7] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::wedge6){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[3];
          new_nodeids[1] = old_nodeids[4];
          new_nodeids[2] = old_nodeids[5];
          new_nodeids[3] = old_nodeids[0];
          new_nodeids[4] = old_nodeids[1];
          new_nodeids[5] = old_nodeids[2];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype == DRT::Element::pyramid5){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[1] = old_nodeids[3];
          new_nodeids[3] = old_nodeids[1];
          // the other nodes can stay the same
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[2] = old_nodeids[2];
          new_nodeids[4] = old_nodeids[4];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else dserror("no rewinding scheme for this type of fluid3");
      }
      // process of rewinding done
      actele->donerewinding_ = true;
      dofillcompleteagain = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  if(dofillcompleteagain) dis.FillComplete(false,false,false);
  
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
