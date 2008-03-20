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

// different implementations
#include "fluid2_genalpha_resVMM.H"

using namespace DRT::UTILS;



/* ----------------------------------------------------------------------
 |                                                            gammi 02/08|

  Depending on the type of the algorithm (the implementation) and the
  element type (tri, quad etc.), the elements allocate common static
  arrays.

  That means that for example all quad4 fluid elements of the stationary
  implementation have a pointer f4 to the same 'implementation class'
  containing all the element arrays for eight noded elements, and all
  tri3 fluid elements of the same problem have a pointer f3 to
  the 'implementation class' containing all the element arrays for the
  3 noded element.

  */

DRT::ELEMENTS::Fluid2GenalphaResVMM* DRT::ELEMENTS::Fluid2::GenalphaResVMM()
{
  switch (NumNode())
  {
  case 4:
  {
    static Fluid2GenalphaResVMM* f4;
    if (f4==NULL)
      f4 = new Fluid2GenalphaResVMM(4);
    return f4;
  }
  case 8:
  {
    static Fluid2GenalphaResVMM* f8;
    if (f8==NULL)
      f8 = new Fluid2GenalphaResVMM(8);
    return f8;
  }
  case 9:
  {
    static Fluid2GenalphaResVMM* f9;
    if (f9==NULL)
      f9 = new Fluid2GenalphaResVMM(9);
    return f9;
  }
  case 3:
  {
    static Fluid2GenalphaResVMM* f3;
    if (f3==NULL)
      f3 = new Fluid2GenalphaResVMM(3);
    return f3;
  }
  case 6:
  {
    static Fluid2GenalphaResVMM* f6;
    if (f6==NULL)
      f6 = new Fluid2GenalphaResVMM(6);
    return f6;
  }
  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 // converts a string into an stabilisation action for this element
 //                                                          gammi 02/08
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::StabilisationAction DRT::ELEMENTS::Fluid2::ConvertStringToStabAction(
  const string& action) const
{
  DRT::ELEMENTS::Fluid2::StabilisationAction act = stabaction_unspecified;

  map<string,StabilisationAction>::const_iterator iter=stabstrtoact_.find(action);

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
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid2::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Fluid2::ActionType act = Fluid2::none;

  // set default value for (at the moment still necessary) control parameter
 bool is_stationary = false;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_fluid_systemmat_and_residual")
	act = Fluid2::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
  {
  	act = Fluid2::calc_fluid_stationary_systemmat_and_residual;
  	is_stationary = true;
  }
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
  {
  	act = Fluid2::calc_fluid_genalpha_sysmat_and_residual;
  }
  else if (action == "time update for subscales")
  {
    act = Fluid2::calc_fluid_genalpha_update_for_subscales;
  }
  else if (action == "time average for subscales and residual")
  {
    act = Fluid2::calc_fluid_genalpha_average_for_subscales_and_residual;
  }
  else
  {

    char errorout[200];
    sprintf(errorout,"Unknown type of action (%s) for Fluid2",action.c_str());

    dserror(errorout);
  }

  // get the material
  RCP<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  switch(act)
  {
    case calc_fluid_stationary_systemmat_and_residual:
    	//no break here at the moment since we use the same code layout
    	//for both stationary and instationary problems controlled
    	//by the value of is_stationary.
    	//It is planned to create a separate implementation class for stationary formulation
    	//as already done in fluid3 whithin a major cleanup of the fluid2 element.
    	// g.bau  11/07
    case calc_fluid_systemmat_and_residual:
    {
      // need current velocity and history vector
      RCP<const Epetra_Vector> vel_pre_np = discretization.GetState("velnp");
      RCP<const Epetra_Vector> hist = discretization.GetState("hist");
      if (vel_pre_np==null || hist==null) dserror("Cannot get state vectors 'velnp' and/or 'hist'");

      // extract local values from the global vectors
      vector<double> my_vel_pre_np(lm.size());
      DRT::UTILS::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);
      vector<double> myhist(lm.size());
      DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

      RCP<const Epetra_Vector> dispnp;
      vector<double> mydispnp;
      RCP<const Epetra_Vector> gridv;
      vector<double> mygridv;

      if (is_ale_)
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

        gridv = discretization.GetState("gridv");
        if (gridv==null) dserror("Cannot get state vectors 'gridv'");
        mygridv.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
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
      const double time = params.get<double>("total time",-1.0);

      // One-step-Theta: timefac = theta*dt
      // BDF2:           timefac = 2/3 * dt
      double timefac = 0;
      if (not is_stationary)
      {
        timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");
      }

      // get flag for fine-scale subgrid viscosity
      string fssgv = params.get<string>("fs subgrid viscosity","No");

      // get Smagorinsky model parameter for fine-scale subgrid viscosity
      ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");
      const double Cs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

      // get fine-scale velocity
      RCP<const Epetra_Vector> fsvelprenp;
      vector<double> myfsvelnp(2*numnode);
      if (fssgv != "No")
      {
        fsvelprenp = discretization.GetState("fsvelnp");
        if (fsvelprenp==null) dserror("Cannot get state vector 'fsvelnp'");
        vector<double> myfsvelprenp(lm.size());
        DRT::UTILS::ExtractMyValues(*fsvelprenp,myfsvelprenp,lm);

        // get from "my_vel_pre_np_fs" only velocity part "myvelnp_fs"
        for (int i=0;i<numnode;++i)
        {
          myfsvelnp[0+(i*2)]=myfsvelprenp[0+(i*3)];
          myfsvelnp[1+(i*2)]=myfsvelprenp[1+(i*3)];
        }
      }
      else
      {
        for (int i=0;i<numnode;++i)
        {
          myfsvelnp[0+(i*2)]=0.0;
          myfsvelnp[1+(i*2)]=0.0;
        }
      }

      // calculate element coefficient matrix and rhs
      f2_sys_mat(lm,myvelnp,myfsvelnp,myprenp,myvhist,mydispnp,mygridv,&elemat1,&elevec1,actmat,time,timefac,fssgv,Cs,is_stationary);

      // This is a very poor way to transport the density to the
      // outside world. Is there a better one?
      params.set("density", actmat->m.fluid->density);
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
      DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

      vector<double> myvelaf(lm.size());
      DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);

      vector<double> myaccam(lm.size());
      DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);

      // create blitz matrix objects
      const int numnode = NumNode();
      blitz::Array<double, 1> eprenp    (  numnode);
      blitz::Array<double, 2> evelnp    (2,numnode,blitz::ColumnMajorArray<2>());
      blitz::Array<double, 2> evelaf    (2,numnode,blitz::ColumnMajorArray<2>());
      blitz::Array<double, 2> eaccam    (2,numnode,blitz::ColumnMajorArray<2>());
      blitz::Array<double, 2> edispnp   (2,numnode,blitz::ColumnMajorArray<2>());
      blitz::Array<double, 2> egridvelaf(2,numnode,blitz::ColumnMajorArray<2>());


      // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
      // Additionally only the 'velocity' components of my_velaf
      // and my_accam are important!
      for (int i=0;i<numnode;++i)
      {
        eprenp(i)   = myvelnp[2+(i*3)];

        evelnp(0,i) = myvelnp[0+(i*3)];
        evelnp(1,i) = myvelnp[1+(i*3)];

        evelaf(0,i) = myvelaf[0+(i*3)];
        evelaf(1,i) = myvelaf[1+(i*3)];

        eaccam(0,i) = myaccam[0+(i*3)];
        eaccam(1,i) = myaccam[1+(i*3)];
      }

      if(is_ale_)
      {
        // get most recent displacements
        RefCountPtr<const Epetra_Vector> dispnp
          =
          discretization.GetState("dispnp");

        // get intermediate grid velocities
        RefCountPtr<const Epetra_Vector> gridvelaf
          =
          discretization.GetState("gridvelaf");

        if (dispnp==null || gridvelaf==null)
        {
          dserror("Cannot get state vectors 'dispnp' and/or 'gridvelaf'");
        }

        vector<double> mydispnp(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

        vector<double> mygridvelaf(lm.size());
        DRT::UTILS::ExtractMyValues(*gridvelaf,mygridvelaf,lm);

        // extract velocity part from "mygridvelaf" and get
        // set element displacements
        for (int i=0;i<numnode;++i)
        {
          egridvelaf(0,i) = mygridvelaf[0+(i*3)];
          egridvelaf(1,i) = mygridvelaf[1+(i*3)];

          edispnp(0,i)    = mydispnp   [0+(i*3)];
          edispnp(1,i)    = mydispnp   [1+(i*3)];
        }
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
      ParameterList& stablist = params.sublist("STABILIZATION");


      // if not available, define map from string to action
      if(stabstrtoact_.empty())
      {
        stabstrtoact_["quasistatic"    ]=subscales_quasistatic;
        stabstrtoact_["time_dependent" ]=subscales_time_dependent;
        stabstrtoact_["no_transient"   ]=inertia_stab_drop;
        stabstrtoact_["yes_transient"  ]=inertia_stab_keep;
        stabstrtoact_["no_pspg"        ]=pstab_assume_inf_sup_stable;
        stabstrtoact_["yes_pspg"       ]=pstab_use_pspg;
        stabstrtoact_["no_supg"        ]=convective_stab_none;
        stabstrtoact_["yes_supg"       ]=convective_stab_supg;
        stabstrtoact_["no_vstab"       ]=viscous_stab_none;
        stabstrtoact_["vstab_gls"      ]=viscous_stab_gls;
        stabstrtoact_["vstab_gls_rhs"  ]=viscous_stab_gls_only_rhs;
        stabstrtoact_["vstab_usfem"    ]=viscous_stab_usfem;
        stabstrtoact_["vstab_usfem_rhs"]=viscous_stab_usfem_only_rhs;
        stabstrtoact_["no_ctab"        ]=continuity_stab_none;
        stabstrtoact_["cstab_qs"       ]=continuity_stab_yes;
        stabstrtoact_["cstab_td"       ]=continuity_stab_td;
        stabstrtoact_["no_cross"       ]=cross_stress_stab_none;
        stabstrtoact_["cross_complete" ]=cross_stress_stab;
        stabstrtoact_["cross_rhs"      ]=cross_stress_stab_only_rhs;
        stabstrtoact_["no_reynolds"    ]=reynolds_stress_stab_none;
        stabstrtoact_["reynolds_rhs"   ]=reynolds_stress_stab_only_rhs;
      }

      StabilisationAction tds      = ConvertStringToStabAction(stablist.get<string>("TDS"));
      StabilisationAction inertia  = ConvertStringToStabAction(stablist.get<string>("TRANSIENT"));
      StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
      StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
      StabilisationAction vstab    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
      StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
      StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
      StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

      // --------------------------------------------------
      // specify what to compute
      const bool compute_elemat = params.get<bool>("compute element matrix");

      // --------------------------------------------------
      // calculate element coefficient matrix
      GenalphaResVMM()->Sysmat(this,
                               elemat1,
                               elevec1,
                               edispnp,
                               egridvelaf,
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
                               compute_elemat
        );

      // This is a very poor way to transport the density to the
      // outside world. Is there a better one?
      params.set("density", actmat->m.fluid->density);
    }
    break;
    case calc_fluid_genalpha_update_for_subscales:
    {
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

      const double dt     = params.get<double>("dt");
      const double gamma  = params.get<double>("gamma");

      sub_acc_old_ = (sub_vel_-sub_vel_old_)/(gamma*dt)
                      -
                      sub_acc_old_*(1.0-gamma)/gamma;

      // most recent subscale velocity becomes the old subscale velocity
      // for the next timestep
      //
      //  ~n   ~n+1
      //  u <- u
      //
      sub_vel_old_=sub_vel_;
    }
    case calc_fluid_genalpha_average_for_subscales_and_residual:
    {
      // nothing at this moment
    }
    break;
    default:
      dserror("Unknown type of action for Fluid2");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid2::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      g.bau 07/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid2 elements, the          |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid2::EvaluateNeumann(ParameterList& params,
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
void DRT::ELEMENTS::Fluid2::f2_sys_mat(vector<int>&              lm,
                                       vector<double>&           evelnp,
                                       vector<double>&           fsevelnp,
                                       vector<double>&           eprenp,
                                       vector<double>&           evhist,
                                       vector<double>&           edispnp,
                                       vector<double>&           egridv,
                                       Epetra_SerialDenseMatrix* sys_mat,
                                       Epetra_SerialDenseVector* residual,
                                       struct _MATERIAL*         material,
                                       double                    time,
                                       double                    timefac,
                                       string                    fssgv,
                                       const double              Cs,
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
  static Epetra_SerialDenseMatrix 	fsvderxy(2,2);
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
  double         		vel_norm, fsvel_norm, pe, re, xi1, xi2, xi;
  vector<double>         	velino(2); /* normed velocity at element centre */
  double         		mk=0.0;
  static vector<double>         velint(2);
  static vector<double>         fsvelint(2);
  static vector<double>         tau(3); // stab parameters
  double                        vart; // artificial subgrid viscosity

  /*------------------------------------------------------- initialise ---*/
    // use one point gauss rule to calculate tau at element center
  GaussRule2D integrationrule_stabili = intrule2D_undefined;
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
  DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
  DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

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

  /*------------------------------------ compute subgrid viscosity nu_art ---*/
  if (fssgv == "artificial_all" || fssgv == "artificial_small")
  {
    if (fssgv == "artificial_small")
    {
      /*--------------------- get fine-scale velocities at element center ---*/
      for (int i=0;i<2;i++)
      {
        fsvelint[i]=0.0;
        for (int j=0;j<iel;j++)
        {
          fsvelint[i] += funct[j]*fsevelnp[i+(2*j)];
        }
      } //end loop over i

      /*----------------------------------------  get fine-scale vel_norm ---*/
      fsvel_norm = sqrt(DSQR(fsvelint[0]) + DSQR(fsvelint[1]));
    }
    /*-----------------------------------------  get all-scale vel_norm ---*/
    else fsvel_norm = vel_norm;

    /*----------------------- compute artificial subgrid viscosity nu_art ---*/
    re = mk * fsvel_norm * hk / visc;     /* convective : viscous forces */
    xi = DMAX(re,1.0);

    vart = (DSQR(hk)*mk*DSQR(fsvel_norm))/(2.0*visc*xi);
  }
  else if (fssgv == "Smagorinsky_all" || fssgv == "Smagorinsky_small" ||
           fssgv == "mixed_Smagorinsky_all" || fssgv == "mixed_Smagorinsky_small")
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                               +-                                 -+ 1
    //                           2   |          / h \           / h \    | -
    //    visc          = (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent              |          \   / ij        \   / ij |
    //                               +-                                 -+
    //                               |                                   |
    //                               +-----------------------------------+
    //                                    'resolved' rate of strain
    //

    double rateofstrain = 0.0;
    {
      /*---------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

      if (fssgv == "Smagorinsky_small" || fssgv == "mixed_Smagorinsky_small")
      {
        /*---- get fine-scale velocity (np,i) derivatives at element center */
        for (int i=0;i<2;i++)
        {
          fsvderxy(0,i)=0.0;
          fsvderxy(1,i)=0.0;
          for (int j=0;j<iel;j++)
          {
            fsvderxy(0,i) += derxy(i,j)*fsevelnp[0+(2*j)];
            fsvderxy(1,i) += derxy(i,j)*fsevelnp[1+(2*j)];
          } /* end of loop over j */
        } /* end of loop over i */
      }
      else
      {
        /*---- get all-scale velocity (np,i) derivatives at element center */
        for (int i=0;i<2;i++)
        {
          fsvderxy(0,i)=0.0;
          fsvderxy(1,i)=0.0;
          for (int j=0;j<iel;j++)
          {
            fsvderxy(0,i) += derxy(i,j)*evelnp[0+(2*j)];
            fsvderxy(1,i) += derxy(i,j)*evelnp[1+(2*j)];
          } /* end of loop over j */
        } /* end of loop over i */
      }

      Epetra_SerialDenseMatrix  epsilon(2,2);

      /*-------------------------- get rate-of-strain tensor epsilon(u) */
      for (int i=0;i<2;i++)
      {
        for (int j=0;j<2;j++)
        {
          epsilon(i,j) = 0.5 * ( fsvderxy(i,j) + fsvderxy(j,i) );
        } /* end of loop over j */
      } /* end of loop over i */

      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<2;mm++)
        {
          rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
        }
      }
      rateofstrain *= 2.0;
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

    vart = Cs * Cs * hk * hk * rateofstrain;

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

      /*----------- get fine-scale velocity (np,i) deriv. at int. point */
      for (int i=0;i<2;i++)
      {
        fsvderxy(0,i)=0.0;
        fsvderxy(1,i)=0.0;
        if (fssgv != "No")
        {
          for (int j=0;j<iel;j++)
          {
            fsvderxy(0,i) += derxy(i,j)*fsevelnp[0+(2*j)];
            fsvderxy(1,i) += derxy(i,j)*fsevelnp[1+(2*j)];
          } /* end of loop over j */
        }
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
     for (int isd=0;isd<2;isd++)
     {
       edeadng[isd] = 0.0;
       for (int inode=0;inode<iel;inode++)
       {
         edeadng[isd]+= bodyforce(isd,inode)*funct[inode];
       }
     }
     edeadng[2] = 0.0;

      // perform integration for entire matrix and rhs
      if(is_stationary==false)
        f2_calmat(*sys_mat,*residual,
                  velint,histvec,gridvelint,press,
                  vderxy,fsvderxy,vderxy2,gradp,
                  funct,tau,vart,
                  derxy,derxy2,edeadng,
                  fac,visc,iel,fssgv,
                  timefac);
      else
        f2_calmat_stationary(*sys_mat,*residual,
                             velint,histvec,gridvelint,press,
                             vderxy,fsvderxy,vderxy2,gradp,
                             funct,tau,vart,
                             derxy,derxy2,edeadng,
                             fac,visc,iel,fssgv);


  } // end of loop over integration points


  return;
} // DRT::ELEMENTS::Fluid2::f2_sys_mat


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
void DRT::ELEMENTS::Fluid2::f2_jaco(const Epetra_SerialDenseMatrix& xyze,
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

} //end of DRT::ELEMENTS::Fluid2::f2_jaco


/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) g.bau 07/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a surface Neumann condition   |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::ELEMENTS::Fluid2::f2_getbodyforce(
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
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative time value
	// indicates an error.
         dserror("Negative time value in body force calculation: time = %f",time);
	// curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
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
} // end of DRT:ELEMENTS:Fluid2:f2_getbodyforce


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
void DRT::ELEMENTS::Fluid2::f2_gder(Epetra_SerialDenseMatrix& derxy,
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
} // end of DRT:ELEMENTS:Fluid2:f2_gder

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

void DRT::ELEMENTS::Fluid2::f2_gder2(const Epetra_SerialDenseMatrix& xyze,
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
} // end of DRT:ELEMENTS:Fluid2:f2_gder2



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
\param **fsvderxy   DOUBLE        (i)   global fine-scale vel derivatives
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

void DRT::ELEMENTS::Fluid2::f2_calmat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    vector<double>&           histvec,
    vector<double>&           gridvint,
    double&   	              press,
    Epetra_SerialDenseMatrix& vderxy,
    Epetra_SerialDenseMatrix& fsvderxy,
    Epetra_SerialDenseMatrix& vderxy2,
    vector<double>&           gradp,
    Epetra_SerialDenseVector& funct,
    vector<double>&           tau,
    const double&             vart,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const vector<double>&     edeadng,
    const double&             fac,
    const double&             visc,
    const int&                iel,
    string                    fssgv,
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

  // parameter for artificial subgrid viscosity
  double vartfac=0.0;
  if (fssgv != "No") vartfac = vart*timefacfac;

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
#define vartfac        vartfac
#define fssgv          fssgv
#define fsvderxy_(i,j) fsvderxy(i,j)

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
#undef vartfac
#undef fssgv
#undef fsvderxy_

  return;
} // end of DRT:ELEMENTS:Fluid2:f2_calmat


void DRT::ELEMENTS::Fluid2::f2_calmat_stationary(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    vector<double>&           velint,
    vector<double>&           histvec,
    vector<double>&           gridvint,
    double&   	              press,
    Epetra_SerialDenseMatrix& vderxy,
    Epetra_SerialDenseMatrix& fsvderxy,
    Epetra_SerialDenseMatrix& vderxy2,
    vector<double>&           gradp,
    Epetra_SerialDenseVector& funct,
    vector<double>&           tau,
    const double&             vart,
    Epetra_SerialDenseMatrix& derxy,
    Epetra_SerialDenseMatrix& derxy2,
    const vector<double>&     edeadng,
    const double&             fac,
    const double&             visc,
    const int&                iel,
    string                    fssgv
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

// parameter for artificial subgrid viscosity
double vartfac=0.0;
if (fssgv != "No") vartfac = vart*fac;

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
#define vartfac        vartfac
#define fssgv          fssgv
#define fsvderxy_(i,j) fsvderxy(i,j)

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
#undef vartfac
#undef fssgv
#undef fsvderxy_

return;
} // end of DRT:ELEMENTS:Fluid2:f2_calmat_stationary

// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::ELEMENTS::Fluid2::is_higher_order_element(
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
int DRT::ELEMENTS::Fluid2Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
