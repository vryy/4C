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
#include "fluid2_impl.H"
#include "fluid2_stationary.H"
#include "fluid2_genalpha_resVMM.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/sutherland_fluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include <Epetra_SerialDenseSolver.h>

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

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action == "calc_fluid_systemmat_and_residual")
    act = Fluid2::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Fluid2::calc_fluid_stationary_systemmat_and_residual;
  else if (action == "calc_fluid_afgenalpha_systemmat_and_residual")
    act = Fluid2::calc_fluid_afgenalpha_systemmat_and_residual;
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
    act = Fluid2::calc_fluid_genalpha_sysmat_and_residual;
  else if (action == "time update for subscales")
    act = Fluid2::calc_fluid_genalpha_update_for_subscales;
  else if (action == "time average for subscales and residual")
    act = Fluid2::calc_fluid_genalpha_average_for_subscales_and_residual;
  else if (action == "get_density")
    act = Fluid2::get_density;
  else
  {

    char errorout[200];
    sprintf(errorout,"Unknown type of action (%s) for Fluid2",action.c_str());

    dserror(errorout);
  }

  // get the material
  RCP<MAT::Material> mat = Material();

  switch(act)
  {
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // the standard one-step-theta + BDF2 implementation as well as
    // generalized-alpha implementation with continuity equation at n+alpha_F
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    case calc_fluid_systemmat_and_residual:
    case calc_fluid_afgenalpha_systemmat_and_residual:
    {
      // if not available, define map from string to action
      if(stabstrtoact_.empty())
      {
        stabstrtoact_["no_pspg"          ]=pstab_assume_inf_sup_stable;
        stabstrtoact_["yes_pspg"         ]=pstab_use_pspg;
        stabstrtoact_["no_supg"          ]=convective_stab_none;
        stabstrtoact_["yes_supg"         ]=convective_stab_supg;
        stabstrtoact_["no_vstab"         ]=viscous_stab_none;
        stabstrtoact_["vstab_gls"        ]=viscous_stab_gls;
        stabstrtoact_["vstab_gls_rhs"    ]=viscous_stab_gls_only_rhs;
        stabstrtoact_["vstab_usfem"      ]=viscous_stab_usfem;
        stabstrtoact_["vstab_usfem_rhs"  ]=viscous_stab_usfem_only_rhs;
        stabstrtoact_["no_cstab"         ]=continuity_stab_none;
        stabstrtoact_["cstab_qs"         ]=continuity_stab_yes;
        stabstrtoact_["no_cross"         ]=cross_stress_stab_none;
        stabstrtoact_["cross_complete"   ]=cross_stress_stab;
        stabstrtoact_["cross_rhs"        ]=cross_stress_stab_only_rhs;
        stabstrtoact_["no_reynolds"      ]=reynolds_stress_stab_none;
        stabstrtoact_["reynolds_rhs"     ]=reynolds_stress_stab_only_rhs;
      }
      return DRT::ELEMENTS::Fluid2ImplInterface::Impl(this)->Evaluate(this,
                                                                      params,
                                                                      discretization,
                                                                      lm,
                                                                      elemat1,
                                                                      elemat2,
                                                                      elevec1,
                                                                      elevec2,
                                                                      elevec3,
                                                                      mat);
    }
    break;
    case calc_fluid_stationary_systemmat_and_residual:
    {
      // if not available, define map from string to action
      if(stabstrtoact_.empty())
      {
        stabstrtoact_["no_pspg"          ]=pstab_assume_inf_sup_stable;
        stabstrtoact_["yes_pspg"         ]=pstab_use_pspg;
        stabstrtoact_["no_supg"          ]=convective_stab_none;
        stabstrtoact_["yes_supg"         ]=convective_stab_supg;
        stabstrtoact_["no_vstab"         ]=viscous_stab_none;
        stabstrtoact_["vstab_gls"        ]=viscous_stab_gls;
        stabstrtoact_["vstab_gls_rhs"    ]=viscous_stab_gls_only_rhs;
        stabstrtoact_["vstab_usfem"      ]=viscous_stab_usfem;
        stabstrtoact_["vstab_usfem_rhs"  ]=viscous_stab_usfem_only_rhs;
        stabstrtoact_["no_cstab"         ]=continuity_stab_none;
        stabstrtoact_["cstab_qs"         ]=continuity_stab_yes;
        stabstrtoact_["no_cross"         ]=cross_stress_stab_none;
        stabstrtoact_["cross_complete"   ]=cross_stress_stab;
        stabstrtoact_["cross_rhs"        ]=cross_stress_stab_only_rhs;
        stabstrtoact_["no_reynolds"      ]=reynolds_stress_stab_none;
        stabstrtoact_["reynolds_rhs"     ]=reynolds_stress_stab_only_rhs;
      }
      return DRT::ELEMENTS::Fluid2StationaryInterface::Impl(this)->Evaluate(this,
                                                                            params,
                                                                            discretization,
                                                                            lm,
                                                                            elemat1,
                                                                            elemat2,
                                                                            elevec1,
                                                                            elevec2,
                                                                            elevec3,
                                                                            mat);
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

      // create SerialDense matrix objects
      const int numnode = NumNode();
      Epetra_SerialDenseVector eprenp    (  numnode);
      Epetra_SerialDenseMatrix evelnp    (2,numnode);
      Epetra_SerialDenseMatrix evelaf    (2,numnode);
      Epetra_SerialDenseMatrix eaccam    (2,numnode);
      Epetra_SerialDenseMatrix edispnp   (2,numnode);
      Epetra_SerialDenseMatrix egridvelaf(2,numnode);


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
      string newtonstr=params.get<string>("Linearisation");

      bool newton = false;
      if(newtonstr=="Newton")
      {
        newton=true;
      }

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
      // Now do the nurbs specific stuff
      std::vector<Epetra_SerialDenseVector> myknots(2);

      // for isogeometric elements
      if(this->Shape()==nurbs4 || this->Shape()==nurbs9)
      {
        DRT::NURBS::NurbsDiscretization* nurbsdis
            =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

        (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());
      }

      // --------------------------------------------------
      // calculate element coefficient matrix
      GenalphaResVMM()->Sysmat(this,
                               myknots,
                               elemat1,
                               elevec1,
                               edispnp,
                               egridvelaf,
                               evelnp,
                               eprenp,
                               eaccam,
                               evelaf,
                               mat,
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
      double dens = 0.0;
      if(mat->MaterialType()== INPAR::MAT::m_fluid)
      {
        MAT::NewtonianFluid* actmat = static_cast<MAT::NewtonianFluid*>(mat.get());
        dens = actmat->Density();
      }
      else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
      {
        MAT::CarreauYasuda* actmat = static_cast<MAT::CarreauYasuda*>(mat.get());
        dens = actmat->Density();
      }
      else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
      {
        MAT::ModPowerLaw* actmat = static_cast<MAT::ModPowerLaw*>(mat.get());
        dens = actmat->Density();
      }
      else
        dserror("no fluid material found");

      params.set("density", dens);
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
      for(int rr=0;rr<spren_.Length();++rr)
      {
	spren_(rr) = sprenp_(rr);
      }
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


      for(int rr=0;rr<2;++rr)
      {
	for(int mm=0;mm<spren_.Length();++mm)
	{
	  saccn_(rr,mm) =
	  (svelnp_(rr,mm)-sveln_(rr,mm))/(gamma*dt)
	    -
	    saccn_(rr,mm)*(1.0-gamma)/gamma;
	}
      }

      // most recent subscale velocity becomes the old subscale velocity
      // for the next timestep
      //
      //  ~n   ~n+1
      //  u <- u
      //
      for(int rr=0;rr<2;++rr)
      {
	for(int mm=0;mm<spren_.Length();++mm)
	{
	  sveln_(rr,mm)=svelnp_(rr,mm);
	}
      }
    }
    case calc_fluid_genalpha_average_for_subscales_and_residual:
    {
      // nothing at this moment
    }
    break;
    case get_density:
    {
      // This is a very poor way to transport the density to the
      // outside world. Is there a better one?
      if(mat->MaterialType()== INPAR::MAT::m_fluid)
      {
        MAT::NewtonianFluid* actmat = static_cast<MAT::NewtonianFluid*>(mat.get());
        params.set("density", actmat->Density());
      }
      else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
      {
        MAT::CarreauYasuda* actmat = static_cast<MAT::CarreauYasuda*>(mat.get());
        params.set("density", actmat->Density());
      }
      else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
      {
        MAT::ModPowerLaw* actmat = static_cast<MAT::ModPowerLaw*>(mat.get());
        params.set("density", actmat->Density());
      }
      else
        dserror("no fluid material found");
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


// get optimal gaussrule for discretization type
GaussRule2D DRT::ELEMENTS::Fluid2::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule2D rule = intrule2D_undefined;
    switch (distype)
    {
    case quad4: case nurbs4:
        rule = intrule_quad_4point;
        break;
    case quad8: case quad9: case nurbs9:
        rule = intrule_quad_9point;
        break;
    case tri3:
        rule = intrule_tri_3point;
        break;
    case tri6:
        rule = intrule_tri_6point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::ELEMENTS::Fluid2::isHigherOrderElement(
  const DRT::Element::DiscretizationType  distype) const
{
  bool hoel = true;
  switch (distype)
  {
    case quad4: case quad8: case quad9: case tri6: case nurbs4: case nurbs9:
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
