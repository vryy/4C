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
#include "fluid3_stationary.H"
#include "fluid3_lin_impl.H"
#include "fluid3_genalpha_resVMM.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include <blitz/array.h>
#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;

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


/*---------------------------------------------------------------------*
|  converts a string into an Action for this element                   |
*----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::ActionType DRT::ELEMENTS::Fluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Fluid3::ActionType act = Fluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = Fluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Fluid3::calc_fluid_stationary_systemmat_and_residual;
  else if (action == "calc_linear_fluid")
    act = Fluid3::calc_linear_fluid;
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
    act = Fluid3::calc_fluid_genalpha_sysmat_and_residual;
  else if (action == "time update for subscales")
    act = Fluid3::calc_fluid_genalpha_update_for_subscales;
  else if (action == "time average for subscales and residual")
    act = Fluid3::calc_fluid_genalpha_average_for_subscales_and_residual;
  else if (action == "calc_fluid_beltrami_error")
    act = Fluid3::calc_fluid_beltrami_error;
  else if (action == "calc_turbulence_statistics")
    act = Fluid3::calc_turbulence_statistics;
  else if (action == "calc_fluid_box_filter")
    act = Fluid3::calc_fluid_box_filter;
  else if (action == "calc_smagorinsky_const")
    act = Fluid3::calc_smagorinsky_const;
  else if (action == "get_density")
    act = Fluid3::get_density;
  else
    dserror("Unknown type of action for Fluid3");
  return act;
}

/*----------------------------------------------------------------------*
 // converts a string into an stabilisation action for this element
 //                                                          gammi 08/07
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::StabilisationAction DRT::ELEMENTS::Fluid3::ConvertStringToStabAction(
  const string& action) const
{
  DRT::ELEMENTS::Fluid3::StabilisationAction act = stabaction_unspecified;

  map<string,StabilisationAction>::const_iterator iter=stabstrtoact_.find(action);

  if (iter != stabstrtoact_.end())
  {
    act = (*iter).second;
  }
  else
  {
    dserror("looking for stab action (%s) not contained in map",action.c_str());
  }
  return act;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3::Evaluate(ParameterList& params,
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
  const DRT::ELEMENTS::Fluid3::ActionType act = convertStringToActionType(action);

  // get the material
  RefCountPtr<MAT::Material> mat = Material();

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_fluid)
    actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();
  else if(mat->MaterialType()== m_carreauyasuda)
    actmat = static_cast<MAT::CarreauYasuda*>(mat.get())->MaterialData();
  else if(mat->MaterialType()== m_modpowerlaw)
    actmat = static_cast<MAT::ModPowerLaw*>(mat.get())->MaterialData();
  else
    dserror("fluid material expected but got type %d", mat->MaterialType());

  switch(act)
  {
      //--------------------------------------------------
      //--------------------------------------------------
      // the standard one-step-theta implementation
      //--------------------------------------------------
      //--------------------------------------------------
      case calc_fluid_systemmat_and_residual:
      {
        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["no_pspg"                ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["yes_pspg"               ]=pstab_use_pspg;
          stabstrtoact_["no_supg"                ]=convective_stab_none;
          stabstrtoact_["yes_supg"               ]=convective_stab_supg;
          stabstrtoact_["no_vstab"               ]=viscous_stab_none;
          stabstrtoact_["vstab_gls"              ]=viscous_stab_gls;
          stabstrtoact_["vstab_gls_rhs"          ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["vstab_usfem"            ]=viscous_stab_usfem;
          stabstrtoact_["vstab_usfem_rhs"        ]=viscous_stab_usfem_only_rhs;
          stabstrtoact_["no_cstab"               ]=continuity_stab_none;
          stabstrtoact_["cstab_qs"               ]=continuity_stab_yes;
          stabstrtoact_["no_cross"               ]=cross_stress_stab_none;
          stabstrtoact_["cross_complete"         ]=cross_stress_stab;
          stabstrtoact_["cross_rhs"              ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no_reynolds"            ]=reynolds_stress_stab_none;
          stabstrtoact_["reynolds_rhs"           ]=reynolds_stress_stab_only_rhs;
          stabstrtoact_["No"                     ]=Fluid3::fssgv_no;
          stabstrtoact_["artificial_all"         ]=Fluid3::fssgv_artificial_all;
          stabstrtoact_["artificial_small"       ]=Fluid3::fssgv_artificial_small;
          stabstrtoact_["Smagorinsky_all"        ]=Fluid3::fssgv_Smagorinsky_all;
          stabstrtoact_["Smagorinsky_small"      ]=Fluid3::fssgv_Smagorinsky_small;
        }

        return DRT::ELEMENTS::Fluid3ImplInterface::Impl(this)->Evaluate(
               this,
               params,
               discretization,
               lm,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3,
               mat,
               actmat);
      }
      break;
      case calc_fluid_stationary_systemmat_and_residual:
      {
        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["no_pspg"                ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["yes_pspg"               ]=pstab_use_pspg;
          stabstrtoact_["no_supg"                ]=convective_stab_none;
          stabstrtoact_["yes_supg"               ]=convective_stab_supg;
          stabstrtoact_["no_vstab"               ]=viscous_stab_none;
          stabstrtoact_["vstab_gls"              ]=viscous_stab_gls;
          stabstrtoact_["vstab_gls_rhs"          ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["vstab_usfem"            ]=viscous_stab_usfem;
          stabstrtoact_["vstab_usfem_rhs"        ]=viscous_stab_usfem_only_rhs;
          stabstrtoact_["no_cstab"               ]=continuity_stab_none;
          stabstrtoact_["cstab_qs"               ]=continuity_stab_yes;
          stabstrtoact_["no_cross"               ]=cross_stress_stab_none;
          stabstrtoact_["cross_complete"         ]=cross_stress_stab;
          stabstrtoact_["cross_rhs"              ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no_reynolds"            ]=reynolds_stress_stab_none;
          stabstrtoact_["reynolds_rhs"           ]=reynolds_stress_stab_only_rhs;
          stabstrtoact_["No"                     ]=Fluid3::fssgv_no;
          stabstrtoact_["artificial_all"         ]=Fluid3::fssgv_artificial_all;
          stabstrtoact_["artificial_small"       ]=Fluid3::fssgv_artificial_small;
          stabstrtoact_["Smagorinsky_all"        ]=Fluid3::fssgv_Smagorinsky_all;
          stabstrtoact_["Smagorinsky_small"      ]=Fluid3::fssgv_Smagorinsky_small;
        }

        return DRT::ELEMENTS::Fluid3StationaryImplInterface::Impl(this)->Evaluate(
               this,
               params,
               discretization,
               lm,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3,
               mat,
               actmat);
      }
      break;
      //--------------------------------------------------
      //--------------------------------------------------
      // the generalized-alpha implementation
      //--------------------------------------------------
      //--------------------------------------------------
      case calc_fluid_genalpha_sysmat_and_residual:
      {
        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["quasistatic"            ]=subscales_quasistatic;
          stabstrtoact_["time_dependent"         ]=subscales_time_dependent;
          stabstrtoact_["no_transient"           ]=inertia_stab_drop;
          stabstrtoact_["yes_transient"          ]=inertia_stab_keep;
          stabstrtoact_["no_pspg"                ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["yes_pspg"               ]=pstab_use_pspg;
          stabstrtoact_["no_supg"                ]=convective_stab_none;
          stabstrtoact_["yes_supg"               ]=convective_stab_supg;
          stabstrtoact_["no_vstab"               ]=viscous_stab_none;
          stabstrtoact_["vstab_gls"              ]=viscous_stab_gls;
          stabstrtoact_["vstab_gls_rhs"          ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["vstab_usfem"            ]=viscous_stab_usfem;
          stabstrtoact_["vstab_usfem_rhs"        ]=viscous_stab_usfem_only_rhs;
          stabstrtoact_["no_cstab"               ]=continuity_stab_none;
          stabstrtoact_["cstab_qs"               ]=continuity_stab_yes;
          stabstrtoact_["cstab_td"               ]=continuity_stab_td;
          stabstrtoact_["no_cross"               ]=cross_stress_stab_none;
          stabstrtoact_["cross_complete"         ]=cross_stress_stab;
          stabstrtoact_["cross_rhs"              ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no_reynolds"            ]=reynolds_stress_stab_none;
          stabstrtoact_["reynolds_rhs"           ]=reynolds_stress_stab_only_rhs;
          stabstrtoact_["No"                     ]=Fluid3::fssgv_no;
          stabstrtoact_["artificial_all"         ]=Fluid3::fssgv_artificial_all;
          stabstrtoact_["artificial_small"       ]=Fluid3::fssgv_artificial_small;
          stabstrtoact_["Smagorinsky_all"        ]=Fluid3::fssgv_Smagorinsky_all;
          stabstrtoact_["Smagorinsky_small"      ]=Fluid3::fssgv_Smagorinsky_small;
        }

        return DRT::ELEMENTS::Fluid3GenalphaResVMMInterface::Impl(this)->Evaluate(
               this,
               params,
               discretization,
               lm,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3,
               mat,
               actmat);
      }
      break;
      case calc_linear_fluid:
      {
        return DRT::ELEMENTS::Fluid3lin_ImplInterface::Impl(this)->Evaluate(
               this,
               params,
               discretization,
               lm,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3,
               mat,
               actmat);
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
          DRT::UTILS::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

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
          DRT::UTILS::ExtractMyValues(*velnp,mysol,lm);

          // integrate mean values
          const DiscretizationType distype = this->Shape();


          switch (distype)
          {
          case DRT::Element::hex8:
          {
            f3_calc_means<8>(discretization,mysol,params);
            break;
          }
          case DRT::Element::hex20:
          {
            f3_calc_means<20>(discretization,mysol,params);
            break;
          }
          case DRT::Element::hex27:
          {
            f3_calc_means<27>(discretization,mysol,params);
            break;
          }
          case DRT::Element::nurbs8:
          {
            f3_calc_means<8>(discretization,mysol,params);
            break;
          }
          case DRT::Element::nurbs27:
          {
            f3_calc_means<27>(discretization,mysol,params);
            break;
          }
          default:
          {
            dserror("Unknown element type for mean value evaluation\n");
          }
          }
        }
      }
      break;
      case calc_fluid_box_filter:
      {
        // --------------------------------------------------
        // extract velocities from the global distributed vectors

        // velocity and pressure values (most recent
        // intermediate solution, i.e. n+alphaF for genalpha
        // and n+1 for one-step-theta)
        RefCountPtr<const Epetra_Vector> vel =
          discretization.GetState("u and p (trial)");

        if (vel==null)
        {
          dserror("Cannot get state vectors 'vel'");
        }

        // extract local values from the global vectors
        vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        // initialise the contribution of this element to the patch volume to zero
        double volume_contribution = 0;

        // integrate the convolution with the box filter function for this element
        // the results are assembled onto the *_hat arrays

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
        case DRT::Element::hex8:
        {
          this->f3_apply_box_filter<8>(myvel,
                                       elevec1.Values(),
                                       elemat1.A(),
                                       elemat2.A(),
                                       volume_contribution);
          break;
        }
        case DRT::Element::tet4:
        {
          this->f3_apply_box_filter<4>(myvel,
                                       elevec1.Values(),
                                       elemat1.A(),
                                       elemat2.A(),
                                       volume_contribution);
          break;
        }
        default:
        {
          dserror("Unknown element type for box filter application\n");
        }
        }

        // hand down the volume contribution to the time integration algorithm
        params.set<double>("volume_contribution",volume_contribution);
      }
      break;
      case calc_smagorinsky_const:
      {
        RefCountPtr<Epetra_MultiVector> filtered_vel                        =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_vel");
        RefCountPtr<Epetra_MultiVector> col_filtered_reynoldsstress         =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_reynoldsstress");
        RefCountPtr<Epetra_MultiVector> col_filtered_modeled_subgrid_stress =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_modeled_subgrid_stress");

        double LijMij   = 0;
        double MijMij   = 0;
        double center   = 0;

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
        case DRT::Element::hex8:
        {
          this->f3_calc_smag_const_LijMij_and_MijMij<8>(
            filtered_vel                       ,
            col_filtered_reynoldsstress        ,
            col_filtered_modeled_subgrid_stress,
            LijMij                             ,
            MijMij                             ,
            center                             );
          break;
        }
        case DRT::Element::tet4:
        {
          this->f3_calc_smag_const_LijMij_and_MijMij<4>(
            filtered_vel                       ,
            col_filtered_reynoldsstress        ,
            col_filtered_modeled_subgrid_stress,
            LijMij                             ,
            MijMij                             ,
            center                             );
          break;
        }
        default:
        {
          dserror("Unknown element type for box filter application\n");
        }
        }

        // set Cs_delta_sq without averaging (only clipping)
        Cs_delta_sq_ = 0.5 * LijMij / MijMij;
        if (Cs_delta_sq_<0)
        {
          Cs_delta_sq_= 0;
        }

        params.set<double>("LijMij",LijMij);
        params.set<double>("MijMij",MijMij);
        params.set<double>("center",center);
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
      break;
      case calc_fluid_genalpha_average_for_subscales_and_residual:
      {
        if(this->Owner() == discretization.Comm().MyPID())
        {

          return DRT::ELEMENTS::Fluid3GenalphaResVMMInterface::Impl(this)->CalcResAvgs(
            this,
            params,
            discretization,
            lm,
            mat,
            actmat);
        }
      }
      break;
      case get_density:
      {
        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        if(mat->MaterialType()== m_fluid)
          params.set("density", actmat->m.fluid->density);
        else if(mat->MaterialType()== m_carreauyasuda)
          params.set("density", actmat->m.carreauyasuda->density);
        else if(mat->MaterialType()== m_modpowerlaw)
          params.set("density", actmat->m.modpowerlaw->density);
        else
          dserror("no fluid material found");
      }
      break;
      default:
        dserror("Unknown type of action for Fluid3");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

// get optimal gaussrule for discretization type
GaussRule3D DRT::ELEMENTS::Fluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule3D_undefined;
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
void DRT::ELEMENTS::Fluid3::f3_int_beltrami_err(
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
  Epetra_SerialDenseMatrix  xjm(3,3);
  Epetra_SerialDenseMatrix  deriv(3,iel);

  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  DRT::Node** nodes = Nodes();
  for(int inode=0;inode<iel;inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode)=x[0];
    xyze(1,inode)=x[1];
    xyze(2,inode)=x[2];
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
    LINALG::FixedSizeSerialDenseMatrix<NSD,NSD>    xjm;

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
    const double det = xjm.Determinant();

    if(det < 0.0)
    {
        printf("\n");
        printf("GLOBAL ELEMENT NO.%i\n",Id());
        printf("NEGATIVE JACOBIAN DETERMINANT: %f\n", det);
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
 |                      _ _           _ _           _ _
 | as well as  numele * u*v, numele * u*w, numele * v*w
 |
 | as well as numele.
 | All results are communicated via the parameter list!
 |
 *---------------------------------------------------------------------*/
template<int iel>
void DRT::ELEMENTS::Fluid3::f3_calc_means(
  DRT::Discretization&      discretization,
  vector<double>&           solution      ,
  ParameterList& 	    params
  )
{
  // get view of solution vector
  LINALG::FixedSizeSerialDenseMatrix<iel,1> sol(&(solution[0]),true);

  // set element data
  const DiscretizationType distype = this->Shape();

  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");

  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  RefCountPtr<vector<double> > planes = params.get<RefCountPtr<vector<double> > >("coordinate vector for hom. planes");

  // get the pointers to the solution vectors
  RefCountPtr<vector<double> > sumu   = params.get<RefCountPtr<vector<double> > >("mean velocity u");
  RefCountPtr<vector<double> > sumv   = params.get<RefCountPtr<vector<double> > >("mean velocity v");
  RefCountPtr<vector<double> > sumw   = params.get<RefCountPtr<vector<double> > >("mean velocity w");
  RefCountPtr<vector<double> > sump   = params.get<RefCountPtr<vector<double> > >("mean pressure p");

  RefCountPtr<vector<double> > sumsqu = params.get<RefCountPtr<vector<double> > >("mean value u^2");
  RefCountPtr<vector<double> > sumsqv = params.get<RefCountPtr<vector<double> > >("mean value v^2");
  RefCountPtr<vector<double> > sumsqw = params.get<RefCountPtr<vector<double> > >("mean value w^2");
  RefCountPtr<vector<double> > sumuv  = params.get<RefCountPtr<vector<double> > >("mean value uv");
  RefCountPtr<vector<double> > sumuw  = params.get<RefCountPtr<vector<double> > >("mean value uw");
  RefCountPtr<vector<double> > sumvw  = params.get<RefCountPtr<vector<double> > >("mean value vw");
  RefCountPtr<vector<double> > sumsqp = params.get<RefCountPtr<vector<double> > >("mean value p^2");

  if(distype == DRT::Element::hex8
     ||
     distype == DRT::Element::hex27
     ||
     distype == DRT::Element::hex20)
  {
    // get node coordinates of element
    LINALG::FixedSizeSerialDenseMatrix<3,iel>  xyze;
    DRT::Node** nodes = Nodes();
    for(int inode=0;inode<iel;inode++)
    {
      const double* x = nodes[inode]->X();
      xyze(0,inode)=x[0];
      xyze(1,inode)=x[1];
      xyze(2,inode)=x[2];
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

      shift=2.0/(static_cast<double>(numplanesinele-1));
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
    LINALG::FixedSizeSerialDenseMatrix<iel,1> funct;

    // get the quad9 gaussrule for the in plane integration
    const IntegrationPoints2D  intpoints = getIntegrationPoints2D(intrule_quad_9point);

    // a hex8 element has two levels, the hex20 and hex27 element have three layers to sample
    // (now we allow even more)
    double layershift=0;

    // loop all levels in element
    for(set<int>::const_iterator id = planesinele.begin();id!=planesinele.end() ;++id)
    {
      // reset temporary values
      double ubar=0;
      double vbar=0;
      double wbar=0;
      double pbar=0;

      double usqbar=0;
      double vsqbar=0;
      double wsqbar=0;
      double uvbar =0;
      double uwbar =0;
      double vwbar =0;
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

#ifdef DEBUG
	// check whether this gausspoint is really inside the desired plane
	{
	  double x[3];
	  x[0]=0;
	  x[1]=0;
	  x[2]=0;
	  for(int inode=0;inode<iel;inode++)
	  {
	    x[0]+=funct(inode)*xyze(0,inode);
	    x[1]+=funct(inode)*xyze(1,inode);
	    x[2]+=funct(inode)*xyze(2,inode);
	  }

	  if(abs(x[normdirect]-(*planes)[*id])>2e-9)
	  {
	    dserror("Mixing up element cut planes during integration");
	  }
	}
#endif

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
          int finode=inode*4;

	  ugp += funct(inode)*sol(finode++);
          vgp += funct(inode)*sol(finode++);
	  wgp += funct(inode)*sol(finode++);
	  pgp += funct(inode)*sol(finode  );
	}

	// add contribution to integral

	double dubar  = ugp*fac;
	double dvbar  = vgp*fac;
	double dwbar  = wgp*fac;
	double dpbar  = pgp*fac;

	ubar   += dubar;
	vbar   += dvbar;
	wbar   += dwbar;
	pbar   += dpbar;

	usqbar += ugp*dubar;
	vsqbar += vgp*dvbar;
	wsqbar += wgp*dwbar;
	uvbar  += ugp*dvbar;
	uwbar  += ugp*dwbar;
	vwbar  += vgp*dwbar;
	psqbar += pgp*dpbar;
      } // end loop integration points

      // add increments from this layer to processor local vectors
      (*sumu  )[*id] += ubar;
      (*sumv  )[*id] += vbar;
      (*sumw  )[*id] += wbar;
      (*sump  )[*id] += pbar;

      (*sumsqu)[*id] += usqbar;
      (*sumsqv)[*id] += vsqbar;
      (*sumsqw)[*id] += wsqbar;
      (*sumuv) [*id] += uvbar;
      (*sumuw) [*id] += uwbar;
      (*sumvw) [*id] += vwbar;
      (*sumsqp)[*id] += psqbar;

      // jump to the next layer in the element.
      // in case of an hex8 element, the two coordinates are -1 and 1(+2)
      // for quadratic elements with three sample planes, we have -1,0(+1),1(+2)

      layershift+=2.0/(static_cast<double>(numplanesinele-1));
    }
  }
  else if(distype == DRT::Element::nurbs8 || distype == DRT::Element::nurbs27)
  {
    // get size of planecoords
    int size = planes->size();

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    if(nurbsdis == NULL)
    {
      dserror("we need a nurbs discretisation for nurbs elements\n");
    }

    // get nurbs dis' element numbers
    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

    // use size of planes and mele to determine number of layers
    int numsublayers=(size-1)/nele_x_mele_x_lele[1];

    // get the knotvector itself
    RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    DRT::Node**   nodes = Nodes();

    // get gid, location in the patch
    int gid = Id();

    vector<int> ele_cart_id;

    int npatch = -1;

    knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);
    if(npatch!=0)
    {
      dserror("expected single patch nurbs problem for calculating means");
    }

    // access elements knot span
    std::vector<blitz::Array<double,1> > eleknots(3);
    knots->GetEleKnots(eleknots,gid);

    // aquire weights from nodes
    LINALG::FixedSizeSerialDenseMatrix<iel,1> weights;

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
	=
	dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights(inode) = cp->W();
    }

    // get shapefunctions, compute all visualisation point positions
    LINALG::FixedSizeSerialDenseMatrix<iel,1> nurbs_shape_funct;

    // there's one additional plane for the last element layer
    int endlayer=0;
    if(ele_cart_id[1]!=nele_x_mele_x_lele[1]-1)
    {
      endlayer=numsublayers;
    }
    else
    {
      endlayer=numsublayers+1;
    }

    // loop layers in element
    for(int rr=0;rr<endlayer;++rr)
    {
      // set gauss point coordinates
      blitz::Array<double, 1> gp(3);

      gp(1)=-1.0+rr*2.0/((double)numsublayers);

      // get the quad9 gaussrule for the in plane integration
      const IntegrationPoints2D  intpoints = getIntegrationPoints2D(intrule_quad_9point);

      // reset temporary values
      double ubar=0;
      double vbar=0;
      double wbar=0;
      double pbar=0;

      double usqbar=0;
      double vsqbar=0;
      double wsqbar=0;
      double uvbar =0;
      double uwbar =0;
      double vwbar =0;
      double psqbar=0;


      // start loop over integration points in layer
      for (int iquad=0;iquad<intpoints.nquad;iquad++)
      {

	// get the other gauss point coordinates
	gp(0)=intpoints.qxg[iquad][0];
	gp(2)=intpoints.qxg[iquad][1];

	// compute the shape function values
	DRT::NURBS::UTILS::nurbs_get_3D_funct
	  (nurbs_shape_funct,
	   gp               ,
	   eleknots         ,
	   weights          ,
	   distype          );

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
	  ugp += nurbs_shape_funct(inode)*sol(inode*4  );
          vgp += nurbs_shape_funct(inode)*sol(inode*4+1);
	  wgp += nurbs_shape_funct(inode)*sol(inode*4+2);
	  pgp += nurbs_shape_funct(inode)*sol(inode*4+3);
	}

	// add contribution to integral
	ubar   += ugp*fac;
	vbar   += vgp*fac;
	wbar   += wgp*fac;
	pbar   += pgp*fac;

	usqbar += ugp*ugp*fac;
	vsqbar += vgp*vgp*fac;
	wsqbar += wgp*wgp*fac;
	uvbar  += ugp*vgp*fac;
	uwbar  += ugp*wgp*fac;
	vwbar  += vgp*wgp*fac;
	psqbar += pgp*pgp*fac;
      } // end loop integration points


      // add increments from this layer to processor local vectors
      (*sumu  )[ele_cart_id[1]*numsublayers+rr] += ubar;
      (*sumv  )[ele_cart_id[1]*numsublayers+rr] += vbar;
      (*sumw  )[ele_cart_id[1]*numsublayers+rr] += wbar;
      (*sump  )[ele_cart_id[1]*numsublayers+rr] += pbar;

      (*sumsqu)[ele_cart_id[1]*numsublayers+rr] += usqbar;
      (*sumsqv)[ele_cart_id[1]*numsublayers+rr] += vsqbar;
      (*sumsqw)[ele_cart_id[1]*numsublayers+rr] += wsqbar;
      (*sumuv) [ele_cart_id[1]*numsublayers+rr] += uvbar;
      (*sumuw) [ele_cart_id[1]*numsublayers+rr] += uwbar;
      (*sumvw) [ele_cart_id[1]*numsublayers+rr] += vwbar;
      (*sumsqp)[ele_cart_id[1]*numsublayers+rr] += psqbar;
    }

  }
  else
  {
    dserror("Unknown element type for mean value evaluation\n");
  }
  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_means

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
template<int iel>
void DRT::ELEMENTS::Fluid3::f3_apply_box_filter(
    vector<double>&           myvel,
    double*                   bvel_hat,
    double*                   breystr_hat,
    double*                   bmodeled_stress_grid_scale_hat,
    double&                   volume
    )
{
  // alloc a fixed size array for nodal velocities
  LINALG::FixedSizeSerialDenseMatrix<3,iel>   evel;

  // wrap matrix objects in fixed-size arrays
  LINALG::FixedSizeSerialDenseMatrix<4*iel,1> myvelvec                     (&(myvel[0])                   ,true);
  LINALG::FixedSizeSerialDenseMatrix<3,1>     vel_hat                      (bvel_hat                      ,true);
  LINALG::FixedSizeSerialDenseMatrix<3,3>     reystr_hat                   (breystr_hat                   ,true);
  LINALG::FixedSizeSerialDenseMatrix<3,3>     modeled_stress_grid_scale_hat(bmodeled_stress_grid_scale_hat,true);

  // split velocity and throw away  pressure, insert into element array
  for (int i=0;i<iel;++i)
  {
    int fi =4*i;

    evel(0,i) = myvelvec(fi++);
    evel(1,i) = myvelvec(fi++);
    evel(2,i) = myvelvec(fi  );
  }

   // number of spatial dimensions is always 3
  const int NSD = 3;

  // set element data
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  LINALG::FixedSizeSerialDenseMatrix<iel,  1>  funct;
  LINALG::FixedSizeSerialDenseMatrix<NSD,NSD>  xjm  ;
  LINALG::FixedSizeSerialDenseMatrix<NSD,NSD>  xji  ;
  LINALG::FixedSizeSerialDenseMatrix<NSD,iel>  deriv;

  // get node coordinates of element
  LINALG::FixedSizeSerialDenseMatrix<NSD,iel> xyze;
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      case DRT::Element::tet10:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        dserror("the box filtering operation is only permitted for linear elements\n");
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);

  // shape functions and derivs at element center
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  const double wquad = intpoints_onepoint.qwgt[0];

  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

  // get Jacobian matrix and determinant
  for (int nn=0;nn<NSD;++nn)
  {
    for (int rr=0;rr<NSD;++rr)
    {
      xjm(nn,rr)=deriv(nn,0)*xyze(rr,0);
      for (int mm=1;mm<iel;++mm)
      {
        xjm(nn,rr)+=deriv(nn,mm)*xyze(rr,mm);
      }
    }
  }

  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);

  //
  //             compute global first derivates
  //
  LINALG::FixedSizeSerialDenseMatrix<NSD,iel> derxy;
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

  */
  xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
  xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
  xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
  xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
  xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
  xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
  xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
  xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
  xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

  // compute global derivates
  for (int nn=0;nn<NSD;++nn)
  {
    for (int rr=0;rr<iel;++rr)
    {
      derxy(nn,rr)=deriv(0,rr)*xji(nn,0);
      for (int mm=1;mm<NSD;++mm)
      {
        derxy(nn,rr)+=deriv(mm,rr)*xji(nn,mm);
      }
    }
  }

  // get velocities (n+alpha_F/1,i) at integration point
  //
  //                   +-----
  //       n+af/1       \                  n+af/1
  //    vel      (x) =   +      N (x) * vel
  //                    /        j         j
  //                   +-----
  //                   node j
  //
  LINALG::FixedSizeSerialDenseMatrix<NSD,1> velint;
  for (int rr=0;rr<NSD;++rr)
  {
    velint(rr)=funct(0)*evel(rr,0);
    for (int mm=1;mm<iel;++mm)
    {
      velint(rr)+=funct(mm)*evel(rr,mm);
    }
  }

  // get velocity (n+alpha_F/1,i) derivatives at integration point
  //
  //       n+af/1      +-----  dN (x)
  //   dvel      (x)    \        k         n+af/1
  //   ------------- =   +     ------ * vel
  //        dx          /        dx        k
  //          j        +-----      j
  //                   node k
  //
  // j : direction of derivative x/y/z
  //
  LINALG::FixedSizeSerialDenseMatrix<NSD,NSD> vderxy;
  for (int nn=0;nn<NSD;++nn)
  {
    for (int rr=0;rr<NSD;++rr)
    {
      vderxy(nn,rr)=derxy(rr,0)*evel(nn,0);
      for (int mm=1;mm<iel;++mm)
      {
        vderxy(nn,rr)+=derxy(rr,mm)*evel(nn,mm);
      }
    }
  }

  /*
                            +-     n+af/1          n+af/1    -+
          / h \       1.0   |  dvel_i    (x)   dvel_j    (x)  |
     eps | u   |    = --- * |  ------------- + -------------  |
          \   / ij    2.0   |       dx              dx        |
                            +-        j               i      -+
  */
  LINALG::FixedSizeSerialDenseMatrix<NSD,NSD> epsilon;

  for (int nn=0;nn<NSD;++nn)
  {
    for (int rr=0;rr<NSD;++rr)
    {
      epsilon(nn,rr)=0.5*(vderxy(nn,rr)+vderxy(rr,nn));
    }
  }

  //
  // modeled part of subgrid scale stresses
  //
  /*    +-                                 -+ 1
        |          / h \           / h \    | -         / h \
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+

        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */

  double rateofstrain = 0;

  for(int rr=0;rr<NSD;rr++)
  {
    for(int mm=0;mm<NSD;mm++)
    {
      rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
    }
  }
  rateofstrain *= 2.0;
  rateofstrain = sqrt(rateofstrain);


  //--------------------------------------------------
  // one point integrations

  // determine contribution to patch volume
  volume = wquad*det;

  for (int rr=0;rr<NSD;++rr)
  {
    double temp=velint(rr)*volume;

    // add contribution to integral over velocities
    vel_hat(rr) += temp;

    // add contribution to integral over reynolds stresses
    for (int nn=0;nn<NSD;++nn)
    {
      reystr_hat(rr,nn) += temp*velint(nn);
    }
  }

  // add contribution to integral over the modeled part of subgrid
  // scale stresses
  double rateofstrain_volume=rateofstrain*volume;
  for (int rr=0;rr<NSD;++rr)
  {
    for (int nn=0;nn<NSD;++nn)
    {
      modeled_stress_grid_scale_hat(rr,nn) += rateofstrain_volume * epsilon(rr,nn);
    }
  }

  return;
} // DRT::ELEMENTS::Fluid3::f3_apply_box_filter

//----------------------------------------------------------------------
// Calculate the quantities LijMij and MijMij, to compare the influence
// of the modeled and resolved stress tensor --- from this relation, Cs
// will be computed
//----------------------------------------------------------------------
template<int iel>
void DRT::ELEMENTS::Fluid3::f3_calc_smag_const_LijMij_and_MijMij(
  RCP<Epetra_MultiVector>& filtered_vel                       ,
  RCP<Epetra_MultiVector>& col_filtered_reynoldsstress        ,
  RCP<Epetra_MultiVector>& col_filtered_modeled_subgrid_stress,
  double&                  LijMij,
  double&                  MijMij,
  double&                  center)
{

  LINALG::FixedSizeSerialDenseMatrix<3,iel> evel_hat                            ;
  LINALG::FixedSizeSerialDenseMatrix<9,iel> ereynoldsstress_hat                 ;
  LINALG::FixedSizeSerialDenseMatrix<9,iel> efiltered_modeled_subgrid_stress_hat;

  for (int nn=0;nn<iel;++nn)
  {
    int lid = (Nodes()[nn])->LID();

    for (int dimi=0;dimi<3;++dimi)
    {
      evel_hat(dimi,nn) = (*((*filtered_vel)(dimi)))[lid];

      for (int dimj=0;dimj<3;++dimj)
      {
        int index=3*dimi+dimj;

        ereynoldsstress_hat(index,nn)
          =(*((*col_filtered_reynoldsstress)(index)))[lid];

        efiltered_modeled_subgrid_stress_hat (index,nn)
          =(*((*col_filtered_modeled_subgrid_stress)(index)))[lid];
      }
    }
  }

  // set element data
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  LINALG::FixedSizeSerialDenseMatrix<iel,1> funct;
  LINALG::FixedSizeSerialDenseMatrix<3,iel> deriv;

  //this will be the y-coordinate of a point in the element interior
  center = 0;

  // get node coordinates of element
  LINALG::FixedSizeSerialDenseMatrix<3,iel> xyze;
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];

    center+=xyze(1,inode);
  }
  center/=iel;


  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
      case DRT::Element::tet10:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points --- i.e. the midpoint
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];

  // shape functions and derivs at element center
  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  LINALG::FixedSizeSerialDenseMatrix<3,3> xjm;
  // get Jacobian matrix and its determinant
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      xjm(nn,rr)=deriv(nn,0)*xyze(rr,0);
      for (int mm=1;mm<iel;++mm)
      {
        xjm(nn,rr)+=deriv(nn,mm)*xyze(rr,mm);
      }
    }
  }
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);

  //
  //             compute global first derivates
  //
  LINALG::FixedSizeSerialDenseMatrix<3,iel> derxy;
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

  */
  LINALG::FixedSizeSerialDenseMatrix<3,3> xji;
  xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
  xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
  xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
  xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
  xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
  xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
  xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
  xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
  xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

  // compute global derivates
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<iel;++rr)
    {
      derxy(nn,rr)=deriv(0,rr)*xji(nn,0);
      for (int mm=1;mm<3;++mm)
      {
        derxy(nn,rr)+=deriv(mm,rr)*xji(nn,mm);
      }
    }
  }

  // get velocities (n+alpha_F/1,i) at integration point
  //
  //                   +-----
  //     ^ n+af/1       \                ^ n+af/1
  //    vel      (x) =   +      N (x) * vel
  //                     /        j         j
  //                    +-----
  //                    node j
  //
  LINALG::FixedSizeSerialDenseMatrix<3,1> velint_hat;
  for (int rr=0;rr<3;++rr)
  {
    velint_hat(rr)=funct(0)*evel_hat(rr,0);
    for (int mm=1;mm<iel;++mm)
    {
      velint_hat(rr)+=funct(mm)*evel_hat(rr,mm);
    }
  }

  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //     ^ n+af/1      +-----  dN (x)
  //   dvel      (x)    \        k       ^ n+af/1
  //   ------------- =   +     ------ * vel
  //       dx           /        dx        k
  //         j         +-----      j
  //                   node k
  //
  // j : direction of derivative x/y/z
  //
  LINALG::FixedSizeSerialDenseMatrix<3,3>  vderxy_hat;

  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      vderxy_hat(nn,rr)=derxy(rr,0)*evel_hat(nn,0);
      for (int mm=1;mm<iel;++mm)
      {
        vderxy_hat(nn,rr)+=derxy(rr,mm)*evel_hat(nn,mm);
      }
    }
  }

  // get filtered reynolds stress (n+alpha_F/1,i) at integration point
  //
  //                        +-----
  //        ^   n+af/1       \                   ^   n+af/1
  //    restress      (x) =   +      N (x) * restress
  //            ij           /        k              k, ij
  //                        +-----
  //                        node k
  //
  LINALG::FixedSizeSerialDenseMatrix<3,3> restress_hat;

  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      int index = 3*nn+rr;
      restress_hat(nn,rr)=funct(0)*ereynoldsstress_hat(index,0);

      for (int mm=1;mm<iel;++mm)
      {
        restress_hat(nn,rr)+=funct(mm)*ereynoldsstress_hat(index,mm);
      }
    }
  }

  //  restress_hat = blitz::sum(funct(k)*ereynoldsstress_hat(i,j,k),k);

  // get filtered modeled subgrid stress (n+alpha_F/1,i) at integration point
  //
  //
  //                   ^                   n+af/1
  //    filtered_modeled_subgrid_stress_hat      (x) =
  //                                       ij
  //
  //            +-----
  //             \                              ^                   n+af/1
  //          =   +      N (x) * filtered_modeled_subgrid_stress_hat
  //             /        k                                         k, ij
  //            +-----
  //            node k
  //
  //
  LINALG::FixedSizeSerialDenseMatrix<3,3> filtered_modeled_subgrid_stress_hat;
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      int index = 3*nn+rr;
      filtered_modeled_subgrid_stress_hat(nn,rr)=funct(0)*efiltered_modeled_subgrid_stress_hat(index,0);

      for (int mm=1;mm<iel;++mm)
      {
        filtered_modeled_subgrid_stress_hat(nn,rr)+=funct(mm)*efiltered_modeled_subgrid_stress_hat(index,mm);
      }
    }
  }

  //filtered_modeled_subgrid_stress_hat = blitz::sum(funct(k)*efiltered_modeled_subgrid_stress_hat(i,j,k),k);


  /*
                            +-   ^ n+af/1        ^   n+af/1    -+
      ^   / h \       1.0   |  dvel_i    (x)   dvel_j      (x)  |
     eps | u   |    = --- * |  ------------- + ---------------  |
          \   / ij    2.0   |       dx              dx          |
                            +-        j               i        -+
  */

  LINALG::FixedSizeSerialDenseMatrix<3,3>  epsilon_hat;
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      epsilon_hat(nn,rr) = 0.5 * ( vderxy_hat(nn,rr) + vderxy_hat(rr,nn) );
    }
  }

  //
  // modeled part of subtestfilter scale stresses
  //
  /*    +-                                 -+ 1
        |      ^   / h \       ^   / h \    | -     ^   / h \
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+

        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */

  double rateofstrain_hat = 0;

  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      rateofstrain_hat += epsilon_hat(rr,mm)*epsilon_hat(rr,mm);
    }
  }
  rateofstrain_hat *= 2.0;
  rateofstrain_hat = sqrt(rateofstrain_hat);

  LINALG::FixedSizeSerialDenseMatrix<3,3> L_ij;
  LINALG::FixedSizeSerialDenseMatrix<3,3> M_ij;

  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      L_ij(rr,mm) = restress_hat(rr,mm) - velint_hat(rr)*velint_hat(mm);
    }
  }

  // this is sqrt(3)
  double filterwidthratio = 1.73;

  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      M_ij(rr,mm) = filtered_modeled_subgrid_stress_hat(rr,mm)
        -
        filterwidthratio*filterwidthratio*rateofstrain_hat*epsilon_hat(rr,mm);
    }
  }

  LijMij =0;
  MijMij =0;
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      LijMij += L_ij(rr,mm)*M_ij(rr,mm);
      MijMij += M_ij(rr,mm)*M_ij(rr,mm);
    }
  }

  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_smag_const_LijMij_and_MijMij


//
// check for higher order derivatives for shape functions
//
bool DRT::ELEMENTS::Fluid3::isHigherOrderElement(
  const DRT::Element::DiscretizationType  distype) const
{
  bool hoel = true;
  switch (distype)
  {
  case hex8: case hex20: case hex27: case tet10: case wedge15: case nurbs8: case nurbs27:
    hoel = true;
    break;
  case tet4: case wedge6: case pyramid5: //!!!TODO:  wedge und pyramid have 2nd derivatives!!!!!!!!!!!!!!!!!!!!!!!!
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
int DRT::ELEMENTS::Fluid3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
