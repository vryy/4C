/*!
\file combust3_evaluate.cpp
\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>

#include "combust3.H"
#include "combust3_sysmat.H"
#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "combust_flamefront.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xdofmapcreation_combust.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>


// converts a std::string into an Action for this element
DRT::ELEMENTS::Combust3::ActionType DRT::ELEMENTS::Combust3::convertStringToActionType(
              const std::string& action) const
{
  DRT::ELEMENTS::Combust3::ActionType act = Combust3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = Combust3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Combust3::calc_fluid_stationary_systemmat_and_residual;
  else if (action == "calc_fluid_beltrami_error")
    act = Combust3::calc_fluid_beltrami_error;
  else if (action == "calc_nitsche_error")
    act = Combust3::calc_nitsche_error;
  else if (action == "calc_void_fraction")
    act = Combust3::calc_void_fraction;
  else if (action == "calc_fluid_box_filter")
    act = Combust3::calc_fluid_box_filter;
  else if (action == "calc_smagorinsky_const")
    act = Combust3::calc_smagorinsky_const;
  else if (action == "store_xfem_info")
    act = Combust3::store_xfem_info;
  else if (action == "get_density")
    act = Combust3::get_density;
  else if (action == "integrate_Shapefunction")
    act = Combust3::integrate_shapefunction;
  else if (action == "reset")
    act = Combust3::reset;
  else if (action == "set_standard_mode")
    act = Combust3::set_standard_mode;
  else
    dserror("Unknown type of action for Combust3");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           g.bau 03/07 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3::Evaluate(Teuchos::ParameterList& params,
                                     DRT::Discretization&      discretization,
                                     std::vector<int>&         lm,
                                     Epetra_SerialDenseMatrix& elemat1,
                                     Epetra_SerialDenseMatrix&,
                                     Epetra_SerialDenseVector& elevec1,
                                     Epetra_SerialDenseVector&,
                                     Epetra_SerialDenseVector&)
{
  // get the action required
  const std::string action(params.get<std::string>("action","none"));
  const DRT::ELEMENTS::Combust3::ActionType act = convertStringToActionType(action);

  // get the list of materials
  const Teuchos::RCP<MAT::Material> material = Material();

  switch(act)
  {
    case get_density:
    {
      std::cout << "Warning! The density is set to 1.0!" << std::endl;
      params.set("density", 1.0);
    }
    break;
    case reset:
    {
      // Reset all information and make element unusable (e.g. it can't answer the numdof question
      // anymore). This way, one can see, if all information is generated correctly or whether
      // something was left from the last nonlinear iteration.
      eleDofManager_ = Teuchos::null;
      ih_ = NULL;
      intersected_ = false;
      splited_ = false;
      bisected_ = false;
      epetra_phinp_ = NULL;
      gradphi_ = NULL;
      gradphi2_ = NULL;
      curvature_ = NULL;
    }
    break;
    case set_standard_mode:
    {
      standard_mode_ = true;
      // reset element dof manager if present
      eleDofManager_ = Teuchos::null;
      ih_ = NULL;
      intersected_ = false;
      splited_ = false;
      bisected_ = false;
      epetra_phinp_ = NULL;
      gradphi_ = NULL;
      gradphi2_ = NULL;
      curvature_ = NULL;
    }
    break;
    case store_xfem_info:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - store_xfem_info");

      // now the element can answer how many (XFEM) dofs it has
      standard_mode_ = false;

      Teuchos::RCP<COMBUST::FlameFront> flamefront = params.get< Teuchos::RCP< COMBUST::FlameFront > >("flamefront",Teuchos::null);

      //--------------------------------------------------
      // find out whether an element is intersected or not
      //--------------------------------------------------
      // remark: initialization call of fluid time integration scheme will also end up here: The initial
      //         flame front has not been incorporated into the fluid field -> no XFEM dofs, yet!
      this->intersected_ = false;
      this->splited_ = false;
      this->bisected_ = false;

      if (flamefront != Teuchos::null) // not the initial call
      {
        // set element variables
        ih_ = flamefront->InterfaceHandle().get();
        epetra_phinp_ = flamefront->Phinp().get();
        if (flamefront->GradPhi()!=Teuchos::null) gradphi_ = flamefront->GradPhi().get();
        if (flamefront->GradPhi2()!=Teuchos::null) gradphi2_ = flamefront->GradPhi2().get();
        if (flamefront->Curvature()!=Teuchos::null) curvature_ = flamefront->Curvature().get();

        this->intersected_ = ih_->ElementIntersected(this->Id());
        this->splited_ = ih_->ElementSplit(this->Id());
        this->bisected_ = ih_->ElementBisected(this->Id());
      }
      else
      {
        // set element variables
        ih_ = NULL;
        epetra_phinp_ = NULL;
        gradphi_ = NULL;
        gradphi2_ = NULL;
        curvature_ = NULL;
      }

      //----------------------------------
      // hand over global dofs to elements
      //----------------------------------
      // get access to global dofman
      const Teuchos::RCP<const XFEM::DofManager> globaldofman = params.get< Teuchos::RCP< XFEM::DofManager > >("dofmanager");

      // create an empty element ansatz map to be filled in the following
      std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_filled;
      // create an empty element ansatz map
      const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty;
      //-------------------------------------------------------------------------
      // build element dof manager according to global dof manager
      // remark: this procedure is closely related to XFEM::createDofMapCombust()
      //-------------------------------------------------------------------------
      // add node dofs, but no element dofs (stress unknowns) for this element
      eleDofManager_ = Teuchos::rcp(new XFEM::ElementDofManager(*this, element_ansatz_empty, *globaldofman));
    }
    break;
    case calc_fluid_systemmat_and_residual:
    case calc_fluid_stationary_systemmat_and_residual:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - calc_fluid_systemmat_and_residual");

      // do no calculation, if not needed
      if (lm.empty()) break;

      // parameter for type of linearization
      const bool newton = params.get<bool>("include reactive terms for linearisation",false);

      // general parameters for two-phase flow and premixed combustion problems
      const INPAR::COMBUST::CombustionType combusttype   = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");
      const INPAR::COMBUST::VelocityJumpType veljumptype = DRT::INPUT::get<INPAR::COMBUST::VelocityJumpType>(params, "veljumptype");
      const INPAR::COMBUST::FluxJumpType fluxjumptype    = DRT::INPUT::get<INPAR::COMBUST::FluxJumpType>(params, "fluxjumptype");
      const INPAR::COMBUST::SmoothGradPhi smoothgradphi  = DRT::INPUT::get<INPAR::COMBUST::SmoothGradPhi>(params, "smoothgradphi");

      const double nitschevel = params.get<double>("nitschevel");
      const double nitschepres = params.get<double>("nitschepres");

      // parameters for premixed combustion problems
      const double flamespeed = params.get<double>("flamespeed");
      const double marksteinlength = params.get<double>("marksteinlength");

      // parameters for two-phase flow problems with surface tension
      // type of surface tension approximation
      const INPAR::COMBUST::SurfaceTensionApprox surftensapprox = DRT::INPUT::get<INPAR::COMBUST::SurfaceTensionApprox>(params, "surftensapprox");
      const bool connected_interface = params.get<bool>("connected_interface");
      const bool smoothed_boundary_integration = params.get<bool>("smoothed_bound_integration");
      const double variablesurftens = params.get<double>("variablesurftens");
      const double second_deriv = params.get<bool>("second_deriv");

      // parameters for Nitsche boundary terms
      const bool nitsche_convflux = params.get<bool>("nitsche_convflux");
      const bool nitsche_convstab = params.get<bool>("nitsche_convstab");
      const bool nitsche_convpenalty = params.get<bool>("nitsche_convpenalty");
      const bool nitsche_mass = params.get<bool>("nitsche_mass");
      const INPAR::COMBUST::WeightType weighttype = params.get<INPAR::COMBUST::WeightType>("weighttype");


      // stabilization terms
      bool pstab = false;
      bool supg = false;
      bool graddiv = false;
      const INPAR::FLUID::StabType stabtype = DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(params.sublist("RESIDUAL-BASED STABILIZATION"),"STABTYPE");
      if (stabtype == INPAR::FLUID::stabtype_residualbased)
      {
        pstab = DRT::INPUT::IntegralValue<int>(params.sublist("RESIDUAL-BASED STABILIZATION"),"PSPG");
        supg  = DRT::INPUT::IntegralValue<int>(params.sublist("RESIDUAL-BASED STABILIZATION"),"SUPG");
        graddiv = DRT::INPUT::IntegralValue<int>(params.sublist("RESIDUAL-BASED STABILIZATION"),"GRAD_DIV");
      }
      else if (stabtype == INPAR::FLUID::stabtype_edgebased)
      {
        // do nothing here
      }
      else
        dserror("Unknown stabilization for combustion problems");

      // stabilization parameter
      const INPAR::FLUID::TauType tautype = DRT::INPUT::IntegralValue<INPAR::FLUID::TauType>(params.sublist("RESIDUAL-BASED STABILIZATION"),"DEFINITION_TAU");
      // check if stabilization parameter definition can be handled by combust3 element
      if (!(tautype == INPAR::FLUID::tau_taylor_hughes_zarins or
            tautype == INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt or
            tautype == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
            tautype == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt or
            tautype == INPAR::FLUID::tau_taylor_hughes_zarins_scaled or
            tautype == INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt or
            tautype == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall or
            tautype == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt or
            tautype == INPAR::FLUID::tau_shakib_hughes_codina or
            tautype == INPAR::FLUID::tau_shakib_hughes_codina_wo_dt or
            tautype == INPAR::FLUID::tau_codina or
            tautype == INPAR::FLUID::tau_codina_wo_dt))
        dserror("unknown type of stabilization parameter definition");

      // time integration parameters
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

      const double time   = params.get<double>("time");
      const double dt     = params.get<double>("dt");
      const double theta  = params.get<double>("theta");
      const double ga_gamma  = params.get<double>("gamma");
      const double ga_alphaF = params.get<double>("alphaF");
      const double ga_alphaM = params.get<double>("alphaM");

      // instationary formulation
      bool instationary = true;
      if (timealgo == INPAR::FLUID::timeint_stationary) instationary = false;
      // generalized alpha time integration scheme
      bool genalpha = false;
      if (timealgo == INPAR::FLUID::timeint_afgenalpha) genalpha = true;

      // smoothed gradient of phi required (surface tension application)
      bool gradphi = true;
      bool gradphi2 = true;
      if (combusttype == INPAR::COMBUST::combusttype_twophaseflow or
          smoothgradphi == INPAR::COMBUST::smooth_grad_phi_none)
      {
        gradphi = false;
        gradphi2 = false;
      }

      // turbulence models
      const INPAR::FLUID::TurbModelAction turbmodel = params.get<INPAR::FLUID::TurbModelAction>("turbmodel");
      const INPAR::FLUID::FineSubgridVisc fssgv= params.get<INPAR::FLUID::FineSubgridVisc>("fssgv");
      double Cs = 0.0;
      double Csgs = 0.0;
      double alpha = 0.0;
      bool CalcN = true;
      double N = 0.0;
      INPAR::FLUID::RefVelocity refvel = INPAR::FLUID::strainrate;
      INPAR::FLUID::RefLength reflength = INPAR::FLUID::cube_edge;
      double c_nu = 0.0;
      bool near_wall_limit = false;
      bool B_gp = false;
      bool mfs_is_conservative = false;

      bool avm3=false;
      if (fssgv!=INPAR::FLUID::no_fssgv
          or turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
      {
        avm3=true;
        //----------------------------------------------------------------------------//
        // get phinp for avm3 preparation (only in this case, this parameter is set!)
        // normally this vector comes from the flamefront, but there is no flamefront
        // at the moment of AVM3 preparation
        // this is a dummy vector containing only 1.0 values
        if (params.get<Teuchos::RCP<Epetra_Vector> >("phinpcol",Teuchos::null)!=Teuchos::null)
          epetra_phinp_=params.get<Teuchos::RCP<Epetra_Vector> >("phinpcol",Teuchos::null).get();
        if (params.get<Teuchos::RCP<COMBUST::InterfaceHandleCombust> >("interfacehandle_",Teuchos::null)!=Teuchos::null)
          ih_ = params.get<Teuchos::RCP<COMBUST::InterfaceHandleCombust> >("interfacehandle_",Teuchos::null).get();

        // prepare avm3
        if (fssgv!=INPAR::FLUID::no_fssgv)
          Cs = params.get<double>("Cs");

        // prepare avm4
        if (turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
        {
          Teuchos::ParameterList& turbmodelparamsmfs = params.sublist("MULTIFRACTAL SUBGRID SCALES");

          // get parameters of model
          Csgs = turbmodelparamsmfs.get<double>("CSGS");

          if (turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
            alpha = 3.0;
          else
            dserror("Unknown filter type for mfs in combust!");

          CalcN = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"CALC_N");

          N = turbmodelparamsmfs.get<double>("N");

          if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "strainrate")
            refvel = INPAR::FLUID::strainrate;
          else if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "resolved")
            refvel = INPAR::FLUID::resolved;
          else if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "fine_scale")
            refvel = INPAR::FLUID::fine_scale;
          else
            dserror("Unknown velocity!");

          if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "cube_edge")
            reflength = INPAR::FLUID::cube_edge;
          else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "sphere_diameter")
            reflength = INPAR::FLUID::sphere_diameter;
          else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "streamlength")
            reflength = INPAR::FLUID::streamlength;
          else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "gradient_based")
            reflength = INPAR::FLUID::gradient_based;
          else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "metric_tensor")
            reflength = INPAR::FLUID::metric_tensor;
          else
            dserror("Unknown length!");

          c_nu = turbmodelparamsmfs.get<double>("C_NU");

          near_wall_limit = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"NEAR_WALL_LIMIT");

          if (turbmodelparamsmfs.get<std::string>("EVALUATION_B") == "element_center")
            B_gp = false;
          else if (turbmodelparamsmfs.get<std::string>("EVALUATION_B") == "integration_point")
            B_gp = true;
          else
            dserror("Unknown evaluation point!");

          if (turbmodelparamsmfs.get<double>("BETA") > 1.0e-9) dserror("No linearization of mfs. Rhs terms only");

          if (turbmodelparamsmfs.get<std::string>("CONVFORM") == "conservative")
            mfs_is_conservative = true;
          else
            mfs_is_conservative = false;
        }
      }

      // extract local (element level) vectors from global state vectors
      DRT::ELEMENTS::Combust3::MyState mystate(
          discretization,
          lm,
          instationary,
          genalpha,
          gradphi,
          gradphi2,
          avm3,
          this,
          epetra_phinp_,
          gradphi_,
          gradphi2_,
          curvature_);

      XFEM::AssemblyType assembly_type=XFEM::standard_assembly;

      // suppress enrichment dofs
      if(params.get<int>("selectedenrichment") != INPAR::COMBUST::selectedenrichment_none)
        assembly_type = XFEM::ComputeAssemblyType(*eleDofManager_, NumNode(), NodeIds());

      // calculate element coefficient matrix and rhs
      COMBUST::callSysmat(assembly_type,
          this, ih_, *eleDofManager_, mystate, elemat1, elevec1,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha,
          combusttype, flamespeed, marksteinlength, nitschevel, nitschepres, surftensapprox, variablesurftens, second_deriv,
          connected_interface,veljumptype,fluxjumptype,smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);
    }
    break;
    case calc_fluid_beltrami_error:
    {
      dserror("Are you sure you want to use this?");
      // add error only for elements which are not ghosted
      if(this->Owner() == discretization.Comm().MyPID())
      {
        // need current velocity and history vector
        Teuchos::RCP<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
        if (vel_pre_np==Teuchos::null)
          dserror("Cannot get state vectors 'velnp'");

        // extract local values from the global vectors
        std::vector<double> my_vel_pre_np(lm.size());
        DRT::UTILS::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

        // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
        const int numnode = NumNode();
        std::vector<double> myprenp(numnode);
        std::vector<double> myvelnp(3*numnode);

        for (int i=0;i<numnode;++i)
        {
          myvelnp[0+(i*3)]=my_vel_pre_np[0+(i*4)];
          myvelnp[1+(i*3)]=my_vel_pre_np[1+(i*4)];
          myvelnp[2+(i*3)]=my_vel_pre_np[2+(i*4)];

          myprenp[i]=my_vel_pre_np[3+(i*4)];
        }

        // integrate beltrami error
        f3_int_beltrami_err(myvelnp,myprenp,material,params);
      }
    }
    break;
    case calc_nitsche_error:
    {
      TEUCHOS_FUNC_TIME_MONITOR("COMBUST3 - evaluate - calc Nitsche errors");

      // add error only for elements which are not ghosted
      if(this->Owner() == discretization.Comm().MyPID())
      {
        // stationary formulation
        const double time = params.get<double>("time");

        // general parameters for two-phase flow and premixed combustion problems
        const INPAR::COMBUST::CombustionType combusttype   = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");
        const INPAR::COMBUST::SmoothGradPhi smoothgradphi  = DRT::INPUT::get<INPAR::COMBUST::SmoothGradPhi>(params, "smoothgradphi");

        const bool smoothed_boundary_integration = params.get<bool>("smoothed_bound_integration");
        // smoothed gradient of phi required (surface tension application)
        bool gradphi = true;
        bool gradphi2 = true;
        if (combusttype == INPAR::COMBUST::combusttype_twophaseflow or
            smoothgradphi == INPAR::COMBUST::smooth_grad_phi_none)
        {
          gradphi = false;
          gradphi2 = false;
        }

        const INPAR::COMBUST::NitscheError NitscheErrorType = DRT::INPUT::get<INPAR::COMBUST::NitscheError>(params, "Nitsche_Compare_Analyt");

        // extract local (element level) vectors from global state vectors
        DRT::ELEMENTS::Combust3::MyState mystate(
            discretization,
            lm,
            false,
            false,
            gradphi,
            gradphi2,
            false,
            this,
            epetra_phinp_,
            gradphi_,
            gradphi2_,
            curvature_);

        // get assembly type
        const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(
            *eleDofManager_, NumNode(), NodeIds());

        // calculate Nitsche norms
        COMBUST::callNitscheErrors(params, NitscheErrorType, assembly_type, this, ih_, *eleDofManager_,
            mystate, material, time, smoothed_boundary_integration, combusttype);
      }
    }
    break;
    case Combust3::integrate_shapefunction:
    {
      // This case integrates the pressures shape functions. It is needed by the
      // Krylov Subspace Projection.

      // extract local (element level) vectors from global state vectors
      // we only need phi and therefore turn everything else off
      DRT::ELEMENTS::Combust3::MyState mystate(
          discretization,
          lm,
          false,
          true,
          false,
          false,
          false,
          this,
          epetra_phinp_,
          gradphi_,
          gradphi2_,
          curvature_);

      const XFEM::AssemblyType assembly_type = XFEM::ComputeAssemblyType(*eleDofManager_, NumNode(), NodeIds());

      // integrate the shape functions.(elemat is a dummy)
      COMBUST::callIntegrateShape(assembly_type, this, ih_, *eleDofManager_, mystate, elemat1, elevec1);
    }
    break;
    case calc_void_fraction:
    {
      // omitting check for 3-Dimensionality, since combust is implemented in 3D only.

      // do nothing if you do not own this element
      if(this->Owner() == discretization.Comm().MyPID())
      {
        const DiscretizationType distype = this->Shape();
        // integrate mean values
        switch(distype)
        {
        case DRT::Element::hex8:
          calc_volume_fraction(discretization, params);
          break;
        default:
          dserror("Unknown element type for mean value evaluation\n");
          break;
        }
      }
    } // end of case calc_void_fraction
    break;
    default:
      dserror("Unknown type of action for Combust3");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3::EvaluateNeumann(Teuchos::ParameterList& params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

// get optimal gaussrule for discretization type
DRT::UTILS::GaussRule3D DRT::ELEMENTS::Combust3::getOptimalGaussrule(const DiscretizationType& distype)
{
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = DRT::UTILS::intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = DRT::UTILS::intrule_hex_27point;
        break;
    case tet4:
        rule = DRT::UTILS::intrule_tet_4point;
        break;
    case tet10:
        rule = DRT::UTILS::intrule_tet_5point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

/*---------------------------------------------------------------------*
 |  calculate error for beltrami test problem               gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3::f3_int_beltrami_err(
    std::vector<double>&      evelnp,
    std::vector<double>&      eprenp,
    Teuchos::RCP<const MAT::Material> material,
    Teuchos::ParameterList&   params
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
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // set constants for analytical solution
  const double t = params.get("total time",-1.0);
  dsassert (t >= 0.0, "beltrami: no total time for error calculation");

  const double a      = M_PI/4.0;
  const double d      = M_PI/2.0;

  // get viscosity
  double  visc = 0.0;
  // just to be sure - actually this has already been checked in Combust3::Evaluate, before
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // use material MAT 3 (first in the material list) for Beltrami flow
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(3);
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  std::cout << "material MAT 3 is used for Beltrami flow" << std::endl;
  const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  visc = mat->Viscosity();

  double         preint;
  std::vector<double> velint  (3);
  std::vector<double> xint    (3);

  std::vector<double> u       (3);

  double         deltap;
  std::vector<double> deltavel(3);

  // gaussian points
  const DRT::UTILS::GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints3D  intpoints(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // declaration of gauss point variables
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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
    LINALG::Matrix<NSD,NSD>    xjm;

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
 | Calculate spatial mean values for channel flow       rasthofer 06/11|
 |                                                      DA wichmann    |
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3::calc_volume_fraction(
  DRT::Discretization&       discretization,
  Teuchos::ParameterList&    params
  )
{
  // set element data
  const DiscretizationType distype = this->Shape();
  if (distype != DRT::Element::hex8)
    dserror("This is for hex8 only");

  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");

  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  Teuchos::RCP<std::vector<double> > planes = params.get<Teuchos::RCP<std::vector<double> > >("coordinate vector for hom. planes");

  // a map to link a material to its index in the vectors
  Teuchos::RCP<std::map<int,int> > matidtoindex = params.get<Teuchos::RCP<std::map<int,int> > >("map materialid to index");
  Teuchos::RCP<const MAT::Material> material = this->Material();
  int matidplus  = -1;
  int matidminus = -1;

  if (INPAR::MAT::m_matlist == material->MaterialType())
  {
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
    matidplus  = matlist->MatID(0);
    matidminus = matlist->MatID(1);
  }
  else
    dserror("COMBUST elements need a material list");


  // get the pointers to the solution vectors
  std::vector<Teuchos::RCP<std::vector<double> > >& sumvol = *(params.get<std::vector<Teuchos::RCP<std::vector<double> > >* >("element volume"));

  // get node coordinates of element
  LINALG::Matrix<3,8>  xyze;
  DRT::Node** nodes = Nodes();
  for(size_t inode=0; inode<8; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode)=x[0];
    xyze(1,inode)=x[1];
    xyze(2,inode)=x[2];
  }

  double min = xyze(normdirect,0);
  double max = xyze(normdirect,0);

  // set maximum and minimum value in wall normal direction
  for(size_t inode=0; inode<8; inode++)
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

  // determine the id of the homogeneous planes intersecting this element
  int eleplaneid = -1;
  for(size_t nplane=0; nplane<planes->size(); ++nplane)
  {
    if (min-2e-9 < (*planes)[nplane] and min+2e-9 > (*planes)[nplane])
    {
      eleplaneid = nplane;
      break;
    }
  }

  if (eleplaneid == 0)
  {
    // this is an element of the lowest element layer. Increase the counter
    // in order to compute the total number of elements in one layer
    int* count = params.get<int*>("count processed elements");

    (*count)++;
  }

  //-------------------------------------------------
  // loop over all integration cells and evaluate on each cell

  const GEO::DomainIntCells& domainintcells = ih_->ElementDomainIntCells(this->Id());

  for (GEO::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
  {
    const bool   indomainplus = cell->getDomainPlus();
    const double cellvol      = cell->VolumeInPhysicalDomain();
    const int    curmatid     = (indomainplus) ? matidplus : matidminus;

    const int phaseindex = matidtoindex->at(curmatid);

    (*(sumvol[phaseindex]))[eleplaneid] += cellvol;
  }
  return;
} // DRT::ELEMENTS::Combust3::calc_volume_fraction

