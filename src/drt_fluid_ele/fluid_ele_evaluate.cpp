/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of actions of fluid element

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/


#include "fluid_ele_factory.H"

#include "fluid_ele.H"
#include "fluid_ele_xwall.H"
#include "fluid_ele_action.H"
#include "fluid_ele_evaluate_utils.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_interface.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_immersed.H"
#include "fluid_ele_parameter_poro.H"
#include "fluid_ele_parameter_xfem.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_parameter_intface.H"
#include "fluid_ele_tds.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_inpar/inpar_fluid.H"

#include "../drt_lib/drt_utils.H"

#include "../drt_opti/topopt_fluidAdjoint3_interface.H"
#include "../drt_opti/topopt_fluidAdjoint3_impl_parameter.H"

#include "../drt_lib/drt_condition_utils.H"


/*
  Depending on the type of action and the element type (tet, hex etc.),
  the elements allocate common static arrays.

  */

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter                         |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidType::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const FLD::Action action = DRT::INPUT::get<FLD::Action>(p, "action");

  if (action == FLD::set_general_fluid_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementGeneralFluidParameter(p, dis.Comm().MyPID());
  }
  else if (action == FLD::set_time_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterTimInt* fldpara =
        DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
    fldpara->SetElementTimeParameter(p);
  }
  else if (action == FLD::set_turbulence_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementTurbulenceParameters(p);
  }
  else if (action == FLD::set_loma_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementLomaParameter(p);
  }
  else if (action == FLD::set_topopt_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementTopoptParameter(p);
  }
  else if (action == FLD::set_general_adjoint_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldpara =
        DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
    fldpara->SetElementGeneralAdjointParameter(p);
  }
  else if (action == FLD::set_adjoint_time_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldpara =
        DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
    fldpara->SetElementAdjointTimeParameter(p);
  }
  else if (action == FLD::set_two_phase_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementTwoPhaseParameter(p);
  }
  else if (action == FLD::set_general_fluid_xfem_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterXFEM* fldpara =
        DRT::ELEMENTS::FluidEleParameterXFEM::Instance();

    fldpara->SetElementGeneralFluidParameter(p, dis.Comm().MyPID());
    fldpara->SetElementTurbulenceParameters(p);
    fldpara->SetElementXFEMParameter(p, dis.Comm().MyPID());
  }

  return;
}


/*----------------------------------------------------------------------*
|  evaluate the element (public)                            g.bau 03/07|
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();

  // get space dimensions
  const int nsd = DRT::UTILS::getDimension(Shape());

  // switch between different physical types as used below
  std::string impltype = "std";
  switch (params.get<int>("Physical Type", INPAR::FLUID::incompressible))
  {
    case INPAR::FLUID::loma:
      impltype = "loma";
      break;
  }

  DRT::ELEMENTS::FluidImmersed* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersed*>(this);
  if (immersedele)  // not a standard immersed element and the node row maps don't know it's nodes
    impltype = "std_immersed";


  DRT::ELEMENTS::FluidXWall* xwallele = dynamic_cast<DRT::ELEMENTS::FluidXWall*>(this);
  if (xwallele)  // not a xwall element and the node row maps don't know it's nodes
    impltype = "xw";

  switch (act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->Evaluate(
              this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for evaluation of off-diagonal matrix block for monolithic
    // low-Mach-number solver
    //-----------------------------------------------------------------------
    case FLD::calc_loma_mono_odblock:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "loma")
          ->Evaluate(this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2,
              elevec3, true);
    }
    break;
    case FLD::calc_turbscatra_statistics:
    {
      if (nsd == 3)
      {
        // do nothing if you do not own this element
        if (this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocity, pressure, and scalar from global
          // distributed vectors
          // --------------------------------------------------
          // velocity/pressure and scalar values (n+1)
          Teuchos::RCP<const Epetra_Vector> velnp =
              discretization.GetState("u and p (n+1,converged)");
          Teuchos::RCP<const Epetra_Vector> scanp =
              discretization.GetState("scalar (n+1,converged)");
          if (velnp == Teuchos::null || scanp == Teuchos::null)
            dserror("Cannot get state vectors 'velnp' and/or 'scanp'");

          // extract local values from the global vectors
          std::vector<double> myvelpre(lm.size());
          std::vector<double> mysca(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp, myvelpre, lm);
          DRT::UTILS::ExtractMyValues(*scanp, mysca, lm);

          // integrate mean values
          const DiscretizationType distype = this->Shape();

          switch (distype)
          {
            case DRT::Element::hex8:
            {
              FLD::f3_calc_scatra_means<8>(this, discretization, myvelpre, mysca, params);
              break;
            }
            case DRT::Element::hex20:
            {
              FLD::f3_calc_scatra_means<20>(this, discretization, myvelpre, mysca, params);
              break;
            }
            case DRT::Element::hex27:
            {
              FLD::f3_calc_scatra_means<27>(this, discretization, myvelpre, mysca, params);
              break;
            }
            default:
            {
              dserror("Unknown element type for turbulent passive scalar mean value evaluation\n");
            }
          }
        }
      }  // end if (nsd == 3)
      else
        dserror("action 'calc_turbscatra_statistics' is a 3D specific action");
    }
    break;
    case FLD::calc_loma_statistics:
    {
      if (nsd == 3)
      {
        // do nothing if you do not own this element
        if (this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocity, pressure, and temperature from
          // global distributed vectors
          // --------------------------------------------------
          // velocity/pressure and scalar values (n+1)
          Teuchos::RCP<const Epetra_Vector> velnp =
              discretization.GetState("u and p (n+1,converged)");
          Teuchos::RCP<const Epetra_Vector> scanp =
              discretization.GetState("scalar (n+1,converged)");
          if (velnp == Teuchos::null || scanp == Teuchos::null)
            dserror("Cannot get state vectors 'velnp' and/or 'scanp'");

          // extract local values from global vectors
          std::vector<double> myvelpre(lm.size());
          std::vector<double> mysca(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp, myvelpre, lm);
          DRT::UTILS::ExtractMyValues(*scanp, mysca, lm);

          // get factor for equation of state
          const double eosfac = params.get<double>("eos factor", 100000.0 / 287.0);

          // integrate mean values
          const DiscretizationType distype = this->Shape();

          switch (distype)
          {
            case DRT::Element::hex8:
            {
              FLD::f3_calc_loma_means<8>(this, discretization, myvelpre, mysca, params, eosfac);
              break;
            }
            case DRT::Element::hex20:
            {
              FLD::f3_calc_loma_means<20>(this, discretization, myvelpre, mysca, params, eosfac);
              break;
            }
            case DRT::Element::hex27:
            {
              FLD::f3_calc_loma_means<27>(this, discretization, myvelpre, mysca, params, eosfac);
              break;
            }
            default:
            {
              dserror("Unknown element type for low-Mach-number mean value evaluation\n");
            }
          }
        }
      }  // end if (nsd == 3)
      else
        dserror("action 'calc_loma_statistics' is a 3D specific action");
    }
    break;
    case FLD::calc_fluid_box_filter:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        int nen = 0;
        if (distype == DRT::Element::hex8)
          nen = 8;
        else if (distype == DRT::Element::tet4)
          nen = 4;
        else
          dserror("not supported");

        // --------------------------------------------------
        // extract velocity and pressure from global
        // distributed vectors
        // --------------------------------------------------
        // velocity, pressure and temperature values (most recent
        // intermediate solution, i.e. n+alphaF for genalpha
        // and n+1 for one-step-theta)
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
        if (vel == Teuchos::null) dserror("Cannot get state vectors 'vel'");
        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);

        std::vector<double> tmp_temp(lm.size());
        std::vector<double> mytemp(nen);
        double thermpress = 0.0;
        // pointer to class FluidEleParameter (access to the general parameter)
        DRT::ELEMENTS::FluidEleParameterStd* fldpara =
            DRT::ELEMENTS::FluidEleParameterStd::Instance();
        if (fldpara->PhysicalType() == INPAR::FLUID::loma)
        {
          Teuchos::RCP<const Epetra_Vector> temp = discretization.GetState("T (trial)");
          if (temp == Teuchos::null) dserror("Cannot get state vectors 'temp'");
          DRT::UTILS::ExtractMyValues(*temp, tmp_temp, lm);

          for (int i = 0; i < nen; i++) mytemp[i] = tmp_temp[nsd + (i * (nsd + 1))];

          // get thermodynamic pressure
          thermpress = params.get<double>("thermpress");
        }
        // initialize the contribution of this element to the patch volume to zero
        double volume_contribution = 0.0;

        // initialize the contributions of this element to the filtered scalar quantities
        double dens_hat = 0.0;
        double dens_strainrate_hat = 0.0;
        // get pointers for vector quantities
        Teuchos::RCP<std::vector<double>> vel_hat =
            params.get<Teuchos::RCP<std::vector<double>>>("vel_hat");
        Teuchos::RCP<std::vector<double>> densvel_hat =
            params.get<Teuchos::RCP<std::vector<double>>>("densvel_hat");
        Teuchos::RCP<std::vector<std::vector<double>>> reynoldsstress_hat =
            params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("reynoldsstress_hat");
        Teuchos::RCP<std::vector<std::vector<double>>> modeled_subgrid_stress =
            params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("modeled_subgrid_stress");
        // Vreman
        double expression_hat = 0.0;
        double alpha2_hat = 0.0;
        Teuchos::RCP<std::vector<std::vector<double>>> strainrate_hat =
            params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("strainrate_hat");
        Teuchos::RCP<std::vector<std::vector<double>>> alphaij_hat =
            params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("alphaij_hat");
        // integrate the convolution with the box filter function for this element
        // the results are assembled onto the *_hat arrays
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            FLD::f3_apply_box_filter<8>(this, fldpara, myvel, mytemp, thermpress, vel_hat,
                densvel_hat, reynoldsstress_hat, modeled_subgrid_stress, volume_contribution,
                dens_hat, dens_strainrate_hat, expression_hat, alpha2_hat, strainrate_hat,
                alphaij_hat);
            break;
          }
          case DRT::Element::tet4:
          {
            FLD::f3_apply_box_filter<4>(this, fldpara, myvel, mytemp, thermpress, vel_hat,
                densvel_hat, reynoldsstress_hat, modeled_subgrid_stress, volume_contribution,
                dens_hat, dens_strainrate_hat, expression_hat, alpha2_hat, strainrate_hat,
                alphaij_hat);
            break;
          }
          default:
          {
            dserror("Unknown element type for box filter application\n");
          }
        }

        // hand down the volume contribution to the time integration algorithm
        params.set<double>("volume_contribution", volume_contribution);
        // as well as the filtered scalar quantities
        params.set<double>("dens_hat", dens_hat);
        params.set<double>("dens_strainrate_hat", dens_strainrate_hat);

        params.set<double>("expression_hat", expression_hat);
        params.set<double>("alpha2_hat", alpha2_hat);

      }  // end if (nsd == 3)
      else
        dserror("action 'calc_fluid_box_filter' is 3D specific action");
    }
    break;
    case FLD::calc_smagorinsky_const:
    {
      if (nsd == 3)
      {
        Teuchos::RCP<Epetra_MultiVector> col_filtered_vel =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_vel");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_reynoldsstress =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_reynoldsstress");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_modeled_subgrid_stress =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_modeled_subgrid_stress");

        // pointer to class FluidEleParameter (access to the general parameter)
        DRT::ELEMENTS::FluidEleParameterStd* fldpara =
            DRT::ELEMENTS::FluidEleParameterStd::Instance();
        // add potential loma specific vectors
        Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel = Teuchos::null;
        Teuchos::RCP<Epetra_Vector> col_filtered_dens = Teuchos::null;
        Teuchos::RCP<Epetra_Vector> col_filtered_dens_strainrate = Teuchos::null;
        if (fldpara->PhysicalType() == INPAR::FLUID::loma)
        {
          col_filtered_dens_vel =
              params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_dens_vel");
          col_filtered_dens = params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_dens");
          col_filtered_dens_strainrate =
              params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_dens_strainrate");
        }

        double LijMij = 0.0;
        double MijMij = 0.0;
        double CI_numerator = 0.0;
        double CI_denominator = 0.0;
        double xcenter = 0.0;
        double ycenter = 0.0;
        double zcenter = 0.0;

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            FLD::f3_calc_smag_const_LijMij_and_MijMij<8>(this, fldpara, col_filtered_vel,
                col_filtered_reynoldsstress, col_filtered_modeled_subgrid_stress,
                col_filtered_dens_vel, col_filtered_dens, col_filtered_dens_strainrate, LijMij,
                MijMij, CI_numerator, CI_denominator, xcenter, ycenter, zcenter);
            break;
          }
          case DRT::Element::tet4:
          {
            FLD::f3_calc_smag_const_LijMij_and_MijMij<4>(this, fldpara, col_filtered_vel,
                col_filtered_reynoldsstress, col_filtered_modeled_subgrid_stress,
                col_filtered_dens_vel, col_filtered_dens, col_filtered_dens_strainrate, LijMij,
                MijMij, CI_numerator, CI_denominator, xcenter, ycenter, zcenter);
            break;
          }
          default:
          {
            dserror("Unknown element type for box filter application\n");
          }
        }

        double Cs_delta_sq = 0.0;
        double Ci_delta_sq = 0.0;
        // set Cs_delta_sq without averaging (only clipping)
        if (abs(MijMij) < 1E-16)
          Cs_delta_sq = 0.0;
        else
          Cs_delta_sq = 0.5 * LijMij / MijMij;
        if (Cs_delta_sq < 0.0)
        {
          Cs_delta_sq = 0.0;
        }

        if (fldpara->PhysicalType() == INPAR::FLUID::loma)
        {
          // and set Ci_delta_sq without averaging (only clipping) for loma
          if (abs(CI_denominator) < 1E-16)
            Ci_delta_sq = 0.0;
          else
            Ci_delta_sq = 0.5 * CI_numerator / CI_denominator;
          if (Ci_delta_sq < 0.0)
          {
            Ci_delta_sq = 0.0;
          }
        }

        // set all values in parameter list
        params.set<double>("ele_Cs_delta_sq", Cs_delta_sq);
        params.set<double>("ele_Ci_delta_sq", Ci_delta_sq);
        params.set<double>("LijMij", LijMij);
        params.set<double>("MijMij", MijMij);
        params.set<double>("CI_numerator", CI_numerator);
        params.set<double>("CI_denominator", CI_denominator);
        params.set<double>("xcenter", xcenter);
        params.set<double>("ycenter", ycenter);
        params.set<double>("zcenter", zcenter);
      }  // end if(nsd == 3)
      else
        dserror("action 'calc_smagorinsky_const' is a 3D specific action");
    }
    break;
    case FLD::calc_vreman_const:
    {
      if (nsd == 3)
      {
        Teuchos::RCP<Epetra_MultiVector> col_filtered_strainrate =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_strainrate");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_alphaij =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_alphaij");
        // pointer to class FluidEleParameter (access to the general parameter)
        Teuchos::RCP<Epetra_Vector> col_filtered_expression = Teuchos::null;
        Teuchos::RCP<Epetra_Vector> col_filtered_alpha2 = Teuchos::null;
        col_filtered_expression =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_expression");
        col_filtered_alpha2 = params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_alpha2");

        double cv_numerator = 0.0;
        double cv_denominator = 0.0;
        double volume = 0.0;
        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            FLD::f3_calc_vreman_const<8>(this, col_filtered_strainrate, col_filtered_alphaij,
                col_filtered_expression, col_filtered_alpha2, cv_numerator, cv_denominator, volume);
            break;
          }
          case DRT::Element::tet4:
          {
            FLD::f3_calc_vreman_const<4>(this, col_filtered_strainrate, col_filtered_alphaij,
                col_filtered_expression, col_filtered_alpha2, cv_numerator, cv_denominator, volume);
            break;
          }
          default:
          {
            dserror("Unknown element type for dynamic vreman application\n");
          }
        }

        elevec1(0) = cv_numerator;
        elevec1(1) = cv_denominator;
      }  // end if(nsd == 3)
      else
        dserror("action 'calc_vreman_const' is a 3D specific action");
    }
    break;
    case FLD::calc_fluid_genalpha_update_for_subscales:
    {
      // time update for time-dependent subgrid-scales
      const double dt = params.get<double>("dt");
      const double gamma = params.get<double>("gamma");
      this->TDS()->Update(dt, gamma);
    }
    break;
    case FLD::calc_model_params_mfsubgr_scales:
    {
      if (nsd == 3)
      {
        // velocity values
        Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
        // fine-scale velocity values
        Teuchos::RCP<const Epetra_Vector> fsvelnp = discretization.GetState("fsvelnp");
        if (velnp == Teuchos::null or fsvelnp == Teuchos::null)
        {
          dserror("Cannot get state vectors");
        }

        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp, myvel, lm);
        std::vector<double> myfsvel(lm.size());
        DRT::UTILS::ExtractMyValues(*fsvelnp, myfsvel, lm);

        // pointer to class FluidEleParameter (access to the general parameter)
        DRT::ELEMENTS::FluidEleParameterStd* fldpara =
            DRT::ELEMENTS::FluidEleParameterStd::Instance();

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            // don't store values of ghosted elements
            if (this->Owner() == discretization.Comm().MyPID())
            {
              FLD::f3_get_mf_params<8, 3, DRT::Element::hex8>(
                  this, fldpara, params, mat, myvel, myfsvel);
            }
            break;
          }
          default:
          {
            dserror("Unknown element type\n");
          }
        }
      }
      else
        dserror("%i D elements does not support calculation of model parameters", nsd);
    }
    break;
    case FLD::calc_mean_Cai:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        int nen = 0;
        if (distype == DRT::Element::hex8)
          nen = 8;
        else
          dserror("not supported");

        // velocity values
        // renamed to "velaf" to be consistent i fluidimplicitintegration.cpp (krank 12/13)
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velaf");
        // scalar values
        Teuchos::RCP<const Epetra_Vector> sca = discretization.GetState("scalar");
        if (vel == Teuchos::null or sca == Teuchos::null)
        {
          dserror("Cannot get state vectors");
        }
        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
        std::vector<double> tmp_sca(lm.size());
        std::vector<double> mysca(nen);
        DRT::UTILS::ExtractMyValues(*sca, tmp_sca, lm);
        for (int i = 0; i < nen; i++) mysca[i] = tmp_sca[nsd + (i * (nsd + 1))];
        // get thermodynamic pressure
        double thermpress = params.get<double>("thermpress at n+alpha_F/n+1", 0.0);

        // pointer to class FluidEleParameter
        DRT::ELEMENTS::FluidEleParameterStd* fldpara =
            DRT::ELEMENTS::FluidEleParameterStd::Instance();

        double Cai = 0.0;
        double vol = 0.0;

        bool is_inflow_ele = false;

        std::vector<DRT::Condition*> myinflowcond;

        // check whether all nodes have a unique inflow condition
        DRT::UTILS::FindElementConditions(this, "TurbulentInflowSection", myinflowcond);
        if (myinflowcond.size() > 1) dserror("More than one inflow condition on one node!");

        if (myinflowcond.size() == 1) is_inflow_ele = true;

        // exclude elemenets of inflow section
        if (not is_inflow_ele)
        {
          switch (distype)
          {
            case DRT::Element::hex8:
            {
              FLD::f3_get_mf_nwc<8, 3, DRT::Element::hex8>(
                  this, fldpara, Cai, vol, myvel, mysca, thermpress);
              break;
            }
            default:
            {
              dserror("Unknown element type\n");
            }
          }
        }

        // hand down the Cai and volume contribution to the time integration algorithm
        params.set<double>("Cai_int", Cai);
        params.set<double>("ele_vol", vol);
      }
      else
        dserror("%i D elements does not support calculation of mean Cai", nsd);
    }
    break;
    case FLD::set_mean_Cai:
    {
      // pointer to class FluidEleParameter
      DRT::ELEMENTS::FluidEleParameterStd* fldpara =
          DRT::ELEMENTS::FluidEleParameterStd::Instance();
      fldpara->SetCsgsPhi(params.get<double>("meanCai"));
    }
    break;
    case FLD::calc_node_normal:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
          case DRT::Element::hex27:
          {
            FLD::ElementNodeNormal<DRT::Element::hex27>(this, params, discretization, lm, elevec1);
            break;
          }
          case DRT::Element::hex20:
          {
            FLD::ElementNodeNormal<DRT::Element::hex20>(this, params, discretization, lm, elevec1);
            break;
          }
          case DRT::Element::hex8:
          {
            FLD::ElementNodeNormal<DRT::Element::hex8>(this, params, discretization, lm, elevec1);
            break;
          }
          case DRT::Element::tet4:
          {
            FLD::ElementNodeNormal<DRT::Element::tet4>(this, params, discretization, lm, elevec1);
            break;
          }
          case DRT::Element::tet10:
          {
            FLD::ElementNodeNormal<DRT::Element::tet10>(this, params, discretization, lm, elevec1);
            break;
          }
          default:
          {
            dserror("Unknown element type for shape function integration\n");
          }
        }
      }
      else
        dserror(
            "action 'calculate node normal' should also work in 2D, but 2D elements are not"
            " added to the template yet. Also it is not tested");
      break;
    }
    case FLD::calc_div_u:
    case FLD::calc_mass_matrix:
    case FLD::calc_fluid_error:
    case FLD::calc_dissipation:
    case FLD::integrate_shape:
    case FLD::calc_divop:
    case FLD::interpolate_velgrad_to_given_point:
    case FLD::interpolate_velocity_to_given_point_immersed:
    case FLD::search_immersed_boundary_elements:
    case FLD::least_squares_matrix_rhs_immersed_boundary:
    case FLD::calc_artificial_velocity_divergence:
    case FLD::interpolate_velocity_to_given_point:
    case FLD::interpolate_pressure_to_given_point:
    case FLD::correct_immersed_fluid_bound_vel:
    case FLD::calc_turbulence_statistics:
    case FLD::xwall_l2_projection:
    case FLD::xwall_calc_mk:
    case FLD::tauw_via_gradient:
    case FLD::velgradient_projection:
    case FLD::presgradient_projection:
    case FLD::calc_velgrad_ele_center:
    case FLD::calc_dt_via_cfl:
    case FLD::calc_mass_flow_periodic_hill:
    case FLD::reset_immersed_ele:
    case FLD::update_immersed_information:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->EvaluateService(
              this, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::set_general_fluid_parameter:
    case FLD::set_time_parameter:
    case FLD::set_turbulence_parameter:
    case FLD::set_loma_parameter:
    case FLD::set_topopt_parameter:
    case FLD::set_general_adjoint_parameter:
    case FLD::set_adjoint_time_parameter:
    case FLD::set_two_phase_parameter:
      //    case FLD::calc_adjoint_neumann: // this is done by the surface elements
      break;
    //-----------------------------------------------------------------------
    // adjoint implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_adjoint_systemmat_and_residual:
    {
      return DRT::ELEMENTS::FluidAdjoint3ImplInterface::Impl(Shape())->Evaluate(
          this, discretization, lm, params, mat, elemat1, elevec1);
      break;
    }
    default:
      dserror("Unknown type of action '%i' for Fluid", act);
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::Fluid::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For fluid elements, the integration   |
 |  integration of volume Neumann conditions (body forces) takes place  |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 | pre-evaluation of FluidIntFaceType class (public)        schott Jun14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidIntFaceType::PreEvaluate(DRT::Discretization& dis,
    Teuchos::ParameterList& p, Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const FLD::Action action = DRT::INPUT::get<FLD::Action>(p, "action");

  if (action == FLD::set_general_face_fluid_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterIntFace* fldintfacepara =
        DRT::ELEMENTS::FluidEleParameterIntFace::Instance();
    fldintfacepara->SetFaceGeneralFluidParameter(p, dis.Comm().MyPID());
  }
  else if (action == FLD::set_general_face_xfem_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterIntFace* fldintfacepara =
        DRT::ELEMENTS::FluidEleParameterIntFace::Instance();
    fldintfacepara->SetFaceGeneralXFEMParameter(p, dis.Comm().MyPID());
  }
  else
    dserror("unknown action type for FluidIntFaceType::PreEvaluate");

  return;
}
