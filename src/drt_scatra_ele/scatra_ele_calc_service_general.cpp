/*----------------------------------------------------------------------*/
/*!

\brief Internal implementation of ScaTra element

\level 1

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"
#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_turbulence.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_immersed_problem/immersed_field_exchange_manager.H"  // for cell migration

#include "../drt_mat/scatra_mat_multiscale.H"

#include "../drt_volmortar/volmortar_shape.H"

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  //(for now) only first dof set considered
  const std::vector<int>& lm = la[0].lm_;
  // determine and evaluate action
  switch (action)
  {
    // calculate global mass matrix
    case SCATRA::calc_mass_matrix:
    {
      // integration points and weights
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        // evaluate values of shape functions and domain integration factor at current integration
        // point
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // loop over dofs
        for (int k = 0; k < numdofpernode_; ++k) CalcMatMass(elemat1_epetra, k, fac, 1.);
      }  // loop over integration points

      break;
    }

    // calculate time derivative for time value t_0
    case SCATRA::calc_initial_time_deriv:
    {
      // calculate matrix and rhs
      CalcInitialTimeDerivative(ele, elemat1_epetra, elevec1_epetra, params, discretization, la);
      break;
    }

    case SCATRA::integrate_shape_functions:
    {
      // calculate integral of shape functions
      const Epetra_IntSerialDenseVector& dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
      IntegrateShapeFunctions(ele, elevec1_epetra, dofids);

      break;
    }

    case SCATRA::integrate_weighted_scalar:
    {
      // calculate integral of scalar
      IntegrateWeightedScalar(params, ele, elevec1_epetra);

      break;
    }

    case SCATRA::calc_flux_domain:
    {
      // get number of dofset associated with velocity related dofs
      const int ndsvel = params.get<int>("ndsvel");

      // get velocity values at nodes
      const Teuchos::RCP<const Epetra_Vector> convel =
          discretization.GetState(ndsvel, "convective velocity field");
      const Teuchos::RCP<const Epetra_Vector> vel =
          discretization.GetState(ndsvel, "velocity field");

      // safety check
      if (convel == Teuchos::null or vel == Teuchos::null) dserror("Cannot get state vector");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel(nsd_ * nen_, -1);
      for (unsigned inode = 0; inode < nen_; ++inode)
        for (unsigned idim = 0; idim < nsd_; ++idim)
          lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // extract local values of (convective) velocity field from global state vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*convel, econvelnp_, lmvel);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*vel, evelnp_, lmvel);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      rotsymmpbc_->RotateMyValuesIfNecessary(econvelnp_);
      rotsymmpbc_->RotateMyValuesIfNecessary(evelnp_);

      // need current values of transported scalar
      // -> extract local values from global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      // access control parameter for flux calculation
      INPAR::SCATRA::FluxType fluxtype = scatrapara_->CalcFluxDomain();
      Teuchos::RCP<std::vector<int>> writefluxids = scatrapara_->WriteFluxIds();

      // we always get an 3D flux vector for each node
      LINALG::Matrix<3, nen_> eflux(true);

      // do a loop for systems of transported scalars
      for (std::vector<int>::iterator it = writefluxids->begin(); it != writefluxids->end(); ++it)
      {
        int k = (*it) - 1;
        // calculate flux vectors for actual scalar
        eflux.Clear();
        CalculateFlux(eflux, ele, fluxtype, k);
        // assembly
        for (unsigned inode = 0; inode < nen_; inode++)
        {
          const int fvi = inode * numdofpernode_ + k;
          elevec1_epetra[fvi] += eflux(0, inode);
          elevec2_epetra[fvi] += eflux(1, inode);
          elevec3_epetra[fvi] += eflux(2, inode);
        }
      }  // loop over numscal

      break;
    }

    case SCATRA::calc_total_and_mean_scalars:
    {
      // get flag for inverting
      bool inverting = params.get<bool>("inverting");

      // need current scalar vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      // calculate scalars and domain integral
      CalculateScalars(ele, elevec1_epetra, inverting);

      break;
    }

    case SCATRA::calc_mean_scalar_time_derivatives:
    {
      CalculateScalarTimeDerivatives(discretization, lm, elevec1_epetra);
      break;
    }

    // calculate filtered fields for calculation of turbulent Prandtl number
    // required for dynamic Smagorinsky model in scatra
    case SCATRA::calc_scatra_box_filter:
    {
      if (nsd_ == 3)
        CalcBoxFilter(ele, params, discretization, la);
      else
        dserror("action 'calc_scatra_box_filter' is 3D specific action");

      break;
    }

    // calculate turbulent prandtl number of dynamic Smagorinsky model
    case SCATRA::calc_turbulent_prandtl_number:
    {
      if (nsd_ == 3)
      {
        // get required quantities, set in dynamic Smagorinsky class
        Teuchos::RCP<Epetra_MultiVector> col_filtered_vel =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_vel");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_dens_vel");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel_temp =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_dens_vel_temp");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_rateofstrain_temp =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_dens_rateofstrain_temp");
        Teuchos::RCP<Epetra_Vector> col_filtered_temp =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_temp");
        Teuchos::RCP<Epetra_Vector> col_filtered_dens =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_dens");
        Teuchos::RCP<Epetra_Vector> col_filtered_dens_temp =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_dens_temp");

        // initialize variables to calculate
        double LkMk = 0.0;
        double MkMk = 0.0;
        double xcenter = 0.0;
        double ycenter = 0.0;
        double zcenter = 0.0;

        // calculate LkMk and MkMk
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            scatra_calc_smag_const_LkMk_and_MkMk(col_filtered_vel, col_filtered_dens_vel,
                col_filtered_dens_vel_temp, col_filtered_dens_rateofstrain_temp, col_filtered_temp,
                col_filtered_dens, col_filtered_dens_temp, LkMk, MkMk, xcenter, ycenter, zcenter,
                ele);
            break;
          }
          default:
          {
            dserror("Unknown element type for box filter application\n");
          }
        }

        // set Prt without averaging (only clipping)
        // calculate inverse of turbulent Prandtl number times (C_s*Delta)^2
        double inv_Prt = 0.0;
        if (abs(MkMk) < 1E-16)
        {
          // std::cout << "warning: abs(MkMk) < 1E-16 -> set inverse of turbulent Prandtl number to
          // zero!"  << std::endl;
          inv_Prt = 0.0;
        }
        else
          inv_Prt = LkMk / MkMk;
        if (inv_Prt < 0.0) inv_Prt = 0.0;

        // set all values in parameter list
        params.set<double>("LkMk", LkMk);
        params.set<double>("MkMk", MkMk);
        params.set<double>("xcenter", xcenter);
        params.set<double>("ycenter", ycenter);
        params.set<double>("zcenter", zcenter);
        params.set<double>("ele_Prt", inv_Prt);
      }
      else
        dserror("action 'calc_turbulent_prandtl_number' is a 3D specific action");

      break;
    }

    case SCATRA::calc_vreman_scatra:
    {
      if (nsd_ == 3)
      {
        Teuchos::RCP<Epetra_MultiVector> col_filtered_phi =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_phi");
        Teuchos::RCP<Epetra_Vector> col_filtered_phi2 =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_phi2");
        Teuchos::RCP<Epetra_Vector> col_filtered_phiexpression =
            params.get<Teuchos::RCP<Epetra_Vector>>("col_filtered_phiexpression");
        Teuchos::RCP<Epetra_MultiVector> col_filtered_alphaijsc =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("col_filtered_alphaijsc");

        // initialize variables to calculate
        double dt_numerator = 0.0;
        double dt_denominator = 0.0;

        // calculate LkMk and MkMk
        switch (distype)
        {
          case DRT::Element::hex8:
          {
            scatra_calc_vreman_dt(col_filtered_phi, col_filtered_phi2, col_filtered_phiexpression,
                col_filtered_alphaijsc, dt_numerator, dt_denominator, ele);
            break;
          }
          default:
          {
            dserror("Unknown element type for vreman scatra application\n");
          }
          break;
        }

        elevec1_epetra(0) = dt_numerator;
        elevec1_epetra(1) = dt_denominator;
      }
      else
        dserror("action 'calc_vreman_scatra' is a 3D specific action");


      break;
    }

    // calculate domain integral, i.e., surface area or volume of domain element
    case SCATRA::calc_domain_integral:
    {
      CalcDomainIntegral(ele, elevec1_epetra);

      break;
    }

    // calculate normalized subgrid-diffusivity matrix
    case SCATRA::calc_subgrid_diffusivity_matrix:
    {
      // calculate mass matrix and rhs
      CalcSubgrDiffMatrix(ele, elemat1_epetra);

      break;
    }

    // calculate mean Cai of multifractal subgrid-scale modeling approach
    case SCATRA::calc_mean_Cai:
    {
      // get number of dofset associated with velocity related dofs
      const int ndsvel = params.get<int>("ndsvel");

      // get convective (velocity - mesh displacement) velocity at nodes
      Teuchos::RCP<const Epetra_Vector> convel =
          discretization.GetState(ndsvel, "convective velocity field");
      if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel(nsd_ * nen_, -1);
      for (unsigned inode = 0; inode < nen_; ++inode)
        for (unsigned idim = 0; idim < nsd_; ++idim)
          lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // extract local values of convective velocity field from global state vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*convel, econvelnp_, lmvel);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      rotsymmpbc_->RotateMyValuesIfNecessary(econvelnp_);

      // get phi for material parameters
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      if (turbparams_->TurbModel() != INPAR::FLUID::multifractal_subgrid_scales)
        dserror("Multifractal_Subgrid_Scales expected");

      double Cai = 0.0;
      double vol = 0.0;
      // calculate Cai and volume, do not include elements of potential inflow section
      if (turbparams_->AdaptCsgsPhi() and turbparams_->Nwl() and (not SCATRA::InflowElement(ele)))
      {
        // use one-point Gauss rule to do calculations at the element center
        DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToStabGaussRule<distype>::rule);
        vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints, 0);

        // adopt integration points and weights for gauss point evaluation of B
        if (turbparams_->BD_Gp())
        {
          const DRT::UTILS::IntPointsAndWeights<nsd_ele_> gauss_intpoints(
              SCATRA::DisTypeToOptGaussRule<distype>::rule);
          intpoints = gauss_intpoints;
        }

        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        {
          const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

          // density at t_(n)
          std::vector<double> densn(numscal_, 1.0);
          // density at t_(n+1) or t_(n+alpha_F)
          std::vector<double> densnp(numscal_, 1.0);
          // density at t_(n+alpha_M)
          std::vector<double> densam(numscal_, 1.0);

          // diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific
          // heat) in case of loma

          diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManager(numscal_));

          // fluid viscosity
          double visc(0.0);

          // set internal variables
          SetInternalVariablesForMatAndRHS();

          // get material
          GetMaterialParams(ele, densn, densnp, densam, visc);

          // get velocity at integration point
          LINALG::Matrix<nsd_, 1> convelint(true);
          convelint.Multiply(econvelnp_, funct_);

          // calculate characteristic element length
          double hk = CalcRefLength(vol, convelint);

          // estimate norm of strain rate
          double strainnorm = GetStrainRate(econvelnp_);
          strainnorm /= sqrt(2.0);

          // get Re from strain rate
          double Re_ele_str = strainnorm * hk * hk * densnp[0] / visc;
          if (Re_ele_str < 0.0) dserror("Something went wrong!");
          // ensure positive values
          if (Re_ele_str < 1.0) Re_ele_str = 1.0;

          // calculate corrected Cai
          //           -3/16
          //  =(1 - (Re)   )
          //
          Cai += (1.0 - pow(Re_ele_str, -3.0 / 16.0)) * fac;
        }
      }

      // hand down the Cai and volume contribution to the time integration algorithm
      params.set<double>("Cai_int", Cai);
      params.set<double>("ele_vol", vol);

      break;
    }

    // calculate dissipation introduced by stabilization and turbulence models
    case SCATRA::calc_dissipation:
    {
      CalcDissipation(params, ele, discretization, la);
      break;
    }

    case SCATRA::calc_integr_grad_reac:
    {
      Teuchos::RCP<const Epetra_Vector> dualphi = discretization.GetState("adjoint phi");
      Teuchos::RCP<const Epetra_Vector> psi = discretization.GetState("psi");
      if (dualphi == Teuchos::null || psi == Teuchos::null)
        dserror("Cannot get state vector 'dual phi' or 'psi' in action calc_integr_grad_reac");
      Epetra_SerialDenseVector mydualphi(nen_);
      Epetra_SerialDenseVector mypsi(nen_);
      DRT::UTILS::ExtractMyValues(*dualphi, mydualphi, lm);
      DRT::UTILS::ExtractMyValues(*psi, mypsi, lm);

      // compute mass matrix
      Epetra_SerialDenseMatrix massmat;
      massmat.Shape(nen_, nen_);
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);
      double area = 0.0;
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);
        area += fac;
        CalcMatMass(massmat, 0, fac, 1.0);
      }

      elevec1_epetra.Multiply('N', 'N', 1.0, massmat, mydualphi, 0.0);
      elevec1_epetra.Multiply('N', 'N', 1.0, massmat, mypsi, 1.0);

      break;
    }

    case SCATRA::calc_integr_grad_diff:
    {
      Teuchos::RCP<const Epetra_Vector> dualphi = discretization.GetState("adjoint phi");
      if (dualphi == Teuchos::null)
        dserror("Cannot get state vector 'dual phi' in action calc_integr_grad_reac");
      Epetra_SerialDenseVector mydualphi(nen_);
      DRT::UTILS::ExtractMyValues(*dualphi, mydualphi, lm);

      // compute stiffness matrix
      Epetra_SerialDenseMatrix stiffmat;
      stiffmat.Shape(nen_, nen_);
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);
      double area = 0.0;
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);
        area += fac;
        CalcMatDiff(stiffmat, 0, fac);
      }

      if (diffmanager_->GetIsotropicDiff(0) == 0.0)
        elevec1_epetra.Scale(0.0);
      else
        elevec1_epetra.Multiply(
            'N', 'N', 1.0 / (diffmanager_->GetIsotropicDiff(0)), stiffmat, mydualphi, 0.0);

      break;
    }

    case SCATRA::recon_gradients_at_nodes:
    {
      const INPAR::SCATRA::L2ProjectionSystemType systemtype =
          params.get<INPAR::SCATRA::L2ProjectionSystemType>(
              "l2 proj system", INPAR::SCATRA::l2_proj_system_std);

      // need current scalar vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      //    CalcGradientAtNodes(ele, elemat1_epetra, elevec1_epetra, elevec2_epetra,
      //    elevec3_epetra);
      CalcGradientAtNodes(ele, elemat1_epetra, elemat2_epetra, systemtype);

      break;
    }

    case SCATRA::calc_grad_ele_center:
    {
      // need current scalar vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      CalcGradientEleCenter(ele, elevec1_epetra, elevec2_epetra);

      break;
    }

    case SCATRA::recon_curvature_at_nodes:
    {
      // need current scalar vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

      Teuchos::RCP<const Epetra_Vector> gradphinp_x = discretization.GetState("gradphinp_x");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'gradphinp_x'");

      Teuchos::RCP<const Epetra_Vector> gradphinp_y = discretization.GetState("gradphinp_y");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'gradphinp_y'");

      Teuchos::RCP<const Epetra_Vector> gradphinp_z = discretization.GetState("gradphinp_z");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'gradphinp_z'");

      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      std::vector<double> mygradphinp_x(lm.size());
      DRT::UTILS::ExtractMyValues(*gradphinp_x, mygradphinp_x, lm);

      std::vector<double> mygradphinp_y(lm.size());
      DRT::UTILS::ExtractMyValues(*gradphinp_y, mygradphinp_y, lm);

      std::vector<double> mygradphinp_z(lm.size());
      DRT::UTILS::ExtractMyValues(*gradphinp_z, mygradphinp_z, lm);

      std::vector<LINALG::Matrix<nen_, nsd_>> egradphinp(numscal_);

      // fill all element arrays
      for (unsigned i = 0; i < nen_; ++i)
      {
        for (int k = 0; k < numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          egradphinp[k](i, 0) = mygradphinp_x[k + (i * numdofpernode_)];
          egradphinp[k](i, 1) = mygradphinp_y[k + (i * numdofpernode_)];
          egradphinp[k](i, 2) = mygradphinp_z[k + (i * numdofpernode_)];
        }
      }  // for i

      CalcCurvatureAtNodes(ele, elemat1_epetra, elevec1_epetra, egradphinp);

      break;
    }

    case SCATRA::calc_mass_center_smoothingfunct:
    {
      double interface_thickness = params.get<double>("INTERFACE_THICKNESS_TPF");

      if (numscal_ > 1)
      {
        std::cout << "#############################################################################"
                     "##############################"
                  << std::endl;
        std::cout << "#                                                 WARNING:                   "
                     "                             #"
                  << std::endl;
        std::cout << "# More scalars than the levelset are transported. Mass center calculations "
                     "have NOT been tested for this. #"
                  << std::endl;
        std::cout << "#                                                                            "
                     "                             # "
                  << std::endl;
        std::cout << "#############################################################################"
                     "##############################"
                  << std::endl;
      }
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // need current scalar vector
        // -> extract local values from the global vectors
        Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
        if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

        // calculate momentum vector and volume for element.
        CalculateMomentumAndVolume(ele, elevec1_epetra, interface_thickness);
      }
      break;
    }

    case SCATRA::calc_integr_pat_rhsvec:
    {
      // extract local values from the global vectors w and phi
      Teuchos::RCP<const Epetra_Vector> rhsvec = discretization.GetState("rhsnodebasedvals");
      if (rhsvec == Teuchos::null) dserror("Cannot get state vector 'rhsnodebasedvals' ");
      Epetra_SerialDenseVector rhs(nen_);
      DRT::UTILS::ExtractMyValues(*rhsvec, rhs, lm);

      Epetra_SerialDenseMatrix mass(nen_, nen_);
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);
        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const double v = fac * funct_(vi);  // no density required here

          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            mass(vi, ui) += v * funct_(ui);
          }
        }
      }

      elevec1_epetra.Multiply('N', 'N', 1.0, mass, rhs, 0.0);

      break;
    }

    case SCATRA::calc_error:
    {
      // check if length suffices
      if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

      // need current solution
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      CalErrorComparedToAnalytSolution(ele, params, elevec1_epetra);

      break;
    }

    case SCATRA::calc_immersed_element_source:
    {
      int scalartoprovidwithsource = 0;
      double segregationconst = params.get<double>("segregation_constant");

      // assembly
      for (unsigned inode = 0; inode < nen_; inode++)
      {
        const int fvi = inode * numdofpernode_ + scalartoprovidwithsource;
        elevec1_epetra[fvi] += segregationconst;
      }

      break;
    }
    case SCATRA::calc_immersed_phi_at_given_point:
    {
      // need current solution
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null)
        dserror("failed to get state phinp for action 'calc_immersed_phi_at_given_point'");

      // extract local values from the global vector
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      break;
    }

    case SCATRA::micro_scale_initialize:
    {
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
      {
        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over all Gauss points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
          // initialize micro scale in multi-scale simulations
          Teuchos::rcp_static_cast<MAT::ScatraMatMultiScale>(ele->Material())
              ->Initialize(ele->Id(), iquad);
      }

      break;
    }

    case SCATRA::micro_scale_prepare_time_step:
    case SCATRA::micro_scale_solve:
    {
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
      {
        // extract state variables at element nodes
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(
            *discretization.GetState("phinp"), ephinp_, lm);

        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over all Gauss points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        {
          // evaluate shape functions at Gauss point
          EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

          // evaluate state variables at Gauss point
          SetInternalVariablesForMatAndRHS();

          if (action == SCATRA::micro_scale_prepare_time_step)
            // prepare time step on micro scale
            Teuchos::rcp_static_cast<MAT::ScatraMatMultiScale>(ele->Material())
                ->PrepareTimeStep(iquad, std::vector<double>(1, scatravarmanager_->Phinp(0)));
          else
          {
            // solve micro scale
            std::vector<double> dummy(1, 0.);
            Teuchos::rcp_static_cast<MAT::ScatraMatMultiScale>(ele->Material())
                ->Evaluate(
                    iquad, std::vector<double>(1, scatravarmanager_->Phinp(0)), dummy[0], dummy);
          }
        }
      }

      break;
    }

    case SCATRA::micro_scale_update:
    {
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
      {
        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over all Gauss points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
          // update multi-scale scalar transport material
          Teuchos::rcp_static_cast<MAT::ScatraMatMultiScale>(ele->Material())->Update(iquad);
      }

      break;
    }

    case SCATRA::micro_scale_output:
    {
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
      {
        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over all Gauss points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
          // create output on micro scale
          Teuchos::rcp_static_cast<MAT::ScatraMatMultiScale>(ele->Material())->Output(iquad);
      }

      break;
    }

    case SCATRA::micro_scale_read_restart:
    {
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
      {
        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // loop over all Gauss points
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
          // read restart on micro scale
          Teuchos::rcp_dynamic_cast<MAT::ScatraMatMultiScale>(ele->Material())->ReadRestart(iquad);
      }

      break;
    }
    /*----------------------------------------------------------------------*
     | computes sources and sinks due to stress fiber assembly / disassembly|
     |                                                          rauch 08/16 |
     *----------------------------------------------------------------------*/
    case SCATRA::calc_cell_mechanotransduction:
    {
      const int NumPhiROCK = DRT::Problem::Instance()
                                 ->CellMigrationParams()
                                 .sublist("SCALAR TRANSPORT DOF IDS")
                                 .get<int>("ROCK");
      const int NumPhiActin = DRT::Problem::Instance()
                                  ->CellMigrationParams()
                                  .sublist("SCALAR TRANSPORT DOF IDS")
                                  .get<int>("ACTIN_MONOMER");
      ;

      if (NumPhiROCK > -1 and NumPhiActin > -1)
      {
        // get immersed manager and pointers to rates
        DRT::ImmersedFieldExchangeManager* immersedmanager =
            DRT::ImmersedFieldExchangeManager::Instance();
        Teuchos::RCP<Epetra_MultiVector> rates = immersedmanager->GetPointerToRates();
        Teuchos::RCP<Epetra_MultiVector> ratesActin = immersedmanager->GetPointerToRatesActin();

        if (rates == Teuchos::null) dserror("rates = Teuchos::null");
        if (ratesActin == Teuchos::null) dserror("ratesActin = Teuchos::null");

        double timefacrhs = scatraparatimint_->TimeFacRhs();
        const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
            SCATRA::DisTypeToOptGaussRule<distype>::rule);

        // gp loop
        for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        {
          double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

          // loop over nodes
          for (unsigned inode = 0; inode < nen_; inode++)
          {
            const int fvi = inode * numdofpernode_ + NumPhiROCK;
            const int fviActin = inode * numdofpernode_ + NumPhiActin;
            int mylid = discretization.ElementColMap()->LID(ele->Id());

            if (mylid == -1) dserror("no corresponding LID found for EleGID %i", ele->Id());

            elevec1_epetra[fvi] +=
                funct_(inode) * (((rates->Pointers())[iquad][mylid])) * fac * timefacrhs;
            elevec1_epetra[fviActin] +=
                funct_(inode) * (((ratesActin->Pointers())[iquad][mylid])) * fac * timefacrhs;

          }  // Node loop
        }    // gp loop
      }      // if valid dof ids are provided
      break;
    }
    /*----------------------------------------------------------------------*
     | computes sources and sinks due to cell protrusion        rauch 08/16 |
     *----------------------------------------------------------------------*/
    case SCATRA::calc_cell_growth_sourcesandsinks:
    {
      // get general stuff
      DRT::Problem* globalproblem = DRT::Problem::Instance();    //< global problem instance
      int numNode = ele->NumNode();                              //< number of nodes of FaceElement
      int numdofpernode = ele->NumDofPerNode(*ele->Nodes()[0]);  //< number of dofs at each node
      if (numNode != 4)
        dserror(
            "Error! numofnodes is not equal 4 !\n"
            "So far, this method is only tested for quad4 surface elements.\n"
            "Review this implementation if you want to use other element types !");

      // get one-step-theta parameter from input file
      static const double theta =
          globalproblem->StructuralDynamicParams().sublist("ONESTEPTHETA").get<double>("THETA");
      if (theta != 1.0)
        dserror(
            "THETA for OST-Time integration not equal 1. This case is not implemented, yet.\n"
            "(1-\theta) part of equations must be implemented. ");

      // get branching limit parameter, i.e. spatial reference lenght(2D)/area(3D) parameter,
      // determining if branching is possible or not
      static const double aBr =
          globalproblem->CellMigrationParams().sublist("PROTRUSION MODULE").get<double>("aBr");
      if (aBr < 0.0) dserror("Invalid Parameter a_Br.");
      // get number of filaments per spatial reference lenght(2D)/area(3D) aBr
      static const int num_fil_on_aBr = globalproblem->CellMigrationParams()
                                            .sublist("PROTRUSION MODULE")
                                            .get<int>("NUM_FIL_ON_aBr");
      if (num_fil_on_aBr < 0) dserror("Invalid Parameter NUM_FIL_ON_aBr.");

      // define limit concentration clim for branching (branching if c_barbedends < clim)
      static const double clim = (const double)num_fil_on_aBr / aBr;

      // frequently used parameters
      static const double k_on = globalproblem->CellMigrationParams()
                                     .sublist("PROTRUSION MODULE")
                                     .get<double>("K_ON");  //< actin polymerisation base rate
      static const double k_br =
          globalproblem->CellMigrationParams()
              .sublist("PROTRUSION MODULE")
              .get<double>("K_BR");  //< branching base rate (assumed to be half as fast as filament
                                     // polymerization)
      if (k_on <= 0.0 or k_br <= 0.0)
        dserror("Invalid Parameters K_ON = %f or K_BR = %f in --CELL DYNAMIC/PROTRUSION MODULE !",
            k_on, k_br);

      // get phi position of concentration from input file
      static const int offset = globalproblem->CellMigrationParams()
                                    .sublist("PROTRUSION MODULE")
                                    .get<int>("NUM_SURF_SCALARS");
      static const int numofcape = globalproblem->CellMigrationParams()
                                       .sublist("PROTRUSION MODULE")
                                       .get<int>("NUMDOF_BARBEDENDS") +
                                   offset;
      static const int numofcam = globalproblem->CellMigrationParams()
                                      .sublist("PROTRUSION MODULE")
                                      .get<int>("NUMDOF_ACTIN") +
                                  offset;
      static const int numofcbr = globalproblem->CellMigrationParams()
                                      .sublist("PROTRUSION MODULE")
                                      .get<int>("NUMDOF_BRANCHES") +
                                  offset;
      if (numofcape == -1 or numofcam == -1 or numofcbr == -1 or offset == -1)
        dserror(
            "Define input parameters NUMDOF_BARBEDENDS , NUMDOF_ACTIN, NUMDOF_BRANCHES, and "
            "NUM_SURF_SCALARS in ---CELL DYNAMIC/PROTRUSION MODULE");

      // get the timefac for the right-hand side
      const double timefacrhs = scatraparatimint_->TimeFacRhs();

      // integration points and weights
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      Teuchos::RCP<const Epetra_MultiVector> phinp = discretization.GetState(0, "phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      Teuchos::RCP<const Epetra_MultiVector> phin = discretization.GetState(0, "phin");
      if (phin == Teuchos::null) dserror("Cannot get state vector 'phin'");


      /////////////////////////////////////////////////////////////////////////////////
      // loop over all integration points
      ////////////////////////////////////////////////////////////////////////////////
      for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
      {
        // get coordinates of current integration point in face element coordinate system --> quad4
        // returns fac := gp_weight * jacobian determinant
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

        // concentrations at n+1
        std::vector<double> myphinp(lm.size());
        DRT::UTILS::ExtractMyValues(*phinp, myphinp, lm);
        // concentrations at n
        std::vector<double> myphin(lm.size());
        DRT::UTILS::ExtractMyValues(*phin, myphin, lm);

        // PHI CONCENTRATION CALCULATIONS
        // phinp concentrations at nodes
        LINALG::Matrix<4, 1> camnp_nd(true);  //< nodal actin monomer concentration at n+1
        LINALG::Matrix<4, 1> cbrn_nd(
            true);  //< nodal branching point (related to Arp2/3) concentration at n

        // phin concentrations at nodes
        LINALG::Matrix<4, 1> capen_nd(true);  //< nodal actin filament pointed end concentration
        LINALG::Matrix<4, 1> camn_nd(true);   //< nodal actin monomer concentration

        // write molecule concentrations from phinp in separate Vectors
        // camnp_nd := actin monomer concentration at surface at element nodes (at time n+1)
        // camn_nd  := actin monomer concentratino at surface at element nodes (at time n)
        // cbrnp_nd := Arp2/3 concentration at surfce at element nodes (at time n+1)
        // capen_nd := pointed end (actin filament) concentration at surface at element nodes (at
        // time n)

        for (int node = 0; node < numNode; node++)
        {
          camnp_nd(node) =
              myphinp[node * numdofpernode + numofcam];  //< nodal actin monomer concentation at n+1

          cbrn_nd(node) =
              myphin[node * numdofpernode + numofcbr];  //< nodal branching point concentration at n
          capen_nd(node) = myphin[node * numdofpernode +
                                  numofcape];  //< nodal filament barbed end concentration at n
          camn_nd(node) =
              myphin[node * numdofpernode + numofcam];  //< nodal actin monomer concentration at n
        }                                               // end element node loop

        // interpolate node values to gp
        // concentration at gauss points
        // camnp_phi := actin monomer concentration at surface at gp from phinp (time n+1)
        // camn_phi  := actin monomer concentration at surface at gp from phin (time n)
        // cbrnp_phi := Arp2/3 concentration at surface at gp from phinp (time n+1)
        // capen_phi := pointed end (actin filament) concentration at surface at gp from phin (time
        // n)

        double camnp_gp = 0.0;
        double cbrn_gp = 0.0;
        double capen_gp = 0.0;
        double camn_gp = 0.0;

        for (int i = 0; i < numNode; i++)
        {
          camnp_gp += funct_(i) * camnp_nd(i);

          cbrn_gp += funct_(i) * cbrn_nd(i);
          capen_gp += funct_(i) * capen_nd(i);
          camn_gp += funct_(i) * camn_nd(i);
        }  // end element node loop to calc gp values

        // actin monomer polymerisation possibility
        // double k_poly = k_on * exp(-traction * delta/ k_bT);
        // NOTE: currently no use of traction dependency -> k_poly = k_on
        double k_poly = k_on;

        // branching probability
        // 0 or 1
        // we evaluate the probability with variables from the last time step
        // this is more stable since this way, we can exclude, that we jump
        // between branching and no branching in the current time step
        double p_br = 0.0;
        if (abs(capen_gp) <= 1e-16)
        {
          p_br = 0.0;
        }
        if (capen_gp < 0.0)
        {
          p_br = 0.0;
        }
        if (capen_gp > 0.0 and abs(capen_gp) > 1e-16 and capen_gp < clim)
        {
          p_br = 1.0;
        }

        //      // DEBUG output
        //      std::cout<<"Ele "<<ele->Id()<<"  gp "<<iquad<<"  capen_gp = "<<capen_gp<<std::endl;
        //      std::cout<<"Ele "<<ele->Id()<<"  gp "<<iquad<<"  camnp_gp = "<<camnp_gp<<std::endl;
        //      std::cout<<"Ele "<<ele->Id()<<"  gp "<<iquad<<"  cbrn_gp  = "<<cbrn_gp<<std::endl;
        //      std::cout<<"Ele "<<ele->Id()<<"  gp "<<iquad<<"  p_br     = "<<p_br<<std::endl;

        // polymerized actin monomers (1. part to actin filaments and 2. part to branches)
        // NOTE !!!!! implemented for case THETA == 1 !!!
        double branching_part = p_br * k_br * cbrn_gp * camnp_gp;
        double cams_dot = (-1.0) * (k_poly * capen_gp * camnp_gp + branching_part * 2.0);

        if (capen_gp < 0.0) cams_dot = 0.0;

        double conc_barbed_ends_gp_dot = 0.0;  //< rate of change of barbed end concentration
        if (abs(p_br) > 1e-14)
        {
          conc_barbed_ends_gp_dot = branching_part;
        }
        // concentration of used arp2/3 complexes for nucleation
        double cbr_dot = (-1.0) * conc_barbed_ends_gp_dot;

        // fill elevec1 for assembly
        for (int node = 0; node < numNode; node++)
        {
          elevec1_epetra[node * numdofpernode + numofcam] +=
              funct_(node) * cams_dot * fac * timefacrhs;

          // no actin filament aging and severing included, yet
          if (conc_barbed_ends_gp_dot < 0.0)
          {
            elevec1_epetra[node * numdofpernode + numofcape] += 0.0;
          }
          else
          {
            elevec1_epetra[node * numdofpernode + numofcape] +=
                funct_(node) * conc_barbed_ends_gp_dot * fac * timefacrhs;
          }
          elevec1_epetra[node * numdofpernode + numofcbr] +=
              funct_(node) * cbr_dot * fac * timefacrhs;
        }  // loop node

      }  // end gp loop

      break;
    }

    case SCATRA::calc_heteroreac_mat_and_rhs:
    {
      //--------------------------------------------------------------------------------
      // extract element based or nodal values
      //--------------------------------------------------------------------------------
      ExtractElementAndNodeValues(ele, params, discretization, la);

      for (int idof = 0; idof < numdofpernode_; idof++)
      {
        // no bodyforce
        bodyforce_[idof].Clear();
      }

      CalcHeteroReacMatAndRHS(ele, elemat1_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::transform_real_to_reference_point:
    {
      // init quantities
      LINALG::Matrix<nsd_, 1> x_real;
      for (unsigned int d = 0; d < nsd_; ++d) x_real(d, 0) = params.get<double*>("point")[d];
      xsi_(0, 0) = 0.0;
      for (unsigned int d = 1; d < nsd_; ++d) xsi_(d, 0) = 0.0;
      int count = 0;
      LINALG::Matrix<nsd_, 1> diff;

      // do the Newton loop
      bool inside = true;
      do
      {
        count++;
        EvalShapeFuncAndDerivsInParameterSpace();
        LINALG::Matrix<nsd_, 1> x_eval;
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          for (unsigned int n = 0; n < nen_; ++n) x_eval(d, 0) += funct_(n, 0) * xyze_(d, n);
          x_eval(d, 0) -= x_real(d, 0);
        }
        diff.MultiplyTN(xij_, x_eval);

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          xsi_(d, 0) -= diff(d, 0);
          if (xsi_(d, 0) > 10.0 || xsi_(d, 0) < -10.0) inside = false;
        }
      } while (count < 20 && diff.Norm1() > 1.0e-10 && inside);

      inside = true;
      for (unsigned int d = 0; d < nsd_; ++d)
        if (xsi_(d, 0) > 1.0 || xsi_(d, 0) < -1.0) inside = false;

      double pointarr[nsd_];
      if (!inside)
      {
        for (unsigned int d = 0; d < nsd_; ++d) pointarr[d] = -123.0;
      }
      else
      {
        for (unsigned int d = 0; d < nsd_; ++d) pointarr[d] = xsi_(d, 0);
      }
      params.set<double*>("point", pointarr);
      params.set<bool>("inside", inside);

      break;
    }

    case SCATRA::evaluate_field_in_point:
    {
      for (unsigned int d = 0; d < nsd_; ++d) xsi_(d, 0) = params.get<double*>("point")[d];

      EvalShapeFuncAndDerivsInParameterSpace();

      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

      if (params.get<int>("numscal") > numscal_)
        dserror("you requested the pointvalue of the %d-th scalar but there is only %d scalars",
            params.get<int>("numscal"), numscal_);

      const double value = funct_.Dot(ephinp_[params.get<int>("numscal")]);

      params.set<double>("value", value);

      break;
    }

    default:
    {
      dserror("Not acting on this action. Forgot implementation?");
      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate service routine                                  fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvaluateService(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  // setup
  if (SetupCalc(ele, discretization) == -1) return 0;

  if (scatrapara_->IsAle())
  {
    // get number of dofset associated with displacement related dofs
    const int ndsdisp = params.get<int>("ndsdisp");

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(nsd_ * nen_, -1);
    for (unsigned inode = 0; inode < nen_; ++inode)
      for (unsigned idim = 0; idim < nsd_; ++idim)
        lmdisp[inode * nsd_ + idim] = la[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*dispnp, edispnp_, lmdisp);

    // add nodal displacements to point coordinates
    UpdateNodeCoordinates();
  }
  else
    edispnp_.Clear();

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params, "action");

  // evaluate action
  EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}


/*-------------------------------------------------------------------------*
 | Element reconstruct grad phi, one deg of freedom (for now) winter 04/14 |
 *-------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcGradientEleCenter(
    const DRT::Element* ele, Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2)
{
  if (distype != DRT::Element::hex8 && distype != DRT::Element::tet4 &&
      distype != DRT::Element::quad4 && distype != DRT::Element::tri3)
    dserror("this is currently only implemented for linear elements");
  // get node coordinates

  //  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  //  // Do ALE specific updates if necessary
  //  if (ele->IsAle())
  //  {
  //    LINALG::Matrix<nsd_,nen_> edisp(true);
  //    ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &edisp, NULL, "disp");
  //
  //    // get new node positions of ALE mesh
  //     xyze_ += edisp;
  //  }

  // evaluate shape functions and derivatives element center
  EvalShapeFuncAndDerivsAtEleCenter();

  if (numdofpernode_ != 1) dserror("Not tested for more than 1 dofset for ScaTra");

  for (int k = 0; k < numdofpernode_; k++)
  {
    LINALG::Matrix<nsd_, 1> grad_phi(true);

    grad_phi.Multiply(derxy_, ephinp_[k]);

    // Compute element vectors. For L2-Projection
    for (unsigned i = 0; i < nsd_; ++i)
    {
      elevec1[i * numdofpernode_ + k] = grad_phi(i, 0);
      // elevec2[node_i*numdofpernode_+k] += funct_(node_i) * fac * grad_phi(1,0);
      // elevec3[node_i*numdofpernode_+k] += funct_(node_i) * fac * grad_phi(2,0);
    }

  }  // loop over degrees of freedom

  // get position of element centroid
  LINALG::Matrix<nsd_, 1> x_centroid(true);
  x_centroid.Multiply(xyze_, funct_);
  for (unsigned i = 0; i < nsd_; ++i)
  {
    elevec2[i] = x_centroid(i);
  }

}  // ScaTraEleCalc::CalcGradientEleCenter

/*-------------------------------------------------------------------------*
 | Element reconstruct grad phi, one deg of freedom (for now) winter 04/14 |
 *-------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcGradientAtNodes(const DRT::Element* ele,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // Loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Get integration factor and evaluate shape func and its derivatives at the integration points.
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);
    // Loop over degrees of freedom
    for (int k = 0; k < numdofpernode_; k++)
    {
      LINALG::Matrix<3, 1> grad_phi(true);

      // Get gradient of phi at Gauss point
      for (unsigned node = 0; node < nen_; node++)
        for (unsigned idim = 0; idim < nsd_; idim++)
          grad_phi(idim, 0) += ephinp_[k](node, 0) * derxy_(idim, node);


      // const double eta_smooth=0.0; //Option for smoothing factor, currently not used.
      // diffmanager_->SetIsotropicDiff(eta_smooth,k);

      // Compute element matrix. For L2-projection
      CalcMatDiff(elemat1, k, fac);
      CalcMatMass(elemat1, k, fac, 1.0);

      // Compute element vectors. For L2-Projection
      for (unsigned node_i = 0; node_i < nen_; node_i++)
      {
        elevec1[node_i * numdofpernode_ + k] += funct_(node_i) * fac * grad_phi(0, 0);
        elevec2[node_i * numdofpernode_ + k] += funct_(node_i) * fac * grad_phi(1, 0);
        elevec3[node_i * numdofpernode_ + k] += funct_(node_i) * fac * grad_phi(2, 0);
      }

    }  // loop over degrees of freedom

  }  // loop over integration points

  return;

}  // ScaTraEleCalc::CalcGradientAtNodes


/*-------------------------------------------------------------------------*
 | Element reconstruct grad phi, one deg of freedom (for now) winter 11/15 |
 *-------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcGradientAtNodes(const DRT::Element* ele,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
    const INPAR::SCATRA::L2ProjectionSystemType& systemtype)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // Loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    /* Get integration factor and evaluate solution shape func and its
     * derivatives at the integration points. */
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);
    // evaluate test functions
    switch (systemtype)
    {
      case INPAR::SCATRA::l2_proj_system_std:
      case INPAR::SCATRA::l2_proj_system_lumped:
      {
        dual_funct_.SetView(funct_.A());
        break;
      }
      case INPAR::SCATRA::l2_proj_system_dual:
      {
        if (diffmanager_->GetIsotropicDiff(0) != 0.0)
          dserror(
              "The dual shape function do not support an "
              "additional diffusive term!");

        VOLMORTAR::UTILS::dual_shape_function<distype>(
            dual_funct_, xsi_.A(), *ele, INPAR::VOLMORTAR::dualquad_no_mod);
        break;
      }
    }

    // Loop over degrees of freedom
    for (int k = 0; k < numdofpernode_; k++)
    {
      LINALG::Matrix<3, 1> grad_phi(true);

      // Get gradient of phi at Gauss point
      for (unsigned node = 0; node < nen_; node++)
        for (unsigned idim = 0; idim < nsd_; idim++)
          grad_phi(idim, 0) += ephinp_[k](node, 0) * derxy_(idim, node);

      // const double eta_smooth=0.0; //Option for smoothing factor, currently not used.
      // diffmanager_->SetIsotropicDiff(eta_smooth,k);

      // Compute element matrix. For L2-projection
      CalcMatDiff(elemat1, k, fac);
      CalcMatMass(elemat1, k, fac, 1.0, funct_, dual_funct_);

      switch (systemtype)
      {
        case INPAR::SCATRA::l2_proj_system_dual:
        {
          for (unsigned vi = 0; vi < nen_; ++vi)
          {
            const int fvi = vi * numdofpernode_ + k;
            // loop all columns
            for (unsigned ui = 0; ui < nen_; ++ui)
            {
              const int fui = ui * numdofpernode_ + k;
              // zeroise off-diagonal entries
              elemat1(fvi, fui) *= static_cast<double>((fui == fvi));
            }
          }
          break;
        }
        case INPAR::SCATRA::l2_proj_system_lumped:
        {
          for (unsigned vi = 0; vi < nen_; ++vi)
          {
            const int fvi = vi * numdofpernode_ + k;

            double sum = 0.0;
            // loop all columns
            for (unsigned ui = 0; ui < nen_; ++ui)
            {
              const int fui = ui * numdofpernode_ + k;
              sum += elemat1(fvi, fui);
              // reset
              elemat1(fvi, fui) = 0.0;
            }
            elemat1(fvi, fvi) = sum;
          }
          break;
        }
        default:
          // do nothing
          break;
      }

      // Compute element vectors. For L2-Projection
      for (unsigned node_i = 0; node_i < nen_; node_i++)
      {
        for (unsigned j = 0; j < nsd_; j++)
        {
          elemat2(node_i, numdofpernode_ * k + j) += dual_funct_(node_i) * fac * grad_phi(j, 0);
        }
      }

    }  // loop over degrees of freedom

  }  // loop over integration points

  return;

}  // ScaTraEleCalc::CalcGradientAtNodes


/*-------------------------------------------------------------------------*
 | Element reconstructed curvature                            winter 04/14 |
 *-------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcCurvatureAtNodes(const DRT::Element* ele,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1,
    const std::vector<LINALG::Matrix<nen_, nsd_>>& egradphinp)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // Loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // Get integration factor and evaluate shape func and its derivatives at the integration points.
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

    // Loop over degrees of freedom
    for (int k = 0; k < numdofpernode_; k++)
    {
      // Get gradient of phi at Gauss point
      LINALG::Matrix<3, 1> grad_phi(true);
      grad_phi.Clear();
      for (unsigned rr = 0; rr < nsd_; rr++)
      {
        for (unsigned inode = 0; inode < nen_; inode++)
          grad_phi(rr) += egradphinp[k](inode, rr) * funct_(inode);
      }


      // get second derivatives of phi at Gauss point
      static LINALG::Matrix<9, 1> grad2_phi;
      grad2_phi.Clear();

      for (unsigned inode = 0; inode < nen_; inode++)
      {
        grad2_phi(0) += derxy_(0, inode) * egradphinp[k](inode, 0);  // ,xx
        grad2_phi(1) += derxy_(1, inode) * egradphinp[k](inode, 1);  // ,yy
        grad2_phi(2) += derxy_(2, inode) * egradphinp[k](inode, 2);  // ,zz
        grad2_phi(3) += derxy_(1, inode) * egradphinp[k](inode, 0);  // ,xy
        grad2_phi(4) += derxy_(2, inode) * egradphinp[k](inode, 0);  // ,xz
        grad2_phi(5) += derxy_(2, inode) * egradphinp[k](inode, 1);  // ,yz
        grad2_phi(6) += derxy_(0, inode) * egradphinp[k](inode, 1);  // ,yx
        grad2_phi(7) += derxy_(0, inode) * egradphinp[k](inode, 2);  // ,zx
        grad2_phi(8) += derxy_(1, inode) * egradphinp[k](inode, 2);  // ,zy
      }

      // get curvature at Gauss point
      double curvature_gp = 0.0;

      double grad_phi_norm = grad_phi.Norm2();
      {
        // check norm of normal gradient
        if (fabs(grad_phi_norm < 1.0E-9))
        {
          std::cout << "grad phi is small -> set to 1.0E9" << grad_phi_norm << std::endl;
          // phi gradient too small -> there must be a local max or min in the level-set field
          // set curvature to a large value
          // note: we have 1.0E12 in CalcCurvature Function
          curvature_gp = 1.0E2;
        }
        else
        {
          double val = grad_phi_norm * grad_phi_norm * grad_phi_norm;
          double invval = 1.0 / val;
          curvature_gp = -invval * (grad_phi(0) * grad_phi(0) * grad2_phi(0) +
                                       grad_phi(1) * grad_phi(1) * grad2_phi(1) +
                                       grad_phi(2) * grad_phi(2) * grad2_phi(2)) -
                         invval * (grad_phi(0) * grad_phi(1) * (grad2_phi(3) + grad2_phi(6)) +
                                      grad_phi(0) * grad_phi(2) * (grad2_phi(4) + grad2_phi(7)) +
                                      grad_phi(1) * grad_phi(2) * (grad2_phi(5) + grad2_phi(8))) +
                         1.0 / grad_phi_norm * (grad2_phi(0) + grad2_phi(1) + grad2_phi(2));

          // *-1.0 because of direction of gradient phi vector, which is minus the normal on the
          // interface pointing from + to -
          curvature_gp *= -1.0;
        }
      }

      // check for norm of grad phi almost zero
      if (fabs(grad_phi_norm < 1.0E-9))
      {
        std::cout << "grad phi is small -> set to 1.0E9" << grad_phi_norm << std::endl;
        grad_phi_norm = 1.0;
      }

      //      const double eta_smooth=0.0; //Option for smoothing factor, currently not used.
      //      diffmanager_->SetIsotropicDiff(eta_smooth,k);

      // Compute element matrix.
      CalcMatDiff(elemat1, k, fac);
      CalcMatMass(elemat1, k, fac, 1.0);

      // Compute rhs vector.
      // element matrix and rhs
      for (unsigned vi = 0; vi < nen_; ++vi)  // loop rows  (test functions)
      {
        elevec1(vi) += fac * funct_(vi) * curvature_gp;
      }

    }  // loop over degrees of freedom

  }  // loop over integration points


  return;

}  // ScaTraEleCalc::CalcCurvatureAtNodes


/*---------------------------------------------------------------------*
 | calculate filtered fields for turbulent Prandtl number   fang 02/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcBoxFilter(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // extract scalar values from global vector
  Teuchos::RCP<const Epetra_Vector> scalar = discretization.GetState("scalar");
  if (scalar == Teuchos::null) dserror("Cannot get scalar!");
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*scalar, ephinp_, la[0].lm_);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get convective (velocity - mesh displacement) velocity at nodes
  Teuchos::RCP<const Epetra_Vector> convel =
      discretization.GetState(ndsvel, "convective velocity field");
  if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

  // determine number of velocity related dofs per node
  const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

  // construct location vector for velocity related dofs
  std::vector<int> lmvel(nsd_ * nen_, -1);
  for (unsigned inode = 0; inode < nen_; ++inode)
    for (unsigned idim = 0; idim < nsd_; ++idim)
      lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

  // extract local values of convective velocity field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*convel, evelnp_, lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  rotsymmpbc_->RotateMyValuesIfNecessary(evelnp_);

  // initialize the contribution of this element to the patch volume to zero
  double volume_contribution = 0.0;
  // initialize the contributions of this element to the filtered scalar quantities
  double dens_hat = 0.0;
  double temp_hat = 0.0;
  double dens_temp_hat = 0.0;
  double phi2_hat = 0.0;
  double phiexpression_hat = 0.0;
  // get pointers for vector quantities
  Teuchos::RCP<std::vector<double>> vel_hat =
      params.get<Teuchos::RCP<std::vector<double>>>("vel_hat");
  Teuchos::RCP<std::vector<double>> densvel_hat =
      params.get<Teuchos::RCP<std::vector<double>>>("densvel_hat");
  Teuchos::RCP<std::vector<double>> densveltemp_hat =
      params.get<Teuchos::RCP<std::vector<double>>>("densveltemp_hat");
  Teuchos::RCP<std::vector<double>> densstraintemp_hat =
      params.get<Teuchos::RCP<std::vector<double>>>("densstraintemp_hat");
  Teuchos::RCP<std::vector<double>> phi_hat =
      params.get<Teuchos::RCP<std::vector<double>>>("phi_hat");
  Teuchos::RCP<std::vector<std::vector<double>>> alphaijsc_hat =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("alphaijsc_hat");
  // integrate the convolution with the box filter function for this element
  // the results are assembled onto the *_hat arrays
  switch (distype)
  {
    case DRT::Element::hex8:
    {
      scatra_apply_box_filter(dens_hat, temp_hat, dens_temp_hat, phi2_hat, phiexpression_hat,
          vel_hat, densvel_hat, densveltemp_hat, densstraintemp_hat, phi_hat, alphaijsc_hat,
          volume_contribution, ele, params);

      break;
    }
    default:
    {
      dserror("Unknown element type for box filter application\n");
      break;
    }
  }

  // hand down the volume contribution to the time integration algorithm
  params.set<double>("volume_contribution", volume_contribution);
  // as well as the filtered scalar quantities
  params.set<double>("dens_hat", dens_hat);
  params.set<double>("temp_hat", temp_hat);
  params.set<double>("dens_temp_hat", dens_temp_hat);
  params.set<double>("phi2_hat", phi2_hat);
  params.set<double>("phiexpression_hat", phiexpression_hat);

  return;
}  // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcBoxFilter


/*-----------------------------------------------------------------------------*
 | calculate mass matrix + rhs for initial time derivative calc.     gjb 03/12 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcInitialTimeDerivative(
    DRT::Element* ele,                    //!< current element
    Epetra_SerialDenseMatrix& emat,       //!< element matrix
    Epetra_SerialDenseVector& erhs,       //!< element residual
    Teuchos::ParameterList& params,       //!< parameter list
    DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la       //!< location array
)
{
  // extract relevant quantities from discretization and parameter list
  ExtractElementAndNodeValues(ele, params, discretization, la);

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

  //------------------------------------------------------------------------------------
  // get material parameters and stabilization parameters (evaluation at element center)
  //------------------------------------------------------------------------------------
  // density at t_(n)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // the stabilisation parameters (one per transported scalar)
  std::vector<double> tau(numscal_, 0.0);

  if (not scatrapara_->MatGP() or not scatrapara_->TauGP())
  {
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);

    if (not scatrapara_->TauGP())
    {
      for (int k = 0; k < numscal_; ++k)  // loop of each transported scalar
      {
        // calculation of stabilization parameter at element center
        CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
            reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k],
            scatravarmanager_->ConVel(k), vol);
      }
    }
  }

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  /*----------------------------------------------------------------------*/
  // element integration loop
  /*----------------------------------------------------------------------*/
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    //------------ get values of variables at integration point
    for (int k = 0; k < numscal_; ++k)  // deal with a system of transported scalars
    {
      // get phi at integration point for all scalars
      const double& phiint = scatravarmanager_->Phinp(k);

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      LINALG::Matrix<nen_, 1> conv = scatravarmanager_->Conv(k);

      // velocity divergence required for conservative form
      double vdiv(0.0);
      if (scatrapara_->IsConservative()) GetDivergence(vdiv, evelnp_);

      // diffusive part used in stabilization terms
      LINALG::Matrix<nen_, 1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
      }

      // calculation of stabilization parameter at integration point
      if (scatrapara_->TauGP())
        CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
            reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k],
            scatravarmanager_->ConVel(k), vol);

      const double fac_tau = fac * tau[k];

      //----------------------------------------------------------------
      // element matrix: transient term
      //----------------------------------------------------------------
      // transient term
      CalcMatMass(emat, k, fac, densam[k]);

      //----------------------------------------------------------------
      // element matrix: stabilization of transient term
      //----------------------------------------------------------------
      // the stabilization term is deactivated in CalcInitialTimeDerivative() on time integrator
      // level
      if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
      {
        // subgrid-scale velocity (dummy)
        LINALG::Matrix<nen_, 1> sgconv(true);
        CalcMatMassStab(emat, k, fac_tau, densam[k], densnp[k], sgconv, diff);

        // remove convective stabilization of inertia term
        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;
          erhs(fvi) += fac_tau * densnp[k] * conv(vi) * densnp[k] * phiint;
        }
      }

      // we solve: (w,dc/dt) = rhs
      // whereas the rhs is based on the standard element evaluation routine
      // including contributions resulting from the time discretization.
      // The contribution from the time discretization has to be removed before solving the system:
      CorrectRHSFromCalcRHSLinMass(erhs, k, fac, densnp[k], phiint);
    }  // loop over each scalar k
  }    // integration loop

  // scale element matrix appropriately to be consistent with scaling of global residual vector
  // computed by AssembleMatAndRHS() routine (see CalcInitialTimeDerivative() routine on time
  // integrator level)
  emat.Scale(scatraparatimint_->TimeFacRhs());

  return;
}  // ScaTraEleCalc::CalcInitialTimeDerivative()


/*----------------------------------------------------------------------*
 |  CorrectRHSFromCalcRHSLinMass                             ehrl 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CorrectRHSFromCalcRHSLinMass(
    Epetra_SerialDenseVector& erhs, const int k, const double fac, const double densnp,
    const double phinp)
{
  // fac->-fac to change sign of rhs
  if (scatraparatimint_->IsIncremental())
    CalcRHSLinMass(erhs, k, 0.0, -fac, 0.0, densnp);
  else
    dserror("Must be incremental!");

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate shape functions over domain (private)           gjb 07/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::IntegrateShapeFunctions(
    const DRT::Element* ele, Epetra_SerialDenseVector& elevec1,
    const Epetra_IntSerialDenseVector& dofids)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // safety check
  if (dofids.M() < numdofpernode_) dserror("Dofids vector is too short. Received not enough flags");

  // loop over integration points
  // this order is not efficient since the integration of the shape functions is always the same for
  // all species
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

    // compute integral of shape functions (only for dofid)
    for (int k = 0; k < numdofpernode_; k++)
    {
      if (dofids[k] >= 0)
      {
        for (unsigned node = 0; node < nen_; node++)
        {
          elevec1[node * numdofpernode_ + k] += funct_(node) * fac;
        }
      }
    }

  }  // loop over integration points

  return;

}  // ScaTraEleCalc::IntegrateShapeFunction

/*----------------------------------------------------------------------*
 |  Integrate weighted scalar                               rauch 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::IntegrateWeightedScalar(
    Teuchos::ParameterList& params, const DRT::Element* ele, Epetra_SerialDenseVector& elevec1)
{
  // extract values from parameter list.
  // these values are set in scatra_timint_ost_endoexocytosis.
  const int scalarid =
      params.get<int>("ScalarID");  //!< dof id of integrated scalar (position in result vector)
  const double scalar = params.get<double>("scalar");  //!< scalar value to be integrated
  const double prefac =
      params.get<double>("user defined prefac");  //!< user defined factor to integral

  // get integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  //////////////////////////////////////
  // loop over integration points
  //////////////////////////////////////
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

    // evaluate element right-hand side vector
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numscal_ + scalarid;
      elevec1[fvi] -= fac * prefac * funct_(vi) * scalar;
    }  // loop over nodes

  }  // loop over integration points

  return;
}  // ScatraEleCalc::IntegrateWeightedScalar

/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalculateFlux(LINALG::Matrix<3, nen_>& flux,
    const DRT::Element* ele, const INPAR::SCATRA::FluxType fluxtype, const int k)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // density at t_(n)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // get material parameters (evaluation at element center)
  if (not scatrapara_->MatGP())
  {
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);
  }

  // integration rule
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    // get material parameters (evaluation at integration point)
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc);

    // get velocity at integration point
    LINALG::Matrix<nsd_, 1> velint(true);
    LINALG::Matrix<nsd_, 1> convelint(true);
    velint.Multiply(evelnp_, funct_);
    convelint.Multiply(econvelnp_, funct_);

    // get gradient of scalar at integration point
    LINALG::Matrix<nsd_, 1> gradphi(true);
    gradphi.Multiply(derxy_, ephinp_[k]);

    // allocate and initialize!
    LINALG::Matrix<nsd_, 1> q(true);

    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case INPAR::SCATRA::flux_total:
        // convective flux contribution
        q.Update(densnp[k] * scatravarmanager_->Phinp(k), convelint);

        // no break statement here!
      case INPAR::SCATRA::flux_diffusive:
        // diffusive flux contribution
        q.Update(-(diffmanager_->GetIsotropicDiff(k)), gradphi, 1.0);

        break;
      default:
        dserror("received illegal flag inside flux evaluation for whole domain");
        break;
    };
    // q at integration point

    // integrate and assemble everything into the "flux" vector
    for (unsigned vi = 0; vi < nen_; vi++)
    {
      for (unsigned idim = 0; idim < nsd_; idim++)
      {
        flux(idim, vi) += fac * funct_(vi) * q(idim);
      }  // idim
    }    // vi

  }  // integration loop

  // set zeros for unused space dimensions
  for (unsigned idim = nsd_; idim < 3; idim++)
  {
    for (unsigned vi = 0; vi < nen_; vi++)
    {
      flux(idim, vi) = 0.0;
    }
  }

  return;
}  // ScaTraCalc::CalculateFlux


/*----------------------------------------------------------------------------------------*
 | calculate domain integral, i.e., surface area or volume of domain element   fang 07/15 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcDomainIntegral(
    const DRT::Element* ele,          //!< the element we are dealing with
    Epetra_SerialDenseVector& scalar  //!< result vector for scalar integral to be computed
)
{
  // initialize variable for domain integral
  double domainintegral(0.);

  // get integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // add contribution from current integration point to domain integral
    for (unsigned vi = 0; vi < nen_; ++vi) domainintegral += funct_(vi) * fac;
  }  // loop over integration points

  // write result into result vector
  scalar(0) = domainintegral;

  return;
}  // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcDomainIntegral


/*----------------------------------------------------------------------*
|  calculate scalar(s) and domain integral                     vg 09/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalculateScalars(
    const DRT::Element* ele, Epetra_SerialDenseVector& scalars, const bool inverting)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integrals of (inverted) scalar(s) and domain
    if (inverting)
    {
      for (unsigned i = 0; i < nen_; i++)
      {
        const double fac_funct_i = fac * funct_(i);
        for (int k = 0; k < numdofpernode_; k++)
        {
          if (std::abs(ephinp_[k](i, 0)) > EPS14)
            scalars[k] += fac_funct_i / ephinp_[k](i, 0);
          else
            dserror("Division by zero");
        }
        // for domain volume
        scalars[numdofpernode_] += fac_funct_i;
      }
    }
    else
    {
      for (unsigned i = 0; i < nen_; i++)
      {
        const double fac_funct_i = fac * funct_(i);
        for (int k = 0; k < numdofpernode_; k++)
        {
          scalars[k] += fac_funct_i * ephinp_[k](i, 0);
        }
        // for domain volume
        scalars[numdofpernode_] += fac_funct_i;
      }
    }
  }  // loop over integration points

  return;
}  // ScaTraEleCalc::CalculateScalars


/*----------------------------------------------------------------------*
 | calculate scalar time derivative(s) and domain integral   fang 03/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalculateScalarTimeDerivatives(
    const DRT::Discretization& discretization,  //!< discretization
    const std::vector<int>& lm,                 //!< location vector
    Epetra_SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
)
{
  // extract scalar time derivatives from global state vector
  const Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
  if (phidtnp == Teuchos::null) dserror("Cannot get state vector \"phidtnp\"!");
  static std::vector<LINALG::Matrix<nen_, 1>> ephidtnp(numscal_);
  DRT::UTILS::ExtractMyValues(*phidtnp, ephidtnp, lm);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integrals of scalar time derivatives
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double fac_funct_vi = fac * funct_(vi);

      for (int k = 0; k < numscal_; ++k) scalars(k) += fac_funct_vi * ephidtnp[k](vi);
    }

    // calculate integral of domain
    scalars(numscal_) += fac;
  }  // loop over integration points

  return;
}  // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalculateScalarTimeDerivatives


/*----------------------------------------------------------------------*
|  calculate momentum vector and minus domain integral          mw 06/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalculateMomentumAndVolume(
    const DRT::Element* ele, Epetra_SerialDenseVector& momandvol, const double interface_thickness)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // coordinates of the current integration point.
    std::vector<double> gpcoord(nsd_);
    // levelset function at gaussian point.
    double ephi_gp = 0.0;
    // fac*funct at GP
    double fac_funct = 0.0;

    for (unsigned i = 0; i < nen_; i++)
    {
      // Levelset function (first scalar stored [0]) at gauss point
      ephi_gp += funct_(i) * ephinp_[0](i, 0);

      // Coordinate * shapefunction to get the coordinate value of the gausspoint.
      for (unsigned idim = 0; idim < nsd_; idim++)
      {
        gpcoord[idim] += funct_(i) * (ele->Nodes()[i]->X()[idim]);
      }

      // Summation of fac*funct_ for volume calculation.
      fac_funct += fac * funct_(i);
    }

    double heavyside_epsilon = 1.0;  // plus side

    // Smoothing function
    if (abs(ephi_gp) <= interface_thickness)
    {
      heavyside_epsilon = 0.5 * (1.0 + ephi_gp / interface_thickness +
                                    1.0 / PI * sin(PI * ephi_gp / interface_thickness));
    }
    else if (ephi_gp < interface_thickness)
    {
      heavyside_epsilon = 0.0;  // minus side
    }


    // add momentum vector and volume
    for (unsigned idim = 0; idim < nsd_; idim++)
    {
      momandvol(idim) += gpcoord[idim] * (1.0 - heavyside_epsilon) * fac_funct;
    }

    momandvol(nsd_) += fac_funct * (1.0 - heavyside_epsilon);

  }  // loop over integration points

  return;
}  // ScaTraEleCalc::CalculateMomentumAndVolume


/*----------------------------------------------------------------------*
 | calculate normalized subgrid-diffusivity matrix              vg 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcSubgrDiffMatrix(
    const DRT::Element* ele, Epetra_SerialDenseMatrix& emat)
{
  /*----------------------------------------------------------------------*/
  // integration loop for one element
  /*----------------------------------------------------------------------*/
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    for (int k = 0; k < numscal_; ++k)
    {
      // set diffusion coeff to 1.0
      diffmanager_->SetIsotropicDiff(1.0, k);

      // calculation of diffusive element matrix
      double timefacfac = scatraparatimint_->TimeFac() * fac;
      CalcMatDiff(emat, k, timefacfac);

      /*subtract SUPG term */
      // emat(fvi,fui) -= taufac*conv(vi)*conv(ui);
    }
  }  // integration loop

  return;
}  // ScaTraImpl::CalcSubgrDiffMatrix


/*----------------------------------------------------------------------------------------*
 | finite difference check on element level (for debugging only) (protected)   fang 10/14 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::FDCheck(DRT::Element* ele,
    Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs,
    Epetra_SerialDenseVector& subgrdiff)
{
  // screen output
  std::cout << "FINITE DIFFERENCE CHECK FOR ELEMENT " << ele->Id();

  // make a copy of state variables to undo perturbations later
  std::vector<LINALG::Matrix<nen_, 1>> ephinp_original(numscal_);
  for (int k = 0; k < numscal_; ++k)
    for (unsigned i = 0; i < nen_; ++i) ephinp_original[k](i, 0) = ephinp_[k](i, 0);

  // generalized-alpha time integration requires a copy of history variables as well
  std::vector<LINALG::Matrix<nen_, 1>> ehist_original(numscal_);
  if (scatraparatimint_->IsGenAlpha())
  {
    for (int k = 0; k < numscal_; ++k)
      for (unsigned i = 0; i < nen_; ++i) ehist_original[k](i, 0) = ehist_[k](i, 0);
  }

  // initialize element matrix and vectors for perturbed state
  Epetra_SerialDenseMatrix emat_dummy(emat);
  Epetra_SerialDenseVector erhs_perturbed(erhs);
  Epetra_SerialDenseVector subgrdiff_dummy(subgrdiff);

  // initialize counter for failed finite difference checks
  unsigned counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over columns of element matrix by first looping over nodes and then over dofs at each node
  for (unsigned inode = 0; inode < nen_; ++inode)
  {
    for (int idof = 0; idof < numdofpernode_; ++idof)
    {
      // number of current column of element matrix
      unsigned col = inode * numdofpernode_ + idof;

      // clear element matrix and vectors for perturbed state
      emat_dummy.Scale(0.0);
      erhs_perturbed.Scale(0.0);
      subgrdiff_dummy.Scale(0.0);

      // fill state vectors with original state variables
      for (int k = 0; k < numscal_; ++k)
        for (unsigned i = 0; i < nen_; ++i) ephinp_[k](i, 0) = ephinp_original[k](i, 0);
      if (scatraparatimint_->IsGenAlpha())
        for (int k = 0; k < numscal_; ++k)
          for (unsigned i = 0; i < nen_; ++i) ehist_[k](i, 0) = ehist_original[k](i, 0);

      // impose perturbation
      if (scatraparatimint_->IsGenAlpha())
      {
        // perturbation of phi(n+alphaF), not of phi(n+1) => scale epsilon by factor alphaF
        ephinp_[idof](inode, 0) += scatraparatimint_->AlphaF() * scatrapara_->FDCheckEps();

        // perturbation of phi(n+alphaF) by alphaF*epsilon corresponds to perturbation of phidtam
        // (stored in ehist_) by alphaM*epsilon/(gamma*dt); note: alphaF/timefac = alphaM/(gamma*dt)
        ehist_[idof](inode, 0) +=
            scatraparatimint_->AlphaF() / scatraparatimint_->TimeFac() * scatrapara_->FDCheckEps();
      }
      else
        ephinp_[idof](inode, 0) += scatrapara_->FDCheckEps();

      // calculate element right-hand side vector for perturbed state
      Sysmat(ele, emat_dummy, erhs_perturbed, subgrdiff_dummy);

      // Now we compare the difference between the current entries in the element matrix
      // and their finite difference approximations according to
      // entries ?= (-erhs_perturbed + erhs_original) / epsilon

      // Note that the element right-hand side equals the negative element residual.
      // To account for errors due to numerical cancellation, we additionally consider
      // entries - erhs_original / epsilon ?= -erhs_perturbed / epsilon

      // Note that we still need to evaluate the first comparison as well. For small entries in the
      // element matrix, the second comparison might yield good agreement in spite of the entries
      // being wrong!
      for (unsigned row = 0; row < static_cast<unsigned>(numdofpernode_ * nen_); ++row)
      {
        // get current entry in original element matrix
        const double entry = emat(row, col);

        // finite difference suggestion (first divide by epsilon and then subtract for better
        // conditioning)
        const double fdval = -erhs_perturbed(row) / scatrapara_->FDCheckEps() +
                             erhs(row) / scatrapara_->FDCheckEps();

        // confirm accuracy of first comparison
        if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in first comparison
        const double abserr1 = entry - fdval;
        if (abs(abserr1) > abs(maxabserr)) maxabserr = abserr1;
        double relerr1(0.);
        if (abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if (abs(fdval) > 1.e-17)
          relerr1 = abserr1 / abs(fdval);
        if (abs(relerr1) > abs(maxrelerr)) maxrelerr = relerr1;

        // evaluate first comparison
        if (abs(relerr1) > scatrapara_->FDCheckTol())
        {
          if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
          std::cout << "emat[" << row << "," << col << "]:  " << entry << "   ";
          std::cout << "finite difference suggestion:  " << fdval << "   ";
          std::cout << "absolute error:  " << abserr1 << "   ";
          std::cout << "relative error:  " << relerr1 << std::endl;

          counter++;
        }

        // first comparison OK
        else
        {
          // left-hand side in second comparison
          const double left = entry - erhs(row) / scatrapara_->FDCheckEps();

          // right-hand side in second comparison
          const double right = -erhs_perturbed(row) / scatrapara_->FDCheckEps();

          // confirm accuracy of second comparison
          if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
            dserror("Finite difference check involves values too close to numerical zero!");

          // absolute and relative errors in second comparison
          const double abserr2 = left - right;
          if (abs(abserr2) > abs(maxabserr)) maxabserr = abserr2;
          double relerr2(0.);
          if (abs(left) > 1.e-17)
            relerr2 = abserr2 / abs(left);
          else if (abs(right) > 1.e-17)
            relerr2 = abserr2 / abs(right);
          if (abs(relerr2) > abs(maxrelerr)) maxrelerr = relerr2;

          // evaluate second comparison
          if (abs(relerr2) > scatrapara_->FDCheckTol())
          {
            if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
            std::cout << "emat[" << row << "," << col << "]-erhs[" << row << "]/eps:  " << left
                      << "   ";
            std::cout << "-erhs_perturbed[" << row << "]/eps:  " << right << "   ";
            std::cout << "absolute error:  " << abserr2 << "   ";
            std::cout << "relative error:  " << relerr2 << std::endl;

            counter++;
          }
        }
      }
    }
  }

  // screen output in case finite difference check is passed
  if (!counter)
    std::cout << " --> PASSED WITH MAXIMUM ABSOLUTE ERROR " << maxabserr
              << " AND MAXIMUM RELATIVE ERROR " << maxrelerr << std::endl;

  // undo perturbations of state variables
  for (int k = 0; k < numscal_; ++k)
    for (unsigned i = 0; i < nen_; ++i) ephinp_[k](i, 0) = ephinp_original[k](i, 0);
  if (scatraparatimint_->IsGenAlpha())
    for (int k = 0; k < numscal_; ++k)
      for (unsigned i = 0; i < nen_; ++i) ehist_[k](i, 0) = ehist_original[k](i, 0);

  return;
}

/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& errors)
{
  if (DRT::INPUT::get<SCATRA::Action>(params, "action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // set constants for analytical solution
  const double t = scatraparatimint_->Time();

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype =
      DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
  switch (errortype)
  {
    case INPAR::SCATRA::calcerror_byfunction:
    {
      const int errorfunctno = params.get<int>("error function number");

      // analytical solution
      double phi_exact(0.0);
      double deltaphi(0.0);
      //! spatial gradient of current scalar value
      LINALG::Matrix<nsd_, 1> gradphi(true);
      LINALG::Matrix<nsd_, 1> gradphi_exact(true);
      LINALG::Matrix<nsd_, 1> deltagradphi(true);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get coordinates at integration point
        // gp reference coordinates
        LINALG::Matrix<nsd_, 1> xyzint(true);
        xyzint.Multiply(xyze_, funct_);

        // function evaluation requires a 3D position vector!!
        double position[3] = {0.0, 0.0, 0.0};

        for (unsigned dim = 0; dim < nsd_; ++dim) position[dim] = xyzint(dim);

        for (int k = 0; k < numdofpernode_; ++k)
        {
          // scalar at integration point at time step n+1
          const double phinp = funct_.Dot(ephinp_[k]);
          // spatial gradient of current scalar value
          gradphi.Multiply(derxy_, ephinp_[k]);

          phi_exact = DRT::Problem::Instance()->Funct(errorfunctno - 1).Evaluate(k, position, t);

          std::vector<double> gradphi_exact_vec = DRT::Problem::Instance()
                                                      ->Funct(errorfunctno - 1)
                                                      .EvaluateSpatialDerivative(k, position, t);

          if (gradphi_exact_vec.size())
          {
            if (nsd_ == nsd_ele_)
              for (unsigned dim = 0; dim < nsd_; ++dim) gradphi_exact(dim) = gradphi_exact_vec[dim];
            else
            {
              // std::cout<<"Warning: Gradient of analytical solution cannot be evaluated correctly
              // for transport on curved surfaces!"<<std::endl;
              gradphi_exact.Clear();
            }
          }
          else
          {
            std::cout << "Warning: Gradient of analytical solution was not evaluated!" << std::endl;
            gradphi_exact.Clear();
          }

          // error at gauss point
          deltaphi = phinp - phi_exact;
          deltagradphi.Update(1.0, gradphi, -1.0, gradphi_exact);

          // 0: delta scalar for L2-error norm
          // 1: delta scalar for H1-error norm
          // 2: analytical scalar for L2 norm
          // 3: analytical scalar for H1 norm

          // the error for the L2 and H1 norms are evaluated at the Gauss point

          // integrate delta scalar for L2-error norm
          errors(k * 4 + 0) += deltaphi * deltaphi * fac;
          // integrate delta scalar for H1-error norm
          errors(k * 4 + 1) += deltaphi * deltaphi * fac;
          // integrate analytical scalar for L2 norm
          errors(k * 4 + 2) += phi_exact * phi_exact * fac;
          // integrate analytical scalar for H1 norm
          errors(k * 4 + 3) += phi_exact * phi_exact * fac;

          // integrate delta scalar derivative for H1-error norm
          errors(k * 4 + 1) += deltagradphi.Dot(deltagradphi) * fac;
          // integrate analytical scalar derivative for H1 norm
          errors(k * 4 + 3) += gradphi_exact.Dot(gradphi_exact) * fac;
        }
      }  // loop over integration points
    }
    break;

    case INPAR::SCATRA::calcerror_spherediffusion:
    {
      // analytical solution
      double phi_exact(0.0);
      double deltaphi(0.0);
      //! spatial gradient of current scalar value
      LINALG::Matrix<nsd_, 1> gradphi(true);
      LINALG::Matrix<nsd_, 1> gradphi_exact(true);
      LINALG::Matrix<nsd_, 1> deltagradphi(true);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get coordinates at integration point
        // gp reference coordinates
        LINALG::Matrix<nsd_, 1> xyzint(true);
        xyzint.Multiply(xyze_, funct_);

        for (int k = 0; k < numscal_; k++)
        {
          const double x = xyzint(0);
          const double y = xyzint(1);
          const double z = xyzint(2);

          // scalar at integration point at time step n+1
          const double phinp = funct_.Dot(ephinp_[k]);
          // spatial gradient of current scalar value
          gradphi.Multiply(derxy_, ephinp_[k]);

          phi_exact = exp(-6 * t) * x * y + 10;

          gradphi_exact(0) = (1.0 - 2.0 * x * x) * y * exp(-6.0 * t);
          gradphi_exact(1) = (1.0 - 2.0 * y * y) * x * exp(-6.0 * t);
          gradphi_exact(2) = -2.0 * x * y * z * exp(-6.0 * t);

          // error at gauss point
          deltaphi = phinp - phi_exact;
          deltagradphi.Update(1.0, gradphi, -1.0, gradphi_exact);

          // 0: delta scalar for L2-error norm
          // 1: delta scalar for H1-error norm
          // 2: analytical scalar for L2 norm
          // 3: analytical scalar for H1 norm

          // the error for the L2 and H1 norms are evaluated at the Gauss point

          // integrate delta scalar for L2-error norm
          errors(k * numscal_ + 0) += deltaphi * deltaphi * fac;
          // integrate delta scalar for H1-error norm
          errors(k * numscal_ + 1) += deltaphi * deltaphi * fac;
          // integrate analytical scalar for L2 norm
          errors(k * numscal_ + 2) += phi_exact * phi_exact * fac;
          // integrate analytical scalar for H1 norm
          errors(k * numscal_ + 3) += phi_exact * phi_exact * fac;

          // integrate delta scalar derivative for H1-error norm
          errors(k * numscal_ + 1) += deltagradphi.Dot(deltagradphi) * fac;
          // integrate analytical scalar derivative for H1 norm
          errors(k * numscal_ + 3) += gradphi_exact.Dot(gradphi_exact) * fac;
        }
      }  // loop over integration points
    }
    break;
    default:
      dserror("Unknown analytical solution!");
      break;
  }  // switch(errortype)

  return;
}  // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalErrorComparedToAnalytSolution


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                  vuong 07/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcHeteroReacMatAndRHS(
    DRT::Element* ele,               ///< the element whose matrix is calculated
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs   ///< element rhs to calculate
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  const double vol = EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_, 0.0);

  if (not scatrapara_->TauGP())
  {
    for (int k = 0; k < numscal_; ++k)  // loop of each transported scalar
    {
      // get velocity at element center
      LINALG::Matrix<nsd_, 1> convelint = scatravarmanager_->ConVel(k);
      // calculation of stabilization parameter at element center
      CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
          reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k], convelint,
          vol);
    }
  }

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP())
  {
    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // loop all scalars
    for (int k = 0; k < numscal_; ++k)  // deal with a system of transported scalars
    {
      // reactive part of the form: (reaction coefficient)*phi
      double rea_phi(0.0);
      rea_phi = densnp[k] * scatravarmanager_->Phinp(k) * reamanager_->GetReaCoeff(k);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      GetRhsInt(rhsint, densnp[k], k);

      double scatrares(0.0);
      // calculate strong residual
      CalcStrongResidual(k, scatrares, densam[k], densnp[k], rea_phi, rhsint, tau[k]);

      if (scatrapara_->TauGP())
      {
        // (re)compute stabilization parameter at integration point, since diffusion may have
        // changed
        CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
            reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k],
            scatravarmanager_->ConVel(k), vol);  // TODO:(Thon) do we really have to do this??
      }

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      // stabilization parameter and integration factors
      const double taufac = tau[k] * fac;
      const double timefacfac = scatraparatimint_->TimeFac() * fac;
      const double timetaufac = scatraparatimint_->TimeFac() * taufac;

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      LINALG::Matrix<nen_, 1> sgconv(true);
      LINALG::Matrix<nen_, 1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
      }

      // including stabilization
      if (reamanager_->Active())
      {
        CalcMatReact(emat, k, timefacfac, timetaufac, taufac, densnp[k], sgconv, diff);
      }

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac = scatraparatimint_->TimeFacRhs() * fac;
      double rhstaufac = scatraparatimint_->TimeFacRhsTau() * taufac;

      ComputeRhsInt(rhsint, densam[k], densnp[k], 0.0);

      RecomputeScatraResForRhs(scatrares, k, diff, densn[k], densnp[k], rea_phi, rhsint);

      //----------------------------------------------------------------
      // standard Galerkin transient, old part of rhs and bodyforce term
      //----------------------------------------------------------------
      CalcRHSHistAndSource(erhs, k, fac, rhsint);

      //----------------------------------------------------------------
      // reactive terms (standard Galerkin and stabilization) on rhs
      //----------------------------------------------------------------

      if (reamanager_->Active())
        CalcRHSReact(erhs, k, rhsfac, rhstaufac, rea_phi, densnp[k], scatrares);

    }  // end loop all scalars
  }    // end loop Gauss points

  return;
}


// template classes

#include "scatra_ele_calc_fwd.hpp"
