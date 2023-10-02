/*----------------------------------------------------------------------*/
/*! \file
\brief quadratic nonlinear tetrahedron
\level 1
*----------------------------------------------------------------------*/
#include "baci_contact_analytical.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_gauss_point_extrapolation.H"
#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_fiber_node.H"
#include "baci_lib_element.H"
#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_prestress_service.H"
#include "baci_lib_utils.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_element_service.H"
#include "baci_so3_prestress.H"
#include "baci_so3_tet10.H"
#include "baci_so3_utils.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_SerialDenseSolver.hpp>

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  EnsureMaterialPostSetup(params);

  SetParamsInterfacePtr(params);

  CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET10, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET10, 1> elevec2(elevec2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET10, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  DRT::ELEMENTS::So_tet10::ActionType act = So_tet10::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_tet10::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_tet10::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_tet10::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_tet10::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_tet10::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_tet10::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = So_tet10::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = So_tet10::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_tet10::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = So_tet10::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = So_tet10::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = So_tet10::calc_struct_reset_all;
  else if (action == "calc_struct_energy")
    act = So_tet10::calc_struct_energy;
  else if (action == "calc_struct_errornorms")
    act = So_tet10::calc_struct_errornorms;
  else if (action == "calc_struct_prestress_update")
    act = So_tet10::prestress_update;
  else if (action == "calc_global_gpstresses_map")
    act = So_tet10::calc_global_gpstresses_map;
  else if (action == "struct_init_gauss_point_data_output")
    act = So_tet10::struct_init_gauss_point_data_output;
  else if (action == "struct_gauss_point_data_output")
    act = So_tet10::struct_gauss_point_data_output;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for So_tet10");

  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (double& i : mydisp) i = 0.0;
      std::vector<double> myres(lm.size());
      for (double& myre : myres) myre = 0.0;
      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &elemat1, nullptr,
          &elevec1, nullptr, &elevec3, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, matptr, nullptr,
          &elevec1, nullptr, &elevec3, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      std::vector<double> mydispmat(lm.size(), 0.0);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> myemat(true);

      so_tet10_nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &myemat, nullptr,
          &elevec1, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      Teuchos::RCP<const Epetra_Vector> acc = discretization.GetState("acceleration");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      if (acc == Teuchos::null) dserror("Cannot get state vectors 'acceleration'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
      std::vector<double> myacc(lm.size());
      DRT::UTILS::ExtractMyValues(*acc, myacc, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &elemat2,
          &elevec1, &elevec2, &elevec3, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);

      if (act == calc_struct_nlnstifflmass) so_tet10_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);
        CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> stress;
        CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> strain;
        auto iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        auto iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        so_tet10_nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr,
            nullptr, nullptr, nullptr, &stress, &strain, params, iostress, iostrain);
        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }

        // output of rotation matrix R with F = U*R
        if (DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_ROT") ==
            true)
        {
          CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> R;
          DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(this, mydisp, R);

          Teuchos::RCP<std::vector<char>> rotdata =
              params.get<Teuchos::RCP<std::vector<char>>>("rotation", Teuchos::null);

          DRT::PackBuffer data;
          AddtoPack(data, R);
          data.StartPacking();
          AddtoPack(data, R);
          std::copy(data().begin(), data().end(), std::back_inserter(*rotdata));
        }
      }
    }
    break;
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    case calc_struct_update_istep:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      Update_element(mydisp, params, Material());
    }
    break;

    case prestress_update:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      switch (pstype_)
      {
        case INPAR::STR::PreStress::mulf:
        {
          // build incremental def gradient for every gauss point
          CORE::LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET10, 9);
          DefGradient(mydisp, gpdefgrd, *prestress_);

          // update deformation gradient and put back to storage
          CORE::LINALG::Matrix<3, 3> deltaF;
          CORE::LINALG::Matrix<3, 3> Fhist;
          CORE::LINALG::Matrix<3, 3> Fnew;
          for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
          {
            prestress_->StoragetoMatrix(gp, deltaF, gpdefgrd);
            prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());
            Fnew.Multiply(deltaF, Fhist);
            prestress_->MatrixtoStorage(gp, Fnew, prestress_->FHistory());
          }

          // push-forward invJ for every gaussian point
          UpdateJacobianMapping(mydisp, *prestress_);

          // Update constraintmixture material
          if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture)
          {
            SolidMaterial()->Update();
          }
          break;
        }
        default:
          dserror(
              "You should either not be here, or the prestressing method you are using is not "
              "implemented for TET10 elements!");
      }
    }
    break;

    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;
    //==================================================================================
    case calc_struct_energy:
    {
      // initialization of internal energy
      double intenergy = 0.0;

      const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_4gp =
          so_tet10_4gp_shapefcts();
      const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_4gp =
          so_tet10_4gp_derivs();
      const static std::vector<double> gpweights_4gp = so_tet10_4gp_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // update element geometry
      CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xcurr;  // current  coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
      UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);
      UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(mydisp, xdisp);
      UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet10, 3>(xrefe, xdisp, xcurr);

      /* =========================================================================*/
      /* ================================================= Loop over Gauss Points */
      /* =========================================================================*/
      CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ;
      // build deformation gradient wrt to material configuration
      // in case of prestressing, build defgrd wrt to last stored configuration
      CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
      for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
      {
        /* get the inverse of the Jacobian matrix which looks like:
        **            [ x_,r  y_,r  z_,r ]^-1
        **     J^-1 = [ x_,s  y_,s  z_,s ]
        **            [ x_,t  y_,t  z_,t ]
        */
        // compute derivatives N_XYZ at gp w.r.t. material coordinates
        // by N_XYZ = J^-1 * N_rst
        N_XYZ.Multiply(invJ_[gp], derivs_4gp[gp]);

        // Gauss weights and Jacobian determinant
        double fac = detJ_[gp] * gpweights_4gp[gp];

        if (::UTILS::PRESTRESS::IsMulf(pstype_))
        {
          // get Jacobian mapping wrt to the stored configuration
          CORE::LINALG::Matrix<3, 3> invJdef;
          prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
          // get derivatives wrt to last spatial configuration
          CORE::LINALG::Matrix<3, 10> N_xyz;
          N_xyz.Multiply(invJdef, derivs_4gp[gp]);

          // build multiplicative incremental defgrd
          defgrd.MultiplyTT(xdisp, N_xyz);
          defgrd(0, 0) += 1.0;
          defgrd(1, 1) += 1.0;
          defgrd(2, 2) += 1.0;

          // get stored old incremental F
          CORE::LINALG::Matrix<3, 3> Fhist;
          prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

          // build total defgrd = delta F * F_old
          CORE::LINALG::Matrix<3, 3> Fnew;
          Fnew.Multiply(defgrd, Fhist);
          defgrd = Fnew;
        }
        else
        {
          defgrd.MultiplyTT(xcurr, N_XYZ);
        }

        // Right Cauchy-Green tensor = F^T * F
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.values(), true);
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);

        // call material for evaluation of strain energy function
        double psi = 0.0;
        SolidMaterial()->StrainEnergy(glstrain, psi, gp, Id());

        // sum up GP contribution to internal energy
        intenergy += fac * psi;
      }

      if (IsParamsInterface())  // new structural time integration
      {
        StrParamsInterface().AddContributionToEnergyType(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        // check length of elevec1
        if (elevec1_epetra.length() < 1) dserror("The given result vector is too short.");

        elevec1_epetra(0) = intenergy;
      }
    }
    break;


    case calc_struct_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - only implemented for SVK material (relevant for energy norm only, L2 and
      //   H1 norms are of course valid for arbitrary materials)
      // - analytical solutions are currently stored in a repository in the CONTACT
      //   namespace, however they could (should?) be moved to a more general location

      // check length of elevec1
      if (elevec1_epetra.length() < 3) dserror("The given result vector is too short.");

      // check material law
      Teuchos::RCP<MAT::Material> mat = Material();

      //******************************************************************
      // only for St.Venant Kirchhoff material
      //******************************************************************
      if (mat->MaterialType() == INPAR::MAT::m_stvenant)
      {
        // declaration of variables
        double l2norm = 0.0;
        double h1norm = 0.0;
        double energynorm = 0.0;

        // shape functions, derivatives and integration weights
        const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> vals =
            so_tet10_11gp_shapefcts();
        const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
            so_tet10_11gp_derivs();
        const static std::vector<double> weights = so_tet10_11gp_weights();

        // get displacements and extract values of this element
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

        // nodal displacement vector
        CORE::LINALG::Matrix<NUMDOF_SOTET10, 1> nodaldisp;
        for (int i = 0; i < NUMDOF_SOTET10; ++i) nodaldisp(i, 0) = mydisp[i];

        // reference geometry (nodal positions)
        CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
        UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);

        // deformation gradient = identity tensor (geometrically linear case!)
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
        for (int i = 0; i < NUMDIM_SOTET10; ++i) defgrd(i, i) = 1;

        // use of 11 GP for the errornorm computation

        //----------------------------------------------------------------
        // loop over all Gauss points
        //----------------------------------------------------------------
        for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
        {
          // Gauss weights and Jacobian determinant
          double fac = detJ_mass_[gp] * weights[gp];

          // Gauss point in reference configuration
          CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> xgp(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int n = 0; n < NUMNOD_SOTET10; ++n) xgp(k, 0) += (vals[gp])(n)*xrefe(n, k);

          //**************************************************************
          // get analytical solution
          CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> uanalyt(true);
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainanalyt(true);
          CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> derivanalyt(true);

          CONTACT::AnalyticalSolutions3D(xgp, uanalyt, strainanalyt, derivanalyt);
          //**************************************************************

          //--------------------------------------------------------------
          // (1) L2 norm
          //--------------------------------------------------------------

          // compute displacements at GP
          CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> ugp(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int n = 0; n < NUMNOD_SOTET10; ++n)
              ugp(k, 0) += (vals[gp])(n)*nodaldisp(NODDOF_SOTET10 * n + k, 0);

          // displacement error
          CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> uerror(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k) uerror(k, 0) = uanalyt(k, 0) - ugp(k, 0);

          // compute GP contribution to L2 error norm
          l2norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (2) H1 norm
          //--------------------------------------------------------------

          // compute derivatives N_XYZ at GP w.r.t. material coordinates
          // by N_XYZ = J^-1 * N_rst
          CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ(true);
          N_XYZ.Multiply(invJ_mass_[gp], derivs[gp]);

          // compute partial derivatives at GP
          CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> derivgp(true);
          for (int l = 0; l < NUMDIM_SOTET10; ++l)
            for (int m = 0; m < NUMDIM_SOTET10; ++m)
              for (int k = 0; k < NUMNOD_SOTET10; ++k)
                derivgp(l, m) += N_XYZ(m, k) * nodaldisp(NODDOF_SOTET10 * k + l, 0);

          // derivative error
          CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> deriverror(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int m = 0; m < NUMDIM_SOTET10; ++m)
              deriverror(k, m) = derivanalyt(k, m) - derivgp(k, m);

          // compute GP contribution to H1 error norm
          h1norm += fac * deriverror.Dot(deriverror);
          h1norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (3) Energy norm
          //--------------------------------------------------------------

          // compute linear B-operator
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
          for (int i = 0; i < NUMNOD_SOTET10; ++i)
          {
            bop(0, NODDOF_SOTET10 * i + 0) = N_XYZ(0, i);
            bop(0, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(0, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(1, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(1, NODDOF_SOTET10 * i + 1) = N_XYZ(1, i);
            bop(1, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 2) = N_XYZ(2, i);

            bop(3, NODDOF_SOTET10 * i + 0) = N_XYZ(1, i);
            bop(3, NODDOF_SOTET10 * i + 1) = N_XYZ(0, i);
            bop(3, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(4, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(4, NODDOF_SOTET10 * i + 1) = N_XYZ(2, i);
            bop(4, NODDOF_SOTET10 * i + 2) = N_XYZ(1, i);
            bop(5, NODDOF_SOTET10 * i + 0) = N_XYZ(2, i);
            bop(5, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(5, NODDOF_SOTET10 * i + 2) = N_XYZ(0, i);
          }

          // compute linear strain at GP
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> straingp(true);
          straingp.Multiply(bop, nodaldisp);

          // strain error
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainerror(true);
          for (int k = 0; k < MAT::NUM_STRESS_3D; ++k)
            strainerror(k, 0) = strainanalyt(k, 0) - straingp(k, 0);

          // compute stress vector and constitutive matrix
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

          SolidMaterial()->Evaluate(&defgrd, &strainerror, params, &stress, &cmat, gp, Id());

          // compute GP contribution to energy error norm
          energynorm += fac * stress.Dot(strainerror);

          // std::cout << "UAnalytical:      " << ugp << std::endl;
          // std::cout << "UDiscrete:        " << uanalyt << std::endl;
          // std::cout << "StrainAnalytical: " << strainanalyt << std::endl;
          // std::cout << "StrainDiscrete:   " << straingp << std::endl;
          // std::cout << "DerivAnalytical:  " << derivanalyt << std::endl;
          // std::cout << "DerivDiscrete:    " << derivgp << std::endl;
        }
        //----------------------------------------------------------------

        // return results
        elevec1_epetra(0) = l2norm;
        elevec1_epetra(1) = h1norm;
        elevec1_epetra(2) = energynorm;
      }
      else
        dserror("ERROR: Error norms only implemented for SVK material");

      ;  // there is nothing to do here at the moment
    }
    break;

      // evaluate stresses and strains at gauss points and store gpstresses in map <EleId,
      // gpstresses >
    case calc_global_gpstresses_map:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstressmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstressmap", Teuchos::null);
        if (gpstressmap == Teuchos::null)
          dserror("no gp stress map available for writing gpstresses");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstrainmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstrainmap", Teuchos::null);
        if (gpstrainmap == Teuchos::null)
          dserror("no gp strain map available for writing gpstrains");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);
        CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> stress;
        CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> strain;
        auto iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        auto iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        // if a linear analysis is desired
        if (kintype_ == INPAR::STR::kinem_linear)
        {
          dserror("Linear case not implemented");
        }


        else
        {
          so_tet10_nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr,
              nullptr, nullptr, nullptr, &stress, &strain, params, iostress, iostrain);
        }
        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstress =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstress->shape(NUMGPT_SOTET10, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET10; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstrain->shape(NUMGPT_SOTET10, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET10; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstrain)(i, j) = strain(i, j);
          }
        }

        // add to map
        (*gpstressmap)[gid] = gpstress;
        (*gpstrainmap)[gid] = gpstrain;

        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }

        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;
    case struct_init_gauss_point_data_output:
    {
      dsassert(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Save number of Gauss of the element for gauss point data output
      StrParamsInterface().GaussPointDataOutputManagerPtr()->AddElementNumberOfGaussPoints(
          NUMGPT_SOTET10);

      // holder for output quantity names and their size
      std::unordered_map<std::string, int> quantities_map{};

      // Ask material for the output quantity names and sizes
      SolidMaterial()->RegisterVtkOutputDataNames(quantities_map);

      // Add quantities to the Gauss point output data manager (if they do not already exist)
      StrParamsInterface().GaussPointDataOutputManagerPtr()->MergeQuantities(quantities_map);
    }
    break;
    case struct_gauss_point_data_output:
    {
      dsassert(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Collection and assembly of gauss point data
      for (const auto& quantity :
          StrParamsInterface().GaussPointDataOutputManagerPtr()->GetQuantities())
      {
        const std::string& quantity_name = quantity.first;
        const int quantity_size = quantity.second;

        // Step 1: Collect the data for each Gauss point for the material
        CORE::LINALG::SerialDenseMatrix gp_data(NUMGPT_SOTET10, quantity_size, true);
        bool data_available = SolidMaterial()->EvaluateVtkOutputData(quantity_name, gp_data);

        // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
        // point)
        if (data_available)
        {
          switch (StrParamsInterface().GaussPointDataOutputManagerPtr()->GetOutputType())
          {
            case INPAR::STR::GaussPointDataOutputType::element_center:
            {
              // compute average of the quantities
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  StrParamsInterface().GaussPointDataOutputManagerPtr()->GetElementCenterData().at(
                      quantity_name);
              CORE::DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::nodes:
            {
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  StrParamsInterface().GaussPointDataOutputManagerPtr()->GetNodalData().at(
                      quantity_name);

              Epetra_IntVector& global_nodal_element_count =
                  *StrParamsInterface().GaussPointDataOutputManagerPtr()->GetNodalDataCount().at(
                      quantity_name);

              static auto gauss_integration =
                  CORE::DRT::UTILS::IntegrationPoints3D(CORE::DRT::UTILS::NumGaussPointsToGaussRule<
                      DRT::Element::DiscretizationType::tet10>(NUMGPT_SOTET10));
              CORE::DRT::UTILS::ExtrapolateGPQuantityToNodesAndAssemble<
                  DRT::Element::DiscretizationType::tet10>(
                  *this, gp_data, *global_data, false, gauss_integration);
              DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::gauss_points:
            {
              std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                  StrParamsInterface().GaussPointDataOutputManagerPtr()->GetGaussPointData().at(
                      quantity_name);
              DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::none:
              dserror(
                  "You specified a Gauss point data output type of none, so you should not end up "
                  "here.");
            default:
              dserror("Unknown Gauss point data output type.");
          }
        }
      }
    }
    break;
    default:
      dserror("Unknown type of action for SolidTet10");
      break;
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = condition.Get<std::vector<int>>("onoff");
  const auto* val = condition.Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOTET10)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOTET10; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = condition.Get<std::vector<int>>("funct");
  CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts =
      so_tet10_4gp_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();
  const static std::vector<double> gpweights = so_tet10_4gp_weights();
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);

  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // compute the Jacobian matrix
    CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> jac;
    jac.Multiply(derivs[gp], xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOTET10; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    double fac = gpweights[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0) ? DRT::Problem::Instance()
                                 ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .Evaluate(xrefegp.A(), time, dim)
                           : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SOTET10; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOTET10 + dim] += shapefcts[gp](nodid) * dim_fac;
        }
      }
    }

  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_tet10::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::InitJacobianMapping()
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs4gp =
      so_tet10_4gp_derivs();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs11gp =
      so_tet10_11gp_derivs();

  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);

  // Initialize for stiffness integration with 4 GPs
  invJ_.resize(NUMGPT_SOTET10);
  detJ_.resize(NUMGPT_SOTET10);
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    invJ_[gp].Multiply(derivs4gp[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    if (::UTILS::PRESTRESS::IsMulfActive(time_, pstype_, pstime_))
      if (!(prestress_->IsInit()))
        prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());
  }

  if (::UTILS::PRESTRESS::IsMulfActive(time_, pstype_, pstime_)) prestress_->IsInit() = true;

  // Initialize for mass integration with 10 GPs

  invJ_mass_.resize(NUMGPT_MASS_SOTET10);
  detJ_mass_.resize(NUMGPT_MASS_SOTET10);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; ++gp)
  {
    invJ_mass_[gp].Multiply(derivs11gp[gp], xrefe);
    detJ_mass_[gp] = invJ_mass_[gp].Invert();
    if (detJ_mass_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_mass_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                             // current displacements
    std::vector<double>* vel,                                              // current velocities
    std::vector<double>* acc,                                              // current accelerations
    std::vector<double>& residual,                                         // current residual displ
    std::vector<double>& dispmat,  // current material displacements
    CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* massmatrix,   // element mass matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET10, 1>* force,       // element internal force vector
    CORE::LINALG::Matrix<NUMDOF_SOTET10, 1>* forceinert,  // element inertial force vector
    CORE::LINALG::Matrix<NUMDOF_SOTET10, 1>* force_str,   // element structural force vector
    CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    CORE::LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  // stress output option
    const INPAR::STR::StrainType iostrain)  // strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_4gp =
      so_tet10_4gp_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_4gp =
      so_tet10_4gp_derivs();
  const static std::vector<double> gpweights_4gp = so_tet10_4gp_weights();
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xcurr;  // current  coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
  UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);
  UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(disp, xdisp);
  UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet10, 3>(xrefe, xdisp, xcurr);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp], derivs_4gp[gp]);
    double detJ = detJ_[gp];

    if (::UTILS::PRESTRESS::IsMulf(pstype_))
    {
      // get Jacobian mapping wrt to the stored configuration
      CORE::LINALG::Matrix<3, 3> invJdef;
      prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      CORE::LINALG::Matrix<3, 10> N_xyz;
      N_xyz.Multiply(invJdef, derivs_4gp[gp]);

      // build multiplicative incremental defgrd
      defgrd.MultiplyTT(xdisp, N_xyz);
      defgrd(0, 0) += 1.0;
      defgrd(1, 1) += 1.0;
      defgrd(2, 2) += 1.0;

      // get stored old incremental F
      CORE::LINALG::Matrix<3, 3> Fhist;
      prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

      // build total defgrd = delta F * F_old
      CORE::LINALG::Matrix<3, 3> Fnew;
      Fnew.Multiply(defgrd, Fhist);
      defgrd = Fnew;
    }
    else
    {
      defgrd.MultiplyTT(xcurr, N_XYZ);
    }

    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.values(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == nullptr) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == nullptr) dserror("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        // inverse of deformation gradient
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> invdefgrd;
        invdefgrd.Invert(defgrd);

        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> temp;
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> euler_almansi;
        temp.Multiply(gl, invdefgrd);
        euler_almansi.MultiplyTN(invdefgrd, temp);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
        break;
    }

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
    for (int i = 0; i < NUMNOD_SOTET10; ++i)
    {
      bop(0, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture ||
        Material()->MaterialType() == INPAR::MAT::m_mixture)
    {
      // gp reference coordinates
      CORE::LINALG::Matrix<NUMNOD_SOTET10, 1> funct(true);
      funct = shapefcts_4gp[gp];
      CORE::LINALG::Matrix<1, NUMDIM_SOTET10> point(true);
      point.MultiplyTN(funct, xrefe);
      params.set("gprefecoord", point);
    }

    UTILS::GetTemperatureForStructuralMaterial<tet10>(shapefcts_4gp[gp], params);
    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) dserror("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) dserror("stress data not available");
        const double detF = defgrd.Determinant();

        CORE::LINALG::Matrix<3, 3> pkstress;
        pkstress(0, 0) = stress(0);
        pkstress(0, 1) = stress(3);
        pkstress(0, 2) = stress(5);
        pkstress(1, 0) = pkstress(0, 1);
        pkstress(1, 1) = stress(1);
        pkstress(1, 2) = stress(4);
        pkstress(2, 0) = pkstress(0, 2);
        pkstress(2, 1) = pkstress(1, 2);
        pkstress(2, 2) = stress(2);

        CORE::LINALG::Matrix<3, 3> temp;
        CORE::LINALG::Matrix<3, 3> cauchystress;
        temp.Multiply(1.0 / detF, defgrd, pkstress);
        cauchystress.MultiplyNT(temp, defgrd);

        (*elestress)(gp, 0) = cauchystress(0, 0);
        (*elestress)(gp, 1) = cauchystress(1, 1);
        (*elestress)(gp, 2) = cauchystress(2, 2);
        (*elestress)(gp, 3) = cauchystress(0, 1);
        (*elestress)(gp, 4) = cauchystress(1, 2);
        (*elestress)(gp, 5) = cauchystress(0, 2);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
        break;
    }

    double detJ_w = detJ * gpweights_4gp[gp];
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update of stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      CORE::LINALG::Matrix<6, NUMDOF_SOTET10> cb;
      cb.Multiply(cmat, bop);
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      CORE::LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.Scale(detJ_w);                       // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);             // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************
    }
  } /* ==================================================== end of Loop over GP */

  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 11 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_11gp =
      so_tet10_11gp_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_11gp =
      so_tet10_11gp_derivs();
  const static std::vector<double> gpweights_11gp = so_tet10_11gp_weights();
  /* ============================================================================*/

  if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
  {
    // since material has only 4 GPs with density information!
    double density = 0.;
    double vol = 0.;

    for (int gp4 = 0; gp4 < NUMGPT_SOTET10; ++gp4)
    {
      vol += gpweights_4gp[gp4];
      density += gpweights_4gp[gp4] * Material()->Density(gp4);
    }
    density /= vol;

    // consistent mass matrix evaluated using a 11-point rule
    for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
    {
      // integrate consistent mass matrix
      double detJ_mass = detJ_mass_[gp];
      double detJ_mass_w = detJ_mass * gpweights_11gp[gp];
      const double factor = detJ_mass_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
      {
        ifactor = shapefcts_11gp[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
        {
          massfactor = shapefcts_11gp[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOTET10 * inod + 0, NUMDIM_SOTET10 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10 * inod + 1, NUMDIM_SOTET10 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10 * inod + 2, NUMDIM_SOTET10 * jnod + 2) += massfactor;
        }
      }

      // check for non constant mass matrix
      if (SolidMaterial()->VaryingDensity())
      {
        /*
         If the density, i.e. the mass matrix, is not constant, a linearization is neccessary.
         In general, the mass matrix can be dependent on the displacements, the velocities and the
         accelerations. We write all the additional terms into the mass matrix, hence, conversion
         from accelerations to velocities and displacements are needed. As those conversions depend
         on the time integration scheme, the factors are set within the respective time integrators
         and read from the parameter list inside the element (this is a little ugly...). */
        CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> nxyz(nxyz_);  // copy!

        double timintfac_dis = 0.0;
        double timintfac_vel = 0.0;
        if (IsParamsInterface())
        {
          timintfac_dis = StrParamsInterface().GetTimIntFactorDisp();
          timintfac_vel = StrParamsInterface().GetTimIntFactorVel();
        }
        else
        {
          timintfac_dis = params.get<double>("timintfac_dis");
          timintfac_vel = params.get<double>("timintfac_vel");
        }
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_disp(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_vel(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass(true);

        /*
        ** B = F . Bl *
        **
        **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
        **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
        **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
        ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
        **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
        **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
        **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
        **      [                                                         ]
        **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
        **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
        **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
        **      [                                                         ]
        **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
        **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
        **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
        */
        // size is 6x12
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
        for (int i = 0; i < NUMNOD_SOTET10; i++)
        {
          bop(0, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * nxyz(i, 0);
          bop(1, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * nxyz(i, 1);
          bop(2, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * nxyz(i, 2);
          /* ~~~ */
          bop(3, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * nxyz(i, 1) + defgrd(0, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * nxyz(i, 1) + defgrd(1, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * nxyz(i, 1) + defgrd(2, 1) * nxyz(i, 0);
          bop(4, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * nxyz(i, 2) + defgrd(0, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * nxyz(i, 2) + defgrd(1, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * nxyz(i, 2) + defgrd(2, 2) * nxyz(i, 1);
          bop(5, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * nxyz(i, 0) + defgrd(0, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * nxyz(i, 0) + defgrd(1, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * nxyz(i, 0) + defgrd(2, 0) * nxyz(i, 2);
        }


        // Right Cauchy-Green tensor = F^T * F
        CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.values(), true);
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        SolidMaterial()->EvaluateNonLinMass(
            &defgrd, &glstrain, params, &linmass_disp, &linmass_vel, gp, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        CORE::LINALG::Matrix<NUMDIM_SOTET10, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
            myacc(idim) += shapefcts_11gp[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOTET10)];

        if (stiffmatrix != nullptr)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          CORE::LINALG::Matrix<1, NUMDOF_SOTET10> cb;
          cb.MultiplyTN(linmass_disp, bop);
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
          {
            double factor = detJ_mass_w * shapefcts_11gp[gp](inod);
            for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
            {
              double massfactor = factor * myacc(idim);
              for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
                for (int jdim = 0; jdim < NUMDIM_SOTET10; ++jdim)
                  (*massmatrix)(inod * NUMDIM_SOTET10 + idim, jnod * NUMDIM_SOTET10 + jdim) +=
                      massfactor * cb(jnod * NUMDIM_SOTET10 + jdim);
            }
          }
        }

        // internal force vector without EAS terms
        if (forceinert != nullptr)
        {
          // integrate nonlinear inertia force term
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
          {
            double forcefactor = shapefcts_11gp[gp](inod) * detJ_mass_w;
            for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
              (*forceinert)(inod * NUMDIM_SOTET10 + idim) += forcefactor * density * myacc(idim);
          }
        }
      }
    }

  }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  return;
}  // DRT::ELEMENTS::So_tet10::SOTET10_nlnstiffmass



/*----------------------------------------------------------------------*
 |  lump mass matrix                                                    |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_lumpmass(
    CORE::LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r = 0; r < (*emass).numRows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fcts at 4 Gauss Points                         |
 *----------------------------------------------------------------------*/
const std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>>
DRT::ELEMENTS::So_tet10::so_tet10_4gp_shapefcts()
{
  static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts(NUMGPT_SOTET10);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const CORE::DRT::UTILS::GaussRule3D gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    CORE::DRT::UTILS::shape_function_3D(shapefcts[gp], r, s, t, DRT::Element::tet10);
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 4 Gauss Points                   |
 *----------------------------------------------------------------------*/
const std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
DRT::ELEMENTS::So_tet10::so_tet10_4gp_derivs()
{
  static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs(NUMGPT_SOTET10);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  for (int gp = 0; gp < NUMGPT_SOTET10; gp++)
  {
    so_tet10_derivs<CORE::DRT::UTILS::GaussRule3D::tet_4point>(derivs[gp], gp);
  }
  derivs_done = true;

  return derivs;
}

template <CORE::DRT::UTILS::GaussRule3D intrule>
void DRT::ELEMENTS::So_tet10::so_tet10_derivs(
    CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>& derivs, const int gp) const
{
  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(intrule);

  const double r = intpoints.qxg[gp][0];
  const double s = intpoints.qxg[gp][1];
  const double t = intpoints.qxg[gp][2];

  CORE::DRT::UTILS::shape_function_3D_deriv1(derivs, r, s, t, DRT::Element::tet10);
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 4 Gauss Points                            |
 *----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_4gp_weights()
{
  static std::vector<double> weights(NUMGPT_SOTET10);
  static bool weights_done = false;
  if (weights_done) return weights;

  const CORE::DRT::UTILS::GaussRule3D gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_SOTET10; gp++) weights[gp] = intpoints.qwgt[gp];
  weights_done = true;

  return weights;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fcts at 10 Gauss Points                        |
 *----------------------------------------------------------------------*/
const std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>>&
DRT::ELEMENTS::So_tet10::so_tet10_11gp_shapefcts()
{
  static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts(NUMGPT_MASS_SOTET10);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const CORE::DRT::UTILS::GaussRule3D gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_11point;
  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    CORE::DRT::UTILS::shape_function_3D(shapefcts[gp], r, s, t, DRT::Element::tet10);
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 10 Gauss Points                  |
 *----------------------------------------------------------------------*/
const std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
DRT::ELEMENTS::So_tet10::so_tet10_11gp_derivs()
{
  static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs(
      NUMGPT_MASS_SOTET10);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
  {
    so_tet10_derivs<CORE::DRT::UTILS::GaussRule3D::tet_11point>(derivs[gp], gp);
  }
  derivs_done = true;

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 10 Gauss Points                           |
 *----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_11gp_weights()
{
  static std::vector<double> weights(NUMGPT_MASS_SOTET10);
  static bool weights_done = false;
  if (weights_done) return weights;

  const CORE::DRT::UTILS::GaussRule3D gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_11point;
  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++) weights[gp] = intpoints.qwgt[gp];
  weights_done = true;

  return weights;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::So_tet10*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet10* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)            |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::DefGradient(const std::vector<double>& disp,
    CORE::LINALG::SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();

  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
  UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(disp, xdisp);
  // update element geometry

  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    CORE::LINALG::Matrix<3, 3> invJdef;
    prestress.StoragetoMatrix(gp, invJdef, prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_xyz;
    N_xyz.Multiply(invJdef, derivs[gp]);

    // build defgrd (independent of xrefe!)
    CORE::LINALG::Matrix<3, 3> defgrd;
    defgrd.MultiplyTT(xdisp, N_xyz);
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;

    prestress.MatrixtoStorage(gp, defgrd, gpdefgrd);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected)          |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::UpdateJacobianMapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();

  // get incremental disp
  CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
  UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(disp, xdisp);

  CORE::LINALG::Matrix<3, 3> invJhist;
  CORE::LINALG::Matrix<3, 3> invJ;
  CORE::LINALG::Matrix<3, 3> defgrd;
  CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_xyz;
  CORE::LINALG::Matrix<3, 3> invJnew;
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp, invJhist, prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist, derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.MultiplyTT(xdisp, N_xyz);
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd, invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp, invJnew, prestress.JHistory());
  }  // for (int gp=0; gp<NUMGPT_SOTET10; ++gp)

  return;
}

/*----------------------------------------------------------------------*
 |  Update material                                                     |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Update_element(std::vector<double>& disp,
    Teuchos::ParameterList& params, const Teuchos::RCP<MAT::Material>& mat)
{
  if (SolidMaterial()->UsesExtendedUpdate())
  {
    // in a first step ommit everything with prestress and EAS!!
    /* ============================================================================*
    ** CONST DERIVATIVES of shape functions for TET_10 with 4 GAUSS POINTS*
    ** ============================================================================*/
    const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
        so_tet10_4gp_derivs();

    // update element geometry

    CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
    CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xcurr;  // current  coord. of element
    CORE::LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
    UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(Nodes(), xrefe);
    UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(disp, xdisp);
    UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet10, 3>(xrefe, xdisp, xcurr);

    /* =========================================================================*/
    /* ================================================= Loop over Gauss Points */
    /* =========================================================================*/

    // build deformation gradient wrt to material configuration
    CORE::LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
    params.set<int>("numgp", static_cast<int>(NUMGPT_SOTET10));
    for (unsigned gp = 0; gp < NUMGPT_SOTET10; ++gp)
    {
      CORE::LINALG::Matrix<1, NUMDIM_SOTET10> point(true);
      point.MultiplyTN(so_tet10_4gp_shapefcts()[gp], xrefe);
      params.set("gprefecoord", point);
      CORE::LINALG::Matrix<3, 10> derivs(false);
      so_tet10_derivs<CORE::DRT::UTILS::GaussRule3D::tet_4point>(derivs, gp);

      UTILS::ComputeDeformationGradient<DRT::Element::tet10>(
          defgrd, kintype_, xdisp, xcurr, invJ_[gp], derivs, pstype_, prestress_, gp);


      // call material update if material = m_growthremodel_elasthyper (calculate and update
      // inelastic deformation gradient)
      if (SolidMaterial()->UsesExtendedUpdate())
      {
        SolidMaterial()->Update(defgrd, gp, params, Id());
      }
    }
  }
  // Update of history for materials
  SolidMaterial()->Update();
}
