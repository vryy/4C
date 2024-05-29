/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluate routines for Solid Hex8 element

\level 1


*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_elements_jacobian.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_gauss_point_extrapolation.hpp"
#include "4C_discretization_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_parameter_timint.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_constraintmixture.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_growthremodel_elasthyper.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_so3_defines.hpp"
#include "4C_so3_element_service.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_hex8_determinant_analysis.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_MultiVector.h>
#include <impl/Kokkos_Traits.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

using VoigtMapping = CORE::LINALG::VOIGT::IndexMappings;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material post_setup() routine has already been called and call it if
  // not
  ensure_material_post_setup(params);

  set_params_interface_ptr(params);

  CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec2(elevec2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  CORE::Elements::ActionType act = CORE::Elements::none;

  if (IsParamsInterface())
  {
    act = params_interface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    act = CORE::Elements::String2ActionType(action);
  }


  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case CORE::Elements::struct_calc_nlnstiff:
    case CORE::Elements::struct_calc_linstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, matptr, nullptr, &elevec1,
          nullptr, &elevec3, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none, INPAR::STR::strain_none);

      break;
    }
    //==================================================================================
    // internal force vector only
    case CORE::Elements::struct_calc_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> myemat(true);

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &myemat, nullptr, &elevec1,
          nullptr, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none, INPAR::STR::strain_none);

      break;
    }
    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case CORE::Elements::struct_calc_nlnstiffmass:
    case CORE::Elements::struct_calc_nlnstifflmass:
    case CORE::Elements::struct_calc_linstiffmass:
    case CORE::Elements::struct_calc_internalinertiaforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      Teuchos::RCP<const Epetra_Vector> acc = discretization.GetState("acceleration");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
      if (acc == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'acceleration'");

      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myvel(lm.size());
      CORE::FE::ExtractMyValues(*vel, myvel, lm);
      std::vector<double> myacc(lm.size());
      CORE::FE::ExtractMyValues(*acc, myacc, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);

      // This matrix is used in the evaluation functions to store the mass matrix. If the action
      // type is CORE::Elements::struct_calc_internalinertiaforce we do not want to actually
      // populate the elemat2 variable, since the inertia terms will be directly added to the right
      // hand side. Therefore, a view is only set in cases where the evaluated mass matrix should
      // also be exported in elemat2.
      CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> mass_matrix_evaluate;
      if (act != CORE::Elements::struct_calc_internalinertiaforce)
        mass_matrix_evaluate.SetView(elemat2);

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      if (act == CORE::Elements::struct_calc_internalinertiaforce)
      {
        nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, nullptr, &mass_matrix_evaluate,
            &elevec1, &elevec2, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
      else  // standard analysis
      {
        nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &mass_matrix_evaluate,
            &elevec1, &elevec2, &elevec3, nullptr, nullptr, nullptr, params,
            INPAR::STR::stress_none, INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
      if (act == CORE::Elements::struct_calc_nlnstifflmass) soh8_lumpmass(&elemat2);

      INPAR::STR::MassLin mass_lin = INPAR::STR::MassLin::ml_none;
      auto modelevaluator_data =
          Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::Data>(ParamsInterfacePtr());
      if (modelevaluator_data != Teuchos::null)
        mass_lin = modelevaluator_data->SDyn().GetMassLinType();
      if (mass_lin == INPAR::STR::MassLin::ml_rotations)
      {
        // In case of Lie group time integration, we need to explicitly add the inertia terms to the
        // force vector, as the global mass matrix is never multiplied with the global acceleration
        // vector.
        CORE::LINALG::Matrix<NUMDOF_SOH8, 1> acceleration(true);
        for (unsigned int i_dof = 0; i_dof < NUMDOF_SOH8; i_dof++)
          acceleration(i_dof) = myacc[i_dof];
        CORE::LINALG::Matrix<NUMDOF_SOH8, 1> internal_inertia(true);
        internal_inertia.Multiply(mass_matrix_evaluate, acceleration);
        elevec2 += internal_inertia;
      }

      break;
    }
    //==================================================================================
    case CORE::Elements::struct_calc_mass_volume:
    {
      // declaration of variables
      double volume_ref = 0.0;
      double volume_mat = 0.0;
      double volume_cur = 0.0;
      double mass_ref = 0.0;
      double mass_mat = 0.0;
      double mass_cur = 0.0;

      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
      }
      // reference and current geometry (nodal positions)
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xref;  // reference coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcur;  // current  coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xmat;  // mat  coord. of element

      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xref(k, 0) = Nodes()[k]->X()[0];
        xref(k, 1) = Nodes()[k]->X()[1];
        xref(k, 2) = Nodes()[k]->X()[2];
        xcur(k, 0) = xref(k, 0) + mydisp[k * NODDOF_SOH8 + 0];
        xcur(k, 1) = xref(k, 1) + mydisp[k * NODDOF_SOH8 + 1];
        xcur(k, 2) = xref(k, 2) + mydisp[k * NODDOF_SOH8 + 2];

        // material displacements for structure with ale
        if (structale_ == true)
        {
          xmat(k, 0) = xref(k, 0) + mydispmat[k * NODDOF_SOH8 + 0];
          xmat(k, 1) = xref(k, 1) + mydispmat[k * NODDOF_SOH8 + 1];
          xmat(k, 2) = xref(k, 2) + mydispmat[k * NODDOF_SOH8 + 2];
        }
      }

      // safety check before the actual evaluation starts
      const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcur);
      if (min_detJ_curr <= 0.0)
      {
        soh8_error_handling(
            min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
        elevec1_epetra(0) = 0.0;
        return 1;
      }

      const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
      const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs =
          soh8_derivs();
      const static std::vector<double> gpweights = soh8_weights();

      // MAT ------------------------
      // build new jacobian mapping with respect to the material configuration
      int err = 0;
      if (structale_ == true)
      {
        err = init_jacobian_mapping(mydispmat);
        if (err)
        {
          // reset class variable before leaving
          init_jacobian_mapping();
          return err;
        }
      }

      std::vector<double> detJmat = detJ_;
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJmat = invJ_;

      err = init_jacobian_mapping(mydisp);
      if (err)
      {
        // reset class variable before leaving
        init_jacobian_mapping();
        return err;
      }

      std::vector<double> detJcur = detJ_;
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJcur = invJ_;

      init_jacobian_mapping();

      std::vector<double> detJref = detJ_;
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJref = invJ_;


      CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
      CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);

      //----------------------------------------------------------------
      // loop over all Gauss points
      //----------------------------------------------------------------
      for (unsigned ip = 0; ip < NUMGPT_SOH8; ++ip)
      {
        const double density = Material()->Density(ip);
        const double wgt = gpweights[ip];

        /*------------------------------------ integration factor  -------*/
        double fac = wgt * detJref[ip];
        volume_ref += fac;
        fac = wgt * detJref[ip] * density;
        mass_ref += fac;

        // MAT ------------------------
        if (structale_)
        {
          fac = wgt * detJmat[ip];
          volume_mat += fac;
          fac = wgt * detJmat[ip] * density;
          mass_mat += fac;

          N_XYZ.Multiply(invJmat[ip], derivs[ip]);
          defgrd.MultiplyTT(xcur, N_XYZ);
          double detFmat = defgrd.Determinant();

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detJcur[ip];
          volume_cur += fac;
          fac = wgt * detJcur[ip] * density * 1 / detFmat;
          mass_cur += fac;
        }
        else
        {
          N_XYZ.Multiply(invJref[ip], derivs[ip]);
          defgrd.MultiplyTT(xcur, N_XYZ);
          double detFref = defgrd.Determinant();

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detJcur[ip];
          volume_cur += fac;
          fac = wgt * detJcur[ip] * density * 1 / detFref;
          mass_cur += fac;
        }
      }

      //----------------------------------------------------------------

      // return results
      if (!structale_)
      {
        volume_mat = volume_ref;
        mass_mat = mass_ref;
      }

      elevec1(0) = volume_ref;
      elevec1(1) = volume_mat;
      elevec1(2) = volume_cur;
      elevec1(3) = mass_ref;
      elevec1(4) = mass_mat;
      elevec1(5) = mass_cur;

      break;
    }

    case CORE::Elements::analyse_jacobian_determinant:
    {
      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);
      // reference and current geometry (nodal positions)
      CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> xcurr;  // current  coord. of element

      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xcurr(0, k) = Nodes()[k]->X()[0] + mydisp[k * NODDOF_SOH8 + 0];
        xcurr(1, k) = Nodes()[k]->X()[1] + mydisp[k * NODDOF_SOH8 + 1];
        xcurr(2, k) = Nodes()[k]->X()[2] + mydisp[k * NODDOF_SOH8 + 2];
      }

      Teuchos::RCP<SoHex8DeterminantAnalysis> det_analyser = SoHex8DeterminantAnalysis::create();
      if (not det_analyser->isValid(xcurr))
        soh8_error_handling(-1.0, params, __LINE__, STR::ELEMENTS::ele_error_determinant_analysis);

      break;
    }

    //==================================================================================
    // recover elementwise stored quantities (e.g. EAS)
    case CORE::Elements::struct_calc_recover:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");

      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW(
            "Cannot get state vectors \"displacement\" "
            "and/or \"residual displacement\"");

      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);

      soh8_recover(lm, myres);
      /* ToDo Probably we have to recover the history information of some special
       * materials as well.                                 hiermeier 04/2016  */

      break;
    }
    //==================================================================================
    // evaluate stresses and strains at gauss points
    case CORE::Elements::struct_calc_stress:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
      Teuchos::RCP<std::vector<char>> straindata = Teuchos::null;
      Teuchos::RCP<std::vector<char>> plstraindata = Teuchos::null;
      INPAR::STR::StressType iostress = INPAR::STR::stress_none;
      INPAR::STR::StrainType iostrain = INPAR::STR::strain_none;
      INPAR::STR::StrainType ioplstrain = INPAR::STR::strain_none;
      if (IsParamsInterface())
      {
        stressdata = str_params_interface().StressDataPtr();
        straindata = str_params_interface().StrainDataPtr();
        plstraindata = str_params_interface().plastic_strain_data_ptr();

        iostress = str_params_interface().GetStressOutputType();
        iostrain = str_params_interface().GetStrainOutputType();
        ioplstrain = str_params_interface().get_plastic_strain_output_type();
      }
      else
      {
        stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        straindata = params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
            params, "iostress", INPAR::STR::stress_none);
        iostrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
            params, "iostrain", INPAR::STR::strain_none);
        // in case of small strain materials calculate plastic strains for post processing
        plstraindata = params.get<Teuchos::RCP<std::vector<char>>>("plstrain", Teuchos::null);
        ioplstrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
            params, "ioplstrain", INPAR::STR::strain_none);
      }
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
      if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");
      if (plstraindata == Teuchos::null) FOUR_C_THROW("Cannot get 'plastic strain' data");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> plstrain;

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr, nullptr,
          nullptr, nullptr, &stress, &strain, &plstrain, params, iostress, iostrain, ioplstrain);

      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, stress);
        data.StartPacking();
        AddtoPack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }

      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, strain);
        data.StartPacking();
        AddtoPack(data, strain);
        std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
      }

      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, plstrain);
        data.StartPacking();
        AddtoPack(data, plstrain);
        std::copy(data().begin(), data().end(), std::back_inserter(*plstraindata));
      }
    }
    break;
    case CORE::Elements::struct_init_gauss_point_data_output:
    {
      FOUR_C_ASSERT(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Save number of Gauss of the element for gauss point data output
      str_params_interface()
          .gauss_point_data_output_manager_ptr()
          ->add_element_number_of_gauss_points(NUMGPT_SOH8);

      // holder for output quantity names and their size
      std::unordered_map<std::string, int> quantities_map{};

      // Ask material for the output quantity names and sizes
      SolidMaterial()->register_output_data_names(quantities_map);

      // Add quantities to the Gauss point output data manager (if they do not already exist)
      str_params_interface().gauss_point_data_output_manager_ptr()->MergeQuantities(quantities_map);
    }
    break;
    case CORE::Elements::struct_gauss_point_data_output:
    {
      FOUR_C_ASSERT(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Collection and assembly of gauss point data
      for (const auto& quantity :
          str_params_interface().gauss_point_data_output_manager_ptr()->GetQuantities())
      {
        const std::string& quantity_name = quantity.first;
        const int quantity_size = quantity.second;

        // Step 1: Collect the data for each Gauss point for the material
        CORE::LINALG::SerialDenseMatrix gp_data(NUMGPT_SOH8, quantity_size, true);
        bool data_available = SolidMaterial()->EvaluateOutputData(quantity_name, gp_data);

        // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
        // point)
        if (data_available)
        {
          switch (str_params_interface().gauss_point_data_output_manager_ptr()->GetOutputType())
          {
            case INPAR::STR::GaussPointDataOutputType::element_center:
            {
              // compute average of the quantities
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  str_params_interface()
                      .gauss_point_data_output_manager_ptr()
                      ->get_element_center_data()
                      .at(quantity_name);
              CORE::FE::AssembleAveragedElementValues(*global_data, gp_data, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::nodes:
            {
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  str_params_interface().gauss_point_data_output_manager_ptr()->GetNodalData().at(
                      quantity_name);

              Epetra_IntVector& global_nodal_element_count =
                  *str_params_interface()
                       .gauss_point_data_output_manager_ptr()
                       ->GetNodalDataCount()
                       .at(quantity_name);

              static auto gauss_integration = CORE::FE::IntegrationPoints3D(
                  CORE::FE::NumGaussPointsToGaussRule<CORE::FE::CellType::hex8>(NUMGPT_SOH8));
              CORE::FE::ExtrapolateGPQuantityToNodesAndAssemble<CORE::FE::CellType::hex8>(
                  *this, gp_data, *global_data, false, gauss_integration);
              DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::gauss_points:
            {
              std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                  str_params_interface()
                      .gauss_point_data_output_manager_ptr()
                      ->GetGaussPointData()
                      .at(quantity_name);
              DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::none:
              FOUR_C_THROW(
                  "You specified a Gauss point data output type of none, so you should not end up "
                  "here.");
            default:
              FOUR_C_THROW("Unknown Gauss point data output type.");
          }
        }
      }
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_eleload:
      FOUR_C_THROW("this method is not supposed to evaluate a load, use evaluate_neumann(...)");
      break;
    //==================================================================================
    case CORE::Elements::struct_calc_fsiload:
      FOUR_C_THROW("Case not yet implemented");
      break;
    //==================================================================================
    case CORE::Elements::struct_calc_update_istep:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      update_element(mydisp, params, Material());
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_reset_istep:
    {
      // restore EAS parameters
      if (eastype_ != soh8_easnone)
      {
        soh8_easrestore();

        // reset EAS internal force
        CORE::LINALG::SerialDenseMatrix* oldfeas = &easdata_.feas;
        oldfeas->putScalar(0.0);
      }
      // Reset of history (if needed)
      SolidMaterial()->reset_step();
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_store_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) FOUR_C_THROW("Provide timestep number to be stored");

      // EAS
      if (eastype_ != soh8_easnone) soh8_easupdate();

      // Material
      SolidMaterial()->StoreHistory(timestep);
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_recover_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) FOUR_C_THROW("Provide timestep number of the timestep to be recovered");

      // EAS
      if (eastype_ != soh8_easnone) soh8_easrestore();

      // Material
      SolidMaterial()->SetHistory(timestep);
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_energy:
    {
      // initialization of internal energy
      double intenergy = 0.0;

      // shape functions and Gauss weights
      const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs =
          soh8_derivs();
      const static std::vector<double> weights = soh8_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      // update element geometry
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;

      UTILS::EvaluateNodalCoordinates<CORE::FE::CellType::hex8, 3>(Nodes(), xrefe);
      UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(mydisp, xdisp);
      UTILS::EvaluateCurrentNodalCoordinates<CORE::FE::CellType::hex8, 3>(xrefe, xdisp, xcurr);

      // safety check before the actual evaluation starts
      const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
      if (min_detJ_curr <= 0.0)
      {
        soh8_error_handling(
            min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
        elevec1_epetra(0) = 0.0;
        return 0;
      }

      // prepare EAS data
      CORE::LINALG::SerialDenseMatrix* alpha = nullptr;              // EAS alphas
      std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
      CORE::LINALG::SerialDenseMatrix M;                             // EAS matrix M at current GP
      double detJ0;                                                  // detJ(origin)
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix
      if (eastype_ != soh8_easnone)
      {
        alpha = &easdata_.alpha;  // get alpha of previous iteration
        soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);
      }

      // loop over all Gauss points
      for (unsigned gp = 0; gp < NUMGPT_SOH8; gp++)
      {
        // Gauss weights and Jacobian determinant
        double fac = detJ_[gp] * weights[gp];

        /* get the inverse of the Jacobian matrix which looks like:
        **            [ x_,r  y_,r  z_,r ]^-1
        **     J^-1 = [ x_,s  y_,s  z_,s ]
        **            [ x_,t  y_,t  z_,t ]
        */
        // compute derivatives N_XYZ at gp w.r.t. material coordinates
        // by N_XYZ = J^-1 * N_rst
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ(true);
        N_XYZ.Multiply(invJ_[gp], derivs[gp]);

        // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(true);

        if (PRESTRESS::IsMulf(pstype_))
        {
          // get Jacobian mapping wrt to the stored configuration
          CORE::LINALG::Matrix<3, 3> invJdef;
          prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
          // get derivatives wrt to last spatial configuration
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
          N_xyz.Multiply(invJdef, derivs[gp]);

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

          // right Cauchy-Green tensor = F^T * F
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
          cauchygreen.MultiplyTN(defgrd, defgrd);

          glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
          glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
          glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
          glstrain(3) = cauchygreen(0, 1);
          glstrain(4) = cauchygreen(1, 2);
          glstrain(5) = cauchygreen(2, 0);
        }
        else if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
        {
          // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
          defgrd.MultiplyTT(xcurr, N_XYZ);

          // right Cauchy-Green tensor = F^T * F
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
          cauchygreen.MultiplyTN(defgrd, defgrd);

          glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
          glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
          glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
          glstrain(3) = cauchygreen(0, 1);
          glstrain(4) = cauchygreen(1, 2);
          glstrain(5) = cauchygreen(2, 0);
        }
        else if (kintype_ == INPAR::STR::KinemType::linear)
        {
          // in kinematically linear analysis the deformation gradient is equal to identity
          // no difference between reference and current state
          for (int i = 0; i < 3; ++i) defgrd(i, i) = 1.0;

          // nodal displacement vector
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1> nodaldisp;
          for (int i = 0; i < NUMDOF_SOH8; ++i) nodaldisp(i, 0) = mydisp[i];
          // compute linear B-operator
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
          for (int i = 0; i < NUMNOD_SOH8; ++i)
          {
            bop(0, NODDOF_SOH8 * i + 0) = N_XYZ(0, i);
            bop(0, NODDOF_SOH8 * i + 1) = 0.0;
            bop(0, NODDOF_SOH8 * i + 2) = 0.0;
            bop(1, NODDOF_SOH8 * i + 0) = 0.0;
            bop(1, NODDOF_SOH8 * i + 1) = N_XYZ(1, i);
            bop(1, NODDOF_SOH8 * i + 2) = 0.0;
            bop(2, NODDOF_SOH8 * i + 0) = 0.0;
            bop(2, NODDOF_SOH8 * i + 1) = 0.0;
            bop(2, NODDOF_SOH8 * i + 2) = N_XYZ(2, i);

            bop(3, NODDOF_SOH8 * i + 0) = N_XYZ(1, i);
            bop(3, NODDOF_SOH8 * i + 1) = N_XYZ(0, i);
            bop(3, NODDOF_SOH8 * i + 2) = 0.0;
            bop(4, NODDOF_SOH8 * i + 0) = 0.0;
            bop(4, NODDOF_SOH8 * i + 1) = N_XYZ(2, i);
            bop(4, NODDOF_SOH8 * i + 2) = N_XYZ(1, i);
            bop(5, NODDOF_SOH8 * i + 0) = N_XYZ(2, i);
            bop(5, NODDOF_SOH8 * i + 1) = 0.0;
            bop(5, NODDOF_SOH8 * i + 2) = N_XYZ(0, i);
          }

          // compute linear strain at GP
          glstrain.Multiply(bop, nodaldisp);
        }
        else
          FOUR_C_THROW("unknown kinematic type for energy calculation");

        // EAS technology: "enhance the strains"  ----------------------------- EAS
        if (eastype_ != soh8_easnone)
        {
          M.shape(MAT::NUM_STRESS_3D, neas_);
          // map local M to global, also enhancement is referred to element origin
          // M = detJ0/detJ T0^{-T} . M
          // CORE::LINALG::SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
          // add enhanced strains = M . alpha to GL strains to "unlock" element
          switch (eastype_)
          {
            case DRT::ELEMENTS::SoHex8::soh8_easfull:
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_easfull>(M.values(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).values());
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
                  1.0, glstrain.A(), 1.0, M.values(), alpha->values());
              break;
            case DRT::ELEMENTS::SoHex8::soh8_easmild:
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_easmild>(M.values(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).values());
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
                  1.0, glstrain.A(), 1.0, M.values(), alpha->values());
              break;
            case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_eassosh8>(
                  M.values(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).values());
              CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
                  1.0, glstrain.A(), 1.0, M.values(), alpha->values());
              break;
            case DRT::ELEMENTS::SoHex8::soh8_easnone:
              break;
            default:
              FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
              break;
          }
        }  // ------------------------------------------------------------------ EAS

        if (defgrd.Determinant() <= 0.0)
        {
          if (IsParamsInterface() and str_params_interface().IsTolerateErrors())
          {
            str_params_interface().SetEleEvalErrorFlag(
                STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
            return 0;
          }
          else
          {
            FOUR_C_THROW("Negative deformation gradient!");
          }
        }

        // call material for evaluation of strain energy function
        double psi = 0.0;
        SolidMaterial()->StrainEnergy(glstrain, psi, gp, Id());

        // sum up GP contribution to internal energy
        intenergy += fac * psi;
      }

      if (IsParamsInterface())  // new structural time integration
      {
        str_params_interface().add_contribution_to_energy_type(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        // check length of elevec1
        if (elevec1_epetra.length() < 1) FOUR_C_THROW("The given result vector is too short.");

        elevec1_epetra(0) = intenergy;
      }
    }
    break;
    //==================================================================================
    case CORE::Elements::multi_calc_dens:
    {
      soh8_homog(params);
    }
    break;
      //==================================================================================
    case CORE::Elements::struct_interpolate_velocity_to_point:
    {
      // get displacements and extract values of this element (set in prepare_fluid_op())
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velocity");

#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (dispnp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");
      if (velnp == Teuchos::null) FOUR_C_THROW("Cannot get state velocity vector");
#endif

      std::vector<double> mydispnp(lm.size());
      CORE::FE::ExtractMyValues(*dispnp, mydispnp, lm);

      std::vector<double> myvelnp(lm.size());
      CORE::FE::ExtractMyValues(*velnp, myvelnp, lm);

      // update element geometry
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>
          xdisp;  // current displacements of element nodes
      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current coord. of element
      UTILS::EvaluateNodalCoordinates<CORE::FE::CellType::hex8, 3>(Nodes(), xrefe);
      UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(mydispnp, xdisp);
      UTILS::EvaluateCurrentNodalCoordinates<CORE::FE::CellType::hex8, 3>(xrefe, xdisp, xcurr);

      // shape functions and derivatives w.r.t. r,s,t
      CORE::LINALG::Matrix<NUMNOD_SOH8, 1> shapefcts;
      CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> deriv1;
      // coordinates of given point in reference coordinates
      CORE::LINALG::Matrix<NUMDIM_SOH8, 1> xsi;
      xsi(0) = elevec2_epetra(0);
      xsi(1) = elevec2_epetra(1);
      xsi(2) = elevec2_epetra(2);
      // evaluate shape functions and derivatives at given point w.r.t r,s,t
      CORE::FE::shape_function<CORE::FE::CellType::hex8>(xsi, shapefcts);
      CORE::FE::shape_function_deriv1<CORE::FE::CellType::hex8>(xsi, deriv1);

      CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> myvelocitynp;
      for (int node = 0; node < NUMNOD_SOH8; ++node)
      {
        for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
        {
          myvelocitynp(node, dim) = myvelnp[node * 3 + dim];
        }
      }

      if (params.get("calculate_velocity", 1))
      {
        //************************************************************************
        // 1.) interpolation of velocity at n+1 to given point
        //************************************************************************

        // give back velocity at given point
        CORE::LINALG::Matrix<3, 1> result;
        result.MultiplyTN(myvelocitynp, shapefcts);
        for (int i = 0; i < 3; ++i) elevec1_epetra(i) = result(i, 0);
      }
      else
      {
        //************************************************************************
        // 2.) calculation of divergence of structural velocity
        //************************************************************************

        // get Jacobian matrix
        // actually compute its transpose....
        /*
                 +-            -+ -1
                 | dx   dy   dz |
                 | --   --   -- |
                 | dr   dr   dr |
                 |              |
                 | dx   dy   dz |
        J^{-T} = | --   --   -- |
                 | ds   ds   ds |
                 |              |
                 | dx   dy   dz |
                 | --   --   -- |
                 | dt   dt   dt |
                 +-            -+
        */

        // global first derivatives of shape functions w.r.t x,y,z (by derxy1 = J^-T * deriv1)
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> derxy1;
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> transJ;
        transJ.Multiply(deriv1, xrefe);
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ(transJ);
        invJ.Invert();
        derxy1.Multiply(invJ, deriv1);

        // build (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * derxy1^T
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrdnp(true);
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrdnp_inv(true);
        defgrdnp.MultiplyTT(xcurr, derxy1);
        defgrdnp_inv.Invert(defgrdnp);

        // evaluate divergence of structural velocity
        double velocitydivergence = 0.0;

        // div(v) = derxy1_A,I (dX_I/dx_i) v_Ai
        for (int I = 0; I < NUMDIM_SOH8; ++I)
        {
          for (int A = 0; A < NUMNOD_SOH8; ++A)
          {
            for (int i = 0; i < NUMDIM_SOH8; ++i)
            {
              velocitydivergence += derxy1(I, A) * (defgrdnp_inv(I, i) * myvelocitynp(A, i));
            }
          }
        }

        // vector to fill
        elevec1_epetra(3) = velocitydivergence;
      }
    }
    break;
    //==================================================================================
    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // allocation and initializiation of these data arrays can only be
    // done in the elements that know the number of EAS parameters
    case CORE::Elements::multi_init_eas:
    {
      soh8_eas_init_multi(params);
    }
    break;
    //==================================================================================
    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // before any microscale simulation, EAS internal data has to be
    // set accordingly
    case CORE::Elements::multi_set_eas:
    {
      soh8_set_eas_multi(params);
    }
    break;
    //==================================================================================
    // read restart of microscale
    case CORE::Elements::multi_readrestart:
    {
      soh8_read_restart_multi();
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_update_prestress:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      switch (pstype_)
      {
        case INPAR::STR::PreStress::mulf:
        {
          // build def gradient for every gauss point
          CORE::LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOH8, 9);
          def_gradient(mydisp, gpdefgrd, *prestress_);

          // update deformation gradient and put back to storage
          CORE::LINALG::Matrix<3, 3> deltaF;
          CORE::LINALG::Matrix<3, 3> Fhist;
          CORE::LINALG::Matrix<3, 3> Fnew;
          for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
          {
            prestress_->StoragetoMatrix(gp, deltaF, gpdefgrd);
            prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());
            Fnew.Multiply(deltaF, Fhist);
            prestress_->MatrixtoStorage(gp, Fnew, prestress_->FHistory());
          }

          // push-forward invJ for every gaussian point
          update_jacobian_mapping(mydisp, *prestress_);

          // Update constraintmixture material
          if (Material()->MaterialType() == CORE::Materials::m_constraintmixture)
          {
            SolidMaterial()->Update();
          }
          break;
        }
        default:
          FOUR_C_THROW(
              "You should either not be here, or the prestressing method you are using is not "
              "implemented for HEX8 elements!");
      }
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss points and store gpstresses in map <EleId, gpstresses
    // >
    case CORE::Elements::struct_calc_global_gpstresses_map:
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
        Teuchos::RCP<std::vector<char>> plstraindata =
            params.get<Teuchos::RCP<std::vector<char>>>("plstrain", Teuchos::null);
        if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
        if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");
        if (plstraindata == Teuchos::null) FOUR_C_THROW("Cannot get 'plastic strain' data");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstressmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstressmap", Teuchos::null);
        if (gpstressmap == Teuchos::null)
          FOUR_C_THROW("no gp stress map available for writing gpstresses");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstrainmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstrainmap", Teuchos::null);
        if (gpstrainmap == Teuchos::null)
          FOUR_C_THROW("no gp strain map available for writing gpstrains");
        std::vector<double> mydisp(lm.size());
        CORE::FE::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        CORE::FE::ExtractMyValues(*res, myres, lm);

        std::vector<double> mydispmat(lm.size(), 0.0);
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
        CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
        CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> plstrain;
        auto iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
            params, "iostress", INPAR::STR::stress_none);
        auto iostrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
            params, "iostrain", INPAR::STR::strain_none);
        auto ioplstrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
            params, "ioplstrain", INPAR::STR::strain_none);

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr, nullptr,
            nullptr, nullptr, &stress, &strain, &plstrain, params, iostress, iostrain, ioplstrain);

        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstress =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstress->shape(NUMGPT_SOH8, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (unsigned i = 0; i < NUMGPT_SOH8; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstrain->shape(NUMGPT_SOH8, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (unsigned i = 0; i < NUMGPT_SOH8; i++)
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
          CORE::COMM::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }

        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }

        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, plstrain);
          data.StartPacking();
          AddtoPack(data, plstrain);
          std::copy(data().begin(), data().end(), std::back_inserter(*plstraindata));
        }
      }
    }
    break;
    //==================================================================================
    case CORE::Elements::struct_calc_predict:
    {
      // do nothing here
      break;
    }
    //==================================================================================
    // create a backup state for all internally stored variables (e.g. EAS)
    case CORE::Elements::struct_create_backup:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res.is_null()) FOUR_C_THROW("Cannot get state vector \"residual displacement\"");

      // extract the part for this element
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);

      soh8_create_eas_backup_state(myres);

      break;
    }
    //==================================================================================
    /* recover internally stored state variables from a previously created backup
     * state (e.g. EAS) */
    case CORE::Elements::struct_recover_from_backup:
    {
      soh8_recover_from_eas_backup_state();

      break;
    }
    default:
      FOUR_C_THROW(
          "Unknown type of action for So_hex8: %s", CORE::Elements::ActionType2String(act).c_str());
      break;
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8::evaluate_neumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  set_params_interface_ptr(params);
  // get values and switches from the condition
  const auto* onoff = &condition.parameters().Get<std::vector<int>>("onoff");
  const auto* val = &condition.parameters().Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = params_interface().GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOH8)
    FOUR_C_THROW("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOH8; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = &condition.parameters().Get<std::vector<int>>("funct");
  CORE::LINALG::Matrix<NUMDIM_SOH8, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOH8; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;


  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  UTILS::EvaluateNodalCoordinates<CORE::FE::CellType::hex8, 3>(Nodes(), xrefe);
  /* ================================================= Loop over Gauss Points */
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // compute the Jacobian matrix
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac;
    jac.Multiply(derivs[gp], xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOH8; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOH8; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = gpweights[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOH8; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0) ? GLOBAL::Problem::Instance()
                                 ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .Evaluate(xrefegp.A(), time, dim)
                           : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SOH8; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOH8 + dim] += shapefcts[gp](nodid) * dim_fac;
        }
      }
    }

  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_hex8::evaluate_neumann

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const double* DRT::ELEMENTS::SoHex8::soh8_get_coordinate_of_gausspoints(const unsigned dim) const
{
  static CORE::LINALG::Matrix<NUMGPT_SOH8, NUMDIM_SOH8> coordinates_of_gps(false);
  static bool init = false;

  if (not init)
  {
    if (gp_rule_.NumPoints() != NUMGPT_SOH8)
      FOUR_C_THROW(
          "Inconsistent number of GPs: "
          "%d != %d",
          gp_rule_.NumPoints(), NUMGPT_SOH8);

    for (unsigned gp = 0; gp < gp_rule_.NumPoints(); ++gp)
      for (unsigned d = 0; d < NUMDIM_SOH8; ++d) coordinates_of_gps(gp, d) = gp_rule_.Point(gp)[d];
    // do it only once
    init = true;
  }

  return &coordinates_of_gps(0, dim);
}

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::init_jacobian_mapping()
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    CORE::Nodes::Node** nodes = Nodes();
    if (!nodes) FOUR_C_THROW("Nodes() returned null pointer");
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(derivs[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) FOUR_C_THROW("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);

    if (PRESTRESS::IsMulfActive(time_, pstype_, pstime_))
    {
      if (!(prestress_->is_init()))
        prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());
    }
  }

  if (PRESTRESS::IsMulfActive(time_, pstype_, pstime_)) prestress_->is_init() = true;
}
/*----------------------------------------------------------------------*
 |  init the element jacobian mapping with respect to the    farah 06/13|
 |  material configuration.                                             |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8::init_jacobian_mapping(std::vector<double>& dispmat)
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xmat;

  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    CORE::Nodes::Node** nodes = Nodes();
    if (!nodes) FOUR_C_THROW("Nodes() returned null pointer");

    xmat(i, 0) = Nodes()[i]->X()[0] + dispmat[i * NODDOF_SOH8 + 0];
    xmat(i, 1) = Nodes()[i]->X()[1] + dispmat[i * NODDOF_SOH8 + 1];
    xmat(i, 2) = Nodes()[i]->X()[2] + dispmat[i * NODDOF_SOH8 + 2];
  }
  invJ_.clear();
  detJ_.clear();
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(derivs[gp], xmat);
    detJ_[gp] = invJ_[gp].Invert();

    if (detJ_[gp] <= 0.0)
    {
      if (IsParamsInterface() and str_params_interface().IsTolerateErrors())
      {
        str_params_interface().SetEleEvalErrorFlag(
            STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
        return 1;
      }
      else
        FOUR_C_THROW("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::SoHex8::soh8_get_min_det_jac_at_corners(
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr) const
{
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> xcurr_t(false);
  xcurr_t.UpdateT(xcurr);
  return CORE::Elements::GetMinimalJacDeterminantAtNodes<CORE::FE::CellType::hex8>(xcurr_t);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_error_handling(const double& det_curr,
    Teuchos::ParameterList& params, const int line_id, const STR::ELEMENTS::EvalErrorFlag flag)
{
  error_handling(det_curr, params, line_id, flag);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_compute_eas_inc(
    const std::vector<double>& residual, CORE::LINALG::SerialDenseMatrix* const eas_inc)
{
  auto* oldKaainv = &easdata_.invKaa;
  auto* oldKda = &easdata_.Kda;
  auto* oldfeas = &easdata_.feas;
  if (!oldKaainv || !oldKda || !oldfeas) FOUR_C_THROW("Missing EAS history data");

  // we need the (residual) displacement at the previous step
  CORE::LINALG::SerialDenseVector res_d_eas(NUMDOF_SOH8);
  for (int i = 0; i < NUMDOF_SOH8; ++i) res_d_eas(i) = residual[i];
  // --- EAS default update ---------------------------
  CORE::LINALG::SerialDenseMatrix eashelp(neas_, 1);
  /*----------- make multiplication eashelp = oldLt * disp_incr[kstep] */
  CORE::LINALG::multiply(eashelp, *oldKda, res_d_eas);
  /*---------------------------------------- add old Rtilde to eashelp */
  eashelp += *oldfeas;
  /*--------- make multiplication alpha_inc = - old Dtildinv * eashelp */
  CORE::LINALG::multiply(*eas_inc, *oldKaainv, eashelp);
  eas_inc->scale(-1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_recover(
    const std::vector<int>& lm, const std::vector<double>& residual)
{
  // for eas
  CORE::LINALG::SerialDenseMatrix* alpha = nullptr;
  CORE::LINALG::SerialDenseMatrix* eas_inc = nullptr;
  // get access to the interface parameters
  const double step_length = str_params_interface().GetStepLength();
  const bool iseas = (eastype_ != soh8_easnone);

  // have eas?
  if (iseas)
  {
    // access general eas history stuff stored in element
    // get alpha of previous iteration
    alpha = &easdata_.alpha;
    // get the old eas increment
    eas_inc = &easdata_.eas_inc;
    if (!alpha || !eas_inc) FOUR_C_THROW("Missing EAS history data (eas_inc and/or alpha)");
  }

  /* if it is a default step, we have to recover the condensed
   * solution vectors */
  if (str_params_interface().IsDefaultStep())
  {
    /* recovery of the enhanced assumed strain increment and
     * update of the eas dofs. */
    if (iseas)
    {
      // first, store the eas state of the previous accepted Newton step
      str_params_interface().sum_into_my_previous_sol_norm(
          NOX::NLN::StatusTest::quantity_eas, neas_, (*alpha)[0], Owner());

      // compute the eas increment
      soh8_compute_eas_inc(residual, eas_inc);

      /*--------------------------- update alpha += step_length * alfa_inc */
      for (int i = 0; i < neas_; ++i) (*alpha)(i, 0) += step_length * (*eas_inc)(i, 0);
    }  // if (iseas)
  }    // if (*isdefault_step_ptr_)
  /* if it is no default step, we can correct the update and the current eas
   * state without the need for any matrix-vector products. */
  else
  {
    // The first step has to be a default step!
    if (old_step_length_ < 0.0) FOUR_C_THROW("The old step length was not defined!");
    /* if this is no full step, we have to adjust the length of the
     * enhanced assumed strain incremental step. */
    if (iseas)
    {
      /* undo the previous step:
       *            alpha_new = alpha_old - old_step * alpha_inc
       * and update the solution variable with the new step length:
       *            alpha_new = alpha_new + new_step * alpha_inc */
      for (int i = 0; i < neas_; ++i)
        (*alpha)(i, 0) += (step_length - old_step_length_) * (*eas_inc)(i, 0);

      //      {
      //        std::cout << "EAS #" << Id() << ":\n";
      //        alpha->Print( std::cout );
      //        eas_inc->Print( std::cout );
      //        std::cout << "\n";
      //      }
    }  // if (nhyb_)
  }    // else

  // save the old step length
  old_step_length_ = step_length;

  // Check if the eas incr is tested and if yes, calculate the element
  // contribution to the norm
  if (iseas)
    str_params_interface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, neas_,
        (*eas_inc)[0], (*alpha)[0], step_length, Owner());

  // the element internal stuff should be up-to-date for now...
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::nlnstiffmass(std::vector<int>& lm,    // location matrix
    std::vector<double>& disp,                                    // current displacements
    std::vector<double>* vel,                                     // current velocities
    std::vector<double>* acc,                                     // current accelerations
    std::vector<double>& residual,                                // current residual displ
    std::vector<double>& dispmat,                                 // current material displacements
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,   // element mass matrix
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,                  // element internal force vector
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* forceinert,             // element inertial force vector
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,              // element structural force vector
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,    // stresses at GP
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,    // strains at GP
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
    Teuchos::ParameterList& params,           // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,    // stress output option
    const INPAR::STR::StrainType iostrain,    // strain output option
    const INPAR::STR::StrainType ioplstrain)  // plastic strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // check for prestressing
  if (PRESTRESS::IsAny(pstype_) && eastype_ != soh8_easnone)
    FOUR_C_THROW("No way you can do mulf or id prestressing with EAS turned on!");

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(false);  // reference coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr(false);  // current  coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp(false);

  UTILS::EvaluateNodalCoordinates<CORE::FE::CellType::hex8, 3>(Nodes(), xrefe);
  UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(disp, xdisp);
  UTILS::EvaluateCurrentNodalCoordinates<CORE::FE::CellType::hex8, 3>(xrefe, xdisp, xcurr);

  // safety check before the actual evaluation starts
  const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
  if (min_detJ_curr <= 0.0)
  {
    soh8_error_handling(
        min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
    return;
  }

  double elediagonallength = 0.0;
  if (!analyticalmaterialtangent_)
    elediagonallength = sqrt(pow(xrefe(0, 0) - xrefe(6, 0), 2) + pow(xrefe(0, 1) - xrefe(6, 1), 2) +
                             pow(xrefe(0, 2) - xrefe(6, 2), 2));

  // we need the (residual) displacement at the previous step
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> nodaldisp;
  for (int i = 0; i < NUMDOF_SOH8; ++i) nodaldisp(i, 0) = disp[i];

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  CORE::LINALG::SerialDenseMatrix* alpha = nullptr;              // EAS alphas
  std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
  CORE::LINALG::SerialDenseMatrix M;                             // EAS matrix M at current GP
  CORE::LINALG::SerialDenseVector feas;                          // EAS portion of internal forces
  CORE::LINALG::SerialDenseMatrix Kaa;                           // EAS matrix Kaa
  CORE::LINALG::SerialDenseMatrix Kda;                           // EAS matrix Kda
  double detJ0;                                                  // detJ(origin)
  CORE::LINALG::SerialDenseMatrix* oldfeas = nullptr;            // EAS history
  CORE::LINALG::SerialDenseMatrix* oldKaainv = nullptr;          // EAS history
  CORE::LINALG::SerialDenseMatrix* oldKda = nullptr;             // EAS history
  CORE::LINALG::SerialDenseMatrix* eas_inc = nullptr;            // EAS increment

  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix

  if (eastype_ != soh8_easnone)
  {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = &easdata_.alpha;  // get alpha of previous iteration
    oldfeas = &easdata_.feas;
    oldKaainv = &easdata_.invKaa;
    oldKda = &easdata_.Kda;
    eas_inc = &easdata_.eas_inc;

    if (!alpha || !oldKaainv || !oldKda || !oldfeas || !eas_inc)
      FOUR_C_THROW("Missing EAS history-data");

    // ============================== DEPRECATED ==============================
    // FixMe deprecated implementation
    if (not IsParamsInterface())
    {
      // we need the (residual) displacement at the previous step
      CORE::LINALG::SerialDenseVector res_d_eas(NUMDOF_SOH8);
      for (int i = 0; i < NUMDOF_SOH8; ++i) res_d_eas(i) = residual[i];

      // this is a line search step, i.e. the direction of the eas increments
      // has been calculated by a Newton step and now it is only scaled
      if (params.isParameter("alpha_ls"))
      {
        double alpha_ls = params.get<double>("alpha_ls");
        // undo step
        eas_inc->scale(-1.);
        alpha->operator+=(*eas_inc);
        // scale increment
        eas_inc->scale(-1. * alpha_ls);
        // add reduced increment
        alpha->operator+=(*eas_inc);
      }
      // add Kda . res_d to feas
      // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
      else
      {
        switch (eastype_)
        {
          case DRT::ELEMENTS::SoHex8::soh8_easfull:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_easfull, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_easfull, soh8_easfull, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_easfull, 1>(1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::SoHex8::soh8_easmild:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_easmild, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_easmild, soh8_easmild, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_easmild, 1>(1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, soh8_eassosh8, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(
                1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::SoHex8::soh8_easnone:
            break;
          default:
            FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }
    }  // if (not IsInterface())
    // ============================== DEPRECATED ==============================

    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.shape(neas_, neas_);

    // EAS matrix K_{d alpha}
    Kda.shape(neas_, NUMDOF_SOH8);

    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);
  }  // -------------------------------------------------------------------- EAS

  // build new jacobian mapping with respect to the material configuration
  if (structale_ == true) init_jacobian_mapping(dispmat);

  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = params.isParameter("cond_rhs_norm");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;

  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    double detJ = detJ_[gp];
    N_XYZ.Multiply(invJ_[gp], derivs[gp]);

    if (PRESTRESS::IsMulf(pstype_))
    {
      // get Jacobian mapping wrt to the stored configuration
      CORE::LINALG::Matrix<3, 3> invJdef;
      prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      CORE::LINALG::Matrix<3, 8> N_xyz;
      N_xyz.Multiply(invJdef, derivs[gp]);

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
    else if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
    {
      // standard kinematically nonlinear analysis
      defgrd.MultiplyTT(xcurr, N_XYZ);
    }
    else
    {
      // in kinematically linear analysis the deformation gradient is equal to identity
      // no difference between reference and current state
      for (int i = 0; i < 3; ++i) defgrd(i, i) = 1.0;
    }

    if (IsParamsInterface())
    {
      double det_defgrd = defgrd.Determinant();
      if (det_defgrd < 0.0)
      {
        if (str_params_interface().IsTolerateErrors())
        {
          str_params_interface().SetEleEvalErrorFlag(
              STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
          stiffmatrix->Clear();
          force->Clear();
          return;
        }
        else
          FOUR_C_THROW("negative deformation gradient determinant");
      }  // if (det_defgrd<0.0)
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
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      bop(0, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.values(), true);
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
    {
      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);
    }
    else
    {
      // build the linearised strain epsilon = B_L . d
      glstrain.Multiply(bop, nodaldisp);
    }

    // deformation gradient consistent with (potentially EAS-modified) GL strains
    // without eas this is equal to the regular defgrd.
    CORE::LINALG::Matrix<3, 3> defgrd_mod(defgrd);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      M.shape(MAT::NUM_STRESS_3D, neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      // CORE::LINALG::SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch (eastype_)
      {
        case DRT::ELEMENTS::SoHex8::soh8_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easfull>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
              1.0, glstrain.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easmild>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
              1.0, glstrain.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_eassosh8>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
              1.0, glstrain.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // calculate deformation gradient consistent with modified GL strain tensor
      if (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->needs_defgrd())
        calc_consistent_defgrd(defgrd, glstrain, defgrd_mod);

      const double det_defgrd_mod = defgrd_mod.Determinant();
      if (det_defgrd_mod <= 0.0)
      {
        soh8_error_handling(det_defgrd_mod, params, __LINE__,
            STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
        return;
      }
    }  // ------------------------------------------------------------------ EAS

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (eastype_ != soh8_easnone)
        {
          FOUR_C_THROW(
              "EA strains are computed with the 'normal' deformation gradient from GL strains, and "
              "not with the deformation gradient that is consistent with EAS!\n"
              "Use the new solid elements instead!");
        }

        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
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
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
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
      case INPAR::STR::strain_log:
      {
        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");

        /// the Eularian logarithmic strain is defined as the natural logarithm of the left stretch
        /// tensor [1,2]: \f[
        ///    e_{log} = e_{hencky} = ln (\mathbf{V}) = \sum_{i=1}^3 (ln \lambda_i) \mathbf{n}_i
        ///    \otimes \mathbf{n}_i
        /// \f]
        ///< h3>References</h3>
        /// <ul>
        /// <li> [1] H. Xiao, Beijing, China, O. T. Bruhns and A. Meyers (1997) Logarithmic strain,
        /// logarithmic spin and logarithmic rate, Eq. 5 <li> [2] Caminero et al. (2011) Modeling
        /// large strain anisotropic elasto-plasticity with logarithmic strain and stress measures,
        /// Eq. 70
        /// </ul>
        ///
        /// \author HdV
        /// \date 08/13

        // eigenvalue decomposition (from elasthyper.cpp)
        CORE::LINALG::Matrix<3, 3> prstr2(true);  // squared principal stretches
        CORE::LINALG::Matrix<3, 1> prstr(true);   // principal stretch
        CORE::LINALG::Matrix<3, 3> prdir(true);   // principal directions
        CORE::LINALG::SYEV(cauchygreen, prstr2, prdir);

        // THE principal stretches
        for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

        // populating the logarithmic strain matrix
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> lnv(true);

        // checking if cauchy green is correctly determined to ensure eigen vectors in correct
        // direction i.e. a flipped eigenvector is also a valid solution C = \sum_{i=1}^3
        // (\lambda_i^2) \mathbf{n}_i \otimes \mathbf{n}_i
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> tempCG(true);

        for (int k = 0; k < 3; ++k)
        {
          double n_00, n_01, n_02, n_11, n_12, n_22 = 0.0;

          n_00 = prdir(0, k) * prdir(0, k);
          n_01 = prdir(0, k) * prdir(1, k);
          n_02 = prdir(0, k) * prdir(2, k);
          n_11 = prdir(1, k) * prdir(1, k);
          n_12 = prdir(1, k) * prdir(2, k);
          n_22 = prdir(2, k) * prdir(2, k);

          // only compute the symmetric components from a single eigenvector,
          // because eigenvalue directions are not consistent (it can be flipped)
          tempCG(0, 0) += (prstr(k)) * (prstr(k)) * n_00;
          tempCG(0, 1) += (prstr(k)) * (prstr(k)) * n_01;
          tempCG(0, 2) += (prstr(k)) * (prstr(k)) * n_02;
          tempCG(1, 0) += (prstr(k)) * (prstr(k)) * n_01;  // symmetry
          tempCG(1, 1) += (prstr(k)) * (prstr(k)) * n_11;
          tempCG(1, 2) += (prstr(k)) * (prstr(k)) * n_12;
          tempCG(2, 0) += (prstr(k)) * (prstr(k)) * n_02;  // symmetry
          tempCG(2, 1) += (prstr(k)) * (prstr(k)) * n_12;  // symmetry
          tempCG(2, 2) += (prstr(k)) * (prstr(k)) * n_22;

          // Computation of the Logarithmic strain tensor

          lnv(0, 0) += (std::log(prstr(k))) * n_00;
          lnv(0, 1) += (std::log(prstr(k))) * n_01;
          lnv(0, 2) += (std::log(prstr(k))) * n_02;
          lnv(1, 0) += (std::log(prstr(k))) * n_01;  // symmetry
          lnv(1, 1) += (std::log(prstr(k))) * n_11;
          lnv(1, 2) += (std::log(prstr(k))) * n_12;
          lnv(2, 0) += (std::log(prstr(k))) * n_02;  // symmetry
          lnv(2, 1) += (std::log(prstr(k))) * n_12;  // symmetry
          lnv(2, 2) += (std::log(prstr(k))) * n_22;
        }

        // compare CG computed with deformation gradient with CG computed
        // with eigenvalues and -vectors to determine/ensure the correct
        // orientation of the eigen vectors
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> diffCG(true);

        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            diffCG(i, j) = cauchygreen(i, j) - tempCG(i, j);
            // the solution to this problem is to evaluate the cauchygreen tensor with tempCG
            // computed with every combination of eigenvector orientations -- up to nine comparisons
            if (diffCG(i, j) > 1e-10)
              FOUR_C_THROW(
                  "eigenvector orientation error with the diffCG giving problems: %10.5e \n BUILD "
                  "SOLUTION TO FIX IT",
                  diffCG(i, j));
          }
        }

        (*elestrain)(gp, 0) = lnv(0, 0);
        (*elestrain)(gp, 1) = lnv(1, 1);
        (*elestrain)(gp, 2) = lnv(2, 2);
        (*elestrain)(gp, 3) = lnv(0, 1);
        (*elestrain)(gp, 4) = lnv(1, 2);
        (*elestrain)(gp, 5) = lnv(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        FOUR_C_THROW("requested strain type not available");
        break;
    }


    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix must be retrieved,
    ** all necessary data must be passed.
    */
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    UTILS::get_temperature_for_structural_material<CORE::FE::CellType::hex8>(shapefcts[gp], params);

    if (Material()->MaterialType() == CORE::Materials::m_constraintmixture ||
        Material()->MaterialType() == CORE::Materials::m_growthremodel_elasthyper ||
        Material()->MaterialType() == CORE::Materials::m_mixture)
    {
      CORE::LINALG::Matrix<NUMDIM_SOH8, 1> point(true);
      soh8_gauss_point_refe_coords(point, xrefe, gp);
      params.set("gp_coords_ref", point);

      // center of element in reference configuration
      point.Clear();
      soh8_element_center_refe_coords(point, xrefe);
      params.set("elecenter_coords_ref", point);
    }

    // if output is requested only active stresses are written.
    params.set<int>("iostress", iostress);

    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());
    so3mat->Evaluate(&defgrd_mod, &glstrain, params, &stress, &cmat, gp, Id());

    // stop if the material evaluation fails
    if (IsParamsInterface() and str_params_interface().IsTolerateErrors())
      if (str_params_interface().GetEleEvalErrorFlag() != STR::ELEMENTS::ele_error_none) return;

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp plastic strains (only in case of plastic strain output)
    switch (ioplstrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (eleplstrain == nullptr) FOUR_C_THROW("plastic strain data not available");
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plglstrain =
            params.get<CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain");
        for (int i = 0; i < 3; ++i) (*eleplstrain)(gp, i) = plglstrain(i);
        for (int i = 3; i < 6; ++i) (*eleplstrain)(gp, i) = 0.5 * plglstrain(i);
        break;
      }
      case INPAR::STR::strain_ea:
      {
        if (eastype_ != soh8_easnone)
        {
          FOUR_C_THROW(
              "EA strains are computed with the 'normal' deformation gradient from GL strains, and "
              "not with the deformation gradient that is consistent with EAS!\n"
              "Use the new solid elements instead!");
        }

        if (eleplstrain == nullptr) FOUR_C_THROW("plastic strain data not available");
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plglstrain =
            params.get<CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain");
        // rewriting Green-Lagrange strains in matrix format
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
        gl(0, 0) = plglstrain(0);
        gl(0, 1) = 0.5 * plglstrain(3);
        gl(0, 2) = 0.5 * plglstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = plglstrain(1);
        gl(1, 2) = 0.5 * plglstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = plglstrain(2);

        // inverse of deformation gradient
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
        temp.Multiply(gl, invdefgrd);
        euler_almansi.MultiplyTN(invdefgrd, temp);

        (*eleplstrain)(gp, 0) = euler_almansi(0, 0);
        (*eleplstrain)(gp, 1) = euler_almansi(1, 1);
        (*eleplstrain)(gp, 2) = euler_almansi(2, 2);
        (*eleplstrain)(gp, 3) = euler_almansi(0, 1);
        (*eleplstrain)(gp, 4) = euler_almansi(1, 2);
        (*eleplstrain)(gp, 5) = euler_almansi(0, 2);
        break;
      }
      case INPAR::STR::strain_none:
        break;
      default:
        FOUR_C_THROW("requested plastic strain type not available");
        break;
    }

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (eastype_ != soh8_easnone)
        {
          FOUR_C_THROW(
              "Cauchy stresses are computed with the 'normal' deformation gradient from 2PK "
              "stresses and not with the deformation gradient that is consistent with EAS!\n"
              "Use the new solid elements instead!");
        }

        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchystress(false);
        p_k2to_cauchy(&stress, &defgrd, &cauchystress);

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
        FOUR_C_THROW("requested stress type not available");
        break;
    }

    const double detJ_w = detJ * gpweights[gp];
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // structural force vector
    if (split_res && force_str != nullptr) force_str->MultiplyTN(detJ_w, bop, stress, 1.);

    // update stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      CORE::LINALG::Matrix<6, NUMDOF_SOH8> cb;
      cb.Multiply(cmat, bop);

      if (analyticalmaterialtangent_)
        stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);  // standard hex8 evaluation
      else
      {
        evaluate_finite_difference_material_tangent(stiffmatrix, stress, disp, detJ_w, detJ, detJ0,
            elediagonallength, bop, cb, N_XYZ, T0invT, M_GP, alpha, M, gp, params);
      }


      if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
      {
        // integrate `geometric' stiffness matrix and add to keu *****************
        CORE::LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
        sfac.Scale(detJ_w);            // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
        std::vector<double> SmB_L(3);  // intermediate Sm.B_L
        // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
        for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
        {
          SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
          SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
          SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
          for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
          {
            double bopstrbop = 0.0;  // intermediate value
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
              bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
            (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
            (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
            (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
          }

        }  // end of integrate `geometric' stiffness******************************
      }

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone)
      {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        CORE::LINALG::SerialDenseMatrix cM(MAT::NUM_STRESS_3D, neas_);  // temporary c . M
        switch (eastype_)
        {
          case DRT::ELEMENTS::SoHex8::soh8_easfull:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_easfull>(cM.values(), cmat.A(), M.values());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D,
                soh8_easfull>(1.0, Kaa, detJ_w, M, cM);
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.values(), detJ_w, M.values(), cb.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.values(), detJ_w, M.values(), stress.A());
            break;
          case DRT::ELEMENTS::SoHex8::soh8_easmild:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_easmild>(cM.values(), cmat.A(), M.values());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D,
                soh8_easmild>(1.0, Kaa, detJ_w, M, cM);
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.values(), detJ_w, M.values(), cb.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.values(), detJ_w, M.values(), stress.A());
            break;
          case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_eassosh8>(cM.values(), cmat.A(), M.values());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
                soh8_eassosh8>(1.0, Kaa, detJ_w, M, cM);
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.values(), detJ_w, M.values(), cb.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.values(), detJ_w, M.values(), stress.A());
            break;
          case DRT::ELEMENTS::SoHex8::soh8_easnone:
            break;
          default:
            FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }  // ---------------------------------------------------------------- EAS
    }

    if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      const double density = Material()->Density(gp);

      // integrate consistent mass matri
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += massfactor;
        }
      }

      // check for non constant mass matrix
      if (so3mat->VaryingDensity())
      {
        /*
         If the density, i.e. the mass matrix, is not constant, a linearization is neccessary.
         In general, the mass matrix can be dependent on the displacements, the velocities and the
         accelerations. We write all the additional terms into the mass matrix, hence, conversion
         from accelerations to velocities and displacements are needed. As those conversions depend
         on the time integration scheme, the factors are set within the respective time integrators
         and read from the parameter list inside the element (this is a little ugly...). */
        double timintfac_dis = 0.0;
        double timintfac_vel = 0.0;
        if (IsParamsInterface())
        {
          timintfac_dis = str_params_interface().GetTimIntFactorDisp();
          timintfac_vel = str_params_interface().GetTimIntFactorVel();
        }
        else
        {
          timintfac_dis = params.get<double>("timintfac_dis");
          timintfac_vel = params.get<double>("timintfac_vel");
        }
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_disp(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_vel(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass(true);

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        so3mat->EvaluateNonLinMass(
            &defgrd, &glstrain, params, &linmass_disp, &linmass_vel, gp, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        CORE::LINALG::Matrix<NUMDIM_SOH8, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
            myacc(idim) += shapefcts[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOH8)];

        if (stiffmatrix != nullptr)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          CORE::LINALG::Matrix<1, NUMDOF_SOH8> cb;
          cb.MultiplyTN(linmass_disp, bop);
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
          {
            double factor = detJ_w * shapefcts[gp](inod);
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
            {
              double massfactor = factor * myacc(idim);
              for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
                for (int jdim = 0; jdim < NUMDIM_SOH8; ++jdim)
                  (*massmatrix)(inod * NUMDIM_SOH8 + idim, jnod * NUMDIM_SOH8 + jdim) +=
                      massfactor * cb(jnod * NUMDIM_SOH8 + jdim);
            }
          }
        }

        // internal force vector without EAS terms
        if (forceinert != nullptr)
        {
          // integrate nonlinear inertia force term
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
          {
            double forcefactor = shapefcts[gp](inod) * detJ_w;
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
              (*forceinert)(inod * NUMDIM_SOH8 + idim) += forcefactor * density * myacc(idim);
          }
        }
      }

    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  // rhs norm of eas equations
  if (eastype_ != soh8_easnone && split_res && force != nullptr)
    // only add for row-map elements
    if (params.get<int>("MyPID") == Owner())
      params.get<double>("cond_rhs_norm") += pow(CORE::LINALG::Norm2(feas), 2.);

  if (force != nullptr && stiffmatrix != nullptr)
  {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ != soh8_easnone)
    {
      // we need the inverse of Kaa. Catch Inf/NaN case
      const double norm1 = Kaa.normOne();
      if (std::isnan(norm1) || std::isinf(norm1) || norm1 == 0.)
      {
        for (int i = 0; i < Kaa.numCols(); ++i)
          for (int j = 0; j < Kaa.numRows(); ++j)
            Kaa(j, i) = std::numeric_limits<double>::quiet_NaN();
      }
      else
      {
        using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
        using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
        Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_inverseKaa;
        solve_for_inverseKaa.setMatrix(Teuchos::rcpFromRef(Kaa));
        solve_for_inverseKaa.invert();
      }

      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas

      CORE::LINALG::SerialDenseMatrix KdaKaa(NUMDOF_SOH8, neas_);  // temporary Kda.Kaa^{-1}
      switch (eastype_)
      {
        case DRT::ELEMENTS::SoHex8::soh8_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_easfull, soh8_easfull>(
              KdaKaa, Kda, Kaa);
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easfull, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.values(), Kda.values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easfull, 1>(
              1.0, force->A(), -1.0, KdaKaa.values(), feas.values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_easmild, soh8_easmild>(
              KdaKaa, Kda, Kaa);
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easmild, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.values(), Kda.values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easmild, 1>(
              1.0, force->A(), -1.0, KdaKaa.values(), feas.values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_eassosh8,
              soh8_eassosh8>(KdaKaa, Kda, Kaa);
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.values(), Kda.values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, 1>(
              1.0, force->A(), -1.0, KdaKaa.values(), feas.values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // store current EAS data in history
      for (int i = 0; i < neas_; ++i)
      {
        for (int j = 0; j < neas_; ++j) (*oldKaainv)(i, j) = Kaa(i, j);
        for (int j = 0; j < NUMDOF_SOH8; ++j) (*oldKda)(i, j) = Kda(i, j);
        (*oldfeas)(i, 0) = feas(i);
      }
    }  // -------------------------------------------------------------------- EAS
  }
  return;
}  // DRT::ELEMENTS::So_hex8::nlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_lumpmass(CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).numRows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Shape fcts at all 8 Gauss Points             maf 05/08|
 *----------------------------------------------------------------------*/
std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> DRT::ELEMENTS::SoHex8::soh8_shapefcts() const
{
  std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts(NUMGPT_SOH8);

  // fill up nodal f at each gp
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    const CORE::LINALG::Matrix<NUMDIM_SOH8, 1> rst_gp(gp_rule_.Point(gp), true);
    CORE::FE::shape_function<CORE::FE::CellType::hex8>(rst_gp, shapefcts[gp]);
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Shape fct derivs at all 8 Gauss Points       maf 05/08|
 *----------------------------------------------------------------------*/
std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> DRT::ELEMENTS::SoHex8::soh8_derivs()
    const
{
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs(NUMGPT_SOH8);

  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    soh8_derivs(derivs[gp], gp);
  }

  return derivs;
}

// Evaluate the derivatives of the shape functions for a specific Gauss point
void DRT::ELEMENTS::SoHex8::soh8_derivs(
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& derivs, const int gp) const
{
  const CORE::LINALG::Matrix<NUMDIM_SOH8, 1> rst_gp(gp_rule_.Point(gp), true);
  CORE::FE::shape_function_deriv1<CORE::FE::CellType::hex8>(rst_gp, derivs);
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Weights at all 8 Gauss Points                maf 05/08|
 *----------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::SoHex8::soh8_weights() const
{
  std::vector<double> weights(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp) weights[gp] = gp_rule_.Weight(gp);

  return weights;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_create_eas_backup_state(const std::vector<double>& displ_incr)
{
  if (eastype_ == soh8_easnone) return;

  // --- create EAS state backup ----------------------------------------------
  {
    const auto* alpha = &easdata_.alpha;
    if (not alpha) FOUR_C_THROW("Can't access the current enhanced strain state.");

    auto* alpha_backup_ptr = &easdata_.alpha_backup;
    if (alpha_backup_ptr)
      *alpha_backup_ptr = *alpha;
    else
      easdata_.alpha_backup = *alpha;
  }

  // --- create EAS increment backup ------------------------------------------
  {
    // compute the current eas increment
    CORE::LINALG::SerialDenseMatrix eas_inc(neas_, 1);
    soh8_compute_eas_inc(displ_incr, &eas_inc);

    auto* eas_inc_backup_ptr = &easdata_.eas_inc_backup;
    if (eas_inc_backup_ptr)
      *eas_inc_backup_ptr = eas_inc;
    else
      easdata_.eas_inc_backup = eas_inc;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::soh8_recover_from_eas_backup_state()
{
  if (eastype_ == soh8_easnone) return;

  CORE::LINALG::SerialDenseMatrix* alpha = nullptr;
  CORE::LINALG::SerialDenseMatrix* eas_inc = nullptr;

  // --- recover state from EAS backup ----------------------------------------
  {
    const auto* alpha_backup = &easdata_.alpha_backup;
    if (not alpha_backup)
      FOUR_C_THROW(
          "Can't access the enhanced strain backup state. Did you "
          "create a backup? See soh8_create_eas_backup_state().");

    alpha = &easdata_.alpha;
    if (not alpha) FOUR_C_THROW("Can't access the enhanced strain state.");

    *alpha = *alpha_backup;
  }

  // --- recover increment from EAS backup ------------------------------------
  {
    const auto* eas_inc_backup = &easdata_.eas_inc_backup;
    if (not eas_inc_backup)
      FOUR_C_THROW(
          "Can't access the enhanced strain increment backup. Did you "
          "create a backup? See soh8_create_eas_backup_state().");

    eas_inc = &easdata_.eas_inc;
    if (not eas_inc) FOUR_C_THROW("Can't access the enhanced strain increment.");

    *eas_inc = *eas_inc_backup;
  }

  // Finally, we have to update the backup state, otherwise a follow-up
  // step length adaption will lead to a wrong eas state.
  {
    old_step_length_ = 0.0;
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 04/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::SoHex8*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8* failed");
    actele->init_jacobian_mapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::def_gradient(const std::vector<double>& disp,
    CORE::LINALG::SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;  // current  coord. of element
  UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(disp, xdisp);

  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    CORE::LINALG::Matrix<3, 3> invJdef;
    prestress.StoragetoMatrix(gp, invJdef, prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
    N_xyz.Multiply(invJdef, derivs[gp]);

    // build defgrd (independent of xrefe!)
    CORE::LINALG::Matrix<3, 3> defgrd(true);
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
    {
      defgrd.MultiplyTT(xdisp, N_xyz);
    }
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;

    prestress.MatrixtoStorage(gp, defgrd, gpdefgrd);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::update_jacobian_mapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();

  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp(false);
  UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(disp, xdisp);

  CORE::LINALG::Matrix<3, 3> invJhist;
  CORE::LINALG::Matrix<3, 3> invJ;
  CORE::LINALG::Matrix<3, 3> defgrd(true);
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
  CORE::LINALG::Matrix<3, 3> invJnew;
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp, invJhist, prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist, derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
    {
      defgrd.MultiplyTT(xdisp, N_xyz);
    }
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd, invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp, invJnew, prestress.JHistory());
  }  // for (int gp=0; gp<NUMGPT_SOH8; ++gp)

  return;
}

/*---------------------------------------------------------------------------------------------*
 |  Update history variables (e.g. remodeling of fiber directions) (protected)      braeu 07/16|
 *---------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::update_element(std::vector<double>& disp,
    Teuchos::ParameterList& params, const Teuchos::RCP<CORE::MAT::Material>& mat)
{
  // Calculate current deformation gradient
  if ((mat->MaterialType() == CORE::Materials::m_constraintmixture) ||
      (mat->MaterialType() == CORE::Materials::m_elasthyper) ||
      (mat->MaterialType() == CORE::Materials::m_growthremodel_elasthyper) ||
      (SolidMaterial()->UsesExtendedUpdate()))
  {
    CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(false);
    CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp(false);
    CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr(false);

    UTILS::EvaluateNodalCoordinates<CORE::FE::CellType::hex8, 3>(Nodes(), xrefe);
    UTILS::EvaluateNodalDisplacements<CORE::FE::CellType::hex8, 3>(disp, xdisp);
    UTILS::EvaluateCurrentNodalCoordinates<CORE::FE::CellType::hex8, 3>(xrefe, xdisp, xcurr);

    /* =========================================================================*/
    /* ================================================= Loop over Gauss Points */
    /* =========================================================================*/
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
    // interpolated values of stress and defgrd for remodeling
    CORE::LINALG::Matrix<3, 3> avg_stress(true);
    CORE::LINALG::Matrix<3, 3> avg_defgrd(true);

    // build deformation gradient wrt to material configuration
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(false);
    params.set<int>("numgp", static_cast<int>(NUMGPT_SOH8));

    // center of element in reference configuration
    CORE::LINALG::Matrix<NUMDIM_SOH8, 1> point(false);
    point.Clear();
    soh8_element_center_refe_coords(point, xrefe);
    params.set("elecenter_coords_ref", point);

    for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
    {
      soh8_gauss_point_refe_coords(point, xrefe, gp);
      params.set("gp_coords_ref", point);
      CORE::LINALG::Matrix<3, 8> derivs(false);
      soh8_derivs(derivs, gp);

      // Compute deformation gradient
      UTILS::compute_deformation_gradient<CORE::FE::CellType::hex8>(
          defgrd, kintype_, xdisp, xcurr, invJ_[gp], derivs, pstype_, prestress_, gp);

      // call material update if material = m_growthremodel_elasthyper (calculate and update
      // inelastic deformation gradient)
      if (SolidMaterial()->UsesExtendedUpdate())
      {
        SolidMaterial()->Update(defgrd, gp, params, Id());
      }
    }  // end loop over gauss points
  }

  // store EAS parameters
  if (eastype_ != soh8_easnone)
  {
    soh8_easupdate();

    // reset EAS internal force
    CORE::LINALG::SerialDenseMatrix* oldfeas = &easdata_.feas;
    oldfeas->putScalar(0.0);
  }
  SolidMaterial()->Update();

  return;
}


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::g_lto_ea(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* euler_almansi)
{
  // e = F^{T-1} . E . F^{-1}

  // rewrite Green-Lagrange strain in tensor notation
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
  gl(0, 0) = (*glstrain)(0);
  gl(0, 1) = 0.5 * (*glstrain)(3);
  gl(0, 2) = 0.5 * (*glstrain)(5);
  gl(1, 0) = gl(0, 1);
  gl(1, 1) = (*glstrain)(1);
  gl(1, 2) = 0.5 * (*glstrain)(4);
  gl(2, 0) = gl(0, 2);
  gl(2, 1) = gl(1, 2);
  gl(2, 2) = (*glstrain)(2);

  // inverse of deformation gradient
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
  invdefgrd.Invert((*defgrd));

  // (3x3) = (3x3) (3x3) (3x3)
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
  temp.Multiply(gl, invdefgrd);
  (*euler_almansi).MultiplyTN(invdefgrd, temp);

}  // GLtoEA()


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::p_k2to_cauchy(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* cauchystress)
{
  // calculate the Jacobi-deterinant
  const double detF = (*defgrd).Determinant();

  // sigma = 1/J . F . S . F^T
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> pkstress;
  pkstress(0, 0) = (*stress)(0);
  pkstress(0, 1) = (*stress)(3);
  pkstress(0, 2) = (*stress)(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = (*stress)(1);
  pkstress(1, 2) = (*stress)(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = (*stress)(2);

  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
  temp.Multiply((1.0 / detF), (*defgrd), pkstress);
  (*cauchystress).MultiplyNT(temp, (*defgrd));

}  // PK2toCauchy()

/*----------------------------------------------------------------------*
 |  Calculate consistent deformation gradient               seitz 04/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::calc_consistent_defgrd(const CORE::LINALG::Matrix<3, 3>& defgrd_disp,
    CORE::LINALG::Matrix<6, 1> glstrain_mod, CORE::LINALG::Matrix<3, 3>& defgrd_mod) const
{
  CORE::LINALG::Matrix<3, 3> R;       // rotation tensor
  CORE::LINALG::Matrix<3, 3> U_mod;   // modified right stretch tensor
  CORE::LINALG::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
  CORE::LINALG::Matrix<3, 3> EW;      // temporarily store eigenvalues
  CORE::LINALG::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
  CORE::LINALG::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

  // ******************************************************************
  // calculate modified right stretch tensor
  // ******************************************************************
  for (int i = 0; i < 3; i++) U_mod(i, i) = 2. * glstrain_mod(i) + 1.;
  U_mod(0, 1) = glstrain_mod(3);
  U_mod(1, 0) = glstrain_mod(3);
  U_mod(1, 2) = glstrain_mod(4);
  U_mod(2, 1) = glstrain_mod(4);
  U_mod(0, 2) = glstrain_mod(5);
  U_mod(2, 0) = glstrain_mod(5);

  CORE::LINALG::SYEV(U_mod, EW, U_mod);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_mod, EW);
  tmp2.MultiplyNT(tmp, U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

  CORE::LINALG::SYEV(U_disp, EW, U_disp);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_disp, EW);
  tmp2.MultiplyNT(tmp, U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(defgrd_disp, U_disp);
  defgrd_mod.Multiply(R, U_mod);

  // you're done here
  return;
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | check the constitutive tensor and/or use the approximation as        |
 | elastic stiffness matrix                                  rauch 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::evaluate_finite_difference_material_tangent(
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress, std::vector<double>& disp,
    const double detJ_w, const double detJ, const double detJ0, const double charelelength,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8>& bop,
    const CORE::LINALG::Matrix<6, NUMDOF_SOH8>& cb,
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& N_XYZ,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& T0invT,
    const std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP,
    const CORE::LINALG::SerialDenseMatrix* alpha, CORE::LINALG::SerialDenseMatrix& M, const int gp,
    Teuchos::ParameterList& params)
{
  // build elastic stiffness matrix directly by finite differences

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());

#ifdef MATERIALFDCHECK
  static CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> stiffmatrix_analytical;
  static CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> stiffmatrix_fd;
  if (gp == 0)
  {
    stiffmatrix_analytical.PutScalar(0.0);
    stiffmatrix_fd.PutScalar(0.0);
  }
  stiffmatrix_analytical.MultiplyTN(detJ_w, bop, cb, 1.0);
#endif

  const double delta = charelelength * 1.0e-08;

  // matrices and vectors
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> cb_fd(true);
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress_fd(true);
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> finitedifference(true);

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current   coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;

  // get nodes
  CORE::Nodes::Node** nodes = Nodes();

  //////////////////////////////////////////////////////////////////////////////
  ////// evaluate partial derivatives of stress (S(d_n+delta) - S(d_n))/delta
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////
  // loop over columns
  /////////////////////////////////////////////
  for (int i = 0; i < NUMDOF_SOH8; ++i)
  {
    // undo disturbance for disp[i-1]
    if (i > 0) disp[i - 1] -= delta;

    // disturb displacements
    disp[i] += delta;

    for (int k = 0; k < NUMNOD_SOH8; ++k)
    {
      const auto& x = nodes[k]->X();
      xrefe(k, 0) = x[0];
      xrefe(k, 1) = x[1];
      xrefe(k, 2) = x[2];

      xcurr(k, 0) = xrefe(k, 0) + disp[k * NODDOF_SOH8 + 0];
      xcurr(k, 1) = xrefe(k, 1) + disp[k * NODDOF_SOH8 + 1];
      xcurr(k, 2) = xrefe(k, 2) + disp[k * NODDOF_SOH8 + 2];
    }

    // build deformation gradient wrt to material configuration

    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd_fd(true);
    defgrd_fd.MultiplyTT(xcurr, N_XYZ);


    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen_fd(true);
    cauchygreen_fd.MultiplyTN(defgrd_fd, defgrd_fd);

    // Green-Lagrange strains matrix E = 0.5 * (cauchygreen_fd - Identity)
    // GL strain vector glstrain_fd={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain_fd(true);
    glstrain_fd(0) = 0.5 * (cauchygreen_fd(0, 0) - 1.0);
    glstrain_fd(1) = 0.5 * (cauchygreen_fd(1, 1) - 1.0);
    glstrain_fd(2) = 0.5 * (cauchygreen_fd(2, 2) - 1.0);
    glstrain_fd(3) = cauchygreen_fd(0, 1);
    glstrain_fd(4) = cauchygreen_fd(1, 2);
    glstrain_fd(5) = cauchygreen_fd(2, 0);

    // deformation gradient consistent with (potentially EAS-modified) GL strains
    // without eas this is equal to the regular defgrd_fd.
    CORE::LINALG::Matrix<3, 3> defgrd_fd_mod(defgrd_fd);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      FOUR_C_THROW("be careful ! fdcheck has not been tested with EAS, yet! ");
      M.shape(MAT::NUM_STRESS_3D, neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      // CORE::LINALG::SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch (eastype_)
      {
        case DRT::ELEMENTS::SoHex8::soh8_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easfull>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
              1.0, glstrain_fd.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easmild>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
              1.0, glstrain_fd.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_eassosh8>(M.values(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).values());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
              1.0, glstrain_fd.A(), 1.0, M.values(), alpha->values());
          break;
        case DRT::ELEMENTS::SoHex8::soh8_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // calculate deformation gradient consistent with modified GL strain tensor
      if (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->needs_defgrd())
        calc_consistent_defgrd(defgrd_fd, glstrain_fd, defgrd_fd_mod);
    }  // ------------------------------------------------------------------ EAS

    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat_fd;
    so3mat->Evaluate(&defgrd_fd_mod, &glstrain_fd, params, &stress_fd, &cmat_fd, gp, Id());

    // finite difference approximation of partial derivative
    //
    //      d S_ij,columnindex
    //
    //
    finitedifference.Update(1.0, stress_fd, -1.0, stress);
    finitedifference.Scale(1.0 / delta);

    /////////////////////////
    // loop over rows
    ////////////////////////
    for (int j = 0; j < MAT::NUM_STRESS_3D; ++j)
    {
      cb_fd(j, i) = finitedifference(j, 0);
    }  // j-loop (rows)


    // reset disp at last loop execution
    if (i == (NUMDOF_SOH8 - 1))
    {
      disp[i] -= delta;

      // reset xcurr
      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        const auto& x = nodes[k]->X();
        xrefe(k, 0) = x[0];
        xrefe(k, 1) = x[1];
        xrefe(k, 2) = x[2];

        xcurr(k, 0) = xrefe(k, 0) + disp[k * NODDOF_SOH8 + 0];
        xcurr(k, 1) = xrefe(k, 1) + disp[k * NODDOF_SOH8 + 1];
        xcurr(k, 2) = xrefe(k, 2) + disp[k * NODDOF_SOH8 + 2];
      }
    }

  }  // i-loop (columns)
  /////////////////////////////// FD LOOP


  ///////////////////////////////////////
  // build approximated stiffness matrix
  ///////////////////////////////////////
  stiffmatrix->MultiplyTN(detJ_w, bop, cb_fd, 1.0);
#ifdef MATERIALFDCHECK
  stiffmatrix_fd.MultiplyTN(detJ_w, bop, cb_fd, 1.0);
#endif

  ///////////////////////////////////////
#ifdef MATERIALFDCHECK
  // after last gp was evaluated
  if (gp == (NUMGPT_SOH8 - 1))
  {
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> errormatrix(true);

    // calc error (subtraction stiffmatrix - stiffmatrix_analytical)
    errormatrix.Update(1.0, stiffmatrix_fd, -1.0, stiffmatrix_analytical);

    for (int i = 0; i < NUMDOF_SOH8; ++i)
    {
      for (int j = 0; j < NUMDOF_SOH8; ++j)
      {
        double relerror = abs(errormatrix(i, j)) / abs((stiffmatrix_analytical)(i, j));
        if (std::min(abs(errormatrix(i, j)), relerror) > delta * 1000.0)
        {
          std::cout << "ELEGID:" << this->Id() << "  gp: " << gp << "  ROW: " << i << "  COL: " << j
                    << "    REL. ERROR: " << relerror
                    << "    ABS. ERROR: " << abs(errormatrix(i, j))
                    << "    stiff. val: " << stiffmatrix_analytical(i, j)
                    << "    approx. val: " << stiffmatrix_fd(i, j) << std::endl;
        }
      }
    }  // check errors
  }    // if last gp of element is reached
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8::get_cauchy_n_dir_and_derivatives_at_xi(
    const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir,
    double& cauchy_n_dir, CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd2,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dn,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dT,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
    double* d_cauchyndir_dc)
{
  FOUR_C_THROW_UNLESS(eastype_ == soh8_easnone && !PRESTRESS::IsMulf(),
      "Evaluation of the Cauchy stress is not possible for EAS-elements or MULF prestressing.");
  if (temp || d_cauchyndir_dT || d2_cauchyndir_dd_dT)
    FOUR_C_THROW("Thermo-elastic Nitsche contact not yet implemented in so hex8");

  cauchy_n_dir = 0.0;

  static CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(true);  // reference coord. of element
  static CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr(true);  // current  coord. of element
  xrefe.Clear();
  xcurr.Clear();
  CORE::Nodes::Node** nodes = Nodes();

  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int d = 0; d < NUMDIM_SOH8; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * NODDOF_SOH8 + d];
    }
  }

  static CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> deriv(true);
  deriv.Clear();
  CORE::FE::shape_function_deriv1<CORE::FE::CellType::hex8>(xi, deriv);

  static CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ(true);
  static CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ(true);
  invJ.Multiply(1.0, deriv, xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, deriv, 0.0);
  static CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static CORE::LINALG::Matrix<9, NUMDOF_SOH8> d_F_dd(true);
  d_F_dd.Clear();
  if (d_cauchyndir_dd || d2_cauchyndir_dd_dn || d2_cauchyndir_dd_ddir || d2_cauchyndir_dd2 ||
      d2_cauchyndir_dd_dxi)
  {
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      d_F_dd(0, NODDOF_SOH8 * i + 0) = N_XYZ(0, i);
      d_F_dd(1, NODDOF_SOH8 * i + 1) = N_XYZ(1, i);
      d_F_dd(2, NODDOF_SOH8 * i + 2) = N_XYZ(2, i);
      d_F_dd(3, NODDOF_SOH8 * i + 0) = N_XYZ(1, i);
      d_F_dd(4, NODDOF_SOH8 * i + 1) = N_XYZ(2, i);
      d_F_dd(5, NODDOF_SOH8 * i + 0) = N_XYZ(2, i);
      d_F_dd(6, NODDOF_SOH8 * i + 1) = N_XYZ(0, i);
      d_F_dd(7, NODDOF_SOH8 * i + 2) = N_XYZ(1, i);
      d_F_dd(8, NODDOF_SOH8 * i + 2) = N_XYZ(0, i);
    }
  }

  static CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF(true);
  static CORE::LINALG::Matrix<9, 9> d2_cauchyndir_dF2(true);
  static CORE::LINALG::Matrix<9, NUMDIM_SOH8> d2_cauchyndir_dF_dn(true);
  static CORE::LINALG::Matrix<9, NUMDIM_SOH8> d2_cauchyndir_dF_ddir(true);

  SolidMaterial()->evaluate_cauchy_n_dir_and_derivatives(defgrd, n, dir, cauchy_n_dir,
      d_cauchyndir_dn, d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2,
      &d2_cauchyndir_dF_dn, &d2_cauchyndir_dF_ddir, -1, Id(), concentration, nullptr, nullptr,
      nullptr);

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->reshape(NUMDOF_SOH8, 1);
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1> d_cauchyndir_dd_mat(d_cauchyndir_dd->values(), true);
    d_cauchyndir_dd_mat.MultiplyTN(1.0, d_F_dd, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dn)
  {
    d2_cauchyndir_dd_dn->reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> d2_cauchyndir_dd_dn_mat(
        d2_cauchyndir_dd_dn->values(), true);
    d2_cauchyndir_dd_dn_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dn, 0.0);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> d2_cauchyndir_dd_ddir_mat(
        d2_cauchyndir_dd_ddir->values(), true);
    d2_cauchyndir_dd_ddir_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_ddir, 0.0);
  }

  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->reshape(NUMDOF_SOH8, NUMDOF_SOH8);
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> d2_cauchyndir_dd2_mat(
        d2_cauchyndir_dd2->values(), true);
    static CORE::LINALG::Matrix<9, NUMDOF_SOH8> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd2_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF2_d_F_dd, 0.0);
  }

  // prepare evaluation of d_cauchyndir_dxi or d2_cauchyndir_dd_dxi
  static CORE::LINALG::Matrix<9, NUMDIM_SOH8> d_F_dxi(true);
  static CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::hex8>::numderiv2,
      NUMNOD_SOH8>
      deriv2(true);
  d_F_dxi.Clear();
  deriv2.Clear();

  if (d_cauchyndir_dxi or d2_cauchyndir_dd_dxi)
  {
    CORE::FE::shape_function_deriv2<CORE::FE::CellType::hex8>(xi, deriv2);

    static CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xXF(true);
    static CORE::LINALG::Matrix<NUMDIM_SOH8,
        CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::hex8>::numderiv2>
        xXFsec(true);
    xXF.Update(1.0, xcurr, 0.0);
    xXF.MultiplyNT(-1.0, xrefe, defgrd, 1.0);
    xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

    for (int a = 0; a < NUMDIM_SOH8; ++a)
    {
      for (int b = 0; b < NUMDIM_SOH8; ++b)
      {
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 0) +=
            xXFsec(a, 0) * invJ(b, 0) + xXFsec(a, 3) * invJ(b, 1) + xXFsec(a, 4) * invJ(b, 2);
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 1) +=
            xXFsec(a, 3) * invJ(b, 0) + xXFsec(a, 1) * invJ(b, 1) + xXFsec(a, 5) * invJ(b, 2);
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 2) +=
            xXFsec(a, 4) * invJ(b, 0) + xXFsec(a, 5) * invJ(b, 1) + xXFsec(a, 2) * invJ(b, 2);
      }
    }
  }

  if (d_cauchyndir_dxi)
  {
    d_cauchyndir_dxi->MultiplyTN(1.0, d_F_dxi, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dxi)
  {
    d2_cauchyndir_dd_dxi->reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> d2_cauchyndir_dd_dxi_mat(
        d2_cauchyndir_dd_dxi->values(), true);

    static CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::hex8>::numderiv2,
        NUMNOD_SOH8>
        deriv2(true);
    deriv2.Clear();
    CORE::FE::shape_function_deriv2<CORE::FE::CellType::hex8>(xi, deriv2);

    static CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::hex8>::numderiv2,
        NUMDIM_SOH8>
        Xsec(true);
    static CORE::LINALG::Matrix<NUMNOD_SOH8,
        CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::hex8>::numderiv2>
        N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static CORE::LINALG::Matrix<9, NUMDOF_SOH8> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d2_cauchyndir_dF2_d_F_dd, d_F_dxi, 0.0);

    static CORE::LINALG::Matrix<9, NUMDIM_SOH8 * NUMDOF_SOH8> d2_F_dxi_dd(true);
    d2_F_dxi_dd.Clear();
    for (int i = 0; i < NUMDIM_SOH8; ++i)
    {
      for (int j = 0; j < NUMDIM_SOH8; ++j)
      {
        for (int k = 0; k < NUMNOD_SOH8; ++k)
        {
          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOH8 * (NODDOF_SOH8 * k + i) + 0) +=
              deriv2(0, k) * invJ(j, 0) + deriv2(3, k) * invJ(j, 1) + deriv2(4, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 0) * invJ(j, 0) - N_XYZ_Xsec(k, 3) * invJ(j, 1) -
              N_XYZ_Xsec(k, 4) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOH8 * (NODDOF_SOH8 * k + i) + 1) +=
              deriv2(3, k) * invJ(j, 0) + deriv2(1, k) * invJ(j, 1) + deriv2(5, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 3) * invJ(j, 0) - N_XYZ_Xsec(k, 1) * invJ(j, 1) -
              N_XYZ_Xsec(k, 5) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOH8 * (NODDOF_SOH8 * k + i) + 2) +=
              deriv2(4, k) * invJ(j, 0) + deriv2(5, k) * invJ(j, 1) + deriv2(2, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 4) * invJ(j, 0) - N_XYZ_Xsec(k, 5) * invJ(j, 1) -
              N_XYZ_Xsec(k, 2) * invJ(j, 2);

          for (int l = 0; l < NUMDIM_SOH8; ++l)
          {
            d2_cauchyndir_dd_dxi_mat(k * 3 + i, l) +=
                d_cauchyndir_dF(VoigtMapping::NonSymToVoigt9(i, j), 0) *
                d2_F_dxi_dd(
                    VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOH8 * (NODDOF_SOH8 * k + i) + l);
          }
        }
      }
    }
  }

  if (d_cauchyndir_dc != nullptr)
  {
    static CORE::LINALG::Matrix<9, 1> d_F_dc(true);
    SolidMaterial()->evaluate_linearization_od(defgrd, *concentration, &d_F_dc);
    *d_cauchyndir_dc = d_cauchyndir_dF.Dot(d_F_dc);
  }
}

FOUR_C_NAMESPACE_CLOSE
