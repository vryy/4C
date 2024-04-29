/*----------------------------------------------------------------------*/
/*! \file
\brief some element evaluate
\level 1


*/
/*----------------------------------------------------------------------*/
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_lib_condition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_aaaraghavanvorp_damage.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_mat_viscoanisotropic.hpp"
#include "4C_mat_viscoelasthyper.hpp"
#include "4C_mat_visconeohooke.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_sh8.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoSh8::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  // get parameter interface
  SetParamsInterfacePtr(params);

  CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec2(elevec2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOH8, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      FOUR_C_THROW("No action supplied");
    else if (action == "calc_struct_linstiff")
      act = ELEMENTS::struct_calc_linstiff;
    else if (action == "calc_struct_nlnstiff")
      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action == "calc_struct_internalforce")
      act = ELEMENTS::struct_calc_internalforce;
    else if (action == "calc_struct_linstiffmass")
      act = ELEMENTS::struct_calc_linstiffmass;
    else if (action == "calc_struct_nlnstiffmass")
      act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action == "calc_struct_nlnstifflmass")
      act = ELEMENTS::struct_calc_nlnstifflmass;
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
    else if (action == "calc_struct_eleload")
      act = ELEMENTS::struct_calc_eleload;
    else if (action == "calc_struct_fsiload")
      act = ELEMENTS::struct_calc_fsiload;
    else if (action == "calc_struct_update_istep")
      act = ELEMENTS::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = ELEMENTS::struct_calc_reset_istep;
    else if (action == "multi_eas_init")
      act = ELEMENTS::multi_init_eas;
    else if (action == "multi_eas_set")
      act = ELEMENTS::multi_set_eas;
    else if (action == "multi_calc_dens")
      act = ELEMENTS::multi_calc_dens;
    else if (action == "multi_readrestart")
      act = ELEMENTS::multi_readrestart;
    else if (action == "calc_stc_matrix")
      act = ELEMENTS::shell_calc_stc_matrix;
    else if (action == "calc_stc_matrix_inverse")
      act = ELEMENTS::shell_calc_stc_matrix_inverse;
    else if (action == "calc_struct_recover")
      act = ELEMENTS::struct_calc_recover;
    else if (action == "calc_struct_energy")
      act = ELEMENTS::struct_calc_energy;
    else if (action == "calc_struct_predict")
      return 0;
    else
      FOUR_C_THROW("Unknown type of action for So_Sh8: %s", action.c_str());
  }

  // what should the element do
  switch (act)
  {
    // linear stiffness
    case ELEMENTS::struct_calc_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (double& i : mydisp) i = 0.0;
      std::vector<double> myres(lm.size());
      for (double& myre : myres) myre = 0.0;
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (eastype_ != DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        sosh8_nlnstiffmass(lm, mydisp, myres, &elemat1, nullptr, &elevec1, nullptr, nullptr,
            nullptr, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
      }
      else
      {
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &elemat1, nullptr, &elevec1,
            nullptr, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
    }
    break;

    // nonlinear stiffness and internal force vector
    case ELEMENTS::struct_calc_nlnstiff:
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
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (eastype_ != DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        sosh8_nlnstiffmass(lm, mydisp, myres, &elemat1, nullptr, &elevec1, &elevec3, nullptr,
            nullptr, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
      }
      else
      {
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &elemat1, nullptr, &elevec1,
            nullptr, &elevec3, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
    }
    break;

    // internal force vector only
    case ELEMENTS::struct_calc_internalforce:
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
      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (eastype_ != DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        sosh8_nlnstiffmass(lm, mydisp, myres, &myemat, nullptr, &elevec1, nullptr, nullptr, nullptr,
            params, INPAR::STR::stress_none, INPAR::STR::strain_none);
      }
      else
      {
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &myemat, nullptr, &elevec1,
            nullptr, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
    }
    break;

    // linear stiffness and consistent mass matrix
    case ELEMENTS::struct_calc_linstiffmass:
      FOUR_C_THROW("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    // nonlinear stiffness, internal force vector, and consistent/lumped mass matrix
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
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

      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (eastype_ != DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        sosh8_nlnstiffmass(lm, mydisp, myres, &elemat1, &elemat2, &elevec1, &elevec3, nullptr,
            nullptr, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
      }
      else
      {
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &elemat1, &elemat2, &elevec1,
            nullptr, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);
      }
      // lump mass
      if (act == ELEMENTS::struct_calc_nlnstifflmass) soh8_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case ELEMENTS::struct_calc_stress:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      Teuchos::RCP<std::vector<char>> stressdata =
          params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
      Teuchos::RCP<std::vector<char>> straindata =
          params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
      Teuchos::RCP<std::vector<char>> plstraindata =
          params.get<Teuchos::RCP<std::vector<char>>>("plstrain", Teuchos::null);

      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'dB_ans_locisplacement'");
      if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get stress 'data'");
      if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get strain 'data'");
      if (plstraindata == Teuchos::null) FOUR_C_THROW("Cannot get plastic strain 'data'");

      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);

      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
      CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> plstrain;

      auto iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
          params, "iostress", INPAR::STR::stress_none);
      auto iostrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
          params, "iostrain", INPAR::STR::strain_none);
      auto ioplstrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
          params, "ioplstrain", INPAR::STR::strain_none);

      // decide whether evaluate 'thin' sosh stiff or 'thick' so_hex8 stiff
      if (eastype_ != DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        sosh8_nlnstiffmass(lm, mydisp, myres, nullptr, nullptr, nullptr, nullptr, &stress, &strain,
            params, iostress, iostrain);
      }
      else
      {
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          CORE::FE::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr, nullptr,
            nullptr, nullptr, &stress, &strain, &plstrain, params, iostress, iostrain, ioplstrain);
      }
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

    case ELEMENTS::struct_calc_eleload:
      FOUR_C_THROW("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    case ELEMENTS::struct_calc_fsiload:
      FOUR_C_THROW("Case not yet implemented");
      break;

    case ELEMENTS::struct_calc_update_istep:
    {
      // update internal EAS parameters
      if (eastype_ == soh8_eassosh8)
      {
        const auto* alpha = &easdata_.alpha;  // Alpha_{n+1}
        auto* alphao = &easdata_.alphao;      // Alpha_n
        // alphao := alpha
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(*alphao, *alpha);

        // store the EAS matrices
        const auto* Kaainv = &easdata_.invKaa;  // Kaa^{-1}_{n+1}
        auto* Kaainvo = &easdata_.invKaao;      // Kaa^{-1}_{n}
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, soh8_eassosh8>(
            *Kaainvo, *Kaainv);

        const auto* Kda = &easdata_.Kda;  // Kda_{n+1}
        auto* Kdao = &easdata_.Kdao;      // Kda_{n}
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, NUMDOF_SOH8>(*Kdao, *Kda);

        // reset EAS internal force
        CORE::LINALG::SerialDenseMatrix* oldfeas = &easdata_.feas;
        oldfeas->putScalar(0.0);
      }
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case ELEMENTS::struct_calc_reset_istep:
    {
      // restore internal EAS parameters
      if (eastype_ == soh8_eassosh8)
      {
        auto* alpha = &easdata_.alpha;          // Alpha_{n+1}
        const auto* alphao = &easdata_.alphao;  // Alpha_n
        // alphao := alphao
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(*alpha, *alphao);

        // store the EAS matrices
        auto* Kaainv = &easdata_.invKaa;          // Kaa^{-1}_{n+1}
        const auto* Kaainvo = &easdata_.invKaao;  // Kaa^{-1}_{n}
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, soh8_eassosh8>(
            *Kaainv, *Kaainvo);

        const auto* Kda = &easdata_.Kda;  // Kda_{n+1}
        auto* Kdao = &easdata_.Kdao;      // Kda_{n}
        CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, NUMDOF_SOH8>(*Kdao, *Kda);

        // reset EAS internal force
        CORE::LINALG::SerialDenseMatrix* oldfeas = &easdata_.feas;
        oldfeas->putScalar(0.0);
      }
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    case ELEMENTS::multi_calc_dens:
    {
      soh8_homog(params);
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // allocation and initializiation of these data arrays can only be
    // done in the elements that know the number of EAS parameters
    case ELEMENTS::multi_init_eas:
    {
      if (eastype_ != soh8_easnone)
      {
        soh8_eas_init_multi(params);
      }
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // before any microscale simulation, EAS internal data has to be
    // set accordingly
    case ELEMENTS::multi_set_eas:
    {
      if (eastype_ != soh8_easnone)
      {
        soh8_set_eas_multi(params);
      }
    }
    break;

    // read restart of microscale
    case ELEMENTS::multi_readrestart:
    {
      Teuchos::RCP<CORE::MAT::Material> mat = Material();

      if (mat->MaterialType() == CORE::Materials::m_struct_multiscale) soh8_read_restart_multi();
    }
    break;

    case ELEMENTS::shell_calc_stc_matrix:
    {
      const auto stc_scaling = CORE::UTILS::GetAsEnum<INPAR::STR::StcScale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        FOUR_C_THROW(
            "Action demands to calculate the STC (Scaled Thickness "
            "Conditiong) matrix, but not suitable scaling has been provided.");
      else
      {
        CalcSTCMatrix(
            elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, false);
      }
    }
    break;
    case ELEMENTS::shell_calc_stc_matrix_inverse:
    {
      const auto stc_scaling = CORE::UTILS::GetAsEnum<INPAR::STR::StcScale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        FOUR_C_THROW(
            "Action demands to calculate the STC (Scaled Thickness "
            "Conditiong) matrix, but not suitable scaling has been provided.");
      else
      {
        CalcSTCMatrix(elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, true);
      }
    }
    break;
    case ELEMENTS::struct_calc_recover:
      SoHex8::Evaluate(params, discretization, lm, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
      break;
    case ELEMENTS::struct_calc_energy:
    {
      if (eastype_ == DRT::ELEMENTS::SoHex8::soh8_easmild)
      {
        SoHex8::Evaluate(params, discretization, lm, elemat1_epetra, elemat2_epetra, elevec1_epetra,
            elevec2_epetra, elevec3_epetra);
        return 0;
      }

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);


      if (IsParamsInterface())  // new structural time integration
      {
        StrParamsInterface().AddContributionToEnergyType(
            sosh8_calc_energy(mydisp, params), STR::internal_energy);
      }
      else  // old structural time integration
      {
        // check length of elevec1
        if (elevec1_epetra.length() < 1) FOUR_C_THROW("The given result vector is too short.");

        elevec1_epetra(0) = sosh8_calc_energy(mydisp, params);
      }

      break;
    }
    case ELEMENTS::struct_calc_predict:
    {
      // do nothing here
      break;
    }
    case ELEMENTS::struct_create_backup:
    case ELEMENTS::struct_recover_from_backup:
    case ELEMENTS::struct_calc_mass_volume:
    case ELEMENTS::analyse_jacobian_determinant:
    {
      SoHex8::Evaluate(params, discretization, lm, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of action for So_sh8: %s | %d",
          ELEMENTS::ActionType2String(act).c_str(), act);
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::SoSh8::sosh8_calc_energy(
    const std::vector<double>& disp, Teuchos::ParameterList& params)
{
  if (PRESTRESS::IsMulf(pstype_)) FOUR_C_THROW("mulf is unsupported for the So_sh8 element!");

  if (kintype_ != INPAR::STR::KinemType::nonlinearTotLag)
    FOUR_C_THROW("Unsupported kinematic type for the So_sh8 element!");

  // initialization of internal energy
  double intenergy = 0.0;

  // shape functions and Gauss weights
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();

  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOH8 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOH8 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOH8 + 2];
  }

  // safety check before the actual evaluation starts
  const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
  if (min_detJ_curr <= 0.0)
  {
    soh8_error_handling(
        min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
    return 0.0;
  }

  // ------------------- EAS-SETUP --------------------------------------------
  /* EAS Technology: declare, intialize, set up, and alpha history */
  // in any case declare variables, sizes etc. only in eascase
  CORE::LINALG::SerialDenseMatrix* alpha = nullptr;              // EAS alphas
  std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, soh8_eassosh8>
      M;         // EAS matrix M at current GP, fixed for sosh8
  double detJ0;  // detJ(origin)
  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix

  switch (eastype_)
  {
    case soh8_eassosh8:
    {
      /*
      ** EAS Update of alphas:
      ** the current alphas are (re-)evaluated out of
      ** Kaa and Kda of previous step to avoid additional element call.
      ** This corresponds to the (innermost) element update loop
      ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
      */
      alpha = &easdata_.alpha;  // get old alpha

      /* evaluation of EAS variables (which are constant for the following):
      ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
      ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
      ** -> T0^{-T}
      */
      soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);

      break;
    }
    case soh8_easnone:
    {
      //      std::cout << "Warning: Solid-Shell8 without EAS" << std::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Solid-Shell8 only with eas_sosh8");
      exit(EXIT_FAILURE);
    }
  }
  // ------------------- END EAS-SETUP ----------------------------------------

  // ------------------- ANS-SETUP --------------------------------------------
  /* ANS Element technology to remedy
   * - transverse-shear locking E_rt and E_st
   * - trapezoidal (curvature-thickness) locking E_tt */
  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  // CORE::LINALG::Matrix<num_ans*num_sp,NUMDOF_SOH8> B_ans_loc(true); //set to 0
  CORE::LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> jac_sps(num_sp);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> jac_cur_sps(num_sp);
  // pointer to derivs evaluated at all sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>* deriv_sp =
      nullptr;  // derivs eval. at all sampling points
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe, xcurr, &deriv_sp, jac_sps, jac_cur_sps, B_ans_loc);

  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  //  static const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
  //  static const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc,
  //  gploc,-gploc}; static const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc,
  //  gploc,-gploc,-gploc, gploc, gploc};
  static const double* r = soh8_get_coordinate_of_gausspoints(0);
  static const double* s = soh8_get_coordinate_of_gausspoints(1);

  // ------------------- END ANS-SETUP ----------------------------------------

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac;
    double detJ = 0.0;
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac_cur;
    double detJ_cur = 0.0;

    if (not sosh8_evaluatejacobians(gp, derivs, xrefe, xcurr, jac, detJ, jac_cur, detJ_cur))
    {
      soh8_error_handling(
          detJ_cur, params, __LINE__, STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
      return 0.0;
    }

    // set up B-Operator in local(parameter) element space including ANS
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop_loc;
    sosh8_get_bop_loc(gp, derivs, jac_cur, r, s, B_ans_loc, bop_loc);

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> TinvT;
    sosh8_evaluateT(jac, TinvT);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
    bop.Multiply(TinvT, bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> lstrain;
    sosh8_get_glstrain_loc(gp, jac_cur, jac, jac_sps, jac_cur_sps, r, s, lstrain);

    // transformation of local glstrains 'back' to global(material) space
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(true);
    glstrain.Multiply(TinvT, lstrain);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
          soh8_eassosh8>(M.A(), detJ0 / detJ, T0invT.A(), M_GP->at(gp).values());
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
          1.0, glstrain.A(), 1.0, M.A(), (*alpha).values());
    }  // ------------------------------------------------------------------ EAS

    const double I3 = sosh8_third_invariant(glstrain);
    if (I3 <= 0.0)
    {
      soh8_error_handling(
          I3, params, __LINE__, STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
      return 0.0;
    }

    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd;
    sosh8_get_deformationgradient(gp, derivs, xcurr, glstrain, defgrd);

    // call material for evaluation of strain energy function
    double psi = 0.0;
    SolidMaterial()->StrainEnergy(glstrain, psi, gp, Id());

    const double detJ_w = detJ * gpweights[gp];

    // sum up GP contribution to internal energy
    intenergy += detJ_w * psi;
  }

  return intenergy;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                       // current displacements
    std::vector<double>& residual,                                   // current residual displ
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,     // element stiffness matrix
    CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,      // element mass matrix
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,      // element internal force vector
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,  // element structural force vector
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  // stress output option
    const INPAR::STR::StrainType iostrain)  // strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOH8 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOH8 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOH8 + 2];
  }

  const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
  if (min_detJ_curr <= 0.0)
  {
    soh8_error_handling(
        min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
    return;
  }

  // -------- EAS-SETUP -------------------------------------------------------
  /* EAS Technology: declare, intialize, set up, and alpha history */
  // in any case declare variables, sizes etc. only in eascase
  CORE::LINALG::SerialDenseMatrix* alpha = nullptr;              // EAS alphas
  std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, soh8_eassosh8>
      M;                                 // EAS matrix M at current GP, fixed for sosh8
  CORE::LINALG::SerialDenseVector feas;  // EAS portion of internal forces
  CORE::LINALG::SerialDenseMatrix Kaa;   // EAS matrix Kaa
  CORE::LINALG::SerialDenseMatrix Kda;   // EAS matrix Kda
  double detJ0;                          // detJ(origin)
  CORE::LINALG::SerialDenseMatrix* oldfeas = nullptr;    // EAS history
  CORE::LINALG::SerialDenseMatrix* oldKaainv = nullptr;  // EAS history
  CORE::LINALG::SerialDenseMatrix* oldKda = nullptr;     // EAS history
  CORE::LINALG::SerialDenseMatrix* eas_inc = nullptr;    // EAS increment

  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix

  switch (eastype_)
  {
    case soh8_eassosh8:
    {
      /*
      ** EAS Update of alphas:
      ** the current alphas are (re-)evaluated out of
      ** Kaa and Kda of previous step to avoid additional element call.
      ** This corresponds to the (innermost) element update loop
      ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
      */
      alpha = &easdata_.alpha;  // get old alpha
      // evaluate current (updated) EAS alphas (from history variables)
      // get stored EAS history
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
        CORE::LINALG::SerialDenseVector res_d(NUMDOF_SOH8);
        for (int i = 0; i < NUMDOF_SOH8; ++i)
        {
          res_d(i) = residual[i];
        }
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
        else
        {
          // add Kda . res_d to feas
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, NUMDOF_SOH8, 1>(
              1.0, *oldfeas, 1.0, *oldKda, res_d);
          // "new" alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, soh8_eassosh8, 1>(
              0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
          CORE::LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(
              1., alpha->values(), 1., eas_inc->values());
        }
      }
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

      break;
    }
    case soh8_easnone:
    {
      //      std::cout << "Warning: Solid-Shell8 without EAS" << std::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Solid-Shell8 only with eas_sosh8");
      exit(EXIT_FAILURE);
    }
  }
  // -------- END EAS-SETUP ---------------------------------------------------

  // -------- ANS-SETUP -------------------------------------------------------
  /* ANS Element technology to remedy
   *  - transverse-shear locking E_rt and E_st
   *  - trapezoidal (curvature-thickness) locking E_tt */
  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  // CORE::LINALG::Matrix<num_ans*num_sp,NUMDOF_SOH8> B_ans_loc(true); //set to 0
  CORE::LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> jac_sps(num_sp);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> jac_cur_sps(num_sp);
  // pointer to derivs evaluated at all sampling points
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>* deriv_sp =
      nullptr;  // derivs eval. at all sampling points
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe, xcurr, &deriv_sp, jac_sps, jac_cur_sps, B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double* r = soh8_get_coordinate_of_gausspoints(0);
  const double* s = soh8_get_coordinate_of_gausspoints(1);
  // -------- END ANS-SETUP ---------------------------------------------------

  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = params.isParameter("cond_rhs_norm");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac;
    double detJ = 0.0;
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac_cur;
    double detJ_cur = 0.0;

    if (not sosh8_evaluatejacobians(gp, derivs, xrefe, xcurr, jac, detJ, jac_cur, detJ_cur))
    {
      soh8_error_handling(
          detJ_cur, params, __LINE__, STR::ELEMENTS::ele_error_negative_det_of_def_gradient);

      if (stiffmatrix) stiffmatrix->Clear();
      if (force) force->Clear();

      return;
    }

    // set up B-Operator in local(parameter) element space including ANS
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop_loc;
    sosh8_get_bop_loc(gp, derivs, jac_cur, r, s, B_ans_loc, bop_loc);

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> TinvT;
    sosh8_evaluateT(jac, TinvT);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
    bop.Multiply(TinvT, bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> lstrain;
    sosh8_get_glstrain_loc(gp, jac_cur, jac, jac_sps, jac_cur_sps, r, s, lstrain);

    // transformation of local glstrains 'back' to global(material) space
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(true);
    glstrain.Multiply(TinvT, lstrain);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
          soh8_eassosh8>(M.A(), detJ0 / detJ, T0invT.A(), M_GP->at(gp).values());
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
          1.0, glstrain.A(), 1.0, M.A(), (*alpha).values());
    }  // ------------------------------------------------------------------ EAS

    const double I3 = sosh8_third_invariant(glstrain);
    if (I3 <= 0.0)
    {
      soh8_error_handling(
          I3, params, __LINE__, STR::ELEMENTS::ele_error_negative_det_of_def_gradient);

      if (stiffmatrix) stiffmatrix->Clear();
      if (force) force->Clear();

      return;
    }


    // return gp strains if necessary
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");

        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);

        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);

        break;
      }
      case INPAR::STR::strain_ea:
      {
        FOUR_C_THROW("no Euler-Almansi strains available for sosh8");
        break;
      }
      case INPAR::STR::strain_none:
      {
        break;
      }
      default:
      {
        FOUR_C_THROW("requested strain option not available");
        exit(EXIT_FAILURE);
      }
    }

    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd;
    sosh8_get_deformationgradient(gp, derivs, xcurr, glstrain, defgrd);
    const double det_defgrd = defgrd.Determinant();
    if (det_defgrd <= 0.0)
    {
      soh8_error_handling(
          det_defgrd, params, __LINE__, STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
      return;
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses if necessary
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
        {
          (*elestress)(gp, i) = stress(i);
        }

        break;
      }
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        sosh8_Cauchy(elestress, gp, defgrd, glstrain, stress);

        break;
      }
      case INPAR::STR::stress_none:
        break;
      default:
        FOUR_C_THROW("requested stress option not available");
        exit(EXIT_FAILURE);
    }

    const double detJ_w = detJ * gpweights[gp];
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }  // if (force!=nullptr)

    // structural force vector
    if (split_res) force_str->MultiplyTN(detJ_w, bop, stress, 1.0);

    // update stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> cb;
      cb.Multiply(cmat, bop);  // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // intergrate `geometric' stiffness matrix and add to keu *****************
      // here also the ANS interpolation comes into play
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> G_ij;
          G_ij(0) = derivs[gp](0, inod) * derivs[gp](0, jnod);  // rr-dir
          G_ij(1) = derivs[gp](1, inod) * derivs[gp](1, jnod);  // ss-dir
          G_ij(3) = derivs[gp](0, inod) * derivs[gp](1, jnod) +
                    derivs[gp](1, inod) * derivs[gp](0, jnod);  // rs-dir

          // do the ANS related stuff if wanted!
          if (anstype_ == anssosh8)
          {
            // ANS modification in tt-dir
            G_ij(2) = 0.25 * (1 - r[gp]) * (1 - s[gp]) * (*deriv_sp)[4](2, inod) *
                          (*deriv_sp)[4](2, jnod) +
                      0.25 * (1 + r[gp]) * (1 - s[gp]) * (*deriv_sp)[5](2, inod) *
                          (*deriv_sp)[5](2, jnod) +
                      0.25 * (1 + r[gp]) * (1 + s[gp]) * (*deriv_sp)[6](2, inod) *
                          (*deriv_sp)[6](2, jnod) +
                      0.25 * (1 - r[gp]) * (1 + s[gp]) * (*deriv_sp)[7](2, inod) *
                          (*deriv_sp)[7](2, jnod);
            // ANS modification in st-dir
            G_ij(4) =
                0.5 * ((1 + r[gp]) * ((*deriv_sp)[1](1, inod) * (*deriv_sp)[1](2, jnod) +
                                         (*deriv_sp)[1](2, inod) * (*deriv_sp)[1](1, jnod)) +
                          (1 - r[gp]) * ((*deriv_sp)[3](1, inod) * (*deriv_sp)[3](2, jnod) +
                                            (*deriv_sp)[3](2, inod) * (*deriv_sp)[3](1, jnod)));
            // ANS modification in rt-dir
            G_ij(5) =
                0.5 * ((1 - s[gp]) * ((*deriv_sp)[0](0, inod) * (*deriv_sp)[0](2, jnod) +
                                         (*deriv_sp)[0](2, inod) * (*deriv_sp)[0](0, jnod)) +
                          (1 + s[gp]) * ((*deriv_sp)[2](0, inod) * (*deriv_sp)[2](2, jnod) +
                                            (*deriv_sp)[2](2, inod) * (*deriv_sp)[2](0, jnod)));
          }
          else if (anstype_ == ansnone)
          {
            G_ij(2) = derivs[gp](2, inod) * derivs[gp](2, jnod);  // tt-dir
            G_ij(4) = derivs[gp](2, inod) * derivs[gp](1, jnod) +
                      derivs[gp](1, inod) * derivs[gp](2, jnod);  // st-dir
            G_ij(5) = derivs[gp](0, inod) * derivs[gp](2, jnod) +
                      derivs[gp](2, inod) * derivs[gp](0, jnod);  // rt-dir
          }
          else
            FOUR_C_THROW("Cannot build geometric stiffness matrix on your ANS-choice!");

          // transformation of local(parameter) space 'back' to global(material) space
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          const double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += Gij;
          (*stiffmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += Gij;
          (*stiffmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += Gij;
        }
      }  // end of intergrate `geometric' stiffness ******************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone)
      {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, soh8_eassosh8> cM;  // temporary c . M
        cM.Multiply(cmat, M);
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
            soh8_eassosh8>(1.0, Kaa.values(), detJ_w, M.A(), cM.A());
        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
            NUMDOF_SOH8>(1.0, Kda.values(), detJ_w, M.A(), cb.A());
        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D, 1>(
            1.0, feas.values(), detJ_w, M.A(), stress.A());
      }  // ------------------------------------------------------------------ EAS
    }    // if (stiffmatrix != nullptr)

    if (massmatrix != nullptr)
    {  // evaluate mass matrix +++++++++++++++++++++++++
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
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
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  // rhs norm of eas equations
  if (eastype_ != soh8_easnone && split_res)
    // only add for row-map elements
    if (params.get<int>("MyPID") == Owner())
      params.get<double>("cond_rhs_norm") += pow(CORE::LINALG::Norm2(feas), 2.);

  if (force != nullptr && stiffmatrix != nullptr)
  {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ == soh8_eassosh8)
    {
      // we need the inverse of Kaa
      using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
      using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
      Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_inverseKaa;
      solve_for_inverseKaa.setMatrix(Teuchos::rcpFromRef(Kaa));
      solve_for_inverseKaa.invert();

      CORE::LINALG::SerialDenseMatrix KdaTKaa(
          NUMDOF_SOH8, soh8_eassosh8);  // temporary Kda^T.Kaa^{-1}
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_eassosh8, soh8_eassosh8>(
          KdaTKaa, Kda, Kaa);
      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, NUMDOF_SOH8>(
          1.0, stiffmatrix->A(), -1.0, KdaTKaa.values(), Kda.values());
      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, 1>(
          1.0, force->A(), -1.0, KdaTKaa.values(), feas.values());

      // store current EAS data in history
      for (int i = 0; i < soh8_eassosh8; ++i)
      {
        for (int j = 0; j < soh8_eassosh8; ++j) (*oldKaainv)(i, j) = Kaa(i, j);
        for (int j = 0; j < NUMDOF_SOH8; ++j) (*oldKda)(i, j) = Kda(i, j);
        (*oldfeas)(i, 0) = feas(i);
      }
    }  // -------------------------------------------------------------------- EAS
  }

  return;
}  // DRT::ELEMENTS::So_sh8::sosh8_nlnstiffmass



/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_anssetup(
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe,  // material element coords
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,  // current element coords
    std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>**
        deriv_sp,  // derivs eval. at all sampling points
    std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>&
        jac_sps,  // jac at all sampling points
    std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>&
        jac_cur_sps,  // current jac at all sampling points
    CORE::LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8>& B_ans_loc) const  // modified B
{
  // static matrix object of derivs at sampling points, kept in memory
  static std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> df_sp(num_sp);
  static bool dfsp_eval;  // flag for re-evaluate everything

  if (dfsp_eval != 0)
  {                      // if true f,df already evaluated
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
  }
  else
  {
    /*====================================================================*/
    /* 8-node hexhedra Solid-Shell node topology
     * and location of sampling points A to H                             */
    /*--------------------------------------------------------------------*/
    /*                      t
     *                      |
     *             4========|================7
     *          // |        |              //||
     *        //   |        |            //  ||
     *      //     |        |   D      //    ||
     *     5=======E=================6       H
     *    ||       |        |        ||      ||
     *    ||   A   |        o--------||-- C -------s
     *    ||       |       /         ||      ||
     *    F        0----- B ---------G ------3
     *    ||     //     /            ||    //
     *    ||   //     /              ||  //
     *    || //     r                ||//
     *     1=========================2
     *
     */
    /*====================================================================*/
    // (r,s,t) gp-locations of sampling points A,B,C,D,E,F,G,H
    // numsp = 8 here set explicitly to allow direct initializing
    //                A,   B,   C,   D,   E,   F,   G,   H
    std::array<double, 8> r = {0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    std::array<double, 8> s = {-1.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0, 1.0};
    std::array<double, 8> t = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i = 0; i < num_sp; ++i)
    {
      // df wrt to r "+0" for each node(0..7) at each sp [i]
      df_sp[i](0, 0) = -(1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 1) = (1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 2) = (1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 3) = -(1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 4) = -(1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 5) = (1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 6) = (1.0 + s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 7) = -(1.0 + s[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to s "+1" for each node(0..7) at each sp [i]
      df_sp[i](1, 0) = -(1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 1) = -(1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 2) = (1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 3) = (1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 4) = -(1.0 - r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 5) = -(1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 6) = (1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 7) = (1.0 - r[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to t "+2" for each node(0..7) at each sp [i]
      df_sp[i](2, 0) = -(1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 1) = -(1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 2) = -(1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 3) = -(1.0 - r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 4) = (1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 5) = (1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 6) = (1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 7) = (1.0 - r[i]) * (1.0 + s[i]) * 0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
    dfsp_eval = true;    // now all arrays are filled statically
  }

  for (int sp = 0; sp < num_sp; ++sp)
  {
    // compute (REFERENCE) Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp], xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp], xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac_cur;
  for (int sp = 0; sp < num_sp; ++sp)
  {
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    jac_cur.Multiply(df_sp[sp], xcurr);

    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_SOH8; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
      {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp * num_ans + 0, inode * 3 + dim) = df_sp[sp](2, inode) * jac_cur(2, dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp * num_ans + 1, inode * 3 + dim) =
            df_sp[sp](1, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(1, dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp * num_ans + 2, inode * 3 + dim) =
            df_sp[sp](0, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(0, dim);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoSh8::sosh8_evaluatejacobians(const unsigned gp,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe,
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_ref, double& detJ_ref,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr, double& detJ_curr) const
{
  sosh8_evaluatejacobian(gp, derivs, xrefe, jac_ref, detJ_ref);
  if (detJ_ref == 0.0)
    FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
  else if (detJ_ref < 0.0)
    FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

  /* compute the CURRENT Jacobian matrix which looks like:
  **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
  **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
  **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
  ** Used to transform the global displacements into parametric space */
  sosh8_evaluatejacobian(gp, derivs, xcurr, jac_curr, detJ_curr);

  return (detJ_curr > 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_evaluatejacobian(const unsigned gp,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& x,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac, double& detJ) const
{
  /* compute the Jacobian matrix which looks like:
  **         [ x_,r  y_,r  z_,r ]
  **     J = [ x_,s  y_,s  z_,s ]
  **         [ x_,t  y_,t  z_,t ]
  */
  jac.Multiply(derivs[gp], x);

  // compute determinant of Jacobian by Sarrus' rule
  detJ = jac.Determinant();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::SoSh8::sosh8_third_invariant(
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain) const
{
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> rcg(glstrain);
  rcg.Scale(2.0);
  for (unsigned i = 0; i < 3; ++i) rcg(i) += 1.0;

  // compute the 3rd invariant, a.k.a. the square product of the det(defGrad)
  const double I3 = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                    0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                    0.25 * rcg(0) * rcg(4) * rcg(4);

  return I3;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_get_bop_loc(const unsigned gp,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr, const double* r,
    const double* s, const CORE::LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8>& B_ans_loc,
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8>& bop_loc) const
{
  // set up B-Operator in local(parameter) element space including ANS
  for (int inode = 0; inode < NUMNOD_SOH8; ++inode)
  {
    for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
    {
      // B_loc_rr = N_r.X_r
      bop_loc(0, inode * 3 + dim) = derivs[gp](0, inode) * jac_curr(0, dim);
      // B_loc_ss = N_s.X_s
      bop_loc(1, inode * 3 + dim) = derivs[gp](1, inode) * jac_curr(1, dim);
      // B_loc_rs = N_r.X_s + N_s.X_r
      bop_loc(3, inode * 3 + dim) =
          derivs[gp](0, inode) * jac_curr(1, dim) + derivs[gp](1, inode) * jac_curr(0, dim);

      // do the ANS related stuff
      switch (anstype_)
      {
        case anssosh8:
        {
          // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
          //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
          //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
          bop_loc(2, inode * 3 + dim) =
              0.25 * (1 - r[gp]) * (1 - s[gp]) * B_ans_loc(0 + 4 * num_ans, inode * 3 + dim)    // E
              + 0.25 * (1 + r[gp]) * (1 - s[gp]) * B_ans_loc(0 + 5 * num_ans, inode * 3 + dim)  // F
              + 0.25 * (1 + r[gp]) * (1 + s[gp]) * B_ans_loc(0 + 6 * num_ans, inode * 3 + dim)  // G
              +
              0.25 * (1 - r[gp]) * (1 + s[gp]) * B_ans_loc(0 + 7 * num_ans, inode * 3 + dim);  // H
          // B_loc_st = interpolation along r of ANS B_loc_st
          //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
          bop_loc(4, inode * 3 + dim) =
              0.5 * (1.0 + r[gp]) * B_ans_loc(1 + 1 * num_ans, inode * 3 + dim)     // B
              + 0.5 * (1.0 - r[gp]) * B_ans_loc(1 + 3 * num_ans, inode * 3 + dim);  // D

          // B_loc_rt = interpolation along s of ANS B_loc_rt
          //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
          bop_loc(5, inode * 3 + dim) =
              0.5 * (1.0 - s[gp]) * B_ans_loc(2 + 0 * num_ans, inode * 3 + dim)     // A
              + 0.5 * (1.0 + s[gp]) * B_ans_loc(2 + 2 * num_ans, inode * 3 + dim);  // C

          break;
        }
        case ansnone:
        {
          // B_loc_tt = N_t.X_t
          bop_loc(2, inode * 3 + dim) = derivs[gp](2, inode) * jac_curr(2, dim);
          // B_loc_st = N_t.X_s + N_s.X_t
          bop_loc(4, inode * 3 + dim) =
              derivs[gp](2, inode) * jac_curr(1, dim) + derivs[gp](1, inode) * jac_curr(2, dim);

          // B_loc_rt = N_r.X_t + N_t.X_r
          bop_loc(5, inode * 3 + dim) =
              derivs[gp](0, inode) * jac_curr(2, dim) + derivs[gp](2, inode) * jac_curr(0, dim);

          break;
        }
        default:
        {
          FOUR_C_THROW("Cannot build bop_loc based on your ANS-choice!");
          exit(EXIT_FAILURE);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_get_glstrain_loc(const unsigned gp,
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr,
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>& jac_sps,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>& jac_cur_sps, const double* r,
    const double* s, CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& lstrain) const
{
  // evaluate glstrains in local(parameter) coords
  // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
  lstrain(0) = 0.5 * (+(jac_curr(0, 0) * jac_curr(0, 0) + jac_curr(0, 1) * jac_curr(0, 1) +
                          jac_curr(0, 2) * jac_curr(0, 2)) -
                         (jac(0, 0) * jac(0, 0) + jac(0, 1) * jac(0, 1) + jac(0, 2) * jac(0, 2)));
  // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
  lstrain(1) = 0.5 * (+(jac_curr(1, 0) * jac_curr(1, 0) + jac_curr(1, 1) * jac_curr(1, 1) +
                          jac_curr(1, 2) * jac_curr(1, 2)) -
                         (jac(1, 0) * jac(1, 0) + jac(1, 1) * jac(1, 1) + jac(1, 2) * jac(1, 2)));
  // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
  lstrain(3) = (+(jac_curr(0, 0) * jac_curr(1, 0) + jac_curr(0, 1) * jac_curr(1, 1) +
                    jac_curr(0, 2) * jac_curr(1, 2)) -
                (jac(0, 0) * jac(1, 0) + jac(0, 1) * jac(1, 1) + jac(0, 2) * jac(1, 2)));


  // do the ANS related stuff if wanted!
  switch (anstype_)
  {
    case anssosh8:
    {
      // ANS modification of strains ************************************** ANS
      double dxdt_A = 0.0;
      double dXdt_A = 0.0;
      double dydt_B = 0.0;
      double dYdt_B = 0.0;
      double dxdt_C = 0.0;
      double dXdt_C = 0.0;
      double dydt_D = 0.0;
      double dYdt_D = 0.0;

      double dzdt_E = 0.0;
      double dZdt_E = 0.0;
      double dzdt_F = 0.0;
      double dZdt_F = 0.0;
      double dzdt_G = 0.0;
      double dZdt_G = 0.0;
      double dzdt_H = 0.0;
      double dZdt_H = 0.0;

      // vector product of rows of jacobians at corresponding sampling point    std::cout <<
      // jac_cur_sps;
      for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
      {
        dxdt_A += jac_cur_sps[0](0, dim) * jac_cur_sps[0](2, dim);  // g_13^A
        dXdt_A += jac_sps[0](0, dim) * jac_sps[0](2, dim);          // G_13^A
        dydt_B += jac_cur_sps[1](1, dim) * jac_cur_sps[1](2, dim);  // g_23^B
        dYdt_B += jac_sps[1](1, dim) * jac_sps[1](2, dim);          // G_23^B
        dxdt_C += jac_cur_sps[2](0, dim) * jac_cur_sps[2](2, dim);  // g_13^C
        dXdt_C += jac_sps[2](0, dim) * jac_sps[2](2, dim);          // G_13^C
        dydt_D += jac_cur_sps[3](1, dim) * jac_cur_sps[3](2, dim);  // g_23^D
        dYdt_D += jac_sps[3](1, dim) * jac_sps[3](2, dim);          // G_23^D

        dzdt_E += jac_cur_sps[4](2, dim) * jac_cur_sps[4](2, dim);
        dZdt_E += jac_sps[4](2, dim) * jac_sps[4](2, dim);
        dzdt_F += jac_cur_sps[5](2, dim) * jac_cur_sps[5](2, dim);
        dZdt_F += jac_sps[5](2, dim) * jac_sps[5](2, dim);
        dzdt_G += jac_cur_sps[6](2, dim) * jac_cur_sps[6](2, dim);
        dZdt_G += jac_sps[6](2, dim) * jac_sps[6](2, dim);
        dzdt_H += jac_cur_sps[7](2, dim) * jac_cur_sps[7](2, dim);
        dZdt_H += jac_sps[7](2, dim) * jac_sps[7](2, dim);
      }
      // E33: remedy of curvature thickness locking
      // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
      lstrain(2) = 0.5 * (0.25 * (1 - r[gp]) * (1 - s[gp]) * (dzdt_E - dZdt_E) +
                             0.25 * (1 + r[gp]) * (1 - s[gp]) * (dzdt_F - dZdt_F) +
                             0.25 * (1 + r[gp]) * (1 + s[gp]) * (dzdt_G - dZdt_G) +
                             0.25 * (1 - r[gp]) * (1 + s[gp]) * (dzdt_H - dZdt_H));
      // E23: remedy of transverse shear locking
      // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
      lstrain(4) = 0.5 * (1 + r[gp]) * (dydt_B - dYdt_B) + 0.5 * (1 - r[gp]) * (dydt_D - dYdt_D);
      // E13: remedy of transverse shear locking
      // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
      lstrain(5) = 0.5 * (1 - s[gp]) * (dxdt_A - dXdt_A) + 0.5 * (1 + s[gp]) * (dxdt_C - dXdt_C);
      // ANS modification of strains ************************************** ANS

      break;
    }
    case ansnone:
    {
      // No ANS!
      // Ett = 0.5 * (dz/dt * dz/dt^T - dZ/dt * dZ/dt^T)
      lstrain(2) =
          0.5 * (+(jac_curr(2, 0) * jac_curr(2, 0) + jac_curr(2, 1) * jac_curr(2, 1) +
                     jac_curr(2, 2) * jac_curr(2, 2)) -
                    (jac(2, 0) * jac(2, 0) + jac(2, 1) * jac(2, 1) + jac(2, 2) * jac(2, 2)));
      // Est = (dz/ds * dy/dt^T - dZ/ds * dY/dt^T)
      lstrain(4) = (+(jac_curr(2, 0) * jac_curr(1, 0) + jac_curr(2, 1) * jac_curr(1, 1) +
                        jac_curr(2, 2) * jac_curr(1, 2)) -
                    (jac(2, 0) * jac(1, 0) + jac(2, 1) * jac(1, 1) + jac(2, 2) * jac(1, 2)));
      // Est = (dz/dr * dx/dt^T - dZ/dr * dX/dt^T)
      lstrain(5) = (+(jac_curr(2, 0) * jac_curr(0, 0) + jac_curr(2, 1) * jac_curr(0, 1) +
                        jac_curr(2, 2) * jac_curr(0, 2)) -
                    (jac(2, 0) * jac(0, 0) + jac(2, 1) * jac(0, 1) + jac(2, 2) * jac(0, 2)));

      break;
    }
    default:
    {
      FOUR_C_THROW("Cannot build local strains based on your ANS-choice!");
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_get_deformationgradient(const unsigned gp,
    const std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
    const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
    CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd) const
{
  /* Caution!! the defgrd can not be modified with ANS to remedy locking
     To get the consistent F a spectral decomposition would be necessary, see sosh8_Cauchy.
     However if one only maps e.g. stresses from current to material configuration,
     I have never noticed any difference to applying just the disp_based F
     which is therefore computed and passed here (no significant add. computation time).  */
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  N_XYZ.Multiply(invJ_[gp], derivs[gp]);
  // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
  defgrd.MultiplyTT(xcurr, N_XYZ);

  // deformation gradient consistent with (potentially EAS-modified) GL strains
  // without eas/ans this is equal to the regular defgrd.
  // This is necessary for material formulations based on the deformation
  // gradient rather than the GL strains.

  // calculate deformation gradient consistent with modified GL strain tensor
  if ((eastype_ != soh8_easnone || anstype_ == ansnone) &&
      (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->NeedsDefgrd()))
    CalcConsistentDefgrd(defgrd, glstrain, defgrd);
}

/*----------------------------------------------------------------------*
 |  evaluate 'T'-transformation matrix )                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_evaluateT(
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac,
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0, 0) = jac(0, 0) * jac(0, 0);
  TinvT(1, 0) = jac(1, 0) * jac(1, 0);
  TinvT(2, 0) = jac(2, 0) * jac(2, 0);
  TinvT(3, 0) = 2 * jac(0, 0) * jac(1, 0);
  TinvT(4, 0) = 2 * jac(1, 0) * jac(2, 0);
  TinvT(5, 0) = 2 * jac(0, 0) * jac(2, 0);

  TinvT(0, 1) = jac(0, 1) * jac(0, 1);
  TinvT(1, 1) = jac(1, 1) * jac(1, 1);
  TinvT(2, 1) = jac(2, 1) * jac(2, 1);
  TinvT(3, 1) = 2 * jac(0, 1) * jac(1, 1);
  TinvT(4, 1) = 2 * jac(1, 1) * jac(2, 1);
  TinvT(5, 1) = 2 * jac(0, 1) * jac(2, 1);

  TinvT(0, 2) = jac(0, 2) * jac(0, 2);
  TinvT(1, 2) = jac(1, 2) * jac(1, 2);
  TinvT(2, 2) = jac(2, 2) * jac(2, 2);
  TinvT(3, 2) = 2 * jac(0, 2) * jac(1, 2);
  TinvT(4, 2) = 2 * jac(1, 2) * jac(2, 2);
  TinvT(5, 2) = 2 * jac(0, 2) * jac(2, 2);

  TinvT(0, 3) = jac(0, 0) * jac(0, 1);
  TinvT(1, 3) = jac(1, 0) * jac(1, 1);
  TinvT(2, 3) = jac(2, 0) * jac(2, 1);
  TinvT(3, 3) = jac(0, 0) * jac(1, 1) + jac(1, 0) * jac(0, 1);
  TinvT(4, 3) = jac(1, 0) * jac(2, 1) + jac(2, 0) * jac(1, 1);
  TinvT(5, 3) = jac(0, 0) * jac(2, 1) + jac(2, 0) * jac(0, 1);


  TinvT(0, 4) = jac(0, 1) * jac(0, 2);
  TinvT(1, 4) = jac(1, 1) * jac(1, 2);
  TinvT(2, 4) = jac(2, 1) * jac(2, 2);
  TinvT(3, 4) = jac(0, 1) * jac(1, 2) + jac(1, 1) * jac(0, 2);
  TinvT(4, 4) = jac(1, 1) * jac(2, 2) + jac(2, 1) * jac(1, 2);
  TinvT(5, 4) = jac(0, 1) * jac(2, 2) + jac(2, 1) * jac(0, 2);

  TinvT(0, 5) = jac(0, 0) * jac(0, 2);
  TinvT(1, 5) = jac(1, 0) * jac(1, 2);
  TinvT(2, 5) = jac(2, 0) * jac(2, 2);
  TinvT(3, 5) = jac(0, 0) * jac(1, 2) + jac(1, 0) * jac(0, 2);
  TinvT(4, 5) = jac(1, 0) * jac(2, 2) + jac(2, 0) * jac(1, 2);
  TinvT(5, 5) = jac(0, 0) * jac(2, 2) + jac(2, 0) * jac(0, 2);

  // now evaluate T^{-T} with solver
  CORE::LINALG::FixedSizeSerialDenseSolver<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, 1>
      solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2 != 0)) FOUR_C_THROW("Inversion of Tinv (Jacobian) failed");
  return;
}

/*----------------------------------------------------------------------*
 |  return Cauchy stress at gp                                 maf 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8::sosh8_Cauchy(
    CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress, const int gp,
    const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress)
{
#if consistent_F
  // double disp1 = defgrd.NormOne();
  // double dispinf = defgrd.NormInf();

  /* to get the consistent (locking-free) F^mod, we need two spectral
   * compositions. First, find R (rotation tensor) from F=RU,
   * then from E^mod = 1/2((U^mod)^2 - 1) find U^mod,
   * and finally F^mod = RU^mod */

  // polar decomposition of displacement based F
  CORE::LINALG::SerialDenseMatrix u(NUMDIM_SOH8, NUMDIM_SOH8);
  CORE::LINALG::SerialDenseMatrix s(NUMDIM_SOH8, NUMDIM_SOH8);
  CORE::LINALG::SerialDenseMatrix v(NUMDIM_SOH8, NUMDIM_SOH8);
  SVD(defgrd, u, s, v);  // Singular Value Decomposition
  CORE::LINALG::SerialDenseMatrix rot(NUMDIM_SOH8, NUMDIM_SOH8);
  CORE::LINALG::multiply(rot, u, v);
  // temp.Multiply('N','N',1.0,v,s,0.0);
  // CORE::LINALG::SerialDenseMatrix stretch_disp(NUMDIM_SOH8,NUMDIM_SOH8);
  // stretch_disp.Multiply('N','T',1.0,temp,v,0.0);
  // defgrd.Multiply('N','N',1.0,rot,stretch_disp,0.0);
  // std::cout << defgrd;

  // get modified squared stretch (U^mod)^2 from glstrain
  CORE::LINALG::SerialDenseMatrix Usq_mod(NUMDIM_SOH8, NUMDIM_SOH8);
  for (int i = 0; i < NUMDIM_SOH8; ++i) Usq_mod(i, i) = 2.0 * glstrain(i) + 1.0;
  // off-diagonal terms are already twice in the Voigt-GLstrain-vector
  Usq_mod(0, 1) = glstrain(3);
  Usq_mod(1, 0) = glstrain(3);
  Usq_mod(1, 2) = glstrain(4);
  Usq_mod(2, 1) = glstrain(4);
  Usq_mod(0, 2) = glstrain(5);
  Usq_mod(2, 0) = glstrain(5);
  // polar decomposition of (U^mod)^2
  SVD(Usq_mod, u, s, v);  // Singular Value Decomposition
  CORE::LINALG::SerialDenseMatrix U_mod(NUMDIM_SOH8, NUMDIM_SOH8);
  for (int i = 0; i < NUMDIM_SOH8; ++i) s(i, i) = sqrt(s(i, i));
  CORE::LINALG::SerialDenseMatrix temp2(NUMDIM_SOH8, NUMDIM_SOH8);
  CORE::LINALG::multiply(temp2, u, s);
  CORE::LINALG::multiply(U_mod, temp2, v);

  // F^mod = RU^mod
  CORE::LINALG::SerialDenseMatrix defgrd_consistent(NUMDIM_SOH8, NUMDIM_SOH8);
  CORE::LINALG::multiply(defgrd_consistent, rot, U_mod);
  defgrd.SetView(defgrd_consistent.A());

  /*
  double mod1 = defgrd.NormOne();
  double modinf = defgrd.NormInf();
  if(((mod1-disp1)/mod1 > 0.03) || ((modinf-dispinf)/modinf > 0.03)){
    std::cout << "difference in F! mod1= " << mod1 << " disp1= " << disp1 << " modinf= " << modinf
  << " dispinf= " << dispinf << std::endl; std::cout << "Fmod" << std::endl << defgrd;
  }
  */
#endif

  double detF = defgrd.Determinant();

  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> pkstress;
  pkstress(0, 0) = stress(0);
  pkstress(0, 1) = stress(3);
  pkstress(0, 2) = stress(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = stress(1);
  pkstress(1, 2) = stress(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = stress(2);

  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchystress;
  CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
  temp.Multiply(1.0 / detF, defgrd, pkstress);
  cauchystress.MultiplyNT(temp, defgrd);

  (*elestress)(gp, 0) = cauchystress(0, 0);
  (*elestress)(gp, 1) = cauchystress(1, 1);
  (*elestress)(gp, 2) = cauchystress(2, 2);
  (*elestress)(gp, 3) = cauchystress(0, 1);
  (*elestress)(gp, 4) = cauchystress(1, 2);
  (*elestress)(gp, 5) = cauchystress(0, 2);

  return;
}

void DRT::ELEMENTS::SoSh8::CalcSTCMatrix(CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>& elemat1,
    const INPAR::STR::StcScale stc_scaling, const int stc_layer, std::vector<int>& lm,
    DRT::Discretization& discretization, bool calcinverse)
{
  /// Compute C based on element aspect ratio
  double stc_fact = 1.0;
  if (stc_scaling == INPAR::STR::stc_currsym)
  {
    // stc_fact = sqrt(sosh8_calcaspectratio());
    stc_fact = sosh8_calcaspectratio();
  }
  else
  {
    // stc_fact = sosh8_calcaspectratio();
    stc_fact = sosh8_calcaspectratio() * sosh8_calcaspectratio();
  }

  // Compute different scaling factors for STC or Inv(STC)
  double factor1 = 0.0;
  double factor2 = 0.0;
  double factor3 = 0.0;
  double factor4 = 0.0;
  if (!calcinverse)
  {
    factor1 = (stc_fact + 1.0) / (2.0 * stc_fact);
    factor2 = (stc_fact - 1.0) / (2.0 * stc_fact);
    factor3 = (1.0 / stc_fact);
    factor4 = (1.0 - 1.0 / stc_fact);
  }
  else
  {
    factor1 = (1.0 + stc_fact) / 2.0;
    factor2 = (1.0 - stc_fact) / 2.0;
    factor3 = stc_fact;
    factor4 = 1 - stc_fact;
  }

  if (stc_scaling == INPAR::STR::stc_curr or stc_scaling == INPAR::STR::stc_currsym)
  {
    CORE::LINALG::Matrix<NUMDOF_SOH8, 1> adjele(true);
    DRT::Node** nodes = Nodes();

    std::vector<DRT::Condition*> cond0;
    std::vector<DRT::Condition*> condFSI0;
    int condnum0 = 1000;    // minimun STCid of layer with nodes 0..3
    bool current0 = false;  // layer with nodes 0..4 to be scaled
    (nodes[0])->GetCondition("STC Layer", cond0);
    (nodes[0])->GetCondition("FSICoupling", condFSI0);
    std::vector<DRT::Condition*> cond1;
    std::vector<DRT::Condition*> condFSI1;
    int condnum1 = 1000;    // minimun STCid of layer with nodes 4..7
    bool current1 = false;  // minimun STCid of layer with nodes 4..7
    (nodes[NUMNOD_SOH8 / 2])->GetCondition("STC Layer", cond1);
    (nodes[NUMNOD_SOH8 / 2])->GetCondition("FSICoupling", condFSI1);

    for (auto& conu : cond0)
    {
      int tmp = conu->Get<int>("ConditionID");
      if (tmp < condnum0) condnum0 = tmp;
    }
    if (condnum0 ==
        stc_layer)  // && (condFSI0.size()==0 or (condFSI0.size()!=0 and condFSI1.size()!=0)))
      current0 = true;


    for (auto& conu : cond1)
    {
      int tmp = conu->Get<int>("ConditionID");
      if (tmp < condnum1) condnum1 = tmp;
    }
    if (condnum1 ==
        stc_layer)  // && (condFSI1.size()==0 or (condFSI0.size()!=0 and condFSI1.size()!=0)))
      current1 = true;


    // both surfaces are to be scaled
    if (current0 and current1)
    {
      // only valid for first round
      if (condnum0 != 1)
        FOUR_C_THROW("STC error: non-initial layer is not connected to a smaller id");
      else
      {
        for (int i = 0; i < NUMNOD_SOH8; i++)
        {
          adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
        }
        for (int ind1 = 0; ind1 < NUMDOF_SOH8 / 2; ind1++)
        {
          elemat1(ind1, ind1) += factor1 / adjele(ind1, 0) * cond0.size();
          elemat1(ind1 + NUMDOF_SOH8 / 2, ind1 + NUMDOF_SOH8 / 2) +=
              factor1 / adjele(ind1 + NUMDOF_SOH8 / 2, 0) * cond1.size();
          elemat1(ind1, ind1 + NUMDOF_SOH8 / 2) += factor2 / adjele(ind1, 0) * cond0.size();
          elemat1(ind1 + NUMDOF_SOH8 / 2, ind1) +=
              factor2 / adjele(ind1 + NUMDOF_SOH8 / 2, 0) * cond1.size();
        }
      }
    }
    // surface with nodes 0..3 is to be scaled
    else if (current0)
    {
      // but not by this element
      if (condnum1 > condnum0)
      {
        for (int i = 0; i < NUMNOD_SOH8; i++)
        {
          adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
        }
        for (int ind1 = NUMDOF_SOH8 / 2; ind1 < NUMDOF_SOH8; ind1++)
        {
          elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
        }
      }
      // this element has to do the whole scaling
      else if (condnum1 <= condnum0)
      {
        for (int i = 0; i < NUMNOD_SOH8; i++)
        {
          adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
        }
        for (int ind1 = 0; ind1 < NUMDOF_SOH8; ind1++)
        {
          if (ind1 < NUMDOF_SOH8 / 2)
          {
            elemat1(ind1, ind1) += factor3 / adjele(ind1, 0) * cond0.size();
            elemat1(ind1, ind1 + NUMDOF_SOH8 / 2) += factor4 / adjele(ind1, 0) * cond0.size();
          }
          else
          {
            elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
          }
        }
      }
    }
    // surface with nodes 4..7 is to be scaled
    else if (current1)
    {
      // but not by this element
      if (condnum0 > condnum1)
      {
        for (int i = 0; i < NUMNOD_SOH8; i++)
        {
          adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
        }
        for (int ind1 = 0; ind1 < NUMDOF_SOH8 / 2; ind1++)
        {
          elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
        }
      }
      // this element has to do the whole scaling
      else if (condnum0 <= condnum1)
      {
        for (int i = 0; i < NUMNOD_SOH8; i++)
        {
          adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
          adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
        }
        for (int ind1 = 0; ind1 < NUMDOF_SOH8; ind1++)
        {
          if (ind1 >= NUMDOF_SOH8 / 2)
          {
            elemat1(ind1, ind1) += factor3 / adjele(ind1, 0) * cond1.size();
            elemat1(ind1, -NUMDOF_SOH8 / 2 + ind1) += factor4 / adjele(ind1, 0) * cond1.size();
          }
          else
          {
            elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
          }
        }
      }
    }
    else
    {
      for (int i = 0; i < NUMNOD_SOH8; i++)
      {
        adjele(NUMDIM_SOH8 * i + 0, 0) = nodes[i]->NumElement();
        adjele(NUMDIM_SOH8 * i + 1, 0) = nodes[i]->NumElement();
        adjele(NUMDIM_SOH8 * i + 2, 0) = nodes[i]->NumElement();
      }
      for (int ind1 = 0; ind1 < NUMDOF_SOH8; ind1++)
      {
        elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
      }
    }
  }
  else
    FOUR_C_THROW("Chosen STC_SCALING not supported!");
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoSh8Type::Initialize(DRT::Discretization& dis)
{
  // sosh8_gmshplotdis(dis);

  int num_morphed_so_hex8_easmild = 0;
  int num_morphed_so_hex8_easnone = 0;

  // Loop through all elements
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::SoSh8*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_sh8* failed");

    if (!actele->nodes_rearranged_)
    {
      bool altered = false;
      switch (actele->thickdir_)
      {
        // check for automatic definition of thickness direction
        case DRT::ELEMENTS::SoSh8::autoj:
        {
          actele->thickdir_ = actele->sosh8_findthickdir();
          altered = true;
          break;
        }
        // check for enforced definition of thickness direction
        case DRT::ELEMENTS::SoSh8::globx:
        {
          CORE::LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(0) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        case DRT::ELEMENTS::SoSh8::globy:
        {
          CORE::LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(1) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        case DRT::ELEMENTS::SoSh8::globz:
        {
          CORE::LINALG::Matrix<NUMDIM_SOH8, 1> thickdirglo(true);
          thickdirglo(2) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        default:
          break;
      }

      if (altered and (actele->thickdir_ != DRT::ELEMENTS::SoSh8::undefined))
      {
        // special element-dependent input of material parameters
        if (actele->Material()->MaterialType() == CORE::Materials::m_viscoanisotropic)
        {
          MAT::ViscoAnisotropic* visco =
              dynamic_cast<MAT::ViscoAnisotropic*>(actele->Material().get());
          visco->Setup(NUMGPT_SOH8, actele->thickvec_);
          if (actele->thickvec_.size() == 0)
            FOUR_C_THROW("zero size thickness vector for element %d", actele->Id());
        }
      }

      int new_nodeids[NUMNOD_SOH8];

      switch (actele->thickdir_)
      {
        case DRT::ELEMENTS::SoSh8::globx:
        case DRT::ELEMENTS::SoSh8::globy:
        case DRT::ELEMENTS::SoSh8::globz:
        {
          FOUR_C_THROW("This should have been replaced by auto(r|s|t)");
          break;
        }
        case DRT::ELEMENTS::SoSh8::autor:
        case DRT::ELEMENTS::SoSh8::enfor:
        {
          // resorting of nodes,
          // such that previous local r-dir is local t-dir afterwards
          new_nodeids[0] = actele->NodeIds()[7];
          new_nodeids[1] = actele->NodeIds()[4];
          new_nodeids[2] = actele->NodeIds()[0];
          new_nodeids[3] = actele->NodeIds()[3];
          new_nodeids[4] = actele->NodeIds()[6];
          new_nodeids[5] = actele->NodeIds()[5];
          new_nodeids[6] = actele->NodeIds()[1];
          new_nodeids[7] = actele->NodeIds()[2];
          //        actele->sosh8_gmshplotlabeledelement(actele->NodeIds());
          //        actele->sosh8_gmshplotlabeledelement(new_nodeids);
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::SoSh8::autos:
        case DRT::ELEMENTS::SoSh8::enfos:
        {
          // resorting of nodes,
          // such that previous local s-dir is local t-dir afterwards
          new_nodeids[0] = actele->NodeIds()[4];
          new_nodeids[1] = actele->NodeIds()[5];
          new_nodeids[2] = actele->NodeIds()[1];
          new_nodeids[3] = actele->NodeIds()[0];
          new_nodeids[4] = actele->NodeIds()[7];
          new_nodeids[5] = actele->NodeIds()[6];
          new_nodeids[6] = actele->NodeIds()[2];
          new_nodeids[7] = actele->NodeIds()[3];
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::SoSh8::autot:
        case DRT::ELEMENTS::SoSh8::enfot:
        {
          // no resorting necessary
          for (int node = 0; node < 8; ++node)
          {
            new_nodeids[node] = actele->NodeIds()[node];
          }
          actele->SetNodeIds(NUMNOD_SOH8, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::SoSh8::undefined:
        {
          if (actele->eastype_ == DRT::ELEMENTS::SoSh8::soh8_eassosh8)
          {
            // here comes plan B: morph So_sh8 to So_hex8
            actele->soh8_reiniteas(DRT::ELEMENTS::SoHex8::soh8_easmild);
            actele->anstype_ = SoSh8::ansnone;
            actele->InitJacobianMapping();
            num_morphed_so_hex8_easmild++;
          }
          else if (actele->eastype_ == DRT::ELEMENTS::SoSh8::soh8_easnone)
          {
            // here comes plan B: morph So_sh8 to So_hex8
            actele->soh8_reiniteas(DRT::ELEMENTS::SoHex8::soh8_easnone);
            actele->anstype_ = SoSh8::ansnone;
            actele->InitJacobianMapping();
            num_morphed_so_hex8_easnone++;
          }
          else if (actele->eastype_ == DRT::ELEMENTS::SoHex8::soh8_easmild)
          {
            // this might happen in post filter (for morped sosh8->soh8)
            actele->soh8_reiniteas(DRT::ELEMENTS::SoHex8::soh8_easmild);
            actele->anstype_ = SoSh8::ansnone;
            actele->InitJacobianMapping();
          }
          else if (actele->eastype_ == DRT::ELEMENTS::SoHex8::soh8_easnone)
          {
            // this might happen in post filter (for morped sosh8->soh8)
            actele->anstype_ = SoSh8::ansnone;
            actele->InitJacobianMapping();
          }
          else
            FOUR_C_THROW("Undefined EAS type");
          break;
        }
        case DRT::ELEMENTS::SoSh8::none:
          break;
        default:
          FOUR_C_THROW("no thickness direction for So_sh8");
      }
      // actele->sosh8_gmshplotlabeledelement(actele->NodeIds());
    }
  }

  if (num_morphed_so_hex8_easmild > 0)
  {
    std::cout << std::endl
              << num_morphed_so_hex8_easmild
              << " Sosh8-Elements have no clear 'thin' direction and have morphed to So_hex8 with "
                 "eas_mild"
              << std::endl;
  }
  if (num_morphed_so_hex8_easnone > 0)
  {
    std::cout << std::endl
              << num_morphed_so_hex8_easnone
              << " Sosh8-Elements have no clear 'thin' direction and have morphed to So_hex8 with "
                 "eas_none"
              << std::endl;
  }

  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false, false, false);

  // loop again to init Jacobian for Sosh8's
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::SoSh8*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_sh8* failed");
    actele->InitJacobianMapping();
  }

  // **************** debug printout ot gmesh **********************************
  // sosh8_gmshplotdis(dis);

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
