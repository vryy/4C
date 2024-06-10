/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 2


*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_mat_service.hpp"
#include "4C_so3_element_service.hpp"
#include "4C_so3_plast_ssn.hpp"
#include "4C_so3_ssn_plast_fwd.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

using VoigtMapping = Core::LinAlg::Voigt::IndexMappings;

/*----------------------------------------------------------------------*
 | evaluate the element (public)                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::Evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material post_setup() routine has already been called and call it if
  // not
  ensure_material_post_setup(params);

  set_params_interface_ptr(params);

  invalid_ele_data();
  if (distype == Core::FE::CellType::nurbs27) get_nurbs_ele_info(&discretization);

  FOUR_C_ASSERT(kintype_ == Inpar::STR::KinemType::nonlinearTotLag,
      "only geometricallly nonlinear formluation for plasticity!");

  // start with "none"
  Core::Elements::ActionType act = Core::Elements::none;
  if (IsParamsInterface())
  {
    act = params_interface().get_action_type();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    act = Core::Elements::String2ActionType(action);
  }

  // what should the element do
  switch (act)
  {
    //============================================================================
    // linear stiffness
    case Core::Elements::struct_calc_linstiff:
    case Core::Elements::struct_calc_linstiffmass:
    {
      FOUR_C_THROW("linear kinematics version out dated");
      break;
    }

    //============================================================================
    // nonlinear stiffness
    case Core::Elements::struct_calc_internalforce:
    {
      // internal force vector
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat1+2, elevec2+3 are not used anyway

      // need current displacement and residual/incremental displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      if ((disp == Teuchos::null))
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(la[0].lm_.size());
      Core::FE::ExtractMyValues(*disp, mydisp, la[0].lm_);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> myemat(true);

      // initialise the vectors
      // Evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      // velocities for TSI problem
      std::vector<double> myvel(0);
      std::vector<double> mytempnp(0);
      if (tsi_)
      {
        if (discretization.HasState(0, "velocity"))
        {
          // get the velocities
          Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0, "velocity");
          if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          Core::FE::ExtractMyValues(*vel, myvel, la[0].lm_);
        }
        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          Core::FE::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, &myemat, nullptr, &elevec1, nullptr, nullptr, params,
          Inpar::STR::stress_none, Inpar::STR::strain_none);


      break;
    }

    //============================================================================
    // nonlinear stiffness
    case Core::Elements::struct_calc_nlnstiff:
    {
      // stiffness
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat1(
          elemat1_epetra.values(), true);
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);

      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(la[0].lm_.size());
      Core::FE::ExtractMyValues(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // Evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      // velocities for TSI problem
      std::vector<double> myvel(0);
      std::vector<double> mytempnp(0);
      if (tsi_)
      {
        if (discretization.HasState(0, "velocity"))
        {
          // get the velocities
          Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0, "velocity");
          if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          Core::FE::ExtractMyValues(*vel, myvel, la[0].lm_);
        }

        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          Core::FE::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, matptr, nullptr, &elevec1, nullptr, nullptr, params,
          Inpar::STR::stress_none, Inpar::STR::strain_none);

      break;
    }  // calc_struct_nlnstiff

    //============================================================================
    // (non)linear stiffness, mass matrix and internal force vector
    case Core::Elements::struct_calc_nlnstiffmass:
    case Core::Elements::struct_calc_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(la[0].lm_.size());
      Core::FE::ExtractMyValues(*disp, mydisp, la[0].lm_);
      // stiffness
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat1(
          elemat1_epetra.values(), true);
      // mass
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat2(
          elemat2_epetra.values(), true);
      // internal force
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);

      // initialise the vectors
      // Evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      // velocities for TSI problem
      std::vector<double> myvel(0);
      std::vector<double> mytempnp(0);
      if (tsi_)
      {
        if (discretization.HasState(0, "velocity"))
        {
          // get the velocities
          Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0, "velocity");
          if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          Core::FE::ExtractMyValues(*vel, myvel, la[0].lm_);
        }
        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          Core::FE::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, &elemat1, &elemat2, &elevec1, nullptr, nullptr, params,
          Inpar::STR::stress_none, Inpar::STR::strain_none);

      if (act == Core::Elements::struct_calc_nlnstifflmass)
        // lump mass matrix
        // we assume #elemat2 is a square matrix
        for (int c = 0; c < elemat2_epetra.numCols(); ++c)  // parse columns
        {
          double d = 0.0;
          for (int r = 0; r < elemat2_epetra.numRows(); ++r)  // parse rows
          {
            d += elemat2(r, c);  // accumulate row entries
            elemat2(r, c) = 0.0;
          }
          elemat2(c, c) = d;  // apply sum of row entries on diagonal
        }
      break;
    }  // calc_struct_nlnstiff(l)mass

    case Core::Elements::struct_calc_stress:
    {
      // elemat1+2,elevec1-3 are not used anyway
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
      Teuchos::RCP<std::vector<char>> straindata = Teuchos::null;
      Inpar::STR::StressType iostress = Inpar::STR::stress_none;
      Inpar::STR::StrainType iostrain = Inpar::STR::strain_none;
      if (IsParamsInterface())
      {
        stressdata = str_params_interface().stress_data_ptr();
        straindata = str_params_interface().strain_data_ptr();

        iostress = str_params_interface().get_stress_output_type();
        iostrain = str_params_interface().get_strain_output_type();
      }
      else
      {
        stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        straindata = params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        iostress = Core::UTILS::GetAsEnum<Inpar::STR::StressType>(
            params, "iostress", Inpar::STR::stress_none);
        iostrain = Core::UTILS::GetAsEnum<Inpar::STR::StrainType>(
            params, "iostrain", Inpar::STR::strain_none);
      }
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
      if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");

      std::vector<double> mydisp((la[0].lm_).size());
      Core::FE::ExtractMyValues(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // Evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      // velocities for TSI problem
      std::vector<double> myvel(0);
      std::vector<double> mytempnp(0);
      std::vector<double> mytempres(0);
      if (tsi_)
      {
        if (discretization.HasState(0, "velocity"))
        {
          // get the velocities
          Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(0, "velocity");
          if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          Core::FE::ExtractMyValues(*vel, myvel, la[0].lm_);
        }
        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          Core::FE::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      Core::LinAlg::Matrix<numgpt_post, numstr_> stress;
      Core::LinAlg::Matrix<numgpt_post, numstr_> strain;

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, nullptr, nullptr, nullptr, &stress, &strain, params,
          iostress, iostrain);
      {
        Core::Communication::PackBuffer data;
        this->AddtoPack(data, stress);
        data.StartPacking();
        this->AddtoPack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }
      {
        Core::Communication::PackBuffer data;
        this->AddtoPack(data, strain);
        data.StartPacking();
        this->AddtoPack(data, strain);
        std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
      }

      break;
    }


    //============================================================================
    // required for predictor TangDis --> can be helpful in compressible case!
    case Core::Elements::struct_calc_reset_istep:
    {
      for (int i = 0; i < numgpt_; i++)
      {
        KbbInv_[i].putScalar(0.0);
        Kbd_[i].putScalar(0.0);
        fbeta_[i].putScalar(0.0);
      }
      break;
    }

    //============================================================================
    case Core::Elements::struct_calc_update_istep:
    {
      // update plastic deformation
      // default: geometrically non-linear analysis with Total Lagrangean approach
      update_plastic_deformation_nln(plspintype_);
      break;
    }  // calc_struct_update_istep

    //============================================================================
    case Core::Elements::struct_calc_stifftemp:  // here we want to build the K_dT block for
                                                 // monolithic TSI
    {
      // stiffness
      Core::LinAlg::Matrix<numdofperelement_, nen_> k_dT(elemat1_epetra.values(), true);

      // calculate matrix block
      nln_kd_t_tsi(&k_dT, params);
    }
    break;

    case Core::Elements::struct_calc_energy:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      std::vector<double> mydisp(la[0].lm_.size());
      Core::FE::ExtractMyValues(*disp, mydisp, la[0].lm_);

      std::vector<double> mytempnp(0);
      if (discretization.HasState(1, "temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
        if (tempnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the dimension of the
        // problem
        const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
        if (la[1].Size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
        Core::FE::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
      }

      double intenergy = calc_int_energy(mydisp, mytempnp, params);

      if (IsParamsInterface())  // new structural time integration
      {
        str_params_interface().add_contribution_to_energy_type(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        FOUR_C_ASSERT(elevec1_epetra.length() < 1, "The given result vector is too short.");
        elevec1_epetra(0) = intenergy;
      }
    }
    break;

    case Core::Elements::struct_calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> myres((la[0].lm_).size());
      Core::FE::ExtractMyValues(*res, myres, la[0].lm_);
      Core::LinAlg::Matrix<nen_ * nsd_, 1> res_d(myres.data(), true);

      std::vector<double> mytempres(0);
      Core::LinAlg::Matrix<nen_, 1> res_t;
      Core::LinAlg::Matrix<nen_, 1>* res_t_ptr = nullptr;
      if (discretization.NumDofSets() > 1)
        if (discretization.HasState(1, "residual temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempres =
              discretization.GetState(1, "residual temperature");
          if (tempres == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'tempres'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempres.resize(((la[0].lm_).size()) / nsd_, 0.0);
          Core::FE::ExtractMyValues(*tempres, mytempres, la[1].lm_);
          res_t = Core::LinAlg::Matrix<nen_, 1>(mytempres.data(), true);
          res_t_ptr = &res_t;
        }

      recover_plasticity_and_eas(&res_d, res_t_ptr);
    }
    break;

    case Core::Elements::struct_calc_predict:
    {
      switch (str_params_interface().get_predictor_type())
      {
        case Inpar::STR::pred_constdis:
        default:
          for (int gp = 0; gp < numgpt_; ++gp) dDp_last_iter_[gp].putScalar(0.0);
          if (eastype_ != soh8p_easnone)
            for (int i = 0; i < neas_; i++) (*alpha_eas_)(i) = (*alpha_eas_last_timestep_)(i);
          break;
        case Inpar::STR::pred_constvel:
        case Inpar::STR::pred_tangdis:
          // do nothing for plasticity
          if (eastype_ != soh8p_easnone)
            for (int i = 0; i < neas_; i++)
              (*alpha_eas_)(i) =
                  0. + (*alpha_eas_last_timestep_)(i) + (*alpha_eas_delta_over_last_timestep_)(i);
          break;
      }
    }
    break;

    case Core::Elements::struct_init_gauss_point_data_output:
    {
      FOUR_C_ASSERT(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");
      // Save number of Gauss of the element for gauss point data output
      str_params_interface()
          .gauss_point_data_output_manager_ptr()
          ->add_element_number_of_gauss_points(numgpt_post);
      // holder for output quantity names and their size
      std::unordered_map<std::string, int> quantities_map{};
      // Ask material for the output quantity names and sizes
      SolidMaterial()->register_output_data_names(quantities_map);
      // Add quantities to the Gauss point output data manager (if they do not already exist)
      str_params_interface().gauss_point_data_output_manager_ptr()->merge_quantities(
          quantities_map);
    }
    break;

    case Core::Elements::struct_gauss_point_data_output:
    {
      FOUR_C_ASSERT(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Collection and assembly of gauss point data
      for (const auto& [quantity_name, quantity_size] :
          str_params_interface().gauss_point_data_output_manager_ptr()->get_quantities())
      {
        // Step 1: Collect the data for each Gauss point for the material
        Core::LinAlg::SerialDenseMatrix gp_data(numgpt_post, quantity_size, true);
        bool data_available = SolidMaterial()->EvaluateOutputData(quantity_name, gp_data);

        // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
        // point)
        if (data_available)
        {
          switch (str_params_interface().gauss_point_data_output_manager_ptr()->get_output_type())
          {
            case Inpar::STR::GaussPointDataOutputType::element_center:
            {
              // compute average of the quantities
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  str_params_interface()
                      .gauss_point_data_output_manager_ptr()
                      ->get_element_center_data()
                      .at(quantity_name);
              Core::FE::AssembleAveragedElementValues<Core::LinAlg::SerialDenseMatrix>(
                  *global_data, gp_data, *this);
              break;
            }
            case Inpar::STR::GaussPointDataOutputType::nodes:
            {
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  str_params_interface().gauss_point_data_output_manager_ptr()->get_nodal_data().at(
                      quantity_name);

              Epetra_IntVector& global_nodal_element_count =
                  *str_params_interface()
                       .gauss_point_data_output_manager_ptr()
                       ->get_nodal_data_count()
                       .at(quantity_name);

              if (distype == Core::FE::CellType::hex8)
              {
                if (quantity_size == 1)
                {
                  Core::LinAlg::Matrix<numgpt_post, 1> gp_data_view(gp_data, true);
                  soh8_expol(gp_data_view, *global_data);
                }
                else if (quantity_size == 9)
                {
                  Core::LinAlg::Matrix<numgpt_post, 9> gp_data_view(gp_data, true);
                  soh8_expol(gp_data_view, *global_data);
                }
              }
              else
                FOUR_C_THROW(
                    "only element centered and Gauss point material internal variables for "
                    "so3_plast");

              Discret::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, *this);
              break;
            }
            case Inpar::STR::GaussPointDataOutputType::gauss_points:
            {
              std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                  str_params_interface()
                      .gauss_point_data_output_manager_ptr()
                      ->get_gauss_point_data()
                      .at(quantity_name);
              Discret::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, *this);
              break;
            }
            case Inpar::STR::GaussPointDataOutputType::none:
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

    //============================================================================
    default:
      FOUR_C_THROW("Unknown type of action for So3Plast");
      break;
  }  // action

  return 0;
}  // Evaluate()



/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::calculate_bop(
    Core::LinAlg::Matrix<numstr_, numdofperelement_>* bop,
    const Core::LinAlg::Matrix<nsd_, nsd_>* defgrd, const Core::LinAlg::Matrix<nsd_, nen_>* N_XYZ,
    const int gp)
{
  // lump mass matrix
  if (bop != nullptr)
  {
    /* non-linear B-operator (may so be called, meaning of B-operator is not so
    **  sharp in the non-linear realm) *
    **   B = F^{i,T} . B_L *
    ** with linear B-operator B_L =  N_XYZ (6x24) = (3x8)
    **
    **   B    =   F  . N_XYZ
    ** (6x24)   (3x3) (3x8)
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
    for (int i = 0; i < nen_; ++i)
    {
      (*bop)(0, numdofpernode_ * i + 0) = (*defgrd)(0, 0) * (*N_XYZ)(0, i);
      (*bop)(0, numdofpernode_ * i + 1) = (*defgrd)(1, 0) * (*N_XYZ)(0, i);
      (*bop)(0, numdofpernode_ * i + 2) = (*defgrd)(2, 0) * (*N_XYZ)(0, i);
      (*bop)(1, numdofpernode_ * i + 0) = (*defgrd)(0, 1) * (*N_XYZ)(1, i);
      (*bop)(1, numdofpernode_ * i + 1) = (*defgrd)(1, 1) * (*N_XYZ)(1, i);
      (*bop)(1, numdofpernode_ * i + 2) = (*defgrd)(2, 1) * (*N_XYZ)(1, i);
      (*bop)(2, numdofpernode_ * i + 0) = (*defgrd)(0, 2) * (*N_XYZ)(2, i);
      (*bop)(2, numdofpernode_ * i + 1) = (*defgrd)(1, 2) * (*N_XYZ)(2, i);
      (*bop)(2, numdofpernode_ * i + 2) = (*defgrd)(2, 2) * (*N_XYZ)(2, i);
      /* ~~~ */
      (*bop)(3, numdofpernode_ * i + 0) =
          (*defgrd)(0, 0) * (*N_XYZ)(1, i) + (*defgrd)(0, 1) * (*N_XYZ)(0, i);
      (*bop)(3, numdofpernode_ * i + 1) =
          (*defgrd)(1, 0) * (*N_XYZ)(1, i) + (*defgrd)(1, 1) * (*N_XYZ)(0, i);
      (*bop)(3, numdofpernode_ * i + 2) =
          (*defgrd)(2, 0) * (*N_XYZ)(1, i) + (*defgrd)(2, 1) * (*N_XYZ)(0, i);
      (*bop)(4, numdofpernode_ * i + 0) =
          (*defgrd)(0, 1) * (*N_XYZ)(2, i) + (*defgrd)(0, 2) * (*N_XYZ)(1, i);
      (*bop)(4, numdofpernode_ * i + 1) =
          (*defgrd)(1, 1) * (*N_XYZ)(2, i) + (*defgrd)(1, 2) * (*N_XYZ)(1, i);
      (*bop)(4, numdofpernode_ * i + 2) =
          (*defgrd)(2, 1) * (*N_XYZ)(2, i) + (*defgrd)(2, 2) * (*N_XYZ)(1, i);
      (*bop)(5, numdofpernode_ * i + 0) =
          (*defgrd)(0, 2) * (*N_XYZ)(0, i) + (*defgrd)(0, 0) * (*N_XYZ)(2, i);
      (*bop)(5, numdofpernode_ * i + 1) =
          (*defgrd)(1, 2) * (*N_XYZ)(0, i) + (*defgrd)(1, 0) * (*N_XYZ)(2, i);
      (*bop)(5, numdofpernode_ * i + 2) =
          (*defgrd)(2, 2) * (*N_XYZ)(0, i) + (*defgrd)(2, 0) * (*N_XYZ)(2, i);
    }
  }
}  // calculate_bop()



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 04/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::So3Plast<distype>::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = condition.parameters().GetIf<std::vector<int>>("onoff");
  const auto* val = condition.parameters().GetIf<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = params_interface().get_total_time();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < nsd_)
    FOUR_C_THROW("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = nsd_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = condition.parameters().GetIf<std::vector<int>>("funct");
  Core::LinAlg::Matrix<nsd_, 1> xrefegp(false);
  const bool havefunct = funct != nullptr && std::any_of(funct->begin(), funct->end(),
                                                 [](const int i) { return i > 0; });

  // update element geometry
  Core::LinAlg::Matrix<nen_, nsd_> xrefe;  // material coord. of element
  Core::Nodes::Node** nodes = Nodes();
  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::LinAlg::Matrix<nen_, 1> shapefunct;
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    Core::LinAlg::Matrix<nsd_, nen_> deriv;
    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);

    // compute the Jacobian matrix
    Core::LinAlg::Matrix<nsd_, nsd_> jac;
    jac.Multiply(deriv, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < nsd_; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < nen_; ++nodid)
          xrefegp(dim) += shapefunct(nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = wgt_[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < nsd_; dim++)
    {
      // function evaluation
      const int functnum = (funct) ? (*funct)[dim] : -1;
      const double functfac =
          (functnum > 0) ? Global::Problem::Instance()
                               ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                               .Evaluate(xrefegp.A(), time, dim)
                         : 1.0;
      const double dim_fac = (*onoff)[dim] * (*val)[dim] * fac * functfac;
      for (int nodid = 0; nodid < nen_; ++nodid)
      {
        elevec1[nodid * nsd_ + dim] += shapefunct(nodid) * dim_fac;
      }
    }

  } /* ==================================================== end of Loop over GP */

  return 0;
}

/*----------------------------------------------------------------------*
 | initialise Jacobian                                      seitz 07/13 |
 | is called once in Initialize() in so3_ssn_plast_eletypes.cpp         |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::init_jacobian_mapping()
{
  return;
}  // init_jacobian_mapping()}


/*----------------------------------------------------------------------*
 | internal force, stiffness and mass for f-bar elements    seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::nln_stiffmass(
    std::vector<double>& disp,         // current displacements
    std::vector<double>& vel,          // current velocities
    std::vector<double>& temperature,  // current temperatures
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*
        stiffmatrix,  // element stiffness matrix
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>* massmatrix,  // element mass matrix
    Core::LinAlg::Matrix<numdofperelement_, 1>* force,      // element internal force vector
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,                         // algorithmic parameters e.g. time
    const Inpar::STR::StressType iostress,                  // stress output option
    const Inpar::STR::StrainType iostrain                   // strain output option
)
{
  const bool eval_tsi = (temperature.size() != 0);
  const bool is_tangDis = str_params_interface().get_predictor_type() == Inpar::STR::pred_tangdis;
  const double dt = str_params_interface().get_delta_time();
  const double timefac_d = str_params_interface().get_tim_int_factor_vel() /
                           str_params_interface().get_tim_int_factor_disp();
  if (timefac_d <= 0 || dt <= 0) FOUR_C_THROW("time integration parameters not provided");

  fill_position_arrays(disp, vel, temperature);

  if (fbar_ || eastype_ != soh8p_easnone)
    evaluate_center();  // deformation gradient at centroid of element

  if (eastype_ != soh8p_easnone) eas_setup();

  // get plastic hyperelastic material
  Mat::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());
  else
  {
    if (tsi_) FOUR_C_THROW("TSI with so3Plast elements only with PlasticElastHyper material");
  }

  // EAS matrix block
  Core::LinAlg::SerialDenseMatrix Kda(numdofperelement_, neas_);
  std::vector<Core::LinAlg::SerialDenseVector> dHda(0);
  if (eastype_ != soh8p_easnone && eval_tsi)
    dHda.resize(numgpt_, Core::LinAlg::SerialDenseVector(neas_));
  // temporary Epetra matrix for this and that
  Core::LinAlg::SerialDenseMatrix tmp;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    invalid_gp_data();

    // shape functions (shapefunct) and their first derivatives (deriv)
    evaluate_shape(xsi_[gp]);
    evaluate_shape_deriv(xsi_[gp]);

    kinematics(gp);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8p_easnone)
    {
      eas_shape(gp);
      eas_enhance_strains();
    }
    // calculate modified deformation gradient
    else if (fbar_)
      setup_fbar_gp();
    else
      set_defgrd_mod() = defgrd();

    output_strains(gp, iostrain, elestrain);

    // Gauss point temperature
    double gp_temp = -1.e12;
    if (eval_tsi) gp_temp = temp().Dot(shape_function());

    // material call *********************************************
    if (plmat != nullptr)
    {
      plmat->EvaluateElast(&defgrd_mod(), &delta_lp(), &set_p_k2(), &set_cmat(), gp, Id());
      if (eval_tsi)
        plmat->evaluate_thermal_stress(&defgrd_mod(), gp_temp, &set_p_k2(), &set_cmat(), gp, Id());
    }
    else
    {
      static Core::LinAlg::Matrix<numstr_, 1> total_glstrain(false);
      total_glstrain(0) = 0.5 * (rcg()(0, 0) - 1.0);
      total_glstrain(1) = 0.5 * (rcg()(1, 1) - 1.0);
      total_glstrain(2) = 0.5 * (rcg()(2, 2) - 1.0);
      total_glstrain(3) = rcg()(0, 1);
      total_glstrain(4) = rcg()(1, 2);
      total_glstrain(5) = rcg()(2, 0);

      SolidMaterial()->Evaluate(
          &defgrd_mod(), &total_glstrain, params, &set_p_k2(), &set_cmat(), gp, Id());
    }
    // material call *********************************************

    output_stress(gp, iostress, elestress);

    // integrate usual internal force and stiffness matrix
    double detJ_w = det_j() * wgt_[gp];
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != nullptr) integrate_force(gp, *force);

    // update stiffness matrix
    if (stiffmatrix != nullptr) integrate_stiff_matrix(gp, *stiffmatrix, Kda);

    if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
      integrate_mass_matrix(gp, *massmatrix);

    if (eval_tsi && (stiffmatrix != nullptr || force != nullptr) && plmat != nullptr)
      integrate_thermo_gp(gp, dHda[gp]);

    // plastic modifications
    if (plmat != nullptr)
      if (!plmat->AllElastic())
        if ((stiffmatrix != nullptr || force != nullptr) && !is_tangDis)
        {
          if (have_plastic_spin())
          {
            if (fbar_)
              condense_plasticity<plspin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), &rc_gvec(), detJ_w, gp, gp_temp, params, force,
                  stiffmatrix, nullptr, nullptr, nullptr, &fbar_fac(), &htensor());
            else if (eastype_ != soh8p_easnone)
              condense_plasticity<plspin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), nullptr, detJ_w, gp, gp_temp, params, force,
                  stiffmatrix, &m_eas(), &Kda, &dHda);
            else
              condense_plasticity<plspin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), nullptr, detJ_w, gp, gp_temp, params, force,
                  stiffmatrix);
          }
          else
          {
            if (fbar_)
              condense_plasticity<zerospin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), &rc_gvec(), detJ_w, gp, gp_temp, params, force,
                  stiffmatrix, nullptr, nullptr, nullptr, &fbar_fac(), &htensor());
            else if (eastype_ != soh8p_easnone)
              condense_plasticity<zerospin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), nullptr, detJ_w, gp, gp_temp, params, force,
                  stiffmatrix, &m_eas(), &Kda, &dHda);
            else
              condense_plasticity<zerospin>(defgrd_mod(), delta_lp(), bop(),
                  &deriv_shape_function_xyz(), nullptr, detJ_w, gp, gp_temp, params, force,
                  stiffmatrix);
          }
        }  // plastic modifications
  }        // gp loop
  invalid_gp_data();

  // Static condensation EAS --> stiff ********************************
  if (stiffmatrix != nullptr && !is_tangDis && eastype_ != soh8p_easnone)
  {
    using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
    using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_inverseKaa;
    solve_for_inverseKaa.setMatrix(KaaInv_);
    solve_for_inverseKaa.invert();

    Core::LinAlg::SerialDenseMatrix kdakaai(numdofperelement_, neas_);
    switch (eastype_)
    {
      case soh8p_easfull:
        Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
            0., kdakaai.values(), 1., Kda.values(), KaaInv_->values());
        if (stiffmatrix != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.values(), Kad_->values());
        if (force != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., force->A(), -1., kdakaai.values(), feas_->values());
        break;
      case soh8p_easmild:
        Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
            0., kdakaai.values(), 1., Kda.values(), KaaInv_->values());
        if (stiffmatrix != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.values(), Kad_->values());
        if (force != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., force->A(), -1., kdakaai.values(), feas_->values());
        break;
      case soh8p_easnone:
        break;
      default:
        FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
        break;
    }

    // TSI with EAS
    if (eval_tsi)
    {
      Core::LinAlg::SerialDenseVector dHdaKaai(neas_);
      switch (eastype_)
      {
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, nen_>(
              0., KdT_eas_->A(), -1., kdakaai.values(), KaT_->values());
          for (int gp = 0; gp < numgpt_; ++gp)
          {
            Core::LinAlg::DenseFunctions::multiply<double, 1,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
                0., dHdaKaai.values(), 1., dHda.at(gp).values(), KaaInv_->values());
            Core::LinAlg::DenseFunctions::multiplyTN<double, numdofperelement_,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
                1., plmat->dHepDissDd(gp).values(), -1., Kad_->values(), dHdaKaai.values());
            Core::LinAlg::DenseFunctions::multiply<double, 1,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
                1., &(plmat->HepDiss(gp)), -1., dHdaKaai.values(), feas_->values());
            Core::LinAlg::DenseFunctions::multiplyTN<double, nen_,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
                0., plmat->dHepDTeas()->at(gp).values(), -1., KaT_->values(), dHdaKaai.values());
          }
          break;
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, nen_>(
              0., KdT_eas_->A(), -1., kdakaai.values(), KaT_->values());
          for (int gp = 0; gp < numgpt_; ++gp)
          {
            Core::LinAlg::DenseFunctions::multiply<double, 1,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
                0., dHdaKaai.values(), 1., dHda.at(gp).values(), KaaInv_->values());
            Core::LinAlg::DenseFunctions::multiplyTN<double, numdofperelement_,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
                1., plmat->dHepDissDd(gp).values(), -1., Kad_->values(), dHdaKaai.values());
            Core::LinAlg::DenseFunctions::multiply<double, 1,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
                1., &(plmat->HepDiss(gp)), -1., dHdaKaai.values(), feas_->values());
            Core::LinAlg::DenseFunctions::multiplyTN<double, nen_,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
                0., plmat->dHepDTeas()->at(gp).values(), -1., KaT_->values(), dHdaKaai.values());
          }
          break;
        case soh8p_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }
  }

  invalid_ele_data();
  return;
}

/*----------------------------------------------------------------------*
 | coupling term k_dT in monolithic TSI                     seitz 06/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::nln_kd_t_tsi(
    Core::LinAlg::Matrix<numdofperelement_, nen_>* k_dT, Teuchos::ParameterList& params)
{
  if (k_dT == nullptr) return;

  // shape functions
  Core::LinAlg::Matrix<nen_, 1> shapefunct(false);

  for (int gp = 0; gp < numgpt_; gp++)
  {
    // shape functions
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    // update linear coupling matrix K_dT
    k_dT->MultiplyNT(1., (*dFintdT_)[gp], shapefunct, 1.);
  }

  // EAS part
  if (eastype_ != soh8p_easnone)
  {
    if (KdT_eas_ == Teuchos::null)
      FOUR_C_THROW("for TSI with EAS the block KdT_eas_ should be acessible here");
    k_dT->Update(1., *KdT_eas_, 1.);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  condense plastic degrees of freedom                     seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
template <int spintype>
void Discret::ELEMENTS::So3Plast<distype>::condense_plasticity(
    const Core::LinAlg::Matrix<nsd_, nsd_>& defgrd, const Core::LinAlg::Matrix<nsd_, nsd_>& deltaLp,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>& bop,
    const Core::LinAlg::Matrix<nsd_, nen_>* N_XYZ, const Core::LinAlg::Matrix<numstr_, 1>* RCG,
    const double detJ_w, const int gp, const double temp, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<numdofperelement_, 1>* force,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>* stiffmatrix,
    const Core::LinAlg::SerialDenseMatrix* M, Core::LinAlg::SerialDenseMatrix* Kda,
    std::vector<Core::LinAlg::SerialDenseVector>* dHda, const double* f_bar_factor,
    const Core::LinAlg::Matrix<numdofperelement_, 1>* htensor)
{
  bool eval_tsi = tsi_ && (temp != -1.e12);

  // get plastic hyperelastic material
  Mat::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());
  else
    FOUR_C_THROW("so3_ssn_plast elements only with PlasticElastHyper material");

  // temporary Epetra matrix for matrix-matrix-matrix products
  Core::LinAlg::SerialDenseMatrix tmp;

  // Nitsche contact
  Core::LinAlg::Matrix<numstr_, 1>* cauchy_ptr = nullptr;
  Core::LinAlg::Matrix<numstr_, spintype + 1> d_cauchy_ddp;
  Core::LinAlg::Matrix<numstr_, spintype + 1>* d_cauchy_ddp_ptr = nullptr;
  Core::LinAlg::Matrix<numstr_, numstr_> d_cauchy_dC;
  Core::LinAlg::Matrix<numstr_, numstr_>* d_cauchy_dC_ptr = nullptr;
  Core::LinAlg::Matrix<numstr_, 9> d_cauchy_dF;
  Core::LinAlg::Matrix<numstr_, 9>* d_cauchy_dF_ptr = nullptr;
  Core::LinAlg::Matrix<numstr_, 1>
      d_cauchy_dT;  // todo: continue with this one; plmat->EvaluatePlast(...) should give you this;
  // however not tested yet
  Core::LinAlg::Matrix<numstr_, 1>* d_cauchy_dT_ptr = nullptr;
  if (is_nitsche_contact_)
  {
    cauchy_ptr = &(cauchy_.at(gp));
    d_cauchy_ddp_ptr = &d_cauchy_ddp;
    d_cauchy_dC_ptr = &d_cauchy_dC;
    d_cauchy_dF_ptr = &d_cauchy_dF;
    if (eval_tsi) d_cauchy_dT_ptr = &d_cauchy_dT;
  }

  // second material call ****************************************************
  Core::LinAlg::Matrix<numstr_, spintype + 1> dpk2ddp;
  Core::LinAlg::Matrix<spintype + 1, 1> ncp;
  Core::LinAlg::Matrix<spintype + 1, spintype + 1> dncpddp;
  Core::LinAlg::Matrix<spintype + 1, numstr_> dncpdc;
  Core::LinAlg::Matrix<spintype + 1, 1> dncpdT;
  Core::LinAlg::Matrix<numstr_, 1> dHdC;
  Core::LinAlg::Matrix<spintype + 1, 1> dHdLp;
  bool active = false;
  bool elast = false;
  bool as_converged = true;
  if (!eval_tsi)
    plmat->EvaluatePlast(&defgrd, &deltaLp, nullptr, params, &dpk2ddp, &ncp, &dncpdc, &dncpddp,
        &active, &elast, &as_converged, gp, nullptr, nullptr, nullptr,
        str_params_interface().get_delta_time(), Id(), cauchy_ptr, d_cauchy_ddp_ptr,
        d_cauchy_dC_ptr, d_cauchy_dF_ptr, d_cauchy_dT_ptr);
  else
    plmat->EvaluatePlast(&defgrd, &deltaLp, &temp, params, &dpk2ddp, &ncp, &dncpdc, &dncpddp,
        &active, &elast, &as_converged, gp, &dncpdT, &dHdC, &dHdLp,
        str_params_interface().get_delta_time(), Id(), cauchy_ptr, d_cauchy_ddp_ptr,
        d_cauchy_dC_ptr, d_cauchy_dF_ptr, d_cauchy_dT_ptr);
  // *************************************************************************

  // Simple matrix do delete the linear dependent row in Voigt-notation.
  // Since all entries are deviatoric, we drop the third entry on the main diagonal.
  Core::LinAlg::Matrix<spintype + 1, spintype> voigt_red(true);
  switch (spintype)
  {
    case zerospin:
      voigt_red(0, 0) = voigt_red(1, 1) = voigt_red(3, 2) = voigt_red(4, 3) = voigt_red(5, 4) = 1.;
      break;
    case plspin:
      voigt_red(0, 0) = voigt_red(1, 1) = voigt_red(3, 2) = voigt_red(4, 3) = voigt_red(5, 4) =
          voigt_red(6, 5) = voigt_red(7, 6) = voigt_red(8, 7) = 1.;
      break;
    default:
      FOUR_C_THROW("Don't know what to do with plastic spin type %d", spintype);
      break;
  }

  // We have to adapt the derivatives to account for the deviatoric structure
  // of delta D^p.
  Core::LinAlg::Matrix<spintype + 1, spintype> dDpdbeta;
  switch (spintype)
  {
    case zerospin:
      dDpdbeta(0, 0) = dDpdbeta(1, 1) = dDpdbeta(5, 4) = dDpdbeta(4, 3) = dDpdbeta(3, 2) = 1.;
      dDpdbeta(2, 0) = dDpdbeta(2, 1) = -1.;
      break;
    case plspin:
      dDpdbeta(0, 0) = dDpdbeta(1, 1) = +1.;
      dDpdbeta(2, 0) = dDpdbeta(2, 1) = -1.;
      dDpdbeta(3, 2) = 1.;
      dDpdbeta(4, 3) = 1.;
      dDpdbeta(5, 4) = 1.;
      dDpdbeta(6, 2) = 1.;
      dDpdbeta(7, 3) = 1.;
      dDpdbeta(8, 4) = 1.;

      dDpdbeta(3, 5) = 1.;
      dDpdbeta(4, 6) = 1.;
      dDpdbeta(5, 7) = 1.;
      dDpdbeta(6, 5) = -1.;
      dDpdbeta(7, 6) = -1.;
      dDpdbeta(8, 7) = -1.;
      break;
    default:
      FOUR_C_THROW("Don't know what to do with plastic spin type %d", spintype);
      break;
  }

  // derivative of internal force w.r.t. beta
  Core::LinAlg::Matrix<numdofperelement_, spintype> kdbeta;
  // derivative of pk2 w.r.t. beta
  Core::LinAlg::Matrix<numstr_, spintype> dpk2db;
  dpk2db.Multiply(dpk2ddp, dDpdbeta);
  if (fbar_)
    kdbeta.MultiplyTN(detJ_w / (*f_bar_factor), bop, dpk2db);
  else
    kdbeta.MultiplyTN(detJ_w, bop, dpk2db);

  Core::LinAlg::Matrix<numstr_, spintype> d_cauchy_db;
  if (is_nitsche_contact_)
  {
    if (N_XYZ == nullptr) FOUR_C_THROW("shape derivative not provided");

    Core::LinAlg::Matrix<9, numdofperelement_> dFdd;
    for (int k = 0; k < nen_; ++k)
    {
      for (int d = 0; d < 3; ++d) dFdd(d, k * nsd_ + d) = (*N_XYZ)(d, k);
      dFdd(3, k * nsd_ + 0) = (*N_XYZ)(1, k);
      dFdd(4, k * nsd_ + 1) = (*N_XYZ)(2, k);
      dFdd(5, k * nsd_ + 0) = (*N_XYZ)(2, k);

      dFdd(6, k * nsd_ + 1) = (*N_XYZ)(0, k);
      dFdd(7, k * nsd_ + 2) = (*N_XYZ)(1, k);
      dFdd(8, k * nsd_ + 2) = (*N_XYZ)(0, k);
    }
    cauchy_deriv_.at(gp).Clear();
    if (fbar_)
    {
      if (RCG == nullptr) FOUR_C_THROW("condense_plasticity(...) requires RCG in case of FBAR");
      Core::LinAlg::Matrix<6, 1> tmp61;
      tmp61.Multiply(.5, d_cauchy_dC, (*RCG));
      cauchy_deriv_.at(gp).MultiplyNT(
          (*f_bar_factor) * (*f_bar_factor) * 2. / 3., tmp61, *htensor, 1.);
      cauchy_deriv_.at(gp).Multiply((*f_bar_factor) * (*f_bar_factor), d_cauchy_dC, bop, 1.);

      Core::LinAlg::Matrix<9, 1> f;
      for (int i = 0; i < 3; ++i) f(i) = defgrd(i, i);
      f(3) = defgrd(0, 1);
      f(4) = defgrd(1, 2);
      f(5) = defgrd(0, 2);
      f(6) = defgrd(1, 0);
      f(7) = defgrd(2, 1);
      f(8) = defgrd(2, 0);

      tmp61.Multiply(1., d_cauchy_dF, f, 0.);
      cauchy_deriv_.at(gp).MultiplyNT(*f_bar_factor / 3., tmp61, *htensor, 1.);
      cauchy_deriv_.at(gp).Multiply(*f_bar_factor, d_cauchy_dF, dFdd, 1.);
    }
    else
    {
      cauchy_deriv_.at(gp).Multiply(1., d_cauchy_dC, bop, 1.);
      cauchy_deriv_.at(gp).Multiply(1., d_cauchy_dF, dFdd, 1.);
    }
    d_cauchy_db.Multiply(d_cauchy_ddp, dDpdbeta);
  }

  if (d_cauchy_dT_ptr)
  {
    cauchy_deriv_T_.at(gp).Clear();
    cauchy_deriv_T_.at(gp).MultiplyNT(1., d_cauchy_dT, shape_function(), 1.);
  }

  // EAS matrix block
  Core::LinAlg::SerialDenseMatrix Kab(neas_, spintype);
  switch (eastype_)
  {
    case soh8p_easnone:
      // do nothing
      break;
    case soh8p_easmild:
      Core::LinAlg::DenseFunctions::multiplyTN<double,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, spintype>(
          1., Kab.values(), detJ_w, M->values(), dpk2db.A());
      break;
    case soh8p_easfull:
      Core::LinAlg::DenseFunctions::multiplyTN<double,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, spintype>(
          1., Kab.values(), detJ_w, M->values(), dpk2db.A());
      break;
    case soh8p_eassosh8:
      Core::LinAlg::DenseFunctions::multiplyTN<double,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, numstr_, spintype>(
          1., Kab.values(), detJ_w, M->values(), dpk2db.A());
      break;
    case soh18p_eassosh18:
      Core::LinAlg::DenseFunctions::multiplyTN<double,
          PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, numstr_, spintype>(
          1., Kab.values(), detJ_w, M->values(), dpk2db.A());
      break;
    default:
      FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
      break;
  }

  // gauss points that require a full condensation
  if (!elast)
  {
    // apply chain rule for kbb block
    Core::LinAlg::Matrix<spintype + 1, spintype> dNCPdb;
    dNCPdb.Multiply(dncpddp, dDpdbeta);
    Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype + 1, spintype>(
        0., KbbInv_[gp].values(), 1., voigt_red.A(), dNCPdb.A());

    // apply chain rule for kbd block
    Core::LinAlg::Matrix<spintype + 1, numdofperelement_> dNCPdd;
    if (fbar_)
    {
      if (RCG == nullptr) FOUR_C_THROW("condense_plasticity(...) requires RCG in case of FBAR");
      Core::LinAlg::Matrix<spintype + 1, 1> tmp61;
      tmp61.Multiply(.5, dncpdc, (*RCG));
      dNCPdd.MultiplyNT((*f_bar_factor) * (*f_bar_factor) * 2. / 3., tmp61, *htensor, 0.);
      dNCPdd.Multiply((*f_bar_factor) * (*f_bar_factor), dncpdc, bop, 1.);
    }
    else
      dNCPdd.Multiply(dncpdc, bop);
    Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype + 1, numdofperelement_>(
        0., Kbd_[gp].values(), 1., voigt_red.A(), dNCPdd.A());

    // EAS block kba
    if (eastype_ != soh8p_easnone)
    {
      Kba_->at(gp).shape(spintype, neas_);
      tmp.shape(spintype + 1, neas_);
      switch (eastype_)
      {
        case soh8p_easnone:
          break;
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              0., tmp.values(), 1., dncpdc.A(), M->values());
          Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype + 1,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              0., Kba_->at(gp).values(), 1., voigt_red.A(), tmp.values());
          break;
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              0., tmp.values(), 1., dncpdc.A(), M->values());
          Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              0., Kba_->at(gp).values(), 1., voigt_red.A(), tmp.values());
          break;
        case soh8p_eassosh8:
          Core::LinAlg::DenseFunctions::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
              0., tmp.values(), 1., dncpdc.A(), M->values());
          Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
              0., Kba_->at(gp).values(), 1., voigt_red.A(), tmp.values());
          break;
        case soh18p_eassosh18:
          Core::LinAlg::DenseFunctions::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas>(
              0., tmp.values(), 1., dncpdc.A(), M->values());
          Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas>(
              0., Kba_->at(gp).values(), 1., voigt_red.A(), tmp.values());
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    // residual
    Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype + 1, 1>(
        0., fbeta_[gp].values(), 1., voigt_red.A(), ncp.A());

    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
    // inverse matrix block [k_beta beta]_ij
    using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
    using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_kbbinv;
    solve_for_kbbinv.setMatrix(Teuchos::rcpFromRef(KbbInv_[gp]));
    int err = solve_for_kbbinv.invert();
    if (err != 0)
    {
      // check, if errors are tolerated or should throw a FOUR_C_THROW
      bool error_tol = false;
      if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error", true);
        return;
      }
      else
        FOUR_C_THROW("Inversion of Kbb failed");
    }
    // temporary  Kdb.Kbb^-1
    Core::LinAlg::Matrix<numdofperelement_, spintype> KdbKbb;
    Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype, spintype>(
        0., KdbKbb.A(), 1., kdbeta.A(), KbbInv_[gp].values());

    // "plastic displacement stiffness"
    // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
    if (stiffmatrix != nullptr)
      Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype,
          numdofperelement_>(1., stiffmatrix->A(), -1., KdbKbb.A(), Kbd_[gp].values());

    // "plastic internal force"
    // plFint = [K_db.K_bb^-1].f_b
    if (force != nullptr)
      Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype, 1>(
          1., force->A(), -1., KdbKbb.A(), fbeta_[gp].values());

    // TSI
    if (eval_tsi)
    {
      // thermal derivative
      (*KbT_).at(gp).size(spintype);
      Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype + 1, 1>(
          0., (*KbT_)[gp].values(), 1., voigt_red.A(), dncpdT.A());

      // condense to K_dT
      Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype, 1>(
          1., (*dFintdT_)[gp].A(), -1., KdbKbb.A(), (*KbT_)[gp].values());

      // Plastic heating and dissipation dC
      Core::LinAlg::Matrix<numdofperelement_, 1> dHepDissDd;
      if (fbar_)
      {
        if (RCG == nullptr) FOUR_C_THROW("condense_plasticity(...) requires RCG in case of FBAR");
        Core::LinAlg::Matrix<1, 1> tmp11;
        tmp11.MultiplyTN(.5, dHdC, (*RCG));
        dHepDissDd.Multiply((*f_bar_factor) * (*f_bar_factor) * 2. / 3., *htensor, tmp11, 0.);
        dHepDissDd.MultiplyTN((*f_bar_factor) * (*f_bar_factor), bop, dHdC, 1.);
      }
      else
        dHepDissDd.MultiplyTN(bop, dHdC);

      // store in material
      Core::LinAlg::DenseFunctions::update<double, numdofperelement_, 1>(
          1., plmat->dHepDissDd(gp).values(), 1., dHepDissDd.A());

      // Plastic heating and dissipation dbeta
      Core::LinAlg::Matrix<spintype, 1> dHepDissDbeta;
      if (spintype != zerospin)
        FOUR_C_THROW("no TSI with plastic spin yet");
      else
        Core::LinAlg::DenseFunctions::multiplyTN<double, zerospin, zerospin + 1, 1>(
            0., dHepDissDbeta.A(), 1., dDpdbeta.A(), dHdLp.A());

      // condense the heating terms
      Core::LinAlg::Matrix<spintype, 1> dHdbKbbi;
      Core::LinAlg::DenseFunctions::multiplyTN<double, spintype, spintype, 1>(
          0., dHdbKbbi.A(), 1., KbbInv_[gp].values(), dHepDissDbeta.A());
      Core::LinAlg::DenseFunctions::multiplyTN<double, 1, spintype, 1>(
          1., &(plmat->HepDiss(gp)), -1., dHdbKbbi.A(), fbeta_[gp].values());
      Core::LinAlg::DenseFunctions::multiplyTN<double, numdofperelement_, spintype, 1>(
          1., plmat->dHepDissDd(gp).values(), -1., Kbd_[gp].values(), dHdbKbbi.A());
      Core::LinAlg::DenseFunctions::multiplyTN<double, 1, spintype, 1>(
          1., &(plmat->dHepDT(gp)), -1., (*KbT_)[gp].values(), dHdbKbbi.A());

      // TSI with EAS
      if (eastype_ != soh8p_easnone)
      {
        // error checks
        if (dHda == nullptr)
          FOUR_C_THROW("dHda is nullptr pointer");
        else if ((int)dHda->size() != numgpt_)
          FOUR_C_THROW("dHda has wrong size");

        switch (eastype_)
        {
          case soh8p_easmild:
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
                1., dHda->at(gp).values(), 1., M->values(), dHdC.values());
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
                1., dHda->at(gp).values(), -1., Kba_->at(gp).values(), dHdbKbbi.values());

            break;
          case soh8p_easfull:
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
                1., dHda->at(gp).values(), 1., M->values(), dHdC.values());
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
                1., dHda->at(gp).values(), -1., Kba_->at(gp).values(), dHdbKbbi.values());
            break;
          case soh8p_eassosh8:
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, numstr_, 1>(
                1., dHda->at(gp).values(), 1., M->values(), dHdC.values());
            Core::LinAlg::DenseFunctions::multiplyTN<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
                1., dHda->at(gp).values(), -1., Kba_->at(gp).values(), dHdbKbbi.values());
            break;
          case soh8p_easnone:
            break;
          default:
            FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }
    }  // TSI

    if (eastype_ != soh8p_easnone)
    {
      // condense plasticity into EAS matrix blocks
      tmp.shape(neas_, spintype);
      switch (eastype_)
      {
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              1., Kda->values(), -1., KdbKbb.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype, spintype>(
              0., tmp.values(), 1., Kab.values(), KbbInv_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype,
              numdofperelement_>(1., Kad_->values(), -1., tmp.values(), Kbd_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              1., KaaInv_->values(), -1., tmp.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
              1., feas_->values(), -1., tmp.values(), fbeta_[gp].values());
          if (eval_tsi)
          {
            Core::LinAlg::SerialDenseMatrix kbTm(spintype, nen_);
            Core::LinAlg::Matrix<nen_, 1> shapefunct;
            Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
            Core::LinAlg::DenseFunctions::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.values(), 1., KbT_->at(gp).values(), shapefunct.values());
            Core::LinAlg::DenseFunctions::multiply<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype, nen_>(
                1., KaT_->values(), -1., tmp.values(), kbTm.values());
          }
          break;
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              1., Kda->values(), -1., KdbKbb.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype, spintype>(
              0., tmp.values(), 1., Kab.values(), KbbInv_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype,
              numdofperelement_>(1., Kad_->values(), -1., tmp.values(), Kbd_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              1., KaaInv_->values(), -1., tmp.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
              1., feas_->values(), -1., tmp.values(), fbeta_[gp].values());
          if (eval_tsi)
          {
            Core::LinAlg::SerialDenseMatrix kbTm(spintype, nen_);
            Core::LinAlg::Matrix<nen_, 1> shapefunct;
            Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
            Core::LinAlg::DenseFunctions::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.values(), 1., KbT_->at(gp).values(), shapefunct.A());
            Core::LinAlg::DenseFunctions::multiply<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype, nen_>(
                1., KaT_->values(), -1., tmp.values(), kbTm.values());
          }
          break;
        case soh8p_eassosh8:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
              1., Kda->values(), -1., KdbKbb.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype, spintype>(
              0., tmp.values(), 1., Kab.values(), KbbInv_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype,
              numdofperelement_>(1., Kad_->values(), -1., tmp.values(), Kbd_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
              1., KaaInv_->values(), -1., tmp.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
              1., feas_->values(), -1., tmp.values(), fbeta_[gp].values());
          if (eval_tsi)
          {
            Core::LinAlg::SerialDenseMatrix kbTm(spintype, nen_);
            Core::LinAlg::Matrix<nen_, 1> shapefunct;
            Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
            Core::LinAlg::DenseFunctions::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.values(), 1., KbT_->at(gp).values(), shapefunct.A());
            Core::LinAlg::DenseFunctions::multiply<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype, nen_>(
                1., KaT_->values(), -1., tmp.values(), kbTm.values());
          }
          break;
        case soh18p_eassosh18:
          Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas>(
              1., Kda->values(), -1., KdbKbb.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype, spintype>(
              0., tmp.values(), 1., Kab.values(), KbbInv_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype,
              numdofperelement_>(1., Kad_->values(), -1., tmp.values(), Kbd_[gp].values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas>(
              1., KaaInv_->values(), -1., tmp.values(), Kba_->at(gp).values());
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype, 1>(
              1., feas_->values(), -1., tmp.values(), fbeta_[gp].values());
          if (eval_tsi)
          {
            Core::LinAlg::SerialDenseMatrix kbTm(spintype, nen_);
            Core::LinAlg::Matrix<nen_, 1> shapefunct;
            Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
            Core::LinAlg::DenseFunctions::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.values(), 1., KbT_->at(gp).values(), shapefunct.A());
            Core::LinAlg::DenseFunctions::multiply<double,
                PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype, nen_>(
                1., KaT_->values(), -1., tmp.values(), kbTm.values());
          }
          break;
        case soh8p_easnone:
          // do nothing
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    if (is_nitsche_contact_)
    {
      Core::LinAlg::Matrix<numstr_, spintype> tmp;
      tmp.Multiply(
          d_cauchy_db, Core::LinAlg::Matrix<spintype, spintype>(KbbInv_.at(gp).values(), true));
      cauchy_.at(gp).Multiply(
          -1., tmp, Core::LinAlg::Matrix<spintype, 1>(fbeta_.at(gp).values(), true), 1.);
      cauchy_deriv_.at(gp).Multiply(-1., tmp,
          Core::LinAlg::Matrix<spintype, numdofperelement_>(Kbd_.at(gp).values(), true), 1.);

      if (d_cauchy_dT_ptr)
      {
        Core::LinAlg::Matrix<6, 1> bla;
        bla.Multiply(-1., tmp, Core::LinAlg::Matrix<spintype, 1>(KbT_->at(gp).values(), true), 1.);
        cauchy_deriv_T_.at(gp).MultiplyNT(1., bla, shape_function(), 1.);
      }
    }
    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
  }

  // simple condensation as kbb is diagonal
  // This destinction is not necessary. However, we use the known diagonal
  // structure to avoid the inversion of the 5x5 matrix kbb.
  else
  {
    if (dDp_last_iter_[gp].normInf() > 0.)
    {
      if (force != nullptr)
        Core::LinAlg::DenseFunctions::multiply<double, numdofperelement_, spintype, 1>(
            1., force->A(), -1., kdbeta.values(), dDp_last_iter_[gp].values());

      switch (eastype_)
      {
        case soh8p_easnone:
          break;
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
              1., feas_->values(), -1., Kab.values(), dDp_last_iter_[gp].values());
          break;
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
              1., feas_->values(), -1., Kab.values(), dDp_last_iter_[gp].values());
          break;
        case soh8p_eassosh8:
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
              1., feas_->values(), -1., Kab.values(), dDp_last_iter_[gp].values());
          break;
        case soh18p_eassosh18:
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, spintype, 1>(
              1., feas_->values(), -1., Kab.values(), dDp_last_iter_[gp].values());
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    KbbInv_[gp].putScalar(0.0);
    for (int i = 0; i < KbbInv_[gp].numRows(); ++i) KbbInv_[gp](i, i) = 1. / (plmat->cpl());
    fbeta_[gp].putScalar(0.0);
    fbeta_[gp] = dDp_last_iter_[gp];
    fbeta_[gp].scale(plmat->cpl());
    Kbd_[gp].putScalar(0.0);
    if (eastype_ != soh8p_easnone) Kba_->at(gp).shape(spintype, neas_);
    if (KbT_ != Teuchos::null) KbT_->at(gp).putScalar(0.0);
  }

  return;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::build_delta_lp(const int gp)
{
  // current plastic flow increment
  set_delta_lp()(0, 0) = dDp_last_iter_[gp](0);
  set_delta_lp()(1, 1) = dDp_last_iter_[gp](1);
  set_delta_lp()(2, 2) = -1.0 * (dDp_last_iter_[gp](0) + dDp_last_iter_[gp](1));
  set_delta_lp()(0, 1) = dDp_last_iter_[gp](2);
  set_delta_lp()(1, 0) = dDp_last_iter_[gp](2);
  set_delta_lp()(1, 2) = dDp_last_iter_[gp](3);
  set_delta_lp()(2, 1) = dDp_last_iter_[gp](3);
  set_delta_lp()(0, 2) = dDp_last_iter_[gp](4);
  set_delta_lp()(2, 0) = dDp_last_iter_[gp](4);
  if (have_plastic_spin())
  {
    set_delta_lp()(0, 1) += dDp_last_iter_[gp](5);
    set_delta_lp()(1, 0) -= dDp_last_iter_[gp](5);
    set_delta_lp()(1, 2) += dDp_last_iter_[gp](6);
    set_delta_lp()(2, 1) -= dDp_last_iter_[gp](6);
    set_delta_lp()(0, 2) += dDp_last_iter_[gp](7);
    set_delta_lp()(2, 0) -= dDp_last_iter_[gp](7);
  }
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::recover_plasticity_and_eas(
    const Core::LinAlg::Matrix<numdofperelement_, 1>* res_d,
    const Core::LinAlg::Matrix<nen_, 1>* res_T)
{
  if (str_params_interface().is_default_step())
  {
    if (eastype_ != soh8p_easnone) recover_eas(res_d, res_T);


    if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    {
      Core::LinAlg::Matrix<nen_, 1> shapefunct;
      double res_t = -1.;
      double* res_t_ptr;
      for (int gp = 0; gp < numgpt_; ++gp)
      {
        if (res_T)
        {
          Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
          res_t = shapefunct.Dot(*res_T);
          res_t_ptr = &res_t;
        }
        else
          res_t_ptr = nullptr;

        // recover plastic variables
        if (have_plastic_spin())
          recover_plasticity<plspin>(res_d, gp, res_t_ptr);
        else
          recover_plasticity<zerospin>(res_d, gp, res_t_ptr);
      }
    }
  }
  else
  {
    //    const double new_step_length   = str_params_interface().get_step_length();
    //
    //    if (eastype_!=soh8p_easnone)
    //      reduce_eas_step(new_step_length,old_step_length_);
    //
    //    if (Material()->MaterialType()==Core::Materials::m_plelasthyper)
    //      for (int gp=0;gp<numgpt_;++gp)
    //        reduce_plasticity_step(new_step_length,old_step_length_,gp);
    //
    //    old_step_length_=new_step_length;
  }
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::recover_eas(
    const Core::LinAlg::Matrix<numdofperelement_, 1>* res_d,
    const Core::LinAlg::Matrix<nen_, 1>* res_T)
{
  if (eastype_ == soh8p_easnone) return;

  const double step_length = str_params_interface().get_step_length();

  // first, store the eas state of the previous accepted Newton step
  str_params_interface().sum_into_my_previous_sol_norm(
      NOX::Nln::StatusTest::quantity_eas, neas_, alpha_eas_->values(), Owner());

  if (str_params_interface().is_default_step()) switch (eastype_)
    {
      case soh8p_easmild:
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numdofperelement_, 1>(
            1.0, feas_->values(), 1.0, Kad_->values(), res_d->A());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, nen_, 1>(
              1., feas_->values(), 1., KaT_->values(), res_T->values());
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        Core::LinAlg::DenseFunctions::update<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_easfull:
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numdofperelement_, 1>(
            1.0, feas_->values(), 1.0, Kad_->values(), res_d->values());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, nen_, 1>(
              1., feas_->values(), 1., KaT_->values(), res_T->values());
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        Core::LinAlg::DenseFunctions::update<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_eassosh8:
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, numdofperelement_, 1>(
            1.0, feas_->values(), 1.0, Kad_->values(), res_d->values());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, nen_, 1>(
              1., feas_->values(), 1., KaT_->values(), res_T->values());
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        Core::LinAlg::DenseFunctions::update<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh18p_eassosh18:
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, numdofperelement_, 1>(
            1.0, feas_->values(), 1.0, Kad_->values(), res_d->values());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          Core::LinAlg::DenseFunctions::multiply<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, nen_, 1>(
              1., feas_->values(), 1., KaT_->values(), res_T->values());
        Core::LinAlg::DenseFunctions::multiply<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        Core::LinAlg::DenseFunctions::update<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_easnone:
        break;
      default:
        FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  else
    FOUR_C_THROW("no line search implemented yet");

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas, neas_,
      alpha_eas_inc_->values(), alpha_eas_->values(), step_length, Owner());

  return;
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
template <int spintype>
void Discret::ELEMENTS::So3Plast<distype>::recover_plasticity(
    const Core::LinAlg::Matrix<numdofperelement_, 1>* res_d, const int gp, const double* res_t)
{
  const double step_length = str_params_interface().get_step_length();

  if (str_params_interface().is_default_step() == false)
    FOUR_C_THROW("no line search implemented yet");

  // first, store the state of the previous accepted Newton step
  str_params_interface().sum_into_my_previous_sol_norm(
      NOX::Nln::StatusTest::quantity_plasticity, spintype, dDp_last_iter_[gp].values(), Owner());

  // temporary Epetra matrix
  Core::LinAlg::SerialDenseVector tmp_v(spintype);
  Core::LinAlg::SerialDenseMatrix tmp_m(spintype, numdofperelement_);

  // first part
  Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype, 1>(
      0., dDp_inc_[gp].values(), -1., KbbInv_[gp].values(), fbeta_[gp].values());

  // second part
  Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype, numdofperelement_>(
      0., tmp_m.values(), 1., KbbInv_[gp].values(), Kbd_[gp].values());
  Core::LinAlg::DenseFunctions::multiply<double, spintype, numdofperelement_, 1>(
      1., dDp_inc_[gp].values(), -1., tmp_m.values(), res_d->values());

  // thermal part
  if (KbT_ != Teuchos::null)
    if (res_t)
    {
      Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype, 1>(
          1., dDp_inc_[gp].values(), -1. * (*res_t), KbbInv_[gp].values(), (*KbT_)[gp].values());
    }

  // EAS part
  if (eastype_ != soh8p_easnone)
  {
    tmp_m.shape(spintype, neas_);
    switch (eastype_)
    {
      case soh8p_easmild:
        Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
            0., tmp_m.values(), 1., KbbInv_[gp].values(), Kba_->at(gp).values());
        Core::LinAlg::DenseFunctions::multiply<double, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
            1., dDp_inc_[gp].values(), -1., tmp_m.values(), alpha_eas_inc_->values());
        break;
      case soh8p_easfull:
        Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
            0., tmp_m.values(), 1., KbbInv_[gp].values(), Kba_->at(gp).values());
        Core::LinAlg::DenseFunctions::multiply<double, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
            1., dDp_inc_[gp].values(), -1., tmp_m.values(), alpha_eas_inc_->values());
        break;
      case soh8p_eassosh8:
        Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas>(
            0., tmp_m.values(), 1., KbbInv_[gp].values(), Kba_->at(gp).values());
        Core::LinAlg::DenseFunctions::multiply<double, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            1., dDp_inc_[gp].values(), -1., tmp_m.values(), alpha_eas_inc_->values());
        break;
      case soh18p_eassosh18:
        Core::LinAlg::DenseFunctions::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas>(
            0., tmp_m.values(), 1., KbbInv_[gp].values(), Kba_->at(gp).values());
        Core::LinAlg::DenseFunctions::multiply<double, spintype,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            1., dDp_inc_[gp].values(), -1., tmp_m.values(), alpha_eas_inc_->values());
        break;
      case soh8p_easnone:
        break;
      default:
        FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }  // EAS part

  Core::LinAlg::DenseFunctions::update<double, spintype, 1>(
      1., dDp_last_iter_[gp], 1., dDp_inc_[gp]);

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_plasticity,
      spintype, dDp_inc_[gp].values(), dDp_inc_[gp].values(), step_length, Owner());
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::reduce_eas_step(
    const double new_step_length, const double old_step_length)
{
  if (eastype_ == soh8p_easnone) return;

  Core::LinAlg::Update(-1., *alpha_eas_inc_, 1., *alpha_eas_);
  alpha_eas_inc_->scale(new_step_length / old_step_length);
  Core::LinAlg::Update(+1., *alpha_eas_inc_, 1., *alpha_eas_);

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas, neas_,
      alpha_eas_inc_->values(), alpha_eas_->values(), new_step_length, Owner());
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::reduce_plasticity_step(
    const double new_step_length, const double old_step_length, const int gp)
{
  Core::LinAlg::Update(-1., dDp_inc_[gp], 1., dDp_last_iter_[gp]);
  dDp_inc_[gp].scale(new_step_length / old_step_length);
  Core::LinAlg::Update(+1., dDp_inc_[gp], 1., dDp_last_iter_[gp]);

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_plasticity,
      plspintype_, dDp_inc_[gp].values(), dDp_inc_[gp].values(), new_step_length, Owner());
}

/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 07/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::update_plastic_deformation_nln(PlSpinType spintype)
{
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
  {
    // loop over all Gauss points
    for (int gp = 0; gp < numgpt_; gp++)
    {
      build_delta_lp(gp);
      static_cast<Mat::PlasticElastHyper*>(Material().get())->UpdateGP(gp, &delta_lp());

      KbbInv_[gp].putScalar(0.0);
      Kbd_[gp].putScalar(0.0);
      fbeta_[gp].putScalar(0.0);
      if (tsi_) (*KbT_)[gp].putScalar(0.0);
    }
  }
  else
  {
    SolidMaterial()->Update();
  }

  if (eastype_ != soh8p_easnone)
  {
    for (int i = 0; i < neas_; i++)
    {
      (*alpha_eas_delta_over_last_timestep_)(i) = (*alpha_eas_)(i) - (*alpha_eas_last_timestep_)(i);
      (*alpha_eas_last_timestep_)(i) = (*alpha_eas_)(i);
    }
    Kad_->putScalar(0.0);
    KaaInv_->putScalar(0.0);
    feas_->putScalar(0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculate internal energy of the element (private)                  |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
double Discret::ELEMENTS::So3Plast<distype>::calc_int_energy(
    std::vector<double>& disp,         // current displacements
    std::vector<double>& temperature,  // current temperatuere
    Teuchos::ParameterList& params)    // strain output option
{
  invalid_ele_data();

  double energy = 0.;

  std::vector<double> bla;
  fill_position_arrays(disp, bla, temperature);

  if (fbar_ || eastype_ != soh8p_easnone)
    evaluate_center();  // deformation gradient at centroid of element

  if (eastype_ != soh8p_easnone) eas_setup();

  // get plastic hyperelastic material
  Mat::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());
  else
    FOUR_C_THROW("elastic strain energy in so3plast elements only for plastic material");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    invalid_gp_data();

    // shape functions (shapefunct) and their first derivatives (deriv)
    evaluate_shape(xsi_[gp]);
    evaluate_shape_deriv(xsi_[gp]);

    kinematics(gp);

    double gp_temp = shape_function().Dot(temp());

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8p_easnone)
    {
      eas_shape(gp);
      eas_enhance_strains();
    }
    // calculate modified deformation gradient
    else if (fbar_)
      setup_fbar_gp();
    else
      set_defgrd_mod() = defgrd();

    const double psi = plmat->StrainEnergyTSI(defgrd_mod(), gp, Id(), gp_temp);

    const double detJ_w = det_j() * wgt_[gp];
    energy += detJ_w * psi;

  }  // gp loop

  return energy;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::get_cauchy_n_dir_and_derivatives_at_xi_elast(
    const Core::LinAlg::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd2,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dT,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dT)
{
  if (distype == Core::FE::CellType::nurbs27) get_nurbs_ele_info();

  if (fbar_ || eastype_ != soh8p_easnone)
    FOUR_C_THROW("cauchy stress not available for fbar or eas elements");

  if (temp || d_cauchyndir_dT || d2_cauchyndir_dd_dT)
    if (!temp || !d_cauchyndir_dT || !d2_cauchyndir_dd_dT)
      FOUR_C_THROW("inconsistent temperature dependency input");
  if (temp && Material()->MaterialType() != Core::Materials::m_plelasthyper)
  {
    FOUR_C_THROW(
        "thermo-mechanical Nitsche contact only with PlasticElastHyper"
        "\nIf you want to do elasticity, set a negative yield stress ;)");
  }

  auto* plmat = dynamic_cast<Mat::PlasticElastHyper*>(Material().get());

  cauchy_n_dir = 0.0;

  static Core::LinAlg::Matrix<nen_, nsd_> xrefe(true);  // reference coord. of element
  static Core::LinAlg::Matrix<nen_, nsd_> xcurr(true);  // current  coord. of element
  static Core::LinAlg::Matrix<nen_, 1> ele_temp(true);
  xrefe.Clear();
  xcurr.Clear();
  ele_temp.Clear();
  Core::Nodes::Node** nodes = Nodes();

  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int d = 0; d < nsd_; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * nsd_ + d];
    }
    if (temp)
      if (!temp->empty()) ele_temp(i) = temp->at(i);
  }

  evaluate_shape(xi);
  evaluate_shape_deriv(xi);

  const double gp_temp = ele_temp.Dot(shape_function());
  double d_cauchyndir_dT_gp = 0.0;
  static Core::LinAlg::Matrix<nsd_, 1> d_T_dxi(true);
  d_T_dxi.Multiply(1.0, deriv_shape_function(), ele_temp, 0.0);

  static Core::LinAlg::Matrix<nsd_, nen_> N_XYZ(true);
  static Core::LinAlg::Matrix<nsd_, nsd_> invJ(true);
  invJ.Multiply(1.0, deriv_shape_function(), xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, deriv_shape_function(), 0.0);
  static Core::LinAlg::Matrix<nsd_, nsd_> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static Core::LinAlg::Matrix<9, numdofperelement_> d_F_dd(true);
  d_F_dd.Clear();
  if (d_cauchyndir_dd || d2_cauchyndir_dd_dn || d2_cauchyndir_dd_ddir || d2_cauchyndir_dd_dxi ||
      d2_cauchyndir_dd_dT)
  {
    for (int i = 0; i < nen_; ++i)
    {
      d_F_dd(0, nsd_ * i + 0) = N_XYZ(0, i);
      d_F_dd(1, nsd_ * i + 1) = N_XYZ(1, i);
      d_F_dd(2, nsd_ * i + 2) = N_XYZ(2, i);
      d_F_dd(3, nsd_ * i + 0) = N_XYZ(1, i);
      d_F_dd(4, nsd_ * i + 1) = N_XYZ(2, i);
      d_F_dd(5, nsd_ * i + 0) = N_XYZ(2, i);
      d_F_dd(6, nsd_ * i + 1) = N_XYZ(0, i);
      d_F_dd(7, nsd_ * i + 2) = N_XYZ(1, i);
      d_F_dd(8, nsd_ * i + 2) = N_XYZ(0, i);
    }
  }

  static Core::LinAlg::Matrix<9, 1> d_cauchyndir_dF(true);
  static Core::LinAlg::Matrix<9, 9> d2_cauchyndir_dF2(true);
  static Core::LinAlg::Matrix<9, nsd_> d2_cauchyndir_dF_dn(true);
  static Core::LinAlg::Matrix<9, nsd_> d2_cauchyndir_dF_ddir(true);
  static Core::LinAlg::Matrix<9, 1> d2_cauchyndir_dF_dT(true);

  if (plmat && temp)
  {
    plmat->evaluate_cauchy_n_dir_and_derivatives(defgrd, n, dir, cauchy_n_dir, d_cauchyndir_dn,
        d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2, &d2_cauchyndir_dF_dn,
        &d2_cauchyndir_dF_ddir, -1, Id(), nullptr, &gp_temp, &d_cauchyndir_dT_gp,
        &d2_cauchyndir_dF_dT);
  }
  else
  {
    SolidMaterial()->evaluate_cauchy_n_dir_and_derivatives(defgrd, n, dir, cauchy_n_dir,
        d_cauchyndir_dn, d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2,
        &d2_cauchyndir_dF_dn, &d2_cauchyndir_dF_ddir, -1, Id(), nullptr, nullptr, nullptr, nullptr);
  }

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->shape(numdofperelement_, 1);
    Core::LinAlg::Matrix<numdofperelement_, 1> d_cauchyndir_dd_mat(d_cauchyndir_dd->values(), true);
    d_cauchyndir_dd_mat.MultiplyTN(1.0, d_F_dd, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dT)
  {
    d2_cauchyndir_dd_dT->shape(numdofperelement_, nen_);
    static Core::LinAlg::Matrix<numdofperelement_, 1> tmp(true);
    tmp.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dT, 0.0);
    Core::LinAlg::Matrix<numdofperelement_, nen_>(d2_cauchyndir_dd_dT->values(), true)
        .MultiplyNT(tmp, shape_function());
  }

  if (d2_cauchyndir_dd_dn)
  {
    d2_cauchyndir_dd_dn->shape(numdofperelement_, nsd_);
    Core::LinAlg::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_dn_mat(
        d2_cauchyndir_dd_dn->values(), true);
    d2_cauchyndir_dd_dn_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dn, 0.0);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->shape(numdofperelement_, nsd_);
    Core::LinAlg::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_ddir_mat(
        d2_cauchyndir_dd_ddir->values(), true);
    d2_cauchyndir_dd_ddir_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_ddir, 0.0);
  }

  if (d_cauchyndir_dT)
  {
    d_cauchyndir_dT->shape(nen_, 1);
    Core::LinAlg::Matrix<nen_, 1>(d_cauchyndir_dT->values(), true)
        .Update(d_cauchyndir_dT_gp, shape_function(), 1.0);
  }

  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->shape(numdofperelement_, numdofperelement_);
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> d2_cauchyndir_dd2_mat(
        d2_cauchyndir_dd2->values(), true);
    static Core::LinAlg::Matrix<9, numdofperelement_> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd2_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF2_d_F_dd, 0.0);
  }

  // prepare evaluation of d_cauchyndir_dxi or d2_cauchyndir_dd_dxi
  static Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, nen_> deriv2(true);
  static Core::LinAlg::Matrix<9, nsd_> d_F_dxi(true);
  deriv2.Clear();
  d_F_dxi.Clear();

  if (d_cauchyndir_dxi or d2_cauchyndir_dd_dxi)
  {
    if (distype == Core::FE::CellType::nurbs27)
    {
      Core::FE::Nurbs::nurbs_get_3D_funct_deriv_deriv2(set_shape_function(),
          set_deriv_shape_function(), deriv2, xi, knots(), weights(), distype);
    }
    else
      Core::FE::shape_function_deriv2<distype>(xi, deriv2);

    static Core::LinAlg::Matrix<nen_, nsd_> xXF(true);
    static Core::LinAlg::Matrix<nsd_, Core::FE::DisTypeToNumDeriv2<distype>::numderiv2> xXFsec(
        true);
    xXF.Update(1.0, xcurr, 0.0);
    xXF.MultiplyNT(-1.0, xrefe, defgrd, 1.0);
    xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

    for (int a = 0; a < nsd_; ++a)
    {
      for (int b = 0; b < nsd_; ++b)
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
    if (temp) d_cauchyndir_dxi->Update(d_cauchyndir_dT_gp, d_T_dxi, 1.0);
  }

  if (d2_cauchyndir_dd_dxi)
  {
    d2_cauchyndir_dd_dxi->shape(numdofperelement_, nsd_);
    Core::LinAlg::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_dxi_mat(
        d2_cauchyndir_dd_dxi->values(), true);

    static Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, nsd_> Xsec(true);
    static Core::LinAlg::Matrix<nen_, 6> N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static Core::LinAlg::Matrix<9, numdofperelement_> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d2_cauchyndir_dF2_d_F_dd, d_F_dxi, 0.0);

    if (temp)
    {
      static Core::LinAlg::Matrix<9, nsd_> tmp(true);
      tmp.MultiplyNT(1.0, d2_cauchyndir_dF_dT, d_T_dxi, 0.0);
      d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d_F_dd, tmp, 1.0);
    }

    static Core::LinAlg::Matrix<9, nsd_ * numdofperelement_> d2_F_dxi_dd(true);
    d2_F_dxi_dd.Clear();
    for (int i = 0; i < nsd_; ++i)
    {
      for (int j = 0; j < nsd_; ++j)
      {
        for (int k = 0; k < nen_; ++k)
        {
          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), numdofpernode_ * (numdofpernode_ * k + i) + 0) +=
              deriv2(0, k) * invJ(j, 0) + deriv2(3, k) * invJ(j, 1) + deriv2(4, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 0) * invJ(j, 0) - N_XYZ_Xsec(k, 3) * invJ(j, 1) -
              N_XYZ_Xsec(k, 4) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), numdofpernode_ * (numdofpernode_ * k + i) + 1) +=
              deriv2(3, k) * invJ(j, 0) + deriv2(1, k) * invJ(j, 1) + deriv2(5, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 3) * invJ(j, 0) - N_XYZ_Xsec(k, 1) * invJ(j, 1) -
              N_XYZ_Xsec(k, 5) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), numdofpernode_ * (numdofpernode_ * k + i) + 2) +=
              deriv2(4, k) * invJ(j, 0) + deriv2(5, k) * invJ(j, 1) + deriv2(2, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 4) * invJ(j, 0) - N_XYZ_Xsec(k, 5) * invJ(j, 1) -
              N_XYZ_Xsec(k, 2) * invJ(j, 2);

          for (int l = 0; l < nsd_; ++l)
          {
            d2_cauchyndir_dd_dxi_mat(k * 3 + i, l) +=
                d_cauchyndir_dF(VoigtMapping::NonSymToVoigt9(i, j), 0) *
                d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                    numdofpernode_ * (numdofpernode_ * k + i) + l);
          }
        }
      }
    }
  }
  invalid_ele_data();
  invalid_gp_data();
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::get_cauchy_n_dir_and_derivatives_at_xi_plast(
    const Core::LinAlg::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd2,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dT,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dT)
{
  if (distype != Core::FE::CellType::hex8 || numgpt_ != 8) FOUR_C_THROW("only for hex8 with 8 gp");
  if (Material()->MaterialType() != Core::Materials::m_plelasthyper)
    FOUR_C_THROW("only PlasticElastHyper materials here");
  if ((int)cauchy_.size() != numgpt_ || (int)cauchy_deriv_.size() != numgpt_)
    FOUR_C_THROW("have you evaluated the cauchy stress???");

  cauchy_n_dir = 0.0;
  if (d_cauchyndir_dxi) d_cauchyndir_dxi->Clear();
  if (d_cauchyndir_dT) d_cauchyndir_dT->shape(nen_, 1);
  if (d2_cauchyndir_dd_dT) d2_cauchyndir_dd_dT->shape(numdofperelement_, nen_);

  Core::LinAlg::Matrix<3, 3> n_dir_dir_n(true);
  n_dir_dir_n.MultiplyNT(.5, n, dir, 1.);
  n_dir_dir_n.MultiplyNT(.5, dir, n, 1.);
  Core::LinAlg::Matrix<6, 1> n_dir_dir_n_v;
  for (int i = 0; i < 3; ++i) n_dir_dir_n_v(i) = n_dir_dir_n(i, i);
  n_dir_dir_n_v(3) = n_dir_dir_n(0, 1) + n_dir_dir_n(1, 0);
  n_dir_dir_n_v(4) = n_dir_dir_n(2, 1) + n_dir_dir_n(1, 2);
  n_dir_dir_n_v(5) = n_dir_dir_n(0, 2) + n_dir_dir_n(2, 0);

  Core::LinAlg::Matrix<3, 1> xi_expol(xi);
  xi_expol.Scale(sqrt(3.));

  Core::LinAlg::Matrix<nen_, 1> shapefunct;
  Core::LinAlg::Matrix<nsd_, nen_> deriv;
  Core::FE::shape_function<distype>(xi_expol, shapefunct);
  Core::FE::shape_function_deriv1<distype>(xi_expol, deriv);

  Core::LinAlg::Matrix<numstr_, 1> cauchy_expol;
  Core::LinAlg::Matrix<numstr_, numdofperelement_> cauchy_deriv_expol;

  Core::LinAlg::Matrix<6, 1> tmp61;
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    cauchy_expol.Update(shapefunct(gp), cauchy_.at(gp), 1.);
    cauchy_deriv_expol.Update(shapefunct(gp), cauchy_deriv_.at(gp), 1.);
    if (d_cauchyndir_dxi)
      for (int d = 0; d < nsd_; ++d)
        (*d_cauchyndir_dxi)(d) += cauchy_.at(gp).Dot(n_dir_dir_n_v) * deriv(d, gp) * sqrt(3.);

    if (d_cauchyndir_dT)
      Core::LinAlg::Matrix<nen_, 1>(d_cauchyndir_dT->values(), true)
          .MultiplyTN(shapefunct(gp), cauchy_deriv_T_.at(gp), n_dir_dir_n_v, 1.);
  }

  cauchy_n_dir = cauchy_expol.Dot(n_dir_dir_n_v);

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->reshape(numdofperelement_, 1);
    Core::LinAlg::Matrix<numdofperelement_, 1>(d_cauchyndir_dd->values(), true)
        .MultiplyTN(cauchy_deriv_expol, n_dir_dir_n_v);
  }
  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->reshape(numdofperelement_, numdofperelement_);
    d2_cauchyndir_dd2->putScalar(0.0);
  }

  Core::LinAlg::Matrix<numstr_, nsd_> d_ndirdirn_v_dn, d_ndirdirn_v_dt;
  for (int i = 0; i < nsd_; ++i)
  {
    for (int j = 0; j < nsd_; ++j)
    {
      for (int a = 0; a < nsd_; ++a)
      {
        d_ndirdirn_v_dn(VoigtMapping::SymToVoigt6(i, j), a) +=
            .5 * ((i == a) * dir(j) + (j == a) * dir(i));
        d_ndirdirn_v_dt(VoigtMapping::SymToVoigt6(i, j), a) +=
            .5 * ((i == a) * n(j) + (j == a) * n(i));
      }
    }
  }
  if (d_cauchyndir_dn) d_cauchyndir_dn->MultiplyTN(d_ndirdirn_v_dn, cauchy_expol);
  if (d_cauchyndir_ddir) d_cauchyndir_ddir->MultiplyTN(d_ndirdirn_v_dt, cauchy_expol);

  if (d2_cauchyndir_dd_dn)
  {
    d2_cauchyndir_dd_dn->reshape(numdofperelement_, nsd_);
    Core::LinAlg::Matrix<numdofperelement_, nsd_>(d2_cauchyndir_dd_dn->values(), true)
        .MultiplyTN(cauchy_deriv_expol, d_ndirdirn_v_dn);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->reshape(numdofperelement_, nsd_);
    Core::LinAlg::Matrix<numdofperelement_, nsd_>(d2_cauchyndir_dd_ddir->values(), true)
        .MultiplyTN(cauchy_deriv_expol, d_ndirdirn_v_dt);
  }
  if (d2_cauchyndir_dd_dxi)
  {
    d2_cauchyndir_dd_dxi->reshape(numdofperelement_, nsd_);
    d2_cauchyndir_dd_dxi->putScalar(0.0);
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::get_cauchy_n_dir_and_derivatives_at_xi(
    const Core::LinAlg::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd2,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dT,
    Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
    double* d_cauchyndir_dc)
{
  if (d_cauchyndir_dc != nullptr) FOUR_C_THROW("Not implemented");

  bool elastic = true;
  auto* plmat = dynamic_cast<Mat::PlasticElastHyper*>(Material().get());
  if (!plmat)
    elastic = true;
  else
    elastic = plmat->AllElastic();

  if (!elastic)
  {
    get_cauchy_n_dir_and_derivatives_at_xi_plast(xi, disp, n, dir, cauchy_n_dir, d_cauchyndir_dd,
        d2_cauchyndir_dd2, d2_cauchyndir_dd_dn, d2_cauchyndir_dd_ddir, d2_cauchyndir_dd_dxi,
        d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dxi, temp, d_cauchyndir_dT,
        d2_cauchyndir_dd_dT);
  }
  else
  {
    get_cauchy_n_dir_and_derivatives_at_xi_elast(xi, disp, n, dir, cauchy_n_dir, d_cauchyndir_dd,
        d2_cauchyndir_dd2, d2_cauchyndir_dd_dn, d2_cauchyndir_dd_ddir, d2_cauchyndir_dd_dxi,
        d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dxi, temp, d_cauchyndir_dT,
        d2_cauchyndir_dd_dT);
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::output_strains(const int gp,
    const Inpar::STR::StrainType iostrain,                 // strain output option
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestrain  // strains at GP
)
{
  // strain output *********************************
  if (eastype_ != soh8p_easnone && numgpt_ != 8)
  {
    // no stress output currently
  }
  else
  {
    switch (iostrain)
    {
      case Inpar::STR::strain_gl:
      {
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        Core::LinAlg::Matrix<numstr_, 1> total_glstrain(false);
        total_glstrain(0) = 0.5 * (rcg()(0, 0) - 1.0);
        total_glstrain(1) = 0.5 * (rcg()(1, 1) - 1.0);
        total_glstrain(2) = 0.5 * (rcg()(2, 2) - 1.0);
        total_glstrain(3) = rcg()(0, 1);
        total_glstrain(4) = rcg()(1, 2);
        total_glstrain(5) = rcg()(2, 0);

        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = total_glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * total_glstrain(i);
      }
      break;
      case Inpar::STR::strain_ea:
      {
        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");

        // inverse of deformation gradient
        Core::LinAlg::Matrix<3, 3> invdefgrd;
        invdefgrd.Invert(defgrd());

        static Core::LinAlg::Matrix<3, 3> tmp1;
        static Core::LinAlg::Matrix<3, 3> total_euler_almansi(true);
        tmp1.MultiplyTN(invdefgrd, invdefgrd);
        total_euler_almansi.MultiplyTN(-1., invdefgrd, tmp1);
        for (int i = 0; i < 3; i++) total_euler_almansi(i, i) = 1.;
        total_euler_almansi.Scale(0.5);

        (*elestrain)(gp, 0) = total_euler_almansi(0, 0);
        (*elestrain)(gp, 1) = total_euler_almansi(1, 1);
        (*elestrain)(gp, 2) = total_euler_almansi(2, 2);
        (*elestrain)(gp, 3) = total_euler_almansi(0, 1);
        (*elestrain)(gp, 4) = total_euler_almansi(1, 2);
        (*elestrain)(gp, 5) = total_euler_almansi(0, 2);
      }
      break;
      case Inpar::STR::strain_none:
        break;
      default:
      {
        FOUR_C_THROW("requested strain type not available");
        break;
      }
    }
  }
  // end of strain output **************************
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::output_stress(const int gp,
    const Inpar::STR::StressType iostress,                 // strain output option
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestress  // strains at GP
)
{
  // return gp stresses
  switch (iostress)
  {
    case Inpar::STR::stress_2pk:
    {
      if (elestress == nullptr) FOUR_C_THROW("stress data not available");
      for (int i = 0; i < numstr_; ++i) (*elestress)(gp, i) = p_k2()(i);
    }
    break;
    case Inpar::STR::stress_cauchy:
    {
      if (elestress == nullptr) FOUR_C_THROW("stress data not available");

      static Core::LinAlg::Matrix<3, 3> pkstress;
      pkstress(0, 0) = p_k2()(0);
      pkstress(0, 1) = p_k2()(3);
      pkstress(0, 2) = p_k2()(5);
      pkstress(1, 0) = pkstress(0, 1);
      pkstress(1, 1) = p_k2()(1);
      pkstress(1, 2) = p_k2()(4);
      pkstress(2, 0) = pkstress(0, 2);
      pkstress(2, 1) = pkstress(1, 2);
      pkstress(2, 2) = p_k2()(2);

      static Core::LinAlg::Matrix<3, 3> cauchystress;
      static Core::LinAlg::Matrix<3, 3> tmp1;
      tmp1.Multiply(1.0 / det_f(), defgrd(), pkstress);
      cauchystress.MultiplyNT(tmp1, defgrd());

      (*elestress)(gp, 0) = cauchystress(0, 0);
      (*elestress)(gp, 1) = cauchystress(1, 1);
      (*elestress)(gp, 2) = cauchystress(2, 2);
      (*elestress)(gp, 3) = cauchystress(0, 1);
      (*elestress)(gp, 4) = cauchystress(1, 2);
      (*elestress)(gp, 5) = cauchystress(0, 2);
    }
    break;
    case Inpar::STR::stress_none:
      break;
    default:
    {
      FOUR_C_THROW("requested stress type not available");
      break;
    }
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::kinematics(const int gp)
{
  // compute Jacobian matrix and determinant
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
  // derivatives of coordinates w.r.t material coordinates xjm_ = dx/ds
  set_inv_j().Multiply(deriv_shape_function(), xrefe());
  // xij_ = ds/dx
  set_det_j() = set_inv_j().Invert();
  if (det_j() < 1.0E-16) FOUR_C_THROW("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", det_j());

  /* get the inverse of the Jacobian matrix which looks like:
   **            [ x_,r  y_,r  z_,r ]^-1
   **     J^-1 = [ x_,s  y_,s  z_,s ]
   **            [ x_,t  y_,t  z_,t ]
   */
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  set_deriv_shape_function_xyz().Multiply(inv_j(), deriv_shape_function());  // (6.21)

  // (material) deformation gradient
  // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
  set_defgrd().MultiplyTT(xcurr(), deriv_shape_function_xyz());

  // inverse deformation gradient and determinant
  set_det_f() = set_inv_defgrd().Invert(defgrd());

  // calcualte total rcg
  set_rcg().MultiplyTN(defgrd(), defgrd());

  // total rcg in strain-like voigt notation
  for (int i = 0; i < 3; i++) set_rc_gvec()(i) = rcg()(i, i);
  set_rc_gvec()(3) = rcg()(0, 1) * 2.;
  set_rc_gvec()(4) = rcg()(1, 2) * 2.;
  set_rc_gvec()(5) = rcg()(0, 2) * 2.;

  // calculate nonlinear B-operator
  calculate_bop(&set_bop(), &defgrd(), &deriv_shape_function_xyz(), gp);

  // build plastic velocity gradient from condensed variables
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper && gp >= 0 && gp < numgpt_)
    build_delta_lp(gp);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::integrate_mass_matrix(
    const int gp, Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>& mass)
{
  const double density = Material()->Density(gp);
  // integrate consistent mass matrix
  const double factor = det_j() * wgt_[gp] * density;
  double ifactor, massfactor;
  for (int inod = 0; inod < nen_; ++inod)
  {
    ifactor = shape_function()(inod) * factor;
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      massfactor = shape_function()(jnod) * ifactor;  // intermediate factor
      mass(3 * inod + 0, 3 * jnod + 0) += massfactor;
      mass(3 * inod + 1, 3 * jnod + 1) += massfactor;
      mass(3 * inod + 2, 3 * jnod + 2) += massfactor;
    }
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::integrate_stiff_matrix(const int gp,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>& stiff,
    Core::LinAlg::SerialDenseMatrix& Kda)
{
  const double detJ_w = det_j() * wgt_[gp];

  // integrate `elastic' and `initial-displacement' stiffness matrix
  // keu = keu + (B^T . C . B) * detJ * w(gp)
  static Core::LinAlg::Matrix<numstr_, numdofperelement_> cb;
  cb.Multiply(cmat(), bop());
  if (fbar_)
    stiff.MultiplyTN(detJ_w * fbar_fac(), bop(), cb, 1.0);
  else
    stiff.MultiplyTN(detJ_w, bop(), cb, 1.0);

  // integrate `geometric' stiffness matrix and add to keu *****************
  Core::LinAlg::Matrix<numstr_, 1> sfac(p_k2());  // auxiliary integrated stress
  if (fbar_)
    sfac.Scale(detJ_w / fbar_fac());  // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
  else
    sfac.Scale(detJ_w);  // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
  double SmB_L[nsd_];    // intermediate Sm.B_L
  // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
  for (int inod = 0; inod < nen_; ++inod)
  {
    SmB_L[0] = sfac(0) * deriv_shape_function_xyz()(0, inod) +
               sfac(3) * deriv_shape_function_xyz()(1, inod) +
               sfac(5) * deriv_shape_function_xyz()(2, inod);
    SmB_L[1] = sfac(3) * deriv_shape_function_xyz()(0, inod) +
               sfac(1) * deriv_shape_function_xyz()(1, inod) +
               sfac(4) * deriv_shape_function_xyz()(2, inod);
    SmB_L[2] = sfac(5) * deriv_shape_function_xyz()(0, inod) +
               sfac(4) * deriv_shape_function_xyz()(1, inod) +
               sfac(2) * deriv_shape_function_xyz()(2, inod);
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      double bopstrbop = 0.0;  // intermediate value
      for (int idim = 0; idim < 3; ++idim)
        bopstrbop += deriv_shape_function_xyz()(idim, jnod) * SmB_L[idim];
      stiff(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
      stiff(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
      stiff(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
    }
  }  // end of integrate `geometric' stiffness******************************

  // integrate additional fbar matrix**************************************
  if (fbar_)
  {
    static Core::LinAlg::Matrix<numstr_, 1> ccg;
    ccg.Multiply(cmat(), rc_gvec());

    static Core::LinAlg::Matrix<numdofperelement_, 1> bopccg(false);  // auxiliary integrated stress
    bopccg.MultiplyTN(detJ_w * fbar_fac() / 3.0, bop(), ccg);

    static Core::LinAlg::Matrix<numdofperelement_, 1> bops(false);  // auxiliary integrated stress
    bops.MultiplyTN(-detJ_w / fbar_fac() / 3.0, bop(), p_k2());
    stiff.MultiplyNT(1., bops, htensor(), 1.);
    stiff.MultiplyNT(1., bopccg, htensor(), 1.);
  }
  // end of integrate additional fbar matrix*****************************

  // EAS technology: integrate matrices --------------------------------- EAS
  if (not(str_params_interface().get_predictor_type() == Inpar::STR::pred_tangdis))
    if (eastype_ != soh8p_easnone)
    {
      // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
      // integrate Kda: Kad += (M^T . cmat . B) * detJ * w(gp)
      // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
      static Core::LinAlg::SerialDenseMatrix cM(numstr_, neas_);  // temporary c . M
      switch (eastype_)
      {
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiply<double, numstr_, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              cM.values(), cmat().A(), m_eas().values());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              1.0, *KaaInv_, detJ_w, m_eas(), cM);
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_,
              numdofperelement_>(1.0, Kad_->values(), detJ_w, m_eas().values(), cb.A());
          Core::LinAlg::DenseFunctions::multiplyTN<double, numdofperelement_, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas>(
              1.0, Kda.values(), detJ_w, cb.A(), m_eas().values());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              1.0, feas_->values(), detJ_w, m_eas().values(), p_k2().A());
          break;
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiply<double, numstr_, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              cM.values(), cmat().A(), m_eas().values());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              1.0, *KaaInv_, detJ_w, m_eas(), cM);
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_,
              numdofperelement_>(1.0, Kad_->values(), detJ_w, m_eas().values(), cb.A());
          Core::LinAlg::DenseFunctions::multiplyTN<double, numdofperelement_, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas>(
              1.0, Kda.values(), detJ_w, cb.A(), m_eas().values());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              1.0, feas_->values(), detJ_w, m_eas().values(), p_k2().A());
          break;
        case soh8p_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }  // ---------------------------------------------------------------- EAS
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::integrate_force(
    const int gp, Core::LinAlg::Matrix<numdofperelement_, 1>& force)
{
  if (fbar_)
    force.MultiplyTN(det_j() * wgt_[gp] / fbar_fac(), bop(), p_k2(), 1.0);
  else
    force.MultiplyTN(det_j() * wgt_[gp], bop(), p_k2(), 1.0);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::integrate_thermo_gp(
    const int gp, Core::LinAlg::SerialDenseVector& dHda)
{
  const double timefac_d = str_params_interface().get_tim_int_factor_vel() /
                           str_params_interface().get_tim_int_factor_disp();
  const double detJ_w = det_j() * wgt_[gp];

  // get plastic hyperelastic material
  Mat::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == Core::Materials::m_plelasthyper)
    plmat = static_cast<Mat::PlasticElastHyper*>(Material().get());
  else
    FOUR_C_THROW("tsi only with m_plelasthyper material type");

  // Gauss point temperature
  const double gp_temp = temp().Dot(shape_function());

  // volumetric part of K_dT = Gough-Joule effect**********************************
  // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccc
  // get the thermal material tangent
  Core::LinAlg::Matrix<numstr_, 1> cTvol(true);
  Core::LinAlg::Matrix<numstr_, numstr_> dcTvoldE;
  plmat->EvaluateCTvol(&defgrd_mod(), &cTvol, &dcTvoldE, gp, Id());
  // end of call material law ccccccccccccccccccccccccccccccccccccccccccccc
  if (fbar_)
    (*dFintdT_)[gp].MultiplyTN(detJ_w / (fbar_fac()), bop(), cTvol, 0.);
  else
    (*dFintdT_)[gp].MultiplyTN(detJ_w, bop(), cTvol, 0.);
  if (eastype_ != soh8p_easnone)
  {
    Core::LinAlg::Matrix<numstr_, nen_> cTm;
    cTm.MultiplyNT(cTvol, shape_function());
    switch (eastype_)
    {
      case soh8p_easfull:
        Core::LinAlg::DenseFunctions::multiplyTN<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, nen_>(
            1., KaT_->values(), detJ_w, m_eas().values(), cTm.A());
        break;
      case soh8p_easmild:
        Core::LinAlg::DenseFunctions::multiplyTN<double,
            PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, nen_>(
            1., KaT_->values(), detJ_w, m_eas().values(), cTm.A());
        break;
      case soh8p_easnone:
        break;
      default:
        FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }

  // elastic heating ******************************************************
  plmat->HepDiss(gp) = 0.;
  plmat->dHepDT(gp) = 0.;
  plmat->dHepDissDd(gp).size(numdofperelement_);
  plmat->dHepDissDd(gp).putScalar(0.0);
  if (eastype_ == soh8p_easnone)
  {
    if (fbar_)
    {
      Core::LinAlg::Matrix<3, 3> defgrd_rate_0;
      defgrd_rate_0.MultiplyTT(xcurr_rate(), deriv_shape_function_xyz_0());

      double he_fac;
      double he_fac_deriv;
      plmat->EvaluateGoughJoule(det_f_0(), gp, Id(), he_fac, he_fac_deriv);

      double fiddfdot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) fiddfdot += inv_defgrd_0()(j, i) * defgrd_rate_0(i, j);

      double j_dot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) j_dot += det_f_0() * inv_defgrd_0()(j, i) * defgrd_rate_0(i, j);

      double He = he_fac * gp_temp * j_dot;

      plmat->HepDiss(gp) = He;

      // derivative of elastic heating w.r.t. temperature *******************
      plmat->dHepDT(gp) = he_fac * j_dot;

      Core::LinAlg::Matrix<numdofperelement_, 1> deriv_jdot_d(true);
      Core::LinAlg::Matrix<numdofperelement_, 1> deriv_j_d(true);

      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_j_d(n) +=
              det_f_0() * inv_defgrd_0()(i, n % 3) * deriv_shape_function_xyz_0()(i, n / 3);

      Core::LinAlg::Matrix<3, 3> tmp;
      tmp.Multiply(inv_defgrd_0(), defgrd_rate_0);
      Core::LinAlg::Matrix<3, 3> tmp2;
      tmp2.Multiply(tmp, inv_defgrd_0());

      deriv_jdot_d.Update(fiddfdot, deriv_j_d, 1.);
      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_jdot_d(n) += det_f_0() * inv_defgrd_0()(i, n % 3) *
                                 deriv_shape_function_xyz_0()(i, n / 3) * timefac_d -
                             det_f_0() * tmp2(i, n % 3) * deriv_shape_function_xyz_0()(i, n / 3);

      Core::LinAlg::Matrix<numdofperelement_, 1> dHedd(true);
      dHedd.Update(gp_temp * he_fac, deriv_jdot_d, 1.);
      dHedd.Update(he_fac_deriv * gp_temp * j_dot, deriv_j_d, 1.);

      Core::LinAlg::DenseFunctions::update<double, numdofperelement_, 1>(
          plmat->dHepDissDd(gp).values(), dHedd.A());
    }
    else
    {
      Core::LinAlg::Matrix<3, 3> defgrd_rate;
      defgrd_rate.MultiplyTT(xcurr_rate(), deriv_shape_function_xyz());

      double he_fac;
      double he_fac_deriv;
      plmat->EvaluateGoughJoule(det_f(), gp, Id(), he_fac, he_fac_deriv);

      double fiddfdot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) fiddfdot += inv_defgrd()(j, i) * defgrd_rate(i, j);

      double j_dot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) j_dot += det_f() * inv_defgrd()(j, i) * defgrd_rate(i, j);

      double He = he_fac * gp_temp * j_dot;

      plmat->HepDiss(gp) = He;

      // derivative of elastic heating w.r.t. temperature *******************
      plmat->dHepDT(gp) = he_fac * j_dot;

      Core::LinAlg::Matrix<numdofperelement_, 1> deriv_jdot_d(true);
      Core::LinAlg::Matrix<numdofperelement_, 1> deriv_j_d(true);

      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_j_d(n) += det_f() * inv_defgrd()(i, n % 3) * deriv_shape_function_xyz()(i, n / 3);

      Core::LinAlg::Matrix<3, 3> tmp;
      tmp.Multiply(inv_defgrd(), defgrd_rate);
      Core::LinAlg::Matrix<3, 3> tmp2;
      tmp2.Multiply(tmp, inv_defgrd());

      deriv_jdot_d.Update(fiddfdot, deriv_j_d, 1.);
      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_jdot_d(n) +=
              det_f() * inv_defgrd()(i, n % 3) * deriv_shape_function_xyz()(i, n / 3) * timefac_d -
              det_f() * tmp2(i, n % 3) * deriv_shape_function_xyz()(i, n / 3);

      Core::LinAlg::Matrix<numdofperelement_, 1> dHedd(true);
      dHedd.Update(gp_temp * he_fac, deriv_jdot_d, 1.);
      dHedd.Update(he_fac_deriv * gp_temp * j_dot, deriv_j_d, 1.);

      Core::LinAlg::DenseFunctions::update<double, numdofperelement_, 1>(
          plmat->dHepDissDd(gp).values(), dHedd.A());
    }
  }
  else
  {
    Core::LinAlg::Matrix<3, 3> defgrd_rate;
    defgrd_rate.MultiplyTT(xcurr_rate(), deriv_shape_function_xyz());

    // Gough-Joule effect
    plmat->HepDiss(gp) = 0.;
    plmat->dHepDT(gp) = 0.;
    plmat->dHepDissDd(gp).size(numdofperelement_);
    // Like this it should be easier to do EAS as well
    Core::LinAlg::Matrix<3, 3> RCGrate;
    RCGrate.MultiplyTN(defgrd_rate, defgrd());
    RCGrate.MultiplyTN(1., defgrd(), defgrd_rate, 1.);
    Core::LinAlg::Matrix<6, 1> RCGrateVec;
    for (int i = 0; i < 3; ++i) RCGrateVec(i, 0) = RCGrate(i, i);
    RCGrateVec(3, 0) = 2. * RCGrate(0, 1);
    RCGrateVec(4, 0) = 2. * RCGrate(1, 2);
    RCGrateVec(5, 0) = 2. * RCGrate(0, 2);

    // enhance the deformation rate
    if (eastype_ != soh8p_easnone)
    {
      Core::LinAlg::SerialDenseVector alpha_dot(neas_);
      switch (eastype_)
      {
        case soh8p_easmild:
          // calculate EAS-rate
          Core::LinAlg::DenseFunctions::update<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
              0., alpha_dot, 1., *alpha_eas_);
          Core::LinAlg::DenseFunctions::update<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., alpha_dot, -1., *alpha_eas_last_timestep_);
          alpha_dot.scale(timefac_d);
          // enhance the strain rate
          // factor 2 because we deal with RCGrate and not GLrate
          Core::LinAlg::DenseFunctions::multiply<double, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., RCGrateVec.A(), 2., m_eas().values(), alpha_dot.values());
          break;
        case soh8p_easfull:
          // calculate EAS-rate
          Core::LinAlg::DenseFunctions::update<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
              0., alpha_dot, 1., *alpha_eas_);
          Core::LinAlg::DenseFunctions::update<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., alpha_dot, -1., *alpha_eas_last_timestep_);
          alpha_dot.scale(timefac_d);
          // enhance the strain rate
          // factor 2 because we deal with RCGrate and not GLrate
          Core::LinAlg::DenseFunctions::multiply<double, numstr_,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., RCGrateVec.A(), 2., m_eas().values(), alpha_dot.values());
          break;
        case soh8p_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }  // enhance the deformation rate

    // heating ************************************************************
    double He = .5 * gp_temp * cTvol.Dot(RCGrateVec);

    plmat->HepDiss(gp) = He;

    // derivative of elastic heating w.r.t. temperature *******************
    plmat->dHepDT(gp) = .5 * cTvol.Dot(RCGrateVec);

    // derivative of elastic heating w.r.t. displacement ******************
    Core::LinAlg::Matrix<numdofperelement_, 1> dHedd(true);
    Core::LinAlg::Matrix<6, nen_ * nsd_> boprate(false);  // (6x24)
    calculate_bop(&boprate, &defgrd_rate, &deriv_shape_function_xyz(), gp);

    Core::LinAlg::Matrix<6, 1> tmp61;
    tmp61.MultiplyTN(dcTvoldE, RCGrateVec);
    dHedd.MultiplyTN(.5 * gp_temp, bop(), tmp61, 1.);

    dHedd.MultiplyTN(timefac_d * gp_temp, bop(), cTvol, 1.);
    dHedd.MultiplyTN(gp_temp, boprate, cTvol, 1.);

    // derivative of elastic heating w.r.t. EAS alphas *******************
    if (eastype_ != soh8p_easnone)
    {
      switch (eastype_)
      {
        case soh8p_easmild:
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              0., dHda.values(), .5 * gp_temp, m_eas().values(), tmp61.A());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              1., dHda.values(), gp_temp * timefac_d, m_eas().values(), cTvol.A());
          break;
        case soh8p_easfull:
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              0., dHda.values(), .5 * gp_temp, m_eas().values(), tmp61.A());
          Core::LinAlg::DenseFunctions::multiplyTN<double,
              PlastEasTypeToNumEas<Discret::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              1., dHda.values(), gp_temp * timefac_d, m_eas().values(), cTvol.A());
          break;
        case soh8p_easnone:
          break;
        default:
          FOUR_C_THROW("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    plmat->dHepDissDd(gp).putScalar(0.0);
    Core::LinAlg::DenseFunctions::update<double, numdofperelement_, 1>(
        plmat->dHepDissDd(gp).values(), dHedd.A());
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::HeatFlux(const std::vector<double>& temperature,
    const std::vector<double>& disp, const Core::LinAlg::Matrix<nsd_, 1>& xi,
    const Core::LinAlg::Matrix<nsd_, 1>& n, double& q, Core::LinAlg::SerialDenseMatrix* dq_dT,
    Core::LinAlg::SerialDenseMatrix* dq_dd, Core::LinAlg::Matrix<nsd_, 1>* dq_dn,
    Core::LinAlg::Matrix<nsd_, 1>* dq_dpxi, Core::LinAlg::SerialDenseMatrix* d2q_dT_dd,
    Core::LinAlg::SerialDenseMatrix* d2q_dT_dn, Core::LinAlg::SerialDenseMatrix* d2q_dT_dpxi)
{
  if (!dq_dT || !dq_dd || !dq_dn || !dq_dpxi || !d2q_dT_dd || !d2q_dT_dn || !d2q_dT_dpxi)
    FOUR_C_THROW("input inconsistent");

  invalid_gp_data();
  invalid_ele_data();

  if (NumMaterial() < 2) FOUR_C_THROW("where's my second material");
  Teuchos::RCP<Mat::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<Mat::FourierIso>(Material(1), true);
  const double k0 = mat_thr->Conductivity();

  std::vector<double> vel(0);
  fill_position_arrays(disp, vel, temperature);

  // shape functions (shapefunct) and their first derivatives (deriv)
  Core::FE::shape_function<distype>(xi, set_shape_function());
  Core::FE::shape_function_deriv1<distype>(xi, set_deriv_shape_function());

  kinematics();

  Core::LinAlg::Matrix<3, 1> GradT;
  GradT.Multiply(deriv_shape_function_xyz(), temp());

  Core::LinAlg::Matrix<3, 1> iFn;
  iFn.Multiply(inv_defgrd(), n);

  Core::LinAlg::Matrix<9, numdofperelement_> dFdd;
  for (int k = 0; k < nen_; ++k)
  {
    for (int d = 0; d < 3; ++d) dFdd(d, k * nsd_ + d) = deriv_shape_function_xyz()(d, k);
    dFdd(3, k * nsd_ + 0) = deriv_shape_function_xyz()(1, k);
    dFdd(4, k * nsd_ + 1) = deriv_shape_function_xyz()(2, k);
    dFdd(5, k * nsd_ + 0) = deriv_shape_function_xyz()(2, k);

    dFdd(6, k * nsd_ + 1) = deriv_shape_function_xyz()(0, k);
    dFdd(7, k * nsd_ + 2) = deriv_shape_function_xyz()(1, k);
    dFdd(8, k * nsd_ + 2) = deriv_shape_function_xyz()(0, k);
  }

  q = -k0 / det_f() * GradT.Dot(iFn);

  Core::LinAlg::Matrix<nen_, 9> dq_dT_dF;

  if (dq_dT)
  {
    dq_dT->shape(nen_, 1);
    for (int n = 0; n < nen_; ++n)
      for (int d = 0; d < nsd_; ++d)
        (*dq_dT)(n, 0) += -k0 / det_f() * deriv_shape_function_xyz()(d, n) * iFn(d);
  }

  if (d2q_dT_dd || dq_dd)
  {
    d2q_dT_dd->shape(nen_, nen_ * nsd_);
    Core::LinAlg::Matrix<nen_, nen_ * nsd_> d2q_dT_dd_m(d2q_dT_dd->values(), true);

    Core::LinAlg::Matrix<3, nen_> tmp;
    tmp.MultiplyTN(inv_defgrd(), deriv_shape_function_xyz());
    Core::LinAlg::Matrix<nen_, 1> tmp2;
    tmp2.MultiplyTN(deriv_shape_function_xyz(), iFn);

    Core::LinAlg::Matrix<9, nen_> dq_dF_v;
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          dq_dF_v(VoigtMapping::NonSymToVoigt9(a, b), c) =
              k0 / det_f() * (tmp2(c) * inv_defgrd()(b, a) + tmp(a, c) * iFn(b));

    d2q_dT_dd_m.MultiplyTN(dq_dF_v, dFdd);

    dq_dd->shape(nen_ * nsd_, 1);
    Core::LinAlg::Matrix<nen_ * nsd_, 1> dq_dd_m(dq_dd->values(), true);
    dq_dd_m.MultiplyTN(d2q_dT_dd_m, temp());
  }

  if (d2q_dT_dn || dq_dn)
  {
    d2q_dT_dn->shape(nen_, nsd_);
    Core::LinAlg::Matrix<nen_, nsd_> d2q_dT_dn_m(d2q_dT_dn->values(), true);
    d2q_dT_dn_m.MultiplyTN(-k0 / det_f(), deriv_shape_function_xyz(), inv_defgrd());

    dq_dn->MultiplyTN(d2q_dT_dn_m, temp());
  }

  if (dq_dpxi || d2q_dT_dpxi)
  {
    d2q_dT_dpxi->shape(nen_, nsd_);
    Core::LinAlg::Matrix<nen_, nsd_> d2q_dT_dpxi_m(d2q_dT_dpxi->values(), true);

    Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, nen_> deriv2;
    Core::FE::shape_function_deriv2<distype>(xi, deriv2);
    const int d2v[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
    Core::LinAlg::Matrix<3, 1> tmp;
    tmp.MultiplyTN(inv_j(), iFn);

    Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<distype>::numderiv2, nen_> tmp3;
    tmp3.Update(deriv2);

    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += -k0 / det_f() * tmp3(d2v[a][b], c) * tmp(b);

    Core::LinAlg::Matrix<nen_, nen_> tmp2;
    tmp2.Multiply(xrefe(), deriv_shape_function_xyz());
    tmp3.Multiply(deriv2, tmp2);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / det_f() * tmp3(d2v[a][b], c) * tmp(b);

    Core::LinAlg::Matrix<nen_, nsd_> xXF(xcurr());
    xXF.MultiplyNT(-1., xrefe(), defgrd(), 1.);

    Core::LinAlg::Matrix<nsd_, Core::FE::DisTypeToNumDeriv2<distype>::numderiv2> xXFsec;
    xXFsec.MultiplyTT(xXF, deriv2);

    Core::LinAlg::Matrix<nsd_, nsd_> tmp4;
    tmp4.MultiplyTN(inv_j(), inv_defgrd());

    Core::LinAlg::Matrix<nsd_, Core::FE::DisTypeToNumDeriv2<distype>::numderiv2> tmp5;
    tmp5.Multiply(tmp4, xXFsec);

    Core::LinAlg::Matrix<nen_, 1> tmp6;
    tmp6.MultiplyTN(deriv_shape_function_xyz(), iFn);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / det_f() * tmp6(c) * tmp5(b, d2v[a][b]);

    Core::LinAlg::Matrix<nsd_, nen_> tmp7;
    tmp7.MultiplyTN(inv_defgrd(), deriv_shape_function_xyz());
    tmp3.MultiplyTN(xXFsec, tmp7);
    tmp.MultiplyTN(inv_j(), iFn);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / det_f() * tmp3(d2v[a][b], c) * tmp(b);

    dq_dpxi->MultiplyTN(d2q_dT_dpxi_m, temp());
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::So3Plast<distype>::get_nurbs_ele_info(Core::FE::Discretization* dis)
{
  if (!IsNurbsElement()) return;

  if (dis == nullptr) dis = Global::Problem::Instance()->GetDis("structure").get();

  dynamic_cast<Discret::Nurbs::NurbsDiscretization*>(dis)->GetKnotVector()->GetEleKnots(
      set_knots(), Id());
  for (int i = 0; i < nen_; ++i)
    set_weights()(i) = dynamic_cast<Discret::Nurbs::ControlPoint*>(Nodes()[i])->W();
}

// template functions
template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::condense_plasticity<5>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::condense_plasticity<8>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>::condense_plasticity<5>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>::condense_plasticity<8>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>::condense_plasticity<5>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>::condense_plasticity<8>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>::condense_plasticity<5>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>::condense_plasticity<8>(
    const Core::LinAlg::Matrix<nsd_, nsd_>&, const Core::LinAlg::Matrix<nsd_, nsd_>&,
    const Core::LinAlg::Matrix<numstr_, numdofperelement_>&,
    const Core::LinAlg::Matrix<nsd_, nen_>*, const Core::LinAlg::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, Core::LinAlg::Matrix<numdofperelement_, 1>*,
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*,
    const Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*,
    std::vector<Core::LinAlg::SerialDenseVector>*, const double*,
    const Core::LinAlg::Matrix<numdofperelement_, 1>*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>::HeatFlux(
    const std::vector<double>&, const std::vector<double>&, const Core::LinAlg::Matrix<nsd_, 1>&,
    const Core::LinAlg::Matrix<nsd_, 1>&, double&, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::Matrix<nsd_, 1>*,
    Core::LinAlg::Matrix<nsd_, 1>*, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*);


template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>::HeatFlux(
    const std::vector<double>&, const std::vector<double>&, const Core::LinAlg::Matrix<nsd_, 1>&,
    const Core::LinAlg::Matrix<nsd_, 1>&, double&, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::Matrix<nsd_, 1>*,
    Core::LinAlg::Matrix<nsd_, 1>*, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>::HeatFlux(
    const std::vector<double>&, const std::vector<double>&, const Core::LinAlg::Matrix<nsd_, 1>&,
    const Core::LinAlg::Matrix<nsd_, 1>&, double&, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::Matrix<nsd_, 1>*,
    Core::LinAlg::Matrix<nsd_, 1>*, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*);

template void Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>::HeatFlux(
    const std::vector<double>&, const std::vector<double>&, const Core::LinAlg::Matrix<nsd_, 1>&,
    const Core::LinAlg::Matrix<nsd_, 1>&, double&, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::Matrix<nsd_, 1>*,
    Core::LinAlg::Matrix<nsd_, 1>*, Core::LinAlg::SerialDenseMatrix*,
    Core::LinAlg::SerialDenseMatrix*, Core::LinAlg::SerialDenseMatrix*);

FOUR_C_NAMESPACE_CLOSE
