/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 2


*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "baci_so3_plast_ssn.H"

#include "baci_so3_ssn_plast_fwd.hpp"

#include "baci_lib_globalproblem.H"
#include "baci_lib_voigt_notation.H"
#include "baci_mat_plasticelasthyper.H"
#include <Epetra_SerialDenseSolver.h>
#include "baci_mat_service.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_structure_new_gauss_point_data_output_manager.H"
#include "baci_nurbs_discret.H"

#include "baci_mat_fourieriso.H"
#include "baci_so3_element_service.H"

#include "baci_discretization_fem_general_utils_gauss_point_postprocess.H"

using VoigtMapping = UTILS::VOIGT::IndexMappings;

/*----------------------------------------------------------------------*
 | evaluate the element (public)                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  SetParamsInterfacePtr(params);

  InvalidEleData();
  if (distype == DRT::Element::nurbs27) GetNurbsEleInfo(&discretization);

  dsassert(kintype_ == INPAR::STR::kinem_nonlinearTotLag,
      "only geometricallly nonlinear formluation for plasticity!");

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
    act = ELEMENTS::String2ActionType(action);
  }

  // what should the element do
  switch (act)
  {
    //============================================================================
    // linear stiffness
    case DRT::ELEMENTS::struct_calc_linstiff:
    case DRT::ELEMENTS::struct_calc_linstiffmass:
    {
      dserror("linear kinematics version out dated");
      break;
    }

    //============================================================================
    // nonlinear stiffness
    case DRT::ELEMENTS::struct_calc_internalforce:
    {
      // internal force vector
      CORE::LINALG::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.A(), true);
      // elemat1+2, elevec2+3 are not used anyway

      // need current displacement and residual/incremental displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      if ((disp == Teuchos::null))
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(la[0].lm_.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      CORE::LINALG::Matrix<numdofperelement_, numdofperelement_> myemat(true);

      // initialise the vectors
      // Evaluate() is called the first time in StructureBaseAlgorithm: at this
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
          if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          DRT::UTILS::ExtractMyValues(*vel, myvel, la[0].lm_);
        }
        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            dserror("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          DRT::UTILS::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, &myemat, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);


      break;
    }

    //============================================================================
    // nonlinear stiffness
    case DRT::ELEMENTS::struct_calc_nlnstiff:
    {
      // stiffness
      CORE::LINALG::Matrix<numdofperelement_, numdofperelement_> elemat1(elemat1_epetra.A(), true);
      CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;
      CORE::LINALG::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.A(), true);

      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(la[0].lm_.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // Evaluate() is called the first time in StructureBaseAlgorithm: at this
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
          if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          DRT::UTILS::ExtractMyValues(*vel, myvel, la[0].lm_);
        }

        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            dserror("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          DRT::UTILS::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, matptr, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);

      break;
    }  // calc_struct_nlnstiff

    //============================================================================
    // (non)linear stiffness, mass matrix and internal force vector
    case DRT::ELEMENTS::struct_calc_nlnstiffmass:
    case DRT::ELEMENTS::struct_calc_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(la[0].lm_.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);
      // stiffness
      CORE::LINALG::Matrix<numdofperelement_, numdofperelement_> elemat1(elemat1_epetra.A(), true);
      // mass
      CORE::LINALG::Matrix<numdofperelement_, numdofperelement_> elemat2(elemat2_epetra.A(), true);
      // internal force
      CORE::LINALG::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.A(), true);

      // initialise the vectors
      // Evaluate() is called the first time in StructureBaseAlgorithm: at this
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
          if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          DRT::UTILS::ExtractMyValues(*vel, myvel, la[0].lm_);
        }
        if (discretization.HasState(1, "temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
          if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            dserror("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
          DRT::UTILS::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
        }
      }

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(mydisp, myvel, mytempnp, &elemat1, &elemat2, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);

      if (act == DRT::ELEMENTS::struct_calc_nlnstifflmass)
        // lump mass matrix
        // we assume #elemat2 is a square matrix
        for (int c = 0; c < elemat2_epetra.N(); ++c)  // parse columns
        {
          double d = 0.0;
          for (int r = 0; r < elemat2_epetra.M(); ++r)  // parse rows
          {
            d += elemat2(r, c);  // accumulate row entries
            elemat2(r, c) = 0.0;
          }
          elemat2(c, c) = d;  // apply sum of row entries on diagonal
        }
      break;
    }  // calc_struct_nlnstiff(l)mass

    case DRT::ELEMENTS::struct_calc_stress:
    {
      // elemat1+2,elevec1-3 are not used anyway

      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
        Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
        Teuchos::RCP<std::vector<char>> straindata = Teuchos::null;
        INPAR::STR::StressType iostress = INPAR::STR::stress_none;
        INPAR::STR::StrainType iostrain = INPAR::STR::strain_none;
        if (IsParamsInterface())
        {
          stressdata = StrParamsInterface().MutableStressDataPtr();
          straindata = StrParamsInterface().MutableStrainDataPtr();

          iostress = StrParamsInterface().GetStressOutputType();
          iostrain = StrParamsInterface().GetStrainOutputType();
        }
        else
        {
          stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
          straindata = params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
          iostress =
              DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
          iostrain =
              DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        }
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");

        std::vector<double> mydisp((la[0].lm_).size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);

        // initialise the vectors
        // Evaluate() is called the first time in StructureBaseAlgorithm: at this
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
            if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
            // extract the velocities
            myvel.resize((la[0].lm_).size());
            DRT::UTILS::ExtractMyValues(*vel, myvel, la[0].lm_);
          }
          if (discretization.HasState(1, "temperature"))
          {
            Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
            if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

            // the temperature field has only one dof per node, disregarded by the dimension of the
            // problem
            const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
            if (la[1].Size() != nen_ * numdofpernode_thr)
              dserror("Location vector length for temperature does not match!");
            // extract the current temperatures
            mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
            DRT::UTILS::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
          }
        }

        CORE::LINALG::Matrix<numgpt_post, numstr_> stress;
        CORE::LINALG::Matrix<numgpt_post, numstr_> strain;

        // default: geometrically non-linear analysis with Total Lagrangean approach
        nln_stiffmass(mydisp, myvel, mytempnp, nullptr, nullptr, nullptr, &stress, &strain, params,
            iostress, iostrain);
        {
          DRT::PackBuffer data;
          this->AddtoPack(data, stress);
          data.StartPacking();
          this->AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          this->AddtoPack(data, strain);
          data.StartPacking();
          this->AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }

      break;
    }


    //============================================================================
    // required for predictor TangDis --> can be helpful in compressible case!
    case DRT::ELEMENTS::struct_calc_reset_istep:
    {
      for (int i = 0; i < numgpt_; i++)
      {
        KbbInv_[i].Scale(0.);
        Kbd_[i].Scale(0.);
        fbeta_[i].Scale(0.);
      }
      break;
    }

    //============================================================================
    case DRT::ELEMENTS::struct_calc_update_istep:
    {
      // update plastic deformation
      // default: geometrically non-linear analysis with Total Lagrangean approach
      UpdatePlasticDeformation_nln(plspintype_);
      break;
    }  // calc_struct_update_istep

    //============================================================================
    case DRT::ELEMENTS::struct_calc_stifftemp:  // here we want to build the K_dT block for
                                                // monolithic TSI
    {
      // stiffness
      CORE::LINALG::Matrix<numdofperelement_, nen_> k_dT(elemat1_epetra.A(), true);

      // calculate matrix block
      nln_kdT_tsi(&k_dT, params);
    }
    break;

    case DRT::ELEMENTS::struct_calc_energy:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      std::vector<double> mydisp(la[0].lm_.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);

      std::vector<double> mytempnp(0);
      if (discretization.HasState(1, "temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1, "temperature");
        if (tempnp == Teuchos::null) dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the dimension of the
        // problem
        const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
        if (la[1].Size() != nen_ * numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempnp.resize(((la[0].lm_).size()) / nsd_, 0.0);
        DRT::UTILS::ExtractMyValues(*tempnp, mytempnp, la[1].lm_);
      }

      double intenergy = CalcIntEnergy(mydisp, mytempnp, params);

      if (IsParamsInterface())  // new structural time integration
      {
        StrParamsInterface().AddContributionToEnergyType(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        dsassert(elevec1_epetra.Length() < 1, "The given result vector is too short.");
        elevec1_epetra(0) = intenergy;
      }
    }
    break;

    case DRT::ELEMENTS::struct_calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> myres((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*res, myres, la[0].lm_);
      CORE::LINALG::Matrix<nen_ * nsd_, 1> res_d(myres.data(), true);

      std::vector<double> mytempres(0);
      CORE::LINALG::Matrix<nen_, 1> res_t;
      CORE::LINALG::Matrix<nen_, 1>* res_t_ptr = nullptr;
      if (discretization.NumDofSets() > 1)
        if (discretization.HasState(1, "residual temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempres =
              discretization.GetState(1, "residual temperature");
          if (tempres == Teuchos::null) dserror("Cannot get state vector 'tempres'");

          // the temperature field has only one dof per node, disregarded by the dimension of the
          // problem
          const int numdofpernode_thr = discretization.NumDof(1, Nodes()[0]);
          if (la[1].Size() != nen_ * numdofpernode_thr)
            dserror("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempres.resize(((la[0].lm_).size()) / nsd_, 0.0);
          DRT::UTILS::ExtractMyValues(*tempres, mytempres, la[1].lm_);
          res_t = CORE::LINALG::Matrix<nen_, 1>(mytempres.data(), true);
          res_t_ptr = &res_t;
        }

      RecoverPlasticityAndEAS(&res_d, res_t_ptr);
    }
    break;

    case DRT::ELEMENTS::struct_calc_predict:
    {
      switch (StrParamsInterface().GetPredictorType())
      {
        case INPAR::STR::pred_constdis:
        default:
          for (int gp = 0; gp < numgpt_; ++gp) dDp_last_iter_[gp].Scale(0.);
          if (eastype_ != soh8p_easnone)
            for (int i = 0; i < neas_; i++) (*alpha_eas_)(i) = (*alpha_eas_last_timestep_)(i);
          break;
        case INPAR::STR::pred_constvel:
        case INPAR::STR::pred_tangdis:
          // do nothing for plasticity
          if (eastype_ != soh8p_easnone)
            for (int i = 0; i < neas_; i++)
              (*alpha_eas_)(i) =
                  0. + (*alpha_eas_last_timestep_)(i) + (*alpha_eas_delta_over_last_timestep_)(i);
          break;
      }
    }
    break;

    case DRT::ELEMENTS::struct_init_gauss_point_data_output:
    {
      dsassert(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");
      // Save number of Gauss of the element for gauss point data output
      StrParamsInterface().MutableGaussPointDataOutputManagerPtr()->AddElementNumberOfGaussPoints(
          numgpt_post);
      // holder for output quantity names and their size
      std::unordered_map<std::string, int> quantities_map{};
      // Ask material for the output quantity names and sizes
      SolidMaterial()->RegisterVtkOutputDataNames(quantities_map);
      // Add quantities to the Gauss point output data manager (if they do not already exist)
      StrParamsInterface().MutableGaussPointDataOutputManagerPtr()->MergeQuantities(quantities_map);
    }
    break;

    case ELEMENTS::struct_gauss_point_data_output:
    {
      dsassert(IsParamsInterface(),
          "This action type should only be called from the new time integration framework!");

      // Collection and assembly of gauss point data
      for (const auto& [quantity_name, quantity_size] :
          StrParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetQuantities())
      {
        // Step 1: Collect the data for each Gauss point for the material
        CORE::LINALG::SerialDenseMatrix gp_data(numgpt_post, quantity_size, true);
        bool data_available = SolidMaterial()->EvaluateVtkOutputData(quantity_name, gp_data);

        // Step 2: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
        // point)
        if (data_available)
        {
          switch (StrParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetOutputType())
          {
            case INPAR::STR::GaussPointDataOutputType::element_center:
            {
              // compute average of the quantities
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  StrParamsInterface()
                      .MutableGaussPointDataOutputManagerPtr()
                      ->GetMutableElementCenterData()
                      .at(quantity_name);
              CORE::DRT::ELEMENTS::AssembleAveragedElementValues<CORE::LINALG::SerialDenseMatrix>(
                  *global_data, gp_data, *this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::nodes:
            {
              Teuchos::RCP<Epetra_MultiVector> global_data =
                  StrParamsInterface()
                      .MutableGaussPointDataOutputManagerPtr()
                      ->GetMutableNodalData()
                      .at(quantity_name);

              Epetra_IntVector& global_nodal_element_count =
                  *StrParamsInterface()
                       .MutableGaussPointDataOutputManagerPtr()
                       ->GetMutableNodalDataCount()
                       .at(quantity_name);

              if (distype == DRT::Element::hex8)
              {
                if (quantity_size == 1)
                {
                  CORE::LINALG::Matrix<numgpt_post, 1> gp_data_view(gp_data, true);
                  soh8_expol(gp_data_view, *global_data);
                }
                else if (quantity_size == 9)
                {
                  CORE::LINALG::Matrix<numgpt_post, 9> gp_data_view(gp_data, true);
                  soh8_expol(gp_data_view, *global_data);
                }
              }
              else
                dserror(
                    "only element centered and Gauss point material internal variables for "
                    "so3_plast");

              DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, this);
              break;
            }
            case INPAR::STR::GaussPointDataOutputType::gauss_points:
            {
              std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
                  StrParamsInterface()
                      .MutableGaussPointDataOutputManagerPtr()
                      ->GetMutableGaussPointData()
                      .at(quantity_name);
              DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, this);
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

    //============================================================================
    default:
      dserror("Unknown type of action for So3_Plast");
      break;
  }  // action

  return 0;
}  // Evaluate()



/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::CalculateBop(
    CORE::LINALG::Matrix<numstr_, numdofperelement_>* bop,
    const CORE::LINALG::Matrix<nsd_, nsd_>* defgrd, const CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ,
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
}  // CalculateBop()



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 04/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = condition.Get<std::vector<int>>("onoff");
  const auto* val = condition.Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < nsd_)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = nsd_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const auto* funct = condition.Get<std::vector<int>>("funct");
  CORE::LINALG::Matrix<nsd_, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < nsd_; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct;

  // update element geometry
  CORE::LINALG::Matrix<nen_, nsd_> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < nen_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    CORE::LINALG::Matrix<nen_, 1> shapefunct;
    CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
    CORE::LINALG::Matrix<nsd_, nen_> deriv;
    CORE::DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp], deriv);

    // compute the Jacobian matrix
    CORE::LINALG::Matrix<nsd_, nsd_> jac;
    jac.Multiply(deriv, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

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
          (functnum > 0) ? DRT::Problem::Instance()
                               ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(functnum - 1)
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::InitJacobianMapping()
{
  return;
}  // InitJacobianMapping()}


/*----------------------------------------------------------------------*
 | internal force, stiffness and mass for f-bar elements    seitz 07/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::nln_stiffmass(
    std::vector<double>& disp,  // current displacements
    std::vector<double>& vel,   // current velocities
    std::vector<double>& temp,  // current temperatures
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*
        stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>* massmatrix,  // element mass matrix
    CORE::LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
    CORE::LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
    CORE::LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,                         // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,                  // stress output option
    const INPAR::STR::StrainType iostrain                   // strain output option
)
{
  const bool eval_tsi = (temp.size() != 0);
  const bool is_tangDis = StrParamsInterface().GetPredictorType() == INPAR::STR::pred_tangdis;
  const double dt = StrParamsInterface().GetDeltaTime();
  const double timefac_d =
      StrParamsInterface().GetTimIntFactorVel() / StrParamsInterface().GetTimIntFactorDisp();
  if (timefac_d <= 0 || dt <= 0) dserror("time integration parameters not provided");

  FillPositionArrays(disp, vel, temp);

  if (fbar_ || eastype_ != soh8p_easnone)
    EvaluateCenter();  // deformation gradient at centroid of element

  if (eastype_ != soh8p_easnone) EasSetup();

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    plmat = static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
  {
    if (tsi_) dserror("TSI with so3Plast elements only with PlasticElastHyper material");
  }

  // EAS matrix block
  Epetra_SerialDenseMatrix Kda(numdofperelement_, neas_);
  std::vector<Epetra_SerialDenseVector> dHda(0);
  if (eastype_ != soh8p_easnone && eval_tsi) dHda.resize(numgpt_, Epetra_SerialDenseVector(neas_));
  // temporary Epetra matrix for this and that
  Epetra_SerialDenseMatrix tmp;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    InvalidGpData();

    // shape functions (shapefunct) and their first derivatives (deriv)
    EvaluateShape(xsi_[gp]);
    EvaluateShapeDeriv(xsi_[gp]);

    Kinematics(gp);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8p_easnone)
    {
      EasShape(gp);
      EasEnhanceStrains();
    }
    // calculate modified deformation gradient
    else if (fbar_)
      SetupFbarGp();
    else
      SetDefgrdMod() = Defgrd();

    OutputStrains(gp, iostrain, elestrain);

    // Gauss point temperature
    double gp_temp = -1.e12;
    if (eval_tsi) gp_temp = Temp().Dot(ShapeFunction());

    // material call *********************************************
    if (plmat != nullptr)
    {
      plmat->EvaluateElast(&DefgrdMod(), &DeltaLp(), &SetPK2(), &SetCmat(), gp, Id());
      if (eval_tsi)
        plmat->EvaluateThermalStress(&DefgrdMod(), gp_temp, &SetPK2(), &SetCmat(), gp, Id());
    }
    else
    {
      static CORE::LINALG::Matrix<numstr_, 1> total_glstrain(false);
      total_glstrain(0) = 0.5 * (RCG()(0, 0) - 1.0);
      total_glstrain(1) = 0.5 * (RCG()(1, 1) - 1.0);
      total_glstrain(2) = 0.5 * (RCG()(2, 2) - 1.0);
      total_glstrain(3) = RCG()(0, 1);
      total_glstrain(4) = RCG()(1, 2);
      total_glstrain(5) = RCG()(2, 0);

      SolidMaterial()->Evaluate(
          &DefgrdMod(), &total_glstrain, params, &SetPK2(), &SetCmat(), gp, Id());
    }
    // material call *********************************************

    OutputStress(gp, iostress, elestress);

    // integrate usual internal force and stiffness matrix
    double detJ_w = DetJ() * wgt_[gp];
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != nullptr) IntegrateForce(gp, *force);

    // update stiffness matrix
    if (stiffmatrix != nullptr) IntegrateStiffMatrix(gp, *stiffmatrix, Kda);

    if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
      IntegrateMassMatrix(gp, *massmatrix);

    if (eval_tsi && (stiffmatrix != nullptr || force != nullptr) && plmat != nullptr)
      IntegrateThermoGp(gp, dHda[gp]);

    // plastic modifications
    if (plmat != nullptr)
      if (!plmat->AllElastic())
        if ((stiffmatrix != nullptr || force != nullptr) && !is_tangDis)
        {
          if (HavePlasticSpin())
          {
            if (fbar_)
              CondensePlasticity<plspin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  &RCGvec(), detJ_w, gp, gp_temp, params, force, stiffmatrix, nullptr, nullptr,
                  nullptr, &FbarFac(), &Htensor());
            else if (eastype_ != soh8p_easnone)
              CondensePlasticity<plspin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  nullptr, detJ_w, gp, gp_temp, params, force, stiffmatrix, &M_eas(), &Kda, &dHda);
            else
              CondensePlasticity<plspin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  nullptr, detJ_w, gp, gp_temp, params, force, stiffmatrix);
          }
          else
          {
            if (fbar_)
              CondensePlasticity<zerospin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  &RCGvec(), detJ_w, gp, gp_temp, params, force, stiffmatrix, nullptr, nullptr,
                  nullptr, &FbarFac(), &Htensor());
            else if (eastype_ != soh8p_easnone)
              CondensePlasticity<zerospin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  nullptr, detJ_w, gp, gp_temp, params, force, stiffmatrix, &M_eas(), &Kda, &dHda);
            else
              CondensePlasticity<zerospin>(DefgrdMod(), DeltaLp(), Bop(), &DerivShapeFunctionXYZ(),
                  nullptr, detJ_w, gp, gp_temp, params, force, stiffmatrix);
          }
        }  // plastic modifications
  }        // gp loop
  InvalidGpData();

  // Static condensation EAS --> stiff ********************************
  if (stiffmatrix != nullptr && !is_tangDis && eastype_ != soh8p_easnone)
  {
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(*KaaInv_);
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix kdakaai(numdofperelement_, neas_);
    switch (eastype_)
    {
      case soh8p_easfull:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
            0., kdakaai.A(), 1., Kda.A(), KaaInv_->A());
        if (stiffmatrix != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.A(), Kad_->A());
        if (force != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., force->A(), -1., kdakaai.A(), feas_->A());
        break;
      case soh8p_easmild:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
            0., kdakaai.A(), 1., Kda.A(), KaaInv_->A());
        if (stiffmatrix != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numdofperelement_>(
              1., stiffmatrix->A(), -1., kdakaai.A(), Kad_->A());
        if (force != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., force->A(), -1., kdakaai.A(), feas_->A());
        break;
      case soh8p_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }

    // TSI with EAS
    if (eval_tsi)
    {
      Epetra_SerialDenseVector dHdaKaai(neas_);
      switch (eastype_)
      {
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, nen_>(
              0., KdT_eas_->A(), -1., kdakaai.A(), KaT_->A());
          for (int gp = 0; gp < numgpt_; ++gp)
          {
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, 1,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
                0., dHdaKaai.A(), 1., dHda.at(gp).A(), KaaInv_->A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
                1., plmat->dHepDissDd(gp).A(), -1., Kad_->A(), dHdaKaai.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, 1,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
                1., &(plmat->HepDiss(gp)), -1., dHdaKaai.A(), feas_->A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, nen_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
                0., plmat->dHepDTeas()->at(gp).A(), -1., KaT_->A(), dHdaKaai.A());
          }
          break;
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, nen_>(
              0., KdT_eas_->A(), -1., kdakaai.A(), KaT_->A());
          for (int gp = 0; gp < numgpt_; ++gp)
          {
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, 1,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
                0., dHdaKaai.A(), 1., dHda.at(gp).A(), KaaInv_->A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
                1., plmat->dHepDissDd(gp).A(), -1., Kad_->A(), dHdaKaai.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double, 1,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
                1., &(plmat->HepDiss(gp)), -1., dHdaKaai.A(), feas_->A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, nen_,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
                0., plmat->dHepDTeas()->at(gp).A(), -1., KaT_->A(), dHdaKaai.A());
          }
          break;
        case soh8p_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }
  }

  InvalidEleData();
  return;
}

/*----------------------------------------------------------------------*
 | coupling term k_dT in monolithic TSI                     seitz 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::nln_kdT_tsi(
    CORE::LINALG::Matrix<numdofperelement_, nen_>* k_dT, Teuchos::ParameterList& params)
{
  if (k_dT == nullptr) return;

  // shape functions
  CORE::LINALG::Matrix<nen_, 1> shapefunct(false);

  for (int gp = 0; gp < numgpt_; gp++)
  {
    // shape functions
    CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
    // update linear coupling matrix K_dT
    k_dT->MultiplyNT(1., (*dFintdT_)[gp], shapefunct, 1.);
  }

  // EAS part
  if (eastype_ != soh8p_easnone)
  {
    if (KdT_eas_ == Teuchos::null)
      dserror("for TSI with EAS the block KdT_eas_ should be acessible here");
    k_dT->Update(1., *KdT_eas_, 1.);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  condense plastic degrees of freedom                     seitz 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <int spintype>
void DRT::ELEMENTS::So3_Plast<distype>::CondensePlasticity(
    const CORE::LINALG::Matrix<nsd_, nsd_>& defgrd, const CORE::LINALG::Matrix<nsd_, nsd_>& deltaLp,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>& bop,
    const CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ, const CORE::LINALG::Matrix<numstr_, 1>* RCG,
    const double detJ_w, const int gp, const double temp, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<numdofperelement_, 1>* force,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>* stiffmatrix,
    const Epetra_SerialDenseMatrix* M, Epetra_SerialDenseMatrix* Kda,
    std::vector<Epetra_SerialDenseVector>* dHda, const double* f_bar_factor,
    const CORE::LINALG::Matrix<numdofperelement_, 1>* htensor)
{
  bool eval_tsi = tsi_ && (temp != -1.e12);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    plmat = static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // temporary Epetra matrix for matrix-matrix-matrix products
  Epetra_SerialDenseMatrix tmp;

  // Nitsche contact
  CORE::LINALG::Matrix<numstr_, 1>* cauchy_ptr = nullptr;
  CORE::LINALG::Matrix<numstr_, spintype + 1> d_cauchy_ddp;
  CORE::LINALG::Matrix<numstr_, spintype + 1>* d_cauchy_ddp_ptr = nullptr;
  CORE::LINALG::Matrix<numstr_, numstr_> d_cauchy_dC;
  CORE::LINALG::Matrix<numstr_, numstr_>* d_cauchy_dC_ptr = nullptr;
  CORE::LINALG::Matrix<numstr_, 9> d_cauchy_dF;
  CORE::LINALG::Matrix<numstr_, 9>* d_cauchy_dF_ptr = nullptr;
  CORE::LINALG::Matrix<numstr_, 1>
      d_cauchy_dT;  // todo: continue with this one; plmat->EvaluatePlast(...) should give you this;
  // however not tested yet
  CORE::LINALG::Matrix<numstr_, 1>* d_cauchy_dT_ptr = nullptr;
  if (is_nitsche_contact_)
  {
    cauchy_ptr = &(cauchy_.at(gp));
    d_cauchy_ddp_ptr = &d_cauchy_ddp;
    d_cauchy_dC_ptr = &d_cauchy_dC;
    d_cauchy_dF_ptr = &d_cauchy_dF;
    if (eval_tsi) d_cauchy_dT_ptr = &d_cauchy_dT;
  }

  // second material call ****************************************************
  CORE::LINALG::Matrix<numstr_, spintype + 1> dpk2ddp;
  CORE::LINALG::Matrix<spintype + 1, 1> ncp;
  CORE::LINALG::Matrix<spintype + 1, spintype + 1> dncpddp;
  CORE::LINALG::Matrix<spintype + 1, numstr_> dncpdc;
  CORE::LINALG::Matrix<spintype + 1, 1> dncpdT;
  CORE::LINALG::Matrix<numstr_, 1> dHdC;
  CORE::LINALG::Matrix<spintype + 1, 1> dHdLp;
  bool active = false;
  bool elast = false;
  bool as_converged = true;
  if (!eval_tsi)
    plmat->EvaluatePlast(&defgrd, &deltaLp, nullptr, params, &dpk2ddp, &ncp, &dncpdc, &dncpddp,
        &active, &elast, &as_converged, gp, nullptr, nullptr, nullptr,
        StrParamsInterface().GetDeltaTime(), Id(), cauchy_ptr, d_cauchy_ddp_ptr, d_cauchy_dC_ptr,
        d_cauchy_dF_ptr, d_cauchy_dT_ptr);
  else
    plmat->EvaluatePlast(&defgrd, &deltaLp, &temp, params, &dpk2ddp, &ncp, &dncpdc, &dncpddp,
        &active, &elast, &as_converged, gp, &dncpdT, &dHdC, &dHdLp,
        StrParamsInterface().GetDeltaTime(), Id(), cauchy_ptr, d_cauchy_ddp_ptr, d_cauchy_dC_ptr,
        d_cauchy_dF_ptr, d_cauchy_dT_ptr);
  // *************************************************************************

  // Simple matrix do delete the linear dependent row in Voigt-notation.
  // Since all entries are deviatoric, we drop the third entry on the main diagonal.
  CORE::LINALG::Matrix<spintype + 1, spintype> voigt_red(true);
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
      dserror("Don't know what to do with plastic spin type %d", spintype);
      break;
  }

  // We have to adapt the derivatives to account for the deviatoric structure
  // of delta D^p.
  CORE::LINALG::Matrix<spintype + 1, spintype> dDpdbeta;
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
      dserror("Don't know what to do with plastic spin type %d", spintype);
      break;
  }

  // derivative of internal force w.r.t. beta
  CORE::LINALG::Matrix<numdofperelement_, spintype> kdbeta;
  // derivative of pk2 w.r.t. beta
  CORE::LINALG::Matrix<numstr_, spintype> dpk2db;
  dpk2db.Multiply(dpk2ddp, dDpdbeta);
  if (fbar_)
    kdbeta.MultiplyTN(detJ_w / (*f_bar_factor), bop, dpk2db);
  else
    kdbeta.MultiplyTN(detJ_w, bop, dpk2db);

  CORE::LINALG::Matrix<numstr_, spintype> d_cauchy_db;
  if (is_nitsche_contact_)
  {
    if (N_XYZ == nullptr) dserror("shape derivative not provided");

    CORE::LINALG::Matrix<9, numdofperelement_> dFdd;
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
      if (RCG == nullptr) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
      CORE::LINALG::Matrix<6, 1> tmp61;
      tmp61.Multiply(.5, d_cauchy_dC, (*RCG));
      cauchy_deriv_.at(gp).MultiplyNT(
          (*f_bar_factor) * (*f_bar_factor) * 2. / 3., tmp61, *htensor, 1.);
      cauchy_deriv_.at(gp).Multiply((*f_bar_factor) * (*f_bar_factor), d_cauchy_dC, bop, 1.);

      CORE::LINALG::Matrix<9, 1> f;
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
    cauchy_deriv_T_.at(gp).MultiplyNT(1., d_cauchy_dT, ShapeFunction(), 1.);
  }

  // EAS matrix block
  Epetra_SerialDenseMatrix Kab(neas_, spintype);
  switch (eastype_)
  {
    case soh8p_easnone:
      // do nothing
      break;
    case soh8p_easmild:
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, spintype>(
          1., Kab.A(), detJ_w, M->A(), dpk2db.A());
      break;
    case soh8p_easfull:
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, spintype>(
          1., Kab.A(), detJ_w, M->A(), dpk2db.A());
      break;
    case soh8p_eassosh8:
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numstr_, spintype>(
          1., Kab.A(), detJ_w, M->A(), dpk2db.A());
      break;
    case soh18p_eassosh18:
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
          PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, numstr_, spintype>(
          1., Kab.A(), detJ_w, M->A(), dpk2db.A());
      break;
    default:
      dserror("Don't know what to do with EAS type %d", eastype_);
      break;
  }

  // gauss points that require a full condensation
  if (!elast)
  {
    // apply chain rule for kbb block
    CORE::LINALG::Matrix<spintype + 1, spintype> dNCPdb;
    dNCPdb.Multiply(dncpddp, dDpdbeta);
    CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype + 1, spintype>(
        0., KbbInv_[gp].A(), 1., voigt_red.A(), dNCPdb.A());

    // apply chain rule for kbd block
    CORE::LINALG::Matrix<spintype + 1, numdofperelement_> dNCPdd;
    if (fbar_)
    {
      if (RCG == nullptr) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
      CORE::LINALG::Matrix<spintype + 1, 1> tmp61;
      tmp61.Multiply(.5, dncpdc, (*RCG));
      dNCPdd.MultiplyNT((*f_bar_factor) * (*f_bar_factor) * 2. / 3., tmp61, *htensor, 0.);
      dNCPdd.Multiply((*f_bar_factor) * (*f_bar_factor), dncpdc, bop, 1.);
    }
    else
      dNCPdd.Multiply(dncpdc, bop);
    CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype + 1, numdofperelement_>(
        0., Kbd_[gp].A(), 1., voigt_red.A(), dNCPdd.A());

    // EAS block kba
    if (eastype_ != soh8p_easnone)
    {
      Kba_->at(gp).Shape(spintype, neas_);
      tmp.Shape(spintype + 1, neas_);
      switch (eastype_)
      {
        case soh8p_easnone:
          break;
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              0., tmp.A(), 1., dncpdc.A(), M->A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype + 1,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              0., Kba_->at(gp).A(), 1., voigt_red.A(), tmp.A());
          break;
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              0., tmp.A(), 1., dncpdc.A(), M->A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              0., Kba_->at(gp).A(), 1., voigt_red.A(), tmp.A());
          break;
        case soh8p_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
              0., tmp.A(), 1., dncpdc.A(), M->A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
              0., Kba_->at(gp).A(), 1., voigt_red.A(), tmp.A());
          break;
        case soh18p_eassosh18:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype + 1, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas>(
              0., tmp.A(), 1., dncpdc.A(), M->A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas>(
              0., Kba_->at(gp).A(), 1., voigt_red.A(), tmp.A());
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    // residual
    CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype + 1, 1>(
        0., fbeta_[gp].A(), 1., voigt_red.A(), ncp.A());

    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
    // inverse matrix block [k_beta beta]_ij
    Epetra_SerialDenseSolver solve_for_kbbinv;
    solve_for_kbbinv.SetMatrix(KbbInv_[gp]);
    int err = solve_for_kbbinv.Invert();
    if (err != 0)
    {
      // check, if errors are tolerated or should throw a dserror
      bool error_tol = false;
      if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error", true);
        return;
      }
      else
        dserror("Inversion of Kbb failed");
    }
    // temporary  Kdb.Kbb^-1
    CORE::LINALG::Matrix<numdofperelement_, spintype> KdbKbb;
    CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype, spintype>(
        0., KdbKbb.A(), 1., kdbeta.A(), KbbInv_[gp].A());

    // "plastic displacement stiffness"
    // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
    if (stiffmatrix != nullptr)
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype,
          numdofperelement_>(1., stiffmatrix->A(), -1., KdbKbb.A(), Kbd_[gp].A());

    // "plastic internal force"
    // plFint = [K_db.K_bb^-1].f_b
    if (force != nullptr)
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype, 1>(
          1., force->A(), -1., KdbKbb.A(), fbeta_[gp].A());

    // TSI
    if (eval_tsi)
    {
      // thermal derivative
      (*KbT_).at(gp).Size(spintype);
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype + 1, 1>(
          0., (*KbT_)[gp].A(), 1., voigt_red.A(), dncpdT.A());

      // condense to K_dT
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype, 1>(
          1., (*dFintdT_)[gp].A(), -1., KdbKbb.A(), (*KbT_)[gp].A());

      // Plastic heating and dissipation dC
      CORE::LINALG::Matrix<numdofperelement_, 1> dHepDissDd;
      if (fbar_)
      {
        if (RCG == nullptr) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
        CORE::LINALG::Matrix<1, 1> tmp11;
        tmp11.MultiplyTN(.5, dHdC, (*RCG));
        dHepDissDd.Multiply((*f_bar_factor) * (*f_bar_factor) * 2. / 3., *htensor, tmp11, 0.);
        dHepDissDd.MultiplyTN((*f_bar_factor) * (*f_bar_factor), bop, dHdC, 1.);
      }
      else
        dHepDissDd.MultiplyTN(bop, dHdC);

      // store in material
      CORE::LINALG::DENSEFUNCTIONS::update<double, numdofperelement_, 1>(
          1., plmat->dHepDissDd(gp).A(), 1., dHepDissDd.A());

      // Plastic heating and dissipation dbeta
      CORE::LINALG::Matrix<spintype, 1> dHepDissDbeta;
      if (spintype != zerospin)
        dserror("no TSI with plastic spin yet");
      else
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, zerospin, zerospin + 1, 1>(
            0., dHepDissDbeta.A(), 1., dDpdbeta.A(), dHdLp.A());

      // condense the heating terms
      CORE::LINALG::Matrix<spintype, 1> dHdbKbbi;
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, spintype, spintype, 1>(
          0., dHdbKbbi.A(), 1., KbbInv_[gp].A(), dHepDissDbeta.A());
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, 1, spintype, 1>(
          1., &(plmat->HepDiss(gp)), -1., dHdbKbbi.A(), fbeta_[gp].A());
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_, spintype, 1>(
          1., plmat->dHepDissDd(gp).A(), -1., Kbd_[gp].A(), dHdbKbbi.A());
      CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, 1, spintype, 1>(
          1., &(plmat->dHepDT(gp)), -1., (*KbT_)[gp].A(), dHdbKbbi.A());

      // TSI with EAS
      if (eastype_ != soh8p_easnone)
      {
        // error checks
        if (dHda == nullptr)
          dserror("dHda is nullptr pointer");
        else if ((int)dHda->size() != numgpt_)
          dserror("dHda has wrong size");

        switch (eastype_)
        {
          case soh8p_easmild:
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
                1., dHda->at(gp).A(), 1., M->A(), dHdC.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
                1., dHda->at(gp).A(), -1., Kba_->at(gp).A(), dHdbKbbi.A());

            break;
          case soh8p_easfull:
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
                1., dHda->at(gp).A(), 1., M->A(), dHdC.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
                1., dHda->at(gp).A(), -1., Kba_->at(gp).A(), dHdbKbbi.A());
            break;
          case soh8p_eassosh8:
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numstr_, 1>(
                1., dHda->at(gp).A(), 1., M->A(), dHdC.A());
            CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
                1., dHda->at(gp).A(), -1., Kba_->at(gp).A(), dHdbKbbi.A());
            break;
          case soh8p_easnone:
            break;
          default:
            dserror("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }
    }  // TSI

    if (eastype_ != soh8p_easnone)
    {
      // condense plasticity into EAS matrix blocks
      tmp.Shape(neas_, spintype);
      switch (eastype_)
      {
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              1., Kda->A(), -1., KdbKbb.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype, spintype>(
              0., tmp.A(), 1., Kab.A(), KbbInv_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype,
              numdofperelement_>(1., Kad_->A(), -1., tmp.A(), Kbd_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              1., KaaInv_->A(), -1., tmp.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
              1., feas_->A(), -1., tmp.A(), fbeta_[gp].A());
          if (eval_tsi)
          {
            Epetra_SerialDenseMatrix kbTm(spintype, nen_);
            CORE::LINALG::Matrix<nen_, 1> shapefunct;
            CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
            CORE::LINALG::DENSEFUNCTIONS::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.A(), 1., KbT_->at(gp).A(), shapefunct.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype, nen_>(
                1., KaT_->A(), -1., tmp.A(), kbTm.A());
          }
          break;
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              1., Kda->A(), -1., KdbKbb.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype, spintype>(
              0., tmp.A(), 1., Kab.A(), KbbInv_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype,
              numdofperelement_>(1., Kad_->A(), -1., tmp.A(), Kbd_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              1., KaaInv_->A(), -1., tmp.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
              1., feas_->A(), -1., tmp.A(), fbeta_[gp].A());
          if (eval_tsi)
          {
            Epetra_SerialDenseMatrix kbTm(spintype, nen_);
            CORE::LINALG::Matrix<nen_, 1> shapefunct;
            CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
            CORE::LINALG::DENSEFUNCTIONS::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.A(), 1., KbT_->at(gp).A(), shapefunct.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype, nen_>(
                1., KaT_->A(), -1., tmp.A(), kbTm.A());
          }
          break;
        case soh8p_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
              1., Kda->A(), -1., KdbKbb.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype, spintype>(
              0., tmp.A(), 1., Kab.A(), KbbInv_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype,
              numdofperelement_>(1., Kad_->A(), -1., tmp.A(), Kbd_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
              1., KaaInv_->A(), -1., tmp.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
              1., feas_->A(), -1., tmp.A(), fbeta_[gp].A());
          if (eval_tsi)
          {
            Epetra_SerialDenseMatrix kbTm(spintype, nen_);
            CORE::LINALG::Matrix<nen_, 1> shapefunct;
            CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
            CORE::LINALG::DENSEFUNCTIONS::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.A(), 1., KbT_->at(gp).A(), shapefunct.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype, nen_>(
                1., KaT_->A(), -1., tmp.A(), kbTm.A());
          }
          break;
        case soh18p_eassosh18:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas>(
              1., Kda->A(), -1., KdbKbb.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype, spintype>(
              0., tmp.A(), 1., Kab.A(), KbbInv_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype,
              numdofperelement_>(1., Kad_->A(), -1., tmp.A(), Kbd_[gp].A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas>(
              1., KaaInv_->A(), -1., tmp.A(), Kba_->at(gp).A());
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype, 1>(
              1., feas_->A(), -1., tmp.A(), fbeta_[gp].A());
          if (eval_tsi)
          {
            Epetra_SerialDenseMatrix kbTm(spintype, nen_);
            CORE::LINALG::Matrix<nen_, 1> shapefunct;
            CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
            CORE::LINALG::DENSEFUNCTIONS::multiplyNT<double, spintype, 1, nen_>(
                0., kbTm.A(), 1., KbT_->at(gp).A(), shapefunct.A());
            CORE::LINALG::DENSEFUNCTIONS::multiply<double,
                PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype, nen_>(
                1., KaT_->A(), -1., tmp.A(), kbTm.A());
          }
          break;
        case soh8p_easnone:
          // do nothing
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    if (is_nitsche_contact_)
    {
      CORE::LINALG::Matrix<numstr_, spintype> tmp;
      tmp.Multiply(d_cauchy_db, CORE::LINALG::Matrix<spintype, spintype>(KbbInv_.at(gp).A(), true));
      cauchy_.at(gp).Multiply(
          -1., tmp, CORE::LINALG::Matrix<spintype, 1>(fbeta_.at(gp).A(), true), 1.);
      cauchy_deriv_.at(gp).Multiply(
          -1., tmp, CORE::LINALG::Matrix<spintype, numdofperelement_>(Kbd_.at(gp).A(), true), 1.);

      if (d_cauchy_dT_ptr)
      {
        CORE::LINALG::Matrix<6, 1> bla;
        bla.Multiply(-1., tmp, CORE::LINALG::Matrix<spintype, 1>(KbT_->at(gp).A(), true), 1.);
        cauchy_deriv_T_.at(gp).MultiplyNT(1., bla, ShapeFunction(), 1.);
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
    if (dDp_last_iter_[gp].NormInf() > 0.)
    {
      if (force != nullptr)
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, numdofperelement_, spintype, 1>(
            1., force->A(), -1., kdbeta.A(), dDp_last_iter_[gp].A());

      switch (eastype_)
      {
        case soh8p_easnone:
          break;
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, spintype, 1>(
              1., feas_->A(), -1., Kab.A(), dDp_last_iter_[gp].A());
          break;
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, spintype, 1>(
              1., feas_->A(), -1., Kab.A(), dDp_last_iter_[gp].A());
          break;
        case soh8p_eassosh8:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, spintype, 1>(
              1., feas_->A(), -1., Kab.A(), dDp_last_iter_[gp].A());
          break;
        case soh18p_eassosh18:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, spintype, 1>(
              1., feas_->A(), -1., Kab.A(), dDp_last_iter_[gp].A());
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    KbbInv_[gp].Scale(0.);
    for (int i = 0; i < KbbInv_[gp].RowDim(); ++i) KbbInv_[gp](i, i) = 1. / (plmat->cpl());
    fbeta_[gp].Scale(0.);
    fbeta_[gp] = dDp_last_iter_[gp];
    fbeta_[gp].Scale(plmat->cpl());
    Kbd_[gp].Scale(0.);
    if (eastype_ != soh8p_easnone) Kba_->at(gp).Shape(spintype, neas_);
    if (KbT_ != Teuchos::null) KbT_->at(gp).Scale(0.);
  }

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::BuildDeltaLp(const int gp)
{
  // current plastic flow increment
  SetDeltaLp()(0, 0) = dDp_last_iter_[gp](0);
  SetDeltaLp()(1, 1) = dDp_last_iter_[gp](1);
  SetDeltaLp()(2, 2) = -1.0 * (dDp_last_iter_[gp](0) + dDp_last_iter_[gp](1));
  SetDeltaLp()(0, 1) = dDp_last_iter_[gp](2);
  SetDeltaLp()(1, 0) = dDp_last_iter_[gp](2);
  SetDeltaLp()(1, 2) = dDp_last_iter_[gp](3);
  SetDeltaLp()(2, 1) = dDp_last_iter_[gp](3);
  SetDeltaLp()(0, 2) = dDp_last_iter_[gp](4);
  SetDeltaLp()(2, 0) = dDp_last_iter_[gp](4);
  if (HavePlasticSpin())
  {
    SetDeltaLp()(0, 1) += dDp_last_iter_[gp](5);
    SetDeltaLp()(1, 0) -= dDp_last_iter_[gp](5);
    SetDeltaLp()(1, 2) += dDp_last_iter_[gp](6);
    SetDeltaLp()(2, 1) -= dDp_last_iter_[gp](6);
    SetDeltaLp()(0, 2) += dDp_last_iter_[gp](7);
    SetDeltaLp()(2, 0) -= dDp_last_iter_[gp](7);
  }
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::RecoverPlasticityAndEAS(
    const CORE::LINALG::Matrix<numdofperelement_, 1>* res_d,
    const CORE::LINALG::Matrix<nen_, 1>* res_T)
{
  if (StrParamsInterface().IsDefaultStep())
  {
    if (eastype_ != soh8p_easnone) RecoverEAS(res_d, res_T);


    if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    {
      CORE::LINALG::Matrix<nen_, 1> shapefunct;
      double res_t = -1.;
      double* res_t_ptr;
      for (int gp = 0; gp < numgpt_; ++gp)
      {
        if (res_T)
        {
          CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
          res_t = shapefunct.Dot(*res_T);
          res_t_ptr = &res_t;
        }
        else
          res_t_ptr = nullptr;

        // recover plastic variables
        if (HavePlasticSpin())
          RecoverPlasticity<plspin>(res_d, gp, res_t_ptr);
        else
          RecoverPlasticity<zerospin>(res_d, gp, res_t_ptr);
      }
    }
  }
  else
  {
    //    const double new_step_length   = StrParamsInterface().GetStepLength();
    //
    //    if (eastype_!=soh8p_easnone)
    //      ReduceEasStep(new_step_length,old_step_length_);
    //
    //    if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    //      for (int gp=0;gp<numgpt_;++gp)
    //        ReducePlasticityStep(new_step_length,old_step_length_,gp);
    //
    //    old_step_length_=new_step_length;
  }
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::RecoverEAS(
    const CORE::LINALG::Matrix<numdofperelement_, 1>* res_d,
    const CORE::LINALG::Matrix<nen_, 1>* res_T)
{
  if (eastype_ == soh8p_easnone) return;

  const double step_length = StrParamsInterface().GetStepLength();

  // first, store the eas state of the previous accepted Newton step
  StrParamsInterface().SumIntoMyPreviousSolNorm(
      NOX::NLN::StatusTest::quantity_eas, neas_, alpha_eas_->A(), Owner());

  if (StrParamsInterface().IsDefaultStep()) switch (eastype_)
    {
      case soh8p_easmild:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numdofperelement_, 1>(
            1.0, feas_->A(), 1.0, Kad_->A(), res_d->A());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, nen_, 1>(
              1., feas_->A(), 1., KaT_->A(), res_T->A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        CORE::LINALG::DENSEFUNCTIONS::update<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_easfull:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numdofperelement_, 1>(
            1.0, feas_->A(), 1.0, Kad_->A(), res_d->A());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, nen_, 1>(
              1., feas_->A(), 1., KaT_->A(), res_T->A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        CORE::LINALG::DENSEFUNCTIONS::update<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_eassosh8:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, numdofperelement_, 1>(
            1.0, feas_->A(), 1.0, Kad_->A(), res_d->A());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, nen_, 1>(
              1., feas_->A(), 1., KaT_->A(), res_T->A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        CORE::LINALG::DENSEFUNCTIONS::update<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh18p_eassosh18:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, numdofperelement_, 1>(
            1.0, feas_->A(), 1.0, Kad_->A(), res_d->A());
        if (KaT_ != Teuchos::null && res_T != nullptr)
          CORE::LINALG::DENSEFUNCTIONS::multiply<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, nen_, 1>(
              1., feas_->A(), 1., KaT_->A(), res_T->A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            0.0, *alpha_eas_inc_, -1.0, *KaaInv_, *feas_);
        CORE::LINALG::DENSEFUNCTIONS::update<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            1.0, *alpha_eas_, 1.0, *alpha_eas_inc_);
        break;
      case soh8p_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  else
    dserror("no line search implemented yet");

  StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, neas_,
      alpha_eas_inc_->A(), alpha_eas_->A(), step_length, Owner());

  return;
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <int spintype>
void DRT::ELEMENTS::So3_Plast<distype>::RecoverPlasticity(
    const CORE::LINALG::Matrix<numdofperelement_, 1>* res_d, const int gp, const double* res_t)
{
  const double step_length = StrParamsInterface().GetStepLength();

  if (StrParamsInterface().IsDefaultStep() == false) dserror("no line search implemented yet");

  // first, store the state of the previous accepted Newton step
  StrParamsInterface().SumIntoMyPreviousSolNorm(
      NOX::NLN::StatusTest::quantity_plasticity, spintype, dDp_last_iter_[gp].A(), Owner());

  // temporary Epetra matrix
  Epetra_SerialDenseVector tmp_v(spintype);
  Epetra_SerialDenseMatrix tmp_m(spintype, numdofperelement_);

  // first part
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype, 1>(
      0., dDp_inc_[gp].A(), -1., KbbInv_[gp].A(), fbeta_[gp].A());

  // second part
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype, numdofperelement_>(
      0., tmp_m.A(), 1., KbbInv_[gp].A(), Kbd_[gp].A());
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, numdofperelement_, 1>(
      1., dDp_inc_[gp].A(), -1., tmp_m.A(), res_d->A());

  // thermal part
  if (KbT_ != Teuchos::null)
    if (res_t)
    {
      CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype, 1>(
          1., dDp_inc_[gp].A(), -1. * (*res_t), KbbInv_[gp].A(), (*KbT_)[gp].A());
    }

  // EAS part
  if (eastype_ != soh8p_easnone)
  {
    tmp_m.Shape(spintype, neas_);
    switch (eastype_)
    {
      case soh8p_easmild:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
            0., tmp_m.A(), 1., KbbInv_[gp].A(), Kba_->at(gp).A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
            1., dDp_inc_[gp].A(), -1., tmp_m.A(), alpha_eas_inc_->A());
        break;
      case soh8p_easfull:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
            0., tmp_m.A(), 1., KbbInv_[gp].A(), Kba_->at(gp).A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
            1., dDp_inc_[gp].A(), -1., tmp_m.A(), alpha_eas_inc_->A());
        break;
      case soh8p_eassosh8:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas>(
            0., tmp_m.A(), 1., KbbInv_[gp].A(), Kba_->at(gp).A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>::neas, 1>(
            1., dDp_inc_[gp].A(), -1., tmp_m.A(), alpha_eas_inc_->A());
        break;
      case soh18p_eassosh18:
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas>(
            0., tmp_m.A(), 1., KbbInv_[gp].A(), Kba_->at(gp).A());
        CORE::LINALG::DENSEFUNCTIONS::multiply<double, spintype,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>::neas, 1>(
            1., dDp_inc_[gp].A(), -1., tmp_m.A(), alpha_eas_inc_->A());
        break;
      case soh8p_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }  // EAS part

  CORE::LINALG::DENSEFUNCTIONS::update<double, spintype, 1>(
      1., dDp_last_iter_[gp], 1., dDp_inc_[gp]);

  StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_plasticity, spintype,
      dDp_inc_[gp].A(), dDp_inc_[gp].A(), step_length, Owner());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::ReduceEasStep(
    const double new_step_length, const double old_step_length)
{
  if (eastype_ == soh8p_easnone) return;

  alpha_eas_->Update(-1., *alpha_eas_inc_, 1.);
  alpha_eas_inc_->Scale(new_step_length / old_step_length);
  alpha_eas_->Update(+1., *alpha_eas_inc_, 1.);

  StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, neas_,
      alpha_eas_inc_->A(), alpha_eas_->A(), new_step_length, Owner());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::ReducePlasticityStep(
    const double new_step_length, const double old_step_length, const int gp)
{
  dDp_last_iter_[gp].Update(-1., dDp_inc_[gp], 1.);
  dDp_inc_[gp].Scale(new_step_length / old_step_length);
  dDp_last_iter_[gp].Update(+1., dDp_inc_[gp], 1.);

  StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_plasticity, plspintype_,
      dDp_inc_[gp].A(), dDp_inc_[gp].A(), new_step_length, Owner());
}

/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 07/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::UpdatePlasticDeformation_nln(PlSpinType spintype)
{
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
  {
    // loop over all Gauss points
    for (int gp = 0; gp < numgpt_; gp++)
    {
      BuildDeltaLp(gp);
      static_cast<MAT::PlasticElastHyper*>(Material().get())->UpdateGP(gp, &DeltaLp());

      KbbInv_[gp].Scale(0.);
      Kbd_[gp].Scale(0.);
      fbeta_[gp].Scale(0.);
      if (tsi_) (*KbT_)[gp].Scale(0.);
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
    Kad_->Scale(0.);
    KaaInv_->Scale(0.);
    feas_->Scale(0.);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculate internal energy of the element (private)                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::So3_Plast<distype>::CalcIntEnergy(
    std::vector<double>& disp,       // current displacements
    std::vector<double>& temp,       // current temperatuere
    Teuchos::ParameterList& params)  // strain output option
{
  InvalidEleData();

  double energy = 0.;

  std::vector<double> bla;
  FillPositionArrays(disp, bla, temp);

  if (fbar_ || eastype_ != soh8p_easnone)
    EvaluateCenter();  // deformation gradient at centroid of element

  if (eastype_ != soh8p_easnone) EasSetup();

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    plmat = static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("elastic strain energy in so3plast elements only for plastic material");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    InvalidGpData();

    // shape functions (shapefunct) and their first derivatives (deriv)
    EvaluateShape(xsi_[gp]);
    EvaluateShapeDeriv(xsi_[gp]);

    Kinematics(gp);

    double gp_temp = ShapeFunction().Dot(Temp());

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8p_easnone)
    {
      EasShape(gp);
      EasEnhanceStrains();
    }
    // calculate modified deformation gradient
    else if (fbar_)
      SetupFbarGp();
    else
      SetDefgrdMod() = Defgrd();

    const double psi = plmat->StrainEnergyTSI(DefgrdMod(), gp, Id(), gp_temp);

    const double detJ_w = DetJ() * wgt_[gp];
    energy += detJ_w * psi;

  }  // gp loop

  return energy;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyNDirAndDerivativesAtXiElast(
    const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Epetra_SerialDenseMatrix* d_cauchyndir_dd,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd2, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd_ddir, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Epetra_SerialDenseMatrix* d_cauchyndir_dT, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dT)
{
  if (distype == DRT::Element::nurbs27) GetNurbsEleInfo();

  if (fbar_ || eastype_ != soh8p_easnone)
    dserror("cauchy stress not available for fbar or eas elements");

  if (temp || d_cauchyndir_dT || d2_cauchyndir_dd_dT)
    if (!temp || !d_cauchyndir_dT || !d2_cauchyndir_dd_dT)
      dserror("inconsistent temperature dependency input");
  if (temp && Material()->MaterialType() != INPAR::MAT::m_plelasthyper)
  {
    dserror(
        "thermo-mechanical Nitsche contact only with PlasticElastHyper"
        "\nIf you want to do elasticity, set a negative yield stress ;)");
  }

  auto* plmat = dynamic_cast<MAT::PlasticElastHyper*>(Material().get());

  cauchy_n_dir = 0.0;

  static CORE::LINALG::Matrix<nen_, nsd_> xrefe(true);  // reference coord. of element
  static CORE::LINALG::Matrix<nen_, nsd_> xcurr(true);  // current  coord. of element
  static CORE::LINALG::Matrix<nen_, 1> ele_temp(true);
  xrefe.Clear();
  xcurr.Clear();
  ele_temp.Clear();
  DRT::Node** nodes = Nodes();

  for (int i = 0; i < nen_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int d = 0; d < nsd_; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * nsd_ + d];
    }
    if (temp)
      if (!temp->empty()) ele_temp(i) = temp->at(i);
  }

  EvaluateShape(xi);
  EvaluateShapeDeriv(xi);

  const double gp_temp = ele_temp.Dot(ShapeFunction());
  double d_cauchyndir_dT_gp = 0.0;
  static CORE::LINALG::Matrix<nsd_, 1> d_T_dxi(true);
  d_T_dxi.Multiply(1.0, DerivShapeFunction(), ele_temp, 0.0);

  static CORE::LINALG::Matrix<nsd_, nen_> N_XYZ(true);
  static CORE::LINALG::Matrix<nsd_, nsd_> invJ(true);
  invJ.Multiply(1.0, DerivShapeFunction(), xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, DerivShapeFunction(), 0.0);
  static CORE::LINALG::Matrix<nsd_, nsd_> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static CORE::LINALG::Matrix<9, numdofperelement_> d_F_dd(true);
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

  static CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF(true);
  static CORE::LINALG::Matrix<9, 9> d2_cauchyndir_dF2(true);
  static CORE::LINALG::Matrix<9, nsd_> d2_cauchyndir_dF_dn(true);
  static CORE::LINALG::Matrix<9, nsd_> d2_cauchyndir_dF_ddir(true);
  static CORE::LINALG::Matrix<9, 1> d2_cauchyndir_dF_dT(true);

  if (plmat && temp)
  {
    plmat->EvaluateCauchyNDirAndDerivatives(defgrd, n, dir, cauchy_n_dir, d_cauchyndir_dn,
        d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2, &d2_cauchyndir_dF_dn,
        &d2_cauchyndir_dF_ddir, -1, Id(), nullptr, &gp_temp, &d_cauchyndir_dT_gp,
        &d2_cauchyndir_dF_dT);
  }
  else
  {
    SolidMaterial()->EvaluateCauchyNDirAndDerivatives(defgrd, n, dir, cauchy_n_dir, d_cauchyndir_dn,
        d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2, &d2_cauchyndir_dF_dn,
        &d2_cauchyndir_dF_ddir, -1, Id(), nullptr, nullptr, nullptr, nullptr);
  }

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->Shape(numdofperelement_, 1);
    CORE::LINALG::Matrix<numdofperelement_, 1> d_cauchyndir_dd_mat(d_cauchyndir_dd->A(), true);
    d_cauchyndir_dd_mat.MultiplyTN(1.0, d_F_dd, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dT)
  {
    d2_cauchyndir_dd_dT->Shape(numdofperelement_, nen_);
    static CORE::LINALG::Matrix<numdofperelement_, 1> tmp(true);
    tmp.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dT, 0.0);
    CORE::LINALG::Matrix<numdofperelement_, nen_>(d2_cauchyndir_dd_dT->A(), true)
        .MultiplyNT(tmp, ShapeFunction());
  }

  if (d2_cauchyndir_dd_dn)
  {
    d2_cauchyndir_dd_dn->Shape(numdofperelement_, nsd_);
    CORE::LINALG::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_dn_mat(
        d2_cauchyndir_dd_dn->A(), true);
    d2_cauchyndir_dd_dn_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dn, 0.0);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->Shape(numdofperelement_, nsd_);
    CORE::LINALG::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_ddir_mat(
        d2_cauchyndir_dd_ddir->A(), true);
    d2_cauchyndir_dd_ddir_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_ddir, 0.0);
  }

  if (d_cauchyndir_dT)
  {
    d_cauchyndir_dT->Shape(nen_, 1);
    CORE::LINALG::Matrix<nen_, 1>(d_cauchyndir_dT->A(), true)
        .Update(d_cauchyndir_dT_gp, ShapeFunction(), 1.0);
  }

  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->Shape(numdofperelement_, numdofperelement_);
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_> d2_cauchyndir_dd2_mat(
        d2_cauchyndir_dd2->A(), true);
    static CORE::LINALG::Matrix<9, numdofperelement_> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd2_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF2_d_F_dd, 0.0);
  }

  // prepare evaluation of d_cauchyndir_dxi or d2_cauchyndir_dd_dxi
  static CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2, nen_>
      deriv2(true);
  static CORE::LINALG::Matrix<9, nsd_> d_F_dxi(true);
  deriv2.Clear();
  d_F_dxi.Clear();

  if (d_cauchyndir_dxi or d2_cauchyndir_dd_dxi)
  {
    if (distype == DRT::Element::nurbs27)
    {
      CORE::DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv_deriv2(
          SetShapeFunction(), SetDerivShapeFunction(), deriv2, xi, Knots(), Weights(), distype);
    }
    else
      CORE::DRT::UTILS::shape_function_deriv2<distype>(xi, deriv2);

    static CORE::LINALG::Matrix<nen_, nsd_> xXF(true);
    static CORE::LINALG::Matrix<nsd_, CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2>
        xXFsec(true);
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
    d2_cauchyndir_dd_dxi->Shape(numdofperelement_, nsd_);
    CORE::LINALG::Matrix<numdofperelement_, nsd_> d2_cauchyndir_dd_dxi_mat(
        d2_cauchyndir_dd_dxi->A(), true);

    static CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2, nsd_>
        Xsec(true);
    static CORE::LINALG::Matrix<nen_, 6> N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static CORE::LINALG::Matrix<9, numdofperelement_> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d2_cauchyndir_dF2_d_F_dd, d_F_dxi, 0.0);

    if (temp)
    {
      static CORE::LINALG::Matrix<9, nsd_> tmp(true);
      tmp.MultiplyNT(1.0, d2_cauchyndir_dF_dT, d_T_dxi, 0.0);
      d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d_F_dd, tmp, 1.0);
    }

    static CORE::LINALG::Matrix<9, nsd_ * numdofperelement_> d2_F_dxi_dd(true);
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
  InvalidEleData();
  InvalidGpData();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyNDirAndDerivativesAtXiPlast(
    const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Epetra_SerialDenseMatrix* d_cauchyndir_dd,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd2, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd_ddir, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Epetra_SerialDenseMatrix* d_cauchyndir_dT, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dT)
{
  if (distype != DRT::Element::hex8 || numgpt_ != 8) dserror("only for hex8 with 8 gp");
  if (Material()->MaterialType() != INPAR::MAT::m_plelasthyper)
    dserror("only PlasticElastHyper materials here");
  if ((int)cauchy_.size() != numgpt_ || (int)cauchy_deriv_.size() != numgpt_)
    dserror("have you evaluated the cauchy stress???");

  cauchy_n_dir = 0.0;
  if (d_cauchyndir_dxi) d_cauchyndir_dxi->Clear();
  if (d_cauchyndir_dT) d_cauchyndir_dT->Shape(nen_, 1);
  if (d2_cauchyndir_dd_dT) d2_cauchyndir_dd_dT->Shape(numdofperelement_, nen_);

  CORE::LINALG::Matrix<3, 3> n_dir_dir_n(true);
  n_dir_dir_n.MultiplyNT(.5, n, dir, 1.);
  n_dir_dir_n.MultiplyNT(.5, dir, n, 1.);
  CORE::LINALG::Matrix<6, 1> n_dir_dir_n_v;
  for (int i = 0; i < 3; ++i) n_dir_dir_n_v(i) = n_dir_dir_n(i, i);
  n_dir_dir_n_v(3) = n_dir_dir_n(0, 1) + n_dir_dir_n(1, 0);
  n_dir_dir_n_v(4) = n_dir_dir_n(2, 1) + n_dir_dir_n(1, 2);
  n_dir_dir_n_v(5) = n_dir_dir_n(0, 2) + n_dir_dir_n(2, 0);

  CORE::LINALG::Matrix<3, 1> xi_expol(xi);
  xi_expol.Scale(sqrt(3.));

  CORE::LINALG::Matrix<nen_, 1> shapefunct;
  CORE::LINALG::Matrix<nsd_, nen_> deriv;
  CORE::DRT::UTILS::shape_function<distype>(xi_expol, shapefunct);
  CORE::DRT::UTILS::shape_function_deriv1<distype>(xi_expol, deriv);

  CORE::LINALG::Matrix<numstr_, 1> cauchy_expol;
  CORE::LINALG::Matrix<numstr_, numdofperelement_> cauchy_deriv_expol;

  CORE::LINALG::Matrix<6, 1> tmp61;
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    cauchy_expol.Update(shapefunct(gp), cauchy_.at(gp), 1.);
    cauchy_deriv_expol.Update(shapefunct(gp), cauchy_deriv_.at(gp), 1.);
    if (d_cauchyndir_dxi)
      for (int d = 0; d < nsd_; ++d)
        (*d_cauchyndir_dxi)(d) += cauchy_.at(gp).Dot(n_dir_dir_n_v) * deriv(d, gp) * sqrt(3.);

    if (d_cauchyndir_dT)
      CORE::LINALG::Matrix<nen_, 1>(d_cauchyndir_dT->A(), true)
          .MultiplyTN(shapefunct(gp), cauchy_deriv_T_.at(gp), n_dir_dir_n_v, 1.);
  }

  cauchy_n_dir = cauchy_expol.Dot(n_dir_dir_n_v);

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->Reshape(numdofperelement_, 1);
    CORE::LINALG::Matrix<numdofperelement_, 1>(d_cauchyndir_dd->A(), true)
        .MultiplyTN(cauchy_deriv_expol, n_dir_dir_n_v);
  }
  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->Reshape(numdofperelement_, numdofperelement_);
    d2_cauchyndir_dd2->Scale(0.);
  }

  CORE::LINALG::Matrix<numstr_, nsd_> d_ndirdirn_v_dn, d_ndirdirn_v_dt;
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
    d2_cauchyndir_dd_dn->Reshape(numdofperelement_, nsd_);
    CORE::LINALG::Matrix<numdofperelement_, nsd_>(d2_cauchyndir_dd_dn->A(), true)
        .MultiplyTN(cauchy_deriv_expol, d_ndirdirn_v_dn);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->Reshape(numdofperelement_, nsd_);
    CORE::LINALG::Matrix<numdofperelement_, nsd_>(d2_cauchyndir_dd_ddir->A(), true)
        .MultiplyTN(cauchy_deriv_expol, d_ndirdirn_v_dt);
  }
  if (d2_cauchyndir_dd_dxi)
  {
    d2_cauchyndir_dd_dxi->Reshape(numdofperelement_, nsd_);
    d2_cauchyndir_dd_dxi->Scale(0.);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyNDirAndDerivativesAtXi(
    const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir,
    double& cauchy_n_dir, Epetra_SerialDenseMatrix* d_cauchyndir_dd,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd2, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dn,
    Epetra_SerialDenseMatrix* d2_cauchyndir_dd_ddir, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    Epetra_SerialDenseMatrix* d_cauchyndir_dT, Epetra_SerialDenseMatrix* d2_cauchyndir_dd_dT,
    const double* concentration, double* d_cauchyndir_dc)
{
  if (d_cauchyndir_dc != nullptr) dserror("Not implemented");

  bool elastic = true;
  auto* plmat = dynamic_cast<MAT::PlasticElastHyper*>(Material().get());
  if (!plmat)
    elastic = true;
  else
    elastic = plmat->AllElastic();

  if (!elastic)
  {
    GetCauchyNDirAndDerivativesAtXiPlast(xi, disp, n, dir, cauchy_n_dir, d_cauchyndir_dd,
        d2_cauchyndir_dd2, d2_cauchyndir_dd_dn, d2_cauchyndir_dd_ddir, d2_cauchyndir_dd_dxi,
        d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dxi, temp, d_cauchyndir_dT,
        d2_cauchyndir_dd_dT);
  }
  else
  {
    GetCauchyNDirAndDerivativesAtXiElast(xi, disp, n, dir, cauchy_n_dir, d_cauchyndir_dd,
        d2_cauchyndir_dd2, d2_cauchyndir_dd_dn, d2_cauchyndir_dd_ddir, d2_cauchyndir_dd_dxi,
        d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dxi, temp, d_cauchyndir_dT,
        d2_cauchyndir_dd_dT);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::OutputStrains(const int gp,
    const INPAR::STR::StrainType iostrain,                 // strain output option
    CORE::LINALG::Matrix<numgpt_post, numstr_>* elestrain  // strains at GP
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
      case INPAR::STR::strain_gl:
      {
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::Matrix<numstr_, 1> total_glstrain(false);
        total_glstrain(0) = 0.5 * (RCG()(0, 0) - 1.0);
        total_glstrain(1) = 0.5 * (RCG()(1, 1) - 1.0);
        total_glstrain(2) = 0.5 * (RCG()(2, 2) - 1.0);
        total_glstrain(3) = RCG()(0, 1);
        total_glstrain(4) = RCG()(1, 2);
        total_glstrain(5) = RCG()(2, 0);

        if (elestrain == nullptr) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = total_glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * total_glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == nullptr) dserror("strain data not available");

        // inverse of deformation gradient
        CORE::LINALG::Matrix<3, 3> invdefgrd;
        invdefgrd.Invert(Defgrd());

        static CORE::LINALG::Matrix<3, 3> tmp1;
        static CORE::LINALG::Matrix<3, 3> total_euler_almansi(true);
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
      case INPAR::STR::strain_none:
        break;
      default:
      {
        dserror("requested strain type not available");
        break;
      }
    }
  }
  // end of strain output **************************
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::OutputStress(const int gp,
    const INPAR::STR::StressType iostress,                 // strain output option
    CORE::LINALG::Matrix<numgpt_post, numstr_>* elestress  // strains at GP
)
{
  // return gp stresses
  switch (iostress)
  {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == nullptr) dserror("stress data not available");
      for (int i = 0; i < numstr_; ++i) (*elestress)(gp, i) = PK2()(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == nullptr) dserror("stress data not available");

      static CORE::LINALG::Matrix<3, 3> pkstress;
      pkstress(0, 0) = PK2()(0);
      pkstress(0, 1) = PK2()(3);
      pkstress(0, 2) = PK2()(5);
      pkstress(1, 0) = pkstress(0, 1);
      pkstress(1, 1) = PK2()(1);
      pkstress(1, 2) = PK2()(4);
      pkstress(2, 0) = pkstress(0, 2);
      pkstress(2, 1) = pkstress(1, 2);
      pkstress(2, 2) = PK2()(2);

      static CORE::LINALG::Matrix<3, 3> cauchystress;
      static CORE::LINALG::Matrix<3, 3> tmp1;
      tmp1.Multiply(1.0 / DetF(), Defgrd(), pkstress);
      cauchystress.MultiplyNT(tmp1, Defgrd());

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
    {
      dserror("requested stress type not available");
      break;
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::Kinematics(const int gp)
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
  SetInvJ().Multiply(DerivShapeFunction(), Xrefe());
  // xij_ = ds/dx
  SetDetJ() = SetInvJ().Invert();
  if (DetJ() < 1.0E-16) dserror("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", DetJ());

  /* get the inverse of the Jacobian matrix which looks like:
   **            [ x_,r  y_,r  z_,r ]^-1
   **     J^-1 = [ x_,s  y_,s  z_,s ]
   **            [ x_,t  y_,t  z_,t ]
   */
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  SetDerivShapeFunctionXYZ().Multiply(InvJ(), DerivShapeFunction());  // (6.21)

  // (material) deformation gradient
  // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
  SetDefgrd().MultiplyTT(Xcurr(), DerivShapeFunctionXYZ());

  // inverse deformation gradient and determinant
  SetDetF() = SetInvDefgrd().Invert(Defgrd());

  // calcualte total rcg
  SetRCG().MultiplyTN(Defgrd(), Defgrd());

  // total rcg in strain-like voigt notation
  for (int i = 0; i < 3; i++) SetRCGvec()(i) = RCG()(i, i);
  SetRCGvec()(3) = RCG()(0, 1) * 2.;
  SetRCGvec()(4) = RCG()(1, 2) * 2.;
  SetRCGvec()(5) = RCG()(0, 2) * 2.;

  // calculate nonlinear B-operator
  CalculateBop(&SetBop(), &Defgrd(), &DerivShapeFunctionXYZ(), gp);

  // build plastic velocity gradient from condensed variables
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper && gp >= 0 && gp < numgpt_)
    BuildDeltaLp(gp);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::IntegrateMassMatrix(
    const int gp, CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>& mass)
{
  const double density = Material()->Density(gp);
  // integrate consistent mass matrix
  const double factor = DetJ() * wgt_[gp] * density;
  double ifactor, massfactor;
  for (int inod = 0; inod < nen_; ++inod)
  {
    ifactor = ShapeFunction()(inod) * factor;
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      massfactor = ShapeFunction()(jnod) * ifactor;  // intermediate factor
      mass(3 * inod + 0, 3 * jnod + 0) += massfactor;
      mass(3 * inod + 1, 3 * jnod + 1) += massfactor;
      mass(3 * inod + 2, 3 * jnod + 2) += massfactor;
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::IntegrateStiffMatrix(const int gp,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>& stiff,
    Epetra_SerialDenseMatrix& Kda)
{
  const double detJ_w = DetJ() * wgt_[gp];

  // integrate `elastic' and `initial-displacement' stiffness matrix
  // keu = keu + (B^T . C . B) * detJ * w(gp)
  static CORE::LINALG::Matrix<numstr_, numdofperelement_> cb;
  cb.Multiply(Cmat(), Bop());
  if (fbar_)
    stiff.MultiplyTN(detJ_w * FbarFac(), Bop(), cb, 1.0);
  else
    stiff.MultiplyTN(detJ_w, Bop(), cb, 1.0);

  // integrate `geometric' stiffness matrix and add to keu *****************
  CORE::LINALG::Matrix<numstr_, 1> sfac(PK2());  // auxiliary integrated stress
  if (fbar_)
    sfac.Scale(detJ_w / FbarFac());  // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
  else
    sfac.Scale(detJ_w);  // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
  double SmB_L[nsd_];    // intermediate Sm.B_L
  // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
  for (int inod = 0; inod < nen_; ++inod)
  {
    SmB_L[0] = sfac(0) * DerivShapeFunctionXYZ()(0, inod) +
               sfac(3) * DerivShapeFunctionXYZ()(1, inod) +
               sfac(5) * DerivShapeFunctionXYZ()(2, inod);
    SmB_L[1] = sfac(3) * DerivShapeFunctionXYZ()(0, inod) +
               sfac(1) * DerivShapeFunctionXYZ()(1, inod) +
               sfac(4) * DerivShapeFunctionXYZ()(2, inod);
    SmB_L[2] = sfac(5) * DerivShapeFunctionXYZ()(0, inod) +
               sfac(4) * DerivShapeFunctionXYZ()(1, inod) +
               sfac(2) * DerivShapeFunctionXYZ()(2, inod);
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      double bopstrbop = 0.0;  // intermediate value
      for (int idim = 0; idim < 3; ++idim)
        bopstrbop += DerivShapeFunctionXYZ()(idim, jnod) * SmB_L[idim];
      stiff(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
      stiff(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
      stiff(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
    }
  }  // end of integrate `geometric' stiffness******************************

  // integrate additional fbar matrix**************************************
  if (fbar_)
  {
    static CORE::LINALG::Matrix<numstr_, 1> ccg;
    ccg.Multiply(Cmat(), RCGvec());

    static CORE::LINALG::Matrix<numdofperelement_, 1> bopccg(false);  // auxiliary integrated stress
    bopccg.MultiplyTN(detJ_w * FbarFac() / 3.0, Bop(), ccg);

    static CORE::LINALG::Matrix<numdofperelement_, 1> bops(false);  // auxiliary integrated stress
    bops.MultiplyTN(-detJ_w / FbarFac() / 3.0, Bop(), PK2());
    stiff.MultiplyNT(1., bops, Htensor(), 1.);
    stiff.MultiplyNT(1., bopccg, Htensor(), 1.);
  }
  // end of integrate additional fbar matrix*****************************

  // EAS technology: integrate matrices --------------------------------- EAS
  if (not(StrParamsInterface().GetPredictorType() == INPAR::STR::pred_tangdis))
    if (eastype_ != soh8p_easnone)
    {
      // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
      // integrate Kda: Kad += (M^T . cmat . B) * detJ * w(gp)
      // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
      static CORE::LINALG::SerialDenseMatrix cM(numstr_, neas_);  // temporary c . M
      switch (eastype_)
      {
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              cM.A(), Cmat().A(), M_eas().A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              1.0, *KaaInv_, detJ_w, M_eas(), cM);
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, numdofperelement_>(
              1.0, Kad_->A(), detJ_w, M_eas().A(), cb.A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas>(
              1.0, Kda.A(), detJ_w, cb.A(), M_eas().A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              1.0, feas_->A(), detJ_w, M_eas().A(), PK2().A());
          break;
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numstr_, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              cM.A(), Cmat().A(), M_eas().A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              1.0, *KaaInv_, detJ_w, M_eas(), cM);
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, numdofperelement_>(
              1.0, Kad_->A(), detJ_w, M_eas().A(), cb.A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, numdofperelement_, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas>(
              1.0, Kda.A(), detJ_w, cb.A(), M_eas().A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              1.0, feas_->A(), detJ_w, M_eas().A(), PK2().A());
          break;
        case soh8p_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }  // ---------------------------------------------------------------- EAS
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::IntegrateForce(
    const int gp, CORE::LINALG::Matrix<numdofperelement_, 1>& force)
{
  if (fbar_)
    force.MultiplyTN(DetJ() * wgt_[gp] / FbarFac(), Bop(), PK2(), 1.0);
  else
    force.MultiplyTN(DetJ() * wgt_[gp], Bop(), PK2(), 1.0);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::IntegrateThermoGp(
    const int gp, Epetra_SerialDenseVector& dHda)
{
  const double timefac_d =
      StrParamsInterface().GetTimIntFactorVel() / StrParamsInterface().GetTimIntFactorDisp();
  const double detJ_w = DetJ() * wgt_[gp];

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = nullptr;
  if (Material()->MaterialType() == INPAR::MAT::m_plelasthyper)
    plmat = static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("tsi only with m_plelasthyper material type");

  // Gauss point temperature
  const double gp_temp = Temp().Dot(ShapeFunction());

  // volumetric part of K_dT = Gough-Joule effect**********************************
  // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccc
  // get the thermal material tangent
  CORE::LINALG::Matrix<numstr_, 1> cTvol(true);
  CORE::LINALG::Matrix<numstr_, numstr_> dcTvoldE;
  plmat->EvaluateCTvol(&DefgrdMod(), &cTvol, &dcTvoldE, gp, Id());
  // end of call material law ccccccccccccccccccccccccccccccccccccccccccccc
  if (fbar_)
    (*dFintdT_)[gp].MultiplyTN(detJ_w / (FbarFac()), Bop(), cTvol, 0.);
  else
    (*dFintdT_)[gp].MultiplyTN(detJ_w, Bop(), cTvol, 0.);
  if (eastype_ != soh8p_easnone)
  {
    CORE::LINALG::Matrix<numstr_, nen_> cTm;
    cTm.MultiplyNT(cTvol, ShapeFunction());
    switch (eastype_)
    {
      case soh8p_easfull:
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, nen_>(
            1., KaT_->A(), detJ_w, M_eas().A(), cTm.A());
        break;
      case soh8p_easmild:
        CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
            PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, nen_>(
            1., KaT_->A(), detJ_w, M_eas().A(), cTm.A());
        break;
      case soh8p_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }

  // elastic heating ******************************************************
  plmat->HepDiss(gp) = 0.;
  plmat->dHepDT(gp) = 0.;
  plmat->dHepDissDd(gp).Size(numdofperelement_);
  plmat->dHepDissDd(gp).Scale(0.);
  if (eastype_ == soh8p_easnone)
  {
    if (fbar_)
    {
      CORE::LINALG::Matrix<3, 3> defgrd_rate_0;
      defgrd_rate_0.MultiplyTT(XcurrRate(), DerivShapeFunctionXYZ_0());

      double he_fac;
      double he_fac_deriv;
      plmat->EvaluateGoughJoule(DetF_0(), gp, Id(), he_fac, he_fac_deriv);

      double fiddfdot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) fiddfdot += InvDefgrd_0()(j, i) * defgrd_rate_0(i, j);

      double j_dot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) j_dot += DetF_0() * InvDefgrd_0()(j, i) * defgrd_rate_0(i, j);

      double He = he_fac * gp_temp * j_dot;

      plmat->HepDiss(gp) = He;

      // derivative of elastic heating w.r.t. temperature *******************
      plmat->dHepDT(gp) = he_fac * j_dot;

      CORE::LINALG::Matrix<numdofperelement_, 1> deriv_jdot_d(true);
      CORE::LINALG::Matrix<numdofperelement_, 1> deriv_j_d(true);

      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_j_d(n) += DetF_0() * InvDefgrd_0()(i, n % 3) * DerivShapeFunctionXYZ_0()(i, n / 3);

      CORE::LINALG::Matrix<3, 3> tmp;
      tmp.Multiply(InvDefgrd_0(), defgrd_rate_0);
      CORE::LINALG::Matrix<3, 3> tmp2;
      tmp2.Multiply(tmp, InvDefgrd_0());

      deriv_jdot_d.Update(fiddfdot, deriv_j_d, 1.);
      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_jdot_d(n) +=
              DetF_0() * InvDefgrd_0()(i, n % 3) * DerivShapeFunctionXYZ_0()(i, n / 3) * timefac_d -
              DetF_0() * tmp2(i, n % 3) * DerivShapeFunctionXYZ_0()(i, n / 3);

      CORE::LINALG::Matrix<numdofperelement_, 1> dHedd(true);
      dHedd.Update(gp_temp * he_fac, deriv_jdot_d, 1.);
      dHedd.Update(he_fac_deriv * gp_temp * j_dot, deriv_j_d, 1.);

      CORE::LINALG::DENSEFUNCTIONS::update<double, numdofperelement_, 1>(
          plmat->dHepDissDd(gp).A(), dHedd.A());
    }
    else
    {
      CORE::LINALG::Matrix<3, 3> defgrd_rate;
      defgrd_rate.MultiplyTT(XcurrRate(), DerivShapeFunctionXYZ());

      double he_fac;
      double he_fac_deriv;
      plmat->EvaluateGoughJoule(DetF(), gp, Id(), he_fac, he_fac_deriv);

      double fiddfdot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) fiddfdot += InvDefgrd()(j, i) * defgrd_rate(i, j);

      double j_dot = 0.;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) j_dot += DetF() * InvDefgrd()(j, i) * defgrd_rate(i, j);

      double He = he_fac * gp_temp * j_dot;

      plmat->HepDiss(gp) = He;

      // derivative of elastic heating w.r.t. temperature *******************
      plmat->dHepDT(gp) = he_fac * j_dot;

      CORE::LINALG::Matrix<numdofperelement_, 1> deriv_jdot_d(true);
      CORE::LINALG::Matrix<numdofperelement_, 1> deriv_j_d(true);

      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_j_d(n) += DetF() * InvDefgrd()(i, n % 3) * DerivShapeFunctionXYZ()(i, n / 3);

      CORE::LINALG::Matrix<3, 3> tmp;
      tmp.Multiply(InvDefgrd(), defgrd_rate);
      CORE::LINALG::Matrix<3, 3> tmp2;
      tmp2.Multiply(tmp, InvDefgrd());

      deriv_jdot_d.Update(fiddfdot, deriv_j_d, 1.);
      for (int i = 0; i < 3; ++i)
        for (int n = 0; n < numdofperelement_; ++n)
          deriv_jdot_d(n) +=
              DetF() * InvDefgrd()(i, n % 3) * DerivShapeFunctionXYZ()(i, n / 3) * timefac_d -
              DetF() * tmp2(i, n % 3) * DerivShapeFunctionXYZ()(i, n / 3);

      CORE::LINALG::Matrix<numdofperelement_, 1> dHedd(true);
      dHedd.Update(gp_temp * he_fac, deriv_jdot_d, 1.);
      dHedd.Update(he_fac_deriv * gp_temp * j_dot, deriv_j_d, 1.);

      CORE::LINALG::DENSEFUNCTIONS::update<double, numdofperelement_, 1>(
          plmat->dHepDissDd(gp).A(), dHedd.A());
    }
  }
  else
  {
    CORE::LINALG::Matrix<3, 3> defgrd_rate;
    defgrd_rate.MultiplyTT(XcurrRate(), DerivShapeFunctionXYZ());

    // Gough-Joule effect
    plmat->HepDiss(gp) = 0.;
    plmat->dHepDT(gp) = 0.;
    plmat->dHepDissDd(gp).Size(numdofperelement_);
    // Like this it should be easier to do EAS as well
    CORE::LINALG::Matrix<3, 3> RCGrate;
    RCGrate.MultiplyTN(defgrd_rate, Defgrd());
    RCGrate.MultiplyTN(1., Defgrd(), defgrd_rate, 1.);
    CORE::LINALG::Matrix<6, 1> RCGrateVec;
    for (int i = 0; i < 3; ++i) RCGrateVec(i, 0) = RCGrate(i, i);
    RCGrateVec(3, 0) = 2. * RCGrate(0, 1);
    RCGrateVec(4, 0) = 2. * RCGrate(1, 2);
    RCGrateVec(5, 0) = 2. * RCGrate(0, 2);

    // enhance the deformation rate
    if (eastype_ != soh8p_easnone)
    {
      Epetra_SerialDenseVector alpha_dot(neas_);
      switch (eastype_)
      {
        case soh8p_easmild:
          // calculate EAS-rate
          CORE::LINALG::DENSEFUNCTIONS::update<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
              0., alpha_dot, 1., *alpha_eas_);
          CORE::LINALG::DENSEFUNCTIONS::update<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., alpha_dot, -1., *alpha_eas_last_timestep_);
          alpha_dot.Scale(timefac_d);
          // enhance the strain rate
          // factor 2 because we deal with RCGrate and not GLrate
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, 1>(
              1., RCGrateVec.A(), 2., M_eas().A(), alpha_dot.A());
          break;
        case soh8p_easfull:
          // calculate EAS-rate
          CORE::LINALG::DENSEFUNCTIONS::update<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
              0., alpha_dot, 1., *alpha_eas_);
          CORE::LINALG::DENSEFUNCTIONS::update<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., alpha_dot, -1., *alpha_eas_last_timestep_);
          alpha_dot.Scale(timefac_d);
          // enhance the strain rate
          // factor 2 because we deal with RCGrate and not GLrate
          CORE::LINALG::DENSEFUNCTIONS::multiply<double, numstr_,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, 1>(
              1., RCGrateVec.A(), 2., M_eas().A(), alpha_dot.A());
          break;
        case soh8p_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }  // enhance the deformation rate

    // heating ************************************************************
    double He = .5 * gp_temp * cTvol.Dot(RCGrateVec);

    plmat->HepDiss(gp) = He;

    // derivative of elastic heating w.r.t. temperature *******************
    plmat->dHepDT(gp) = .5 * cTvol.Dot(RCGrateVec);

    // derivative of elastic heating w.r.t. displacement ******************
    CORE::LINALG::Matrix<numdofperelement_, 1> dHedd(true);
    CORE::LINALG::Matrix<6, nen_ * nsd_> boprate(false);  // (6x24)
    CalculateBop(&boprate, &defgrd_rate, &DerivShapeFunctionXYZ(), gp);

    CORE::LINALG::Matrix<6, 1> tmp61;
    tmp61.MultiplyTN(dcTvoldE, RCGrateVec);
    dHedd.MultiplyTN(.5 * gp_temp, Bop(), tmp61, 1.);

    dHedd.MultiplyTN(timefac_d * gp_temp, Bop(), cTvol, 1.);
    dHedd.MultiplyTN(gp_temp, boprate, cTvol, 1.);

    // derivative of elastic heating w.r.t. EAS alphas *******************
    if (eastype_ != soh8p_easnone)
    {
      switch (eastype_)
      {
        case soh8p_easmild:
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              0., dHda.A(), .5 * gp_temp, M_eas().A(), tmp61.A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>::neas, numstr_, 1>(
              1., dHda.A(), gp_temp * timefac_d, M_eas().A(), cTvol.A());
          break;
        case soh8p_easfull:
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              0., dHda.A(), .5 * gp_temp, M_eas().A(), tmp61.A());
          CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double,
              PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>::neas, numstr_, 1>(
              1., dHda.A(), gp_temp * timefac_d, M_eas().A(), cTvol.A());
          break;
        case soh8p_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }
    }

    plmat->dHepDissDd(gp).Scale(0.);
    CORE::LINALG::DENSEFUNCTIONS::update<double, numdofperelement_, 1>(
        plmat->dHepDissDd(gp).A(), dHedd.A());
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::HeatFlux(const std::vector<double>& temp,
    const std::vector<double>& disp, const CORE::LINALG::Matrix<nsd_, 1>& xi,
    const CORE::LINALG::Matrix<nsd_, 1>& n, double& q, Epetra_SerialDenseMatrix* dq_dT,
    Epetra_SerialDenseMatrix* dq_dd, CORE::LINALG::Matrix<nsd_, 1>* dq_dn,
    CORE::LINALG::Matrix<nsd_, 1>* dq_dpxi, Epetra_SerialDenseMatrix* d2q_dT_dd,
    Epetra_SerialDenseMatrix* d2q_dT_dn, Epetra_SerialDenseMatrix* d2q_dT_dpxi)
{
  if (!dq_dT || !dq_dd || !dq_dn || !dq_dpxi || !d2q_dT_dd || !d2q_dT_dn || !d2q_dT_dpxi)
    dserror("input inconsistent");

  InvalidGpData();
  InvalidEleData();

  if (NumMaterial() < 2) dserror("where's my second material");
  Teuchos::RCP<MAT::FourierIso> mat_thr =
      Teuchos::rcp_dynamic_cast<MAT::FourierIso>(Material(1), true);
  const double k0 = mat_thr->Conductivity();

  std::vector<double> vel(0);
  FillPositionArrays(disp, vel, temp);

  // shape functions (shapefunct) and their first derivatives (deriv)
  CORE::DRT::UTILS::shape_function<distype>(xi, SetShapeFunction());
  CORE::DRT::UTILS::shape_function_deriv1<distype>(xi, SetDerivShapeFunction());

  Kinematics();

  CORE::LINALG::Matrix<3, 1> GradT;
  GradT.Multiply(DerivShapeFunctionXYZ(), Temp());

  CORE::LINALG::Matrix<3, 1> iFn;
  iFn.Multiply(InvDefgrd(), n);

  CORE::LINALG::Matrix<9, numdofperelement_> dFdd;
  for (int k = 0; k < nen_; ++k)
  {
    for (int d = 0; d < 3; ++d) dFdd(d, k * nsd_ + d) = DerivShapeFunctionXYZ()(d, k);
    dFdd(3, k * nsd_ + 0) = DerivShapeFunctionXYZ()(1, k);
    dFdd(4, k * nsd_ + 1) = DerivShapeFunctionXYZ()(2, k);
    dFdd(5, k * nsd_ + 0) = DerivShapeFunctionXYZ()(2, k);

    dFdd(6, k * nsd_ + 1) = DerivShapeFunctionXYZ()(0, k);
    dFdd(7, k * nsd_ + 2) = DerivShapeFunctionXYZ()(1, k);
    dFdd(8, k * nsd_ + 2) = DerivShapeFunctionXYZ()(0, k);
  }

  q = -k0 / DetF() * GradT.Dot(iFn);

  CORE::LINALG::Matrix<nen_, 9> dq_dT_dF;

  if (dq_dT)
  {
    dq_dT->Shape(nen_, 1);
    for (int n = 0; n < nen_; ++n)
      for (int d = 0; d < nsd_; ++d)
        (*dq_dT)(n, 0) += -k0 / DetF() * DerivShapeFunctionXYZ()(d, n) * iFn(d);
  }

  if (d2q_dT_dd || dq_dd)
  {
    d2q_dT_dd->Shape(nen_, nen_ * nsd_);
    CORE::LINALG::Matrix<nen_, nen_ * nsd_> d2q_dT_dd_m(d2q_dT_dd->A(), true);

    CORE::LINALG::Matrix<3, nen_> tmp;
    tmp.MultiplyTN(InvDefgrd(), DerivShapeFunctionXYZ());
    CORE::LINALG::Matrix<nen_, 1> tmp2;
    tmp2.MultiplyTN(DerivShapeFunctionXYZ(), iFn);

    CORE::LINALG::Matrix<9, nen_> dq_dF_v;
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          dq_dF_v(VoigtMapping::NonSymToVoigt9(a, b), c) =
              k0 / DetF() * (tmp2(c) * InvDefgrd()(b, a) + tmp(a, c) * iFn(b));

    d2q_dT_dd_m.MultiplyTN(dq_dF_v, dFdd);

    dq_dd->Shape(nen_ * nsd_, 1);
    CORE::LINALG::Matrix<nen_ * nsd_, 1> dq_dd_m(dq_dd->A(), true);
    dq_dd_m.MultiplyTN(d2q_dT_dd_m, Temp());
  }

  if (d2q_dT_dn || dq_dn)
  {
    d2q_dT_dn->Shape(nen_, nsd_);
    CORE::LINALG::Matrix<nen_, nsd_> d2q_dT_dn_m(d2q_dT_dn->A(), true);
    d2q_dT_dn_m.MultiplyTN(-k0 / DetF(), DerivShapeFunctionXYZ(), InvDefgrd());

    dq_dn->MultiplyTN(d2q_dT_dn_m, Temp());
  }

  if (dq_dpxi || d2q_dT_dpxi)
  {
    d2q_dT_dpxi->Shape(nen_, nsd_);
    CORE::LINALG::Matrix<nen_, nsd_> d2q_dT_dpxi_m(d2q_dT_dpxi->A(), true);

    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2, nen_> deriv2;
    CORE::DRT::UTILS::shape_function_deriv2<distype>(xi, deriv2);
    const int d2v[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
    CORE::LINALG::Matrix<3, 1> tmp;
    tmp.MultiplyTN(InvJ(), iFn);

    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2, nen_> tmp3;
    tmp3.Update(deriv2);

    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += -k0 / DetF() * tmp3(d2v[a][b], c) * tmp(b);

    CORE::LINALG::Matrix<nen_, nen_> tmp2;
    tmp2.Multiply(Xrefe(), DerivShapeFunctionXYZ());
    tmp3.Multiply(deriv2, tmp2);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / DetF() * tmp3(d2v[a][b], c) * tmp(b);

    CORE::LINALG::Matrix<nen_, nsd_> xXF(Xcurr());
    xXF.MultiplyNT(-1., Xrefe(), Defgrd(), 1.);

    CORE::LINALG::Matrix<nsd_, CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2> xXFsec;
    xXFsec.MultiplyTT(xXF, deriv2);

    CORE::LINALG::Matrix<nsd_, nsd_> tmp4;
    tmp4.MultiplyTN(InvJ(), InvDefgrd());

    CORE::LINALG::Matrix<nsd_, CORE::DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2> tmp5;
    tmp5.Multiply(tmp4, xXFsec);

    CORE::LINALG::Matrix<nen_, 1> tmp6;
    tmp6.MultiplyTN(DerivShapeFunctionXYZ(), iFn);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / DetF() * tmp6(c) * tmp5(b, d2v[a][b]);

    CORE::LINALG::Matrix<nsd_, nen_> tmp7;
    tmp7.MultiplyTN(InvDefgrd(), DerivShapeFunctionXYZ());
    tmp3.MultiplyTN(xXFsec, tmp7);
    tmp.MultiplyTN(InvJ(), iFn);
    for (int a = 0; a < nsd_; ++a)
      for (int b = 0; b < nsd_; ++b)
        for (int c = 0; c < nen_; ++c)
          (*d2q_dT_dpxi)(c, a) += +k0 / DetF() * tmp3(d2v[a][b], c) * tmp(b);

    dq_dpxi->MultiplyTN(d2q_dT_dpxi_m, Temp());
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetNurbsEleInfo(DRT::Discretization* dis)
{
  if (!IsNurbsElement()) return;

  if (dis == nullptr) dis = DRT::Problem::Instance()->GetDis("structure").get();

  dynamic_cast<DRT::NURBS::NurbsDiscretization*>(dis)->GetKnotVector()->GetEleKnots(
      SetKnots(), Id());
  for (int i = 0; i < nen_; ++i)
    SetWeights()(i) = dynamic_cast<DRT::NURBS::ControlPoint*>(Nodes()[i])->W();
}

// template functions
template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::CondensePlasticity<5>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::CondensePlasticity<8>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex18>::CondensePlasticity<5>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex18>::CondensePlasticity<8>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>::CondensePlasticity<5>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>::CondensePlasticity<8>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::nurbs27>::CondensePlasticity<5>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::nurbs27>::CondensePlasticity<8>(
    const CORE::LINALG::Matrix<nsd_, nsd_>&, const CORE::LINALG::Matrix<nsd_, nsd_>&,
    const CORE::LINALG::Matrix<numstr_, numdofperelement_>&,
    const CORE::LINALG::Matrix<nsd_, nen_>*, const CORE::LINALG::Matrix<numstr_, 1>*, const double,
    const int, const double, Teuchos::ParameterList&, CORE::LINALG::Matrix<numdofperelement_, 1>*,
    CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*, const Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, std::vector<Epetra_SerialDenseVector>*, const double*,
    const CORE::LINALG::Matrix<numdofperelement_, 1>*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>::HeatFlux(const std::vector<double>&,
    const std::vector<double>&, const CORE::LINALG::Matrix<nsd_, 1>&,
    const CORE::LINALG::Matrix<nsd_, 1>&, double&, Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, CORE::LINALG::Matrix<nsd_, 1>*, CORE::LINALG::Matrix<nsd_, 1>*,
    Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*);


template void DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>::HeatFlux(const std::vector<double>&,
    const std::vector<double>&, const CORE::LINALG::Matrix<nsd_, 1>&,
    const CORE::LINALG::Matrix<nsd_, 1>&, double&, Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, CORE::LINALG::Matrix<nsd_, 1>*, CORE::LINALG::Matrix<nsd_, 1>*,
    Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::tet4>::HeatFlux(const std::vector<double>&,
    const std::vector<double>&, const CORE::LINALG::Matrix<nsd_, 1>&,
    const CORE::LINALG::Matrix<nsd_, 1>&, double&, Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, CORE::LINALG::Matrix<nsd_, 1>*, CORE::LINALG::Matrix<nsd_, 1>*,
    Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*);

template void DRT::ELEMENTS::So3_Plast<DRT::Element::nurbs27>::HeatFlux(const std::vector<double>&,
    const std::vector<double>&, const CORE::LINALG::Matrix<nsd_, 1>&,
    const CORE::LINALG::Matrix<nsd_, 1>&, double&, Epetra_SerialDenseMatrix*,
    Epetra_SerialDenseMatrix*, CORE::LINALG::Matrix<nsd_, 1>*, CORE::LINALG::Matrix<nsd_, 1>*,
    Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*, Epetra_SerialDenseMatrix*);
