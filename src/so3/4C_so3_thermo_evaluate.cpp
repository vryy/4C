// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_mat_thermoplasticlinelast.hpp"
#include "4C_mat_trait_thermo_solid.hpp"
#include "4C_so3_thermo.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | pre-evaluate the element (public)                         dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::pre_evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la)
{
  // if the coupling variables are required before evaluate() is called the 1st
  // time
  // here for Robinson's material
  if (la.size() > 1)
  {
    // the temperature field has only one dof per node, disregarded by the
    // dimension of the problem
    const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);

    if (discretization.has_state(1, "temperature"))
    {
      if (la[1].size() != nen_ * numdofpernode_thr)
      {
        FOUR_C_THROW(
            "Location vector length for temperatures does not match!\n"
            "la[1].Size()= %i\tnen_*numdofpernode_thr= %i",
            la[1].size(), nen_ * numdofpernode_thr);
      }
      // check if you can get the temperature state
      std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
          discretization.get_state(1, "temperature");
      if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'.");

      // extract local values of the global vectors
      std::shared_ptr<std::vector<double>> nodaltempnp =
          std::make_shared<std::vector<double>>(la[1].lm_.size());
      Core::FE::extract_my_values(*tempnp, *nodaltempnp, la[1].lm_);

      // now set the current temperature vector in the parameter list
      params.set<std::shared_ptr<std::vector<double>>>("nodal_tempnp", nodaltempnp);
    }
  }  // initial temperature dependence
}  // pre_evaluate()


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::Elements::So3Thermo<So3Ele, distype>::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // set the pointer to the parameter list in element
  So3Ele::set_params_interface_ptr(params);

  // what actions are available
  // (action == "calc_struct_stifftemp")
  // (action == "calc_struct_stress")
  // (action == default)

  // start with "none"
  typename So3Thermo::ActionType act = So3Thermo::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_struct_stifftemp")
    act = So3Thermo::calc_struct_stifftemp;
  else if (action == "calc_struct_stress")
    act = So3Thermo::calc_struct_stress;

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // coupling terms K_dT in stiffness matrix K^{TSI} for monolithic TSI
    case So3Thermo::calc_struct_stifftemp:
    {
      if (la.size() > 1)
      {
        evaluate_coupl_with_thr(params, discretization, la, elemat1_epetra, elemat2_epetra,
            elevec1_epetra, elevec2_epetra, elevec3_epetra);
      }
      break;
    }

    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      // here you can pass, e.g. for Robinson's material the current temperature
      // T_n+1 needed to calculate e.g. the young's modulus E(T_n+1)
      pre_evaluate(params, discretization, la);

      // call the purely structural methods
      So3Ele::evaluate(params, discretization,
          la[0].lm_,  // only the first column, i.e. the structural field is passed
          elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      // add the temperature-dependent terms to the structural field, i.e.
      // it's a TSI problem
      if (la.size() > 1)
      {
        evaluate_coupl_with_thr(params, discretization,
            la,  // coupled TSI is considered, i.e. pass the completed location array
            elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
      }
      break;
    }  // default

  }  // action

  return 0;
}  // evaluate()


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 08/12 |
 | here is the action for the coupling to the thermal field             |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::Elements::So3Thermo<So3Ele, distype>::evaluate_coupl_with_thr(
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
    Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  ActionType act = none;

  // get the required action for coupling with the thermal field
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_struct_internalforce")
    act = calc_struct_internalforce;
  else if (action == "calc_struct_nlnstiff")
    act = calc_struct_nlnstiff;
  else if (action == "calc_struct_nlnstiffmass")
    act = calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stifftemp")
    act = calc_struct_stifftemp;
  else if (action == "calc_struct_stress")
    act = calc_struct_stress;
  else if (action == "calc_struct_reset_istep")
    act = calc_struct_reset_istep;
  else if (action == "calc_struct_update_istep")
    act = calc_struct_update_istep;
  else if (action == "calc_struct_energy")
    act = calc_struct_energy;
  else if (action == "calc_struct_predict")
    return 0;
  else if (action == "calc_struct_recover")
    return 0;
  else
    FOUR_C_THROW("Unknown type of action for So3_Thermo: %s", action.c_str());

  // what should the element do
  switch (act)
  {
    //============================================================================
    // internal force vector for TSI only
    case calc_struct_internalforce:
    {
      // internal force vector
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat1+2, elevec2+3 are not used anyway

      // need current displacement and residual/incremental displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
          discretization.get_state(0, "displacement");
      std::shared_ptr<const Core::LinAlg::Vector<double>> res =
          discretization.get_state(0, "residual displacement");

      if ((disp == nullptr) or (res == nullptr))
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");

      // build the location vector only for the structure field
      std::vector<double> mydisp((la[0].lm_).size());
      Core::FE::extract_my_values(*disp, mydisp, la[0].lm_);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> myemat(true);

      // initialise the vectors
      // evaluate() is called the first time in Thermo::BaseAlgorithm: at this stage the
      // coupling field is not yet known. Pass coupling vectors filled with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);

      // need current temperature state, call the temperature discretization
      // disassemble temperature
      if (discretization.has_state(1, "temperature"))
      {
        // check if you can get the temperature state
        std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
            discretization.get_state(1, "temperature");
        if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
        if (la[1].size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");
        // extract the current temperatures
        Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);

        // default: geometrically non-linear analysis with Total Lagrangean approach
        if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::nonlinearTotLag)
        {
          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat1(
              elemat1_epetra.values(), true);


          nln_stifffint_tsi(la,          // location array
              discretization,            // discr
              mydisp,                    // current displacements
              mytempnp,                  // current temperature
              nullptr,                   // element stiffness matrix
              &elevec1,                  // element internal force vector
              nullptr,                   // stresses at GP
              params,                    // algorithmic parameters e.g. time
              Inpar::Solid::stress_none  // stress output option
          );

        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::nonlinearTotLag)

        // geometric Inpar::Solid::KinemType::linear
        else if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::linear)
        {
          // calculate the THERMOmechanical term for fint
          lin_fint_tsi(la, mydisp, mytempnp, &elevec1, nullptr, params, Inpar::Solid::stress_none);
        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::linear)
      }
      break;
    }  // calc_struct_internalforce

    //============================================================================
    // (non)linear stiffness for TSI
    case calc_struct_nlnstiff:
    {
      // internal force vector
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat2, elevec2+3 are not used anyway
      // elemat1 only for geometrically nonlinear analysis

      // need current displacement and residual/incremental displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
          discretization.get_state(0, "displacement");

      if (disp == nullptr) FOUR_C_THROW("Cannot get state vectors 'displacement' ");

      // build the location vector only for the structure field
      std::vector<double> mydisp((la[0].lm_).size());
      Core::FE::extract_my_values(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // evaluate() is called the first time in Thermo::BaseAlgorithm: at this stage the
      // coupling field is not yet known. Pass coupling vectors filled with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);

      // need current temperature state, call the temperature discretization
      // disassemble temperature
      if (discretization.has_state(1, "temperature"))
      {
        // check if you can get the temperature state
        std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
            discretization.get_state(1, "temperature");
        if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
        if (la[1].size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");
        // extract the current temperatures
        Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);

        // default: geometrically non-linear analysis with Total Lagrangean approach
        if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::nonlinearTotLag)
        {
          // stiffness
          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat1(
              elemat1_epetra.values(), true);

          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>* matptr = nullptr;
          if (elemat1.is_initialized()) matptr = &elemat1;

          nln_stifffint_tsi(la,          // location array
              discretization,            // discr
              mydisp,                    // current displacements
              mytempnp,                  // current temperature
              matptr,                    // element stiffness matrix
              &elevec1,                  // element internal force vector
              nullptr,                   // stresses at GP
              params,                    // algorithmic parameters e.g. time
              Inpar::Solid::stress_none  // stress output option
          );

        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::nonlinearTotLag)

        // geometric linear
        else if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::linear)
        {
          // calculate the THERMOmechanical term for fint
          lin_fint_tsi(la, mydisp, mytempnp, &elevec1, nullptr, params, Inpar::Solid::stress_none);
        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::linear)
      }

      break;
    }  // calc_struct_nlnstiff

    //============================================================================
    // (non)linear stiffness, mass matrix and internal force vector for TSI
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // internal force vector
      Core::LinAlg::Matrix<numdofperelement_, 1> elevec1(elevec1_epetra.values(), true);
      // elevec2+3 and elemat2 are not used anyway,
      // elemat1 only for geometrically nonlinear analysis

      // need current displacement and residual/incremental displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
          discretization.get_state(0, "displacement");

      if (disp == nullptr) FOUR_C_THROW("Cannot get state vectors 'displacement'");

      // build the location vector only for the structure field
      std::vector<double> mydisp((la[0].lm_).size());
      Core::FE::extract_my_values(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);

      // need current temperature state, call the temperature discretization
      // disassemble temperature
      if (discretization.has_state(1, "temperature"))
      {
        // check if you can get the temperature state
        std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
            discretization.get_state(1, "temperature");
        if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
        if (la[1].size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");
        // extract the current temperatures
        Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);

        // default: geometrically non-linear analysis with Total Lagrangean approach
        if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::nonlinearTotLag)
        {
          // stiffness
          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_> elemat1(
              elemat1_epetra.values(), true);

          nln_stifffint_tsi(la,          // location array
              discretization,            // discr
              mydisp,                    // current displacements
              mytempnp,                  // current temperature
              &elemat1,                  // element stiffness matrix
              &elevec1,                  // element internal force vector
              nullptr,                   // stresses at GP
              params,                    // algorithmic parameters e.g. time
              Inpar::Solid::stress_none  // stress output option
          );

        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::nonlinearTotLag)

        // geometric linear
        else if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::linear)
        {
          // build the current temperature vector
          Core::LinAlg::Matrix<nen_ * numdofpernode_, 1> etemp(&(mytempnp[1]), true);  // view only!
          // calculate the THERMOmechanical term for fint
          lin_fint_tsi(la, mydisp, mytempnp, &elevec1, nullptr, params, Inpar::Solid::stress_none);
        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::linear)
      }

      break;
    }  // calc_struct_nlnstiff(l)mass

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // elemat1+2,elevec1-3 are not used anyway
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
          discretization.get_state(0, "displacement");
      if (disp == nullptr) FOUR_C_THROW("Cannot get state vectors 'displacement'");

      std::vector<double> mydisp((la[0].lm_).size());
      Core::FE::extract_my_values(*disp, mydisp, la[0].lm_);

      std::shared_ptr<std::vector<char>> couplstressdata;
      Inpar::Solid::StressType iocouplstress;
      if (this->is_params_interface())
      {
        couplstressdata = this->str_params_interface().coupling_stress_data_ptr();
        iocouplstress = this->str_params_interface().get_coupling_stress_output_type();
      }
      else
      {
        couplstressdata = params.get<std::shared_ptr<std::vector<char>>>("couplstress", nullptr);
        iocouplstress =
            params.get<Inpar::Solid::StressType>("iocouplstress", Inpar::Solid::stress_none);
      }

      // get the temperature dependent stress
      Core::LinAlg::Matrix<numgpt_post, numstr_> couplstress(true);

      // initialise the vectors
      // evaluate() is called the first time in Thermo::BaseAlgorithm: at this stage the
      // coupling field is not yet known. Pass coupling vectors filled with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);

      // need current temperature state,
      // call the temperature discretization: thermo equates 2nd dofset
      // disassemble temperature
      if (discretization.has_state(1, "temperature"))
      {
        // check if you can get the temperature state
        std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
            discretization.get_state(1, "temperature");
        if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
        if (la[1].size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");

        // extract the current temperatures
        Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);

        // default: geometrically non-linear analysis with Total Lagrangean approach
        if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::nonlinearTotLag)
        {
#ifdef TSIASOUTPUT
          std::cout << "thermal stress" << couplstress << std::endl;
          std::cout << "iocouplstress = " << iocouplstress << std::endl;
#endif

          // calculate the thermal stress
          nln_stifffint_tsi(la,  // location array
              discretization,    // discr
              mydisp,            // current displacements
              mytempnp,          // current temperature
              nullptr,           // element stiffness matrix
              nullptr,           // element internal force vector
              &couplstress,      // stresses at GP
              params,            // algorithmic parameters e.g. time
              iocouplstress      // stress output option
          );


#ifdef TSIASOUTPUT
          std::cout << "thermal stress" << couplstress << std::endl;
#endif

        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::nonlinearTotLag)

        // geometric linear
        else if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::linear)
        {
          // purely structural method, this is the coupled routine, i.e., a 2nd
          // discretisation exists, i.e., --> we always have a temperature state

          // calculate the THERMOmechanical term for fint: temperature stresses
          lin_fint_tsi(la, mydisp, mytempnp, nullptr, &couplstress, params, iocouplstress);
        }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::linear)

#ifdef TSIASOUTPUT
        std::cout << "thermal stress" << couplstress << std::endl;
#endif
      }

      // total stress is the sum of the mechanical stress and the thermal stress
      // stress = stress_d + stress_T
      //        stress.update(1.0,couplstress,1.0);
      // --> so far the addition of s_d and s_T was realised here
      // ==> from now on: we fill 2 different vectors (stressdata,couplstressdata)
      //     which are used in the post processing.
      //     --> advantage: different numbers of Gauss points for the stress
      //         and couplstress are possible
      //         --> important e.g. in case of Tet4 (s_d: 1GP, s_T: 5GP)
      //             --> for s_T we use the library intrepid
      //     --> in ParaView you can visualise the mechanical and the thermal
      //         stresses separately
      //     --> to get the total stress you have to calculate both vectors
      //         within ParaView using programmable filters

      // pack the data for postprocessing
      {
        Core::Communication::PackBuffer data;
        add_to_pack(data, couplstress);
        std::copy(data().begin(), data().end(), std::back_inserter(*couplstressdata));
      }

      break;
    }  // calc_struct_stress

    //============================================================================
    // required for predictor TangDis --> can be helpful in compressible case!
    case calc_struct_reset_istep:
    {
      // do nothing, actual implementation is in So3Ele
      break;
    }

    //============================================================================
    case calc_struct_update_istep:
    {
      auto thermoSolid = std::dynamic_pointer_cast<Mat::Trait::ThermoSolid>(material());
      if (thermoSolid != nullptr)
      {
        std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);
        if (discretization.has_state(1, "temperature"))
        {
          // check if you can get the temperature state
          std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
              discretization.get_state(1, "temperature");
          if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the
          // dimension of the problem
          const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
          if (la[1].size() != nen_ * numdofpernode_thr)
            FOUR_C_THROW("Location vector length for temperature does not match!");
          // extract the current temperatures
          Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);
        }

        // vector of the current element temperatures
        Core::LinAlg::Matrix<nen_, 1> etemp(mytempnp.data());
        Core::LinAlg::Matrix<nen_, 1> shapefunct;

        // --------------------------------------------------
        for (int gp = 0; gp < numgpt_; ++gp)
        {
          Core::FE::shape_function<distype>(xsi_[gp], shapefunct);

          // product of shapefunctions and element temperatures
          Core::LinAlg::Matrix<1, 1> NT(false);
          NT.multiply_tn(shapefunct, etemp);

          thermoSolid->reinit(nullptr, nullptr, NT(0), map_my_gp_to_so_hex8(gp));
          thermoSolid->commit_current_state();
        }
      }

      break;
    }  // calc_struct_update_istep

    //============================================================================
    // coupling term k_dT of stiffness matrix for monolithic TSI
    case calc_struct_stifftemp:
    {
      // mechanical-thermal system matrix
      Core::LinAlg::Matrix<numdofperelement_, nen_> stiffmatrix_kdT(elemat1_epetra.values(), true);
      // elemat2,elevec1-3 are not used anyway
      // need current displacement and residual/incremental displacements
      std::shared_ptr<const Core::LinAlg::Vector<double>> disp =
          discretization.get_state(0, "displacement");
      if (disp == nullptr) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp((la[0].lm_).size());
      // build the location vector only for the structure field
      Core::FE::extract_my_values(*disp, mydisp, la[0].lm_);

      // initialise the vectors
      // evaluate() is called the first time in structure_base_algorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      std::vector<double> mytempnp(((la[0].lm_).size()) / nsd_, 0.0);

      // need current temperature state, call the temperature discretization
      // disassemble temperature
      if (discretization.has_state(1, "temperature"))
      {
        // check if you can get the temperature state
        std::shared_ptr<const Core::LinAlg::Vector<double>> tempnp =
            discretization.get_state(1, "temperature");
        if (tempnp == nullptr) FOUR_C_THROW("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_thr = discretization.num_dof(1, nodes()[0]);
        if (la[1].size() != nen_ * numdofpernode_thr)
          FOUR_C_THROW("Location vector length for temperature does not match!");
        // extract the current temperatures
        Core::FE::extract_my_values(*tempnp, mytempnp, la[1].lm_);
      }
      // default: geometrically non-linear analysis with Total Lagrangean approach
      if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::nonlinearTotLag)
      {
        // calculate the mechanical-thermal sub matrix k_dT of K_TSI
        nln_kd_t_tsi(la, discretization, mydisp, mytempnp, &stiffmatrix_kdT, params);

      }  // (So3Ele::KinematicType() == nonlinear)

      // geometric linear
      else if (So3Ele::kinematic_type() == Inpar::Solid::KinemType::linear)
      {
        // calculate the mechanical-thermal sub matrix k_dT of K_TSI
        lin_kd_t_tsi(la, mydisp, mytempnp, &stiffmatrix_kdT, params);
      }  // (So3Ele::KinematicType() == Inpar::Solid::KinemType::linear)

      break;
    }  // calc_struct_stifftemp

    // strain energy: the solid element knows what to do...
    case calc_struct_energy:
    {
      So3Ele::evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
    //============================================================================
    default:
      FOUR_C_THROW("Unknown type of action for So3_Thermo");
      break;
  }  // action

  return 0;
}  // evaluate_coupl_with_thr()


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 05/10 |
 | contribution to r_d (private)                                        |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::lin_fint_tsi(
    Core::Elements::LocationArray& la,                      // location array
    std::vector<double>& disp,                              // current displacements
    std::vector<double>& temp,                              // current temperature
    Core::LinAlg::Matrix<numdofperelement_, 1>* force,      // element internal force vector
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
    Teuchos::ParameterList& params,                         // algorithmic parameters e.g. time
    const Inpar::Solid::StressType iostress                 // stress output option
)
{
  // update element geometry hex8, 3D: (8x3)
  Core::LinAlg::Matrix<nen_, nsd_> xrefe;  // X, material coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurr;  // x, current  coord. of element
  // vector of the current element temperatures
  Core::LinAlg::Matrix<nen_, 1> etemp;


  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * numdofpernode_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * numdofpernode_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * numdofpernode_ + 2];

    etemp(i, 0) = temp[i + 0];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  Core::LinAlg::Matrix<nsd_, nen_> N_XYZ;
  // build deformation gradient w.r.t. to material configuration
  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(true);
  // shape functions and their first derivatives
  Core::LinAlg::Matrix<nen_, 1> shapefunct;
  Core::LinAlg::Matrix<nsd_, nen_> deriv;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.multiply(invJ_[gp], deriv);  // (6.21)
    double detJ = detJ_[gp];           // (6.22)

    // geometrically linear, i.e. reference == current state, i.e. F == I
    // set to initial state (defgrd == identity)
    for (int i = 0; i < 3; ++i) defgrd(i, i) = 1.0;

    // calculate the linear B-operator
    Core::LinAlg::Matrix<numstr_, numdofperelement_> boplin;
    calculate_boplin(&boplin, &N_XYZ);

    // call material law

    // product of shapefunctions and element temperatures for couplstress
    // N_T . T
    Core::LinAlg::Matrix<1, 1> NT(false);
    NT.multiply_tn(shapefunct, etemp);
    // scalar-valued current element temperature T_{n+1}
    // temperature-dependent material parameters, i.e. E(T), pass T_{n+1}
    double scalartemp = NT(0, 0);
    // insert T_{n+1} into parameter list
    params.set<double>("temperature", scalartemp);

    // calculate the stress part dependent on the temperature in the material
    Core::LinAlg::Matrix<numstr_, 1> ctemp(true);
    Core::LinAlg::Matrix<numstr_, 1> couplstress(true);
    Core::LinAlg::Matrix<numstr_, numstr_> cmat(true);
    Core::LinAlg::Matrix<numstr_, 1> glstrain(true);

    // insert strain increment into parameter list
    // calculate iterative strains
    Core::LinAlg::Matrix<numstr_, 1> straininc(true);
    params.set<Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>>("straininc", straininc);

    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) couplstress = C . Delta T
    // do not call the material for Robinson's material
    if (material()->material_type() != Core::Materials::m_vp_robinson)
      materialize(&couplstress, &ctemp, &NT, &cmat, &glstrain, params);

    // end of call material law

    // return gp stresses
    switch (iostress)
    {
      case Inpar::Solid::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");

        for (int i = 0; i < numstr_; ++i) (*elestress)(gp, i) = couplstress(i);
        break;
      }
      case Inpar::Solid::stress_cauchy:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");

        // push forward of material stress to the spatial configuration
        Core::LinAlg::Matrix<nsd_, nsd_> cauchycouplstress;
        p_k2to_cauchy(&couplstress, &defgrd, &cauchycouplstress);

        (*elestress)(gp, 0) = cauchycouplstress(0, 0);
        (*elestress)(gp, 1) = cauchycouplstress(1, 1);
        (*elestress)(gp, 2) = cauchycouplstress(2, 2);
        (*elestress)(gp, 3) = cauchycouplstress(0, 1);
        (*elestress)(gp, 4) = cauchycouplstress(1, 2);
        (*elestress)(gp, 5) = cauchycouplstress(0, 2);
        break;
      }
      case Inpar::Solid::stress_none:
        break;
      default:
        FOUR_C_THROW("requested stress type not available");
        break;
    }

    // integrate internal force vector r_d
    // f = f + (B^T . sigma_temp) * detJ * w(gp)
    if (force != nullptr)
    {
      // old implementation hex8_thermo double detJ_w = detJ*gpweights[gp];
      double detJ_w = detJ * intpoints_.weight(gp);  // gpweights[gp];
      force->multiply_tn(detJ_w, boplin, couplstress, 1.0);
    }  // if (force != nullptr)

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}  // lin_fint_tsi()


/*----------------------------------------------------------------------*
 | evaluate only the mechanical-thermal stiffness term       dano 03/11 |
 | for monolithic TSI, contribution to k_dT (private)                   |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::lin_kd_t_tsi(Core::Elements::LocationArray& la,
    std::vector<double>& disp,                                       // current displacement
    std::vector<double>& temp,                                       // current temperatures
    Core::LinAlg::Matrix<numdofperelement_, nen_>* stiffmatrix_kdT,  // (nsd_*nen_ x nen_)
    Teuchos::ParameterList& params)
{
  // update element geometry (8x3)
  Core::LinAlg::Matrix<nen_, nsd_> xrefe(false);  // X, material coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurr(false);  // x, current  coord. of element
  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * numdofpernode_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * numdofpernode_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * numdofpernode_ + 2];
  }

  // get current element displacements in vector notation --> use for elestrain
  Core::LinAlg::Matrix<numdofperelement_, 1> edisp(false);
  for (int i = 0; i < numdofperelement_; i++)
  {
    edisp(i, 0) = disp[i + 0];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  Core::LinAlg::Matrix<nsd_, nen_> N_XYZ;
  // shape functions and their first derivatives
  Core::LinAlg::Matrix<nen_, 1> shapefunct;
  Core::LinAlg::Matrix<nsd_, nen_> deriv;

  // ------------------------------------------------ initialise material
  // get the thermal material tangent
  Core::LinAlg::Matrix<numstr_, 1> ctemp(true);

  // for the stiffness matrix k_dT
  Core::LinAlg::Matrix<nen_, 1> etemp(false);

  // get the temperature vector
  for (int i = 0; i < nen_; ++i) etemp(i, 0) = temp[i];

  // initialize matrices required for k_dT
  Core::LinAlg::Matrix<numdofperelement_, 1> Bstress_T(true);
  Core::LinAlg::Matrix<numdofperelement_, 1> Bcouplstress_T(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.multiply(invJ_[gp], deriv);  // (6.21)
    double detJ = detJ_[gp];           // (6.22)

    // calculate the linear B-operator B_L = N_XYZ
    Core::LinAlg::Matrix<numstr_, numdofperelement_> boplin;
    calculate_boplin(&boplin, &N_XYZ);

    // call material law

    std::shared_ptr<Mat::Trait::ThermoSolid> thermoSolidMaterial =
        std::dynamic_pointer_cast<Mat::Trait::ThermoSolid>(material());
    if (thermoSolidMaterial != nullptr)
    {
      // calculate the nodal strains: strain = B . d
      Core::LinAlg::Matrix<numstr_, 1> strain(false);
      strain.multiply(boplin, edisp);

      Core::LinAlg::Matrix<1, 1> NT(false);
      NT.multiply_tn(shapefunct, etemp);  // (1x1)
      thermoSolidMaterial->reinit(nullptr, &strain, NT(0), map_my_gp_to_so_hex8(gp));
      // full thermal derivative of stress wrt to scalar temperature (needs to be post-multiplied
      // with shape functions)

      // this element will be removed anyway: create dummy matrices for defgrd and gl strain
      Core::LinAlg::Matrix<3, 3> defgrd = Core::LinAlg::identity_matrix<3>();
      Core::LinAlg::Matrix<6, 1> glstrain(true);

      params.set("temperature", NT(0));
      ctemp = thermoSolidMaterial->evaluate_d_stress_d_scalar(
          defgrd, glstrain, params, map_my_gp_to_so_hex8(gp), id());
    }
    // get thermal material tangent
    else
    {
      compute_ctemp(&ctemp, params);
    }
    // end of call material law

    double detJ_w = detJ * intpoints_.weight(gp);
    // update linear coupling matrix K_dT
    if (stiffmatrix_kdT != nullptr)
    {
      // C_T . N_T
      Core::LinAlg::Matrix<numstr_, nen_> cn(false);
      cn.multiply_nt(ctemp, shapefunct);  // (6x8)=(6x1)(1x8)
      // integrate stiffness term
      // k_dT = k_dT + (B^T . C_T . N_T) * detJ * w(gp)
      stiffmatrix_kdT->multiply_tn(detJ_w, boplin, cn, 1.0);

      // in case of temperature-dependent Young's modulus, additional term for
      // coupling stiffness matrix k_dT
      {
        // k_dT += B_d^T . dC/dT  . B_d . d . N_T
        stiffmatrix_kdT->multiply_nt(detJ_w, Bstress_T, shapefunct, 1.0);

        // k_dT += B_d^T . dC_T/dT  . N_T . T . N_T
        // (24x8)                          (24x1)        (8x1)
        stiffmatrix_kdT->multiply_nt(detJ_w, Bcouplstress_T, shapefunct, 1.0);
      }

      // Be careful: scaling with time factor is done in tsi_monolithic!!
    }  // (stiffmatrix_kdT != nullptr)
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}  // lin_kdT_tsi()


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 03/10 |
 | originally by maf 04/07  (private)                                   |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::nln_stifffint_tsi(
    Core::Elements::LocationArray& la,         // location array
    Core::FE::Discretization& discretization,  ///< discretisation to extract knot vector
    std::vector<double>& disp,                 // current displacements
    std::vector<double>& temp,                 // current temperature
    Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*
        stiffmatrix,                                        // element stiffness matrix
    Core::LinAlg::Matrix<numdofperelement_, 1>* force,      // element internal force vector
    Core::LinAlg::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
    Teuchos::ParameterList& params,                         // algorithmic parameters e.g. time
    const Inpar::Solid::StressType iostress                 // stress output option
)
{
  // update element geometry hex8, 3D: (8x3)
  Core::LinAlg::Matrix<nen_, nsd_> xrefe;  // X, material coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurr;  // x, current coord. of element
  // vector of the current element temperatures
  Core::LinAlg::Matrix<nen_, 1> etemp;

  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * numdofpernode_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * numdofpernode_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * numdofpernode_ + 2];

    etemp(i, 0) = temp[i + 0];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  Core::LinAlg::Matrix<nsd_, nen_> N_XYZ(false);
  // build deformation gradient w.r.t. to material configuration
  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(false);
  // shape functions and their first derivatives
  Core::LinAlg::Matrix<nen_, 1> shapefunct(false);
  Core::LinAlg::Matrix<nsd_, nen_> deriv(false);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);


    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.multiply(invJ_[gp], deriv);  // (6.21)
    double detJ = detJ_[gp];           // (6.22)

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.multiply_tt(xcurr, N_XYZ);

    // right Cauchy-Green tensor = F^T . F
    Core::LinAlg::Matrix<nsd_, nsd_> cauchygreen(false);
    cauchygreen.multiply_tn(defgrd, defgrd);

    // inverse of right Cauchy-Green tensor = F^{-1} . F^{-T}
    Core::LinAlg::Matrix<nsd_, nsd_> Cinv(false);
    Cinv.invert(cauchygreen);
    Core::LinAlg::Matrix<numstr_, 1> Cinv_vct(false);
    Cinv_vct(0) = Cinv(0, 0);
    Cinv_vct(1) = Cinv(1, 1);
    Cinv_vct(2) = Cinv(2, 2);
    Cinv_vct(3) = Cinv(0, 1);
    Cinv_vct(4) = Cinv(1, 2);
    Cinv_vct(5) = Cinv(2, 0);

    // calculate linear B-operator
    Core::LinAlg::Matrix<numstr_, numdofperelement_> boplin(false);
    calculate_boplin(&boplin, &N_XYZ);

    // calculate nonlinear B-operator
    Core::LinAlg::Matrix<numstr_, numdofperelement_> bop(false);
    calculate_bop(&bop, &defgrd, &N_XYZ);

    // temperature
    // described as a matrix (for stress calculation): NT = N_T . T
    Core::LinAlg::Matrix<1, 1> NT(false);
    NT.multiply_tn(shapefunct, etemp);
    // scalar-valued current element temperature T_{n+1}
    // temperature-dependent material parameters, i.e. E(T), pass T_{n+1}
    // insert T_{n+1} into parameter list
    params.set<double>("temperature", NT(0, 0));

    // call material law

    // calculate the stress part dependent on the temperature in the material
    Core::LinAlg::Matrix<numstr_, 1> ctemp(true);
    Core::LinAlg::Matrix<numstr_, 1> couplstress(true);
    Core::LinAlg::Matrix<numstr_, numstr_> cmat_T(true);
    Core::LinAlg::Matrix<numstr_, 1> glstrain(true);

    // insert strain increment into parameter list
    // calculate iterative strains
    Core::LinAlg::Matrix<numstr_, 1> straininc(true);
    params.set<Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>>("straininc", straininc);
    // insert matrices into parameter list which are only required for thrplasthyperelast
    params.set<Core::LinAlg::Matrix<nsd_, nsd_>>("defgrd", defgrd);
    params.set<Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>>("Cinv_vct", Cinv_vct);

    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) couplstress = C . Delta T
    // do not call the material for Robinson's material
    if (material()->material_type() != Core::Materials::m_vp_robinson)
      materialize(&couplstress, &ctemp, &NT, &cmat_T, &glstrain, params);

    // end of call material law

    // return gp stresses
    switch (iostress)
    {
      case Inpar::Solid::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        for (int i = 0; i < numstr_; ++i) (*elestress)(gp, i) = couplstress(i);
        break;
      }
      case Inpar::Solid::stress_cauchy:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");

        // push forward of material stress to the spatial configuration
        // sigma = 1/J . F . S_temp . F^T
        Core::LinAlg::Matrix<nsd_, nsd_> cauchycouplstress(false);
        p_k2to_cauchy(&couplstress, &defgrd, &cauchycouplstress);

        (*elestress)(gp, 0) = cauchycouplstress(0, 0);
        (*elestress)(gp, 1) = cauchycouplstress(1, 1);
        (*elestress)(gp, 2) = cauchycouplstress(2, 2);
        (*elestress)(gp, 3) = cauchycouplstress(0, 1);
        (*elestress)(gp, 4) = cauchycouplstress(1, 2);
        (*elestress)(gp, 5) = cauchycouplstress(0, 2);
        break;
      }
      case Inpar::Solid::stress_none:
        break;
      default:
        FOUR_C_THROW("requested stress type not available");
        break;
    }

    // integrate internal force vector r_d
    // f = f + (B^T . sigma_temp) . detJ . w(gp)
    double detJ_w = detJ * intpoints_.weight(gp);
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) . detJ . w(gp)
      force->multiply_tn(detJ_w, bop, couplstress, 1.0);
    }

    // update stiffness matrix k_dd
    if (stiffmatrix != nullptr)
    {
      // integrate temperature-dependent `elastic' and `initial-displacement'
      // stiffness matrix
      // keu = keu + (B^T . Cmat_T . B) . detJ . w(gp)
      // Neo-Hookean type: dC_T/dd = m . Delta T . (-1) . ( Cinv boeppel Cinv )_{abcd}
      // St.Venant Kirchhoff: dC_T/dd == 0
      // with ( Cinv boeppel Cinv )_{abcd} = 1/2 . ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )
      Core::LinAlg::Matrix<numstr_, numdofperelement_> cb(false);  // cb(6x24) // cmattemp (6x6)
      cb.multiply(cmat_T, bop);
      stiffmatrix->multiply_tn(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************

      // kgeo += ( B_L^T . B_L . sigma_temp) . detJ . w(gp)
      // (B_L^T . sigma . B_L) = (24x6)(6x1)(6x24)
      // --> size of matrices do not fit --> multiply component-by-component
      // with linear B-operator B_L = Ni,Xj, see NiliFEM-Skript (6.20)

      Core::LinAlg::Matrix<numstr_, 1> sfac(couplstress);  // auxiliary integrated stress
      // detJ . w(gp) . [S11,S22,S33,S12=S21,S23=S32,S13=S31]
      sfac.scale(detJ_w);
      // intermediate sigma_temp . B_L (6x1).(6x24)
      std::vector<double> StempB_L(3);
      for (int inod = 0; inod < nen_; ++inod)
      {
        // (3x1) = (6x1) (6x24)
        // S11*N_XYZ(1,i)+S23*N_XYZ(2,i)+S12*N_XYZ(3,i)
        StempB_L[0] =
            sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        // S23*N_XYZ(1,i)+S22*N_XYZ(2,i)+S13*N_XYZ(3,i)
        StempB_L[1] =
            sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        // S12*N_XYZ(1,i)+S13*N_XYZ(2,i)+S33*N_XYZ(3,i)
        StempB_L[2] =
            sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        // (B_L^T . sigma . B_L) = (24x6)(6x24)
        for (int jnod = 0; jnod < nen_; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < nsd_; ++idim)
          {
            // double     (3x8)                 (3x1)
            bopstrbop += N_XYZ(idim, jnod) * StempB_L[idim];
          }
          // (24x24)
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************
    }    // fill k_dd
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}  // nln_stifffint_tsi()


/*----------------------------------------------------------------------*
 | evaluate only the mechanical-thermal stiffness term       dano 11/12 |
 | for monolithic TSI, contribution to k_dT (private)                   |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::nln_kd_t_tsi(Core::Elements::LocationArray& la,
    Core::FE::Discretization& discretization,  ///< discretisation to extract knot vector
    std::vector<double>& disp,                 // current displacement
    std::vector<double>& temp,                 // current temperature
    Core::LinAlg::Matrix<numdofperelement_, nen_>* stiffmatrix_kdT,  // (nsd_*nen_ x nen_)
    Teuchos::ParameterList& params)
{
  // update element geometry (8x3)
  Core::LinAlg::Matrix<nen_, nsd_> xrefe(false);  // X, material coord. of element
  Core::LinAlg::Matrix<nen_, nsd_> xcurr(false);  // x, current  coord. of element
  for (int i = 0; i < nen_; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * numdofpernode_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * numdofpernode_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * numdofpernode_ + 2];
  }

  // ------------------------------------------------ initialise material
  // get the thermal material tangent
  Core::LinAlg::Matrix<numstr_, 1> ctemp(true);

  // ------------------------------ temperature-dependent Young's modulus
  // if young's modulus is temperature-dependent, E(T) additional terms arise
  // for the stiffness matrix k_dT
  Core::LinAlg::Matrix<nen_, 1> etemp(true);

  // get the temperature vector
  for (int i = 0; i < nen_; ++i) etemp(i, 0) = temp[i];

  // initialise matrices required for k_dT
  Core::LinAlg::Matrix<numdofperelement_, 1> Bstress_T(true);
  Core::LinAlg::Matrix<numdofperelement_, 1> Bcouplstress_T(true);

  // shape functions and their first derivatives
  Core::LinAlg::Matrix<nen_, 1> shapefunct(false);
  Core::LinAlg::Matrix<nsd_, nen_> deriv(false);
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  Core::LinAlg::Matrix<nsd_, nen_> N_XYZ(false);
  // build deformation gradient w.r.t. to material configuration
  Core::LinAlg::Matrix<nsd_, nsd_> defgrd(false);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    Core::FE::shape_function<distype>(xsi_[gp], shapefunct);
    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 . N_rst
    N_XYZ.multiply(invJ_[gp], deriv);  // (6.21)
    double detJ = detJ_[gp];           // (6.22)

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T . N_XYZ^T
    defgrd.multiply_tt(xcurr, N_XYZ);

    // right Cauchy-Green tensor = F^T . F
    Core::LinAlg::Matrix<3, 3> cauchygreen;
    cauchygreen.multiply_tn(defgrd, defgrd);
    // initialise inverse of right Cauchy-Green tensor = F^{-1} . F^{-T}
    Core::LinAlg::Matrix<numstr_, 1> Cinv_vct(true);

    // calculate nonlinear B-operator
    Core::LinAlg::Matrix<numstr_, numdofperelement_> bop(false);
    calculate_bop(&bop, &defgrd, &N_XYZ);

    // call material law

    std::shared_ptr<Mat::Trait::ThermoSolid> thermoSolidMaterial =
        std::dynamic_pointer_cast<Mat::Trait::ThermoSolid>(material());
    if (thermoSolidMaterial != nullptr)
    {
      Core::LinAlg::Matrix<numstr_, 1> glstrain(false);
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);  // Voigt notation
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);

      Core::LinAlg::Matrix<1, 1> NT(false);
      NT.multiply_tn(shapefunct, etemp);  // (1x1)
      thermoSolidMaterial->reinit(nullptr, &glstrain, NT(0), map_my_gp_to_so_hex8(gp));
      // full thermal derivative of stress wrt to scalar temperature (needs to be post-multiplied
      // with shape functions)
      params.set("temperature", NT(0));
      ctemp = thermoSolidMaterial->evaluate_d_stress_d_scalar(
          defgrd, glstrain, params, map_my_gp_to_so_hex8(gp), id());
    }
    else if (material()->material_type() == Core::Materials::m_thermoplhyperelast)
    {
      // inverse of Right Cauchy-Green tensor = F^{-1} . F^{-T}
      Core::LinAlg::Matrix<nsd_, nsd_> Cinv(false);
      Cinv.invert(cauchygreen);
      Cinv_vct(0) = Cinv(0, 0);
      Cinv_vct(1) = Cinv(1, 1);
      Cinv_vct(2) = Cinv(2, 2);
      Cinv_vct(3) = Cinv(0, 1);
      Cinv_vct(4) = Cinv(1, 2);
      Cinv_vct(5) = Cinv(2, 0);

      params.set<Core::LinAlg::Matrix<nsd_, nsd_>>("defgrd", defgrd);
      params.set<Core::LinAlg::Matrix<numstr_, 1>>("Cinv_vct", Cinv_vct);

      std::shared_ptr<Mat::ThermoPlasticHyperElast> thermoplhyperelast =
          std::dynamic_pointer_cast<Mat::ThermoPlasticHyperElast>(material());
      // get thermal material tangent
      thermoplhyperelast->setup_cthermo(ctemp, defgrd.determinant(), Cinv_vct);
    }  // m_thermoplhyperelast

    // get thermal material tangent
    else
    {
      compute_ctemp(&ctemp, params);
    }

    // end of call material law

    double detJ_w = detJ * intpoints_.weight(gp);
    // update linear coupling matrix K_dT
    if (stiffmatrix_kdT != nullptr)
    {
      // C_temp . N_temp
      Core::LinAlg::Matrix<numstr_, nen_> cn(false);
      cn.multiply_nt(ctemp, shapefunct);  // (6x8)=(6x1)(1x8)
      // integrate stiffness term
      // k_dT = k_dT + (B^T . C_T . N_temp) . detJ . w(gp)
      stiffmatrix_kdT->multiply_tn(detJ_w, bop, cn, 1.0);

      // in case of temperature-dependent Young's modulus, additional term for
      // coupling stiffness matrix k_dT
      if ((material()->material_type() == Core::Materials::m_thermostvenant))
      {
        // k_dT += B_d^T . stress_T . N_T
        stiffmatrix_kdT->multiply_nt(detJ_w, Bstress_T, shapefunct, 1.0);

        // k_dT += B_d^T . couplstress_T . N_T
        // (24x8)                          (24x1)        (8x1)
        stiffmatrix_kdT->multiply_nt(detJ_w, Bcouplstress_T, shapefunct, 1.0);

        // Be careful: scaling with time factor is done in tsi_monolithic!!
      }  // m_thermostvenant
    }
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

}  // nln_kdT_tsi()


/*----------------------------------------------------------------------*
 | material law with temperature part for So3_thermo         dano 05/10 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::materialize(
    Core::LinAlg::Matrix<numstr_, 1>* couplstress,  // temperature-dependent stress part
    Core::LinAlg::Matrix<numstr_, 1>* ctemp,        // temperature-dependent material tangent
    Core::LinAlg::Matrix<1, 1>* Ntemp,              // temperature of element
    Core::LinAlg::Matrix<numstr_, numstr_>* cmat,   // (mechanical) material tangent
    Core::LinAlg::Matrix<numstr_, 1>* glstrain,     // Green-Lagrange strain tensor
    Teuchos::ParameterList& params                  // parameter
)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!couplstress) FOUR_C_THROW("No stress vector supplied");
#endif

  // backwards compatibility: not all materials use the new interface
  std::shared_ptr<Mat::Trait::ThermoSolid> thermoSolidMaterial =
      std::dynamic_pointer_cast<Mat::Trait::ThermoSolid>(material());
  if (thermoSolidMaterial != nullptr)
  {
    // new interface already includes temperature terms in stress evaluation
    // no additive splitting into stress and couplstress -> nothing to do here
    return;
  }

  // All materials that have a pure Core::LinAlg::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  std::shared_ptr<Core::Mat::Material> mat = material();
  switch (mat->material_type())
  {
    // small strain von Mises thermoelastoplastic material
    case Core::Materials::m_thermopllinelast:
    {
      std::shared_ptr<Mat::ThermoPlasticLinElast> thrpllinelast =
          std::dynamic_pointer_cast<Mat::ThermoPlasticLinElast>(material());
      thrpllinelast->evaluate(*Ntemp, *ctemp, *couplstress);
      return;
      break;
    }
    // visco-plastic Robinson's material
    case Core::Materials::m_vp_robinson:
    {
      // no temperature-dependent stress terms
      return;
      break;
    }
    // thermo-hyperelasto-plastic material
    case Core::Materials::m_thermoplhyperelast:
    {
      std::shared_ptr<Mat::ThermoPlasticHyperElast> thermoplhyperelast =
          std::dynamic_pointer_cast<Mat::ThermoPlasticHyperElast>(material());
      thermoplhyperelast->evaluate(*Ntemp, *ctemp, *cmat, *couplstress, params);
      return;
      break;
    }

    default:
      FOUR_C_THROW("Unknown type of temperature dependent material");
      break;
  }  // switch (mat->material_type())
}  // Materialize()


/*----------------------------------------------------------------------*
 | get the constant temperature fraction for couplstress      dano 05/10 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::compute_ctemp(
    Core::LinAlg::Matrix<numstr_, 1>* ctemp, Teuchos::ParameterList& params)
{
  switch (material()->material_type())
  {
    // thermo st.venant-kirchhoff-material
    case Core::Materials::m_thermostvenant:
    {
      FOUR_C_THROW("Ctemp call no longer valid for ThermoStVenantKirchhoff. Fix implementation.");
    }
    // small strain von Mises thermoelastoplastic material
    case Core::Materials::m_thermopllinelast:
    {
      std::shared_ptr<Mat::ThermoPlasticLinElast> thrpllinelast =
          std::dynamic_pointer_cast<Mat::ThermoPlasticLinElast>(material());
      return thrpllinelast->setup_cthermo(*ctemp);
      break;
    }
    // visco-plastic Robinson's material
    case Core::Materials::m_vp_robinson:
    {
      // so far: do nothing, because the displacement-dependent coupling term
      // is neglected
      return;
      break;
    }
    // thermo-hyperelasto-plastic material
    case Core::Materials::m_thermoplhyperelast:
    {
      std::shared_ptr<Mat::ThermoPlasticHyperElast> thermoplhyperelast =
          std::dynamic_pointer_cast<Mat::ThermoPlasticHyperElast>(material());
      Core::LinAlg::Matrix<3, 3> defgrd = params.get<Core::LinAlg::Matrix<3, 3>>("defgrd");
      Core::LinAlg::Matrix<6, 1> Cinv = params.get<Core::LinAlg::Matrix<6, 1>>("Cinv_vct");
      thermoplhyperelast->setup_cthermo(*ctemp, defgrd.determinant(), Cinv);
      return;
      break;
    }
    default:
      FOUR_C_THROW("Cannot ask material for the temperature-dependent material tangent");
      break;
  }  // switch (mat->material_type())

}  // Ctemp()


/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                        dano 11/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::calculate_bop(
    Core::LinAlg::Matrix<numstr_, numdofperelement_>* bop, Core::LinAlg::Matrix<nsd_, nsd_>* defgrd,
    Core::LinAlg::Matrix<nsd_, nen_>* N_XYZ)
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
 | calculate the nonlinear B-operator in vector notation     dano 11/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::calculate_bop_vec(
    Core::LinAlg::Matrix<1, numdofperelement_>& bopvec, Core::LinAlg::Matrix<nsd_, nsd_>& defgrd,
    Core::LinAlg::Matrix<nsd_, nen_>& N_XYZ)
{
  // ---------------------------------------------- build F as vector 9x1
  // F != F^T, i.e. Voigt notation NOT admissible
  // F (3x3) --> (9x1)
  Core::LinAlg::Matrix<nsd_ * nsd_, 1> defgrd_vec(false);
  defgrd_vec(0) = defgrd(0, 0);
  defgrd_vec(1) = defgrd(0, 1);
  defgrd_vec(2) = defgrd(0, 2);
  defgrd_vec(3) = defgrd(1, 0);
  defgrd_vec(4) = defgrd(1, 1);
  defgrd_vec(5) = defgrd(1, 2);
  defgrd_vec(6) = defgrd(2, 0);
  defgrd_vec(7) = defgrd(2, 1);
  defgrd_vec(8) = defgrd(2, 2);

  // ------------------------ build N_X operator (w.r.t. material config)
  // N_XYZ (3x8) --> 9x24
  Core::LinAlg::Matrix<9, numdofperelement_> N_X(true);  // set to zero
  for (int i = 0; i < nen_; ++i)
  {
    N_X(0, 3 * i + 0) = N_XYZ(0, i);
    N_X(1, 3 * i + 1) = N_XYZ(0, i);
    N_X(2, 3 * i + 2) = N_XYZ(0, i);

    N_X(3, 3 * i + 0) = N_XYZ(1, i);
    N_X(4, 3 * i + 1) = N_XYZ(1, i);
    N_X(5, 3 * i + 2) = N_XYZ(1, i);

    N_X(6, 3 * i + 0) = N_XYZ(2, i);
    N_X(7, 3 * i + 1) = N_XYZ(2, i);
    N_X(8, 3 * i + 2) = N_XYZ(2, i);
  }
  bopvec.multiply_tn(1.0, defgrd_vec, N_X);

}  // Calculate (1x24) B-Operator


/*----------------------------------------------------------------------*
 | calculate the linear B-operator                           dano 11/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::calculate_boplin(
    Core::LinAlg::Matrix<numstr_, numdofperelement_>* boplin,
    Core::LinAlg::Matrix<nsd_, nen_>* N_XYZ)
{
  // lump mass matrix
  if (boplin != nullptr)
  {
    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    for (int i = 0; i < nen_; ++i)
    {
      (*boplin)(0, numdofpernode_ * i + 0) = (*N_XYZ)(0, i);
      (*boplin)(0, numdofpernode_ * i + 1) = 0.0;
      (*boplin)(0, numdofpernode_ * i + 2) = 0.0;
      (*boplin)(1, numdofpernode_ * i + 0) = 0.0;
      (*boplin)(1, numdofpernode_ * i + 1) = (*N_XYZ)(1, i);
      (*boplin)(1, numdofpernode_ * i + 2) = 0.0;
      (*boplin)(2, numdofpernode_ * i + 0) = 0.0;
      (*boplin)(2, numdofpernode_ * i + 1) = 0.0;
      (*boplin)(2, numdofpernode_ * i + 2) = (*N_XYZ)(2, i);
      /* ~~~ */
      (*boplin)(3, numdofpernode_ * i + 0) = (*N_XYZ)(1, i);
      (*boplin)(3, numdofpernode_ * i + 1) = (*N_XYZ)(0, i);
      (*boplin)(3, numdofpernode_ * i + 2) = 0.0;
      (*boplin)(4, numdofpernode_ * i + 0) = 0.0;
      (*boplin)(4, numdofpernode_ * i + 1) = (*N_XYZ)(2, i);
      (*boplin)(4, numdofpernode_ * i + 2) = (*N_XYZ)(1, i);
      (*boplin)(5, numdofpernode_ * i + 0) = (*N_XYZ)(2, i);
      (*boplin)(5, numdofpernode_ * i + 1) = 0.0;
      (*boplin)(5, numdofpernode_ * i + 2) = (*N_XYZ)(0, i);
    }
  }
}  // calculate_boplin()


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::p_k2to_cauchy(
    Core::LinAlg::Matrix<numstr_, 1>* stress, Core::LinAlg::Matrix<nsd_, nsd_>* defgrd,
    Core::LinAlg::Matrix<nsd_, nsd_>* cauchystress)
{
  // calculate the Jacobi-deterinant
  const double detF = (*defgrd).determinant();

  // sigma = 1/J . F . S . F^T
  Core::LinAlg::Matrix<nsd_, nsd_> pkstress(false);
  pkstress(0, 0) = (*stress)(0);
  pkstress(0, 1) = (*stress)(3);
  pkstress(0, 2) = (*stress)(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = (*stress)(1);
  pkstress(1, 2) = (*stress)(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = (*stress)(2);

  Core::LinAlg::Matrix<nsd_, nsd_> temp(false);
  temp.multiply((1.0 / detF), (*defgrd), pkstress);
  (*cauchystress).multiply_nt(temp, (*defgrd));

}  // PK2toCauchy()


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::g_lto_ea(
    Core::LinAlg::Matrix<numstr_, 1>* glstrain, Core::LinAlg::Matrix<nsd_, nsd_>* defgrd,
    Core::LinAlg::Matrix<nsd_, nsd_>* euler_almansi)
{
  // e = F^{T-1} . E . F^{-1}

  // rewrite Green-Lagrange strain in tensor notation
  Core::LinAlg::Matrix<nsd_, nsd_> gl;
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
  Core::LinAlg::Matrix<nsd_, nsd_> invdefgrd;
  invdefgrd.invert(*defgrd);

  // (3x3) = (3x3) (3x3) (3x3)
  Core::LinAlg::Matrix<nsd_, nsd_> temp(false);
  temp.multiply(gl, invdefgrd);
  (*euler_almansi).multiply_tn(invdefgrd, temp);

}  // GLtoEdata()


/*----------------------------------------------------------------------*
 | initialise Jacobian                                       dano 08/12 |
 | is called once in initialize() in so3_thermo_eletypes.cpp            |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::init_jacobian_mapping_special_for_tsi_elements(
    Core::FE::Discretization& dis)
{
  // get the material coordinates
  Core::LinAlg::Matrix<nen_, nsd_> xrefe;
  for (int i = 0; i < nen_; ++i)
  {
    xrefe(i, 0) = nodes()[i]->x()[0];
    xrefe(i, 1) = nodes()[i]->x()[1];
    xrefe(i, 2) = nodes()[i]->x()[2];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  // initialise the derivatives of the shape functions
  Core::LinAlg::Matrix<nsd_, nen_> deriv;
  Core::LinAlg::Matrix<nen_, 1> funct;


  // coordinates of the current integration point (xsi_)
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // get the coordinates of Gauss points, here use intrepid
    const double* gpcoord = intpoints_.point(gp);
    for (int idim = 0; idim < nsd_; idim++)
    {
      xsi_[gp](idim) = gpcoord[idim];
    }
    // first derivatives of shape functions (deriv)

    Core::FE::shape_function_deriv1<distype>(xsi_[gp], deriv);

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
    // derivatives of coordinates w.r.t. material coordinates xjm_ = dx/ds
    invJ_[gp].multiply(deriv, xrefe);
    // xij_ = ds/dx
    detJ_[gp] = invJ_[gp].invert();
    if (detJ_[gp] < 1.0E-16) FOUR_C_THROW("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", detJ_[gp]);
  }  // end gp loop
}  // init_jacobian_mapping()

/*----------------------------------------------------------------------------------*
 | map the GP ordering as defined in Intrepid to the So_hex8 ordering  proell 05/18 |
 *----------------------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::Elements::So3Thermo<So3Ele, distype>::map_my_gp_to_so_hex8(int myGp)
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      const std::array<int, NUMGPT_SOH8> hex8_map = {6, 7, 5, 4, 2, 3, 1, 0};
      return hex8_map[myGp];
    }
    default:
      return myGp;
  }
}


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE

// --- explicit instantiations --- //
#include "4C_so3_thermo_fwd.hpp"
