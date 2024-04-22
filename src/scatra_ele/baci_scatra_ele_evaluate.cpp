/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation methods of scatra element

\level 1


*----------------------------------------------------------------------*/

#include "baci_lib_discret.hpp"
#include "baci_mat_elchmat.hpp"
#include "baci_scatra_ele.hpp"
#include "baci_scatra_ele_action.hpp"
#include "baci_scatra_ele_factory.hpp"
#include "baci_scatra_ele_interface.hpp"
#include "baci_scatra_ele_parameter_boundary.hpp"
#include "baci_scatra_ele_parameter_elch.hpp"
#include "baci_scatra_ele_parameter_elch_diffcond.hpp"
#include "baci_scatra_ele_parameter_lsreinit.hpp"
#include "baci_scatra_ele_parameter_std.hpp"
#include "baci_scatra_ele_parameter_timint.hpp"
#include "baci_scatra_ele_parameter_turbulence.hpp"

FOUR_C_NAMESPACE_OPEN



/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter                         |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportType::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const auto action = Teuchos::getIntegralValue<SCATRA::Action>(p, "action");

  switch (action)
  {
    case SCATRA::Action::set_general_scatra_parameter:
    {
      ScaTraEleParameterStd::Instance(dis.Name())->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_nodeset_parameter:
    {
      ScaTraEleParameterStd::Instance(dis.Name())->SetNodesetParameters(p);

      break;
    }

    case SCATRA::Action::set_turbulence_scatra_parameter:
    {
      ScaTraEleParameterTurbulence::Instance(dis.Name())->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_time_parameter:
    {
      ScaTraEleParameterTimInt::Instance(dis.Name())->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_mean_Cai:
    {
      ScaTraEleParameterTurbulence::Instance(dis.Name())->SetCsgsPhi(p.get<double>("meanCai"));

      break;
    }

    case SCATRA::Action::set_lsreinit_scatra_parameter:
    {
      // set general parameters first
      ScaTraEleParameterStd::Instance(dis.Name())->SetParameters(p);

      // set additional, problem-dependent parameters
      ScaTraEleParameterLsReinit::Instance(dis.Name())->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_elch_scatra_parameter:
    {
      // set general parameters first
      ScaTraEleParameterStd::Instance(dis.Name())->SetParameters(p);

      // set additional, problem-dependent parameters
      ScaTraEleParameterElch::Instance(dis.Name())->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_scatra_ele_boundary_parameter:
    {
      // set additional, problem-dependent parameters
      ScaTraEleParameterBoundary::Instance("scatra")->SetParameters(p);

      break;
    }

    case SCATRA::Action::set_diffcond_scatra_parameter:
    {
      // set general parameters first
      ScaTraEleParameterStd::Instance(dis.Name())->SetParameters(p);

      // set additional, problem-dependent parameters
      ScaTraEleParameterElch::Instance(dis.Name())->SetParameters(p);
      ScaTraEleParameterElchDiffCond::Instance(dis.Name())->SetParameters(p);

      break;
    }

    default:
      // do nothing in all other cases
      break;
  }
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 01/09|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  FOUR_C_THROW("not implemented. Use the Evaluate() method with Location Array instead!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 01/09|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  // perform additional operations specific to implementation type
  switch (impltype_)
  {
    case INPAR::SCATRA::impltype_elch_diffcond:
    case INPAR::SCATRA::impltype_elch_diffcond_multiscale:
    case INPAR::SCATRA::impltype_elch_diffcond_thermo:
    case INPAR::SCATRA::impltype_elch_electrode:
    case INPAR::SCATRA::impltype_elch_electrode_growth:
    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    case INPAR::SCATRA::impltype_elch_NP:
    case INPAR::SCATRA::impltype_elch_scl:
    {
      // adapt number of transported scalars for electrochemistry problems
      numscal -= 1;

      // get the material of the first element
      // we assume here, that the material is equal for all elements in this discretization
      Teuchos::RCP<MAT::Material> material = Material();
      if (material->MaterialType() == INPAR::MAT::m_elchmat)
      {
        const auto* actmat = static_cast<const MAT::ElchMat*>(material.get());

        numscal = actmat->NumScal();
      }

      break;
    }
    case INPAR::SCATRA::impltype_levelset:
    case INPAR::SCATRA::impltype_lsreinit:
    {
      // decide whether reinitialization is active or not
      if (not params.get<bool>("solve reinit eq", false))
        impltype_ = INPAR::SCATRA::impltype_levelset;
      else
        impltype_ = INPAR::SCATRA::impltype_lsreinit;

      break;
    }

    case INPAR::SCATRA::impltype_std:
    case INPAR::SCATRA::impltype_thermo_elch_electrode:
    case INPAR::SCATRA::impltype_thermo_elch_diffcond:
    case INPAR::SCATRA::impltype_advreac:
    case INPAR::SCATRA::impltype_refconcreac:
    case INPAR::SCATRA::impltype_chemo:
    case INPAR::SCATRA::impltype_chemoreac:
    case INPAR::SCATRA::impltype_aniso:
    case INPAR::SCATRA::impltype_cardiac_monodomain:
    case INPAR::SCATRA::impltype_loma:
    case INPAR::SCATRA::impltype_poro:
    case INPAR::SCATRA::impltype_pororeac:
    case INPAR::SCATRA::impltype_pororeacECM:
    case INPAR::SCATRA::impltype_multipororeac:
    case INPAR::SCATRA::impltype_one_d_artery:
    case INPAR::SCATRA::impltype_no_physics:
      // do nothing in these cases
      break;

    default:
    {
      // other implementation types are invalid
      FOUR_C_THROW("Invalid implementation type!");
      break;
    }
  }

  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<SCATRA::Action>(params, "action");
  switch (action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Transport element)
    case SCATRA::Action::calc_mat_and_rhs:
    case SCATRA::Action::calc_rhs:
    {
      return ScaTraFactory::ProvideImpl(
          Shape(), impltype_, numdofpernode, numscal, discretization.Name())
          ->Evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }

    case SCATRA::Action::calc_scatra_mono_odblock_fluid:
    case SCATRA::Action::calc_scatra_mono_odblock_mesh:
    case SCATRA::Action::calc_scatra_mono_odblock_scatrathermo:
    case SCATRA::Action::calc_scatra_mono_odblock_thermoscatra:
    {
      return ScaTraFactory::ProvideImpl(
          Shape(), impltype_, numdofpernode, numscal, discretization.Name())
          ->EvaluateOD(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }

    case SCATRA::Action::check_scatra_element_parameter:
    case SCATRA::Action::calc_initial_time_deriv:
    case SCATRA::Action::integrate_shape_functions:
    case SCATRA::Action::calc_flux_domain:
    case SCATRA::Action::calc_total_and_mean_scalars:
    case SCATRA::Action::calc_mean_scalar_time_derivatives:
    case SCATRA::Action::calc_domain_and_bodyforce:
    case SCATRA::Action::calc_scatra_box_filter:
    case SCATRA::Action::calc_turbulent_prandtl_number:
    case SCATRA::Action::calc_vreman_scatra:
    case SCATRA::Action::calc_subgrid_diffusivity_matrix:
    case SCATRA::Action::get_material_parameters:
    case SCATRA::Action::calc_mean_Cai:
    case SCATRA::Action::calc_dissipation:
    case SCATRA::Action::calc_mat_and_rhs_lsreinit_correction_step:
    case SCATRA::Action::calc_node_based_reinit_velocity:
    case SCATRA::Action::calc_domain_integral:
    case SCATRA::Action::calc_error:
    case SCATRA::Action::calc_elch_conductivity:
    case SCATRA::Action::calc_elch_electrode_soc_and_c_rate:
    case SCATRA::Action::calc_elch_elctrode_mean_concentration:
    case SCATRA::Action::calc_elch_domain_kinetics:
    case SCATRA::Action::calc_mass_center_smoothingfunct:
    case SCATRA::Action::get_material_internal_state:
    case SCATRA::Action::set_material_internal_state:
    case SCATRA::Action::get_material_ionic_currents:
    case SCATRA::Action::time_update_material:
    case SCATRA::Action::calc_immersed_element_source:
    case SCATRA::Action::calc_elch_boundary_kinetics_point:
    case SCATRA::Action::micro_scale_initialize:
    case SCATRA::Action::micro_scale_prepare_time_step:
    case SCATRA::Action::micro_scale_solve:
    case SCATRA::Action::micro_scale_update:
    case SCATRA::Action::micro_scale_output:
    case SCATRA::Action::micro_scale_read_restart:
    case SCATRA::Action::micro_scale_set_time:
    case SCATRA::Action::calc_heteroreac_mat_and_rhs:
    case SCATRA::Action::calc_mass_matrix:
    case SCATRA::Action::transform_real_to_reference_point:
    case SCATRA::Action::evaluate_field_in_point:
    {
      return ScaTraFactory::ProvideImpl(
          Shape(), impltype_, numdofpernode, numscal, discretization.Name())
          ->EvaluateService(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case SCATRA::Action::set_general_scatra_parameter:
    case SCATRA::Action::set_turbulence_scatra_parameter:
    case SCATRA::Action::set_time_parameter:
    case SCATRA::Action::set_mean_Cai:
    case SCATRA::Action::set_nodeset_parameter:
    case SCATRA::Action::set_lsreinit_scatra_parameter:
    case SCATRA::Action::set_elch_scatra_parameter:
    case SCATRA::Action::set_scatra_ele_boundary_parameter:
    case SCATRA::Action::set_diffcond_scatra_parameter:
      // these actions have already been evaluated during element pre-evaluate
      break;

    default:
    {
      FOUR_C_THROW("Unknown type of action '%i' for ScaTra", action);
      break;
    }
  }  // switch(action)

  return 0;
}  // DRT::ELEMENTS::Transport::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                        gjb 01/09|
 |                                                                      |
 |  The function is just a dummy. For the transport elements, the       |
 |  integration of the volume Neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilization terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
