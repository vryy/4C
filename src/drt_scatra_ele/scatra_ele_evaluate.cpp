/*!
\file scatra_element_evaluate.cpp
\brief

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/

#include "scatra_ele.H"

#include "scatra_ele_factory.H"
#include "scatra_ele_interface.H"
#include "scatra_ele_action.H"

#include "scatra_ele_parameter.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_lsreinit.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_elch_diffcond.H"

#include "scatra_ele_calc_utils.H"

#include "../drt_inpar/inpar_scatra.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_mat/elchmat.H"

#include "../drt_lib/drt_discret.H"



/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter                         |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportType::PreEvaluate(DRT::Discretization&               dis,
                                            Teuchos::ParameterList&               p,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix1,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector1,
                                            Teuchos::RCP<Epetra_Vector>           systemvector2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector3)
{
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(p,"action");

  switch(action)
  {
  case SCATRA::set_general_scatra_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterStd::Instance()->SetElementGeneralParameters(p);

    break;
  }

  case SCATRA::set_turbulence_scatra_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterStd::Instance()->SetElementTurbulenceParameters(p);

    break;
  }

  case SCATRA::set_time_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance()->SetElementTimeParameter(p);

    break;
  }

  case SCATRA::set_mean_Cai:
  {
    DRT::ELEMENTS::ScaTraEleParameterStd::Instance()->SetCsgsPhi(p.get<double>("meanCai"));

    break;
  }

  case SCATRA::set_lsreinit_scatra_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterLsReinit* scatrapara = DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance();
    // set general parameters first
    scatrapara->SetElementGeneralParameters(p);
    // set additional problem-dependent parameters
    scatrapara->SetElementLsReinitScaTraParameter(p);

    break;
  }

  case SCATRA::set_elch_scatra_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterElch::Instance()->SetElementGeneralParameters(p);
    break;
  }

  case SCATRA::set_diffcond_scatra_parameter:
  {
    DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance()->SetElementGeneralParameters(p);
    break;
  }

  default:
    // do nothing in all other cases
    break;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 01/09|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  // perform additional operations specific to implementation type
  switch(impltype_)
  {
  case INPAR::SCATRA::impltype_elch_diffcond:
  case INPAR::SCATRA::impltype_elch_electrode:
  case INPAR::SCATRA::impltype_elch_NP:
  {
    // adapt number of transported scalars for electrochemistry problems
    numscal -= 1;

    // get the material of the first element
    // we assume here, that the material is equal for all elements in this discretization
    Teuchos::RCP<MAT::Material> material = Material();
    if (material->MaterialType() == INPAR::MAT::m_elchmat)
    {
      const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(material.get());

      numscal = actmat->NumScal();
    }

    break;
  }

  case INPAR::SCATRA::impltype_levelset:
  case INPAR::SCATRA::impltype_lsreinit:
  {
    // decide whether reinitialization is active or not
    if(not params.get<bool>("solve reinit eq",false))
      impltype_ = INPAR::SCATRA::impltype_levelset;
    else
      impltype_ = INPAR::SCATRA::impltype_lsreinit;

    break;
  }

  case INPAR::SCATRA::impltype_std:
  case INPAR::SCATRA::impltype_advreac:
  case INPAR::SCATRA::impltype_aniso:
  case INPAR::SCATRA::impltype_cardiac_monodomain:
  case INPAR::SCATRA::impltype_loma:
  case INPAR::SCATRA::impltype_poro:
  case INPAR::SCATRA::impltype_pororeac:
    // do nothing in these cases
    break;

  default:
  {
    // other implementation types are invalid
    dserror("Invalid implementation type!");
    break;
  }
  }

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch(action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Transport element)
    case SCATRA::calc_mat_and_rhs:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(),impltype_,numdofpernode,numscal)->Evaluate(
              this,
              params,
              discretization,
              lm,
              elemat1,
              elemat2,
              elevec1,
              elevec2,
              elevec3
              );
    }
    break;
    case SCATRA::calc_scatra_mono_odblock_fluid:
    case SCATRA::calc_scatra_mono_odblock_mesh:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(),impltype_,numdofpernode,numscal)->EvaluateOD(
              this,
              params,
              discretization,
              lm,
              elemat1,
              elemat2,
              elevec1,
              elevec2,
              elevec3
              );
    }
    break;
    case SCATRA::check_scatra_element_parameter:
    case SCATRA::calc_initial_time_deriv:
    case SCATRA::integrate_shape_functions:
    case SCATRA::calc_flux_domain:
    case SCATRA::calc_mean_scalars:
    case SCATRA::calc_domain_and_bodyforce:
    case SCATRA::calc_scatra_box_filter:
    case SCATRA::calc_turbulent_prandtl_number:
    case SCATRA::calc_vreman_scatra:
    case SCATRA::calc_subgrid_diffusivity_matrix:
    case SCATRA::get_material_parameters:
    case SCATRA::calc_mean_Cai:
    case SCATRA::calc_dissipation:
    case SCATRA::calc_mat_and_rhs_lsreinit_correction_step:
    case SCATRA::calc_node_based_reinit_velocity:
    case SCATRA::calc_error:
    case SCATRA::calc_elch_conductivity:
    case SCATRA::calc_elch_electrode_soc:
    case SCATRA::calc_integr_grad_reac:
    case SCATRA::recon_gradients_at_nodes:
    case SCATRA::recon_curvature_at_nodes:
    case SCATRA::calc_mass_center_smoothingfunct:
    case SCATRA::get_material_internal_state:
    case SCATRA::set_material_internal_state:
    case SCATRA::get_material_ionic_currents:
    case SCATRA::time_update_material:
    case SCATRA::calc_integr_pat_rhsvec:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(),impltype_,numdofpernode,numscal)->EvaluateService(
               this,
               params,
               discretization,
               lm,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3);
      break;
    }
    case SCATRA::set_general_scatra_parameter:
    case SCATRA::set_turbulence_scatra_parameter:
    case SCATRA::set_time_parameter:
    case SCATRA::set_mean_Cai:
    case SCATRA::set_lsreinit_scatra_parameter:
    case SCATRA::set_elch_scatra_parameter:
    case SCATRA::set_diffcond_scatra_parameter:
      break;
    default:
      dserror("Unknown type of action '%i' for ScaTra", action);
      break;
  } // switch(action)

  return 0;
} //DRT::ELEMENTS::Transport::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                        gjb 01/09|
 |                                                                      |
 |  The function is just a dummy. For the transport elements, the       |
 |  integration of the volume Neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilization terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


