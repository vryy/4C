/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation methods of the porofluidmultiphase element

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_lib_discret.hpp"
#include "4C_porofluidmultiphase_ele.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"
#include "4C_porofluidmultiphase_ele_factory.hpp"
#include "4C_porofluidmultiphase_ele_interface.hpp"
#include "4C_porofluidmultiphase_ele_parameter.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhase::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the element and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));

  // check for the action parameter
  const POROFLUIDMULTIPHASE::Action action =
      CORE::UTILS::GetAsEnum<POROFLUIDMULTIPHASE::Action>(params, "action");
  switch (action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only PoroFluidMultiPhase element)
    case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
    case POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat:
    case POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat:
    case POROFLUIDMULTIPHASE::calc_error:
    case POROFLUIDMULTIPHASE::calc_pres_and_sat:
    case POROFLUIDMULTIPHASE::calc_solidpressure:
    case POROFLUIDMULTIPHASE::calc_porosity:
    case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
    case POROFLUIDMULTIPHASE::calc_phase_velocities:
    case POROFLUIDMULTIPHASE::calc_initial_time_deriv:
    case POROFLUIDMULTIPHASE::calc_valid_dofs:
    case POROFLUIDMULTIPHASE::calc_domain_integrals:
    {
      std::vector<CORE::LINALG::SerialDenseMatrix*> elemat(2);
      elemat[0] = &elemat1;
      elemat[1] = &elemat2;

      std::vector<CORE::LINALG::SerialDenseVector*> elevec(3);
      elevec[0] = &elevec1;
      elevec[1] = &elevec2;
      elevec[2] = &elevec3;

      return DRT::ELEMENTS::PoroFluidMultiPhaseFactory::ProvideImpl(
          Shape(), numdofpernode, discretization.Name())
          ->Evaluate(this, params, discretization, la, elemat, elevec);
      break;
    }
    case POROFLUIDMULTIPHASE::set_timestep_parameter:
    case POROFLUIDMULTIPHASE::set_general_parameter:
      // these actions have already been evaluated during element pre-evaluate
      break;
    default:
    {
      FOUR_C_THROW("Unknown type of action '%i' for PoroFluidMultiPhase", action);
      break;
    }
  }  // switch(action)

  return 0;
}  // DRT::ELEMENTS::PoroFluidMultiPhase::Evaluate


/*----------------------------------------------------------------------*
 |  dummy                                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhase::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("EvaluateNeumann for PoroFluidMultiPhase  not yet implemented!");
  //    The function is just a dummy. For PoroFluidMultiPhase elements, the integration
  //    integration of volume Neumann conditions (body forces) takes place
  //    in the element. We need it there for potential stabilisation terms!
  return 0;
}

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter             vuong 08/16 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseType::PreEvaluate(DRT::Discretization& dis,
    Teuchos::ParameterList& p, Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const POROFLUIDMULTIPHASE::Action action =
      CORE::UTILS::GetAsEnum<POROFLUIDMULTIPHASE::Action>(p, "action");

  switch (action)
  {
    case POROFLUIDMULTIPHASE::set_general_parameter:
    {
      DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(dis.Name())->SetGeneralParameters(p);

      break;
    }

    case POROFLUIDMULTIPHASE::set_timestep_parameter:
    {
      DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(dis.Name())
          ->SetTimeStepParameters(p);

      break;
    }
    default:
      // do nothing in all other cases
      break;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
