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

#include "scatra_ele_impl_utils.H"

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

  if (action == SCATRA::set_general_scatra_parameter)
  {
    DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara = DRT::ELEMENTS::ScaTraEleParameterStd::Instance();
    scatrapara->SetElementGeneralScaTraParameter(p,dis.Comm().MyPID());
  }
  if (action == SCATRA::set_turbulence_scatra_parameter)
  {
    DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara = DRT::ELEMENTS::ScaTraEleParameterStd::Instance();
    scatrapara->SetElementTurbulenceParameter(p);
  }
  else if (action == SCATRA::set_time_parameter)
  {
    DRT::ELEMENTS::ScaTraEleParameterTimInt* scatrapara = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance();
    scatrapara->SetElementTimeParameter(p);
  }
  else if (action == SCATRA::set_mean_Cai)
  {
    const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(p, "scatratype");
    if (scatratype == INPAR::SCATRA::scatratype_undefined)
      dserror("Element parameter SCATRATYPE has not been set!");
    switch(scatratype)
    {
    case INPAR::SCATRA::scatratype_condif:
    case INPAR::SCATRA::scatratype_loma:
    {
      DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara = DRT::ELEMENTS::ScaTraEleParameterStd::Instance();
      scatrapara->SetCsgsPhi(p.get<double>("meanCai"));
      break;
    }
    default:
        dserror("set your scatratype for mfs here");
        break;
    }
  }
  else if (action == SCATRA::set_lsreinit_scatra_parameter)
  {
    DRT::ELEMENTS::ScaTraEleParameterLsReinit* scatrapara = DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance();
    // set general parameters first
    scatrapara->SetElementGeneralScaTraParameter(p,dis.Comm().MyPID());
    // set additional problem-dependent parameters
    scatrapara->SetElementLsReinitScaTraParameter(p);
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
  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = this->NumDofPerNode(*(this->Nodes()[0]));
  int numscal = numdofpernode;
  if (SCATRA::IsElchProblem(scatratype))
  {
    numscal -= 1;

    // get the material of the first element
    // we assume here, that the material is equal for all elements in this discretization
    Teuchos::RCP<MAT::Material> material = this->Material();
    if (material->MaterialType() == INPAR::MAT::m_elchmat)
    {
      const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(material.get());

      if (actmat->Current())
        numscal -= DRT::UTILS::getDimension(this->Shape());
    }
  }

  // switch between different physical types as used below
  std::string impltype = "std";
  switch(scatratype)
  {
  case INPAR::SCATRA::scatratype_condif:     impltype = "std";     break;
  case INPAR::SCATRA::scatratype_loma:       impltype = "loma";    break;
  case INPAR::SCATRA::scatratype_elch:       impltype = "elch";    break;
  case INPAR::SCATRA::scatratype_levelset:
  {
    if (not params.get<bool>("solve reinit eq",false))
      impltype = "std";
    else
      impltype = "lsreinit";
    break;
  }
  case INPAR::SCATRA::scatratype_poro:       impltype = "poro";    break;
  default: dserror("Unknown scatratype for calc_mat_and_rhs!");    break;
  }

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch(action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Transport element)
    case SCATRA::calc_mat_and_rhs:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(), impltype, numdofpernode, numscal)->Evaluate(
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
    case SCATRA::calc_initial_time_deriv:
    case SCATRA::integrate_shape_functions:
    case SCATRA::calc_flux_domain:
    case SCATRA::calc_mean_scalars:
    case SCATRA::calc_domain_and_bodyforce:
    case SCATRA::calc_scatra_box_filter:
    case SCATRA::calc_turbulent_prandtl_number:
    case SCATRA::calc_vreman_scatra:
    case SCATRA::calc_subgrid_diffusivity_matrix:
    case SCATRA::calc_mean_Cai:
    case SCATRA::calc_dissipation:
    case SCATRA::calc_mat_and_rhs_lsreinit_correction_step:
    case SCATRA::calc_node_based_reinit_velocity:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(), impltype, numdofpernode, numscal)->EvaluateService(
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
      break;
    default:
      dserror("Unknown type of action '%i' for ScaTra", action);
      break;
  } // end of switch(act)

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


