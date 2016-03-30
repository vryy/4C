/*!-------------------------------------------------------------------------*\
 * \file meshfree_scatra_cell_evaluate.cpp
 *
 * \brief evaluation of a meshfree scatra cell
 *
 * <pre>
 * \maintainer Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "meshfree_scatra_cell.H"       // eh klar...
#include "meshfree_scatra_cell_interface.H"       // DRT::ELEMENTS::MeshfreeScaTraImplInterface

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_parameter_std.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../drt_scatra_ele/scatra_ele_factory.H"
#include "../drt_scatra_ele/scatra_ele_calc_utils.H"

#include "../drt_inpar/inpar_scatra.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_mat/material.H"
#include "../drt_mat/elchmat.H"

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter               nis Dec13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeTransportType::PreEvaluate(
  DRT::Discretization&                 dis,
  Teuchos::ParameterList&              p,
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
  Teuchos::RCP<Epetra_Vector>          systemvector1,
  Teuchos::RCP<Epetra_Vector>          systemvector2,
  Teuchos::RCP<Epetra_Vector>          systemvector3)
{
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(p,"action");

  if (action == SCATRA::set_general_scatra_parameter)
    DRT::ELEMENTS::ScaTraEleParameterStd::Instance(dis.Name())->SetParameters(p);
  else if (action == SCATRA::set_time_parameter)
    DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(dis.Name())->SetParameters(p);

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        hoermann 08/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeTransport::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  dserror("not implemented. Use the Evaluate() method with Location Array instead!");
  return -1;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             nis Mar12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeTransport::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    LocationArray&            la,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the csomputations
  const int numdofpernode = this->NumDofPerNode(*(this->Nodes()[0]));
  int numscal = numdofpernode;

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch(action){
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Transport element)
    case SCATRA::calc_mat_and_rhs:
    {
      return DRT::ELEMENTS::ScaTraFactory::ProvideMeshfreeImpl(Shape(), impltype_, numdofpernode, numscal, discretization.Name() )->Evaluate(
              this,
              params,
              discretization,
              la,
              elemat1,
              elemat2,
              elevec1,
              elevec2,
              elevec3
              );
      break;
    }
    case SCATRA::integrate_shape_functions:
//    {
//      return DRT::ELEMENTS::ScaTraFactory::ProvideImpl(Shape(), impltype_, numdofpernode, numscal)->EvaluateService(
//               this,
//               params,
//               discretization,
//               lm,
//               elemat1,
//               elemat2,
//               elevec1,
//               elevec2,
//               elevec3);
//      break;
//    }
    case SCATRA::calc_initial_time_deriv:
    case SCATRA::calc_flux_domain:
    case SCATRA::calc_total_and_mean_scalars:
    case SCATRA::calc_domain_and_bodyforce:
    case SCATRA::calc_scatra_box_filter:
    case SCATRA::calc_turbulent_prandtl_number:
    case SCATRA::calc_vreman_scatra:
    case SCATRA::calc_subgrid_diffusivity_matrix:
    case SCATRA::calc_mean_Cai:
    case SCATRA::calc_dissipation:
    case SCATRA::calc_mat_and_rhs_lsreinit_correction_step:
    case SCATRA::calc_node_based_reinit_velocity:
    case SCATRA::set_general_scatra_parameter:
    case SCATRA::set_turbulence_scatra_parameter:
    case SCATRA::set_time_parameter:
    case SCATRA::set_mean_Cai:
    case SCATRA::set_lsreinit_scatra_parameter:
      break;
    default:
      dserror("Unknown type of action '%i' for MeshfreeScaTra", action);
      break;
  } // end of switch(act)

  return(0);
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                       nis Mar12 |
 |                                                                      |
 |  Not decided on an implementation of Neumann BC for meshfree schemes |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeTransport::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}
