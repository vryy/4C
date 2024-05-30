/*----------------------------------------------------------------------------*/
/*! \file
\brief evaluate methods of the 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_w1_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Scatra::pre_evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la)
{
  const int numnode = num_node();

  if (la.Size() > 1)
  {
    //  dofs per node of second dofset
    const int numdofpernode = NumDofPerNode(1, *(Nodes()[0]), discretization.Name());

    if (la[1].Size() != numnode * numdofpernode)
      FOUR_C_THROW("calc_struct_nlnstiff: Location vector length for velocities does not match!");

    if (discretization.HasState(1, "scalarfield"))
    {
      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState(1, "scalarfield");

      if (phinp == Teuchos::null) FOUR_C_THROW("pre_evaluate: Cannot get state vector 'phinp' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myphi =
          Teuchos::rcp(new std::vector<double>(la[1].lm_.size()));
      CORE::FE::ExtractMyValues(*phinp, *myphi, la[1].lm_);

      double meanphi = 0.0;
      for (int i = 0; i < numnode; ++i)
      {
        meanphi += (*myphi)[i] / numnode;
      }
      params.set<double>("scalar", meanphi);
    }
    // Get pointer for scatra material in the same element
    Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;
    scatradis = GLOBAL::Problem::Instance()->GetDis("scatra");
    CORE::Elements::Element* scatraele = scatradis->gElement(Id());
    Teuchos::RCP<CORE::MAT::Material> scatramat =
        Teuchos::rcp_dynamic_cast<CORE::MAT::Material>(scatraele->Material());
    params.set<Teuchos::RCP<CORE::MAT::Material>>("scatramat", scatramat);
  }
  CORE::LINALG::Matrix<2, 1> xrefe(true);
  CORE::Nodes::Node** nodes = Nodes();
  for (int i = 0; i < numnode; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(0) += x[0] / numnode;
    xrefe(1) += x[1] / numnode;
  }
  params.set("elecenter_coords_ref", xrefe);
}
/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1Scatra::my_evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1Scatra::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  set_params_interface_ptr(params);

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
    if (action == "none") FOUR_C_THROW("No action supplied");
  }

  // what should the element do
  switch (act)
  {
      //==================================================================================
      // coupling terms in force-vector and stiffness matrix
      //  case Wall1_Scatra::calc_struct_multidofsetcoupling:
      //  {
      //
      //    MyEvaluate(params,
      //                      discretization,
      //                      la,
      //                      elemat1_epetra,
      //                      elemat2_epetra,
      //                      elevec1_epetra,
      //                      elevec2_epetra,
      //                      elevec3_epetra);
      //  }
      //  break;
      //  case Wall1_Scatra::postprocess_stress:
      //  {
      //    Wall1::Evaluate(params,
      //                          discretization,
      //                          la[0].lm_,
      //                          elemat1_epetra,
      //                          elemat2_epetra,
      //                          elevec1_epetra,
      //                          elevec2_epetra,
      //                          elevec3_epetra);
      //  }
      //  break;
    /*case Wall1_Scatra::calc_struct_update_istep:
    {
      so3_ele::Evaluate(params,
                        discretization,
                        la[0].lm_,
                        elemat1_epetra,
                        elemat2_epetra,
                        elevec1_epetra,
                        elevec2_epetra,
                        elevec3_epetra);
    }
    break;*/
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating

      pre_evaluate(params, discretization, la);

      Wall1::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      my_evaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // action

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
