/*----------------------------------------------------------------------------*/
/*! \file
\brief evaluate methods of the 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "wall1_scatra.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_mat/material.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::PreEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  const int numnode = NumNode();

  if (la.Size() > 1)
  {
    //  dofs per node of second dofset
    const int numdofpernode = NumDofPerNode(1, *(Nodes()[0]), discretization.Name());

    if (la[1].Size() != numnode * numdofpernode)
      dserror("calc_struct_nlnstiff: Location vector length for velocities does not match!");

    if (discretization.HasState(1, "scalarfield"))
    {
      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState(1, "scalarfield");

      if (phinp == Teuchos::null) dserror("pre_evaluate: Cannot get state vector 'phinp' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myphi =
          Teuchos::rcp(new std::vector<double>(la[1].lm_.size()));
      DRT::UTILS::ExtractMyValues(*phinp, *myphi, la[1].lm_);

      double meanphi = 0.0;
      for (int i = 0; i < numnode; ++i)
      {
        meanphi += (*myphi)[i] / numnode;
      }
      params.set<double>("scalar", meanphi);
    }
    // Get pointer for scatra material in the same element
    Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;
    scatradis = DRT::Problem::Instance()->GetDis("scatra");
    DRT::Element* scatraele = scatradis->gElement(Id());
    Teuchos::RCP<MAT::Material> scatramat =
        Teuchos::rcp_dynamic_cast<MAT::Material>(scatraele->Material());
    params.set<Teuchos::RCP<MAT::Material>>("scatramat", scatramat);
  }
  Teuchos::RCP<std::vector<double>> xrefe = Teuchos::rcp(new std::vector<double>(2));
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnode; ++i)
  {
    const double* x = nodes[i]->X();
    (*xrefe)[0] += x[0] / numnode;
    (*xrefe)[1] += x[1] / numnode;
  }
  params.set<Teuchos::RCP<std::vector<double>>>("position", xrefe);
  return;
}
/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1_Scatra::MyEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1_Scatra::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  SetParamsInterfacePtr(params);

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
      dserror("No action supplied");
    else if (action == "postprocess_stress")
      act = ELEMENTS::struct_postprocess_stress;
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

      PreEvaluate(params, discretization, la);

      Wall1::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // action

  return 0;
}
