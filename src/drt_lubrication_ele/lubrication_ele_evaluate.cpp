/*!
\file lubrication_ele_evaluate.cpp
\brief

\level 3

\maintainer Alexander Seitz

*/

#include "../drt_lubrication_ele/lubrication_ele_action.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lubrication_ele/lubrication_ele.H"
#include "../drt_lubrication_ele/lubrication_ele_factory.H"
#include "../drt_lubrication_ele/lubrication_ele_interface.H"
#include "../drt_lubrication_ele/lubrication_ele_parameter.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Lubrication::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    LocationArray&            la,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{

  // check for the action parameter
  const LUBRICATION::Action action = DRT::INPUT::get<LUBRICATION::Action>(params,"action");
  switch(action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Lubrication element)
    case LUBRICATION::calc_mat_and_rhs:
    {
      return DRT::ELEMENTS::LubricationFactory::ProvideImpl(Shape(),discretization.Name())->Evaluate(
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
    case LUBRICATION::calc_lubrication_coupltang:
    {
      return DRT::ELEMENTS::LubricationFactory::ProvideImpl(Shape(),discretization.Name())->EvaluateEHLMon(
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
    case LUBRICATION::calc_error:
    case LUBRICATION::calc_mean_pressures:
    {
      return DRT::ELEMENTS::LubricationFactory::ProvideImpl(Shape(),discretization.Name())->EvaluateService(
               this,
               params,
               discretization,
               la,
               elemat1,
               elemat2,
               elevec1,
               elevec2,
               elevec3);
      break;
    }
    case LUBRICATION::set_time_parameter:
    case LUBRICATION::set_general_lubrication_parameter:
    // these actions have already been evaluated during element pre-evaluate
      break;
    default:
    {
      dserror("Unknown type of action '%i' for Lubrication", action);
      break;
    }
  } // switch(action)

  return 0;
} //DRT::ELEMENTS::Lubrication::Evaluate


/*----------------------------------------------------------------------*
 |  dummy                                                   wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Lubrication::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
//    The function is just a dummy. For Lubrication elements, the integration
//    integration of volume Neumann conditions (body forces) takes place
//    in the element. We need it there for potential stabilisation terms! (wirtz)
  return 0;
}

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter             wirtz 10/15 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::LubricationType::PreEvaluate(DRT::Discretization&               dis,
                                            Teuchos::ParameterList&               p,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix1,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector1,
                                            Teuchos::RCP<Epetra_Vector>           systemvector2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector3)
{
  const LUBRICATION::Action action = DRT::INPUT::get<LUBRICATION::Action>(p,"action");

  switch(action)
  {
  case LUBRICATION::set_general_lubrication_parameter:
  {
    DRT::ELEMENTS::LubricationEleParameter::Instance(dis.Name())->SetGeneralParameters(p);

    break;
  }

  case LUBRICATION::set_time_parameter:
  {
    DRT::ELEMENTS::LubricationEleParameter::Instance(dis.Name())->SetTimeParameters(p);

    break;
  }
  default:
    // do nothing in all other cases
    break;
  }

  return;
}

