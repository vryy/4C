/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


*/

#include "4C_lib_discret.hpp"
#include "4C_lubrication_ele.hpp"
#include "4C_lubrication_ele_action.hpp"
#include "4C_lubrication_ele_factory.hpp"
#include "4C_lubrication_ele_interface.hpp"
#include "4C_lubrication_ele_parameter.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // check for the action parameter
  const LUBRICATION::Action action = Core::UTILS::GetAsEnum<LUBRICATION::Action>(params, "action");
  switch (action)
  {
    // all physics-related stuff is included in the implementation class(es) that can
    // be used in principle inside any element (at the moment: only Lubrication element)
    case LUBRICATION::calc_mat_and_rhs:
    {
      return Discret::ELEMENTS::LubricationFactory::ProvideImpl(Shape(), discretization.Name())
          ->Evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case LUBRICATION::calc_lubrication_coupltang:
    {
      return Discret::ELEMENTS::LubricationFactory::ProvideImpl(Shape(), discretization.Name())
          ->EvaluateEHLMon(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case LUBRICATION::calc_error:
    case LUBRICATION::calc_mean_pressures:
    {
      return Discret::ELEMENTS::LubricationFactory::ProvideImpl(Shape(), discretization.Name())
          ->EvaluateService(
              this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case LUBRICATION::set_time_parameter:
    case LUBRICATION::set_general_lubrication_parameter:
      // these actions have already been evaluated during element pre-evaluate
      break;
    default:
    {
      FOUR_C_THROW("Unknown type of action '%i' for Lubrication", action);
      break;
    }
  }  // switch(action)

  return 0;
}  // Discret::ELEMENTS::Lubrication::Evaluate


/*----------------------------------------------------------------------*
 |  dummy                                                   wirtz 10/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Lubrication::evaluate_neumann(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  //    The function is just a dummy. For Lubrication elements, the integration
  //    integration of volume Neumann conditions (body forces) takes place
  //    in the element. We need it there for potential stabilisation terms! (wirtz)
  return 0;
}

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter             wirtz 10/15 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationType::pre_evaluate(Discret::Discretization& dis,
    Teuchos::ParameterList& p, Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const LUBRICATION::Action action = Core::UTILS::GetAsEnum<LUBRICATION::Action>(p, "action");

  switch (action)
  {
    case LUBRICATION::set_general_lubrication_parameter:
    {
      Discret::ELEMENTS::LubricationEleParameter::Instance(dis.Name())->set_general_parameters(p);

      break;
    }

    case LUBRICATION::set_time_parameter:
    {
      Discret::ELEMENTS::LubricationEleParameter::Instance(dis.Name())->SetTimeParameters(p);

      break;
    }
    default:
      // do nothing in all other cases
      break;
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
