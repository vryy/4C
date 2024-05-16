/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale3.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_boundary_integration.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_element_integration_select.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3SurfaceImplInterface* DRT::ELEMENTS::Ale3SurfaceImplInterface::Impl(
    DRT::ELEMENTS::Ale3Surface* ele)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return DRT::ELEMENTS::Ale3SurfaceImpl<CORE::FE::CellType::quad4>::Instance(
          CORE::UTILS::SingletonAction::create);
    }
    default:
      FOUR_C_THROW("shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
      break;
  }
  return nullptr;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::Ale3SurfaceImpl<distype>* DRT::ELEMENTS::Ale3SurfaceImpl<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::Ale3SurfaceImpl<distype>>(
            new DRT::ELEMENTS::Ale3SurfaceImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale3Surface::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  const Ale3::ActionType act = CORE::UTILS::GetAsEnum<Ale3::ActionType>(params, "action");

  switch (act)
  {
    case Ale3::ba_calc_ale_node_normal:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      dispnp = discretization.GetState("dispnp");

      if (dispnp != Teuchos::null)
      {
        mydispnp.resize(lm.size());
        CORE::FE::ExtractMyValues(*dispnp, mydispnp, lm);
      }

      Ale3SurfaceImplInterface::Impl(this)->ElementNodeNormal(
          this, params, discretization, lm, elevec1, mydispnp);

      break;
    }
    default:
      FOUR_C_THROW("Unknown type of action '%i' for Ale3Surface", act);
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::Ale3Surface::Evaluate

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale3Surface::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
inline void DRT::ELEMENTS::Ale3SurfaceImpl<distype>::ElementNodeNormal(Ale3Surface* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, std::vector<double>& mydispnp)
{
  CORE::FE::ElementNodeNormal<distype>(funct_, deriv_, fac_, unitnormal_, drs_, xsi_, xyze_, ele,
      discretization, elevec1, mydispnp, false, true);
}

FOUR_C_NAMESPACE_CLOSE
