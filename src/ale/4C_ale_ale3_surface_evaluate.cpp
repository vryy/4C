/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale3.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale3SurfaceImplInterface* Discret::ELEMENTS::Ale3SurfaceImplInterface::Impl(
    Discret::ELEMENTS::Ale3Surface* ele)
{
  switch (ele->Shape())
  {
    case Core::FE::CellType::quad4:
    {
      return Discret::ELEMENTS::Ale3SurfaceImpl<Core::FE::CellType::quad4>::Instance(
          Core::UTILS::SingletonAction::create);
    }
    default:
      FOUR_C_THROW("shape %d (%d nodes) not supported", ele->Shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::Ale3SurfaceImpl<distype>* Discret::ELEMENTS::Ale3SurfaceImpl<distype>::Instance(
    Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::Ale3SurfaceImpl<distype>>(
            new Discret::ELEMENTS::Ale3SurfaceImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Discret::ELEMENTS::Ale3Surface::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  const Ale3::ActionType act = Core::UTILS::GetAsEnum<Ale3::ActionType>(params, "action");

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
        Core::FE::ExtractMyValues(*dispnp, mydispnp, lm);
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
}  // end of Discret::ELEMENTS::Ale3Surface::Evaluate

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Discret::ELEMENTS::Ale3Surface::evaluate_neumann(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
inline void Discret::ELEMENTS::Ale3SurfaceImpl<distype>::ElementNodeNormal(Ale3Surface* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, std::vector<double>& mydispnp)
{
  Core::FE::ElementNodeNormal<distype>(funct_, deriv_, fac_, unitnormal_, drs_, xsi_, xyze_, ele,
      discretization, elevec1, mydispnp, false, true);
}

FOUR_C_NAMESPACE_CLOSE
