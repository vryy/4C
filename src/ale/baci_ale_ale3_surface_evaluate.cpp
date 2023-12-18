/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_ale3.H"
#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_geometry_position_array.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_lib_discret.H"
#include "baci_lib_element_integration_select.H"
#include "baci_lib_utils.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3Surface_Impl_Interface* DRT::ELEMENTS::Ale3Surface_Impl_Interface::Impl(
    DRT::ELEMENTS::Ale3Surface* ele)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return DRT::ELEMENTS::Ale3Surface_Impl<CORE::FE::CellType::quad4>::Instance(
          CORE::UTILS::SingletonAction::create);
    }
    default:
      dserror("shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
      break;
  }
  return nullptr;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::Ale3Surface_Impl<distype>* DRT::ELEMENTS::Ale3Surface_Impl<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::Ale3Surface_Impl<distype>>(
            new DRT::ELEMENTS::Ale3Surface_Impl<distype>());
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
  const Ale3::ActionType act = DRT::INPUT::get<Ale3::ActionType>(params, "action");

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
        DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);
      }

      Ale3Surface_Impl_Interface::Impl(this)->ElementNodeNormal(
          this, params, discretization, lm, elevec1, mydispnp);

      break;
    }
    default:
      dserror("Unknown type of action '%i' for Ale3Surface", act);
      break;
  }  // end of switch(act)

  return 0;
}  // end of DRT::ELEMENTS::Ale3Surface::Evaluate

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale3Surface::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
inline void DRT::ELEMENTS::Ale3Surface_Impl<distype>::ElementNodeNormal(Ale3Surface* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, std::vector<double>& mydispnp)
{
  CORE::FE::ElementNodeNormal<distype>(funct_, deriv_, fac_, unitnormal_, drs_, xsi_, xyze_, ele,
      discretization, elevec1, mydispnp, false, true);
}

BACI_NAMESPACE_CLOSE
