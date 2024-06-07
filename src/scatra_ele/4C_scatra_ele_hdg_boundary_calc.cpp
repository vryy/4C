/*----------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG boundary elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#include "4C_scatra_ele_hdg_boundary_calc.hpp"

#include "4C_fem_general_node.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_hdg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraHDGBoundaryImplInterface*
Discret::ELEMENTS::ScaTraHDGBoundaryImplInterface::Impl(const Core::Elements::Element* ele)
{
  switch (ele->Shape())
  {
    case Core::FE::CellType::quad4:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::quad4>::Instance();
    }
    case Core::FE::CellType::quad8:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::quad8>::Instance();
    }
    case Core::FE::CellType::quad9:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::quad9>::Instance();
    }
    case Core::FE::CellType::tri3:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::tri3>::Instance();
    }
    case Core::FE::CellType::tri6:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::tri6>::Instance();
    }
    case Core::FE::CellType::line2:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::line2>::Instance();
    }
    case Core::FE::CellType::line3:
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::line3>::Instance();
    }
    case Core::FE::CellType::nurbs2:  // 1D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::nurbs2>::Instance();
    }
    case Core::FE::CellType::nurbs3:  // 1D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::nurbs3>::Instance();
    }
    case Core::FE::CellType::nurbs4:  // 2D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::nurbs4>::Instance();
    }
    case Core::FE::CellType::nurbs9:  // 2D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<Core::FE::CellType::nurbs9>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>*
Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::Instance(Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>>(
            new Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::ScaTraHDGBoundaryImpl()
    : xyze_(true),
      funct_(true),
      deriv_(true),
      unitnormal_(true),
      velint_(true),
      drs_(0.0),
      fac_(0.0)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::evaluate_neumann(
    Discret::ELEMENTS::ScaTraHDGBoundary* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra)
{
  Core::LinAlg::SerialDenseVector dummy_vec2, dummy_vec3;
  Core::LinAlg::SerialDenseMatrix dummy_mat2;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::project_neumann_field, params);

  const int* nodeids = ele->NodeIds();

  Core::Elements::Element* parent = ele->parent_element();
  Teuchos::RCP<Core::Elements::FaceElement>* faces = parent->Faces();
  bool same = false;
  for (int i = 0; i < parent->NumFace(); ++i)
  {
    const int* nodeidsfaces = faces[i]->NodeIds();

    if (faces[i]->num_node() != ele->num_node()) break;

    for (int j = 0; j < ele->num_node(); ++j)
    {
      if (nodeidsfaces[j] == nodeids[j])
        same = true;
      else
      {
        same = false;
        break;
      }
    }
    if (same == true)
    {
      // i is the number we were searching for!!!!
      params.set<int>("face", i);
      ele->parent_element()->Evaluate(params, discretization, la, elemat1_epetra, dummy_mat2,
          elevec1_epetra, dummy_vec2, dummy_vec3);
      // break;
    }
  }
  if (same == false && (faces[0]->num_node() != ele->num_node()))
    FOUR_C_THROW("Neumann boundary condition implemented only for surface elements");

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
