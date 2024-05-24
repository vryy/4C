/*----------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG boundary elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#include "4C_scatra_ele_hdg_boundary_calc.hpp"

#include "4C_lib_node.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_hdg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGBoundaryImplInterface* DRT::ELEMENTS::ScaTraHDGBoundaryImplInterface::Impl(
    const DRT::Element* ele)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::quad4>::Instance();
    }
    case CORE::FE::CellType::quad8:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::quad8>::Instance();
    }
    case CORE::FE::CellType::quad9:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::quad9>::Instance();
    }
    case CORE::FE::CellType::tri3:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::tri3>::Instance();
    }
    case CORE::FE::CellType::tri6:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::tri6>::Instance();
    }
    case CORE::FE::CellType::line2:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::line2>::Instance();
    }
    case CORE::FE::CellType::line3:
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::line3>::Instance();
    }
    case CORE::FE::CellType::nurbs2:  // 1D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::nurbs2>::Instance();
    }
    case CORE::FE::CellType::nurbs3:  // 1D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::nurbs3>::Instance();
    }
    case CORE::FE::CellType::nurbs4:  // 2D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::nurbs4>::Instance();
    }
    case CORE::FE::CellType::nurbs9:  // 2D nurbs boundary element
    {
      return ScaTraHDGBoundaryImpl<CORE::FE::CellType::nurbs9>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>*
DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::Instance(CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>>(
            new DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::ScaTraHDGBoundaryImpl()
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
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::evaluate_neumann(
    DRT::ELEMENTS::ScaTraHDGBoundary* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra)
{
  CORE::LINALG::SerialDenseVector dummy_vec2, dummy_vec3;
  CORE::LINALG::SerialDenseMatrix dummy_mat2;

  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::project_neumann_field, params);

  const int* nodeids = ele->NodeIds();

  DRT::Element* parent = ele->parent_element();
  Teuchos::RCP<DRT::FaceElement>* faces = parent->Faces();
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
