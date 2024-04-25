/*-----------------------------------------------------------*/
/*! \file

\brief evaluation routines for the fluid poro boundary element


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_boundary_calc.hpp"
#include "4C_fluid_ele_boundary_factory.hpp"
#include "4C_fluid_ele_poro.hpp"
#include "4C_inpar_fluid.hpp"

FOUR_C_NAMESPACE_OPEN

int DRT::ELEMENTS::FluidPoroBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const auto act = CORE::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");

  // switch between different physical types as used below
  std::string impltype = "poro";
  switch (params.get<int>("Physical Type", INPAR::FLUID::poro))
  {
    case INPAR::FLUID::poro:
      impltype = "poro";
      break;
    case INPAR::FLUID::poro_p1:
      impltype = "poro_p1";
      break;
    default:
      FOUR_C_THROW("invalid physical type for porous fluid!");
      break;
  }

  switch (act)
  {
    case FLD::calc_flowrate:
    case FLD::no_penetration:
    case FLD::no_penetrationIDs:
    case FLD::poro_boundary:
    case FLD::poro_prescoupl:
    case FLD::poro_splitnopenetration:
    case FLD::poro_splitnopenetration_OD:
    case FLD::poro_splitnopenetration_ODdisp:
    case FLD::poro_splitnopenetration_ODpres:
    case FLD::fpsi_coupling:
    {
      DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), impltype)
          ->EvaluateAction(
              this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    default:  // call standard fluid boundary element
    {
      FluidBoundary::Evaluate(
          params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
  }

  return 0;
}

void DRT::ELEMENTS::FluidPoroBoundary::LocationVector(const Discretization& dis, LocationArray& la,
    bool doDirichlet, const std::string& condstring, Teuchos::ParameterList& params) const
{
  // get the action required
  const auto act = CORE::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");
  switch (act)
  {
    case FLD::poro_boundary:
    case FLD::fpsi_coupling:
    case FLD::calc_flowrate:
    case FLD::poro_splitnopenetration_ODdisp:
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      ParentElement()->LocationVector(dis, la, doDirichlet);
      break;
    case FLD::poro_splitnopenetration:
    case FLD::poro_splitnopenetration_OD:
    {
      // This is a hack ...
      // Remove pressure dofs from location vector!!!
      // call standard fluid boundary element
      FluidBoundary::LocationVector(dis, la, doDirichlet, condstring, params);
      int dim = la[0].stride_[0] - 1;
      // extract velocity dofs from first dofset
      for (int i = NumNode(); i > 0; --i)
      {
        la[0].lm_.erase(la[0].lm_.begin() + (i - 1) * (dim + 1) + dim);
        la[0].lmowner_.erase(la[0].lmowner_.begin() + (i - 1) * (dim + 1) + dim);
        la[0].stride_[i - 1] = dim;
      }
    }
    break;
    case FLD::ba_none:
      FOUR_C_THROW("No action supplied");
      break;
    default:
      // call standard fluid boundary element
      FluidBoundary::LocationVector(dis, la, doDirichlet, condstring, params);
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
