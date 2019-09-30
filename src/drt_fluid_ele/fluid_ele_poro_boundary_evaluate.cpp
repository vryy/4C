/*-----------------------------------------------------------*/
/*! \file

\brief evaluation routines for the fluid poro boundary element

\maintainer Johannes Kremheller

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_ele_poro.H"
#include "fluid_ele_action.H"
#include "fluid_ele_boundary_calc.H"
#include "fluid_ele_boundary_factory.H"

#include "../drt_inpar/inpar_fluid.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            vuong 11/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidPoroBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params, "action");

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
      dserror("invalid physical type for porous fluid!");
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
  }  // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element      (public) vuong 11/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoroBoundary::LocationVector(const Discretization& dis, LocationArray& la,
    bool doDirichlet, const std::string& condstring, Teuchos::ParameterList& params) const
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params, "action");
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
      // Todo: this is a hack ...
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
      dserror("No action supplied");
      break;
    default:
      // call standard fluid boundary element
      FluidBoundary::LocationVector(dis, la, doDirichlet, condstring, params);
      break;
  }
  return;
}
