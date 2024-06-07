/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluate boundary conditions for scalar transport problems

\level 2

 *----------------------------------------------------------------------*/
#include "4C_fem_discretization.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_boundary_calc.hpp"
#include "4C_scatra_ele_boundary_factory.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportBoundary::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  FOUR_C_THROW("not implemented. Use the Evaluate() method with Location Array instead!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportBoundary::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  // perform additional operations specific to implementation type
  switch (parent_element()->ImplType())
  {
    case Inpar::ScaTra::impltype_elch_diffcond:
    case Inpar::ScaTra::impltype_elch_diffcond_thermo:
    case Inpar::ScaTra::impltype_elch_electrode:
    case Inpar::ScaTra::impltype_elch_electrode_growth:
    case Inpar::ScaTra::impltype_elch_electrode_thermo:
    case Inpar::ScaTra::impltype_elch_NP:
    {
      // adapt number of transported scalars for electrochemistry problems
      numscal -= 1;

      // get the material of the first element
      // we assume here, that the material is equal for all elements in this discretization
      // get the parent element including its material
      Teuchos::RCP<Core::Mat::Material> material = parent_element()->Material();
      if (material->MaterialType() == Core::Materials::m_elchmat)
        numscal = static_cast<const Mat::ElchMat*>(material.get())->NumScal();

      break;
    }

    case Inpar::ScaTra::impltype_std:
    case Inpar::ScaTra::impltype_advreac:
    case Inpar::ScaTra::impltype_refconcreac:
    case Inpar::ScaTra::impltype_chemo:
    case Inpar::ScaTra::impltype_chemoreac:
    case Inpar::ScaTra::impltype_aniso:
    case Inpar::ScaTra::impltype_cardiac_monodomain:
    case Inpar::ScaTra::impltype_levelset:
    case Inpar::ScaTra::impltype_loma:
    case Inpar::ScaTra::impltype_poro:
    case Inpar::ScaTra::impltype_pororeac:
    case Inpar::ScaTra::impltype_thermo_elch_diffcond:
    case Inpar::ScaTra::impltype_thermo_elch_electrode:
    case Inpar::ScaTra::impltype_multipororeac:
      // do nothing in these cases
      break;

    default:
    {
      // other implementation types are invalid
      FOUR_C_THROW("Invalid implementation type!");
      break;
    }
  }

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return Discret::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(
      this, parent_element()->ImplType(), numdofpernode, numscal, discretization.Name())
      ->Evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition on boundary element   fang 01/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportBoundary::evaluate_neumann(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // add Neumann boundary condition to parameter list
  params.set<Core::Conditions::Condition*>("condition", &condition);

  LocationArray la(discretization.NumDofSets());
  Core::Elements::Element::LocationVector(discretization, la, false);

  // evaluate boundary element
  return Evaluate(params, discretization, la, *elemat1, *elemat1, elevec1, elevec1, elevec1);
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::TransportBoundary::LocationVector(const Discretization& dis,
    LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<ScaTra::BoundaryAction>(params, "action");
  switch (action)
  {
    case ScaTra::BoundaryAction::calc_weak_Dirichlet:
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      parent_element()->LocationVector(dis, la, doDirichlet);
      break;
    default:
      Core::Elements::Element::LocationVector(dis, la, doDirichlet);
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
