/*!----------------------------------------------------------------------
\file scatra_element_boundary_evaluate.cpp
\brief

Evaluate boundary conditions for scalar transport problems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *----------------------------------------------------------------------*/
#include "../drt_mat/elchmat.H"

#include "scatra_ele_action.H"
#include "scatra_ele_boundary_calc.H"
#include "scatra_ele_boundary_factory.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  if(scatratype == INPAR::SCATRA::scatratype_elch)
  {
    numscal -= 1;

    // get the material of the first element
    // we assume here, that the material is equal for all elements in this discretization
    // get the parent element including its material
    Teuchos::RCP<MAT::Material> material = ParentElement()->Material();
    if (material->MaterialType() == INPAR::MAT::m_elchmat)
      numscal = static_cast<const MAT::ElchMat*>(material.get())->NumScal();
  }

  // determine implementation type
  INPAR::SCATRA::ImplType impltype = INPAR::SCATRA::impltype_undefined;
  switch(scatratype)
  {
  case INPAR::SCATRA::scatratype_condif:             impltype = INPAR::SCATRA::impltype_std;                break;
  case INPAR::SCATRA::scatratype_advreac:            impltype = INPAR::SCATRA::impltype_advreac;            break;
  case INPAR::SCATRA::scatratype_anisotrop:          impltype = INPAR::SCATRA::impltype_aniso;              break;
  case INPAR::SCATRA::scatratype_cardiac_monodomain: impltype = INPAR::SCATRA::impltype_cardiac_monodomain; break;
  case INPAR::SCATRA::scatratype_levelset:           impltype = INPAR::SCATRA::impltype_levelset;           break;
  case INPAR::SCATRA::scatratype_loma:               impltype = INPAR::SCATRA::impltype_loma;               break;
  case INPAR::SCATRA::scatratype_elch:
  {
    // At this point, an instance of the singleton parameter class ScaTraEleParameterElch already exists
    DRT::ELEMENTS::ScaTraEleParameterElch* elchpara = DRT::ELEMENTS::ScaTraEleParameterElch::Instance();

    if(ParentElement()->Material()->MaterialType() == INPAR::MAT::m_electrode)
      impltype = INPAR::SCATRA::impltype_elch_electrode;
    else if(elchpara->ElchType() == INPAR::ELCH::elchtype_diffcond)
      impltype = INPAR::SCATRA::impltype_elch_diffcond;
    else if(elchpara->ElchType() == INPAR::ELCH::elchtype_nernst_planck)
      impltype = INPAR::SCATRA::impltype_elch_NP;
    // else: impltype just remains undefined
    break;
  }
  case INPAR::SCATRA::scatratype_poro:     impltype = INPAR::SCATRA::impltype_poro;     break;
  case INPAR::SCATRA::scatratype_pororeac: impltype = INPAR::SCATRA::impltype_pororeac; break;
  default:
  {
    dserror("Unknown scatratype for boundary element evaluation!");
    break;
  }
  }

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(this,impltype,numdofpernode,numscal)->EvaluateAction(
      this,
      params,
      discretization,
      lm,
      elemat1,
      elemat2,
      elevec1,
      elevec2,
      elevec3
      );
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition on boundary element   fang 01/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  // add Neumann boundary condition to parameter list
  params.set<DRT::Condition*>("condition",&condition);

  // evaluate boundary element
  return Evaluate(params,discretization,lm,*elemat1,*elemat1,elevec1,elevec1,elevec1);
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::LocationVector(
    const Discretization&   dis,
    LocationArray&          la,
    bool                    doDirichlet,
    const std::string&      condstring,
    Teuchos::ParameterList& params
    ) const
{
  // check for the action parameter
  const SCATRA::BoundaryAction action = DRT::INPUT::get<SCATRA::BoundaryAction>(params,"action");
  switch (action)
  {
  case SCATRA::bd_calc_weak_Dirichlet:
    // special cases: the boundary element assembles also into
    // the inner dofs of its parent element
    // note: using these actions, the element will get the parent location vector
    //       as input in the respective evaluate routines
    ParentElement()->LocationVector(dis,la,doDirichlet);
    break;
  default:
    DRT::Element::LocationVector(dis,la,doDirichlet);
    break;
  }
  return;
}
