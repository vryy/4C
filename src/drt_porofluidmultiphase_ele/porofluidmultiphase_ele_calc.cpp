/*----------------------------------------------------------------------*/
/*!
 \brief implementation of the evaluation routines of the porofluidmultiphase element

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "porofluidmultiphase_ele_calc.H"
#include "porofluid_phasemanager.H"
#include "porofluid_variablemanager.H"
#include "porofluid_evaluator.H"
#include "porofluidmultiphase_ele_parameter.H"

#include "../drt_mat/fluidporo_multiphase.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::PoroFluidMultiPhaseEleCalc(
    const int numdofpernode, const std::string& disname)
    : ele_(NULL),
      totalnumdofpernode_(numdofpernode),
      numfluidphases_(0),
      para_(DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(
          disname)),  // standard parameter list
      xsi_(true),
      xyze0_(true),
      funct_(true),
      deriv_(true),
      deriv2_(true),
      derxy_(true),
      derxy2_(true),
      xjm_(true),
      xij_(true),
      det_(0.0),
      J_(0.0),
      phasemanager_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 | singleton access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>*
DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Instance(const int numdofpernode,
    const std::string& disname, const PoroFluidMultiPhaseEleCalc* delete_me)
{
  static std::map<std::pair<std::string, int>, PoroFluidMultiPhaseEleCalc<distype>*> instances;

  std::pair<std::string, int> key(disname, numdofpernode);

  if (delete_me == NULL)
  {
    if (instances.find(key) == instances.end())
      instances[key] = new PoroFluidMultiPhaseEleCalc<distype>(numdofpernode, disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::pair<std::string, int>,
             PoroFluidMultiPhaseEleCalc<distype>*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}

/*----------------------------------------------------------------------*
 | singleton destruction                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Done()
{
  // delete instance
  Instance(0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate  routine                                         vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, std::vector<Epetra_SerialDenseMatrix*>& elemat,
    std::vector<Epetra_SerialDenseVector*>& elevec)
{
  // check for the action parameter
  const POROFLUIDMULTIPHASE::Action action =
      DRT::INPUT::get<POROFLUIDMULTIPHASE::Action>(params, "action");

  // setup
  if (SetupCalc(ele, discretization, action) == -1) return -1;

  // extract element based or nodal values
  ExtractElementAndNodeValues(ele, params, discretization, la);

  // evaluate action
  EvaluateAction(ele, params, discretization, action, la, elemat, elevec);

  return 0;
}

/*----------------------------------------------------------------------*
 | evaluate action                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const POROFLUIDMULTIPHASE::Action& action, DRT::Element::LocationArray& la,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec)
{
  // determine and evaluate action
  switch (action)
  {
    case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
    case POROFLUIDMULTIPHASE::calc_initial_time_deriv:
    case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
    case POROFLUIDMULTIPHASE::calc_domain_integrals:
    {
      // loop over gauss points and evaluate terms (standard call)
      GaussPointLoop(ele, elemat, elevec, discretization, la);
      break;
    }
    case POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat:
    {
      // loop over gauss points and evaluate off-diagonal terms
      GaussPointLoopODStruct(ele, elemat, elevec, discretization, la);
      break;
    }
    case POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat:
    {
      // loop over gauss points and evaluate off-diagonal terms
      GaussPointLoopODScatra(ele, elemat, elevec, discretization, la);
      break;
    }
    case POROFLUIDMULTIPHASE::calc_pres_and_sat:
    {
      // loop over nodes and evaluate terms
      NodeLoop(ele, elemat, elevec, discretization, la,
          false  // no Jacobian at the node needed
      );
      break;
    }
    case POROFLUIDMULTIPHASE::calc_porosity:
    case POROFLUIDMULTIPHASE::calc_solidpressure:
    {
      // loop over nodes and evaluate terms
      NodeLoop(ele, elemat, elevec, discretization, la,
          true  // Jacobian at the node needed since porosity and (pressure in case of volume
                // fractions )depends on J
      );
      break;
    }
    case POROFLUIDMULTIPHASE::calc_valid_dofs:
    {
      // evaluate the element
      EvaluateOnlyElement(ele, elemat, elevec, discretization, la);
      break;
    }
    default:
    {
      dserror("Not acting on this action. Forgot implementation?");
      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoop(DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // prepare gauss point evaluation
  PrepareGaussPointLoop(ele);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      POROFLUIDMULTIPHASE::ELEUTILS::DisTypeToOptGaussRule<distype>::rule);

  // start loop over gauss points
  GaussPointLoop(intpoints, ele, elemat, elevec, discretization, la);

  return;
}

/*----------------------------------------------------------------------*
| calculate off-diagonal fluid-struct-coupling matrix  kremheller 03/17 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoopODStruct(DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // prepare gauss point evaluation
  PrepareGaussPointLoop(ele);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      POROFLUIDMULTIPHASE::ELEUTILS::DisTypeToOptGaussRule<distype>::rule);

  // start loop over gauss points
  GaussPointLoopODStruct(intpoints, ele, elemat, elevec, discretization, la);

  return;
}

/*----------------------------------------------------------------------*
| calculate off-diagonal fluid-scatra-coupling matrix  kremheller 06/17 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoopODScatra(DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // prepare gauss point evaluation
  PrepareGaussPointLoop(ele);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(
      POROFLUIDMULTIPHASE::ELEUTILS::DisTypeToOptGaussRule<distype>::rule);

  // start loop over gauss points
  GaussPointLoopODScatra(intpoints, ele, elemat, elevec, discretization, la);

  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::PrepareGaussPointLoop(DRT::Element* ele)
{
  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoop(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // start the loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    variablemanager_->EvaluateGPVariables(funct_, derxy_);

    // set the current gauss point state in the manager
    phasemanager_->EvaluateGPState(J_, *variablemanager_);

    // stabilization parameter and integration factors
    // const double taufac     = tau[k]*fac;
    const double timefacfac = para_->TimeFac() * fac;
    // const double timetaufac = para_->TimeFac()*taufac;

    double rhsfac = para_->TimeFacRhs() * fac;

    //----------------------------------------------------------------
    // 1) element matrix
    //----------------------------------------------------------------
    evaluator_->EvaluateMatrix(elemat, funct_, derxy_, totalnumdofpernode_, *phasemanager_,
        *variablemanager_, timefacfac, fac);

    //----------------------------------------------------------------
    // 2) element right hand side
    //----------------------------------------------------------------

    evaluator_->EvaluateVector(elevec, funct_, derxy_, totalnumdofpernode_, *phasemanager_,
        *variablemanager_, rhsfac, fac);

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
  }

  return;
}

/*----------------------------------------------------------------------*
| calculate off-diagonal fluid-struct-coupling matrix  kremheller 03/17 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoopODStruct(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // start the loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    variablemanager_->EvaluateGPVariables(funct_, derxy_);

    // set the current gauss point state in the manager
    phasemanager_->EvaluateGPState(J_, *variablemanager_);

    // stabilization parameter and integration factors
    // const double taufac     = tau[k]*fac;
    double rhsfac = para_->TimeFacRhs() * fac;
    // const double timetaufac = para_->TimeFac()*taufac;

    //----------------------------------------------------------------
    // 1) element off-diagonal matrix
    //----------------------------------------------------------------
    evaluator_->EvaluateMatrixODStruct(elemat, funct_, deriv_, derxy_, xjm_, totalnumdofpernode_,
        *phasemanager_, *variablemanager_, rhsfac, fac, det_);

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
  }

  return;
}

/*----------------------------------------------------------------------*
| calculate off-diagonal fluid-scatra-coupling matrix  kremheller 03/17 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoopODScatra(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // start the loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    variablemanager_->EvaluateGPVariables(funct_, derxy_);

    // set the current gauss point state in the manager
    phasemanager_->EvaluateGPState(J_, *variablemanager_);

    // stabilization parameter and integration factors
    // const double taufac     = tau[k]*fac;
    double rhsfac = para_->TimeFacRhs() * fac;

    // const double timetaufac = para_->TimeFac()*taufac;

    //----------------------------------------------------------------
    // 1) element off-diagonal matrix
    //----------------------------------------------------------------
    evaluator_->EvaluateMatrixODScatra(elemat, funct_, derxy_, totalnumdofpernode_, *phasemanager_,
        *variablemanager_, rhsfac, fac);

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | loop over nodes for evaluation                                  vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::NodeLoop(DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    const bool jacobian_needed)
{
  // Note: this is an evaluation at every node (e.g. for post processing)
  //       Hence, the shape function matrix funct and the derivative matrix derxy
  //       are set in such a way that the standard routines can be used

  // dummy derivative matrix as zero
  derxy_.Clear();

  for (int inode = 0; inode < nen_; ++inode)
  {
    // dummy shape function matrix (set to zero)
    funct_.Clear();
    // set value at current node to 1
    funct_(inode) = 1.0;

    if (jacobian_needed)
    {
      ComputeJacobianAtNode(inode);
    }
    else
      J_ = 1.0;

    // evaluate needed variables
    variablemanager_->EvaluateGPVariables(funct_, derxy_);

    // set the current node state as GP state in the manager
    phasemanager_->EvaluateGPState(J_, *variablemanager_);

    //----------------------------------------------------------------
    // 1) element matrix
    //----------------------------------------------------------------
    evaluator_->EvaluateMatrix(
        elemat, funct_, derxy_, totalnumdofpernode_, *phasemanager_, *variablemanager_, 1.0, 1.0);

    //----------------------------------------------------------------
    // 2) element vector
    //----------------------------------------------------------------
    evaluator_->EvaluateVector(
        elevec, funct_, derxy_, totalnumdofpernode_, *phasemanager_, *variablemanager_, 1.0, 1.0);

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
  }
}

/*-----------------------------------------------------------------------------*
 | loop over nodes for evaluation                                  vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvaluateOnlyElement(DRT::Element* ele,
    std::vector<Epetra_SerialDenseMatrix*>& elemat, std::vector<Epetra_SerialDenseVector*>& elevec,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  // Note: this is an evaluation of the element where no Gauss point data and function values (from
  //       variable and phase-managers are required

  // dummy derivative matrix as zero
  derxy_.Clear();
  funct_.Clear();

  //----------------------------------------------------------------
  // 1) element matrix
  //----------------------------------------------------------------
  evaluator_->EvaluateMatrix(
      elemat, funct_, derxy_, totalnumdofpernode_, *phasemanager_, *variablemanager_, 1.0, 1.0);

  //----------------------------------------------------------------
  // 2) element vector
  //----------------------------------------------------------------
  evaluator_->EvaluateVector(
      elevec, funct_, derxy_, totalnumdofpernode_, *phasemanager_, *variablemanager_, 1.0, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | setup element evaluation                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::SetupCalc(DRT::Element* ele,
    DRT::Discretization& discretization, const POROFLUIDMULTIPHASE::Action& action)
{
  // get element coordinates
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, xyze0_);

  // set current coordinates to initial coordinates
  // the displacements will be added later in ExtractElementAndNodeValues() for the moving mesh case
  xyze_ = xyze0_;

  // set element
  ele_ = ele;

  Teuchos::RCP<MAT::Material> mat = ele->Material();
  if (mat->MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
      mat->MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("PoroFluidMultiPhase element got unsupported material type %d", mat->MaterialType());

  Teuchos::RCP<MAT::FluidPoroMultiPhase> actmat =
      Teuchos::rcp_static_cast<MAT::FluidPoroMultiPhase>(mat);
  if (actmat == Teuchos::null) dserror("cast failed");
  numfluidphases_ = actmat->NumFluidPhases();

  // Note:
  // here the phase manager, the variable manager and the evaluator classes are
  // rebuild here, allocating new memory. This is actually not necessary.
  // If this proves to be a performance issue, then the Create*** methods
  // could be adapted to have static members, which are only built once.

  // build the phase manager
  phasemanager_ = DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
      *para_, nsd_, ele->Material()->MaterialType(), action, totalnumdofpernode_, numfluidphases_);
  // setup the manager
  phasemanager_->Setup(ele);

  // rebuild the phase manager
  variablemanager_ =
      DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<nsd_, nen_>::CreateVariableManager(
          *para_, action, actmat, totalnumdofpernode_, numfluidphases_);

  // build the evaluator
  evaluator_ = DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<nsd_, nen_>::CreateEvaluator(
      *para_, action, totalnumdofpernode_, numfluidphases_, *phasemanager_);

  return 0;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values                     vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  variablemanager_->ExtractElementAndNodeValues(*ele, discretization, la, xyze_);
  return;
}


/*------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     vuong 08/16 |
 *------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int iquad                                          ///< id of current Gauss point
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim = 0; idim < nsd_; idim++) xsi_(idim) = gpcoord[idim];

  det_ = EvalShapeFuncAndDerivsInParameterSpace();

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", ele_->Id(), det_);

  // compute global spatial derivatives
  derxy_.Multiply(xij_, deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderiv_)
  {
    // get global second derivatives
    DRT::UTILS::gder2<distype, nen_>(xjm_, derxy_, deriv2_, xyze_, derxy2_);
  }
  else
    derxy2_.Clear();

  //------------------------get determinant of Jacobian dX / ds
  // transposed jacobian "dX/ds"
  LINALG::Matrix<nsd_, nsd_> xjm0;
  xjm0.MultiplyNT(deriv_, xyze0_);

  // inverse of transposed jacobian "ds/dX"
  const double det0 = xjm0.Determinant();

  // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds)
  // )^-1
  J_ = det_ / det0;

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad] * det_;

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

}  // PoroFluidMultiPhaseEleCalc::EvalShapeFuncAndDerivsAtIntPoint

/*--------------------------------------------------------------------------*
 | evaluate shape functions and derivatives in parameter space   vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvalShapeFuncAndDerivsInParameterSpace()
{
  double det = 0.0;

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_, funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);
  if (use2ndderiv_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_, deriv2_);
  }

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */

  xjm_.MultiplyNT(deriv_, xyze_);
  det = xij_.Invert(xjm_);

  return det;
}

/*--------------------------------------------------------------------------*
 | Compute Jacobian at node 'inode'                        kremheller 04/17 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::ComputeJacobianAtNode(const int inode)
{
  // get parameter space coordinates of current node
  LINALG::Matrix<3, 1> myXi = DRT::UTILS::getNodeCoordinates(inode, distype);
  for (int idim = 0; idim < nsd_; idim++) xsi_(idim) = myXi(idim);

  det_ = EvalShapeFuncAndDerivsInParameterSpace();

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", ele_->Id(), det_);

  // compute global spatial derivatives
  derxy_.Multiply(xij_, deriv_);

  //------------------------get determinant of Jacobian dX / ds
  // transposed jacobian "dX/ds"
  LINALG::Matrix<nsd_, nsd_> xjm0;
  xjm0.MultiplyNT(deriv_, xyze0_);

  // inverse of transposed jacobian "ds/dX"
  const double det0 = xjm0.Determinant();

  // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds)
  // )^-1
  J_ = det_ / det0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2,2>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2,3>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3,2>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tri3>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad4>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad9>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex8>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tet10>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::nurbs27>;
