/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_integrator.cpp

\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for
       augmented contact)

\level 2

\maintainer Matthias Mayr

\date Apr 28, 2014

*/
/*---------------------------------------------------------------------*/
#include "contact_augmented_integrator.H"
#include "contact_integrator_utils.H"
#include "contact_aug_element_utils.H"

#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_serialdensevector.H"

#include <Epetra_Map.h>
#include <Teuchos_TimeMonitor.hpp>

// define and initialize static member
CONTACT::INTEGRATOR::UniqueProjInfoPair CONTACT::AUG::IntegrationWrapper::projInfo_(0);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegrationWrapper::IntegrationWrapper(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm)
    : CONTACT::CoIntegrator::CoIntegrator(params, eletype, comm), integrator_(NULL)

{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlane called without specific shape "
        "function defined!");

  // check for problem dimension
  dsassert(Dim() == 3, "ERROR: 3D integration method called for non-3D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of "
        "MortarElement pair!");
  if (cell == Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  if (ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlane supports no Dual shape functions for the "
        "augmented Lagrange solving strategy!");

  GlobalTimeMonitor* timer_ptr = cparams_ptr->GetTimer<GlobalTimeID>(0);

  timer_ptr->start(GlobalTimeID::IntegrateDerivCell3DAuxPlane);
  integrator_ = IntegratorGeneric::Create(Dim(), sele.Shape(), mele.Shape(), *cparams_ptr, this);
  integrator_->IntegrateDerivCell3DAuxPlane(sele, mele, *cell, auxn);
  timer_ptr->stop(GlobalTimeID::IntegrateDerivCell3DAuxPlane);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivEle3D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror(
        "ERROR: IntegrateDerivCell3DAuxPlane called without specific shape "
        "function defined!");

  // check for problem dimension
  dsassert(Dim() == 3, "ERROR: 3D integration method called for non-3D problem");

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

  // check input data
  for (unsigned test = 0; test < meles.size(); ++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror(
          "ERROR: IntegrateDerivCell3D called on a wrong type of "
          "MortarElement pair!");
  }

  // contact with wear
  if (wearlaw_ != INPAR::WEAR::wear_none) dserror("Wear is not supported!");

  // Boundary Segmentation check -- HasProj()-check
  //  *boundary_ele = BoundarySegmCheck3D(sele,meles);
  *boundary_ele = false;

  GlobalTimeMonitor* timer_ptr = cparams_ptr->GetTimer<GlobalTimeID>(0);
  timer_ptr->start(GlobalTimeID::IntegrateDerivEle3D);

  *proj = INTEGRATOR::FindFeasibleMasterElements(sele, meles, boundary_ele, *this, projInfo_);

  for (auto& info_pair : projInfo_)
  {
    MORTAR::MortarElement& mele = *(info_pair.first);
    integrator_ = IntegratorGeneric::Create(Dim(), sele.Shape(), mele.Shape(), *cparams_ptr, this);
    integrator_->Evaluate(sele, mele, *boundary_ele, info_pair.second);
  }

  timer_ptr->stop(GlobalTimeID::IntegrateDerivEle3D);

  Epetra_Vector* sele_times = cparams_ptr->Get<Epetra_Vector>(0);
  const int slid = sele_times->Map().LID(sele.Id());
  if (slid == -1)
    dserror("Couldn't find the current slave element GID #%d on proc #%d.", sele.Id(),
        sele_times->Map().Comm().MyPID());
  (*sele_times)[slid] += timer_ptr->getLastTimeIncr();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSlaveElement(MORTAR::MortarElement& sele,
    const Epetra_Comm& comm, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr, true);

  IntegrateDerivSlaveElement(sele, comm, cparams_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSlaveElement(MORTAR::MortarElement& sele,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  integrator_ = IntegratorGeneric::Create(Dim(), sele.Shape(), sele.Shape(), *cparams_ptr, this);
  integrator_->IntegrateDerivSlaveElement(sele);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSegment2D(MORTAR::MortarElement& sele,
    double& sxia, double& sxib, MORTAR::MortarElement& mele, double& mxia, double& mxib,
    const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape() == MORTAR::MortarElement::line3 ||
      ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: Petrov-Galerkin / quadratic FE interpolation not yet implemented.");

  // check for problem dimension
  dsassert(Dim() == 2, "ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia < -1.0) || (sxib > 1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia < -1.0) || (mxib > 1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  GlobalTimeMonitor* timer_ptr = cparams_ptr->GetTimer<GlobalTimeID>(0);
  timer_ptr->start(GlobalTimeID::IntegrateDerivSegment2D);

  integrator_ = IntegratorGeneric::Create(Dim(), sele.Shape(), mele.Shape(), *cparams_ptr, this);
  integrator_->IntegrateDerivSegment2D(sele, sxia, sxib, mele, mxia, mxib);

  timer_ptr->stop(GlobalTimeID::IntegrateDerivSegment2D);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivEle2D(MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  // check for problem dimension
  if (Dim() != 2) dserror("ERROR: 2D integration method called for non-2D problem");

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // check input data
  for (int i = 0; i < (int)meles.size(); ++i)
  {
    if ((!sele.IsSlave()) || (meles[i]->IsSlave()))
      dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  }

  // number of nodes (slave) and problem dimension
  const int nrow = sele.NumNode();

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k = 0; k < nrow; ++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  GlobalTimeMonitor* timer_ptr = cparams_ptr->GetTimer<GlobalTimeID>(0);
  timer_ptr->start(GlobalTimeID::IntegrateDerivEle2D);

  // Boundary Segmentation check -- HasProj()-check
  if (IntegrationType() == INPAR::MORTAR::inttype_elements_BS)
    *boundary_ele = BoundarySegmCheck2D(sele, meles);

  if (*boundary_ele == false || IntegrationType() == INPAR::MORTAR::inttype_elements)
  {
    INTEGRATOR::FindFeasibleMasterElements(sele, meles, *this, projInfo_);

    for (auto& info_pair : projInfo_)
    {
      MORTAR::MortarElement& mele = *(info_pair.first);
      integrator_ =
          IntegratorGeneric::Create(Dim(), sele.Shape(), mele.Shape(), *cparams_ptr, this);
      integrator_->Evaluate(sele, mele, false, info_pair.second);
    }
  }  // boundary_ele check

  timer_ptr->stop(GlobalTimeID::IntegrateDerivEle2D);

  Epetra_Vector* sele_times = cparams_ptr->Get<Epetra_Vector>(0);
  const int slid = sele_times->Map().LID(sele.Id());
  if (slid == -1)
    dserror("Couldn't find the current slave element GID #%d on proc #%d.", sele.Id(),
        sele_times->Map().Comm().MyPID());
  (*sele_times)[slid] += timer_ptr->getLastTimeIncr();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create(int probdim,
    DRT::Element::DiscretizationType slavetype, DRT::Element::DiscretizationType mastertype,
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper)
{
  switch (probdim)
  {
    case 2:
      return Create2D(slavetype, mastertype, cparams, wrapper);
    case 3:
      return Create3D(slavetype, mastertype, cparams, wrapper);
    default:
      dserror("Unsupported problem dimension %d", probdim);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create2D(
    DRT::Element::DiscretizationType slavetype, DRT::Element::DiscretizationType mastertype,
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper)
{
  switch (slavetype)
  {
    case DRT::Element::line2:
      return Create2D<DRT::Element::line2>(mastertype, cparams, wrapper);
    case DRT::Element::nurbs2:
      return Create2D<DRT::Element::nurbs2>(mastertype, cparams, wrapper);
    case DRT::Element::nurbs3:
      return Create2D<DRT::Element::nurbs3>(mastertype, cparams, wrapper);
    default:
      dserror("Unsupported slave element type %d|\"%s\"", slavetype,
          DRT::DistypeToString(slavetype).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType slavetype>
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create2D(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper)
{
  switch (mastertype)
  {
    case DRT::Element::line2:
      return Create2D<slavetype, DRT::Element::line2>(cparams, wrapper);
    case DRT::Element::nurbs2:
      return Create2D<slavetype, DRT::Element::nurbs2>(cparams, wrapper);
    case DRT::Element::nurbs3:
      return Create2D<slavetype, DRT::Element::nurbs3>(cparams, wrapper);
    default:
      dserror("Unsupported master element type %d|\"%s\"", mastertype,
          DRT::DistypeToString(mastertype).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType slavetype, DRT::Element::DiscretizationType mastertype>
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create2D(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper)
{
  const enum INPAR::CONTACT::VariationalApproach var_type = cparams.GetVariationalApproachType();

  switch (var_type)
  {
    case INPAR::CONTACT::var_incomplete:
    {
      typedef DebugIncompleteIntPolicy<2, slavetype, mastertype> incomplete_policy;
      return Integrator<2, slavetype, mastertype, incomplete_policy>::Instance(&cparams, wrapper);
    }
    case INPAR::CONTACT::var_complete:
    {
      typedef DebugCompleteIntPolicy<2, slavetype, mastertype> complete_policy;
      return Integrator<2, slavetype, mastertype, complete_policy>::Instance(&cparams, wrapper);
    }
    default:
    {
      dserror("Unknown variational approach! (var_type= \"%s\" | %d)",
          INPAR::CONTACT::VariationalApproach2String(var_type).c_str(), var_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create3D(
    DRT::Element::DiscretizationType slavetype, DRT::Element::DiscretizationType mastertype,
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper)
{
  switch (slavetype)
  {
    case DRT::Element::quad4:
      return Create3D<DRT::Element::quad4>(mastertype, cparams, wrapper);
    case DRT::Element::tri3:
      return Create3D<DRT::Element::tri3>(mastertype, cparams, wrapper);
    case DRT::Element::nurbs4:
      return Create3D<DRT::Element::nurbs4>(mastertype, cparams, wrapper);
    case DRT::Element::nurbs9:
      return Create3D<DRT::Element::nurbs9>(mastertype, cparams, wrapper);
    default:
      dserror("Unsupported slave element type %d|\"%s\"", DRT::DistypeToString(mastertype).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType slavetype>
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create3D(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper)
{
  switch (mastertype)
  {
    case DRT::Element::quad4:
      return Create3D<slavetype, DRT::Element::quad4>(cparams, wrapper);
    case DRT::Element::tri3:
      return Create3D<slavetype, DRT::Element::tri3>(cparams, wrapper);
    case DRT::Element::nurbs4:
      return Create3D<slavetype, DRT::Element::nurbs4>(cparams, wrapper);
    case DRT::Element::nurbs9:
      return Create3D<slavetype, DRT::Element::nurbs9>(cparams, wrapper);
    default:
      dserror(
          "Unsupported master element type %d|\"%s\"", DRT::DistypeToString(mastertype).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType slavetype, DRT::Element::DiscretizationType mastertype>
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create3D(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper)
{
  const enum INPAR::CONTACT::VariationalApproach var_type = cparams.GetVariationalApproachType();

  switch (var_type)
  {
    case INPAR::CONTACT::var_incomplete:
    {
      typedef DebugIncompleteIntPolicy<3, slavetype, mastertype> incomplete_policy;
      return Integrator<3, slavetype, mastertype, incomplete_policy>::Instance(&cparams, wrapper);
    }
    case INPAR::CONTACT::var_complete:
    {
      typedef DebugCompleteIntPolicy<3, slavetype, mastertype> complete_policy;
      return Integrator<3, slavetype, mastertype, complete_policy>::Instance(&cparams, wrapper);
    }
    default:
    {
      dserror("Unknown variational approach! (var_type= \"%s\" | %d)",
          INPAR::CONTACT::VariationalApproach2String(var_type).c_str(), var_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>*
CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::Instance(
    CONTACT::ParamsInterface* cparams, CONTACT::CoIntegrator* wrapper, const bool delete_me)
{
  static Integrator<probdim, slavetype, mastertype, IntPolicy>* instance = NULL;

  if (delete_me)
  {
    if (instance)
    {
      instance->IntPolicy::timer_.write(std::cout);
      delete instance;
    }

    return NULL;
  }

  if (not instance)
  {
    instance = new Integrator<probdim, slavetype, mastertype, IntPolicy>();
  }

  instance->Init(cparams, wrapper);
  instance->IntPolicy::timer_.setComm(&wrapper->Comm());

  return instance;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::Done()
{
  Instance(NULL, NULL, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::Integrator()
    : IntegratorGeneric(), IntPolicy()
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::Integrator(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator& wrapper)
    : IntegratorGeneric(cparams, wrapper), IntPolicy()
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::IntegrateDerivSegment2D(
    MORTAR::MortarElement& sele, double& sxia, double& sxib, MORTAR::MortarElement& mele,
    double& mxia, double& mxib)
{
  dserror(
      "Deprecated method! The segmented based integration is no longer "
      "supported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype,
    IntPolicy>::IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, MORTAR::IntCell& cell, double* auxn)
{
  dserror(
      "Deprecated method! The segmented based integration is no longer "
      "supported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype,
    IntPolicy>::IntegrateDerivSlaveElement(MORTAR::MortarElement& sele)
{
  // set evaluator
  const enum MORTAR::ActionType action = CParams().GetActionType();
  SetEvaluator(action);

  for (int gp = 0; gp < this->Wrapper().nGP(); ++gp)
  {
    const double eta[2] = {this->Wrapper().Coordinate(gp, 0), this->Wrapper().Coordinate(gp, 1)};
    const double wgt = this->Wrapper().Weight(gp);

    // get Gauss point in slave element coordinates
    const double sxi[2] = {eta[0], eta[1]};
    const LINALG::Matrix<2, 1> sxi_mat(sxi, true);

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(this->ShapeFcn(), sxi, lmval_, lmderiv_, my::SLAVENUMNODE, true);

    // evaluate shape function and derivative values (on slave element)
    shape_function_and_deriv1<slavetype>(sele, sxi_mat, sval_, sderiv_);

    // integrate the slave jacobian
    const double jac = sele.Jacobian(sxi);

    // evaluate the convective slave base vectors
    LINALG::Matrix<3, 2> stau;
    sele.Metrics(sxi, &stau(0, 0), &stau(0, 1));

    // evaluate the slave Jacobian 1-st order derivative
    evaluator_->Deriv_Jacobian(sele, sxi, sderiv_, stau);

    // *** SLAVE NODES ****************************************************
    // compute the tributary area
    GP_AugA(sele, lmval_, wgt, jac);

    // compute 1-st order derivative of the tributary area
    Get_Deriv1st_AugA(sele, lmval_, wgt, jac, derivjac_);

    // compute 2-nd order derivative of the tributary area
    evaluator_->Get_Deriv2nd_AugA(sele, lmval_, wgt, deriv2ndjac_);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::Evaluate(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, bool boundary_ele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo)
{
  if (this->Wrapper().IntegrationType() != INPAR::MORTAR::inttype_elements)
    dserror("How did you come here?");

  const enum MORTAR::ActionType action = CParams().GetActionType();

  // set the evaluator: 1-st derivatives only, or 1-st AND 2-nd derivatives
  SetEvaluator(action);

  // choose the integration scheme
  switch (action)
  {
    case MORTAR::eval_static_constraint_rhs:
    {
      IntegrateWeightedGap(sele, mele, boundary_ele, projInfo);
      break;
    }
    case MORTAR::eval_force_stiff:
    case MORTAR::eval_force:
    {
      IntegrateDerivEle(sele, mele, boundary_ele, projInfo);
      break;
    }
    case MORTAR::eval_wgap_gradient_error:
    {
      IntegrateWeightedGapGradientError(sele, mele, boundary_ele, projInfo);
      break;
    }
    default:
    {
      dserror("Unconsidered ActionType = %d | \"%s\" ", action,
          MORTAR::ActionType2String(action).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::SetEvaluator(
    const enum MORTAR::ActionType action)
{
  switch (action)
  {
    case MORTAR::eval_static_constraint_rhs:
    {
      /* do nothing, since no derivatives have to be evaluated */
      break;
    }
    case MORTAR::eval_force:
    case MORTAR::eval_wgap_gradient_error:
    {
      if (evaluator_.is_null() or evaluator_->GetType() != Evaluator::Type::deriv1st_only)
        evaluator_ = Teuchos::rcp(new EvaluatorDeriv1stOnly(*this));

      //      static int count = 0;
      //      std::cout << "eval_force = " << ++count << std::endl;
      break;
    }
    case MORTAR::eval_force_stiff:
    {
      if (evaluator_.is_null() or evaluator_->GetType() != Evaluator::Type::full)
        evaluator_ = Teuchos::rcp(new EvaluatorFull(*this));

      //      static int count = 0;
      //      std::cout << "eval_force_stiff = " << ++count << std::endl;
      break;
    }
    default:
    {
      dserror("Unconsidered ActionType = %d | \"%s\" ", action,
          MORTAR::ActionType2String(action).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::IntegrateDerivEle(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, bool boundary_ele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo)
{
  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);

  const int linsize = GetLinSize(sele);

  {
    // get the gausspoints of this slave / master element pair
    const unsigned num_gps = projInfo.gaussPoints_.size();

    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    HardReset(linsize);

    for (my::gp_id_ = 0; static_cast<unsigned>(my::gp_id_) < num_gps; ++my::gp_id_)
    {
      const int gp = projInfo.gaussPoints_[my::gp_id_];

      // coordinates and weight
      const double eta[2] = {this->Wrapper().Coordinate(gp, 0), this->Wrapper().Coordinate(gp, 1)};
      const double wgt = this->Wrapper().Weight(gp) * projInfo.scaling_[my::gp_id_];

      // get Gauss point in slave element coordinates
      const double sxi[2] = {eta[0], eta[1]};
      const LINALG::Matrix<2, 1> sxi_mat(sxi, true);

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(this->ShapeFcn(), sxi, lmval_, lmderiv_, my::SLAVENUMNODE, true);

      // evaluate trace space shape functions (on both elements)
      shape_function_and_deriv1<slavetype>(sele, sxi_mat, sval_, sderiv_);

      // evaluate the convective slave base vectors
      LINALG::Matrix<3, 2> stau;
      sele.Metrics(sxi, &stau(0, 0), &stau(0, 1));

      // evaluate the two Jacobians (int. cell and slave element)
      const double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian 1-st and 2-nd order derivatives
      //      Deriv_Jacobian( sele, sxi, sderiv_, stau, derivjac_, &deriv2ndjac_ );
      evaluator_->Deriv_Jacobian(sele, sxi, sderiv_, stau);

      const double uniqueProjalpha = projInfo.uniqueProjAlpha_[my::gp_id_];
      const LINALG::Matrix<2, 1>& uniqueMxi = projInfo.uniqueMxi_[my::gp_id_];

      mele.GetNodalCoords(mcoord_);

      // get mval
      shape_function_and_deriv1_and_deriv2<mastertype>(mele, uniqueMxi, mval_, mderiv_, mderiv2nd_);

      // evaluate the convective master base vectors
      LINALG::Matrix<3, 2> mtau;
      mele.Metrics(uniqueMxi.A(), &mtau(0, 0), &mtau(0, 1));

      // evaluate the GP master coordinate 1-st and 2-nd order derivatives
      evaluator_->Deriv_MXiGP(
          sele, mele, &sxi[0], uniqueMxi.A(), uniqueProjalpha, sval_, mval_, mderiv_, mtau);

      //**********************************************************************
      // evaluate at GP and lin char. quantities
      //**********************************************************************
      // calculate the averaged normal + derivative at gp level
      GP_Normal_DerivNormal(
          sele, sval_, &gpn_[0], dn_non_unit_, ddn_non_unit_, dn_unit_, ddn_unit_);

      // integrate scaling factor kappa
      GP_kappa(sele, lmval_, wgt, jacslave);

      // integrate the inner integral relating to the first order derivative of
      // the discrete normal gap for later usage (for all found slave nodes)
      IntPolicy::Get_Deriv1st_GapN(
          sele, mele, sval_, mval_, &gpn_[0], mtau, dmxigp_, deriv_gapn_sl_, deriv_gapn_ma_);

      // evaluate normal gap (split into slave and master contributions)
      double gapn_sl = 0.0;
      double gapn_ma = 0.0;
      GapN(sele, mele, sval_, mval_, gpn_, gapn_sl, gapn_ma);

      // evaluate the weighted gap (slave / master)
      GP_WGap(sele, lmval_, gapn_sl, gapn_ma, wgt, jacslave);

      // 1-st order derivative of the weighted gap (variation)
      IntPolicy::Get_Deriv1st_WGap(
          sele, lmval_, gapn_sl, gapn_ma, wgt, jacslave, derivjac_, deriv_gapn_sl_, deriv_gapn_ma_);

      // 1-st order derivative of the weighted gap (necessary for the
      // linearization of the constraint equations in case of the complete AND
      // incomplete variational approach)
      IntPolicy::Get_Deriv1st_WGap_Complete(linsize, sele, mele, sval_, mval_, lmval_, &gpn_[0],
          mtau, dmxigp_, gapn_sl, gapn_ma, wgt, jacslave, derivjac_);

      IntPolicy::Get_Debug(sele, lmval_, gapn_sl, gapn_ma, wgt, jacslave, gpn_, uniqueMxi.A());

      IntPolicy::Get_Deriv1st_Debug(sele, lmval_, sval_, sderiv_, stau, derivjac_, dmxigp_,
          dn_unit_, deriv_gapn_sl_, gapn_sl, wgt, jacslave);

      switch (CParams().GetActionType())
      {
        case MORTAR::eval_force_stiff:
        {
          Get_Deriv1st_Kappa(sele, lmval_, wgt, derivjac_);

          Get_Deriv2nd_Kappa(sele, lmval_, wgt, deriv2ndjac_);

          IntPolicy::Get_Deriv2nd_WGap(sele, mele, sval_, mval_, lmval_, mderiv_, mderiv2nd_, mtau,
              gpn_, wgt, gapn_sl, gapn_ma, jacslave, derivjac_, deriv2ndjac_, dmxigp_, ddmxigp_,
              dn_unit_, ddn_unit_, deriv_gapn_sl_, deriv_gapn_ma_);

          IntPolicy::Get_Deriv2nd_Debug(sele, lmval_, sval_, sderiv_, stau, derivjac_,
              deriv_gapn_sl_, deriv2ndjac_, ddmxigp_, dn_unit_, ddn_unit_, gapn_sl, wgt, jacslave);

          break;
        }
        default:
          // do nothing
          break;
      }

      WeakReset(linsize);
    }  // GP-loop

    IntPolicy::CompleteNodeData(sele);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::IntegrateWeightedGap(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele, bool boundary_ele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo)
{
  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);

  const int linsize = GetLinSize(sele);

  {
    // get the gausspoints of this slave / master element pair
    const unsigned num_gps = projInfo.gaussPoints_.size();

    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    HardReset(linsize);

    for (my::gp_id_ = 0; static_cast<unsigned>(my::gp_id_) < num_gps; ++my::gp_id_)
    {
      const int gp = projInfo.gaussPoints_[my::gp_id_];

      // coordinates and weight
      const double eta[2] = {this->Wrapper().Coordinate(gp, 0), this->Wrapper().Coordinate(gp, 1)};
      const double wgt = this->Wrapper().Weight(gp) * projInfo.scaling_[my::gp_id_];

      // get Gauss point in slave element coordinates
      const double sxi[2] = {eta[0], eta[1]};
      const LINALG::Matrix<2, 1> sxi_mat(sxi, true);

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(this->ShapeFcn(), sxi, lmval_, lmderiv_, my::SLAVENUMNODE, true);

      // evaluate trace space shape functions (on both elements)
      shape_function_and_deriv1<slavetype>(sele, sxi_mat, sval_, sderiv_);

      // evaluate the two Jacobians (int. cell and slave element)
      const double jacslave = sele.Jacobian(sxi);

      const LINALG::Matrix<2, 1>& uniqueMxi = projInfo.uniqueMxi_[my::gp_id_];

      mele.GetNodalCoords(mcoord_);

      // get mval and mderiv1
      shape_function_and_deriv1<mastertype>(mele, uniqueMxi, mval_, mderiv_);

      // integrate scaling factor kappa
      GP_kappa(sele, lmval_, wgt, jacslave);

      // calculate the averaged unified GP normal
      GP_Normal(sele, sval_, &gpn_[0]);

      // evaluate normal gap (split into slave and master contributions)
      double gapn_sl = 0.0;
      double gapn_ma = 0.0;
      GapN(sele, mele, sval_, mval_, gpn_, gapn_sl, gapn_ma);

      // evaluate the weighted gap (slave / master)
      GP_WGap(sele, lmval_, gapn_sl, gapn_ma, wgt, jacslave);

      WeakReset(linsize);
    }  // GP-loop
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype,
    IntPolicy>::IntegrateWeightedGapGradientError(MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, bool boundary_ele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo)
{
  // access unordered maps
  std::unordered_map<int, Deriv1stMap>* grad_error_ma_ptr =
      this->CParams().template GetUnorderedMap<int, Deriv1stMap>(0);
  std::unordered_map<int, Deriv1stMap>* grad_error_jac_ptr =
      this->CParams().template GetUnorderedMap<int, Deriv1stMap>(1);

  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);

  const int linsize = GetLinSize(sele);

  std::vector<unsigned> active_nlids;
  active_nlids.reserve(my::SLAVENUMNODE);
  ExtractActiveSlaveNodeLIDs(active_nlids, sele);

  {
    // get the gausspoints of this slave / master element pair
    const unsigned num_gps = projInfo.gaussPoints_.size();

    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    HardReset(linsize);

    for (my::gp_id_ = 0; static_cast<unsigned>(my::gp_id_) < num_gps; ++my::gp_id_)
    {
      const int gp = projInfo.gaussPoints_[my::gp_id_];

      // coordinates and weight
      const double eta[2] = {this->Wrapper().Coordinate(gp, 0), this->Wrapper().Coordinate(gp, 1)};
      const double wgt = this->Wrapper().Weight(gp) * projInfo.scaling_[my::gp_id_];

      // get Gauss point in slave element coordinates
      const double sxi[2] = {eta[0], eta[1]};
      const LINALG::Matrix<2, 1> sxi_mat(sxi, true);

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(this->ShapeFcn(), sxi, lmval_, lmderiv_, my::SLAVENUMNODE, true);

      // evaluate trace space shape functions (on both elements)
      shape_function_and_deriv1<slavetype>(sele, sxi_mat, sval_, sderiv_);

      // evaluate the convective slave base vectors
      LINALG::Matrix<3, 2> stau;
      sele.Metrics(sxi, &stau(0, 0), &stau(0, 1));

      // evaluate the two Jacobians (int. cell and slave element)
      const double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian 1-st and 2-nd order derivatives
      evaluator_->Deriv_Jacobian(sele, sxi, sderiv_, stau);

      const double uniqueProjalpha = projInfo.uniqueProjAlpha_[my::gp_id_];
      const LINALG::Matrix<2, 1>& uniqueMxi = projInfo.uniqueMxi_[my::gp_id_];

      mele.GetNodalCoords(mcoord_);

      // get mval and mderiv1
      shape_function_and_deriv1<mastertype>(mele, uniqueMxi, mval_, mderiv_);

      // evaluate the convective master base vectors
      LINALG::Matrix<3, 2> mtau;
      mele.Metrics(uniqueMxi.A(), &mtau(0, 0), &mtau(0, 1));

      // evaluate the GP master coordinate 1-st and 2-nd order derivatives
      evaluator_->Deriv_MXiGP(
          sele, mele, &sxi[0], uniqueMxi.A(), uniqueProjalpha, sval_, mval_, mderiv_, mtau);

      //**********************************************************************
      // evaluate at GP and lin char. quantities
      //**********************************************************************
      // calculate the averaged normal + derivative at gp level
      GP_Normal_DerivNormal(
          sele, sval_, &gpn_[0], dn_non_unit_, ddn_non_unit_, dn_unit_, ddn_unit_);

      // integrate the inner integral relating to the first order derivative of
      // the discrete normal gap for later usage (for all found slave nodes)
      IntPolicy::Get_Deriv1st_GapN(
          sele, mele, sval_, mval_, &gpn_[0], mtau, dmxigp_, deriv_gapn_sl_, deriv_gapn_ma_);

      // evaluate normal gap (split into slave and master contributions)
      double gapn_sl = 0.0;
      double gapn_ma = 0.0;
      GapN(sele, mele, sval_, mval_, gpn_, gapn_sl, gapn_ma);

      IntPolicy::Get_Deriv1st_WGapN_Error(sele, active_nlids, lmval_, gpn_, gapn_sl, gapn_ma, wgt,
          jacslave, derivjac_, mtau, dmxigp_, deriv_gapn_ma_, *grad_error_ma_ptr,
          *grad_error_jac_ptr);

      WeakReset(linsize);
    }  // GP-loop
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
int CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::GetLinSize(
    MORTAR::MortarElement& sele) const
{
  int linsize = 0;
  const DRT::Node* const* mynodes = sele.Nodes();
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    const CoNode& cnode = static_cast<const CoNode&>(*mynodes[i]);
    linsize += cnode.GetLinsize();
  }

  return linsize;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype,
    IntPolicy>::ExtractActiveSlaveNodeLIDs(std::vector<unsigned>& active_nlids,
    const MORTAR::MortarElement& sele) const
{
  const Epetra_Map* active_snode_row_map = this->CParams().template Get<Epetra_Map>(1);

  const int* nodeids = sele.NodeIds();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    if (active_snode_row_map->LID(nodeids[i]) != -1)
    {
      active_nlids.push_back(i);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::HardReset(
    const unsigned linsize)
{
  GEN_DATA::reset(my::SLAVEDIM, 0, dsxigp_);

  GEN_DATA::reset(my::MASTERDIM, linsize + my::MASTERNUMNODE * probdim, dmxigp_);
  GEN_DATA::reset(linsize + my::MASTERNUMNODE * probdim, dalpha_);
  GEN_DATA::reset(my::MASTERDIM, linsize + my::MASTERNUMNODE * probdim, ddmxigp_);

  std::fill(gpn_, gpn_ + 3, 0.0);
  GEN_DATA::reset(probdim, linsize + probdim * my::MASTERNUMNODE, dn_non_unit_);
  GEN_DATA::reset(probdim, linsize + probdim * my::MASTERNUMNODE, ddn_non_unit_);
  GEN_DATA::reset(probdim, linsize + probdim * my::MASTERNUMNODE, dn_unit_);
  GEN_DATA::reset(probdim, linsize + probdim * my::MASTERNUMNODE, ddn_unit_);

  GEN_DATA::reset(probdim * my::SLAVENUMNODE, deriv_gapn_sl_);
  GEN_DATA::reset(linsize + probdim * my::MASTERNUMNODE, deriv_gapn_ma_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype, class IntPolicy>
void CONTACT::AUG::Integrator<probdim, slavetype, mastertype, IntPolicy>::WeakReset(
    const unsigned linsize)
{
  GEN_DATA::reset(my::SLAVEDIM, 0, dsxigp_);

  GEN_DATA::weak_reset(dmxigp_);
  GEN_DATA::weak_reset(dalpha_);
  GEN_DATA::weak_reset(ddmxigp_);

  std::fill(gpn_, gpn_ + 3, 0.0);
  GEN_DATA::weak_reset(dn_non_unit_);
  GEN_DATA::weak_reset(ddn_non_unit_);
  GEN_DATA::reset(probdim, linsize + probdim * my::MASTERNUMNODE, dn_unit_);
  GEN_DATA::weak_reset(ddn_unit_);

  GEN_DATA::reset(probdim * my::SLAVENUMNODE, deriv_gapn_sl_);
  GEN_DATA::reset(linsize + probdim * my::MASTERNUMNODE, deriv_gapn_ma_);
}

/*----------------------------------------------------------------------------*/
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::line2>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::line2, DRT::Element::line2>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs2>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs2, DRT::Element::nurbs2>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs2, DRT::Element::nurbs3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs3>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs3, DRT::Element::nurbs3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::nurbs3, DRT::Element::nurbs2>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4, DRT::Element::quad4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4, DRT::Element::tri3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4, DRT::Element::nurbs4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4, DRT::Element::nurbs9>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3, DRT::Element::quad4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3, DRT::Element::tri3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3, DRT::Element::nurbs4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3, DRT::Element::nurbs9>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs4>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs4, DRT::Element::nurbs4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs4, DRT::Element::quad4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs4, DRT::Element::tri3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs4, DRT::Element::nurbs9>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs9>(
    DRT::Element::DiscretizationType mastertype, CONTACT::ParamsInterface& cparams,
    CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs9, DRT::Element::nurbs9>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs9, DRT::Element::quad4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs9, DRT::Element::tri3>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);
template CONTACT::AUG::IntegratorGeneric*
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::nurbs9, DRT::Element::nurbs4>(
    CONTACT::ParamsInterface& cparams, CONTACT::CoIntegrator* wrapper);

#include "contact_augmented_integrator_list.H"
