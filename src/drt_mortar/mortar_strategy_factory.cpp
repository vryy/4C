/*---------------------------------------------------------------------*/
/*!
\file mortar_strategy_factory.cpp

\brief Base class for the CONTACT/MESHTYING factories.

\level 3

\maintainer Michael Hiermeier

\date Feb 4, 2016

*/
/*---------------------------------------------------------------------*/
#include "mortar_strategy_factory.H"

#include "mortar_element.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_knotvector.H"
#include "../drt_nurbs_discret/drt_control_point.H"

#include <Epetra_SerialDenseVector.h>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::STRATEGY::Factory::Factory()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null),
      discret_ptr_(Teuchos::null),
      comm_ptr_(Teuchos::null),
      dim_(-1)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr)
{
  // call Setup() after Init()
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::Setup()
{
  CheckInit();

  //  get the underlying discretization
  discret_ptr_ = gstate_ptr_->GetMutableDiscret();

  // get a copy of the underlying structural communicator
  comm_ptr_ = Teuchos::rcp(discret_ptr_->Comm().Clone());

  // get the problem dimension
  dim_ = DRT::Problem::Instance()->NDim();

  // since this is an abstract class the setup flag stays false.
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& MORTAR::STRATEGY::Factory::GState()
    const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::DiscretizationInterface& MORTAR::STRATEGY::Factory::Discret()
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::DiscretizationInterface& MORTAR::STRATEGY::Factory::Discret() const
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_Comm& MORTAR::STRATEGY::Factory::Comm()
{
  CheckInitSetup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::STRATEGY::Factory::Comm() const
{
  CheckInitSetup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Comm> MORTAR::STRATEGY::Factory::CommPtr()
{
  CheckInitSetup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Comm> MORTAR::STRATEGY::Factory::CommPtr() const
{
  CheckInitSetup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const int& MORTAR::STRATEGY::Factory::Dim() const
{
  if (dim_==-1)
    dserror("Call the STR::MODELEVEALUATOR::Setup() routine first to "
        "set the problem dimension variable!");
  return dim_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::PrepareNURBSElement(
    const DRT::DiscretizationInterface& discret,
    Teuchos::RCP<DRT::Element> ele,
    Teuchos::RCP<MORTAR::MortarElement> cele) const
{
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discret));
  if (nurbsdis==NULL)
    dserror("Dynamic cast failed!");

  Teuchos::RCP<const DRT::NURBS::Knotvector> knots =
      nurbsdis->GetKnotVector();
  std::vector<Epetra_SerialDenseVector> parentknots(Dim());
  std::vector<Epetra_SerialDenseVector> mortarknots(Dim() - 1);

  double normalfac = 0.0;
  Teuchos::RCP<DRT::FaceElement> faceele =
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);

  bool zero_size = knots->GetBoundaryEleAndParentKnots(
      parentknots,
      mortarknots,
      normalfac,
      faceele->ParentMasterElement()->Id(),
      faceele->FaceMasterNumber());

  // store nurbs specific data to node
  cele->ZeroSized() = zero_size;
  cele->Knots()     = mortarknots;
  cele->NormalFac() = normalfac;

  return;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::PrepareNURBSNode(const DRT::Node* node,
    Teuchos::RCP<MORTAR::MortarNode> mnode) const
{
  const DRT::NURBS::ControlPoint* cp =
      dynamic_cast<const DRT::NURBS::ControlPoint*>(node);

  mnode->NurbsW() = cp->W();

  return;
}
