/*---------------------------------------------------------------------*/
/*!
\file mortar_strategy_factory.cpp

\brief Base class for the CONTACT/MESHTYING factories.

\level 3

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#include "mortar_strategy_factory.H"

#include "mortar_element.H"
#include "mortar_interface.H"

#include "../drt_inpar/inpar_contact.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_nurbs_discret/drt_knotvector.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include <Epetra_SerialDenseVector.h>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::STRATEGY::Factory::Factory()
    : isinit_(false),
      issetup_(false),
      discret_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
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
void MORTAR::STRATEGY::Factory::Init(Teuchos::RCP<DRT::DiscretizationInterface> dis)
{
  // call Setup() after Init()
  issetup_ = false;

  gstate_ptr_ = Teuchos::null;

  discret_ptr_ = dis;

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::Setup()
{
  CheckInit();

  //  get the underlying discretization
  if (gstate_ptr_ != Teuchos::null) discret_ptr_ = gstate_ptr_->GetMutableDiscret();

  // get a copy of the underlying structural communicator
  comm_ptr_ = Teuchos::rcp(discret_ptr_->Comm().Clone());

  // get the problem dimension
  dim_ = DRT::Problem::Instance()->NDim();

  // Note: Since this is an abstract class, the setup flag stays false.

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckInit() const
{
  if (not IsInit()) dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& MORTAR::STRATEGY::Factory::GState() const
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
  if (dim_ == -1)
    dserror(
        "Call the STR::MODELEVEALUATOR::Setup() routine first to "
        "set the problem dimension variable!");
  return dim_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckDimension() const
{
  if (Dim() != 2 && Dim() != 3) dserror("Mortar meshtying/contact problems must be 2D or 3D");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::PrepareNURBSElement(const DRT::DiscretizationInterface& discret,
    Teuchos::RCP<DRT::Element> ele, Teuchos::RCP<MORTAR::MortarElement> cele) const
{
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discret));
  if (nurbsdis == NULL) dserror("Dynamic cast failed!");

  Teuchos::RCP<const DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();
  std::vector<Epetra_SerialDenseVector> parentknots(Dim());
  std::vector<Epetra_SerialDenseVector> mortarknots(Dim() - 1);

  double normalfac = 0.0;
  Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele, true);
  if (faceele.is_null()) dserror("Cast to FaceElement failed!");

  bool zero_size = knots->GetBoundaryEleAndParentKnots(parentknots, mortarknots, normalfac,
      faceele->ParentMasterElement()->Id(), faceele->FaceMasterNumber());

  // store nurbs specific data to node
  cele->ZeroSized() = zero_size;
  cele->Knots() = mortarknots;
  cele->NormalFac() = normalfac;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::PrepareNURBSNode(
    const DRT::Node* node, Teuchos::RCP<MORTAR::MortarNode> mnode) const
{
  const DRT::NURBS::ControlPoint* cp = dynamic_cast<const DRT::NURBS::ControlPoint*>(node);

  mnode->NurbsW() = cp->W();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::BuildSearchTree(
    const std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const
{
  for (unsigned i = 0; i < interfaces.size(); ++i) interfaces[i]->CreateSearchTree();

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::PrintStrategyBanner(
    const enum INPAR::CONTACT::SolvingStrategy soltype)
{
  // some parameters
  const Teuchos::ParameterList& smortar = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar, "LM_SHAPEFCN");
  INPAR::CONTACT::SystemType systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact, "SYSTEM");
  INPAR::MORTAR::AlgorithmType algorithm =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar, "ALGORITHM");
  bool smoothing = DRT::INPUT::IntegralValue<int>(scontact, "DISCR_SMOOTHING");
  bool nonSmoothGeometries = DRT::INPUT::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

  if (nonSmoothGeometries)
  {
    if (soltype == INPAR::CONTACT::solution_lagmult)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Lagrange Multiplier Strategy =============================\n";
      IO::cout << "===== NONSMOOTH - GEOMETRIES ===================================\n";
      IO::cout << "================================================================\n\n";
    }
    else
      dserror("ERROR: Invalid system type for contact/meshtying interface smoothing");
  }
  else if (smoothing)
  {
    if (soltype == INPAR::CONTACT::solution_lagmult)
    {
      IO::cout << "================================================================\n";
      IO::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      IO::cout << "================================================================\n\n";
      IO::cout << "================================================================\n";
      IO::cout << "===== Interface smoothing approach with     ====================\n";
      IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
      IO::cout << "===== (Saddle point formulation) ===============================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (soltype == INPAR::CONTACT::solution_penalty)
    {
      IO::cout << "================================================================\n";
      IO::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      IO::cout << "================================================================\n\n";
      IO::cout << "================================================================\n";
      IO::cout << "===== Interface smoothing approach with     ====================\n";
      IO::cout << "===== Standard Penalty strategy             ====================\n";
      IO::cout << "===== (Pure displacement formulation)===========================\n";
      IO::cout << "================================================================\n\n";
    }
    else
      dserror("ERROR: Invalid system type for contact/meshtying interface smoothing");
  }
  else
  {
    if (algorithm == INPAR::MORTAR::algorithm_mortar)
    {
      // saddle point formulation
      if (systype == INPAR::CONTACT::system_saddlepoint)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult &&
            shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_combo)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Combination of different Solving Strategies ==============\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_augmented)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Augmented Lagrange strategy ==============================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_std_lagrange)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange strategy ===============================\n";
          IO::cout << "===== Derived from the Augmented formulation ===================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_steepest_ascent)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Steepest Ascent strategy =================================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_xcontact &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Extended contact strategy ================================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }

      // condensed formulation
      else if (systype == INPAR::CONTACT::system_condensed ||
               systype == INPAR::CONTACT::system_condensed_lagmult)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_standard &&
                 DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(smortar, "LM_QUAD") ==
                     INPAR::MORTAR::lagmult_const)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== const Lagrange multiplier strategy =======================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }
    }
    else if (algorithm == INPAR::MORTAR::algorithm_nts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Node-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_lts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Line-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_stl)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Segment-To-Line approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_gpts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      IO::cout << "================================================================\n\n";
    }
    // invalid system type
    else
      dserror("ERROR: Invalid system type for contact/meshtying");
  }
  return;
}
