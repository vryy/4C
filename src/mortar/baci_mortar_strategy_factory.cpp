/*---------------------------------------------------------------------*/
/*! \file
\brief Base class for the CONTACT/MESHTYING factories.

\level 3

*/
/*---------------------------------------------------------------------*/
#include "baci_mortar_strategy_factory.H"

#include "baci_io.H"
#include "baci_io_pstream.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_mortar_element.H"
#include "baci_mortar_interface.H"
#include "baci_nurbs_discret.H"
#include "baci_nurbs_discret_control_point.H"
#include "baci_nurbs_discret_knotvector.H"
#include "baci_structure_new_timint_basedataglobalstate.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::STRATEGY::Factory::Factory()
    : discret_ptr_(Teuchos::null),
      isinit_(false),
      issetup_(false),
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
void MORTAR::STRATEGY::Factory::Init(Teuchos::RCP<DRT::Discretization> dis)
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
  if (gstate_ptr_ != Teuchos::null) discret_ptr_ = gstate_ptr_->GetDiscret();

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
DRT::Discretization& MORTAR::STRATEGY::Factory::Discret()
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& MORTAR::STRATEGY::Factory::Discret() const
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
void MORTAR::STRATEGY::Factory::PrepareNURBSElement(const DRT::Discretization& discret,
    Teuchos::RCP<DRT::Element> ele, Teuchos::RCP<MORTAR::MortarElement> cele) const
{
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discret));
  if (nurbsdis == nullptr) dserror("Dynamic cast failed!");

  Teuchos::RCP<const DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();
  std::vector<CORE::LINALG::SerialDenseVector> parentknots(Dim());
  std::vector<CORE::LINALG::SerialDenseVector> mortarknots(Dim() - 1);

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
    const DRT::Node* node, Teuchos::RCP<MORTAR::Node> mnode) const
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
      INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar, "LM_SHAPEFCN");
  INPAR::CONTACT::SystemType systype =
      INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact, "SYSTEM");
  INPAR::MORTAR::AlgorithmType algorithm =
      INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar, "ALGORITHM");
  bool nonSmoothGeometries = INPUT::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

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
      dserror("Invalid system type for contact/meshtying interface smoothing");
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
        else
          dserror("Invalid strategy or shape function type for contact/meshtying");
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
                 INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(smortar, "LM_QUAD") ==
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
          dserror("Invalid strategy or shape function type for contact/meshtying");
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
      dserror("Invalid system type for contact/meshtying");
  }
  return;
}

BACI_NAMESPACE_CLOSE
