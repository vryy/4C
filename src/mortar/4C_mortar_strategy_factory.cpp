/*---------------------------------------------------------------------*/
/*! \file
\brief Base class for the CONTACT/MESHTYING factories.

\level 3

*/
/*---------------------------------------------------------------------*/
#include "4C_mortar_strategy_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_nurbs_discret.hpp"
#include "4C_nurbs_discret_control_point.hpp"
#include "4C_nurbs_discret_knotvector.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

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
  check_init();

  //  get the underlying discretization
  if (gstate_ptr_ != Teuchos::null) discret_ptr_ = gstate_ptr_->GetDiscret();

  // get a copy of the underlying structural communicator
  comm_ptr_ = Teuchos::rcp(discret_ptr_->Comm().Clone());

  // get the problem dimension
  dim_ = GLOBAL::Problem::Instance()->NDim();

  // Note: Since this is an abstract class, the setup flag stays false.

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::check_init_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::check_init() const
{
  if (not is_init()) FOUR_C_THROW("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& MORTAR::STRATEGY::Factory::g_state() const
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& MORTAR::STRATEGY::Factory::discret()
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& MORTAR::STRATEGY::Factory::discret() const
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_Comm& MORTAR::STRATEGY::Factory::comm()
{
  check_init_setup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::STRATEGY::Factory::comm() const
{
  check_init_setup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Comm> MORTAR::STRATEGY::Factory::comm_ptr()
{
  check_init_setup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Comm> MORTAR::STRATEGY::Factory::comm_ptr() const
{
  check_init_setup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const int& MORTAR::STRATEGY::Factory::dim() const
{
  if (dim_ == -1)
    FOUR_C_THROW(
        "Call the STR::MODELEVEALUATOR::Setup() routine first to "
        "set the problem dimension variable!");
  return dim_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::CheckDimension() const
{
  if (dim() != 2 && dim() != 3) FOUR_C_THROW("Mortar meshtying/contact problems must be 2D or 3D");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::prepare_nurbs_element(const DRT::Discretization& discret,
    Teuchos::RCP<DRT::Element> ele, Teuchos::RCP<MORTAR::Element> cele) const
{
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discret));
  if (nurbsdis == nullptr) FOUR_C_THROW("Dynamic cast failed!");

  Teuchos::RCP<const DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();
  std::vector<CORE::LINALG::SerialDenseVector> parentknots(dim());
  std::vector<CORE::LINALG::SerialDenseVector> mortarknots(dim() - 1);

  double normalfac = 0.0;
  Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele, true);
  if (faceele.is_null()) FOUR_C_THROW("Cast to FaceElement failed!");

  bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots, normalfac,
      faceele->ParentMasterElement()->Id(), faceele->FaceMasterNumber());

  // store nurbs specific data to node
  cele->ZeroSized() = zero_size;
  cele->Knots() = mortarknots;
  cele->NormalFac() = normalfac;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::prepare_nurbs_node(
    const DRT::Node* node, Teuchos::RCP<MORTAR::Node> mnode) const
{
  const DRT::NURBS::ControlPoint* cp = dynamic_cast<const DRT::NURBS::ControlPoint*>(node);

  mnode->NurbsW() = cp->W();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::Factory::BuildSearchTree(
    const std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces) const
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
  const Teuchos::ParameterList& smortar = GLOBAL::Problem::Instance()->mortar_coupling_params();
  const Teuchos::ParameterList& scontact = GLOBAL::Problem::Instance()->contact_dynamic_params();
  INPAR::MORTAR::ShapeFcn shapefcn =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar, "LM_SHAPEFCN");
  INPAR::CONTACT::SystemType systype =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(scontact, "SYSTEM");
  INPAR::MORTAR::AlgorithmType algorithm =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar, "ALGORITHM");
  bool nonSmoothGeometries = CORE::UTILS::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

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
      FOUR_C_THROW("Invalid system type for contact/meshtying interface smoothing");
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
          FOUR_C_THROW("Invalid strategy or shape function type for contact/meshtying");
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
                 CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(smortar, "LM_QUAD") ==
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
          FOUR_C_THROW("Invalid strategy or shape function type for contact/meshtying");
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
      FOUR_C_THROW("Invalid system type for contact/meshtying");
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
