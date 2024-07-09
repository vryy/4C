/*---------------------------------------------------------------------*/
/*! \file
\brief Base class for the CONTACT/MESHTYING factories.

\level 3

*/
/*---------------------------------------------------------------------*/
#include "4C_mortar_strategy_factory.hpp"

#include "4C_fem_nurbs_discretization.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mortar::STRATEGY::Factory::Factory()
    : discret_ptr_(Teuchos::null),
      isinit_(false),
      issetup_(false),
      comm_ptr_(Teuchos::null),
      dim_(-1)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::init(Teuchos::RCP<Core::FE::Discretization> dis)
{
  // call setup() after init()
  issetup_ = false;

  discret_ptr_ = dis;

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::setup()
{
  check_init();

  // get a copy of the underlying structural communicator
  comm_ptr_ = Teuchos::rcp(discret_ptr_->get_comm().Clone());

  // get the problem dimension
  dim_ = Global::Problem::instance()->n_dim();

  // Note: Since this is an abstract class, the setup flag stays false.

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::check_init_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::check_init() const
{
  if (not is_init()) FOUR_C_THROW("Call init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::Discretization& Mortar::STRATEGY::Factory::discret()
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Core::FE::Discretization& Mortar::STRATEGY::Factory::discret() const
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_Comm& Mortar::STRATEGY::Factory::get_comm()
{
  check_init_setup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& Mortar::STRATEGY::Factory::get_comm() const
{
  check_init_setup();
  return *comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Comm> Mortar::STRATEGY::Factory::comm_ptr()
{
  check_init_setup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Comm> Mortar::STRATEGY::Factory::comm_ptr() const
{
  check_init_setup();
  return comm_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const int& Mortar::STRATEGY::Factory::n_dim() const
{
  if (dim_ == -1)
    FOUR_C_THROW(
        "Call the Solid::MODELEVEALUATOR::setup() routine first to "
        "set the problem dimension variable!");
  return dim_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::check_dimension() const
{
  if (n_dim() != 2 && n_dim() != 3)
    FOUR_C_THROW("Mortar meshtying/contact problems must be 2D or 3D");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::prepare_nurbs_element(const Core::FE::Discretization& discret,
    Teuchos::RCP<Core::Elements::Element> ele, Teuchos::RCP<Mortar::Element> cele) const
{
  const Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<const Core::FE::Nurbs::NurbsDiscretization*>(&(discret));
  if (nurbsdis == nullptr) FOUR_C_THROW("Dynamic cast failed!");

  Teuchos::RCP<const Core::FE::Nurbs::Knotvector> knots = nurbsdis->get_knot_vector();
  std::vector<Core::LinAlg::SerialDenseVector> parentknots(n_dim());
  std::vector<Core::LinAlg::SerialDenseVector> mortarknots(n_dim() - 1);

  double normalfac = 0.0;
  Teuchos::RCP<Core::Elements::FaceElement> faceele =
      Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
  if (faceele.is_null()) FOUR_C_THROW("Cast to FaceElement failed!");

  bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots, normalfac,
      faceele->parent_master_element()->id(), faceele->face_master_number());

  // store nurbs specific data to node
  cele->zero_sized() = zero_size;
  cele->knots() = mortarknots;
  cele->normal_fac() = normalfac;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::prepare_nurbs_node(
    const Core::Nodes::Node* node, Teuchos::RCP<Mortar::Node> mnode) const
{
  const Core::FE::Nurbs::ControlPoint* cp =
      dynamic_cast<const Core::FE::Nurbs::ControlPoint*>(node);

  mnode->nurbs_w() = cp->w();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::build_search_tree(
    const std::vector<Teuchos::RCP<Mortar::Interface>>& interfaces) const
{
  for (unsigned i = 0; i < interfaces.size(); ++i) interfaces[i]->create_search_tree();

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::STRATEGY::Factory::print_strategy_banner(
    const enum Inpar::CONTACT::SolvingStrategy soltype)
{
  // some parameters
  const Teuchos::ParameterList& smortar = Global::Problem::instance()->mortar_coupling_params();
  const Teuchos::ParameterList& scontact = Global::Problem::instance()->contact_dynamic_params();
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(smortar, "LM_SHAPEFCN");
  Inpar::CONTACT::SystemType systype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(scontact, "SYSTEM");
  Inpar::Mortar::AlgorithmType algorithm =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(smortar, "ALGORITHM");
  bool nonSmoothGeometries = Core::UTILS::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

  if (nonSmoothGeometries)
  {
    if (soltype == Inpar::CONTACT::solution_lagmult)
    {
      Core::IO::cout << "================================================================\n";
      Core::IO::cout << "===== Lagrange Multiplier Strategy =============================\n";
      Core::IO::cout << "===== NONSMOOTH - GEOMETRIES ===================================\n";
      Core::IO::cout << "================================================================\n\n";
    }
    else
      FOUR_C_THROW("Invalid system type for contact/meshtying interface smoothing");
  }
  else
  {
    if (algorithm == Inpar::Mortar::algorithm_mortar)
    {
      // saddle point formulation
      if (systype == Inpar::CONTACT::system_saddlepoint)
      {
        if (soltype == Inpar::CONTACT::solution_lagmult &&
            shapefcn == Inpar::Mortar::shape_standard)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
          Core::IO::cout << "===== (Saddle point formulation) ===============================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_lagmult &&
                 shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          Core::IO::cout << "===== (Saddle point formulation) ===============================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_lagmult &&
                 shapefcn == Inpar::Mortar::shape_petrovgalerkin)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          Core::IO::cout << "===== (Saddle point formulation) ===============================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_penalty &&
                 shapefcn == Inpar::Mortar::shape_standard)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Standard Penalty strategy ================================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_penalty &&
                 shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Penalty strategy ====================================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_uzawa &&
                 shapefcn == Inpar::Mortar::shape_standard)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_uzawa && shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else
          FOUR_C_THROW("Invalid strategy or shape function type for contact/meshtying");
      }

      // condensed formulation
      else if (systype == Inpar::CONTACT::system_condensed ||
               systype == Inpar::CONTACT::system_condensed_lagmult)
      {
        if (soltype == Inpar::CONTACT::solution_lagmult && shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          Core::IO::cout << "===== (Condensed formulation) ==================================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_lagmult &&
                 shapefcn == Inpar::Mortar::shape_standard &&
                 Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(smortar, "LM_QUAD") ==
                     Inpar::Mortar::lagmult_const)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== const Lagrange multiplier strategy =======================\n";
          Core::IO::cout << "===== (Condensed formulation) ==================================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_lagmult &&
                 shapefcn == Inpar::Mortar::shape_petrovgalerkin)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          Core::IO::cout << "===== (Condensed formulation) ==================================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_penalty &&
                 shapefcn == Inpar::Mortar::shape_standard)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Standard Penalty strategy ================================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_penalty &&
                 shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Penalty strategy ====================================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_uzawa &&
                 shapefcn == Inpar::Mortar::shape_standard)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else if (soltype == Inpar::CONTACT::solution_uzawa && shapefcn == Inpar::Mortar::shape_dual)
        {
          Core::IO::cout << "================================================================\n";
          Core::IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          Core::IO::cout << "===== (Pure displacement formulation) ==========================\n";
          Core::IO::cout << "================================================================\n\n";
        }
        else
          FOUR_C_THROW("Invalid strategy or shape function type for contact/meshtying");
      }
    }
    else if (algorithm == Inpar::Mortar::algorithm_nts)
    {
      Core::IO::cout << "================================================================\n";
      Core::IO::cout << "===== Node-To-Segment approach =================================\n";
      Core::IO::cout << "================================================================\n\n";
    }
    else if (algorithm == Inpar::Mortar::algorithm_lts)
    {
      Core::IO::cout << "================================================================\n";
      Core::IO::cout << "===== Line-To-Segment approach =================================\n";
      Core::IO::cout << "================================================================\n\n";
    }
    else if (algorithm == Inpar::Mortar::algorithm_stl)
    {
      Core::IO::cout << "================================================================\n";
      Core::IO::cout << "===== Segment-To-Line approach =================================\n";
      Core::IO::cout << "================================================================\n\n";
    }
    else if (algorithm == Inpar::Mortar::algorithm_gpts)
    {
      Core::IO::cout << "================================================================\n";
      Core::IO::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      Core::IO::cout << "================================================================\n\n";
    }
    // invalid system type
    else
      FOUR_C_THROW("Invalid system type for contact/meshtying");
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
