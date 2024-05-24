/*----------------------------------------------------------------------------*/
/*! \file
\brief GP projector template

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_contact_aug_projector.hpp"

#include "4C_contact_aug_element_utils.hpp"
#include "4C_linalg_gauss.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_element.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get(const unsigned probdim,
    CORE::FE::CellType ref_type, CORE::FE::CellType tar_type, const bool debug)
{
  switch (probdim)
  {
    case 2:
      return get2_d(ref_type, tar_type, debug);
    case 3:
      return get3_d(ref_type, tar_type, debug);
    default:
      FOUR_C_THROW("Unsupported problem dimension! (probdim=%d)", probdim);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::get2_d(
    CORE::FE::CellType ref_type, CORE::FE::CellType tar_type, const bool debug)
{
  switch (ref_type)
  {
    case CORE::FE::CellType::line2:
      return get2_d<CORE::FE::CellType::line2>(tar_type, debug);
    case CORE::FE::CellType::nurbs2:
      return get2_d<CORE::FE::CellType::nurbs2>(tar_type, debug);
    case CORE::FE::CellType::nurbs3:
      return get2_d<CORE::FE::CellType::nurbs3>(tar_type, debug);
    default:
      FOUR_C_THROW("Unsupported reference-type %s.", CORE::FE::CellTypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <CORE::FE::CellType ref_type>
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::get2_d(
    CORE::FE::CellType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case CORE::FE::CellType::line2:
      if (debug) return Projector<ProjDebugger, 2, ref_type, CORE::FE::CellType::line2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, CORE::FE::CellType::line2>::Instance();
    case CORE::FE::CellType::nurbs2:
      if (debug)
        return Projector<ProjDebugger, 2, ref_type, CORE::FE::CellType::nurbs2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, CORE::FE::CellType::nurbs2>::Instance();
    case CORE::FE::CellType::nurbs3:
      if (debug)
        return Projector<ProjDebugger, 2, ref_type, CORE::FE::CellType::nurbs3>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, CORE::FE::CellType::nurbs3>::Instance();
    default:
      FOUR_C_THROW("Unsupported target-type %s.", CORE::FE::CellTypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::get3_d(
    CORE::FE::CellType ref_type, CORE::FE::CellType tar_type, const bool debug)
{
  switch (ref_type)
  {
    case CORE::FE::CellType::quad4:
      return get3_d<CORE::FE::CellType::quad4>(tar_type, debug);
    case CORE::FE::CellType::tri3:
      return get3_d<CORE::FE::CellType::tri3>(tar_type, debug);
    case CORE::FE::CellType::nurbs4:
      return get3_d<CORE::FE::CellType::nurbs4>(tar_type, debug);
    case CORE::FE::CellType::nurbs9:
      return get3_d<CORE::FE::CellType::nurbs9>(tar_type, debug);
    default:
      FOUR_C_THROW("Unsupported reference-type %s.", CORE::FE::CellTypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <CORE::FE::CellType ref_type>
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::get3_d(
    CORE::FE::CellType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case CORE::FE::CellType::quad4:
      if (debug) return Projector<ProjDebugger, 3, ref_type, CORE::FE::CellType::quad4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, CORE::FE::CellType::quad4>::Instance();
    case CORE::FE::CellType::tri3:
      if (debug) return Projector<ProjDebugger, 3, ref_type, CORE::FE::CellType::tri3>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, CORE::FE::CellType::tri3>::Instance();
    case CORE::FE::CellType::nurbs4:
      if (debug)
        return Projector<ProjDebugger, 3, ref_type, CORE::FE::CellType::nurbs4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, CORE::FE::CellType::nurbs4>::Instance();
    case CORE::FE::CellType::nurbs9:
      if (debug)
        return Projector<ProjDebugger, 3, ref_type, CORE::FE::CellType::nurbs9>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, CORE::FE::CellType::nurbs9>::Instance();
    default:
      FOUR_C_THROW("Unsupported target-type %s.", CORE::FE::CellTypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
CONTACT::AUG::ProjectorBase*
CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::Instance()
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Projector<DebugPolicy, probdim, ref_type, tar_type>>(
            new Projector<DebugPolicy, probdim, ref_type, tar_type>);
      });

  return singleton_owner.Instance(CORE::UTILS::SingletonAction::create);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
void CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::Setup()
{
  ref_val_.putScalar(0.0);
  x_ref_.putScalar(0.0);
  n_ref_.putScalar(0.0);
  iter_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::operator()(
    MORTAR::Element& ref_ele, const double* ref_xi, MORTAR::Element& target_ele, double* target_xi,
    double& alpha)
{
  Setup();

  const DRT::Node* const* ref_nodes = ref_ele.Nodes();

  const CORE::LINALG::Matrix<REF_DIM, 1> rxi(ref_xi, true);

  shape_function<ref_type>(ref_ele, rxi, ref_val_);

  for (unsigned i = 0; i < REF_NUMNODES; ++i)
  {
    const MORTAR::Node& mnode = static_cast<const MORTAR::Node&>(*ref_nodes[i]);
    const double* x_ref_node = mnode.xspatial();
    const double* n_ref_node = mnode.MoData().n();

    for (unsigned j = 0; j < probdim; ++j)
    {
      x_ref_(j, 0) += x_ref_node[j] * ref_val_(i, 0);
      n_ref_(j, 0) += n_ref_node[j] * ref_val_(i, 0);
    }
  }

  // get the spatial coordinates of the target element
  const DRT::Node* const* tnodes = target_ele.Nodes();

  for (unsigned i = 0; i < TAR_NUMNODES; ++i)
  {
    const MORTAR::Node& mnode = static_cast<const MORTAR::Node&>(*tnodes[i]);
    const double* tar_x = mnode.xspatial();

    std::copy(tar_x, tar_x + probdim, &tar_coords_(0, i));
  }

  // initial value
  CORE::LINALG::Matrix<TAR_DIM, 1> txi(target_xi, true);
  alpha = 0.0;

  CORE::LINALG::Matrix<TAR_DIM, 1> txi_center(false);
  CORE::FE::getLocalCenterPosition<TAR_DIM>(tar_type, txi_center);
  std::copy(txi_center.A(), txi_center.A() + TAR_DIM, txi.A());

  rhs_gp(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha);
  DebugPolicy::writeVector(std::cout, probdim, rhs_.A(), "Rhs");

  const double ref_rhs_nrm2 = std::max(1.0, rhs_.Norm2());
  double dx_nrm2 = 0.0;
  double rhs_nrm2 = 0.0;
  bool is_parallel_proj = false;

  for (; iter_ < MORTARMAXITER; ++iter_)
  {
    l_mat_gp(lmat_, tar_deriv1_, target_ele, tar_coords_, target_xi, n_ref_);

    rhs_.Scale(-1.0);
    const double det = CORE::LINALG::gaussElimination<true, probdim>(lmat_, rhs_, dx_);

    // safety check
    if (std::abs(det) < 1.0e-12)
    {
      is_parallel_proj = true;
      std::cout << "WARNING: GPProjection parallel to master element:\n"
                << "Determinant:           " << det << "\n"
                << "Reference element GID: " << ref_ele.Id() << "\n"
                << "Target element GID:    " << target_ele.Id() << "\n"
                << "GP will be skipped for this master element!" << std::endl;
      break;
    }

    // update
    for (unsigned i = 0; i < TAR_DIM; ++i) txi(i, 0) += dx_(i, 0);
    alpha += dx_(probdim - 1, 0);

    DebugPolicy::writeVector(std::cout, TAR_DIM, txi.A(), "txi");

    // new right-hand side
    if (not rhs_gp(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha))
    {
      DebugPolicy::writeVector(std::cout, probdim, rhs_.A(), "Rhs (failed)");
      //      std::cout << "ShapeFunction evaluation failed @ [" << txi(0,0) << ", " <<
      //          txi(1,0) << "]" << std::endl;
      iter_ = MORTARMAXITER;
      break;
    }
    DebugPolicy::writeVector(std::cout, probdim, rhs_.A(), "Rhs");

    rhs_nrm2 = rhs_.Norm2();

    if (rhs_nrm2 <= MORTARCONVTOL * ref_rhs_nrm2)
    {
      dx_nrm2 = dx_.Norm2();

      rel_sol_tolerance_ = alpha * alpha;
      for (unsigned i = 0; i < TAR_DIM; ++i) rel_sol_tolerance_ += txi(i, 0) * txi(i, 0);

      rel_sol_tolerance_ = std::max(1.0, std::sqrt(rel_sol_tolerance_));
      rel_sol_tolerance_ *= MORTARCONVTOL;

      if (dx_nrm2 <= rel_sol_tolerance_) break;
    }
  }

  // Newton iteration unconverged
  if (iter_ == MORTARMAXITER or is_parallel_proj)
  {
    std::fill(txi.A(), txi.A() + TAR_DIM, 1.0e12);

    return false;
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::rhs_gp(
    CORE::LINALG::Matrix<probdim, 1>& rhs, const CORE::LINALG::Matrix<probdim, 1>& x_ref,
    const CORE::LINALG::Matrix<probdim, 1>& n_ref, MORTAR::Element& target_ele,
    const CORE::LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
    const double& alpha) const
{
  CORE::LINALG::Matrix<probdim, 1> x_tar(false);
  const bool status = get_global_position<tar_type>(target_ele, tar_coords, tar_xi, x_tar);

  // evaluate right hand side
  for (unsigned i = 0; i < probdim; ++i)
    rhs(i, 0) = x_tar(i, 0) - alpha * n_ref(i, 0) - x_ref(i, 0);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
template <CORE::FE::CellType type, unsigned numnodes>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::get_global_position(
    MORTAR::Element& ele, const CORE::LINALG::Matrix<probdim, numnodes>& coords, const double* xi,
    CORE::LINALG::Matrix<probdim, 1>& pos) const
{
  const unsigned dim = CORE::FE::dim<type>;
  const CORE::LINALG::Matrix<dim, 1> mat_xi(xi, true);

  CORE::LINALG::Matrix<numnodes, 1> val(true);
  const bool status = shape_function<type>(ele, mat_xi, val);

  pos.Multiply(coords, val);
  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
    CORE::FE::CellType tar_type>
void CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::l_mat_gp(
    CORE::LINALG::Matrix<probdim, probdim>& lmat,
    CORE::LINALG::Matrix<TAR_DIM, TAR_NUMNODES>& tar_deriv1, MORTAR::Element& tar_ele,
    const CORE::LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
    const CORE::LINALG::Matrix<probdim, 1>& n_ref) const
{
  const CORE::LINALG::Matrix<TAR_DIM, 1> txi(tar_xi, true);
  std::fill(lmat.A(), lmat.A() + 2 * probdim, 0.0);

  shape_function_deriv1<tar_type>(tar_ele, txi, tar_deriv1);

  for (unsigned i = 0; i < TAR_NUMNODES; ++i)
    for (unsigned j = 0; j < TAR_DIM; ++j)
      for (unsigned k = 0; k < probdim; ++k) lmat(k, j) += tar_deriv1(j, i) * tar_coords(k, i);

  for (unsigned i = 0; i < probdim; ++i) lmat(i, probdim - 1) = -n_ref(i, 0);
}


template CONTACT::AUG::ProjectorBase*
CONTACT::AUG::ProjectorBase::get2_d<CORE::FE::CellType::line2>(
    CORE::FE::CellType tar_type, const bool debug);
template CONTACT::AUG::ProjectorBase*
CONTACT::AUG::ProjectorBase::get3_d<CORE::FE::CellType::quad4>(
    CORE::FE::CellType tar_type, const bool debug);
template CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::get3_d<CORE::FE::CellType::tri3>(
    CORE::FE::CellType tar_type, const bool debug);

// standard discretization types
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 2, CORE::FE::CellType::line2,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 2,
    CORE::FE::CellType::line2, CORE::FE::CellType::line2>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::quad4, CORE::FE::CellType::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::quad4>;

// pure NURBS discretization types
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs2>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3,
    CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs4>;

FOUR_C_NAMESPACE_CLOSE
