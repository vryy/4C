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
CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::Get(const unsigned probdim,
    Core::FE::CellType ref_type, Core::FE::CellType tar_type, const bool debug)
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
CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::get2_d(
    Core::FE::CellType ref_type, Core::FE::CellType tar_type, const bool debug)
{
  switch (ref_type)
  {
    case Core::FE::CellType::line2:
      return get2_d<Core::FE::CellType::line2>(tar_type, debug);
    case Core::FE::CellType::nurbs2:
      return get2_d<Core::FE::CellType::nurbs2>(tar_type, debug);
    case Core::FE::CellType::nurbs3:
      return get2_d<Core::FE::CellType::nurbs3>(tar_type, debug);
    default:
      FOUR_C_THROW("Unsupported reference-type %s.", Core::FE::CellTypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType ref_type>
CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::get2_d(
    Core::FE::CellType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case Core::FE::CellType::line2:
      if (debug) return Projector<ProjDebugger, 2, ref_type, Core::FE::CellType::line2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, Core::FE::CellType::line2>::Instance();
    case Core::FE::CellType::nurbs2:
      if (debug)
        return Projector<ProjDebugger, 2, ref_type, Core::FE::CellType::nurbs2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, Core::FE::CellType::nurbs2>::Instance();
    case Core::FE::CellType::nurbs3:
      if (debug)
        return Projector<ProjDebugger, 2, ref_type, Core::FE::CellType::nurbs3>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, Core::FE::CellType::nurbs3>::Instance();
    default:
      FOUR_C_THROW("Unsupported target-type %s.", Core::FE::CellTypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::get3_d(
    Core::FE::CellType ref_type, Core::FE::CellType tar_type, const bool debug)
{
  switch (ref_type)
  {
    case Core::FE::CellType::quad4:
      return get3_d<Core::FE::CellType::quad4>(tar_type, debug);
    case Core::FE::CellType::tri3:
      return get3_d<Core::FE::CellType::tri3>(tar_type, debug);
    case Core::FE::CellType::nurbs4:
      return get3_d<Core::FE::CellType::nurbs4>(tar_type, debug);
    case Core::FE::CellType::nurbs9:
      return get3_d<Core::FE::CellType::nurbs9>(tar_type, debug);
    default:
      FOUR_C_THROW("Unsupported reference-type %s.", Core::FE::CellTypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType ref_type>
CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::get3_d(
    Core::FE::CellType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case Core::FE::CellType::quad4:
      if (debug) return Projector<ProjDebugger, 3, ref_type, Core::FE::CellType::quad4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, Core::FE::CellType::quad4>::Instance();
    case Core::FE::CellType::tri3:
      if (debug) return Projector<ProjDebugger, 3, ref_type, Core::FE::CellType::tri3>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, Core::FE::CellType::tri3>::Instance();
    case Core::FE::CellType::nurbs4:
      if (debug)
        return Projector<ProjDebugger, 3, ref_type, Core::FE::CellType::nurbs4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, Core::FE::CellType::nurbs4>::Instance();
    case Core::FE::CellType::nurbs9:
      if (debug)
        return Projector<ProjDebugger, 3, ref_type, Core::FE::CellType::nurbs9>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, Core::FE::CellType::nurbs9>::Instance();
    default:
      FOUR_C_THROW("Unsupported target-type %s.", Core::FE::CellTypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
CONTACT::Aug::ProjectorBase*
CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::Instance()
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Projector<DebugPolicy, probdim, ref_type, tar_type>>(
            new Projector<DebugPolicy, probdim, ref_type, tar_type>);
      });

  return singleton_owner.Instance(Core::UTILS::SingletonAction::create);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
void CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::setup()
{
  ref_val_.putScalar(0.0);
  x_ref_.putScalar(0.0);
  n_ref_.putScalar(0.0);
  iter_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
bool CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::operator()(
    Mortar::Element& ref_ele, const double* ref_xi, Mortar::Element& target_ele, double* target_xi,
    double& alpha)
{
  setup();

  const Core::Nodes::Node* const* ref_nodes = ref_ele.Nodes();

  const Core::LinAlg::Matrix<REF_DIM, 1> rxi(ref_xi, true);

  shape_function<ref_type>(ref_ele, rxi, ref_val_);

  for (unsigned i = 0; i < REF_NUMNODES; ++i)
  {
    const Mortar::Node& mnode = static_cast<const Mortar::Node&>(*ref_nodes[i]);
    const double* x_ref_node = mnode.xspatial();
    const double* n_ref_node = mnode.MoData().n();

    for (unsigned j = 0; j < probdim; ++j)
    {
      x_ref_(j, 0) += x_ref_node[j] * ref_val_(i, 0);
      n_ref_(j, 0) += n_ref_node[j] * ref_val_(i, 0);
    }
  }

  // get the spatial coordinates of the target element
  const Core::Nodes::Node* const* tnodes = target_ele.Nodes();

  for (unsigned i = 0; i < TAR_NUMNODES; ++i)
  {
    const Mortar::Node& mnode = static_cast<const Mortar::Node&>(*tnodes[i]);
    const double* tar_x = mnode.xspatial();

    std::copy(tar_x, tar_x + probdim, &tar_coords_(0, i));
  }

  // initial value
  Core::LinAlg::Matrix<TAR_DIM, 1> txi(target_xi, true);
  alpha = 0.0;

  Core::LinAlg::Matrix<TAR_DIM, 1> txi_center(false);
  Core::FE::getLocalCenterPosition<TAR_DIM>(tar_type, txi_center);
  std::copy(txi_center.A(), txi_center.A() + TAR_DIM, txi.A());

  rhs_gp(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha);
  DebugPolicy::write_vector(std::cout, probdim, rhs_.A(), "Rhs");

  const double ref_rhs_nrm2 = std::max(1.0, rhs_.Norm2());
  double dx_nrm2 = 0.0;
  double rhs_nrm2 = 0.0;
  bool is_parallel_proj = false;

  for (; iter_ < MORTARMAXITER; ++iter_)
  {
    l_mat_gp(lmat_, tar_deriv1_, target_ele, tar_coords_, target_xi, n_ref_);

    rhs_.Scale(-1.0);
    const double det = Core::LinAlg::gaussElimination<true, probdim>(lmat_, rhs_, dx_);

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

    DebugPolicy::write_vector(std::cout, TAR_DIM, txi.A(), "txi");

    // new right-hand side
    if (not rhs_gp(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha))
    {
      DebugPolicy::write_vector(std::cout, probdim, rhs_.A(), "Rhs (failed)");
      //      std::cout << "shape_function evaluation failed @ [" << txi(0,0) << ", " <<
      //          txi(1,0) << "]" << std::endl;
      iter_ = MORTARMAXITER;
      break;
    }
    DebugPolicy::write_vector(std::cout, probdim, rhs_.A(), "Rhs");

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
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
bool CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::rhs_gp(
    Core::LinAlg::Matrix<probdim, 1>& rhs, const Core::LinAlg::Matrix<probdim, 1>& x_ref,
    const Core::LinAlg::Matrix<probdim, 1>& n_ref, Mortar::Element& target_ele,
    const Core::LinAlg::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
    const double& alpha) const
{
  Core::LinAlg::Matrix<probdim, 1> x_tar(false);
  const bool status = get_global_position<tar_type>(target_ele, tar_coords, tar_xi, x_tar);

  // evaluate right hand side
  for (unsigned i = 0; i < probdim; ++i)
    rhs(i, 0) = x_tar(i, 0) - alpha * n_ref(i, 0) - x_ref(i, 0);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
template <Core::FE::CellType type, unsigned numnodes>
bool CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::get_global_position(
    Mortar::Element& ele, const Core::LinAlg::Matrix<probdim, numnodes>& coords, const double* xi,
    Core::LinAlg::Matrix<probdim, 1>& pos) const
{
  const unsigned dim = Core::FE::dim<type>;
  const Core::LinAlg::Matrix<dim, 1> mat_xi(xi, true);

  Core::LinAlg::Matrix<numnodes, 1> val(true);
  const bool status = shape_function<type>(ele, mat_xi, val);

  pos.Multiply(coords, val);
  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
    Core::FE::CellType tar_type>
void CONTACT::Aug::Projector<DebugPolicy, probdim, ref_type, tar_type>::l_mat_gp(
    Core::LinAlg::Matrix<probdim, probdim>& lmat,
    Core::LinAlg::Matrix<TAR_DIM, TAR_NUMNODES>& tar_deriv1, Mortar::Element& tar_ele,
    const Core::LinAlg::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
    const Core::LinAlg::Matrix<probdim, 1>& n_ref) const
{
  const Core::LinAlg::Matrix<TAR_DIM, 1> txi(tar_xi, true);
  std::fill(lmat.A(), lmat.A() + 2 * probdim, 0.0);

  shape_function_deriv1<tar_type>(tar_ele, txi, tar_deriv1);

  for (unsigned i = 0; i < TAR_NUMNODES; ++i)
    for (unsigned j = 0; j < TAR_DIM; ++j)
      for (unsigned k = 0; k < probdim; ++k) lmat(k, j) += tar_deriv1(j, i) * tar_coords(k, i);

  for (unsigned i = 0; i < probdim; ++i) lmat(i, probdim - 1) = -n_ref(i, 0);
}


template CONTACT::Aug::ProjectorBase*
CONTACT::Aug::ProjectorBase::get2_d<Core::FE::CellType::line2>(
    Core::FE::CellType tar_type, const bool debug);
template CONTACT::Aug::ProjectorBase*
CONTACT::Aug::ProjectorBase::get3_d<Core::FE::CellType::quad4>(
    Core::FE::CellType tar_type, const bool debug);
template CONTACT::Aug::ProjectorBase* CONTACT::Aug::ProjectorBase::get3_d<Core::FE::CellType::tri3>(
    Core::FE::CellType tar_type, const bool debug);

// standard discretization types
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 2, Core::FE::CellType::line2,
    Core::FE::CellType::line2>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 2,
    Core::FE::CellType::line2, Core::FE::CellType::line2>;

template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::quad4,
    Core::FE::CellType::quad4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::quad4, Core::FE::CellType::quad4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::quad4,
    Core::FE::CellType::tri3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::quad4, Core::FE::CellType::tri3>;

template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::tri3,
    Core::FE::CellType::tri3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3, Core::FE::CellType::tri3,
    Core::FE::CellType::tri3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::tri3,
    Core::FE::CellType::quad4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3, Core::FE::CellType::tri3,
    Core::FE::CellType::quad4>;

// pure NURBS discretization types
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs2,
    Core::FE::CellType::nurbs2>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs2>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs3,
    Core::FE::CellType::nurbs3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs2,
    Core::FE::CellType::nurbs3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs3>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs3,
    Core::FE::CellType::nurbs2>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs2>;

template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs4,
    Core::FE::CellType::nurbs4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs9,
    Core::FE::CellType::nurbs9>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs9>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs4,
    Core::FE::CellType::nurbs9>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs9>;
template class CONTACT::Aug::Projector<CONTACT::Aug::ProjDebugger, 3, Core::FE::CellType::nurbs9,
    Core::FE::CellType::nurbs4>;
template class CONTACT::Aug::Projector<CONTACT::Aug::EmptyProjDebugger, 3,
    Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs4>;

FOUR_C_NAMESPACE_CLOSE
