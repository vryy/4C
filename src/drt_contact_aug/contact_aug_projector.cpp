/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_projector.cpp

\brief GP projector template

\level 3

\maintainer Michael Hiermeier
\date Jun 23, 2017

*/
/*----------------------------------------------------------------------------*/

#include "contact_aug_projector.H"
#include "../drt_mortar/mortar_element.H"
#include "contact_aug_element_utils.H"
#include "../linalg/linalg_gauss.H"
#include "../drt_mortar/mortar_defines.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get(const unsigned probdim,
    DRT::Element::DiscretizationType ref_type, DRT::Element::DiscretizationType tar_type,
    const bool debug)
{
  switch (probdim)
  {
    case 2:
      return Get2D(ref_type, tar_type, debug);
    case 3:
      return Get3D(ref_type, tar_type, debug);
    default:
      dserror("Unsupported problem dimension! (probdim=%d)", probdim);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get2D(
    DRT::Element::DiscretizationType ref_type, DRT::Element::DiscretizationType tar_type,
    const bool debug)
{
  switch (ref_type)
  {
    case DRT::Element::line2:
      return Get2D<DRT::Element::line2>(tar_type, debug);
    case DRT::Element::nurbs2:
      return Get2D<DRT::Element::nurbs2>(tar_type, debug);
    case DRT::Element::nurbs3:
      return Get2D<DRT::Element::nurbs3>(tar_type, debug);
    default:
      dserror("Unsupported reference-type %s.", DRT::DistypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType ref_type>
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get2D(
    DRT::Element::DiscretizationType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case DRT::Element::line2:
      if (debug) return Projector<ProjDebugger, 2, ref_type, DRT::Element::line2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, DRT::Element::line2>::Instance();
    case DRT::Element::nurbs2:
      if (debug) return Projector<ProjDebugger, 2, ref_type, DRT::Element::nurbs2>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, DRT::Element::nurbs2>::Instance();
    case DRT::Element::nurbs3:
      if (debug) return Projector<ProjDebugger, 2, ref_type, DRT::Element::nurbs3>::Instance();
      return Projector<EmptyProjDebugger, 2, ref_type, DRT::Element::nurbs3>::Instance();
    default:
      dserror("Unsupported target-type %s.", DRT::DistypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get3D(
    DRT::Element::DiscretizationType ref_type, DRT::Element::DiscretizationType tar_type,
    const bool debug)
{
  switch (ref_type)
  {
    case DRT::Element::quad4:
      return Get3D<DRT::Element::quad4>(tar_type, debug);
    case DRT::Element::tri3:
      return Get3D<DRT::Element::tri3>(tar_type, debug);
    case DRT::Element::nurbs4:
      return Get3D<DRT::Element::nurbs4>(tar_type, debug);
    case DRT::Element::nurbs9:
      return Get3D<DRT::Element::nurbs9>(tar_type, debug);
    default:
      dserror("Unsupported reference-type %s.", DRT::DistypeToString(ref_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType ref_type>
CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get3D(
    DRT::Element::DiscretizationType tar_type, const bool debug)
{
  switch (tar_type)
  {
    case DRT::Element::quad4:
      if (debug) return Projector<ProjDebugger, 3, ref_type, DRT::Element::quad4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, DRT::Element::quad4>::Instance();
    case DRT::Element::tri3:
      if (debug) return Projector<ProjDebugger, 3, ref_type, DRT::Element::tri3>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, DRT::Element::tri3>::Instance();
    case DRT::Element::nurbs4:
      if (debug) return Projector<ProjDebugger, 3, ref_type, DRT::Element::nurbs4>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, DRT::Element::nurbs4>::Instance();
    case DRT::Element::nurbs9:
      if (debug) return Projector<ProjDebugger, 3, ref_type, DRT::Element::nurbs9>::Instance();
      return Projector<EmptyProjDebugger, 3, ref_type, DRT::Element::nurbs9>::Instance();
    default:
      dserror("Unsupported target-type %s.", DRT::DistypeToString(tar_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
CONTACT::AUG::ProjectorBase*
CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::Instance(bool delete_me)
{
  static ProjectorBase* instance = NULL;

  if (delete_me)
  {
    if (instance) delete instance;

    return NULL;
  }

  if (not instance) instance = new Projector<DebugPolicy, probdim, ref_type, tar_type>;

  return instance;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
void CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::Done()
{
  Instance(true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
void CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::Setup()
{
  ref_val_.Scale(0.0);
  x_ref_.Scale(0.0);
  n_ref_.Scale(0.0);
  iter_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::operator()(
    MORTAR::MortarElement& ref_ele, const double* ref_xi, MORTAR::MortarElement& target_ele,
    double* target_xi, double& alpha)
{
  Setup();

  const DRT::Node* const* ref_nodes = ref_ele.Nodes();

  const LINALG::Matrix<REF_DIM, 1> rxi(ref_xi, true);

  shape_function<ref_type>(ref_ele, rxi, ref_val_);

  for (unsigned i = 0; i < REF_NUMNODES; ++i)
  {
    const MORTAR::MortarNode& mnode = static_cast<const MORTAR::MortarNode&>(*ref_nodes[i]);
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
    const MORTAR::MortarNode& mnode = static_cast<const MORTAR::MortarNode&>(*tnodes[i]);
    const double* tar_x = mnode.xspatial();

    std::copy(tar_x, tar_x + probdim, &tar_coords_(0, i));
  }

  // initial value
  LINALG::Matrix<TAR_DIM, 1> txi(target_xi, true);
  alpha = 0.0;

  LINALG::Matrix<TAR_DIM, 1> txi_center(false);
  DRT::UTILS::getLocalCenterPosition<TAR_DIM>(tar_type, txi_center);
  std::copy(txi_center.A(), txi_center.A() + TAR_DIM, txi.A());

  RhsGP(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha);
  DebugPolicy::writeVector(std::cout, probdim, rhs_.A(), "Rhs");

  const double ref_rhs_nrm2 = std::max(1.0, rhs_.Norm2());
  double dx_nrm2 = 0.0;
  double rhs_nrm2 = 0.0;
  bool is_parallel_proj = false;

  for (; iter_ < MORTARMAXITER; ++iter_)
  {
    LMatGP(lmat_, tar_deriv1_, target_ele, tar_coords_, target_xi, n_ref_);

    rhs_.Scale(-1.0);
    const double det = LINALG::gaussElimination<true, probdim>(lmat_, rhs_, dx_);

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
    if (not RhsGP(rhs_, x_ref_, n_ref_, target_ele, tar_coords_, target_xi, alpha))
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
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::RhsGP(
    LINALG::Matrix<probdim, 1>& rhs, const LINALG::Matrix<probdim, 1>& x_ref,
    const LINALG::Matrix<probdim, 1>& n_ref, MORTAR::MortarElement& target_ele,
    const LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
    const double& alpha) const
{
  LINALG::Matrix<probdim, 1> x_tar(false);
  const bool status = GetGlobalPosition<tar_type>(target_ele, tar_coords, tar_xi, x_tar);

  // evaluate right hand side
  for (unsigned i = 0; i < probdim; ++i)
    rhs(i, 0) = x_tar(i, 0) - alpha * n_ref(i, 0) - x_ref(i, 0);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
template <DRT::Element::DiscretizationType type, unsigned numnodes>
bool CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::GetGlobalPosition(
    MORTAR::MortarElement& ele, const LINALG::Matrix<probdim, numnodes>& coords, const double* xi,
    LINALG::Matrix<probdim, 1>& pos) const
{
  const unsigned dim = DRT::UTILS::DisTypeToDim<type>::dim;
  const LINALG::Matrix<dim, 1> mat_xi(xi, true);

  LINALG::Matrix<numnodes, 1> val(true);
  const bool status = shape_function<type>(ele, mat_xi, val);

  pos.Multiply(coords, val);
  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class DebugPolicy, unsigned probdim, DRT::Element::DiscretizationType ref_type,
    DRT::Element::DiscretizationType tar_type>
void CONTACT::AUG::Projector<DebugPolicy, probdim, ref_type, tar_type>::LMatGP(
    LINALG::Matrix<probdim, probdim>& lmat, LINALG::Matrix<TAR_DIM, TAR_NUMNODES>& tar_deriv1,
    MORTAR::MortarElement& tar_ele, const LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords,
    const double* tar_xi, const LINALG::Matrix<probdim, 1>& n_ref) const
{
  const LINALG::Matrix<TAR_DIM, 1> txi(tar_xi, true);
  std::fill(lmat.A(), lmat.A() + 2 * probdim, 0.0);

  shape_function_deriv1<tar_type>(tar_ele, txi, tar_deriv1);

  for (unsigned i = 0; i < TAR_NUMNODES; ++i)
    for (unsigned j = 0; j < TAR_DIM; ++j)
      for (unsigned k = 0; k < probdim; ++k) lmat(k, j) += tar_deriv1(j, i) * tar_coords(k, i);

  for (unsigned i = 0; i < probdim; ++i) lmat(i, probdim - 1) = -n_ref(i, 0);
}


template CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get2D<DRT::Element::line2>(
    DRT::Element::DiscretizationType tar_type, const bool debug);
template CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get3D<DRT::Element::quad4>(
    DRT::Element::DiscretizationType tar_type, const bool debug);
template CONTACT::AUG::ProjectorBase* CONTACT::AUG::ProjectorBase::Get3D<DRT::Element::tri3>(
    DRT::Element::DiscretizationType tar_type, const bool debug);

// standard discretization types
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 2, DRT::Element::line2,
    DRT::Element::line2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 2, DRT::Element::line2,
    DRT::Element::line2>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::quad4,
    DRT::Element::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::quad4,
    DRT::Element::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::quad4,
    DRT::Element::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::quad4,
    DRT::Element::tri3>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::tri3,
    DRT::Element::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::tri3,
    DRT::Element::tri3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::tri3,
    DRT::Element::quad4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::tri3,
    DRT::Element::quad4>;

// pure NURBS discretization types
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs2,
    DRT::Element::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs2,
    DRT::Element::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs3,
    DRT::Element::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs3,
    DRT::Element::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs2,
    DRT::Element::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs2,
    DRT::Element::nurbs3>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs3,
    DRT::Element::nurbs2>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs3,
    DRT::Element::nurbs2>;

template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs4,
    DRT::Element::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs4,
    DRT::Element::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs9,
    DRT::Element::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs9,
    DRT::Element::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs4,
    DRT::Element::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs4,
    DRT::Element::nurbs9>;
template class CONTACT::AUG::Projector<CONTACT::AUG::ProjDebugger, 3, DRT::Element::nurbs9,
    DRT::Element::nurbs4>;
template class CONTACT::AUG::Projector<CONTACT::AUG::EmptyProjDebugger, 3, DRT::Element::nurbs9,
    DRT::Element::nurbs4>;
