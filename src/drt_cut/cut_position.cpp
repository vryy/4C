/*----------------------------------------------------------------------------*/
/*!
\file cut_position.cpp

\brief Class to check whether a point lies inside an element, and also to
       compute local coordinates of a point with respect to the element.
       The standard as well as the embedded cases are treated here. Feel
       free to extend the content for your needs.

\level 2

\maintainer Michael Hiermeier
*/
/*----------------------------------------------------------------------------*/

#include "cut_position.H"
#include "cut_boundingbox.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create(
    const Element& element, const Point& point)
{
  const PositionFactory factory;
  return factory.CreatePosition(element, point);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create(
    const Element& element, const LINALG::Matrix<rdim, 1>& xyz)
{
  const PositionFactory factory;
  if (rdim < factory.ProbDim())
    dserror(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.ProbDim());
  return factory.CreatePosition(element, xyz.A());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim, unsigned cdim, unsigned rdim_2>
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create(const LINALG::Matrix<rdim, cdim>& xyze,
    const LINALG::Matrix<rdim_2, 1>& xyz, const DRT::Element::DiscretizationType& distype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.ProbDim();
  const unsigned num_nodes_ele = DRT::UTILS::getNumberOfElementNodes(distype);

  if (rdim < probdim or cdim != num_nodes_ele)
    dserror(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.M(), xyze.N());

  const double* xyze_ptr = xyze.A();
  Epetra_SerialDenseMatrix xyze_eptra;
  if (rdim > probdim)
  {
    xyze_eptra.Shape(probdim, num_nodes_ele);
    FixMatrixShape(xyze, xyze_eptra);
    xyze_ptr = xyze_eptra.A();
  }

  if (rdim_2 < probdim)
    dserror(
        "Dimension mismatch of xyz! \n"
        "expected input: %d x 1 (rows x cols)\n"
        "received input : %d x 1 (rows x cols)",
        probdim, rdim_2);

  return factory.CreatePosition(xyze_ptr, xyz.A(), distype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create(const Epetra_SerialDenseMatrix& xyze,
    const LINALG::Matrix<rdim, 1>& xyz, const DRT::Element::DiscretizationType& distype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.ProbDim();
  const unsigned num_nodes_ele = DRT::UTILS::getNumberOfElementNodes(distype);

  if (static_cast<unsigned>(xyze.M()) < probdim or static_cast<unsigned>(xyze.N()) != num_nodes_ele)
    dserror(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.M(), xyze.N());

  const double* xyze_ptr = xyze.A();
  Epetra_SerialDenseMatrix xyze_eptra;
  if (static_cast<unsigned>(xyze.M()) > probdim)
  {
    xyze_eptra.Shape(probdim, num_nodes_ele);
    FixMatrixShape(xyze, xyze_eptra);
    xyze_ptr = xyze_eptra.A();
  }

  if (xyz.M() < probdim)
    dserror(
        "Dimension mismatch of xyz! \n"
        "expected input: %d x 1 (rows x cols)\n"
        "received input : %d x 1 (rows x cols)",
        probdim, xyz.M());

  return factory.CreatePosition(xyze_ptr, xyz.A(), distype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create(const std::vector<Node*> nodes,
    const LINALG::Matrix<rdim, 1>& xyz, DRT::Element::DiscretizationType distype)
{
  const PositionFactory factory;
  if (rdim < factory.ProbDim())
    dserror(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.ProbDim());
  return factory.CreatePosition(nodes, xyz.A(), distype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned numNodesElement,
    unsigned dim>
void GEO::CUT::PositionGeneric<probdim, eletype, numNodesElement, dim>::ConstructBoundingBox()
{
  bbside_ = Teuchos::rcp(BoundingBox::Create());

  for (unsigned i = 0; i < numNodesElement; ++i)
  {
    LINALG::Matrix<3, 1> x1(&this->xyze_(0, i), true);
    bbside_->AddPoint(x1);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned numNodesElement,
    unsigned dim>
bool GEO::CUT::ComputePosition<probdim, eletype, numNodesElement, dim>::Compute(const double& Tol)
{
  /* If the given point is outside the element bounding box, no need
   * to perform the complex calculations */
  //  if( not this->bbside_->Within( 1.0, &this->px_(0,0) ) )
  //  {
  //    this->pos_status_ = Position::position_outside_of_bbox;
  //    return false;
  //  }

  //      KERNEL::DebugComputePosition<probdim,eletype> cp( this->xsi_ );
  KERNEL::ComputePosition<probdim, eletype> cp(this->xsi_);
  this->pos_status_ =
      (cp(this->xyze_, this->px_) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cp.GetTolerance();

  if (this->pos_status_ == Position::position_valid) return WithinLimitsTol(Tol);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned numNodesElement,
    unsigned dim>
bool GEO::CUT::ComputeEmbeddedPosition<probdim, eletype, numNodesElement,
    dim>::IsGivenPointWithinElement()
{
  // If the given point is outside the side's bounding box, no need to perform
  // the complex calculations
  if (not this->bbside_->Within(1.0, &this->px_(0, 0)))
  {
    this->pos_status_ = Position::position_outside_of_bbox;
    return false;
  }

  // try to compute the local coordinates and the distance using the Newton scheme
  KERNEL::ComputeDistance<probdim, eletype> cd(xsi_aug_);
  // KERNEL::DebugComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );

  double dist = 0.0;
  this->pos_status_ =
      (cd(this->xyze_, this->px_, dist) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cd.GetTolerance();

  return (this->pos_status_ == Position::position_valid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned numNodesElement,
    unsigned dim>
bool GEO::CUT::ComputeEmbeddedPosition<probdim, eletype, numNodesElement, dim>::Compute(
    const double& Tol, const bool& allow_dist)
{
  xsi_aug_ = 0.0;
  /* If the given point is outside the side's bounding box, no need
   * to perform the complex calculations */
  if (not allow_dist)
  {
    if (not this->bbside_->Within(1.0, &this->px_(0, 0)))
    {
      this->pos_status_ = Position::position_outside_of_bbox;
      return false;
    }
  }

  // try to compute the local coordinates and the distance using the Newton scheme
  KERNEL::ComputeDistance<probdim, eletype> cd(xsi_aug_);
  //  KERNEL::DebugComputeDistance< probdim, eletype > cd( xsi_aug_ );

  double dist = 0.0;
  this->pos_status_ =
      (cd(this->xyze_, this->px_, dist) ? Position::position_valid : Position::position_invalid);

  this->compute_tolerance_ = cd.GetTolerance();

  // if newton did not converge, try some other checks
  if (this->pos_status_ != Position::position_valid)
  {
    // if the surface has no area zero area
    if (cd.ZeroArea())
    {
      // Comment: If there is no area an the surface is a line it would be possible to calculate the
      // distance between the lines,
      // but we would still loose the sign of the distance (if we need this we would still need to
      // add something like Position::position_absdistance_valid) - for this case the normal is just
      // not defined anymore
      this->pos_status_ = Position::position_zero_area;
      return false;
    }
    else  // this branch is a little bit unsafe --> write a warning when this case is used
    {
      // Check the l2-norm of the actual reached residual. If it is smaller than 1.0,
      // we will trust the distance value
      // --> this is quite critical as we cannot guarantee that the precicion of the distance
      // anymore (but only option to get a signed distance)
      double res_nrm2 = cd.GetResidualL2Norm();
      if (res_nrm2 < 1.0)
      {
        std::cout << "==| WARNING: ComputeDistance didn't converge, but residual is smaller than 1 "
                     "(calculated distance is used) |=="
                  << std::endl;
        this->pos_status_ = Position::position_distance_valid;
      }

      // as the local coordinates didn't converge
      return false;
    }
  }

  // Std Return in case Compute Distance converged!
  return WithinLimitsTol(Tol, allow_dist);
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::PositionFactory::PositionFactory() : probdim_(DRT::Problem::Instance()->NDim()) {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::PositionFactory::CreatePosition(
    const Element& element, const Point& point) const
{
  DRT::Element::DiscretizationType distype = element.Shape();

  switch (distype)
  {
    case DRT::Element::line2:
      return CreateConcretePosition<DRT::Element::line2>(element, point);
    case DRT::Element::tri3:
      return CreateConcretePosition<DRT::Element::tri3>(element, point);
    case DRT::Element::tri6:
      return CreateConcretePosition<DRT::Element::tri6>(element, point);
    case DRT::Element::quad4:
      return CreateConcretePosition<DRT::Element::quad4>(element, point);
    case DRT::Element::quad8:
      return CreateConcretePosition<DRT::Element::quad8>(element, point);
    case DRT::Element::quad9:
      return CreateConcretePosition<DRT::Element::quad9>(element, point);
    case DRT::Element::hex8:
      return CreateConcretePosition<DRT::Element::hex8>(element, point);
    case DRT::Element::hex20:
      return CreateConcretePosition<DRT::Element::hex20>(element, point);
    case DRT::Element::tet4:
      return CreateConcretePosition<DRT::Element::tet4>(element, point);
    case DRT::Element::pyramid5:
      return CreateConcretePosition<DRT::Element::pyramid5>(element, point);
    case DRT::Element::wedge6:
      return CreateConcretePosition<DRT::Element::wedge6>(element, point);
    default:
      dserror("Unsupported distype = %s", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::PositionFactory::CreatePosition(
    const Element& element, const double* xyz) const
{
  DRT::Element::DiscretizationType distype = element.Shape();

  switch (distype)
  {
    case DRT::Element::line2:
      return CreateConcretePosition<DRT::Element::line2>(element, xyz);
    case DRT::Element::tri3:
      return CreateConcretePosition<DRT::Element::tri3>(element, xyz);
    case DRT::Element::tri6:
      return CreateConcretePosition<DRT::Element::tri6>(element, xyz);
    case DRT::Element::quad4:
      return CreateConcretePosition<DRT::Element::quad4>(element, xyz);
    case DRT::Element::quad8:
      return CreateConcretePosition<DRT::Element::quad8>(element, xyz);
    case DRT::Element::quad9:
      return CreateConcretePosition<DRT::Element::quad9>(element, xyz);
    case DRT::Element::hex8:
      return CreateConcretePosition<DRT::Element::hex8>(element, xyz);
    case DRT::Element::hex20:
      return CreateConcretePosition<DRT::Element::hex20>(element, xyz);
    case DRT::Element::tet4:
      return CreateConcretePosition<DRT::Element::tet4>(element, xyz);
    case DRT::Element::pyramid5:
      return CreateConcretePosition<DRT::Element::pyramid5>(element, xyz);
    case DRT::Element::wedge6:
      return CreateConcretePosition<DRT::Element::wedge6>(element, xyz);
    default:
      dserror("Unsupported distype = %s", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::PositionFactory::CreatePosition(
    const double* xyze, const double* xyz, const DRT::Element::DiscretizationType& distype) const
{
  switch (distype)
  {
    case DRT::Element::line2:
      return CreateConcretePosition<DRT::Element::line2>(xyze, xyz);
    case DRT::Element::tri3:
      return CreateConcretePosition<DRT::Element::tri3>(xyze, xyz);
    case DRT::Element::tri6:
      return CreateConcretePosition<DRT::Element::tri6>(xyze, xyz);
    case DRT::Element::quad4:
      return CreateConcretePosition<DRT::Element::quad4>(xyze, xyz);
    case DRT::Element::quad8:
      return CreateConcretePosition<DRT::Element::quad8>(xyze, xyz);
    case DRT::Element::quad9:
      return CreateConcretePosition<DRT::Element::quad9>(xyze, xyz);
    case DRT::Element::hex8:
      return CreateConcretePosition<DRT::Element::hex8>(xyze, xyz);
    case DRT::Element::hex20:
      return CreateConcretePosition<DRT::Element::hex20>(xyze, xyz);
    case DRT::Element::tet4:
      return CreateConcretePosition<DRT::Element::tet4>(xyze, xyz);
    case DRT::Element::pyramid5:
      return CreateConcretePosition<DRT::Element::pyramid5>(xyze, xyz);
    case DRT::Element::wedge6:
      return CreateConcretePosition<DRT::Element::wedge6>(xyze, xyz);
    default:
      dserror("Unsupported distype = %s", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Position> GEO::CUT::PositionFactory::CreatePosition(
    const std::vector<GEO::CUT::Node*> nodes, const double* xyz,
    DRT::Element::DiscretizationType distype) const
{
  if (distype == DRT::Element::dis_none)
  {
    plain_element_set elements;
    FindCommonElements(nodes, elements);
    if (elements.size() != 1)
      dserror("Couldn't find a unique element corresponding to the given nodes.");

    distype = elements[0]->Shape();
  }

  switch (distype)
  {
    case DRT::Element::line2:
      return CreateConcretePosition<DRT::Element::line2>(nodes, xyz);
    case DRT::Element::tri3:
      return CreateConcretePosition<DRT::Element::tri3>(nodes, xyz);
    case DRT::Element::tri6:
      return CreateConcretePosition<DRT::Element::tri6>(nodes, xyz);
    case DRT::Element::quad4:
      return CreateConcretePosition<DRT::Element::quad4>(nodes, xyz);
    case DRT::Element::quad8:
      return CreateConcretePosition<DRT::Element::quad8>(nodes, xyz);
    case DRT::Element::quad9:
      return CreateConcretePosition<DRT::Element::quad9>(nodes, xyz);
    case DRT::Element::hex8:
      return CreateConcretePosition<DRT::Element::hex8>(nodes, xyz);
    case DRT::Element::hex20:
      return CreateConcretePosition<DRT::Element::hex20>(nodes, xyz);
    case DRT::Element::tet4:
      return CreateConcretePosition<DRT::Element::tet4>(nodes, xyz);
    case DRT::Element::pyramid5:
      return CreateConcretePosition<DRT::Element::pyramid5>(nodes, xyz);
    case DRT::Element::wedge6:
      return CreateConcretePosition<DRT::Element::wedge6>(nodes, xyz);
    default:
      dserror("Unsupported distype = %s", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}


template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<2>(
    const Element& element, const LINALG::Matrix<2, 1>& xyz);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3>(
    const Element& element, const LINALG::Matrix<3, 1>& xyz);

template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 3, 3>(
    const LINALG::Matrix<3, 3>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 6, 3>(
    const LINALG::Matrix<3, 6>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 4, 3>(
    const LINALG::Matrix<3, 4>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 8, 3>(
    const LINALG::Matrix<3, 8>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 9, 3>(
    const LINALG::Matrix<3, 9>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3, 2, 3>(
    const LINALG::Matrix<3, 2>& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<2, 2, 2>(
    const LINALG::Matrix<2, 2>& xyze, const LINALG::Matrix<2, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);

template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3>(
    const Epetra_SerialDenseMatrix& xyze, const LINALG::Matrix<3, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<2>(
    const Epetra_SerialDenseMatrix& xyze, const LINALG::Matrix<2, 1>& xyz,
    const DRT::Element::DiscretizationType& distype);


template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<3>(
    const std::vector<Node*> nodes, const LINALG::Matrix<3, 1>& xyz,
    DRT::Element::DiscretizationType distype);
template Teuchos::RCP<GEO::CUT::Position> GEO::CUT::Position::Create<2>(
    const std::vector<Node*> nodes, const LINALG::Matrix<2, 1>& xyz,
    DRT::Element::DiscretizationType distype);

/* --- ComputeEmbeddedPosition --- */
// embedded element types
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::tri3>;
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::tri6>;
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::quad4>;
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::quad8>;
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::quad9>;
template class GEO::CUT::ComputeEmbeddedPosition<2, DRT::Element::line2>;
template class GEO::CUT::ComputeEmbeddedPosition<3, DRT::Element::line2>;

// non-embedded element types for the embedded case (only necessary due to compiler problems)
// template class GEO::CUT::ComputeEmbeddedPosition<3,DRT::Element::hex16>;
// template class GEO::CUT::ComputeEmbeddedPosition<3,DRT::Element::hex18>;
// template class GEO::CUT::ComputeEmbeddedPosition<3,DRT::Element::hex27>;
// template class GEO::CUT::ComputeEmbeddedPosition<3,DRT::Element::tet10>;
// template class GEO::CUT::ComputeEmbeddedPosition<3,DRT::Element::wedge15>;

/* --- ComputePosition --- */
// non-embedded cases (only)
template class GEO::CUT::ComputePosition<1, DRT::Element::line2>;

template class GEO::CUT::ComputePosition<2, DRT::Element::tri3>;
template class GEO::CUT::ComputePosition<2, DRT::Element::tri6>;
template class GEO::CUT::ComputePosition<2, DRT::Element::quad4>;
template class GEO::CUT::ComputePosition<2, DRT::Element::quad8>;
template class GEO::CUT::ComputePosition<2, DRT::Element::quad9>;

template class GEO::CUT::ComputePosition<3, DRT::Element::tet4>;
template class GEO::CUT::ComputePosition<3, DRT::Element::tet10>;
template class GEO::CUT::ComputePosition<3, DRT::Element::hex8>;
template class GEO::CUT::ComputePosition<3, DRT::Element::hex16>;
template class GEO::CUT::ComputePosition<3, DRT::Element::hex18>;
template class GEO::CUT::ComputePosition<3, DRT::Element::hex20>;
template class GEO::CUT::ComputePosition<3, DRT::Element::hex27>;
template class GEO::CUT::ComputePosition<3, DRT::Element::pyramid5>;
template class GEO::CUT::ComputePosition<3, DRT::Element::wedge6>;
template class GEO::CUT::ComputePosition<3, DRT::Element::wedge15>;

/* --- PositionGeneric --- */
// embedded cases
template class GEO::CUT::PositionGeneric<3, DRT::Element::tri3>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::tri6>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::quad4>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::quad8>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::quad9>;
template class GEO::CUT::PositionGeneric<2, DRT::Element::line2>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::line2>;

// non-embedded cases
template class GEO::CUT::PositionGeneric<1, DRT::Element::line2>;

template class GEO::CUT::PositionGeneric<2, DRT::Element::tri3>;
template class GEO::CUT::PositionGeneric<2, DRT::Element::tri6>;
template class GEO::CUT::PositionGeneric<2, DRT::Element::quad4>;
template class GEO::CUT::PositionGeneric<2, DRT::Element::quad8>;
template class GEO::CUT::PositionGeneric<2, DRT::Element::quad9>;

template class GEO::CUT::PositionGeneric<3, DRT::Element::tet4>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::tet10>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::hex8>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::hex16>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::hex18>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::hex20>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::hex27>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::pyramid5>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::wedge6>;
template class GEO::CUT::PositionGeneric<3, DRT::Element::wedge15>;

// unused / impossible cases (only necessary due to compiler problems)
// template class GEO::CUT::PositionGeneric<2,DRT::Element::hex8>;
// template class GEO::CUT::PositionGeneric<2,DRT::Element::tet4>;
// template class GEO::CUT::PositionGeneric<2,DRT::Element::pyramid5>;
// template class GEO::CUT::PositionGeneric<2,DRT::Element::wedge6>;
