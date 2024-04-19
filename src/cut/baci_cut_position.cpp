/*----------------------------------------------------------------------------*/
/*! \file

\brief Class to check whether a point lies inside an element, and also to
       compute local coordinates of a point with respect to the element.
       The standard as well as the embedded cases are treated here. Feel
       free to extend the content for your needs.

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_cut_position.hpp"

#include "baci_cut_boundingbox.hpp"
#include "baci_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create(
    const Element& element, const Point& point, INPAR::CUT::CutFloattype floattype)
{
  const PositionFactory factory;
  return factory.CreatePosition(element, point, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create(const Element& element,
    const CORE::LINALG::Matrix<rdim, 1>& xyz, INPAR::CUT::CutFloattype floattype)
{
  const PositionFactory factory;
  if (rdim < factory.ProbDim())
    FOUR_C_THROW(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.ProbDim());
  return factory.CreatePosition(element, xyz.A(), floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim, unsigned cdim, unsigned rdim_2>
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create(
    const CORE::LINALG::Matrix<rdim, cdim>& xyze, const CORE::LINALG::Matrix<rdim_2, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.ProbDim();
  const unsigned num_nodes_ele = CORE::FE::getNumberOfElementNodes(distype);

  if (rdim < probdim or cdim != num_nodes_ele)
    FOUR_C_THROW(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.M(), xyze.N());

  const double* xyze_ptr = xyze.A();
  CORE::LINALG::SerialDenseMatrix xyze_eptra;
  if (rdim > probdim)
  {
    xyze_eptra.shape(probdim, num_nodes_ele);
    FixMatrixShape(xyze, xyze_eptra);
    xyze_ptr = xyze_eptra.values();
  }

  if (rdim_2 < probdim)
    FOUR_C_THROW(
        "Dimension mismatch of xyz! \n"
        "expected input: %d x 1 (rows x cols)\n"
        "received input : %d x 1 (rows x cols)",
        probdim, rdim_2);

  return factory.CreatePosition(xyze_ptr, xyz.A(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create(
    const CORE::LINALG::SerialDenseMatrix& xyze, const CORE::LINALG::Matrix<rdim, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.ProbDim();
  const unsigned num_nodes_ele = CORE::FE::getNumberOfElementNodes(distype);

  if (static_cast<unsigned>(xyze.numRows()) < probdim or
      static_cast<unsigned>(xyze.numCols()) != num_nodes_ele)
    FOUR_C_THROW(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.numRows(), xyze.numCols());

  const double* xyze_ptr = xyze.values();
  CORE::LINALG::SerialDenseMatrix xyze_eptra;
  if (static_cast<unsigned>(xyze.numRows()) > probdim)
  {
    xyze_eptra.shape(probdim, num_nodes_ele);
    FixMatrixShape(xyze, xyze_eptra);
    xyze_ptr = xyze_eptra.values();
  }

  if (xyz.numRows() < probdim)
    FOUR_C_THROW(
        "Dimension mismatch of xyz! \n"
        "expected input: %d x 1 (rows x cols)\n"
        "received input : %d x 1 (rows x cols)",
        probdim, xyz.numRows());

  return factory.CreatePosition(xyze_ptr, xyz.A(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create(
    const std::vector<Node*> nodes, const CORE::LINALG::Matrix<rdim, 1>& xyz,
    CORE::FE::CellType distype, INPAR::CUT::CutFloattype floattype)
{
  const PositionFactory factory;
  if (rdim < factory.ProbDim())
    FOUR_C_THROW(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.ProbDim());
  return factory.CreatePosition(nodes, xyz.A(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType eletype, unsigned numNodesElement, unsigned dim,
    INPAR::CUT::CutFloattype floattype>
void CORE::GEO::CUT::PositionGeneric<probdim, eletype, numNodesElement, dim,
    floattype>::ConstructBoundingBox()
{
  bbside_ = Teuchos::rcp(BoundingBox::Create());

  for (unsigned i = 0; i < numNodesElement; ++i)
  {
    CORE::LINALG::Matrix<3, 1> x1(&this->xyze_(0, i), true);
    bbside_->AddPoint(x1);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType eletype, unsigned numNodesElement, unsigned dim,
    INPAR::CUT::CutFloattype floattype>
bool CORE::GEO::CUT::ComputePosition<probdim, eletype, numNodesElement, dim, floattype>::Compute(
    const double& Tol)
{
  /* If the given point is outside the element bounding box, no need
   * to perform the complex calculations */
  //  if( not this->bbside_->Within( 1.0, &this->px_(0,0) ) )
  //  {
  //    this->pos_status_ = Position::position_outside_of_bbox;
  //    return false;
  //  }

  // no cln used for compute position
  KERNEL::ComputePosition<probdim, eletype, floattype> cp(this->xsi_);
  //  KERNEL::DebugComputePosition<probdim,eletype,floattype> cp( this->xsi_ );
  this->pos_status_ =
      (cp(this->xyze_, this->px_) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cp.GetTolerance();

  if (this->pos_status_ == Position::position_valid) return WithinLimitsTol(Tol);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType eletype, unsigned numNodesElement, unsigned dim,
    INPAR::CUT::CutFloattype floattype>
bool CORE::GEO::CUT::ComputeEmbeddedPosition<probdim, eletype, numNodesElement, dim,
    floattype>::IsGivenPointWithinElement()
{
  // If the given point is outside the side's bounding box, no need to perform
  // the complex calculations
  if (not this->bbside_->Within(1.0, &this->px_(0, 0)))
  {
    this->pos_status_ = Position::position_outside_of_bbox;
    return false;
  }

  // try to compute the local coordinates and the distance using the Newton scheme
  KERNEL::ComputeDistance<probdim, eletype, (floattype == INPAR::CUT::floattype_cln)> cd(xsi_aug_);
  //  KERNEL::DebugComputeDistance< probdim, eletype,(floattype==INPAR::CUT::floattype_cln) > cd(
  //  xsi_aug_ );

  double dist = 0.0;
  this->pos_status_ =
      (cd(this->xyze_, this->px_, dist) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cd.GetTolerance();

  return (this->pos_status_ == Position::position_valid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType eletype, unsigned numNodesElement, unsigned dim,
    INPAR::CUT::CutFloattype floattype>
bool CORE::GEO::CUT::ComputeEmbeddedPosition<probdim, eletype, numNodesElement, dim,
    floattype>::Compute(const double& Tol, const bool& allow_dist)
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
  KERNEL::ComputeDistance<probdim, eletype, (floattype == INPAR::CUT::floattype_cln)> cd(xsi_aug_);
  //  KERNEL::DebugComputeDistance< probdim, eletype,(floattype==INPAR::CUT::floattype_cln) > cd(
  //  xsi_aug_ );

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
CORE::GEO::CUT::PositionFactory::PositionFactory() : probdim_(GLOBAL::Problem::Instance()->NDim())
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::PositionFactory::CreatePosition(
    const Element& element, const Point& point, INPAR::CUT::CutFloattype floattype) const
{
  CORE::FE::CellType distype = element.Shape();

  switch (distype)
  {
    case CORE::FE::CellType::line2:
      return CreateConcretePosition<CORE::FE::CellType::line2>(element, point, floattype);
    case CORE::FE::CellType::tri3:
      return CreateConcretePosition<CORE::FE::CellType::tri3>(element, point, floattype);
    case CORE::FE::CellType::tri6:
      return CreateConcretePosition<CORE::FE::CellType::tri6>(element, point, floattype);
    case CORE::FE::CellType::quad4:
      return CreateConcretePosition<CORE::FE::CellType::quad4>(element, point, floattype);
    case CORE::FE::CellType::quad8:
      return CreateConcretePosition<CORE::FE::CellType::quad8>(element, point, floattype);
    case CORE::FE::CellType::quad9:
      return CreateConcretePosition<CORE::FE::CellType::quad9>(element, point, floattype);
    case CORE::FE::CellType::hex8:
      return CreateConcretePosition<CORE::FE::CellType::hex8>(element, point, floattype);
    case CORE::FE::CellType::hex20:
      return CreateConcretePosition<CORE::FE::CellType::hex20>(element, point, floattype);
    case CORE::FE::CellType::tet4:
      return CreateConcretePosition<CORE::FE::CellType::tet4>(element, point, floattype);
    case CORE::FE::CellType::pyramid5:
      return CreateConcretePosition<CORE::FE::CellType::pyramid5>(element, point, floattype);
    case CORE::FE::CellType::wedge6:
      return CreateConcretePosition<CORE::FE::CellType::wedge6>(element, point, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", CORE::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::PositionFactory::CreatePosition(
    const Element& element, const double* xyz, INPAR::CUT::CutFloattype floattype) const
{
  CORE::FE::CellType distype = element.Shape();

  switch (distype)
  {
    case CORE::FE::CellType::line2:
      return CreateConcretePosition<CORE::FE::CellType::line2>(element, xyz, floattype);
    case CORE::FE::CellType::tri3:
      return CreateConcretePosition<CORE::FE::CellType::tri3>(element, xyz, floattype);
    case CORE::FE::CellType::tri6:
      return CreateConcretePosition<CORE::FE::CellType::tri6>(element, xyz, floattype);
    case CORE::FE::CellType::quad4:
      return CreateConcretePosition<CORE::FE::CellType::quad4>(element, xyz, floattype);
    case CORE::FE::CellType::quad8:
      return CreateConcretePosition<CORE::FE::CellType::quad8>(element, xyz, floattype);
    case CORE::FE::CellType::quad9:
      return CreateConcretePosition<CORE::FE::CellType::quad9>(element, xyz, floattype);
    case CORE::FE::CellType::hex8:
      return CreateConcretePosition<CORE::FE::CellType::hex8>(element, xyz, floattype);
    case CORE::FE::CellType::hex20:
      return CreateConcretePosition<CORE::FE::CellType::hex20>(element, xyz, floattype);
    case CORE::FE::CellType::tet4:
      return CreateConcretePosition<CORE::FE::CellType::tet4>(element, xyz, floattype);
    case CORE::FE::CellType::pyramid5:
      return CreateConcretePosition<CORE::FE::CellType::pyramid5>(element, xyz, floattype);
    case CORE::FE::CellType::wedge6:
      return CreateConcretePosition<CORE::FE::CellType::wedge6>(element, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", CORE::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::PositionFactory::CreatePosition(
    const double* xyze, const double* xyz, const CORE::FE::CellType& distype,
    INPAR::CUT::CutFloattype floattype) const
{
  switch (distype)
  {
    case CORE::FE::CellType::line2:
      return CreateConcretePosition<CORE::FE::CellType::line2>(xyze, xyz, floattype);
    case CORE::FE::CellType::tri3:
      return CreateConcretePosition<CORE::FE::CellType::tri3>(xyze, xyz, floattype);
    case CORE::FE::CellType::tri6:
      return CreateConcretePosition<CORE::FE::CellType::tri6>(xyze, xyz, floattype);
    case CORE::FE::CellType::quad4:
      return CreateConcretePosition<CORE::FE::CellType::quad4>(xyze, xyz, floattype);
    case CORE::FE::CellType::quad8:
      return CreateConcretePosition<CORE::FE::CellType::quad8>(xyze, xyz, floattype);
    case CORE::FE::CellType::quad9:
      return CreateConcretePosition<CORE::FE::CellType::quad9>(xyze, xyz, floattype);
    case CORE::FE::CellType::hex8:
      return CreateConcretePosition<CORE::FE::CellType::hex8>(xyze, xyz, floattype);
    case CORE::FE::CellType::hex20:
      return CreateConcretePosition<CORE::FE::CellType::hex20>(xyze, xyz, floattype);
    case CORE::FE::CellType::tet4:
      return CreateConcretePosition<CORE::FE::CellType::tet4>(xyze, xyz, floattype);
    case CORE::FE::CellType::pyramid5:
      return CreateConcretePosition<CORE::FE::CellType::pyramid5>(xyze, xyz, floattype);
    case CORE::FE::CellType::wedge6:
      return CreateConcretePosition<CORE::FE::CellType::wedge6>(xyze, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", CORE::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::PositionFactory::CreatePosition(
    const std::vector<CORE::GEO::CUT::Node*> nodes, const double* xyz, CORE::FE::CellType distype,
    INPAR::CUT::CutFloattype floattype) const
{
  if (distype == CORE::FE::CellType::dis_none)
  {
    plain_element_set elements;
    FindCommonElements(nodes, elements);
    if (elements.size() != 1)
      FOUR_C_THROW("Couldn't find a unique element corresponding to the given nodes.");

    distype = elements[0]->Shape();
  }

  switch (distype)
  {
    case CORE::FE::CellType::line2:
      return CreateConcretePosition<CORE::FE::CellType::line2>(nodes, xyz, floattype);
    case CORE::FE::CellType::tri3:
      return CreateConcretePosition<CORE::FE::CellType::tri3>(nodes, xyz, floattype);
    case CORE::FE::CellType::tri6:
      return CreateConcretePosition<CORE::FE::CellType::tri6>(nodes, xyz, floattype);
    case CORE::FE::CellType::quad4:
      return CreateConcretePosition<CORE::FE::CellType::quad4>(nodes, xyz, floattype);
    case CORE::FE::CellType::quad8:
      return CreateConcretePosition<CORE::FE::CellType::quad8>(nodes, xyz, floattype);
    case CORE::FE::CellType::quad9:
      return CreateConcretePosition<CORE::FE::CellType::quad9>(nodes, xyz, floattype);
    case CORE::FE::CellType::hex8:
      return CreateConcretePosition<CORE::FE::CellType::hex8>(nodes, xyz, floattype);
    case CORE::FE::CellType::hex20:
      return CreateConcretePosition<CORE::FE::CellType::hex20>(nodes, xyz, floattype);
    case CORE::FE::CellType::tet4:
      return CreateConcretePosition<CORE::FE::CellType::tet4>(nodes, xyz, floattype);
    case CORE::FE::CellType::pyramid5:
      return CreateConcretePosition<CORE::FE::CellType::pyramid5>(nodes, xyz, floattype);
    case CORE::FE::CellType::wedge6:
      return CreateConcretePosition<CORE::FE::CellType::wedge6>(nodes, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", CORE::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CUT::CutFloattype CORE::GEO::CUT::PositionFactory::UsePosFloattype(
    INPAR::CUT::CutFloattype floattype)
{
  if (general_pos_floattype_ != INPAR::CUT::floattype_none)
    return general_pos_floattype_;
  else
    return floattype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CUT::CutFloattype CORE::GEO::CUT::PositionFactory::UseDistFloattype(
    INPAR::CUT::CutFloattype floattype)
{
  if (general_dist_floattype_ != INPAR::CUT::floattype_none)
    return general_dist_floattype_;
  else
    return floattype;
}

INPAR::CUT::CutFloattype CORE::GEO::CUT::PositionFactory::general_pos_floattype_ =
    INPAR::CUT::floattype_none;
INPAR::CUT::CutFloattype CORE::GEO::CUT::PositionFactory::general_dist_floattype_ =
    INPAR::CUT::floattype_none;

template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<2>(
    const Element& element, const CORE::LINALG::Matrix<2, 1>& xyz,
    INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3>(
    const Element& element, const CORE::LINALG::Matrix<3, 1>& xyz,
    INPAR::CUT::CutFloattype floattype);

template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 3, 3>(
    const CORE::LINALG::Matrix<3, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 6, 3>(
    const CORE::LINALG::Matrix<3, 6>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 4, 3>(
    const CORE::LINALG::Matrix<3, 4>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 8, 3>(
    const CORE::LINALG::Matrix<3, 8>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 9, 3>(
    const CORE::LINALG::Matrix<3, 9>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3, 2, 3>(
    const CORE::LINALG::Matrix<3, 2>& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<2, 2, 2>(
    const CORE::LINALG::Matrix<2, 2>& xyze, const CORE::LINALG::Matrix<2, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);

template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3>(
    const CORE::LINALG::SerialDenseMatrix& xyze, const CORE::LINALG::Matrix<3, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<2>(
    const CORE::LINALG::SerialDenseMatrix& xyze, const CORE::LINALG::Matrix<2, 1>& xyz,
    const CORE::FE::CellType& distype, INPAR::CUT::CutFloattype floattype);


template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<3>(
    const std::vector<Node*> nodes, const CORE::LINALG::Matrix<3, 1>& xyz,
    CORE::FE::CellType distype, INPAR::CUT::CutFloattype floattype);
template Teuchos::RCP<CORE::GEO::CUT::Position> CORE::GEO::CUT::Position::Create<2>(
    const std::vector<Node*> nodes, const CORE::LINALG::Matrix<2, 1>& xyz,
    CORE::FE::CellType distype, INPAR::CUT::CutFloattype floattype);

/* --- ComputeEmbeddedPosition --- */
// embedded element types
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::tri3>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::tri6>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad4>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad8>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad9>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<2, CORE::FE::CellType::line2>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::line2>;

// non-embedded element types for the embedded case (only necessary due to compiler problems)
// template class
// CORE::GEO::CUT::ComputeEmbeddedPosition<3,CORE::FE::CellType::hex16>; template
// class CORE::GEO::CUT::ComputeEmbeddedPosition<3,CORE::FE::CellType::hex18>;
// template class
// CORE::GEO::CUT::ComputeEmbeddedPosition<3,CORE::FE::CellType::hex27>; template
// class CORE::GEO::CUT::ComputeEmbeddedPosition<3,CORE::FE::CellType::tet10>;
// template class
// CORE::GEO::CUT::ComputeEmbeddedPosition<3,CORE::FE::CellType::wedge15>;

template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::tri3,
    CORE::FE::num_nodes<CORE::FE::CellType::tri3>, CORE::FE::dim<CORE::FE::CellType::tri3>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::tri6,
    CORE::FE::num_nodes<CORE::FE::CellType::tri6>, CORE::FE::dim<CORE::FE::CellType::tri6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad4,
    CORE::FE::num_nodes<CORE::FE::CellType::quad4>, CORE::FE::dim<CORE::FE::CellType::quad4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad8,
    CORE::FE::num_nodes<CORE::FE::CellType::quad8>, CORE::FE::dim<CORE::FE::CellType::quad8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::quad9,
    CORE::FE::num_nodes<CORE::FE::CellType::quad9>, CORE::FE::dim<CORE::FE::CellType::quad9>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<2, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputeEmbeddedPosition<3, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;
/* --- ComputePosition --- */
// non-embedded cases (only)
template class CORE::GEO::CUT::ComputePosition<1, CORE::FE::CellType::line2>;

template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::tri3>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::tri6>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad4>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad8>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad9>;

template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::tet4>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::tet10>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex8>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex16>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex18>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex20>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex27>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::pyramid5>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::wedge6>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::wedge15>;

template class CORE::GEO::CUT::ComputePosition<1, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;

template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::tri3,
    CORE::FE::num_nodes<CORE::FE::CellType::tri3>, CORE::FE::dim<CORE::FE::CellType::tri3>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::tri6,
    CORE::FE::num_nodes<CORE::FE::CellType::tri6>, CORE::FE::dim<CORE::FE::CellType::tri6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad4,
    CORE::FE::num_nodes<CORE::FE::CellType::quad4>, CORE::FE::dim<CORE::FE::CellType::quad4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad8,
    CORE::FE::num_nodes<CORE::FE::CellType::quad8>, CORE::FE::dim<CORE::FE::CellType::quad8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<2, CORE::FE::CellType::quad9,
    CORE::FE::num_nodes<CORE::FE::CellType::quad9>, CORE::FE::dim<CORE::FE::CellType::quad9>,
    INPAR::CUT::floattype_cln>;

template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::tet4,
    CORE::FE::num_nodes<CORE::FE::CellType::tet4>, CORE::FE::dim<CORE::FE::CellType::tet4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::tet10,
    CORE::FE::num_nodes<CORE::FE::CellType::tet10>, CORE::FE::dim<CORE::FE::CellType::tet10>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex8,
    CORE::FE::num_nodes<CORE::FE::CellType::hex8>, CORE::FE::dim<CORE::FE::CellType::hex8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex16,
    CORE::FE::num_nodes<CORE::FE::CellType::hex16>, CORE::FE::dim<CORE::FE::CellType::hex16>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex18,
    CORE::FE::num_nodes<CORE::FE::CellType::hex18>, CORE::FE::dim<CORE::FE::CellType::hex18>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex20,
    CORE::FE::num_nodes<CORE::FE::CellType::hex20>, CORE::FE::dim<CORE::FE::CellType::hex20>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::hex27,
    CORE::FE::num_nodes<CORE::FE::CellType::hex27>, CORE::FE::dim<CORE::FE::CellType::hex27>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::pyramid5,
    CORE::FE::num_nodes<CORE::FE::CellType::pyramid5>, CORE::FE::dim<CORE::FE::CellType::pyramid5>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::wedge6,
    CORE::FE::num_nodes<CORE::FE::CellType::wedge6>, CORE::FE::dim<CORE::FE::CellType::wedge6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::ComputePosition<3, CORE::FE::CellType::wedge15,
    CORE::FE::num_nodes<CORE::FE::CellType::wedge15>, CORE::FE::dim<CORE::FE::CellType::wedge15>,
    INPAR::CUT::floattype_cln>;

/* --- PositionGeneric --- */
// embedded cases
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tri3>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tri6>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad4>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad8>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad9>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::line2>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::line2>;

// non-embedded cases
template class CORE::GEO::CUT::PositionGeneric<1, CORE::FE::CellType::line2>;

template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::tri3>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::tri6>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad4>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad8>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad9>;

template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tet4>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tet10>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex8>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex16>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex18>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex20>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex27>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::pyramid5>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::wedge6>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::wedge15>;

// unused / impossible cases (only necessary due to compiler problems)
// template class CORE::GEO::CUT::PositionGeneric<2,CORE::FE::CellType::hex8>;
// template class CORE::GEO::CUT::PositionGeneric<2,CORE::FE::CellType::tet4>;
// template class CORE::GEO::CUT::PositionGeneric<2,CORE::FE::CellType::pyramid5>;
// template class CORE::GEO::CUT::PositionGeneric<2,CORE::FE::CellType::wedge6>;

// embedded cases
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tri3,
    CORE::FE::num_nodes<CORE::FE::CellType::tri3>, CORE::FE::dim<CORE::FE::CellType::tri3>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tri6,
    CORE::FE::num_nodes<CORE::FE::CellType::tri6>, CORE::FE::dim<CORE::FE::CellType::tri6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad4,
    CORE::FE::num_nodes<CORE::FE::CellType::quad4>, CORE::FE::dim<CORE::FE::CellType::quad4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad8,
    CORE::FE::num_nodes<CORE::FE::CellType::quad8>, CORE::FE::dim<CORE::FE::CellType::quad8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::quad9,
    CORE::FE::num_nodes<CORE::FE::CellType::quad9>, CORE::FE::dim<CORE::FE::CellType::quad9>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;

// non-embedded cases
template class CORE::GEO::CUT::PositionGeneric<1, CORE::FE::CellType::line2,
    CORE::FE::num_nodes<CORE::FE::CellType::line2>, CORE::FE::dim<CORE::FE::CellType::line2>,
    INPAR::CUT::floattype_cln>;

template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::tri3,
    CORE::FE::num_nodes<CORE::FE::CellType::tri3>, CORE::FE::dim<CORE::FE::CellType::tri3>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::tri6,
    CORE::FE::num_nodes<CORE::FE::CellType::tri6>, CORE::FE::dim<CORE::FE::CellType::tri6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad4,
    CORE::FE::num_nodes<CORE::FE::CellType::quad4>, CORE::FE::dim<CORE::FE::CellType::quad4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad8,
    CORE::FE::num_nodes<CORE::FE::CellType::quad8>, CORE::FE::dim<CORE::FE::CellType::quad8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<2, CORE::FE::CellType::quad9,
    CORE::FE::num_nodes<CORE::FE::CellType::quad9>, CORE::FE::dim<CORE::FE::CellType::quad9>,
    INPAR::CUT::floattype_cln>;

template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tet4,
    CORE::FE::num_nodes<CORE::FE::CellType::tet4>, CORE::FE::dim<CORE::FE::CellType::tet4>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::tet10,
    CORE::FE::num_nodes<CORE::FE::CellType::tet10>, CORE::FE::dim<CORE::FE::CellType::tet10>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex8,
    CORE::FE::num_nodes<CORE::FE::CellType::hex8>, CORE::FE::dim<CORE::FE::CellType::hex8>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex16,
    CORE::FE::num_nodes<CORE::FE::CellType::hex16>, CORE::FE::dim<CORE::FE::CellType::hex16>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex18,
    CORE::FE::num_nodes<CORE::FE::CellType::hex18>, CORE::FE::dim<CORE::FE::CellType::hex18>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex20,
    CORE::FE::num_nodes<CORE::FE::CellType::hex20>, CORE::FE::dim<CORE::FE::CellType::hex20>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::hex27,
    CORE::FE::num_nodes<CORE::FE::CellType::hex27>, CORE::FE::dim<CORE::FE::CellType::hex27>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::pyramid5,
    CORE::FE::num_nodes<CORE::FE::CellType::pyramid5>, CORE::FE::dim<CORE::FE::CellType::pyramid5>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::wedge6,
    CORE::FE::num_nodes<CORE::FE::CellType::wedge6>, CORE::FE::dim<CORE::FE::CellType::wedge6>,
    INPAR::CUT::floattype_cln>;
template class CORE::GEO::CUT::PositionGeneric<3, CORE::FE::CellType::wedge15,
    CORE::FE::num_nodes<CORE::FE::CellType::wedge15>, CORE::FE::dim<CORE::FE::CellType::wedge15>,
    INPAR::CUT::floattype_cln>;

FOUR_C_NAMESPACE_CLOSE
