/*----------------------------------------------------------------------------*/
/*! \file

\brief Class to check whether a point lies inside an element, and also to
       compute local coordinates of a point with respect to the element.
       The standard as well as the embedded cases are treated here. Feel
       free to extend the content for your needs.

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_cut_position.hpp"

#include "4C_cut_boundingbox.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create(
    const Element& element, const Point& point, Core::Geo::Cut::CutFloatType floattype)
{
  const PositionFactory factory;
  return factory.create_position(element, point, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create(const Element& element,
    const Core::LinAlg::Matrix<rdim, 1>& xyz, Core::Geo::Cut::CutFloatType floattype)
{
  const PositionFactory factory;
  if (rdim < factory.n_prob_dim())
    FOUR_C_THROW(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.n_prob_dim());
  return factory.create_position(element, xyz.data(), floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim, unsigned cdim, unsigned rdim_2>
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create(
    const Core::LinAlg::Matrix<rdim, cdim>& xyze, const Core::LinAlg::Matrix<rdim_2, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.n_prob_dim();
  const unsigned num_nodes_ele = Core::FE::getNumberOfElementNodes(distype);

  if (rdim < probdim or cdim != num_nodes_ele)
    FOUR_C_THROW(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.m(), xyze.n());

  const double* xyze_ptr = xyze.data();
  Core::LinAlg::SerialDenseMatrix xyze_eptra;
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

  return factory.create_position(xyze_ptr, xyz.data(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create(
    const Core::LinAlg::SerialDenseMatrix& xyze, const Core::LinAlg::Matrix<rdim, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype)
{
  const PositionFactory factory;
  const unsigned probdim = factory.n_prob_dim();
  const unsigned num_nodes_ele = Core::FE::getNumberOfElementNodes(distype);

  if (static_cast<unsigned>(xyze.numRows()) < probdim or
      static_cast<unsigned>(xyze.numCols()) != num_nodes_ele)
    FOUR_C_THROW(
        "Dimension mismatch of xyze! \n"
        "expected input: %d x %d (rows x cols)\n"
        "received input : %d x %d (rows x cols)",
        probdim, num_nodes_ele, xyze.numRows(), xyze.numCols());

  const double* xyze_ptr = xyze.values();
  Core::LinAlg::SerialDenseMatrix xyze_eptra;
  if (static_cast<unsigned>(xyze.numRows()) > probdim)
  {
    xyze_eptra.shape(probdim, num_nodes_ele);
    FixMatrixShape(xyze, xyze_eptra);
    xyze_ptr = xyze_eptra.values();
  }

  if (xyz.num_rows() < probdim)
    FOUR_C_THROW(
        "Dimension mismatch of xyz! \n"
        "expected input: %d x 1 (rows x cols)\n"
        "received input : %d x 1 (rows x cols)",
        probdim, xyz.num_rows());

  return factory.create_position(xyze_ptr, xyz.data(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned rdim>
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create(
    const std::vector<Node*> nodes, const Core::LinAlg::Matrix<rdim, 1>& xyz,
    Core::FE::CellType distype, Core::Geo::Cut::CutFloatType floattype)
{
  const PositionFactory factory;
  if (rdim < factory.n_prob_dim())
    FOUR_C_THROW(
        "The given point has the wrong row dimension!\n"
        "rdim < prodbim <--> %d < %d",
        rdim, factory.n_prob_dim());
  return factory.create_position(nodes, xyz.data(), distype, floattype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, Core::FE::CellType eletype, unsigned num_nodes_element, unsigned dim,
    Core::Geo::Cut::CutFloatType floattype>
void Core::Geo::Cut::PositionGeneric<probdim, eletype, num_nodes_element, dim,
    floattype>::construct_bounding_box()
{
  bbside_ = Teuchos::rcp(BoundingBox::create());

  for (unsigned i = 0; i < num_nodes_element; ++i)
  {
    Core::LinAlg::Matrix<3, 1> x1(&this->xyze_(0, i), true);
    bbside_->add_point(x1);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, Core::FE::CellType eletype, unsigned num_nodes_element, unsigned dim,
    Core::Geo::Cut::CutFloatType floattype>
bool Core::Geo::Cut::ComputePosition<probdim, eletype, num_nodes_element, dim, floattype>::compute(
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
  Kernel::ComputePosition<probdim, eletype, floattype> cp(this->xsi_);
  //  Kernel::DebugComputePosition<probdim,eletype,floattype> cp( this->xsi_ );
  this->pos_status_ =
      (cp(this->xyze_, this->px_) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cp.get_tolerance();

  if (this->pos_status_ == Position::position_valid) return within_limits_tol(Tol);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, Core::FE::CellType eletype, unsigned num_nodes_element, unsigned dim,
    Core::Geo::Cut::CutFloatType floattype>
bool Core::Geo::Cut::ComputeEmbeddedPosition<probdim, eletype, num_nodes_element, dim,
    floattype>::is_given_point_within_element()
{
  // If the given point is outside the side's bounding box, no need to perform
  // the complex calculations
  if (not this->bbside_->within(1.0, &this->px_(0, 0)))
  {
    this->pos_status_ = Position::position_outside_of_bbox;
    return false;
  }

  // try to compute the local coordinates and the distance using the Newton scheme
  Kernel::ComputeDistance<probdim, eletype, (floattype == floattype_cln)> cd(xsi_aug_);
  //  Kernel::DebugComputeDistance< probdim, eletype,(floattype==floattype_cln) > cd(
  //  xsi_aug_ );

  double dist = 0.0;
  this->pos_status_ =
      (cd(this->xyze_, this->px_, dist) ? Position::position_valid : Position::position_invalid);
  this->compute_tolerance_ = cd.get_tolerance();

  return (this->pos_status_ == Position::position_valid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, Core::FE::CellType eletype, unsigned num_nodes_element, unsigned dim,
    Core::Geo::Cut::CutFloatType floattype>
bool Core::Geo::Cut::ComputeEmbeddedPosition<probdim, eletype, num_nodes_element, dim,
    floattype>::compute(const double& Tol, const bool& allow_dist)
{
  xsi_aug_ = 0.0;
  /* If the given point is outside the side's bounding box, no need
   * to perform the complex calculations */
  if (not allow_dist)
  {
    if (not this->bbside_->within(1.0, &this->px_(0, 0)))
    {
      this->pos_status_ = Position::position_outside_of_bbox;
      return false;
    }
  }

  // try to compute the local coordinates and the distance using the Newton scheme
  Kernel::ComputeDistance<probdim, eletype, (floattype == floattype_cln)> cd(xsi_aug_);
  //  Kernel::DebugComputeDistance< probdim, eletype,(floattype==floattype_cln) > cd(
  //  xsi_aug_ );

  double dist = 0.0;
  this->pos_status_ =
      (cd(this->xyze_, this->px_, dist) ? Position::position_valid : Position::position_invalid);

  this->compute_tolerance_ = cd.get_tolerance();

  // if newton did not converge, try some other checks
  if (this->pos_status_ != Position::position_valid)
  {
    // if the surface has no area zero area
    if (cd.zero_area())
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
      double res_nrm2 = cd.get_residual_l2_norm();
      if (res_nrm2 < 1.0)
      {
        std::cout
            << "==| WARNING: compute_distance didn't converge, but residual is smaller than 1 "
               "(calculated distance is used) |=="
            << std::endl;
        this->pos_status_ = Position::position_distance_valid;
      }

      // as the local coordinates didn't converge
      return false;
    }
  }

  // Std Return in case Compute Distance converged!
  return within_limits_tol(Tol, allow_dist);
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::PositionFactory::PositionFactory() : probdim_(Global::Problem::instance()->n_dim())
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::PositionFactory::create_position(
    const Element& element, const Point& point, Core::Geo::Cut::CutFloatType floattype) const
{
  Core::FE::CellType distype = element.shape();

  switch (distype)
  {
    case Core::FE::CellType::line2:
      return create_concrete_position<Core::FE::CellType::line2>(element, point, floattype);
    case Core::FE::CellType::tri3:
      return create_concrete_position<Core::FE::CellType::tri3>(element, point, floattype);
    case Core::FE::CellType::tri6:
      return create_concrete_position<Core::FE::CellType::tri6>(element, point, floattype);
    case Core::FE::CellType::quad4:
      return create_concrete_position<Core::FE::CellType::quad4>(element, point, floattype);
    case Core::FE::CellType::quad8:
      return create_concrete_position<Core::FE::CellType::quad8>(element, point, floattype);
    case Core::FE::CellType::quad9:
      return create_concrete_position<Core::FE::CellType::quad9>(element, point, floattype);
    case Core::FE::CellType::hex8:
      return create_concrete_position<Core::FE::CellType::hex8>(element, point, floattype);
    case Core::FE::CellType::hex20:
      return create_concrete_position<Core::FE::CellType::hex20>(element, point, floattype);
    case Core::FE::CellType::tet4:
      return create_concrete_position<Core::FE::CellType::tet4>(element, point, floattype);
    case Core::FE::CellType::pyramid5:
      return create_concrete_position<Core::FE::CellType::pyramid5>(element, point, floattype);
    case Core::FE::CellType::wedge6:
      return create_concrete_position<Core::FE::CellType::wedge6>(element, point, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::PositionFactory::create_position(
    const Element& element, const double* xyz, Core::Geo::Cut::CutFloatType floattype) const
{
  Core::FE::CellType distype = element.shape();

  switch (distype)
  {
    case Core::FE::CellType::line2:
      return create_concrete_position<Core::FE::CellType::line2>(element, xyz, floattype);
    case Core::FE::CellType::tri3:
      return create_concrete_position<Core::FE::CellType::tri3>(element, xyz, floattype);
    case Core::FE::CellType::tri6:
      return create_concrete_position<Core::FE::CellType::tri6>(element, xyz, floattype);
    case Core::FE::CellType::quad4:
      return create_concrete_position<Core::FE::CellType::quad4>(element, xyz, floattype);
    case Core::FE::CellType::quad8:
      return create_concrete_position<Core::FE::CellType::quad8>(element, xyz, floattype);
    case Core::FE::CellType::quad9:
      return create_concrete_position<Core::FE::CellType::quad9>(element, xyz, floattype);
    case Core::FE::CellType::hex8:
      return create_concrete_position<Core::FE::CellType::hex8>(element, xyz, floattype);
    case Core::FE::CellType::hex20:
      return create_concrete_position<Core::FE::CellType::hex20>(element, xyz, floattype);
    case Core::FE::CellType::tet4:
      return create_concrete_position<Core::FE::CellType::tet4>(element, xyz, floattype);
    case Core::FE::CellType::pyramid5:
      return create_concrete_position<Core::FE::CellType::pyramid5>(element, xyz, floattype);
    case Core::FE::CellType::wedge6:
      return create_concrete_position<Core::FE::CellType::wedge6>(element, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::PositionFactory::create_position(
    const double* xyze, const double* xyz, const Core::FE::CellType& distype,
    Core::Geo::Cut::CutFloatType floattype) const
{
  switch (distype)
  {
    case Core::FE::CellType::line2:
      return create_concrete_position<Core::FE::CellType::line2>(xyze, xyz, floattype);
    case Core::FE::CellType::tri3:
      return create_concrete_position<Core::FE::CellType::tri3>(xyze, xyz, floattype);
    case Core::FE::CellType::tri6:
      return create_concrete_position<Core::FE::CellType::tri6>(xyze, xyz, floattype);
    case Core::FE::CellType::quad4:
      return create_concrete_position<Core::FE::CellType::quad4>(xyze, xyz, floattype);
    case Core::FE::CellType::quad8:
      return create_concrete_position<Core::FE::CellType::quad8>(xyze, xyz, floattype);
    case Core::FE::CellType::quad9:
      return create_concrete_position<Core::FE::CellType::quad9>(xyze, xyz, floattype);
    case Core::FE::CellType::hex8:
      return create_concrete_position<Core::FE::CellType::hex8>(xyze, xyz, floattype);
    case Core::FE::CellType::hex20:
      return create_concrete_position<Core::FE::CellType::hex20>(xyze, xyz, floattype);
    case Core::FE::CellType::tet4:
      return create_concrete_position<Core::FE::CellType::tet4>(xyze, xyz, floattype);
    case Core::FE::CellType::pyramid5:
      return create_concrete_position<Core::FE::CellType::pyramid5>(xyze, xyz, floattype);
    case Core::FE::CellType::wedge6:
      return create_concrete_position<Core::FE::CellType::wedge6>(xyze, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::PositionFactory::create_position(
    const std::vector<Core::Geo::Cut::Node*> nodes, const double* xyz, Core::FE::CellType distype,
    Core::Geo::Cut::CutFloatType floattype) const
{
  if (distype == Core::FE::CellType::dis_none)
  {
    plain_element_set elements;
    FindCommonElements(nodes, elements);
    if (elements.size() != 1)
      FOUR_C_THROW("Couldn't find a unique element corresponding to the given nodes.");

    distype = elements[0]->shape();
  }

  switch (distype)
  {
    case Core::FE::CellType::line2:
      return create_concrete_position<Core::FE::CellType::line2>(nodes, xyz, floattype);
    case Core::FE::CellType::tri3:
      return create_concrete_position<Core::FE::CellType::tri3>(nodes, xyz, floattype);
    case Core::FE::CellType::tri6:
      return create_concrete_position<Core::FE::CellType::tri6>(nodes, xyz, floattype);
    case Core::FE::CellType::quad4:
      return create_concrete_position<Core::FE::CellType::quad4>(nodes, xyz, floattype);
    case Core::FE::CellType::quad8:
      return create_concrete_position<Core::FE::CellType::quad8>(nodes, xyz, floattype);
    case Core::FE::CellType::quad9:
      return create_concrete_position<Core::FE::CellType::quad9>(nodes, xyz, floattype);
    case Core::FE::CellType::hex8:
      return create_concrete_position<Core::FE::CellType::hex8>(nodes, xyz, floattype);
    case Core::FE::CellType::hex20:
      return create_concrete_position<Core::FE::CellType::hex20>(nodes, xyz, floattype);
    case Core::FE::CellType::tet4:
      return create_concrete_position<Core::FE::CellType::tet4>(nodes, xyz, floattype);
    case Core::FE::CellType::pyramid5:
      return create_concrete_position<Core::FE::CellType::pyramid5>(nodes, xyz, floattype);
    case Core::FE::CellType::wedge6:
      return create_concrete_position<Core::FE::CellType::wedge6>(nodes, xyz, floattype);
    default:
      FOUR_C_THROW("Unsupported distype = %s", Core::FE::CellTypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::CutFloatType Core::Geo::Cut::PositionFactory::use_pos_floattype(
    Core::Geo::Cut::CutFloatType floattype)
{
  if (general_pos_floattype_ != floattype_none)
    return general_pos_floattype_;
  else
    return floattype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::CutFloatType Core::Geo::Cut::PositionFactory::use_dist_floattype(
    Core::Geo::Cut::CutFloatType floattype)
{
  if (general_dist_floattype_ != floattype_none)
    return general_dist_floattype_;
  else
    return floattype;
}

Core::Geo::Cut::CutFloatType Core::Geo::Cut::PositionFactory::general_pos_floattype_ =
    floattype_none;
Core::Geo::Cut::CutFloatType Core::Geo::Cut::PositionFactory::general_dist_floattype_ =
    floattype_none;

template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<2>(
    const Element& element, const Core::LinAlg::Matrix<2, 1>& xyz,
    Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3>(
    const Element& element, const Core::LinAlg::Matrix<3, 1>& xyz,
    Core::Geo::Cut::CutFloatType floattype);

template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 3, 3>(
    const Core::LinAlg::Matrix<3, 3>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 6, 3>(
    const Core::LinAlg::Matrix<3, 6>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 4, 3>(
    const Core::LinAlg::Matrix<3, 4>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 8, 3>(
    const Core::LinAlg::Matrix<3, 8>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 9, 3>(
    const Core::LinAlg::Matrix<3, 9>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3, 2, 3>(
    const Core::LinAlg::Matrix<3, 2>& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<2, 2, 2>(
    const Core::LinAlg::Matrix<2, 2>& xyze, const Core::LinAlg::Matrix<2, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);

template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3>(
    const Core::LinAlg::SerialDenseMatrix& xyze, const Core::LinAlg::Matrix<3, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<2>(
    const Core::LinAlg::SerialDenseMatrix& xyze, const Core::LinAlg::Matrix<2, 1>& xyz,
    const Core::FE::CellType& distype, Core::Geo::Cut::CutFloatType floattype);


template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<3>(
    const std::vector<Node*> nodes, const Core::LinAlg::Matrix<3, 1>& xyz,
    Core::FE::CellType distype, Core::Geo::Cut::CutFloatType floattype);
template Teuchos::RCP<Core::Geo::Cut::Position> Core::Geo::Cut::Position::create<2>(
    const std::vector<Node*> nodes, const Core::LinAlg::Matrix<2, 1>& xyz,
    Core::FE::CellType distype, Core::Geo::Cut::CutFloatType floattype);

/* --- ComputeEmbeddedPosition --- */
// embedded element types
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::tri3>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::tri6>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad4>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad8>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad9>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<2, Core::FE::CellType::line2>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::line2>;

// non-embedded element types for the embedded case (only necessary due to compiler problems)
// template class
// Core::Geo::Cut::ComputeEmbeddedPosition<3,Core::FE::CellType::hex16>; template
// class Core::Geo::Cut::ComputeEmbeddedPosition<3,Core::FE::CellType::hex18>;
// template class
// Core::Geo::Cut::ComputeEmbeddedPosition<3,Core::FE::CellType::hex27>; template
// class Core::Geo::Cut::ComputeEmbeddedPosition<3,Core::FE::CellType::tet10>;
// template class
// Core::Geo::Cut::ComputeEmbeddedPosition<3,Core::FE::CellType::wedge15>;

template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::tri3,
    Core::FE::num_nodes<Core::FE::CellType::tri3>, Core::FE::dim<Core::FE::CellType::tri3>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::tri6,
    Core::FE::num_nodes<Core::FE::CellType::tri6>, Core::FE::dim<Core::FE::CellType::tri6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad4,
    Core::FE::num_nodes<Core::FE::CellType::quad4>, Core::FE::dim<Core::FE::CellType::quad4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad8,
    Core::FE::num_nodes<Core::FE::CellType::quad8>, Core::FE::dim<Core::FE::CellType::quad8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::quad9,
    Core::FE::num_nodes<Core::FE::CellType::quad9>, Core::FE::dim<Core::FE::CellType::quad9>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<2, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputeEmbeddedPosition<3, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;
/* --- ComputePosition --- */
// non-embedded cases (only)
template class Core::Geo::Cut::ComputePosition<1, Core::FE::CellType::line2>;

template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::tri3>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::tri6>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad4>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad8>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad9>;

template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::tet4>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::tet10>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex8>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex16>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex18>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex20>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex27>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::pyramid5>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::wedge6>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::wedge15>;

template class Core::Geo::Cut::ComputePosition<1, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;

template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::tri3,
    Core::FE::num_nodes<Core::FE::CellType::tri3>, Core::FE::dim<Core::FE::CellType::tri3>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::tri6,
    Core::FE::num_nodes<Core::FE::CellType::tri6>, Core::FE::dim<Core::FE::CellType::tri6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad4,
    Core::FE::num_nodes<Core::FE::CellType::quad4>, Core::FE::dim<Core::FE::CellType::quad4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad8,
    Core::FE::num_nodes<Core::FE::CellType::quad8>, Core::FE::dim<Core::FE::CellType::quad8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<2, Core::FE::CellType::quad9,
    Core::FE::num_nodes<Core::FE::CellType::quad9>, Core::FE::dim<Core::FE::CellType::quad9>,
    Core::Geo::Cut::floattype_cln>;

template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::tet4,
    Core::FE::num_nodes<Core::FE::CellType::tet4>, Core::FE::dim<Core::FE::CellType::tet4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::tet10,
    Core::FE::num_nodes<Core::FE::CellType::tet10>, Core::FE::dim<Core::FE::CellType::tet10>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex8,
    Core::FE::num_nodes<Core::FE::CellType::hex8>, Core::FE::dim<Core::FE::CellType::hex8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex16,
    Core::FE::num_nodes<Core::FE::CellType::hex16>, Core::FE::dim<Core::FE::CellType::hex16>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex18,
    Core::FE::num_nodes<Core::FE::CellType::hex18>, Core::FE::dim<Core::FE::CellType::hex18>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex20,
    Core::FE::num_nodes<Core::FE::CellType::hex20>, Core::FE::dim<Core::FE::CellType::hex20>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::hex27,
    Core::FE::num_nodes<Core::FE::CellType::hex27>, Core::FE::dim<Core::FE::CellType::hex27>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::pyramid5,
    Core::FE::num_nodes<Core::FE::CellType::pyramid5>, Core::FE::dim<Core::FE::CellType::pyramid5>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::wedge6,
    Core::FE::num_nodes<Core::FE::CellType::wedge6>, Core::FE::dim<Core::FE::CellType::wedge6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::ComputePosition<3, Core::FE::CellType::wedge15,
    Core::FE::num_nodes<Core::FE::CellType::wedge15>, Core::FE::dim<Core::FE::CellType::wedge15>,
    Core::Geo::Cut::floattype_cln>;

/* --- PositionGeneric --- */
// embedded cases
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tri3>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tri6>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad4>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad8>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad9>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::line2>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::line2>;

// non-embedded cases
template class Core::Geo::Cut::PositionGeneric<1, Core::FE::CellType::line2>;

template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::tri3>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::tri6>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad4>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad8>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad9>;

template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tet4>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tet10>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex8>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex16>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex18>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex20>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex27>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::pyramid5>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::wedge6>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::wedge15>;

// unused / impossible cases (only necessary due to compiler problems)
// template class Core::Geo::Cut::PositionGeneric<2,Core::FE::CellType::hex8>;
// template class Core::Geo::Cut::PositionGeneric<2,Core::FE::CellType::tet4>;
// template class Core::Geo::Cut::PositionGeneric<2,Core::FE::CellType::pyramid5>;
// template class Core::Geo::Cut::PositionGeneric<2,Core::FE::CellType::wedge6>;

// embedded cases
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tri3,
    Core::FE::num_nodes<Core::FE::CellType::tri3>, Core::FE::dim<Core::FE::CellType::tri3>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tri6,
    Core::FE::num_nodes<Core::FE::CellType::tri6>, Core::FE::dim<Core::FE::CellType::tri6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad4,
    Core::FE::num_nodes<Core::FE::CellType::quad4>, Core::FE::dim<Core::FE::CellType::quad4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad8,
    Core::FE::num_nodes<Core::FE::CellType::quad8>, Core::FE::dim<Core::FE::CellType::quad8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::quad9,
    Core::FE::num_nodes<Core::FE::CellType::quad9>, Core::FE::dim<Core::FE::CellType::quad9>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;

// non-embedded cases
template class Core::Geo::Cut::PositionGeneric<1, Core::FE::CellType::line2,
    Core::FE::num_nodes<Core::FE::CellType::line2>, Core::FE::dim<Core::FE::CellType::line2>,
    Core::Geo::Cut::floattype_cln>;

template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::tri3,
    Core::FE::num_nodes<Core::FE::CellType::tri3>, Core::FE::dim<Core::FE::CellType::tri3>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::tri6,
    Core::FE::num_nodes<Core::FE::CellType::tri6>, Core::FE::dim<Core::FE::CellType::tri6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad4,
    Core::FE::num_nodes<Core::FE::CellType::quad4>, Core::FE::dim<Core::FE::CellType::quad4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad8,
    Core::FE::num_nodes<Core::FE::CellType::quad8>, Core::FE::dim<Core::FE::CellType::quad8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<2, Core::FE::CellType::quad9,
    Core::FE::num_nodes<Core::FE::CellType::quad9>, Core::FE::dim<Core::FE::CellType::quad9>,
    Core::Geo::Cut::floattype_cln>;

template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tet4,
    Core::FE::num_nodes<Core::FE::CellType::tet4>, Core::FE::dim<Core::FE::CellType::tet4>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::tet10,
    Core::FE::num_nodes<Core::FE::CellType::tet10>, Core::FE::dim<Core::FE::CellType::tet10>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex8,
    Core::FE::num_nodes<Core::FE::CellType::hex8>, Core::FE::dim<Core::FE::CellType::hex8>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex16,
    Core::FE::num_nodes<Core::FE::CellType::hex16>, Core::FE::dim<Core::FE::CellType::hex16>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex18,
    Core::FE::num_nodes<Core::FE::CellType::hex18>, Core::FE::dim<Core::FE::CellType::hex18>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex20,
    Core::FE::num_nodes<Core::FE::CellType::hex20>, Core::FE::dim<Core::FE::CellType::hex20>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::hex27,
    Core::FE::num_nodes<Core::FE::CellType::hex27>, Core::FE::dim<Core::FE::CellType::hex27>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::pyramid5,
    Core::FE::num_nodes<Core::FE::CellType::pyramid5>, Core::FE::dim<Core::FE::CellType::pyramid5>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::wedge6,
    Core::FE::num_nodes<Core::FE::CellType::wedge6>, Core::FE::dim<Core::FE::CellType::wedge6>,
    Core::Geo::Cut::floattype_cln>;
template class Core::Geo::Cut::PositionGeneric<3, Core::FE::CellType::wedge15,
    Core::FE::num_nodes<Core::FE::CellType::wedge15>, Core::FE::dim<Core::FE::CellType::wedge15>,
    Core::Geo::Cut::floattype_cln>;

FOUR_C_NAMESPACE_CLOSE
