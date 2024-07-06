/*----------------------------------------------------------------------*/
/*! \file

\brief stores data about nearest object in oct tree

\level 3

*----------------------------------------------------------------------*/

#include "4C_fem_geometry_searchtree_nearestobject.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Geo::NearestObject::NearestObject()
    : object_type_(NOTYPE_OBJECT), node_id_(-1), line_id_(-1), surf_id_(-1), label_(-1)
{
  physcoord_.put_scalar(0.0);
  return;
}


Core::Geo::NearestObject::NearestObject(const Core::Geo::NearestObject& old)
    : object_type_(old.object_type_),
      node_id_(old.node_id_),
      line_id_(old.line_id_),
      surf_id_(old.surf_id_),
      label_(old.label_),
      physcoord_(old.physcoord_)
{
  return;
}


Core::Geo::NearestObject& Core::Geo::NearestObject::operator=(const Core::Geo::NearestObject& old)
{
  object_type_ = old.object_type_;
  node_id_ = old.node_id_;
  line_id_ = old.line_id_;
  surf_id_ = old.surf_id_;
  label_ = old.label_;
  physcoord_ = old.physcoord_;
  return *this;
}


void Core::Geo::NearestObject::clear()
{
  object_type_ = NOTYPE_OBJECT;
  node_id_ = -1;
  line_id_ = -1;
  surf_id_ = -1;
  label_ = -1;
  physcoord_.put_scalar(0.0);
  return;
}


void Core::Geo::NearestObject::set_node_object_type(
    const int nodeId, const int label, const Core::LinAlg::Matrix<3, 1>& physcoord)
{
  object_type_ = NODE_OBJECT;
  node_id_ = nodeId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  line_id_ = -1;
  surf_id_ = -1;
}


void Core::Geo::NearestObject::set_line_object_type(const int lineId, const int surfId,
    const int label, const Core::LinAlg::Matrix<3, 1>& physcoord)
{
  object_type_ = LINE_OBJECT;
  line_id_ = lineId;
  surf_id_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  node_id_ = -1;
}


void Core::Geo::NearestObject::set_surface_object_type(
    const int surfId, const int label, const Core::LinAlg::Matrix<3, 1>& physcoord)
{
  object_type_ = SURFACE_OBJECT;
  surf_id_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  node_id_ = -1;
  line_id_ = -1;
}

FOUR_C_NAMESPACE_CLOSE
