/*----------------------------------------------------------------------*/
/*! \file

\brief stores data about nearest object in oct tree

\level 3

*----------------------------------------------------------------------*/

#include "baci_discretization_geometry_searchtree_nearestobject.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::GEO::NearestObject::NearestObject()
    : object_type_(NOTYPE_OBJECT), node_id_(-1), line_id_(-1), surf_id_(-1), label_(-1)
{
  physcoord_.PutScalar(0.0);
  return;
}


CORE::GEO::NearestObject::NearestObject(const CORE::GEO::NearestObject& old)
    : object_type_(old.object_type_),
      node_id_(old.node_id_),
      line_id_(old.line_id_),
      surf_id_(old.surf_id_),
      label_(old.label_),
      physcoord_(old.physcoord_)
{
  return;
}


CORE::GEO::NearestObject& CORE::GEO::NearestObject::operator=(const CORE::GEO::NearestObject& old)
{
  object_type_ = old.object_type_;
  node_id_ = old.node_id_;
  line_id_ = old.line_id_;
  surf_id_ = old.surf_id_;
  label_ = old.label_;
  physcoord_ = old.physcoord_;
  return *this;
}


void CORE::GEO::NearestObject::clear()
{
  object_type_ = NOTYPE_OBJECT;
  node_id_ = -1;
  line_id_ = -1;
  surf_id_ = -1;
  label_ = -1;
  physcoord_.PutScalar(0.0);
  return;
}


void CORE::GEO::NearestObject::setNodeObjectType(
    const int nodeId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
{
  object_type_ = NODE_OBJECT;
  node_id_ = nodeId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  line_id_ = -1;
  surf_id_ = -1;
}


void CORE::GEO::NearestObject::setLineObjectType(const int lineId, const int surfId,
    const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
{
  object_type_ = LINE_OBJECT;
  line_id_ = lineId;
  surf_id_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  node_id_ = -1;
}


void CORE::GEO::NearestObject::setSurfaceObjectType(
    const int surfId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
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
