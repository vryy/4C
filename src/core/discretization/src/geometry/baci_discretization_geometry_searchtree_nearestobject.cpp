/*----------------------------------------------------------------------*/
/*! \file

\brief stores data about nearest object in oct tree

\level 3

*----------------------------------------------------------------------*/

#include "baci_discretization_geometry_searchtree_nearestobject.hpp"

BACI_NAMESPACE_OPEN


CORE::GEO::NearestObject::NearestObject()
    : objectType_(NOTYPE_OBJECT), nodeId_(-1), lineId_(-1), surfId_(-1), label_(-1)
{
  physcoord_.PutScalar(0.0);
  return;
}


CORE::GEO::NearestObject::NearestObject(const CORE::GEO::NearestObject& old)
    : objectType_(old.objectType_),
      nodeId_(old.nodeId_),
      lineId_(old.lineId_),
      surfId_(old.surfId_),
      label_(old.label_),
      physcoord_(old.physcoord_)
{
  return;
}


CORE::GEO::NearestObject& CORE::GEO::NearestObject::operator=(const CORE::GEO::NearestObject& old)
{
  objectType_ = old.objectType_;
  nodeId_ = old.nodeId_;
  lineId_ = old.lineId_;
  surfId_ = old.surfId_;
  label_ = old.label_;
  physcoord_ = old.physcoord_;
  return *this;
}


void CORE::GEO::NearestObject::clear()
{
  objectType_ = NOTYPE_OBJECT;
  nodeId_ = -1;
  lineId_ = -1;
  surfId_ = -1;
  label_ = -1;
  physcoord_.PutScalar(0.0);
  return;
}


void CORE::GEO::NearestObject::setNodeObjectType(
    const int nodeId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = NODE_OBJECT;
  nodeId_ = nodeId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  lineId_ = -1;
  surfId_ = -1;
}


void CORE::GEO::NearestObject::setLineObjectType(const int lineId, const int surfId,
    const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = LINE_OBJECT;
  lineId_ = lineId;
  surfId_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  nodeId_ = -1;
}


void CORE::GEO::NearestObject::setSurfaceObjectType(
    const int surfId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = SURFACE_OBJECT;
  surfId_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  nodeId_ = -1;
  lineId_ = -1;
}

BACI_NAMESPACE_CLOSE
