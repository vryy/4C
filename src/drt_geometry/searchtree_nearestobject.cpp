/*----------------------------------------------------------------------*/
/*! \file

\brief stores data about nearest object in oct tree

\level 3

\maintainer Martin Kronbichler
 */


#include "../drt_geometry/searchtree_nearestobject.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 08/08|
 *----------------------------------------------------------------------*/
GEO::NearestObject::NearestObject()
    : objectType_(NOTYPE_OBJECT), nodeId_(-1), lineId_(-1), surfId_(-1), label_(-1)
{
  physcoord_.PutScalar(0.0);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       u.may 08/08|
 *----------------------------------------------------------------------*/
GEO::NearestObject::NearestObject(const GEO::NearestObject& old)
    : objectType_(old.objectType_),
      nodeId_(old.nodeId_),
      lineId_(old.lineId_),
      surfId_(old.surfId_),
      label_(old.label_),
      physcoord_(old.physcoord_)
{
  return;
}



/*----------------------------------------------------------------------*
 * assignment operator                                       u.may 08/08|
 *----------------------------------------------------------------------*/
GEO::NearestObject& GEO::NearestObject::operator=(const GEO::NearestObject& old)
{
  objectType_ = old.objectType_;
  nodeId_ = old.nodeId_;
  lineId_ = old.lineId_;
  surfId_ = old.surfId_;
  label_ = old.label_;
  physcoord_ = old.physcoord_;
  return *this;
}


/*----------------------------------------------------------------------*
 * clear nearest object                                      u.may 05/09|
 *----------------------------------------------------------------------*/
void GEO::NearestObject::clear()
{
  objectType_ = NOTYPE_OBJECT;
  nodeId_ = -1;
  lineId_ = -1;
  surfId_ = -1;
  label_ = -1;
  physcoord_.PutScalar(0.0);
  return;
}


/*----------------------------------------------------------------------*
 |  set node object type                                     u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::NearestObject::setNodeObjectType(
    const int nodeId, const int label, const LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = NODE_OBJECT;
  nodeId_ = nodeId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  lineId_ = -1;
  surfId_ = -1;
}


/*----------------------------------------------------------------------*
 |  set line object type                                     u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::NearestObject::setLineObjectType(
    const int lineId, const int surfId, const int label, const LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = LINE_OBJECT;
  lineId_ = lineId;
  surfId_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  nodeId_ = -1;
}


/*----------------------------------------------------------------------*
 |  set surface object type                                  u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::NearestObject::setSurfaceObjectType(
    const int surfId, const int label, const LINALG::Matrix<3, 1>& physcoord)
{
  objectType_ = SURFACE_OBJECT;
  surfId_ = surfId;
  label_ = label;
  physcoord_ = physcoord;

  // reset unused variables
  nodeId_ = -1;
  lineId_ = -1;
}
