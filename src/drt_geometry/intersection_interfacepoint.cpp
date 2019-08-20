/*----------------------------------------------------------------------*/
/*! \file

\brief  InterfacePoint stores and delivers all data a point lying on the
        intersection interface has to know

\level 2

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#include "intersection_interfacepoint.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::InterfacePoint::InterfacePoint() : pType_(NOTYPE), nnode_(0), nline_(0), nsurf_(0)
{
  nodeId_ = -1;
  lineId_.clear();
  surfId_.clear();
  coord_.Clear();
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::InterfacePoint::InterfacePoint(GEO::pointType& pType, int nodeId, std::vector<int>& lineId,
    std::vector<int>& surfId, LINALG::Matrix<3, 1>& coordinates)
    : pType_(pType), nodeId_(nodeId), lineId_(lineId), surfId_(surfId), coord_(coordinates)
{
  setNodeLineSurfNumbers(pType);
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::InterfacePoint::InterfacePoint(const GEO::InterfacePoint& old)
    : pType_(old.pType_),
      nnode_(old.nnode_),
      nline_(old.nline_),
      nsurf_(old.nsurf_),
      nodeId_(old.nodeId_),
      lineId_(old.lineId_),
      surfId_(old.surfId_),
      coord_(old.coord_)
{
  return;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::InterfacePoint::~InterfacePoint() { return; }


/*----------------------------------------------------------------------*
 * assignment operator                                       u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::InterfacePoint& GEO::InterfacePoint::operator=(const GEO::InterfacePoint& point)
{
  pType_ = point.pType_;
  nnode_ = point.nnode_;
  nline_ = point.nline_;
  nsurf_ = point.nsurf_;
  nodeId_ = point.nodeId_;
  lineId_ = point.lineId_;
  surfId_ = point.surfId_;
  coord_ = point.coord_;
  return *this;
}


/*----------------------------------------------------------------------*
 |  set xfem number of nodes, lines, surfaces the interface  u.may 07/08|
 |  point is lying on according to point type                           |
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setNodeLineSurfNumbers(const GEO::pointType pType)
{
  switch (pType)
  {
    case INTERNAL:
    {
      nnode_ = 0;
      nline_ = 0;
      nsurf_ = 0;
      break;
    }
    case SURFACE:
    {
      nnode_ = 0;
      nline_ = 0;
      nsurf_ = 1;
      break;
    }
    case LINE:
    {
      nnode_ = 0;
      nline_ = 1;
      nsurf_ = 2;
      break;
    }
    case NODE:
    {
      nnode_ = 1;
      nline_ = 3;
      nsurf_ = 3;
      break;
    }
    default:
      dserror("incorrect point type");
  }
}


/*----------------------------------------------------------------------*
 |  set point type the interface point                       u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setPointType(const GEO::pointType pType)
{
  pType_ = pType;
  setNodeLineSurfNumbers(pType);
}


/*----------------------------------------------------------------------*
 |  set xfem node ids the interface point is lying on        u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setNodeId(const int nodeId)
{
  if (nnode_ != 1) dserror("point type is not correct (nodeId)");

  nodeId_ = nodeId;
}


/*----------------------------------------------------------------------*
 |  set xfem line ids the interface point is lying on        u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setLineId(const std::vector<int>& lineId)
{
  if (nline_ != (int)lineId.size()) dserror("point type is not correct (lineId)");

  if (!lineId_.empty()) lineId_.clear();

  lineId_ = lineId;
}


/*----------------------------------------------------------------------*
 |  set xfem surface ids the interface point is lying on     u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setSurfaceId(const std::vector<int>& surfId)
{
  if (nsurf_ != (int)surfId.size()) dserror("point type is not correct (surfId)");

  if (!surfId_.empty()) surfId_.clear();

  surfId_ = surfId;
}


/*----------------------------------------------------------------------*
 |  set coordinates of the interface point                   u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setCoord(const LINALG::Matrix<3, 1>& coordinates)
{
  if (!coordinates.IsInitialized()) dserror("dimension of coordinates is not correct");

  coord_ = coordinates;
}



/*----------------------------------------------------------------------*
 |  set X-coordinates of the interface point                 u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setCoordX(const double coordX)
{
  if (!coord_.IsInitialized()) dserror("coordinates not yet initialized");

  coord_(0) = coordX;
}



/*----------------------------------------------------------------------*
 |  set Y-coordinates of the interface point                 u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setCoordY(const double coordY)
{
  if (!coord_.IsInitialized()) dserror("coordinates not yet initialized");

  coord_(1) = coordY;
}



/*----------------------------------------------------------------------*
 |  set Z-coordinates of the interface point                 u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setCoordZ(const double coordZ)
{
  if (!coord_.IsInitialized()) dserror("coordinates not yet initialized");

  coord_(2) = coordZ;
}


/*----------------------------------------------------------------------*
 |  set single coordinates of the interface point            u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::InterfacePoint::setSingleCoord(const int index, const double coord)
{
  if (!coord_.IsInitialized()) dserror("coordinates not yet initialized");
  if (index > 2) dserror("index out of range");

  coord_(index) = coord;
}
