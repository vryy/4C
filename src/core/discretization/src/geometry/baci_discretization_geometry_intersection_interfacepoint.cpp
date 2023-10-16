/*----------------------------------------------------------------------*/
/*! \file

\brief  InterfacePoint stores and delivers all data a point lying on the
        intersection interface has to know

\level 2

*----------------------------------------------------------------------*/


#include "baci_discretization_geometry_intersection_interfacepoint.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
CORE::GEO::InterfacePoint::InterfacePoint() : pType_(NOTYPE), nnode_(0), nline_(0), nsurf_(0)
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
CORE::GEO::InterfacePoint::InterfacePoint(CORE::GEO::pointType& pType, int nodeId,
    std::vector<int>& lineId, std::vector<int>& surfId, CORE::LINALG::Matrix<3, 1>& coordinates)
    : pType_(pType), nodeId_(nodeId), lineId_(lineId), surfId_(surfId), coord_(coordinates)
{
  setNodeLineSurfNumbers(pType);
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
CORE::GEO::InterfacePoint::InterfacePoint(const CORE::GEO::InterfacePoint& old)
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
CORE::GEO::InterfacePoint::~InterfacePoint() { return; }


/*----------------------------------------------------------------------*
 * assignment operator                                       u.may 07/08|
 *----------------------------------------------------------------------*/
CORE::GEO::InterfacePoint& CORE::GEO::InterfacePoint::operator=(
    const CORE::GEO::InterfacePoint& point)
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
void CORE::GEO::InterfacePoint::setNodeLineSurfNumbers(const CORE::GEO::pointType pType)
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
