/*!----------------------------------------------------------------------
\file intersection_interfacepoint.cpp

\brief  InterfacePoint stores and delivers all data a point lying on the 
        intersection interface has to know 

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/


#ifdef CCADISCRET

#include "intersection_interfacepoint.H"


#ifdef PARALLEL
#include <mpi.h>
#endif



/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
XFEM::InterfacePoint::InterfacePoint():
pType_(NOTYPE),
nnode_(0),
nline_(0),
nsurf_(0)
{
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
XFEM::InterfacePoint::InterfacePoint(
  XFEM::pointType     pType,
  int                 nodeId,
  std::vector<int>    lineId,
  std::vector<int>    surfId,
  std::vector<double> coordinates
  ):
pType_(pType),
nodeId_(nodeId),
lineId_(lineId),
surfId_(surfId),
coord_(coordinates)
{
  setNodeLineSurfNumbers(pType);
  return;
}

 
/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 07/08|
 *----------------------------------------------------------------------*/
XFEM::InterfacePoint::InterfacePoint(
    const XFEM::InterfacePoint& old) :
pType_(old.pType_),
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
 |  set xfem number of nodes, lines, surfaces the interface  u.may 07/08|
 |  point is lying on according to point type                           |                                          
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setNodeLineSurfNumbers(
  const XFEM::pointType pType)
{
  switch(pType)
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
void XFEM::InterfacePoint::setPointType(
      const XFEM::pointType   pType
      )
{
  pType_ = pType;
  setNodeLineSurfNumbers(pType);
}
    

/*----------------------------------------------------------------------*
 |  set xfem node ids the interface point is lying on        u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setNodeId(
      const int    nodeId
      )
{
  if(nnode_ != 1)
    dserror("point type is not correct (nodeId)");
  
  nodeId_ = nodeId;
}
       

/*----------------------------------------------------------------------*
 |  set xfem line ids the interface point is lying on        u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setLineId(
      const std::vector<int>&    lineId
      )
{
  if(nline_ != (int) lineId.size())
    dserror("point type is not correct (lineId)");
  
  if(!lineId_.empty()) 
    lineId_.clear();
  
  lineId_ = lineId;
} 
  

/*----------------------------------------------------------------------*
 |  set xfem surface ids the interface point is lying on     u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setSurfaceId(
      const std::vector<int>&    surfId
      )
{
  if(nsurf_ != (int) surfId.size())
    dserror("point type is not correct (surfId)");
    
  if(!surfId_.empty()) 
      surfId_.clear();
  
  surfId_ = surfId;
}  
     

/*----------------------------------------------------------------------*
 |  set coordinates of the interface point                   u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setCoord(
      const std::vector<double>&    coordinates
      )
{
  if((int) coordinates.size() != 3)
    dserror("dimension of coordinates is not correct");
  
  if(!coord_.empty()) 
    coord_.clear();
  
  coord_ = coordinates;
} 



/*----------------------------------------------------------------------*
 |  set X-coordinates of the interface point                 u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setCoordX(
      const double    coordX
      )
{
  if((int) coord_.size() != 3)
    dserror("coordinates not yet initialized");
  
  coord_[0] = coordX;
} 



/*----------------------------------------------------------------------*
 |  set Y-coordinates of the interface point                 u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setCoordY(
      const double    coordY
      )
{
  if((int) coord_.size() != 3)
    dserror("coordinates not yet initialized");
  
  coord_[1] = coordY;
} 



/*----------------------------------------------------------------------*
 |  set Z-coordinates of the interface point                 u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setCoordZ(
      const double    coordZ
      )
{
  if((int) coord_.size() != 3)
    dserror("coordinates not yet initialized");
  
  coord_[2] = coordZ;
} 


/*----------------------------------------------------------------------*
 |  set single coordinates of the interface point            u.may 07/08|                                         
 *----------------------------------------------------------------------*/
void XFEM::InterfacePoint::setSingleCoord(
      const int	      index,
      const double    coord
      )
{
  if((int) coord_.size() != 3)
    dserror("coordinates not yet initialized");
  if(index > 2)
    dserror("index out of range");

  coord_[index] = coord;
} 



#endif  // #ifdef CCADISCRET

