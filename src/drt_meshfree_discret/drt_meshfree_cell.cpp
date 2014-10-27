/*!--------------------------------------------------------------------------
\file drt_meshfree_cell.cpp
\brief

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*--------------------------------------------------------------------------*/

#include "drt_meshfree_cell.H"
#include "drt_meshfree_node.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_node.H"

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::Cell(int id, int owner)
  : DRT::Element(id,owner),  // necessary due to virtual inheritance from DRT::Element
    DRT::MESHFREE::MeshfreeBin(id,owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::Cell(const DRT::MESHFREE::Cell& old)
  : DRT::Element(old),  // necessary due to virtual inheritance from DRT::Element
    DRT::MESHFREE::MeshfreeBin(old),
    pointid_(old.pointid_),
    point_(old.point_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::~Cell()
{
  return;
}


/*--------------------------------------------------------------------------*
 |  << operator                                                   nis Jan12 |
 *--------------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::Cell& cell)
{
  cell.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Print(std::ostream& os) const
{
  os << std::setw(6) << Id() << " Owner " << std::setw(3) << Owner() << " ";

  const int npoint = NumPoint();
  const int* pointids = PointIds();
  if (npoint > 0)
  {
    os << " Points ";
    for (int i=0; i<npoint; ++i) os << std::setw(6) << pointids[i] << " ";
  }

  const int nnode = NumNode();
  const int* nodeids = NodeIds();
  if (nnode > 0)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << std::setw(6) << nodeids[i] << " ";
  }

  return;
}

/*--------------------------------------------------------------------------*
 |  set point numbers to element                          (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::SetPointIds(const int npoint, const int* points)
{
  pointid_.resize(npoint);
  for (int i=0; i<npoint; ++i) pointid_[i] = points[i];
  point_.resize(0);
  return;
}

/*--------------------------------------------------------------------------*
 |  set point numbers to element                          (public) nis Jan12 |
 *--------------------------------------------------------------------------*/

void DRT::MESHFREE::Cell::SetPointIds(const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  linedef->ExtractIntVector(distype,pointid_);
  for (unsigned i=0; i<pointid_.size(); ++i)
    pointid_[i] -= 1;
  point_.resize(0);
}

/*----------------------------------------------------------------------*
 |  Build point pointers (protected)                           nis Jan12 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::Cell::BuildPointPointers(std::map<int,Teuchos::RCP<DRT::Node> >& points)
{
  int        npoint   = NumPoint();
  const int* pointids = PointIds();
  point_.resize(npoint);
  for (int i=0; i<npoint; ++i)
  {
    std::map<int,Teuchos::RCP<DRT::Node> >::const_iterator curr = points.find(pointids[i]);
    // this point is not on this proc
    if (curr==points.end()) dserror("Meshfree cell %d cannot find point %d",Id(),pointids[i]);
    else
      point_[i] = curr->second.get();
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Build point pointers (protected)                           nis Jan12 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::Cell::BuildPointPointers(DRT::Node** points)
{
  point_.resize(NumPoint());
  for (int i=0; i<NumPoint(); ++i) point_[i] = points[i];
  return true;
}
/*----------------------------------------------------------------------*
 |  Pack data  (public)                                       nis Jan12 |
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::MESHFREE::MeshfreeBin
  DRT::MESHFREE::MeshfreeBin::Pack(data);
  // add vector pointid_
  AddtoPack(data,pointid_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data  (public)                                     nis Jan12 |
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class DRT::MESHFREE::MeshfreeBin
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::MESHFREE::MeshfreeBin::Unpack(basedata);
  // extract pointid_
  ExtractfromPack(position,data,pointid_);
  // point_ is NOT communicated
  point_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  Coordinates of cell center computed by point position      nis Jan14 |
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::CenterAndMaxRadius(
  LINALG::Matrix<3,1>&  center,
  double& max_radius)
{
  // make sure center is clear
  center.Clear();
  max_radius = 0.0;

  // compute cell center
  const int numpoints = pointid_.size();
  for (int i=0; i<numpoints; ++i)
  {
    const LINALG::Matrix<3,1> cpoint_xyz(const_cast<double*>(point_[i]->X()),true);
    center.Update(1.0,cpoint_xyz,1.0);
  }
  center.Scale(1.0/numpoints);

  // compute maximum radius
  LINALG::Matrix<3,1> diff;
  for (int i=0; i<numpoints; ++i)
  {
    const LINALG::Matrix<3,1> cpoint_xyz(const_cast<double*>(point_[i]->X()),true);
    diff.Update(1.0,center,1.0,cpoint_xyz);
    const double dist = diff.Norm2();
    if (dist>max_radius)
      max_radius = dist;
  }

  return;
}
