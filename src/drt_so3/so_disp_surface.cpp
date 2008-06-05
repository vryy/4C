/*!----------------------------------------------------------------------
\file so_disp_surface.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Element* parent,
                              const int lsurface) :
DRT::Element(id,element_sodispsurface,owner),
parent_(parent),
lsurface_(lsurface),
data_()
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)  for sending SoDispSurfaces                 umay 10/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(int id, int owner) :
DRT::Element(id,element_sodispsurface,owner),
parent_(),
lsurface_(),
data_()
{

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(const DRT::ELEMENTS::SoDispSurface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoDispSurface::Clone() const
{
  DRT::ELEMENTS::SoDispSurface* newelement = new DRT::ELEMENTS::SoDispSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::SoDispSurface::Shape() const
{
  switch (NumNode())
  {
  case 3: return tri3;
  case 4: return quad4;
  case 6: return tri6;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Pack(vector<char>& data) const
{
  data.resize(0);
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::~SoDispSurface()
{
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::SoDispSurface::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements
  return DRT::UTILS::ElementBoundaryFactory<SoDispLine,SoDispSurface>(DRT::UTILS::buildLines,this);
}
  

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Print(ostream& os) const
{
  os << "SoDispSurface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
