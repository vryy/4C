/*!----------------------------------------------------------------------
\file so_weg6.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET

#include "so_weg6.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

using namespace DRT::Utils;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::So_weg6::So_weg6(int id, int owner) :
DRT::Element(id,element_so_weg6,owner),
data_()
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::So_weg6::So_weg6(const DRT::Elements::So_weg6& old) :
DRT::Element(old),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::So_weg6::Clone() const
{
  DRT::Elements::So_weg6* newelement = new DRT::Elements::So_weg6(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::So_weg6::Shape() const
{
  return wedge6;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::So_weg6::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // stresstype_
  AddtoPack(data,stresstype_);
  // kintype_
  AddtoPack(data,kintype_);
  // rewind flags
  AddtoPack(data,donerewinding_);
  AddtoPack(data,rewind_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::So_weg6::Unpack(const vector<char>& data)
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
  // stresstype_
  ExtractfromPack(position,data,stresstype_);
  // kintype_
  ExtractfromPack(position,data,kintype_);
  // rewinding flags
  ExtractfromPack(position,data,donerewinding_);
  ExtractfromPack(position,data,rewind_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::So_weg6::~So_weg6()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_weg6::Print(ostream& os) const
{
  os << "So_weg6 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return So_weg6Register (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::So_weg6::ElementRegister() const
{
  return rcp(new DRT::Elements::Sow6Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::So_weg6::Volumes()
{
  volume_.resize(1);
  return 0;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::So_weg6::Surfaces()
{
  const int nsurf = NumSurface();
  surfaces_.resize(nsurf);
  surfaceptrs_.resize(nsurf);
  // first the 3 quad surfaces (#0..2)
  for (int qisurf = 0; qisurf < 3; ++qisurf) 
  {
        const int nnode_surf = 4;
        const int surfid = qisurf;
        int nodeids[nnode_surf];
        DRT::Node* nodes[nnode_surf];
        for (int qinode = 0; qinode < nnode_surf; ++qinode) {
          nodeids[qinode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[surfid][qinode]];
          nodes[qinode] = Nodes()[eleNodeNumbering_wedge15_quadsurfaces[surfid][qinode]];
        }
        surfaces_[qisurf] = rcp(new DRT::Elements::Sow6Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[qisurf] = surfaces_[qisurf].get();
  };
  // then the tri's...
  for (int tisurf = 0; tisurf < 2; ++tisurf) 
  {
      const int nnode_surf = 3;
      const int surfid = tisurf + 3;
      int nodeids[nnode_surf];
      DRT::Node* nodes[nnode_surf];
      for (int tinode = 0; tinode < nnode_surf; ++tinode) {
        nodeids[tinode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[surfid][tinode]];
        nodes[tinode] = Nodes()[eleNodeNumbering_wedge15_trisurfaces[surfid][tinode]];
      }
      surfaces_[surfid] = rcp(new DRT::Elements::Sow6Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
      surfaceptrs_[surfid] = surfaces_[surfid].get();
  };
  return (DRT::Element**)(&(surfaceptrs_[0]));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::So_weg6::Lines()
{
  dserror("So_weg6 lines not yet implemented");
  const int nline = NumLine();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  return (DRT::Element**)(&(lineptrs_[0]));
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Sow6Register::Sow6Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Sow6Register::Sow6Register(
                               const DRT::Elements::Sow6Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Elements::Sow6Register* DRT::Elements::Sow6Register::Clone() const
{
  return new DRT::Elements::Sow6Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Sow6Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Sow6Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Sow6Register::~Sow6Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Sow6Register::Print(ostream& os) const
{
  os << "Sow6Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
