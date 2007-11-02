/*!----------------------------------------------------------------------
\file so_disp.cpp
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

#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

using namespace DRT::Utils;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::SoDisp::SoDisp(int id, int owner) :
DRT::Element(id,element_sodisp,owner),
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
DRT::Elements::SoDisp::SoDisp(const DRT::Elements::SoDisp& old) :
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
DRT::Element* DRT::Elements::SoDisp::Clone() const
{
  DRT::Elements::SoDisp* newelement = new DRT::Elements::SoDisp(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::SoDisp::Shape() const
{
    switch (NumNode())
    {
    case  4: return tet4;
    case  8: return hex8;
    case 10: return tet10;
    case 20: return hex20;
    case 27: return hex27;
    case  6: return wedge6;
    case 15: return wedge15;
    default:
      dserror("unexpected number of nodes %d", NumNode());
    }
    return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDisp::Pack(vector<char>& data) const
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
void DRT::Elements::SoDisp::Unpack(const vector<char>& data)
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
DRT::Elements::SoDisp::~SoDisp()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDisp::Print(ostream& os) const
{
  os << "SoDisp ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return SoDispRegister (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::SoDisp::ElementRegister() const
{
  return rcp(new DRT::Elements::SoDispRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::SoDisp::Volumes()
{
  volume_.resize(1);
  return 0;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::SoDisp::Surfaces()
{
  
    const DiscretizationType distype = Shape();
    const int nsurf = NumSurface();
    surfaces_.resize(nsurf);
    surfaceptrs_.resize(nsurf);

    switch (distype)
    {
    case tet4:
        CreateSurfacesTet(nsurf, 3);
        break;
    case tet10:
        CreateSurfacesTet(nsurf, 6);
        break;
    case hex8:
        CreateSurfacesHex(nsurf, 4);
        break;
    case hex20:
        CreateSurfacesHex(nsurf, 8);
        break;
    case hex27:
        CreateSurfacesHex(nsurf, 9);
        break;
    case wedge6:
        CreateSurfacesWegde6(nsurf);
        break;
    case wedge15:
        CreateSurfacesWegde15(nsurf);
        break;
    default:
        dserror("distype not supported");
    }
    return (DRT::Element**)(&(surfaceptrs_[0]));
}

// support for above
void DRT::Elements::SoDisp::CreateSurfacesTet(const int& nsurf,
                                              const int& nnode)
{
    for(int isurf=0;isurf<nsurf;isurf++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_tet10_surfaces[isurf][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_tet10_surfaces[isurf][inode]];
        }
        surfaces_[isurf] = rcp(new DRT::Elements::SoDispSurface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}


// support for above
void DRT::Elements::SoDisp::CreateSurfacesHex(const int& nsurf,
                                               const int& nnode)
{
    for(int isurf=0;isurf<nsurf;isurf++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_hex27_surfaces[isurf][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_hex27_surfaces[isurf][inode]];
        }
        surfaces_[isurf] = rcp(new DRT::Elements::SoDispSurface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}

// support for above
void DRT::Elements::SoDisp::CreateSurfacesWegde6(const int& nsurf)
{
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
          surfaces_[qisurf] = rcp(new DRT::Elements::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
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
        surfaces_[surfid] = rcp(new DRT::Elements::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[surfid] = surfaces_[surfid].get();
    };
}

// support for above
void DRT::Elements::SoDisp::CreateSurfacesWegde15(const int& nsurf)
{
    // first the 3 quad surfaces (#0..2)
    for (int qisurf = 0; qisurf < 3; ++qisurf) 
    {
      const int nnode_surf = 8;
      int nodeids[nnode_surf];
      DRT::Node* nodes[nnode_surf];
      for (int qinode = 0; qinode < nnode_surf; ++qinode) {
        nodeids[qinode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[qisurf][qinode]];
        nodes[qinode] = Nodes()[eleNodeNumbering_wedge15_quadsurfaces[qisurf][qinode]];
      }
      surfaces_[qisurf] = rcp(new DRT::Elements::SoDispSurface(qisurf,Owner(),nnode_surf,nodeids,nodes,this,qisurf));
      surfaceptrs_[qisurf] = surfaces_[qisurf].get();
    };
    // then the tri's...
    for (int tisurf = 0; tisurf < 2; ++tisurf) 
    {
        const int nnode_surf = 6;
        const int surfid = tisurf + 3;
        int nodeids[nnode_surf];
        DRT::Node* nodes[nnode_surf];
        for (int tinode = 0; tinode < nnode_surf; ++tinode) {
          nodeids[tinode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[tisurf][tinode]];
          nodes[tinode] = Nodes()[eleNodeNumbering_wedge15_trisurfaces[tisurf][tinode]];
        }
        surfaces_[surfid] = rcp(new DRT::Elements::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[surfid] = surfaces_[surfid].get();
    };
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::SoDisp::Lines()
{
  dserror("SoDisp lines not yet implemented");
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
DRT::Elements::SoDispRegister::SoDispRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispRegister::SoDispRegister(
                               const DRT::Elements::SoDispRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispRegister* DRT::Elements::SoDispRegister::Clone() const
{
  return new DRT::Elements::SoDispRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispRegister::Pack(vector<char>& data) const
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
void DRT::Elements::SoDispRegister::Unpack(const vector<char>& data)
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
DRT::Elements::SoDispRegister::~SoDispRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispRegister::Print(ostream& os) const
{
  os << "SoDispRegister ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
