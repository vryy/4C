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
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(int id, int owner) :
DRT::Element(id,element_sodisp,owner),
data_()
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  donerewinding_ = false;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(const DRT::ELEMENTS::SoDisp& old) :
DRT::Element(old),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_),
donerewinding_(old.donerewinding_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoDisp::Clone() const
{
  DRT::ELEMENTS::SoDisp* newelement = new DRT::ELEMENTS::SoDisp(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::SoDisp::Shape() const
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
void DRT::ELEMENTS::SoDisp::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::SoDisp::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::SoDisp::~SoDisp()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Print(ostream& os) const
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
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::SoDisp::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::SoDispRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Volumes()
{
  volume_.resize(1);
  return 0;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Surfaces()
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
void DRT::ELEMENTS::SoDisp::CreateSurfacesTet(const int& nsurf,
                                              const int& nnode)
{
    for(int isurf=0;isurf<nsurf;isurf++)
    {
        vector<int> nodeids(nnode);
        vector<DRT::Node*> nodes(nnode);

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_tet10_surfaces[isurf][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_tet10_surfaces[isurf][inode]];
        }
        surfaces_[isurf] = rcp(new DRT::ELEMENTS::SoDispSurface(isurf,Owner(),nnode,&nodeids[0],&nodes[0],this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}


// support for above
void DRT::ELEMENTS::SoDisp::CreateSurfacesHex(const int& nsurf,
                                               const int& nnode)
{
    for(int isurf=0;isurf<nsurf;isurf++)
    {
        vector<int> nodeids(nnode);
        vector<DRT::Node*> nodes(nnode);

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_hex27_surfaces[isurf][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_hex27_surfaces[isurf][inode]];
        }
        surfaces_[isurf] = rcp(new DRT::ELEMENTS::SoDispSurface(isurf,Owner(),nnode,&nodeids[0],&nodes[0],this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}

// support for above
void DRT::ELEMENTS::SoDisp::CreateSurfacesWegde6(const int& nsurf)
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
          surfaces_[qisurf] = rcp(new DRT::ELEMENTS::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
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
        surfaces_[surfid] = rcp(new DRT::ELEMENTS::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[surfid] = surfaces_[surfid].get();
    };
}

// support for above
void DRT::ELEMENTS::SoDisp::CreateSurfacesWegde15(const int& nsurf)
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
      surfaces_[qisurf] = rcp(new DRT::ELEMENTS::SoDispSurface(qisurf,Owner(),nnode_surf,nodeids,nodes,this,qisurf));
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
        surfaces_[surfid] = rcp(new DRT::ELEMENTS::SoDispSurface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[surfid] = surfaces_[surfid].get();
    };
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDisp::Lines()
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
DRT::ELEMENTS::SoDispRegister::SoDispRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister::SoDispRegister(
                               const DRT::ELEMENTS::SoDispRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispRegister* DRT::ELEMENTS::SoDispRegister::Clone() const
{
  return new DRT::ELEMENTS::SoDispRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispRegister::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::SoDispRegister::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::SoDispRegister::~SoDispRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispRegister::Print(ostream& os) const
{
  os << "SoDispRegister ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDispRegister::Initialize(DRT::Discretization& dis)
{
  bool dofillcompleteagain = false;
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_sodisp) continue;
    DRT::ELEMENTS::SoDisp* actele = dynamic_cast<DRT::ELEMENTS::SoDisp*>(dis.lColElement(i));
    if (!actele) dserror("cast to SoDisp* failed");
    
    const DRT::Element::DiscretizationType distype = actele->Shape();
    bool possiblytorewind = false;
    switch(distype)
    {
    case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        possiblytorewind = true;
        break;
    case DRT::Element::tet4: case DRT::Element::tet10:
        possiblytorewind = true;
        break;
    case DRT::Element::wedge6: case DRT::Element::wedge15:
        possiblytorewind = true;
        break;
    case DRT::Element::pyramid5:
        possiblytorewind = true;
        break;
    default:
        dserror("invalid discretization type for fluid3");
    }
    
    if ( (possiblytorewind) && (!actele->donerewinding_) ) {
      const bool rewind = checkRewinding3D(actele);

      if (rewind) {
        if (distype==DRT::Element::tet4){
          const int iel = actele->NumNode();
          vector<int> new_nodeids(iel);
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[1] = old_nodeids[2];
          new_nodeids[2] = old_nodeids[1];
          new_nodeids[3] = old_nodeids[3];
          actele->SetNodeIds(iel, &new_nodeids[0]);
        }
        else if (distype==DRT::Element::hex8){
          const int iel = actele->NumNode();
          vector<int> new_nodeids(iel);
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[4];
          new_nodeids[1] = old_nodeids[5];
          new_nodeids[2] = old_nodeids[6];
          new_nodeids[3] = old_nodeids[7];
          new_nodeids[4] = old_nodeids[0];
          new_nodeids[5] = old_nodeids[1];
          new_nodeids[6] = old_nodeids[2];
          new_nodeids[7] = old_nodeids[3];
          actele->SetNodeIds(iel, &new_nodeids[0]);
        }
        else if (distype==DRT::Element::wedge6){
          const int iel = actele->NumNode();
          vector<int> new_nodeids(iel);
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[3];
          new_nodeids[1] = old_nodeids[4];
          new_nodeids[2] = old_nodeids[5];
          new_nodeids[3] = old_nodeids[0];
          new_nodeids[4] = old_nodeids[1];
          new_nodeids[5] = old_nodeids[2];
          actele->SetNodeIds(iel, &new_nodeids[0]);
        }
        else if (distype == DRT::Element::pyramid5){
          const int iel = actele->NumNode();
          vector<int> new_nodeids(iel);
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[1] = old_nodeids[3];
          new_nodeids[3] = old_nodeids[1];
          // the other nodes can stay the same
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[2] = old_nodeids[2];
          new_nodeids[4] = old_nodeids[4];
          actele->SetNodeIds(iel, &new_nodeids[0]);
        }
        else dserror("no rewinding scheme for this type of fluid3");
      }
      // process of rewinding done
      actele->donerewinding_ = true;
      dofillcompleteagain = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  if(dofillcompleteagain) dis.FillComplete(false,false,false);
  
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
