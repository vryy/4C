#ifdef D_ALE
#ifdef CCADISCRET

#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


using namespace DRT::Utils;

DRT::Elements::Ale3::Ale3(int id, int owner)
  : DRT::Element(id,element_ale3,owner),
    data_()
{
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
}


DRT::Elements::Ale3::Ale3(const DRT::Elements::Ale3& old)
  : DRT::Element(old),
    data_(old.data_),
    surfaces_(old.surfaces_),
    surfaceptrs_(old.surfaceptrs_)
{
  return;
}


DRT::Element* DRT::Elements::Ale3::Clone() const
{
  DRT::Elements::Ale3* newelement = new DRT::Elements::Ale3(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::Elements::Ale3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


void DRT::Elements::Ale3::Pack(vector<char>& data) const
{
  data.resize(0);
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // Gaussrule
  //AddtoPack(data,gaussrule_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);
}


void DRT::Elements::Ale3::Unpack(const vector<char>& data)
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
  // Gaussrule
  //ExtractfromPack(position,data,gaussrule_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


DRT::Elements::Ale3::~Ale3()
{
}


void DRT::Elements::Ale3::Print(ostream& os) const
{
  os << "Ale3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


RefCountPtr<DRT::ElementRegister> DRT::Elements::Ale3::ElementRegister() const
{
  return rcp(new DRT::Elements::Ale3Register(Type()));
}


//
// get vector of surfaces
//
DRT::Element** DRT::Elements::Ale3::Surfaces()
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
      CreateSurfaceWedge(nsurf, 6);
        break;
    case wedge15:
      CreateSurfaceWedge(nsurf, 15);
        break;
    case pyramid5:
      CreateSurfacesPyramid();
      break;
    default:
        dserror("distype not supported");
    }
    return (DRT::Element**)(&(surfaceptrs_[0]));
}



// support for above
void DRT::Elements::Ale3::CreateSurfacesTet(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::Ale3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}


// support for above
void DRT::Elements::Ale3::CreateSurfacesHex(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::Ale3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}

void DRT::Elements::Ale3::CreateSurfaceWedge(const int& nsurf, const int& wedgetype)
{
  int nnode, trisurfacenodes, quadsurfacenodes;

  switch (wedgetype){
  case 6:
    trisurfacenodes=3;
    quadsurfacenodes=4;
    for (int isurf=0; isurf<nsurf; isurf++)
    {
      int nodeids[nnode];
      DRT::Node* nodes[nnode];

      if (isurf <=2)
        nnode=trisurfacenodes;
      else
        nnode=quadsurfacenodes;


      for (int inode=0; inode<nnode;inode++)
      {
        if (isurf < 2){
        nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
      nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
    }
    else
      nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
      nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
  }
       surfaces_[isurf] = rcp(new DRT::Elements::Ale3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
    break;

  case 15:
    trisurfacenodes=6;
    quadsurfacenodes=8;

    for (int isurf=0; isurf<nsurf; isurf++)
    {

      int nodeids[nnode];
      DRT::Node* nodes[nnode];

      if (isurf <=2)
        nnode=trisurfacenodes;
      else
        nnode=quadsurfacenodes;


      for (int inode=0; inode<nnode;inode++)
      {
        if (isurf < 2){
        nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
      nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
    }
    else
      nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
      nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
  }
       surfaces_[isurf] = rcp(new DRT::Elements::Ale3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
    break;

}
}

void DRT::Elements::Ale3::CreateSurfacesPyramid()
{
  // Quad surface

        const int nnode_surf = 4;
        const int surfid = 0;
        int nodeids[nnode_surf];
        DRT::Node* nodes[nnode_surf];
        for (int qinode = 0; qinode < nnode_surf; qinode++) {
          nodeids[qinode] = NodeIds()[eleNodeNumbering_hex27_surfaces[surfid][qinode]];
          nodes[qinode] = Nodes()[eleNodeNumbering_hex27_surfaces[surfid][qinode]];
        }
        surfaces_[surfid] = rcp(new DRT::Elements::Ale3Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[surfid] = surfaces_[surfid].get();

  // tri surfaces
  for (int tisurf = 0; tisurf < 4; tisurf++)
  {
      const int nnode_surf = 3;
      const int surfid = tisurf+1;
      int nodeids[nnode_surf];
      DRT::Node* nodes[nnode_surf];
      for (int tinode = 0; tinode < nnode_surf; tinode++) {
        nodeids[tinode] = NodeIds()[eleNodeNumbering_pyramid5_trisurfaces[tisurf][tinode]];
        nodes[tinode] = Nodes()[eleNodeNumbering_pyramid5_trisurfaces[tisurf][tinode]];
      }
      surfaces_[surfid] = rcp(new DRT::Elements::Ale3Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
      surfaceptrs_[surfid] = surfaces_[surfid].get();

  }
}

DRT::Element** DRT::Elements::Ale3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Ale3 element itself
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


DRT::Elements::Ale3Register::Ale3Register(DRT::Element::ElementType etype)
  : ElementRegister(etype)
{
}


DRT::Elements::Ale3Register::Ale3Register(const DRT::Elements::Ale3Register& old)
  : ElementRegister(old)
{
}


DRT::Elements::Ale3Register* DRT::Elements::Ale3Register::Clone() const
{
  return new DRT::Elements::Ale3Register(*this);
}


void DRT::Elements::Ale3Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);
}


void DRT::Elements::Ale3Register::Unpack(const vector<char>& data)
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
}


DRT::Elements::Ale3Register::~Ale3Register()
{
}


void DRT::Elements::Ale3Register::Print(ostream& os) const
{
  os << "Ale3Register ";
  ElementRegister::Print(os);
}


#endif
#endif
