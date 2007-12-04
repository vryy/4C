/*!----------------------------------------------------------------------
\file fluid3.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*/
// map to convert strings tao actions (stabilisation)
/*----------------------------------------------------------------------*/
map<string,DRT::Elements::Fluid3::StabilisationAction> DRT::Elements::Fluid3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3::Fluid3(int id, int owner) :
DRT::Element(id,element_fluid3,owner),
is_ale_(false),
data_()
{
    gaussrule_ = intrule_hex_27point;
    surfaces_.resize(0);
    surfaceptrs_.resize(0);
    lines_.resize(0);
    lineptrs_.resize(0);

    sub_acc_old_.resize(0,0);
    sub_vel_.resize(0,0);
    sub_vel_old_.resize(0,0);
    sub_pre_.resize(0);
    sub_pre_old_.resize(0);

    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3::Fluid3(const DRT::Elements::Fluid3& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
is_ale_(old.is_ale_),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
    gaussrule_ = old.gaussrule_;
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid3 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Fluid3::Clone() const
{
  DRT::Elements::Fluid3* newelement = new DRT::Elements::Fluid3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Fluid3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::Pack(vector<char>& data) const
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
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  // is_ale_
  AddtoPack(data,is_ale_);
  // rewinding bools
  AddtoPack(data,rewind_);
  AddtoPack(data,donerewinding_);

  // history variables
  AddtoPack(data,sub_acc_old_.extent(blitz::firstDim));
  AddtoPack(data,sub_acc_old_.extent(blitz::secondDim));
  
  int size = sub_acc_old_.extent(blitz::firstDim)
             *sub_acc_old_.extent(blitz::secondDim)
             *sizeof(double);
  AddtoPack(data,sub_acc_old_.data(),size);
  AddtoPack(data,sub_vel_.data()    ,size);
  AddtoPack(data,sub_vel_old_.data(),size);

  size = sub_acc_old_.extent(blitz::secondDim)*sizeof(double);
  AddtoPack(data,sub_pre_.data()    ,size);
  AddtoPack(data,sub_pre_old_.data(),size);
    
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // Gaussrule
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = GaussRule3D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);
  // rewinding bools
  ExtractfromPack(position,data,rewind_); 
  ExtractfromPack(position,data,donerewinding_);


  // history variables (subscale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position,data,firstdim);
    ExtractfromPack(position,data,secondim);

    sub_acc_old_.resize(firstdim,secondim);
    sub_vel_    .resize(firstdim,secondim);
    sub_vel_old_.resize(firstdim,secondim);
 
    int size = firstdim*secondim*sizeof(double);
       
    ExtractfromPack(position,data,&(sub_acc_old_.data()[0]),size);
    ExtractfromPack(position,data,&(sub_vel_.data()[0])    ,size);
    ExtractfromPack(position,data,&(sub_vel_old_.data()[0]),size);

    sub_pre_    .resize(secondim);
    sub_pre_old_.resize(secondim);

    ExtractfromPack(position,data,&(sub_pre_.data()[0])    ,secondim*sizeof(double));
    ExtractfromPack(position,data,&(sub_pre_old_.data()[0]),secondim*sizeof(double));

    
  }
  
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3::~Fluid3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3::Print(ostream& os) const
{
  os << "Fluid3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Fluid3::ElementRegister() const
{
  return rcp(new DRT::Elements::Fluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Fluid3::Lines()
{
    const DiscretizationType distype = Shape();
    const int nline = NumLine();
    lines_.resize(nline);
    lineptrs_.resize(nline);

    switch (distype)
    {
    case tet4:
        CreateLinesTet(nline, 2);
        break;
    case tet10:
        CreateLinesTet(nline, 3);
        break;
    case hex8:
        CreateLinesHex(nline, 2);
        break;
    case hex20: case hex27:
        CreateLinesHex(nline, 3);
        break;
    default:
        dserror("distype not supported");
    }
    return (DRT::Element**)(&(lineptrs_[0]));
}

// support for above
void DRT::Elements::Fluid3::CreateLinesTet(const int& nline,
                                           const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_tet10_lines[iline][inode]];
             nodes[inode]   = Nodes()[eleNodeNumbering_tet10_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Fluid3Line(iline,Owner(),nnode,nodeids,nodes,NULL,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}


// support for above
void DRT::Elements::Fluid3::CreateLinesHex(const int& nline,
                                           const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_hex27_lines[iline][inode]];
             nodes[inode]   = Nodes()[eleNodeNumbering_hex27_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Fluid3Line(iline,Owner(),nnode,nodeids,nodes,NULL,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Fluid3::Surfaces()
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
void DRT::Elements::Fluid3::CreateSurfacesTet(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::Fluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}


// support for above
void DRT::Elements::Fluid3::CreateSurfacesHex(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::Fluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}

void DRT::Elements::Fluid3::CreateSurfaceWedge(
        const int& nsurf,
        const int& wedgetype)
{
    switch (wedgetype)
    {
    case 6:
    {
        const int trisurfacenodes=3;   // tri3
        const int quadsurfacenodes=4;  // quad4
        for (int isurf=0; isurf<nsurf; isurf++)
        {
            int nnode = 0;
            if (isurf <=2)
                nnode=trisurfacenodes;
            else
                nnode=quadsurfacenodes;

            int nodeids[nnode];
            DRT::Node* nodes[nnode];

            for (int inode=0; inode<nnode; inode++)
            {
                if (isurf < 2)
                {
                    nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
                    nodes[  inode] = Nodes()[ eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
                }
                else
                {
                    nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
                    nodes[  inode] = Nodes()[ eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
                }
                surfaces_[isurf] = rcp(new DRT::Elements::Fluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
                surfaceptrs_[isurf] = surfaces_[isurf].get();
            }
        }
        break;
    }
    case 15:
    {
        const int trisurfacenodes=6;   // tri6
        const int quadsurfacenodes=8;  // quad8

        for (int isurf=0; isurf<nsurf; isurf++)
        {
            int nnode = 0;
            if (isurf <=2)
                nnode=trisurfacenodes;
            else
                nnode=quadsurfacenodes;

            int nodeids[nnode];
            DRT::Node* nodes[nnode];
            for (int inode=0; inode<nnode; inode++)
            {
                if (isurf < 2)
                {
                    nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
                    nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_trisurfaces[isurf][inode]];
                }
                else
                {
                    nodeids[inode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
                    nodes[inode] = Nodes()[ eleNodeNumbering_wedge15_quadsurfaces[isurf][inode]];
                }
                surfaces_[isurf] = rcp(new DRT::Elements::Fluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
                surfaceptrs_[isurf] = surfaces_[isurf].get();
            }
        }
        break;
    }
    default:
        dserror("incorrect wedgetype");
    } // end switch wedge type
}

void DRT::Elements::Fluid3::CreateSurfacesPyramid()
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
        surfaces_[surfid] = rcp(new DRT::Elements::Fluid3Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
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
      surfaces_[surfid] = rcp(new DRT::Elements::Fluid3Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
      surfaceptrs_[surfid] = surfaces_[surfid].get();

  }
}

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Fluid3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Fluid3 element itself
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Register::Fluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Register::Fluid3Register(
                               const DRT::Elements::Fluid3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Register* DRT::Elements::Fluid3Register::Clone() const
{
  return new DRT::Elements::Fluid3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Register::Pack(vector<char>& data) const
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
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Register::~Fluid3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Register::Print(ostream& os) const
{
  os << "Fluid3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
