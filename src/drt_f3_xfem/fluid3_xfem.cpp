/*!----------------------------------------------------------------------
\file fluid3_xfem.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3::XFluid3(int id, int owner) :
DRT::Element(id,element_xfluid3,owner),
material_(0),
is_ale_(false),
data_()
{
    gaussrule_ = hex_27point;
    surfaces_.resize(0);
    surfaceptrs_.resize(0);
    lines_.resize(0);
    lineptrs_.resize(0);
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3::XFluid3(const DRT::Elements::XFluid3& old) :
DRT::Element(old),
material_(old.material_),
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
DRT::Element* DRT::Elements::XFluid3::Clone() const
{
    DRT::Elements::XFluid3* newelement = new DRT::Elements::XFluid3(*this);
    return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::XFluid3::Shape() const
{
    switch (NumNode())
    {
    case  4: return tet4;
    case  8: return hex8;
    case 10: return tet10;
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
void DRT::Elements::XFluid3::Pack(vector<char>& data) const
{
    data.resize(0);

    // pack type of this instance of ParObject
    int type = UniqueParObjectId();
    AddtoPack(data,type);
    // add base class Element
    vector<char> basedata(0);
    Element::Pack(basedata);
    AddtoPack(data,basedata);
    // ngp_
    AddtoPack(data,gaussrule_);
    // material_
    AddtoPack(data,material_);
    // is_ale_
    AddtoPack(data,is_ale_);
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
void DRT::Elements::XFluid3::Unpack(const vector<char>& data)
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
  
    ExtractfromPack(position,data,gaussrule_);
    ExtractfromPack(position,data,material_);
    ExtractfromPack(position,data,is_ale_);
    // data_
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
DRT::Elements::XFluid3::~XFluid3()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::Print(ostream& os) const
{
    os << "XFluid3 ";
    Element::Print(os);
    cout << endl;
    cout << data_;
    return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::XFluid3::ElementRegister() const
{
    return rcp(new DRT::Elements::XFluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/ 
DRT::Element** DRT::Elements::XFluid3::Lines()
{
    const DiscretizationType distype = Shape();
	const int nline = NumLine();
	lines_.resize(nline);
	lineptrs_.resize(nline);
	
    switch (distype)
    {
    case tet4:
        CreateLinesTet(6, 2);
        break;
    case tet10:
        CreateLinesTet(6, 3);
        break;
    case hex8:
        CreateLinesHex(12, 2);
        break;
	case hex20: case hex27:
		CreateLinesHex(12, 3);
        break;
    default:
        dserror("distype not supported");
    }
  	return (DRT::Element**)(&(lineptrs_[0]));
}

// support for above
void DRT::Elements::XFluid3::CreateLinesTet(const int& nline,
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
        lines_[iline] = rcp(new DRT::Elements::XFluid3Line(iline,Owner(),nnode,nodeids,nodes,NULL,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}        


// support for above
void DRT::Elements::XFluid3::CreateLinesHex(const int& nline,
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
        lines_[iline] = rcp(new DRT::Elements::XFluid3Line(iline,Owner(),nnode,nodeids,nodes,NULL,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}   

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::XFluid3::Surfaces()
{
    const DiscretizationType distype = Shape(); 
    const int nsurf = NumSurface();
  	surfaces_.resize(nsurf);
  	surfaceptrs_.resize(nsurf);
    
    switch (distype)
    {
    case tet4:
        CreateSurfacesTet(4, 3);
        break;
    case tet10:
        CreateSurfacesTet(4, 6);
        break;
  	case hex8:
        CreateSurfacesHex(6, 4);
        break;
    case hex20: case hex27:
        CreateSurfacesHex(6, 9);
        break;
  	default: 
        dserror("distype not supported");
    }
  	return (DRT::Element**)(&(surfaceptrs_[0]));
}



// support for above
void DRT::Elements::XFluid3::CreateSurfacesTet(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}        


// support for above
void DRT::Elements::XFluid3::CreateSurfacesHex(const int& nsurf,
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
        surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),nnode,nodeids,nodes,this,isurf));
        surfaceptrs_[isurf] = surfaces_[isurf].get();
    }
}   


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::XFluid3::Volumes()
{
    volume_.resize(1);
    volume_[0] = this;     //points to XFluid3 element itself
    return &volume_[0];
}
 
   


  

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register::XFluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
    return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register::XFluid3Register(
                               const DRT::Elements::XFluid3Register& old) :
ElementRegister(old)
{
    return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Register* DRT::Elements::XFluid3Register::Clone() const
{
  return new DRT::Elements::XFluid3Register(*this);
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Register::Pack(vector<char>& data) const
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
void DRT::Elements::XFluid3Register::Unpack(const vector<char>& data)
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
DRT::Elements::XFluid3Register::~XFluid3Register()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Register::Print(ostream& os) const
{
    os << "XFluid3Register ";
    ElementRegister::Print(os);
    return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
