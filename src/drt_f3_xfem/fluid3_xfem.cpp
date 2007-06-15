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
	int nodeids[3];
	DRT::Node* nodes[3];
	
    switch (distype)
    {
    case tet4:
        for(int iline=0;iline<nline;iline++)
        {
            for (int inode=0;inode<2;inode++)
            {
                nodeids[inode] = NodeIds()[tet10_lines[iline][inode]];
                nodes[inode] = Nodes()[tet10_lines[iline][inode]];
            }
            lines_[iline] = rcp(new DRT::Elements::XFluid3Line(iline,Owner(),2,nodeids,nodes,NULL,this,iline));
            lineptrs_[iline] = lines_[iline].get();
        }
        break;
    case tet10:
        for(int iline=0;iline<nline;iline++)
        {
            for (int inode=0;inode<3;inode++)
            {
                nodeids[inode] = NodeIds()[tet10_lines[iline][inode]];
                nodes[inode] = Nodes()[tet10_lines[iline][inode]];
            }
            lines_[iline] = rcp(new DRT::Elements::XFluid3Line(iline,Owner(),3,nodeids,nodes,NULL,this,iline));
            lineptrs_[iline] = lines_[iline].get();
        }
        break;
    case hex8:
        for(int iline=0;iline<nline;iline++)
        {
	  	    nodeids[0] = NodeIds()[hex27_lines[iline][0]];
	  	    nodeids[1] = NodeIds()[hex27_lines[iline][1]];
	  	    nodes[0] = Nodes()[0];
	  	    nodes[1] = Nodes()[1];
	  	    lines_[iline] = rcp(new DRT::Elements::XFluid3Line(iline,Owner(),2,nodeids,nodes,NULL,this,iline));
	  	    lineptrs_[iline] = lines_[iline].get();
        }    
        break;
	case hex27:
		for(int iline=0;iline<nline;iline++)
        {
    		nodeids[0] = NodeIds()[hex27_lines[iline][0]];
    	  	nodeids[1] = NodeIds()[hex27_lines[iline][1]];
    	  	nodeids[2] = NodeIds()[hex27_lines[iline][2]];
    	  	nodes[0] = Nodes()[0];
    	  	nodes[1] = Nodes()[1];
    	  	nodes[2] = Nodes()[2];
	  	    lines_[iline] =	rcp(new DRT::Elements::XFluid3Line(iline,Owner(),3,nodeids,nodes,NULL,this,iline));
	  	    lineptrs_[iline] = lines_[iline].get();
        }
        break;
    default:
        dserror("distype not supported");
    }
  	return (DRT::Element**)(&(lineptrs_[0]));
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
  	int nodeids[9];
  	DRT::Node* nodes[9];
    
    switch (distype)
    {
    case tet4:
        for (int isurf=0;isurf<nsurf;isurf++)
        {
            for (int inode=0;inode<3;inode++)
            {
                nodeids[inode] = NodeIds()[tet10_surfaces[isurf][inode]];
                nodes[inode] = Nodes()[tet10_surfaces[isurf][inode]];
            }
            surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),3,nodeids,nodes,this,isurf));
            surfaceptrs_[isurf] = surfaces_[isurf].get();
        }
        break;
    case tet10:
        for (int isurf=0;isurf<nsurf;isurf++)
        {
            for (int inode=0;inode<6;inode++)
            {
                nodeids[inode] = NodeIds()[tet10_surfaces[isurf][inode]];
                nodes[inode] = Nodes()[tet10_surfaces[isurf][inode]];
            }
            surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),3,nodeids,nodes,this,isurf));
            surfaceptrs_[isurf] = surfaces_[isurf].get();
        }
        break;
  	case hex8:
    	for (int isurf=0;isurf<nsurf;isurf++)
    	{
            for (int inode=0;inode<4;inode++)
            {
                nodeids[inode] = NodeIds()[hex27_surfaces[isurf][inode]];
                nodes[inode] = Nodes()[hex27_surfaces[isurf][inode]];
            }
            surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),4,nodeids,nodes,this,isurf));
            surfaceptrs_[isurf] = surfaces_[isurf].get();
        }
        break;
    case hex27:
        for (int isurf=0;isurf<nsurf;isurf++)
        {
            for (int inode=0;inode<9;inode++)
            {
                nodeids[inode] = NodeIds()[hex27_surfaces[isurf][inode]];
                nodes[inode] = Nodes()[hex27_surfaces[isurf][inode]];
            }
            surfaces_[isurf] = rcp(new DRT::Elements::XFluid3Surface(isurf,Owner(),9,nodeids,nodes,this,isurf));
            surfaceptrs_[isurf] = surfaces_[isurf].get();
	
  	     }
         break;
  	default: 
        dserror("distype not supported");
    }

  	return (DRT::Element**)(&(surfaceptrs_[0]));
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

// update!!
const int DRT::Elements::XFluid3::hex27_surfaces[6][9] = {{0,  3,  2,  1, 11, 10,  9,  8, 20},
                                                          {4,  5,  6,  7, 12, 13, 14, 15, 21},
                                                          {0,  1,  5,  4,  8, 17, 12, 16, 22},
                                                          {2,  3,  7,  6, 10, 19, 14, 18, 23},
                                                          {0,  4,  7,  3, 16, 15, 19, 11, 24},
                                                          {1,  2,  6,  5,  9, 18, 13, 17, 25}};
// update!!
const int DRT::Elements::XFluid3::hex27_lines[12][3] = {{0,  1,  8},
                                                        {1,  2,  9},
                                                        {2,  3, 10},
                                                        {3,  0, 11},
                                                        {0,  4, 16},
                                                        {1,  5, 17},
                                                        {2,  6, 18},
                                                        {3,  7, 19},
                                                        {4,  5, 12},
                                                        {5,  6, 13},
                                                        {6,  7, 14},
                                                        {7,  4, 15}};

// update!!
const int DRT::Elements::XFluid3::tet10_surfaces[4][6] = {{0,  1,  8,  0,  1,  8},
                                                          {1,  2,  9,  0,  1,  8},
                                                          {2,  3, 10,  0,  1,  8},
                                                          {3,  0, 11,  0,  1,  8}};

// update!!
const int DRT::Elements::XFluid3::tet10_lines[6][3] = {{0,  1,  8},
                                                       {1,  2,  9},
                                                       {2,  3, 10},
                                                       {3,  0, 11},
                                                       {0,  4, 16},
                                                       {1,  5, 17}};

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
