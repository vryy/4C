/*!----------------------------------------------------------------------
\file wall1.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::Wall1(int id, int owner) :
DRT::Element(id,element_wall1,owner),
data_(),
material_(0),
thickness_(0.0)
{
  ngp_[0] = ngp_[1] = 0;
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::Wall1(const DRT::Elements::Wall1& old) :
DRT::Element(old),
data_(old.data_),
material_(old.material_),
thickness_(old.thickness_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  for (int i=0; i<2; ++i) ngp_[i] = old.ngp_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Wall1::Clone() const
{
  DRT::Elements::Wall1* newelement = new DRT::Elements::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Wall1::Shape() const
{
  switch (NumNode())
  {
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  case 3: return tri3;
  case 6: return tri6;

  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // material_
  AddtoPack(data,material_);
  //thickness
  AddtoPack(data,thickness_);
  // ngp_
  AddtoPack(data,ngp_,2*sizeof(int));
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Unpack(const vector<char>& data)
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
  // material_
  ExtractfromPack(position,data,material_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,2*sizeof(int));
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1::~Wall1()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::Print(ostream& os) const
{
  os << "Wall1 ";
  Element::Print(os);
  os << " ngp_: " << ngp_[0] << " " << ngp_[1] << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Wall1Register (public)              mgit 03/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Wall1::ElementRegister() const
{
  return rcp(new DRT::Elements::Wall1Register(Type()));
}



/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 07/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Wall1::Lines()
{
    const DiscretizationType distype = Shape();
    const int nline   = NumLine();
    lines_.resize(nline);
    lineptrs_.resize(nline);

    switch (distype)
    {
    case tri3:
        CreateLinesTri(nline, 2);
        break;
    case tri6:
        CreateLinesTri(nline, 3);
        break;
    case quad4:
        CreateLinesQuad(nline, 2);
        break;
    case quad8:
        CreateLinesQuad(nline, 3);
        break;
    case quad9:
        CreateLinesQuad(nline, 3);
        break;
    default:
        dserror("distype not supported");
    }

    return (DRT::Element**)(&(lineptrs_[0]));
}

void DRT::Elements::Wall1::CreateLinesTri(const int& nline,
                                         const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::Utils::eleNodeNumbering_tri6_lines[iline][inode]];
             nodes[inode]   = Nodes()[DRT::Utils::eleNodeNumbering_tri6_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Wall1Line(iline,Owner(),nnode,nodeids,nodes,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}        

void DRT::Elements::Wall1::CreateLinesQuad(const int& nline,
                                          const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];	
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::Utils::eleNodeNumbering_quad9_lines[iline][inode]];
             nodes[inode]   = Nodes()[DRT::Utils::eleNodeNumbering_quad9_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Wall1Line(iline,Owner(),nnode,nodeids,nodes,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}    


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Wall1::Surfaces()
{
  surfaces_.resize(1);
  surfaces_[0] = this;
  return &surfaces_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::Wall1Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::Wall1Register(
                               const DRT::Elements::Wall1Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register* DRT::Elements::Wall1Register::Clone() const
{
  return new DRT::Elements::Wall1Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Pack(vector<char>& data) const
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
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Register::~Wall1Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Register::Print(ostream& os) const
{
  os << "Wall1Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::Elements::Wall1Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
