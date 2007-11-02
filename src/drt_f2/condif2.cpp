/*!----------------------------------------------------------------------
\file condif2.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "condif2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::Utils;


/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 05/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2::Condif2(int id, int owner) :
DRT::Element(id,element_condif2,owner),
data_()
{
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          vg 05/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2::Condif2(const DRT::Elements::Condif2& old) :
DRT::Element(old),
data_(old.data_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Condif2 and return pointer to it (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Condif2::Clone() const
{
  DRT::Elements::Condif2* newelement = new DRT::Elements::Condif2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Condif2::Shape() const
{
  switch (NumNode())
  {
  case  3: return tri3;
  case  4: return quad4;
  case  6: return tri6;
  case  8: return quad8;
  case  9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::Pack(vector<char>& data) const
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
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::Unpack(const vector<char>& data)
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
 |  dtor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2::~Condif2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                 vg 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2::Print(ostream& os) const
{
  os << "Condif2 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Condif2Register (public)                vg 05/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Condif2::ElementRegister() const
{
  return rcp(new DRT::Elements::Condif2Register(Type()));
}

//
// get vector of lines
//
DRT::Element** DRT::Elements::Condif2::Lines()
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

void DRT::Elements::Condif2::CreateLinesTri(const int& nline,
                                           const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_tri6_lines[iline][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_tri6_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Condif2Line(iline,Owner(),nnode,nodeids,nodes,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}

void DRT::Elements::Condif2::CreateLinesQuad(const int& nline,
                                            const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_quad9_lines[iline][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_quad9_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::Condif2Line(iline,Owner(),nnode,nodeids,nodes,this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)                  vg 05/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Condif2::Surfaces()
{
  surface_.resize(1);
  surface_[0] = this; //points to Condif2 element itself
  return &surface_[0];
}


GaussRule2D DRT::Elements::Condif2::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule2D rule;
    switch (distype)
    {
    case quad4:
        rule = intrule_quad_4point;
        break;
    case quad8: case quad9:
        rule = intrule_quad_9point;
        break;
    case tri3:
        rule = intrule_tri_3point;
        break;
    case tri6:
        rule = intrule_tri_6point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2Register::Condif2Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          vg 05/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2Register::Condif2Register(
                               const DRT::Elements::Condif2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2Register* DRT::Elements::Condif2Register::Clone() const
{
  return new DRT::Elements::Condif2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2Register::Pack(vector<char>& data) const
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
 |                                                             vg 05/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                               vg 05/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Condif2Register::~Condif2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                              vg 05/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Condif2Register::Print(ostream& os) const
{
  os << "Condif2Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
