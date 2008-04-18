/*!----------------------------------------------------------------------
\file bele3.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "bele3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(int id, int owner) :
DRT::Element(id,element_fluid2,owner),
data_()
{
  gaussrule_ = intrule2D_undefined;
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::Bele3(const DRT::ELEMENTS::Bele3& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
data_(old.data_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  gaussrule_ = old.gaussrule_;
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Bele3 and return pointer to it (public) |
 |                                                          gammi 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele3::Clone() const
{
  DRT::ELEMENTS::Bele3* newelement = new DRT::ELEMENTS::Bele3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele3::Shape() const
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
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Pack(vector<char>& data) const
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

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Unpack(const vector<char>& data)
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
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = GaussRule2D(gausrule_integer); //explicit conversion from integer to enum

  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3::~Bele3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3::Print(ostream& os) const
{
  os << "Bele3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Bele3Register (public)              gammi 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Bele3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Bele3Register(Type()));
}

//
// get vector of lines
//
DRT::Element** DRT::ELEMENTS::Bele3::Lines()
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

void DRT::ELEMENTS::Bele3::CreateLinesTri(const int& nline,
                                           const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        vector<int> nodeids(nnode);
        vector<DRT::Node*> nodes(nnode);

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_tri6_lines[iline][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_tri6_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::ELEMENTS::Bele3Line(iline,Owner(),nnode,&nodeids[0],&nodes[0],this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}

void DRT::ELEMENTS::Bele3::CreateLinesQuad(const int& nline,
                                            const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        vector<int> nodeids(nnode);
        vector<DRT::Node*> nodes(nnode);

        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[eleNodeNumbering_quad9_lines[iline][inode]];
             nodes[inode]   = Nodes()[  eleNodeNumbering_quad9_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::ELEMENTS::Bele3Line(iline,Owner(),nnode,&nodeids[0],&nodes[0],this,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Bele3::Surfaces()
{
  surface_.resize(1);
  surface_[0] = this; //points to Bele3 element itself
  return &surface_[0];
}


GaussRule2D DRT::ELEMENTS::Bele3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule2D rule = intrule2D_undefined;
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
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register::Bele3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register::Bele3Register(
                               const DRT::ELEMENTS::Bele3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Register* DRT::ELEMENTS::Bele3Register::Clone() const
{
  return new DRT::ELEMENTS::Bele3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Bele3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Bele3Register::~Bele3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Register::Print(ostream& os) const
{
  os << "Bele3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
