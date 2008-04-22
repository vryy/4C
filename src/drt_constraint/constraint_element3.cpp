/*!----------------------------------------------------------------------
\file constraint_element3.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraint_element3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(int id, int owner) :
DRT::Element(id,element_constraintelement3,owner),
data_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(const DRT::ELEMENTS::ConstraintElement3& old) :
DRT::Element(old),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of ConstraintElement3 and return pointer to it (public) |
 |                                                          gammi 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ConstraintElement3::Clone() const
{
  DRT::ELEMENTS::ConstraintElement3* newelement = new DRT::ELEMENTS::ConstraintElement3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::ConstraintElement3::Shape() const
{
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3::Pack(vector<char>& data) const
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
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::~ConstraintElement3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3::Print(ostream& os) const
{
  os << "ConstraintElement3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return ConstraintElement3Register (public)              gammi 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::ConstraintElement3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::ConstraintElement3Register(Type()));
}

//
// get vector of lines
//
DRT::Element** DRT::ELEMENTS::ConstraintElement3::Lines()
{
    dserror("noe, kriegst nix!");
    surface_.resize(1);
    surface_[0] = this; //points to ConstraintElement3 element itself
    return &surface_[0];
}

//void DRT::ELEMENTS::ConstraintElement3::CreateLinesTri(const int& nline,
//                                           const int& nnode)
//{
//    for(int iline=0;iline<nline;iline++)
//    {
//        vector<int> nodeids(nnode);
//        vector<DRT::Node*> nodes(nnode);
//
//        for (int inode=0;inode<nnode;inode++)
//        {
//             nodeids[inode] = NodeIds()[eleNodeNumbering_tri6_lines[iline][inode]];
//             nodes[inode]   = Nodes()[  eleNodeNumbering_tri6_lines[iline][inode]];
//        }
//        lines_[iline] = rcp(new DRT::ELEMENTS::ConstraintElement3Line(iline,Owner(),nnode,&nodeids[0],&nodes[0],this,iline));
//        lineptrs_[iline] = lines_[iline].get();
//    }
//}

//void DRT::ELEMENTS::ConstraintElement3::CreateLinesQuad(const int& nline,
//                                            const int& nnode)
//{
//    for(int iline=0;iline<nline;iline++)
//    {
//        vector<int> nodeids(nnode);
//        vector<DRT::Node*> nodes(nnode);
//
//        for (int inode=0;inode<nnode;inode++)
//        {
//             nodeids[inode] = NodeIds()[eleNodeNumbering_quad9_lines[iline][inode]];
//             nodes[inode]   = Nodes()[  eleNodeNumbering_quad9_lines[iline][inode]];
//        }
//        lines_[iline] = rcp(new DRT::ELEMENTS::ConstraintElement3Line(iline,Owner(),nnode,&nodeids[0],&nodes[0],this,iline));
//        lineptrs_[iline] = lines_[iline].get();
//    }
//}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::ConstraintElement3::Surfaces()
{
  surface_.resize(1);
  surface_[0] = this; //points to ConstraintElement3 element itself
  return &surface_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register::ConstraintElement3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register::ConstraintElement3Register(
                               const DRT::ELEMENTS::ConstraintElement3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register* DRT::ELEMENTS::ConstraintElement3Register::Clone() const
{
  return new DRT::ELEMENTS::ConstraintElement3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::ConstraintElement3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::ConstraintElement3Register::~ConstraintElement3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3Register::Print(ostream& os) const
{
  os << "ConstraintElement3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
