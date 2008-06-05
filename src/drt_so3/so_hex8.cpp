/*!----------------------------------------------------------------------
\file so_hex8.cpp
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

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(int id, int owner) :
DRT::Element(id,element_so_hex8,owner),
data_()
{
  kintype_ = soh8_totlag;
  eastype_ = soh8_easnone;
  neas_ = 0;
#if defined(PRESTRESS) || defined(POSTSTRESS)
  glprestrain_ = rcp(new Epetra_SerialDenseMatrix(NUMGPT_SOH8,NUMSTR_SOH8));
#endif
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (int i=0; i<NUMGPT_SOH8; ++i)
    invJ_[i].Shape(3,3);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(const DRT::ELEMENTS::So_hex8& old) :
DRT::Element(old),
kintype_(old.kintype_),
eastype_(old.eastype_),
neas_(old.neas_),
data_(old.data_),
detJ_(old.detJ_)
{
#if defined(PRESTRESS) || defined(POSTSTRESS)
  glprestrain_ = rcp(new Epetra_SerialDenseMatrix(*(old.glprestrain_)));
#endif
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
    invJ_[i] = old.invJ_[i];
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex8::Clone() const
{
  DRT::ELEMENTS::So_hex8* newelement = new DRT::ELEMENTS::So_hex8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_hex8::Shape() const
{
  return hex8;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // kintype_
  AddtoPack(data,kintype_);
  // eastype_
  AddtoPack(data,eastype_);
  // neas_
  AddtoPack(data,neas_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  // glprestrain_
  AddtoPack(data,*glprestrain_);
#endif

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Unpack(const vector<char>& data)
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
  // kintype_
  ExtractfromPack(position,data,kintype_);
  // eastype_
  ExtractfromPack(position,data,eastype_);
  // neas_
  ExtractfromPack(position,data,neas_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);
#if defined(PRESTRESS) || defined(POSTSTRESS)
  // glprestrain_
  ExtractfromPack(position,data,*glprestrain_);
#endif
  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size;
  ExtractfromPack(position,data,size);
  invJ_.resize(size);
  for (int i=0; i<size; ++i)
  {
    invJ_[i].Shape(0,0);
    ExtractfromPack(position,data,invJ_[i]);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::~So_hex8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Print(ostream& os) const
{
  os << "So_hex8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_expol(Epetra_SerialDenseMatrix& stresses,
                                        Epetra_SerialDenseMatrix& nodalstresses)
{
  static Epetra_SerialDenseMatrix expol(NUMNOD_SOH8,NUMGPT_SOH8);
  static bool isfilled;

  if (isfilled==true)
  {
    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);
  }
  else
  {
    double sq3=sqrt(3);
    expol(0,0)=1.25+0.75*sq3;
    expol(0,1)=-0.25-0.25*sq3;
    expol(0,2)=-0.25+0.25*sq3;
    expol(0,3)=-0.25-0.25*sq3;
    expol(0,4)=-0.25-0.25*sq3;
    expol(0,5)=-0.25+0.25*sq3;
    expol(0,6)=1.25-0.75*sq3;
    expol(0,7)=-0.25+0.25*sq3;
    expol(1,1)=1.25+0.75*sq3;
    expol(1,2)=-0.25-0.25*sq3;
    expol(1,3)=-0.25+0.25*sq3;
    expol(1,4)=-0.25+0.25*sq3;
    expol(1,5)=-0.25-0.25*sq3;
    expol(1,6)=-0.25+0.25*sq3;
    expol(1,7)=1.25-0.75*sq3;
    expol(2,2)=1.25+0.75*sq3;
    expol(2,3)=-0.25-0.25*sq3;
    expol(2,4)=1.25-0.75*sq3;
    expol(2,5)=-0.25+0.25*sq3;
    expol(2,6)=-0.25-0.25*sq3;
    expol(2,7)=-0.25+0.25*sq3;
    expol(3,3)=1.25+0.75*sq3;
    expol(3,4)=-0.25+0.25*sq3;
    expol(3,5)=1.25-0.75*sq3;
    expol(3,6)=-0.25+0.25*sq3;
    expol(3,7)=-0.25-0.25*sq3;
    expol(4,4)=1.25+0.75*sq3;
    expol(4,5)=-0.25-0.25*sq3;
    expol(4,6)=-0.25+0.25*sq3;
    expol(4,7)=-0.25-0.25*sq3;
    expol(5,5)=1.25+0.75*sq3;
    expol(5,6)=-0.25-0.25*sq3;
    expol(5,7)=-0.25+0.25*sq3;
    expol(6,6)=1.25+0.75*sq3;
    expol(6,7)=-0.25-0.25*sq3;
    expol(7,7)=1.25+0.75*sq3;

    for (int i=0;i<NUMNOD_SOH8;++i)
    {
      for (int j=0;j<i;++j)
      {
        expol(i,j)=expol(j,i);
      }
    }

    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);

    isfilled = true;
  }
}

/*----------------------------------------------------------------------*
 |  allocate and return So_hex8Register (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_hex8::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Soh8Register(Type()));
}

  /*====================================================================*/
  /* 8-node hexhedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (r,s,t) of nodes
   * of biunit cube [-1,1]x[-1,1]x[-1,1]
   *  8-node hexahedron: node 0,1,...,7
   *                      t
   *                      |
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------s
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      r                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Surfaces()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisNames(map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);

//  // element fiber direction vector
//  if (fiberdirection_.size()!=0)
//  {
//    string fibervecname = "FiberVec";
//    names[fibervecname] = 3;
//  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisData(name,data);

//  // these are the names so_hex8 recognizes, do nothing for everything else
//  if (name != "FiberVec") return;
//
//  // check sizes
//  if ((name == "FiberVec") && (data.size()!=fiberdirection_.size()))
//    dserror("FiberVec size mismatch ");
//
//  if (name == "FiberVec"){
//      data = fiberdirection_;
//  }
//  else dserror("weirdo impossible case????");

  return;
}



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register::Soh8Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register::Soh8Register(
                               const DRT::ELEMENTS::Soh8Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Register* DRT::ELEMENTS::Soh8Register::Clone() const
{
  return new DRT::ELEMENTS::Soh8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Soh8Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Soh8Register::~Soh8Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Register::Print(ostream& os) const
{
  os << "Soh8Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
