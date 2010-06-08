/*!----------------------------------------------------------------------
\file torsion3.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TORSION3
#ifdef CCADISCRET

#include "torsion3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_fixedsizematrix.H"

DRT::ELEMENTS::Torsion3Type DRT::ELEMENTS::Torsion3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Torsion3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Torsion3* object = new DRT::ELEMENTS::Torsion3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Torsion3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="TORSION3" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::Torsion3(id,owner));
    return ele;
  }
  return Teuchos::null;
}


DRT::ELEMENTS::Torsion3RegisterType DRT::ELEMENTS::Torsion3RegisterType::instance_;


DRT::ParObject* DRT::ELEMENTS::Torsion3RegisterType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Torsion3Register* object =
    new DRT::ELEMENTS::Torsion3Register(DRT::Element::element_torsion3);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(int id, int owner) :
DRT::Element(id,element_torsion3,owner),
data_(),
springconstant_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(const DRT::ELEMENTS::Torsion3& old) :
 DRT::Element(old),
 data_(old.data_),
 springconstant_(old.springconstant_)
{
  return;
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Torsion3 and return pointer to it (public)|
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Torsion3::Clone() const
{
  DRT::ELEMENTS::Torsion3* newelement = new DRT::ELEMENTS::Torsion3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::~Torsion3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Print(ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Torsion3Register (public)             cyron 02/10|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Torsion3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Torsion3Register(Type()));
}


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Torsion3::Shape() const
{
  return line3;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  AddtoPack(data,springconstant_);
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  ExtractfromPack(position,data,springconstant_);
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Torsion3::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}


//------------- class Torsion3Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::Torsion3Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::Torsion3Register(
                               const DRT::ELEMENTS::Torsion3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register* DRT::ELEMENTS::Torsion3Register::Clone() const
{
  return new DRT::ELEMENTS::Torsion3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Pack(vector<char>& data) const
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


/*-----------------------------------------------------------------------*
 |  Unpack data (public)                                      cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::~Torsion3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Print(ostream& os) const
{
  os << "Torsion3Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Torsion3Register::Initialize(DRT::Discretization& dis)
{
  //no special initilization required for this element

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION3

