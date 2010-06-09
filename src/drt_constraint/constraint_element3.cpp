/*!----------------------------------------------------------------------
\file constraint_element3.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "constraint_element3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"

using namespace DRT::UTILS;


DRT::ELEMENTS::ConstraintElement3Type DRT::ELEMENTS::ConstraintElement3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::ConstraintElement3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::ConstraintElement3* object = new DRT::ELEMENTS::ConstraintElement3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ConstraintElement3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="CONSTRELE3" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::ConstraintElement3(id,owner));
    return ele;
  }
  return Teuchos::null;
}


void DRT::ELEMENTS::ConstraintElement3Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
}

void DRT::ELEMENTS::ConstraintElement3Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
}


DRT::ELEMENTS::ConstraintElement3RegisterType DRT::ELEMENTS::ConstraintElement3RegisterType::instance_;


DRT::ParObject* DRT::ELEMENTS::ConstraintElement3RegisterType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::ConstraintElement3Register* object =
    new DRT::ELEMENTS::ConstraintElement3Register(DRT::Element::element_constraintelement3);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(int id, int owner) :
DRT::Element(id, owner),
data_()
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::ConstraintElement3(const DRT::ELEMENTS::ConstraintElement3& old) :
DRT::Element(old),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ConstraintElement3::Clone() const
{
  DRT::ELEMENTS::ConstraintElement3* newelement = new DRT::ELEMENTS::ConstraintElement3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
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
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
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

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3::~ConstraintElement3()
{
  return;
}


/*----------------------------------------------------------------------*
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
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register::ConstraintElement3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register::ConstraintElement3Register(
                               const DRT::ELEMENTS::ConstraintElement3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register* DRT::ELEMENTS::ConstraintElement3Register::Clone() const
{
  return new DRT::ELEMENTS::ConstraintElement3Register(*this);
}

/*----------------------------------------------------------------------*
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
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3Register::Unpack(const vector<char>& data)
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
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ConstraintElement3Register::~ConstraintElement3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ConstraintElement3Register::Print(ostream& os) const
{
  os << "ConstraintElement3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
