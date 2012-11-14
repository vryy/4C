
/*!----------------------------------------------------------------------
\file red_acinus.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
            </pre>

*----------------------------------------------------------------------*/

#include "red_airway.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;

DRT::ELEMENTS::RedAcinusType DRT::ELEMENTS::RedAcinusType::instance_;


DRT::ParObject* DRT::ELEMENTS::RedAcinusType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::RedAcinus* object = new DRT::ELEMENTS::RedAcinus(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="RED_ACINUS" )
  {
    Teuchos::RCP<DRT::Element> ele =  Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele =  Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id,owner));
  return ele;
}


void DRT::ELEMENTS::RedAcinusType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["RED_ACINUS"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedString("TYPE")
    .AddNamedDouble("AcinusVolume")
    .AddNamedDouble("AlveolarDuctVolume")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(int id, int owner) :
DRT::Element(id,owner),
data_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(const DRT::ELEMENTS::RedAcinus& old) :
DRT::Element(old),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAcinus and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedAcinus::Clone() const
{
  DRT::ELEMENTS::RedAcinus* newelement = new DRT::ELEMENTS::RedAcinus(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::RedAcinus::Shape() const
{
  switch (NumNode())
  {
  case  2: return line2;
  case  3: return line3;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  Element::Pack(data);

  AddtoPack(data,elemType_);
  AddtoPack(data,resistance_);

  map<std::string,double>::const_iterator it;

  AddtoPack(data,(int)(elemParams_.size()));
  for (it = elemParams_.begin(); it!= elemParams_.end(); it++)
  {
    AddtoPack(data,it->first);
    AddtoPack(data,it->second);
  }

  AddtoPack(data,generation_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);

  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  ExtractfromPack(position,data,elemType_);
  ExtractfromPack(position,data,resistance_);
  map<std::string,double> it;
  int n = 0;

  ExtractfromPack(position,data,n);

  for (int i = 0; i<n; i++)
  {
    std::string name;
    double val;
    ExtractfromPack(position,data,name);
    ExtractfromPack(position,data,val);
    elemParams_[name] = val;
  }

  // extract generation
  ExtractfromPack(position,data,generation_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                           ismail 01/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::~RedAcinus()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Print(ostream& os) const
{
  os << "RedAcinus ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::VisNames(std::map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);

#if 0
  // see whether we have additional data for visualization in our container
  ostringstream temp;
  temp << 1;

  // in flow of volumetric flow profile
  string name = "flow_in";
  names.insert(std::pair<string,int>(name,1));

  // out flow of volumetric flow profile
  name = "flow_out";
  names.insert(std::pair<string,int>(name,1));
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if(DRT::Element::VisData(name,data))
    return true;

  return false;
}



/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::getParams(std::string name, double & var)
{

  map<std::string,double>::iterator it;
  it = elemParams_.find(name);
  if (it == elemParams_.end())
  {
    dserror ("[%s] is not found with in the element variables",name.c_str());
    exit(1);
  }
  var = elemParams_[name];

}

/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 03/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::getParams(std::string name, int & var)
{

  if (name == "Generation")
  {
    var = generation_; 
  }
  else
  {
    dserror ("[%s] is not found with in the element INT variables",name.c_str());
    exit(1);
  }

}

