/*!----------------------------------------------------------------------
\file so_hex8_poro.cpp
\brief

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "so_hex8_poro.H"
//#include "../drt_lib/drt_discret.H"
//#include "../drt_lib/drt_utils.H"
//#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"


DRT::ELEMENTS::So_hex8PoroType DRT::ELEMENTS::So_hex8PoroType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_hex8PoroType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_hex8_poro* object = new DRT::ELEMENTS::So_hex8_poro(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDH8PORO" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex8_poro(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_hex8_poro(id,owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}


void DRT::ELEMENTS::So_hex8PoroType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}


void DRT::ELEMENTS::So_hex8PoroType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8PORO"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddNamedString("EAS")
    .AddOptionalNamedDoubleVector("RAD",3)
    .AddOptionalNamedDoubleVector("AXI",3)
    .AddOptionalNamedDoubleVector("CIR",3)
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8_poro::So_hex8_poro(int id, int owner) :
DRT::Element(id,owner),
DRT::ELEMENTS::So_hex8(id,owner),
DRT::ELEMENTS::So3_Poro<DRT::Element::hex8>(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8_poro::So_hex8_poro(const DRT::ELEMENTS::So_hex8_poro& old) :
    DRT::Element(old),
    DRT::ELEMENTS::So_hex8(old),
    DRT::ELEMENTS::So3_Poro<DRT::Element::hex8>(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_hex8_poro::Clone() const
{
  DRT::ELEMENTS::So_hex8_poro* newelement = new DRT::ELEMENTS::So_hex8_poro(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8_poro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_hex8 Element
  DRT::ELEMENTS::So_hex8::Pack(data);
  // add base class So3_poro Element
  DRT::ELEMENTS::So3_Poro<DRT::Element::hex8>::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8_poro::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class So_hex8 Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::So_hex8::Unpack(basedata);
  // extract base class So3_poro Element
  basedata.clear();
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::So3_Poro<DRT::Element::hex8>::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            vuong 03/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8_poro::~So_hex8_poro()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8_poro::Print(ostream& os) const
{
  os << "So_hex8_poro ";
  Element::Print(os);
  cout << endl;
  cout << So_hex8::data_;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8_poro::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  return So_hex8::ReadElement(eletype, eledistype, linedef);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8_poro::NumDofPerNode(const unsigned nds,
                                              const DRT::Node& node) const
{
  if (nds==1)
    return So3_Poro<DRT::Element::hex8>::NumDofPerNode(nds, node);

  return So_hex8::NumDofPerNode(node);
};

#endif
