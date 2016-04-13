/*!----------------------------------------------------------------------
\file shell8.cpp

\maintainer Michael Gee

*----------------------------------------------------------------------*/
#include "shell8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"


DRT::ELEMENTS::Shell8Type DRT::ELEMENTS::Shell8Type::instance_;


DRT::ELEMENTS::Shell8Type& DRT::ELEMENTS::Shell8Type::Instance()
{
  return instance_;
}


DRT::ParObject* DRT::ELEMENTS::Shell8Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Shell8* object = new DRT::ELEMENTS::Shell8(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell8Type::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SHELL8" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Shell8(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell8Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Shell8(id,owner));
  return ele;
}


void DRT::ELEMENTS::Shell8Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 6;
  dimns = 6;
  nv = 6;
}

void DRT::ELEMENTS::Shell8Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Shell8Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SHELL8"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedDouble("SDC")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedDouble("SDC")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedDouble("SDC")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedDouble("SDC")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedDouble("THICK")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_TRI")
    .AddNamedString("FORCES")
    .AddNamedString("EAS")
    .AddString("EAS2")
    .AddString("EAS3")
    .AddString("EAS4")
    .AddString("EAS5")
    .AddNamedString("ANS")
    .AddNamedDouble("SDC")
    ;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell8LineType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new Shell8Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::Shell8(int id, int owner) :
DRT::Element(id,owner),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0),
data_(),
old_step_length_(-1.0),
interface_ptr_(Teuchos::null)
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::Shell8(const DRT::ELEMENTS::Shell8& old) :
DRT::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_),
data_(old.data_),
old_step_length_(old.old_step_length_),
interface_ptr_(old.interface_ptr_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Shell8::Clone() const
{
  DRT::ELEMENTS::Shell8* newelement = new DRT::ELEMENTS::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Shell8::Shape() const
{
  switch (NumNode())
  {
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // forcetype_
  AddtoPack(data,forcetype_);
  // thickness_
  AddtoPack(data,thickness_);
  // ngp_
  AddtoPack(data,ngp_,3*sizeof(int));
  // ngptri_
  AddtoPack(data,ngptri_);
  // nhyb_
  AddtoPack(data,nhyb_);
  // eas_
  AddtoPack(data,eas_,5*sizeof(int));
  // ans_
  AddtoPack(data,ans_);
  // sdc_
  AddtoPack(data,sdc_);
  // material_
  AddtoPack(data,material_);
  // old_step_length_
  AddtoPack(data,old_step_length_);
  // data_
  AddtoPack(data,data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // forcetype_
  forcetype_ = static_cast<ForceType>( ExtractInt(position,data) );
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeof(int));
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // old_step_length_
  ExtractfromPack(position,data,old_step_length_);
  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Print(std::ostream& os) const
{
  os << "Shell8 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mwgee 01/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Shell8::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Shell8Line,Shell8>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mwgee 01/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Shell8::Surfaces()
{
  std::vector<Teuchos::RCP<Element> > surfaces(1);
  surfaces[0]= Teuchos::rcp(this, false);
  return surfaces;
}
