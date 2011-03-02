/*!----------------------------------------------------------------------
\file xfluid3.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xdiff3.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;

DRT::ELEMENTS::XDiff3Type DRT::ELEMENTS::XDiff3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::XDiff3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::XDiff3* object = new DRT::ELEMENTS::XDiff3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::XDiff3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="XDIFF3" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::XDiff3(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::XDiff3Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::XDiff3(id,owner));
  return ele;
}


void DRT::ELEMENTS::XDiff3Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 1;
  dimns = 1;
  nv = 1;
}

void DRT::ELEMENTS::XDiff3Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{

}

void DRT::ELEMENTS::XDiff3Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["XDIFF3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    ;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::XDiff3SurfaceType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new XDiff3Surface( id, owner ) );
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::XDiff3LineType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new XDiff3Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
// map to convert strings to actions (stabilization)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::XDiff3::StabilisationAction> DRT::ELEMENTS::XDiff3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3::XDiff3(int id, int owner) :
DRT::Element(id,owner),
eleDofManager_(Teuchos::null),
eleDofManager_uncondensed_(Teuchos::null),
output_mode_(false)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3::XDiff3(const DRT::ELEMENTS::XDiff3& old) :
DRT::Element(old),
eleDofManager_(old.eleDofManager_),
eleDofManager_uncondensed_(old.eleDofManager_uncondensed_),
output_mode_(old.output_mode_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of XDiff3 and return pointer to it (public)|
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::XDiff3::Clone() const
{
  return new DRT::ELEMENTS::XDiff3(*this);
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::XDiff3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  case  0: return dis_none;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XDiff3::Pack(DRT::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  AddtoPack(data,output_mode_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XDiff3::Unpack(const std::vector<char>& data)
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

  ExtractfromPack(position,data,output_mode_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3::~XDiff3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XDiff3::Print(ostream& os) const
{
  os << "XDiff3 ";
  if (output_mode_)
    os << "(outputmode=true)";
  Element::Print(os);
  cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::XDiff3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<XDiff3Line,XDiff3>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::XDiff3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<XDiff3Surface,XDiff3>(DRT::UTILS::buildSurfaces,this);
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::XDiff3::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


/*----------------------------------------------------------------------*
 |  constructor
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3::MyState::MyState(
    const DRT::Discretization&      discret,
    const std::vector<int>&         lm,
    const bool                      instat
    ) :
      instationary(instat)
{
  DRT::UTILS::ExtractMyValues(*discret.GetState("velnp"),velnp,lm);
  if (instat)
  {
    DRT::UTILS::ExtractMyValues(*discret.GetState("veln") ,veln ,lm);
    DRT::UTILS::ExtractMyValues(*discret.GetState("velnm"),velnm,lm);
    DRT::UTILS::ExtractMyValues(*discret.GetState("accn") ,accn ,lm);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3::DLMInfo::DLMInfo(const int nd, const int na)
: oldKaainv_(LINALG::SerialDenseMatrix(na,na,true)),
  oldKad_(LINALG::SerialDenseMatrix(na,nd,true)),
  oldfa_(LINALG::SerialDenseVector(na,true)),
  stressdofs_(LINALG::SerialDenseVector(na,true))
{
  return;
}



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
