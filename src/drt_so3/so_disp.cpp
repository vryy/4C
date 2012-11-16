/*!----------------------------------------------------------------------
\file so_disp.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include "so_disp.H"
#include "so_line.H"
#include "so_surface.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::SoDispType DRT::ELEMENTS::SoDispType::instance_;


DRT::ParObject* DRT::ELEMENTS::SoDispType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::SoDisp* object =
    new DRT::ELEMENTS::SoDisp(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoDispType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLID3" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoDisp(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoDispType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoDisp(id,owner));
  return ele;
}


void DRT::ELEMENTS::SoDispType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::SoDispType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::SoDispType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLID3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    .AddNamedInt("GP_PYRAMID")
    .AddNamedInt("GP_TET")
    .AddNamedString("GP_ALT")
    .AddNamedString("KINEM")
    ;

}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(int id, int owner) :
DRT::Element(id,owner),
kintype_(sodisp_totlag),
stresstype_(sodisp_stress_none),
gaussrule_(DRT::UTILS::intrule3D_undefined),
numnod_disp_(-1),
numdof_disp_(-1),
numgpt_disp_(-1)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::SoDisp(const DRT::ELEMENTS::SoDisp& old) :
DRT::Element(old),
kintype_(old.kintype_),
stresstype_(old.stresstype_),
gaussrule_(old.gaussrule_),
numnod_disp_(old.numnod_disp_),
numdof_disp_(old.numdof_disp_),
numgpt_disp_(old.numgpt_disp_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoDisp::Clone() const
{
  DRT::ELEMENTS::SoDisp* newelement = new DRT::ELEMENTS::SoDisp(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::SoDisp::Shape() const
{
    switch (NumNode())
    {
    case  4: return tet4;
    case  8: return hex8;
    case 10: return tet10;
    case 20: return hex20;
    case 27: return hex27;
    case  6: return wedge6;
    case 15: return wedge15;
    case  5: return pyramid5;
    default:
      dserror("unexpected number of nodes %d", NumNode());
    }
    return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  AddtoPack(data,stresstype_);
  AddtoPack(data,kintype_);
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  AddtoPack(data,numnod_disp_);
  AddtoPack(data,numdof_disp_);
  AddtoPack(data,numgpt_disp_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Unpack(const std::vector<char>& data)
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

  stresstype_ = static_cast<StressType>( ExtractInt(position,data) );
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );

  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule3D(gausrule_integer); //explicit conversion from integer to enum

  ExtractfromPack(position,data,numnod_disp_);
  ExtractfromPack(position,data,numdof_disp_);
  ExtractfromPack(position,data,numgpt_disp_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::~SoDisp()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::Print(ostream& os) const
{
  os << "SoDisp ";
  Element::Print(os);
  cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::SoDisp::Lines()
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
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::SoDisp::Surfaces()
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
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::SoDisp::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}


