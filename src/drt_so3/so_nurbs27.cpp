/*!----------------------------------------------------------------------
\file so_nurbs27.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "so_line.H"
#include "so_surface.H"
#include "so_nurbs27.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::NURBS::So_nurbs27Type DRT::ELEMENTS::NURBS::So_nurbs27Type::instance_;


DRT::ParObject* DRT::ELEMENTS::NURBS::So_nurbs27Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NURBS::So_nurbs27* object = new DRT::ELEMENTS::NURBS::So_nurbs27(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::So_nurbs27Type::Create( const string eletype,
                                                                         const string eledistype,
                                                                         const int id,
                                                                         const int owner )
{
  if ( eletype=="SONURBS27" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NURBS::So_nurbs27(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::So_nurbs27Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NURBS::So_nurbs27(id,owner));
  return ele;
}


void DRT::ELEMENTS::NURBS::So_nurbs27Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::NURBS::So_nurbs27Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::NURBS::So_nurbs27Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SONURBS27"];

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedIntVector("GP",3)
    ;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::So_nurbs27::So_nurbs27(int id, int owner) :
DRT::Element(id,owner),
data_()
{
  kintype_ = sonurbs27_totlag;

  invJ_.resize(NUMGPT_SONURBS27, LINALG::Matrix<NUMDIM_SONURBS27,NUMDIM_SONURBS27>(true));
  detJ_.resize(NUMGPT_SONURBS27, 0.0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::So_nurbs27::So_nurbs27(const DRT::ELEMENTS::NURBS::So_nurbs27& old) :
DRT::Element(old),
kintype_(old.kintype_),
data_   (old.data_   ),
detJ_   (old.detJ_   )
{
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::NURBS::So_nurbs27::Clone() const
{
  DRT::ELEMENTS::NURBS::So_nurbs27* newelement = new DRT::ELEMENTS::NURBS::So_nurbs27(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::So_nurbs27::Shape() const
{
  return nurbs27;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data,kintype_);
  // data_
  AddtoPack(data,data_);

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
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::Unpack(const vector<char>& data)
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
  // kintype_
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SONURBS27,NUMDIM_SONURBS27>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::So_nurbs27::~So_nurbs27()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::Print(ostream& os) const
{
  os << "So_nurbs27 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                           |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NURBS::So_nurbs27::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                                      |
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NURBS::So_nurbs27::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,So_nurbs27>(DRT::UTILS::buildSurfaces,this);

}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NURBS::So_nurbs27::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)


  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,So_nurbs27>(DRT::UTILS::buildLines,this);
}


