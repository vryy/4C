/*----------------------------------------------------------------------*/
/*!
\file reynolds_ele.cpp

\brief Reynolds elements

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15270
</pre>
 */
/*----------------------------------------------------------------------*/

#include "reynolds_ele.H"

#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

DRT::ELEMENTS::ReynoldsType DRT::ELEMENTS::ReynoldsType::instance_;

DRT::ELEMENTS::ReynoldsType& DRT::ELEMENTS::ReynoldsType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::ReynoldsType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Reynolds* object =
    new DRT::ELEMENTS::Reynolds(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ReynoldsType::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="REYNOLDS" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Reynolds(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ReynoldsType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Reynolds(id,owner));
  return ele;
}


void DRT::ELEMENTS::ReynoldsType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;

}

void DRT::ELEMENTS::ReynoldsType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluidDNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::ReynoldsType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["REYNOLDS"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    ;

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    ;

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    ;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::ReynoldsBoundaryType DRT::ELEMENTS::ReynoldsBoundaryType::instance_;

DRT::ELEMENTS::ReynoldsBoundaryType& DRT::ELEMENTS::ReynoldsBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::ReynoldsBoundaryType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new ReynoldsBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Reynolds::Reynolds(int id, int owner) :
DRT::Element(id,owner),
distype_(dis_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Reynolds::Reynolds(const DRT::ELEMENTS::Reynolds& old) :
DRT::Element(old),
distype_(old.distype_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Reynolds and return pointer to it        |
 |                                                 (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Reynolds::Clone() const
{
  DRT::ELEMENTS::Reynolds* newelement = new DRT::ELEMENTS::Reynolds(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Reynolds element                     (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Reynolds::Shape() const
{
  return distype_;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Reynolds::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  Element::Pack(data);

  // add internal data
  AddtoPack(data,distype_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Reynolds::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  // extract internal data
  distype_ = static_cast<DiscretizationType>(ExtractInt(position,data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)         wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Reynolds::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Reynolds::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)        wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Reynolds::NumVolume() const
{
  return DRT::UTILS::getNumberOfElementVolumes(distype_);
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Reynolds::~Reynolds()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Reynolds::Print(std::ostream& os) const
{
  os << "Reynolds element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << DRT::DistypeToString(distype_) << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Reynolds::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumLine() > 1) // 3D and 2D
    return DRT::UTILS::ElementBoundaryFactory<ReynoldsBoundary,Reynolds>(DRT::UTILS::buildLines,this);
  else
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > lines(1);
    lines[0]= Teuchos::rcp(this, false);
    return lines;
  }

}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                         wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Reynolds::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  if (NumSurface() > 1) // 3D
    return DRT::UTILS::ElementBoundaryFactory<ReynoldsBoundary,Reynolds>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    // 1D
    dserror("Surfaces() for 1D-Reynolds element not implemented");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)               wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Reynolds::Volumes()
{
  dserror("Volumes() for 1D-/2D-Reynolds element not implemented");
  return DRT::Element::Volumes();
}


/*----------------------------------------------------------------------*
 | read element input                                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Reynolds::ReadElement(
    const std::string& eletype,
    const std::string& distype,
    DRT::INPUT::LineDefinition* linedef
    )
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // set discretization type
  SetDisType(DRT::StringToDistype(distype));

  return true;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ReynoldsBoundary::ReynoldsBoundary(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Reynolds* parent,
                              const int lbeleid) :
DRT::FaceElement(id,owner)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent,lbeleid);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ReynoldsBoundary::ReynoldsBoundary(const DRT::ELEMENTS::ReynoldsBoundary& old) :
DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::ReynoldsBoundary::Clone() const
{
  DRT::ELEMENTS::ReynoldsBoundary* newelement = new DRT::ELEMENTS::ReynoldsBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                   (public) wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::ReynoldsBoundary::Shape() const
{
  return DRT::UTILS::getShapeOfBoundaryElement(NumNode(), ParentElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                      wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsBoundary::Pack(DRT::PackBuffer& data) const
{
  dserror("This ReynoldsBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsBoundary::Unpack(const std::vector<char>& data)
{
  dserror("This ReynoldsBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                           wirtz 10/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ReynoldsBoundary::~ReynoldsBoundary()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             wirtz 10/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ReynoldsBoundary::Print(std::ostream& os) const
{
  os << "ReynoldsBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  "<<Shape()<< std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)      wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ReynoldsBoundary::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)  wirtz 10/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ReynoldsBoundary::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::ReynoldsBoundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of ReynoldsBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                            wirtz 10/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::ReynoldsBoundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of ReynoldsBoundary not implemented");
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}
