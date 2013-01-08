/*!----------------------------------------------------------------------
\file fluid_ele.cpp
\brief

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>

*----------------------------------------------------------------------*/

#include "fluid_ele.H"
#include "fluid_ele_tds.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"


DRT::ELEMENTS::FluidType DRT::ELEMENTS::FluidType::instance_;


DRT::ParObject* DRT::ELEMENTS::FluidType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Fluid* object = new DRT::ELEMENTS::Fluid(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidType::Create(const string  eletype,
                                                             const string  eledistype,
                                                             const int     id,
                                                             const int     owner)
{
  if ( eletype=="FLUID3" )
  {
      return Teuchos::rcp(new DRT::ELEMENTS::Fluid(id,owner));
  }
  else if ( eletype=="FLUID2" )
  {
      return Teuchos::rcp(new DRT::ELEMENTS::Fluid(id,owner));
  }
  else if (eletype=="FLUID")
  {
      return Teuchos::rcp(new DRT::ELEMENTS::Fluid(id,owner));
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::Fluid(id,owner));
}


void DRT::ELEMENTS::FluidType::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf-1;
  np = 1;
}


void DRT::ELEMENTS::FluidType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluidDNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::FluidType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["FLUID3"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX20"]
    .AddIntVector("HEX20",20)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["HEX27"]
    .AddIntVector("HEX27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["WEDGE15"]
    .AddIntVector("WEDGE15",15)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["PYRAMID5"]
    .AddIntVector("PYRAMID5",5)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS8"]
    .AddIntVector("NURBS8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS27"]
    .AddIntVector("NURBS27",27)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  // 2D elements
  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  std::map<std::string,DRT::INPUT::LineDefinition>& defs2D = definitions["FLUID2"];

   defs2D["QUAD4"]
     .AddIntVector("QUAD4",4)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["QUAD8"]
     .AddIntVector("QUAD8",8)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["QUAD9"]
     .AddIntVector("QUAD9",9)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["TRI3"]
     .AddIntVector("TRI3",3)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["TRI6"]
     .AddIntVector("TRI6",6)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["NURBS4"]
     .AddIntVector("NURBS4",4)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["NURBS9"]
     .AddIntVector("NURBS9",9)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   defs2D["THQ9"]
     .AddIntVector("THQ9",9)
     .AddNamedInt("MAT")
     .AddNamedString("NA")
     ;

   std::map<std::string,DRT::INPUT::LineDefinition>& defsgeneral = definitions["FLUID"];

     defsgeneral["HEX8"]
       .AddIntVector("HEX8",8)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["HEX20"]
       .AddIntVector("HEX20",20)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["HEX27"]
       .AddIntVector("HEX27",27)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["TET4"]
       .AddIntVector("TET4",4)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["TET10"]
       .AddIntVector("TET10",10)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["WEDGE6"]
       .AddIntVector("WEDGE6",6)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["WEDGE15"]
       .AddIntVector("WEDGE15",15)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["PYRAMID5"]
       .AddIntVector("PYRAMID5",5)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["NURBS8"]
       .AddIntVector("NURBS8",8)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["NURBS27"]
       .AddIntVector("NURBS27",27)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     // 2D elements
     defsgeneral["QUAD4"]
       .AddIntVector("QUAD4",4)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["QUAD8"]
       .AddIntVector("QUAD8",8)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["QUAD9"]
       .AddIntVector("QUAD9",9)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["TRI3"]
       .AddIntVector("TRI3",3)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["TRI6"]
       .AddIntVector("TRI6",6)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["NURBS4"]
       .AddIntVector("NURBS4",4)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;

     defsgeneral["NURBS9"]
       .AddIntVector("NURBS9",9)
       .AddNamedInt("MAT")
       .AddNamedString("NA")
       ;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid::Fluid(int id, int owner) :
DRT::Element(id,owner),
is_ale_(false)
{
    distype_= dis_none;
    tds_=Teuchos::null;
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid::Fluid(const DRT::ELEMENTS::Fluid& old) :
DRT::Element(old             ),
distype_    (old.distype_    ),
is_ale_     (old.is_ale_     )
{
  tds_ = Teuchos::null;
  if (old.tds_ != Teuchos::null)
    dserror("Clone() method for deep copying tds_ not yet implemented!");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid::Clone() const
{
  DRT::ELEMENTS::Fluid* newelement = new DRT::ELEMENTS::Fluid(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // is_ale_
  AddtoPack(data,is_ale_);
  // Discretisation type
  AddtoPack(data,distype_);

  // time-dependent subgrid scales
  bool is_tds(false);
  if (tds_!= Teuchos::null)
  {
    is_tds = true;
    AddtoPack(data,is_tds);
    tds_->Pack(data);
  }
  else
  {
    AddtoPack(data,is_tds);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid::Unpack(const std::vector<char>& data)
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
  // is_ale_
  is_ale_ = ExtractInt(position,data);
  // distype
  distype_ = static_cast<DiscretizationType>( ExtractInt(position,data) );

  // time-dependent subgrid scales
  bool is_tds = ExtractInt(position,data);
  if (is_tds)
  {
    tds_=Teuchos::rcp(new FLD::TDSEleData());
    std::vector<char> pbtest;
    ExtractfromPack(position,data,pbtest);
    if (pbtest.size() == 0)
      dserror("Seems no TDS data available");
    tds_->Unpack(pbtest);
  }
  else
    tds_=Teuchos::null;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid::~Fluid()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid::Print(ostream& os) const
{
  os << "Fluid ";
  Element::Print(os);
  //cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                 ae  02/010|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Fluid::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<FluidBoundary,Fluid>(DRT::UTILS::buildLines,this);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          ehrl  02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Fluid::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<FluidBoundary,Fluid>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1) // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<RCP<Element> > surfaces(1);
    surfaces[0]= Teuchos::rcp(this, false);
    return surfaces;
  }
  else  // 1D elements
  {
    dserror("Surfaces() does not exist for 1D-element ");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                 ehrl 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Fluid::Volumes()
{
  if (NumVolume()==1) // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<RCP<Element> > volumes(1);
    volumes[0]= Teuchos::rcp(this, false);
    return volumes;
  }
  else //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get internal faces element (public)                     schott 03/12|
 *----------------------------------------------------------------------*/
RCP<DRT::Element> DRT::ELEMENTS::Fluid::CreateInternalFaces( DRT::Element* parent_slave,   //!< parent slave fluid3 element
                                                              int nnode,                    //!< number of surface nodes
                                                              const int* nodeids,           //!< node ids of surface element
                                                              DRT::Node** nodes,            //!< nodes of surface element
                                                              const int lsurface_master     //!< local surface number w.r.t master parent element
                                                             )
{
  // dynamic cast for slave parent element
  DRT::ELEMENTS::Fluid * slave_pele = dynamic_cast<DRT::ELEMENTS::Fluid *>( parent_slave );


  // insert both parent elements
  return DRT::UTILS::ElementIntFaceFactory<FluidIntFace,Fluid>(-1,             // internal face element id
                                                                -1,             // owner of internal face element
                                                                nnode,
                                                                nodeids,
                                                                nodes,
                                                                this,           // master parent element
                                                                slave_pele,     // slave parent element
                                                                lsurface_master // local surface number w.r.t master parent element
                                                                );

}


int DRT::ELEMENTS::Fluid::NumDofPerNode(const unsigned nds, const DRT::Node& node) const
{
  if (nds==1)
  {
    // what's the current problem type?
    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    switch (probtype)
    {
      case prb_poroelast:
      case prb_poroscatra:
      {
        return DRT::Problem::Instance()->NDim();
        break;
      }
      default: // scalar transport
      {
        return 1;
        break;
      }
    }
  }
  else if(nds==0)
    return NumDofPerNode(node);
  else
  {
    dserror("invalid number of dof sets");
    return -1;
  }
}


/*----------------------------------------------------------------------*
 |  activate time dependent subgrid scales (public)      gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid::ActivateTDS(int nquad,int nsd, double** saccn, double** sveln, double** svelnp)
{

  if (tds_ == Teuchos::null)
    tds_ = Teuchos::rcp(new FLD::TDSEleData());

  tds_->ActivateTDS(nquad, nsd, saccn, sveln, svelnp);

}
