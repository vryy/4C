/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_poro.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "fluid_ele_poro.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_inpar/inpar_fluid.H"

DRT::ELEMENTS::FluidPoroEleType DRT::ELEMENTS::FluidPoroEleType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::FluidPoroEleType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::FluidPoro* object = new DRT::ELEMENTS::FluidPoro(-1,-1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroEleType::Create(const std::string  eletype,
                                                             const std::string  eledistype,
                                                             const int     id,
                                                             const int     owner)
{
  if (eletype=="FLUIDPORO")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::FluidPoro(id,owner));
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroEleType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::FluidPoro(id,owner));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoroEleType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_fluid;
  FluidType::SetupElementDefinition(definitions_fluid);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_fluid =
      definitions_fluid["FLUID"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["FLUIDPORO"];

  //3D
  defs["HEX8"]     = defs_fluid["HEX8"];
  defs["HEX20"]    = defs_fluid["HEX20"];
  defs["HEX27"]    = defs_fluid["HEX27"];
  defs["TET4"]     = defs_fluid["TET4"];
  defs["TET10"]    = defs_fluid["TET10"];
  defs["WEDGE6"]   = defs_fluid["WEDGE6"];
  defs["WEDGE15"]  = defs_fluid["WEDGE15"];
  defs["PYRAMID5"] = defs_fluid["PYRAMID5"];
  defs["NURBS8"]   = defs_fluid["NURBS8"];
  defs["NURBS27"]  = defs_fluid["NURBS27"];

  //2D
  defs["QUAD4"]    = defs_fluid["QUAD4"];
  defs["QUAD8"]    = defs_fluid["QUAD8"];
  defs["QUAD9"]    = defs_fluid["QUAD9"];
  defs["TRI3"]     = defs_fluid["TRI3"];
  defs["TRI6"]     = defs_fluid["TRI6"];
  defs["NURBS4"]   = defs_fluid["NURBS4"];
  defs["NURBS9"]   = defs_fluid["NURBS9"];
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoro::FluidPoro(int id, int owner) :
Fluid(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoro::FluidPoro(const DRT::ELEMENTS::FluidPoro& old) :
Fluid(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                                     |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidPoro::Clone() const
{
  DRT::ELEMENTS::FluidPoro* newelement = new DRT::ELEMENTS::FluidPoro(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  Fluid::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::ExtractfromPack(position,data,basedata);
  Fluid::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                           |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::FluidPoro::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<FluidBoundary,FluidPoro>(DRT::UTILS::buildLines,this);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
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
 |  get vector of surfaces (public)                                     |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::FluidPoro::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<FluidBoundary,FluidPoro>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1) // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element> > surfaces(1);
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
 |  get vector of volumes (length 1) (public)                            |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::FluidPoro::Volumes()
{
  if (NumVolume()==1) // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element> > volumes(1);
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
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoro::Print(std::ostream& os) const
{
  os << "FluidPoro ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidPoro::NumDofPerNode(const unsigned nds, const DRT::Node& node, const std::string disname) const
{
  if (nds==0)
  {
    return Fluid::NumDofPerNode(node);
  }
  else if (nds==1)
  {
    //get number of dimensions
    const int nsd = DRT::UTILS::getDimension(Shape());
    //get fluid dynamic parameters
    const Teuchos::ParameterList& params= DRT::Problem::Instance()->FluidDynamicParams();
    //get physical type from parameter list
    INPAR::FLUID::PhysicalType physicaltype = DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(params,"PHYSICAL_TYPE");

    switch(physicaltype)
    {
    case INPAR::FLUID::poro:
    case INPAR::FLUID::poro_p2:
      //second dofset is porous structure -> nsd dofs per node
      return nsd;
      break;
    case INPAR::FLUID::poro_p1:
      //second dofset is porous structure (+porosity dof) -> nsd+1 dofs per node
      return nsd+1;
      break;
    default:
      dserror("invalid fluid physical type for porous media");
      return -1;
      break;
    }
  }
  else if(nds==2)
  {
    // what's the current problem type?
    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    switch (probtype)
    {
      case prb_poroscatra:
      {
        //scalar transport -> 1 dof per node (assuming only one scalar being transported)
        return 1;
        break;
      }
      default:
      {
        dserror("invalid number of dofsets (3) for this problem type");
        return -1;
        break;
      }
    }
  }
  else
  {
    dserror("invalid number of dof sets");
    return -1;
  }
  return -1;
}
