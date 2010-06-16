/*!----------------------------------------------------------------------
\file fluid3.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;

DRT::ELEMENTS::Fluid3Type DRT::ELEMENTS::Fluid3Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Fluid3Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Fluid3* object = new DRT::ELEMENTS::Fluid3(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Fluid3Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="FLUID3" )
  {
    if ( eledistype!="NURBS8" and eledistype!="NURBS27" )
    {
      return rcp(new DRT::ELEMENTS::Fluid3(id,owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Fluid3Type::Create( const int id, const int owner )
{
  return rcp(new DRT::ELEMENTS::Fluid3(id,owner));
}


void DRT::ELEMENTS::Fluid3Type::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf-1;
  np = 1;
}


void DRT::ELEMENTS::Fluid3Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluid3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Fluid3Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
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
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Fluid3BoundaryType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new Fluid3Boundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
// map to convert strings to actions (stabilization)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::Fluid3::StabilisationAction> DRT::ELEMENTS::Fluid3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(int id, int owner) :
DRT::Element(id,owner),
is_ale_(false),
data_()
{
    distype_= dis_none;

    Cs_delta_sq_=0;

    saccn_ .Shape(0,0);
    svelnp_.Shape(0,0);
    sveln_ .Shape(0,0);

    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(const DRT::ELEMENTS::Fluid3& old) :
DRT::Element(old             ),
distype_    (old.distype_    ),
is_ale_     (old.is_ale_     ),
data_       (old.data_       ),
Cs_delta_sq_(old.Cs_delta_sq_),
saccn_      (old.saccn_      ),
svelnp_     (old.svelnp_     ),
sveln_      (old.sveln_      )
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid3 and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid3::Clone() const
{
  DRT::ELEMENTS::Fluid3* newelement = new DRT::ELEMENTS::Fluid3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // is_ale_
  AddtoPack(data,is_ale_);
  // Cs_delta_sq_, the Smagorinsky constant for the dynamic Smagorinsky model
  AddtoPack(data,Cs_delta_sq_);
  // Discretisation type
  AddtoPack(data,distype_);

  // history variables
  AddtoPack(data,saccn_.M());
  AddtoPack(data,saccn_.N());

  int size = saccn_.M()*saccn_.N()*sizeof(double);

  AddtoPack(data,saccn_ .A(),size);
  AddtoPack(data,svelnp_.A(),size);
  AddtoPack(data,sveln_ .A(),size);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Unpack(const vector<char>& data)
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
  // is_ale_
  ExtractfromPack(position,data,is_ale_);
  // extract Cs_delta_sq_, the Smagorinsky constant for the dynamic
  // Smagorinsky model
  ExtractfromPack(position,data,Cs_delta_sq_);
  // distype
  ExtractfromPack(position,data,distype_);

  // history variables (subscale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position,data,firstdim);
    ExtractfromPack(position,data,secondim);

    saccn_ .Shape(firstdim,secondim);
    svelnp_.Shape(firstdim,secondim);
    sveln_ .Shape(firstdim,secondim);

    int size = firstdim*secondim*sizeof(double);

    ExtractfromPack(position,data,&(saccn_ .A()[0]),size);
    ExtractfromPack(position,data,&(svelnp_.A()[0]),size);
    ExtractfromPack(position,data,&(sveln_ .A()[0]),size);
  }

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::~Fluid3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Print(ostream& os) const
{
  os << "Fluid3 ";
  Element::Print(os);
  //cout << endl;
  cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                 ae  02/010|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine()>1) // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<Fluid3Boundary,Fluid3>(DRT::UTILS::buildLines,this);
  }
  else if (NumLine()==1) // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    vector<RCP<Element> > surfaces(1);
    surfaces[0]= rcp(this, false);
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)   // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<Fluid3Boundary,Fluid3>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1) // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
  {
    // 2D (we return the element itself)
    vector<RCP<Element> > surfaces(1);
    surfaces[0]= rcp(this, false);
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Volumes()
{
  if (NumVolume()==1) // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    vector<RCP<Element> > volumes(1);
    volumes[0]= rcp(this, false);
    return volumes;
  }
  else //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  activate time dependend subscales (public)           gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::ActivateTDS(int nquad,int nsd)
   {
     if(saccn_.M() != nsd
        ||
        saccn_.N() != nquad)
     {
       saccn_ .Shape(nsd,nquad);
       memset(saccn_.A() ,0,nsd*nquad*sizeof(double));

       sveln_ .Shape(nsd,nquad);
       memset(sveln_.A() ,0,nsd*nquad*sizeof(double));

       svelnp_.Shape(nsd,nquad);
       memset(svelnp_.A(),0,nsd*nquad*sizeof(double));
     }
     return;
   }


/*----------------------------------------------------------------------*
 |  activate time dependend subscales (public)           gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::UpdateSvelnpInOneDirection(
    const double  fac1   ,
    const double  fac2   ,
    const double  fac3   ,
    const double  resM   ,
    const double  alphaF ,
    const int     dim    ,
    const int     iquad  ,
    double&       svelaf
    )
{

    /*
        ~n+1           ~n           ~ n            n+1
        u    =  fac1 * u  + fac2 * acc  -fac3 * res
         (i)

    */

  svelnp_(dim,iquad)=
    fac1*sveln_(dim,iquad)
    +
    fac2*saccn_(dim,iquad)
    -
    fac3*resM;

  /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

  */
  svelaf=
    alphaF      *svelnp_(dim,iquad)
    +
    (1.0-alphaF)*sveln_ (dim,iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  activate time dependend subscales (public)           gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::UpdateSvelnpInOneDirection(
    const double  fac1   ,
    const double  fac2   ,
    const double  fac3   ,
    const double  resM   ,
    const double  alphaF ,
    const int     dim    ,
    const int     iquad  ,
    double&       svelnp,
    double&       svelaf
    )
    {
      UpdateSvelnpInOneDirection(fac1   ,
                                 fac2   ,
                                 fac3   ,
                                 resM   ,
                                 alphaF ,
                                 dim    ,
                                 iquad  ,
                                 svelaf );

      svelnp=svelnp_(dim,iquad);

      return;
    }


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
