/*!----------------------------------------------------------------------
\file truss2.cpp
\brief two dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "truss2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Truss2Type DRT::ELEMENTS::Truss2Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Truss2Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Truss2* object = new DRT::ELEMENTS::Truss2(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss2Type::Create( const std::string eletype,
                                                              const std::string eledistype,
                                                              const int    id,
                                                              const int    owner )
{
  if ( eletype=="TRUSS2" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss2(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss2Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss2(id,owner));
  return ele;
}


void DRT::ELEMENTS::Truss2Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 2;
  dimns = 3;
}

void DRT::ELEMENTS::Truss2Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
}

void DRT::ELEMENTS::Truss2Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["TRUSS2"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedString("KINEM")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2::Truss2(int id, int owner) :
DRT::Element(id,owner),
data_(),
isinit_(false),
material_(0),
lrefe_(0),
crosssec_(0),
kintype_(tr2_totlag),

//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2::Truss2(const DRT::ELEMENTS::Truss2& old) :
 DRT::Element(old),
 data_(old.data_),
 isinit_(old.isinit_),
 X_(old.X_),
 material_(old.material_),
 lrefe_(old.lrefe_),
 crosssec_(old.crosssec_),
 kintype_(old. kintype_),
 gaussrule_(old.gaussrule_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Truss2 and return pointer to it (public) |
 |                                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss2::Clone() const
{
  DRT::ELEMENTS::Truss2* newelement = new DRT::ELEMENTS::Truss2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2::~Truss2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Print(ostream& os) const
{
  os << "Truss2 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Truss2::Shape() const
{
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  //whether element has already been initialized
  AddtoPack(data,isinit_);
  //nodal reference coordinates
  AddtoPack(data,X_);
  //material type
  AddtoPack(data,material_);
  //reference length
  AddtoPack(data,lrefe_);
  //cross section
  AddtoPack(data,crosssec_);
  // gaussrule_
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  //kinematic type
  AddtoPack(data,kintype_);
  AddtoPack(data,data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Unpack(const std::vector<char>& data)
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
  //whether element has already been initialized
  isinit_ = ExtractInt(position,data);
  //nodal reference coordinates
  ExtractfromPack(position,data,X_);
  //material type
  ExtractfromPack(position,data,material_);
  //reference length
  ExtractfromPack(position,data,lrefe_);
  //cross section
  ExtractfromPack(position,data,crosssec_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  // kinematic type
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );
  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Truss2::Lines()
{
  std::vector<RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

void DRT::ELEMENTS::Truss2::SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe)
{
  /*this method initialized geometric variables of the element; such an initialization can only be done one time when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this
   * method will initialize the geometric variables iff the class variable isinit_ == false and afterwards set this variable to
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/
  if(!isinit_)
  {
    isinit_ = true;

    //setting reference coordinates
    X_ = xrefe;

    //length in reference configuration
    lrefe_ = std::pow(pow(X_(2)-X_(0),2)+pow(X_(3)-X_(1),2),0.5);
  }

  return;
}




int DRT::ELEMENTS::Truss2Type::Initialize(DRT::Discretization& dis)
{
  //reference node positions
  LINALG::Matrix<4,1> xrefe;

  //setting beam reference director correctly
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {
    //in case that current element is not a truss2 element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    //if we get so far current element is a truss2 element and  we get a pointer at it
    DRT::ELEMENTS::Truss2* currele = dynamic_cast<DRT::ELEMENTS::Truss2*>(dis.lColElement(i));
    if (!currele) dserror("cast to Truss2* failed");

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k=0; k<2; k++) //element has two nodes
        for(int l= 0; l < 2; l++)
          xrefe(k*2 + l) = currele->Nodes()[k]->X()[l];
    }

    currele->SetUpReferenceGeometry(xrefe);


  } //for (int i=0; i<dis_.NumMyColElements(); ++i)


  return 0;
}


