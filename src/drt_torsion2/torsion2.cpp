/*!----------------------------------------------------------------------
\file torsion2.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "torsion2.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Torsion2Type DRT::ELEMENTS::Torsion2Type::instance_;

DRT::ParObject* DRT::ELEMENTS::Torsion2Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Torsion2* object = new DRT::ELEMENTS::Torsion2(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Torsion2Type::Create( const std::string eletype,
                                                                const std::string eledistype,
                                                                const int    id,
                                                                const int    owner )
{
  if ( eletype=="TORSION2" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Torsion2(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Torsion2Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Torsion2(id,owner));
  return ele;
}


void DRT::ELEMENTS::Torsion2Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 2;
  dimns = 3;
}

void DRT::ELEMENTS::Torsion2Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure2DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Torsion2Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["TORSION2"];

  defs["LINE3"]
    .AddIntVector("LINE3",3)
    .AddNamedInt("MAT")
    .AddNamedString("BENDINGPOTENTIAL")
    ;

  defs["LIN3"]
    .AddIntVector("LIN3",3)
    .AddNamedInt("MAT")
    .AddNamedString("BENDINGPOTENTIAL")
    ;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion2::Torsion2(int id, int owner) :
DRT::Element(id,owner),
data_(),
isinit_(false),
theta_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion2::Torsion2(const DRT::ELEMENTS::Torsion2& old) :
 DRT::Element(old),
 data_(old.data_),
 isinit_(old.isinit_),
 theta_(old.theta_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Torsion2 and return pointer to it (public)|
 |                                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Torsion2::Clone() const
{
  DRT::ELEMENTS::Torsion2* newelement = new DRT::ELEMENTS::Torsion2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion2::~Torsion2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion2::Print(std::ostream& os) const
{

}


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Torsion2::Shape() const
{
  return line3;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion2::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  Element::Pack(data);
  AddtoPack(data,isinit_);
   AddtoPack(data,theta_);
  AddtoPack(data,bendingpotential_);
  AddtoPack(data,data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion2::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  isinit_ = ExtractInt(position,data);
  ExtractfromPack(position,data,theta_);
  bendingpotential_ = static_cast<BendingPotential>( ExtractInt(position,data) );
  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Torsion2::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

void DRT::ELEMENTS::Torsion2::SetUpReferenceGeometry(const LINALG::Matrix<6,1>& xrefe)
{
  /*this method initialized geometric variables of the element; such an initialization can only be done one time when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this
   * method will initialize the geometric variables iff the class variable isinit_ == false and afterwards set this variable to
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/
  if(!isinit_)
  {
    isinit_ = true;

    //calculation of the reference angle between the trusses
    theta_=acos( ( (xrefe(2)-xrefe(0))*(xrefe(4)-xrefe(2)) + (xrefe(3)-xrefe(1))*(xrefe(5)-xrefe(3)) )
       / sqrt( pow(xrefe(2)-xrefe(0),2)+pow(xrefe(3)-xrefe(1),2) ) / sqrt( pow(xrefe(4)-xrefe(2),2)+pow(xrefe(5)-xrefe(3),2) ) );


  if(( ( (xrefe(2)-xrefe(0))*(xrefe(4)-xrefe(2)) + (xrefe(3)-xrefe(1))*(xrefe(5)-xrefe(3)) )
         / sqrt( pow(xrefe(2)-xrefe(0),2)+pow(xrefe(3)-xrefe(1),2) )
         / sqrt( pow(xrefe(4)-xrefe(2),2)+pow(xrefe(5)-xrefe(3),2) ) ) > 1){
    if(( ( (xrefe(2)-xrefe(0))*(xrefe(4)-xrefe(2)) + (xrefe(3)-xrefe(1))*(xrefe(5)-xrefe(3)) )
           / sqrt( pow(xrefe(2)-xrefe(0),2)+pow(xrefe(3)-xrefe(1),2) )
           / sqrt( pow(xrefe(4)-xrefe(2),2)+pow(xrefe(5)-xrefe(3),2) ) -1  )<10e-7){
      theta_=0;
    }
    else
      std::cout<<"\n Fehler bei der Winkelberechnung";
  }


    //cross product to determine the algebraic sign of the reference angle
    if (( (xrefe(2)-xrefe(0))*(xrefe(5)-xrefe(3)) - (xrefe(3)-xrefe(1))*(xrefe(4)-xrefe(2)) ) <0)
       theta_*=-1;

  }

  return;
}


int DRT::ELEMENTS::Torsion2Type::Initialize(DRT::Discretization& dis)
{
  //reference node positions
  LINALG::Matrix<6,1> xrefe;

  //setting beam reference director correctly
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {
    //in case that current element is not a torsion2 element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    //if we get so far current element is a beam3 element and  we get a pointer at it
    DRT::ELEMENTS::Torsion2* currele = dynamic_cast<DRT::ELEMENTS::Torsion2*>(dis.lColElement(i));
    if (!currele) dserror("cast to Torsion2* failed");

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k=0; k<3; ++k) //element has three nodes
        for(int l= 0; l < 2; ++l)
          xrefe(k*2 + l) = currele->Nodes()[k]->X()[l];
    }

    currele->SetUpReferenceGeometry(xrefe);


  } //for (int i=0; i<dis_.NumMyColElements(); ++i)


  return 0;
}


