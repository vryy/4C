/*!----------------------------------------------------------------------
\file beam2.cpp
\brief two dimensional nonlinear corotational Timoshenko beam element

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Beam2Type DRT::ELEMENTS::Beam2Type::instance_;

DRT::ParObject* DRT::ELEMENTS::Beam2Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Beam2* object = new DRT::ELEMENTS::Beam2(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam2Type::Create( const std::string eletype,
                                                             const std::string eledistype,
                                                             const int         id,
                                                             const int         owner )
{
  if ( eletype=="BEAM2" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam2(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Beam2Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Beam2(id,owner));
  return ele;
}


void DRT::ELEMENTS::Beam2Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 3;
  nv = 3;
}

void DRT::ELEMENTS::Beam2Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeBeam2DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Beam2Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["BEAM2"];

  defs["LINE2"]
    .AddIntVector("LINE2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;

  defs["LIN2"]
    .AddIntVector("LIN2",2)
    .AddNamedInt("MAT")
    .AddNamedDouble("CROSS")
    .AddNamedDouble("SHEARCORR")
    .AddNamedDouble("INERMOM")
    ;
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::Beam2(int id, int owner) :
DRT::Element(id,owner),
isinit_(false),
lrefe_(0),
crosssec_(0),
crosssecshear_(0),
mominer_(0),
numperiodsnew_(0),
numperiodsold_(0),
numperiodsconv_(0),
alphanew_(0),
alphaold_(0),
alphaconv_(0),
alpha0_(0),
//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::Beam2(const DRT::ELEMENTS::Beam2& old) :
DRT::Element(old),
isinit_(old.isinit_),
lrefe_(old.lrefe_),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
mominer_(old.mominer_),
numperiodsnew_(old.numperiodsnew_),
numperiodsold_(old.numperiodsold_),
numperiodsconv_(old.numperiodsconv_),
alphanew_(old.alphanew_),
alphaold_(old.alphaold_),
alphaconv_(old.alphaconv_),
alpha0_(old.alpha0_),
jacobi_(old.jacobi_),
jacobimass_(old.jacobimass_),
jacobinode_(old.jacobinode_),
gaussrule_(old.gaussrule_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam2 and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam2::Clone() const
{
  DRT::ELEMENTS::Beam2* newelement = new DRT::ELEMENTS::Beam2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::~Beam2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Print(std::ostream& os) const
{
  os << "Beam2 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam2::Shape() const
{
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  AddtoPack(data,isinit_);
  AddtoPack(data,lrefe_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack(data,mominer_);
  AddtoPack(data,numperiodsnew_);
  AddtoPack(data,numperiodsold_);
  AddtoPack(data,numperiodsconv_);
  AddtoPack(data,alphanew_);
  AddtoPack(data,alphaold_);
  AddtoPack(data,alphaconv_);
  AddtoPack(data,alpha0_);
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobimass_);
  AddtoPack(data,jacobinode_);

  // gaussrule_
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Unpack(const std::vector<char>& data)
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
  ExtractfromPack(position,data,lrefe_);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  ExtractfromPack(position,data,mominer_);
  ExtractfromPack(position,data,numperiodsnew_);
  ExtractfromPack(position,data,numperiodsold_);
  ExtractfromPack(position,data,numperiodsconv_);
  ExtractfromPack(position,data,alphanew_);
  ExtractfromPack(position,data,alphaold_);
  ExtractfromPack(position,data,alphaconv_);
  ExtractfromPack(position,data,alpha0_);
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);

  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::Beam2::Lines()
{
  std::vector<Teuchos::RCP<Element> > lines(1);
  lines[0]= Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam2::MyGaussRule(int nnode, IntegrationType integrationtype)
{
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule1D_undefined;

  switch(integrationtype)
  {
    case gaussexactintegration:
    {
      gaussrule = DRT::UTILS::intrule_line_2point;
      break;
    }
    case gaussunderintegration:
    {
      gaussrule =  DRT::UTILS::intrule_line_1point;
      break;
    }
    case lobattointegration:
    {
      gaussrule =  DRT::UTILS::intrule_line_lobatto2point;
      break;
    }
    default:
      dserror("unknown type of integration");
  }

  return gaussrule;
}

//sets up element reference geomtry for reference nodal position vector xrefe (may be used also after simulation start)
void DRT::ELEMENTS::Beam2::SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe)
{
  /*this method initializes geometric variables of the element; such an initialization can only be done one time when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this
   * method will initialize the geometric variables iff the class variable isinit_ == false and afterwards set this variable to
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/
  if(!isinit_)
  {
    isinit_ = true;

    //length in reference configuration
    lrefe_  = std::pow( pow(xrefe(3)-xrefe(1),2) + pow(xrefe(2)-xrefe(0),2) , 0.5 );

    //resize jacobi_, jacobimass_, jacobinode_  so they can each store for each Gauss point jacobi determinant lrefe_ / 2.0

    //underintegration in elasticity -> nnode - 1 Gauss points required
    jacobi_.resize(1,lrefe_ / 2.0);
    //exact integration of mass matrix -> nnode Gauss point required
    jacobimass_.resize(2,lrefe_ / 2.0);
    //for nodal quadrature as many Gauss points as nodes required
    jacobinode_.resize(2,lrefe_ / 2.0);


    // beta is the rotation angle out of x-axis in a x-y-plane in reference configuration
    double cos_alpha0 = (xrefe(2)-xrefe(0))/lrefe_;
    double sin_alpha0 = (xrefe(3)-xrefe(1))/lrefe_;

    //we calculate beta in a range between -pi < beta <= pi
    if (cos_alpha0 >= 0)
      alpha0_ = asin(sin_alpha0);
    else
    { if (sin_alpha0 >= 0)
        alpha0_ =  acos(cos_alpha0);
      else
        alpha0_ = -acos(cos_alpha0);
     }

    //initially the absolute rotation of the element frame equals the reference rotation alpha0_
    alphanew_  = alpha0_;
    alphaold_  = alpha0_;
    alphaconv_ = alpha0_;

    /*the angle alpha0_ is exactly the angle beta gained from the coordinate positions by evaluation
     * of the sine- and cosine-functions without adding or substracting any multiple of 2*PI*/
    numperiodsnew_  = 0;
    numperiodsold_  = 0;
    numperiodsconv_ = 0;

  }

  return;
} //DRT::ELEMENTS::Beam3::SetUpReferenceGeometry()



int DRT::ELEMENTS::Beam2Type::Initialize(DRT::Discretization& dis)
{

  //reference node position
  LINALG::Matrix<4,1> xrefe;

  //setting up geometric variables for beam3 elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {
    //in case that current element is not a beam2 element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->ElementType() != *this) continue;

    //if we get so far current element is a beam2 element and  we get a pointer at it
    DRT::ELEMENTS::Beam2* currele = dynamic_cast<DRT::ELEMENTS::Beam2*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam2* failed");

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

  } //for (int num=0; num<dis_.NumMyColElements(); ++num)

  return 0;
}


