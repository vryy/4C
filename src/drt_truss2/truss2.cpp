/*!----------------------------------------------------------------------
\file truss2.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TRUSS2
#ifdef CCADISCRET

#include "truss2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2::Truss2(int id, int owner) :
DRT::Element(id,element_truss2,owner),
data_(),
isinit_(false),
material_(0),
lrefe_(0),
crosssec_(0),
kintype_(tr3_totlag),

//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
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
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss2::Clone() const
{
  DRT::ELEMENTS::Truss2* newelement = new DRT::ELEMENTS::Truss2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2::~Truss2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Print(ostream& os) const
{
  os << "Truss2 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Truss2Register (public)               cyron 08/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Truss2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Truss2Register(Type()));
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
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
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
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  //whether element has already been initialized
  ExtractfromPack(position,data,isinit_);
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
  ExtractfromPack(position,data,kintype_);
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Truss2::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}

void DRT::ELEMENTS::Truss2::SetUpReferenceGeometry(const LINALG::Matrix<6,1>& xrefe)
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
    lrefe_ = pow(pow(X_(3)-X_(0),2)+pow(X_(4)-X_(1),2)+pow(X_(5)-X_(2),2),0.5);
  }

  return;
}




//------------- class Truss2Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2Register::Truss2Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2Register::Truss2Register(
                               const DRT::ELEMENTS::Truss2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2Register* DRT::ELEMENTS::Truss2Register::Clone() const
{
  return new DRT::ELEMENTS::Truss2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*-----------------------------------------------------------------------*
 |  Unpack data (public)                                      cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss2Register::~Truss2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss2Register::Print(ostream& os) const
{
  os << "Truss2Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Truss2Register::Initialize(DRT::Discretization& dis)
{		
  //reference node positions
  LINALG::Matrix<6,1> xrefe;

  //setting beam reference director correctly
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {
    //in case that current element is not a beam3 element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(i)->Type() != DRT::Element::element_truss2) continue;

    //if we get so far current element is a beam3 element and  we get a pointer at it
    DRT::ELEMENTS::Truss2* currele = dynamic_cast<DRT::ELEMENTS::Truss2*>(dis.lColElement(i));
    if (!currele) dserror("cast to Truss2* failed");

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k=0; k<2; k++) //element has two nodes
        for(int l= 0; l < 3; l++)
          xrefe(k*3 + l) = currele->Nodes()[k]->X()[l];
    }

    currele->SetUpReferenceGeometry(xrefe);


  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

	
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS2
