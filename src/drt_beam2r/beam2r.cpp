/*!----------------------------------------------------------------------
\file beam2.cpp
\brief two dimensional nonlinear beam element using Reissner`s theory.
\According to Crisfield Non-linear finite element analysis of solids and structures Vol.1 section 7.4
<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2R
#ifdef CCADISCRET

#include "beam2r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(int id, int owner) :
DRT::Element(id,element_beam2r,owner),
data_(),
isinit_(false),
lrefe_(0),
crosssec_(0),
crosssecshear_(0),
mominer_(0),
thetaav0_(0),
//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(const DRT::ELEMENTS::Beam2r& old) :
DRT::Element(old),
data_(old.data_),
isinit_(old.isinit_),
lrefe_(old.lrefe_),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
mominer_(old.mominer_),
thetaav0_(old.thetaav0_),
gaussrule_(old.gaussrule_)

{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam2r and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam2r::Clone() const
{
  DRT::ELEMENTS::Beam2r* newelement = new DRT::ELEMENTS::Beam2r(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::~Beam2r()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Print(ostream& os) const
{
  os << "Beam2r ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Beam2rRegister (public)               cyron 01/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Beam2r::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Beam2rRegister(Type()));
}



/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam2r::Shape() const
{
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Pack(vector<char>& data) const
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
  //reference length
  AddtoPack(data,lrefe_);
  //cross section
  AddtoPack(data,crosssec_);
   //cross section with shear correction
  AddtoPack(data,crosssecshear_);
  //moment of inertia of area
  AddtoPack(data,mominer_);
  // average absolute angle in reference configuration
  AddtoPack(data,thetaav0_);
  // gaussrule_
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::Unpack(const vector<char>& data)
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
  //reference length
  ExtractfromPack(position,data,lrefe_);
  //cross section
  ExtractfromPack(position,data,crosssec_);
   //cross section with shear correction
  ExtractfromPack(position,data,crosssecshear_);
  //moment of inertia of area
  ExtractfromPack(position,data,mominer_);
  //average absolute angle in reference configuration
  ExtractfromPack(position,data,thetaav0_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam2r::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}

//sets up element reference geomtry for reference nodal position vector xrefe (may be used also after simulation start)
void DRT::ELEMENTS::Beam2r::SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe)
{
  /*this method initializes geometric variables of the element; such an initialization can only be done once when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this 
   * method will initialize the geometric variables if the class variable isinit_ == false and afterwards set this variable to 
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/ 
  if(!isinit_)
  {
    isinit_ = true;

    //length in reference configuration
    lrefe_  = pow( pow(xrefe(3)-xrefe(1),2) + pow(xrefe(2)-xrefe(0),2) , 0.5 );

    /*note: thetaav0 is here chosen as if the beam has no curvature in reference configuration.
     * This is an assumption but we can prove that this leads to no stresses in reference configuration. 
     * The current thetaav can later be computed from thataav0 and the increments of theta1 and theta2 */
    //thetaav0 in reference configuration
    double cos_thetaav0 = (xrefe(2)-xrefe(0))/lrefe_;
    double sin_thetaav0 = (xrefe(3)-xrefe(1))/lrefe_;

    //we calculate thetaav0 in a range between -pi < thetaav0 <= pi, Crisfield Vol. 1 (7.60)
    if (cos_thetaav0 >= 0)
      thetaav0_ = asin(sin_thetaav0);
    else
    { if (sin_thetaav0 >= 0)
    	thetaav0_ =  acos(cos_thetaav0);
      else
    	  thetaav0_ = -acos(cos_thetaav0);
    }
    
  }

  return;
} //DRT::ELEMENTS::Beam2r::SetUpReferenceGeometry()



//------------- class Beam2rRegister: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2rRegister::Beam2rRegister(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2rRegister::Beam2rRegister(
                               const DRT::ELEMENTS::Beam2rRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2rRegister* DRT::ELEMENTS::Beam2rRegister::Clone() const
{
  return new DRT::ELEMENTS::Beam2rRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2rRegister::Pack(vector<char>& data) const
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
 |  Unpack data (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2rRegister::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2rRegister::~Beam2rRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2rRegister::Print(ostream& os) const
{
  os << "Beam2rRegister ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Beam2rRegister::Initialize(DRT::Discretization& dis)
{
  
  //reference node position
  LINALG::Matrix<4,1> xrefe;
 
  //setting up geometric variables for beam2r elements
  for (int num=0; num<  dis.NumMyColElements(); ++num)
  {    
    //in case that current element is not a beam2r element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->Type() != DRT::Element::element_beam2r) continue;
    
    //if we get so far current element is a beam2r element and  we get a pointer at it
    DRT::ELEMENTS::Beam2r* currele = dynamic_cast<DRT::ELEMENTS::Beam2r*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam2r* failed");
    
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


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R
