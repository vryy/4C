/*!----------------------------------------------------------------------
\file beam2.cpp
\brief two dimensional nonlinear corotational Timoshenko beam element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

#include "beam2.H"
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
DRT::ELEMENTS::Beam2::Beam2(int id, int owner) :
DRT::Element(id,element_beam2,owner),
data_(),
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
data_(old.data_),
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
void DRT::ELEMENTS::Beam2::Print(ostream& os) const
{
  os << "Beam2 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Beam2Register (public)               cyron 01/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Beam2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Beam2Register(Type()));
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
void DRT::ELEMENTS::Beam2::Pack(vector<char>& data) const
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
  //number of 2*PI rotations of element frame compared to angle gained from sine- and cosine evaluatoin
  AddtoPack(data,numperiodsnew_);
  AddtoPack(data,numperiodsold_);
  AddtoPack(data,numperiodsconv_);
  //absolute angle of element frame
  AddtoPack(data,alphanew_);
  AddtoPack(data,alphaold_);
  AddtoPack(data,alphaconv_);
  //angle relative to x-axis in reference configuration
  AddtoPack(data,alpha0_);
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
void DRT::ELEMENTS::Beam2::Unpack(const vector<char>& data)
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
  //number of 2*PI rotations of element frame compared to angle gained from sine- and cosine evaluatoin
  ExtractfromPack(position,data,numperiodsnew_);
  ExtractfromPack(position,data,numperiodsold_);
  ExtractfromPack(position,data,numperiodsconv_);
  //absolute angle of element frame
  ExtractfromPack(position,data,alphanew_);
  ExtractfromPack(position,data,alphaold_);
  ExtractfromPack(position,data,alphaconv_);
  //angle relative to x-axis in reference configuration
  ExtractfromPack(position,data,alpha0_);
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam2::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}



//------------- class Beam2Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register::Beam2Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register::Beam2Register(
                               const DRT::ELEMENTS::Beam2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register* DRT::ELEMENTS::Beam2Register::Clone() const
{
  return new DRT::ELEMENTS::Beam2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Beam2Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Beam2Register::~Beam2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2Register::Print(ostream& os) const
{
  os << "Beam2Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Beam2Register::Initialize(DRT::Discretization& dis)
{
  LINALG::SerialDenseMatrix xrefe;
  xrefe.Shape(2,2);

  //setting beam reference director correctly
  for (int i=0; i<dis.NumMyColElements(); ++i)
    {
      //in case that current element is not a beam2 element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(i)->Type() != DRT::Element::element_beam2) continue;

      //if we get so far current element is a beam2 element and  we get a pointer at it
      DRT::ELEMENTS::Beam2* currele = dynamic_cast<DRT::ELEMENTS::Beam2*>(dis.lColElement(i));
      if (!currele) dserror("cast to Beam2* failed");
      
      /*the following part initializes geometric variables of the element; such an initialization can only be done one time when the element is
       * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this 
       * method will initialize the geometric variables iff the class variable isinit_ == false and afterwards set this variable to 
       * isinit_ = true; if this part is called and finds alreday isinit_ == true it will just do nothing*/    
      if(!currele->isinit_)
      {
        currele->isinit_ = true;

        //getting element's reference coordinates
        for (int k=0; k<2; ++k) //element has two nodes
          {
            xrefe(0,k) = currele->Nodes()[k]->X()[0];
            xrefe(1,k) = currele->Nodes()[k]->X()[1];
          }
  
        //length in reference configuration
        currele->lrefe_  = pow( pow(xrefe(0,1)-xrefe(0,0),2) + pow(xrefe(1,1)-xrefe(1,0),2) , 0.5 );
  
        // beta is the rotation angle out of x-axis in a x-y-plane in reference configuration
        double cos_alpha0 = (xrefe(0,1)-xrefe(0,0))/currele->lrefe_;
        double sin_alpha0 = (xrefe(1,1)-xrefe(1,0))/currele->lrefe_;
  
        //we calculate beta in a range between -pi < beta <= pi
        if (cos_alpha0 >= 0)
        	currele->alpha0_ = asin(sin_alpha0);
        else
        {	if (sin_alpha0 >= 0)
            currele->alpha0_ =  acos(cos_alpha0);
          else
            currele->alpha0_ = -acos(cos_alpha0);
         }
        
        //initially the absolute rotation of the element frame equals the reference rotation alpha0_
        currele->alphanew_  = currele->alpha0_;
        currele->alphaold_  = currele->alpha0_;
        currele->alphaconv_ = currele->alpha0_;
        
        /*the angle alpha0_ is exactly the angle beta gained from the coordinate positions by evaluation
         * of the sine- and cosine-functions without adding or substracting any multiple of 2*PI*/
        currele->numperiodsnew_  = 0;
        currele->numperiodsold_  = 0;
        currele->numperiodsconv_ = 0;
      }
    } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
