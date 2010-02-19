/*!----------------------------------------------------------------------
\file torsion3.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TORSION3
#ifdef CCADISCRET

#include "torsion3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(int id, int owner) :
DRT::Element(id,element_torsion3,owner),
data_(),
isinit_(false),
theta_(0.0),
springconstant_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(const DRT::ELEMENTS::Torsion3& old) :
 DRT::Element(old),
 data_(old.data_),
 isinit_(old.isinit_),
 theta_(old.theta_),
 springconstant_(old.springconstant_)
{
  return;
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Torsion3 and return pointer to it (public)|
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Torsion3::Clone() const
{
  DRT::ELEMENTS::Torsion3* newelement = new DRT::ELEMENTS::Torsion3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::~Torsion3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Print(ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Torsion3Register (public)             cyron 02/10|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Torsion3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Torsion3Register(Type()));
}


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Torsion3::Shape() const
{
  return line3;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  AddtoPack(data,isinit_);
  AddtoPack(data,theta_);
  AddtoPack(data,springconstant_);
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  ExtractfromPack(position,data,isinit_);
  ExtractfromPack(position,data,theta_);
  ExtractfromPack(position,data,springconstant_);
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Torsion3::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}

void DRT::ELEMENTS::Torsion3::SetUpReferenceGeometry(const LINALG::Matrix<9,1>& xrefe)
{   
  /*this method initialized geometric variables of the element; such an initialization can only be done one time when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this 
   * method will initialize the geometric variables iff the class variable isinit_ == false and afterwards set this variable to 
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/ 
  if(!isinit_)
  {
    isinit_ = true;
    
    //reference length of vectors 1-->2  und 2-->3
    LINALG::Matrix<2,1> lrefe;
    for (int j=0; j<2; ++j)
      lrefe(j) = sqrt( pow(xrefe(3+j*3)-xrefe(0+j*3),2) + pow(xrefe(4+j*3)-xrefe(1+j*3),2) + pow(xrefe(5+j*3)-xrefe(2+j*3),2) );
      
    //reference angle theta
    double dotprod=0.0;
    for (int j=0; j<3; ++j)
      dotprod +=  (xrefe(3+j)-xrefe(0+j)) * (xrefe(6+j)-xrefe(3+j));    
    
    theta_=acos(dotprod/lrefe(0)/lrefe(1));
    
  if( (dotprod/lrefe(0)/lrefe(1)) > 1){
    if( (dotprod/lrefe(0)/lrefe(1)-1)<10e-7){
      theta_=0;
    }
    else
      std::cout<<"\n error in computation of angle";
  }   

  }
 
  return;
}




//------------- class Torsion3Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::Torsion3Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::Torsion3Register(
                               const DRT::ELEMENTS::Torsion3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register* DRT::ELEMENTS::Torsion3Register::Clone() const
{
  return new DRT::ELEMENTS::Torsion3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Pack(vector<char>& data) const
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
 |  Unpack data (public)                                      cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3Register::~Torsion3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3Register::Print(ostream& os) const
{
  os << "Torsion3Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Torsion3Register::Initialize(DRT::Discretization& dis)
{     
  //reference node positions
  LINALG::Matrix<9,1> xrefe;

  //setting beam reference director correctly
  for (int i=0; i<  dis.NumMyColElements(); ++i)
  {    
    //in case that current element is not a beam3 element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(i)->Type() != DRT::Element::element_torsion3) continue;
    
    //if we get so far current element is a beam3 element and  we get a pointer at it
    DRT::ELEMENTS::Torsion3* currele = dynamic_cast<DRT::ELEMENTS::Torsion3*>(dis.lColElement(i));
    if (!currele) dserror("cast to Torsion3* failed");
    
    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {   
      for (int k=0; k<3; k++) //element has three nodes
        for(int l= 0; l < 3; l++)
          xrefe(k*3 + l) = currele->Nodes()[k]->X()[l];
    }
 
    currele->SetUpReferenceGeometry(xrefe);

  
  } //for (int i=0; i<dis_.NumMyColElements(); ++i)

  
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION3

