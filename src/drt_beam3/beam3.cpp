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
#ifdef D_BEAM3
#ifdef CCADISCRET

#include "beam3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(int id, int owner) :
DRT::Element(id,element_beam3,owner),
data_(),
material_(0),
lrefe_(0),
timenew_(0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
lumpedflag_(0),
thermalenergy_(0),
Arbeit_N(0),
Arbeit_M(0),
Arbeit_Q(0),
x_verschiebung(0),
//since lines_ is a vector it calls its constructor automatically -> 
//no initialization of lines_ necessary here

//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(const DRT::ELEMENTS::Beam3& old) :
DRT::Element(old),
data_(old.data_),
material_(old.material_),
lrefe_(old.lrefe_),
Told_(old.Told_),
Tnew_(old.Tnew_),
Tmid_(old.Tmid_),
curvold_(old.curvold_),
curvnew_(old.curvnew_),
alphaold_(old.alphaold_),
alphanew_(old.alphanew_),
timenew(old,timenew),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
Iyy_(old.Iyy_),
Izz_(old.Izz_),
Irr_(old.Irr_),
lumpedflag_(old.lumpedflag_),
thermalenergy_(old.thermalenergy_),
Arbeit_N(old.Arbeit_N),
Arbeit_M(old.Arbeit_M),
Arbeit_Q(old.Arbeit_Q),
x_verschiebung(old.x_verschiebung),
lines_(old.lines_),
gaussrule_(old.gaussrule_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3 and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3::Clone() const
{
  DRT::ELEMENTS::Beam3* newelement = new DRT::ELEMENTS::Beam3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::~Beam3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::Print(ostream& os) const
{
  os << "Beam3 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Beam3Register (public)               cyron 01/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Beam3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Beam3Register(Type()));
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3::Shape() const
{
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  //material type
  AddtoPack(data,material_);
  //reference length
  AddtoPack(data,lrefe_);
  //central coordinate triad and related data
  AddtoPack(data,Told_);
  AddtoPack(data,Tnew_);
  AddtoPack(data,Tmid_);
  AddtoPack(data,curvold_);
  AddtoPack(data,curvnew_);
  AddtoPack(data,alphaold_);
  AddtoPack(data,alphanew_);
  AddtoPack(data,timenew);
  //cross section
  AddtoPack(data,crosssec_);
   //cross section with shear correction
  AddtoPack(data,crosssecshear_);
  //moments of inertia of area
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,Irr_);
  //flag determining if consistent or lumped mass matrix
  AddtoPack(data,lumpedflag_);
  //thermal energy responsible for statistical forces
  AddtoPack(data,thermalenergy_);
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
void DRT::ELEMENTS::Beam3::Unpack(const vector<char>& data)
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
  //material type
  ExtractfromPack(position,data,material_);
  //reference length
  ExtractfromPack(position,data,lrefe_);
  //central coordinate triad and related data
  ExtractfromPack(position,data,Told_);
  ExtractfromPack(position,data,Tnew_);
  ExtractfromPack(position,data,Tmid_);
  ExtractfromPack(position,data,curvold_);
  ExtractfromPack(position,data,curvnew_);
  ExtractfromPack(position,data,alphaold_);
  ExtractfromPack(position,data,alphanew_);
  ExtractfromPack(position,data,timenew);
  //cross section
  ExtractfromPack(position,data,crosssec_);
  //cross section with shear correction
  ExtractfromPack(position,data,crosssecshear_);
  //moments of inertia of area
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,Irr_);
  //flag determining if consistent or lumped mass matrix
  ExtractfromPack(position,data,lumpedflag_);
  //thermal energy responsible for statistical forces
  ExtractfromPack(position,data,thermalenergy_);
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
DRT::Element** DRT::ELEMENTS::Beam3::Lines()
{
  lines_.resize(1);
  lines_[0] = this;
  return &lines_[0];
}




//------------- class Beam3Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Register::Beam3Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Register::Beam3Register(
                               const DRT::ELEMENTS::Beam3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Register* DRT::ELEMENTS::Beam3Register::Clone() const
{
  return new DRT::ELEMENTS::Beam3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Beam3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Beam3Register::~Beam3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Register::Print(ostream& os) const
{
  os << "Beam3Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Beam3Register::Initialize(DRT::Discretization& dis)
{		
  //random generator for seeding only (necessary for thermal noise)
  ranlib::Normal<double> seedgenerator(0,1);
  seedgenerator.seed((unsigned int)std::time(0));
  
  //variable for nodal point coordinates in reference configuration
  LINALG::SerialDenseMatrix xrefe;
  xrefe.Shape(3,2);
  
  //setting beam reference director correctly
  for (int i=0; i<dis.NumMyColElements(); ++i)
    {
      //in case that current element is not a beam3 element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(i)->Type() != DRT::Element::element_beam2) continue;
      
      //if we get so far current element is a beam3 element and  we get a pointer at it
      DRT::ELEMENTS::Beam3* currele = dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(i));
      if (!currele) dserror("cast to Beam3* failed");
      
      //getting element's reference coordinates     
      for (int k=0; k<2; ++k) //element has two nodes
        {
          xrefe(0,k) = currele->Nodes()[k]->X()[0];
          xrefe(1,k) = currele->Nodes()[k]->X()[1];
          xrefe(2,k) = currele->Nodes()[k]->X()[2];
        }
      //Initializing member matrices
      Told_.Shape(3,3);
      Tnew_.Shape(3,3);
      Tmid_.Shape(3,3);
      alphaold_.Shape(3,1);
      alphanew_.Shape(3,1);
      curvold_.Shape(3,1);
      curvnew_.Shape(3,1);
      
      //length in reference configuration
      currele->lrefe_ = pow(pow(xrefe(0,1)-xrefe(0,0),2)+pow(xrefe(1,1)-xrefe(1,0),2)+pow(xrefe(2,1)-xrefe(2,0),2),0.5);  
      
      /*initial triad Told_ = [t1,t2,t3] with t1-axis equals beam axis and t2 and t3 are principal axes of the moment of inertia 
       * of area tensor */
      for (int line; line<3; line++)	
      {
	      Told_(line,0) = xrefe(line,1)-xrefe(line,0)/currele->lrefe_;
      }
      /*in the following two more or less arbitrary axes t2 and t3 are calculated in order to complete the triad Told_, which 
       * works in case of a rotationally symmetric crosssection; in case of different kinds of crosssections one has to mo-
       * dify the following code lines in such a way that t2 and t3 are still the principal axes related with Iyy_ and Izz_*/
      
      //seeting t2(0)=0 and calculating other elements by setting scalar product t1 o t2 to zero
      double lin1norm = pow(pow(Told_(1,0),2)+pow(Told_(2,0),2),0.5);
      Told_(0,1) =  0;
      Told_(1,1) =  Told_(1,0)/lin1norm;
      Told_(2,1) = -Told_(2,0)/lin1norm;
      //calculating t3 by crossproduct t1 x t2
      Told_(0,2) = -Told_(1,0)*Told_(2,1)-Told_(2,0)*Told_(1,1);
      Told_(1,2) =    Told_(0,0)*Told_(2,1);
      Told_(2,2) =    Told_(0,0)*Told_(1,1);
       
      
    } //for (int i=0; i<dis_.NumMyColElements(); ++i)
	
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
