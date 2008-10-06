/*!----------------------------------------------------------------------
\file beam3.cpp
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
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(int id, int owner) :
DRT::Element(id,element_beam3,owner),
data_(),
material_(0),
lrefe_(0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
zeta_(0),
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
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 curvconv_(old.curvconv_),
 curvold_(old.curvold_),
 curvnew_(old.curvnew_),
 betaplusalphaconv_(old.betaplusalphaconv_),
 betaplusalphaold_(old.betaplusalphaold_),
 betaplusalphanew_(old.betaplusalphanew_),
 betaminusalphaconv_(old.betaminusalphaconv_),
 betaminusalphaold_(old.betaminusalphaold_),
 betaminusalphanew_(old.betaminusalphanew_),
 crosssec_(old.crosssec_),
 crosssecshear_(old.crosssecshear_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 zeta_(old.zeta_),
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
  AddtoPack(data,Qconv_);
  AddtoPack(data,Qold_);
  AddtoPack(data,Qnew_);
  AddtoPack(data,curvconv_);
  AddtoPack(data,curvold_);
  AddtoPack(data,curvnew_);  
  AddtoPack(data,betaplusalphaconv_);
  AddtoPack(data,betaplusalphaold_);
  AddtoPack(data,betaplusalphanew_);
  AddtoPack(data,betaminusalphaconv_);
  AddtoPack(data,betaminusalphaold_);
  AddtoPack(data,betaminusalphanew_);
  //cross section
  AddtoPack(data,crosssec_);
   //cross section with shear correction
  AddtoPack(data,crosssecshear_);
  //moments of inertia of area
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,Irr_);
  //viscosity of surrounding fluid
  AddtoPack(data,zeta_);
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
  ExtractfromPack(position,data,Qconv_);
  ExtractfromPack(position,data,Qold_);
  ExtractfromPack(position,data,Qnew_);
  ExtractfromPack(position,data,curvconv_);
  ExtractfromPack(position,data,curvold_);
  ExtractfromPack(position,data,curvnew_); 
  ExtractfromPack(position,data,betaplusalphaconv_);
  ExtractfromPack(position,data,betaplusalphaold_);
  ExtractfromPack(position,data,betaplusalphanew_);
  ExtractfromPack(position,data,betaminusalphaconv_);
  ExtractfromPack(position,data,betaminusalphaold_);
  ExtractfromPack(position,data,betaminusalphanew_);   
  //cross section
  ExtractfromPack(position,data,crosssec_);
  //cross section with shear correction
  ExtractfromPack(position,data,crosssecshear_);
  //moments of inertia of area
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,Irr_);
  //viscosity of surrounding fluid
  ExtractfromPack(position,data,zeta_);
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam3::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
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


  //variable for nodal point coordinates in reference configuration
  BlitzMat3x2 xrefe;
  //center triad in reference configuration
  BlitzMat3x3 Tref;
  
  //setting beam reference director correctly
  for (int num=0; num<  dis.NumMyColElements(); ++num)
    {    
      //in case that current element is not a beam3 element there is nothing to do and we go back
      //to the head of the loop
      if (dis.lColElement(num)->Type() != DRT::Element::element_beam3) continue;
      
      //if we get so far current element is a beam3 element and  we get a pointer at it
      DRT::ELEMENTS::Beam3* currele = dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(num));
      if (!currele) dserror("cast to Beam3* failed");
      
      //getting element's reference coordinates     
      for (int k=0; k<2; ++k) //element has two nodes
        {
          xrefe(0,k) = currele->Nodes()[k]->X()[0];
          xrefe(1,k) = currele->Nodes()[k]->X()[1];
          xrefe(2,k) = currele->Nodes()[k]->X()[2];
        }
      
      //length in reference configuration
      currele->lrefe_ = pow(pow(xrefe(0,1)-xrefe(0,0),2)+pow(xrefe(1,1)-xrefe(1,0),2)+pow(xrefe(2,1)-xrefe(2,0),2),0.5);  
      
      /*initial triad Tref = [t1,t2,t3] with t1-axis equals beam axis and t2 and t3 are principal axes of the moment of inertia 
       * of area tensor */
      for (int k=0; k<3; k++)	
      {
        Tref(k,0) = ( xrefe(k,1)-xrefe(k,0) )/currele->lrefe_;        
      }
      /*in the following two more or less arbitrary axes t2 and t3 are calculated in order to complete the triad Told_, which 
       * works in case of a rotationally symmetric crosssection; in case of different kinds of crosssections one has to mo-
       * dify the following code lines in such a way that t2 and t3 are still the principal axes related with Iyy_ and Izz_*/
      
      //t2 is calculated as a unit vector in the x2x3-plane orthogonal to t1
      Tref(0,1) = 0;
      //if t1 is a unit direction vector for the x1-axis t2 is set to a unit direction vector of the x2-axis
      if (Tref(1,0) == 0 && Tref(2,0) == 0)
      {    
        Tref(1,1) = 1;
        Tref(2,1) = 0;
      }
      //otherwise t2 is calculated from the scalar product with t1
      else
      { 
        //seeting t2(0)=0 and calculating other elements by setting scalar product t1 o t2 to zero
        double lin1norm = pow(pow(Tref(1,0),2)+pow(Tref(2,0),2),0.5);
        Tref(1,1) =  Tref(1,0)/lin1norm;
        Tref(2,1) = -Tref(2,0)/lin1norm;
      }
   
      //calculating t3 by crossproduct t1 x t2
      Tref(0,2) = -Tref(1,0)*Tref(2,1)-Tref(2,0)*Tref(1,1);
      Tref(1,2) =  Tref(0,0)*Tref(2,1);
      Tref(2,2) =  Tref(0,0)*Tref(1,1);
    
      /*the center triad in reference configuration is stored as a quaternion whose equivalent would be the rotation
       * from the identity matrix into the reference configuration all operations are performed according to
       * Crisfield, Vol. 2, section 16.10 and the there described Spurrier's algorithm*/
      double trace = Tref(0,0) + Tref(1,1) + Tref(2,2);
      if(trace>Tref(0,0)  && trace>Tref(1,1) &&trace>Tref(2,2))
      {
        currele->Qconv_(3) = 0.5 * pow(1 + trace, 0.5);
        currele->Qconv_(0) = (Tref(1,2) - Tref(2,1)) / (4*currele->Qconv_(3));
        currele->Qconv_(1) = (Tref(2,0) - Tref(0,2)) / (4*currele->Qconv_(3));
        currele->Qconv_(2) = (Tref(0,1) - Tref(1,0)) / (4*currele->Qconv_(3));
      }
      else if(Tref(0,0)>Tref(1,1) && Tref(0,0)>Tref(2,2))
      {
        for(int i = 0 ; i<3 ; i++)
        {
          int k = (i+1)% 3;
          int j = (i+2)% 3;
          int l;
          
          currele->Qconv_(i) = pow(0.5*Tref(i,i) + 0.25*(1 - trace) , 0.5);
          
          l = j;        
          currele->Qconv_(l) = 0.25*(Tref(l,i) + Tref(i,l)) / currele->Qconv_(i);
          
          l = k;        
          currele->Qconv_(l) = 0.25*(Tref(l,i) + Tref(i,l)) / currele->Qconv_(i);
          
          currele->Qconv_(3) = 0.25*(Tref(k,j) - Tref(j,k)) / currele->Qconv_(i);
        }
      }
  
      currele->Qold_ = currele->Qconv_;
      currele->Qnew_ = currele->Qconv_;
   
      //the here employed beam element does not need data about the current position of the nodal directors so that
      //initilization of those can be skipped (the nodal displacements handeled in beam3_evaluate.cpp are not the actual angles,
      //but only the differences between actual angles and angles in reference configuration, respectively. Thus the
      //director orientation in reference configuration cancels out and can be assumed to be zero without loss of 
      //generality
      for (int k=0; k<3; k++) 
      {
        currele->curvconv_(k) = 0;
        currele->curvold_(k) = 0;
        currele->curvnew_(k) = 0;
        currele->betaplusalphaconv_(k)  = 0;
        currele->betaplusalphaold_(k)  = 0;   
        currele->betaplusalphanew_(k)  = 0;
        currele->betaminusalphaconv_(k) = 0;
        currele->betaminusalphaold_(k) = 0;
        currele->betaminusalphaold_(k) = 0;
      }
            
    } //for (int num=0; num<dis_.NumMyColElements(); ++num)
	
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
