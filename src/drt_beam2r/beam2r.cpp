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
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H" // for shape functions

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(int id, int owner) :
DRT::Element(id,element_beam2r,owner),
data_(),
crosssec_(0),
crosssecshear_(0),
gaussrule_(DRT::UTILS::intrule1D_undefined),
isinit_(false),
lrefe_(0),
mominer_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2r::Beam2r(const DRT::ELEMENTS::Beam2r& old) :
DRT::Element(old),
data_(old.data_),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
gaussrule_(old.gaussrule_),
isinit_(old.isinit_),
lrefe_(old.lrefe_),
mominer_(old.mominer_),
alpha_(old.alpha_),
alphamass_(old.alphamass_),
floc_(old.floc_),
theta0_(old.theta0_)
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
 |  Checks the number of nodes and returns the  			   (public) |
 |  DiscretizationType                                      cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam2r::Shape() const
{
  int numnodes = NumNode();
  switch(numnodes)
  {
  	case 2:
  			return line2;
  			break;
  	case 3:
  			return line3;
  			break;
  	case 4:
  	  		return line4;
  	  		break;
  	case 5:
  	  		return line5;
  	  		break;
  	default:
  			dserror("Only Line2, Line3 and Line4 elements are implemented.");
  			break;
  	
  }
  return dis_none;
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
  
  //add all class variables of beam2r element
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_);
  AddtoPack(data,gaussrule_);
  AddtoPack(data,isinit_);
  AddtoPack(data,lrefe_);
  AddtoPack(data,mominer_);
   
  for(int i=0; i<NumNode()-1;i++)
    AddtoPack(data,alpha_[i]);
  
  for(int i=0; i<NumNode()-1;i++)
    AddtoPack(data,alphamass_[i]);
  
  for(int i=0; i<NumNode();i++)
    AddtoPack(data,floc_[i]);
  
  for(int i=0; i<NumNode()-1;i++)
    AddtoPack(data,theta0_[i]);
   
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
 
  //extract all class variables of beam2r element
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  
  ExtractfromPack(position,data,isinit_);
  ExtractfromPack(position,data,lrefe_);
  ExtractfromPack(position,data,mominer_);
 
  alpha_.resize(NumNode()-1);
  for (int i=0;i<NumNode()-1;i++)
    ExtractfromPack(position,data,alpha_[i]);
  
  alphamass_.resize(NumNode()-1);
  for (int i=0;i<NumNode()-1;i++)
    ExtractfromPack(position,data,alphamass_[i]);
  
  floc_.resize(NumNode());
  for (int i=0;i<NumNode();i++)
    ExtractfromPack(position,data,floc_[i]);
  
  theta0_.resize(NumNode()-1);
  for (int i=0;i<NumNode()-1;i++)
    ExtractfromPack(position,data,theta0_[i]);
 
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08   |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam2r::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}


/*----------------------------------------------------------------------*
 |  sets up element reference geometry for reference nodal position 	|
 |  vector xrefe (may be used also after simulation start)  cyron 01/08 |
 *----------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam2r::SetUpReferenceGeometry(const vector<double>& xrefe)
{
  /*this method initializes geometric variables of the element; such an initialization can only be done once when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this 
   * method will initialize the geometric variables if the class variable isinit_ == false and afterwards set this variable to 
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/ 
  if(!isinit_)
  {

  	floc_.resize(nnode);//set vector size for local stochastic forces
  	
  	isinit_ = true;
    
	  //create Matrix for the derivates of the shapefunctions at the GP
	  LINALG::Matrix<1,nnode> shapefuncderiv;
      
	  //Get DiscretizationType
	  DRT::Element::DiscretizationType distype = Shape();
	  
	  //Get the applied integrationpoints
	  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule_);
	  
	  //Loop through all GPs and calculate alpha and theta0
	  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
	  {
    	
		  //Get position xi of GP
		  const double xi = gausspoints.qxg[numgp][0];
		  
		  //Get derivatives of shapefunctions at GP
		  DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);
		  			  
		  double dxdxi_gp(0.0),dzdxi_gp(0.0);
    	
		  //calculate dx/dxi and dz/dxi
		  for(int i=0; i<nnode; i++)
		  {
			  dxdxi_gp+=shapefuncderiv(i)*xrefe[2*i];
			  dzdxi_gp+=shapefuncderiv(i)*xrefe[2*i+1];
		  }//for(int i=0; i<nnode; i++)
		  //Store length factor for every GP
		  //note: the length factor alpha replaces the determinant and refers by definition always to the reference configuration
		  alpha_[numgp]= pow(pow( dxdxi_gp ,2.0) + pow( dzdxi_gp ,2.0) ,0.5);
    	
		  /*calculate sin and cos theta0 for each gausspoint for a stress-free-reference-configuration
		   *the formulas are derived from Crisfield Vol.1 (7.132) and (7.133) for no strain
		   * Therfore we assume that we have no strain at each GP in ref. config. 
		   */
		  double cos_theta0 = alpha_[numgp]*dxdxi_gp/(pow(dxdxi_gp,2)+pow(dzdxi_gp,2));
    	
		  double sin_theta0 = alpha_[numgp]*dzdxi_gp/(pow(dxdxi_gp,2)+pow(dzdxi_gp,2));
    	  
		  //we calculate thetaav0 in a range between -pi < thetaav0 <= pi, Crisfield Vol. 1 (7.60)
		  if (cos_theta0 >= 0)
			  theta0_[numgp] = asin(sin_theta0);
		  else
		  {   
			  if (sin_theta0 >= 0)
				  theta0_[numgp] = acos(cos_theta0);
			  else
				  theta0_[numgp] = -acos(cos_theta0);
		  }
		  
		  /* Here we force the triad to point in positive xi direction at every gausspoint.
		   * We compare vector t1 [Crisfield Vol. 1 (7.110)] with the tangent on our element 
		   * in positive xi direction via scalar product. If the difference is more than +/-90 degrees
		   * we turn the triad with 180 degrees. Our stress free reference configuration is
		   * unaffected by this modification.
		   */
		  double check = 0.0;
		  
		  check = (cos(theta0_[numgp])*dxdxi_gp+sin(theta0_[numgp])*dzdxi_gp);
		  
		  if (check<0)
		  {
			  if(theta0_[numgp]> 0 )
			  theta0_[numgp]-=3.141592653589;
			  else
		      theta0_[numgp]+=3.141592653589;	  
		  }
	  
	  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)
	  
	  //Now we get the integrationfactor alphamass_ for a complete integration of the massmatrix
	  
	  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode);
	  
	  //Get the applied integrationpoints
	  DRT::UTILS::IntegrationPoints1D gausspointsmass(gaussrule_);
	  
	  //Loop through all GPs and calculate alpha and theta0
	  for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
	  {
    	
		  //Get position xi of GP
		  const double xi = gausspointsmass.qxg[numgp][0];
		  
		  //Get derivatives of shapefunctions at GP
		  DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);
		  			  
		  double dxdxi_gp(0.0),dzdxi_gp(0.0);
    	
		  //calculate dx/dxi and dz/dxi
		  for(int i=0; i<nnode; i++)
		  {
			  dxdxi_gp+=shapefuncderiv(i)*xrefe[2*i];
			  dzdxi_gp+=shapefuncderiv(i)*xrefe[2*i+1];
		  }//for(int i=0; i<nnode; i++)
		  //Store length factor for every GP
		  //note: the length factor alpha replaces the determinant and refers by definition always to the reference configuration
		  alphamass_[numgp]= pow(pow( dxdxi_gp ,2.0) + pow( dzdxi_gp ,2.0) ,0.5);
    	
	  }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
	 
	  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);
	 
  }//if(!isinit_)

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
/*-----------------------------------------------------------------------*
 | Initialize (public) Setting up geometric variables for beam2r elements|
 *-----------------------------------------------------------------------*/ 
int DRT::ELEMENTS::Beam2rRegister::Initialize(DRT::Discretization& dis)
{

 //loop through all elements
 for (int num=0; num<  dis.NumMyColElements(); ++num)
  {    
    //in case that current element is not a beam2r element there is nothing to do and we go back
    //to the head of the loop
    if (dis.lColElement(num)->Type() != DRT::Element::element_beam2r) continue;
    //if we get so far current element is a beam2r element and  we get a pointer at it
    DRT::ELEMENTS::Beam2r* currele = dynamic_cast<DRT::ELEMENTS::Beam2r*>(dis.lColElement(num));
    if (!currele) dserror("cast to Beam2r* failed");
    
    //reference node position
    vector<double> xrefe;
    
    int nnode= currele->NumNode();
    //resize xrefe for the number of coordinates we need to store
    xrefe.resize(2*nnode);

    //getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {   
      for (int k=0; k<nnode; k++) //element has k nodes
        for(int l= 0; l < 2; l++)// element node has two coordinates x and z
          xrefe[k*2 + l] = currele->Nodes()[k]->X()[l];
    }
    //resize alpha_, alphamass_ and theta0_ so they can each store 1 value at each GP
    //it can be proven, that we will always need one GP less than nodes zu underintagrate in order to get rid of shear locking
    currele->alpha_.resize(nnode-1);
    currele->alphamass_.resize(nnode);
    currele->theta0_.resize(nnode-1);
    currele->floc_.resize(nnode);
    
    //SetUpReferenceGeometry is a templated function
    switch(nnode)
    {
  		case 2:  		
  		{	
  			currele->SetUpReferenceGeometry<2>(xrefe);
  			break;
  		}
  		case 3:
  		{
  			currele->SetUpReferenceGeometry<3>(xrefe);
  			break;
  		}
  		case 4:
  		{
  			currele->SetUpReferenceGeometry<4>(xrefe);
  			break;
  		} 
  		case 5:
  		{
  			currele->SetUpReferenceGeometry<5>(xrefe);
  			break;
  		}   		
  		default:
  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");	
  	} 
 
  } //for (int num=0; num<dis_.NumMyColElements(); ++num)

  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R
