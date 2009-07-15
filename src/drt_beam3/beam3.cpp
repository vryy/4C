/*!----------------------------------------------------------------------
\file beam3.cpp
\brief three dimensional nonlinear corotational Timoshenko beam element

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
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(int id, int owner) :
DRT::Element(id,element_beam3,owner),
data_(),
isinit_(false),
lrefe_(0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
gaussrule_(DRT::UTILS::intrule1D_undefined)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3::Beam3(const DRT::ELEMENTS::Beam3& old) :
 DRT::Element(old),
 data_(old.data_),
 isinit_(old.isinit_),
 lrefe_(old.lrefe_),
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 curvconv_(old.curvconv_),
 curvold_(old.curvold_),
 curvnew_(old.curvnew_),
 thetaconv_(old.thetaconv_),
 thetaold_(old.thetaold_),
 thetanew_(old.thetanew_),
 thetaprimeconv_(old.thetaprimeconv_),
 thetaprimeold_(old.thetaprimeold_),
 thetaprimenew_(old.thetaprimenew_),
 crosssec_(old.crosssec_),
 crosssecshear_(old.crosssecshear_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 alpha_(old.alpha_),
 alphamass_(old.alphamass_),
 floc_(old.floc_),
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
  //whether element has already been initialized
  AddtoPack(data,isinit_);
  //reference length
  AddtoPack(data,lrefe_);
  //central coordinate triad and related data
  AddtoPack(data,Qconv_);
  AddtoPack(data,Qold_);
  AddtoPack(data,Qnew_);
  AddtoPack(data,curvconv_);
  AddtoPack(data,curvold_);
  AddtoPack(data,curvnew_);
  AddtoPack(data,thetaconv_);
  AddtoPack(data,thetaold_);
  AddtoPack(data,thetanew_);
  AddtoPack(data,thetaprimeconv_);
  AddtoPack(data,thetaprimeold_);
  AddtoPack(data,thetaprimenew_);
  //cross section
  AddtoPack(data,crosssec_);
   //cross section with shear correction
  AddtoPack(data,crosssecshear_);
  //moments of inertia of area
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack(data,Irr_);
  //alpha for underintegration
  AddtoPack(data,alpha_);
  //alpha for complete integration
  AddtoPack(data,alphamass_);
  //stochastic forces
  AddtoPack(data,floc_);
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
  //whether element has already been initialized
  ExtractfromPack(position,data,isinit_);
  //reference length
  ExtractfromPack(position,data,lrefe_);
  //central coordinate triad and related data
  ExtractfromPack(position,data,Qconv_);
  ExtractfromPack(position,data,Qold_);
  ExtractfromPack(position,data,Qnew_);
  ExtractfromPack(position,data,curvconv_);
  ExtractfromPack(position,data,curvold_);
  ExtractfromPack(position,data,curvnew_);
  ExtractfromPack(position,data,thetaconv_);
  ExtractfromPack(position,data,thetaold_);
  ExtractfromPack(position,data,thetanew_);
  ExtractfromPack(position,data,thetaprimeconv_);
  ExtractfromPack(position,data,thetaprimeold_);
  ExtractfromPack(position,data,thetaprimenew_);
  //cross section
  ExtractfromPack(position,data,crosssec_);
  //cross section with shear correction
  ExtractfromPack(position,data,crosssecshear_);
  //moments of inertia of area
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_);
  ExtractfromPack(position,data,Irr_);
  //alpha for underintegration
  ExtractfromPack(position,data,alpha_);
  //alpha for complete integration
  ExtractfromPack(position,data,alphamass_);
  //stochastic forces
  ExtractfromPack(position,data,floc_);
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


/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   cyron 10/08|
 *----------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::SetUpReferenceGeometry(const vector<double>& xrefe,const vector<double>& rotrefe)
{
  /*this method initialized geometric variables of the element; such an initialization can only be done once when the element is
   * generated and never again (especially not in the frame of a restart); to make sure that this requirement is not violated this
   * method will initialize the geometric variables if the class variable isinit_ == false and afterwards set this variable to
   * isinit_ = true; if this method is called and finds alreday isinit_ == true it will just do nothing*/

  if(!isinit_)
  {
    isinit_ = true;
    
  //Set the applied Gaussrule ( It can be proven that we need 1 GP less than nodes to integrate exact )
  //note: we use a static cast for the enumeration here cf. Practical C++ Programming p.185
  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);

  //create Matrix for the derivates of the shapefunctions at the GP
	LINALG::Matrix<1,nnode> shapefuncderiv;
	
	//create Matrix for the shapefunctions at the GP
	LINALG::Matrix<1,nnode> funct;
	
	//Get DiscretizationType
	DRT::Element::DiscretizationType distype = Shape();
	
	//Get the applied integrationpoints
	DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule_);
	
	//Loop through all GPs and calculate alpha the triads at the GPs
	for(int numgp=0; numgp < gausspoints.nquad; numgp++)
	{
	  	
		//Get position xi of GP
		const double xi = gausspoints.qxg[numgp][0];
		
		//Get derivatives of shapefunctions at GP
		DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);
		
		//Get shapefunctions at GP
		DRT::UTILS::shape_function_1D(funct,xi,distype);
	
		//triad in reference configuration at GP
		LINALG::Matrix<3,3> Tref;
	
	    //length in reference configuration
	    lrefe_ = pow(pow(xrefe[0]-xrefe[3*nnode-3],2)+pow(xrefe[1]-xrefe[3*nnode-2],2)+pow(xrefe[2]-xrefe[3*nnode-1],2),0.5);
	
	    /*initial triad Tref = [t1,t2,t3] is set in a way for which we don`t have strains in reference configuration*/
	    LINALG::Matrix<3,1> dxdxi;
	
	    dxdxi.Clear();
	    thetaconv_[numgp].Clear();
	    thetaprimeconv_[numgp].Clear();
	
	    //calculate vector dxdxi
	    for(int node=0; node<nnode; node++)
	    {
	    	for(int dof=0; dof<3 ; dof++)
	    	{
		    	dxdxi(dof) += shapefuncderiv(node) * xrefe[3*node+dof];
			    thetaconv_[numgp](dof) += funct(node) * rotrefe[3*node+dof];
			    thetaprimeconv_[numgp](dof) += shapefuncderiv(node) * rotrefe[3*node+dof]; 		
	    	}//for(int dof=0; dof<3 ; dof++)
	    }//for(int node=0; node<nnode; node++)
	
	    //Store length factor for every GP
	    //note: the length factor alpha replaces the determinant and refers to the reference configuration by definition
	    alpha_[numgp]= pow(pow( dxdxi(0) ,2.0) + pow( dxdxi(1) ,2.0) + pow(dxdxi(2) ,2.0) ,0.5);	

	    for (int k=0; k<3; k++)
	    {
	  		//t1 axis points in positive direction along xi and is a unit vector
	  		Tref(k,0)=dxdxi(k)/alpha_[numgp];
	    }
	
	    //t2 is a unit vector in the x2x3-plane orthogonal to t1
	    Tref(0,1) = 0;
	    //if t1 is parallel to the x1-axis t2 is set parallel to the x2-axis
	    if (Tref(1,0) == 0 && Tref(2,0) == 0)
	    {
	         Tref(1,1) = 1;
	         Tref(2,1) = 0;
	    }
	
	    //otherwise t2 is calculated from the scalar product with t1
	    else
	    {
	        //setting t2(0)=0 and calculating other elements by setting scalar product t1 o t2 to zero
	        double lin1norm = pow(pow(Tref(1,0),2)+pow(Tref(2,0),2),0.5);
	        Tref(1,1) = -Tref(2,0)/lin1norm;
	        Tref(2,1) =  Tref(1,0)/lin1norm;
	    }	
	           	
        //calculating t3 by crossproduct t1 x t2
        Tref(0,2) = Tref(1,0)*Tref(2,1)-Tref(1,1)*Tref(2,0);
        Tref(1,2) = Tref(2,0)*Tref(0,1)-Tref(2,1)*Tref(0,0);
        Tref(2,2) = Tref(0,0)*Tref(1,1)-Tref(0,1)*Tref(1,0);

        /*the center triad in reference configuration is stored as a quaternion whose equivalent would be the rotation
        * from the identity matrix into the reference configuration*/
        triadtoquaternion(Tref,Qconv_[numgp]);

        //the here employed beam element does not need data about the current position of the nodal directors so that
        //initilization of those can be skipped (the nodal displacements handeled in beam3_evaluate.cpp are not the current angles,
        //but only the differences between current angles and angles in reference configuration, respectively. Thus the
        //director orientation in reference configuration cancels out and can be assumed to be zero without loss of
        //generality

        curvconv_[numgp].Clear();
    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)
	
	//Now all triads have been calculated and Qold_ and Qnew_ can be updated
    Qold_ = Qconv_;
    Qnew_ = Qconv_;

    curvold_ = curvconv_;
    curvnew_ = curvconv_;

    thetaold_ = thetaconv_;
    thetanew_ = thetaconv_;

    thetaprimeold_ = thetaprimeconv_;
    thetaprimenew_ = thetaprimeconv_;

	//Now we get the integrationfactor alphamass_ for a complete integration of the massmatrix therefor we increase the gaussrule by 1
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
	
	    LINALG::Matrix<3,1> dxdximass;
	
	    dxdximass.Clear();
	    //calculate dx/dxi and dz/dxi
	    for(int node=0; node<nnode; node++)
	    {
		    dxdximass(0)+=shapefuncderiv(node)*xrefe[3*node];
		    dxdximass(1)+=shapefuncderiv(node)*xrefe[3*node+1];
		    dxdximass(2)+=shapefuncderiv(node)*xrefe[3*node+2];
		
	    }//for(int node=0; node<nnode; node++)
	
	    //Store length factor for every GP
	    //note: the length factor alpha replaces the determinant and refers by definition always to the reference configuration
	    alphamass_[numgp]= pow(pow( dxdximass(0) ,2.0) + pow( dxdximass(1) ,2.0) + pow(dxdximass(2) ,2.0) ,0.5);	
	
	}//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)

	//reset gaussrule_ for further calculations
	gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);
	
	return;

  }//if(!isinit_)

}//DRT::ELEMENTS::Beam3::SetUpReferenceGeometry()

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

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3Register::Initialize(DRT::Discretization& dis)
{		
	  //setting up geometric variables for beam3 elements
	
	  for (int num=0; num<  dis.NumMyColElements(); ++num)
	  {
	    //in case that current element is not a beam3 element there is nothing to do and we go back
	    //to the head of the loop
	    if (dis.lColElement(num)->Type() != DRT::Element::element_beam3) continue;
	
	    //if we get so far current element is a beam3 element and  we get a pointer at it
	    DRT::ELEMENTS::Beam3* currele = dynamic_cast<DRT::ELEMENTS::Beam3*>(dis.lColElement(num));
	    if (!currele) dserror("cast to Beam3* failed");
	
	    //reference node position
	    vector<double> xrefe;
	    vector<double> rotrefe;
	    const int nnode= currele->NumNode();
	
	    //resize xrefe for the number of coordinates we need to store
	    xrefe.resize(3*nnode);
	    rotrefe.resize(3*nnode);
	
	    //getting element's nodal coordinates and treating them as reference configuration
	    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
	      dserror("Cannot get nodes in order to compute reference configuration'");
	    else
	    {
	      for (int node=0; node<nnode; node++) //element has k nodes
	        for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
	        {
	        	xrefe[node*3 + dof] = currele->Nodes()[node]->X()[dof];
	        	rotrefe[node*3 + dof]= 0.0;
	        }	
	    }
	
	    //resize alpha_, alphamass_ and theta0_ so they can each store 1 value at each GP
	    currele->alpha_.resize(nnode-1);
	    currele->alphamass_.resize(nnode);
	    currele->Qconv_.resize((nnode-1));
	    currele->Qold_.resize((nnode-1));
	    currele->Qnew_.resize((nnode-1));
	    currele->curvconv_.resize((nnode-1));
	    currele->curvold_.resize((nnode-1));
	    currele->curvnew_.resize((nnode-1));
	    currele->thetaconv_.resize((nnode-1));
	    currele->thetaold_.resize((nnode-1));
	    currele->thetanew_.resize((nnode-1));
	    currele->thetaprimeconv_.resize((nnode-1));
	    currele->thetaprimeold_.resize((nnode-1));
	    currele->thetaprimenew_.resize((nnode-1));

	    //SetUpReferenceGeometry is a templated function
	    switch(nnode)
	    {
	  		case 2:  		
	  		{	
	  			currele->SetUpReferenceGeometry<2>(xrefe,rotrefe);
	  			break;
	  		}
	  		case 3:
	  		{
	  			currele->SetUpReferenceGeometry<3>(xrefe,rotrefe);
	  			break;
	  		}
	  		case 4:
	  		{
	  			currele->SetUpReferenceGeometry<4>(xrefe,rotrefe);
	  			break;
	  		}
	  		case 5:
	  		{
	  			currele->SetUpReferenceGeometry<5>(xrefe,rotrefe);
	  			break;
	  		}   		
	  		default:
	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");	
	  	}
	
	  } //for (int num=0; num<dis_.NumMyColElements(); ++num)

	  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
