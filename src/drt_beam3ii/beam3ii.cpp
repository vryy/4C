/*!----------------------------------------------------------------------
\file beam3ii.cpp
\brief three dimensional nonlinear corotational Timoshenko beam element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3II
#ifdef CCADISCRET

#include "beam3ii.H"
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
DRT::ELEMENTS::Beam3ii::Beam3ii(int id, int owner) :
DRT::Element(id,element_beam3ii,owner),
isinit_(false),
nodeI_(0),
nodeJ_(0),
crosssec_(0),
crosssecshear_(0),
Iyy_(0),
Izz_(0),
Irr_(0),
jacobi_(0),
jacobimass_(0),
jacobinode_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ii::Beam3ii(const DRT::ELEMENTS::Beam3ii& old) :
 DRT::Element(old),
 isinit_(old.isinit_),
 Qconv_(old.Qconv_),
 Qold_(old.Qold_),
 Qnew_(old.Qnew_),
 dispthetaconv_(old.dispthetaconv_),
 dispthetaold_(old.dispthetaold_),
 nodeI_(old.nodeI_),
 nodeJ_(old.nodeJ_),
 crosssec_(old.crosssec_),
 crosssecshear_(old.crosssecshear_),
 Iyy_(old.Iyy_),
 Izz_(old.Izz_),
 Irr_(old.Irr_),
 jacobi_(old.jacobi_),
 jacobimass_(old.jacobimass_),
 jacobinode_(old.jacobinode_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam3ii and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam3ii::Clone() const
{
  DRT::ELEMENTS::Beam3ii* newelement = new DRT::ELEMENTS::Beam3ii(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3ii::~Beam3ii()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::Print(ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Beam3iiRegister (public)               cyron 01/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Beam3ii::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Beam3iiRegister(Type()));
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam3ii::Shape() const
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
void DRT::ELEMENTS::Beam3ii::Pack(vector<char>& data) const
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
  AddtoPack(data,jacobi_);
  AddtoPack(data,jacobimass_);
  AddtoPack(data,jacobinode_);
  AddtoPack(data,nodeI_);
  AddtoPack(data,nodeJ_);
  AddtoPack(data,crosssec_);
  AddtoPack(data,crosssecshear_); 
  AddtoPack(data,isinit_);
  AddtoPack(data,Irr_);
  AddtoPack(data,Iyy_);
  AddtoPack(data,Izz_);
  AddtoPack<4,1>(data,Qconv_);
  AddtoPack<4,1>(data,Qnew_);
  AddtoPack<4,1>(data,Qold_);
  AddtoPack<3,1>(data,dispthetaconv_);
  AddtoPack<3,1>(data,dispthetaold_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::Unpack(const vector<char>& data)
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
  
  
  //extract all class variables of beam3ii element 
  ExtractfromPack(position,data,jacobi_);
  ExtractfromPack(position,data,jacobimass_);
  ExtractfromPack(position,data,jacobinode_);
  ExtractfromPack(position,data,nodeI_);
  ExtractfromPack(position,data,nodeJ_);
  ExtractfromPack(position,data,crosssec_);
  ExtractfromPack(position,data,crosssecshear_);
  ExtractfromPack(position,data,isinit_);
  ExtractfromPack(position,data,Irr_);
  ExtractfromPack(position,data,Iyy_);
  ExtractfromPack(position,data,Izz_); 
  ExtractfromPack<4,1>(position,data,Qconv_);
  ExtractfromPack<4,1>(position,data,Qnew_);
  ExtractfromPack<4,1>(position,data,Qold_);
  ExtractfromPack<3,1>(position,data,dispthetaconv_);
  ExtractfromPack<3,1>(position,data,dispthetaold_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Beam3ii::Lines()
{
  vector<RCP<Element> > lines(1);
  lines[0]= rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Beam3ii::MyGaussRule(int nnode, IntegrationType integrationtype)
{
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule1D_undefined;
  
  switch(nnode)
  {
    case 2:
    {     
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
      break;
    }
    case 3:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_3point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_2point;
          break;
        }
        case lobattointegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_lobatto3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 4:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_4point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 5:
    {
      switch(integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_5point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule =  DRT::UTILS::intrule_line_4point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }
  
  return gaussrule;
}

/*----------------------------------------------------------------------*
 | sets up geometric data from current nodal position as reference
 | position; this method can be used by the register class or when ever
 | a new beam element is generated for which some reference configuration
 | has to be stored; prerequesite for applying this method is that the
 | element nodes are already known (public)                   cyron 10/08|
 *----------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3ii::SetUpReferenceGeometry(const vector<double>& xrefe,const vector<double>& rotrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
   *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
   *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
   *the method and is true the element is initialized another time with respective xrefe and rotrefe;
   *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a 
   *second initilization in principle (e.g. for periodic boundary conditions*/

  if(!isinit_ || secondinit)
  {
    isinit_ = true;
    
    /*first the nodes for the reference triad \Lambda_r of the element are chosen according to eq. (6.2), Crisfield 1999;
     *note that the first node of the element in BACI is node 0 so that we need -1 in the end to convert from the notation
     *in Crisfield 1999 to the BACI convention*/
    nodeI_ = (int)floor(0.5*(NumNode()+1)) - 1;
    nodeJ_ = (int)floor(0.5*(NumNode()+2)) - 1;
    
    
    //resize and initialized STL vectors for rotational displacements so that they can store one value at each node
    dispthetaconv_.resize(nnode);
    dispthetaold_.resize(nnode);
    dispthetanew_.resize(nnode);
    for(int i=0; i<nnode; i++)
      for(int j=0; j<3; j++)
      {
        dispthetaconv_[i](j) = 0;
        dispthetaold_[i](j) = 0;
        dispthetanew_[i](j) = 0;
      }
      
    
    //resize STL vectors for Jacobi determinants so that they can store one value at each Gauss point
    jacobi_.resize(nnode-1);
    jacobimass_.resize(nnode);
    jacobinode_.resize(nnode);

    //create Matrix for the derivates of the shapefunctions at the GP
    LINALG::Matrix<1,nnode> shapefuncderiv;
    
    //create Matrix for the shapefunctions at the GP
    LINALG::Matrix<1,nnode> funct;
    
    //derivative of curve in physical space with respect to curve parameter xi \in [-1;1] on element level
    LINALG::Matrix<3,1> drdxi;
    
    //Get DiscretizationType
    DRT::Element::DiscretizationType distype = Shape();
    
    //Get the applied integrationpoints for underintegration
    DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));
    
    //Loop through all GPs and calculate jacobi the triads at the GPs
    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {   
      //Get position xi of GP
      const double xi = gausspoints.qxg[numgp][0];
      
      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);
      
      //Get shapefunctions at GP
      DRT::UTILS::shape_function_1D(funct,xi,distype);

      drdxi.Clear();

      //calculate vector dxdxi
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<3 ; dof++)
          drdxi(dof) += shapefuncderiv(node) * xrefe[3*node+dof];

      //Store Jacobi determinant with respect to reference configuration
      jacobi_[numgp]= drdxi.Norm2();
      
    }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    

    
    //Get the applied integrationpoints for exact integration of mass matrix
    DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode,gaussexactintegration));
        
    //Loop through all GPs and calculate jacobi and theta0
    for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
    {
        
      //Get position xi of GP
      const double xi = gausspointsmass.qxg[numgp][0];
      
      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

      drdxi.Clear();
      //calculate dx/dxi and dz/dxi
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<3; dof++)
          drdxi(dof)+=shapefuncderiv(node)*xrefe[3*node+dof];
      

      //Store Jacobi determinant with respect to reference configuration
      jacobimass_[numgp]= drdxi.Norm2();  
    
    }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
    
    
    //compute Jacobi determinant at gauss points for Lobatto quadrature (i.e. at nodes)
    for(int numgp=0; numgp< nnode; numgp++)
    {  
      //Get position xi of nodes
      const double xi = -1.0 + 2*numgp / (nnode - 1);
      
      //Get derivatives of shapefunctions at GP
      DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

      drdxi.Clear();
      //calculate dx/dxi and dz/dxi
      for(int node=0; node<nnode; node++)
        for(int dof=0; dof<3; dof++)
          drdxi(dof)+=shapefuncderiv(node)*xrefe[3*node+dof];
      
      //Store Jacobi determinant for each node (Jacobi determinant refers by definition always to the reference configuration)
      jacobinode_[numgp]= drdxi.Norm2(); 
    
    }//for(int numgp=0; numgp< nnode; numgp++)
    
  }//if(!isinit_)

	return;
	
}//DRT::ELEMENTS::Beam3ii::SetUpReferenceGeometry()

//------------- class Beam3iiRegister: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3iiRegister::Beam3iiRegister(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3iiRegister::Beam3iiRegister(
                               const DRT::ELEMENTS::Beam3iiRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3iiRegister* DRT::ELEMENTS::Beam3iiRegister::Clone() const
{
  return new DRT::ELEMENTS::Beam3iiRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3iiRegister::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Beam3iiRegister::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Beam3iiRegister::~Beam3iiRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3iiRegister::Print(ostream& os) const
{
  os << "Beam3iiRegister ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  Initialize (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3iiRegister::Initialize(DRT::Discretization& dis)
{		
	  //setting up geometric variables for beam3ii elements
	
	  for (int num=0; num<  dis.NumMyColElements(); ++num)
	  {
	    //in case that current element is not a beam3ii element there is nothing to do and we go back
	    //to the head of the loop
	    if (dis.lColElement(num)->Type() != DRT::Element::element_beam3ii) continue;
	
	    //if we get so far current element is a beam3ii element and  we get a pointer at it
	    DRT::ELEMENTS::Beam3ii* currele = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(dis.lColElement(num));
	    if (!currele) dserror("cast to Beam3ii* failed");
	
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
