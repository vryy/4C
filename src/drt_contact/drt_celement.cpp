/*!----------------------------------------------------------------------
\file drt_celement.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_celement.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement::CElement(int id, ElementType etype, int owner, 
                    const DRT::Element::DiscretizationType& shape, 
                    const int numnode,
                    const int* nodeids,
                    const bool isslave) :
DRT::Element(id,etype,owner),
shape_(shape),
isslave_(isslave)
{
  SetNodeIds(numnode,nodeids);
  RefArea()=0.0;
  Area()=RefArea();
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement::CElement(const CONTACT::CElement& old) :
DRT::Element(old),
shape_(old.shape_),
isslave_(old.isslave_),
refarea_(old.refarea_),
area_(old.area_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CElement* CONTACT::CElement::Clone() const
{
  CONTACT::CElement* newele = new CONTACT::CElement(*this);
  return newele;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CElement& element)
{
  element.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Print(ostream& os) const
{
  os << "Contact Element ";
  DRT::Element::Print(os);
  if (isslave_) os << " Slave  Side ";
  else          os << " Master Side ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Element
  vector<char> basedata(0);
  DRT::Element::Pack(basedata);
  AddtoPack(data,basedata);
  // add shape_
  AddtoPack(data,shape_);
  // add isslave_
  AddtoPack(data,isslave_);
  // add refarea_
  AddtoPack(data,refarea_);
  // add area_
  AddtoPack(data,area_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Element::Unpack(basedata);
  // shape_
  ExtractfromPack(position,data,shape_);
  // isslave_
  ExtractfromPack(position,data,isslave_);
  // refarea_
  ExtractfromPack(position,data,refarea_);
  // area_
  ExtractfromPack(position,data,area_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate element (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
int CONTACT::CElement::Evaluate(ParameterList& params,
                                DRT::Discretization&      discretization,
                                vector<int>&              lm,
                                Epetra_SerialDenseMatrix& elemat1,
                                Epetra_SerialDenseMatrix& elemat2,
                                Epetra_SerialDenseVector& elevec1,
                                Epetra_SerialDenseVector& elevec2,
                                Epetra_SerialDenseVector& elevec3)
{
  dserror("CONTACT::CElement::Evaluate not yet impl.");
  return -1;
}

/*----------------------------------------------------------------------*
 |  Get local coordinates for local node id                   popp 12/07|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::LocalCoordinatesOfNode(int lid, double* xi)
{
	// 2D linear and quadratic case
	if (Shape()==line2 || Shape()==line3)
	{
		if (lid==0)
			xi[0]=-1.0;
		else if (lid==1)
			xi[0]= 1.0;
		else if (lid==2)
			xi[0]= 0.0;
		else
			dserror("ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range",lid,Id());
		
		// we are in the 2D case here!
		xi[1]=0.0;
	}
	
	// unknown case (3D not implemented yet)
	else
		dserror("ERROR: LocalCoordinatesOfNode called for unknown element type");
		
	return true;
}

/*----------------------------------------------------------------------*
 |  Get local numbering for global node id                    popp 12/07|
 *----------------------------------------------------------------------*/
int CONTACT::CElement::GetLocalNodeId(int nid)
{
	int lid = -1;	
	
	// look for global ID nid in element's modes
	for (int i=0;i<NumNode();++i)
	{
		if ( *(NodeIds()+i) == nid )
		{
			lid=i;
			break;
		}
	}
	
	if (lid<0)
		dserror("ERROR: Cannot find node % in segment %", nid, Id());
	
	return lid;
}

/*----------------------------------------------------------------------*
 |  Build element normal at node                              popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::BuildNormalAtNode(int nid, vector<double>& n)
{
	if (n.size()!=3) n.resize(3);
	
	// find this node in my list of nodes and get local numbering
	int lid = GetLocalNodeId(nid);
	
	// get local coordinates for this node
	double xi[2];
	LocalCoordinatesOfNode(lid,xi);
	
	// build an outward unit normal at xi and return it
	ComputeNormalAtXi(xi,n);

	return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal at loc. coord. xi                  popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::ComputeNormalAtXi(double* xi, vector<double>& n)
{
	int nnodes = NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
	vector<double> g(3);
	/*
	// test dual shape function values and derivatives at xi
	EvaluateShape_Dual1D(xi, val, deriv, nnodes);
	if (Shape()==line2)
		cout << "DualShapeFct:" << endl << val[0] << "\t" << val[1] << endl;
	if (Shape()==line3)
		cout << "DualShapeFct:" << endl << val[0] << "\t" << val[1] << "\t" << val[2] << endl;
	*/
	// get shape function values and derivatives at xi
	EvaluateShape_1D(xi, val, deriv, nnodes);

	// get coordinates of element nodes
	LINALG::SerialDenseMatrix coord(3,nnodes);
	coord = GetNodalCoords();
	
	// build basis vector g
	for (int i=0;i<nnodes;++i)
	{
		g[0]+=deriv[i]*coord(0,i);
		g[1]+=deriv[i]*coord(1,i);
		g[2]+=deriv[i]*coord(2,i);
	}
	
	// only 2D case implemented, normal easy to compute from g
	n[0] = g[1];
	n[1] =-g[0];
	n[2] = 0.0;
	
	double length = sqrt(n[0]*n[0]+n[1]*n[1]);
	if (length==0.0)
		dserror("ERROR: ComputeNormalAtXi computes normal of length zero!");
	
	// create unit normal (division by length)
	for (int i=0;i<3;++i)
		n[i]/=length;
	
	return;
}
/*----------------------------------------------------------------------*
 |  Get nodal coordinates of the element                      popp 01/08|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix CONTACT::CElement::GetNodalCoords()
{
	int nnodes = NumNode();
	DRT::Node** mynodes = Nodes();
	LINALG::SerialDenseMatrix coord(3,nnodes);
	if (!mynodes)
		dserror("ERROR: GetNodalCoords: Null pointer!");
	
	for (int i=0;i<nnodes;++i)
	{
		CNode* mycnode = static_cast<CNode*> (mynodes[i]);
		if (!mycnode)
			dserror("ERROR: GetNodalCoords: Null pointer!");
		coord(0,i) = mycnode->xspatial()[0];
		coord(1,i) = mycnode->xspatial()[1];
		coord(2,i) = mycnode->xspatial()[2];
	}
	
	return coord;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant - LINEAR / QUAD 1D          popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::Jacobian_1D(const vector<double>& val,
																	 		const vector<double>& deriv,
  																		const LINALG::SerialDenseMatrix& coord)
{
	double jac = 0.0;
	
	// 2D linear case (2noded line element)
	if (Shape()==line2)
		jac = Area()/2;
	
	// 2D quadratic case (3noded line element)
	else if (Shape()==line3)
	{
		double g[3] = {0.0, 0.0, 0.0};
		for (int i=0;i<static_cast<int>(val.size());++i)
		{
			g[0] += deriv[i]*coord(0,i);
			g[1] += deriv[i]*coord(1,i);
			g[2] += deriv[i]*coord(2,i);
		}
		jac = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
	}
	
	// unknown case
	else
		dserror("ERROR: Jacobian_1D called for unknown element type!");
		
	return jac;
}

/*----------------------------------------------------------------------*
 |  Compute length (in 3D area) of the element                popp 12/07|
 *----------------------------------------------------------------------*/
double CONTACT::CElement::ComputeArea()
{
	double area = 0.0;
	
	// 2D linear case (2noded line element)
	if (Shape()==line2)
	{
		// no integration necessary (constant Jacobian)
		int nnodes = NumNode();
		LINALG::SerialDenseMatrix coord(3,nnodes);
		coord = GetNodalCoords();
		
		// build vector between the two nodes
		double tang[3] = {0.0, 0.0, 0.0};
		for (int k=0;k<3;++k)
		{
			tang[k]=coord(k,1)-coord(k,0);
		}
		area=sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]);
	}
	
	// 2D quadratic case (3noded line element)
	else if (Shape()==line3)
	{
		// Gauss quadrature
		CONTACT::Integrator integrator(CONTACT_NGP,true);
		int nnodes = NumNode();
		double detg = 0.0;
		vector<double> val(nnodes);
		vector<double> deriv(nnodes);
			
		LINALG::SerialDenseMatrix coord(3,nnodes);
		coord = GetNodalCoords();
		
			// loop over all Gauss points, build Jacobian and compute area
		for (int j=0;j<integrator.nGP();++j)
		{
			double gpc[2] = {integrator.Coordinate(j), 0.0};
			EvaluateShape_1D(gpc, val, deriv, nnodes);
			detg = Jacobian_1D(val,deriv,coord);			
			area+= integrator.Weight(j)*detg;
		}	
	}
	
	// other cases (3D) not implemented yet
	else
		dserror("ERROR: Area computation not implemented for this type of CElement");
	
	return area;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions - LINEAR / QUAD 1D               popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::EvaluateShape_1D(const double* xi, vector<double>& val,
																				 vector<double>& deriv, const int valdim)
{
	if (!xi)
		dserror("ERROR: EvaluateShape_1D called with xi=NULL");
	
	// 2D linear case (2noded line element)
	if ((valdim==2)&& (Shape()==line2))
	{
		val[0] = 0.5*(1-xi[0]);
		val[1] = 0.5*(1+xi[0]);
		
		deriv[0] = -0.5;
		deriv[1] =  0.5;
	}

	// 2D quadratic case (3noded line element)
	else if ((valdim==3) && (Shape()==line3))
	{
		val[0] = 0.5*xi[0]*(xi[0]-1);
		val[1] = 0.5*xi[0]*(xi[0]+1);
		val[2] = (1-xi[0])*(1+xi[0]);
			
		deriv[0] = xi[0]-0.5;
		deriv[1] = xi[0]+0.5;
		deriv[2] = -2*xi[0];
	}
	
	// unknown case
	else
		dserror("ERROR: EvaluateShape_1D called for unknown CElement type");

	return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate dual shape functions - QUAD 1D                   popp 12/07|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::EvaluateShape_Dual1D(const double* xi, vector<double>& val,
																						 vector<double>& deriv, const int valdim)
{
	if (!xi)
		dserror("ERROR: EvaluateShape_Dual1D called with xi=NULL");
	
	// 2D linear case (2noded line element)
	if ((valdim==2) && (Shape()==line2))
	{
			val[0] = 0.5*(1-3*xi[0]);
			val[1] = 0.5*(1+3*xi[0]);
			
			deriv[0] = -1.5;
			deriv[1] =  1.5;
	}
	
	// 2D quadratic case (3noded line element)
	else if ((valdim==3) && (Shape()==line3))
	{
		// establish fundamental data	
		int nnodes = NumNode();
		double detg = 0.0;
	
		LINALG::SerialDenseMatrix coord(3,nnodes);
		coord = GetNodalCoords();
		
		// compute entries to bi-ortho matrices Me/De with Gauss quadrature
		CONTACT::Integrator integrator(CONTACT_NGP,true);
		
		Epetra_SerialDenseMatrix Me(nnodes,nnodes);
		Epetra_SerialDenseMatrix De(nnodes,nnodes);
		
		for (int i=0;i<integrator.nGP();++i)
		{
			double gpc[2] = {integrator.Coordinate(i), 0.0};
			EvaluateShape_1D(gpc, val, deriv, nnodes);
			detg = Jacobian_1D(val,deriv,coord);
			
			for (int j=0;j<nnodes;++j)
				for (int k=0;k<nnodes;++k)
				{
					Me(j,k)+=integrator.Weight(i)*val[j]*val[k]*detg;
					De(j,k)+=(j==k)*integrator.Weight(i)*val[j]*detg;
				}	
		}
		
		// invert bi-ortho matrix Me
		LINALG::SymmetricInverse(Me,nnodes);
		
		// get solution matrix with dual parameters
		Epetra_SerialDenseMatrix Ae(nnodes,nnodes);
		Ae.Multiply('N','N',1.0,De,Me,0.0);
		
		// evaluate dual shape functions at loc. coord. xi
		// need standard shape functions at xi first
		EvaluateShape_1D(xi, val, deriv, nnodes);
		
		vector<double> valtemp(nnodes);
		vector<double> derivtemp(nnodes);
		for (int i=0;i<nnodes;++i)
		{
			valtemp[i]=0.0;
			derivtemp[i]=0.0;
			for (int j=0;j<nnodes;++j)
			{
				valtemp[i]+=Ae(i,j)*val[j];
				derivtemp[i]+=Ae(i,j)*deriv[j];
			}
		}
		val=valtemp;
		deriv=derivtemp;
	}
	
	// unknown case
	else
		dserror("ERROR: EvaluateShape_Dual1D called for unknown element type");
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords                  popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::LocalToGlobal(const double* xi, double* globcoord,
																			bool inttype)
{
	// check input
	if (!xi)
		dserror("ERROR: LocalToGlobal called with xi=NULL");
	if (!globcoord)
		dserror("ERROR: LocalToGlobal called with globcoord=NULL");
	if (Shape()!=line3 && Shape()!=line2)
		dserror("ERROR: LocalToGlobal called for CEl type != line3");
	
	// collect fundamental data
	int nnodes = NumNode();
	DRT::Node** mynodes = Nodes();
	if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");
	LINALG::SerialDenseMatrix coord(3,nnodes);
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
	
	// Evaluate shape, get nodal coords  and interpolate global coords
	EvaluateShape_1D(xi, val, deriv, nnodes);
	for (int i=0;i<3;++i)
		globcoord[i]=0.0;
	for (int i=0;i<nnodes;++i)
	{
		CNode* mycnode = static_cast<CNode*> (mynodes[i]);
		if (!mycnode) dserror("ERROR: LocalToGlobal: Null pointer!");
		coord(0,i) = mycnode->xspatial()[0];
		coord(1,i) = mycnode->xspatial()[1];
		coord(2,i) = mycnode->xspatial()[2];
		
		if (inttype)
		{
			// use shape function values for interpolation
			globcoord[0]+=val[i]*coord(0,i);
			globcoord[1]+=val[i]*coord(1,i);
			globcoord[2]+=val[i]*coord(2,i);
		}
		else
		{
			// use shape function derivatives for interpolation
			globcoord[0]+=deriv[i]*coord(0,i);
			globcoord[1]+=deriv[i]*coord(1,i);
			globcoord[2]+=deriv[i]*coord(2,i);
		}
	}
	
	return true;
}

#endif  // #ifdef CCADISCRET
